module mod_particle_sorting
  use iso_fortran_env, only: int32, real64, int8
  use mod_particles,   only: ParticleSet
  use omp_lib,         only: omp_get_max_threads
  implicit none
  private

  public :: sort_particles_by_cell
  public :: redistribute_species_by_cell
  public :: compute_particle_cell_ids
  public :: cell_index_from_position
  public :: check_particles_are_sorted
  public :: check_cell_indexing
  public :: get_cell_particle_range

  ! --- persistent scratch buffers for sort_particles_by_cell -----------
  ! Buffers only grow (doubling, same convention as ParticleSet%ensure_capacity)
  ! and are never deallocated mid-run, avoiding a large heap
  ! allocate+deallocate on every single sort.
  !
  ! part%sp is not staged/reordered here at all, see the note in
  ! sort_particles_by_cell below - it's provably constant across every
  ! particle in a given ParticleSet, so sorting it is a no-op.
  !
  ! The 7 core real fields (x,y,z,vx,vy,vz,w - together ~69% of this
  ! routine's original cache-miss cost) are staged through ONE packed
  ! AoS record (RealFields7) instead of 7 separate SoA arrays: the
  ! reorder scatter writes to a target position that's essentially
  ! random relative to the source index (that's the point of sorting),
  ! so with 7 separate arrays every write is an independent cache-miss-
  ! prone random write into a different array; packed into one 56-byte
  ! record (7 reals, no padding - unlike a mixed-type record, an
  ! all-real64 type needs no alignment padding) the scatter touches one
  ! cache line per particle instead of up to seven. All fields the same
  ! size/kind is deliberate: the smaller int32/int8 metadata fields
  ! (cell_id, flag_dead, flag_cex - together ~12% of the original cost)
  ! are cheap enough, and different alignment enough, to leave as plain
  ! SoA arrays rather than risk padding bloat by folding them in too.
  !
  ! IMPORTANT: unpack this back into part%x/y/z/vx/vy/vz/w with ONE
  ! explicit loop that reads staging(pos) once and writes all 7 fields
  ! before moving to the next pos (see sort_particles_by_cell) - NOT as
  ! 7 separate whole-array-section statements
  ! (part%pv(1,1:np)=staging(1:np)%x, etc). The latter was tried first and
  ! measured as a net regression (sorting ms roughly doubled on the
  ! matched benchmark): each field-extraction statement independently
  ! re-scans the ENTIRE staging array (which is far larger than any
  ! cache) to pull out just its one field, so 7 separate statements
  ! cost ~7x the DRAM traffic of a single combined pass. Reading the
  ! whole packed record once per iteration and writing all 7
  ! destination arrays immediately keeps total traffic proportional to
  ! one pass, matching the original unpack's cost profile.
  !
  ! mod_state%sort_particles_local parallelizes over (ptype,iproc) with
  ! !$omp parallel do, so these are declared threadprivate: every OpenMP
  ! thread gets and keeps its own persistent copy, growing independently,
  ! with no risk of two threads stomping on a shared buffer.
  type :: RealFields7
    real(real64) :: x, y, z, vx, vy, vz, w
  end type RealFields7

  integer(int32), save :: scratch_cap = 0_int32
  integer(int32), save :: cell_scratch_cap = 0_int32

  integer(int32),     allocatable, save :: next_slot(:)
  integer(int32),     allocatable, save :: cell_id_new(:)
  type(RealFields7),  allocatable, save :: staging(:)
  integer(int8),      allocatable, save :: flag_dead_new(:)
  integer(int32),     allocatable, save :: flag_cex_new(:)

  !$omp threadprivate(scratch_cap, cell_scratch_cap, next_slot, cell_id_new, &
  !$omp                staging, flag_dead_new, flag_cex_new)

  ! --- shared (NOT threadprivate) scratch for redistribute_species_by_cell
  ! One staging buffer for a whole species across all iproc, plus per-cell
  ! global counts/starts. Grow-only, like the per-thread buffers above.
  integer(int32), save :: gscratch_cap = 0_int32
  integer(int32), save :: gcell_cap    = 0_int32
  type(RealFields7), allocatable, save :: gstage(:)
  integer(int32),    allocatable, save :: gcell(:)
  integer(int8),     allocatable, save :: gdead(:)
  integer(int32),    allocatable, save :: gcex(:)
  integer(int32),    allocatable, save :: gstart(:)

contains

  subroutine ensure_sort_scratch(np)
    integer(int32), intent(in) :: np
    integer(int32) :: new_cap

    if (np <= scratch_cap) return

    new_cap = max(np, max(1_int32, 2_int32*scratch_cap))

    if (allocated(staging)) then
      deallocate(staging, cell_id_new, flag_dead_new, flag_cex_new)
    end if

    allocate(staging(new_cap))
    allocate(cell_id_new(new_cap))
    allocate(flag_dead_new(new_cap))
    allocate(flag_cex_new(new_cap))

    scratch_cap = new_cap
  end subroutine ensure_sort_scratch

  subroutine ensure_sort_cell_scratch(ncells)
    integer(int32), intent(in) :: ncells

    if (ncells <= cell_scratch_cap) return

    if (allocated(next_slot)) deallocate(next_slot)
    allocate(next_slot(ncells))
    cell_scratch_cap = ncells
  end subroutine ensure_sort_cell_scratch

  pure integer(int32) function cell_index_from_position(x, y, z, h, n) result(ic)
    !=============================================================
    ! Return the 1D physical cell index corresponding to particle
    ! position, using the same convention as the legacy code:
    !
    !   ic = (ix-1) + n(1)*((iy-1) + n(2)*(iz-1)) + 1
    !
    ! where:
    !   ix = INT(x/hx) + 1
    !   iy = INT(y/hy) + 1
    !   iz = INT(z/hz) + 1
    !
    ! Valid cell ids are therefore:
    !   1 ... n(1)*n(2)*n(3)
    !
    ! Particle positions are clamped into the physical cell range.
    !=============================================================
    real(real64),   intent(in) :: x, y, z
    real(real64),   intent(in) :: h(3)
    integer(int32), intent(in) :: n(3)

    integer(int32) :: ix, iy, iz

    ix = int(x / h(1), int32) + 1_int32
    iy = int(y / h(2), int32) + 1_int32
    iz = int(z / h(3), int32) + 1_int32

    ix = max(1_int32, min(n(1), ix))
    iy = max(1_int32, min(n(2), iy))
    iz = max(1_int32, min(n(3), iz))

    ic = (ix - 1_int32) + n(1) * ((iy - 1_int32) + n(2) * (iz - 1_int32)) + 1_int32
  end function cell_index_from_position


  subroutine compute_particle_cell_ids(part, n, h)
    !=============================================================
    ! Compute cell_id(i) for each active particle in ParticleSet.
    ! Also fills:
    !   - cell_count(ic): number of particles in cell ic
    !   - cell_start(ic): first particle index of cell ic
    !
    ! This routine does NOT reorder the particle arrays.
    !=============================================================
    type(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)

    integer(int32) :: i, ic, ncells, lo, hi

    ncells = n(1) * n(2) * n(3)

    call part%ensure_cell_storage(ncells)

    ! Only cells cell_lo..cell_hi can be nonzero (see ParticleSet).
    call clear_cell_range(part)

    if (part%n <= 0_int32) return

    lo = huge(1_int32)
    hi = 0_int32
    do i = 1, part%n
      ic = cell_index_from_position(part%pv(1,i), part%pv(2,i), part%pv(3,i), h, n)
      part%cell_id(i)     = ic
      part%cell_count(ic) = part%cell_count(ic) + 1_int32
      lo = min(lo, ic)
      hi = max(hi, ic)
    end do
    part%cell_lo = lo
    part%cell_hi = hi

    ! cell_start is only meaningful where cell_count > 0; outside the
    ! occupied range it stays 0.
    part%cell_start(lo) = 1_int32
    do ic = lo+1, hi
      part%cell_start(ic) = part%cell_start(ic-1) + part%cell_count(ic-1)
    end do
  end subroutine compute_particle_cell_ids


  subroutine clear_cell_range(part)
    ! Zero cell_count/cell_start over the set's occupied range, leaving
    ! both all-zero, and mark the range empty.
    type(ParticleSet), intent(inout) :: part

    if (part%cell_lo <= part%cell_hi) then
      part%cell_count(part%cell_lo:part%cell_hi) = 0_int32
      part%cell_start(part%cell_lo:part%cell_hi) = 0_int32
    end if
    part%cell_lo = 1_int32
    part%cell_hi = 0_int32
  end subroutine clear_cell_range


  subroutine sort_particles_by_cell(part, n, h)
    !=============================================================
    ! Reorder particles so that particles belonging to the same cell
    ! are contiguous in memory.
    !
    ! On exit:
    !   - x/y/z, vx/vy/vz, w, sp are sorted by cell
    !   - cell_id(i) is sorted and nondecreasing
    !   - cell_count(ic) is valid
    !   - cell_start(ic) is valid
    !=============================================================
    type(ParticleSet), intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)

    integer(int32) :: np, ncells
    integer(int32) :: i, ic, pos

    if (.not. allocated(part%pv)) return

    np     = part%n
    ncells = n(1) * n(2) * n(3)

    call part%ensure_cell_storage(ncells)
    call compute_particle_cell_ids(part, n, h)

    if (np <= 1_int32) return

    call ensure_sort_scratch(np)
    call ensure_sort_cell_scratch(ncells)

    next_slot(1:ncells) = part%cell_start(1:ncells)

    ! part%sp is deliberately NOT staged/reordered here: since this
    ! ParticleSet holds exactly one species (indexed by ptype at the
    ! caller), part%sp(i) == the same species id for every i, always -
    ! every write site sets it to that same constant (mod_injection.f90,
    ! mod_restart.f90, mod_particleBC.f90/mod_particleMover.f90's
    ! secondary-electron creation all write sp(k)=ptype/species_id, and
    ! it is never read anywhere else to branch on). Reordering an array
    ! where every element is already identical is a no-op, so skipping it
    ! here saves a real cost (was ~10% of this routine's cache misses)
    ! for zero behavioral difference.
    do i = 1, np
      ic  = part%cell_id(i)
      pos = next_slot(ic)

      staging(pos)%x   = part%pv(1,i)
      staging(pos)%y   = part%pv(2,i)
      staging(pos)%z   = part%pv(3,i)
      staging(pos)%vx  = part%pv(4,i)
      staging(pos)%vy  = part%pv(5,i)
      staging(pos)%vz  = part%pv(6,i)
      staging(pos)%w   = part%w(i)
      cell_id_new(pos) = ic

      flag_dead_new(pos) = part%flag_dead(i)
      flag_cex_new(pos)  = part%flag_cex(i)

      next_slot(ic) = next_slot(ic) + 1_int32
    end do

    ! Unpack as ONE loop reading staging(pos) once per iteration - see the
    ! module header comment for why this must not be split into separate
    ! whole-array-section statements per field.
    do pos = 1, np
      part%pv(1,pos)  = staging(pos)%x
      part%pv(2,pos)  = staging(pos)%y
      part%pv(3,pos)  = staging(pos)%z
      part%pv(4,pos) = staging(pos)%vx
      part%pv(5,pos) = staging(pos)%vy
      part%pv(6,pos) = staging(pos)%vz
      part%w(pos)  = staging(pos)%w
    end do

    part%cell_id(1:np)   = cell_id_new(1:np)
    part%flag_dead(1:np) = flag_dead_new(1:np)
    part%flag_cex(1:np)  = flag_cex_new(1:np)

  end subroutine sort_particles_by_cell


  subroutine redistribute_species_by_cell(parts, n, h, species_id)
    !=============================================================
    ! Global counting sort of ONE species across all iproc, then deal
    ! the sorted particles back out in contiguous, equal-size chunks:
    ! iproc k gets global positions (k-1)*chunk+1 .. k*chunk. Since the
    ! cell index runs x fastest, then y, then z, each iproc ends up
    ! owning a thin z-slab of the domain - same as legacy part_sorting
    ! (Src/sorting.f90). A thread's field gathers (mover) and charge
    ! scatter (deposit) then touch only its slab of E/B/np instead of
    ! the whole grid, which keeps them in L2 and is what lets legacy's
    ! push scale ~2x with hyperthreading.
    !
    ! On exit each parts(k) is sorted by cell with valid cell_id,
    ! cell_count and cell_start (as after sort_particles_by_cell), and
    ! per-species counts per iproc are balanced to within one particle.
    !
    ! Must be called OUTSIDE any parallel region (opens its own).
    !=============================================================
    type(ParticleSet), intent(inout) :: parts(:)
    integer(int32),    intent(in)    :: n(3)
    real(real64),      intent(in)    :: h(3)
    integer(int32),    intent(in)    :: species_id

    integer(int32) :: nproc, ncells, ntot, chunk
    integer(int32) :: p, ic, i, g, g0, g1, m, acc, n_old
    integer(int32) :: plo(size(parts)), phi(size(parts))
    integer(int32) :: nblk, b, b0, b1, i0, i1

    nproc  = size(parts)
    ncells = n(1) * n(2) * n(3)

    ntot = 0_int32
    do p = 1, nproc
      if (allocated(parts(p)%pv)) ntot = ntot + max(0_int32, parts(p)%n)
    end do
    if (ntot <= 0_int32) return

    chunk = (ntot + nproc - 1_int32) / nproc

    if (ntot > gscratch_cap) then
      if (allocated(gstage)) deallocate(gstage, gcell, gdead, gcex)
      gscratch_cap = max(ntot, 2_int32*gscratch_cap)
      allocate(gstage(gscratch_cap), gcell(gscratch_cap), &
               gdead(gscratch_cap), gcex(gscratch_cap))
    end if
    if (ncells > gcell_cap) then
      if (allocated(gstart)) deallocate(gstart)
      allocate(gstart(ncells))
      gcell_cap = ncells
    end if

    ! 1. Per-iproc cell ids and per-cell counts (cell_start from this is
    !    overwritten below).
    !$omp parallel do private(p) schedule(static)
    do p = 1, nproc
      if (allocated(parts(p)%pv)) then
        call compute_particle_cell_ids(parts(p), n, h)
      else
        call parts(p)%ensure_cell_storage(ncells)
        parts(p)%cell_count = 0_int32
      end if
    end do
    !$omp end parallel do

    ! Occupied cell range of each set (from step 1). Sets are only read
    ! inside their own range - everything outside it is zero.
    do p = 1, nproc
      plo(p) = parts(p)%cell_lo
      phi(p) = parts(p)%cell_hi
    end do

    ! 2. Global per-cell totals, then exclusive prefix -> global start.
    !    Parallel over contiguous cell blocks; each block adds, in iproc
    !    order, the part of every set's range that overlaps it. (A per-cell
    !    loop over p here was outer-loop-vectorized by ifx 2025.0 into an
    !    AVX-512 gather with dest == index register, which faults SIGILL.)
    nblk = max(1_int32, min(ncells, 8_int32*int(omp_get_max_threads(), int32)))
    !$omp parallel do private(b,b0,b1,p,i0,i1) schedule(static)
    do b = 1, nblk
      b0 = int((int(ncells,8)*(b-1))/nblk, int32) + 1_int32
      b1 = int((int(ncells,8)*b)/nblk, int32)
      gstart(b0:b1) = 0_int32
      do p = 1, nproc
        i0 = max(b0, plo(p)); i1 = min(b1, phi(p))
        if (i0 <= i1) gstart(i0:i1) = gstart(i0:i1) + parts(p)%cell_count(i0:i1)
      end do
    end do
    !$omp end parallel do

    call exclusive_prefix_sum(gstart, ncells)

    ! 3. Per-(cell,iproc) write cursors, stored in each set's cell_start
    !    (only inside its range; rebuilt in step 6). gstart is advanced in
    !    place as the running cursor (not needed afterwards); iproc order
    !    within a cell keeps the result independent of thread scheduling.
    !$omp parallel do private(b,b0,b1,p,i0,i1) schedule(static)
    do b = 1, nblk
      b0 = int((int(ncells,8)*(b-1))/nblk, int32) + 1_int32
      b1 = int((int(ncells,8)*b)/nblk, int32)
      do p = 1, nproc
        i0 = max(b0, plo(p)); i1 = min(b1, phi(p))
        if (i0 <= i1) then
          parts(p)%cell_start(i0:i1) = gstart(i0:i1)
          gstart(i0:i1) = gstart(i0:i1) + parts(p)%cell_count(i0:i1)
        end if
      end do
    end do
    !$omp end parallel do

    ! 4. Scatter every particle to its global sorted slot. Each iproc owns
    !    its cursors, so no atomics are needed.
    !$omp parallel do private(p,i,ic,g) schedule(static)
    do p = 1, nproc
      if (.not. allocated(parts(p)%pv)) cycle
      do i = 1, parts(p)%n
        ic = parts(p)%cell_id(i)
        g  = parts(p)%cell_start(ic)
        parts(p)%cell_start(ic) = g + 1_int32

        gstage(g)%x  = parts(p)%pv(1,i)
        gstage(g)%y  = parts(p)%pv(2,i)
        gstage(g)%z  = parts(p)%pv(3,i)
        gstage(g)%vx = parts(p)%pv(4,i)
        gstage(g)%vy = parts(p)%pv(5,i)
        gstage(g)%vz = parts(p)%pv(6,i)
        gstage(g)%w  = parts(p)%w(i)
        gcell(g) = ic
        gdead(g) = parts(p)%flag_dead(i)
        gcex(g)  = parts(p)%flag_cex(i)
      end do
    end do
    !$omp end parallel do

    ! 5+6. Each iproc copies in its own contiguous chunk (so it first-
    !      touches its own memory) and rebuilds its cell_count/cell_start.
    !$omp parallel do private(p,g0,g1,m,i,g,ic,n_old) schedule(static)
    do p = 1, nproc
      g0 = (p - 1_int32) * chunk + 1_int32
      g1 = min(p * chunk, ntot)
      m  = max(0_int32, g1 - g0 + 1_int32)

      n_old = 0_int32
      if (allocated(parts(p)%pv)) n_old = parts(p)%n
      call parts(p)%ensure_capacity(m)
      parts(p)%n = m
      ! Back to all-zero cell lists (only the step-1 range was touched),
      ! then rebuild over the new, contiguous range.
      call clear_cell_range(parts(p))
      ! Slots this set vacates keep zeroed flags (invariant for slots > n).
      if (n_old > m) then
        parts(p)%flag_dead(m+1:n_old) = 0_int8
        parts(p)%flag_cex(m+1:n_old)  = 0_int32
      end if

      do i = 1, m
        g = g0 + i - 1_int32
        parts(p)%pv(1,i) = gstage(g)%x
        parts(p)%pv(2,i) = gstage(g)%y
        parts(p)%pv(3,i) = gstage(g)%z
        parts(p)%pv(4,i) = gstage(g)%vx
        parts(p)%pv(5,i) = gstage(g)%vy
        parts(p)%pv(6,i) = gstage(g)%vz
        parts(p)%w(i)    = gstage(g)%w
        parts(p)%sp(i)   = species_id
        ic = gcell(g)
        parts(p)%cell_id(i)   = ic
        parts(p)%flag_dead(i) = gdead(g)
        parts(p)%flag_cex(i)  = gcex(g)
        parts(p)%cell_count(ic) = parts(p)%cell_count(ic) + 1_int32
      end do

      if (m > 0_int32) then
        ! Chunk is sorted by cell, so its range is first..last cell.
        parts(p)%cell_lo = gcell(g0)
        parts(p)%cell_hi = gcell(g1)
        parts(p)%cell_start(parts(p)%cell_lo) = 1_int32
        do ic = parts(p)%cell_lo+1, parts(p)%cell_hi
          parts(p)%cell_start(ic) = parts(p)%cell_start(ic-1) + parts(p)%cell_count(ic-1)
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine redistribute_species_by_cell


  subroutine exclusive_prefix_sum(a, na)
    ! In place: a(i) <- 1 + sum(a(1:i-1)) (1-based starts). Blocked so all
    ! threads share the pass instead of one thread walking all ncells.
    integer(int32), intent(inout) :: a(:)
    integer(int32), intent(in)    :: na

    integer(int32), allocatable :: bsum(:)
    integer(int32) :: nb, b, i, i0, i1, acc, t

    nb = max(1, min(na, 8*omp_get_max_threads()))
    allocate(bsum(0:nb))
    bsum = 0_int32

    !$omp parallel do private(b,i,i0,i1,acc) schedule(static)
    do b = 1, nb
      i0 = int((int(na,8)*(b-1))/nb, int32) + 1_int32
      i1 = int((int(na,8)*b)/nb, int32)
      acc = 0_int32
      do i = i0, i1
        acc = acc + a(i)
      end do
      bsum(b) = acc
    end do
    !$omp end parallel do

    acc = 1_int32
    do b = 1, nb
      t = bsum(b)
      bsum(b) = acc
      acc = acc + t
    end do

    !$omp parallel do private(b,i,i0,i1,acc,t) schedule(static)
    do b = 1, nb
      i0 = int((int(na,8)*(b-1))/nb, int32) + 1_int32
      i1 = int((int(na,8)*b)/nb, int32)
      acc = bsum(b)
      do i = i0, i1
        t = a(i)
        a(i) = acc
        acc = acc + t
      end do
    end do
    !$omp end parallel do

    deallocate(bsum)
  end subroutine exclusive_prefix_sum

  logical function check_particles_are_sorted(part) result(ok)
    !=============================================================
    ! Return .true. if cell_id(:) is nondecreasing over active
    ! particles.
    !=============================================================
    type(ParticleSet), intent(in) :: part
    integer(int32) :: i

    ok = .true.

    if (.not. allocated(part%pv)) return
    if (part%n <= 1_int32) return

    if (.not. allocated(part%cell_id)) then
      ok = .false.
      return
    end if

    do i = 2, part%n
      if (part%cell_id(i) < part%cell_id(i-1)) then
        ok = .false.
        return
      end if
    end do
  end function check_particles_are_sorted


  logical function check_cell_indexing(part, n) result(ok)
    !=============================================================
    ! Strong consistency check for the Plist-equivalent metadata:
    !
    !   cell_id(i)
    !   cell_count(ic)
    !   cell_start(ic)
    !
    ! Checks:
    !   1) all cell ids lie in [1,ncells]
    !   2) cell_count matches a direct scan over cell_id
    !   3) cell_start/cell_count define valid contiguous ranges
    !   4) each particle in a cell range really belongs to that cell
    !=============================================================
    type(ParticleSet), intent(in) :: part
    integer(int32),     intent(in) :: n(3)

    integer(int32) :: ncells
    integer(int32) :: icell, i, i0, i1
    integer(int32), allocatable :: counts_scan(:)

    ok = .true.

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    if (.not. allocated(part%cell_id)) then
      ok = .false.
      return
    end if
    if (.not. allocated(part%cell_start)) then
      ok = .false.
      return
    end if
    if (.not. allocated(part%cell_count)) then
      ok = .false.
      return
    end if

    ncells = n(1) * n(2) * n(3)

    if (size(part%cell_count) /= ncells) then
      ok = .false.
      return
    end if

    if (size(part%cell_start) /= ncells) then
      ok = .false.
      return
    end if

    allocate(counts_scan(ncells))
    counts_scan = 0_int32

    do i = 1, part%n
      if (part%cell_id(i) < 1_int32 .or. part%cell_id(i) > ncells) then
        ok = .false.
        deallocate(counts_scan)
        return
      end if
      counts_scan(part%cell_id(i)) = counts_scan(part%cell_id(i)) + 1_int32
    end do

    do icell = 1, ncells
      if (part%cell_count(icell) /= counts_scan(icell)) then
        ok = .false.
        deallocate(counts_scan)
        return
      end if
    end do

    do icell = 1, ncells
      if (part%cell_count(icell) > 0_int32) then
        i0 = part%cell_start(icell)
        i1 = i0 + part%cell_count(icell) - 1_int32

        if (i0 < 1_int32 .or. i1 > part%n) then
          ok = .false.
          deallocate(counts_scan)
          return
        end if

        do i = i0, i1
          if (part%cell_id(i) /= icell) then
            ok = .false.
            deallocate(counts_scan)
            return
          end if
        end do
      end if
    end do

    deallocate(counts_scan)
  end function check_cell_indexing


  subroutine get_cell_particle_range(part, icell, i0, i1, count)
    !=============================================================
    ! Helper for collisions:
    ! return the contiguous particle range for one cell.
    !
    ! If the cell is empty:
    !   count = 0
    !   i0 = 0
    !   i1 = -1
    !=============================================================
    type(ParticleSet), intent(in)  :: part
    integer(int32),     intent(in)  :: icell
    integer(int32),     intent(out) :: i0, i1, count

    count = 0_int32
    i0    = 0_int32
    i1    = -1_int32

    if (.not. allocated(part%cell_count)) return
    if (.not. allocated(part%cell_start)) return

    if (icell < 1_int32 .or. icell > size(part%cell_count)) return

    count = part%cell_count(icell)
    if (count <= 0_int32) return

    i0 = part%cell_start(icell)
    i1 = i0 + count - 1_int32
  end subroutine get_cell_particle_range

end module mod_particle_sorting