module mod_collisionsGwenael
  !===========================================================================
  ! mod_collisionsGwenael
  !
  ! Goal: reproduce the LEGACY collision algorithm (Src/collisions.f90,
  ! subroutines `collisions` + `collision_OMP`, G. Fubiani) on top of the
  ! new modular architecture (ParticleSet, cell-sorted neighbour lists,
  ! mod_RNG), instead of the simplified/approximate re-implementations that
  ! currently live in mod_MCCcollisions / mod_BMCCcollisions.
  !
  ! This module replaces BOTH of those: a projectile species loops over
  ! *all* of its reaction channels in a single per-particle pass, whether
  ! the target is an untracked background neutral (legacy "MC" branch,
  ! Maxwellian-sampled target velocity, fixed reference density) or a
  ! tracked charged species (legacy "DSMC" branch, real target particle
  ! found via a cell-neighbour search, real local density from the
  ! deposited density grid). A single null-collision draw per projectile
  ! decides whether *any* reaction occurs and, if so, which one - exactly
  ! as in collision_OMP.
  !
  ! Known, deliberate departures from collision_OMP:
  !   1) No 2D source/sink diagnostic maps (ss2D_xy/xz/yz): these are not
  !      part of the current architecture's collision call signature, and
  !      nothing downstream currently consumes them.
  !   2) No true multi-MPI-rank reduction of particle counts (the existing
  !      MCC/BMCC modules in this branch don't do this either): nu_max is
  !      computed from this rank's particle counts only.
  !   3) The legacy `indexx` ascending sort of collision frequencies before
  !      the cumulative null-collision walk is dropped: summing a fixed
  !      (unsorted) list of non-negative rates and comparing against a
  !      single uniform draw gives the exact same selection PROBABILITIES
  !      regardless of summation order, so the sort has no statistical
  !      effect - it only existed to ever so slightly reduce floating-point
  !      round-off accumulation. It is reproduced here as a plain ordered
  !      sum over p_ncol(ptype).
  !   4) sig_Eex(:,2) ("dE_heavy", extra energy to heavy dissociation
  !      products) is read but intentionally NOT applied: grep of
  !      Src/collisions.f90 shows it is also unused there (only
  !      sig_Eex(:,ind_Eth) is ever read). The current mod_BMCCcollisions /
  !      mod_collisionProducts DO apply a dE_heavy kick; that is a
  !      deviation from legacy that this module does not reproduce.
  !
  ! Faithfully reproduced legacy mechanics:
  !   - Per ptype: global nu_max(ptype) = sum_icol np_mx(ttype)*sigv_mx,
  !     clipped to nu_uplim(ptype); Pmax = nu_max*ns_coll*dt; total trial
  !     count Nc_tmp = INT(dNc) + stochastic round-up from a SINGLE
  !     ran2(iseed(1)) draw (same odd legacy detail: always thread/proc 1's
  !     stream). Every OMP thread (iproc) then runs the SAME Nc_tmp trials,
  !     using a corrected nu_max_OMP(iproc) = Nc_tmp/(n_total*ns_coll*dt)
  !     so the realised trial count matches exactly.
  !   - np_mx(ttype): for a tracked (charged) ttype, the running MAX of the
  !     deposited density grid np_red over the physical domain (legacy
  !     calc_rho.f90); for an untracked background neutral, ni0(ttype).
  !   - Per trial: incident particle ip drawn uniformly from the FROZEN
  !     particle count (np_tot is only applied at the very end, exactly
  !     like legacy's deferred `np_tot = np_tot + np_add`). For each
  !     reaction channel of ptype, if the target species has not yet been
  !     "visited" in this trial, draw/find ONE target (Maxwellian sample
  !     for background neutrals, or a real particle from a neighbouring
  !     cell - own cell, then +-x, +-y, +-z, searched across all iproc
  !     exactly like legacy's Plist/tproc walk - for tracked targets), and
  !     cache its velocity/relative-energy for every later channel sharing
  !     the same target species within the SAME trial (legacy
  !     flag_ttype/vr/Ekr/vx_cm caching).
  !   - Single ran2 draw decides both the null-collision test AND (by
  !     walking the same cumulative sum) which channel fires.
  !   - Byproduct generation: equal polar-angle spacing in the CM frame,
  !     random azimuth, then a single random CM-frame rotation shared by
  !     all byproducts (legacy NOTES (2)); reactant slots are reused by
  !     byproducts of the same species when possible, otherwise new
  !     particles are appended; un-reused reactants are marked
  !     `flag_dead` (NOT removed immediately - removal/compaction is the
  !     job of mod_particleBC, exactly as for legacy's "kill done in
  !     mover").
  !   - Charge-exchange channels (rt==4): velocity hand-off only, no
  !     particle creation/destruction.
  !===========================================================================

  use iso_fortran_env,          only: int32, int8, real64
  use mod_particles,            only: ParticleSet
  use mod_RNG,                  only: ran2, load_gauss
  use omp_lib,                  only: omp_get_max_threads
  use mod_collisionDiagnostics, only: count_rxn, init_debug_diagnostics, record_np_mx, &
                                       record_target_search, record_npt

  implicit none
  private
  public :: perform_collisions_gwenael

  ! Print np_mx every nsav_npmx collision calls (ns_coll=10 calls per step,
  ! so nsav_npmx=100 → prints every 1000 steps, matching nsav=1000 in legacy)
  integer(int32), save :: call_count_npmx = 0_int32
  integer(int32), parameter :: nsav_npmx = 100_int32

  real(real64), parameter :: QE_ABS = 1.602176634e-19_real64
  real(real64), parameter :: PI_R8  = 3.14159265358979323846_real64

contains

  subroutine perform_collisions_gwenael( &
      part, n, h, ntype_tracked, ntype_all, mass, Ti, Nm, p_ncol, sig_list, &
      col_info, sigv_mx, sig, sig_Er, sig_Eex, ni0, ns_coll, dt, nu_uplim, &
      iseed, mpi_rank, Pcoll, dom_volume, np_red, bcnd, &
      ix_plane, iy_plane, iz_plane, sour_xy, sour_xz, sour_yz, &
      sink_xy, sink_xz, sink_yz, mom_loss)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in)    :: n(3)
    real(real64),   intent(in)    :: h(3)
    integer(int32), intent(in)    :: ntype_tracked, ntype_all
    real(real64),   intent(in)    :: mass(:), Ti(:), Nm(:)
    integer(int32), intent(in)    :: p_ncol(:)
    integer(int32), intent(in)    :: sig_list(:,:), col_info(:,:)
    real(real64),   intent(in)    :: sigv_mx(:,:), sig(:,:), sig_Er(:), sig_Eex(:,:)
    real(real64),   intent(in)    :: ni0(:)
    integer(int32), intent(in)    :: ns_coll
    real(real64),   intent(in)    :: dt, nu_uplim(:)
    integer(int32), intent(inout) :: iseed(:)
    integer(int32), intent(in)    :: mpi_rank
    real(real64),   intent(inout) :: Pcoll(:,:)
    ! Momentum-conservation diagnostic - vector counterpart of Pcoll,
    ! shaped (3, ntype, nproc): mom_loss(:,ptype,iproc).
    real(real64),   intent(inout) :: mom_loss(:,:,:)
    real(real64),   intent(in)    :: dom_volume
    ! Deposited density grid (legacy np_red), ghosted like the field
    ! arrays: (0:n(1)+2, 0:n(2)+2, 0:n(3)+2, ntype_tracked).
    real(real64),   intent(in)    :: np_red(0:,0:,0:,:)
    ! Boundary-condition flag grid, same ghosting; -1 marks an interior
    ! (physical, non-wall/non-ghost) cell - legacy calc_rho.f90 only takes
    ! np_mx's running max over bcnd==-1 cells.
    integer(int32), intent(in)    :: bcnd(0:,0:,0:)
    ! Legacy ss2D_xy/xz/yz equivalent: reaction-driven particle
    ! production (sour)/destruction (sink) per plane, per tracked species,
    ! per iproc (see mod_state.f90 sour_avg_xy/etc for why per-iproc).
    integer(int32), intent(in)    :: ix_plane, iy_plane, iz_plane
    real(real64),   intent(inout) :: sour_xy(0:,0:,:,:), sour_xz(0:,0:,:,:), sour_yz(0:,0:,:,:)
    real(real64),   intent(inout) :: sink_xy(0:,0:,:,:), sink_xz(0:,0:,:,:), sink_yz(0:,0:,:,:)

    integer(int32) :: ptype, iproc, nproc
    integer(int32) :: npt_sig
    integer(int32) :: Nc_tmp
    real(real64)   :: nu_max, Pmax, sum_np_tot, dNc, rnd
    real(real64), allocatable :: np_mx(:)
    integer(int32), allocatable :: n_add(:,:)   ! deferred particle-creation counters (species,iproc)

    integer(int32) :: s

    ! Precomputed cross-thread cell-target lookup for find_charged_target -
    ! see this subroutine's call to precompute_cell_target_info below for why.
    integer(int32) :: ncells
    integer(int32), allocatable :: total_cell_count(:,:)  ! (ncells, ntype_tracked)
    integer(int32), allocatable :: first_proc(:,:)        ! (ncells, ntype_tracked)

    nproc   = int(size(part,2), int32)
    npt_sig = count_valid_pts(sig_Er)
    if (npt_sig < 2_int32) return

    allocate(np_mx(ntype_all))
    allocate(n_add(ntype_tracked, nproc))
    n_add = 0_int32

    ! find_charged_target used to re-sum part(ttype,:)%cell_count LIVE
    ! (cell_count_all_procs, an O(nproc) loop) on every single trial that
    ! needs a charged target, plus an identical O(nproc) "which proc holds
    ! it" search right after - and Nc_tmp (trials/thread) measured in the
    ! ~800K-900K range on the ITER case, x nproc threads x up to 7 cells
    ! checked per trial (home + 6 neighbour fallback) = order a BILLION
    ! redundant scalar sums per collision call. Every tracked target
    ! species' per-cell/per-proc particle counts are frozen for this whole
    ! call (any reactions are deferred via n_add, applied only after every
    ! ptype's trial loop finishes - see this module's header comment,
    ! point 4 in "Faithfully reproduced legacy mechanics"), so the sum is
    ! call-invariant and only needs computing ONCE, as one vectorized pass
    ! per (species, thread) instead of a scalar loop re-run per trial.
    ncells = int(n(1), int32) * int(n(2), int32) * int(n(3), int32)
    allocate(total_cell_count(ncells, ntype_tracked))
    allocate(first_proc(ncells, ntype_tracked))
    call precompute_cell_target_info(part, ntype_tracked, nproc, ncells, &
                                      total_cell_count, first_proc)

    call init_debug_diagnostics(ntype_tracked)

    ! --- legacy np_mx(ttype): running max of deposited density over
    !     INTERIOR cells only (bcnd==-1), legacy calc_rho.f90 loop bounds
    !     1:n(i)+1 - not the full 1:n(i) physical range; ghost/wall cells
    !     are excluded either way via the bcnd mask. Fixed (ni0) reference
    !     density for background neutrals.
    do s = 1_int32, ntype_all
      if (s <= ntype_tracked) then
        np_mx(s) = 0.0_real64
        block
          integer(int32) :: ix, iy, iz
          do iz = 1_int32, n(3)+1_int32
            do iy = 1_int32, n(2)+1_int32
              do ix = 1_int32, n(1)+1_int32
                if (bcnd(ix,iy,iz) == -1_int32) then
                  np_mx(s) = max(np_mx(s), np_red(ix,iy,iz,s))
                end if
              end do
            end do
          end do
        end block
        call record_np_mx(s, np_mx(s))
      else
        np_mx(s) = ni0(s)
      end if
    end do

    ! --- print np_mx every nsav steps for direct comparison with legacy ---
    if (mpi_rank == 0 .and. mod(call_count_npmx, nsav_npmx) == 0) then
      do s = 1_int32, ntype_all
        write(*,'(a,i2,a,es12.4)') ' MOD_NPMX ttype=', s, ' np_mx=', np_mx(s)
      end do
    end if
    call_count_npmx = call_count_npmx + 1_int32

    !-------------------------------------------------------------------
    ! Outer loop over projectile species (serial: Nc_tmp must be the
    ! same, single value seen by every OMP thread for this ptype, exactly
    ! like legacy's two-phase structure).
    !-------------------------------------------------------------------
    do ptype = 1_int32, ntype_tracked

      if (p_ncol(ptype) <= 0_int32)    cycle
      if (mass(ptype)   <= 0.0_real64) cycle

      ! --- nu_max(ptype), clipped to nu_uplim
      nu_max = 0.0_real64
      do s = 1_int32, p_ncol(ptype)
        block
          integer(int32) :: ind_col, n_re, ttype
          ind_col = sig_list(ptype, s)
          n_re    = col_info(ind_col, 1)
          ttype   = col_info(ind_col, 2 + n_re)
          if (ttype >= 1_int32 .and. ttype <= ntype_all) then
            nu_max = nu_max + np_mx(ttype) * sigv_mx(ptype, ind_col)
          end if
        end block
      end do

      if (ptype <= size(nu_uplim)) nu_max = min(nu_max, nu_uplim(ptype))
      if (nu_max <= 0.0_real64) cycle

      Pmax = nu_max * real(ns_coll, real64) * dt
      if (Pmax > 1.0_real64) then
        write(*,*) 'mod_collisionsGwenael: collision probability > 1 for ptype=', ptype
        write(*,*) 'Pmax=', Pmax, ' - clipping to 1 (legacy would abort here)'
        Pmax = 1.0_real64
      end if

      ! --- total particle count for ptype, this rank (no MPI reduction)
      sum_np_tot = 0.0_real64
      do iproc = 1_int32, nproc
        if (allocated(part(ptype,iproc)%pv)) sum_np_tot = sum_np_tot + real(part(ptype,iproc)%n, real64)
      end do

      ! --- per-thread trial quota Nc_tmp, legacy stochastic round-up
      !     (always drawn from iseed(1), matching the legacy oddity)
      dNc = sum_np_tot * Pmax / real(nproc, real64)
      Nc_tmp = int(dNc, int32)
      rnd = ran2(iseed(1))
      if (rnd <= (dNc - real(Nc_tmp, real64))) Nc_tmp = Nc_tmp + 1_int32

      if (Nc_tmp <= 0_int32) cycle

      !-----------------------------------------------------------------
      ! Per-particle collision determination, parallel over iproc.
      ! Every thread runs the SAME number of trials (Nc_tmp), each drawn
      ! from its OWN (frozen) particle count - matching legacy exactly.
      !-----------------------------------------------------------------
      !$omp parallel do default(shared) private(iproc) schedule(static)
      do iproc = 1_int32, nproc

        if (.not. allocated(part(ptype,iproc)%pv)) cycle

        block
          integer(int32) :: n_total_local
          real(real64)   :: nu_max_OMP_local

          n_total_local = part(ptype,iproc)%n
          if (n_total_local <= 0_int32) cycle

          nu_max_OMP_local = real(Nc_tmp, real64) / real(n_total_local, real64) / &
                              (real(ns_coll, real64) * dt)
          if (nu_max_OMP_local <= 0.0_real64) cycle

          call run_trials_for_thread( &
              part, n, h, ptype, iproc, nproc, n_total_local, Nc_tmp, &
              nu_max_OMP_local, ntype_tracked, ntype_all, mass, Ti, Nm, &
              p_ncol, sig_list, col_info, sig, sig_Er, sig_Eex, ni0, &
              np_red, np_mx, iseed(iproc), Pcoll, n_add, &
              ix_plane, iy_plane, iz_plane, sour_xy, sour_xz, sour_yz, &
              sink_xy, sink_xz, sink_yz, mom_loss, &
              ncells, total_cell_count, first_proc)
        end block

      end do
      !$omp end parallel do

    end do ! ptype

    ! --- apply deferred particle creation (legacy np_tot = np_tot + np_add)
    do ptype = 1_int32, ntype_tracked
      do iproc = 1_int32, nproc
        if (n_add(ptype,iproc) > 0_int32) then
          part(ptype,iproc)%n = part(ptype,iproc)%n + n_add(ptype,iproc)
        end if
      end do
    end do

  end subroutine perform_collisions_gwenael


  !=========================================================================
  ! Run the Nc_tmp trials owned by one (ptype,iproc) thread.
  ! Mirrors legacy collision_OMP's `do ic=1,Nc_tmp` loop body.
  !=========================================================================
  subroutine run_trials_for_thread( &
      part, n, h, ptype, iproc, nproc, n_total, Nc_tmp, nu_max_OMP, &
      ntype_tracked, ntype_all, mass, Ti, Nm, p_ncol, sig_list, col_info, &
      sig, sig_Er, sig_Eex, ni0, np_red, np_mx, iseed_local, Pcoll, n_add, &
      ix_plane, iy_plane, iz_plane, sour_xy, sour_xz, sour_yz, &
      sink_xy, sink_xz, sink_yz, mom_loss, &
      ncells, total_cell_count, first_proc)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in)    :: n(3)
    real(real64),   intent(in)    :: h(3)
    integer(int32), intent(in)    :: ptype, iproc, nproc, n_total, Nc_tmp
    real(real64),   intent(in)    :: nu_max_OMP
    integer(int32), intent(in)    :: ntype_tracked, ntype_all
    real(real64),   intent(in)    :: mass(:), Ti(:), Nm(:)
    integer(int32), intent(in)    :: p_ncol(:)
    integer(int32), intent(in)    :: sig_list(:,:), col_info(:,:)
    real(real64),   intent(in)    :: sig(:,:), sig_Er(:), sig_Eex(:,:)
    real(real64),   intent(in)    :: ni0(:)
    real(real64),   intent(in)    :: np_red(0:,0:,0:,:)
    real(real64),   intent(in)    :: np_mx(:)
    integer(int32), intent(inout) :: iseed_local
    real(real64),   intent(inout) :: Pcoll(:,:)
    real(real64),   intent(inout) :: mom_loss(:,:,:)
    integer(int32), intent(inout) :: n_add(:,:)
    integer(int32), intent(in)    :: ix_plane, iy_plane, iz_plane
    real(real64),   intent(inout) :: sour_xy(0:,0:,:,:), sour_xz(0:,0:,:,:), sour_yz(0:,0:,:,:)
    real(real64),   intent(inout) :: sink_xy(0:,0:,:,:), sink_xz(0:,0:,:,:), sink_yz(0:,0:,:,:)
    ! Precomputed cross-thread cell-target lookup - see
    ! perform_collisions_gwenael's precompute_cell_target_info call.
    integer(int32), intent(in)    :: ncells
    integer(int32), intent(in)    :: total_cell_count(ncells,ntype_tracked)
    integer(int32), intent(in)    :: first_proc(ncells,ntype_tracked)

    integer(int32) :: ic, ip, icol, ind_col, n_re, ttype, c_ind
    integer(int32) :: npt_sig
    real(real64)   :: vx1, vy1, vz1
    real(real64)   :: rnd1, sum_nu, cumul

    ! Per-trial, per-target-species cache (indexed 1:ntype_all), mirroring
    ! legacy's flag_ttype / vr / Ekr / mu / vx_cm / sav_it / sav_tproc.
    logical        :: visited(ntype_all)
    logical        :: search_failed(ntype_all)
    real(real64)   :: cache_vr(ntype_all), cache_Ekr(ntype_all), cache_mu(ntype_all)
    real(real64)   :: cache_vxcm(ntype_all), cache_vycm(ntype_all), cache_vzcm(ntype_all)
    real(real64)   :: cache_tvx(ntype_all), cache_tvy(ntype_all), cache_tvz(ntype_all)
    integer(int32) :: cache_it(ntype_all), cache_itproc(ntype_all)

    real(real64)   :: nu_store(size(sig_list,2))
    real(real64)   :: vz_sav
    logical        :: dbg_sample

    npt_sig = count_valid_pts(sig_Er)
    vz_sav  = 0.0_real64   ! legacy: reset once per (ptype,iproc) call, persists across ic

    do ic = 1_int32, Nc_tmp

      ! Debug instrumentation (mod_collisionDiagnostics) is expensive
      ! (atomic increments from every OMP thread) - subsample 1-in-1000
      ! trials using the trial counter itself, never an extra ran2() draw,
      ! so the physics RNG stream is completely unaffected by whether
      ! instrumentation is active.
      dbg_sample = (mod(ic, 1000_int32) == 0_int32)

      ! --- incident particle, drawn from the frozen count
      rnd1 = ran2(iseed_local)
      ip = int(real(n_total, real64) * rnd1, int32) + 1_int32
      if (ip > n_total) ip = n_total

      vx1 = part(ptype,iproc)%pv(4,ip)
      vy1 = part(ptype,iproc)%pv(5,ip)
      vz1 = part(ptype,iproc)%pv(6,ip)

      visited       = .false.
      search_failed = .false.
      sum_nu        = 0.0_real64
      nu_store      = 0.0_real64

      do icol = 1_int32, p_ncol(ptype)

        ind_col = sig_list(ptype, icol)
        n_re    = col_info(ind_col, 1)
        ttype   = col_info(ind_col, 2 + n_re)

        if (ttype < 1_int32 .or. ttype > ntype_all) cycle
        if (search_failed(ttype)) cycle

        if (.not. visited(ttype)) then

          if (ttype > ntype_tracked) then
            ! ----- background neutral target (legacy "MC" branch) -----
            call sample_neutral_target( &
                ttype, mass, Ti, iseed_local, vz_sav, &
                cache_tvx(ttype), cache_tvy(ttype), cache_tvz(ttype))
            cache_it(ttype)     = 0_int32
            cache_itproc(ttype) = 0_int32
          else
            ! ----- tracked charged target (legacy "DSMC" branch) -----
            call find_charged_target( &
                part, ttype, nproc, n, h, &
                part(ptype,iproc)%pv(1,ip), part(ptype,iproc)%pv(2,ip), part(ptype,iproc)%pv(3,ip), &
                iseed_local, cache_it(ttype), cache_itproc(ttype), &
                ncells, total_cell_count(:,ttype), first_proc(:,ttype))

            if (cache_it(ttype) <= 0_int32) then
              search_failed(ttype) = .true.
              if (dbg_sample) call record_target_search(ttype, .false.)
              cycle
            end if
            if (dbg_sample) call record_target_search(ttype, .true.)

            cache_tvx(ttype) = part(ttype, cache_itproc(ttype))%pv(4,cache_it(ttype))
            cache_tvy(ttype) = part(ttype, cache_itproc(ttype))%pv(5,cache_it(ttype))
            cache_tvz(ttype) = part(ttype, cache_itproc(ttype))%pv(6,cache_it(ttype))
          end if

          cache_mu(ttype) = abs(mass(ptype))*abs(mass(ttype)) / &
                            (abs(mass(ptype)) + abs(mass(ttype)))

          cache_vr(ttype) = sqrt( (vx1-cache_tvx(ttype))**2 + &
                                  (vy1-cache_tvy(ttype))**2 + &
                                  (vz1-cache_tvz(ttype))**2 )

          cache_Ekr(ttype) = 0.5_real64*cache_mu(ttype)*cache_vr(ttype)**2 / QE_ABS

          cache_vxcm(ttype) = (abs(mass(ptype))*vx1 + abs(mass(ttype))*cache_tvx(ttype)) / &
                              (abs(mass(ptype)) + abs(mass(ttype)))
          cache_vycm(ttype) = (abs(mass(ptype))*vy1 + abs(mass(ttype))*cache_tvy(ttype)) / &
                              (abs(mass(ptype)) + abs(mass(ttype)))
          cache_vzcm(ttype) = (abs(mass(ptype))*vz1 + abs(mass(ttype))*cache_tvz(ttype)) / &
                              (abs(mass(ptype)) + abs(mass(ttype)))

          visited(ttype) = .true.

        end if

        block
          real(real64) :: sig_p, np_t

          sig_p = interp_sigma(cache_Ekr(ttype), sig, ind_col, sig_Er, npt_sig)

          if (ttype > ntype_tracked) then
            np_t = ni0(ttype)
          else
            np_t = local_density_8pt(np_red, ttype, &
                       part(ptype,iproc)%pv(1,ip), part(ptype,iproc)%pv(2,ip), &
                       part(ptype,iproc)%pv(3,ip), n, h)
            if (dbg_sample) call record_npt(ttype, np_t, np_mx(ttype))
          end if

          nu_store(icol) = np_t * sig_p * cache_vr(ttype) / nu_max_OMP
          sum_nu = sum_nu + nu_store(icol)
        end block

      end do ! icol

      ! --- single draw: null-collision test AND channel selection
      rnd1 = ran2(iseed_local)
      if (rnd1 > sum_nu) cycle  ! null collision

      cumul  = 0.0_real64
      c_ind  = 0_int32
      ttype  = 0_int32
      do icol = 1_int32, p_ncol(ptype)
        if (nu_store(icol) <= 0.0_real64) cycle
        cumul = cumul + nu_store(icol)
        if (rnd1 <= cumul) then
          ind_col = sig_list(ptype, icol)
          n_re    = col_info(ind_col, 1)
          ttype   = col_info(ind_col, 2 + n_re)
          c_ind   = ind_col
          exit
        end if
      end do

      if (c_ind <= 0_int32) cycle

      call count_rxn(c_ind)

      call apply_reaction_products( &
          part, n, h, n_add, ptype, iproc, ip, ttype, &
          cache_tvx(ttype), cache_tvy(ttype), cache_tvz(ttype), &
          cache_it(ttype), cache_itproc(ttype), ntype_tracked, &
          c_ind, col_info, sig_Eex, mass, Nm, Pcoll, iseed_local, &
          ix_plane, iy_plane, iz_plane, sour_xy, sour_xz, sour_yz, &
          sink_xy, sink_xz, sink_yz, mom_loss)

    end do ! ic

  end subroutine run_trials_for_thread


  !=========================================================================
  ! Apply the products of the chosen reaction c_ind.
  ! Mirrors legacy collision_OMP lines ~745-1003 (all-but-CEX + CEX).
  !=========================================================================
  subroutine apply_reaction_products( &
      part, n, h, n_add, ptype, iproc, ip, ttype, tvx, tvy, tvz, it, itproc, &
      ntype_tracked, c_ind, col_info, sig_Eex, mass, Nm, Pcoll, iseed_local, &
      ix_plane, iy_plane, iz_plane, sour_xy, sour_xz, sour_yz, &
      sink_xy, sink_xz, sink_yz, mom_loss)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in)    :: n(3)
    real(real64),   intent(in)    :: h(3)
    integer(int32), intent(inout) :: n_add(:,:)
    integer(int32), intent(in)    :: ptype, iproc, ip, ttype
    real(real64),   intent(in)    :: tvx, tvy, tvz
    integer(int32), intent(in)    :: it, itproc
    integer(int32), intent(in)    :: ntype_tracked
    integer(int32), intent(in)    :: c_ind
    integer(int32), intent(in)    :: col_info(:,:)
    real(real64),   intent(in)    :: sig_Eex(:,:), mass(:), Nm(:)
    real(real64),   intent(inout) :: Pcoll(:,:)
    integer(int32), intent(inout) :: iseed_local
    integer(int32), intent(in)    :: ix_plane, iy_plane, iz_plane
    real(real64),   intent(inout) :: sour_xy(0:,0:,:,:), sour_xz(0:,0:,:,:), sour_yz(0:,0:,:,:)
    real(real64),   intent(inout) :: sink_xy(0:,0:,:,:), sink_xz(0:,0:,:,:), sink_yz(0:,0:,:,:)
    ! Momentum-conservation diagnostic - vector counterpart of Pcoll,
    ! shaped (3, ntype, nproc): mom_loss(:,ptype,iproc).
    real(real64),   intent(inout) :: mom_loss(:,:,:)

    integer(int32) :: n_re, n_by, rt, i_by, btype, ib, ibproc
    integer(int32) :: flag_proj_kept, flag_targ_kept
    real(real64)   :: vx1, vy1, vz1
    real(real64)   :: mu, vr, sum_mass_inv, Eth, Erel_after
    real(real64)   :: vx_cm, vy_cm, vz_cm
    real(real64)   :: th, th_add, phi, cos_th, sin_th, cos_th_s, phi_s
    real(real64)   :: ex, ey, ez, ex1, ey1, ez1
    real(real64)   :: rnd1, rnd2, vp, v2old, v2new
    real(real64)   :: bvx_old, bvy_old, bvz_old

    n_re = col_info(c_ind, 1)
    n_by = col_info(c_ind, 2)
    rt   = col_info(c_ind, 2 + n_re + n_by + 1)

    !-----------------------------------------------------------------
    ! Charge exchange: velocity hand-off only (legacy lines 953-1003).
    ! No particle creation/destruction.
    !-----------------------------------------------------------------
    if (rt == 4_int32) then

      ! Projectile takes on the target's (sampled or real) velocity.
      vx1 = part(ptype,iproc)%pv(4,ip)
      vy1 = part(ptype,iproc)%pv(5,ip)
      vz1 = part(ptype,iproc)%pv(6,ip)
      v2old = vx1*vx1 + vy1*vy1 + vz1*vz1

      part(ptype,iproc)%pv(4,ip) = tvx
      part(ptype,iproc)%pv(5,ip) = tvy
      part(ptype,iproc)%pv(6,ip) = tvz

      v2new = tvx*tvx + tvy*tvy + tvz*tvz

      Pcoll(ptype,iproc) = Pcoll(ptype,iproc) + &
          0.5_real64*Nm(ptype)*abs(mass(ptype))*(v2new - v2old)

      mom_loss(1,ptype,iproc) = mom_loss(1,ptype,iproc) + Nm(ptype)*abs(mass(ptype))*(tvx - vx1)
      mom_loss(2,ptype,iproc) = mom_loss(2,ptype,iproc) + Nm(ptype)*abs(mass(ptype))*(tvy - vy1)
      mom_loss(3,ptype,iproc) = mom_loss(3,ptype,iproc) + Nm(ptype)*abs(mass(ptype))*(tvz - vz1)

      ! If the target was itself a tracked charged particle, it would, in
      ! the general legacy formulation, receive the projectile's old
      ! velocity. None of the present chemistry exercises CEX between two
      ! tracked species (all CEX reactions here target a background
      ! neutral), so this branch is intentionally not exercised; tracked
      ! ttype CEX targets are left untouched rather than guessed at.
      return

    end if

    !-----------------------------------------------------------------
    ! All other reaction types (elastic, excitation, ionization,
    ! dissociation, recombination, ...).
    !-----------------------------------------------------------------
    vx1 = part(ptype,iproc)%pv(4,ip)
    vy1 = part(ptype,iproc)%pv(5,ip)
    vz1 = part(ptype,iproc)%pv(6,ip)

    mu = abs(mass(ptype))*abs(mass(ttype)) / (abs(mass(ptype)) + abs(mass(ttype)))
    vr = sqrt((vx1-tvx)**2 + (vy1-tvy)**2 + (vz1-tvz)**2)

    vx_cm = (abs(mass(ptype))*vx1 + abs(mass(ttype))*tvx) / (abs(mass(ptype)) + abs(mass(ttype)))
    vy_cm = (abs(mass(ptype))*vy1 + abs(mass(ttype))*tvy) / (abs(mass(ptype)) + abs(mass(ttype)))
    vz_cm = (abs(mass(ptype))*vz1 + abs(mass(ttype))*tvz) / (abs(mass(ptype)) + abs(mass(ttype)))

    Eth = sig_Eex(c_ind, 1)   ! sig_Eex(:,2) "dE_heavy" intentionally unused - see module header

    Erel_after = 0.5_real64*mu*vr*vr - Eth*QE_ABS
    if (Erel_after < 0.0_real64) return   ! legacy aborts; we treat as a failed/null event

    sum_mass_inv = 0.0_real64
    do i_by = 1_int32, n_by
      btype = col_info(c_ind, 2+n_re+i_by)
      if (mass(btype) > 0.0_real64) sum_mass_inv = sum_mass_inv + 1.0_real64/abs(mass(btype))
    end do
    if (sum_mass_inv <= 0.0_real64) return

    rnd1 = ran2(iseed_local); rnd2 = ran2(iseed_local)
    th_add = 2.0_real64*PI_R8/real(n_by, real64)
    cos_th = 1.0_real64 - 2.0_real64*rnd1
    th     = acos(cos_th)
    phi    = 2.0_real64*PI_R8*rnd2

    rnd1 = ran2(iseed_local); rnd2 = ran2(iseed_local)
    cos_th_s = 1.0_real64 - 2.0_real64*rnd1
    phi_s    = 2.0_real64*PI_R8*rnd2

    flag_proj_kept = 0_int32
    flag_targ_kept = 0_int32

    do i_by = 1_int32, n_by

      btype = col_info(c_ind, 2+n_re+i_by)

      cos_th = cos(th)
      sin_th = sin(th)
      th = th + th_add

      ex = cos_th
      ey = sin_th*sin(phi)
      ez = sin_th*cos(phi)
      call scatter_legacy(ex1, ey1, ez1, ex, ey, ez, cos_th_s, phi_s)

      if (mass(btype) <= 0.0_real64) cycle   ! untracked/negative-mass product: discard

      ! Elastic/excitation: target keeps its identity in place (legacy
      ! `if (rt==1 .or. rt==3) ... ib=ip / ib=it`); for other reaction
      ! types the first matching reactant slot is reused, subsequent
      ! same-species byproducts are appended as new particles.
      if (btype == ptype .and. flag_proj_kept == 0_int32) then
        ib = ip; ibproc = iproc
        flag_proj_kept = 1_int32
      else if (btype == ttype .and. ttype <= ntype_tracked .and. flag_targ_kept == 0_int32) then
        ib = it; ibproc = itproc
        flag_targ_kept = 1_int32
      else if (btype <= ntype_tracked) then
        call append_deferred(part, n_add, btype, iproc, &
                              part(ptype,iproc)%pv(1,ip), part(ptype,iproc)%pv(2,ip), &
                              part(ptype,iproc)%pv(3,ip), ib)
        ibproc = iproc
      else
        cycle   ! btype > ntype_tracked: untracked background-neutral byproduct, no ParticleSet to write into
      end if

      vp = sqrt(2.0_real64*(Erel_after/sum_mass_inv) / abs(mass(btype))**2)

      bvx_old = part(btype,ibproc)%pv(4,ib)
      bvy_old = part(btype,ibproc)%pv(5,ib)
      bvz_old = part(btype,ibproc)%pv(6,ib)
      v2old = bvx_old*bvx_old + bvy_old*bvy_old + bvz_old*bvz_old

      part(btype,ibproc)%pv(4,ib) = vx_cm + vp*ex1
      part(btype,ibproc)%pv(5,ib) = vy_cm + vp*ey1
      part(btype,ibproc)%pv(6,ib) = vz_cm + vp*ez1

      v2new = part(btype,ibproc)%pv(4,ib)**2 + part(btype,ibproc)%pv(5,ib)**2 + &
              part(btype,ibproc)%pv(6,ib)**2

      Pcoll(btype,ibproc) = Pcoll(btype,ibproc) + &
          0.5_real64*Nm(btype)*abs(mass(btype))*(v2new - v2old)

      mom_loss(1,btype,ibproc) = mom_loss(1,btype,ibproc) + &
          Nm(btype)*abs(mass(btype))*(part(btype,ibproc)%pv(4,ib) - bvx_old)
      mom_loss(2,btype,ibproc) = mom_loss(2,btype,ibproc) + &
          Nm(btype)*abs(mass(btype))*(part(btype,ibproc)%pv(5,ib) - bvy_old)
      mom_loss(3,btype,ibproc) = mom_loss(3,btype,ibproc) + &
          Nm(btype)*abs(mass(btype))*(part(btype,ibproc)%pv(6,ib) - bvz_old)

      ! Legacy ss2D(1,...) equivalent: this reaction produced one btype
      ! macroparticle here (whether ib is a freshly-appended slot or a
      ! reused reactant slot doesn't change the physical production rate).
      call deposit_plane_event( &
          part(btype,ibproc)%pv(1,ib), part(btype,ibproc)%pv(2,ib), part(btype,ibproc)%pv(3,ib), &
          n, h, ix_plane, iy_plane, iz_plane, btype, ibproc, Nm(btype), &
          sour_xy, sour_xz, sour_yz)

    end do

    ! Reactants that did not survive as a byproduct slot are marked dead
    ! (compaction happens later, in mod_particleBC - exactly like legacy's
    ! "kill done in mover"). Legacy ss2D(2,...) equivalent: deposit the
    ! destruction at the reactant's own (unchanged) position.
    if (flag_proj_kept == 0_int32) then
      if (allocated(part(ptype,iproc)%flag_dead)) part(ptype,iproc)%flag_dead(ip) = 1_int8
      call deposit_plane_event( &
          part(ptype,iproc)%pv(1,ip), part(ptype,iproc)%pv(2,ip), part(ptype,iproc)%pv(3,ip), &
          n, h, ix_plane, iy_plane, iz_plane, ptype, iproc, Nm(ptype), &
          sink_xy, sink_xz, sink_yz)
    end if
    if (flag_targ_kept == 0_int32 .and. ttype <= ntype_tracked .and. it > 0_int32) then
      if (allocated(part(ttype,itproc)%flag_dead)) part(ttype,itproc)%flag_dead(it) = 1_int8
      call deposit_plane_event( &
          part(ttype,itproc)%pv(1,it), part(ttype,itproc)%pv(2,it), part(ttype,itproc)%pv(3,it), &
          n, h, ix_plane, iy_plane, iz_plane, ttype, itproc, Nm(ttype), &
          sink_xy, sink_xz, sink_yz)
    end if

  end subroutine apply_reaction_products


  !=========================================================================
  ! Deposit one reaction-driven production/destruction event into whichever
  ! of the three output planes it is nearest to. Nearest-cell (not the
  ! moments code's bilinear/plane-straddle interpolation) - this fires from
  ! inside the hot per-trial collision loop, so a cheaper deposition was
  ! chosen deliberately; it is not legacy's exact ki(1-4) bilinear formula.
  !=========================================================================
  subroutine deposit_plane_event(x, y, z, n, h, ix_plane, iy_plane, iz_plane, &
                                  species, iproc, weight, arr_xy, arr_xz, arr_yz)
    real(real64),   intent(in)    :: x, y, z
    integer(int32), intent(in)    :: n(3)
    real(real64),   intent(in)    :: h(3)
    integer(int32), intent(in)    :: ix_plane, iy_plane, iz_plane
    integer(int32), intent(in)    :: species, iproc
    real(real64),   intent(in)    :: weight
    real(real64),   intent(inout) :: arr_xy(0:,0:,:,:), arr_xz(0:,0:,:,:), arr_yz(0:,0:,:,:)

    integer(int32) :: ix, iy, iz

    ix = int(x / h(1), int32) + 1_int32
    iy = int(y / h(2), int32) + 1_int32
    iz = int(z / h(3), int32) + 1_int32

    if (iz == iz_plane .and. ix >= 0_int32 .and. ix <= n(1)+2_int32 .and. &
                              iy >= 0_int32 .and. iy <= n(2)+2_int32) then
      arr_xy(ix,iy,species,iproc) = arr_xy(ix,iy,species,iproc) + weight
    end if
    if (iy == iy_plane .and. ix >= 0_int32 .and. ix <= n(1)+2_int32 .and. &
                              iz >= 0_int32 .and. iz <= n(3)+2_int32) then
      arr_xz(ix,iz,species,iproc) = arr_xz(ix,iz,species,iproc) + weight
    end if
    if (ix == ix_plane .and. iy >= 0_int32 .and. iy <= n(2)+2_int32 .and. &
                              iz >= 0_int32 .and. iz <= n(3)+2_int32) then
      arr_yz(iy,iz,species,iproc) = arr_yz(iy,iz,species,iproc) + weight
    end if
  end subroutine deposit_plane_event


  !=========================================================================
  ! Maxwellian sample of a background-neutral target velocity, replicating
  ! legacy's vz_sav Box-Muller-pair caching trick exactly (load_gauss
  ! produces two Gaussian numbers per call; the second one is cached and
  ! reused - WITHOUT rescaling - as the z-component of the NEXT neutral
  ! sample drawn in this thread's trial loop, whatever species it is for).
  !=========================================================================
  subroutine sample_neutral_target(ttype, mass, Ti, iseed_local, vz_sav, vx2, vy2, vz2)
    integer(int32), intent(in)    :: ttype
    real(real64),   intent(in)    :: mass(:), Ti(:)
    integer(int32), intent(inout) :: iseed_local
    real(real64),   intent(inout) :: vz_sav
    real(real64),   intent(out)   :: vx2, vy2, vz2

    real(real64) :: vt0, rnd(2)

    vt0 = sqrt(2.0_real64*QE_ABS*Ti(ttype) / abs(mass(ttype)))

    rnd(1) = ran2(iseed_local); rnd(2) = ran2(iseed_local)
    call load_gauss(vx2, vy2, vt0, rnd)

    if (vz_sav == 0.0_real64) then
      rnd(1) = ran2(iseed_local); rnd(2) = ran2(iseed_local)
      call load_gauss(vz2, vz_sav, vt0, rnd)
    else
      vz2 = vz_sav
      vz_sav = 0.0_real64
    end if
  end subroutine sample_neutral_target


  !=========================================================================
  ! Find a real target particle of tracked species ttype near position
  ! (xp,yp,zp): own cell first, then the 6 face-neighbours (+x,-x,+y,-y,
  ! +z,-z, in that order, stopping at the first non-empty ring) - exactly
  ! the legacy Plist search order - summed across all iproc, then a
  ! uniformly-random pick within whichever iproc actually holds particles
  ! in the resolved cell (legacy's tproc walk).
  !=========================================================================
  subroutine find_charged_target(part, ttype, nproc, n, h, xp, yp, zp, iseed_local, it, itproc, &
                                  ncells, total_cell_count, first_proc)
    type(ParticleSet), intent(in)    :: part(:,:)
    integer(int32),    intent(in)    :: ttype, nproc
    integer(int32),    intent(in)    :: n(3)
    real(real64),      intent(in)    :: h(3)
    real(real64),      intent(in)    :: xp, yp, zp
    integer(int32),    intent(inout) :: iseed_local
    integer(int32),    intent(out)   :: it, itproc
    ! Precomputed cross-thread cell-target lookup (this ttype's column of
    ! perform_collisions_gwenael's arrays) - replaces two separate O(nproc)
    ! live scans (cell_count_all_procs + the "which proc" search below)
    ! with O(1) lookups. See perform_collisions_gwenael's
    ! precompute_cell_target_info call for why this is safe (target
    ! particle counts are frozen for the whole collision call). nproc is
    ! no longer used in this body but kept in the signature to match the
    ! call site unchanged.
    integer(int32),    intent(in)    :: ncells
    integer(int32),    intent(in)    :: total_cell_count(ncells)
    integer(int32),    intent(in)    :: first_proc(ncells)

    integer(int32) :: ix, iy, iz, ici, ict, idir
    integer(int32) :: cnt_total, cnt_proc, target_proc
    real(real64)   :: rnd

    it     = 0_int32
    itproc = 0_int32

    ix = int(xp/h(1), int32) + 1_int32
    iy = int(yp/h(2), int32) + 1_int32
    iz = int(zp/h(3), int32) + 1_int32
    ix = max(1_int32, min(n(1), ix))
    iy = max(1_int32, min(n(2), iy))
    iz = max(1_int32, min(n(3), iz))

    ici = (ix-1_int32) + n(1)*((iy-1_int32) + n(2)*(iz-1_int32)) + 1_int32

    ict = ici
    cnt_total = total_cell_count(ict)

    if (cnt_total <= 0_int32) then
      do idir = 1_int32, 6_int32
        select case (idir)
        case (1_int32); if (ix <  n(1)) ict = ici + 1_int32
        case (2_int32); if (ix >  1_int32) ict = ici - 1_int32
        case (3_int32); if (iy <  n(2)) ict = ici + n(1)
        case (4_int32); if (iy >  1_int32) ict = ici - n(1)
        case (5_int32); if (iz <  n(3)) ict = ici + n(1)*n(2)
        case (6_int32); if (iz >  1_int32) ict = ici - n(1)*n(2)
        end select
        cnt_total = total_cell_count(ict)
        if (cnt_total > 0_int32) exit
      end do
    end if

    if (cnt_total <= 0_int32) return  ! search_failed for this ttype this trial

    ! Which iproc holds the particle - precomputed as the same "first proc
    ! in 1..nproc order with cell_count>0" legacy's cnt_proc loop picked.
    target_proc = first_proc(ict)
    if (target_proc <= 0_int32) return

    cnt_proc = part(ttype,target_proc)%cell_count(ict)
    rnd = ran2(iseed_local)
    it = part(ttype,target_proc)%cell_start(ict) + min(cnt_proc-1_int32, int(real(cnt_proc,real64)*rnd, int32))
    itproc = target_proc

  end subroutine find_charged_target


  subroutine precompute_cell_target_info(part, ntype_tracked, nproc, ncells, &
                                          total_cell_count, first_proc)
    ! Once-per-call replacement for what used to be cell_count_all_procs
    ! (an O(nproc) live re-sum) plus find_charged_target's own O(nproc)
    ! "which proc holds it" search, both called from inside the per-trial
    ! loop - see perform_collisions_gwenael's call site for the full
    ! rationale and the "frozen during this call" argument that makes this
    ! valid. Vectorized whole-array ops (Fortran array-section syntax), not
    ! a scalar ncells-deep loop, so the compiler can SIMD this - the total
    ! element count (ntype_tracked * nproc * ncells) is itself large for
    ! ITER-scale grids, but each element op is a cheap add/compare, and
    ! this only runs ONCE per collision call versus the up to ~10^9 scalar,
    ! branchy, function-call-per-lookup ops the live version needed.
    type(ParticleSet), intent(in)  :: part(:,:)
    integer(int32),    intent(in)  :: ntype_tracked, nproc, ncells
    integer(int32),    intent(out) :: total_cell_count(ncells,ntype_tracked)
    integer(int32),    intent(out) :: first_proc(ncells,ntype_tracked)

    integer(int32) :: s, jproc, nblk, b, b0, b1, i0, i1
    integer(int32) :: clo(ntype_tracked,nproc), chi(ntype_tracked,nproc)

    ! Each set's cell_count is zero outside cell_lo..cell_hi (see
    ! ParticleSet); after the per-species redistribution that range is one
    ! z-slab. Reading only inside it - and running over cells in parallel -
    ! replaces a serial ntype*nproc*ncells whole-array pass. Same result:
    ! same sums, and first_proc is still the lowest jproc with a particle.
    do s = 1_int32, ntype_tracked
      do jproc = 1_int32, nproc
        clo(s,jproc) = 1_int32
        chi(s,jproc) = 0_int32
        if (.not. allocated(part(s,jproc)%pv)) cycle
        if (.not. allocated(part(s,jproc)%cell_count)) cycle
        if (size(part(s,jproc)%cell_count) < ncells) cycle
        clo(s,jproc) = max(1_int32, part(s,jproc)%cell_lo)
        chi(s,jproc) = min(ncells, part(s,jproc)%cell_hi)
      end do
    end do

    ! Parallel over contiguous cell blocks, each adding every set's
    ! overlapping slice in jproc order (contiguous array sections). A
    ! per-cell loop over jproc here would be outer-loop-vectorized by ifx
    ! 2025.0 into an AVX-512 gather that faults SIGILL - see
    ! redistribute_species_by_cell.
    nblk = max(1_int32, min(ncells, 8_int32*int(omp_get_max_threads(), int32)))
    !$omp parallel do private(b,b0,b1,s,jproc,i0,i1) schedule(static)
    do b = 1_int32, nblk
      b0 = int((int(ncells,8)*(b-1))/nblk, int32) + 1_int32
      b1 = int((int(ncells,8)*b)/nblk, int32)
      do s = 1_int32, ntype_tracked
        total_cell_count(b0:b1,s) = 0_int32
        first_proc(b0:b1,s)       = 0_int32
        do jproc = 1_int32, nproc
          i0 = max(b0, clo(s,jproc)); i1 = min(b1, chi(s,jproc))
          if (i0 > i1) cycle
          total_cell_count(i0:i1,s) = total_cell_count(i0:i1,s) + part(s,jproc)%cell_count(i0:i1)
          where (first_proc(i0:i1,s) == 0_int32 .and. part(s,jproc)%cell_count(i0:i1) > 0_int32)
            first_proc(i0:i1,s) = jproc
          end where
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine precompute_cell_target_info


  !=========================================================================
  ! Legacy's "DSMC" local density: an UNWEIGHTED average of the deposited
  ! density at the 8 grid points surrounding the projectile (legacy
  ! collision_OMP lines ~530-533) - not a true trilinear interpolation.
  !=========================================================================
  real(real64) function local_density_8pt(np_red, ttype, xp, yp, zp, n, h) result(np_t)
    real(real64),   intent(in) :: np_red(0:,0:,0:,:)
    integer(int32), intent(in) :: ttype
    real(real64),   intent(in) :: xp, yp, zp
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)

    integer(int32) :: ix, iy, iz

    ix = int(xp/h(1), int32) + 1_int32
    iy = int(yp/h(2), int32) + 1_int32
    iz = int(zp/h(3), int32) + 1_int32
    ix = max(1_int32, min(n(1), ix))
    iy = max(1_int32, min(n(2), iy))
    iz = max(1_int32, min(n(3), iz))

    np_t = 0.125_real64*( &
        np_red(ix,  iy,  iz,  ttype) + np_red(ix+1,iy,  iz,  ttype) + &
        np_red(ix,  iy+1,iz,  ttype) + np_red(ix+1,iy+1,iz,  ttype) + &
        np_red(ix,  iy,  iz+1,ttype) + np_red(ix+1,iy,  iz+1,ttype) + &
        np_red(ix,  iy+1,iz+1,ttype) + np_red(ix+1,iy+1,iz+1,ttype) )
  end function local_density_8pt


  !=========================================================================
  ! Append a new particle of species btype (at the position of the
  ! triggering projectile) into a DEFERRED staging slot: the array is
  ! grown if needed, but part(btype,iproc)%n is NOT bumped here - that
  ! happens once, for every species, after all collisions for this step
  ! are done (legacy's deferred np_tot = np_tot + np_add).
  !=========================================================================
  subroutine append_deferred(part, n_add, btype, iproc, xp, yp, zp, ib)
    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32),    intent(inout) :: n_add(:,:)
    integer(int32),    intent(in)    :: btype, iproc
    real(real64),      intent(in)    :: xp, yp, zp
    integer(int32),    intent(out)   :: ib

    ib = part(btype,iproc)%n + n_add(btype,iproc) + 1_int32

    if (ib > part(btype,iproc)%nmax) then
      call part(btype,iproc)%ensure_capacity(ib)
    end if

    part(btype,iproc)%pv(1,ib) = xp
    part(btype,iproc)%pv(2,ib) = yp
    part(btype,iproc)%pv(3,ib) = zp
    part(btype,iproc)%pv(4,ib) = 0.0_real64
    part(btype,iproc)%pv(5,ib) = 0.0_real64
    part(btype,iproc)%pv(6,ib) = 0.0_real64
    if (allocated(part(btype,iproc)%flag_dead)) part(btype,iproc)%flag_dead(ib) = 0_int8
    if (allocated(part(btype,iproc)%flag_cex))  part(btype,iproc)%flag_cex(ib)  = 0_int32

    n_add(btype,iproc) = n_add(btype,iproc) + 1_int32
  end subroutine append_deferred


  !=========================================================================
  ! Rotate vector (vx,vy,vz) by polar angle costheta and azimuthal angle
  ! phi - verbatim port of legacy's `scatter` subroutine.
  !=========================================================================
  subroutine scatter_legacy(vx1, vy1, vz1, vx, vy, vz, costheta, phi)
    real(real64), intent(out) :: vx1, vy1, vz1
    real(real64), intent(in)  :: vx, vy, vz, costheta, phi

    real(real64) :: sintheta, cosphi, sinphi, v, vv

    sintheta = sqrt(max(1.0_real64 - costheta**2, 0.0_real64))
    sinphi = sin(phi)
    cosphi = cos(phi)
    v = sqrt(vx**2 + vy**2 + vz**2)

    if (v == 0.0_real64) then
      vx1 = 0.0_real64; vy1 = 0.0_real64; vz1 = 0.0_real64
      return
    end if

    if (abs(vy) > abs(vz)) then
      vv  = sqrt(vx**2 + vy**2)
      vx1 = vx*costheta + (vy*v*sinphi + vx*vz*cosphi)/vv*sintheta
      vy1 = vy*costheta + (-vx*v*sinphi + vy*vz*cosphi)/vv*sintheta
      vz1 = vz*costheta - vv*cosphi*sintheta
    else
      vv  = sqrt(vx**2 + vz**2)
      vx1 = vx*costheta + (vz*v*sinphi - vy*vx*cosphi)/vv*sintheta
      vy1 = vy*costheta + vv*cosphi*sintheta
      vz1 = vz*costheta - (vx*v*sinphi + vy*vz*cosphi)/vv*sintheta
    end if
  end subroutine scatter_legacy


  !=========================================================================
  ! Bisection + linear interpolation of cross section sig(:,ind_col) on
  ! the shared energy grid sig_Er - equivalent to legacy's 5-step
  ! dichotomy followed by a local linear scan (same end result, fewer
  ! lines); clamped at the table edges like the legacy code's implicit
  ! behaviour at ipt_L=1 / ipt_R=npt.
  !=========================================================================
  real(real64) function interp_sigma(E, sig, ind_col, sig_Er, npt) result(sig_p)
    real(real64),   intent(in) :: E, sig(:,:), sig_Er(:)
    integer(int32), intent(in) :: ind_col, npt
    integer(int32) :: lo, hi, mid

    if (npt < 2_int32) then; sig_p = 0.0_real64; return; end if
    if (E <= sig_Er(1))   then; sig_p = sig(1,ind_col);   return; end if
    if (E >= sig_Er(npt)) then; sig_p = sig(npt,ind_col); return; end if

    lo = 1_int32; hi = npt
    do while (hi - lo > 1_int32)
      mid = (lo + hi) / 2_int32
      if (E <= sig_Er(mid)) then; hi = mid; else; lo = mid; end if
    end do

    sig_p = sig(lo,ind_col) + &
        (E - sig_Er(lo))*(sig(hi,ind_col) - sig(lo,ind_col)) / (sig_Er(hi) - sig_Er(lo))
  end function interp_sigma


  integer(int32) function count_valid_pts(sig_Er) result(npt)
    real(real64), intent(in) :: sig_Er(:)
    integer(int32) :: ipt
    npt = 0_int32
    do ipt = 1_int32, int(size(sig_Er), int32)
      if (sig_Er(ipt) > 0.0_real64) npt = ipt
    end do
  end function count_valid_pts

end module mod_collisionsGwenael