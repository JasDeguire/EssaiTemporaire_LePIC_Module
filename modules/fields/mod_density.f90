module mod_density
  !!
  !! Fast density reduction / charge-density builder
  !!
  !! Main performance change vs the previous modular version:
  !!   - No full-size temporary np_work allocation/copy every timestep.
  !!   - Periodic stitching is applied directly to np_thread, exactly where the
  !!     old reduce_species_density already did it.
  !!   - build_rho_from_np now reads np_red directly instead of copying it.
  !!
  !! This keeps the same numerical convention as your current modular code:
  !!   reduce_species_density:
  !!      np_thread -> periodic density stitching -> np_red
  !!   build_rho_from_np:
  !!      rho = - sum_s charge(s) * np_red(s)
  !!
  use iso_fortran_env, only: real64, int32
  use mpi
  use mod_constants, only: qe, eps0

  implicit none
  private

  public :: reduce_species_density
  public :: reduce_density_and_rho
  public :: sync_species_density
  public :: clear_np_planes
  public :: build_rho_from_np
  public :: build_rho_from_np_thread
  public :: density_max_per_species
  public :: average_species_density

contains

  subroutine reduce_species_density(n, bcnd, np_thread, ntype, nproc, mpi_comm, np_red)
    ! np_thread -> periodic stitching -> np_red, MPI-summed across ranks.
    ! Used where the full per-species density must be globally valid
    ! right away (initial load, restart). The per-step path uses
    ! reduce_density_and_rho + a deferred sync_species_density instead.
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc, mpi_comm
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
    real(real64),   intent(out)   :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)

    ! Important: this is intentionally in-place.
    ! In your timestep, np_thread is rebuilt by deposition before this is called,
    ! so modifying its periodic/ghost planes here is equivalent to the old
    ! np_work = np_thread copy, but avoids the large allocation/copy.
    call apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)
    call sum_thread_density(n, np_thread, ntype, nproc, np_red)
    call sync_species_density(n, ntype, mpi_comm, np_red)
  end subroutine reduce_species_density


  subroutine reduce_density_and_rho(n, bcnd, np_thread, ntype, nproc, mpi_comm, charge, &
                                    np_red, rho, zlo, zhi, zend)
    ! Per-step path, mirroring legacy calc_rho: rho is built from this
    ! rank's thread sum and only rho (one grid) is MPI-summed. np_red is
    ! left RANK-LOCAL - the caller must run sync_species_density on it
    ! before anything reads it (legacy likewise only reduces the
    ! per-species density, dens_red, on the steps that need it). Saves
    ! allreducing ntype grids every step (~30 ms/step at ITER 2x16).
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc, mpi_comm
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
    real(real64),   intent(in)    :: charge(ntype)
    real(real64),   intent(out)   :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(out)   :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    ! Optional per-(ptype,iproc) z-extents of np_thread (see
    ! np_plane_live); planes outside them are known to be zero and are
    ! skipped in the sum.
    integer(int32), intent(in), optional :: zlo(ntype,nproc), zhi(ntype,nproc), zend(ntype,nproc)

    integer :: ierr, mpi_size

    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    call apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)
    call sum_thread_density(n, np_thread, ntype, nproc, np_red, charge, rho, zlo, zhi, zend)

    if (mpi_size > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, rho, (n(1)+3)*(n(2)+3)*(n(3)+3), &
          MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
    end if
  end subroutine reduce_density_and_rho


  subroutine sync_species_density(n, ntype, mpi_comm, np_red)
    ! MPI-sum a rank-local np_red in place (no-op on a single rank).
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, mpi_comm
    real(real64),   intent(inout) :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)

    integer :: ierr, mpi_size

    call MPI_Comm_size(mpi_comm, mpi_size, ierr)
    if (mpi_size > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, np_red, &
          (n(1)+3)*(n(2)+3)*(n(3)+3)*ntype, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
    end if
  end subroutine sync_species_density


  subroutine sum_thread_density(n, np_thread, ntype, nproc, np_red, charge, rho, zlo, zhi, zend)
    ! np_red = sum over iproc of np_thread on interior points 1..n+1,
    ! zero on the ghost shell. If charge/rho are present, also
    ! rho = -sum_s charge(s)*np_red(s) (same operation order as
    ! build_rho_from_np). Works a row (fixed iy,iz) at a time with iproc
    ! outside the contiguous ix loop, so every access is unit-stride and
    ! vectorizable and the row being accumulated stays in L1 - the
    ! previous per-point loop over iproc did nproc strided scalar loads per
    ! point, and the full-array zero before it ran on one thread.
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype, nproc
    real(real64),   intent(in)  :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
    real(real64),   intent(out) :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(in),  optional :: charge(ntype)
    real(real64),   intent(out), optional :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in),  optional :: zlo(ntype,nproc), zhi(ntype,nproc), zend(ntype,nproc)

    integer :: ix, iy, iz, ptype, iproc, nx1
    logical :: do_rho, ranged

    do_rho = present(rho) .and. present(charge)
    ranged = present(zlo) .and. present(zhi) .and. present(zend)
    nx1 = n(1) + 1

    !$omp parallel do collapse(2) private(iy,iz,ix,ptype,iproc) schedule(static) default(shared)
    do iz = 0, n(3)+2
      do iy = 0, n(2)+2
        if (iz < 1 .or. iz > n(3)+1 .or. iy < 1 .or. iy > n(2)+1) then
          do ptype = 1, ntype
            np_red(:,iy,iz,ptype) = 0.0_real64
          end do
          if (do_rho) rho(:,iy,iz) = 0.0_real64
          cycle
        end if

        do ptype = 1, ntype
          np_red(0,iy,iz,ptype)      = 0.0_real64
          np_red(nx1+1,iy,iz,ptype)  = 0.0_real64
          ! Start from zero and add only live planes. Bit-identical to the
          ! full sum: every skipped term is exactly 0, and x + 0 == x.
          do ix = 1, nx1
            np_red(ix,iy,iz,ptype) = 0.0_real64
          end do
          do iproc = 1, nproc
            if (ranged) then
              if (.not. np_plane_live(iz, zlo(ptype,iproc), zhi(ptype,iproc), &
                                      zend(ptype,iproc), n(3))) cycle
            end if
            do ix = 1, nx1
              np_red(ix,iy,iz,ptype) = np_red(ix,iy,iz,ptype) + np_thread(ix,iy,iz,ptype,iproc)
            end do
          end do
        end do

        if (do_rho) then
          rho(0,iy,iz)     = 0.0_real64
          rho(nx1+1,iy,iz) = 0.0_real64
          do ix = 1, nx1
            rho(ix,iy,iz) = 0.0_real64
          end do
          do ptype = 1, ntype
            do ix = 1, nx1
              rho(ix,iy,iz) = rho(ix,iy,iz) - charge(ptype) * np_red(ix,iy,iz,ptype)
            end do
          end do
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine sum_thread_density


  pure logical function np_plane_live(iz, lo, hi, zend, n3) result(live)
    ! True if plane iz of an np_thread slice may be nonzero, given its
    ! tracked extent: planes lo..hi, plus the end blocks 0:2 and n3:n3+2
    ! when zend /= 0 (what apply_periodic_density_bc's z-stitching can
    ! write when the deposit reached near either end).
    integer, intent(in)        :: iz
    integer(int32), intent(in) :: lo, hi, zend, n3

    live = (iz >= lo .and. iz <= hi)
    if (.not. live .and. zend /= 0_int32) live = (iz <= 2 .or. iz >= n3)
  end function np_plane_live


  subroutine clear_np_planes(np_local, n, lo, hi, zend)
    ! Zero exactly the planes np_plane_live reports for (lo,hi,zend),
    ! leaving the whole slice zero (everything else already was).
    integer(int32), intent(in)    :: n(3)
    real(real64),   intent(inout) :: np_local(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)    :: lo, hi, zend

    integer :: iz

    do iz = 0, n(3)+2
      if (np_plane_live(iz, lo, hi, zend, n(3))) np_local(:,:,iz) = 0.0_real64
    end do
  end subroutine clear_np_planes


  subroutine build_rho_from_np(n, np_red, charge, ntype, rho, bcnd, flag_pbc)
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype
    real(real64),   intent(in)  :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(in)  :: charge(ntype)
    real(real64),   intent(out) :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)  :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)  :: flag_pbc

    integer :: ix, iy, iz, ptype
    real(real64) :: r

    ! Keep current modular behavior: no second periodic correction here.
    ! The previous file had that correction commented out after copying np_red
    ! into np_work. Therefore we read np_red directly.
    rho = 0.0_real64

    !$omp parallel do collapse(3) private(ptype,r) schedule(static) default(shared)
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          r = 0.0_real64
          do ptype = 1, ntype
            r = r - charge(ptype) * np_red(ix,iy,iz,ptype)
          end do
          rho(ix,iy,iz) = r
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine build_rho_from_np


  subroutine build_rho_from_np_thread(n, np_thread, charge, ntype, nproc, rho, bcnd, flag_pbc)
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
    real(real64),   intent(in)    :: charge(ntype)
    real(real64),   intent(out)   :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)    :: flag_pbc

    integer :: ix, iy, iz, ptype, iproc
    real(real64) :: r

    ! This routine is not used in your current mod_simulation timestep path,
    ! but it is kept as a fast drop-in equivalent.
    if (flag_pbc == 1_int32) then
      call apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)
    end if

    rho = 0.0_real64

    !$omp parallel do collapse(3) private(iproc,ptype,r) schedule(static) default(shared)
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          r = 0.0_real64
          do iproc = 1, nproc
            do ptype = 1, ntype
              r = r - charge(ptype) * np_thread(ix,iy,iz,ptype,iproc)
            end do
          end do
          rho(ix,iy,iz) = r
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine build_rho_from_np_thread


  subroutine density_max_per_species(n, bcnd, np_red, ntype, np_mx)
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype
    integer,        intent(in)  :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(in)  :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(out) :: np_mx(ntype)

    integer :: ix, iy, iz, ptype
    real(real64) :: local_max

    np_mx = 0.0_real64

    do ptype = 1, ntype
      local_max = 0.0_real64

      !$omp parallel do collapse(3) reduction(max:local_max) schedule(static) default(shared)
      do iz = 1, n(3)+1
        do iy = 1, n(2)+1
          do ix = 1, n(1)+1
            local_max = max(local_max, np_red(ix,iy,iz,ptype))
          end do
        end do
      end do
      !$omp end parallel do

      np_mx(ptype) = local_max
    end do

  end subroutine density_max_per_species


  ! Spatial average of np_red(:,:,:,ptype) over the physical node grid.
  ! np_red is populated on nodes 1..n(dim)+1 in each dimension; when a
  ! dimension is periodic, apply_periodic_density_bc makes node n(dim)+1
  ! a duplicate of node 1, so that node is excluded here to avoid double
  ! counting the seam plane.
  real(real64) function average_species_density(n, np_red, ntype, ptype, flag_pbc, flag_pbcz)
    integer(int32), intent(in) :: n(3)
    integer,        intent(in) :: ntype, ptype
    real(real64),   intent(in) :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    integer(int32), intent(in) :: flag_pbc, flag_pbcz

    integer :: ny_uniq, nz_uniq

    ny_uniq = merge(n(2), n(2)+1, flag_pbc  == 1_int32)
    nz_uniq = merge(n(3), n(3)+1, flag_pbcz == 1_int32)

    average_species_density = sum(np_red(1:n(1)+1, 1:ny_uniq, 1:nz_uniq, ptype)) / &
        real((n(1)+1) * ny_uniq * nz_uniq, real64)
  end function average_species_density


  subroutine apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)

    integer :: ix, iy, iz, iproc, ptype

    ! ------------------------------------------------------------
    ! Legacy periodic density stitching in y.
    ! Written with explicit ptype loop to avoid array-section temporaries.
    ! ------------------------------------------------------------
    !$omp parallel do collapse(3) private(ptype) schedule(static) default(shared)
    do iproc = 1, nproc
      do iz = 1, n(3)+1
        do ix = 1, n(1)+1
          if (bcnd(ix,1,iz) == 0) then
            do ptype = 1, ntype
              np_thread(ix,1,iz,ptype,iproc) = 0.5_real64 * ( &
                   np_thread(ix,1,iz,ptype,iproc) + np_thread(ix,n(2)+1,iz,ptype,iproc) )
              np_thread(ix,0,iz,ptype,iproc) = np_thread(ix,n(2),iz,ptype,iproc)
            end do
          end if

          if (bcnd(ix,n(2)+1,iz) == 0) then
            do ptype = 1, ntype
              np_thread(ix,n(2)+2,iz,ptype,iproc) = np_thread(ix,2,iz,ptype,iproc)
              np_thread(ix,n(2)+1,iz,ptype,iproc) = np_thread(ix,1,iz,ptype,iproc)
            end do
          end if
        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------
    ! Legacy periodic density stitching in z.
    ! Written with explicit ptype loop to avoid array-section temporaries.
    ! ------------------------------------------------------------
    !$omp parallel do collapse(3) private(ptype) schedule(static) default(shared)
    do iproc = 1, nproc
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          if (bcnd(ix,iy,1) == 0) then
            do ptype = 1, ntype
              np_thread(ix,iy,1,ptype,iproc) = 0.5_real64 * ( &
                   np_thread(ix,iy,1,ptype,iproc) + np_thread(ix,iy,n(3)+1,ptype,iproc) )
              np_thread(ix,iy,0,ptype,iproc) = np_thread(ix,iy,n(3),ptype,iproc)
            end do
          end if

          if (bcnd(ix,iy,n(3)+1) == 0) then
            do ptype = 1, ntype
              np_thread(ix,iy,n(3)+2,ptype,iproc) = np_thread(ix,iy,2,ptype,iproc)
              np_thread(ix,iy,n(3)+1,ptype,iproc) = np_thread(ix,iy,1,ptype,iproc)
            end do
          end if
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_periodic_density_bc

end module mod_density
