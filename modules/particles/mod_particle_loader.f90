module mod_particle_loader
  use iso_fortran_env, only: real64, int32
  use omp_lib
  use mod_config,    only: Config
  use mod_particles, only: ParticleSet
  use mod_rng,       only: ran2
  use mod_utils,     only: stop_calculation
  implicit none
  private
  public :: load_particles_modular

contains

  subroutine load_particles_modular(cfg, mpi_rank, mpi_size, n, h, bcnd, kq, vt0, Nm, ni0, iseed, ntype_trk, part, np_thread)

    type(Config),   intent(in)    :: cfg
    integer,        intent(in)    :: mpi_rank, mpi_size
    integer(int32), intent(in)    :: n(3)
    real(real64),   intent(in)    :: h(3)
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(in)    :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(in)    :: vt0(:)
    real(real64),   intent(in)    :: Nm(:)
    real(real64),   intent(in)    :: ni0(:)
    integer(int32), intent(inout) :: iseed(:)
    integer(int32), intent(in)    :: ntype_trk

    type(ParticleSet), intent(inout) :: part(:,:)

    real(real64), intent(out) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype_trk,size(part,2))

    integer(int32) :: nproc
    integer(int32) :: np_cell, n_cell, ix_load
    integer(int32) :: jmax, nmax
    real(real64)   :: x_load, ymax, zmax

    integer :: iproc, ptype, j, k
    integer :: ix, iy, iz

    real(real64) :: x, y, z, vx, vy, vt
    real(real64) :: rnd(2)
    real(real64) :: vz_sav(ntype_trk)

    real(real64) :: px, py, pz, kp, ki(8)

    integer, allocatable      :: np_tot(:,:)
    real(real64), allocatable :: vxp(:,:,:,:)

    integer(int32) :: iseed_omp

    nproc = omp_get_max_threads()

    if (size(part,1) /= ntype_trk) error stop "load_particles_modular: wrong part dim 1"
    if (size(part,2) < nproc)      error stop "load_particles_modular: wrong part dim 2"
    if (size(iseed)  < nproc)      error stop "load_particles_modular: wrong iseed size"
    if (size(np_thread,5) < nproc) error stop "load_particles_modular: wrong np_thread dim 5"

    np_cell = int(cfg%np_cell, int32)

    x_load = cfg%x_load
    if (x_load <= 0.0_real64) then
      x_load = real(n(1),real64) * h(1)
    else
      x_load = min(x_load, real(n(1),real64) * h(1))
    end if

    ymax = real(n(2),real64) * h(2)
    zmax = real(n(3),real64) * h(3)

    ix_load = int(floor(x_load / h(1)), int32) + 1_int32
    if (ix_load < 1_int32)      ix_load = 1_int32
    if (ix_load > n(1)+1_int32) ix_load = n(1)+1_int32

    ! Legacy-compatible count of loadable cells
    n_cell = 0_int32
    do iz = 1, n(3)
      do iy = 1, n(2)
        do ix = 1, ix_load - 1
          if ( bcnd(ix,iy,iz)       == -1 .or. &
               bcnd(ix+1,iy,iz)     == -1 .or. &
               bcnd(ix+1,iy+1,iz)   == -1 .or. &
               bcnd(ix,iy+1,iz)     == -1 .or. &
               bcnd(ix,iy,iz+1)     == -1 .or. &
               bcnd(ix+1,iy,iz+1)   == -1 .or. &
               bcnd(ix+1,iy+1,iz+1) == -1 .or. &
               bcnd(ix,iy+1,iz+1)   == -1 ) then
            n_cell = n_cell + 1_int32
          end if
        end do
      end do
    end do

    jmax = nint(real(np_cell*n_cell, real64) / real(max(1, mpi_size*nproc), real64))
    nmax = int(1.3_real64 * real(max(1,jmax),real64), int32) + 32_int32

    allocate(vxp(6, nmax, ntype_trk, nproc))
    allocate(np_tot(ntype_trk, nproc))

    vxp       = 0.0_real64
    np_thread = 0.0_real64
    np_tot    = 0

    if (mpi_rank == 0) then
      write(*,*) " "
      write(*,'(a)') "Loading particles..."
    end if

    !$omp parallel  &
    !$omp private(iproc,iseed_omp,ptype,j,k,ix,iy,iz,x,y,z,vx,vy,vt,rnd,px,py,pz,kp,ki,vz_sav)

    iproc = omp_get_thread_num() + 1
    if (iproc > nproc) error stop 'iproc > nproc in load_particles_modular'

    iseed_omp = iseed(iproc)

    np_tot(:,iproc) = 0
    np_thread(:,:,:,:,iproc) = 0.0_real64
    vz_sav = 0.0_real64

    do j = 1, jmax

70    continue
      rnd(1) = ran2(iseed_omp)
      x = rnd(1) * x_load
      rnd(1) = ran2(iseed_omp)
      y = rnd(1) * ymax
      rnd(1) = ran2(iseed_omp)
      z = rnd(1) * zmax

      ix = int(x / h(1)) + 1
      iy = int(y / h(2)) + 1
      iz = int(z / h(3)) + 1

      if ( bcnd(ix,iy,iz)       >= 1 .and. &
           bcnd(ix+1,iy,iz)     >= 1 .and. &
           bcnd(ix+1,iy+1,iz)   >= 1 .and. &
           bcnd(ix,iy+1,iz)     >= 1 .and. &
           bcnd(ix,iy,iz+1)     >= 1 .and. &
           bcnd(ix+1,iy,iz+1)   >= 1 .and. &
           bcnd(ix+1,iy+1,iz+1) >= 1 .and. &
           bcnd(ix,iy+1,iz+1)   >= 1 ) then
        goto 70
      end if

      do ptype = 1, ntype_trk

        rnd(1) = ran2(iseed_omp)
        if (rnd(1) > ni0(ptype)) cycle

        np_tot(ptype,iproc) = np_tot(ptype,iproc) + 1
        k = np_tot(ptype,iproc)

        if (k > nmax) then
          print*, 'k > nmax in load_particles_modular'
          call stop_calculation
        end if

        vt = vt0(ptype)

        vxp(1,k,ptype,iproc) = x
        vxp(2,k,ptype,iproc) = y
        vxp(3,k,ptype,iproc) = z

        rnd(1) = ran2(iseed_omp)
        rnd(2) = ran2(iseed_omp)
        call load_gauss(vx,vy,vt,rnd)
        vxp(4,k,ptype,iproc) = vx
        vxp(5,k,ptype,iproc) = vy

        ! Exact legacy behavior
        if (vz_sav(ptype) == 0.0_real64) then
          rnd(1) = ran2(iseed_omp)
          rnd(2) = ran2(iseed_omp)
          call load_gauss(vx,vy,vt,rnd)
          vxp(6,k,ptype,iproc) = vx
          vz_sav(ptype) = vy
        else
          vxp(6,k,ptype,iproc) = vz_sav(ptype)
          vz_sav(ptype) = 0.0_real64
        end if

        px = (ix*h(1) - x)/h(1)
        py = (iy*h(2) - y)/h(2)
        pz = (iz*h(3) - z)/h(3)

        kp = Nm(ptype)/(h(1)*h(2)*h(3))

        ki(1)= kp*px*py*pz
        ki(2)= kp*(1.d0-px)*py*pz
        ki(3)= kp*(1.d0-px)*(1.d0-py)*pz
        ki(4)= kp*px*(1.d0-py)*pz
        ki(5)= kp*px*py*(1.d0-pz)
        ki(6)= kp*(1.d0-px)*py*(1.d0-pz)
        ki(7)= kp*(1.d0-px)*(1.d0-py)*(1.d0-pz)
        ki(8)= kp*px*(1.d0-py)*(1.d0-pz)

        np_thread(ix,iy,iz,ptype,iproc)       = np_thread(ix,iy,iz,ptype,iproc)       + kq(ix,iy,iz)*ki(1)
        np_thread(ix+1,iy,iz,ptype,iproc)     = np_thread(ix+1,iy,iz,ptype,iproc)     + kq(ix+1,iy,iz)*ki(2)
        np_thread(ix+1,iy+1,iz,ptype,iproc)   = np_thread(ix+1,iy+1,iz,ptype,iproc)   + kq(ix+1,iy+1,iz)*ki(3)
        np_thread(ix,iy+1,iz,ptype,iproc)     = np_thread(ix,iy+1,iz,ptype,iproc)     + kq(ix,iy+1,iz)*ki(4)
        np_thread(ix,iy,iz+1,ptype,iproc)     = np_thread(ix,iy,iz+1,ptype,iproc)     + kq(ix,iy,iz+1)*ki(5)
        np_thread(ix+1,iy,iz+1,ptype,iproc)   = np_thread(ix+1,iy,iz+1,ptype,iproc)   + kq(ix+1,iy,iz+1)*ki(6)
        np_thread(ix+1,iy+1,iz+1,ptype,iproc) = np_thread(ix+1,iy+1,iz+1,ptype,iproc) + kq(ix+1,iy+1,iz+1)*ki(7)
        np_thread(ix,iy+1,iz+1,ptype,iproc)   = np_thread(ix,iy+1,iz+1,ptype,iproc)   + kq(ix,iy+1,iz+1)*ki(8)

      end do
    end do

    iseed(iproc) = iseed_omp

    !$omp end parallel

    do iproc = 1, nproc
      do ptype = 1, ntype_trk
        call part(ptype,iproc)%from_vxp(vxp(:,:,ptype,iproc), int(np_tot(ptype,iproc),int32), int(ptype,int32))
      end do
    end do

    if (mpi_rank == 0) then
      do ptype = 1, ntype_trk
        write(*,'(a,i0,a,i0)') "species ", ptype, ": ", sum(np_tot(ptype,1:nproc))
      end do
    end if

    deallocate(vxp, np_tot)

  end subroutine load_particles_modular


  subroutine load_gauss(vx, vy, vt, rnd)
    use mod_constants, only: pi
    real(real64), intent(out) :: vx, vy
    real(real64), intent(in)  :: vt
    real(real64), intent(in)  :: rnd(2)
    real(real64) :: theta, vp

    vp    = vt * sqrt(-log(1.0_real64 - rnd(1)))
    theta = 2.0_real64 * pi * rnd(2)
    vx    = vp * cos(theta)
    vy    = vp * sin(theta)
  end subroutine load_gauss

end module mod_particle_loader