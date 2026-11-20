module mod_planeMoments
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  use mod_constants,   only: qe
  use mod_simParams,   only: SimParams
  implicit none
  private

  public :: compute_particle_plane_moments_species

contains

  subroutine compute_particle_plane_moments_species(part, nproc, n, h, ptype, mass_species, &
                                                    ix_plane, iy_plane, iz_plane,             &
                                                    data_xy, data_xz, data_yz, params)
    ! Plane moments (density, mean velocity, temperature) of one species
    ! on the xy/xz/yz diagnostic planes.
    !
    ! Parallel over iproc: each iproc accumulates into its own slice of
    ! a_xy/a_xz/a_yz (7 sums per plane point: weight, w*v(1:3),
    ! w*v(1:3)**2), the slices are summed, then the per-cell
    ! normalisation below runs as before. This loop used to run on one
    ! thread over every particle of the rank (~370 ms per call at ITER
    ! 2x16). Only the summation order differs from the serial version.
    type(ParticleSet), intent(in)    :: part(:)
    integer(int32),    intent(in)    :: nproc
    integer(int32),    intent(in)    :: n(3)
    integer(int32),    intent(in)    :: ptype
    integer(int32),    intent(in)    :: ix_plane, iy_plane, iz_plane
    real(real64),      intent(in)    :: h(3)
    real(real64),      intent(in)    :: mass_species
    real(real64),      intent(inout) :: data_xy(5,0:n(1)+2,0:n(2)+2)
    real(real64),      intent(inout) :: data_xz(5,0:n(1)+2,0:n(3)+2)
    real(real64),      intent(inout) :: data_yz(5,0:n(2)+2,0:n(3)+2)
    type(SimParams),   intent(in)    :: params
    integer(int32), parameter :: np_avg=1, Tp_avg=2, u1_avg=3, u2_avg=4, u3_avg=5
    integer(int32) :: iproc, i, k
    integer(int32) :: ix, iy, iz
    real(real64)   :: x, y, z, vx, vy, vz, half_dt
    real(real64)   :: px, py, pz
    real(real64)   :: pplane
    real(real64)   :: wx1, wx2, wy1, wy2, wz1, wz2
    real(real64)   :: w11, w21, w22, w12
    real(real64)   :: q(7)
    real(real64)   :: cnt, u2, v2mean, thermal_v2
    real(real64), allocatable :: a_xy(:,:,:,:), a_xz(:,:,:,:), a_yz(:,:,:,:)
    real(real64), allocatable :: t_xy(:,:,:),   t_xz(:,:,:),   t_yz(:,:,:)

    allocate(a_xy(7,0:n(1)+2,0:n(2)+2,nproc), a_xz(7,0:n(1)+2,0:n(3)+2,nproc), &
             a_yz(7,0:n(2)+2,0:n(3)+2,nproc))
    allocate(t_xy(7,0:n(1)+2,0:n(2)+2), t_xz(7,0:n(1)+2,0:n(3)+2), t_yz(7,0:n(2)+2,0:n(3)+2))

    half_dt = params%dt * 0.5_real64

    !$omp parallel do private(iproc,i,x,y,z,vx,vy,vz,ix,iy,iz,px,py,pz,pplane, &
    !$omp&   wx1,wx2,wy1,wy2,wz1,wz2,w11,w21,w22,w12,q) schedule(static)
    do iproc = 1, nproc
      ! Zeroed by the owning thread (first touch).
      a_xy(:,:,:,iproc) = 0.0_real64
      a_xz(:,:,:,iproc) = 0.0_real64
      a_yz(:,:,:,iproc) = 0.0_real64

      if (.not. allocated(part(iproc)%pv)) cycle
      if (part(iproc)%n <= 0_int32) cycle

      do i = 1, part(iproc)%n
        x = part(iproc)%pv(1,i) - part(iproc)%pv(4,i) * half_dt
        y = part(iproc)%pv(2,i) - part(iproc)%pv(5,i) * half_dt
        z = part(iproc)%pv(3,i) - part(iproc)%pv(6,i) * half_dt
        vx = part(iproc)%pv(4,i)
        vy = part(iproc)%pv(5,i)
        vz = part(iproc)%pv(6,i)

        ix = int(x / h(1), int32) + 1_int32
        iy = int(y / h(2), int32) + 1_int32
        iz = int(z / h(3), int32) + 1_int32
        if (ix < 1_int32 .or. ix > n(1)) cycle
        if (iy < 1_int32 .or. iy > n(2)) cycle
        if (iz < 1_int32 .or. iz > n(3)) cycle

        ! Skip the weight math for the vast majority of particles, which
        ! are on none of the three planes.
        if (iz /= iz_plane .and. iz /= iz_plane-1_int32 .and. &
            iy /= iy_plane .and. iy /= iy_plane-1_int32 .and. &
            ix /= ix_plane .and. ix /= ix_plane-1_int32) cycle

        px = (real(ix, real64)*h(1) - x) / h(1)
        py = (real(iy, real64)*h(2) - y) / h(2)
        pz = (real(iz, real64)*h(3) - z) / h(3)
        px = max(0.0_real64, min(1.0_real64, px))
        py = max(0.0_real64, min(1.0_real64, py))
        pz = max(0.0_real64, min(1.0_real64, pz))
        wx1 = px
        wx2 = 1.0_real64 - px
        wy1 = py
        wy2 = 1.0_real64 - py
        wz1 = pz
        wz2 = 1.0_real64 - pz

        q = [1.0_real64, vx, vy, vz, vx*vx, vy*vy, vz*vz]

        if (iz == iz_plane .or. iz == iz_plane-1_int32) then
          if (iz == iz_plane) then
            pplane = pz
          else
            pplane = 1.0_real64 - pz
          end if
          w11 = pplane * wx1 * wy1
          w21 = pplane * wx2 * wy1
          w22 = pplane * wx2 * wy2
          w12 = pplane * wx1 * wy2
          do k = 1, 7
            a_xy(k,ix  ,iy  ,iproc) = a_xy(k,ix  ,iy  ,iproc) + w11*q(k)
            a_xy(k,ix+1,iy  ,iproc) = a_xy(k,ix+1,iy  ,iproc) + w21*q(k)
            a_xy(k,ix+1,iy+1,iproc) = a_xy(k,ix+1,iy+1,iproc) + w22*q(k)
            a_xy(k,ix  ,iy+1,iproc) = a_xy(k,ix  ,iy+1,iproc) + w12*q(k)
          end do
        end if

        if (iy == iy_plane .or. iy == iy_plane-1_int32) then
          if (iy == iy_plane) then
            pplane = py
          else
            pplane = 1.0_real64 - py
          end if
          w11 = pplane * wx1 * wz1
          w21 = pplane * wx2 * wz1
          w22 = pplane * wx2 * wz2
          w12 = pplane * wx1 * wz2
          do k = 1, 7
            a_xz(k,ix  ,iz  ,iproc) = a_xz(k,ix  ,iz  ,iproc) + w11*q(k)
            a_xz(k,ix+1,iz  ,iproc) = a_xz(k,ix+1,iz  ,iproc) + w21*q(k)
            a_xz(k,ix+1,iz+1,iproc) = a_xz(k,ix+1,iz+1,iproc) + w22*q(k)
            a_xz(k,ix  ,iz+1,iproc) = a_xz(k,ix  ,iz+1,iproc) + w12*q(k)
          end do
        end if

        if (ix == ix_plane .or. ix == ix_plane-1_int32) then
          if (ix == ix_plane) then
            pplane = px
          else
            pplane = 1.0_real64 - px
          end if
          w11 = pplane * wy1 * wz1
          w21 = pplane * wy2 * wz1
          w22 = pplane * wy2 * wz2
          w12 = pplane * wy1 * wz2
          do k = 1, 7
            a_yz(k,iy  ,iz  ,iproc) = a_yz(k,iy  ,iz  ,iproc) + w11*q(k)
            a_yz(k,iy+1,iz  ,iproc) = a_yz(k,iy+1,iz  ,iproc) + w21*q(k)
            a_yz(k,iy+1,iz+1,iproc) = a_yz(k,iy+1,iz+1,iproc) + w22*q(k)
            a_yz(k,iy  ,iz+1,iproc) = a_yz(k,iy  ,iz+1,iproc) + w12*q(k)
          end do
        end if
      end do
    end do
    !$omp end parallel do

    ! Sum the per-iproc slices.
    !$omp parallel do private(iy,iproc) schedule(static)
    do iy = 0, n(2)+2
      t_xy(:,:,iy) = a_xy(:,:,iy,1)
      do iproc = 2, nproc
        t_xy(:,:,iy) = t_xy(:,:,iy) + a_xy(:,:,iy,iproc)
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(iz,iproc) schedule(static)
    do iz = 0, n(3)+2
      t_xz(:,:,iz) = a_xz(:,:,iz,1)
      t_yz(:,:,iz) = a_yz(:,:,iz,1)
      do iproc = 2, nproc
        t_xz(:,:,iz) = t_xz(:,:,iz) + a_xz(:,:,iz,iproc)
        t_yz(:,:,iz) = t_yz(:,:,iz) + a_yz(:,:,iz,iproc)
      end do
    end do
    !$omp end parallel do

    ! Normalise. Velocity sums are added into data_*(u1:u3) (as the
    ! serial version accumulated straight into them) and divided by the
    ! weight only where it is positive.
    do iy = 0, n(2)+2
      do ix = 0, n(1)+2
        data_xy(u1_avg:u3_avg,ix,iy) = data_xy(u1_avg:u3_avg,ix,iy) + t_xy(2:4,ix,iy)
        cnt = t_xy(1,ix,iy)
        if (cnt > 0.0_real64) then
          data_xy(np_avg,ix,iy) = cnt
          data_xy(u1_avg,ix,iy) = data_xy(u1_avg,ix,iy) / cnt
          data_xy(u2_avg,ix,iy) = data_xy(u2_avg,ix,iy) / cnt
          data_xy(u3_avg,ix,iy) = data_xy(u3_avg,ix,iy) / cnt
          u2 = data_xy(u1_avg,ix,iy)**2 + data_xy(u2_avg,ix,iy)**2 + data_xy(u3_avg,ix,iy)**2
          v2mean = (t_xy(5,ix,iy) + t_xy(6,ix,iy) + t_xy(7,ix,iy)) / cnt
          thermal_v2 = max(0.0_real64, v2mean - u2)
          data_xy(Tp_avg,ix,iy) = data_xy(Tp_avg,ix,iy) + &
                        mass_species * thermal_v2 / (3.0_real64 * qe)
        end if
      end do
    end do
    do iz = 0, n(3)+2
      do ix = 0, n(1)+2
        data_xz(u1_avg:u3_avg,ix,iz) = data_xz(u1_avg:u3_avg,ix,iz) + t_xz(2:4,ix,iz)
        cnt = t_xz(1,ix,iz)
        if (cnt > 0.0_real64) then
          data_xz(np_avg,ix,iz) = cnt
          data_xz(u1_avg,ix,iz) = data_xz(u1_avg,ix,iz) / cnt
          data_xz(u2_avg,ix,iz) = data_xz(u2_avg,ix,iz) / cnt
          data_xz(u3_avg,ix,iz) = data_xz(u3_avg,ix,iz) / cnt
          u2 = data_xz(u1_avg,ix,iz)**2 + data_xz(u2_avg,ix,iz)**2 + data_xz(u3_avg,ix,iz)**2
          v2mean = (t_xz(5,ix,iz) + t_xz(6,ix,iz) + t_xz(7,ix,iz)) / cnt
          thermal_v2 = max(0.0_real64, v2mean - u2)
          data_xz(Tp_avg,ix,iz) = data_xz(Tp_avg,ix,iz) + &
                      mass_species * thermal_v2 / (3.0_real64 * qe)
        end if
      end do
    end do
    do iz = 0, n(3)+2
      do iy = 0, n(2)+2
        data_yz(u1_avg:u3_avg,iy,iz) = data_yz(u1_avg:u3_avg,iy,iz) + t_yz(2:4,iy,iz)
        cnt = t_yz(1,iy,iz)
        if (cnt > 0.0_real64) then
          data_yz(np_avg,iy,iz) = cnt
          data_yz(u1_avg,iy,iz) = data_yz(u1_avg,iy,iz) / cnt
          data_yz(u2_avg,iy,iz) = data_yz(u2_avg,iy,iz) / cnt
          data_yz(u3_avg,iy,iz) = data_yz(u3_avg,iy,iz) / cnt
          u2 = data_yz(u1_avg,iy,iz)**2 + data_yz(u2_avg,iy,iz)**2 + data_yz(u3_avg,iy,iz)**2
          v2mean = (t_yz(5,iy,iz) + t_yz(6,iy,iz) + t_yz(7,iy,iz)) / cnt
          thermal_v2 = max(0.0_real64, v2mean - u2)
          data_yz(Tp_avg,iy,iz) = data_yz(Tp_avg,iy,iz) + &
                  mass_species * thermal_v2 / (3.0_real64 * qe)
        end if
      end do
    end do

    deallocate(a_xy, a_xz, a_yz, t_xy, t_xz, t_yz)
  end subroutine compute_particle_plane_moments_species
end module mod_planeMoments