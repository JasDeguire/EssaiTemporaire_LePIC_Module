module mod_chargeDeposition
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  implicit none
  private

  public :: clear_np_thread
  public :: deposit_particle_set_to_np_thread
  public :: HeatRegionTally

  ! Keep this .false. for production speed.
  ! Switch to .true. only while debugging bad particle positions.
  logical, parameter :: DEPOSITION_SAFETY_CHECKS = .false.

  ! Electron heating-region tally (Nh, sum_dEk), accumulated inside the
  ! deposit loop instead of a separate pass over all electrons - legacy
  ! does the same inside its charge_deposition. The separate pass cost
  ! ~70 ms/step at ITER scale (2x16), mostly memory traffic re-reading pv
  ! plus false sharing on the per-iproc Nh/sum_dEk accumulators.
  ! ixl..ek_coef are inputs describing the region; Nh/sum_dEk accumulate.
  type :: HeatRegionTally
    integer(int32) :: ixl = 0_int32, ixr = -1_int32
    integer(int32) :: flag_circxh = 0_int32, flag_ahp = 0_int32
    real(real64)   :: R2 = 0.0_real64, yc = 0.0_real64, zc = 0.0_real64
    real(real64)   :: ek_coef = 0.0_real64      ! 0.5*Nm*mass
    ! Skin-depth weighting w = exp((r-R)/delta) (delta <= 0: w = 1);
    ! sum_dEk and sum_w are then sums of w*Ek and w.
    real(real64)   :: R = 0.0_real64, delta = 0.0_real64
    integer(int32) :: Nh = 0_int32
    real(real64)   :: sum_dEk = 0.0_real64
    real(real64)   :: sum_w = 0.0_real64
  end type HeatRegionTally

contains

  subroutine clear_np_thread(n, ntype, nproc, np_thread)
    integer(int32), intent(in)    :: n(3)
    integer(int32), intent(in)    :: ntype, nproc
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)

    np_thread = 0.0_real64
  end subroutine clear_np_thread


  subroutine deposit_particle_set_to_np_thread(part, n, h, kq, Nm_species, np_local, heat, &
                                               iz_lo, iz_hi)
    ! Fast production version of charge deposition.
    ! Same deposition convention as the previous modular code, but the
    ! expensive debug/error guards are compile-time disabled by default.
    ! If heat is present, also tallies live particles in the heating
    ! region into heat%Nh/heat%sum_dEk (see HeatRegionTally).
    ! If iz_lo/iz_hi are present, returns the range of z-planes this call
    ! wrote into np_local (iz_lo > iz_hi if it wrote nothing).
    type(ParticleSet), intent(in)    :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: Nm_species
    real(real64),       intent(inout) :: np_local(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    type(HeatRegionTally), intent(inout), optional :: heat
    integer(int32),     intent(out),   optional :: iz_lo, iz_hi

    integer(int32) :: zlo_loc, zhi_loc

    logical        :: do_heat, in_heat
    integer(int32) :: ixh, hixl, hixr, hcirc, hahp, nh_loc
    real(real64)   :: hR2, hyc, hzc, hek, r2, sdek_loc
    real(real64)   :: hR, hdelta, wh, sw_loc

    integer(int32) :: i
    integer(int32) :: ix, iy, iz
    real(real64)   :: x, y, z
    real(real64)   :: px, py, pz
    real(real64)   :: wx2, wy2, wz2
    real(real64)   :: k1
    real(real64)   :: w1, w2, w3, w4, w5, w6, w7, w8
    real(real64)   :: xmax_loc, ymax_loc, zmax_loc
    logical        :: has_dead

    if (present(iz_lo)) iz_lo = huge(1_int32)
    if (present(iz_hi)) iz_hi = -huge(1_int32)
    zlo_loc = huge(1_int32)
    zhi_loc = -huge(1_int32)

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    xmax_loc = real(n(1), real64) * h(1)
    ymax_loc = real(n(2), real64) * h(2)
    zmax_loc = real(n(3), real64) * h(3)

    k1 = Nm_species / (h(1) * h(2) * h(3))
    has_dead = allocated(part%flag_dead)

    ! Heat-tally parameters copied to locals so the loop has no aliasing
    ! concerns and the tally accumulates in registers.
    do_heat  = present(heat)
    nh_loc   = 0_int32
    sdek_loc = 0.0_real64
    sw_loc   = 0.0_real64
    hixl = 0_int32; hixr = -1_int32; hcirc = 0_int32; hahp = 0_int32
    hR2 = 0.0_real64; hyc = 0.0_real64; hzc = 0.0_real64; hek = 0.0_real64
    hR = 0.0_real64; hdelta = 0.0_real64
    if (do_heat) then
      hixl = heat%ixl;  hixr = heat%ixr
      hcirc = heat%flag_circxh;  hahp = heat%flag_ahp
      hR2 = heat%R2;  hyc = heat%yc;  hzc = heat%zc;  hek = heat%ek_coef
      hR = heat%R;  hdelta = heat%delta
    end if

    do i = 1, part%n
      if (has_dead) then
        if (part%flag_dead(i) /= 0) cycle
      end if

      x = part%pv(1,i)
      y = part%pv(2,i)
      z = part%pv(3,i)

      ! Heating-region tally: before the roundoff guard below, since the
      ! separate pass this replaces counted every live particle.
      if (do_heat) then
        ixh = int(x / h(1), int32) + 1_int32
        if (ixh >= hixl .and. ixh <= hixr) then
          in_heat = .true.
          if (hcirc == 1_int32) then
            r2 = (y-hyc)**2 + (z-hzc)**2
            if (hahp == 0_int32) then
              in_heat = (r2 <= hR2)
            else
              in_heat = (r2 >= hR2)
            end if
          end if
          if (in_heat) then
            wh = 1.0_real64
            if (hdelta > 0.0_real64) wh = exp((sqrt((y-hyc)**2 + (z-hzc)**2) - hR)/hdelta)
            sdek_loc = sdek_loc + wh * hek * (part%pv(4,i)**2 + part%pv(5,i)**2 + part%pv(6,i)**2)
            sw_loc   = sw_loc + wh
            nh_loc   = nh_loc + 1_int32
          end if
        end if
      end if

      ! Cheap guard kept for particles slightly outside due to roundoff.
      ! This matches the previous behavior where tiny excursions were skipped.
      if (x < -1.0e-12_real64 .or. x > xmax_loc + 1.0e-12_real64) cycle
      if (y < -1.0e-12_real64 .or. y > ymax_loc + 1.0e-12_real64) cycle
      if (z < -1.0e-12_real64 .or. z > zmax_loc + 1.0e-12_real64) cycle

      if (DEPOSITION_SAFETY_CHECKS) then
        if (.not. (x == x .and. y == y .and. z == z)) then
          write(*,*) 'NaN particle position in deposition'
          write(*,*) 'i = ', i
          write(*,*) 'x,y,z = ', x, y, z
          error stop 'deposit_particle_set_to_np_thread: NaN position'
        end if

        if (x < 0.0_real64 .or. x > xmax_loc .or. &
            y < 0.0_real64 .or. y > ymax_loc .or. &
            z < 0.0_real64 .or. z > zmax_loc) then
          write(*,*) 'Out-of-range particle position in deposition'
          write(*,*) 'i = ', i
          write(*,*) 'x,y,z = ', x, y, z
          write(*,*) 'xmax,ymax,zmax = ', xmax_loc, ymax_loc, zmax_loc
          error stop 'deposit_particle_set_to_np_thread: particle out of bounds'
        end if
      end if

      ix = int(x / h(1), int32) + 1_int32
      iy = int(y / h(2), int32) + 1_int32
      iz = int(z / h(3), int32) + 1_int32
      zlo_loc = min(zlo_loc, iz)
      zhi_loc = max(zhi_loc, iz + 1_int32)

      if (DEPOSITION_SAFETY_CHECKS) then
        if (ix < 1_int32 .or. ix > n(1)+1_int32 .or. &
            iy < 1_int32 .or. iy > n(2)+1_int32 .or. &
            iz < 1_int32 .or. iz > n(3)+1_int32) then
          write(*,*) 'Bad deposition cell index'
          write(*,*) 'i = ', i
          write(*,*) 'ix,iy,iz = ', ix, iy, iz
          write(*,*) 'x,y,z = ', x, y, z
          write(*,*) 'n = ', n
          error stop 'deposit_particle_set_to_np_thread: invalid cell index'
        end if
      end if

      px = (real(ix, real64) * h(1) - x) / h(1)
      py = (real(iy, real64) * h(2) - y) / h(2)
      pz = (real(iz, real64) * h(3) - z) / h(3)

      wx2 = 1.0_real64 - px
      wy2 = 1.0_real64 - py
      wz2 = 1.0_real64 - pz

      w1 = k1 * px  * py  * pz
      w2 = k1 * wx2 * py  * pz
      w3 = k1 * wx2 * wy2 * pz
      w4 = k1 * px  * wy2 * pz
      w5 = k1 * px  * py  * wz2
      w6 = k1 * wx2 * py  * wz2
      w7 = k1 * wx2 * wy2 * wz2
      w8 = k1 * px  * wy2 * wz2

      np_local(ix  ,iy  ,iz  ) = np_local(ix  ,iy  ,iz  ) + kq(ix  ,iy  ,iz  ) * w1
      np_local(ix+1,iy  ,iz  ) = np_local(ix+1,iy  ,iz  ) + kq(ix+1,iy  ,iz  ) * w2
      np_local(ix+1,iy+1,iz  ) = np_local(ix+1,iy+1,iz  ) + kq(ix+1,iy+1,iz  ) * w3
      np_local(ix  ,iy+1,iz  ) = np_local(ix  ,iy+1,iz  ) + kq(ix  ,iy+1,iz  ) * w4
      np_local(ix  ,iy  ,iz+1) = np_local(ix  ,iy  ,iz+1) + kq(ix  ,iy  ,iz+1) * w5
      np_local(ix+1,iy  ,iz+1) = np_local(ix+1,iy  ,iz+1) + kq(ix+1,iy  ,iz+1) * w6
      np_local(ix+1,iy+1,iz+1) = np_local(ix+1,iy+1,iz+1) + kq(ix+1,iy+1,iz+1) * w7
      np_local(ix  ,iy+1,iz+1) = np_local(ix  ,iy+1,iz+1) + kq(ix  ,iy+1,iz+1) * w8
    end do

    if (present(iz_lo)) iz_lo = zlo_loc
    if (present(iz_hi)) iz_hi = zhi_loc

    if (do_heat) then
      heat%Nh      = heat%Nh + nh_loc
      heat%sum_dEk = heat%sum_dEk + sdek_loc
      heat%sum_w   = heat%sum_w + sw_loc
    end if
  end subroutine deposit_particle_set_to_np_thread

end module mod_chargeDeposition