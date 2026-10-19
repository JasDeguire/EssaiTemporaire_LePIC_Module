module mod_heating
  use iso_fortran_env, only: int32, real64, int8
  use mod_particles,   only: ParticleSet
  use mod_rng,         only: ran2
  implicit none
  private

  public :: apply_electron_heating

contains

  subroutine apply_electron_heating(part, h, vt, iseed, nudt, xl_pow, xr_pow, &
                                    flag_circxh, flag_ahp, R_ahp, delta_heat,   &
                                    ymax, zmax, Nm_e, mass_e, p_loss_heating)
    class(ParticleSet), intent(inout) :: part
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: vt
    integer(int32),     intent(inout) :: iseed
    real(real64),       intent(in)    :: nudt
    real(real64),       intent(in)    :: xl_pow, xr_pow
    integer(int32),     intent(in)    :: flag_circxh, flag_ahp
    
    real(real64),       intent(in)    :: R_ahp, ymax, zmax
    ! > 0: re-sampling probability scaled by exp((r-R_ahp)/delta_heat)
    real(real64),       intent(in)    :: delta_heat
    real(real64),       intent(in)    :: Nm_e, mass_e
    real(real64),       intent(inout) :: p_loss_heating

    integer(int32) :: i, ix
    real(real64)   :: xp_new, yp_new, zp_new
    real(real64)   :: rnd(2), vz_sav, v2old, wh
    integer(int32)            :: ixl_pow, ixr_pow
    

    ixl_pow = int(xl_pow / h(1), int32) + 1_int32
    ixr_pow = int(xr_pow / h(1), int32) + 1_int32

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    vz_sav = 0.0_real64

    do i = 1, part%n

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) /= 0_int8) cycle
      end if

      xp_new = part%pv(1,i)
      ix = int(xp_new / h(1), int32) + 1_int32

      if (ix < ixl_pow .or. ix > ixr_pow) cycle

      if (flag_circxh == 1_int32) then
        yp_new = part%pv(2,i)
        zp_new = part%pv(3,i)

        if (flag_ahp == 0_int32) then
          if ( ((yp_new - ymax/2.0_real64)**2 + (zp_new - zmax/2.0_real64)**2) > R_ahp**2 ) cycle
        else
          if ( ((yp_new - ymax/2.0_real64)**2 + (zp_new - zmax/2.0_real64)**2) < R_ahp**2 ) cycle
        end if
      end if

      wh = 1.0_real64
      if (delta_heat > 0.0_real64) then
        wh = exp((sqrt((part%pv(2,i) - ymax/2.0_real64)**2 + &
                       (part%pv(3,i) - zmax/2.0_real64)**2) - R_ahp)/delta_heat)
      end if

      rnd(1) = ran2(iseed)
      if (rnd(1) > nudt*wh) cycle

      v2old = part%pv(4,i)*part%pv(4,i) + part%pv(5,i)*part%pv(5,i) + part%pv(6,i)*part%pv(6,i)

      rnd(1) = ran2(iseed)
      rnd(2) = ran2(iseed)
      call load_gauss_local(part%pv(4,i), part%pv(5,i), vt, rnd)

      if (vz_sav == 0.0_real64) then
        rnd(1) = ran2(iseed)
        rnd(2) = ran2(iseed)
        call load_gauss_local(part%pv(6,i), vz_sav, vt, rnd)
      else
        part%pv(6,i) = vz_sav
        vz_sav = 0.0_real64
      end if

      p_loss_heating = p_loss_heating + 0.5_real64 * Nm_e * mass_e * &
            (part%pv(4,i)*part%pv(4,i) + part%pv(5,i)*part%pv(5,i) + part%pv(6,i)*part%pv(6,i) - v2old)

    end do

  end subroutine apply_electron_heating


  subroutine load_gauss_local(vx, vy, vt, rnd)
    use mod_constants, only: pi
    real(real64), intent(out) :: vx, vy
    real(real64), intent(in)  :: vt
    real(real64), intent(in)  :: rnd(2)
    real(real64) :: theta, vp

    vp    = vt * sqrt(-log(1.0_real64 - rnd(1)))
    theta = 2.0_real64 * pi * rnd(2)
    vx    = vp * cos(theta)
    vy    = vp * sin(theta)
  end subroutine load_gauss_local

end module mod_heating