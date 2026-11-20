module mod_rng

  use iso_fortran_env,  only: int32, real64
  use mod_constants,    only: pi
  implicit none
  private

  public :: ran2
  public :: seed_initialization
  public :: load_gauss

contains


  subroutine seed_initialization(iseed, mpi_rank, iproc)
    integer(int32), intent(out) :: iseed
    integer,        intent(in)  :: mpi_rank
    integer(int32), intent(in)  :: iproc

    ! EXACT legacy formula:
    iseed = 123456_int32 * iproc * int(10*mpi_rank + 1, int32)

    if (iseed == 0_int32) iseed = 1_int32
  end subroutine seed_initialization


  !=========================================================
  ! Legacy ran2 RNG (exact same as original)
  !=========================================================
  real(real64) function ran2(irand)
    integer(int32), intent(inout) :: irand

    integer(int32), parameter :: ia = 16807_int32
    integer(int32), parameter :: im = 2147483647_int32
    integer(int32), parameter :: iq = 127773_int32
    integer(int32), parameter :: ir = 2836_int32
    real(real64),   parameter :: am = 1.0_real64 / real(im, real64)

    integer(int32) :: k

    k = irand / iq
    irand = ia * (irand - k * iq) - ir * k
    if (irand < 0_int32) irand = irand + im

    ran2 = am * real(irand, real64)

  end function ran2


  ! Gaussian velocity sampling in the transverse plane.
  ! Given two uniform random numbers rnd(1:2) in (0,1), returns
  ! two velocity components (vx, vy) from a Maxwellian at thermal speed vt.
  subroutine load_gauss(vx, vy, vt, rnd)
    real(real64), intent(out) :: vx, vy
    real(real64), intent(in)  :: vt
    real(real64), intent(in)  :: rnd(2)
    real(real64) :: vp, theta

    vp    = vt * sqrt(-log(1.0_real64 - rnd(1)))
    theta = 2.0_real64 * pi * rnd(2)
    vx    = vp * cos(theta)
    vy    = vp * sin(theta)
  end subroutine load_gauss

end module mod_rng