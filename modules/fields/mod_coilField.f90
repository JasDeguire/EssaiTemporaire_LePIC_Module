module mod_coilField
  ! RF mode 'Coil' (cfg%flag_coil_ant==1): prescribed inductive field of the
  ! real planar coil, Phase 1 of "A better model for the planar antenna".
  !
  ! Each turn j (radius a_j, axial position x_j, relative current w_j) is a
  ! circular loop coaxial with the domain's x axis through the y-z centre.
  ! Its vector potential per ampere is the closed form (Jackson 3rd ed. 5.37)
  !
  !   A_loop(r,x) = mu0/(pi k) sqrt(a/r) [ (1 - k^2/2) K(k) - E(k) ],
  !   k^2 = 4 a r / ((a+r)^2 + (x-x_j)^2),
  !
  ! summed once over the turns on a 2D (r,x) table at startup. The table
  ! holds the vacuum A_hat, and B_hat = curl(A_hat theta_hat) by central
  ! differences of psi = r*A_hat:
  !
  !   B_hat_r = -(1/r) d(psi)/dx,     B_hat_x = (1/r) d(psi)/dr,
  !
  ! so the discrete div B (same central operators) vanishes to round-off.
  ! Plasma shielding is the assumed factor g(x) = exp(-(x-x0)/delta) past the
  ! window x0 (delta = State%gams, updated every RF period), applied at
  ! evaluation time so the table never needs rebuilding:
  !
  !   A = g A_hat,  B_r = g (B_hat_r + A_hat/delta),  B_x = g B_hat_x.
  !
  ! The mover then uses E_theta = omega I A sin(wt), B_RF = I B cos(wt)
  ! (A = I A_hat cos(wt), E = -dA/dt), with I = cfg%E0_RF rescaled toward
  ! Pabs by update_rf_convergence exactly as E0_RF is for the other modes.
  use iso_fortran_env, only: int32, real64
  use mod_constants,   only: pi
  implicit none
  private

  public :: CoilField

  real(real64), parameter :: mu0 = 4.0e-7_real64*pi

  type :: CoilField
    integer(int32) :: nturn = 0_int32
    real(real64), allocatable :: a(:), xc(:), w(:)
    ! Table nodes r_j = j*dr (j=0..nr), x_i = i*dx (i=0..nx), per ampere
    ! of coil current: A_hat (T m/A), B_hat_r and B_hat_x (T/A).
    integer(int32) :: nr = 0_int32, nx = 0_int32
    real(real64)   :: dr = 0.0_real64, dx = 0.0_real64
    real(real64), allocatable :: Ahat(:,:), Brhat(:,:), Bxhat(:,:)
  contains
    procedure :: build
    procedure :: eval
    procedure :: self_check
    procedure :: write_table
  end type CoilField

contains

  pure subroutine ellip_KE(m, K, E)
    ! Complete elliptic integrals of the first and second kind, parameter
    ! m = k^2 in [0,1), by the arithmetic-geometric mean.
    real(real64), intent(in)  :: m
    real(real64), intent(out) :: K, E
    real(real64) :: a, b, c, an, p, s
    integer :: it

    a = 1.0_real64
    b = sqrt(1.0_real64 - m)
    s = 0.5_real64*m
    p = 0.5_real64
    do it = 1, 40
      if (abs(a-b) <= 1.0e-15_real64*a) exit
      c  = 0.5_real64*(a - b)
      an = 0.5_real64*(a + b)
      b  = sqrt(a*b)
      a  = an
      p  = 2.0_real64*p
      s  = s + p*c*c
    end do
    K = pi/(2.0_real64*a)
    E = K*(1.0_real64 - s)
  end subroutine ellip_KE


  pure real(real64) function A_loop(r, x, a, xc) result(Aval)
    ! Vector potential (T m) of a unit-current loop of radius a at x = xc.
    real(real64), intent(in) :: r, x, a, xc
    real(real64) :: m, f, K, E

    if (r <= 0.0_real64) then
      Aval = 0.0_real64
      return
    end if
    m = 4.0_real64*a*r / ((a + r)**2 + (x - xc)**2)
    m = min(m, 1.0_real64 - 1.0e-14_real64)
    if (m < 1.0e-3_real64) then
      ! (1-m/2)K - E cancels to O(m^2): use its series instead.
      f = 0.5_real64*pi * m*m/16.0_real64 * &
          (1.0_real64 + 0.75_real64*m + 75.0_real64/128.0_real64*m*m)
    else
      call ellip_KE(m, K, E)
      f = (1.0_real64 - 0.5_real64*m)*K - E
    end if
    Aval = mu0/(pi*sqrt(m)) * sqrt(a/r) * f
  end function A_loop


  subroutine build(self, a, xc, w, h, xmax, ymax, zmax)
    ! Tabulate A_hat, B_hat on (r,x) at half the grid spacing, covering the
    ! whole box (r up to the y-z corner). psi gets one extra node on every
    ! side so all table nodes use central differences.
    class(CoilField), intent(inout) :: self
    real(real64),     intent(in)    :: a(:), xc(:), w(:)
    real(real64),     intent(in)    :: h(3), xmax, ymax, zmax

    real(real64), allocatable :: psi(:,:)
    real(real64) :: r, x, rmax
    integer(int32) :: i, j, t

    self%nturn = int(size(a), int32)
    self%a  = a
    self%xc = xc
    self%w  = w

    self%dr = 0.5_real64*min(h(2), h(3))
    self%dx = 0.5_real64*h(1)
    rmax    = sqrt((0.5_real64*ymax)**2 + (0.5_real64*zmax)**2)
    self%nr = int(ceiling(rmax/self%dr), int32) + 1_int32
    self%nx = int(ceiling(xmax/self%dx), int32) + 1_int32

    allocate(psi(0:self%nr+1, -1:self%nx+1))
    !$omp parallel do private(i,j,t,r,x) schedule(static)
    do i = -1, self%nx+1
      x = real(i, real64)*self%dx
      do j = 0, self%nr+1
        r = real(j, real64)*self%dr
        psi(j,i) = 0.0_real64
        do t = 1, self%nturn
          psi(j,i) = psi(j,i) + self%w(t)*r*A_loop(r, x, self%a(t), self%xc(t))
        end do
      end do
    end do
    !$omp end parallel do

    if (allocated(self%Ahat)) deallocate(self%Ahat, self%Brhat, self%Bxhat)
    allocate(self%Ahat (0:self%nr, 0:self%nx), &
             self%Brhat(0:self%nr, 0:self%nx), &
             self%Bxhat(0:self%nr, 0:self%nx))
    do i = 0, self%nx
      ! Axis: A and B_r vanish, B_x = lim (1/r) d(psi)/dr = 2 psi(dr)/dr^2.
      self%Ahat (0,i) = 0.0_real64
      self%Brhat(0,i) = 0.0_real64
      self%Bxhat(0,i) = 2.0_real64*psi(1,i)/self%dr**2
      do j = 1, self%nr
        r = real(j, real64)*self%dr
        self%Ahat (j,i) = psi(j,i)/r
        self%Brhat(j,i) = -(psi(j,i+1) - psi(j,i-1))/(2.0_real64*self%dx*r)
        self%Bxhat(j,i) =  (psi(j+1,i) - psi(j-1,i))/(2.0_real64*self%dr*r)
      end do
    end do
    deallocate(psi)
  end subroutine build


  pure subroutine eval(self, xp, rp, x0, delta, A, Br, Bx)
    ! Bilinear (r,x) interpolation of the table, shielding g(x) applied -
    ! per ampere of coil current, amplitudes only (no time factor).
    class(CoilField), intent(in)  :: self
    real(real64),     intent(in)  :: xp, rp, x0, delta
    real(real64),     intent(out) :: A, Br, Bx

    integer(int32) :: j, i
    real(real64)   :: tr, tx, w00, w10, w01, w11, g

    tr = rp/self%dr
    tx = xp/self%dx
    j  = min(max(int(tr, int32), 0_int32), self%nr-1_int32)
    i  = min(max(int(tx, int32), 0_int32), self%nx-1_int32)
    tr = min(max(tr - real(j, real64), 0.0_real64), 1.0_real64)
    tx = min(max(tx - real(i, real64), 0.0_real64), 1.0_real64)

    w00 = (1.0_real64-tr)*(1.0_real64-tx)
    w10 = tr*(1.0_real64-tx)
    w01 = (1.0_real64-tr)*tx
    w11 = tr*tx

    A  = w00*self%Ahat (j,i) + w10*self%Ahat (j+1,i) + w01*self%Ahat (j,i+1) + w11*self%Ahat (j+1,i+1)
    Br = w00*self%Brhat(j,i) + w10*self%Brhat(j+1,i) + w01*self%Brhat(j,i+1) + w11*self%Brhat(j+1,i+1)
    Bx = w00*self%Bxhat(j,i) + w10*self%Bxhat(j+1,i) + w01*self%Bxhat(j,i+1) + w11*self%Bxhat(j+1,i+1)

    if (xp > x0) then
      g  = exp(-(xp - x0)/delta)
      Br = g*(Br + A/delta)
      Bx = g*Bx
      A  = g*A
    end if
  end subroutine eval


  subroutine self_check(self, mpi_rank)
    ! Vacuum checks of the table (slide "Checks maths"), printed on rank 0:
    !   1. on-axis B_x vs the exact sum of mu0 w a^2 / (2 (a^2+(x-x_j)^2)^1.5)
    !   2. far field of the loop formula vs a magnetic dipole, at R = 100 a
    !   3. discrete div B with the same central operators, to round-off
    class(CoilField), intent(in) :: self
    integer,          intent(in) :: mpi_rank

    real(real64) :: x, Bex, err_ax, Bmax_ax, Rfar, th, Afar, Adip, err_dip
    real(real64) :: r, divB, div_max, B_max
    integer(int32) :: i, j, t

    if (mpi_rank /= 0) return

    err_ax  = 0.0_real64
    Bmax_ax = 0.0_real64
    do i = 0, self%nx
      x   = real(i, real64)*self%dx
      Bex = 0.0_real64
      do t = 1, self%nturn
        Bex = Bex + self%w(t)*mu0*self%a(t)**2 / &
              (2.0_real64*(self%a(t)**2 + (x - self%xc(t))**2)**1.5_real64)
      end do
      err_ax  = max(err_ax, abs(self%Bxhat(0,i) - Bex))
      Bmax_ax = max(Bmax_ax, abs(Bex))
    end do
    err_ax = err_ax/max(Bmax_ax, tiny(1.0_real64))

    Rfar = 100.0_real64*self%a(1)
    th   = 0.25_real64*pi
    Afar = A_loop(Rfar*sin(th), self%xc(1) + Rfar*cos(th), self%a(1), self%xc(1))
    Adip = mu0*self%a(1)**2*sin(th)/(4.0_real64*Rfar**2)
    err_dip = abs(Afar/Adip - 1.0_real64)

    div_max = 0.0_real64
    B_max   = maxval(abs(self%Bxhat)) + maxval(abs(self%Brhat))
    do i = 1, self%nx-1
      do j = 1, self%nr-1
        r    = real(j, real64)*self%dr
        divB = (real(j+1, real64)*self%dr*self%Brhat(j+1,i) - &
                real(j-1, real64)*self%dr*self%Brhat(j-1,i))/(2.0_real64*self%dr*r) + &
               (self%Bxhat(j,i+1) - self%Bxhat(j,i-1))/(2.0_real64*self%dx)
        div_max = max(div_max, abs(divB))
      end do
    end do
    div_max = div_max*min(self%dr, self%dx)/max(B_max, tiny(1.0_real64))

    write(*,'(1x,a,i0,a,f8.4,a,f8.4,a,i0,a,i0)') 'RF antenna: Coil, ', self%nturn, &
      ' turns, table dr(cm)=', self%dr*1.0e2_real64, ', dx(cm)=', self%dx*1.0e2_real64, &
      ', nr=', self%nr, ', nx=', self%nx
    write(*,'(1x,a,es10.3,a,es10.3)') 'Coil vacuum field: B_x on axis per ampere (T/A) max=', &
      Bmax_ax, ', rel. error vs exact=', err_ax
    write(*,'(1x,a,es10.3,a)') 'Coil vacuum field: loop A vs dipole at R=100a, rel. diff=', &
      err_dip, ' (expected ~1e-4)'
    write(*,'(1x,a,es10.3)') 'Coil vacuum field: max |div B| h/|B|max=', div_max
  end subroutine self_check


  subroutine write_table(self, fname, mpi_rank)
    ! Vacuum table (no shielding) for plotting: one line per (r,x) node.
    class(CoilField), intent(in) :: self
    character(len=*), intent(in) :: fname
    integer,          intent(in) :: mpi_rank
    integer :: u
    integer(int32) :: i, j

    if (mpi_rank /= 0) return
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(a)') '# r(cm), x(cm), A_hat(T m/A), B_hat_r(T/A), B_hat_x(T/A) - vacuum, per ampere'
    do i = 0, self%nx
      do j = 0, self%nr
        write(u,'(5(1x,es14.6))') real(j, real64)*self%dr*1.0e2_real64, &
          real(i, real64)*self%dx*1.0e2_real64, self%Ahat(j,i), self%Brhat(j,i), self%Bxhat(j,i)
      end do
      write(u,*)
    end do
    close(u)
  end subroutine write_table

end module mod_coilField
