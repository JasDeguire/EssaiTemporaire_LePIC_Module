module mod_CoulombCollisions
  !=========================================================================
  ! Coulomb small-angle scattering (Nanbu 2000) ported from the standalone
  ! coulomb.f90 onto the modular ParticleSet flat-array architecture.
  !
  ! Reference:
  !   K. Nanbu & S. Yonemura, J. Comp. Phys. 145, 639 (1998)
  !   K. Nanbu, IEEE Trans. Plasma Sci. 28, 971 (2000)
  !
  ! Public entry point:
  !   perform_coulomb_collisions -- called from mod_simulation at the same
  !   cadence as MC collisions. Does three things per call:
  !     1. Compute species temperatures Ti(ptype) from global particle
  !        velocity moments (reduces over all iproc).
  !     2. Compute the Nanbu Coulomb parameter Aab and pair counts Nc
  !        for every (ptype,ttype) interaction, using the global Ti and
  !        the input electron density n_e.
  !     3. Run the OMP-parallel per-iproc collision loop: random permute
  !        of live particle indices, then Nanbu scattering for each pair.
  !
  ! Structural choices vs. the standalone:
  !   - Particles are NOT copied into a global vxp array. Velocities are
  !     read/written directly from part(ptype,iproc)%vx/vy/vz.
  !   - Collision pairs are formed within each iproc's own particle set
  !     (same OMP-thread semantics as the existing MC collision module).
  !     This is statistically equivalent to a fully mixed global pairing
  !     provided each iproc holds enough particles (>> 1).
  !   - Nc is computed locally per iproc from local particle counts, then
  !     corrected for the odd-particle edge case exactly as in the standalone.
  !   - The global Fisher-Yates permutation with !$OMP CRITICAL is replaced
  !     by a per-iproc in-place shuffle (no critical sections needed since
  !     each thread owns its own index array).
  !   - Ti is computed globally (reduction over all iproc, all species)
  !     before the OMP collision loop, so the Coulomb logarithm uses
  !     thermodynamically consistent temperatures rather than local
  !     per-thread estimates.
  !=========================================================================
  use iso_fortran_env, only: int8, int32, real64
  use mod_particles,   only: ParticleSet
  use mod_constants,   only: qe, eps0, pi
  use mod_RNG,         only: ran2
  implicit none
  private

  public :: perform_coulomb_collisions

  ! Maximum number of species for fixed-size local arrays inside the
  ! per-iproc collision loop (avoids heap allocation on the hot path).
  integer(int32), parameter :: NTYPE_MAX = 10_int32

contains

  !=========================================================================
  ! Top-level entry point.
  ! n_e    : global average electron number density [m^-3], used for the
  !          Debye length. Pass the domain-averaged np_red value for ptype=1,
  !          or cfg%ni0(1)*Nm(1) as a target density.
  ! Nm     : macroparticle weights, one per species (all same value in ITER)
  ! mass   : species masses [kg]
  ! charge : species charges [C]
  !=========================================================================
  subroutine perform_coulomb_collisions( &
      part, ntype, nproc, mass, charge, Nm, n_e, dt, iseed)

    type(ParticleSet), intent(inout) :: part(ntype, nproc)
    integer(int32),    intent(in)    :: ntype, nproc
    real(real64),      intent(in)    :: mass(ntype)
    real(real64),      intent(in)    :: charge(ntype)
    real(real64),      intent(in)    :: Nm(ntype)
    real(real64),      intent(in)    :: n_e      ! [m^-3] global electron density
    real(real64),      intent(in)    :: dt       ! [s]
    integer(int32),    intent(inout) :: iseed(nproc)

    integer(int32) :: ptype, ttype, iproc
    integer(int32) :: N_global(NTYPE_MAX)
    real(real64)   :: Ti(NTYPE_MAX)             ! kinetic temperature [eV]
    real(real64)   :: Aab(NTYPE_MAX, NTYPE_MAX) ! Nanbu Coulomb parameter [m^3/s^3 units]
    real(real64)   :: n_tot                     ! total number density [m^-3]

    if (ntype > NTYPE_MAX) then
      write(*,*) 'perform_coulomb_collisions: ntype=', ntype, &
                 ' exceeds NTYPE_MAX=', NTYPE_MAX, '; increase NTYPE_MAX.'
      error stop 'perform_coulomb_collisions: NTYPE_MAX too small'
    end if

    ! --- Step 1: global particle counts and temperatures ---
    call compute_global_moments(part, ntype, nproc, mass, N_global, Ti)

    ! --- Step 2: n_tot, Aab matrix ---
    n_tot = compute_n_tot(N_global, Nm, n_e, ntype)

    do ptype = 1_int32, ntype
      do ttype = 1_int32, ntype
        if (abs(charge(ptype)) < 1.0e-30_real64 .or. &
            abs(charge(ttype)) < 1.0e-30_real64) then
          Aab(ptype, ttype) = 0.0_real64
          cycle
        end if
        Aab(ptype, ttype) = compute_Aab( &
          charge(ptype), charge(ttype), mass(ptype), mass(ttype), &
          Ti(ptype), Ti(ttype), n_tot, n_e)
      end do
    end do

    ! --- Step 3: per-iproc collision loops ---
    !$omp parallel do private(iproc) schedule(static)
    do iproc = 1_int32, nproc
      call collide_one_thread( &
        part, ntype, nproc, iproc, mass, charge, Nm, &
        Ti, Aab, n_tot, dt, iseed(iproc))
    end do
    !$omp end parallel do

  end subroutine perform_coulomb_collisions


  !=========================================================================
  ! Compute global particle counts N_global(ptype) and kinetic temperatures
  ! Ti(ptype) [eV] from the velocity moments, reducing over all iproc.
  ! Dead particles (flag_dead /= 0) are excluded.
  !=========================================================================
  subroutine compute_global_moments(part, ntype, nproc, mass, N_global, Ti)
    type(ParticleSet), intent(in)  :: part(ntype, nproc)
    integer(int32),    intent(in)  :: ntype, nproc
    real(real64),      intent(in)  :: mass(ntype)
    integer(int32),    intent(out) :: N_global(NTYPE_MAX)
    real(real64),      intent(out) :: Ti(NTYPE_MAX)

    real(real64) :: v2sum(NTYPE_MAX)
    integer(int32) :: ptype, iproc, i, nlive

    N_global = 0_int32
    v2sum    = 0.0_real64

    do ptype = 1_int32, ntype
      do iproc = 1_int32, nproc
        if (.not. allocated(part(ptype,iproc)%pv)) cycle
        nlive = 0_int32
        do i = 1_int32, part(ptype,iproc)%n
          if (part(ptype,iproc)%flag_dead(i) /= 0_int8) cycle
          nlive = nlive + 1_int32
          v2sum(ptype) = v2sum(ptype) + &
            part(ptype,iproc)%pv(4,i)**2 + &
            part(ptype,iproc)%pv(5,i)**2 + &
            part(ptype,iproc)%pv(6,i)**2
        end do
        N_global(ptype) = N_global(ptype) + nlive
      end do
    end do

    ! Ti [eV] = m * <v^2> / (3 * qe)
    do ptype = 1_int32, ntype
      if (N_global(ptype) > 0_int32) then
        Ti(ptype) = mass(ptype) * (v2sum(ptype) / real(N_global(ptype), real64)) &
                    / (3.0_real64 * qe)
      else
        Ti(ptype) = 0.0_real64
      end if
    end do
  end subroutine compute_global_moments


  !=========================================================================
  ! Compute n_tot = sum of all species number densities [m^-3].
  ! n_e and N_global(1) give the electron Nm-weighted density; the other
  ! species are scaled by their own Nm weight.
  ! If all species share the same Nm, this simplifies to n_e * Ntot/Ne.
  !=========================================================================
  real(real64) function compute_n_tot(N_global, Nm, n_e, ntype)
    integer(int32), intent(in) :: N_global(NTYPE_MAX)
    real(real64),   intent(in) :: Nm(ntype)
    real(real64),   intent(in) :: n_e
    integer(int32), intent(in) :: ntype

    integer(int32) :: ptype
    real(real64)   :: Ne_weighted

    Ne_weighted = real(N_global(1), real64) * Nm(1)
    compute_n_tot = 0.0_real64
    if (Ne_weighted <= 0.0_real64) return

    do ptype = 1_int32, ntype
      compute_n_tot = compute_n_tot + &
        real(N_global(ptype), real64) * Nm(ptype) / Ne_weighted * n_e
    end do
  end function compute_n_tot


  !=========================================================================
  ! Compute Aab = (n_tot/(4*pi)) * (q_a*q_b/(eps0*mu))^2 * ln(Lambda)
  ! with the Debye length computed from n_e and Ti(1) (electrons).
  !=========================================================================
  real(real64) function compute_Aab( &
      qa, qb, ma, mb, Ta_eV, Tb_eV, n_tot, n_e)
    real(real64), intent(in) :: qa, qb, ma, mb
    real(real64), intent(in) :: Ta_eV, Tb_eV  ! [eV]
    real(real64), intent(in) :: n_tot, n_e

    real(real64) :: mu, g2avg, b_min, lbde, CoulombLog

    mu    = ma * mb / (ma + mb)
    g2avg = 3.0_real64 * qe * (Ta_eV/ma + Tb_eV/mb)

    if (g2avg <= 0.0_real64 .or. n_e <= 0.0_real64) then
      compute_Aab = 0.0_real64
      return
    end if

    b_min       = abs(qa*qb) / (4.0_real64 * pi * eps0 * mu * g2avg)
    lbde        = sqrt((eps0 * Ta_eV) / (n_e * abs(qa)))
    CoulombLog  = log(lbde / b_min)
    if (CoulombLog <= 1.0_real64) CoulombLog = 1.0_real64  ! floor

    compute_Aab = (n_tot / (4.0_real64 * pi)) * &
                  (qa * qb / (eps0 * mu))**2 * CoulombLog
  end function compute_Aab


  !=========================================================================
  ! Per-iproc collision loop, called from inside !$omp parallel do.
  ! Mirrors coulomb_collisions() in the standalone but operates directly
  ! on part(ptype,iproc) flat arrays with no global vxp copy.
  !=========================================================================
  subroutine collide_one_thread( &
      part, ntype, nproc, iproc, mass, charge, Nm, &
      Ti, Aab, n_tot, dt, iseed)

    type(ParticleSet), intent(inout) :: part(ntype, nproc)
    integer(int32),    intent(in)    :: ntype, nproc, iproc
    real(real64),      intent(in)    :: mass(ntype), charge(ntype), Nm(ntype)
    real(real64),      intent(in)    :: Ti(NTYPE_MAX)
    real(real64),      intent(in)    :: Aab(NTYPE_MAX, NTYPE_MAX)
    real(real64),      intent(in)    :: n_tot, dt
    integer(int32),    intent(inout) :: iseed

    ! Local arrays sized to NTYPE_MAX for stack allocation
    integer(int32) :: N_local(NTYPE_MAX)
    integer(int32) :: N_local_tot
    integer(int32) :: Nc_local(NTYPE_MAX, NTYPE_MAX)
    integer(int32) :: sum_Nc(NTYPE_MAX)
    integer(int32), allocatable :: ind(:,:)  ! ind(max_n, ntype): live particle indices

    integer(int32) :: ptype, ttype, ip, ii, itg, max_n
    integer(int32) :: Nc_pt, imin, imax, flag_odd
    real(real64)   :: dtp, rnd

    ! --- Count live particles per species in this iproc ---
    N_local = 0_int32
    do ptype = 1_int32, ntype
      if (.not. allocated(part(ptype,iproc)%pv)) cycle
      N_local(ptype) = count_live(part(ptype,iproc))
    end do
    N_local_tot = sum(N_local(1:ntype))
    if (N_local_tot <= 1_int32) return

    ! --- Build per-species live-particle index arrays ---
    max_n = maxval(N_local(1:ntype))
    if (max_n <= 0_int32) return
    allocate(ind(max_n, ntype))
    ind = 0_int32
    do ptype = 1_int32, ntype
      if (N_local(ptype) <= 0_int32) cycle
      call build_live_index(part(ptype,iproc), N_local(ptype), ind(:,ptype))
    end do

    ! --- Compute local Nc and correct for odd numbers (Nanbu 2000 §III) ---
    Nc_local = 0_int32
    sum_Nc   = 0_int32
    do ptype = 1_int32, ntype
      do ttype = 1_int32, ntype
        if (Nc_local(ttype,ptype) == 0_int32) then
          if (N_local_tot > 0_int32) then
            Nc_local(ptype,ttype) = &
              int(real(N_local(ptype),real64) * real(N_local(ttype),real64) &
                  / real(N_local_tot, real64), int32)
          end if
        end if
        sum_Nc(ptype) = sum_Nc(ptype) + max(Nc_local(ptype,ttype), Nc_local(ttype,ptype))
      end do
    end do

    do ptype = 1_int32, ntype
      if (sum_Nc(ptype) /= N_local(ptype)) then
        Nc_local(ptype,ptype) = N_local(ptype) - sum_Nc(ptype) + Nc_local(ptype,ptype)
      end if
    end do

    ! --- Random permutation of each species' index array ---
    do ptype = 1_int32, ntype
      if (N_local(ptype) > 1_int32) &
        call fisher_yates(ind(1:N_local(ptype), ptype), N_local(ptype), iseed)
    end do

    ! --- Collision loops ---
    sum_Nc = 0_int32
    do ptype = 1_int32, ntype
      do ttype = 1_int32, ntype

        Nc_pt = Nc_local(ptype, ttype)
        if (Nc_pt <= 0_int32) then
          sum_Nc(ptype) = sum_Nc(ptype) + max(Nc_local(ptype,ttype), Nc_local(ttype,ptype))
          cycle
        end if
        if (Aab(ptype,ttype) == 0.0_real64) then
          sum_Nc(ptype) = sum_Nc(ptype) + max(Nc_local(ptype,ttype), Nc_local(ttype,ptype))
          cycle
        end if

        imin = sum_Nc(ptype) + 1_int32
        imax = sum_Nc(ptype) + Nc_pt

        dtp      = dt
        flag_odd = 0_int32

        if (ptype == ttype) then
          ! Odd-particle correction (like-particle pairs)
          if (mod(imax - sum_Nc(ptype), 2_int32) /= 0_int32) then
            imax     = imax + 1_int32
            dtp      = real(Nc_pt, real64) / real(Nc_pt + 1_int32, real64) * dt
            flag_odd = 1_int32
          end if
          imax = sum_Nc(ptype) + (imax - sum_Nc(ptype)) / 2_int32
        end if

        do ip = imin, imax

          if (ip > N_local(ptype) .or. ip > max_n) exit

          ii = ind(ip, ptype)

          if (ptype /= ttype) then
            if (ip > N_local(ttype)) exit
            itg = ind(ip, ttype)
          else
            ! Like-particle: second half of the shuffled index array
            if (ip + imax > N_local(ptype)) exit
            itg = ind(ip + imax - sum_Nc(ptype), ttype)
            ! Particle #1 collides twice when N is odd
            if (flag_odd == 1_int32 .and. ip == imax) then
              rnd = ran2(iseed)
              itg = ind(max(1_int32, int(rnd * real(imax - sum_Nc(ptype) - 1_int32, real64), int32) + 1_int32), ttype)
            end if
          end if

          if (ii <= 0_int32 .or. itg <= 0_int32) cycle
          if (ii > part(ptype,iproc)%n .or. itg > part(ttype,iproc)%n) cycle

          call nanbu_scatter( &
            part(ptype,iproc), part(ttype,iproc), &
            ii, itg, mass(ptype), mass(ttype), Aab(ptype,ttype), dtp, iseed)

        end do

        sum_Nc(ptype) = sum_Nc(ptype) + max(Nc_local(ptype,ttype), Nc_local(ttype,ptype))
      end do
    end do

    deallocate(ind)
  end subroutine collide_one_thread


  !=========================================================================
  ! Nanbu scattering kernel (coulomb_scattering in standalone).
  ! Updates vx/vy/vz of particles ii (ptype) and itg (ttype) in-place.
  !=========================================================================
  subroutine nanbu_scatter(pa, pb, ii, itg, ma, mb, Aab, dt, iseed)
    type(ParticleSet), intent(inout) :: pa, pb
    integer(int32),    intent(in)    :: ii, itg
    real(real64),      intent(in)    :: ma, mb, Aab, dt
    integer(int32),    intent(inout) :: iseed

    real(real64) :: g(3), h(3)
    real(real64) :: g_perp, g2_perp, g_norm, phi, tau
    real(real64) :: A_tau, cos_chi, sin_chi
    real(real64) :: mass_ratio_a, mass_ratio_b
    real(real64) :: rnd, dv(3)

    g(1) = pa%pv(4,ii) - pb%pv(4,itg)
    g(2) = pa%pv(5,ii) - pb%pv(5,itg)
    g(3) = pa%pv(6,ii) - pb%pv(6,itg)

    g2_perp = g(2)**2 + g(3)**2
    g_norm  = sqrt(g(1)**2 + g2_perp)
    if (g_norm < 1.0e-30_real64) return  ! no relative motion, skip

    g_perp  = sqrt(g2_perp)
    tau     = Aab * dt / g_norm**3

    rnd = ran2(iseed)
    phi = 2.0_real64 * pi * rnd

    if (g_perp > 1.0e-30_real64) then
      h(1) =  g_perp * cos(phi)
      h(2) = -(g(1)*g(2)*cos(phi) + g_norm*g(3)*sin(phi)) / g_perp
      h(3) = -(g(1)*g(3)*cos(phi) - g_norm*g(2)*sin(phi)) / g_perp
    else
      ! g is along x: use degenerate h
      h(1) = 0.0_real64
      h(2) = g_norm * cos(phi)
      h(3) = g_norm * sin(phi)
    end if

    rnd = ran2(iseed)
    cos_chi = sample_cos_chi(tau, rnd)
    sin_chi = sqrt(max(0.0_real64, 1.0_real64 - cos_chi**2))

    mass_ratio_b = mb / (ma + mb)
    mass_ratio_a = ma / (ma + mb)

    dv(1) = g(1)*(1.0_real64 - cos_chi) + h(1)*sin_chi
    dv(2) = g(2)*(1.0_real64 - cos_chi) + h(2)*sin_chi
    dv(3) = g(3)*(1.0_real64 - cos_chi) + h(3)*sin_chi

    pa%pv(4,ii) = pa%pv(4,ii) - mass_ratio_b * dv(1)
    pa%pv(5,ii) = pa%pv(5,ii) - mass_ratio_b * dv(2)
    pa%pv(6,ii) = pa%pv(6,ii) - mass_ratio_b * dv(3)

    pb%pv(4,itg) = pb%pv(4,itg) + mass_ratio_a * dv(1)
    pb%pv(5,itg) = pb%pv(5,itg) + mass_ratio_a * dv(2)
    pb%pv(6,itg) = pb%pv(6,itg) + mass_ratio_a * dv(3)

  end subroutine nanbu_scatter


  !=========================================================================
  ! Sample cos(chi) from the Nanbu distribution for a given tau = Aab*dt/g^3.
  ! Three regimes as in the standalone (Nanbu 2000 eq. 13):
  !   tau <= 0.1   : small-angle approx  cos_chi = max(1 + tau*ln(U), -1)
  !   0.1 < tau <= 6: full solution via Newton-Raphson for A(tau)
  !   tau > 6      : isotropic          cos_chi = 2U - 1
  !=========================================================================
  real(real64) function sample_cos_chi(tau, u)
    real(real64), intent(in) :: tau, u   ! u ~ Uniform(0,1)

    real(real64) :: A_tau

    if (tau <= 0.1_real64) then
      sample_cos_chi = max(1.0_real64 + tau * log(u), -1.0_real64)
    else if (tau <= 6.0_real64) then
      A_tau = newton_raphson_A(tau)
      sample_cos_chi = (A_tau + log(u + (1.0_real64 - u) * exp(-2.0_real64*A_tau))) / A_tau
    else
      sample_cos_chi = 2.0_real64*u - 1.0_real64
    end if
  end function sample_cos_chi


  !=========================================================================
  ! Newton-Raphson solver for A(tau): coth(A) - 1/A = exp(-tau).
  ! Direct port of NewtRaph() from the standalone.
  !=========================================================================
  real(real64) function newton_raphson_A(tau)
    real(real64), intent(in) :: tau

    real(real64), parameter :: tol = 1.0e-8_real64
    integer(int32), parameter :: max_iter = 50_int32
    real(real64) :: A, f_val, df_val, et
    integer(int32) :: i

    A  = 10.0_real64   ! initial guess — works well for tau in [0.1, 6]
    et = exp(-tau)

    do i = 1_int32, max_iter
      f_val  = 1.0_real64/tanh(A) - 1.0_real64/A - et
      df_val = -1.0_real64/sinh(A)**2 + 1.0_real64/A**2
      if (abs(f_val) < tol) exit
      if (abs(df_val) < 1.0e-30_real64) exit
      A = A - f_val / df_val
      if (A <= 0.0_real64) A = 5.0_real64   ! guard against negative A
    end do

    newton_raphson_A = A
  end function newton_raphson_A


  !=========================================================================
  ! Count live particles (flag_dead == 0) in a ParticleSet.
  !=========================================================================
  integer(int32) function count_live(p)
    type(ParticleSet), intent(in) :: p
    integer(int32) :: i
    count_live = 0_int32
    if (.not. allocated(p%pv)) return
    do i = 1_int32, p%n
      if (p%flag_dead(i) == 0_int8) count_live = count_live + 1_int32
    end do
  end function count_live


  !=========================================================================
  ! Build an index array of the indices of live particles in p.
  ! ind(1:nlive) on exit contains the original array indices of live particles.
  !=========================================================================
  subroutine build_live_index(p, nlive, ind)
    type(ParticleSet), intent(in)  :: p
    integer(int32),    intent(in)  :: nlive
    integer(int32),    intent(out) :: ind(nlive)
    integer(int32) :: i, k
    k = 0_int32
    do i = 1_int32, p%n
      if (p%flag_dead(i) /= 0_int8) cycle
      k = k + 1_int32
      if (k > nlive) exit
      ind(k) = i
    end do
  end subroutine build_live_index


  !=========================================================================
  ! Fisher-Yates in-place shuffle of integer array a(1:n).
  ! No critical sections — called per-thread on a thread-private array.
  !=========================================================================
  subroutine fisher_yates(a, n, iseed)
    integer(int32), intent(inout) :: a(n)
    integer(int32), intent(in)    :: n
    integer(int32), intent(inout) :: iseed
    integer(int32) :: i, j, tmp
    real(real64)   :: rnd
    do i = n, 2_int32, -1_int32
      rnd = ran2(iseed)
      j   = 1_int32 + int(rnd * real(i, real64), int32)
      j   = min(j, i)
      tmp  = a(i)
      a(i) = a(j)
      a(j) = tmp
    end do
  end subroutine fisher_yates

end module mod_CoulombCollisions
