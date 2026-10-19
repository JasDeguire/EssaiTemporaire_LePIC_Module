module mod_injection
  !=========================================================================
  ! Port of legacy Src/part_injection.f90 and Src/part_flux_injection.f90
  ! onto the modular ParticleSet flat-array architecture.
  !
  ! Two public entry points:
  !
  !   inject_particles_volume:
  !     Port of part_injection(). Volume injection of electron+ion pairs.
  !     Called every step (ns_inj=1 hardcoded in legacy, same here) when
  !     flag_inj==1. Runs INSIDE the !$omp parallel do over iproc — each
  !     thread injects into its own part(ptype,iproc), exactly as legacy
  !     calls part_injection() inside the !$OMP PARALLEL block.
  !     Supports all four legacy injection modes via opt_inj:
  !       +/-1 : uniform (opt>0) or cosine (opt<0) spatial distribution
  !        +3  : electron beam (shifted Maxwellian flux, disk geometry)
  !        +4  : cathode emission (half-Maxwellian, flat disk, +/-Z)
  !        +/-2: N_inj counter driven (particle count = ion wall losses)
  !
  !   inject_flux_particles:
  !     Port of part_flux_injection() + load_flux_OMP(). Surface flux
  !     injection of negative ions off the extraction electrode (jne>0).
  !     Called outside the OMP region, fires its own internal !$omp parallel
  !     do, every step (ns_flx=1, same as legacy). Only the jne>0 path is
  !     implemented (flag_thr thruster path is present but guarded; the
  !     user confirmed only jne>0 is used).
  !
  ! Source diagnostic arrays (sour_xy, sour_xz, sour_yz, sour_fx_yz) are
  ! passed in from State, where they are allocated — see note in
  ! mod_state.f90 about their allocation.
  !=========================================================================
  use iso_fortran_env,  only: int8, int32, int64, real64
  use mod_particles,    only: ParticleSet
  use mod_config,       only: Config
  use mod_RNG,          only: ran2
  use mod_constants,    only: pi, qe
  use omp_lib
  implicit none
  private

  public :: inject_particles_volume
  public :: inject_flux_particles

contains

  !=========================================================================
  ! Volume injection — port of part_injection()
  ! Must be called from INSIDE an !$omp parallel do over iproc.
  ! Each call handles one (iproc) thread's injection quota.
  !=========================================================================
  subroutine inject_particles_volume( &
      part, iproc, cfg, n, h, bcnd, phi, &
      vt0, Nm, mass, ni0, charge, &
      ntype, nproc, dt, &
      zg_sec, n_cath, dir_sec_inout, &
      tag_neg, tag_beam, &
      flag_pbc, flag_pbcz, &
      iseed, N_inj, &
      sour_xy, sour_xz, sour_yz, &
      iz_pl, ix_pl, &
      P_loss, mom_loss )

    type(ParticleSet), intent(inout) :: part(ntype, nproc)
    integer(int32),    intent(in)    :: iproc
    type(Config),      intent(in)    :: cfg
    integer(int32),    intent(in)    :: n(3)
    real(real64),      intent(in)    :: h(3)
    integer(int32),    intent(in)    :: bcnd(0:n(1)+2, 0:n(2)+2, 0:n(3)+2)
    real(real64),      intent(in)    :: phi(0:n(1)+2, 0:n(2)+2, 0:n(3)+2)
    real(real64),      intent(in)    :: vt0(ntype)
    real(real64),      intent(in)    :: Nm(ntype)
    real(real64),      intent(in)    :: mass(ntype)
    real(real64),      intent(in)    :: ni0(ntype)
    real(real64),      intent(in)    :: charge(ntype)
    integer(int32),    intent(in)    :: ntype, nproc
    real(real64),      intent(in)    :: dt
    real(real64),      intent(in)    :: zg_sec(2)
    integer(int32),    intent(in)    :: n_cath
    integer(int32),    intent(inout) :: dir_sec_inout
    integer(int32),    intent(in)    :: tag_neg, tag_beam
    integer(int32),    intent(in)    :: flag_pbc, flag_pbcz
    integer(int32),    intent(inout) :: iseed
    integer(int32),    intent(inout) :: N_inj(ntype, nproc)
    integer(int32),    intent(inout) :: sour_xy(0:n(1)+2, 0:n(2)+2, ntype, nproc)
    integer(int32),    intent(inout) :: sour_xz(0:n(1)+2, 0:n(3)+2, ntype, nproc)
    integer(int32),    intent(inout) :: sour_yz(0:n(2)+2, 0:n(3)+2, ntype, nproc)
    integer(int32),    intent(in)    :: iz_pl, ix_pl
    real(real64),      intent(inout) :: P_loss(4, ntype, nproc)
    ! Momentum-conservation diagnostic - vector counterpart of P_loss(4,:,:),
    ! filled alongside the kinetic-energy term below (no momentum
    ! equivalent for the e-beam potential-energy term further down -
    ! potential energy has no associated momentum vector).
    real(real64),      intent(inout) :: mom_loss(3, 4, ntype, nproc)

    integer(int32) :: i, k, ptype, N_inj_tmp, ng_sec, ix, iy, iz
    real(real64)   :: rnd(3), x, y, z, vx, vy, vz, vt, vz_sav(ntype)
    real(real64)   :: xm, ym, zm, dx, dy, dz, Rb, vb, vmax, fmax, dN_inj
    real(real64)   :: z_tmp, dt_tmp, phip
    real(real64)   :: px, py, pz, ki(8)
    integer(int32) :: mpi_size_1  ! nproc_mpi, hardcoded 1 for current single-MPI runs

    mpi_size_1 = 1_int32  ! adjust if MPI is enabled

    vz_sav = 0.0_real64
    xm = (cfg%xl_pow + cfg%xr_pow) / 2.0_real64
    dx = cfg%xr_pow - cfg%xl_pow
    ym = (cfg%yl_pow + cfg%yr_pow) / 2.0_real64
    dy = cfg%yr_pow - cfg%yl_pow
    zm = (cfg%zl_pow + cfg%zr_pow) / 2.0_real64
    dz = cfg%zr_pow - cfg%zl_pow
    Rb = min(dx, dy) / 2.0_real64

    ng_sec = 1_int32
    if (dir_sec_inout == -1_int32) ng_sec = 2_int32

    vb   = 0.0_real64
    fmax = 0.0_real64
    if (abs(cfg%opt_inj) == 3_int32 .or. abs(cfg%opt_inj) == 4_int32) then
      vb   = sqrt(2.0_real64 * qe * cfg%THm / abs(mass(1)))
      vt   = vt0(1)
      vmax = (vb + sqrt(vb**2 + 2.0_real64*vt**2)) / 2.0_real64
      fmax = vmax * exp(-(vmax - vb)**2 / vt**2)
    end if

    if (abs(cfg%opt_inj) /= 2_int32) then
      dN_inj   = cfg%I_inj * dt / (qe * Nm(1)) / real(mpi_size_1,real64) / real(nproc,real64)
      N_inj_tmp = int(dN_inj, int32)
      rnd(1) = ran2(iseed)
      if (rnd(1) <= (dN_inj - real(N_inj_tmp,real64))) N_inj_tmp = N_inj_tmp + 1_int32
    else
      N_inj_tmp = N_inj(2, iproc)  ! match ion wall-loss count (ptype=2)
    end if

    do i = 1_int32, N_inj_tmp

70    rnd(1) = ran2(iseed)
      rnd(2) = ran2(iseed)
      rnd(3) = ran2(iseed)

      if (cfg%opt_inj > 0_int32) then
        x = -rnd(1)*dx + cfg%xr_pow
        y = -rnd(2)*dy + cfg%yr_pow
        z = -rnd(3)*dz + cfg%zr_pow
      else
        x = xm + (dx/pi) * asin(2.0_real64*rnd(1) - 1.0_real64)
        if (flag_pbc == 0 .or. &
            (flag_pbc == 1 .and. flag_pbcz == 1)) then
          y = ym + (dy/pi) * asin(2.0_real64*rnd(2) - 1.0_real64)
        else
          y = -rnd(2)*dy + cfg%yr_pow
        end if
        if (flag_pbc == 0) then
          z = zm + (real(n(3),real64)*h(3)/pi) * asin(2.0_real64*rnd(3) - 1.0_real64)
        else
          z = -rnd(3)*dz + cfg%zr_pow
        end if
      end if

      if (abs(cfg%opt_inj) == 3_int32 .or. abs(cfg%opt_inj) == 4_int32) then
        if (((x-xm)**2 + (y-ym)**2) > Rb**2) goto 70
      end if

      ix = floor(x / h(1)) + 1_int32
      iy = floor(y / h(2)) + 1_int32
      iz = floor(z / h(3)) + 1_int32

      ! Guard against a candidate position landing outside the bcnd array
      ! (e.g. a heating-region bound in conditions.inp - xl_pow/xr_pow/
      ! yl_pow/yr_pow/zl_pow/zr_pow - that extends past the actual domain
      ! box). Legacy's part_injection.f90 has no equivalent check and reads
      ! out-of-bounds silently (usually landing in adjacent heap memory
      ! without crashing); resample instead of risking a segfault here.
      if (ix < 0_int32 .or. ix > n(1)+1_int32 .or. &
          iy < 0_int32 .or. iy > n(2)+1_int32 .or. &
          iz < 0_int32 .or. iz > n(3)+1_int32) goto 70

      ! Only inject in open domain (not inside solid boundaries)
      if (bcnd(ix,iy,iz) >= 1 .and. bcnd(ix+1,iy,iz) >= 1 .and. &
          bcnd(ix+1,iy+1,iz) >= 1 .and. bcnd(ix,iy+1,iz) >= 1 .and. &
          bcnd(ix,iy,iz+1) >= 1 .and. bcnd(ix+1,iy,iz+1) >= 1 .and. &
          bcnd(ix+1,iy+1,iz+1) >= 1 .and. bcnd(ix,iy+1,iz+1) >= 1) goto 70

      do ptype = 1_int32, ntype

        ! e-beam or cathode: inject only electron species
        if (abs(cfg%opt_inj) == 3_int32 .or. abs(cfg%opt_inj) == 4_int32) then
          if (tag_beam == 0_int32) then
            if (ptype >= 2_int32) cycle
          else
            if (ptype /= tag_beam) cycle
          end if

          if (n_cath == 2_int32 .and. abs(cfg%opt_inj) == 4_int32) then
            dir_sec_inout = dir_sec_inout + 1_int32
            if (dir_sec_inout > 2_int32) dir_sec_inout = 1_int32
            ng_sec = dir_sec_inout
          end if

          z = zg_sec(ng_sec)
        end if

        ! Probabilistic species selection via ni0 ratio
        rnd(1) = ran2(iseed)
        if (rnd(1) > ni0(ptype)) cycle

        part(ptype, iproc)%n = part(ptype, iproc)%n + 1_int32
        k = part(ptype, iproc)%n
        call part(ptype, iproc)%ensure_capacity(k)

        vt = vt0(ptype)

        if (abs(cfg%opt_inj) == 3_int32) then
          ! Shifted Maxwellian e-beam flux
          z_tmp = z
80        call shifted_maxwellian_flux_local(vz, vb, vt, fmax, iseed)
          rnd(1) = ran2(iseed)
          rnd(2) = ran2(iseed)
          call load_gauss_2d(vx, vy, vt, rnd(1), rnd(2))
          rnd(1) = ran2(iseed)
          dt_tmp = rnd(1) * dt
          z = z_tmp + vz * dt_tmp
          if (z < 0.0_real64 .or. z > real(n(3),real64)*h(3)) goto 80

        else if (abs(cfg%opt_inj) == 4_int32) then
          ! Half-Maxwellian cathode emission along Z
          z_tmp = z
          vb = sqrt(2.0_real64 * qe * cfg%THm / abs(mass(1)))
85        rnd(1) = ran2(iseed)
          vz = real(dir_sec_inout, real64) * vb * sqrt(-log(1.0_real64 - rnd(1)))
          rnd(1) = ran2(iseed)
          rnd(2) = ran2(iseed)
          call load_gauss_2d(vx, vy, vb, rnd(1), rnd(2))
          rnd(1) = ran2(iseed)
          dt_tmp = rnd(1) * dt
          z = z_tmp + vz * dt_tmp
          if (z < 0.0_real64 .or. z > real(n(3),real64)*h(3)) goto 85

        else
          ! Gaussian thermal distribution
          rnd(1) = ran2(iseed)
          rnd(2) = ran2(iseed)
          call load_gauss_2d(vx, vy, vt, rnd(1), rnd(2))
          if (vz_sav(ptype) == 0.0_real64) then
            rnd(1) = ran2(iseed)
            rnd(2) = ran2(iseed)
            call load_gauss_2d(vx, vy, vt, rnd(1), rnd(2))
            vz = vx
            vz_sav(ptype) = vy
          else
            vz = vz_sav(ptype)
            vz_sav(ptype) = 0.0_real64
          end if
        end if

        part(ptype,iproc)%pv(1,k)  = x
        part(ptype,iproc)%pv(2,k)  = y
        part(ptype,iproc)%pv(3,k)  = z
        part(ptype,iproc)%pv(4,k) = vx
        part(ptype,iproc)%pv(5,k) = vy
        part(ptype,iproc)%pv(6,k) = vz
        part(ptype,iproc)%w(k)  = 1.0_real64
        part(ptype,iproc)%sp(k) = ptype
        part(ptype,iproc)%flag_dead(k) = 0_int8
        part(ptype,iproc)%flag_cex(k)  = 0_int32

        ix = floor(x / h(1)) + 1_int32
        iy = floor(y / h(2)) + 1_int32
        iz = floor(z / h(3)) + 1_int32

        if (iz == iz_pl) &
          sour_xy(ix,iy,ptype,iproc) = sour_xy(ix,iy,ptype,iproc) + 1_int32
        if (iy == n(2)/2+1) &
          sour_xz(ix,iz,ptype,iproc) = sour_xz(ix,iz,ptype,iproc) + 1_int32
        if (ix == ix_pl) &
          sour_yz(iy,iz,ptype,iproc) = sour_yz(iy,iz,ptype,iproc) + 1_int32

        P_loss(4,ptype,iproc) = P_loss(4,ptype,iproc) + &
          0.5_real64 * Nm(ptype) * mass(ptype) * (vx*vx + vy*vy + vz*vz)

        mom_loss(1,4,ptype,iproc) = mom_loss(1,4,ptype,iproc) + Nm(ptype)*mass(ptype)*vx
        mom_loss(2,4,ptype,iproc) = mom_loss(2,4,ptype,iproc) + Nm(ptype)*mass(ptype)*vy
        mom_loss(3,4,ptype,iproc) = mom_loss(3,4,ptype,iproc) + Nm(ptype)*mass(ptype)*vz

        ! For e-beam: also account for potential energy at injection point
        if (abs(cfg%opt_inj) == 3_int32) then
          px = (real(ix,real64)*h(1) - x) / h(1)
          py = (real(iy,real64)*h(2) - y) / h(2)
          pz = (real(iz,real64)*h(3) - z) / h(3)
          ki(1) = px*py*pz
          ki(2) = (1.0_real64-px)*py*pz
          ki(3) = (1.0_real64-px)*(1.0_real64-py)*pz
          ki(4) = px*(1.0_real64-py)*pz
          ki(5) = px*py*(1.0_real64-pz)
          ki(6) = (1.0_real64-px)*py*(1.0_real64-pz)
          ki(7) = (1.0_real64-px)*(1.0_real64-py)*(1.0_real64-pz)
          ki(8) = px*(1.0_real64-py)*(1.0_real64-pz)
          phip = ki(1)*phi(ix,iy,iz)     + ki(2)*phi(ix+1,iy,iz)   + &
                 ki(3)*phi(ix+1,iy+1,iz) + ki(4)*phi(ix,iy+1,iz)   + &
                 ki(5)*phi(ix,iy,iz+1)   + ki(6)*phi(ix+1,iy,iz+1) + &
                 ki(7)*phi(ix+1,iy+1,iz+1) + ki(8)*phi(ix,iy+1,iz+1)
          P_loss(4,ptype,iproc) = P_loss(4,ptype,iproc) + &
            Nm(ptype) * charge(ptype) * phip
        end if

      end do ! ptype
    end do ! N_inj_tmp

  end subroutine inject_particles_volume


  !=========================================================================
  ! Surface flux injection — port of part_flux_injection() + load_flux_OMP()
  ! Called OUTSIDE any OMP region; spawns its own !$omp parallel do.
  ! Only the jne>0 negative-ion path is implemented (flag_thr==0 assumed).
  !=========================================================================
  subroutine inject_flux_particles( &
      part, cfg, n, h, bcnd, &
      Nm, mass, &
      ntype, nproc, mpi_rank, mpi_size, dt, istep, nsav, &
      xg1, Lgy, Lgz, ymax, zmax, Sg, &
      tag_neg, tag_neu, &
      iseed, &
      sour_xy, sour_xz, sour_fx_yz, &
      iz_pl, ix_pl, &
      P_loss, N_flx, mom_loss )

    type(ParticleSet), intent(inout) :: part(ntype, nproc)
    type(Config),      intent(in)    :: cfg
    integer(int32),    intent(in)    :: n(3)
    real(real64),      intent(in)    :: h(3)
    integer(int32),    intent(in)    :: bcnd(0:n(1)+2, 0:n(2)+2, 0:n(3)+2)
    real(real64),      intent(in)    :: Nm(ntype)
    real(real64),      intent(in)    :: mass(ntype)
    integer(int32),    intent(in)    :: ntype, nproc, mpi_rank, mpi_size
    real(real64),      intent(in)    :: dt
    integer(int32),    intent(in)    :: istep, nsav
    real(real64),      intent(in)    :: xg1, Lgy, Lgz, ymax, zmax, Sg
    integer(int32),    intent(in)    :: tag_neg, tag_neu
    integer(int32),    intent(inout) :: iseed(nproc)
    integer(int32),    intent(inout) :: sour_xy(0:n(1)+2, 0:n(2)+2, ntype, nproc)
    integer(int32),    intent(inout) :: sour_xz(0:n(1)+2, 0:n(3)+2, ntype, nproc)
    integer(int32),    intent(inout) :: sour_fx_yz(0:n(2)+2, 0:n(3)+2, ntype, nproc)
    integer(int32),    intent(in)    :: iz_pl, ix_pl
    real(real64),      intent(inout) :: P_loss(4, ntype, nproc)
    integer(int32),    intent(inout) :: N_flx(ntype, nproc)
    real(real64),      intent(inout) :: mom_loss(3, 4, ntype, nproc)

    integer(int32) :: iproc, ptype, Nh, Nh_tmp_shared
    real(real64)   :: jH, dNh, rnd, vt_neg

    if (cfg%jne <= 0.0_real64) return

    ptype = tag_neg
    if (ptype <= 0_int32) return

    jH    = cfg%jne * 10.0_real64  ! A/m^2 (jne input in mA/cm^2)
    dNh   = jH * Sg * dt / (Nm(ptype)*qe) / real(mpi_size,real64) / real(nproc,real64)
    Nh_tmp_shared = int(dNh, int32)
    rnd = ran2(iseed(1))
    if (rnd <= (dNh - real(Nh_tmp_shared,real64))) Nh_tmp_shared = Nh_tmp_shared + 1_int32
    N_flx(ptype,:) = Nh_tmp_shared

    ! Thermal velocity of H- from THm temperature
    vt_neg = sqrt(2.0_real64 * qe * cfg%THm / abs(mass(ptype)))

    !$omp parallel do private(iproc, Nh) schedule(static)
    do iproc = 1_int32, nproc
      Nh = N_flx(ptype, iproc)
      if (Nh <= 0_int32) cycle

      call load_flux_omp_local( &
        part       = part(ptype, iproc), &
        n          = n, &
        h          = h, &
        bcnd       = bcnd, &
        Nh         = Nh, &
        ptype      = ptype, &
        iproc      = iproc, &
        nproc      = nproc, &
        iseed      = iseed(iproc), &
        Nm_p       = Nm(ptype), &
        mass_p     = mass(ptype), &
        vt         = vt_neg, &
        dt         = dt, &
        xg1        = xg1, &
        Lgy        = Lgy, &
        Lgz        = Lgz, &
        ymax       = ymax, &
        zmax       = zmax, &
        ix_inj     = int(xg1/h(1),int32) + 1_int32, &
        iz_pl      = iz_pl, &
        tag_neu    = tag_neu, &
        ntype_s    = ntype, &
        nproc_s    = nproc, &
        sour_xy    = sour_xy, &
        sour_xz    = sour_xz, &
        sour_fx_yz = sour_fx_yz, &
        P_loss_inj = P_loss(4, ptype, iproc), &
        mom_loss_inj = mom_loss(:, 4, ptype, iproc) )
    end do
    !$omp end parallel do

    ! Source term diagnostic output — matches legacy's MOD(istep,nsav)==1 write
    if (mpi_rank == 0 .and. mod(istep, nsav) == 1_int32) then
      call write_sour_fx_yz(n(2), n(3), ntype, nproc, sour_fx_yz, ptype, &
                             istep, dt, h, Nm(ptype))
    end if

    if (istep == 1_int32 .and. mpi_rank == 0) then
      write(*,'(a,f8.3,a,f6.1)') &
        ' Negative ion current off extraction electrode: IH-(A)= ', &
        dNh * real(nproc*mpi_size,real64) * qe * Nm(ptype) / dt, &
        ', jH-(A/m2)= ', jH
    end if

  end subroutine inject_flux_particles


  !=========================================================================
  ! Internal: load_flux_OMP equivalent for one iproc's H- injection
  !=========================================================================
  subroutine load_flux_omp_local( &
      part, n, h, bcnd, Nh, ptype, iproc, nproc, iseed, &
      Nm_p, mass_p, vt, dt, &
      xg1, Lgy, Lgz, ymax, zmax, ix_inj, iz_pl, &
      tag_neu, ntype_s, nproc_s, &
      sour_xy, sour_xz, sour_fx_yz, P_loss_inj, mom_loss_inj )

    type(ParticleSet), intent(inout) :: part
    integer(int32),    intent(in)    :: n(3)
    real(real64),      intent(in)    :: h(3)
    integer(int32),    intent(in)    :: bcnd(0:n(1)+2, 0:n(2)+2, 0:n(3)+2)
    integer(int32),    intent(in)    :: Nh, ptype, iproc, nproc
    integer(int32),    intent(inout) :: iseed
    real(real64),      intent(in)    :: Nm_p, mass_p, vt, dt
    real(real64),      intent(in)    :: xg1, Lgy, Lgz, ymax, zmax
    integer(int32),    intent(in)    :: ix_inj, iz_pl, tag_neu, ntype_s, nproc_s
    integer(int32),    intent(inout) :: sour_xy(0:n(1)+2, 0:n(2)+2, ntype_s, nproc_s)
    integer(int32),    intent(inout) :: sour_xz(0:n(1)+2, 0:n(3)+2, ntype_s, nproc_s)
    integer(int32),    intent(inout) :: sour_fx_yz(0:n(2)+2, 0:n(3)+2, ntype_s, nproc_s)
    real(real64),      intent(inout) :: P_loss_inj
    real(real64),      intent(inout) :: mom_loss_inj(3)

    integer(int32) :: i, k, iy, iz
    real(real64)   :: rnd(2), y_inj, z_inj, vx, vy, vz, x_inj, dt_tmp

    x_inj = xg1

    do i = 1_int32, Nh

90    rnd(1) = ran2(iseed)
      y_inj = (ymax - Lgy) / 2.0_real64 + Lgy * rnd(1)
      iy    = int(y_inj / h(2), int32) + 1_int32

      rnd(1) = ran2(iseed)
      z_inj = (zmax - Lgz) / 2.0_real64 + Lgz * rnd(1)
      iz    = int(z_inj / h(3), int32) + 1_int32

      ! Must land on wall surface (all 8 surrounding cells solid/wall)
      if (.not. (bcnd(ix_inj,iy,iz) >= 1 .and. &
                 bcnd(ix_inj,iy+1,iz) >= 1 .and. &
                 bcnd(ix_inj,iy,iz+1) >= 1 .and. &
                 bcnd(ix_inj,iy+1,iz+1) >= 1 .and. &
                 bcnd(ix_inj+1,iy,iz) >= 1 .and. &
                 bcnd(ix_inj+1,iy+1,iz) >= 1 .and. &
                 bcnd(ix_inj+1,iy,iz+1) >= 1 .and. &
                 bcnd(ix_inj+1,iy+1,iz+1) >= 1)) goto 90

      part%n = part%n + 1_int32
      k = part%n
      call part%ensure_capacity(k)

      ! H- injected toward -x with half-Maxwellian flux
      rnd(1) = ran2(iseed)
      vx = -vt * sqrt(-log(1.0_real64 - rnd(1)))   ! toward LHS (negative x)

      rnd(1) = ran2(iseed)
      rnd(2) = ran2(iseed)
      call load_gauss_2d(vy, vz, vt, rnd(1), rnd(2))

      ! Spread the flux over ns_flx*dt (ns_flx=1 so just dt)
      rnd(1) = ran2(iseed)
      dt_tmp = rnd(1) * dt

      part%pv(1,k)  = x_inj + vx * dt_tmp
      part%pv(2,k)  = y_inj
      part%pv(3,k)  = z_inj
      part%pv(4,k) = vx
      part%pv(5,k) = vy
      part%pv(6,k) = vz
      part%w(k)  = 1.0_real64
      part%sp(k) = ptype
      part%flag_dead(k) = 0_int8
      part%flag_cex(k)  = tag_neu   ! H- originated from neutral H

      sour_fx_yz(iy, iz, ptype, iproc) = sour_fx_yz(iy, iz, ptype, iproc) + 1_int32
      if (iz == iz_pl) &
        sour_xy(ix_inj, iy, ptype, iproc) = sour_xy(ix_inj, iy, ptype, iproc) + 1_int32
      if (iy == n(2)/2+1) &
        sour_xz(ix_inj, iz, ptype, iproc) = sour_xz(ix_inj, iz, ptype, iproc) + 1_int32

      P_loss_inj = P_loss_inj + 0.5_real64 * Nm_p * mass_p * (vx*vx + vy*vy + vz*vz)

      mom_loss_inj(1) = mom_loss_inj(1) + Nm_p * mass_p * vx
      mom_loss_inj(2) = mom_loss_inj(2) + Nm_p * mass_p * vy
      mom_loss_inj(3) = mom_loss_inj(3) + Nm_p * mass_p * vz

    end do

  end subroutine load_flux_omp_local


  !=========================================================================
  ! Write sour_fx_yz diagnostic — port of the MOD(istep,nsav)==1 block
  ! in part_flux_injection()
  !=========================================================================
  subroutine write_sour_fx_yz(ny, nz, ntype, nproc, sour_fx_yz, ptype, &
                               istep, dt, h, Nm_p)
    integer(int32), intent(in) :: ny, nz, ntype, nproc, ptype, istep
    integer(int32), intent(in) :: sour_fx_yz(0:ny+2, 0:nz+2, ntype, nproc)
    real(real64),   intent(in) :: dt, h(3), Nm_p

    integer(int32) :: iy, iz
    integer(int32) :: sour_sum(0:ny+2, 0:nz+2)
    real(real64)   :: k_norm
    character(1)   :: pnum

    sour_sum = sum(sour_fx_yz(:,:,ptype,:), dim=3)

    k_norm = Nm_p / (h(1)*h(2)*h(3))
    write(pnum,'(i1)') ptype

    open(14, file='./Output/Output_2D/sour_fx'//pnum//'_yz.mco', status='REPLACE')
    write(14,*) ny, nz
    do iz = nz+1, 1, -1
      write(14,'(800(e18.6,1x))') &
        (real(sour_sum(iy,iz),real64)*k_norm / (real(istep,real64)*dt), iy=1,ny+1)
    end do
    close(14)

  end subroutine write_sour_fx_yz


  !=========================================================================
  ! Helpers
  !=========================================================================
  subroutine shifted_maxwellian_flux_local(v, vb, vt, fmax, iseed)
    real(real64),   intent(out)   :: v
    real(real64),   intent(in)    :: vb, vt, fmax
    integer(int32), intent(inout) :: iseed
    real(real64) :: vm, vp, f, rnd(2)

    vm = max(0.0_real64, vb - 4.0_real64*vt)
    vp = vb + 4.0_real64*vt

90  rnd(1) = ran2(iseed)
    rnd(2) = ran2(iseed)
    v = vm + rnd(1)*(vp - vm)
    f = v * exp(-(v - vb)**2 / vt**2)
    if (rnd(2) > f/fmax) goto 90
  end subroutine shifted_maxwellian_flux_local


  subroutine load_gauss_2d(vx, vy, vt, rnd1, rnd2)
    real(real64), intent(out) :: vx, vy
    real(real64), intent(in)  :: vt, rnd1, rnd2
    real(real64) :: vp, theta
    vp    = vt * sqrt(-log(1.0_real64 - rnd1))
    theta = 2.0_real64 * pi * rnd2
    vx    = vp * cos(theta)
    vy    = vp * sin(theta)
  end subroutine load_gauss_2d

end module mod_injection