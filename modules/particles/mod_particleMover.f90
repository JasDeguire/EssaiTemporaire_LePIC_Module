module mod_particleMover
  use iso_fortran_env, only: int32, real64, int8
  use mod_particles,   only: ParticleSet
  use mod_constants,   only: qe, pi
  use mod_rng,         only: ran2, load_gauss
  use mod_particleBC,  only: particle_is_lost, SeeParams
  implicit none
  private

  public :: move_and_bc_electrostatic
  public :: move_and_bc_electrostatic_fast
  public :: move_and_bc_boris
  public :: move_and_bc_boris_fast
  public :: interpolate_E_trilinear
  public :: gather_E_energy_conserving

contains

  pure integer(int32) function clamp_index(i, ilo, ihi) result(ic)
    integer(int32), intent(in) :: i, ilo, ihi
    ic = max(ilo, min(ihi, i))
  end function clamp_index


  pure subroutine interpolate_E_trilinear(xp, yp, zp, n, h, E, Exp, Eyp, Ezp)
    real(real64),   intent(in)  :: xp, yp, zp
    integer(int32), intent(in)  :: n(3)
    real(real64),   intent(in)  :: h(3)
    real(real64),   intent(in)  :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(out) :: Exp, Eyp, Ezp

    integer(int32) :: ix, iy, iz
    real(real64)   :: px, py, pz
    real(real64)   :: wx2, wy2, wz2
    real(real64)   :: w1, w2, w3, w4, w5, w6, w7, w8

    ix = int(xp / h(1), int32) + 1_int32
    iy = int(yp / h(2), int32) + 1_int32
    iz = int(zp / h(3), int32) + 1_int32

    ix = clamp_index(ix, 0_int32, n(1)+1_int32)
    iy = clamp_index(iy, 0_int32, n(2)+1_int32)
    iz = clamp_index(iz, 0_int32, n(3)+1_int32)

    px = (real(ix, real64)*h(1) - xp) / h(1)
    py = (real(iy, real64)*h(2) - yp) / h(2)
    pz = (real(iz, real64)*h(3) - zp) / h(3)

    wx2 = 1.0_real64 - px
    wy2 = 1.0_real64 - py
    wz2 = 1.0_real64 - pz

    w1 = px  * py  * pz
    w2 = wx2 * py  * pz
    w3 = wx2 * wy2 * pz
    w4 = px  * wy2 * pz
    w5 = px  * py  * wz2
    w6 = wx2 * py  * wz2
    w7 = wx2 * wy2 * wz2
    w8 = px  * wy2 * wz2

    Exp = w1*E(1,ix  ,iy  ,iz  ) + w2*E(1,ix+1,iy  ,iz  ) + &
          w3*E(1,ix+1,iy+1,iz  ) + w4*E(1,ix  ,iy+1,iz  ) + &
          w5*E(1,ix  ,iy  ,iz+1) + w6*E(1,ix+1,iy  ,iz+1) + &
          w7*E(1,ix+1,iy+1,iz+1) + w8*E(1,ix  ,iy+1,iz+1)

    Eyp = w1*E(2,ix  ,iy  ,iz  ) + w2*E(2,ix+1,iy  ,iz  ) + &
          w3*E(2,ix+1,iy+1,iz  ) + w4*E(2,ix  ,iy+1,iz  ) + &
          w5*E(2,ix  ,iy  ,iz+1) + w6*E(2,ix+1,iy  ,iz+1) + &
          w7*E(2,ix+1,iy+1,iz+1) + w8*E(2,ix  ,iy+1,iz+1)

    Ezp = w1*E(3,ix  ,iy  ,iz  ) + w2*E(3,ix+1,iy  ,iz  ) + &
          w3*E(3,ix+1,iy+1,iz  ) + w4*E(3,ix  ,iy+1,iz  ) + &
          w5*E(3,ix  ,iy  ,iz+1) + w6*E(3,ix+1,iy  ,iz+1) + &
          w7*E(3,ix+1,iy+1,iz+1) + w8*E(3,ix  ,iy+1,iz+1)
  end subroutine interpolate_E_trilinear


  pure subroutine gather_E_energy_conserving(xp, yp, zp, n, h, E, Exp, Eyp, Ezp)
    ! Energy-conserving gather (Powis & Kaganovich, Phys. Plasmas 31, 023901
    ! (2024), Eq. 10): each E component uses the nearest-grid-point (S0,
    ! piecewise-constant) weight in its own normal direction, and ordinary
    ! linear weight in the other two (tangential) directions. This is the
    ! exact coordinate-by-coordinate derivative of the same trilinear
    ! potential-energy shape function used for charge deposition - the
    ! property that makes the scheme energy-conserving.
    !
    ! "Nearest face" in a direction turns out to be exactly the cell the
    ! particle is already in: cell ix spans [(ix-1)*h(1), ix*h(1)), and its
    ! own center (ix-0.5)*h(1) - which is where calc_Efield_energy_conserving
    ! (mod_electricField.f90) stores E(1,ix,...) - is strictly nearer to
    ! every point in that span than either neighboring cell's center. So the
    ! normal direction just reuses the same cell index ix/iy/iz as the
    ! tangential (bilinear) directions do, with no separate rounding.
    real(real64),   intent(in)  :: xp, yp, zp
    integer(int32), intent(in)  :: n(3)
    real(real64),   intent(in)  :: h(3)
    real(real64),   intent(in)  :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(out) :: Exp, Eyp, Ezp

    integer(int32) :: ix, iy, iz
    real(real64)   :: px, py, pz
    real(real64)   :: wx2, wy2, wz2

    ix = int(xp / h(1), int32) + 1_int32
    iy = int(yp / h(2), int32) + 1_int32
    iz = int(zp / h(3), int32) + 1_int32

    ix = clamp_index(ix, 0_int32, n(1)+1_int32)
    iy = clamp_index(iy, 0_int32, n(2)+1_int32)
    iz = clamp_index(iz, 0_int32, n(3)+1_int32)

    px = (real(ix, real64)*h(1) - xp) / h(1)
    py = (real(iy, real64)*h(2) - yp) / h(2)
    pz = (real(iz, real64)*h(3) - zp) / h(3)

    wx2 = 1.0_real64 - px
    wy2 = 1.0_real64 - py
    wz2 = 1.0_real64 - pz

    ! E_x: nearest (fixed ix), bilinear in y,z
    Exp = py*pz*E(1,ix,iy,iz) + py*wz2*E(1,ix,iy,iz+1) + &
          wy2*pz*E(1,ix,iy+1,iz) + wy2*wz2*E(1,ix,iy+1,iz+1)

    ! E_y: nearest (fixed iy), bilinear in x,z
    Eyp = px*pz*E(2,ix,iy,iz) + px*wz2*E(2,ix,iy,iz+1) + &
          wx2*pz*E(2,ix+1,iy,iz) + wx2*wz2*E(2,ix+1,iy,iz+1)

    ! E_z: nearest (fixed iz), bilinear in x,y
    Ezp = px*py*E(3,ix,iy,iz) + px*wy2*E(3,ix,iy+1,iz) + &
          wx2*py*E(3,ix+1,iy,iz) + wx2*wy2*E(3,ix+1,iy+1,iz)
  end subroutine gather_E_energy_conserving


  pure subroutine compute_rf_field(xp, yp, zp, time, ymax, zmax, R_ahp, gams, &
                                    E0_RF, omega_RF, flag_planar_ant, x0, &
                                    EAy, EAz, EAth, theta)
    ! RF antenna (inductive) heating field, azimuthal about the domain's
    ! y-z center. Two selectable geometries (cfg%flag_planar_ant, set from
    ! the 'Wrapped' vs 'Planar' RF mode in mod_readConditions.f90):
    !
    !   flag_planar_ant==0 (wall-coil, default): coil wrapped around the
    !   chamber wall at radius R_ahp. Field peaks at R_ahp and decays
    !   exponentially INWARD (skin depth gams) as rp decreases; zero for
    !   rp > R_ahp. Ported from legacy Src/part_expmover.f90:147-155.
    !
    !   flag_planar_ant==1 (planar-coil): coil under a dielectric window
    !   at x=x0 (=xl_pow). Radial term is the bounded Faraday's-law
    !   profile for a uniform time-varying B_z confined to r<=R_ahp, the
    !   window radius (cfg%R_win if given, else yr_pow-ymax/2)
    !   (E ~ r inside, E ~ R_ahp/r outside - continuous at rp=R_ahp,
    !   unlike a naive r/R_ahp that grows without bound past the coil's
    !   own extent). Skin depth now decays AXIALLY away from the window,
    !   not radially - the induced field penetrates into the plasma bulk
    !   along x, not toward the axis.
    !
    ! Called from both movers below, once for the field-add (before the
    ! velocity update) and reusing its EAth/theta outputs for the
    ! power-accumulation call after the position update.
    real(real64),   intent(in)  :: xp, yp, zp, time, ymax, zmax, R_ahp, gams
    real(real64),   intent(in)  :: E0_RF, omega_RF
    integer(int32), intent(in)  :: flag_planar_ant
    real(real64),   intent(in)  :: x0
    real(real64),   intent(out) :: EAy, EAz, EAth, theta
    real(real64) :: rp, radial_factor

    rp    = sqrt((yp-ymax/2.0_real64)**2 + (zp-zmax/2.0_real64)**2)
    theta = atan2(zp-zmax/2.0_real64, yp-ymax/2.0_real64)
    if (theta < 0.0_real64) theta = theta + 2.0_real64*pi

    if (flag_planar_ant == 1_int32) then
       if (rp <= R_ahp) then
          radial_factor = rp / R_ahp
       else
          radial_factor = R_ahp / rp
       end if
       EAth = E0_RF * radial_factor * exp(-(xp-x0)/gams) * cos(omega_RF*time)
    else if (rp <= R_ahp) then
       EAth = E0_RF * (rp/R_ahp) * exp((rp-R_ahp)/gams) * cos(omega_RF*time)
    else
       ! Outside the coil radius: no field. The skin-depth profile grows as
       ! exp((rp-R_ahp)/gams) past the coil, so without this cutoff any
       ! plasma at rp > R_ahp inside the xl_pow..xr_pow slab (e.g. the
       ! expansion-chamber side of the ixr_pow cell) sees runaway fields.
       ! Legacy Src/part_expmover.f90 has no cutoff.
       EAth = 0.0_real64
    end if

    EAy  = -EAth*sin(theta)
    EAz  =  EAth*cos(theta)
  end subroutine compute_rf_field


  subroutine move_and_bc_electrostatic( part, n, h, E, q, m, dt, &
                                         use_energy_conserving, &
                                         bcnd, wall_cell, xmax, ymax, zmax, flag_pbc, flag_nmn, &
                                         ptype, tag_neg, flag_die, dtype, qmacro, &
                                         sum_q_xz_local, sum_q_yz_local, p_mac_boundary, &
                                         P_loss_wall, Nm_species, see, part_electrons, &
                                         iseed, P_loss_see, &
                                         flag_RFant, ixl_pow, ixr_pow, R_ahp, gams, &
                                         E0_RF, omega_RF, time, P_RF_local, &
                                         flag_planar_ant, x0, &
                                         mom_loss_wall, mom_loss_see, mom_RF_local, &
                                         P_loss_coll, mom_loss_coll )
    ! Fused electrostatic push + boundary-condition/SEE pass: one read and
    ! one write of each live particle's position/velocity instead of two
    ! (push writing back into part%x/vx/..., then a separate BC pass
    ! re-reading that same data fresh from memory). Ported line-for-line
    ! from move_particles_electrostatic + apply_particle_bc (mod_particleBC.f90);
    ! any change here should be mirrored there and vice versa if the two
    ! non-fused routines are still needed elsewhere.
    !
    ! use_energy_conserving selects gather_E_energy_conserving (above)
    ! instead of the inlined trilinear gather below, matching whichever of
    ! calc_Efield_modular/calc_Efield_energy_conserving (mod_electricField.f90)
    ! populated E this step - see cfg%push_scheme (mod_config.f90).
    type(ParticleSet),  intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: q, m, dt
    logical,            intent(in)    :: use_energy_conserving

    integer(int32),     intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    logical(1),         intent(in)    :: wall_cell(0:n(1)+1,0:n(2)+1,0:n(3)+1)
    real(real64),       intent(in)    :: xmax, ymax, zmax
    integer(int32),     intent(in)    :: flag_pbc, flag_nmn, ptype, tag_neg
    integer(int32),     intent(in)    :: flag_die
    integer(int32),     intent(in)    :: dtype(:)
    real(real64),       intent(in)    :: qmacro
    real(real64),       intent(inout) :: sum_q_xz_local(0:n(1)+2,0:n(3)+2)
    real(real64),       intent(inout) :: sum_q_yz_local(2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(inout) :: p_mac_boundary(:,:)
    real(real64),       intent(inout) :: P_loss_wall
    real(real64),       intent(in)    :: Nm_species
    type(SeeParams),    intent(in)    :: see
    type(ParticleSet),  intent(inout), optional :: part_electrons
    integer(int32),     intent(inout) :: iseed
    real(real64),       intent(inout) :: P_loss_see

    ! RF antenna (inductive) heating - see compute_rf_field above.
    ! ixl_pow/ixr_pow: heating-slab grid-x bounds (inclusive), gate on the
    ! pre-push ix, matching legacy Src/part_expmover.f90:142.
    integer(int32),     intent(in)    :: flag_RFant, ixl_pow, ixr_pow
    real(real64),       intent(in)    :: R_ahp, gams, E0_RF, omega_RF, time
    real(real64),       intent(inout) :: P_RF_local
    integer(int32),     intent(in)    :: flag_planar_ant
    real(real64),       intent(in)    :: x0
    ! Momentum-conservation diagnostic - see mod_particleBC.f90's
    ! apply_particle_bc header comment for what mom_loss_wall/mom_loss_see/
    ! P_loss_coll+mom_loss_coll are for. mom_RF_local is the RF-antenna
    ! impulse vector counterpart of P_RF_local (q*E_RF*dt, exact - no v
    ! projection needed the way the power term needs v_theta).
    real(real64),       intent(inout) :: mom_loss_wall(3)
    real(real64),       intent(inout) :: mom_loss_see(3)
    real(real64),       intent(inout) :: mom_RF_local(3)
    real(real64),       intent(inout) :: P_loss_coll
    real(real64),       intent(inout) :: mom_loss_coll(3)

    integer(int32) :: i, i_shift, i_see, ip_sec, n_sec
    integer(int32) :: ix, iy, iz
    integer(int32) :: flag_lost, np_lost
    integer(int32) :: igrid, d_ind

    real(real64) :: qmdt
    real(real64) :: xp, yp, zp
    real(real64) :: px, py, pz
    real(real64) :: wx2, wy2, wz2
    real(real64) :: w1, w2, w3, w4, w5, w6, w7, w8
    real(real64) :: Exp, Eyp, Ezp

    real(real64) :: xp_new, yp_new, zp_new
    real(real64) :: vpx_new, vpy_new, vpz_new
    real(real64) :: ki4(4)
    real(real64) :: Ek_eV, Ek_J
    real(real64) :: rnd(2)
    real(real64) :: vx_sec, vy_sec, vz_sec

    logical :: do_see

    logical      :: in_rf_region
    real(real64) :: EAy, EAz, EAth, theta_rf, v_theta

    logical :: has_w, has_sp, has_flag_dead, has_flag_cex

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    qmdt   = dt*q/m
    ! part_electrons is OPTIONAL (not just always-passed-but-unused): when
    ! this call pushes the electron species itself, the caller
    ! (mod_state.f90) omits it rather than passing self%part(1,iproc) - the
    ! same actual storage as `part` in that case. Passing the same
    ! ParticleSet to two INTENT(INOUT) dummies is non-conforming argument
    ! association (see mod_restart.f90's E0_RF comment for the same class
    ! of issue elsewhere), and even where q<0 makes do_see provably false
    ! at runtime, the compiler can't see that - it has no choice but to
    ! assume part and part_electrons might overlap and re-load part's
    ! several allocatable-array descriptors every iteration instead of
    ! hoisting them once outside this loop. Making the argument itself
    ! absent for that call removes the aliasing possibility at the
    ! language level, not just the runtime branch.
    do_see = present(part_electrons) .and. (see%gam_sec > 0.0_real64) .and. (q > 0.0_real64)

    ! Same class of forced re-check as do_see above: part%w/sp/flag_dead/
    ! flag_cex's allocation status is a ParticleSet-level property, fixed
    ! for this whole call, but calls inside the loop below (compute_rf_field,
    ! part_electrons%ensure_capacity, etc.) keep the compiler from proving
    ! that and it re-reads each array descriptor's alloc flag on every
    ! surviving particle instead. perf confirmed these four checks as
    ! measurably hot (~5% combined, ITER case) - hoisting them here changes
    ! nothing about which branch executes for any given particle, only
    ! removes the redundant re-check.
    has_w         = allocated(part%w)
    has_sp        = allocated(part%sp)
    has_flag_dead = allocated(part%flag_dead)
    has_flag_cex  = allocated(part%flag_cex)

    np_lost = 0_int32

    do i = 1, part%n

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) /= 0_int8) then
          np_lost = np_lost + 1_int32
          ! Remove the energy/momentum this particle carried when a
          ! collision flagged it dead last step (mirrors legacy
          ! Src/part_expmover.f90:70-84) - otherwise it stays counted in
          ! Pcoll/mom_loss(:,3,:) after this step quietly discards it.
          P_loss_coll = P_loss_coll - Nm_species * 0.5_real64 * m * &
              (part%pv(4,i)*part%pv(4,i) + part%pv(5,i)*part%pv(5,i) + part%pv(6,i)*part%pv(6,i))
          mom_loss_coll(1) = mom_loss_coll(1) - Nm_species * m * part%pv(4,i)
          mom_loss_coll(2) = mom_loss_coll(2) - Nm_species * m * part%pv(5,i)
          mom_loss_coll(3) = mom_loss_coll(3) - Nm_species * m * part%pv(6,i)
          cycle
        end if
      end if

      ! ---- push ----
      xp = part%pv(1,i)
      yp = part%pv(2,i)
      zp = part%pv(3,i)

      ix = int(xp / h(1), int32) + 1_int32
      iy = int(yp / h(2), int32) + 1_int32
      iz = int(zp / h(3), int32) + 1_int32

      ix = max(0_int32, min(n(1)+1_int32, ix))
      iy = max(0_int32, min(n(2)+1_int32, iy))
      iz = max(0_int32, min(n(3)+1_int32, iz))

      px = (real(ix, real64)*h(1) - xp) / h(1)
      py = (real(iy, real64)*h(2) - yp) / h(2)
      pz = (real(iz, real64)*h(3) - zp) / h(3)

      wx2 = 1.0_real64 - px
      wy2 = 1.0_real64 - py
      wz2 = 1.0_real64 - pz

      w1 = px  * py  * pz
      w2 = wx2 * py  * pz
      w3 = wx2 * wy2 * pz
      w4 = px  * wy2 * pz
      w5 = px  * py  * wz2
      w6 = wx2 * py  * wz2
      w7 = wx2 * wy2 * wz2
      w8 = px  * wy2 * wz2

      if (use_energy_conserving) then
        call gather_E_energy_conserving(xp, yp, zp, n, h, E, Exp, Eyp, Ezp)
      else
        Exp = w1*E(1,ix  ,iy  ,iz  ) + w2*E(1,ix+1,iy  ,iz  ) + &
              w3*E(1,ix+1,iy+1,iz  ) + w4*E(1,ix  ,iy+1,iz  ) + &
              w5*E(1,ix  ,iy  ,iz+1) + w6*E(1,ix+1,iy  ,iz+1) + &
              w7*E(1,ix+1,iy+1,iz+1) + w8*E(1,ix  ,iy+1,iz+1)

        Eyp = w1*E(2,ix  ,iy  ,iz  ) + w2*E(2,ix+1,iy  ,iz  ) + &
              w3*E(2,ix+1,iy+1,iz  ) + w4*E(2,ix  ,iy+1,iz  ) + &
              w5*E(2,ix  ,iy  ,iz+1) + w6*E(2,ix+1,iy  ,iz+1) + &
              w7*E(2,ix+1,iy+1,iz+1) + w8*E(2,ix  ,iy+1,iz+1)

        Ezp = w1*E(3,ix  ,iy  ,iz  ) + w2*E(3,ix+1,iy  ,iz  ) + &
              w3*E(3,ix+1,iy+1,iz  ) + w4*E(3,ix  ,iy+1,iz  ) + &
              w5*E(3,ix  ,iy  ,iz+1) + w6*E(3,ix+1,iy  ,iz+1) + &
              w7*E(3,ix+1,iy+1,iz+1) + w8*E(3,ix  ,iy+1,iz+1)
      end if

      ! ---- RF antenna field-add (flag_RFant==1 only) ----
      in_rf_region = .false.
      if (flag_RFant == 1_int32) then
        if (ix >= ixl_pow .and. ix <= ixr_pow) then
          in_rf_region = .true.
          call compute_rf_field(xp, yp, zp, time, ymax, zmax, R_ahp, gams, &
                                 E0_RF, omega_RF, flag_planar_ant, x0, &
                                 EAy, EAz, EAth, theta_rf)
          Eyp = Eyp + EAy
          Ezp = Ezp + EAz
        end if
      end if

      vpx_new = part%pv(4,i) + qmdt*Exp
      vpy_new = part%pv(5,i) + qmdt*Eyp
      vpz_new = part%pv(6,i) + qmdt*Ezp

      xp_new = xp + dt*vpx_new
      yp_new = yp + dt*vpy_new
      zp_new = zp + dt*vpz_new

      ! ---- RF antenna power absorbed (uses pre-update part%vy/vz(i)) ----
      if (in_rf_region) then
        v_theta = -0.5_real64*(part%pv(5,i)+vpy_new)*sin(theta_rf) + &
                   0.5_real64*(part%pv(6,i)+vpz_new)*cos(theta_rf)
        P_RF_local = P_RF_local + Nm_species*q*EAth*v_theta*dt

        ! Impulse delivered by the RF field alone (exact - no v-projection
        ! needed the way the power term above needs v_theta, since impulse
        ! is linear in E). No x-component: the prescribed RF field only has
        ! azimuthal (y,z) components (EAy,EAz above).
        mom_RF_local(2) = mom_RF_local(2) + Nm_species*q*EAy*dt
        mom_RF_local(3) = mom_RF_local(3) + Nm_species*q*EAz*dt
      end if

      ! ---- boundary conditions / SEE, operating on the just-pushed state ----
      ix = floor(xp_new / h(1)) + 1_int32
      iy = floor(yp_new / h(2)) + 1_int32
      iz = floor(zp_new / h(3)) + 1_int32

      if (ix < 0_int32)      ix = 0_int32
      if (ix > n(1)+1_int32) ix = n(1)+1_int32
      if (iy < 0_int32)      iy = 0_int32
      if (iy > n(2)+1_int32) iy = n(2)+1_int32
      if (iz < 0_int32)      iz = 0_int32
      if (iz > n(3)+1_int32) iz = n(3)+1_int32

      flag_lost = 0_int32

      if (particle_is_lost(wall_cell, ix, iy, iz, n)) flag_lost = 1_int32

      if (ptype == tag_neg) then
        if (xp_new < 0.0_real64 .and. flag_nmn == 1_int32) flag_lost = 2_int32
      end if

      if (flag_lost >= 1_int32) then

        igrid = bcnd(ix,iy,iz)

        if (flag_lost == 2_int32 .and. ptype == tag_neg) igrid = 0_int32

        if (flag_die == 1_int32 .and. igrid > 0_int32) then
          if (dtype(igrid) > 1_int32) then

            if (flag_pbc == 1_int32) then
              if (zp_new >= zmax) then
                zp_new = zp_new - zmax
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
              if (zp_new <= 0.0_real64) then
                zp_new = zmax + zp_new
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
            end if

            if (ix < 0_int32)      ix = 0_int32
            if (ix > n(1)+1_int32) ix = n(1)+1_int32
            if (iy < 0_int32)      iy = 0_int32
            if (iy > n(2)+1_int32) iy = n(2)+1_int32
            if (iz < 0_int32)      iz = 0_int32
            if (iz > n(3)+1_int32) iz = n(3)+1_int32

            pz = (real(iz,real64)*h(3) - zp_new) / h(3)

            if (dtype(igrid) == 2_int32) then
              px = (real(ix,real64)*h(1) - xp_new) / h(1)

              ki4(1) = qmacro * px               * pz
              ki4(2) = qmacro * (1.0_real64-px) * pz
              ki4(3) = qmacro * (1.0_real64-px) * (1.0_real64-pz)
              ki4(4) = qmacro * px               * (1.0_real64-pz)

              sum_q_xz_local(ix  ,iz  ) = sum_q_xz_local(ix  ,iz  ) + ki4(1)
              sum_q_xz_local(ix+1,iz  ) = sum_q_xz_local(ix+1,iz  ) + ki4(2)
              sum_q_xz_local(ix+1,iz+1) = sum_q_xz_local(ix+1,iz+1) + ki4(3)
              sum_q_xz_local(ix  ,iz+1) = sum_q_xz_local(ix  ,iz+1) + ki4(4)
            end if

            if (dtype(igrid) == 3_int32 .or. dtype(igrid) == 4_int32) then
              py = (real(iy,real64)*h(2) - yp_new) / h(2)

              ki4(1) = qmacro * py               * pz
              ki4(2) = qmacro * (1.0_real64-py) * pz
              ki4(3) = qmacro * (1.0_real64-py) * (1.0_real64-pz)
              ki4(4) = qmacro * py               * (1.0_real64-pz)

              d_ind = dtype(igrid)

              sum_q_yz_local(d_ind-2,iy  ,iz  ) = sum_q_yz_local(d_ind-2,iy  ,iz  ) + ki4(1)
              sum_q_yz_local(d_ind-2,iy+1,iz  ) = sum_q_yz_local(d_ind-2,iy+1,iz  ) + ki4(2)
              sum_q_yz_local(d_ind-2,iy+1,iz+1) = sum_q_yz_local(d_ind-2,iy+1,iz+1) + ki4(3)
              sum_q_yz_local(d_ind-2,iy  ,iz+1) = sum_q_yz_local(d_ind-2,iy  ,iz+1) + ki4(4)
            end if

          end if
        end if

        ! Wall diagnostics
        if (igrid < 0_int32) igrid = 0_int32

        Ek_eV = 0.5_real64 * m * &
          (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new) / qe

        p_mac_boundary(1,igrid) = p_mac_boundary(1,igrid) + qmacro
        p_mac_boundary(2,igrid) = p_mac_boundary(2,igrid) + abs(qmacro) * Ek_eV

        Ek_J = 0.5_real64 * m * &
            (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new)
        P_loss_wall = P_loss_wall + Nm_species * Ek_J

        mom_loss_wall(1) = mom_loss_wall(1) + Nm_species * m * vpx_new
        mom_loss_wall(2) = mom_loss_wall(2) + Nm_species * m * vpy_new
        mom_loss_wall(3) = mom_loss_wall(3) + Nm_species * m * vpz_new

        ! Secondary electron emission (positive ions hitting igrid_sec)
        if (do_see .and. igrid == see%igrid_sec) then
          rnd(1) = ran2(iseed)
          n_sec  = int(see%gam_sec, int32)
          if (rnd(1) <= (see%gam_sec - real(n_sec, real64))) n_sec = n_sec + 1_int32

          if (n_sec > 0_int32) then
            ! No lock needed: this routine's only caller (state%advance_particles_local,
            ! mod_state.f90) parallelizes over iproc alone, one thread owning
            ! this iproc's part_electrons = part(1,iproc) for the whole call.
            call part_electrons%ensure_capacity(part_electrons%n + n_sec)
            do ip_sec = 1, n_sec
              rnd(1) = ran2(iseed)
              vz_sec = -sign(1.0_real64, vpz_new) * see%vt_sec * sqrt(-log(1.0_real64 - rnd(1)))
              rnd(1) = ran2(iseed)
              rnd(2) = ran2(iseed)
              call load_gauss(vx_sec, vy_sec, see%vt_sec, rnd)

              part_electrons%n    = part_electrons%n + 1_int32
              i_see               = part_electrons%n
              part_electrons%pv(1,i_see)  = xp_new
              part_electrons%pv(2,i_see)  = yp_new
              part_electrons%pv(3,i_see)  = merge(see%zg_sec(1), see%zg_sec(2), vpz_new < 0.0_real64)
              part_electrons%pv(4,i_see) = vx_sec
              part_electrons%pv(5,i_see) = vy_sec
              part_electrons%pv(6,i_see) = vz_sec
              if (allocated(part_electrons%flag_dead)) part_electrons%flag_dead(i_see) = 0_int8
              if (allocated(part_electrons%flag_cex))  part_electrons%flag_cex(i_see)  = 0_int32

              P_loss_see = P_loss_see + 0.5_real64 * see%Nm_e * &
                  (vx_sec*vx_sec + vy_sec*vy_sec + vz_sec*vz_sec)

              mom_loss_see(1) = mom_loss_see(1) + see%Nm_e * see%mass_e * vx_sec
              mom_loss_see(2) = mom_loss_see(2) + see%Nm_e * see%mass_e * vy_sec
              mom_loss_see(3) = mom_loss_see(3) + see%Nm_e * see%mass_e * vz_sec
            end do
          end if
        end if

        np_lost = np_lost + 1_int32
        cycle
      end if

      ! Survivors: Neumann reflection on LHS only
      if (flag_nmn == 1_int32) then
        if (xp_new <= 0.0_real64) then
          xp_new  = -xp_new
          vpx_new = -vpx_new
        end if
      end if

      if (flag_pbc == 1_int32) then
        if (yp_new >= ymax) yp_new = yp_new - ymax
        if (yp_new <= 0.0_real64) yp_new = ymax + yp_new
        if (zp_new >= zmax) zp_new = zp_new - zmax
        if (zp_new <= 0.0_real64) zp_new = zmax + zp_new
      end if

      i_shift = i - np_lost

      part%pv(1,i_shift)  = xp_new
      part%pv(2,i_shift)  = yp_new
      part%pv(3,i_shift)  = zp_new
      part%pv(4,i_shift) = vpx_new
      part%pv(5,i_shift) = vpy_new
      part%pv(6,i_shift) = vpz_new

      if (has_w)         part%w(i_shift)         = part%w(i)
      if (has_sp)        part%sp(i_shift)        = part%sp(i)
      if (has_flag_dead) part%flag_dead(i_shift) = 0_int8
      if (has_flag_cex)  part%flag_cex(i_shift)  = part%flag_cex(i)

    end do

    part%n = part%n - np_lost
    if (part%n < 0_int32) part%n = 0_int32

    ! Keep the invariant "slots past n have zeroed flags" by clearing only
    ! the np_lost slots compaction just vacated. This used to zero the
    ! whole spare capacity plus cell_id/cell_count/cell_start (~19 MB of
    ! ncells-sized arrays per species per iproc) on every call - pure
    ! memory-bandwidth cost that dominated the push at high thread counts
    ! and is why it barely gained from hyperthreading. The cell lists are
    ! only read by collisions, right after sort_particles_local rebuilds
    ! them, so invalidating them here bought nothing.
    if (allocated(part%flag_dead)) then
      part%flag_dead(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int8
    end if
    if (allocated(part%flag_cex)) then
      part%flag_cex(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int32
    end if

  end subroutine move_and_bc_electrostatic


  subroutine move_and_bc_electrostatic_fast( part, n, h, E, q, m, dt, &
                                              bcnd, wall_cell, xmax, ymax, zmax, flag_pbc, flag_nmn, &
                                              ptype, tag_neg, flag_die, dtype, qmacro, &
                                              sum_q_xz_local, sum_q_yz_local, p_mac_boundary, &
                                              P_loss_wall, Nm_species, see, part_electrons, &
                                              iseed, P_loss_see, &
                                              mom_loss_wall, mom_loss_see, &
                                              P_loss_coll, mom_loss_coll )
    ! Lean fast path for move_and_bc_electrostatic, used whenever the
    ! energy-conserving pusher and RF-antenna heating are both off (the
    ! common case - both are opt-in features). Identical to the full
    ! routine with those two branches removed from the hot per-particle
    ! loop entirely (not just short-circuited): a local throwaway
    ! experiment measured this as a ~25% mover_push reduction, evidently
    ! from branch-density/ILP effects near the cache-miss-heavy boundary
    ! check, not from the branches' own cost. The BC/SEE tail below is a
    ! deliberate duplicate of move_and_bc_electrostatic's, not a shared
    ! subroutine call - an earlier attempt to extract it into a shared
    ! subroutine (relying on -ipo to inline it back) produced a 3-5x
    ! regression instead of the expected win, apparently from real
    ! per-call overhead that this file's -ipo build wasn't eliminating.
    ! Any future BC/SEE fix must be mirrored in both routines (and in
    ! move_and_bc_boris, which has its own copy of the same block).
    type(ParticleSet),  intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(in)    :: q, m, dt

    integer(int32),     intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    logical(1),         intent(in)    :: wall_cell(0:n(1)+1,0:n(2)+1,0:n(3)+1)
    real(real64),       intent(in)    :: xmax, ymax, zmax
    integer(int32),     intent(in)    :: flag_pbc, flag_nmn, ptype, tag_neg
    integer(int32),     intent(in)    :: flag_die
    integer(int32),     intent(in)    :: dtype(:)
    real(real64),       intent(in)    :: qmacro
    real(real64),       intent(inout) :: sum_q_xz_local(0:n(1)+2,0:n(3)+2)
    real(real64),       intent(inout) :: sum_q_yz_local(2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(inout) :: p_mac_boundary(:,:)
    real(real64),       intent(inout) :: P_loss_wall
    real(real64),       intent(in)    :: Nm_species
    type(SeeParams),    intent(in)    :: see
    type(ParticleSet),  intent(inout), optional :: part_electrons
    integer(int32),     intent(inout) :: iseed
    real(real64),       intent(inout) :: P_loss_see
    ! Momentum-conservation diagnostic - see move_and_bc_electrostatic's
    ! header comment.
    real(real64),       intent(inout) :: mom_loss_wall(3)
    real(real64),       intent(inout) :: mom_loss_see(3)
    real(real64),       intent(inout) :: P_loss_coll
    real(real64),       intent(inout) :: mom_loss_coll(3)

    integer(int32) :: i, i_shift, i_see, ip_sec, n_sec
    integer(int32) :: ix, iy, iz
    integer(int32) :: flag_lost, np_lost
    integer(int32) :: igrid, d_ind

    real(real64) :: qmdt
    real(real64) :: xp, yp, zp
    real(real64) :: px, py, pz
    real(real64) :: wx2, wy2, wz2
    real(real64) :: w1, w2, w3, w4, w5, w6, w7, w8
    real(real64) :: Exp, Eyp, Ezp

    real(real64) :: xp_new, yp_new, zp_new
    real(real64) :: vpx_new, vpy_new, vpz_new
    real(real64) :: ki4(4)
    real(real64) :: Ek_eV, Ek_J
    real(real64) :: rnd(2)
    real(real64) :: vx_sec, vy_sec, vz_sec

    logical :: do_see

    logical :: has_w, has_sp, has_flag_dead, has_flag_cex

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    qmdt   = dt*q/m
    ! part_electrons is OPTIONAL - see move_and_bc_electrostatic's header
    ! comment on why the electron-species call omits it entirely rather
    ! than passing self%part(1,iproc) aliased against `part`.
    do_see = present(part_electrons) .and. (see%gam_sec > 0.0_real64) .and. (q > 0.0_real64)

    ! See move_and_bc_electrostatic's header comment on has_w/has_sp/
    ! has_flag_dead/has_flag_cex.
    has_w         = allocated(part%w)
    has_sp        = allocated(part%sp)
    has_flag_dead = allocated(part%flag_dead)
    has_flag_cex  = allocated(part%flag_cex)

    np_lost = 0_int32

    do i = 1, part%n

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) /= 0_int8) then
          np_lost = np_lost + 1_int32
          ! Remove the energy/momentum this particle carried when a
          ! collision flagged it dead last step (mirrors legacy
          ! Src/part_expmover.f90:70-84) - otherwise it stays counted in
          ! Pcoll/mom_loss(:,3,:) after this step quietly discards it.
          P_loss_coll = P_loss_coll - Nm_species * 0.5_real64 * m * &
              (part%pv(4,i)*part%pv(4,i) + part%pv(5,i)*part%pv(5,i) + part%pv(6,i)*part%pv(6,i))
          mom_loss_coll(1) = mom_loss_coll(1) - Nm_species * m * part%pv(4,i)
          mom_loss_coll(2) = mom_loss_coll(2) - Nm_species * m * part%pv(5,i)
          mom_loss_coll(3) = mom_loss_coll(3) - Nm_species * m * part%pv(6,i)
          cycle
        end if
      end if

      ! ---- push ----
      xp = part%pv(1,i)
      yp = part%pv(2,i)
      zp = part%pv(3,i)

      ix = int(xp / h(1), int32) + 1_int32
      iy = int(yp / h(2), int32) + 1_int32
      iz = int(zp / h(3), int32) + 1_int32

      ix = max(0_int32, min(n(1)+1_int32, ix))
      iy = max(0_int32, min(n(2)+1_int32, iy))
      iz = max(0_int32, min(n(3)+1_int32, iz))

      px = (real(ix, real64)*h(1) - xp) / h(1)
      py = (real(iy, real64)*h(2) - yp) / h(2)
      pz = (real(iz, real64)*h(3) - zp) / h(3)

      wx2 = 1.0_real64 - px
      wy2 = 1.0_real64 - py
      wz2 = 1.0_real64 - pz

      w1 = px  * py  * pz
      w2 = wx2 * py  * pz
      w3 = wx2 * wy2 * pz
      w4 = px  * wy2 * pz
      w5 = px  * py  * wz2
      w6 = wx2 * py  * wz2
      w7 = wx2 * wy2 * wz2
      w8 = px  * wy2 * wz2

      Exp = w1*E(1,ix  ,iy  ,iz  ) + w2*E(1,ix+1,iy  ,iz  ) + &
            w3*E(1,ix+1,iy+1,iz  ) + w4*E(1,ix  ,iy+1,iz  ) + &
            w5*E(1,ix  ,iy  ,iz+1) + w6*E(1,ix+1,iy  ,iz+1) + &
            w7*E(1,ix+1,iy+1,iz+1) + w8*E(1,ix  ,iy+1,iz+1)

      Eyp = w1*E(2,ix  ,iy  ,iz  ) + w2*E(2,ix+1,iy  ,iz  ) + &
            w3*E(2,ix+1,iy+1,iz  ) + w4*E(2,ix  ,iy+1,iz  ) + &
            w5*E(2,ix  ,iy  ,iz+1) + w6*E(2,ix+1,iy  ,iz+1) + &
            w7*E(2,ix+1,iy+1,iz+1) + w8*E(2,ix  ,iy+1,iz+1)

      Ezp = w1*E(3,ix  ,iy  ,iz  ) + w2*E(3,ix+1,iy  ,iz  ) + &
            w3*E(3,ix+1,iy+1,iz  ) + w4*E(3,ix  ,iy+1,iz  ) + &
            w5*E(3,ix  ,iy  ,iz+1) + w6*E(3,ix+1,iy  ,iz+1) + &
            w7*E(3,ix+1,iy+1,iz+1) + w8*E(3,ix  ,iy+1,iz+1)

      vpx_new = part%pv(4,i) + qmdt*Exp
      vpy_new = part%pv(5,i) + qmdt*Eyp
      vpz_new = part%pv(6,i) + qmdt*Ezp

      xp_new = xp + dt*vpx_new
      yp_new = yp + dt*vpy_new
      zp_new = zp + dt*vpz_new

      ! ---- boundary conditions / SEE, operating on the just-pushed state ----
      ix = floor(xp_new / h(1)) + 1_int32
      iy = floor(yp_new / h(2)) + 1_int32
      iz = floor(zp_new / h(3)) + 1_int32

      if (ix < 0_int32)      ix = 0_int32
      if (ix > n(1)+1_int32) ix = n(1)+1_int32
      if (iy < 0_int32)      iy = 0_int32
      if (iy > n(2)+1_int32) iy = n(2)+1_int32
      if (iz < 0_int32)      iz = 0_int32
      if (iz > n(3)+1_int32) iz = n(3)+1_int32

      flag_lost = 0_int32

      if (particle_is_lost(wall_cell, ix, iy, iz, n)) flag_lost = 1_int32

      if (ptype == tag_neg) then
        if (xp_new < 0.0_real64 .and. flag_nmn == 1_int32) flag_lost = 2_int32
      end if

      if (flag_lost >= 1_int32) then

        igrid = bcnd(ix,iy,iz)

        if (flag_lost == 2_int32 .and. ptype == tag_neg) igrid = 0_int32

        if (flag_die == 1_int32 .and. igrid > 0_int32) then
          if (dtype(igrid) > 1_int32) then

            if (flag_pbc == 1_int32) then
              if (zp_new >= zmax) then
                zp_new = zp_new - zmax
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
              if (zp_new <= 0.0_real64) then
                zp_new = zmax + zp_new
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
            end if

            if (ix < 0_int32)      ix = 0_int32
            if (ix > n(1)+1_int32) ix = n(1)+1_int32
            if (iy < 0_int32)      iy = 0_int32
            if (iy > n(2)+1_int32) iy = n(2)+1_int32
            if (iz < 0_int32)      iz = 0_int32
            if (iz > n(3)+1_int32) iz = n(3)+1_int32

            pz = (real(iz,real64)*h(3) - zp_new) / h(3)

            if (dtype(igrid) == 2_int32) then
              px = (real(ix,real64)*h(1) - xp_new) / h(1)

              ki4(1) = qmacro * px               * pz
              ki4(2) = qmacro * (1.0_real64-px) * pz
              ki4(3) = qmacro * (1.0_real64-px) * (1.0_real64-pz)
              ki4(4) = qmacro * px               * (1.0_real64-pz)

              sum_q_xz_local(ix  ,iz  ) = sum_q_xz_local(ix  ,iz  ) + ki4(1)
              sum_q_xz_local(ix+1,iz  ) = sum_q_xz_local(ix+1,iz  ) + ki4(2)
              sum_q_xz_local(ix+1,iz+1) = sum_q_xz_local(ix+1,iz+1) + ki4(3)
              sum_q_xz_local(ix  ,iz+1) = sum_q_xz_local(ix  ,iz+1) + ki4(4)
            end if

            if (dtype(igrid) == 3_int32 .or. dtype(igrid) == 4_int32) then
              py = (real(iy,real64)*h(2) - yp_new) / h(2)

              ki4(1) = qmacro * py               * pz
              ki4(2) = qmacro * (1.0_real64-py) * pz
              ki4(3) = qmacro * (1.0_real64-py) * (1.0_real64-pz)
              ki4(4) = qmacro * py               * (1.0_real64-pz)

              d_ind = dtype(igrid)

              sum_q_yz_local(d_ind-2,iy  ,iz  ) = sum_q_yz_local(d_ind-2,iy  ,iz  ) + ki4(1)
              sum_q_yz_local(d_ind-2,iy+1,iz  ) = sum_q_yz_local(d_ind-2,iy+1,iz  ) + ki4(2)
              sum_q_yz_local(d_ind-2,iy+1,iz+1) = sum_q_yz_local(d_ind-2,iy+1,iz+1) + ki4(3)
              sum_q_yz_local(d_ind-2,iy  ,iz+1) = sum_q_yz_local(d_ind-2,iy  ,iz+1) + ki4(4)
            end if

          end if
        end if

        ! Wall diagnostics
        if (igrid < 0_int32) igrid = 0_int32

        Ek_eV = 0.5_real64 * m * &
          (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new) / qe

        p_mac_boundary(1,igrid) = p_mac_boundary(1,igrid) + qmacro
        p_mac_boundary(2,igrid) = p_mac_boundary(2,igrid) + abs(qmacro) * Ek_eV

        Ek_J = 0.5_real64 * m * &
            (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new)
        P_loss_wall = P_loss_wall + Nm_species * Ek_J

        mom_loss_wall(1) = mom_loss_wall(1) + Nm_species * m * vpx_new
        mom_loss_wall(2) = mom_loss_wall(2) + Nm_species * m * vpy_new
        mom_loss_wall(3) = mom_loss_wall(3) + Nm_species * m * vpz_new

        ! Secondary electron emission (positive ions hitting igrid_sec)
        if (do_see .and. igrid == see%igrid_sec) then
          rnd(1) = ran2(iseed)
          n_sec  = int(see%gam_sec, int32)
          if (rnd(1) <= (see%gam_sec - real(n_sec, real64))) n_sec = n_sec + 1_int32

          if (n_sec > 0_int32) then
            ! No lock needed: this routine's only caller (state%advance_particles_local,
            ! mod_state.f90) parallelizes over iproc alone, one thread owning
            ! this iproc's part_electrons = part(1,iproc) for the whole call.
            call part_electrons%ensure_capacity(part_electrons%n + n_sec)
            do ip_sec = 1, n_sec
              rnd(1) = ran2(iseed)
              vz_sec = -sign(1.0_real64, vpz_new) * see%vt_sec * sqrt(-log(1.0_real64 - rnd(1)))
              rnd(1) = ran2(iseed)
              rnd(2) = ran2(iseed)
              call load_gauss(vx_sec, vy_sec, see%vt_sec, rnd)

              part_electrons%n    = part_electrons%n + 1_int32
              i_see               = part_electrons%n
              part_electrons%pv(1,i_see)  = xp_new
              part_electrons%pv(2,i_see)  = yp_new
              part_electrons%pv(3,i_see)  = merge(see%zg_sec(1), see%zg_sec(2), vpz_new < 0.0_real64)
              part_electrons%pv(4,i_see) = vx_sec
              part_electrons%pv(5,i_see) = vy_sec
              part_electrons%pv(6,i_see) = vz_sec
              if (allocated(part_electrons%flag_dead)) part_electrons%flag_dead(i_see) = 0_int8
              if (allocated(part_electrons%flag_cex))  part_electrons%flag_cex(i_see)  = 0_int32

              P_loss_see = P_loss_see + 0.5_real64 * see%Nm_e * &
                  (vx_sec*vx_sec + vy_sec*vy_sec + vz_sec*vz_sec)

              mom_loss_see(1) = mom_loss_see(1) + see%Nm_e * see%mass_e * vx_sec
              mom_loss_see(2) = mom_loss_see(2) + see%Nm_e * see%mass_e * vy_sec
              mom_loss_see(3) = mom_loss_see(3) + see%Nm_e * see%mass_e * vz_sec
            end do
          end if
        end if

        np_lost = np_lost + 1_int32
        cycle
      end if

      ! Survivors: Neumann reflection on LHS only
      if (flag_nmn == 1_int32) then
        if (xp_new <= 0.0_real64) then
          xp_new  = -xp_new
          vpx_new = -vpx_new
        end if
      end if

      if (flag_pbc == 1_int32) then
        if (yp_new >= ymax) yp_new = yp_new - ymax
        if (yp_new <= 0.0_real64) yp_new = ymax + yp_new
        if (zp_new >= zmax) zp_new = zp_new - zmax
        if (zp_new <= 0.0_real64) zp_new = zmax + zp_new
      end if

      i_shift = i - np_lost

      part%pv(1,i_shift)  = xp_new
      part%pv(2,i_shift)  = yp_new
      part%pv(3,i_shift)  = zp_new
      part%pv(4,i_shift) = vpx_new
      part%pv(5,i_shift) = vpy_new
      part%pv(6,i_shift) = vpz_new

      if (has_w)         part%w(i_shift)         = part%w(i)
      if (has_sp)        part%sp(i_shift)        = part%sp(i)
      if (has_flag_dead) part%flag_dead(i_shift) = 0_int8
      if (has_flag_cex)  part%flag_cex(i_shift)  = part%flag_cex(i)

    end do

    part%n = part%n - np_lost
    if (part%n < 0_int32) part%n = 0_int32

    ! Keep the invariant "slots past n have zeroed flags" by clearing only
    ! the np_lost slots compaction just vacated. This used to zero the
    ! whole spare capacity plus cell_id/cell_count/cell_start (~19 MB of
    ! ncells-sized arrays per species per iproc) on every call - pure
    ! memory-bandwidth cost that dominated the push at high thread counts
    ! and is why it barely gained from hyperthreading. The cell lists are
    ! only read by collisions, right after sort_particles_local rebuilds
    ! them, so invalidating them here bought nothing.
    if (allocated(part%flag_dead)) then
      part%flag_dead(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int8
    end if
    if (allocated(part%flag_cex)) then
      part%flag_cex(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int32
    end if

  end subroutine move_and_bc_electrostatic_fast


  subroutine move_and_bc_boris( part, n, h, E, n_B, h_B, Bi, q, m, dt, &
                                 use_energy_conserving, &
                                 bcnd, wall_cell, xmax, ymax, zmax, flag_pbc, flag_nmn, &
                                 ptype, tag_neg, flag_die, dtype, qmacro, &
                                 sum_q_xz_local, sum_q_yz_local, p_mac_boundary, &
                                 P_loss_wall, Nm_species, see, part_electrons, &
                                 iseed, P_loss_see, &
                                 flag_RFant, ixl_pow, ixr_pow, R_ahp, gams, &
                                 E0_RF, omega_RF, time, P_RF_local, &
                                 flag_planar_ant, x0, &
                                 mom_loss_wall, mom_loss_see, mom_RF_local, &
                                 P_loss_coll, mom_loss_coll, &
                                 n_boris_used, n_boris_total )
    ! Same fusion as move_and_bc_electrostatic, but with the Boris push
    ! (move_particles_boris) in place of the electrostatic one. See that
    ! routine's header comment for the rationale; the BC/SEE block below
    ! is identical to it (and to apply_particle_bc in mod_particleBC.f90).
    ! use_energy_conserving: see move_and_bc_electrostatic - only affects
    ! how E is gathered, the magnetic (Boris rotation) part is unchanged.
    type(ParticleSet),  intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32),     intent(in)    :: n_B(3)
    real(real64),       intent(in)    :: h_B(3)
    real(real64),       intent(in)    :: Bi(4,0:n_B(1)+2,0:n_B(2)+2,0:n_B(3)+2)
    real(real64),       intent(in)    :: q, m, dt
    logical,            intent(in)    :: use_energy_conserving

    integer(int32),     intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    logical(1),         intent(in)    :: wall_cell(0:n(1)+1,0:n(2)+1,0:n(3)+1)
    real(real64),       intent(in)    :: xmax, ymax, zmax
    integer(int32),     intent(in)    :: flag_pbc, flag_nmn, ptype, tag_neg
    integer(int32),     intent(in)    :: flag_die
    integer(int32),     intent(in)    :: dtype(:)
    real(real64),       intent(in)    :: qmacro
    real(real64),       intent(inout) :: sum_q_xz_local(0:n(1)+2,0:n(3)+2)
    real(real64),       intent(inout) :: sum_q_yz_local(2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(inout) :: p_mac_boundary(:,:)
    real(real64),       intent(inout) :: P_loss_wall
    real(real64),       intent(in)    :: Nm_species
    type(SeeParams),    intent(in)    :: see
    type(ParticleSet),  intent(inout), optional :: part_electrons
    integer(int32),     intent(inout) :: iseed
    real(real64),       intent(inout) :: P_loss_see

    ! RF antenna (inductive) heating - see compute_rf_field above.
    integer(int32),     intent(in)    :: flag_RFant, ixl_pow, ixr_pow
    real(real64),       intent(in)    :: R_ahp, gams, E0_RF, omega_RF, time
    real(real64),       intent(inout) :: P_RF_local
    integer(int32),     intent(in)    :: flag_planar_ant
    real(real64),       intent(in)    :: x0
    ! Momentum-conservation diagnostic - see move_and_bc_electrostatic's
    ! header comment.
    real(real64),       intent(inout) :: mom_loss_wall(3)
    real(real64),       intent(inout) :: mom_loss_see(3)
    real(real64),       intent(inout) :: mom_RF_local(3)
    real(real64),       intent(inout) :: P_loss_coll
    real(real64),       intent(inout) :: mom_loss_coll(3)
    ! TEMPORARY diagnostic (see mod_state.f90's advance_particles_local and
    ! mod_simulation.f90's diagnostic print): counts how many live particles
    ! this call actually ran the full Boris rotation for vs. how many hit
    ! the do_boris==.false. negligible-|B| fallback below, to check whether
    ! that skip is firing at all on a given case before trusting it as a
    ! performance explanation. Remove once that question is answered.
    integer(int32),     intent(out)   :: n_boris_used, n_boris_total

    integer(int32) :: i, i_shift, i_see, ip_sec, n_sec
    integer(int32) :: ix, iy, iz
    integer(int32) :: ixB, iyB, izB
    integer(int32) :: flag_lost, np_lost
    integer(int32) :: igrid, d_ind

    real(real64) :: qm2dt
    real(real64) :: xp, yp, zp
    real(real64) :: px, py, pz
    real(real64) :: wx2, wy2, wz2
    real(real64) :: w1, w2, w3, w4, w5, w6, w7, w8

    real(real64) :: pxB, pyB, pzB
    real(real64) :: wx2B, wy2B, wz2B
    real(real64) :: b1, b2, b3, b4, b5, b6, b7, b8

    real(real64) :: Exp, Eyp, Ezp
    real(real64) :: Bpx, Bpy, Bpz

    real(real64) :: vminus_x, vminus_y, vminus_z
    real(real64) :: vprime_x, vprime_y, vprime_z
    real(real64) :: vplus_x, vplus_y, vplus_z

    real(real64) :: tx, ty, tz
    real(real64) :: sx, sy, sz
    real(real64) :: t2, inv_denom

    logical :: uniform_B, same_grid_B, do_boris

    real(real64) :: xp_new, yp_new, zp_new
    real(real64) :: vpx_new, vpy_new, vpz_new
    real(real64) :: ki4(4)
    real(real64) :: Ek_eV, Ek_J
    real(real64) :: rnd(2)
    real(real64) :: vx_sec, vy_sec, vz_sec

    logical :: do_see

    logical      :: in_rf_region
    real(real64) :: EAy, EAz, EAth, theta_rf, v_theta

    logical :: has_w, has_sp, has_flag_dead, has_flag_cex

    n_boris_used  = 0_int32
    n_boris_total = 0_int32

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    qm2dt  = 0.5_real64 * dt * q / m
    ! part_electrons is OPTIONAL - see move_and_bc_electrostatic's header
    ! comment on why the electron-species call omits it entirely rather
    ! than passing self%part(1,iproc) aliased against `part`.
    do_see = present(part_electrons) .and. (see%gam_sec > 0.0_real64) .and. (q > 0.0_real64)

    ! See move_and_bc_electrostatic's header comment on has_w/has_sp/
    ! has_flag_dead/has_flag_cex.
    has_w         = allocated(part%w)
    has_sp        = allocated(part%sp)
    has_flag_dead = allocated(part%flag_dead)
    has_flag_cex  = allocated(part%flag_cex)

    uniform_B = (n_B(1) == 1_int32 .and. n_B(2) == 1_int32 .and. n_B(3) == 1_int32)

    same_grid_B = .false.
    if (.not. uniform_B) then
      same_grid_B = all(n_B == n) .and. &
                    abs(h_B(1)-h(1)) < 1.0e-14_real64 .and. &
                    abs(h_B(2)-h(2)) < 1.0e-14_real64 .and. &
                    abs(h_B(3)-h(3)) < 1.0e-14_real64
    end if

    np_lost = 0_int32

    do i = 1_int32, part%n

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) /= 0_int8) then
          np_lost = np_lost + 1_int32
          ! Remove the energy/momentum this particle carried when a
          ! collision flagged it dead last step (mirrors legacy
          ! Src/part_expmover.f90:70-84) - otherwise it stays counted in
          ! Pcoll/mom_loss(:,3,:) after this step quietly discards it.
          P_loss_coll = P_loss_coll - Nm_species * 0.5_real64 * m * &
              (part%pv(4,i)*part%pv(4,i) + part%pv(5,i)*part%pv(5,i) + part%pv(6,i)*part%pv(6,i))
          mom_loss_coll(1) = mom_loss_coll(1) - Nm_species * m * part%pv(4,i)
          mom_loss_coll(2) = mom_loss_coll(2) - Nm_species * m * part%pv(5,i)
          mom_loss_coll(3) = mom_loss_coll(3) - Nm_species * m * part%pv(6,i)
          cycle
        end if
      end if

      ! ---- push ----
      xp = part%pv(1,i)
      yp = part%pv(2,i)
      zp = part%pv(3,i)

      ix = int(xp / h(1), int32) + 1_int32
      iy = int(yp / h(2), int32) + 1_int32
      iz = int(zp / h(3), int32) + 1_int32

      ix = max(0_int32, min(n(1)+1_int32, ix))
      iy = max(0_int32, min(n(2)+1_int32, iy))
      iz = max(0_int32, min(n(3)+1_int32, iz))

      px = (real(ix, real64)*h(1) - xp) / h(1)
      py = (real(iy, real64)*h(2) - yp) / h(2)
      pz = (real(iz, real64)*h(3) - zp) / h(3)

      wx2 = 1.0_real64 - px
      wy2 = 1.0_real64 - py
      wz2 = 1.0_real64 - pz

      w1 = px  * py  * pz
      w2 = wx2 * py  * pz
      w3 = wx2 * wy2 * pz
      w4 = px  * wy2 * pz
      w5 = px  * py  * wz2
      w6 = wx2 * py  * wz2
      w7 = wx2 * wy2 * wz2
      w8 = px  * wy2 * wz2

      if (use_energy_conserving) then
        call gather_E_energy_conserving(xp, yp, zp, n, h, E, Exp, Eyp, Ezp)
      else
        Exp = w1*E(1,ix  ,iy  ,iz  ) + w2*E(1,ix+1,iy  ,iz  ) + &
              w3*E(1,ix+1,iy+1,iz  ) + w4*E(1,ix  ,iy+1,iz  ) + &
              w5*E(1,ix  ,iy  ,iz+1) + w6*E(1,ix+1,iy  ,iz+1) + &
              w7*E(1,ix+1,iy+1,iz+1) + w8*E(1,ix  ,iy+1,iz+1)

        Eyp = w1*E(2,ix  ,iy  ,iz  ) + w2*E(2,ix+1,iy  ,iz  ) + &
              w3*E(2,ix+1,iy+1,iz  ) + w4*E(2,ix  ,iy+1,iz  ) + &
              w5*E(2,ix  ,iy  ,iz+1) + w6*E(2,ix+1,iy  ,iz+1) + &
              w7*E(2,ix+1,iy+1,iz+1) + w8*E(2,ix  ,iy+1,iz+1)

        Ezp = w1*E(3,ix  ,iy  ,iz  ) + w2*E(3,ix+1,iy  ,iz  ) + &
              w3*E(3,ix+1,iy+1,iz  ) + w4*E(3,ix  ,iy+1,iz  ) + &
              w5*E(3,ix  ,iy  ,iz+1) + w6*E(3,ix+1,iy  ,iz+1) + &
              w7*E(3,ix+1,iy+1,iz+1) + w8*E(3,ix  ,iy+1,iz+1)
      end if

      ! ---- RF antenna field-add (flag_RFant==1 only) ----
      in_rf_region = .false.
      if (flag_RFant == 1_int32) then
        if (ix >= ixl_pow .and. ix <= ixr_pow) then
          in_rf_region = .true.
          call compute_rf_field(xp, yp, zp, time, ymax, zmax, R_ahp, gams, &
                                 E0_RF, omega_RF, flag_planar_ant, x0, &
                                 EAy, EAz, EAth, theta_rf)
          Eyp = Eyp + EAy
          Ezp = Ezp + EAz
        end if
      end if

      ! ------------------------------------------------------------
      ! B-field interpolation - skipped where |B| is negligible.
      ! Bi(4,:,:,:) = |B|[T] (mod_magneticField.f90) is precomputed
      ! but was previously never read here, so this routine always
      ! paid for a B gather + full Boris rotation even in field-free
      ! regions. Legacy gates the same way (Bi(4,ix,iy,iz).gt.1.d-5,
      ! i.e. >0.1G - Src/part_expmover.f90:158,214) and falls back to
      ! the cheap electrostatic-only push below otherwise. A uniform
      ! field (n_B all ==1, no spatial map) has no negligible region
      ! to skip, so it is never gated - matching legacy, which also
      ! always treats its "constant B-field" case as active.
      ! ------------------------------------------------------------
      do_boris = .true.

      if (uniform_B) then

        Bpx = Bi(1,1,1,1)
        Bpy = Bi(2,1,1,1)
        Bpz = Bi(3,1,1,1)

      else if (same_grid_B) then

        if (Bi(4,ix,iy,iz) > 1.0e-5_real64) then
          Bpx = w1*Bi(1,ix  ,iy  ,iz  ) + w2*Bi(1,ix+1,iy  ,iz  ) + &
                w3*Bi(1,ix+1,iy+1,iz  ) + w4*Bi(1,ix  ,iy+1,iz  ) + &
                w5*Bi(1,ix  ,iy  ,iz+1) + w6*Bi(1,ix+1,iy  ,iz+1) + &
                w7*Bi(1,ix+1,iy+1,iz+1) + w8*Bi(1,ix  ,iy+1,iz+1)

          Bpy = w1*Bi(2,ix  ,iy  ,iz  ) + w2*Bi(2,ix+1,iy  ,iz  ) + &
                w3*Bi(2,ix+1,iy+1,iz  ) + w4*Bi(2,ix  ,iy+1,iz  ) + &
                w5*Bi(2,ix  ,iy  ,iz+1) + w6*Bi(2,ix+1,iy  ,iz+1) + &
                w7*Bi(2,ix+1,iy+1,iz+1) + w8*Bi(2,ix  ,iy+1,iz+1)

          Bpz = w1*Bi(3,ix  ,iy  ,iz  ) + w2*Bi(3,ix+1,iy  ,iz  ) + &
                w3*Bi(3,ix+1,iy+1,iz  ) + w4*Bi(3,ix  ,iy+1,iz  ) + &
                w5*Bi(3,ix  ,iy  ,iz+1) + w6*Bi(3,ix+1,iy  ,iz+1) + &
                w7*Bi(3,ix+1,iy+1,iz+1) + w8*Bi(3,ix  ,iy+1,iz+1)
        else
          do_boris = .false.
        end if

      else

        ixB = int(xp / h_B(1), int32) + 1_int32
        iyB = int(yp / h_B(2), int32) + 1_int32
        izB = int(zp / h_B(3), int32) + 1_int32

        ixB = max(0_int32, min(n_B(1)+1_int32, ixB))
        iyB = max(0_int32, min(n_B(2)+1_int32, iyB))
        izB = max(0_int32, min(n_B(3)+1_int32, izB))

        if (Bi(4,ixB,iyB,izB) > 1.0e-5_real64) then
          pxB = (real(ixB, real64)*h_B(1) - xp) / h_B(1)
          pyB = (real(iyB, real64)*h_B(2) - yp) / h_B(2)
          pzB = (real(izB, real64)*h_B(3) - zp) / h_B(3)

          wx2B = 1.0_real64 - pxB
          wy2B = 1.0_real64 - pyB
          wz2B = 1.0_real64 - pzB

          b1 = pxB  * pyB  * pzB
          b2 = wx2B * pyB  * pzB
          b3 = wx2B * wy2B * pzB
          b4 = pxB  * wy2B * pzB
          b5 = pxB  * pyB  * wz2B
          b6 = wx2B * pyB  * wz2B
          b7 = wx2B * wy2B * wz2B
          b8 = pxB  * wy2B * wz2B

          Bpx = b1*Bi(1,ixB  ,iyB  ,izB  ) + b2*Bi(1,ixB+1,iyB  ,izB  ) + &
                b3*Bi(1,ixB+1,iyB+1,izB  ) + b4*Bi(1,ixB  ,iyB+1,izB  ) + &
                b5*Bi(1,ixB  ,iyB  ,izB+1) + b6*Bi(1,ixB+1,iyB  ,izB+1) + &
                b7*Bi(1,ixB+1,iyB+1,izB+1) + b8*Bi(1,ixB  ,iyB+1,izB+1)

          Bpy = b1*Bi(2,ixB  ,iyB  ,izB  ) + b2*Bi(2,ixB+1,iyB  ,izB  ) + &
                b3*Bi(2,ixB+1,iyB+1,izB  ) + b4*Bi(2,ixB  ,iyB+1,izB  ) + &
                b5*Bi(2,ixB  ,iyB  ,izB+1) + b6*Bi(2,ixB+1,iyB  ,izB+1) + &
                b7*Bi(2,ixB+1,iyB+1,izB+1) + b8*Bi(2,ixB  ,iyB+1,izB+1)

          Bpz = b1*Bi(3,ixB  ,iyB  ,izB  ) + b2*Bi(3,ixB+1,iyB  ,izB  ) + &
                b3*Bi(3,ixB+1,iyB+1,izB  ) + b4*Bi(3,ixB  ,iyB+1,izB  ) + &
                b5*Bi(3,ixB  ,iyB  ,izB+1) + b6*Bi(3,ixB+1,iyB  ,izB+1) + &
                b7*Bi(3,ixB+1,iyB+1,izB+1) + b8*Bi(3,ixB  ,iyB+1,izB+1)
        else
          do_boris = .false.
        end if

      end if

      n_boris_total = n_boris_total + 1_int32
      if (do_boris) n_boris_used = n_boris_used + 1_int32

      ! ------------------------------------------------------------
      ! Boris push without sqrt(|B|^2) where B is present; otherwise
      ! the same electrostatic-only leapfrog update used by
      ! move_and_bc_electrostatic_fast above (mirrors legacy's own
      ! fallback, Src/part_expmover.f90:235-239).
      ! ------------------------------------------------------------
      if (do_boris) then

        vminus_x = part%pv(4,i) + qm2dt*Exp
        vminus_y = part%pv(5,i) + qm2dt*Eyp
        vminus_z = part%pv(6,i) + qm2dt*Ezp

        tx = qm2dt * Bpx
        ty = qm2dt * Bpy
        tz = qm2dt * Bpz

        t2        = tx*tx + ty*ty + tz*tz
        inv_denom = 1.0_real64 / (1.0_real64 + t2)

        sx = 2.0_real64 * tx * inv_denom
        sy = 2.0_real64 * ty * inv_denom
        sz = 2.0_real64 * tz * inv_denom

        vprime_x = vminus_x + (vminus_y*tz - vminus_z*ty)
        vprime_y = vminus_y + (vminus_z*tx - vminus_x*tz)
        vprime_z = vminus_z + (vminus_x*ty - vminus_y*tx)

        vplus_x = vminus_x + (vprime_y*sz - vprime_z*sy)
        vplus_y = vminus_y + (vprime_z*sx - vprime_x*sz)
        vplus_z = vminus_z + (vprime_x*sy - vprime_y*sx)

        vpx_new = vplus_x + qm2dt*Exp
        vpy_new = vplus_y + qm2dt*Eyp
        vpz_new = vplus_z + qm2dt*Ezp

      else

        vpx_new = part%pv(4,i) + 2.0_real64*qm2dt*Exp
        vpy_new = part%pv(5,i) + 2.0_real64*qm2dt*Eyp
        vpz_new = part%pv(6,i) + 2.0_real64*qm2dt*Ezp

      end if

      xp_new = xp + dt*vpx_new
      yp_new = yp + dt*vpy_new
      zp_new = zp + dt*vpz_new

      ! ---- RF antenna power absorbed (uses pre-update part%vy/vz(i)) ----
      if (in_rf_region) then
        v_theta = -0.5_real64*(part%pv(5,i)+vpy_new)*sin(theta_rf) + &
                   0.5_real64*(part%pv(6,i)+vpz_new)*cos(theta_rf)
        P_RF_local = P_RF_local + Nm_species*q*EAth*v_theta*dt

        ! Impulse delivered by the RF field alone (exact - no v-projection
        ! needed the way the power term above needs v_theta, since impulse
        ! is linear in E). No x-component: the prescribed RF field only has
        ! azimuthal (y,z) components (EAy,EAz above).
        mom_RF_local(2) = mom_RF_local(2) + Nm_species*q*EAy*dt
        mom_RF_local(3) = mom_RF_local(3) + Nm_species*q*EAz*dt
      end if

      ! ---- boundary conditions / SEE, operating on the just-pushed state ----
      ix = floor(xp_new / h(1)) + 1_int32
      iy = floor(yp_new / h(2)) + 1_int32
      iz = floor(zp_new / h(3)) + 1_int32

      if (ix < 0_int32)      ix = 0_int32
      if (ix > n(1)+1_int32) ix = n(1)+1_int32
      if (iy < 0_int32)      iy = 0_int32
      if (iy > n(2)+1_int32) iy = n(2)+1_int32
      if (iz < 0_int32)      iz = 0_int32
      if (iz > n(3)+1_int32) iz = n(3)+1_int32

      flag_lost = 0_int32

      if (particle_is_lost(wall_cell, ix, iy, iz, n)) flag_lost = 1_int32

      if (ptype == tag_neg) then
        if (xp_new < 0.0_real64 .and. flag_nmn == 1_int32) flag_lost = 2_int32
      end if

      if (flag_lost >= 1_int32) then

        igrid = bcnd(ix,iy,iz)

        if (flag_lost == 2_int32 .and. ptype == tag_neg) igrid = 0_int32

        if (flag_die == 1_int32 .and. igrid > 0_int32) then
          if (dtype(igrid) > 1_int32) then

            if (flag_pbc == 1_int32) then
              if (zp_new >= zmax) then
                zp_new = zp_new - zmax
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
              if (zp_new <= 0.0_real64) then
                zp_new = zmax + zp_new
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
            end if

            if (ix < 0_int32)      ix = 0_int32
            if (ix > n(1)+1_int32) ix = n(1)+1_int32
            if (iy < 0_int32)      iy = 0_int32
            if (iy > n(2)+1_int32) iy = n(2)+1_int32
            if (iz < 0_int32)      iz = 0_int32
            if (iz > n(3)+1_int32) iz = n(3)+1_int32

            pz = (real(iz,real64)*h(3) - zp_new) / h(3)

            if (dtype(igrid) == 2_int32) then
              px = (real(ix,real64)*h(1) - xp_new) / h(1)

              ki4(1) = qmacro * px               * pz
              ki4(2) = qmacro * (1.0_real64-px) * pz
              ki4(3) = qmacro * (1.0_real64-px) * (1.0_real64-pz)
              ki4(4) = qmacro * px               * (1.0_real64-pz)

              sum_q_xz_local(ix  ,iz  ) = sum_q_xz_local(ix  ,iz  ) + ki4(1)
              sum_q_xz_local(ix+1,iz  ) = sum_q_xz_local(ix+1,iz  ) + ki4(2)
              sum_q_xz_local(ix+1,iz+1) = sum_q_xz_local(ix+1,iz+1) + ki4(3)
              sum_q_xz_local(ix  ,iz+1) = sum_q_xz_local(ix  ,iz+1) + ki4(4)
            end if

            if (dtype(igrid) == 3_int32 .or. dtype(igrid) == 4_int32) then
              py = (real(iy,real64)*h(2) - yp_new) / h(2)

              ki4(1) = qmacro * py               * pz
              ki4(2) = qmacro * (1.0_real64-py) * pz
              ki4(3) = qmacro * (1.0_real64-py) * (1.0_real64-pz)
              ki4(4) = qmacro * py               * (1.0_real64-pz)

              d_ind = dtype(igrid)

              sum_q_yz_local(d_ind-2,iy  ,iz  ) = sum_q_yz_local(d_ind-2,iy  ,iz  ) + ki4(1)
              sum_q_yz_local(d_ind-2,iy+1,iz  ) = sum_q_yz_local(d_ind-2,iy+1,iz  ) + ki4(2)
              sum_q_yz_local(d_ind-2,iy+1,iz+1) = sum_q_yz_local(d_ind-2,iy+1,iz+1) + ki4(3)
              sum_q_yz_local(d_ind-2,iy  ,iz+1) = sum_q_yz_local(d_ind-2,iy  ,iz+1) + ki4(4)
            end if

          end if
        end if

        ! Wall diagnostics
        if (igrid < 0_int32) igrid = 0_int32

        Ek_eV = 0.5_real64 * m * &
          (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new) / qe

        p_mac_boundary(1,igrid) = p_mac_boundary(1,igrid) + qmacro
        p_mac_boundary(2,igrid) = p_mac_boundary(2,igrid) + abs(qmacro) * Ek_eV

        Ek_J = 0.5_real64 * m * &
            (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new)
        P_loss_wall = P_loss_wall + Nm_species * Ek_J

        mom_loss_wall(1) = mom_loss_wall(1) + Nm_species * m * vpx_new
        mom_loss_wall(2) = mom_loss_wall(2) + Nm_species * m * vpy_new
        mom_loss_wall(3) = mom_loss_wall(3) + Nm_species * m * vpz_new

        ! Secondary electron emission (positive ions hitting igrid_sec)
        if (do_see .and. igrid == see%igrid_sec) then
          rnd(1) = ran2(iseed)
          n_sec  = int(see%gam_sec, int32)
          if (rnd(1) <= (see%gam_sec - real(n_sec, real64))) n_sec = n_sec + 1_int32

          if (n_sec > 0_int32) then
            ! No lock needed: this routine's only caller (state%advance_particles_local,
            ! mod_state.f90) parallelizes over iproc alone, one thread owning
            ! this iproc's part_electrons = part(1,iproc) for the whole call.
            call part_electrons%ensure_capacity(part_electrons%n + n_sec)
            do ip_sec = 1, n_sec
              rnd(1) = ran2(iseed)
              vz_sec = -sign(1.0_real64, vpz_new) * see%vt_sec * sqrt(-log(1.0_real64 - rnd(1)))
              rnd(1) = ran2(iseed)
              rnd(2) = ran2(iseed)
              call load_gauss(vx_sec, vy_sec, see%vt_sec, rnd)

              part_electrons%n    = part_electrons%n + 1_int32
              i_see               = part_electrons%n
              part_electrons%pv(1,i_see)  = xp_new
              part_electrons%pv(2,i_see)  = yp_new
              part_electrons%pv(3,i_see)  = merge(see%zg_sec(1), see%zg_sec(2), vpz_new < 0.0_real64)
              part_electrons%pv(4,i_see) = vx_sec
              part_electrons%pv(5,i_see) = vy_sec
              part_electrons%pv(6,i_see) = vz_sec
              if (allocated(part_electrons%flag_dead)) part_electrons%flag_dead(i_see) = 0_int8
              if (allocated(part_electrons%flag_cex))  part_electrons%flag_cex(i_see)  = 0_int32

              P_loss_see = P_loss_see + 0.5_real64 * see%Nm_e * &
                  (vx_sec*vx_sec + vy_sec*vy_sec + vz_sec*vz_sec)

              mom_loss_see(1) = mom_loss_see(1) + see%Nm_e * see%mass_e * vx_sec
              mom_loss_see(2) = mom_loss_see(2) + see%Nm_e * see%mass_e * vy_sec
              mom_loss_see(3) = mom_loss_see(3) + see%Nm_e * see%mass_e * vz_sec
            end do
          end if
        end if

        np_lost = np_lost + 1_int32
        cycle
      end if

      ! Survivors: Neumann reflection on LHS only
      if (flag_nmn == 1_int32) then
        if (xp_new <= 0.0_real64) then
          xp_new  = -xp_new
          vpx_new = -vpx_new
        end if
      end if

      if (flag_pbc == 1_int32) then
        if (yp_new >= ymax) yp_new = yp_new - ymax
        if (yp_new <= 0.0_real64) yp_new = ymax + yp_new
        if (zp_new >= zmax) zp_new = zp_new - zmax
        if (zp_new <= 0.0_real64) zp_new = zmax + zp_new
      end if

      i_shift = i - np_lost

      part%pv(1,i_shift)  = xp_new
      part%pv(2,i_shift)  = yp_new
      part%pv(3,i_shift)  = zp_new
      part%pv(4,i_shift) = vpx_new
      part%pv(5,i_shift) = vpy_new
      part%pv(6,i_shift) = vpz_new

      if (has_w)         part%w(i_shift)         = part%w(i)
      if (has_sp)        part%sp(i_shift)        = part%sp(i)
      if (has_flag_dead) part%flag_dead(i_shift) = 0_int8
      if (has_flag_cex)  part%flag_cex(i_shift)  = part%flag_cex(i)

    end do

    part%n = part%n - np_lost
    if (part%n < 0_int32) part%n = 0_int32

    ! Keep the invariant "slots past n have zeroed flags" by clearing only
    ! the np_lost slots compaction just vacated. This used to zero the
    ! whole spare capacity plus cell_id/cell_count/cell_start (~19 MB of
    ! ncells-sized arrays per species per iproc) on every call - pure
    ! memory-bandwidth cost that dominated the push at high thread counts
    ! and is why it barely gained from hyperthreading. The cell lists are
    ! only read by collisions, right after sort_particles_local rebuilds
    ! them, so invalidating them here bought nothing.
    if (allocated(part%flag_dead)) then
      part%flag_dead(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int8
    end if
    if (allocated(part%flag_cex)) then
      part%flag_cex(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int32
    end if

  end subroutine move_and_bc_boris


  subroutine move_and_bc_boris_fast( part, n, h, E, n_B, h_B, Bi, q, m, dt, &
                                      bcnd, wall_cell, xmax, ymax, zmax, flag_pbc, flag_nmn, &
                                      ptype, tag_neg, flag_die, dtype, qmacro, &
                                      sum_q_xz_local, sum_q_yz_local, p_mac_boundary, &
                                      P_loss_wall, Nm_species, see, part_electrons, &
                                      iseed, P_loss_see, &
                                      mom_loss_wall, mom_loss_see, &
                                      P_loss_coll, mom_loss_coll, &
                                      n_boris_used, n_boris_total )
    ! Lean fast path for move_and_bc_boris, used whenever the energy-
    ! conserving pusher and RF-antenna heating are both off (mirrors
    ! move_and_bc_electrostatic_fast - see that routine's header comment
    ! for the measured win and, critically, why the BC/SEE tail below is a
    ! deliberate duplicate rather than a shared subroutine call: extracting
    ! it once already caused a 3-5x regression from real per-call overhead
    ! this file's -ipo build did not eliminate). A throwaway compile-time-
    ! constant-false test on this exact ITER case (RF and SEE both already
    ! inactive there) measured ~4-5% mover_push reduction before this
    ! routine existed - same effect confirmed here. do_see is NOT dropped
    ! (unlike RF): it's a per-species-config toggle independent of RF/
    ! energy-conserving mode and move_and_bc_electrostatic_fast keeps it
    ! too, so this mirrors that precedent rather than assuming SEE is
    ! always off just because this ITER case has it off.
    type(ParticleSet),  intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
    real(real64),       intent(in)    :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32),     intent(in)    :: n_B(3)
    real(real64),       intent(in)    :: h_B(3)
    real(real64),       intent(in)    :: Bi(4,0:n_B(1)+2,0:n_B(2)+2,0:n_B(3)+2)
    real(real64),       intent(in)    :: q, m, dt

    integer(int32),     intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    logical(1),         intent(in)    :: wall_cell(0:n(1)+1,0:n(2)+1,0:n(3)+1)
    real(real64),       intent(in)    :: xmax, ymax, zmax
    integer(int32),     intent(in)    :: flag_pbc, flag_nmn, ptype, tag_neg
    integer(int32),     intent(in)    :: flag_die
    integer(int32),     intent(in)    :: dtype(:)
    real(real64),       intent(in)    :: qmacro
    real(real64),       intent(inout) :: sum_q_xz_local(0:n(1)+2,0:n(3)+2)
    real(real64),       intent(inout) :: sum_q_yz_local(2,0:n(2)+2,0:n(3)+2)
    real(real64),       intent(inout) :: p_mac_boundary(:,:)
    real(real64),       intent(inout) :: P_loss_wall
    real(real64),       intent(in)    :: Nm_species
    type(SeeParams),    intent(in)    :: see
    type(ParticleSet),  intent(inout), optional :: part_electrons
    integer(int32),     intent(inout) :: iseed
    real(real64),       intent(inout) :: P_loss_see
    ! Momentum-conservation diagnostic - see move_and_bc_electrostatic's
    ! header comment.
    real(real64),       intent(inout) :: mom_loss_wall(3)
    real(real64),       intent(inout) :: mom_loss_see(3)
    real(real64),       intent(inout) :: P_loss_coll
    real(real64),       intent(inout) :: mom_loss_coll(3)
    ! TEMPORARY diagnostic - see move_and_bc_boris's header comment.
    integer(int32),     intent(out)   :: n_boris_used, n_boris_total

    integer(int32) :: i, i_shift, i_see, ip_sec, n_sec
    integer(int32) :: ix, iy, iz
    integer(int32) :: ixB, iyB, izB
    integer(int32) :: flag_lost, np_lost
    integer(int32) :: igrid, d_ind

    real(real64) :: qm2dt
    real(real64) :: xp, yp, zp
    real(real64) :: px, py, pz
    real(real64) :: wx2, wy2, wz2
    real(real64) :: w1, w2, w3, w4, w5, w6, w7, w8

    real(real64) :: pxB, pyB, pzB
    real(real64) :: wx2B, wy2B, wz2B
    real(real64) :: b1, b2, b3, b4, b5, b6, b7, b8

    real(real64) :: Exp, Eyp, Ezp
    real(real64) :: Bpx, Bpy, Bpz

    real(real64) :: vminus_x, vminus_y, vminus_z
    real(real64) :: vprime_x, vprime_y, vprime_z
    real(real64) :: vplus_x, vplus_y, vplus_z

    real(real64) :: tx, ty, tz
    real(real64) :: sx, sy, sz
    real(real64) :: t2, inv_denom

    logical :: uniform_B, same_grid_B, do_boris

    real(real64) :: xp_new, yp_new, zp_new
    real(real64) :: vpx_new, vpy_new, vpz_new
    real(real64) :: ki4(4)
    real(real64) :: Ek_eV, Ek_J
    real(real64) :: rnd(2)
    real(real64) :: vx_sec, vy_sec, vz_sec

    logical :: do_see

    logical :: has_w, has_sp, has_flag_dead, has_flag_cex

    n_boris_used  = 0_int32
    n_boris_total = 0_int32

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    qm2dt  = 0.5_real64 * dt * q / m
    ! part_electrons is OPTIONAL - see move_and_bc_electrostatic's header
    ! comment on why the electron-species call omits it entirely rather
    ! than passing self%part(1,iproc) aliased against `part`.
    do_see = present(part_electrons) .and. (see%gam_sec > 0.0_real64) .and. (q > 0.0_real64)

    ! See move_and_bc_electrostatic's header comment on has_w/has_sp/
    ! has_flag_dead/has_flag_cex.
    has_w         = allocated(part%w)
    has_sp        = allocated(part%sp)
    has_flag_dead = allocated(part%flag_dead)
    has_flag_cex  = allocated(part%flag_cex)

    uniform_B = (n_B(1) == 1_int32 .and. n_B(2) == 1_int32 .and. n_B(3) == 1_int32)

    same_grid_B = .false.
    if (.not. uniform_B) then
      same_grid_B = all(n_B == n) .and. &
                    abs(h_B(1)-h(1)) < 1.0e-14_real64 .and. &
                    abs(h_B(2)-h(2)) < 1.0e-14_real64 .and. &
                    abs(h_B(3)-h(3)) < 1.0e-14_real64
    end if

    np_lost = 0_int32

    do i = 1_int32, part%n

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) /= 0_int8) then
          np_lost = np_lost + 1_int32
          ! Remove the energy/momentum this particle carried when a
          ! collision flagged it dead last step (mirrors legacy
          ! Src/part_expmover.f90:70-84) - otherwise it stays counted in
          ! Pcoll/mom_loss(:,3,:) after this step quietly discards it.
          P_loss_coll = P_loss_coll - Nm_species * 0.5_real64 * m * &
              (part%pv(4,i)*part%pv(4,i) + part%pv(5,i)*part%pv(5,i) + part%pv(6,i)*part%pv(6,i))
          mom_loss_coll(1) = mom_loss_coll(1) - Nm_species * m * part%pv(4,i)
          mom_loss_coll(2) = mom_loss_coll(2) - Nm_species * m * part%pv(5,i)
          mom_loss_coll(3) = mom_loss_coll(3) - Nm_species * m * part%pv(6,i)
          cycle
        end if
      end if

      ! ---- push ----
      xp = part%pv(1,i)
      yp = part%pv(2,i)
      zp = part%pv(3,i)

      ix = int(xp / h(1), int32) + 1_int32
      iy = int(yp / h(2), int32) + 1_int32
      iz = int(zp / h(3), int32) + 1_int32

      ix = max(0_int32, min(n(1)+1_int32, ix))
      iy = max(0_int32, min(n(2)+1_int32, iy))
      iz = max(0_int32, min(n(3)+1_int32, iz))

      px = (real(ix, real64)*h(1) - xp) / h(1)
      py = (real(iy, real64)*h(2) - yp) / h(2)
      pz = (real(iz, real64)*h(3) - zp) / h(3)

      wx2 = 1.0_real64 - px
      wy2 = 1.0_real64 - py
      wz2 = 1.0_real64 - pz

      w1 = px  * py  * pz
      w2 = wx2 * py  * pz
      w3 = wx2 * wy2 * pz
      w4 = px  * wy2 * pz
      w5 = px  * py  * wz2
      w6 = wx2 * py  * wz2
      w7 = wx2 * wy2 * wz2
      w8 = px  * wy2 * wz2

      Exp = w1*E(1,ix  ,iy  ,iz  ) + w2*E(1,ix+1,iy  ,iz  ) + &
            w3*E(1,ix+1,iy+1,iz  ) + w4*E(1,ix  ,iy+1,iz  ) + &
            w5*E(1,ix  ,iy  ,iz+1) + w6*E(1,ix+1,iy  ,iz+1) + &
            w7*E(1,ix+1,iy+1,iz+1) + w8*E(1,ix  ,iy+1,iz+1)

      Eyp = w1*E(2,ix  ,iy  ,iz  ) + w2*E(2,ix+1,iy  ,iz  ) + &
            w3*E(2,ix+1,iy+1,iz  ) + w4*E(2,ix  ,iy+1,iz  ) + &
            w5*E(2,ix  ,iy  ,iz+1) + w6*E(2,ix+1,iy  ,iz+1) + &
            w7*E(2,ix+1,iy+1,iz+1) + w8*E(2,ix  ,iy+1,iz+1)

      Ezp = w1*E(3,ix  ,iy  ,iz  ) + w2*E(3,ix+1,iy  ,iz  ) + &
            w3*E(3,ix+1,iy+1,iz  ) + w4*E(3,ix  ,iy+1,iz  ) + &
            w5*E(3,ix  ,iy  ,iz+1) + w6*E(3,ix+1,iy  ,iz+1) + &
            w7*E(3,ix+1,iy+1,iz+1) + w8*E(3,ix  ,iy+1,iz+1)

      ! ------------------------------------------------------------
      ! B-field interpolation - skipped where |B| is negligible.
      ! See move_and_bc_boris's header comment on this gate.
      ! ------------------------------------------------------------
      do_boris = .true.

      if (uniform_B) then

        Bpx = Bi(1,1,1,1)
        Bpy = Bi(2,1,1,1)
        Bpz = Bi(3,1,1,1)

      else if (same_grid_B) then

        if (Bi(4,ix,iy,iz) > 1.0e-5_real64) then
          Bpx = w1*Bi(1,ix  ,iy  ,iz  ) + w2*Bi(1,ix+1,iy  ,iz  ) + &
                w3*Bi(1,ix+1,iy+1,iz  ) + w4*Bi(1,ix  ,iy+1,iz  ) + &
                w5*Bi(1,ix  ,iy  ,iz+1) + w6*Bi(1,ix+1,iy  ,iz+1) + &
                w7*Bi(1,ix+1,iy+1,iz+1) + w8*Bi(1,ix  ,iy+1,iz+1)

          Bpy = w1*Bi(2,ix  ,iy  ,iz  ) + w2*Bi(2,ix+1,iy  ,iz  ) + &
                w3*Bi(2,ix+1,iy+1,iz  ) + w4*Bi(2,ix  ,iy+1,iz  ) + &
                w5*Bi(2,ix  ,iy  ,iz+1) + w6*Bi(2,ix+1,iy  ,iz+1) + &
                w7*Bi(2,ix+1,iy+1,iz+1) + w8*Bi(2,ix  ,iy+1,iz+1)

          Bpz = w1*Bi(3,ix  ,iy  ,iz  ) + w2*Bi(3,ix+1,iy  ,iz  ) + &
                w3*Bi(3,ix+1,iy+1,iz  ) + w4*Bi(3,ix  ,iy+1,iz  ) + &
                w5*Bi(3,ix  ,iy  ,iz+1) + w6*Bi(3,ix+1,iy  ,iz+1) + &
                w7*Bi(3,ix+1,iy+1,iz+1) + w8*Bi(3,ix  ,iy+1,iz+1)
        else
          do_boris = .false.
        end if

      else

        ixB = int(xp / h_B(1), int32) + 1_int32
        iyB = int(yp / h_B(2), int32) + 1_int32
        izB = int(zp / h_B(3), int32) + 1_int32

        ixB = max(0_int32, min(n_B(1)+1_int32, ixB))
        iyB = max(0_int32, min(n_B(2)+1_int32, iyB))
        izB = max(0_int32, min(n_B(3)+1_int32, izB))

        if (Bi(4,ixB,iyB,izB) > 1.0e-5_real64) then
          pxB = (real(ixB, real64)*h_B(1) - xp) / h_B(1)
          pyB = (real(iyB, real64)*h_B(2) - yp) / h_B(2)
          pzB = (real(izB, real64)*h_B(3) - zp) / h_B(3)

          wx2B = 1.0_real64 - pxB
          wy2B = 1.0_real64 - pyB
          wz2B = 1.0_real64 - pzB

          b1 = pxB  * pyB  * pzB
          b2 = wx2B * pyB  * pzB
          b3 = wx2B * wy2B * pzB
          b4 = pxB  * wy2B * pzB
          b5 = pxB  * pyB  * wz2B
          b6 = wx2B * pyB  * wz2B
          b7 = wx2B * wy2B * wz2B
          b8 = pxB  * wy2B * wz2B

          Bpx = b1*Bi(1,ixB  ,iyB  ,izB  ) + b2*Bi(1,ixB+1,iyB  ,izB  ) + &
                b3*Bi(1,ixB+1,iyB+1,izB  ) + b4*Bi(1,ixB  ,iyB+1,izB  ) + &
                b5*Bi(1,ixB  ,iyB  ,izB+1) + b6*Bi(1,ixB+1,iyB  ,izB+1) + &
                b7*Bi(1,ixB+1,iyB+1,izB+1) + b8*Bi(1,ixB  ,iyB+1,izB+1)

          Bpy = b1*Bi(2,ixB  ,iyB  ,izB  ) + b2*Bi(2,ixB+1,iyB  ,izB  ) + &
                b3*Bi(2,ixB+1,iyB+1,izB  ) + b4*Bi(2,ixB  ,iyB+1,izB  ) + &
                b5*Bi(2,ixB  ,iyB  ,izB+1) + b6*Bi(2,ixB+1,iyB  ,izB+1) + &
                b7*Bi(2,ixB+1,iyB+1,izB+1) + b8*Bi(2,ixB  ,iyB+1,izB+1)

          Bpz = b1*Bi(3,ixB  ,iyB  ,izB  ) + b2*Bi(3,ixB+1,iyB  ,izB  ) + &
                b3*Bi(3,ixB+1,iyB+1,izB  ) + b4*Bi(3,ixB  ,iyB+1,izB  ) + &
                b5*Bi(3,ixB  ,iyB  ,izB+1) + b6*Bi(3,ixB+1,iyB  ,izB+1) + &
                b7*Bi(3,ixB+1,iyB+1,izB+1) + b8*Bi(3,ixB  ,iyB+1,izB+1)
        else
          do_boris = .false.
        end if

      end if

      n_boris_total = n_boris_total + 1_int32
      if (do_boris) n_boris_used = n_boris_used + 1_int32

      ! ------------------------------------------------------------
      ! Boris push without sqrt(|B|^2) where B is present; otherwise
      ! the same electrostatic-only leapfrog update used by
      ! move_and_bc_electrostatic_fast above (mirrors legacy's own
      ! fallback, Src/part_expmover.f90:235-239).
      ! ------------------------------------------------------------
      if (do_boris) then

        vminus_x = part%pv(4,i) + qm2dt*Exp
        vminus_y = part%pv(5,i) + qm2dt*Eyp
        vminus_z = part%pv(6,i) + qm2dt*Ezp

        tx = qm2dt * Bpx
        ty = qm2dt * Bpy
        tz = qm2dt * Bpz

        t2        = tx*tx + ty*ty + tz*tz
        inv_denom = 1.0_real64 / (1.0_real64 + t2)

        sx = 2.0_real64 * tx * inv_denom
        sy = 2.0_real64 * ty * inv_denom
        sz = 2.0_real64 * tz * inv_denom

        vprime_x = vminus_x + (vminus_y*tz - vminus_z*ty)
        vprime_y = vminus_y + (vminus_z*tx - vminus_x*tz)
        vprime_z = vminus_z + (vminus_x*ty - vminus_y*tx)

        vplus_x = vminus_x + (vprime_y*sz - vprime_z*sy)
        vplus_y = vminus_y + (vprime_z*sx - vprime_x*sz)
        vplus_z = vminus_z + (vprime_x*sy - vprime_y*sx)

        vpx_new = vplus_x + qm2dt*Exp
        vpy_new = vplus_y + qm2dt*Eyp
        vpz_new = vplus_z + qm2dt*Ezp

      else

        vpx_new = part%pv(4,i) + 2.0_real64*qm2dt*Exp
        vpy_new = part%pv(5,i) + 2.0_real64*qm2dt*Eyp
        vpz_new = part%pv(6,i) + 2.0_real64*qm2dt*Ezp

      end if

      xp_new = xp + dt*vpx_new
      yp_new = yp + dt*vpy_new
      zp_new = zp + dt*vpz_new

      ! ---- boundary conditions / SEE, operating on the just-pushed state ----
      ix = floor(xp_new / h(1)) + 1_int32
      iy = floor(yp_new / h(2)) + 1_int32
      iz = floor(zp_new / h(3)) + 1_int32

      if (ix < 0_int32)      ix = 0_int32
      if (ix > n(1)+1_int32) ix = n(1)+1_int32
      if (iy < 0_int32)      iy = 0_int32
      if (iy > n(2)+1_int32) iy = n(2)+1_int32
      if (iz < 0_int32)      iz = 0_int32
      if (iz > n(3)+1_int32) iz = n(3)+1_int32

      flag_lost = 0_int32

      if (particle_is_lost(wall_cell, ix, iy, iz, n)) flag_lost = 1_int32

      if (ptype == tag_neg) then
        if (xp_new < 0.0_real64 .and. flag_nmn == 1_int32) flag_lost = 2_int32
      end if

      if (flag_lost >= 1_int32) then

        igrid = bcnd(ix,iy,iz)

        if (flag_lost == 2_int32 .and. ptype == tag_neg) igrid = 0_int32

        if (flag_die == 1_int32 .and. igrid > 0_int32) then
          if (dtype(igrid) > 1_int32) then

            if (flag_pbc == 1_int32) then
              if (zp_new >= zmax) then
                zp_new = zp_new - zmax
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
              if (zp_new <= 0.0_real64) then
                zp_new = zmax + zp_new
                iz = floor(zp_new / h(3), int32) + 1_int32
              end if
            end if

            if (ix < 0_int32)      ix = 0_int32
            if (ix > n(1)+1_int32) ix = n(1)+1_int32
            if (iy < 0_int32)      iy = 0_int32
            if (iy > n(2)+1_int32) iy = n(2)+1_int32
            if (iz < 0_int32)      iz = 0_int32
            if (iz > n(3)+1_int32) iz = n(3)+1_int32

            pz = (real(iz,real64)*h(3) - zp_new) / h(3)

            if (dtype(igrid) == 2_int32) then
              px = (real(ix,real64)*h(1) - xp_new) / h(1)

              ki4(1) = qmacro * px               * pz
              ki4(2) = qmacro * (1.0_real64-px) * pz
              ki4(3) = qmacro * (1.0_real64-px) * (1.0_real64-pz)
              ki4(4) = qmacro * px               * (1.0_real64-pz)

              sum_q_xz_local(ix  ,iz  ) = sum_q_xz_local(ix  ,iz  ) + ki4(1)
              sum_q_xz_local(ix+1,iz  ) = sum_q_xz_local(ix+1,iz  ) + ki4(2)
              sum_q_xz_local(ix+1,iz+1) = sum_q_xz_local(ix+1,iz+1) + ki4(3)
              sum_q_xz_local(ix  ,iz+1) = sum_q_xz_local(ix  ,iz+1) + ki4(4)
            end if

            if (dtype(igrid) == 3_int32 .or. dtype(igrid) == 4_int32) then
              py = (real(iy,real64)*h(2) - yp_new) / h(2)

              ki4(1) = qmacro * py               * pz
              ki4(2) = qmacro * (1.0_real64-py) * pz
              ki4(3) = qmacro * (1.0_real64-py) * (1.0_real64-pz)
              ki4(4) = qmacro * py               * (1.0_real64-pz)

              d_ind = dtype(igrid)

              sum_q_yz_local(d_ind-2,iy  ,iz  ) = sum_q_yz_local(d_ind-2,iy  ,iz  ) + ki4(1)
              sum_q_yz_local(d_ind-2,iy+1,iz  ) = sum_q_yz_local(d_ind-2,iy+1,iz  ) + ki4(2)
              sum_q_yz_local(d_ind-2,iy+1,iz+1) = sum_q_yz_local(d_ind-2,iy+1,iz+1) + ki4(3)
              sum_q_yz_local(d_ind-2,iy  ,iz+1) = sum_q_yz_local(d_ind-2,iy  ,iz+1) + ki4(4)
            end if

          end if
        end if

        ! Wall diagnostics
        if (igrid < 0_int32) igrid = 0_int32

        Ek_eV = 0.5_real64 * m * &
          (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new) / qe

        p_mac_boundary(1,igrid) = p_mac_boundary(1,igrid) + qmacro
        p_mac_boundary(2,igrid) = p_mac_boundary(2,igrid) + abs(qmacro) * Ek_eV

        Ek_J = 0.5_real64 * m * &
            (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new)
        P_loss_wall = P_loss_wall + Nm_species * Ek_J

        mom_loss_wall(1) = mom_loss_wall(1) + Nm_species * m * vpx_new
        mom_loss_wall(2) = mom_loss_wall(2) + Nm_species * m * vpy_new
        mom_loss_wall(3) = mom_loss_wall(3) + Nm_species * m * vpz_new

        ! Secondary electron emission (positive ions hitting igrid_sec)
        if (do_see .and. igrid == see%igrid_sec) then
          rnd(1) = ran2(iseed)
          n_sec  = int(see%gam_sec, int32)
          if (rnd(1) <= (see%gam_sec - real(n_sec, real64))) n_sec = n_sec + 1_int32

          if (n_sec > 0_int32) then
            ! No lock needed: this routine's only caller (state%advance_particles_local,
            ! mod_state.f90) parallelizes over iproc alone, one thread owning
            ! this iproc's part_electrons = part(1,iproc) for the whole call.
            call part_electrons%ensure_capacity(part_electrons%n + n_sec)
            do ip_sec = 1, n_sec
              rnd(1) = ran2(iseed)
              vz_sec = -sign(1.0_real64, vpz_new) * see%vt_sec * sqrt(-log(1.0_real64 - rnd(1)))
              rnd(1) = ran2(iseed)
              rnd(2) = ran2(iseed)
              call load_gauss(vx_sec, vy_sec, see%vt_sec, rnd)

              part_electrons%n    = part_electrons%n + 1_int32
              i_see               = part_electrons%n
              part_electrons%pv(1,i_see)  = xp_new
              part_electrons%pv(2,i_see)  = yp_new
              part_electrons%pv(3,i_see)  = merge(see%zg_sec(1), see%zg_sec(2), vpz_new < 0.0_real64)
              part_electrons%pv(4,i_see) = vx_sec
              part_electrons%pv(5,i_see) = vy_sec
              part_electrons%pv(6,i_see) = vz_sec
              if (allocated(part_electrons%flag_dead)) part_electrons%flag_dead(i_see) = 0_int8
              if (allocated(part_electrons%flag_cex))  part_electrons%flag_cex(i_see)  = 0_int32

              P_loss_see = P_loss_see + 0.5_real64 * see%Nm_e * &
                  (vx_sec*vx_sec + vy_sec*vy_sec + vz_sec*vz_sec)

              mom_loss_see(1) = mom_loss_see(1) + see%Nm_e * see%mass_e * vx_sec
              mom_loss_see(2) = mom_loss_see(2) + see%Nm_e * see%mass_e * vy_sec
              mom_loss_see(3) = mom_loss_see(3) + see%Nm_e * see%mass_e * vz_sec
            end do
          end if
        end if

        np_lost = np_lost + 1_int32
        cycle
      end if

      ! Survivors: Neumann reflection on LHS only
      if (flag_nmn == 1_int32) then
        if (xp_new <= 0.0_real64) then
          xp_new  = -xp_new
          vpx_new = -vpx_new
        end if
      end if

      if (flag_pbc == 1_int32) then
        if (yp_new >= ymax) yp_new = yp_new - ymax
        if (yp_new <= 0.0_real64) yp_new = ymax + yp_new
        if (zp_new >= zmax) zp_new = zp_new - zmax
        if (zp_new <= 0.0_real64) zp_new = zmax + zp_new
      end if

      i_shift = i - np_lost

      part%pv(1,i_shift)  = xp_new
      part%pv(2,i_shift)  = yp_new
      part%pv(3,i_shift)  = zp_new
      part%pv(4,i_shift) = vpx_new
      part%pv(5,i_shift) = vpy_new
      part%pv(6,i_shift) = vpz_new

      if (has_w)         part%w(i_shift)         = part%w(i)
      if (has_sp)        part%sp(i_shift)        = part%sp(i)
      if (has_flag_dead) part%flag_dead(i_shift) = 0_int8
      if (has_flag_cex)  part%flag_cex(i_shift)  = part%flag_cex(i)

    end do

    part%n = part%n - np_lost
    if (part%n < 0_int32) part%n = 0_int32

    ! Keep the invariant "slots past n have zeroed flags" by clearing only
    ! the np_lost slots compaction just vacated. This used to zero the
    ! whole spare capacity plus cell_id/cell_count/cell_start (~19 MB of
    ! ncells-sized arrays per species per iproc) on every call - pure
    ! memory-bandwidth cost that dominated the push at high thread counts
    ! and is why it barely gained from hyperthreading. The cell lists are
    ! only read by collisions, right after sort_particles_local rebuilds
    ! them, so invalidating them here bought nothing.
    if (allocated(part%flag_dead)) then
      part%flag_dead(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int8
    end if
    if (allocated(part%flag_cex)) then
      part%flag_cex(part%n+1:min(part%nmax, part%n+np_lost)) = 0_int32
    end if

  end subroutine move_and_bc_boris_fast

end module mod_particleMover