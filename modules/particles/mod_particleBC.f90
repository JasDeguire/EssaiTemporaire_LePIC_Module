module mod_particleBC
  use iso_fortran_env,  only: int32, int8, real64
  use mod_particles,    only: ParticleSet
  use mod_constants,    only: qe
  use mod_rng,          only: ran2, load_gauss

  implicit none
  private

  public :: apply_particle_bc
  public :: particle_is_lost
  public :: SeeParams

  ! Secondary electron emission parameters.
  ! Set gam_sec <= 0 (or charge_species <= 0) to disable SEE for a species.
  ! vt_sec = sqrt(2*qe*|THm|/|m_e|) — precomputed by the caller.
  type :: SeeParams
    real(real64)   :: gam_sec   = 0.0_real64
    integer(int32) :: igrid_sec = 0_int32
    real(real64)   :: zg_sec(2) = 0.0_real64
    real(real64)   :: vt_sec    = 0.0_real64
    real(real64)   :: Nm_e      = 0.0_real64
    real(real64)   :: mass_e    = 0.0_real64
  end type SeeParams

contains

  logical function particle_is_lost(wall_cell, ix, iy, iz, n) result(is_lost)
    integer(int32), intent(in) :: ix, iy, iz
    integer(int32), intent(in) :: n(3)
    logical(1),     intent(in) :: wall_cell(0:n(1)+1,0:n(2)+1,0:n(3)+1)

    ! wall_cell(ix,iy,iz) is precomputed once (mod_boundary.f90::build_boundary)
    ! as the exact 8-corner-of-bcnd-all->=1 predicate this used to compute here
    ! on every call -- bcnd is static after init, so this is a cache lookup now.
    is_lost = wall_cell(ix,iy,iz)
  end function particle_is_lost


  subroutine apply_particle_bc( part, n, h, bcnd, wall_cell, xmax, ymax, zmax, &
                                flag_pbc, flag_nmn, ptype, tag_neg,  &
                                flag_die, dtype, qmacro,              &
                                sum_q_xz_local, sum_q_yz_local,       &
                                p_mac_boundary, mass_species,          &
                                P_loss_wall, Nm_species,               &
                                charge_species, see,                   &
                                part_electrons, iseed, P_loss_see,     &
                                mom_loss_wall, mom_loss_see,           &
                                P_loss_coll, mom_loss_coll )

    type(ParticleSet),  intent(inout) :: part
    integer(int32),     intent(in)    :: n(3)
    real(real64),       intent(in)    :: h(3)
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
    real(real64),       intent(in)    :: mass_species
    real(real64),       intent(inout) :: P_loss_wall
    real(real64),       intent(in)    :: Nm_species
    real(real64),       intent(in)    :: charge_species
    type(SeeParams),    intent(in)    :: see
    type(ParticleSet),  intent(inout) :: part_electrons
    integer(int32),     intent(inout) :: iseed
    real(real64),       intent(inout) :: P_loss_see
    ! Momentum-conservation diagnostic - vector counterparts of P_loss_wall/
    ! P_loss_see, filled at the same call sites below. mom_loss_coll/
    ! P_loss_coll additionally close a pre-existing gap: a particle a
    ! collision flagged dead last step is discarded below with no energy/
    ! momentum bookkeeping at all (unlike legacy Src/part_expmover.f90:70-84,
    ! which subtracts it from P_loss(3)) - fixed here so Pcoll/mom_loss(:,3,:)
    ! actually balance instead of quietly leaking the dead particle's energy
    ! and momentum out of the books.
    real(real64),       intent(inout) :: mom_loss_wall(3)
    real(real64),       intent(inout) :: mom_loss_see(3)
    real(real64),       intent(inout) :: P_loss_coll
    real(real64),       intent(inout) :: mom_loss_coll(3)

    integer(int32) :: i, i_shift, i_see, ip_sec, n_sec
    integer(int32) :: ix, iy, iz
    integer(int32) :: flag_lost, np_lost
    integer(int32) :: igrid, d_ind

    real(real64) :: xp_new, yp_new, zp_new
    real(real64) :: vpx_new, vpy_new, vpz_new
    real(real64) :: px, py, pz
    real(real64) :: ki4(4)
    real(real64) :: Ek_eV, Ek_J
    real(real64) :: rnd(2)
    real(real64) :: vx_sec, vy_sec, vz_sec

    logical :: do_see

    if (.not. allocated(part%pv)) return
    if (part%n <= 0_int32) return

    do_see = (see%gam_sec > 0.0_real64) .and. (charge_species > 0.0_real64)

    np_lost = 0_int32

    do i = 1, part%n

      xp_new  = part%pv(1,i)
      yp_new  = part%pv(2,i)
      zp_new  = part%pv(3,i)
      vpx_new = part%pv(4,i)
      vpy_new = part%pv(5,i)
      vpz_new = part%pv(6,i)

      if (allocated(part%flag_dead)) then
        if (part%flag_dead(i) == 1_int8) then
          np_lost = np_lost + 1_int32
          P_loss_coll = P_loss_coll - Nm_species * 0.5_real64 * mass_species * &
              (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new)
          mom_loss_coll(1) = mom_loss_coll(1) - Nm_species * mass_species * vpx_new
          mom_loss_coll(2) = mom_loss_coll(2) - Nm_species * mass_species * vpy_new
          mom_loss_coll(3) = mom_loss_coll(3) - Nm_species * mass_species * vpz_new
          cycle
        end if
      end if

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

        Ek_eV = 0.5_real64 * mass_species * &
          (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new) / qe

        p_mac_boundary(1,igrid) = p_mac_boundary(1,igrid) + qmacro
        p_mac_boundary(2,igrid) = p_mac_boundary(2,igrid) + abs(qmacro) * Ek_eV

        Ek_J = 0.5_real64 * mass_species * &
            (vpx_new*vpx_new + vpy_new*vpy_new + vpz_new*vpz_new)
        P_loss_wall = P_loss_wall + Nm_species * Ek_J

        mom_loss_wall(1) = mom_loss_wall(1) + Nm_species * mass_species * vpx_new
        mom_loss_wall(2) = mom_loss_wall(2) + Nm_species * mass_species * vpy_new
        mom_loss_wall(3) = mom_loss_wall(3) + Nm_species * mass_species * vpz_new

        ! Secondary electron emission (positive ions hitting igrid_sec)
        if (do_see .and. igrid == see%igrid_sec) then
          rnd(1) = ran2(iseed)
          n_sec  = int(see%gam_sec, int32)
          if (rnd(1) <= (see%gam_sec - real(n_sec, real64))) n_sec = n_sec + 1_int32

          if (n_sec > 0_int32) then
            ! No lock needed: apply_particle_bc's only caller
            ! (state%advance_particles_local, mod_state.f90) parallelizes
            ! over iproc alone, one thread owning this iproc's
            ! part_electrons = part(1,iproc) for the whole call - so no
            ! other thread can ever touch it concurrently. A prior version
            ! of that caller parallelized over (ptype,iproc) collapsed
            ! together, which really could race here across ptypes of the
            ! same iproc; this critical section is a leftover from that.
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

      if (allocated(part%w))         part%w(i_shift)         = part%w(i)
      if (allocated(part%sp))        part%sp(i_shift)        = part%sp(i)
      if (allocated(part%flag_dead)) part%flag_dead(i_shift) = 0_int8
      if (allocated(part%flag_cex))  part%flag_cex(i_shift)  = part%flag_cex(i)

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

  end subroutine apply_particle_bc

end module mod_particleBC
