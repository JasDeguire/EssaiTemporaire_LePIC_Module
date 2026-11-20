module mod_collisions
  use iso_fortran_env,    only: int32, real64
  use mod_particles,      only: ParticleSet
  use mod_collisionDiagnostics, only: init_rxn_counts
  use mod_collisionsGwenael, only: perform_collisions_gwenael
  use mod_CoulombCollisions, only: perform_coulomb_collisions

  implicit none
  private
  public :: perform_collisions_step, perform_coulomb_step

contains

  subroutine perform_collisions_step( &
      part, n, h, ntype_tracked, ntype_all, mass, charge, Ti, Nm, p_ncol, sig_list, col_info, &
      sigv_mx, sig, sig_Er, sig_Eex, ni0, ns_coll, dt, nu_uplim, iseed, &
      mpi_rank, Pcoll, dom_volume, np_red, bcnd, &
      ix_plane, iy_plane, iz_plane, sour_xy, sour_xz, sour_yz, &
      sink_xy, sink_xz, sink_yz, mom_loss)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
    integer(int32), intent(in) :: ntype_tracked, ntype_all
    real(real64), intent(in) :: mass(:), charge(:), Ti(:), Nm(:)
    integer(int32), intent(in) :: p_ncol(:)
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sigv_mx(:,:), sig(:,:), sig_Er(:), sig_Eex(:,:)
    real(real64), intent(in) :: ni0(:)
    integer(int32), intent(in) :: ns_coll
    real(real64), intent(in) :: dt, nu_uplim(:)
    integer(int32), intent(inout) :: iseed(:)
    integer(int32), intent(in) :: mpi_rank
    real(real64), intent(inout) :: Pcoll(:,:)
    real(real64), intent(inout) :: mom_loss(:,:,:)
    real(real64), intent(in) :: dom_volume
    real(real64), intent(in) :: np_red(0:,0:,0:,:)
    integer(int32), intent(in) :: bcnd(0:,0:,0:)
    integer(int32), intent(in) :: ix_plane, iy_plane, iz_plane
    real(real64), intent(inout) :: sour_xy(0:,0:,:,:), sour_xz(0:,0:,:,:), sour_yz(0:,0:,:,:)
    real(real64), intent(inout) :: sink_xy(0:,0:,:,:), sink_xz(0:,0:,:,:), sink_yz(0:,0:,:,:)

    call init_rxn_counts(int(size(sigv_mx,2), int32))

    call perform_collisions_gwenael( &
        part          = part, &
        n             = n, &
        h             = h, &
        ntype_tracked = ntype_tracked, &
        ntype_all     = ntype_all, &
        mass          = mass, &
        Ti            = Ti, &
        Nm            = Nm, &
        p_ncol        = p_ncol, &
        sig_list      = sig_list, &
        col_info      = col_info, &
        sigv_mx       = sigv_mx, &
        sig           = sig, &
        sig_Er        = sig_Er, &
        sig_Eex       = sig_Eex, &
        ni0           = ni0, &
        ns_coll       = ns_coll, &
        dt            = dt, &
        nu_uplim      = nu_uplim, &
        iseed         = iseed, &
        mpi_rank      = mpi_rank, &
        Pcoll         = Pcoll, &
        mom_loss      = mom_loss, &
        dom_volume    = dom_volume, &
        np_red        = np_red, &
        bcnd          = bcnd, &
        ix_plane      = ix_plane, &
        iy_plane      = iy_plane, &
        iz_plane      = iz_plane, &
        sour_xy       = sour_xy, &
        sour_xz       = sour_xz, &
        sour_yz       = sour_yz, &
        sink_xy       = sink_xy, &
        sink_xz       = sink_xz, &
        sink_yz       = sink_yz)

  end subroutine perform_collisions_step


  ! Coulomb small-angle scattering (Nanbu 2000), with its own cadence
  ! independent of MC collisions. dt must already be scaled by ns_coulomb
  ! (i.e. pass dt * nb_step_coulomb from the caller).
  subroutine perform_coulomb_step(part, ntype, nproc, mass, charge, Nm, n_e, dt, iseed)
    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32),    intent(in)    :: ntype, nproc
    real(real64),      intent(in)    :: mass(:), charge(:), Nm(:)
    real(real64),      intent(in)    :: n_e    ! [m^-3] global average electron density
    real(real64),      intent(in)    :: dt     ! effective dt = base_dt * nb_step_coulomb
    integer(int32),    intent(inout) :: iseed(:)

    call perform_coulomb_collisions( &
        part   = part, &
        ntype  = ntype, &
        nproc  = nproc, &
        mass   = mass(1:ntype), &
        charge = charge(1:ntype), &
        Nm     = Nm(1:ntype), &
        n_e    = n_e, &
        dt     = dt, &
        iseed  = iseed)

  end subroutine perform_coulomb_step

end module mod_collisions
