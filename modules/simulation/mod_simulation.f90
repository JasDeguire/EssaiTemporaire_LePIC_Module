module mod_simulation
  use iso_fortran_env, only: int32, real64, int8, int64

  use mod_state,                 only: State, advance_particles_local
  use mod_density,               only: reduce_species_density, reduce_density_and_rho, &
                                       build_rho_from_np, average_species_density
  use mod_restart,               only: write_restart_modular
  use mod_injection,             only: inject_particles_volume, inject_flux_particles
  use mod_PoissonSolver_legacy,  only: solve_poisson_legacy, print_poisson_breakdown, reset_poisson_breakdown
  use mod_electricField,         only: calc_Efield_modular, calc_Efield_energy_conserving
  use mod_output_2d,             only: write_density_planes, write_scalar_planes, &
                                       write_vector_component_planes, &
                                       write_plane_xy_scalar_2d, write_plane_xz_scalar_2d, &
                                       write_plane_yz_scalar_2d, &
                                       write_plane_xy_flux_2d, write_plane_xz_flux_2d, &
                                       write_plane_yz_flux_2d
  use mod_constants,             only: eps0, qe
  use mod_collisionDiagnostics, only: print_rxn_counts, reset_rxn_counts, &
                                       print_debug_diagnostics, reset_debug_diagnostics
  use mod_collisions, only: perform_collisions_step, perform_coulomb_step
  use mod_intro,      only: print_introduction
  use mpi

  implicit none
  private
  public :: Simulation

  logical, parameter :: DEBUG_SIMULATION = .false.

  ! TRIED, MEASURED, REJECTED (keep .false.): sort particles by cell every
  ! step, immediately before the fused mover+BC+deposit region, instead of
  ! only right before collisions (mod_state%sort_particles_local's other
  ! call site, ~line 627 below). Rationale that motivated the experiment:
  ! deposit_particle_set_to_np_thread's charge scatter and the mover's E/B
  ! gather (mod_particleMover.f90) both touch the field grid at
  ! particle-position-dependent, effectively random cell indices unless
  ! particles are already cell-sorted, and most steps run that scatter/
  ! gather on whatever order injection/the previous push left particles in.
  !
  ! Measured on the ITER example case, MPI=2 x OMP=16, 32-core dev box:
  ! sorting every step roughly DOUBLED average per-step wall time (~670ms
  ! -> ~1015-1818ms on non-collision steps). Worse, each individual sort
  ! call got ~3x more expensive too (97-180ms once per 10 steps at
  ! baseline vs 336-520ms every step here) - not just paid more often -
  ! for reasons not fully root-caused (mover_push itself looked flat to
  ! slightly better, ~268ms vs baseline's ~340ms, consistent with the
  ! locality theory being directionally right, but nowhere near enough to
  ! offset the added sort cost). Left in, disabled, so nobody re-tries
  ! this exact experiment without knowing the result.
  logical, parameter :: SORT_EVERY_STEP = .false.

  type :: Simulation
    type(State) :: state
    !type(CollisionWorkspace) :: coll_ws

    real(real64) :: t_Erho    = 0.0_real64
    real(real64) :: t_poisson = 0.0_real64
    real(real64) :: t_sort    = 0.0_real64
    real(real64) :: t_avg     = 0.0_real64
    real(real64) :: t_MC      = 0.0_real64
    real(real64) :: t_mover   = 0.0_real64
    real(real64) :: t_bck     = 0.0_real64
    real(real64) :: t_total   = 0.0_real64
    integer(int32) :: t_count = 0_int32

    ! Sub-timers breaking down mover and E/rho (deposit) into their actual
    ! component phases, added to localize the remaining mover/E-rho cost
    ! after class(ParticleSet)->type(ParticleSet) had no measurable effect.
    real(real64) :: t_mover_push = 0.0_real64
    real(real64) :: t_mover_bc   = 0.0_real64
    real(real64) :: t_dep_clear  = 0.0_real64
    real(real64) :: t_dep_loop   = 0.0_real64
    real(real64) :: t_dep_reduce = 0.0_real64

    ! Min/avg across iproc of t_mover_push/t_dep_loop (t_mover_push and
    ! t_dep_loop above already hold the max across iproc, since that's what
    ! bounds the fused OMP region's wall time). Only used to check OMP load
    ! imbalance across iproc - see advance_particles_local (mod_state.f90).
    real(real64) :: t_mover_push_min = 0.0_real64
    real(real64) :: t_mover_push_avg = 0.0_real64
    real(real64) :: t_dep_loop_min   = 0.0_real64
    real(real64) :: t_dep_loop_avg   = 0.0_real64

    ! TEMPORARY diagnostic - see move_and_bc_boris's header comment
    ! (mod_particleMover.f90): tallies how many live particles ran the full
    ! Boris rotation vs. hit the negligible-|B| electrostatic fallback,
    ! summed over the steps between diagnostic prints (same
    ! accumulate-then-reset pattern as t_mover_push etc above).
    ! int64: summed over ~all particles x nsav steps, which passes 2**31
    ! (e.g. ~4.5M electrons/rank x 1000 steps at ITER 8x24).
    integer(int64) :: n_boris_used  = 0_int64
    integer(int64) :: n_boris_total = 0_int64

    ! Toggles DATA.BAK/ <-> DATA.BAK2/ on each backup write, matching
    ! legacy's flag_wrt (so a crash mid-write never destroys the only
    ! valid backup - the other directory still has the previous one).
    integer(int32) :: flag_wrt = 0_int32

  contains
    procedure :: init
    procedure :: build_initial_fields
    procedure :: write_initial_diagnostics
    procedure :: deposit_all_particles
    procedure :: collisions_step
    procedure :: coulomb_step
    procedure :: output_step
    procedure :: reset_2d_averages
    procedure :: advance_one_step
    procedure :: print_debug_state
    procedure :: print_diagnostics
    procedure :: run
    procedure :: finalize
  end type Simulation

contains

  subroutine init(self, comm_in)
    class(Simulation), intent(inout) :: self
    integer, intent(in) :: comm_in
    integer :: rank, ierr

    ! Print the banner/license once, from rank 0 only, before anything else.
    call MPI_Comm_rank(comm_in, rank, ierr)
    if (rank == 0) call print_introduction()

    call self%state%init(comm_in)

    ! The energy-conserving push_scheme (mod_config.f90) supports plain wall
    ! boundaries, y/z-periodic (flag_pbc/flag_pbcz - see the ghost fix-up in
    ! calc_Efield_energy_conserving, mod_electricField.f90), and dielectric
    ! (flag_die) walls. Dielectric needs no special-casing in the E-field
    ! calc or particle mover: apply_dielectric_bc_to_phi (mod_state.f90)
    ! folds each dielectric node's accumulated-charge potential into phi
    ! BEFORE the Poisson solve/E-field calc run, every step, regardless of
    ! push_scheme - so by the time calc_Efield_energy_conserving differences
    ! phi across faces, a dielectric node's phi is just another already-
    ! correct nodal value, indistinguishable from a fixed Dirichlet wall's.
    ! The face-centered gather (gather_E_energy_conserving,
    ! mod_particleMover.f90) and the dielectric surface-charge deposit
    ! (move_and_bc_electrostatic's flag_die branch) are likewise shared,
    ! unbranched code for both push schemes. Neumann (flag_nmn) is the one
    ! still-unimplemented case: calc_Efield_modular gives bcnd==-2 nodes a
    ! dedicated one-sided ghost formula that calc_Efield_energy_conserving
    ! has no counterpart for. Fail fast at startup rather than silently
    ! produce wrong physics for that untested combination.
    if (trim(self%state%cfg%push_scheme) == 'energy') then
      if (self%state%dom%flag_nmn == 1) then
        if (self%state%mpi_rank == 0) then
          write(*,'(a)') "ERROR: cfg%push_scheme = 'energy' does not support " // &
            "Neumann boundaries (flag_nmn=0 required). Plain walls, " // &
            "periodic (flag_pbc/flag_pbcz), and dielectric (flag_die) are supported."
        end if
        error stop 1
      end if
    end if
  end subroutine init


  subroutine build_initial_fields(self)
    class(Simulation), intent(inout) :: self

    call build_rho_from_np( &
         n        = int(self%state%dom%n, int32), &
         np_red   = self%state%fld%np, &
         charge   = self%state%chem%charge(1:self%state%ntype), &
         ntype    = int(self%state%ntype), &
         rho      = self%state%fld%rho, &
         bcnd     = self%state%dom%bcnd, &
         flag_pbc = int(self%state%dom%flag_pbc, int32))

    call self%state%apply_dielectric_bc_to_phi()

    call solve_poisson_legacy( &
         pdec        = self%state%pdec, &
         phi_global  = self%state%fld%phi, &
         bcnd_global = self%state%dom%bcnd, &
         rhs_global  = self%state%fld%rho(0:self%state%dom%n(1)+1, &
                                          0:self%state%dom%n(2)+1, &
                                          0:self%state%dom%n(3)+1), &
         h           = self%state%dom%h, &
         n_in        = int(self%state%dom%n, int32), &
         ncycl       = 20000, &
         eps         = self%state%cfg%eps, &
         omega       = self%state%cfg%omega, &
         ng          = self%state%cfg%ng, &
         flag_pbc_in = self%state%dom%flag_pbc, &
         flag_nmn_in = self%state%dom%flag_nmn )

    if (trim(self%state%cfg%push_scheme) == 'energy') then
      call calc_Efield_energy_conserving( &
           n         = int(self%state%dom%n, int32), &
           h         = self%state%dom%h, &
           phi       = self%state%fld%phi, &
           E         = self%state%fld%E, &
           bcnd      = self%state%dom%bcnd, &
           flag_pbc  = int(self%state%dom%flag_pbc, int32), &
           flag_pbcz = int(self%state%dom%flag_pbcz, int32) )
    else
      call calc_Efield_modular( &
           n    = int(self%state%dom%n, int32), &
           h    = self%state%dom%h, &
           phi  = self%state%fld%phi, &
           E    = self%state%fld%E, &
           bcnd = self%state%dom%bcnd )
    end if
  end subroutine build_initial_fields


  subroutine write_initial_diagnostics(self)
    class(Simulation), intent(in) :: self

    character(len=128) :: prefix
    character(len=10)  :: s
    integer(int32) :: i
    integer(int32) :: iy_plane_phi, iy_plane_density, iy_plane_E

    if (self%state%mpi_rank /= 0) return

    iy_plane_density = int(self%state%dom%n(2)/2 + 1, int32)
    iy_plane_phi     = int(self%state%dom%n(2)/2,     int32)
    iy_plane_E       = int(self%state%dom%n(2)/2 + 1, int32)

    do i = 1, self%state%ntype
      write(s,'(i0)') i
      prefix = './Output/Output_2D/n' // trim(s)

      call write_density_planes( &
        np       = self%state%fld%np, &
        n        = int(self%state%dom%n, int32), &
        ptype    = i, &
        ix_plane = self%state%params%ix_plot_plane, &
        iy_plane = iy_plane_density, &
        iz_plane = self%state%params%iz_plot_plane, &
        every    = 1_int32, &
        prefix   = prefix )
    end do

    call write_scalar_planes( &
      f        = self%state%fld%phi, &
      n        = int(self%state%dom%n, int32), &
      ix_plane = self%state%params%ix_plot_plane, &
      iy_plane = iy_plane_phi, &
      iz_plane = self%state%params%iz_plot_plane, &
      every    = 1_int32, &
      prefix   = './Output/Output_2D/phi1' )

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      1_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, './Output/Output_2D/Ex')

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      2_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, './Output/Output_2D/Ey')

    call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                      3_int32, self%state%params%ix_plot_plane, &
                                      iy_plane_E, self%state%params%iz_plot_plane, &
                                      1_int32, './Output/Output_2D/Ez')
  end subroutine write_initial_diagnostics


  subroutine deposit_all_particles(self)
    class(Simulation), intent(inout) :: self

    real(real64)   :: dt0, dt1

    ! The per-thread clear+deposit loop that used to live here now runs
    ! inside state%advance_particles_local(), fused into the same OMP
    ! region as the push and BC/SEE pass (see mod_state.f90) instead of
    ! its own separate fork - one less parallel region per step. t_dep_clear
    ! and t_dep_loop are kept at zero (folded into t_mover_push, where that
    ! fused region is now timed) rather than removed, so old log-parsing
    ! scripts expecting those fields don't break.
    self%t_dep_clear = self%t_dep_clear
    self%t_dep_loop  = self%t_dep_loop

    ! Builds next step's rho here, right after the deposit (nothing else
    ! writes rho or np in between), and MPI-sums only rho. fld%np is left
    ! rank-local until a reader calls state%sync_np() - see
    ! reduce_density_and_rho.
    dt0 = MPI_Wtime()
    call reduce_density_and_rho( &
      n         = int(self%state%dom%n, int32), &
      bcnd      = self%state%dom%bcnd, &
      np_thread = self%state%np_thread, &
      ntype     = int(self%state%ntype), &
      nproc     = int(self%state%nproc), &
      mpi_comm  = self%state%comm, &
      charge    = self%state%chem%charge(1:self%state%ntype), &
      np_red    = self%state%fld%np, &
      rho       = self%state%fld%rho, &
      zlo       = self%state%np_zlo, &
      zhi       = self%state%np_zhi, &
      zend      = self%state%np_zend )
    self%state%np_synced = .false.
    dt1 = MPI_Wtime()
    self%t_dep_reduce = self%t_dep_reduce + (dt1 - dt0)
  end subroutine deposit_all_particles


  subroutine collisions_step(self)
    class(Simulation), intent(inout) :: self

    call self%state%sync_np()

    call perform_collisions_step( &
        part          = self%state%part, &
        n             = self%state%dom%n, &
        h             = self%state%dom%h, &
        ntype_tracked = self%state%ntype, &
        ntype_all     = self%state%rxn%ntype, &
        mass          = self%state%chem%mass(1:self%state%rxn%ntype), &
        charge        = self%state%chem%charge(1:self%state%rxn%ntype), &
        Ti            = self%state%chem%Ti(1:self%state%rxn%ntype), &
        Nm            = self%state%params%Nm(1:self%state%rxn%ntype), &
        p_ncol        = self%state%chem%p_ncol(1:self%state%rxn%ntype), &
        sig_list      = self%state%rxn%sig_list, &
        col_info      = self%state%rxn%col_info, &
        sigv_mx       = self%state%rxn%sigv_mx, &
        sig           = self%state%rxn%sig, &
        sig_Er        = self%state%rxn%sig_Er, &
        sig_Eex       = self%state%rxn%sig_Eex, &
        ni0           = self%state%chem%ni0(1:self%state%rxn%ntype), &
        ns_coll       = self%state%params%nb_step_collisions, &
        dt            = self%state%params%dt, &
        nu_uplim      = self%state%params%nu_uplim(1:self%state%rxn%ntype), &
        iseed         = self%state%params%iseed, &
        mpi_rank      = int(self%state%mpi_rank, int32), &
        Pcoll         = self%state%P_loss(3,:,:), &
        mom_loss      = self%state%mom_loss(:,3,:,:), &
        dom_volume    = self%state%dom%h(1) * self%state%dom%h(2) * self%state%dom%h(3) * &
                        real(self%state%dom%n(1)*self%state%dom%n(2)*self%state%dom%n(3), real64), &
        np_red        = self%state%fld%np, &
        bcnd          = self%state%dom%bcnd, &
        ix_plane      = self%state%params%ix_plot_plane, &
        iy_plane      = int(self%state%dom%n(2)/2 + 1, int32), &
        iz_plane      = self%state%params%iz_plot_plane, &
        sour_xy       = self%state%sour_avg_xy, &
        sour_xz       = self%state%sour_avg_xz, &
        sour_yz       = self%state%sour_avg_yz, &
        sink_xy       = self%state%sink_avg_xy, &
        sink_xz       = self%state%sink_avg_xz, &
        sink_yz       = self%state%sink_avg_yz)

  end subroutine collisions_step


  subroutine coulomb_step(self)
    class(Simulation), intent(inout) :: self

    call self%state%sync_np()

    call perform_coulomb_step( &
        part   = self%state%part, &
        ntype  = self%state%ntype, &
        nproc  = self%state%nproc, &
        mass   = self%state%chem%mass(1:self%state%ntype), &
        charge = self%state%chem%charge(1:self%state%ntype), &
        Nm     = self%state%params%Nm(1:self%state%ntype), &
        n_e    = average_species_density( &
                   int(self%state%dom%n, int32), self%state%fld%np, &
                   int(self%state%ntype, int32), 1_int32, &
                   int(self%state%dom%flag_pbc, int32), int(self%state%dom%flag_pbcz, int32)), &
        dt     = self%state%params%dt * real(self%state%params%nb_step_coulomb, real64), &
        iseed  = self%state%params%iseed)

  end subroutine coulomb_step


  subroutine reset_2d_averages(self)
    class(Simulation), intent(inout) :: self

    if (allocated(self%state%np_avg_xy))  self%state%np_avg_xy  = 0.0_real64
    if (allocated(self%state%np_avg_xz))  self%state%np_avg_xz  = 0.0_real64
    if (allocated(self%state%np_avg_yz))  self%state%np_avg_yz  = 0.0_real64

    if (allocated(self%state%phi_avg_xy)) self%state%phi_avg_xy = 0.0_real64
    if (allocated(self%state%phi_avg_xz)) self%state%phi_avg_xz = 0.0_real64
    if (allocated(self%state%phi_avg_yz)) self%state%phi_avg_yz = 0.0_real64

    if (allocated(self%state%data_pavg_xy)) self%state%data_pavg_xy = 0.0_real64
    if (allocated(self%state%data_pavg_xz)) self%state%data_pavg_xz = 0.0_real64
    if (allocated(self%state%data_pavg_yz)) self%state%data_pavg_yz = 0.0_real64

    if (allocated(self%state%sour_avg_xy)) self%state%sour_avg_xy = 0.0_real64
    if (allocated(self%state%sour_avg_xz)) self%state%sour_avg_xz = 0.0_real64
    if (allocated(self%state%sour_avg_yz)) self%state%sour_avg_yz = 0.0_real64
    if (allocated(self%state%sink_avg_xy)) self%state%sink_avg_xy = 0.0_real64
    if (allocated(self%state%sink_avg_xz)) self%state%sink_avg_xz = 0.0_real64
    if (allocated(self%state%sink_avg_yz)) self%state%sink_avg_yz = 0.0_real64

    if (allocated(self%state%avg3D)) self%state%avg3D = 0.0_real64

    self%state%cnt_avg = 0_int32
  end subroutine reset_2d_averages


  subroutine output_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32), intent(in) :: istep

    integer(int32) :: i
    integer(int32) :: ix_plane, iy_plane_E, iz_plane
    character(len=256) :: prefix
    character(len=16)  :: sstep, sspecies
    real(real64) :: avg_factor, dr_cell, cell_volume, dt_interval
    real(real64), allocatable :: tmp_xy(:,:), tmp_xz(:,:), tmp_yz(:,:)
    real(real64), allocatable :: ang_xy(:,:), ang_xz(:,:), ang_yz(:,:)

    if (self%state%mpi_rank /= 0) return

    ix_plane   = self%state%params%ix_plot_plane
    iz_plane   = self%state%params%iz_plot_plane
    iy_plane_E = int(self%state%dom%n(2)/2 + 1, int32)

    avg_factor  = real(max(1_int32, self%state%cnt_avg), real64)
    cell_volume = self%state%dom%h(1) * self%state%dom%h(2) * self%state%dom%h(3)
    dt_interval = real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt

    allocate(tmp_xy(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2))
    allocate(tmp_xz(0:self%state%dom%n(1)+2,0:self%state%dom%n(3)+2))
    allocate(tmp_yz(0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))
    allocate(ang_xy(0:self%state%dom%n(1)+2,0:self%state%dom%n(2)+2))
    allocate(ang_xz(0:self%state%dom%n(1)+2,0:self%state%dom%n(3)+2))
    allocate(ang_yz(0:self%state%dom%n(2)+2,0:self%state%dom%n(3)+2))

    write(sstep,'(i0)') istep

    ! --- Debye length ratio (legacy dr_*.mco): dr_cell/lambda_De, from
    ! species 1 (electron) T/n moments - dividing by data_pavg_xy(1,...)
    ! rather than np_avg_xy since legacy computes this from the SAME
    ! per-particle-moment n/T pair used for the u/j outputs below, not the
    ! grid-deposited density used for the n<i> files above. Guarded exactly
    ! like legacy (only take the ratio where lambda_De > 0) to avoid a
    ! divide-by-zero where a cell has no sampled electrons this step.
    dr_cell = (self%state%dom%h(1) * self%state%dom%h(2) * self%state%dom%h(3)) ** (1.0_real64/3.0_real64)

    tmp_xy = 0.0_real64
    where (self%state%data_pavg_xy(1,:,:,1) > 0.0_real64)
      tmp_xy = sqrt( eps0 * self%state%data_pavg_xy(2,:,:,1) / (self%state%data_pavg_xy(1,:,:,1) * qe) )
    end where
    where (tmp_xy > 0.0_real64) tmp_xy = dr_cell / tmp_xy

    tmp_xz = 0.0_real64
    where (self%state%data_pavg_xz(1,:,:,1) > 0.0_real64)
      tmp_xz = sqrt( eps0 * self%state%data_pavg_xz(2,:,:,1) / (self%state%data_pavg_xz(1,:,:,1) * qe) )
    end where
    where (tmp_xz > 0.0_real64) tmp_xz = dr_cell / tmp_xz

    tmp_yz = 0.0_real64
    where (self%state%data_pavg_yz(1,:,:,1) > 0.0_real64)
      tmp_yz = sqrt( eps0 * self%state%data_pavg_yz(2,:,:,1) / (self%state%data_pavg_yz(1,:,:,1) * qe) )
    end where
    where (tmp_yz > 0.0_real64) tmp_yz = dr_cell / tmp_yz

    prefix = './Output/Output_2D/it' // trim(sstep) // '_dr'

    call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
    call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
    call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

    do i = 1, self%state%ntype
      write(sspecies,'(i0)') i

      tmp_xy = self%state%np_avg_xy(:,:,i) / avg_factor
      tmp_xz = self%state%np_avg_xz(:,:,i) / avg_factor
      tmp_yz = self%state%np_avg_yz(:,:,i) / avg_factor

      prefix = './Output/Output_2D/it' // trim(sstep) // '_n' // trim(sspecies)

      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

      tmp_xy = self%state%data_pavg_xy(2,:,:,i) / avg_factor
      tmp_xz = self%state%data_pavg_xz(2,:,:,i) / avg_factor
      tmp_yz = self%state%data_pavg_yz(2,:,:,i) / avg_factor

      prefix = './Output/Output_2D/it' // trim(sstep) // '_T' // trim(sspecies)

      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

      ! --- raw mean velocity components (legacy u<i>x/y/z_*.mco). Not
      ! divided by avg_factor: unlike np_avg_xy/data_pavg(Tp_avg), these
      ! (data_pavg indices 3-5) are overwritten - not summed - by each
      ! compute_particle_plane_moments_species call, so they already hold
      ! a single normalized snapshot (see mod_planeMoments.f90).
      prefix = './Output/Output_2D/it' // trim(sstep) // '_u' // trim(sspecies) // 'x'
      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', self%state%data_pavg_xy(3,:,:,i), int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', self%state%data_pavg_xz(3,:,:,i), int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', self%state%data_pavg_yz(3,:,:,i), int(self%state%dom%n, int32), 1_int32)

      prefix = './Output/Output_2D/it' // trim(sstep) // '_u' // trim(sspecies) // 'y'
      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', self%state%data_pavg_xy(4,:,:,i), int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', self%state%data_pavg_xz(4,:,:,i), int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', self%state%data_pavg_yz(4,:,:,i), int(self%state%dom%n, int32), 1_int32)

      prefix = './Output/Output_2D/it' // trim(sstep) // '_u' // trim(sspecies) // 'z'
      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', self%state%data_pavg_xy(5,:,:,i), int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', self%state%data_pavg_xz(5,:,:,i), int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', self%state%data_pavg_yz(5,:,:,i), int(self%state%dom%n, int32), 1_int32)

      ! --- flux magnitude + direction (legacy j<i>_*.mco): n*|u| plus a
      ! "vector" block of the in-plane flow angle. Uses the per-particle
      ! moment density (index 1), the same source legacy's own j uses -
      ! not np_avg_xy, which is the separate grid-deposited density used
      ! for n<i>_*.mco above.
      tmp_xy = self%state%data_pavg_xy(1,:,:,i) * sqrt( self%state%data_pavg_xy(3,:,:,i)**2 &
                                                        + self%state%data_pavg_xy(4,:,:,i)**2 &
                                                        + self%state%data_pavg_xy(5,:,:,i)**2 )
      ang_xy = atan2( self%state%data_pavg_xy(4,:,:,i), self%state%data_pavg_xy(3,:,:,i) )

      tmp_xz = self%state%data_pavg_xz(1,:,:,i) * sqrt( self%state%data_pavg_xz(3,:,:,i)**2 &
                                                        + self%state%data_pavg_xz(4,:,:,i)**2 &
                                                        + self%state%data_pavg_xz(5,:,:,i)**2 )
      ang_xz = atan2( self%state%data_pavg_xz(5,:,:,i), self%state%data_pavg_xz(3,:,:,i) )

      tmp_yz = self%state%data_pavg_yz(1,:,:,i) * sqrt( self%state%data_pavg_yz(3,:,:,i)**2 &
                                                        + self%state%data_pavg_yz(4,:,:,i)**2 &
                                                        + self%state%data_pavg_yz(5,:,:,i)**2 )
      ang_yz = atan2( self%state%data_pavg_yz(5,:,:,i), self%state%data_pavg_yz(4,:,:,i) )

      prefix = './Output/Output_2D/it' // trim(sstep) // '_j' // trim(sspecies)
      call write_plane_xy_flux_2d(trim(prefix)//'_xy.mco', tmp_xy, ang_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_flux_2d(trim(prefix)//'_xz.mco', tmp_xz, ang_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_flux_2d(trim(prefix)//'_yz.mco', tmp_yz, ang_yz, int(self%state%dom%n, int32), 1_int32)

      ! --- reaction-driven source/sink rate maps (legacy sour<i>/sink<i>_
      ! *.mco, plt_src==1 branch): sum the per-iproc accumulators from
      ! mod_collisionsGwenael's deposit_plane_event, then convert the
      ! interval's accumulated particle count into a rate density
      ! (particles / m^3 / s), same (nsav-1)*dt convention as Iw/Pw above.
      tmp_xy = sum(self%state%sour_avg_xy(:,:,i,:), dim=3) / (cell_volume * dt_interval)
      tmp_xz = sum(self%state%sour_avg_xz(:,:,i,:), dim=3) / (cell_volume * dt_interval)
      tmp_yz = sum(self%state%sour_avg_yz(:,:,i,:), dim=3) / (cell_volume * dt_interval)

      prefix = './Output/Output_2D/it' // trim(sstep) // '_sour' // trim(sspecies)
      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

      tmp_xy = sum(self%state%sink_avg_xy(:,:,i,:), dim=3) / (cell_volume * dt_interval)
      tmp_xz = sum(self%state%sink_avg_xz(:,:,i,:), dim=3) / (cell_volume * dt_interval)
      tmp_yz = sum(self%state%sink_avg_yz(:,:,i,:), dim=3) / (cell_volume * dt_interval)

      prefix = './Output/Output_2D/it' // trim(sstep) // '_sink' // trim(sspecies)
      call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
      call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
      call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)
    end do

    tmp_xy = self%state%phi_avg_xy / avg_factor
    tmp_xz = self%state%phi_avg_xz / avg_factor
    tmp_yz = self%state%phi_avg_yz / avg_factor

    prefix = './Output/Output_2D/it' // trim(sstep) // '_phi'

    call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
    call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
    call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)

    if (allocated(self%state%avg3D)) then
      ! flag_avg3D==1: time-averaged E planes and the averaged 3D phi /
      ! species-2 density map, as legacy write_data.f90 and main.f90 do.
      block
        character(len=2), parameter :: ecomp(3) = ['Ex', 'Ey', 'Ez']
        integer(int32) :: comp, n1, n2, n3
        integer :: u3d
        real(real64), allocatable :: map3D(:,:,:,:)

        n1 = int(self%state%dom%n(1), int32)
        n2 = int(self%state%dom%n(2), int32)
        n3 = int(self%state%dom%n(3), int32)

        do comp = 1, 3
          tmp_xy = 0.0_real64
          tmp_xz = 0.0_real64
          tmp_yz = 0.0_real64
          tmp_xy(1:n1+1,1:n2+1) = self%state%avg3D(comp,:,:,iz_plane)   / avg_factor
          tmp_xz(1:n1+1,1:n3+1) = self%state%avg3D(comp,:,iy_plane_E,:) / avg_factor
          tmp_yz(1:n2+1,1:n3+1) = self%state%avg3D(comp,ix_plane,:,:)   / avg_factor

          prefix = './Output/Output_2D/it' // trim(sstep) // '_' // ecomp(comp)
          call write_plane_xy_scalar_2d(trim(prefix)//'_xy.mco', tmp_xy, int(self%state%dom%n, int32), 1_int32)
          call write_plane_xz_scalar_2d(trim(prefix)//'_xz.mco', tmp_xz, int(self%state%dom%n, int32), 1_int32)
          call write_plane_yz_scalar_2d(trim(prefix)//'_yz.mco', tmp_yz, int(self%state%dom%n, int32), 1_int32)
        end do

        ! Same record layout as legacy DATA/phi_n_3D.dat, so
        ! post_analysis/phi_n_3D/convert.f90 reads it unchanged.
        write(*,*) 'Saving 3D potentiel and density maps ...'
        open(newunit=u3d, file='./Output/phi_n_3D.dat', form='unformatted', &
             status='replace', action='write')
        write(u3d) self%state%dom%n(1), self%state%dom%n(2), self%state%dom%n(3)
        ! Explicit heap buffer: the inline expression makes the compiler
        ! build a stack temporary (~39 MB at ITER size), which segfaults
        ! under the default 8 MB stack limit.
        allocate(map3D(2, n1+1, n2+1, n3+1))
        map3D = self%state%avg3D(4:5,:,:,:) / avg_factor
        write(u3d) map3D
        deallocate(map3D)
        close(u3d)
      end block
    else
      prefix = './Output/Output_2D/it' // trim(sstep) // '_Ex'
      call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                        1_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

      prefix = './Output/Output_2D/it' // trim(sstep) // '_Ey'
      call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                        2_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)

      prefix = './Output/Output_2D/it' // trim(sstep) // '_Ez'
      call write_vector_component_planes(self%state%fld%E, int(self%state%dom%n, int32), &
                                        3_int32, ix_plane, iy_plane_E, iz_plane, 1_int32, prefix)
    end if

    ! --- legacy phi_Te_ne_cntr.dat: domain-center phi/Te/ne (species 1 =
    ! electrons) plus ne sampled along x. Same data_pavg source as the
    ! dr/u/j outputs above (not np_avg_xy), matching legacy's own internal
    ! consistency. Written here (not print_diagnostics) because phi_avg_xy/
    ! data_pavg_xy are only freshly populated for THIS interval at this
    ! point - print_diagnostics(istep) runs before advance_one_step(istep),
    ! i.e. before this interval's compute_plane_moments_local/
    ! accumulate_2d_averages, so it would see last interval's already-reset
    ! (zero) values instead.
    block
      integer :: ucntr, k
      integer(int32) :: nx, iy_c
      real(real64) :: simulation_time
      logical :: file_exists

      nx  = int(self%state%dom%n(1), int32)
      iy_c = int(self%state%dom%n(2)/2_int32 + 1_int32, int32)
      simulation_time = self%state%params%dt * real(istep, real64)

      inquire(file='./Output/phi_Te_ne_cntr.dat', exist=file_exists)
      open(newunit=ucntr, file='./Output/phi_Te_ne_cntr.dat', status='unknown', &
           position='append', action='write')
      if (.not. file_exists) then
        write(ucntr,'(a)') '# Time(s), phi(xm/2,ym/2,z_pl), Te(xm/2,ym/2,z_pl), ne at ym/2, z_pl and '// &
            'xm/2, xm/10, 2xm/10, 3xm/10, 4xm/10, 6xm/10, 7xm/10, 8xm/10, 9xm/10'
      end if
      write(ucntr,'(20(1x,es16.8))') simulation_time, &
          self%state%phi_avg_xy(nx/2_int32+1_int32, iy_c) / avg_factor, &
          self%state%data_pavg_xy(2, nx/2_int32+1_int32, iy_c, 1) / avg_factor, &
          self%state%data_pavg_xy(1, nx/2_int32+1_int32, iy_c, 1) / avg_factor, &
          ( self%state%data_pavg_xy(1, nint(real(k,real64)*real(nx,real64)/10.0_real64)+1_int32, iy_c, 1) &
              / avg_factor, k = 1_int32, 4_int32 ), &
          ( self%state%data_pavg_xy(1, nint(real(k,real64)*real(nx,real64)/10.0_real64)+1_int32, iy_c, 1) &
              / avg_factor, k = 6_int32, 9_int32 )
      close(ucntr)
    end block

    if (DEBUG_SIMULATION) then
      write(*,*) "MOD output_step: istep, cnt_avg = ", istep, self%state%cnt_avg
    end if

    deallocate(tmp_xy, tmp_xz, tmp_yz)
    deallocate(ang_xy, ang_xz, ang_yz)
  end subroutine output_step


  subroutine print_debug_state(self, istep, label)
    class(Simulation), intent(in) :: self
    integer(int32), intent(in) :: istep
    character(len=*), intent(in) :: label

    if (.not. DEBUG_SIMULATION) return
    if (self%state%mpi_rank /= 0) return

    write(*,*) "MOD DEBUG ", trim(label), " istep = ", istep
    write(*,*) "rho min/max/sum = ", minval(self%state%fld%rho), &
                                  maxval(self%state%fld%rho), &
                                  sum(self%state%fld%rho)
    write(*,*) "phi min/max/sum = ", minval(self%state%fld%phi), &
                                  maxval(self%state%fld%phi), &
                                  sum(self%state%fld%phi)
  end subroutine print_debug_state


  subroutine advance_one_step(self, istep)
    class(Simulation), intent(inout) :: self
    integer(int32), intent(in) :: istep

    real(real64)   :: vt_heat
    real(real64)   :: t0, t1, tstep0
    real(real64)   :: t_move_dbg, t_bc_dbg, t_deposit_dbg
    real(real64)   :: t_move_min_dbg, t_move_avg_dbg, t_deposit_min_dbg, t_deposit_avg_dbg
    integer(int32) :: n_boris_used_dbg, n_boris_total_dbg
    integer(int32) :: iproc_h

    tstep0 = MPI_Wtime()

    ! RF antenna heating: periodic convergence (skin depth + E0_RF rescale),
    ! ported from legacy Src/main.f90:887-933. Placed here, before
    ! build_rho_from_np below, so it reads self%state%fld%np as left by the
    ! PREVIOUS step's deposit_all_particles - mirrors legacy, whose
    ! equivalent block (main.f90:859) also runs before that iteration's
    ! calc_rho (main.f90:970), reading the same "previous period" density
    ! and already-accumulated P_RF.
    ! Cadence is a count of accumulated steps, not mod(istep,ns_RF): after a
    ! restart istep starts at the backup step, which is generally not a
    ! multiple of ns_RF, and a partial first window would under-read P_RF.
    if (self%state%cfg%flag_RFant == 1_int32) then
      if (self%state%nstep_RF >= self%state%params%ns_RF) then
        call self%state%update_rf_convergence(istep)
      end if
      ! This step's mover call (below) accumulates into P_RF.
      self%state%nstep_RF = self%state%nstep_RF + 1_int32
    end if

    ! E/rho
    !
    ! rho is already built: by the PREVIOUS step's deposit_all_particles
    ! (reduce_density_and_rho), or for istep==1 by build_initial_fields.
    ! Nothing writes rho in between, so rebuilding it here from fld%np
    ! (as this used to) would only recompute it - and would need the full
    ! per-species MPI sync that the per-step path now skips.

    ! Poisson + E field
    t0 = MPI_Wtime()
    call self%state%apply_dielectric_bc_to_phi()

    call solve_poisson_legacy( &
      pdec        = self%state%pdec, &
      phi_global  = self%state%fld%phi, &
      bcnd_global = self%state%dom%bcnd, &
      rhs_global  = self%state%fld%rho(0:self%state%dom%n(1)+1, &
                                      0:self%state%dom%n(2)+1, &
                                      0:self%state%dom%n(3)+1), &
      h           = self%state%dom%h, &
      n_in        = int(self%state%dom%n, int32), &
      ncycl       = 20000, &
      eps         = self%state%cfg%eps, &
      omega       = self%state%cfg%omega, &
      ng          = self%state%cfg%ng, &
      flag_pbc_in = self%state%dom%flag_pbc, &
      flag_nmn_in = self%state%dom%flag_nmn )

    if (trim(self%state%cfg%push_scheme) == 'energy') then
      call calc_Efield_energy_conserving( &
        n         = int(self%state%dom%n, int32), &
        h         = self%state%dom%h, &
        phi       = self%state%fld%phi, &
        E         = self%state%fld%E, &
        bcnd      = self%state%dom%bcnd, &
        flag_pbc  = int(self%state%dom%flag_pbc, int32), &
        flag_pbcz = int(self%state%dom%flag_pbcz, int32) )
    else
      call calc_Efield_modular( &
        n    = int(self%state%dom%n, int32), &
        h    = self%state%dom%h, &
        phi  = self%state%fld%phi, &
        E    = self%state%fld%E, &
        bcnd = self%state%dom%bcnd )
    end if
    t1 = MPI_Wtime()
    self%t_poisson = self%t_poisson + (t1 - t0)

    ! Averaging
    t0 = MPI_Wtime()
    if (self%state%params%nb_step_averaging > 0_int32) then
      if (mod(istep, self%state%params%nb_step_averaging) == 1_int32) then
        call self%state%compute_plane_moments_local()
        call self%state%accumulate_2d_averages()
      end if
    end if
    t1 = MPI_Wtime()
    self%t_avg = self%t_avg + (t1 - t0)

    ! Reset dielectric charge counters
    if (allocated(self%state%sum_q_xz)) self%state%sum_q_xz = 0.0_real64
    if (allocated(self%state%sum_q_yz)) self%state%sum_q_yz = 0.0_real64

    ! Electron heating — fires BEFORE the push, matching legacy's sequence:
    !   [outside OMP] read Nh/sum_dEk from previous deposit → compute vt
    !   [OMP]         reset Nh/sum_dEk; call eheating(vt); push; deposit
    ! In modular the deposit/accumulation stays after BC (below), but the
    ! vt computation and heating application happen here, before the push.
    if (self%state%cfg%Pabs > 0.0_real64 .and. self%state%cfg%flag_heat == 1) then
      if (self%state%params%nb_step_heating > 0_int32) then
        if (mod(istep, self%state%params%nb_step_heating) == 0_int32) then
          call self%state%update_heating_vt(vt_heat)
          if (self%state%cfg%flag_inj == 1) vt_heat = self%state%params%vt0(1)
          if (allocated(self%state%Nh))      self%state%Nh      = 0_int32
          if (allocated(self%state%sum_dEk)) self%state%sum_dEk = 0.0_real64
          call self%state%apply_electron_heating_local(vt_heat)
        end if
      end if
    end if

    ! EXPERIMENTAL (SORT_EVERY_STEP): see flag declaration above. Sorts
    ! immediately before the fused push+deposit region below, on whatever
    ! order the previous step (injection/BC) left particles in, so that
    ! region's field gather (mover) and charge scatter (deposit) both run
    ! on cell-sorted particles instead of arbitrary order. Timed into the
    ! same t_sort accumulator as the collision-triggered sort below, so
    ! this shows up in the existing "sorting" line of the diagnostic
    ! print without adding a new counter.
    if (SORT_EVERY_STEP) then
      t0 = MPI_Wtime()
      call self%state%sort_particles_local()
      t1 = MPI_Wtime()
      self%t_sort = self%t_sort + (t1 - t0)
    end if

    ! Mover + particle BC/SEE + charge deposition, fused into a single OMP
    ! region (see state%advance_particles_local in mod_state.f90) instead
    ! of three separate parallel regions. advance_particles_local returns
    ! its own internal (per-thread, max-across-iproc) split of that region
    ! into push/BC/deposit time, since the outer MPI_Wtime() bracket here
    ! can no longer separate them now that they share one fork. t_mover_bc
    ! and t_dep_loop - previously kept at zero after the fusion - now carry
    ! that real BC time and deposit-loop time respectively, restoring the
    ! "mover / deposit timing breakdown" print block below to real numbers.
    t0 = MPI_Wtime()
    call advance_particles_local(self%state, istep, t_move_dbg, t_bc_dbg, t_deposit_dbg, &
                                  t_move_min_dbg, t_move_avg_dbg, &
                                  t_deposit_min_dbg, t_deposit_avg_dbg, &
                                  n_boris_used_dbg, n_boris_total_dbg)
    t1 = MPI_Wtime()
    self%n_boris_used  = self%n_boris_used  + int(n_boris_used_dbg, int64)
    self%n_boris_total = self%n_boris_total + int(n_boris_total_dbg, int64)
    self%t_mover_push = self%t_mover_push + t_move_dbg
    self%t_mover_bc   = self%t_mover_bc   + t_bc_dbg
    self%t_dep_loop   = self%t_dep_loop   + t_deposit_dbg

    self%t_mover_push_min = self%t_mover_push_min + t_move_min_dbg
    self%t_mover_push_avg = self%t_mover_push_avg + t_move_avg_dbg
    self%t_dep_loop_min   = self%t_dep_loop_min   + t_deposit_min_dbg
    self%t_dep_loop_avg   = self%t_dep_loop_avg   + t_deposit_avg_dbg

    ! t_mover kept as the combined push+BC+deposit total for the existing
    ! top-level summary line (comparable to legacy's fused mover timer,
    ! which also spans push+BC/SEE+deposit - see Src/main.f90 ctime(6)).
    self%t_mover = self%t_mover_push + self%t_mover_bc + self%t_dep_loop

    ! Nh/sum_dEk (electron heating-region tally) are now accumulated
    ! inside advance_particles_local's ptype==1 deposit, above - one
    ! less full pass over the electrons per step (see HeatRegionTally in
    ! mod_chargeDeposition.f90). It overwrites Nh(iproc)/sum_dEk(iproc)
    ! every step, so no reset is needed here any more.

    ! Reduce this step's per-thread deposit (already computed inside
    ! advance_particles_local, above) into np_red BEFORE collisions (legacy
    ! ordering: calc_rho is called before collisions in Src/main.f90, so
    ! the density collisions read - np_red/np_mx - is always synchronized
    ! with the SAME particle positions collisions are about to act on,
    ! never a step stale. This single reduction still also feeds next
    ! iteration's E/rho/Poisson, same as before - nothing else deposits
    ! again before then.)
    t0 = MPI_Wtime()
    call self%deposit_all_particles()
    t1 = MPI_Wtime()
    self%t_Erho = self%t_Erho + (t1 - t0)

    ! MC collisions (every nb_step_collisions steps)
    t0 = MPI_Wtime()
    ! Fires at the END of steps nb, 2*nb, ... (mod == 0). That is the same
    ! point, on the same density, as legacy's collisions at the START of
    ! steps nb+1, 2*nb+1, ... (Src/main.f90: MOD(it,ns_coll).eq.1, before
    ! that step's push). Using mod == 1 here ran them one step later than
    ! legacy, so a heating step (mod(istep,4)==0) could directly follow a
    ! collision step and build vt from a heating-region tally taken before
    ! the collisions changed the electron population: those heating steps
    ! delivered ~10% more than Pabs (+2..5% on average, and a denser
    ! plasma than legacy).
    if (self%state%params%nb_step_collisions > 0_int32) then
      if (mod(istep, self%state%params%nb_step_collisions) == 0_int32) then

        ! Legacy-like: build Plist/cell lists from current post-mover particles
        call self%state%sort_particles_local()

        t1 = MPI_Wtime()
        self%t_sort = self%t_sort + (t1 - t0)

        t0 = MPI_Wtime()
        call self%collisions_step()

      end if
    end if
    t1 = MPI_Wtime()
    self%t_MC = self%t_MC + (t1 - t0)

    ! Coulomb collisions (every nb_step_coulomb steps, independent of MC)
    if (self%state%cfg%flag_coulomb == 1_int32) then
      if (self%state%params%nb_step_coulomb > 0_int32) then
        if (mod(istep, self%state%params%nb_step_coulomb) == 1_int32) then
          call self%coulomb_step()
        end if
      end if
    end if

    ! Volume particle injection (legacy: part_injection, inside OMP per iproc)
    ! Fires every step when flag_inj==1 (ns_inj=1 hardcoded, same as legacy).
    if (self%state%cfg%flag_inj == 1_int32) then
      !$omp parallel do private(iproc_h) schedule(static)
      do iproc_h = 1, self%state%nproc
        call inject_particles_volume( &
          part         = self%state%part, &
          iproc        = iproc_h, &
          cfg          = self%state%cfg, &
          n            = int(self%state%dom%n, int32), &
          h            = self%state%dom%h, &
          bcnd         = self%state%dom%bcnd, &
          phi          = self%state%fld%phi, &
          vt0          = self%state%params%vt0, &
          Nm           = self%state%params%Nm, &
          mass         = self%state%chem%mass, &
          ni0          = self%state%chem%ni0, &
          charge       = self%state%chem%charge, &
          ntype        = self%state%ntype, &
          nproc        = self%state%nproc, &
          dt           = self%state%params%dt, &
          zg_sec       = self%state%dom%zg_sec, &
          n_cath       = self%state%dom%n_cath, &
          dir_sec_inout = self%state%dom%dir_sec, &
          tag_neg      = int(self%state%chem%tag_neg,  int32), &
          tag_beam     = int(self%state%chem%tag_beam, int32), &
          flag_pbc     = int(self%state%dom%flag_pbc,  int32), &
          flag_pbcz    = int(self%state%dom%flag_pbcz, int32), &
          iseed        = self%state%params%iseed(iproc_h), &
          N_inj        = self%state%N_inj, &
          sour_xy      = self%state%sour_xy, &
          sour_xz      = self%state%sour_xz, &
          sour_yz      = self%state%sour_yz, &
          iz_pl        = self%state%params%iz_plot_plane, &
          ix_pl        = self%state%params%ix_plot_plane, &
          P_loss       = self%state%P_loss, &
          mom_loss     = self%state%mom_loss )
      end do
      !$omp end parallel do
      self%state%N_inj = 0_int32
    end if

    ! Flux injection of negative ions off extraction electrode
    ! (legacy: part_flux_injection, every step when jne>0)
    if (self%state%cfg%jne > 0.0_real64) then
      call inject_flux_particles( &
        part         = self%state%part, &
        cfg          = self%state%cfg, &
        n            = int(self%state%dom%n, int32), &
        h            = self%state%dom%h, &
        bcnd         = self%state%dom%bcnd, &
        Nm           = self%state%params%Nm, &
        mass         = self%state%chem%mass, &
        ntype        = self%state%ntype, &
        nproc        = self%state%nproc, &
        mpi_rank     = self%state%mpi_rank, &
        mpi_size     = self%state%mpi_size, &
        dt           = self%state%params%dt, &
        istep        = istep, &
        nsav         = self%state%cfg%nsav, &
        xg1          = self%state%dom%xg1, &
        Lgy          = self%state%dom%Lgy, &
        Lgz          = self%state%dom%Lgz, &
        ymax         = self%state%dom%ymax, &
        zmax         = self%state%dom%zmax, &
        Sg           = self%state%dom%Sg, &
        tag_neg      = int(self%state%chem%tag_neg, int32), &
        tag_neu      = int(self%state%chem%tag_neu, int32), &
        iseed        = self%state%params%iseed, &
        sour_xy      = self%state%sour_xy, &
        sour_xz      = self%state%sour_xz, &
        sour_fx_yz   = self%state%sour_fx_yz, &
        iz_pl        = self%state%params%iz_plot_plane, &
        ix_pl        = self%state%params%ix_plot_plane, &
        P_loss       = self%state%P_loss, &
        N_flx        = self%state%N_flx, &
        mom_loss     = self%state%mom_loss )
      self%state%N_flx = 0_int32
    end if

    ! Output/write
    t0 = MPI_Wtime()
    if (mod(istep, self%state%cfg%nsav) == 1_int32) then
      call self%output_step(istep)
      call self%reset_2d_averages()
    end if
    t1 = MPI_Wtime()
    self%t_avg = self%t_avg + (t1 - t0)

    ! Backup (legacy: it.gt.1 .and. MOD(it,nbak).eq.1)
    t0 = MPI_Wtime()
    if (self%state%cfg%nbak > 0_int32) then
      if (istep > 1_int32 .and. mod(istep, self%state%cfg%nbak) == 1_int32) then
        if (self%state%mpi_rank == 0_int32) write(*,*) 'Backing up simulation data ...'
        call write_restart_modular( &
          mpi_rank = self%state%mpi_rank, &
          ntype    = self%state%ntype, &
          nproc    = self%state%nproc, &
          tag_neg  = int(self%state%chem%tag_neg, int32), &
          part     = self%state%part, &
          time     = real(istep, real64) * self%state%params%dt, &
          flag_wrt = self%flag_wrt, &
          flag_RFant = self%state%cfg%flag_RFant, &
          E0_RF      = self%state%cfg%E0_RF )
        if (self%state%mpi_rank == 0_int32) write(*,*) 'Done!'
      end if
    end if
    t1 = MPI_Wtime()
    self%t_bck = self%t_bck + (t1 - t0)

    ! Total wall time
    t1 = MPI_Wtime()
    self%t_total = self%t_total + (t1 - tstep0)
    self%t_count = self%t_count + 1_int32

  end subroutine advance_one_step


  subroutine run(self, nsteps)
    class(Simulation), intent(inout) :: self
    integer, intent(in) :: nsteps

    integer(int32) :: istep, istep_offset, istep_end
    integer(int32) :: ptype
    character(len=8) :: pnum

    if (self%state%mpi_rank == 0) then
      write(*,*) " "
      write(*,*) " "
      write(*,"(a)") ">>> Entering the PIC loop <<<"
      write(*,*) " "
    end if

    ! Fresh (non-restart) run: clear out any stale scalar time-series files
    ! left over from a previous run in the same Output/ directory, so
    ! print_diagnostics/output_step's own file_exists check (below) sees a
    ! genuinely-absent file and starts each one anew (fresh header, time
    ! series from t=0) instead of silently appending a new run's data onto
    ! an old one. A real restart (flag_restart>0) must NOT do this - those
    ! files are expected to already hold the original run's history and are
    ! meant to just keep growing.
    if (self%state%cfg%flag_restart == 0_int32 .and. self%state%mpi_rank == 0_int32) then
      call reset_output_file('./Output/Ptot.dat')
      call reset_output_file('./Output/Mtot.dat')
      call reset_output_file('./Output/phi_Te_ne_cntr.dat')
      do ptype = 1_int32, self%state%ntype
        write(pnum,'(i0)') ptype
        call reset_output_file('./Output/Iw'//trim(pnum)//'.dat')
        call reset_output_file('./Output/Pw'//trim(pnum)//'.dat')
      end do
    end if

    ! On restart, self%state%time is set from the backup file's stored
    ! physical time (= istep_original * dt). Convert back to a step offset
    ! so istep numbering, timing diagnostics, and backup cadence all
    ! continue from where the original run left off instead of resetting
    ! to 1. For a fresh run self%state%time == 0 so istep_offset == 0.
    istep_offset = nint(self%state%time / self%state%params%dt, int32)
    istep_end    = istep_offset + int(nsteps, int32)

    do istep = istep_offset + 1_int32, istep_end

      if (mod(istep, self%state%cfg%nsav) == 1_int32) then
        call self%print_diagnostics(istep)
      end if

      call self%advance_one_step(istep)

    end do
  end subroutine run


  ! Deletes filename if it exists; a no-op otherwise. Used at the start of a
  ! fresh (non-restart) run to clear stale Output/ scalar time-series files
  ! from a previous run before the append-based writes further down
  ! (print_diagnostics, output_step) start appending to them.
  subroutine reset_output_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit
    logical :: file_exists

    inquire(file=filename, exist=file_exists)
    if (file_exists) then
      open(newunit=unit, file=filename, status='old', action='write')
      close(unit, status='delete')
    end if
  end subroutine reset_output_file


  subroutine finalize(self)
    class(Simulation), intent(inout) :: self

    call self%state%finalize()
  end subroutine finalize

  subroutine print_diagnostics(self,istep)
    class(Simulation), intent(inout) :: self
    integer(int32), intent(in)   :: istep
    integer(int32)               :: ptype
    real(real64)                 :: simulation_time
    real(real64) :: Pwall, Pabs, Pcoll, Pinj
    real(real64) :: Iw1, Iw2
    real(real64) :: Ekw1, Ekw2
    integer :: ierr
    integer(int32) :: nspecies_print
    integer(int32), allocatable :: npart_global(:)
    ! self%state%P_loss/p_mac/part(:,:) hold only THIS MPI rank's local
    ! subdomain, so the sum(...)'s below (unchanged from before - they
    ! already summed across this rank's own OMP threads/iproc) are still
    ! only a per-rank-local total. raw(:) packs those local totals so ONE
    ! MPI_Allreduce can turn them into true global totals before anything
    ! is printed - same two-step "reduce across MPI, THEN rank-0 prints"
    ! order legacy uses for its own equivalent diagnostics (sum_np_tot's
    ! MPI_ALLREDUCE followed by "if(mpi_rank.eq.0) write(*,101) ...",
    ! Src/main.f90). Without this, each rank would print only its own
    ! local slice of particles/power/current, not the physically
    ! meaningful whole-domain total.
    real(real64) :: raw(9)
    ! Per-(species,grid-region) current/power, legacy's Iw<i>.dat/Pw<i>.dat -
    ! same two-step MPI reduction as raw(:) above, but keeping the full
    ! (ntype,2,0:ngrid) shape instead of collapsing to scalars.
    real(real64), allocatable :: p_mac_global(:,:,:)
    integer(int32) :: igrid, kn_start

    ! Momentum-conservation diagnostic - vector (x,y,z) counterpart of
    ! Pwall/Pabs/Pcoll/Pinj, same two-step MPI reduction. raw_mom groups
    ! by component: (1:3)=wall, (4:6)=abs/RF, (7:9)=coll, (10:12)=inj -
    ! mirrors self%state%mom_loss's (axis,slot,ptype,iproc) layout.
    ! mom_total is the instantaneous Sum(Nm*mass*v) over every live tracked
    ! macroparticle: printed/written alongside Fwall/Fabs/Fcoll/Finj so the
    ! balance Fwall*dtp + Fabs*dtp + Fcoll*dtp + Finj =~= Delta(mom_total)
    ! between two prints can actually be checked, not just eyeballed -
    ! same sign convention as Pwall/Pabs/Pcoll/Pinj (Fwall negated like
    ! Pwall, so all four terms read as "signed contribution to the
    ! plasma's momentum" and a plain sum over one period equals the
    ! change in mom_total if the books are closed).
    real(real64) :: raw_mom(12)
    real(real64) :: mom_total(3)
    real(real64) :: Fwall_mom(3), Fabs_mom(3), Fcoll_mom(3), Finj_mom(3)
    integer(int32) :: iproc_m, ip_m

    nspecies_print   = self%state%rxn%ntype - self%state%rxn%n_neu
    simulation_time  = self%state%params%dt * real(istep, real64)
    allocate(npart_global(nspecies_print))
    do ptype = 1_int32, nspecies_print
      npart_global(ptype) = sum(self%state%part(ptype,:)%n)
    end do

    allocate(p_mac_global(self%state%ntype, 2_int32, 0:self%state%cfg%ngrid))
    p_mac_global = sum(self%state%p_mac, dim=4)

    raw(1) = sum(self%state%P_loss(1,:,:))
    raw(2) = sum(self%state%P_loss(2,:,:))
    raw(3) = sum(self%state%P_loss(3,:,:))
    raw(4) = sum(self%state%P_loss(4,:,:))
    raw(5) = sum(self%state%p_mac(1,1,:,:))
    raw(6) = sum(self%state%p_mac(2:self%state%ntype,1,:,:))
    raw(7) = sum(self%state%p_mac(1,2,:,:))
    raw(8) = sum(self%state%p_mac(2:self%state%ntype,2,:,:))
    ! sum(abs(x_i)) over the global particle set equals the sum, across
    ! ranks, of each rank's own local sum(abs(its particles)) - abs() is
    ! applied per-element before any summing, so this reduces the same
    ! way as the plain sums above (unlike abs(sum(x_i)), which is NOT
    ! separable across ranks - Ekw1 below instead takes abs() of raw(5)
    ! AFTER raw(5) is already the global sum, matching how Iw1/Iw2 treat
    ! their own sums).
    raw(9) = sum(abs(self%state%p_mac(2:self%state%ntype,1,:,:)))

    raw_mom(1)  = sum(self%state%mom_loss(1,1,:,:))
    raw_mom(2)  = sum(self%state%mom_loss(2,1,:,:))
    raw_mom(3)  = sum(self%state%mom_loss(3,1,:,:))
    raw_mom(4)  = sum(self%state%mom_loss(1,2,:,:))
    raw_mom(5)  = sum(self%state%mom_loss(2,2,:,:))
    raw_mom(6)  = sum(self%state%mom_loss(3,2,:,:))
    raw_mom(7)  = sum(self%state%mom_loss(1,3,:,:))
    raw_mom(8)  = sum(self%state%mom_loss(2,3,:,:))
    raw_mom(9)  = sum(self%state%mom_loss(3,3,:,:))
    raw_mom(10) = sum(self%state%mom_loss(1,4,:,:))
    raw_mom(11) = sum(self%state%mom_loss(2,4,:,:))
    raw_mom(12) = sum(self%state%mom_loss(3,4,:,:))

    ! Instantaneous total momentum of every live tracked macroparticle on
    ! this rank (mass<=0 marks an untracked background species with no
    ! ParticleSet storage - same guard the mover uses to skip the push).
    mom_total = 0.0_real64
    do ptype = 1_int32, self%state%ntype
      if (self%state%chem%mass(ptype) <= 0.0_real64) cycle
      do iproc_m = 1_int32, self%state%nproc
        if (.not. allocated(self%state%part(ptype,iproc_m)%pv)) cycle
        do ip_m = 1_int32, self%state%part(ptype,iproc_m)%n
          if (allocated(self%state%part(ptype,iproc_m)%flag_dead)) then
            if (self%state%part(ptype,iproc_m)%flag_dead(ip_m) /= 0_int8) cycle
          end if
          mom_total(1) = mom_total(1) + self%state%params%Nm(ptype) * self%state%chem%mass(ptype) * &
              self%state%part(ptype,iproc_m)%pv(4,ip_m)
          mom_total(2) = mom_total(2) + self%state%params%Nm(ptype) * self%state%chem%mass(ptype) * &
              self%state%part(ptype,iproc_m)%pv(5,ip_m)
          mom_total(3) = mom_total(3) + self%state%params%Nm(ptype) * self%state%chem%mass(ptype) * &
              self%state%part(ptype,iproc_m)%pv(6,ip_m)
        end do
      end do
    end do

    if (self%state%mpi_size > 1_int32) then
      call MPI_Allreduce(MPI_IN_PLACE, npart_global, nspecies_print, &
                          MPI_INTEGER, MPI_SUM, self%state%comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, raw, 9_int32, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, self%state%comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, raw_mom, 12_int32, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, self%state%comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, mom_total, 3_int32, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, self%state%comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, p_mac_global, &
                          size(p_mac_global, kind=int32), &
                          MPI_DOUBLE_PRECISION, MPI_SUM, self%state%comm, ierr)
    end if

    ! Same normalization Iw1/Iw2 already apply to raw(5)/raw(6) below -
    ! p_mac is already in physical (current/power) units at accumulation
    ! time (see mod_state.f90), just averaged over the save interval here.
    p_mac_global = p_mac_global / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Pabs  = raw(2) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Pcoll = raw(3) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Pinj  = raw(4) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Pwall = -raw(1) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Iw1 = raw(5) / &
          (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Iw2 = raw(6) / &
          (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Ekw1 = 0.0_real64
    if (abs(raw(5)) > 0.0_real64) Ekw1 = raw(7) / abs(raw(5))

    Ekw2 = 0.0_real64
    if (raw(9) > 0.0_real64) Ekw2 = raw(8) / raw(9)

    Fwall_mom = -raw_mom(1:3) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Fabs_mom = raw_mom(4:6) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Fcoll_mom = raw_mom(7:9) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    Finj_mom = raw_mom(10:12) / &
        (real(self%state%cfg%nsav - 1_int32, real64) * self%state%params%dt)

    ! Everything below this point is output only (no further physics
    ! state is touched) - restrict it to rank 0, matching legacy's own
    ! convention of a single rank-0-only diagnostic print instead of one
    ! copy per MPI rank. print_rxn_counts/print_debug_diagnostics/
    ! print_poisson_breakdown are similarly rank-0-only calls now; their
    ! matching reset_* calls still run on every rank below (unconditional,
    ! after this block), since every rank must keep zeroing its own local
    ! accumulators each period regardless of who prints.
    if (self%state%mpi_rank == 0_int32) then

    call print_rxn_counts()
    call print_debug_diagnostics()

    write(*,'(a)') " "
    write(*,'(a)') " "

    write(*,'(a,i0,a,F8.3,a)') &
      ' TIME STEP ', istep, ' --> ', &
      self%state%params%dt*istep*1.d6, " us"

    write(*,'(a)') " -------------------------"
    write(*,'(a)') " Particles diagnostics: "

    do ptype = 1_int32, nspecies_print
      write(*,'(a,a,a,i0)') &
        " npart ", self%state%chem%pname(ptype) , &
        " = ", npart_global(ptype)
    end do

    write(*,'(a,i12,a)') " ===== MODULAR DIAGNOSTIC it = ", istep, " ====="
    write(*,'(a,3ES24.15)') " phi min/max/sum = ", &
      minval(self%state%fld%phi), &
      maxval(self%state%fld%phi), &
      sum(self%state%fld%phi)
    write(*,'(a)') " ==========================================="

    write(*,'(a)') " -------------------------"

    write(*,"(a)") " Power diagnostics: "

    write(*,'(A,ES10.2)') ' Pwall (W)  = ', Pwall
    write(*,'(A,ES10.2)') ' Pabs  (W)  = ', Pabs
    write(*,'(A,ES10.2)') ' Pcoll (W)  = ', Pcoll
    write(*,'(A,ES10.2)') ' Pinj  (W)  = ', Pinj

    write(*,'(a,ES10.2,ES10.2)') ' I_w  (A)  = ', Iw1, Iw2
    write(*,'(a,ES10.2,ES10.2)') ' Ek_w (eV) = ', Ekw1, Ekw2

    write(*,'(a)') " -------------------------"
    write(*,"(a)") " Momentum diagnostics (x,y,z): "
    write(*,'(A,3ES10.2)') ' Fwall (N)  = ', Fwall_mom
    write(*,'(A,3ES10.2)') ' Fabs  (N)  = ', Fabs_mom
    write(*,'(A,3ES10.2)') ' Fcoll (N)  = ', Fcoll_mom
    write(*,'(A,3ES10.2)') ' Finj (kg.m/s) = ', Finj_mom
    write(*,'(A,3ES10.2)') ' mom_total (kg.m/s) = ', mom_total

    ! --- legacy DATA/*.dat scalar time series, kept in Output/ instead
    ! (this is the modular tree, not legacy's DATA/ layout). I_inj and
    ! Vgrd(igrid_sec) (legacy's RF-antenna injected current and per-grid
    ! bias voltage) are not tracked anywhere in modular yet, so those two
    ! legacy columns are omitted here rather than zero-padded - a script
    ! expecting legacy's exact 11-column Ptot.dat will need updating.
    block
      integer :: ufile
      logical :: file_exists

      inquire(file='./Output/Ptot.dat', exist=file_exists)
      open(newunit=ufile, file='./Output/Ptot.dat', status='unknown', &
           position='append', action='write')
      if (.not. file_exists) then
        write(ufile,'(a)') '# Time (s), Pwall (W), Pabs, Pcoll, Pinj, I_mw (A), I_pw, Ek_ew (eV), Ek_iw'
      end if
      write(ufile,'(20(1x,es16.8))') simulation_time, Pwall, Pabs, Pcoll, Pinj, Iw1, Iw2, Ekw1, Ekw2
      close(ufile)
    end block

    ! --- momentum-conservation diagnostic: same signed-flux convention as
    ! Ptot.dat above (Fwall negated like Pwall, Finj not time-normalized
    ! like Pinj), plus mom_total so Fwall*dtp+Fabs*dtp+Fcoll*dtp+Finj can be
    ! checked against the actual Delta(mom_total) between two rows.
    block
      integer :: ufile
      logical :: file_exists

      inquire(file='./Output/Mtot.dat', exist=file_exists)
      open(newunit=ufile, file='./Output/Mtot.dat', status='unknown', &
           position='append', action='write')
      if (.not. file_exists) then
        write(ufile,'(a)') '# Time (s), Fwall_x, Fwall_y, Fwall_z (N), '// &
            'Fabs_x, Fabs_y, Fabs_z (N), Fcoll_x, Fcoll_y, Fcoll_z (N), '// &
            'Finj_x, Finj_y, Finj_z (kg.m/s), mom_total_x, mom_total_y, mom_total_z (kg.m/s)'
      end if
      write(ufile,'(20(1x,es16.8))') simulation_time, &
          Fwall_mom, Fabs_mom, Fcoll_mom, Finj_mom, mom_total
      close(ufile)
    end block

    ! --- legacy Iw<i>.dat/Pw<i>.dat: per-species, per-grid-region
    ! current/power time series (p_mac_global reduced/normalized above).
    kn_start = 1_int32
    if (self%state%dom%flag_nmn == 1_int32) kn_start = 0_int32

    do ptype = 1_int32, self%state%ntype
      block
        integer :: uIw, uPw
        character(len=8) :: pnum
        write(pnum,'(i0)') ptype
        open(newunit=uIw, file='./Output/Iw'//trim(pnum)//'.dat', status='unknown', &
             position='append', action='write')
        open(newunit=uPw, file='./Output/Pw'//trim(pnum)//'.dat', status='unknown', &
             position='append', action='write')
        write(uIw,'(20(1x,es16.8))') simulation_time, &
            ( p_mac_global(ptype,1,igrid), igrid=kn_start, self%state%cfg%ngrid )
        write(uPw,'(20(1x,es16.8))') simulation_time, &
            ( p_mac_global(ptype,2,igrid), igrid=kn_start, self%state%cfg%ngrid )
        close(uIw)
        close(uPw)
      end block
    end do

    write(*,'(a)') " -------------------------"

    if (self%t_count > 0_int32) then
      write(*,'(a)') " Timing diagnostics: "

      write(*,'(A,F7.1)') ' E/rho      (ms) = ', &
        1000.0_real64*self%t_Erho/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' poisson    (ms) = ', &
        1000.0_real64*self%t_poisson/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' sorting    (ms) = ', &
        1000.0_real64*self%t_sort/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' avg/write  (ms) = ', &
        1000.0_real64*self%t_avg/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' MC         (ms) = ', &
        1000.0_real64*self%t_MC/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' mover      (ms) = ', &
        1000.0_real64*self%t_mover/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' bck        (ms) = ', &
        1000.0_real64*self%t_bck/real(self%t_count,real64)

      write(*,'(A,F7.1)') ' total      (ms) = ', &
        1000.0_real64*self%t_total/real(self%t_count,real64)

      write(*,'(a)') " -------------------------"
    end if

    call print_poisson_breakdown(self%t_count)

    if (self%t_count > 0_int32) then
      write(*,'(a)') " ----- mover / deposit timing breakdown -----"
      write(*,'(a,f8.2,a,f8.2)') &
        "  mover_push(ms)=", 1000.0_real64*self%t_mover_push/real(self%t_count,real64), &
        "  mover_bc(ms)=",   1000.0_real64*self%t_mover_bc/real(self%t_count,real64)

      ! TEMPORARY diagnostic - see move_and_bc_boris's header comment
      ! (mod_particleMover.f90). Only meaningful when the Boris pusher ran
      ! at all this window (n_boris_total==0 otherwise, e.g. an
      ! electrostatic-only case).
      if (self%n_boris_total > 0_int64) then
        write(*,'(a,i0,a,i0,a,f6.2,a)') &
          "  boris used=", self%n_boris_used, " / ", self%n_boris_total, &
          "  (", 100.0_real64*real(self%n_boris_used,real64)/real(self%n_boris_total,real64), &
          "% ran full Boris rotation, rest hit the negligible-|B| fallback)"
      end if
      write(*,'(a,f8.2,a,f8.2,a,f8.2)') &
        "  dep_clear(ms)=",  1000.0_real64*self%t_dep_clear/real(self%t_count,real64), &
        "  dep_loop(ms)=",   1000.0_real64*self%t_dep_loop/real(self%t_count,real64), &
        "  dep_reduce(ms)=", 1000.0_real64*self%t_dep_reduce/real(self%t_count,real64)

      ! OMP load-imbalance check across iproc: min/avg/max of the same
      ! per-thread timings mover_push(ms)/dep_loop(ms) above already
      ! reduce with max. imbalance = max/avg; close to 1 means threads are
      ! evenly loaded, so the mover_push(ms)/dep_loop(ms) cost above is
      ! genuine per-particle work, not a few slow iproc dragging the max up.
      block
        real(real64) :: push_min_ms, push_avg_ms, push_max_ms, push_imbalance
        real(real64) :: dep_min_ms, dep_avg_ms, dep_max_ms, dep_imbalance

        push_min_ms = 1000.0_real64*self%t_mover_push_min/real(self%t_count,real64)
        push_avg_ms = 1000.0_real64*self%t_mover_push_avg/real(self%t_count,real64)
        push_max_ms = 1000.0_real64*self%t_mover_push/real(self%t_count,real64)
        push_imbalance = 0.0_real64
        if (push_avg_ms > 0.0_real64) push_imbalance = push_max_ms / push_avg_ms

        dep_min_ms = 1000.0_real64*self%t_dep_loop_min/real(self%t_count,real64)
        dep_avg_ms = 1000.0_real64*self%t_dep_loop_avg/real(self%t_count,real64)
        dep_max_ms = 1000.0_real64*self%t_dep_loop/real(self%t_count,real64)
        dep_imbalance = 0.0_real64
        if (dep_avg_ms > 0.0_real64) dep_imbalance = dep_max_ms / dep_avg_ms

        write(*,'(a,i0,a)') " ----- OMP load imbalance across iproc (nproc=", self%state%nproc, &
          ", min/avg/max, max/avg ratio) -----"
        write(*,'(a,f8.2,a,f8.2,a,f8.2,a,f6.3)') &
          "  mover_push(ms) min=", push_min_ms, " avg=", push_avg_ms, &
          " max=", push_max_ms, " max/avg=", push_imbalance
        write(*,'(a,f8.2,a,f8.2,a,f8.2,a,f6.3)') &
          "  dep_loop(ms)   min=", dep_min_ms, " avg=", dep_avg_ms, &
          " max=", dep_max_ms, " max/avg=", dep_imbalance

        ! Per-iproc live particle counts (species 1) alongside the timing
        ! ratios above - min/avg/max alone can't distinguish "one thread
        ! genuinely owns ~0 particles" (self%part(ptype,iproc)%n <= 0, the
        ! mover_push cycle guard in advance_particles_local) from "every
        ! thread has particles but one is just slower". Printed every call
        ! since it's O(nproc), negligible next to the timing block above.
        block
          integer(int32) :: iproc_dbg
          ! Written value-by-value (not accumulated into one fixed-length
          ! internal-write buffer) since nproc can reach 192 on some
          ! machines - a bounded buffer overflowed here ("output statement
          ! overflows record") the first time this was tried at nproc=32.
          write(*,'(a)',advance='no') " npart[iproc] (species 1) ="
          do iproc_dbg = 1, self%state%nproc
            write(*,'(1x,i0)',advance='no') self%state%part(1,iproc_dbg)%n
          end do
          write(*,*)
        end block
      end block

      write(*,'(a)') " ---------------------------------------------"
    end if




    write(*,'(a)') "  "
    write(*,*) "  "

    end if ! self%state%mpi_rank == 0

    deallocate(npart_global)
    deallocate(p_mac_global)

    ! Reset legacy-style power accumulators after printing - unconditional
    ! on every rank (not just rank 0), so each rank's own local
    ! accumulators start clean for the next reporting period regardless
    ! of which rank did the printing above.
    call reset_rxn_counts()
    call reset_debug_diagnostics()
    if (allocated(self%state%P_loss) .or. allocated(self%state%p_mac)) then
      self%state%P_loss = 0.0_real64
      self%state%p_mac  = 0.0_real64
    end if
    if (allocated(self%state%mom_loss)) self%state%mom_loss = 0.0_real64

    ! Reset timers
    self%t_Erho    = 0.0_real64
    self%t_poisson = 0.0_real64
    call reset_poisson_breakdown()
    self%t_sort    = 0.0_real64
    self%t_avg     = 0.0_real64
    self%t_MC      = 0.0_real64
    self%t_mover   = 0.0_real64
    self%t_mover_push = 0.0_real64
    self%t_mover_bc   = 0.0_real64
    self%t_dep_clear  = 0.0_real64
    self%t_dep_loop   = 0.0_real64
    self%t_dep_reduce = 0.0_real64
    self%t_mover_push_min = 0.0_real64
    self%t_mover_push_avg = 0.0_real64
    self%t_dep_loop_min   = 0.0_real64
    self%t_dep_loop_avg   = 0.0_real64
    self%n_boris_used  = 0_int64
    self%n_boris_total = 0_int64
    self%t_bck     = 0.0_real64
    self%t_total   = 0.0_real64
    self%t_count   = 0_int32

  end subroutine print_diagnostics

end module mod_simulation