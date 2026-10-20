module mod_state
  use iso_fortran_env, only: int32, real64
  use omp_lib, only: omp_get_max_threads, omp_get_num_procs, omp_get_wtime
  use mpi

  use mod_config,         only: Config
  use mod_readConditions, only: read_input
  use mod_utils,          only: stop_calculation

  use mod_constants,      only: eps0, qe, c, pi
  use mod_domain,         only: Domain
  use mod_boundary,       only: build_boundary

  use mod_fields,         only: Fields
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_charge_weights, only: build_kq
  use mod_debug_checks,   only: checkpoint_poisson_decomp, checkpoint_kq
  use mod_density,        only: reduce_species_density, sync_species_density, clear_np_planes

  use mod_chemistryState, only: ChemistryState
  use mod_reactionsDB,    only: ReactionsDB
  use mod_part_info,      only: npart
  use mod_particles,      only: ParticleSet
  use mod_particle_loader,  only: load_particles_modular
  use mod_restart,          only: restart_particles_modular
  use mod_particle_sorting, only: sort_particles_by_cell, redistribute_species_by_cell, check_particles_are_sorted, check_cell_indexing
  use mod_particleMover,    only: move_and_bc_electrostatic, move_and_bc_electrostatic_fast, &
                                   move_and_bc_boris, move_and_bc_boris_fast
  use mod_particleBC,       only: apply_particle_bc, SeeParams
  use mod_chargeDeposition, only: deposit_particle_set_to_np_thread, HeatRegionTally
  use mod_magneticField,    only: MagneticField
  use mod_simParams,        only: SimParams
  use mod_heating,          only: apply_electron_heating
  use mod_planeMoments,     only: compute_particle_plane_moments_species

  implicit none
  private
  public :: State
  public :: advance_particles_local

  ! Production defaults. Keep expensive consistency checks available in the
  ! source, but do not run them in normal legacy-comparison/production runs.
  logical, parameter :: DEBUG_INIT_CHECKS = .false.
  logical, parameter :: VALIDATE_SORTING  = .false.

  type :: State
    type(Config)         :: cfg
    type(Domain)         :: dom
    type(Fields)         :: fld
    type(PoissonDecomp)  :: pdec

    ! Set from the .bak file's stored time on restart (cfg%flag_restart>0);
    ! left at 0 for a fresh load. Not currently consumed by the main step
    ! loop's own istep-based time bookkeeping - exposed here for callers
    ! that want to resynchronize elapsed-time reporting after a restart.
    real(real64) :: time = 0.0_real64

    type(ChemistryState) :: chem
    type(ReactionsDB)    :: rxn
    type(MagneticField)  :: magField
    type(SimParams)      :: params
    type(ParticleSet), allocatable :: part(:,:)

    real(real64), allocatable :: np_thread(:,:,:,:,:)

    ! Per-(ptype,iproc) z-extent of np_thread that may be nonzero: planes
    ! np_zlo..np_zhi, plus the end planes 0:2 and n3:n3+2 if np_zend /= 0
    ! (see np_plane_live in mod_density). Every plane outside it is
    ! exactly zero. Since each iproc owns a z-slab of particles
    ! (redistribute_species_by_cell), this lets the per-step clear and
    ! the thread reduction touch only that slab instead of the whole grid.
    integer(int32), allocatable :: np_zlo(:,:), np_zhi(:,:), np_zend(:,:)

    ! .false. while fld%np holds only this rank's density (the per-step
    ! deposit leaves it rank-local and MPI-sums only rho); sync_np() makes
    ! it global. Every reader of fld%np must call sync_np() first.
    logical :: np_synced = .true.
    real(real64), allocatable :: sum_q_xz(:,:,:,:)
    real(real64), allocatable :: sum_q_yz(:,:,:,:,:)

    integer :: mpi_rank = -1
    integer :: mpi_size = -1
    integer :: comm     = MPI_COMM_NULL

    integer(int32) :: ntype = 0
    integer(int32) :: nproc = 1
    integer(int32), allocatable :: Nh(:)
    real(real64), allocatable :: sum_dEk(:)
    ! Heating weights sum (= real(Nh) unless cfg%delta_heat > 0)
    real(real64), allocatable :: sum_wh(:)

    real(real64), allocatable :: P_loss(:,:,:)
    real(real64), allocatable :: p_mac(:,:,:,:)

    ! Momentum-conservation diagnostic - vector (x,y,z) counterpart of
    ! P_loss, same (slot,ptype,iproc) convention: mom_loss(:,1,:,:) wall,
    ! mom_loss(:,2,:,:) RF-antenna impulse (folded in from mom_RF the same
    ! way P_RF folds into P_loss(2,:,:)), mom_loss(:,3,:,:) collisions
    ! (incl. the dead-particle correction below), mom_loss(:,4,:,:)
    ! injection+SEE. See print_diagnostics (mod_simulation.f90) for how
    ! this is reduced/printed/reset, mirroring Pwall/Pabs/Pcoll/Pinj.
    real(real64), allocatable :: mom_loss(:,:,:,:)
    real(real64), allocatable :: mom_RF(:,:,:)

    ! RF antenna (inductive) heating - flag_RFant==1 only. gams: dynamic
    ! skin depth (m), refined every ns_RF steps by update_rf_convergence
    ! from the actual local electron density. P_RF(ptype,iproc): power
    ! absorbed via the RF field, accumulated every step by the mover and
    ! folded into P_loss(2,:,:)/reset every ns_RF steps.
    real(real64) :: gams = 0.0_real64
    real(real64), allocatable :: P_RF(:,:)
    ! Number of mover steps accumulated in P_RF since the last
    ! update_rf_convergence. Drives the convergence cadence instead of
    ! mod(istep,ns_RF), so a restart (istep offset not a multiple of
    ! ns_RF) still averages over a full RF period.
    integer(int32) :: nstep_RF = 0_int32

    ! Injection source diagnostic arrays — same shape as legacy:
    !   sour_xy/xz/yz: volume injection (part_injection)
    !   sour_fx_yz:    flux injection off extraction electrode (part_flux_injection)
    !   N_inj: per-(ptype,iproc) injected particle counter (opt_inj==2 mode)
    !   N_flx: per-(ptype,iproc) flux counter (H- wall-loss driven injection)
    integer(int32), allocatable :: sour_xy(:,:,:,:)
    integer(int32), allocatable :: sour_xz(:,:,:,:)
    integer(int32), allocatable :: sour_yz(:,:,:,:)
    integer(int32), allocatable :: sour_fx_yz(:,:,:,:)
    integer(int32), allocatable :: N_inj(:,:)
    integer(int32), allocatable :: N_flx(:,:)

    real(real64), allocatable :: data_pavg_xy(:,:,:,:)
    real(real64), allocatable :: data_pavg_xz(:,:,:,:)
    real(real64), allocatable :: data_pavg_yz(:,:,:,:)

    real(real64), allocatable :: np_avg_xy(:,:,:)
    real(real64), allocatable :: np_avg_xz(:,:,:)
    real(real64), allocatable :: np_avg_yz(:,:,:)

    real(real64), allocatable :: phi_avg_xy(:,:)
    real(real64), allocatable :: phi_avg_xz(:,:)
    real(real64), allocatable :: phi_avg_yz(:,:)

    ! 3D time averages (legacy avg3D, cfg%flag_avg3D==1 only):
    ! (1:3,...) = E, (4,...) = phi, (5,...) = density of species 2, on the
    ! physical nodes 1..n+1. Summed at the same cadence and with the same
    ! cnt_avg as the 2D plane averages; written by output_step (averaged
    ! E planes + phi_n_3D.dat). Rank 0 only - every rank holds the same
    ! global E/phi/np at that point and only rank 0 writes, so the other
    ! ranks would just duplicate ~100 MB and the accumulation.
    real(real64), allocatable :: avg3D(:,:,:,:)

    ! Reaction-driven particle production/destruction per plane, per
    ! tracked species, per OMP thread (iproc) - legacy's plt_src==1
    ! ss2D_xy/xz/yz. Dimensioned per-iproc (like Pcoll/p_mac) because they
    ! are written concurrently from perform_collisions_gwenael's
    ! "!$omp parallel do ... iproc" loop; summed over iproc only at
    ! output/write time. Accumulate continuously across steps, reset only
    ! in reset_2d_averages (same cadence as data_pavg_xy).
    real(real64), allocatable :: sour_avg_xy(:,:,:,:)
    real(real64), allocatable :: sour_avg_xz(:,:,:,:)
    real(real64), allocatable :: sour_avg_yz(:,:,:,:)
    real(real64), allocatable :: sink_avg_xy(:,:,:,:)
    real(real64), allocatable :: sink_avg_xz(:,:,:,:)
    real(real64), allocatable :: sink_avg_yz(:,:,:,:)

    integer(int32) :: cnt_avg = 0_int32

  contains
    procedure :: init
    procedure :: init_chemistry
    procedure :: init_domain_and_boundary
    procedure :: init_particles

    procedure :: sort_particles_local
    procedure :: apply_electron_heating_local
    procedure :: compute_heating_region_moments
    procedure :: update_heating_skin_depth
    procedure :: write_skin_depth_record
    procedure :: update_heating_vt
    procedure :: update_rf_convergence
    procedure :: compute_plane_moments_local
    procedure :: apply_dielectric_bc_to_phi
    procedure :: accumulate_2d_averages
    procedure :: sync_np

    procedure :: finalize_particles_only
    procedure :: finalize
  end type State

contains

  subroutine init(self, comm_in)
    class(State), intent(inout) :: self
    integer,      intent(in)    :: comm_in
    integer :: ierr

    self%comm = comm_in
    call MPI_Comm_rank(self%comm, self%mpi_rank, ierr)
    call MPI_Comm_size(self%comm, self%mpi_size, ierr)

    call read_input(self%cfg, self%mpi_rank)

    eps0 = self%cfg%k_eps0 * 8.854187817d-12
    if (self%mpi_rank == 0 .and. self%cfg%k_eps0 /= 1.0_real64) then
      write(*,*) 'eps0 HAS BEEN RE-SCALED: k_eps0 = ', self%cfg%k_eps0
    end if

    call self%dom%init_from_config(self%cfg)
    call self%dom%allocate_masks_domain()

    call self%fld%allocate_from_domain(self%dom)
    call self%fld%zero()

    call self%init_chemistry()
    call self%init_domain_and_boundary()

    call self%params%build(self%cfg, self%dom, self%chem, self%rxn, self%magField, self%mpi_rank)

    self%ntype = int(self%rxn%ntype - self%rxn%n_neu, int32)
    if (self%ntype < 1_int32) self%ntype = 1_int32

    ! Matches legacy (Src/main.f90:143): particle-domain parallelism always
    ! sizes to the real OMP runtime thread count, not the input file's
    ! "# of OMP threads" field (cfg%omp_rank_max), which is a separate
    ! legacy knob used only for restart redistribution math (see
    ! mod_restart.f90's nproc_written). Sizing self%nproc from
    ! cfg%omp_rank_max instead left most threads idle in the mover/deposit
    ! loops (or crashed outright) whenever that stale input-file value
    ! didn't match OMP_NUM_THREADS.
    self%nproc = max(1_int32, omp_get_max_threads())
    call self%params%init_seeds(self%nproc, self%mpi_rank)
    call self%params%print_summary(self%mpi_rank, self%cfg%nsav)

    ! RF antenna heating: R_ahp (antenna radius) must be set here whenever
    ! flag_RFant==1, independent of flag_circxh - apply_electron_heating_local's
    ! lazy R_ahp assignment never runs under RF-antenna mode since that
    ! forces flag_heat=0 (mod_readConditions.f90). Mirrors legacy
    ! Src/main.f90:312.
    if (self%cfg%flag_circxh == 1 .or. self%cfg%flag_RFant == 1) then
      self%cfg%R_ahp = self%cfg%yr_pow - self%dom%ymax / 2.0_real64
    end if

    ! Skin-depth heating needs the circular heating cross-section (inner
    ! side) for r and R_ahp to mean anything.
    if (self%cfg%flag_heat == 1 .and. self%cfg%delta_heat > 0.0_real64) then
      if (self%cfg%flag_circxh /= 1 .or. self%cfg%flag_ahp /= 0) then
        if (self%mpi_rank == 0) write(*,'(a)') 'Skin-depth heating (delta_heat > 0) ' // &
          'needs a circular cross-section spanning the heating region, please correct ...'
        call stop_calculation
      end if
      if (self%mpi_rank == 0) write(*,'(1x,a,f7.3,a,f7.3)') &
        'Skin-depth heating: w(r) = exp((r-R)/delta), R(cm)=', self%cfg%R_ahp*1.0e2_real64, &
        ', delta(cm)=', self%cfg%delta_heat*1.0e2_real64
      if (self%cfg%flag_delta_sc == 1 .and. self%mpi_rank == 0) write(*,'(1x,a,i0,a)') &
        'Skin depth recomputed every ', self%params%ns_RF, &
        ' steps from <ne> over R-delta <= r <= R (delta above = initial guess)'
    end if
    ! Planar antenna: an explicit window radius (may be smaller than the
    ! domain cross-section) overrides the yr_pow-derived radius.
    if (self%cfg%flag_planar_ant == 1 .and. self%cfg%R_win > 0.0_real64) then
      self%cfg%R_ahp = self%cfg%R_win
    end if

    if (self%cfg%flag_RFant == 1) then
      ! Initial nominal skin depth from the startup wp - mirrors legacy
      ! Src/main.f90:661-667. Refined every ns_RF steps by
      ! update_rf_convergence (called from mod_simulation.f90's
      ! advance_one_step).
      self%gams = c / self%params%wp / sqrt(self%cfg%k_eps0)

      if (self%mpi_rank == 0) then
        write(*,'(1x,a,i0,a,f6.2)') &
          '# of time steps per RF period=', self%params%ns_RF, &
          ', RF skin depth gam(cm)=', self%gams*1.0e2_real64
        if (self%cfg%flag_planar_ant == 1) then
          write(*,'(1x,a,f7.3,a,f7.3)') 'RF antenna: Planar, window radius (cm)=', &
            self%cfg%R_ahp*1.0e2_real64, ', window at x (cm)=', &
            max(0.0_real64, self%cfg%xl_pow)*1.0e2_real64
        else
          write(*,'(1x,a,f7.3)') 'RF antenna: Wrapped, coil radius (cm)=', &
            self%cfg%R_ahp*1.0e2_real64
        end if
      end if
      if (self%cfg%R_ahp <= 0.0_real64) then
        if (self%mpi_rank == 0) write(*,'(a)') 'RF antenna radius must be > 0 ' // &
          '(set yr_pow > ymax/2, or R_win for Planar), please correct ...'
        call stop_calculation
      end if

      ! Round nsav to a multiple of the RF period - legacy Src/main.f90:666.
      self%cfg%nsav = self%params%ns_RF * max(1_int32, &
           nint(real(self%cfg%nsav, real64) / real(self%params%ns_RF, real64), int32))
    end if

    allocate(self%np_thread(0:self%dom%n(1)+2, &
                            0:self%dom%n(2)+2, &
                            0:self%dom%n(3)+2, &
                            self%ntype, self%nproc))
    self%np_thread = 0.0_real64

    ! Start as "whole grid may be nonzero": the initial load and restart
    ! deposits write np_thread without tracking extents; the first
    ! per-step deposit clears everything and narrows these.
    allocate(self%np_zlo(self%ntype, self%nproc), self%np_zhi(self%ntype, self%nproc), &
             self%np_zend(self%ntype, self%nproc))
    self%np_zlo  = 0_int32
    self%np_zhi  = int(self%dom%n(3), int32) + 2_int32
    self%np_zend = 1_int32

    call self%init_particles()

    allocate(self%sum_q_xz(0:self%dom%n(1)+2, &
                           0:self%dom%n(3)+2, &
                           self%ntype, self%nproc))
    self%sum_q_xz = 0.0_real64

    allocate(self%sum_q_yz(2, &
                           0:self%dom%n(2)+2, &
                           0:self%dom%n(3)+2, &
                           self%ntype, self%nproc))
    self%sum_q_yz = 0.0_real64

    allocate(self%P_loss(4, self%ntype, self%nproc))
    self%P_loss = 0.0_real64

    allocate(self%mom_loss(3, 4, self%ntype, self%nproc))
    self%mom_loss = 0.0_real64

    allocate(self%P_RF(self%ntype, self%nproc))
    self%P_RF = 0.0_real64

    allocate(self%mom_RF(3, self%ntype, self%nproc))
    self%mom_RF = 0.0_real64

    allocate(self%p_mac(self%ntype, 2, 0:self%cfg%ngrid, self%nproc))
    self%p_mac = 0.0_real64

    ! Injection source arrays and counters
    associate(nx => self%dom%n(1), ny => self%dom%n(2), nz => self%dom%n(3), &
              nt => self%ntype, np => self%nproc)
      allocate(self%sour_xy(0:nx+2, 0:ny+2, nt, np));   self%sour_xy    = 0_int32
      allocate(self%sour_xz(0:nx+2, 0:nz+2, nt, np));   self%sour_xz    = 0_int32
      allocate(self%sour_yz(0:ny+2, 0:nz+2, nt, np));   self%sour_yz    = 0_int32
      allocate(self%sour_fx_yz(0:ny+2, 0:nz+2, nt, np)); self%sour_fx_yz = 0_int32
      allocate(self%N_inj(nt, np));  self%N_inj = 0_int32
      allocate(self%N_flx(nt, np));  self%N_flx = 0_int32
    end associate

    allocate(self%Nh(self%nproc))
    allocate(self%sum_dEk(self%nproc))
    allocate(self%sum_wh(self%nproc))
    self%Nh      = 0_int32
    self%sum_dEk = 0.0_real64
    self%sum_wh  = 0.0_real64

    allocate(self%np_avg_xy(0:self%dom%n(1)+2,0:self%dom%n(2)+2,self%ntype))
    allocate(self%np_avg_xz(0:self%dom%n(1)+2,0:self%dom%n(3)+2,self%ntype))
    allocate(self%np_avg_yz(0:self%dom%n(2)+2,0:self%dom%n(3)+2,self%ntype))

    allocate(self%phi_avg_xy(0:self%dom%n(1)+2,0:self%dom%n(2)+2))
    allocate(self%phi_avg_xz(0:self%dom%n(1)+2,0:self%dom%n(3)+2))
    allocate(self%phi_avg_yz(0:self%dom%n(2)+2,0:self%dom%n(3)+2))

    self%np_avg_xy  = 0.0_real64
    self%np_avg_xz  = 0.0_real64
    self%np_avg_yz  = 0.0_real64
    self%phi_avg_xy = 0.0_real64
    self%phi_avg_xz = 0.0_real64
    self%phi_avg_yz = 0.0_real64
    self%cnt_avg    = 0_int32

    if (self%cfg%flag_avg3D == 1 .and. self%mpi_rank == 0) then
      allocate(self%avg3D(5, 1:self%dom%n(1)+1, 1:self%dom%n(2)+1, 1:self%dom%n(3)+1))
      self%avg3D = 0.0_real64
    end if

    allocate(self%data_pavg_xy(5, 0:self%dom%n(1)+2, 0:self%dom%n(2)+2, self%ntype))
    allocate(self%data_pavg_xz(5, 0:self%dom%n(1)+2, 0:self%dom%n(3)+2, self%ntype))
    allocate(self%data_pavg_yz(5, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2, self%ntype))

    self%data_pavg_xy = 0.0_real64
    self%data_pavg_xz = 0.0_real64
    self%data_pavg_yz = 0.0_real64

    allocate(self%sour_avg_xy(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, self%ntype, self%nproc))
    allocate(self%sour_avg_xz(0:self%dom%n(1)+2, 0:self%dom%n(3)+2, self%ntype, self%nproc))
    allocate(self%sour_avg_yz(0:self%dom%n(2)+2, 0:self%dom%n(3)+2, self%ntype, self%nproc))
    allocate(self%sink_avg_xy(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, self%ntype, self%nproc))
    allocate(self%sink_avg_xz(0:self%dom%n(1)+2, 0:self%dom%n(3)+2, self%ntype, self%nproc))
    allocate(self%sink_avg_yz(0:self%dom%n(2)+2, 0:self%dom%n(3)+2, self%ntype, self%nproc))

    self%sour_avg_xy = 0.0_real64
    self%sour_avg_xz = 0.0_real64
    self%sour_avg_yz = 0.0_real64
    self%sink_avg_xy = 0.0_real64
    self%sink_avg_xz = 0.0_real64
    self%sink_avg_yz = 0.0_real64

  end subroutine init


  subroutine init_chemistry(self)
    class(State), intent(inout) :: self

    call self%chem%init(npart, self%rxn%ncol_mx, self%rxn%npt_mx)
    self%chem%ngas = self%cfg%ngas
    call self%rxn%load(self%chem, trim(self%cfg%rname), self%mpi_rank)

    if (self%mpi_rank == 0) then
      write(*,'(i0,a,i0,a,i0,a)') self%rxn%ntype, ' species:  (', &
           self%rxn%ntype-self%rxn%n_neu, ' charged and ', self%rxn%n_neu, ' neutral)'
      write(*,'(*(a,1x))') self%chem%pname
      write(*,'(a,i0)') 'Total number of reactions: ', self%chem%ncol
      write(*,*) '  '
    end if
  end subroutine init_chemistry


  subroutine init_domain_and_boundary(self)
    class(State), intent(inout) :: self
    integer :: ierr
    real(real64), allocatable :: phi_rt(:,:,:)
    real(real64) :: lmax, gmax

    call build_boundary(self%dom, self%cfg, self%fld, self%mpi_rank)
    call self%magField%build_from_cfg(self%cfg, self%dom, self%mpi_rank)
    call self%magField%write_macho_planes('./Output/Output_2D', 1, self%mpi_rank)

    call self%pdec%init(int(self%dom%n(1),int32), int(self%dom%n(2),int32), int(self%dom%n(3),int32), self%comm)
    call self%pdec%scatter_from_global(self%fld%phi, self%dom%bcnd)

    if (DEBUG_INIT_CHECKS) then
      call checkpoint_poisson_decomp(self%mpi_rank, self%comm, &
                                     int(self%dom%n(3),int32), self%pdec%k0, self%pdec%m, &
                                     self%pdec%phi_dom, self%pdec%bcnd_dom)
    end if

    call build_kq(self%dom%bcnd, self%fld%kq)

    if (DEBUG_INIT_CHECKS) then
      call checkpoint_kq(self%fld%kq, self%comm, self%mpi_rank, 'after build_kq', self%pdec)

      allocate(phi_rt(0:self%dom%n(1)+2, 0:self%dom%n(2)+2, 0:self%dom%n(3)+2))
      phi_rt = -999.0_real64

      call self%pdec%gather_phi_to_global(phi_rt)

      lmax = maxval(abs(phi_rt - self%fld%phi))
      call MPI_Allreduce(lmax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)

      deallocate(phi_rt)
    end if

    if (self%mpi_rank == 0) then
      write(*,'(a)') '  '
      write(*,'(a)') 'Building boundary...'
      write(*,'(a,3(i0,1x))')       'n = ', self%dom%n(1), self%dom%n(2), self%dom%n(3)
      write(*,'(a,3(1p,e12.4,1x))') 'h(m) = ', self%dom%h(1), self%dom%h(2), self%dom%h(3)
      write(*,'(a,3(1p,e12.4,1x))') 'box(m) = ', self%dom%xmax, self%dom%ymax, self%dom%zmax
      write(*,'(a,1p,e12.4)')       'Sg(m^2) = ', self%dom%Sg
      write(*,'(a,4(i0,1x))')       'flags(pbc,pbcz,nmn,die) = ', &
           self%dom%flag_pbc, self%dom%flag_pbcz, self%dom%flag_nmn, self%dom%flag_die
      write(*,*) '  '
    end if
  end subroutine init_domain_and_boundary


  subroutine init_particles(self)
    class(State), intent(inout) :: self
    integer(int32) :: ntype_trk
    integer(int32) :: ptype, iproc
    real(real64)   :: E0_RF_restart

    ntype_trk = self%ntype

    if (allocated(self%part)) then
      call self%finalize_particles_only()
    end if
    allocate(self%part(self%ntype, self%nproc))

    call self%fld%allocate_species_density(self%dom, self%ntype)

    self%np_thread = 0.0_real64

    if (self%cfg%flag_restart > 0_int32) then
      call restart_particles_modular( &
            cfg       = self%cfg, &
            mpi_rank  = self%mpi_rank, &
            mpi_size  = self%mpi_size, &
            n         = int(self%dom%n, int32), &
            h         = self%dom%h, &
            kq        = self%fld%kq, &
            Nm        = self%params%Nm, &
            ntype     = ntype_trk, &
            nproc     = self%nproc, &
            tag_neg   = int(self%chem%tag_neg, int32), &
            iseed     = self%params%iseed, &
            part      = self%part, &
            np_thread = self%np_thread, &
            time_out  = self%time, &
            E0_RF_out = E0_RF_restart )
      ! Applied after the call returns, not passed as an alias of
      ! self%cfg%E0_RF directly - cfg (intent in, whole struct) and
      ! E0_RF_out (intent out) would otherwise be overlapping actual
      ! arguments, which is non-conforming per the Fortran standard's
      ! argument-association rules.
      self%cfg%E0_RF = E0_RF_restart
    else
      call load_particles_modular( &
            cfg       = self%cfg, &
            mpi_rank  = self%mpi_rank, &
            mpi_size  = self%mpi_size, &
            n         = int(self%dom%n, int32), &
            h         = self%dom%h, &
            bcnd      = self%dom%bcnd, &
            kq        = self%fld%kq, &
            vt0       = self%params%vt0, &
            Nm        = self%params%Nm, &
            ni0       = self%chem%ni0(1:self%ntype), &
            iseed     = self%params%iseed, &
            ntype_trk = ntype_trk, &
            part      = self%part, &
            np_thread = self%np_thread )
    end if
    call reduce_species_density( &
         n         = int(self%dom%n, int32), &
         bcnd      = self%dom%bcnd, &
         np_thread = self%np_thread, &
         ntype     = int(self%ntype), &
         nproc     = int(self%nproc), &
         mpi_comm  = self%comm, &
         np_red    = self%fld%np )

    ! Preallocate extra particle capacity once, before the PIC loop.
    ! This allows collision-created particles to be appended without
    ! reallocating inside the Monte Carlo collision kernel.
    do iproc = 1_int32, self%nproc
      do ptype = 1_int32, self%ntype
        call self%part(ptype,iproc)%ensure_capacity( &
          self%part(ptype,iproc)%n + max(100000_int32, self%part(ptype,iproc)%n / 2_int32) )
      end do
    end do

  end subroutine init_particles


  subroutine sort_particles_local(self)
    class(State), intent(inout) :: self
    integer(int32) :: ptype, iproc
    logical :: ok_sorted, ok_cells

    if (.not. allocated(self%part)) return

    ! Global per-species sort + redistribution: each iproc ends up with a
    ! contiguous run of cells (a z-slab), as in legacy part_sorting. This
    ! replaced a per-iproc sort that left every thread's particles spread
    ! over the whole domain, so every thread's field gathers and deposit
    ! swept the full E/B/np grids - see redistribute_species_by_cell.
    do ptype = 1, self%ntype
      call redistribute_species_by_cell(self%part(ptype,:), int(self%dom%n, int32), &
                                        self%dom%h, ptype)
    end do

    if (VALIDATE_SORTING) then
      do ptype = 1, self%ntype
        do iproc = 1, self%nproc
          if (.not. allocated(self%part(ptype,iproc)%pv)) cycle
          if (self%part(ptype,iproc)%n <= 1_int32) cycle

          ok_sorted = check_particles_are_sorted(self%part(ptype,iproc))
          if (.not. ok_sorted) then
            write(*,'(a,3(i0,1x))') 'Sorting failed on rank, ptype, iproc = ', &
                 self%mpi_rank, ptype, iproc
            error stop 'mod_state%sort_particles_local: particle sorting failed'
          end if

          ok_cells = check_cell_indexing(self%part(ptype,iproc), int(self%dom%n, int32))
          if (.not. ok_cells) then
            write(*,'(a,3(i0,1x))') 'Cell indexing failed on rank, ptype, iproc = ', &
                 self%mpi_rank, ptype, iproc
            error stop 'mod_state%sort_particles_local: cell indexing failed'
          end if
        end do
      end do
    end if
  end subroutine sort_particles_local


  subroutine apply_electron_heating_local(self, vt)
    class(State), intent(inout) :: self
    real(real64), intent(in)    :: vt
    integer(int32) :: iproc
    real(real64)   :: Ploss_before

    if (.not. allocated(self%part)) return
    if (self%ntype < 1) return
    if (.not. allocated(self%P_loss)) return
    if (self%cfg%flag_circxh == 1) self%cfg%R_ahp = self%cfg%yr_pow - self%dom%ymax/2.0_real64

    Ploss_before = sum(self%P_loss(2,1,:))

    !$omp parallel do private(iproc) schedule(static)
    do iproc = 1, self%nproc
      if (.not. allocated(self%part(1,iproc)%pv)) cycle
      if (self%part(1,iproc)%n <= 0_int32) cycle

      call apply_electron_heating( &
          part           = self%part(1,iproc), &
          h              = self%dom%h, &
          vt             = vt, &
          iseed          = self%params%iseed(iproc), &
          nudt           = self%params%nudt, &
          xl_pow         = self%cfg%xl_pow, &
          xr_pow         = self%cfg%xr_pow, &
          flag_circxh    = self%cfg%flag_circxh, &
          flag_ahp       = self%cfg%flag_ahp, &
          R_ahp          = self%cfg%R_ahp, &
          delta_heat     = self%cfg%delta_heat, &
          ymax           = self%dom%ymax, &
          zmax           = self%dom%zmax, &
          Nm_e           = self%params%Nm(1), &
          mass_e         = self%chem%mass(1), &
          p_loss_heating = self%P_loss(2,1,iproc) )
    end do
    !$omp end parallel do

    !if (self%mpi_rank == 0) &
     ! write(*,'(a,es12.4)') ' HEAT_MOD_dEk=', sum(self%P_loss(2,1,:)) - Ploss_before

  end subroutine apply_electron_heating_local


  subroutine apply_dielectric_bc_to_phi(self)
    use iso_fortran_env, only: int32, real64
    class(State), intent(inout) :: self

    integer(int32) :: ix, iy, iz, ptype, iproc
    integer(int32) :: igrid
    real(real64)   :: Ca
    real(real64)   :: qsum_xz, qsum_yz_1, qsum_yz_2

    if (self%dom%flag_die /= 1) return
    if (.not. allocated(self%sum_q_xz)) return
    if (.not. allocated(self%sum_q_yz)) return

    Ca = self%params%Ca_cell

    do iz = 0, self%dom%n(3)+2
      do iy = 0, self%dom%n(2)+2
        do ix = 0, self%dom%n(1)+2

          igrid = self%dom%bcnd(ix,iy,iz)
          if (igrid <= 0) cycle

          do ptype = 1, self%ntype

            if (self%dom%dtype(igrid) == 2) then
              qsum_xz = 0.0_real64
              do iproc = 1, self%nproc
                qsum_xz = qsum_xz + self%sum_q_xz(ix,iz,ptype,iproc)
              end do
              self%fld%phi(ix,iy,iz) = self%fld%phi(ix,iy,iz) + qsum_xz / (2.0_real64 * Ca)
            end if

            if (self%dom%dtype(igrid) == 3) then
              qsum_yz_1 = 0.0_real64
              do iproc = 1, self%nproc
                qsum_yz_1 = qsum_yz_1 + self%sum_q_yz(1,iy,iz,ptype,iproc)
              end do
              self%fld%phi(ix,iy,iz) = self%fld%phi(ix,iy,iz) + qsum_yz_1 / Ca
            end if

            if (self%dom%dtype(igrid) == 4) then
              qsum_yz_2 = 0.0_real64
              do iproc = 1, self%nproc
                qsum_yz_2 = qsum_yz_2 + self%sum_q_yz(2,iy,iz,ptype,iproc)
              end do
              self%fld%phi(ix,iy,iz) = self%fld%phi(ix,iy,iz) + qsum_yz_2 / Ca
            end if

          end do
        end do
      end do
    end do
  end subroutine apply_dielectric_bc_to_phi


  subroutine advance_particles_local(self, istep, t_move_out, t_bc_out, t_deposit_out, &
                                      t_move_min_out, t_move_avg_out, &
                                      t_deposit_min_out, t_deposit_avg_out, &
                                      n_boris_used_out, n_boris_total_out)
    ! Fused push + boundary-condition/SEE + charge-deposition pass.
    !
    ! Not a type-bound procedure (called as advance_particles_local(self%state, ...)
    ! from mod_simulation.f90, not self%state%advance_particles_local(...)) so that
    ! self below can be type(State) instead of class(State) - a type-bound
    ! procedure's passed-object argument is required by the standard to be
    ! polymorphic. State is never extended, so the polymorphism bought nothing
    ! but prevented the compiler from proving self%np_thread(...) below is a
    ! plain contiguous array; profiling showed the per-step clear of that
    ! array (__intel_avx_rep_memset) as the single largest self-time cost in
    ! the whole run, disproportionately larger than the equivalent clear in
    ! legacy's plain (non-derived-type) np(...) array.
    !
    ! Parallelizes over iproc only (no collapse with ptype), one thread
    ! owning every species of its own iproc for the whole pass - matching
    ! legacy's single fused OMP region in Src/main.f90 instead of three
    ! separate parallel regions (mover, BC, deposit), each paying its own
    ! fork/join and re-walking the particle arrays cold. This also removes
    ! a latent cross-thread race the old collapse(2) loops had: SEE writes
    ! into part(1,iproc) and draws from iseed(iproc) (intent(inout)) keyed
    ! only by iproc, so two threads processing different ptypes of the
    ! SAME iproc under collapse(2) static chunking could hit those
    ! concurrently once ntype>1.
    !
    ! Deposition runs in a second ptype loop, after every species has been
    ! pushed and BC'd for this iproc - not folded into the same ptype pass
    ! - because a secondary electron created while processing another
    ! species' BC this step must still land in this step's electron
    ! deposit (mirrors legacy's separate push-then-deposit loop order).
    !
    ! t_move_out/t_bc_out/t_deposit_out: diagnostic wall-clock split of the
    ! fused region, since the outer MPI_Wtime() bracket around the whole
    ! call (mod_simulation.f90) can no longer tell push apart from deposit
    ! now that they share one parallel region. Each thread times its own
    ! move_and_bc_*/apply_particle_bc call (summed per iproc) and its own
    ! deposit loop with omp_get_wtime(), into per-iproc slots it alone
    ! writes to; the max across iproc (the slowest thread sets the region's
    ! wall time) is returned. Purely diagnostic - does not change the
    ! computation. t_bc_out is always 0: push and BC used to be two
    ! separate calls timed separately, but move_and_bc_electrostatic/boris
    ! (mod_particleMover.f90) now fuse them into one per-particle loop, so
    ! there is nothing left to time apart from t_move_out - kept as a
    ! separate output (rather than removed) so mod_simulation.f90's
    ! t_mover_bc field and print line don't need touching again.
    !
    ! t_move_min_out/t_move_avg_out and t_deposit_min_out/t_deposit_avg_out:
    ! min and mean across iproc of the same per-thread timings that
    ! t_move_out/t_deposit_out already reduce with maxval. Added to check
    ! OMP load imbalance across iproc directly (max/avg ratio close to 1
    ! means threads are evenly loaded; a large ratio means a few iproc are
    ! doing most of the work and the mover cost is imbalance, not raw
    ! per-particle cost) instead of guessing from the max alone.
    type(State), intent(inout) :: self
    integer(int32), intent(in) :: istep
    real(real64), intent(out) :: t_move_out, t_bc_out, t_deposit_out
    real(real64), intent(out) :: t_move_min_out, t_move_avg_out
    real(real64), intent(out) :: t_deposit_min_out, t_deposit_avg_out
    ! TEMPORARY diagnostic - see move_and_bc_boris's header comment
    ! (mod_particleMover.f90) and mod_simulation.f90's diagnostic print.
    integer(int32), intent(out) :: n_boris_used_out, n_boris_total_out

    integer(int32)  :: ptype, iproc
    integer(int32)  :: tag_neg_local, ispec
    real(real64)    :: q_species, m_species, dt_local, qmacro
    real(real64)    :: tw0, tw1
    logical         :: use_boris
    logical         :: use_energy_conserving
    logical         :: use_fast_electrostatic
    logical         :: use_fast_boris
    type(SeeParams) :: see
    real(real64), allocatable :: t_move_thread(:), t_bc_thread(:), t_deposit_thread(:)
    integer(int32), allocatable :: n_boris_used_thread(:), n_boris_total_thread(:)
    integer(int32)  :: n_boris_used_call, n_boris_total_call

    ! RF antenna heating: quantities constant for the whole step, computed
    ! once here rather than per-particle inside the movers (mirrors legacy
    ! computing omega_RF once per part_mover call).
    real(real64)   :: time_rf, omega_RF_local
    integer(int32) :: ixl_pow_rf, ixr_pow_rf
    real(real64)   :: x0_rf

    ! Electron heating-region tally (Nh/sum_dEk), done inside the ptype==1
    ! deposit below instead of a separate pass afterwards - see
    ! HeatRegionTally in mod_chargeDeposition.f90.
    logical               :: do_heat_tally
    type(HeatRegionTally) :: heat_tpl, heat_loc

    t_move_out        = 0.0_real64
    t_bc_out          = 0.0_real64
    t_deposit_out     = 0.0_real64
    t_move_min_out    = 0.0_real64
    t_move_avg_out    = 0.0_real64
    t_deposit_min_out = 0.0_real64
    t_deposit_avg_out = 0.0_real64
    n_boris_used_out  = 0_int32
    n_boris_total_out = 0_int32

    if (.not. allocated(self%part)) return

    allocate(t_move_thread(self%nproc), t_bc_thread(self%nproc), t_deposit_thread(self%nproc))
    t_move_thread    = 0.0_real64
    t_bc_thread      = 0.0_real64
    t_deposit_thread = 0.0_real64

    allocate(n_boris_used_thread(self%nproc), n_boris_total_thread(self%nproc))
    n_boris_used_thread  = 0_int32
    n_boris_total_thread = 0_int32

    dt_local = self%params%dt
    use_energy_conserving = (trim(self%cfg%push_scheme) == 'energy')
    ! Fast electrostatic path (mod_particleMover.f90): skips both this and
    ! RF-antenna heating entirely rather than branching on them per-particle
    ! - see move_and_bc_electrostatic_fast's header comment.
    use_fast_electrostatic = (.not. use_energy_conserving) .and. &
                              (self%cfg%flag_RFant /= 1_int32)
    ! Same condition, Boris side - see move_and_bc_boris_fast's header
    ! comment.
    use_fast_boris = (.not. use_energy_conserving) .and. &
                      (self%cfg%flag_RFant /= 1_int32)

    time_rf        = 0.0_real64
    omega_RF_local = 0.0_real64
    ixl_pow_rf     = 0_int32
    ixr_pow_rf     = 0_int32
    x0_rf          = 0.0_real64
    if (self%cfg%flag_RFant == 1_int32) then
      time_rf        = real(istep, real64) * dt_local
      omega_RF_local = 2.0_real64 * pi * self%cfg%f0_RF
      ! Clamp to the domain box - see update_rf_convergence's header comment.
      ixl_pow_rf     = int(max(0.0_real64, self%cfg%xl_pow) / self%dom%h(1), int32) + 1_int32
      ixr_pow_rf     = int(min(self%dom%xmax, self%cfg%xr_pow) / self%dom%h(1), int32) + 1_int32
      ! Planar-coil geometry only: window position, same clamp as ixl_pow_rf.
      x0_rf          = max(0.0_real64, self%cfg%xl_pow)
    end if

    tag_neg_local = -1_int32
    do ispec = 1_int32, self%ntype
      if (self%chem%charge(ispec) < 0.0_real64 .and. ispec /= 1_int32) then
        tag_neg_local = ispec
        exit
      end if
    end do

    see%gam_sec   = real(self%cfg%gam_sec, real64)
    see%igrid_sec = int(self%cfg%igrid_sec, int32)
    see%zg_sec    = self%dom%zg_sec
    see%Nm_e      = self%params%Nm(1)
    see%mass_e    = self%chem%mass(1)
    see%vt_sec    = 0.0_real64
    if (self%cfg%THm /= 0.0_real64 .and. self%chem%mass(1) /= 0.0_real64) then
      see%vt_sec = sqrt(2.0_real64 * qe * abs(self%cfg%THm) / abs(self%chem%mass(1)))
    end if

    do_heat_tally = self%cfg%Pabs > 0.0_real64 .and. self%cfg%flag_heat == 1 .and. &
                    allocated(self%Nh) .and. allocated(self%sum_dEk)
    if (do_heat_tally) then
      heat_tpl%ixl         = int(self%cfg%xl_pow / self%dom%h(1), int32) + 1_int32
      heat_tpl%ixr         = int(self%cfg%xr_pow / self%dom%h(1), int32) + 1_int32
      heat_tpl%flag_circxh = int(self%cfg%flag_circxh, int32)
      heat_tpl%flag_ahp    = int(self%cfg%flag_ahp, int32)
      heat_tpl%R2          = self%cfg%R_ahp**2
      heat_tpl%yc          = self%dom%ymax / 2.0_real64
      heat_tpl%zc          = self%dom%zmax / 2.0_real64
      heat_tpl%ek_coef     = 0.5_real64 * self%params%Nm(1) * self%chem%mass(1)
      heat_tpl%R           = self%cfg%R_ahp
      heat_tpl%delta       = self%cfg%delta_heat
      heat_tpl%Nh          = 0_int32
      heat_tpl%sum_dEk     = 0.0_real64
      heat_tpl%sum_w       = 0.0_real64
    end if

    !$omp parallel do private(iproc,ptype,q_species,m_species,use_boris,qmacro,tw0,tw1, &
    !$omp&                    n_boris_used_call,n_boris_total_call,heat_loc) schedule(static)
    do iproc = 1, self%nproc

      do ptype = 1, self%ntype
        if (.not. allocated(self%part(ptype,iproc)%pv)) cycle
        if (self%part(ptype,iproc)%n <= 0_int32) cycle

        q_species = self%chem%charge(ptype)
        m_species = self%chem%mass(ptype)
        qmacro    = self%params%Nm(ptype) * q_species

        tw0 = omp_get_wtime()
        if (m_species > 0.0_real64) then
          ! Push+BC/SEE fused into one per-particle loop (move_and_bc_*,
          ! mod_particleMover.f90) instead of two full passes over the
          ! same position/velocity arrays - see that module for why
          ! (measured ~76ms push + ~56ms BC as separate passes).
          use_boris = (self%cfg%flag_B == 1) .and. &
                      .not. (self%cfg%flag_B_pos == 1 .and. q_species > 0.0_real64)

          if (use_boris) then
            ! ptype==1 (electrons) is split out rather than always passing
            ! part_electrons=self%part(1,iproc): for that species it's the
            ! same actual storage as `part` above, and do_see is provably
            ! false anyway since electrons carry q<0 - see mod_particleMover.f90's
            ! header comment on move_and_bc_electrostatic for why the
            ! compiler can't use that runtime fact on its own to stop
            ! assuming the two INTENT(INOUT) ParticleSets might alias.
            !
            ! use_fast_boris picks move_and_bc_boris_fast (RF-antenna and
            ! energy-conserving branches removed entirely, not just runtime-
            ! skipped) over the general move_and_bc_boris whenever both are
            ! off for this run - see move_and_bc_boris_fast's header comment.
            if (use_fast_boris) then
              if (ptype == 1_int32) then
                call move_and_bc_boris_fast( &
                    part           = self%part(ptype,iproc), &
                    n              = int(self%dom%n, int32), &
                    h              = self%dom%h, &
                    E              = self%fld%E, &
                    n_B            = self%magField%n_B, &
                    h_B            = self%magField%h_B, &
                    Bi             = self%magField%Bi, &
                    q              = q_species, &
                    m              = m_species, &
                    dt             = dt_local, &
                    bcnd           = self%dom%bcnd, &
                    wall_cell      = self%dom%wall_cell, &
                    xmax           = self%dom%xmax, &
                    ymax           = self%dom%ymax, &
                    zmax           = self%dom%zmax, &
                    flag_pbc       = int(self%dom%flag_pbc, int32), &
                    flag_nmn       = int(self%dom%flag_nmn, int32), &
                    ptype          = ptype, &
                    tag_neg        = tag_neg_local, &
                    flag_die       = int(self%dom%flag_die, int32), &
                    dtype          = self%dom%dtype, &
                    qmacro         = qmacro, &
                    sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                    sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                    p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                    P_loss_wall    = self%P_loss(1,ptype,iproc), &
                    Nm_species     = self%params%Nm(ptype), &
                    see            = see, &
                    iseed          = self%params%iseed(iproc), &
                    P_loss_see     = self%P_loss(4,1,iproc), &
                    mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                    mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                    P_loss_coll    = self%P_loss(3,ptype,iproc), &
                    mom_loss_coll  = self%mom_loss(:,3,ptype,iproc), &
                    n_boris_used   = n_boris_used_call, &
                    n_boris_total  = n_boris_total_call )
                n_boris_used_thread(iproc)  = n_boris_used_thread(iproc)  + n_boris_used_call
                n_boris_total_thread(iproc) = n_boris_total_thread(iproc) + n_boris_total_call
              else
                call move_and_bc_boris_fast( &
                    part           = self%part(ptype,iproc), &
                    n              = int(self%dom%n, int32), &
                    h              = self%dom%h, &
                    E              = self%fld%E, &
                    n_B            = self%magField%n_B, &
                    h_B            = self%magField%h_B, &
                    Bi             = self%magField%Bi, &
                    q              = q_species, &
                    m              = m_species, &
                    dt             = dt_local, &
                    bcnd           = self%dom%bcnd, &
                    wall_cell      = self%dom%wall_cell, &
                    xmax           = self%dom%xmax, &
                    ymax           = self%dom%ymax, &
                    zmax           = self%dom%zmax, &
                    flag_pbc       = int(self%dom%flag_pbc, int32), &
                    flag_nmn       = int(self%dom%flag_nmn, int32), &
                    ptype          = ptype, &
                    tag_neg        = tag_neg_local, &
                    flag_die       = int(self%dom%flag_die, int32), &
                    dtype          = self%dom%dtype, &
                    qmacro         = qmacro, &
                    sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                    sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                    p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                    P_loss_wall    = self%P_loss(1,ptype,iproc), &
                    Nm_species     = self%params%Nm(ptype), &
                    see            = see, &
                    part_electrons = self%part(1,iproc), &
                    iseed          = self%params%iseed(iproc), &
                    P_loss_see     = self%P_loss(4,1,iproc), &
                    mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                    mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                    P_loss_coll    = self%P_loss(3,ptype,iproc), &
                    mom_loss_coll  = self%mom_loss(:,3,ptype,iproc), &
                    n_boris_used   = n_boris_used_call, &
                    n_boris_total  = n_boris_total_call )
                n_boris_used_thread(iproc)  = n_boris_used_thread(iproc)  + n_boris_used_call
                n_boris_total_thread(iproc) = n_boris_total_thread(iproc) + n_boris_total_call
              end if
            else
              if (ptype == 1_int32) then
                call move_and_bc_boris( &
                    part           = self%part(ptype,iproc), &
                    n              = int(self%dom%n, int32), &
                    h              = self%dom%h, &
                    E              = self%fld%E, &
                    n_B            = self%magField%n_B, &
                    h_B            = self%magField%h_B, &
                    Bi             = self%magField%Bi, &
                    q              = q_species, &
                    m              = m_species, &
                    dt             = dt_local, &
                    use_energy_conserving = use_energy_conserving, &
                    bcnd           = self%dom%bcnd, &
                    wall_cell      = self%dom%wall_cell, &
                    xmax           = self%dom%xmax, &
                    ymax           = self%dom%ymax, &
                    zmax           = self%dom%zmax, &
                    flag_pbc       = int(self%dom%flag_pbc, int32), &
                    flag_nmn       = int(self%dom%flag_nmn, int32), &
                    ptype          = ptype, &
                    tag_neg        = tag_neg_local, &
                    flag_die       = int(self%dom%flag_die, int32), &
                    dtype          = self%dom%dtype, &
                    qmacro         = qmacro, &
                    sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                    sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                    p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                    P_loss_wall    = self%P_loss(1,ptype,iproc), &
                    Nm_species     = self%params%Nm(ptype), &
                    see            = see, &
                    iseed          = self%params%iseed(iproc), &
                    P_loss_see     = self%P_loss(4,1,iproc), &
                    flag_RFant     = self%cfg%flag_RFant, &
                    ixl_pow        = ixl_pow_rf, &
                    ixr_pow        = ixr_pow_rf, &
                    R_ahp          = self%cfg%R_ahp, &
                    gams           = self%gams, &
                    E0_RF          = self%cfg%E0_RF, &
                    omega_RF       = omega_RF_local, &
                    time           = time_rf, &
                    P_RF_local     = self%P_RF(ptype,iproc), &
                    flag_planar_ant = self%cfg%flag_planar_ant, &
                    x0             = x0_rf, &
                    mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                    mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                    mom_RF_local   = self%mom_RF(:,ptype,iproc), &
                    P_loss_coll    = self%P_loss(3,ptype,iproc), &
                    mom_loss_coll  = self%mom_loss(:,3,ptype,iproc), &
                    n_boris_used   = n_boris_used_call, &
                    n_boris_total  = n_boris_total_call )
                n_boris_used_thread(iproc)  = n_boris_used_thread(iproc)  + n_boris_used_call
                n_boris_total_thread(iproc) = n_boris_total_thread(iproc) + n_boris_total_call
              else
                call move_and_bc_boris( &
                    part           = self%part(ptype,iproc), &
                    n              = int(self%dom%n, int32), &
                    h              = self%dom%h, &
                    E              = self%fld%E, &
                    n_B            = self%magField%n_B, &
                    h_B            = self%magField%h_B, &
                    Bi             = self%magField%Bi, &
                    q              = q_species, &
                    m              = m_species, &
                    dt             = dt_local, &
                    use_energy_conserving = use_energy_conserving, &
                    bcnd           = self%dom%bcnd, &
                    wall_cell      = self%dom%wall_cell, &
                    xmax           = self%dom%xmax, &
                    ymax           = self%dom%ymax, &
                    zmax           = self%dom%zmax, &
                    flag_pbc       = int(self%dom%flag_pbc, int32), &
                    flag_nmn       = int(self%dom%flag_nmn, int32), &
                    ptype          = ptype, &
                    tag_neg        = tag_neg_local, &
                    flag_die       = int(self%dom%flag_die, int32), &
                    dtype          = self%dom%dtype, &
                    qmacro         = qmacro, &
                    sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                    sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                    p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                    P_loss_wall    = self%P_loss(1,ptype,iproc), &
                    Nm_species     = self%params%Nm(ptype), &
                    see            = see, &
                    part_electrons = self%part(1,iproc), &
                    iseed          = self%params%iseed(iproc), &
                    P_loss_see     = self%P_loss(4,1,iproc), &
                    flag_RFant     = self%cfg%flag_RFant, &
                    ixl_pow        = ixl_pow_rf, &
                    ixr_pow        = ixr_pow_rf, &
                    R_ahp          = self%cfg%R_ahp, &
                    gams           = self%gams, &
                    E0_RF          = self%cfg%E0_RF, &
                    omega_RF       = omega_RF_local, &
                    time           = time_rf, &
                    P_RF_local     = self%P_RF(ptype,iproc), &
                    flag_planar_ant = self%cfg%flag_planar_ant, &
                    x0             = x0_rf, &
                    mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                    mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                    mom_RF_local   = self%mom_RF(:,ptype,iproc), &
                    P_loss_coll    = self%P_loss(3,ptype,iproc), &
                    mom_loss_coll  = self%mom_loss(:,3,ptype,iproc), &
                    n_boris_used   = n_boris_used_call, &
                    n_boris_total  = n_boris_total_call )
                n_boris_used_thread(iproc)  = n_boris_used_thread(iproc)  + n_boris_used_call
                n_boris_total_thread(iproc) = n_boris_total_thread(iproc) + n_boris_total_call
              end if
            end if
          else if (use_fast_electrostatic) then
            ! ptype==1 split out - see the move_and_bc_boris call above.
            if (ptype == 1_int32) then
              call move_and_bc_electrostatic_fast( &
                  part           = self%part(ptype,iproc), &
                  n              = int(self%dom%n, int32), &
                  h              = self%dom%h, &
                  E              = self%fld%E, &
                  q              = q_species, &
                  m              = m_species, &
                  dt             = dt_local, &
                  bcnd           = self%dom%bcnd, &
                  wall_cell      = self%dom%wall_cell, &
                  xmax           = self%dom%xmax, &
                  ymax           = self%dom%ymax, &
                  zmax           = self%dom%zmax, &
                  flag_pbc       = int(self%dom%flag_pbc, int32), &
                  flag_nmn       = int(self%dom%flag_nmn, int32), &
                  ptype          = ptype, &
                  tag_neg        = tag_neg_local, &
                  flag_die       = int(self%dom%flag_die, int32), &
                  dtype          = self%dom%dtype, &
                  qmacro         = qmacro, &
                  sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                  sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                  p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                  P_loss_wall    = self%P_loss(1,ptype,iproc), &
                  Nm_species     = self%params%Nm(ptype), &
                  see            = see, &
                  iseed          = self%params%iseed(iproc), &
                  P_loss_see     = self%P_loss(4,1,iproc), &
                  mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                  mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                  P_loss_coll    = self%P_loss(3,ptype,iproc), &
                  mom_loss_coll  = self%mom_loss(:,3,ptype,iproc) )
            else
              call move_and_bc_electrostatic_fast( &
                  part           = self%part(ptype,iproc), &
                  n              = int(self%dom%n, int32), &
                  h              = self%dom%h, &
                  E              = self%fld%E, &
                  q              = q_species, &
                  m              = m_species, &
                  dt             = dt_local, &
                  bcnd           = self%dom%bcnd, &
                  wall_cell      = self%dom%wall_cell, &
                  xmax           = self%dom%xmax, &
                  ymax           = self%dom%ymax, &
                  zmax           = self%dom%zmax, &
                  flag_pbc       = int(self%dom%flag_pbc, int32), &
                  flag_nmn       = int(self%dom%flag_nmn, int32), &
                  ptype          = ptype, &
                  tag_neg        = tag_neg_local, &
                  flag_die       = int(self%dom%flag_die, int32), &
                  dtype          = self%dom%dtype, &
                  qmacro         = qmacro, &
                  sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                  sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                  p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                  P_loss_wall    = self%P_loss(1,ptype,iproc), &
                  Nm_species     = self%params%Nm(ptype), &
                  see            = see, &
                  part_electrons = self%part(1,iproc), &
                  iseed          = self%params%iseed(iproc), &
                  P_loss_see     = self%P_loss(4,1,iproc), &
                  mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                  mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                  P_loss_coll    = self%P_loss(3,ptype,iproc), &
                  mom_loss_coll  = self%mom_loss(:,3,ptype,iproc) )
            end if
          else
            ! ptype==1 split out - see the move_and_bc_boris call above.
            if (ptype == 1_int32) then
              call move_and_bc_electrostatic( &
                  part           = self%part(ptype,iproc), &
                  n              = int(self%dom%n, int32), &
                  h              = self%dom%h, &
                  E              = self%fld%E, &
                  q              = q_species, &
                  m              = m_species, &
                  dt             = dt_local, &
                  use_energy_conserving = use_energy_conserving, &
                  bcnd           = self%dom%bcnd, &
                  wall_cell      = self%dom%wall_cell, &
                  xmax           = self%dom%xmax, &
                  ymax           = self%dom%ymax, &
                  zmax           = self%dom%zmax, &
                  flag_pbc       = int(self%dom%flag_pbc, int32), &
                  flag_nmn       = int(self%dom%flag_nmn, int32), &
                  ptype          = ptype, &
                  tag_neg        = tag_neg_local, &
                  flag_die       = int(self%dom%flag_die, int32), &
                  dtype          = self%dom%dtype, &
                  qmacro         = qmacro, &
                  sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                  sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                  p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                  P_loss_wall    = self%P_loss(1,ptype,iproc), &
                  Nm_species     = self%params%Nm(ptype), &
                  see            = see, &
                  iseed          = self%params%iseed(iproc), &
                  P_loss_see     = self%P_loss(4,1,iproc), &
                  flag_RFant     = self%cfg%flag_RFant, &
                  ixl_pow        = ixl_pow_rf, &
                  ixr_pow        = ixr_pow_rf, &
                  R_ahp          = self%cfg%R_ahp, &
                  gams           = self%gams, &
                  E0_RF          = self%cfg%E0_RF, &
                  omega_RF       = omega_RF_local, &
                  time           = time_rf, &
                  P_RF_local     = self%P_RF(ptype,iproc), &
                  flag_planar_ant = self%cfg%flag_planar_ant, &
                  x0             = x0_rf, &
                  mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                  mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                  mom_RF_local   = self%mom_RF(:,ptype,iproc), &
                  P_loss_coll    = self%P_loss(3,ptype,iproc), &
                  mom_loss_coll  = self%mom_loss(:,3,ptype,iproc) )
            else
              call move_and_bc_electrostatic( &
                  part           = self%part(ptype,iproc), &
                  n              = int(self%dom%n, int32), &
                  h              = self%dom%h, &
                  E              = self%fld%E, &
                  q              = q_species, &
                  m              = m_species, &
                  dt             = dt_local, &
                  use_energy_conserving = use_energy_conserving, &
                  bcnd           = self%dom%bcnd, &
                  wall_cell      = self%dom%wall_cell, &
                  xmax           = self%dom%xmax, &
                  ymax           = self%dom%ymax, &
                  zmax           = self%dom%zmax, &
                  flag_pbc       = int(self%dom%flag_pbc, int32), &
                  flag_nmn       = int(self%dom%flag_nmn, int32), &
                  ptype          = ptype, &
                  tag_neg        = tag_neg_local, &
                  flag_die       = int(self%dom%flag_die, int32), &
                  dtype          = self%dom%dtype, &
                  qmacro         = qmacro, &
                  sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
                  sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
                  p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
                  P_loss_wall    = self%P_loss(1,ptype,iproc), &
                  Nm_species     = self%params%Nm(ptype), &
                  see            = see, &
                  part_electrons = self%part(1,iproc), &
                  iseed          = self%params%iseed(iproc), &
                  P_loss_see     = self%P_loss(4,1,iproc), &
                  flag_RFant     = self%cfg%flag_RFant, &
                  ixl_pow        = ixl_pow_rf, &
                  ixr_pow        = ixr_pow_rf, &
                  R_ahp          = self%cfg%R_ahp, &
                  gams           = self%gams, &
                  E0_RF          = self%cfg%E0_RF, &
                  omega_RF       = omega_RF_local, &
                  time           = time_rf, &
                  P_RF_local     = self%P_RF(ptype,iproc), &
                  flag_planar_ant = self%cfg%flag_planar_ant, &
                  x0             = x0_rf, &
                  mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
                  mom_loss_see   = self%mom_loss(:,4,1,iproc), &
                  mom_RF_local   = self%mom_RF(:,ptype,iproc), &
                  P_loss_coll    = self%P_loss(3,ptype,iproc), &
                  mom_loss_coll  = self%mom_loss(:,3,ptype,iproc) )
            end if
          end if
        else
          ! m_species<=0: no push (same guard as before this fusion),
          ! BC/SEE still applies - fall back to the standalone routine
          ! since move_and_bc_* always includes the push.
          call apply_particle_bc( &
              part           = self%part(ptype,iproc), &
              n              = int(self%dom%n, int32), &
              h              = self%dom%h, &
              bcnd           = self%dom%bcnd, &
              wall_cell      = self%dom%wall_cell, &
              xmax           = self%dom%xmax, &
              ymax           = self%dom%ymax, &
              zmax           = self%dom%zmax, &
              flag_pbc       = int(self%dom%flag_pbc, int32), &
              flag_nmn       = int(self%dom%flag_nmn, int32), &
              ptype          = ptype, &
              tag_neg        = tag_neg_local, &
              flag_die       = int(self%dom%flag_die, int32), &
              dtype          = self%dom%dtype, &
              qmacro         = qmacro, &
              sum_q_xz_local = self%sum_q_xz(:,:,ptype,iproc), &
              sum_q_yz_local = self%sum_q_yz(:,:,:,ptype,iproc), &
              p_mac_boundary = self%p_mac(ptype,:,:,iproc), &
              mass_species   = m_species, &
              P_loss_wall    = self%P_loss(1,ptype,iproc), &
              Nm_species     = self%params%Nm(ptype), &
              charge_species = q_species, &
              see            = see, &
              part_electrons = self%part(1,iproc), &
              iseed          = self%params%iseed(iproc), &
              P_loss_see     = self%P_loss(4,1,iproc), &
              mom_loss_wall  = self%mom_loss(:,1,ptype,iproc), &
              mom_loss_see   = self%mom_loss(:,4,1,iproc), &
              P_loss_coll    = self%P_loss(3,ptype,iproc), &
              mom_loss_coll  = self%mom_loss(:,3,ptype,iproc) )
        end if
        tw1 = omp_get_wtime()
        t_move_thread(iproc) = t_move_thread(iproc) + (tw1 - tw0)
      end do

      tw0 = omp_get_wtime()
      if (do_heat_tally) heat_loc = heat_tpl
      do ptype = 1, self%ntype
        ! Clear only the planes the previous deposit (+ periodic BC) could
        ! have made nonzero, then mark empty until this deposit reports.
        call clear_np_planes(self%np_thread(:,:,:,ptype,iproc), int(self%dom%n, int32), &
                             self%np_zlo(ptype,iproc), self%np_zhi(ptype,iproc), &
                             self%np_zend(ptype,iproc))
        self%np_zlo(ptype,iproc)  = huge(1_int32)
        self%np_zhi(ptype,iproc)  = -huge(1_int32)
        self%np_zend(ptype,iproc) = 0_int32

        if (.not. allocated(self%part(ptype,iproc)%pv)) cycle
        if (self%part(ptype,iproc)%n <= 0_int32) cycle

        if (ptype == 1_int32 .and. do_heat_tally) then
          call deposit_particle_set_to_np_thread( &
            part       = self%part(ptype,iproc), &
            n          = int(self%dom%n, int32), &
            h          = self%dom%h, &
            kq         = self%fld%kq, &
            Nm_species = self%params%Nm(ptype), &
            np_local   = self%np_thread(:,:,:,ptype,iproc), &
            heat       = heat_loc, &
            iz_lo      = self%np_zlo(ptype,iproc), &
            iz_hi      = self%np_zhi(ptype,iproc) )
        else
          call deposit_particle_set_to_np_thread( &
            part       = self%part(ptype,iproc), &
            n          = int(self%dom%n, int32), &
            h          = self%dom%h, &
            kq         = self%fld%kq, &
            Nm_species = self%params%Nm(ptype), &
            np_local   = self%np_thread(:,:,:,ptype,iproc), &
            iz_lo      = self%np_zlo(ptype,iproc), &
            iz_hi      = self%np_zhi(ptype,iproc) )
        end if

        ! apply_periodic_density_bc's z-stitching copies between planes
        ! {1,2,n3,n3+1} and {0,1,n3+1,n3+2}: if this deposit reached
        ! within two planes of either end, both end blocks may turn nonzero.
        if (self%np_zlo(ptype,iproc) <= 2_int32 .or. &
            self%np_zhi(ptype,iproc) >= int(self%dom%n(3), int32)) then
          self%np_zend(ptype,iproc) = 1_int32
        end if
      end do
      ! One write per thread (no per-particle updates of shared arrays).
      if (do_heat_tally) then
        self%Nh(iproc)      = heat_loc%Nh
        self%sum_dEk(iproc) = heat_loc%sum_dEk
        self%sum_wh(iproc)  = heat_loc%sum_w
      end if
      tw1 = omp_get_wtime()
      t_deposit_thread(iproc) = tw1 - tw0
    end do
    !$omp end parallel do

    t_move_out    = maxval(t_move_thread)
    t_bc_out      = maxval(t_bc_thread)
    t_deposit_out = maxval(t_deposit_thread)

    t_move_min_out    = minval(t_move_thread)
    t_move_avg_out    = sum(t_move_thread) / real(self%nproc, real64)
    t_deposit_min_out = minval(t_deposit_thread)
    t_deposit_avg_out = sum(t_deposit_thread) / real(self%nproc, real64)

    n_boris_used_out  = sum(n_boris_used_thread)
    n_boris_total_out = sum(n_boris_total_thread)

    deallocate(t_move_thread, t_bc_thread, t_deposit_thread)
    deallocate(n_boris_used_thread, n_boris_total_thread)
  end subroutine advance_particles_local


  subroutine compute_heating_region_moments(self)
    class(State), intent(inout) :: self
    integer(int32) :: iproc, i, ix
    real(real64)   :: x, y, z, v2
    real(real64)   :: ymax_half, zmax_half

    if (.not. allocated(self%part)) return
    if (.not. allocated(self%Nh)) return
    if (.not. allocated(self%sum_dEk)) return

    self%Nh      = 0_int32
    self%sum_dEk = 0.0_real64

    ymax_half = self%dom%ymax / 2.0_real64
    zmax_half = self%dom%zmax / 2.0_real64

    !$omp parallel do private(iproc,i,ix,x,y,z,v2) schedule(static)
    do iproc = 1, self%nproc
      if (.not. allocated(self%part(1,iproc)%pv)) cycle
      if (self%part(1,iproc)%n <= 0_int32) cycle

      do i = 1, self%part(1,iproc)%n
        x  = self%part(1,iproc)%pv(1,i)
        ix = int(x / self%dom%h(1), int32) + 1_int32

        if (ix < self%cfg%xl_pow/self%dom%h(1) .or. ix > self%cfg%xr_pow/self%dom%h(1)) cycle

        if (self%cfg%flag_circxh == 1) then
          y = self%part(1,iproc)%pv(2,i)
          z = self%part(1,iproc)%pv(3,i)

          if (self%cfg%flag_ahp == 0) then
            if (((y-ymax_half)**2 + (z-zmax_half)**2) > self%cfg%R_ahp**2) cycle
          else
            if (((y-ymax_half)**2 + (z-zmax_half)**2) < self%cfg%R_ahp**2) cycle
          end if
        end if

        v2 = self%part(1,iproc)%pv(4,i)**2 + &
             self%part(1,iproc)%pv(5,i)**2 + &
             self%part(1,iproc)%pv(6,i)**2

        self%sum_dEk(iproc) = self%sum_dEk(iproc) + &
            0.5_real64 * self%params%Nm(1) * self%chem%mass(1) * v2

        self%Nh(iproc) = self%Nh(iproc) + 1_int32
      end do
    end do
    !$omp end parallel do
  end subroutine compute_heating_region_moments


  subroutine update_heating_vt(self, vt_heat)
    use mpi
    class(State), intent(inout) :: self
    real(real64), intent(out)   :: vt_heat

    integer :: ierr
    integer(int32) :: sum_Nh, sum_Nh_global
    real(real64)   :: sum_dEk_tot, sum_dEk_global
    real(real64)   :: sum_w, sum_w_global
    real(real64)   :: Te

    vt_heat = self%params%vt0(1)

    if (.not. allocated(self%Nh)) return
    if (.not. allocated(self%sum_dEk)) return

    sum_Nh       = sum(self%Nh)
    sum_dEk_tot  = sum(self%sum_dEk)
    sum_w        = sum(self%sum_wh)

    sum_Nh_global  = sum_Nh
    sum_dEk_global = sum_dEk_tot
    sum_w_global   = sum_w

    if (self%mpi_size > 1) then
      call MPI_Allreduce(sum_Nh, sum_Nh_global, 1, MPI_INTEGER, MPI_SUM, self%comm, ierr)
      call MPI_Allreduce(sum_dEk_tot, sum_dEk_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      call MPI_Allreduce(sum_w, sum_w_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
    end if

    if (sum_Nh_global <= 0 .or. sum_w_global <= 0.0_real64) return

    ! Electron i is re-sampled with probability nudt*w_i, so the expected
    ! energy gain is nudt*sum_i w_i*(1.5*Nm*qe*Te - Ek_i); setting it to
    ! Pabs*dt_heat gives Te below (w = 1: sum_w = Nh, the legacy formula).
    Te = (2.0_real64/3.0_real64) * &
         (sum_dEk_global + self%cfg%Pabs / self%cfg%nu_h) / &
         (self%params%Nm(1) * qe * sum_w_global)

    if (Te > 0.0_real64) then
      vt_heat = sqrt(2.0_real64 * qe * Te / abs(self%chem%mass(1)))
    end if

    !if (self%mpi_rank == 0) &
     ! write(*,'(a,es12.4,a,i10,a,es12.4)') &
      !' HEAT_MOD vt=', vt_heat, ' Nh=', sum_Nh_global, ' sum_dEk=', sum_dEk_global

  end subroutine update_heating_vt


  subroutine update_heating_skin_depth(self, istep)
    ! Self-consistent skin depth of the Maxwellian heating profile
    ! w(r) = exp((r-R_ahp)/delta_heat): every RF period, delta_heat relaxes
    ! toward c/wp of the electron density averaged over the previous skin
    ! layer R_ahp-delta_heat <= r <= R_ahp of the heating region. Same
    ! grid loop as update_rf_convergence below, but over the layer the
    ! field reaches rather than the whole cross-section.
    class(State),   intent(inout) :: self
    integer(int32), intent(in)    :: istep

    ! Under-relaxation: delta can feed back on the layer it averages over
    real(real64), parameter :: alpha = 0.3_real64
    integer(int32) :: ix, iy, iz, cnt
    integer(int32) :: ixl, ixr, iyl, iyr, izl, izr
    real(real64)   :: yp, zp, r2, rin2, ymax_half, zmax_half
    real(real64)   :: np_layer, wp_layer, delta_calc, delta_min

    call self%sync_np()

    ymax_half = self%dom%ymax / 2.0_real64
    zmax_half = self%dom%zmax / 2.0_real64
    rin2 = max(0.0_real64, self%cfg%R_ahp - self%cfg%delta_heat)**2

    ixl = int(max(0.0_real64, self%cfg%xl_pow)/self%dom%h(1), int32) + 1_int32
    ixr = int(min(self%dom%xmax, self%cfg%xr_pow)/self%dom%h(1), int32) + 1_int32
    iyl = int(max(0.0_real64, ymax_half - self%cfg%R_ahp)/self%dom%h(2), int32) + 1_int32
    iyr = int(min(self%dom%ymax, ymax_half + self%cfg%R_ahp)/self%dom%h(2), int32) + 1_int32
    izl = int(max(0.0_real64, zmax_half - self%cfg%R_ahp)/self%dom%h(3), int32) + 1_int32
    izr = int(min(self%dom%zmax, zmax_half + self%cfg%R_ahp)/self%dom%h(3), int32) + 1_int32

    np_layer = 0.0_real64
    cnt      = 0_int32
    !$omp parallel do collapse(2) private(iy,iz,ix,yp,zp,r2) &
    !$omp&  reduction(+:np_layer,cnt) schedule(static)
    do iz = izl, izr
      do iy = iyl, iyr
        yp = real(iy-1, real64)*self%dom%h(2)
        zp = real(iz-1, real64)*self%dom%h(3)
        r2 = (yp-ymax_half)**2 + (zp-zmax_half)**2
        if (r2 >= rin2 .and. r2 <= self%cfg%R_ahp**2) then
          do ix = ixl, ixr
            np_layer = np_layer + self%fld%np(ix,iy,iz,1)
            cnt      = cnt + 1_int32
          end do
        end if
      end do
    end do
    !$omp end parallel do

    self%nstep_RF = 0_int32
    if (cnt <= 0_int32) return
    np_layer = np_layer / real(cnt, real64)
    if (np_layer <= 0.0_real64) return

    ! Physical eps0 (undo k_eps0), as in update_rf_convergence
    wp_layer   = sqrt( np_layer * self%chem%charge(1)**2 / &
                       ((eps0/self%cfg%k_eps0) * self%chem%mass(1)) )
    delta_calc = c / wp_layer

    ! Bounded between ~2 cells and the heating radius (beyond R the
    ! profile is essentially uniform anyway)
    delta_min = 2.0_real64 * max(self%dom%h(2), self%dom%h(3))
    self%cfg%delta_heat = (1.0_real64 - alpha)*self%cfg%delta_heat + alpha*delta_calc
    self%cfg%delta_heat = min(max(self%cfg%delta_heat, delta_min), self%cfg%R_ahp)

    if (self%mpi_rank == 0) then
      write(*,'(1x,a,i0,a,es10.2,a,f7.3,a,f7.3)') 'it=', istep, &
        ', <ne>_layer(/m^3)=', np_layer, ', c/wp(cm)=', delta_calc*1.0e2_real64, &
        ', heating skin depth delta(cm)=', self%cfg%delta_heat*1.0e2_real64
    end if
    call self%write_skin_depth_record(istep, np_layer, delta_calc, self%cfg%delta_heat, &
                                      0.0_real64, 0.0_real64)
  end subroutine update_heating_skin_depth


  subroutine update_rf_convergence(self, istep)
    ! RF antenna heating: periodic (every ns_RF steps) convergence step -
    ! recompute the skin depth from the actual local electron density, and
    ! rescale E0_RF toward the target absorbed power cfg%Pabs. Ported from
    ! legacy Src/main.f90:887-933. Follows the MPI-reduction idiom of
    ! update_heating_vt above and the OMP-reduction-over-grid idiom of
    ! compute_heating_region_moments above.
    use mpi
    class(State),   intent(inout) :: self
    integer(int32), intent(in)    :: istep

    integer :: ierr
    integer(int32) :: ix, iy, iz, cnt_RFcels
    integer(int32) :: ixl, ixr, iyl, iyr, izl, izr
    real(real64)   :: xl_c, xr_c, yl_c, yr_c, zl_c, zr_c
    real(real64)   :: yp, zp, ymax_half, zmax_half
    real(real64)   :: np_avgRF, wp_rf, P_RFtot, P_RFtot_local

    ! Clamp the heating-region bounds to the domain box before turning them
    ! into grid indices - cfg%xl_pow etc. are user input and, by convention,
    ! often set larger than the box to mean "the whole domain" (elsewhere
    ! they're only ever used as comparison bounds, where an oversized value
    ! is harmless). Used directly as array indices below, an unclamped
    ! oversized bound reads self%fld%np far out of bounds. Mirrors legacy
    ! Src/main.f90:680-692.
    call self%sync_np()

    xl_c = max(0.0_real64, self%cfg%xl_pow)
    xr_c = min(self%dom%xmax, self%cfg%xr_pow)
    yl_c = max(0.0_real64, self%cfg%yl_pow)
    yr_c = min(self%dom%ymax, self%cfg%yr_pow)
    zl_c = max(0.0_real64, self%cfg%zl_pow)
    zr_c = min(self%dom%zmax, self%cfg%zr_pow)
    ! Planar antenna: average over the window disk (r <= R_ahp about the
    ! y-z centre), whose radius need not match the yl_pow..zr_pow box.
    if (self%cfg%flag_planar_ant == 1_int32) then
      yl_c = max(0.0_real64,     self%dom%ymax/2.0_real64 - self%cfg%R_ahp)
      yr_c = min(self%dom%ymax, self%dom%ymax/2.0_real64 + self%cfg%R_ahp)
      zl_c = max(0.0_real64,     self%dom%zmax/2.0_real64 - self%cfg%R_ahp)
      zr_c = min(self%dom%zmax, self%dom%zmax/2.0_real64 + self%cfg%R_ahp)
    end if

    ixl = int(xl_c/self%dom%h(1), int32) + 1_int32
    ixr = int(xr_c/self%dom%h(1), int32) + 1_int32
    iyl = int(yl_c/self%dom%h(2), int32) + 1_int32
    iyr = int(yr_c/self%dom%h(2), int32) + 1_int32
    izl = int(zl_c/self%dom%h(3), int32) + 1_int32
    izr = int(zr_c/self%dom%h(3), int32) + 1_int32

    ymax_half = self%dom%ymax / 2.0_real64
    zmax_half = self%dom%zmax / 2.0_real64

    ! <np> over the heating volume, restricted to the circular
    ! cross-section of radius R_ahp - legacy Src/main.f90:888-909. Every
    ! rank holds the full self%fld%np grid (already MPI-Allreduced by
    ! reduce_species_density) and the same cfg-derived loop bounds, so no
    ! MPI reduction is needed here - matches legacy, which also doesn't
    ! reduce cnt_RFcels.
    np_avgRF   = 0.0_real64
    cnt_RFcels = 0_int32
    !$omp parallel do collapse(2) private(iy,iz,ix,yp,zp) &
    !$omp&  reduction(+:np_avgRF,cnt_RFcels) schedule(static)
    do iz = izl, izr
      do iy = iyl, iyr
        yp = real(iy-1, real64)*self%dom%h(2)
        zp = real(iz-1, real64)*self%dom%h(3)
        if (((yp-ymax_half)**2 + (zp-zmax_half)**2) <= self%cfg%R_ahp**2) then
          do ix = ixl, ixr
            np_avgRF   = np_avgRF + self%fld%np(ix,iy,iz,1)
            cnt_RFcels = cnt_RFcels + 1_int32
          end do
        end if
      end do
    end do
    !$omp end parallel do

    if (cnt_RFcels <= 0_int32) return
    np_avgRF = np_avgRF / real(cnt_RFcels, real64)

    ! wp/gams from the true physical eps0 (undo the k_eps0 rescale applied
    ! at init(), line ~135) - matches legacy Src/main.f90:910-912 exactly.
    wp_rf = sqrt( np_avgRF * self%chem%charge(1)**2 / &
                  ((eps0/self%cfg%k_eps0) * self%chem%mass(1)) )
    self%gams = c / wp_rf

    ! Reduce this period's absorbed power (accumulated by the mover into
    ! self%P_RF every step since the last call) across MPI ranks.
    P_RFtot_local = sum(self%P_RF(1:self%ntype,1:self%nproc)) / &
                    real(max(1_int32, self%nstep_RF), real64) / self%params%dt
    P_RFtot = P_RFtot_local
    if (self%mpi_size > 1) then
      call MPI_Allreduce(P_RFtot_local, P_RFtot, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, self%comm, ierr)
    end if

    ! Fold into the same power-loss counter print_diagnostics already
    ! reports as "Pabs (W)" - reuses P_loss(2,:,:), matching legacy's own
    ! reuse (Src/main.f90:923).
    self%P_loss(2,:,:)   = self%P_loss(2,:,:)   + self%P_RF(:,:)
    self%mom_loss(:,2,:,:) = self%mom_loss(:,2,:,:) + self%mom_RF(:,:,:)

    ! Iterative E0_RF toward the target cfg%Pabs.
    if (P_RFtot > 0.0_real64) then
      self%cfg%E0_RF = self%cfg%E0_RF * sqrt(self%cfg%Pabs / P_RFtot)
    end if
    if (self%mpi_size > 1) then
      call MPI_Bcast(self%cfg%E0_RF, 1, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    end if

    if (self%mpi_rank == 0) then
      write(*,'(1x,a,i0,a,es10.2,a,es10.2,a,es10.2,a,f6.2)') &
        'it=', istep, ', E0(V/m)=', self%cfg%E0_RF, &
        ', <P_RF>(W)=', P_RFtot, ', <np>(/m^3)=', np_avgRF, &
        ', skin depth gam(cm)=', self%gams*1.0e2_real64
    end if
    call self%write_skin_depth_record(istep, np_avgRF, self%gams, self%gams, &
                                      self%cfg%E0_RF, P_RFtot)

    self%P_RF     = 0.0_real64
    self%mom_RF   = 0.0_real64
    self%nstep_RF = 0_int32
  end subroutine update_rf_convergence


  subroutine write_skin_depth_record(self, istep, ne, delta_calc, delta, E0, P)
    ! Appends one line to Output/skin_depth.dat (rank 0) at every skin-depth
    ! update, i.e. every ns_RF steps: the value then holds until the next
    ! line. Written by update_rf_convergence (RF antenna: delta = gams,
    ! E0_RF and absorbed power filled) and update_heating_skin_depth
    ! (self-consistent Maxwellian heating: E0 = P = 0). Cleared at the
    ! start of a fresh run by Simulation%run, appended to on restart.
    class(State),   intent(in) :: self
    integer(int32), intent(in) :: istep
    real(real64),   intent(in) :: ne, delta_calc, delta, E0, P

    integer :: u
    logical :: file_exists

    if (self%mpi_rank /= 0) return
    inquire(file='./Output/skin_depth.dat', exist=file_exists)
    open(newunit=u, file='./Output/skin_depth.dat', status='unknown', &
         position='append', action='write')
    if (.not. file_exists) then
      write(u,'(a)') '# Time(s), it, <ne>(m^-3) over the skin layer (heating) or heated ' // &
          'volume (antenna), c/wp(cm), delta used(cm), E0_RF(V/m), P_RF(W)'
    end if
    write(u,'(1x,es16.8,1x,i12,5(1x,es16.8))') self%params%dt*real(istep, real64), istep, &
        ne, delta_calc*1.0e2_real64, delta*1.0e2_real64, E0, P
    close(u)
  end subroutine write_skin_depth_record


  subroutine compute_plane_moments_local(self)
    class(State), intent(inout) :: self
    integer(int32) :: ptype

    if (.not. allocated(self%part)) return

    do ptype = 1, self%ntype
      call compute_particle_plane_moments_species( &
          part         = self%part(ptype,:), &
          nproc        = self%nproc, &
          n            = int(self%dom%n, int32), &
          h            = self%dom%h, &
          ptype        = ptype, &
          mass_species = self%chem%mass(ptype), &
          ix_plane     = self%params%ix_plot_plane , &
          iy_plane     = int(self%dom%n(2)/2 +1, int32), &
          iz_plane     = self%params%iz_plot_plane, &
          data_xy      = self%data_pavg_xy(:,:,:,ptype), &
          data_xz      = self%data_pavg_xz(:,:,:,ptype), &
          data_yz      = self%data_pavg_yz(:,:,:,ptype), &
          params       = self%params)
    end do
  end subroutine compute_plane_moments_local

  subroutine sync_np(self)
    ! Collective: must be reached by every rank (np_synced changes in
    ! lockstep on all ranks, so the Allreduce below is always matched).
    class(State), intent(inout) :: self

    if (self%np_synced) return
    call sync_species_density(int(self%dom%n, int32), int(self%ntype), self%comm, self%fld%np)
    self%np_synced = .true.
  end subroutine sync_np


  subroutine accumulate_2d_averages(self)
    class(State), intent(inout) :: self

    integer(int32) :: ptype
    integer(int32) :: ix_plane, iy_density, iy_phi, iz_plane
    integer(int32) :: ny
    integer(int32) :: ix, iy, iz, ispec

    call self%sync_np()

    ny = self%dom%n(2)

    ix_plane   = self%params%ix_plot_plane
    iy_density = int(ny/2 + 1, int32)
    iy_phi     = int(ny/2,     int32)
    iz_plane   = self%params%iz_plot_plane

    do ptype = 1, self%ntype
      self%np_avg_xy(:,:,ptype) = self%np_avg_xy(:,:,ptype) + &
          self%fld%np(:,:,iz_plane,ptype)

      self%np_avg_xz(:,:,ptype) = self%np_avg_xz(:,:,ptype) + &
          self%fld%np(:,iy_density,:,ptype)

      self%np_avg_yz(:,:,ptype) = self%np_avg_yz(:,:,ptype) + &
          self%fld%np(ix_plane,:,:,ptype)
    end do

    self%phi_avg_xy(:,:) = self%phi_avg_xy(:,:) + self%fld%phi(:,:,iz_plane)
    self%phi_avg_xz(:,:) = self%phi_avg_xz(:,:) + self%fld%phi(:,iy_phi,:)
    self%phi_avg_yz(:,:) = self%phi_avg_yz(:,:) + self%fld%phi(ix_plane,:,:)

    ! Legacy main.f90 avg3D block: E, phi and the species-2 density. fld%np
    ! is already rank-summed by the sync_np() above (legacy's dens_red).
    if (allocated(self%avg3D)) then
      ispec = min(2_int32, self%ntype)
      !$omp parallel do private(ix,iy,iz) schedule(static)
      do iz = 1, self%dom%n(3)+1
        do iy = 1, self%dom%n(2)+1
          do ix = 1, self%dom%n(1)+1
            self%avg3D(1:3,ix,iy,iz) = self%avg3D(1:3,ix,iy,iz) + self%fld%E(1:3,ix,iy,iz)
            self%avg3D(4,ix,iy,iz)   = self%avg3D(4,ix,iy,iz)   + self%fld%phi(ix,iy,iz)
            self%avg3D(5,ix,iy,iz)   = self%avg3D(5,ix,iy,iz)   + self%fld%np(ix,iy,iz,ispec)
          end do
        end do
      end do
      !$omp end parallel do
    end if

    self%cnt_avg = self%cnt_avg + 1_int32

  end subroutine accumulate_2d_averages


  subroutine finalize_particles_only(self)
    class(State), intent(inout) :: self
    integer(int32) :: ptype, iproc

    if (.not. allocated(self%part)) return

    do ptype = 1, size(self%part,1)
      do iproc = 1, size(self%part,2)
        call self%part(ptype,iproc)%destroy()
      end do
    end do

    deallocate(self%part)
  end subroutine finalize_particles_only


  subroutine finalize(self)
    class(State), intent(inout) :: self

    call self%pdec%destroy()
    call self%fld%destroy()
    call self%rxn%destroy()
    call self%magField%destroy()
    call self%finalize_particles_only()

    if (allocated(self%np_thread))    deallocate(self%np_thread)
    if (allocated(self%np_zlo))       deallocate(self%np_zlo, self%np_zhi, self%np_zend)
    if (allocated(self%sum_q_xz))     deallocate(self%sum_q_xz)
    if (allocated(self%sum_q_yz))     deallocate(self%sum_q_yz)
    if (allocated(self%P_loss))       deallocate(self%P_loss)
    if (allocated(self%mom_loss))     deallocate(self%mom_loss)
    if (allocated(self%P_RF))         deallocate(self%P_RF)
    if (allocated(self%mom_RF))       deallocate(self%mom_RF)
    if (allocated(self%p_mac))        deallocate(self%p_mac)
    if (allocated(self%sour_xy))      deallocate(self%sour_xy)
    if (allocated(self%sour_xz))      deallocate(self%sour_xz)
    if (allocated(self%sour_yz))      deallocate(self%sour_yz)
    if (allocated(self%sour_fx_yz))   deallocate(self%sour_fx_yz)
    if (allocated(self%N_inj))        deallocate(self%N_inj)
    if (allocated(self%N_flx))        deallocate(self%N_flx)
    if (allocated(self%Nh))           deallocate(self%Nh)
    if (allocated(self%sum_dEk))      deallocate(self%sum_dEk)
    if (allocated(self%np_avg_xy))    deallocate(self%np_avg_xy)
    if (allocated(self%np_avg_xz))    deallocate(self%np_avg_xz)
    if (allocated(self%np_avg_yz))    deallocate(self%np_avg_yz)
    if (allocated(self%phi_avg_xy))   deallocate(self%phi_avg_xy)
    if (allocated(self%phi_avg_xz))   deallocate(self%phi_avg_xz)
    if (allocated(self%phi_avg_yz))   deallocate(self%phi_avg_yz)
    if (allocated(self%avg3D))        deallocate(self%avg3D)
    if (allocated(self%data_pavg_xy)) deallocate(self%data_pavg_xy)
    if (allocated(self%data_pavg_xz)) deallocate(self%data_pavg_xz)
    if (allocated(self%data_pavg_yz)) deallocate(self%data_pavg_yz)
    if (allocated(self%sour_avg_xy))  deallocate(self%sour_avg_xy)
    if (allocated(self%sour_avg_xz))  deallocate(self%sour_avg_xz)
    if (allocated(self%sour_avg_yz))  deallocate(self%sour_avg_yz)
    if (allocated(self%sink_avg_xy))  deallocate(self%sink_avg_xy)
    if (allocated(self%sink_avg_xz))  deallocate(self%sink_avg_xz)
    if (allocated(self%sink_avg_yz))  deallocate(self%sink_avg_yz)

    if (allocated(self%params%iseed)) deallocate(self%params%iseed)
  end subroutine finalize

end module mod_state
