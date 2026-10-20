module mod_config
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: Config, NBMAX, nm_rg, npart

  integer, parameter :: NBMAX = 5
  integer, parameter :: nm_rg = 5
  integer, parameter :: npart = 32

  type :: Config
    integer           :: n(3)      = 0
    real(real64)      :: tmax      = 0.0_real64
    real(real64)      :: xl_pow    = 0.0_real64, xr_pow = 0.0_real64
    real(real64)      :: yl_pow    = 0.0_real64, yr_pow = 0.0_real64
    integer           :: flag_circxh = 0
    real(real64)      :: R_ahp    = 0.0_real64
    real(real64)      :: zl_pow    = 0.0_real64, zr_pow = 0.0_real64
    integer           :: nsav      = 0
    real(real64)      :: eps       = 0.0_real64
    real(real64)      :: omega     = 0.0_real64
    real(real64)      :: kt        = 0.0_real64
    character(len=20) :: rname     = ""
    ! Number of PIC time steps to run (sim%run). Read from the first line of
    ! conditions.inp, right after rname - REQUIRED (must be > 0). legacy's
    ! Src/read_input.f90 reads only rname from that line and ignores the rest.
    integer           :: nsteps    = 0
    integer           :: ngrid     = 0
    integer           :: ng        = 0

    integer           :: n_B(3)    = 1
    integer           :: flag_avg3D = 1
    real(real64)      :: np_dup    = 1.0_real64

    real(real64)      :: xl_rg(nm_rg) = 0.0_real64
    real(real64)      :: xr_rg(nm_rg) = 0.0_real64
    real(real64)      :: yl_rg     = 0.0_real64, yr_rg = 0.0_real64
    real(real64)      :: zl_rg     = 0.0_real64, zr_rg = 0.0_real64

    real(real64)      :: I_inj     = 0.0_real64
    real(real64)      :: Ca        = 0.0_real64

    real(real64)      :: phi0_RF   = 0.0_real64, f0_RF = 0.0_real64
    real(real64)      :: phi1_RF   = 0.0_real64, f1_RF = 0.0_real64
    real(real64)      :: E0_RF     = 0.0_real64
    ! Planar antenna only: radius (m) of the dielectric window/coil, centred
    ! on the y-z axis. <= 0 falls back to R_ahp = yr_pow - ymax/2.
    real(real64)      :: R_win     = 0.0_real64
    ! Maxwellian heating only: skin depth (m) of a wall-localized heating
    ! profile w(r) = exp((r-R_ahp)/delta_heat) inside the circular heating
    ! cross-section. <= 0 keeps the uniform heating (w = 1).
    real(real64)      :: delta_heat = 0.0_real64
    ! 1: delta_heat is recomputed every RF period (1/f0_RF) as c/wp of the
    ! electron density averaged over the previous skin layer (input < 0,
    ! |input| = initial guess). See State%update_heating_skin_depth.
    integer           :: flag_delta_sc = 0

    real(real64)      :: Ti(npart)     = 0.0_real64
    real(real64)      :: x_load    = 0.0_real64
    integer           :: flag_ahp  = 0

    integer           :: flag_B_pos = 0
    integer           :: flag_B     = 0
    integer           :: nB         = 0

    character(len=1)  :: B_file(NBMAX)  = " "
    character(len=20) :: B_name(NBMAX)  = ""
    character(len=2)  :: B_info(NBMAX)  = ""
    real(real64)      :: B_scale(NBMAX) = 0.0_real64
    real(real64)      :: B0(NBMAX)      = 0.0_real64
    real(real64)      :: dL(NBMAX)      = 0.0_real64
    real(real64)      :: x0(NBMAX)      = 0.0_real64
    real(real64)      :: y0(NBMAX)      = 0.0_real64
    real(real64)      :: z0(NBMAX)      = 0.0_real64

    real(real64)      :: Pabs      = 0.0_real64
    real(real64)      :: nu_h      = 0.0_real64
    integer           :: opt_inj   = 0
    integer           :: flag_heat = 0
    integer           :: flag_inj     = 0
    integer           :: flag_coulomb = 0  ! 1 = enable Coulomb collisions (Nanbu 2000)

    ! Selects the E-field calc + gather scheme used by advance_particles_local
    ! (mod_state.f90) and calc_Efield_* (mod_electricField.f90):
    !   'momentum' - node-centered E, same trilinear shape function for
    !                deposit and gather (default, matches legacy/existing
    !                behavior exactly).
    !   'energy'   - explicit energy-conserving scheme (Powis & Kaganovich,
    !                Phys. Plasmas 31, 023901 (2024)): face-centered E,
    !                nearest-grid-point gather normal to each component,
    !                linear gather tangentially. Supports plain wall,
    !                y/z-periodic (flag_pbc/flag_pbcz), and dielectric
    !                (flag_die) boundaries. Neumann (flag_nmn) is not yet
    !                supported - see mod_simulation.f90's push_scheme guard.
    ! Read from conditions.inp (mod_readConditions.f90), the last line
    ! before END - REQUIRED, not optional: a conditions.inp missing this
    ! line makes the read hit EOF, which is a hard failure (err=999 ->
    ! stop_calculation), not a fallback to the default below. Every
    ! conditions.inp (including archived per-case copies under input_dir/)
    ! must have this line. legacy's Src/read_input.f90 has no matching read
    ! for it, but tolerates its presence fine - its own END-search loop
    ! skips any non-END line, same as it already does for the
    ! Coulomb-collisions line just above this one in the file.
    character(len=20) :: push_scheme = 'momentum'

    real(real64)      :: n0        = 0.0_real64
    real(real64)      :: ngas      = 0.0_real64
    integer           :: np_cell   = 0

    integer           :: nbak      = 0
    integer           :: cnt_plt   = 0
    real(real64)      :: tseq_init = 0.0_real64
    real(real64)      :: tseq_final= 0.0_real64
    real(real64)      :: tseq      = 0.0_real64

    real(real64)      :: k_eps0    = 0.0_real64
    integer           :: flag_convP= 0

    integer           :: flag_pdf  = 0
    integer           :: ptype_pdf = 0
    integer           :: n_rg      = 0

    integer           :: omp_rank_max = 1
    integer           :: mpi_rank_max = 1
    integer           :: flag_restart = 0

    real(real64)      :: Lhy = 0.0_real64
    real(real64)      :: Lhz = 0.0_real64
    integer           :: nhy = 0
    integer           :: nhz = 0
    integer           :: ind_g = 0

    real(real64)      :: jne = 0.0_real64
    real(real64)      :: THm = 0.0_real64
    integer           :: num_grd = 0

    real(real64)      :: gam_sec  = 0.0_real64
    integer           :: igrid_sec = 0

    integer           :: flag_RFpot = 0
    integer           :: flag_RFant = 0
    integer           :: flag_planar_ant = 0
    integer           :: flag_grd   = 0
  end type Config

end module mod_config
