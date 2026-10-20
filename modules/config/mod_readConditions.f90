module mod_readConditions
  use iso_fortran_env, only: real64
  use mod_utils
  use mod_config, only: Config, nm_rg, npart, NBMAX
  implicit none
  private
  public :: read_input

  interface read_input
    module procedure read_input_cfg
    module procedure read_input_legacy
  end interface read_input
 contains
  ! --------------------------------------------------------------------
  ! Legacy kernel: EXACT signature, EXACT logic.
  ! Only edits: remove include/COMMON, rely on mod_part_info globals.
  ! --------------------------------------------------------------------
    subroutine read_input_legacy(n,tmax,xl_pow,xr_pow,yl_pow,yr_pow,zl_pow,zr_pow,&
        nsav,eps,omega,kt,rname,ngrid,ng,xl_rg,xr_rg,yl_rg,yr_rg,zl_rg,zr_rg,&
        I_inj,Ca,mpi_rank,flag_avg3D,np_dup,n_B,phi0_RF,f0_RF,phi1_RF,f1_RF)

        implicit none
        integer,      intent(out) :: n(3), ngrid, ng, nsav, flag_avg3D, n_B(3)
        integer,      intent(in)  :: mpi_rank
        real(real64), intent(out) :: tmax, eps, omega, kt
        real(real64), intent(out) :: xl_pow, xr_pow, yl_pow, yr_pow, zl_pow, zr_pow
        real(real64), intent(out) :: xl_rg(nm_rg), xr_rg(nm_rg), yl_rg, yr_rg, zl_rg, zr_rg
        real(real64), intent(out) :: I_inj, Ca, np_dup
        character(len=*), intent(out) :: rname
        real(real64), intent(out) :: phi0_RF, f0_RF, phi1_RF, f1_RF

        type(Config) :: cfg

        call read_input_cfg(cfg, mpi_rank)

        ! copy to legacy outputs
        n = cfg%n
        tmax = cfg%tmax
        xl_pow = cfg%xl_pow; xr_pow = cfg%xr_pow
        yl_pow = cfg%yl_pow; yr_pow = cfg%yr_pow
        zl_pow = cfg%zl_pow; zr_pow = cfg%zr_pow


        nsav = cfg%nsav
        eps = cfg%eps
        omega = cfg%omega
        kt = cfg%kt
        rname = cfg%rname
        ngrid = cfg%ngrid
        ng = cfg%ng
        xl_rg = cfg%xl_rg
        xr_rg = cfg%xr_rg
        yl_rg = cfg%yl_rg; yr_rg = cfg%yr_rg
        zl_rg = cfg%zl_rg; zr_rg = cfg%zr_rg
        I_inj = cfg%I_inj
        Ca = cfg%Ca
        flag_avg3D = cfg%flag_avg3D
        np_dup = cfg%np_dup
        n_B = cfg%n_B
        phi0_RF = cfg%phi0_RF; f0_RF = cfg%f0_RF
        phi1_RF = cfg%phi1_RF; f1_RF = cfg%f1_RF
    end subroutine read_input_legacy



  ! --------------------------------------------------------------------
  ! OOP entry point:
  ! call read_input(cfg, mpi_rank)
  !
  ! For now: call legacy kernel, then copy the results + globals into cfg.
  ! This is the safest "first step" (no physics changes).
  ! --------------------------------------------------------------------
    subroutine read_input_cfg(cfg, mpi_rank)
        implicit none
        type(Config), intent(inout) :: cfg
        integer,      intent(in)    :: mpi_rank

        integer :: iB, flag_read, i_rg
        character(len=3) :: end_file
        character(len=256) :: line1
        character(len=1) :: ans
        character(len=256) :: line_rf, line_heat
        character(len=16)  :: rf_mode
        integer :: ios

        ! locals for signature outputs (still stored in cfg)
        integer :: flag_avg3D

        ! Init things that must reset each read
        cfg%xl_rg = 0.0_real64
        cfg%xr_rg = 0.0_real64
        cfg%flag_grd = 0

        open(10,file='./input_dir/conditions.inp')

        ! Read line 1 as a whole record first so a missing nsteps is an
        ! error instead of silently consuming a value from line 2.
        read(10,'(a)',end=999) line1
        read(line1,*,err=998,end=998) cfg%rname, cfg%nsteps
        if (cfg%nsteps <= 0) goto 998
        read(10,*,end=999) cfg%Ti(1)

        read(10,*,err=999) cfg%ng
        read(10,*,err=999) cfg%eps, cfg%omega
        read(10,*,err=999) cfg%n(1), cfg%n(2), cfg%n(3)

        read(10,*,err=999) cfg%xl_pow, cfg%xr_pow, cfg%yl_pow, cfg%yr_pow, cfg%zl_pow, cfg%zr_pow
        cfg%xl_pow = cfg%xl_pow*1.0e-2_real64
        cfg%xr_pow = cfg%xr_pow*1.0e-2_real64
        cfg%yl_pow = cfg%yl_pow*1.0e-2_real64
        cfg%yr_pow = cfg%yr_pow*1.0e-2_real64
        cfg%zl_pow = cfg%zl_pow*1.0e-2_real64
        cfg%zr_pow = cfg%zr_pow*1.0e-2_real64

        read(10,*,err=999) cfg%x_load
        cfg%x_load = cfg%x_load*1.0e-2_real64

        cfg%flag_ahp = 0
        if (cfg%yr_pow < 0.0_real64) then
            cfg%flag_ahp = 1
            cfg%yr_pow = abs(cfg%yr_pow)
        end if

        read(10,*,err=999) cfg%nB
        cfg%flag_B_pos = 0
        if (cfg%nB < 0) cfg%flag_B_pos = 1
        cfg%nB = abs(cfg%nB)

        read(10,*,err=999) cfg%n_B(1), cfg%n_B(2), cfg%n_B(3)

        cfg%flag_B = 0
        do iB=1,cfg%nB
            read(10,*,err=999) cfg%B_file(iB), cfg%B_name(iB), cfg%B_scale(iB), cfg%B0(iB), &
                            cfg%B_info(iB), cfg%dL(iB), cfg%x0(iB), cfg%y0(iB), cfg%z0(iB)
            cfg%B0(iB) = cfg%B0(iB)*1.0e-4_real64
            if (cfg%B_file(iB) /= 's' .and. cfg%B_file(iB) /= 'S') cfg%flag_B = 1
        end do

        if (cfg%n_B(1) <= 1 .or. cfg%n_B(2) <= 1 .or. cfg%n_B(3) <= 1 .or. cfg%flag_B == 0 .or. cfg%nB == 0) cfg%n_B = 1
        if (cfg%flag_B == 0) cfg%nB = 0

        cfg%dL = cfg%dL*1.0e-2_real64
        cfg%x0 = cfg%x0*1.0e-2_real64
        cfg%y0 = cfg%y0*1.0e-2_real64
        cfg%z0 = cfg%z0*1.0e-2_real64

        ! Pabs, nu_h [, delta_heat(cm)]: same one-record trick as the RF line
        ! below, so an absent optional value can't eat the next line.
        read(10,'(a)',err=999,end=999) line_heat
        read(line_heat,*,iostat=ios) cfg%Pabs, cfg%nu_h, cfg%delta_heat
        if (ios /= 0) then
            cfg%delta_heat = 0.0_real64
            read(line_heat,*,err=999) cfg%Pabs, cfg%nu_h
        end if
        cfg%flag_delta_sc = 0
        if (cfg%delta_heat < 0.0_real64) cfg%flag_delta_sc = 1
        cfg%delta_heat = abs(cfg%delta_heat)*1.0e-2_real64
        cfg%flag_heat = 1
        if (cfg%nu_h < 0.0_real64) cfg%flag_heat = 0

        read(10,*,err=999) cfg%I_inj, cfg%opt_inj
        cfg%flag_inj = 1
        if (cfg%I_inj < 0.0_real64) cfg%flag_inj = 0

        read(10,*,err=999) cfg%tmax
        read(10,*,err=999) cfg%kt
        read(10,*,err=999) cfg%n0
        read(10,*,err=999) cfg%ngas
        read(10,*,err=999) cfg%np_cell
        read(10,*,err=999) cfg%ngrid

        flag_avg3D = 1
        if (cfg%ngrid < 0) then
            flag_avg3D = 0
            cfg%ngrid  = abs(cfg%ngrid)
        end if
        cfg%flag_avg3D = flag_avg3D

        read(10,*,err=999) cfg%nsav
        read(10,*,err=999) cfg%nbak, cfg%cnt_plt, cfg%tseq_init, cfg%tseq_final
        cfg%tseq = abs(cfg%tseq_final - cfg%tseq_init)/real(cfg%cnt_plt, real64)

        read(10,*,err=999) cfg%k_eps0
        cfg%flag_convP = 0
        if (cfg%k_eps0 < 0.0_real64) then
            cfg%flag_convP = 1
            cfg%k_eps0 = abs(cfg%k_eps0)
        end if

        read(10,*,err=999) ans, cfg%ptype_pdf, cfg%yl_rg, cfg%yr_rg, cfg%zl_rg, cfg%zr_rg, cfg%n_rg, &
                            (cfg%xl_rg(i_rg), cfg%xr_rg(i_rg), i_rg=1,cfg%n_rg)

        cfg%zl_rg = cfg%zl_rg*1.0e-2_real64
        cfg%zr_rg = cfg%zr_rg*1.0e-2_real64
        cfg%yl_rg = cfg%yl_rg*1.0e-2_real64
        cfg%yr_rg = cfg%yr_rg*1.0e-2_real64
        cfg%xl_rg = cfg%xl_rg*1.0e-2_real64
        cfg%xr_rg = cfg%xr_rg*1.0e-2_real64

        if (mpi_rank == 0 .and. cfg%n_rg > nm_rg) then
            print*, 'Warning: nm_rg is too small, please correct ...'
            call stop_calculation
        end if

        cfg%flag_pdf = 0
        if (ans=='y' .or. ans=='Y') cfg%flag_pdf = 1

        read(10,*,err=999) ans, cfg%omp_rank_max, cfg%mpi_rank_max, cfg%np_dup
        cfg%flag_restart = 0
        if (ans=='y' .or. ans=='Y') cfg%flag_restart = 1
        if (ans=='e' .or. ans=='E') cfg%flag_restart = 2
        if (cfg%flag_restart == 0) cfg%np_dup = 1.0_real64

        read(10,*,err=999) ans, cfg%Lhy, cfg%Lhz, cfg%nhy, cfg%nhz, cfg%ind_g
        cfg%flag_grd = 0
        if (ans=='y' .or. ans=='Y') cfg%flag_grd = 1

        read(10,*,err=999) cfg%jne, cfg%THm, cfg%num_grd, cfg%Ca

        ! gam_sec, igrid_sec, RF mode, phi0_RF, f0_RF, phi1_RF, f1_RF [, E0_RF [, R_win]]
        ! Read as one record so the optional trailing values can't consume
        ! tokens from the next line: try 9 values, then 8, then 7 (a trailing
        ! '!' comment makes the longer reads fail and fall through).
        read(10,'(a)',err=999,end=999) line_rf
        read(line_rf,*,iostat=ios) cfg%gam_sec, cfg%igrid_sec, rf_mode, cfg%phi0_RF, &
                                   cfg%f0_RF, cfg%phi1_RF, cfg%f1_RF, cfg%E0_RF, cfg%R_win
        if (ios /= 0) then
            cfg%R_win = 0.0_real64
            read(line_rf,*,iostat=ios) cfg%gam_sec, cfg%igrid_sec, rf_mode, cfg%phi0_RF, &
                                       cfg%f0_RF, cfg%phi1_RF, cfg%f1_RF, cfg%E0_RF
        end if
        if (ios /= 0) then
            cfg%E0_RF = 0.0_real64
            read(line_rf,*,err=999) cfg%gam_sec, cfg%igrid_sec, rf_mode, cfg%phi0_RF, &
                                    cfg%f0_RF, cfg%phi1_RF, cfg%f1_RF
        end if
        cfg%R_win = cfg%R_win*1.0e-2_real64

        ! RF mode (case-insensitive; legacy single letters still accepted):
        !   'No'        : no RF
        !   'Electrode' : RF potential on electrode igrid_sec (legacy 'D', modular 'y')
        !   'Wrapped'   : inductive coil wrapped around the driver wall (legacy 'A').
        !                 Field peaks at R_ahp and decays radially inward (skin depth).
        !   'Planar'    : inductive planar coil under a dielectric window at x=xl_pow
        !                 (legacy-modular 'P'). Field decays axially away from the
        !                 window; radial profile bounded at the window radius R_win.
        cfg%flag_RFpot      = 0
        cfg%flag_RFant      = 0
        cfg%flag_planar_ant = 0
        select case (to_lower(rf_mode))
        case ('no', 'n', 'none', 'off')
        case ('electrode', 'y', 'd')
            cfg%flag_RFpot = 1
        case ('wrapped', 'a')
            cfg%flag_RFant = 1
        case ('planar', 'p')
            cfg%flag_RFant      = 1
            cfg%flag_planar_ant = 1
        case default
            if (mpi_rank == 0) write(*,'(a)') 'Unknown RF mode "' // trim(rf_mode) // &
                '" in conditions.inp (expected No, Electrode, Wrapped or Planar)'
            call stop_calculation
        end select

        if (cfg%flag_RFpot == 1) then
            if (mpi_rank == 0) write(*,'(a)') 'RF mode "Electrode" is not implemented in ' // &
                'LePIC+ yet (legacy Src/main.f90:972-982 was not ported), please correct ...'
            call stop_calculation
        end if

        if (cfg%flag_RFant == 1) then
            ! RF antenna (inductive) heating - Src/read_input.f90:128-134
            cfg%flag_inj  = 0
            cfg%I_inj     = 0.0_real64
            cfg%flag_heat = 0
            ! E0_RF is rescaled multiplicatively toward Pabs and ns_RF = 1/(f0_RF*dt),
            ! so all three must be positive.
            if (cfg%f0_RF <= 0.0_real64 .or. cfg%E0_RF <= 0.0_real64 .or. cfg%Pabs <= 0.0_real64) then
                if (mpi_rank == 0) write(*,'(a)') 'RF antenna needs f0_RF > 0, E0_RF > 0 ' // &
                    '(initial guess) and Pabs > 0 (target power), please correct ...'
                call stop_calculation
            end if
        end if

        if (cfg%flag_delta_sc == 1 .and. cfg%flag_heat == 1 .and. cfg%f0_RF <= 0.0_real64) then
            if (mpi_rank == 0) write(*,'(a)') 'Self-consistent skin depth (delta < 0) needs ' // &
                'f0 > 0 on the RF line (update period 1/f0), please correct ...'
            call stop_calculation
        end if

        if (cfg%R_win > 0.0_real64 .and. cfg%flag_planar_ant == 0 .and. mpi_rank == 0) &
            write(*,'(a)') 'Warning: R_win is only used by the Planar RF antenna, ignored ...'

        if (cfg%I_inj > 0.0_real64 .and. cfg%gam_sec > 0.0_real64) then
            print*, 'Warning: gam_sec>0 is incompatible with I_inj>0, please correct ...'
            call stop_calculation
        end if

        read(10,*,err=999) ans
        cfg%flag_coulomb = 0
        if (ans=='y' .or. ans=='Y') cfg%flag_coulomb = 1

        ! Modular-only line - legacy's Src/read_input.f90 has no matching
        ! read, but tolerates it fine: its own END-search loop just skips
        ! any line that isn't END, same as it already does for the
        ! Coulomb-collisions line above. 'momentum' (default, unchanged
        ! behavior) or 'energy' (Powis & Kaganovich EC-PIC scheme) - see
        ! cfg%push_scheme's declaration comment (mod_config.f90) and the
        ! push_scheme guard in mod_simulation.f90's init().
        read(10,*,err=999) cfg%push_scheme

        flag_read = 0
        do while (flag_read == 0)
            read(10,*,err=999) end_file
            if (end_file=='END' .or. end_file=='end') flag_read = 1
        end do

        if (mpi_rank == 0) write(*,"(a)"), 'Input file was read correctly ...'
        close(10)
        return

        999 if (mpi_rank == 0) write(*,"(a)"), 'Input file was not read correctly!'
        call stop_calculation
        return

        998 if (mpi_rank == 0) write(*,"(a)") 'conditions.inp line 1 must be: ' // &
            '<reactions file> <number of time steps (> 0)>'
        call stop_calculation
    end subroutine read_input_cfg

    pure function to_lower(s) result(r)
        character(len=*), intent(in) :: s
        character(len=len_trim(s))   :: r
        integer :: i, ic
        r = s(1:len_trim(s))
        do i = 1, len(r)
            ic = iachar(r(i:i))
            if (ic >= iachar('A') .and. ic <= iachar('Z')) r(i:i) = achar(ic + 32)
        end do
    end function to_lower

end module mod_readConditions
