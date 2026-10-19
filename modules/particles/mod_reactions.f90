module mod_reactions
  use iso_fortran_env, only: real64
  use mpi
  use mod_utils,          only: stop_calculation, indexx
  use mod_chemistryState, only: ChemistryState, ind_nre, ind_nby, ind_Eth, ind_dE
  use mod_constants
  implicit none
  private
  public :: read_reactions

contains

  subroutine read_reactions(chem, sig, sig_Er, sig_list, sig_Eex, ncol_mx, sig_type, &
                            npt_mx, rname, col_info, scol_rank, scol_info, sigv_mx, &
                            ni0, ntype, n_neu, mpi_rank)

    type(ChemistryState), intent(inout) :: chem
    integer,              intent(in)    :: ncol_mx, npt_mx, mpi_rank
    character(len=*),     intent(in)    :: rname

    ! NOTE: keep arrays external for now (fastest path)
    real(real64), intent(inout) :: sig_Er(:)                 ! (npt_mx)
    real(real64), intent(inout) :: sig(:,:)                  ! (npt_mx, ncol_mx)
    integer,      intent(inout) :: sig_list(:,:)             ! (npart, ncol_mx)
    real(real64), intent(inout) :: sig_Eex(:,:)              ! (ncol_mx, 2)
    integer,      intent(inout) :: sig_type(:)               ! (ncol_mx)
    integer,      intent(inout) :: col_info(:,:)             ! (ncol_mx, 10)
    integer,      intent(inout) :: scol_rank(:,:)            ! (npart, ncol_mx)
    real(real64), intent(inout) :: scol_info(:,:)            ! (ncol_mx, 4)
    real(real64), intent(inout) :: sigv_mx(:,:)              ! (npart, ncol_mx)
    real(real64), intent(inout) :: ni0(:)                    ! (npart)

    integer,      intent(out)   :: ntype, n_neu

    ! locals
    integer :: ipt, i, j, icol, ptype, flag_pa, flag_re, i_re
    integer :: r_sy, l_sy, n_re, n_by, ind_label, iscol
    integer :: sig_npt(ncol_mx), indx(ncol_mx)
    real(real64) :: sort_arr(ncol_mx), sig_scale(2), Ps
    real(real64), allocatable :: sig_tmp(:,:,:)
    character(len=4)  :: name
    character(len=50) :: rlabel
    character(len=6)  :: tname(2)
    character(len=20) :: rtype
    integer, parameter :: jmax = 10000

    ! -------- sanity checks --------
    if (chem%npart <= 0) then
      if (mpi_rank == 0) print*, "ChemistryState not initialized: chem%npart <= 0"
      call stop_calculation
    end if

    if (size(sig_Er) < npt_mx) call stop_calculation
    if (size(sig,1) < npt_mx .or. size(sig,2) < ncol_mx) call stop_calculation
    if (size(sig_list,1) < chem%npart .or. size(sig_list,2) < ncol_mx) call stop_calculation
    if (size(sig_type) < ncol_mx) call stop_calculation
    if (size(sig_Eex,1) < ncol_mx .or. size(sig_Eex,2) < 2) call stop_calculation
    if (size(col_info,1) < ncol_mx .or. size(col_info,2) < 10) call stop_calculation
    if (size(scol_rank,1) < chem%npart .or. size(scol_rank,2) < ncol_mx) call stop_calculation
    if (size(scol_info,1) < ncol_mx .or. size(scol_info,2) < 4) call stop_calculation
    if (size(sigv_mx,1) < chem%npart .or. size(sigv_mx,2) < ncol_mx) call stop_calculation
    if (size(ni0) < chem%npart) call stop_calculation

    if (.not.allocated(chem%pname))  call stop_calculation
    if (.not.allocated(chem%mass))   call stop_calculation
    if (.not.allocated(chem%charge)) call stop_calculation
    if (.not.allocated(chem%Ti))     call stop_calculation
    if (.not.allocated(chem%ni0))    call stop_calculation
    if (.not.allocated(chem%p_ncol)) call stop_calculation
    if (.not.allocated(chem%p_nscol))call stop_calculation

    allocate(sig_tmp(npt_mx, ncol_mx, 2))

    ! Open files
    open(10, file='./input_dir/'//trim(rname))
    open(11, file='./Output/particles.out')

    ! Init arrays
    do icol=1,ncol_mx
      sig_npt(icol)  = 0
      sig_type(icol) = 0
      do ipt=1,npt_mx
        sig_Er(ipt)   = 0.0_real64
        sig(ipt,icol) = 0.0_real64
        sig_tmp(ipt,icol,1) = 0.0_real64
        sig_tmp(ipt,icol,2) = 0.0_real64
      end do
    end do

    icol     = 0
    flag_pa  = 0
    flag_re  = 0
    col_info = 0

    chem%p_ncol  = 0
    chem%p_nscol = 0
    iscol    = 0
    sort_arr = 0.0_real64
    indx     = 0

    chem%tag_neg  = 0
    chem%tag_neu  = 0
    chem%tag_beam = 0
    n_neu = 0

    ! Electron
    chem%pname(1)  = '[e]'
    chem%charge(1) = -qe
    chem%mass(1)   = 9.10938188e-31_real64
    ni0(1)         = 1.0_real64
    ntype          = 1

    ! ----------------
    ! Scan input file
    ! ----------------
    do j=1,jmax

19    read(10,*,END=101,ERR=100) name

      if (trim(name) == 'IONS' .or. trim(name) == 'Ions' .or. trim(name) == 'ions') then
        flag_pa = 1
        if (mpi_rank==0) then
          write(11,*) ' '
          write(11,*) 'Ions:'
        end if
        goto 20
      end if

      if (trim(name) == 'NEUT' .or. trim(name) == 'Neut' .or. trim(name) == 'neut') then
        flag_pa = 1
        if (mpi_rank==0) then
          write(11,*) ' '
          write(11,*) 'Neutrals:'
        end if
        goto 20
      end if

      if (trim(name) == 'REAC' .or. trim(name) == 'Reac' .or. trim(name) == 'reac') then
        if (flag_pa==0) goto 101
        icol = icol + 1
        if (mpi_rank==0 .and. icol>ncol_mx) then
          print*, 'icol > ncol_mx, please correct.'
          print*, 'Abort simulation ...'
          call stop_calculation
        end if
        sig_type(icol) = 1
        goto 21
      end if

      if (trim(name) == 'SURF' .or. trim(name) == 'Surf' .or. trim(name) == 'surf') then
        if (flag_pa==0) goto 101
        goto 22
      end if

      goto 19

      ! ---------------------------
      ! Ions / Neutrals definition
      ! ---------------------------
20    do while (trim(name) /= '----')
        read(10,*,END=101,ERR=100) name
      end do

200   read(10,fmt='(A6)',END=101,ERR=100,advance='no') tname(1)

      if (tname(1)(1:2) == '--') then
        if (mpi_rank==0 .and. ntype>chem%npart) then
          print*, 'Warning ntype>npart, please correct ...'
          call stop_calculation
        end if
        goto 99
      else
        ntype = ntype + 1
        if (ntype > chem%npart) then
          if (mpi_rank==0) print*, 'Warning ntype>npart, please correct ...'
          call stop_calculation
        end if

        chem%pname(ntype) = tname(1)
        if (mpi_rank==0) write(11,*) trim(chem%pname(ntype))

        if (trim(chem%pname(ntype)) == '[H-]') chem%tag_neg  = ntype
        if (trim(chem%pname(ntype)) == '[H]')  chem%tag_neu  = ntype
        if (trim(chem%pname(ntype)) == '[eb]') chem%tag_beam = ntype

        read(10,*,END=101,ERR=100) chem%mass(ntype), chem%charge(ntype), chem%Ti(ntype), ni0(ntype)

        if (mpi_rank==0 .and. ni0(ntype)>1.0_real64) then
          print*, 'Warning (n/n0)>1 for ', trim(chem%pname(ntype))
          print*, 'Please correct ...'
          call stop_calculation
        end if

        chem%charge(ntype) = chem%charge(ntype) * qe
        chem%mass(ntype)   = chem%mass(ntype)   * amu

        if (chem%charge(ntype) == 0.0_real64) then
          n_neu = n_neu + 1
          ni0(ntype) = ni0(ntype) * chem%ngas
          chem%mass(ntype) = -abs(chem%mass(ntype))
        end if

        goto 200
      end if

      ! -------------
      ! Reactions
      ! -------------
21    read(10,*,END=101,ERR=100) rlabel

      if (flag_re==0) then
        if (mpi_rank==0) then
          write(11,*) ' '
          write(11,*) 'Reactions:'
        end if
        flag_re = 1
      end if

      ind_label = 2
      l_sy = 0
      r_sy = 0

      do i_re=1,50
        if (rlabel(i_re:i_re) == '[') l_sy = i_re
        if (rlabel(i_re:i_re) == ']') r_sy = i_re
        if (rlabel(i_re:i_re) == '>') then
          n_re = ind_label - 2
          col_info(icol, ind_nre) = n_re
        end if

        if (r_sy > l_sy) then
          do ptype=1,(ntype+1)
            if (mpi_rank==0 .and. ptype==(ntype+1)) then
              write(*,500) icol
500           format(' Warning: unknown particle in reaction #',1x,i2)
              print*, 'Please correct ...'
              call stop_calculation
            end if
            if (rlabel(l_sy:r_sy) == chem%pname(ptype)) exit
          end do
          ind_label = ind_label + 1
          col_info(icol, ind_label) = ptype
          l_sy = 0
          r_sy = 0
        end if
      end do

      n_by = ind_label - n_re - 2

      col_info(icol, ind_nby) = n_by

      sig_scale = 0.0_real64
      read(10,*,END=101,ERR=100) sig_Eex(icol, ind_Eth), sig_Eex(icol, ind_dE)
      read(10,*,END=101,ERR=100) sig_scale(1), sig_scale(2)
      read(10,*,END=101,ERR=100) rtype

      sig_Eex(icol, ind_Eth) = sig_Eex(icol, ind_Eth) * sig_scale(1)

      if (trim(rtype)=='COLLISION'      .or. trim(rtype)=='collision')      col_info(icol, ind_nby+1+n_re+n_by) = 1
      if (trim(rtype)=='IONIZATION'     .or. trim(rtype)=='ionization')     col_info(icol, ind_nby+1+n_re+n_by) = 2
      if (trim(rtype)=='EXCITATION'     .or. trim(rtype)=='excitation')     col_info(icol, ind_nby+1+n_re+n_by) = 3
      if (trim(rtype)=='CHARGEEXCHANGE' .or. trim(rtype)=='chargeexchange') col_info(icol, ind_nby+1+n_re+n_by) = 4
      if (trim(rtype)=='DISSOCIATION'   .or. trim(rtype)=='dissociation')   col_info(icol, ind_nby+1+n_re+n_by) = 5

      if (mpi_rank==0 .and. col_info(icol, ind_nby+1+n_re+n_by) == 0) then
        print*, 'Warning: unknown reaction type in reaction #', icol
        print*, 'Please correct ...'
        call stop_calculation
      end if

      ! number of byproduct electrons
      do i_re = ind_nby+1+n_re, ind_nby+n_re+n_by
        if (col_info(icol, i_re) == 1) then
          col_info(icol, ind_nby+1+n_re+n_by+1) = col_info(icol, ind_nby+1+n_re+n_by+1) + 1
        end if
      end do

      do while (trim(name) /= '----')
        read(10,*,END=101,ERR=100) name
      end do

      if (mpi_rank==0) write(11,'(i3,1x,a30,1x,f8.2)') icol, rlabel, sig_Eex(icol, ind_Eth)

      do ipt=1,npt_mx
        if (mpi_rank==0 .and. ipt==npt_mx) then
          print*, 'Maximum number of data points in cross section reached'
          print*, 'Abort simulation ...'
          call stop_calculation
        end if

        read(10,fmt='(A)',advance='no') name
        if (trim(name) == '----') then
          read(10,fmt='(A)') name
          exit
        else
          backspace(10)
        end if

        sig_npt(icol) = sig_npt(icol) + 1
        read(10,*) sig_tmp(ipt,icol,1), sig_tmp(ipt,icol,2)
        sig_tmp(ipt,icol,1) = sig_tmp(ipt,icol,1) * sig_scale(1)
        sig_tmp(ipt,icol,2) = sig_tmp(ipt,icol,2) * sig_scale(2)
      end do

      goto 99

      ! ----------------
      ! Surface reactions
      ! ----------------
22    do while (trim(name) /= '----')
        read(10,*,END=101,ERR=100) name
      end do

      if (mpi_rank==0) then
        write(11,*) ' '
        write(11,*) 'Surface reactions:'
      end if

220   iscol = iscol + 1
      read(10,fmt='(A6)',END=101,ERR=100,advance='no') tname(1)

      if (tname(1)(1:2) == '--') then
        chem%nscol = iscol - 1
        call indexx(chem%nscol, sort_arr, indx)

        do iscol=1,chem%nscol
          ptype = int(scol_info(indx(iscol),1))
          chem%p_nscol(ptype) = chem%p_nscol(ptype) + 1
          scol_rank(ptype, chem%p_nscol(ptype)) = indx(iscol)
        end do
        goto 99
      else
        read(10,*,END=101,ERR=100) tname(2), scol_info(iscol,3), scol_info(iscol,4)

        if (mpi_rank==0) write(11,*) trim(tname(1)), trim(tname(2))

        do i=1,2
          do ptype=1,(ntype+1)
            if (mpi_rank==0 .and. ptype==(ntype+1)) then
              write(*,501) iscol
501           format(' Warning: unknown particle found in surface reaction #',1x,i2)
              print*, 'Please correct ...'
              call stop_calculation
            end if
            if (tname(i) == chem%pname(ptype)) exit
          end do
          scol_info(iscol,i) = real(ptype, real64)
        end do

        sort_arr(iscol) = scol_info(iscol,3)
        goto 220
      end if

99    continue
    end do

    if (mpi_rank==0) then
      print*, '# of collisions > maximum allowed, please correct'
      print*, 'Abort calculation ...'
      call stop_calculation
    end if

100 continue
    if (mpi_rank==0) then
      print*, 'An error occured while reading cross section #', icol
      print*, 'Abort calculation ...'
      call stop_calculation
    end if

101 continue
    if (mpi_rank==0 .and. flag_pa==0) then
      print*, 'Warning: "PARTICLE" section must be placed first, please correct ...'
      call stop_calculation
    end if

    if (mpi_rank==0) write(*,"(a)"), 'Gas chemistry read correctly (end of file reached)'

    chem%ncol  = icol
    chem%ntype = ntype
    chem%n_neu = n_neu

    if (chem%ngas == 0.0_real64) chem%ncol = 0

    if (mpi_rank == 0) write(*,*) " "
    if (mpi_rank == 0) write(*,"(a)") "Initializing chemistry/reactions..."

    if (chem%ncol == 0) then
      if (mpi_rank==0) print*, 'Collision-less mode is set'
    end if

    call ordering(chem, sig, sig_Er, sig_tmp, sig_npt, sig_list, col_info, sig_Eex, sigv_mx, ncol_mx, npt_mx, ntype, mpi_rank)

    deallocate(sig_tmp)

    ! Wall probability check
    do ptype=1,ntype
      Ps = 0.0_real64
      do icol=1,chem%p_nscol(ptype)
        Ps = Ps + scol_info(scol_rank(ptype,icol),3)
      end do
      if (mpi_rank==0 .and. Ps > 1.0_real64) then
        print*, 'Total wall collision probability greater than 1 for ', trim(chem%pname(ptype))
        print*, 'Please correct ...'
        call stop_calculation
      end if
    end do

    chem%ni0(1:chem%npart) = ni0(1:chem%npart)

    close(10)
    close(11)
  end subroutine read_reactions


  subroutine ordering(chem, sig, sig_Er, sig_tmp, sig_npt, sig_list, col_info, &
                      sig_Eex, sigv_mx, ncol_mx, npt_mx, ntype, mpi_rank)
    type(ChemistryState), intent(inout) :: chem
    integer, intent(in) :: ncol_mx, npt_mx, ntype, mpi_rank

    real(real64), intent(inout) :: sig_Er(:)
    real(real64), intent(inout) :: sig(:,:)
    real(real64), intent(inout) :: sig_tmp(:,:,:)
    integer,      intent(in)    :: sig_npt(:)
    integer,      intent(inout) :: sig_list(:,:)
    integer,      intent(in)    :: col_info(:,:)
    real(real64), intent(in)    :: sig_Eex(:,:)
    real(real64), intent(inout) :: sigv_mx(:,:)

    integer :: ipt, ipt_new, npt, icol, ptype, length_sort
    integer :: rtype1, rtype2
    real(real64) :: mu, vr, Ek_L, Ek_R, sig_L, sig_R, sig_p, Ekp, sort_Ek_sav, Eth
    real(real64), allocatable :: sort_Ek(:), sort_Ek_tmp(:)
    integer,      allocatable :: indx(:)
    integer,      allocatable :: ind_col(:)
    character(len=1) :: pnum

    sigv_mx  = 0.0_real64
    sig_list = 0

    allocate(ind_col(ntype))
    ind_col = 0

    length_sort = sum(sig_npt(1:chem%ncol))
    if (length_sort <= 0) then
      chem%sig_npt_mx = 0
      return
    end if

    allocate(sort_Ek(length_sort), sort_Ek_tmp(length_sort), indx(length_sort))

    ipt_new = 0
    do icol=1,chem%ncol
      do ipt=1,sig_npt(icol)
        ipt_new = ipt_new + 1
        sort_Ek(ipt_new) = sig_tmp(ipt,icol,1)
      end do
    end do

    call indexx(length_sort, sort_Ek, indx)

    ipt_new = 0
    sort_Ek_sav = 1.0d10
    do ipt=1,length_sort
      if (sort_Ek(indx(ipt)) /= sort_Ek_sav) then
        ipt_new = ipt_new + 1
        sort_Ek_tmp(ipt_new) = sort_Ek(indx(ipt))
        sort_Ek_sav = sort_Ek_tmp(ipt_new)
      end if
    end do

    chem%sig_npt_mx = ipt_new

    if (mpi_rank==0 .and. chem%sig_npt_mx >= npt_mx) then
      print*, 'sig_npt_mx > npt_mx, please correct'
      print*, 'sig_npt_mx=', chem%sig_npt_mx
      print*, 'Abort calculation ...'
      call stop_calculation
    end if

    sig_Er(1:chem%sig_npt_mx) = sort_Ek_tmp(1:chem%sig_npt_mx)

    deallocate(sort_Ek, sort_Ek_tmp)

    ! Interpolate sig onto the common energy grid
    do icol=1,chem%ncol
      npt = sig_npt(icol)
      if (npt < 2) cycle
      Eth = sig_Eex(icol, ind_Eth)

      do ipt=1,chem%sig_npt_mx
        Ekp = sig_Er(ipt)

        do ipt_new=1,npt-1
          Ek_L  = sig_tmp(ipt_new,icol,1)
          Ek_R  = sig_tmp(ipt_new+1,icol,1)
          sig_L = sig_tmp(ipt_new,icol,2)
          sig_R = sig_tmp(ipt_new+1,icol,2)

          if (Ekp < Ek_L) then
            sig_p = sig_L
            exit
          end if

          if (Ekp >= Ek_L .and. Ekp <= Ek_R) then
            sig_p = sig_L + (Ekp - Ek_L) * (sig_R - sig_L) / (Ek_R - Ek_L)
            exit
          end if
        end do

        if (sig_L==0.0_real64 .and. sig_R>0.0_real64) then
          if (Eth > Ek_R) then
            print'(" Eth greater than threshold energy in reaction #",1x,i3)', icol
            print*, 'Please correct...'
            call stop_calculation
          end if
        end if

        sig(ipt,icol) = sig_p
      end do
    end do

    ! Build per-ptype reaction lists and max(sig*v)
    do icol=1,chem%ncol
      rtype1 = col_info(icol, ind_nby+1)
      ind_col(rtype1) = ind_col(rtype1) + 1
      sig_list(rtype1, ind_col(rtype1)) = icol

      do ipt=1,chem%sig_npt_mx
        rtype1 = col_info(icol, ind_nby+1)
        rtype2 = col_info(icol, ind_nby+2)

        mu = abs(chem%mass(rtype1))*abs(chem%mass(rtype2)) / (abs(chem%mass(rtype1))+abs(chem%mass(rtype2)))
        vr = sqrt(2.0_real64 * sig_Er(ipt) * qe / mu)

        if (sig_Er(ipt) >= 1.0_real64) then
          sigv_mx(rtype1,icol) = max(sigv_mx(rtype1,icol), vr * sig(ipt,icol))
        end if
      end do
    end do

    chem%p_ncol(1:ntype) = ind_col(1:ntype)

    ! Dump reactions for debugging/plots (legacy behavior)
    if (mpi_rank==0) then
      do ptype=1,ntype
        write(pnum,'(i1)') ptype
        open(10,file='./Output/reactions.'//pnum)
        do icol=1,chem%p_ncol(ptype)
          do ipt=1,chem%sig_npt_mx
            write(10,'(1x,f12.4,1x,es12.4)') sig_Er(ipt), sig(ipt, sig_list(ptype,icol))
          end do
          write(10,*) ' '
        end do
        close(10)
      end do
    end if

    deallocate(indx, ind_col)
  end subroutine ordering

end module mod_reactions
