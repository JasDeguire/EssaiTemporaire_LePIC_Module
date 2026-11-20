module mod_restart
  !===========================================================================
  ! Port of legacy Src/restart.f90 + the backup-write block of Src/main.f90
  ! (the `it.gt.1 .and. MOD(it,nbak).eq.1` section) onto the modular
  ! ParticleSet architecture.
  !
  ! Two pieces:
  !   - write_restart_modular : periodic backup write, one unformatted file
  !     per MPI rank, alternating Backups/DATA.BAK/ and Backups/DATA.BAK2/ exactly like
  !     legacy (so a crash mid-write never destroys the only valid backup).
  !   - restart_particles_modular : the read-back side, supporting both
  !     flag_restart=1 (same MPI/OMP rank count as the run that wrote the
  !     files) and flag_restart=2 ("external" files, written by a run with
  !     a different rank count, redistributed round-robin across the
  !     current ranks - same algorithm as legacy's restart()), followed by
  !     np_dup macroparticle duplication (np_dup>1) or random thinning
  !     (0<np_dup<1).
  !
  ! File format: matches legacy's particles*.bak layout exactly (time,
  ! np_tot(ntype,nproc), then per iproc per ptype a vxp(6,n) block, then
  ! flag_cex for the negative-ion species if present) WITH ONE DELIBERATE
  ! ADDITION: a leading `nproc_written` integer.
  !
  ! Why: in legacy, the "current" thread count (nproc, from
  ! omp_get_max_threads()) and the "original" thread count used only for
  ! restart redistribution math (omp_rank_max, read from the input file)
  ! are two genuinely independent variables. In this modular codebase,
  ! mod_state.f90 sets self%nproc DIRECTLY from cfg%omp_rank_max (the same
  ! input-file field) - so that field can no longer carry "the original
  ! run's thread count" separately from "this run's thread count" the way
  ! legacy intends for a flag_restart=2 cross-config restart. Recording
  ! the writer's own nproc directly in the file header is the minimal,
  ! self-contained way to preserve full restart-with-different-rank-count
  ! capability without changing the input file format or the existing
  ! cfg%omp_rank_max -> self%nproc convention used everywhere else.
  ! cfg%mpi_rank_max is NOT affected by this and is used exactly as
  ! legacy uses it (self%mpi_size comes from the real MPI runtime, not
  ! from cfg%mpi_rank_max, so that field is never repurposed).
  !
  ! Not ported: legacy's sum_dEk/Nh re-accumulation inside restart_OMP
  ! (heating-region energy bookkeeping for the adaptive electron
  ! temperature). This is a deliberate scope reduction - those
  ! accumulators are naturally rebuilt by the heating routine over the
  ! following steps regardless, so skipping their one-time restart
  ! recompute only affects the very first few heating updates after a
  ! restart, not the long-run physics.
  !
  ! SECOND DELIBERATE ADDITION: when flag_RFant==1, a scalar E0_RF record
  ! immediately after `time`, matching legacy (Src/restart.f90:63,87),
  ! which also stores time+E0_RF and nothing else RF-specific. Gated on
  ! flag_RFant so the file format is byte-for-byte unchanged for every
  ! non-RF run (the vast majority, and the only ones with any pre-existing
  ! local backups to stay compatible with - RF antenna heating didn't
  ! exist in this modular codebase before this record was added, so there
  ! is no old-format RF backup to preserve compatibility with). gams (skin
  ! depth) is deliberately NOT persisted, here or in legacy: it's fully
  ! re-derived from the local electron density on the first
  ! update_rf_convergence call after restart, and that density is already
  ! reconstructed by this module's existing post-restart deposit below.
  !===========================================================================
  use iso_fortran_env,       only: int8, int32, real64
  use mod_particles,         only: ParticleSet
  use mod_config,            only: Config
  use mod_chargeDeposition,  only: deposit_particle_set_to_np_thread
  use mod_RNG,                only: ran2
  implicit none
  private

  public :: write_restart_modular
  public :: restart_particles_modular

contains

  !=========================================================================
  ! WRITE side - one call per backup step, same cadence as legacy
  ! (it.gt.1 .and. MOD(it,nbak).eq.1), called from mod_simulation.f90.
  ! Alternates Backups/DATA.BAK/ <-> Backups/DATA.BAK2/ via the caller-owned flag_wrt.
  !=========================================================================
  subroutine write_restart_modular(mpi_rank, ntype, nproc, tag_neg, part, time, flag_wrt, &
                                    flag_RFant, E0_RF)
    integer,           intent(in)    :: mpi_rank
    integer(int32),    intent(in)    :: ntype, nproc
    integer(int32),    intent(in)    :: tag_neg
    type(ParticleSet), intent(in)    :: part(:,:)
    real(real64),      intent(in)    :: time
    integer(int32),    intent(inout) :: flag_wrt
    integer(int32),    intent(in)    :: flag_RFant
    real(real64),      intent(in)    :: E0_RF

    integer        :: unit_no, ptype, iproc, npar, ios
    character(len=64) :: fname
    character(len=16) :: rank_str
    real(real64), allocatable :: vxp_buf(:,:)
    integer(int32), allocatable :: np_tot(:,:)

    unit_no = 41 + mpi_rank
    write(rank_str, '(I0)') mpi_rank

    if (flag_wrt == 0_int32) then
      fname = 'Backups/DATA.BAK/particles' // trim(rank_str) // '.bak'
      flag_wrt = 1_int32
    else
      fname = 'Backups/DATA.BAK2/particles' // trim(rank_str) // '.bak'
      flag_wrt = 0_int32
    end if

    allocate(np_tot(ntype, nproc))
    do iproc = 1_int32, nproc
      do ptype = 1_int32, ntype
        np_tot(ptype, iproc) = part(ptype, iproc)%n
      end do
    end do

    open(unit_no, file=trim(fname), form='UNFORMATTED', status='REPLACE', iostat=ios)
    if (ios /= 0_int32) then
      write(*,*) 'write_restart_modular: failed to open ', trim(fname), ' (iostat=', ios, ')'
      write(*,*) 'Does the directory exist? Expected: Backups/DATA.BAK/ and Backups/DATA.BAK2/', &
                 ' relative to the run''s working directory.'
      error stop 'write_restart_modular: cannot open backup file for writing'
    end if

    write(unit_no) nproc       ! header: this writer's own thread count - see module header
    write(unit_no) time
    if (flag_RFant == 1_int32) write(unit_no) E0_RF
    write(unit_no) np_tot(1:ntype, 1:nproc)

    do iproc = 1_int32, nproc
      do ptype = 1_int32, ntype
        npar = np_tot(ptype, iproc)
        if (npar > 0_int32) then
          allocate(vxp_buf(6, npar))
          call pset_to_vxp(part(ptype, iproc), npar, vxp_buf)
          write(unit_no) vxp_buf
          deallocate(vxp_buf)
        else
          write(unit_no) reshape([real(real64) ::], [6, 0])
        end if
      end do

      if (tag_neg > 0_int32) then
        npar = np_tot(tag_neg, iproc)
        if (npar > 0_int32) then
          write(unit_no) part(tag_neg, iproc)%flag_cex(1:npar)
        else
          write(unit_no) [integer(int32) ::]
        end if
      end if
    end do

    close(unit_no)
    deallocate(np_tot)
  end subroutine write_restart_modular


  !=========================================================================
  ! READ side - called once from mod_state%init_particles when
  ! cfg%flag_restart > 0, in place of load_particles_modular.
  !=========================================================================
  subroutine restart_particles_modular( &
      cfg, mpi_rank, mpi_size, n, h, kq, Nm, ntype, nproc, tag_neg, &
      iseed, part, np_thread, time_out, E0_RF_out)

    type(Config),      intent(in)    :: cfg
    integer,           intent(in)    :: mpi_rank, mpi_size
    integer(int32),    intent(in)    :: n(3)
    real(real64),      intent(in)    :: h(3)
    real(real64),      intent(in)    :: kq(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),      intent(in)    :: Nm(:)
    integer(int32),    intent(in)    :: ntype, nproc
    integer(int32),    intent(in)    :: tag_neg
    integer(int32),    intent(inout) :: iseed(:)
    type(ParticleSet), intent(inout) :: part(:,:)
    real(real64),      intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
    real(real64),      intent(out)   :: time_out
    real(real64),      intent(out)   :: E0_RF_out

    integer(int32) :: ptype, iproc

    do iproc = 1_int32, nproc
      do ptype = 1_int32, ntype
        call part(ptype, iproc)%clear()
      end do
    end do

    if (cfg%flag_restart == 1_int32) then
      call read_same_rank_count(cfg, mpi_rank, ntype, nproc, tag_neg, part, time_out, E0_RF_out)
    else
      call read_external_ranks(cfg, mpi_rank, mpi_size, ntype, nproc, tag_neg, part, time_out, E0_RF_out)
    end if

    call apply_np_dup(cfg%np_dup, ntype, nproc, iseed, part)

    ! Density deposit (legacy: restart_OMP), reusing the existing
    ! production deposition kernel rather than re-deriving the CIC
    ! weights here.
    np_thread = 0.0_real64
    !$omp parallel do collapse(2) private(ptype,iproc) schedule(static)
    do ptype = 1_int32, ntype
      do iproc = 1_int32, nproc
        if (part(ptype, iproc)%n <= 0_int32) cycle
        call deposit_particle_set_to_np_thread( &
          part       = part(ptype, iproc), &
          n          = n, &
          h          = h, &
          kq         = kq, &
          Nm_species = Nm(ptype), &
          np_local   = np_thread(:,:,:,ptype,iproc) )
      end do
    end do
    !$omp end parallel do

  end subroutine restart_particles_modular


  !=========================================================================
  ! flag_restart == 1 : same MPI/OMP rank count as the writer. Each rank
  ! reads only its own Backups/DATA.BAK/particles<mpi_rank>.bak.
  !=========================================================================
  subroutine read_same_rank_count(cfg, mpi_rank, ntype, nproc, tag_neg, part, time_out, E0_RF_out)
    type(Config),      intent(in)    :: cfg
    integer,           intent(in)    :: mpi_rank
    integer(int32),    intent(in)    :: ntype, nproc
    integer(int32),    intent(in)    :: tag_neg
    type(ParticleSet), intent(inout) :: part(:,:)
    real(real64),      intent(out)   :: time_out
    real(real64),      intent(out)   :: E0_RF_out

    integer        :: unit_no, ptype, iproc, npar, nproc_written
    character(len=64) :: fname
    character(len=16) :: rank_str
    real(real64), allocatable :: vxp_buf(:,:)
    integer(int32), allocatable :: np_tot(:,:)
    integer(int32), allocatable :: cex_buf(:)

    unit_no = 41 + mpi_rank
    write(rank_str, '(I0)') mpi_rank
    fname = 'Backups/DATA.BAK/particles' // trim(rank_str) // '.bak'

    open(unit_no, file=trim(fname), form='UNFORMATTED', status='OLD', err=999)
    read(unit_no, err=999) nproc_written
    read(unit_no, err=999) time_out

    E0_RF_out = cfg%E0_RF
    if (cfg%flag_RFant == 1_int32) then
      read(unit_no, err=999) E0_RF_out
    end if

    if (nproc_written /= nproc) then
      write(*,*) 'restart_particles_modular: backup file was written with nproc=', &
                 nproc_written, ' but this run uses nproc=', nproc
      write(*,*) 'flag_restart=1 assumes the SAME rank count; use flag_restart=2', &
                 ' (external/Backups/SAV_DATA) for a different rank count.'
      error stop 'restart_particles_modular: rank-count mismatch for flag_restart=1'
    end if

    allocate(np_tot(ntype, nproc))
    read(unit_no, err=999) np_tot(1:ntype, 1:nproc)

    do iproc = 1_int32, nproc
      do ptype = 1_int32, ntype
        npar = np_tot(ptype, iproc)
        allocate(vxp_buf(6, max(npar,1_int32)))
        read(unit_no, err=999) vxp_buf(:, 1:npar)
        if (npar > 0_int32) call part(ptype, iproc)%from_vxp(vxp_buf(:,1:npar), npar, ptype)
        deallocate(vxp_buf)
      end do

      if (tag_neg > 0_int32) then
        npar = np_tot(tag_neg, iproc)
        allocate(cex_buf(max(npar,1_int32)))
        read(unit_no, err=999) cex_buf(1:npar)
        if (npar > 0_int32) part(tag_neg, iproc)%flag_cex(1:npar) = cex_buf(1:npar)
        deallocate(cex_buf)
      end if
    end do

    close(unit_no)
    deallocate(np_tot)
    return

999 continue
    write(*,*) 'restart_particles_modular: failed to read ', trim(fname)
    error stop 'restart_particles_modular: backup file read error (flag_restart=1)'
  end subroutine read_same_rank_count


  !=========================================================================
  ! flag_restart == 2 : "external" files, written by a run with a
  ! (possibly) different MPI rank count and/or OMP thread count, read
  ! from Backups/SAV_DATA/ and redistributed round-robin into the current
  ! part(ptype,iproc) arrays - same algorithm as legacy's restart().
  !=========================================================================
  subroutine read_external_ranks(cfg, mpi_rank, mpi_size, ntype, nproc, tag_neg, part, time_out, E0_RF_out)
    type(Config),      intent(in)    :: cfg
    integer,           intent(in)    :: mpi_rank, mpi_size
    integer(int32),    intent(in)    :: ntype, nproc
    integer(int32),    intent(in)    :: tag_neg
    type(ParticleSet), intent(inout) :: part(:,:)
    real(real64),      intent(out)   :: time_out
    real(real64),      intent(out)   :: E0_RF_out

    integer        :: unit_no, irank, mpi_rank_div
    integer        :: ptype, iproc_orig, iproc_tmp
    integer(int32) :: nproc_written
    integer(int32) :: n1, n2, npar
    character(len=64) :: fname
    character(len=16) :: rank_str
    real(real64), allocatable :: vxp_buf(:,:)
    integer(int32), allocatable :: np_tot_ext(:,:)
    integer(int32), allocatable :: cex_buf(:)
    integer(int32) :: n1_tag, n2_tag

    if (mod(int(cfg%mpi_rank_max,int32), int(mpi_size,int32)) /= 0_int32) then
      write(*,*) 'restart_particles_modular: mpi_rank_max=', cfg%mpi_rank_max, &
                 ' is not a multiple of the current MPI size=', mpi_size
      error stop 'restart_particles_modular: bad mpi_rank_max for flag_restart=2'
    end if
    mpi_rank_div = int(cfg%mpi_rank_max, int32) / mpi_size

    time_out  = 0.0_real64
    E0_RF_out = cfg%E0_RF
    unit_no = 41 + mpi_rank

    do irank = mpi_rank*mpi_rank_div, (mpi_rank+1)*mpi_rank_div - 1_int32

      write(rank_str, '(I0)') irank
      fname = 'Backups/SAV_DATA/particles' // trim(rank_str) // '.bak'
      write(*,*) 'Reading external file ', trim(fname), ', MPI rank=', mpi_rank

      open(unit_no, file=trim(fname), form='UNFORMATTED', status='OLD', err=998)
      read(unit_no, err=998) nproc_written
      read(unit_no, err=998) time_out
      if (cfg%flag_RFant == 1_int32) then
        read(unit_no, err=998) E0_RF_out
      end if

      allocate(np_tot_ext(ntype, nproc_written))
      read(unit_no, err=998) np_tot_ext(1:ntype, 1:nproc_written)

      iproc_tmp = 0
      do iproc_orig = 1, int(nproc_written)

        if (nproc_written == nproc) then
          iproc_tmp = iproc_orig
        else
          iproc_tmp = iproc_tmp + 1
          if (iproc_tmp > nproc) iproc_tmp = 1
        end if

        n1_tag = 0_int32
        n2_tag = 0_int32

        do ptype = 1_int32, ntype
          npar = np_tot_ext(ptype, iproc_orig)
          n1 = part(ptype, iproc_tmp)%n
          n2 = n1 + npar

          allocate(vxp_buf(6, max(npar,1_int32)))
          read(unit_no, err=998) vxp_buf(:, 1:npar)

          if (npar > 0_int32) then
            call part(ptype, iproc_tmp)%ensure_capacity(n2)
            call append_vxp(part(ptype, iproc_tmp), vxp_buf(:,1:npar), n1, npar, ptype)
          end if
          deallocate(vxp_buf)

          if (ptype == tag_neg) then
            n1_tag = n1
            n2_tag = n2
          end if
        end do

        if (tag_neg > 0_int32) then
          npar = n2_tag - n1_tag
          allocate(cex_buf(max(npar,1_int32)))
          read(unit_no, err=998) cex_buf(1:npar)
          if (npar > 0_int32) part(tag_neg, iproc_tmp)%flag_cex(n1_tag+1:n2_tag) = cex_buf(1:npar)
          deallocate(cex_buf)
        end if

      end do

      close(unit_no)
      deallocate(np_tot_ext)
    end do

    return

998 continue
    write(*,*) 'restart_particles_modular: failed to read external file ', trim(fname)
    error stop 'restart_particles_modular: backup file read error (flag_restart=2)'
  end subroutine read_external_ranks


  !=========================================================================
  ! np_dup macroparticle duplication (np_dup > 1, integer factor, exact
  ! copies appended) or random thinning (0 < np_dup < 1, each particle
  ! kept independently with probability np_dup) - same as legacy.
  !=========================================================================
  subroutine apply_np_dup(np_dup, ntype, nproc, iseed, part)
    real(real64),      intent(in)    :: np_dup
    integer(int32),    intent(in)    :: ntype, nproc
    integer(int32),    intent(inout) :: iseed(:)
    type(ParticleSet), intent(inout) :: part(:,:)

    integer(int32) :: ptype, iproc, i, k, n0, n_keep
    real(real64)   :: rnd

    if (abs(np_dup - 1.0_real64) < 1.0e-12_real64) return

    if (abs(np_dup) > 1.0_real64) then
      ! --- duplicate: append (round(|np_dup|)-1) extra copies of the
      !     originally-loaded range, exactly like legacy.
      do iproc = 1_int32, nproc
        do ptype = 1_int32, ntype
          n0 = part(ptype, iproc)%n
          if (n0 <= 0_int32) cycle
          call part(ptype, iproc)%ensure_capacity(n0 * nint(abs(np_dup)))
          do k = 1_int32, nint(abs(np_dup)) - 1_int32
            part(ptype,iproc)%pv(1,k*n0+1:(k+1)*n0)  = part(ptype,iproc)%pv(1,1:n0)
            part(ptype,iproc)%pv(2,k*n0+1:(k+1)*n0)  = part(ptype,iproc)%pv(2,1:n0)
            part(ptype,iproc)%pv(3,k*n0+1:(k+1)*n0)  = part(ptype,iproc)%pv(3,1:n0)
            part(ptype,iproc)%pv(4,k*n0+1:(k+1)*n0) = part(ptype,iproc)%pv(4,1:n0)
            part(ptype,iproc)%pv(5,k*n0+1:(k+1)*n0) = part(ptype,iproc)%pv(5,1:n0)
            part(ptype,iproc)%pv(6,k*n0+1:(k+1)*n0) = part(ptype,iproc)%pv(6,1:n0)
            part(ptype,iproc)%w(k*n0+1:(k+1)*n0)  = part(ptype,iproc)%w(1:n0)
            part(ptype,iproc)%sp(k*n0+1:(k+1)*n0) = part(ptype,iproc)%sp(1:n0)
            part(ptype,iproc)%flag_dead(k*n0+1:(k+1)*n0) = part(ptype,iproc)%flag_dead(1:n0)
            part(ptype,iproc)%flag_cex(k*n0+1:(k+1)*n0)  = part(ptype,iproc)%flag_cex(1:n0)
          end do
          part(ptype, iproc)%n = n0 * nint(abs(np_dup))
        end do
      end do
    else
      ! --- thin: keep each particle independently with probability np_dup,
      !     compacting in place (same ran2 stream convention as legacy:
      !     one thread/iproc's own iseed entry, drawn in particle order).
      do iproc = 1_int32, nproc
        do ptype = 1_int32, ntype
          n0 = part(ptype, iproc)%n
          if (n0 <= 0_int32) cycle
          n_keep = 0_int32
          do i = 1_int32, n0
            rnd = ran2(iseed(iproc))
            if (rnd <= abs(np_dup)) then
              n_keep = n_keep + 1_int32
              if (n_keep /= i) then
                part(ptype,iproc)%pv(1,n_keep)  = part(ptype,iproc)%pv(1,i)
                part(ptype,iproc)%pv(2,n_keep)  = part(ptype,iproc)%pv(2,i)
                part(ptype,iproc)%pv(3,n_keep)  = part(ptype,iproc)%pv(3,i)
                part(ptype,iproc)%pv(4,n_keep) = part(ptype,iproc)%pv(4,i)
                part(ptype,iproc)%pv(5,n_keep) = part(ptype,iproc)%pv(5,i)
                part(ptype,iproc)%pv(6,n_keep) = part(ptype,iproc)%pv(6,i)
                part(ptype,iproc)%w(n_keep)  = part(ptype,iproc)%w(i)
                part(ptype,iproc)%sp(n_keep) = part(ptype,iproc)%sp(i)
                part(ptype,iproc)%flag_dead(n_keep) = part(ptype,iproc)%flag_dead(i)
                part(ptype,iproc)%flag_cex(n_keep)  = part(ptype,iproc)%flag_cex(i)
              end if
            end if
          end do
          part(ptype, iproc)%n = n_keep
        end do
      end do
    end if
  end subroutine apply_np_dup


  !=========================================================================
  ! Small helpers
  !=========================================================================
  subroutine pset_to_vxp(part, npar, vxp_buf)
    type(ParticleSet), intent(in)  :: part
    integer(int32),    intent(in)  :: npar
    real(real64),      intent(out) :: vxp_buf(6, npar)
    integer(int32) :: i
    do i = 1_int32, npar
      vxp_buf(1,i) = part%pv(1,i)
      vxp_buf(2,i) = part%pv(2,i)
      vxp_buf(3,i) = part%pv(3,i)
      vxp_buf(4,i) = part%pv(4,i)
      vxp_buf(5,i) = part%pv(5,i)
      vxp_buf(6,i) = part%pv(6,i)
    end do
  end subroutine pset_to_vxp

  ! Append npar particles from vxp_buf into part, starting right after the
  ! n1 particles already present (part%n is left at n1+npar on exit).
  ! Caller must have already called ensure_capacity(n1+npar).
  subroutine append_vxp(part, vxp_buf, n1, npar, species_id)
    type(ParticleSet), intent(inout) :: part
    real(real64),      intent(in)    :: vxp_buf(6, npar)
    integer(int32),    intent(in)    :: n1, npar, species_id
    integer(int32) :: i
    do i = 1_int32, npar
      part%pv(1,n1+i)  = vxp_buf(1,i)
      part%pv(2,n1+i)  = vxp_buf(2,i)
      part%pv(3,n1+i)  = vxp_buf(3,i)
      part%pv(4,n1+i) = vxp_buf(4,i)
      part%pv(5,n1+i) = vxp_buf(5,i)
      part%pv(6,n1+i) = vxp_buf(6,i)
      part%w(n1+i)  = 1.0_real64
      part%sp(n1+i) = species_id
      part%flag_dead(n1+i) = 0_int8
      part%flag_cex(n1+i)  = 0_int32
    end do
    part%n = n1 + npar
  end subroutine append_vxp

end module mod_restart