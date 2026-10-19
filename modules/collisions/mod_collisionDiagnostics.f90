module mod_collisionDiagnostics
  use iso_fortran_env, only: int32, int64, real64
  implicit none
  private
  public :: init_rxn_counts, count_rxn, print_rxn_counts, reset_rxn_counts
  public :: init_debug_diagnostics, record_np_mx, record_target_search, record_npt
  public :: print_debug_diagnostics, reset_debug_diagnostics

  integer(int64), allocatable, save :: rxn_count(:)

  ! --- per-tracked-species debug accumulators (mod_collisionsGwenael
  !     instrumentation, see that module's header). All reset on the same
  !     cadence as rxn_count via reset_debug_diagnostics().
  real(real64),   allocatable, save :: dbg_np_mx(:)           ! latest np_mx(ttype) snapshot
  integer(int64), allocatable, save :: dbg_search_ok(:)        ! find_charged_target successes
  integer(int64), allocatable, save :: dbg_search_fail(:)      ! find_charged_target failures
  integer(int64), allocatable, save :: dbg_npt_count(:)        ! number of np_t samples taken
  real(real64),   allocatable, save :: dbg_npt_sum(:)          ! running sum of np_t samples
  real(real64),   allocatable, save :: dbg_npt_max(:)          ! max np_t sample seen
  integer(int64), allocatable, save :: dbg_npt_exceeds_npmx(:) ! count of np_t > np_mx(ttype)

contains

  subroutine init_rxn_counts(ncol)
    integer(int32), intent(in) :: ncol

    if (.not. allocated(rxn_count)) then
      allocate(rxn_count(ncol))
      rxn_count = 0_int64
    end if
  end subroutine init_rxn_counts

  subroutine count_rxn(c_ind)
    integer(int32), intent(in) :: c_ind
    if (.not. allocated(rxn_count)) return
    if (c_ind < 1 .or. c_ind > size(rxn_count)) return
    !$omp atomic
    rxn_count(c_ind) = rxn_count(c_ind) + 1_int64
  end subroutine count_rxn

  subroutine print_rxn_counts()
    integer :: i
    if (.not. allocated(rxn_count)) return
    write(*,'(a)') " Accepted collision reactions:"
    do i = 1, size(rxn_count)
      if (rxn_count(i) > 0_int64) then
        write(*,'(a,i4,a,i12)') " reaction ", i, " accepted = ", rxn_count(i)
      end if
    end do
  end subroutine print_rxn_counts

  subroutine reset_rxn_counts()
    if (allocated(rxn_count)) rxn_count = 0_int64
  end subroutine reset_rxn_counts

  !=========================================================================
  ! Debug instrumentation for mod_collisionsGwenael, tracking H2+ (and
  ! every other tracked species) target-density behaviour so a growing
  ! new-vs-legacy discrepancy can be localized empirically instead of by
  ! code-reading alone.
  !=========================================================================

  subroutine init_debug_diagnostics(ntype_tracked)
    integer(int32), intent(in) :: ntype_tracked
    if (.not. allocated(dbg_np_mx)) then
      allocate(dbg_np_mx(ntype_tracked));            dbg_np_mx = 0.0_real64
      allocate(dbg_search_ok(ntype_tracked));        dbg_search_ok = 0_int64
      allocate(dbg_search_fail(ntype_tracked));      dbg_search_fail = 0_int64
      allocate(dbg_npt_count(ntype_tracked));        dbg_npt_count = 0_int64
      allocate(dbg_npt_sum(ntype_tracked));          dbg_npt_sum = 0.0_real64
      allocate(dbg_npt_max(ntype_tracked));          dbg_npt_max = 0.0_real64
      allocate(dbg_npt_exceeds_npmx(ntype_tracked)); dbg_npt_exceeds_npmx = 0_int64
    end if
  end subroutine init_debug_diagnostics

  subroutine record_np_mx(ttype, val)
    integer(int32), intent(in) :: ttype
    real(real64),   intent(in) :: val
    if (.not. allocated(dbg_np_mx)) return
    if (ttype < 1 .or. ttype > size(dbg_np_mx)) return
    dbg_np_mx(ttype) = val   ! latest snapshot, not accumulated - written once per ttype per call, no race
  end subroutine record_np_mx

  subroutine record_target_search(ttype, success)
    integer(int32), intent(in) :: ttype
    logical,        intent(in) :: success
    if (.not. allocated(dbg_search_ok)) return
    if (ttype < 1 .or. ttype > size(dbg_search_ok)) return
    if (success) then
      !$omp atomic
      dbg_search_ok(ttype) = dbg_search_ok(ttype) + 1_int64
    else
      !$omp atomic
      dbg_search_fail(ttype) = dbg_search_fail(ttype) + 1_int64
    end if
  end subroutine record_target_search

  subroutine record_npt(ttype, npt, npmx)
    integer(int32), intent(in) :: ttype
    real(real64),   intent(in) :: npt, npmx
    if (.not. allocated(dbg_npt_count)) return
    if (ttype < 1 .or. ttype > size(dbg_npt_count)) return
    !$omp atomic
    dbg_npt_count(ttype) = dbg_npt_count(ttype) + 1_int64
    !$omp atomic
    dbg_npt_sum(ttype) = dbg_npt_sum(ttype) + npt
    !$omp critical (dbg_npt_max_update)
    if (npt > dbg_npt_max(ttype)) dbg_npt_max(ttype) = npt
    !$omp end critical (dbg_npt_max_update)
    if (npt > npmx) then
      !$omp atomic
      dbg_npt_exceeds_npmx(ttype) = dbg_npt_exceeds_npmx(ttype) + 1_int64
    end if
  end subroutine record_npt

  subroutine print_debug_diagnostics()
    integer :: s
    real(real64) :: avg
    if (.not. allocated(dbg_np_mx)) return
    write(*,'(a)') " ----- mod_collisionsGwenael debug diagnostics -----"
    do s = 1, size(dbg_np_mx)
      if (dbg_search_ok(s) == 0_int64 .and. dbg_search_fail(s) == 0_int64) cycle
      avg = 0.0_real64
      if (dbg_npt_count(s) > 0_int64) avg = dbg_npt_sum(s) / real(dbg_npt_count(s), real64)
      write(*,'(a,i3,a,es12.4,a,i10,a,i10,a,es12.4,a,es12.4,a,i10)') &
        " ttype=", s, &
        " np_mx=", dbg_np_mx(s), &
        " search_ok=", dbg_search_ok(s), &
        " search_fail=", dbg_search_fail(s), &
        " np_t_avg=", avg, &
        " np_t_max=", dbg_npt_max(s), &
        " np_t>np_mx count=", dbg_npt_exceeds_npmx(s)
    end do
    write(*,'(a)') " ----------------------------------------------------"
  end subroutine print_debug_diagnostics

  subroutine reset_debug_diagnostics()
    if (.not. allocated(dbg_np_mx)) return
    dbg_search_ok = 0_int64
    dbg_search_fail = 0_int64
    dbg_npt_count = 0_int64
    dbg_npt_sum = 0.0_real64
    dbg_npt_max = 0.0_real64
    dbg_npt_exceeds_npmx = 0_int64
    ! dbg_np_mx is a snapshot, not accumulated - left alone, will just be
    ! overwritten by the next call to record_np_mx
  end subroutine reset_debug_diagnostics

end module mod_collisionDiagnostics