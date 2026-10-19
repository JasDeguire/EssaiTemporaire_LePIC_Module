module mod_PoissonSolver_legacy
  use iso_fortran_env, only: real64, int32
  use mpi
  use mod_poisson_decomp, only: PoissonDecomp
  use mod_part_info, only: flag_pbc, flag_nmn
  implicit none
  private
  public :: solve_poisson_legacy
  public :: print_poisson_breakdown, reset_poisson_breakdown

  ! --- internal timing breakdown, to localize where time actually goes
  ! inside solve_poisson_legacy (scatter to local decomposition, the
  ! pdesolver multigrid call itself, gather back to global) - added
  ! because eliminating the per-call temp-array allocation/copy did not
  ! move the measured "poisson (ms)" total, so the cost must be one of
  ! these three, not the part that was changed.
  real(real64), save :: t_scatter = 0.0_real64
  real(real64), save :: t_solve   = 0.0_real64
  real(real64), save :: t_gather  = 0.0_real64
  real(real64), save :: sum_k_it  = 0.0_real64
  real(real64), save :: sum_res   = 0.0_real64
  real(real64), save :: max_k_it  = 0.0_real64

  interface
    subroutine pdesolver(u,b,bcnd,h,n,ncycl,eps,omega,k,ktot,res,ng,rank,nproc)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in)    :: n(3), ng, ncycl, rank, nproc
      integer, intent(inout) :: k
      real(real64), intent(inout) :: h(3)
      real(real64), intent(inout) :: u(0:n(1)+2,0:n(2)+2,-1:n(3)/nproc+2)
      real(real64), intent(inout) :: b(0:n(1)+1,0:n(2)+1,0:n(3)/nproc+1)
      integer,      intent(inout) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)/nproc+2)
      real(real64), intent(in)    :: eps, omega
      real(real64), intent(inout) :: ktot, res
    end subroutine pdesolver
  end interface

contains

  subroutine solve_poisson_legacy(pdec, phi_global, bcnd_global, rhs_global, h, n_in, &
                                  ncycl, eps, omega, ng, flag_pbc_in, flag_nmn_in)
    type(PoissonDecomp), intent(inout) :: pdec
    real(real64), intent(inout)        :: phi_global(0:,0:,0:)
    integer(int32), intent(in)         :: bcnd_global(0:,0:,0:)
    real(real64), intent(in)           :: rhs_global(0:,0:,0:)
    real(real64), intent(in)           :: h(3)
    integer(int32), intent(in)         :: n_in(3)
    integer, intent(in)                :: ncycl, ng
    real(real64), intent(in)           :: eps, omega
    integer, intent(in)                :: flag_pbc_in, flag_nmn_in

    integer :: k_it
    real(real64) :: ktot, res
    integer :: n_leg(3)
    real(real64) :: h_leg(3)
    real(real64) :: tA, tB

    flag_pbc = flag_pbc_in
    flag_nmn = flag_nmn_in

    tA = MPI_Wtime()
    call pdec%scatter_from_global(phi_global, bcnd_global)
    call pdec%scatter_rhs_from_global(rhs_global)
    tB = MPI_Wtime()
    t_scatter = t_scatter + (tB - tA)

    n_leg = int(n_in, kind=4)
    h_leg = h

    ! Legacy zeroed the whole u_mg array (including the extra z=-1 ghost
    ! layer) before every solve; with phi_dom now persistent across calls
    ! (see mod_PoissonDecomposition%init), only that one 2D ghost layer
    ! needs resetting here - everything else (phi_dom(:,:,0:m+2),
    ! rhs_dom, bcnd_dom) was just freshly written by the scatter calls
    ! above, so no separate copy into a temporary is needed at all.
    pdec%phi_dom(:,:,-1) = 0.0_real64

    k_it = 0
    ktot = 0.0_real64
    res  = 0.0_real64

    tA = MPI_Wtime()
    call pdesolver(pdec%phi_dom, pdec%rhs_dom, pdec%bcnd_dom, h_leg, n_leg, &
                   ncycl, eps, omega, k_it, ktot, res, ng, pdec%rank, pdec%nproc)
    tB = MPI_Wtime()
    t_solve = t_solve + (tB - tA)
    sum_k_it = sum_k_it + real(k_it, real64)
    sum_res  = sum_res + res
    max_k_it = max(max_k_it, real(k_it, real64))

    tA = MPI_Wtime()
    call pdec%gather_phi_to_global(phi_global)
    tB = MPI_Wtime()
    t_gather = t_gather + (tB - tA)
  end subroutine solve_poisson_legacy


  subroutine print_poisson_breakdown(t_count)
    integer(int32), intent(in) :: t_count
    real(real64) :: denom
    if (t_count <= 0_int32) return
    denom = real(t_count, real64)
    write(*,'(a)') " ----- poisson timing breakdown -----"
    write(*,'(a,f8.2,a,f8.2,a,f8.2)') &
      "  scatter(ms)=", 1000.0_real64*t_scatter/denom, &
      "  solve(ms)=",   1000.0_real64*t_solve/denom, &
      "  gather(ms)=",  1000.0_real64*t_gather/denom
    write(*,'(a,f8.2,a,f8.2,a,es10.2)') &
      "  k_it_avg=", sum_k_it/denom, &
      "  k_it_max=", max_k_it, &
      "  res_avg=", sum_res/denom
    write(*,'(a)') " -------------------------------------"
  end subroutine print_poisson_breakdown


  subroutine reset_poisson_breakdown()
    t_scatter = 0.0_real64
    t_solve   = 0.0_real64
    t_gather  = 0.0_real64
    sum_k_it  = 0.0_real64
    sum_res   = 0.0_real64
    max_k_it  = 0.0_real64
  end subroutine reset_poisson_breakdown

end module mod_PoissonSolver_legacy