module mod_poisson_decomp
  use iso_fortran_env, only: int32, real64
  use mpi
  implicit none
  private
  public :: PoissonDecomp

  type :: PoissonDecomp
    ! Global sizes
    integer(int32) :: nx=0, ny=0, nz=0

    ! MPI
    integer :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nproc = -1

    ! Z-slab decomposition info
    integer(int32) :: m = 0
    integer(int32) :: k0 = 0
    integer(int32) :: k1 = 0
    integer(int32) :: kl = 0
    integer(int32) :: kr = 0

    ! Allgatherv bookkeeping
    integer, allocatable :: recvcounts_phi(:), displs_phi(:)
    integer, allocatable :: recvcounts_rhs(:), displs_rhs(:)

    ! Local arrays
    real(real64),   allocatable :: phi_dom(:,:,:)
    real(real64),   allocatable :: rhs_dom(:,:,:)
    integer(int32), allocatable :: bcnd_dom(:,:,:)

  contains
    procedure :: init
    procedure :: destroy
    procedure :: scatter_from_global
    procedure :: gather_phi_to_global
    procedure :: gather_rhs_to_global
    procedure :: scatter_rhs_from_global
  end type PoissonDecomp

contains

  subroutine init(self, nx, ny, nz, comm)
    class(PoissonDecomp), intent(inout) :: self
    integer(int32), intent(in) :: nx, ny, nz
    integer,        intent(in) :: comm
    integer :: ierr, r
    integer(int32) :: planes_r
    integer :: off

    self%nx = nx; self%ny = ny; self%nz = nz
    self%comm = comm
    call MPI_Comm_rank(comm, self%rank, ierr)
    call MPI_Comm_size(comm, self%nproc, ierr)

    if (mod(nz, self%nproc) /= 0) then
      if (self%rank == 0) write(*,*) "ERROR: nz not divisible by nproc (legacy requires this)."
      error stop "PoissonDecomp%init: nz % nproc != 0"
    end if

    self%m  = nz / self%nproc
    self%k0 = self%rank * self%m
    self%k1 = self%k0 + self%m - 1

    ! phi_dom's z range starts at -1, not 0: pdesolver (legacy multigrid)
    ! expects one extra ghost layer below the local domain
    ! (u(0:n1+2,0:n2+2,-1:m+2)). Allocating it here lets
    ! mod_PoissonSolver_legacy pass phi_dom straight through with no
    ! temporary copy.
    allocate(self%phi_dom(0:nx+2, 0:ny+2, -1:self%m+2))
    allocate(self%rhs_dom(0:nx+1, 0:ny+1, 0:self%m+1))
    allocate(self%bcnd_dom(0:nx+2, 0:ny+2, 0:self%m+2))

    self%phi_dom  = 0.0_real64
    self%rhs_dom  = 0.0_real64
    self%bcnd_dom = 0_int32

    allocate(self%recvcounts_phi(0:self%nproc-1), self%displs_phi(0:self%nproc-1))
    allocate(self%recvcounts_rhs(0:self%nproc-1), self%displs_rhs(0:self%nproc-1))

    ! ---- phi counts/displs (nx+3)*(ny+3) per plane ----
    off = 0
    do r=0, self%nproc-1
      planes_r = self%m
      if (r == 0)            planes_r = planes_r + 1
      if (r == self%nproc-1) planes_r = planes_r + 2

      self%recvcounts_phi(r) = int(planes_r, kind(off)) * (self%nx+3) * (self%ny+3)
      self%displs_phi(r)     = off
      off = off + self%recvcounts_phi(r)
    end do

    ! ---- rhs counts/displs (nx+2)*(ny+2) per plane ----
    off = 0
    do r=0, self%nproc-1
      planes_r = self%m
      if (r == 0)            planes_r = planes_r + 1   ! include global z=0
      if (r == self%nproc-1) planes_r = planes_r + 1   ! include global z=nz+1

      self%recvcounts_rhs(r) = int(planes_r, kind(off)) * (self%nx+2) * (self%ny+2)
      self%displs_rhs(r)     = off
      off = off + self%recvcounts_rhs(r)
    end do

    ! Send range (local indices) used in gather
    self%kl = 1
    if (self%rank == 0) self%kl = 0

    self%kr = self%m
    if (self%rank == self%nproc-1) self%kr = self%m + 2
  end subroutine init


  subroutine scatter_from_global(self, phi, bcnd)
    class(PoissonDecomp), intent(inout) :: self
    real(real64), intent(in) :: phi(0:,0:,0:)
    integer,      intent(in) :: bcnd(0:,0:,0:)

    integer :: iz0, iz1, k

    iz0 = self%k0
    iz1 = self%k0 + self%m + 2

    if (iz0 < lbound(phi,3) .or. iz1 > ubound(phi,3)) then
      write(*,*) 'scatter_from_global bounds error'
      write(*,*) 'iz0, iz1 = ', iz0, iz1
      write(*,*) 'phi z bounds = ', lbound(phi,3), ubound(phi,3)
      error stop
    end if

    ! Plane-parallel copies (were single-threaded array assignments).
    !$omp parallel do private(k) schedule(static)
    do k = 0, self%m+2
      self%phi_dom(:,:,k)  = phi(:,:,iz0+k)
      self%bcnd_dom(:,:,k) = int(bcnd(:,:,iz0+k), int32)
    end do
    !$omp end parallel do
  end subroutine scatter_from_global




  subroutine scatter_rhs_from_global(self, rhs)
    class(PoissonDecomp), intent(inout) :: self
    real(real64), intent(in) :: rhs(0:,0:,0:)

    integer :: iz0, iz1, k

    iz0 = self%k0
    iz1 = self%k0 + self%m + 1

    if (iz0 < lbound(rhs,3) .or. iz1 > ubound(rhs,3)) then
      write(*,*) 'scatter_rhs_from_global bounds error'
      write(*,*) 'iz0, iz1 = ', iz0, iz1
      write(*,*) 'rhs z bounds = ', lbound(rhs,3), ubound(rhs,3)
      error stop
    end if

    !$omp parallel do private(k) schedule(static)
    do k = 0, self%m+1
      self%rhs_dom(:,:,k) = rhs(:,:,iz0+k)
    end do
    !$omp end parallel do
  end subroutine scatter_rhs_from_global



  subroutine gather_phi_to_global(self, phi)
    class(PoissonDecomp), intent(inout) :: self
    real(real64), intent(inout) :: phi(0:self%nx+2,0:self%ny+2,0:self%nz+2)
    integer :: ierr
    integer :: sendcount

    sendcount = (self%kr - self%kl + 1) * (self%nx+3) * (self%ny+3)

    call MPI_Allgatherv( &
      self%phi_dom(0:self%nx+2, 0:self%ny+2, self%kl:self%kr), sendcount, MPI_DOUBLE_PRECISION, &
      phi, self%recvcounts_phi, self%displs_phi, MPI_DOUBLE_PRECISION, self%comm, ierr)
  end subroutine gather_phi_to_global


  subroutine gather_rhs_to_global(self, rhs)
    class(PoissonDecomp), intent(inout) :: self
    real(real64), intent(inout) :: rhs(0:self%nx+1,0:self%ny+1,0:self%nz+1)
    integer :: ierr
    integer :: sendcount

    sendcount = (self%kr - self%kl + 1) * (self%nx+2) * (self%ny+2)

    call MPI_Allgatherv( &
      self%rhs_dom(0:self%nx+1,0:self%ny+1,self%kl:self%kr), sendcount, MPI_DOUBLE_PRECISION, &
      rhs, self%recvcounts_rhs, self%displs_rhs, MPI_DOUBLE_PRECISION, self%comm, ierr)
  end subroutine gather_rhs_to_global


  subroutine destroy(self)
    class(PoissonDecomp), intent(inout) :: self
    if (allocated(self%phi_dom))        deallocate(self%phi_dom)
    if (allocated(self%rhs_dom))        deallocate(self%rhs_dom)
    if (allocated(self%bcnd_dom))       deallocate(self%bcnd_dom)
    if (allocated(self%recvcounts_phi)) deallocate(self%recvcounts_phi)
    if (allocated(self%displs_phi))     deallocate(self%displs_phi)
    if (allocated(self%recvcounts_rhs)) deallocate(self%recvcounts_rhs)
    if (allocated(self%displs_rhs))     deallocate(self%displs_rhs)
  end subroutine destroy

end module mod_poisson_decomp