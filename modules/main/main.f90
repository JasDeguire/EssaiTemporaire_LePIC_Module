program main
  use mpi
  use mod_simulation, only: Simulation
  implicit none

  type(Simulation) :: sim
  integer :: ierr

  call MPI_Init(ierr)

  call sim%init(MPI_COMM_WORLD)
  call sim%run(sim%state%cfg%nsteps)
  call sim%finalize()

  call MPI_Finalize(ierr)
end program main
