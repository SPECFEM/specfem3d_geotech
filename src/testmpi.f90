! this program is intended to test Fortran MPI compiler
program testmpi
use mpi
implicit none
integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL
integer :: ier

call MPI_INIT(ier)
call MPI_BARRIER(MPI_COMM_WORLD,ier)
call MPI_FINALIZE(ier)

end program testmpi
