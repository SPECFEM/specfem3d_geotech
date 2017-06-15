! routines requied for MPI proecsses
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module mpi_library
use mpi
contains

! start MPI processes
subroutine start_process(ismpi,ounit)
use global,only:myrank,nproc
implicit none
logical,intent(out) :: ismpi
integer,intent(in) :: ounit
integer :: errcode
ismpi=.true. ! parallel
call MPI_INIT(errcode)
if(errcode /= 0) call mpierror('ERROR: cannot initialize MPI!',errcode,ounit)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find processor ID (rank)!',errcode,ounit)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find number of processors!',errcode,ounit)
return
end subroutine start_process
!=======================================================

! close all MPI processes
subroutine close_process()
implicit none
integer :: errcode
call MPI_FINALIZE(errcode)
stop
return
end subroutine close_process
!=======================================================

! MPI error
subroutine mpierror(errtag,errcode,ounit)
implicit none
character(len=*) :: errtag
integer :: errcode,ounit
write(ounit,'(a,a,i4,a)')errtag,' (MPI ERROR code:',errcode,')!'
stop
end subroutine mpierror
!=======================================================

! syncronize all MPI processes
subroutine sync_process()
implicit none
integer :: errcode

call MPI_BARRIER(MPI_COMM_WORLD,errcode)

end subroutine sync_process
!=======================================================

! write error and stop
subroutine error_stop(errtag,ounit,myrank)
implicit none
character(*),intent(in) :: errtag
integer,intent(in) :: ounit,myrank
integer :: errcode,ierr
if(myrank==0)write(ounit,'(a)')errtag
call close_process
! stop all the MPI processes, and exit
write(ounit,'(a)')'aborting MPI...'
call MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
stop
end subroutine error_stop
!=======================================================

! get processor tag
function proc_tag() result(ptag)
use global,only:myrank,nproc
implicit none
character(len=20) :: format_str,ptag

write(format_str,*)ceiling(log10(real(nproc)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

write(ptag,fmt=format_str)'_proc',myrank

return
end function
!=======================================================

end module mpi_library
