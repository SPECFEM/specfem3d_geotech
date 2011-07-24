! routines requied to mimic similar MPI routines found in mpi_library.f90
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module serial_library
contains

subroutine start_process(myid,nproc,ounit)
implicit none
integer,intent(out) :: myid,nproc
integer,intent(in) :: ounit
integer :: errcode
myid=1; nproc=1
return
end subroutine start_process
!=======================================================

subroutine close_process()
implicit none
integer :: errcode
stop
return
end subroutine close_process
!=======================================================

subroutine sync_process()
implicit none
integer :: errcode
return
end subroutine sync_process
!=======================================================

! write error and stop
subroutine error_stop(errtag,ounit,myid)
implicit none
character(*),intent(in) :: errtag
integer,intent(in) :: ounit,myid
if(myid==1)write(ounit,'(a)')errtag
stop
end subroutine error_stop
!=======================================================

end module serial_library
