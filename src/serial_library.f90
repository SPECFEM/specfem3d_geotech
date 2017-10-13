! routines requied to mimic similar MPI routines found in mpi_library.f90 and
! others
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module serial_library
use set_precision

contains

!-------------------------------------------------------------------------------
subroutine start_process()
use global,only:ismpi,myrank,nproc
implicit none
ismpi=.false. ! serial
myrank=0; nproc=1
return
end subroutine start_process
!===============================================================================

subroutine close_process()
implicit none
stop
return
end subroutine close_process
!===============================================================================

subroutine sync_process()
implicit none
return
end subroutine sync_process
!===============================================================================

! write error and stop
subroutine control_error(errtag)
use global,only:myrank,stdout
implicit none
character(*),intent(in) :: errtag
if(myrank==0)write(stdout,'(a)')errtag
stop
end subroutine control_error
!===============================================================================

! get processor tag
function proc_tag() result(ptag)
character(len=20) :: ptag

ptag=''

return
end function
!===============================================================================

subroutine prepare_ghost()
implicit none
return
end subroutine prepare_ghost
!===============================================================================

subroutine prepare_ghost_gdof()
implicit none
return
end subroutine prepare_ghost_gdof
!===============================================================================

subroutine modify_ghost()
implicit none
return
end subroutine modify_ghost
!===============================================================================

subroutine assemble_ghosts(neq,array,array_g)
implicit none
integer,intent(in) :: neq
real(kind=kreal),dimension(0:neq),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g
array_g=array
return
end subroutine assemble_ghosts
!===============================================================================

! this subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts_nodal(nndof,array,array_g)
use global,only:nnode
implicit none
integer,intent(in) :: nndof
! number of active ghost partitions for a node
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(nndof,nnode),intent(out) :: array_g
array_g=array
return
end subroutine assemble_ghosts_nodal
!===============================================================================

! this subroutine counts the active ghost partitions for each node on the
! interfaces.
! logical flag representing whether the nodes in the interfaces are intact or
! void has to be communicated across the processors
subroutine count_active_nghosts()
implicit none
return
end subroutine count_active_nghosts
!===============================================================================

! this subroutine distributes the excavation loads discarded by a processors due
! to the special geoemtry partition. it will not distribute if the load is used
! within the partition
subroutine distribute2ghosts(gdof,nndof,neq,array,array_g)
use global,only:nnode
implicit none
integer,intent(in) :: nndof,neq
integer,dimension(nndof,nnode),intent(in) :: gdof ! global degree of freedom
! number of active ghost partitions for a node
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g

integer :: i,j,igdof
real(kind=kreal),parameter :: zero=0.0_kreal
! store nodal values to gdof locations
do j=1,nnode
  do i=1,nndof
    igdof=gdof(i,j)
    array_g(igdof)=array(i,j)
  enddo
enddo
array_g(0)=zero
return
end subroutine distribute2ghosts
!===============================================================================

! deallocate ghost variables
subroutine free_ghost()
implicit none
return
end subroutine free_ghost
!===============================================================================

end module serial_library
!===============================================================================
