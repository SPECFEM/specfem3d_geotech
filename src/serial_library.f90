! routines requied to mimic similar MPI routines found in mpi_library.f90 and
! others
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module serial_library
use set_precision
contains

subroutine start_process(ismpi,ounit)
use global,only:myrank,nproc
implicit none
logical,intent(out) :: ismpi
integer,intent(in) :: ounit
ismpi=.false. ! serial
myrank=0; nproc=1
return
end subroutine start_process
!=======================================================

subroutine close_process()
implicit none
stop
return
end subroutine close_process
!=======================================================

subroutine sync_process()
implicit none
return
end subroutine sync_process
!=======================================================

! write error and stop
subroutine error_stop(errtag,ounit,myrank)
implicit none
character(*),intent(in) :: errtag
integer,intent(in) :: ounit,myrank
if(myrank==1)write(ounit,'(a)')errtag
stop
end subroutine error_stop
!=======================================================

! get processor tag
function proc_tag() result(ptag)
character(len=20) :: ptag

ptag=''

return
end function
!=======================================================

subroutine prepare_ghost()
implicit none
return
end subroutine prepare_ghost
!=======================================================

subroutine modify_ghost(isnode)
use global,only:nnode,nndof
implicit none
logical,intent(in) :: isnode(nnode)
return
end subroutine modify_ghost
!=======================================================

subroutine assemble_ghosts(neq,array,array_g)
implicit none
integer,intent(in) :: neq
real(kind=kreal),dimension(0:neq),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g
array_g=array
return
end subroutine assemble_ghosts
!=======================================================

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
!===========================================

! this subroutine counts the active ghost partitions for each node on the
! interfaces.
! logical flag representing whether the nodes in the interfaces are intact or
! void has to be communicated across the processors
subroutine count_active_nghosts(ngpart_node)
use global,only:nnode
implicit none
! number of active ghost partitions for a node
integer,dimension(nnode),intent(out) :: ngpart_node
! only the interfacial nodes can be saved for the storage (TODO)
ngpart_node=0
return
end subroutine count_active_nghosts
!===========================================

! this subroutine distributes the excavation loads discarded by a processors due
! to the special geoemtry partition. it will not distribute if the load is used
! within the partition
subroutine distribute2ghosts(gdof,nndof,neq,ngpart_node, &
array,array_g)
use global,only:nnode
implicit none
integer,intent(in) :: nndof,neq
integer,dimension(nndof,nnode),intent(in) :: gdof ! global degree of freedom
! number of active ghost partitions for a node
integer,dimension(nnode),intent(in) :: ngpart_node
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
!===========================================

! deallocate ghost variables
subroutine free_ghost()
implicit none
return
end subroutine free_ghost
!===========================================

end module serial_library

