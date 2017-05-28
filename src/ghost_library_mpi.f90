! this module contains the routines to process parallel communications across
! the ghost partitions
! TODO:
!   - better to avoid using allocatable component of derived type variable which
!     is fortran 95 standard? How much influence will be on the performance with
!     repeated allocate and deallocate
! REVISION:
!   HNG, APR 15,2011; HNG, Apr 09,2010
module ghost_library_mpi
use set_precision_mpi
use mpi_library
use global,only : myrank,nproc
private :: sort_array_coord,rank_buffers,swap_all_buffers
integer :: ngpart,maxngnode

! ghost partitions
type ghost_partition
  integer :: id,mindex,nnode
  integer,dimension(:),allocatable :: node,gdof
  logical,dimension(:),allocatable :: isnode
  ! whether the node on the interface is intact (.true.) or void (.false.)
end type ghost_partition
type(ghost_partition),dimension(:),allocatable :: gpart

contains

!-----------------------------------------------------------
subroutine prepare_ghost()
use global,only:ndim,nnode,nndof,ngllx,nglly,ngllz,g_num,g_coord,gdof,gfile,   &
part_path,stdout
use math_library, only : quick_sort
use math_library_mpi, only : maxvec

implicit none
integer :: istat

integer,dimension(8) :: ign,jgn,kgn ! ith, jth, and kth GLL indices of node
integer :: ig0,ig1,jg0,jg1,kg0,kg1
integer :: i_g,j_g,k_g
integer,dimension(6,4) :: node_face ! local node numbers in each face
integer,dimension(12,2) :: node_edge ! local node numbers in each edge

character(len=20) :: format_str
character(len=80) :: fname

integer :: mpartid,maxngpart ! partition ID
integer :: i_elmt,i_gpart,gpartid,melmt,ngelmt,igllp,inode,maxngelmt
integer :: etype,eid
integer :: ncount,new_ncount,ngllxy,maxngll_face
logical,dimension(nnode) :: switch_node
integer,dimension(:),allocatable :: itmp_array ! temporary integer array
double precision,dimension(:),allocatable :: xp,yp,zp
integer,dimension(3) :: ngll_vec ! (/ngllx,nglly,ngllz/) in ascending order
character(len=250) :: errtag

! find maximum ngll points in a face
! this is needed only for memory allocation
ngll_vec=quick_sort((/ngllx,nglly,ngllz/),3)
maxngll_face=ngll_vec(3)*ngll_vec(2)

! local node numbering in each face CUBIT/EXODUS convention
node_face(1,:)=(/1,2,6,5/) ! front
node_face(2,:)=(/2,3,7,6/) ! right
node_face(3,:)=(/4,3,7,8/) ! back
node_face(4,:)=(/1,4,8,5/) ! left
node_face(5,:)=(/1,2,3,4/) ! bottom
node_face(6,:)=(/5,6,7,8/) ! top

! local node numbering in each edge CUBIT/EXODUS convention
! bottom edges
node_edge(1,:)=(/1,2/);  node_edge(2,:)=(/2,3/)
node_edge(3,:)=(/3,4/);  node_edge(4,:)=(/4,1/)
! side edges
node_edge(5,:)=(/1,5/);  node_edge(6,:)=(/2,6/)
node_edge(7,:)=(/3,7/);  node_edge(8,:)=(/4,8/)
! top edges
node_edge(9,:)=(/5,6/);  node_edge(10,:)=(/6,7/)
node_edge(11,:)=(/7,8/); node_edge(12,:)=(/8,5/)

! GLL indices for nodes of a hexahedral mesh
! ith
ign(1)=1; ign(4)=1; ign(5)=1; ign(8)=1
ign(2)=ngllx; ign(3)=ngllx; ign(6)=ngllx; ign(7)=ngllx
! jth
jgn(1)=1; jgn(2)=1; jgn(5)=1; jgn(6)=1
jgn(4)=nglly; jgn(3)=nglly; jgn(8)=nglly; jgn(7)=nglly
! kth
kgn(1)=1; kgn(2)=1; kgn(4)=1; kgn(3)=1
kgn(5)=ngllz; kgn(6)=ngllz; kgn(8)=ngllz; kgn(7)=ngllz

ngllxy=ngllx*nglly

! open appropriate ghost file
write(format_str,*)ceiling(log10(real(nproc)+1))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
write(fname, fmt=format_str)trim(part_path)//trim(gfile)//'_proc',myrank
open(unit=11,file=trim(fname),status='old',action='read',iostat = istat)
if( istat /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  call error_stop(errtag,stdout,myrank)
endif

read(11,*) ! skip 1 line
read(11,*)mpartid ! master partition ID

if(mpartid/=myrank)then
  write(errtag,*)'ERROR: wrong gpart file partition ',mpartid,' !'
  call error_stop(errtag,stdout,myrank)
endif

read(11,*) ! skip 1 line
read(11,*)ngpart,maxngpart,maxngelmt

! allocate gpart
allocate(gpart(ngpart))
gpart(1:ngpart)%nnode=0
gpart(1:ngpart)%id=-1
gpart(1:ngpart)%mindex=-1
!allocate(cghost(ngpart)[*]) ! it will not work here, why?!:(

do i_gpart=1,ngpart ! ghost partitions loop
  read(11,*) ! skip 1 line
  read(11,*)gpartid,gpart(i_gpart)%mindex !gpart_mindex(i_gpart)
  gpartid=gpartid+1 ! index coarray starts from 1

  read(11,*) ! skip 1 line
  read(11,*)ngelmt
  allocate(itmp_array(ngelmt*maxngll_face)) ! I don't know why this doesn't work!
  itmp_array=-1
  switch_node=.false.
  ncount=0
  do i_elmt=1,ngelmt ! ghost elements loop
    read(11,*)melmt,etype,eid

    ! initialize
    ig0=-1; ig1=-1
    jg0=-1; jg1=-1
    kg0=-1; kg1=-1
    ! find range of GLL indices
    if(etype==1)then ! node
      ig0=ign(eid); ig1=ign(eid)
      jg0=jgn(eid); jg1=jgn(eid)
      kg0=kgn(eid); kg1=kgn(eid)
    elseif(etype==2)then ! edge
      ig0=minval(ign(node_edge(eid,:))); ig1=maxval(ign(node_edge(eid,:)))
      jg0=minval(jgn(node_edge(eid,:))); jg1=maxval(jgn(node_edge(eid,:)))
      kg0=minval(kgn(node_edge(eid,:))); kg1=maxval(kgn(node_edge(eid,:)))
    elseif(etype==4)then ! face
      ig0=minval(ign(node_face(eid,:))); ig1=maxval(ign(node_face(eid,:)))
      jg0=minval(jgn(node_face(eid,:))); jg1=maxval(jgn(node_face(eid,:)))
      kg0=minval(kgn(node_face(eid,:))); kg1=maxval(kgn(node_face(eid,:)))
    else
      write(errtag,*)'ERROR: wrong etype:',etype,' for ghost partition ',mpartid,'!'
      call error_stop(errtag,stdout,myrank)
    endif

    do k_g=kg0,kg1
      do j_g=jg0,jg1
        do i_g=ig0,ig1
          igllp=(k_g-1)*ngllxy+(j_g-1)*ngllx+i_g
          inode=g_num(igllp,melmt)
          if(.not.switch_node(inode))then
            ncount=ncount+1
            itmp_array(ncount)=inode
            switch_node(inode)=.true.
          endif
        enddo
      enddo
    enddo
  enddo ! do i_elmt

  gpart(i_gpart)%nnode=ncount
  gpart(i_gpart)%id=gpartid
  allocate(gpart(i_gpart)%node(ncount),gpart(i_gpart)%isnode(ncount))
  gpart(i_gpart)%isnode=.true. ! all nodes are active
  allocate(gpart(i_gpart)%gdof(ncount*nndof))

  gpart(i_gpart)%node=itmp_array(1:ncount)

  !order nodal array to match with ghost partitions
  !extract coordinates
  allocate(xp(ncount),yp(ncount),zp(ncount))

  xp=g_coord(1,itmp_array(1:ncount))
  yp=g_coord(2,itmp_array(1:ncount))
  zp=g_coord(3,itmp_array(1:ncount))
  deallocate(itmp_array)

  !sort nodes
  call sort_array_coord(ndim,ncount,xp,yp,zp,gpart(i_gpart)%node,new_ncount)
  !call sort_array_coord(ndim,ncount,xp,yp,zp,gpart_node(i_gpart,1:ncount),    &
  !new_ncount)
  deallocate(xp,yp,zp)
  if(ncount/=new_ncount)then
    write(errtag,*)'ERROR: number of ghost nodes mismatched after sorting!'
    call error_stop(errtag,stdout,myrank)
  endif

  ! find ghost gdof
  gpart(i_gpart)%gdof=reshape(gdof(:,gpart(i_gpart)%node),(/ncount*nndof/))

enddo ! do i_gpart
close(11)
!stop
maxngnode=maxvec(gpart(1:ngpart)%nnode)

return
end subroutine prepare_ghost
!=======================================================

! modify ghost gdof based on the modified gdof
subroutine modify_ghost(isnode)
use global,only:nnode,nndof,gdof

implicit none
logical,intent(in) :: isnode(nnode)

integer :: i,i_gpart,ncount

! modify ghost gdof based on the modified gdof
do i_gpart=1,ngpart
  ncount=gpart(i_gpart)%nnode
  gpart(i_gpart)%gdof=reshape(gdof(:,gpart(i_gpart)%node),(/ncount*nndof/))
  do i=1,ncount
    gpart(i_gpart)%isnode(i)=isnode(gpart(i_gpart)%node(i))
  enddo
enddo ! do i_gpart
return
end subroutine modify_ghost
!=======================================================

! this subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts(neq,array,array_g)
use global,only:nndof
!use math_library, only : maxscal_par
use mpi
implicit none
integer,intent(in) :: neq
real(kind=kreal),dimension(0:neq),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g
real(kind=kreal),dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
real(kind=kreal),parameter :: zero=0.0_kreal

integer :: ierr,i_gpart,ncount
array_g=array
send_array=zero; recv_array=zero
!do i_gpart=1,ngpart
!  ncount=gpart(i_gpart)%nnode*nndof
  ! array to send
!  send_array(1:ncount,i_gpart)=array(gpart(i_gpart)%gdof)
!enddo
!call sync_process()
do i_gpart=1,ngpart
  ncount=gpart(i_gpart)%nnode*nndof
  ! array to send
  send_array(1:ncount,i_gpart)=array(gpart(i_gpart)%gdof)
  ! send
  call MPI_ISSEND(send_array(1,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%id-1,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(1,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%id-1,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ncount=gpart(i_gpart)%nnode*nndof
  array_g(gpart(i_gpart)%gdof)=array_g(gpart(i_gpart)%gdof)+ &
  recv_array(1:ncount,i_gpart)
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
array_g(0)=zero
return
end subroutine assemble_ghosts
!=======================================================

! this subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts_nodal(nndof,array,array_g)
use mpi
use global,only:nnode
implicit none
integer,intent(in) :: nndof
! number of active ghost partitions for a node
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(nndof,nnode),intent(out) :: array_g

!logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
real(kind=kreal),dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
!integer,dimension(nnode) :: ngpart_node ! number of ghost partition for a node
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(:,gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(1,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%id-1,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(1,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%id-1,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  array_g(:,gpart(i_gpart)%node)=array_g(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/nndof,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_nodal
!=======================================================

! this subroutine counts the active ghost partitions for each node on the
! interfaces.
! logical flag representing whether the nodes in the interfaces are intact or
! void has to be communicated across the processors
subroutine count_active_nghosts(ngpart_node)
use mpi
use global,only:nnode
implicit none
! number of active ghost partitions for a node
integer,dimension(nnode),intent(out) :: ngpart_node
! only the interfacial nodes can be saved for the storage (TODO)

logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
integer :: i,ierr,i_gpart,ignode,ngnode

ngpart_node=0
lsend_array=.true.; lrecv_array=.true.
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode

  ! store array-to-send
  lsend_array(1:ngnode,i_gpart)=gpart(i_gpart)%isnode(1:ngnode)

  ! send
  call MPI_ISSEND(lsend_array(1,i_gpart),ngnode,MPI_LOGICAL,gpart(i_gpart)%id-1,&
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(lrecv_array(1,i_gpart),ngnode,MPI_LOGICAL,gpart(i_gpart)%id-1,&
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! count active partitons along the interfaces
do i_gpart=1,ngpart
  do i=1,gpart(i_gpart)%nnode
    ignode=gpart(i_gpart)%node(i)
    if(lrecv_array(i,i_gpart))ngpart_node(ignode)=ngpart_node(ignode)+1
  enddo
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()

return
end subroutine count_active_nghosts
!=======================================================

! this subroutine distributes the excavation loads discarded by a processors due
! to the special geoemtry partition. it will not distribute if the load is used
! within the partition
subroutine distribute2ghosts(gdof,nndof,neq,ngpart_node,array,array_g)
use mpi
use global,only:nnode
implicit none
integer,intent(in) :: nndof,neq
integer,dimension(nndof,nnode),intent(in) :: gdof ! global degree of freedom
! number of active ghost partitions for a node
integer,dimension(nnode),intent(in) :: ngpart_node
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g

real(kind=kreal),dimension(nndof,nnode) :: tarray
!logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
real(kind=kreal),dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
real(kind=kreal),dimension(nndof,maxngnode) :: garray
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
!integer,dimension(nnode) :: ngpart_node ! number of ghost partition for a node
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: i,j,ierr,i_gpart,igdof,ignode,ncount,ngnode

array_g=zero
tarray=array; send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  ! store array-to-send in a garray
  garray=zero
  garray(:,1:ngnode)=array(:,gpart(i_gpart)%node)
  ! make appropriate correction so that only the discarded loads will equally be
  ! distributed
  do j=1,ngnode
    ignode=gpart(i_gpart)%node(j)
    if(ngpart_node(ignode)==0)then
      ! this node is dead
      garray(:,j)=zero !stop 'strange!'
      cycle
    endif
    do i=1,nndof
      if(gdof(i,ignode)/=0)then
        ! do not distribute if the load is already used in the partition
        garray(i,j)=zero !garray(i,j)/real(ngpart_node(ignode)+1,kreal) !
        !tarray(i,ignode)=garray(i,j)
      else
        ! distribute equally to ghost partitions if the load is discarded by the
        ! partition.
        !if(ngpart_node(ignode)==0)stop 'strange!'
        garray(i,j)=garray(i,j)/real(ngpart_node(ignode),kreal)
        tarray(i,ignode)=zero
      endif
    enddo
  enddo
  send_array(1:ncount,i_gpart)=reshape(garray(:,1:ngnode),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(1,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%id-1,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(1,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%id-1,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  tarray(:,gpart(i_gpart)%node)=tarray(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/nndof,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()

! store nodal values to gdof locations
do j=1,nnode
  do i=1,nndof
    igdof=gdof(i,j)
    array_g(igdof)=tarray(i,j)
  enddo
enddo
array_g(0)=zero
return
end subroutine distribute2ghosts
!===========================================

! deallocate ghost variables
subroutine free_ghost()
implicit none
integer :: i
do i=1,ngpart
  deallocate(gpart(i)%node,gpart(i)%gdof)
enddo
deallocate(gpart)
return
end subroutine free_ghost
!===========================================

! routines below are imported and modified from SPECFEM3D

! subroutines to sort MPI buffers to assemble between chunks
subroutine sort_array_coord(ndim,npoint,x,y,z,ibool,nglob)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points
!
! returns: sorted indexing array (ibool),  reordering array (iglob) &
! number of global points (nglob)

 use math_constants, only : zerotol
 implicit none

 integer,intent(in) :: ndim,npoint
 double precision,dimension(npoint),intent(in) :: x,y,z
 integer,dimension(npoint),intent(inout) :: ibool
 integer,intent(out) :: nglob

 integer,dimension(npoint) :: iglob
 integer,dimension(npoint) :: iloc,ind,ninseg
 logical,dimension(npoint) :: ifseg

 integer,dimension(npoint) :: iwork(npoint)
 double precision,dimension(npoint) :: work(npoint)

 integer :: i,j
 integer :: nseg,ioff,iseg,ig
 double precision :: xtol

! establish initial pointers
 do i=1,npoint
   iloc(i)=i
 enddo

! define a tolerance, normalized radius is 1., so let's use a small value
 xtol = zerotol !SMALLVAL_TOL

 ifseg(:)=.false.

 nseg=1
 ifseg(1)=.true.
 ninseg(1)=npoint

 do j=1,NDIM

! sort within each segment
  ioff=1
  do iseg=1,nseg
    if(j == 1) then

      call rank_buffers(x(ioff),ind,ninseg(iseg))

    else if(j == 2) then

      call rank_buffers(y(ioff),ind,ninseg(iseg))

    else

      call rank_buffers(z(ioff),ind,ninseg(iseg))

    endif

    call swap_all_buffers(ibool(ioff),iloc(ioff), &
            x(ioff),y(ioff),z(ioff),iwork,work,ind,ninseg(iseg))

    ioff=ioff+ninseg(iseg)
  enddo

! check for jumps in current coordinate
  if(j == 1) then
    do i=2,npoint
      if(dabs(x(i)-x(i-1)) > xtol) ifseg(i)=.true.
    enddo
  else if(j == 2) then
    do i=2,npoint
      if(dabs(y(i)-y(i-1)) > xtol) ifseg(i)=.true.
    enddo
  else
    do i=2,npoint
      if(dabs(z(i)-z(i-1)) > xtol) ifseg(i)=.true.
    enddo
  endif

! count up number of different segments
  nseg=0
  do i=1,npoint
    if(ifseg(i)) then
      nseg=nseg+1
      ninseg(nseg)=1
    else
      ninseg(nseg)=ninseg(nseg)+1
    endif
  enddo
  enddo

! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npoint
    if(ifseg(i)) ig=ig+1
    iglob(iloc(i))=ig
  enddo

  nglob=ig

  end subroutine sort_array_coord
!===========================================

! -------------------- library for sorting routine ------------------

! sorting routines put here in same file to allow for inlining

  subroutine rank_buffers(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer :: n
  double precision :: A(n)
  integer :: IND(n)

  integer :: i,j,l,ir,indx
  double precision :: q

  do j=1,n
    IND(j)=j
  enddo

  if(n == 1) return

  L=n/2+1
  ir=n
  100 CONTINUE
   IF(l>1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   ENDIF
   i=l
   j=l+l
  200    CONTINUE
   IF(J <= IR) THEN
      IF(J < IR) THEN
         IF(A(IND(j)) < A(IND(j+1))) j=j+1
      ENDIF
      IF (q < A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      ENDIF
   goto 200
   ENDIF
   IND(I)=INDX
  goto 100
  end subroutine rank_buffers
!===========================================

  subroutine swap_all_buffers(IA,IB,A,B,C,IW,W,ind,n)
!
! swap arrays IA, IB, A, B and C according to addressing in array IND
!
  implicit none

  integer :: n

  integer :: IND(n)
  integer :: IA(n),IB(n),IW(n)
  double precision :: A(n),B(n),C(n),W(n)

  integer :: i

  do i=1,n
    W(i)=A(i)
    IW(i)=IA(i)
  enddo

  do i=1,n
    A(i)=W(ind(i))
    IA(i)=IW(ind(i))
  enddo

  do i=1,n
    W(i)=B(i)
    IW(i)=IB(i)
  enddo

  do i=1,n
    B(i)=W(ind(i))
    IB(i)=IW(ind(i))
  enddo

  do i=1,n
    W(i)=C(i)
  enddo

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all_buffers
!===========================================

end module ghost_library_mpi


