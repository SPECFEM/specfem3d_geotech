! this module contains the routines to partition mesh
! most of routines in this library are copied and modified from the original
! SPECFEM3D package (Komatitsch and Tromp 1999, Peter et al. 2011)
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module partmesh_library

implicit none

! Useful kind types
integer,parameter :: kreal=selected_real_kind(15)
integer,parameter :: short = SELECTED_INT_KIND(4), long = SELECTED_INT_KIND(18)

! Number of nodes per elements.
integer,parameter  :: ESIZE = 8

! Number of faces per element.
integer,parameter  :: nfaces = 6

! very large and very small values
real(kind=kreal),parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9 !double precision

! acoustic-elastic load balancing:
! assumes that elastic at least ~6 times more expensive than acoustic
integer,parameter :: ACOUSTIC_LOAD = 1
integer,parameter :: ELASTIC_LOAD = 4

character(len=1),parameter :: CR=achar(13) ! carriage return to overwrite
integer :: nelmt_interface
integer,dimension(:),allocatable :: elmt_interface,part_interface
integer,dimension(:),allocatable :: node_npart
integer,dimension(:,:),allocatable :: connect_interface
type nodes
  integer,dimension(:),allocatable :: part
end type nodes
type(nodes),dimension(:),allocatable :: node

!  include './constants_decompose_mesh_SCOTCH.h'

contains

!-----------------------------------------------
! Creating dual graph (adjacency is defined by 'ncommonnodes' between two elements).
!-----------------------------------------------
subroutine mesh2dual_ncommonnodes(nelmnts,nnodes,nsize, sup_neighbour,elmnts,xadj,adjncy,nnodes_elmnts, &
nodes_elmnts,max_neighbour,ncommonnodes)

integer,intent(in) :: nelmnts !integer(long)
integer,intent(in) :: nnodes
integer,intent(in) :: nsize !integer(long)
integer,intent(in) :: sup_neighbour !integer(long)
integer,dimension(0:esize*nelmnts-1),intent(in) :: elmnts    

integer,dimension(0:nelmnts) :: xadj
integer,dimension(0:sup_neighbour*nelmnts-1) :: adjncy
integer,dimension(0:nnodes-1) :: nnodes_elmnts
integer,dimension(0:nsize*nnodes-1) :: nodes_elmnts
integer,intent(out) :: max_neighbour
integer,intent(in) :: ncommonnodes

! local parameters
integer :: i,j,k,l,m,nb_edges
logical :: is_neighbour
integer :: num_node,n
integer :: elem_base,elem_target
integer :: connectivity


! initializes
xadj(:) = 0
adjncy(:) = 0
nnodes_elmnts(:) = 0
nodes_elmnts(:) = 0
nb_edges = 0

! list of elements per node
do i = 0, esize*nelmnts-1
    nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
    nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
end do

! checking which elements are neighbours ('ncommonnodes' criteria)
do j = 0, nnodes-1
    do k = 0, nnodes_elmnts(j)-1
      do l = k+1, nnodes_elmnts(j)-1

          connectivity = 0
          elem_base = nodes_elmnts(k+j*nsize)
          elem_target = nodes_elmnts(l+j*nsize)
          do n = 1, esize
            num_node = elmnts(esize*elem_base+n-1)
            do m = 0, nnodes_elmnts(num_node)-1
                if ( nodes_elmnts(m+num_node*nsize) == elem_target ) then
                  connectivity = connectivity + 1
                end if
            end do
          end do

          if ( connectivity >=  ncommonnodes) then

            is_neighbour = .false.

            do m = 0, xadj(nodes_elmnts(k+j*nsize))
                if ( .not.is_neighbour ) then
                  if ( adjncy(nodes_elmnts(k+j*nsize)*sup_neighbour+m) == nodes_elmnts(l+j*nsize) ) then
                      is_neighbour = .true.

                  end if
                end if
            end do
            if ( .not.is_neighbour ) then
                adjncy(nodes_elmnts(k+j*nsize)*sup_neighbour+xadj(nodes_elmnts(k+j*nsize))) = nodes_elmnts(l+j*nsize)
                
                xadj(nodes_elmnts(k+j*nsize)) = xadj(nodes_elmnts(k+j*nsize)) + 1
                if (xadj(nodes_elmnts(k+j*nsize))>sup_neighbour) stop 'ERROR : too much neighbours per element, modify the mesh.'

                adjncy(nodes_elmnts(l+j*nsize)*sup_neighbour+xadj(nodes_elmnts(l+j*nsize))) = nodes_elmnts(k+j*nsize)

                xadj(nodes_elmnts(l+j*nsize)) = xadj(nodes_elmnts(l+j*nsize)) + 1
                if (xadj(nodes_elmnts(l+j*nsize))>sup_neighbour) stop 'ERROR : too much neighbours per element, modify the mesh.'
            end if
          end if
      end do
    end do
end do

max_neighbour = maxval(xadj)

! making adjacency arrays compact (to be used for partitioning)
do i = 0, nelmnts-1
    k = xadj(i)
    xadj(i) = nb_edges
    do j = 0, k-1
      adjncy(nb_edges) = adjncy(i*sup_neighbour+j)
      nb_edges = nb_edges + 1
    end do
end do

xadj(nelmnts) = nb_edges
end subroutine mesh2dual_ncommonnodes
!=======================================================


!--------------------------------------------------
! construct local numbering for the elements in each partition
!--------------------------------------------------
subroutine Construct_glob2loc_elmnts(nelmnts, part, glob2loc_elmnts,npart)

! include './constants_decompose_mesh_SCOTCH.h'

integer,intent(in) :: nelmnts !integer(long)
integer,dimension(0:nelmnts-1),intent(in) :: part
integer,dimension(:),pointer :: glob2loc_elmnts

integer :: num_glob,num_part,npart
integer,dimension(0:npart-1) :: num_loc

! allocates local numbering array
allocate(glob2loc_elmnts(0:nelmnts-1))

! initializes number of local points per partition
do num_part = 0, npart-1
    num_loc(num_part) = 0
end do

! local numbering
do num_glob = 0, nelmnts-1
    ! gets partition
    num_part = part(num_glob)
    ! increments local numbering of elements (starting with 0,1,2,...)
    glob2loc_elmnts(num_glob) = num_loc(num_part)
    num_loc(num_part) = num_loc(num_part) + 1
end do

end subroutine Construct_glob2loc_elmnts
!=======================================================


!--------------------------------------------------
! construct local numbering for the nodes in each partition
!--------------------------------------------------
subroutine Construct_glob2loc_nodes(nelmnts, nnodes, nsize, nnodes_elmnts, nodes_elmnts, part, &
      glob2loc_nodes_npart, glob2loc_nodes_parts, glob2loc_nodes,npart)

! include './constants_decompose_mesh_SCOTCH.h'

integer,intent(in) :: nelmnts,nsize !integer(long)
integer,intent(in) :: nnodes
integer,dimension(0:nelmnts-1),intent(in) :: part
integer,dimension(0:nnodes-1),intent(in) :: nnodes_elmnts
integer,dimension(0:nsize*nnodes-1),intent(in) :: nodes_elmnts
integer,dimension(:),pointer :: glob2loc_nodes_npart
integer,dimension(:),pointer :: glob2loc_nodes_parts
integer,dimension(:),pointer :: glob2loc_nodes

integer :: num_node
integer :: el
integer :: num_part
integer :: size_glob2loc_nodes,npart
integer,dimension(0:npart-1) :: parts_node
integer,dimension(0:npart-1) :: num_parts

allocate(glob2loc_nodes_npart(0:nnodes))

size_glob2loc_nodes = 0
parts_node(:) = 0

do num_node = 0, nnodes-1
    glob2loc_nodes_npart(num_node) = size_glob2loc_nodes
    do el = 0, nnodes_elmnts(num_node)-1
      parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1

    end do

    do num_part = 0, npart-1
      if ( parts_node(num_part) == 1 ) then
          size_glob2loc_nodes = size_glob2loc_nodes + 1
          parts_node(num_part) = 0

      end if
    end do

end do

glob2loc_nodes_npart(nnodes) = size_glob2loc_nodes

allocate(glob2loc_nodes_parts(0:glob2loc_nodes_npart(nnodes)-1))
allocate(glob2loc_nodes(0:glob2loc_nodes_npart(nnodes)-1))

glob2loc_nodes(0) = 0

parts_node(:) = 0
num_parts(:) = 0
size_glob2loc_nodes = 0


do num_node = 0, nnodes-1
    do el = 0, nnodes_elmnts(num_node)-1
      parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1

    end do
    do num_part = 0, npart-1

      if ( parts_node(num_part) == 1 ) then
          glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
          glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
          size_glob2loc_nodes = size_glob2loc_nodes + 1
          num_parts(num_part) = num_parts(num_part) + 1
          parts_node(num_part) = 0
      end if

    end do
end do
end subroutine Construct_glob2loc_nodes
!=======================================================

subroutine write_ssbc(out_path,inp_path,bcfile,nelmt,part,npart,glob2loc_elmt)
! this subroutine writes the partitioned side-set boundary conditions
integer,intent(in) :: nelmt !integer(long)

integer,dimension(nelmt),intent(in) :: part ! numbering starts from 1 only in this routine
integer,intent(in) :: npart
integer,dimension(:),intent(in):: glob2loc_elmt ! numbering starts from 1 only in this routine
character(len=256),intent(in) :: inp_path,out_path
character(len=80),intent(in) :: bcfile
character(len=20) :: format_str,format_str1
character(len=80) :: out_fname

integer,dimension(npart) :: mpart_icount,mpart_nelmt
integer :: bc_nelmt ! number of BC elements 
integer,dimension(:,:),allocatable :: bc_elmt ! fist row = element ID, second row = face ID
integer :: i_elmt,i_part,ipart,istat

type master_partition
integer,dimension(:),allocatable :: iloc ! index location of element in the partition 
end type master_partition
type(master_partition),dimension(npart) :: mpart

! open BC file
open(unit=16,file=trim(inp_path)//trim(bcfile),status='old',action='read',iostat=istat)
if(istat/=0) then
print*,'error file open:',trim(bcfile)
stop
endif
read(16,*)bc_nelmt
if(bc_nelmt>0)then
allocate(bc_elmt(2,bc_nelmt))
read(16,*,iostat=istat)bc_elmt
if(istat/=0)then
  write(*,*)'ERROR: cannot read BC information!'
  stop
endif
endif
close(16)
! find number of BC elements in each partitions
mpart_nelmt=0
do i_elmt=1,bc_nelmt
mpart_nelmt(part(bc_elmt(1,i_elmt)))=mpart_nelmt(part(bc_elmt(1,i_elmt)))+1
enddo
! allocate derived type variables
do i_part=1,npart
allocate(mpart(i_part)%iloc(mpart_nelmt(i_part)))
enddo
mpart_icount=0
! partition BC
do i_elmt=1,bc_nelmt
ipart=part(bc_elmt(1,i_elmt)) ! partion
mpart_icount(ipart)=mpart_icount(ipart)+1  
mpart(ipart)%iloc(mpart_icount(ipart))=i_elmt
enddo

bc_elmt(1,:)=glob2loc_elmt(bc_elmt(1,:)) ! local element numbering in the partition

! format string for file name
write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

! format string for element ID and face ID
write(format_str1,*)ceiling(log10(real(maxval(bc_elmt(1,:)))+1.))
format_str1='(i'//trim(adjustl(format_str1))//'i2)' ! i2 is sufficient for face id or node id 

! write BCs in each partition
do i_part=1,npart
! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(bcfile)//'_proc',i_part-1
open(unit=16,file=trim(out_fname),status='replace',action='write',iostat=istat)
if(istat/=0) then
  print*,'ERROR: cannot open file "'//trim(out_fname)//'!'
  stop
endif
write(16,*)mpart_nelmt(i_part)
write(16,format_str1)bc_elmt(:,mpart(i_part)%iloc)
close(16)
enddo ! i_part
deallocate(bc_elmt)
! deallocate derived type variables
do i_part=1,npart
deallocate(mpart(i_part)%iloc)
enddo

end subroutine write_ssbc
!=======================================================

subroutine write_traction(out_path,inp_path,trfile,nelmt,part,npart,glob2loc_elmt)
! this subroutine writes the partitioned traction boundary conditions
integer,intent(in) :: nelmt !integer(long)

integer,dimension(nelmt),intent(in) :: part ! numbering starts from 1 only in this routine
integer,intent(in) :: npart
integer,dimension(:),intent(in):: glob2loc_elmt ! numbering starts from 1 only in this routine
character(len=256),intent(in) :: inp_path,out_path
character(len=80),intent(in) :: trfile
character(len=20) :: format_str,format_str1
character(len=80) :: out_fname


integer,dimension(npart) :: mpart_icount,mpart_nelmt
integer :: tr_nelmt ! number of traction elements 
integer,dimension(:,:),allocatable :: tr_elmt ! fist row = element ID, second row = entity ID
integer :: ios,i_elmt,i_part,itrac,ipart,istat,ntrac,tractype

integer :: iaxis
real(kind=kreal) :: q0(3),q1(3),x1,x2
logical :: ispart(npart)

type master_partition
integer,dimension(:),allocatable :: iloc ! index location of element in the partition 
end type master_partition
type(master_partition),dimension(npart) :: mpart

! open traction file
open(unit=16,file=trim(inp_path)//trim(trfile),status='old',action='read',iostat=istat)
if(istat/=0) then
print*,'error file open:',trim(trfile)
stop
endif
!read(16,*)ntrac
itrac=0
do ! itrac=1,ntrac
read(16,*,iostat=ios)tractype  
if(ios/=0)exit
itrac=itrac+1
if(tractype==0)then ! point loading
  read(16,*)q0
elseif(tractype==1)then ! uniform loading
  read(16,*)q0
elseif(tractype==2)then ! linear loading
  read(16,*)iaxis,x1,x2,q0,q1
else
  write(*,*)'ERROR: unsupported traction type:',tractype
  stop
endif
read(16,*)tr_nelmt
if(tr_nelmt>0)then
  allocate(tr_elmt(2,tr_nelmt))
  read(16,*,iostat=istat)tr_elmt
  if(istat/=0)then
    write(*,*)'ERROR: cannot read traction information!'
    stop
  endif
endif

! find number of traction elements in each partitions
mpart_nelmt=0
do i_elmt=1,tr_nelmt
  mpart_nelmt(part(tr_elmt(1,i_elmt)))=mpart_nelmt(part(tr_elmt(1,i_elmt)))+1
enddo
! allocate derived type variables
ispart=.false.
do i_part=1,npart
  if(mpart_nelmt(i_part)>0)then
    ispart(i_part)=.true.
  endif
  allocate(mpart(i_part)%iloc(mpart_nelmt(i_part)))
enddo
mpart_icount=0
! partition BC
do i_elmt=1,tr_nelmt
  ipart=part(tr_elmt(1,i_elmt)) ! partion
  mpart_icount(ipart)=mpart_icount(ipart)+1  
  mpart(ipart)%iloc(mpart_icount(ipart))=i_elmt
enddo

tr_elmt(1,:)=glob2loc_elmt(tr_elmt(1,:)) ! local element numbering in the partition

! format string for file name
write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

! format string for element ID and face ID
write(format_str1,*)ceiling(log10(real(maxval(tr_elmt(1,:)))+1.))
format_str1='(i'//trim(adjustl(format_str1))//'i2)' ! i2 is sufficient for face id or node id 
  
! write BCs in each partition
do i_part=1,npart    
  ! open output file
  write(out_fname, fmt=format_str)trim(out_path)//trim(trfile)//'_proc',i_part-1
  if(itrac==1)then
    open(unit=17,file=trim(out_fname),status='replace',action='write',iostat=istat)
  else
    open(unit=17,file=trim(out_fname),status='old',action='write',position='append',iostat=istat)
  endif
  if(istat/=0) then
    print*,'ERROR: cannot open file "'//trim(out_fname)//'!'
    stop
  endif
  if(ispart(i_part))then
    write(17,*)tractype
    if(tractype==0)then
      write(17,'(3(f12.6,1x))')q0
    elseif(tractype==1)then
      write(17,'(3(f12.6,1x))')q0
    elseif(tractype==2)then
      write(17,'(i1,1x,2(f12.6,1x),6(f12.6,1x))')iaxis,x1,x2,q0,q1
    endif
    write(17,*)mpart_nelmt(i_part)
    write(17,format_str1)tr_elmt(:,mpart(i_part)%iloc)
  endif
  close(17)
enddo ! i_part
deallocate(tr_elmt)
! deallocate derived type variables
do i_part=1,npart
  deallocate(mpart(i_part)%iloc)
enddo  
enddo ! itrac
close(16)

write(*,'(a,i3,a)',advance='no')'ntrac=',itrac,'...'

end subroutine write_traction
!=======================================================

subroutine find_interface(nelmt,nnode,part,connect,npart)
! this subroutine find the elements on the interfaces. This routine only works for 8-noded
! hexahedral elements
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
!TODO
! - some optimization is possible
! - after finding interfacial elements we can only store/process the information 
!   only for those elements
implicit none
integer,intent(in) :: nelmt,nnode
integer,dimension(nelmt),intent(in) :: part ! numbering starts from 1 only in this routine
integer,dimension(esize,nelmt),intent(in) :: connect    
integer,intent(in) :: npart

integer :: i,icount,istat,j
logical,dimension(:,:),allocatable :: node_part

allocate(node_part(nnode,npart))
allocate(node_npart(nnode))

node_part=.false.  
do i=1,nelmt
node_part(connect(:,i),part(i))=.true.
enddo

! determine the number of partitions for all nodes
node_npart=0
allocate(node(nnode))
do i=1,nnode 
node_npart(i)=count(node_part(i,:))
allocate(node(i)%part(node_npart(i)))
icount=0
do j=1,npart
  if(node_part(i,j))then
    icount=icount+1
    node(i)%part(icount)=j
  endif
enddo
enddo 
!print*,'hello:',maxval(node_part)
!stop
deallocate(node_part)

! count the number of elements in the interface  
nelmt_interface=0
do i=1,nelmt
if(maxval(node_npart(connect(:,i)))>1)nelmt_interface=nelmt_interface+1
enddo
!max_share=maxval(node_npart)-1
print*,'total elements on the interface: ',nelmt_interface
!stop

! number of neighboring element for the interface elements
allocate(elmt_interface(nelmt_interface))

icount=0
do i=1,nelmt
if(maxval(node_npart(connect(:,i)))>1)then
  icount=icount+1
!    nelmt_share(icount)=maxval(node_npart(connect(:,i)))-1
  elmt_interface(icount)=i
endif
enddo
allocate(connect_interface(esize,nelmt_interface))
connect_interface=connect(:,elmt_interface)
allocate(part_interface(nelmt_interface))
part_interface=part(elmt_interface)
!print*,'SUCCESS'
return
end subroutine find_interface
!=======================================================

subroutine detect_ghost(out_phead,nelmt,nnode,npart,max_neighbour,glob2loc_elmt)
! this subroutine detect all the ghost partitions/elements/components and write
! in a separate file for each master partition. This routine only works for 8-noded
! hexahedral elements
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
!TODO
! - why there is error for 2 partitions? See the explanation in max_nelmt in the code below
! - some optimization is possible
! - after finding interfacial elements we can only store/process the information 
!   only for those elements
implicit none
integer,parameter :: nenode=8,nedge=12,nface=6 ! number of nodes, edges, and faces per element 
integer,parameter :: nnode_face=4,nnode_edge=2 ! number of nodes per face and per edge

integer,intent(in) :: nelmt,nnode
!integer,dimension(nelmt),intent(in) :: part ! numbering starts from 1 only in this routine
!integer,dimension(esize,nelmt),intent(in) :: connect    
integer,intent(in) :: npart,max_neighbour
integer,dimension(:),intent(in) :: glob2loc_elmt

integer :: i,ig,j,k,i_part,j_part
!integer,dimension(nnode) :: node_npart
integer :: icount !,nelmt_interface
!integer,dimension(:),allocatable :: elmt_interface
integer,dimension(6,4) :: node_face ! local node numbers in each face
integer,dimension(6,4) :: edge_face ! local edge numbers in each face
integer,dimension(12,2) :: node_edge ! local node numbers in each edge
integer :: ielmt,inode,ipart,jelmt,ie,je  
character(len=256),intent(in) :: out_phead
character(len=20) :: format_str
character(len=80) :: out_fname  

integer :: istat
integer,dimension(:),allocatable :: list_gpart
integer :: ind,ind_gpart(1)
integer :: max_nelmt,max_ngpart
integer :: ecomp ! element componnet (1: node, 2: edge, 3: face)
logical,dimension(:,:),allocatable :: mpart_gpart !node_part,
logical :: elmt_face(6),elmt_edge(12),elmt_node(8)

integer :: ngelmt
integer :: gpartid,myindex(1)
integer,allocatable :: nghost(:),adj_ind(:,:),adj_id1(:,:),adj_id2(:,:), &
adj_nelmt(:),adj_elmt(:,:),eadjid(:) !,adj_elmt(:,:),adj1(:,:),adj2(:,:)
integer,dimension(nnode) :: indadj ! i think this can be declared even smaller down to 2*nenode see comments on the code
integer :: ncom,ie_interface,je_interface
integer,allocatable :: gelmt(:),ginterface(:),gadjid(:),adj_iface(:),adj_iedge(:),adj_inode(:)
integer :: adj_nface,adj_nedge,adj_nnode,ia

! ghost partition  
type ghost_partition
integer :: nelmt
!integer,dimension(4) :: gorder,morder
integer,dimension(:),pointer :: melmt,ecomp,meid ! master element,entity type,master element's entity id
integer,dimension(:),pointer :: gelmt,geid ! ghost element,ghost element's entity id
end type ghost_partition

! master partition  
type partition
integer :: ielmt,nelmt
integer,dimension(:),pointer :: elmt
integer,dimension(:),pointer :: eid_interface
integer :: igpart,ngpart
integer,dimension(:),pointer :: gpartid ! ghost partition id
type(ghost_partition),dimension(:),pointer :: gpart ! ghost partition
end type partition
type(partition),dimension(npart) :: mpart ! master partition

!type nodes
!  integer,dimension(:),allocatable :: part
!end type nodes
!type(nodes),dimension(:),allocatable :: node

! local node numbering in each face CUBIT/EXODUS convention  
node_face(1,:)=(/1,2,6,5/)
node_face(2,:)=(/2,3,7,6/)
node_face(3,:)=(/4,3,7,8/)
node_face(4,:)=(/1,4,8,5/)
node_face(5,:)=(/1,2,3,4/)
node_face(6,:)=(/5,6,7,8/)

! local node numbering in each edge CUBIT/EXODUS convention
node_edge(1,:)=(/1,2/)
node_edge(2,:)=(/2,3/)
node_edge(3,:)=(/3,4/)
node_edge(4,:)=(/4,1/)
node_edge(5,:)=(/1,5/)
node_edge(6,:)=(/2,6/)
node_edge(7,:)=(/3,7/)
node_edge(8,:)=(/4,8/)
node_edge(9,:)=(/5,6/)
node_edge(10,:)=(/6,7/)
node_edge(11,:)=(/7,8/)
node_edge(12,:)=(/8,5/)

! edge id in the face
edge_face(1,1)=1; edge_face(1,2)=6; edge_face(1,3)=-9; edge_face(1,4)=-5;
edge_face(2,1)=2; edge_face(2,2)=7; edge_face(2,3)=-10; edge_face(2,4)=-6;
edge_face(3,1)=3; edge_face(3,2)=8; edge_face(3,3)=-11; edge_face(3,4)=-7;
edge_face(4,1)=4; edge_face(4,2)=5; edge_face(4,3)=-12; edge_face(4,4)=-8;
edge_face(5,1)=1; edge_face(5,2)=2; edge_face(5,3)=3; edge_face(5,4)=4;
edge_face(6,1)=9; edge_face(6,2)=10; edge_face(6,3)=11; edge_face(6,4)=12;

!  allocate(node_part(nnode,npart))
!  
!  node_part=.false.  
!  do i=1,nelmt
!    node_part(connect(:,i),part(i))=.true.
!  enddo
!  
!  ! determine the number of partitions for all nodes
!  node_npart=0
!  allocate(node(nnode))
!  do i=1,nnode 
!    node_npart(i)=count(node_part(i,:))
!    allocate(node(i)%part(node_npart(i)))
!    icount=0
!    do j=1,npart
!      if(node_part(i,j))then
!        icount=icount+1
!        node(i)%part(icount)=j
!      endif
!    enddo
!  enddo 
  
!  deallocate(node_part)
!  !print*,minval(node_npart),maxval(node_npart)
!  !stop

!  !node_npart=0
!  ! determine the number of partitions for all nodes
!  !do i_part=1,npart
!  !  do i=1,nelmt
!  !    if(part(i)==i_part)node_npart(connect(:,i))=node_npart(connect(:,i))+1
!  !  enddo
!  !enddo

!  ! count the number of elements in the interface  
!  nelmt_interface=0
!  do i=1,nelmt
!    if(maxval(node_npart(connect(:,i)))>1)nelmt_interface=nelmt_interface+1
!  enddo
!  !max_share=maxval(node_npart)-1
!  print*,'total elements in the interface: ',nelmt_interface
!  !stop

!  ! number of neighboring element for the interface elements
!  allocate(elmt_interface(nelmt_interface))

!  icount=0
!  do i=1,nelmt
!    if(maxval(node_npart(connect(:,i)))>1)then
!      icount=icount+1
!  !    nelmt_share(icount)=maxval(node_npart(connect(:,i)))-1
!      elmt_interface(icount)=i
!    endif
!  enddo

! find the adjacent elements and sharing ID
write(*,'(a)',advance='no')' finding adjacent elements...'
allocate(adj_elmt(nelmt_interface,max_neighbour),adj_ind(nelmt_interface,max_neighbour),adj_nelmt(nelmt_interface), &
adj_id1(nelmt_interface,max_neighbour),adj_id2(nelmt_interface,max_neighbour),stat=istat)

if(istat/=0)then
  print*,'ERROR: insufficient memory!'
  stop
endif  
!adj_elmt=0
indadj=0 ! adjacency indicator, an element value = 2 indicates the adjacency
adj_nelmt=0; adj_elmt=0
adj_ind=0; adj_id1=0; adj_id2=0
do i=1,nelmt_interface-1   
  !ielmt=elmt_interface(i) 
  do j=i+1,nelmt_interface          
    !jelmt=elmt_interface(j)
    
    ! sorting 2 array using the rank as index we can use very small array instead of indadj(nnode)
    ! test
    ! connect_adj(1:nenode)=connect_interface(:,i);
    ! connect_adj(nenode+1:2*nenode)=connect_interface(:,j);
    ! print*,connect_adj
    ! stop
    ! test
    
    ! we are only interested in the part of i and j
    indadj(connect_interface(:,i))=0
    indadj(connect_interface(:,j))=0      
    indadj(connect_interface(:,i))=indadj(connect_interface(:,i))+1
    indadj(connect_interface(:,j))=indadj(connect_interface(:,j))+1
    ncom=sum(indadj(connect_interface(:,i)))-nenode ! number of common nodes 
          
    if(ncom==0)then
      ! no common node
      cycle
    elseif(ncom==1)then
      ! common node
      !adj_elmt(i,j)=1
      !adj_elmt(j,i)=1
      !adj1(i,j)=get_nodeid(indadj(connect_interface(:,i))==2)
      !adj2(i,j)=get_nodeid(indadj(connect_interface(:,j))==2)
      
      adj_nelmt(i)=adj_nelmt(i)+1
      
      adj_id1(i,adj_nelmt(i))=get_nodeid(indadj(connect_interface(:,i))==2)
      adj_id2(i,adj_nelmt(i))=get_nodeid(indadj(connect_interface(:,j))==2)
      
      adj_ind(i,adj_nelmt(i))=1
      adj_elmt(i,adj_nelmt(i))=j
      
      adj_nelmt(j)=adj_nelmt(j)+1
      adj_ind(j,adj_nelmt(j))=1
      adj_elmt(j,adj_nelmt(j))=i
    elseif(ncom==2)then
      ! common edge
      !adj_elmt(i,j)=2
      !adj_elmt(j,i)=2
      !adj1(i,j)=get_edgeid(indadj(connect_interface(:,i))==2,node_edge)
      !adj2(i,j)=get_edgeid(indadj(connect_interface(:,j))==2,node_edge)
      
      adj_nelmt(i)=adj_nelmt(i)+1
      
      adj_id1(i,adj_nelmt(i))=get_edgeid(indadj(connect_interface(:,i))==2,node_edge)
      adj_id2(i,adj_nelmt(i))=get_edgeid(indadj(connect_interface(:,j))==2,node_edge)
      
      adj_ind(i,adj_nelmt(i))=2
      adj_elmt(i,adj_nelmt(i))=j
      
      adj_nelmt(j)=adj_nelmt(j)+1
      adj_ind(j,adj_nelmt(j))=2
      adj_elmt(j,adj_nelmt(j))=i
    elseif(ncom==4)then
      ! common face
      !adj_elmt(i,j)=4
      !adj_elmt(j,i)=4
      !adj1(i,j)=get_faceid(indadj(connect_interface(:,i))==2,node_face)        
      !adj2(i,j)=get_faceid(indadj(connect_interface(:,j))==2,node_face)
      
      adj_nelmt(i)=adj_nelmt(i)+1
      
      adj_id1(i,adj_nelmt(i))=get_faceid(indadj(connect_interface(:,i))==2,node_face)        
      adj_id2(i,adj_nelmt(i))=get_faceid(indadj(connect_interface(:,j))==2,node_face)
      
      adj_ind(i,adj_nelmt(i))=4
      adj_elmt(i,adj_nelmt(i))=j
      
      adj_nelmt(j)=adj_nelmt(j)+1
      adj_ind(j,adj_nelmt(j))=4 
      adj_elmt(j,adj_nelmt(j))=i         
    else
      write(*,*)'ERROR: wrong number of common nodes!',ncom
      stop
    endif      
  enddo
  !print*,i,nelmt_interface
enddo 
write(*,*)'complete!'

allocate(mpart_gpart(npart,npart))
mpart_gpart=.false.
do i=1,nelmt_interface
  !ielmt=elmt_interface(i)
  do j=1,nenode
    inode=connect_interface(j,i)    
    do k=1,node_npart(inode)
      mpart_gpart(part_interface(i),node(inode)%part(k))=.true.
    enddo
  enddo
enddo

do i=1,nnode
  deallocate(node(i)%part)
enddo
deallocate(node,node_npart,connect_interface)

do i_part=1,npart
  mpart_gpart(i_part,i_part)=.false. ! itself cannot be its ghost
enddo
allocate(nghost(npart))
! count ghost partitions of each master partition
do i_part=1,npart
  mpart(i_part)%ngpart=count(mpart_gpart(i_part,:))
  nghost(i_part)=mpart(i_part)%ngpart
enddo
max_ngpart=maxval(mpart(1:npart)%ngpart) ! maximum number of ghost partitions
  
deallocate(mpart_gpart) ! mpart_gpart seems that this can be utilized for efficiency later, but I don't know how to.
! count interfacial elements of each partition
mpart(1:npart)%nelmt=0  
!mpart(1:npart)%ngpart=0
do i=1,nelmt_interface
  !ielmt=elmt_interface(i)
  ipart=part_interface(i)
  ! count elements   
  mpart(ipart)%nelmt=mpart(ipart)%nelmt+1
enddo

max_nelmt=maxval(mpart(1:npart)%nelmt) ! maximum number of element in each partition
! beause there are several repeating elements this number may be larger for the allocation
! therefore it gives bound error for npart=2, we can actually increase this number by some logical values
! to solve the problem. correct solution would be to count the acual sharing entities of all elements!!!

print*,'max ngpart: ',max_ngpart,' max_nelmt: ',max_nelmt

allocate(gelmt(max_nelmt),ginterface(max_nelmt),gadjid(max_nelmt))
allocate(adj_iface(max_nelmt),adj_iedge(max_nelmt),adj_inode(max_nelmt))

! allocate type variable gpart
do i_part=1,npart    
  allocate(mpart(i_part)%elmt(mpart(i_part)%nelmt))
  allocate(mpart(i_part)%eid_interface(mpart(i_part)%nelmt))
enddo

mpart(1:npart)%ielmt=0
! partition interfacial elements
do i=1,nelmt_interface
  ielmt=elmt_interface(i)
  ipart=part_interface(i)
  ! count elements
  mpart(ipart)%ielmt=mpart(ipart)%ielmt+1
  mpart(ipart)%elmt(mpart(ipart)%ielmt)=ielmt
  mpart(ipart)%eid_interface(mpart(ipart)%ielmt)=i
enddo
deallocate(part_interface)
if(.not.all(mpart(1:npart)%ielmt==mpart(1:npart)%nelmt))then
  print*,'ERROR: total number of interfacial elements mismatched!'
  stop
endif

! allocate type variable gpart
do i_part=1,npart    
  allocate(mpart(i_part)%gpart(mpart(i_part)%ngpart))
  allocate(mpart(i_part)%gpartid(mpart(i_part)%ngpart))
  mpart(i_part)%gpartid=0    
  do j_part=1,mpart(i_part)%ngpart
    ! obviously max_nelmt is a waste of memory. there should be some other ways,
    ! but this should not be a problem, because this routine can be run after deallocating
    ! most of the arrays. In case of the elements which shares more than one components
    ! (e.g., some times more than 1 face are shared), there might be segmentation fault.
    ! In this case we can take some larger value than max_nelmt for allocation. However,
    ! this problem should be rare.
    allocate(mpart(i_part)%gpart(j_part)%melmt(max_nelmt))
    allocate(mpart(i_part)%gpart(j_part)%gelmt(max_nelmt))
    allocate(mpart(i_part)%gpart(j_part)%ecomp(max_nelmt))
    allocate(mpart(i_part)%gpart(j_part)%meid(max_nelmt))
    allocate(mpart(i_part)%gpart(j_part)%geid(max_nelmt))
    mpart(i_part)%gpart(j_part)%nelmt=0
  enddo  
enddo
!print*,max_nelmt
allocate(list_gpart(max_ngpart),stat=istat)
if(istat/=0)then
  print*,'ERROR: out of memory!'
  stop
endif

! write master and ghost partitions
write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,a,i'//trim(adjustl(format_str))//',a,i'//trim(adjustl(format_str))//')'
!write(*,*)' '
allocate(eadjid(nelmt_interface))
!inorder=(/1,2,3,4/) ! default local node ordering
mpart(1:npart)%igpart=0 ! counter for unique gpart  
!ipart=0
! find neighbouring entities lying on the different partition
ipart_loop: do i_part=1,npart-1
  write(*,fmt=format_str,advance='no')CR,' matching ghost entities - partition: ',i_part,'/',npart-1
  jpart_loop: do j_part=i_part+1,npart
    !ipart=ipart+1
    !ipart=0
    do i=1,mpart(i_part)%nelmt
      ielmt=mpart(i_part)%elmt(i)
      ie_interface=mpart(i_part)%eid_interface(i)
      elmt_face=.false.; elmt_edge=.false.; elmt_node=.false. 
      
      !print*,adj_elmt(ie_interface,
      ngelmt=mpart(j_part)%nelmt
      gelmt(1:ngelmt)=mpart(j_part)%elmt(:)
      ginterface(1:ngelmt)=mpart(j_part)%eid_interface(:)
      
      gadjid=0
      !print*,size(adjid),ie_interface,minval(ginterface(1:ngelmt)),maxval(ginterface(1:ngelmt))
      eadjid=0 ! elemental adjacency indicators
      eadjid(adj_elmt(ie_interface,1:adj_nelmt(ie_interface)))= &
      adj_ind(ie_interface,1:adj_nelmt(ie_interface))
      !do ig=1,adj_nelmt(ie_interface)
      !  ind=find(ginterface(1:ngelmt)==adjelmt(ie_interface,ig))
      !  if(ind>0)gadjid(ind)=adjind(ie_interface,ig)
      !enddo
      gadjid(1:ngelmt)=eadjid(ginterface(1:ngelmt)) !adj_elmt(ie_interface,ginterface(1:ngelmt))
      !print*,gadjid(1:ngelmt)
      !print*,'hi',adjind(ie_interface,:)
      !print*,'oh0',ginterface(1:ngelmt)
      !print*,'oh1',adjelmt(ie_interface,:)
      !if(all(eadjid(ginterface(1:ngelmt))/=gadjid(1:ngelmt)))then
      !print*,'hello0',eadjid(ginterface(1:ngelmt))
      !print*,'hello1',gadjid(1:ngelmt)
      !endif
      !stop
      ! count faces/edges/nodes
      adj_nface=0; adj_nedge=0; adj_nnode=0
      do ia=1,ngelmt !mpart(j_part)%nelmt
        if(gadjid(ia)==4)then
          adj_nface=adj_nface+1
          adj_iface(adj_nface)=ia
        elseif(gadjid(ia)==2)then
          adj_nedge=adj_nedge+1
          adj_iedge(adj_nedge)=ia
        elseif(gadjid(ia)==1)then
          adj_nnode=adj_nnode+1
          adj_inode(adj_nnode)=ia
        endif
      enddo
      ! list jpart elements
      ! print*,'hi',adj_nface,adj_nedge,adj_nnode
      ! print*,'hi'
      
      ! we check in the order face, edge, and node to exclude redundant sharing (low level repeating) 
      ! check for face/s
      do ia=1,adj_nface
        je_interface=ginterface(adj_iface(ia))
        jelmt=gelmt(adj_iface(ia))
        ecomp=4
        if(ie_interface<je_interface)then
          !ecomp=adj_elmt(ie_interface,je_interface)
          ind=find(adj_elmt(ie_interface,:)==je_interface)
          ie=adj_id1(ie_interface,ind) !adj1(ie_interface,je_interface)
          je=adj_id2(ie_interface,ind) !adj2(ie_interface,je_interface)
        elseif(ie_interface>je_interface)then
          !ecomp=adj_elmt(je_interface,ie_interface)
          ind=find(adj_elmt(je_interface,:)==ie_interface)
          je=adj_id1(je_interface,ind) !adj1(je_interface,ie_interface)
          ie=adj_id2(je_interface,ind) !adj2(je_interface,ie_interface)
        else
          write(*,*)'ERROR: wrong adjacency for face!'
          stop
        endif
        if(elmt_face(ie))cycle
        call set_partition(i_part,j_part,ielmt,jelmt,ecomp,ie,je) !,inorder,jnorder)
        !elmt_face(ie)=.true.
        !elmt_edge(abs(edge_face(ie,:)))=.true.
        !elmt_node(node_face(ie,:))=.true.
        
        ! ielmt is the ghost of jelmt                
        call set_partition(j_part,i_part,jelmt,ielmt,ecomp,je,ie) !,jnorder,inorder)
        elmt_face(ie)=.true.
        elmt_edge(abs(edge_face(ie,:)))=.true.
        elmt_node(node_face(ie,:))=.true.
      enddo
      
      ! chek for edge/s
      do ia=1,adj_nedge
        je_interface=ginterface(adj_iedge(ia))
        jelmt=gelmt(adj_iedge(ia))
        ecomp=2
        if(ie_interface<je_interface)then          
          !ecomp=adj_elmt(ie_interface,je_interface)
          ind=find(adj_elmt(ie_interface,:)==je_interface)
          ie=adj_id1(ie_interface,ind) !adj1(ie_interface,je_interface)
          je=adj_id2(ie_interface,ind) !adj2(ie_interface,je_interface)
        elseif(ie_interface>je_interface)then
          !ecomp=adj_elmt(je_interface,ie_interface)
          ind=find(adj_elmt(je_interface,:)==ie_interface)
          je=adj_id1(je_interface,ind) !adj1(je_interface,ie_interface)
          ie=adj_id2(je_interface,ind) !adj2(je_interface,ie_interface)
        else
          write(*,*)'ERROR: wrong adjacency for edge!'
          stop
        endif
        if(elmt_edge(ie))cycle
        call set_partition(i_part,j_part,ielmt,jelmt,ecomp,ie,je) !,inorder,jnorder)
        !elmt_face(ie)=.true.
        !elmt_edge(abs(edge_face(ie,:)))=.true.
        !elmt_node(node_face(ie,:))=.true.
        
        ! ielmt is the ghost of jelmt                
        call set_partition(j_part,i_part,jelmt,ielmt,ecomp,je,ie) !,jnorder,inorder)
        elmt_edge(ie)=.true.
        elmt_node(node_edge(ie,:))=.true.
      enddo
      
      ! check for node/s
      do ia=1,adj_nnode  
        je_interface=ginterface(adj_inode(ia))
        jelmt=gelmt(adj_inode(ia))
        ecomp=1
        if(ie_interface<je_interface)then
          !ecomp=adj_elmt(ie_interface,je_interface)
          ind=find(adj_elmt(ie_interface,:)==je_interface)
          ie=adj_id1(ie_interface,ind) !adj1(ie_interface,je_interface)
          je=adj_id2(ie_interface,ind) !adj2(ie_interface,je_interface)
        elseif(ie_interface>je_interface)then
          !ecomp=adj_elmt(je_interface,ie_interface)
          ind=find(adj_elmt(je_interface,:)==ie_interface)
          je=adj_id1(je_interface,ind) !adj1(je_interface,ie_interface)
          ie=adj_id2(je_interface,ind) !adj2(je_interface,ie_interface)
        else
          write(*,*)'ERROR: wrong adjacency for node!'
          stop
        endif
        if(elmt_node(ie))cycle
        call set_partition(i_part,j_part,ielmt,jelmt,ecomp,ie,je) !,inorder,jnorder)
                  
        ! ielmt is the ghost of jelmt                
        call set_partition(j_part,i_part,jelmt,ielmt,ecomp,je,ie) !,jnorder,inorder)
        elmt_node(ie)=.true.
      enddo
    enddo ! i
    !print*,'total:',ipart
    !print*,mpart(i_part)%nelmt,mpart(j_part)%nelmt
    !stop
  enddo jpart_loop ! j_part
enddo ipart_loop ! i_part  
!print*,'total',ipart
if(.not.all(mpart(1:npart)%igpart==mpart(1:npart)%ngpart))then
  print*,'ERROR: total number of ghost partitions mismatched!'    
  print*,'Total number of ghost partitions'
  print*,mpart(1:npart)%ngpart
  print*,'Counted number of ghost partitions'
  print*,mpart(1:npart)%igpart
  stop
endif
deallocate(adj_elmt,adj_id1,adj_id2,adj_ind,adj_nelmt,eadjid)
deallocate(elmt_interface)
deallocate(list_gpart)

write(*,*)' complete!'
write(*,'(a)',advance='no')' saving ghost files...'

! write master and ghost partitions
write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'  

! write interfaces in each partition
do i_part=1,npart    

  ! open output file
  write(out_fname, fmt=format_str)trim(out_phead)//'ghost_proc',i_part-1
  open(unit=16,file=trim(out_fname),status='replace',action='write',iostat=istat)
  if(istat/=0) then
    print*,'error file open:',trim(out_fname)
    stop
  endif
  write(16,*)'this partition'
  write(16,*)i_part-1
  write(16,*)'number of ghost partitions, max ngpart, max number of elements'
  !print*,mpart(i_part)%gpart(1:mpart(i_part)%ngpart)%nelmt
  !print*,'max',maxval(mpart(i_part)%gpart(1:mpart(i_part)%ngpart)%nelmt)
  write(16,*)mpart(i_part)%ngpart,max_ngpart,maxval(mpart(i_part)%gpart(1:mpart(i_part)%ngpart)%nelmt)
  do j_part=1,mpart(i_part)%ngpart
    gpartid=mpart(i_part)%gpartid(j_part)
    myindex=maxloc(mpart(gpartid)%gpartid,logical(mpart(gpartid)%gpartid==i_part,1))
    if(myindex(1)==0)then
      print*,'ERROR: master partition not found in its ghost partition!'
      stop
    endif
    write(16,*)'ghost partition ',j_part,' , my index there'
    write(16,*)mpart(i_part)%gpartid(j_part)-1,myindex(1)   
    write(16,*)'number of ghost elements in this ghost partition'     
    !ngelmt=mpart(i_part)%gpart(j_part)%nelmt
    write(16,*)mpart(i_part)%gpart(j_part)%nelmt !,sum(mpart(i_part)%gpart(j_part)%ecomp(1:ngelmt))
    do i=1,mpart(i_part)%gpart(j_part)%nelmt
      ecomp=mpart(i_part)%gpart(j_part)%ecomp(i)        
      write(16,*)glob2loc_elmt(mpart(i_part)%gpart(j_part)%melmt(i)),mpart(i_part)%gpart(j_part)%ecomp(i), &
        mpart(i_part)%gpart(j_part)%meid(i) !,mpart(i_part)%gpart(j_part)%morder(1:ecomp), &
        !glob2loc_elmt(mpart(i_part)%gpart(j_part)%gelmt(i)), &
        !mpart(i_part)%gpart(j_part)%geid(i),mpart(i_part)%gpart(j_part)%gorder(1:ecomp)
    enddo      
  enddo
  close(16)

enddo
write(*,*)' complete!' 

! deallocate type variable gpart
do i_part=1,npart       
  do j_part=1,mpart(i_part)%ngpart
    deallocate(mpart(i_part)%gpart(j_part)%melmt)
    deallocate(mpart(i_part)%gpart(j_part)%gelmt)
    deallocate(mpart(i_part)%gpart(j_part)%ecomp)
    deallocate(mpart(i_part)%gpart(j_part)%meid)
    deallocate(mpart(i_part)%gpart(j_part)%geid)
    
  enddo  
enddo
do i_part=1,npart    
  deallocate(mpart(i_part)%gpart)
  deallocate(mpart(i_part)%gpartid)
  deallocate(mpart(i_part)%elmt)
  deallocate(mpart(i_part)%eid_interface)
enddo
!deallocate(mpart)
!=======================================================
  
contains
!--------------------------------------------------
! this function finds the node ID of a true node in isnode
function find(isnode) result(ind)
implicit none
logical,intent(in) :: isnode(:)
integer :: i,ind,nlen
nlen=size(isnode)
ind=0
do i=1,nlen
  if(isnode(i))then
    ind=i
    return
  endif
enddo  
end function find
!--------------------------------------------------

! this function finds the node ID of a true node in isnode
function get_nodeid(isnode) result(id)
implicit none
logical,dimension(nenode),intent(in) :: isnode
integer :: i,id
do i=1,nenode
  if(isnode(i))then
    id=i
    return
  endif
enddo
write(*,*)'ERROR: no common node for node ID!'
stop   
end function get_nodeid
!--------------------------------------------------

! this function finds the edge ID according to node_edge of a set of true nodes in isnode
function get_edgeid(isnode,node_edge) result(id)
implicit none
logical,dimension(nenode),intent(in) :: isnode
integer :: i,id
integer,dimension(nedge,2),intent(in) :: node_edge ! local node numbers in each edge

do i=1,nedge
  if(all(isnode(node_edge(i,:))))then
    id=i
    return
  endif
enddo
write(*,*)'ERROR: no common node for edge ID!'
stop   
end function get_edgeid
!-------------------------------------------------- 
! this function finds the face ID according to node_face of a set of true nodes in isnode
function get_faceid(isnode,node_face) result(id)
implicit none
logical,dimension(nenode),intent(in) :: isnode
integer :: i,id
integer,dimension(nface,4),intent(in) :: node_face ! local node numbers in each face

do i=1,nface
  if(all(isnode(node_face(i,:))))then
    id=i
    return
  endif
enddo
write(*,*)'ERROR: no common node for face ID!'
stop   
end function get_faceid
!--------------------------------------------------   

! this subroutine set all the members of each master-ghost partition pair
subroutine set_partition(master,ghost,melmt,gelmt,ecomp,meid,geid) !,morder,gorder)

implicit none
integer,intent(in) :: master,ghost ! master and ghost partitions
integer,intent(in) :: melmt,gelmt,ecomp,meid,geid
!integer,dimension(4),intent(in) :: morder,gorder

! find/determine the index of the ghost partition
if (mpart(master)%igpart==0)then ! first time count
  mpart(master)%igpart=1
  ind_gpart(1)=1
  mpart(master)%gpartid(1)=ghost        
else ! already counted
  list_gpart=-1    
  list_gpart(1:mpart(master)%igpart)=mpart(master)%gpartid(1:mpart(master)%igpart) ! list of existing ghost partitions        
  ind_gpart=maxloc(list_gpart,logical(list_gpart==ghost,1)) ! search if this ghost partion is already in the list and find the index
  
  if(ind_gpart(1)==0)then ! not found, i.e., this is a new gpart         
    mpart(master)%igpart=mpart(master)%igpart+1
    ind_gpart(1)=mpart(master)%igpart
    mpart(master)%gpartid(ind_gpart(1))=ghost                
  endif
endif

! count ghost elements in the ghost partition
  
!if(master==1 .and. ghost==2 )print*,'first0 nelmt: ',mpart(master)%gpart(ind_gpart(1))%nelmt,ind_gpart(1), &
!mpart(master)%gpartid(ind_gpart(1))
!print*,'hi0',mpart(master)%gpart(ind_gpart(1))%nelmt
!stop
mpart(master)%gpart(ind_gpart(1))%nelmt=mpart(master)%gpart(ind_gpart(1))%nelmt+1         
ind=mpart(master)%gpart(ind_gpart(1))%nelmt
!if(master==1 .and. ghost==2 )print*,'first1 nelmt: ',mpart(master)%gpart(ind_gpart(1))%nelmt,ind_gpart(1), &
!mpart(master)%gpartid(ind_gpart(1))

! set all ghost partition members
!print*,'hi1',master,ind_gpart(1),ind
mpart(master)%gpart(ind_gpart(1))%melmt(ind)=melmt
mpart(master)%gpart(ind_gpart(1))%gelmt(ind)=gelmt
mpart(master)%gpart(ind_gpart(1))%ecomp(ind)=ecomp
mpart(master)%gpart(ind_gpart(1))%meid(ind)=meid
mpart(master)%gpart(ind_gpart(1))%geid(ind)=geid
!mpart(master)%gpart(ind_gpart(1))%morder=morder
!mpart(master)%gpart(ind_gpart(1))%gorder=gorder
end subroutine set_partition
!--------------------------------------------------  

end subroutine detect_ghost
!=======================================================

!--------------------------------------------------
! Construct interfaces between each partitions.
! Two adjacent elements in distinct partitions make an entry in array tab_interfaces :
! 1/ first element, 2/ second element, 3/ number of common nodes, 4/ first node,
! 5/ second node, if relevant.

! interface ignores acoustic and elastic elements 

! Elements with undefined material are considered as elastic elements.
!--------------------------------------------------
subroutine Construct_interfaces(nelmnts,sup_neighbour,part, &
elmnts,xadj,adjncy,tab_interfaces, tab_size_interfaces, &
ninterfaces,npart)

integer,intent(in) :: nelmnts,sup_neighbour !integer(long)
integer,dimension(0:nelmnts-1),intent(in) :: part
integer,dimension(0:esize*nelmnts-1),intent(in) :: elmnts
integer,dimension(0:nelmnts),intent(in) :: xadj
integer,dimension(0:sup_neighbour*nelmnts-1),intent(in) :: adjncy
integer,dimension(:),pointer :: tab_size_interfaces, tab_interfaces
integer,intent(out) :: ninterfaces

integer,intent(in) :: npart

! local parameters  
integer :: num_part,num_part_bis,el,el_adj,num_interface, &
num_edge,ncommon_nodes,num_node, num_node_bis
integer :: i,j

! counts number of interfaces between partitions
ninterfaces = 0
do  i = 0, npart-1
    do j = i+1, npart-1
      ninterfaces = ninterfaces + 1
    end do
end do

print*,npart,ninterfaces

allocate(tab_size_interfaces(0:ninterfaces))
tab_size_interfaces(:) = 0

num_interface = 0
num_edge = 0

! determines acoustic/elastic elements based upon given vs velocities
! and counts same elements for each interface
do num_part = 0, npart-1
    do num_part_bis = num_part+1, npart-1
      do el = 0, nelmnts-1
          if ( part(el) == num_part ) then                
            ! looks at all neighbor elements 
            do el_adj = xadj(el), xadj(el+1)-1
                ! adds element if neighbor element lies in next partition
                if ( part(adjncy(el_adj)) == num_part_bis ) then
                  num_edge = num_edge + 1
                end if

            end do
          end if
      end do
      ! stores number of elements at interface
      tab_size_interfaces(num_interface+1) = tab_size_interfaces(num_interface) + num_edge
      num_edge = 0
      num_interface = num_interface + 1

    end do
end do


! stores element indices for elements from above search at each interface
num_interface = 0
num_edge = 0

allocate(tab_interfaces(0:(tab_size_interfaces(ninterfaces)*7-1)))
tab_interfaces(:) = 0

do num_part = 0, npart-1
    do num_part_bis = num_part+1, npart-1
      do el = 0, nelmnts-1
          if ( part(el) == num_part ) then
            do el_adj = xadj(el), xadj(el+1)-1
                ! adds element if in adjacent partition                    
                if ( part(adjncy(el_adj)) == num_part_bis ) then
                  tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+0) = el
                  tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+1) = adjncy(el_adj)
                  ncommon_nodes = 0
                  do num_node = 0, esize-1
                      do num_node_bis = 0, esize-1
                        if ( elmnts(el*esize+num_node) == elmnts(adjncy(el_adj)*esize+num_node_bis) ) then
                            tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+3+ncommon_nodes) &
                                = elmnts(el*esize+num_node)
                            ncommon_nodes = ncommon_nodes + 1
                        end if
                      end do
                  end do
                  if ( ncommon_nodes > 0 ) then
                      tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+2) = ncommon_nodes
                  else
                      print *, "Error while building interfaces!", ncommon_nodes
                  end if
                  num_edge = num_edge + 1
                end if
            end do
          end if

      end do
      num_edge = 0
      num_interface = num_interface + 1
    end do
end do

end subroutine Construct_interfaces
!=======================================================

!--------------------------------------------------
! Construct interfaces between each partitions.
! Two adjacent elements in distinct partitions make an entry in array tab_interfaces :
! 1/ first element, 2/ second element, 3/ number of common nodes, 4/ first node,
! 5/ second node, if relevant.

! No interface between acoustic and elastic elements.

! Elements with undefined material are considered as elastic elements.
!--------------------------------------------------
subroutine Construct_interfaces_no_ac_el_sep(nelmnts, &
sup_neighbour,part,elmnts,xadj,adjncy,tab_interfaces, &
tab_size_interfaces,ninterfaces,nb_materials,cs_material, &
num_material,npart)

integer,intent(in) :: nelmnts,sup_neighbour !integer(long)
integer,dimension(0:nelmnts-1),intent(in) :: part
integer,dimension(0:esize*nelmnts-1),intent(in) :: elmnts
integer,dimension(0:nelmnts),intent(in) :: xadj
integer,dimension(0:sup_neighbour*nelmnts-1),intent(in) :: adjncy
integer,dimension(:),pointer :: tab_size_interfaces, tab_interfaces
integer,intent(out) :: ninterfaces
integer,dimension(1:nelmnts),intent(in) :: num_material
integer,intent(in) :: nb_materials,npart
! vs velocities
real(kind=kreal),dimension(1:nb_materials),intent(in) :: cs_material !double precision



! local parameters  
integer :: num_part,num_part_bis,el,el_adj,num_interface, &
num_edge,ncommon_nodes,num_node,num_node_bis
integer :: i,j
logical :: is_acoustic_el,is_acoustic_el_adj

! counts number of interfaces between partitions
ninterfaces = 0
do  i = 0, npart-1
    do j = i+1, npart-1
      ninterfaces = ninterfaces + 1
    end do
end do

allocate(tab_size_interfaces(0:ninterfaces))
tab_size_interfaces(:) = 0

num_interface = 0
num_edge = 0

! determines acoustic/elastic elements based upon given vs velocities
! and counts same elements for each interface
do num_part = 0, npart-1
    do num_part_bis = num_part+1, npart-1
      do el = 0, nelmnts-1
          if ( part(el) == num_part ) then
            ! determines whether element is acoustic or not
            if(num_material(el+1) > 0) then
                if ( cs_material(num_material(el+1)) < TINYVAL) then
                  is_acoustic_el = .true.
                else
                  is_acoustic_el = .false.
                end if
            else
                is_acoustic_el = .false.
            end if
            ! looks at all neighbor elements 
            do el_adj = xadj(el), xadj(el+1)-1
                ! determines whether neighbor element is acoustic or not
                if(num_material(adjncy(el_adj)+1) > 0) then
                  if ( cs_material(num_material(adjncy(el_adj)+1)) < TINYVAL) then
                      is_acoustic_el_adj = .true.
                  else
                      is_acoustic_el_adj = .false.
                  end if
                else
                  is_acoustic_el_adj = .false.
                end if
                ! adds element if neighbor element has same material acoustic/not-acoustic and lies in next partition
                if ( (part(adjncy(el_adj)) == num_part_bis) .and. (is_acoustic_el .eqv. is_acoustic_el_adj) ) then
                  num_edge = num_edge + 1
                end if
            end do
          end if
      end do
      ! stores number of elements at interface
      tab_size_interfaces(num_interface+1) = tab_size_interfaces(num_interface) + num_edge
      num_edge = 0
      num_interface = num_interface + 1

    end do
end do


! stores element indices for elements from above search at each interface
num_interface = 0
num_edge = 0

allocate(tab_interfaces(0:(tab_size_interfaces(ninterfaces)*7-1)))
tab_interfaces(:) = 0

do num_part = 0, npart-1
    do num_part_bis = num_part+1, npart-1
      do el = 0, nelmnts-1
          if ( part(el) == num_part ) then
            if(num_material(el+1) > 0) then
                if ( cs_material(num_material(el+1)) < TINYVAL) then
                  is_acoustic_el = .true.
                else
                  is_acoustic_el = .false.
                end if
            else
                is_acoustic_el = .false.
            end if
            do el_adj = xadj(el), xadj(el+1)-1
                if(num_material(adjncy(el_adj)+1) > 0) then
                  if ( cs_material(num_material(adjncy(el_adj)+1)) < TINYVAL) then
                      is_acoustic_el_adj = .true.
                  else
                      is_acoustic_el_adj = .false.
                  end if
                else
                  is_acoustic_el_adj = .false.
                end if
                if ( (part(adjncy(el_adj)) == num_part_bis) .and. (is_acoustic_el .eqv. is_acoustic_el_adj) ) then
                  tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+0) = el
                  tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+1) = adjncy(el_adj)
                  ncommon_nodes = 0
                  do num_node = 0, esize-1
                      do num_node_bis = 0, esize-1
                        if ( elmnts(el*esize+num_node) == elmnts(adjncy(el_adj)*esize+num_node_bis) ) then
                            tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+3+ncommon_nodes) &
                                = elmnts(el*esize+num_node)
                            ncommon_nodes = ncommon_nodes + 1
                        end if
                      end do
                  end do
                  if ( ncommon_nodes > 0 ) then
                      tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+2) = ncommon_nodes
                  else
                      print *, "Error while building interfaces!", ncommon_nodes
                  end if
                  num_edge = num_edge + 1
                end if
            end do
          end if

      end do
      num_edge = 0
      num_interface = num_interface + 1
    end do
end do
end subroutine Construct_interfaces_no_ac_el_sep
!=======================================================


!--------------------------------------------------
! Write nodes (their coordinates) pertaining to iproc partition in the corresponding Database
!--------------------------------------------------
subroutine write_glob2loc_nodes_database(out_path,xfile,yfile,zfile,iproc,npgeo, &
nodes_coords,glob2loc_nodes_npart,glob2loc_nodes_parts, &
glob2loc_nodes,nnodes,npart)

!integer,intent(in) :: IIN_database
integer,intent(in) :: nnodes,iproc,npart !, num_phase
integer,intent(inout) :: npgeo

real(kind=kreal),dimension(3,nnodes) :: nodes_coords !double precision
integer,dimension(:),pointer :: glob2loc_nodes_npart
integer,dimension(:),pointer :: glob2loc_nodes_parts
integer,dimension(:),pointer :: glob2loc_nodes

integer :: i,istat,j

character(len=256),intent(in) :: out_path
character(len=80),intent(in) :: xfile,yfile,zfile
character(len=20) :: format_str
character(len=80) :: out_fname

write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
!print*,format_str
!print*,trim(out_path)//trim(xfile)//'_proc',iproc
! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(xfile)//'_proc',iproc
open(unit=16,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(yfile)//'_proc',iproc
open(unit=17,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(zfile)//'_proc',iproc
open(unit=18,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

!if ( num_phase == 1 ) then
! counts number of points in partition
    npgeo = 0
    do i = 0, nnodes-1
      do j = glob2loc_nodes_npart(i), glob2loc_nodes_npart(i+1)-1
          if ( glob2loc_nodes_parts(j) == iproc ) then
            npgeo = npgeo + 1

          end if

      end do
    end do
write(16,*)npgeo
write(17,*)npgeo
write(18,*)npgeo

!write(out_fname, fmt=format_str)trim(out_path)//'gnode_proc',iproc
!open(unit=19,file=trim(out_fname),&
!  status='unknown', action='write', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'error file open:',trim(out_fname)
!  stop
!endif
!write(19,*)npgeo
!else
! writes out point coordinates
    do i = 0, nnodes-1
      do j = glob2loc_nodes_npart(i), glob2loc_nodes_npart(i+1)-1
          if ( glob2loc_nodes_parts(j) == iproc ) then
            write(16,'(i10,f25.15)') glob2loc_nodes(j)+1, nodes_coords(1,i+1)
            write(17,'(i10,f25.15)') glob2loc_nodes(j)+1, nodes_coords(2,i+1)
            write(18,'(i10,f25.15)') glob2loc_nodes(j)+1, nodes_coords(3,i+1)
!            write(19,'(i10)')i+1
          end if
      end do
    end do
!end if
close(16)
close(17)
close(18)
!close(19)

end subroutine Write_glob2loc_nodes_database
!=======================================================

!--------------------------------------------------
! Write material properties in the Database
!--------------------------------------------------
subroutine write_material_properties_database(out_path,matfile, &
count_def_mat,count_undef_mat,mat_domain,mat_prop, &
undef_mat_domain,undef_mat_prop,nwmat,waterid,iproc,npart) 

integer,intent(in) :: count_def_mat,count_undef_mat,iproc,npart
integer,dimension(count_def_mat) :: mat_domain
real(kind=kreal),dimension(6,count_def_mat)  :: mat_prop !double precision
integer,dimension(count_undef_mat) :: undef_mat_domain
character(len=30),dimension(6,count_undef_mat) :: undef_mat_prop
integer,intent(in) :: nwmat
integer,intent(in) :: waterid(nwmat)
integer :: i,istat
character(len=256),intent(in) :: out_path
character(len=80),intent(in) :: matfile
character(len=20) :: format_str
character(len=80) :: out_fname

write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(matfile)//'_proc',iproc
open(unit=16,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

write(16,*)'# material properties (domain,id,gamma,ym,nu,phi,coh,psi)'   
write(16,*)  count_def_mat !,count_undef_mat 
do i = 1, count_def_mat
  ! database material definition
  !
  ! format:  #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id     
  !
  ! (note that this order of the properties is different than the input in nummaterial_velocity_file)
  !
    write(16,*) i,mat_domain(i),mat_prop(1,i), mat_prop(2,i), mat_prop(3,i), &
                        mat_prop(4,i), mat_prop(5,i), mat_prop(6,i)
end do
!do i = 1, count_undef_mat
!   write(IIN_database,*) undef_mat_domain(i),trim(undef_mat_prop(1,i)),trim(undef_mat_prop(2,i)), &
!                        trim(undef_mat_prop(3,i)),trim(undef_mat_prop(4,i)), &
!                        trim(undef_mat_prop(5,i)),trim(undef_mat_prop(6,i)),trim(undef_mat_prop(7,i))
!end do
! write water properties if any
write(16,*)nwmat
do i=1,nwmat
  !print*,nwmat,i,waterid(i)      
  write(16,*)waterid(i) 
enddo
close(16)

end subroutine  write_material_properties_database
!===========================================

!--------------------------------------------------
! Write elements on boundaries (and their four nodes on boundaries) pertaining to iproc partition in the corresponding Database
!--------------------------------------------------
subroutine write_boundaries_database(out_path,uxfile,uyfile,uzfile,iproc, nelmnts, nspec2D_xmin, nspec2D_xmax, &
                      nspec2D_ymin, nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                      ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                      ibelm_ymax, ibelm_bottom, ibelm_top, &
                      nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                      nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, & 
                      glob2loc_elmnts, glob2loc_nodes_npart, &
                      glob2loc_nodes_parts, glob2loc_nodes, part,npart )
    
integer,intent(in) :: iproc,npart
integer,intent(in) :: nelmnts !integer(long)
integer,intent(in) :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin, nspec2D_ymax,nspec2D_bottom,nspec2D_top
integer,dimension(nspec2D_xmin),intent(in) :: ibelm_xmin
integer,dimension(nspec2D_xmax),intent(in) :: ibelm_xmax
integer,dimension(nspec2D_ymin),intent(in) :: ibelm_ymin
integer,dimension(nspec2D_ymax),intent(in) :: ibelm_ymax
integer,dimension(nspec2D_bottom),intent(in) :: ibelm_bottom
integer,dimension(nspec2D_top),intent(in) :: ibelm_top 

integer,dimension(nspec2D_xmin),intent(in) :: nodes_ibelm_xmin
integer,dimension(nspec2D_xmax),intent(in) :: nodes_ibelm_xmax
integer,dimension(nspec2D_ymin),intent(in) :: nodes_ibelm_ymin
integer,dimension(4,nspec2D_ymax),intent(in) :: nodes_ibelm_ymax
integer,dimension(4,nspec2D_bottom),intent(in) :: nodes_ibelm_bottom
integer,dimension(4,nspec2D_top),intent(in) :: nodes_ibelm_top    
integer,dimension(:),pointer :: glob2loc_elmnts
integer,dimension(:),pointer :: glob2loc_nodes_npart
integer,dimension(:),pointer :: glob2loc_nodes_parts
integer,dimension(:),pointer :: glob2loc_nodes
integer,dimension(1:nelmnts) :: part

! local parameters  
integer :: i
integer :: loc_node1,loc_node2,loc_node3,loc_node4
integer :: loc_nspec2D_xmin,loc_nspec2D_xmax,loc_nspec2D_ymin, &
            loc_nspec2D_ymax,loc_nspec2D_bottom,loc_nspec2D_top

integer :: istat
!character(len=80) :: prname
character(len=256) :: out_path,uxfile,uyfile,uzfile
character(len=20) :: format_str
character(len=80) :: out_fname

write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(uxfile)//'_proc',iproc
open(unit=16,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

! counts number of elements for boundary at xmin, xmax, ymin, ymax, bottom, top in this partition
loc_nspec2D_xmin = 0
do i=1,nspec2D_xmin  
    if(part(ibelm_xmin(i)) == iproc) then
      loc_nspec2D_xmin = loc_nspec2D_xmin + 1
    end if
end do
write(16,*)loc_nspec2D_xmin

! outputs element index and element node indices
! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
!          (this is assigned by CUBIT, if this changes the following indexing must be changed as well)
!          while glob2loc_elmnts(.) is shifted from 0 to nspec-1  thus 
!          we need to have the arg of glob2loc_elmnts start at 0 ==> glob2loc_nodes(ibelm_** -1 )
do i=1,nspec2D_xmin  
    if(part(ibelm_xmin(i)) == iproc) then
      write(16,*) glob2loc_elmnts(ibelm_xmin(i)-1)+1, nodes_ibelm_xmin(i)  
    end if
end do
close(16)

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(uyfile)//'_proc',iproc
open(unit=16,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

loc_nspec2D_xmax = 0
do i=1,nspec2D_xmax  
    if(part(ibelm_xmax(i)) == iproc) then
      loc_nspec2D_xmax = loc_nspec2D_xmax + 1
    end if
end do
write(16,*)loc_nspec2D_xmax

do i=1,nspec2D_xmax     
    if(part(ibelm_xmax(i)) == iproc) then         
      write(16,*) glob2loc_elmnts(ibelm_xmax(i)-1)+1, nodes_ibelm_xmax(i)  
    end if
end do
close(16)

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(uzfile)//'_proc',iproc
open(unit=16,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif
loc_nspec2D_ymin = 0
do i=1,nspec2D_ymin  
    if(part(ibelm_ymin(i)) == iproc) then
      loc_nspec2D_ymin = loc_nspec2D_ymin + 1
    end if
end do
write(16,*)loc_nspec2D_ymin   

do i=1,nspec2D_ymin     
    if(part(ibelm_ymin(i)) == iproc) then          
      write(16,*) glob2loc_elmnts(ibelm_ymin(i)-1)+1, nodes_ibelm_ymin(i)  
    end if
end do
close(16)   

end subroutine write_boundaries_database
!===========================================

!--------------------------------------------------
! Write elements (their nodes) pertaining to iproc partition in the corresponding Database
!--------------------------------------------------
subroutine write_partition_database(out_path,confile,idfile,iproc, nspec, nelmnts, elmnts, &
  glob2loc_elmnts, glob2loc_nodes_npart,glob2loc_nodes_parts, glob2loc_nodes, &
  part,matid,ngnod,npart)

!    include './constants_decompose_mesh_SCOTCH.h'

!integer,intent(in) :: IIN_database
integer,intent(in) :: iproc,npart !num_phase
integer,intent(in) :: nelmnts !
integer,intent(inout) :: nspec
integer,dimension(0:nelmnts-1) :: part
integer,dimension(0:esize*nelmnts-1) :: elmnts
integer,dimension(:),pointer :: glob2loc_elmnts
integer,dimension(2,nelmnts) :: matid
integer,dimension(:),pointer :: glob2loc_nodes_npart
integer,dimension(:),pointer :: glob2loc_nodes_parts
integer,dimension(:),pointer :: glob2loc_nodes
integer,intent(in) :: ngnod

integer :: i,istat,j,k
integer,dimension(0:ngnod-1) :: loc_nodes
character(len=256),intent(in) :: out_path
character(len=80),intent(in) :: confile,idfile
character(len=20) :: format_str
character(len=80) :: out_fname

write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(confile)//'_proc',iproc
open(unit=16,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

! open output file
write(out_fname, fmt=format_str)trim(out_path)//trim(idfile)//'_proc',iproc
open(unit=17,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif



!if ( num_phase == 1 ) then
! counts number of spectral elements in this partition
    nspec = 0
    do i = 0, nelmnts-1
      if ( part(i) == iproc ) then
          nspec = nspec + 1
      end if
    end do
    write(16,*)nspec
    write(17,*)nspec
    
    !write(out_fname, fmt=format_str)trim(out_phead)//'elmt_proc',iproc
    ! open(unit=19,file=trim(out_fname),&
    !   status='unknown', action='write', form='formatted', iostat = istat)
    ! if( istat /= 0 ) then
    !   print*,'error file open:',trim(out_fname)
    !   stop
    ! endif
    ! write(19,*)nspec

!else
! writes out element corner indices
    do i = 0, nelmnts-1
      if ( part(i) == iproc ) then

          do j = 0, ngnod-1
            do k = glob2loc_nodes_npart(elmnts(i*ngnod+j)), glob2loc_nodes_npart(elmnts(i*ngnod+j)+1)-1

                if ( glob2loc_nodes_parts(k) == iproc ) then
                  loc_nodes(j) = glob2loc_nodes(k)
                end if
            end do

          end do
          
          ! format:
          ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id8
          !write(16,*) glob2loc_elmnts(i)+1, matid(1,i+1), matid(2,i+1),(loc_nodes(k)+1, k=0,ngnod-1)
          write(16,*) glob2loc_elmnts(i)+1, (loc_nodes(k)+1, k=0,ngnod-1)
          write(17,*) glob2loc_elmnts(i)+1, matid(2,i+1)
          !write(19,'(i10)')i+1
      end if
    end do
!end if
close(16)
close(17)
!close(19)
end subroutine write_partition_database
!=======================================================


!--------------------------------------------------
! Write interfaces (element and common nodes) pertaining to iproc partition in the corresponding Database
!--------------------------------------------------
subroutine write_interfaces_database(out_phead,tab_interfaces, tab_size_interfaces, iproc, ninterfaces, &
      my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, glob2loc_nodes_npart, glob2loc_nodes_parts, &
      glob2loc_nodes, npart)

!    include './constants_decompose_mesh_SCOTCH.h'

!   integer,intent(in) :: IIN_database
integer,intent(in) :: iproc
integer,intent(in) :: ninterfaces,npart
integer,intent(inout) :: my_ninterface
integer,dimension(:),pointer :: tab_size_interfaces
integer,dimension(:),pointer :: tab_interfaces
integer,dimension(0:ninterfaces-1),intent(inout)  :: my_interfaces
integer,dimension(0:ninterfaces-1),intent(inout)  :: my_nb_interfaces
integer,dimension(:),pointer :: glob2loc_elmnts
integer,dimension(:),pointer :: glob2loc_nodes_npart
integer,dimension(:),pointer :: glob2loc_nodes_parts
integer,dimension(:),pointer :: glob2loc_nodes

integer,dimension(4) :: local_nodes
integer :: local_elmnt
!integer :: num_phase

integer :: i,j,k,l
integer :: num_interface

integer :: count_faces
integer :: istat
!character(len=80) :: prname
character(len=256) :: out_phead
character(len=20) :: format_str
character(len=80) :: out_fname    

write(format_str,*)ceiling(log10(real(npart)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

! open output file
write(out_fname, fmt=format_str)trim(out_phead)//'interface_proc',iproc
open(unit=16,file=trim(out_fname),&
  status='unknown', action='write', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'error file open:',trim(out_fname)
  stop
endif

num_interface = 0

!if ( num_phase == 1 ) then
! counts number of interfaces to neighbouring partitions
    my_interfaces(:) = 0
    my_nb_interfaces(:) = 0
  
    ! double loops over all partitions
    do i = 0, npart-1       
      do j = i+1, npart-1
          ! only counts if specified partition (iproc) appears and interface elements increment
          if ( (tab_size_interfaces(num_interface) < tab_size_interfaces(num_interface+1)) .and. &
              (i == iproc .or. j == iproc) ) then
            ! sets flag  
            my_interfaces(num_interface) = 1
            ! sets number of elements on interface
            my_nb_interfaces(num_interface) = tab_size_interfaces(num_interface+1) - tab_size_interfaces(num_interface)
          end if
          num_interface = num_interface + 1
      end do
    end do
    my_ninterface = sum(my_interfaces(:))

    ! writes out MPI interfaces elements
    write(16,*) my_ninterface, maxval(my_nb_interfaces)
    !print*,'hi'

!else
! writes out MPI interface elements
  do i = 0, npart-1
      do j = i+1, npart-1
        if ( my_interfaces(num_interface) == 1 ) then
            if ( i == iproc ) then
              write(16,*) j, my_nb_interfaces(num_interface)
            else
              write(16,*) i, my_nb_interfaces(num_interface)
            end if
            
            count_faces = 0
            do k = tab_size_interfaces(num_interface), tab_size_interfaces(num_interface+1)-1
              !print*,'hello',i,iproc
              !stop
              if ( i == iproc ) then
                  local_elmnt = glob2loc_elmnts(tab_interfaces(k*7+0))+1
              else
                  !print*,'what!',k,size(tab_size_interfaces)
                  !print*,minval(tab_size_interfaces),maxval(tab_size_interfaces),num_interface
                  !print*,size(tab_interfaces),(k*7+1)
                  !print*,tab_interfaces(k*7+1),size(glob2loc_elmnts)
                  !stop
                  local_elmnt = glob2loc_elmnts(tab_interfaces(k*7+1))+1
                  
              end if

!!$                  if ( tab_interfaces(k*7+2) == 1 ) then
!!$                     do l = glob2loc_nodes_npart(tab_interfaces(k*7+3)), &
!!$                          glob2loc_nodes_npart(tab_interfaces(k*7+3)+1)-1
!!$                        if ( glob2loc_nodes_parts(l) == iproc ) then
!!$                           local_nodes(1) = glob2loc_nodes(l)+1
!!$                        end if
!!$                     end do
!!$
!!$                     write(16,*) local_elmnt, tab_interfaces(k*7+2), local_nodes(1), -1
!!$                  else
!!$                     if ( tab_interfaces(k*7+2) == 2 ) then
!!$                        do l = glob2loc_nodes_npart(tab_interfaces(k*7+3)), &
!!$                             glob2loc_nodes_npart(tab_interfaces(k*7+3)+1)-1
!!$                           if ( glob2loc_nodes_parts(l) == iproc ) then
!!$                              local_nodes(1) = glob2loc_nodes(l)+1
!!$                           end if
!!$                        end do
!!$                        do l = glob2loc_nodes_npart(tab_interfaces(k*7+4)), &
!!$                           glob2loc_nodes_npart(tab_interfaces(k*7+4)+1)-1
!!$                           if ( glob2loc_nodes_parts(l) == iproc ) then
!!$                              local_nodes(2) = glob2loc_nodes(l)+1
!!$                           end if
!!$                        end do
!!$                        write(16,*) local_elmnt, tab_interfaces(k*7+2), local_nodes(1), local_nodes(2)
!!$                     else
!!$                        write(16,*) "erreur_write_interface_", tab_interfaces(k*7+2)
!!$                     end if
!!$                  end if
              print*,tab_interfaces(k*7+2)
              select case (tab_interfaces(k*7+2))
              case (1)
                  ! single point element
                  do l = glob2loc_nodes_npart(tab_interfaces(k*7+3)), &
                      glob2loc_nodes_npart(tab_interfaces(k*7+3)+1)-1
                    if ( glob2loc_nodes_parts(l) == iproc ) then
                        local_nodes(1) = glob2loc_nodes(l)+1
                    end if
                  end do
                  write(16,*) local_elmnt, tab_interfaces(k*7+2), local_nodes(1), -1, -1, -1
              case (2)
                  ! edge element
                  do l = glob2loc_nodes_npart(tab_interfaces(k*7+3)), &
                      glob2loc_nodes_npart(tab_interfaces(k*7+3)+1)-1
                    if ( glob2loc_nodes_parts(l) == iproc ) then
                        local_nodes(1) = glob2loc_nodes(l)+1
                    end if
                  end do
                  do l = glob2loc_nodes_npart(tab_interfaces(k*7+4)), &
                      glob2loc_nodes_npart(tab_interfaces(k*7+4)+1)-1
                    if ( glob2loc_nodes_parts(l) == iproc ) then
                        local_nodes(2) = glob2loc_nodes(l)+1
                    end if
                  end do
                  write(16,*) local_elmnt, tab_interfaces(k*7+2), local_nodes(1), local_nodes(2), -1, -1
              case (4)
                  ! face element
                  print*,'hi1'
                  count_faces = count_faces + 1
                  do l = glob2loc_nodes_npart(tab_interfaces(k*7+3)), &
                      glob2loc_nodes_npart(tab_interfaces(k*7+3)+1)-1
                    if ( glob2loc_nodes_parts(l) == iproc ) then
                        local_nodes(1) = glob2loc_nodes(l)+1
                    end if
                  end do
                  do l = glob2loc_nodes_npart(tab_interfaces(k*7+4)), &
                      glob2loc_nodes_npart(tab_interfaces(k*7+4)+1)-1
                    if ( glob2loc_nodes_parts(l) == iproc ) then
                        local_nodes(2) = glob2loc_nodes(l)+1
                    end if
                  end do
                  do l = glob2loc_nodes_npart(tab_interfaces(k*7+5)), &
                      glob2loc_nodes_npart(tab_interfaces(k*7+5)+1)-1
                    if ( glob2loc_nodes_parts(l) == iproc ) then
                        local_nodes(3) = glob2loc_nodes(l)+1
                    end if
                  end do
                  do l = glob2loc_nodes_npart(tab_interfaces(k*7+6)), &
                      glob2loc_nodes_npart(tab_interfaces(k*7+6)+1)-1
                    if ( glob2loc_nodes_parts(l) == iproc ) then
                        local_nodes(4) = glob2loc_nodes(l)+1
                    end if
                  end do
                  write(16,*) local_elmnt, tab_interfaces(k*7+2), &
                      local_nodes(1), local_nodes(2),local_nodes(3), local_nodes(4)
              case default
                  print *, "error in write_interfaces_database!", tab_interfaces(k*7+2), iproc
              end select
            end do
      
            ! outputs infos
            !print*,'  partition MPI interface:',iproc,num_interface
            !print*,'    element faces: ',count_faces

        end if

        num_interface = num_interface + 1
      end do
  end do

!end if
close(16)
end subroutine write_interfaces_database
!=======================================================

!--------------------------------------------------
! Write elements on surface boundaries (and their four nodes on boundaries) 
! pertaining to iproc partition in the corresponding Database
!--------------------------------------------------
subroutine write_moho_surface_database(IIN_database, iproc, nelmnts, & 
                      glob2loc_elmnts, glob2loc_nodes_npart, &
                      glob2loc_nodes_parts, glob2loc_nodes, part, &
                      nspec2D_moho,ibelm_moho,nodes_ibelm_moho)
    
integer,intent(in) :: IIN_database
integer,intent(in) :: iproc
integer,intent(in) :: nelmnts !integer(long)

integer,dimension(:),pointer :: glob2loc_elmnts
integer,dimension(:),pointer :: glob2loc_nodes_npart
integer,dimension(:),pointer :: glob2loc_nodes_parts
integer,dimension(:),pointer :: glob2loc_nodes
integer,dimension(1:nelmnts) :: part

integer,intent(in) :: nspec2D_moho
integer,dimension(nspec2D_moho),intent(in) :: ibelm_moho
integer,dimension(4,nspec2D_moho),intent(in) :: nodes_ibelm_moho

integer :: i,j
integer :: loc_node1,loc_node2,loc_node3,loc_node4
integer :: loc_nspec2D_moho
  
! counts number of elements for moho surface in this partition
! optional moho
loc_nspec2D_moho = 0
do i=1,nspec2D_moho
    if(part(ibelm_moho(i)) == iproc) then
      loc_nspec2D_moho = loc_nspec2D_moho + 1
    end if
end do
! checks if anything to do
if( loc_nspec2D_moho == 0 ) return

! format: #surface_id, #number of elements
write(IIN_database,*) 7, loc_nspec2D_moho

! outputs element index and element node indices
! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
!          (this is assigned by CUBIT, if this changes the following indexing must be changed as well)
!          while glob2loc_elmnts(.) is shifted from 0 to nspec-1  thus 
!          we need to have the arg of glob2loc_elmnts start at 0 ==> glob2loc_nodes(ibelm_** -1 )

! optional moho
do i=1,nspec2D_moho    
    if(part(ibelm_moho(i)) == iproc) then
      do j = glob2loc_nodes_npart(nodes_ibelm_moho(1,i)-1), glob2loc_nodes_npart(nodes_ibelm_moho(1,i))-1
          if (glob2loc_nodes_parts(j) == iproc ) then
            loc_node1 = glob2loc_nodes(j)+1
          end if
      end do
      do j = glob2loc_nodes_npart(nodes_ibelm_moho(2,i)-1), glob2loc_nodes_npart(nodes_ibelm_moho(2,i))-1
          if (glob2loc_nodes_parts(j) == iproc ) then
            loc_node2 = glob2loc_nodes(j)+1
          end if
      end do
      do j = glob2loc_nodes_npart(nodes_ibelm_moho(3,i)-1), glob2loc_nodes_npart(nodes_ibelm_moho(3,i))-1
          if (glob2loc_nodes_parts(j) == iproc ) then
            loc_node3 = glob2loc_nodes(j)+1
          end if
      end do
      do j = glob2loc_nodes_npart(nodes_ibelm_moho(4,i)-1), glob2loc_nodes_npart(nodes_ibelm_moho(4,i))-1
          if (glob2loc_nodes_parts(j) == iproc ) then
            loc_node4 = glob2loc_nodes(j)+1
          end if
      end do
      write(IIN_database,*) glob2loc_elmnts(ibelm_moho(i)-1)+1, loc_node1, loc_node2, loc_node3, loc_node4  
    end if

end do
end subroutine write_moho_surface_database
!=======================================================


!--------------------------------------------------
! loading : sets weights for acoustic/elastic elements to account for different 
!               expensive calculations in specfem simulations
!--------------------------------------------------

subroutine acoustic_elastic_load (elmnts_load,nelmnts,nb_materials,num_material,mat_prop)
!
! note: 
!   acoustic material = domainID 1  (stored in mat_prop(6,..) )
!   elastic material    = domainID 2
!
implicit none

integer,intent(in) :: nelmnts !integer(long)
integer,intent(in) :: nb_materials

! load weights
integer,dimension(1:nelmnts),intent(out) :: elmnts_load

! materials  
integer,dimension(1:nelmnts),intent(in) :: num_material
real(kind=kreal),dimension(6,nb_materials),intent(in) :: mat_prop !double precision

! local parameters
logical,dimension(nb_materials) :: is_acoustic,is_elastic    
integer :: i,el

! sets acoustic/elastic flags for materials
is_acoustic(:) = .false.
is_elastic(:) = .false.
do i = 1, nb_materials
    ! acoustic material has idomain_id 1
    if (int(mat_prop(6,i)) == 1 ) then
      is_acoustic(i) = .true.
    endif
    ! elastic material has idomain_id 2
    if (int(mat_prop(6,i)) == 2 ) then
      is_elastic(i) = .true.
    endif
enddo

! sets weights for elements
do el = 0, nelmnts-1
  ! acoustic element (cheap)
  if ( is_acoustic(num_material(el+1)) ) then
    elmnts_load(el+1) = elmnts_load(el+1)*ACOUSTIC_LOAD
  endif
  ! elastic element (expensive)
  if ( is_elastic(num_material(el+1)) ) then
    elmnts_load(el+1) = elmnts_load(el+1)*ELASTIC_LOAD
  endif
enddo
end subroutine acoustic_elastic_load
!=======================================================

!--------------------------------------------------
! Repartitioning : two coupled acoustic/elastic elements are transfered to the same partition
!--------------------------------------------------

subroutine acoustic_elastic_repartitioning (nelmnts, nnodes, elmnts, &
                      nb_materials, num_material, mat_prop, &
                      sup_neighbour, nsize, &
                      nproc, part, nfaces_coupled, faces_coupled)

implicit none

integer,intent(in) :: nelmnts !integer(long)
integer,intent(in) :: nnodes,nproc,nb_materials
integer,intent(in) :: sup_neighbour,nsize !integer(long)

integer,dimension(1:nelmnts),intent(in) :: num_material

real(kind=kreal),dimension(6,nb_materials),intent(in) :: mat_prop !double precision

integer,dimension(0:nelmnts-1) :: part
integer,dimension(0:esize*nelmnts-1) :: elmnts

integer,intent(out) :: nfaces_coupled
integer,dimension(:,:),pointer :: faces_coupled


logical,dimension(nb_materials) :: is_acoustic,is_elastic

! neighbors
integer,dimension(:),allocatable :: xadj
integer,dimension(:),allocatable :: adjncy
integer,dimension(:),allocatable :: nnodes_elmnts
integer,dimension(:),allocatable :: nodes_elmnts
integer :: max_neighbour        

integer :: i,iface
integer :: el,el_adj
logical :: is_repartitioned

! sets acoustic/elastic flags for materials
is_acoustic(:) = .false.
is_elastic(:) = .false.
do i = 1, nb_materials
    if (int(mat_prop(6,i)) == 1 ) then
      is_acoustic(i) = .true.
    endif
    if (int(mat_prop(6,i)) == 2 ) then
      is_elastic(i) = .true.
    endif
enddo

! gets neighbors by 4 common nodes (face)
allocate(xadj(0:nelmnts))
allocate(adjncy(0:sup_neighbour*nelmnts-1))
allocate(nnodes_elmnts(0:nnodes-1))
allocate(nodes_elmnts(0:nsize*nnodes-1))
!call mesh2dual_ncommonnodes(nelmnts, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts,4)
call mesh2dual_ncommonnodes(nelmnts, nnodes, nsize, sup_neighbour, &
                            elmnts, xadj, adjncy, nnodes_elmnts, &
                            nodes_elmnts, max_neighbour, 4)

! counts coupled elements
nfaces_coupled = 0
do el = 0, nelmnts-1
    if ( is_acoustic(num_material(el+1)) ) then
      do el_adj = xadj(el), xadj(el+1) - 1
          if ( is_elastic(num_material(adjncy(el_adj)+1)) ) then
            nfaces_coupled = nfaces_coupled + 1
          endif
      enddo
    endif
enddo

! coupled elements
allocate(faces_coupled(2,nfaces_coupled))

! stores elements indices
nfaces_coupled = 0
do el = 0, nelmnts-1
    if ( is_acoustic(num_material(el+1)) ) then
      do el_adj = xadj(el), xadj(el+1) - 1
          if ( is_elastic(num_material(adjncy(el_adj)+1)) ) then
            nfaces_coupled = nfaces_coupled + 1
            faces_coupled(1,nfaces_coupled) = el
            faces_coupled(2,nfaces_coupled) = adjncy(el_adj)
          endif
      enddo
    endif
enddo

! puts coupled elements into same partition
do i = 1, nfaces_coupled*nproc
    is_repartitioned = .false.
    do iface = 1, nfaces_coupled
      if ( part(faces_coupled(1,iface)) /= part(faces_coupled(2,iface)) ) then
          if ( part(faces_coupled(1,iface)) < part(faces_coupled(2,iface)) ) then
            part(faces_coupled(2,iface)) = part(faces_coupled(1,iface))
          else
            part(faces_coupled(1,iface)) = part(faces_coupled(2,iface))
          endif
          is_repartitioned = .true.
      endif
    enddo
    if ( .not. is_repartitioned ) then
      exit
    endif
enddo
end subroutine acoustic_elastic_repartitioning
!=======================================================

!--------------------------------------------------
! Repartitioning : two coupled moho surface elements are transfered to the same partition
!--------------------------------------------------
subroutine moho_surface_repartitioning (nelmnts,nnodes,elmnts, &
sup_neighbour,nsize,nproc,part,nspec2D_moho,ibelm_moho, &
nodes_ibelm_moho)

implicit none

! number of (spectral) elements  ( <-> nspec )
integer,intent(in) :: nelmnts !integer(long)

! number of (global) nodes, number or processes
integer,intent(in) :: nnodes,nproc 

! maximum number of neighours and max number of elements-that-contain-the-same-node
integer,intent(in) :: sup_neighbour,nsize !integer(long)
    
! partition index on each element
integer,dimension(0:nelmnts-1) :: part

! mesh element indexing
! ( elmnts(esize,nspec) )
integer,dimension(0:esize*nelmnts-1) :: elmnts

! moho surface
integer,intent(in) :: nspec2D_moho
integer,dimension(nspec2D_moho),intent(in) :: ibelm_moho
integer,dimension(4,nspec2D_moho),intent(in) :: nodes_ibelm_moho

! local parameters
integer :: nfaces_coupled
integer,dimension(:,:),pointer :: faces_coupled

logical,dimension(:),allocatable :: is_moho,node_is_moho

! for neighbors
integer,dimension(:),allocatable :: xadj
integer,dimension(:),allocatable :: adjncy
integer,dimension(:),allocatable :: nnodes_elmnts
integer,dimension(:),allocatable :: nodes_elmnts
integer :: max_neighbour        

integer :: i,j,iface,inode,ispec2D,counter
integer :: el,el_adj
logical :: is_repartitioned

! temporary flag arrays
allocate( is_moho(0:nelmnts-1)) ! element ids start from 0
allocate( node_is_moho(0:nnodes-1) ) ! node ids start from 0
is_moho(:) = .false.
node_is_moho(:) = .false.
    
! sets moho flags for known elements
do ispec2D = 1, nspec2D_moho
  ! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
  el = ibelm_moho(ispec2D) - 1
  is_moho(el) = .true.  
  
  ! sets node flags      
  do j=1,4
    ! note: assumes that node indices in nodes_ibelm_* arrays are in the range from 1 to nodes
    inode = nodes_ibelm_moho(j,ispec2D) - 1
    node_is_moho(inode) = .true.
  enddo
enddo

! checks if element has moho surface 
do el = 0, nelmnts-1
  if( is_moho(el) ) cycle
  
  ! loops over all element corners         
  counter = 0   
  do i=0,esize-1
    ! note: assumes that node indices in elmnts array are in the range from 0 to nodes-1
    inode = elmnts(el*esize+i)
    if( node_is_moho(inode) ) counter = counter + 1  
  enddo
  
  ! sets flag if it has a surface
  if( counter == 4 ) is_moho(el) = .true.
enddo

! statistics output
counter = 0
do el=0, nelmnts-1
  if ( is_moho(el) ) counter = counter + 1
enddo
print*,'  moho elements = ',counter

! gets neighbors by 4 common nodes (face)
allocate(xadj(0:nelmnts)) ! contains number of adjacent elements (neighbours)
allocate(adjncy(0:sup_neighbour*nelmnts-1)) ! contains all element id indices of adjacent elements
allocate(nnodes_elmnts(0:nnodes-1))
allocate(nodes_elmnts(0:nsize*nnodes-1))

call mesh2dual_ncommonnodes(nelmnts, nnodes, nsize, sup_neighbour, &
                    elmnts, xadj, adjncy, nnodes_elmnts, &
                    nodes_elmnts, max_neighbour, 4)

! counts coupled elements
nfaces_coupled = 0
do el = 0, nelmnts-1
    if ( is_moho(el) ) then
      do el_adj = xadj(el), xadj(el+1) - 1
        ! increments counter if it contains face
        if( is_moho(adjncy(el_adj)) ) nfaces_coupled = nfaces_coupled + 1
      enddo
    endif
enddo

! coupled elements
allocate(faces_coupled(2,nfaces_coupled))

! stores elements indices
nfaces_coupled = 0
do el = 0, nelmnts-1
    if ( is_moho(el) ) then
      do el_adj = xadj(el), xadj(el+1) - 1
          if ( is_moho(adjncy(el_adj)) ) then
            nfaces_coupled = nfaces_coupled + 1
            faces_coupled(1,nfaces_coupled) = el
            faces_coupled(2,nfaces_coupled) = adjncy(el_adj)
          endif
      enddo
    endif
enddo

! puts coupled elements into same partition
do i = 1, nfaces_coupled*nproc
    is_repartitioned = .false.
    do iface = 1, nfaces_coupled
      if ( part(faces_coupled(1,iface)) /= part(faces_coupled(2,iface)) ) then
          ! coupled moho elements are in different partitions
          if ( part(faces_coupled(1,iface)) < part(faces_coupled(2,iface)) ) then
            part(faces_coupled(2,iface)) = part(faces_coupled(1,iface))
          else
            part(faces_coupled(1,iface)) = part(faces_coupled(2,iface))
          endif
          is_repartitioned = .true.
      endif
    enddo
    if ( .not. is_repartitioned ) then
      exit
    endif
enddo
end subroutine moho_surface_repartitioning
!=======================================================

function rank_array(x,n) result(irank)
integer,intent(in) :: n ! size of the vector data x
integer,dimension(n) :: x ! data vector to sort
integer :: temp
integer :: i,j,inum
integer,dimension(n) :: irank,xnew

do i = 2, n
j = i - 1
temp = x(i)
do while (j>=1 .and. x(j)>temp)
  x(j+1) = x(j)
  j = j - 1
end do
x(j+1) = temp
end do
xnew=x
inum=1
irank(1)=1
do i=2,n
if(x(i)/=x(i-1))inum=inum+1
irank(i)=inum
enddo

end function rank_array
!=======================================================

function sort(x,n) result(xnew)
integer,intent(in) :: n ! size of the vector data x
integer,dimension(n) :: x ! data vector to sort
integer :: temp
integer :: i,j
integer,dimension(n) :: xnew

do i = 2, n
j = i - 1
temp = x(i)
do while (j>=1 .and. x(j)>temp)
  x(j+1) = x(j)
  j = j - 1
end do
x(j+1) = temp
end do
xnew=x
end function sort
!=======================================================

subroutine insertion_sort(x,n)
integer,intent(in) :: n ! size of the vector data x
real,intent(inout),dimension(n) :: x ! data vector to sort
real :: temp
integer :: i,j

do i = 2, n
j = i - 1
temp = x(i)
do while (j>=1 .and. x(j)>temp)
  x(j+1) = x(j)
  j = j - 1
end do
x(j+1) = temp
end do
end subroutine insertion_sort
!=======================================================

function is_equal(x1,x2,n) result(flag)
integer,intent(in) :: n ! size of the vectors
integer,dimension(n),intent(in) :: x1,x2 ! data vectors to compare
logical :: flag
integer :: i
flag=.true. 
do i = 1,n
if(x1(i)/=x2(i))then
  flag=.false.
  return
endif
enddo
return
end function is_equal
!=======================================================
end module partmesh_library

