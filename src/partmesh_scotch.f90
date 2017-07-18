! this module contains the routines to partition mesh.
! this library is copied and modified from the original
! SPECFEM3D package (Komatitsch and Tromp 1999, Peter et al. 2011)
! REVISION
!   HNG, Dec 09,2011; HNG, Jul 12,2011; HNG, Apr 09,2010
module partmesh_scotch
use global
use partmesh_library
use string_library,only:str2int,str2real
implicit none

!include './scotchf.h'
include "ptscotchf.h"

! number of partitions
integer :: npart ! e.g. 4 for partitioning for 4 CPUs or 4 processes

! mesh arrays
integer,dimension(:),allocatable :: part

integer,dimension(:),allocatable :: xadj
integer,dimension(:),allocatable :: adjncy
integer,dimension(:),allocatable :: nnodeelmt
integer,dimension(:),allocatable :: nodeselmt
integer,dimension(:),allocatable :: elmt_weight

integer,dimension(:),allocatable :: glob2loc_elmt
integer,dimension(:),allocatable :: glob2loc_nodes_npart
integer,dimension(:),allocatable :: glob2loc_nodes_parts
integer,dimension(:),allocatable :: glob2loc_nodes

integer,dimension(:),allocatable :: tab_size_interfaces, tab_interfaces
integer,dimension(:),allocatable :: my_interfaces
integer,dimension(:),allocatable :: my_nb_interfaces
integer ::  ninterfaces
integer :: my_ninterface

! max number of elements that contain the same node.
integer :: nsize !integer(long).
integer :: nb_edges

integer :: max_neighbour ! Real maximum number of neighbours per element
! integer(long). Majoration of the maximum number of neighbours per element
integer :: sup_neighbour

integer :: nnode_loc, nelmt_loc
integer :: num_elmnt, num_node, num_mat

! viscoelastic material
integer :: nvisco
integer,allocatable :: imat_visco(:)
real(kind=kreal),allocatable :: muratio(:),viscosity(:)
! water table information
integer :: nwmat
integer,allocatable :: waterid(:)

!! boundaries
!integer :: nelmt2D_xmin,nelmt2D_xmax,nelmt2D_ymin,nelmt2D_ymax,nelmt2D_bottom, &
!nelmt2D_top
!integer,dimension(:),allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,          &
!ibelm_ymax,ibelm_bottom,ibelm_top
!integer,dimension(:),allocatable :: nodes_ibelm_xmin,nodes_ibelm_xmax,         &
!nodes_ibelm_ymin
!integer,dimension(:,:),allocatable :: nodes_ibelm_ymax,nodes_ibelm_bottom,     &
!nodes_ibelm_top
!
!! moho surface (optional)
!integer :: nelmt2D_moho
!integer,dimension(:),allocatable :: ibelm_moho
!integer,dimension(:,:),allocatable :: nodes_ibelm_moho

character(len=60) :: lineword(9)
character(len=256) :: prname

!logical,dimension(:),allocatable :: mask_nodeselmt
!integer,dimension(:),allocatable :: used_nodeselmt

real(kind=kreal),dimension(SCOTCH_GRAPHDIM) :: scotchgraph !double precision
real(kind=kreal),dimension(SCOTCH_STRATDIM) :: scotchstrat !double precision
!character(len=256),parameter :: scotch_strategy='b{job=t,map=t,poli=S,sep=h{pass=30}}'

!pll
character(len=60),allocatable :: mfile(:)
real(kind=kreal),dimension(:,:),allocatable :: mat_prop ! double precision
!integer,dimension(:),allocatable :: mat_domain
integer,dimension(:),allocatable :: undef_mat_domain
integer :: count_def_mat,count_undef_mat,imat
character (len=30),dimension(:,:),allocatable :: undef_mat_prop

! default mesh file directory
character(len=256) :: out_phead ! output path and header

contains

! reads in mesh files
subroutine read_mesh_files
implicit none
integer :: i,ios,istat

character(len=80) :: inp_fname

print*,'matfile:',trim(adjustl(matfile))
! sets number of nodes per element
ngnode = esize

!! reads node coordinates
!! open file to read
!inp_fname=trim(inp_path)//trim(xfile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'could not open file:',trim(inp_fname)
!  stop
!endif
!read(98,*)nnode
!allocate(g_coord(3,nnode))
!read(98,*)g_coord(1,:)
!!write(*,'(f25.15)')g_coord(1,3)
!!stop
!close(98)
!
!! open file to read
!inp_fname=trim(inp_path)//trim(yfile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'could not open file:',trim(inp_fname)
!  stop
!endif
!read(98,*) nnode
!read(98,*)g_coord(2,:)
!close(98)
!
!! open file to read
!inp_fname=trim(inp_path)//trim(zfile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'could not open file:',trim(inp_fname)
!  stop
!endif
!read(98,*) nnode
!read(98,*)g_coord(3,:)
!close(98)
!
!print*, 'total nodes:', nnode
!
!! reads mesh elements indexing
!!(CUBIT calls this the connectivity, guess in the sense that it connects with
!! the points index in the global coordinate file "g_coord_file"; it doesn't
!! tell you which point is connected with others)
!! open file to read
!inp_fname=trim(inp_path)//trim(confile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'could not open file:',trim(inp_fname)
!  stop
!endif
!
!read(98,*) nelmt
!allocate(g_num(esize,nelmt))
!do ispec = 1, nelmt
!  ! format: # element_id  #id_node1 ... #id_node8
!
!  ! note: be aware that here we can have different node ordering for a cube element;
!  !          the ordering from Cubit files might not be consistent for multiple volumes, or uneven, unstructured grids
!  !
!  !          guess here it assumes that spectral elements ordering is like first at the bottom of the element, anticlock-wise, i.e.
!  !             point 1 = (0,0,0), point 2 = (0,1,0), point 3 = (1,1,0), point 4 = (1,0,0)
!  !          then top (positive z-direction) of element
!  !             point 5 = (0,0,1), point 6 = (0,1,1), point 7 = (1,1,1), point 8 = (1,0,1)
!
!  !read(98,*) num_elmnt, g_num(5,num_elmnt), g_num(1,num_elmnt),g_num(4,num_elmnt), g_num(8,num_elmnt), &
!  !      g_num(6,num_elmnt), g_num(2,num_elmnt), g_num(3,num_elmnt), g_num(7,num_elmnt)
!
!  read(98,*)g_num(:,ispec)
!
!  !outputs info for each element to see ordering
!  !print*,'ispec: ',ispec
!  !print*,'  ',num_elmnt, g_num(5,num_elmnt), g_num(1,num_elmnt),g_num(4,num_elmnt), g_num(8,num_elmnt), &
!  !      g_num(6,num_elmnt), g_num(2,num_elmnt), g_num(3,num_elmnt), g_num(7,num_elmnt)
!  !print*,'elem:',num_elmnt
!  !do i=1,8
!  !  print*,' i ',i,'val :',g_num(i,num_elmnt),&
!  !    g_coord(1,g_num(i,num_elmnt)),g_coord(2,g_num(i,num_elmnt)),g_coord(3,g_num(i,num_elmnt))
!  !enddo
!  !print*
!
!enddo
!close(98)
!print*, 'total elements:', nelmt
!
!! reads material associations
!! open file to read
!inp_fname=trim(inp_path)//trim(idfile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'could not open file:',trim(inp_fname)
!  stop
!endif
!
!read(98,*) nelmt
!!print*,nelmt; stop
!allocate(mat_id(nelmt))
!
!do i_spec = 1, nelmt
!  ! format: # id_element #flag
!  ! note: we assume elements are sorted in materials_file
!  read(98,*)ispec,mat_id(ispec);
!enddo
!close(98)

! reads material definitions
!
! note: format of nummaterial_velocity_file must be
!
! #(1)material_domain_id #(2)material_id  #(3)gam  #(4)ym   #(5)nu   #(6)phi  #(7)anisotropy_flag
!
! where
!     material_domain_id : 1=elastic / 2=acoustic / 3=poroelastic
!     material_id               : number of material/volume
!     gam                           : density
!     ym                             : P-velocity
!     nu                             : S-velocity
!     phi                      : 0=no attenuation/1=IATTENUATION_SEDIMENTS_40, 2=..., 13=IATTENUATION_BEDROCK
!     anisotropy_flag        : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90
!count_def_mat = 0
count_undef_mat = 0
! open file to read
inp_fname=trim(inp_path)//trim(adjustl(matfile))
open(unit=98, file=trim(inp_fname),&
  status='old', form='formatted', iostat = istat)
if( istat /= 0 ) then
  print*,'could not open file:',trim(inp_fname)
  stop
endif

! note: format #material_domain_id #material_id #...
!read(98,*) ! skip line
!read(98,*,iostat=istat) idummy,num_mat
!print *,'materials:'
! counts materials (defined/undefined)
!do while (istat == 0)
!   print*, '  num_mat = ',num_mat
!   if(num_mat /= -1) then
!      count_def_mat = count_def_mat + 1
!   else
!      count_undef_mat = count_undef_mat + 1
!   end if
!   read(98,*,iostat=istat) idummy,num_mat
!enddo
!close(98)
read(98,*) ! skip 1 line
read(98,*)count_def_mat
print*, '  defined = ',count_def_mat, 'undefined = ',count_undef_mat
! check with material flags
if( count_def_mat > 0 .and. maxval(mat_id) > count_def_mat ) then
  print*,'error material definitions:'
  print*,'  materials associated in materials_file:',maxval(mat_id)
  print*,'  bigger than defined materials in nummaterial_velocity_file:',count_def_mat
  stop 'error materials'
endif
allocate(mfile(count_def_mat))
mfile=""
allocate(mat_prop(6,count_def_mat),mat_domain(count_def_mat),type_blk(count_def_mat))
allocate(undef_mat_prop(6,count_undef_mat),undef_mat_domain(count_def_mat))
! reads in defined material properties
!open(unit=98, file=inp_path(1:len_trim(inp_path))//trim(matfile), &
!      status='old', form='formatted')
!read(98,*) ! skip ONE line
do imat=1,count_def_mat
    ! material definitions
    !
    ! format: note that we save the arguments in a slightly different order in
    ! mat_prop(:,:)
    ! #(6) material_domain_id #(0) material_id  #(1) gam #(2) ym #(3) nu #(4)
    ! phi #(5) anisotropy_flag
    !
    ! idomain_id,gam,ym,nu,phi,coh,psi
    !read(98,*) num_mat,mat_domain(num_mat),mat_prop(:,num_mat)
    !read(98,*) num_mat, mat_prop(1,num_mat),mat_prop(2,num_mat),&
    !           mat_prop(3,num_mat),mat_prop(4,num_mat),mat_prop(5,num_mat)
    !mat_domian(num_mat) = idomain_id
    !mat_prop(1,num_mat) = gam
    !mat_prop(2,num_mat) = ym
    !mat_prop(3,num_mat) = nu
    !mat_prop(4,num_mat) = phi
    !mat_prop(5,num_mat) = coh
    !mat_prop(6,num_mat) = psi
    lineword=""
    read(98,*,iostat=ios)lineword
    num_mat=str2int(lineword(1))
    mat_domain(num_mat)=str2int(lineword(2))
    type_blk(num_mat)=str2int(lineword(3))
    if(type_blk(num_mat).eq.0)then
      ! block material properties
      mat_prop(1,num_mat)=str2real(lineword(4))
      mat_prop(2,num_mat)=str2real(lineword(5))
      mat_prop(3,num_mat)=str2real(lineword(6))
      mat_prop(4,num_mat)=str2real(lineword(7))
      mat_prop(5,num_mat)=str2real(lineword(8))
      mat_prop(6,num_mat)=str2real(lineword(9))
    elseif(type_blk(num_mat).eq.-1)then
      ! tomographic model defined on regular structured grid
      mfile(num_mat)=trim(lineword(4))
    else
      print*,'ERROR: type_blk:',type_blk(num_mat),' is unsupported!'
      stop
    endif

    if(num_mat < 0 .or. num_mat > count_def_mat)then
      print*,"ERROR: Invalid nummaterial_velocity_file file."
      stop
    endif

enddo
! reads in undefined material properties
do imat=1,count_undef_mat
  read(98,'(i6,i6,6a30)')num_mat,undef_mat_domain(num_mat),                    &
  undef_mat_prop(:,num_mat)
  !undef_mat_prop(7,imat),undef_mat_prop(6,imat),undef_mat_prop(1,imat),       &
  !undef_mat_prop(2,imat),undef_mat_prop(3,imat),undef_mat_prop(4,imat),       &
  !undef_mat_prop(5,imat)
enddo
! viscoelastic properties
nvisco=count(mat_domain.eq.VISCOELASTIC_DOMAIN)
allocate(imat_visco(nvisco),muratio(nvisco),viscosity(nvisco))
do i=1,nvisco
  read(98,*)imat_visco(i),muratio(i),viscosity(i)
enddo

! water properties
nwmat=0
read(98,*,iostat=ios)nwmat
allocate(waterid(nwmat))
do i=1,nwmat
  read(98,*)waterid(i)
enddo
close(98)

!! reads in absorbing boundary files
!! open file to read
!inp_fname=trim(inp_path)//trim(uxfile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'WARNING: could not open file:',trim(inp_fname)
!  !stop
!endif
!
!if( istat /= 0 ) then
!  nelmt2D_xmin = 0
!else
!  read(98,*) nelmt2D_xmin
!endif
!allocate(ibelm_xmin(nelmt2D_xmin))
!allocate(nodes_ibelm_xmin(nelmt2D_xmin))
!do ispec2D = 1,nelmt2D_xmin
!  ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
!  ! note: ordering for CUBIT seems such that the normal of the face points outward of the element the face belongs to;
!  !         in other words, nodes are in increasing order such that when looking from within the element outwards,
!  !         they are ordered clockwise
!  !
!  !          doesn't necessarily have to start on top-rear, then bottom-rear, bottom-front, and finally top-front i.e.:
!  !          point 1 = (0,1,1), point 2 = (0,1,0), point 3 = (0,0,0), point 4 = (0,0,1)
!  read(98,*) ibelm_xmin(ispec2D), nodes_ibelm_xmin(ispec2D)
!
!  !outputs info for each element for check of ordering
!  !print*,'ispec2d:',ispec2d
!  !print*,'  xmin:', ibelm_xmin(ispec2D), nodes_ibelm_xmin(1,ispec2D), nodes_ibelm_xmin(2,ispec2D), &
!  !      nodes_ibelm_xmin(3,ispec2D), nodes_ibelm_xmin(4,ispec2D)
!  !do i=1,4
!  !  print*,'i',i,'val:',ibelm_xmin(ispec2d),g_coord(1,nodes_ibelm_xmin(i,ispec2D)), &
!  !      g_coord(2,nodes_ibelm_xmin(i,ispec2D)),g_coord(3,nodes_ibelm_xmin(i,ispec2D))
!  !enddo
!  !print*
!enddo
!close(98)
!print*, 'absorbing boundaries:'
!print*, '  nelmt2D_xmin = ', nelmt2D_xmin
!
!! reads in absorbing boundary files
!inp_fname=trim(inp_path)//trim(uyfile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'WARNING: could not open file:',trim(inp_fname)
!  !stop
!endif
!
!if( istat /= 0 ) then
!  nelmt2D_xmax = 0
!else
!  read(98,*) nelmt2D_xmax
!endif
!allocate(ibelm_xmax(nelmt2D_xmax))
!allocate(nodes_ibelm_xmax(nelmt2D_xmax))
!do ispec2D = 1,nelmt2D_xmax
!  ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
!  read(98,*) ibelm_xmax(ispec2D), nodes_ibelm_xmax(ispec2D)
!enddo
!close(98)
!print*, '  nelmt2D_xmax = ', nelmt2D_xmax
!
!! reads in absorbing boundary files
!inp_fname=trim(inp_path)//trim(uzfile)
!open(unit=98, file=trim(inp_fname),&
!  status='old', form='formatted', iostat = istat)
!if( istat /= 0 ) then
!  print*,'WARNING: could not open file:',trim(inp_fname)
!  !stop
!endif
!
!if( istat /= 0 ) then
!  nelmt2D_ymin = 0
!else
!  read(98,*) nelmt2D_ymin
!endif
!allocate(ibelm_ymin(nelmt2D_ymin))
!allocate(nodes_ibelm_ymin(nelmt2D_ymin))
!do ispec2D = 1,nelmt2D_ymin
!  ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
!  read(98,*) ibelm_ymin(ispec2D), nodes_ibelm_ymin(ispec2D)
!enddo
!close(98)
!print*, '  nelmt2D_ymin = ', nelmt2D_ymin
end subroutine read_mesh_files
!=======================================================

! checks valence of nodes
subroutine check_valence
implicit none
integer :: ispec,inode
logical,dimension(:),allocatable :: mask_nodeselmt
integer,dimension(:),allocatable :: used_nodeselmt

allocate(mask_nodeselmt(nnode))
allocate(used_nodeselmt(nnode))
mask_nodeselmt(:) = .false.
used_nodeselmt(:) = 0
do ispec = 1, nelmt
  do inode = 1, ESIZE
    mask_nodeselmt(g_num(inode,ispec)) = .true.
    used_nodeselmt(g_num(inode,ispec)) = used_nodeselmt(g_num(inode,ispec)) + 1
  enddo
enddo
print *, 'nodes valence: '
print *, '  min = ',minval(used_nodeselmt(:)),'max = ', maxval(used_nodeselmt(:))
do inode = 1, nnode
  if (.not. mask_nodeselmt(inode)) then
    stop 'ERROR : nodes not used.'
  endif
enddo
nsize = maxval(used_nodeselmt(:)) ! max number of element sharing a node
sup_neighbour = ngnode * nsize - (ngnode + (ngnode/2 - 1)*nfaces)
print*, '  nsize = ',nsize, 'sup_neighbour = ', sup_neighbour

deallocate(mask_nodeselmt,used_nodeselmt)
end subroutine check_valence
!=======================================================

! divides model into partitions using scotch library functions
subroutine scotch_partitioning
implicit none
integer :: istat

g_num = g_num - 1

! determines maximum neighbors based on 1 common node
allocate(xadj(1:nelmt+1))
allocate(adjncy(1:sup_neighbour*nelmt))
allocate(nnodeelmt(1:nnode))
allocate(nodeselmt(1:nsize*nnode))
call mesh2dual_ncommonnode(nelmt,nnode,nsize,sup_neighbour,g_num,xadj,      &
adjncy,nnodeelmt,nodeselmt,max_neighbour,1)
print*, 'mesh2dual: '
print*, '  max_neighbour = ',max_neighbour
print*, '  sup_neighbour = ', sup_neighbour
!print*,xadj
!print*,adjncy
nb_edges = xadj(nelmt+1)

! allocates & initializes partioning of elements
allocate(part(1:nelmt))
part(:) = -1

! initializes
! elements load array
allocate(elmt_weight(1:nelmt))

! uniform load
elmt_weight(:) = 1

! in case of acoustic/elastic simulation, weights elements accordingly
call acoustic_elastic_load(elmt_weight,nelmt,count_def_mat,mat_id,mat_domain)

! SCOTCH partitioning
call scotchfstratinit (scotchstrat(1), istat)
  if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot initialize strat'
endif

! no need to use this for default strategy
!call scotchfstratgraphmap (scotchstrat(1), trim(scotch_strategy), istat)
!  if (istat /= 0) then
!    stop 'ERROR : MAIN : Cannot build strat'
!endif

call scotchfgraphinit (scotchgraph(1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot initialize graph'
endif

! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
! arguments: #(1) graph_structure       #(2) baseval(either 0/1)    #(3) number_of_vertices
!                    #(4) adjacency_index_array         #(5) adjacency_end_index_array (optional)
!                    #(6) vertex_load_array (optional) #(7) vertex_label_array
!                    #(7) number_of_arcs                    #(8) adjacency_array
!                    #(9) arc_load_array (optional)      #(10) istator
call scotchfgraphbuild (scotchgraph(1), 0, nelmt, &
                      xadj(1), xadj(1), &
                      elmt_weight (1), xadj (1), &
                      nb_edges, adjncy(1), &
                      adjncy(1), istat)

! w/out element load, but adjacency array
!call scotchfgraphbuild (scotchgraph (1), 0, nelmt, &
!                      xadj (1), xadj (1), &
!                      xadj (1), xadj (1), &
!                      nb_edges, adjncy (1), &
!                      adjncy (1), istat)


if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot build graph'
endif

call scotchfgraphcheck (scotchgraph (1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Invalid check'
endif

call scotchfgraphpart (scotchgraph (1), npart, scotchstrat(1),part(1),istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot part graph'
endif

call scotchfgraphexit (scotchgraph (1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot destroy graph'
endif

call scotchfstratexit (scotchstrat(1), istat)
if (istat /= 0) then
    stop 'ERROR : MAIN : Cannot destroy strat'
endif


! re-partitioning puts acoustic-elastic coupled elements into same partition
!  integer  :: nfaces_coupled
!  integer, dimension(:,:), allocatable  :: faces_coupled
!    call acoustic_elastic_repartitioning (nelmt, nnode, g_num, &
!                   count_def_mat, mat_id , mat_prop, &
!                   sup_neighbour, nsize, &
!                   npart, part, nfaces_coupled, faces_coupled)


! local number of each element for each partition
call Construct_glob2loc_elmt(nelmt, part, glob2loc_elmt,npart)

! local number of each node for each partition
call Construct_glob2loc_nodes(nelmt,nnode,nsize,nnodeelmt,nodeselmt,   &
part,glob2loc_nodes_npart,glob2loc_nodes_parts,glob2loc_nodes,npart)

! mpi interfaces

! call detect_ghost(out_phead,nelmt,nnode,part+1,g_num+1,npart,              &
! glob2loc_elmt(0:nelmt-1)+1) ! I need all indices starting from 1 not 0
! acoustic/elastic boundaries WILL BE SEPARATED into different MPI partitions
! TODO:WARNING:call below may not be necessary
call Construct_interfaces(nelmt,sup_neighbour,part,g_num,xadj,adjncy,         &
tab_interfaces,tab_size_interfaces,ninterfaces,npart)
!or: uncomment if you want acoustic/elastic boundaries NOT to be separated into
!different MPI partitions
!call Construct_interfaces_no_ac_el_sep(nelmt,sup_neighbour,part,g_num,xadj,  &
!adjncy,tab_interfaces,tab_size_interfaces,ninterfaces,count_def_mat,          &
!mat_prop(3,:),mat_id,npart)
end subroutine scotch_partitioning
!=======================================================

! writes out new Databases files for each partition
subroutine write_mesh_databases
implicit none
integer :: ipart
!character(len=20) :: format_str
!character(len=80) :: out_fname

allocate(my_interfaces(0:ninterfaces-1))
allocate(my_nb_interfaces(0:ninterfaces-1))

!write(format_str,*)ceiling(log10(real(npart)+1))
!format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//',a)'
!print*,format_str; stop

! writes out Database file for each partition
do ipart = 0, npart-1

  ! opens output file
  !write(out_fname, fmt=format_str)trim(out_phead)//'proc',ipart,'_database'
  !print*,trim(out_fname); stop
  !open(unit=15,file=trim(out_fname),&
  !     status='unknown', action='write', form='formatted', iostat = istat)
  !if( istat /= 0 ) then
  ! print*,'error file open:',trim(out_fname)
  ! stop
  !endif

  ! gets number of nodes
  call write_glob2loc_nodes_database(out_path,coordfile,ipart,         &
  nnode_loc,g_coord,glob2loc_nodes_npart,glob2loc_nodes_parts,           &
  glob2loc_nodes,nnode,npart)

  ! gets number of spectral elements
  call write_partition_database(out_path,confile,idfile,ipart,nelmt_loc,nelmt, &
  g_num,glob2loc_elmt,glob2loc_nodes_npart,glob2loc_nodes_parts,            &
  glob2loc_nodes,part,mat_id,ngnode,npart)

  ! writes out node coordinate locations
  !write(15,*) nnode_loc

  !call write_glob2loc_nodes_database(15, ipart, nnode_loc, g_coord,&
  !  glob2loc_nodes_npart, glob2loc_nodes_parts, &
  !  glob2loc_nodes, nnode, 2)

  call write_material_properties_database(out_path,matfile,count_def_mat,      &
  mat_domain,type_blk,mfile,mat_prop,nvisco,  &
  imat_visco,muratio,viscosity,nwmat,waterid,ipart,npart)

  ! writes out spectral element indices
  !write(15,*) nelmt_loc

  !call write_partition_database(15, ipart, nelmt_loc, nelmt, g_num, &
  !  glob2loc_elmt, glob2loc_nodes_npart, &
  !  glob2loc_nodes_parts, glob2loc_nodes, part, mat_id, ngnode, 2)

  ! writes out absorbing/free-surface boundaries
  !call write_boundaries_database(out_phead,ipart, nelmt, nelmt2D_xmin, nelmt2D_xmax, nelmt2D_ymin, &
  !  nelmt2D_ymax, nelmt2D_bottom, nelmt2D_top, &
  !  ibelm_xmin, ibelm_xmax, ibelm_ymin, &
  !  ibelm_ymax, ibelm_bottom, ibelm_top, &
  !  nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
  !  nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
  !  glob2loc_elmt,part,npart)

  ! write MPI interfaces
  !call Write_interfaces_database(out_phead,tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
  !  my_ninterface, my_interfaces, my_nb_interfaces, &
  !  glob2loc_elmt, glob2loc_nodes_npart, glob2loc_nodes_parts, &
  !  glob2loc_nodes, npart)
enddo

deallocate(mat_id,g_coord,xadj,adjncy,nnodeelmt,nodeselmt,elmt_weight)
deallocate(tab_size_interfaces,tab_interfaces,my_interfaces,my_nb_interfaces)

deallocate(mfile,mat_prop,mat_domain)
deallocate(undef_mat_prop,undef_mat_domain)
if(nwmat>0)deallocate(waterid)

! write BC files
write(*,'(a)',advance='no')'writing bc files...'
call write_ssbc(out_path,inp_path,uxfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
call write_ssbc(out_path,inp_path,uyfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
call write_ssbc(out_path,inp_path,uzfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
write(*,'(a)')'complete!'
if(istraction)then
  ! write traction files
  write(*,'(a)',advance='no')'writing traction files...'
  call write_traction(out_path,inp_path,trfile,part+1,npart,glob2loc_elmt(0:nelmt-1)+1)
  write(*,'(a)')'complete!'
endif
deallocate(glob2loc_nodes_npart,glob2loc_nodes_parts,glob2loc_nodes)

!deallocate(ibelm_xmin,nodes_ibelm_xmin)
!deallocate(ibelm_xmax,nodes_ibelm_xmax)
!deallocate(ibelm_ymin,nodes_ibelm_ymin)

write(*,'(a)',advance='yes')'finding interfaces...'
call find_interface(nelmt,nnode,part+1,g_num+1,npart) ! I need all indices starting from 1 not 0
write(*,'(a)')'complete!'
deallocate(g_num,part)
write(*,'(a)',advance='yes')'writing interfaces...'
call detect_ghost(out_phead,nnode,npart,max_neighbour,                   &
glob2loc_elmt(0:nelmt-1)+1) ! I need all indices starting from 1 not 0
write(*,'(a)')'complete!'
deallocate(glob2loc_elmt)

write(*,*)'number of partitions: ',npart
!print*, 'finished successfully'
!write(*,*)'-----------------------------------'
end subroutine write_mesh_databases
!=======================================================

end module partmesh_scotch

