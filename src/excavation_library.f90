! this module contains the routines requred to analyse the excavation topology
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module excavation_library
use set_precision
implicit none
logical,allocatable :: isnode_intact(:)
! number of active ghost partitions for a node
integer,allocatable :: ngpart_node(:)

contains
!-------------------------------------------------------------------------------
! this subroutine analyzes the excavation and determine the elements in the
! intact and void regions
subroutine intact_void_elmt(nexcavid,excavid,ismat,nelmt_intact,nelmt_void,    &
elmt_intact,elmt_void)
use global,only:mat_id,nelmt,nnode,g_num,nmatblk
implicit none
integer,intent(in) :: nexcavid,nelmt_intact,nelmt_void
integer,intent(in) :: excavid(nexcavid)
logical,intent(in) :: ismat(nmatblk)
integer,intent(out) :: elmt_intact(nelmt_intact),elmt_void(nelmt_void)
integer :: i_elmt,ielmt_intact,ielmt_void
logical :: ismat_off(nmatblk)

ismat_off=.false.
ismat_off(excavid)=.true.

! find intact and void elements
ielmt_intact=0; ielmt_void=0
do i_elmt=1,nelmt
  if(ismat_off(mat_id(i_elmt)))then !(mat_id(i_elmt)==excavid)then
    ielmt_void=ielmt_void+1
    elmt_void(ielmt_void)=i_elmt
  elseif(ismat(mat_id(i_elmt)))then
    ielmt_intact=ielmt_intact+1
    elmt_intact(ielmt_intact)=i_elmt
  endif
enddo
if(ielmt_void/=nelmt_void .or. ielmt_intact/=nelmt_intact)then
  write(*,'(/,a)')'ERROR: counted intact/void elements mismatch!'
  stop
endif

! find intact and void nodes
isnode_intact=.false.
do i_elmt=1,nelmt_intact
  isnode_intact(g_num(:,elmt_intact(i_elmt)))=.true.
enddo
return
end subroutine intact_void_elmt
!===============================================================================

! this subroutine analyzes the excavation and determines the nodes
! in the intact and void regions
subroutine intact_void_node(isnode,nnode_intact,nnode_void,node_intact,        &
node_void,nmir)
use global,only:nnode
implicit none
logical,intent(in) :: isnode(nnode)
integer,intent(in) ::nnode_intact,nnode_void
integer,intent(out) :: node_intact(nnode_intact),node_void(nnode_void),        &
nmir(nnode)
integer :: i_node,inode_intact,inode_void

! find intact and void nodes
inode_intact=0; inode_void=0
do i_node=1,nnode
  if(isnode(i_node))then
    inode_intact=inode_intact+1
    node_intact(inode_intact)=i_node
    nmir(i_node)=inode_intact
  else
    inode_void=inode_void+1
    node_void(inode_void)=i_node
  endif
enddo
if(inode_intact/=nnode_intact .and. inode_void/=nnode_void)then
  write(*,'(/,a)')'ERROR: counted intact/void nodes mismatch!'
  stop
endif
return
end subroutine intact_void_node
!===============================================================================

! this subroutine modifies the gdof array setting all the freedoms corresponding
! to the dead nodes as fixed
subroutine modify_gdof(gdof,nnode_void,node_void,neq)
use global,only:nndof,nnode
implicit none
integer,intent(inout) :: gdof(nndof,nnode)
integer,intent(in) :: nnode_void,node_void(nnode_void)
integer,intent(out) :: neq
integer :: i,j

! set global degrees of freedoms of dead nodes to 0
gdof(:,node_void)=0
! compute modified gdof
neq=0
do j=1,nnode
  do i=1,nndof
    if(gdof(i,j)/=0)then
      neq=neq+1
      gdof(i,j)=neq
    endif
  enddo
enddo
return
end subroutine modify_gdof
!===============================================================================

! this subroutine corrects the nodal valency removing the dead-element nodes
subroutine correct_nvalency(node_valency,nelmt_void,gnum_void)
use global,only : ngll,nnode
implicit none
integer,intent(inout) :: node_valency(nnode)
integer,intent(in) :: nelmt_void
!integer,intent(in) :: gnum_void(ngll,nelmt_void)
integer,intent(in) :: gnum_void(:,:)
integer :: i_elmt,num(ngll)

! correct node valency subtracting dead-element nodes
do i_elmt=1,nelmt_void
  !ielmt=elmt_void(i_elmt)
  num=gnum_void(:,i_elmt)
  node_valency(num)=node_valency(num)-1
enddo
return
end subroutine correct_nvalency
!===============================================================================

! this subroutine computes the excavation loads. this load consists of both
! gravity and stress load
subroutine excavation_load(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,dshape_hex8,  &
dlagrange_gll,gll_weights,stress_elmt,extload)
use global,only:ndim,nedof,nenod,ngll,ngnode,nst,g_coord,gam_blk
use preprocess,only:compute_bmat
use math_library,only:determinant,invert
implicit none
integer,intent(in) :: nelmt,neq,gnod(8) ! nelmt (only void elements)
integer,intent(in) :: g_num(nenod,nelmt),gdof_elmt(nedof,nelmt),mat_id(nelmt)
! only void elements
real(kreal),intent(in) :: dshape_hex8(ndim,ngnode,ngll),                        &
dlagrange_gll(ndim,ngll,ngll),gll_weights(ngll)
real(kreal),intent(in) :: stress_elmt(nst,ngll,nelmt)
real(kreal),intent(inout) :: extload(0:neq)

real(kreal) :: detjac,zero=0.0_kreal
real(kreal) :: coord(ngnode,ndim),jac(ndim,ndim),deriv(ndim,nenod),             &
bmat(nst,nedof),eld(nedof),bload(nedof),eload(nedof),sigma(nst)
integer :: egdof(nedof),num(nenod)
integer :: i,idof,i_elmt

! compute excavation load
do i_elmt=1,nelmt
  bload=zero; eld=zero
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  idof=0
  do i=1,ngll
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)
    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat(deriv,bmat)
    sigma=stress_elmt(:,i,i_elmt)
    eload=MATMUL(sigma,bmat)
    bload=bload+eload*detjac*gll_weights(i)
    ! since interpolation funtions are orthogonal we compute only nonzero terms
    idof=idof+3
    eld(idof)=eld(idof)+detjac*gll_weights(i)
    !eld(3:nedof:3)=eld(3:nedof:3)+lagrange_gll(i,:)*detjac*gll_weights(i)
  end do ! i=1,ngll
  extload(egdof)=extload(egdof)+eld*gam_blk(mat_id(i_elmt))+bload
enddo
extload(0)=zero
return
end subroutine excavation_load
!===============================================================================

! this subroutine computes the excavation loads. this load consists of both
! gravity and stress load
subroutine excavation_load_nodal(nelmt,gnod,g_num,mat_id,dshape_hex8,          &
dlagrange_gll,gll_weights,stress_elmt,excavload)
use global,only:ndim,nedof,nndof,nenod,ngll,ngnode,nnode,nst,g_coord,gam_blk
use preprocess,only:compute_bmat
use math_library,only:determinant,invert
implicit none
integer,intent(in) :: nelmt,gnod(8) ! nelmt (only void elements)
!integer,intent(in) :: g_num(nenod,nelmt),mat_id(nelmt)
integer,intent(in) :: g_num(:,:),mat_id(:)
real(kreal),intent(in) :: dshape_hex8(ndim,ngnode,ngll),                        &
dlagrange_gll(ndim,ngll,ngll),gll_weights(ngll)
real(kreal),intent(in) :: stress_elmt(nst,ngll,nelmt)
real(kreal),intent(out) :: excavload(nndof,nnode) ! nnode include all nodes

real(kreal) :: detjac,zero=0.0_kreal
real(kreal) :: coord(ngnode,ndim),jac(ndim,ndim),deriv(ndim,nenod), &
bmat(nst,nedof),eld(nedof),bload(nedof),eload(nedof),tload(nedof),sigma(nst)
integer :: num(nenod)
integer :: i,idof,i_elmt

excavload=zero
! compute excavation load
do i_elmt=1,nelmt
  bload=zero; eld=zero; tload=zero
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  idof=0
  do i=1,ngll
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)
    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat(deriv,bmat)
    sigma=stress_elmt(:,i,i_elmt)
    eload=MATMUL(sigma,bmat)
    bload=bload+eload*detjac*gll_weights(i)
    ! since interpolation funtions are orthogonal we compute only nonzero terms
    idof=idof+3
    eld(idof)=eld(idof)+detjac*gll_weights(i)
    !eld(3:nedof:3)=eld(3:nedof:3)+lagrange_gll(i,:)*detjac*gll_weights(i)
  end do ! i=1,ngll
  tload=eld*gam_blk(mat_id(i_elmt))+bload
  !excavload(:,num)=excavload(:,num)+reshape(tload,(/nndof,ngll/))
  do i=1,nndof
    excavload(i,num)=excavload(i,num)+tload(i:nedof:nndof)
  enddo
  !extload(egdof)=extload(egdof)+eld*gam_blk(mat_id(i_elmt))+bload
enddo
!extload(0)=zero
return
end subroutine excavation_load_nodal
!===============================================================================
end module excavation_library
!===============================================================================
