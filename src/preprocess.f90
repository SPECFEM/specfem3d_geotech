! this module contains preprocessing library routines
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Jul 07,2011
module preprocess

contains

!-------------------------------------------------------------------------------
! this subroutine computes strain-displacement matrix (B)
subroutine compute_bmat(deriv,bmat)
use set_precision
implicit none
real(kind=kreal),intent(in) :: deriv(:,:)
real(kind=kreal),intent(out) :: bmat(:,:)
integer :: n,nod,nst,j1,j2,j3
real(kind=kreal) :: dx,dy,dz

! check size
nst=ubound(bmat,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!',nst
  stop
endif
nod=ubound(deriv,2)

!bmat=0.0_kreal
!DO m=1,nod
!  n=3*m !j3
!  k=n-1 !j2
!  l=k-1 !j1
!  dx=deriv(1,m)
!  dy=deriv(2,m)
!  dz=deriv(3,m)
!  bmat(1,l)=dx
!  bmat(4,k)=dx
!  bmat(6,n)=dx
!  bmat(2,k)=dy
!  bmat(4,l)=dy
!  bmat(5,n)=dy
!  bmat(3,n)=dz
!  bmat(5,k)=dz
!  bmat(6,l)=dz
!end DO

bmat=0.0_kreal
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  dx=deriv(1,n)
  dy=deriv(2,n)
  dz=deriv(3,n)

  bmat(1,j1)=dx
  bmat(2,j2)=dy
  bmat(3,j3)=dz
  bmat(4,j1)=dy
  bmat(4,j2)=dx
  bmat(5,j2)=dz
  bmat(5,j3)=dy
  bmat(6,j1)=dz
  bmat(6,j3)=dx
enddo
!print*,'Bmat discrepancy:',maxval(abs(bmat-bbmat))
return
end subroutine compute_bmat
!===============================================================================

! this subrotine computes the stiffness matrix, diagonal preconditioner, and
! body loads (gravity and pseudostatic loads) optionally
! TODO: optional precoditioner,optional assembly of stiffness
subroutine stiffness_bodyload(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,gam,nu,ym, &
dshape_hex8,dlagrange_gll,gll_weights,storekm,dprecon,extload,gravity,pseudoeq)
use set_precision
use global,only:ndim,nst,ngll,nedof,nenod,ngnode,g_coord,eqkx,eqky,eqkz,nmatblk
use math_library,only:determinant,invert
use elastic,only:compute_cmat
implicit none
integer,intent(in) :: nelmt,neq,gnod(8)
!nelmt (only intact elements)
!integer,intent(in) :: g_num(nenod,nelmt),gdof_elmt(nedof,nelmt),mat_id(nelmt)
! only intact elements
integer,intent(in) :: g_num(:,:),gdof_elmt(:,:),mat_id(:)
! only intact elements
real(kind=kreal),intent(in) :: gam(nmatblk),nu(nmatblk),ym(nmatblk)
real(kind=kreal),intent(in) :: dshape_hex8(ndim,ngnode,ngll),                   &
dlagrange_gll(ndim,ngll,ngll),gll_weights(ngll)
real(kind=kreal),intent(out) :: storekm(nedof,nedof,nelmt),dprecon(0:neq)
real(kind=kreal),intent(inout),optional :: extload(0:neq)
logical,intent(in),optional :: gravity,pseudoeq

real(kind=kreal) :: detjac,zero=0.0_kreal
real(kind=kreal) :: cmat(nst,nst),coord(ngnode,ndim),jac(ndim,ndim),           &
deriv(ndim,nenod),bmat(nst,nedof),eld(nedof),eqload(nedof),km(nedof,nedof)
integer :: egdof(nedof),num(nenod)
integer :: i,idof,i_elmt,k

if(present(extload).and.(.not.present(gravity) .or. .not.present(pseudoeq)))then
  write(*,'(/,a)')'ERROR: both "gravity" and "pseudoeq" must be defined for &
  &"extload"!'
  stop
endif

! compute stiffness matrices
storekm=zero; dprecon=zero
do i_elmt=1,nelmt
  call compute_cmat(cmat,ym(mat_id(i_elmt)),nu(mat_id(i_elmt)))
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  km=zero; eld=zero; eqload=zero
  idof=0
  do i=1,ngll
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)

    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat(deriv,bmat)
    km=km+matmul(matmul(transpose(bmat),cmat),bmat)*detjac*gll_weights(i)
    idof=idof+3
    ! interpolation functions are orthogonal, hence it is simple
    eld(idof)=eld(idof)+detjac*gll_weights(i)
    !eld(3:nedof:3)=eld(3:nedof:3)+lagrange_gll(i,:)*detjac*gll_weights(i)
  end do ! i=1,ngll
  storekm(:,:,i_elmt)=km
  do k=1,nedof
    dprecon(egdof(k))=dprecon(egdof(k))+km(k,k)
  end do

  if(.not.present(extload))cycle
  ! compute body loads
  ! gravity load and add to extload
  if(gravity)extload(egdof)=extload(egdof)-eld*gam(mat_id(i_elmt))
  if(pseudoeq)then
    ! compute pseudostatic earthquake loads and add to extload
    eqload(1:nedof:3)=eqkx*eld(3:nedof:3)
    eqload(2:nedof:3)=eqky*eld(3:nedof:3)
    eqload(3:nedof:3)=eqkz*eld(3:nedof:3)
    extload(egdof)=extload(egdof)+eqload*gam(mat_id(i_elmt)) ! KN
  endif
end do ! i_elmt=1,nelmt
!write(*,*)'complete!'
dprecon(0)=zero
if(present(extload))extload(0)=zero
end subroutine stiffness_bodyload
!===============================================================================
end module preprocess
!===============================================================================

