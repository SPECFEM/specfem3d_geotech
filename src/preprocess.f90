! this module contains preprocessing library routines
! REVISION:
!   HNG, Jul 07,2011
module preprocess
! interface to auxilliary routines
interface
! this subroutine computes B matrix
subroutine compute_bmat(bmat,deriv)
 use set_precision
 implicit none
 real(kind=kreal),intent(in)::deriv(:,:)
 real(kind=kreal),intent(out)::bmat(:,:)
end subroutine compute_bmat

! this subroutine computes C matrix (elasticity matrix)
subroutine compute_cmat(cmat,e,v)
 use set_precision
 implicit none
 real(kind=kreal),intent(in)::e,v
 real(kind=kreal),intent(out)::cmat(:,:)
end subroutine compute_cmat
end interface

contains
! this subrotine computes the stiffness matrix, diagonal preconditioner, and
! body loads (gravity and pseudostatic loads) optionally
! TODO: optional precoditioner,optional assembly of stiffness
subroutine stiffness_bodyload(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,gam,nu,ym, &
dshape_hex8,dlagrange_gll,gll_weights,storkm,dprecon,extload,gravity,pseudoeq)
use set_precision
use global,only:ndim,nst,ngll,nedof,nenod,ngnod,g_coord,eqkx,eqky,eqkz,nmat
use math_library,only:determinant,invert
implicit none
integer,intent(in) :: nelmt,neq,gnod(8) ! nelmt (only intact elements)
integer,intent(in) :: g_num(nenod,nelmt),gdof_elmt(nedof,nelmt),mat_id(nelmt) ! only intact elements
real(kind=kreal),intent(in) :: gam(nmat),nu(nmat),ym(nmat)
real(kind=kreal),intent(in) :: dshape_hex8(ndim,ngnod,ngll),                   &
dlagrange_gll(ndim,ngll,ngll),gll_weights(ngll)
real(kind=kreal),intent(out) :: storkm(nedof,nedof,nelmt),dprecon(0:neq)
real(kind=kreal),intent(inout),optional :: extload(0:neq)
logical,intent(in),optional :: gravity,pseudoeq

real(kind=kreal) :: detjac,zero=0.0_kreal
real(kind=kreal) :: cmat(nst,nst),coord(ngnod,ndim),jac(ndim,ndim),deriv(ndim,nenod), &
bmat(nst,nedof),eld(nedof),eqload(nedof),km(nedof,nedof)
integer :: egdof(nedof),num(nenod)
integer :: i,idof,i_elmt,k

if(present(extload).and.(.not.present(gravity) .or. .not.present(pseudoeq)))then
  write(*,'(/,a)')'ERROR: both "gravity" and "pseudoeq" must be defined for "extload"!'
  stop
endif

! compute stiffness matrices
storkm=zero; dprecon=zero
!----element stiffness integration, storage and preconditioner----
do i_elmt=1,nelmt
  call compute_cmat(cmat,ym(mat_id(i_elmt)),nu(mat_id(i_elmt)))
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod))) !transpose(g_coord(:,num(1:ngnod)))
  !print*,ngllx,nglly,ngllz,ngll,nndof,nenod ,gdof(:,g_num(:,ielmt))
  egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,ielmt)),(/nndof*nenod/)) !g=g_g(:,ielmt)
  km=zero; eld=zero; eqload=zero
  idof=0
  do i=1,ngll
    !call shape_function(fun,gll_points(i))
    ! compute Jacobian at GLL point using 20 noded element
    !call shape_derivative(der,gll_points(:,i))
    jac=matmul(dshape_hex8(:,:,i),coord) !jac=matmul(der,coord)
    detjac=determinant(jac)
    call invert(jac)

    deriv=matmul(jac,dlagrange_gll(:,i,:)) ! use der for gll
    call compute_bmat(bmat,deriv) !!! gll bmat matrix
    km=km+matmul(matmul(transpose(bmat),cmat),bmat)*detjac*gll_weights(i)
    idof=idof+3
    ! interpolation functions are orthogonal, hence it is simple
    eld(idof)=eld(idof)+detjac*gll_weights(i)
    !eld(3:nedof:3)=eld(3:nedof:3)+lagrange_gll(i,:)*detjac*gll_weights(i)
  end do ! i=1,ngll
  storkm(:,:,i_elmt)=km
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
end module preprocess

