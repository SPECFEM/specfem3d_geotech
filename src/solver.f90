! collection of solvers
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module solver
use set_precision
use global, only : cg_maxiter,cg_tol,g_num,nedof
use math_constants, only : zero

contains

! diagonally preconditioned conjuate-gradient solver
subroutine pcg_solver(myid,ngpart,maxngnode,neq,nelmt,k,u,f,dprecon,gdof_elmt,cg_iter,errcode,errtag)
implicit none
integer,intent(in) :: myid,ngpart,maxngnode,neq,nelmt ! nelmt (for intact) may not be same as global nelmt 
real(kind=kreal),dimension(nedof,nedof,nelmt),intent(in) :: k ! only for intact elements
real(kind=kreal),dimension(0:neq),intent(inout) :: u
real(kind=kreal),dimension(0:neq),intent(in) :: f,dprecon
!integer,dimension(nndof,nnode),intent(in) :: gdof
integer,dimension(nedof,nelmt),intent(in) :: gdof_elmt ! only for intact elements
integer,intent(out) :: cg_iter
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz
real(kind=kreal),dimension(0:neq) :: kp,p,r,z
real(kind=kreal),dimension(nedof,nedof) :: km

errtag="ERROR: unknown!"
errcode=-1

!---PCG solver
kp=zero
if(maxval(abs(u)).gt.zero)then 
  do i_elmt=1,nelmt   
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/)) 
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,u(egdof))
  end do
  kp(0)=zero
endif
r=f-kp
z=dprecon*r

p=z
!----pcg iteration----
pcg: do cg_iter=1,cg_maxiter
  kp=zero
  do i_elmt=1,nelmt      
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
    km=k(:,:,i_elmt)   
    kp(egdof)=kp(egdof)+matmul(km,p(egdof))
  end do
  kp(0)=zero
  
  rz=dot_product(r,z)
  alpha=rz/dot_product(p,kp)
  u=u+alpha*p
  
  if(abs(alpha)*maxval(abs(p))/maxval(abs(u)).le.cg_tol)then
    errcode=0
    return
  endif
  
  r=r-alpha*kp  
  z=dprecon*r
  beta=dot_product(r,z)/rz
  p=z+beta*p
  !write(*,'(i3,f25.18,f25.18,f25.18)')cg_iter,alpha,beta,rz
  
end do pcg
write(errtag,'(a)')'ERROR: PCG solver doesn''t converge!'
return
end subroutine pcg_solver
!============================================

end module solver
