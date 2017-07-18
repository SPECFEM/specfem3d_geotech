! collection of solvers
! AUTHOR
!   Hom Nath Gharti
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module solver_mpi
use set_precision
use global, only : cg_maxiter,cg_tol,g_num,nedof,nndof,nnode
use math_constants, only : zero
!use math_library !, only : dot_product_par,maxvec
use math_library_mpi
use ghost_library_mpi

contains

! diagonally preconditioned conjuate-gradient solver
subroutine pcg_solver(neq,nelmt,k,u_g,f,dprecon_g,gdof_elmt,cg_iter,errcode,   &
errtag)
!use math_library
implicit none
integer,intent(in) :: neq,nelmt
! nelmt (for intact) may not be same as global nelmt
real(kind=kreal),dimension(nedof,nedof,nelmt),intent(in) :: k
! only for intact elements
real(kind=kreal),dimension(0:neq),intent(inout) :: u_g
real(kind=kreal),dimension(0:neq),intent(in) :: f,dprecon_g
!integer,dimension(nedof,nelmt),intent(in) :: gdof_elmt
integer,dimension(:,:),intent(in) :: gdof_elmt
! only for intact elements
integer,intent(out) :: cg_iter
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz
real(kind=kreal),dimension(0:neq) :: kp,p,p_g,r,z,z_g
real(kind=kreal),dimension(nedof,nedof) :: km

errtag="ERROR: unknown!"
errcode=-1

kp=zero
if(maxval(abs(u_g)).gt.zero)then
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt)
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,u_g(egdof))
  enddo
  kp(0)=zero
endif
r=f-kp
z=dprecon_g*r

call assemble_ghosts(neq,z,z_g)

p=z
! pcg iteration
pcg: do cg_iter=1,cg_maxiter
  call assemble_ghosts(neq,p,p_g)
  kp=zero
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt)
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,p_g(egdof))
  enddo
  kp(0)=zero

  rz=dot_product_par(r,z_g)
  alpha=rz/dot_product_par(p_g,kp)
  u_g=u_g+alpha*p_g

  if(abs(alpha)*maxvec(abs(p_g))/maxvec(abs(u_g)).le.cg_tol)then
    errcode=0
    return
  endif
  !if(myrank==0)print*,abs(alpha)*maxvec(abs(p_g))/maxvec(abs(u_g))
  r=r-alpha*kp
  z=dprecon_g*r
  call assemble_ghosts(neq,z,z_g) !,gdof)
  beta=dot_product_par(r,z_g)/rz
  p=z+beta*p
  !if(myrank==0)write(*,'(i3,f25.18,f25.18,f25.18)')cg_iter,alpha,beta,rz
enddo pcg
write(errtag,'(a)')'ERROR: PCG solver doesn''t converge!'
return
end subroutine pcg_solver
!===============================================================================

end module solver_mpi
!===============================================================================
