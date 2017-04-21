! this module contains math constants
! math parameters
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module math_constants
use set_precision
implicit none
real(kind=kreal),parameter :: zero=0.0_kreal,half=0.5_kreal,one=1.0_kreal,     &
two=2.0_kreal
real(kind=kreal),parameter :: pi=3.141592653589793_kreal
real(kind=kreal),parameter :: deg2rad=pi/180.0_kreal,rad2deg=180.0_kreal/pi

! tolerance value for zero
real(kind=kreal),parameter :: inftol=1.0e32_kreal,zerotol = 1.0e-12_kreal
end module math_constants
!=======================================================

! this module coatins math routines
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module math_library
use set_precision

contains
function get_normal(x1,x2,x3) result(nx)
real(kind=kreal),dimension(3),intent(in) :: x1,x2,x3
real(kind=kreal),dimension(3) :: nx
real(kind=kreal),dimension(3) :: v1,v2
real(kind=kreal) :: norm

! two vectors
v1=x2-x1
v2=x3-x1

! cross product
nx(1)=v1(2)*v2(3)-v2(2)*v1(3)
nx(2)=v2(1)*v1(3)-v1(1)*v2(3)
nx(3)=v1(1)*v2(2)-v2(1)*v1(2)
norm=sqrt(sum(nx**2))
if(norm<=0.0_kreal)then
  write(*,*)'ERROR: undefined normal!'
  stop
endif
! unit normal
nx=nx/norm
return
end function get_normal
!=======================================================

function norm(x) result(l2n)
!
! this function calculates the l2 norm of vector x
!
implicit none
real(kind=kreal),intent(in) :: x(:)
real(kind=kreal)::l2n
l2n=sqrt(sum(x**2))
return
end function norm
!=======================================================

recursive function factorial(n) result(nfact)
implicit none
integer, intent(in) :: n
integer :: nfact
if(n > 0) then
  nfact = n * factorial(n-1)
  return
elseif (n==0)then
  nfact = 1
else
  write(*,*)'ERROR: undefined factorial!'
  stop
end if
end function factorial
!=======================================================

! this function returns the determinant of a 1x1, 2x2 or 3x3
! matrix.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
function determinant(xmat)result(det)
implicit none
real(kind=kreal),intent(in) :: xmat(:,:)
real(kind=kreal)::det
integer::n
n=ubound(xmat,1)
select case(n)
case(1)
  det=1.0_kreal
case(2)
  det=xmat(1,1)*xmat(2,2)-xmat(1,2)*xmat(2,1)
case(3)
  det=xmat(1,1)*(xmat(2,2)*xmat(3,3)-xmat(3,2)*xmat(2,3))
  det=det-xmat(1,2)*(xmat(2,1)*xmat(3,3)-xmat(3,1)*xmat(2,3))
  det=det+xmat(1,3)*(xmat(2,1)*xmat(3,2)-xmat(3,1)*xmat(2,2))
case default
  write(*,*)'ERROR: determinant of ',n,'X',n,' matrix not supported!'
end select
return
end function determinant
!=======================================================

! this subroutine inverts a small square matrix onto itself.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine invert(xmat)
implicit none
real(kind=kreal),intent(in out)::xmat(:,:)
real(kind=kreal)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
integer::ndim,i,k
ndim=ubound(xmat,1)
if(ndim==2)then
  det=xmat(1,1)*xmat(2,2)-xmat(1,2)*xmat(2,1)
  j11=xmat(1,1)
  xmat(1,1)=xmat(2,2)
  xmat(2,2)=j11
  xmat(1,2)=-xmat(1,2)
  xmat(2,1)=-xmat(2,1)
  xmat=xmat/det
else if(ndim==3)then
  det=xmat(1,1)*(xmat(2,2)*xmat(3,3)-xmat(3,2)*xmat(2,3))
  det=det-xmat(1,2)*(xmat(2,1)*xmat(3,3)-xmat(3,1)*xmat(2,3))
  det=det+xmat(1,3)*(xmat(2,1)*xmat(3,2)-xmat(3,1)*xmat(2,2))
  j11=xmat(2,2)*xmat(3,3)-xmat(3,2)*xmat(2,3)
  j21=-xmat(2,1)*xmat(3,3)+xmat(3,1)*xmat(2,3)
  j31=xmat(2,1)*xmat(3,2)-xmat(3,1)*xmat(2,2)
  j12=-xmat(1,2)*xmat(3,3)+xmat(3,2)*xmat(1,3)
  j22=xmat(1,1)*xmat(3,3)-xmat(3,1)*xmat(1,3)
  j32=-xmat(1,1)*xmat(3,2)+xmat(3,1)*xmat(1,2)
  j13=xmat(1,2)*xmat(2,3)-xmat(2,2)*xmat(1,3)
  j23=-xmat(1,1)*xmat(2,3)+xmat(2,1)*xmat(1,3)
  j33=xmat(1,1)*xmat(2,2)-xmat(2,1)*xmat(1,2)

  xmat(1,1)=j11; xmat(1,2)=j12; xmat(1,3)=j13
  xmat(2,1)=j21; xmat(2,2)=j22; xmat(2,3)=j23
  xmat(3,1)=j31; xmat(3,2)=j32; xmat(3,3)=j33
  xmat=xmat/det
else
  do k=1,ndim
    con=xmat(k,k)
    xmat(k,k)=1.0_kreal
    xmat(k,:)=xmat(k,:)/con
    do i=1,ndim
      if(i/=k)then
        con=xmat(i,k)
        xmat(i,k)=0.0_kreal
        xmat(i,:)=xmat(i,:)-xmat(k,:)*con
      end if
    end do
  end do
end if
return
end subroutine invert
!=======================================================

! this subroutine forms the stress invariants 3d.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine stress_invariant(stress,sigm,dsbar,theta)
implicit none
real(kind=kreal),intent(in)::stress(:)
real(kind=kreal),intent(out),optional::sigm,dsbar,theta
real(kind=kreal)::sx,sy,sz,txy,dx,dy,dz,xj3,sine,s1,s2,s3,s4,s5,s6,ds1,ds2,  &
ds3,d2,d3,sq3,zero=0.0_kreal,small=1.e-12_kreal,one=1.0_kreal,two=2.0_kreal, &
three=3.0_kreal,six=6.0_kreal,thpt5=13.5_kreal
integer :: nst

! check size
nst=ubound(stress,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!'
  stop
endif

sq3=sqrt(three)
s1=stress(1)
s2=stress(2)
s3=stress(3)
s4=stress(4)
s5=stress(5)
s6=stress(6)
sigm=(s1+s2+s3)/three
d2=((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/six+s4*s4+s5*s5+s6*s6

if(d2<small)d2=small ! special case of hydrostatic pressure or just at the tip

ds1=s1-sigm
ds2=s2-sigm
ds3=s3-sigm
d3=ds1*ds2*ds3-ds1*s5*s5-ds2*s6*s6-ds3*s4*s4+two*s4*s5*s6
dsbar=sq3*sqrt(d2)
if(dsbar<small)then
  theta=zero
else
  sine=-three*sq3*d3/(two*sqrt(d2)**3)
  if(sine>=one)sine=one
  if(sine<-one)sine=-one
  theta=asin(sine)/three
end if
return
end subroutine stress_invariant
!=======================================================

! quick sort of integer list
function quick_sort(x,n) result(xnew)
integer,intent(in) :: n ! size of the vector data x
integer, dimension(n) :: x ! data vector to sort
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
end function quick_sort
!=======================================================

! quick sort of real list
function rquick_sort(x,n) result(xnew)
integer,intent(in) :: n ! size of the vector data x
real(kind=kreal), dimension(n) :: x ! data vector to sort
real(kind=kreal) :: temp
integer :: i,j
real(kind=kreal),dimension(n) :: xnew

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
end function rquick_sort
!=======================================================

! insertion sort of integer list
subroutine insertion_sort(x,n)
integer,intent(in) :: n ! size of the vector data x
real, intent(inout), dimension(n) :: x ! data vector to sort
real :: temp
integer :: i, j

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

end module math_library
