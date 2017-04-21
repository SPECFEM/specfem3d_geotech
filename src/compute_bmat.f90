! This subroutine forms the strain-displacement matrix (bmat)
! REFERENCE:
!  copied and modified from
!  Smith and Griffiths (2004): Programming the finite element method
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
SUBROUTINE compute_bmat(bmat,deriv)
use set_precision
implicit none
real(kind=kreal),intent(in) :: deriv(:,:)
real(kind=kreal),intent(OUT) :: bmat(:,:)
integer :: k,l,m,n,nod,nst
real(kind=kreal) :: x,y,z

! check size
nst=ubound(bmat,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!'
  stop
endif
nod=ubound(deriv,2)

bmat=0.0_kreal
DO m=1,nod
  n=3*m
  k=n-1
  l=k-1
  x=deriv(1,m)
  y=deriv(2,m)
  z=deriv(3,m)
  bmat(1,l)=x
  bmat(4,k)=x
  bmat(6,n)=x
  bmat(2,k)=y
  bmat(4,l)=y
  bmat(5,n)=y
  bmat(3,n)=z
  bmat(5,k)=z
  bmat(6,l)=z
END DO
RETURN
END SUBROUTINE compute_bmat

