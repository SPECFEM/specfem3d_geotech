! This subroutine forms the strain-displacement matrix (bmat)
! in 2D (ih=3 or 4) or 3D (ih=6)
! REFERENCE:
!  copied and modified from
!  Smith and Griffiths (2004): Programming the finite element method
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
SUBROUTINE compute_bmat(bmat,deriv)
use set_precision
 IMPLICIT NONE
 REAL(kind=kreal),INTENT(IN)::deriv(:,:)
 REAL(kind=kreal),INTENT(OUT)::bmat(:,:)
 INTEGER::k,l,m,n,ih,nod
 REAL(kind=kreal)::x,y,z
 bmat=0.0_kreal
 ih=UBOUND(bmat,1)
 nod=UBOUND(deriv,2)
 SELECT CASE (ih)
 CASE(3,4)
   DO m=1,nod
     k=2*m
     l=k-1
     x=deriv(1,m)
     y=deriv(2,m)
     bmat(1,l)=x
     bmat(3,k)=x
     bmat(2,k)=y
     bmat(3,l)=y
   END DO
 CASE(6)
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
 CASE DEFAULT
   WRITE(*,*)'ERROR: wrong dimension for "nst" in bmat matrix!'
 END SELECT
RETURN
END SUBROUTINE compute_bmat

