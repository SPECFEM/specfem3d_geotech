! this subroutine returns the elastic matrix for ih=3 (plane strain),
! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6 (3D)
! REFERENCE: 
!  copied and modified from 
!  Smith and Griffiths (2004): Programming the finite element method
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
subroutine compute_cmat(cmat,e,v)
use set_precision
 implicit none
 real(kind=kreal),intent(in)::e,v
 real(kind=kreal),intent(out)::cmat(:,:)
 real(kind=kreal)::v1,v2,c,vv,zero=0.0_kreal,pt5=0.5_kreal,one=1.0_kreal,two=2.0_kreal
 integer::i,ih
 cmat=zero  
 ih=ubound(cmat,1)
 v1=one-v
 c=e/((one+v)*(one-two*v))
 select case(ih)
 case(3)
   cmat(1,1)=v1*c
   cmat(2,2)=v1*c
   cmat(1,2)=v*c
   cmat(2,1)=v*c
   cmat(3,3)=pt5*c*(one-two*v)
 case(4)
   cmat(1,1)=v1*c
   cmat(2,2)=v1*c
   cmat(4,4)=v1*c
   cmat(3,3)=pt5*c*(one-two*v) 
   cmat(1,2)=v*c
   cmat(2,1)=v*c
   cmat(1,4)=v*c
   cmat(4,1)=v*c
   cmat(2,4)=v*c
   cmat(4,2)=v*c
 case(6)
   v2=v/(one-v)
   vv=(one-two*v)/(one-v)*pt5
   do i=1,3
     cmat(i,i)=one
   end do
   do i=4,6
     cmat(i,i)=vv
   end do
   cmat(1,2)=v2
   cmat(2,1)=v2
   cmat(1,3)=v2
   cmat(3,1)=v2
   cmat(2,3)=v2
   cmat(3,2)=v2
   cmat=cmat*e/(two*(one+v)*vv)
 case default
   write(*,*)'ERROR: wrong size for "cmat" matrix!'
 end select
return
end subroutine compute_cmat    
