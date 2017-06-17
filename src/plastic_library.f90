! this module contains the routine for Mohr-Coulomb plasticity, Viscoplastic
! algorithm, and strength reduction techniques
! AUTHOR
!   Hom Nath Gharti
! REVISION
!   HNG, Jul 07,2011; HNG, Apr 09,2010
module plastic_library
use set_precision
use math_constants
use conversion_constants,only:DEG2RAD,RAD2DEG

contains

!-------------------------------------------------------------------------------
! this function computes the stable pseudo-time step
! for viscoplastic algorithm (Cormeau 1975, Smith and Griffiths 2003)
function dt_viscoplas(nmat,nuf,phif,ymf,ismat) result(dt_min)
implicit none
integer,intent(in) :: nmat
! friction angle,Poisson's ratio, strength reduction factor
real(kind=kreal),intent(in) :: ymf(nmat),phif(nmat),nuf(nmat) ! phif in degrees
logical,optional,intent(in) :: ismat(nmat)
real(kind=kreal) :: dt_min
real(kind=kreal) :: dt,snphi
real(kind=kreal),parameter :: r4=4.0_kreal
integer :: i_mat
logical :: ismat_on(nmat)
ismat_on=.true.
if(present(ismat))ismat_on=ismat
! compute minimum pseudo-time step for viscoplasticity
dt_min=inftol
do i_mat=1,nmat
  if(.not.ismat_on(i_mat))cycle
  snphi=sin(phif(i_mat)*deg2rad)
  dt=r4*(one+nuf(i_mat))*(one-two*nuf(i_mat))/(ymf(i_mat)*(one-two*nuf(i_mat)+ &
  snphi**2))
  if(dt<dt_min)dt_min=dt
end do
end function dt_viscoplas
!===============================================================================

subroutine strength_reduction(srf,phinu,nmat,coh,nu,phi,psi,cohf,nuf,phif,psif,&
istat)
implicit none
real(kind=kreal),intent(in) :: srf
logical :: phinu
integer,intent(in) :: nmat
real(kind=kreal),intent(in) :: coh(nmat),nu(nmat),phi(nmat),psi(nmat)
! phi,psi in degrees
real(kind=kreal),intent(out) :: cohf(nmat),nuf(nmat),phif(nmat),psif(nmat)
! phif,psif in degrees
! status whether the material properties has changed
integer,intent(out) :: istat

integer :: i_mat
real(kind=kreal) :: beta_nuphi,omtnu,snphi,snphif,tnphi,tnpsi

cohf=coh; nuf=nu; phif=phi; psif=psi
istat=0

! strength reduction
if(srf/=one)then
  do i_mat=1,nmat
    tnphi=tan(phi(i_mat)*deg2rad)
    phif(i_mat)=atan(tnphi/srf)*rad2deg
    tnpsi=tan(psi(i_mat)*deg2rad)
    psif(i_mat)=atan(tnpsi/srf)*rad2deg
    cohf(i_mat)=coh(i_mat)/srf
  enddo
endif

if(.not.phinu)return
! correction for phi-nu inequality sin(phi)>=1-2*nu
! Reference: Zheng et al 2005, IJNME
! currently only the Poisson's ratio is corrected
do i_mat=1,nmat
  omtnu=one-two*nu(i_mat)
  snphif=sin(phif(i_mat)*deg2rad)
  if(snphif<omtnu)then
    snphi=sin(phi(i_mat)*deg2rad)
    beta_nuphi=snphi/omtnu
    if(beta_nuphi<one)beta_nuphi=one+zerotol
    nuf(i_mat)=half*(one-snphif/beta_nuphi)
    istat=1 ! material properties has changed
  endif
enddo
return
end subroutine strength_reduction
!===============================================================================

! this subroutine calculates the value of the yield function
! for a mohr-coulomb material (phi in degrees, theta in radians).
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine mohcouf(phi,c,sigm,dsbar,theta,f)
implicit none
real(kind=kreal),intent(in)::phi,c,sigm,dsbar,theta
real(kind=kreal),intent(out)::f
real(kind=kreal)::phir,snph,csph,csth,snth,r3=3.0_kreal

phir=phi*deg2rad
snph=sin(phir)
csph=cos(phir)
csth=cos(theta)
snth=sin(theta)
f=snph*sigm+dsbar*(csth/sqrt(r3)-snth*snph/r3)-c*csph
return
end subroutine mohcouf
!===============================================================================

! this subroutine forms the derivatives of a mohr-coulomb potential
! function with respect to the three stress invariants
! (psi in degrees, theta in radians).
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine mohcouq(psi,dsbar,theta,dq1,dq2,dq3)
implicit none
real(kind=kreal),intent(in)::psi,dsbar,theta
real(kind=kreal),intent(out)::dq1,dq2,dq3
real(kind=kreal)::psir,snth,snps,sq3,c1,csth,cs3th,tn3th,tnth,pt49=0.49_kreal,&
pt5=0.5_kreal,r3=3.0_kreal

psir=psi*deg2rad
snth=sin(theta)
snps=sin(psir)
sq3=sqrt(r3)
dq1=snps

if(abs(snth).gt.pt49)then
  c1=one
  if(snth.lt.zero)c1=-one
  dq2=(sq3*pt5-c1*snps*pt5/sq3)*sq3*pt5/dsbar
  dq3=zero
else
  csth=cos(theta)
  cs3th=cos(r3*theta)
  tn3th=tan(r3*theta)
  tnth=snth/csth
  dq2=sq3*csth/dsbar*((one+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*pt5
  dq3=pt5*r3*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar)
end if
return
end subroutine mohcouq
!===============================================================================

! this subroutine forms the derivatives of the invariants with respect to
! stress in 3d.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine formm(stress,m1,m2,m3)
implicit none
real(kind=kreal),intent(in)::stress(:)
real(kind=kreal),intent(out)::m1(:,:),m2(:,:),m3(:,:)
real(kind=kreal)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,  &
r3=3.0_kreal,r6=6.0_kreal
integer::nst,i,j

nst=ubound(stress,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!'
  stop
endif

sx=stress(1);  sy=stress(2);  sz=stress(3)
txy=stress(4); tyz=stress(5); tzx=stress(6)
sigm=(sx+sy+sz)/r3
dx=sx-sigm; dy=sy-sigm; dz=sz-sigm

m1=zero; m2=zero
m1(1:3,1:3)=one/(r3*sigm)
do i=1,3
  m2(i,i)=two
  m2(i+3,i+3)=r6
end do
m2(1,2)=-one
m2(1,3)=-one
m2(2,3)=-one
m3(1,1)=dx
m3(1,2)=dz
m3(1,3)=dy
m3(1,4)=txy
m3(1,5)=-two*tyz
m3(1,6)=tzx
m3(2,2)=dy
m3(2,3)=dx
m3(2,4)=txy
m3(2,5)=tyz
m3(2,6)=-two*tzx
m3(3,3)=dz
m3(3,4)=-two*txy
m3(3,5)=tyz
m3(3,6)=tzx
m3(4,4)=-r3*dz
m3(4,5)=r3*tzx
m3(4,6)=r3*tyz
m3(5,5)=-r3*dx
m3(5,6)=r3*txy
m3(6,6)=-r3*dy
do i=1,6
  do j=i+1,6
    m1(j,i)=m1(i,j)
    m2(j,i)=m2(i,j)
    m3(j,i)=m3(i,j)
  end do
end do
m1=m1/r3; m2=m2/r3; m3=m3/r3
return
end subroutine formm
!===============================================================================

end module plastic_library
!===============================================================================
