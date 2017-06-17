! This module contains elastic routines
! AUTHOR
!   Hom Nath Gharti
! REVISION
!   HNG, Jul 18,2012; Jul 12,2011; HNG, Apr 09,2010
module elastic
use set_precision
contains
! this subroutine computes the elastic matrix (3D)
! REFERENCE:
!  Zienkiewicz, Taylor, and Zhu (2005): "The finite element method: its basis &
!  fundamentals", Sixth edition, pp195
subroutine compute_cmat(cmat,E,v)
!use set_precision
implicit none
real(kind=kreal),intent(in) :: E,v ! young's modulus, poisson's ratio
real(kind=kreal),intent(out) :: cmat(:,:) ! elasticity matrix

real(kind=kreal) :: zero=0.0_kreal,half=0.5_kreal,one=1.0_kreal,two=2.0_kreal
real(kind=kreal) :: Efact,onepv,onemv,onem2v
integer :: nst

! check size
nst=ubound(cmat,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!'
  stop
endif

cmat=zero

onepv=one+v; onemv=one-v; onem2v=one-two*v
Efact=E/(onepv*onem2v)
cmat(1,1)=onemv*Efact; cmat(1,2)=v*Efact  ; cmat(1,3)=cmat(1,2)
cmat(2,1)=cmat(1,2)  ; cmat(2,2)=cmat(1,1); cmat(2,3)=cmat(1,2)
cmat(3,1)=cmat(1,3)  ; cmat(3,2)=cmat(2,3); cmat(3,3)=cmat(1,1)
cmat(4,4)=half*onem2v*Efact;
cmat(5,5)=cmat(4,4)
cmat(6,6)=cmat(4,4)

return
end subroutine compute_cmat
!===============================================================================

! this subroutine returns the elastic matrix (3D)
! REFERENCE:
!  copied and modified from
!  Smith and Griffiths (2004): Programming the finite element method
subroutine compute_cmatOLD(cmat,e,v)
implicit none
real(kind=kreal),intent(in)::e,v
real(kind=kreal),intent(out)::cmat(:,:)
real(kind=kreal)::v1,v2,c,vv,zero=0.0_kreal,pt5=0.5_kreal,one=1.0_kreal,       &
two=2.0_kreal
integer::i,nst

! check size
nst=ubound(cmat,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!'
  stop
endif

cmat=zero
v1=one-v
c=e/((one+v)*(one-two*v))

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
return
end subroutine compute_cmatOLD
!===============================================================================

subroutine compute_cmat_elastic(bulkmod,shearmod,cmat)
implicit none
real(kind=kreal),intent(in) :: bulkmod,shearmod !K,mu
real(kind=kreal),intent(out) :: cmat(:,:)

integer :: nst
real(kind=kreal) :: twomuthird

! check size
nst=ubound(cmat,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!'
  stop
endif

twomuthird=2.0_kreal*shearmod/3.0_kreal

cmat=0.0_kreal

cmat(1,1) = bulkmod+2.0_kreal*twomuthird ! C1111
cmat(1,2) = bulkmod-twomuthird     ! C1122
cmat(1,3) = cmat(1,2)              ! C1133
!cmat(1,4) = 0                     ! C1112
!cmat(1,5) = 0                     ! C1123
!cmat(1,6) = 0                     ! C1113

cmat(2,1) = cmat(1,2)              ! C2211
cmat(2,2) = cmat(1,1)              ! C2222
cmat(2,3) = cmat(1,2)              ! C2233
!cmat(2,4) = 0;                    ! C2212
!cmat(2,5) = 0;                    ! C2223
!cmat(2,6) = 0;                    ! C2213

cmat(3,1) = cmat(1,3)              ! C3311
cmat(3,2) = cmat(2,3)              ! C3322
cmat(3,3) = cmat(1,1)              ! C3333
!cmat(3,4) = 0;                    ! C3312
!cmat(3,5) = 0;                    ! C3323
!cmat(3,6) = 0;                    ! C3313

!cmat(4,1) = 0;                    ! C1211
!cmat(4,2 = 0;                     ! C1222
!cmat(4,3) = 0;                    ! C1233
cmat(4,4) = shearmod               ! C1212
!cmat(4,5) = 0;                    ! C1223
!cmat(4,6) = 0;                    ! C1213

!cmat(5,1) = 0;                    ! C2311
!cmat(5,2) = 0;                    ! C2322
!cmat(5,3) = 0;                    ! C2333
!cmat(5,4) = 0;                    ! C2312
cmat(5,5) = cmat(4,4)              ! C2323
!cmat(5,6) = 0                     ! C2313

!cmat(6,1) = 0;                    ! C1311
!cmat(6,2) = 0;                    ! C1322
!cmat(6,3) = 0;                    ! C1333
!cmat(6,4) = 0;                    ! C1312
!cmat(6,5) = 0;                    ! C1323
cmat(6,6) = cmat(4,4)              ! C1313

return
end subroutine compute_cmat_elastic
!===========================================

!! compute elastic stress
!subroutine compute_stress_elastic(nelmt,nedofu,imapu,gnod,g_num,gdof_elmt,elmt_bulkmod,       &
!elmt_shearmod,nodalu,dshape_hex8,dlagrange_gll,gll_weights,islast,stress_elas,strain_elas,bodyload)
!use global,only:g_coord,ndim,nst,ngnod,ngll,nedof,nnode,nndof,idofu
!use math_library,only:determinant,invert
!use weakform_residual,only:compute_rmat_stress
!implicit none
!integer,intent(in) :: nelmt,nedofu,imapu(nedofu),gnod(ngnod),g_num(ngll,nelmt),gdof_elmt(nedof,nelmt)
!real(kind=kreal),intent(in) :: elmt_bulkmod(ngll,nelmt),elmt_shearmod(ngll,nelmt)
!real(kind=kreal),intent(in) :: nodalu(nndof,nnode),dshape_hex8(:,:,:),dlagrange_gll(:,:,:),gll_weights(:)
!logical,intent(in) :: islast
!real(kind=kreal),intent(inout) :: stress_elas(:,:,:),strain_elas(:,:,:),bodyload(:)
!
!real(kind=kreal) :: detjac
!real(kind=kreal) :: bload(nedofu),bmat(nst,nedofu),eload(nedofu),cmat(nst,nst),coord(ngnod,ndim), &
!dinterpf(ndim,ngll),eld(nedofu),jac(ndim,ndim),epst(nst),sigma(nst),bulkmod(ngll),shearmod(ngll)
!integer :: egdofu(nedofu),num(ngll)
!integer :: i,ielmt
!
!do ielmt=1,nelmt
!  num=g_num(:,ielmt)
!  coord=transpose(g_coord(:,num(gnod))) !transpose(g_coord(:,num(1:ngnod)))
!  egdofu=gdof_elmt(imapu,ielmt) !reshape(gdof(:,g_num(:,ielmt)),(/nedof/)) !g=g_g(:,i_elmt)
!
!  bulkmod=elmt_bulkmod(:,ielmt)
!  shearmod=elmt_shearmod(:,ielmt)
!
!  !eld=reshape(nodalu(imapu,g_num(:,ielmt)),(/nedofu/)) !x(egdof)
!  eld=reshape(nodalu(idofu,g_num(:,ielmt)),(/nedofu/)) !x(egdof)
!  bload=0.0_kreal
!  do i=1,ngll ! loop over integration points
!    call compute_cmat_elastic(bulkmod(i),shearmod(i),cmat)
!    !call shape_derivative(der,gll_points(:,i))
!    jac=matmul(dshape_hex8(:,:,i),coord) !jac=matmul(der,coord)
!    detjac=determinant(jac)
!    call invert(jac)
!
!    dinterpf=matmul(jac,dlagrange_gll(:,i,:))
!    call compute_rmat_stress(dinterpf,bmat)
!    !call compute_bmat(dinterpf,bmat)
!    epst=matmul(bmat,eld)
!    !strain_local(:,i,ielmt)=epst !+strain_local(:,i,ielmt)
!    sigma=matmul(cmat,epst)
!    !print*,'Elastic stress',maxval(abs(sigma)); stop
!    stress_elas(:,i,ielmt)=sigma !+stress_elmt(:,i,ielmt)
!
!    if(islast)then
!      strain_elas(:,i,ielmt)=epst
!    endif
!
!    eload=matmul(sigma,bmat)
!    bload=bload+eload*detjac*gll_weights(i)
!  enddo ! i
!  !----compute the total bodyloads vector----
!  bodyload(egdofu)=bodyload(egdofu)+bload
!enddo
!return
!end subroutine compute_stress_elastic
!!===========================================
!
!subroutine overburden_stress(nelmt,g_num,mat_id,z0,p0,k0,stress_elmt)
!use global,only:nenod,ngll,nmatblk,nst,g_coord,rho
!implicit none
!integer,intent(in) :: nelmt
!integer,intent(in) :: g_num(nenod,nelmt),mat_id(nelmt)
!real(kind=kreal),intent(in) :: z0,p0,k0
!real(kind=kreal),intent(inout) :: stress_elmt(nst,ngll,nelmt)
!
!real(kind=kreal) :: z,szz,gam(nmatblk)
!integer :: i,i_elmt,num(nenod)
!
!gam=9.81_kreal*rho
!do i_elmt=1,nelmt
!  num=g_num(:,i_elmt)
!
!  do i=1,ngll ! loop over integration points
!    z=g_coord(3,num(i))
!    szz=p0-gam(mat_id(i_elmt))*abs(z-z0) ! compression -
!    stress_elmt(3,i,i_elmt)=szz
!    stress_elmt(1,i,i_elmt)=k0*szz
!    stress_elmt(2,i,i_elmt)=k0*szz
!  end do ! i GLL
!end do ! i_elmt
!return
!end subroutine overburden_stress
!!===========================================
!
!! TODO: is it possible to compute stress_elmt only for intact elements just in case?
!! it seems that the subarray cannot be a receiving array
!! this subroutine computes elastic stress from the known displacement
!subroutine elastic_stress(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,dshape_hex8,dlagrange_gll,x,stress_elmt)
!use global,only:ndim,nedof,nenod,ngnod,ngll,nst,g_coord,ym,nu
!!use preprocess,only:compute_bmat,compute_cmat
!use math_library,only:determinant,invert
!use weakform_residual,only:compute_rmat_stress
!implicit none
!integer,intent(in) :: nelmt,neq,gnod(8)
!integer,intent(in) :: g_num(nenod,nelmt),gdof_elmt(nedof,nelmt),mat_id(nelmt)
!real(kind=kreal),intent(in) :: dshape_hex8(ndim,ngnod,ngll),dlagrange_gll(ndim,ngll,ngll),x(0:neq)
!real(kind=kreal),intent(inout) :: stress_elmt(nst,ngll,nelmt)
!
!real(kind=kreal) :: detjac,zero=0.0_kreal
!real(kind=kreal) :: cmat(nst,nst),coord(ngnod,ndim),jac(ndim,ndim),dinterpf(ndim,nenod), &
!bmat(nst,nedof),eld(nedof),eps(nst),sigma(nst)
!integer :: egdof(nedof),num(nenod)
!integer :: i,i_elmt
!
!do i_elmt=1,nelmt
!  call compute_cmatOLD(cmat,ym(mat_id(i_elmt)),nu(mat_id(i_elmt)))
!  num=g_num(:,i_elmt)
!  coord=transpose(g_coord(:,num(gnod))) !transpose(g_coord(:,num(1:ngnod)))
!  egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/)) !g=g_g(:,i_elmt)
!  eld=x(egdof)
!
!  do i=1,ngll ! loop over integration points
!    !call shape_derivative(der,gll_points(:,i))
!    jac=matmul(dshape_hex8(:,:,i),coord) !jac=matmul(der,coord)
!    detjac=determinant(jac)
!    call invert(jac)!
!
!    dinterpf=matmul(jac,dlagrange_gll(:,i,:))
!    call compute_rmat_stress(bmat,dinterpf)
!    !call compute_bmat(deriv,bmat)
!    eps=matmul(bmat,eld)
!    sigma=matmul(cmat,eps)
!    stress_elmt(:,i,i_elmt)=sigma
!  end do ! i GLL
!end do ! i_elmt
!return
!end subroutine elastic_stress
!!===========================================
!
!subroutine elastic_stress_intact(nelmt_intact,neq,gnod,elmt_intact,g_num,gdof_elmt, &
!mat_id,dshape_hex8,dlagrange_gll,x,stress_elmt)
!use global,only:ndim,nedof,nelmt,nenod,ngnod,ngll,nst,g_coord,ym,nu
!!use preprocess,only:compute_bmat,compute_cmat
!use math_library,only:determinant,invert
!use weakform_residual,only:compute_rmat_stress
!implicit none
!integer,intent(in) :: nelmt_intact,neq,gnod(8)
!integer,intent(in) :: elmt_intact(nelmt_intact),g_num(nenod,nelmt_intact), &
!gdof_elmt(nedof,nelmt_intact),mat_id(nelmt_intact)
!real(kind=kreal),intent(in) :: dshape_hex8(ndim,ngnod,ngll),dlagrange_gll(ndim,ngll,ngll),x(0:neq)
!real(kind=kreal),intent(inout) :: stress_elmt(nst,ngll,nelmt)
!
!real(kind=kreal) :: detjac,zero=0.0_kreal
!real(kind=kreal) :: cmat(nst,nst),coord(ngnod,ndim),jac(ndim,ndim),dinterpf(ndim,nenod), &
!bmat(nst,nedof),eld(nedof),eps(nst),sigma(nst)
!integer :: egdof(nedof),num(nenod)
!integer :: i,i_elmt,ielmt
!
!do i_elmt=1,nelmt_intact
!  ielmt=elmt_intact(i_elmt)
!  call compute_cmatOLD(cmat,ym(mat_id(i_elmt)),nu(mat_id(i_elmt)))
!  num=g_num(:,i_elmt)
!  coord=transpose(g_coord(:,num(gnod))) !transpose(g_coord(:,num(1:ngnod)))
!  egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/)) !g=g_g(:,i_elmt)
!  eld=x(egdof)
!  !print*,egdof
!  !stop
!  do i=1,ngll ! loop over integration points
!    !call shape_derivative(der,gll_points(:,i))
!    jac=matmul(dshape_hex8(:,:,i),coord) !jac=matmul(der,coord)
!    detjac=determinant(jac)
!    call invert(jac)!
!
!    dinterpf=matmul(jac,dlagrange_gll(:,i,:))
!    call compute_rmat_stress(bmat,dinterpf)
!    !call compute_bmat(deriv,bmat)
!    eps=matmul(bmat,eld)
!    sigma=matmul(cmat,eps)
!    stress_elmt(:,i,ielmt)=stress_elmt(:,i,ielmt)+sigma
!    !if(i_elmt==1.and.i==1)then
!    !      print*,'s',sigma
!    !      print*,'e',eld
!    !      !print*,'ev',evpt
!    !     stop
!    !    endif
!  end do ! i GLL
!end do ! i_elmt
!!print*,size(stress_elmt)
!return
!end subroutine elastic_stress_intact
!!===========================================
end module elastic
!===============================================================================
