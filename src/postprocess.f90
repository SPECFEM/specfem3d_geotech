! this module contains post processing library routines
! AUTHOR
!   Hom Nath Gharti
module postprocess
use set_precision

contains

!-------------------------------------------------------------------------------
subroutine overburden_stress(nelmt,g_num,mat_id,z0,p0,k0,stress_elmt)
use global,only:nenod,ngll,nst,g_coord,gam_blk
implicit none
integer,intent(in) :: nelmt
integer,intent(in) :: g_num(nenod,nelmt),mat_id(nelmt)
real(kind=kreal),intent(in) :: z0,p0,k0
real(kind=kreal),intent(inout) :: stress_elmt(nst,ngll,nelmt)

real(kind=kreal) :: z,szz
integer :: i,i_elmt,num(nenod)

do i_elmt=1,nelmt
  num=g_num(:,i_elmt)

  do i=1,ngll ! loop over integration points
    z=g_coord(3,num(i))
    szz=p0-gam_blk(mat_id(i_elmt))*abs(z-z0) ! compression -
    stress_elmt(3,i,i_elmt)=szz
    stress_elmt(1,i,i_elmt)=k0*szz
    stress_elmt(2,i,i_elmt)=k0*szz
  end do ! i GLL
end do ! i_elmt
return
end subroutine overburden_stress
!===============================================================================

! TODO: is it possible to compute stress_elmt only for intact elements just in
! case?
! it seems that the subarray cannot be a receiving array
! this subroutine computes elastic stress from the known displacement
subroutine elastic_stress(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,dshape_hex8,   &
dlagrange_gll,x,stress_elmt)
use global,only:ndim,nedof,nenod,ngnode,ngll,nst,g_coord,ym_blk,nu_blk
use elastic,only:compute_cmat
use preprocess,only:compute_bmat
use math_library,only:determinant,invert
implicit none
integer,intent(in) :: nelmt,neq,gnod(8)
integer,intent(in) :: g_num(nenod,nelmt),gdof_elmt(nedof,nelmt),mat_id(nelmt)
real(kind=kreal),intent(in) :: dshape_hex8(ndim,ngnode,ngll),                   &
dlagrange_gll(ndim,ngll,ngll),x(0:neq)
real(kind=kreal),intent(inout) :: stress_elmt(nst,ngll,nelmt)

real(kind=kreal) :: detjac
real(kind=kreal) :: cmat(nst,nst),coord(ngnode,ndim),jac(ndim,ndim),            &
deriv(ndim,nenod),bmat(nst,nedof),eld(nedof),eps(nst),sigma(nst)
integer :: egdof(nedof),num(nenod)
integer :: i,i_elmt

do i_elmt=1,nelmt
  call compute_cmat(cmat,ym_blk(mat_id(i_elmt)),nu_blk(mat_id(i_elmt)))
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  eld=x(egdof)

  do i=1,ngll ! loop over integration points
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)!

    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat(deriv,bmat)
    eps=matmul(bmat,eld)
    sigma=matmul(cmat,eps)
    stress_elmt(:,i,i_elmt)=sigma
  end do ! i GLL
end do ! i_elmt
return
end subroutine elastic_stress
!===============================================================================

subroutine elastic_stress_intact(nelmt_intact,neq,gnod,elmt_intact,g_num,      &
gdof_elmt,mat_id,dshape_hex8,dlagrange_gll,x,stress_intact)
use global,only:ndim,nedof,nenod,ngnode,ngll,nst,g_coord,ym_blk,nu_blk
use elastic,only:compute_cmat
use preprocess,only:compute_bmat
use math_library,only:determinant,invert
implicit none
integer,intent(in) :: nelmt_intact,neq,gnod(8)
!integer,intent(in) :: elmt_intact(nelmt_intact),g_num(nenod,nelmt_intact), &
!gdof_elmt(nedof,nelmt_intact),mat_id(nelmt_intact)
integer,intent(in) :: elmt_intact(:),g_num(:,:),gdof_elmt(:,:),mat_id(:)
real(kind=kreal),intent(in) :: dshape_hex8(ndim,ngnode,ngll),                  &
dlagrange_gll(ndim,ngll,ngll),x(0:neq)
real(kind=kreal),intent(inout) :: stress_intact(nst,ngll,nelmt_intact)

real(kind=kreal) :: detjac
real(kind=kreal) :: cmat(nst,nst),coord(ngnode,ndim),jac(ndim,ndim),           &
deriv(ndim,nenod),bmat(nst,nedof),eld(nedof),eps(nst),sigma(nst)
integer :: egdof(nedof),num(nenod)
integer :: i,i_elmt,ielmt

do i_elmt=1,nelmt_intact
  ielmt=elmt_intact(i_elmt)
  call compute_cmat(cmat,ym_blk(mat_id(i_elmt)),nu_blk(mat_id(i_elmt)))
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  eld=x(egdof)
  do i=1,ngll ! loop over integration points
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)!

    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat(deriv,bmat)
    eps=matmul(bmat,eld)
    sigma=matmul(cmat,eps)
    stress_intact(:,i,ielmt)=stress_intact(:,i,ielmt)+sigma
  end do ! i GLL
end do ! i_elmt
return
end subroutine elastic_stress_intact
!===============================================================================

! this routine save data to files Ensight Gold format
! TODO: make it optional
subroutine save_data(ptail,format_str,istep,nnode,nodalu,scf,vmeps,            &
stress_global)
use global,only:nst,out_path,file_head,savedata
use math_constants
use math_library,only:rquick_sort,stress_invariant
use visual
implicit none
character(len=20),intent(in) :: format_str,ptail
integer,intent(in) :: istep,nnode
real(kind=kreal),intent(in) :: nodalu(3,nnode),scf(nnode),vmeps(nnode),        &
stress_global(nst,nnode)

integer :: i_node,npart
real(kind=kreal) :: psigma(3,nnode),nsigma(nnode),taumax(nnode)
real(kind=kreal) :: dsbar,lode_theta,sigm,two_third=2.0_kreal/3.0_kreal
character(len=250) :: out_fname
character(len=80) :: destag ! this must be 80 characters long

if(savedata%psigma.or.savedata%nsigma.or.savedata%maxtau)then
  ! compute principal stresses
  psigma=zero
  do i_node=1,nnode
    call stress_invariant(stress_global(:,i_node),sigm,dsbar,lode_theta)
    psigma(:,i_node)=sigm+two_third*dsbar* &
          sin((/ lode_theta-two_third*pi,lode_theta,lode_theta+two_third*pi /))
    ! put in ascending order (compresson (-) as a major principal stress!)
    psigma(:,i_node)=rquick_sort(psigma(:,i_node),3)
  enddo
endif
! store initial principal stresses
!if(istep==0)psigma0=psigma
! compute major principal stress factor
!scf=one
!if(istep>0)then
!  where(psigma0(1,:)==zero)
!    scf=inftol
!  elsewhere
!      scf=psigma(1,:)/psigma0(1,:)
!  end where
!endif

if(savedata%maxtau)then
  ! compute maximum shear stress
  taumax=zero
  do i_node=1,nnode
    taumax(i_node)=half*abs(psigma(1,i_node)-psigma(3,i_node))
  enddo
endif

if(savedata%nsigma)then
  ! compute normal stress
  nsigma=zero
  do i_node=1,nnode
    nsigma(i_node)=half*(psigma(1,i_node)+psigma(3,i_node))
  enddo
endif

if(savedata%disp)then
  ! write displacement vector
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',     &
  istep,trim(ptail)//'.dis'
  npart=1;
  destag='Displacement field'
  call write_ensight_pernodeVECAS(out_fname,destag,npart,3,nnode,real(nodalu))
endif

if(savedata%stress)then
  ! write stress tensor
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',     &
  istep,trim(ptail)//'.sig'
  npart=1;
  destag='Effective stress tensor'
  call write_ensight_pernodeVECAS(out_fname,destag,npart,6,nnode,              &
  real(stress_global))
endif

if(savedata%psigma)then
  ! write principal stress
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',     &
  istep,trim(ptail)//'.psig'
  npart=1;
  destag='principal stresses'
  call write_ensight_pernodeVECAS(out_fname,destag,npart,3,nnode,real(psigma))
endif

if(savedata%scf)then
  ! write stress concentration factor
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',     &
  istep,trim(ptail)//'.scf'
  npart=1;
  destag='principal stress concentration factor'
  call write_ensight_pernodeSCALAS(out_fname,destag,npart,nnode,real(scf))
endif

if(savedata%maxtau)then
  ! write maximum shear stress
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',     &
  istep,trim(ptail)//'.mtau'
  npart=1;
  destag='maximum shear stress'
  call write_ensight_pernodeSCALAS(out_fname,destag,npart,nnode,real(taumax))
endif

if(savedata%nsigma)then
  ! write normal stress on maximum-shear-stress plane
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',     &
  istep,trim(ptail)//'.nsig'
  npart=1;
  destag='normal stress'
  call write_ensight_pernodeSCALAS(out_fname,destag,npart,nnode,real(nsigma))
endif

if(savedata%vmeps)then
  ! write accumulated effective plastic strain
  write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',     &
  istep,trim(ptail)//'.eps'
  npart=1;
  destag='von Mises effective plastic strain'
  call write_ensight_pernodeSCALAS(out_fname,destag,npart,nnode,real(vmeps))
endif
end subroutine save_data
!===============================================================================
end module postprocess
!===============================================================================
