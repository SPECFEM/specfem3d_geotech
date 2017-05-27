! include 'license.txt'
! this is a main routine for multistage excavation
! this program was originally based on the book "Programming the finite element
! method" Smith and Griffiths (2004)
! REVISION:
!   HNG, Aug 25,2011; HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
subroutine semexcav3d(ismpi,myid,gnod,sum_file,ptail,format_str)
! import necessary libraries
use global
use string_library, only : parse_file
use math_constants
use gll_library
!use mesh_spec
use shape_library
use math_library
use preprocess
!use gauss_library
use excavation
use plastic_library
#if (USE_MPI)
use mpi_library
use ghost_library_mpi
use math_library_mpi
use solver_mpi
#else
use serial_library
use math_library_serial
use solver
#endif
use visual
use postprocess

implicit none
logical,intent(in) :: ismpi
integer,intent(in) :: myid
integer,intent(in) :: gnod(8)
character(len=250),intent(in) :: sum_file
character(len=20),intent(in) :: ptail,format_str

integer :: funit,i,ios,istat,j,k,neq
integer :: i_elmt,i_node,i_inc,i_srf,i_excav,ielmt,imat,inode
real(kind=kreal),parameter :: r3=3.0_kreal,two_third=two/r3
real(kind=kreal) :: detjac,dq1,dq2,dq3,dsbar,dt,f,fmax,lode_theta,phifr, &
sf,sigm

real(kind=kreal) :: uerr,umax,uxmax
integer :: cg_iter,cg_tot,nl_iter,nl_tot
logical :: nl_isconv ! logical variable to check convergence of
! nonlinear (NL) iterations

real(kind=kreal) :: cmat(nst,nst),devp(nst),eps(nst),erate(nst),evp(nst),      &
flow(nst,nst),m1(nst,nst),m2(nst,nst),m3(nst,nst),effsigma(nst),sigma(nst)
! dynamic arrays
integer,allocatable::gdof(:,:),gdof_elmt(:,:),num(:),node_valency(:)
! factored parameters
real(kind=kreal),allocatable :: cohf(:),nuf(:),phif(:),psif(:),ymf(:)
real(kind=kreal),allocatable::bodyload(:),bmat(:,:),bload(:),coord(:,:),      &
der(:,:),deriv(:,:),dprecon(:),eld(:),eload(:),evpt(:,:,:),excavload(:,:),    &
extload(:),jac(:,:),load(:),nodalu(:,:),storkm(:,:,:),oldx(:),x(:),           &
stress_local(:,:,:),stress_global(:,:),scf(:),vmeps(:)
!,psigma(:,:),psigma0(:,:),taumax(:),nsigma(:)
integer,allocatable :: egdof(:) ! elemental global degree of freedom

integer :: map2exodus(8),ngllxy,node_hex8(8)
real(kind=kreal),allocatable :: dshape_hex8(:,:,:)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal
!double precision
real(kind=kreal),allocatable :: xigll(:),wxgll(:) !double precision
real(kind=kreal),allocatable :: etagll(:),wygll(:) !double precision
real(kind=kreal),allocatable :: zetagll(:),wzgll(:) !double precision
real(kind=kreal),allocatable :: gll_weights(:),gll_points(:,:)
real(kind=kreal),allocatable :: lagrange_gll(:,:),dlagrange_gll(:,:,:)

character(len=250) :: out_fname
character(len=20), parameter :: wild_char='********************'
character(len=20) :: ensight_etype
character(len=80) :: buffer,destag ! this must be 80 characters long
character(len=250) :: geo_file !,sum_file
integer :: npart

logical :: gravity,pseudoeq ! gravity load and pseudostatic load
real(kind=kreal),allocatable :: wpressure(:) ! water pressure
logical,allocatable :: submerged_node(:)

! excavation
integer :: id0,id1
integer :: nelmt_intact,nelmt_void
integer,allocatable :: elmt_intact(:),elmt_void(:)
integer :: nnode_intact,nnode_void
integer,allocatable :: nmir(:),node_intact(:),node_void(:)
logical,allocatable :: ismat(:),isnode(:)

integer :: ipart !,myid,nproc
integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode
integer :: tot_neq,max_neq,min_neq
! number of active ghost partitions for a node
integer,allocatable :: ngpart_node(:)
character(len=250) :: errtag ! error message
integer :: errcode

errtag=""; errcode=-1

ipart=myid-1 ! partition id starts from 0

allocate(ismat(nmat))
ismat=.true.

ngllxy=ngllx*nglly
ensight_etype='hexa8'

! map sequential node numbering to exodus/cubit order for 8-noded hexahedra
map2exodus=(/ 1,2,4,3,5,6,8,7 /)

! apply displacement boundary conditions
if(myid==1)write(stdout,'(a)',advance='no')'applying BC...'
allocate(gdof(nndof,nnode),stat=istat)
if (istat/=0)then
  write(stdout,*)'ERROR: cannot allocate memory!'
  stop
endif
gdof=1
call apply_bc(ismpi,myid,nproc,gdof,neq,errcode,errtag)
if(errcode/=0)call error_stop(errtag,stdout,myid)
if(myid==1)write(stdout,*)'complete!'
!-------------------------------------

allocate(isnode(nnode),num(nenod),evpt(nst,ngll,nelmt),coord(ngnod,ndim),      &
jac(ndim,ndim),der(ndim,ngnod),deriv(ndim,nenod),bmat(nst,nedof),eld(nedof),   &
bload(nedof),eload(nedof),nodalu(nndof,nnode),egdof(nedof),stat=istat)
if (istat/=0)then
  write(stdout,*)'ERROR: cannot allocate memory!'
  stop
endif

tot_neq=sumscal(neq); max_neq=maxscal(neq); min_neq=minscal(neq)
if(myid==1)then
  write(stdout,*)'degrees of freedoms => total:',tot_neq,' max:',max_neq,      &
  ' min:',min_neq
endif

! get gll points and weights
allocate(xigll(ngllx),wxgll(ngllx),etagll(nglly),wygll(nglly),zetagll(ngllz),  &
wzgll(ngllz))
call zwgljd(xigll,wxgll,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(etagll,wygll,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(zetagll,wzgll,ngllz,jacobi_alpha,jacobi_beta)

! get derivatives of shape functions for 8-noded hex
allocate(dshape_hex8(ndim,ngnod,ngll))
call dshape_function_hex8(ngnod,ngllx,nglly,ngllz,xigll,etagll,zetagll,   &
dshape_hex8)
deallocate(xigll,wxgll,etagll,wygll,zetagll,wzgll)
! compute gauss-lobatto-legendre quadrature information
allocate(gll_weights(ngll),gll_points(ndim,ngll),lagrange_gll(ngll,ngll),      &
dlagrange_gll(ndim,ngll,ngll))
call gll_quadrature(ndim,ngllx,nglly,ngllz,ngll,gll_points,gll_weights,        &
lagrange_gll,dlagrange_gll)
!--------------------------------

! store elemental global degrees of freedoms from nodal gdof
! this removes the repeated use of reshape later but it has larger size than gdof!!!
allocate(gdof_elmt(nedof,nelmt))
gdof_elmt=0
do i_elmt=1,nelmt
  gdof_elmt(:,i_elmt)=reshape(gdof(:,g_num(:,i_elmt)),(/nedof/)) !g=g_g(:,i_elmt)
enddo
!-------------------------------

if(myid==1)write(stdout,'(a)',advance='no')'preprocessing...'

! initalize intact elements and nodes
nnode_intact=nnode; nelmt_intact=nelmt
nnode_void=0; nelmt_void=0
allocate(node_intact(nnode_intact),elmt_intact(nelmt_intact))
node_intact=(/ (i,i=1,nnode) /)
elmt_intact=(/ (i,i=1,nelmt) /)

allocate(stress_local(nst,ngll,nelmt))
! compute initial stress assuming elastic domain
stress_local=zero

if(s0_type==0)then
  ! compute initial stress using SEM itself

  allocate(extload(0:neq),x(0:neq),dprecon(0:neq),storkm(nedof,nedof,          &
  nelmt_intact),stat=istat) ! elastic(0:neq),
  if (istat/=0)then
    write(stdout,*)'ERROR: cannot allocate memory!'
    stop
  endif
  extload=zero; gravity=.true.; pseudoeq=.false.
  call stiffness_bodyload(nelmt_intact,neq,gnod,g_num(:,elmt_intact),          &
  gdof_elmt(:,elmt_intact),mat_id(elmt_intact),gam,nu,ym,dshape_hex8,          &
  dlagrange_gll,gll_weights,storkm,dprecon,extload,gravity,       &
  pseudoeq)

  !print*,minval(dprecon),maxval(dprecon)
  !print*,minval(extload),maxval(extload)
  !print*,minval(storkm),maxval(storkm)
  !stop
  if(myid==1)write(stdout,*)'complete!'
  !-------------------------------

  ! apply traction boundary conditions
  if(istraction)then
    if(myid==1)write(*,'(a)',advance='no')'applying traction...'
    call apply_traction(ismpi,myid,nproc,gnod,gdof,neq,extload,errcode,errtag)
    if(errcode/=0)call error_stop(errtag,stdout,myid)
    if(myid==1)write(*,*)'complete!'
  endif
  !-------------------------------

  ! compute water pressure
  if(iswater)then
    if(myid==1)write(stdout,'(a)',advance='no')'computing water pressure...'
    allocate(wpressure(nnode),submerged_node(nnode))
    call compute_pressure(wpressure,submerged_node,errcode,   &
    errtag)
    if(errcode/=0)call error_stop(errtag,stdout,myid)
    ! write pore pressure file

    ! open Ensight Gold data file to store data
    out_fname=trim(out_path)//trim(file_head)//trim(ptail)//'.por'
    npart=1;
    destag='Pore pressure'
    call write_ensight_pernode(out_fname,destag,npart,1,nnode,real(wpressure))
    if(myid==1)write(stdout,*)'complete!'
  endif
  !-------------------------------

  if(myid==1)write(stdout,'(a)')'--------------------------------------------'

  ! prepare ghost partitions for the communication
  call prepare_ghost(myid,nproc,gdof)

  ! assemble from ghost partitions
  call assemble_ghosts(nndof,neq,dprecon,dprecon)
  !print*,minval(dprecon),maxval(dprecon)
  !print*,minval(storkm),maxval(storkm)
  !stop
  dprecon(1:)=one/dprecon(1:); dprecon(0)=zero

  ! compute displacement due to graviy loading to compute initial stress
  x=zero
  call pcg_solver(neq,nelmt_intact,storkm,x,extload,     &
  dprecon,gdof_elmt(:,elmt_intact),cg_iter,errcode,errtag)
  if(errcode/=0)call error_stop(errtag,stdout,myid)
  x(0)=zero

  call elastic_stress(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,dshape_hex8,       &
  dlagrange_gll,x,stress_local)
  deallocate(extload,dprecon,x,storkm)
elseif(s0_type==1)then
  ! compute initial stress using simple relation for overburden pressure
  call overburden_stress(nelmt,g_num,mat_id,z_datum,s0_datum,epk0,stress_local)
else
  write(stdout,*)'ERROR: s0_type:',s0_type,' not supported!'
  stop
endif
!-------------------------------

allocate(stress_global(nst,nnode),vmeps(nnode))
allocate(scf(nnode)) !,psigma(ndim,nnode),psigma0(ndim,nnode),taumax(nnode),nsigma(nnode))
allocate(nmir(nnode),node_valency(nnode))

! compute node valency only once
node_valency=0
do i_elmt=1,nelmt_intact
  ielmt=elmt_intact(i_elmt)
  num=g_num(:,ielmt)
  node_valency(num)=node_valency(num)+1
enddo

! open summary file
!sum_file = trim(out_path)//trim(file_head)//'_summary'//trim(ptail)
open(unit=10,file=trim(sum_file),status='old',position='append',action='write',&
iostat=ios)
write(10,*)'CG_MAXITER, CG_TOL, NL_MAXITER, NL_TOL'
write(10,*)cg_maxiter,cg_tol,nl_maxiter,nl_tol
write(10,*)'Number of SRFs'
write(10,*)nsrf

allocate(cohf(nmat),nuf(nmat),phif(nmat),psif(nmat),ymf(nmat))

if(myid==1)then
  write(stdout,'(a,e12.4,a,i5)')'CG_TOL:',cg_tol,' CG_MAXITER:',cg_maxiter
  write(stdout,'(a,e12.4,a,i5)')'NL_TOL:',nl_tol,' NL_MAXITER:',nl_maxiter
  write(stdout,'(a)',advance='no')'SRFs:'
  do i_srf=1,nsrf
    write(stdout,'(f7.4,1x)',advance='no')srf(i_srf)
  enddo
  write(stdout,*)
  write(stdout,'(a,i3)')'nexcav:',nexcav
  write(stdout,'(/,a)')'--------------------------------------------'
endif

srf_loop: do i_srf=1,nsrf
write(10,*)'SRF'
write(10,*)srf(i_srf)
write(10,*)'Number of excavation stages'
write(10,*)nexcav
write(10,*)'STEP, CGITER, NLITER, UXMAX, UMAX, fmax'
close(10)

allocate(excavload(nndof,nnode))
allocate(ngpart_node(nnode))

! excavation-stage loop
nodalu=zero
excavation_stage: do i_excav=0,nexcav

  vmeps=zero; scf=inftol
  if(i_excav>0)then
  !deallocate(load,bodyload,oldx,extload,x,dprecon,stat=istat)
  if(myid==1)write(stdout,'(/,a,i4)')'excavation stage:',i_excav

  ! find appropriate indices in excavid
  id0=(i_excav-1)*nexcavid(i_excav)+1; id1=id0+nexcavid(i_excav)-1

  ! disable excavated material id
  ismat(excavid(id0:id1))=off !ismat(excavid(i_excav))=off
  !print*,id0,id1
  !stop
  ! count intact and void elements after excavation
  nelmt_void=0
  do i=id0,id1
    nelmt_void=nelmt_void+count(mat_id==excavid(i))
  enddo
  !nelmt_void=count(mat_id==excavid(i_excav))
  nelmt_intact=nelmt_intact-nelmt_void
  if(nelmt_intact<=0)then
    write(stdout,'(/,a)')'WARNING: no intact elements left!'
    exit excavation_stage
  endif
  allocate(elmt_intact(nelmt_intact),elmt_void(nelmt_void))
  ! find intact and void elements after excavation
  call intact_void_elmt(nexcavid(i_excav),excavid(id0:id1),ismat,nelmt_intact, &
  nelmt_void,elmt_intact,elmt_void,isnode)

  ! count intact and void nodes after excavation
  nnode_intact=count(isnode)
  nnode_void=nnode-nnode_intact
  allocate(node_intact(nnode_intact),node_void(nnode_void))
  ! find intact and void nodes after excavation
  call intact_void_node(isnode,nnode_intact,nnode_void,node_intact,node_void,  &
  nmir)

  tot_nelmt=sumscal(nelmt_intact); tot_nnode=sumscal(nnode_intact)
  max_nelmt=maxscal(nelmt_intact); max_nnode=maxscal(nnode_intact)
  min_nelmt=minscal(nelmt_intact); min_nnode=minscal(nnode_intact)
  if(myid==1)then
    write(stdout,*)'intact elements => total:',tot_nelmt,' max:',max_nelmt,' min:',min_nelmt
    write(stdout,*)'intact nodes    => total:',tot_nnode,' max:',max_nnode,' min:',min_nnode
  endif

  !write(stdout,'(a,i10)')' total intact nodes:',nnode_intact
  !write(stdout,'(a,i10)')' total intact elements:',nelmt_intact

  ! correct node valency subtracting dead element-nodes
  call correct_nvalency(node_valency,nelmt_void,g_num(:,elmt_void))

  ! modify global-degrees-of-freedom (gdof)
  call modify_gdof(gdof,nnode_void,node_void,neq)

  tot_neq=sumscal(neq); max_neq=maxscal(neq); min_neq=minscal(neq)
  if(myid==1)then
    write(stdout,*)'degrees of freedoms => total:',tot_neq,' max:',max_neq,' min:',min_neq
  endif
  !write(stdout,'(a,i10)')' total degrees of freedoms:',neq

  ! store elemental global degrees of freedoms from nodal gdof
  ! this removes the repeated use of reshape later but it has larger size than gdof!!!
  gdof_elmt=0
  do i_elmt=1,nelmt
    gdof_elmt(:,i_elmt)=reshape(gdof(:,g_num(:,i_elmt)),(/nedof/)) !g=g_g(:,i_elmt)
  enddo

  ! write geo file for this stage
  ! open Ensight Gold geo file to store mesh data
  write(geo_file,fmt=format_str)trim(out_path)//trim(file_head)//'_step',i_excav,trim(ptail)//'.geo'
  npart=1
  destag='unstructured meshes'
  call write_ensight_geocoord(geo_file,destag,npart,nnode_intact,real(g_coord(:,node_intact)),funit)

  ! writes element information
  buffer=ensight_etype
  write(funit)buffer
  write(funit)nelmt_intact*(ngllx-1)*(nglly-1)*(ngllz-1)

  ! do not substract 1 for ensight file
  do i_elmt=1,nelmt_intact
    ielmt=elmt_intact(i_elmt)
    do k=1,ngllz-1
      do j=1,nglly-1
        do i=1,ngllx-1
          ! corner nodes in a sequential numbering
          node_hex8(1)=(k-1)*ngllxy+(j-1)*ngllx+i
          node_hex8(2)=node_hex8(1)+1

          node_hex8(3)=node_hex8(1)+ngllx
          node_hex8(4)=node_hex8(3)+1

          node_hex8(5)=node_hex8(1)+ngllxy
          node_hex8(6)=node_hex8(5)+1

          node_hex8(7)=node_hex8(5)+ngllx
          node_hex8(8)=node_hex8(7)+1
          ! map to exodus/cubit numbering and write
          write(funit)nmir(g_num(node_hex8(map2exodus),ielmt))
        enddo
      enddo
    enddo
  enddo
  close(funit)

  ! reallocate those arrays whose size depend on the neq
  allocate(load(0:neq),bodyload(0:neq),extload(0:neq),oldx(0:neq),x(0:neq), & !,extload1(0:neq)
  dprecon(0:neq),storkm(nedof,nedof,nelmt_intact),stat=istat) ! elastic(0:neq),
  if (istat/=0)then
    write(stdout,*)'ERROR: cannot allocate memory!'
    stop
  endif

  ! modify ghost partitions after excavation
  !call prepare_ghost(myid,nproc,gdof)
  call modify_ghost(myid,gdof,isnode)
  call count_active_nghosts(ngpart_node)

  excavload=zero; extload=zero; ! extload1=zero

  ! compute excavation load at gdofs
  !call excavation_load(nelmt_void,neq,gnod,g_num(:,elmt_void),gdof_elmt(:,elmt_void), &
  !mat_id(elmt_void),dshape_hex8,dlagrange_gll,gll_weights, &
  !stress_local(:,:,elmt_void),extload)

  ! compute excavation load at nodes
  call excavation_load_nodal(nelmt_void,gnod,g_num(:,elmt_void), &
  mat_id(elmt_void),dshape_hex8,dlagrange_gll,gll_weights, &
  stress_local(:,:,elmt_void),excavload)

  ! if the excavation load is discarded by the partition (it can happens due to
  ! the special combination of partition geometry and excavation geoemtry) it
  ! should be distributed equally to the active sharing partitions.
  call distribute2ghosts(gdof,nndof,neq,ngpart_node,excavload,extload)
!#else
!  ! store nodal values to gdof locations
!  do j=1,nnode
!    do i=1,nndof
!      igdof=gdof(i,j)
!      extload(igdof)=extload(igdof)+excavload(i,j)
!    enddo
!  enddo
!#endif
  !extload=extload1
  extload(0)=zero; !extload1(0)=zero

  ! strength reduction
  call strength_reduction(srf(i_srf),phinu,nmat,coh,nu,phi,psi,cohf,nuf,phif,  &
  psif,istat)

  ! compute minimum pseudo-time step for viscoplasticity
  dt=dt_viscoplas(nmat,nuf,phif,ym,ismat)

  ! compute stiffness matrix
  gravity=.false.; pseudoeq=.false.
  call stiffness_bodyload(nelmt_intact,neq,gnod,g_num(:,elmt_intact),          &
  gdof_elmt(:,elmt_intact),mat_id(elmt_intact),gam,nuf,ym,dshape_hex8,         &
  dlagrange_gll,gll_weights,storkm,dprecon)!,extload,gravity,pseudoeq)

  ! assemble from ghost partitions
  call assemble_ghosts(nndof,neq,dprecon,dprecon)
  dprecon(0)=zero; dprecon(1:)=one/dprecon(1:)

  ! find global dt
  dt=minscal(dt)

  cg_tot=0; nl_tot=0
  ! load incremental loop
  if(myid==1)write(stdout,'(a,i10)')' total load increments:',ninc
  extload=extload/ninc
  load_increment: do i_inc=1,ninc
  !stress_local(:,:,elmt_void)=zero

  bodyload=zero; evpt=zero

  x=zero; oldx=zero
  ! plastic iteration loop
  plastic: do nl_iter=1,nl_maxiter
    fmax=zero

    load=extload+bodyload
    load(0)=zero

    ! pcg solver
    !x=zero
    call pcg_solver(neq,nelmt_intact,storkm,x,load,      &
    dprecon,gdof_elmt(:,elmt_intact),cg_iter,errcode,errtag)
    if(errcode/=0)call error_stop(errtag,stdout,myid)
    cg_tot=cg_tot+cg_iter
    x(0)=zero

    if(allelastic)then
      !print*,size(stress_local)
      call elastic_stress_intact(nelmt_intact,neq,gnod,elmt_intact,            &
      g_num(:,elmt_intact),gdof_elmt(:,elmt_intact),mat_id(elmt_intact),       &
      dshape_hex8,dlagrange_gll,x,stress_local(:,:,:))

      exit plastic
    endif

    ! check plastic convergence
    uerr=maxvec(abs(x-oldx))/maxvec(abs(x))
    oldx=x
    nl_isconv=uerr.le.nl_tol

    ! compute stress and check failure
    do i_elmt=1,nelmt_intact
      ielmt=elmt_intact(i_elmt)

      imat=mat_id(ielmt)
      !tnph=tan(phi(mat_id(ielmt))*pi/r180)
      !phif=atan(tnph/srf(i_srf))*r180/pi
      !tnps=tan(psi(mat_id(ielmt))*pi/r180)
      !psif=atan(tnps/srf(i_srf))*r180/pi
      !cf=coh(mat_id(ielmt))/srf(i_srf)
      call compute_cmat(cmat,ym(imat),nuf(imat))
      num=g_num(:,ielmt)
      coord=transpose(g_coord(:,num(gnod))) !transpose(g_coord(:,num(1:ngnod)))
      egdof=gdof_elmt(:,ielmt) !reshape(gdof(:,g_num(:,ielmt)),(/nedof/))
      eld=x(egdof)
      !print*,egdof
      !stop
      bload=zero
      do i=1,ngll ! loop over integration points
        jac=matmul(dshape_hex8(:,:,i),coord)
        detjac=determinant(jac)
        call invert(jac)

        deriv=matmul(jac,dlagrange_gll(:,i,:))
        call compute_bmat(bmat,deriv)
        eps=matmul(bmat,eld)
        eps=eps-evpt(:,i,ielmt)
        sigma=matmul(cmat,eps)

        ! compute effective pressure
        effsigma=sigma
        if(iswater)then
          if(submerged_node(num(i)))then
            ! water pressure is compressive (negative)
            effsigma(1:3)=effsigma(1:3)+wpressure(num(i))
          endif
        endif

        effsigma=effsigma+stress_local(:,i,ielmt)
        call stress_invariant(effsigma,sigm,dsbar,lode_theta)
        ! check whether yield is violated
        call mohcouf(phif(imat),cohf(imat),sigm,dsbar,lode_theta,f)
        if(f>fmax)fmax=f

        if(f>=zero)then !.or.(nl_isconv.or.nl_iter==nl_maxiter))then
          call mohcouq(psif(imat),dsbar,lode_theta,dq1,dq2,dq3)
          call formm(effsigma,m1,m2,m3)
          !if(dsbar<=zerotol)print*,m1*dq1+m2*dq2+m3*dq3
          flow=f*(m1*dq1+m2*dq2+m3*dq3)

          erate=matmul(flow,effsigma)
          evp=erate*dt
          evpt(:,i,ielmt)=evpt(:,i,ielmt)+evp
          devp=matmul(cmat,evp)
          ! if not converged we need body load for next iteration
          if(.not.nl_isconv .and. nl_iter/=nl_maxiter)then
            !devp(1:3)=devp(1:3)-wpressure(num(i))
            eload=matmul(devp,bmat)
            bload=bload+eload*detjac*gll_weights(i)
          end if
        end if
        if(nl_isconv.or.nl_iter==nl_maxiter)then
          devp=sigma
          ! compute von Mises effective plastic strain
          vmeps(num(i))=vmeps(num(i))+sqrt(two_third*dot_product(evpt(:,i,ielmt),evpt(:,i,ielmt)))
          ! update stresses
          stress_local(:,i,ielmt)=effsigma
          phifr=phif(imat)*deg2rad !atan(tnph/srf(i_srf))
          sf=(sigm*sin(phifr)-cohf(imat)*cos(phifr))/(-dsbar*(cos(lode_theta)/ &
          sqrt(r3)-sin(lode_theta)*sin(phifr)/r3))
          if(sf<scf(num(i)))scf(num(i))=sf
        endif
        !if(i_elmt==1.and.i==1)then
        !  print*,'s',sigma
        !  print*,'e',eld
        !  print*,'ev',maxval(abs(evpt))
        !  stop
        !endif
      end do ! i_gll

      if(nl_isconv .or. nl_iter==nl_maxiter)cycle
      ! compute total body load vector
      bodyload(egdof)=bodyload(egdof)+bload
    end do ! i_elmt
    bodyload(0)=zero
    fmax=maxscal(fmax)
    uxmax=maxvec(abs(x))
    if(myid==1)then
      write(stdout,'(a,a,i4,a,i4,a,f12.6,a,f12.6,a,f12.6)',advance='no')CR,    &
      ' ninc:',i_inc,' nl_iter:',nl_iter,' f_max:',fmax,' uerr:',uerr,' umax:',&
      uxmax
    endif

    if(nl_isconv.or.nl_iter==nl_maxiter)exit
  end do plastic ! plastic iteration
  if(nl_iter>=nl_maxiter .and. .not.nl_isconv)then
    write(stdout,*)'WARNING: nonconvergence in nonlinear iterations!'
    write(stdout,*)'desired tolerance:',nl_tol,' achieved tolerance:',uerr
  endif
  nl_tot=nl_tot+nl_iter
  !if(myid==1)print*,cg_tot,nl_tot
  ! nodal displacement
  do i=1,nndof
    do j=1,nnode
      if(gdof(i,j)/=0)then
        nodalu(i,j)=nodalu(i,j)+x(gdof(i,j)) !-elastic(gdof(i,j))
      endif
    enddo
  enddo
  enddo load_increment ! load increment loop
  deallocate(load,bodyload,extload,oldx,x,dprecon,storkm,stat=istat)
  if (istat/=0)then
    write(stdout,*)'ERROR: cannot deallocate memory!'
    stop
  endif

  ! write summary
  uxmax=maxvec(abs(reshape(nodalu(:,node_intact),(/nndof*nnode_intact/))))
  umax=maxvec(sqrt(nodalu(1,node_intact)*nodalu(1,node_intact)+ &
  nodalu(2,node_intact)*nodalu(2,node_intact)+nodalu(3,node_intact)*           &
  nodalu(3,node_intact)))
  open(10,file=trim(sum_file),status='old',position='append',action='write')
  write(10,*)i_excav,cg_tot,nl_tot,uxmax,umax,fmax
  close(10)

  ! compute average effective strain
  do i_node=1,nnode_intact
    inode=node_intact(i_node)
    vmeps(inode)=vmeps(inode)/real(node_valency(inode),kreal)
  enddo

  endif ! if(i_excav>0)

  ! compute stress_global
  stress_global=zero
  do i_elmt=1,nelmt_intact
    ielmt=elmt_intact(i_elmt)
    num=g_num(:,ielmt)
    stress_global(:,num)=stress_global(:,num)+stress_local(:,:,ielmt)
  enddo

  ! compute average stress at sharing nodes
  do i_node=1,nnode_intact
    inode=node_intact(i_node)
    stress_global(:,inode)=stress_global(:,inode)/real(node_valency(inode),    &
    kreal)
  enddo

  call save_data(ptail,format_str,i_excav,nnode_intact,           &
  nodalu(:,node_intact),scf(node_intact),                 &
  vmeps(node_intact),stress_global(:,node_intact))

  ! deallocate those variables whose size depend on changing geometry
  deallocate(elmt_intact,node_intact,stat=istat)
  if(i_excav>0)deallocate(elmt_void,node_void)

  !call sync_process
  if(nl_iter==nl_maxiter)exit

enddo excavation_stage ! i_excav time stepping loop
enddo srf_loop ! i_srf safety factor loop
deallocate(mat_id,gam,ym,coh,nu,phi,psi,srf)
deallocate(excavload,g_coord,g_num,isnode,nmir)
call free_ghost()
!-----------------------------------

return
end subroutine semexcav3d
!===========================================


