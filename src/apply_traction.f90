! TODO: make use of the orthogonality
! this routine read and applies the traction specified in the traction file
! AUTHOR
!   Hom Nath Gharti
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
subroutine apply_traction(ismpi,gnod,neq,load,errcode,errtag)
use global
use math_constants
use shape_library,only : dshape_function_quad4
use gll_library,only : gll_quadrature2d,zwgljd
implicit none
logical,intent(in) :: ismpi
integer,intent(in) :: gnod(ngnod)
integer,intent(in) :: neq
real(kind=kreal),intent(inout) :: load(0:neq)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: ngdof(ndim)
integer :: i,j,k,i1,i2,i3,i4,i5,i6,inod,i_face,i_gll,ios,iaxis
integer :: ielmt,iface,nface,inode
integer :: num(nenod),ngllxy,ngllyz,ngllzx,maxngll2d
integer :: nfdof,nfgll
integer,allocatable :: fgdof(:) ! face global degrees of freedom
integer :: tractype,count_trac
logical :: trac_stat

real(kind=kreal) :: coord(ndim,4)
real(kind=kreal) :: x,x1,x2,detjac
real(kind=kreal),dimension(ndim) :: face_normal,dx_dxi,dx_deta,dq_dx,q,q1,q2
real(kind=kreal),allocatable :: ftracload(:) ! face traction load

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal
real(kind=kreal),dimension(ngllx) :: xigll,wxgll !double precision
real(kind=kreal),dimension(nglly) :: etagll,wygll !double precision
real(kind=kreal),dimension(ngllz) :: zetagll,wzgll !double precision

real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4,dshape_quad4_xy, &
dshape_quad4_yz,dshape_quad4_zx
real(kind=kreal),dimension(:),allocatable :: gll_weights,gll_weights_xy,       &
gll_weights_yz,gll_weights_zx
real(kind=kreal),dimension(:,:),allocatable :: gll_points_xy,gll_points_yz,    &
gll_points_zx
real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll,lagrange_gll_xy,   &
lagrange_gll_yz,lagrange_gll_zx
real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll,                &
dlagrange_gll_xy,dlagrange_gll_yz,dlagrange_gll_zx

! face sign or normal orientation (outward +, inward -)
real(kind=kreal) :: fsign(6)
character(len=20) :: format_str,ptail
character(len=80) :: fname
character(len=80) :: data_path

type faces
  integer,allocatable :: nod(:) !ngllx*nglly) !,allocatable :: nod(:)
  integer :: gnod(4) ! corner nodes only
end type faces
type (faces) :: face(6)

errtag="ERROR: unknown!"
errcode=-1
! set data path
if(ismpi)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

if(ismpi)then
  write(format_str,*)ceiling(log10(real(nproc)+1))
  format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str)) &
  //')'
  write(ptail, fmt=format_str)'_proc',myrank
else
  ptail=""
endif

fname=trim(data_path)//trim(trfile)//trim(ptail)
open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "',trim(fname),'" cannot be opened!'
  return
endif

ngllxy=ngllx*nglly
ngllyz=nglly*ngllz
ngllzx=ngllz*ngllx
maxngll2d=max(ngllxy,ngllyz,ngllzx)

allocate(face(1)%nod(ngllzx),face(3)%nod(ngllzx))
allocate(face(2)%nod(ngllyz),face(4)%nod(ngllyz))
allocate(face(5)%nod(ngllxy),face(6)%nod(ngllxy))
! local node numbers for the faces. faces are numbered in exodus/CUBIT
! convention,but the node numbers follow the incremental indicial convention
inod=0
i1=0; i2=0; i3=0; i4=0; i5=0; i6=0
do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx
      inod=inod+1
      if (i==1)then
        ! face 4
        i4=i4+1
        face(4)%nod(i4)=inod
      endif
      if (i==ngllx)then
        ! face 2
        i2=i2+1
        face(2)%nod(i2)=inod
      endif
      if (j==1)then
        ! face 1
        i1=i1+1
        face(1)%nod(i1)=inod
      endif
      if (j==nglly)then
        ! face 3
        i3=i3+1
        face(3)%nod(i3)=inod
      endif
      if (k==1)then
        ! face 5
        i5=i5+1
        face(5)%nod(i5)=inod
      endif
      if (k==ngllz)then
        ! face 6
        i6=i6+1
        face(6)%nod(i6)=inod
      endif
    enddo
  enddo
enddo
! find corners nodes
do i_face=1,6 ! there are 6 faces in a hexahedron
  if(i_face==1 .or. i_face==3)then ! ZX plane
    face(i_face)%gnod(1)=face(i_face)%nod(1)
    face(i_face)%gnod(2)=face(i_face)%nod(ngllx)
    face(i_face)%gnod(3)=face(i_face)%nod(ngllzx)
    face(i_face)%gnod(4)=face(i_face)%nod(ngllzx-ngllx+1)
  elseif(i_face==2 .or. i_face==4)then ! YZ plane
    face(i_face)%gnod(1)=face(i_face)%nod(1)
    face(i_face)%gnod(2)=face(i_face)%nod(nglly)
    face(i_face)%gnod(3)=face(i_face)%nod(ngllyz)
    face(i_face)%gnod(4)=face(i_face)%nod(ngllyz-nglly+1)
  elseif(i_face==5 .or. i_face==6)then ! XY plane
    face(i_face)%gnod(1)=face(i_face)%nod(1)
    face(i_face)%gnod(2)=face(i_face)%nod(ngllx)
    face(i_face)%gnod(3)=face(i_face)%nod(ngllxy)
    face(i_face)%gnod(4)=face(i_face)%nod(ngllxy-ngllx+1)
  else
    write(errtag,'(a)')'ERROR: wrong face ID for traction!'
    return
  endif
enddo

! orientation of the normals
fsign(1)=one
fsign(2)=one
fsign(6)=one
fsign(3)=-one
fsign(4)=-one
fsign(5)=-one

! compute GLL points and weights
call zwgljd(xigll,wxgll,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(etagll,wygll,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(zetagll,wzgll,ngllz,jacobi_alpha,jacobi_beta)

allocate(fgdof(ndim*maxngll2d),ftracload(ndim*maxngll2d))
allocate(dshape_quad4(2,4,maxngll2d),dshape_quad4_xy(2,4,ngllxy),              &
dshape_quad4_yz(2,4,ngllyz),dshape_quad4_zx(2,4,ngllzx))
allocate(gll_weights_xy(ngllxy),gll_points_xy(2,ngllxy),                       &
lagrange_gll_xy(ngllxy,ngllxy),dlagrange_gll_xy(2,ngllxy,ngllxy))
allocate(gll_weights_yz(ngllyz),gll_points_yz(2,ngllyz),                       &
lagrange_gll_yz(ngllyz,ngllyz),dlagrange_gll_yz(2,ngllyz,ngllyz))
allocate(gll_weights_zx(ngllzx),gll_points_zx(2,ngllzx),                       &
lagrange_gll_zx(ngllzx,ngllzx),dlagrange_gll_zx(2,ngllzx,ngllzx))
allocate(gll_weights(maxngll2d),lagrange_gll(maxngll2d,maxngll2d),             &
dlagrange_gll(2,maxngll2d,maxngll2d))
call dshape_function_quad4(4,ngllx,nglly,xigll,etagll,dshape_quad4_xy)
call dshape_function_quad4(4,nglly,ngllz,etagll,zetagll,dshape_quad4_yz)
call dshape_function_quad4(4,ngllz,ngllx,zetagll,xigll,dshape_quad4_zx)

call gll_quadrature2d(2,ngllx,nglly,ngllxy,gll_points_xy,gll_weights_xy,       &
lagrange_gll_xy,dlagrange_gll_xy)
call gll_quadrature2d(2,nglly,ngllz,ngllyz,gll_points_yz,gll_weights_yz,       &
lagrange_gll_yz,dlagrange_gll_yz)
call gll_quadrature2d(2,ngllz,ngllx,ngllzx,gll_points_zx,gll_weights_zx,       &
lagrange_gll_zx,dlagrange_gll_zx)

!read(11,*)ntrac
trac_stat=.true. ! necessary for empty trfile
count_trac=0
traction: do ! i_trac=1,ntrac
  read(11,*,iostat=ios)tractype
  if(ios/=0)exit traction
  count_trac=count_trac+1
  trac_stat=.false.

  if(tractype==0)then ! point loading
    read(11,*)q ! vector
    read(11,*)nface ! number of points
    do i_face=1,nface

      read(11,*)ielmt,inode
      ngdof=gdof(:,g_num(gnod(inode),ielmt))
      load(ngdof)=load(ngdof)+q
    enddo
    trac_stat=.true.
  elseif(tractype==1)then ! uniform loading
    read(11,*)q ! vector
    read(11,*)nface
    do i_face=1,nface
      read(11,*)ielmt,iface
      if(iface==1 .or. iface==3)then
        nfgll=ngllzx
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_zx
        gll_weights(1:nfgll)=gll_weights_zx
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
      elseif(iface==2 .or. iface==4)then
        nfgll=ngllyz
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_yz
        gll_weights(1:nfgll)=gll_weights_yz
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
      elseif(iface==5 .or. iface==6)then
        nfgll=ngllxy
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
        gll_weights(1:nfgll)=gll_weights_xy
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
      else
        write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        exit traction
      endif
      nfdof=nfgll*ndim

      num=g_num(:,ielmt)
      coord=g_coord(:,num(face(iface)%gnod))
      fgdof(1:nfdof)=reshape(gdof(:,g_num(face(iface)%nod,ielmt)),(/nfdof/))

      ftracload=zero
      ! compute numerical integration
      do i_gll=1,nfgll
        ! compute d(area)
        dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
        dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))
        ! Normal
        face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
        face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
        face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

        detjac=sqrt(dot_product(face_normal,face_normal))
        face_normal=fsign(iface)*face_normal/detjac
        !stop

        ! TODO:for constant q this can be computed only once!!
        ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
        q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(1) !only in X direction
        ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
        q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(2) !only in Y direction
        ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
        q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(3) !only in Z direction
      enddo
      load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
    enddo
    trac_stat=.true.
  elseif(tractype==2)then ! linearly distributed loading
    read(11,*)iaxis,x1,x2,q1,q2
    ! q1 and q2 are vectors, x1 and x2 can be any coordinates
    dq_dx=(q2-q1)/(x2-x1)
    read(11,*)nface
    do i_face=1,nface
      read(11,*)ielmt,iface
      if(iface==1 .or. iface==3)then
        nfgll=ngllzx
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_zx
        gll_weights(1:nfgll)=gll_weights_zx
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
      elseif(iface==2 .or. iface==4)then
        nfgll=ngllyz
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_yz
        gll_weights(1:nfgll)=gll_weights_yz
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
      elseif(iface==5 .or. iface==6)then
        nfgll=ngllzx
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
        gll_weights(1:nfgll)=gll_weights_xy
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
      else
        write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        exit traction
      endif
      nfdof=nfgll*ndim

      num=g_num(:,ielmt)
      coord=g_coord(:,num(face(iface)%gnod))
      fgdof(1:nfdof)=reshape(gdof(:,g_num(face(iface)%nod,ielmt)),(/nfdof/))
      ftracload=zero
      ! compute numerical integration
      do i_gll=1,nfgll
        x=g_coord(iaxis,num(face(iface)%nod(i_gll)))
        q=q1+dq_dx*(x-x1) ! vector of nodal values

        ! compute two vectors dx_dxi and dx_deta
        dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
        dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

        ! Normal = (dx_sxi x dx_deta)
        face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
        face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
        face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

        detjac=sqrt(dot_product(face_normal,face_normal))
        face_normal=fsign(iface)*face_normal/detjac

        ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
        q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(1) !only in X direction
        ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
        q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(2) !only in Y direction
        ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
        q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(3) !only in Z direction
     enddo
     load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
   enddo
   trac_stat=.true.
  else
    write(errtag,'(a)')'ERROR: traction type ',tractype,' not supported!'
    exit traction
  endif

enddo traction
close(11)

deallocate(dshape_quad4,dshape_quad4_xy,dshape_quad4_yz,dshape_quad4_zx)
deallocate(gll_weights_xy,gll_points_xy,lagrange_gll_xy,dlagrange_gll_xy)
deallocate(gll_weights_yz,gll_points_yz,lagrange_gll_yz,dlagrange_gll_yz)
deallocate(gll_weights_zx,gll_points_zx,lagrange_gll_zx,dlagrange_gll_zx)
deallocate(gll_weights,lagrange_gll,dlagrange_gll)
do i=1,6
  deallocate(face(i)%nod)
enddo
if(.not.trac_stat)then
  write(errtag,'(a)')'ERROR: all tractions cannot be read!'
  return
endif

errcode=0

return

end subroutine apply_traction
!===============================================================================
