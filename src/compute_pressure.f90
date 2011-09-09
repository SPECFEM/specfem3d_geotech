! this subroutine computes hydrostatic pressure at staurated nodes given the freesurface profile
! z always + up
! REVISION:
!   HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO:
!   partially saturated
subroutine compute_pressure(ismpi,myid,nproc,wpressure,submerged_node,errcode,errtag)
use global
use math_library,only:determinant
implicit none
logical,intent(in) :: ismpi
integer,intent(in) :: myid,nproc
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
real(kind=kreal),dimension(nnode) :: wpressure ! water pressure
logical,dimension(nnode),intent(out) :: submerged_node
real(kind=kreal) :: A,B,C,D ! coefficients of a plane
real(kind=kreal) :: v1(3),v2(3),xf(3,4),xp(3),xmin,xmax,ymin,ymax,zmin,zmax
real(kind=kreal) :: rx1,rx2,z1,z2,z
integer :: rdir !(1: X-axis, 2: Y-axis, 3: Z-axis)
integer :: i,j,k,i1,i2,i3,i4,i5,i6,inod,i_elmt,i_face,i_node,ios
integer :: ngllxy,ngllyz,ngllzx
! water surface segments
integer :: i_wsurf,nwsurf ! number of water table surfaces
! unit weight of water
real(kind=kreal),parameter :: zero=0.0_kreal,gamw=9.81_kreal !KN/m3

character(len=20) :: format_str,ptail
character(len=250) :: fname
character(len=150) :: data_path
integer :: ipart ! partition ID

logical :: wsurf_mesh

type water_surface
  integer :: rdir,stype
  real(kind=kreal) :: rx1,rx2,z1,z2 ! z-coordinates of water surface
  integer :: nface
  integer,allocatable :: ielmt(:),iface(:)
end type water_surface
type(water_surface),allocatable :: wsurf(:)
! faces of a hexahedron
type faces
  integer,allocatable :: nod(:)
  integer :: gnod(4) ! 4 corner nodes only
end type faces
type (faces) :: face(6) ! hexahedral element has 6 faces

errtag="ERROR: unknown!"
errcode=-1
! set data path wsfile always stored in inp_path because no need to partition
!if(ismpi)then
!  data_path=trim(part_path)
!else
  data_path=trim(inp_path)
!endif

wsurf_mesh=.false.

ipart=myid-1 ! partition ID starts from 0
! all processors read only one water file from the input folder
!if(ismpi)then
!  write(format_str,*)ceiling(log10(real(nproc)+1))
!  format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
!  write(ptail, fmt=format_str)'_proc',ipart
!else
  ptail=""
!endif

ngllxy=ngllx*nglly
ngllyz=nglly*ngllz
ngllzx=ngllz*ngllx

! segment just below can be computed only once!! it is computed in apply_bc,compute_pressure, and apply_traction!!!
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

! find submerged nodes
submerged_node=.false.
do i_elmt=1,nelmt
  if(water(mat_id(i_elmt)))then
    submerged_node(g_num(:,i_elmt))=.true.
  endif
enddo

fname=trim(data_path)//trim(wsfile)
open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "',trim(fname),'" cannot be opened!'
  return
endif
read(11,*)nwsurf
allocate(wsurf(nwsurf))
do i_wsurf=1,nwsurf
  read(11,*)wsurf(i_wsurf)%stype
  if(wsurf(i_wsurf)%stype==0)then ! sweep horizontal line all across (z=constant)
    read(11,*)wsurf(i_wsurf)%rdir,wsurf(i_wsurf)%rx1,wsurf(i_wsurf)%rx2,wsurf(i_wsurf)%z1
  elseif(wsurf(i_wsurf)%stype==1)then ! sweep inclined line all across (z=linear)
    read(11,*)wsurf(i_wsurf)%rdir,wsurf(i_wsurf)%rx1,wsurf(i_wsurf)%rx2,wsurf(i_wsurf)%z1,wsurf(i_wsurf)%z2
  elseif(wsurf(i_wsurf)%stype==2)then ! meshed surface attached with the model
    read(11,*)wsurf(i_wsurf)%nface
    allocate(wsurf(i_wsurf)%ielmt(wsurf(i_wsurf)%nface))
    allocate(wsurf(i_wsurf)%iface(wsurf(i_wsurf)%nface))
    wsurf_mesh=.true.
    do i_face=1,wsurf(i_wsurf)%nface
      read(11,*)wsurf(i_wsurf)%ielmt(i_face),wsurf(i_wsurf)%iface(i_face)
    enddo
  endif
enddo
close(11)

! compute distance to the free surface
wpressure=zero
nodal: do i_node=1,nnode
  if(submerged_node(i_node))then
    xp=g_coord(:,i_node) ! coordinates of the node

    do i_wsurf=1,nwsurf
      if(wsurf(i_wsurf)%stype==0)then ! sweep horizontal line all across (z=constant)
        rdir=wsurf(i_wsurf)%rdir
        rx1=wsurf(i_wsurf)%rx1
        rx2=wsurf(i_wsurf)%rx2
        z=wsurf(i_wsurf)%z1
        !print*,rdir,rx1,rx2,z,xp(3)
        !stop
        if(xp(rdir)>=rx1 .and. xp(rdir)<=rx2 .and. xp(3)<=z)then
          ! point lies below this water surface
          ! compute pressure
          if(z>xp(3))then
            wpressure(i_node)=gamw*(z-xp(3))
          endif
          cycle nodal
        endif
      elseif(wsurf(i_wsurf)%stype==1)then ! sweep inclined line all across (z=linear)
        rdir=wsurf(i_wsurf)%rdir
        rx1=wsurf(i_wsurf)%rx1
        rx2=wsurf(i_wsurf)%rx2
        z1=wsurf(i_wsurf)%z1
        z2=wsurf(i_wsurf)%z2

        if(xp(rdir)>=min(rx1,rx2) .and. xp(rdir)<=max(rx1,rx2) .and. xp(3)<=max(z1,z2))then
          ! point lies below this water surface
          !compute z
          z=z1+(z2-z1)*(xp(rdir)-rx1)/(rx2-rx1)
          !print*,rdir,rx1,rx2,z1,z2,z,xp(rdir),xp(3)
          !stop
          ! compute pressure
          if(z>xp(3))then
            wpressure(i_node)=gamw*(z-xp(3))
          endif
          cycle nodal
        endif
      elseif(wsurf(i_wsurf)%stype==2)then ! meshed surface in the model

        do i_face=1,wsurf(i_wsurf)%nface
          xf=g_coord(:,g_num(face(wsurf(i_wsurf)%iface(i_face))%gnod,wsurf(i_wsurf)%ielmt(i_face))) ! coordinates vector of the face
          xmin=minval(xf(1,:)); xmax=maxval(xf(1,:))
          ymin=minval(xf(2,:)); ymax=maxval(xf(2,:))
          zmin=minval(xf(3,:)); zmax=maxval(xf(3,:))

          if(xp(1)>=xmin .and. xp(1)<=xmax .and. xp(2)>=ymin .and. xp(2)<=ymax .and. xp(3)<=zmax)then
            ! find equation of the plane
            ! Ax+By+Cz+D=0
            ! two vectors
            v1=xf(:,2)-xf(:,1)
            v2=xf(:,3)-xf(:,1)
            ! cross product
            A=v1(2)*v2(3)-v2(2)*v1(3)
            B=v2(1)*v1(3)-v1(1)*v2(3)
            C=v1(1)*v2(2)-v2(1)*v1(2)
            D=determinant(xf(:,1:3))

            ! find Z-coordinates on the plane just above the point (xp)
            if (C==zero)then
              write(*,*)'WARNING: free surface face is vertical!'
              z=maxval(xf(3,:))
            else
              z=-(A*xp(1)+B*xp(2)+D)/C
            endif
            if (z>maxval(xf(3,:)) .or. z<minval(xf(3,:)))then
              write(errtag,'(a)')'ERROR: free surface cannot be determined!'
              return
            endif

            ! compute pressure
            if(z>xp(3))then
              wpressure(i_node)=gamw*(z-xp(3))
            endif
            cycle nodal
          endif
        enddo
      endif
    enddo
  endif
enddo nodal
! deallocate variables
do i=1,6
  deallocate(face(i)%nod)
enddo
if(wsurf_mesh)then
  do i=1,nwsurf
    deallocate(wsurf(i)%ielmt,wsurf(i)%iface)
  enddo
endif
deallocate(wsurf)

errcode=0
return
end subroutine compute_pressure
!=======================================================

