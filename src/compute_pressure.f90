! this subroutine computes hydrostatic pressure at staurated nodes given the
! freesurface profile
! z always + up
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO:
!   partially saturated
subroutine compute_pressure(wpressure,submerged_node,errcode,errtag)
use global
use math_library,only:determinant
use element,only:hexface
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
real(kind=kreal),dimension(nnode) :: wpressure ! water pressure
logical,dimension(nnode),intent(out) :: submerged_node
real(kind=kreal) :: A,B,C,D ! coefficients of a plane
real(kind=kreal) :: v1(3),v2(3),xf(3,4),xp(3),xmin,xmax,ymin,ymax,zmin,zmax
real(kind=kreal) :: rx1,rx2,z1,z2,z
integer :: rdir !(1: X-axis, 2: Y-axis, 3: Z-axis)
integer :: i,i_elmt,i_face,i_node,ios
! water surface segments
integer :: i_wsurf,nwsurf ! number of water table surfaces
! unit weight of water
real(kind=kreal),parameter :: zero=0.0_kreal,gamw=9.81_kreal !KN/m3

character(len=250) :: fname
character(len=150) :: data_path

logical :: wsurf_mesh

type water_surface
  integer :: rdir,stype
  real(kind=kreal) :: rx1,rx2,z1,z2 ! z-coordinates of water surface
  integer :: nface
  integer,allocatable :: ielmt(:),iface(:)
end type water_surface
type(water_surface),allocatable :: wsurf(:)

errtag="ERROR: unknown!"
errcode=-1

! set data path wsfile always stored in inp_path because no need to partition
data_path=trim(inp_path)

wsurf_mesh=.false.

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
  write(errtag,'(a)')'ERROR: input file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)nwsurf
allocate(wsurf(nwsurf))
do i_wsurf=1,nwsurf
  read(11,*)wsurf(i_wsurf)%stype
  if(wsurf(i_wsurf)%stype==0)then
    ! sweep horizontal line all across (z=constant)
    read(11,*)wsurf(i_wsurf)%rdir,wsurf(i_wsurf)%rx1,wsurf(i_wsurf)%rx2,       &
    wsurf(i_wsurf)%z1
  elseif(wsurf(i_wsurf)%stype==1)then
    ! sweep inclined line all across (z=linear)
    read(11,*)wsurf(i_wsurf)%rdir,wsurf(i_wsurf)%rx1,wsurf(i_wsurf)%rx2,       &
    wsurf(i_wsurf)%z1,wsurf(i_wsurf)%z2
  elseif(wsurf(i_wsurf)%stype==2)then
    ! meshed surface attached with the model
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
      if(wsurf(i_wsurf)%stype==0)then
      ! sweep horizontal line all across (z=constant)
        rdir=wsurf(i_wsurf)%rdir
        rx1=wsurf(i_wsurf)%rx1
        rx2=wsurf(i_wsurf)%rx2
        z=wsurf(i_wsurf)%z1
        if(xp(rdir)>=rx1 .and. xp(rdir)<=rx2 .and. xp(3)<=z)then
          ! point lies below this water surface
          ! compute pressure
          if(z>xp(3))then
            wpressure(i_node)=gamw*(z-xp(3))
          endif
          cycle nodal
        endif
      elseif(wsurf(i_wsurf)%stype==1)then
      ! sweep inclined line all across (z=linear)
        rdir=wsurf(i_wsurf)%rdir
        rx1=wsurf(i_wsurf)%rx1
        rx2=wsurf(i_wsurf)%rx2
        z1=wsurf(i_wsurf)%z1
        z2=wsurf(i_wsurf)%z2

        if(xp(rdir)>=min(rx1,rx2) .and. xp(rdir)<=max(rx1,rx2) .and.           &
        xp(3)<=max(z1,z2))then
          ! point lies below this water surface
          !compute z
          z=z1+(z2-z1)*(xp(rdir)-rx1)/(rx2-rx1)
          ! compute pressure
          if(z>xp(3))then
            wpressure(i_node)=gamw*(z-xp(3))
          endif
          cycle nodal
        endif
      elseif(wsurf(i_wsurf)%stype==2)then ! meshed surface in the model

        do i_face=1,wsurf(i_wsurf)%nface
          xf=g_coord(:,g_num(hexface(wsurf(i_wsurf)%iface(i_face))%gnode,          &
          wsurf(i_wsurf)%ielmt(i_face))) ! coordinates vector of the face
          xmin=minval(xf(1,:)); xmax=maxval(xf(1,:))
          ymin=minval(xf(2,:)); ymax=maxval(xf(2,:))
          zmin=minval(xf(3,:)); zmax=maxval(xf(3,:))

          if(xp(1)>=xmin .and. xp(1)<=xmax .and. xp(2)>=ymin .and. xp(2)<=ymax &
          .and. xp(3)<=zmax)then
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
if(wsurf_mesh)then
  do i=1,nwsurf
    deallocate(wsurf(i)%ielmt,wsurf(i)%iface)
  enddo
endif
deallocate(wsurf)

errcode=0
return
end subroutine compute_pressure
!===============================================================================

