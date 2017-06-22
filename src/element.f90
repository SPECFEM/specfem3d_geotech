! AUTHOR
!   Hom Nath Gharti
module element
use set_precision
implicit none
! map sequential node numbering to exodus/cubit order for 8-noded hexahedra
integer,parameter :: map2exodus(8)=(/ 1,2,4,3,5,6,8,7 /)
integer :: hex8_gnode(8)
real(kind=kreal) :: hexface_sign(6)
! face sign or normal orientation (outward +, inward -)

type hex_faces
  integer,allocatable :: node(:)
  integer :: gnode(4) ! geometric (corner) nodes only
end type hex_faces
type (hex_faces) :: hexface(6)

contains

!-------------------------------------------------------------------------------
subroutine prepare_hex(errcode,errtag)
use global,only:ngllx,ngllz,ngll,ngllxy
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

errtag="ERROR: unknown!"
errcode=-1

! geometrical nodes (corner nodes) in EXODUS/CUBIT order
! bottom nodes
hex8_gnode(1)=1;
hex8_gnode(2)=ngllx
hex8_gnode(3)=ngllxy;
hex8_gnode(4)=hex8_gnode(3)-ngllx+1
! top nodes
hex8_gnode(5)=(ngllz-1)*ngllxy+1;
hex8_gnode(6)=hex8_gnode(5)+ngllx-1
hex8_gnode(7)=ngll;
hex8_gnode(8)=hex8_gnode(7)-ngllx+1

errcode=0
return

end subroutine prepare_hex
!===============================================================================

subroutine prepare_hexface(errcode,errtag)
use global,only:ngllx,nglly,ngllz,ngllxy,ngllyz,ngllzx
use math_constants,only:ONE
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i,i1,i2,i3,i4,i5,i6,j,k
integer :: i_face,inode

errtag="ERROR: unknown!"
errcode=-1

allocate(hexface(1)%node(ngllzx),hexface(3)%node(ngllzx))
allocate(hexface(2)%node(ngllyz),hexface(4)%node(ngllyz))
allocate(hexface(5)%node(ngllxy),hexface(6)%node(ngllxy))

! local node numbers for the faces (faces are numbered in exodus/CUBIT
! convention)
inode=0
i1=0; i2=0; i3=0; i4=0; i5=0; i6=0
do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx
      inode=inode+1
      if (i==1)then
        ! face 4
        i4=i4+1
        hexface(4)%node(i4)=inode
      endif
      if (i==ngllx)then
        ! face 2
        i2=i2+1
        hexface(2)%node(i2)=inode
      endif
      if (j==1)then
        ! face 1
        i1=i1+1
        hexface(1)%node(i1)=inode
      endif
      if (j==nglly)then
        ! face 3
        i3=i3+1
        hexface(3)%node(i3)=inode
      endif
      if (k==1)then
        ! face 5
        i5=i5+1
        hexface(5)%node(i5)=inode
      endif
      if (k==ngllz)then
        ! face 6
        i6=i6+1
        hexface(6)%node(i6)=inode
      endif
    enddo
  enddo
enddo

! find geometric corners nodes
do i_face=1,6 ! there are 6 faces in a hexahedron
  if(i_face==1 .or. i_face==3)then ! ZX plane
    hexface(i_face)%gnode(1)=hexface(i_face)%node(1)
    hexface(i_face)%gnode(2)=hexface(i_face)%node(ngllx)
    hexface(i_face)%gnode(3)=hexface(i_face)%node(ngllzx)
    hexface(i_face)%gnode(4)=hexface(i_face)%node(ngllzx-ngllx+1)
  elseif(i_face==2 .or. i_face==4)then ! YZ plane
    hexface(i_face)%gnode(1)=hexface(i_face)%node(1)
    hexface(i_face)%gnode(2)=hexface(i_face)%node(nglly)
    hexface(i_face)%gnode(3)=hexface(i_face)%node(ngllyz)
    hexface(i_face)%gnode(4)=hexface(i_face)%node(ngllyz-nglly+1)
  elseif(i_face==5 .or. i_face==6)then ! XY plane
    hexface(i_face)%gnode(1)=hexface(i_face)%node(1)
    hexface(i_face)%gnode(2)=hexface(i_face)%node(ngllx)
    hexface(i_face)%gnode(3)=hexface(i_face)%node(ngllxy)
    hexface(i_face)%gnode(4)=hexface(i_face)%node(ngllxy-ngllx+1)
  else
    write(errtag,'(a)')'ERROR: wrong face ID for traction!'
    return
  endif
enddo

! orientation of the normals
hexface_sign(1)=one
hexface_sign(2)=one
hexface_sign(6)=one
hexface_sign(3)=-one
hexface_sign(4)=-one
hexface_sign(5)=-one

errcode=0
return

end subroutine prepare_hexface
!===============================================================================

subroutine cleanup_hexface(errcode,errtag)
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i

errtag="ERROR: unknown!"
errcode=-1

do i=1,6
  deallocate(hexface(i)%node)
enddo

errcode=0

return

end subroutine cleanup_hexface

end module element
!===============================================================================
