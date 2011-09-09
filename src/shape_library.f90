! this module contains routines to compute shape functions and their derivatives
! for hexhedral and quadrilateral elements
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module shape_library
use set_precision
contains
! this subroutines computes the shape fucntions at gll
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
subroutine shape_function_hex8(ndim,ngnod,ngllx,nglly,ngllz,xigll,etagll,      &
zetagll,shape_hex8)
use set_precision
use math_constants
implicit none

integer,intent(in) :: ndim,ngnod,ngllx,nglly,ngllz

! gauss-lobatto-legendre points of integration
real(kind=kreal) :: xigll(ngllx)
real(kind=kreal) :: etagll(nglly)
real(kind=kreal) :: zetagll(ngllz)

! 3d shape functions
real(kind=kreal) :: shape_hex8(ngnod,ngllx,nglly,ngllz)

integer :: i,j,k,i_gnod

! location of the nodes of the 3d quadrilateral elements
!real(kind=kreal) xi,eta,gamma
real(kind=kreal) :: xip,xim,etap,etam,zetap,zetam

! for checking the 3d shape functions
real(kind=kreal) :: sum_shape

real(kind=kreal), parameter :: one_eighth = 0.125d0

! check that the parameter file is correct
if(ngnod /= 8)then
  write(*,*)'ERROR: elements must have 8 geometrical nodes!'
  stop
endif

! compute shape functions
! case of a 3d 8-node element (dhatt-touzot p. 115)
do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx

      !xi = xigll(i)
      !eta = etagll(j)
      !gamma = zetagll(k)

      xip = one + xigll(i)
      xim = one - xigll(i)

      etap = one + etagll(j)
      etam = one - etagll(j)

      zetap = one + zetagll(k)
      zetam = one - zetagll(k)

      shape_hex8(1,i,j,k) = one_eighth*xim*etam*zetam
      shape_hex8(2,i,j,k) = one_eighth*xip*etam*zetam
      shape_hex8(3,i,j,k) = one_eighth*xip*etap*zetam
      shape_hex8(4,i,j,k) = one_eighth*xim*etap*zetam
      shape_hex8(5,i,j,k) = one_eighth*xim*etam*zetap
      shape_hex8(6,i,j,k) = one_eighth*xip*etam*zetap
      shape_hex8(7,i,j,k) = one_eighth*xip*etap*zetap
      shape_hex8(8,i,j,k) = one_eighth*xim*etap*zetap
    enddo
  enddo
enddo

! check the shape functions and their derivatives

do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx

      sum_shape = zero

      do i_gnod=1,ngnod
        sum_shape = sum_shape + shape_hex8(i_gnod,i,j,k)
      enddo

      ! sum of shape functions should be one
      if(abs(sum_shape-one) >  zerotol)then
        write(*,*)'ERROR: error shape functions!'
        stop
      endif
    enddo
  enddo
enddo

end subroutine shape_function_hex8
!=======================================================

! this subroutines computes derivatives of the shape fucntions at gll
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
subroutine dshape_function_hex8(ndim,ngnod,ngllx,nglly,ngllz,xigll,etagll,     &
zetagll,dshape_hex8)
use set_precision
use math_constants
implicit none
integer,intent(in) :: ndim,ngnod,ngllx,nglly,ngllz

! gauss-lobatto-legendre points of integration
real(kind=kreal) :: xigll(ngllx)
real(kind=kreal) :: etagll(nglly)
real(kind=kreal) :: zetagll(ngllz)

! derivatives of the 3d shape functions
real(kind=kreal) :: dshape_hex8(ndim,ngnod,ngllx*nglly*ngllz)

integer :: i,j,k,i_gnod
integer :: igll,ngll

! location of the nodes of the 3d quadrilateral elements
!real(kind=kreal) :: xi,eta,gamma
real(kind=kreal) :: xip,xim,etap,etam,zetap,zetam

! for checking the 3d shape functions
real(kind=kreal) :: sum_dshapexi,sum_dshapeeta,sum_dshapezeta

real(kind=kreal), parameter :: one_eighth = 0.125_kreal

! check that the parameter file is correct
if(ngnod /= 8)then
  write(*,*)'ERROR: elements must have 8 geometrical nodes!'
  stop
endif

ngll=ngllx*nglly*ngllz

! compute the derivatives of 3d shape functions
igll=0
do k=1,ngllz
  zetap = one + zetagll(k)
  zetam = one - zetagll(k)
  do j=1,nglly
    etap = one + etagll(j)
    etam = one - etagll(j)
    do i=1,ngllx
      igll=igll+1

      !xi = xigll(i)
      !eta = etagll(j)
      !gamma = zetagll(k)

      xip = one + xigll(i)
      xim = one - xigll(i)

      dshape_hex8(1,1,igll) = - one_eighth*etam*zetam
      dshape_hex8(1,2,igll) = one_eighth*etam*zetam
      dshape_hex8(1,3,igll) = one_eighth*etap*zetam
      dshape_hex8(1,4,igll) = - one_eighth*etap*zetam
      dshape_hex8(1,5,igll) = - one_eighth*etam*zetap
      dshape_hex8(1,6,igll) = one_eighth*etam*zetap
      dshape_hex8(1,7,igll) = one_eighth*etap*zetap
      dshape_hex8(1,8,igll) = - one_eighth*etap*zetap

      dshape_hex8(2,1,igll) = - one_eighth*xim*zetam
      dshape_hex8(2,2,igll) = - one_eighth*xip*zetam
      dshape_hex8(2,3,igll) = one_eighth*xip*zetam
      dshape_hex8(2,4,igll) = one_eighth*xim*zetam
      dshape_hex8(2,5,igll) = - one_eighth*xim*zetap
      dshape_hex8(2,6,igll) = - one_eighth*xip*zetap
      dshape_hex8(2,7,igll) = one_eighth*xip*zetap
      dshape_hex8(2,8,igll) = one_eighth*xim*zetap

      dshape_hex8(3,1,igll) = - one_eighth*xim*etam
      dshape_hex8(3,2,igll) = - one_eighth*xip*etam
      dshape_hex8(3,3,igll) = - one_eighth*xip*etap
      dshape_hex8(3,4,igll) = - one_eighth*xim*etap
      dshape_hex8(3,5,igll) = one_eighth*xim*etam
      dshape_hex8(3,6,igll) = one_eighth*xip*etam
      dshape_hex8(3,7,igll) = one_eighth*xip*etap
      dshape_hex8(3,8,igll) = one_eighth*xim*etap

    enddo
  enddo
enddo

! check the shape functions and their derivatives

do i=1,ngll
      sum_dshapexi = zero
      sum_dshapeeta = zero
      sum_dshapezeta = zero

      do i_gnod=1,ngnod
        sum_dshapexi = sum_dshapexi + dshape_hex8(1,i_gnod,i)
        sum_dshapeeta = sum_dshapeeta + dshape_hex8(2,i_gnod,i)
        sum_dshapezeta = sum_dshapezeta + dshape_hex8(3,i_gnod,i)
      enddo

      ! sum of derivative of shape functions should be zero
      if(abs(sum_dshapexi) >  zerotol)then
        write(*,*)'ERROR: derivative xi shape functions!'
        stop
      endif
      if(abs(sum_dshapeeta) >  zerotol)then
        write(*,*)'ERROR: derivative eta shape functions!'
        stop
      endif
      if(abs(sum_dshapezeta) >  zerotol)then
        write(*,*)'ERROR: derivative gamma shape functions!'
        stop
      endif
enddo

end subroutine dshape_function_hex8
!=======================================================

! this subroutines computes derivatives of the shape fucntions at gll
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
subroutine dshape_function_quad4(ngnod2d,ngllx,nglly,xigll,etagll,dshape_quad4)
use set_precision
use math_constants
implicit none
integer,intent(in) :: ngnod2d,ngllx,nglly

! gauss-lobatto-legendre points of integration
real(kind=kreal) :: xigll(ngllx)
real(kind=kreal) :: etagll(nglly)

! derivatives of the 3d shape functions
real(kind=kreal) :: dshape_quad4(2,ngnod2d,ngllx*nglly)

integer :: i,j,i_gnod
integer :: igll,ngll

! location of the nodes of the 2d quadrilateral elements
!real(kind=kreal) :: xi,eta
real(kind=kreal) :: xip,xim,etap,etam

! for checking the 3d shape functions
real(kind=kreal) :: sum_dshapexi,sum_dshapeeta

real(kind=kreal), parameter :: one_fourth = 0.25_kreal

! check that the parameter file is correct
if(ngnod2d /= 4)then
  write(*,*)'ERROR: element faces must have 4 geometrical nodes!'
  stop
endif

ngll=ngllx*nglly

! compute the derivatives of 2d shape functions
igll=0
do j=1,nglly
  etap = one + etagll(j)
  etam = one - etagll(j)
  do i=1,ngllx
    igll=igll+1

    xip = one + xigll(i)
    xim = one - xigll(i)

    ! corner nodes
    !shape_quad4(1,igll) = one_fourth*xim*etam
    !shape_quad4(2,igll) = one_fourth*xip*etam
    !shape_quad4(3,igll) = one_fourth*xip*etap
    !shape_quad4(4,igll) = one_fourth*xim*etap

    dshape_quad4(1,1,igll) = -one_fourth*etam
    dshape_quad4(1,2,igll) = one_fourth*etam
    dshape_quad4(1,3,igll) = one_fourth*etap
    dshape_quad4(1,4,igll) = -one_fourth*etap

    dshape_quad4(2,1,igll) = -one_fourth*xim
    dshape_quad4(2,2,igll) = -one_fourth*xip
    dshape_quad4(2,3,igll) = one_fourth*xip
    dshape_quad4(2,4,igll) = one_fourth*xim

    enddo
  enddo

  ! check the shape functions and their derivatives
  do i=1,ngll
    sum_dshapexi = zero
    sum_dshapeeta = zero

    do i_gnod=1,ngnod2d
      sum_dshapexi = sum_dshapexi + dshape_quad4(1,i_gnod,i)
      sum_dshapeeta = sum_dshapeeta + dshape_quad4(2,i_gnod,i)
    enddo

    ! sum of derivative of shape functions should be zero
    if(abs(sum_dshapexi) >  zerotol)then
      write(*,*)'ERROR: derivative xi shape functions!'
      stop
    endif
    if(abs(sum_dshapeeta) >  zerotol)then
      write(*,*)'ERROR: derivative eta shape functions!'
      stop
    endif
  enddo

end subroutine dshape_function_quad4
!=======================================================
end module shape_library
