module integration
use set_precision
implicit none

real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4_xy,              &
dshape_quad4_yz,dshape_quad4_zx
real(kind=kreal),dimension(:),allocatable :: gll_weights_xy,                   &
gll_weights_yz,gll_weights_zx
real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll_xy,                &
lagrange_gll_yz,lagrange_gll_zx
real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll_xy,             &
dlagrange_gll_yz,dlagrange_gll_zx

contains

!-------------------------------------------------------------------------------
subroutine prepare_integration2d
use global,only:ngllx,nglly,ngllz,ngllxy,ngllyz,ngllzx
use gll_library,only:gll_quadrature2d,zwgljd
use shape_library,only:dshape_function_quad4
implicit none
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal
real(kind=kreal),dimension(ngllx) :: xigll,wxgll !double precision
real(kind=kreal),dimension(nglly) :: etagll,wygll !double precision
real(kind=kreal),dimension(ngllz) :: zetagll,wzgll !double precision
real(kind=kreal),dimension(:,:),allocatable :: gll_points_xy,gll_points_yz,    &
gll_points_zx

! compute GLL points and weights
call zwgljd(xigll,wxgll,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(etagll,wygll,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(zetagll,wzgll,ngllz,jacobi_alpha,jacobi_beta)

allocate(dshape_quad4_xy(2,4,ngllxy),dshape_quad4_yz(2,4,ngllyz),              &
dshape_quad4_zx(2,4,ngllzx))
allocate(gll_weights_xy(ngllxy),gll_points_xy(2,ngllxy),                       &
lagrange_gll_xy(ngllxy,ngllxy),dlagrange_gll_xy(2,ngllxy,ngllxy))
allocate(gll_weights_yz(ngllyz),gll_points_yz(2,ngllyz),                       &
lagrange_gll_yz(ngllyz,ngllyz),dlagrange_gll_yz(2,ngllyz,ngllyz))
allocate(gll_weights_zx(ngllzx),gll_points_zx(2,ngllzx),                       &
lagrange_gll_zx(ngllzx,ngllzx),dlagrange_gll_zx(2,ngllzx,ngllzx))

call dshape_function_quad4(4,ngllx,nglly,xigll,etagll,dshape_quad4_xy)
call dshape_function_quad4(4,nglly,ngllz,etagll,zetagll,dshape_quad4_yz)
call dshape_function_quad4(4,ngllz,ngllx,zetagll,xigll,dshape_quad4_zx)

call gll_quadrature2d(2,ngllx,nglly,ngllxy,gll_points_xy,gll_weights_xy,       &
lagrange_gll_xy,dlagrange_gll_xy)
call gll_quadrature2d(2,nglly,ngllz,ngllyz,gll_points_yz,gll_weights_yz,       &
lagrange_gll_yz,dlagrange_gll_yz)
call gll_quadrature2d(2,ngllz,ngllx,ngllzx,gll_points_zx,gll_weights_zx,       &
lagrange_gll_zx,dlagrange_gll_zx)

deallocate(gll_points_xy,gll_points_yz,gll_points_zx)

end subroutine prepare_integration2d
!===============================================================================

end module integration
!===============================================================================
