! this module contains global parameters and variables
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Jul 07,2011; HNG, Apr 09,2010
!  precision parameters
module set_precision
implicit none
integer,parameter :: kreal=8 !selected_real_kind(15)
end module set_precision
!===============================================================================

! this module contains math constants and parameters
module math_constants
use set_precision
implicit none
real(kind=kreal),parameter :: ZERO=0.0_kreal,HALF=0.5_kreal,ONE=1.0_kreal,     &
two=2.0_kreal,THREE=3.0_kreal,FOUR_THIRD=4.0_kreal/THREE
real(kind=kreal),parameter :: PI=3.141592653589793_kreal

! tolerance value for zero
real(kind=kreal),parameter :: INFTOL=1.0e32_kreal,ZEROTOL = 1.0e-12_kreal

! gravitational constant: G ( m^3 kg^{-1} s^{-2} )
! source: 2014 CODATA recommended values
! http://www.physics.nist.gov/cgi-bin/cuu/Value?bg
real(kind=kreal),parameter :: GRAVCONS=6.67408e-11_kreal
end module math_constants
!===============================================================================

! this module contains conversion parameters for some physical quantities
module conversion_constants
use set_precision
use math_constants,only:PI
implicit none
! CGI to SI units
real(kind=kreal),parameter :: M2KM=1e-3_kreal
real(kind=kreal),parameter :: CGI2SI_MOMENT=1e-7_kreal
! SI to CGI units
real(kind=kreal),parameter :: KM2M=1e+3_kreal
real(kind=kreal),parameter :: SI2CGI_MOMENT=1e+7_kreal
real(kind=kreal),parameter :: DEG2RAD=PI/180.0_kreal,RAD2DEG=180.0_kreal/PI
end module conversion_constants
!===============================================================================

! global parameters/variables
module global
use set_precision
implicit none
logical,parameter :: off=.false., on=.true.
character(len=3) :: method
integer,parameter :: ndim=3
integer,parameter :: nndof=3 ! number of nodal degree of freedoms - ux, uy, uz
! number of gauss-lobatto-legendre points
integer :: ngllx,nglly,ngllz,ngll
integer,parameter :: ng=8 ! number of gauss points for FEM

! rank or ID of this processor (starts from 0), number of processors
integer :: myrank,nproc

integer :: nenod, nedof ! number of elemental nodes, number of elemental
!degree of freedoms -> nedof=nndof*nenod
integer :: ngnode ! number of geometrical nodes. usually, for FEM ngnod=nenod
integer :: nnode,nelmt
integer,allocatable :: mat_id(:)

integer,allocatable :: g_num(:,:),gdof(:,:)
! acceleration due to gravity
real(kind=kreal),parameter :: mtokm=0.001_kreal,agrav=9.81_kreal
real(kind=kreal),allocatable :: g_coord(:,:)
integer :: nmatblk
integer,allocatable :: mat_domain(:),type_blk(:)
real(kind=kreal),allocatable :: gam_blk(:),rho_blk(:),ym_blk(:),nu_blk(:),     &
coh_blk(:),phi_blk(:),psi_blk(:)
character(len=60),allocatable :: mfile_blk(:)
real(kind=kreal),allocatable :: bulkmod_blk(:),shearmod_blk(:)
logical,allocatable :: water(:)
integer,parameter :: ELASTIC_DOMAIN=1,VISCOELASTIC_DOMAIN=11

logical :: allelastic,iseqload,iswater,istraction,phinu
! pseudostatic coefficients for earthquake loading eqkh=ah/g, eqkv=av/g
real(kind=kreal) :: eqkx,eqky,eqkz
! where ah and av are horizontal and vertical pseduostatic accelerations

integer,parameter :: nst=6 ! number of unique stress components
! sx,sy,sz,tauxy,tauyz,tauzx

character(len=250) :: file_head,inp_path,out_path,part_path
! displacement BC, ghost, traction, and water surface files
character(len=250) :: confile,idfile
character(len=250),dimension(3) :: coordfile
character(len=250) :: matfile,uxfile,uyfile,uzfile,gfile,trfile,wsfile
integer :: cg_maxiter,nl_maxiter,nexcav,ninc,nsrf,ntstep
real(kind=kreal) :: cg_tol,nl_tol
! Excavation ID (regions), nunber of excavation IDs (regions) in each stage
integer,allocatable :: excavid(:),nexcavid(:)
real(kind=kreal),allocatable :: srf(:) ! strength reduction factors

! initial stress
logical :: isstress0,usek0 ! use k0 to compute horizontal stress also for
! s0_type=0
real(kind=kreal) :: s0_type ! initial stress type
! 0: by default compute by SEM,
! 1: simple overburden pressure use s0+gamma*z
!    only for horizontal and homogeneous, and
! 2: read from file
real(kind=kreal) :: z_datum,s0_datum,epk0
! z-coordinate at free surface, stress at free surface, and
! earth-pressure coefficient

! save options
type savedata_options
  logical :: disp,stress,porep,psigma,maxtau,nsigma,scf,vmeps
end type savedata_options
type(savedata_options) :: savedata

! others
character(len=1),parameter :: CR=achar(13) ! carriage return to overwrite
!previous line

!character(len=80) :: phead
integer :: stdout=6

end module global
!===============================================================================
