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
integer :: ngnod ! number of geometrical nodes. usually, for FEM ngnod=nenod
integer :: nnode,nelmt
integer,allocatable :: mat_id(:)
integer,allocatable :: g_num(:,:)
integer,allocatable :: gdof(:,:)
! acceleration due to gravity
real(kind=kreal),parameter :: mtokm=0.001_kreal,agrav=9.81_kreal
real(kind=kreal),allocatable :: g_coord(:,:)
integer :: nmat
real(kind=kreal),allocatable :: gam(:),rho(:),ym(:),nu(:),coh(:),phi(:),psi(:)
logical,allocatable :: water(:)

logical :: allelastic,iseqload,iswater,istraction,phinu
! pseudostatic coefficients for earthquake loading eqkh=ah/g, eqkv=av/g
real(kind=kreal) :: eqkx,eqky,eqkz
! where ah and av are horizontal and vertical pseduostatic accelerations

integer,parameter :: nst=6 ! number of unique stress components
! sx,sy,sz,tauxy,tauyz,tauzx

character(len=250) :: file_head,inp_path,out_path,part_path
! displacement BC, ghost, traction, and water surface files
character(len=250) :: uxfile,uyfile,uzfile,gfile,trfile,wsfile
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
