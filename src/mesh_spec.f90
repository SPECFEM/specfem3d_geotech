! this module contains routines to create spectral elements from the eight-noded
! hexahedral elements
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module mesh_spec

use set_precision
private :: rank, get_global,get_global_indirect_addressing, swap_all

contains

! this subroutine convert all hexahedral meshes (8-noded) to spectral elements
! of arbitrary order defined by ngllx, nglly, and ngllz
subroutine hex2spec(ndim,ngnod,nelmt,nnode,ngllx,nglly,ngllz,errcode,errtag)
use global,only : g_coord,g_num
use shape_library,only : shape_function_hex8
use gll_library, only : zwgljd
!use math_constants

implicit none
integer,intent(in) :: ndim,ngnod,ngllx,nglly,ngllz
integer,intent(inout) :: nelmt,nnode
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: ngll
real(kind=kreal),parameter :: jacobi_alpha=0.0d0,jacobi_beta=0.0d0,zero=0.0d0 !double precision
integer :: i,i_elmt,i_gnod,j,k
real(kind=kreal),dimension(:),allocatable :: xstore,ystore,zstore
real(kind=kreal),dimension(ngnod,ngllx,nglly,ngllz) :: shape_hex8

real(kind=kreal) :: xgll,ygll,zgll
real(kind=kreal),dimension(ngllx) :: xigll,wxgll !double precision
real(kind=kreal),dimension(nglly) :: etagll,wygll !double precision
real(kind=kreal),dimension(ngllz) :: zetagll,wzgll !double precision
real(kind=kreal) :: xmin,xmax

integer :: ipoint,npoint
integer :: istat
integer :: ienode,inode

integer, dimension(:), allocatable :: iglob

errtag="ERROR: unknown!"
errcode=-1

ngll=ngllx*nglly*ngllz

xmin=minval(g_coord(1,:))
xmax=maxval(g_coord(1,:))

npoint=nelmt*(ngllx*nglly*ngllz)
allocate(xstore(npoint),ystore(npoint),zstore(npoint),stat=istat)
if(istat/=0)then
  write(errtag,'(a)')'ERROR: cannot allocate memory!'
  return
endif
!xigll=0.0d0
! get gll points and weights
call zwgljd(xigll,wxgll,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(etagll,wygll,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(zetagll,wzgll,ngllz,jacobi_alpha,jacobi_beta)

! if number of points is odd, the middle abscissa is exactly zero
if(mod(ngllx,2) /= 0) xigll((ngllx-1)/2+1) = zero
if(mod(nglly,2) /= 0) etagll((nglly-1)/2+1) = zero
if(mod(ngllz,2) /= 0) zetagll((ngllz-1)/2+1) = zero

! get shape function for 8-noded hex
call shape_function_hex8(ngnod,ngllx,nglly,ngllz,xigll,etagll,zetagll,shape_hex8)

! compute coordinates all local gll points
xstore=zero
ystore=zero
zstore=zero

ipoint=0
do i_elmt=1,nelmt
  do k=1,ngllz
    do j=1,nglly
      do i=1,ngllx
        xgll = zero
        ygll = zero
        zgll = zero

        do i_gnod=1,ngnod
          xgll = xgll + shape_hex8(i_gnod,i,j,k)*g_coord(1,g_num(i_gnod,i_elmt))
          ygll = ygll + shape_hex8(i_gnod,i,j,k)*g_coord(2,g_num(i_gnod,i_elmt))
          zgll = zgll + shape_hex8(i_gnod,i,j,k)*g_coord(3,g_num(i_gnod,i_elmt))
        enddo

        ipoint=ipoint+1

        xstore(ipoint) = xgll
        ystore(ipoint) = ygll
        zstore(ipoint) = zgll

      enddo
    enddo
  enddo
enddo
deallocate(g_coord,g_num) ! no longer need these

!open(11,file='test.vtk',status='replace',action='write')

!write(11,'(a)')'# vtk DataFile Version 2.0'
!write(11,'(a)')'test'
!write(11,'(a)')'ASCII'
!write(11,'(a)')'DATASET POLYDATA'
!write(11,'(a,i7,a)')'POINTS',npoint,' float'
!do i=1,npoint
!  write(11,*)xstore(i),ystore(i),zstore(i)
!enddo
!close(11)

allocate(iglob(npoint))

! gets ibool indexing from local (gll points) to global points
call get_global(ndim,xstore,ystore,zstore,iglob,nnode,npoint,xmin,xmax)

! now we got the new number of nodes (nnode)
!- we can create a new indirect addressing to reduce cache misses
call get_global_indirect_addressing(nnode,npoint,iglob)

allocate(g_coord(3,nnode),g_num(ngll,nelmt)) ! alocate with new number of nodes
ipoint=0
do i_elmt=1,nelmt
  ienode=0
  do k=1,ngllz
    do j=1,nglly
      do i=1,ngllx
        ienode=ienode+1
        ipoint=ipoint+1
        !if(ifseg(ipoint))then
        inode=iglob(ipoint) !ibool(i,j,k,i_elmt) ! iglob(locval(ipoint))
        g_num(ienode,i_elmt)=inode
        g_coord(1,inode)=xstore(ipoint)
        g_coord(2,inode)=ystore(ipoint)
        g_coord(3,inode)=zstore(ipoint)
        !endif
      enddo
    enddo
  enddo
enddo

!open(11,file='test1.vtk',status='replace',action='write')

!write(11,'(a)')'# vtk DataFile Version 2.0'
!write(11,'(a)')'test'
!write(11,'(a)')'ASCII'
!write(11,'(a)')'DATASET POLYDATA'
!write(11,'(a,i7,a)')'POINTS',nnode,' float'
!do i=1,nnode
!  write(11,*)g_coord(1,i),g_coord(2,i),g_coord(3,i)
!enddo
!close(11)
deallocate(iglob,xstore,ystore,zstore)

errcode=0
return

end subroutine hex2spec
!============================================

subroutine get_global(ndim,xold,yold,zold,iglob,nnode,npoint,xmin,xmax)

! this routine must be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by paul f. fischer

! leave the sorting subroutines in the same source file to allow for inlining

implicit none
integer :: ndim,npoint,nnode
integer :: iglob(npoint)
integer :: iloc(npoint) ! at first this is the location (indices) of all points, at the end
!it becomes the orderd location (indices) according to the sorted rank
logical :: ifseg(npoint) ! initailly all are false, afterward only the unique
!ponts become true
real(kind=kreal),intent(in) :: xold(npoint),yold(npoint),zold(npoint) !double precision
real(kind=kreal) :: xp(npoint),yp(npoint),zp(npoint) ! double precision. initially these are the
!original coordinates of all points, at the end these become the coordinates of
!ordered points
real(kind=kreal) :: xmin,xmax

integer :: i,j,ier
integer :: nseg,ioff,iseg,ig

integer, dimension(:), allocatable :: ind,ninseg,iwork
real(kind=kreal), dimension(:), allocatable :: work !double precision

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
real(kind=kreal) :: smalltol !double precision

xp=xold
yp=yold
zp=zold

! define geometrical tolerance based upon typical size of the model
smalltol = 1.e-10_kreal * abs(xmax - xmin)

! dynamically allocate arrays
  allocate(ind(npoint), &
          ninseg(npoint), &
          iwork(npoint), &
          work(npoint),stat=ier)
  if( ier /= 0 )then
    write(*,*)'ERROR: error allocating arrays!'
    stop
  endif

! establish initial pointers
  do i=1,npoint
    iloc(i)=i
  enddo

  ifseg=.false.

  nseg=1
  ifseg(1)=.true.
  ninseg(1)=npoint

  do j=1,ndim

! sort within each segment
    ioff=1
    do iseg=1,nseg
      if(j == 1) then
        call rank(xp(ioff),ind,ninseg(iseg))
      else if(j == 2) then
        call rank(yp(ioff),ind,ninseg(iseg))
      else
        call rank(zp(ioff),ind,ninseg(iseg))
      endif
      call swap_all(iloc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))
      ioff=ioff+ninseg(iseg)
    enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
    if(j == 1) then
      do i=2,npoint
        if(abs(xp(i)-xp(i-1)) > smalltol) ifseg(i)=.true.
      enddo
    else if(j == 2) then
      do i=2,npoint
        if(abs(yp(i)-yp(i-1)) > smalltol) ifseg(i)=.true.
      enddo
    else
      do i=2,npoint
        if(abs(zp(i)-zp(i-1)) > smalltol) ifseg(i)=.true.
      enddo
    endif

! count up number of different segments
    nseg=0
    do i=1,npoint
      if(ifseg(i)) then
        nseg=nseg+1
        ninseg(nseg)=1
      else
        ninseg(nseg)=ninseg(nseg)+1
      endif
    enddo
  enddo ! j=1,ndim

! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npoint
    if(ifseg(i)) ig=ig+1
    iglob(iloc(i))=ig
  enddo

  nnode=ig

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine get_global
!===========================================

! sorting routines put in same file to allow for inlining

  subroutine rank(a,ind,n)
!
! use heap sort (numerical recipes)
!
  implicit none

  integer :: n
  real(kind=kreal) :: a(n) !double precision
  integer :: ind(n)

  integer :: i,j,l,ir,indx
  real(kind=kreal) :: q !double precision

  do j=1,n
   ind(j)=j
  enddo

  if (n == 1) return

  l=n/2+1
  ir=n
  100 continue
   if (l>1) then
      l=l-1
      indx=ind(l)
      q=a(indx)
   else
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   endif
   i=l
   j=l+l
  200    continue
   if (j <= ir) then
      if (j<ir) then
         if ( a(ind(j))<a(ind(j+1)) ) j=j+1
      endif
      if (q<a(ind(j))) then
         ind(i)=ind(j)
         i=j
         j=j+j
      else
         j=ir+1
      endif
   goto 200
   endif
   ind(i)=indx
  goto 100

end subroutine rank
!===========================================

subroutine swap_all(ia,a,b,c,iw,w,ind,n)
!
! swap arrays ia, a, b and c according to addressing in array ind
!
  implicit none

  integer :: n

  integer :: ind(n)
  integer :: ia(n),iw(n)
  real(kind=kreal) :: a(n),b(n),c(n),w(n) !double precision

  integer :: i

  iw(:) = ia(:)
  w(:) = a(:)

  do i=1,n
    ia(i)=iw(ind(i))
    a(i)=w(ind(i))
  enddo

  w(:) = b(:)

  do i=1,n
    b(i)=w(ind(i))
  enddo

  w(:) = c(:)

  do i=1,n
    c(i)=w(ind(i))
  enddo

end subroutine swap_all
!===========================================

subroutine get_global_indirect_addressing(nnode,npoint,ibool)
!- we can create a new indirect addressing to reduce cache misses
! (put into this subroutine but compiler keeps on complaining that it can't vectorize loops...)

implicit none

integer :: nnode,npoint
integer, dimension(npoint) :: ibool

! mask to sort ibool
integer, dimension(nnode) :: mask_ibool
integer, dimension(npoint) :: copy_ibool_ori
integer :: inumber
integer:: i_point

mask_ibool = -1
copy_ibool_ori = ibool
! reduces misses
inumber = 0
do i_point=1,npoint
  if(mask_ibool(copy_ibool_ori(i_point)) == -1) then
    inumber = inumber + 1
    ibool(i_point) = inumber
    mask_ibool(copy_ibool_ori(i_point)) = inumber
  else
    ! use an existing point created previously
    ibool(i_point) = mask_ibool(copy_ibool_ori(i_point))
  endif
enddo
return
end subroutine get_global_indirect_addressing
!===========================================

end module mesh_spec
