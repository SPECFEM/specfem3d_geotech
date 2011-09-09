! MPI math library
module math_library_mpi
use math_constants
use set_precision_mpi

private :: iminscal,fminscal
private :: iminvec,fminvec
private :: isumscal,fsumscal
private :: imaxscal,fmaxscal
private :: imaxvec,fmaxvec

! global sum of a scalar in all processors
interface sumscal
  module procedure isumscal
  module procedure fsumscal
end interface

! global maximum of a scalar in all processors
interface minscal
  module procedure iminscal
  module procedure fminscal
end interface

! global maximum of a scalar in all processors
interface maxscal
  module procedure imaxscal
  module procedure fmaxscal
end interface

! global maximum of a vector in all processors
interface maxvec
  module procedure imaxvec
  module procedure fmaxvec
end interface

! global minimum of a scalar in all processors
interface minvec
  module procedure iminvec
  module procedure fminvec
end interface
contains
!=======================================================
!=======================================================

function iminscal(scal) result(gmin)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmin
integer :: ierr

call MPI_ALLREDUCE(scal,gmin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

return
end function iminscal
!=======================================================

function fminscal(scal) result(gmin)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmin
integer :: ierr

call MPI_ALLREDUCE(scal,gmin,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)

return
end function fminscal
!=======================================================

function imaxscal(scal) result(gmax)
!
! this finds a global maximum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmax
integer :: ierr

call MPI_ALLREDUCE(scal,gmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

return
end function imaxscal
!=======================================================

function fmaxscal(scal) result(gmax)
!
! this finds a global maximum of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmax
integer :: ierr

call MPI_ALLREDUCE(scal,gmax,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)

return
end function fmaxscal
!=======================================================

function imaxvec(vec) result(gmax)
implicit none
integer,intent(in)::vec(:)
integer :: lmax,gmax ! local and global
integer :: ierr

lmax=maxval(vec)

call MPI_ALLREDUCE(lmax,gmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

return
end function imaxvec
!=======================================================

function fmaxvec(vec) result(gmax)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmax,gmax ! local and global
integer :: ierr

lmax=maxval(vec)

call MPI_ALLREDUCE(lmax,gmax,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)

return
end function fmaxvec
!=======================================================

function iminvec(vec) result(gmin)
implicit none
integer,intent(in)::vec(:)
integer :: lmin,gmin ! local and global
integer :: ierr

lmin=minval(vec)

call MPI_ALLREDUCE(lmin,gmin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

return
end function iminvec
!=======================================================

function fminvec(vec) result(gmin)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmin,gmin ! local and global
integer :: ierr

lmin=minval(vec)

call MPI_ALLREDUCE(lmin,gmin,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)

return
end function fminvec
!=======================================================

function isumscal(scal) result(gsum)
!
! this finds a global summation of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gsum
integer :: ierr

call MPI_ALLREDUCE(scal,gsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

return
end function isumscal
!=======================================================

function fsumscal(scal) result(gsum)
!
! this finds a global summation of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gsum
integer :: ierr

call MPI_ALLREDUCE(scal,gsum,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)

return
end function fsumscal
!=======================================================

function dot_product_par(vec1,vec2) result(gdot)
!
! this finds global dot product of two vectors across the processors
!
implicit none
real(kind=kreal),intent(in)::vec1(:),vec2(:)
real(kind=kreal) :: ldot,gdot
integer :: ierr

! find local dot
ldot=dot_product(vec1,vec2)
call MPI_ALLREDUCE(ldot,gdot,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)

return
end function dot_product_par
!=======================================================

end module math_library_mpi
