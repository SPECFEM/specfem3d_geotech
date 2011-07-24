! MPI math library
module math_library_mpi
use math_constants
use set_precision_mpi

private :: iminvec_par,fminvec_par
private :: isumscal_par,fsumscal_par
private :: imaxscal_par,fmaxscal_par
private :: imaxvec_par,fmaxvec_par

! global sum of a scalar in all processors
interface sumscal_par
  module procedure isumscal_par
  module procedure fsumscal_par
end interface

! global maximum of a scalar in all processors
interface minscal_par
  module procedure iminscal_par
  module procedure fminscal_par
end interface

! global maximum of a scalar in all processors
interface maxscal_par
  module procedure imaxscal_par
  module procedure fmaxscal_par
end interface

! global maximum of a vector in all processors
interface maxvec_par
  module procedure imaxvec_par
  module procedure fmaxvec_par
end interface

! global minimum of a scalar in all processors
interface minvec_par
  module procedure iminvec_par
  module procedure fminvec_par
end interface
contains
!=======================================================
!=======================================================

function iminscal_par(scal) result(gmin)
!
! this finds a summation of a scalar across the processors
!
implicit none 
integer,intent(in)::scal 
integer :: gmin
integer :: ierr
  
call MPI_ALLREDUCE(scal,gmin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
 
return
end function iminscal_par
!=======================================================

function fminscal_par(scal) result(gmin)
!
! this finds a summation of a scalar across the processors
!
implicit none 
real(kind=kreal),intent(in)::scal 
real(kind=kreal) :: gmin
integer :: ierr
  
call MPI_ALLREDUCE(scal,gmin,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
 
return
end function fminscal_par
!=======================================================

function imaxscal_par(scal) result(gmax)
!
! this finds a summation of a scalar across the processors
!
implicit none 
integer,intent(in)::scal 
integer :: gmax
integer :: ierr
  
call MPI_ALLREDUCE(scal,gmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
 
return
end function imaxscal_par
!=======================================================

function fmaxscal_par(scal) result(gmax)
!
! this finds a summation of a scalar across the processors
!
implicit none 
real(kind=kreal),intent(in)::scal 
real(kind=kreal) :: gmax
integer :: ierr
  
call MPI_ALLREDUCE(scal,gmax,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
 
return
end function fmaxscal_par
!=======================================================

function imaxvec_par(vec) result(gmax)
implicit none
integer,intent(in)::vec(:)
integer :: lmax,gmax ! local and global
integer :: ierr

lmax=maxval(vec)

call MPI_ALLREDUCE(lmax,gmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

return
end function imaxvec_par
!=======================================================

function fmaxvec_par(vec) result(gmax)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmax,gmax ! local and global
integer :: ierr

lmax=maxval(vec)

call MPI_ALLREDUCE(lmax,gmax,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)

return
end function fmaxvec_par
!=======================================================

function iminvec_par(vec) result(gmin)
implicit none
integer,intent(in)::vec(:)
integer :: lmin,gmin ! local and global
integer :: ierr

lmin=minval(vec)

call MPI_ALLREDUCE(lmin,gmin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

return
end function iminvec_par
!=======================================================

function fminvec_par(vec) result(gmin)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmin,gmin ! local and global
integer :: ierr

lmin=minval(vec)

call MPI_ALLREDUCE(lmin,gmin,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)

return
end function fminvec_par
!=======================================================

function isumscal_par(scal) result(gsum)
!
! this finds a summation of a scalar across the processors
!
implicit none 
integer,intent(in)::scal 
integer :: gsum
integer :: ierr
  
call MPI_ALLREDUCE(scal,gsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
 
return
end function isumscal_par
!=======================================================

function fsumscal_par(scal) result(gsum)
!
! this finds a summation of a scalar across the processors
!
implicit none 
real(kind=kreal),intent(in)::scal 
real(kind=kreal) :: gsum
integer :: ierr
  
call MPI_ALLREDUCE(scal,gsum,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
 
return
end function fsumscal_par
!=======================================================

function dot_product_par(vec1,vec2) result(gdot)
!
! this finds dot product of two vectors across the processors
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
