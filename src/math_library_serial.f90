! Serial math library this equivalent to math_library_mpi
module math_library_serial
use math_constants

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
! this finds a summation of a scalar across the processors
!
implicit none 
integer,intent(in)::scal 
integer :: gmin
integer :: ierr
  
gmin=scal
 
return
end function iminscal
!=======================================================

function fminscal(scal) result(gmin)
!
! this finds a summation of a scalar across the processors
!
implicit none 
real(kind=kreal),intent(in)::scal 
real(kind=kreal) :: gmin
integer :: ierr
  
gmin=scal
 
return
end function fminscal
!=======================================================

function imaxscal(scal) result(gmax)
!
! this finds a summation of a scalar across the processors
!
implicit none 
integer,intent(in)::scal 
integer :: gmax
integer :: ierr
  
gmax=scal
 
return
end function imaxscal
!=======================================================

function fmaxscal(scal) result(gmax)
!
! this finds a summation of a scalar across the processors
!
implicit none 
real(kind=kreal),intent(in)::scal 
real(kind=kreal) :: gmax
integer :: ierr
  
gmax=scal
 
return
end function fmaxscal
!=======================================================

function imaxvec(vec) result(gmax)
implicit none
integer,intent(in)::vec(:)
integer :: lmax,gmax ! local and global
integer :: ierr

lmax=maxval(vec)

gmax=lmax

return
end function imaxvec
!=======================================================

function fmaxvec(vec) result(gmax)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmax,gmax ! local and global
integer :: ierr

lmax=maxval(vec)

gmax=lmax

return
end function fmaxvec
!=======================================================

function iminvec(vec) result(gmin)
implicit none
integer,intent(in)::vec(:)
integer :: lmin,gmin ! local and global
integer :: ierr

lmin=minval(vec)

gmin=lmin

return
end function iminvec
!=======================================================

function fminvec(vec) result(gmin)
implicit none
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmin,gmin ! local and global
integer :: ierr

lmin=minval(vec)

gmin=lmin

return
end function fminvec
!=======================================================

function isumscal(scal) result(gsum)
!
! this finds a summation of a scalar across the processors
!
implicit none 
integer,intent(in)::scal 
integer :: gsum
integer :: ierr
  
gsum=scal
 
return
end function isumscal
!=======================================================

function fsumscal(scal) result(gsum)
!
! this finds a summation of a scalar across the processors
!
implicit none 
real(kind=kreal),intent(in)::scal 
real(kind=kreal) :: gsum
integer :: ierr
  
gsum=scal
 
return
end function fsumscal
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

gdot=ldot 
  
return
end function dot_product_par
!=======================================================

end module math_library_serial
