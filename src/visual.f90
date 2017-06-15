! this module contains visualization routines
! REVISION:
!  HNG, Mar 11,2011; HNG, Apr 09,2010
! TODO:
!  define the array in main program such that entire array can be written
!  without transpose in this routine (done!)
module visual
contains
! this subroutine writes an ensight geofile only upto coordinates and returns
! the file unit to the calling program so that the calling program can writes
! the remaining part (connectivity) of the geo file and close it.
subroutine write_ensight_geocoord(out_fname,destag,npart,nnode,coord,funit)
character(len=250),intent(in) :: out_fname
character(len=80),intent(in) :: destag
integer,intent(in) :: npart,nnode
real,dimension(3,nnode),intent(in) :: coord !3D
integer,intent(out) :: funit

character(len=80) :: buffer ! this must be 80 characters long
integer :: i,ios

funit=11
open(unit=funit,file=trim(out_fname),access='stream',form='unformatted',       &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(funit)buffer
buffer='Created by write_ensight Routine'
write(funit)buffer
buffer='semfem3d'
write(funit)buffer
buffer='node id off'
write(funit)buffer
buffer='element id off'
write(funit)buffer
!call write_string(extent_str//char(0),fd)
!do j=1,3
!  do i=1,2
!    call write_float(real(extent(i,j)),fd)
!  enddo
!enddo
buffer='part'
write(funit)buffer
write(funit)npart
buffer=destag !'unstructured meshes'
write(funit)buffer
buffer='coordinates'
write(funit)buffer
write(funit)nnode
do i=1,3
  write(11)coord(i,:)
enddo
return
end subroutine write_ensight_geocoord
!============================================

! this subroutine writes an ensight geo file that consists of the mesh
! information
subroutine write_ensight_geo(out_fname,etype,destag,npart,nelmt,nnode,coord,   &
connect)
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
character(len=80),intent(in) :: destag
integer,intent(in) :: npart,nelmt,nnode
real,dimension(3,nnode),intent(in) :: coord !3D
integer,dimension(8,nelmt),intent(in) :: connect !hexahedral elements

character(len=80) :: buffer ! this must be 80 characters long
integer :: i,ios
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(11)buffer
buffer='Created by write_ensight Routine'
write(11)buffer
buffer='semfem3d'
write(11)buffer
buffer='node id off'
write(11)buffer
buffer='element id off'
write(11)buffer
!call write_string(extent_str//char(0),fd)
!do j=1,3
!  do i=1,2
!    call write_float(real(extent(i,j)),fd)
!  enddo
!enddo
buffer='part'
write(11)buffer
write(11)npart
buffer=destag !'unstructured meshes'
write(11)buffer
buffer='coordinates'
write(11)buffer
write(11)nnode
do i=1,3
  write(11)coord(i,:)
enddo
! writes element information
buffer=etype
write(11)buffer
write(11)nelmt

! do not substract 1 for ensight file
write(11)connect
close(11)
return
end subroutine write_ensight_geo
!============================================

! this subroutines writes ensight gold per-node variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_pernode(out_fname,destag,npart,ncomp,n,var)
character(len=250),intent(in) :: out_fname
character(len=80),intent(in) :: destag
integer,intent(in) :: npart,ncomp,n
real,dimension(ncomp,n),intent(in) :: var

character(len=80) :: buffer ! this must be 80 characters long
integer :: i,ios

if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(/,a)')'ERROR: invalid ncomp for ensight per_node variable!'
  print*,ncomp
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=destag
write(11)buffer
buffer='part'
write(11)buffer
write(11)npart
buffer='coordinates'
write(11)buffer
do i=1,ncomp
  write(11)var(i,:)
enddo

close(11)
return
end subroutine write_ensight_pernode
!============================================

end module visual
