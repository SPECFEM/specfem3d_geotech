! this subroutine applies displacement boundary conditions and determines
! global degrees of freedom
! AUTHOR
!   Hom Nath Gharti
! REVISION
!   HNG, Jul 12,2011; ; HNG, Apr 09,2010
subroutine apply_bc(neq,errcode,errtag)
use global
use math_constants,only:ZERO
use element,only:hexface
implicit none
integer,intent(out) :: neq
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: i,ios,j
integer :: bctype,ielmt,iface,nelpart,i_elpart
real(kind=kreal) :: val
character(len=250) :: fname
character(len=150) :: data_path

errtag="ERROR: unknown!"
errcode=-1
! set data path
if(ismpi)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

fname=trim(data_path)//trim(uxfile)//trim(ptail)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
bcux: do
  read(11,*,iostat=ios)bctype,val
  if(ios/=0)exit
  if(val/=zero)then
    write(errtag,*)'ERROR: nonzero displacement BC not implemented!'
    return
  endif
  if(bctype==0)then ! point
    write(errtag,*)'ERROR: nodal displacement BC not implemented!'
    return
  elseif(bctype==1)then ! edge
    write(errtag,*)'ERROR: edge displacement BC not implemented!'
    return
  elseif(bctype==2)then ! face
    read(11,*)nelpart
    do i_elpart=1,nelpart
      read(11,*)ielmt,iface ! This will read a line and proceed to next line
      gdof(1,g_num(hexface(iface)%node,ielmt))=0
    enddo
  else
    write(errtag,*)'ERROR: undefined displacement BC type ux!',bctype
    return
  endif
enddo bcux
close(11)

fname=trim(data_path)//trim(uyfile)//trim(ptail)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

bcuy: do
  read(11,*,iostat=ios)bctype,val
  if(ios/=0)exit
  if(val/=zero)then
    write(errtag,*)'ERROR: nonzero displacement BC not implemented!'
    return
  endif
  if(bctype==0)then ! point
    write(errtag,*)'ERROR: nodal displacement BC not implemented!'
    return
  elseif(bctype==1)then ! edge
    write(errtag,*)'ERROR: edge displacement BC not implemented!'
    return
  elseif(bctype==2)then ! face
    read(11,*)nelpart
    do i_elpart=1,nelpart
      read(11,*)ielmt,iface ! This will read a line and proceed to next line
      gdof(2,g_num(hexface(iface)%node,ielmt))=0
    enddo
  else
    write(errtag,*)'ERROR: undefined displacement BC type ux!',bctype
    return
  endif
enddo bcuy
close(11)

fname=trim(data_path)//trim(uzfile)//trim(ptail)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

bcuz: do
  read(11,*,iostat=ios)bctype,val
  if(ios/=0)exit
  if(val/=zero)then
    write(errtag,*)'ERROR: nonzero displacement BC not implemented!'
    return
  endif
  if(bctype==0)then ! point
    write(errtag,*)'ERROR: nodal displacement BC not implemented!'
    return
  elseif(bctype==1)then ! edge
    write(errtag,*)'ERROR: edge displacement BC not implemented!'
    return
  elseif(bctype==2)then ! face
    read(11,*)nelpart
    do i_elpart=1,nelpart
      read(11,*)ielmt,iface ! This will read a line and proceed to next line
      gdof(3,g_num(hexface(iface)%node,ielmt))=0
    enddo
  else
    write(errtag,*)'ERROR: undefined displacement BC type ux!',bctype
    return
  endif
enddo bcuz
close(11)

! compute modified gdof
neq=0
do j=1,ubound(gdof,2)
  do i=1,ubound(gdof,1)
    if(gdof(i,j)/=0)then
      neq=neq+1
      gdof(i,j)=neq
    endif
  enddo
enddo

! compute nodal to global
errcode=0
return

end subroutine apply_bc
!===============================================================================
