! driver routine for partitioning the mesh
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
! TODO:
!   - check why segmentation fault occurs for g95
!   - remove writing unwanted files
!   - can be optimised further the detect_ghost routine
program partmesh
!use set_precision
use string_library
use global
use serial_library
use partmesh_scotch!,only: npart,out_phead,read_mesh_files,check_valence,      &
!scotch_partitioning,write_mesh_databases
use input
implicit none
integer :: istat
character(len=256) :: inp_fname
character(len=80) :: path,ext
real :: cpu_tstart,cpu_tend

logical :: ismpi
character(len=250) :: errtag ! error message
integer :: errcode

myrank=0; nproc=1;
errtag=""; errcode=-1

ismpi=.false.
!----input and initialisation----
if (command_argument_count() <= 0) then
  write(*,*)'ERROR: no input file!'
  stop
endif

! get input file name
call get_command_argument(1, inp_fname)

call parse_file(inp_fname,path,file_head,ext)

! open file to read
open(unit=11,file=trim(inp_fname),status='old', action='read',iostat=istat)
if (istat /= 0)then
  write(*,*)'ERROR: input file "'//trim(inp_fname)//'" cannot be opened!'
  stop
endif
call read_input(ismpi,inp_fname,errcode,errtag,.true.)

npart=nproc
!bc_stat=-1
!traction_stat=-1
!mesh_stat=-1
!material_stat=-1
!
out_path=trim(part_path)
out_phead=trim(out_path)//trim(file_head)//'_'
!
!istraction=.false.
!write(*,'(a)',advance='no')'reading main input file...'
!do
!  read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
!  if (ios/=0)exit
!  ! check for blank and comment line
!  if (isblank(line) .or. iscomment(line,'#'))cycle
!
!  ! look for line continuation
!  tag=trim(line)
!  call last_char(line,tmp_char,ind)
!  if (tmp_char=='&')then
!    slen=len(line)
!    tag=trim(line(1:ind-1))
!    read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
!    tag=trim(tag)//trim(line)
!  endif
!  call first_token(tag,token)
!
!  ! read pre information
!  if (trim(token)=='preinfo:')then
!    preinfo_stat=-1;
!    call split_string(tag,',',args,narg)
!    npart=get_integer('nproc',args,narg);
!    if(npart<=1)then
!      write(*,*)'ERROR: number of processors should be greater than 1!'
!      stop
!    endif
!    call seek_string('inp_path',strval,args,narg)
!    if (.not. isblank(strval))inp_path=trim(strval)
!    slen=len_trim(inp_path)
!    if(inp_path(slen:slen)/='/')inp_path=trim(inp_path)//'/'
!    call seek_string('part_path',strval,args,narg)
!    if (.not. isblank(strval))out_path=trim(strval)
!    slen=len_trim(out_path)
!    if(out_path(slen:slen)/='/')out_path=trim(out_path)//'/'
!    out_phead=trim(out_path)//trim(file_head)//'_'
!
!    preinfo_stat=0
!    cycle
!  endif
!  ! read mesh information
!  if (trim(token)=='mesh:')then
!    call split_string(tag,',',args,narg)
!    xfile=get_string('xfile',args,narg)
!    yfile=get_string('yfile',args,narg)
!    zfile=get_string('zfile',args,narg)
!    confile=get_string('confile',args,narg)
!    idfile=get_string('idfile',args,narg)
!
!    mesh_stat=0
!    cycle
!  endif
!
!  ! read bc information
!  if (trim(token)=='bc:')then
!    call split_string(tag,',',args,narg)
!    uxfile=get_string('uxfile',args,narg)
!    uyfile=get_string('uyfile',args,narg)
!    uzfile=get_string('uzfile',args,narg)
!
!    bc_stat=0
!    cycle
!  endif
!
!  ! read traction information
!  if (trim(token)=='traction:')then
!    call split_string(tag,',',args,narg)
!    trfile=get_string('trfile',args,narg)
!    traction_stat=0
!    istraction=.true.
!    !print*,trfile
!    cycle
!  endif
!
!  ! read material list
!  if (trim(token)=='material:')then
!    material_stat=-1
!    call split_string(tag,',',args,narg)
!    matfile=get_string('matfile',args,narg)
!
!    material_stat=0
!    cycle
!  endif
!enddo ! do
!close(11)
!
!! check for material list
!!if (mat_count/=nmatblk)then
!!        write(*,'(/,a)')'ERROR: number of materials doesn''t match with total number!'
!!  stop
!!endif
!
!! check input status
!if (preinfo_stat /= 0)then
!  write(*,'(a)')'ERROR: cannot read pre information! make sure the line with "preinfo:" token is correct.'
!  stop
!endif
!
!! check input status
!if (mesh_stat /= 0)then
!  write(*,'(a)')'ERROR: cannot read mesh information! make sure the line with "mesh:" token is correct.'
!  stop
!endif
!
!! check output status
!if (bc_stat /= 0)then
!  write(*,'(a)')'ERROR: cannot read BC information! make sure the line with "bc:" token is correct.'
!  stop
!endif
!
!! check material status
!if (material_stat /= 0)then
!  write(*,'(a)')'ERROR: cannot read material information! make sure the line with "material:" token is correct.'
!  stop
!endif
!
!write(*,*)'complete!'
write(*,*)'partition directory: ',trim(out_path)

! starting timer
call cpu_time(cpu_tstart)

! reads in (CUBIT) mesh files: mesh_file,nodes_coord_file, ...
call read_mesh_files()

! checks valence of nodes
  call check_valence()

! partitions mesh
  call scotch_partitioning()

! writes out database files
  call write_mesh_databases()

  write(*,*)'mesh partitioning finished successfully!'

! elapsed time
  call cpu_time(cpu_tend)
  write(*,*)'total elapsed time: ',cpu_tend-cpu_tstart,'s'
  write(*,*)'-----------------------------------------'

end program partmesh
!=======================================================

! parse file name and stop path, file head, and extension
subroutine parse_file(fname,path,head,ext)
character(len=*),intent(in) :: fname
character(len=*),intent(out) :: path,head,ext
integer :: i,ipath,iext,slen

slen=len(fname)

! set default output
path=""
head=""
ext=""

if(len_trim(fname)==0)stop

! find dot position
iext=slen+1
do i=slen,1,-1
  if(fname(i:i)==".")then
    iext=i
    exit
  endif
enddo

! find slash position
ipath=0
do i=iext,1,-1
  if(fname(i:i)=="/" .or. fname(i:i)=="\")then
    ipath=i
    exit
  endif
enddo

! determine path, head, and extension
head=fname(ipath+1:iext-1)
path=fname(1:ipath)
ext=fname(iext+1:slen)
stop
end subroutine parse_file
!=======================================================

