! this is a main program SPECFEM3D_GEOTECH
! REVISION:
!   HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
program semgeotech
! import necessary libraries
use global
use string_library, only : parse_file
!use math_constants
use mesh_spec
#if (USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
use visual

implicit none
integer :: funit,i,ios,j,k
integer :: i_elmt

integer :: gnod(8),map2exodus(8),ngllxy,node_hex8(8)

character(len=250) :: arg1,inp_fname,prog
character(len=250) :: path
character(len=20), parameter :: wild_char='********************'
character(len=20) :: ensight_etype
character(len=80) :: buffer,destag ! this must be 80 characters long
character(len=20) :: ext,format_str,ptail
character(len=250) :: case_file,geo_file,sum_file
integer :: npart,nt,tinc,tstart,twidth,ts ! ts: time set for ensight gold

real(kind=kreal) :: cpu_tstart,cpu_tend,telap,max_telap,mean_telap

logical :: ismpi !.true. : MPI, .false. : serial
integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode

character(len=250) :: errtag ! error message
integer :: errcode
logical :: isopen ! flag to check whether the file is opened

myrank=0; nproc=1;
errtag=""; errcode=-1

call start_process(ismpi,stdout)

call get_command_argument(0, prog)
!----input and initialisation----
if (command_argument_count() <= 0) then
  call error_stop('ERROR: no input file!',stdout,myrank)
endif

call get_command_argument(1, arg1)
if(trim(arg1)==('--help'))then
  if(myrank==0)then
    write(stdout,'(a)')'Usage: '//trim(prog)//' [Options] [input_file]'
    write(stdout,'(a)')'Options:'
    write(stdout,'(a)')'    --help        : Display this information.'
    write(stdout,'(a)')'    --version     : Display version information.'
  endif
  !call sync_process
  call close_process()
elseif(trim(arg1)==('--version'))then
  if(myrank==0)then
    write(stdout,'(a)')'SPECFEM3D_GEOTECH 1.2 Beta'
    write(stdout,'(a)')'This is free software; see the source for copying '
    write(stdout,'(a)')'conditions.  There is NO warranty; not even for '
    write(stdout,'(a)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
  endif
  !call sync_process
  call close_process()
endif

! starting timer
call cpu_time(cpu_tstart)

! get input file name
call get_command_argument(1, inp_fname)

! read input data
call read_input(ismpi,myrank,inp_fname,errcode,errtag)
if(errcode/=0)call error_stop(errtag,stdout,myrank)
!call sync_process()

tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)
if(myrank==0)then
write(stdout,*)'elements => total:',tot_nelmt,' max:',max_nelmt,' min:',min_nelmt
write(stdout,*)'nodes    => total:',tot_nnode,' max:',max_nnode,' min:',min_nnode
endif

if (trim(method)/='sem')then
  write(errtag,'(a)')'ERROR: wrong input for sem3d!'
  call error_stop(errtag,stdout,myrank)
endif

call parse_file(inp_fname,path,file_head,ext)

! get processor tag
ptail=proc_tag()

ensight_etype='hexa8'
ts=1 ! time set
tstart=1; tinc=1
if(nexcav==0)then
  nt=nsrf
else
  nt=nexcav+1 ! include 0 excavation stage (i.e., initial)
  tstart=0
endif
twidth=ceiling(log10(real(nt)+1.))

! write original meshes
! write Ensight gold .case file
case_file=trim(out_path)//trim(file_head)//'_original'//trim(ptail)//'.case'
geo_file=trim(file_head)//'_original'//trim(ptail)//'.geo'

open(unit=11,file=trim(case_file),status='replace',action='write',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(case_file)//'" cannot be opened!'
  call error_stop(errtag,stdout,myrank)
endif

write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  ensight gold'

write(11,'(a)')'GEOMETRY'
write(11,'(a,a,/)')'model:    ',trim(geo_file)

close(11)

! write original mesh with 8-noded hexahedra
! open Ensight Gold geo file to store mesh data
geo_file=trim(out_path)//trim(geo_file)
npart=1
destag='unstructured meshes'
call write_ensight_geo(geo_file,ensight_etype,destag,npart,nelmt,nnode,        &
real(g_coord),g_num)

! create spectral elements
if(myrank==0)write(stdout,'(a)',advance='no')'creating spectral elements...'
call hex2spec(ndim,ngnod,nelmt,nnode,ngllx,nglly,ngllz,errcode,errtag)
if(errcode/=0)call error_stop(errtag,stdout,myrank)
if(myrank==0)write(stdout,*)'complete!'

tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)
if(myrank==0)then
write(stdout,*)'elements => total:',tot_nelmt,' max:',max_nelmt,' min:',min_nelmt
write(stdout,*)'nodes    => total:',tot_nnode,' max:',max_nnode,' min:',min_nnode
endif

!call sync_process

!stop
nenod=ngll !(ngllx*nglly*ngllz) ! number of elemental nodes (nodes per element)
! number of elemental degrees of freedom
nedof=nndof*nenod

ngllxy=ngllx*nglly

! geometrical nodes (corner nodes) in EXODUS/CUBIT order
! bottom nodes
gnod(1)=1;
gnod(2)=ngllx
gnod(3)=ngllxy;
gnod(4)=gnod(3)-ngllx+1
! top nodes
gnod(5)=(ngllz-1)*ngllxy+1;
gnod(6)=gnod(5)+ngllx-1
gnod(7)=ngll;
gnod(8)=gnod(7)-ngllx+1

! map sequential node numbering to exodus/cubit order for 8-noded hexahedra
map2exodus=(/ 1,2,4,3,5,6,8,7 /)

case_file=trim(out_path)//trim(file_head)//trim(ptail)//'.case'
open(unit=11,file=trim(case_file),status='replace',action='write',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(case_file)//'" cannot be opened!'
  call error_stop(errtag,stdout,myrank)
endif

write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  ensight gold'

write(11,'(a)')'GEOMETRY'
!write(11,'(a,a,/)')'model:    ',trim(geo_file)
if(nexcav==0)then
  write(11,'(a,a/)')'model:    ',trim(file_head)//trim(ptail)//'.geo'
else
  write(11,'(a,i10,a,a/)')'model:    ',ts,' ',trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.geo'
endif

write(11,'(a)')'VARIABLE'
!write(11,'(a,i10,a,a,a,a,/)')'vector per node: ',ts,' displacement ', &
!trim(file_head)//'_step'//wild_char(1:twidth)//'.dis'
!write(11,'(a,i10,a,a,a,a,/)')'vector per node: ',ts,' principal_stress ', &
!trim(file_head)//'_step'//wild_char(1:twidth)//'.sig'

if(savedata%disp)then
  write(11,'(a,i10,a,a,a,a,/)')'vector per node: ',ts,' ','displacement',' ',  &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.dis'
endif
if(savedata%stress)then
  write(11,'(a,i10,a,a,a,a,/)')'tensor symm per node: ',ts,' ','stress',' ',   &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.sig'
endif
if(savedata%psigma)then
  write(11,'(a,i10,a,a,a,a,/)')'vector per node: ',ts,' ','principal_stress',' ', &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.psig'
endif
if(savedata%porep)then
   write(11,'(a,a,a,a,a,/)')'scalar per node: ',' ','pore_pressure',' ', &
   trim(file_head)//trim(ptail)//'.por'
endif
if(savedata%vmeps)then
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','plastic_strain',' ', &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.eps'
endif
if(savedata%scf)then
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','stress_concentration_factor',' ', &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.scf'
endif
if(savedata%maxtau)then
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','max_shear_stress',' ', &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.mtau'
endif
if(savedata%nsigma)then
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','normal_stress',' ', &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.nsig'
endif
write(11,'(a)')'TIME'
write(11,'(a,i10)')'time set:',ts
write(11,'(a,i10)')'number of steps:',nt
write(11,'(a,i10)')'filename start number:',tstart
write(11,'(a,i10)')'filename increment:',tinc
write(11,'(a)',advance='no')'time values: '

if(nexcav==0)then
  do i=1,nt
    write(11,'(e12.5)',advance='yes')srf(i)
  enddo
else
  do i=0,nt-1
    write(11,'(e12.5)',advance='yes')real(i)
  enddo
endif
close(11)

! Format string
write(format_str,*)twidth
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//',a)'

! write geo file for inital stage (original)
! open Ensight Gold geo file to store mesh data
if(nexcav==0)then
  write(geo_file,fmt=format_str)trim(out_path)//trim(file_head)//trim(ptail)//'.geo'
else
  write(geo_file,fmt=format_str)trim(out_path)//trim(file_head)//'_step',0,trim(ptail)//'.geo'
endif
npart=1
destag='unstructured meshes'
call write_ensight_geocoord(geo_file,destag,npart,nnode,real(g_coord),funit)

! writes element information
buffer=ensight_etype
write(funit)buffer
write(funit)nelmt*(ngllx-1)*(nglly-1)*(ngllz-1)

! do not substract 1 for ensight file
do i_elmt=1,nelmt
  do k=1,ngllz-1
    do j=1,nglly-1
      do i=1,ngllx-1
        ! corner nodes in a sequential numbering
        node_hex8(1)=(k-1)*ngllxy+(j-1)*ngllx+i
        node_hex8(2)=node_hex8(1)+1

        node_hex8(3)=node_hex8(1)+ngllx
        node_hex8(4)=node_hex8(3)+1

        node_hex8(5)=node_hex8(1)+ngllxy
        node_hex8(6)=node_hex8(5)+1

        node_hex8(7)=node_hex8(5)+ngllx
        node_hex8(8)=node_hex8(7)+1
        ! map to exodus/cubit numbering and write
        write(funit)g_num(node_hex8(map2exodus),i_elmt)
      enddo
    enddo
  enddo
enddo
close(funit)

! open summary file
sum_file = trim(out_path)//trim(file_head)//'_summary'//trim(ptail)
open(unit=10,file=trim(sum_file),status='replace',action='write',iostat=ios)
write(10,*)'--------------------------------------------'
write(10,*)'Result summary produced by SPECFEM3D_GEOTECH'
write(10,*)'--------------------------------------------'
close(10)

if(myrank==0)write(stdout,'(a)')'--------------------------------------------'

! call main routines
if(nexcav==0)then
  ! slope stability
  call semslope3d(ismpi,gnod,sum_file,ptail,format_str)
else
  ! excavation
  call semexcav3d(ismpi,gnod,sum_file,ptail,format_str)
endif
!-----------------------------------

! compute elapsed time
call cpu_time(cpu_tend)
telap=cpu_tend-cpu_tstart
max_telap=maxscal(telap)
mean_telap=sumscal(telap)/real(nproc,kreal)

write(format_str,*)ceiling(log10(real(max_telap)+1.))+5 ! 1 . and 4 decimals
format_str='(3(f'//trim(adjustl(format_str))//'.4,1X))'
open(10,file=trim(sum_file),status='old',position='append',action='write')
write(10,*)'ELAPSED TIME, MAX ELAPSED TIME, MEAN ELAPSED TIME'
write(10,fmt=format_str)telap,max_telap,mean_telap
close(10)
!-----------------------------------

if(myrank==0)then
  write(stdout,*) ! write new line
  write(stdout,'(a)')'--------------------------------------------'
  inquire(stdout,opened=isopen)
  if(isopen)close(stdout)
endif

call sync_process
call close_process()

end program semgeotech
!===========================================

