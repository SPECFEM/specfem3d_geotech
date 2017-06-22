! this is a main program SPECFEM3D_GEOTECH
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
program semgeotech
! import necessary libraries
use global
use string_library, only : parse_file
use input
use mesh_spec
use element
use integration
#if (USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
use visual
! main drive routines
use slope
use excavation
implicit none
integer :: funit,i,ios,j,k
integer :: i_elmt

integer :: gnum_hex(8),node_hex8(8)

character(len=250) :: arg1,inp_fname,prog
character(len=250) :: path
character(len=20), parameter :: wild_char='********************'
character(len=20) :: ensight_etype
character(len=80) :: buffer,destag ! this must be 80 characters long
character(len=20) :: ext,format_str
character(len=250) :: case_file,geo_file,sum_file
integer :: npart,nt,tinc,tstart,twidth,ts ! ts: time set for ensight gold

real(kind=kreal) :: cpu_tstart,cpu_tend,telap,max_telap,mean_telap

integer :: tot_nelmt,max_nelmt,min_nelmt,tot_nnode,max_nnode,min_nnode

character(len=250) :: errtag ! error message
integer :: errcode
logical :: isfile,isopen ! flag to check whether the file is opened

myrank=0; nproc=1;
errtag=""; errcode=-1

call start_process()

call get_command_argument(0, prog)
if (command_argument_count() <= 0) then
  call error_stop('ERROR: no input file!')
endif

call get_command_argument(1, arg1)
if(trim(arg1)==('--help'))then
  if(myrank==0)then
    write(stdout,'(a)')'Usage:'
    write(stdout,'(a)')'For information:'
    write(stdout,'(a)')'  '//trim(prog)//' [Options]'
    write(stdout,'(a)')'  Options:'
    write(stdout,'(a)')'    --help        : Display this information.'
    write(stdout,'(a)')'    --version     : Display version information.'
    if(trim(prog).eq.'./bin/semgeotech')then
      write(stdout,'(a)')'For a serial run:'
      write(stdout,'(a)')'  '//trim(prog)//' [input_file]'
    elseif(trim(prog).eq.'./bin/psemgeotech')then
      write(stdout,'(a)')'For a parallel run:'
      write(stdout,'(a)')'  mpirun -n [NP] p'//trim(prog)//' [input_file]'
    endif
    write(stdout,'(a)')'See doc/manual_SPECFEM3D_GEOTECH.pdf for details.'
  endif
  !call sync_process
  call close_process()
elseif(trim(arg1)==('--version'))then
  if(myrank==0)then
    write(stdout,'(a)')'SPECFEM3D_GEOTECH 1.2.0'
    write(stdout,'(a)')'This is free software; see the source for copying '
    write(stdout,'(a)')'conditions.  There is NO warranty; not even for '
    write(stdout,'(a)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
  endif
  !call sync_process
  call close_process()
endif

! starting timer
call cpu_time(cpu_tstart)

! get processor tag
ptail=proc_tag()

! get input file name
call get_command_argument(1, inp_fname)

if(myrank==0)write(stdout,*)'Input file name: "',trim(inp_fname),'"'
inquire(file=trim(inp_fname),exist=isfile)
if(.not.isfile)then
 if(myrank==0)then
   write(stdout,*)'Input file: "',trim(inp_fname),'" doesn''t exist!'
 endif
 call close_process()
endif
! read input data
call read_input(inp_fname,errcode,errtag)
if(errcode/=0)call error_stop(errtag)

tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)
if(myrank==0)then
write(stdout,*)'elements => total:',tot_nelmt,' max:',&
max_nelmt,' min:',min_nelmt
write(stdout,*)'nodes    => total:',tot_nnode,' max:',&
max_nnode,' min:',min_nnode
endif

if (trim(method)/='sem')then
  write(errtag,'(a)')'ERROR: wrong input for sem3d!'
  call error_stop(errtag)
endif

call parse_file(inp_fname,path,file_head,ext)

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
  call error_stop(errtag)
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
call hex2spec(ndim,ngnode,nelmt,nnode,ngllx,nglly,ngllz,errcode,errtag)
if(errcode/=0)call error_stop(errtag)
if(myrank==0)write(stdout,*)'complete!'

tot_nelmt=sumscal(nelmt); tot_nnode=sumscal(nnode)
max_nelmt=maxscal(nelmt); max_nnode=maxscal(nnode)
min_nelmt=minscal(nelmt); min_nnode=minscal(nnode)
if(myrank==0)then
write(stdout,*)'elements => total:',tot_nelmt,' max:',max_nelmt,' min:',       &
min_nelmt
write(stdout,*)'nodes    => total:',tot_nnode,' max:',max_nnode,' min:',       &
min_nnode
endif

nenod=ngll !(ngllx*nglly*ngllz) ! number of elemental nodes (nodes per element)
! number of elemental degrees of freedom
nedof=nndof*nenod

ngllxy=ngllx*nglly
ngllyz=nglly*ngllz
ngllzx=ngllz*ngllx

maxngll2d=max(ngllxy,ngllyz,ngllzx)

! prepare
call prepare_hex(errcode,errtag)
call prepare_hexface(errcode,errtag)


case_file=trim(out_path)//trim(file_head)//trim(ptail)//'.case'
open(unit=11,file=trim(case_file),status='replace',action='write',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(case_file)//'" cannot be opened!'
  call error_stop(errtag)
endif

write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  ensight gold'

write(11,'(a)')'GEOMETRY'
if(nexcav==0)then
  write(11,'(a,a/)')'model:    ',trim(file_head)//trim(ptail)//'.geo'
else
  write(11,'(a,i10,a,a/)')'model:    ',ts,' ',trim(file_head)//'_step'//      &
  wild_char(1:twidth)//trim(ptail)//'.geo'
endif

write(11,'(a)')'VARIABLE'

if(savedata%disp)then
  write(11,'(a,i10,a,a,a,a,/)')'vector per node: ',ts,' ','displacement',' ', &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.dis'
endif
if(savedata%stress)then
  write(11,'(a,i10,a,a,a,a,/)')'tensor symm per node: ',ts,' ','stress',' ',  &
  trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.sig'
endif
if(savedata%psigma)then
  write(11,'(a,i10,a,a,a,a,/)')'vector per node: ',ts,' ','principal_stress',  &
  ' ',trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.psig'
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
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ',                     &
  'stress_concentration_factor',' ',trim(file_head)//'_step'//                 &
  wild_char(1:twidth)//trim(ptail)//'.scf'
endif
if(savedata%maxtau)then
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ','max_shear_stress',  &
  ' ',trim(file_head)//'_step'//wild_char(1:twidth)//trim(ptail)//'.mtau'
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
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))// &
',a)'

! write geo file for inital stage (original)
! open Ensight Gold geo file to store mesh data
if(nexcav==0)then
  write(geo_file,fmt=format_str)trim(out_path)//trim(file_head)//trim(ptail)// &
  '.geo'
else
  write(geo_file,fmt=format_str)trim(out_path)//trim(file_head)//'_step',0,    &
  trim(ptail)//'.geo'
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
        gnum_hex=g_num(node_hex8(map2exodus),i_elmt)
        write(funit)gnum_hex
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
  call slope3d(sum_file,format_str)
else
  ! excavation
  call excavation3d(sum_file,format_str)
endif
!-------------------------------------------------------------------------------

! clean up
call cleanup_hexface(errcode,errtag)
call cleanup_integration2d(errcode,errtag)

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
!-------------------------------------------------------------------------------

if(myrank==0)then
  write(stdout,*) ! write new line
  write(stdout,'(a)')'--------------------------------------------'
  inquire(stdout,opened=isopen)
  if(isopen)close(stdout)
endif

call sync_process
call close_process()

end program semgeotech
!===============================================================================

