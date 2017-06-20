! this program writes Ensight SOS file
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Aug 17,2010
! COMPILE:
!   gfortran write_sos.f90 -o write_sos
! USE:
!   ./write_sos [input file]
program write_sos
implicit none

character(len=80) :: inp_fname,stemp1,stemp2,stemp3,case_head,file_head,sos_path
character(len=80),allocatable,dimension(:) :: server_name,server_exec,data_path
character(len=20) :: format_str
integer :: i_proc,ios,nproc,proc_count,slen

! input and initialisation
if (command_argument_count() <= 0) then
  write(*,'(/,a)')'ERROR: no input file!'
  stop
endif

write(*,'(a)',advance='no')'running...'

! get input file name
call get_command_argument(1, inp_fname)

! open file to read
open(unit=11,file=trim(inp_fname),status='old', action='read',iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: input file "'//trim(inp_fname)//'" cannot be opened!'
  stop
endif

read(11,*) ! skip 1 line
read(11,*)file_head
read(11,*) ! skip 1 line
read(11,*)sos_path
slen=len_trim(sos_path)
if(sos_path(slen:slen)/='/')sos_path=trim(sos_path)//'/'
read(11,*) ! skip 1 line
read(11,*)nproc
allocate(server_name(nproc),server_exec(nproc),data_path(nproc))
read(11,*) ! skip 1 line

proc_count=0
do
  read(11,*,iostat=ios)stemp1,stemp2,stemp3
  if (ios/=0.or.proc_count==nproc)exit
  proc_count=proc_count+1
  server_name(proc_count)=stemp1
  server_exec(proc_count)=stemp2
  data_path(proc_count)=stemp3
enddo
close(11)
if (proc_count<1)then
  write(*,'(/,a)')'ERROR: no information for processor/s!'
  stop
endif

if(proc_count==1 .and. proc_count<nproc)then
  server_name(2:nproc)=server_name(1)
  server_exec(2:nproc)=server_exec(1)
  data_path(2:nproc)=data_path(1)
endif

write(format_str,*)ceiling(log10(real(nproc)+1))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

! open SOS file for original data
! open file to read
open(unit=11,file=trim(sos_path)//trim(file_head)//'_original.sos',            &
status='replace', action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: input file "'//trim(file_head)//&
  '_original.sos" cannot be opened!'
  stop
endif

! open SOS file for new data
! open file to read
open(unit=12,file=trim(sos_path)//trim(file_head)//'.sos',status='replace',    &
action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: input file "'//trim(file_head)//&
  '.sos" cannot be opened!'
  stop
endif
write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  master_server gold'

write(11,'(a)')'SERVERS'
write(11,'(a,i2,/)')'number of servers:    ',nproc

write(12,'(a)')'FORMAT'
write(12,'(a,/)')'type:  master_server gold'

write(12,'(a)')'SERVERS'
write(12,'(a,i2,/)')'number of servers:    ',nproc

! loop over output slices
do i_proc=1,nproc
  ! original data
  write(11,'(a,i2)')'#Server ',i_proc
  write(11,'(a)')'machine id: '//trim(server_name(i_proc))
  write(11,'(a)')'executable: '//trim(server_exec(i_proc))
  write(11,'(a)')'#login id: '
  write(11,'(a)')'data_path: '//trim(data_path(i_proc))

  ! new data
  write(12,'(a,i2)')'#Server ',i_proc
  write(12,'(a)')'machine id: '//trim(server_name(i_proc))
  write(12,'(a)')'executable: '//trim(server_exec(i_proc))
  write(12,'(a)')'#login id: '
  write(12,'(a)')'data_path: '//trim(data_path(i_proc))

  ! file header
  write(case_head,fmt=format_str)trim(file_head)//'_proc',i_proc

  write(11,'(a,/)')'casefile: '//trim(case_head)//'_original.case'

  write(12,'(a,/)')'casefile: '//trim(case_head)//'.case'

enddo
close(11)
close(12)
write(*,*)'complete!'

end program write_sos
!===============================================================================
