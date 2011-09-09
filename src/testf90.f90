! this program is intended to test the Fortran compiler
program testf90
implicit none
character(len=250) :: prog
integer :: i,ios

! ghost partitions
type derived_type
  integer,dimension(:),allocatable :: acomp
end type derived_type
type(derived_type),dimension(:),allocatable :: derived

call get_command_argument(0, prog)
!----input and initialisation----
if (command_argument_count() <= 0) then
  write(*,*)'ERROR: no input file!'
endif

allocate(derived(10))
do i=1,10
  allocate(derived(i)%acomp(10))
enddo

open(unit=10,file='testf90.txt',access='stream',form='unformatted',       &
status='replace',action='write',iostat=ios)
close(10)
do i=1,10
  deallocate(derived(i)%acomp)
enddo
deallocate(derived)

end program testf90
