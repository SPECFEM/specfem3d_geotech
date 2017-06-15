! this subroutine reads the input information from a structured ASCII text file
! AUTHOR
!   Hom Nath Gharti
! REVISION:
!   HNG, Jul 07,2011; HNG, Apr 09,2010
! TODO:
!   - prompt warning or error for unknown argument/s
subroutine read_input(ismpi,inp_fname,errcode,errtag)
use global
use math_constants,only:zero,zerotol
use string_library
implicit none

integer :: i
character(len=*),intent(in) :: inp_fname
logical,intent(in) :: ismpi
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
character(len=250) :: line
character(len=800) ::tag
character(len=250) :: strval,token
character(len=1) :: tmp_char
character(len=250),dimension(50) :: args
character(len=250) :: confile,idfile,matfile
character(len=250),dimension(3) :: coordfile
integer :: id,ind,ios,narg,slen

integer :: bc_stat,preinfo_stat,mesh_stat,material_stat,control_stat,         &
eqload_stat,stress0_stat,traction_stat,water_stat,save_stat
integer :: mat_count,nwmat
integer :: ielmt,i_node,inode,imat,mat_domain,tmp_nelmt,tmp_nnode

character(len=20) :: format_str,ptail
character(len=250) :: fname
character(len=250) :: data_path,mat_path

integer :: nproc_inp ! partition ID
integer :: iselastic,ismatpart,istat,ival,issave,nexcavid_all
integer,allocatable :: ivect(:)
real(kind=kreal) :: rval
real(kind=kreal),allocatable :: rvect(:)

errtag="ERROR: unknown!"
errcode=-1

if(ismpi)then
  write(format_str,*)ceiling(log10(real(nproc)+1))
  format_str='(a,i'//trim(adjustl(format_str))//'.'//                         &
  trim(adjustl(format_str))//')'
  write(ptail, fmt=format_str)'_proc',myrank
else
  ptail=""
endif

! reading main input information
if(myrank==0)write(*,'(a)',advance='no')'reading main input file...'
preinfo_stat=-1
bc_stat=-1
traction_stat=0
eqload_stat=0
water_stat=0
mesh_stat=-1
material_stat=-1
stress0_stat=-1
control_stat=-1
save_stat=0

iselastic=0
ismatpart=1
allelastic=.false.
istraction=.false.
iseqload=.false.
iswater=.false.
isstress0=.false.
usek0=.false.
phinu=.false.

savedata%disp=.false.
savedata%stress=.false.
savedata%porep=.false.
savedata%psigma=.false.
savedata%maxtau=.false.
savedata%nsigma=.false.
savedata%scf=.false.

s0_type=0 ! by default compute initial stress using SEM

mat_count=0

! default value
method='sem'
if(ismpi)then
  inp_path='./partition/'
  part_path='./partition/'
  mat_path='./partition/'
else
  inp_path='./input/'
  mat_path='./input/'
endif
out_path='./output/'

eqkx=0.0_kreal
eqky=0.0_kreal
eqkz=0.0_kreal

! default CG and NL parameters
cg_tol=zerotol; cg_maxiter=100
nl_tol=zerotol; nl_maxiter=100

nexcav=0
nsrf=1 ! number of strength reduction factors
ninc=1 ! number of load increments

! open file to read
open(unit=11,file=trim(inp_fname),status='old', action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "'//trim(inp_fname)//                 &
  '" cannot be opened!'
  return
endif
do
  ! This will read a line and proceed to next line
  read(11,'(a)',iostat=ios)line
  if (ios/=0)exit
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle

  ! look for line continuation
  tag=trim(line)
  call last_char(line,tmp_char,ind)
  if (tmp_char=='&')then
    slen=len(line)
    tag=trim(line(1:ind-1))
    ! This will read a line and proceed to next line
    read(11,'(a)',iostat=ios)line
    tag=trim(tag)//trim(line)
  endif
  call first_token(tag,token)

  ! read pre information
  if (trim(token)=='preinfo:')then
    if(preinfo_stat==1)then
      write(errtag,*)'ERROR: copy of line type preinfo: not permitted!'
      return
    endif
    preinfo_stat=-1;
    call split_string(tag,',',args,narg)
    if(ismpi)then
      nproc_inp=get_integer('nproc',args,narg);
      if(nproc_inp/=nproc)then
        write(errtag,*)'ERROR: number of processors and images must be equal!'
        return
      endif
    endif
    call seek_string('method',strval,args,narg)
    if (.not. isblank(strval))method=trim(strval)
    !method=get_string('method',args,narg)
    if (trim(method)=='sem')then
      ngllx=get_integer('ngllx',args,narg);
      nglly=get_integer('nglly',args,narg);
      ngllz=get_integer('ngllz',args,narg);
      ngll=ngllx*nglly*ngllz ! total GLL points
    endif
    nenod=get_integer('nenod',args,narg);
    ! number of elemental degrees of freedom
    ! nedof=nndof*nenod
    if (method=='sem')then
      ! number of geometrical nodes
      ngnod=get_integer('ngnod',args,narg);
    elseif(method=='fem')then
      ngnod=nenod ! default number of geometrical nodes
    else
      write(errtag,*)'ERROR: wrong value for method!'
      return
    endif
    call seek_string('inp_path',strval,args,narg)
    if (.not. isblank(strval))inp_path=trim(strval)
    slen=len_trim(inp_path)
    if(inp_path(slen:slen)/='/')inp_path=trim(inp_path)//'/'
    if(ismpi)then
      call seek_string('part_path',strval,args,narg)
      if (.not. isblank(strval))part_path=trim(strval)
      slen=len_trim(part_path)
      if(part_path(slen:slen)/='/')part_path=trim(part_path)//'/'
    endif
    call seek_string('out_path',strval,args,narg)
    if (.not. isblank(strval))out_path=trim(strval)
    slen=len_trim(out_path)
    if(out_path(slen:slen)/='/')out_path=trim(out_path)//'/'

    preinfo_stat=1
    cycle
  endif
  ! read mesh information
  if (trim(token)=='mesh:')then
    if(mesh_stat==1)then
      write(errtag,*)'ERROR: copy of line type mesh: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    coordfile(1)=get_string('xfile',args,narg)
    coordfile(2)=get_string('yfile',args,narg)
    coordfile(3)=get_string('zfile',args,narg)
    confile=get_string('confile',args,narg)
    idfile=get_string('idfile',args,narg)
    if(ismpi)gfile=get_string('gfile',args,narg)

    mesh_stat=1
    cycle
  endif

  ! read bc information
  if (trim(token)=='bc:')then
    if(bc_stat==1)then
      write(errtag,*)'ERROR: copy of line type bc: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    uxfile=get_string('uxfile',args,narg)
    uyfile=get_string('uyfile',args,narg)
    uzfile=get_string('uzfile',args,narg)

    bc_stat=1
    cycle
  endif

  ! read initial stress information
  if (trim(token)=='stress0:')then
    if(stress0_stat==1)then
      write(errtag,*)'ERROR: copy of line type stress0: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    isstress0=.true.
    call seek_integer('type',ival,args,narg,istat)
    if(istat==0)s0_type=ival
    if(s0_type==1)then
      z_datum=get_real('z0',args,narg)
      s0_datum=get_real('s0',args,narg)
    endif
    call seek_real('k0',rval,args,narg,istat)
    if(istat==0)epk0=rval
    call seek_integer('usek0',ival,args,narg,istat)
    if(istat==0)usek0=.true.
    stress0_stat=1
    cycle
  endif

  ! read traction information
  if (trim(token)=='traction:')then
    if(traction_stat==1)then
      write(errtag,*)'ERROR: copy of line type traction: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    trfile=get_string('trfile',args,narg)
    traction_stat=1
    istraction=.true.
    cycle
  endif

  ! read material list
  if (trim(token)=='material:')then
    if(material_stat==1)then
      write(errtag,*)'ERROR: copy of line type material: not permitted!'
      return
    endif
    material_stat=-1
    call split_string(tag,',',args,narg)
    call seek_integer('ispart',ival,args,narg,istat)
    if(istat==0)ismatpart=ival
    ! if not partitioned default location is ./input
    if(ismatpart==0)mat_path='./input/'
    call seek_string('matpath',strval,args,narg)
    if (.not. isblank(strval))mat_path=trim(strval)
    slen=len_trim(mat_path)
    if(mat_path(slen:slen)/='/')mat_path=trim(mat_path)//'/'
    matfile=get_string('matfile',args,narg)
    call seek_integer('allelastic',iselastic,args,narg,istat)
    if(istat==0 .and. iselastic==1)allelastic=.true.

    material_stat=1
    cycle
  endif

  ! read earthquake loading information
  if (trim(token)=='eqload:')then
    if(eqload_stat==1)then
      write(errtag,*)'ERROR: copy of line type eqload: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    eqkx=get_real('eqkx',args,narg)
    eqky=get_real('eqky',args,narg)
    eqkz=get_real('eqkz',args,narg)
    eqload_stat=1
    iseqload=.true.
    cycle
  endif

  ! read water surface information
  if (trim(token)=='water:')then
    if(water_stat==1)then
      write(errtag,*)'ERROR: copy of line type water: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    wsfile=get_string('wsfile',args,narg)
    water_stat=1
    iswater=.true.
    cycle
  endif

  ! read control information
  if (trim(token)=='control:')then
    if(control_stat==1)then
      write(errtag,*)'ERROR: copy of line type control: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    call seek_real('cg_tol',rval,args,narg,istat)
    if(istat==0)cg_tol=rval
    call seek_integer('cg_maxiter',ival,args,narg,istat)
    if(istat==0)cg_maxiter=ival
    call seek_real('nl_tol',rval,args,narg,istat)
    if(istat==0)nl_tol=rval
    call seek_integer('nl_maxiter',ival,args,narg,istat)
    if(istat==0)nl_maxiter=ival
    call seek_integer('nsrf',ival,args,narg,istat)
    if(istat==0)nsrf=ival
    allocate(srf(nsrf),rvect(nsrf))
    srf=1.0_kreal
    call seek_real_vect('srf',rvect,nsrf,args,narg,istat)
    if(nsrf>1 .and. istat/=0)then
      write(errtag,'(a)')'ERROR: argument "srf" not found or insufficient list&
      & for "srf"!'
      return
    endif
    if(istat==0)srf=rvect
    deallocate(rvect)
    if(minval(srf)<=zero)then
      write(errtag,'(a)')'ERROR: "srf" must be positive!'
      return
    endif
    ! excavation
    call seek_integer('nexcav',ival,args,narg,istat)
    if(istat==0)nexcav=ival
    if(nexcav>0)then
      allocate(nexcavid(nexcav))
      nexcavid=1 ! default -> 1 region in each stage
      allocate(ivect(nexcav))
      call seek_integer_vect('nexcavid',ivect,nexcav,args,narg,istat)
      if(istat==0)nexcavid=ivect
      deallocate(ivect)
      nexcavid_all=sum(nexcavid)
      allocate(excavid(nexcavid_all))
      excavid=get_integer_vect('excavid',nexcavid_all,args,narg)
    endif
    if(nsrf>1 .and. nexcav>1)then
      write(errtag,'(a)')'ERROR: cannot run slope stabiliy and excavation&
      & simultaneously!'
      return
    endif
    !---------------------------

    call seek_integer('ninc',ival,args,narg,istat)
    if(istat==0)ninc=ival
    call seek_integer('phinu',ival,args,narg,istat)
    if(istat==0.and.ival==1)phinu=.true.
    control_stat=1
    cycle
  endif

  ! read save options
  if (trim(token)=='save:')then
    if(save_stat==1)then
      write(errtag,*)'ERROR: copy of line type save: not permitted!'
      return
    endif
    save_stat=-1
    call split_string(tag,',',args,narg)
    call seek_integer('disp',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%disp=.true.
    call seek_integer('stress',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%stress=.true.
    call seek_integer('porep',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%porep=.true.
    call seek_integer('psigma',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%psigma=.true.
    call seek_integer('nsigma',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%nsigma=.true.
    call seek_integer('scf',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%scf=.true.
    call seek_integer('vmeps',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%vmeps=.true.

    save_stat=1
    cycle
  endif

  write(errtag,'(a)')'ERROR: invalid line type: "'//trim(token)//'"!'
  return


enddo ! do
if(.not.iswater)savedata%porep=.false.

! check input status
if (preinfo_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read pre information! make sure the line&
  & with "preinfo:" token is correct.'
  return
endif

! check input status
if (mesh_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read mesh information! make sure the line&
  & with "mesh:" token is correct.'
  return
endif

! check output status
if (bc_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read BC information! make sure the line&
  & with "bc:" token is correct.'
  return
endif

! check material status
if (material_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read material information! make sure the&
  & line with "material:" token is correct.'
  return
endif

! check control status
if (control_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read control information! make sure the&
  & line with "control:" token is correct.'
  return
endif
if(myrank==0)write(*,*)'complete!'
!-------------------------------------------------------------------------------

! set data path
if(ismpi)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

if(myrank==0)write(*,'(a)',advance='no')'reading mesh & material properties...'
! read coordinates information
do i=1,ndim
  fname=trim(data_path)//trim(coordfile(i))//trim(ptail)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  read(11,*)tmp_nnode
  if(i==1)then
    nnode=tmp_nnode
    allocate(g_coord(ndim,nnode))
  endif
  if(tmp_nnode/=nnode)then
    write(errtag,'(a)')'ERROR: total number of nodes mismatch!'
    return
  endif
  if(ismpi)then
    do i_node=1,nnode
      read(11,*)inode,g_coord(i,inode)
    enddo
  else
    do i_node=1,nnode
      read(11,*)g_coord(i,i_node)
    enddo
  endif
enddo
close(11)
! read connectivity
fname=trim(data_path)//trim(confile)//trim(ptail)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)nelmt
allocate(g_num(nenod,nelmt))
if(ismpi)then
  do i=1,nelmt
    read(11,*)ielmt,g_num(:,ielmt)
  enddo
else
  do i=1,nelmt
    read(11,*)g_num(:,i)
  enddo
endif
close(11)

! read material id
fname=trim(data_path)//trim(idfile)//trim(ptail)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)tmp_nelmt
if(tmp_nelmt/=nelmt)then
  write(errtag,'(a)')'ERROR: total number of elements mismatch!'
  return
endif
allocate(mat_id(nelmt))
if(ismpi)then
  do i=1,nelmt
    read(11,*)ielmt,mat_id(ielmt)
  enddo
else
  do i=1,nelmt
    read(11,*)mat_id(i)
  enddo
endif
close(11)

! read material lists
if(ismatpart==0)then ! material file not partitioned
  fname=trim(mat_path)//trim(matfile)
elseif(ismatpart==1)then ! material file partitioned
  fname=trim(mat_path)//trim(matfile)//trim(ptail)
else
  write(errtag,'(a)')'ERROR: illegal ismatpart value!'
  return
endif
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)
read(11,*)nmat
allocate(gam(nmat),rho(nmat),ym(nmat),coh(nmat),nu(nmat),phi(nmat),psi(nmat),  &
water(nmat))
do i=1,nmat
  read(11,*)imat,mat_domain,gam(i),ym(i),nu(i),phi(i),coh(i),psi(i)
enddo
if(minval(mat_id)<1 .or. maxval(mat_id)>nmat)then
  write(errtag,'(a)')'ERROR: material IDs must be consistent with the defined&
  & material regions!'
  return
endif
water=.false.
if(iswater)then
  read(11,*,iostat=ios)nwmat
  if( ios /= 0 ) then
    write(errtag,'(a)')'ERROR: water IDs cannot be read from material list!'
    return
  endif
  do i=1,nwmat
    read(11,*)id
    water(id)=.true.
  enddo
endif
rho=gam/9.81_kreal ! Kg/m3

! read bc
errcode=0
if(myrank==0)write(*,*)'complete!'

end subroutine read_input
!===============================================================================
