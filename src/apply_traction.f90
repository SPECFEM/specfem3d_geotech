! TODO: make use of the orthogonality
! this routine read and applies the traction specified in the traction file
! AUTHOR
!   Hom Nath Gharti
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
subroutine apply_traction(gnod,neq,load,errcode,errtag)
use global
use math_constants
use element,only:hexface,hexface_sign
use integration,only:dshape_quad4_xy,dshape_quad4_yz,dshape_quad4_zx,          &
                     gll_weights_xy,gll_weights_yz,gll_weights_zx,             &
                     lagrange_gll_xy,lagrange_gll_yz,lagrange_gll_zx!,         &
                     !dlagrange_gll_xy,dlagrange_gll_yz,dlagrange_gll_zx
implicit none
integer,intent(in) :: gnod(ngnode)
integer,intent(in) :: neq
real(kind=kreal),intent(inout) :: load(0:neq)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: ngdof(ndim)
integer :: i_face,i_gll,ios,iaxis
integer :: ielmt,iface,nface,inode
integer :: num(nenod)
integer :: nfdof,nfgll
integer,allocatable :: fgdof(:) ! face global degrees of freedom
integer :: tractype,count_trac
logical :: trac_stat

real(kind=kreal) :: coord(ndim,4)
real(kind=kreal) :: x,x1,x2,detjac
real(kind=kreal),dimension(ndim) :: face_normal,dx_dxi,dx_deta,dq_dx,q,q1,q2
real(kind=kreal),allocatable :: ftracload(:) ! face traction load

real(kind=kreal),allocatable :: dshape_quad4(:,:,:),gll_weights(:),            &
lagrange_gll(:,:),dlagrange_gll(:,:,:)

character(len=80) :: fname
character(len=80) :: data_path

errtag="ERROR: unknown!"
errcode=-1
! set data path
if(ismpi)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

fname=trim(data_path)//trim(trfile)//trim(ptail)
open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "',trim(fname),'" cannot be opened!'
  return
endif

!read(11,*)ntrac
trac_stat=.true. ! necessary for empty trfile
count_trac=0
traction: do ! i_trac=1,ntrac
  read(11,*,iostat=ios)tractype
  if(ios/=0)exit traction
  count_trac=count_trac+1
  trac_stat=.false.

  if(tractype==0)then ! point loading
    read(11,*)q ! vector
    read(11,*)nface ! number of points
    do i_face=1,nface

      read(11,*)ielmt,inode
      ngdof=gdof(:,g_num(gnod(inode),ielmt))
      load(ngdof)=load(ngdof)+q
    enddo
    trac_stat=.true.
  elseif(tractype==1)then ! uniform loading
    read(11,*)q ! vector
    read(11,*)nface
    do i_face=1,nface
      read(11,*)ielmt,iface
      if(iface==1 .or. iface==3)then
        nfgll=ngllzx
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_zx
        gll_weights(1:nfgll)=gll_weights_zx
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
      elseif(iface==2 .or. iface==4)then
        nfgll=ngllyz
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_yz
        gll_weights(1:nfgll)=gll_weights_yz
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
      elseif(iface==5 .or. iface==6)then
        nfgll=ngllxy
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
        gll_weights(1:nfgll)=gll_weights_xy
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
      else
        write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        exit traction
      endif
      nfdof=nfgll*ndim

      num=g_num(:,ielmt)
      coord=g_coord(:,num(hexface(iface)%gnode))
      fgdof(1:nfdof)=reshape(gdof(:,g_num(hexface(iface)%node,ielmt)),(/nfdof/))

      ftracload=zero
      ! compute numerical integration
      do i_gll=1,nfgll
        ! compute d(area)
        dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
        dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))
        ! Normal
        face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
        face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
        face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

        detjac=sqrt(dot_product(face_normal,face_normal))
        face_normal=hexface_sign(iface)*face_normal/detjac
        !stop

        ! TODO:for constant q this can be computed only once!!
        ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
        q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(1) !only in X direction
        ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
        q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(2) !only in Y direction
        ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
        q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(3) !only in Z direction
      enddo
      load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
    enddo
    trac_stat=.true.
  elseif(tractype==2)then ! linearly distributed loading
    read(11,*)iaxis,x1,x2,q1,q2
    ! q1 and q2 are vectors, x1 and x2 can be any coordinates
    dq_dx=(q2-q1)/(x2-x1)
    read(11,*)nface
    do i_face=1,nface
      read(11,*)ielmt,iface
      if(iface==1 .or. iface==3)then
        nfgll=ngllzx
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_zx
        gll_weights(1:nfgll)=gll_weights_zx
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
      elseif(iface==2 .or. iface==4)then
        nfgll=ngllyz
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_yz
        gll_weights(1:nfgll)=gll_weights_yz
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
      elseif(iface==5 .or. iface==6)then
        nfgll=ngllzx
        lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
        gll_weights(1:nfgll)=gll_weights_xy
        dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
      else
        write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        exit traction
      endif
      nfdof=nfgll*ndim

      num=g_num(:,ielmt)
      coord=g_coord(:,num(hexface(iface)%gnode))
      fgdof(1:nfdof)=reshape(gdof(:,g_num(hexface(iface)%node,ielmt)),(/nfdof/))
      ftracload=zero
      ! compute numerical integration
      do i_gll=1,nfgll
        x=g_coord(iaxis,num(hexface(iface)%node(i_gll)))
        q=q1+dq_dx*(x-x1) ! vector of nodal values

        ! compute two vectors dx_dxi and dx_deta
        dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
        dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

        ! Normal = (dx_sxi x dx_deta)
        face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
        face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
        face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

        detjac=sqrt(dot_product(face_normal,face_normal))
        face_normal=hexface_sign(iface)*face_normal/detjac

        ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
        q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(1) !only in X direction
        ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
        q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(2) !only in Y direction
        ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
        q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
        ! *face_normal(3) !only in Z direction
     enddo
     load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
   enddo
   trac_stat=.true.
  else
    write(errtag,'(a)')'ERROR: traction type ',tractype,' not supported!'
    exit traction
  endif

enddo traction
close(11)

deallocate(dshape_quad4)
deallocate(gll_weights,lagrange_gll,dlagrange_gll)
if(.not.trac_stat)then
  write(errtag,'(a)')'ERROR: all tractions cannot be read!'
  return
endif

errcode=0

return

end subroutine apply_traction
!===============================================================================
