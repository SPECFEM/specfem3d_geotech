! this module contains routines to compute Gauss-Legendre-Lobatto quadrature
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module gll_library
use set_precision

private :: endw1,endw2,gammaf,jacg

contains
! this subroutine computes GLL quadrature points and weights for 3D
subroutine gll_quadrature(ndim,ngllx,nglly,ngllz,ngll,gll_points,gll_weights,lagrange_gll,dlagrange_gll)
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll
real(kind=kreal),dimension(ngll),intent(out) :: gll_weights
real(kind=kreal),dimension(ndim,ngll),intent(out) :: gll_points
real(kind=kreal),dimension(ngll,ngll),intent(out) :: lagrange_gll
real(kind=kreal),dimension(ndim,ngll,ngll),intent(out) :: dlagrange_gll

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal
integer :: i,ii,j,k,n
real(kind=kreal) :: xi,eta,zeta
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx ! gll points and weights
real(kind=kreal),dimension(nglly) :: gllpy,gllwy ! gll points and weights
real(kind=kreal),dimension(ngllz) :: gllpz,gllwz ! gll points and weights
real(kind=kreal),dimension(ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(nglly) :: lagrange_y,lagrange_dy
real(kind=kreal),dimension(ngllz) :: lagrange_z,lagrange_dz

! compute everything in indexed order

! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz=ngll, need to call only once
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

n=0
do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx
      n=n+1
      ! integration points
      gll_points(1,n)=gllpx(i)
      gll_points(2,n)=gllpy(j)
      gll_points(3,n)=gllpz(k)

      ! integration weights
      gll_weights(n)=gllwx(i)*gllwy(j)*gllwz(k)
    enddo
  enddo
enddo

do ii=1,ngll ! ngllx*nglly*ngllz
  xi=gll_points(1,ii)
  eta=gll_points(2,ii)
  zeta=gll_points(3,ii)

  ! compute 1d lagrange polynomials on GLL points
  ! this can also be computed in a simple manner due to the orthogonality
  call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
  call lagrange1dGLL(nglly,gllpy,eta,lagrange_y,lagrange_dy)
  call lagrange1dGLL(ngllz,gllpz,zeta,lagrange_z,lagrange_dz)
  n=0
  do k=1,ngllz
    do j=1,nglly
      do i=1,ngllx
        n=n+1
        lagrange_gll(ii,n)=lagrange_x(i)*lagrange_y(j)*lagrange_z(k)
        dlagrange_gll(1,ii,n)=lagrange_dx(i)*lagrange_y(j)*lagrange_z(k)
        dlagrange_gll(2,ii,n)=lagrange_x(i)*lagrange_dy(j)*lagrange_z(k)
        dlagrange_gll(3,ii,n)=lagrange_x(i)*lagrange_y(j)*lagrange_dz(k)
      enddo
    enddo
  enddo
enddo
return
end subroutine gll_quadrature
!===========================================

! this subroutine computes GLL quadrature points and weights for 2D
subroutine gll_quadrature2d(ndim,ngllx,nglly,ngll,gll_points2d,gll_weights2d,lagrange_gll2d,dlagrange_gll2d)
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngll
real(kind=kreal),dimension(ngll),intent(out) :: gll_weights2d
real(kind=kreal),dimension(ndim,ngll),intent(out) :: gll_points2d
real(kind=kreal),dimension(ngll,ngll),intent(out) :: lagrange_gll2d
real(kind=kreal),dimension(ndim,ngll,ngll),intent(out) :: dlagrange_gll2d

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal
integer :: i,ii,j,n
real(kind=kreal) :: xi,eta !,zeta
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx ! gll points and weights
real(kind=kreal),dimension(nglly) :: gllpy,gllwy ! gll points and weights
real(kind=kreal),dimension(ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(nglly) :: lagrange_y,lagrange_dy

! compute everything in indexed order

! gll points and weights (source: http://mathworld.wolfram.com/lobattoquadrature.html)
!gllp(1)=-1.0_kreal   ; gllw(1)=1.0_kreal/3.0_kreal
!gllp(2)= 0.0_kreal   ; gllw(2)=4.0_kreal/3.0_kreal
!gllp(3)= 1.0_kreal   ; gllw(3)=gllw(1)

! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz=ngll, need to call only once
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)

n=0
do j=1,nglly
  do i=1,ngllx
    n=n+1
    ! integration points
    gll_points2d(1,n)=gllpx(i)
    gll_points2d(2,n)=gllpy(j)

    ! integration weights
    gll_weights2d(n)=gllwx(i)*gllwy(j)
  enddo
enddo

do ii=1,ngll ! ngllx*nglly
  xi=gll_points2d(1,ii)
  eta=gll_points2d(2,ii)

  ! compute 1d lagrange polynomials
  ! this can also be computed in a simple manner due to the orthogonality
  call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
  call lagrange1dGLL(nglly,gllpy,eta,lagrange_y,lagrange_dy)
  n=0
  do j=1,nglly
    do i=1,ngllx
      n=n+1
      lagrange_gll2d(ii,n)=lagrange_x(i)*lagrange_y(j)
      dlagrange_gll2d(1,ii,n)=lagrange_dx(i)*lagrange_y(j)
      dlagrange_gll2d(2,ii,n)=lagrange_x(i)*lagrange_dy(j)
    enddo
  enddo
enddo

return
end subroutine gll_quadrature2d
!===========================================

! this subroutine computes GLL quadrature points and weights for 1D
subroutine gll_quadrature1d(ndim,ngllx,ngll,gll_points1d,gll_weights1d,lagrange_gll1d,dlagrange_gll1d)
implicit none
integer,intent(in) :: ndim,ngllx,ngll
real(kind=kreal),dimension(ngll),intent(out) :: gll_weights1d
real(kind=kreal),dimension(ndim,ngll),intent(out) :: gll_points1d
real(kind=kreal),dimension(ngll,ngll),intent(out) :: lagrange_gll1d
real(kind=kreal),dimension(ndim,ngll,ngll),intent(out) :: dlagrange_gll1d

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal
integer :: i,ii,n
real(kind=kreal) :: xi
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx ! gll points and weights
real(kind=kreal),dimension(ngllx) :: lagrange_x,lagrange_dx

! compute everything in indexed order

! gll points and weights (source: http://mathworld.wolfram.com/lobattoquadrature.html)
!gllp(1)=-1.0_kreal   ; gllw(1)=1.0_kreal/3.0_kreal
!gllp(2)= 0.0_kreal   ; gllw(2)=4.0_kreal/3.0_kreal
!gllp(3)= 1.0_kreal   ; gllw(3)=gllw(1)

! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz=ngll, need to call only once
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)

n=0
do i=1,ngllx
  n=n+1
  ! integration points
  gll_points1d(1,n)=gllpx(i)

  ! integration weights
  gll_weights1d(n)=gllwx(i)
enddo

do ii=1,ngll ! ngllx
  xi=gll_points1d(1,ii)

  ! compute 1d lagrange polynomials
  ! this can also be computed in a simple manner due to the orthogonality
  call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)

  n=0
  do i=1,ngllx
    n=n+1
    lagrange_gll1d(ii,n)=lagrange_x(i)
    dlagrange_gll1d(1,ii,n)=lagrange_dx(i)
  enddo
enddo

return
end subroutine gll_quadrature1d
!===========================================

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
subroutine lagrange1dGEN(nenod,xi,phi,dphi_dxi)
implicit none
integer,intent(in) :: nenod ! number of nodes in an 1d element
integer :: i,j,k
real(kind=kreal),intent(in) :: xi ! point where to calculate lagrange function and
!its derivative
real(kind=kreal),dimension(nenod),intent(out) :: phi,dphi_dxi
real(kind=kreal),dimension(nenod) :: xii,term,dterm,sum_term
real(kind=kreal) :: dx

! compute natural coordnates
dx=2.0_kreal/real((nenod-1),kreal)! length = 2.0 as xi is taken -1 to +1
do i=1,nenod
  ! coordinates when origin is in the left
  xii(i)=real((i-1),kreal)*dx
enddo

! origin is tranformed to mid point
xii=xii-1.0_kreal

do i=1,nenod
  k=0
  phi(i)=1.0_kreal
  do j=1,nenod
    if(j/=i)then
      k=k+1
      term(k)=(xi-xii(j))/(xii(i)-xii(j))
      dterm(k)=1.0_kreal/(xii(i)-xii(j)) ! derivative of the term wrt xi

      phi(i)=phi(i)*(xi-xii(j))/(xii(i)-xii(j))
    endif
  enddo

  sum_term=1.0_kreal
  do j=1,nenod-1
    do k=1,nenod-1
      if(k==j)then
        sum_term(j)=sum_term(j)*dterm(k)
      else
        sum_term(j)=sum_term(j)*term(k)
      endif
    enddo
  enddo
  dphi_dxi(i)=0.0_kreal
  do j=1,nenod-1
    dphi_dxi(i)=dphi_dxi(i)+sum_term(j)
  enddo
enddo

return
end subroutine lagrange1dGEN
!===========================================

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
subroutine lagrange1dGLL(nenod,xii,xi,phi,dphi_dxi)
implicit none
integer,intent(in) :: nenod ! number of nodes in an 1d element
real(kind=kreal),dimension(nenod),intent(in) :: xii
real(kind=kreal),intent(in) :: xi ! point where to calculate lagrange function and
!its derivative
real(kind=kreal),dimension(nenod),intent(out) :: phi,dphi_dxi

integer :: i,j,k
real(kind=kreal),dimension(nenod) :: term,dterm,sum_term

do i=1,nenod
  k=0
  phi(i)=1.0_kreal
  do j=1,nenod
    if(j/=i)then
      k=k+1
      term(k)=(xi-xii(j))/(xii(i)-xii(j))
      dterm(k)=1.0_kreal/(xii(i)-xii(j)) ! derivative of the term wrt xi

      phi(i)=phi(i)*(xi-xii(j))/(xii(i)-xii(j))
    endif
  enddo

  sum_term=1.0_kreal
  do j=1,nenod-1
    do k=1,nenod-1
      if(k==j)then
        sum_term(j)=sum_term(j)*dterm(k)
      else
        sum_term(j)=sum_term(j)*term(k)
      endif
    enddo
  enddo
  dphi_dxi(i)=0.0_kreal
  do j=1,nenod-1
    dphi_dxi(i)=dphi_dxi(i)+sum_term(j)
  enddo
enddo

return
end subroutine lagrange1dGLL
!===========================================

!===========================================
!
!  Library to compute the Gauss-Lobatto-Legendre points and weights
!  Based on Gauss-Lobatto routines from M.I.T.
!  Department of Mechanical Engineering
!
!===========================================

real(kind=kreal) function endw1(n,alpha,beta) !double precision

implicit none

integer n
real(kind=kreal) alpha,beta !double precision

real(kind=kreal), parameter :: zero=0._kreal,one=1._kreal,two=2._kreal,three=3._kreal,four=4._kreal !double precision
real(kind=kreal) apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3 !double precision
!double precision, external :: gammaf
integer i

f3 = zero
apb   = alpha+beta
if (n == 0) then
  endw1 = zero
  return
endif
f1   = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
f1   = f1*(apb+two)*two**(apb+two)/two
if (n == 1) then
  endw1 = f1
  return
endif
fint1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
fint1 = fint1*two**(apb+two)
fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
fint2 = fint2*two**(apb+three)
f2    = (-two*(beta+two)*fint1 + (apb+four)*fint2) * (apb+three)/four
if (n == 2) then
  endw1 = f2
  return
endif
do i=3,n
  di   = dble(i-1)
  abn  = alpha+beta+di
  abnn = abn+di
  a1   = -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
  a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
  a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
  f3   =  -(a2*f2+a1*f1)/a3
  f1   = f2
  f2   = f3
enddo
endw1  = f3

end function endw1
!=======================================================================

real(kind=kreal) function endw2(n,alpha,beta) !double precision

implicit none

integer n
real(kind=kreal) alpha,beta !double precision

real(kind=kreal), parameter :: zero=0._kreal,one=1._kreal,two=2._kreal,three=3._kreal,four=4._kreal !double precision
real(kind=kreal) apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3 !double precision
!real(kind=kreal), external :: gammaf
integer i

apb   = alpha+beta
f3 = zero
if (n == 0) then
  endw2 = zero
  return
endif
f1   = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
f1   = f1*(apb+two)*two**(apb+two)/two
if (n == 1) then
  endw2 = f1
  return
endif
fint1 = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
fint1 = fint1*two**(apb+two)
fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
fint2 = fint2*two**(apb+three)
f2    = (two*(alpha+two)*fint1 - (apb+four)*fint2) * (apb+three)/four
if (n == 2) then
  endw2 = f2
  return
endif
do i=3,n
  di   = dble(i-1)
  abn  = alpha+beta+di
  abnn = abn+di
  a1   =  -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
  a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
  a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
  f3   =  -(a2*f2+a1*f1)/a3
  f1   = f2
  f2   = f3
enddo
endw2  = f3

end function endw2

!
!=======================================================================
!

real(kind=kreal) function gammaf (x) !double precision

implicit none

real(kind=kreal), parameter :: pi = 3.141592653589793_kreal !double precision

real(kind=kreal) x !double precision

real(kind=kreal), parameter :: half=0.5_kreal,one=1._kreal,two=2._kreal !double precision

gammaf = one

if (x == -half) gammaf = -two*sqrt(pi)
if (x ==  half) gammaf =  sqrt(pi)
if (x ==  one ) gammaf =  one
if (x ==  two ) gammaf =  one
if (x ==  1.5_kreal) gammaf =  sqrt(pi)/2._kreal
if (x ==  2.5_kreal) gammaf =  1.5_kreal*sqrt(pi)/2._kreal
if (x ==  3.5_kreal) gammaf =  2.5_kreal*1.5_kreal*sqrt(pi)/2._kreal
if (x ==  3._kreal ) gammaf =  2._kreal
if (x ==  4._kreal ) gammaf = 6._kreal
if (x ==  5._kreal ) gammaf = 24._kreal
if (x ==  6._kreal ) gammaf = 120._kreal

end function gammaf

!
!=====================================================================
!

subroutine jacg (xjac,np,alpha,beta)

!=======================================================================
!
! computes np Gauss points, which are the zeros of the
! Jacobi polynomial with parameters alpha and beta
!
!                  .alpha = beta =  0.0  ->  Legendre points
!                  .alpha = beta = -0.5  ->  Chebyshev points
!
!=======================================================================

implicit none

integer np
real(kind=kreal) alpha,beta !double precision
real(kind=kreal) xjac(np) !double precision

integer k,j,i,jmin,jm,n
real(kind=kreal) xlast,dth,x,x1,x2,recsum,delx,xmin,swap !double precision
real(kind=kreal) p,pd,pm1,pdm1,pm2,pdm2 !double precision

integer, parameter :: K_MAX_ITER = 10
real(kind=kreal), parameter :: zero = 0._kreal, eps = 1.0e-12_kreal !double precision

pm1 = zero
pm2 = zero
pdm1 = zero
pdm2 = zero

xlast = 0._kreal
n   = np-1
dth = 4._kreal*atan(1._kreal)/(2._kreal*dble(n)+2._kreal)
p = 0._kreal
pd = 0._kreal
jmin = 0
do j=1,np
  if(j == 1) then
    x = cos((2._kreal*(dble(j)-1._kreal)+1._kreal)*dth)
  else
    x1 = cos((2._kreal*(dble(j)-1._kreal)+1._kreal)*dth)
    x2 = xlast
    x  = (x1+x2)/2._kreal
  endif
  do k=1,K_MAX_ITER
    call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
    recsum = 0._kreal
    jm = j-1
    do i=1,jm
        recsum = recsum+1._kreal/(x-xjac(np-i+1))
    enddo
    delx = -p/(pd-recsum*p)
    x    = x+delx
    if(abs(delx) < eps) goto 31
  enddo
31  continue
  xjac(np-j+1) = x
  xlast        = x
enddo
do i=1,np
  xmin = 2._kreal
  do j=i,np
    if(xjac(j) < xmin) then
        xmin = xjac(j)
        jmin = j
    endif
  enddo
  if(jmin /= i) then
    swap = xjac(i)
    xjac(i) = xjac(jmin)
    xjac(jmin) = swap
  endif
enddo

end subroutine jacg

!
!=====================================================================
!

subroutine jacobf (poly,pder,polym1,pderm1,polym2,pderm2,n,alp,bet,x)

!=======================================================================
!
! Computes the Jacobi polynomial of degree n and its derivative at x
!
!=======================================================================

implicit none

real(kind=kreal) poly,pder,polym1,pderm1,polym2,pderm2,alp,bet,x !double precision
integer n

real(kind=kreal) apb,polyl,pderl,dk,a1,a2,b3,a3,a4,polyn,pdern,psave,pdsave !double precision
integer k

apb  = alp+bet
poly = 1._kreal
pder = 0._kreal
psave = 0._kreal
pdsave = 0._kreal

if (n == 0) return

polyl = poly
pderl = pder
poly  = (alp-bet+(apb+2._kreal)*x)/2._kreal
pder  = (apb+2._kreal)/2._kreal
if (n == 1) return

do k=2,n
  dk = dble(k)
  a1 = 2._kreal*dk*(dk+apb)*(2._kreal*dk+apb-2._kreal)
  a2 = (2._kreal*dk+apb-1._kreal)*(alp**2-bet**2)
  b3 = (2._kreal*dk+apb-2._kreal)
  a3 = b3*(b3+1._kreal)*(b3+2._kreal)
  a4 = 2._kreal*(dk+alp-1._kreal)*(dk+bet-1._kreal)*(2._kreal*dk+apb)
  polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
  pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
  psave  = polyl
  pdsave = pderl
  polyl  = poly
  poly   = polyn
  pderl  = pder
  pder   = pdern
enddo

polym1 = polyl
pderm1 = pderl
polym2 = psave
pderm2 = pdsave

end subroutine jacobf

!
!------------------------------------------------------------------------
!

real(kind=kreal) FUNCTION PNDLEG (Z,N) !double precision

!------------------------------------------------------------------------
!
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!------------------------------------------------------------------------
implicit none

real(kind=kreal) z !double precision
integer n

real(kind=kreal) P1,P2,P1D,P2D,P3D,FK,P3 !double precision
integer k

P1   = 1._kreal
P2   = Z
P1D  = 0._kreal
P2D  = 1._kreal
P3D  = 1._kreal

do K = 1, N-1
  FK  = dble(K)
  P3  = ((2._kreal*FK+1._kreal)*Z*P2 - FK*P1)/(FK+1._kreal)
  P3D = ((2._kreal*FK+1._kreal)*P2 + (2._kreal*FK+1._kreal)*Z*P2D - FK*P1D) / (FK+1._kreal)
  P1  = P2
  P2  = P3
  P1D = P2D
  P2D = P3D
enddo

PNDLEG = P3D

end function pndleg

!
!------------------------------------------------------------------------
!

real(kind=kreal) FUNCTION PNLEG (Z,N) !double precision

!------------------------------------------------------------------------
!
!     Compute the value of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!------------------------------------------------------------------------
implicit none

real(kind=kreal) z !double precision
integer n

real(kind=kreal) P1,P2,P3,FK !double precision
integer k

P1   = 1._kreal
P2   = Z
P3   = P2

do K = 1, N-1
  FK  = dble(K)
  P3  = ((2._kreal*FK+1._kreal)*Z*P2 - FK*P1)/(FK+1._kreal)
  P1  = P2
  P2  = P3
enddo

PNLEG = P3

end function pnleg

!
!------------------------------------------------------------------------
!

real(kind=kreal) function pnormj (n,alpha,beta) !double precision

implicit none

real(kind=kreal) alpha,beta !double precision
integer n

real(kind=kreal) one,two,dn,const,prod,dindx,frac !double precision
!real(kind=kreal), external :: gammaf
integer i

one   = 1._kreal
two   = 2._kreal
dn    = dble(n)
const = alpha+beta+one

if (n <= 1) then
  prod   = gammaf(dn+alpha)*gammaf(dn+beta)
  prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
  pnormj = prod * two**const/(two*dn+const)
  return
endif

prod  = gammaf(alpha+one)*gammaf(beta+one)
prod  = prod/(two*(one+const)*gammaf(const+one))
prod  = prod*(one+alpha)*(two+alpha)
prod  = prod*(one+beta)*(two+beta)

do i=3,n
  dindx = dble(i)
  frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
  prod  = prod*frac
enddo

pnormj = prod * two**const/(two*dn+const)

end function pnormj

!
!------------------------------------------------------------------------
!

subroutine zwgjd(z,w,np,alpha,beta)

!=======================================================================
!
!     Z w g j d : Generate np Gauss-Jacobi points and weights
!                 associated with Jacobi polynomial of degree n = np-1
!
!     Note : Coefficients alpha and beta must be greater than -1.
!     ----
!=======================================================================

implicit none

real(kind=kreal), parameter :: zero=0._kreal,one=1._kreal,two=2._kreal !double precision

integer np
real(kind=kreal) z(np),w(np) !double precision
real(kind=kreal) alpha,beta !double precision

integer n,np1,np2,i
real(kind=kreal) p,pd,pm1,pdm1,pm2,pdm2 !double precision
real(kind=kreal) apb,dnp1,dnp2,fac1,fac2,fac3,fnorm,rcoef !double precision
!real(kind=kreal), external :: gammaf,pnormj

pd = zero
pm1 = zero
pm2 = zero
pdm1 = zero
pdm2 = zero

n    = np-1
apb  = alpha+beta
p    = zero
pdm1 = zero

if(np <= 0)then
  write(*,*)'ERROR: number of Gauss points < 1!'
  stop
endif

if((alpha <= -one) .or. (beta <= -one))then
  write(*,*)'ERROR: alpha and beta must be greater than -1!'
  stop
endif

if (np == 1) then
  z(1) = (beta-alpha)/(apb+two)
  w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two) * two**(apb+one)
  return
endif

call jacg(z,np,alpha,beta)

np1   = n+1
np2   = n+2
dnp1  = dble(np1)
dnp2  = dble(np2)
fac1  = dnp1+alpha+beta+one
fac2  = fac1+dnp1
fac3  = fac2+one
fnorm = pnormj(np1,alpha,beta)
rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
do i=1,np
  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
  w(i) = -rcoef/(p*pdm1)
enddo

end subroutine zwgjd

!
!------------------------------------------------------------------------
!

subroutine zwgljd(z,w,np,alpha,beta)

!=======================================================================
!
!     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
!     -----------   weights associated with Jacobi polynomials of degree
!                   n = np-1.
!
!     Note : alpha and beta coefficients must be greater than -1.
!            Legendre polynomials are special case of Jacobi polynomials
!            just by setting alpha and beta to 0.
!
!=======================================================================

implicit none

real(kind=kreal), parameter :: zero=0._kreal,one=1._kreal,two=2._kreal !double precision

integer np
real(kind=kreal) alpha,beta !double precision
real(kind=kreal) z(np), w(np) !double precision

integer n,nm1,i
real(kind=kreal) p,pd,pm1,pdm1,pm2,pdm2 !double precision
real(kind=kreal) alpg,betg !double precision
!real(kind=kreal), external :: endw1,endw2

p = zero
pm1 = zero
pm2 = zero
pdm1 = zero
pdm2 = zero

n   = np-1
nm1 = n-1
pd  = zero

if(np <= 1)then
  write(*,*)'ERROR: number of Gauss-Lobatto points < 2!'
  stop
endif

! with spectral elements, use at least 3 points
if(np < 3)then
  write(*,*)'WARNING: number of Gauss-Lobatto points < 3!'
  !stop
endif
!if (np <= 2) stop 'minimum number of Gauss-Lobatto points for the SEM is 3'

if((alpha <= -one) .or. (beta <= -one))then
  write(*,*)'ERROR: alpha and beta must be greater than -1!'
  stop
endif

if (nm1 > 0) then
  alpg  = alpha+one
  betg  = beta+one
  call zwgjd(z(2),w(2),nm1,alpg,betg)
endif

z(1)  = - one
z(np) =  one

do i=2,np-1
  w(i) = w(i)/(one-z(i)**2)
enddo

call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
w(1)  = endw1(n,alpha,beta)/(two*pd)
call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
w(np) = endw2(n,alpha,beta)/(two*pd)

end subroutine zwgljd
end module gll_library


