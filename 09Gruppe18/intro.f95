  module potpar
  double precision :: V0,R0,a0
  double precision :: V1,R1,a1
  end module
  
  program levinson

  use potpar

  implicit none
  integer*4, parameter :: Nrmax=2000
  integer*4, parameter :: NEmax=2000
  double precision, parameter :: mn=939.56563
  double precision, parameter :: amu=931.49432
  double precision, parameter :: hbarc=197.327053

  integer*4 :: l,Nr,NE,nb,ie,i,j,ir,x1,x2
  double precision :: A,mA,Anorm,Emax,dE,Rmax,rr
  double precision :: mred,mhb2,dr,E,k2,r,pi
  double precision :: bes1,neu1,bes2,neu2
  double precision :: y(0:Nrmax),delta(0:NEmax),bs(0:Nrmax,5)
  double precision :: z(0:Nrmax)
  double precision :: Vpot,omega,k
  double precision :: E1,cc1,E2,cc2,Ec,cc
!
  pi=4.d0*atan(1.d0)

! read integration parameters
  open(1,file='data.inp',form='formatted')
  read(1,*) A         ! mass number
  read(1,*) mA        ! rel. atomic mass of target nucleus
  read(1,*) l         ! orbital angular momentum quantum number
  read(1,*) Anorm     ! start normalization of wave function
  read(1,*) Emax,dE   ! maximal energy, energy step in MeV
  read(1,*) Rmax,Nr   ! maximum integration in fm, number of mesh points
  read(1,*) V0,rr,a0  ! SW potential parameters in MeV, fm, fm
  read(1,*) V1,R1,a1  ! parameters of additional Gauss potential in MeV,fm,fm**2
!
! derived quantities
  mred=mn*mA*amu/(mA*amu+mn)  ! reduced mass in MeV
  mhb2=2.d0*mred/hbarc**2   ! 2.*mred/hbarc**2 in fm**-2
  dr=min(Rmax/Nr,dfloat(Nrmax)) ! radiale Schrittweite in fm
  R0=rr*A**(1.d0/3.d0)      ! half density radius of potential in fm module potpar

!---------------------------------------------------------------------------------------
!Startwerte deklarieren
z(0)=0
z(1)=0.1
x1=1980
x2=1981
j=0

do ie=2000,1,-1
 do i=1,Nrmax-1
  z(i+1)=2*z(i)-z(i-1)-(dr**2)*omega(i*dr,l,dble(ie*de),mhb2)*z(i)
 end do
!
 do i=0,Nrmax
  y(i)=(1+((dr**2)/12)*omega(i*dr,l,dble(ie*de),mhb2))*z(i)
!  write(*,*) i*dr, y(i),Vpot(i*dr)
 end do

k=sqrt(ie*de*mhb2)
 call bessel(k*x1*dr,l,bes1,neu1)
 call bessel(k*x2*dr,l,bes2,neu2)
 delta(ie)=datan(-(y(x1)*bes2-y(x2)*bes1)/(y(x1)*neu2-y(x2)*neu1))
  if(ie .lt. 1998) then
  if(delta(ie+1)-delta(ie) .gt. 0.75*(4.*datan(1.d0))) then
  j=j+1
  end if
  end if
write(*,*) ie*de,delta(ie)+j*4.*atan(1.d0)
end do

end

!---------------------------------------------------------------------------------------
  double precision function omega(x,l,E,mhb2)
!
! determines the omega function at x

  use potpar

  implicit none
  double precision :: x,E,mhb2,vpot
  integer*4 :: l
  
  Vpot=-V0/(1.d0+exp((x-R0)/a0))
  Vpot=Vpot+V1*exp(-(x-R1)**2/a1)

  omega = -(l*(l+1))/(x**2)+(E-Vpot)*mhb2

  return
  end

!-------------------------------------------------------------------------------------
  subroutine bessel(x,l,bes,neu)
!
! evaluates the Bessel and Neumann functions in 
! Riccati form for given l at the argument x 
!   j0=sin(x),n0=cos(x)
!
  implicit none
!
  integer*4 :: l
  double precision :: x,bes,neu
!
  bes=0.d0
  neu=0.d0
  if(l.eq.0) then
    bes=sin(x)
    neu=cos(x)
  end if
!
  if(l.eq.1) then
    bes=sin(x)/x-cos(x)
    neu=cos(x)/x+sin(x)
  end if
!
  if(l.eq.2) then
    bes=sin(x)*(3/x**2-1.d0)-3.d0*cos(x)/x
    neu=cos(x)*(3/x**2-1.d0)+3.d0*sin(x)/x
  end if 
!
  return
  end

!------------------------------------------------------------------------------------
 double precision function Vpot(x)
!
! determines the potential at x

  use potpar

  implicit none
  double precision :: x

  Vpot=-V0/(1.d0+exp((x-R0)/a0))
  Vpot=Vpot+V1*exp(-(x-R1)**2/a1)
  return
  end

