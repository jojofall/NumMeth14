  module pendpar
!-----------------------------------------------------------
!  Module pendpar contains the parameters of the driven
!  pendulum
!----------------------------------------------------------- 
    double precision :: eta,om0,om2
    double precision :: fd,omd
  end module pendpar

  program exampleA
!
  use pendpar
!
  implicit none
!
  integer,parameter :: Nt=9000  ! Anzahl Zeitschritte
  integer,parameter :: No=500   ! Anzahl der Schritte in omega
  double precision :: dt        ! Zeitschritt in s
  double precision,parameter :: gg=9.80665 ! m/s**2  ! g
  double precision,parameter :: hh=6.626076e-34 ! Js ! h
!
  integer :: i,Npi,io,ne,m
  double precision :: Ft,om
  double precision :: y(0:2),z(0:2)
  double precision :: pi,lp,mp
  double precision :: hom0,vE,alf0,En,sum
!
  pi=4.*atan(1.d0)
!
! input of parameters of the pendulum
  lp=0.2484057 ! Länge des Pendels in m
  mp=1.e-2     ! Masse des Pendels in kg
  eta=2.00     ! 1/s
  fd=2.0       ! reduzierte Kraft 1/s**2
  omd= 5.5     ! Treiber Kreisfrequenz 1/s
  ne=0
! derived quantities
!  om2=gg/lp        ! Quadrat der Eigenkreisfrequenz
!  om0=sqrt(gg/lp)  ! Eigenkreisfrequenz
!  dt=pi/(om0*30.)       ! Schrittweite 1/60. der Periode in s
!  hom0=hh*om0/(2.*pi)   ! hbar*om0
!  vE=mp*lp**2/(2.*hom0)           ! 
!  alf0=sqrt((ne+0.5)/(vE*om0**2)) ! qmax maximale Auslenkung
