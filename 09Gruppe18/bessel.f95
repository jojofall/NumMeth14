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
