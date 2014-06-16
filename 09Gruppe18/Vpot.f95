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


!-----------------------------------------


