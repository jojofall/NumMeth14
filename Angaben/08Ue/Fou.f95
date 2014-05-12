!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Einfacher Einschub für Fouriertransformation !
!                                              !
! Es wird die Groesse                          !
!   Ft(om)=(1/pi)*int_0^(dt*Nt) cos(om*t)*z(t) !
!                                              !
! Methode: Trapezregel                         !
!                                              !
! Annahme:                                     !
!   t,z(t) steht auf File 2 in folgender Form  !
!       0     z(0)                             !
!      dt    z(dt)                             !  
!     2*dt  z(2*dt)                            !
!     3*dt  z(3*dt)                            !
!      ...    ...                              !
!                                              !
! Ausgabe: auf file 3 = 'Fou.d'                !
!    om,Ft,Ft**2                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! fourier transformation
  open(3,file='Fou.d',form='formatted')
  rewind 3
  do io=1,No
    rewind 2
    om=io*0.2
    read(2,2000)z
    Ft=cos(om*z(0))*z(1)/2.d0
    do i=1,Nt
      read(2,2000)z
      Ft=Ft+cos(om*z(0))*z(1)
    end do
    Ft=Ft-cos(om*z(0))*z(1)/2.d0
    Ft=dt*Ft/pi
    write(3,2000)om,Ft,Ft*Ft
  end do
