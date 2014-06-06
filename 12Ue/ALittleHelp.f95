!dimesnion of the image
  integer Nimax,Njmax,
  parameter (Nimax=254,Njmax=333)


!read in image
  do i=1,Nimax
     do j=1,Njmax
        read(*,'(3I5)') dummyi,dummyj,MR(i,j)
     enddo
     read (*,*)
  enddo

!write out image
   do i=1,Nimax
      do j=1,Njmax
            write(32,'(3I5)') i,j,mag(i,j)
      enddo
      write (32,*)
   enddo

!mean values and standard deviation
   integer MW(5)
!                  BG  WM  GM  CSF  SB
   parameter (MW=(/30,426,602,1223,167/))
   integer SIGMA(5)
   parameter (SIGMA=(/30,59,102,307,69/))
 
!random tissue in 1..5:
        newmag=INT(rand()*5.d0)+1 

!new boundary condiition (period possible since BG->BG)
         !/* Check periodic boundary conditions */
              if (Inew .le. 0) then
                 Inew = Nimax
              else 
                 if(Inew .gt. Nimax ) then 
                    Inew = 1
                 endif
              endif
              if (Jnew .le. 0) then
                 Jnew = Njmax
              else 
                 if(Jnew .gt. Njmax) then 
                    Jnew = 1
                 endif
              endif

!check for floating point overflow
           if ((Enew-Eold)/T.gt.20.d0) then
              r=0.d0
           else
              if ((Enew-Eold)/T.lt.0.05d0) then
                 r=1.d0
              else
                 r=exp(-(Enew-Eold)/T)
              end if
           end if
