module ising2d
   implicit none
   ! Turn on or off the periodic boundary conditions
   logical, parameter :: periodic = .true.
   ! Shape of the array small_spins(-V_BOUND:V_BOUND,-H_BOUND,H_BOUND). This corresponds to 
   ! the number of neighbours to each side (H_bound) and above/below (V_BOUND) the site (i,j)
   ! you want to include in your interaction term. If you only want to have nearest neighbours
   ! these values should be set to 1. See the function get_econf for the use of the small_spins array.
   integer,parameter                :: V_BOUND=1,H_BOUND=1                         

   ! The smallest kind of intergers
   integer,parameter                :: short = SELECTED_INT_KIND(1)                
   ! all_spins is an array that stores the spins (the periodicity of this array is toggled by the periodic flag above). 
   ! The periodicity is implemented as an extra boundary layer around your spin box. Make sure 
   ! that you use the subroutine set_spin when you want to update the spins in this array, as this
   ! subroutine automatically takes care of the periodicity.
   integer(KIND=short), allocatable :: all_spins(:,:)                              
   ! array storing the initial spin configuration. To be used in the MRI exercise.
   integer(KIND=short), allocatable :: initial_spins(:,:)                           
   ! Precomputed energy differences and acceptance probabilities maped by the energy index econf.
   real,allocatable                 :: cached_delta_e(:), cached_accept_prob(:)    

contains
   subroutine ising2d_init(width,height,temperature,total_econfs)
      !< Allocates the (global) spin, and the (private) delta_es and accepts arrays. 
      !< Initializes the all_spins array to 0.
      !< Precomputes the energy differences between spin configurations. 
      !< The energy differences are labeled by an energy index of your choice. (See get_econf)
      !< Precomputes the acceptance probability given the energy index by calling the function accept_prob.
      integer,intent(in)      :: width,height   ! The height and the width of the spin box
      real(KIND=8),intent(in) :: temperature    ! The temperature of the system
      integer,intent(in)      :: total_econfs   ! The total number of energy configurations (See get_econf)
      integer                 :: er             ! Error flag
      integer                 :: econf          ! loop variable (runs over the energy configurations)
      real(KIND=8)            :: delta_etmp     ! Temporary variable for the energy difference

      ! The periodicity is implemented by padding the borders of the all_spins array with extra spins that are 
      ! copies of the spins in all_spins(1:width,1:height). Note that one should use the subroutine set_spin
      ! to update the all_spins array!
      allocate(all_spins(-V_BOUND+1:height+V_BOUND,-H_BOUND+1:width+H_BOUND),&
               initial_spins(-V_BOUND+1:height+V_BOUND,-H_BOUND+1:width+H_BOUND),&
               cached_delta_e(total_econfs),cached_accept_prob(total_econfs),stat=er)
      if (er .ne. 0) stop 'Error: Could not allocate the spin array'
      all_spins = 0
      cached_delta_e = 0
      cached_accept_prob = 0
      
      ! Precompute the energy differences and the acceptance probability
      do econf = 1,total_econfs
         ! Calculate the energy difference for energy configuration econf.
         delta_etmp = delta_e(econf)
         ! Store the energy difference in cached_delta_e 
         cached_delta_e(econf) = delta_etmp
         ! Calculate the corresponding acceptance probability
         cached_accept_prob(econf) = accept_prob(temperature,delta_etmp)
      enddo

      return
   end subroutine
   
   subroutine ising2d_end()
      !< Deallocates the arrays allocated by ising_init for a clean exit.
      deallocate(all_spins,initial_spins,cached_delta_e,cached_accept_prob)
      return
   end subroutine

   integer pure function get_econf(small_spins,newspin)
   !< Maps the spin configuration and the trial spin to an energy index
   ! small_spins is a matrix of the neighbouring spins. Note that the center of small_spins is at (0,0).
   ! The shape of the small_spins array is by default (-1:1,-1:1), i.e. a 3x3 matrix.
   integer(KIND=short),intent(in) :: small_spins(-V_BOUND:V_BOUND,-H_BOUND:H_BOUND)  
   integer(KIND=short),intent(in) :: newspin   ! The new (trial) value for the spin in the center of small_spins.

      ! #####################################################################################
      ! 0) Here you should add a few lines of code which 
      !    *quickly* maps the energy difference between the old and the new spin configuration
      !    to a (unique) energy index. You will have to come up with a good way to 
      !    index the energies yourself based on the form of the spin interaction.
      !    Whenever you update this routine, make sure that you also update the way you 
      !    calculate the total number of energy configurations (total_econfs).
      !    NOTE 1: The spins outside of your box/image will have the value 0 
      !          (unless you have periodic boundary conditions)
      !    NOTE 2: Make sure that the energy index starts at 1 (e.g. by using an appropriate offset).
      ! #####################################################################################

      get_econf = 

      return
   end function

   real(KIND=8) pure function delta_e(econf)
   !< Calculates the energy difference corresponding to the energy index econf (given by the function get_econf)
      integer,intent(in) :: econf 

      ! #####################################################################################
      ! 1) Here you should calculate the energy difference corresponding to your energy index.
      !    These energies are cached, so speed is not an issue.
      ! #####################################################################################

      delta_e = 

      return
   end function

   real pure function accept_prob(temperature,delta_e)
      real(KIND=8), intent(in) :: temperature
      real(KIND=8), intent(in) :: delta_e      ! Change in energy

      ! #####################################################################################
      ! 2) Here you should calculate the acceptance probability given the temperature and  
      !    the change in energy. The probabilities are cached, so speed is not an issue.
      ! #####################################################################################

      accept_prob = 

      return
   end function

   subroutine set_spin(i,j,newspin,height,width)
      !< Set the spin with index (i,j) to newspin, honouring periodic boundary conditions
      integer, intent(in)              :: i,j
      integer(KIND=short), intent(in)  :: newspin
      integer, intent(in)              :: height, width   ! The height and the width of the spin box
      
      ! Change the spin in the all_spins array
      all_spins(i,j) = newspin

      ! Make sure that the copies at the boundaries of all_spins are also updated
      if (periodic) then
         if(i .le. V_BOUND) then
             all_spins(i+height,j) = newspin
             if(j .le. H_BOUND) all_spins(height+i,j+width) = newspin
             if(j .gt. width-H_BOUND) all_spins(i+height,j-width) = newspin
         elseif(i .gt. height-V_BOUND) then
             all_spins(i-height,j) = newspin
             if(j .le. H_BOUND) all_spins(i-height,j+width) = newspin
             if(j .gt. width-H_BOUND) all_spins(i-height,j-width) = newspin
         endif
         if(j .le. H_BOUND) then
             all_spins(i,j+width) = newspin
         elseif(j .gt. width-H_BOUND) then
             all_spins(i,j-width) = newspin
         endif
      endif

      return
   end subroutine

   subroutine check(spin_array,height,width)
      integer(KIND=short),intent(in) :: spin_array(-V_BOUND+1:width+V_BOUND,-H_BOUND+1:height+H_BOUND)
      integer,intent(in)             :: height, width ! Height and width of the spin box

      if (periodic) then
         ! Checks that the spin board is consistent
         if(any(spin_array(1-V_BOUND:0,:).ne.spin_array(height+1-V_BOUND:height,:))) then
            stop "INTERNAL ERROR: NORTH BORDER"
         elseif(any(spin_array(height+1:height+V_BOUND,:).ne.spin_array(1:V_BOUND,:))) then
            stop "INTERNAL ERROR: SOUTH BORDER"
         elseif(any(spin_array(:,1-H_BOUND:0).ne.spin_array(:,width+1-H_BOUND:width))) then
            stop "INTERNAL ERROR: WEST BORDER"
         elseif(any(spin_array(:,width+1:width+H_BOUND).ne.spin_array(:,1:H_BOUND))) then
            stop "INTERNAL ERROR: EAST BORDER"
         endif
      else
         ! Checks that the spin board is set to zero
         if(any(spin_array(1-V_BOUND:0,:).ne. 0 )) then
            stop "INTERNAL ERROR: NORTH BORDER"
         elseif(any(spin_array(:,1-H_BOUND:0).ne.0)) then
            stop "INTERNAL ERROR: WEST BORDER"
         elseif(any(spin_array(height+1:height+V_BOUND,:).ne.0)) then
            stop "INTERNAL ERROR: SOUTH BORDER"
         elseif(any(spin_array(:,width+1:width+H_BOUND).ne.0)) then
            stop "INTERNAL ERROR: EAST BORDER"
         endif
      endif

      return
   end subroutine

   subroutine print_all_spins(ounit,info,height,width,clear)
        !< Prints the current spin configuration to the console for debugging
        integer, intent(in) :: ounit,height,width
        real,intent(in)     :: info(:)
        logical, intent(in) :: clear
        character, parameter :: CHARS(-1:4) = (/ '.', ' ', 'O', 'x', 'o', 'M' /)
        integer :: i,j,nrows

        ! We are actually only printing a reasonable part of the all_spins array.
        nrows = min(height+2*V_BOUND,23)
        ! Write an ansi escape character to clear the screen
        if(clear) write(ounit,"(A,'[',I2.2,'A')",advance="no") char(27),nrows+1
        ! Write the all_spins array, including the boundary!
        do i = 1-V_BOUND,nrows-V_BOUND
            do j = 1-H_BOUND,min(width+H_BOUND,60)
                write (ounit,"(A1)",advance="no") CHARS(all_spins(i,j))
            enddo
            write (ounit,*)
        enddo
        ! Write some extra information
        write(ounit,*) "Extra info:",info
        return
    end subroutine

end module

program ising_mc
   use ising2d
   implicit none
   integer,parameter    :: ounit = 9, iunit = 10         ! Set the units for output and input
   integer              :: warmups                       ! Number of warmup steps *per spin*
   integer              :: measurements                  ! Number of measurement steps *per spin*
   integer              :: width                         ! Width of the spin box
   integer              :: height                        ! Height of the spin box
   integer(KIND=short),allocatable  :: small_spins(:,:)  ! Array containing the neighbouring spins
   real(KIND=8)         :: temperature                   ! Temperature of the system
   logical              :: read_from_file                ! Read the initial spin configuration from file

   integer,allocatable  :: i_rnds(:), j_rnds(:)          ! Arrays containing random positions
   real,allocatable     :: p_rnds(:)                     ! Array containing the Monte Carlo dice values
   real                 :: p                             ! temporary variable for the Monte Carlo dice
   integer              :: total_econfs                  ! Total number of energy configurations
   integer              :: econf                         ! Energy configuration index
   integer              :: i, j, k, l                    ! loop counters
   integer              :: imeasure, itry, ioffset       ! more loop counters
   integer(KIND=short)  :: newspin                       ! temporary variable containing the value of the trial spin
   integer              :: system_size                   ! system size = width*height

   integer              :: er                            ! Error flag
   ! Input
   integer           :: nargs                            ! Number of input arguments
   character(len=32) :: sval                             ! Temporary string 
   character(len=32) :: filename                         ! File name string

   ! Open the output unit
   open(ounit,file="out")

   ! Initialize: Read the input form the input arguments
   !(You could also read this directly from an input file)
   nargs = command_argument_count()
   call get_command_argument(1, sval)
   read (sval,*) height
   call get_command_argument(2, sval)
   read (sval,*) width
   call get_command_argument(3, sval)
   read (sval,*) temperature
   call get_command_argument(4, sval)
   read (sval,*) warmups      ! Warmup steps *per spin*
   call get_command_argument(5, sval)
   read (sval,*) measurements ! Measurement steps *per spin*
   read_from_file = nargs .ge. 6
   if (read_from_file) call get_command_argument(6, filename)

   ! Determine the system size
   system_size = height*width

   ! Write some info to the output
   write(ounit,"('Grid size:',I4,' by',I4,', temperature:',F9.3)") width, height, temperature
   write (ounit,*) ""

   ! #####################################################################################
   ! 2) Here you should calculate the total number of (unique) energy configurations. 
   !    In other words, total_econfs should be set to the largest index that your function 
   !    get_econf will return.
   ! #####################################################################################
   total_econfs = 

   ! Initialize the ising2d module
   call ising2d_init(width,height,temperature,total_econfs)

   ! Read the starting configuration from file.
   if (read_from_file) then
      ! Read the starting configuration from file
      open(iunit,file=TRIM(ADJUSTL(filename)),status='old')
      do i = 1,height
         do j= 1,width
            read(iunit,*) k,l,newspin
            if (i .ne. k .or. j .ne. l) stop 'Error: File format error! Check your input file.'
            call set_spin(i,j,newspin,height,width) ! Store the values in the all_spins array
         enddo
         read(iunit,*)
      enddo
   else
      ! Initialize all the spins to spin up
      do i = 1,height
         do j= 1,width
            call set_spin(i,j,1_short,height,width) ! Store the value in the all_spins array
         enddo
      enddo
   endif
   ! Store the configuration for future use in the MRI exercise (initial_spins is stored in the ising2d module)
   initial_spins = all_spins
   ! Check that no simple mistakes were made..
   call check(initial_spins,height,width)

   ! Allocate the small temporary spin configuration array, and the array for
   ! the preconputed random numbers and positions.
   allocate(small_spins(-V_BOUND:V_BOUND,-H_BOUND:H_BOUND),&
            p_rnds((warmups+measurements)*system_size),&
            i_rnds((warmups+measurements)*system_size),&
            j_rnds((warmups+measurements)*system_size),stat=er)
   if (er .ne. 0) stop 'Error: Could not allocate p_rnds'

   ! Precompute the random sites (i,j) and the Monte Carlo dice (p_rnds) for speed
   call random_number(p_rnds) ! Fills the array p_rnds with uniformly distributed random numbers in the interval [0,1]
   ! Calculate i and j
   ! #####################################################################################
   ! 4) MRI optional: Try to roughly exclude the trial points (i_rnds,j_rnds) outside of the skull. 
   !    Remember to also set the intensity of these points to zero in the all_spins and
   !    initial_spins array (using the set_spin subroutine).
   ! #####################################################################################
   p_rnds = p_rnds*width
   j_rnds = int(p_rnds)+1
   p_rnds = (p_rnds - (j_rnds -1))*height
   i_rnds = int(p_rnds)+1
   ! Get new random numbers for the dice
   call random_number(p_rnds)
   ! #####################################################################################
   ! 5) MRI: Rewrite the acceptance condition to the form: f(p)*T > delta_e1 + delta_e2
   !         where delta_e1 is the energy difference from the first term and delta_e2 is 
   !         the energy difference from the second term.
   !         Precompute f(p) here and store it in p_rnds, i.e. p_rnds = f(p_rnds).
   ! #####################################################################################
   


   ! Warmup loop
   do itry = 1,warmups*system_size
      ! Get the new random position and Monte Carlo dice
      i = i_rnds(itry)
      j = j_rnds(itry)
      p = p_rnds(itry)

      ! Copy the neighbouring spins to the small_spins array
      small_spins = all_spins(i-V_BOUND:i+V_BOUND,j-H_BOUND:j+H_BOUND)

      ! #####################################################################################
      ! 5) Here you should write a simple warmup loop. 
      !    Hint: Use get_econf and cached_accept_prob.
      !    MRI: You also need to include the intensity (z) from the initial_spins array.
      !         Hint: Use get_econf and cached_delta_e for delta_e1, and compute delta_e2 as
      !               delta_e2 = z^2*A + z*B + C, where A, B, and C are small (precomputed) arrays.
      !               Use the acceptance condition you obtained at point 5.
      !    Make sure that you use the subroutine set_spin when you update the spin configuration.
      ! #####################################################################################


   enddo


   ! Measurement/SA loop
   avgfactor = 1d0/measurements
   ioffset = 1+warmups*system_size ! Number of random numbers already used in the warmup loop
   ! Measure your configurations "measurements" times.
   do imeasure = 1,measurements
      ! Make the inner loop run over the system size 
      do itry = ioffset+(imeasure-1)*system_size,ioffset+imeasure*system_size
         ! Get the new random position and Monte Carlo dice
         i = i_rnds(itry)
         j = j_rnds(itry)
         p = p_rnds(itry)

         ! Copy the neighbouring spins to the small_spins array
         small_spins = all_spins(i-V_BOUND:i+V_BOUND,j-H_BOUND:j+H_BOUND)

         ! #####################################################################################
         ! 6) Here you should write the loop in which updates your spins.
         !    Ising: You also need to *update* the current magnetization when a flip is accepted.
         !    MRI: You also need to include the intensity (z) from the initial_spins array.
         !         Hint: Use get_econf and cached_delta_e for delta_e1, and compute delta_e2 as
         !               delta_e2 = z^2*A + z*B + C, where A, B, and C are small (precomputed) arrays.
         !               Use the acceptance condition you obtained at point 5.
         !    Make sure that you use the subroutine set_spin when you update the spin configuration.
         ! #####################################################################################

      enddo
      ! #####################################################################################
      ! 7) Ising: Use the current magnetization to construct the expectation values of the binder cumulant.
      !    Watch out for overflow errors (convert to dble before you apply the exponent)!
      !    MRI: Here you can adjust your temperature according to your SA scheme
      ! #####################################################################################

   enddo
   ! Check that no simple mistakes were made..
   call check(all_spins,height,width)

   ! Open the data output unit
   close(ounit)
   open(ounit,file="ising2d.dat",status='unknown',position='append')
   ! #####################################################################################
   ! 8) Here you should append your new data to the ising2d.dat file
   ! #####################################################################################
   write(ounit,*) 

   close(ounit)
   call ising2d_end()

   return
end program
