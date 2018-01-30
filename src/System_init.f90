! =================================================================================================
!    Routine to initialize the system
!  To do (or think) :
!  - Format file for input_file and event library
!    
! =================================================================================================

    subroutine Init_system( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer :: i,id
      integer                           :: nsites, sys_dim, u0, ios, node_state, nevent, nstep
      real                              :: temp
      character (len=50)                :: input_file, string, input_event, algo_txt
      logical                           :: EOF

      character (len=1)  :: delims
      CHARACTER (len=100), dimension (50) :: args
      integer :: nargs

!   ...We must read input file
      print*, " - To do :: Init_system "
      print*, " - To do :: Read_input "

      sys_dim = 2
      nsites = 10

!  ::: Read input_file
      input_file = "input_KMC.dat"
      open( newunit=u0, file=trim(input_file), iostat=ios )
      if ( ios /= 0 ) call print_error( "input_file does't open!!" )
      !write (*,*) trim(input_file)," open...",u0,ios

      EOF = .false.
      delims = " "
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         if ( args(1) == "system_dimension" ) read ( args(2), '(i5)' ) sys_dim
         if ( args(1) == "Nber_node_1dim" )   read ( args(2), '(i5)' ) nsites
         if ( args(1) == "node_prop" )        read ( args(2), '(i5)' ) node_state
         if ( args(1) == "input_event" )      read ( args(2), '(a)'  ) input_event
         if ( args(1) == "algorithm" )        read ( args(2), '(a)'  ) algo_txt
         if ( args(1) == "temperature" )      read ( args(2), '(f6.6)' ) temp
         if ( args(1) == "nstep" )            read ( args(2), '(I6)' ) nstep
         if ( args(1) == "freq_write" )       read ( args(2), '(I6)' ) struc% freq_write
      enddo
      write (*,*) " LECTURE ",sys_dim, nsites, node_state, input_event
      close( u0 )

      call struc%builder( sys_dim, nsites, algorithm=algo_txt, temperature=temp, step=nstep )
      call struc%print_type

!  ::: Lecture of state and rate of each node
      open( newunit=u0, file=trim(input_event), iostat=ios )
      if ( ios /= 0 ) call print_error( "input_event does't open!!" )

      EOF = .false.
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)
         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call struc% event% builder( nevent )
            exit
         endif
      enddo
      do i = 1,struc% event% nevent
         read (u0,*) id, struc% event% init_state(id), struc% event% final_state(id),   &
                         struc% event% ebarrier(id), struc% event% de(id)
         write (*,*) id, struc% event% init_state(id), struc% event% final_state(id),   &
                         struc% event% ebarrier(id), struc% event% de(id)
      enddo

!  ::: Initialization of node State
!      This step depend on system...
      call distribution_state ( struc, 0.05 ) ! This routine is specific for vacancies diffusion
      call neig_list( struc )

    end subroutine Init_system
! =================================================================================================

    subroutine distribution_state ( struc, per100 )
      use derived_types
      use random
      implicit none 
      type( KMC_type ), intent( inout ) :: struc
      integer                           :: i
      real, intent( in )                :: per100
      real                              :: rnd

      call set_random_seed()

      do i = 1,struc% tot_sites
         call random_number( rnd ) 
         if ( rnd <= 1 - per100 ) then
            struc% site( i ) = 1
         else
            struc% site( i ) = 0
         endif
      enddo

    end subroutine distribution_state
! =================================================================================================

    subroutine neig_list( struc )
      use derived_types
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer :: i,j, x, nx, y, ny
      nx = struc% nsites(1)
      ny = struc% nsites(2)
!
!  :: System 2D 
      do i = 1,struc% tot_sites 
         struc% nneig( i ) = 4

         x = MODULO(i,nx)           ! ** begin at 1 because i begin at 1
         y = i / nx
         struc% neig(1,i) = i - 1
         struc% neig(2,i) = i + 1
         struc% neig(3,i) = i - nx
         struc% neig(4,i) = i + nx

!       -- Peridoc Bound Condition
         if (x == 1)      struc% neig(1,i) = i + (nx - 1)
         if (x == nx)     struc% neig(2,i) = i - (nx - 1)
         if (y == 0)      struc% neig(3,i) = i + ny*(nx - 1)
         if (y == ny - 1) struc% neig(4,i) = i - ny*(nx - 1) 

         if (i < 4) &
           write (*,*) i,x,y,nx, " neig ", (struc% neig(j,i),j=1,4)
      enddo

    end subroutine neig_list
! =================================================================================================



















