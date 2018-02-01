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

!  ::: Read input_file
      call read_input( struc )

      call struc% printer

      call read_event( struc )
      call struc% event% printer

!  ::: Initialization of node State
!      This step depend on system...
      call distribution_state ( struc, 0.0 ) ! This routine is specific for vacancies diffusion

      call neig_list( struc )

    end subroutine Init_system
! =================================================================================================

    subroutine read_input( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer                           :: ny,nz
      integer                           :: nsites, sys_dim, u0, ios, node_state, nstep
      real                              :: temp
      character (len=50)                :: string, input_event, algo_txt
      character (len=1)                 :: px,py,pz
      logical                           :: EOF

      character (len=1)  :: delims
      CHARACTER (len=100), dimension (50) :: args
      integer :: nargs

      sys_dim = 2
      nsites = 10
      struc% nprop = 0

      !input_file = "input_KMC.dat"
      open( newunit=u0, file=trim(struc% input_file), iostat=ios )
      if ( ios /= 0 ) call print_error( "input_file does't open!!" )
      !write (*,*) trim(input_file)," open...",u0,ios

      EOF = .false.
      delims = " "
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         if ( args(1) == "system_dimension" ) read ( args(2), '(i5)' ) sys_dim
         if ( args(1) == "Nber_node_1dim".or.  &
              args(1) == "nsite_x" )          read ( args(2), '(i5)' ) nsites
         if ( args(1) == "nsite_y" )          read ( args(2), '(i5)' ) ny
         if ( args(1) == "nsite_z" )          read ( args(2), '(i5)' ) nz

         if ( args(1) == "node_prop" )        read ( args(2), '(i5)' ) node_state
         if ( args(1) == "input_event" )      read ( args(2), '(a)'  ) struc% input_event
         if ( args(1) == "algorithm" )        read ( args(2), '(a)'  ) algo_txt
         if ( args(1) == "temperature" )      read ( args(2), '(f6.6)' ) temp
         if ( args(1) == "nstep" )            read ( args(2), '(I6)' ) nstep
         if ( args(1) == "freq_write" )       read ( args(2), '(I6)' ) struc% freq_write
         if ( args(1) == "calc_properties" )  read ( args(2), '(I6)' ) struc% nprop
         if ( args(1) == "bound_condition" )  then
            read ( args(2), '(a)' ) px
            if ( trim(px) == "p") struc% period_x = .true.
            read ( args(3), '(a)' ) py
            if ( trim(py) == "p") struc% period_y = .true.
            read ( args(4), '(a)' ) pz
            if ( trim(pz) == "p") struc% period_z = .true.
         endif
      enddo
      write (*,*) " LECTURE ",sys_dim, nsites, node_state, input_event
      close( u0 )

      if ( sys_dim == 2 )   &
        call struc%builder( sys_dim, nsites, algorithm=algo_txt, temperature=temp, step=nstep )
      if ( sys_dim == 3 )   &
        call struc%builder( sys_dim, nsites, size_y=ny, size_z=nz, algorithm=algo_txt, temperature=temp, step=nstep )

    end subroutine read_input
! =================================================================================================

!   subroutine distribution_state ( struc, per100 )
!     use derived_types
!     use random
!     implicit none 
!     type( KMC_type ), intent( inout ) :: struc
!     integer                           :: i
!     real, intent( in )                :: per100
!     real                              :: rnd

!     call set_random_seed()

!     do i = 1,struc% tot_sites
!        call random_number( rnd ) 
!        if ( rnd <= 1 - per100 ) then
!           struc% site( i ) = 1
!        else
!           struc% site( i ) = 0
!        endif
!     enddo

!   end subroutine distribution_state
! =================================================================================================

