
  module derived_types
    implicit none

    private
    public :: constructor, destructor, print_state

    type, public :: KMC_type
      logical                            :: bavard
      integer                            :: nsites
      real                               :: sum_rate, t
      integer, dimension(:), allocatable :: site
      real, dimension(:), allocatable    :: rate
     contains
      procedure :: print => print_kmc_type
    end type KMC_type

    contains

    subroutine constructor( this, a )
      implicit none
      type( KMC_type ) :: this
      integer, intent( in ) :: a
      this% t = 0
      this% sum_rate = 0
      this% nsites = a
      this% bavard = .false.
      allocate( this% site(a), this% rate(a) )
      if ( .not.allocated(this% site) ) print*, " CONSTRUCTOR problem..."
      if ( .not.allocated(this% rate) ) print*, " CONSTRUCTOR problem..."
    end subroutine constructor
! ............................................................................
    subroutine destructor ( this )
      implicit none
      type( KMC_type ) :: this
      deallocate( this% site, this% rate )
      if ( allocated(this% site) ) print*, " DESTRUCTOR problem ..."
      if ( allocated(this% rate) ) print*, " DESTRUCTOR problem ..."
    end subroutine destructor
! ............................................................................
    subroutine print_kmc_type ( this )
      implicit none
      type( KMC_type ) :: this
      print*, " - To do :: print_kmc_type "
    end subroutine print_kmc_type
! ............................................................................
    subroutine print_state( this )
      implicit none
      type( KMC_type ) :: this
      print*, " - To do :: print_state "
    end subroutine print_state
  end module derived_types
! ============================================================================

  module KMC_routine
    use derived_types
    implicit none


  contains

    subroutine Init_system( struc )
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      integer :: nsite
!   ...We must read input file
      print*, " - To do :: Init_system "
      print*, " - To do :: Read_input "
      nsite = 10
      call constructor( struc, nsite )
    end subroutine Init_system( struc )
! .......................................................
    subroutine rate_sum_calc( struct )
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      print*, " - To do :: rate_sum_calc "
    end subroutine rate_sum_calc
! .......................................................
    integer function choose_event( struc )
      implicit none
      type( KMC_type ) :: struc
      print*, " - To do :: choose_event "
    end function choose_event
! .......................................................
    subroutine event_applied( struc )
      implicit none
      type( KMC_type ) :: struc
      print*, " - To do :: Event_applied "
    end subroutine event_applied
! .......................................................
    subroutine time_increment( struc )
      implicit none
      type( KMC_type ) :: struc
      print*, " - To do :: analyse "
    end subroutine
! .......................................................
    subroutine analyse( struc )
      implicit none
      type( KMC_type ) :: struc
      print*, " - To do :: analyse "
    end subroutine analyse

  end module KMC_routine






