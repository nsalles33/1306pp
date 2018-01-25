! =================================================================================================
!     Kinetic Monte Carlo Kernel
! -------------------------------------------------------------------------------------------------
!  Nicolas Salles, chercheur Precaire
!  Miha Gunde, future Toulousain
!  Code inspiration: Comput. Mehtods Appl. Mech. Engrg. 197 (2008) 3386-3398
! =================================================================================================
!  The idea its to write a kernel enough general to 
!  connect a wide range of problem.
!
!  System I want:
!  - Vacuncies diffusion : 2D and 3D
!  - Atomic deposition ( 3D )
!  - Flip atomic nature ( 2D ) 
!  
!  Modularity:
!  - read input file ( input_KMC.dat ) to Initialize the system
!  - connection with shared library (.so) to introduce other king of system... ( GRAAL )
!  - Accecpt differents kind of KMC algorithms ( Gillepsie & other )
! 
!  NOTE: 
!  - I try to design the program in object oriented... may be is not convinient... I try
!  
! =================================================================================================
  module derived_types
    implicit none

    private 
    public :: constructor, destructor, print_state

    type, public :: KMC_type
      logical                            :: bavard
      integer                            :: nsites
      real                               :: sum_rate, t
      integer, dimension(*), allocatable :: site
      real, dimension(*), allocatable    :: rate
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

  Program Kernel
!
  use derived_type
  use KMC_routine
  implicit none
!
  integer          :: nstep, ievent
  type( KMC_type ) :: sys 
!

! -=::: SYSTEM INITIALIZATION :::=-
  call Init_system( sys )

! -=::: SYSTEM EVOLUTION :::=-
  nstep = 0
  do while ( nstep <= sys% max_step )
     nstep = nstep + 1
     call print_state( sys )
     call rate_sum_calc( sys )
     ievent = choose_event( sys )
     call Event_Applied( sys, ievent )
     call Time_increment( sys )
  enddo
!
  call print_state( sys )

! -=::: ANALYSE & CONCLUSION :::=-
  call analyse( sys )
  call print_conclusion( sys )


! -=::: FINILIZATION :::=-
  call destructor( sys )

  end program Kernel











