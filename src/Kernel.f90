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

  Program Kernel
!
  use derived_types
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
  !call print_conclusion( sys )


! -=::: FINILIZATION :::=-
  call sys% destroy

  end program Kernel











