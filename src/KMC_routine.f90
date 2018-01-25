
  module KMC_routine
    use 
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






