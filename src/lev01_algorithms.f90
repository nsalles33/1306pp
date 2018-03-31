!
!   All Kinetic Monte Carlo Algorithms Available
!   - Residence Time Algorithms ( BKL )
!   - Gillepsie
! =================================================================================================

  subroutine Algorithm( struc )
    use iso_c_binding
    use derived_types
    implicit none
    !
    type( KMC_type ) :: struc
    character (len=3) :: temp1
    character (len=9) :: temp2
    !
    !
    temp1 = trim(struc% algorithm)
    !if ( trim(struc% algorithm) == 'BKL' ) call algo_bkl( struc )
    if ( temp1 == 'BKL' ) call algo_bkl( struc )
    !
    temp2 = trim(struc% algorithm)
    !if ( trim(struc% algorithm) == 'Gillepsie' ) call algo_Gillepsie( struc )
    if ( temp2 == 'Gillepsie' ) call algo_Gillepsie( struc )
    !
  end subroutine Algorithm
! =================================================================================================

  subroutine algo_bkl( struc )
    use derived_types
    use lib_hub
    implicit none

    type( KMC_type ) :: struc
    integer :: isite, ievent

    !print*, " -= BKL Algorithm =- "

    call rate_sum_calc( struc )

    call choose_event_hub( struc, isite, ievent )

    if ( struc% conv == 1 ) return

    call event_applied_hub( struc, isite, ievent )

    call Time_increment( struc )
    !stop " in algo_BKL..."

  end subroutine algo_bkl
! .................................................................................................

    subroutine time_increment( struc )
      use iso_c_binding
      use derived_types
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      !character (len = 3, kind=c_char)  :: tmp1
      !character (len = 9, kind=c_char)  :: tmp2
      real( c_double )                  :: rdn
      !
      !tmp1 = struc% algorithm
      !tmp2 = struc% algorithm
      !print*, " - To do :: Time_increment ", struc% algorithm, tmp1,tmp2
      !
      !if ( trim(struc% algorithm) == "BKL" ) then
         call random_number( rdn )
         struc% rand_time = - log( rdn )/struc% sum_rate
         struc% time = struc% time + struc% rand_time
      !else if ( trim(struc% algorithm) == "Gillepsie" ) then
      !else if ( tmp2 == "Gillepsie" ) then
      !   print*, " time_increment",struc% rand_time
      !   struc% time = struc% time + struc% rand_time
      !endif
      !
      !stop " Time_increment..."
      !
    end subroutine
! =================================================================================================

  subroutine algo_gillepsie( struc )
    use derived_types
    use lib_hub
    implicit none
    !
    type( KMC_type ) :: struc
    integer :: ievent, isite
    !
    call rate_sum_calc( struc )
    !
    call choose_event_hub( struc, isite, ievent )
    !
    if ( struc% conv == 1 ) return
    !
    call event_applied_hub( struc, isite, ievent )
    !
    !call Time_increment( struc )
    struc% time = struc% time + struc% rand_time
    !
  end subroutine algo_gillepsie
! =================================================================================================
