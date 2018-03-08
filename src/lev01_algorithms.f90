!
!   All Kinetic Monte Carlo Algorithms Available
!   - Residence Time Algorithms ( BKL )
!   - Gillepsie
! =================================================================================================

  subroutine Algorithm( struc )
    use derived_types
    implicit none

    type( KMC_type ) :: struc
    
    if ( trim(struc% algorithm) == 'BKL' ) call algo_BKL( struc )
    if ( trim(struc% algorithm) == 'Gillepsie' ) call algo_Gillepsie( struc )

  end subroutine Algorithm
! =================================================================================================

  subroutine algo_bkl( struc )
    use derived_types
    use KMC_routine

#ifdef SHARED_LIB
    use dlopen_lib
#endif
    implicit none

    type( KMC_type ) :: struc
    integer :: isite, ievent
#ifdef SHARED_LIB
    procedure( choose_event ), bind( C ), pointer :: choose_event_proc
    procedure( event_applied ), bind( C ), pointer :: event_applied_proc
#endif

    !print*, " -= BKL Algorithm =- "

    call rate_sum_calc( struc )

#ifdef SHARED_LIB
    call c_f_procpointer( proc_choose_event, choose_event_proc )
    call choose_event_proc( struc, isite, ievent )
    call c_f_procpointer( proc_event_applied, event_applied_proc )
    call event_applied_proc( struc, isite, ievent )
#else
    call choose_event( struc, isite, ievent )
    call event_applied( struc, isite, ievent )
#endif
    call Time_increment( struc )

  end subroutine algo_bkl
! .................................................................................................

    subroutine time_increment( struc )
      use iso_c_binding
      use derived_types
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      real( c_double )                  :: rdn

!      print*, " - To do :: Time_increment "

      call random_number( rdn )
      struc% rand_time = - log( rdn )/struc% sum_rate
      struc% time = struc% time + struc% rand_time

    end subroutine
! =================================================================================================

  subroutine algo_gillepsie( struc )
    use derived_types
    use KMC_routine
#ifdef SHARED_LIB
    use dlopen_lib
#endif
    implicit none

    type( KMC_type ) :: struc
    integer :: ievent, isite
#ifdef SHARED_LIB
    procedure( choose_event ), bind( C ), pointer :: choose_event_proc
    procedure( event_applied ), bind( C ), pointer :: event_applied_proc
#endif

#ifdef SHARED_LIB
    call c_f_procpointer( proc_choose_event, choose_event_proc )
    call choose_event_proc( struc, isite, ievent )
    call c_f_procpointer( proc_event_applied, event_applied_proc )
    call event_applied_proc( struc, isite, ievent )
#else
    call choose_event( struc, isite, ievent )
    call event_applied( struc, isite, ievent )
#endif
    call Time_increment( struc )

  end subroutine algo_gillepsie
! =================================================================================================
