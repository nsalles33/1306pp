!
!    Routine to calculate the rate of each event and each site  
!
!
    subroutine rate_sum_calc( struc )
      use iso_c_binding
      use derived_types

#ifdef SHARED_LIB
      use dlopen_lib
#endif

      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int )                  :: i
      real( c_double )                  :: sum_rate

      real( c_double ), dimension(:), pointer :: rate

#ifdef SHARED_LIB
      procedure( event_rate_calc ), bind( C ), pointer :: event_rate_calc_proc

      call c_f_procpointer( proc_event_rate_calc, event_rate_calc_proc )
      call event_rate_calc_proc( struc )
#else
      call event_rate_calc( struc )
#endif
!
!  -- Sum over all the site
!
      sum_rate = 0
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      do i = 1,struc% tot_sites
         sum_rate = sum_rate + rate( i )
      enddo
      struc% sum_rate = sum_rate
!
!      write (*,*) " DONE "
!
    end subroutine rate_sum_calc
! ..................................................................................................





