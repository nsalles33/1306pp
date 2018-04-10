
! =================================================================================================

  module hidden_table
    use iso_c_binding
    implicit none
    integer( c_int ), dimension(:), target, allocatable :: h_i_state,    &
                                                           h_f_state,    &  
                                                           h_site,       &
                                                           h_nneig,      &
                                                           h_nevt,       &
                                                           h_neig,       &
                                                           h_event_site, &
                                                           h_spec
    !
    real( c_double ), dimension(:), target, allocatable :: h_rate,       &
                                                           h_prop,       &
                                                           h_f0,         &
                                                           h_ebarrier,   &
                                                           h_de,         &
                                                           h_event_rate, &
                                                           h_ppressure,  &
                                                           h_masse
    !
    real( c_double ), dimension(:,:), target, allocatable :: h_ebond  
    !
  end module hidden_table
! =================================================================================================
  !
  module derived_types
    use iso_c_binding
    use hidden_table
    !use errors
    implicit none

    !private 
    integer( c_int ), parameter :: nvois = 10
    !real( c_double ), parameter :: kb = 1.38065e-23 ! J/K
    real( c_double ), parameter :: kb = 8.61733e-5     ! eV/K
    real( c_double ), parameter :: na = 6.02214e23     ! atomes 
    real( c_double ), parameter :: pi = 4.0*atan(1.0)  ! 3.14...
    real( c_double ), parameter :: u = 1.66e-27        ! kg
    real( c_double ), parameter :: ang = 1e-10         ! m 
    !public :: print_state

! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public, bind( C ) :: event_type  ! This type will can be modulate by the user!!
      integer( c_int ) :: nevent, &
                          nbond,  &
                          nchem_react
      !
      type( c_ptr )    :: ptr_i_state,  &
                          ptr_f_state,  &
                          ptr_f0,       &
                          ptr_ebarrier, &
                          ptr_de,       &
                          ptr_ebond
    end type event_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!
! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public, bind( C ) :: KMC_type
      integer( c_int )                               :: bavard,      &
                                                        conv      
      integer( c_int ), dimension( 3 )               :: period,      &
                                                        nsites
      !
      character ( len=50,kind=c_char )               :: algorithm,   &
                                                        input_file,  &
                                                        input_event, &
                                                        libname,     &
                                                        init_mod
      !
      integer( c_int )                               :: tot_sites,   &
                                                        max_step,    &
                                                        sys_dim,     &
                                                        freq_write,  &
                                                        nprop,       &
                                                        node_state,  &
                                                        nspec,       &
                                                        npressure
      !
      real( c_double )                                :: sum_rate,    &
                                                         rand_rate,   &
                                                         time,        &
                                                         rand_time,   &
                                                         max_time,    &
                                                         temp, kt,    &
                                                         per100, f0,  &
                                                         scale,       &
                                                         step
      !
      type( c_ptr )  ::  ptr_site, ptr_nneig, ptr_nevt, ptr_neig, ptr_event_site,    &
                         ptr_spec, ptr_rate, ptr_prop, ptr_event_rate, ptr_pressure, &
                         ptr_masse
      !
      type( event_type )                   :: event
      !
    end type KMC_type
    !
  end module derived_types
! =================================================================================================
