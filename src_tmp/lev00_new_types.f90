
! =================================================================================================

  module hidden_table
    use iso_c_binding
    implicit none
    integer( c_int ), dimension(:), target, allocatable :: h_i_state, h_f_state, h_site, h_nneig, h_nevt
    integer( c_int ), dimension(:,:), target, allocatable :: h_neig, h_event_site
    real( c_double ), dimension(:), target, allocatable :: h_rate, h_prop, h_ebarrier, h_de
    real( c_double ), dimension(:,:), target, allocatable :: h_event_rate, h_ebond  
    !
  end module hidden_table
! =================================================================================================

  module derived_types
    use iso_c_binding
    use hidden_table
!    use errors
    implicit none

!    private 
!    real, parameter :: kb = 1.38065e-12 ! J/K
    real( c_double ), parameter :: kb = 8.61733e-7  ! eV/K
!    public :: print_state

! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public, bind( C ) :: event_type  ! This type will can be modulate by the user!!
      integer( c_int ) :: nevent, nbond
      type( c_ptr )    :: ptr_i_state,  &
                          ptr_f_state,  &
                          ptr_ebarrier, &
                          ptr_de,       &
                          ptr_ebond
    end type event_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!
! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public, bind( C ) :: KMC_type
      integer( c_int )                               :: bavard      
      integer( c_int ), dimension( 3 )               :: period,      &
                                                        nsites
      character ( len=50,kind=c_char )               :: algorithm,   &
                                                        input_file,  &
                                                        input_event, &
                                                        libname
      integer( c_int )                               :: tot_sites,   &
                                                        max_step,    &
                                                        sys_dim,     &
                                                        freq_write,  &
                                                        nprop,       &
                                                        node_state
      real( c_double )                                :: sum_rate,    &
                                                        rand_rate,   &
                                                        time,        &
                                                        rand_time,   &
                                                        max_time,    &
                                                        temp, kt,    &
                                                        per100, f0
    !  character (len=25), dimension(:), allocatable :: txtprop
      !integer( c_int ), dimension(3)                 :: nsites
      !
      type( c_ptr )  ::  ptr_site, ptr_nneig, ptr_nevt, ptr_neig, ptr_event_site, &
                         ptr_rate, ptr_prop, ptr_event_rate
      !
      type( event_type )                   :: event
      !
    end type KMC_type
    !
  end module derived_types
! =================================================================================================
