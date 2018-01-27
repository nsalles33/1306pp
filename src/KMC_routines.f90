
! =================================================================================================
  module errors
    implicit none
  contains
      subroutine print_error( text )
        character (len=*), intent( in ) :: text
        print*, trim(text)
        stop " ===== FATAL ERROR!!"
      end subroutine print_error
  end module errors
! =================================================================================================
  module random
    implicit none
  contains
! ............................................................................
    subroutine set_random_seed()
     ! ----- setting a random seed, based on current time -----
     INTEGER :: i_seed
     INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
     INTEGER, DIMENSION(1:8) :: dt_seed

     ! ----- Set up random seed portably -----
     CALL RANDOM_SEED(size=i_seed)
     ALLOCATE(a_seed(1:i_seed))
     CALL RANDOM_SEED(get=a_seed)
     CALL DATE_AND_TIME(values=dt_seed)
     a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
     CALL RANDOM_SEED(put=a_seed)
     DEALLOCATE(a_seed)
    end subroutine set_random_seed
! ............................................................................

  end module random
! =================================================================================================
  module derived_types
    use errors
    implicit none

    private 
    public :: print_state

! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public :: event_type  ! This type will can be modulate by the user!!
      integer                              :: nevent
      integer, dimension(:), allocatable   :: init_state, final_state, ebarrier, de
     contains
      procedure :: builder => constructor_event_type 
      procedure :: destroy => destructor_event_type
    end type event_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public :: KMC_type
      logical                            :: bavard
      integer                            :: tot_sites,   &
                                            max_step,    &
                                            sys_dim
      real                               :: sum_rate,    &
                                            time,        &
                                            max_time
      integer, dimension(3)              :: nsites
      integer, dimension(:), allocatable :: site
      real, dimension(:), allocatable    :: rate
      type( event_type )                 :: event
     contains
      procedure                          :: builder    => constructor_kmc_type
      procedure                          :: destroy    => destructor_kmc_type
      procedure                          :: print_type => print_kmc_type
    end type KMC_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

  CONTAINS ! %%%%%%%%%%%%%%%%%%%%%%%%%

! ............................................................................
    subroutine constructor_kmc_type( this, dim_, size1_, size2_, size3_, step_, time_ )
      implicit none
      class( KMC_type ) :: this
      integer, intent( in ) :: size1_, dim_
      integer, optional, intent( in ) :: size2_, size3_, step_, time_
      print*, " Enter in KMC_type constructor..."
      this% time = 0
      this% sum_rate = 0
      this% bavard = .false.

!  ::: max step & time
      this% max_step = 1
      this% max_time = 1
      if ( present(step_) ) this% max_step = step_
      if ( present(time_) ) this% max_time = time_

!  ::: System Size
      if ( dim_ > 3 .or. dim_ == 0 ) &
         call print_error( " BAD SYSTEM DIMENSION " )
      this% sys_dim = dim_
      this% nsites(:) = 1
      this% nsites(1) = size1_
      if ( present(size2_) ) this% nsites(2) = size2_
      if ( .not.present(size2_).and.dim_ == 2) this% nsites(2) = size1_
      if ( present(size3_) ) this% nsites(3) = size3_
      if ( dim_ == 3.and.(.not.present(size2_)).and.(.not.present(size3_)) )  &
        call print_error( "System 3D => Ly and Lz must be declared " )
      this% tot_sites = this% nsites(1)*this% nsites(2)*this% nsites(3)

!  ::: Table allocation 
      allocate( this% site(this% tot_sites), this% rate(this% tot_sites) )
      if ( .not.allocated(this% site) ) call print_error( " KMC_type => CONSTRUCTOR problem..." )
      if ( .not.allocated(this% rate) ) call print_error( " KMC_type => CONSTRUCTOR problem..." )

      print*, " KMC_type Constructor DONE"

    end subroutine constructor_kmc_type
! ............................................................................

    subroutine destructor_kmc_type ( this )
      implicit none
      class( KMC_type ) :: this
!  :: destroy event_type
      call this% event% destroy  
      deallocate( this% site, this% rate )
      if ( allocated(this% site) ) print*, " KMC_type => DESTRUCTOR problem ..."
      if ( allocated(this% rate) ) print*, " KMC_type => DESTRUCTOR problem ..."
    end subroutine destructor_kmc_type
! ............................................................................

    subroutine constructor_event_type( this, n )
      implicit none
      class( event_type )   :: this
      integer, intent( in ) :: n
      print*, " Enter in EVENT_type constructor..."
      this% nevent = n
      allocate( this% Init_state(n), this% Final_state(n), this% Ebarrier(n), this% dE(n) )
      if ( .not.allocated(this% init_state) .or. .not.allocated(this% final_state) .or.   &
           .not.allocated(this% ebarrier) .or. .not.allocated(this% de) )                 &
         call print_error( " EVENT_type => CONSTRUCTOR problem " )
      print*, " EVENT_type constructor DONE "
    end subroutine constructor_event_type
! ............................................................................

    subroutine destructor_event_type( this )
      implicit none
      class( event_type ) :: this
      deallocate( this% init_state, this% final_state, this% ebarrier, this% de)
      if ( allocated(this% init_state).or.allocated(this% final_state).or.  &
           allocated(this% ebarrier).or.allocated(this% de))                &
         call print_error( "EVENT_type => DESTRUCTOR problem " )

    end subroutine destructor_event_type
! ............................................................................

    subroutine print_kmc_type ( this )
      implicit none
      class( KMC_type ) :: this
      integer           :: i

      print*, " - To do :: print_kmc_type "
      write (*,*) " ======== SYSTEM PARAMETERS ======= "
      write (*,*) " MAX KMC steps     : ", this% max_step
      write (*,*) " MAX KMC time (s?) : ", this% max_time
      write (*,*) " System Dimension  : ", this% sys_dim
      do i = 1,this%sys_dim
         write (*,*) " Nber of Node : ", this% nsites(i)
      enddo
      write (*,*) " System Size (nodes) : ", this% tot_sites
      write (*,*) " ================================== "

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
    use errors
    implicit none


  contains


    subroutine rate_sum_calc( struc )
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
    subroutine event_applied( struc, event )
      implicit none
      type( KMC_type ) :: struc
      integer, intent( in ) :: event
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
! .......................................................
    subroutine print_conclusion( struc )
      implicit none
      type( KMC_type ) :: struc
      print*, " - To do :: print_cnclusion "
    end subroutine print_conclusion

  end module KMC_routine






