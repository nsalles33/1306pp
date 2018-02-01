
! =================================================================================================
  module errors
    implicit none
  contains
! ............................................................................
      subroutine print_error( text )
        character (len=*), intent( in ) :: text
        print*, " ERROR: ",trim(text)
        stop " ===== FATAL ERROR!!"
      end subroutine print_error
! ............................................................................
      subroutine warning( text )
        character (len=*), intent( in ) :: text
        print*, " WARNING: ",trim(text)
      end subroutine warning
! ............................................................................
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
!    real, parameter :: kb = 1.38065e-12 ! J/K
    real, parameter :: kb = 8.61733e-7  ! eV/K
    public :: print_state

! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public :: event_type  ! This type will can be modulate by the user!!
      integer                              :: nevent
      integer, dimension(:), allocatable   :: init_state,   &
                                              final_state
      real, dimension(:), allocatable      :: ebarrier,     &
                                              de
     contains
      procedure :: builder => constructor_event_type 
      procedure :: destroy => destructor_event_type
      procedure :: printer => print_event
    end type event_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!
! '''''''''''''''''''''''''''' NEW TYPE '''''''''''''''''''''''''''''''''''''''
    type, public :: KMC_type
      logical                              :: bavard,      &
                                              period_x,    &
                                              period_y,    &
                                              period_z
      character (len=25)                   :: algorithm,   &
                                              input_file,  &
                                              input_event
      integer                              :: tot_sites,   &
                                              max_step,    &
                                              sys_dim,     &
                                              freq_write,  &
                                              nprop
      real                                 :: sum_rate,    &
                                              rand_rate,   &
                                              time,        &
                                              rand_time,   &
                                              max_time,    &
                                              temp, kt, f0
    !  character (len=25), dimension(:), allocatable :: txtprop
      integer, dimension(3)                         :: nsites
      integer, dimension(:), allocatable            :: site,        &
                                                       env_site,    &
                                                       nneig,       &
                                                       nevt
      integer, dimension(:,:), allocatable          :: neig,        &
                                                       event_site
      real, dimension(:), allocatable      :: rate,        &
                                              prop        
      real, dimension(:,:), allocatable    :: event_rate
      type( event_type )                   :: event
     contains
      procedure                            :: builder    => constructor_kmc_type
      procedure                            :: destroy    => destructor_kmc_type
      procedure                            :: printer    => print_kmc_type
      procedure                            :: init_table => init_table
    end type KMC_type
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

  CONTAINS ! %%%%%%%%%%%%%%%%%%%%%%%%%

! ............................................................................
! ............................ KMC_T BUILDER .................................
! ............................................................................
    subroutine constructor_kmc_type( this, dim_, size_x, size_y, size_z, step, time, algorithm, temperature )
      implicit none
      class( KMC_type ) :: this
      integer, intent( in ) :: size_x, dim_
      integer,           optional, intent( in ) :: size_y, size_z, step
      character (len=*), optional, intent( in ) :: algorithm 
      real,              optional, intent( in ) :: temperature, time

      integer :: n !,i,jn
      print*, " Enter in KMC_type constructor..."

      !this% period_x = .false.
      !this% period_y = .false.
      !this% period_z = .false.

      this% time = 0
      this% sum_rate = 0
      this% bavard = .false.
      this% algorithm = trim(algorithm)

!  ::: Temperature
!      print*, 'Temperature..'
      if ( present(temperature) ) then
         this% temp = temperature
      else
         this% temp = 300.0
      endif
      this% kt = 1.0 ;! kb*this% temp
      this% f0 = 1.0 ;! 1e-12

!  ::: max step & time
!      print*, 'Step  and Time...'
      this% max_step = 1
      this% max_time = 1.0
      if ( present(step) ) this% max_step = step
      if ( present(time) ) this% max_time = time

!  ::: System Size
!      print*, 'System Size...'
      if ( dim_ > 3 .or. dim_ == 0 ) &
         call print_error( " BAD SYSTEM DIMENSION " )
      this% sys_dim = dim_
      this% nsites(:) = 1
      this% nsites(1) = size_x

      if ( present(size_y) ) this% nsites(2) = size_y
      if ( .not.present(size_y).and.dim_ == 2) this% nsites(2) = size_x
      if ( present(size_z) ) this% nsites(3) = size_z

      if ( dim_ == 3.and.(.not.present(size_y)).and.(.not.present(size_z)) )  &
        call print_error( "System 3D => Ly and Lz must be declared " )

      this% tot_sites = this% nsites(1)*this% nsites(2)*this% nsites(3)
      n = this% tot_sites

!  ::: Table allocation 
!      print*, 'Table allocation...'
      allocate( this% site(n), this% rate(n), this% nneig(n), this% neig(10,n), this% event_rate(10,n),  &
                this% nevt(n), this% event_site(10,n) )
      if ( .not.allocated(this% site)  .or. .not.allocated(this% rate) .or.   &
           .not.allocated(this% nneig) .or. .not.allocated(this% neig) .or.   &
           .not.allocated(this% event_rate).or. .not.allocated(this% nevt) .or.&
           .not.allocated(this% event_site) )                                 &
         call print_error( " KMC_type => CONSTRUCTOR problem..." )
!
!  ::: Properties table
!      print*, 'Prop Table allocation...'
      if ( this% nprop /= 0 ) allocate( this% prop( this% nprop ) ) !, this% txtprop(this% nprop) )
!
      print*, " KMC_type Constructor DONE"

    end subroutine constructor_kmc_type
! ............................................................................
!
    subroutine destructor_kmc_type ( this )
      implicit none
      class( KMC_type ) :: this
!  :: destroy event_type
      call this% event% destroy  

      deallocate( this% site, this% rate, this% nneig, this% neig, this% event_rate,   &
                  this% nevt, this% event_site )

      if ( allocated(this% site).or.allocated(this% rate).or.       &
           allocated(this% nneig).or.allocated(this% neig).or.      &
           allocated(this% event_rate).or.allocated(this% nevt).or. &
           allocated(this% event_site) )   &
         print*, " KMC_type => DESTRUCTOR problem ..."

    end subroutine destructor_kmc_type
! ............................................................................
!
    subroutine Init_table( this )
      implicit none
      class( KMC_type ) :: this
      integer :: i, jn
!  ::: Rate Initialization at 0.
      do i = 1,this% tot_sites
         this% rate( i ) = 0.0
         do jn = 1,this% nneig( i )
            this% event_rate( jn, i ) = 0.0
         enddo
      enddo
!  ::: Properties Table initialization 
      do i = 1, this% nprop
         this% prop( i ) = 0.0
      enddo
    end subroutine Init_table
! ............................................................................
!
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
!
    subroutine destructor_event_type( this )
      implicit none
      class( event_type ) :: this
      deallocate( this% init_state, this% final_state, this% ebarrier, this% de)
      if ( allocated(this% init_state).or.allocated(this% final_state).or.  &
           allocated(this% ebarrier).or.allocated(this% de))                &
         call print_error( "EVENT_type => DESTRUCTOR problem " )

    end subroutine destructor_event_type
! ............................................................................
!
    subroutine print_kmc_type ( this )
      implicit none
      class( KMC_type ) :: this
      integer           :: i

!      print*, " - To do :: print_kmc_type "
      write (*,*) " ======== SYSTEM PARAMETERS ======= "
      write (*,*) " Algorithm         : ", trim(this% algorithm)
      write (*,*) " MAX KMC steps     : ", this% max_step
      write (*,*) " MAX KMC time (s?) : ", this% max_time
      write (*,*) " System Dimension  : ", this% sys_dim
      do i = 1,this%sys_dim
         write (*,*) " Nber of Node : ", this% nsites(i)
      enddo
      write (*,*) " System Size (nodes) : ", this% tot_sites
      write (*,*) " Freq Write          : ", this% freq_write
      write (*,*) " Bound Condition     : ", this% period_x, this% period_y, this% period_z
      write (*,*) " ================================== "

    end subroutine print_kmc_type
! ............................................................................
    subroutine print_event( this ) 
      implicit none
      class( event_type ) :: this
      integer :: i, id, nevt2

      nevt2 = this% nevent/2

      write (*,*) " =========== EVENT LIB ========= "
      do i = 1,nevt2
         id = i
         write (*,*) " Evt:",id,this% init_state( id ), this% final_state( id ), &
                                this% ebarrier( id ), this% de( id )
         id = i + nevt2
         write (*,*) " Inv:",id,this% init_state( id ), this% final_state( id ), &
                                this% ebarrier( id ), this% de( id )
      enddo

    end subroutine print_event
! ............................................................................
!
    subroutine print_state( this, step, u0 )
      implicit none
      type( KMC_type ) :: this
      integer, intent ( in ) :: step, u0
      integer :: ios, i, x, y, z, nx, ny, nxy
      
      if ( MODULO( step, this% freq_write ) /= 0 ) return 

      nx = this% nsites(1) 
      ny = this% nsites(2)
      nxy = nx*ny 

!  ::: Write Configuration
!    -- 3D cubic
      write (u0,fmt='(1x,I6)',iostat=ios) this% tot_sites
        if (ios /=0) write (*,*) " Problem write => state_file.xyz  "

      write (u0,*,iostat=ios) step, this% sys_dim,"DIM", (this% nsites(i),i=1,this% sys_dim)
        if (ios /=0) write (*,*) " Problem write => state_file.xyz  "

      do i = 0,this%tot_sites - 1
         x = MODULO( i, nx ) 
         y = MODULO( i/nx, ny )
         z = i / nxy
         write (u0,'(1x,I2,3(2x,I4))') this% site(i+1), x, y, z 
      enddo

!  ::: Write Energetic Statistic
      if ( step == 0 )  then
        write (*,*) "   Time   |  rand_rate    |    rand_time    " !,( " | ",trim(this% txtprop(i)),i=1,this%nprop )
        write (*,*) " ----------------------------------------------------------------- "
      endif
      write (*,'(4(2x,E10.4))') this% time, this% rand_rate, this% rand_time, this% prop(1) 
      write (100,'(4(2x,E10.4))') this% time, this% rand_rate, this% rand_time, this% prop(1) 

    end subroutine print_state

  end module derived_types
! ============================================================================
!
  module KMC_routine
    use derived_types
    use errors
    implicit none

  contains

! .......................................................
!    subroutine analyse( struc )
!      implicit none
!      type( KMC_type ) :: struc
!!      print*, " - To do :: analyse "
!    end subroutine analyse
! .......................................................
    subroutine print_conclusion( struc )
      implicit none
      type( KMC_type ) :: struc
!      print*, " - To do :: print_cnclusion "
    end subroutine print_conclusion

  end module KMC_routine






