

    subroutine read_event( obj )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: obj

      character (len=50)                :: string
      integer                           :: u0, i, id, ios, nevent
      logical                           :: EOF

      character (len=1)  :: delims
      CHARACTER (len=100), dimension (50) :: args
      integer :: nargs

      !  ::: Lecture of state and rate of each node
      open( newunit=u0, file=trim(obj% input_event), iostat=ios )
      if ( ios /= 0 ) call print_error( "input_event does't open!!" )

      EOF = .false.
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
!         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)
         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call obj% event% builder( 2*nevent )
            exit
         endif
      enddo
!
      do i = 1,obj% event% nevent/2
         read (u0,*) id, obj% event% init_state(id), obj% event% final_state(id),   &
                         obj% event% ebarrier(id), obj% event% de(id)
         write (*,*) id, obj% event% init_state(id), obj% event% final_state(id),   &
                         obj% event% ebarrier(id), obj% event% de(id)

!  ::: Event inverse...
         obj% event% init_state( id + nevent ) = obj% event% final_state( id )         
         obj% event% final_state( id + nevent ) = obj% event% init_state( id )
         obj% event% ebarrier( id + nevent ) =  obj% event% ebarrier( id ) - obj% event% de( id )         
         obj% event% de( id + nevent ) = - obj% event% de( id )
      enddo

      close( u0 )

    end subroutine read_event
! .................................................................................................
    
    subroutine distribution_state ( struc, per100 )
      use derived_types
      use random
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      real, intent( in )                :: per100
     
      integer                           :: i  
      real                              :: rnd

      call set_random_seed()

      do i = 1,struc% tot_sites
         call random_number( rnd )
         if ( rnd <= 1 - per100 ) then
            struc% site( i ) = 1
         else
            struc% site( i ) = 0
         endif
      enddo

    end subroutine distribution_state
! .................................................................................................

    subroutine event_rate_calc( obj )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: obj

      integer :: i, id, ievt, jevt, jn
      real    :: kt, f0

     ! integer, dimension( obj% tot_sites ) :: nevt
     ! integer, dimension( 10, obj% tot_sites ) :: evt_site

      kt = obj% kt   
      f0 = obj% f0

      do i = 1,obj% tot_sites
         do jn = 1,obj% nneig( i )
            obj% event_rate( jn, i ) = 0.0
         enddo
      enddo
!
! ::: We have 1 event by site : Adsorption or Desorption
      do i = 1,obj% tot_sites
         obj% rate( i ) = 0.0
         id = obj% site( i )

         ievt = 0
         do jevt = 1,obj% event% nevent
            if ( obj% event% init_state( jevt ) == id ) then
               ievt = ievt + 1
               obj% event_site( ievt, i ) = jevt
            endif
         enddo
         obj% nevt( i ) = ievt

         do jevt = 1, obj% nevt( i )
            obj% event_rate( jevt, i ) = obj% event% ebarrier( obj% event_site(jevt, i) ) !f0*exp( - obj% event% ebarrier(evt_site(jevt, i))/kt )
            obj% rate( i ) = obj% rate( i ) + obj% event_rate( jevt, i )
         enddo

         !write (*,*) " *** ", i, f0, obj% rate( i )

      enddo
!
    end subroutine event_rate_calc
! .................................................................................................

    subroutine choose_event( struc, isite, ievent )
      use derived_types
      use random
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      integer, intent( inout ) :: isite, ievent

      integer :: i, jn
      real :: rsum, rdn, rrdn

!      print*, " - To do :: choose_event "
      isite = 0 ; ievent = 0
      call random_number( rdn )
      rrdn = rdn*struc% sum_rate
      struc% rand_rate = rrdn
!      
      rsum = 0.0
      do i = 1,struc% tot_sites
      !
         do jn = 1,struc% nevt( i )
      !  !
            ievent = jn
            rsum = rsum + struc% event_rate( jn, i )
!            write (*,*) "rsum:",i ,jn ,rsum, struc% event_rate( jn, i )
            if ( rsum > rrdn ) exit
      !  !
         enddo
      !
         isite = i
         if ( rsum > rrdn ) exit
      enddo
      if ( rsum <= rrdn ) write (*,*) " PB choose event...", rsum,struc% sum_rate
!      write (*,*) struc% sum_rate,rdn,rrdn,rsum
!      write (*,'(2(a,1x,i6))') " We Choose on site ", isite," event ",ievent
!      write (*,*) struc% site( isite ), struc% site( struc% neig( ievent,isite ) )

    end subroutine choose_event
! .................................................................................................

    subroutine event_applied( struc, is, jn )
      use derived_types
      use errors
      implicit none
      type( KMC_type )      :: struc
      integer, intent( in ) :: is,   &  ! Site selected
                               jn       ! "ievent" of site "is" selected 
      integer :: jevt

      jevt = struc% event_site( jn, is )

      if ( struc% site( is ) /= struc% event% init_state( jevt ) ) then
         call warning( " Problem site state not correspond to initial state event ")
!         write (*,*) is,struc% site( is ), (struc% neig(j,is),struc% site( struc%neig(j,is) ),j=1,struc% nneig(is))
!         write (*,*) is,struc% rate( is ), (struc% neig(j,is),struc% event_rate( j,is ),j=1,struc% nneig(is))
      endif
      
      struc% site( is ) = struc% event% final_state( jevt )

    end subroutine event_applied
! .................................................................................................

    subroutine analyse( obj ) 
      use derived_types
      implicit none

      type( KMC_type ), intent( inout ) :: obj
      integer                        :: i

      !obj% txtprop( 1 ) = "  Cover (%) "
      obj% prop( 1 ) = 0.0

      do i = 1, obj% tot_sites
         if ( obj% site( i ) == 1 ) obj% prop( 1 ) = obj% prop(1) + 1
      enddo

      obj% prop(1) = obj% prop(1) / obj% tot_sites

    end subroutine analyse
! .................................................................................................
















