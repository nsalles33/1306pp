
!  module diff_vac
!    use derived_types
!    implicit none

!   contains
!
!  USER routine:
!  - subroutine read_event( obj ) 
!  - subroutine distribution_state( obj, frac )
!  - event_rate_calc( obj )
!
! .................................................................................................

    subroutine read_event( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      character (len=50)                :: string
      integer                           :: u0, i, id, ios, nevent
      logical                           :: EOF

      character (len=1)  :: delims
      CHARACTER (len=100), dimension (50) :: args
      integer :: nargs

      !  ::: Lecture of state and rate of each node
      open( newunit=u0, file=trim(struc% input_event), iostat=ios )
      if ( ios /= 0 ) call print_error( "input_event does't open!!" )

      EOF = .false.
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)
         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call struc% event% builder( nevent )
            exit
         endif
      enddo
!
      do i = 1,struc% event% nevent
         read (u0,*) id, struc% event% init_state(id), struc% event% final_state(id),   &
                         struc% event% ebarrier(id), struc% event% de(id)
         write (*,*) id, struc% event% init_state(id), struc% event% final_state(id),   &
                         struc% event% ebarrier(id), struc% event% de(id)
      enddo

      close( u0 )
!
    end subroutine read_event
! .................................................................................................

    subroutine distribution_state ( struc, per100 )
      use derived_types
      use random
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      integer                           :: i
      real, intent( in )                :: per100
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

    subroutine event_rate_calc( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer :: i, jn, j, k, kn, dbd, l
      real :: kt, ebd, f0
      logical :: see

      integer, dimension( struc% tot_sites ) :: bd, ddb, v

      see = .true.
      kt = struc% kt
      ebd = struc% event% de( 1 )
      f0 = struc% f0
!
      do i = 1,struc% tot_sites
         struc% rate( i ) = 0.0
         do jn = 1,struc% nneig( i )
            struc% event_rate( jn, i ) = 0.0
         enddo
      enddo
!
!  -- We define the event rate with the environment 
!  1) How many dandling bond per vacancy
      do i = 1,struc% tot_sites
         struc% rate( i ) = 0.0

         if ( struc% site( i ) == 1 ) cycle

!    -- Dandling bond...
         ddb( i ) = 0 ; v( i ) = 0 ; bd( i ) = 0
         do jn = 1,struc% nneig( i )
            j = struc% neig( jn, i )
            if (j == 0) write (*,*) "PB system_rate: ",i,jn,j
            struc% event_rate( jn, i ) = 0.0
            if ( struc% site( j ) == 1 ) ddb( i ) = ddb( i ) + 1
            if ( struc% site( j ) == 0 ) v( i ) = v( i ) + 1
         enddo
         if ( ddb( i ) == 0 ) cycle
!
!    -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
!       ---------------------------------------------------------
!       dbd( 0 ) = v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2
!       ---------------------------------------------------------
!         write (*,*) i, struc% nneig( i ), bd( i ), ddb( i ), v( i )
         do jn = 1,struc% nneig( i )  ! *** neighbour is 1
            j = struc% neig( jn, i )

            if ( struc% site( j ) == 0 ) cycle

            ddb( j ) = 0 ; bd( j ) = 0 ; v( i ) = 0
            do kn = 1,struc% nneig(j)
               k = struc% neig( kn, j )
               if (k <= 0.or.k>struc%tot_sites)  &
                 write (*,*) k,"PB system_rate",i,jn,j,(l,struc% neig(l,j),l=1,struc% nneig(j))
               if ( struc% site( k ) == 0 ) ddb( j ) = ddb( j ) + 1
               if ( struc% site( k ) == 1 ) bd( j ) = bd( j ) + 1
            enddo

            dbd = v( i ) + bd( j ) - ddb( i ) - ddb( j ) + 2
            struc% event_rate( jn, i ) = f0*exp( - dbd*ebd/kt )
            struc% rate( i ) = struc% rate( i ) + struc% event_rate( jn, i )

!            write (*,*) j, struc% nneig( j ), bd( j ), ddb( j ), v( j )
!            write (*,*) " dbd :",v( i ),'+',bd( j ),'-',ddb( i ),'-',ddb( j ),'+ 2 =',dbd
!            write (*,*) i, jn, j, dbd,f0,kt,ebd, -dbd*ebd/kt, struc% event_rate( jn, i )

         enddo
!         write (*,*) " *** ", i, f0, ebd, struc% rate( i )
!

      enddo

    end subroutine event_rate_calc
! ..................................................................................................

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
         do jn = 1,struc% nneig( i )
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
! ..................................................................................................
    
    subroutine event_applied( struc, is, jn )
      implicit none
      type( KMC_type ) :: struc
      integer, intent( in ) :: is, jn
      integer :: j
!      print*, " - To do :: Event_applied "

!  ::: Flip 0 -> 1
      if ( struc% site( is ) == 1 ) then
         call warning( " site Flip Problem 0 = 1 ")
         write (*,*) is,struc% site( is ), (struc% neig(j,is),struc% site( struc%neig(j,is) ),j=1,struc% nneig(is))
         write (*,*) is,struc% rate( is ), (struc% neig(j,is),struc% event_rate( j,is ),j=1,struc% nneig(is))
      endif
      struc% site( is ) = 1
      j = struc% neig( jn, is )
      if ( struc% site( j ) == 0 ) call warning( " site Flip Problem 1 = 0 ")
      struc% site( j ) = 0

    end subroutine event_applied    
! ..................................................................................................


!  end module diff_vac
