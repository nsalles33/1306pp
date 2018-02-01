!
!    Routine to calculate the rate of each event and each site  
!
!
    subroutine rate_sum_calc( struc )
      use derived_types
      !use diff_vac
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      integer :: i
      real :: sum_rate


!      print*, " Dans rate_sum_calc..."

      call event_rate_calc( struc )
!
!  -- Sum over all the site
!
      sum_rate = 0
      do i = 1,struc% tot_sites
         sum_rate = sum_rate + struc% rate( i )
      enddo
      struc% sum_rate = sum_rate
!
!      write (*,*) " DONE "
!
    end subroutine rate_sum_calc
! ..................................................................................................

!   subroutine event_rate_calc( struc )
!     use derived_types
!     use errors
!     implicit none

!     type( KMC_type ), intent( inout ) :: struc
!     integer :: i, jn, j, k, kn, dbd, l
!     real :: kt, ebd, f0, sum_rate
!     logical :: see

!     integer, dimension( struc% tot_sites ) :: bd, ddb, v

!     see = .true.
!     kt = struc% kt
!     ebd = struc% event% de( 1 )
!     f0 = struc% f0
!
!     do i = 1,struc% tot_sites
!        struc% rate( i ) = 0.0
!        do jn = 1,struc% nneig( i )
!           struc% event_rate( jn, i ) = 0.0
!        enddo
!     enddo
!
!  -- We define the event rate with the environment 
!  1) How many dandling bond per vacancy
!     do i = 1,struc% tot_sites
!        struc% rate( i ) = 0.0

!        if ( struc% site( i ) == 1 ) cycle

!    -- Dandling bond...
!        ddb( i ) = 0 ; v( i ) = 0 ; bd( i ) = 0
!        do jn = 1,struc% nneig( i )
!           j = struc% neig( jn, i )
!           if (j == 0) write (*,*) "PB system_rate: ",i,jn,j
!           struc% event_rate( jn, i ) = 0.0
!           if ( struc% site( j ) == 1 ) ddb( i ) = ddb( i ) + 1
!           if ( struc% site( j ) == 0 ) v( i ) = v( i ) + 1
!        enddo
!        if ( ddb( i ) == 0 ) cycle
!
!    -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
!       ---------------------------------------------------------
!       dbd( 0 ) = v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2
!       ---------------------------------------------------------
!         write (*,*) i, struc% nneig( i ), bd( i ), ddb( i ), v( i )
!        do jn = 1,struc% nneig( i )  ! *** neighbour is 1
!           j = struc% neig( jn, i )

!           if ( struc% site( j ) == 0 ) cycle

!           ddb( j ) = 0 ; bd( j ) = 0 ; v( i ) = 0
!           do kn = 1,struc% nneig(j)
!              k = struc% neig( kn, j )
!              if (k <= 0.or.k>struc%tot_sites)  &
!                write (*,*) k,"PB system_rate",i,jn,j,(l,struc% neig(l,j),l=1,struc% nneig(j))
!              if ( struc% site( k ) == 0 ) ddb( j ) = ddb( j ) + 1
!              if ( struc% site( k ) == 1 ) bd( j ) = bd( j ) + 1
!           enddo

!           dbd = v( i ) + bd( j ) - ddb( i ) - ddb( j ) + 2
!           struc% event_rate( jn, i ) = f0*exp( - dbd*ebd/kt )
!           struc% rate( i ) = struc% rate( i ) + struc% event_rate( jn, i )

!            write (*,*) j, struc% nneig( j ), bd( j ), ddb( j ), v( j )
!            write (*,*) " dbd :",v( i ),'+',bd( j ),'-',ddb( i ),'-',ddb( j ),'+ 2 =',dbd
!            write (*,*) i, jn, j, dbd,f0,kt,ebd, -dbd*ebd/kt, struc% event_rate( jn, i )

!        enddo
!         write (*,*) " *** ", i, f0, ebd, struc% rate( i )
!

!     enddo

!   end subroutine event_rate_calc
! ..................................................................................................





