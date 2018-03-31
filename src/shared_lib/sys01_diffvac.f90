
  module variables
    use iso_c_binding
    implicit none
    !
    integer( c_int ), dimension(:), pointer :: site, nneig, neig, nevt, spec, event_site
    !
    real( c_double ), dimension(:), pointer :: rate, prop, event_rate
    !
    ! ::: EVENT
    integer( c_int ), dimension(:), pointer :: init_state, final_state
    real( c_double ), dimension(:), pointer :: ebarrier, de
    real( c_double ), dimension(:,:), pointer :: ebond
    !
  end module variables

!
!  USER routine:
!  - subroutine read_event( obj ) 
!  - subroutine event_rate_calc( obj )
!  - subroutine 
!
! .................................................................................................

    subroutine read_event( struc ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use errors
      use variables
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      character (len=500,kind=c_char)    :: string
      integer( c_int )                  :: u0, i, id, ibd, jbd, ios, nevent, nbond
      logical                           :: EOF

      character (len=1,kind=c_char)  :: delims = " "
      CHARACTER (len=100,kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs

      !integer( c_int ), dimension(:), pointer :: init_state, final_state
      !real( c_double ), dimension(:), pointer :: ebarrier, de
      !real( c_double ), dimension(:,:), pointer :: ebond

!  ::: Lecture of state and rate of each node
      print*, " INPUT_EVENT: ",trim(struc% input_event)
      open( newunit=u0, file=trim(struc% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )
      !
      EOF = .false.
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)

         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call builder_event_type( struc% event, nevent )
            exit
         endif

      enddo
      !
      call link_int1_ptr( struc% event% ptr_i_state, init_state, struc% event% nevent)
      call link_int1_ptr( struc% event% ptr_f_state, final_state, struc% event% nevent)
      call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier, struc% event% nevent)
      call link_real1_ptr( struc% event% ptr_de, de, struc% event% nevent)
      !
      do i = 1,struc% event% nevent
         read (u0,*) id, init_state(id), final_state(id), ebarrier(id), de(id)
         write (*,*) id, init_state(id), final_state(id), ebarrier(id), de(id)
      enddo
      !
      EOF = .false.
      nbond = 0
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)
         
         if ( args(1) == "Energy_Bond" ) then
            read( args(2), '(i5)' ) nbond
            exit
         endif
      enddo
      !
      if ( nbond /= 0 ) then
         call link_real2_ptr( struc% event% ptr_ebond, ebond, struc% event% nbond, struc% event% nbond)
         do i = 1,nbond
            call read_line( u0, string, EOF )
            call parse( trim(string), delims, args, nargs ) 
            read( args(2), '(i5)') ibd 
            read( args(3), '(i5)') jbd 
            read( args(4), '(f10.5)') ebond(ibd,jbd)
            write (*,*) i, ibd, jbd, ebond( ibd, jbd )
         enddo
      endif
      !
      close( u0 )
!     !
    end subroutine read_event
! .................................................................................................

    subroutine event_rate_calc( struc ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use errors
      use variables
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ) :: i, j0, nj, jn, j, k, k0, kn, dbd, l, evt0
      real( c_double ) :: kt, ebd, f0, eb

      integer, dimension( struc% tot_sites ) :: bd, ddb, v

      !integer( c_int ), dimension(:), pointer :: site, nneig, neig
      !integer( c_int ), dimension(:,:), pointer :: neig
      !real( c_double ), dimension(:), pointer :: de, rate, event_rate, ebarrier
      !real( c_double ), dimension(:,:), pointer :: ebond
      call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig, neig, nvois*struc% tot_sites )
      call link_real1_ptr( struc% event% ptr_de, de, struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier, struc% event% nevent )
      call link_real2_ptr( struc% event% ptr_ebond, ebond, struc% event% nbond, struc% event% nbond )
      !
      !kt = struc% kt
      kt = -ebond(1,1)
      ebd = ebond(1,1)
      eb = 0.0 !ebarrier(1)
      f0 = struc% f0
      !
      rate = 0.0
      event_rate = 0.0
      !
      ! -- We define the event rate with the environment 
      ! 1) How many dandling bond per vacancy
      do i = 1,struc% tot_sites
         rate( i ) = 0.0
         evt0 = (i-1)*nvois + 1
         do jn = 1,nvois
            event_rate( evt0 + jn ) = 0.0
         enddo
         !
         if ( site( i ) == 1 ) cycle
         !
         ! -- Dandling bond...
         ddb( i ) = 0 ; v( i ) = 0 ; bd( i ) = 0
         j0 = nneig( i )
         evt0 = j0  !! *** Conection between neighbour and event
         nj = neig( j0 )
         !
         do jn = 1,nj
            j = neig( j0 + jn )
            if (j == 0) write (*,*) "PB system_rate: ",i,jn,j
            !event_rate( evt0 + jn ) = 0.0
            if ( site( j ) == 1 ) ddb( i ) = ddb( i ) + 1
            if ( site( j ) == 0 ) v( i ) = v( i ) + 1
         enddo
         !
         if ( ddb( i ) == 0 ) cycle
         !
         ! -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
         !  ---------------------------------------------------------
         !  dbd( 0 ) = v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2
         !  ---------------------------------------------------------
         !write (*,*) i, struc% nneig( i ), bd( i ), ddb( i ), v( i )
         do jn = 1,neig( j0 )  ! *** neighbour is 1
            j = neig( j0 + jn )
            !
            if ( site( j ) == 0 ) cycle
            !
            ddb( j ) = 0 ; bd( j ) = 0 ; v( i ) = 0
            k0 = nneig( j )
            do kn = 1,neig( k0 )
               k = neig( k0 + kn )
               if (k <= 0.or.k > struc% tot_sites)  &
                 write (*,*) k,"PB system_rate",i,jn,j,(l,neig(k0+l),l=1,neig(k0))
               if ( site( k ) == 0 ) ddb( j ) = ddb( j ) + 1
               if ( site( k ) == 1 ) bd( j ) = bd( j ) + 1
            enddo
            !
            dbd = v( i ) + bd( j ) - ddb( i ) - ddb( j ) + 2
            if ( - dbd <= 0 ) then
               event_rate( j0 + jn ) = f0*exp( - eb/kt )
            elseif ( - dbd > 0 ) then
               event_rate( j0 + jn ) = f0*exp( - (eb + (-dbd)*ebd)/kt )
            endif
            rate( i ) = rate( i ) + event_rate( j0 + jn )
            !
            !write (*,*) j, nneig( j ), bd( j ), ddb( j ), v( j )
            !write (*,*) " dbd :",v( i ),'+',bd( j ),'-',ddb( i ),'-',ddb( j ),'+ 2 =',dbd
            !write (*,*) i, jn, j, f0,kt,dbd, dbd*ebd/kt, event_rate( j0 + jn )
            !
         enddo
         !write (*,*) " *** ", i, f0, ebd, rate( i )
         !
      enddo
      !
    end subroutine event_rate_calc
! ..................................................................................................

    subroutine choose_event( struc, isite, ievent ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use random
      use variables
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ), intent( inout ) :: isite, ievent

      integer( c_int ) :: i, j0, jn
      real( c_double ) :: rsum, rdn, rrdn

      !integer( c_int ), dimension(:), pointer :: nneig, neig
      !real( c_double ), dimension(:), pointer :: event_rate
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig, neig, nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      !
      !print*, " - To do :: choose_event "
      isite = 0 ; ievent = 0
      call random_number( rdn )
      rrdn = rdn*struc% sum_rate
      struc% rand_rate = rrdn
      !      
      rsum = 0.0
      do i = 1,struc% tot_sites
         !
         j0 = nneig( i )
         !write(*,*) i,j0,(jn,neig( j0+jn ),jn=0,4)
         do jn = 1,neig( j0 )
         !  !
            ievent = jn
            rsum = rsum + event_rate( j0 + jn )
            !write (*,*) "rsum:",i,j0, jn ,rsum, event_rate( j0+jn )
            if ( rsum > rrdn ) exit
         !  !
         enddo
         !
         isite = i
         if ( rsum > rrdn ) exit
      enddo
      if ( rsum <= rrdn ) write (*,*) " PB choose event...", rsum,struc% sum_rate
      !write (*,*) struc% sum_rate,rdn,rrdn,rsum
      !write (*,'(2(a,1x,i6))') " We Choose on site ", isite," event ",ievent
      !write (*,*) struc% site( isite ), struc% site( struc% neig( ievent,isite ) )
      !
    end subroutine choose_event
! ..................................................................................................
    
    subroutine event_applied( struc, is, jn ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use variables
      implicit none
      type( KMC_type ) :: struc
      integer( c_int ), intent( in ) :: is, jn
      integer( c_int ) :: j, j0
      !
      !integer( c_int ), dimension(:), pointer :: site, nneig, neig
      !real( c_double ), dimension(:), pointer :: rate, event_rate 
      !
      call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig, neig, nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      !
      !  ::: Flip 0 -> 1
      !
      j0 = nneig( is )
      if ( site( is ) == 1 ) then
         call warning( " site Flip Problem 0 = 1 ")
         write (*,*) is, site( is ), ( neig(j0+j), site( neig(j0+j) ),j=1,neig(j0) )
         write (*,*) is, rate( is ), ( neig(j0+j), event_rate( j0+j ),j=1,neig(j0) )
      endif
      site( is ) = 1
      j = neig( j0 + jn )
      if ( site( j ) == 0 ) call warning( " site Flip Problem 1 = 0 ")
      site( j ) = 0
      !
    end subroutine event_applied    
! ..................................................................................................
!
    subroutine analyse( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use variables
      implicit none
      !
      type( KMC_type ) :: obj
      integer( c_int ) :: i, j, j0, jv, k, k0, kv, ngp, gpv, nc, mixgp, nvac, maxsize, max_gp
      !
      integer( c_int ), dimension( obj% tot_sites ) :: gp, histo
      integer( c_int ), dimension( -1:int(obj% tot_sites/2) ) :: clster
      !
      !integer( c_int ), dimension(:), pointer :: site, nneig, neig
      !real( c_double ), dimension(:), pointer :: prop
      call link_int1_ptr( obj% ptr_site, site, obj% tot_sites )
      call link_int1_ptr( obj% ptr_nneig, nneig, obj% tot_sites )
      call link_int1_ptr( obj% ptr_neig, neig, nvois*obj% tot_sites )
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )
      !
      !
      ngp = 0 ; gp = 0 ; max_gp = 0
      clster = 0 ; nvac = 0
      do i = 1,obj% tot_sites
         if ( site( i ) == 1 ) cycle
         !
         nvac = nvac + 1
         nc = 0
         gp( i ) = -1
         gpv = 0 ; mixgp = 0
         j0 = nneig( i )
         do jv = 1,neig( j0 )
            j = neig( j0 + jv )
            !
            if ( site( j ) == 1 ) cycle
            !
            if ( gp( j ) /= 0.and.gpv == 0 ) then
               gpv = gp( j )
            elseif ( gp( j ) /= 0.and.gpv /= 0.and.gp(j) /= gpv) then
               mixgp = gp( j )
            endif
            !
            nc = nc + 1
            !
         enddo
         !write (*,*) i,nc,ngp,gpv,mixgp
         !
         if ( nc /= 0.and.gpv /= 0.and.mixgp == 0 ) then
            gp( i ) = gpv
         elseif ( nc /= 0.and.gpv == 0 ) then
            ngp = ngp + 1
            gp( i ) = ngp
            do jv = 1,neig( j0 )
               j = neig( j0 + jv )
               if ( site( j ) == 0 ) gp( j ) = gp( i )
            enddo
         elseif ( mixgp /= 0 ) then
            !
            gp( i ) = min( mixgp, gpv )
            do j = 1,i-1
               if ( site( j ) == 1 ) cycle
               if ( gp( j ) == mixgp.or.gp( j ) == gpv ) then
                  gp( j ) = gp( i )
                  k0 = nneig( j )
                  do kv = 1,neig( k0 )
                     k = neig( k0 + kv )
                     if ( site( k ) == 0 ) gp( k ) = gp( i )
                  enddo
               endif
            enddo
            !
         endif
         !
         clster( gp(i) ) = clster( gp(i) ) + 1
         max_gp = max( max_gp, gp(i) )
         !write (*,*) i, gp(i), ngp, clster( gp(i) ), nvac
         !
         !if (i < obj% nsites(1).or. i > obj%nsites(2)*(obj%nsites(1)-1)+0) &
         !  write (*,*) i," neig ", nneig(i),(neig(nneig(i)+j),j=0,4)
      enddo
      !
      nvac = 0
      histo = 0 ; maxsize = 1
      histo( 1 ) = clster( -1 )
      !
      do i = 1,max_gp
         nvac = nvac + clster( i )
         histo( clster(i) ) = histo( clster(i) ) + 1
         maxsize = max( maxsize, clster(i) )
      enddo
      !
      prop( 1 ) = 0
      do i = 1,maxsize
         prop( 1 ) = prop( 1 ) + real(i*histo( i ))/real( clster(-1) + max_gp )
!         write (*,*) prop(1),i,histo(i), clster(-1) , max_gp
      enddo
      !
    end subroutine analyse
! ..................................................................................................
! ..................................................................................................

!  end module diff_vac







