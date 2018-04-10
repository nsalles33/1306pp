
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
    real( c_double ), dimension(:), pointer :: ebarrier, de, f0
    real( c_double ), dimension(:,:), pointer :: ebond
    !
    ! ::: Local update
    integer( c_int ) :: save_site, save_gp(25)
    !
    ! ::: Brownian Mvt
    integer( c_int ), dimension(3) :: rold, pos
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
      !
      type( KMC_type ), intent( inout ) :: struc
      character (len=500,kind=c_char)    :: string
      integer( c_int )                  :: u0, i, id, ibd, jbd, ios, nevent, nbond
      logical                           :: EOF
      !
      character (len=1,kind=c_char)  :: delims = " "
      CHARACTER (len=100,kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs
      !
      !  ::: INIT save site
      !
      save_site = 0 
      rold = 0
      !
      !  ::: Lecture of state and rate of each node
      !
      print*, " INPUT_EVENT: ",trim(struc% input_event)
      open( newunit=u0, file=trim(struc% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )
      !
      EOF = .false.
      do while ( .not.EOF )
         !
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         !
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)
         !
         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call builder_event_type( struc% event, nevent )
            !
            call link_int1_ptr( struc% event% ptr_i_state, init_state, struc% event% nevent)
            call link_int1_ptr( struc% event% ptr_f_state, final_state, struc% event% nevent)
            call link_real1_ptr( struc% event% ptr_f0, f0, struc% event% nevent)
            call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier, struc% event% nevent)
            call link_real1_ptr( struc% event% ptr_de, de, struc% event% nevent)
            !
            do i = 1,struc% event% nevent
               read (u0,*) id, init_state(id), final_state(id), f0(id), ebarrier(id), de(id)
               write (*,*) id, init_state(id), final_state(id), f0(id), ebarrier(id), de(id)
            enddo
            !
         endif
         !
         !
         if ( args(1) == "Energy_Bond" ) then
            !
            read( args(2), '(i5)' ) nbond
            !
            if ( nbond /= 0 ) then
               !
               call link_real2_ptr( struc% event% ptr_ebond, ebond, struc% event% nbond, struc% event% nbond)
               !
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
         endif
         !
      enddo
      !
      !
      close( u0 )
      !
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
      integer( c_int ) :: i, is, j0, nj, jn, j, k, k0, kn, dbd, l, evt0, ifirst, ilast
      real( c_double ) :: kt, ebd, eb

      integer, dimension( struc% tot_sites ) :: bd, ddb, v, ilist

      call link_int1_ptr( struc% ptr_site,             site,       struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig,            nneig,      struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig,             neig,       nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate,            rate,       struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate,      event_rate, nvois*struc% tot_sites )
      !
      call link_real1_ptr( struc% event% ptr_f0,       f0,         struc% event% nevent )
      call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier,   struc% event% nevent )
      call link_real1_ptr( struc% event% ptr_de,       de,         struc% event% nevent )
      call link_real2_ptr( struc% event% ptr_ebond,    ebond,      struc% event% nbond, struc% event% nbond )
      !
      !kt = struc% kt
      kt = -ebond(1,1)/4.0
      ebd = ebond(1,1)
      eb = ebarrier(1)
      !f0 = struc% f0
      !
      if ( save_site /= 0 ) then
         ifirst = 1
         ilast = save_site
         do i = ifirst,ilast
            ilist( i ) = save_gp(i)
         enddo
      else
         ifirst = 1
         ilast = struc% tot_sites
         do i = ifirst,ilast
            ilist( i ) = i
         enddo
      endif
      !
      !rate = 0.0
      !event_rate = 0.0
      !
      ! -- We define the event rate with the environment 
      ! 1) How many dandling bond per vacancy
      !do i = 1,struc% tot_sites
      do i = ifirst,ilast
         !
         is = ilist( i )
         !if ( save_site /= 0 ) is = save_site
         !
         !
         rate( is ) = 0.0
         !evt0 = ( is - 1 )*nvois + 1
         evt0 = nneig( is )
         do jn = 1,neig( evt0 )
            event_rate( evt0 + jn ) = 0.0
         enddo
         !
         if ( site( is ) == 1 ) cycle
         !
         ! -- Dandling bond...
         !
         ddb( is ) = 0 ; v( is ) = 0 ; bd( is ) = 0
         j0 = nneig( is )
         evt0 = j0  !! *** Conection between neighbour and event
         nj = neig( j0 )
         !
         do jn = 1,nj
            j = neig( j0 + jn )
            if (j == 0) write (*,*) "PB system_rate: ", is, jn, j
            !event_rate( evt0 + jn ) = 0.0
            if ( site( j ) > 0 ) ddb( is ) = ddb( is ) + 1
            if ( site( j ) == 0 ) v( is ) = v( is ) + 1
         enddo
         !
         if ( ddb( is ) == 0 ) cycle
         !
         ! -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
         !  ---------------------------------------------------------
         !  dbd( 0 ) = v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2
         !  ---------------------------------------------------------
         !write (*,*) i, struc% nneig( i ), bd( i ), ddb( i ), v( i )
         !
         do jn = 1,nj !neig( j0 )  ! *** neighbour is 1
            j = neig( j0 + jn )
            !
            if ( site( j ) == 0 ) cycle
            !
            ddb( j ) = 0 ; bd( j ) = 0 !; v( i ) = 0
            k0 = nneig( j )
            do kn = 1,neig( k0 )
               k = neig( k0 + kn )
               if (k <= 0.or.k > struc% tot_sites)  &
                 write (*,*) k,"PB system_rate",i,jn,j,(l,neig(k0+l),l=1,neig(k0))
               if ( site( k ) == 0 ) ddb( j ) = ddb( j ) + 1
               if ( site( k ) == 1 ) bd( j ) = bd( j ) + 1
            enddo
            !
            dbd = v( is ) + bd( j ) - ddb( is ) - ddb( j ) + 2
            if ( - dbd <= 0 ) then
               event_rate( j0 + jn ) = f0(1)*exp( - eb/kt )
            elseif ( - dbd > 0 ) then
               event_rate( j0 + jn ) = f0(1)*exp( - (eb + (-dbd)*ebd)/kt )
            endif
            rate( is ) = rate( is ) + event_rate( j0 + jn )
            !
            !write (*,*) j, nneig( j ), bd( j ), ddb( j ), v( j )
            !write (*,*) " dbd :",v( i ),'+',bd( j ),'-',ddb( i ),'-',ddb( j ),'+ 2 =',dbd
            !write (*,*) is, j0, jn, j, site( j ), f0(1), kt, (-dbd)*ebd/kt, event_rate( j0 + jn )
            !
         enddo
         !
         !write (*,'(a,1x,i8,3(1x,e10.2))') " *** ", i, f0(1), ebd, rate( i )
         !write (*,*) " *** ", is, f0(1), ebd, rate( is )
         !
         !if ( save_site /= 0 ) exit
         !
      enddo
      !
      !stop "Event_rate_calc..."
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

      call link_int1_ptr( struc% ptr_site,        site,       struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig,       nneig,      struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig,        neig,       nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate,       rate,       struc% tot_sites )
      !
      !
      isite = 0 ; ievent = 0
      call random_number( rdn )
      rrdn = rdn*struc% sum_rate
      struc% rand_rate = rrdn
      !      
      rsum = 0.0
      do i = 1,struc% tot_sites
         !
         !if ( site( i ) == 1 ) cycle
         j0 = nneig( i )
         !write(*,*) i,j0,(jn,neig( j0+jn ),jn=0,4)
         do jn = 1,neig( j0 )
         !  !
            ievent = jn
            rsum = rsum + event_rate( j0 + jn )
            !write (*,*) "rsum:",i, j0, jn ,rsum, event_rate( j0 + jn ), rate( i )
            if ( rsum > rrdn ) exit
         !  !
         enddo
         !
         isite = i
         if ( rsum > rrdn ) exit
      enddo
      !
      if ( rsum <= rrdn ) then
         write (*,*) " PB choose event...", rsum, rrdn, struc% sum_rate
         stop " EXIT..."
      endif
      !
      !write (*,*) struc% sum_rate,rdn,rrdn,rsum
      !write (*,'(2(a,1x,i6))') " We Choose on site ", isite," event ",ievent
      !write (*,*) site( isite ), site( neig( j0 + ievent ) )
      !
    end subroutine choose_event
! ..................................................................................................
    
    subroutine event_applied( struc, is, jn ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use variables
      implicit none
      !
      type( KMC_type )               :: struc
      integer( c_int ), intent( in ) :: is, jn
      integer( c_int )               :: j, jv, j0, k, k0, kv,l, o0, ov, &
                                        igp, see, nx, ny, nz, nxy
      !
      integer( c_int ), dimension(3) :: x
      !
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
      !
      site( is ) = 1
      !
      ! ::: Empty the table
      rate( is ) = 0.0
      do jv = 1,neig( j0 )
         event_rate( j0 + jv ) = 0.0   
      enddo
      !
      j = neig( j0 + jn )
      if ( site( j ) == 0 ) call warning( " site Flip Problem 1 = 0 ")
      !
      site( j ) = 0
      !
      !
      ! ::: Deplacement :::
      !
      nx = struc% nsites(1)
      ny = struc% nsites(2)
      nz = struc% nsites(3)
      nxy = nx*ny
      !
      x( 1 ) = MODULO( j-1, nx )
      x( 2 ) = MODULO( (j-1) / nx, ny )
      x( 3 ) = 1
      if ( struc% sys_dim == 3 ) x( 3 ) = (j-1) / nxy
      !
      !print*,"pos",pos
      pos = pos + (x - rold)**2 
      !print*, " Deplacment: "
      !print*,rold
      !print*,x   
      !print*,"pos",pos
      !
      !
      ! ::: Save local neig
      j0 = nneig( j )
      save_gp(1) = j
      igp = 1
      do jv = 1,neig( j0 )
         !
         ! ::: Verif list
         see = 0
         do l = 1, igp
            if ( save_gp( l ) == neig( j0 + jv ) ) see = 1
         enddo
         if ( see == 1 ) cycle
         !
         ! ::: Include site
         igp = igp + 1
         k = neig( j0 + jv )
         save_gp( igp ) = k
         k0 = nneig( k )
         do kv = 1,neig( k0 )
            !
            ! ::: Verif list
            see = 0
            do l = 1, igp
               if ( save_gp( l ) == neig( k0 + kv ) ) then
                  see = 1
                  exit
               endif
            enddo
            if ( see == 1 ) cycle
            !
            ! ::: Include site
            igp = igp + 1
            save_gp( igp ) = neig( k0 + kv ) 
            o0 = neig( k0 + kv ) 
            o0 = nneig( o0 )
            !
!           do ov = 1,neig( o0 )
!              !
!              ! ::: Verif list
!              see = 0
!              do l = 1, igp
!                 if ( save_gp( l ) == neig( o0 + ov ) ) then
!                    see = 1
!                    exit
!                 endif
!              enddo
!              if ( see == 1 ) cycle
!              !
!              igp = igp + 1
!              save_gp( igp ) = neig( o0 + ov )
!              !
!           enddo
            !
         enddo
      enddo
      !
      save_site = igp
      !do j = 1,save_site
      !   write (*,*)  j, save_gp( j )
      !enddo
      !
      !save_site = 0
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
      integer( c_int ) :: i, j, j0, jv, k, k0, kv,   &
                          ngp, gpv, nc, mixgp, nvac, &
                          maxsize, max_gp, nx, ny, nz, nxy
      real( c_double ) :: x, y, z
      !
      integer( c_int ), dimension( obj% tot_sites ) :: gp, histo
      integer( c_int ), dimension( -1:int(obj% tot_sites/2) ) :: clster
      !
      call link_int1_ptr( obj% ptr_site, site, obj% tot_sites )
      call link_int1_ptr( obj% ptr_nneig, nneig, obj% tot_sites )
      call link_int1_ptr( obj% ptr_neig, neig, nvois*obj% tot_sites )
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )
      !
      !
      ! ::: Cluster size :::
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
         !write (*,*) prop(1),i,histo(i), clster(-1) , max_gp
      enddo
      !
      !
      if ( obj% nprop == 1 ) return
      !
      !
      !
      ! ::: Variance compute :::
      !
      do i = 1,obj% tot_sites
         if ( site( i ) == 0 ) then
            j = i
            exit
         endif
      enddo
      !
      if ( sum(rold) == 0 ) then
         print*, " First step: init R_old..."
         nx = obj% nsites(1)
         ny = obj% nsites(2)
         nz = obj% nsites(3)
         nxy = nx*ny
         !
         rold( 1 ) = MODULO( j-1, nx )
         rold( 2 ) = MODULO( (j-1) / nx, ny )
         rold( 3 ) = 1
         if ( obj% sys_dim == 3 ) rold( 3 ) = (j-1) / nxy
         !
         pos = 0
         print*,rold,pos
         !
      else
         !
         x = pos(1); y = pos(2); z = pos(3)
         !print*, (x**2 + y**2 + z**2), obj% step
         prop( 2 ) = sqrt( (x**2 + y**2 + z**2)/real( obj% step ) )
         !
      endif
      !
      !
      
      !
    end subroutine analyse
! ..................................................................................................
! ..................................................................................................

!  end module diff_vac







