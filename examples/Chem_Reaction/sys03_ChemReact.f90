
  module user_system
    use iso_c_binding
    use derived_types
    use sub_new_types
    use errors
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
   contains
    !
    !  USER routine:
    !  - subroutine read_event( obj ) 
    !  - subroutine event_rate_calc( obj )
    !  - subroutine choose_event( obj, is, jn )
    !  - subroutine event_applied( obj, is, jn )
    !  - subroutine analyse( obj )
    !
! .................................................................................................
    !
    subroutine alloc_chem_event_type( this, nspec )
      implicit none
      type( event_type )   :: this
      integer( c_int ), intent ( in ) :: nspec
      integer( c_int ) :: n

      print*, " Enter in EVENT_type constructor..."
      !this% nevent = nevt
      n = this% nevent + this% nchem_react*nspec
      print*, " Allocation size :", n, this% nevent, this% nchem_react, nspec
      !
      allocate( h_i_state(n), h_f_state(n), h_ebarrier(n), h_dE(n))
      if ( .not.allocated(h_i_state) .or.  &
           .not.allocated(h_f_state) .or.  &
           .not.allocated(h_ebarrier) .or. &
           .not.allocated(h_de) )          &
         call error( " EVENT_type => CONSTRUCTOR problem " )
      !
      if ( this% nbond /= 0 ) then
         allocate( h_ebond(this% nbond, this% nbond) )
         if ( .not.allocated(h_ebond) )          &
            call error( " EVENT_type => CONSTRUCTOR bond problem " )
      endif
      !
      !if ( nspec /= 0 ) then
      !   allocate( h_spec( this% nspec ) )
      !   if ( .not.allocated(h_init_spec) )      &
      !      call error( " EVENT_type => CONSTRUCTOR init_spec problem " )
      !endif
      !
      this% ptr_i_state      = c_loc( h_i_state(1) )
      this% ptr_f_state      = c_loc( h_f_state(1) )
      this% ptr_ebarrier     = c_loc( h_ebarrier(1) )
      this% ptr_de           = c_loc( h_de(1) )
      if ( this% nbond /= 0 )  &
         this% ptr_ebond     = c_loc( h_ebond(1,1) )
      !if ( this% nspec /= 0 ) &
      !   this% ptr_init_spec = c_loc( h_init_spec(1) )
      !
      print*, " EVENT_type constructor DONE "
      !
    end subroutine alloc_chem_event_type
! .................................................................................................
!
    subroutine read_event( struc )
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      character (len=50,kind=c_char)    :: string
      integer( c_int )                  :: u0, i, j, id, ios, nspec
      logical                           :: EOF

      integer( c_int ), dimension(:), pointer :: spec

      character (len=1,kind=c_char)  :: delims
      CHARACTER (len=100,kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs
      delims = " "
      !
      call link_int1_ptr( struc% ptr_spec, spec,        struc% nspec)
      !
      !  ::: Initialize event parameter
      !
      struc% event% nevent      = 0
      struc% event% nbond       = 0
      struc% event% nchem_react = 0
      !
      !  ::: Lecture of state and rate of each node
      !
      print*, " INPUT_EVENT: ",trim(struc% input_event)
      open( newunit=u0, file=trim(struc% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )
      !
      EOF = .false.
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         !write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)
         !
         if ( args(1) == "Init_species" ) then
            read( args(2), '(i5)' ) nspec
             if ( nspec /= struc% nspec) call error( " PB Number of Init species..." )
            write (*,*) trim(args(1)), struc% nspec
            do i = 1,struc% nspec
               read (u0,*) args(1)
               read( args(1), '(i5)' ) spec( i )
               print*, "spec", i, "value", spec(i)
            enddo
         endif
         !
         if ( args(1) == "Number_of_reaction" ) then
            read( args(2), '(i5)' ) struc% event% nchem_react
            call alloc_chem_event_type( struc% event, struc% nspec )
            exit
         endif
         !
      enddo
      !
      id = struc% event% nchem_react*struc% nspec
      call link_int1_ptr( struc% event% ptr_i_state,   init_state,  id )
      call link_int1_ptr( struc% event% ptr_f_state,   final_state, id )
      call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier,    id )
      call link_real1_ptr( struc% event% ptr_de,       de,          id )
      !
      do i = 1,struc% event% nchem_react
         id = ( i - 1 )*struc% nspec
         read( u0, * ) ios, ( init_state( id + j ), j=1,struc% nspec ), ( final_state( id + j ), j=1,struc% nspec ), &
                           ebarrier( i ), de( i )
         write ( *, * ) id, ( init_state( id + j ), j=1,struc% nspec ), ( final_state( id + j ), j=1,struc% nspec ), &
                             ebarrier( i ), de( i )
      enddo 
      !
      close( u0 )
      !stop "read_event"
      !
    end subroutine read_event
! .................................................................................................

    subroutine event_rate_calc( obj )
      implicit none
      type( kmc_type ), intent (inout) :: obj
      !
      integer( c_int ) :: is, ir, j, n, evt( obj% event% nchem_react )
      real( c_double ) :: kt, h
      !
      !integer( c_int ), dimension(:), pointer :: site, nevt, init_state, event_site
      !real( c_double ), dimension(:), pointer :: ebarrier
      !
      n = obj% nspec*obj% event% nchem_react
      call link_int1_ptr( obj% ptr_site,             site,        obj% tot_sites )
      call link_int1_ptr( obj% ptr_nevt,             nevt,        obj% tot_sites )
      call link_int1_ptr( obj% event% ptr_i_state,   init_state,  n )
      call link_real1_ptr( obj% event% ptr_ebarrier, ebarrier,    n )
      call link_int1_ptr( obj% ptr_event_site,       event_site,  nvois*obj% tot_sites )
      call link_real1_ptr( obj% ptr_event_rate,      event_rate,  nvois*obj% tot_sites )
      !
      kt = obj% kt
      !f0 = obj% f0
      !
      !rate = 0.0
      event_rate = 0.0
      !
      ! 1) est ce possible de faire la reaction
      ! 2) calcul du temps de la reaction
      ! 3) Conservation du temps le plus bas
      !
      nevt = 0
      do is = 1,obj% tot_sites
         !
         if ( site( is ) == 0 ) cycle
         !
         evt = 0
         do ir = 1,obj% event% nchem_react
            !
            evt( ir ) = 1
            do j = 1,obj% nspec 
               if ( spec( j ) < init_state( ir - 1 + j ) ) evt( ir ) = 0
               !print*, j, site( j ), ir - 1 + j, init_state( ir - 1 + j ), nevt( ir )
            enddo
            !
            if ( evt( ir ) == 0 ) cycle
            nevt( is ) = nevt( is ) + evt( ir )
            !
            h = 0.0
            if ( ir == 1 ) h = real( spec(1) )*real( spec(2) )
            if ( ir == 2.and.spec(2) > 0) h = real( spec(2)*( spec(2) - 1 ) )/2.0
            event_rate( ir ) = h*ebarrier( ir )
            !print*, ir," chem react :", event_rate( ir ), evt( ir ), h, ebarrier( ir ), nevt( is )
            !
         enddo
         !
      enddo
      !
      !stop " event_rate_calc..."
      !
    end subroutine event_rate_calc
! .................................................................................................

    !real( c_double ) function Chem_combination( ievent ) result ( h )    
 !! 
   ! end function Chem_combination
! .................................................................................................

    subroutine choose_event( obj, is, jn )
      implicit none
      type( kmc_type ) :: obj
      integer( c_int ), intent( out ) :: is, jn
      !
      integer( c_int ) :: i, nber, isite
      real( c_double ) :: rdn, tmin, t
      !
      !integer( c_int ), dimension(:), pointer :: nevt
      !real( c_double ), dimension(:), pointer :: event_rate
      !
      tmin = 1000
      is = 0; jn = 0; nber = 0
      obj% rand_time = 0.0
      !
      call link_int1_ptr( obj% ptr_nevt,             nevt,        obj% tot_sites )
      call link_real1_ptr( obj% ptr_event_rate,      event_rate,  nvois*obj% tot_sites )
      !
      do isite = 1, obj% tot_sites
         !
         do i = 1,obj% event% nchem_react
            !
            if ( nevt( isite ) == 0 ) cycle
            nber = nber + nevt( isite )
            !
            call random_number( rdn )
            t = -log( rdn )/event_rate( i )
            !print*, i, event_rate( i ), rdn, t
            !
            if ( t < tmin ) then
               tmin = t
               is = i
            endif
            !
         enddo
         !
         if ( nber /= 0 ) then
            obj% rand_time = tmin
         else
            obj% conv = 1
         endif
         !
      enddo
      !print*, " event : ",is, event_rate( is ), obj% rand_time
      !stop " Choose_event..."
      !
    end subroutine choose_event
! .................................................................................................

    subroutine event_applied( obj, is , jn )
      implicit none
      type( kmc_type ) :: obj
      integer( c_int ), intent( in ) :: is, jn
      integer( c_int ) :: i, id, n
      !
      n = obj% nspec*obj% event% nchem_react
      call link_int1_ptr( obj% ptr_spec,             spec,        obj% nspec )
      call link_int1_ptr( obj% event% ptr_f_state,   final_state, n )
      !
      id =  (is - 1 )*obj%nspec
      do i = 1,obj% nspec
         spec( i ) = spec( i ) + final_state( id + i )
         !print*, is, id , final_state( id + i ), i, spec( i )
      enddo
      !
      !stop " event_applied..."
      !
    end subroutine event_applied
! .................................................................................................

    subroutine analyse( obj )
      implicit none
      type( kmc_type ) :: obj
      integer( c_int ) :: i
      !
      call link_int1_ptr( obj% ptr_spec, spec, obj% nspec )
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )
      !
      if ( obj% nprop /= obj% nspec ) call error( " ERROR change analyse routine... nprop/=nspec" )
      !
      do i = 1,obj% nprop
         prop( i ) = spec( i ) 
         !print*, " analyse", i, prop(i)
      enddo
      !
    end subroutine
! .................................................................................................

  end module user_system
! =================================================================================================










