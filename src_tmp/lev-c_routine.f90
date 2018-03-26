  module c_routine
   use iso_c_binding

   contains
!     subroutine error( text ) bind( C )
!       use iso_c_binding
!       implicit none
!       character ( kind=c_char ), dimension(*), intent( in ) :: text
!       integer( c_int ) :: len, i
!       len = 0
!       do
!          if ( text(len+1) == C_NULL_CHAR ) exit
!          len = len + 1
!       enddo
!       print*, " ERROR: ",(text(i),i=1,len)
!       stop " ===== FATAL ERROR!!"
!     end subroutine error

    subroutine builder_event_type( this, n ) bind( C )
      use iso_c_binding
      use derived_types
      implicit none
      type( event_type )   :: this
      integer( c_int ), intent( in ) :: n

      print*, " Enter in EVENT_type constructor..."
      this% nevent = n

!      allocate( h_i_state(n), h_f_state(n), h_Ebarrier(n), h_dE(n) )
!      if ( .not.allocated(h_i_state) .or. .not.allocated(h_f_state) .or.   &
!           .not.allocated(h_ebarrier) .or. .not.allocated(h_de) )                 &
!         call error( " EVENT_type => CONSTRUCTOR problem " )

!    this% ptr_i_state = c_loc( h_i_state(1) )
!    this% ptr_f_state = c_loc( h_f_state(1) )
!    this% ptr_ebarrier = c_loc( h_ebarrier(1) )
!    this% ptr_de = c_loc( h_de(1) )

!      print*, " EVENT_type constructor DONE "

    end subroutine builder_event_type
  end module c_routine
! ............................................................................
