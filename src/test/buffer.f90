
module buffer
  use iso_c_binding
  implicit none
  integer( c_int ), dimension(:), allocatable, target :: a

contains

  subroutine fix_a( ptr_c, n ) bind( C )
    implicit none
    type( c_ptr ), intent( inout ) :: ptr_c
    integer( c_int ), intent( in ) :: n
    allocate( a(n) )
    ptr_c = c_loc( a(1) )
  end subroutine fix_a

  subroutine link_ptr( ptr_c, ptr_f, n ) bind( C )
    implicit none
    type( c_ptr ) :: ptr_c
    integer( c_int ), intent( in ) :: n
    integer( c_int ), dimension(:), pointer :: ptr_f
    call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
  end subroutine link_ptr

end module buffer
