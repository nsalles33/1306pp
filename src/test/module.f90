

module type_def
  use iso_c_binding
  implicit none

  type, public, bind( C ) :: sub_t
    integer( c_int ) :: nevt
    real( c_double ), dimension( 10 ) :: tab
  end type


  type, public, bind( C ) :: test_t
    integer( c_int ) :: n
    integer( c_int ), dimension(3) :: p
    character( kind=c_char, len=20 ) :: libname,  &
                                        funcname
    type( c_ptr ) :: ptr
    type( sub_t ) :: evt
  end type

end module type_def
! ...................................................................................
