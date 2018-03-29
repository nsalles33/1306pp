#include <stdio.h>
#include <stdlib.h>

/* --------- VARIABLE DECLARATION... */

struct sub_t {
  int nevt;
  double tab[10];
};

struct test_t {
  int n;
  int p[3];
  char libname[20], funcname[20];
  int *ptr;
  struct sub_t evt;
};

/* --------- FUNCTION DECLARATION... */

//void plus1( struct test_t * );

/* --------- FUNCTION DEFINITION... */

void plus1( struct test_t *truc ) {

  int i;

  printf( " libname: %s \n", truc->libname);
  printf( " funcname: %s \n", truc->funcname);
  printf( " n = %d\n",truc->n);
  printf( " p[0] %d p[1] %d p[2] %d\n", truc->p[0], truc->p[1], truc->p[2] );

  for ( i = 0; i < truc->n; i++ ) {
    truc->evt.tab[i] = (double)truc->ptr[i];
  }

  for ( i = 0; i < truc->n; i++ ) {
    truc->ptr[i] = truc->ptr[i] + 1;
    printf( " %f c-ptr+1 %d\n", truc->evt.tab[i], truc->ptr[i] );
  }

}

/*
  subroutine plus1( truc ) bind( C )
    use iso_c_binding
    use type_def
    implicit none
    type( test_t ), intent( inout ) :: truc
    integer( c_int ) :: i
    integer( c_int ), dimension(:), pointer :: ptr_a
    call link_ptr( truc% ptr, ptr_a, truc% n )
    do i = 1,truc% n
       ptr_a( i ) = ptr_a( i ) + 1
       write(*,*) 'ptr+1:',ptr_a(i)
    enddo
  end subroutine plus1
*/

