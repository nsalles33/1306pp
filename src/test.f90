program test
implicit none

integer :: i,j,n

j = 10
do i = 0,50 
   n = modulo(i,j)
   write (*,*) i, "%",j,n
enddo

end program test
