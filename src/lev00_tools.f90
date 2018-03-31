! =================================================================================================
!     The Routine TOOLS BOX for everything 
!  In my tools box there are:
!  - read_line
!  - parse
!  - ...
! =================================================================================================
!
      subroutine read_line(fd, line, end_of_file)
      !--------------------
      ! read a line, makes possible to use # for comment lines, skips empty lines, 
      !  is pretty much a copy from QE.
      !---------
      ! fd ==> file descriptor
      ! line ==> what it reads
        implicit none
        integer, intent(in) :: fd
        integer             :: ios
        character(len=*), intent(out) :: line
        logical, optional, intent(out) :: end_of_file
        logical :: tend

        tend = .false.
101     read(fd,fmt='(A256)',END=111, iostat=ios) line
        if (ios /= 0) print*, " Reading Problem..."
        if(line == ' ' .or. line(1:1) == '#') go to 101
        go to 105
111     tend = .true.
        go to 105
105   continue

        if( present(end_of_file)) then
          end_of_file = tend
        endif
      end subroutine read_line
! .............................................................................

      subroutine parse(instring, delims, args, nargs)
!
      implicit none
!
      CHARACTER (len=*), intent (in) :: instring
      character (len=1)              :: delims,delim2
      character (len=500)            :: strgtmp
      CHARACTER (len=100), dimension (50), INTENT(INOUT) :: args
      integer, intent (out)                              :: nargs
      integer                                            :: indexs,leng,nmax,test
      character                                          :: tab

     ! allocate(args(nargs))
      args(:) = " "
      delim2 = '"'
      
      tab = char(9)
      !write (*,*) 'length of tab', len(tab)
      strgtmp = TRIM(adjustl(instring))
      !write (*,*) 'in ',strgtmp

      nargs = 0
      nmax = size(args)
      do
        leng = len(trim(strgtmp))
        !write (*,*) 'do',nargs+1,trim(strgtmp),':o'
        if (leng == 0) exit
        nargs = nargs + 1
        !write (*,*) 'length',leng,len(instring)
        if ( strgtmp(1:1) == delim2 ) then
           indexs = SCAN(strgtmp(1:leng),delim2,.true.)
           !write (*,*) "Element trouver ",strgtmp(1:1)," second ",indexs,leng
           args(nargs) = strgtmp(2:indexs-1)
           strgtmp = trim(adjustl(strgtmp(indexs+1:leng)))
           cycle
        elseif ( strgtmp(1:1) == tab ) then
           !write (*,*) " Tab trouver..."
           indexs = VERIFY(strgtmp(2:leng),tab)
           !write (*,*) " prochain non tab ", indexs
           strgtmp = trim(adjustl(strgtmp(indexs+1:leng)))
           !cycle
        !elseif ( SCAN(strgtmp,delims) > SCAN(strgtmp,tab) ) then
           !write (*,*) " Tab avant space"    
        endif
        indexs = SCAN(strgtmp,delims)
        test = SCAN(strgtmp,char(9))
        if ( test /= 0.and.test < indexs) indexs = test
        args(nargs)= strgtmp(1:indexs-1)
        strgtmp = trim(adjustl(strgtmp(indexs+1:leng)))
        if (nargs == nmax) exit
      enddo
      !
      !write (*,*) nargs,'out ',(trim(args(i)),i=1,nargs)
      end subroutine parse
! .............................................................................
 !
 !  subroutine link_int1_ptr( ptr_c, ptr_f, n ) bind( C )
 !    use iso_c_binding
 !    implicit none
 !    type( c_ptr ) :: ptr_c
 !    integer( c_int ), intent( in ) :: n
 !    integer( c_int ), dimension(:), pointer :: ptr_f
 !    call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
 !  end subroutine link_int1_ptr
!!............................................................................

 !  subroutine link_int2_ptr( ptr_c, ptr_f, n, m ) bind( C )
 !    use iso_c_binding
 !    implicit none
 !    type( c_ptr ) :: ptr_c
 !    integer( c_int ), intent( in ) :: n, m
 !    integer( c_int ), dimension(:,:), pointer :: ptr_f
 !    call c_f_pointer( ptr_c, ptr_f, shape=[n,m] )
 !  end subroutine link_int2_ptr
!!............................................................................

 !  subroutine link_real1_ptr( ptr_c, ptr_f, n ) bind( C )
 !    use iso_c_binding
 !    implicit none
 !    type( c_ptr ) :: ptr_c
 !    integer( c_int ), intent( in ) :: n
 !    real( c_double ), dimension(:), pointer :: ptr_f
 !    call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
 !  end subroutine link_real1_ptr
!!............................................................................

 !  subroutine link_real2_ptr( ptr_c, ptr_f, n, m ) bind( C )
 !    use iso_c_binding
 !    implicit none
 !    type( c_ptr ) :: ptr_c
 !    integer( c_int ), intent( in ) :: n,m
 !    real( c_double ), dimension(:,:), pointer :: ptr_f
 !    call c_f_pointer( ptr_c, ptr_f, shape=[ n,m ] )
 !  end subroutine link_real2_ptr
! ............................................................................



! .............................................................................

      subroutine sort1( n, vect )
!
!  Sorting by insertion => N^2 routine
!
        implicit none
        integer, intent ( in ):: n
        integer :: i, j, a
        integer, dimension(n), intent( inout ) :: vect

        do j = 2,n
           a = vect( j )
           do i = j - 1,1,-1
              if ( vect(i) <= a ) goto 10
              vect( i + 1 ) = vect(i)
           enddo
           i = 0
10         vect( i + 1 ) = a
        enddo

      end subroutine sort1
! .............................................................................
SUBROUTINE HPSORT(N,RA)
  real RA(N)
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
  else
    RRA=RA(IR)
    RA(IR)=RA(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  goto 10
END


! .............................................................................





















