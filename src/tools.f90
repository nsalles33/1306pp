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
      character (len=1) :: delims,delim2
      character (len=100) :: strgtmp
      CHARACTER (len=100), dimension (50), INTENT(INOUT) :: args
      integer, intent (out) :: nargs
      integer :: indexs,leng,nmax,test
      character :: tab

     ! allocate(args(nargs))
      args(:) = " "
      delim2 = '"'
      
      tab = char(9)
!      write (*,*) 'length of tab', len(tab)
      strgtmp = TRIM(adjustl(instring))
!      write (*,*) 'in ',strgtmp

      nargs = 0
      nmax = size(args)
      do
        leng = len(trim(strgtmp))
!        write (*,*) 'do',nargs+1,trim(strgtmp),':o'
        if (leng == 0) exit
        nargs = nargs + 1
!        write (*,*) 'length',leng,len(instring)
        if ( strgtmp(1:1) == delim2 ) then
           indexs = SCAN(strgtmp(1:leng),delim2,.true.)
!           write (*,*) "Element trouver ",strgtmp(1:1)," second ",indexs,leng
           args(nargs) = strgtmp(2:indexs-1)
           strgtmp = trim(adjustl(strgtmp(indexs+1:leng)))
           cycle
        elseif ( strgtmp(1:1) == tab ) then
!           write (*,*) " Tab trouver..."
           indexs = VERIFY(strgtmp(2:leng),tab)
!           write (*,*) " prochain non tab ", indexs
           strgtmp = trim(adjustl(strgtmp(indexs+1:leng)))
!           cycle
     !   elseif ( SCAN(strgtmp,delims) > SCAN(strgtmp,tab) ) then
     !      write (*,*) " Tab avant space"    
        endif
        indexs = SCAN(strgtmp,delims)
        test = SCAN(strgtmp,char(9))
        if ( test /= 0.and.test < indexs) indexs = test
        args(nargs)= strgtmp(1:indexs-1)
        strgtmp = trim(adjustl(strgtmp(indexs+1:leng)))
        if (nargs == nmax) exit
      enddo
!
!      write (*,*) nargs,'out ',(trim(args(i)),i=1,nargs)
      end subroutine parse
! .............................................................................

