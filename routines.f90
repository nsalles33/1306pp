module routines
 implicit none
 public

 contains


  subroutine read_line(fd, line, end_of_file)
  !--------------------
  ! read a line, makes possible to use # for comment lines, skips empty lines, 
  !  is pretty much a copy from QE.
  !
   implicit none
   integer, intent(in) :: fd
   character(len=*), intent(out) :: line
   logical, optional, intent(out) :: end_of_file
   logical :: tend

   tend = .false.
101   read(fd,fmt='(A256)',END=111) line
      if(line == ' ' .or. line(1:1) == '#') go to 101
      go to 105
111   tend = .true.
      go to 105
105   continue

      if( present(end_of_file)) then
        end_of_file = tend
      endif
  end subroutine read_line


  subroutine get_nsites(fd_sites,nsites)
  !---------------------
  ! get the total number of sites
  !---------------------
   implicit none
   integer, intent(in) :: fd_sites
   integer, intent(out) :: nsites
   logical :: eof
   character(len=64) :: line

    eof=.false.
    do while (.not.eof)
      call read_line(fd_sites,line,eof)
      line = trim(adjustl(line))
      if(line(1:6)=='nsites') then
        line=trim(adjustl(line(7:)))
        if (line(1:1)=='=') line=trim(adjustl(line(2:)))
        read(line,*) nsites
        eof=.true.
      endif
    end do
    rewind(fd_sites)
  end subroutine get_nsites


  subroutine get_nevt(fd_events,nevt)
  !-------------------
  ! get total number of events
  !-------------------
   implicit none
   integer, intent(in) :: fd_events
   integer, intent(out) :: nevt
   logical :: eof
   character(len=64) :: line

   eof = .false.
   do while ( .not. eof )
     call read_line( fd_events, line, eof )
     line = trim ( adjustl ( line ) )
     if ( line(1:4) == 'nevt' ) then
       line = trim ( adjustl ( line(5:) ) )
       if (line(1:1)=='=') line=trim(adjustl(line(2:)))
       read(line,*) nevt
       eof = .true.
     endif
   end do
   rewind(fd_events)
  end subroutine get_nevt


  subroutine get_hash_prob(fd,hash1,hash2,prob,nevt)
  !-----------------------
  ! parse through the event input file and extract the total number of events.
  ! It is given by 'nevt #', probably in the first line.
  !
  ! In the second run, Parse through events input and extract the initial and
  !  final hash of an event, and
  ! the probability of this event. The line with this info should be right under the
  ! line with '@number' of an event, keep this number as ievt (index of event).
  !-----------------------
   implicit none
   integer, intent(in) :: fd
   integer, optional, intent(out) :: nevt
   integer, allocatable, intent(out) :: hash1(:), hash2(:)
   real, allocatable, intent(out) :: prob(:)
   integer :: ievt
   character(len=256) :: line
   logical :: eof
 
   eof = .false.
   do while ( .not. eof )
     call read_line( fd, line, eof )
     line = trim ( adjustl ( line ) )   !! to get rid of surrounding spaces
     if ( line(1:4) == 'nevt' ) then
       line = trim ( adjustl ( line(5:) ) )  !! get rid of 'nevt'
       if (line(1:1)=='=') line=trim(adjustl(line(2:)))
       read(line,*) ievt  !! get number from character
       eof = .true.
     endif
   end do
   rewind(fd)
   
   if( present(nevt)) nevt=ievt

   allocate(hash1(1:nevt))
   allocate(hash2(1:nevt))
   allocate(prob(1:nevt))

   eof = .false.
   do while ( .not. eof )
     call read_line( fd, line, eof )
     line = trim ( adjustl ( line ) ) 
     if ( line(1:1) == '@' ) then
       line = trim ( adjustl ( line(2:) ) )   !! to get rid of '@' and possible spaces
       read(line,*) ievt
       read(fd,*) hash1(ievt), hash2(ievt), prob(ievt)
     endif
   end do
   rewind(fd)
  end subroutine get_hash_prob


  subroutine get_ev_coord( fd, ev_idx, ev_init_nat, ev_init_typ, ev_init_coord, &
                                       ev_final_nat, ev_final_typ, ev_final_coord )
  !-------------------------------------
  ! extract the initial and final coordinates of the chosen event
  !-------------------------------------
   implicit none
   integer, intent(in) :: fd
   integer, intent(in) :: ev_idx
   character(len=256) :: line
   logical :: eof
   integer :: ievt, i
   integer, intent(out) :: ev_init_nat, ev_final_nat
   integer, allocatable, intent(out) :: ev_init_typ(:), ev_final_typ(:)
   real, allocatable, intent(out) :: ev_init_coord(:,:), ev_final_coord(:,:)

   eof=.false.
   do while (.not. eof)
     call read_line(fd,line,eof)
     line = trim ( adjustl (line) )
     if ( line(1:1) == '@' ) read(line(2:),*) ievt
     !!! get the wanted event given by ev_idx
     if (ievt == ev_idx) then
       do while (.not.eof)
         call read_line(fd,line,eof)
         line = trim ( adjustl (line) )
         if ( line(1:13) =='begin initial' ) then
           !!! read initial configuration
           read(fd,*) ev_init_nat
           allocate(ev_init_typ(1:ev_init_nat))
           allocate(ev_init_coord(1:ev_init_nat,1:3))
           do i=1,ev_init_nat
             read(fd,*) ev_init_typ(i),ev_init_coord(i,1), &
                           ev_init_coord(i,2), ev_init_coord(i,3)
           end do
         elseif( line(1:11)=='begin final') then
           !!! read final configuration
           read(fd,*) ev_final_nat
           allocate(ev_final_typ(1:ev_final_nat))
           allocate(ev_final_coord(1:ev_final_nat,1:3))
           do i=1,ev_final_nat
             read(fd,*) ev_final_typ(i), ev_final_coord(i,1),&
                           ev_final_coord(i,2), ev_final_coord(i,3)
           end do
           !!! finished reading all necessary, break the loop
           eof=.true.
         endif
       end do
     endif
   end do
   rewind(fd)
  end subroutine get_ev_coord


  subroutine periodic(c)
  !--------------------------------
  ! periodic boundary conditions in "crystal" coordinates
  !--------------------------------
   implicit none
   real :: c
   if(c < (-0.5)) c = c + 1.0
   if(c >= 0.5) c = c - 1.0

   return
  end subroutine periodic


  subroutine set_random_seed()
   ! ----- setting a random seed, based on current time -----
   INTEGER :: i_seed
   INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
   INTEGER, DIMENSION(1:8) :: dt_seed

   ! ----- Set up random seed portably -----
   CALL RANDOM_SEED(size=i_seed)
   ALLOCATE(a_seed(1:i_seed))
   CALL RANDOM_SEED(get=a_seed)
   CALL DATE_AND_TIME(values=dt_seed)
   a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
   CALL RANDOM_SEED(put=a_seed)
   DEALLOCATE(a_seed)
  end subroutine set_random_seed


  subroutine choose_p(G,d,rnd,idx)
  !-----------------------------------------
  ! choose some event from the list of probabilities which is summed up
  ! G1--G2--G3----G4-----G5-----G6------------G7
  ! "should the probabilities be sorted into some order?"
  ! ----------------------------------------
  ! G   ==> vector containing the probability values of events
  ! d   ==> dimension of vector G
  ! rnd ==> input random number
  ! idx ==> output index of the event chosen
  ! -------------
  ! k    ==> dummy sum of the smaller G values
  ! rnd1 ==> rnd scaled by sum(G), for testing remove this and put rnd intent(inout)
  !-----------------------------------------
   implicit none
   integer, intent(in) :: d
   real, dimension(d), intent(in) :: G
   real, intent(in) :: rnd
   integer, intent(out) :: idx

   real :: k,rnd1
   integer :: i
   
!! do we normalize the G vector by sum(G)? or it's ok like this? 
!! now it's opposite, the rnd is scaled by sum(G)... 

   k=0
   rnd1=rnd*sum(G)
   idx=0
   do i=1,d
     if ( k < rnd1 ) then
       idx = idx + 1
       k = k + G(i)
     else
       exit
     endif
   end do
  end subroutine choose_p




!!! ------------------------- !!!
!                               !
!    3-dimensional routines     ! 
!                               !
!!! ------------------------- !!!

  subroutine get_tf3D(r2, r1, A)
  !------------------------
  ! get the transformation matrix from r2 to r1 such that
  !       A r2 = r1
  !------------------------
   implicit none
   real, dimension(3), intent(in) :: r2
   real, dimension(3), intent(in) :: r1
   real, dimension(3,3), intent(out) :: A

   A(1,1) = r1(1)*r2(1)
   A(1,2) = r1(1)*r2(2)
   A(1,3) = r1(1)*r2(3)
   A(2,1) = r1(2)*r2(1)
   A(2,2) = r1(2)*r2(2)
   A(2,3) = r1(2)*r2(3)
   A(3,1) = r1(3)*r2(1)
   A(3,2) = r1(3)*r2(2)
   A(3,3) = r1(3)*r2(3)
 
  end subroutine get_tf3D

 
  subroutine get_tf3D_short(r1,r2,A)
  !--------------------------
  ! Transformation matrix from r1 to r2
  ! A r1 = r2
  !--------------------------
   implicit none
   real, dimension(3), intent(in) :: r1, r2
   real, dimension(3,3), intent(out) :: A

   A = spread (r2, dim=2, ncopies=3 ) * spread( r1, dim=1, ncopies=3 )
  end subroutine get_tf3D_short
 

  subroutine read_latvecs3D(fd,latvecs)
   implicit none
   integer :: i
   integer, intent(in) :: fd
   real, dimension(3,3), intent(out) :: latvecs
   
   do i=1,3
    read(fd,*) latvecs(i,1), latvecs(i,2), latvecs(i,3)
   end do
  end subroutine read_latvecs3D


  subroutine read_sites3D(fd,site_type, site_coords, Nsites)
  !---------------------
  ! read input file of the 3D sites, file structure:
  ! 
  ! integer(number of sites)
  ! integer(site_type) coord_x coord_y coord_z
  !
   implicit none
   integer, intent(in) :: fd  !! file descriptor (unit) number
   integer :: i
   integer, intent(out) :: Nsites
   integer, allocatable, intent(out) :: site_type(:)
   real, allocatable, intent(out) :: site_coords(:,:)

   read(fd,*) Nsites
   allocate( site_type( 1:Nsites ) )
   allocate( site_coords( 1:Nsites, 1:3 ) )
   do i=1, Nsites
    read(fd,*) site_type(i), site_coords(i,1), site_coords(i,2), site_coords(i,3)
   end do
  end subroutine read_sites3D


  subroutine read_sites3D_new(fd_sites, nsites,site_hash, coord)
   implicit none
   integer, intent(in) :: fd_sites
   integer, allocatable, intent(out) :: site_hash(:)
   real, allocatable, intent(out) :: coord(:,:)
   integer, intent(in) :: nsites
   integer :: isite
   character(len=64) :: line
   logical :: eof
   
   if (.not. allocated(site_hash)) allocate(site_hash(1:nsites))
   if (.not. allocated(coord)) allocate(coord(1:nsites,1:3))
    eof=.false.
    do while (.not.eof)
      call read_line(fd_sites,line,eof)
      line=trim(adjustl(line))
      if (line(1:5)=='begin') then
        do isite=1,nsites
          read(fd_sites,*) site_hash(isite), coord(isite,1), coord(isite,2), coord(isite,3)
        end do
        eof=.true.
       endif
    end do
    rewind(fd_sites)
  end subroutine read_sites3D_new


  subroutine pbcdist3D(A,B,Lx,Ly,Lz,dist)
  !--------------------------------
  ! distance between two points in 3-dimensions
  ! A, B         ==> vectors of two points
  ! Lx, Ly, Lz   ==> x,y,z sizes of the box
  ! dist         ==> output distance
  !--------------------------------
   implicit none
   real, dimension(3), intent(in) :: A, B
   real, intent(in) :: Lx, Ly, Lz
   real, intent(out) :: dist
   real :: sqr
 
    sqr = ( B(1) - A(1) - Lx*nint(( B(1) - A(1)) / Lx ))**2 + &
          ( B(2) - A(2) - Ly*nint(( B(2) - A(2)) / Ly ))**2 + & 
          ( B(3) - A(3) - Lz*nint(( B(3) - A(3)) / Lz ))**2
    dist = sqr**0.5 
  end subroutine pbcdist3D
  
  
  subroutine cart_to_crist3D(xpp,ct)
  !----------------------------
  ! cartesian to crystallographic coordinates transform, in 3-dimension
  ! v_crist = B^-1 * R_cart; where B is the matrix formed by unit cell vectors
  ! --------
  ! xpp(3)      ==> input vector of position in cartesian
  ! ct(3,3)     ==> conversion matrix, vectors of the Bravais lattice
  !----------------------------
  ! bt(3,3) ==> inverse matrix of ct, used locally
  ! xc(3)   ==> copy of xpp, used locally
  ! detct   ==> determinant of ct, used locally
  !
   implicit none
   real, dimension(3), intent(inout) :: xpp
   real, dimension(3,3), intent(in) :: ct

   real,dimension(3) :: xc
   real :: detct
   real, dimension(3,3) :: bt

       bt(:,:)=0.0
       xc(:) = 0.0       

  ! -----------------------------------------------
  !  inverse matrix of ct(:,:)
  !------------------------------------------------
       detct=ct(1,1)*ct(2,2)*ct(3,3)+&
             ct(1,2)*ct(2,3)*ct(3,1)+&
             ct(2,1)*ct(3,2)*ct(1,3)&
            -ct(1,3)*ct(2,2)*ct(3,1)&
            -ct(3,2)*ct(2,3)*ct(1,1)&
            -ct(1,2)*ct(2,1)*ct(3,3)

       bt(1,1)= ct(2,2)*ct(3,3)-ct(2,3)*ct(3,2)
       bt(1,2)=-(ct(1,2)*ct(3,3)-ct(1,3)*ct(3,2))
       bt(1,3)= ct(1,2)*ct(2,3)-ct(1,3)*ct(2,2)
       bt(2,1)=-(ct(2,1)*ct(3,3)-ct(2,3)*ct(3,1))
       bt(2,2)= ct(1,1)*ct(3,3)-ct(3,1)*ct(1,3)
       bt(2,3)=-(ct(1,1)*ct(2,3)-ct(1,3)*ct(2,1))
       bt(3,1)= ct(2,1)*ct(3,2)-ct(2,2)*ct(3,1)
       bt(3,2)= -(ct(1,1)*ct(3,2)-ct(1,2)*ct(3,1))
       bt(3,3)= ct(1,1)*ct(2,2)-ct(2,1)*ct(1,2)
  !------------------------------------------------

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))/detct
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))/detct
       xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))/detct

       xpp(:) = xc(:)

       return
  end subroutine cart_to_crist3D


  subroutine crist_to_cart3D(xpp,bt)
  !--------------------------------
  ! crystallographic to cartesian transformation in 3-dimensions
  ! R_cart = B * v_crist; where B is the matrix formed by cell vectors horizontally
  ! -----------
  ! xpp(3)    ==> input vector in crystallographic, output vector in cartesian
  ! bt(3,3)   ==> input conversion matrix, vectors of the Bravais lattice
  !-----
  ! xc(3)   ==> local vector
  !
   implicit none
   real, dimension(3), intent(inout) :: xpp
   real, dimension(3,3), intent(in) :: bt
   real, dimension(3) :: xc

       xc(:) = 0.0 

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))
       xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))

       xpp(:)=0.0

       xpp(:) = xc(:)

       return
  end subroutine crist_to_cart3D





!!! -------------------------- !!!
!                                !
!    2-dimensional routines      !
!                                !
!!! -------------------------- !!!


  subroutine read_latvecs2D(fd,latvecs)
   implicit none
   integer :: i
   integer, intent(in) :: fd
   real, dimension(2,2), intent(out) :: latvecs
   
   do i=1,2
    read(fd,*) latvecs(i,1), latvecs(i,2)
   end do
  end subroutine read_latvecs2D


  subroutine read_sites2D(fd,site_type, site_coords)
  !---------------------
  ! read input file of the 2D sites, file structure:
  ! 
  ! integer(number of sites)
  ! integer(site_type) coord_x coord_y
  !
   implicit none
   integer, intent(in) :: fd  !! file descriptor (unit) number
   integer :: i, Nsites
   integer, allocatable, intent(out) :: site_type(:)
   real, allocatable, intent(out) :: site_coords(:,:)

   read(fd,*) Nsites
   allocate( site_type( 1:Nsites ) )
   allocate( site_coords( 1:Nsites, 1:2 ) )
   do i=1, Nsites
    read(fd,*) site_type(i), site_coords(i,1), site_coords(i,2)
   end do
  end subroutine read_sites2D


  subroutine pbcdist2D(A,B,Lx,Ly,dist)
  !--------------------------------
  ! distance between two points in 2-dimensions
  ! A, B     ==> vectors of two points
  ! Lx, Ly   ==> x,y size of the box
  ! dist     ==> output distance
  !--------------------------------
   implicit none
   real, dimension(2), intent(in) :: A, B
   real, intent(in) :: Lx, Ly
   real, intent(out) :: dist
   real :: sqr
 
    sqr = ( B(1) - A(1) - Lx*nint(( B(1) - A(1)) / Lx ))**2 + &
          ( B(2) - A(2) - Ly*nint(( B(2) - A(2)) / Ly ))**2  
    dist = sqr**0.5 
  end subroutine pbcdist2D


  subroutine cart_to_crist2D(xpp,ct)
  !----------------------------
  ! cartesian to crystallographic coordinates transform, in 2-dimension
  ! v_crist = B^-1 * R_cart; where B is the matrix formed by unit cell vectors
  ! --------
  ! xpp(2)      ==> input vector of position in cartesian
  ! ct(2,2)     ==> conversion matrix, vectors of the Bravais lattice
  !----------------------------
  ! bt(2,2) ==> inverse matrix of ct, used locally
  ! xc(2)   ==> copy of xpp, used locally
  ! detct   ==> determinant of ct, used locally
  !
   implicit none
   real, dimension(2), intent(inout) :: xpp
   real, dimension(2,2), intent(in) :: ct

   real,dimension(2) :: xc
   real :: detct
   real, dimension(2,2) :: bt

       bt(:,:)=0.0
       xc(:) = 0.0       

  ! -----------------------------------------------
  !  inverse matrix of ct(:,:)
  !------------------------------------------------
       detct=ct(1,1)*ct(2,2) - ct(1,2)*ct(2,1)


       bt(1,1)= ct(2,2)
       bt(1,2)=-ct(1,2)
       bt(2,1)=-ct(2,1)
       bt(2,2)= ct(1,1)
  !------------------------------------------------

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1))/detct
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2))/detct

       xpp(:) = xc(:)

       return
  end subroutine cart_to_crist2D


  subroutine crist_to_cart2D(xpp,bt)
  !--------------------------------
  ! crystallographic to cartesian transformation in 2-dimensions 
  ! R_cart = B * v_crist; where B is the matrix formed by cell vectors horizontally
  ! -----------
  ! xpp(2)    ==> input vector in crystallographic, output vector in cartesian
  ! bt(2,2)   ==> input conversion matrix, vectors of the Bravais lattice
  !-----
  ! xc(2)   ==> local vector
  !
   implicit none
   real, dimension(2), intent(inout) :: xpp
   real, dimension(2,2), intent(in) :: bt
   real, dimension(2) :: xc

       xc(:) = 0.0   


       xc(1) = xpp(1)*bt(1,1)+xpp(2)*bt(2,1)
       xc(2) = xpp(1)*bt(1,2)+xpp(2)*bt(2,2)

       xpp(:)=0.0

       xpp(:) = xc(:)

       return
  end subroutine crist_to_cart2D


  subroutine get_tf2D(r1,r2,A)
  !-------------------
  ! get the transfer matrix from r1 to r2
  !   A r1 = r2
  !-------------------
   implicit none
   real, dimension(2), intent(in) :: r1,r2
   real, dimension(2,2), intent(out) :: A

    A = spread(r2,dim=2,ncopies=2) * spread(r1,dim=1,ncopies=2)
  end subroutine get_tf2D

end module routines
