module routines
 implicit none
 public

 contains

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

  subroutine read_latvecs(fd,latvecs)
   implicit none
   integer :: i
   integer, intent(in) :: fd
   real, dimension(3,3), intent(out) :: latvecs
   
   do i=1,3
    read(fd,*) latvecs(i,1), latvecs(i,2), latvecs(i,3)
   end do
  end subroutine read_latvecs

  subroutine read_sites(fd,site_type, site_coords)
   implicit none
   integer, intent(in) :: fd  !! file descriptor (unit) number
   integer :: i, Nsites
   integer, allocatable, intent(out) :: site_type(:)
   real, allocatable, intent(out) :: site_coords(:,:)

   read(fd,*) Nsites
   allocate( site_type( 1:Nsites ) )
   allocate( site_coords( 1:Nsites, 1:3 ) )
   do i=1, Nsites
    read(fd,*) site_type(i), site_coords(i,1), site_coords(i,2), site_coords(i,3)
   end do
  end subroutine read_sites

  subroutine choose_p(G,d,rnd,idx)
  !-----------------------------------------
  ! choose some event from the list of probabilities which is summed up
  ! G1--G2--G3----G4-----G5-----G6------------G7
  ! "should they be sorted into some order?"
  ! ----------------------------------------
  ! G   ==> vector containing the probability values of events
  ! d   ==> dimension of vector G
  ! rnd ==> input random number
  ! idx ==> output index of the event chosen
  !-----------------------------------------
   implicit none
   integer, intent(in) :: d
   real, dimension(d), intent(in) :: G
   real, intent(inout) :: rnd
   integer, intent(out) :: idx

   real :: k
   integer :: i
   
   k=0
   rnd=rnd*sum(G)
   idx=0
   do i=1,d
     if ( k < rnd ) then
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


end module routines
