program main_test
 use routines
 real :: r
 real, dimension(3) :: A, B
 real :: Lx, Ly, Lz, dist, rnd, k
 real, dimension(3,3) :: BT,latvecs
 integer :: j,i,cnt,fd
 real, dimension(6) :: G,G1
 real, allocatable :: coords(:,:)
 integer, allocatable :: typ(:)

! --- test periodic -----
write(*,*) 'test periodic(r)'
 r=1.3
 write(*,*) ' r =',r
 call periodic(r)
 write(*,*) ' periodic(r) =',r
write(*,*)
!-----

! --- test pbcdist3D ----
write(*,*) 'test pbcdist3D'
 dist=0.0
 A(:)=0.0
 B(1)=1.0
 B(2)=1.0
 B(3)=1.0
 Lx=0.8
 Ly=5.0
 Lz=5.0
 call pbcdist3D(A,B,Lx,Ly,Lz,dist)
 write(*,'(A,F3.1,A,F3.1,A,F3.1,A,F8.4)') ' distance between (0, 0, 0) and (1, 1, 1) &
            in a box size ',Lx,' x ',Ly,' x ' ,Lz,'  is ',dist
write(*,*)
!------------------------

! --- test pbcdist2D ----
write(*,*) 'test pbcdist2D'
 dist=0.0
 A(:)=0.0
 B(1)=1.0
 B(2)=1.0
 Lx=0.3
 Ly=0.2
 call pbcdist2D(A,B,Lx,Ly,dist)
 write(*,'(A,F3.1,A,F3.1,A,F8.4)') ' distance between (0, 0) and (1, 1) &
            in a box size ',Lx,' x ',Ly,'  is ',dist
write(*,*)
!------------------------

! -- test crist_to_cart3D ----
write(*,*) 'test crist_to_cart'
BT(1,1)=-0.5
BT(1,2)=0.5
BT(1,3)=0.5
BT(2,1)=0.5
BT(2,2)=0.5
BT(2,3)=0.5
BT(3,1)=-0.5
BT(3,2)=-0.5
BT(3,3)=0.5
A(1)=0.0
A(2)=1.0
A(3)=1.0
write(*,'(A,F5.2,F5.2,F5.2)') ' A in crist',A(:)
call crist_to_cart3D(A,BT)
write(*,'(A,F5.2,F5.2,F5.2)') ' A in cart ',A(:)
call cart_to_crist3D(A,BT)
write(*,'(A,F5.2,F5.2,F5.2)') ' A in crist',A(:)
!-------------------------

write(*,*)
call set_random_seed()
G(1)=0.1
G(2)=0.1
G(3)=0.1
G(4)=0.2
G(5)=0.4
G(6)=0.4

!call random_number(rnd)
!rnd=rnd*sum(G)
!k=0
!cnt=0
!write(*,*) 'rnd',rnd
!do i=1,size(G)
!  if (k<rnd) then
!  cnt=cnt+1
!  k=k+G(i)
!  write(*,*) 'k',k
! else
!  exit
! endif
!end do
!write(*,*) G
!write(*,*) cnt,G(cnt)



call random_number(rnd)
write(*,*)
call choose_p(G,6,rnd,cnt)
write(*,*) rnd
write(*,*) cnt

!!!!!!!!!------------
write(*,*)
fd=901
open(unit=fd,file='sites.in',status='old')
call read_latvecs(fd,latvecs)
write(*,*) 'latvecs:'
do i=1,3
 write(*,*) latvecs(i,:)
end do

call read_sites(fd,typ,coords)
write(*,*) 'types',typ
write(*,*) 'coords:'
do i=1,size(coords)/3
  write(*,*) coords(i,:)
end do

end program main_test
