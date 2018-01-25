program hashtest
 use routines
 integer :: nevt, ievt, hash, fd, isite, nsites, k, idx, this_site, i
 integer, allocatable :: init_hash(:), final_hash(:)
 integer, allocatable :: site_hash(:), ev_site(:)
 real :: c, rnd
 real, allocatable :: init_pos1(:,:), final_pos1(:,:), init_pos2(:,:), final_pos2(:,:)
 real, allocatable ::  coords(:,:)
 real, allocatable ::  coord(:,:)
 real, allocatable :: PP(:), G(:)
 logical :: eof
 character(len=256) :: line
 integer :: some_nmbr, ev_init_nat, ev_final_nat
 integer, allocatable :: hash1(:), hash2(:),ev_init_typ(:),ev_final_typ(:)
 real, allocatable :: prob(:), ev_init_coord(:,:),ev_final_coord(:,:)
 real, allocatable :: prob_M(:)
 real, dimension(3,3) :: A
 real, dimension(2,2) :: B
 real, dimension(2) :: r, v

 open(unit=444,file='events.in',status='old')
 open(unit=111, file='site.in',status='old')

 call set_random_seed()

! write(*,*) nevt
! write(*,*) prob(:)
! write(*,*) 'en'

! call read_sites3D(111,site_hash,coords,nsites) 

! write(999,*) nsites
! write(999,*)
! do isite=1,nsites 
!   write(999,*) site_hash(isite), coords(isite,1), coords(isite,2), coords(isite,3)
! end do
 

!!! get number of sites and number of events
 call get_nsites(111,nsites)
 call get_nevt(444,nevt)
 write(*,*) 'nsites',nsites,'nevt',nevt

!!! allocate coord matrix
 allocate(coord(1:nsites,1:3))

!!! allocate site_hash vector
 allocate(site_hash(1:nsites))

!!! allocate prob_M vector
 allocate(prob_M(1:nsites*nevt))

!!! allocate event site vector for keeping track of sites
 allocate(ev_site(1:nsites*nevt))

!!! read the sites and hashes
 call read_sites3D_new(111,nsites,site_hash,coord)

!!! get the hashes and probabilities of events
 call get_hash_prob(444,hash1,hash2,prob,nevt)

!!! get all events with the same initial hash as site_hash(isite), and put their
!!! probabilities into vector prob_M( prob ), probabilities for all events on all sites.
!!! at the same time, make another vector which keeps track of which site it is, 
!!! e.g. ev_site( 1, 1, 1, 2, 2, 3, 5 ) says the first 3 events are on site 1, second two
!!! on site 2, the sixth event is on site 3, and the seventh on site 5.
 prob_M(:) = 0.0
 ev_site(:) = 0
 k=1
 do isite=1,nsites
   do ievt=1,nevt
     if ( hash1(ievt)==site_hash(isite) ) then
       prob_M(k) = prob(ievt)
       ev_site(k) = isite
       k = k+1
     endif
   end do
 end do
 
 do i=1,size(prob_M)
   write(*,'(A4,F5.2)') 'prob',prob_M(i)
   write(*,'(A2,I2)') 'on',ev_site(i)
 end do
 

 call random_number(rnd)
 call random_number(rnd)
 call choose_p(prob_M,size(prob_M),rnd,idx)

 write(*,*) 'chosen event',idx,'with',prob_M(idx)
 write(*,*) 'whish should be at site',ev_site(idx)
write(*,*)
 call choose_p(prob,size(prob),rnd,idx)
 write(*,*) 'chosen event',idx
! write(*,*) init_hash(idx), final_hash(idx), PP(idx)
! write(*,*) init_pos1(idx,1),init_pos1(idx,2),init_pos1(idx,3)
! write(*,*) final_pos1(idx,1),final_pos1(idx,2),final_pos1(idx,3)
 
 call get_ev_coord(444,idx, ev_init_nat, ev_init_typ, ev_init_coord,&
                            ev_final_nat,ev_final_typ,ev_final_coord)
 write(*,*) ev_init_typ(:)
 write(*,*) ev_init_coord(1,:)
 write(*,*) ev_init_coord(2,:)
 write(*,*) ev_final_coord(1,:)
 write(*,*) ev_final_coord(2,:)

 r = (/ 1.0, 0.5 /) 
 v = (/ 0.4, 0.2 /)
 write(*,*)
 B= spread(r,dim=2,ncopies=2)*spread(v,dim=1,ncopies=2)
 do i=1,2
   write(*,'(F5.2,F5.2)') B(i,1), B(i,2)
 end do
 write(*,*)

! call get_tf3D(ev_final_coord(1,:),ev_init_coord(1,:),A)
! do i=1,3
!  write(*,'(F5.2,F5.2,F5.2)') A(i,:)
! end do
! write(*,*) 'init',ev_init_coord(1,:)
! write(*,*) 'final',ev_final_coord(1,:)
! r=matmul(A,ev_final_coord(1,:))
! write(*,*) 'A applied to final = initial',r
! 
! call get_tf3D_short(ev_final_coord(1,:),ev_init_coord(1,:),A)
! do i=1,3
! write(*,'(F5.2,F5.2,F5.2)') A(i,:)
! end do
! write(*,*) 'init',ev_init_coord(1,:)
! write(*,*) 'final',ev_final_coord(1,:)
! r=matmul(A,ev_final_coord(1,:))
! write(*,*) 'A applied to final = initial',r
! 
! 
! call get_tf3D(ev_init_coord(1,:),ev_final_coord(1,:),A)
! do i=1,3
!  write(*,'(F5.2,F5.2,F5.2)') A(i,:)
! end do
! write(*,*) 'init',ev_init_coord(1,:)
! write(*,*) 'final',ev_final_coord(1,:)
! r=matmul(A,ev_init_coord(1,:))
! write(*,*) 'A applied to initial = final',r
! 
! call get_tf3D_short(ev_init_coord(1,:),ev_final_coord(1,:),A)
! do i=1,3
!  write(*,'(F5.2,F5.2,F5.2)') A(i,:)
! end do
! write(*,*) 'init',ev_init_coord(1,:)
! write(*,*) 'final',ev_final_coord(1,:)
! r=matmul(A,ev_init_coord(1,:))
! write(*,*) 'A applied to initial = final',r

 write(*,*) coords(this_site,:)
 do isite=1,nsites
   if (coords(isite,1)==coords(this_site,1)+init_pos2(idx,1) &
       .and. coords(isite,2)==coords(this_site,2)+init_pos2(idx,2) &
       .and. coords(isite,3)==coords(this_site,3)+init_pos2(idx,3)) then 
     write(*,*) init_pos2(idx,:),isite
     coords(isite,:)=coords(isite,:)+(final_pos2(idx,:)-init_pos2(idx,:))
   endif
 coords(this_site,:)=coords(this_site,:) + (final_pos1(idx,:)-init_pos1(idx,:))
 end do

 write(999,*) nsites
 write(999,*)
 do isite=1,nsites 
   write(999,*) site_hash(isite), coords(isite,1), coords(isite,2), coords(isite,3)
 end do
end program
