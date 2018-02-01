!
!     Routine to create the list of neighbour
!
    subroutine neig_list( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc

      if ( struc% sys_dim == 2 ) call cubic_2D( struc )  
      if ( struc% sys_dim == 3 ) call cubic_3D( struc )  

    end subroutine neig_list
! =================================================================================================

    subroutine cubic_2D( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer :: i, j, x, nx, y, ny, nxy, id

      nx = struc% nsites(1)
      ny = struc% nsites(2)
      nxy = nx*ny
      write(*,*) " Neig_list DONE",nx,ny,ny*(nx-1)

      do i = 0,struc% tot_sites - 1
         id = i + 1
         struc% nneig( id ) = 4

         x = MODULO( i, nx )           ! ** begin at 1 because i begin at 1
         y = i / nx
         if ( x /= 0 )      struc% neig(1,id) = id - 1
         if ( x /= nx - 1 ) struc% neig(2,id) = id + 1
         if ( y /= 0 )      struc% neig(3,id) = id - nx
         if ( y /= ny - 1 ) struc% neig(4,id) = id + nx

!       -- Peridoc Bound Condition
         if (struc% period_x.and.x == 0)      struc% neig(1,id) = id + (nx - 1)
         if (struc% period_x.and.x == nx - 1) struc% neig(2,id) = id - (nx - 1)
         if (struc% period_y.and.y == 0)      struc% neig(3,id) = id + ny*(nx - 1)
         if (struc% period_y.and.y == ny - 1) struc% neig(4,id) = id - ny*(nx - 1)

         do j = 1,struc% nneig( id )
            if ( struc% neig( j, id ) == 0 ) then
               call warning( "neig_list : index neig = 0 " )
               write (*,*) " neig_list :",id,j,struc% neig(j,id)
            endif
         enddo

         !if (i < 4) &
         !  write (*,*) id,x,y,nx, " neig ", (struc% neig(j,id),j=1,4)
      enddo


    end subroutine cubic_2D
! =================================================================================================
    subroutine cubic_3D( struc )
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer :: i,x, nx, y, ny, z, nz, nxy, id, jn, tmp
      integer, dimension (:), allocatable :: atmp

      nx = struc% nsites(1)
      ny = struc% nsites(2)
      nz = struc% nsites(3)
      nxy = nx*ny
      write(*,*) " Neig_list DONE",nx,ny,ny*(nx-1),nz

      !
!  :: System 3D
        do i = 0, struc% tot_sites - 1
           id = i + 1
           struc% nneig( id ) = 6
           struc% neig( :, id ) = 0

           x = MODULO( i, nx )
           y = MODULO( i / nx, ny )
           z = i / nxy

           if ( x /= 0 )      struc% neig(1,id) = id - 1
           if ( x /= nx - 1 ) struc% neig(2,id) = id + 1
           if ( y /= 0 )      struc% neig(3,id) = id - nx
           if ( y /= ny - 1 ) struc% neig(4,id) = id + nx
           if ( z /= 0 )      struc% neig(5,id) = id - nxy
           if ( z /= nz - 1 ) struc% neig(6,id) = id + nxy
!
!         -- Periodic Bound Condition
           if ( struc% period_x ) then
              if ( x == 0 )      struc% neig(1,id) = id + (nx - 1)
              if ( x == nx - 1 ) struc% neig(2,id) = id - (nx - 1)
           endif
           if ( struc% period_y ) then
              if ( y == 0 )      struc% neig(3,id) = id + ny*(nx - 1)
              if ( y == ny - 1 ) struc% neig(4,id) = id - ny*(nx - 1)
           endif
           if ( struc% period_z ) then
              if ( z == 0 )      struc% neig(5,id) = id + nxy*(nz - 1)
              if ( z == nz - 1 ) struc% neig(6,id) = id - nxy*(nz - 1)
           endif

           atmp = pack( struc% neig(:,id), struc% neig(:,id) /= 0 )
           tmp = size( atmp )
           do jn = 1,struc% nneig( id )
              if (jn <= tmp) then
                 struc% neig( jn, id ) = atmp( jn )
              else
                 struc% neig( jn, id ) = 0
              endif
           enddo
           struc% nneig( id ) = tmp

         !  write (*,*) id,x,y,z,nxy,(nz-1),struc% nneig(id), "neig ", (struc% neig(j,id),j=1,struc% nneig(id))
        enddo
  
    end subroutine cubic_3D

! =================================================================================================





