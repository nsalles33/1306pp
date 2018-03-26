#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* ------- VARIABLE DECLARATION... 
   ------------------------------- */
struct event_type {

  int nevent, nbond;

  int *ptr_i_state, 
      *ptr_f_state;

  double *ptr_ebarrier, 
         *ptr_de,
         **ptr_ebond;
};

struct kmc_type {

  int bavard;
  int period[3], nsites[3];

  char algorithm[50], 
       input_file[50], 
       input_event[50], 
       libname[50];

  int tot_sites, 
      max_step, 
      sys_dim, 
      freq_write, 
      nprop,
      node_state;

  double sum_rate, 
         rand_rate, 
         time, 
         rand_time, 
         max_time, 
         temp, 
         kt, 
         per100, 
         f0;
  
  int     *ptr_site, *ptr_nneig, *ptr_nevt, **ptr_neig, **ptr_event_site;
  double  *ptr_rate, *ptr_prop, **ptr_event_rate;

  struct event_type event;
};

extern int *__hidden_table_MOD_h_i_state, *__hidden_table_MOD_h_f_state;
extern double *__hidden_table_MOD_h_ebarrier, *__hidden_table_MOD_h_de, **__hidden_table_MOD_h_ebond;


/* ------- FUNCTION DECLARATION...
   ------------------------------- */

void builder_event_type( struct event_type *, int * );


/* ------- FUNCTION DEFINITION... 
   ------------------------------ */

char ** parsing( char *string, int len ) {
  int i = 0;
  char **pch;
  pch = malloc( len*sizeof(char) );
  pch[i] = strtok( string, "\t " );
  while ( pch[i] != NULL ) {
   // printf( "sparse %d %s \n",i, pch[i] );
    i++;
    pch[i] = strtok( NULL, "\t " );
  }
  return pch;
}
// ...........................................................

int read_line( FILE *fp, struct event_type *event ) {

  char string[100], **word;
  char *b = string;
  size_t bufsize = 100;
  size_t nword;

  int nevent, nbond;

  nword = getline(&b, &bufsize, fp );
  //printf( "Pb de lecture... %d %s %d\n",nword, &string[0], strcmp(&string[0],"#"));

  while ( strcmp(&string[0],"#") >= 0 ) {
     //printf( "read input_event %zu char %s | %d \n", nword, string, strcmp(&string[0],"#") );

     word = parsing( string, bufsize );
     if ( !strcmp(word[0],"Number_of_event") ) {
        nevent = 2*atoi( word[1] );
        printf( " nevent %d\n ",nevent/2 );
        builder_event_type( event, &nevent );
        //break;
        return nevent;
     }
     if ( !strcmp( word[0], "Energy_Bond") ) {
        nbond = atoi( word[1] );
        //break;
        return nbond;
     }
     nword = getline(&b, &bufsize, fp );
     //if ( !strcmp(&string[0],"#") ) exit(1);
  }

  return 0;
}	
// ...........................................................

void builder_event_type( struct event_type *event, int *nevent ) {
  int int_size = sizeof( int );
  int double_size = sizeof( double );

  printf( " Enter in EVENT_type constructor...\n" );
  event->nevent = *nevent;

  __hidden_table_MOD_h_i_state = malloc( *nevent*int_size );
  __hidden_table_MOD_h_f_state = malloc( *nevent*int_size );
  __hidden_table_MOD_h_ebarrier = malloc( *nevent*double_size );
  __hidden_table_MOD_h_de = malloc( *nevent*double_size );

  __hidden_table_MOD_h_ebond = malloc( event->nbond*sizeof( double* ) );
  for (int i = 0; i < event->nbond; i++ ) 
    __hidden_table_MOD_h_ebond[i] = malloc( event->nbond*sizeof( double** ) );

  event->ptr_i_state = __hidden_table_MOD_h_i_state;
  event->ptr_f_state = __hidden_table_MOD_h_f_state;
  event->ptr_ebarrier = __hidden_table_MOD_h_ebarrier;
  event->ptr_de = __hidden_table_MOD_h_de;
  event->ptr_ebond = __hidden_table_MOD_h_ebond;

/*  for (int i = 0; i < event->nbond; i++ ) {
    event->ptr_ebond[i][i] = 0.0;
    printf( "ptr_ebond %f \n", event->ptr_ebond[i][i] );
  }
*/
  printf( " EVENT_type constructor DONE \n" );

}
// ...........................................................

void read_event( struct kmc_type *struc ) {

  char string[100], **word;
  char *b = string;
  int i, id, nevent, nbond, ibd, jbd;
  size_t bufsize = 100;
  size_t nword;

  FILE *fp;

  printf( " === PRINT STRUCTURE KMC ===\n");
  printf( " Periodicity %d %d %d \n",struc->period[0], struc->period[1], struc->period[2] );
  printf( " nsites %d %d %d\n", struc->nsites[0], struc->nsites[1], struc->nsites[2] );
  printf( " INPUT EVENT %s \n", struc->input_event ); 
  printf( " ALGORITHM   %s \n", struc->algorithm );
  printf( " INPUT FILE  %s \n", struc->input_file );
  printf( " LIBNAME     %s \n", struc->libname );

  for ( id = 0; id < 3; id ++) {
     printf( "%d site %d ", id, struc->ptr_site[id] ); //, struc->event.ptr_ebond[1][1] );
  }

  fp = fopen( struc->input_event, "r" );
  if ( fp == NULL ) perror( " PB open input_event...\n" );

  nword = getline(&b, &bufsize, fp );
  //printf( "Pb de lecture... %d %s %d\n",nword, &string[0], strcmp(&string[0],"#"));

  while ( strcmp(&string[0],"#") >= 0 ) {
     //printf( "read input_event %zu char %s | %d \n", nword, string, strcmp(&string[0],"#") );

     word = parsing( string, bufsize );
     if ( !strcmp(word[0],"Number_of_event") ) { 
        nevent = 2*atoi( word[1] );
        printf( " nevent %d\n ",nevent/2 );
        builder_event_type( &struc->event, &nevent );
        break;
     }
     nword = getline(&b, &bufsize, fp );
     //if ( !strcmp(&string[0],"#") ) exit(1);

  }
//  exit(1);
  /* ==== read the event ==== */
  nword = getline(&b, &bufsize, fp );
  //word = parsing( string, bufsize );
  int n = struc->event.nevent / 2 ;
  for ( id = 0; id < struc->event.nevent/2; id ++ ) {
     word = parsing( string, bufsize );

     struc->event.ptr_i_state[ id ] = atoi( word[1] );           
     struc->event.ptr_f_state[ id ] = atoi( word[2] );           
     struc->event.ptr_ebarrier[ id ] = atof( word[3] );
     struc->event.ptr_de[ id ] = atof( word[4] );
     printf( " read %d event %d %d %f %f \n", id, struc->event.ptr_i_state[ id ], struc->event.ptr_f_state[ id ], struc->event.ptr_ebarrier[ id ], struc->event.ptr_de[ id ] ); 

     struc->event.ptr_i_state[ id + n ] = struc->event.ptr_f_state[ id ]; 
     struc->event.ptr_f_state[ id + n ] = struc->event.ptr_i_state[ id ];
     struc->event.ptr_ebarrier[ id + n ] = struc->event.ptr_ebarrier[ id ] - struc->event.ptr_de[ id ];
     struc->event.ptr_de[ id + n ] = - struc->event.ptr_de[ id ];
     printf( " read %d event %d %d %f %f \n", id+n, struc->event.ptr_i_state[ id + n ], struc->event.ptr_f_state[ id + n ], struc->event.ptr_ebarrier[ id + n ], struc->event.ptr_de[ id + n ] ); 

     nword = getline(&b, &bufsize, fp );
    // word = parsing( string, bufsize );
  }

  /* ==== read the energy ==== */
  //nword = getline(&b, &bufsize, fp );
  while ( strcmp(&string[0],"#") >= 0 ) {
//    printf( "entre dans le while-2\n");
    word = parsing( string, bufsize );
//    printf( " %s %s \n",word[0], word[1] );
    if( !strcmp( word[0], "Energy_Bond") ) {
      nbond = atoi( word[1] );
//      printf( " dans energy_bond %s %d \n",word[0],nbond );
      break;
    }	    
    nword = getline(&b, &bufsize, fp );
  }    

  for ( id = 0; id < nbond; id++ ) {
     nword = getline(&b, &bufsize, fp );
     word = parsing( string, bufsize );
//     printf( " ebond %s %s %s %s --\n", word[0], word[1], word[2], word[3] );
     ibd = atoi( word[1] );
     jbd = atoi( word[2] );
//     printf( " %d %d --\n", ibd, jbd );
     struc->event.ptr_ebond[ibd][jbd] = atof( word[3] );
//     printf( " %d %d %f --\n", ibd, jbd, struc->event.ptr_ebond[ibd][jbd] );
  }	  

  fclose( fp );

 // exit(1);

}
/* ======
    subroutine read_event( struc ) bind( C )
      use iso_c_binding
      use derived_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      character (len=50,kind=c_char)    :: string
      integer( c_int )                  :: u0, i, id, ios, nevent
      logical                           :: EOF

      character (len=1,kind=c_char)  :: delims
      CHARACTER (len=100,kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs

      integer( c_int ), dimension(:), pointer :: init_state, final_state
      real( c_double ), dimension(:), pointer :: ebarrier, de

!  ::: Lecture of state and rate of each node
      print*, " INPUT_EVENT: ",trim(struc% input_event)
      open( newunit=u0, file=trim(struc% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )

      EOF = .false.
      do while ( .not.EOF )
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)

         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call builder_event_type( struc% event, nevent )
            exit
         endif

      enddo
!
      call link_int1_ptr( struc% event% ptr_i_state, init_state, struc% event% nevent)
      call link_int1_ptr( struc% event% ptr_f_state, final_state, struc% event% nevent)
      call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier, struc% event% nevent)
      call link_real1_ptr( struc% event% ptr_de, de, struc% event% nevent)
!
      do i = 1,struc% event% nevent
         read (u0,*) id, init_state(id), final_state(id), ebarrier(id), de(id)
         write (*,*) id, init_state(id), final_state(id), ebarrier(id), de(id)
      enddo

      close( u0 )
!
    end subroutine read_event
*/
// .................................................................................................

//void event_rate_calc( struct kmc_type *struc ) {
//  printf( " Hello in event_rate_calc...");
//}

void event_rate_calc( struct kmc_type *struc ) {

  int i, jn, j, kn, dbd, l, ntot;
  double kt, ebd, f0;

  int *bd, *ddb, *v;
  ntot = struc->tot_sites;
  bd = malloc( ntot*sizeof( int ) );
  ddb = malloc( ntot*sizeof( int ) );
  v = malloc( ntot*sizeof( int ) );


  kt = struc->kt;
  ebd = struc->event.ptr_ebond[ 1 ][ 1 ]; 
  f0 = struc-> f0;
  
  printf( " %d\n", struc->ptr_neig[0][0]);
  for ( i = 0; i < struc->tot_sites; i++ ) {
    struc->ptr_rate[i] = 0.0;
    for ( jn = 0; jn < struc->ptr_nneig[i]; jn++ ) {
      printf( " %d neig[%d]  \n",i,jn); //,struc->ptr_neig[jn][i] );
      //struc->ptr_event_rate[jn][i] = 0.0;
    }
      exit(1);
  }
  printf( " Hello in event_rate_calc...\n");
  exit(1);


//
//!  -- We define the event rate with the environment 
//!  1) How many dandling bond per vacancy  
  for ( i = 0; i < struc->tot_sites; i++ ) {

    if ( struc->ptr_site[ i ] == 1 ) continue;
    
//  -- Dandling bond
    ddb[i] = 0; v[i] = 0; bd[i] = 0;
    for ( jn = 0; jn < struc->ptr_nneig[i]; jn++ ) {
      j = struc->ptr_neig[jn][i];
      if ( j == 0 ) printf( " PB system_rate: %d %d %d\n", i,jn,j );
      struc->ptr_event_rate[jn][i] = 0.0;
      if ( struc->ptr_site[j] == 1 ) ddb[i]++;
      if ( struc->ptr_site[j] == 0 ) v[i]++;
    }
    if ( ddb[i] == 0 ) continue;
//
//    -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
//      ---------------------------------------------------------
//       dbd( 0 ) = v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2
//      ---------------------------------------------------------

     for ( jn = 0; jn < struc->ptr_nneig[i]; jn++ ) {
       j = struc->ptr_neig[jn][i];

       if ( struc->ptr_site[j] == 0 ) continue;
     }

  }
}

/*
!
!  -- We define the event rate with the environment 
!  1) How many dandling bond per vacancy
      do i = 1,struc% tot_sites
         rate( i ) = 0.0

         if ( site( i ) == 1 ) cycle

!    -- Dandling bond...
         ddb( i ) = 0 ; v( i ) = 0 ; bd( i ) = 0
         do jn = 1,nneig( i )
            j = neig( jn, i )
            if (j == 0) write (*,*) "PB system_rate: ",i,jn,j
            event_rate( jn, i ) = 0.0
            if ( site( j ) == 1 ) ddb( i ) = ddb( i ) + 1
            if ( site( j ) == 0 ) v( i ) = v( i ) + 1
         enddo
         if ( ddb( i ) == 0 ) cycle
!
!    -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
!       ---------------------------------------------------------
!       dbd( 0 ) = v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2
!       ---------------------------------------------------------
!         write (*,*) i, struc% nneig( i ), bd( i ), ddb( i ), v( i )
         do jn = 1,nneig( i )  ! *** neighbour is 1
            j = neig( jn, i )

            if ( site( j ) == 0 ) cycle

            ddb( j ) = 0 ; bd( j ) = 0 ; v( i ) = 0
            do kn = 1,nneig(j)
               k = neig( kn, j )
               if (k <= 0.or.k > struc% tot_sites)  &
                 write (*,*) k,"PB system_rate",i,jn,j,(l,neig(l,j),l=1,nneig(j))
               if ( site( k ) == 0 ) ddb( j ) = ddb( j ) + 1
               if ( site( k ) == 1 ) bd( j ) = bd( j ) + 1
            enddo

            dbd = v( i ) + bd( j ) - ddb( i ) - ddb( j ) + 2
            event_rate( jn, i ) = f0*exp( - dbd*ebd/kt )
            rate( i ) = rate( i ) + event_rate( jn, i )

!            write (*,*) j, struc% nneig( j ), bd( j ), ddb( j ), v( j )
!            write (*,*) " dbd :",v( i ),'+',bd( j ),'-',ddb( i ),'-',ddb( j ),'+ 2 =',dbd
!            write (*,*) i, jn, j, dbd,f0,kt,ebd, -dbd*ebd/kt, struc% event_rate( jn, i )

         enddo
!         write (*,*) " *** ", i, f0, ebd, struc% rate( i )
!

      enddo

    end subroutine event_rate_calc
*/
// ..................................................................................................

void choose_event( struct kmc_type struc, int isite, int ievent ) {

}

/*
    subroutine choose_event( struc, isite, ievent ) bind( C )
      use iso_c_binding
      use derived_types
      use random
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ), intent( inout ) :: isite, ievent

      integer( c_int ) :: i, jn
      real( c_double ) :: rsum, rdn, rrdn

      integer( c_int ), dimension(:), pointer :: nneig
      real( c_double ), dimension(:,:), pointer :: event_rate
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_real2_ptr( struc% ptr_event_rate, event_rate, 10, struc% tot_sites )

!      print*, " - To do :: choose_event "
      isite = 0 ; ievent = 0
      call random_number( rdn )
      rrdn = rdn*struc% sum_rate
      struc% rand_rate = rrdn
!      
      rsum = 0.0
      do i = 1,struc% tot_sites
      !
         do jn = 1,nneig( i )
      !  !
            ievent = jn
            rsum = rsum + event_rate( jn, i )
!            write (*,*) "rsum:",i ,jn ,rsum, struc% event_rate( jn, i )
            if ( rsum > rrdn ) exit
      !  !
         enddo
      !
         isite = i
         if ( rsum > rrdn ) exit
      enddo
      if ( rsum <= rrdn ) write (*,*) " PB choose event...", rsum,struc% sum_rate
!      write (*,*) struc% sum_rate,rdn,rrdn,rsum
!      write (*,'(2(a,1x,i6))') " We Choose on site ", isite," event ",ievent
!      write (*,*) struc% site( isite ), struc% site( struc% neig( ievent,isite ) )

    end subroutine choose_event
*/
// ..................................................................................................

void event_applied( struct kmc_type struc, int is, int jn ) {

}

/*  
    subroutine event_applied( struc, is, jn ) bind( C )
      use iso_c_binding
      use derived_types
      implicit none
      type( KMC_type ) :: struc
      integer( c_int ), intent( in ) :: is, jn
      integer( c_int ) :: j
!      print*, " - To do :: Event_applied "

      integer( c_int ), dimension(:), pointer :: site, nneig
      integer( c_int ), dimension(:,:), pointer :: neig
      real( c_double ), dimension(:), pointer :: rate 
      real( c_double ), dimension(:,:), pointer :: event_rate 
      call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int2_ptr( struc% ptr_neig, neig, 10, struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      call link_real2_ptr( struc% ptr_event_rate, event_rate, 10, struc% tot_sites )

!  ::: Flip 0 -> 1
      if ( site( is ) == 1 ) then
         call warning( " site Flip Problem 0 = 1 ")
         write (*,*) is, site( is ), ( neig(j,is), site( neig(j,is) ),j=1,nneig(is) )
         write (*,*) is, rate( is ), ( neig(j,is), event_rate( j,is ),j=1,nneig(is) )
      endif
      site( is ) = 1
      j = neig( jn, is )
      if ( site( j ) == 0 ) call warning( " site Flip Problem 1 = 0 ")
      site( j ) = 0

    end subroutine event_applied    
*/
// ..................................................................................................

void analyse( struct kmc_type obj ) {

}

/*
    subroutine analyse( obj ) bind( C )
      use iso_c_binding
      use derived_types
      implicit none

      type( KMC_type ) :: obj
      integer( c_int ) :: i, j, jv, k, kv, ngp, gpv, nc, mixgp, nvac, maxsize, max_gp

      integer( c_int ), dimension( obj% tot_sites ) :: gp, histo
      integer( c_int ), dimension( -1:int(obj% tot_sites/2)) :: clster
      !
      integer( c_int ), dimension(:), pointer :: site, nneig
      integer( c_int ), dimension(:,:), pointer :: neig
      real( c_double ), dimension(:), pointer :: prop
      call link_int1_ptr( obj% ptr_site, site, obj% tot_sites )
      call link_int1_ptr( obj% ptr_nneig, nneig, obj% tot_sites )
      call link_int2_ptr( obj% ptr_neig, neig, 10, obj% tot_sites )
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )
      !
      ngp = 0 ; gp = 0 ; max_gp = 0
      clster = 0 ; nvac = 0
      do i = 1,obj% tot_sites
         if ( site( i ) == 1 ) cycle
         !
         nvac = nvac + 1
         nc = 0
         gp( i ) = -1
         gpv = 0 ; mixgp = 0
         do jv = 1,nneig( i )
            j = neig( jv, i )
            !
            if ( site( j ) == 1 ) cycle
            !if ( gp( j ) == 0 ) gpn = 1 
            !
            if ( gp( j ) /= 0.and.gpv == 0 ) then
               gpv = gp( j )
            elseif ( gp( j ) /= 0.and.gpv /= 0.and.gp(j) /= gpv) then
               mixgp = gp( j )
            endif
            !
            nc = nc + 1
            !gp( j ) = ngp           
            !
         enddo
       !  write (*,*) i,nc,ngp,gpv,mixgp
         !
         if ( nc /= 0.and.gpv /= 0.and.mixgp == 0 ) then
            gp( i ) = gpv
         elseif ( nc /= 0.and.gpv == 0 ) then
            ngp = ngp + 1
            gp( i ) = ngp
            do jv = 1,nneig( i )
               j = neig( jv, i )
               if ( site( j ) == 0 ) gp( j ) = gp( i )
            enddo
         elseif ( mixgp /= 0 ) then
            !
            gp( i ) = min( mixgp, gpv )
            do j = 1,i-1
               if ( site( j ) == 1 ) cycle
               if ( gp( j ) == mixgp.or.gp( j ) == gpv ) then
                  gp( j ) = gp( i )
                  do kv = 1,nneig( j )
                     k = neig( kv, j )
                     if ( site( k ) == 0 ) gp( k ) = gp( i )
                  enddo
               endif
            enddo
            !
         endif
         !
         clster( gp(i) ) = clster( gp(i) ) + 1
         max_gp = max( max_gp, gp(i) )
      !   write (*,*) i, gp(i), ngp, clster( gp(i) ), nvac
         !
      enddo
      !
      nvac = 0
      histo = 0 ; maxsize = 1
      histo( 1 ) = clster( -1 )
      do i = 1,max_gp
         nvac = nvac + clster( i )
         histo( clster(i) ) = histo( clster(i) ) + 1
         maxsize = max( maxsize, clster(i) )
      enddo
      prop( 1 ) = 0
      do i = 1,maxsize
         prop( 1 ) = prop( 1 ) + real(i*histo( i ))/real( clster(-1) + max_gp )
!         write (*,*) prop(1),i,histo(i), clster(-1) , max_gp
      enddo
      !
    end subroutine analyse
! ..................................................................................................
! ..................................................................................................

!  end module diff_vac

==== */





