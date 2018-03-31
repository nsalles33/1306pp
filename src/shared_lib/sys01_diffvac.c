#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* ------- VARIABLE DECLARATION... 
   ------------------------------- */
struct event_type {

  int nevent, nbond, nchem_react;

  int *ptr_i_state, 
      *ptr_f_state;

  double *ptr_ebarrier, 
         *ptr_de,
         **ptr_ebond;
};

struct kmc_type {

  int bavard, conv;
  int period[3], nsites[3];

  char algorithm[50], 
       input_file[50], 
       input_event[50], 
       libname[50],
       init_mod[50];

  int tot_sites, 
      max_step, 
      sys_dim, 
      freq_write, 
      nprop,
      node_state,
      nspec;

  double sum_rate, 
         rand_rate, 
         time, 
         rand_time, 
         max_time, 
         temp, 
         kt, 
         per100, 
         f0;
  
  int     *ptr_site, *ptr_nneig, *ptr_nevt, *ptr_neig, *ptr_event_site, *ptr_spec;
  double  *ptr_rate, *ptr_prop, *ptr_event_rate;

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
// ......................................................................

double random_number( double min, double max ) {
   return min + (double)rand() / ( (double)(RAND_MAX + 1.0)*( max - min) );
}
// ......................................................................
// ......................................................................

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

   srand( time(NULL) );
 // exit(1);

}
// .................................................................................................

void event_rate_calc( struct kmc_type *struc ) {

  int i, j0, jn, j, k0, kn, k, nvj, nvk, dbd, l, ntot;
  double kt, eb, ebd, f0;

  int *bd, *ddb, *v;
  ntot = struc->tot_sites;
  bd = malloc( ntot*sizeof( int ) );
  ddb = malloc( ntot*sizeof( int ) );
  v = malloc( ntot*sizeof( int ) );


  kt = struc->kt;
  eb = struc->event.ptr_ebarrier[ 0 ];
  ebd = struc->event.ptr_ebond[ 1 ][ 1 ]; 
  f0 = struc-> f0;
  
  //printf( " %d\n", struc->ptr_neig[0][0]);
  for ( i = 0; i < struc->tot_sites; i++ ) {
    struc->ptr_rate[ i ] = 0.0;
    j0 = struc->ptr_nneig[ i ] - 1;
    for ( jn = 0; jn <= struc->ptr_neig[ j0 ]; jn++ ) {
      //printf( " %d %d neig[%d]  %d \n",i,j0,j0+jn,struc->ptr_neig[ j0 + jn ] );
      struc->ptr_event_rate[ j0+jn ] = 0.0;
    }
    //exit(1);
  }
  //printf( " Hello in event_rate_calc...\n");
  //exit(1);


//
//!  -- We define the event rate with the environment 
//!  1) How many dandling bond per vacancy  
  for ( i = 0; i < struc->tot_sites; i++ ) {

    if ( struc->ptr_site[ i ] == 1 ) continue;
    
//  -- Dandling bond
    ddb[i] = 0; v[i] = 0; bd[i] = 0;
    j0 = struc->ptr_nneig[ i ] - 1;
    nvj = struc->ptr_neig[ j0 ];
    for ( jn = 1; jn <= nvj ; jn++ ) {
      j = struc->ptr_neig[ j0 + jn ] - 1;
      if ( j< 0 ) printf( " PB system_rate: %d %d %d\n", i,jn,j );
      struc->ptr_event_rate[ j0 + jn ] = 0.0;
      if ( struc->ptr_site[ j ] == 1 ) ddb[i]++;
      if ( struc->ptr_site[ j ] == 0 ) v[i]++;
      //printf( " %d %d ddb %d v %d\n", i, j, ddb[i], v[i]);
    }
    //printf( " %d ddb %d v %d\n", i, ddb[i], v[i]);
    //exit(1);
    if ( ddb[i] == 0 ) continue;
//
//    -- Dandling bond change between site(i) = 0 and site(j) = 1 ...
//      ---------------------------------------------------------
//       dbd( 0 ) = - ( v( 0 ) + bd( 1 ) - ddb( 0 ) - ddb( 1 ) + 2 )
//      ---------------------------------------------------------
     j0 = struc->ptr_nneig[ i ] - 1;
     nvj = struc->ptr_neig[ j0 ];
     for ( jn = 1; jn <= nvj; jn++ ) {
         j = struc->ptr_neig[ j0 + jn ] - 1;

         if ( struc->ptr_site[j] == 0 ) continue;

         ddb[ j ] = 0; bd[ j ] = 0; v[ j ] = 0;
         k0 = struc->ptr_nneig[ j ] - 1;
         nvk = struc->ptr_neig[ k0 ];
         for ( kn = 1; kn <= nvk; kn++) {
             k = struc->ptr_neig[ k0 + kn ] - 1;
             //if ( k < 0 || k >= struc->tot_sites ) 
             //  printf( " %d PB system_rate %d %d %d \n ",k,i,jn,j );
             if ( struc->ptr_site[ k ] == 0 ) ddb[ j ] += 1;  
             if ( struc->ptr_site[ k ] == 1 ) bd[ j ] += 1;  
         }

         dbd =  ddb[ j ] + ddb[ i ] - v[ i ] - bd[ j ] - 2;
         if ( dbd <= 0 ) {
           struc->ptr_event_rate[ j0 + jn ] = f0*exp( -eb/kt );
         } else if (dbd > 0 ) {
           struc->ptr_event_rate[ j0 + jn ] = f0*exp( -(eb + dbd*ebd )/kt );
         }
         struc->ptr_rate[ i ] += struc->ptr_event_rate[ j0 + jn ];

         //printf( " %d + %d - %d - %d - 2 = %d \n",ddb[ j ], ddb[ i ], v[ i ], bd[ j ], dbd);
         //printf( " eb = %f | ebd = %f \n", eb, ebd );
         //printf( " %d %d %d event_rate[ %d ] %f\n",i,j,jn,j0+jn,struc->ptr_event_rate[ j0 + jn] );
     }
     //exit(1);
     //printf( " %d rate %f\n",i,struc->ptr_rate[ i ] ); 

  }

  free( ddb );
  free( v );
  free( bd );

  //exit(1);

}
// ....................................................................................
// ..................................................................................................

void choose_event( struct kmc_type *struc, int *isite, int *ievent ) {
   
   int i, j0, jn, j, nvj;

   *isite = 0;
   *ievent = 0;
   
   //srand( time(NULL) );

   double rdn = random_number( 0.0, 1.0 );
   double rrdn = rdn*struc->sum_rate;
   double *evt_rate = struc->ptr_event_rate;
   struc-> rand_rate = rrdn;

   //printf( " rdn %f rrdn %f\n", rdn, rrdn );

   double rsum = 0.0;
   for (i = 0; i < struc->tot_sites; i++) {
      
      j0 = struc->ptr_nneig[ i ] - 1;
      nvj = struc->ptr_neig[ j0 ];
      for ( jn = 1; jn <= nvj; jn++) {
          *ievent = jn;
          rsum += evt_rate[ j0 + jn ];
          if ( rsum > rrdn ) break;
      }
      *isite = i;
      if ( rsum > rrdn ) break;
   }
   if ( rsum <= rrdn ) 
     printf( " PB choose_event... %f %f %f\n", rrdn, rsum, struc->sum_rate);
   //printf( " We choose is %d  jn %d\n", *isite, *ievent );
   //exit(1);
}
// ..................................................................................................

void event_applied( struct kmc_type *struc, int *is, int *jn ) {

  int j, j0;
  
  //printf( " %d %d \n", *is, *jn);
  //exit(1);
  if ( struc->ptr_site[ *is ] == 1 ) {
    printf( " site Flip Problem 0 = 1 : %d %d \n", *is, struc->ptr_site[ *is ]);
  }
  struc->ptr_site[ *is ] = 1;
  j0 = struc->ptr_nneig[ *is ] - 1;
  j = struc->ptr_neig[ j0 + *jn ] - 1;
  if ( struc->ptr_site[ j ] == 0 ) 
    printf( " site Flip Problem 1 = 0 : %d %d \n", j, struc->ptr_site[ j ]);
  struc->ptr_site[ j ] = 0  ;
}
// ..................................................................................................

void analyse( struct kmc_type *obj ) {

  int i,j, j0, jv, jn,k0, k, kv, nvj, nvk, ngp, gpv, nc, mixgp, nvac, maxsize, max_gp;

  int gp[obj->tot_sites], histo[obj->tot_sites];
  int clster[obj->tot_sites];

  for ( i = 0; i < obj->tot_sites; i++) {
    clster[ i ] = 0;
    gp[ i ] = 0;
  }

  ngp = 1; max_gp = 0; nvac = 0;
  for ( i = 0; i < obj->tot_sites; i++ ) {

      if ( obj->ptr_site[ i ] == 1 ) continue;

      nvac += 1;
      nc = 0;
      gp[ i ] = 1;
      gpv = 0; mixgp = 0;
      j0 = obj->ptr_nneig[ i ] - 1;
      nvj = obj->ptr_neig[ j0 ];
      for ( jn = 1; jn <= nvj; jn++ ) {
          j = obj->ptr_neig[ j0 + jn ] - 1;
          
          if ( obj->ptr_site[ j ] == 1 ) continue; 

          if ( gp[ j ] != 0 && gpv == 0 ) {
             gpv = gp[ j ];
          } else if ( gp[ j ] != 0 && gpv != 0 && gp[ j ] != gpv ) {
             mixgp = gp[ j ];
          }

          nc += 1;
      } // ----

      if ( nc != 0 && gpv != 0 && mixgp == 0 ) {
         gp[ i ] = gpv;
      } else if ( nc != 0 && gpv == 0 ) {
         ngp += 1;
         gp[ i ] = ngp;
         j0 = obj->ptr_nneig[ i ] - 1;
         nvj = obj->ptr_neig[ j0 ];
         for ( jv = 1; jv <= nvj; jv++ ) {
             j = obj->ptr_neig[ j0 + jv ] - 1;
             if ( obj->ptr_site[ j ] == 0 ) gp[ j ] = gp[ i ];
         }
      } else if ( mixgp != 0 ) {
         gp[ i ] = fmin( mixgp, gpv );
         for ( j = 0; j < i; j++ ) {
             if ( obj->ptr_site[ j ] == 1 ) continue;
             if ( gp[ j ] == mixgp || gp[ j ] == gpv ) {
                gp[ j ] = gp[ i ];
                k0 = obj->ptr_nneig[ j ] - 1;
		nvk = obj->ptr_neig[ k0 ];
                for ( kv = 1; kv <= nvk; kv++ ) {
		    k = obj->ptr_neig[ k0 + kv ] - 1;
                    if ( obj->ptr_site[ k ] == 0 ) gp[ k ] = gp[ i ];
                } 
             }
         } 

      } // --- else if

      clster[ gp[ i ] ] += 1;
      max_gp = fmax( max_gp, gp[ i ] );
      histo[ i ] = 0; 
      
  } // ---- for "i"

  nvac = 0;
  maxsize = 1;
  for ( i = 0; i < max_gp; i++ ) {
      nvac += clster[ i ];
      histo[ clster[i] ] += 1;
      maxsize = fmax( maxsize, clster[i] ); 
  }
  obj->ptr_prop[ 0 ] = 0.;
  for ( i = 0; i < maxsize; i++ ) {
      obj->ptr_prop[ 0 ] += (double)( (i+1)*histo[i] )/(double)( clster[1] + max_gp );
  }
}

//! ..................................................................................................






