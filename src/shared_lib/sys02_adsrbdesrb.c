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
         **ptr_ebond,
         *ptr_spec;
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
// .................................................................................
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
//.................................................................................

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

}
// .................................................................................................

void event_rate_calc( struct kmc_type *obj ) {
  
  int i, id, ievt, jevt, j0, jv, nvj;
  double kt, f0, eb, ebd;

  kt = obj->kt;
  f0 = obj->f0;

  for ( i= 0; i < obj->tot_sites; i++ ) {
      obj->ptr_rate[ i ] = 0.0;
      j0 = obj->ptr_nneig[ i ] - 1;
      nvj = obj->ptr_neig[ j0 ];
      for ( jv = 1; jv <= nvj; jv++ ) 
          obj->ptr_event_rate[ j0 + jv ];
  }

  for ( i= 0; i < obj->tot_sites; i++ ) {

      id   = obj->ptr_site[ i ];
      ievt = 0;
      j0   = obj->ptr_nneig[ i ] - 1;
      nvj = obj->ptr_neig[ j0 ];
      for ( jevt = 0; jevt < obj->event.nevent; jevt++ ) {
          if ( obj->event.ptr_i_state[ jevt ] == id ) {
             ievt = ievt + 1;
	     obj->ptr_event_site[ j0 + ievt ] = jevt;
          }
      }
      obj->ptr_event_site[ j0 ] = ievt;
      obj->ptr_nevt[ i ] = ievt;

      for ( jevt = 1; jevt <= obj->ptr_nevt[ i ]; jevt++ ) {
          ievt = obj->ptr_event_site[ j0 + jevt ];
          obj->ptr_event_rate[ j0 + jevt ] = obj->event.ptr_ebarrier[ ievt ]; //f0*exp( - obj->event.ebarrier[ ievt ]/kt )
          obj->ptr_rate[ i ] += obj->ptr_event_rate[ j0 + jevt ];
          //printf( " %d event %d %d %f\n", i,jevt, ievt, obj->ptr_event_rate[ j0 + jevt ] );
      }
      //printf( " *** %d %f\n", i, obj->ptr_rate[ i ] );

  } // -- For "i"

  //exit(1);
}
// .................................................................................................

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
// .................................................................................................

void event_applied( struct kmc_type *obj, int *is, int *jn ) {

  int jevt, j0;

  j0 = obj->ptr_nneig[ *is ] - 1;
  jevt = obj->ptr_event_site[ j0 + *jn ];

  if ( obj->ptr_site[ *is ] != obj->event.ptr_i_state[ jevt ] ) 
   printf( " Problem site state not correspond to initial state event\n" );

  obj->ptr_site[ *is ] = obj->event.ptr_f_state[ jevt ];
  //exit(1);
}
// .................................................................................................

void analyse( struct kmc_type *obj ) {

  int i;
  obj->ptr_prop[ 0 ] = 0.0;
  for ( i = 0; i < obj->tot_sites; i++ ) 
      if ( obj->ptr_site[ i ] == 1 ) obj->ptr_prop[ 0 ] += 1;

  obj->ptr_prop[ 0 ] /= (double)obj->tot_sites;

}
// .................................................................................................
















