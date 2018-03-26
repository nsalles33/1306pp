#include <string.h>
#include <stdio.h>
#include <stdlib.h>

char** parsing( char*, int );

void main() {

  int i = 0;
  char str[] ="- This, a sample string.";
  char **mot;
  
  mot = parsing( str );
  //parsing( str );

  while ( mot[i] != NULL ) {
    printf( " mot[%d] : %s\n", i++, mot[i] );
  }

}

char ** parsing( char *string, int len ) {
  int i = 0;
  char **pch;
  pch = malloc( len*sizeof(char) );
  pch[i] = strtok( string, "#\t " );
  while ( pch[i] != NULL ) {
    printf( "sparse %s \n",pch[i] );
    i++;
    pch[i] = strtok( NULL, "#\t " );
  }
  return pch;
}

