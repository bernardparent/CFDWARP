#include <exm.h>
#include "share.h"
#include <stdio.h>
#include <stdlib.h>

void test_pdma ( long numlines ) {
  EXM_pdmaline_t *pdma1, *pdma2;
  long line, cnt;
  double sum;

  srand ( 105 );
  pdma1 = ( EXM_pdmaline_t * ) malloc ( ( numlines ) * sizeof ( EXM_pdmaline_t ) );
  pdma2 = ( EXM_pdmaline_t * ) malloc ( ( numlines ) * sizeof ( EXM_pdmaline_t ) );

  for ( line = 0; line < numlines; line++ ) {
    pdma1[line].val[0] = ( double ) rand (  );
    pdma1[line].val[1] = ( double ) rand (  );
    pdma1[line].val[2] = ( double ) rand (  );
    pdma1[line].val[3] = ( double ) rand (  );
    pdma1[line].val[4] = ( double ) rand (  );
    pdma1[line].val[5] = ( double ) rand (  );
  }

  /* the following shouldn't be necessary for the pdma to work */
  if ( FALSE ) {
    pdma1[0].val[0] = 0.0;
    pdma1[0].val[1] = 0.0;
    pdma1[1].val[0] = 0.0;
    pdma1[numlines - 1].val[4] = 0.0;
    pdma1[numlines - 1].val[3] = 0.0;
    pdma1[numlines - 2].val[4] = 0.0;
  }

  for ( line = 0; line < numlines; line++ ) {
    for ( cnt = 0; cnt < 6; cnt++ )
      pdma2[line].val[cnt] = pdma1[line].val[cnt];
  }

  EXM_solve_PDMA ( pdma2, numlines );

  for ( line = 0; line < numlines; line++ ) {
    sum = pdma1[line].val[2] * pdma2[line + 0].val[5];
    if ( line >= 2 )
      sum += pdma1[line].val[0] * pdma2[line - 2].val[5];
    if ( line >= 1 )
      sum += pdma1[line].val[1] * pdma2[line - 1].val[5];
    if ( line <= numlines - 3 )
      sum += pdma1[line].val[4] * pdma2[line + 2].val[5];
    if ( line <= numlines - 2 )
      sum += pdma1[line].val[3] * pdma2[line + 1].val[5];
    printf ( "%ld  %E   %E \n", line, pdma1[line].val[5], sum );

  }

  free ( pdma1 );
  free ( pdma2 );

}


void test_long(void){
  long nm1,n;
  printf("\nTesting long type..\n");
  nm1=1;
  n=1;
  do {
    nm1=n;
    n*=2;
    printf("%ld ",n);
  } while(n>nm1);
  printf("\n");
}


void test_long_long(void){
  long long nm1,n;
  printf("\nTesting long long type..\n");
  nm1=1;
  n=1;
  do {
    nm1=n;
    n*=2;
    printf("%lld ",n);
  } while(n>nm1);
  printf("\n");
}


int main ( int argc, char **argv ) {
  test_pdma ( 20 );
  test_long ();
  test_long_long ();
  return ( EXIT_SUCCESS );
}
