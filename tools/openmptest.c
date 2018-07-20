#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#ifndef round
#define round(a) (floor(a+0.5e0))
#endif

#define sqr(a)   ((a)*(a))
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))
#define pi  3.14159265358979323846

#define nf 12
#define imax 100
#define jmax 100
#define kmax 130

typedef double sqmat_t[nf][nf];

/* this function calls wfprintf with the same arguments
   as it is given and exits.*/
void fatal_error ( const char *formatstr, ... ) {
  va_list ap;
  char *newstr;
  newstr = ( char * ) malloc ( 10000 * sizeof ( char ) );
  fprintf ( stderr, "\n\n" );
  va_start ( ap, formatstr );
  vsprintf ( newstr, formatstr, ap );
  va_end ( ap );
  fprintf ( stderr, "%s", newstr );
  free ( newstr );
  fprintf ( stderr, "\nFatal error. Exiting.\n\n" );
  exit ( EXIT_FAILURE );
}

void set_matrix_to_identity ( sqmat_t A ) {
  long row, col;

  for ( row = 0; row < nf; row++ ) {
    for ( col = 0; col < nf; col++ ) {
      A[row][col] = 0.0e0;
    }
  }
  for ( row = 0; row < nf; row++ ) {
    A[row][row] = 1.0e0;
  }
}

void invert_matrix ( sqmat_t mattmp, sqmat_t mat2 ) {
  long row, row2, col;
  double multfact;
  sqmat_t mat1;
  for ( row = 0; row < nf; row++ ) {
    for ( col = 0; col < nf; col++ ) {
      mat1[row][col] = mattmp[row][col];
    }
  }
  set_matrix_to_identity ( mat2 );
  for ( row = 0; row < nf; row++ ) {
    for ( row2 = 0; row2 < nf; row2++ ) {
      if ( row2 != row ) {
        if ( mat1[row][row] == 0.0 )
          fatal_error ( "Matrix can not be inverted.\n" );
        multfact = -mat1[row2][row] / ( mat1[row][row] );

        for ( col = 0; col < nf; col++ ) {
          mat1[row2][col] = mat1[row2][col] + mat1[row][col] * multfact;
          mat2[row2][col] = mat2[row2][col] + mat2[row][col] * multfact;
        }
      }
    }
  }
  for ( row = 0; row < nf; row++ ) {
    for ( col = 0; col < nf; col++ ) {
      if ( mat1[row][row] == 0.0 )
        fatal_error ( "Matrix can not be inverted.\n" );
      mat2[row][col] = mat2[row][col] / ( mat1[row][row] );
    }
    mat1[row][row] = 1.0e0;
  }
}

void invert_matrix_row ( long j, double *sum ) {
  long i, k, row, col;
  sqmat_t mat[imax], matinv[imax];

  for ( k = 0; k < kmax; k++ ) {
    for ( i = 0; i < imax; i++ ) {
      // init mat
      for ( row = 0; row < nf; row++ ) {
        for ( col = 0; col < nf; col++ ) {
          mat[i][row][col] =
            ( ( double ) row * ( double ) col ) / ( 1.0 + ( double ) i + ( double ) j + ( double ) k );
          if ( row == col )
            mat[i][row][col] = 10.0;
        }
      }
      // invert mat
      invert_matrix ( mat[i], matinv[i] );
      // invert mat again
      invert_matrix ( matinv[i], mat[i] );
    }

    *sum = 0.0;
    for ( i = 0; i < imax; i++ ) {
      for ( row = 0; row < nf; row++ ) {
        for ( col = 0; col < nf; col++ ) {
          if ( row == col )
            *sum += mat[i][row][col];
        }
      }
    }
  }
}

void simplest_loop (  ) {
  long j;
  double sum[jmax];

#pragma omp parallel for private(j)
  for ( j = 0; j < jmax; j++ ) {
    invert_matrix_row ( j, &( sum[j] ) );
  }

  printf ( "sum=%E\n", sum[5] );
}

void simple_loop (  ) {
  long j, numthreads, thisthread, js, je;
  double sum[jmax];

#pragma omp parallel   default(shared) private(thisthread,numthreads,js,je,j)
  {
    thisthread = omp_get_thread_num (  );
    numthreads = omp_get_num_threads (  );
    js = round ( jmax * thisthread / numthreads );
    je = round ( jmax * ( thisthread + 1 ) / numthreads ) - 1;
    printf ( "jmax=%d thisthread=%ld numthreads=%ld js=%ld je=%ld\n", jmax, thisthread, numthreads, js, je );
    for ( j = js; j <= je; j++ ) {
      invert_matrix_row ( j, &( sum[j] ) );
    }
  }
  printf ( "sum=%E\n", sum[5] );
}

int main ( int argc, char *argv[] ) {
  //simplest_loop();
  simple_loop (  );
  return ( EXIT_SUCCESS );
}
