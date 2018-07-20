/* reference: G Chiocchia, Exact solutions to transonic and supersonic flows. Technical report AR-211, AGARD, 1985. */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <exm.h>

#define nj 25  //25, 49, 97
#define ni nj
#define MIRROR FALSE            /* mirror the solution with respect to the x axis? */

#define kmax 1.5
#define kmin 0.7
#define dq 1e-10
#define q0 0.5
#define gamma 1.4
#define Rgas 286.0

#define postfilename "ringleb.dat"
#define gridfilename "ringleb.grid"
#define initfilename "ringleb.init"

void Findxyfromkq ( double k, double q, double *x, double *y ) {
  double a, rho, J;
  a = sqrt ( 1.0 - ( gamma - 1 ) / 2.0 * sqr ( q ) );
  rho = pow ( a, 2.0 / ( gamma - 1.0 ) );
  J =
    1.0 / a + 1.0 / 3.0 / powint ( a, 3 ) + 1.0 / 5.0 / powint ( a,
                                                                 5 ) - 0.5 * log ( ( 1.0 + a ) / ( 1.0 -
                                                                                                   a ) );
  *x = 1.0 / 2.0 / rho * ( 2.0 / k / k - 1.0 / q / q ) - J / 2.0;
  *y = 1.0 / k / rho / q * sqrt ( 1.0 - sqr ( q / k ) );
  if ( 1.0 - sqr ( q / k ) < 0.0 )
    printf ( "q=%E  k=%E\n", q, k );
}

/* needed to create a more or less uniform mesh spacing along i */
double spacingfunction_i2 ( double x ) {
  double y, A, B, C, D;
  A = 6.61229818085967073309E-04;
  B = 1.28288906000451397027E-01;
  C = 2.77647338125862397362E+00;
  D = -1.90476228725907792416E+00;
  y = A + B * x + C * x * x + D * x * x * x;

  return ( y );
}

/* needed to create a more or less uniform mesh spacing along i */
double spacingfunction_i ( double x ) {
  double y;
  y = spacingfunction_i2 ( x * 0.994 ) / spacingfunction_i2 ( 0.994 );
  return ( y );
}

int main (  ) {
  double T, L, q, k, a, rho, P, x, y, x1, y1, x2, y2, vx, vy;
  long i, j, i2, i2max;
  FILE *postfile;

  FILE *gridfile, *initfile;

  gridfile = fopen ( gridfilename, "w" );
  initfile = fopen ( initfilename, "w" );
  if ( MIRROR )
    i2max = ni * 2 - 1;
  else
    i2max = ni;
  fprintf ( gridfile, "Grid2D(\n  is=1; ie=%ld; js=1; je=%ld;\n  Size(is,js,ie,je);\n", ( long ) i2max,
            ( long ) nj );

  postfile = fopen ( postfilename, "w" );
  fprintf ( postfile, "VARIABLES= \"X\", \"Y\",\n" );
  fprintf ( postfile, "\"V[0]\",\n" );
  fprintf ( postfile, "\"V[1]\",\n" );
  fprintf ( postfile, "\"M[0]\",\n" );
  fprintf ( postfile, "\"M[1]\",\n" );
  fprintf ( postfile, "\"rho\",\n" );
  fprintf ( postfile, "\"P\",\n" );
  fprintf ( postfile, "\"T\",\n" );
  fprintf ( postfile, "\"a\",\n" );
  fprintf ( postfile, "ZONE I=%ld, J=%ld F=POINT\n", ( long ) i2max, ( long ) nj );

  printf ( "Writing ringleb data to..\nGridfile: %s\nInit file: %s\nTecplot data file: %s\n", gridfilename,
           initfilename, postfilename );
  for ( j = 0; j < nj; j++ ) {
    for ( i2 = 0; i2 < i2max; i2++ ) {
      if ( i2 > ni - 2 )
        i = ni - ( i2 - ni ) - 2;
      else
        i = i2;
      if ( i >= 0 ) {
        k = kmin + ( kmax - kmin ) * ( double ) j / ( double ) ( nj - 1 );
        q = q0 + ( k - q0 ) * spacingfunction_i ( ( double ) i / ( double ) ( ni - 1 ) );
        Findxyfromkq ( k, q, &x, &y );
        if ( q + dq <= k ) {
          Findxyfromkq ( k, q, &x1, &y1 );
          Findxyfromkq ( k, q + dq, &x2, &y2 );
        } else {
          if ( q - dq > k )
            printf ( "problem here q=%E  dq=%E  k=%E   q-dq=%E  spacingfunction_i=%E\n", q, dq, k, q - dq,
                     spacingfunction_i ( ( double ) i / ( double ) ( ni - 1 ) ) );
          Findxyfromkq ( k, q - dq, &x1, &y1 );
          Findxyfromkq ( k, q, &x2, &y2 );
        }
        a = sqrt ( 1.0 - ( gamma - 1 ) / 2.0 * sqr ( q ) );
        rho = pow ( a, 2.0 / ( gamma - 1.0 ) );
        P = 1.0 / gamma * pow ( a, 2.0 * gamma / ( gamma - 1.0 ) );
//      J=1.0/a+1.0/3.0/powint(a,3)+1.0/5.0/powint(a,5) -0.5*log((1.0+a)/(1.0-a));
        L = sqrt ( sqr ( x2 - x1 ) + sqr ( y2 - y1 ) );
        vx = ( x2 - x1 ) / L * q;
        vy = ( y2 - y1 ) / L * q;
        if ( i != i2 ) {
          y = -y;
          vx = -vx;
        }
        //if (j==0) printf("i=%ld   i2=%ld  y=%E  vx=%E\n",i,i2,y,vx);
        T = P / rho / Rgas;
        fprintf ( postfile, "%E %E %E %E %E %E %E %E %E %E\n", x, y, vx, vy, vx / a, vy / a, rho, P, T, a );
        fprintf ( gridfile, "  Point(%ld,%ld,%E,%E);\n", 1 + i2, 1 + j, x, y );
        fprintf ( initfile, "    Region(%ld,%ld,%ld,%ld,INIT_TYPE1,%E,%E,%E,%E);\n", 1 + i2, 1 + j, 1 + i2, 1 + j, vx,
                  vy, P, T );

      }
    }
  }

  fprintf ( gridfile, ");\n" );

  fclose ( postfile );
  fclose ( gridfile );
  fclose ( initfile );
  printf ( "[done]\n" );
  return ( EXIT_SUCCESS );
}
