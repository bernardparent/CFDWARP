#include "share.h"
#include <stdio.h>
#include <assert.h>
#include <exm.h>

#define FACT_PRECISION 10000

int chkarg ( int argc, char **argv, char *arg ) {
  int cnt, tmp;
  tmp = 0;
  for ( cnt = 1; cnt < argc; cnt++ ) {
    if ( strcmp ( argv[cnt], arg ) == 0 ) {
      tmp = cnt;
    }
  }
  return ( tmp );
}

double phi ( double M ) {
  double tmp;
  tmp = atan ( 1.0 / sqrt ( M * M - 1.0 ) );
  return ( tmp );
}

int main ( int argc, char **argv ) {
  double delta, r, xramp, yramp, L, H, alpha, T2overT1;
  long numsteps;
  double ddelta, gamma, dM, delta_end, delta_start, M, M1, x, y, xcowl, ycowl;
  long cnt;
  bool VALIDOPTIONS;

  /* take care of the command line options */
  VALIDOPTIONS = TRUE;
  if ( chkarg ( argc, argv, "-d" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-d" ) + 1], "%lg", &delta_end );
    delta_end = delta_end * pi / 180.0e0;
  } else
    VALIDOPTIONS = FALSE;
  if ( chkarg ( argc, argv, "-M" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-M" ) + 1], "%lg", &M );
  } else
    VALIDOPTIONS = FALSE;
  if ( chkarg ( argc, argv, "-g" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-g" ) + 1], "%lg", &gamma );
  } else
    VALIDOPTIONS = FALSE;
  if ( chkarg ( argc, argv, "-s" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-s" ) + 1], "%ld", &numsteps );
  } else
    VALIDOPTIONS = FALSE;
  if ( chkarg ( argc, argv, "-x" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-x" ) + 1], "%lg", &xcowl );
  } else
    VALIDOPTIONS = FALSE;
  if ( chkarg ( argc, argv, "-y" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-y" ) + 1], "%lg", &ycowl );
  } else
    VALIDOPTIONS = FALSE;
  if ( chkarg ( argc, argv, "-H" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-H" ) + 1], "%lg", &H );
  } else
    VALIDOPTIONS = FALSE;

  if ( !VALIDOPTIONS ) {
    fprintf ( stderr, "\nAllowed and Required Flags:\n\n"
              "Flag   \tArg                        \tArg Type \tRequired? \n"
              "---------------------------------------------------------------\n"
              "-x     \t<x position of cowl>        \tdouble   \tY\n"
              "-y     \t<y position of cowl>        \tdouble   \tY\n"
              "-H     \t<height of fan>             \tdouble   \tY\n"
              "-M     \t<Mach number>               \tdouble   \tY\n"
              "-d     \t<flow turning angle in deg> \tdouble   \tY\n"
              "-s     \t<number of steps>           \tint      \tY\n"
              "-g     \t<gamma>                     \tdouble   \tY\n\n" );

  } else {
    M1 = M;
    alpha = asin( 1.0 / M );
    L = H / tan( alpha );
    numsteps = numsteps * FACT_PRECISION + 1;
    delta_start = rad ( 0.0 );
    xramp = xcowl - L;
    yramp = ycowl + ( xramp - xcowl ) * tan ( delta_start + asin ( 1.0 / M ) );
    //printf("%E %E %E\n",xramp,xcowl,L);
    r = ( ycowl - yramp - ( xcowl - xramp ) * tan ( delta_start ) )
      / ( sin ( delta_start + phi ( M ) ) - cos ( delta_start + phi ( M ) ) * tan ( delta_start ) );
    delta = delta_start;
    ddelta = ( delta_end - delta_start ) / numsteps;
    printf ( "x (m)            y (m)\n");
    for ( cnt = 0; cnt < numsteps; cnt++ ) {
      dM = -ddelta * M * ( 1.0 + ( gamma - 1.0 ) / 2.0 * M * M ) / sqrt ( M * M - 1.0 );
      r = r * sin ( phi ( M ) - ddelta ) / sin ( pi - phi ( M + dM ) );
      M += dM;
      delta += ddelta;
      x = xcowl - r * cos ( delta + phi ( M ) );
      y = ycowl - r * sin ( delta + phi ( M ) );
      if ( cnt % FACT_PRECISION == 0 )
        printf ( "%E    %E\n", x, y );
    }

    T2overT1=(1.0 + (gamma - 1.0)/2.0 * M1*M1) / (1.0 + (gamma - 1.0)/2.0 * M*M);
    printf ( "\nM1    = %E", M1 );
    printf ( "\nM2    = %E", M );
    printf ( "\nP2/P1 = %E", pow( T2overT1 , gamma/(gamma - 1.0)) );
    printf ( "\nT2/T1 = %E\n", T2overT1 );
  }

  return ( 0 );
}
