#include "share.h"
#include <stdio.h>
#include <string.h>

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

int main ( int argc, char **argv ) {
  double Re, x, M1, M2, M, x1, x2, Re1, Re2;
  bool Re_flag, x_flag, M_flag, validOptions, TURB, COMP;
  long cnt, numsteps;

  /* take care of the command line options */
  validOptions = FALSE;
  Re_flag = FALSE;
  x_flag = FALSE;
  M_flag = FALSE;

  if ( chkarg ( argc, argv, "-x" ) != 0 ) {
    x_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-x" ) + 1], "%lg", &x1 );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-x" ) + 2, argc - 1 )], "%lg", &x2 ) != 1 )
      x2 = x1;
  }
  if ( chkarg ( argc, argv, "-Re" ) != 0 ) {
    Re_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-Re" ) + 1], "%lg", &Re1 );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-Re" ) + 2, argc - 1 )], "%lg", &Re2 ) != 1 )
      Re2 = Re1;
  }
  if ( chkarg ( argc, argv, "-M" ) != 0 ) {
    M_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-M" ) + 1], "%lg", &M1 );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-M" ) + 2, argc - 1 )], "%lg", &M2 ) != 1 )
      M2 = M1;
  }
  if ( chkarg ( argc, argv, "-t" ) != 0 )
    TURB = TRUE;
  else
    TURB = FALSE;
  if ( chkarg ( argc, argv, "-c" ) != 0 )
    COMP = TRUE;
  else
    COMP = FALSE;
  if ( chkarg ( argc, argv, "-s" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-s" ) + 1], "%ld", &numsteps );
  } else {
    numsteps = 1;
  }
  if ( Re_flag && x_flag ) {
    validOptions = TRUE;
    if ( COMP && !M_flag )
      validOptions = FALSE;
  }

  if ( !validOptions ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag  \tArg                              \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-Re   \t<Reynolds number per m [1/m]>    \tdouble   \tY\n"
              "-x    \t<distance from leading edge>     \tdouble   \tY\n"
              "-M    \t<Mach number>                    \tdouble   \tY if comp\n"
              "-c    \t<compressible? none>             \tnone     \tN\n"
              "-s    \t<number of steps? 10>            \tinteger  \tN\n"
              "-t    \t<turbulent? none>                \tnone     \tN\n\n" );

  } else {
    if ( TURB )
      printf ( "Turbulent-" );
    else
      printf ( "Laminar-" );
    if ( COMP )
      printf ( "Compressible \n" );
    else
      printf ( "Incompressible \n" );
    printf ( "\n" );
    printf
      ( "x[m]           Re[1/m]        M              Cf             delta[m]       VanDriest      y/yplus[m]\n" );
    for ( cnt = 0; cnt < numsteps; cnt++ ) {
      M = cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( M2 - M1 ) + M1;
      Re = cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( Re2 - Re1 ) + Re1;
      x = cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( x2 - x1 ) + x1;
      printf ( "%E   %E   %E   %E   %E   %E   %E\n",
               x,
               Re,
               M,
               SkinFrictionCoeff ( M, Re, x, TURB, COMP ),
               BdryLayerThickness ( M, Re, x, TURB, COMP ),
               VanDriestIICorrel ( M ),
               1.0 / ( Re * sqrt ( SkinFrictionCoeff ( M, Re, x, TURB, COMP ) / 2.0 ) ) );
    }

  }

  return ( 0 );
}
