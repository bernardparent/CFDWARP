#include "share.h"
#include <stdio.h>

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
  double MC1, MC2, P1, P2, T1, T2, R1, R2, gamma1, gamma2, U1, U2, UC, U1b, U2b, U1a, U2a, fact;
  bool U1_flag, gamma1_flag, P1_flag, T1_flag, R1_flag;
  bool U2_flag;
  bool validOptions, noroot;
  long numsteps, step;
  double dshear, betac, betad;

  /* take care of the command line options */
  validOptions = FALSE;
  U1_flag = FALSE;
  U2_flag = FALSE;
  gamma1_flag = FALSE;
  P1_flag = FALSE;
  T1_flag = FALSE;
  R1_flag = FALSE;

  if ( chkarg ( argc, argv, "-U1" ) != 0 ) {
    U1_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-U1" ) + 1], "%lg", &U1a );
    if ( !( sscanf ( argv[chkarg ( argc, argv, "-U1" ) + 2], "%lg", &U1b ) == 1 ) )
      U1b = U1a;
  }
  if ( chkarg ( argc, argv, "-U2" ) != 0 ) {
    U2_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-U2" ) + 1], "%lg", &U2a );
    if ( !( sscanf ( argv[chkarg ( argc, argv, "-U2" ) + 2], "%lg", &U2b ) == 1 ) )
      U2b = U2a;
  }
  if ( chkarg ( argc, argv, "-g1" ) != 0 ) {
    gamma1_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-g1" ) + 1], "%lg", &gamma1 );
    gamma2 = gamma1;
  }
  if ( chkarg ( argc, argv, "-g2" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-g2" ) + 1], "%lg", &gamma2 );
  }
  if ( chkarg ( argc, argv, "-P1" ) != 0 ) {
    P1_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-P1" ) + 1], "%lg", &P1 );
    P2 = P1;
  }
  if ( chkarg ( argc, argv, "-P2" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-P2" ) + 1], "%lg", &P2 );
  }
  if ( chkarg ( argc, argv, "-T1" ) != 0 ) {
    T1_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-T1" ) + 1], "%lg", &T1 );
    T2 = T1;
  }
  if ( chkarg ( argc, argv, "-T2" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-T2" ) + 1], "%lg", &T2 );
  }
  if ( chkarg ( argc, argv, "-R1" ) != 0 ) {
    R1_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-R1" ) + 1], "%lg", &R1 );
    R2 = R1;
  }
  if ( chkarg ( argc, argv, "-R2" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-R2" ) + 1], "%lg", &R2 );
  }
  if ( chkarg ( argc, argv, "-s" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-s" ) + 1], "%ld", &numsteps );
  } else {
    numsteps = 1;
  }

  if ( U1_flag && U2_flag && R1_flag && P1_flag && T1_flag && gamma1_flag ) {
    validOptions = TRUE;
  }

  if ( !validOptions ) {
    fprintf ( stderr, "\nAllowed and Required Flags:\n\n"
              "Flag  \tArg                              \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-U1   \t<speed of jet 1 [m/s]>           \tdouble   \tY\n"
              "-U2   \t<speed of jet 2 [m/s]>           \tdouble   \tY\n"
              "-P1   \t<pressure of jet 1 [pa]>         \tdouble   \tY\n"
              "-P2   \t<pressure of jet 2 [pa]>         \tdouble   \tN\n"
              "-T1   \t<temperature of jet 1 [K]>       \tdouble   \tY\n"
              "-T2   \t<temperature of jet 2 [K]>       \tdouble   \tN\n"
              "-R1   \t<gas constant of jet 1 [J/kgK]>  \tdouble   \tY\n"
              "-R2   \t<gas constant of jet 2 [J/kgK]>  \tdouble   \tN\n"
              "-g1   \t<ratio of S-H of jet 1>          \tdouble   \tY\n"
              "-g1   \t<ratio of S-H of jet 1>          \tdouble   \tY\n"
              "-s    \t<number of steps>                \tdouble   \tN\n\n" );

  } else {
    numsteps++;
    printf
      ( "U1 [m/s]      U2 [m/s]      MC1           MC2           betac         betad         growth [m/m]\n" );
    for ( step = 0; step <= numsteps; step++ ) {
      fact = 1.0e0 - ( double ) ( step ) / ( double ) ( max ( 1.0e0, numsteps - 1.0e0 ) );
      U1 = fact * U1a + ( 1.0e0 - fact ) * U1b;
      U2 = fact * U2a + ( 1.0e0 - fact ) * U2b;
      FindConvMachNumber ( gamma1, R1, P1, T1, U1, gamma2, R2, P2, T2, U2, &MC1, &MC2, &UC, &noroot );
      if ( noroot ) {
        printf ( "no root could be found.. \n" );
      } else {

        FindShearLayerGrowth ( gamma1, R1, P1, T1, U1,
                               gamma2, R2, P2, T2, U2, UC, MC1, MC2, &dshear, &betac, &betad );
        printf ( "%E  %E  %E  %E  %E  %E  %E\n", U1, U2, MC1, MC2, betac, betad, dshear );

      }

    }
  }

  return ( 0 );
}
