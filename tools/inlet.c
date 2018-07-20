#include <exm.h>
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
  bool validOptions, l_flag, m_flag;
  double U1, U2, U3, T1, T2, T3, gamma, Rgas, Pdynamic, L, P1, P2, P3, Ms, M1, M2, M3, W1,
    mdot, altitude, W2, W3;

  /* take care of the command line options */
  validOptions = FALSE;
  m_flag = FALSE;
  l_flag = FALSE;

  if ( chkarg ( argc, argv, "-m" ) != 0 ) {
    m_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-m" ) + 1], "%lg", &M1 );
  }
  if ( chkarg ( argc, argv, "-l" ) != 0 ) {
    l_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-l" ) + 1], "%lg", &L );
  }
  if ( m_flag && l_flag ) {
    validOptions = TRUE;
  }

  if ( !validOptions ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag  \tArg                              \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-m   \t<flight Mach number [m/s]>       \tdouble   \tY\n"
              "-l   \t<length of the inlet [m]>        \tdouble   \tY\n" "\n\n" );

  } else {
    T3 = 900.0e0;
    gamma = 1.4e0;
    Rgas = 287.06e0;
    Pdynamic = 67032.0e0;
    FindExtCompInletProps ( M1, L, Rgas,
                            T3, Pdynamic, gamma,
                            &T1, &T2, &P1, &P2, &P3, &Ms, &M2, &M3, &W1, &mdot, &altitude );

    U1 = M1 * sqrt ( gamma * Rgas * T1 );
    U2 = M2 * sqrt ( gamma * Rgas * T2 );
    U3 = M3 * sqrt ( gamma * Rgas * T3 );
    W2 = W1 * ( P1 / T1 * U1 ) / ( P2 / T2 * U2 );
    W3 = W2 * ( P2 / T2 * U2 ) / ( P3 / T3 * U3 );

    printf ( "Altitude [m]  "
             "Ms            "
             "mdot [kg/ms]  "
             "W1 [m]        "
             "W2 [m]        "
             "W3 [m]        "
             "M1            "
             "M2            "
             "M3            "
             "P1 [Pa]       "
             "P2 [Pa]       " "P3 [Pa]       " "T1 [K]        " "T2 [K]        " "T3 [K]        \n" );
    printf ( "%E  %E  %E  %E  %E  %E  %E  %E  %E  %E  %E  %E  %E  %E  %E\n",
             altitude, Ms, mdot, W1, W2, W3, M1, M2, M3, P1, P2, P3, T1, T2, T3 );
  }

  return ( 0 );
}
