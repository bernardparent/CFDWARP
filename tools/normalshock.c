#include "share.h"
#include <stdio.h>
#include <assert.h>

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
  double gamma_1, gamma_2, gamma, M1, M1_1, M1_2, M2, T2overT1, P2overP1, U2overU1, Pstag2overPstag1;
  bool gamma_flag, mach_flag;
  bool validOptions;
  long cnt, numsteps;

  /* take care of the command line options */
  validOptions = FALSE;
  gamma_flag = FALSE;
  mach_flag = FALSE;
  if ( chkarg ( argc, argv, "-m" ) != 0 ) {
    mach_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-m" ) + 1], "%lg", &M1_1 );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-m" ) + 2, argc - 1 )], "%lg", &M1_2 ) != 1 )
      M1_2 = M1_1;
  }
  if ( chkarg ( argc, argv, "-g" ) != 0 ) {
    gamma_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-g" ) + 1], "%lg", &gamma_1 );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-g" ) + 2, argc - 1 )], "%lg", &gamma_2 ) != 1 )
      gamma_2 = gamma_1;
  }
  if ( chkarg ( argc, argv, "-s" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-s" ) + 1], "%ld", &numsteps );
  } else {
    numsteps = 1;
  }
  if ( mach_flag && gamma_flag ) {
    validOptions = TRUE;
  }
  if ( !validOptions ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag  \tArg                 \tArg Type \tRequired? \n"
              "---------------------------------------------------------------\n"
              "-m    \t<Mach number>       \tdouble   \tY\n"
              "-g    \t<gamma>             \tdouble   \tY\n"
              "-s    \t<number of steps>   \tint      \tN\n" "\n" );
  } else {
    printf ( "M1            "
             "gamma         "
             "M2            " "T2 over T1    " "P2 over P1    " "U2 over U1    " "Pstag2 over Pstag1\n" );
    for ( cnt = 0; cnt < numsteps; cnt++ ) {
      M1 = cnt / ( double ) ( max ( numsteps - 1, 1 ) ) * ( M1_2 - M1_1 ) + M1_1;
      assert ( M1 >= 1.0e0 );
      gamma = cnt / ( double ) ( max ( numsteps - 1, 1 ) ) * ( gamma_2 - gamma_1 ) + gamma_1;
      FindNormalShockProps ( M1, gamma, &M2, &P2overP1, &T2overT1, &U2overU1, &Pstag2overPstag1 );
      printf ( "%E  %E  %E  %E  %E  %E  %E\n",
               M1, gamma, M2, T2overT1, P2overP1, U2overU1, Pstag2overPstag1 );
    }
  }
  return ( 0 );
}
