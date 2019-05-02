// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2001 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
  double delta, delta_1, delta_2, gamma, gamma_1, gamma_2, M1, M1_1, M1_2, M2, T2overT1, P2overP1, phi;
  bool delta_flag, gamma_flag, mach_flag;
  bool validOptions, shockvalid;
  long cnt, numsteps;

  /* take care of the command line options */
  validOptions = FALSE;
  delta_flag = FALSE;
  gamma_flag = FALSE;
  mach_flag = FALSE;
  if ( chkarg ( argc, argv, "-d" ) != 0 ) {
    delta_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-d" ) + 1], "%lg", &delta_1 );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-d" ) + 2, argc - 1 )], "%lg", &delta_2 ) != 1 )
      delta_2 = delta_1;
    delta_1 = delta_1 * pi / 180.0e0;
    delta_2 = delta_2 * pi / 180.0e0;
  }
  if ( chkarg ( argc, argv, "-M" ) != 0 ) {
    mach_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-M" ) + 1], "%lg", &M1_1 );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-M" ) + 2, argc - 1 )], "%lg", &M1_2 ) != 1 )
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
  if ( delta_flag && mach_flag && gamma_flag ) {
    validOptions = TRUE;
  }

  if ( !validOptions ) {
    fprintf ( stderr, "\nAllowed and Required Flags:\n\n"
              "Flag  \tArg                        \tArg Type \tRequired? \n"
              "---------------------------------------------------------------\n"
              "-M    \t<Mach number>              \tdouble   \tY\n"
              "-d    \t<flow turning angle in deg>\tdouble   \tY\n"
              "-s    \t<number of steps>          \tint      \tN\n"
              "-g    \t<gamma>                    \tdouble   \tY\n\n" );

  } else {
    printf ( "M1           "
             "gamma        "
             "delta (deg)  "
             "M2           "
             "T2/T1        " 
             "P2/P1        " 
             "q2/q1        " 
             "phi (deg)\n" );
    for ( cnt = 0; cnt < numsteps; cnt++ ) {
      M1 = cnt / ( double ) ( max ( numsteps - 1, 1 ) ) * ( M1_2 - M1_1 ) + M1_1;
      assert ( M1 >= 1.0e0 );
      gamma = cnt / ( double ) ( max ( numsteps - 1, 1 ) ) * ( gamma_2 - gamma_1 ) + gamma_1;
      delta = cnt / ( double ) ( max ( numsteps - 1, 1 ) ) * ( delta_2 - delta_1 ) + delta_1;
      FindWedgeShockProps ( gamma, delta, M1, &M2, &T2overT1, &P2overP1, &phi, &shockvalid );
      if ( shockvalid ) {
        printf ( "%E %E %E %E %E %E %E %E\n",
                 M1, gamma, delta/pi*180.0, M2, T2overT1,
                 P2overP1, M2 / M1 * sqrt ( T2overT1 ), 180.0e0 * phi / pi );
      } else {
        printf ( "Couldn't find root...\n" );
      }
    }
  }

  return ( 0 );
}
