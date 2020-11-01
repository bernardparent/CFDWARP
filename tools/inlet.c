// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2001 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
   of conditions and the following disclaimer in the documentation and/or other
   materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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

  if ( chkarg ( argc, argv, "-M" ) != 0 ) {
    m_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-M" ) + 1], "%lg", &M1 );
  }
  if ( chkarg ( argc, argv, "-L" ) != 0 ) {
    l_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-L" ) + 1], "%lg", &L );
  }
  if ( chkarg ( argc, argv, "-Pdyn" ) != 0 ) {
    l_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-Pdyn" ) + 1], "%lg", &Pdynamic );
  }
  if ( chkarg ( argc, argv, "-T3" ) != 0 ) {
    l_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-T3" ) + 1], "%lg", &T3 );
  }
  if ( m_flag && l_flag ) {
    validOptions = TRUE;
  }

  if ( !validOptions ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag  \tArg                              \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-M   \t<flight Mach number [m/s]>       \tdouble   \tY\n"
              "-Pdyn\t<flight dynamic pressure [Pa]>   \tdouble   \tY\n"
              "-T3  \t<temperature at inlet exit [K]>  \tdouble   \tY\n"
              "-L   \t<length of the inlet [m]>        \tdouble   \tY\n" "\n\n" );

  } else {
    gamma = 1.4e0;
    Rgas = 287.06e0;
    FindExtCompInletProps ( M1, L, Rgas,
                            T3, Pdynamic, gamma,
                            &T1, &T2, &P1, &P2, &P3, &Ms, &M2, &M3, &W1, &mdot, &altitude );

    U1 = M1 * sqrt ( gamma * Rgas * T1 );
    U2 = M2 * sqrt ( gamma * Rgas * T2 );
    U3 = M3 * sqrt ( gamma * Rgas * T3 );
    W2 = W1 * ( P1 / T1 * U1 ) / ( P2 / T2 * U2 );
    W3 = W2 * ( P2 / T2 * U2 ) / ( P3 / T3 * U3 );

    printf ( "Altitude [m]  : %E\n"
             "Ms            : %E\n"
             "mdot [kg/ms]  : %E\n"
             "W1 [m]        : %E\n"
             "W2 [m]        : %E\n"
             "W3 [m]        : %E\n"
             "M1            : %E\n"
             "M2            : %E\n"
             "M3            : %E\n"
             "U1 [m/s]      : %E\n"
             "U2 [m/s]      : %E\n"
             "U3 [m/s]      : %E\n"
             "P1 [Pa]       : %E\n"
             "P2 [Pa]       : %E\n" 
             "P3 [Pa]       : %E\n" 
             "T1 [K]        : %E\n" 
             "T2 [K]        : %E\n" 
             "T3 [K]        : %E\n" ,
             altitude, Ms, mdot, W1, W2, W3, M1, M2, M3, U1, U2, U3, P1, P2, P3, T1, T2, T3 );
  }

  return ( 0 );
}
