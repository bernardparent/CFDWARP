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
#include <stdlib.h>

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
  double T[3], P[3], delta[3], M[3], phi[3], H[3], x[4], y[4], U[3];
  double L, gamma, Rgas, Pdynamic, Mshock, mdot, altitude, Tratio, Pratio, tmp;
  bool shockvalid;
  long cnt;

  T[2] = 900.0e0;
  Rgas = 287.06e0;
  Pdynamic = 67032.0e0;
  delta[0] = 0.0;
  phi[0] = 0.0;

  if ( chkarg ( argc, argv, "-d1" ) != 0 && chkarg ( argc, argv, "-d2" ) != 0
       && chkarg ( argc, argv, "-l" ) != 0 && chkarg ( argc, argv, "-m" ) != 0
       && chkarg ( argc, argv, "-g" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-d1" ) + 1], "%lg", &delta[1] );
    delta[1] = delta[1] * pi / 180.0e0;
    sscanf ( argv[chkarg ( argc, argv, "-d2" ) + 1], "%lg", &delta[2] );
    delta[2] = delta[2] * pi / 180.0e0;
    sscanf ( argv[chkarg ( argc, argv, "-l" ) + 1], "%lg", &L );
    sscanf ( argv[chkarg ( argc, argv, "-m" ) + 1], "%lg", &M[0] );
    sscanf ( argv[chkarg ( argc, argv, "-g" ) + 1], "%lg", &gamma );
  } else {
    fprintf ( stderr, "\nAllowed and Required Flags:\n\n"
              "Flag  \tArg                                \tArg Type \tRequired? \n"
              "---------------------------------------------------------------\n"
              "-m    \t<Mach number>                      \tdouble   \tY\n"
              "-d1   \t<first flow turning angle in deg>  \tdouble   \tY\n"
              "-d2   \t<first flow turning angle in deg>  \tdouble   \tY\n"
              "-g    \t<gamma>                            \tdouble   \tY\n"
              "-l    \t<length>                           \tdouble   \tY\n\n" );

    exit ( EXIT_FAILURE );
  }

  FindExtCompInletProps ( M[0], L, Rgas, T[2], Pdynamic, gamma,
                          &T[0], &T[1], &P[0], &P[1], &P[2], &Mshock, &M[1], &M[2], &H[0], &mdot, &altitude );
  for ( cnt = 0; cnt < 2; cnt++ ) {
    FindWedgeShockProps ( gamma, delta[cnt + 1], M[cnt], &M[cnt + 1], &Tratio, &Pratio,
                          &phi[cnt + 1], &shockvalid );
    if ( !shockvalid )
      fprintf ( stderr, "problem, shock not valid\n" );
    T[cnt + 1] = T[cnt] * Tratio;
    P[cnt + 1] = P[cnt] * Pratio;
    H[cnt + 1] = H[cnt] * Tratio / Pratio;
  }
  for ( cnt = 0; cnt < 3; cnt++ ) {
    U[cnt] = M[cnt] * sqrt ( gamma * Rgas * T[cnt] );
  }
  x[0] = 0.0e0;
  y[0] = 0.0e0;

  x[1] =
    ( tan ( phi[1] ) * L - tan ( phi[2] + delta[1] ) * L ) / ( -tan ( phi[2] + delta[1] ) +
                                                               tan ( delta[1] ) );
  y[1] = tan ( delta[1] ) * x[1];

  H[0] = tan ( phi[1] ) * L;
  for ( cnt = 0; cnt < 2; cnt++ ) {
    H[cnt + 1] = H[cnt] / ( P[cnt + 1] / T[cnt + 1] * U[cnt + 1] ) * ( P[cnt] / T[cnt] * U[cnt] );
  }

  tmp = ( H[2] / tan ( phi[2] - delta[2] ) );
  x[2] = cos ( delta[2] + delta[1] ) * tmp + x[1];
  y[2] = sin ( delta[2] + delta[1] ) * tmp + y[1];
  x[3] = L;
  y[3] = tan ( phi[1] ) * L;

  printf ( "\nPosition of the compression ramp kinks, with last point the position of the cowl\n" );
  for ( cnt = 0; cnt < 4; cnt++ )
    printf ( "x[%ld]=%E     y[%ld]=%E\n", cnt, x[cnt], cnt, y[cnt] );

  printf ( "\nPressure, temperature, Mach number and Flow Speed in each region, starting with ambiant\n" );
  for ( cnt = 0; cnt < 3; cnt++ )
    printf ( "P[%ld]=%E  T[%ld]=%E  M[%ld]=%E  U[%ld]=%E\n", cnt, P[cnt], cnt, T[cnt], cnt, M[cnt], cnt,
             U[cnt] );

  printf ( "\nKink turning angle (in degrees) and angle of the shock (in degrees)\n" );
  for ( cnt = 1; cnt < 3; cnt++ )
    printf ( "delta[%ld]=%E  phi[%ld]=%E  \n", cnt, delta[cnt] * 180.0 / pi, cnt, phi[cnt] * 180.0 / pi );

  printf ( "\n" );
  return ( EXIT_SUCCESS );
}
