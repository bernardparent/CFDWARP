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

/* Program to determine the properties of hydrogen
 * and air for an injector tube flow located
 * in ZONE 2 of an external compression shcramjet,
 * flying at Mach MFlight. Pressure of hydrogen is set
 * equal to the pressure of the air in zone 2.
 */

#include <exm.h>
#include "share.h"
#include <stdio.h>

#define RgasH2 4124.0e0
#define RgasAir 287.0e0
#define gamma   1.4e0           /* same for air and hydrogen */

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

typedef struct {
  double PH2, TstagH2, mdotH2, AcsH2;
} arg_t;

static double _errfunct ( void *arg, double MH2 ) {
  double TstagH2, mdotH2, PH2, AcsH2;
  double TH2, UH2, tmp;

  TstagH2 = ( ( arg_t * ) arg )->TstagH2;
  mdotH2 = ( ( arg_t * ) arg )->mdotH2;
  AcsH2 = ( ( arg_t * ) arg )->AcsH2;
  PH2 = ( ( arg_t * ) arg )->PH2;

  TH2 = TstagH2 / ( 1.0e0 + ( gamma - 1.0e0 ) / 2.0e0 * MH2 * MH2 );
  UH2 = MH2 * sqrt ( gamma * RgasH2 * TH2 );
  tmp = UH2 * AcsH2 * PH2 / RgasH2 / TH2 - mdotH2;
  return ( tmp );
}

int main ( int argc, char **argv ) {
  bool validOptions, l_flag, m_flag, a_flag, t_flag, phi_flag;
  double U1, U2, T1, T2, T3, Rgas, Pdynamic, L, P1, P2, P3, Ms, M1, M2, M3, W1,
    mdot, altitude, W2, AcsH2_start, AcsH2_end;
  double AcsH2, TstagH2, MH2, mdot_air, mdot_H2, TH2, rhoH2, phi,
    UH2, TstagH2_start, TstagH2_end, L_start, L_end, phi_start, phi_end, M1_start, M1_end;
  long IFLAG, numsteps, cnt;
  arg_t arg;
  double MCair, MCH2, UC;
  bool noroot;

  /* take care of the command line options */
  validOptions = FALSE;
  phi_flag = FALSE;
  m_flag = FALSE;
  l_flag = FALSE;
  t_flag = FALSE;
  a_flag = FALSE;

  if ( chkarg ( argc, argv, "-m" ) != 0 ) {
    m_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-m" ) + 1], "%lg", &M1_start );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-m" ) + 2, argc - 1 )], "%lg", &M1_end ) != 1 )
      M1_end = M1_start;
  }
  if ( chkarg ( argc, argv, "-l" ) != 0 ) {
    l_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-l" ) + 1], "%lg", &L_start );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-l" ) + 2, argc - 1 )], "%lg", &L_end ) != 1 )
      L_end = L_start;
  }
  if ( chkarg ( argc, argv, "-t" ) != 0 ) {
    t_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-t" ) + 1], "%lg", &TstagH2_start );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-t" ) + 2, argc - 1 )], "%lg", &TstagH2_end ) != 1 )
      TstagH2_end = TstagH2_start;
  }
  if ( chkarg ( argc, argv, "-p" ) != 0 ) {
    phi_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-p" ) + 1], "%lg", &phi_start );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-p" ) + 2, argc - 1 )], "%lg", &phi_end ) != 1 )
      phi_end = phi_start;
  }
  if ( chkarg ( argc, argv, "-a" ) != 0 ) {
    a_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-a" ) + 1], "%lg", &AcsH2_start );
    if ( sscanf ( argv[min ( chkarg ( argc, argv, "-a" ) + 2, argc - 1 )], "%lg", &AcsH2_end ) != 1 )
      AcsH2_end = AcsH2_start;
  }
  if ( chkarg ( argc, argv, "-s" ) != 0 ) {
    sscanf ( argv[chkarg ( argc, argv, "-s" ) + 1], "%ld", &numsteps );
  } else {
    numsteps = 1;
  }
  if ( m_flag && l_flag && t_flag && a_flag && phi_flag ) {
    validOptions = TRUE;
  }

  if ( !validOptions ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag  \tArg                                        \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-m   \t<flight Mach number [m/s]>                  \tdouble   \tY\n"
              "-p   \t<global equivalence ratio>                  \tdouble   \tY\n"
              "-a   \t<fuel cross-section area per unit depth [m] \tdouble   \tY\n"
              "-t   \t<stag T of hydrogen jet [K]>                \tdouble   \tY\n"
              "-l   \t<length of the shcramjet inlet [m]>         \tdouble   \tY\n"
              "-s   \t<number of steps>                           \tlong     \tY\n" "\n\n" );

  } else {
    printf ( " Mfli "
             "   phi "
             "     H "
             "     L "
             "     P " "    Mc " " Uair " " Tair " "  Mair " "  UH2 " "  TH2 " "   MH2 " " ToH2 " "\n" );
    printf ( "      "
             "       "
             "     m "
             "     m "
             "    Pa " "       " "  m/s " "    K " "       " "  m/s  " "   K  " "      " "    K " "\n" );

    for ( cnt = 0; cnt < numsteps; cnt++ ) {
      AcsH2 = cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( AcsH2_end - AcsH2_start ) + AcsH2_start;
      TstagH2 =
        cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( TstagH2_end - TstagH2_start ) + TstagH2_start;
      phi = cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( phi_end - phi_start ) + phi_start;
      L = cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( L_end - L_start ) + L_start;
      M1 = cnt / ( double ) ( max ( 1, numsteps - 1 ) ) * ( M1_end - M1_start ) + M1_start;
      T3 = 900.0e0;
      Rgas = RgasAir;
      Pdynamic = 67032.0e0;
      FindExtCompInletProps ( M1, L, Rgas,
                              T3, Pdynamic, gamma,
                              &T1, &T2, &P1, &P2, &P3, &Ms, &M2, &M3, &W1, &mdot, &altitude );

      U1 = M1 * sqrt ( gamma * Rgas * T1 );
      U2 = M2 * sqrt ( gamma * Rgas * T2 );
      W2 = W1 * ( P1 / T1 * U1 ) / ( P2 / T2 * U2 );

      mdot_air = mdot;
      mdot_H2 = phi * 0.235e0 * mdot_air * 2.0e0 * 2.01588e0 / 31.9988e0;

      arg.PH2 = P2;
      arg.TstagH2 = TstagH2;
      arg.mdotH2 = mdot_H2;
      arg.AcsH2 = AcsH2;
      MH2 = EXM_find_root_zero_in ( &( _errfunct ), &arg, 1.1e0, 20.0e0, 1.0e-12, 1.0e-12, &IFLAG );
      if ( IFLAG == 4 ) {
        printf ( "no root\n" );
      } else {
        TH2 = TstagH2 / ( 1.0e0 + ( gamma - 1.0e0 ) / 2.0e0 * MH2 * MH2 );
        rhoH2 = P2 / RgasH2 / TH2;
        UH2 = MH2 * sqrt ( gamma * RgasH2 * TH2 );
        mdot_H2 = rhoH2 * UH2 * AcsH2;

        FindConvMachNumber ( gamma, RgasAir, P2, T2, U2,
                             gamma, RgasH2, P2, TH2, UH2, &MCair, &MCH2, &UC, &noroot );
        printf ( "%5.1f %6.3f %6.3f %6.3f %6.0f %6.2f %5.0f %5.0f %6.2f %5.0f %5.0f %6.2f %5.0f\n",
                 M1, phi, W2, L, P2, MCH2, U2, T2, M2, UH2, TH2, MH2, TstagH2 );
      }
    }
  }
  return ( 0 );
}
