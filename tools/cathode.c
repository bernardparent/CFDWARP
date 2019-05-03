// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2013 Bernard Parent

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
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#define echarge 1.60217646E-19  /* elementary charge in Coulomb */
#define kB 1.3806503e-23        /* Boltzman constant */
#define epsilon0 8.854188E-12   /* permittivity of free space */
#define ebar 2.7182818284

#define PLASMA_AIR 0
#define PLASMA_N2 1
#define PLASMA_N2lowE 2

typedef struct {
  double gamma, Nn, Ti, j;
  long plasmatype;
} arg_t;

/* the following is taken from Raizer, Gas Discharge Physics
       see p. 56 for A and B
*/
double _A ( long plasmatype ) {
  double A;
  A = 0.0;
  if ( plasmatype == PLASMA_AIR ) {
    A = 15.0;                   /* 1 / cm Torr */
  }
  if ( plasmatype == PLASMA_N2 ) {
    A = 12.0;                   /* 1 / cm Torr */
  }
  if ( plasmatype == PLASMA_N2lowE ) {
    A = 8.8;                    /* 1 / cm Torr */
  }
  return ( A );
}

/* the following is taken from Raizer, Gas Discharge Physics
       see p. 56 for A and B
*/
double _B ( long plasmatype ) {
  double B;
  B = 0.0;
  if ( plasmatype == PLASMA_AIR ) {
    B = 365.0;                  /* Vstar / cm Torr */
  }
  if ( plasmatype == PLASMA_N2 ) {
    B = 342.0;                  /* Vstar / cm Torr */
  }
  if ( plasmatype == PLASMA_N2lowE ) {
    B = 275.0;                  /* Vstar / cm Torr */
  }
  return ( B );
}

/* find the ion mobility at the cathode in units of m2/Vs */
double _mui_min ( long plasmatype, double Nn, double Ti, double Estar ) {
  double mui;
  mui = 0.0;
  Estar = max ( 1e-99, fabs ( Estar ) );
  if ( plasmatype == PLASMA_AIR ) {
    mui = 1.0 / Nn * min ( 8.32e22 / sqrt ( Ti ), 2.13e12 / sqrt ( Estar ) );
  }
  if ( plasmatype == PLASMA_N2 || plasmatype == PLASMA_N2lowE ) {
    mui = 1.0 / Nn * min ( 7.45e22 / sqrt ( Ti ), 1.76e12 / sqrt ( Estar ) );
  }
  return ( mui );
}

/* find the ion mobility at the edge of the sheath (near quasi-neutral region) in units of m2/Vs */
double _mui_max ( long plasmatype, double Nn, double Ti, double Estar ) {
  double mui;
  mui = _mui_min ( plasmatype, Nn, Ti, 0.0 );
  return ( mui );
}

/* find the ion mobility in units of m2/Vs */
double _mui ( long plasmatype, double Nn, double Ti, double Estar ) {
  double mui, mui1, mui2;
  mui = 0.0;
  mui1 = _mui_min ( plasmatype, Nn, Ti, Estar );
  mui2 = _mui_max ( plasmatype, Nn, Ti, Estar );
  /* use a geometric average between the mobility at the cathode and the mobility at the sheath edge */
  mui = 2.0 / ( 1.0 / mui1 + 1.0 / mui2 );
  return ( mui );
}

/* the following is taken from Raizer, Gas Discharge Physics
    see p. 180 for Vtilde, Etilde, dtilde, jtilde, and jmin
    see p. 134 for Pdmin, EoverPmin, Vmin
    see p. 56 for A and B

   valid range: 31.1e-20<E/N<249e-20 Vm2   

*/

double _jres ( void *arg, double dtilde ) {
  double jres, Pdmin, Vmin, P, jtilde, mui, jmin;
  double gamma, Nn, Ti, j, Estar, Etilde;
  long plasmatype;
  gamma = ( ( arg_t * ) arg )->gamma;
  Nn = ( ( arg_t * ) arg )->Nn;
  Ti = ( ( arg_t * ) arg )->Ti;
  j = ( ( arg_t * ) arg )->j;
  plasmatype = ( ( arg_t * ) arg )->plasmatype;
  P = Nn * Ti * kB / 133.0;     /* Torr */
  Pdmin = ebar / _A ( plasmatype ) * ( 1.0 / gamma + 1.0 );     /* Torr cm */
  Vmin = ebar * _B ( plasmatype ) / _A ( plasmatype ) * log ( 1.0 / gamma + 1.0 );      /* Vstar */

  jtilde = 1.0 / dtilde / sqr ( 1.0 + log ( dtilde ) );

  Etilde = 1.0 / ( 1.0 + log ( dtilde ) );
  Estar = Etilde * _B ( plasmatype ) * P * 100.0 / Nn;

  mui = _mui ( plasmatype, Nn, Ti, Estar ) * 1e4;       /* cm2 / Vstar s */
  jmin = sqr ( P ) * ( 1.0 + gamma ) / 9e11 * mui * P * sqr ( Vmin ) / 4.0 / pi / Pdmin / Pdmin / Pdmin;

  jres = j - jtilde * jmin * 1e4;
//  printf("dtilde=%E gamma=%E Nn=%E Ti=%E  jres=%E\n",dtilde,gamma,Nn,Ti,jres);

  return ( jres );
}

/* find the cathode sheath thickness in m */
double _d ( long plasmatype, double Ti, double Nn, double gamma, double j ) {
  double d, dtilde, relerr, abserr, dtildemin, dtildemax, Pdmin, P;
  long IFLAG;
  arg_t arg;

/*
  d=0.01; 
  dd=1e-8; 
  do {
    res= (j-_j(Ti,Nn,gamma,d)) ;
    dddj=dd/(_j(Ti,Nn,gamma,d+dd)-_j(Ti,Nn,gamma,d));
    d=max(0.01,d+dd);
    dd=dddj*res;
    printf("d=%E  j=%E  res=%E\n",d,j,res);
  } while(fabs(res/j)>1e-10);
*/

  arg.gamma = gamma;
  arg.Nn = Nn;
  arg.Ti = Ti;
  arg.j = j;
  arg.plasmatype = plasmatype;
  dtildemin = 1.0 / ebar;
  dtildemax = 2000.0;
  relerr = 1e-10;
  abserr = 1e-100;
  dtilde = EXM_find_root_zero_in ( &_jres, &arg, dtildemin, dtildemax, relerr, abserr, &IFLAG );
  if ( IFLAG != 1 )
    printf ( "a root for dtilde has not been found (IFLAG=%ld)\n", IFLAG );
  P = Nn * Ti * kB / 133.0;     /* Torr */
  Pdmin = ebar / _A ( plasmatype ) * ( 1.0 / gamma + 1.0 );     /* Torr cm */

  d = dtilde * Pdmin / P * 0.01;

  return ( d );
}

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
  double j, jtilde, gamma, Pdmin, EoverPmin, Vmin, dtilde, d, P, Vtilde, Etilde, Ti, Nn;
  bool Ti_flag, gamma_flag, P_flag, j_flag;
  long plasmatype;
  char *strtmp;

  Ti = 300;                     /* K */
  P = 101300;                   /* Pa */
  j = 0.5;                      /* A/m2 */
  gamma = 0.1;
  strtmp = ( char * ) malloc ( 0 );
  /* take care of the command line options */
  Ti_flag = FALSE;
  gamma_flag = FALSE;
  P_flag = FALSE;
  j_flag = FALSE;

  if ( chkarg ( argc, argv, "-Ti" ) != 0 ) {
    Ti_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-Ti" ) + 1], "%lg", &Ti );
  }

  if ( chkarg ( argc, argv, "-P" ) != 0 ) {
    P_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-P" ) + 1], "%lg", &P );
  }

  if ( chkarg ( argc, argv, "-j" ) != 0 ) {
    j_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-j" ) + 1], "%lg", &j );
  }

  if ( chkarg ( argc, argv, "-g" ) != 0 ) {
    gamma_flag = TRUE;
    sscanf ( argv[chkarg ( argc, argv, "-g" ) + 1], "%lg", &gamma );
  }

  plasmatype = -1;
  if ( chkarg ( argc, argv, "-type" ) != 0 ) {

    strtmp = ( char * ) realloc ( strtmp,
                                  ( 2 +
                                    ( long ) strlen ( argv[chkarg ( argc, argv, "-type" ) + 1] ) ) *
                                  sizeof ( char ) );
    strcpy ( strtmp, argv[chkarg ( argc, argv, "-type" ) + 1] );
    if ( strcmp ( strtmp, "AIR" ) == 0 )
      plasmatype = PLASMA_AIR;
    if ( strcmp ( strtmp, "N2" ) == 0 )
      plasmatype = PLASMA_N2;
    if ( strcmp ( strtmp, "N2lowE" ) == 0 )
      plasmatype = PLASMA_N2lowE;
    if ( plasmatype == -1 ) {
      printf ( "wrong plasma type :%s:\n", strtmp );
      exit ( EXIT_FAILURE );
    }
  }

  if ( !( Ti_flag && P_flag && j_flag && gamma_flag ) ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag   \tArg                                    \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-g     \t<electron emission coefficient gamma>  \tdouble   \tY\n"
              "-Ti    \t<ion temperature [K]>                  \tdouble   \tY\n"
              "-P     \t<plasma pressure [Pa]>                 \tdouble   \tY\n"
              "-j     \t<current density [A/m2]>               \tdouble   \tY\n"
              "-type  \t<plasma type, AIR, N2, N2lowE>         \tstring   \tY\n" "\n" );

  } else {
    Nn = P / kB / Ti;
    d = _d ( plasmatype, Ti, Nn, gamma, j ) * 100.0;    /* cm */

    P = P / 133.0;              /* Torr */

    Pdmin = ebar / _A ( plasmatype ) * ( 1.0 / gamma + 1.0 );   /* Torr cm */
    EoverPmin = _B ( plasmatype );      /* Vstar / cm Torr */
    Vmin = ebar * _B ( plasmatype ) / _A ( plasmatype ) * log ( 1.0 / gamma + 1.0 );    /* Vstar */

    dtilde = P * d / Pdmin;
    Vtilde = dtilde / ( 1.0 + log ( dtilde ) );
    Etilde = 1.0 / ( 1.0 + log ( dtilde ) );
    jtilde = 1.0 / dtilde / sqr ( 1.0 + log ( dtilde ) );
    printf ( "\n" );
    printf ( "Type of Plasma: " );
    if ( plasmatype == PLASMA_N2 )
      printf ( "N2 (100<EoverP<600 Vstar/cmTorr)\n" );
    if ( plasmatype == PLASMA_N2lowE )
      printf ( "N2 (27<EoverP<200 Vstar/cmTorr)\n" );
    if ( plasmatype == PLASMA_AIR )
      printf ( "AIR (100<EoverP<800 Vstar/cmTorr)\n" );
    printf ( "Type of Discharge: " );
    if ( jtilde > 1.0 )
      printf ( "High Current Discharge\n" );
    else
      printf ( "Low Current Discharge\n" );
    printf ( "\n" );
    printf ( "mui (at the cathode):                  %E m2/Vs\n",
             _mui_min ( plasmatype, Nn, Ti, Etilde * EoverPmin * P * 100.0 / Nn ) );
    printf ( "mui (at the sheath edge):              %E m2/Vs\n",
             _mui_max ( plasmatype, Nn, Ti, Etilde * EoverPmin * P * 100.0 / Nn ) );
    printf ( "mui (geometric average):               %E m2/Vs\n",
             _mui ( plasmatype, Nn, Ti, Etilde * EoverPmin * P * 100.0 / Nn ) );
    printf ( "d:                                     %E m \n", d / 100.0 );
    printf ( "Pressure:                              %E torr\n", P );
    printf ( "Voltage drop in cathode sheath:        %E Vstar\n", Vtilde * Vmin );
    printf ( "Current in cathode sheath:             %E A/m2\n", j );
    printf ( "Electric field at the cathode:         %E Vstar/m\n", Etilde * EoverPmin * P * 100.0 );
    printf ( "Reduced Electric field at the cathode: %E Vm2\n", Etilde * EoverPmin * P * 100.0 / Nn );
    printf ( "E over P at the cathode:               %E Vstar/cmTorr\n", EoverPmin );
    printf ( "jtilde:                                %E \n", jtilde );

  }
  free ( strtmp );
  return ( EXIT_SUCCESS );
}
