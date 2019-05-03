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

#include "share.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

static void InitArray ( double *array, long num, ... ) {
  va_list ap;
  long cnt;

  va_start ( ap, num );
  for ( cnt = 0; cnt < num; cnt++ ) {
    array[cnt] = va_arg ( ap, double );
  }
  va_end ( ap );
}

/**************************************************************
 Turbulent Shear Layer growth predictions
 **************************************************************/

typedef struct {
  double gamma1g, R1g, P1g, T1g, U1g;
  double gamma2g, R2g, P2g, T2g, U2g;
} arg_t;

static double _errfunct ( void *arg, double UC ) {
  double tmp;
  tmp =
    ( ( arg_t * ) arg )->P1g * pow ( 1.0e0 +
                                     ( ( ( arg_t * ) arg )->gamma1g -
                                       1.0e0 ) / 2.0e0 * sqr ( ( UC - ( ( arg_t * ) arg )->U1g )
                                                               / sqrt ( ( ( arg_t * ) arg )->gamma1g *
                                                                        ( ( arg_t * ) arg )->R1g *
                                                                        ( ( arg_t * ) arg )->T1g ) ),
                                     ( ( arg_t * ) arg )->gamma1g / ( ( ( arg_t * ) arg )->gamma1g - 1.0e0 ) )
    - ( ( arg_t * ) arg )->P2g * pow ( 1.0e0 +
                                       ( ( ( arg_t * ) arg )->gamma2g -
                                         1.0e0 ) / 2.0e0 * sqr ( ( UC - ( ( arg_t * ) arg )->U2g )
                                                                 / sqrt ( ( ( arg_t * ) arg )->gamma2g *
                                                                          ( ( arg_t * ) arg )->R2g *
                                                                          ( ( arg_t * ) arg )->T2g ) ),
                                       ( ( arg_t * ) arg )->gamma2g / ( ( ( arg_t * ) arg )->gamma2g -
                                                                        1.0e0 ) );
  return ( tmp );
}

void FindConvMachNumber ( double gamma1, double R1, double P1, double T1, double U1,
                          double gamma2, double R2, double P2, double T2, double U2,
                          double *MC1, double *MC2, double *UC, bool * noroot ) {
  long IFLAG;
  arg_t arg;

  arg.gamma1g = gamma1;
  arg.R1g = R1;
  arg.P1g = P1;
  arg.T1g = T1;
  arg.U1g = U1;
  arg.gamma2g = gamma2;
  arg.R2g = R2;
  arg.P2g = P2;
  arg.T2g = T2;
  arg.U2g = U2;

  *UC = EXM_find_root_zero_in ( &( _errfunct ), &arg, min ( U1, U2 ), max ( U1, U2 ),
                                1.0e-12, 1.0e-12, &IFLAG );
  *MC1 = fabs ( *UC - U1 ) / sqrt ( gamma1 * R1 * T1 );
  *MC2 = fabs ( *UC - U2 ) / sqrt ( gamma2 * R2 * T2 );
  *noroot = FALSE;
  if ( IFLAG == 4 )
    *noroot = TRUE;
/*  printf("%E\n",_errfunct(*UC)); */
}

void FindShearLayerGrowth ( double gamma1, double R1, double P1, double T1, double U1,
                            double gamma2, double R2, double P2, double T2, double U2,
                            double UC, double MC1, double MC2, double *dshear,
                            double *betac, double *betad ) {
  double r, s, dsheardx1, dsheardx2;
  if ( U2 < U1 ) {
    r = U2 / U1;
    s = sqrt ( P2 / ( R2 * T2 ) / P1 * ( R1 * T1 ) );
  } else {
    r = U1 / U2;
    s = sqrt ( P1 / ( R1 * T1 ) / P2 * ( R2 * T2 ) );
  }
  *betad = 1.0e0 - ( 1.0e0 - s ) / ( 1.0e0 + s ) / ( 1.0e0 + 2.9e0 * ( 1.0e0 + r ) / ( 1.0e0 - r ) );
  *betac = 0.2e0 + 0.8e0 * exp ( -3.0e0 * sqr ( ( MC1 + MC2 ) / 2.0e0 ) );
  dsheardx1 = 0.17 * ( *betad ) * ( *betac ) * ( 1 - r ) * ( 1 + s ) / ( 1 + r * s );
  dsheardx2 = 0.17 * ( *betad ) * ( *betac ) * ( 1 - r ) * ( 1 + s ) / ( 1 + r * s );
  *dshear = ( dsheardx1 + dsheardx2 ) / 2.0e0;
}

/****************************************************************
 One Dimensional Gas Dynamic theory to find a normal shock
 ****************************************************************/

/* from shock angle phi (in radians) find the props after the shock */
void FindWedgeShockPropsFromPhi ( double gam,
                                  double phi, double M1,
                                  double *T2overT1, double *P2overP1,
                                  double *Mx2, double *My2, bool * shockvalid ) {
  double Ms1, Ms2, Mp1, Mp2;
  Ms1 = M1 * cos ( pi * 0.5e0 - phi );
  Mp1 = M1 * sin ( pi * 0.5e0 - phi );
  if ( Ms1 > 1.0e0 ) {
    *shockvalid = TRUE;
    *P2overP1 = ( 2.0e0 * gam / ( gam + 1.0e0 ) * Ms1 * Ms1 - ( gam - 1.0e0 ) / ( gam + 1.0e0 ) );
    *T2overT1 = ( *P2overP1 ) * ( ( gam - 1.0e0 ) / ( gam + 1.0e0 ) + 2.0e0 / ( gam + 1.0e0 ) / Ms1 / Ms1 );
    Ms2 =
      sqrt ( ( Ms1 * Ms1 + 2.0e0 / ( gam - 1.0e0 ) ) / ( 2.0e0 * gam / ( gam - 1.0e0 ) * Ms1 * Ms1 -
                                                         1.0e0 ) );
    Mp2 = Mp1 * sqrt ( 1.0e0 / ( *T2overT1 ) );
    *Mx2 = cos ( 0.5e0 * pi - phi ) * Ms2 + cos ( phi ) * Mp2;
    *My2 = -sin ( 0.5e0 * pi - phi ) * Ms2 + sin ( phi ) * Mp2;
  } else {
    *shockvalid = FALSE;
  }
}

/* from wedge angle del (in radians) find the props after the shock 
   from a phi limit of phi1<phi<phi2*/
void FindWedgeShockProps ( double gamma, double del,
                           double M1,
                           double *M2, double *T2overT1, double *P2overP1, double *phiout, bool * valid ) {
  bool shockvalid, shockvalidold;
  double delnew, err, errold, phi, dphi, TOL, phi1, phi2, Mx2, My2;
  long cnt, maxiter;

  phi1 = 0.0e0;
  phi2 = pi * 0.5e0;
  TOL = 1.0e-10;
  maxiter = 10000;
  dphi = ( phi2 - phi1 ) / 10000.0e0;
  phi = phi1 - 3.0e0 * dphi;
  shockvalid = FALSE;
  cnt = 0;
  err = 1.0e99;
  Mx2 = 0.0;                    // to avoid compiler warning
  My2 = 0.0;                    // to avoid compiler warning
  do {
    phi = phi + dphi;
    cnt++;
    shockvalidold = shockvalid;
    FindWedgeShockPropsFromPhi ( gamma, phi, M1, T2overT1, P2overP1, &Mx2, &My2, &shockvalid );
    errold = err;
    delnew = atan ( My2 / Mx2 );
    err = delnew - del;
    if ( !shockvalid )
      err = 1.0e99;
    if ( ( err * errold < 0.0e0 || err == 0.0e0 ) && shockvalid && shockvalidold ) {
      phi = phi + 2.0e0 * dphi;
      dphi = -dphi * 0.3e0;
      shockvalid = FALSE;
    }
  } while ( fabs ( err ) > TOL && cnt < maxiter );
  *valid = TRUE;
  *M2 = sqrt ( Mx2 * Mx2 + My2 * My2 );
  if ( cnt == maxiter )
    *valid = FALSE;
  *phiout = phi;
}

/*****************************************************************
 Boundary Layer thickness and wall shear stress
 *****************************************************************/

/* adiabatic wall only!! and M<6*/
double VanDriestIICorrel ( double M ) {
  double minimum, maximum, Mmax, tmp;
  minimum = 0.28e0;
  maximum = 1.0e0;
  Mmax = 6.0e0;
  tmp = minimum + ( maximum - minimum ) * exp ( -3.0e0 * sqr ( M / Mmax ) );
  return ( tmp );
}

double SkinFrictionCoeff ( double M, double Re, double x, bool TURB, bool COMP ) {
  double tmp;
  if ( COMP ) {
    if ( TURB ) {
      /* the first part (before the VanDriestII) is taken from
         a standard textbook; see Introduction to Fluid Mechanics
         by Fox/Macdonald for example. It can be used confidently
         for incompressible flows from Re=5.0E5 to Re=1.0E7. 
         the second part is the Van Driest II correlation applicable
         till Mach 6 approximately to transform the incompressible
         boundary layer Cf to a compressible one. */
      tmp = 0.0594e0 / pow ( Re * x, 0.2e0 ) * VanDriestIICorrel ( M );
    } else {
      /* ?????? */
      /* cold flat plate --> approximate <-- below Mach 20 */
      tmp = ( 0.730e0 - fabs ( M ) / 70.0e0 ) / pow ( Re * x, 0.5e0 );
    }
  } else {
    if ( TURB ) {
      tmp = 0.0594e0 / pow ( Re * x, 0.2e0 );
    } else {
      tmp = 0.730e0 / pow ( Re * x, 0.5e0 );
    }
  }
  return ( tmp );
}

double BdryLayerThickness ( double M, double Re, double x, bool TURB, bool COMP ) {
  double tmp;
  if ( COMP ) {
    if ( TURB ) {
      /* ?????? */
      tmp = 0.0e0;
    } else {
      /* ?????? */
      /* over a cold flat plate only!! --> approx <--- below Mach 20 */
      tmp = ( 5.48e0 + fabs ( M ) * 1.2e0 ) * x / sqrt ( Re * x );
    }
  } else {
    if ( TURB ) {
      tmp = 0.382e0 * x / pow ( Re * x, 0.2e0 );
    } else {
      tmp = 5.48e0 * x / sqrt ( Re * x );
    }
  }
  return ( tmp );
}

/*****************************************************************
 Normal Shock Properties
 *****************************************************************/

void FindNormalShockProps ( double M1, double gam, double *M2, double *P2overP1,
                            double *T2overT1, double *U2overU1, double *Pstag2overPstag1 ) {
  double one, two;
  one = 1.0e0;
  two = 2.0e0;
  *P2overP1 = two * gam / ( gam + one ) * sqr ( M1 ) - ( gam - one ) / ( gam + one );
  *T2overT1 = ( two * gam / ( gam + one ) * sqr ( M1 ) - ( gam - one ) / ( gam + one ) )
    * ( ( gam - one ) / ( gam + one ) + two / ( gam + one ) / sqr ( M1 ) );
  *U2overU1 = ( two + ( gam - one ) * sqr ( M1 ) ) / ( ( gam + one ) * sqr ( M1 ) );
  *Pstag2overPstag1 = pow ( pow ( ( gam - one ) / ( gam + one ) + two / ( gam + one ) / sqr ( M1 ), gam )
                            * ( two * gam / ( gam + one ) * sqr ( M1 ) - ( gam - one ) / ( gam + one ) )
                            , -one / ( gam - one ) );
  *M2 = sqrt ( ( sqr ( M1 ) + two / ( gam - one ) ) / ( two * gam / ( gam - one ) * sqr ( M1 ) - one ) );
}

/*****************************************************************
 US standard atmosphere
*****************************************************************/

#define nn 38

void FindTandHfromPstdatm ( double P, double *H, double *T ) {
  double Patm[nn], Tatm[nn], hatm[nn];
  long i, nn1, nn2;
  double dp, dp1, fact1, fact2;

  InitArray ( Patm, nn,
              1.0610e0, 1.0000e0, 0.9421e0, 0.8870e0, 0.8345e0, 0.7846e0,
              0.7372e0, 0.6920e0, 0.6492e0, 0.6085e0, 0.5700e0, 0.5334e0,
              0.4660e0, 0.4057e0, 0.3519e0, 0.3040e0, 0.2615e0, 0.2240e0,
              0.1915e0, 0.1636e0, 0.1399e0, 0.1195e0, 0.1022e0, 8.734e-2,
              7.466e-2, 6.383e-2, 5.457e-2, 3.995e-2, 2.933e-2, 2.160e-2,
              1.595e-2, 1.181e-2, 2.834e-3, 7.874e-4, 2.217e-4, 5.448e-5, 1.023e-5, 1.622e-6 );
  InitArray ( Tatm, nn,
              291.4e0, 288.2e0, 284.9e0, 281.7e0, 278.4e0, 275.2e0,
              271.9e0, 268.7e0, 265.4e0, 262.2e0, 258.9e0, 255.7e0,
              249.2e0, 242.7e0, 236.2e0, 229.7e0, 223.3e0, 216.8e0,
              216.7e0, 216.7e0, 216.7e0, 216.7e0, 216.7e0, 216.7e0,
              216.7e0, 216.7e0, 216.7e0, 218.6e0, 220.6e0, 222.5e0,
              224.5e0, 226.5e0, 250.4e0, 270.7e0, 255.8e0, 219.7e0, 180.7e0, 180.7e0 );

  for ( i = 0; i < nn; i++ )
    Patm[i] = Patm[i] * 1.01325e5;
  hatm[0] = -500.0e0;
  for ( i = 1; i < 12; i++ )
    hatm[i] = hatm[i - 1] + 500.0e0;
  for ( i = 12; i < 27; i++ )
    hatm[i] = hatm[i - 1] + 1000.0e0;
  for ( i = 27; i < 32; i++ )
    hatm[i] = hatm[i - 1] + 2000.0e0;
  for ( i = 32; i < 38; i++ )
    hatm[i] = hatm[i - 1] + 10000.0e0;

  nn1 = -1;
  nn2 = -1;
  for ( i = 0; i < nn - 1; i++ ) {
    if ( ( Patm[i] >= P ) && ( Patm[i + 1] < P ) ) {
      nn1 = i;
      nn2 = i + 1;
    }
  }

  if ( nn1 != -1 ) {
    dp = Patm[nn1] - Patm[nn2];
    dp1 = Patm[nn1] - P;
    fact2 = dp1 / dp;
    fact1 = 1.0e0 - fact2;
    *T = fact1 * Tatm[nn1] + fact2 * Tatm[nn2];
    *H = fact1 * hatm[nn1] + fact2 * hatm[nn2];
  } else {
    fprintf ( stderr, "Pressure not in US std atm Range\n" );
    exit ( 1 );
  }
}

#undef nn

/*****************************************************************
 Inlet subroutines
 *****************************************************************/

typedef struct {
  double T2overT1;
} arg2_t;

static double _errfunct2 ( void *arg, double Ms ) {
  double tmp, M2, P2overP1, T2overT1, U2overU1, Pstag2overPstag1;

  FindNormalShockProps ( Ms, 1.4e0, &M2, &P2overP1, &T2overT1, &U2overU1, &Pstag2overPstag1 );
  tmp = T2overT1 - ( ( arg2_t * ) arg )->T2overT1;
  return ( tmp );
}

void FindMsExtCompInlet ( double T1, double T3, double gam, double *Ms ) {
  long IFLAG;
  arg2_t arg;
  arg.T2overT1 = sqrt ( T3 / T1 );
  *Ms = EXM_find_root_zero_in ( &( _errfunct2 ), &arg, 1.0e0, 20.0e0, 1.0e-12, 1.0e-12, &IFLAG );
  if ( IFLAG == 4 ) {
    fprintf ( stderr, "couldn't find root in FindMsShcramjet, exiting..\n" );
    exit ( 1 );
  }
}

void FindMachAfterShock ( double M1, double M1n, double gam, double T1, double T2, double *M2 ) {
  double M1t2, M2n2, M2t2;
  M1t2 = sqr ( M1 ) - sqr ( M1n );
  M2n2 = ( sqr ( M1n ) + 2.0e0 / ( gam - 1.0e0 ) ) / ( 2.0e0 * gam / ( gam - 1.0e0 ) * sqr ( M1n ) - 1.0e0 );
  M2t2 = M1t2 * T1 / T2;
  *M2 = sqrt ( M2t2 + M2n2 );
}

void FindExtCompInletProps ( double M1, double L, double Rgas,
                             double T3, double Pdyn, double gam,
                             double *T1, double *T2, double *P1,
                             double *P2, double *P3, double *Ms,
                             double *M2, double *M3, double *W1, double *mdot, double *altitude ) {
  double xi;
  double U1;

  *P1 = 2.0e0 * Pdyn / gam / pow ( M1, 2 );
  FindTandHfromPstdatm ( *P1, altitude, T1 );
  FindMsExtCompInlet ( *T1, T3, gam, Ms );
  xi = ( 2.0e0 * gam / ( gam + 1.0e0 ) * sqr ( *Ms ) - ( gam - 1.0e0 ) / ( gam + 1.0e0 ) );
  *P2 = xi * ( *P1 );
  *P3 = xi * ( *P2 );
  ( *T2 ) = ( *T1 ) * ( 2.0e0 * gam / ( gam + 1.0e0 ) * sqr ( *Ms ) - ( gam - 1.0e0 ) / ( gam + 1.0e0 ) )
    * ( ( gam - 1.0e0 ) / ( gam + 1.0e0 ) + 2.0e0 / ( gam + 1.0e0 ) / sqr ( *Ms ) );
  FindMachAfterShock ( M1, *Ms, gam, *T1, *T2, M2 );
  FindMachAfterShock ( *M2, *Ms, gam, *T2, T3, M3 );
  U1 = M1 * sqrt ( gam * Rgas * ( *T1 ) );

  *W1 = L * tan ( asin ( *Ms / M1 ) );
  *mdot = ( *P1 ) * ( U1 ) / Rgas / ( *T1 ) * ( *W1 );
}
