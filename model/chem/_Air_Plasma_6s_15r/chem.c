// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2016 Bernard Parent

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

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
//#include <model/.active/model_eef.h>

#define nr 15

#define ATTACHMENT TRUE
#define TOWNSEND TRUE
#define TOWNSEND_IMPLICIT FALSE

#define Estarmin 1e-40

  /* Remember: 
     by definition
     W[specO2]=qi*Nk[specN2]/N  
     where 
     then
     Wp[13]=kf[13]*Nk[specN2]=qi*Nk[specN2]/N
     or 
     kf[13]=qi/N

   */



double _k1a ( double EoverN ) {
  double k, theta;
  theta = log ( EoverN );
  k = exp ( -0.0105809 * sqr ( theta ) - 2.40411e-75 * pow ( theta, 46.0 ) );
  return ( k );
}

double _k1b ( double EoverN ) {
  double k, theta;
  theta = log ( EoverN );
  k = exp ( -0.0102785 * sqr ( theta ) - 2.42260e-75 * pow ( theta, 46.0 ) );

  return ( k );
}

void find_W ( spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long spec;
  double kf[nr];
  double Wp[nr];
  double Nk[ns];
  spec_t calM;
  double N;

  Estar = max ( Estarmin, Estar );

  for ( spec = 0; spec < ns; spec++ )
    W[spec] = 0.0;
  /* find the gas temperature and the electron temperature */
  /* find the mole fraction Xk in moles/cm^3 */
  for ( spec = 0; spec < ns; spec++ ) {
    calM[spec] = _calM ( spec );
    Nk[spec] = rhok[spec] / calM[spec] * 1e-6 * calA;
  }
  N = 0.0;
  for ( spec = 0; spec < ns; spec++ )
    N += Nk[spec];              /* N is in 1/cm3 */

  /* set the reaction rates for each reaction */
  if ( TOWNSEND ) {
    kf[0] = _k1a ( Estar );     /* 1a */
    kf[1] = _k1b ( Estar );     /* 1b */
  } else {
    kf[0] = 0.0;
    kf[1] = 0.0;
  }
  kf[2] = 2.0e-7 * pow ( 300.0 / Te, 0.7 );     /* 2a */
  kf[3] = 2.8e-7 * pow ( 300.0 / Te, 0.5 );     /* 2b */
  kf[4] = 2.0e-7 * pow ( 300.0 / T, 0.5 );      /* 3a */
  kf[5] = kf[4];                /* 3b */
  kf[6] = 2.0e-25 * pow ( 300.0 / T, 2.5 );     /* 4a */
  kf[7] = kf[6];                /* 4b */
  kf[8] = kf[6];                /* 4c */
  kf[9] = kf[6];                /* 4d */
  if ( ATTACHMENT ) {
    kf[10] = 1.4e-29 * 300.0 / Te * exp ( -600.0 / T ) * exp ( 700.0 * ( Te - T ) / Te / T );   /* 5a */
    kf[11] = 1.07e-31 * sqr ( 300.0 / Te ) * exp ( -70.0 / T ) * exp ( 1500.0 * ( Te - T ) / Te / T );  /* 5b */
  } else {
    kf[10] = 0.0;
    kf[11] = 0.0;
  }
  kf[12] = 8.6e-10 * exp ( -6030.0 / T ) * ( 1.0 - exp ( -1570.0 / T ) );       /* 6 */
  kf[13] = 2.0e11 * Qbeam / N;  /* 7a *//* Qbeam is in watts per m3, N in 1/cm3 */
  kf[14] = 1.8e11 * Qbeam / N;  /* 7b */

  /* reaction 7b */
/* N2  ->  N2+  +  e- */
  Wp[14] = kf[14] * Nk[specN2];
  W[specN2] += -Wp[14];
  W[specN2plus] += +Wp[14];
  W[speceminus] += +Wp[14];

  /* reaction 7a */
/* O2  ->  O2+  +  e-  */
  Wp[13] = kf[13] * Nk[specO2];
  W[specO2] += -Wp[13];
  W[specO2plus] += +Wp[13];
  W[speceminus] += +Wp[13];

  /* reaction 6 */
/* O2  +  O2-  ->  2*O2   +   e-   */
  Wp[12] = kf[12] * Nk[specO2] * Nk[specO2minus];
  W[specO2] += +Wp[12];
  W[specO2minus] += -Wp[12];
  W[speceminus] += +Wp[12];

  /* reaction 5b */
/* O2  +  N2  +  e-   ->  N2  +  O2-  */
  Wp[11] = kf[11] * Nk[speceminus] * Nk[specO2] * Nk[specN2];
  W[specO2] += -Wp[11];
  W[specO2minus] += +Wp[11];
  W[speceminus] += -Wp[11];

  /* reaction 5a */
/*  2*O2  +  e-  ->  O2  +  O2-   */
  Wp[10] = kf[10] * Nk[speceminus] * sqr ( Nk[specO2] );
  W[specO2] += -Wp[10];
  W[specO2minus] += +Wp[10];
  W[speceminus] += -Wp[10];

  /* reaction 4d */
/* O2  +  O2+   +  O2-  ->  3*O2   */
  Wp[9] = kf[9] * Nk[specO2] * Nk[specO2plus] * Nk[specO2minus];
  W[specO2] += +2.0 * Wp[9];
  W[specO2plus] += -Wp[9];
  W[specO2minus] += -Wp[9];

  /* reaction 4c */
/* O2  +   N2+  +  O2-   ->  2*O2  + N2  */
  Wp[8] = kf[8] * Nk[specO2] * Nk[specN2plus] * Nk[specO2minus];
  W[specO2] += +Wp[8];
  W[specN2] += +Wp[8];
  W[specN2plus] += -Wp[8];
  W[specO2minus] += -Wp[8];

  /* reaction 4b */
/*  N2  +  O2+  +  O2-  -> 2*O2  + N2  */
  Wp[7] = kf[7] * Nk[specN2] * Nk[specO2plus] * Nk[specO2minus];
  W[specO2] += +2.0 * Wp[7];
  W[specO2plus] += -Wp[7];
  W[specO2minus] += -Wp[7];

  /* reaction 4a */
/*  N2  +  N2+  +  O2-   ->  O2  + 2*N2 */
  Wp[6] = kf[6] * Nk[specN2] * Nk[specN2plus] * Nk[specO2minus];
  W[specO2] += +Wp[6];
  W[specN2] += +Wp[6];
  W[specN2plus] += -Wp[6];
  W[specO2minus] += -Wp[6];

  /* reaction 3b */
/* O2-  +  O2+  ->  2*O2 */
  Wp[5] = kf[5] * Nk[specO2plus] * Nk[specO2minus];
  W[specO2] += +2.0 * Wp[5];
  W[specO2plus] += -Wp[5];
  W[specO2minus] += -Wp[5];

  /* reaction 3a */
/*  O2-  +  N2+  ->  O2  +  N2 */
  Wp[4] = kf[4] * Nk[specN2plus] * Nk[specO2minus];
  W[specO2] += +Wp[4];
  W[specN2] += +Wp[4];
  W[specN2plus] += -Wp[4];
  W[specO2minus] += -Wp[4];

  /* reaction 2b */
/* e-  +  N2+  ->  N2  */
  Wp[3] = kf[3] * Nk[speceminus] * Nk[specN2plus];
  W[specN2] += +Wp[3];
  W[specN2plus] += -Wp[3];
  W[speceminus] += -Wp[3];

  /* reaction 2a */
/* e-  +  O2+  -> O2  */
  Wp[2] = kf[2] * Nk[speceminus] * Nk[specO2plus];
  W[specO2] += +Wp[2];
  W[specO2plus] += -Wp[2];
  W[speceminus] += -Wp[2];

  /* reaction 1b */
/* e-  +  O2  ->  O2+   +  2*e-  */
  Wp[1] = kf[1] * Nk[speceminus] * Nk[specO2];
  W[specO2] += -Wp[1];
  W[specO2plus] += +Wp[1];
  W[speceminus] += +Wp[1];

  /* reaction 1a */
/* e-  +  N2  ->  N2+  +  2*e-  */
  Wp[0] = kf[0] * Nk[speceminus] * Nk[specN2];
  W[specN2] += -Wp[0];
  W[specN2plus] += +Wp[0];
  W[speceminus] += +Wp[0];

  for ( spec = 0; spec < ns; spec++ )
    W[spec] = W[spec] / calA * calM[spec] * 1.0e6;
}

void find_dW_dx ( spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, r, s, spec;           /* counters */
  spec_t calM;
  double Nk[ns];
  double dkdTe, dkdT, dkdQb;
  double N;
  double kf[nr];
  double sigma, theta1, theta2, dktownsenddNk, ktownsend;

  Estar = max ( Estarmin, Estar );

  /* first, initialize all derivatives to zero */
  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }

  /* find the mole fraction Xk in moles/cm^3 */
  for ( spec = 0; spec < ns; spec++ ) {
    calM[spec] = _calM ( spec );
    Nk[spec] = rhok[spec] / calM[spec] * 1e-6 * calA;
  }
  N = 0.0;
  for ( spec = 0; spec < ns; spec++ )
    N += Nk[spec];              /* N is in 1/cm3 */

  /* set the reaction rates for each reaction */

  /* reaction 1a */
/* e-  +  N2  ->  N2+  +  2*e-  */

//  Wp[0]=kf[0]*Nk[ns]*Nk[specN2];
//  W[ns]+=+Wp[0];
//  W[specN2]+=-Wp[0];
//  W[specN2plus]+=+Wp[0];
//  W[speceminus]+=+Wp[0];

  if ( TOWNSEND && TOWNSEND_IMPLICIT ) {
    theta1 = 0.0105809;
    theta2 = 2.40411e-75;
    ktownsend = exp ( -theta1 * sqr ( log ( Estar ) ) - theta2 * pow ( log ( Estar ), 46.0 ) ) * 1e-6;
    sigma = 0.0;
    for ( spec = 0; spec < ns; spec++ )
      sigma += fabs ( _C ( spec ) ) * mu[spec] * Nk[spec] * 1e6;

    for ( spec = 0; spec < ns; spec++ ) {
      dktownsenddNk =
        ktownsend * fabs ( _C ( spec ) ) * mu[spec] * ( 2.0 * theta1 / notzero ( sigma,1e-99 ) * log ( Estar ) +
                                                        46.0 * theta2 / notzero ( sigma,1e-99 ) *
                                                        pow ( log ( Estar ), 45.0 ) );
      dWdrhok[specN2][spec] += -dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
      dWdrhok[specN2plus][spec] += dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
      dWdrhok[speceminus][spec] += dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
    }

  }

  /* reaction 1b */
/* e-  +  O2  ->  O2+   +  2*e-  */

//  Wp[1]=kf[1]*Nk[ns]*Nk[specO2];
//  W[specO2]+=-Wp[1];
//  W[specO2plus]+=+Wp[1];
//  W[speceminus]+=+Wp[1];

  if ( TOWNSEND && TOWNSEND_IMPLICIT ) {
    theta1 = 0.0102785;
    theta2 = 2.42260e-75;
    ktownsend = exp ( -theta1 * sqr ( log ( Estar ) ) - theta2 * pow ( log ( Estar ), 46.0 ) ) * 1e-6;
    sigma = 0.0;
    for ( spec = 0; spec < ns; spec++ )
      sigma += fabs ( _C ( spec ) ) * mu[spec] * Nk[spec] * 1e6;

    for ( spec = 0; spec < ns; spec++ ) {
      dktownsenddNk =
        ktownsend * fabs ( _C ( spec ) ) * mu[spec] * ( 2.0 * theta1 / notzero ( sigma,1e-99 ) * log ( Estar ) +
                                                        46.0 * theta2 / notzero ( sigma,1e-99 ) *
                                                        pow ( log ( Estar ), 45.0 ) );
      dWdrhok[specO2][spec] += -dktownsenddNk * Nk[specO2] * Nk[speceminus] * 1e12;
      dWdrhok[specO2plus][spec] += dktownsenddNk * Nk[specO2] * Nk[speceminus] * 1e12;
      dWdrhok[speceminus][spec] += dktownsenddNk * Nk[specO2] * Nk[speceminus] * 1e12;
    }

  }

  /* reaction 2a */
/* e-  +  O2+  -> O2  */
  kf[2] = 2.0e-7 * pow ( 300.0 / Te, 0.7 );     /* 2a */

//  Wp[2]=kf[2]*Nk[speceminus]*Nk[specO2plus];
//  W[specO2]+=+Wp[2];
//  W[specO2plus]+=-Wp[2];
//  W[speceminus]+=-Wp[2];

  dWdrhok[specO2][speceminus] += +kf[2] * Nk[specO2plus];
  dWdrhok[specO2][specO2plus] += +kf[2] * Nk[speceminus];

  dWdrhok[specO2plus][speceminus] += -kf[2] * Nk[specO2plus];
  dWdrhok[specO2plus][specO2plus] += -kf[2] * Nk[speceminus];

  dWdrhok[speceminus][speceminus] += -kf[2] * Nk[specO2plus];
  dWdrhok[speceminus][specO2plus] += -kf[2] * Nk[speceminus];

  dkdTe = -0.7 * 2.0e-7 * pow ( Te / 300.0, -1.7 ) / 300.0 * Nk[speceminus] * Nk[specO2plus];
  dWdTe[speceminus] += -dkdTe;
  dWdTe[specO2] += +dkdTe;
  dWdTe[specO2plus] += -dkdTe;

  /* reaction 2b */
/* e-  +  N2+  ->  N2  */

  kf[3] = 2.8e-7 * pow ( 300.0 / Te, 0.5 );     /* 2b */

//  Wp[3]=kf[3]*Nk[speceminus]*Nk[specN2plus];
//  W[specN2]+=+Wp[3];
//  W[specN2plus]+=-Wp[3];
//  W[speceminus]+=-Wp[3];

  dWdrhok[specN2][speceminus] += +kf[3] * Nk[specN2plus];
  dWdrhok[specN2][specN2plus] += +kf[3] * Nk[speceminus];
  dWdrhok[specN2plus][speceminus] += -kf[3] * Nk[specN2plus];
  dWdrhok[specN2plus][specN2plus] += -kf[3] * Nk[speceminus];
  dWdrhok[speceminus][speceminus] += -kf[3] * Nk[specN2plus];
  dWdrhok[speceminus][specN2plus] += -kf[3] * Nk[speceminus];

  dkdTe = -0.5 * 2.8e-7 * pow ( Te / 300.0, -1.5 ) / 300.0 * Nk[speceminus] * Nk[specN2plus];
  dWdTe[specN2] += +dkdTe;
  dWdTe[specN2plus] += -dkdTe;
  dWdTe[speceminus] += -dkdTe;

  /* reaction 3a */
/*  O2-  +  N2+  ->  O2  +  N2 */

  kf[4] = 2.0e-7 * pow ( 300.0 / T, 0.5 );      /* 3a */
//  Wp[4]=kf[4]*Nk[specN2plus]*Nk[specO2minus];
//  W[specO2]+=+Wp[4];
//  W[specN2]+=+Wp[4];
//  W[specN2plus]+=-Wp[4];
//  W[specO2minus]+=-Wp[4];

  dWdrhok[specO2][specN2plus] += +kf[4] * Nk[specO2minus];
  dWdrhok[specO2][specO2minus] += +kf[4] * Nk[specN2plus];
  dWdrhok[specN2][specN2plus] += +kf[4] * Nk[specO2minus];
  dWdrhok[specN2][specO2minus] += +kf[4] * Nk[specN2plus];
  dWdrhok[specN2plus][specN2plus] += -kf[4] * Nk[specO2minus];
  dWdrhok[specN2plus][specO2minus] += -kf[4] * Nk[specN2plus];
  dWdrhok[specO2minus][specN2plus] += -kf[4] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2minus] += -kf[4] * Nk[specN2plus];
  dkdT = -0.5 * 2.0e-7 * pow ( T / 300.0, -1.5 ) / 300.0 * Nk[specN2plus] * Nk[specO2minus];
  dWdT[specO2] += +dkdT;
  dWdT[specN2] += +dkdT;
  dWdT[specN2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 3b */
/* O2-  +  O2+  ->  2*O2 */

  kf[5] = kf[4];
//  Wp[5]=kf[5]*Nk[specO2plus]*Nk[specO2minus];
//  W[specO2]+=+2.0*Wp[5];
//  W[specO2plus]+=-Wp[5];
//  W[specO2minus]+=-Wp[5];
  dWdrhok[specO2][specO2plus] += +2.0 * kf[5] * Nk[specO2minus];
  dWdrhok[specO2][specO2minus] += +2.0 * kf[5] * Nk[specO2plus];
  dWdrhok[specO2plus][specO2plus] += -kf[5] * Nk[specO2minus];
  dWdrhok[specO2plus][specO2minus] += -kf[5] * Nk[specO2plus];
  dWdrhok[specO2minus][specO2plus] += -kf[5] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2minus] += -kf[5] * Nk[specO2plus];
  dkdT = -0.5 * 2.0e-7 * pow ( T / 300.0, -1.5 ) / 300.0 * Nk[specO2plus] * Nk[specO2minus];
  dWdT[specO2] += +2.0 * dkdT;
  dWdT[specO2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4a */
/*  N2  +  N2+  +  O2-   ->  O2  + 2*N2 */

  kf[6] = 2.0e-25 * pow ( 300.0 / T, 2.5 );     /* 4a */
//  Wp[6]=kf[6]*Nk[specN2]*Nk[specN2plus]*Nk[specO2minus];
//  W[specO2]+=+Wp[6];
//  W[specN2]+=+Wp[6];
//  W[specN2plus]+=-Wp[6];
//  W[specO2minus]+=-Wp[6];

  dWdrhok[specO2][specN2] += +kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specO2][specN2plus] += +kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdrhok[specO2][specO2minus] += +kf[6] * Nk[specN2] * Nk[specN2plus];
  dWdrhok[specN2][specN2] += +kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specN2][specN2plus] += +kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdrhok[specN2][specO2minus] += +kf[6] * Nk[specN2] * Nk[specN2plus];
  dWdrhok[specN2plus][specN2] += -kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specN2plus][specN2plus] += -kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdrhok[specN2plus][specO2minus] += -kf[6] * Nk[specN2] * Nk[specN2plus];
  dWdrhok[specO2minus][specN2] += -kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specO2minus][specN2plus] += -kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2minus] += -kf[6] * Nk[specN2] * Nk[specN2plus];
  dkdT = -2.5 * kf[6] / T * Nk[specN2] * Nk[specN2plus] * Nk[specO2minus];
  dWdT[specO2] += +dkdT;
  dWdT[specN2] += +dkdT;
  dWdT[specN2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4b */
/*  N2  +  O2+  +  O2-  -> 2*O2  + N2  */
  kf[7] = kf[6];                /* 4b */
//  Wp[7]=kf[7]*Nk[specN2]*Nk[specO2plus]*Nk[specO2minus];
//  W[specO2]+=+2.0*Wp[7];
//  W[specO2plus]+=-Wp[7];
//  W[specO2minus]+=-Wp[7];
  dWdrhok[specO2][specN2] += +2.0 * kf[7] * Nk[specO2plus] * Nk[specO2minus];
  dWdrhok[specO2][specO2plus] += +2.0 * kf[7] * Nk[specN2] * Nk[specO2minus];
  dWdrhok[specO2][specO2minus] += +2.0 * kf[7] * Nk[specN2] * Nk[specO2plus];
  dWdrhok[specO2plus][specN2] += -kf[7] * Nk[specO2plus] * Nk[specO2minus];
  dWdrhok[specO2plus][specO2plus] += -kf[7] * Nk[specN2] * Nk[specO2minus];
  dWdrhok[specO2plus][specO2minus] += -kf[7] * Nk[specN2] * Nk[specO2plus];
  dWdrhok[specO2minus][specN2] += -kf[7] * Nk[specO2plus] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2plus] += -kf[7] * Nk[specN2] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2minus] += -kf[7] * Nk[specN2] * Nk[specO2plus];
  dkdT = -2.5 * kf[7] / T * Nk[specN2] * Nk[specO2plus] * Nk[specO2minus];
  dWdT[specO2] += +2.0 * dkdT;  /* 2.0e-25*pow(300.0/T,2.5) */
  dWdT[specO2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4c */
/* O2  +   N2+  +  O2-   ->  2*O2  + N2  */
  kf[8] = kf[6];                /* 4c */
//  Wp[8]=kf[8]*Nk[specO2]*Nk[specN2plus]*Nk[specO2minus];
//  W[specO2]+=+Wp[8];
//  W[specN2]+=+Wp[8];
//  W[specN2plus]+=-Wp[8];
//  W[specO2minus]+=-Wp[8];

  dWdrhok[specO2][specO2] += +kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specO2][specN2plus] += +kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdrhok[specO2][specO2minus] += +kf[8] * Nk[specO2] * Nk[specN2plus];
  dWdrhok[specN2][specO2] += +kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specN2][specN2plus] += +kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdrhok[specN2][specO2minus] += +kf[8] * Nk[specO2] * Nk[specN2plus];
  dWdrhok[specN2plus][specO2] += -kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specN2plus][specN2plus] += -kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdrhok[specN2plus][specO2minus] += -kf[8] * Nk[specO2] * Nk[specN2plus];
  dWdrhok[specO2minus][specO2] += -kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdrhok[specO2minus][specN2plus] += -kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2minus] += -kf[8] * Nk[specO2] * Nk[specN2plus];
  dkdT = -2.5 * kf[8] / T * Nk[specO2] * Nk[specN2plus] * Nk[specO2minus];
  dWdT[specO2] += +dkdT;
  dWdT[specN2] += +dkdT;
  dWdT[specN2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4d */
/* O2  +  O2+   +  O2-  ->  3*O2   */
  kf[9] = kf[6];                /* 4d */
//  Wp[9]=kf[9]*Nk[specO2]*Nk[specO2plus]*Nk[specO2minus];
//  W[specO2]+=+2.0*Wp[9];
//  W[specO2plus]+=-Wp[9];
//  W[specO2minus]+=-Wp[9];

  dWdrhok[specO2][specO2] += +2.0 * kf[9] * Nk[specO2plus] * Nk[specO2minus];
  dWdrhok[specO2][specO2plus] += +2.0 * kf[9] * Nk[specO2] * Nk[specO2minus];
  dWdrhok[specO2][specO2minus] += +2.0 * kf[9] * Nk[specO2] * Nk[specO2plus];
  dWdrhok[specO2plus][specO2] += -kf[9] * Nk[specO2plus] * Nk[specO2minus];
  dWdrhok[specO2plus][specO2plus] += -kf[9] * Nk[specO2] * Nk[specO2minus];
  dWdrhok[specO2plus][specO2minus] += -kf[9] * Nk[specO2] * Nk[specO2plus];
  dWdrhok[specO2minus][specO2] += -kf[9] * Nk[specO2plus] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2plus] += -kf[9] * Nk[specO2] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2minus] += -kf[9] * Nk[specO2] * Nk[specO2plus];
  dkdT = -2.5 * kf[9] / T * Nk[specO2] * Nk[specO2plus] * Nk[specO2minus];
  dWdT[specO2] += +2.0 * dkdT;  /* 2.0e-25*pow(300.0/T,2.5) */
  dWdT[specO2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 5a */
/*  2*O2  +  e-  ->  O2  +  O2-   */
  if ( ATTACHMENT ) {
    kf[10] = 1.4e-29 * 300.0 / Te * exp ( -600.0 / T ) * exp ( 700.0 * ( Te - T ) / Te / T );   /* 5a */
  } else {
    kf[10] = 0.0;
  }
//  Wp[10]=kf[10]*Nk[speceminus]*sqr(Nk[specO2]);
//  W[specO2]+=-Wp[10];
//  W[specO2minus]+=+Wp[10];
//  W[speceminus]+=-Wp[10];

  dWdrhok[specO2][speceminus] += -kf[10] * sqr ( Nk[specO2] );
  dWdrhok[specO2][specO2] += -kf[10] * Nk[speceminus] * 2.0 * Nk[specO2];
  dWdrhok[specO2minus][speceminus] += +kf[10] * sqr ( Nk[specO2] );
  dWdrhok[specO2minus][specO2] += +kf[10] * Nk[speceminus] * 2.0 * Nk[specO2];
  dWdrhok[speceminus][speceminus] += -kf[10] * sqr ( Nk[specO2] );
  dWdrhok[speceminus][specO2] += -kf[10] * Nk[speceminus] * 2.0 * Nk[specO2];

  dkdTe = ( -1.0 / Te + 700.0 / sqr ( Te ) ) * kf[10] * Nk[speceminus] * sqr ( Nk[specO2] );
  dkdT =
    ( 600.0 / sqr ( T ) - 700.0 / Te / T -
      700.0 * ( Te - T ) / Te / sqr ( T ) ) * kf[10] * Nk[speceminus] * sqr ( Nk[specO2] );
  dWdTe[speceminus] += -dkdTe;
  dWdTe[specO2] += -dkdTe;
  dWdTe[specO2minus] += +dkdTe;
  dWdT[speceminus] += -dkdT;
  dWdT[specO2] += -dkdT;
  dWdT[specO2minus] += +dkdT;

  /* reaction 5b */
/* O2  +  N2  +  e-   ->  N2  +  O2-  */
  if ( ATTACHMENT ) {
    kf[11] = 1.07e-31 * sqr ( 300.0 / Te ) * exp ( -70.0 / T ) * exp ( 1500.0 * ( Te - T ) / Te / T );  /* 5b */
  } else {
    kf[11] = 0.0;
  }
//  Wp[11]=kf[11]*Nk[speceminus]*Nk[specO2]*Nk[specN2];
//  W[specO2]+=-Wp[11];
//  W[specO2minus]+=+Wp[11];
//  W[speceminus]+=-Wp[11];

  dWdrhok[specO2][speceminus] += -kf[11] * Nk[specO2] * Nk[specN2];
  dWdrhok[specO2][specO2] += -kf[11] * Nk[speceminus] * Nk[specN2];
  dWdrhok[specO2][specN2] += -kf[11] * Nk[speceminus] * Nk[specO2];

  dWdrhok[specO2minus][speceminus] += +kf[11] * Nk[specO2] * Nk[specN2];
  dWdrhok[specO2minus][specO2] += +kf[11] * Nk[speceminus] * Nk[specN2];
  dWdrhok[specO2minus][specN2] += +kf[11] * Nk[speceminus] * Nk[specO2];

  dWdrhok[speceminus][speceminus] += -kf[11] * Nk[specO2] * Nk[specN2];
  dWdrhok[speceminus][specO2] += -kf[11] * Nk[speceminus] * Nk[specN2];
  dWdrhok[speceminus][specN2] += -kf[11] * Nk[speceminus] * Nk[specO2];

  dkdTe = ( -2.0 / Te + 1500.0 / sqr ( Te ) ) * kf[11] * Nk[speceminus] * Nk[specO2] * Nk[specN2];
  dkdT = ( 70.0 / sqr ( T ) - 1500.0 / sqr ( T ) ) * kf[11] * Nk[speceminus] * Nk[specO2] * Nk[specN2];
  dWdTe[specO2] += -dkdTe;
  dWdTe[specO2minus] += +dkdTe;
  dWdTe[speceminus] += -dkdTe;
  dWdT[specO2] += -dkdT;
  dWdT[specO2minus] += +dkdT;
  dWdT[speceminus] += -dkdT;

  /* reaction 6 */
/* O2  +  O2-  ->  2*O2   +   e-   */
  kf[12] = 8.6e-10 * exp ( -6030.0 / T ) * ( 1.0 - exp ( -1570.0 / T ) );       /* 6 */
//  Wp[12]=kf[12]*Nk[specO2]*Nk[specO2minus];
//  W[specO2]+=+Wp[12];
//  W[specO2minus]+=-Wp[12];
//  W[speceminus]+=+Wp[12];
  dWdrhok[specO2][specO2] += +kf[12] * Nk[specO2minus];
  dWdrhok[specO2][specO2minus] += +kf[12] * Nk[specO2];
  dWdrhok[specO2minus][specO2] += -kf[12] * Nk[specO2minus];
  dWdrhok[specO2minus][specO2minus] += -kf[12] * Nk[specO2];
  dWdrhok[speceminus][specO2] += +kf[12] * Nk[specO2minus];
  dWdrhok[speceminus][specO2minus] += +kf[12] * Nk[specO2];

  dkdT =
    ( kf[12] * 6030.0 / sqr ( T ) - 8.6e-10 * exp ( -6030.0 / T ) * exp ( -1570.0 / T ) * 1570.0 / sqr ( T )
     ) * Nk[specO2] * Nk[specO2minus];

  dWdT[specO2] += +dkdT;
  dWdT[specO2minus] += -dkdT;
  dWdT[speceminus] += +dkdT;

  /* reaction 7a */
/* O2  ->  O2+  +  e-  */
  kf[13] = 2.0e11 * Qbeam / N;  /* 7a *//* Qbeam is in watts per m3, N in 1/cm3 */
//  Wp[13]=kf[13]*Nk[specO2];
//  W[specO2]+=-Wp[13];
//  W[specO2plus]+=+Wp[13];  
//  W[speceminus]+=+Wp[13];  

  dWdrhok[specO2][specO2] += -kf[13];
  dWdrhok[specO2plus][specO2] += +kf[13];
  dWdrhok[speceminus][specO2] += +kf[13];
  for ( s = 0; s < ns; s++ ) {
    dWdrhok[specO2][s] += +kf[13] * Nk[specO2] / N;
    dWdrhok[specO2plus][s] += -kf[13] * Nk[specO2] / N;
    dWdrhok[speceminus][s] += -kf[13] * Nk[specO2] / N;
  }
  dkdQb = 2.0e11 / N * Nk[specO2];
  dWdQbeam[specO2] += -dkdQb;
  dWdQbeam[specO2plus] += +dkdQb;
  dWdQbeam[speceminus] += +dkdQb;

  /* reaction 7b */
/* N2  ->  N2+  +  e- */
  kf[14] = 1.8e11 * Qbeam / N;  /* 7b */
//  Wp[14]=kf[14]*Nk[specN2];
//  W[specN2]+=-Wp[14];
//  W[specN2plus]+=+Wp[14];  
//  W[speceminus]+=+Wp[14];  
  dWdrhok[specN2][specN2] += -kf[14];
  dWdrhok[specN2plus][specN2] += +kf[14];
  dWdrhok[speceminus][specN2] += +kf[14];
  for ( s = 0; s < ns; s++ ) {
    dWdrhok[specN2][s] += +kf[14] * Nk[specN2] / N;
    dWdrhok[specN2plus][s] += -kf[14] * Nk[specN2] / N;
    dWdrhok[speceminus][s] += -kf[14] * Nk[specN2] / N;
  }
  dkdQb = 1.8e11 / N * Nk[specN2];
  dWdQbeam[specN2] += -dkdQb;
  dWdQbeam[specN2plus] += +dkdQb;
  dWdQbeam[speceminus] += +dkdQb;

  for ( s = 0; s < ns; s++ ) {
    for ( r = 0; r < ns; r++ ) {
      dWdrhok[s][r] = dWdrhok[s][r] * calM[s] / calM[r];
    }
    dWdT[s] = dWdT[s] / calA * calM[s] * 1.0e6;
    dWdTe[s] = dWdTe[s] / calA * calM[s] * 1.0e6;
    dWdQbeam[s] = dWdQbeam[s] / calA * calM[s] * 1.0e6;
  }

}
