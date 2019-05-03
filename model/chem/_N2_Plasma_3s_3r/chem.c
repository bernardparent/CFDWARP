// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015 Bernard Parent

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

#define nr 3

#define ATTACHMENT TRUE
#define TOWNSEND TRUE
#define TOWNSEND_SEMI_IMPLICIT FALSE    //jacobian found through "constant current" approach
#define TOWNSEND_IMPLICIT FALSE  //jacobian found through exact linearization

#define Estarmin 1e-40


double _k1a ( double EoverN ) {
  double k, theta;
  theta = log ( EoverN );
  k = exp ( -0.0105809 * sqr ( theta ) - 2.40411e-75 * pow ( theta, 46.0 ) );
  return ( k );
}

double _dk1adN ( double EoverN, spec_t Nk ) {
  double theta, dk1adN, k1a, N;
  long spec;
  theta = log ( EoverN );
  k1a = exp ( -0.0105809 * sqr ( theta ) - 2.40411e-75 * pow ( theta, 46.0 ) );
  N = 0.0;
  for ( spec = 0; spec < ns; spec++ )
    N += Nk[spec];
  dk1adN =
    k1a * ( -0.0105809 * 2.0 * theta -
            2.40411e-75 * 46.0 * pow ( theta, 45.0 ) ) * ( 1.0 / EoverN ) * ( -EoverN / N );
  return ( dk1adN );
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
  } else {
    kf[0] = 0.0;
  }

  kf[1] = 1.8e11 * Qbeam / N;   /* 7b */
  kf[2] = 2.8e-7 * pow ( 300.0 / Te, 0.5 );     /* 2b */

  /* reaction 7b */
//       N2  N2+ e-  
// MR=  {1,  0,  0 },  /* 7b */
// MP=  {0,  1,  1 },  /* 7b */
  Wp[1] = kf[1] * Nk[specN2];
  W[specN2] += -Wp[1];
  W[specN2plus] += +Wp[1];
  W[speceminus] += +Wp[1];

  /* reaction 2b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 2b */
// MP=  {0,  1,  0,  0,  0,  0,  0,  0 },  /* 2b */
  Wp[2] = kf[2] * Nk[speceminus] * Nk[specN2plus];
  W[specN2] += +Wp[2];
  W[specN2plus] += -Wp[2];
  W[speceminus] += -Wp[2];

  /* reaction 1a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  1,  0,  0,  0,  0,  0,  1   },  /* 1a */
// MP=  {0,  0,  0,  0,  0,  1,  0,  2  },  /* 1a */
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
  spec2_t dWdNk;
  double Nk[ns];
  double dkdTe, dkdQb;
  double N, dkdN;
  double kf[nr];
  double theta1, theta2, sigma, dktownsenddNk, ktownsend;

  /* first, initialize all derivatives to zero */
  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdNk[s][k] = 0.0;
    }
  }

  /* make sure Estar is not zero */
  Estar = max ( Estarmin, Estar );

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
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  1,  0,  0,  0,  0,  0,  1   },  /* 1a */
// MP=  {0,  0,  0,  0,  0,  1,  0,  2  },  /* 1a */
//  Wp[0]=kf[0]*Nk[speceminus]*Nk[specN2];
//  W[specN2]+=-Wp[0];
//  W[specN2plus]+=+Wp[0];
//  W[speceminus]+=+Wp[0];

  if ( TOWNSEND && TOWNSEND_IMPLICIT ) {
    kf[0] = _k1a ( Estar );
    dkdN = _dk1adN ( Estar, Nk );
    dWdNk[specN2][speceminus] += -kf[0] * Nk[specN2];
    dWdNk[specN2][specN2] += -kf[0] * Nk[speceminus];

    dWdNk[specN2plus][speceminus] += kf[0] * Nk[specN2];
    dWdNk[specN2plus][specN2] += kf[0] * Nk[speceminus];

    dWdNk[speceminus][speceminus] += kf[0] * Nk[specN2];
    dWdNk[speceminus][specN2] += kf[0] * Nk[speceminus];

    for ( spec = 0; spec < ns; spec++ ) {
      dWdNk[speceminus][spec] += dkdN * Nk[speceminus] * Nk[specN2];
      dWdNk[specN2][spec] -= dkdN * Nk[speceminus] * Nk[specN2];
      dWdNk[specN2plus][spec] += dkdN * Nk[speceminus] * Nk[specN2];
    }
  }

  if ( TOWNSEND && TOWNSEND_SEMI_IMPLICIT ) {

    theta1 = 0.0105809;
    theta2 = 2.40411e-75;
    ktownsend = exp ( -theta1 * sqr ( log ( Estar ) ) - theta2 * pow ( log ( Estar ), 46.0 ) ) * 1e-6;
    sigma = 0.0;
    for ( spec = 0; spec < ns; spec++ )
      sigma += fabs ( _C ( spec ) ) * mu[spec] * Nk[spec] * 1e6;

    for ( spec = 0; spec < ns; spec++ ) {
      dktownsenddNk =
        ktownsend * fabs ( _C ( spec ) ) * mu[spec] * ( 2.0 * theta1 / notzero ( sigma,1e-99 ) * log ( Estar ) +
                                                        46.0 * theta2 / notzero ( sigma ,1e-99) *
                                                        pow ( log ( Estar ), 45.0 ) );
      dWdNk[specN2][spec] += -dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
      dWdNk[specN2plus][spec] += dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
      dWdNk[speceminus][spec] += dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
    }
  }

  /* reaction 2b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 2b */
// MP=  {0,  1,  0,  0,  0,  0,  0,  0 },  /* 2b */

  kf[2] = 2.8e-7 * pow ( 300.0 / Te, 0.5 );     /* 2b */

//  Wp[2]=kf[2]*Nk[speceminus]*Nk[specN2plus];
//  W[specN2]+=+Wp[2];
//  W[specN2plus]+=-Wp[2];
//  W[speceminus]+=-Wp[2];

  dWdNk[specN2][speceminus] += +kf[2] * Nk[specN2plus];
  dWdNk[specN2][specN2plus] += +kf[2] * Nk[speceminus];
  dWdNk[specN2plus][speceminus] += -kf[2] * Nk[specN2plus];
  dWdNk[specN2plus][specN2plus] += -kf[2] * Nk[speceminus];
  dWdNk[speceminus][speceminus] += -kf[2] * Nk[specN2plus];
  dWdNk[speceminus][specN2plus] += -kf[2] * Nk[speceminus];

  dkdTe = -0.5 * 2.8e-7 * pow ( Te / 300.0, -1.5 ) / 300.0 * Nk[speceminus] * Nk[specN2plus];
  dWdTe[specN2] += +dkdTe;
  dWdTe[specN2plus] += -dkdTe;
  dWdTe[speceminus] += -dkdTe;

  /* reaction 7b */
//       O2  N2  O   N   O2+ N2+ O2- e-  
// MR=  {0,  1,  0,  0,  0,  0,  0,  0 },  /* 7b */
// MP=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 7b */
  kf[1] = 1.8e11 * Qbeam / N;   /* 7b */
//  Wp[14]=kf[1]*Nk[specN2];
//  W[specN2]+=-Wp[1];
//  W[specN2plus]+=+Wp[1];  
//  W[speceminus]+=+Wp[1];  

  dWdNk[specN2][specN2] += -kf[1];
  dWdNk[specN2plus][specN2] += +kf[1];
  dWdNk[speceminus][specN2] += +kf[1];
  for ( s = 0; s < ns; s++ ) {
    dWdNk[specN2][s] += +kf[1] * Nk[specN2] / N;
    dWdNk[specN2plus][s] += -kf[1] * Nk[specN2] / N;
    dWdNk[speceminus][s] += -kf[1] * Nk[specN2] / N;
  }
  dkdQb = 1.8e11 / N * Nk[specN2];
  dWdQbeam[specN2] += -dkdQb;
  dWdQbeam[specN2plus] += +dkdQb;
  dWdQbeam[speceminus] += +dkdQb;

  for ( s = 0; s < ns; s++ ) {
    for ( r = 0; r < ns; r++ ) {
      dWdrhok[s][r] = dWdNk[s][r] * calM[s] / calM[r];
    }
    dWdT[s] = dWdT[s] / calA * calM[s] * 1.0e6;
    dWdTe[s] = dWdTe[s] / calA * calM[s] * 1.0e6;
    dWdQbeam[s] = dWdQbeam[s] / calA * calM[s] * 1.0e6;
  }

}
