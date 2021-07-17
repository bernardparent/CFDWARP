// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2005-2018 Bernard Parent
Copyright 2021 Prasanna Thoguluva Rajendran
All rights reserved.


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
#include <model/share/chem_share.h>

#define nr 28

#define DISSOCIATION TRUE
#define ATTACHMENT TRUE
#define TOWNSEND TRUE

#define TOWNSEND_SEMI_IMPLICIT FALSE    //jacobian found through "constant current" approach
#define TOWNSEND_IMPLICIT FALSE  //jacobian found through exact linearization

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

double _dk1bdN ( double EoverN, spec_t Nk ) {
  double theta, dk1bdN, k1b, N;
  long spec;
  theta = log ( EoverN );
  k1b = exp ( -0.0102785 * sqr ( theta ) - 2.42260e-75 * pow ( theta, 46.0 ) );
  N = 0.0;
  for ( spec = 0; spec < ns; spec++ )
    N += Nk[spec];
  dk1bdN =
    k1b * ( -0.0102785 * 2.0 * theta -
            2.42260e-75 * 46.0 * pow ( theta, 45.0 ) ) * ( 1.0 / EoverN ) * ( -EoverN / N );
  return ( dk1bdN );
}



void find_W_Macheret2007old ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long spec, cnt;
  double kf[nr];
  double Wp[nr];
  double Nk[ns];
  spec_t calM;
  double N, tmp;

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
  if ( DISSOCIATION ) {
    tmp = exp ( -59380.0 / T ) * ( 1.0 - exp ( -2240.0 / T ) );
    kf[15] = 3.7e-8 * tmp;      /* 8a */
    kf[16] = 9.3e-9 * tmp;      /* 8b */
    kf[17] = 1.3e-7 * tmp;      /* 8c */
    tmp = exp ( -113200.0 / T ) * ( 1.0 - exp ( -3354.0 / T ) );
    kf[18] = 5.0e-8 * tmp;      /* 8d */
    kf[19] = 5.0e-8 * tmp;      /* 8e */
    kf[20] = 1.1e-7 * tmp;      /* 8f */
    kf[21] = 2.45e-31 * pow ( T, -0.63 );       /* 9a */
    kf[22] = 2.76e-34 * exp ( 720.0 / T );      /* 9b */
    kf[23] = 8.8e-31 * pow ( T, -0.63 );        /* 9c */
    tmp = 8.27e-34 * exp ( 500.0 / T );
    kf[24] = tmp;               /* 9d */
    kf[25] = tmp;               /* 9e */
    kf[26] = tmp;               /* 9f */
    kf[27] = tmp;               /* 9g */
  } else {
    for ( cnt = 15; cnt < 28; cnt++ ) {
      kf[cnt] = 0.0;
    }
  }

// reaction 9g 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  0,  3,  0,  0,  0,  0},  /* 9g */
//MP= 0,  1,  0,  1,  0,  0,  0,  0},  /* 9g */
  Wp[27] = kf[27] * Nk[specN] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[27];
  W[specN] += -2.0 * Wp[27];

// reaction 9f 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  1,  2,  0,  0,  0,  0},  /* 9f */
//MP= 0,  1,  1,  0,  0,  0,  0,  0},  /* 9f */
  Wp[26] = kf[26] * Nk[specO] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[26];
  W[specN] += -2.0 * Wp[26];

// reaction 9e 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  0,  2,  0,  0,  0,  0},  /* 9e */
//MP= 0,  2,  0,  0,  0,  0,  0,  0},  /* 9e */
  Wp[25] = kf[25] * Nk[specN2] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[25];
  W[specN] += -2.0 * Wp[25];

// reaction 9d 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  0,  2,  0,  0,  0,  0},  /* 9d */
//MP= 1,  1,  0,  0,  0,  0,  0,  0},  /* 9d */
  Wp[24] = kf[24] * Nk[specO2] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[24];
  W[specN] += -2.0 * Wp[24];

// reaction 9c 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  3,  0,  0,  0,  0,  0},  /* 9c */
//MP= 1,  0,  1,  0,  0,  0,  0,  0},  /* 9c */
  Wp[23] = kf[23] * Nk[specO] * Nk[specO] * Nk[specO];
  W[specO2] += +Wp[23];
  W[specO] += -2.0 * Wp[23];

// reaction 9b 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  2,  0,  0,  0,  0,  0},  /* 9b */
//MP= 1,  1,  0,  0,  0,  0,  0,  0},  /* 9b */
  Wp[22] = kf[22] * Nk[specN2] * Nk[specO] * Nk[specO];
  W[specO2] += +Wp[22];
  W[specO] += -2.0 * Wp[22];

// reaction 9a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  2,  0,  0,  0,  0,  0},  /* 9a */
//MP= 2,  0,  0,  0,  0,  0,  0,  0},  /* 9a */
  Wp[21] = kf[21] * Nk[specO2] * Nk[specO] * Nk[specO];
  W[specO2] += +Wp[21];
  W[specO] += -2.0 * Wp[21];

// reaction 8f 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  1,  0,  0,  0,  0,  0},  /* 8f */
//MP= 0,  0,  1,  2,  0,  0,  0,  0},  /* 8f */
  Wp[20] = kf[20] * Nk[specN2] * Nk[specO];
  W[specN2] += -Wp[20];
  W[specN] += +2.0 * Wp[20];

// reaction 8e 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  2,  0,  0,  0,  0,  0,  0},  /* 8e */
//MP= 0,  1,  0,  2,  0,  0,  0,  0},  /* 8e */
  Wp[19] = kf[19] * Nk[specN2] * Nk[specN2];
  W[specN2] += -Wp[19];
  W[specN] += +2.0 * Wp[19];

// reaction 8d 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  1,  0,  0,  0,  0,  0,  0},  /* 8d */
//MP= 1,  0,  0,  2,  0,  0,  0,  0},  /* 8d */
  Wp[18] = kf[18] * Nk[specO2] * Nk[specN2];
  W[specN2] += -Wp[18];
  W[specN] += +2.0 * Wp[18];

// reaction 8c 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  1,  0,  0,  0,  0,  0},  /* 8c */
//MP= 0,  0,  3,  0,  0,  0,  0,  0},  /* 8c */
  Wp[17] = kf[17] * Nk[specO2] * Nk[specO];
  W[specO2] += -Wp[17];
  W[specO] += +2.0 * Wp[17];

// reaction 8b 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  1,  0,  0,  0,  0,  0,  0},  /* 8b */
//MP= 0,  1,  2,  0,  0,  0,  0,  0},  /* 8b */
  Wp[16] = kf[16] * Nk[specO2] * Nk[specN2];
  W[specO2] += -Wp[16];
  W[specO] += +2.0 * Wp[16];

// reaction 8a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 2,  0,  0,  0,  0,  0,  0,  0},  /* 8a */
//MP= 1,  0,  2,  0,  0,  0,  0,  0},  /* 8a */
  Wp[15] = kf[15] * Nk[specO2] * Nk[specO2];
  W[specO2] += -Wp[15];
  W[specO] += +2.0 * Wp[15];

  /* reaction 7b */
//       O2  N2  O   N   O2+ N2+ O2- e-  
// MR=  {0,  1,  0,  0,  0,  0,  0,  0 },  /* 7b */
// MP=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 7b */
  Wp[14] = kf[14] * Nk[specN2];
  W[specN2] += -Wp[14];
  W[specN2plus] += +Wp[14];
  W[speceminus] += +Wp[14];

  /* reaction 7a */
//        O2  N2  O   N   O2+ N2+ O2- e-    
// MR=  { 1,  0,  0,  0,  0,  0,  0,  0 },  /* 7a */
// MP=  { 0,  0,  0,  0,  1,  0,  0,  1 },  /* 7a */
  Wp[13] = kf[13] * Nk[specO2];
  W[specO2] += -Wp[13];
  W[specO2plus] += +Wp[13];
  W[speceminus] += +Wp[13];

  /* reaction 6 */
//    {  O2  N2  O   N   O2+ N2+ O2- e-    
// MR=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 6  */
// MP=  {2,  0,  0,  0,  0,  0,  0,  1 },  /* 6  */
  Wp[12] = kf[12] * Nk[specO2] * Nk[specO2minus];
  W[specO2] += +Wp[12];
  W[specO2minus] += -Wp[12];
  W[speceminus] += +Wp[12];

  /* reaction 5b */
//    {/* O2  N2  O   N   O2+ N2+ O2- e-    */
// MR=  { 1,  1,  0,  0,  0,  0,  0,  1 },  /* 5b */
// MP=  { 0,  1,  0,  0,  0,  0,  1,  0 },  /* 5b */
  Wp[11] = kf[11] * Nk[speceminus] * Nk[specO2] * Nk[specN2];
  W[specO2] += -Wp[11];
  W[specO2minus] += +Wp[11];
  W[speceminus] += -Wp[11];

  /* reaction 5a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-    */
// MR=  {2,  0,  0,  0,  0,  0,  0,  1 },  /* 5a */
// MP=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 5a */
  Wp[10] = kf[10] * Nk[speceminus] * sqr ( Nk[specO2] );
  W[specO2] += -Wp[10];
  W[specO2minus] += +Wp[10];
  W[speceminus] += -Wp[10];

  /* reaction 4d */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {1,  0,  0,  0,  1,  0,  1,  0 },  /* 4d */
// MP=  {3,  0,  0,  0,  0,  0,  0,  0 },  /* 4d */
  Wp[9] = kf[9] * Nk[specO2] * Nk[specO2plus] * Nk[specO2minus];
  W[specO2] += +2.0 * Wp[9];
  W[specO2plus] += -Wp[9];
  W[specO2minus] += -Wp[9];

  /* reaction 4c */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {1,  0,  0,  0,  0,  1,  1,  0 },  /* 4c */
// MP=  {2,  1,  0,  0,  0,  0,  0,  0 },  /* 4c */
  Wp[8] = kf[8] * Nk[specO2] * Nk[specN2plus] * Nk[specO2minus];
  W[specO2] += +Wp[8];
  W[specN2] += +Wp[8];
  W[specN2plus] += -Wp[8];
  W[specO2minus] += -Wp[8];

  /* reaction 4b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  1,  0,  0,  1,  0,  1,  0 },  /* 4b */
// MP=  {2,  1,  0,  0,  0,  0,  0,  0 },  /* 4b */
  Wp[7] = kf[7] * Nk[specN2] * Nk[specO2plus] * Nk[specO2minus];
  W[specO2] += +2.0 * Wp[7];
  W[specO2plus] += -Wp[7];
  W[specO2minus] += -Wp[7];

  /* reaction 4a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  1,  0,  0,  0,  1,  1,  0 },  /* 4a */
// MP=  {1,  2,  0,  0,  0,  0,  0,  0 },  /* 4a */
  Wp[6] = kf[6] * Nk[specN2] * Nk[specN2plus] * Nk[specO2minus];
  W[specO2] += +Wp[6];
  W[specN2] += +Wp[6];
  W[specN2plus] += -Wp[6];
  W[specO2minus] += -Wp[6];

  /* reaction 3b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  1,  0,  1,  0 },  /* 3b */
// MP=  {2,  0,  0,  0,  0,  0,  0,  0 },  /* 3b */
  Wp[5] = kf[5] * Nk[specO2plus] * Nk[specO2minus];
  W[specO2] += +2.0 * Wp[5];
  W[specO2plus] += -Wp[5];
  W[specO2minus] += -Wp[5];

  /* reaction 3a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  0,  1,  1,  0 },  /* 3a */
// MP=  {1,  1,  0,  0,  0,  0,  0,  0 },  /* 3a */
  Wp[4] = kf[4] * Nk[specN2plus] * Nk[specO2minus];
  W[specO2] += +Wp[4];
  W[specN2] += +Wp[4];
  W[specN2plus] += -Wp[4];
  W[specO2minus] += -Wp[4];

  /* reaction 2b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 2b */
// MP=  {0,  0,  0,  2,  0,  0,  0,  0 },  /* 2b */
  Wp[3] = kf[3] * Nk[speceminus] * Nk[specN2plus];
  W[specN] += +2.0 * Wp[3];
  W[specN2plus] += -Wp[3];
  W[speceminus] += -Wp[3];

  /* reaction 2a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  1,  0,  0,  1 },  /* 2a */
// MP=  {0,  0,  2,  0,  0,  0,  0,  0 },  /* 2a */
  Wp[2] = kf[2] * Nk[speceminus] * Nk[specO2plus];
  W[specO] += +2.0 * Wp[2];
  W[specO2plus] += -Wp[2];
  W[speceminus] += -Wp[2];

  /* reaction 1b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {1,  0,  0,  0,  0,  0,  0,  1 },  /* 1b */
// MP=  {0,  0,  0,  0,  1,  0,  0,  2 },  /* 1b */
  Wp[1] = kf[1] * Nk[speceminus] * Nk[specO2];
  W[specO2] += -Wp[1];
  W[specO2plus] += +Wp[1];
  W[speceminus] += +Wp[1];

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

void find_dW_dx_Macheret2007old ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, r, s, spec;           /* counters */
  spec_t calM;
  spec2_t dWdNk;
  double Nk[ns];
  double dkdTe, dkdT, dkdQb;
  double N, dkdN;
  double kf[nr];
  double tmp, dtmpdT, sigma, theta1, theta2, ktownsend, dktownsenddNk;

  Estar = max ( Estarmin, Estar );

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
//  W[speceminus]+=+Wp[0];
//  W[specN2]+=-Wp[0];
//  W[specN2plus]+=+Wp[0];

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
                                                        46.0 * theta2 / notzero ( sigma,1e-99 ) *
                                                        pow ( log ( Estar ), 45.0 ) );
      dWdNk[specN2][spec] += -dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
      dWdNk[specN2plus][spec] += dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
      dWdNk[speceminus][spec] += dktownsenddNk * Nk[specN2] * Nk[speceminus] * 1e12;
    }

  }

  /* reaction 1b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {1,  0,  0,  0,  0,  0,  0,  1 },  /* 1b */
// MP=  {0,  0,  0,  0,  1,  0,  0,  2 },  /* 1b */
//  Wp[1]=kf[1]*Nk[speceminus]*Nk[specO2];
//  W[specO2]+=-Wp[1];
//  W[specO2plus]+=+Wp[1];
//  W[speceminus]+=+Wp[1];

  if ( TOWNSEND && TOWNSEND_IMPLICIT ) {
    kf[1] = _k1b ( Estar );
    dkdN = _dk1bdN ( Estar, Nk );

    dWdNk[specO2][speceminus] += -kf[1] * Nk[specO2];
    dWdNk[specO2][specO2] += -kf[1] * Nk[speceminus];

    dWdNk[specO2plus][speceminus] += kf[1] * Nk[specO2];
    dWdNk[specO2plus][specO2] += kf[1] * Nk[speceminus];

    dWdNk[speceminus][speceminus] += kf[1] * Nk[specO2];
    dWdNk[speceminus][specO2] += kf[1] * Nk[speceminus];

    for ( spec = 0; spec < ns; spec++ ) {
      dWdNk[specO2][spec] -= dkdN * Nk[speceminus] * Nk[specO2];
      dWdNk[specO2plus][spec] += dkdN * Nk[speceminus] * Nk[specO2];
      dWdNk[speceminus][spec] += dkdN * Nk[speceminus] * Nk[specO2];
    }

  }

  if ( TOWNSEND && TOWNSEND_SEMI_IMPLICIT ) {

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
      dWdNk[specO2][spec] += -dktownsenddNk * Nk[specO2] * Nk[speceminus] * 1e12;
      dWdNk[specO2plus][spec] += dktownsenddNk * Nk[specO2] * Nk[speceminus] * 1e12;
      dWdNk[speceminus][spec] += dktownsenddNk * Nk[specO2] * Nk[speceminus] * 1e12;
    }

  }

  /* reaction 2a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  1,  0,  0,  1 },  /* 2a */
// MP=  {0,  0,  2,  0,  0,  0,  0,  0 },  /* 2a */
  kf[2] = 2.0e-7 * pow ( 300.0 / Te, 0.7 );     /* 2a */

//  Wp[2]=kf[2]*Nk[speceminus]*Nk[specO2plus];
//  W[specO]+=+2.0*Wp[2];
//  W[specO2plus]+=-Wp[2];
//  W[speceminus]+=-Wp[2];

  dWdNk[specO][speceminus] += +2.0 * kf[2] * Nk[specO2plus];
  dWdNk[specO][specO2plus] += +2.0 * kf[2] * Nk[speceminus];

  dWdNk[specO2plus][speceminus] += -kf[2] * Nk[specO2plus];
  dWdNk[specO2plus][specO2plus] += -kf[2] * Nk[speceminus];

  dWdNk[speceminus][speceminus] += -kf[2] * Nk[specO2plus];
  dWdNk[speceminus][specO2plus] += -kf[2] * Nk[speceminus];

  dkdTe = -0.7 * 2.0e-7 * pow ( Te / 300.0, -1.7 ) / 300.0 * Nk[speceminus] * Nk[specO2plus];
  dWdTe[speceminus] += -dkdTe;
  dWdTe[specO] += +2.0 * dkdTe;
  dWdTe[specO2plus] += -dkdTe;

  /* reaction 2b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 2b */
// MP=  {0,  0,  0,  2,  0,  0,  0,  0 },  /* 2b */

  kf[3] = 2.8e-7 * pow ( 300.0 / Te, 0.5 );     /* 2b */

//  Wp[3]=kf[3]*Nk[speceminus]*Nk[specN2plus];
//  W[specN]+=+2.0*Wp[3];
//  W[specN2plus]+=-Wp[3];
//  W[speceminus]+=-Wp[3];

  dWdNk[specN][speceminus] += +2.0 * kf[3] * Nk[specN2plus];
  dWdNk[specN][specN2plus] += +2.0 * kf[3] * Nk[speceminus];
  dWdNk[specN2plus][speceminus] += -kf[3] * Nk[specN2plus];
  dWdNk[specN2plus][specN2plus] += -kf[3] * Nk[speceminus];
  dWdNk[speceminus][speceminus] += -kf[3] * Nk[specN2plus];
  dWdNk[speceminus][specN2plus] += -kf[3] * Nk[speceminus];

  dkdTe = -0.5 * 2.8e-7 * pow ( Te / 300.0, -1.5 ) / 300.0 * Nk[speceminus] * Nk[specN2plus];
  dWdTe[specN] += +2.0 * dkdTe;
  dWdTe[specN2plus] += -dkdTe;
  dWdTe[speceminus] += -dkdTe;

  /* reaction 3a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  0,  1,  1,  0 },  /* 3a */
// MP=  {1,  1,  0,  0,  0,  0,  0,  0 },  /* 3a */

  kf[4] = 2.0e-7 * pow ( 300.0 / T, 0.5 );      /* 3a */
//  Wp[4]=kf[4]*Nk[specN2plus]*Nk[specO2minus];
//  W[specO2]+=+Wp[4];
//  W[specN2]+=+Wp[4];
//  W[specN2plus]+=-Wp[4];
//  W[specO2minus]+=-Wp[4];

  dWdNk[specO2][specN2plus] += +kf[4] * Nk[specO2minus];
  dWdNk[specO2][specO2minus] += +kf[4] * Nk[specN2plus];
  dWdNk[specN2][specN2plus] += +kf[4] * Nk[specO2minus];
  dWdNk[specN2][specO2minus] += +kf[4] * Nk[specN2plus];
  dWdNk[specN2plus][specN2plus] += -kf[4] * Nk[specO2minus];
  dWdNk[specN2plus][specO2minus] += -kf[4] * Nk[specN2plus];
  dWdNk[specO2minus][specN2plus] += -kf[4] * Nk[specO2minus];
  dWdNk[specO2minus][specO2minus] += -kf[4] * Nk[specN2plus];
  dkdT = -0.5 * 2.0e-7 * pow ( T / 300.0, -1.5 ) / 300.0 * Nk[specN2plus] * Nk[specO2minus];
  dWdT[specO2] += +dkdT;
  dWdT[specN2] += +dkdT;
  dWdT[specN2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 3b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  0,  0,  0,  1,  0,  1,  0 },  /* 3b */
// MP=  {2,  0,  0,  0,  0,  0,  0,  0 },  /* 3b */

  kf[5] = kf[4];
//  Wp[5]=kf[5]*Nk[specO2plus]*Nk[specO2minus];
//  W[specO2]+=+2.0*Wp[5];
//  W[specO2plus]+=-Wp[5];
//  W[specO2minus]+=-Wp[5];
  dWdNk[specO2][specO2plus] += +2.0 * kf[5] * Nk[specO2minus];
  dWdNk[specO2][specO2minus] += +2.0 * kf[5] * Nk[specO2plus];
  dWdNk[specO2plus][specO2plus] += -kf[5] * Nk[specO2minus];
  dWdNk[specO2plus][specO2minus] += -kf[5] * Nk[specO2plus];
  dWdNk[specO2minus][specO2plus] += -kf[5] * Nk[specO2minus];
  dWdNk[specO2minus][specO2minus] += -kf[5] * Nk[specO2plus];
  dkdT = -0.5 * 2.0e-7 * pow ( T / 300.0, -1.5 ) / 300.0 * Nk[specO2plus] * Nk[specO2minus];
  dWdT[specO2] += +2.0 * dkdT;
  dWdT[specO2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  1,  0,  0,  0,  1,  1,  0 },  /* 4a */
// MP=  {1,  2,  0,  0,  0,  0,  0,  0 },  /* 4a */
  kf[6] = 2.0e-25 * pow ( 300.0 / T, 2.5 );     /* 4a */
//  Wp[6]=kf[6]*Nk[specN2]*Nk[specN2plus]*Nk[specO2minus];
//  W[specO2]+=+Wp[6];
//  W[specN2]+=+Wp[6];
//  W[specN2plus]+=-Wp[6];
//  W[specO2minus]+=-Wp[6];

  dWdNk[specO2][specN2] += +kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specO2][specN2plus] += +kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdNk[specO2][specO2minus] += +kf[6] * Nk[specN2] * Nk[specN2plus];
  dWdNk[specN2][specN2] += +kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specN2][specN2plus] += +kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdNk[specN2][specO2minus] += +kf[6] * Nk[specN2] * Nk[specN2plus];
  dWdNk[specN2plus][specN2] += -kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specN2plus][specN2plus] += -kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdNk[specN2plus][specO2minus] += -kf[6] * Nk[specN2] * Nk[specN2plus];
  dWdNk[specO2minus][specN2] += -kf[6] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specO2minus][specN2plus] += -kf[6] * Nk[specN2] * Nk[specO2minus];
  dWdNk[specO2minus][specO2minus] += -kf[6] * Nk[specN2] * Nk[specN2plus];
  dkdT = -2.5 * kf[6] / T * Nk[specN2] * Nk[specN2plus] * Nk[specO2minus];
  dWdT[specO2] += +dkdT;
  dWdT[specN2] += +dkdT;
  dWdT[specN2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4b */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {0,  1,  0,  0,  1,  0,  1,  0 },  /* 4b */
// MP=  {2,  1,  0,  0,  0,  0,  0,  0 },  /* 4b */
  kf[7] = kf[6];                /* 4b */
//  Wp[7]=kf[7]*Nk[specN2]*Nk[specO2plus]*Nk[specO2minus];
//  W[specO2]+=+2.0*Wp[7];
//  W[specO2plus]+=-Wp[7];
//  W[specO2minus]+=-Wp[7];
  dWdNk[specO2][specN2] += +2.0 * kf[7] * Nk[specO2plus] * Nk[specO2minus];
  dWdNk[specO2][specO2plus] += +2.0 * kf[7] * Nk[specN2] * Nk[specO2minus];
  dWdNk[specO2][specO2minus] += +2.0 * kf[7] * Nk[specN2] * Nk[specO2plus];
  dWdNk[specO2plus][specN2] += -kf[7] * Nk[specO2plus] * Nk[specO2minus];
  dWdNk[specO2plus][specO2plus] += -kf[7] * Nk[specN2] * Nk[specO2minus];
  dWdNk[specO2plus][specO2minus] += -kf[7] * Nk[specN2] * Nk[specO2plus];
  dWdNk[specO2minus][specN2] += -kf[7] * Nk[specO2plus] * Nk[specO2minus];
  dWdNk[specO2minus][specO2plus] += -kf[7] * Nk[specN2] * Nk[specO2minus];
  dWdNk[specO2minus][specO2minus] += -kf[7] * Nk[specN2] * Nk[specO2plus];
  dkdT = -2.5 * kf[7] / T * Nk[specN2] * Nk[specO2plus] * Nk[specO2minus];
  dWdT[specO2] += +2.0 * dkdT;  /* 2.0e-25*pow(300.0/T,2.5) */
  dWdT[specO2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4c */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {1,  0,  0,  0,  0,  1,  1,  0 },  /* 4c */
// MP=  {2,  1,  0,  0,  0,  0,  0,  0 },  /* 4c */
  kf[8] = kf[6];                /* 4c */
//  Wp[8]=kf[8]*Nk[specO2]*Nk[specN2plus]*Nk[specO2minus];
//  W[specO2]+=+Wp[8];
//  W[specN2]+=+Wp[8];
//  W[specN2plus]+=-Wp[8];
//  W[specO2minus]+=-Wp[8];

  dWdNk[specO2][specO2] += +kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specO2][specN2plus] += +kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdNk[specO2][specO2minus] += +kf[8] * Nk[specO2] * Nk[specN2plus];
  dWdNk[specN2][specO2] += +kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specN2][specN2plus] += +kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdNk[specN2][specO2minus] += +kf[8] * Nk[specO2] * Nk[specN2plus];
  dWdNk[specN2plus][specO2] += -kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specN2plus][specN2plus] += -kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdNk[specN2plus][specO2minus] += -kf[8] * Nk[specO2] * Nk[specN2plus];
  dWdNk[specO2minus][specO2] += -kf[8] * Nk[specN2plus] * Nk[specO2minus];
  dWdNk[specO2minus][specN2plus] += -kf[8] * Nk[specO2] * Nk[specO2minus];
  dWdNk[specO2minus][specO2minus] += -kf[8] * Nk[specO2] * Nk[specN2plus];
  dkdT = -2.5 * kf[8] / T * Nk[specO2] * Nk[specN2plus] * Nk[specO2minus];
  dWdT[specO2] += +dkdT;
  dWdT[specN2] += +dkdT;
  dWdT[specN2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 4d */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-     */
// MR=  {1,  0,  0,  0,  1,  0,  1,  0 },  /* 4d */
// MP=  {3,  0,  0,  0,  0,  0,  0,  0 },  /* 4d */
  kf[9] = kf[6];                /* 4d */
//  Wp[9]=kf[9]*Nk[specO2]*Nk[specO2plus]*Nk[specO2minus];
//  W[specO2]+=+2.0*Wp[9];
//  W[specO2plus]+=-Wp[9];
//  W[specO2minus]+=-Wp[9];

  dWdNk[specO2][specO2] += +2.0 * kf[9] * Nk[specO2plus] * Nk[specO2minus];
  dWdNk[specO2][specO2plus] += +2.0 * kf[9] * Nk[specO2] * Nk[specO2minus];
  dWdNk[specO2][specO2minus] += +2.0 * kf[9] * Nk[specO2] * Nk[specO2plus];
  dWdNk[specO2plus][specO2] += -kf[9] * Nk[specO2plus] * Nk[specO2minus];
  dWdNk[specO2plus][specO2plus] += -kf[9] * Nk[specO2] * Nk[specO2minus];
  dWdNk[specO2plus][specO2minus] += -kf[9] * Nk[specO2] * Nk[specO2plus];
  dWdNk[specO2minus][specO2] += -kf[9] * Nk[specO2plus] * Nk[specO2minus];
  dWdNk[specO2minus][specO2plus] += -kf[9] * Nk[specO2] * Nk[specO2minus];
  dWdNk[specO2minus][specO2minus] += -kf[9] * Nk[specO2] * Nk[specO2plus];
  dkdT = -2.5 * kf[9] / T * Nk[specO2] * Nk[specO2plus] * Nk[specO2minus];
  dWdT[specO2] += +2.0 * dkdT;  /* 2.0e-25*pow(300.0/T,2.5) */
  dWdT[specO2plus] += -dkdT;
  dWdT[specO2minus] += -dkdT;

  /* reaction 5a */
//    {/*O2  N2  O   N   O2+ N2+ O2- e-    */
// MR=  {2,  0,  0,  0,  0,  0,  0,  1 },  /* 5a */
// MP=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 5a */

  if ( ATTACHMENT ) {
    kf[10] = 1.4e-29 * 300.0 / Te * exp ( -600.0 / T ) * exp ( 700.0 * ( Te - T ) / Te / T );   /* 5a */
  } else {
    kf[10] = 0.0;
  }
//  Wp[10]=kf[10]*Nk[speceminus]*sqr(Nk[specO2]);
//  W[specO2]+=-Wp[10];
//  W[specO2minus]+=+Wp[10];
//  W[speceminus]+=-Wp[10];

  dWdNk[specO2][speceminus] += -kf[10] * sqr ( Nk[specO2] );
  dWdNk[specO2][specO2] += -kf[10] * Nk[speceminus] * 2.0 * Nk[specO2];
  dWdNk[specO2minus][speceminus] += +kf[10] * sqr ( Nk[specO2] );
  dWdNk[specO2minus][specO2] += +kf[10] * Nk[speceminus] * 2.0 * Nk[specO2];
  dWdNk[speceminus][speceminus] += -kf[10] * sqr ( Nk[specO2] );
  dWdNk[speceminus][specO2] += -kf[10] * Nk[speceminus] * 2.0 * Nk[specO2];

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
//    {/* O2  N2  O   N   O2+ N2+ O2- e-    */
// MR=  { 1,  1,  0,  0,  0,  0,  0,  1 },  /* 5b */
// MP=  { 0,  1,  0,  0,  0,  0,  1,  0 },  /* 5b */
  if ( ATTACHMENT ) {
    kf[11] = 1.07e-31 * sqr ( 300.0 / Te ) * exp ( -70.0 / T ) * exp ( 1500.0 * ( Te - T ) / Te / T );  /* 5b */
  } else {
    kf[11] = 0.0;
  }
//  Wp[11]=kf[11]*Nk[speceminus]*Nk[specO2]*Nk[specN2];
//  W[specO2]+=-Wp[11];
//  W[specO2minus]+=+Wp[11];
//  W[speceminus]+=-Wp[11];

  dWdNk[specO2][speceminus] += -kf[11] * Nk[specO2] * Nk[specN2];
  dWdNk[specO2][specO2] += -kf[11] * Nk[speceminus] * Nk[specN2];
  dWdNk[specO2][specN2] += -kf[11] * Nk[speceminus] * Nk[specO2];

  dWdNk[specO2minus][speceminus] += +kf[11] * Nk[specO2] * Nk[specN2];
  dWdNk[specO2minus][specO2] += +kf[11] * Nk[speceminus] * Nk[specN2];
  dWdNk[specO2minus][specN2] += +kf[11] * Nk[speceminus] * Nk[specO2];

  dWdNk[speceminus][speceminus] += -kf[11] * Nk[specO2] * Nk[specN2];
  dWdNk[speceminus][specO2] += -kf[11] * Nk[speceminus] * Nk[specN2];
  dWdNk[speceminus][specN2] += -kf[11] * Nk[speceminus] * Nk[specO2];

  dkdTe = ( -2.0 / Te + 1500.0 / sqr ( Te ) ) * kf[11] * Nk[speceminus] * Nk[specO2] * Nk[specN2];
  dkdT = ( 70.0 / sqr ( T ) - 1500.0 / sqr ( T ) ) * kf[11] * Nk[speceminus] * Nk[specO2] * Nk[specN2];
  dWdTe[specO2] += -dkdTe;
  dWdTe[specO2minus] += +dkdTe;
  dWdTe[speceminus] += -dkdTe;
  dWdT[specO2] += -dkdT;
  dWdT[specO2minus] += +dkdT;
  dWdT[speceminus] += -dkdT;

  /* reaction 6 */
//    {  O2  N2  O   N   O2+ N2+ O2- e-    
// MR=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 6  */
// MP=  {2,  0,  0,  0,  0,  0,  0,  1 },  /* 6  */
  kf[12] = 8.6e-10 * exp ( -6030.0 / T ) * ( 1.0 - exp ( -1570.0 / T ) );       /* 6 */
//  Wp[12]=kf[12]*Nk[specO2]*Nk[specO2minus];
//  W[specO2]+=+Wp[12];
//  W[specO2minus]+=-Wp[12];
//  W[speceminus]+=+Wp[12];
  dWdNk[specO2][specO2] += +kf[12] * Nk[specO2minus];
  dWdNk[specO2][specO2minus] += +kf[12] * Nk[specO2];
  dWdNk[specO2minus][specO2] += -kf[12] * Nk[specO2minus];
  dWdNk[specO2minus][specO2minus] += -kf[12] * Nk[specO2];
  dWdNk[speceminus][specO2] += +kf[12] * Nk[specO2minus];
  dWdNk[speceminus][specO2minus] += +kf[12] * Nk[specO2];

  dkdT =
    ( kf[12] * 6030.0 / sqr ( T ) - 8.6e-10 * exp ( -6030.0 / T ) * exp ( -1570.0 / T ) * 1570.0 / sqr ( T )
     ) * Nk[specO2] * Nk[specO2minus];

  dWdT[specO2] += +dkdT;
  dWdT[specO2minus] += -dkdT;
  dWdT[speceminus] += +dkdT;

  /* reaction 7a */
//        O2  N2  O   N   O2+ N2+ O2- e-    
// MR=  { 1,  0,  0,  0,  0,  0,  0,  0 },  /* 7a */
// MP=  { 0,  0,  0,  0,  1,  0,  0,  1 },  /* 7a */
  kf[13] = 2.0e11 * Qbeam / N;  /* 7a *//* Qbeam is in watts per m3, N in 1/cm3 */
//  Wp[13]=kf[13]*Nk[specO2];
//  W[specO2]+=-Wp[13];
//  W[specO2plus]+=+Wp[13];  
//  W[speceminus]+=+Wp[13];  

  dWdNk[specO2][specO2] += -kf[13];
  dWdNk[specO2plus][specO2] += +kf[13];
  dWdNk[speceminus][specO2] += +kf[13];
  for ( s = 0; s < ns; s++ ) {
    dWdNk[specO2][s] += +kf[13] * Nk[specO2] / N;
    dWdNk[specO2plus][s] += -kf[13] * Nk[specO2] / N;
    dWdNk[speceminus][s] += -kf[13] * Nk[specO2] / N;
  }
  dkdQb = 2.0e11 / N * Nk[specO2];
  dWdQbeam[specO2] += -dkdQb;
  dWdQbeam[specO2plus] += +dkdQb;
  dWdQbeam[speceminus] += +dkdQb;

  /* reaction 7b */
//       O2  N2  O   N   O2+ N2+ O2- e-  
// MR=  {0,  1,  0,  0,  0,  0,  0,  0 },  /* 7b */
// MP=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 7b */
  kf[14] = 1.8e11 * Qbeam / N;  /* 7b */
//  Wp[14]=kf[14]*Nk[specN2];
//  W[specN2]+=-Wp[14];
//  W[specN2plus]+=+Wp[14];  
//  W[speceminus]+=+Wp[14];  
  dWdNk[specN2][specN2] += -kf[14];
  dWdNk[specN2plus][specN2] += +kf[14];
  dWdNk[speceminus][specN2] += +kf[14];
  for ( s = 0; s < ns; s++ ) {
    dWdNk[specN2][s] += +kf[14] * Nk[specN2] / N;
    dWdNk[specN2plus][s] += -kf[14] * Nk[specN2] / N;
    dWdNk[speceminus][s] += -kf[14] * Nk[specN2] / N;
  }
  dkdQb = 1.8e11 / N * Nk[specN2];
  dWdQbeam[specN2] += -dkdQb;
  dWdQbeam[specN2plus] += +dkdQb;
  dWdQbeam[speceminus] += +dkdQb;

// reaction 8a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 2,  0,  0,  0,  0,  0,  0,  0},  /* 8a */
//MP= 1,  0,  2,  0,  0,  0,  0,  0},  /* 8a */
//  Wp[15]=kf[15]*Nk[specO2]*Nk[specO2];
//  W[specO2]+=-Wp[15];
//  W[specO]+=+2.0*Wp[15];

  tmp = exp ( -59380.0 / T ) * ( 1.0 - exp ( -2240.0 / T ) );
  dtmpdT = 59380.0 / sqr ( T ) * tmp - exp ( -59380.0 / T ) * 2240.0 / sqr ( T ) * exp ( -2240.0 / T );
  kf[15] = 3.7e-8 * tmp;        /* 8a */
  dWdNk[specO2][specO2] += -kf[15] * Nk[specO2];
  dWdNk[specO2][specO2] += -kf[15] * Nk[specO2];
  dWdNk[specO][specO2] += +2.0 * kf[15] * Nk[specO2];
  dWdNk[specO][specO2] += +2.0 * kf[15] * Nk[specO2];
  dWdT[specO2] += -dtmpdT * 3.7e-8 * Nk[specO2] * Nk[specO2];
  dWdT[specO] += +2.0 * dtmpdT * 3.7e-8 * Nk[specO2] * Nk[specO2];

// reaction 8b 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  1,  0,  0,  0,  0,  0,  0},  /* 8b */
//MP= 0,  1,  2,  0,  0,  0,  0,  0},  /* 8b */
//  Wp[16]=kf[16]*Nk[specO2]*Nk[specN2];
//  W[specO2]+=-Wp[16];
//  W[specO]+=+2.0*Wp[16];
//  tmp=exp(-59380.0/T)*(1.0-exp(-2240.0/T));  
  kf[16] = 9.3e-9 * tmp;        /* 8b */
  dWdNk[specO2][specO2] += -kf[16] * Nk[specN2];
  dWdNk[specO2][specN2] += -kf[16] * Nk[specO2];
  dWdNk[specO][specO2] += +2.0 * kf[16] * Nk[specN2];
  dWdNk[specO][specN2] += +2.0 * kf[16] * Nk[specO2];
  dWdT[specO2] += -dtmpdT * 9.3e-9 * Nk[specO2] * Nk[specN2];
  dWdT[specO] += +2.0 * dtmpdT * 9.3e-9 * Nk[specO2] * Nk[specN2];

// reaction 8c 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  1,  0,  0,  0,  0,  0},  /* 8c */
//MP= 0,  0,  3,  0,  0,  0,  0,  0},  /* 8c */
//  Wp[17]=kf[17]*Nk[specO2]*Nk[specO];
//  W[specO2]+=-Wp[17];
//  W[specO]+=+2.0*Wp[17];
//  tmp=exp(-59380.0/T)*(1.0-exp(-2240.0/T));  
  kf[17] = 1.3e-7 * tmp;        /* 8c */
  dWdNk[specO2][specO2] += -kf[17] * Nk[specO];
  dWdNk[specO2][specO] += -kf[17] * Nk[specO2];
  dWdNk[specO][specO2] += +2.0 * kf[17] * Nk[specO];
  dWdNk[specO][specO] += +2.0 * kf[17] * Nk[specO2];
  dWdT[specO2] += -dtmpdT * 1.3e-7 * Nk[specO2] * Nk[specO];
  dWdT[specO] += +2.0 * dtmpdT * 1.3e-7 * Nk[specO2] * Nk[specO];

// reaction 8d 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  1,  0,  0,  0,  0,  0,  0},  /* 8d */
//MP= 1,  0,  0,  2,  0,  0,  0,  0},  /* 8d */
//  Wp[18]=kf[18]*Nk[specO2]*Nk[specN2];
//  W[specN2]+=-Wp[18];
//  W[specN]+=+2.0*Wp[18];
  tmp = exp ( -113200.0 / T ) * ( 1.0 - exp ( -3354.0 / T ) );
  dtmpdT = 113200.0 / sqr ( T ) * tmp - exp ( -113200.0 / T ) * 3354.0 / sqr ( T ) * exp ( -3354.0 / T );
  kf[18] = 5.0e-8 * tmp;        /* 8d */
  dWdNk[specN2][specO2] += -kf[18] * Nk[specN2];
  dWdNk[specN2][specN2] += -kf[18] * Nk[specO2];
  dWdNk[specN][specO2] += +2.0 * kf[18] * Nk[specN2];
  dWdNk[specN][specN2] += +2.0 * kf[18] * Nk[specO2];
  dWdT[specN2] += -dtmpdT * 5.0e-8 * Nk[specO2] * Nk[specN2];
  dWdT[specN] += +2.0 * dtmpdT * 5.0e-8 * Nk[specO2] * Nk[specN2];

// reaction 8e 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  2,  0,  0,  0,  0,  0,  0},  /* 8e */
//MP= 0,  1,  0,  2,  0,  0,  0,  0},  /* 8e */
//  Wp[19]=kf[19]*Nk[specN2]*Nk[specN2];
//  W[specN2]+=-Wp[19];
//  W[specN]+=+2.0*Wp[19];
//  tmp=exp(-113200.0/T)*(1.0-exp(-3354.0/T));
  kf[19] = 5.0e-8 * tmp;        /* 8e */
  dWdNk[specN2][specN2] += -kf[19] * Nk[specN2];
  dWdNk[specN2][specN2] += -kf[19] * Nk[specN2];
  dWdNk[specN][specN2] += +2.0 * kf[19] * Nk[specN2];
  dWdNk[specN][specN2] += +2.0 * kf[19] * Nk[specN2];
  dWdT[specN2] += -dtmpdT * 5.0e-8 * Nk[specN2] * Nk[specN2];
  dWdT[specN] += +2.0 * dtmpdT * 5.0e-8 * Nk[specN2] * Nk[specN2];

// reaction 8f 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  1,  0,  0,  0,  0,  0},  /* 8f */
//MP= 0,  0,  1,  2,  0,  0,  0,  0},  /* 8f */
//  Wp[20]=kf[20]*Nk[specN2]*Nk[specO];
//  W[specN2]+=-Wp[20];
//  W[specN]+=+2.0*Wp[20];
//  tmp=exp(-113200.0/T)*(1.0-exp(-3354.0/T));
  kf[20] = 1.1e-7 * tmp;        /* 8f */
  dWdNk[specN2][specN2] += -kf[20] * Nk[specO];
  dWdNk[specN2][specO] += -kf[20] * Nk[specN2];
  dWdNk[specN][specN2] += +2.0 * kf[20] * Nk[specO];
  dWdNk[specN][specO] += +2.0 * kf[20] * Nk[specN2];
  dWdT[specN2] += -dtmpdT * 1.1e-7 * Nk[specN2] * Nk[specO];
  dWdT[specN] += +2.0 * dtmpdT * 1.1e-7 * Nk[specN2] * Nk[specO];

  if ( ATTACHMENT ) {
// reaction 9a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  2,  0,  0,  0,  0,  0},  /* 9a */
//MP= 2,  0,  0,  0,  0,  0,  0,  0},  /* 9a */
//  Wp[21]=kf[21]*Nk[specO2]*Nk[specO]*Nk[specO];
//  W[specO2]+=+Wp[21];
//  W[specO]+=-2.0*Wp[21];
    kf[21] = 2.45e-31 * pow ( T, -0.63 );       /* 9a */
    dtmpdT = -0.63 * pow ( T, -1.63 );
    dWdNk[specO2][specO2] += +kf[21] * Nk[specO] * Nk[specO];
    dWdNk[specO2][specO] += +kf[21] * Nk[specO2] * Nk[specO];
    dWdNk[specO2][specO] += +kf[21] * Nk[specO2] * Nk[specO];
    dWdNk[specO][specO2] += -2.0 * kf[21] * Nk[specO] * Nk[specO];
    dWdNk[specO][specO] += -2.0 * kf[21] * Nk[specO2] * Nk[specO];
    dWdNk[specO][specO] += -2.0 * kf[21] * Nk[specO2] * Nk[specO];
    dWdT[specO2] += +dtmpdT * 2.45e-31 * Nk[specO2] * Nk[specO] * Nk[specO];
    dWdT[specO] += -2.0 * dtmpdT * 2.45e-31 * Nk[specO2] * Nk[specO] * Nk[specO];

// reaction 9b 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  2,  0,  0,  0,  0,  0},  /* 9b */
//MP= 1,  1,  0,  0,  0,  0,  0,  0},  /* 9b */
//  Wp[22]=kf[22]*Nk[specN2]*Nk[specO]*Nk[specO];
//  W[specO2]+=+Wp[22];
//  W[specO]+=-2.0*Wp[22];
    kf[22] = 2.76e-34 * exp ( 720.0 / T );      /* 9b */
    dtmpdT = -720.0 / sqr ( T ) * exp ( 720.0 / T );
    dWdNk[specO2][specN2] += +kf[22] * Nk[specO] * Nk[specO];
    dWdNk[specO2][specO] += +kf[22] * Nk[specN2] * Nk[specO];
    dWdNk[specO2][specO] += +kf[22] * Nk[specN2] * Nk[specO];
    dWdNk[specO][specN2] += -2.0 * kf[22] * Nk[specO] * Nk[specO];
    dWdNk[specO][specO] += -2.0 * kf[22] * Nk[specN2] * Nk[specO];
    dWdNk[specO][specO] += -2.0 * kf[22] * Nk[specN2] * Nk[specO];
    dWdT[specO2] += +dtmpdT * 2.76e-34 * Nk[specN2] * Nk[specO] * Nk[specO];
    dWdT[specO] += -2.0 * dtmpdT * 2.76e-34 * Nk[specN2] * Nk[specO] * Nk[specO];

// reaction 9c 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  3,  0,  0,  0,  0,  0},  /* 9c */
//MP= 1,  0,  1,  0,  0,  0,  0,  0},  /* 9c */
//  Wp[23]=kf[23]*Nk[specO]*Nk[specO]*Nk[specO];
//  W[specO2]+=+Wp[23];
//  W[specO]+=-2.0*Wp[23];
    kf[23] = 8.8e-31 * pow ( T, -0.63 );        /* 9c */
    dtmpdT = -0.63 * pow ( T, -1.63 );
    dWdNk[specO2][specO] += +kf[23] * Nk[specO] * Nk[specO];
    dWdNk[specO2][specO] += +kf[23] * Nk[specO] * Nk[specO];
    dWdNk[specO2][specO] += +kf[23] * Nk[specO] * Nk[specO];
    dWdNk[specO][specO] += -2.0 * kf[23] * Nk[specO] * Nk[specO];
    dWdNk[specO][specO] += -2.0 * kf[23] * Nk[specO] * Nk[specO];
    dWdNk[specO][specO] += -2.0 * kf[23] * Nk[specO] * Nk[specO];
    dWdT[specO2] += +dtmpdT * 8.8e-31 * Nk[specO] * Nk[specO] * Nk[specO];
    dWdT[specO] += -2.0 * dtmpdT * 8.8e-31 * Nk[specO] * Nk[specO] * Nk[specO];

// reaction 9d 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  0,  2,  0,  0,  0,  0},  /* 9d */
//MP= 1,  1,  0,  0,  0,  0,  0,  0},  /* 9d */
//  Wp[24]=kf[24]*Nk[specO2]*Nk[specN]*Nk[specN];
//  W[specN2]+=+Wp[24];
//  W[specN]+=-2.0*Wp[24];
    tmp = 8.27e-34 * exp ( 500.0 / T );
    dtmpdT = -8.27e-34 * 500.0 / sqr ( T ) * exp ( 500.0 / T );
    kf[24] = tmp;               /* 9d */
    dWdNk[specN2][specO2] += +kf[24] * Nk[specN] * Nk[specN];
    dWdNk[specN2][specN] += +kf[24] * Nk[specO2] * Nk[specN];
    dWdNk[specN2][specN] += +kf[24] * Nk[specO2] * Nk[specN];
    dWdNk[specN][specO2] += -2.0 * kf[24] * Nk[specN] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[24] * Nk[specO2] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[24] * Nk[specO2] * Nk[specN];
    dWdT[specN2] += +dtmpdT * Nk[specO2] * Nk[specN] * Nk[specN];
    dWdT[specN] += -2.0 * dtmpdT * Nk[specO2] * Nk[specN] * Nk[specN];

// reaction 9e 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  0,  2,  0,  0,  0,  0},  /* 9e */
//MP= 0,  2,  0,  0,  0,  0,  0,  0},  /* 9e */
//  Wp[25]=kf[25]*Nk[specN2]*Nk[specN]*Nk[specN];
//  W[specN2]+=+Wp[25];
//  W[specN]+=-2.0*Wp[25];
//  tmp=8.27e-34*exp(500.0/T);
    kf[25] = tmp;               /* 9e */
    dWdNk[specN2][specN2] += +kf[25] * Nk[specN] * Nk[specN];
    dWdNk[specN2][specN] += +kf[25] * Nk[specN2] * Nk[specN];
    dWdNk[specN2][specN] += +kf[25] * Nk[specN2] * Nk[specN];
    dWdNk[specN][specN2] += -2.0 * kf[25] * Nk[specN] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[25] * Nk[specN2] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[25] * Nk[specN2] * Nk[specN];
    dWdT[specN2] += +dtmpdT * Nk[specN2] * Nk[specN] * Nk[specN];
    dWdT[specN] += -2.0 * dtmpdT * Nk[specN2] * Nk[specN] * Nk[specN];

// reaction 9f 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  1,  2,  0,  0,  0,  0},  /* 9f */
//MP= 0,  1,  1,  0,  0,  0,  0,  0},  /* 9f */
//  Wp[26]=kf[26]*Nk[specO]*Nk[specN]*Nk[specN];
//  W[specN2]+=+Wp[26];
//  W[specN]+=-2.0*Wp[26];
//  tmp=8.27e-34*exp(500.0/T);
    kf[26] = tmp;               /* 9f */
    dWdNk[specN2][specO] += +kf[26] * Nk[specN] * Nk[specN];
    dWdNk[specN2][specN] += +kf[26] * Nk[specO] * Nk[specN];
    dWdNk[specN2][specN] += +kf[26] * Nk[specO] * Nk[specN];
    dWdNk[specN][specO] += -2.0 * kf[26] * Nk[specN] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[26] * Nk[specO] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[26] * Nk[specO] * Nk[specN];
    dWdT[specN2] += +dtmpdT * Nk[specO] * Nk[specN] * Nk[specN];
    dWdT[specN] += -2.0 * dtmpdT * Nk[specO] * Nk[specN] * Nk[specN];

// reaction 9g 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  0,  3,  0,  0,  0,  0},  /* 9g */
//MP= 0,  1,  0,  1,  0,  0,  0,  0},  /* 9g */
//  Wp[27]=kf[27]*Nk[specN]*Nk[specN]*Nk[specN];
//  W[specN2]+=+Wp[27];
//  W[specN]+=-2.0*Wp[27];
//  tmp=8.27e-34*exp(500.0/T);
    kf[27] = tmp;               /* 9g */
    dWdNk[specN2][specN] += +kf[27] * Nk[specN] * Nk[specN];
    dWdNk[specN2][specN] += +kf[27] * Nk[specN] * Nk[specN];
    dWdNk[specN2][specN] += +kf[27] * Nk[specN] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[27] * Nk[specN] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[27] * Nk[specN] * Nk[specN];
    dWdNk[specN][specN] += -2.0 * kf[27] * Nk[specN] * Nk[specN];
    dWdT[specN2] += +dtmpdT * Nk[specN] * Nk[specN] * Nk[specN];
    dWdT[specN] += -2.0 * dtmpdT * Nk[specN] * Nk[specN] * Nk[specN];
  }

  for ( s = 0; s < ns; s++ ) {
    for ( r = 0; r < ns; r++ ) {
      dWdrhok[s][r] = dWdNk[s][r] * calM[s] / calM[r];
    }
    dWdT[s] = dWdT[s] / calA * calM[s] * 1.0e6;
    dWdTe[s] = dWdTe[s] / calA * calM[s] * 1.0e6;
    dWdQbeam[s] = dWdQbeam[s] / calA * calM[s] * 1.0e6;
  }

}




void find_Qei_Macheret2007old(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  double theta;
  
  *Qei=0.0;  
  theta=log(Estar);

  if (TOWNSEND)  {
    add_to_Qei(specN2, exp(-0.0105809*sqr(theta)-2.40411e-75*pow(theta,46.0)), rhok, Qei);
    add_to_Qei(specO2, exp(-0.0102785*sqr(theta)-2.42260e-75*pow(theta,46.0)), rhok, Qei);
  }
}



void find_dQei_dx_Macheret2007old(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  double theta;
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
  theta=log(Estar);

  if (TOWNSEND)  {
    add_to_dQei(specN2, exp(-0.0105809*sqr(theta)-2.40411e-75*pow(theta,46.0)), 0.0,  rhok, dQeidrhok, dQeidTe);
    add_to_dQei(specO2, exp(-0.0102785*sqr(theta)-2.42260e-75*pow(theta,46.0)), 0.0,  rhok, dQeidrhok, dQeidTe);
  }
}


