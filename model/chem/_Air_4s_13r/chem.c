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

#define nr 13

void find_W ( spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long spec;
  double kf[nr];
  double Wp[nr];
  double Nk[ns];
  spec_t calM;
  double N, tmp;

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

  tmp = exp ( -59380.0 / T ) * ( 1.0 - exp ( -2240.0 / T ) );
  kf[0] = 3.7e-8 * tmp;         /* 8a */
  kf[1] = 9.3e-9 * tmp;         /* 8b */
  kf[2] = 1.3e-7 * tmp;         /* 8c */
  tmp = exp ( -113200.0 / T ) * ( 1.0 - exp ( -3354.0 / T ) );
  kf[3] = 5.0e-8 * tmp;         /* 8d */
  kf[4] = 5.0e-8 * tmp;         /* 8e */
  kf[5] = 1.1e-7 * tmp;         /* 8f */
  kf[6] = 2.45e-31 * pow ( T, -0.63 );  /* 9a */
  kf[7] = 2.76e-34 * exp ( 720.0 / T ); /* 9b */
  kf[8] = 8.8e-31 * pow ( T, -0.63 );   /* 9c */
  tmp = 8.27e-34 * exp ( 500.0 / T );
  kf[9] = tmp;                  /* 9d */
  kf[10] = tmp;                 /* 9e */
  kf[11] = tmp;                 /* 9f */
  kf[12] = tmp;                 /* 9g */

// reaction 9g 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  0,  3,  0,  0,  0,  0},  /* 9g */
//MP= 0,  1,  0,  1,  0,  0,  0,  0},  /* 9g */
  Wp[12] = kf[12] * Nk[specN] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[12];
  W[specN] += -2.0 * Wp[12];

// reaction 9f 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  1,  2,  0,  0,  0,  0},  /* 9f */
//MP= 0,  1,  1,  0,  0,  0,  0,  0},  /* 9f */
  Wp[11] = kf[11] * Nk[specO] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[11];
  W[specN] += -2.0 * Wp[11];

// reaction 9e 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  0,  2,  0,  0,  0,  0},  /* 9e */
//MP= 0,  2,  0,  0,  0,  0,  0,  0},  /* 9e */
  Wp[10] = kf[10] * Nk[specN2] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[10];
  W[specN] += -2.0 * Wp[10];

// reaction 9d 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  0,  2,  0,  0,  0,  0},  /* 9d */
//MP= 1,  1,  0,  0,  0,  0,  0,  0},  /* 9d */
  Wp[9] = kf[9] * Nk[specO2] * Nk[specN] * Nk[specN];
  W[specN2] += +Wp[9];
  W[specN] += -2.0 * Wp[9];

// reaction 9c 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  3,  0,  0,  0,  0,  0},  /* 9c */
//MP= 1,  0,  1,  0,  0,  0,  0,  0},  /* 9c */
  Wp[8] = kf[8] * Nk[specO] * Nk[specO] * Nk[specO];
  W[specO2] += +Wp[8];
  W[specO] += -2.0 * Wp[8];

// reaction 9b 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  2,  0,  0,  0,  0,  0},  /* 9b */
//MP= 1,  1,  0,  0,  0,  0,  0,  0},  /* 9b */
  Wp[7] = kf[7] * Nk[specN2] * Nk[specO] * Nk[specO];
  W[specO2] += +Wp[7];
  W[specO] += -2.0 * Wp[7];

// reaction 9a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  2,  0,  0,  0,  0,  0},  /* 9a */
//MP= 2,  0,  0,  0,  0,  0,  0,  0},  /* 9a */
  Wp[6] = kf[6] * Nk[specO2] * Nk[specO] * Nk[specO];
  W[specO2] += +Wp[6];
  W[specO] += -2.0 * Wp[6];

// reaction 8f 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  1,  0,  0,  0,  0,  0},  /* 8f */
//MP= 0,  0,  1,  2,  0,  0,  0,  0},  /* 8f */
  Wp[5] = kf[5] * Nk[specN2] * Nk[specO];
  W[specN2] += -Wp[5];
  W[specN] += +2.0 * Wp[5];

// reaction 8e 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  2,  0,  0,  0,  0,  0,  0},  /* 8e */
//MP= 0,  1,  0,  2,  0,  0,  0,  0},  /* 8e */
  Wp[4] = kf[4] * Nk[specN2] * Nk[specN2];
  W[specN2] += -Wp[4];
  W[specN] += +2.0 * Wp[4];

// reaction 8d 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  1,  0,  0,  0,  0,  0,  0},  /* 8d */
//MP= 1,  0,  0,  2,  0,  0,  0,  0},  /* 8d */
  Wp[3] = kf[3] * Nk[specO2] * Nk[specN2];
  W[specN2] += -Wp[3];
  W[specN] += +2.0 * Wp[3];

// reaction 8c 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  1,  0,  0,  0,  0,  0},  /* 8c */
//MP= 0,  0,  3,  0,  0,  0,  0,  0},  /* 8c */
  Wp[2] = kf[2] * Nk[specO2] * Nk[specO];
  W[specO2] += -Wp[2];
  W[specO] += +2.0 * Wp[2];

// reaction 8b 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  1,  0,  0,  0,  0,  0,  0},  /* 8b */
//MP= 0,  1,  2,  0,  0,  0,  0,  0},  /* 8b */
  Wp[1] = kf[1] * Nk[specO2] * Nk[specN2];
  W[specO2] += -Wp[1];
  W[specO] += +2.0 * Wp[1];

// reaction 8a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 2,  0,  0,  0,  0,  0,  0,  0},  /* 8a */
//MP= 1,  0,  2,  0,  0,  0,  0,  0},  /* 8a */
  Wp[0] = kf[0] * Nk[specO2] * Nk[specO2];
  W[specO2] += -Wp[0];
  W[specO] += +2.0 * Wp[0];

  for ( spec = 0; spec < ns; spec++ )
    W[spec] = W[spec] / calA * calM[spec] * 1.0e6;
}

void find_dW_dx ( spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, r, s, spec;           /* counters */
  spec_t calM;
  double Nk[ns];
  double N;
  double kf[nr];
  double tmp, dtmpdT;

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

// reaction 8a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 2,  0,  0,  0,  0,  0,  0,  0},  /* 8a */
//MP= 1,  0,  2,  0,  0,  0,  0,  0},  /* 8a */
//  Wp[15]=kf[15]*Nk[specO2]*Nk[specO2];
//  W[specO2]+=-Wp[15];
//  W[specO]+=+2.0*Wp[15];

  tmp = exp ( -59380.0 / T ) * ( 1.0 - exp ( -2240.0 / T ) );
  dtmpdT = 59380.0 / sqr ( T ) * tmp - exp ( -59380.0 / T ) * 2240.0 / sqr ( T ) * exp ( -2240.0 / T );
  kf[0] = 3.7e-8 * tmp;         /* 8a */
  dWdrhok[specO2][specO2] += -kf[0] * Nk[specO2];
  dWdrhok[specO2][specO2] += -kf[0] * Nk[specO2];
  dWdrhok[specO][specO2] += +2.0 * kf[0] * Nk[specO2];
  dWdrhok[specO][specO2] += +2.0 * kf[0] * Nk[specO2];
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
  kf[1] = 9.3e-9 * tmp;         /* 8b */
  dWdrhok[specO2][specO2] += -kf[1] * Nk[specN2];
  dWdrhok[specO2][specN2] += -kf[1] * Nk[specO2];
  dWdrhok[specO][specO2] += +2.0 * kf[1] * Nk[specN2];
  dWdrhok[specO][specN2] += +2.0 * kf[1] * Nk[specO2];
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
  kf[2] = 1.3e-7 * tmp;         /* 8c */
  dWdrhok[specO2][specO2] += -kf[2] * Nk[specO];
  dWdrhok[specO2][specO] += -kf[2] * Nk[specO2];
  dWdrhok[specO][specO2] += +2.0 * kf[2] * Nk[specO];
  dWdrhok[specO][specO] += +2.0 * kf[2] * Nk[specO2];
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
  kf[3] = 5.0e-8 * tmp;         /* 8d */
  dWdrhok[specN2][specO2] += -kf[3] * Nk[specN2];
  dWdrhok[specN2][specN2] += -kf[3] * Nk[specO2];
  dWdrhok[specN][specO2] += +2.0 * kf[3] * Nk[specN2];
  dWdrhok[specN][specN2] += +2.0 * kf[3] * Nk[specO2];
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
  kf[4] = 5.0e-8 * tmp;         /* 8e */
  dWdrhok[specN2][specN2] += -kf[4] * Nk[specN2];
  dWdrhok[specN2][specN2] += -kf[4] * Nk[specN2];
  dWdrhok[specN][specN2] += +2.0 * kf[4] * Nk[specN2];
  dWdrhok[specN][specN2] += +2.0 * kf[4] * Nk[specN2];
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
  kf[5] = 1.1e-7 * tmp;         /* 8f */
  dWdrhok[specN2][specN2] += -kf[5] * Nk[specO];
  dWdrhok[specN2][specO] += -kf[5] * Nk[specN2];
  dWdrhok[specN][specN2] += +2.0 * kf[5] * Nk[specO];
  dWdrhok[specN][specO] += +2.0 * kf[5] * Nk[specN2];
  dWdT[specN2] += -dtmpdT * 1.1e-7 * Nk[specN2] * Nk[specO];
  dWdT[specN] += +2.0 * dtmpdT * 1.1e-7 * Nk[specN2] * Nk[specO];

// reaction 9a 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 1,  0,  2,  0,  0,  0,  0,  0},  /* 9a */
//MP= 2,  0,  0,  0,  0,  0,  0,  0},  /* 9a */
//  Wp[21]=kf[21]*Nk[specO2]*Nk[specO]*Nk[specO];
//  W[specO2]+=+Wp[21];
//  W[specO]+=-2.0*Wp[21];
  kf[6] = 2.45e-31 * pow ( T, -0.63 );  /* 9a */
  dtmpdT = -0.63 * pow ( T, -1.63 );
  dWdrhok[specO2][specO2] += +kf[6] * Nk[specO] * Nk[specO];
  dWdrhok[specO2][specO] += +kf[6] * Nk[specO2] * Nk[specO];
  dWdrhok[specO2][specO] += +kf[6] * Nk[specO2] * Nk[specO];
  dWdrhok[specO][specO2] += -2.0 * kf[6] * Nk[specO] * Nk[specO];
  dWdrhok[specO][specO] += -2.0 * kf[6] * Nk[specO2] * Nk[specO];
  dWdrhok[specO][specO] += -2.0 * kf[6] * Nk[specO2] * Nk[specO];
  dWdT[specO2] += +dtmpdT * 2.45e-31 * Nk[specO2] * Nk[specO] * Nk[specO];
  dWdT[specO] += -2.0 * dtmpdT * 2.45e-31 * Nk[specO2] * Nk[specO] * Nk[specO];

// reaction 9b 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  1,  2,  0,  0,  0,  0,  0},  /* 9b */
//MP= 1,  1,  0,  0,  0,  0,  0,  0},  /* 9b */
//  Wp[22]=kf[22]*Nk[specN2]*Nk[specO]*Nk[specO];
//  W[specO2]+=+Wp[22];
//  W[specO]+=-2.0*Wp[22];
  kf[7] = 2.76e-34 * exp ( 720.0 / T ); /* 9b */
  dtmpdT = -720.0 / sqr ( T ) * exp ( 720.0 / T );
  dWdrhok[specO2][specN2] += +kf[7] * Nk[specO] * Nk[specO];
  dWdrhok[specO2][specO] += +kf[7] * Nk[specN2] * Nk[specO];
  dWdrhok[specO2][specO] += +kf[7] * Nk[specN2] * Nk[specO];
  dWdrhok[specO][specN2] += -2.0 * kf[7] * Nk[specO] * Nk[specO];
  dWdrhok[specO][specO] += -2.0 * kf[7] * Nk[specN2] * Nk[specO];
  dWdrhok[specO][specO] += -2.0 * kf[7] * Nk[specN2] * Nk[specO];
  dWdT[specO2] += +dtmpdT * 2.76e-34 * Nk[specN2] * Nk[specO] * Nk[specO];
  dWdT[specO] += -2.0 * dtmpdT * 2.76e-34 * Nk[specN2] * Nk[specO] * Nk[specO];

// reaction 9c 
//    O2  N2  O   N   O2+ N2+ O2- e-   
//MR= 0,  0,  3,  0,  0,  0,  0,  0},  /* 9c */
//MP= 1,  0,  1,  0,  0,  0,  0,  0},  /* 9c */
//  Wp[23]=kf[23]*Nk[specO]*Nk[specO]*Nk[specO];
//  W[specO2]+=+Wp[23];
//  W[specO]+=-2.0*Wp[23];
  kf[8] = 8.8e-31 * pow ( T, -0.63 );   /* 9c */
  dtmpdT = -0.63 * pow ( T, -1.63 );
  dWdrhok[specO2][specO] += +kf[8] * Nk[specO] * Nk[specO];
  dWdrhok[specO2][specO] += +kf[8] * Nk[specO] * Nk[specO];
  dWdrhok[specO2][specO] += +kf[8] * Nk[specO] * Nk[specO];
  dWdrhok[specO][specO] += -2.0 * kf[8] * Nk[specO] * Nk[specO];
  dWdrhok[specO][specO] += -2.0 * kf[8] * Nk[specO] * Nk[specO];
  dWdrhok[specO][specO] += -2.0 * kf[8] * Nk[specO] * Nk[specO];
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
  kf[9] = tmp;                  /* 9d */
  dWdrhok[specN2][specO2] += +kf[9] * Nk[specN] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[9] * Nk[specO2] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[9] * Nk[specO2] * Nk[specN];
  dWdrhok[specN][specO2] += -2.0 * kf[9] * Nk[specN] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[9] * Nk[specO2] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[9] * Nk[specO2] * Nk[specN];
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
  kf[10] = tmp;                 /* 9e */
  dWdrhok[specN2][specN2] += +kf[10] * Nk[specN] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[10] * Nk[specN2] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[10] * Nk[specN2] * Nk[specN];
  dWdrhok[specN][specN2] += -2.0 * kf[10] * Nk[specN] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[10] * Nk[specN2] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[10] * Nk[specN2] * Nk[specN];
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
  kf[11] = tmp;                 /* 9f */
  dWdrhok[specN2][specO] += +kf[11] * Nk[specN] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[11] * Nk[specO] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[11] * Nk[specO] * Nk[specN];
  dWdrhok[specN][specO] += -2.0 * kf[11] * Nk[specN] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[11] * Nk[specO] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[11] * Nk[specO] * Nk[specN];
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
  kf[12] = tmp;                 /* 9g */
  dWdrhok[specN2][specN] += +kf[12] * Nk[specN] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[12] * Nk[specN] * Nk[specN];
  dWdrhok[specN2][specN] += +kf[12] * Nk[specN] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[12] * Nk[specN] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[12] * Nk[specN] * Nk[specN];
  dWdrhok[specN][specN] += -2.0 * kf[12] * Nk[specN] * Nk[specN];
  dWdT[specN2] += +dtmpdT * Nk[specN] * Nk[specN] * Nk[specN];
  dWdT[specN] += -2.0 * dtmpdT * Nk[specN] * Nk[specN] * Nk[specN];

  for ( s = 0; s < ns; s++ ) {
    for ( r = 0; r < ns; r++ ) {
      dWdrhok[s][r] = dWdrhok[s][r] * calM[s] / calM[r];
    }
    dWdT[s] = dWdT[s] / calA * calM[s] * 1.0e6;
    dWdTe[s] = dWdTe[s] / calA * calM[s] * 1.0e6;
    dWdQbeam[s] = dWdQbeam[s] / calA * calM[s] * 1.0e6;
  }

}





void find_Qei(spec_t rhok, double Estar, double Te, double *Qei){  
  *Qei=0.0;
}


void find_dQei_dx(spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
}

