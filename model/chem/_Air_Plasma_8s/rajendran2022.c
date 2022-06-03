// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2021 Prasanna Thoguluva Rajendran

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

#include <model/thermo/_thermo.h>
#include <model/share/chem_share.h>



/* set all reactions to true except for testing purposes */
const static bool REACTION[23]=
  {
   TRUE, /* reaction 0 */
   TRUE, /* reaction 1 */
   TRUE, /* reaction 2 */
   TRUE, /* reaction 3 */
   TRUE, /* reaction 4 */
   TRUE, /* reaction 5 */
   TRUE, /* reaction 6 */
   TRUE, /* reaction 7 */
   TRUE, /* reaction 8 */
   TRUE, /* reaction 9 */
   TRUE, /* reaction 10 */
   TRUE, /* reaction 11 */
   TRUE, /* reaction 12 */
   TRUE, /* reaction 13 */
   TRUE, /* reaction 14 */
   TRUE, /* reaction 15 */
   TRUE, /* reaction 16 */
   TRUE, /* reaction 17 */
   TRUE, /* reaction 18 */
   TRUE, /* reaction 19 */ 
   TRUE, /* reaction 20 */
   TRUE, /* reaction 21 */
   TRUE, /* reaction 22 */
  };

#define specEND -1

const static long specM1[]=
  {
   specN2, specO2, specEND
  };

const static long specM2[]=
  {
   specN2, specO2, specN, specO, specEND
  };



void find_W_Rajendran2022 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                       double Qbeam, spec_t W ) {
  double tmp,Nt,R;
  long k;
  spec_t N,X;

  Nt = 0.0;
  R=1.9872;
  for ( k = 0; k < ns; k++ ) {
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    Nt += N[k];
  }
  

  if (REACTION[1]){
    add_to_W_fw_2r3p ( specN2, speceminus,   specN2plus, speceminus, speceminus, 2.6291e-8 * calA, 0.142, 536330.0*R, Te, X, W);
  }
  
  if (REACTION[2]) {
    add_to_W_fw_2r3p ( specO2, speceminus,   specO2plus, speceminus, speceminus,   3.869e-8 * calA, 0.1166, 567360.0*R, Te, X, W);
  }

  if (REACTION[3]){
    add_to_W_2r2p ( specO2plus, speceminus,   specO, specO,  2.0e-7 * pow (300.0/Te, 0.7), N, W);
  }

  if (REACTION[4]){
    add_to_W_2r2p ( specN2plus, speceminus,   specN, specN,  2.8e-7 * pow (300.0/Te, 0.5), N, W);
  }
  
  tmp = 2.0e-7 * pow (300.0/T, 0.5);
  if (REACTION[5]){
    add_to_W_2r2p ( specO2minus, specN2plus,   specO2, specN2,  tmp, N, W);
  }

  if (REACTION[6]){
    add_to_W_2r2p ( specO2minus, specO2plus,   specO2, specO2,  tmp, N, W);
  }
  
  tmp = 2.0e-25 * pow (300.0/T, 2.5);
  if (REACTION[7]){
    for (k=0; specM1[k]!=specEND; k++){
      add_to_W_3r3p ( specO2minus, specN2plus, specM1[k],   specO2, specN2, specM1[k], tmp, N, W);
    }
  }

  if (REACTION[8]){
    for (k=0; specM1[k]!=specEND; k++){
      add_to_W_3r3p ( specO2minus, specO2plus, specM1[k],   specO2, specO2, specM1[k], tmp, N, W);
    }
  }
  
  if (REACTION[9]){
    add_to_W_3r2p ( speceminus, specO2, specO2,   specO2minus, specO2,  1.4e-29 * 300.0 / Te * exp ( -600.0 / T ) * exp ( 700.0 * ( Te - T ) / Te / T ), N, W);
  }
  
  if (REACTION[10]){
    add_to_W_3r2p ( speceminus, specO2, specN2,   specO2minus, specN2,  1.07e-31 * sqr ( 300.0 / Te ) * exp ( -70.0 / T ) * exp ( 1500.0 * ( Te - T ) / Te / T ), N, W);
  }
  
  if (REACTION[11]){
    add_to_W_2r3p ( specO2minus, specO2,   speceminus, specO2, specO2,  8.6e-10 * exp ( -6030.0 / T ) * ( 1.0 - exp ( -1570.0 / T ) ), N, W);
  }

  if (REACTION[12]){
    tmp = 2.0e11 * Qbeam / Nt * N[specO2]; /* Qbeam is in watts per m3, Nt in 1/cm3 */
    W[specO2] += -tmp * _calM(specO2) / calA * 1.0e6;
    W[specO2plus] += +tmp * _calM(specO2plus) / calA * 1.0e6;
    W[speceminus] += +tmp * _calM(speceminus) / calA * 1.0e6;
  }
  
  if (REACTION[13]){
    tmp = 1.8e11 * Qbeam / Nt * N[specN2];
    W[specN2] += -tmp * _calM(specN2) / calA * 1.0e6;
    W[specN2plus] += +tmp * _calM(specN2plus) / calA * 1.0e6;
    W[speceminus] += +tmp * _calM(speceminus) / calA * 1.0e6;
  }
  
  tmp = exp ( -59380.0 / T ) * ( 1.0 - exp ( -2240.0 / T ) );
  if (REACTION[14]){
    add_to_W_2r3p ( specO2, specO2,   specO, specO, specO2,  3.7e-8 * tmp, N, W);
  }

  if (REACTION[15]){
    add_to_W_2r3p ( specO2, specN2,   specO, specO, specN2,  9.3e-9 * tmp, N, W);
  }

  if (REACTION[16]){
    add_to_W_2r3p ( specO2, specO,   specO, specO, specO,  1.3e-7 * tmp, N, W);
  }

  tmp = exp ( -113200.0 / T ) * ( 1.0 - exp ( -3354.0 / T ) );
  if (REACTION[17]){
    for (k=0; specM1[k]!=specEND; k++){
      add_to_W_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 5.0e-8 * tmp, N, W);
    }
  }

  if (REACTION[18]){
    add_to_W_2r3p ( specN2, specO,   specN, specN, specO, 1.1e-7 * tmp, N, W);
  }
  
  tmp = pow ( T, -0.63 );
  if (REACTION[19]){
    add_to_W_3r2p ( specO, specO, specO2,   specO2, specO2, 2.45e-31 * tmp, N, W);
  }

  if (REACTION[20]){
    add_to_W_3r2p ( specO, specO, specN2,   specO2, specN2, 2.76e-34 * exp ( 720.0 / T ), N, W);
  }

  if (REACTION[21]){
    add_to_W_3r2p ( specO, specO, specO,   specO2, specO, 8.8e-31 * tmp, N, W);
  }

  tmp = 8.27e-34 * exp ( 500.0 / T );
  if (REACTION[22]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_W_3r2p ( specN, specN, specM2[k],   specN2, specM2[k], tmp, N, W);
    }
  }

}



/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx_Rajendran2022 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    
  spec_t N,X;
  double Nt,R,kf,tmp,dkfdTe,dkfdT,dkfdTv,dtmpdT,dkfdQb;

  /* initialize all derivatives to zero */
  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }
  
  /* find properties needed by add_to_dW* functions in proper units */
  Nt = 0.0;
  R=1.9872;
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA; /* particules/cm^3 */
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    Nt += N[k];
  }


  if (REACTION[1]){
    add_to_dW_fw_2r3p ( specN2, speceminus,   specN2plus, speceminus, speceminus, 2.6291e-8 * calA, 0.142, 536330.0*R, Te, X, dWdTe, dWdrhok);
  }
  
  if (REACTION[2]) {
    add_to_dW_fw_2r3p ( specO2, speceminus,   specO2plus, speceminus, speceminus,   3.869e-8 * calA, 0.1166, 567360.0*R, Te, X, dWdTe, dWdrhok);
  }
  
  if (REACTION[3]){
    kf=2.0e-7 * pow (300.0/Te, 0.7);
    dkfdTe = -0.7 * 2.0e-7 * pow ( Te / 300.0, -1.7 ) / 300.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r2p ( specO2plus, speceminus,   specO, specO,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[4]){
    kf=2.8e-7 * pow (300.0/Te, 0.5);
    dkfdTe = -0.5 * 2.8e-7 * pow ( Te / 300.0, -1.5 ) / 300.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r2p ( specN2plus, speceminus,   specN, specN,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  tmp = 2.0e-7 * pow (300.0/T, 0.5);
  if (REACTION[5]){
    kf=tmp;
    dkfdTe = 0.0;
    dkfdT= -0.5 * 2.0e-7 * pow ( T / 300.0, -1.5 ) / 300.0;
    dkfdTv=0.0;
    add_to_dW_2r2p ( specO2minus, specN2plus,   specO2, specN2,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[6]){
    kf=tmp;
    dkfdTe = 0.0;
    dkfdT= -0.5 * 2.0e-7 * pow ( T / 300.0, -1.5 ) / 300.0;
    dkfdTv=0.0;
    add_to_dW_2r2p ( specO2plus, specO2minus,   specO2, specO2,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  tmp = 2.0e-25 * pow (300.0/T, 2.5);
  if (REACTION[7]){
    kf=tmp;
    dkfdTe = 0.0;
    dkfdT= -2.5 * kf / T;
    dkfdTv=0.0;
    for (k=0; specM1[k]!=specEND; k++){
      add_to_dW_3r3p ( specO2minus, specN2plus, specM1[k],   specO2, specN2, specM1[k], kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[8]){
    kf=tmp;
    dkfdTe = 0.0;
    dkfdT= -2.5 * kf / T;
    dkfdTv=0.0;
    for (k=0; specM1[k]!=specEND; k++){
      add_to_dW_3r3p ( specO2minus, specO2plus, specM1[k],   specO2, specO2, specM1[k], kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }
  
  if (REACTION[9]){
    kf=1.4e-29 * 300.0 / Te * exp ( -600.0 / T ) * exp ( 700.0 * ( Te - T ) / Te / T );
    dkfdTe = ( -1.0 / Te + 700.0 / sqr ( Te ) ) * kf;
    dkfdT= ( 600.0 / sqr ( T ) - 700.0 / Te / T - 700.0 * ( Te - T ) / Te / sqr ( T ) ) * kf;
    dkfdTv=0.0;
    add_to_dW_3r2p ( speceminus, specO2, specO2,   specO2minus, specO2,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (REACTION[10]){
    kf=1.07e-31 * sqr ( 300.0 / Te ) * exp ( -70.0 / T ) * exp ( 1500.0 * ( Te - T ) / Te / T );
    dkfdTe = ( -2.0 / Te + 1500.0 / sqr ( Te ) ) * kf;
    dkfdT= ( 70.0 / sqr ( T ) - 1500.0 / sqr ( T ) ) * kf;
    dkfdTv=0.0;
    add_to_dW_3r2p ( speceminus, specO2, specN2,   specO2minus, specN2,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (REACTION[11]){
    kf=8.6e-10 * exp ( -6030.0 / T ) * ( 1.0 - exp ( -1570.0 / T ) );
    dkfdTe = 0.0;
    dkfdT= kf * 6030.0 / sqr ( T ) - 8.6e-10 * exp ( -6030.0 / T ) * exp ( -1570.0 / T ) * 1570.0 / sqr ( T );
    dkfdTv=0.0;
    add_to_dW_2r3p ( specO2minus, specO2,   speceminus, specO2, specO2,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[12]){
    tmp = 2.0e11 * Qbeam / Nt;  /* Qbeam is in watts per m3, N in 1/cm3 */

    dWdrhok[specO2][specO2] += -tmp * _calM(specO2) / _calM(specO2);
    dWdrhok[specO2plus][specO2] += +tmp * _calM(specO2plus) / _calM(specO2);
    dWdrhok[speceminus][specO2] += +tmp * _calM(speceminus) / _calM(specO2);
    for ( s = 0; s < ns; s++ ) {
      dWdrhok[specO2][s] += +tmp * N[specO2] / Nt * _calM(specO2) / _calM(s);
      dWdrhok[specO2plus][s] += -tmp * N[specO2] / Nt * _calM(specO2plus) / _calM(s);
      dWdrhok[speceminus][s] += -tmp * N[specO2] / Nt * _calM(speceminus) / _calM(s);
    }
    dkfdQb = 2.0e11 / Nt * N[specO2] / calA * 1.0e6;
    dWdQbeam[specO2] += -dkfdQb * _calM(specO2);
    dWdQbeam[specO2plus] += +dkfdQb * _calM(specO2plus);
    dWdQbeam[speceminus] += +dkfdQb * _calM(speceminus);
  }
  
  if (REACTION[13]){
    tmp = 1.8e11 * Qbeam / Nt;  /* Qbeam is in watts per m3, N in 1/cm3 */

    dWdrhok[specN2][specN2] += -tmp * _calM(specN2) / _calM(specN2);
    dWdrhok[specN2plus][specN2] += +tmp * _calM(specN2plus) / _calM(specN2);
    dWdrhok[speceminus][specN2] += +tmp * _calM(speceminus) / _calM(specN2);
    for ( s = 0; s < ns; s++ ) {
      dWdrhok[specN2][s] += +tmp * N[specN2] / Nt * _calM(specN2) / _calM(s);
      dWdrhok[specN2plus][s] += -tmp * N[specN2] / Nt * _calM(specN2plus) / _calM(s);
      dWdrhok[speceminus][s] += -tmp * N[specN2] / Nt * _calM(speceminus) / _calM(s);
    }
    dkfdQb = 1.8e11 / Nt * N[specN2] / calA * 1.0e6;
    dWdQbeam[specN2] += -dkfdQb * _calM(specN2);
    dWdQbeam[specN2plus] += +dkfdQb * _calM(specN2plus);
    dWdQbeam[speceminus] += +dkfdQb * _calM(speceminus);
  }
  
  tmp = exp ( -59380.0 / T ) * ( 1.0 - exp ( -2240.0 / T ) );
  dtmpdT = 59380.0 / sqr ( T ) * tmp - exp ( -59380.0 / T ) * 2240.0 / sqr ( T ) * exp ( -2240.0 / T );
  if (REACTION[14]){
    kf=3.7e-8 * tmp;
    dkfdTe = 0.0;
    dkfdT= 3.7e-8 * dtmpdT;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specO2, specO2,   specO, specO, specO2,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[15]){
    kf=9.3e-9 * tmp;
    dkfdTe = 0.0;
    dkfdT= 9.3e-9 * dtmpdT;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specO2, specN2,   specO, specO, specN2,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[16]){
    kf=1.3e-7 * tmp;
    dkfdTe = 0.0;
    dkfdT= 1.3e-7 * dtmpdT;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specO2, specO,   specO, specO, specO,  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  tmp = exp ( -113200.0 / T ) * ( 1.0 - exp ( -3354.0 / T ) );
  dtmpdT = 113200.0 / sqr ( T ) * tmp - exp ( -113200.0 / T ) * 3354.0 / sqr ( T ) * exp ( -3354.0 / T );
  if (REACTION[17]){
    kf=5.0e-8 * tmp;
    dkfdTe = 0.0;
    dkfdT= 5.0e-8 * dtmpdT;
    dkfdTv=0.0;
    for (k=0; specM1[k]!=specEND; k++){
      add_to_dW_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[18]){
    kf=1.1e-7 * tmp;
    dkfdTe = 0.0;
    dkfdT= 1.1e-7 * dtmpdT;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specN2, specO,   specN, specN, specO, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  tmp = pow ( T, -0.63 );
  if (REACTION[19]){
    kf=2.45e-31 * tmp;
    dkfdTe = 0.0;
    dkfdT= 2.45e-31 * -0.63 * pow ( T, -1.63 );
    dkfdTv=0.0;
    add_to_dW_3r2p ( specO, specO, specO2,   specO2, specO2, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[20]){
    kf=2.76e-34 * exp ( 720.0 / T );
    dkfdTe = 0.0;
    dkfdT= 2.76e-34 * -720.0 / sqr ( T ) * exp ( 720.0 / T );
    dkfdTv=0.0;
    add_to_dW_3r2p ( specO, specO, specN2,   specO2, specN2, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[21]){
    kf=8.8e-31 * tmp;
    dkfdTe = 0.0;
    dkfdT= 8.8e-31 * -0.63 * pow ( T, -1.63 );
    dkfdTv=0.0;
    add_to_dW_3r2p ( specO, specO, specO,   specO2, specO, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  tmp = 8.27e-34 * exp ( 500.0 / T );
  if (REACTION[22]){
    kf=tmp;
    dkfdTe = 0.0;
    dkfdT= -8.27e-34 * 500.0 / sqr ( T ) * exp ( 500.0 / T );
    dkfdTv=0.0;
    for (k=0; specM2[k]!=specEND; k++){
      add_to_dW_3r2p ( specN, specN, specM2[k],   specN2, specM2[k], kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

}



void find_Qei_Rajendran2022(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  *Qei=0.0;  

  if (REACTION[1])
    add_to_Qei(specN2, 2.6291e-8 * pow( Te, 0.142 ) * exp ( -536330.0 / Te ), rhok, Qei);
  if (REACTION[2]) 
    add_to_Qei(specO2, 3.869e-8 * pow( Te, 0.1166) * exp ( -567360.0 / Te ), rhok, Qei);
}



void find_dQei_dx_Rajendran2022(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  ;

  if (REACTION[1])
    add_to_dQei(specN2, 2.6291e-8 * pow( Te, 0.142 ) * exp ( -536330.0 / Te ), 0.0,  rhok, dQeidrhok, dQeidTe);
  if (REACTION[2]) 
    add_to_dQei(specO2, 3.869e-8 * pow( Te, 0.1166) * exp ( -567360.0 / Te ), 0.0,  rhok, dQeidrhok, dQeidTe);
}
