// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2018,2021 Bernard Parent

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

#define TEMAX_TOWNSEND 60000.0

const static bool REACTION[27]=
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
   TRUE, /* reaction 23 */
   TRUE, /* reaction 24 */
   TRUE, /* reaction 25 */
   TRUE, /* reaction 26 */
  };


#define specEND -1

const static long specM1[]=
  {
   specN, specNO, specEND
  };

const static long specM2[]=
  {
   specO, specNO, specO2, specEND
  };

const static long specM3[]=
  {
   specO2, specN2, specEND
  };

const static long specM4[]=
  {
   specO, specN, specNO, specEND
  };



void find_W_DunnKang1973 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;
  spec_t X;
  double R;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
  }

  R=1.9872;

  if (REACTION[1]){
    for (k=0; specM1[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM1[k], specO, specO, specM1[k], 3.6E18, -1.0, 59500.0*R, T, X, W );
      add_to_W_fw_3r2p ( specO, specO, specM1[k], specO2, specM1[k], 3.0E15, -0.5, 0.0*R, T, X, W );
    }
  }

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specN2, specM2[k], specN, specN, specM2[k], 1.9E17, -0.5, 113000.0*R, T, X, W );
      add_to_W_fw_3r2p ( specN, specN, specM2[k], specN2, specM2[k],  1.1E16, -0.5, 0.0*R, T, X, W );
    }
  }

  if (REACTION[3]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM3[k], specN, specO, specM3[k], 3.9e20, -1.5, 75500.0*R, T, X, W );
      add_to_W_fw_3r2p ( specN, specO, specM3[k], specNO, specM3[k], 1.0E20, -1.5, 0.0*R, T, X, W );
    }
  }

  if (REACTION[4]){
    add_to_W_fw_2r2p ( specO, specNO, specN, specO2, 3.2E9, 1.0, 19700.0*R, T, X, W );
    add_to_W_fw_2r2p ( specN, specO2, specO, specNO, 1.3e10, 1.0, 3580.0*R, T, X, W );
  }

  if (REACTION[5]){
    add_to_W_fw_2r2p ( specO, specN2, specN, specNO, 7.0E13, 0.0, 38000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specN, specNO, specO, specN2, 1.56E13, 0.0, 0.0, T, X, W );
  }

  if (REACTION[6]){
    add_to_W_fw_2r3p ( specN, specN2, specN, specN, specN, 4.085E22, -1.5, 113000.0*R, T, X, W );
    add_to_W_fw_3r2p ( specN, specN, specN, specN, specN2, 2.27e21, -1.5, 0.0, T, X, W );
  }

  if (REACTION[7]){
    add_to_W_fw_2r2p ( specO, specN, specNOplus, speceminus, 1.4E6, 1.5, 31900.0*R, T, X, W );
    add_to_W_fw_2r2p ( specNOplus, speceminus, specO, specN, 6.7e21, -1.5, 0.0, Te, X, W );
  }

  if (REACTION[8]) {
    add_to_W_fw_2r3p ( specO, speceminus, specOplus, speceminus, speceminus, 3.6E31, -2.91, 158000.0*R, min(TEMAX_TOWNSEND,Te), X, W );
    add_to_W_fw_3r2p ( specOplus, speceminus, speceminus, specO, speceminus, 2.2E40, -4.5, 0.0, Te, X, W );
  }

  if (REACTION[9]){
    add_to_W_fw_2r3p ( specN, speceminus, specNplus, speceminus, speceminus, 1.1E32, -3.14, 169000.0*R, min(TEMAX_TOWNSEND,Te), X, W );
    add_to_W_fw_3r2p ( specNplus, speceminus, speceminus, specN, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, W );
  }

  if (REACTION[10]){
    add_to_W_fw_2r2p ( specO, specO, specO2plus, speceminus, 1.6E17, -0.98, 80800.0*R, T, X, W );
    add_to_W_fw_2r2p ( specO2plus, speceminus, specO, specO, 2.2E40, -4.50, 0.0, Te, X, W );
  }

  if (REACTION[11]){
    add_to_W_fw_2r2p ( specO, specO2plus, specO2, specOplus, 2.92E18, -1.11, 28000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specO2, specOplus, specO, specO2plus, 7.8e11, 0.5, 0.0, T, X, W );
  }

  if (REACTION[12]){
    add_to_W_fw_2r2p ( specN2, specNplus, specN, specN2plus, 2.02e11, 0.81, 13000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specN, specN2plus, specN2, specNplus, 7.8e11, 0.5, 0.0, T, X, W );
  }

  if (REACTION[13]){
    add_to_W_fw_2r2p ( specN, specN, specN2plus, speceminus, 1.4e13, 0.0, 67800.0*R, T, X, W );
    add_to_W_fw_2r2p ( specN2plus, speceminus, specN, specN, 1.5e22, -1.5, 0.0, Te, X, W );
  }

  if (REACTION[14]){
    add_to_W_fw_2r3p ( specO2, specN2, specNO, specNOplus, speceminus, 1.38E20, -1.84, 141000.0*R, T, X, W );
    add_to_W_fw_3r2p ( specNO, specNOplus, speceminus, specO2, specN2, 1.0E24, -2.5, 0.0, Te, X, W );
  }

  if (REACTION[15]){
    add_to_W_fw_2r3p ( specNO, specN2, specNOplus, speceminus, specN2, 2.2E15, -0.35, 108000.0*R, T, X, W );
    add_to_W_fw_3r2p ( specNOplus, speceminus, specN2, specNO, specN2, 2.2E26, -2.5, 0.0, Te, X, W );
  }

  if (REACTION[16]){
    add_to_W_fw_2r2p ( specO, specNOplus, specNO, specOplus, 3.63E15, -0.6, 50800.0*R, T, X, W );
    add_to_W_fw_2r2p ( specNO, specOplus, specO, specNOplus, 1.5E13, 0.0, 0.0, T, X, W );
  }

  if (REACTION[17]){
    add_to_W_fw_2r2p ( specN2, specOplus, specO, specN2plus, 3.4E19, -2.0, 23000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specO, specN2plus, specN2, specOplus, 2.48E19, -2.2, 0.0, T, X, W );
  }

  if (REACTION[18]){
    add_to_W_fw_2r2p ( specN, specNOplus, specNO, specNplus, 1.0E19, -0.93, 61000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specNO, specNplus, specN, specNOplus, 4.8e14, 0.0, 0.0, T, X, W );
  }

  if (REACTION[19]){
    add_to_W_fw_2r2p ( specO2, specNOplus, specNO, specO2plus, 1.8E15, 0.17, 33000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specNO, specO2plus, specO2, specNOplus, 1.8E13, 0.5, 0.0, T, X, W );
  }

  if (REACTION[20]){
    add_to_W_fw_2r2p ( specO, specNOplus, specO2, specNplus, 1.34E13, 0.31, 77270.0*R, T, X, W );
    add_to_W_fw_2r2p ( specO2, specNplus, specO, specNOplus, 1e14, 0.0, 0.0, T, X, W );
  }

  if (REACTION[21]){
    add_to_W_fw_2r3p ( specNO, specO2, specNOplus, speceminus, specO2, 8.8e15, -0.35, 108000.0*R, T, X, W );
    add_to_W_fw_3r2p ( specNOplus, speceminus, specO2, specNO, specO2, 8.8e26, -2.5, 0.0, Te, X, W );
  }

  if (REACTION[22]){
    add_to_W_fw_2r3p ( specO2, specO, specO, specO, specO, 9E19, -1.0, 59500.0*R, T, X, W );
    add_to_W_fw_3r2p ( specO, specO, specO, specO2, specO, 7.5E16, -0.5, 0.0, T, X, W );
  }

  if (REACTION[23]){
    add_to_W_fw_2r3p ( specO2, specO2, specO, specO, specO2, 3.24E19, -1.0, 59500.0*R, T, X, W );
    add_to_W_fw_3r2p ( specO, specO, specO2, specO2, specO2, 2.7e16, -0.5, 0.0, T, X, W );
  }

  if (REACTION[24]){
    add_to_W_fw_2r3p ( specO2, specN2, specO, specO, specN2, 7.2E18, -1.0, 59500.0*R, T, X, W );
    add_to_W_fw_3r2p ( specO, specO, specN2, specO2, specN2, 6.0e15, -0.5, 0.0, T, X, W );
  }

  if (REACTION[25]){
    add_to_W_fw_2r3p ( specN2, specN2, specN, specN, specN2, 4.7E17, -0.5, 113000.0*R, T, X, W );
    add_to_W_fw_3r2p ( specN, specN, specN2, specN2, specN2, 2.72e16, -0.5, 0.0, T, X, W );
  }

  if (REACTION[26]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM4[k], specN, specO, specM4[k], 7.8E20, -1.5, 75500.0*R, T, X, W );
      add_to_W_fw_3r2p ( specN, specO, specM4[k], specNO, specM4[k], 2.0E20, -1.5, 0.0, T, X, W );
    }
  }

}


void find_dW_dx_DunnKang1973 ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */
  spec_t X;
  double R;
  
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }

  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }

  R=1.9872;


  if (REACTION[1]){
    for (k=0; specM1[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM1[k], specO, specO, specM1[k], 3.6E18, -1.0, 59500.0*R, T, X, dWdT, dWdrhok );
      add_to_dW_fw_3r2p ( specO, specO, specM1[k], specO2, specM1[k], 3.0E15, -0.5, 0.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specN2, specM2[k], specN, specN, specM2[k], 1.9E17, -0.5, 113000.0*R, T, X, dWdT, dWdrhok );
      add_to_dW_fw_3r2p ( specN, specN, specM2[k], specN2, specM2[k],  1.1E16, -0.5, 0.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[3]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM3[k], specN, specO, specM3[k], 3.9e20, -1.5, 75500.0*R, T, X, dWdT, dWdrhok );
      add_to_dW_fw_3r2p ( specN, specO, specM3[k], specNO, specM3[k], 1.0E20, -1.5, 0.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[4]){
    add_to_dW_fw_2r2p ( specO, specNO, specN, specO2, 3.2E9, 1.0, 19700.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specN, specO2, specO, specNO, 1.3e10, 1.0, 3580.0*R, T, X, dWdT, dWdrhok );
  }

  if (REACTION[5]){
    add_to_dW_fw_2r2p ( specO, specN2, specN, specNO, 7.0E13, 0.0, 38000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specN, specNO, specO, specN2, 1.56E13, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[6]){
    add_to_dW_fw_2r3p ( specN, specN2, specN, specN, specN, 4.085E22, -1.5, 113000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specN, specN, specN, specN, specN2, 2.27e21, -1.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[7]){
    add_to_dW_fw_2r2p ( specO, specN, specNOplus, speceminus, 1.4E6, 1.5, 31900.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specNOplus, speceminus, specO, specN, 6.7e21, -1.5, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[8]) {
    add_to_dW_fw_2r3p ( specO, speceminus, specOplus, speceminus, speceminus, 3.6E31, -2.91, 158000.0*R, min(TEMAX_TOWNSEND,Te), X, dWdTe, dWdrhok );
    add_to_dW_fw_3r2p ( specOplus, speceminus, speceminus, specO, speceminus, 2.2E40, -4.5, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[9]){
    add_to_dW_fw_2r3p ( specN, speceminus, specNplus, speceminus, speceminus, 1.1E32, -3.14, 169000.0*R, min(TEMAX_TOWNSEND,Te), X, dWdTe, dWdrhok );
    add_to_dW_fw_3r2p ( specNplus, speceminus, speceminus, specN, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[10]){
    add_to_dW_fw_2r2p ( specO, specO, specO2plus, speceminus, 1.6E17, -0.98, 80800.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specO2plus, speceminus, specO, specO, 2.2E40, -4.50, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[11]){
    add_to_dW_fw_2r2p ( specO, specO2plus, specO2, specOplus, 2.92E18, -1.11, 28000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specO2, specOplus, specO, specO2plus, 7.8e11, 0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[12]){
    add_to_dW_fw_2r2p ( specN2, specNplus, specN, specN2plus, 2.02e11, 0.81, 13000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specN, specN2plus, specN2, specNplus, 7.8e11, 0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[13]){
    add_to_dW_fw_2r2p ( specN, specN, specN2plus, speceminus, 1.4e13, 0.0, 67800.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specN2plus, speceminus, specN, specN, 1.5e22, -1.5, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[14]){
    add_to_dW_fw_2r3p ( specO2, specN2, specNO, specNOplus, speceminus, 1.38E20, -1.84, 141000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specNO, specNOplus, speceminus, specO2, specN2, 1.0E24, -2.5, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[15]){
    add_to_dW_fw_2r3p ( specNO, specN2, specNOplus, speceminus, specN2, 2.2E15, -0.35, 108000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specNOplus, speceminus, specN2, specNO, specN2, 2.2E26, -2.5, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[16]){
    add_to_dW_fw_2r2p ( specO, specNOplus, specNO, specOplus, 3.63E15, -0.6, 50800.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specNO, specOplus, specO, specNOplus, 1.5E13, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[17]){
    add_to_dW_fw_2r2p ( specN2, specOplus, specO, specN2plus, 3.4E19, -2.0, 23000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specO, specN2plus, specN2, specOplus, 2.48E19, -2.2, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[18]){
    add_to_dW_fw_2r2p ( specN, specNOplus, specNO, specNplus, 1.0E19, -0.93, 61000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specNO, specNplus, specN, specNOplus, 4.8e14, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[19]){
    add_to_dW_fw_2r2p ( specO2, specNOplus, specNO, specO2plus, 1.8E15, 0.17, 33000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specNO, specO2plus, specO2, specNOplus, 1.8E13, 0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[20]){
    add_to_dW_fw_2r2p ( specO, specNOplus, specO2, specNplus, 1.34E13, 0.31, 77270.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_2r2p ( specO2, specNplus, specO, specNOplus, 1e14, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[21]){
    add_to_dW_fw_2r3p ( specNO, specO2, specNOplus, speceminus, specO2, 8.8e15, -0.35, 108000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specNOplus, speceminus, specO2, specNO, specO2, 8.8e26, -2.5, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[22]){
    add_to_dW_fw_2r3p ( specO2, specO, specO, specO, specO, 9E19, -1.0, 59500.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specO, specO, specO, specO2, specO, 7.5E16, -0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[23]){
    add_to_dW_fw_2r3p ( specO2, specO2, specO, specO, specO2, 3.24E19, -1.0, 59500.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specO, specO, specO2, specO2, specO2, 2.7e16, -0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[24]){
    add_to_dW_fw_2r3p ( specO2, specN2, specO, specO, specN2, 7.2E18, -1.0, 59500.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specO, specO, specN2, specO2, specN2, 6.0e15, -0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[25]){
    add_to_dW_fw_2r3p ( specN2, specN2, specN, specN, specN2, 4.7E17, -0.5, 113000.0*R, T, X, dWdT, dWdrhok );
    add_to_dW_fw_3r2p ( specN, specN, specN2, specN2, specN2, 2.72e16, -0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[26]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM4[k], specN, specO, specM4[k], 7.8E20, -1.5, 75500.0*R, T, X, dWdT, dWdrhok );
      add_to_dW_fw_3r2p ( specN, specO, specM4[k], specNO, specM4[k], 2.0E20, -1.5, 0.0, T, X, dWdT, dWdrhok );
    }
  }



}


void find_Qei_DunnKang1973(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){

    Te=min(TEMAX_TOWNSEND,Te);
    if (REACTION[8]) 
      add_to_Qei(specO, 3.6e31/sqr(calA)*pow(Te,-2.91)*exp(-158000.0/Te), rhok, Qei);
    if (REACTION[9]) 
      add_to_Qei(specN, 1.1e32/sqr(calA)*pow(Te,-3.14)*exp(-169000.0/Te), rhok, Qei);
 
}



void find_dQei_dx_DunnKang1973(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
    Te=min(TEMAX_TOWNSEND,Te);
    if (REACTION[8]) 
      add_to_dQei(specO, 3.6e31/sqr(calA)*pow(Te,-2.91)*exp(-158000.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[9]) 
      add_to_dQei(specN, 1.1e32/sqr(calA)*pow(Te,-3.14)*exp(-169000.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);

}
