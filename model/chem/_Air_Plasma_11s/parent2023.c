// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Prasanna Thoguluva Rajendran
Copyright 2022 Bernard Parent

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

#define EXPON_TTV 1.0
#define EXPON_TTE 0.0
#define EXPON_TVTE 0.0

/* set all reactions to true except for testing purposes */
const static bool REACTION[29]=
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
   TRUE, /* reaction 27 */
   TRUE, /* reaction 28 */
  };

#define specEND -1

const static long specM1[]=
  {
   specN, specO, specNplus, specOplus, specEND
  };

const static long specM2[]=
  {
   specN2, specO2, specNO, specN2plus, specO2plus, specNOplus, specEND
  };

const static long specM3[]=
  {
   specN, specO, specNO, specNplus, specOplus,  specEND
  };

const static long specM4[]=
  {
   specN2, specO2, specN2plus, specO2plus, specNOplus,  specEND
  };


void find_W_Parent2023 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                       double Qbeam, spec_t W ) {
  double R,TTv,TvTe,TTe;
  long k;
  spec_t X;

  /* find properties needed by add_to_W* functions */
  R=Rchem;
  TTv=pow(T,EXPON_TTV)*pow(Tv,1.0-EXPON_TTV);
  TTe=pow(T,EXPON_TTE)*pow(Te,1.0-EXPON_TTE);
  TvTe=pow(Tv,EXPON_TVTE)*pow(Te,1.0-EXPON_TVTE);

  
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
  }


  if (REACTION[1]) {
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, T, X, W );
    }
  }
      

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, T, X, W );
    }
  }

  if (REACTION[3]){
    add_to_W_fw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 3.0e24, -1.6, 113200.0*R, TvTe, X, W );
    add_to_W_bw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 3.0e24, -1.6, 113200.0*R, TTe, X, W );
  }

  if (REACTION[4]){
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, T, X, W );
    }
  }

  if (REACTION[5]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, T, X, W );
    }
  }

  if (REACTION[6]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, T, X, W );
    }
  }

  if (REACTION[7]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, T, X, W );
    }
  }

  if (REACTION[8]){
    add_to_W_fwbw_2r2p ( specNO, specO,   specN, specO2, 8.4e12, 0.0, 19400.0*R, T, X, W );
  }

  if (REACTION[9]){
    add_to_W_fwbw_2r2p ( specN2, specO,   specNO, specN, 5.7e12, 0.42, 42938.0*R, T, X, W );
  }

  if (REACTION[10]){
    add_to_W_fw_2r2p ( specN, specO,   specNOplus, speceminus, 5.3e12, 0.0, 32000.0*R, T, X, W );
    add_to_W_bw_2r2p ( specN, specO,   specNOplus, speceminus, 5.3e12, 0.0, 32000.0*R, Te, X, W );
  }

  if (REACTION[11]){
    add_to_W_fw_2r2p ( specO, specO,   specO2plus, speceminus, 1.1e13, 0.0, 81200.0*R, T, X, W );
    add_to_W_bw_2r2p ( specO, specO,   specO2plus, speceminus, 1.1e13, 0.0, 81200.0*R, Te, X, W );
  }

  if (REACTION[12]){
    add_to_W_fw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, T, X, W );
    add_to_W_bw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, Te, X, W );
  }

  if (REACTION[13]){
    add_to_W_fwbw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.0e12, 0.5, 77200.0*R, T, X, W );
  }

  if (REACTION[14]){
    add_to_W_fwbw_2r2p ( specNplus, specN2,   specN2plus, specN, 1.0e12, 0.5, 12200.0*R, T, X, W );
  }

  if (REACTION[15]){
    add_to_W_fwbw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.7e13, 0.14, 28600.0*R, T, X, W );
  }

  if (REACTION[16]){
    add_to_W_fwbw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.4e5, 1.90, 26600.0*R, T, X, W );
  }

  if (REACTION[17]){
    add_to_W_fwbw_2r2p ( specO2plus, specN2,   specN2plus, specO2, 9.9e12, 0.00, 40700.0*R, T, X, W );
  }

  if (REACTION[18]){
    add_to_W_fwbw_2r2p ( specO2plus, specO,   specOplus, specO2, 4.0e12, -0.09, 18000.0*R, T, X, W );
  }

  if (REACTION[19]){
    add_to_W_fwbw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.4e13, -1.08, 12800.0*R, T, X, W );
  }

  if (REACTION[20]){
    add_to_W_fwbw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.4e13, 0.41, 32600.0*R, T, X, W );
  }

  if (REACTION[21]){
    add_to_W_fwbw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.2e12, 0.29, 48600.0*R, T, X, W );
  }

  if (REACTION[22]){
    add_to_W_fwbw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.1e11, 0.36, 22800.0*R, T, X, W );
  }

  if (REACTION[23]){
    add_to_W_fwbw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.2e13, 0.00, 35500.0*R, T, X, W );
  }

  if (REACTION[24]) {
    add_to_W_fwbw_2r3p ( specO, speceminus, specOplus, speceminus, speceminus, 6.37e16, 0.0029, 477190.0*R, Te, X, W );
  }

  if (REACTION[25]){
    add_to_W_fwbw_2r3p ( specN, speceminus, specNplus, speceminus, speceminus, 1.06e18, -0.2072, 629700.0*R, Te, X, W );
  }

  if (REACTION[26]){
    add_to_W_fwbw_2r3p ( specO2, speceminus, specO2plus, speceminus, speceminus, 2.33e16, 0.1166, 567360.0*R, Te, X, W );
  }

  if (REACTION[27]){
    add_to_W_fwbw_2r3p ( specN2, speceminus, specN2plus, speceminus, speceminus, 1.58e16, 0.1420, 536330.0*R, Te, X, W );
  }

  if (REACTION[28]){
    add_to_W_fwbw_2r3p ( specNO, speceminus, specNOplus, speceminus, speceminus, 5.63e18, -0.2607, 686030.0*R, Te, X, W );
  }


}



/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx_Parent2023 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s, spec;                    
  double TTv,TTe,TvTe,R;
  spec_t dWdTTv,dWdTTe,dWdTvTe;
  spec_t X;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  
  R=Rchem;
  TTv=pow(T,EXPON_TTV)*pow(Tv,1.0-EXPON_TTV);
  TTe=pow(T,EXPON_TTE)*pow(Te,1.0-EXPON_TTE);
  TvTe=pow(Tv,EXPON_TVTE)*pow(Te,1.0-EXPON_TVTE);
  /* initialize all derivatives to zero */
  for ( s = 0; s < ns; s++ ) {
    dWdTTv[s] = 0.0;
    dWdTTe[s] = 0.0;
    dWdTvTe[s] = 0.0;
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }



  if (REACTION[1]) {
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, T, X, dWdT, dWdrhok );
    }
  }
      

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[3]){
    add_to_dW_fw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 3.0e24, -1.6, 113200.0*R, TvTe, X, dWdTvTe, dWdrhok );
    add_to_dW_bw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 3.0e24, -1.6, 113200.0*R, TTe, X, dWdTTe, dWdrhok );
  }

  if (REACTION[4]){
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[5]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[6]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[7]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, T, X, dWdT, dWdrhok );
    }
  }



  if (REACTION[8]){
    add_to_dW_fwbw_2r2p ( specNO, specO,   specN, specO2, 8.4e12, 0.0, 19400.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[9]){
    add_to_dW_fwbw_2r2p ( specN2, specO,   specNO, specN, 5.7e12, 0.42, 42938.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[10]){
    add_to_dW_fw_2r2p ( specN, specO,   specNOplus, speceminus, 5.3e12, 0.0, 32000.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_bw_2r2p ( specN, specO,   specNOplus, speceminus, 5.3e12, 0.0, 32000.0*R, Te, X, dWdTe, dWdrhok);
  }

  if (REACTION[11]){
    add_to_dW_fw_2r2p ( specO, specO,   specO2plus, speceminus, 1.1e13, 0.0, 81200.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_bw_2r2p ( specO, specO,   specO2plus, speceminus, 1.1e13, 0.0, 81200.0*R, Te, X, dWdTe, dWdrhok);
  }

  if (REACTION[12]){
    add_to_dW_fw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_bw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, Te, X, dWdTe, dWdrhok);
  }


  if (REACTION[13]){
    add_to_dW_fwbw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.0e12, 0.5, 77200.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[14]){
    add_to_dW_fwbw_2r2p ( specNplus, specN2,   specN2plus, specN, 1.0e12, 0.5, 12200.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[15]){
    add_to_dW_fwbw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.7e13, 0.14, 28600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[16]){
    add_to_dW_fwbw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.4e5, 1.90, 26600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[17]){
    add_to_dW_fwbw_2r2p ( specO2plus, specN2,   specN2plus, specO2, 9.9e12, 0.00, 40700.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[18]){
    add_to_dW_fwbw_2r2p ( specO2plus, specO,   specOplus, specO2, 4.0e12, -0.09, 18000.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[19]){
    add_to_dW_fwbw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.4e13, -1.08, 12800.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[20]){
    add_to_dW_fwbw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.4e13, 0.41, 32600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[21]){
    add_to_dW_fwbw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.2e12, 0.29, 48600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[22]){
    add_to_dW_fwbw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.1e11, 0.36, 22800.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[23]){
    add_to_dW_fwbw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.2e13, 0.00, 35500.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[24]) {
    add_to_dW_fwbw_2r3p ( specO, speceminus, specOplus, speceminus, speceminus, 6.37e16, 0.0029, 477190.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[25]){
    add_to_dW_fwbw_2r3p ( specN, speceminus, specNplus, speceminus, speceminus, 1.06e18, -0.2072, 629700.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[26]){
    add_to_dW_fwbw_2r3p ( specO2, speceminus, specO2plus, speceminus, speceminus, 2.33e16, 0.1166, 567360.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[27]){
    add_to_dW_fwbw_2r3p ( specN2, speceminus, specN2plus, speceminus, speceminus, 1.58e16, 0.1420, 536330.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[28]){
    add_to_dW_fwbw_2r3p ( specNO, speceminus, specNOplus, speceminus, speceminus, 5.63e18, -0.2607, 686030.0*R, Te, X, dWdTe, dWdrhok );
  }



  for (spec=0; spec<ns; spec++){
    dWdT[spec]+=dWdTTv[spec]*EXPON_TTV*TTv/T;
    dWdTv[spec]+=dWdTTv[spec]*(1.0-EXPON_TTV)*TTv/Tv;

    dWdT[spec]+=dWdTTe[spec]*EXPON_TTE*TTe/T;
    dWdTe[spec]+=dWdTTe[spec]*(1.0-EXPON_TTE)*TTe/Te;

    dWdTv[spec]+=dWdTvTe[spec]*EXPON_TVTE*TvTe/Tv;
    dWdTe[spec]+=dWdTvTe[spec]*(1.0-EXPON_TVTE)*TvTe/Te;

  }
  
}



void find_Qei_Parent2023(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){

    if (REACTION[24]) 
      add_to_Qei(specO,_ionizationpot(specO), 6.37e16/calA*pow(Te,0.0029)*exp(-477190.0/Te), rhok, Qei);
    if (REACTION[25]) 
      add_to_Qei(specN,_ionizationpot(specN), 1.06e18/calA*pow(Te,-0.2072)*exp(-629700.0/Te), rhok, Qei);
    if (REACTION[26]) 
      add_to_Qei(specO2,_ionizationpot(specO2), 2.33e16/calA*pow(Te,0.1166)*exp(-567360.0/Te), rhok, Qei);
    if (REACTION[27]) 
      add_to_Qei(specN2,_ionizationpot(specN2), 1.58e16/calA*pow(Te,0.1420)*exp(-536330.0/Te), rhok, Qei);
    if (REACTION[28]) 
      add_to_Qei(specNO,_ionizationpot(specNO), 5.63e18/calA*pow(Te,-0.2607)*exp(-686030.0/Te), rhok, Qei);
 
}



void find_dQei_dx_Parent2023(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){

    if (REACTION[24]) 
      add_to_dQei(specO,_ionizationpot(specO), 6.37e16/calA*pow(Te,0.0029)*exp(-477190.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[25]) 
      add_to_dQei(specN,_ionizationpot(specN), 1.06e18/calA*pow(Te,-0.2072)*exp(-629700.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[26]) 
      add_to_dQei(specO2,_ionizationpot(specO2), 2.33e16/calA*pow(Te,0.1166)*exp(-567360.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[27]) 
      add_to_dQei(specN2,_ionizationpot(specN2), 1.58e16/calA*pow(Te,0.1420)*exp(-536330.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[28]) 
      add_to_dQei(specNO,_ionizationpot(specNO), 5.63e18/calA*pow(Te,-0.2607)*exp(-686030.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);

}
