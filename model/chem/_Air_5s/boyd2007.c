// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2020 Bernard Parent

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
const static bool REACTION[9]=
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
  };

#define specEND -1

const static long specM1[]=
  {
   specN, specO, specEND
  };

const static long specM2[]=
  {
   specN2, specO2, specNO, specEND
  };

const static long specM3[]=
  {
   specN, specO, specNO, specEND
  };

const static long specM4[]=
  {
   specN2, specO2, specEND
  };


void find_W_Boyd2007 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                       double Qbeam, spec_t W ) {
  double R,TTv;
  long k;
  spec_t X;

  /* find properties needed by add_to_W* functions */
  R=Rchem;
  TTv=sqrt(Tv*T);
  
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
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, T, X, W );
    }
  }

  if (REACTION[4]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, T, X, W );
    }
  }

  if (REACTION[5]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, T, X, W );
    }
  }

  if (REACTION[6]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, T, X, W );
    }
  }

  if (REACTION[7]){
    add_to_W_fwbw_2r2p ( specNO, specO,   specN, specO2, 8.4e12, 0.0, 19400.0*R, T, X, W );
  }

  if (REACTION[8]){
    add_to_W_fwbw_2r2p ( specN2, specO,   specNO, specN, 5.7e12, 0.42, 42938.0*R, T, X, W );
  }



}



/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx_Boyd2007 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
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
  TTv=sqrt(T*Tv);
  TTe=sqrt(T*Te);
  TvTe=sqrt(Tv*Te);
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
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[4]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[5]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[6]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, T, X, dWdT, dWdrhok );
    }
  }



  if (REACTION[7]){
    add_to_dW_fwbw_2r2p ( specNO, specO,   specN, specO2, 8.4e12, 0.0, 19400.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[8]){
    add_to_dW_fwbw_2r2p ( specN2, specO,   specNO, specN, 5.7e12, 0.42, 42938.0*R, T, X, dWdT, dWdrhok);
  }



  for (spec=0; spec<ns; spec++){
    dWdT[spec]+=dWdTTv[spec]*0.5/TTv*Tv;
    dWdTv[spec]+=dWdTTv[spec]*0.5/TTv*T;

    dWdT[spec]+=dWdTTe[spec]*0.5/TTe*Te;
    dWdTe[spec]+=dWdTTe[spec]*0.5/TTe*T;

    dWdTv[spec]+=dWdTvTe[spec]*0.5/TvTe*Te;
    dWdTe[spec]+=dWdTvTe[spec]*0.5/TvTe*Tv;
  }
  
}
