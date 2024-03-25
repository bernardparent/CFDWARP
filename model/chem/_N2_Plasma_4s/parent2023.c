// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Prasanna Thoguluva Rajendran
Copyright 2022,2023 Bernard Parent

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
   FALSE, /* reaction 0 */
   TRUE, /* reaction 1 */
   TRUE, /* reaction 2 */
   TRUE, /* reaction 3 */
   TRUE, /* reaction 4 */
   TRUE, /* reaction 5 */
   TRUE, /* reaction 6 */
  };

#define specEND -1

const static long specM1[]=
  {
   specN,  specEND
  };

const static long specM2[]=
  {
   specN2,  specN2plus, specEND
  };



void find_W_Parent2023 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                       double Qbeam, spec_t W ) {
  double tmp,Nt,R,TTv,TvTe,TTe;
  long k;
  spec_t X,N;

  /* find properties needed by add_to_W* functions */
  R=Rchem;
  TTv=pow(T,EXPON_TTV)*pow(Tv,1.0-EXPON_TTV);
  TTe=pow(T,EXPON_TTE)*pow(Te,1.0-EXPON_TTE);
  TvTe=pow(Tv,EXPON_TVTE)*pow(Te,1.0-EXPON_TVTE);

  
  Nt = 0.0;
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
    Nt += N[k];
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
    add_to_W_fw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, T, X, W );
    add_to_W_bw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, Te, X, W );
  }

  if (REACTION[5]){
    add_to_W_fwbw_2r3p ( specN2, speceminus, specN2plus, speceminus, speceminus, 1.58e16, 0.1420, 536330.0*R, Te, X, W );
  }

  if (REACTION[6]){
    tmp = 1.8e11 * Qbeam / Nt * N[specN2];
    W[specN2] += -tmp * _calM(specN2) / calA * 1.0e6;
    W[specN2plus] += +tmp * _calM(specN2plus) / calA * 1.0e6;
    W[speceminus] += +tmp * _calM(speceminus) / calA * 1.0e6;
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
  double dkfdQb,tmp,Nt,TTv,TTe,TvTe,R;
  spec_t dWdTTv,dWdTTe,dWdTvTe;
  spec_t X,N;

  Nt = 0.0;
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
    Nt += N[k];
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
    add_to_dW_fw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_bw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, Te, X, dWdTe, dWdrhok);
  }


  if (REACTION[5]){
    add_to_dW_fwbw_2r3p ( specN2, speceminus, specN2plus, speceminus, speceminus, 1.58e16, 0.1420, 536330.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[6]){
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

    if (REACTION[5]) 
      add_to_Qei(gl,Te,specN2,_ionizationpot(specN2), 1.58e16/calA*pow(Te,0.1420)*exp(-536330.0/Te), rhok, Qei);
 
}



void find_dQei_dx_Parent2023(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){

    if (REACTION[5]) 
      add_to_dQei(gl,Te,specN2,_ionizationpot(specN2), 1.58e16/calA*pow(Te,0.1420)*exp(-536330.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);

}
