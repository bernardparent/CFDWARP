// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2023 Ajjay Omprakas

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
#include <model/share/model_share.h>


const static bool REACTION[20]=
  {
   FALSE, /* reaction 0 */
   TRUE, /* reaction 1 */
   TRUE, /* reaction 2 */
   TRUE, /* reaction 3 */
   TRUE, /* reaction 4 */
   TRUE, /* reaction 5 */
   TRUE, /* reaction 6 */
   TRUE, /* reaction 7 */
   TRUE, /* reaction 8 */
  };



void find_W_bavafa2008 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;
  spec_t X, N;
  double R;
  double theta;

  //Estar=3.6e-21;
  theta = log(Estar);

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
    W[k] = 0.0;
  }

  R=Rchem;
  
  
  if (REACTION[1]){
    add_to_W_2r3p ( specNH3, speceminus, specNH3plus, speceminus, speceminus, _averaged_rate(np,gl,1,min( exp ( 1.2346e-14 * pow ( theta, 9.0 ) + 0.3115 * theta ), 1.5640e-07 )), N, W );
  }

  if (REACTION[2]){
    add_to_W_2r4p ( specNH3, speceminus, specNH2plus, specH, speceminus, speceminus, _averaged_rate(np,gl,2,min( exp ( 0.4853 * theta + 5.5515e-131 * pow( theta, 79.0 ) ), 2.4490e-10 )), N, W );
  }

  if (REACTION[3]){
    add_to_W_2r3p ( specNH3, speceminus, specNH2, specH, speceminus, _averaged_rate(np,gl,3,min( exp ( 2.0636e-07 * pow( theta, 5.0 ) - 4.3548e-11 * pow( theta, 7.0 ) ), 1.0400e-08 * exp( -1.1290e+19 * Estar ) + 2.0460e-09 * exp( -3.6270e05 * Estar ) )), N, W );
  }

  if (REACTION[4]){
    add_to_W_2r4p ( specNH3, speceminus, specNH, specH, specH, speceminus, _averaged_rate(np,gl,4,min( exp ( 1.4602e-07 * pow ( theta, 5.0 ) - 1.3747e-21 * pow( theta, 13.0 ) ), 1.200e-08 )), N, W );
  }

  if (REACTION[5]){
    add_to_W_fw_2r1p ( specNH2, specH, specNH3, 1.0000E-13*calA, 0.0, 0.0*R, T, X, W );
  }

  if (REACTION[6]){
    add_to_W_fw_2r2p ( specNH2, specNH2, specNH, specNH3, 1.4000E-12*calA, 0.0, 0.0*R, T, X, W );
  }


  if (REACTION[7]){
    add_to_W_fw_2r2p ( specNH3, specNH3plus,   specNH2, specNH4plus, 2.2000E-09*calA, 0.0, 0.0*R, T, X, W );
  }


  if (REACTION[8]) {
    add_to_W_fw_2r2p ( specNH3, specNH2plus, specNH2, specNH3plus, 1.0000E-09*calA, 0.0, 0.0*R, T, X, W );
  }
  
}


void find_dW_dx_bavafa2008 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */
  spec_t X, N;
  double R;
  double dkfdTe,dkfdT,dkfdTv,kf;
  double theta;
  
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA;
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

  //Estar=3.6e-21;
  /* find properties needed by add_to_dW* functions in proper units */
  R=Rchem;
  theta = log(Estar);
  
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA; /* particules/cm^3 */
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }


  if (REACTION[1]){
    kf = min( exp ( 1.2346e-14 * pow ( theta, 9.0 ) + 0.3115 * theta ), 1.5640e-07 );
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specNH3, speceminus, specNH3plus, speceminus, speceminus, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }

  if (REACTION[2]){
    kf = min( exp ( 0.4853 * theta + 5.5515e-131 * pow( theta, 79.0 ) ), 2.4490e-10 );
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r4p ( specNH3, speceminus, specNH2plus, specH, speceminus, speceminus, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }

  if (REACTION[3]){
    kf = min( exp ( 2.0636e-07 * pow( theta, 5.0 ) - 4.3548e-11 * pow( theta, 7.0 ) ), 1.0400e-08 * exp( -1.1290e+19 * Estar ) + 2.0460e-09 * exp( -3.6270e05 * Estar ) );
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specNH3, speceminus, specNH2, specH, speceminus, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }

  if (REACTION[4]){
    kf = min( exp ( 1.4602e-07 * pow ( theta, 5.0 ) - 1.3747e-21 * pow( theta, 13.0 ) ), 1.200e-08 );
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r4p ( specNH3, speceminus, specNH, specH, specH, speceminus, kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }

  if (REACTION[5]){
    add_to_dW_fw_2r1p ( specNH2, specH, specNH3, 1.0000E-13*calA, 0.0, 0.0*R, T, X, dWdT, dWdrhok );
  }

  if (REACTION[6]){
    add_to_dW_fw_2r2p ( specNH2, specNH2, specNH, specNH3, 1.4000E-12*calA, 0.0, 0.0*R, T, X, dWdT, dWdrhok );
  }

  if (REACTION[7]){
    add_to_dW_fw_2r2p ( specNH3, specNH3plus,   specNH2, specNH4plus, 2.2000E-09*calA, 0.0, 0.0*R, T, X, dWdT, dWdrhok );
  }


  if (REACTION[8]) {
    add_to_dW_fw_2r2p ( specNH3, specNH2plus, specNH2, specNH3plus, 1.0000E-09*calA, 0.0, 0.0*R, T, X, dWdT, dWdrhok );
  }
  
}


void find_Qei_bavafa2008(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
 
  double theta;
  *Qei=0.0;  
  theta=log(Estar);

  if (REACTION[1]) {
    add_to_Qei(specNH3, _ionizationpot(specNH3), min( exp ( 1.2346e-14 * pow ( theta, 9.0 ) + 0.3115 * theta ), 1.5640e-07 ), rhok, Qei);
  }
  if (REACTION[2]) {
    add_to_Qei(specNH3, _ionizationpot(specNH3), min( exp ( 0.4853 * theta + 5.5515e-131 * pow( theta, 79.0 ) ), 2.4490e-10 ), rhok, Qei);
  }   

}



void find_dQei_dx_bavafa2008(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  
  double theta;
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
  theta=log(Estar);

  if (REACTION[1]) {
    add_to_dQei(specNH3, _ionizationpot(specNH3), min( exp ( 1.2346e-14 * pow ( theta, 9.0 ) + 0.3115 * theta ), 1.5640e-07 ), 0.0,  rhok, dQeidrhok, dQeidTe);
  }
  if (REACTION[2]) {
    add_to_dQei(specNH3, _ionizationpot(specNH3), min( exp ( 0.4853 * theta + 5.5515e-131 * pow( theta, 79.0 ) ), 2.4490e-10 ), 0.0,  rhok, dQeidrhok, dQeidTe);
  }

}


