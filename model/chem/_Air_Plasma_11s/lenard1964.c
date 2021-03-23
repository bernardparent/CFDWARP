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
const static bool REACTION[21]=
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
  };

#define specEND -1

const static long specM1[]=
  {
   specN, specN2, specNO, specNOplus, specEND
  };

const static long specM2[]=
  {
   specO, specO2, specNO, specNOplus, specEND
  };

const static long specM3[]=
  {
   specO, specO2, specN, specN2, specNOplus, specEND
  };



void find_W_Lenard1964 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  double N[ns];
  double R,kf;
  long k;
  
      
  /* find properties needed by add_to_W* functions */
  R=1.9872;
  for ( k = 0; k < ns; k++ ) {
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
  }


  if (REACTION[1])
    add_to_W_2r3p ( specO2, specO2,   specO, specO, specO2,  _kf_Arrhenius(2, 2.3e19, -1.0, 59400.0*R, T) , N, W);
    
  if (REACTION[2])
    add_to_W_2r3p ( specO2, specO,   specO, specO, specO,    _kf_Arrhenius(2, 8.5e19, -1.0, 59400.0*R, T) , N, W);
  
  if (REACTION[3]){
    kf=_kf_Arrhenius(2, 3.0e18, -1.0, 59400.0*R, T);
    for (k=0; specM1[k]!=specEND; k++)
      add_to_W_2r3p ( specO2, specM1[k],   specO, specO, specM1[k],  kf   , N, W); 
  }
  
  if (REACTION[4])
    add_to_W_2r3p ( specN2, specN2,   specN, specN, specN2,    _kf_Arrhenius(2, 3.8e19, -1.0, 113200.0*R, T) , N, W);

  if (REACTION[5])
    add_to_W_2r3p ( specN2, specN,   specN, specN, specN,    _kf_Arrhenius(2, 1.3e20, -1.0, 113200.0*R, T) , N, W);

  if (REACTION[6]){
    kf=_kf_Arrhenius(2, 1.9e19, -1.0, 113200.0*R, T);
    for (k=0; specM2[k]!=specEND; k++)
      add_to_W_2r3p ( specN2, specM2[k],   specN, specN, specM2[k],   kf , N, W);
  }
  
  if (REACTION[7]){
    kf=_kf_Arrhenius(2, 2.4e17, -0.5, 75500.0*R, T);
    for (k=0; specM3[k]!=specEND; k++)
      add_to_W_2r3p ( specNO, specM3[k],     specN, specO, specM3[k],     kf , N, W);
  }

  if (REACTION[8])
    add_to_W_2r2p ( specO, specN2,   specN, specNO,    _kf_Arrhenius(2, 6.8e13, 0.0, 37750.0*R, T) , N, W);

  if (REACTION[9])
    add_to_W_2r2p ( specO, specNO,   specN, specO2,    _kf_Arrhenius(2, 4.3e7, 1.5, 19100.0*R, T) , N, W);

  if (REACTION[10])
    add_to_W_2r2p ( specN, specO,    specNOplus, speceminus,    _kf_Arrhenius(2, 1.3e8, 1.0, 31900.0*R, T) , N, W);

  if (REACTION[11])
    add_to_W_3r2p ( specO, specO, specO2,   specO2, specO2,    _kf_Arrhenius(3, 1.9e16, -0.5, 0.0, T) , N, W);

  if (REACTION[12])
    add_to_W_3r2p ( specO, specO, specO,   specO2, specO,    _kf_Arrhenius(3, 7.1e16, -0.5, 0.0, T) , N, W);

  if (REACTION[13]){
    kf=_kf_Arrhenius(3, 2.5e15, -0.5, 0.0, T);
    for (k=0; specM1[k]!=specEND; k++)
      add_to_W_3r2p ( specO, specO, specM1[k],   specO2, specM1[k],    kf , N, W);
  }

  if (REACTION[14])
    add_to_W_3r2p ( specN, specN, specN2,   specN2, specN2,    _kf_Arrhenius(3, 2.0e18, -1.0, 0.0, T) , N, W);

  if (REACTION[15])
    add_to_W_3r2p ( specN, specN, specN,   specN2, specN,    _kf_Arrhenius(3, 7.0e18, -1.0, 0.0, T) , N, W);

  if (REACTION[16]){
    kf=_kf_Arrhenius(3, 1e18, -1.0, 0.0, T);
    for (k=0; specM2[k]!=specEND; k++)
      add_to_W_3r2p ( specN, specN, specM2[k],   specN2, specM2[k],    kf , N, W);
  }     

  if (REACTION[17]){
    kf=_kf_Arrhenius(3, 6e16, -0.5, 0.0, T);
    for (k=0; specM3[k]!=specEND; k++)
      add_to_W_3r2p ( specN, specO, specM3[k],   specNO, specM3[k],    kf , N, W);
  }     

  if (REACTION[18])
    add_to_W_2r2p ( specN, specNO,   specO, specN2,    _kf_Arrhenius(2, 1.5e13, 0.0, 0.0, T) , N, W);

  if (REACTION[19])
    add_to_W_2r2p ( specN, specO2,   specO, specNO,    _kf_Arrhenius(2, 1.8e8, 1.5, 3020.0*R, T) , N, W);

  if (REACTION[20])
    add_to_W_2r2p ( speceminus, specNOplus,   specN, specO,    _kf_Arrhenius(2, 2e19, -1.0, 0.0, Te) , N, W);


}






/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx_Lenard1964 ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, 
                  double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, 
                  spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    
  spec_t N;
  double R,kf,dkfdTe,dkfdT,dkfdTv;
  
  
  R=1.9872;
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
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA;
  }


  if (REACTION[1]) {
    kf=_kf_Arrhenius(2, 2.3e19, -1.0, 59400.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 2.3e19, -1.0, 59400.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_2r3p ( specO2, specO2,   specO, specO, specO2,  kf , N,  dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
    
  if (REACTION[2]) {
    kf=_kf_Arrhenius(2, 8.5e19, -1.0, 59400.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 8.5e19, -1.0, 59400.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_2r3p ( specO2, specO,   specO, specO, specO,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (REACTION[3]){
    kf=_kf_Arrhenius(2, 3.0e18, -1.0, 59400.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 3.0e18, -1.0, 59400.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    for (k=0; specM1[k]!=specEND; k++)
      add_to_dW_2r3p ( specO2, specM1[k],   specO, specO, specM1[k],  kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe); 
  }
  
  if (REACTION[4]) {
    kf=_kf_Arrhenius(2, 3.8e19, -1.0, 113200.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 3.8e19, -1.0, 113200.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_2r3p ( specN2, specN2,   specN, specN, specN2,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[5]) {
    kf=_kf_Arrhenius(2, 1.3e20, -1.0, 113200.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 1.3e20, -1.0, 113200.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_2r3p ( specN2, specN,   specN, specN, specN,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[6]){
    kf=_kf_Arrhenius(2, 1.9e19, -1.0, 113200.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 1.9e19, -1.0, 113200.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    for (k=0; specM2[k]!=specEND; k++)
      add_to_dW_2r3p ( specN2, specM2[k],   specN, specN, specM2[k],   kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (REACTION[7]){
    kf=_kf_Arrhenius(2, 2.4e17, -0.5, 75500.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 2.4e17, -0.5, 75500.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    for (k=0; specM3[k]!=specEND; k++)
      add_to_dW_2r3p ( specNO, specM3[k],     specN, specO, specM3[k],     kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[8]){
    kf=_kf_Arrhenius(2, 6.8e13, 0.0, 37750.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 6.8e13, 0.0, 37750.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_2r2p ( specO, specN2,   specN, specNO,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[9]){
    kf=_kf_Arrhenius(2, 4.3e7, 1.5, 19100.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 4.3e7, 1.5, 19100.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_2r2p ( specO, specNO,   specN, specO2,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[10]){
    kf=_kf_Arrhenius(2, 1.3e8, 1.0, 31900.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 1.3e8, 1.0, 31900.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_2r2p ( specN, specO,    specNOplus, speceminus,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[11]){
    kf=_kf_Arrhenius(3, 1.9e16, -0.5, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(3, 1.9e16, -0.5, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_3r2p ( specO, specO, specO2,   specO2, specO2,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[12]){
    kf=_kf_Arrhenius(3, 7.1e16, -0.5, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(3, 7.1e16, -0.5, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;
    add_to_dW_3r2p ( specO, specO, specO,   specO2, specO,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[13]){
    kf=_kf_Arrhenius(3, 2.5e15, -0.5, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(3, 2.5e15, -0.5, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;    
    for (k=0; specM1[k]!=specEND; k++)
      add_to_dW_3r2p ( specO, specO, specM1[k],   specO2, specM1[k],    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[14]){
    kf=_kf_Arrhenius(3, 2.0e18, -1.0, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(3, 2.0e18, -1.0, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;        
    add_to_dW_3r2p ( specN, specN, specN2,   specN2, specN2,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[15]){
    kf=_kf_Arrhenius(3, 7.0e18, -1.0, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(3, 7.0e18, -1.0, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;            
    add_to_dW_3r2p ( specN, specN, specN,   specN2, specN,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[16]){
    kf=_kf_Arrhenius(3, 1e18, -1.0, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(3, 1e18, -1.0, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;            
    for (k=0; specM2[k]!=specEND; k++)
      add_to_dW_3r2p ( specN, specN, specM2[k],   specN2, specM2[k],    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }     

  if (REACTION[17]){
    kf=_kf_Arrhenius(3, 6e16, -0.5, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(3, 6e16, -0.5, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;                
    for (k=0; specM3[k]!=specEND; k++)
      add_to_dW_3r2p ( specN, specO, specM3[k],   specNO, specM3[k],    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }     

  if (REACTION[18]){
    kf=_kf_Arrhenius(2, 1.5e13, 0.0, 0.0, T);
    dkfdT=_dkfdT_Arrhenius(2, 1.5e13, 0.0, 0.0, T);
    dkfdTv=0.0;
    dkfdTe=0.0;                
    add_to_dW_2r2p ( specN, specNO,   specO, specN2,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[19]){
    kf=_kf_Arrhenius(2, 1.8e8, 1.5, 3020.0*R, T);
    dkfdT=_dkfdT_Arrhenius(2, 1.8e8, 1.5, 3020.0*R, T);
    dkfdTv=0.0;
    dkfdTe=0.0;                
    add_to_dW_2r2p ( specN, specO2,   specO, specNO,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[20]){
    kf=_kf_Arrhenius(2, 2e19, -1.0, 0.0, Te);
    dkfdTe=_dkfdT_Arrhenius(2, 2e19, -1.0, 0.0, Te);
    dkfdT=0.0;
    dkfdTv=0.0;      
    add_to_dW_2r2p ( speceminus, specNOplus,   specN, specO,    kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  
}




void find_We_Lenard1964 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, 
                          double *We_create, double *We_destroy ) {
  long k;
  spec_t N;
  double R;

  R=1.9872;
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA;
  }
  *We_create=0.0;
  *We_destroy=0.0;
  if (REACTION[10])
    *We_create+=_kf_Arrhenius(2, 1.3e8, 1.0, 31900.0*R, Tv)*N[specN]*N[specO]/ calA * _calM(speceminus) * 1.0e6;
  if (REACTION[20])
    *We_destroy+=_kf_Arrhenius(2, 2e19, -1.0, 0.0, Te)*N[speceminus]*N[specNOplus]/ calA * _calM(speceminus) * 1.0e6;
}


void find_Qei_Lenard1964(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  *Qei=0.0;
}



void find_dQei_dx_Lenard1964(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;
}
