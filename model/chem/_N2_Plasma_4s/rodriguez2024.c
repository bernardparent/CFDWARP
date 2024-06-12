// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2021 Prasanna Thoguluva Rajendran
Copyright 2024 Felipe Martin Rodriguez Fuentes

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

#define Estarmin 1e-40

/* set all reactions to true except for testing purposes */
const static bool REACTION[6]=
  {
   TRUE, /* reaction 0 */
   TRUE, /* reaction 1 */
   TRUE, /* reaction 2 */
   TRUE, /* reaction 3 */
   TRUE, /* reaction 4 */
   TRUE, /* reaction 5 */
  };

#define specEND -1

const static long specM1[]=
  {
   specN2,  specEND
  };

const static long specM2[]=
  {
   specN2,  specN,  specEND
  };


static double _kf1(double Te){
  /* N2 ionization */
  double kf1;
  double Te_control[] = 
  { 
    9.46979690785538,
    9.51956710917422,
    9.84281034346154,
    10.1366387167620,
    10.7022762534868,
    11.6533000019413,
    11.8041027220493,
    12.6902842369606,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -41.1651192144546,
    -37.0450701265803,
    -30.9051454150272,
    -28.6391530514683,
    -25.2959688327970,
    -19.9307259125084,
    -19.3648566793161,
    -17.2892786324613,
    -15.7126305428502
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf1 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf1 ));
}


static double _kf2(double Te){
  /* N2+ recombination */
  double kf2;
  double Te_control[] = 
  { 
    2.96134635039148,
    5.71525385948086,
    9.17412056693579,
    9.86376375536634,
    12.3695988451580,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -15.1096633905206,
    -15.6250837536637,
    -16.9584559651989,
    -17.4850816245454,
    -18.7960420527258,
    -20.7232658369464
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf2 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf2 ));
}



void find_W_Rodriguez2024 (gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                       double Qbeam, spec_t W ) {
  double tmp,Nt;
  long k;
  spec_t N;

  Nt = 0.0;
  for ( k = 0; k < ns; k++ ) {
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
    Nt += N[k];
  }
  
  Estar = max ( Estarmin, Estar );
  

  if (REACTION[1]){
    add_to_W_2r3p ( specN2, speceminus,   specN2plus, speceminus, speceminus, _kf1(Te), N, W );
  }
  

  if (REACTION[2]){
    add_to_W_2r2p ( speceminus, specN2plus, specN, specN, _kf2(Te), N, W );
  }
  
  
  if (REACTION[3]){
    tmp = 1.8e11 * Qbeam / Nt * N[specN2];
    W[specN2] += -tmp * _calM(specN2) / calA * 1.0e6;
    W[specN2plus] += +tmp * _calM(specN2plus) / calA * 1.0e6;
    W[speceminus] += +tmp * _calM(speceminus) / calA * 1.0e6;
  }
  

  tmp = exp ( -113200.0 / T ) * ( 1.0 - exp ( -3354.0 / T ) );
  if (REACTION[4]){
    for (k=0; specM1[k]!=specEND; k++){
      add_to_W_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 5.0e-8 * tmp, N, W);
    }
  }


  tmp = 8.27e-34 * exp ( 500.0 / T );
  if (REACTION[5]){
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

void find_dW_dx_Rodriguez2024 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    
  spec_t N;
  double Nt,kf,tmp,dkfdTe,dkfdT,dkfdTv,dtmpdT,dkfdQb;
  //double Estar_from_Te,Te_from_Estar;

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
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA; /* particules/cm^3 */
    Nt += N[k];
  }
                  
    
  /*find_EoverN_from_Te(Te, &Estar_from_Te);
  //if (Estar_from_Te>Estar) Estar=Estar_from_Te;   //??? needs to be verified
  find_Te_from_EoverN(Estar, &Te_from_Estar);
  Te_from_Estar=max(300.0,Te_from_Estar); 
  Estar = max ( Estarmin, Estar );
  theta = log ( Estar );
  */

  if (REACTION[1]){
    dkfdTe = 0.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specN2, speceminus, specN2plus, speceminus, speceminus, _kf1(Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  

  if (REACTION[2]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r2p ( specN2plus, speceminus, specN, specN, _kf2(Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }
  

  
  if (REACTION[3]){
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
  


  tmp = exp ( -113200.0 / T ) * ( 1.0 - exp ( -3354.0 / T ) );
  dtmpdT = 113200.0 / sqr ( T ) * tmp - exp ( -113200.0 / T ) * 3354.0 / sqr ( T ) * exp ( -3354.0 / T );
  if (REACTION[4]){
    kf=5.0e-8 * tmp;
    dkfdTe = 0.0;
    dkfdT= 5.0e-8 * dtmpdT;
    dkfdTv=0.0;
    for (k=0; specM1[k]!=specEND; k++){
      add_to_dW_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }


  tmp = 8.27e-34 * exp ( 500.0 / T );
  if (REACTION[5]){
    kf=tmp;
    dkfdTe = 0.0;
    dkfdT= -8.27e-34 * 500.0 / sqr ( T ) * exp ( 500.0 / T );
    dkfdTv=0.0;
    for (k=0; specM2[k]!=specEND; k++){
      add_to_dW_3r2p ( specN, specN, specM2[k],   specN2, specM2[k], kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

}



void find_Qei_Rodriguez2024( gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  
  *Qei=0.0;  

  if (REACTION[1])
    add_to_Qei(gl,Te,specN2,_ionizationpot(specN2), _kf1(Te), rhok, Qei);
}



void find_dQei_dx_Rodriguez2024( gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  

  if (REACTION[1])
    add_to_dQei(gl,Te,specN2,_ionizationpot(specN2), _kf1(Te), 0.0,  rhok, dQeidrhok, dQeidTe);
}
