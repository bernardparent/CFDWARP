/*
SPDX-License-Identifier: BSD-2-Clause

Copyright 2020 Aaron Trinh

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


#define TOWNSEND TRUE
#define TOWNSEND_IMPLICIT FALSE
#define Estarmin 1e-40


/* set all reactions to true except for testing purposes */
const static bool REACTION[44]=
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
   TRUE, /* reaction 29 */
   TRUE, /* reaction 30 */
   TRUE, /* reaction 31 */
   TRUE, /* reaction 32 */
   TRUE, /* reaction 33 */
   TRUE, /* reaction 34 */
   TRUE, /* reaction 35 */
   TRUE, /* reaction 36 */
   TRUE, /* reaction 37 */
   TRUE, /* reaction 38 */
   TRUE, /* reaction 39 */
   TRUE, /* reaction 40 */
   TRUE, /* reaction 41 */
   TRUE, /* reaction 42 */
   TRUE  /* reaction 43 */
  };

#define specEND -1

const static long specM1[]=
  {
    specAr, specC, specN, specO, specC2, specN2, specO2, specCN, specCO, specNO, specCO2, specNCO,
    specArplus, specCplus, specNplus, specOplus, specC2plus, specN2plus, specO2plus, specCNplus,
    specCOplus, specNOplus, specEND
  };
  
const static long specM2[]=
  {
    specAr, specC2, specN2, specO2, specCN, specCO, specNO, specCO2, specEND
  };
  
const static long specM3[]=
  {
    specC, specN, specO, specEND
  };

const static long specM4[]=
  {
    specAr, specC2, specN2, specO2, specCN, specCO, specEND
  };
  
const static long specM5[]=
  {
    specC, specN, specO, specNO, specCO2, specEND
  };
  
const static long specM6[]=
  {
    specC2, specN2, specO2, specCN, specCO, specNO, specCO2, specEND
  };



void write_model_chem_template(FILE **controlfile){
}


void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
}


void find_W (gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W) {
  double N[ns];
  double R,theta,kf,TTv,TvTe,TTe;
  long k;
  spec_t X;  

  /* find properties needed by add_to_W* functions */
  R = 1.9872;
  TTv=sqrt(Tv*T);
  TTe=sqrt(T*Te);
  TvTe=sqrt(Tv*Te);
  
  for (k = 0; k < ns; k++) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1.0e-6 * calA;  /* particules/cm^3 */
  }
  Estar = max (Estarmin, Estar);
  theta = log (Estar);

  
  if (REACTION[1])
    for (k = 0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specC2, specM1[k], specC, specC, specM1[k], 3.7e14, 0.0, 69900.0*R, TTv, X, W);
      add_to_W_bw_2r3p (specC2, specM1[k], specC, specC, specM1[k], 3.7e14, 0.0, 69900.0*R, T, X, W);   
    }
  
  if (REACTION[2])
    for (k = 0; specM2[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specN2, specM2[k], specN, specN, specM2[k], 7.0e21, -1.60, 69900.0*R, TTv, X, W);   
      add_to_W_bw_2r3p (specN2, specM2[k], specN, specN, specM2[k], 7.0e21, -1.60, 69900.0*R, T, X, W);
    }
  
  if (REACTION[3])    
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specN2, specM3[k], specN, specN, specM3[k], 3.0e22, -1.60, 113200.0*R, TTv, X, W); 
      add_to_W_bw_2r3p (specN2, specM3[k], specN, specN, specM3[k], 3.0e22, -1.60, 113200.0*R, T, X, W);  
    }
  
  if (REACTION[4]) {
    add_to_W_fw_2r3p (specN2, speceminus, specN, specN, speceminus, 1.2e25, -1.60, 113200.0*R, TvTe, X, W);   
    add_to_W_bw_2r3p (specN2, speceminus, specN, specN, speceminus, 1.2e25, -1.60, 113200.0*R, TTe, X, W); 
  }
  
  if (REACTION[5])
    for (k = 0; specM2[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specO2, specM2[k], specO, specO, specM2[k], 2.0e21, -1.50, 59750.0*R, TTv, X, W);  
      add_to_W_bw_2r3p (specO2, specM2[k], specO, specO, specM2[k], 2.0e21, -1.50, 59750.0*R, T, X, W);
    }
  
  if (REACTION[6])
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specO2, specM3[k], specO, specO, specM3[k], 1.0e22, -1.50, 59750.0*R, TTv, X, W);
      add_to_W_bw_2r3p (specO2, specM3[k], specO, specO, specM3[k], 1.0e22, -1.50, 59750.0*R, T, X, W);   
    }
  
  if (REACTION[7])
    for (k = 0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specCN, specM1[k], specC, specN, specM1[k], 2.5e14, 0.00, 71000.0*R, TTv, X, W);   
      add_to_W_bw_2r3p (specCN, specM1[k], specC, specN, specM1[k], 2.5e14, 0.00, 71000.0*R, T, X, W);  
    }
  
  if (REACTION[8]) {
    add_to_W_fw_2r3p (specCO, specAr, specC, specO, specAr, 2.3e19, -1.00, 129000.0*R, T, X, W);
    add_to_W_bw_2r3p (specCO, specAr, specC, specO, specAr, 2.3e19, -1.00, 129000.0*R, T, X, W);
  }

  if (REACTION[9])
    for (k = 0; specM2[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specCO, specM2[k], specC, specO, specM2[k], 3.4e20, -1.00, 129000.0*R, TTv, X, W);   
      add_to_W_bw_2r3p (specCO, specM2[k], specC, specO, specM2[k], 3.4e20, -1.00, 129000.0*R, T, X, W);
    }
  
  if (REACTION[10])
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specCO, specM3[k], specC, specO, specM3[k], 2.3e20, -1.00, 129000.0*R, TTv, X, W);  
      add_to_W_bw_2r3p (specCO, specM3[k], specC, specO, specM3[k], 2.3e20, -1.00, 129000.0*R, T, X, W);
    } 
  
  if (REACTION[11])
    for (k = 0; specM4[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specNO, specM4[k], specN, specO, specM4[k], 5.0e15, 0.00, 75500.0*R, TTv, X, W);   
      add_to_W_bw_2r3p (specNO, specM4[k], specN, specO, specM4[k], 5.0e15, 0.00, 75500.0*R, T, X, W);
    }
  
  if (REACTION[12])
    for (k = 0; specM5[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specNO, specM5[k], specN, specO, specM5[k], 1.1e17, 0.00, 75500.0*R, TTv, X, W);   
      add_to_W_bw_2r3p (specNO, specM5[k], specN, specO, specM5[k], 1.1e17, 0.00, 75500.0*R, T, X, W);
    }
  
  if (REACTION[13]) {
      add_to_W_fw_2r3p (specCO2, specAr, specCO, specO, specAr, 6.9e20, -1.50, 63275.0*R, T, X, W);
      add_to_W_bw_2r3p (specCO2, specAr, specCO, specO, specAr, 6.9e20, -1.50, 63275.0*R, T, X, W);
  }
  
  if (REACTION[14])
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specCO2, specM3[k], specCO, specO, specM3[k], 1.4e22, -1.00, 63275.0*R, TTv, X, W);  
      add_to_W_bw_2r3p (specCO2, specM3[k], specCO, specO, specM3[k], 1.4e22, -1.00, 63275.0*R, T, X, W);
    } 
  
  if (REACTION[15])
    for (k = 0; specM6[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specCO2, specM6[k], specCO, specO, specM6[k], 6.9e21, -1.00, 63275.0*R, TTv, X, W);   
      add_to_W_bw_2r3p (specCO2, specM6[k], specCO, specO, specM6[k], 6.9e21, -1.00, 63275.0*R, T, X, W); 
    }
  
  if (REACTION[16])
    for (k = 0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p (specNCO, specM1[k], specCO, specN, specM1[k], 6.3e16, -0.50, 24000.0*R, TTv, X, W);   
      add_to_W_bw_2r3p (specNCO, specM1[k], specCO, specN, specM1[k], 6.3e16, -0.50, 24000.0*R, T, X, W);
    }
    
  if (REACTION[17])
    add_to_W_fwbw_2r2p (specNO, specO, specN, specO2, 8.4e12, 0.00, 19450.0*R, T, X, W);

  if (REACTION[18])
    add_to_W_fwbw_2r2p (specN2, specO, specNO, specN, 6.4e17, -1.00, 38370.0*R, T, X, W);

  if (REACTION[19])
    add_to_W_fwbw_2r2p (specCO, specO, specC, specO2, 3.9e13, -0.18, 69200.0*R, T, X, W);
    
  if (REACTION[20])
    add_to_W_fwbw_2r2p (specCO, specC, specC2, specO, 2.0e17, -1.00, 58000.0*R, T, X, W);
    
  if (REACTION[21])
    add_to_W_fwbw_2r2p (specCO, specN, specCN, specO, 1.0e14, 0.00, 38600.0*R, T, X, W);

  if (REACTION[22])
    add_to_W_fwbw_2r2p (specN2, specC, specCN, specN, 1.1e14, -0.11, 23200.0*R, T, X, W);

  if (REACTION[23])
    add_to_W_fwbw_2r2p (specCN, specO, specNO, specC, 1.6e13, 0.10, 14600.0*R, T, X, W);
    
  if (REACTION[24])
    add_to_W_fwbw_2r2p (specCN, specC, specC2, specN, 5.0e13, 0.00, 13000.0*R, T, X, W);
    
  if (REACTION[25])
    add_to_W_fwbw_2r2p (specCO2, specO, specO2, specCO, 2.1e13, 0.00, 27800.0*R, T, X, W);

  if (REACTION[26])
    add_to_W_fwbw_2r2p (specCN, specO2, specNCO, specO, 6.6e12, 0.00, -200.0*R, T, X, W);

  if (REACTION[27])
    add_to_W_fwbw_2r2p (specCN, specCO2, specNCO, specCO, 4.0e14, 0.00, 19200.0*R, T, X, W);
    
  if (REACTION[28])
    add_to_W_fwbw_2r2p (specCN, specNO, specNCO, specN, 1.0e14, 0.00, 21200.0*R, T, X, W);
    
  if (REACTION[29])
    add_to_W_fwbw_2r2p (specCO, specNO, specNCO, specO, 3.8e17, -0.873, 51600.0*R, T, X, W);

  if (REACTION[30])
    add_to_W_fwbw_2r2p (specCN, specCO, specNCO, specC, 1.5e16, -0.487, 65800.0*R, T, X, W);
    
  if (REACTION[31]) {
    add_to_W_fwbw_2r2p (specN, specO, specNOplus, speceminus, 8.8e8, 1.00, 31900.0*R, T, X, W);
    add_to_W_fwbw_2r2p (specN, specO, specNOplus, speceminus, 8.8e8, 1.00, 31900.0*R, TvTe, X, W);
  }
    
  if (REACTION[32]) {
    add_to_W_fwbw_2r2p (specO, specO, specO2plus, speceminus, 7.1e2, 2.70, 80600.0*R, T, X, W);
    add_to_W_fwbw_2r2p (specO, specO, specO2plus, speceminus, 7.1e2, 2.70, 80600.0*R, TvTe, X, W);
  }
    
  if (REACTION[33]) {
    add_to_W_fwbw_2r2p (specC, specO, specCOplus, speceminus, 8.8e8, 1.00, 33100.0*R, T, X, W);
    add_to_W_fwbw_2r2p (specC, specO, specCOplus, speceminus, 8.8e8, 1.00, 33100.0*R, TvTe, X, W);
  }

  if (REACTION[34])
    add_to_W_fwbw_2r2p (specNOplus, specC, specNO, specCplus, 1.0e13, 0.00, 23200.0*R, T, X, W);

  if (REACTION[35])
    add_to_W_fwbw_2r2p (specO2plus, specO, specOplus, specO2, 4.0e12, -0.09, 18000.0*R, T, X, W);
    
  if (REACTION[36])
    add_to_W_fwbw_2r2p (specNOplus, specN, specOplus, specN2, 3.4e13, -1.08, 12800.0*R, T, X, W);
    
  if (REACTION[37])
    add_to_W_fwbw_2r2p (specNOplus, specO, specO2plus, specN, 7.2e12, 0.29, 48600.0*R, T, X, W);

  if (REACTION[38])
    add_to_W_fwbw_2r2p (specCO, specCplus, specCOplus, specC, 1.0e13, 0.00, 31400.0*R, T, X, W);

  if (REACTION[39])
    add_to_W_fwbw_2r2p (specO2, specCplus, specO2plus, specC, 1.0e13, 0.00, 9400.0*R, T, X, W);
    
  if (REACTION[40])
    add_to_W_fwbw_2r3p (specC, speceminus, specCplus, speceminus, speceminus, 3.9e33, -3.78, 130700.0*R, Te, X, W);
    
  if (REACTION[41])
    add_to_W_fwbw_2r3p (specO, speceminus, specOplus, speceminus, speceminus, 3.9e33, -3.78, 158500.0*R, Te, X, W);

  if (REACTION[42]) {
    kf=_kf_Arrhenius(2, 1.07e11, -0.52, 0.0, Te);
    add_to_W_2r1p (specOplus, speceminus, specO, kf, N, W);
  }
    
  if (REACTION[43]) {
    kf=_kf_Arrhenius(2, 2.02e11, -0.46, 0.0, Te);
    add_to_W_2r1p (specCplus, speceminus, specC, kf, N, W);
  }
  
}


/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx (gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam) {
  long k, s, spec;                    
  spec_t dWdTTv,dWdTTe,dWdTvTe;
  spec_t X;
  spec_t N;
  double R,theta,kf,TTv,TTe,TvTe,dkfdTe,dkfdT,dkfdTv,dTadT,dTadTv,Ta,dWdTa;
  
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  
  R=Rchem;
  TTv=sqrt(T*Tv);
  TTe=sqrt(T*Te);
  TvTe=sqrt(Tv*Te);
  /* initialize all derivatives to zero */
  for (s = 0; s < ns; s++) {
    dWdTTv[s] = 0.0;
    dWdTTe[s] = 0.0;
    dWdTvTe[s] = 0.0;
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for (k = 0; k < ns; k++) {
      dWdrhok[s][k] = 0.0;
    }
  }

  /* find properties needed by add_to_dW* functions in proper units */
  for (k = 0; k < ns; k++) {
    N[k] = rhok[k] / _calM ( k ) * 1.0e-6 * calA;
  }
  Estar = max (Estarmin, Estar);
  theta = log (Estar);

  if (REACTION[1])
    for (k = 0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specC2, specM1[k], specC, specC, specM1[k], 3.7e14, 0.0, 69900.0*R, TTv, X, dWdTTv, dWdrhok);   
      add_to_dW_bw_2r3p (specC2, specM1[k], specC, specC, specM1[k], 3.7e14, 0.0, 69900.0*R, T, X, dWdT, dWdrhok);
    }
  
  if (REACTION[2])
    for (k = 0; specM2[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specN2, specM2[k], specN, specN, specM2[k], 7.0e21, -1.60, 69900.0*R, TTv, X, dWdTTv, dWdrhok);   
      add_to_dW_bw_2r3p (specN2, specM2[k], specN, specN, specM2[k], 7.0e21, -1.60, 69900.0*R, T, X, dWdT, dWdrhok);   
    }
  
  if (REACTION[3])    
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specN2, specM3[k], specN, specN, specM3[k], 3.0e22, -1.60, 113200.0*R, TTv, X, dWdTTv, dWdrhok);   
      add_to_dW_bw_2r3p (specN2, specM3[k], specN, specN, specM3[k], 3.0e22, -1.60, 113200.0*R, T, X, dWdT, dWdrhok);   
    }
  
  if (REACTION[4]) {
    add_to_dW_fw_2r3p (specN2, speceminus, specN, specN, speceminus, 1.2e25, -1.60, 113200.0*R, TvTe, X, dWdTvTe, dWdrhok);   
    add_to_dW_bw_2r3p (specN2, speceminus, specN, specN, speceminus, 1.2e25, -1.60, 113200.0*R, TTe, X, dWdTTe, dWdrhok); 
  }
  
  if (REACTION[5])
    for (k = 0; specM2[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specO2, specM2[k], specO, specO, specM2[k], 2.0e21, -1.50, 59750.0*R, TTv, X, dWdTTv, dWdrhok);  
      add_to_dW_bw_2r3p (specO2, specM2[k], specO, specO, specM2[k], 2.0e21, -1.50, 59750.0*R, T, X, dWdT, dWdrhok); 
    }
  
  if (REACTION[6])
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specO2, specM3[k], specO, specO, specM3[k], 1.0e22, -1.50, 59750.0*R, TTv, X, dWdTTv, dWdrhok);  
      add_to_dW_bw_2r3p (specO2, specM3[k], specO, specO, specM3[k], 1.0e22, -1.50, 59750.0*R, T, X, dWdT, dWdrhok);  
    } 
  
  if (REACTION[7])
    for (k = 0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specCN, specM1[k], specC, specN, specM1[k], 2.5e14, 0.00, 71000.0*R, TTv, X, dWdTTv, dWdrhok); 
      add_to_dW_bw_2r3p (specCN, specM1[k], specC, specN, specM1[k], 2.5e14, 0.00, 71000.0*R, T, X, dWdT, dWdrhok); 
    }  
  
  if (REACTION[8])
    add_to_dW_fwbw_2r3p (specCO, specAr, specC, specO, specAr, 2.3e19, -1.00, 129000.0*R, T, X, dWdT, dWdrhok);
  
  if (REACTION[9])
    for (k = 0; specM2[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specCO, specM2[k], specC, specO, specM2[k], 3.4e20, -1.00, 129000.0*R, TTv, X, dWdTTv, dWdrhok);   
      add_to_dW_bw_2r3p (specCO, specM2[k], specC, specO, specM2[k], 3.4e20, -1.00, 129000.0*R, T, X, dWdT, dWdrhok);  
    }
  
  if (REACTION[10])
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specCO, specM3[k], specC, specO, specM3[k], 2.3e20, -1.00, 129000.0*R, TTv, X, dWdTTv, dWdrhok);   
      add_to_dW_bw_2r3p (specCO, specM3[k], specC, specO, specM3[k], 2.3e20, -1.00, 129000.0*R, T, X, dWdT, dWdrhok); 
    }
  
  if (REACTION[11])
    for (k = 0; specM4[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specNO, specM4[k], specN, specO, specM4[k], 5.0e15, 0.00, 75500.0*R, TTv, X, dWdTTv, dWdrhok);
      add_to_dW_bw_2r3p (specNO, specM4[k], specN, specO, specM4[k], 5.0e15, 0.00, 75500.0*R, T, X, dWdT, dWdrhok);
    }
  
  if (REACTION[12])
    for (k = 0; specM5[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specNO, specM5[k], specN, specO, specM5[k], 1.1e17, 0.00, 75500.0*R, TTv, X, dWdTTv, dWdrhok);
      add_to_dW_bw_2r3p (specNO, specM5[k], specN, specO, specM5[k], 1.1e17, 0.00, 75500.0*R, T, X, dWdT, dWdrhok);
    }
  
  if (REACTION[13])
    add_to_dW_fwbw_2r3p (specCO2, specAr, specCO, specO, specAr, 6.9e20, -1.50, 63275.0*R, T, X, dWdT, dWdrhok);
  
  if (REACTION[14])
    for (k = 0; specM3[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specCO2, specM3[k], specCO, specO, specM3[k], 1.4e22, -1.00, 63275.0*R, TTv, X, dWdTTv, dWdrhok);  
      add_to_dW_bw_2r3p (specCO2, specM3[k], specCO, specO, specM3[k], 1.4e22, -1.00, 63275.0*R, T, X, dWdT, dWdrhok); 
    } 
  
  if (REACTION[15])
    for (k = 0; specM6[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specCO2, specM6[k], specCO, specO, specM6[k], 6.9e21, -1.00, 63275.0*R, TTv, X, dWdTTv, dWdrhok);  
      add_to_dW_bw_2r3p (specCO2, specM6[k], specCO, specO, specM6[k], 6.9e21, -1.00, 63275.0*R, T, X, dWdT, dWdrhok);  
    } 
  
  if (REACTION[16])
    for (k = 0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p (specNCO, specM1[k], specCO, specN, specM1[k], 6.3e16, -0.50, 24000.0*R, TTv, X, dWdTTv, dWdrhok);
      add_to_dW_bw_2r3p (specNCO, specM1[k], specCO, specN, specM1[k], 6.3e16, -0.50, 24000.0*R, T, X, dWdT, dWdrhok);
    }
  
  if (REACTION[17])
    add_to_dW_fwbw_2r2p (specNO, specO, specN, specO2, 8.4e12, 0.00, 19450.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[18])
    add_to_dW_fwbw_2r2p (specN2, specO, specNO, specN, 6.4e17, -1.00, 38370.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[19])
    add_to_dW_fwbw_2r2p (specCO, specO, specC, specO2, 3.9e13, -0.18, 69200.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[20])
    add_to_dW_fwbw_2r2p (specCO, specC, specC2, specO, 2.0e17, -1.00, 58000.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[21])
    add_to_dW_fwbw_2r2p (specCO, specN, specCN, specO, 1.0e14, 0.00, 38600.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[22])
    add_to_dW_fwbw_2r2p (specN2, specC, specCN, specN, 1.1e14, -0.11, 23200.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[23])
    add_to_dW_fwbw_2r2p (specCN, specO, specNO, specC, 1.6e13, 0.10, 14600.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[24])
    add_to_dW_fwbw_2r2p (specCN, specC, specC2, specN, 5.0e13, 0.00, 13000.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[25])
    add_to_dW_fwbw_2r2p (specCO2, specO, specO2, specCO, 2.1e13, 0.00, 27800.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[26])
    add_to_dW_fwbw_2r2p (specCN, specO2, specNCO, specO, 6.6e12, 0.00, -200.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[27])
    add_to_dW_fwbw_2r2p (specCN, specCO2, specNCO, specCO, 4.0e14, 0.00, 19200.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[28])
    add_to_dW_fwbw_2r2p (specCN, specNO, specNCO, specN, 1.0e14, 0.00, 21200.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[29])
    add_to_dW_fwbw_2r2p (specCO, specNO, specNCO, specO, 3.8e17, -0.873, 51600.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[30])
    add_to_dW_fwbw_2r2p (specCN, specCO, specNCO, specC, 1.5e16, -0.487, 65800.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[31]) {
    add_to_dW_fw_2r2p (specN, specO, specNOplus, speceminus, 8.8e8, 1.00, 31900.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_bw_2r2p (specN, specO, specNOplus, speceminus, 8.8e8, 1.00, 31900.0*R, TvTe, X, dWdTvTe, dWdrhok);
  }
    
  if (REACTION[32]) {
    add_to_dW_fw_2r2p (specO, specO, specO2plus, speceminus, 7.1e2, 2.70, 80600.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_bw_2r2p (specO, specO, specO2plus, speceminus, 7.1e2, 2.70, 80600.0*R, TvTe, X, dWdTvTe, dWdrhok);
  }
    
  if (REACTION[33]) {
    add_to_dW_fw_2r2p (specC, specO, specCOplus, speceminus, 8.8e8, 1.00, 33100.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_bw_2r2p (specC, specO, specCOplus, speceminus, 8.8e8, 1.00, 33100.0*R, TvTe, X, dWdTvTe, dWdrhok);
  }

  if (REACTION[34])
    add_to_dW_fwbw_2r2p (specNOplus, specC, specNO, specCplus, 1.0e13, 0.00, 23200.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[35])
    add_to_dW_fwbw_2r2p (specO2plus, specO, specOplus, specO2, 4.0e12, -0.09, 18000.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[36])
    add_to_dW_fwbw_2r2p (specNOplus, specN, specOplus, specN2, 3.4e13, -1.08, 12800.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[37])
    add_to_dW_fwbw_2r2p (specNOplus, specO, specO2plus, specN, 7.2e12, 0.29, 48600.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[38])
    add_to_dW_fwbw_2r2p (specCO, specCplus, specCOplus, specC, 1.0e13, 0.00, 31400.0*R, T, X, dWdT, dWdrhok);

  if (REACTION[39])
    add_to_dW_fwbw_2r2p (specO2, specCplus, specO2plus, specC, 1.0e13, 0.00, 9400.0*R, T, X, dWdT, dWdrhok);
    
  if (REACTION[40])
    add_to_dW_fwbw_2r3p (specC, speceminus, specCplus, speceminus, speceminus, 3.9e33, -3.78, 130700.0*R, Te, X, dWdTe, dWdrhok);
    
  if (REACTION[41])
    add_to_dW_fwbw_2r3p (specO, speceminus, specOplus, speceminus, speceminus, 3.9e33, -3.78, 158500.0*R, Te, X, dWdTe, dWdrhok);

  if (REACTION[42]) {
    kf=_kf_Arrhenius(2, 1.07e11, -0.52, 0.0, Te);
    dkfdT=0.0;
    dkfdTv=0.0;
    dkfdTe=_dkfdT_Arrhenius(2, 1.07e11, -0.52, 0.0, Te);
    add_to_dW_2r1p ( specOplus, speceminus,   specO,   kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
    
  if (REACTION[43]) {
    kf=_kf_Arrhenius(2, 2.02e11, -0.46, 0.0, Te);
    dkfdT=0.0;
    dkfdTv=0.0;
    dkfdTe=_dkfdT_Arrhenius(2, 2.02e11, -0.46, 0.0, Te);
    add_to_dW_2r1p ( specCplus, speceminus,   specC,   kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
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


void find_Qei(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  *Qei=0.0;
}



void find_dQei_dx(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  

}


  



