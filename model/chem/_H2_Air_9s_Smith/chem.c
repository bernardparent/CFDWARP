// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2020 Ajjay Omprakas

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
#include <model/share/chem_share.h>




const static bool REACTION[28]=
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
  };
  
  #define specEND -1
  
  const static long specM1[]=
  {
    specO2, specH2O, specN2, specEND
  };
  
  const static long specM2[]=
  {
     specH2, specH2O, specN2, specEND
  };

void write_model_chem_template(FILE **controlfile){
}


void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
}


void find_W ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;
  spec_t X;
   
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
  }


  if (REACTION[1])
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_W_fwbw_3r2p ( specO, specO, specM2[k], specO2, specM2[k], 1.2E17, -1.0, 0.0, T, X, W);
  }
  
  if (REACTION[2])
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_W_fwbw_3r2p ( specO, specH, specM2[k], specOH, specM2[k], 5.0E17, -1.0, 0.0, T, X, W);
  }
  
  if (REACTION[3])
    add_to_W_fwbw_2r2p ( specO, specH2, specOH, specH, 3.87E04, 2.7, 6260.0, T, X, W);

  if (REACTION[4])
    add_to_W_fwbw_2r2p ( specO, specHO2, specOH, specO2, 2.0E13, 0.0, 0.0, T, X, W);

  if (REACTION[5])
    add_to_W_fwbw_2r2p ( specO, specH2O2, specOH, specHO2, 9.63E6, 2.0, 4000.0, T, X, W);

  if (REACTION[6])
  {
    for(k=0; specM1[k]!=specEND ; k++)
     add_to_W_fwbw_3r2p ( specH, specO2, specM1[k], specHO2, specM1[k], 2.8E18, -0.862, 0.0, T, X, W);
  }
  
  if (REACTION[7])
    add_to_W_fwbw_3r2p ( specH, specO2, specO2, specHO2, specO2, 2.8E19, -1.24, 0.0, T, X, W);

  if (REACTION[8])
    add_to_W_fwbw_3r2p ( specH, specO2, specH2O, specHO2, specH2O, 11.26E18, -0.76, 0.0, T, X, W);

  if (REACTION[9])
    add_to_W_fwbw_3r2p ( specH, specO2, specO2, specHO2, specO2, 2.6E19, -1.24, 0.0, T, X, W);

  if (REACTION[10])
    add_to_W_fwbw_2r2p ( specH, specO2, specOH, specO, 2.65E16, -0.6707, 17041.0, T, X, W);

  if (REACTION[11])
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_W_fwbw_3r2p ( specH, specH, specM2[k], specH2, specM2[k], 1.0E18, -1.0, 0.0, T, X, W);
  }
  if (REACTION[12])
    add_to_W_fwbw_3r2p ( specH, specH,
                             specH2, specH2, specH2, 9.0E16, -0.6, 0.0, T, X, W);

  if (REACTION[13])
    add_to_W_fwbw_3r2p ( specH, specH,
                             specH2O, specH2, specH2O, 6.0E19, -1.25, 0.0, T, X, W);

  if (REACTION[14])
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_W_fwbw_3r2p ( specH, specOH, specM2[k], specH2O, specM2[k], 2.2E22, -2.0, 0.0, T, X, W);
  }
  
  if (REACTION[15])
    add_to_W_fwbw_2r2p ( specH, specHO2, specH2O, specO, 3.97E12, 0.0, 671.0, T, X, W);

  if (REACTION[16])
    add_to_W_fwbw_2r2p ( specH, specHO2, specOH, specOH, 0.84e14, 0.0, 635.0, T, X, W );

  if (REACTION[17])
    add_to_W_fwbw_2r2p ( specH, specH2O2, specHO2, specH2, 1.2e07, 2.0, 5200.0, T, X, W );

  if (REACTION[18])
    add_to_W_fwbw_2r2p ( specH, specH2O2, specOH, specH2O, 1.0E13, 0.0, 3600.0, T, X, W);

  if (REACTION[19])
    add_to_W_fwbw_2r2p ( specOH, specH2, specH, specH2O, 2.16E08, 1.510, 3430.0, T, X, W );

  if (REACTION[20])
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_W_fwbw_3r2p ( specOH, specOH, specM2[k], specH2O2, specM2[k], 7.4E13, -0.371, 0.0, T, X, W);
  }
  
  if (REACTION[21])
    add_to_W_fwbw_2r2p ( specOH, specOH, specH2O, specO, 3.57E04, 2.4, -2100.0, T, X, W);

  if (REACTION[22])
    add_to_W_fwbw_2r2p ( specOH, specHO2, specH2O, specO2, 1.45E13, 0.0, -500.0, T, X, W);

  if (REACTION[23])
    add_to_W_fwbw_2r2p ( specOH, specH2O2, specH2O, specHO2, 2.0E12, 0.0, 427.0, T, X, W);

  if (REACTION[24])
    add_to_W_fwbw_2r2p ( specOH, specH2O2, specH2O, specHO2, 1.7E18, 0.0, 29410.0, T, X, W );

  if (REACTION[25])
    add_to_W_fwbw_2r2p ( specHO2, specHO2, specO2, specH2O2, 1.3E11, 0.0, -1630.0, T, X, W);

  if (REACTION[26])
    add_to_W_fwbw_2r2p ( specHO2, specHO2, specO2, specH2O2, 4.2E14, 0.0, 12000.0, T, X, W);

  if (REACTION[27])
    add_to_W_fwbw_2r2p ( specOH, specHO2, specO2, specH2O, 0.5E16, 0.0, 17330.0, T, X, W);

    
}

void find_dW_dx ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */
  spec_t X;

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


  if (REACTION[1]) 
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_dW_fwbw_3r2p ( specO, specO, specM2[k], specO2, specM2[k], 1.2E17, -1.0, 0.0, T, X, dWdT, dWdrhok );
  }
  
  if (REACTION[2]) 
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_dW_fwbw_3r2p ( specO, specH, specM2[k], specOH, specM2[k], 5.0E17, -1.0, 0.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[3]) 
    add_to_dW_fwbw_2r2p ( specO, specH2, specOH, specH, 3.87E04, 2.7, 6260.0, T, X, dWdT, dWdrhok );

  if (REACTION[4])  
    add_to_dW_fwbw_2r2p ( specO, specHO2, specOH, specO2, 2.0E13, 0.0, 0.0, T, X, dWdT, dWdrhok );

  if (REACTION[5]) 
    add_to_dW_fwbw_2r2p ( specO, specH2O2, specOH, specHO2, 9.63E6, 2.0, 4000.0, T, X, dWdT, dWdrhok );

  if (REACTION[6])
  {
    for(k=0; specM1[k]!=specEND ; k++) 
     add_to_dW_fwbw_3r2p ( specH, specO2, specM1[k], specHO2, specM1[k], 2.8E18, -0.862, 0.0, T, X, dWdT, dWdrhok );
  }
  
  if (REACTION[7]) 
    add_to_dW_fwbw_3r2p ( specH, specO2, specO2, specHO2, specO2, 2.8E19, -1.24, 0.0, T, X, dWdT, dWdrhok );

  if (REACTION[8]) 
    add_to_dW_fwbw_3r2p ( specH, specO2, specH2O, specHO2, specH2O, 11.26E18, -0.76, 0.0, T, X, dWdT, dWdrhok );

  if (REACTION[9]) 
    add_to_dW_fwbw_3r2p ( specH, specO2, specO2, specHO2, specO2, 2.6E19, -1.24, 0.0, T, X, dWdT, dWdrhok );

  if (REACTION[10])  
    add_to_dW_fwbw_2r2p ( specH, specO2, specOH, specO, 2.65E16, -0.6707, 17041.0, T, X, dWdT, dWdrhok );

  if (REACTION[11]) 
  {
    for(k=0; specM2[k]!=specEND ; k++) 
     add_to_dW_fwbw_3r2p ( specH, specH, specM2[k], specH2, specM2[k], 1.0E18, -1.0, 0.0, T, X, dWdT, dWdrhok );
  }
  if (REACTION[12]) 
    add_to_dW_fwbw_3r2p ( specH, specH,
                             specH2, specH2, specH2, 9.0E16, -0.6, 0.0, T, X, dWdT, dWdrhok );

  if (REACTION[13]) 
    add_to_dW_fwbw_3r2p ( specH, specH,
                             specH2O, specH2, specH2O, 6.0E19, -1.25, 0.0, T, X, dWdT, dWdrhok );

  if (REACTION[14])
  {
    for(k=0; specM2[k]!=specEND ; k++)  
    add_to_dW_fwbw_3r2p ( specH, specOH, specM2[k], specH2O, specM2[k], 2.2E22, -2.0, 0.0, T, X, dWdT, dWdrhok );
  }
  
  if (REACTION[15]) 
    add_to_dW_fwbw_2r2p ( specH, specHO2, specH2O, specO, 3.97E12, 0.0, 671.0, T, X, dWdT, dWdrhok );

  if (REACTION[16]) 
    add_to_dW_fwbw_2r2p ( specH, specHO2, specOH, specOH, 0.84e14, 0.0, 635.0, T, X, dWdT, dWdrhok );

  if (REACTION[17]) 
    add_to_dW_fwbw_2r2p ( specH, specH2O2, specHO2, specH2, 1.2e07, 2.0, 5200.0, T, X, dWdT, dWdrhok );

  if (REACTION[18]) 
    add_to_dW_fwbw_2r2p ( specH, specH2O2, specOH, specH2O, 1.0E13, 0.0, 3600.0, T, X, dWdT, dWdrhok );

  if (REACTION[19]) 
    add_to_dW_fwbw_2r2p ( specOH, specH2, specH, specH2O, 2.16E08, 1.510, 3430.0, T, X, dWdT, dWdrhok );

  if (REACTION[20]) 
  {
    for(k=0; specM2[k]!=specEND ; k++)
     add_to_dW_fwbw_3r2p ( specOH, specOH, specM2[k], specH2O2, specM2[k], 7.4E13, -0.371, 0.0, T, X, dWdT, dWdrhok );
  }
  if (REACTION[21]) 
    add_to_dW_fwbw_2r2p ( specOH, specOH, specH2O, specO, 3.57E04, 2.4, -2100.0, T, X, dWdT, dWdrhok );

  if (REACTION[22]) 
    add_to_dW_fwbw_2r2p ( specOH, specHO2, specH2O, specO2, 1.45E13, 0.0, -500.0, T, X, dWdT, dWdrhok );

  if (REACTION[23]) 
    add_to_dW_fwbw_2r2p ( specOH, specH2O2, specH2O, specHO2, 2.0E12, 0.0, 427.0, T, X, dWdT, dWdrhok );

  if (REACTION[24]) 
    add_to_dW_fwbw_2r2p ( specOH, specH2O2, specH2O, specHO2, 1.7E18, 0.0, 29410.0, T, X, dWdT, dWdrhok );

  if (REACTION[25]) 
    add_to_dW_fwbw_2r2p ( specHO2, specHO2, specO2, specH2O2, 1.3E11, 0.0, -1630.0, T, X, dWdT, dWdrhok );

  if (REACTION[26]) 
    add_to_dW_fwbw_2r2p ( specHO2, specHO2, specO2, specH2O2, 4.2E14, 0.0, 12000.0, T, X, dWdT, dWdrhok );

  if (REACTION[27]) 
    add_to_dW_fwbw_2r2p ( specOH, specHO2, specO2, specH2O, 0.5E16, 0.0, 17330.0, T, X,  dWdT, dWdrhok );

}




void find_Qei(spec_t rhok, double Estar, double Te, double *Qei){  
  *Qei=0.0;
}



void find_dQei_dx(spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
}
