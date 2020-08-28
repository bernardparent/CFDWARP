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


const static bool REACTION[34]=
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
  };

#define specEND -1

const static long specM[]=
  {
   specH2, specH2O, specO2, specH, specO, specOH, specHO2, specH2O2, specN2, specN, specNO, specHNO, specNO2, specEND
  };

const static long eta1[]=
  {
   1.0, 6.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };

const static long eta2[]=
  {
   2.0, 6.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };

const static long eta3[]=
  {
   1.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };

const static long eta4[]=
  {
   2.0, 16.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };

const static long eta5[]=
  {
   1.0, 15.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };

const static long eta6[]=
  {
   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };



void find_W_Jachimowski1988 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;
  spec_t X;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
  }


  if (REACTION[1]) {
    add_to_W_fwbw_2r2p ( specH2, specO2,   specOH, specOH,  1.7e13, 0.0, 48000.0, T, X, W );
  }

  if (REACTION[2]) {
    add_to_W_fwbw_2r2p ( specH, specO2,   specOH, specO,  2.6e14, 0.0, 16800.0, T, X, W );
  }

  if (REACTION[3]) {
    add_to_W_fwbw_2r2p ( specO, specH2,   specOH, specH,  1.8e10, 1.0, 8900.0, T, X, W );
  }

  if (REACTION[4]) {
    add_to_W_fwbw_2r2p ( specOH, specH2,   specH2O, specH,  2.2e13, 0.0, 5150.0, T, X, W );
  }

  if (REACTION[5]) {
    add_to_W_fwbw_2r2p ( specOH, specOH,   specH2O, specO,  6.3e12, 0.0, 1090.0, T, X, W );
  }

  if (REACTION[6]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_3r2p ( specH, specOH, specM[k],   specH2O, specM[k], eta1[k]*2.2e22, -2.0, 0.0, T, X, W );
    }
  }

  if (REACTION[7]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_3r2p ( specH, specH, specM[k],   specH2, specM[k], eta2[k]*6.4e17, -1.0, 0.0, T, X, W );
    }
  }

  if (REACTION[8]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_3r2p ( specH, specO, specM[k],   specOH, specM[k], eta3[k]*6.0e16, -0.6, 0.0, T, X, W );
    }
  }

  if (REACTION[9]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_3r2p ( specH, specO2, specM[k],   specHO2, specM[k], eta4[k]*2.1e15, 0.0, -1000.0, T, X, W );
    }
  }

  if (REACTION[10]) {
    add_to_W_fwbw_2r2p ( specHO2, specH,   specH2, specO2,  1.3e13, 0.0, 0.0, T, X, W );
  }

  if (REACTION[11]) {
    add_to_W_fwbw_2r2p ( specHO2, specH,   specOH, specOH,  1.4e14, 0.0, 1080.0, T, X, W );
  }

  if (REACTION[12]) {
    add_to_W_fwbw_2r2p ( specHO2, specH,   specH2O, specO,  1.0e13, 0.0, 1080.0, T, X, W );
  }

  if (REACTION[13]) {
    add_to_W_fwbw_2r2p ( specHO2, specO,   specO2, specOH,  1.5e13, 0.0, 950.0, T, X, W );
  }

  if (REACTION[14]) {
    add_to_W_fwbw_2r2p ( specHO2, specOH,   specH2O, specO2,  8e12, 0.0, 0.0, T, X, W );
  }

  if (REACTION[15]) {
    add_to_W_fwbw_2r2p ( specHO2, specHO2,   specH2O2, specO2,  2e12, 0.0, 0.0, T, X, W );
  }

  if (REACTION[16]) {
    add_to_W_fwbw_2r2p ( specH, specH2O2,   specH2, specHO2,  1.4e12, 0.0, 3600.0, T, X, W );
  }

  if (REACTION[17]) {
    add_to_W_fwbw_2r2p ( specO, specH2O2,   specOH, specHO2,  1.4e13, 0.0, 6400.0, T, X, W );
  }

  if (REACTION[18]) {
    add_to_W_fwbw_2r2p ( specOH, specH2O2,   specH2O, specHO2,  6.1e12, 0.0, 1430.0, T, X, W );
  }

  if (REACTION[19]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_2r3p ( specM[k], specH2O2,   specOH, specOH, specM[k], eta5[k]*1.2e17, 0.0, 45500.0, T, X, W );
    }
  }

  if (REACTION[20]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_3r2p ( specO, specO, specM[k],    specO2,  specM[k], eta6[k]*6.0e17, 0.0, -1800.0, T, X, W );
    }
  }

  if (REACTION[21]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_3r2p ( specN, specN, specM[k],    specN2,  specM[k], eta6[k]*2.8e17, -0.75, 0.0, T, X, W );
    }
  }

  if (REACTION[22]) {
    add_to_W_fwbw_2r2p ( specN, specO2,   specNO, specO,  6.4e9, 1.0, 6300.0, T, X, W );
  }

  if (REACTION[23]) {
    add_to_W_fwbw_2r2p ( specN, specNO,   specN2, specO,  1.6e13, 0.0, 0.0, T, X, W );
  }

  if (REACTION[24]) {
    add_to_W_fwbw_2r2p ( specN, specOH,   specNO, specH,  6.3e11, 0.5, 0.0, T, X, W );
  }

  if (REACTION[25]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_3r2p ( specH, specNO, specM[k],    specHNO,  specM[k], eta6[k]*5.4e15, 0.0, -600.0, T, X, W );
    }
  }

  if (REACTION[26]) {
    add_to_W_fwbw_2r2p ( specH, specHNO,   specNO, specH2,  4.8e12, 0.0, 0.0, T, X, W );
  }

  if (REACTION[27]) {
    add_to_W_fwbw_2r2p ( specO, specHNO,   specNO, specOH,  5e11, 0.5, 0.0, T, X, W );
  }

  if (REACTION[28]) {
    add_to_W_fwbw_2r2p ( specOH, specHNO,   specNO, specH2O,  3.6e13, 0.0, 0.0, T, X, W );
  }

  if (REACTION[29]) {
    add_to_W_fwbw_2r2p ( specHO2, specHNO,   specNO, specH2O2,  2e12, 0.0, 0.0, T, X, W );
  }

  if (REACTION[30]) {
    add_to_W_fwbw_2r2p ( specHO2, specNO,   specNO2, specOH,  3.4e12, 0.0, -260.0, T, X, W );
  }

  if (REACTION[31]) {
    add_to_W_fwbw_2r2p ( specH, specNO2,   specNO, specOH,  3.5e14, 0.0, 1500.0, T, X, W );
  }

  if (REACTION[32]) {
    add_to_W_fwbw_2r2p ( specO, specNO2,   specNO, specO2,  1e13, 0.0, 600.0, T, X, W );
  }

  if (REACTION[33]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_W_fwbw_2r3p ( specM[k], specNO2,     specNO,  specO, specM[k], eta6[k]*1.16e16, 0.0, 66000.0, T, X, W );
    }
  }
}


void find_dW_dx_Jachimowski1988 ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
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





  if (REACTION[1]) {
    add_to_dW_fwbw_2r2p ( specH2, specO2,   specOH, specOH,  1.7e13, 0.0, 48000.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[2]) {
    add_to_dW_fwbw_2r2p ( specH, specO2,   specOH, specO,  2.6e14, 0.0, 16800.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[3]) {
    add_to_dW_fwbw_2r2p ( specO, specH2,   specOH, specH,  1.8e10, 1.0, 8900.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[4]) {
    add_to_dW_fwbw_2r2p ( specOH, specH2,   specH2O, specH,  2.2e13, 0.0, 5150.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[5]) {
    add_to_dW_fwbw_2r2p ( specOH, specOH,   specH2O, specO,  6.3e12, 0.0, 1090.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[6]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_3r2p ( specH, specOH, specM[k],   specH2O, specM[k], eta1[k]*2.2e22, -2.0, 0.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[7]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_3r2p ( specH, specH, specM[k],   specH2, specM[k], eta2[k]*6.4e17, -1.0, 0.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[8]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_3r2p ( specH, specO, specM[k],   specOH, specM[k], eta3[k]*6.0e16, -0.6, 0.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[9]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_3r2p ( specH, specO2, specM[k],   specHO2, specM[k], eta4[k]*2.1e15, 0.0, -1000.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[10]) {
    add_to_dW_fwbw_2r2p ( specHO2, specH,   specH2, specO2,  1.3e13, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[11]) {
    add_to_dW_fwbw_2r2p ( specHO2, specH,   specOH, specOH,  1.4e14, 0.0, 1080.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[12]) {
    add_to_dW_fwbw_2r2p ( specHO2, specH,   specH2O, specO,  1.0e13, 0.0, 1080.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[13]) {
    add_to_dW_fwbw_2r2p ( specHO2, specO,   specO2, specOH,  1.5e13, 0.0, 950.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[14]) {
    add_to_dW_fwbw_2r2p ( specHO2, specOH,   specH2O, specO2,  8e12, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[15]) {
    add_to_dW_fwbw_2r2p ( specHO2, specHO2,   specH2O2, specO2,  2e12, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[16]) {
    add_to_dW_fwbw_2r2p ( specH, specH2O2,   specH2, specHO2,  1.4e12, 0.0, 3600.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[17]) {
    add_to_dW_fwbw_2r2p ( specO, specH2O2,   specOH, specHO2,  1.4e13, 0.0, 6400.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[18]) {
    add_to_dW_fwbw_2r2p ( specOH, specH2O2,   specH2O, specHO2,  6.1e12, 0.0, 1430.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[19]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_2r3p ( specM[k], specH2O2,   specOH, specOH, specM[k], eta5[k]*1.2e17, 0.0, 45500.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[20]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_3r2p ( specO, specO, specM[k],    specO2,  specM[k], eta6[k]*6.0e17, 0.0, -1800.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[21]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_3r2p ( specN, specN, specM[k],    specN2,  specM[k], eta6[k]*2.8e17, -0.75, 0.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[22]) {
    add_to_dW_fwbw_2r2p ( specN, specO2,   specNO, specO,  6.4e9, 1.0, 6300.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[23]) {
    add_to_dW_fwbw_2r2p ( specN, specNO,   specN2, specO,  1.6e13, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[24]) {
    add_to_dW_fwbw_2r2p ( specN, specOH,   specNO, specH,  6.3e11, 0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[25]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_3r2p ( specH, specNO, specM[k],    specHNO,  specM[k], eta6[k]*5.4e15, 0.0, -600.0, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[26]) {
    add_to_dW_fwbw_2r2p ( specH, specHNO,   specNO, specH2,  4.8e12, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[27]) {
    add_to_dW_fwbw_2r2p ( specO, specHNO,   specNO, specOH,  5e11, 0.5, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[28]) {
    add_to_dW_fwbw_2r2p ( specOH, specHNO,   specNO, specH2O,  3.6e13, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[29]) {
    add_to_dW_fwbw_2r2p ( specHO2, specHNO,   specNO, specH2O2,  2e12, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[30]) {
    add_to_dW_fwbw_2r2p ( specHO2, specNO,   specNO2, specOH,  3.4e12, 0.0, -260.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[31]) {
    add_to_dW_fwbw_2r2p ( specH, specNO2,   specNO, specOH,  3.5e14, 0.0, 1500.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[32]) {
    add_to_dW_fwbw_2r2p ( specO, specNO2,   specNO, specO2,  1e13, 0.0, 600.0, T, X, dWdT, dWdrhok );
  }

  if (REACTION[33]) {
    for (k=0; specM[k]!=specEND; k++) {
      add_to_dW_fwbw_2r3p ( specM[k], specNO2,     specNO,  specO, specM[k], eta6[k]*1.16e16, 0.0, 66000.0, T, X, dWdT, dWdrhok );
    }
  }


}
