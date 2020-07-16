/*
SPDX-License-Identifier: BSD-2-Clause

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


#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
#include <model/share/chem_share.h> 


#define TOWNSEND TRUE
#define TOWNSEND_IMPLICIT FALSE
#define Estarmin 1e-40


/* set all reactions to true except for testing purposes */
const static bool REACTION[43]=
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
  };

#define specEND -1


void find_W (spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W) {
  double N[ns];
  double theta,kf;
  long k;
  double R;
  

  /* find properties needed by add_to_W* functions */
  R = 1.9872;
  for (k = 0; k < ns; k++) {
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1.0e-6 * calA;  /* particules/cm^3 */
  }
  Estar = max (Estarmin, Estar);
  theta = log (Estar);

  
}


/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx (spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam) {
  long k, s;                    
  spec_t N;
  double R,theta,kf,dkfdTe,dkfdT,dkfdTv;
  
  
  R=1.9872;
  /* initialize all derivatives to zero */
  for (s = 0; s < ns; s++) {
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


  
}



