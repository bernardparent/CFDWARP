// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2018 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/share/chem_share.h>

#define Tmin 500.0

void find_W ( spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;
  spec_t Gs, X;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    Gs[k] = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * _calM ( k );    /* J/mole */
    W[k] = 0.0;
  }

  if ( T > Tmin ) {

    // reaction #1
    add_to_W_Arrhenius2r3p ( specO2, specN, specO, specO, specN, 3.6E18, -1.0, 118800.0, T, X, Gs, W );
    // reaction #2
    add_to_W_Arrhenius2r3p ( specO2, specNO, specO, specO, specNO, 3.6E18, -1.0, 118800.0, T, X, Gs, W );
    // reaction #3 
    add_to_W_Arrhenius2r3p ( specN2, specO, specN, specN, specO, 1.9E17, -0.5, 226000.0, T, X, Gs, W );
    // reaction #4 
    add_to_W_Arrhenius2r3p ( specN2, specNO, specN, specN, specNO, 1.9E17, -0.5, 226000.0, T, X, Gs, W );
    // reaction #5 
    add_to_W_Arrhenius2r3p ( specN2, specO2, specN, specN, specO2, 1.9E17, -0.5, 226000.0, T, X, Gs, W );
    // reaction #6 
    add_to_W_Arrhenius2r3p ( specNO, specO2, specN, specO, specO2, 3.9E20, -1.5, 151000.0, T, X, Gs, W );
    // reaction #7 
    add_to_W_Arrhenius2r3p ( specNO, specN2, specN, specO, specN2, 3.9E20, -1.5, 151000.0, T, X, Gs, W );
    // reaction #8 
    add_to_W_Arrhenius2r2p ( specO, specNO, specN, specO2, 3.2E9, 1.0, 39400.0, T, X, Gs, W );
    // reaction #9 
    add_to_W_Arrhenius2r2p ( specO, specN2, specN, specNO, 7.0E13, 0.0, 76000.0, T, X, Gs, W );
    // reaction #10 
    add_to_W_Arrhenius2r3p ( specN, specN2, specN, specN, specN, 4.085E22, -1.5, 226000.0, T, X, Gs, W );

    // reaction #11 
    add_to_W_Arrhenius2r2p ( specO, specN, specNOplus, speceminus, 1.4E6, 1.5, 63800.0, T, X, Gs, W );
    // reaction #12 
    add_to_W_Arrhenius2r3p ( specO, speceminus,
                             specOplus, speceminus, speceminus, 3.6E31, -2.91, 316000.0, T, X, Gs, W );
    // reaction #13 
    add_to_W_Arrhenius2r3p ( specN, speceminus,
                             specNplus, speceminus, speceminus, 1.1E32, -3.14, 338000.0, T, X, Gs, W );
    // reaction #14 
    add_to_W_Arrhenius2r2p ( specO, specO, specO2plus, speceminus, 1.6E17, -0.98, 161600.0, T, X, Gs, W );
    // reaction #15 
    add_to_W_Arrhenius2r2p ( specO, specO2plus, specO2, specOplus, 2.92E18, -1.11, 56000.0, T, X, Gs, W );
    // reaction #16 
    add_to_W_Arrhenius2r2p ( specN2, specNplus, specN, specN2plus, 2.02e11, 0.81, 26000.0, T, X, Gs, W );
    // reaction #17 
    add_to_W_Arrhenius2r2p ( specN, specN, specN2plus, speceminus, 1.4e13, 0.0, 135600.0, T, X, Gs, W );
    // reaction #18 
    add_to_W_Arrhenius2r2p ( specO, specNOplus, specNO, specOplus, 3.63E15, -0.6, 101600.0, T, X, Gs, W );
    // reaction #19 
    add_to_W_Arrhenius2r2p ( specN2, specOplus, specO, specN2plus, 3.4E19, -2.0, 46000.0, T, X, Gs, W );
    // reaction #20 
    add_to_W_Arrhenius2r2p ( specN, specNOplus, specNO, specNplus, 1.0E19, -0.93, 122000.0, T, X, Gs, W );
    // reaction #21 
    add_to_W_Arrhenius2r2p ( specO2, specNOplus, specNO, specO2plus, 1.8E15, 0.17, 66000.0, T, X, Gs, W );
    // reaction #22 
    add_to_W_Arrhenius2r2p ( specO, specNOplus, specO2, specNplus, 1.34E13, 0.31, 154540.0, T, X, Gs, W );
    // reaction #23 
    add_to_W_Arrhenius2r3p ( specO2, specO, specO, specO, specO, 9E19, -1.0, 119000.0, T, X, Gs, W );

    // reaction #24 
    add_to_W_Arrhenius2r3p ( specO2, specO2, specO, specO, specO2, 3.24E19, -1.0, 119000.0, T, X, Gs, W );
    // reaction #25 
    add_to_W_Arrhenius2r3p ( specO2, specN2, specO, specO, specN2, 7.2E18, -1.0, 119000.0, T, X, Gs, W );
    // reaction #26 
    add_to_W_Arrhenius2r3p ( specN2, specN2, specN, specN, specN2, 4.7E17, -0.5, 226000.0, T, X, Gs, W );
    // reaction #27 
    add_to_W_Arrhenius2r3p ( specNO, specO, specN, specO, specO, 7.8E20, -1.5, 151000.0, T, X, Gs, W );
    // reaction #28 
    add_to_W_Arrhenius2r3p ( specNO, specN, specO, specN, specN, 7.8E20, -1.5, 151000.0, T, X, Gs, W );
    // reaction #29 
    add_to_W_Arrhenius2r3p ( specNO, specNO, specN, specO, specNO, 7.8E20, -1.5, 151000.0, T, X, Gs, W );
    // reaction #30 
    add_to_W_Arrhenius2r3p ( specO2, specN2,
                             specNO, specNOplus, speceminus, 1.38E20, -1.84, 282000.0, T, X, Gs, W );
    // reaction #31 
    add_to_W_Arrhenius2r3p ( specNO, specN2,
                             specNOplus, speceminus, specN2, 2.2E15, -0.35, 216000.0, T, X, Gs, W );
  }
}

void find_dW_dx ( spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */
  spec_t Gs, X, dGsdT;
  spec2_t dWdX;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    Gs[k] = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * _calM ( k );    /* J/mole */
    dGsdT[k] = ( _cpk_from_T ( k, T ) - _sk_from_T ( k, T ) - T * _dsk_dT_from_T ( k, T )
       ) * _calM ( k );
  }

  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
      dWdX[s][k] = 0.0;
    }
  }

  if ( T > Tmin ) {

    // reaction #1 
    add_to_dW_Arrhenius2r3p ( specO2, specN,
                              specO, specO, specN, 3.6E18, -1.0, 118800.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #2 
    add_to_dW_Arrhenius2r3p ( specO2, specNO,
                              specO, specO, specNO, 3.6E18, -1.0, 118800.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #3 
    add_to_dW_Arrhenius2r3p ( specN2, specO,
                              specN, specN, specO, 1.9E17, -0.5, 226000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #4  
    add_to_dW_Arrhenius2r3p ( specN2, specNO,
                              specN, specN, specNO, 1.9E17, -0.5, 226000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #5 
    add_to_dW_Arrhenius2r3p ( specN2, specO2,
                              specN, specN, specO2, 1.9E17, -0.5, 226000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #6 
    add_to_dW_Arrhenius2r3p ( specNO, specO2,
                              specN, specO, specO2, 3.9E20, -1.5, 151000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #7 
    add_to_dW_Arrhenius2r3p ( specNO, specN2,
                              specN, specO, specN2, 3.9E20, -1.5, 151000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #8 
    add_to_dW_Arrhenius2r2p ( specO, specNO,
                              specN, specO2, 3.2E9, 1.0, 39400.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #9 
    add_to_dW_Arrhenius2r2p ( specO, specN2,
                              specN, specNO, 7.0E13, 0.0, 76000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #10  
    add_to_dW_Arrhenius2r3p ( specN, specN2,
                              specN, specN, specN, 4.085E22, -1.5, 226000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #11 
    add_to_dW_Arrhenius2r2p ( specO, specN,
                              specNOplus, speceminus, 1.4E6, 1.5, 63800.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #12 
    add_to_dW_Arrhenius2r3p ( specO, speceminus,
                              specOplus, speceminus, speceminus,
                              3.6E31, -2.91, 316000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #13 
    add_to_dW_Arrhenius2r3p ( specN, speceminus,
                              specNplus, speceminus, speceminus,
                              1.1E32, -3.14, 338000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #14 
    add_to_dW_Arrhenius2r2p ( specO, specO,
                              specO2plus, speceminus, 1.6E17, -0.98, 161600.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #15 
    add_to_dW_Arrhenius2r2p ( specO, specO2plus,
                              specO2, specOplus, 2.92E18, -1.11, 56000.0, T, X, Gs, dGsdT, dWdT, dWdX );

    // reaction #16 
    add_to_dW_Arrhenius2r2p ( specN2, specNplus,
                              specN, specN2plus, 2.02e11, 0.81, 26000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #17 
    add_to_dW_Arrhenius2r2p ( specN, specN,
                              specN2plus, speceminus, 1.4e13, 0.0, 135600.0, T, X, Gs, dGsdT, dWdT, dWdX );

    // reaction #18 
    add_to_dW_Arrhenius2r2p ( specO, specNOplus,
                              specNO, specOplus, 3.63E15, -0.6, 101600.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #19 
    add_to_dW_Arrhenius2r2p ( specN2, specOplus,
                              specO, specN2plus, 3.4E19, -2.0, 46000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #20 
    add_to_dW_Arrhenius2r2p ( specN, specNOplus,
                              specNO, specNplus, 1.0E19, -0.93, 122000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #21 
    add_to_dW_Arrhenius2r2p ( specO2, specNOplus,
                              specNO, specO2plus, 1.8E15, 0.17, 66000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #22 
    add_to_dW_Arrhenius2r2p ( specO, specNOplus,
                              specO2, specNplus, 1.34E13, 0.31, 154540.0, T, X, Gs, dGsdT, dWdT, dWdX );

    // reaction #23 
    add_to_dW_Arrhenius2r3p ( specO2, specO,
                              specO, specO, specO, 9E19, -1.0, 119000.0, T, X, Gs, dGsdT, dWdT, dWdX );

    // reaction #24 
    add_to_dW_Arrhenius2r3p ( specO2, specO2,
                              specO, specO, specO2, 3.24E19, -1.0, 119000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #25 
    add_to_dW_Arrhenius2r3p ( specO2, specN2,
                              specO, specO, specN2, 7.2E18, -1.0, 119000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #26 
    add_to_dW_Arrhenius2r3p ( specN2, specN2,
                              specN, specN, specN2, 4.7E17, -0.5, 226000.0, T, X, Gs, dGsdT, dWdT, dWdX );

    // reaction #27 
    add_to_dW_Arrhenius2r3p ( specNO, specO,
                              specN, specO, specO, 7.8E20, -1.5, 151000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #28 
    add_to_dW_Arrhenius2r3p ( specNO, specN,
                              specO, specN, specN, 7.8E20, -1.5, 151000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #29 
    add_to_dW_Arrhenius2r3p ( specNO, specNO,
                              specN, specO, specNO, 7.8E20, -1.5, 151000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #30 
    add_to_dW_Arrhenius2r3p ( specO2, specN2,
                              specNO, specNOplus, speceminus,
                              1.38E20, -1.84, 282000.0, T, X, Gs, dGsdT, dWdT, dWdX );
    // reaction #31 
    add_to_dW_Arrhenius2r3p ( specNO, specN2,
                              specNOplus, speceminus, specN2,
                              2.2E15, -0.35, 216000.0, T, X, Gs, dGsdT, dWdT, dWdX );

    for ( s = 0; s < ns; s++ ) {
      for ( k = 0; k < ns; k++ ) {
        dWdrhok[s][k] = dWdX[s][k] / _calM ( k ) * 1.0e-06;
      }
    }

  }
}
