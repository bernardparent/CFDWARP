// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2008, 2018 Bernard Parent

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

/* taken from Kundu, K.P., Penko, P.F., VanOverbeke, T.J., "A Practical Mechanism for Computing Combustion in Gas Turbine Engines," 35th Joint Propulsion Conference and Exhibit, AIAA Paper 99-2218, 1999. */

#define Xmin 1e-40
#define Tmin 750.0              //T needs to be clipped to prevent instabilities (WARNING: SHOULD LOOK INTO THIS)

void find_W ( spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;
  spec_t X, Gs;
  double brr, frr;

  for ( k = 0; k < ns; k++ )
    W[k] = 0.0;

  if (T>Tmin){
  for ( k = 0; k < ns; k++ ) {
    X[k] = max ( -1e6, rhok[k] / _calM ( k ) * 1.0e-06 );       /* moles per cm3 */
    Gs[k] = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * _calM ( k );    /* J/moles */
  }

  X[specH2] = max ( Xmin, X[specH2] );
  X[specN2] = max ( Xmin, X[specN2] );
  X[specC3H8] = max ( Xmin, X[specC3H8] );
  X[specO] = max ( Xmin, X[specO] );


//      Rgasconstant=1.987192004e0 cal/mol/K; 
// forward reaction 1: N2+C3H8 -> 3CH + 5H + N2
  frr =
    4.50e+10 * exp ( -30000.0e0 / ( 1.987192004e0 * T ) ) * pow ( X[specN2], 0.8 ) * pow ( X[specC3H8], 0.8 );
  W[specCH] += +3.0 * frr;
  W[specH] += +5.0 * frr;
  W[specC3H8] -= +1.0 * frr;

// forward reaction 1: N2+C12H23 -> 12CH + 11H + N2
//      frr  = 5.50e+10 * exp(-30000.0e0/(1.987192004e0*T))*pow(X[specN2],0.8)*pow(X[specC12H23],0.8);
//      W[specCH] += +12.0*frr; 
//      W[specH] += +11.0*frr; 
//      W[specC3H8]-= +1.0*frr; 

// forward reaction 2: CH + H2 + N2 -> 2NH + CH
  frr =
    1.00e+16 * exp ( -78000.0e0 / ( 1.987192004e0 * T ) ) * sqr ( X[specCH] ) * pow ( X[specH2],
                                                                                      0.1 ) * X[specN2];
  W[specNH] += +2.0 * frr;
  W[specH2] -= +1.0 * frr;
  W[specN2] -= +1.0 * frr;

// backward reaction 2: CH + 2NH -> H2 + N2 + CH
  brr = 1.95e+15 * exp ( -53900.0e0 / ( 1.987192004e0 * T ) ) * X[specCH] * sqr ( X[specNH] );
  W[specH2] += +1.0 * brr;
  W[specN2] += +1.0 * brr;
  W[specNH] -= +2.0 * brr;

// forward reaction 3: O + N2 + HO2 -> 2NO + H + O
  frr =
    1.95e+08 * exp ( -41900.0e0 / ( 1.987192004e0 * T ) ) * sqrt ( T ) * pow ( X[specO],
                                                                               0.1 ) * sqrt ( X[specN2] ) *
    X[specHO2];
  W[specNO] += 2.0 * frr;
  W[specH] += 1.0 * frr;
  W[specN2] -= 1.0 * frr;
  W[specHO2] -= 1.0 * frr;

// backward reaction 3: 2NO + H -> N2 + HO2
  brr = 1.25e+10 * exp ( -6000.0e0 / ( 1.987192004e0 * T ) ) * X[specNO] * X[specH];
  W[specN2] += 1.0 * brr;
  W[specHO2] += 1.0 * brr;
  W[specNO] -= 2.0 * brr;
  W[specH] -= 1.0 * brr;

// forward reaction 4: N2 + CO + HO2 -> CO2 + OH + N2
  frr =
    3.50e+13 * exp ( -22934.0e0 / ( 1.987192004e0 * T ) ) * pow ( X[specN2], 0.1 ) * X[specCO] * X[specHO2];
  W[specCO2] += 1.0 * frr;
  W[specOH] += 1.0 * frr;
  W[specCO] -= 1.0 * frr;
  W[specHO2] -= 1.0 * frr;

// forward reaction 5: N2 + O2 -> 2O + N2
  frr = 1.00e+18 * exp ( -122239.0e0 / ( 1.987192004e0 * T ) ) * X[specN2] * X[specO2];
  W[specO] += 2.0 * frr;
  W[specO2] -= 1.0 * frr;

// backward reaction 5: H2 + 2O -> O2 + H2
  brr = 1.00e+18 * X[specH2] * sqr ( X[specO] );
  W[specO] -= 2.0 * brr;
  W[specO2] += 1.0 * brr;

  for ( k = 0; k < ns; k++ )
    W[k] *= 1.0e+06 * _calM ( k );

  // reaction #6: H2 + OH -> H2O + H
  add_to_W_fwbw_2r2p ( specH2, specOH, specH2O, specH, 1.17e+11, 1.3, 3626.0e0, T, X, Gs, W );

  // reaction #7: H2 + O -> H + OH
  add_to_W_fwbw_2r2p ( specH2, specO, specH, specOH, 2.50e+15, 0.0, 6000.0e0, T, X, Gs, W );

  // reaction #8: H + O2 -> O + OH
  add_to_W_fwbw_2r2p ( specH, specO2, specO, specOH, 4.00e+14, 0.0, 18000.0e0, T, X, Gs, W );

  // reaction #9: H2 + H + H -> H2+ H2
  add_to_W_fwbw_3r2p ( specH2, specH, specH, specH2, specH2, 4.00e+20, -1.0, 0.0e0, T, X, Gs, W );

  // reaction #10: H + O2 -> HO2
  add_to_W_fwbw_2r1p ( specH, specO2, specHO2, 1.00e+15, -0.87, 0.0e0, T, X, Gs, W );

  // reaction #11: H + HO2 -> H2 + O2
  add_to_W_fwbw_2r2p ( specH, specHO2, specH2, specO2, 1.50e+14, 0.0e0, 0.0e0, T, X, Gs, W );

  // reaction #12: O + HO2 -> OH + O2
  add_to_W_fwbw_2r2p ( specO, specHO2, specOH, specO2, 2.50e+13, 0.0e0, 0.0e0, T, X, Gs, W );

  // reaction #13: CO + OH -> CO2 + H
  add_to_W_fwbw_2r2p ( specCO, specOH, specCO2, specH, 1.51e+07, 1.3e0, -758.0, T, X, Gs, W );

  // reaction #14: N2 + CH+ CH -> C2H2 + N2 
  add_to_W_fwbw_3r2p ( specN2, specCH, specCH, specC2H2, specN2, 1.00e+18, 0.0e0, -758.0, T, X, Gs, W );

  // reaction #15: C2H2 + O2 -> CO + CO + H2 
  add_to_W_fwbw_2r3p ( specC2H2, specO2, specCO, specCO, specH2, 3.00e+16, 0.0e0, 19000.0, T, X, Gs, W );

  // reaction #16: CH + O -> CO + H 
  add_to_W_fwbw_2r2p ( specCH, specO, specCO, specH, 1.00e+12, 0.7e0, 0.0, T, X, Gs, W );

  // reaction #17: CH + OH -> CO + H2 
  add_to_W_fwbw_2r2p ( specCH, specOH, specCO, specH2, 1.00e+13, 0.0e0, 0.0, T, X, Gs, W );

  // reaction #18: CH + NO -> NH + CO 
  add_to_W_fwbw_2r2p ( specCH, specNO, specNH, specCO, 1.00e+11, 0.0e0, 0.0, T, X, Gs, W );

  // reaction #19: N2 + O -> N + NO 
  add_to_W_fwbw_2r2p ( specN2, specO, specN, specNO, 9.00e+13, 0.0e0, 75000.0, T, X, Gs, W );

  // reaction #20: N + O2 -> NO + O 
  add_to_W_fwbw_2r2p ( specN, specO2, specNO, specO, 6.30e+09, 1.0e0, 6300.0, T, X, Gs, W );
  // reaction #21: NO + H -> N + OH 
  add_to_W_fwbw_2r2p ( specNO, specH, specN, specOH, 1.00e+12, 0.0e0, 0.0, T, X, Gs, W );

  // reaction #22: NH + O -> NO + H 
  add_to_W_fwbw_2r2p ( specNH, specO, specNO, specH, 2.5e+4, 2.64e0, 0.0, T, X, Gs, W );

  // reaction #23:  NH + NO -> N2 + OH
  add_to_W_fwbw_2r2p ( specNH, specNO, specN2, specOH, 2.00e+15, -0.8, 0.0, T, X, Gs, W );

  }
}

void find_dW_dx ( spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, p, r, s;              /* counters */
  long row, col;
  double rho, frr, dfrrdT, brr, dbrrdT;
  double XM, tmp1, tmp2;
  double dWdX[ns][ns];
  spec_t X, w, Gs, dGsdT;


  for ( r = 0; r < ns; r++ ) {
    for ( s = 0; s < ns; s++ ) {
      dWdrhok[r][s] = 0.0;
    }
    dWdT[r] = 0.0;
    dWdTe[r] = 0.0;
    dWdTv[r] = 0.0;
    dWdQbeam[r] = 0.0;
  }

  if (T>Tmin){

  /* Get values for temperature, species composition, and density */

  rho = 0.0;
  for ( k = 0; k < ns; k++ )
    rho += rhok[k];

  for ( s = 0; s < ns; s++ ) {
    w[s] = rhok[s] / rho;
    X[s] = max ( -1e6, w[s] * rho / _calM ( s ) * 1.0e-06 );    /* moles per cm3 */
  }

  X[specH2] = max ( Xmin, X[specH2] );
  X[specN2] = max ( Xmin, X[specN2] );
  X[specC3H8] = max ( Xmin, X[specC3H8] );
  X[specO] = max ( Xmin, X[specO] );

  XM = 0.0;
  for ( k = 0; k < ns; k++ )
    XM += X[k];

  for ( k = 0; k < ns; k++ ) {
    dWdT[k] = 0.0e0;
    Gs[k] = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * _calM ( k );    /* J/moles */
    dGsdT[k] = ( _cpk_from_T ( k, T ) - T * _dsk_dT_from_T ( k, T ) - _sk_from_T ( k, T ) ) * _calM ( k );
  }

  for ( row = 0; row < ns; row++ ) {
    for ( col = 0; col < ns; col++ ) {
      dWdX[row][col] = 0.0e0;
    }
  }

//      Rgasconstant=1.987192004e0 cal/mol/K; 
// forward reaction 1: N2+C3H8 -> 3CH + 5H + N2
  frr = 4.50e+10 * exp ( -30000.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT = frr * 30000.0e0 / 1.987192004e0 / sqr ( T );
  tmp1 = pow ( X[specN2], 0.8 );
  tmp2 = pow ( X[specC3H8], 0.8 );
  dWdT[specCH] += +3.0 * dfrrdT * tmp1 * tmp2;
  dWdT[specH] += +5.0 * dfrrdT * tmp1 * tmp2;
  dWdT[specC3H8] -= +1.0 * dfrrdT * tmp1 * tmp2;

  dWdX[specCH][specN2] += +0.8 * 3.0 * frr * tmp1 / X[specN2] * tmp2;
  dWdX[specH][specN2] += +0.8 * 5.0 * frr * tmp1 / X[specN2] * tmp2;
  dWdX[specC3H8][specN2] -= +0.8 * 1.0 * frr * tmp1 / X[specN2] * tmp2;

  dWdX[specCH][specC3H8] += +0.8 * 3.0 * frr * tmp1 * tmp2 / X[specC3H8];
  dWdX[specH][specC3H8] += +0.8 * 5.0 * frr * tmp1 * tmp2 / X[specC3H8];
  dWdX[specC3H8][specC3H8] -= +0.8 * 1.0 * frr * tmp1 * tmp2 / X[specC3H8];

// forward reaction 1: N2+C12H23 -> 12CH + 11H + N2
//      frr  = 5.50e+10 * exp(-30000.0e0/(1.987192004e0*T))*pow(X[specN2],0.8)*pow(X[specC12H23],0.8);
//      dfrrdT=0.0;
//      dWdT[specCH] += +12.0*frr; 
//      dWdT[specH] += +11.0*frr; 
//      dWdT[specC3H8]-= +1.0*frr; 

// forward reaction 2: CH + H2 + N2 -> 2NH + CH
  frr = 1.00e+16 * exp ( -78000.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT = frr * 78000.0e0 / 1.987192004e0 / sqr ( T );
  tmp1 = pow ( X[specH2], 0.1 );
  dWdT[specNH] += +2.0 * dfrrdT * sqr ( X[specCH] ) * tmp1 * X[specN2];
  dWdT[specH2] -= +1.0 * dfrrdT * sqr ( X[specCH] ) * tmp1 * X[specN2];
  dWdT[specN2] -= +1.0 * dfrrdT * sqr ( X[specCH] ) * tmp1 * X[specN2];

  dWdX[specNH][specCH] += +2.0 * 2.0 * frr * ( X[specCH] ) * tmp1 * X[specN2];
  dWdX[specH2][specCH] -= +2.0 * 1.0 * frr * ( X[specCH] ) * tmp1 * X[specN2];
  dWdX[specN2][specCH] -= +2.0 * 1.0 * frr * ( X[specCH] ) * tmp1 * X[specN2];

  dWdX[specNH][specH2] += +0.1 * 2.0 * frr * sqr ( X[specCH] ) * tmp1 / X[specH2] * X[specN2];
  dWdX[specH2][specH2] -= +0.1 * 1.0 * frr * sqr ( X[specCH] ) * tmp1 / X[specH2] * X[specN2];
  dWdX[specN2][specH2] -= +0.1 * 1.0 * frr * sqr ( X[specCH] ) * tmp1 / X[specH2] * X[specN2];

  dWdX[specNH][specN2] += +2.0 * frr * sqr ( X[specCH] ) * tmp1;
  dWdX[specH2][specN2] -= +1.0 * frr * sqr ( X[specCH] ) * tmp1;
  dWdX[specN2][specN2] -= +1.0 * frr * sqr ( X[specCH] ) * tmp1;

// backward reaction 2: CH + 2NH -> H2 + N2 + CH
  brr = 1.95e+15 * exp ( -53900.0e0 / ( 1.987192004e0 * T ) );
  dbrrdT = brr * 53900.0e0 / 1.987192004e0 / sqr ( T );

  dWdT[specH2] += +1.0 * dbrrdT * X[specCH] * sqr ( X[specNH] );
  dWdT[specN2] += +1.0 * dbrrdT * X[specCH] * sqr ( X[specNH] );
  dWdT[specNH] -= +2.0 * dbrrdT * X[specCH] * sqr ( X[specNH] );

  dWdX[specH2][specCH] += +1.0 * brr * sqr ( X[specNH] );
  dWdX[specN2][specCH] += +1.0 * brr * sqr ( X[specNH] );
  dWdX[specNH][specCH] -= +2.0 * brr * sqr ( X[specNH] );

  dWdX[specH2][specNH] += +2.0 * 1.0 * brr * X[specCH] * X[specNH];
  dWdX[specN2][specNH] += +2.0 * 1.0 * brr * X[specCH] * X[specNH];
  dWdX[specNH][specNH] -= +2.0 * 2.0 * brr * X[specCH] * X[specNH];

// forward reaction 3: O + N2 + HO2 -> 2NO + H + O
  frr = 1.95e+08 * exp ( -41900.0e0 / ( 1.987192004e0 * T ) ) * sqrt ( T );
  dfrrdT = frr * 41900.0e0 / 1.987192004e0 / sqr ( T ) + 0.5 * frr / T;
  tmp1 = pow ( X[specO], 0.1 );
  tmp2 = sqrt ( X[specN2] );

  dWdT[specNO] += 2.0 * dfrrdT * tmp1 * tmp2 * X[specHO2];
  dWdT[specH] += 1.0 * dfrrdT * tmp1 * tmp2 * X[specHO2];
  dWdT[specN2] -= 1.0 * dfrrdT * tmp1 * tmp2 * X[specHO2];
  dWdT[specHO2] -= 1.0 * dfrrdT * tmp1 * tmp2 * X[specHO2];

  dWdX[specNO][specO] += 0.1 * 2.0 * frr * tmp1 / X[specO] * tmp2 * X[specHO2];
  dWdX[specH][specO] += 0.1 * 1.0 * frr * tmp1 / X[specO] * tmp2 * X[specHO2];
  dWdX[specN2][specO] -= 0.1 * 1.0 * frr * tmp1 / X[specO] * tmp2 * X[specHO2];
  dWdX[specHO2][specO] -= 0.1 * 1.0 * frr * tmp1 / X[specO] * tmp2 * X[specHO2];

  dWdX[specNO][specN2] += 0.5 * 2.0 * frr * tmp1 / tmp2 * X[specHO2];
  dWdX[specH][specN2] += 0.5 * 1.0 * frr * tmp1 / tmp2 * X[specHO2];
  dWdX[specN2][specN2] -= 0.5 * 1.0 * frr * tmp1 / tmp2 * X[specHO2];
  dWdX[specHO2][specN2] -= 0.5 * 1.0 * frr * tmp1 / tmp2 * X[specHO2];

  dWdX[specNO][specHO2] += 2.0 * frr * tmp1 * tmp2;
  dWdX[specH][specHO2] += 1.0 * frr * tmp1 * tmp2;
  dWdX[specN2][specHO2] -= 1.0 * frr * tmp1 * tmp2;
  dWdX[specHO2][specHO2] -= 1.0 * frr * tmp1 * tmp2;

// backward reaction 3: 2NO + H -> N2 + HO2
  brr = 1.25e+10 * exp ( -6000.0e0 / ( 1.987192004e0 * T ) );
  dbrrdT = brr * 6000.0e0 / 1.987192004e0 / sqr ( T );

  dWdT[specN2] += 1.0 * dbrrdT * X[specNO] * X[specH];
  dWdT[specHO2] += 1.0 * dbrrdT * X[specNO] * X[specH];
  dWdT[specNO] -= 2.0 * dbrrdT * X[specNO] * X[specH];
  dWdT[specH] -= 1.0 * dbrrdT * X[specNO] * X[specH];

  dWdX[specN2][specNO] += 1.0 * brr * X[specH];
  dWdX[specHO2][specNO] += 1.0 * brr * X[specH];
  dWdX[specNO][specNO] -= 2.0 * brr * X[specH];
  dWdX[specH][specNO] -= 1.0 * brr * X[specH];

  dWdX[specN2][specH] += 1.0 * brr * X[specNO];
  dWdX[specHO2][specH] += 1.0 * brr * X[specNO];
  dWdX[specNO][specH] -= 2.0 * brr * X[specNO];
  dWdX[specH][specH] -= 1.0 * brr * X[specNO];

// forward reaction 4: N2 + CO + HO2 -> CO2 + OH + N2
  frr = 3.50e+13 * exp ( -22934.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT = frr * 22934.0e0 / 1.987192004e0 / sqr ( T );
  tmp1 = pow ( X[specN2], 0.1 );

  dWdT[specCO2] += 1.0 * dfrrdT * tmp1 * X[specCO] * X[specHO2];
  dWdT[specOH] += 1.0 * dfrrdT * tmp1 * X[specCO] * X[specHO2];
  dWdT[specCO] -= 1.0 * dfrrdT * tmp1 * X[specCO] * X[specHO2];
  dWdT[specHO2] -= 1.0 * dfrrdT * tmp1 * X[specCO] * X[specHO2];

  dWdX[specCO2][specN2] += 0.1 * 1.0 * frr * tmp1 / X[specN2] * X[specCO] * X[specHO2];
  dWdX[specOH][specN2] += 0.1 * 1.0 * frr * tmp1 / X[specN2] * X[specCO] * X[specHO2];
  dWdX[specCO][specN2] -= 0.1 * 1.0 * frr * tmp1 / X[specN2] * X[specCO] * X[specHO2];
  dWdX[specHO2][specN2] -= 0.1 * 1.0 * frr * tmp1 / X[specN2] * X[specCO] * X[specHO2];

  dWdX[specCO2][specCO] += 1.0 * frr * tmp1 * X[specHO2];
  dWdX[specOH][specCO] += 1.0 * frr * tmp1 * X[specHO2];
  dWdX[specCO][specCO] -= 1.0 * frr * tmp1 * X[specHO2];
  dWdX[specHO2][specCO] -= 1.0 * frr * tmp1 * X[specHO2];

  dWdX[specCO2][specHO2] += 1.0 * frr * tmp1 * X[specCO];
  dWdX[specOH][specHO2] += 1.0 * frr * tmp1 * X[specCO];
  dWdX[specCO][specHO2] -= 1.0 * frr * tmp1 * X[specCO];
  dWdX[specHO2][specHO2] -= 1.0 * frr * tmp1 * X[specCO];

// forward reaction 5: N2 + O2 -> 2O + N2
  frr = 1.00e+18 * exp ( -122239.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT = frr * 122239.0e0 / 1.987192004e0 / sqr ( T );

  dWdT[specO] += 2.0 * dfrrdT * X[specN2] * X[specO2];
  dWdT[specO2] -= 1.0 * dfrrdT * X[specN2] * X[specO2];

  dWdX[specO][specN2] += 2.0 * frr * X[specO2];
  dWdX[specO2][specN2] -= 1.0 * frr * X[specO2];

  dWdX[specO][specO2] += 2.0 * frr * X[specN2];
  dWdX[specO2][specO2] -= 1.0 * frr * X[specN2];

// backward reaction 5: H2 + 2O -> O2 + H2
  brr = 1.00e+18;
  dbrrdT = 0.0;

  dWdT[specO] -= 2.0 * dbrrdT * X[specH2] * sqr ( X[specO] );
  dWdT[specO2] += 1.0 * dbrrdT * X[specH2] * sqr ( X[specO] );

  dWdX[specO][specH2] -= 2.0 * brr * sqr ( X[specO] );
  dWdX[specO2][specH2] += 1.0 * brr * sqr ( X[specO] );

  dWdX[specO][specO] -= 2.0 * 2.0 * brr * X[specH2] * X[specO];
  dWdX[specO2][specO] += 2.0 * 1.0 * brr * X[specH2] * X[specO];

  for ( k = 0; k < ns; k++ ) {
    for ( p = 0; p < ns; p++ ) {
      dWdrhok[k][p] = dWdX[k][p] / _calM ( p ) * _calM ( k );
      dWdX[k][p] = 0.0;
    }
    dWdT[k] = dWdT[k] * _calM ( k ) * 1E6;
  }

  // reaction #6: H2 + OH -> H2O + H
  add_to_dW_fwbw_2r2p ( specH2, specOH,
                            specH2O, specH, 1.17e+11, 1.3, 3626.0e0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #7: H2 + O -> H + OH
  add_to_dW_fwbw_2r2p ( specH2, specO,
                            specH, specOH, 2.50e+15, 0.0, 6000.0e0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #8: H + O2 -> O + OH
  add_to_dW_fwbw_2r2p ( specH, specO2,
                            specO, specOH, 4.00e+14, 0.0, 18000.0e0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #9: H2 + H + H -> 2H2
  add_to_dW_fwbw_3r2p ( specH2, specH, specH,
                            specH2, specH2, 4.00e+20, -1.0, 0.0e0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #10: H + O2 -> HO2
  add_to_dW_fwbw_2r1p ( specH, specO2, specHO2, 1.00e+15, -0.87, 0.0e0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #11: H + HO2 -> H2 + O2
  add_to_dW_fwbw_2r2p ( specH, specHO2,
                            specH2, specO2, 1.50e+14, 0.0e0, 0.0e0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #12: O + HO2 -> OH + O2
  add_to_dW_fwbw_2r2p ( specO, specHO2,
                            specOH, specO2, 2.50e+13, 0.0e0, 0.0e0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #13: CO + OH -> CO2 + H
  add_to_dW_fwbw_2r2p ( specCO, specOH,
                            specCO2, specH, 1.51e+07, 1.3e0, -758.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #14: N2 + CH+ CH -> C2H2 + N2 
  add_to_dW_fwbw_3r2p ( specN2, specCH, specCH,
                            specC2H2, specN2, 1.00e+18, 0.0e0, -758.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #15: C2H2 + O2 -> CO + CO + H2 
  add_to_dW_fwbw_2r3p ( specC2H2, specO2,
                            specCO, specCO, specH2, 3.00e+16, 0.0e0, 19000.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #16: CH + O -> CO + H 
  add_to_dW_fwbw_2r2p ( specCH, specO, specCO, specH, 1.00e+12, 0.7e0, 0.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #17: CH + OH -> CO + H2 
  add_to_dW_fwbw_2r2p ( specCH, specOH,
                            specCO, specH2, 1.00e+13, 0.0e0, 0.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #18: CH + NO -> NH + CO 
  add_to_dW_fwbw_2r2p ( specCH, specNO,
                            specNH, specCO, 1.00e+11, 0.0e0, 0.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #19: N2 + O -> N + NO 
  add_to_dW_fwbw_2r2p ( specN2, specO,
                            specN, specNO, 9.00e+13, 0.0e0, 75000.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #20: N + O2 -> NO + O 
  add_to_dW_fwbw_2r2p ( specN, specO2,
                            specNO, specO, 6.30e+09, 1.0e0, 6300.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #21: NO + H -> N + OH 
  add_to_dW_fwbw_2r2p ( specNO, specH, specN, specOH, 1.00e+12, 0.0e0, 0.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #22: NH + O -> NO + H 
  add_to_dW_fwbw_2r2p ( specNH, specO, specNO, specH, 2.5e+4, 2.64e0, 0.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  // reaction #23:  NH + NO -> N2 + OH
  add_to_dW_fwbw_2r2p ( specNH, specNO,
                            specN2, specOH, 2.00e+15, -0.8, 0.0, T, X, Gs, dGsdT, dWdT, dWdrhok );

  }
}
