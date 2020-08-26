// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2001 Giovanni Fusina
Copyright 2002 Timothy Hui

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

#define maxreact 21

#define Kcmin 1.0e-90

typedef double react_t[maxreact];


void find_W_Jachimowski1988old ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k, r, s;
  double X[ns];
  double w[ns], rho, frr[maxreact], Kc[maxreact];
  double Gs[ns], dGs[maxreact];
  double third, third6, third7, third8, third9, third19;

  /* Get values for temperature, species composition, and density */

  rho = 0.0;
  for ( k = 0; k < ns; k++ )
    rho += rhok[k];

  for ( s = 0; s < ns; s++ ) {
    w[s] = rhok[s] / rho;
  }

  /* Calculate properties for each species */

  /* Initialize matrices to store values for Gibbs free energy [J/mol] */

  for ( k = 0; k < ns; k++ ) {
    Gs[k] = 0.0e0;
  }

  /* Initialize matrices to store values for forward reaction rate,
     difference in Gibbs free energy */

  for ( k = 0; k < ns; k++ ) {
    X[k] = w[k] * rho / _calM ( k ) * 1.0e-06;
    Gs[k] = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * _calM ( k );
  }

  for ( r = 0; r < maxreact; r++ ) {
    frr[r] = 0.0e0;
    dGs[r] = 0.0e0;
  }

  /* Calculate difference in Gibbs for each reaction */

  dGs[1] = 2.0e0 * Gs[specOH] - ( Gs[specH2] + Gs[specO2] );
  dGs[2] = ( Gs[specOH] + Gs[specO] ) - ( Gs[specH] + Gs[specO2] );
  dGs[3] = ( Gs[specOH] + Gs[specH] ) - ( Gs[specO] + Gs[specH2] );
  dGs[4] = ( Gs[specH2O] + Gs[specH] ) - ( Gs[specOH] + Gs[specH2] );
  dGs[5] = ( Gs[specH2O] + Gs[specO] ) - ( 2.0e0 * Gs[specOH] );
  dGs[6] = Gs[specH2O] - ( Gs[specH] + Gs[specOH] );
  dGs[7] = Gs[specH2] - 2.0e0 * Gs[specH];
  dGs[8] = Gs[specOH] - ( Gs[specH] + Gs[specO] );
  dGs[9] = Gs[specHO2] - ( Gs[specH] + Gs[specO2] );
  dGs[10] = ( Gs[specH2] + Gs[specO2] ) - ( Gs[specHO2] + Gs[specH] );
  dGs[11] = 2.0e0 * Gs[specOH] - ( Gs[specHO2] + Gs[specH] );
  dGs[12] = ( Gs[specH2O] + Gs[specO] ) - ( Gs[specHO2] + Gs[specH] );
  dGs[13] = ( Gs[specO2] + Gs[specOH] ) - ( Gs[specHO2] + Gs[specO] );
  dGs[14] = ( Gs[specH2O] + Gs[specO2] ) - ( Gs[specHO2] + Gs[specOH] );
  dGs[15] = ( Gs[specH2O2] + Gs[specO2] ) - 2.0e0 * Gs[specHO2];
  dGs[16] = ( Gs[specH2] + Gs[specHO2] ) - ( Gs[specH] + Gs[specH2O2] );
  dGs[17] = ( Gs[specOH] + Gs[specHO2] ) - ( Gs[specO] + Gs[specH2O2] );
  dGs[18] = ( Gs[specH2O] + Gs[specHO2] ) - ( Gs[specOH] + Gs[specH2O2] );
  dGs[19] = 2.0e0 * Gs[specOH] - Gs[specH2O2];
  dGs[20] = Gs[specO2] - 2.0e0 * Gs[specO];

  /* Determine Kp (storing results in Kc since this is ultimately
     what is sought here). [(atm) ** (n-m)] */

  for ( r = 0; r < maxreact; r++ ) {
    Kc[r] = 0.0e0;
  }

  for ( r = 1; r < maxreact; r++ ) {
    Kc[r] = max ( Kcmin, exp ( -dGs[r] / ( calR * T ) ) );
  }

  /* 
     determine Kc by multiplying Kp by 
     ((RuT/100000*100**3)**(prod-react))
     [Kc] = (cm**3/mol) ** (prod-react) 

     Note: WARP version prior Feb09 2006 has the form
     ((RuT/101325*100**3)**(prod-react))
     the change from 101325 to 100000 is due to the 
     implementation of NASA new thermo-coefficient data set
   */

  Kc[6] = Kc[6] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[7] = Kc[7] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[8] = Kc[8] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[9] = Kc[9] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[19] = Kc[19] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[20] = Kc[20] * ( 1E6 * calR / THERMO_P_REF * T );

  /* Calculate forward reaction rate, here a different R is used
     where R = 1.987192004 cal/g-K. [(cm**3/mol)**(m-1)/s] */

  frr[1] = 1.70e+13 * exp ( -48000.0e0 / ( 1.987192004e0 * T ) );
  frr[2] = 2.60e+14 * exp ( -16800.0e0 / ( 1.987192004e0 * T ) );
  frr[3] = 1.80e+10 * exp ( -8900.0e0 / ( 1.987192004e0 * T ) ) * T;
  frr[4] = 2.20e+13 * exp ( -5150.0e0 / ( 1.987192004e0 * T ) );
  frr[5] = 6.30e+12 * exp ( -1090.0e0 / ( 1.987192004e0 * T ) );
  frr[6] = 2.20e+22 * pow ( T, ( -2.0e0 ) );
  frr[7] = 6.40e+17 / T;
  frr[8] = 6.00e+16 * pow ( T, ( -0.6e0 ) );
  frr[9] = 2.10e+15 * exp ( 1000.0e0 / ( 1.987192004e0 * T ) );
  frr[10] = 1.30e+13;
  frr[11] = 1.40e+14 * exp ( -1080.0e0 / ( 1.987192004e0 * T ) );
  frr[12] = 1.00e+13 * exp ( -1080.0e0 / ( 1.987192004e0 * T ) );
  frr[13] = 1.50e+13 * exp ( -950.0e0 / ( 1.987192004e0 * T ) );
  frr[14] = 8.00e+12;
  frr[15] = 2.00e+12;
  frr[16] = 1.40e+12 * exp ( -3600.0e0 / ( 1.987192004e0 * T ) );
  frr[17] = 1.40e+13 * exp ( -6400.0e0 / ( 1.987192004e0 * T ) );
  frr[18] = 6.10e+12 * exp ( -1430.0e0 / ( 1.987192004e0 * T ) );
  frr[19] = 1.20e+17 * exp ( -45500.0e0 / ( 1.987192004e0 * T ) );
  frr[20] = 6.00e+17 * exp ( 1800.0e0 / ( 1.987192004e0 * T ) );

  /* Calculate third body efficiencies */

  third = 0.0e0;
  for ( k = 0; k < ns; k++ ) {
    third = third + X[k];
  }
  third6 = third + 5.0e0 * X[specH2O];
  third7 = third + X[specH2] + 5.0e0 * X[specH2O];
  third8 = third + 4.0e0 * X[specH2O];
  third9 = third + X[specH2] + 15.0e0 * X[specH2O];
  third19 = third + 14.0e0 * X[specH2O];

/* Calculate the species production rates [kg/s] */

  for ( s = 0; s < ns; s++ ) {
    W[s] = 0.0e0;
  }

  W[specH2] = ( -frr[1] * ( X[specH2] * X[specO2] - 1.0e0 / Kc[1] * ( X[specOH] * X[specOH] ) )
                - frr[3] * ( X[specO] * X[specH2] - 1.0e0 / Kc[3] * ( X[specOH] * X[specH] ) )
                - frr[4] * ( X[specOH] * X[specH2] - 1.0e0 / Kc[4] * ( X[specH2O] * X[specH] ) )
                + frr[7] * ( X[specH] * X[specH] - 1.0e0 / Kc[7] * ( X[specH2] ) ) * third7
                + frr[10] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[10] * ( X[specH2] * X[specO2] ) )
                + frr[16] * ( X[specH] * X[specH2O2] - 1.0e0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
     ) * 1.0e+06 * _calM ( specH2 );

  W[specO2] = ( -frr[1] * ( X[specH2] * X[specO2] - 1.0e0 / Kc[1] * ( X[specOH] * X[specOH] ) )
                - frr[2] * ( X[specH] * X[specO2] - 1.0e0 / Kc[2] * ( X[specOH] * X[specO] ) )
                - frr[9] * ( X[specH] * X[specO2] - 1.0e0 / Kc[9] * ( X[specHO2] ) ) * third9
                + frr[10] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[10] * ( X[specH2] * X[specO2] ) )
                + frr[13] * ( X[specHO2] * X[specO] - 1.0e0 / Kc[13] * ( X[specO2] * X[specOH] ) )
                + frr[14] * ( X[specHO2] * X[specOH] - 1.0e0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                + frr[15] * ( X[specHO2] * X[specHO2] - 1.0e0 / Kc[15] * ( X[specH2O2] * X[specO2] ) )
                + frr[20] * ( X[specO] * X[specO] -
                              1.0e0 / Kc[20] * ( X[specO2] ) ) * third ) * 1.0e+06 * _calM ( specO2 );

  W[specH] = ( -frr[2] * ( X[specH] * X[specO2] - 1.0e0 / Kc[2] * ( X[specOH] * X[specO] ) )
               + frr[3] * ( X[specO] * X[specH2] - 1.0e0 / Kc[3] * ( X[specOH] * X[specH] ) )
               + frr[4] * ( X[specOH] * X[specH2] - 1.0e0 / Kc[4] * ( X[specH2O] * X[specH] ) )
               - frr[6] * ( X[specH] * X[specOH] - 1.0e0 / Kc[6] * ( X[specH2O] ) ) * third6
               - 2.0e0 * frr[7] * ( X[specH] * X[specH] - 1.0e0 / Kc[7] * ( X[specH2] ) ) * third7
               - frr[8] * ( X[specH] * X[specO] - 1.0e0 / Kc[8] * ( X[specOH] ) ) * third8
               - frr[9] * ( X[specH] * X[specO2] - 1.0e0 / Kc[9] * ( X[specHO2] ) ) * third9
               - frr[10] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[10] * ( X[specH2] * X[specO2] ) )
               - frr[11] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[11] * ( X[specOH] * X[specOH] ) )
               - frr[12] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[12] * ( X[specH2O] * X[specO] ) )
               - frr[16] * ( X[specH] * X[specH2O2] - 1.0e0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
     ) * 1.0e+06 * _calM ( specH );

  W[specO] = ( +frr[2] * ( X[specH] * X[specO2] - 1.0e0 / Kc[2] * ( X[specOH] * X[specO] ) )
               - frr[3] * ( X[specO] * X[specH2] - 1.0e0 / Kc[3] * ( X[specOH] * X[specH] ) )
               + frr[5] * ( X[specOH] * X[specOH] - 1.0e0 / Kc[5] * ( X[specH2O] * X[specO] ) )
               - frr[8] * ( X[specH] * X[specO] - 1.0e0 / Kc[8] * ( X[specOH] ) ) * third8
               + frr[12] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[12] * ( X[specH2O] * X[specO] ) )
               - frr[13] * ( X[specHO2] * X[specO] - 1.0e0 / Kc[13] * ( X[specO2] * X[specOH] ) )
               - frr[17] * ( X[specO] * X[specH2O2] - 1.0e0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
               - 2.0e0 * frr[20] * ( X[specO] * X[specO] -
                                     1.0e0 / Kc[20] * ( X[specO2] ) ) * third ) * 1.0e+06 * _calM ( specO );

  W[specOH] = ( +2.0e0 * frr[1] * ( X[specH2] * X[specO2] - 1.0e0 / Kc[1] * ( X[specOH] * X[specOH] ) )
                + frr[2] * ( X[specH] * X[specO2] - 1.0e0 / Kc[2] * ( X[specOH] * X[specO] ) )
                + frr[3] * ( X[specO] * X[specH2] - 1.0e0 / Kc[3] * ( X[specOH] * X[specH] ) )
                - frr[4] * ( X[specOH] * X[specH2] - 1.0e0 / Kc[4] * ( X[specH2O] * X[specH] ) )
                - 2.0e0 * frr[5] * ( X[specOH] * X[specOH] - 1.0e0 / Kc[5] * ( X[specH2O] * X[specO] ) )
                - frr[6] * ( X[specH] * X[specOH] - 1.0e0 / Kc[6] * ( X[specH2O] ) ) * third6
                + frr[8] * ( X[specH] * X[specO] - 1.0e0 / Kc[8] * ( X[specOH] ) ) * third8
                + 2.0e0 * frr[11] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[11] * ( X[specOH] * X[specOH] ) )
                + frr[13] * ( X[specHO2] * X[specO] - 1.0e0 / Kc[13] * ( X[specO2] * X[specOH] ) )
                - frr[14] * ( X[specHO2] * X[specOH] - 1.0e0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                + frr[17] * ( X[specO] * X[specH2O2] - 1.0e0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
                - frr[18] * ( X[specOH] * X[specH2O2] - 1.0e0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
                + 2.0e0 * frr[19] * ( X[specH2O2] -
                                      1.0e0 / Kc[19] * ( X[specOH] * X[specOH] ) ) * third19 ) * 1.0e+06 *
    _calM ( specOH );

  W[specH2O] = ( +frr[4] * ( X[specOH] * X[specH2] - 1.0e0 / Kc[4] * ( X[specH2O] * X[specH] ) )
                 + frr[5] * ( X[specOH] * X[specOH] - 1.0e0 / Kc[5] * ( X[specH2O] * X[specO] ) )
                 + frr[6] * ( X[specH] * X[specOH] - 1.0e0 / Kc[6] * ( X[specH2O] ) ) * third6
                 + frr[12] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[12] * ( X[specH2O] * X[specO] ) )
                 + frr[14] * ( X[specHO2] * X[specOH] - 1.0e0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                 + frr[18] * ( X[specOH] * X[specH2O2] - 1.0e0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
     ) * 1.0e+06 * _calM ( specH2O );

  W[specHO2] = ( +frr[9] * ( X[specH] * X[specO2] - 1.0e0 / Kc[9] * ( X[specHO2] ) ) * third9
                 - frr[10] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[10] * ( X[specH2] * X[specO2] ) )
                 - frr[11] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[11] * ( X[specOH] * X[specOH] ) )
                 - frr[12] * ( X[specHO2] * X[specH] - 1.0e0 / Kc[12] * ( X[specH2O] * X[specO] ) )
                 - frr[13] * ( X[specHO2] * X[specO] - 1.0e0 / Kc[13] * ( X[specO2] * X[specOH] ) )
                 - frr[14] * ( X[specHO2] * X[specOH] - 1.0e0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                 - 2.0e0 * frr[15] * ( X[specHO2] * X[specHO2] -
                                       1.0e0 / Kc[15] * ( X[specH2O2] * X[specO2] ) )
                 + frr[16] * ( X[specH] * X[specH2O2] - 1.0e0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
                 + frr[17] * ( X[specO] * X[specH2O2] - 1.0e0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
                 + frr[18] * ( X[specOH] * X[specH2O2] - 1.0e0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
     ) * 1.0e+06 * _calM ( specHO2 );

  W[specH2O2] = ( +frr[15] * ( X[specHO2] * X[specHO2] - 1.0e0 / Kc[15] * ( X[specH2O2] * X[specO2] ) )
                  - frr[16] * ( X[specH] * X[specH2O2] - 1.0e0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
                  - frr[17] * ( X[specO] * X[specH2O2] - 1.0e0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
                  - frr[18] * ( X[specOH] * X[specH2O2] - 1.0e0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
                  - frr[19] * ( X[specH2O2] -
                                1.0e0 / Kc[19] * ( X[specOH] * X[specOH] ) ) * third19 ) * 1.0e+06 *
    _calM ( specH2O2 );
}

/***********************************************************************
 This function calculates the analytical Jacobian for the chemical
 source term.
************************************************************************/

void find_dW_dx_Jachimowski1988old ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, p, r, s;              /* counters */
  long row, col;
  double rho;
  react_t dGs;
  double third, third6, third7, third8, third9, third19;
  double temp6, temp7, temp8, temp9, temp19, temp19b, temp20;
  double dGsdT[ns];
  double dWsdX[ns][ns];
  react_t ddGsdT;
  react_t Kc, dKcdT, frr, dfrrdT;
  spec_t X, w, Gs;

  for ( r = 0; r < ns; r++ ) {
    for ( s = 0; s < ns; s++ ) {
      dWdrhok[r][s] = 0.0;
    }
    dWdT[r] = 0.0;
    dWdTe[r] = 0.0;
    dWdTv[r] = 0.0;
    dWdQbeam[r] = 0.0;
  }

  /* Get values for temperature, species composition, and density */

  rho = 0.0;
  for ( k = 0; k < ns; k++ )
    rho += rhok[k];

  for ( s = 0; s < ns; s++ ) {
    w[s] = rhok[s] / rho;
  }

  /* Calculate properties for each species */

  /* Initialize matrices to store values for Gibbs free energy [J/mol]
     and specific heat at constant volume [J/g-K] and dWi/dT */

  for ( k = 0; k < ns; k++ ) {
    Gs[k] = 0.0e0;
    dWdT[k] = 0.0e0;
  }

  /* Initialize matrices to store values for forward reaction rate
     difference in Gibbs free energy, rate of change of forward
     reaction rate wrt temperature, and rate of change of equilibrium
     constant wrt temperature. */

  for ( r = 0; r < maxreact; r++ ) {
    frr[r] = 0.0e0;
    dGs[r] = 0.0e0;
    dfrrdT[r] = 0.0e0;
    dKcdT[r] = 0.0e0;
  }

  /* Calculate concentrations [mol/cm**3], Gibbs free energy [J/mol],
     and specific heat at constant volume [J/g-K] */

  for ( k = 0; k < ns; k++ ) {
    X[k] = w[k] * rho / _calM ( k ) * 1.0e-06;
    Gs[k] = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * _calM ( k );
    dGsdT[k] = ( _cpk_from_T ( k, T ) - _sk_from_T ( k, T ) - T * _dsk_dT_from_T ( k, T )

       ) * _calM ( k );
  }

  /* Calculate difference in Gibbs for each reaction */

  dGs[1] = 2.0e0 * Gs[specOH] - ( Gs[specH2] + Gs[specO2] );
  dGs[2] = ( Gs[specOH] + Gs[specO] ) - ( Gs[specH] + Gs[specO2] );
  dGs[3] = ( Gs[specOH] + Gs[specH] ) - ( Gs[specO] + Gs[specH2] );
  dGs[4] = ( Gs[specH2O] + Gs[specH] ) - ( Gs[specOH] + Gs[specH2] );
  dGs[5] = ( Gs[specH2O] + Gs[specO] ) - 2.0e0 * Gs[specOH];
  dGs[6] = Gs[specH2O] - ( Gs[specH] + Gs[specOH] );
  dGs[7] = Gs[specH2] - 2.0e0 * Gs[specH];
  dGs[8] = Gs[specOH] - ( Gs[specH] + Gs[specO] );
  dGs[9] = Gs[specHO2] - ( Gs[specH] + Gs[specO2] );
  dGs[10] = ( Gs[specH2] + Gs[specO2] ) - ( Gs[specHO2] + Gs[specH] );
  dGs[11] = 2.0e0 * Gs[specOH] - ( Gs[specHO2] + Gs[specH] );
  dGs[12] = ( Gs[specH2O] + Gs[specO] ) - ( Gs[specHO2] + Gs[specH] );
  dGs[13] = ( Gs[specO2] + Gs[specOH] ) - ( Gs[specHO2] + Gs[specO] );
  dGs[14] = ( Gs[specH2O] + Gs[specO2] ) - ( Gs[specHO2] + Gs[specOH] );
  dGs[15] = ( Gs[specH2O2] + Gs[specO2] ) - 2.0e0 * Gs[specHO2];
  dGs[16] = ( Gs[specH2] + Gs[specHO2] ) - ( Gs[specH] + Gs[specH2O2] );
  dGs[17] = ( Gs[specOH] + Gs[specHO2] ) - ( Gs[specO] + Gs[specH2O2] );
  dGs[18] = ( Gs[specH2O] + Gs[specHO2] ) - ( Gs[specOH] + Gs[specH2O2] );
  dGs[19] = 2.0e0 * Gs[specOH] - Gs[specH2O2];
  dGs[20] = Gs[specO2] - 2.0e0 * Gs[specO];

  ddGsdT[1] = 2.0e0 * dGsdT[specOH] - ( dGsdT[specH2] + dGsdT[specO2] );
  ddGsdT[2] = ( dGsdT[specOH] + dGsdT[specO] ) - ( dGsdT[specH] + dGsdT[specO2] );
  ddGsdT[3] = ( dGsdT[specOH] + dGsdT[specH] ) - ( dGsdT[specO] + dGsdT[specH2] );
  ddGsdT[4] = ( dGsdT[specH2O] + dGsdT[specH] ) - ( dGsdT[specOH] + dGsdT[specH2] );
  ddGsdT[5] = ( dGsdT[specH2O] + dGsdT[specO] ) - 2.0e0 * dGsdT[specOH];
  ddGsdT[6] = dGsdT[specH2O] - ( dGsdT[specH] + dGsdT[specOH] );
  ddGsdT[7] = dGsdT[specH2] - 2.0e0 * dGsdT[specH];
  ddGsdT[8] = dGsdT[specOH] - ( dGsdT[specH] + dGsdT[specO] );
  ddGsdT[9] = dGsdT[specHO2] - ( dGsdT[specH] + dGsdT[specO2] );
  ddGsdT[10] = ( dGsdT[specH2] + dGsdT[specO2] ) - ( dGsdT[specHO2] + dGsdT[specH] );
  ddGsdT[11] = 2.0e0 * dGsdT[specOH] - ( dGsdT[specHO2] + dGsdT[specH] );
  ddGsdT[12] = ( dGsdT[specH2O] + dGsdT[specO] ) - ( dGsdT[specHO2] + dGsdT[specH] );
  ddGsdT[13] = ( dGsdT[specO2] + dGsdT[specOH] ) - ( dGsdT[specHO2] + dGsdT[specO] );
  ddGsdT[14] = ( dGsdT[specH2O] + dGsdT[specO2] ) - ( dGsdT[specHO2] + dGsdT[specOH] );
  ddGsdT[15] = ( dGsdT[specH2O2] + dGsdT[specO2] ) - 2.0e0 * dGsdT[specHO2];
  ddGsdT[16] = ( dGsdT[specH2] + dGsdT[specHO2] ) - ( dGsdT[specH] + dGsdT[specH2O2] );
  ddGsdT[17] = ( dGsdT[specOH] + dGsdT[specHO2] ) - ( dGsdT[specO] + dGsdT[specH2O2] );
  ddGsdT[18] = ( dGsdT[specH2O] + dGsdT[specHO2] ) - ( dGsdT[specOH] + dGsdT[specH2O2] );
  ddGsdT[19] = 2.0e0 * dGsdT[specOH] - dGsdT[specH2O2];
  ddGsdT[20] = dGsdT[specO2] - 2.0e0 * dGsdT[specO];

  /* Determine Kp (storing results in Kc since this is ultimately
     what is sought here). [(atm) ** (n-m)] */

  for ( r = 0; r < maxreact; r++ ) {
    Kc[r] = 0.0e0;
    dKcdT[r] = 0.0e0;
  }

  for ( r = 1; r < maxreact; r++ ) {
    Kc[r] = max ( Kcmin, exp ( -dGs[r] / ( calR * T ) ) );
    dKcdT[r] = Kc[r] * ( -ddGsdT[r] / ( calR * T ) + dGs[r] / ( calR * sqr ( T ) ) );
  }

  /* 
     determine Kc by multiplying Kp by 
     ((RuT/100000*100**3)**(prod-react))
     [Kc] = (cm**3/mol) ** (prod-react) 

     Note: WARP version prior Feb09 2006 has the form
     ((RuT/101325*100**3)**(prod-react))
     the change from 101325 to 100000 is due to the 
     implementation of NASA new thermo-coefficient data set
   */

  Kc[6] = Kc[6] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[7] = Kc[7] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[8] = Kc[8] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[9] = Kc[9] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[19] = Kc[19] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[20] = Kc[20] * ( 1E6 * calR / THERMO_P_REF * T );

  dKcdT[6] = dKcdT[6] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[6] / T;
  dKcdT[7] = dKcdT[7] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[7] / T;
  dKcdT[8] = dKcdT[8] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[8] / T;
  dKcdT[9] = dKcdT[9] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[9] / T;
  dKcdT[19] = dKcdT[19] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[19] / T;
  dKcdT[20] = dKcdT[20] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[20] / T;

  /* Calculate forward reaction rate, here a different R is used
     where R = 1.987192004 cal/g-K. [(cm**3/mol)**(m-1)/s] */

  frr[1] = 1.70e+13 * exp ( -48000.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[1] = frr[1] * ( 48000.0e0 / ( 1.987192004e0 * T * T ) );

  frr[2] = 2.60e+14 * exp ( -16800.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[2] = frr[2] * ( 16800.0e0 / ( 1.987192004e0 * T * T ) );

  frr[3] = 1.80e+10 * exp ( -8900.0e0 / ( 1.987192004e0 * T ) ) * T;
  dfrrdT[3] = frr[3] / T + frr[3] * ( 8900.0e0 / ( 1.987192004e0 * T * T ) );

  frr[4] = 2.20e+13 * exp ( -5150.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[4] = frr[4] * ( 5150.0e0 / ( 1.987192004e0 * T * T ) );

  frr[5] = 6.30e+12 * exp ( -1090.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[5] = frr[5] * ( 1090.0e0 / ( 1.987192004e0 * T * T ) );

  frr[6] = 2.20e+22 * pow ( T, ( -2.0e0 ) );
  dfrrdT[6] = -2.0e0 * 2.20e+22 / ( T * T * T );

  frr[7] = 6.40e+17 / T;
  dfrrdT[7] = -6.40e+17 / ( T * T );

  frr[8] = 6.00e+16 * pow ( T, ( -0.6e0 ) );
  dfrrdT[8] = -0.6e0 * 6.00e+16 * pow ( T, ( -1.6e0 ) );

  frr[9] = 2.10e+15 * exp ( 1000.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[9] = -frr[9] * ( 1000.0e0 / ( 1.987192004e0 * T * T ) );

  frr[10] = 1.30e+13;
  dfrrdT[10] = 0.0e0;

  frr[11] = 1.40e+14 * exp ( -1080.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[11] = frr[11] * ( 1080.0e0 / ( 1.987192004e0 * T * T ) );

  frr[12] = 1.00e+13 * exp ( -1080.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[12] = frr[12] * ( 1080.0e0 / ( 1.987192004e0 * T * T ) );

  frr[13] = 1.50e+13 * exp ( -950.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[13] = frr[13] * ( 950.0e0 / ( 1.987192004e0 * T * T ) );

  frr[14] = 8.00e+12;
  dfrrdT[14] = 0.0e0;

  frr[15] = 2.00e+12;
  dfrrdT[15] = 0.0e0;

  frr[16] = 1.40e+12 * exp ( -3600.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[16] = frr[16] * ( 3600.0e0 / ( 1.987192004e0 * T * T ) );

  frr[17] = 1.40e+13 * exp ( -6400.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[17] = frr[17] * ( 6400.0e0 / ( 1.987192004e0 * T * T ) );

  frr[18] = 6.10e+12 * exp ( -1430.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[18] = frr[18] * ( 1430.0e0 / ( 1.987192004e0 * T * T ) );

  frr[19] = 1.20e+17 * exp ( -45500.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[19] = frr[19] * ( 45500.0e0 / ( 1.987192004e0 * T * T ) );

  frr[20] = 6.00e+17 * exp ( 1800.0e0 / ( 1.987192004e0 * T ) );
  dfrrdT[20] = -frr[20] * ( 1800.0e0 / ( 1.987192004e0 * T * T ) );

  /* Calculate third body efficiencies */

  third = 0.0e0;
  for ( k = 0; k < ns; k++ ) {
    third = third + X[k];
  }
  third6 = third + 5.0e0 * X[specH2O];
  third7 = third + X[specH2] + 5.0e0 * X[specH2O];
  third8 = third + 4.0e0 * X[specH2O];
  third9 = third + X[specH2] + 15.0e0 * X[specH2O];
  third19 = third + 14.0e0 * X[specH2O];

/*
	W[specH2] = (-frr[1]  * (X[specH2]*X[specO2] - 1.0e0/Kc[1]*(X[specOH]*X[specOH]))
          	-frr[3]  * (X[specO]*X[specH2] - 1.0e0/Kc[3]*(X[specOH]*X[specH]))
           	-frr[4]  * (X[specOH]*X[specH2] - 1.0e0/Kc[4]*(X[specH2O]*X[specH]))
           	+frr[7]  * (X[specH]*X[specH] - 1.0e0/Kc[7]*(X[specH2]))*third7
           	+frr[10] * (X[specHO2]*X[specH] - 1.0e0/Kc[10]*(X[specH2]*X[specO2]))
           	+frr[16] * (X[specH]*X[specH2O2] - 1.0e0/Kc[16]*(X[specH2]*X[specHO2]))
           )*1.0e+06*_calM(specH2);
*/

  dWdT[specH2] = ( -dfrrdT[1] * ( X[specH2] * X[specO2] - 1.0 / Kc[1] * ( X[specOH] * X[specOH] ) )
                   - dfrrdT[3] * ( X[specO] * X[specH2] - 1.0 / Kc[3] * ( X[specOH] * X[specH] ) )
                   - dfrrdT[4] * ( X[specOH] * X[specH2] - 1.0 / Kc[4] * ( X[specH2O] * X[specH] ) )
                   + dfrrdT[7] * ( X[specH] * X[specH] - 1.0 / Kc[7] * ( X[specH2] ) ) * third7
                   + dfrrdT[10] * ( X[specHO2] * X[specH] - 1.0 / Kc[10] * ( X[specH2] * X[specO2] ) )
                   + dfrrdT[16] * ( X[specH] * X[specH2O2] - 1.0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
     );

  dWdT[specO2] = ( -dfrrdT[1] * ( X[specH2] * X[specO2] - 1.0 / Kc[1] * ( X[specOH] * X[specOH] ) )
                   - dfrrdT[2] * ( X[specH] * X[specO2] - 1.0 / Kc[2] * ( X[specOH] * X[specO] ) )
                   - dfrrdT[9] * ( X[specH] * X[specO2] - 1.0 / Kc[9] * ( X[specHO2] ) ) * third9
                   + dfrrdT[10] * ( X[specHO2] * X[specH] - 1.0 / Kc[10] * ( X[specH2] * X[specO2] ) )
                   + dfrrdT[13] * ( X[specHO2] * X[specO] - 1.0 / Kc[13] * ( X[specO2] * X[specOH] ) )
                   + dfrrdT[14] * ( X[specHO2] * X[specOH] - 1.0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                   + dfrrdT[15] * ( X[specHO2] * X[specHO2] - 1.0 / Kc[15] * ( X[specH2O2] * X[specO2] ) )
                   + dfrrdT[20] * ( X[specO] * X[specO] - 1.0 / Kc[20] * ( X[specO2] ) ) * third );

  dWdT[specH] = ( -dfrrdT[2] * ( X[specH] * X[specO2] - 1.0 / Kc[2] * ( X[specOH] * X[specO] ) )
                  + dfrrdT[3] * ( X[specO] * X[specH2] - 1.0 / Kc[3] * ( X[specOH] * X[specH] ) )
                  + dfrrdT[4] * ( X[specOH] * X[specH2] - 1.0 / Kc[4] * ( X[specH2O] * X[specH] ) )
                  - dfrrdT[6] * ( X[specH] * X[specOH] - 1.0 / Kc[6] * ( X[specH2O] ) ) * third6
                  - 2.0e0 * dfrrdT[7] * ( X[specH] * X[specH] - 1.0 / Kc[7] * ( X[specH2] ) ) * third7
                  - dfrrdT[8] * ( X[specH] * X[specO] - 1.0 / Kc[8] * ( X[specOH] ) ) * third8
                  - dfrrdT[9] * ( X[specH] * X[specO2] - 1.0 / Kc[9] * ( X[specHO2] ) ) * third9
                  - dfrrdT[10] * ( X[specHO2] * X[specH] - 1.0 / Kc[10] * ( X[specH2] * X[specO2] ) )
                  - dfrrdT[11] * ( X[specHO2] * X[specH] - 1.0 / Kc[11] * ( X[specOH] * X[specOH] ) )
                  - dfrrdT[12] * ( X[specHO2] * X[specH] - 1.0 / Kc[12] * ( X[specH2O] * X[specO] ) )
                  - dfrrdT[16] * ( X[specH] * X[specH2O2] - 1.0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
     );

  dWdT[specO] = ( +dfrrdT[2] * ( X[specH] * X[specO2] - 1.0 / Kc[2] * ( X[specOH] * X[specO] ) )
                  - dfrrdT[3] * ( X[specO] * X[specH2] - 1.0 / Kc[3] * ( X[specOH] * X[specH] ) )
                  + dfrrdT[5] * ( X[specOH] * X[specOH] - 1.0 / Kc[5] * ( X[specH2O] * X[specO] ) )
                  - dfrrdT[8] * ( X[specH] * X[specO] - 1.0 / Kc[8] * ( X[specOH] ) ) * third8
                  + dfrrdT[12] * ( X[specHO2] * X[specH] - 1.0 / Kc[12] * ( X[specH2O] * X[specO] ) )
                  - dfrrdT[13] * ( X[specHO2] * X[specO] - 1.0 / Kc[13] * ( X[specO2] * X[specOH] ) )
                  - dfrrdT[17] * ( X[specO] * X[specH2O2] - 1.0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
                  - 2.0e0 * dfrrdT[20] * ( X[specO] * X[specO] - 1.0 / Kc[20] * ( X[specO2] ) ) * third );

  dWdT[specOH] = ( +2.0e0 * dfrrdT[1] * ( X[specH2] * X[specO2] - 1.0 / Kc[1] * ( X[specOH] * X[specOH] ) )
                   + dfrrdT[2] * ( X[specH] * X[specO2] - 1.0 / Kc[2] * ( X[specOH] * X[specO] ) )
                   + dfrrdT[3] * ( X[specO] * X[specH2] - 1.0 / Kc[3] * ( X[specOH] * X[specH] ) )
                   - dfrrdT[4] * ( X[specOH] * X[specH2] - 1.0 / Kc[4] * ( X[specH2O] * X[specH] ) )
                   - 2.0e0 * dfrrdT[5] * ( X[specOH] * X[specOH] - 1.0 / Kc[5] * ( X[specH2O] * X[specO] ) )
                   - dfrrdT[6] * ( X[specH] * X[specOH] - 1.0 / Kc[6] * ( X[specH2O] ) ) * third6
                   + dfrrdT[8] * ( X[specH] * X[specO] - 1.0 / Kc[8] * ( X[specOH] ) ) * third8
                   + 2.0e0 * dfrrdT[11] * ( X[specHO2] * X[specH] - 1.0 / Kc[11] * ( X[specOH] * X[specOH] ) )
                   + dfrrdT[13] * ( X[specHO2] * X[specO] - 1.0 / Kc[13] * ( X[specO2] * X[specOH] ) )
                   - dfrrdT[14] * ( X[specHO2] * X[specOH] - 1.0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                   + dfrrdT[17] * ( X[specO] * X[specH2O2] - 1.0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
                   - dfrrdT[18] * ( X[specOH] * X[specH2O2] - 1.0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
                   + 2.0e0 * dfrrdT[19] * ( X[specH2O2] -
                                            1.0 / Kc[19] * ( X[specOH] * X[specOH] ) ) * third19 );

  dWdT[specH2O] = ( +dfrrdT[4] * ( X[specOH] * X[specH2] - 1.0 / Kc[4] * ( X[specH2O] * X[specH] ) )
                    + dfrrdT[5] * ( X[specOH] * X[specOH] - 1.0 / Kc[5] * ( X[specH2O] * X[specO] ) )
                    + dfrrdT[6] * ( X[specH] * X[specOH] - 1.0 / Kc[6] * ( X[specH2O] ) ) * third6
                    + dfrrdT[12] * ( X[specHO2] * X[specH] - 1.0 / Kc[12] * ( X[specH2O] * X[specO] ) )
                    + dfrrdT[14] * ( X[specHO2] * X[specOH] - 1.0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                    + dfrrdT[18] * ( X[specOH] * X[specH2O2] - 1.0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
     );

  dWdT[specHO2] = ( +dfrrdT[9] * ( X[specH] * X[specO2] - 1.0 / Kc[9] * ( X[specHO2] ) ) * third9
                    - dfrrdT[10] * ( X[specHO2] * X[specH] - 1.0 / Kc[10] * ( X[specH2] * X[specO2] ) )
                    - dfrrdT[11] * ( X[specHO2] * X[specH] - 1.0 / Kc[11] * ( X[specOH] * X[specOH] ) )
                    - dfrrdT[12] * ( X[specHO2] * X[specH] - 1.0 / Kc[12] * ( X[specH2O] * X[specO] ) )
                    - dfrrdT[13] * ( X[specHO2] * X[specO] - 1.0 / Kc[13] * ( X[specO2] * X[specOH] ) )
                    - dfrrdT[14] * ( X[specHO2] * X[specOH] - 1.0 / Kc[14] * ( X[specH2O] * X[specO2] ) )
                    - 2.0e0 * dfrrdT[15] * ( X[specHO2] * X[specHO2] -
                                             1.0 / Kc[15] * ( X[specH2O2] * X[specO2] ) )
                    + dfrrdT[16] * ( X[specH] * X[specH2O2] - 1.0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
                    + dfrrdT[17] * ( X[specO] * X[specH2O2] - 1.0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
                    + dfrrdT[18] * ( X[specOH] * X[specH2O2] - 1.0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
     );

  dWdT[specH2O2] = ( +dfrrdT[15] * ( X[specHO2] * X[specHO2] - 1.0 / Kc[15] * ( X[specH2O2] * X[specO2] ) )
                     - dfrrdT[16] * ( X[specH] * X[specH2O2] - 1.0 / Kc[16] * ( X[specH2] * X[specHO2] ) )
                     - dfrrdT[17] * ( X[specO] * X[specH2O2] - 1.0 / Kc[17] * ( X[specOH] * X[specHO2] ) )
                     - dfrrdT[18] * ( X[specOH] * X[specH2O2] - 1.0 / Kc[18] * ( X[specH2O] * X[specHO2] ) )
                     - dfrrdT[19] * ( X[specH2O2] - 1.0 / Kc[19] * ( X[specOH] * X[specOH] ) ) * third19 );

  dWdT[specN2] = 0.0e0;

/*--------------------------dWi/dKcj*dKcj/dT---------------------------*/
  /* Calculate dWi/dKcj*dKcj/dT */

/*
	W[specH2] = (-frr[1]  * (X[specH2]*X[specO2] - 1.0e0/Kc[1]*(X[specOH]*X[specOH]))
          	-frr[3]  * (X[specO]*X[specH2] - 1.0e0/Kc[3]*(X[specOH]*X[specH]))
           	-frr[4]  * (X[specOH]*X[specH2] - 1.0e0/Kc[4]*(X[specH2O]*X[specH]))
           	+frr[7]  * (X[specH]*X[specH] - 1.0e0/Kc[7]*(X[specH2]))*third7
           	+frr[10] * (X[specHO2]*X[specH] - 1.0e0/Kc[10]*(X[specH2]*X[specO2]))
           	+frr[16] * (X[specH]*X[specH2O2] - 1.0e0/Kc[16]*(X[specH2]*X[specHO2]))
           )*1.0e+06*_calM(specH2);
*/

  dWdT[specH2] += ( -frr[1] / sqr ( Kc[1] ) * ( X[specOH] * X[specOH] ) * dKcdT[1]
                    - frr[3] / sqr ( Kc[3] ) * ( X[specOH] * X[specH] ) * dKcdT[3]
                    - frr[4] / sqr ( Kc[4] ) * ( X[specH2O] * X[specH] ) * dKcdT[4]
                    + frr[7] / sqr ( Kc[7] ) * ( X[specH2] ) * third7 * dKcdT[7]
                    + frr[10] / sqr ( Kc[10] ) * ( X[specH2] * X[specO2] ) * dKcdT[10]
                    + frr[16] / sqr ( Kc[16] ) * ( X[specH2] * X[specHO2] ) * dKcdT[16]
     );

  dWdT[specO2] += ( -frr[1] / sqr ( Kc[1] ) * ( X[specOH] * X[specOH] ) * dKcdT[1]
                    - frr[2] / sqr ( Kc[2] ) * ( X[specOH] * X[specO] ) * dKcdT[2]
                    - frr[9] / sqr ( Kc[9] ) * ( X[specHO2] ) * third9 * dKcdT[9]
                    + frr[10] / sqr ( Kc[10] ) * ( X[specH2] * X[specO2] ) * dKcdT[10]
                    + frr[13] / sqr ( Kc[13] ) * ( X[specO2] * X[specOH] ) * dKcdT[13]
                    + frr[14] / sqr ( Kc[14] ) * ( X[specH2O] * X[specO2] ) * dKcdT[14]
                    + frr[15] / sqr ( Kc[15] ) * ( X[specH2O2] * X[specO2] ) * dKcdT[15]
                    + frr[20] / sqr ( Kc[20] ) * ( X[specO2] ) * third * dKcdT[20]
     );

  dWdT[specH] += ( -frr[2] / sqr ( Kc[2] ) * ( X[specOH] * X[specO] ) * dKcdT[2]
                   + frr[3] / sqr ( Kc[3] ) * ( X[specOH] * X[specH] ) * dKcdT[3]
                   + frr[4] / sqr ( Kc[4] ) * ( X[specH2O] * X[specH] ) * dKcdT[4]
                   - frr[6] / sqr ( Kc[6] ) * ( X[specH2O] ) * third6 * dKcdT[6]
                   - 2.0e0 * frr[7] / sqr ( Kc[7] ) * ( X[specH2] ) * third7 * dKcdT[7]
                   - frr[8] / sqr ( Kc[8] ) * ( X[specOH] ) * third8 * dKcdT[8]
                   - frr[9] / sqr ( Kc[9] ) * ( X[specHO2] ) * third9 * dKcdT[9]
                   - frr[10] / sqr ( Kc[10] ) * ( X[specH2] * X[specO2] ) * dKcdT[10]
                   - frr[11] / sqr ( Kc[11] ) * ( X[specOH] * X[specOH] ) * dKcdT[11]
                   - frr[12] / sqr ( Kc[12] ) * ( X[specH2O] * X[specO] ) * dKcdT[12]
                   - frr[16] / sqr ( Kc[16] ) * ( X[specH2] * X[specHO2] ) * dKcdT[16]
     );

  dWdT[specO] += ( +frr[2] / sqr ( Kc[2] ) * ( X[specOH] * X[specO] ) * dKcdT[2]
                   - frr[3] / sqr ( Kc[3] ) * ( X[specOH] * X[specH] ) * dKcdT[3]
                   + frr[5] / sqr ( Kc[5] ) * ( X[specH2O] * X[specO] ) * dKcdT[5]
                   - frr[8] / sqr ( Kc[8] ) * ( X[specOH] ) * third8 * dKcdT[8]
                   + frr[12] / sqr ( Kc[12] ) * ( X[specH2O] * X[specO] ) * dKcdT[12]
                   - frr[13] / sqr ( Kc[13] ) * ( X[specO2] * X[specOH] ) * dKcdT[13]
                   - frr[17] / sqr ( Kc[17] ) * ( X[specOH] * X[specHO2] ) * dKcdT[17]
                   - 2.0e0 * frr[20] / sqr ( Kc[20] ) * ( X[specO2] ) * third * dKcdT[20]
     );

  dWdT[specOH] += ( +2.0e0 * frr[1] / sqr ( Kc[1] ) * ( X[specOH] * X[specOH] ) * dKcdT[1]
                    + frr[2] / sqr ( Kc[2] ) * ( X[specOH] * X[specO] ) * dKcdT[2]
                    + frr[3] / sqr ( Kc[3] ) * ( X[specOH] * X[specH] ) * dKcdT[3]
                    - frr[4] / sqr ( Kc[4] ) * ( X[specH2O] * X[specH] ) * dKcdT[4]
                    - 2.0e0 * frr[5] / sqr ( Kc[5] ) * ( X[specH2O] * X[specO] ) * dKcdT[5]
                    - frr[6] / sqr ( Kc[6] ) * ( X[specH2O] ) * third6 * dKcdT[6]
                    + frr[8] / sqr ( Kc[8] ) * ( X[specOH] ) * third8 * dKcdT[8]
                    + 2.0e0 * frr[11] / sqr ( Kc[11] ) * ( X[specOH] * X[specOH] ) * dKcdT[11]
                    + frr[13] / sqr ( Kc[13] ) * ( X[specO2] * X[specOH] ) * dKcdT[13]
                    - frr[14] / sqr ( Kc[14] ) * ( X[specH2O] * X[specO2] ) * dKcdT[14]
                    + frr[17] / sqr ( Kc[17] ) * ( X[specOH] * X[specHO2] ) * dKcdT[17]
                    - frr[18] / sqr ( Kc[18] ) * ( X[specH2O] * X[specHO2] ) * dKcdT[18]
                    + 2.0e0 * frr[19] / sqr ( Kc[19] ) * ( X[specOH] * X[specOH] ) * third19 * dKcdT[19]
     );

  dWdT[specH2O] += ( +frr[4] / sqr ( Kc[4] ) * ( X[specH2O] * X[specH] ) * dKcdT[4]
                     + frr[5] / sqr ( Kc[5] ) * ( X[specH2O] * X[specO] ) * dKcdT[5]
                     + frr[6] / sqr ( Kc[6] ) * ( X[specH2O] ) * third6 * dKcdT[6]
                     + frr[12] / sqr ( Kc[12] ) * ( X[specH2O] * X[specO] ) * dKcdT[12]
                     + frr[14] / sqr ( Kc[14] ) * ( X[specH2O] * X[specO2] ) * dKcdT[14]
                     + frr[18] / sqr ( Kc[18] ) * ( X[specH2O] * X[specHO2] ) * dKcdT[18]
     );

  dWdT[specHO2] += ( +frr[9] / sqr ( Kc[9] ) * ( X[specHO2] ) * third9 * dKcdT[9]
                     - frr[10] / sqr ( Kc[10] ) * ( X[specH2] * X[specO2] ) * dKcdT[10]
                     - frr[11] / sqr ( Kc[11] ) * ( X[specOH] * X[specOH] ) * dKcdT[11]
                     - frr[12] / sqr ( Kc[12] ) * ( X[specH2O] * X[specO] ) * dKcdT[12]
                     - frr[13] / sqr ( Kc[13] ) * ( X[specO2] * X[specOH] ) * dKcdT[13]
                     - frr[14] / sqr ( Kc[14] ) * ( X[specH2O] * X[specO2] ) * dKcdT[14]
                     - 2.0e0 * frr[15] / sqr ( Kc[15] ) * ( X[specH2O2] * X[specO2] ) * dKcdT[15]
                     + frr[16] / sqr ( Kc[16] ) * ( X[specH2] * X[specHO2] ) * dKcdT[16]
                     + frr[17] / sqr ( Kc[17] ) * ( X[specOH] * X[specHO2] ) * dKcdT[17]
                     + frr[18] / sqr ( Kc[18] ) * ( X[specH2O] * X[specHO2] ) * dKcdT[18]
     );

  dWdT[specH2O2] += ( +frr[15] / sqr ( Kc[15] ) * ( X[specH2O2] * X[specO2] ) * dKcdT[15]
                      - frr[16] / sqr ( Kc[16] ) * ( X[specH2] * X[specHO2] ) * dKcdT[16]
                      - frr[17] / sqr ( Kc[17] ) * ( X[specOH] * X[specHO2] ) * dKcdT[17]
                      - frr[18] / sqr ( Kc[18] ) * ( X[specH2O] * X[specHO2] ) * dKcdT[18]
                      - frr[19] / sqr ( Kc[19] ) * ( X[specOH] * X[specOH] ) * third19 * dKcdT[19]
     );
  dWdT[specN2] = 0.0;

  for ( k = 0; k < ns; k++ )
    dWdT[k] = dWdT[k] * 1.0e6 * _calM ( k );

/*------------------------------------------------------------------------*/
  /* Now that dWi/dT is found, let's determine the derivatives required
     in order to calculate the Jacobian for the terms dWi/drhoj.  This
     term is composed of the following partial derivatives:

     dWi/drhoj = dWi/dXj*dXj/drhoj + dWi/dT*dT/drhoj

     Let's first determine dWi/dXj */

/*-----------------------------dWi/dXj------------------------------------*/
  /* Initialize the matrix the store values for dWsdX */

  for ( row = 0; row < ns; row++ ) {
    for ( col = 0; col < ns; col++ ) {
      dWsdX[row][col] = 0.0e0;
    }
  }

  temp6 = frr[6] * ( X[specH] * X[specOH] - 1.0e0 / Kc[6] * ( X[specH2O] ) );
  temp7 = frr[7] * ( X[specH] * X[specH] - 1.0e0 / Kc[7] * ( X[specH2] ) );
  temp8 = frr[8] * ( X[specH] * X[specO] - 1.0e0 / Kc[8] * ( X[specOH] ) );
  temp9 = frr[9] * ( X[specH] * X[specO2] - 1.0e0 / Kc[9] * ( X[specHO2] ) );
  temp20 = frr[20] * ( X[specO] * X[specO] - 1.0e0 / Kc[20] * ( X[specO2] ) );

  dWsdX[specH2][specH2] = ( -frr[1] * X[specO2] - frr[3] * X[specO] - frr[4] * X[specOH]
                            - frr[7] / Kc[7] * third7 + temp7 * 2.0e0
                            - frr[10] / Kc[10] * X[specO2] - frr[16] / Kc[16] * X[specHO2]
     );

  dWsdX[specH2][specO2] = ( -frr[1] * X[specH2] + temp7 - frr[10] / Kc[10] * X[specH2]
     );

  dWsdX[specH2][specH] = ( +frr[3] / Kc[3] * X[specOH] + frr[4] / Kc[4] * X[specH2O]
                           + temp7 + 2.0e0 * frr[7] * X[specH] * third7
                           + frr[10] * X[specHO2] + frr[16] * X[specH2O2]
     );

  dWsdX[specH2][specO] = ( -frr[3] * X[specH2] + temp7 );

  dWsdX[specH2][specOH] = ( +frr[1] / Kc[1] * 2.0e0 * X[specOH] + frr[3] / Kc[3] * X[specH]
                            - frr[4] * X[specH2] + temp7 );

  dWsdX[specH2][specH2O] = ( +frr[4] / Kc[4] * X[specH] + temp7 * 6.0e0 );

  dWsdX[specH2][specHO2] = ( temp7 + frr[10] * X[specH] - frr[16] / Kc[16] * X[specH2]
     );

  dWsdX[specH2][specH2O2] = ( temp7 + frr[16] * X[specH]
     );

  dWsdX[specH2][specN2] = ( temp7 );

  dWsdX[specO2][specH2] = ( -frr[1] * X[specO2] - temp9 * 2.0 - frr[10] / Kc[10] * X[specO2]
                            + temp20 );

  dWsdX[specO2][specO2] = ( -frr[1] * X[specH2] - frr[2] * X[specH]
                            - frr[9] * X[specH] * third9 - temp9
                            - frr[10] / Kc[10] * X[specH2] - frr[13] / Kc[13] * X[specOH]
                            - frr[14] / Kc[14] * X[specH2O] - frr[15] / Kc[15] * X[specH2O2]
                            - frr[20] / Kc[20] * third + temp20 );

  dWsdX[specO2][specH] = ( -frr[2] * X[specO2] - frr[9] * X[specO2] * third9 - temp9
                           + frr[10] * X[specHO2] + temp20 );

  dWsdX[specO2][specO] = ( +frr[2] / Kc[2] * X[specOH] - temp9 + frr[13] * X[specHO2]
                           + 2.0e0 * frr[20] * X[specO] * third + temp20 );

  dWsdX[specO2][specOH] = ( +frr[1] / Kc[1] * 2.0e0 * X[specOH] + frr[2] / Kc[2] * X[specO] - temp9
                            - frr[13] / Kc[13] * X[specO2] + frr[14] * X[specHO2] + temp20 );

  dWsdX[specO2][specH2O] = ( -temp9 * 16 - frr[14] / Kc[14] * X[specO2] + temp20 );

  dWsdX[specO2][specHO2] = ( +frr[9] / Kc[9] * third9 - temp9 + frr[10] * X[specH]
                             + frr[13] * X[specO] + frr[14] * X[specOH] + frr[15] * 2.0e0 * X[specHO2] +
                             temp20 );

  dWsdX[specO2][specH2O2] = ( -temp9 - frr[15] / Kc[15] * X[specO2] + temp20 );

  dWsdX[specO2][specN2] = ( -temp9 + temp20 );

  dWsdX[specH][specH2] = ( +frr[3] * X[specO] + frr[4] * X[specOH]
                           + 2.0e0 * frr[7] / Kc[7] * third7 + frr[10] / Kc[10] * X[specO2] +
                           frr[16] / Kc[16] * X[specHO2]
                           - temp6 - temp7 * 4.0 - temp8 - temp9 * 2.0 );

  dWsdX[specH][specO2] = ( -frr[2] * X[specH] - temp6 - temp7 * 2.0 - temp8 - frr[9] * X[specH]
                           * third9 - temp9 + frr[10] / Kc[10] * X[specH2]
     );

  dWsdX[specH][specH] = ( -frr[2] * X[specO2] - frr[3] / Kc[3] * X[specOH]
                          - frr[4] / Kc[4] * X[specH2O] - frr[6] * X[specOH] * third6 - temp6
                          - 2.0e0 * frr[7] * 2.0e0 * X[specH] * third7 - temp7 * 2.0 -
                          frr[8] * X[specO] * third8 - temp8 - frr[9] * X[specO2] * third9 - temp9 -
                          frr[10] * X[specHO2]
                          - frr[11] * X[specHO2] - frr[12] * X[specHO2] - frr[16] * X[specH2O2]
     );

  dWsdX[specH][specO] = ( +frr[2] / Kc[2] * X[specOH] + frr[3] * X[specH2] - temp6 - temp7 * 2.0
                          - frr[8] * X[specH] * third8 - temp8 - temp9 + frr[12] / Kc[12] * X[specH2O]
     );

  dWsdX[specH][specOH] = ( +frr[2] / Kc[2] * X[specO] - frr[3] / Kc[3] * X[specH] + frr[4] * X[specH2]
                           - frr[6] * X[specH] * third6 - temp6 - temp7 * 2.0
                           + frr[8] / Kc[8] * third8 - temp8 - temp9 + frr[11] / Kc[11] * 2.0e0 * X[specOH]
     );

  dWsdX[specH][specH2O] = ( -frr[4] / Kc[4] * X[specH] + frr[6] / Kc[6] * third6
                            - temp6 * 6.0e0 - temp7 * 12.0e0 - temp8 * 5.0e0 - temp9 * 16.0e0
                            + frr[12] / Kc[12] * X[specO]
     );

  dWsdX[specH][specHO2] = ( -temp6 - temp7 * 2.0 - temp8 + frr[9] / Kc[9] * third9
                            - temp9 - frr[10] * X[specH] - frr[11] * X[specH]
                            - frr[12] * X[specH] + frr[16] / Kc[16] * X[specH2]
     );

  dWsdX[specH][specH2O2] = ( -temp6 - temp7 * 2.0 - temp8 - temp9 - frr[16] * X[specH]
     );

  dWsdX[specH][specN2] = ( -temp6 - temp7 * 2.0 - temp8 - temp9 );

  temp20 = 2.0e0 * frr[20] * ( X[specO] * X[specO] - 1.0e0 / Kc[20] * ( X[specO2] ) );

  dWsdX[specO][specH2] = ( -frr[3] * X[specO] - temp8 - temp20 );

  dWsdX[specO][specO2] = ( +frr[2] * X[specH] - temp8 + frr[13] / Kc[13] * X[specOH]
                           + 2.0e0 * frr[20] / Kc[20] * third - temp20 );

  dWsdX[specO][specH] = ( +frr[2] * X[specO2] + frr[3] / Kc[3] * X[specOH]
                          - frr[8] * X[specO] * third8 - temp8 + frr[12] * X[specHO2] - temp20 );

  dWsdX[specO][specO] = ( -frr[2] / Kc[2] * X[specOH] - frr[3] * X[specH2]
                          - frr[5] / Kc[5] * X[specH2O] - frr[8] * X[specH] * third8 - temp8
                          - frr[12] / Kc[12] * X[specH2O] - frr[13] * X[specHO2]
                          - frr[17] * X[specH2O2] - 2.0e0 * frr[20] * 2.0e0 * X[specO] * third - temp20 );

  dWsdX[specO][specOH] = ( -frr[2] / Kc[2] * X[specO] + frr[3] / Kc[3] * X[specH]
                           + 2.0e0 * frr[5] * X[specOH] + frr[8] / Kc[8] * third8 - temp8
                           + frr[13] / Kc[13] * X[specO2] + frr[17] / Kc[17] * X[specHO2]
                           - temp20 );

  dWsdX[specO][specH2O] = ( -frr[5] / Kc[5] * X[specO] - temp8 * 5 - frr[12] / Kc[12] * X[specO] - temp20 );

  dWsdX[specO][specHO2] = ( -temp8 + frr[12] * X[specH] - frr[13] * X[specO]
                            + frr[17] / Kc[17] * X[specOH] - temp20 );

  dWsdX[specO][specH2O2] = ( -temp8 - frr[17] * X[specO] - temp20 );

  dWsdX[specO][specN2] = ( -temp8 - temp20 );

  temp19 = 2.0e0 * frr[19] * ( X[specH2O2] - 1.0 / Kc[19] * ( X[specOH] * X[specOH] ) );

  dWsdX[specOH][specH2] = ( +2.0e0 * frr[1] * X[specO2] + frr[3] * X[specO] - frr[4] * X[specOH]
                            - temp6 + temp8 + temp19 );

  dWsdX[specOH][specO2] = ( +2.0e0 * frr[1] * X[specH2] + frr[2] * X[specH] - temp6 + temp8
                            - frr[13] / Kc[13] * X[specOH] + frr[14] / Kc[14] * X[specH2O]
                            + temp19 );

  dWsdX[specOH][specH] = ( +frr[2] * X[specO2] - frr[3] / Kc[3] * X[specOH]
                           + frr[4] / Kc[4] * X[specH2O] - frr[6] * X[specOH] * third6 - temp6
                           + frr[8] * X[specO] * third8 + temp8 + 2.0e0 * frr[11] * X[specHO2] + temp19 );

  dWsdX[specOH][specO] = ( -frr[2] / Kc[2] * X[specOH] + frr[3] * X[specH2]
                           + 2.0e0 * frr[5] / Kc[5] * X[specH2O] - temp6 + frr[8] * X[specH] * third8 + temp8
                           + frr[13] * X[specHO2] + frr[17] * X[specH2O2] + temp19 );

  dWsdX[specOH][specOH] = ( -2.0e0 * frr[1] / Kc[1] * 2.0e0 * X[specOH] - frr[2] / Kc[2] * X[specO]
                            - frr[3] / Kc[3] * X[specH] - frr[4] * X[specH2] -
                            2.0e0 * frr[5] * 2.0e0 * X[specOH]
                            - frr[6] * X[specH] * third6 - temp6 - frr[8] / Kc[8] * third8 + temp8 -
                            2.0e0 * frr[11] / Kc[11] * 2.0e0 * X[specOH] - frr[13] / Kc[13] * X[specO2]
                            - frr[14] * X[specHO2] - frr[17] / Kc[17] * X[specHO2] - frr[18] * X[specH2O2]
                            - 2.0e0 * frr[19] / Kc[19] * 2.0e0 * X[specOH] * third19 + temp19 );

  dWsdX[specOH][specH2O] = ( +frr[4] / Kc[4] * X[specH] + 2.0e0 * frr[5] / Kc[5] * X[specO]
                             + frr[6] / Kc[6] * third6 - temp6 * 6.0e0 + temp8 * 5.0e0
                             + frr[14] / Kc[14] * X[specO2] + frr[18] / Kc[18] * X[specHO2]
                             + temp19 * 15.0e0 );

  dWsdX[specOH][specHO2] = ( -temp6 + temp8 + 2.0e0 * frr[11] * X[specH] + frr[13] * X[specO]
                             - frr[14] * X[specOH] - frr[17] / Kc[17] * X[specOH]
                             + frr[18] / Kc[18] * X[specH2O] + temp19 );

  dWsdX[specOH][specH2O2] = ( -temp6 + temp8 + frr[17] * X[specO] - frr[18] * X[specOH]
                              + 2.0e0 * frr[19] * third19 + temp19 );

  dWsdX[specOH][specN2] = ( -temp6 + temp8 + temp19 );

  dWsdX[specH2O][specH2] = ( +frr[4] * X[specOH] + temp6 );

  dWsdX[specH2O][specO2] = ( temp6 - frr[14] / Kc[14] * X[specH2O]
     );

  dWsdX[specH2O][specH] = ( -frr[4] / Kc[4] * X[specH2O] + frr[6] * X[specOH] * third6
                            + temp6 + frr[12] * X[specHO2]
     );

  dWsdX[specH2O][specO] = ( -frr[5] / Kc[5] * X[specH2O] + temp6 - frr[12] / Kc[12] * X[specH2O]
     );

  dWsdX[specH2O][specOH] = ( +frr[4] * X[specH2] + frr[5] * 2.0e0 * X[specOH]
                             + frr[6] * X[specH] * third6 + temp6 + frr[14] * X[specHO2]
                             + frr[18] * X[specH2O2]
     );

  dWsdX[specH2O][specH2O] = ( -frr[4] / Kc[4] * X[specH] - frr[5] / Kc[5] * X[specO]
                              - frr[6] / Kc[6] * third6 + temp6 * 6
                              - frr[12] / Kc[12] * X[specO] - frr[14] / Kc[14] * X[specO2]
                              - frr[18] / Kc[18] * X[specHO2]
     );

  dWsdX[specH2O][specHO2] = ( temp6 + frr[12] * X[specH] + frr[14] * X[specOH]
                              - frr[18] / Kc[18] * X[specH2O]
     );

  dWsdX[specH2O][specH2O2] = ( temp6 + frr[18] * X[specOH]
     );

  dWsdX[specH2O][specN2] = ( temp6 );

  dWsdX[specHO2][specH2] = ( temp9 * 2.0 + frr[10] / Kc[10] * X[specO2]
                             - frr[16] / Kc[16] * X[specHO2]
     );

  dWsdX[specHO2][specO2] = ( +frr[9] * X[specH] * third9 + temp9 + frr[10] / Kc[10] * X[specH2]
                             + frr[13] / Kc[13] * X[specOH] + frr[14] / Kc[14] * X[specH2O]
                             + frr[15] / Kc[15] * X[specH2O2] * 2.0e0 );

  dWsdX[specHO2][specH] = ( +frr[9] * X[specO2] * third9 + temp9 - frr[10] * X[specHO2]
                            - frr[11] * X[specHO2] - frr[12] * X[specHO2]
                            + frr[16] * X[specH2O2]
     );

  dWsdX[specHO2][specO] = ( temp9 + frr[12] / Kc[12] * X[specH2O] - frr[13] * X[specHO2]
                            + frr[17] * X[specH2O2]
     );

  dWsdX[specHO2][specOH] = ( temp9 + frr[11] / Kc[11] * 2.0e0 * X[specOH]
                             + frr[13] / Kc[13] * X[specO2] - frr[14] * X[specHO2]
                             - frr[17] / Kc[17] * X[specHO2] + frr[18] * X[specH2O2]
     );

  dWsdX[specHO2][specH2O] = ( temp9 * 16 + frr[12] / Kc[12] * X[specO] + frr[14] / Kc[14] * X[specO2]
                              - frr[18] / Kc[18] * X[specHO2]
     );

  dWsdX[specHO2][specHO2] = ( -frr[9] / Kc[9] * third9 + temp9 - frr[10] * X[specH]
                              - frr[11] * X[specH] - frr[12] * X[specH] - frr[13] * X[specO]
                              - frr[14] * X[specOH] - frr[15] * 4.0e0 * X[specHO2] -
                              frr[16] / Kc[16] * X[specH2]
                              - frr[17] / Kc[17] * X[specOH] - frr[18] / Kc[18] * X[specH2O]
     );

  dWsdX[specHO2][specH2O2] = ( temp9 + frr[15] / Kc[15] * X[specO2] * 2.0e0 + frr[16] * X[specH]
                               + frr[17] * X[specO] + frr[18] * X[specOH]
     );

  dWsdX[specHO2][specN2] = ( temp9 );

  temp19b = frr[19] * ( X[specH2O2] - 1.0 / Kc[19] * ( X[specOH] * X[specOH] ) );

  dWsdX[specH2O2][specH2] = ( +frr[16] / Kc[16] * X[specHO2] - temp19b );

  dWsdX[specH2O2][specO2] = ( -frr[15] / Kc[15] * X[specH2O2] - temp19b );

  dWsdX[specH2O2][specH] = ( -frr[16] * X[specH2O2] - temp19b );

  dWsdX[specH2O2][specO] = ( -frr[17] * X[specH2O2] - temp19b );

  dWsdX[specH2O2][specOH] = ( +frr[17] / Kc[17] * X[specHO2] - frr[18] * X[specH2O2]
                              + frr[19] / Kc[19] * 2.0e0 * X[specOH] * third19 - temp19b );

  dWsdX[specH2O2][specH2O] = ( +frr[18] / Kc[18] * X[specHO2] - temp19b * 15.0 );

  dWsdX[specH2O2][specHO2] = ( +frr[15] * 2.0e0 * X[specHO2] + frr[16] / Kc[16] * X[specH2]
                               + frr[17] / Kc[17] * X[specOH] + frr[18] / Kc[18] * X[specH2O]
                               - temp19b );

  dWsdX[specH2O2][specH2O2] = ( -frr[15] / Kc[15] * X[specO2] - frr[16] * X[specH] - frr[17] * X[specO]
                                - frr[18] * X[specOH] - frr[19] * third19 - temp19b );

  dWsdX[specH2O2][specN2] = ( -temp19b );

  for ( k = 0; k < ns; k++ ) {
    for ( p = 0; p < ns; p++ ) {
      dWdrhok[k][p] = dWsdX[k][p] / _calM ( p ) * _calM ( k );
    }
  }

}
