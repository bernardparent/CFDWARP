// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2005-2007 Jean P. Sislian

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
#include <model/.active/model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
#include <model/fluid/_fluid.h>

#define Tmin 400.0              //needed to prevent some reaction terms from being too small leading to blowups

#define nr 52
/* Note: ns is specified in chem.hh */

typedef double react_t[nr];


void write_model_chem_template(FILE **controlfile){
}


void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
}

/*const static long mR[nr][ns+1]=
             {*//* O    O2   H2   H2O  H    OH   HO2  CO   CO2  CH3  CH4  H2O2 CHO  CH2O CH3O C2H3 C2H4 C2H5 C2H6 N2 */
                                                                                                                                                                                                                        /* {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}, *//* 1 */
/* }; */

/* const static long mP[nr][ns+1]=
               {*//* O    O2   H2   H2O  H    OH   HO2  CO   CO2  CH3  CH4  H2O2 CHO  CH2O CH3O C2H3 C2H4 C2H5 C2H6 N2 */
                                                                                                                                                                                                                        /* {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}, *//* 1 */
/* }; */

  /* Remember: 
     by definition
     W[0]=qi*nk[1]/N  
     where 
     then
     Wp[13]=kf[13]*nk[1]=qi*nk[1]/N
     or 
     kf[13]=qi/N

   */

/*
This is the chemical file of Yungster and Rabinowitz Methane/Air Combustion Model
Reference: 
Yungster,S.,Rabinowitz,M.J.,"Computation of Shock-Induced Combustion Using a Detailed
Methane-Air Mechanism," Journal of Propulsion and Power, Volume 10, No. 5, pp609~617, 1994
*/

/*
The rate reactions are expressed by the empirical Arrhenius form
A (Cf in the code) = pre-exponential factor
E/Ru (E in the code) = activation energy/universal gas constant
b (n_f in the code) = temperature exponent

Reaction 33, 48, 49 are pressure-dependent
The extra coefficients are stored in arry number 53, 54, 55
*/

static double Cf[nr + 3] = {
  1.59E+17, 3.87E+04, 2.16E+08, 2.10E+08, 6.40E+17, 8.40E+21, 7.00E+17, 1.50E+14, 2.50E+13, 2.00E+13,
  6.02E+13, 1.00E+17, 1.22E+07, 3.01E+14, 7.23E+13, 3.00E+13, 1.00E+14, 4.20E+12, 1.86E+17, 1.26E+08,
  3.50E+13, 7.23E+05, 1.00E+14, 8.91E-13, 5.00E+16, 8.43E+13, 8.00E+12, 4.30E+13, 5.20E+13, 2.28E+13,
  3.20E+11, 4.90E+12, 7.05E+16, 7.80E+06, 1.90E+09, 5.60E+12, 1.50E+06, 4.60E+12, 2.00E+13, 5.00E+12,
  4.28E-13, 1.00E+14, 3.98E+12, 3.16E+11, 3.00E+13, 3.00E+13, 2.00E+12, 4.97E+10, 7.10E+25, 5.40E+02,
  2.20E+07, 5.50E-01, 1.19E+35, 6.24E+39, 2.23E+61
};

static double n_f[nr + 3] = {
  -0.927E+0, 2.700E+0, 1.510E+0, 1.400E+0, -1.000E+0, -2.000E+0, -0.800E+0, 0.000E+0, 0.000E+0, 0.000E+0,
  0.000E+0, 0.000E+0, 1.350E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, -1.000E+0, 1.620E+0,
  0.000E+0, 2.460E+0, 0.000E+0, 7.400E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0,
  0.500E+0, 0.000E+0, -0.558E+0, 2.110E+0, 1.440E+0, 0.000E+0, 2.130E+0, 0.000E+0, 0.000E+0, 0.000E+0,
  7.600E+0, 0.000E+0, 0.000E+0, 0.700E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.730E+0, -2.792E+0, 3.500E+0,
  1.900E+0, 4.000E+0, -4.911E+0, -6.800E+0, -11.992E+0
};

static double E[nr + 3] = {
  8491.28E+0, 3151.16E+0, 1725.92E+0, -199.65E+0, 0.0E+0, 0.0E+0, 0.0E+0, 505.15E+0, 348.79E+0, 0.0E+0,
  0.0E+0, 22851.89E+0, -317.52E+0, 1515.44E+0, 0.0E+0, 0.0E+0, 0.0E+0, 0.0E+0, 8551.42E+0, 1094.49E+0,
  1768.01E+0, -488.31E+0, 20085.61E+0, -481.09E+0, 38487.40E+0, 0.0E+0, 0.0E+0, 15503.20E+0,
  17559.87E+0, 0.0E+0,
  0.0E+0, 5905.41E+0, 52788.95E+0, 3896.85E+0, 4365.91E+0, 28179.99E+0, 1226.79E+0, 9056.57E+0,
  0.0E+0, 0.0E+0,
  -1775.23E+0, 12628.68E+0, -120.27E+0, 4029.15E+0, 1503.41E+0, 0.0E+0, 2513.71E+0, 18549.48E+0,
  46866.83E+0, 2621.95E+0,
  565.28E+0, 4173.48E+0, 53829.68E+0, 21384.56E+0, 49895.06E+0
};

static double a[3] = {
  0.555E+0, 0.667E+0, 0.805E+0
};

static double b[3] = {
  405.62E+0, 653.88E+0, 302.71E+0
};

static double w[3] = {
  4580.87E+0, 8733.74E+0, 10730.56E+0
};

/* species concentration [mol/cm^3] */
double _Xk ( long k, spec_t rhok ) {
  double tmp;
  tmp = rhok[k] / _calM ( k ) * 1.0e-6;
  return ( tmp );
}

/* gibbs energy of species [J/mol] */
double _Gk ( long k, double T ) {
  double tmp;
  tmp = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * _calM ( k );
  return ( tmp );
}

/* derivative of gibbs energy of species k wrt temperature [J/(mol K)] */
double _dGkdT ( long k, double T ) {
  double tmp;
  tmp = ( _cpk_from_T ( k, T ) - T * _dsk_dT_from_T ( k, T ) - _sk_from_T ( k, T ) ) * _calM ( k );
  return ( tmp );
}

/* 
forward reaction rate based on the modified Arrhenius equation 
[{(cm^3/mol)^x}/s] where x is a whole number multiple depending on the 
reaction under consideration 

reaction 33, 48, 49 are "pressure" (aka. concentration) dependent 
and require a different treatment

Note: Equation 15 of Yungster,S.;Rabinowitz, M.J.,"Computation of Shock-Induced
Combustion Using a Detailed Methane-Air Mechanism," is incorrect.
the correct formula is Fc^Xt instead of Fc*Xt
*/

double _kf ( long r, double third, double T ) {
  double tmp, kf_inf, kf_o, Fc, pr, xt;
  kf_o = 0.0E+0;
  Fc = 0.0E+0;
  if ( r == 32 || r == 47 || r == 48 ) {
    kf_inf = Cf[r] * pow ( T, n_f[r] ) * exp ( -E[r] / T );
    if ( r == 32 ) {
      kf_o = Cf[52] * pow ( T, n_f[52] ) * exp ( -E[52] / T );
      Fc = a[0] * exp ( -b[0] / T ) + ( 1.0E+0 - a[0] ) * exp ( -w[0] / T );
    }
    if ( r == 47 ) {
      kf_o = Cf[53] * pow ( T, n_f[53] ) * exp ( -E[53] / T );
      Fc = a[1] * exp ( -b[1] / T ) + ( 1.0E+0 - a[1] ) * exp ( -w[1] / T );
    }
    if ( r == 48 ) {
      kf_o = Cf[54] * pow ( T, n_f[54] ) * exp ( -E[54] / T );
      Fc = a[2] * exp ( -b[2] / T ) + ( 1.0E+0 - a[2] ) * exp ( -w[2] / T );
    }

    pr = third * kf_o / kf_inf;
    xt = 1.0E+0 / ( 1.0E+0 + log10 ( pr ) * log10 ( pr ) );
    tmp = kf_inf * pr / ( 1.0E+0 + pr ) * pow ( Fc, xt );
  } else {
    tmp = Cf[r] * pow ( T, n_f[r] ) * exp ( -E[r] / T );
  }
  return ( tmp );
}

double _dkfdT ( long r, double third, double T ) {
  double tmp, kf_o, Fc, pr, xt, xia, xib, xic, xiCf, xin_f, xiE, dkf_odT, dprdT, dxtdT, dFcdT;
  kf_o = 0.0E+0;
  Fc = 0.0E+0;
  if ( r == 32 || r == 47 || r == 48 ) {
    if ( r == 32 ) {
      xiCf = Cf[52];
      xin_f = n_f[52];
      xiE = E[52];
      xia = a[0];
      xib = b[0];
      xic = w[0];
    }
    if ( r == 47 ) {
      xiCf = Cf[53];
      xin_f = n_f[53];
      xiE = E[53];
      xia = a[1];
      xib = b[1];
      xic = w[1];
    }
    if ( r == 48 ) {
      xiCf = Cf[54];
      xin_f = n_f[54];
      xiE = E[54];
      xia = a[2];
      xib = b[2];
      xic = w[2];
    }

    kf_o = xiCf * pow ( T, xin_f ) * exp ( -xiE / T );
    dkf_odT = xiCf * exp ( -xiE / T ) * pow ( T, xin_f - 2.0 ) * ( xiE + xin_f * T );

    Fc = xia * exp ( -xib / T ) + ( 1.0E+0 - xia ) * exp ( -xic / T );
    dFcdT = +xia * xib * exp ( -xib / T ) / sqr ( T ) + xic * ( 1.0E+0 - xia ) * exp ( -xic / T ) / sqr ( T );

    pr = third * kf_o / ( Cf[r] * pow ( T, n_f[r] ) * exp ( -E[r] / T ) );
    dprdT = third * dkf_odT / ( Cf[r] * pow ( T, n_f[r] ) * exp ( -E[r] / T ) )
      + third * kf_o / Cf[r] * ( -exp ( E[r] / T ) * pow ( T, -n_f[r] - 2.0 ) * ( E[r] + n_f[r] * T ) );

    xt = 1.0E+0 / ( 1.0E+0 + log10 ( pr ) * log10 ( pr ) );
    dxtdT =
      -1.0E+0 / sqr ( 1.0E+0 + log10 ( pr ) * log10 ( pr ) ) * 2.0 * log10 ( pr ) / log ( 10.0 ) / pr * dprdT;

    //original non-derived: tmp=third*kf_o/(1.0E+0+pr)*pow(Fc,xt);

    tmp = third * dkf_odT / ( 1.0E+0 + pr ) * pow ( Fc, xt )
      - third * kf_o / sqr ( 1.0E+0 + pr ) * pow ( Fc, xt ) * dprdT
      + third * kf_o / ( 1.0E+0 + pr ) * pow ( Fc, xt ) * ( xt / Fc * dFcdT + log ( Fc ) * dxtdT );
  } else {
    tmp = Cf[r] * n_f[r] * pow ( T, n_f[r] - 1.0 ) * exp ( -E[r] / T )
      + Cf[r] * pow ( T, n_f[r] ) * exp ( -E[r] / T ) * E[r] / sqr ( T );
  }
  return ( tmp );
}

/* 
derivative of the forward reaction rate wrt the density of species k 

Notes: Since reaction 33, 48, and 49 are "pressure" (aka. concentration) dependent,
numerical interpretation is required
*/

double _dkfdT_old ( long r, double third, double T ) {
  double tmp;
  double Tp, Tm, kf_p, kf_m, kf;
  double kf_inf, kf_o, Fc, pr, xt;
  long cnt;

  Tp = T + 1.0;
  Tm = T - 1.0;
  kf_p = 0.0;
  kf_m = 0.0;
  kf_o = 0.0E+0;
  Fc = 0.0E+0;
  if ( r == 32 || r == 47 || r == 48 ) {
    for ( cnt = 1; cnt < 3; cnt++ ) {
      if ( cnt == 1 )
        T = Tm;
      if ( cnt == 2 )
        T = Tp;
      kf_inf = Cf[r] * pow ( T, n_f[r] ) * exp ( -E[r] / T );
      if ( r == 32 ) {
        kf_o = Cf[52] * pow ( T, n_f[52] ) * exp ( -E[52] / T );
        Fc = a[0] * exp ( -b[0] / T ) + ( 1.0E+0 - a[0] ) * exp ( -w[0] / T );
      }
      if ( r == 47 ) {
        kf_o = Cf[53] * pow ( T, n_f[53] ) * exp ( -E[53] / T );
        Fc = a[1] * exp ( -b[1] / T ) + ( 1.0E+0 - a[1] ) * exp ( -w[1] / T );
      }
      if ( r == 48 ) {
        kf_o = Cf[54] * pow ( T, n_f[54] ) * exp ( -E[54] / T );
        Fc = a[2] * exp ( -b[2] / T ) + ( 1.0E+0 - a[2] ) * exp ( -w[2] / T );
      }

      pr = third * kf_o / kf_inf;
      xt = 1.0E+0 / ( 1.0E+0 + log10 ( pr ) * log10 ( pr ) );
      kf = kf_inf * pr / ( 1.0E+0 + pr ) * pow ( Fc, xt );

      if ( cnt == 1 )
        kf_m = kf;
      if ( cnt == 2 )
        kf_p = kf;
    }
    tmp = ( kf_p - kf_m ) / ( Tp - Tm );
  } else {
    tmp = _kf ( r, third, T ) * ( n_f[r] / T + E[r] / ( T * T ) );
  }
  return ( tmp );
}

double _dkfdthird ( long r, double third, double T ) {
  double tmp, kf_inf, kf_o, Fc, pr, xt;
  kf_o = 0.0E+0;
  Fc = 0.0E+0;
  if ( r == 32 || r == 47 || r == 48 ) {
    kf_inf = Cf[r] * pow ( T, n_f[r] ) * exp ( -E[r] / T );
    if ( r == 32 ) {
      kf_o = Cf[52] * pow ( T, n_f[52] ) * exp ( -E[52] / T );
      Fc = a[0] * exp ( -b[0] / T ) + ( 1.0E+0 - a[0] ) * exp ( -w[0] / T );
    }
    if ( r == 47 ) {
      kf_o = Cf[53] * pow ( T, n_f[53] ) * exp ( -E[53] / T );
      Fc = a[1] * exp ( -b[1] / T ) + ( 1.0E+0 - a[1] ) * exp ( -w[1] / T );
    }
    if ( r == 48 ) {
      kf_o = Cf[54] * pow ( T, n_f[54] ) * exp ( -E[54] / T );
      Fc = a[2] * exp ( -b[2] / T ) + ( 1.0E+0 - a[2] ) * exp ( -w[2] / T );
    }

    pr = third * kf_o / kf_inf;
    xt = 1.0E+0 / ( 1.0E+0 + log10 ( pr ) * log10 ( pr ) );
    //tmp=kf_inf*pr/(1.0E0+pr)*pow(Fc,xt);
    tmp = kf_inf * kf_o / kf_inf / ( 1.0E0 + pr ) * pow ( Fc, xt )
      - kf_inf * kf_o / kf_inf * pr / sqr ( 1.0e0 + pr ) * pow ( Fc, xt )
      - kf_inf * pr / ( 1.0E0 + pr ) * pow ( Fc,
                                             xt ) * log ( Fc ) / sqr ( 1.0E+0 +
                                                                       log10 ( pr ) * log10 ( pr ) ) * 2.0 *
      log10 ( pr ) / pr / log ( 10.0 ) * kf_o / kf_inf;

  } else {
    tmp = 0.0;
  }
  return ( tmp );
}



void find_W ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {

  long spec;
  long r, k;
  double X[ns], Gs[ns];
  double dGs[nr], Kc[nr], kf[nr], react[nr];
  double third, third5, third11, third13, third18, third24, third41;

  T = max ( T, Tmin );

  for ( spec = 0; spec < ns; spec++ ) {
    W[spec] = 0.0E+00;
  }

  third = 0.0E+0;
  for ( k = 0; k < ns; k++ ) {
    X[k] = _Xk ( k, rhok );
    Gs[k] = _Gk ( k, T );
    third = third + X[k];
  }

/*
   O[0]    O2[1]    H2[2]   H2O[3]     H[4]    OH[5]   HO2[6]    CO[7]   CO2[8]  CH3[9]   
CH4[10] H2O2[11]  CHO[12] CH2O[13] CH3O[14] C2H3[15] C2H4[16] C2H5[17] C2H6[18]  N2[19]
*/

  third5 =
    third + 0.9E+0 * X[2] + 1.6E+0 * X[1] + 1.6E+0 * X[19] + 8.5E+0 * X[3] + 1.6E+0 * X[7] + 1.6E+0 * X[8];
  third11 =
    third + 1.9E+0 * X[2] + 0.2E+0 * X[1] + 0.2E+0 * X[19] + 17.5E+0 * X[3] + 1.1E+0 * X[7] + 3.3E+0 * X[8];
  third13 = third + 11.0E+0 * X[1] + 1.0E+0 * X[19] + 2.0E+0 * X[7] + 6.0E+0 * X[8];
  third18 = third + 0.87E+0 * X[2] + 7.12E+0 * X[3];
  third24 =
    third + 1.9E+0 * X[2] + 0.2E+0 * X[1] + 0.2E+0 * X[19] + 17.5E+0 * X[3] + 1.1E+0 * X[7] + 3.3E+0 * X[8];
  third41 =
    third + 1.9E+0 * X[2] + 0.2E+0 * X[1] + 0.2E+0 * X[19] + 17.5E+0 * X[3] + 1.1E+0 * X[7] + 3.3E+0 * X[8];

/* 
Calculate the difference between the gibbs free energy between the 
products and reactants of each reaction
*/

  dGs[0] = Gs[5] + Gs[0] - ( Gs[4] + Gs[1] );
  dGs[1] = Gs[5] + Gs[4] - ( Gs[0] + Gs[2] );
  dGs[2] = Gs[3] + Gs[4] - ( Gs[5] + Gs[2] );
  dGs[3] = Gs[0] + Gs[3] - ( Gs[5] + Gs[5] );
  dGs[4] = Gs[2] - ( Gs[4] + Gs[4] );
  dGs[5] = Gs[3] - ( Gs[4] + Gs[5] );
  dGs[6] = Gs[6] - ( Gs[4] + Gs[1] );
  dGs[7] = Gs[5] + Gs[5] - ( Gs[6] + Gs[4] );
  dGs[8] = Gs[2] + Gs[1] - ( Gs[6] + Gs[4] );
  dGs[9] = Gs[1] + Gs[5] - ( Gs[6] + Gs[0] );
  dGs[10] = Gs[3] + Gs[1] - ( Gs[6] + Gs[5] );
  dGs[11] = Gs[5] + Gs[5] - Gs[11];
  dGs[12] = Gs[8] + Gs[4] - ( Gs[7] + Gs[5] );
  dGs[13] = Gs[8] - ( Gs[7] + Gs[0] );
  dGs[14] = Gs[7] + Gs[2] - ( Gs[12] + Gs[4] );
  dGs[15] = Gs[7] + Gs[5] - ( Gs[12] + Gs[0] );
  dGs[16] = Gs[7] + Gs[3] - ( Gs[12] + Gs[5] );
  dGs[17] = Gs[7] + Gs[6] - ( Gs[12] + Gs[1] );
  dGs[18] = Gs[7] + Gs[4] - Gs[12];
  dGs[19] = Gs[12] + Gs[2] - ( Gs[13] + Gs[4] );
  dGs[20] = Gs[12] + Gs[5] - ( Gs[13] + Gs[0] );
  dGs[21] = Gs[12] + Gs[3] - ( Gs[13] + Gs[5] );
  dGs[22] = Gs[12] + Gs[6] - ( Gs[13] + Gs[1] );
  dGs[23] = Gs[12] + Gs[10] - ( Gs[13] + Gs[9] );
  dGs[24] = Gs[12] + Gs[4] - Gs[13];
  dGs[25] = Gs[13] + Gs[4] - ( Gs[9] + Gs[0] );
  dGs[26] = Gs[13] + Gs[2] - ( Gs[9] + Gs[5] );
  dGs[27] = Gs[14] + Gs[0] - ( Gs[9] + Gs[1] );
  dGs[28] = Gs[13] + Gs[5] - ( Gs[9] + Gs[1] );
  dGs[29] = Gs[14] + Gs[5] - ( Gs[9] + Gs[6] );
  dGs[30] = Gs[10] + Gs[7] - ( Gs[9] + Gs[12] );
  dGs[31] = Gs[17] + Gs[4] - ( Gs[9] + Gs[9] );
  dGs[32] = Gs[9] + Gs[4] - Gs[10];
  dGs[33] = Gs[9] + Gs[2] - ( Gs[10] + Gs[4] );
  dGs[34] = Gs[9] + Gs[5] - ( Gs[10] + Gs[0] );
  dGs[35] = Gs[9] + Gs[6] - ( Gs[10] + Gs[1] );
  dGs[36] = Gs[9] + Gs[3] - ( Gs[10] + Gs[5] );
  dGs[37] = Gs[9] + Gs[11] - ( Gs[10] + Gs[6] );
  dGs[38] = Gs[13] + Gs[2] - ( Gs[14] + Gs[4] );
  dGs[39] = Gs[13] + Gs[3] - ( Gs[14] + Gs[5] );
  dGs[40] = Gs[13] + Gs[6] - ( Gs[14] + Gs[1] );
  dGs[41] = Gs[13] + Gs[4] - Gs[14];
  dGs[42] = Gs[13] + Gs[12] - ( Gs[15] + Gs[1] );
  dGs[43] = Gs[15] + Gs[2] - ( Gs[16] + Gs[4] );
  dGs[44] = Gs[15] + Gs[3] - ( Gs[16] + Gs[5] );
  dGs[45] = Gs[16] + Gs[2] - ( Gs[17] + Gs[4] );
  dGs[46] = Gs[16] + Gs[6] - ( Gs[17] + Gs[1] );
  dGs[47] = Gs[16] + Gs[4] - Gs[17];
  dGs[48] = Gs[9] + Gs[9] - Gs[18];
  dGs[49] = Gs[17] + Gs[2] - ( Gs[18] + Gs[4] );
  dGs[50] = Gs[17] + Gs[3] - ( Gs[18] + Gs[5] );
  dGs[51] = Gs[17] + Gs[10] - ( Gs[18] + Gs[9] );

/* 
determine Kp (storing results in Kc since this is ultimately what is 
sought here) 
*/

  for ( r = 0; r < nr; r++ ) {
    Kc[r] = exp ( -dGs[r] / ( calR * T ) );
    kf[r] = _kf ( r, third, T );
  }

/* 
determine Kc 
Kc = Kp /((RuT/100000*100**3)**(prod-react))
[Kc] = (cm**3/mol) ** (prod-react) 

Note: WARP version prior Feb09 2006 has the form
Kc = Kp /((RuT/101325*100**3)**(prod-react))
the change from 101325 to 100000 is due to the 
implementation of NASA new thermo-coefficient data set
*/

  Kc[4] = Kc[4] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[5] = Kc[5] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[6] = Kc[6] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[11] = Kc[11] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[13] = Kc[13] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[18] = Kc[18] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[24] = Kc[24] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[32] = Kc[32] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[41] = Kc[41] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[47] = Kc[47] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[48] = Kc[48] / ( 1E6 * calR / THERMO_P_REF * T );

/* calculate the overall rate of individual reactions */

  react[0] = kf[0] * ( X[4] * X[1] - ( X[5] * X[0] ) / Kc[0] );
  react[1] = kf[1] * ( X[0] * X[2] - ( X[5] * X[4] ) / Kc[1] );
  react[2] = kf[2] * ( X[5] * X[2] - ( X[3] * X[4] ) / Kc[2] );
  react[3] = kf[3] * ( X[5] * X[5] - ( X[0] * X[3] ) / Kc[3] );
  react[4] = kf[4] * ( X[4] * X[4] - ( X[2] ) / Kc[4] ) * third;
  react[5] = kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * third5;
  react[6] = kf[6] * ( X[4] * X[1] - ( X[6] ) / Kc[6] ) * third;
  react[7] = kf[7] * ( X[6] * X[4] - ( X[5] * X[5] ) / Kc[7] );
  react[8] = kf[8] * ( X[6] * X[4] - ( X[2] * X[1] ) / Kc[8] );
  react[9] = kf[9] * ( X[6] * X[0] - ( X[1] * X[5] ) / Kc[9] );
  react[10] = kf[10] * ( X[6] * X[5] - ( X[3] * X[1] ) / Kc[10] );
  react[11] = kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * third11;
  react[12] = kf[12] * ( X[7] * X[5] - ( X[8] * X[4] ) / Kc[12] );
  react[13] = kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * third13;
  react[14] = kf[14] * ( X[12] * X[4] - ( X[7] * X[2] ) / Kc[14] );
  react[15] = kf[15] * ( X[12] * X[0] - ( X[7] * X[5] ) / Kc[15] );
  react[16] = kf[16] * ( X[12] * X[5] - ( X[7] * X[3] ) / Kc[16] );
  react[17] = kf[17] * ( X[12] * X[1] - ( X[7] * X[6] ) / Kc[17] );
  react[18] = kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * third18;
  react[19] = kf[19] * ( X[13] * X[4] - ( X[12] * X[2] ) / Kc[19] );
  react[20] = kf[20] * ( X[13] * X[0] - ( X[12] * X[5] ) / Kc[20] );
  react[21] = kf[21] * ( X[13] * X[5] - ( X[12] * X[3] ) / Kc[21] );
  react[22] = kf[22] * ( X[13] * X[1] - ( X[12] * X[6] ) / Kc[22] );
  react[23] = kf[23] * ( X[13] * X[9] - ( X[12] * X[10] ) / Kc[23] );
  react[24] = kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * third24;
  react[25] = kf[25] * ( X[9] * X[0] - ( X[13] * X[4] ) / Kc[25] );
  react[26] = kf[26] * ( X[9] * X[5] - ( X[13] * X[2] ) / Kc[26] );
  react[27] = kf[27] * ( X[9] * X[1] - ( X[14] * X[0] ) / Kc[27] );
  react[28] = kf[28] * ( X[9] * X[1] - ( X[13] * X[5] ) / Kc[28] );
  react[29] = kf[29] * ( X[9] * X[6] - ( X[14] * X[5] ) / Kc[29] );
  react[30] = kf[30] * ( X[9] * X[12] - ( X[10] * X[7] ) / Kc[30] );
  react[31] = kf[31] * ( X[9] * X[9] - ( X[17] * X[4] ) / Kc[31] );
  react[32] = kf[32] * ( X[10] - ( X[9] * X[4] ) / Kc[32] );
  react[33] = kf[33] * ( X[10] * X[4] - ( X[9] * X[2] ) / Kc[33] );
  react[34] = kf[34] * ( X[10] * X[0] - ( X[9] * X[5] ) / Kc[34] );
  react[35] = kf[35] * ( X[10] * X[1] - ( X[9] * X[6] ) / Kc[35] );
  react[36] = kf[36] * ( X[10] * X[5] - ( X[9] * X[3] ) / Kc[36] );
  react[37] = kf[37] * ( X[10] * X[6] - ( X[9] * X[11] ) / Kc[37] );
  react[38] = kf[38] * ( X[14] * X[4] - ( X[13] * X[2] ) / Kc[38] );
  react[39] = kf[39] * ( X[14] * X[5] - ( X[13] * X[3] ) / Kc[39] );
  react[40] = kf[40] * ( X[14] * X[1] - ( X[13] * X[6] ) / Kc[40] );
  react[41] = kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * third41;
  react[42] = kf[42] * ( X[15] * X[1] - ( X[13] * X[12] ) / Kc[42] );
  react[43] = kf[43] * ( X[16] * X[4] - ( X[15] * X[2] ) / Kc[43] );
  react[44] = kf[44] * ( X[16] * X[5] - ( X[15] * X[3] ) / Kc[44] );
  react[45] = kf[45] * ( X[17] * X[4] - ( X[16] * X[2] ) / Kc[45] );
  react[46] = kf[46] * ( X[17] * X[1] - ( X[16] * X[6] ) / Kc[46] );
  react[47] = kf[47] * ( X[17] - ( X[16] * X[4] ) / Kc[47] );
  react[48] = kf[48] * ( X[18] - ( X[9] * X[9] ) / Kc[48] );
  react[49] = kf[49] * ( X[18] * X[4] - ( X[17] * X[2] ) / Kc[49] );
  react[50] = kf[50] * ( X[18] * X[5] - ( X[17] * X[3] ) / Kc[50] );
  react[51] = kf[51] * ( X[18] * X[9] - ( X[17] * X[10] ) / Kc[51] );

/* calculate the chemical source term [kg_k/(m^3 s)] */
  W[0] =
    ( react[0] - react[1] + react[3] - react[9] - react[13] - react[15] - react[20] - react[25] + react[27]
      - react[34] )
    * 1.0E+6 * _calM ( 0 );
  W[1] =
    ( -react[0] - react[6] + react[8] + react[9] + react[10] - react[17] - react[22] - react[27] - react[28]
      - react[35] - react[40] - react[42] - react[46] )
    * 1.0E+6 * _calM ( 1 );
  W[2] =
    ( -react[1] - react[2] + react[4] + react[8] + react[14] + react[19] + react[26] + react[33] + react[38]
      + react[43] + react[45] + react[49] )
    * 1.0E+6 * _calM ( 2 );
  W[3] =
    ( react[2] + react[3] + react[5] + react[10] + react[16] + react[21] + react[36] + react[39] + react[44]
      + react[50] )
    * 1.0E+6 * _calM ( 3 );
  W[4] = ( -react[0] + react[1] + react[2] - react[4] - react[4] - react[5] - react[6] - react[7] - react[8]
           + react[12] - react[14] + react[18] - react[19] + react[24] + react[25] + react[31] + react[32] -
           react[33]
           - react[38] + react[41] - react[43] - react[45] + react[47] - react[49] )
    * 1.0E+6 * _calM ( 4 );
  W[5] = ( react[0] + react[1] - react[2] - react[3] - react[3] - react[5] + react[7] + react[7] + react[9]
           - react[10] + react[11] + react[11] - react[12] + react[15] - react[16] + react[20] - react[21] -
           react[26]
           + react[28] + react[29] + react[34] - react[36] - react[39] - react[44] - react[50] )
    * 1.0E+6 * _calM ( 5 );
  W[6] =
    ( react[6] - react[7] - react[8] - react[9] - react[10] + react[17] + react[22] - react[29] + react[35]
      - react[37] + react[40] + react[46] )
    * 1.0E+6 * _calM ( 6 );
  W[7] = ( -react[12] - react[13] + react[14] + react[15] + react[16] + react[17] + react[18] + react[30] )
    * 1.0E+6 * _calM ( 7 );
  W[8] = ( react[12] + react[13] )
    * 1.0E+6 * _calM ( 8 );
  W[9] =
    ( -react[23] - react[25] - react[26] - react[27] - react[28] - react[29] - react[30] - react[31] -
      react[31]
      + react[32] + react[33] + react[34] + react[35] + react[36] + react[37] + react[48] + react[48] -
      react[51] )
    * 1.0E+6 * _calM ( 9 );
  W[10] =
    ( react[23] + react[30] - react[32] - react[33] - react[34] - react[35] - react[36] - react[37] +
      react[51] )
    * 1.0E+6 * _calM ( 10 );
  W[11] = ( -react[11] + react[37] )
    * 1.0E+6 * _calM ( 11 );
  W[12] =
    ( -react[14] - react[15] - react[16] - react[17] - react[18] + react[19] + react[20] + react[21] +
      react[22]
      + react[23] + react[24] - react[30] + react[42] )
    * 1.0E+6 * _calM ( 12 );
  W[13] =
    ( -react[19] - react[20] - react[21] - react[22] - react[23] - react[24] + react[25] + react[26] +
      react[28]
      + react[38] + react[39] + react[40] + react[41] + react[42] )
    * 1.0E+6 * _calM ( 13 );
  W[14] = ( react[27] + react[29] - react[38] - react[39] - react[40] - react[41] )
    * 1.0E+6 * _calM ( 14 );
  W[15] = ( -react[42] + react[43] + react[44] )
    * 1.0E+6 * _calM ( 15 );
  W[16] = ( -react[43] - react[44] + react[45] + react[46] + react[47] )
    * 1.0E+6 * _calM ( 16 );
  W[17] = ( react[31] - react[45] - react[46] - react[47] + react[49] + react[50] + react[51] )
    * 1.0E+6 * _calM ( 17 );
  W[18] = ( -react[48] - react[49] - react[50] - react[51] )
    * 1.0E+6 * _calM ( 18 );

/*  for(k=0; k<ns; k++){
    W[k]=W[k]*_Omega(np,gl);
  }*/

}

void FinddKcdT ( double T, react_t Kc, react_t dKcdT ) {
  long k, r;
  react_t dGs;
  spec_t Gs, dGkdT;

  for ( k = 0; k < ns; k++ ) {
    Gs[k] = _Gk ( k, T );
    dGkdT[k] = _dGkdT ( k, T );
  }
  dGs[0] = Gs[5] + Gs[0] - ( Gs[4] + Gs[1] );
  dGs[1] = Gs[5] + Gs[4] - ( Gs[0] + Gs[2] );
  dGs[2] = Gs[3] + Gs[4] - ( Gs[5] + Gs[2] );
  dGs[3] = Gs[0] + Gs[3] - ( Gs[5] + Gs[5] );
  dGs[4] = Gs[2] - ( Gs[4] + Gs[4] );
  dGs[5] = Gs[3] - ( Gs[4] + Gs[5] );
  dGs[6] = Gs[6] - ( Gs[4] + Gs[1] );
  dGs[7] = Gs[5] + Gs[5] - ( Gs[6] + Gs[4] );
  dGs[8] = Gs[2] + Gs[1] - ( Gs[6] + Gs[4] );
  dGs[9] = Gs[1] + Gs[5] - ( Gs[6] + Gs[0] );
  dGs[10] = Gs[3] + Gs[1] - ( Gs[6] + Gs[5] );
  dGs[11] = Gs[5] + Gs[5] - Gs[11];
  dGs[12] = Gs[8] + Gs[4] - ( Gs[7] + Gs[5] );
  dGs[13] = Gs[8] - ( Gs[7] + Gs[0] );
  dGs[14] = Gs[7] + Gs[2] - ( Gs[12] + Gs[4] );
  dGs[15] = Gs[7] + Gs[5] - ( Gs[12] + Gs[0] );
  dGs[16] = Gs[7] + Gs[3] - ( Gs[12] + Gs[5] );
  dGs[17] = Gs[7] + Gs[6] - ( Gs[12] + Gs[1] );
  dGs[18] = Gs[7] + Gs[4] - Gs[12];
  dGs[19] = Gs[12] + Gs[2] - ( Gs[13] + Gs[4] );
  dGs[20] = Gs[12] + Gs[5] - ( Gs[13] + Gs[0] );
  dGs[21] = Gs[12] + Gs[3] - ( Gs[13] + Gs[5] );
  dGs[22] = Gs[12] + Gs[6] - ( Gs[13] + Gs[1] );
  dGs[23] = Gs[12] + Gs[10] - ( Gs[13] + Gs[9] );
  dGs[24] = Gs[12] + Gs[4] - Gs[13];
  dGs[25] = Gs[13] + Gs[4] - ( Gs[9] + Gs[0] );
  dGs[26] = Gs[13] + Gs[2] - ( Gs[9] + Gs[5] );
  dGs[27] = Gs[14] + Gs[0] - ( Gs[9] + Gs[1] );
  dGs[28] = Gs[13] + Gs[5] - ( Gs[9] + Gs[1] );
  dGs[29] = Gs[14] + Gs[5] - ( Gs[9] + Gs[6] );
  dGs[30] = Gs[10] + Gs[7] - ( Gs[9] + Gs[12] );
  dGs[31] = Gs[17] + Gs[4] - ( Gs[9] + Gs[9] );
  dGs[32] = Gs[9] + Gs[4] - Gs[10];
  dGs[33] = Gs[9] + Gs[2] - ( Gs[10] + Gs[4] );
  dGs[34] = Gs[9] + Gs[5] - ( Gs[10] + Gs[0] );
  dGs[35] = Gs[9] + Gs[6] - ( Gs[10] + Gs[1] );
  dGs[36] = Gs[9] + Gs[3] - ( Gs[10] + Gs[5] );
  dGs[37] = Gs[9] + Gs[11] - ( Gs[10] + Gs[6] );
  dGs[38] = Gs[13] + Gs[2] - ( Gs[14] + Gs[4] );
  dGs[39] = Gs[13] + Gs[3] - ( Gs[14] + Gs[5] );
  dGs[40] = Gs[13] + Gs[6] - ( Gs[14] + Gs[1] );
  dGs[41] = Gs[13] + Gs[4] - Gs[14];
  dGs[42] = Gs[13] + Gs[12] - ( Gs[15] + Gs[1] );
  dGs[43] = Gs[15] + Gs[2] - ( Gs[16] + Gs[4] );
  dGs[44] = Gs[15] + Gs[3] - ( Gs[16] + Gs[5] );
  dGs[45] = Gs[16] + Gs[2] - ( Gs[17] + Gs[4] );
  dGs[46] = Gs[16] + Gs[6] - ( Gs[17] + Gs[1] );
  dGs[47] = Gs[16] + Gs[4] - Gs[17];
  dGs[48] = Gs[9] + Gs[9] - Gs[18];
  dGs[49] = Gs[17] + Gs[2] - ( Gs[18] + Gs[4] );
  dGs[50] = Gs[17] + Gs[3] - ( Gs[18] + Gs[5] );
  dGs[51] = Gs[17] + Gs[10] - ( Gs[18] + Gs[9] );

  for ( r = 0; r < nr; r++ ) {
    Kc[r] = exp ( -dGs[r] / ( calR * T ) );
  }

/* find the derivative of each source term wrt temperature */
  for ( r = 0; r < nr; r++ ) {
    dKcdT[r] = ( dGs[r] / ( calR * T * T ) );
  }

  dKcdT[0] += -( dGkdT[5] + dGkdT[0] - ( dGkdT[4] + dGkdT[1] ) ) / calR / T;
  dKcdT[1] += -( dGkdT[5] + dGkdT[4] - ( dGkdT[0] + dGkdT[2] ) ) / calR / T;
  dKcdT[2] += -( dGkdT[3] + dGkdT[4] - ( dGkdT[5] + dGkdT[2] ) ) / calR / T;
  dKcdT[3] += -( dGkdT[0] + dGkdT[3] - ( dGkdT[5] + dGkdT[5] ) ) / calR / T;
  dKcdT[4] += -( dGkdT[2] - ( dGkdT[4] + dGkdT[4] ) ) / calR / T;
  dKcdT[5] += -( dGkdT[3] - ( dGkdT[4] + dGkdT[5] ) ) / calR / T;
  dKcdT[6] += -( dGkdT[6] - ( dGkdT[4] + dGkdT[1] ) ) / calR / T;
  dKcdT[7] += -( dGkdT[5] + dGkdT[5] - ( dGkdT[6] + dGkdT[4] ) ) / calR / T;
  dKcdT[8] += -( dGkdT[2] + dGkdT[1] - ( dGkdT[6] + dGkdT[4] ) ) / calR / T;
  dKcdT[9] += -( dGkdT[1] + dGkdT[5] - ( dGkdT[6] + dGkdT[0] ) ) / calR / T;
  dKcdT[10] += -( dGkdT[3] + dGkdT[1] - ( dGkdT[6] + dGkdT[5] ) ) / calR / T;
  dKcdT[11] += -( dGkdT[5] + dGkdT[5] - ( dGkdT[11] ) ) / calR / T;
  dKcdT[12] += -( dGkdT[8] + dGkdT[4] - ( dGkdT[7] + dGkdT[5] ) ) / calR / T;
  dKcdT[13] += -( dGkdT[8] - ( dGkdT[7] + dGkdT[0] ) ) / calR / T;
  dKcdT[14] += -( dGkdT[7] + dGkdT[2] - ( dGkdT[12] + dGkdT[4] ) ) / calR / T;
  dKcdT[15] += -( dGkdT[7] + dGkdT[5] - ( dGkdT[12] + dGkdT[0] ) ) / calR / T;
  dKcdT[16] += -( dGkdT[7] + dGkdT[3] - ( dGkdT[12] + dGkdT[5] ) ) / calR / T;
  dKcdT[17] += -( dGkdT[7] + dGkdT[6] - ( dGkdT[12] + dGkdT[1] ) ) / calR / T;
  dKcdT[18] += -( dGkdT[7] + dGkdT[4] - ( dGkdT[12] ) ) / calR / T;
  dKcdT[19] += -( dGkdT[12] + dGkdT[2] - ( dGkdT[13] + dGkdT[4] ) ) / calR / T;
  dKcdT[20] += -( dGkdT[12] + dGkdT[5] - ( dGkdT[13] + dGkdT[0] ) ) / calR / T;
  dKcdT[21] += -( dGkdT[12] + dGkdT[3] - ( dGkdT[13] + dGkdT[5] ) ) / calR / T;
  dKcdT[22] += -( dGkdT[12] + dGkdT[6] - ( dGkdT[13] + dGkdT[1] ) ) / calR / T;
  dKcdT[23] += -( dGkdT[12] + dGkdT[10] - ( dGkdT[13] + dGkdT[9] ) ) / calR / T;
  dKcdT[24] += -( dGkdT[12] + dGkdT[4] - ( dGkdT[13] ) ) / calR / T;
  dKcdT[25] += -( dGkdT[13] + dGkdT[4] - ( dGkdT[9] + dGkdT[0] ) ) / calR / T;
  dKcdT[26] += -( dGkdT[13] + dGkdT[2] - ( dGkdT[9] + dGkdT[5] ) ) / calR / T;
  dKcdT[27] += -( dGkdT[14] + dGkdT[0] - ( dGkdT[9] + dGkdT[1] ) ) / calR / T;
  dKcdT[28] += -( dGkdT[13] + dGkdT[5] - ( dGkdT[9] + dGkdT[1] ) ) / calR / T;
  dKcdT[29] += -( dGkdT[14] + dGkdT[5] - ( dGkdT[9] + dGkdT[6] ) ) / calR / T;
  dKcdT[30] += -( dGkdT[10] + dGkdT[7] - ( dGkdT[9] + dGkdT[12] ) ) / calR / T;
  dKcdT[31] += -( dGkdT[17] + dGkdT[4] - ( dGkdT[9] + dGkdT[9] ) ) / calR / T;
  dKcdT[32] += -( dGkdT[9] + dGkdT[4] - ( dGkdT[10] ) ) / calR / T;
  dKcdT[33] += -( dGkdT[9] + dGkdT[2] - ( dGkdT[10] + dGkdT[4] ) ) / calR / T;
  dKcdT[34] += -( dGkdT[9] + dGkdT[5] - ( dGkdT[10] + dGkdT[0] ) ) / calR / T;
  dKcdT[35] += -( dGkdT[9] + dGkdT[6] - ( dGkdT[10] + dGkdT[1] ) ) / calR / T;
  dKcdT[36] += -( dGkdT[9] + dGkdT[3] - ( dGkdT[10] + dGkdT[5] ) ) / calR / T;
  dKcdT[37] += -( dGkdT[9] + dGkdT[11] - ( dGkdT[10] + dGkdT[6] ) ) / calR / T;
  dKcdT[38] += -( dGkdT[13] + dGkdT[2] - ( dGkdT[14] + dGkdT[4] ) ) / calR / T;
  dKcdT[39] += -( dGkdT[13] + dGkdT[3] - ( dGkdT[14] + dGkdT[5] ) ) / calR / T;
  dKcdT[40] += -( dGkdT[13] + dGkdT[6] - ( dGkdT[14] + dGkdT[1] ) ) / calR / T;
  dKcdT[41] += -( dGkdT[13] + dGkdT[4] - ( dGkdT[14] ) ) / calR / T;
  dKcdT[42] += -( dGkdT[13] + dGkdT[12] - ( dGkdT[15] + dGkdT[1] ) ) / calR / T;
  dKcdT[43] += -( dGkdT[15] + dGkdT[2] - ( dGkdT[16] + dGkdT[4] ) ) / calR / T;
  dKcdT[44] += -( dGkdT[15] + dGkdT[3] - ( dGkdT[16] + dGkdT[5] ) ) / calR / T;
  dKcdT[45] += -( dGkdT[16] + dGkdT[2] - ( dGkdT[17] + dGkdT[4] ) ) / calR / T;
  dKcdT[46] += -( dGkdT[16] + dGkdT[6] - ( dGkdT[17] + dGkdT[1] ) ) / calR / T;
  dKcdT[47] += -( dGkdT[16] + dGkdT[4] - ( dGkdT[17] ) ) / calR / T;
  dKcdT[48] += -( dGkdT[9] + dGkdT[9] - ( dGkdT[18] ) ) / calR / T;
  dKcdT[49] += -( dGkdT[17] + dGkdT[2] - ( dGkdT[18] + dGkdT[4] ) ) / calR / T;
  dKcdT[50] += -( dGkdT[17] + dGkdT[3] - ( dGkdT[18] + dGkdT[5] ) ) / calR / T;
  dKcdT[51] += -( dGkdT[17] + dGkdT[10] - ( dGkdT[18] + dGkdT[9] ) ) / calR / T;

  for ( r = 0; r < nr; r++ ) {
    dKcdT[r] = dKcdT[r] * Kc[r];
  }

  dKcdT[4] = dKcdT[4] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[4] * 1E6 * calR / THERMO_P_REF;
  dKcdT[5] = dKcdT[5] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[5] * 1E6 * calR / THERMO_P_REF;
  dKcdT[6] = dKcdT[6] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[6] * 1E6 * calR / THERMO_P_REF;
  dKcdT[11] =
    dKcdT[11] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[11] / ( 1E6 * calR / THERMO_P_REF * sqr ( T ) );
  dKcdT[13] = dKcdT[13] * ( 1E6 * calR / THERMO_P_REF * T ) + Kc[13] * 1E6 * calR / THERMO_P_REF;
  dKcdT[18] =
    dKcdT[18] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[18] / ( 1E6 * calR / THERMO_P_REF * sqr ( T ) );
  dKcdT[24] =
    dKcdT[24] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[24] / ( 1E6 * calR / THERMO_P_REF * sqr ( T ) );
  dKcdT[32] =
    dKcdT[32] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[32] / ( 1E6 * calR / THERMO_P_REF * sqr ( T ) );
  dKcdT[41] =
    dKcdT[41] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[41] / ( 1E6 * calR / THERMO_P_REF * sqr ( T ) );
  dKcdT[47] =
    dKcdT[47] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[47] / ( 1E6 * calR / THERMO_P_REF * sqr ( T ) );
  dKcdT[48] =
    dKcdT[48] / ( 1E6 * calR / THERMO_P_REF * T ) - Kc[48] / ( 1E6 * calR / THERMO_P_REF * sqr ( T ) );

  Kc[4] = Kc[4] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[5] = Kc[5] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[6] = Kc[6] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[11] = Kc[11] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[13] = Kc[13] * ( 1E6 * calR / THERMO_P_REF * T );
  Kc[18] = Kc[18] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[24] = Kc[24] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[32] = Kc[32] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[41] = Kc[41] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[47] = Kc[47] / ( 1E6 * calR / THERMO_P_REF * T );
  Kc[48] = Kc[48] / ( 1E6 * calR / THERMO_P_REF * T );

}

void find_dW_dx ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long r, k, p;
  double kf[nr], kb[nr], react[nr];
  react_t Kc, dKcdT;
  double X[ns];
  double dWsdX[ns][ns];
  double third, third5, third11, third13, third18, third24, third41;

  T = max ( T, Tmin );

  /* initialize derivatives to zero */
  for ( r = 0; r < ns; r++ ) {
    dWdT[r] = 0.0;
    dWdTe[r] = 0.0;
    dWdTv[r] = 0.0;
    dWdQbeam[r] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[r][k] = 0.0;
    }
  }

  third = 0.0;
  for ( k = 0; k < ns; k++ ) {
    X[k] = _Xk ( k, rhok );
    third = third + X[k];
  }

  third5 =
    third + 0.9E+0 * X[2] + 1.6E+0 * X[1] + 1.6E+0 * X[19] + 8.5E+0 * X[3] + 1.6E+0 * X[7] + 1.6E+0 * X[8];
  third11 =
    third + 1.9E+0 * X[2] + 0.2E+0 * X[1] + 0.2E+0 * X[19] + 17.5E+0 * X[3] + 1.1E+0 * X[7] + 3.3E+0 * X[8];
  third13 = third + 11.0E+0 * X[1] + 1.0E+0 * X[19] + 2.0E+0 * X[7] + 6.0E+0 * X[8];
  third18 = third + 0.87E+0 * X[2] + 7.12E+0 * X[3];
  third24 =
    third + 1.9E+0 * X[2] + 0.2E+0 * X[1] + 0.2E+0 * X[19] + 17.5E+0 * X[3] + 1.1E+0 * X[7] + 3.3E+0 * X[8];
  third41 =
    third + 1.9E+0 * X[2] + 0.2E+0 * X[1] + 0.2E+0 * X[19] + 17.5E+0 * X[3] + 1.1E+0 * X[7] + 3.3E+0 * X[8];

  FinddKcdT ( T, Kc, dKcdT );

  for ( r = 0; r < nr; r++ ) {
    kf[r] = _kf ( r, third, T );
    kb[r] = kf[r] / Kc[r];
  }

  react[0] = _dkfdT ( 0, third, T ) * ( X[4] * X[1] - ( X[5] * X[0] ) / Kc[0] )
    + kf[0] / sqr ( Kc[0] ) * X[5] * X[0] * dKcdT[0];
  react[1] = _dkfdT ( 1, third, T ) * ( X[0] * X[2] - ( X[5] * X[4] ) / Kc[1] )
    + kf[1] / sqr ( Kc[1] ) * X[5] * X[4] * dKcdT[1];
  react[2] = _dkfdT ( 2, third, T ) * ( X[5] * X[2] - ( X[3] * X[4] ) / Kc[2] )
    + kf[2] / sqr ( Kc[2] ) * X[3] * X[4] * dKcdT[2];
  react[3] = _dkfdT ( 3, third, T ) * ( X[5] * X[5] - ( X[0] * X[3] ) / Kc[3] )
    + kf[3] / sqr ( Kc[3] ) * X[0] * X[3] * dKcdT[3];
  react[4] = _dkfdT ( 4, third, T ) * ( X[4] * X[4] - ( X[2] ) / Kc[4] ) * third
    + kf[4] / sqr ( Kc[4] ) * X[2] * dKcdT[4] * third;
  react[5] = _dkfdT ( 5, third, T ) * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * third5
    + kf[5] / sqr ( Kc[5] ) * X[3] * dKcdT[5] * third5;
  react[6] = _dkfdT ( 6, third, T ) * ( X[4] * X[1] - ( X[6] ) / Kc[6] ) * third
    + kf[6] / sqr ( Kc[6] ) * X[6] * dKcdT[6] * third;
  react[7] = _dkfdT ( 7, third, T ) * ( X[6] * X[4] - ( X[5] * X[5] ) / Kc[7] )
    + kf[7] / sqr ( Kc[7] ) * X[5] * X[5] * dKcdT[7];
  react[8] = _dkfdT ( 8, third, T ) * ( X[6] * X[4] - ( X[2] * X[1] ) / Kc[8] )
    + kf[8] / sqr ( Kc[8] ) * X[2] * X[1] * dKcdT[8];
  react[9] = _dkfdT ( 9, third, T ) * ( X[6] * X[0] - ( X[1] * X[5] ) / Kc[9] )
    + kf[9] / sqr ( Kc[9] ) * X[1] * X[5] * dKcdT[9];
  react[10] = _dkfdT ( 10, third, T ) * ( X[6] * X[5] - ( X[3] * X[1] ) / Kc[10] )
    + kf[10] / sqr ( Kc[10] ) * X[3] * X[1] * dKcdT[10];
  react[11] = _dkfdT ( 11, third, T ) * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * third11
    + kf[11] / sqr ( Kc[11] ) * X[5] * X[5] * dKcdT[11] * third11;
  react[12] = _dkfdT ( 12, third, T ) * ( X[7] * X[5] - ( X[8] * X[4] ) / Kc[12] )
    + kf[12] / sqr ( Kc[12] ) * X[8] * X[4] * dKcdT[12];
  react[13] = _dkfdT ( 13, third, T ) * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * third13
    + kf[13] / sqr ( Kc[13] ) * X[8] * dKcdT[13] * third13;
  react[14] = _dkfdT ( 14, third, T ) * ( X[12] * X[4] - ( X[7] * X[2] ) / Kc[14] )
    + kf[14] / sqr ( Kc[14] ) * X[7] * X[2] * dKcdT[14];
  react[15] = _dkfdT ( 15, third, T ) * ( X[12] * X[0] - ( X[7] * X[5] ) / Kc[15] )
    + kf[15] / sqr ( Kc[15] ) * X[7] * X[5] * dKcdT[15];
  react[16] = _dkfdT ( 16, third, T ) * ( X[12] * X[5] - ( X[7] * X[3] ) / Kc[16] )
    + kf[16] / sqr ( Kc[16] ) * X[7] * X[3] * dKcdT[16];
  react[17] = _dkfdT ( 17, third, T ) * ( X[12] * X[1] - ( X[7] * X[6] ) / Kc[17] )
    + kf[17] / sqr ( Kc[17] ) * X[7] * X[6] * dKcdT[17];
  react[18] = _dkfdT ( 18, third, T ) * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * third18
    + kf[18] / sqr ( Kc[18] ) * X[7] * X[4] * dKcdT[18] * third18;
  react[19] = _dkfdT ( 19, third, T ) * ( X[13] * X[4] - ( X[12] * X[2] ) / Kc[19] )
    + kf[19] / sqr ( Kc[19] ) * X[12] * X[2] * dKcdT[19];
  react[20] = _dkfdT ( 20, third, T ) * ( X[13] * X[0] - ( X[12] * X[5] ) / Kc[20] )
    + kf[20] / sqr ( Kc[20] ) * X[12] * X[5] * dKcdT[20];
  react[21] = _dkfdT ( 21, third, T ) * ( X[13] * X[5] - ( X[12] * X[3] ) / Kc[21] )
    + kf[21] / sqr ( Kc[21] ) * X[12] * X[3] * dKcdT[21];
  react[22] = _dkfdT ( 22, third, T ) * ( X[13] * X[1] - ( X[12] * X[6] ) / Kc[22] )
    + kf[22] / sqr ( Kc[22] ) * X[12] * X[6] * dKcdT[22];
  react[23] = _dkfdT ( 23, third, T ) * ( X[13] * X[9] - ( X[12] * X[10] ) / Kc[23] )
    + kf[23] / sqr ( Kc[23] ) * X[12] * X[10] * dKcdT[23];
  react[24] = _dkfdT ( 24, third, T ) * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * third24
    + kf[24] / sqr ( Kc[24] ) * X[12] * X[4] * dKcdT[24] * third24;
  react[25] = _dkfdT ( 25, third, T ) * ( X[9] * X[0] - ( X[13] * X[4] ) / Kc[25] )
    + kf[25] / sqr ( Kc[25] ) * X[13] * X[4] * dKcdT[25];
  react[26] = _dkfdT ( 26, third, T ) * ( X[9] * X[5] - ( X[13] * X[2] ) / Kc[26] )
    + kf[26] / sqr ( Kc[26] ) * X[13] * X[2] * dKcdT[26];
  react[27] = _dkfdT ( 27, third, T ) * ( X[9] * X[1] - ( X[14] * X[0] ) / Kc[27] )
    + kf[27] / sqr ( Kc[27] ) * X[14] * X[0] * dKcdT[27];
  react[28] = _dkfdT ( 28, third, T ) * ( X[9] * X[1] - ( X[13] * X[5] ) / Kc[28] )
    + kf[28] / sqr ( Kc[28] ) * X[13] * X[5] * dKcdT[28];
  react[29] = _dkfdT ( 29, third, T ) * ( X[9] * X[6] - ( X[14] * X[5] ) / Kc[29] )
    + kf[29] / sqr ( Kc[29] ) * X[14] * X[5] * dKcdT[29];
  react[30] = _dkfdT ( 30, third, T ) * ( X[9] * X[12] - ( X[10] * X[7] ) / Kc[30] )
    + kf[30] / sqr ( Kc[30] ) * X[10] * X[7] * dKcdT[30];
  react[31] = _dkfdT ( 31, third, T ) * ( X[9] * X[9] - ( X[17] * X[4] ) / Kc[31] )
    + kf[31] / sqr ( Kc[31] ) * X[17] * X[4] * dKcdT[31];
  react[32] = _dkfdT ( 32, third, T ) * ( X[10] - ( X[9] * X[4] ) / Kc[32] )
    + kf[32] / sqr ( Kc[32] ) * X[9] * X[4] * dKcdT[32];
  react[33] = _dkfdT ( 33, third, T ) * ( X[10] * X[4] - ( X[9] * X[2] ) / Kc[33] )
    + kf[33] / sqr ( Kc[33] ) * X[9] * X[2] * dKcdT[33];
  react[34] = _dkfdT ( 34, third, T ) * ( X[10] * X[0] - ( X[9] * X[5] ) / Kc[34] )
    + kf[34] / sqr ( Kc[34] ) * X[9] * X[5] * dKcdT[34];
  react[35] = _dkfdT ( 35, third, T ) * ( X[10] * X[1] - ( X[9] * X[6] ) / Kc[35] )
    + kf[35] / sqr ( Kc[35] ) * X[9] * X[6] * dKcdT[35];
  react[36] = _dkfdT ( 36, third, T ) * ( X[10] * X[5] - ( X[9] * X[3] ) / Kc[36] )
    + kf[36] / sqr ( Kc[36] ) * X[9] * X[3] * dKcdT[36];
  react[37] = _dkfdT ( 37, third, T ) * ( X[10] * X[6] - ( X[9] * X[11] ) / Kc[37] )
    + kf[37] / sqr ( Kc[37] ) * X[9] * X[11] * dKcdT[37];
  react[38] = _dkfdT ( 38, third, T ) * ( X[14] * X[4] - ( X[13] * X[2] ) / Kc[38] )
    + kf[38] / sqr ( Kc[38] ) * X[13] * X[2] * dKcdT[38];
  react[39] = _dkfdT ( 39, third, T ) * ( X[14] * X[5] - ( X[13] * X[3] ) / Kc[39] )
    + kf[39] / sqr ( Kc[39] ) * X[13] * X[3] * dKcdT[39];
  react[40] = _dkfdT ( 40, third, T ) * ( X[14] * X[1] - ( X[13] * X[6] ) / Kc[40] )
    + kf[40] / sqr ( Kc[40] ) * X[13] * X[6] * dKcdT[40];
  react[41] = _dkfdT ( 41, third, T ) * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * third41
    + kf[41] / sqr ( Kc[41] ) * X[13] * X[4] * dKcdT[41] * third41;
  react[42] = _dkfdT ( 42, third, T ) * ( X[15] * X[1] - ( X[13] * X[12] ) / Kc[42] )
    + kf[42] / sqr ( Kc[42] ) * X[13] * X[12] * dKcdT[42];
  react[43] = _dkfdT ( 43, third, T ) * ( X[16] * X[4] - ( X[15] * X[2] ) / Kc[43] )
    + kf[43] / sqr ( Kc[43] ) * X[15] * X[2] * dKcdT[43];
  react[44] = _dkfdT ( 44, third, T ) * ( X[16] * X[5] - ( X[15] * X[3] ) / Kc[44] )
    + kf[44] / sqr ( Kc[44] ) * X[15] * X[3] * dKcdT[44];
  react[45] = _dkfdT ( 45, third, T ) * ( X[17] * X[4] - ( X[16] * X[2] ) / Kc[45] )
    + kf[45] / sqr ( Kc[45] ) * X[16] * X[2] * dKcdT[45];
  react[46] = _dkfdT ( 46, third, T ) * ( X[17] * X[1] - ( X[16] * X[6] ) / Kc[46] )
    + kf[46] / sqr ( Kc[46] ) * X[16] * X[6] * dKcdT[46];
  react[47] = _dkfdT ( 47, third, T ) * ( X[17] - ( X[16] * X[4] ) / Kc[47] )
    + kf[47] / sqr ( Kc[47] ) * X[16] * X[4] * dKcdT[47];
  react[48] = _dkfdT ( 48, third, T ) * ( X[18] - ( X[9] * X[9] ) / Kc[48] )
    + kf[48] / sqr ( Kc[48] ) * X[9] * X[9] * dKcdT[48];
  react[49] = _dkfdT ( 49, third, T ) * ( X[18] * X[4] - ( X[17] * X[2] ) / Kc[49] )
    + kf[49] / sqr ( Kc[49] ) * X[17] * X[2] * dKcdT[49];
  react[50] = _dkfdT ( 50, third, T ) * ( X[18] * X[5] - ( X[17] * X[3] ) / Kc[50] )
    + kf[50] / sqr ( Kc[50] ) * X[17] * X[3] * dKcdT[50];
  react[51] = _dkfdT ( 51, third, T ) * ( X[18] * X[9] - ( X[17] * X[10] ) / Kc[51] )
    + kf[51] / sqr ( Kc[51] ) * X[17] * X[10] * dKcdT[51];

  dWdT[0] =
    ( react[0] - react[1] + react[3] - react[9] - react[13] - react[15] - react[20] - react[25] + react[27]
      - react[34] );
  dWdT[1] =
    ( -react[0] - react[6] + react[8] + react[9] + react[10] - react[17] - react[22] - react[27] - react[28]
      - react[35] - react[40] - react[42] - react[46] );
  dWdT[2] =
    ( -react[1] - react[2] + react[4] + react[8] + react[14] + react[19] + react[26] + react[33] + react[38]
      + react[43] + react[45] + react[49] );
  dWdT[3] =
    ( react[2] + react[3] + react[5] + react[10] + react[16] + react[21] + react[36] + react[39] + react[44]
      + react[50] );
  dWdT[4] =
    ( -react[0] + react[1] + react[2] - react[4] - react[4] - react[5] - react[6] - react[7] - react[8]
      + react[12] - react[14] + react[18] - react[19] + react[24] + react[25] + react[31] + react[32] -
      react[33]
      - react[38] + react[41] - react[43] - react[45] + react[47] - react[49] );
  dWdT[5] = ( react[0] + react[1] - react[2] - react[3] - react[3] - react[5] + react[7] + react[7] + react[9]
              - react[10] + react[11] + react[11] - react[12] + react[15] - react[16] + react[20] -
              react[21] - react[26]
              + react[28] + react[29] + react[34] - react[36] - react[39] - react[44] - react[50] );
  dWdT[6] =
    ( react[6] - react[7] - react[8] - react[9] - react[10] + react[17] + react[22] - react[29] + react[35]
      - react[37] + react[40] + react[46] );
  dWdT[7] =
    ( -react[12] - react[13] + react[14] + react[15] + react[16] + react[17] + react[18] + react[30] );
  dWdT[8] = ( react[12] + react[13] );
  dWdT[9] =
    ( -react[23] - react[25] - react[26] - react[27] - react[28] - react[29] - react[30] - react[31] -
      react[31]
      + react[32] + react[33] + react[34] + react[35] + react[36] + react[37] + react[48] + react[48] -
      react[51] );
  dWdT[10] =
    ( react[23] + react[30] - react[32] - react[33] - react[34] - react[35] - react[36] - react[37] +
      react[51] );
  dWdT[11] = ( -react[11] + react[37] );
  dWdT[12] =
    ( -react[14] - react[15] - react[16] - react[17] - react[18] + react[19] + react[20] + react[21] +
      react[22]
      + react[23] + react[24] - react[30] + react[42] );
  dWdT[13] =
    ( -react[19] - react[20] - react[21] - react[22] - react[23] - react[24] + react[25] + react[26] +
      react[28]
      + react[38] + react[39] + react[40] + react[41] + react[42] );
  dWdT[14] = ( react[27] + react[29] - react[38] - react[39] - react[40] - react[41] );
  dWdT[15] = ( -react[42] + react[43] + react[44] );
  dWdT[16] = ( -react[43] - react[44] + react[45] + react[46] + react[47] );
  dWdT[17] = ( react[31] - react[45] - react[46] - react[47] + react[49] + react[50] + react[51] );
  dWdT[18] = ( -react[48] - react[49] - react[50] - react[51] );
  dWdT[19] = 0.0E+0;

  for ( k = 0; k < ns; k++ )
    dWdT[k] = dWdT[k] * 1.0e6 * _calM ( k );

/*
  W[0] = ( react[0] -react[1] +react[3] -react[9] -react[13]-react[15]-react[20]-react[25]+react[27]
          -react[34]) 
        *1.0E+6*_calM(0);
  W[1] = (-react[0] -react[6] +react[8] +react[9] +react[10]-react[17]-react[22]-react[27]-react[28] 
          -react[35]-react[40]-react[42]-react[46])
        *1.0E+6*_calM(1);
  W[2] = (-react[1] -react[2] +react[4] +react[8] +react[14]+react[19]+react[26]+react[33]+react[38]
          +react[43]+react[45]+react[49]) 
        *1.0E+6*_calM(2);
  W[3] = ( react[2] +react[3] +react[5] +react[10]+react[16]+react[21]+react[36]+react[39]+react[44]
          +react[50]) 
        *1.0E+6*_calM(3);
  W[4] = (-react[0] +react[1] +react[2] -react[4] -react[4] -react[5] -react[6] -react[7] -react[8] 
          +react[12]-react[14]+react[18]-react[19]+react[24]+react[25]+react[31]+react[32]-react[33] 
          -react[38]+react[41]-react[43]-react[45]+react[47]-react[49]) 
        *1.0E+6*_calM(4);
  W[5] = ( react[0] +react[1] -react[2] -react[3] -react[3] -react[5] +react[7] +react[7] +react[9] 
          -react[10]+react[11]+react[11]-react[12]+react[15]-react[16]+react[20]-react[21]-react[26] 
          +react[28]+react[29]+react[34]-react[36]-react[39]-react[44]-react[50]) 
        *1.0E+6*_calM(5);
  W[6] = ( react[6] -react[7] -react[8] -react[9] -react[10]+react[17]+react[22]-react[29]+react[35]
          -react[37]+react[40]+react[46]) 
        *1.0E+6*_calM(6);
  W[7] = (-react[12]-react[13]+react[14]+react[15]+react[16]+react[17]+react[18]+react[30]) 
        *1.0E+6*_calM(7);
  W[8] = ( react[12]+react[13]) 
        *1.0E+6*_calM(8);
  W[9] = (-react[23]-react[25]-react[26]-react[27]-react[28]-react[29]-react[30]-react[31]-react[31]
          +react[32]+react[33]+react[34]+react[35]+react[36]+react[37]+react[48]+react[48]-react[51])
        *1.0E+6*_calM(9);
  W[10]= ( react[23]+react[30]-react[32]-react[33]-react[34]-react[35]-react[36]-react[37]+react[51])
        *1.0E+6*_calM(10);
  W[11]= (-react[11]+react[37]) 
        *1.0E+6*_calM(11);
  W[12]= (-react[14]-react[15]-react[16]-react[17]-react[18]+react[19]+react[20]+react[21]+react[22]
          +react[23]+react[24]-react[30]+react[42]) 
        *1.0E+6*_calM(12);
  W[13]= (-react[19]-react[20]-react[21]-react[22]-react[23]-react[24]+react[25]+react[26]+react[28] 
          +react[38]+react[39]+react[40]+react[41]+react[42]) 
        *1.0E+6*_calM(13);
  W[14]= ( react[27]+react[29]-react[38]-react[39]-react[40]-react[41])
        *1.0E+6*_calM(14);
  W[15]= (-react[42]+react[43]+react[44]) 
        *1.0E+6*_calM(15);
  W[16]= (-react[43]-react[44]+react[45]+react[46]+react[47]) 
        *1.0E+6*_calM(16);
  W[17]= ( react[31]-react[45]-react[46]-react[47]+react[49]+react[50]+react[51]) 
        *1.0E+6*_calM(17);
  W[18]= (-react[48]-react[49]-react[50]-react[51]) 
        *1.0E+6*_calM(18);

  react[0]  = kf[0] *(X[4] *X[1] -(X[5] *X[0] )/Kc[0]);
  react[1]  = kf[1] *(X[0] *X[2] -(X[5] *X[4] )/Kc[1]);
  react[2]  = kf[2] *(X[5] *X[2] -(X[3] *X[4] )/Kc[2]);  
  react[3]  = kf[3] *(X[5] *X[5] -(X[0] *X[3] )/Kc[3]);  
  react[4]  = kf[4] *(X[4] *X[4] -(X[2]       )/Kc[4])*third; 
  react[5]  = kf[5] *(X[4] *X[5] -(X[3]       )/Kc[5])*third5;  
  react[6]  = kf[6] *(X[4] *X[1] -(X[6]       )/Kc[6])*third;  
  react[7]  = kf[7] *(X[6] *X[4] -(X[5] *X[5] )/Kc[7]);  
  react[8]  = kf[8] *(X[6] *X[4] -(X[2] *X[1] )/Kc[8]);  
  react[9]  = kf[9] *(X[6] *X[0] -(X[1] *X[5] )/Kc[9]);  
  react[10] = kf[10]*(X[6] *X[5] -(X[3] *X[1] )/Kc[10]);  
  react[11] = kf[11]*(X[11]      -(X[5] *X[5] )/Kc[11])*third11;  
  react[12] = kf[12]*(X[7] *X[5] -(X[8] *X[4] )/Kc[12]);  
  react[13] = kf[13]*(X[7] *X[0] -(X[8]       )/Kc[13])*third13;  
  react[14] = kf[14]*(X[12]*X[4] -(X[7] *X[2] )/Kc[14]);  
  react[15] = kf[15]*(X[12]*X[0] -(X[7] *X[5] )/Kc[15]);  
  react[16] = kf[16]*(X[12]*X[5] -(X[7] *X[3] )/Kc[16]);  
  react[17] = kf[17]*(X[12]*X[1] -(X[7] *X[6] )/Kc[17]);  
  react[18] = kf[18]*(X[12]      -(X[7] *X[4] )/Kc[18])*third18;  
  react[19] = kf[19]*(X[13]*X[4] -(X[12]*X[2] )/Kc[19]);  
  react[20] = kf[20]*(X[13]*X[0] -(X[12]*X[5] )/Kc[20]);    
  react[21] = kf[21]*(X[13]*X[5] -(X[12]*X[3] )/Kc[21]);  
  react[22] = kf[22]*(X[13]*X[1] -(X[12]*X[6] )/Kc[22]);  
  react[23] = kf[23]*(X[13]*X[9] -(X[12]*X[10])/Kc[23]);  
  react[24] = kf[24]*(X[13]      -(X[12]*X[4] )/Kc[24])*third24;  
  react[25] = kf[25]*(X[9] *X[0] -(X[13]*X[4] )/Kc[25]);  
  react[26] = kf[26]*(X[9] *X[5] -(X[13]*X[2] )/Kc[26]);  
  react[27] = kf[27]*(X[9] *X[1] -(X[14]*X[0] )/Kc[27]);  
  react[28] = kf[28]*(X[9] *X[1] -(X[13]*X[5] )/Kc[28]);  
  react[29] = kf[29]*(X[9] *X[6] -(X[14]*X[5] )/Kc[29]);  
  react[30] = kf[30]*(X[9] *X[12]-(X[10]*X[7] )/Kc[30]);  
  react[31] = kf[31]*(X[9] *X[9] -(X[17]*X[4] )/Kc[31]);    
  react[32] = kf[32]*(X[10]      -(X[9] *X[4] )/Kc[32]);
  react[33] = kf[33]*(X[10]*X[4] -(X[9] *X[2] )/Kc[33]);
  react[34] = kf[34]*(X[10]*X[0] -(X[9] *X[5] )/Kc[34]);
  react[35] = kf[35]*(X[10]*X[1] -(X[9] *X[6] )/Kc[35]);
  react[36] = kf[36]*(X[10]*X[5] -(X[9] *X[3] )/Kc[36]);
  react[37] = kf[37]*(X[10]*X[6] -(X[9] *X[11])/Kc[37]);    
  react[38] = kf[38]*(X[14]*X[4] -(X[13]*X[2] )/Kc[38]);  
  react[39] = kf[39]*(X[14]*X[5] -(X[13]*X[3] )/Kc[39]);  
  react[40] = kf[40]*(X[14]*X[1] -(X[13]*X[6] )/Kc[40]);  
  react[41] = kf[41]*(X[14]      -(X[13]*X[4] )/Kc[41])*third41;  
  react[42] = kf[42]*(X[15]*X[1] -(X[13]*X[12])/Kc[42]);  
  react[43] = kf[43]*(X[16]*X[4] -(X[15]*X[2] )/Kc[43]);  
  react[44] = kf[44]*(X[16]*X[5] -(X[15]*X[3] )/Kc[44]);  
  react[45] = kf[45]*(X[17]*X[4] -(X[16]*X[2] )/Kc[45]);  
  react[46] = kf[46]*(X[17]*X[1] -(X[16]*X[6] )/Kc[46]);  
  react[47] = kf[47]*(X[17]      -(X[16]*X[4] )/Kc[47]);  
  react[48] = kf[48]*(X[18]      -(X[9] *X[9] )/Kc[48]);  
  react[49] = kf[49]*(X[18]*X[4] -(X[17]*X[2] )/Kc[49]);  
  react[50] = kf[50]*(X[18]*X[5] -(X[17]*X[3] )/Kc[50]);  
  react[51] = kf[51]*(X[18]*X[9] -(X[17]*X[10])/Kc[51]);

  third5  = third+ 0.9E+0*X[2]+ 1.6E+0*X[1]+ 1.6E+0*X[19]+ 8.5E+0*X[3]+1.6E+0*X[7]+1.6E+0*X[8];
  third11 = third+ 1.9E+0*X[2]+ 0.2E+0*X[1]+ 0.2E+0*X[19]+17.5E+0*X[3]+1.1E+0*X[7]+3.3E+0*X[8];
  third13 = third+11.0E+0*X[1]+ 1.0E+0*X[19]+2.0E+0*X[7]+  6.0E+0*X[8];
  third18 = third+0.87E+0*X[2]+7.12E+0*X[3];
  third24 = third+ 1.9E+0*X[2]+ 0.2E+0*X[1]+ 0.2E+0*X[19]+17.5E+0*X[3]+1.1E+0*X[7]+3.3E+0*X[8];
  third41 = third+ 1.9E+0*X[2]+ 0.2E+0*X[1]+ 0.2E+0*X[19]+17.5E+0*X[3]+1.1E+0*X[7]+3.3E+0*X[8];

*/

/*

  W[0] = ( react[0] -react[1] +react[3] -react[9] -react[13]-react[15]-react[20]-react[25]+react[27]
          -react[34]) 
        *1.0E+6*_calM(0);

???? 32,47,48 double _dkfdthird(long r, double third, double T){
*/

/* find the derivative of each source term wrt concentration */

  for ( k = 0; k < ns; k++ ) {
    for ( p = 0; p < ns; p++ ) {
      dWsdX[k][p] = 0.0e0;
    }
  }

  dWsdX[0][0] =
    -kf[0] * X[5] / Kc[0] - kf[1] * X[2] - kf[3] * X[3] / Kc[3] - kf[9] * X[6] - kf[13] * X[7] * third13 -
    kf[15] * X[12] - kf[20] * X[13] - kf[25] * X[9] - kf[27] * X[14] / Kc[27] - kf[34] * X[10];

  dWsdX[0][1] =
    +kf[0] * X[4] + kf[9] * X[5] / Kc[9] + kf[27] * X[9] - kf[13] * ( X[7] * X[0] -
                                                                      ( X[8] ) / Kc[13] ) * 11.0e0;

  dWsdX[0][2] = ( -1.0 * kf[1] * X[0] );

  dWsdX[0][3] = -1.0 * kf[3] * X[0] / Kc[3];

  dWsdX[0][4] = ( 1.0 * kf[0] * X[1] + 1.0 * kf[1] * X[5] / Kc[1] + 1.0 * kf[25] * X[13] / Kc[25] );

  dWsdX[0][5] = ( -1.0 * kb[0] * X[0] + 1.0 * kb[1] * X[4] + 2.0 * kf[3] * X[5]
                  + 1.0 * kb[9] * X[1] + 1.0 * kb[15] * X[7] + 1.0 * kb[20] * X[12]
                  + 1.0 * kb[34] * X[9] );

  dWsdX[0][6] = ( -1.0 * kf[9] * X[0] );

  dWsdX[0][7] =
    ( -1.0 * kf[13] * X[0] * third13 + 1.0 * kb[15] * X[5] ) - kf[13] * ( X[7] * X[0] -
                                                                          ( X[8] ) / Kc[13] ) * 2.0e0;

  dWsdX[0][8] = ( 1.0 * kb[13] * third13 ) - kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 6.0e0;

  dWsdX[0][9] = ( -1.0 * kf[25] * X[0] + 1.0 * kf[27] * X[1] + 1.0 * kb[34] * X[5] );

  dWsdX[0][10] = ( -1.0 * kf[34] * X[0] );

  dWsdX[0][12] = ( -1.0 * kf[15] * X[0] + 1.0 * kb[20] * X[5] );

  dWsdX[0][13] = ( -1.0 * kf[20] * X[0] + 1.0 * kb[25] * X[4] );

  dWsdX[0][14] = ( -1.0 * kb[27] * X[0] );

  dWsdX[0][19] = -kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 1.0e0;

  for ( k = 0; k < ns; k++ ) {
    dWsdX[0][k] += ( -1.0 * kf[13] * X[7] * X[0] + 1.0 * kb[13] * X[8] );
  }

  dWsdX[1][0] = ( 1.0 * kb[0] * X[5] + 1.0 * kf[9] * X[6] + 1.0 * kb[27] * X[14] );

  dWsdX[1][1] =
    ( -1.0 * kf[0] * X[4] - 1.0 * kf[6] * X[4] * third - 1.0 * kb[8] * X[2] - 1.0 * kb[9] * X[5] -
      1.0 * kb[10] * X[3]
      - 1.0 * kf[17] * X[12] - 1.0 * kf[22] * X[13] - 1.0 * kf[27] * X[9] - 1.0 * kf[28] * X[9] -
      1.0 * kf[35] * X[10]
      - 1.0 * kf[40] * X[14] - 1.0 * kf[42] * X[15] - 1.0 * kf[46] * X[17] );

  dWsdX[1][2] = ( -1.0 * kb[8] * X[1] );

  dWsdX[1][3] = ( -1.0 * kb[10] * X[1] );

  dWsdX[1][4] = ( -1.0 * kf[0] * X[1] - 1.0 * kf[6] * X[1] * third + 1.0 * kf[8] * X[6] );

  dWsdX[1][5] = ( 1.0 * kb[0] * X[0] - 1.0 * kb[9] * X[1] + 1.0 * kf[10] * X[6] + 1.0 * kb[28] * X[13] );

  dWsdX[1][6] =
    ( 1.0 * kb[6] * third + 1.0 * kf[8] * X[4] + 1.0 * kf[9] * X[0] + 1.0 * kf[10] * X[5] +
      1.0 * kb[17] * X[7]
      + 1.0 * kb[22] * X[12] + 1.0 * kb[35] * X[9] + 1.0 * kb[40] * X[13] + 1.0 * kb[46] * X[16] );

  dWsdX[1][7] = ( 1.0 * kb[17] * X[6] );

  dWsdX[1][9] = ( -1.0 * kf[27] * X[1] - 1.0 * kf[28] * X[1] + 1.0 * kb[35] * X[6] );

  dWsdX[1][10] = ( -1.0 * kf[35] * X[1] );

  dWsdX[1][12] = ( -1.0 * kf[17] * X[1] + 1.0 * kb[22] * X[6] + 1.0 * kb[42] * X[13] );

  dWsdX[1][13] = ( -1.0 * kf[22] * X[1] + 1.0 * kb[28] * X[5] + 1.0 * kb[40] * X[6] + 1.0 * kb[42] * X[12] );

  dWsdX[1][14] = ( 1.0 * kb[27] * X[0] - 1.0 * kf[40] * X[1] );

  dWsdX[1][15] = ( -1.0 * kf[42] * X[1] );

  dWsdX[1][16] = ( 1.0 * kb[46] * X[6] );

  dWsdX[1][17] = ( -1.0 * kf[46] * X[1] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[1][k] += ( -1.0 * kf[6] * X[4] * X[1] + 1.0 * kb[6] * X[6] );
  }

  dWsdX[2][0] = ( -1.0 * kf[1] * X[2] );

  dWsdX[2][1] = ( -1.0 * kb[8] * X[2] );

  dWsdX[2][2] =
    ( -1.0 * kf[1] * X[0] - 1.0 * kf[2] * X[5] - 1.0 * kb[4] * third - 1.0 * kb[8] * X[1] -
      1.0 * kb[14] * X[7]
      - 1.0 * kb[19] * X[12] - 1.0 * kb[26] * X[13] - 1.0 * kb[33] * X[9] - 1.0 * kb[38] * X[13]
      - 1.0 * kb[43] * X[15] - 1.0 * kb[45] * X[16] - 1.0 * kb[49] * X[17] );

  dWsdX[2][3] = ( 1.0 * kb[2] * X[4] );

  dWsdX[2][4] =
    ( 1.0 * kb[1] * X[5] + 1.0 * kb[2] * X[3] + 2.0 * kf[4] * X[4] * third + 1.0 * kf[8] * X[6] +
      1.0 * kf[14] * X[12]
      + 1.0 * kf[19] * X[13] + 1.0 * kf[33] * X[10] + 1.0 * kf[38] * X[14] + 1.0 * kf[43] * X[16]
      + 1.0 * kf[45] * X[17] + 1.0 * kf[49] * X[18] );

  dWsdX[2][5] = ( 1.0 * kb[1] * X[4] - 1.0 * kf[2] * X[2] + 1.0 * kf[26] * X[9] );

  dWsdX[2][6] = ( 1.0 * kf[8] * X[4] );

  dWsdX[2][7] = ( -1.0 * kb[14] * X[2] );

  dWsdX[2][9] = ( 1.0 * kf[26] * X[5] - 1.0 * kb[33] * X[2] );

  dWsdX[2][10] = ( 1.0 * kf[33] * X[4] );

  dWsdX[2][12] = ( 1.0 * kf[14] * X[4] - 1.0 * kb[19] * X[2] );

  dWsdX[2][13] = ( 1.0 * kf[19] * X[4] - 1.0 * kb[26] * X[2] - 1.0 * kb[38] * X[2] );

  dWsdX[2][14] = ( 1.0 * kf[38] * X[4] );

  dWsdX[2][15] = ( -1.0 * kb[43] * X[2] );

  dWsdX[2][16] = ( 1.0 * kf[43] * X[4] - 1.0 * kb[45] * X[2] );

  dWsdX[2][17] = ( 1.0 * kf[45] * X[4] - 1.0 * kb[49] * X[2] );

  dWsdX[2][18] = ( 1.0 * kf[49] * X[4] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[2][k] += ( 1.0 * kf[4] * X[4] * X[4] - 1.0 * kb[4] * X[2] );
  }

  dWsdX[3][0] = ( -1.0 * kb[3] * X[3] );

  dWsdX[3][1] = ( -1.0 * kf[10] * X[3] / Kc[10] ) + 1.6e0 * kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] );

  dWsdX[3][2] = ( 1.0 * kf[2] * X[5] ) + 0.9e0 * kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] );

  dWsdX[3][3] =
    ( -1.0 * kb[2] * X[4] - 1.0 * kb[3] * X[0] - 1.0 * kb[5] * third5 - 1.0 * kb[10] * X[1] -
      1.0 * kb[16] * X[7]
      - 1.0 * kb[21] * X[12] - 1.0 * kb[36] * X[9] - 1.0 * kb[39] * X[13] - 1.0 * kb[44] * X[15]
      - 1.0 * kb[50] * X[17] )
    + 8.5e0 * kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] );

  dWsdX[3][4] = ( -1.0 * kb[2] * X[3] + 1.0 * kf[5] * X[5] * third5 );

  dWsdX[3][5] = ( 1.0 * kf[2] * X[2] + 2.0 * kf[3] * X[5] + 1.0 * kf[5] * X[4] * third5 + 1.0 * kf[10] * X[6]
                  + 1.0 * kf[16] * X[12] + 1.0 * kf[21] * X[13] + 1.0 * kf[36] * X[10] + 1.0 * kf[39] * X[14]
                  + 1.0 * kf[44] * X[16] + 1.0 * kf[50] * X[18] );

  dWsdX[3][6] = ( 1.0 * kf[10] * X[5] );

  dWsdX[3][7] = ( -1.0 * kb[16] * X[3] ) + 1.6e0 * kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] );

  dWsdX[3][8] = +1.6e0 * kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] );

  dWsdX[3][9] = ( -1.0 * kb[36] * X[3] );

  dWsdX[3][10] = ( 1.0 * kf[36] * X[5] );

  dWsdX[3][12] = ( 1.0 * kf[16] * X[5] - 1.0 * kb[21] * X[3] );

  dWsdX[3][13] = ( 1.0 * kf[21] * X[5] - 1.0 * kb[39] * X[3] );

  dWsdX[3][14] = ( 1.0 * kf[39] * X[5] );

  dWsdX[3][15] = ( -1.0 * kb[44] * X[3] );

  dWsdX[3][16] = ( 1.0 * kf[44] * X[5] );

  dWsdX[3][17] = ( -1.0 * kb[50] * X[3] );

  dWsdX[3][18] = ( 1.0 * kf[50] * X[5] );

  dWsdX[3][19] = +1.6e0 * kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[3][k] += ( 1.0 * kf[5] * X[4] * X[5] - 1.0 * kb[5] * X[3] );
  }

  dWsdX[4][0] = ( 1.0 * kb[0] * X[5] + 1.0 * kf[1] * X[2] + 1.0 * kf[25] * X[9] );

  dWsdX[4][1] = ( -1.0 * kf[0] * X[4] - 1.0 * kf[6] * X[4] * third + 1.0 * kb[8] * X[2] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 0.2e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 0.2e0;

  dWsdX[4][2] =
    ( 1.0 * kf[1] * X[0] + 1.0 * kf[2] * X[5] + 2.0 * kb[4] * third + 1.0 * kb[8] * X[1] + 1.0 * kb[14] * X[7]
      + 1.0 * kb[19] * X[12] + 1.0 * kb[33] * X[9] + 1.0 * kb[38] * X[13] + 1.0 * kb[43] * X[15]
      + 1.0 * kb[45] * X[16] + 1.0 * kb[49] * X[17] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 0.9e0
    + kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * 0.87e0
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 1.9e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 1.9e0;

  dWsdX[4][3] = ( -1.0 * kb[2] * X[4] + 1.0 * kb[5] * third5 )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 8.5e0
    + kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * 7.12e0
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 17.5e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 17.5e0;

  dWsdX[4][4] =
    ( -1.0 * kf[0] * X[1] - 1.0 * kb[1] * X[5] - 1.0 * kb[2] * X[3] - 4.0 * kf[4] * X[4] * third
      - 1.0 * kf[5] * X[5] * third5 - 1.0 * kf[6] * X[1] * third - 1.0 * kf[7] * X[6] - 1.0 * kf[8] * X[6]
      - 1.0 * kb[12] * X[8] - 1.0 * kf[14] * X[12] - 1.0 * kb[18] * X[7] * third18 - 1.0 * kf[19] * X[13]
      - 1.0 * kb[24] * X[12] * third24 - 1.0 * kb[25] * X[13] - 1.0 * kb[31] * X[17] - 1.0 * kb[32] * X[9]
      - 1.0 * kf[33] * X[10] - 1.0 * kf[38] * X[14] - 1.0 * kb[41] * X[13] * third41 - 1.0 * kf[43] * X[16]
      - 1.0 * kf[45] * X[17] - 1.0 * kb[47] * X[16] - 1.0 * kf[49] * X[18] );

  dWsdX[4][5] =
    ( 1.0 * kb[0] * X[0] - 1.0 * kb[1] * X[4] + 1.0 * kf[2] * X[2] - 1.0 * kf[5] * X[4] * third5 +
      2.0 * kb[7] * X[5]
      + 1.0 * kf[12] * X[7] );

  dWsdX[4][6] = ( 1.0 * kb[6] * third - 1.0 * kf[7] * X[4] - 1.0 * kf[8] * X[4] );

  dWsdX[4][7] =
    ( 1.0 * kf[12] * X[5] + 1.0 * kb[14] * X[2] -
      1.0 * kb[18] * X[4] * ( third + 0.87E0 * X[2] + 7.12E0 * X[3] ) )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0 + kf[24] * ( X[13] -
                                                                      ( X[12] * X[4] ) / Kc[24] ) * 1.1e0 +
    kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 1.1e0;

  dWsdX[4][8] = ( -1.0 * kb[12] * X[4] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 3.3e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 3.3e0;

  dWsdX[4][9] = ( 1.0 * kf[25] * X[0] + 2.0 * kf[31] * X[9] - 1.0 * kb[32] * X[4] + 1.0 * kb[33] * X[2] );

  dWsdX[4][10] = ( 1.0 * kf[32] - 1.0 * kf[33] * X[4] );

  dWsdX[4][12] =
    ( -1.0 * kf[14] * X[4] + 1.0 * kf[18] * third18 + 1.0 * kb[19] * X[2] - 1.0 * kb[24] * X[4] * third24 );

  dWsdX[4][13] = ( -1.0 * kf[19] * X[4] + 1.0 * kf[24] * third24 - 1.0 * kb[25] * X[4] + 1.0 * kb[38] * X[2]
                   - 1.0 * kb[41] * X[4] * third41 );

  dWsdX[4][14] = ( -1.0 * kf[38] * X[4] + 1.0 * kf[41] * third41 );

  dWsdX[4][15] = ( 1.0 * kb[43] * X[2] );

  dWsdX[4][16] = ( -1.0 * kf[43] * X[4] + 1.0 * kb[45] * X[2] - 1.0 * kb[47] * X[4] );

  dWsdX[4][17] = ( -1.0 * kb[31] * X[4] - 1.0 * kf[45] * X[4] + 1.0 * kf[47] + 1.0 * kb[49] * X[2] );

  dWsdX[4][18] = ( -1.0 * kf[49] * X[4] );

  dWsdX[4][19] = -kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 0.2e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 0.2e0;

  for ( k = 0; k < ns; k++ ) {
    dWsdX[4][k] +=
      -kf[4] * ( X[4] * X[4] - ( X[2] ) / Kc[4] ) * 2.0 - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] )
      - kf[6] * ( X[4] * X[1] - ( X[6] ) / Kc[6] )
      + kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] )
      + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] )
      + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] )
      + ( X[10] - ( X[9] * X[4] ) / Kc[32] ) * _dkfdthird ( 32, third, T )
      + ( X[17] - ( X[16] * X[4] ) / Kc[47] ) * _dkfdthird ( 47, third, T );
  }

  dWsdX[5][0] =
    ( -1.0 * kb[0] * X[5] + 1.0 * kf[1] * X[2] + 2.0 * kb[3] * X[3] + 1.0 * kf[9] * X[6] +
      1.0 * kf[15] * X[12]
      + 1.0 * kf[20] * X[13] + 1.0 * kf[34] * X[10] );

  dWsdX[5][1] = ( 1.0 * kf[0] * X[4] - 1.0 * kb[9] * X[5] + 1.0 * kb[10] * X[3] + 1.0 * kf[28] * X[9] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0
    + kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 2.0e0 * 0.2e0;

  dWsdX[5][2] = ( 1.0 * kf[1] * X[0] - 1.0 * kf[2] * X[5] + 1.0 * kb[26] * X[13] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 0.9e0
    + kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 2.0 * 1.9e0;

  dWsdX[5][3] =
    ( 1.0 * kb[2] * X[4] + 2.0 * kb[3] * X[0] + 1.0 * kb[5] * third5 + 1.0 * kb[10] * X[1] +
      1.0 * kb[16] * X[7]
      + 1.0 * kb[21] * X[12] + 1.0 * kb[36] * X[9] + 1.0 * kb[39] * X[13] + 1.0 * kb[44] * X[15]
      + 1.0 * kb[50] * X[17] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 8.5e0
    + kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 2.0 * 17.5e0;

  dWsdX[5][4] =
    ( 1.0 * kf[0] * X[1] - 1.0 * kb[1] * X[5] + 1.0 * kb[2] * X[3] - 1.0 * kf[5] * X[5] * third5 +
      2.0 * kf[7] * X[6]
      + 1.0 * kb[12] * X[8] );

  dWsdX[5][5] =
    ( -1.0 * kb[0] * X[0] - 1.0 * kb[1] * X[4] - 1.0 * kf[2] * X[2] - 4.0 * kf[3] * X[5] -
      1.0 * kf[5] * X[4] * third5 - 4.0 * kb[7] * X[5] - 1.0 * kb[9] * X[1] - 1.0 * kf[10] * X[6] -
      4.0 * kb[11] * X[5] * third11 - 1.0 * kf[12] * X[7] - 1.0 * kb[15] * X[7] - 1.0 * kf[16] * X[12] -
      1.0 * kb[20] * X[12] - 1.0 * kf[21] * X[13]
      - 1.0 * kf[26] * X[9] - 1.0 * kb[28] * X[13] - 1.0 * kb[29] * X[14] - 1.0 * kb[34] * X[9] -
      1.0 * kf[36] * X[10]
      - 1.0 * kf[39] * X[14] - 1.0 * kf[44] * X[16] - 1.0 * kf[50] * X[18] );

  dWsdX[5][6] = ( 2.0 * kf[7] * X[4] + 1.0 * kf[9] * X[0] - 1.0 * kf[10] * X[5] + 1.0 * kf[29] * X[9] );

  dWsdX[5][7] = ( -1.0 * kf[12] * X[5] - 1.0 * kb[15] * X[5] + 1.0 * kb[16] * X[3] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0
    + kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 2.0 * 1.1e0;

  dWsdX[5][8] = ( 1.0 * kb[12] * X[4] )
    - kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0
    + kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 2.0 * 3.3e0;

  dWsdX[5][9] =
    ( -1.0 * kf[26] * X[5] + 1.0 * kf[28] * X[1] + 1.0 * kf[29] * X[6] - 1.0 * kb[34] * X[5] +
      1.0 * kb[36] * X[3] );

  dWsdX[5][10] = ( 1.0 * kf[34] * X[0] - 1.0 * kf[36] * X[5] );

  dWsdX[5][11] = ( 2.0 * kf[11] * third11 );

  dWsdX[5][12] = ( 1.0 * kf[15] * X[0] - 1.0 * kf[16] * X[5] - 1.0 * kb[20] * X[5] + 1.0 * kb[21] * X[3] );

  dWsdX[5][13] =
    ( 1.0 * kf[20] * X[0] - 1.0 * kf[21] * X[5] + 1.0 * kb[26] * X[2] - 1.0 * kb[28] * X[5] +
      1.0 * kb[39] * X[3] );

  dWsdX[5][14] = ( -1.0 * kb[29] * X[5] - 1.0 * kf[39] * X[5] );

  dWsdX[5][15] = ( 1.0 * kb[44] * X[3] );

  dWsdX[5][16] = ( -1.0 * kf[44] * X[5] );

  dWsdX[5][17] = ( 1.0 * kb[50] * X[3] );

  dWsdX[5][18] = ( -1.0 * kf[50] * X[5] );

  dWsdX[5][19] = -kf[5] * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * 1.6e0
    + kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 2.0 * 0.2e0;

  for ( k = 0; k < ns; k++ ) {
    dWsdX[5][k] +=
      ( -1.0 * kf[5] * X[4] * X[5] + kb[5] * X[3] + 2.0 * kf[11] * X[11] - 2.0 * kb[11] * X[5] * X[5] );
  }

  dWsdX[6][0] = ( -1.0 * kf[9] * X[6] );

  dWsdX[6][1] =
    ( 1.0 * kf[6] * X[4] * third + 1.0 * kb[8] * X[2] + 1.0 * kb[9] * X[5] + 1.0 * kb[10] * X[3] +
      1.0 * kf[17] * X[12]
      + 1.0 * kf[22] * X[13] + 1.0 * kf[35] * X[10] + 1.0 * kf[40] * X[14] + 1.0 * kf[46] * X[17] );

  dWsdX[6][2] = ( 1.0 * kb[8] * X[1] );

  dWsdX[6][3] = ( 1.0 * kb[10] * X[1] );

  dWsdX[6][4] = ( 1.0 * kf[6] * X[1] * third - 1.0 * kf[7] * X[6] - 1.0 * kf[8] * X[6] );

  dWsdX[6][5] = ( 2.0 * kb[7] * X[5] + 1.0 * kb[9] * X[1] - 1.0 * kf[10] * X[6] + 1.0 * kb[29] * X[14] );

  dWsdX[6][6] =
    ( -1.0 * kb[6] * third - 1.0 * kf[7] * X[4] - 1.0 * kf[8] * X[4] - 1.0 * kf[9] * X[0] -
      1.0 * kf[10] * X[5]
      - 1.0 * kb[17] * X[7] - 1.0 * kb[22] * X[12] - 1.0 * kf[29] * X[9] - 1.0 * kb[35] * X[9] -
      1.0 * kf[37] * X[10]
      - 1.0 * kb[40] * X[13] - 1.0 * kb[46] * X[16] );

  dWsdX[6][7] = ( -1.0 * kb[17] * X[6] );

  dWsdX[6][9] = ( -1.0 * kf[29] * X[6] - 1.0 * kb[35] * X[6] + 1.0 * kb[37] * X[11] );

  dWsdX[6][10] = ( 1.0 * kf[35] * X[1] - 1.0 * kf[37] * X[6] );

  dWsdX[6][11] = ( 1.0 * kb[37] * X[9] );

  dWsdX[6][12] = ( 1.0 * kf[17] * X[1] - 1.0 * kb[22] * X[6] );

  dWsdX[6][13] = ( 1.0 * kf[22] * X[1] - 1.0 * kb[40] * X[6] );

  dWsdX[6][14] = ( 1.0 * kb[29] * X[5] + 1.0 * kf[40] * X[1] );

  dWsdX[6][16] = ( -1.0 * kb[46] * X[6] );

  dWsdX[6][17] = ( 1.0 * kf[46] * X[1] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[6][k] += ( 1.0 * kf[6] * X[4] * X[1] - 1.0 * kb[6] * X[6] );
  }

  dWsdX[7][0] = ( -1.0 * kf[13] * X[7] * third13 + 1.0 * kf[15] * X[12] );

  dWsdX[7][1] = ( 1.0 * kf[17] * X[12] )
    - kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 11.0e0;

  dWsdX[7][2] = ( -1.0 * kb[14] * X[7] )
    + kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * 0.87e0;

  dWsdX[7][3] = ( -1.0 * kb[16] * X[7] )
    + kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * 7.12e0;

  dWsdX[7][4] = ( 1.0 * kb[12] * X[8] + 1.0 * kf[14] * X[12] - 1.0 * kb[18] * X[7] * third18 );

  dWsdX[7][5] = ( -1.0 * kf[12] * X[7] - 1.0 * kb[15] * X[7] + 1.0 * kf[16] * X[12] );

  dWsdX[7][6] = ( -1.0 * kb[17] * X[7] );

  dWsdX[7][7] =
    ( -1.0 * kf[12] * X[5] - 1.0 * kf[13] * X[0] * third13 - 1.0 * kb[14] * X[2] - 1.0 * kb[15] * X[5]
      - 1.0 * kb[16] * X[3] - 1.0 * kb[17] * X[6] - 1.0 * kb[18] * X[4] * third18 - 1.0 * kb[30] * X[10] )
    - kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 2.0e0;

  dWsdX[7][8] = ( 1.0 * kb[12] * X[4] + 1.0 * kb[13] * third13 )
    - kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 6.0e0;

  dWsdX[7][9] = ( 1.0 * kf[30] * X[12] );

  dWsdX[7][10] = ( -1.0 * kb[30] * X[7] );

  dWsdX[7][12] = ( 1.0 * kf[14] * X[4] + 1.0 * kf[15] * X[0] + 1.0 * kf[16] * X[5] + 1.0 * kf[17] * X[1]
                   + 1.0 * kf[18] * third18 + 1.0 * kf[30] * X[9] );
  dWsdX[7][19] = -kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[7][k] +=
      ( -1.0e0 * kf[13] * X[7] * X[0] + 1.0 * kb[13] * X[8] + 1.0 * kf[18] * X[12] -
        1.0 * kb[18] * X[7] * X[4] );
  }

  dWsdX[8][0] = ( 1.0 * kf[13] * X[7] * third13 );

  dWsdX[8][1] = kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 11.0e0;

  dWsdX[8][4] = ( -1.0 * kb[12] * X[8] );

  dWsdX[8][5] = ( 1.0 * kf[12] * X[7] );

  dWsdX[8][7] = ( 1.0 * kf[12] * X[5] + 1.0 * kf[13] * X[0] * third13 )
    + kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 2.0e0;

  dWsdX[8][8] = ( -1.0 * kb[12] * X[4] - 1.0 * kb[13] * third13 )
    + kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * 6.0e0;

  dWsdX[8][19] = +kf[13] * ( X[7] * X[0] - ( X[8] ) / Kc[13] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[8][k] += ( 1.0e0 * kf[13] * X[7] * X[0] - 1.0 * kb[13] * X[8] );
  }

  dWsdX[9][0] = ( -1.0 * kf[25] * X[9] + 1.0 * kb[27] * X[14] + 1.0 * kf[34] * X[10] );

  dWsdX[9][1] = ( -1.0 * kf[27] * X[9] - 1.0 * kf[28] * X[9] + 1.0 * kf[35] * X[10] );

  dWsdX[9][2] = ( 1.0 * kb[26] * X[13] - 1.0 * kb[33] * X[9] );

  dWsdX[9][3] = ( -1.0 * kb[36] * X[9] );

  dWsdX[9][4] = ( 1.0 * kb[25] * X[13] + 2.0 * kb[31] * X[17] - 1.0 * kb[32] * X[9] + 1.0 * kf[33] * X[10] );

  dWsdX[9][5] =
    ( -1.0 * kf[26] * X[9] + 1.0 * kb[28] * X[13] + 1.0 * kb[29] * X[14] - 1.0 * kb[34] * X[9] +
      1.0 * kf[36] * X[10] );

  dWsdX[9][6] = ( -1.0 * kf[29] * X[9] - 1.0 * kb[35] * X[9] + 1.0 * kf[37] * X[10] );

  dWsdX[9][7] = ( 1.0 * kb[30] * X[10] );

  dWsdX[9][9] =
    ( -1.0 * kf[23] * X[13] - 1.0 * kf[25] * X[0] - 1.0 * kf[26] * X[5] - 1.0 * kf[27] * X[1] -
      1.0 * kf[28] * X[1]
      - 1.0 * kf[29] * X[6] - 1.0 * kf[30] * X[12] - 4.0 * kf[31] * X[9] - 1.0 * kb[32] * X[4] -
      1.0 * kb[33] * X[2]
      - 1.0 * kb[34] * X[5] - 1.0 * kb[35] * X[6] - 1.0 * kb[36] * X[3] - 1.0 * kb[37] * X[11] -
      4.0 * kb[48] * X[9]
      - 1.0 * kf[51] * X[18] );

  dWsdX[9][10] =
    ( 1.0 * kb[23] * X[12] + 1.0 * kb[30] * X[7] + 1.0 * kf[32] + 1.0 * kf[33] * X[4] + 1.0 * kf[34] * X[0]
      + 1.0 * kf[35] * X[1] + 1.0 * kf[36] * X[5] + 1.0 * kf[37] * X[6] + 1.0 * kb[51] * X[17] );

  dWsdX[9][11] = ( -1.0 * kb[37] * X[9] );

  dWsdX[9][12] = ( 1.0 * kb[23] * X[10] - 1.0 * kf[30] * X[9] );

  dWsdX[9][13] = ( -1.0 * kf[23] * X[9] + 1.0 * kb[25] * X[4] + 1.0 * kb[26] * X[2] + 1.0 * kb[28] * X[5] );

  dWsdX[9][14] = ( 1.0 * kb[27] * X[0] + 1.0 * kb[29] * X[5] );

  dWsdX[9][17] = ( 2.0 * kb[31] * X[4] + 1.0 * kb[51] * X[10] );

  dWsdX[9][18] = ( 2.0 * kf[48] - 1.0 * kf[51] * X[9] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[9][k] += _dkfdthird ( 32, third, T ) * ( X[10] - ( X[9] * X[4] ) / Kc[32] )
      + _dkfdthird ( 48, third, T ) * ( X[18] - ( X[9] * X[9] ) / Kc[48] ) * 2.0;
  }

  dWsdX[10][0] = ( -1.0 * kf[34] * X[10] );

  dWsdX[10][1] = ( -1.0 * kf[35] * X[10] );

  dWsdX[10][2] = ( 1.0 * kb[33] * X[9] );

  dWsdX[10][3] = ( 1.0 * kb[36] * X[9] );

  dWsdX[10][4] = ( 1.0 * kb[32] * X[9] - 1.0 * kf[33] * X[10] );

  dWsdX[10][5] = ( 1.0 * kb[34] * X[9] - 1.0 * kf[36] * X[10] );

  dWsdX[10][6] = ( 1.0 * kb[35] * X[9] - 1.0 * kf[37] * X[10] );

  dWsdX[10][7] = ( -1.0 * kb[30] * X[10] );

  dWsdX[10][9] =
    ( 1.0 * kf[23] * X[13] + 1.0 * kf[30] * X[12] + 1.0 * kb[32] * X[4] + 1.0 * kb[33] * X[2] +
      1.0 * kb[34] * X[5]
      + 1.0 * kb[35] * X[6] + 1.0 * kb[36] * X[3] + 1.0 * kb[37] * X[11] + 1.0 * kf[51] * X[18] );

  dWsdX[10][10] =
    ( -1.0 * kb[23] * X[12] - 1.0 * kb[30] * X[7] - 1.0 * kf[32] - 1.0 * kf[33] * X[4] - 1.0 * kf[34] * X[0]
      - 1.0 * kf[35] * X[1] - 1.0 * kf[36] * X[5] - 1.0 * kf[37] * X[6] - 1.0 * kb[51] * X[17] );

  dWsdX[10][11] = ( 1.0 * kb[37] * X[9] );

  dWsdX[10][12] = ( -1.0 * kb[23] * X[10] + 1.0 * kf[30] * X[9] );

  dWsdX[10][13] = ( 1.0 * kf[23] * X[9] );

  dWsdX[10][17] = ( -1.0 * kb[51] * X[10] );

  dWsdX[10][18] = ( 1.0 * kf[51] * X[9] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[10][k] += -_dkfdthird ( 32, third, T ) * ( X[10] - ( X[9] * X[4] ) / Kc[32] );
  }

  dWsdX[11][1] = -kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 0.2e0;

  dWsdX[11][2] = -kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 1.9e0;

  dWsdX[11][3] = -kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 17.5e0;

  dWsdX[11][5] = ( 2.0 * kb[11] * X[5] * third11 );

  dWsdX[11][6] = ( 1.0 * kf[37] * X[10] );

  dWsdX[11][7] = -kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 1.1e0;

  dWsdX[11][8] = -kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 3.3e0;

  dWsdX[11][9] = ( -1.0 * kb[37] * X[11] );

  dWsdX[11][10] = ( 1.0 * kf[37] * X[6] );

  dWsdX[11][11] = ( -1.0 * kf[11] * third11 - 1.0 * kb[37] * X[9] );

  dWsdX[11][19] = -kf[11] * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * 0.2e0;

  for ( k = 0; k < ns; k++ ) {
    dWsdX[11][k] += ( -1.0 * kf[11] * X[11] + 1.0 * kb[11] * X[5] * X[5] );
  }

  dWsdX[12][0] = ( -1.0 * kf[15] * X[12] + 1.0 * kf[20] * X[13] );

  dWsdX[12][1] = ( -1.0 * kf[17] * X[12] + 1.0 * kf[22] * X[13] + 1.0 * kf[42] * X[15] )
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 0.2e0;

  dWsdX[12][2] = ( 1.0 * kb[14] * X[7] - 1.0 * kb[19] * X[12] )
    - kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * 0.87e0
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 1.9e0;

  dWsdX[12][3] = ( 1.0 * kb[16] * X[7] - 1.0 * kb[21] * X[12] )
    - kf[18] * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * 7.12e0
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 17.5e0;

  dWsdX[12][4] =
    ( -1.0 * kf[14] * X[12] + 1.0 * kb[18] * X[7] * third18 + 1.0 * kf[19] * X[13] -
      1.0 * kb[24] * X[12] * third24 );

  dWsdX[12][5] = ( 1.0 * kb[15] * X[7] - 1.0 * kf[16] * X[12] - 1.0 * kb[20] * X[12] + 1.0 * kf[21] * X[13] );

  dWsdX[12][6] = ( 1.0 * kb[17] * X[7] - 1.0 * kb[22] * X[12] );

  dWsdX[12][7] = ( 1.0 * kb[14] * X[2] + 1.0 * kb[15] * X[5] + 1.0 * kb[16] * X[3] + 1.0 * kb[17] * X[6]
                   + 1.0 * kb[18] * X[4] * third18 + 1.0 * kb[30] * X[10] )
    + kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 1.1;

  dWsdX[12][8] = +kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 3.3;

  dWsdX[12][9] = ( 1.0 * kf[23] * X[13] - 1.0 * kf[30] * X[12] );

  dWsdX[12][10] = ( -1.0 * kb[23] * X[12] + 1.0 * kb[30] * X[7] );

  dWsdX[12][12] =
    ( -1.0 * kf[14] * X[4] - 1.0 * kf[15] * X[0] - 1.0 * kf[16] * X[5] - 1.0 * kf[17] * X[1] -
      1.0 * kf[18] * third18 - 1.0 * kb[19] * X[2] - 1.0 * kb[20] * X[5] - 1.0 * kb[21] * X[3] -
      1.0 * kb[22] * X[6] - 1.0 * kb[23] * X[10]
      - 1.0 * kb[24] * X[4] * third24 - 1.0 * kf[30] * X[9] - 1.0 * kb[42] * X[13] );

  dWsdX[12][13] =
    ( 1.0 * kf[19] * X[4] + 1.0 * kf[20] * X[0] + 1.0 * kf[21] * X[5] + 1.0 * kf[22] * X[1] +
      1.0 * kf[23] * X[9]
      + 1.0 * kf[24] * third24 - 1.0 * kb[42] * X[12] );

  dWsdX[12][15] = ( 1.0 * kf[42] * X[1] );

  dWsdX[12][19] = +kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 0.2e0;

  for ( k = 0; k < ns; k++ ) {
    dWsdX[12][k] +=
      ( -1.0 * kf[18] * X[12] + 1.0 * kb[18] * X[7] * X[4] + 1.0 * kf[24] * X[13] -
        1.0 * kb[24] * X[12] * X[4] );
  }

  dWsdX[13][0] = ( -1.0 * kf[20] * X[13] + 1.0 * kf[25] * X[9] );

  dWsdX[13][1] = ( -1.0 * kf[22] * X[13] + 1.0 * kf[28] * X[9] + 1.0 * kf[40] * X[14] + 1.0 * kf[42] * X[15] )
    - kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 0.2e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 0.2e0;

  dWsdX[13][2] = ( 1.0 * kb[19] * X[12] - 1.0 * kb[26] * X[13] - 1.0 * kb[38] * X[13] )
    - kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 1.9e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 1.9e0;

  dWsdX[13][3] = ( 1.0 * kb[21] * X[12] - 1.0 * kb[39] * X[13] )
    - kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 17.5e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 17.5e0;

  dWsdX[13][4] =
    ( -1.0 * kf[19] * X[13] + 1.0 * kb[24] * X[12] * third24 - 1.0 * kb[25] * X[13] + 1.0 * kf[38] * X[14]
      - 1.0 * kb[41] * X[13] * third41 );

  dWsdX[13][5] = ( 1.0 * kb[20] * X[12] - 1.0 * kf[21] * X[13] + 1.0 * kf[26] * X[9] - 1.0 * kb[28] * X[13]
                   + 1.0 * kf[39] * X[14] );

  dWsdX[13][6] = ( 1.0 * kb[22] * X[12] - 1.0 * kb[40] * X[13] );

  dWsdX[13][7] = -kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 1.1e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 1.1e0;

  dWsdX[13][8] = -kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 3.3e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 3.3e0;

  dWsdX[13][9] = ( -1.0 * kf[23] * X[13] + 1.0 * kf[25] * X[0] + 1.0 * kf[26] * X[5] + 1.0 * kf[28] * X[1] );

  dWsdX[13][10] = ( 1.0 * kb[23] * X[12] );

  dWsdX[13][12] =
    ( 1.0 * kb[19] * X[2] + 1.0 * kb[20] * X[5] + 1.0 * kb[21] * X[3] + 1.0 * kb[22] * X[6] +
      1.0 * kb[23] * X[10]
      + 1.0 * kb[24] * X[4] * third24 - 1.0 * kb[42] * X[13] );

  dWsdX[13][13] =
    ( -1.0 * kf[19] * X[4] - 1.0 * kf[20] * X[0] - 1.0 * kf[21] * X[5] - 1.0 * kf[22] * X[1] -
      1.0 * kf[23] * X[9]
      - 1.0 * kf[24] * third24 - 1.0 * kb[25] * X[4] - 1.0 * kb[26] * X[2] - 1.0 * kb[28] * X[5] -
      1.0 * kb[38] * X[2]
      - 1.0 * kb[39] * X[3] - 1.0 * kb[40] * X[6] - 1.0 * kb[41] * X[4] * third41 - 1.0 * kb[42] * X[12] );

  dWsdX[13][14] =
    ( 1.0 * kf[38] * X[4] + 1.0 * kf[39] * X[5] + 1.0 * kf[40] * X[1] + 1.0 * kf[41] * third41 );

  dWsdX[13][15] = ( 1.0 * kf[42] * X[1] );

  dWsdX[13][19] = -kf[24] * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * 0.2e0
    + kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 0.2e0;

  for ( k = 0; k < ns; k++ ) {
    dWsdX[13][k] +=
      ( -1.0 * kf[24] * X[13] + 1.0 * kb[24] * X[12] * X[4] + 1.0 * kf[41] * X[14] - kb[41] * X[13] * X[4] );
  }

  dWsdX[14][0] = ( -1.0 * kb[27] * X[14] );

  dWsdX[14][1] = ( 1.0 * kf[27] * X[9] - 1.0 * kf[40] * X[14] )
    - kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 0.2e0;

  dWsdX[14][2] = ( 1.0 * kb[38] * X[13] )
    - kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 1.9e0;

  dWsdX[14][3] = ( 1.0 * kb[39] * X[13] )
    - kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 17.5e0;

  dWsdX[14][4] = ( -1.0 * kf[38] * X[14] + 1.0 * kb[41] * X[13] * third41 );

  dWsdX[14][5] = ( -1.0 * kb[29] * X[14] - 1.0 * kf[39] * X[14] );

  dWsdX[14][6] = ( 1.0 * kf[29] * X[9] + 1.0 * kb[40] * X[13] );

  dWsdX[14][7] = -kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 1.1e0;

  dWsdX[14][8] = -kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 3.3e0;

  dWsdX[14][9] = ( 1.0 * kf[27] * X[1] + 1.0 * kf[29] * X[6] );

  dWsdX[14][13] =
    ( 1.0 * kb[38] * X[2] + 1.0 * kb[39] * X[3] + 1.0 * kb[40] * X[6] + 1.0 * kb[41] * X[4] * third41 );

  dWsdX[14][14] =
    ( -1.0 * kb[27] * X[0] - 1.0 * kb[29] * X[5] - 1.0 * kf[38] * X[4] - 1.0 * kf[39] * X[5] -
      1.0 * kf[40] * X[1]
      - 1.0 * kf[41] * third41 );

  dWsdX[14][19] = -kf[41] * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * 0.2e0;

  for ( k = 0; k < ns; k++ ) {
    dWsdX[14][k] += ( -1.0 * kf[41] * X[14] + 1.0 * kb[41] * X[13] * X[4] );
  }

  dWsdX[15][1] = ( -1.0 * kf[42] * X[15] );

  dWsdX[15][2] = ( -1.0 * kb[43] * X[15] );

  dWsdX[15][3] = ( -1.0 * kb[44] * X[15] );

  dWsdX[15][4] = ( 1.0 * kf[43] * X[16] );

  dWsdX[15][5] = ( 1.0 * kf[44] * X[16] );

  dWsdX[15][12] = ( 1.0 * kb[42] * X[13] );

  dWsdX[15][13] = ( 1.0 * kb[42] * X[12] );

  dWsdX[15][15] = ( -1.0 * kf[42] * X[1] - 1.0 * kb[43] * X[2] - 1.0 * kb[44] * X[3] );

  dWsdX[15][16] = ( 1.0 * kf[43] * X[4] + 1.0 * kf[44] * X[5] );

  dWsdX[16][1] = ( 1.0 * kf[46] * X[17] );

  dWsdX[16][2] = ( 1.0 * kb[43] * X[15] - 1.0 * kb[45] * X[16] );

  dWsdX[16][3] = ( 1.0 * kb[44] * X[15] );

  dWsdX[16][4] = ( -1.0 * kf[43] * X[16] + 1.0 * kf[45] * X[17] - 1.0 * kb[47] * X[16] );

  dWsdX[16][5] = ( -1.0 * kf[44] * X[16] );

  dWsdX[16][6] = ( -1.0 * kb[46] * X[16] );

  dWsdX[16][15] = ( 1.0 * kb[43] * X[2] + 1.0 * kb[44] * X[3] );

  dWsdX[16][16] =
    ( -1.0 * kf[43] * X[4] - 1.0 * kf[44] * X[5] - 1.0 * kb[45] * X[2] - 1.0 * kb[46] * X[6] -
      1.0 * kb[47] * X[4] );

  dWsdX[16][17] = ( 1.0 * kf[45] * X[4] + 1.0 * kf[46] * X[1] + 1.0 * kf[47] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[16][k] += ( X[17] - ( X[16] * X[4] ) / Kc[47] ) * _dkfdthird ( 47, third, T );
  }

  dWsdX[17][1] = ( -1.0 * kf[46] * X[17] );

  dWsdX[17][2] = ( 1.0 * kb[45] * X[16] - 1.0 * kb[49] * X[17] );

  dWsdX[17][3] = ( -1.0 * kb[50] * X[17] );

  dWsdX[17][4] =
    ( -1.0 * kb[31] * X[17] - 1.0 * kf[45] * X[17] + 1.0 * kb[47] * X[16] + 1.0 * kf[49] * X[18] );

  dWsdX[17][5] = ( 1.0 * kf[50] * X[18] );

  dWsdX[17][6] = ( 1.0 * kb[46] * X[16] );

  dWsdX[17][9] = ( 2.0 * kf[31] * X[9] + 1.0 * kf[51] * X[18] );

  dWsdX[17][10] = ( -1.0 * kb[51] * X[17] );

  dWsdX[17][16] = ( 1.0 * kb[45] * X[2] + 1.0 * kb[46] * X[6] + 1.0 * kb[47] * X[4] );

  dWsdX[17][17] =
    ( -1.0 * kb[31] * X[4] - 1.0 * kf[45] * X[4] - 1.0 * kf[46] * X[1] - 1.0 * kf[47] - 1.0 * kb[49] * X[2]
      - 1.0 * kb[50] * X[3] - 1.0 * kb[51] * X[10] );

  dWsdX[17][18] = ( 1.0 * kf[49] * X[4] + 1.0 * kf[50] * X[5] + 1.0 * kf[51] * X[9] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[17][k] += -( X[17] - ( X[16] * X[4] ) / Kc[47] ) * _dkfdthird ( 47, third, T );
  }

  dWsdX[18][2] = ( 1.0 * kb[49] * X[17] );

  dWsdX[18][3] = ( 1.0 * kb[50] * X[17] );

  dWsdX[18][4] = ( -1.0 * kf[49] * X[18] );

  dWsdX[18][5] = ( -1.0 * kf[50] * X[18] );

  dWsdX[18][9] = ( 2.0 * kb[48] * X[9] - 1.0 * kf[51] * X[18] );

  dWsdX[18][10] = ( 1.0 * kb[51] * X[17] );

  dWsdX[18][17] = ( 1.0 * kb[49] * X[2] + 1.0 * kb[50] * X[3] + 1.0 * kb[51] * X[10] );

  dWsdX[18][18] = ( -1.0 * kf[48] - 1.0 * kf[49] * X[4] - 1.0 * kf[50] * X[5] - 1.0 * kf[51] * X[9] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[18][k] += -( X[18] - ( X[9] * X[9] ) / Kc[48] ) * _dkfdthird ( 48, third, T );
  }

  for ( k = 0; k < ns; k++ ) {
    for ( p = 0; p < ns; p++ ) {
      dWdrhok[k][p] = dWsdX[k][p] / _calM ( p ) * _calM ( k );
    }
  }

}




void find_Qei(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){  
  *Qei=0.0;
}


void find_dQei_dx(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
}

