#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/.active/model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
#include <src/common.h>         /* only needed for display_matrix() function */

#define nr 52                   /* Note: ns is specified in chem.hh */

/*
chemSim.w is the simplified version of chem.w
this file calculates _dkfdT of reaction 33, 48, and 49 from kf_inf only
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

static double n_f[nr + 3] = {
  -0.927E+0, 2.700E+0, 1.510E+0, 1.400E+0, -1.000E+0, -2.000E+0, -0.800E+0, 0.000E+0, 0.000E+0, 0.000E+0,
  0.000E+0, 0.000E+0, 1.350E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, -1.000E+0, 1.620E+0,
  0.000E+0, 2.460E+0, 0.000E+0, 7.400E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.000E+0,
  0.500E+0, 0.000E+0, -0.558E+0, 2.110E+0, 1.440E+0, 0.000E+0, 2.130E+0, 0.000E+0, 0.000E+0, 0.000E+0,
  7.600E+0, 0.000E+0, 0.000E+0, 0.700E+0, 0.000E+0, 0.000E+0, 0.000E+0, 0.730E+0, -2.792E+0, 3.500E+0,
  1.900E+0, 4.000E+0, -4.911E+0, -6.800E+0, -11.992E+0
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

void write_chem_template ( FILE ** controlfile ) {

}

void read_and_init_chem_actions ( char *action, char **argum, SOAP_codex_t * codex ) {

}

/* species concentration [mol/cm^3] */
double _Xk ( long k, np_t np ) {
  double tmp;
  tmp = _w ( np, k ) * _rho ( np ) / SpeciesMolWeight ( k ) * 1.0e-6;
  return ( tmp );
}

/* gibbs energy of species [J/mol] */
double _Gk ( long k, np_t np ) {
  double tmp, T;
  T = _T ( np );
  tmp = ( _hk_from_T ( k, T ) - T * _sk_from_T ( k, T ) ) * SpeciesMolWeight ( k );
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

double _kf ( long r, double third, np_t np ) {
  double tmp, T, kf_inf, kf_o, Fc, pr, xt;
  T = _T ( np );
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

/* derivative of gibbs energy of species k wrt temperature [J/(mol K)] */
double _dGkdT ( long k, np_t np ) {
  double tmp, T;
  T = _T ( np );
  tmp = ( _cpk_from_T ( k, T ) - _dTskdT_from_T ( k, T ) ) * SpeciesMolWeight ( k );
  return ( tmp );
}

/* 
derivative of the forward reaction rate wrt the density of species k 

Notes: for reaction 33, 48, and 49, take derivative of Kf_inf only (simplification)
*/

double _dkfdT ( long r, double third, np_t np ) {
  double tmp, T;
  T = _T ( np );
  tmp = _kf ( r, third, np ) * ( n_f[r] / T + E[r] / ( T * T ) );
  return ( tmp );
}

/* 
derivative of the concentration input in [mol/cm^3] of species k 
wrt the density [kg/m^3] of species k 
*/
double _dXdrhok ( long k ) {
  double tmp;
  tmp = 1.0e-6 / SpeciesMolWeight ( k );
  return ( tmp );
}

void find_Schem ( np_t np, gl_t * gl, flux_t S ) {
  long flux;
  long r, k;
  double T;
  double X[ns], Gs[ns];
  double dGs[nr], Kc[nr], kf[nr], react[nr];
  double third, third5, third11, third13, third18, third24, third41;

  T = _T ( np );
  for ( flux = 0; flux < nf; flux++ ) {
    S[flux] = 0.0E+00;
  }

  third = 0.0E+0;
  for ( k = 0; k < ns; k++ ) {
    X[k] = _Xk ( k, np );
    Gs[k] = _Gk ( k, np );
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
    kf[r] = _kf ( r, third, np );
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
  S[0] =
    ( react[0] - react[1] + react[3] - react[9] - react[13] - react[15] - react[20] - react[25] + react[27]
      - react[34] )
    * 1.0E+6 * SpeciesMolWeight ( 0 );
  S[1] =
    ( -react[0] - react[6] + react[8] + react[9] + react[10] - react[17] - react[22] - react[27] - react[28]
      - react[35] - react[40] - react[42] - react[46] )
    * 1.0E+6 * SpeciesMolWeight ( 1 );
  S[2] =
    ( -react[1] - react[2] + react[4] + react[8] + react[14] + react[19] + react[26] + react[33] + react[38]
      + react[43] + react[45] + react[49] )
    * 1.0E+6 * SpeciesMolWeight ( 2 );
  S[3] =
    ( react[2] + react[3] + react[5] + react[10] + react[16] + react[21] + react[36] + react[39] + react[44]
      + react[50] )
    * 1.0E+6 * SpeciesMolWeight ( 3 );
  S[4] = ( -react[0] + react[1] + react[2] - react[4] - react[4] - react[5] - react[6] - react[7] - react[8]
           + react[12] - react[14] + react[18] - react[19] + react[24] + react[25] + react[31] + react[32] -
           react[33]
           - react[38] + react[41] - react[43] - react[45] + react[47] - react[49] )
    * 1.0E+6 * SpeciesMolWeight ( 4 );
  S[5] = ( react[0] + react[1] - react[2] - react[3] - react[3] - react[5] + react[7] + react[7] + react[9]
           - react[10] + react[11] + react[11] - react[12] + react[15] - react[16] + react[20] - react[21] -
           react[26]
           + react[28] + react[29] + react[34] - react[36] - react[39] - react[44] - react[50] )
    * 1.0E+6 * SpeciesMolWeight ( 5 );
  S[6] =
    ( react[6] - react[7] - react[8] - react[9] - react[10] + react[17] + react[22] - react[29] + react[35]
      - react[37] + react[40] + react[46] )
    * 1.0E+6 * SpeciesMolWeight ( 6 );
  S[7] = ( -react[12] - react[13] + react[14] + react[15] + react[16] + react[17] + react[18] + react[30] )
    * 1.0E+6 * SpeciesMolWeight ( 7 );
  S[8] = ( react[12] + react[13] )
    * 1.0E+6 * SpeciesMolWeight ( 8 );
  S[9] =
    ( -react[23] - react[25] - react[26] - react[27] - react[28] - react[29] - react[30] - react[31] -
      react[31]
      + react[32] + react[33] + react[34] + react[35] + react[36] + react[37] + react[48] + react[48] -
      react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 9 );
  S[10] =
    ( react[23] + react[30] - react[32] - react[33] - react[34] - react[35] - react[36] - react[37] +
      react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 10 );
  S[11] = ( -react[11] + react[37] )
    * 1.0E+6 * SpeciesMolWeight ( 11 );
  S[12] =
    ( -react[14] - react[15] - react[16] - react[17] - react[18] + react[19] + react[20] + react[21] +
      react[22]
      + react[23] + react[24] - react[30] + react[42] )
    * 1.0E+6 * SpeciesMolWeight ( 12 );
  S[13] =
    ( -react[19] - react[20] - react[21] - react[22] - react[23] - react[24] + react[25] + react[26] +
      react[28]
      + react[38] + react[39] + react[40] + react[41] + react[42] )
    * 1.0E+6 * SpeciesMolWeight ( 13 );
  S[14] = ( react[27] + react[29] - react[38] - react[39] - react[40] - react[41] )
    * 1.0E+6 * SpeciesMolWeight ( 14 );
  S[15] = ( -react[42] + react[43] + react[44] )
    * 1.0E+6 * SpeciesMolWeight ( 15 );
  S[16] = ( -react[43] - react[44] + react[45] + react[46] + react[47] )
    * 1.0E+6 * SpeciesMolWeight ( 16 );
  S[17] = ( react[31] - react[45] - react[46] - react[47] + react[49] + react[50] + react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 17 );
  S[18] = ( -react[48] - react[49] - react[50] - react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 18 );

  for ( k = 0; k < ns; k++ ) {
    S[k] = S[k] * _Omega ( np, gl );
  }

}

/*end of find_Schem */

void find_dSchem_dU ( np_t np, gl_t * gl, sqmat_t C ) {
  long flux, flux2;
  long r, k, p;
  double T;
  double Kc[nr], kf[nr], kb[nr], dGs[nr], react[nr], dKcdT[nr];
  double X[ns], Gs[ns], dGsdT[ns], dWsdT[ns];
  double dWsdX[ns][ns];
  double dTdrhoE;
  double third, third5, third11, third13, third18, third24, third41;
  spec_t dTdrhok;
  dim_t dTdrhoV;

  T = _T ( np );
  for ( flux = 0; flux < nf; flux++ ) {
    for ( flux2 = 0; flux2 < nf; flux2++ ) {
      C[flux][flux2] = 0.0e0;
    }
  }

  third = 0;
  for ( k = 0; k < ns; k++ ) {
    X[k] = _Xk ( k, np );
    Gs[k] = _Gk ( k, np );
    dGsdT[k] = _dGkdT ( k, np );
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

/* Calculate the difference between the gibbs free energy between the products 
   and reactants of each reaction */

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
determine Kc by multiplying Kp by 
((RuT/100000*100**3)**(prod-react))
[Kc] = (cm**3/mol) ** (prod-react) 

Note: WARP version prior Feb09 2006 has the form
((RuT/101325*100**3)**(prod-react))
the change from 101325 to 100000 is due to the 
implementation of NASA new thermo-coefficient data set
*/

  for ( r = 0; r < nr; r++ ) {
    Kc[r] = exp ( -dGs[r] / ( calR * T ) );
  }

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

  for ( r = 0; r < nr; r++ ) {
    kf[r] = _kf ( r, third, np );
    kb[r] = kf[r] / Kc[r];
  }

/* find the derivative of each source term wrt temperature */

  react[0] = _dkfdT ( 0, third, np ) * ( X[4] * X[1] - ( X[5] * X[0] ) / Kc[0] );
  react[1] = _dkfdT ( 1, third, np ) * ( X[0] * X[2] - ( X[5] * X[4] ) / Kc[1] );
  react[2] = _dkfdT ( 2, third, np ) * ( X[5] * X[2] - ( X[3] * X[4] ) / Kc[2] );
  react[3] = _dkfdT ( 3, third, np ) * ( X[5] * X[5] - ( X[0] * X[3] ) / Kc[3] );
  react[4] = _dkfdT ( 4, third, np ) * ( X[4] * X[4] - ( X[2] ) / Kc[4] ) * third;
  react[5] = _dkfdT ( 5, third, np ) * ( X[4] * X[5] - ( X[3] ) / Kc[5] ) * third5;
  react[6] = _dkfdT ( 6, third, np ) * ( X[4] * X[1] - ( X[6] ) / Kc[6] ) * third;
  react[7] = _dkfdT ( 7, third, np ) * ( X[6] * X[4] - ( X[5] * X[5] ) / Kc[7] );
  react[8] = _dkfdT ( 8, third, np ) * ( X[6] * X[4] - ( X[2] * X[1] ) / Kc[8] );
  react[9] = _dkfdT ( 9, third, np ) * ( X[6] * X[0] - ( X[1] * X[5] ) / Kc[9] );
  react[10] = _dkfdT ( 10, third, np ) * ( X[6] * X[5] - ( X[3] * X[1] ) / Kc[10] );
  react[11] = _dkfdT ( 11, third, np ) * ( X[11] - ( X[5] * X[5] ) / Kc[11] ) * third11;
  react[12] = _dkfdT ( 12, third, np ) * ( X[7] * X[5] - ( X[8] * X[4] ) / Kc[12] );
  react[13] = _dkfdT ( 13, third, np ) * ( X[7] * X[0] - ( X[8] ) / Kc[13] ) * third13;
  react[14] = _dkfdT ( 14, third, np ) * ( X[12] * X[4] - ( X[7] * X[2] ) / Kc[14] );
  react[15] = _dkfdT ( 15, third, np ) * ( X[12] * X[0] - ( X[7] * X[5] ) / Kc[15] );
  react[16] = _dkfdT ( 16, third, np ) * ( X[12] * X[5] - ( X[7] * X[3] ) / Kc[16] );
  react[17] = _dkfdT ( 17, third, np ) * ( X[12] * X[1] - ( X[7] * X[6] ) / Kc[17] );
  react[18] = _dkfdT ( 18, third, np ) * ( X[12] - ( X[7] * X[4] ) / Kc[18] ) * third18;
  react[19] = _dkfdT ( 19, third, np ) * ( X[13] * X[4] - ( X[12] * X[2] ) / Kc[19] );
  react[20] = _dkfdT ( 20, third, np ) * ( X[13] * X[0] - ( X[12] * X[5] ) / Kc[20] );
  react[21] = _dkfdT ( 21, third, np ) * ( X[13] * X[5] - ( X[12] * X[3] ) / Kc[21] );
  react[22] = _dkfdT ( 22, third, np ) * ( X[13] * X[1] - ( X[12] * X[6] ) / Kc[22] );
  react[23] = _dkfdT ( 23, third, np ) * ( X[13] * X[9] - ( X[12] * X[10] ) / Kc[23] );
  react[24] = _dkfdT ( 24, third, np ) * ( X[13] - ( X[12] * X[4] ) / Kc[24] ) * third24;
  react[25] = _dkfdT ( 25, third, np ) * ( X[9] * X[0] - ( X[13] * X[4] ) / Kc[25] );
  react[26] = _dkfdT ( 26, third, np ) * ( X[9] * X[5] - ( X[13] * X[2] ) / Kc[26] );
  react[27] = _dkfdT ( 27, third, np ) * ( X[9] * X[1] - ( X[14] * X[0] ) / Kc[27] );
  react[28] = _dkfdT ( 28, third, np ) * ( X[9] * X[1] - ( X[13] * X[5] ) / Kc[28] );
  react[29] = _dkfdT ( 29, third, np ) * ( X[9] * X[6] - ( X[14] * X[5] ) / Kc[29] );
  react[30] = _dkfdT ( 30, third, np ) * ( X[9] * X[12] - ( X[10] * X[7] ) / Kc[30] );
  react[31] = _dkfdT ( 31, third, np ) * ( X[9] * X[9] - ( X[17] * X[4] ) / Kc[31] );
  react[32] = _dkfdT ( 32, third, np ) * ( X[10] - ( X[9] * X[4] ) / Kc[32] );
  react[33] = _dkfdT ( 33, third, np ) * ( X[10] * X[4] - ( X[9] * X[2] ) / Kc[33] );
  react[34] = _dkfdT ( 34, third, np ) * ( X[10] * X[0] - ( X[9] * X[5] ) / Kc[34] );
  react[35] = _dkfdT ( 35, third, np ) * ( X[10] * X[1] - ( X[9] * X[6] ) / Kc[35] );
  react[36] = _dkfdT ( 36, third, np ) * ( X[10] * X[5] - ( X[9] * X[3] ) / Kc[36] );
  react[37] = _dkfdT ( 37, third, np ) * ( X[10] * X[6] - ( X[9] * X[11] ) / Kc[37] );
  react[38] = _dkfdT ( 38, third, np ) * ( X[14] * X[4] - ( X[13] * X[2] ) / Kc[38] );
  react[39] = _dkfdT ( 39, third, np ) * ( X[14] * X[5] - ( X[13] * X[3] ) / Kc[39] );
  react[40] = _dkfdT ( 40, third, np ) * ( X[14] * X[1] - ( X[13] * X[6] ) / Kc[40] );
  react[41] = _dkfdT ( 41, third, np ) * ( X[14] - ( X[13] * X[4] ) / Kc[41] ) * third41;
  react[42] = _dkfdT ( 42, third, np ) * ( X[15] * X[1] - ( X[13] * X[12] ) / Kc[42] );
  react[43] = _dkfdT ( 43, third, np ) * ( X[16] * X[4] - ( X[15] * X[2] ) / Kc[43] );
  react[44] = _dkfdT ( 44, third, np ) * ( X[16] * X[5] - ( X[15] * X[3] ) / Kc[44] );
  react[45] = _dkfdT ( 45, third, np ) * ( X[17] * X[4] - ( X[16] * X[2] ) / Kc[45] );
  react[46] = _dkfdT ( 46, third, np ) * ( X[17] * X[1] - ( X[16] * X[6] ) / Kc[46] );
  react[47] = _dkfdT ( 47, third, np ) * ( X[17] - ( X[16] * X[4] ) / Kc[47] );
  react[48] = _dkfdT ( 48, third, np ) * ( X[18] - ( X[9] * X[9] ) / Kc[48] );
  react[49] = _dkfdT ( 49, third, np ) * ( X[18] * X[4] - ( X[17] * X[2] ) / Kc[49] );
  react[50] = _dkfdT ( 50, third, np ) * ( X[18] * X[5] - ( X[17] * X[3] ) / Kc[50] );
  react[51] = _dkfdT ( 51, third, np ) * ( X[18] * X[9] - ( X[17] * X[10] ) / Kc[51] );

  dWsdT[0] =
    ( react[0] - react[1] + react[3] - react[9] - react[13] - react[15] - react[20] - react[25] + react[27]
      - react[34] )
    * 1.0E+6 * SpeciesMolWeight ( 0 );
  dWsdT[1] =
    ( -react[0] - react[6] + react[8] + react[9] + react[10] - react[17] - react[22] - react[27] - react[28]
      - react[35] - react[40] - react[42] - react[46] )
    * 1.0E+6 * SpeciesMolWeight ( 1 );
  dWsdT[2] =
    ( -react[1] - react[2] + react[4] + react[8] + react[14] + react[19] + react[26] + react[33] + react[38]
      + react[43] + react[45] + react[49] )
    * 1.0E+6 * SpeciesMolWeight ( 2 );
  dWsdT[3] =
    ( react[2] + react[3] + react[5] + react[10] + react[16] + react[21] + react[36] + react[39] + react[44]
      + react[50] )
    * 1.0E+6 * SpeciesMolWeight ( 3 );
  dWsdT[4] =
    ( -react[0] + react[1] + react[2] - react[4] - react[4] - react[5] - react[6] - react[7] - react[8]
      + react[12] - react[14] + react[18] - react[19] + react[24] + react[25] + react[31] + react[32] -
      react[33]
      - react[38] + react[41] - react[43] - react[45] + react[47] - react[49] )
    * 1.0E+6 * SpeciesMolWeight ( 4 );
  dWsdT[5] =
    ( react[0] + react[1] - react[2] - react[3] - react[3] - react[5] + react[7] + react[7] + react[9]
      - react[10] + react[11] + react[11] - react[12] + react[15] - react[16] + react[20] - react[21] -
      react[26]
      + react[28] + react[29] + react[34] - react[36] - react[39] - react[44] - react[50] )
    * 1.0E+6 * SpeciesMolWeight ( 5 );
  dWsdT[6] =
    ( react[6] - react[7] - react[8] - react[9] - react[10] + react[17] + react[22] - react[29] + react[35]
      - react[37] + react[40] + react[46] )
    * 1.0E+6 * SpeciesMolWeight ( 6 );
  dWsdT[7] =
    ( -react[12] - react[13] + react[14] + react[15] + react[16] + react[17] + react[18] + react[30] )
    * 1.0E+6 * SpeciesMolWeight ( 7 );
  dWsdT[8] = ( react[12] + react[13] )
    * 1.0E+6 * SpeciesMolWeight ( 8 );
  dWsdT[9] =
    ( -react[23] - react[25] - react[26] - react[27] - react[28] - react[29] - react[30] - react[31] -
      react[31]
      + react[32] + react[33] + react[34] + react[35] + react[36] + react[37] + react[48] + react[48] -
      react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 9 );
  dWsdT[10] =
    ( react[23] + react[30] - react[32] - react[33] - react[34] - react[35] - react[36] - react[37] +
      react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 10 );
  dWsdT[11] = ( -react[11] + react[37] )
    * 1.0E+6 * SpeciesMolWeight ( 11 );
  dWsdT[12] =
    ( -react[14] - react[15] - react[16] - react[17] - react[18] + react[19] + react[20] + react[21] +
      react[22]
      + react[23] + react[24] - react[30] + react[42] )
    * 1.0E+6 * SpeciesMolWeight ( 12 );
  dWsdT[13] =
    ( -react[19] - react[20] - react[21] - react[22] - react[23] - react[24] + react[25] + react[26] +
      react[28]
      + react[38] + react[39] + react[40] + react[41] + react[42] )
    * 1.0E+6 * SpeciesMolWeight ( 13 );
  dWsdT[14] = ( react[27] + react[29] - react[38] - react[39] - react[40] - react[41] )
    * 1.0E+6 * SpeciesMolWeight ( 14 );
  dWsdT[15] = ( -react[42] + react[43] + react[44] )
    * 1.0E+6 * SpeciesMolWeight ( 15 );
  dWsdT[16] = ( -react[43] - react[44] + react[45] + react[46] + react[47] )
    * 1.0E+6 * SpeciesMolWeight ( 16 );
  dWsdT[17] = ( react[31] - react[45] - react[46] - react[47] + react[49] + react[50] + react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 17 );
  dWsdT[18] = ( -react[48] - react[49] - react[50] - react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 18 );
  dWsdT[19] = 0.0E+0;

/*******************************************/

  dKcdT[0] = dGs[0] / ( calR * T * T ) - ( dGsdT[5] + dGsdT[0] - ( dGsdT[4] + dGsdT[1] ) ) / calR / T;
  dKcdT[1] = dGs[1] / ( calR * T * T ) - ( dGsdT[5] + dGsdT[4] - ( dGsdT[0] + dGsdT[2] ) ) / calR / T;
  dKcdT[2] = dGs[2] / ( calR * T * T ) - ( dGsdT[3] + dGsdT[4] - ( dGsdT[5] + dGsdT[2] ) ) / calR / T;
  dKcdT[3] = dGs[3] / ( calR * T * T ) - ( dGsdT[0] + dGsdT[3] - ( dGsdT[5] + dGsdT[5] ) ) / calR / T;
  dKcdT[4] = 1.0E+0 / T + dGs[4] / ( calR * T * T ) - ( dGsdT[2] - ( dGsdT[4] + dGsdT[4] ) ) / calR / T;
  dKcdT[5] = 1.0E+0 / T + dGs[5] / ( calR * T * T ) - ( dGsdT[3] - ( dGsdT[4] + dGsdT[5] ) ) / calR / T;
  dKcdT[6] = 1.0E+0 / T + dGs[6] / ( calR * T * T ) - ( dGsdT[6] - ( dGsdT[4] + dGsdT[1] ) ) / calR / T;
  dKcdT[7] = dGs[7] / ( calR * T * T ) - ( dGsdT[5] + dGsdT[5] - ( dGsdT[6] + dGsdT[4] ) ) / calR / T;
  dKcdT[8] = dGs[8] / ( calR * T * T ) - ( dGsdT[2] + dGsdT[1] - ( dGsdT[6] + dGsdT[4] ) ) / calR / T;
  dKcdT[9] = dGs[9] / ( calR * T * T ) - ( dGsdT[1] + dGsdT[5] - ( dGsdT[6] + dGsdT[0] ) ) / calR / T;
  dKcdT[10] = dGs[10] / ( calR * T * T ) - ( dGsdT[3] + dGsdT[1] - ( dGsdT[6] + dGsdT[5] ) ) / calR / T;
  dKcdT[11] = -1.0E+0 / T + dGs[11] / ( calR * T * T ) - ( dGsdT[5] + dGsdT[5] - ( dGsdT[11] ) ) / calR / T;
  dKcdT[12] = dGs[12] / ( calR * T * T ) - ( dGsdT[8] + dGsdT[4] - ( dGsdT[7] + dGsdT[5] ) ) / calR / T;
  dKcdT[13] = 1.0E+0 / T + dGs[13] / ( calR * T * T ) - ( dGsdT[8] - ( dGsdT[7] + dGsdT[0] ) ) / calR / T;
  dKcdT[14] = dGs[14] / ( calR * T * T ) - ( dGsdT[7] + dGsdT[2] - ( dGsdT[12] + dGsdT[4] ) ) / calR / T;
  dKcdT[15] = dGs[15] / ( calR * T * T ) - ( dGsdT[7] + dGsdT[5] - ( dGsdT[12] + dGsdT[0] ) ) / calR / T;
  dKcdT[16] = dGs[16] / ( calR * T * T ) - ( dGsdT[7] + dGsdT[3] - ( dGsdT[12] + dGsdT[5] ) ) / calR / T;
  dKcdT[17] = dGs[17] / ( calR * T * T ) - ( dGsdT[7] + dGsdT[6] - ( dGsdT[12] + dGsdT[1] ) ) / calR / T;
  dKcdT[18] = -1.0E+0 / T + dGs[18] / ( calR * T * T ) - ( dGsdT[7] + dGsdT[4] - ( dGsdT[12] ) ) / calR / T;
  dKcdT[19] = dGs[19] / ( calR * T * T ) - ( dGsdT[12] + dGsdT[2] - ( dGsdT[13] + dGsdT[4] ) ) / calR / T;
  dKcdT[20] = dGs[20] / ( calR * T * T ) - ( dGsdT[12] + dGsdT[5] - ( dGsdT[13] + dGsdT[0] ) ) / calR / T;
  dKcdT[21] = dGs[21] / ( calR * T * T ) - ( dGsdT[12] + dGsdT[3] - ( dGsdT[13] + dGsdT[5] ) ) / calR / T;
  dKcdT[22] = dGs[22] / ( calR * T * T ) - ( dGsdT[12] + dGsdT[6] - ( dGsdT[13] + dGsdT[1] ) ) / calR / T;
  dKcdT[23] = dGs[23] / ( calR * T * T ) - ( dGsdT[12] + dGsdT[10] - ( dGsdT[13] + dGsdT[9] ) ) / calR / T;
  dKcdT[24] = -1.0E+0 / T + dGs[24] / ( calR * T * T ) - ( dGsdT[12] + dGsdT[4] - ( dGsdT[13] ) ) / calR / T;
  dKcdT[25] = dGs[25] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[4] - ( dGsdT[9] + dGsdT[0] ) ) / calR / T;
  dKcdT[26] = dGs[26] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[2] - ( dGsdT[9] + dGsdT[5] ) ) / calR / T;
  dKcdT[27] = dGs[27] / ( calR * T * T ) - ( dGsdT[14] + dGsdT[0] - ( dGsdT[9] + dGsdT[1] ) ) / calR / T;
  dKcdT[28] = dGs[28] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[5] - ( dGsdT[9] + dGsdT[1] ) ) / calR / T;
  dKcdT[29] = dGs[29] / ( calR * T * T ) - ( dGsdT[14] + dGsdT[5] - ( dGsdT[9] + dGsdT[6] ) ) / calR / T;
  dKcdT[30] = dGs[30] / ( calR * T * T ) - ( dGsdT[10] + dGsdT[7] - ( dGsdT[9] + dGsdT[12] ) ) / calR / T;
  dKcdT[31] = dGs[31] / ( calR * T * T ) - ( dGsdT[17] + dGsdT[4] - ( dGsdT[9] + dGsdT[9] ) ) / calR / T;
  dKcdT[32] = -1.0E+0 / T + dGs[32] / ( calR * T * T ) - ( dGsdT[9] + dGsdT[4] - ( dGsdT[10] ) ) / calR / T;
  dKcdT[33] = dGs[33] / ( calR * T * T ) - ( dGsdT[9] + dGsdT[2] - ( dGsdT[10] + dGsdT[4] ) ) / calR / T;
  dKcdT[34] = dGs[34] / ( calR * T * T ) - ( dGsdT[9] + dGsdT[5] - ( dGsdT[10] + dGsdT[0] ) ) / calR / T;
  dKcdT[35] = dGs[35] / ( calR * T * T ) - ( dGsdT[9] + dGsdT[6] - ( dGsdT[10] + dGsdT[1] ) ) / calR / T;
  dKcdT[36] = dGs[36] / ( calR * T * T ) - ( dGsdT[9] + dGsdT[3] - ( dGsdT[10] + dGsdT[5] ) ) / calR / T;
  dKcdT[37] = dGs[37] / ( calR * T * T ) - ( dGsdT[9] + dGsdT[11] - ( dGsdT[10] + dGsdT[6] ) ) / calR / T;
  dKcdT[38] = dGs[38] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[2] - ( dGsdT[14] + dGsdT[4] ) ) / calR / T;
  dKcdT[39] = dGs[39] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[3] - ( dGsdT[14] + dGsdT[5] ) ) / calR / T;
  dKcdT[40] = dGs[40] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[6] - ( dGsdT[14] + dGsdT[1] ) ) / calR / T;
  dKcdT[41] = -1.0E+0 / T + dGs[41] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[4] - ( dGsdT[14] ) ) / calR / T;
  dKcdT[42] = dGs[42] / ( calR * T * T ) - ( dGsdT[13] + dGsdT[12] - ( dGsdT[15] + dGsdT[1] ) ) / calR / T;
  dKcdT[43] = dGs[43] / ( calR * T * T ) - ( dGsdT[15] + dGsdT[2] - ( dGsdT[16] + dGsdT[4] ) ) / calR / T;
  dKcdT[44] = dGs[44] / ( calR * T * T ) - ( dGsdT[15] + dGsdT[3] - ( dGsdT[16] + dGsdT[5] ) ) / calR / T;
  dKcdT[45] = dGs[45] / ( calR * T * T ) - ( dGsdT[16] + dGsdT[2] - ( dGsdT[17] + dGsdT[4] ) ) / calR / T;
  dKcdT[46] = dGs[46] / ( calR * T * T ) - ( dGsdT[16] + dGsdT[6] - ( dGsdT[17] + dGsdT[1] ) ) / calR / T;
  dKcdT[47] = -1.0E+0 / T + dGs[47] / ( calR * T * T ) - ( dGsdT[16] + dGsdT[4] - ( dGsdT[17] ) ) / calR / T;
  dKcdT[48] = -1.0E+0 / T + dGs[48] / ( calR * T * T ) - ( dGsdT[9] + dGsdT[9] - ( dGsdT[18] ) ) / calR / T;
  dKcdT[49] = dGs[49] / ( calR * T * T ) - ( dGsdT[17] + dGsdT[2] - ( dGsdT[18] + dGsdT[4] ) ) / calR / T;
  dKcdT[50] = dGs[50] / ( calR * T * T ) - ( dGsdT[17] + dGsdT[3] - ( dGsdT[18] + dGsdT[5] ) ) / calR / T;
  dKcdT[51] = dGs[51] / ( calR * T * T ) - ( dGsdT[17] + dGsdT[10] - ( dGsdT[18] + dGsdT[9] ) ) / calR / T;

/***********************************/

  react[0] = kf[0] / Kc[0] * X[5] * X[0] * dKcdT[0];
  react[1] = kf[1] / Kc[1] * X[5] * X[4] * dKcdT[1];
  react[2] = kf[2] / Kc[2] * X[3] * X[4] * dKcdT[2];
  react[3] = kf[3] / Kc[3] * X[0] * X[3] * dKcdT[3];
  react[4] = kf[4] / Kc[4] * X[2] * dKcdT[4] * third;
  react[5] = kf[5] / Kc[5] * X[3] * dKcdT[5] * third5;
  react[6] = kf[6] / Kc[6] * X[6] * dKcdT[6] * third;
  react[7] = kf[7] / Kc[7] * X[5] * X[5] * dKcdT[7];
  react[8] = kf[8] / Kc[8] * X[2] * X[1] * dKcdT[8];
  react[9] = kf[9] / Kc[9] * X[1] * X[5] * dKcdT[9];
  react[10] = kf[10] / Kc[10] * X[3] * X[1] * dKcdT[10];
  react[11] = kf[11] / Kc[11] * X[5] * X[5] * dKcdT[11] * third11;
  react[12] = kf[12] / Kc[12] * X[8] * X[4] * dKcdT[12];
  react[13] = kf[13] / Kc[13] * X[8] * dKcdT[13] * third13;
  react[14] = kf[14] / Kc[14] * X[7] * X[2] * dKcdT[14];
  react[15] = kf[15] / Kc[15] * X[7] * X[5] * dKcdT[15];
  react[16] = kf[16] / Kc[16] * X[7] * X[3] * dKcdT[16];
  react[17] = kf[17] / Kc[17] * X[7] * X[6] * dKcdT[17];
  react[18] = kf[18] / Kc[18] * X[7] * X[4] * dKcdT[18] * third18;
  react[19] = kf[19] / Kc[19] * X[12] * X[2] * dKcdT[19];
  react[20] = kf[20] / Kc[20] * X[12] * X[5] * dKcdT[20];
  react[21] = kf[21] / Kc[21] * X[12] * X[3] * dKcdT[21];
  react[22] = kf[22] / Kc[22] * X[12] * X[6] * dKcdT[22];
  react[23] = kf[23] / Kc[23] * X[12] * X[10] * dKcdT[23];
  react[24] = kf[24] / Kc[24] * X[12] * X[4] * dKcdT[24] * third24;
  react[25] = kf[25] / Kc[25] * X[13] * X[4] * dKcdT[25];
  react[26] = kf[26] / Kc[26] * X[13] * X[2] * dKcdT[26];
  react[27] = kf[27] / Kc[27] * X[14] * X[0] * dKcdT[27];
  react[28] = kf[28] / Kc[28] * X[13] * X[5] * dKcdT[28];
  react[29] = kf[29] / Kc[29] * X[14] * X[5] * dKcdT[29];
  react[30] = kf[30] / Kc[30] * X[10] * X[7] * dKcdT[30];
  react[31] = kf[31] / Kc[31] * X[17] * X[4] * dKcdT[31];
  react[32] = kf[32] / Kc[32] * X[9] * X[4] * dKcdT[32];
  react[33] = kf[33] / Kc[33] * X[9] * X[2] * dKcdT[33];
  react[34] = kf[34] / Kc[34] * X[9] * X[5] * dKcdT[34];
  react[35] = kf[35] / Kc[35] * X[9] * X[6] * dKcdT[35];
  react[36] = kf[36] / Kc[36] * X[9] * X[3] * dKcdT[36];
  react[37] = kf[37] / Kc[37] * X[9] * X[11] * dKcdT[37];
  react[38] = kf[38] / Kc[38] * X[13] * X[2] * dKcdT[38];
  react[39] = kf[39] / Kc[39] * X[13] * X[3] * dKcdT[39];
  react[40] = kf[40] / Kc[40] * X[13] * X[6] * dKcdT[40];
  react[41] = kf[41] / Kc[41] * X[13] * X[4] * dKcdT[41] * third41;
  react[42] = kf[42] / Kc[42] * X[13] * X[12] * dKcdT[42];
  react[43] = kf[43] / Kc[43] * X[15] * X[2] * dKcdT[43];
  react[44] = kf[44] / Kc[44] * X[15] * X[3] * dKcdT[44];
  react[45] = kf[45] / Kc[45] * X[16] * X[2] * dKcdT[45];
  react[46] = kf[46] / Kc[46] * X[16] * X[6] * dKcdT[46];
  react[47] = kf[47] / Kc[47] * X[16] * X[4] * dKcdT[47];
  react[48] = kf[48] / Kc[48] * X[9] * X[9] * dKcdT[48];
  react[49] = kf[49] / Kc[49] * X[17] * X[2] * dKcdT[49];
  react[50] = kf[50] / Kc[50] * X[17] * X[3] * dKcdT[50];
  react[51] = kf[51] / Kc[51] * X[17] * X[10] * dKcdT[51];

/* get the final answer */

  dWsdT[0] =
    dWsdT[0] + ( react[0] - react[1] + react[3] - react[9] - react[13] - react[15] - react[20] - react[25]
                 + react[27] - react[34] )
    * 1.0E+6 * SpeciesMolWeight ( 0 );
  dWsdT[1] =
    dWsdT[1] + ( -react[0] - react[6] + react[8] + react[9] + react[10] - react[17] - react[22] - react[27]
                 - react[28] - react[35] - react[40] - react[42] - react[46] )
    * 1.0E+6 * SpeciesMolWeight ( 1 );
  dWsdT[2] =
    dWsdT[2] + ( -react[1] - react[2] + react[4] + react[8] + react[14] + react[19] + react[26] + react[33]
                 + react[38] + react[43] + react[45] + react[49] )
    * 1.0E+6 * SpeciesMolWeight ( 2 );
  dWsdT[3] =
    dWsdT[3] + ( react[2] + react[3] + react[5] + react[10] + react[16] + react[21] + react[36] + react[39]
                 + react[44] + react[50] )
    * 1.0E+6 * SpeciesMolWeight ( 3 );
  dWsdT[4] =
    dWsdT[4] + ( -react[0] + react[1] + react[2] - react[4] - react[4] - react[5] - react[6] - react[7]
                 - react[8] + react[12] - react[14] + react[18] - react[19] + react[24] + react[25] +
                 react[31]
                 + react[32] - react[33] - react[38] + react[41] - react[43] - react[45] + react[47] -
                 react[49] )
    * 1.0E+6 * SpeciesMolWeight ( 4 );
  dWsdT[5] =
    dWsdT[5] + ( react[0] + react[1] - react[2] - react[3] - react[3] - react[5] + react[7] + react[7]
                 + react[9] - react[10] + react[11] + react[11] - react[12] + react[15] - react[16] +
                 react[20]
                 - react[21] - react[26] + react[28] + react[29] + react[34] - react[36] - react[39] -
                 react[44]
                 - react[50] )
    * 1.0E+6 * SpeciesMolWeight ( 5 );
  dWsdT[6] =
    dWsdT[6] + ( react[6] - react[7] - react[8] - react[9] - react[10] + react[17] + react[22] - react[29]
                 + react[35] - react[37] + react[40] + react[46] )
    * 1.0E+6 * SpeciesMolWeight ( 6 );
  dWsdT[7] =
    dWsdT[7] + ( -react[12] - react[13] + react[14] + react[15] + react[16] + react[17] + react[18] +
                 react[30] )
    * 1.0E+6 * SpeciesMolWeight ( 7 );
  dWsdT[8] = dWsdT[8] + ( react[12] + react[13] )
    * 1.0E+6 * SpeciesMolWeight ( 8 );
  dWsdT[9] =
    dWsdT[9] + ( -react[23] - react[25] - react[26] - react[27] - react[28] - react[29] - react[30] -
                 react[31]
                 - react[31] + react[32] + react[33] + react[34] + react[35] + react[36] + react[37] +
                 react[48]
                 + react[48] - react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 9 );
  dWsdT[10] =
    dWsdT[10] + ( react[23] + react[30] - react[32] - react[33] - react[34] - react[35] - react[36] -
                  react[37]
                  + react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 10 );
  dWsdT[11] = dWsdT[11] + ( -react[11] + react[37] )
    * 1.0E+6 * SpeciesMolWeight ( 11 );
  dWsdT[12] =
    dWsdT[12] + ( -react[14] - react[15] - react[16] - react[17] - react[18] + react[19] + react[20] +
                  react[21]
                  + react[22] + react[23] + react[24] - react[30] + react[42] )
    * 1.0E+6 * SpeciesMolWeight ( 12 );
  dWsdT[13] =
    dWsdT[13] + ( -react[19] - react[20] - react[21] - react[22] - react[23] - react[24] + react[25] +
                  react[26]
                  + react[28] + react[38] + react[39] + react[40] + react[41] + react[42] )
    * 1.0E+6 * SpeciesMolWeight ( 13 );
  dWsdT[14] = dWsdT[14] + ( react[27] + react[29] - react[38] - react[39] - react[40] - react[41] )
    * 1.0E+6 * SpeciesMolWeight ( 14 );
  dWsdT[15] = dWsdT[15] + ( -react[42] + react[43] + react[44] )
    * 1.0E+6 * SpeciesMolWeight ( 15 );
  dWsdT[16] = dWsdT[16] + ( -react[43] - react[44] + react[45] + react[46] + react[47] )
    * 1.0E+6 * SpeciesMolWeight ( 16 );
  dWsdT[17] =
    dWsdT[17] + ( react[31] - react[45] - react[46] - react[47] + react[49] + react[50] + react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 17 );
  dWsdT[18] = dWsdT[18] + ( -react[48] - react[49] - react[50] - react[51] )
    * 1.0E+6 * SpeciesMolWeight ( 18 );
  dWsdT[19] = 0.0E+0;

/* find the derivative of each source term wrt concentration */

  for ( k = 0; k < ns; k++ ) {
    for ( p = 0; p < ns; p++ ) {
      dWsdX[k][p] = 0.0e0;
    }
  }

  dWsdX[0][0] = 1.0E+6 * SpeciesMolWeight ( 0 ) *
    ( -1.0 * kb[0] * X[5] - 1.0 * kf[1] * X[2] - 1.0 * kb[3] * X[3] - 1.0 * kf[9] * X[6] -
      1.0 * kf[13] * X[7] * third13 - 1.0 * kf[15] * X[12] - 1.0 * kf[20] * X[13] - 1.0 * kf[25] * X[9] -
      1.0 * kb[27] * X[14] - 1.0 * kf[34] * X[10] );

  dWsdX[0][1] = 1.0E+6 * SpeciesMolWeight ( 0 ) *
    ( 1.0 * kf[0] * X[4] + 1.0 * kb[9] * X[5] + 1.0 * kf[27] * X[9] );

  dWsdX[0][2] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kf[1] * X[0] );

  dWsdX[0][3] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kb[3] * X[0] );

  dWsdX[0][4] = 1.0E+6 * SpeciesMolWeight ( 0 ) *
    ( 1.0 * kf[0] * X[1] + 1.0 * kb[1] * X[5] + 1.0 * kb[25] * X[13] );

  dWsdX[0][5] = 1.0E+6 * SpeciesMolWeight ( 0 ) *
    ( -1.0 * kb[0] * X[0] + 1.0 * kb[1] * X[4] + 2.0 * kf[3] * X[5]
      + 1.0 * kb[9] * X[1] + 1.0 * kb[15] * X[7] + 1.0 * kb[20] * X[12]
      + 1.0 * kb[34] * X[9] );

  dWsdX[0][6] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kf[9] * X[0] );

  dWsdX[0][7] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kf[13] * X[0] * third13 + 1.0 * kb[15] * X[5] );

  dWsdX[0][8] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( 1.0 * kb[13] * third13 );

  dWsdX[0][9] = 1.0E+6 * SpeciesMolWeight ( 0 ) *
    ( -1.0 * kf[25] * X[0] + 1.0 * kf[27] * X[1] + 1.0 * kb[34] * X[5] );

  dWsdX[0][10] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kf[34] * X[0] );

  dWsdX[0][12] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kf[15] * X[0] + 1.0 * kb[20] * X[5] );

  dWsdX[0][13] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kf[20] * X[0] + 1.0 * kb[25] * X[4] );

  dWsdX[0][14] = 1.0E+6 * SpeciesMolWeight ( 0 ) * ( -1.0 * kb[27] * X[0] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[0][k] = dWsdX[0][k] + 1.0E+6 * SpeciesMolWeight ( 0 ) *
      ( -1.0 * kf[13] * X[7] * X[0] + 1.0 * kb[13] * X[8] );
  }

  dWsdX[1][0] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( 1.0 * kb[0] * X[5] + 1.0 * kf[9] * X[6] + 1.0 * kb[27] * X[14] );

  dWsdX[1][1] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( -1.0 * kf[0] * X[4] - 1.0 * kf[6] * X[4] * third - 1.0 * kb[8] * X[2] - 1.0 * kb[9] * X[5] -
      1.0 * kb[10] * X[3]
      - 1.0 * kf[17] * X[12] - 1.0 * kf[22] * X[13] - 1.0 * kf[27] * X[9] - 1.0 * kf[28] * X[9] -
      1.0 * kf[35] * X[10]
      - 1.0 * kf[40] * X[14] - 1.0 * kf[42] * X[15] - 1.0 * kf[46] * X[17] );

  dWsdX[1][2] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( -1.0 * kb[8] * X[1] );

  dWsdX[1][3] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( -1.0 * kb[10] * X[1] );

  dWsdX[1][4] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( -1.0 * kf[0] * X[1] - 1.0 * kf[6] * X[1] * third + 1.0 * kf[8] * X[6] );

  dWsdX[1][5] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( 1.0 * kb[0] * X[0] - 1.0 * kb[9] * X[1] + 1.0 * kf[10] * X[6] + 1.0 * kb[28] * X[13] );

  dWsdX[1][6] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( 1.0 * kb[6] * third + 1.0 * kf[8] * X[4] + 1.0 * kf[9] * X[0] + 1.0 * kf[10] * X[5] +
      1.0 * kb[17] * X[7]
      + 1.0 * kb[22] * X[12] + 1.0 * kb[35] * X[9] + 1.0 * kb[40] * X[13] + 1.0 * kb[46] * X[16] );

  dWsdX[1][7] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( 1.0 * kb[17] * X[6] );

  dWsdX[1][9] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( -1.0 * kf[27] * X[1] - 1.0 * kf[28] * X[1] + 1.0 * kb[35] * X[6] );

  dWsdX[1][10] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( -1.0 * kf[35] * X[1] );

  dWsdX[1][12] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( -1.0 * kf[17] * X[1] + 1.0 * kb[22] * X[6] + 1.0 * kb[42] * X[13] );

  dWsdX[1][13] = 1.0E+6 * SpeciesMolWeight ( 1 ) *
    ( -1.0 * kf[22] * X[1] + 1.0 * kb[28] * X[5] + 1.0 * kb[40] * X[6] + 1.0 * kb[42] * X[12] );

  dWsdX[1][14] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( 1.0 * kb[27] * X[0] - 1.0 * kf[40] * X[1] );

  dWsdX[1][15] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( -1.0 * kf[42] * X[1] );

  dWsdX[1][16] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( 1.0 * kb[46] * X[6] );

  dWsdX[1][17] = 1.0E+6 * SpeciesMolWeight ( 1 ) * ( -1.0 * kf[46] * X[1] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[1][k] = dWsdX[1][k] + 1.0E+6 * SpeciesMolWeight ( 1 ) *
      ( -1.0 * kf[6] * X[4] * X[1] + 1.0 * kb[6] * X[6] );
  }

  dWsdX[2][0] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( -1.0 * kf[1] * X[2] );

  dWsdX[2][1] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( -1.0 * kb[8] * X[2] );

  dWsdX[2][2] = 1.0E+6 * SpeciesMolWeight ( 2 ) *
    ( -1.0 * kf[1] * X[0] - 1.0 * kf[2] * X[5] - 1.0 * kb[4] * third - 1.0 * kb[8] * X[1] -
      1.0 * kb[14] * X[7]
      - 1.0 * kb[19] * X[12] - 1.0 * kb[26] * X[13] - 1.0 * kb[33] * X[9] - 1.0 * kb[38] * X[13]
      - 1.0 * kb[43] * X[15] - 1.0 * kb[45] * X[16] - 1.0 * kb[49] * X[17] );

  dWsdX[2][3] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kb[2] * X[4] );

  dWsdX[2][4] = 1.0E+6 * SpeciesMolWeight ( 2 ) *
    ( 1.0 * kb[1] * X[5] + 1.0 * kb[2] * X[3] + 2.0 * kf[4] * X[4] * third + 1.0 * kf[8] * X[6] +
      1.0 * kf[14] * X[12]
      + 1.0 * kf[19] * X[13] + 1.0 * kf[33] * X[10] + 1.0 * kf[38] * X[14] + 1.0 * kf[43] * X[16]
      + 1.0 * kf[45] * X[17] + 1.0 * kf[49] * X[18] );

  dWsdX[2][5] = 1.0E+6 * SpeciesMolWeight ( 2 ) *
    ( 1.0 * kb[1] * X[4] - 1.0 * kf[2] * X[2] + 1.0 * kf[26] * X[9] );

  dWsdX[2][6] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[8] * X[4] );

  dWsdX[2][7] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( -1.0 * kb[14] * X[2] );

  dWsdX[2][9] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[26] * X[5] - 1.0 * kb[33] * X[2] );

  dWsdX[2][10] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[33] * X[4] );

  dWsdX[2][12] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[14] * X[4] - 1.0 * kb[19] * X[2] );

  dWsdX[2][13] = 1.0E+6 * SpeciesMolWeight ( 2 ) *
    ( 1.0 * kf[19] * X[4] - 1.0 * kb[26] * X[2] - 1.0 * kb[38] * X[2] );

  dWsdX[2][14] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[38] * X[4] );

  dWsdX[2][15] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( -1.0 * kb[43] * X[2] );

  dWsdX[2][16] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[43] * X[4] - 1.0 * kb[45] * X[2] );

  dWsdX[2][17] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[45] * X[4] - 1.0 * kb[49] * X[2] );

  dWsdX[2][18] = 1.0E+6 * SpeciesMolWeight ( 2 ) * ( 1.0 * kf[49] * X[4] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[2][k] = dWsdX[2][k] + 1.0E+6 * SpeciesMolWeight ( 2 ) *
      ( 1.0 * kf[4] * X[4] * X[4] - 1.0 * kb[4] * X[2] );
  }

  dWsdX[3][0] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( -1.0 * kb[3] * X[3] );

  dWsdX[3][1] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( -1.0 * kb[10] * X[3] );

  dWsdX[3][2] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[2] * X[5] );

  dWsdX[3][3] = 1.0E+6 * SpeciesMolWeight ( 3 ) *
    ( -1.0 * kb[2] * X[4] - 1.0 * kb[3] * X[0] - 1.0 * kb[5] * third5 - 1.0 * kb[10] * X[1] -
      1.0 * kb[16] * X[7]
      - 1.0 * kb[21] * X[12] - 1.0 * kb[36] * X[9] - 1.0 * kb[39] * X[13] - 1.0 * kb[44] * X[15]
      - 1.0 * kb[50] * X[17] );

  dWsdX[3][4] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( -1.0 * kb[2] * X[3] + 1.0 * kf[5] * X[5] * third5 );

  dWsdX[3][5] = 1.0E+6 * SpeciesMolWeight ( 3 ) *
    ( 1.0 * kf[2] * X[2] + 2.0 * kf[3] * X[5] + 1.0 * kf[5] * X[4] * third5 + 1.0 * kf[10] * X[6]
      + 1.0 * kf[16] * X[12] + 1.0 * kf[21] * X[13] + 1.0 * kf[36] * X[10] + 1.0 * kf[39] * X[14]
      + 1.0 * kf[44] * X[16] + 1.0 * kf[50] * X[18] );

  dWsdX[3][6] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[10] * X[5] );

  dWsdX[3][7] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( -1.0 * kb[16] * X[3] );

  dWsdX[3][9] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( -1.0 * kb[36] * X[3] );

  dWsdX[3][10] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[36] * X[5] );

  dWsdX[3][12] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[16] * X[5] - 1.0 * kb[21] * X[3] );

  dWsdX[3][13] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[21] * X[5] - 1.0 * kb[39] * X[3] );

  dWsdX[3][14] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[39] * X[5] );

  dWsdX[3][15] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( -1.0 * kb[44] * X[3] );

  dWsdX[3][16] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[44] * X[5] );

  dWsdX[3][17] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( -1.0 * kb[50] * X[3] );

  dWsdX[3][18] = 1.0E+6 * SpeciesMolWeight ( 3 ) * ( 1.0 * kf[50] * X[5] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[3][k] = dWsdX[3][k] + 1.0E+6 * SpeciesMolWeight ( 3 ) *
      ( 1.0 * kf[5] * X[4] * X[5] - 1.0 * kb[5] * X[3] );
  }

  dWsdX[4][0] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( 1.0 * kb[0] * X[5] + 1.0 * kf[1] * X[2] + 1.0 * kf[25] * X[9] );

  dWsdX[4][1] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( -1.0 * kf[0] * X[4] - 1.0 * kf[6] * X[4] * third + 1.0 * kb[8] * X[2] );

  dWsdX[4][2] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( 1.0 * kf[1] * X[0] + 1.0 * kf[2] * X[5] + 2.0 * kb[4] * third + 1.0 * kb[8] * X[1] + 1.0 * kb[14] * X[7]
      + 1.0 * kb[19] * X[12] + 1.0 * kb[33] * X[9] + 1.0 * kb[38] * X[13] + 1.0 * kb[43] * X[15]
      + 1.0 * kb[45] * X[16] + 1.0 * kb[49] * X[17] );

  dWsdX[4][3] = 1.0E+6 * SpeciesMolWeight ( 4 ) * ( -1.0 * kb[2] * X[4] + 1.0 * kb[5] * third5 );

  dWsdX[4][4] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( -1.0 * kf[0] * X[1] - 1.0 * kb[1] * X[5] - 1.0 * kb[2] * X[3] - 4.0 * kf[4] * X[4] * third
      - 1.0 * kf[5] * X[5] * third5 - 1.0 * kf[6] * X[1] * third - 1.0 * kf[7] * X[6] - 1.0 * kf[8] * X[6]
      - 1.0 * kb[12] * X[8] - 1.0 * kf[14] * X[12] - 1.0 * kb[18] * X[7] * third18 - 1.0 * kf[19] * X[13]
      - 1.0 * kb[24] * X[12] * third24 - 1.0 * kb[25] * X[13] - 1.0 * kb[31] * X[17] - 1.0 * kb[32] * X[9]
      - 1.0 * kf[33] * X[10] - 1.0 * kf[38] * X[14] - 1.0 * kb[41] * X[13] * third41 - 1.0 * kf[43] * X[16]
      - 1.0 * kf[45] * X[17] - 1.0 * kb[47] * X[16] - 1.0 * kf[49] * X[18] );

  dWsdX[4][5] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( 1.0 * kb[0] * X[0] - 1.0 * kb[1] * X[4] + 1.0 * kf[2] * X[2] - 1.0 * kf[5] * X[4] * third5 +
      2.0 * kb[7] * X[5]
      + 1.0 * kf[12] * X[7] );

  dWsdX[4][6] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( 1.0 * kb[6] * third - 1.0 * kf[7] * X[4] - 1.0 * kf[8] * X[4] );

  dWsdX[4][7] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( 1.0 * kf[12] * X[5] + 1.0 * kb[14] * X[2] - 1.0 * kb[18] * X[4] * third18 );

  dWsdX[4][8] = 1.0E+6 * SpeciesMolWeight ( 4 ) * ( -1.0 * kb[12] * X[4] );

  dWsdX[4][9] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( 1.0 * kf[25] * X[0] + 2.0 * kf[31] * X[9] - 1.0 * kb[32] * X[4] + 1.0 * kb[33] * X[2] );

  dWsdX[4][10] = 1.0E+6 * SpeciesMolWeight ( 4 ) * ( 1.0 * kf[32] - 1.0 * kf[33] * X[4] );

  dWsdX[4][12] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( -1.0 * kf[14] * X[4] + 1.0 * kf[18] * third18 + 1.0 * kb[19] * X[2] - 1.0 * kb[24] * X[4] * third24 );

  dWsdX[4][13] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( -1.0 * kf[19] * X[4] + 1.0 * kf[24] * third24 - 1.0 * kb[25] * X[4] + 1.0 * kb[38] * X[2]
      - 1.0 * kb[41] * X[4] * third41 );

  dWsdX[4][14] = 1.0E+6 * SpeciesMolWeight ( 4 ) * ( -1.0 * kf[38] * X[4] + 1.0 * kf[41] * third41 );

  dWsdX[4][15] = 1.0E+6 * SpeciesMolWeight ( 4 ) * ( 1.0 * kb[43] * X[2] );

  dWsdX[4][16] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( -1.0 * kf[43] * X[4] + 1.0 * kb[45] * X[2] - 1.0 * kb[47] * X[4] );

  dWsdX[4][17] = 1.0E+6 * SpeciesMolWeight ( 4 ) *
    ( -1.0 * kb[31] * X[4] - 1.0 * kf[45] * X[4] + 1.0 * kf[47] + 1.0 * kb[49] * X[2] );

  dWsdX[4][18] = 1.0E+6 * SpeciesMolWeight ( 4 ) * ( -1.0 * kf[49] * X[4] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[4][k] = dWsdX[4][k] + 1.0E+6 * SpeciesMolWeight ( 4 ) *
      ( -2.0 * kf[4] * X[4] * X[4] + 2.0 * kb[4] * X[2] - 1.0 * kf[5] * X[4] * X[5] + 1.0 * kb[5] * X[3]
        - 1.0 * kf[6] * X[4] * X[1] + 1.0 * kb[6] * X[6] + 1.0 * kf[18] * X[12] - 1.0 * kb[18] * X[7] * X[4]
        + 1.0 * kf[24] * X[13] - 1.0 * kb[24] * X[12] * X[4] + kf[41] * X[14] - 1.0 * kb[41] * X[13] * X[4] );
  }

  dWsdX[5][0] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( -1.0 * kb[0] * X[5] + 1.0 * kf[1] * X[2] + 2.0 * kb[3] * X[3] + 1.0 * kf[9] * X[6] +
      1.0 * kf[15] * X[12]
      + 1.0 * kf[20] * X[13] + 1.0 * kf[34] * X[10] );

  dWsdX[5][1] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( 1.0 * kf[0] * X[4] - 1.0 * kb[9] * X[5] + 1.0 * kb[10] * X[3] + 1.0 * kf[28] * X[9] );

  dWsdX[5][2] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( 1.0 * kf[1] * X[0] - 1.0 * kf[2] * X[5] + 1.0 * kb[26] * X[13] );

  dWsdX[5][3] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( 1.0 * kb[2] * X[4] + 2.0 * kb[3] * X[0] + 1.0 * kb[5] * third5 + 1.0 * kb[10] * X[1] +
      1.0 * kb[16] * X[7]
      + 1.0 * kb[21] * X[12] + 1.0 * kb[36] * X[9] + 1.0 * kb[39] * X[13] + 1.0 * kb[44] * X[15]
      + 1.0 * kb[50] * X[17] );

  dWsdX[5][4] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( 1.0 * kf[0] * X[1] - 1.0 * kb[1] * X[5] + 1.0 * kb[2] * X[3] - 1.0 * kf[5] * X[5] * third5 +
      2.0 * kf[7] * X[6]
      + 1.0 * kb[12] * X[8] );

  dWsdX[5][5] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( -1.0 * kb[0] * X[0] - 1.0 * kb[1] * X[4] - 1.0 * kf[2] * X[2] - 4.0 * kf[3] * X[5] -
      1.0 * kf[5] * X[4] * third5 - 4.0 * kb[7] * X[5] - 1.0 * kb[9] * X[1] - 1.0 * kf[10] * X[6] -
      4.0 * kb[11] * X[5] * third11 - 1.0 * kf[12] * X[7] - 1.0 * kb[15] * X[7] - 1.0 * kf[16] * X[12] -
      1.0 * kb[20] * X[12] - 1.0 * kf[21] * X[13]
      - 1.0 * kf[26] * X[9] - 1.0 * kb[28] * X[13] - 1.0 * kb[29] * X[14] - 1.0 * kb[34] * X[9] -
      1.0 * kf[36] * X[10]
      - 1.0 * kf[39] * X[14] - 1.0 * kf[44] * X[16] - 1.0 * kf[50] * X[18] );

  dWsdX[5][6] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( 2.0 * kf[7] * X[4] + 1.0 * kf[9] * X[0] - 1.0 * kf[10] * X[5] + 1.0 * kf[29] * X[9] );

  dWsdX[5][7] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( -1.0 * kf[12] * X[5] - 1.0 * kb[15] * X[5] + 1.0 * kb[16] * X[3] );

  dWsdX[5][8] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( 1.0 * kb[12] * X[4] );

  dWsdX[5][9] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( -1.0 * kf[26] * X[5] + 1.0 * kf[28] * X[1] + 1.0 * kf[29] * X[6] - 1.0 * kb[34] * X[5] +
      1.0 * kb[36] * X[3] );

  dWsdX[5][10] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( 1.0 * kf[34] * X[0] - 1.0 * kf[36] * X[5] );

  dWsdX[5][11] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( 2.0 * kf[11] * third11 );

  dWsdX[5][12] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( 1.0 * kf[15] * X[0] - 1.0 * kf[16] * X[5] - 1.0 * kb[20] * X[5] + 1.0 * kb[21] * X[3] );

  dWsdX[5][13] = 1.0E+6 * SpeciesMolWeight ( 5 ) *
    ( 1.0 * kf[20] * X[0] - 1.0 * kf[21] * X[5] + 1.0 * kb[26] * X[2] - 1.0 * kb[28] * X[5] +
      1.0 * kb[39] * X[3] );

  dWsdX[5][14] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( -1.0 * kb[29] * X[5] - 1.0 * kf[39] * X[5] );

  dWsdX[5][15] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( 1.0 * kb[44] * X[3] );

  dWsdX[5][16] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( -1.0 * kf[44] * X[5] );

  dWsdX[5][17] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( 1.0 * kb[50] * X[3] );

  dWsdX[5][18] = 1.0E+6 * SpeciesMolWeight ( 5 ) * ( -1.0 * kf[50] * X[5] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[5][k] = dWsdX[5][k] + 1.0E+6 * SpeciesMolWeight ( 5 ) *
      ( -1.0 * kf[5] * X[4] * X[5] + kb[5] * X[3] + 2.0 * kf[11] * X[11] - 2.0 * kb[11] * X[5] * X[5] );
  }

  dWsdX[6][0] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( -1.0 * kf[9] * X[6] );

  dWsdX[6][1] = 1.0E+6 * SpeciesMolWeight ( 6 ) *
    ( 1.0 * kf[6] * X[4] * third + 1.0 * kb[8] * X[2] + 1.0 * kb[9] * X[5] + 1.0 * kb[10] * X[3] +
      1.0 * kf[17] * X[12]
      + 1.0 * kf[22] * X[13] + 1.0 * kf[35] * X[10] + 1.0 * kf[40] * X[14] + 1.0 * kf[46] * X[17] );

  dWsdX[6][2] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kb[8] * X[1] );

  dWsdX[6][3] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kb[10] * X[1] );

  dWsdX[6][4] = 1.0E+6 * SpeciesMolWeight ( 6 ) *
    ( 1.0 * kf[6] * X[1] * third - 1.0 * kf[7] * X[6] - 1.0 * kf[8] * X[6] );

  dWsdX[6][5] = 1.0E+6 * SpeciesMolWeight ( 6 ) *
    ( 2.0 * kb[7] * X[5] + 1.0 * kb[9] * X[1] - 1.0 * kf[10] * X[6] + 1.0 * kb[29] * X[14] );

  dWsdX[6][6] = 1.0E+6 * SpeciesMolWeight ( 6 ) *
    ( -1.0 * kb[6] * third - 1.0 * kf[7] * X[4] - 1.0 * kf[8] * X[4] - 1.0 * kf[9] * X[0] -
      1.0 * kf[10] * X[5]
      - 1.0 * kb[17] * X[7] - 1.0 * kb[22] * X[12] - 1.0 * kf[29] * X[9] - 1.0 * kb[35] * X[9] -
      1.0 * kf[37] * X[10]
      - 1.0 * kb[40] * X[13] - 1.0 * kb[46] * X[16] );

  dWsdX[6][7] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( -1.0 * kb[17] * X[6] );

  dWsdX[6][9] = 1.0E+6 * SpeciesMolWeight ( 6 ) *
    ( -1.0 * kf[29] * X[6] - 1.0 * kb[35] * X[6] + 1.0 * kb[37] * X[11] );

  dWsdX[6][10] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kf[35] * X[1] - 1.0 * kf[37] * X[6] );

  dWsdX[6][11] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kb[37] * X[9] );

  dWsdX[6][12] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kf[17] * X[1] - 1.0 * kb[22] * X[6] );

  dWsdX[6][13] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kf[22] * X[1] - 1.0 * kb[40] * X[6] );

  dWsdX[6][14] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kb[29] * X[5] + 1.0 * kf[40] * X[1] );

  dWsdX[6][16] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( -1.0 * kb[46] * X[6] );

  dWsdX[6][17] = 1.0E+6 * SpeciesMolWeight ( 6 ) * ( 1.0 * kf[46] * X[1] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[6][k] = dWsdX[6][k] + 1.0E+6 * SpeciesMolWeight ( 6 ) *
      ( 1.0 * kf[6] * X[4] * X[1] - 1.0 * kb[6] * X[6] );
  }

  dWsdX[7][0] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( -1.0 * kf[13] * X[7] * third13 + 1.0 * kf[15] * X[12] );

  dWsdX[7][1] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( 1.0 * kf[17] * X[12] );

  dWsdX[7][2] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( -1.0 * kb[14] * X[7] );

  dWsdX[7][3] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( -1.0 * kb[16] * X[7] );

  dWsdX[7][4] = 1.0E+6 * SpeciesMolWeight ( 7 ) *
    ( 1.0 * kb[12] * X[8] + 1.0 * kf[14] * X[12] - 1.0 * kb[18] * X[7] * third18 );

  dWsdX[7][5] = 1.0E+6 * SpeciesMolWeight ( 7 ) *
    ( -1.0 * kf[12] * X[7] - 1.0 * kb[15] * X[7] + 1.0 * kf[16] * X[12] );

  dWsdX[7][6] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( -1.0 * kb[17] * X[7] );

  dWsdX[7][7] = 1.0E+6 * SpeciesMolWeight ( 7 ) *
    ( -1.0 * kf[12] * X[5] - 1.0 * kf[13] * X[0] * third13 - 1.0 * kb[14] * X[2] - 1.0 * kb[15] * X[5]
      - 1.0 * kb[16] * X[3] - 1.0 * kb[17] * X[6] - 1.0 * kb[18] * X[4] * third18 - 1.0 * kb[30] * X[10] );

  dWsdX[7][8] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( 1.0 * kb[12] * X[4] + 1.0 * kb[13] * third13 );

  dWsdX[7][9] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( 1.0 * kf[30] * X[12] );

  dWsdX[7][10] = 1.0E+6 * SpeciesMolWeight ( 7 ) * ( -1.0 * kb[30] * X[7] );

  dWsdX[7][12] = 1.0E+6 * SpeciesMolWeight ( 7 ) *
    ( 1.0 * kf[14] * X[4] + 1.0 * kf[15] * X[0] + 1.0 * kf[16] * X[5] + 1.0 * kf[17] * X[1]
      + 1.0 * kf[18] * third18 + 1.0 * kf[30] * X[9] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[7][k] = dWsdX[7][k] + 1.0E+6 * SpeciesMolWeight ( 7 ) *
      ( -1.0e0 * kf[13] * X[7] * X[0] + 1.0 * kb[13] * X[8] + 1.0 * kf[18] * X[12] -
        1.0 * kb[18] * X[7] * X[4] );
  }

  dWsdX[8][0] = 1.0E+6 * SpeciesMolWeight ( 8 ) * ( 1.0 * kf[13] * X[7] * third13 );

  dWsdX[8][4] = 1.0E+6 * SpeciesMolWeight ( 8 ) * ( -1.0 * kb[12] * X[8] );

  dWsdX[8][5] = 1.0E+6 * SpeciesMolWeight ( 8 ) * ( 1.0 * kf[12] * X[7] );

  dWsdX[8][7] = 1.0E+6 * SpeciesMolWeight ( 8 ) * ( 1.0 * kf[12] * X[5] + 1.0 * kf[13] * X[0] * third13 );

  dWsdX[8][8] = 1.0E+6 * SpeciesMolWeight ( 8 ) * ( -1.0 * kb[12] * X[4] - 1.0 * kb[13] * third13 );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[8][k] = dWsdX[8][k] + 1.0E+6 * SpeciesMolWeight ( 8 ) *
      ( 1.0e0 * kf[13] * X[7] * X[0] - 1.0 * kb[13] * X[8] );
  }

  dWsdX[9][0] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( -1.0 * kf[25] * X[9] + 1.0 * kb[27] * X[14] + 1.0 * kf[34] * X[10] );

  dWsdX[9][1] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( -1.0 * kf[27] * X[9] - 1.0 * kf[28] * X[9] + 1.0 * kf[35] * X[10] );

  dWsdX[9][2] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( 1.0 * kb[26] * X[13] - 1.0 * kb[33] * X[9] );

  dWsdX[9][3] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( -1.0 * kb[36] * X[9] );

  dWsdX[9][4] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( 1.0 * kb[25] * X[13] + 2.0 * kb[31] * X[17] - 1.0 * kb[32] * X[9] + 1.0 * kf[33] * X[10] );

  dWsdX[9][5] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( -1.0 * kf[26] * X[9] + 1.0 * kb[28] * X[13] + 1.0 * kb[29] * X[14] - 1.0 * kb[34] * X[9] +
      1.0 * kf[36] * X[10] );

  dWsdX[9][6] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( -1.0 * kf[29] * X[9] - 1.0 * kb[35] * X[9] + 1.0 * kf[37] * X[10] );

  dWsdX[9][7] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( 1.0 * kb[30] * X[10] );

  dWsdX[9][9] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( -1.0 * kf[23] * X[13] - 1.0 * kf[25] * X[0] - 1.0 * kf[26] * X[5] - 1.0 * kf[27] * X[1] -
      1.0 * kf[28] * X[1]
      - 1.0 * kf[29] * X[6] - 1.0 * kf[30] * X[12] - 4.0 * kf[31] * X[9] - 1.0 * kb[32] * X[4] -
      1.0 * kb[33] * X[2]
      - 1.0 * kb[34] * X[5] - 1.0 * kb[35] * X[6] - 1.0 * kb[36] * X[3] - 1.0 * kb[37] * X[11] -
      4.0 * kb[48] * X[9]
      - 1.0 * kf[51] * X[18] );

  dWsdX[9][10] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( 1.0 * kb[23] * X[12] + 1.0 * kb[30] * X[7] + 1.0 * kf[32] + 1.0 * kf[33] * X[4] + 1.0 * kf[34] * X[0]
      + 1.0 * kf[35] * X[1] + 1.0 * kf[36] * X[5] + 1.0 * kf[37] * X[6] + 1.0 * kb[51] * X[17] );

  dWsdX[9][11] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( -1.0 * kb[37] * X[9] );

  dWsdX[9][12] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( 1.0 * kb[23] * X[10] - 1.0 * kf[30] * X[9] );

  dWsdX[9][13] = 1.0E+6 * SpeciesMolWeight ( 9 ) *
    ( -1.0 * kf[23] * X[9] + 1.0 * kb[25] * X[4] + 1.0 * kb[26] * X[2] + 1.0 * kb[28] * X[5] );

  dWsdX[9][14] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( 1.0 * kb[27] * X[0] + 1.0 * kb[29] * X[5] );

  dWsdX[9][17] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( 2.0 * kb[31] * X[4] + 1.0 * kb[51] * X[10] );

  dWsdX[9][18] = 1.0E+6 * SpeciesMolWeight ( 9 ) * ( 2.0 * kf[48] - 1.0 * kf[51] * X[9] );

  dWsdX[10][0] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( -1.0 * kf[34] * X[10] );

  dWsdX[10][1] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( -1.0 * kf[35] * X[10] );

  dWsdX[10][2] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kb[33] * X[9] );

  dWsdX[10][3] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kb[36] * X[9] );

  dWsdX[10][4] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kb[32] * X[9] - 1.0 * kf[33] * X[10] );

  dWsdX[10][5] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kb[34] * X[9] - 1.0 * kf[36] * X[10] );

  dWsdX[10][6] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kb[35] * X[9] - 1.0 * kf[37] * X[10] );

  dWsdX[10][7] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( -1.0 * kb[30] * X[10] );

  dWsdX[10][9] = 1.0E+6 * SpeciesMolWeight ( 10 ) *
    ( 1.0 * kf[23] * X[13] + 1.0 * kf[30] * X[12] + 1.0 * kb[32] * X[4] + 1.0 * kb[33] * X[2] +
      1.0 * kb[34] * X[5]
      + 1.0 * kb[35] * X[6] + 1.0 * kb[36] * X[3] + 1.0 * kb[37] * X[11] + 1.0 * kf[51] * X[18] );

  dWsdX[10][10] = 1.0E+6 * SpeciesMolWeight ( 10 ) *
    ( -1.0 * kb[23] * X[12] - 1.0 * kb[30] * X[7] - 1.0 * kf[32] - 1.0 * kf[33] * X[4] - 1.0 * kf[34] * X[0]
      - 1.0 * kf[35] * X[1] - 1.0 * kf[36] * X[5] - 1.0 * kf[37] * X[6] - 1.0 * kb[51] * X[17] );

  dWsdX[10][11] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kb[37] * X[9] );

  dWsdX[10][12] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( -1.0 * kb[23] * X[10] + 1.0 * kf[30] * X[9] );

  dWsdX[10][13] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kf[23] * X[9] );

  dWsdX[10][17] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( -1.0 * kb[51] * X[10] );

  dWsdX[10][18] = 1.0E+6 * SpeciesMolWeight ( 10 ) * ( 1.0 * kf[51] * X[9] );

  dWsdX[11][5] = 1.0E+6 * SpeciesMolWeight ( 11 ) * ( 2.0 * kb[11] * X[5] * third11 );

  dWsdX[11][6] = 1.0E+6 * SpeciesMolWeight ( 11 ) * ( 1.0 * kf[37] * X[10] );

  dWsdX[11][9] = 1.0E+6 * SpeciesMolWeight ( 11 ) * ( -1.0 * kb[37] * X[11] );

  dWsdX[11][10] = 1.0E+6 * SpeciesMolWeight ( 11 ) * ( 1.0 * kf[37] * X[6] );

  dWsdX[11][11] = 1.0E+6 * SpeciesMolWeight ( 11 ) * ( -1.0 * kf[11] * third11 - 1.0 * kb[37] * X[9] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[11][k] = dWsdX[11][k] + 1.0E+6 * SpeciesMolWeight ( 11 ) *
      ( -1.0 * kf[11] * X[11] + 1.0 * kb[11] * X[5] * X[5] );
  }

  dWsdX[12][0] = 1.0E+6 * SpeciesMolWeight ( 12 ) * ( -1.0 * kf[15] * X[12] + 1.0 * kf[20] * X[13] );

  dWsdX[12][1] = 1.0E+6 * SpeciesMolWeight ( 12 ) *
    ( -1.0 * kf[17] * X[12] + 1.0 * kf[22] * X[13] + 1.0 * kf[42] * X[15] );

  dWsdX[12][2] = 1.0E+6 * SpeciesMolWeight ( 12 ) * ( 1.0 * kb[14] * X[7] - 1.0 * kb[19] * X[12] );

  dWsdX[12][3] = 1.0E+6 * SpeciesMolWeight ( 12 ) * ( 1.0 * kb[16] * X[7] - 1.0 * kb[21] * X[12] );

  dWsdX[12][4] = 1.0E+6 * SpeciesMolWeight ( 12 ) *
    ( -1.0 * kf[14] * X[12] + 1.0 * kb[18] * X[7] * third18 + 1.0 * kf[19] * X[13] -
      1.0 * kb[24] * X[12] * third24 );

  dWsdX[12][5] = 1.0E+6 * SpeciesMolWeight ( 12 ) *
    ( 1.0 * kb[15] * X[7] - 1.0 * kf[16] * X[12] - 1.0 * kb[20] * X[12] + 1.0 * kf[21] * X[13] );

  dWsdX[12][6] = 1.0E+6 * SpeciesMolWeight ( 12 ) * ( 1.0 * kb[17] * X[7] - 1.0 * kb[22] * X[12] );

  dWsdX[12][7] = 1.0E+6 * SpeciesMolWeight ( 12 ) *
    ( 1.0 * kb[14] * X[2] + 1.0 * kb[15] * X[5] + 1.0 * kb[16] * X[3] + 1.0 * kb[17] * X[6]
      + 1.0 * kb[18] * X[4] * third18 + 1.0 * kb[30] * X[10] );

  dWsdX[12][9] = 1.0E+6 * SpeciesMolWeight ( 12 ) * ( 1.0 * kf[23] * X[13] - 1.0 * kf[30] * X[12] );

  dWsdX[12][10] = 1.0E+6 * SpeciesMolWeight ( 12 ) * ( -1.0 * kb[23] * X[12] + 1.0 * kb[30] * X[7] );

  dWsdX[12][12] = 1.0E+6 * SpeciesMolWeight ( 12 ) *
    ( -1.0 * kf[14] * X[4] - 1.0 * kf[15] * X[0] - 1.0 * kf[16] * X[5] - 1.0 * kf[17] * X[1] -
      1.0 * kf[18] * third18 - 1.0 * kb[19] * X[2] - 1.0 * kb[20] * X[5] - 1.0 * kb[21] * X[3] -
      1.0 * kb[22] * X[6] - 1.0 * kb[23] * X[10]
      - 1.0 * kb[24] * X[4] * third24 - 1.0 * kf[30] * X[9] - 1.0 * kb[42] * X[13] );

  dWsdX[12][13] = 1.0E+6 * SpeciesMolWeight ( 12 ) *
    ( 1.0 * kf[19] * X[4] + 1.0 * kf[20] * X[0] + 1.0 * kf[21] * X[5] + 1.0 * kf[22] * X[1] +
      1.0 * kf[23] * X[9]
      + 1.0 * kf[24] * third24 - 1.0 * kb[42] * X[12] );

  dWsdX[12][15] = 1.0E+6 * SpeciesMolWeight ( 12 ) * ( 1.0 * kf[42] * X[1] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[12][k] = dWsdX[12][k] + 1.0E+6 * SpeciesMolWeight ( 12 ) *
      ( -1.0 * kf[18] * X[12] + 1.0 * kb[18] * X[7] * X[4] + 1.0 * kf[24] * X[13] -
        1.0 * kb[24] * X[12] * X[4] );
  }

  dWsdX[13][0] = 1.0E+6 * SpeciesMolWeight ( 13 ) * ( -1.0 * kf[20] * X[13] + 1.0 * kf[25] * X[9] );

  dWsdX[13][1] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( -1.0 * kf[22] * X[13] + 1.0 * kf[28] * X[9] + 1.0 * kf[40] * X[14] + 1.0 * kf[42] * X[15] );

  dWsdX[13][2] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( 1.0 * kb[19] * X[12] - 1.0 * kb[26] * X[13] - 1.0 * kb[38] * X[13] );

  dWsdX[13][3] = 1.0E+6 * SpeciesMolWeight ( 13 ) * ( 1.0 * kb[21] * X[12] - 1.0 * kb[39] * X[13] );

  dWsdX[13][4] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( -1.0 * kf[19] * X[13] + 1.0 * kb[24] * X[12] * third24 - 1.0 * kb[25] * X[13] + 1.0 * kf[38] * X[14]
      - 1.0 * kb[41] * X[13] * third41 );

  dWsdX[13][5] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( 1.0 * kb[20] * X[12] - 1.0 * kf[21] * X[13] + 1.0 * kf[26] * X[9] - 1.0 * kb[28] * X[13]
      + 1.0 * kf[39] * X[14] );

  dWsdX[13][6] = 1.0E+6 * SpeciesMolWeight ( 13 ) * ( 1.0 * kb[22] * X[12] - 1.0 * kb[40] * X[13] );

  dWsdX[13][9] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( -1.0 * kf[23] * X[13] + 1.0 * kf[25] * X[0] + 1.0 * kf[26] * X[5] + 1.0 * kf[28] * X[1] );

  dWsdX[13][10] = 1.0E+6 * SpeciesMolWeight ( 13 ) * ( 1.0 * kb[23] * X[12] );

  dWsdX[13][12] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( 1.0 * kb[19] * X[2] + 1.0 * kb[20] * X[5] + 1.0 * kb[21] * X[3] + 1.0 * kb[22] * X[6] +
      1.0 * kb[23] * X[10]
      + 1.0 * kb[24] * X[4] * third24 - 1.0 * kb[42] * X[13] );

  dWsdX[13][13] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( -1.0 * kf[19] * X[4] - 1.0 * kf[20] * X[0] - 1.0 * kf[21] * X[5] - 1.0 * kf[22] * X[1] -
      1.0 * kf[23] * X[9]
      - 1.0 * kf[24] * third24 - 1.0 * kb[25] * X[4] - 1.0 * kb[26] * X[2] - 1.0 * kb[28] * X[5] -
      1.0 * kb[38] * X[2]
      - 1.0 * kb[39] * X[3] - 1.0 * kb[40] * X[6] - 1.0 * kb[41] * X[4] * third41 - 1.0 * kb[42] * X[12] );

  dWsdX[13][14] = 1.0E+6 * SpeciesMolWeight ( 13 ) *
    ( 1.0 * kf[38] * X[4] + 1.0 * kf[39] * X[5] + 1.0 * kf[40] * X[1] + 1.0 * kf[41] * third41 );

  dWsdX[13][15] = 1.0E+6 * SpeciesMolWeight ( 13 ) * ( 1.0 * kf[42] * X[1] );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[13][k] = dWsdX[13][k] + 1.0E+6 * SpeciesMolWeight ( 13 ) *
      ( -1.0 * kf[24] * X[13] + 1.0 * kb[24] * X[12] * X[4] + 1.0 * kf[41] * X[14] - kb[41] * X[13] * X[4] );
  }

  dWsdX[14][0] = 1.0E+6 * SpeciesMolWeight ( 14 ) * ( -1.0 * kb[27] * X[14] );

  dWsdX[14][1] = 1.0E+6 * SpeciesMolWeight ( 14 ) * ( 1.0 * kf[27] * X[9] - 1.0 * kf[40] * X[14] );

  dWsdX[14][2] = 1.0E+6 * SpeciesMolWeight ( 14 ) * ( 1.0 * kb[38] * X[13] );

  dWsdX[14][3] = 1.0E+6 * SpeciesMolWeight ( 14 ) * ( 1.0 * kb[39] * X[13] );

  dWsdX[14][4] = 1.0E+6 * SpeciesMolWeight ( 14 ) *
    ( -1.0 * kf[38] * X[14] + 1.0 * kb[41] * X[13] * third41 );

  dWsdX[14][5] = 1.0E+6 * SpeciesMolWeight ( 14 ) * ( -1.0 * kb[29] * X[14] - 1.0 * kf[39] * X[14] );

  dWsdX[14][6] = 1.0E+6 * SpeciesMolWeight ( 14 ) * ( 1.0 * kf[29] * X[9] + 1.0 * kb[40] * X[13] );

  dWsdX[14][9] = 1.0E+6 * SpeciesMolWeight ( 14 ) * ( 1.0 * kf[27] * X[1] + 1.0 * kf[29] * X[6] );

  dWsdX[14][13] = 1.0E+6 * SpeciesMolWeight ( 14 ) *
    ( 1.0 * kb[38] * X[2] + 1.0 * kb[39] * X[3] + 1.0 * kb[40] * X[6] + 1.0 * kb[41] * X[4] * third41 );

  dWsdX[14][14] = 1.0E+6 * SpeciesMolWeight ( 14 ) *
    ( -1.0 * kb[27] * X[0] - 1.0 * kb[29] * X[5] - 1.0 * kf[38] * X[4] - 1.0 * kf[39] * X[5] -
      1.0 * kf[40] * X[1]
      - 1.0 * kf[41] * third41 );

  for ( k = 0; k < ns; k++ ) {
    dWsdX[14][k] = dWsdX[14][k] + 1.0E+6 * SpeciesMolWeight ( 14 ) *
      ( -1.0 * kf[41] * X[14] + 1.0 * kb[41] * X[13] * X[4] );
  }

  dWsdX[15][1] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( -1.0 * kf[42] * X[15] );

  dWsdX[15][2] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( -1.0 * kb[43] * X[15] );

  dWsdX[15][3] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( -1.0 * kb[44] * X[15] );

  dWsdX[15][4] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( 1.0 * kf[43] * X[16] );

  dWsdX[15][5] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( 1.0 * kf[44] * X[16] );

  dWsdX[15][12] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( 1.0 * kb[42] * X[13] );

  dWsdX[15][13] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( 1.0 * kb[42] * X[12] );

  dWsdX[15][15] = 1.0E+6 * SpeciesMolWeight ( 15 ) *
    ( -1.0 * kf[42] * X[1] - 1.0 * kb[43] * X[2] - 1.0 * kb[44] * X[3] );

  dWsdX[15][16] = 1.0E+6 * SpeciesMolWeight ( 15 ) * ( 1.0 * kf[43] * X[4] + 1.0 * kf[44] * X[5] );

  dWsdX[16][1] = 1.0E+6 * SpeciesMolWeight ( 16 ) * ( 1.0 * kf[46] * X[17] );

  dWsdX[16][2] = 1.0E+6 * SpeciesMolWeight ( 16 ) * ( 1.0 * kb[43] * X[15] - 1.0 * kb[45] * X[16] );

  dWsdX[16][3] = 1.0E+6 * SpeciesMolWeight ( 16 ) * ( 1.0 * kb[44] * X[15] );

  dWsdX[16][4] = 1.0E+6 * SpeciesMolWeight ( 16 ) *
    ( -1.0 * kf[43] * X[16] + 1.0 * kf[45] * X[17] - 1.0 * kb[47] * X[16] );

  dWsdX[16][5] = 1.0E+6 * SpeciesMolWeight ( 16 ) * ( -1.0 * kf[44] * X[16] );

  dWsdX[16][6] = 1.0E+6 * SpeciesMolWeight ( 16 ) * ( -1.0 * kb[46] * X[16] );

  dWsdX[16][15] = 1.0E+6 * SpeciesMolWeight ( 16 ) * ( 1.0 * kb[43] * X[2] + 1.0 * kb[44] * X[3] );

  dWsdX[16][16] = 1.0E+6 * SpeciesMolWeight ( 16 ) *
    ( -1.0 * kf[43] * X[4] - 1.0 * kf[44] * X[5] - 1.0 * kb[45] * X[2] - 1.0 * kb[46] * X[6] -
      1.0 * kb[47] * X[4] );

  dWsdX[16][17] = 1.0E+6 * SpeciesMolWeight ( 16 ) *
    ( 1.0 * kf[45] * X[4] + 1.0 * kf[46] * X[1] + 1.0 * kf[47] );

  dWsdX[17][1] = 1.0E+6 * SpeciesMolWeight ( 17 ) * ( -1.0 * kf[46] * X[17] );

  dWsdX[17][2] = 1.0E+6 * SpeciesMolWeight ( 17 ) * ( 1.0 * kb[45] * X[16] - 1.0 * kb[49] * X[17] );

  dWsdX[17][3] = 1.0E+6 * SpeciesMolWeight ( 17 ) * ( -1.0 * kb[50] * X[17] );

  dWsdX[17][4] = 1.0E+6 * SpeciesMolWeight ( 17 ) *
    ( -1.0 * kb[31] * X[17] - 1.0 * kf[45] * X[17] + 1.0 * kb[47] * X[16] + 1.0 * kf[49] * X[18] );

  dWsdX[17][5] = 1.0E+6 * SpeciesMolWeight ( 17 ) * ( 1.0 * kf[50] * X[18] );

  dWsdX[17][6] = 1.0E+6 * SpeciesMolWeight ( 17 ) * ( 1.0 * kb[46] * X[16] );

  dWsdX[17][9] = 1.0E+6 * SpeciesMolWeight ( 17 ) * ( 2.0 * kf[31] * X[9] + 1.0 * kf[51] * X[18] );

  dWsdX[17][10] = 1.0E+6 * SpeciesMolWeight ( 17 ) * ( -1.0 * kb[51] * X[17] );

  dWsdX[17][16] = 1.0E+6 * SpeciesMolWeight ( 17 ) *
    ( 1.0 * kb[45] * X[2] + 1.0 * kb[46] * X[6] + 1.0 * kb[47] * X[4] );

  dWsdX[17][17] = 1.0E+6 * SpeciesMolWeight ( 17 ) *
    ( -1.0 * kb[31] * X[4] - 1.0 * kf[45] * X[4] - 1.0 * kf[46] * X[1] - 1.0 * kf[47] - 1.0 * kb[49] * X[2]
      - 1.0 * kb[50] * X[3] - 1.0 * kb[51] * X[10] );

  dWsdX[17][18] = 1.0E+6 * SpeciesMolWeight ( 17 ) *
    ( 1.0 * kf[49] * X[4] + 1.0 * kf[50] * X[5] + 1.0 * kf[51] * X[9] );

  dWsdX[18][2] = 1.0E+6 * SpeciesMolWeight ( 18 ) * ( 1.0 * kb[49] * X[17] );

  dWsdX[18][3] = 1.0E+6 * SpeciesMolWeight ( 18 ) * ( 1.0 * kb[50] * X[17] );

  dWsdX[18][4] = 1.0E+6 * SpeciesMolWeight ( 18 ) * ( -1.0 * kf[49] * X[18] );

  dWsdX[18][5] = 1.0E+6 * SpeciesMolWeight ( 18 ) * ( -1.0 * kf[50] * X[18] );

  dWsdX[18][9] = 1.0E+6 * SpeciesMolWeight ( 18 ) * ( 2.0 * kb[48] * X[9] - 1.0 * kf[51] * X[18] );

  dWsdX[18][10] = 1.0E+6 * SpeciesMolWeight ( 18 ) * ( 1.0 * kb[51] * X[17] );

  dWsdX[18][17] = 1.0E+6 * SpeciesMolWeight ( 18 ) *
    ( 1.0 * kb[49] * X[2] + 1.0 * kb[50] * X[3] + 1.0 * kb[51] * X[10] );

  dWsdX[18][18] = 1.0E+6 * SpeciesMolWeight ( 18 ) *
    ( -1.0 * kf[48] - 1.0 * kf[49] * X[4] - 1.0 * kf[50] * X[5] - 1.0 * kf[51] * X[9] );

/* assemble the terms for the chemical Jacobian */
  find_dT_dx ( np, gl, &dTdrhoE, dTdrhok, dTdrhoV );
  for ( k = 0; k < ns; k++ ) {
    /* derivative of S[k] wrt to rhok[ns] */
    for ( p = 0; p < ns; p++ ) {
      C[k][p] = ( dWsdX[k][p] * _dXdrhok ( p ) + dWsdT[k] * dTdrhok[p] );
    }
    /* derivative of S[k] wrt to rhov[nd] */
    for ( p = ns; p < ns + nd; p++ ) {
      C[k][p] = ( dWsdT[k] * dTdrhoV[p - ns] );
    }
    /* derivative of S[k] wrt to rhoE */
    C[k][ns + nd] = ( dWsdT[k] * dTdrhoE );
    /* derivative of S[k] wrt to rhok */
    C[k][ns + nd + 1] = ( -dTdrhoE * dWsdT[k] );
  }

  /* This section checks the analytical with the numerical Jacobian */
  /*
     find_numerical_jacobian_2(&(find_Schem), np, C2);
     for (k=0; k<ns+nd+1; k++){
     for (p=0; p<ns+nd+1; p++){
     printf("\nAna C[%ld][%ld] = %.7e\tNum C2[%ld][%ld] = %.7e\tDiff = %.8e",k,p,C[k][p],k,p,C2[k][p],
     (C[k][p]-C2[k][p])/C2[k][p]);
     }
     }
   */

}                               /*end of find_dSchem_dU */
