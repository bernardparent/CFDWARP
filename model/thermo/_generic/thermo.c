// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2018,2020 Bernard Parent
Copyright 2020 Aaron Trinh
Copyright 2001 Jason Etele
Copyright 2000 Giovanni Fusina


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
#include <soap.h>

#define T_accuracy 1.0e-12
  #define ns_c 54

#define RN2 296.8E0
#define Thetav 3353.0E0
#define evlim 0.1
#define Tvlim (Thetav/log(RN2*Thetav/evlim+1.0))

// temperature overlap of polynomials for cpk, hk, sk dskdT (in Kelvins) 
#define dToverlap 10.0e0
// minimum range of temperature within one set to cpk, hk polynomial
// make sure dTrangemin>dToverlap
#define dTrangemin 100.0e0
// maximum temperature in Kelvin used to determine the viscosity, thermal conductivity, and mass diffusion coefficients
#define MAX_T_FOR_NUK_ETA_KAPPA_POLYNOMIALS 12000.0  

typedef char speciesname_t[ns_c];

const static double Pd[3][5]=
 {
   {2.3527333E+0, -1.3589968E+0, 5.2202460E-1, -9.4262883E-2, 6.4354629E-3},
   {1.2660308E+0, -1.6441443E-1, 2.2945928E-2, -1.6324168E-3, 4.5833672E-5},
   {8.5263337E-1, -1.3552911E-2, 2.6162080E-4, -2.4647654E-6, 8.6538568E-9}
 };

const static double Pe[3][3]=
 {
   {1.1077725E+0, -9.4802344E-3, +1.6918277E-3},
   {1.0871429E+0, +3.1964282E-3, -8.9285689E-5},
   {1.1059000E+0, +6.5136364E-4, -3.4090910E-6}
 };

const static speciesname_t speciesname[ns_c]=
  {
   "e-",
   "O2",
   "N2",
   "O",
   "N",
   "O2+",
   "N2+",
   "O+",
   "N+",
   "O2-",
   "O-",
   "O3",
   "N2(A)",
   "NO",
   "NO+",
   "H2",
   "H2O",
   "H",
   "OH",
   "HO2",
   "CO",
   "CO2",
   "CH3",
   "CH4",
   "H2O2",
   "CHO",
   "CH2O",
   "CH3O",
   "C2H3",
   "C2H4",
   "C2H5",
   "C2H6",
   "CH",
   "NH",
   "C2H2",
   "C3H8",
   "C12H23",
   "He",
   "Air",
   "H2+",
   "Cs",
   "Cs+",
   "Ar",
   "C",
   "C2",
   "CN",
   "NCO",
   "Ar+",
   "C+",
   "C2+",
   "CN+",
   "CO+",
   "HNO",
   "NO2"
  };


const static long numatoms[ns_c]=
  {
   0,  /* e- */
   2,  /* O2 */
   2,  /* N2 */
   1,  /* O */
   1,  /* N */
   2,  /* O2+ */
   2,  /* N2+ */
   1,  /* O+ */
   1,  /* N+ */
   2,  /* O2- */
   1,  /* O- */
   3,  /* O3 */
   2,  /* N2(A) */
   2,  /* NO */
   2,  /* NO+ */
   2,  /* H2 */
   3,  /*H2O */
   1,  /*H */
   2,  /*OH*/
   3,  /*HO2*/
   2,  /*CO*/
   3,  /*CO2*/
   4,  /* CH3*/
   5,  /* CH4 */
   4,  /*H2O2*/
   3,  /*CHO */
   4,  /*CH2O*/
   5,  /*CH3O*/
   5,  /*C2H3*/
   6,  /*C2H4*/
   7,  /*C2H5*/
   8,  /*C2H6*/
   2,  /*CH*/  
   2,  /*NH*/
   4,  /*C2H2*/
   11, /*C3H8*/
   35, /*C12H23*/
   1,  /*He*/
   2,  /*Air*/
   2,  /*H2+*/
   1,  /*Cs*/
   1,  /*Cs+*/
   1,  /*Ar*/
   1,  /*C*/
   2,  /*C2*/
   2,  /*CN*/
   3,  /*NCO*/
   1,  /*Ar+*/
   1,  /*C+*/ 
   2,  /*C2+*/
   2,  /*CN+*/
   2,  /*CO+*/
   3,  /*HNO*/
   3   /*NO2*/
  }; 



/*
transport properties: Lennard-Jones Potential Parameters: epsilon & sigma

Primary Reference: all species except those not mentioned in the secondary reference
Svehla, R.A,"Estimated Viscosities and Thermal Conductivites of Gases at high temperature,"
	NASA TR R-132, 1962
CFDWARP/model/thermo/_generic/ref/NASA-TR-R-132.pdf

Secondary Reference: HO2, HCO, C2H3, C2H5, C4H8, CH2O, CH3, CH3O, HNO, NO2
Sandia National Laboratories, "A Fortran Computer Code Package for the Evaluation of Gas Phase
	Multicomponent Transport Properties," SAND86-8246, 1986
Sandia National Laboratories: SAND86-8246 Update/Revision, 1998
CFDWARP/model/thermo/_generic/ref/SAND86-8246.pdf

the transport properties of CHO, H2CO, C4H8O, C6H12O, C8H16, and C12H24
are calculated as           HCO, CH2O, C4H8, C2H5OC2H5, C6H12, and C6H12 respectively

*/

/*The epsilon parameter values below are sourced from one of two tables above and include
 * Boltzmann's constant inside already. Thus the units are in Kelvin. (epsilon/k_b)
 * 
 * Units of sigma: A/10 (dA, deciangstrom)
 * 
 * */

const static double Peps[ns_c]=
  {
   106.7e0,     /* e- */ /* !! unknown value: fixed to the one of O */
   106.7e0,     /* O2 */
   71.4e0,      /* N2 */
   106.7e0,     /* O */
   71.4e0,      /* N */
   106.7e0,     /* O2+ */ /* !! unknown value: fixed to the one of O2 */
   71.4e0,      /* N2+ */ /* !! unknown value: fixed to the one of N2 */
   106.7e0,     /* O+ */  /* !! unknown value: fixed to the one of O */
   71.4e0,      /* N+ */  /* !! unknown value: fixed to the one of N */
   106.7e0,     /* O2- */ /* !! unknown value: fixed to the one of O2 */
   106.7e0,     /* O- */  /* !! unknown value: fixed to the one of O */
   106.7,       /* O3 */  /* !! unknown value: fixed to the one of O2 */
   71.4,        /* N2(A) */  /* !! unknown value: fixed to the one of N2 */
   116.7e0,     /* NO */
   116.7e0,     /* NO+ */  /* !! unknown value: fixed to the one of NO */
   59.7e0,      /* H2 */
   809.1e0,     /*H2O */
   37.0e0,      /*H */
   79.8e0,      /*OH*/
   107.4e0,     /*HO2*/
   91.7e0,      /*CO*/
   195.2e0,     /*CO2*/
   144.0e0,     /* CH3*/
   148.6e0,     /* CH4 */
   289.3E+0,    /*H2O2*/
   498.0E+0,    /*CHO */
   498.0E+0,    /*CH2O*/
   417.0E+0,    /*CH3O*/
   209.0E+0,    /*C2H3*/
   224.7E+0,    /*C2H4*/
   252.3E+0,    /*C2H5*/
   215.7E+0,    /*C2H6*/
   68.6E+0,     //CH from NASA TR R-132.pdf
   65.3E+0,     //NH from NASA TR R-132.pdf
   231.8E+0,    //C2H2
   237.1E+0,    //C3H8 from NASA TR R-132.pdf
   297.1E+0,    //C12H23 from NASA TR R-132.pdf
   10.22E+0,    //He
   78.46e0,     /* Air -> obtained from N2 and O2 assuming 4:1 ratio */
   59.7e0,      /* H2+ (unknown!, fixed to the one of H2) */
   71.4,        /* Cs (unknown!, fixed to one of N2) */
   71.4,        /* Cs+ (unknown!, fixed to one of N2*/
   93.3,        /*Ar*/
   30.6,        /*C*/ 
   78.8,        /*C2*/ 
   75.0,        /*CN*/
   232.4,       /*NCO*/  
   93.3,        /*Ar+*/  /* !! unknown value: fixed to the one of Ar */
   30.6,        /*C+*/  /* !! unknown value: fixed to the one of C */
   78.8,        /*C2+*/  /* !! unknown value: fixed to the one of C2 */
   75.0,        /*CN+*/  /* !! unknown value: fixed to the one of CN */
   91.7,        /*CO+*/  /* !! unknown value: fixed to the one of CO */
   116.7,        /*HNO*/
   200.0,        /*NO2*/
  };

const static double Psig[ns_c]=
  {
   0.34e0,      /* e- */   /* !! unknown value: fixed to the one of O */
   0.3467e0,    /* O2 */
   0.3798e0,    /* N2 */
   0.34e0,      /* O */
   0.3298e0,    /* N */
   0.84e0,      /* O2+ */  /* !! determined from the mobility at 300K */
   1.22e0,      /* N2+ */  /* !! determined from the mobility at 300K */
   0.34e0,      /* O+ */   /* !! unknown value: fixed to the one of O */
   0.3298e0,    /* N+ */   /* !! unknown value: fixed to the one of N */
   0.99e0,      /* O2- */  /* !! determined from the mobility at 300K */
   0.34e0,      /* O- */  /* !! unknown value: fixed to the one of O */
   0.3467e0,    /* O3 */   /* !! unknown value: fixed to the one of O2 */
   0.3798e0,    /* N2(A) */  /* !! unknown value: fixed to the one of N2 */
   0.3492e0,    /* NO */   
   0.3492e0,    /* NO+ */  /* !! unknown value: fixed to the one of NO */
   0.2827E+0,   /* H2 */
   0.2641E+0,   /*H2O */
   0.2708E+0,   /*H */
   0.3147E+0,   /*OH*/
   0.3458E+0,   /*HO2*/
   0.3690E+0,   /*CO*/
   0.3941E+0,   /*CO2*/
   0.3800E+0,   /* CH3*/
   0.3758E+0,   /* CH4 */
   0.4196E+0,   /*H2O2*/
   0.3590E+0,   /*CHO */
   0.3590E+0,   /*CH2O*/
   0.3690E+0,   /*CH3O*/
   0.4100E+0,   /*C2H3*/
   0.4163E+0,   /*C2H4*/
   0.4302E+0,   /*C2H5*/
   0.4443E+0,   /*C2H6*/
   0.3370E+0,   //CH from NASA TR R-132.pdf
   0.3312E+0,   //NH from NASA TR R-132.pdf
   0.4033E+0,   //C2H2
   0.5118E+0,   //C3H8 from NASA TR R-132.pdf
   0.6182E+0,   //C12H23 from NASA TR R-132.pdf     
   0.2551E+0,   //He
   0.3732e0,    /* Air -> obtained from N2 and O2 assuming a 4:1 ratio */
   0.2827E+0,   /* H2+ !! fixed to the one of H2 */
   0.3567E+0,   /* Cs, unknown!, fixed to the one of Na */
   0.3567E+0,   /* Cs+, unknown!, fixed to the one of Na */
   0.3542E+0,   /*Ar*/
   0.3385E+0,   /*C*/
   0.3913E+0,   /*C2*/
   0.3856E+0,   /*CN*/
   0.3828E+0,   /*NCO*/
   0.3543E+0,   /* Ar+ !! fixed to the one of Ar */   
   0.3385E+0,   /*C+ !! fixed to the one of C */
   0.3913E+0,   /*C2+* !! fixed to the one of C2 */
   0.3856E+0,   /*CN+* !! fixed to the one of CN */
   0.3690E+0,   /*CO+* !! fixed to the one of CO */  
   0.3492E+0,   /*HNO*/
   0.2500E+0,   /*NO2*/
  };


/*
Molecular Weight Reference:
NASA RP-1311 volume I & II, 1994 & 1996

Units: kg/mol

*/     

const static double calM[ns_c]=
  {
   5.4858E-7,     /* e- */
   31.9988e-3,    /* O2 */
   28.0134e-3,    /* N2 */
   15.9994e-3,    /* O */
   14.0067e-3,    /* N */
   31.99825e-3,   /* O2+ */
   28.01285e-3,   /* N2+ */
   15.99885e-3,   /* O+ */
   14.00615e-3,   /* N+ */
   31.99935e-3,   /* O2- */   
   15.99995e-3,   /* O- */    
   47.99820e-3,   /* O3 */    
   28.01340e-3,   /* N2(A) */    
   30.00610e-3,   /* NO */    
   30.00555e-3,   /* NO+ */    
   2.01588E-3,    /* H2 */
   18.01528E-3,   /*H2O */
   1.00794E-3,    /*H */
   17.00734E-3,   /*OH*/
   33.00674E-3,   /*HO2*/
   28.01010E-3,   /*CO*/
   44.00950E-3,   /*CO2*/
   15.03452E-3,   /* CH3*/
   16.04246E-3,   /* CH4 */
   34.01468E-3,   /*H2O2*/
   29.01804E-3,   /*CHO */
   30.02598E-3,   /*CH2O*/
   31.03392E-3,   /*CH3O*/
   27.04522E-3,   /*C2H3*/
   28.05316E-3,   /*C2H4*/
   29.06110E-3,   /*C2H5*/
   30.06904E-3,   /*C2H6*/
   13.01864E-3,   //CH
   15.01464E-3,   //NH
   26.03728E-3,   //C2H2
   44.09562E-3,   //C3H8
   167.31102E-3,  //C12H23
   4.002602E-3,   //He
   28.96512E-3,   // Air  
   2.01533E-3,    /* H2+ */
   132.90545E-3,  /* Cs */
   132.90490E-3,  /* Cs+ */
   39.94800E-3,    /*Ar*/ 
   12.10170E-3,    /*C*/
   24.20340E-3,    /*C2*/
   26.01740E-3,    /*CN*/
   42.01680E-3,    /*NCO*/
   39.94745E-3,    /*Ar+*/
   12.10115E-3,    /*C+*/
   24.20285E-3,    /*C2+*/
   26.01685E-3,    /*CN+*/
   28.10045E-3,    /*CO+*/
   31.01400E-3,    /*HNO*/
   46.00550E-3,    /*NO2*/
  };
  

/* the species number of charges */
const static long ck[ns_c]=
  {
   -1, /* e- */
   0,  /* O2 */
   0,  /* N2 */
   0,  /* O */
   0,  /* N */
   +1, /* O2+ */
   +1, /* N2+ */
   +1, /* O+ */
   +1, /* N+ */
   -1, /* O2- */   
   -1, /* O- */    
   0,  /* O3 */    
   0,  /* N2(A) */    
   0,  /* NO */    
   +1, /* NO+ */    
   0,  /* H2 */
   0,  /*H2O */
   0,  /*H */
   0,  /*OH*/
   0,  /*HO2*/
   0,  /*CO*/
   0,  /*CO2*/
   0,  /* CH3*/
   0,  /* CH4 */
   0,  /*H2O2*/
   0,  /*CHO */
   0,  /*CH2O*/
   0,  /*CH3O*/
   0,  /*C2H3*/
   0,  /*C2H4*/
   0,  /*C2H5*/
   0,  /*C2H6*/
   0,  /*CH*/
   0,  /*NH*/
   0,  /*C2H2*/
   0,  /*C3H8*/
   0,  /*C12H23*/
   0,  /*He */
   0,  /* Air */
   +1, /* H2+ */
   0,  /* Cs */
   +1, /* Cs+ */
   0,  /*Ar*/
   0,  /*C*/
   0,  /*C2*/
   0,  /*CN*/
   0,  /*NCO*/
   +1, /*Ar+*/
   +1, /*C+*/
   +1, /*C2+*/
   +1, /*CN+*/
   +1, /*CO+*/
   0,  /*HNO*/
   0,  /*NO2*/
  };
  
  
  
const static double Pa[ns_c][3][11]=
  {

/* These polynomials come from 
	McBride, B. "NASA Glenn Coefficients for Calculating Thermodynamic Properties
	of Individual Species" NASA-TP-2002-211556, september 2002
  CFDWARP/model/thermo/_generic/ref/mcbride.2002.pdf
*/

/* species e-
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    

    {
      {
       +298.15e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +0.0e0,           /* a1 */
       +0.0e0,           /* a2 */
       +2.5e0,           /* a3 */
       +0.0e0,           /* a4 */
       +0.0e0,           /* a5 */
       +0.0e0,           /* a6 */
       +0.0e0,           /* a7 */
       -7.45375e+2,      /* b1 */
       -1.172081224e+1   /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +0.0e0,           /* a1 */
       +0.0e0,           /* a2 */
       +2.5e0,           /* a3 */
       +0.0e0,           /* a4 */
       +0.0e0,           /* a5 */
       +0.0e0,           /* a6 */
       +0.0e0,           /* a7 */
       -7.45375e+2,      /* b1 */
       -1.172081224e+1   /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +0.0e0,           /* a1 */
       +0.0e0,           /* a2 */
       +2.5e0,           /* a3 */
       +0.0e0,           /* a4 */
       +0.0e0,           /* a5 */
       +0.0e0,           /* a6 */
       +0.0e0,           /* a7 */
       -7.45375e+2,      /* b1 */
       -1.172081224e+1   /* b2 */        
      }
    },

    
/* species O2
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -3.425563420e+04, /* a1 */
       +4.847000970e+02, /* a2 */
       +1.119010961e+00, /* a3 */
       +4.293889240e-03, /* a4 */
       -6.836300520e-07, /* a5 */
       -2.023372700e-09, /* a6 */
       +1.039040018e-12, /* a7 */
       -3.391454870e+03, /* b1 */
       +1.849699470e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -1.037939022e+06, /* a1 */
       +2.344830282e+03, /* a2 */
       +1.819732036e+00, /* a3 */
       +1.267847582e-03, /* a4 */
       -2.188067988e-07, /* a5 */
       +2.053719572e-11, /* a6 */
       -8.193467050e-16, /* a7 */
       -1.689010929e+04, /* b1 */
       +1.738716506e+01  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +4.975294300e+08, /* a1 */
       -2.866106874e+05, /* a2 */
       +6.690352250e+01, /* a3 */
       -6.169959020e-03, /* a4 */
       +3.016396027e-07, /* a5 */
       -7.421416600e-12, /* a6 */
       +7.278175770e-17, /* a7 */
       +2.293554027e+06, /* b1 */
       -5.530621610e+02  /* b2 */        
      }
    },
    
  
/* species N2
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +2.210371497e+04, /* a1 */
       -3.818461820e+02, /* a2 */
       +6.082738360e+00, /* a3 */
       -8.530914410e-03, /* a4 */
       +1.384646189e-05, /* a5 */
       -9.625793620e-09, /* a6 */
       +2.519705809e-12, /* a7 */
       +7.108460860e+02, /* b1 */
       -1.076003744e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +5.877124060e+05, /* a1 */
       -2.239249073e+03, /* a2 */
       +6.066949220e+00, /* a3 */
       -6.139685500e-04, /* a4 */
       +1.491806679e-07, /* a5 */
       -1.923105485e-11, /* a6 */
       +1.061954386e-15, /* a7 */
       +1.283210415e+04, /* b1 */
       -1.586640027e+01  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +8.310139160e+08, /* a1 */
       -6.420733540e+05, /* a2 */
       +2.020264635e+02, /* a3 */
       -3.065092046e-02, /* a4 */
       +2.486903333e-06, /* a5 */
       -9.705954110e-11, /* a6 */
       +1.437538881e-15, /* a7 */
       +4.938707040e+06, /* b1 */
       -1.672099740e+03  /* b2 */        
      }
    },
  

/* species O
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -7.953611300e+03, /* a1 */
       +1.607177787e+02, /* a2 */
       +1.966226438e+00, /* a3 */
       +1.013670310e-03, /* a4 */
       -1.110415423e-06, /* a5 */
       +6.517507500e-10, /* a6 */
       -1.584779251e-13, /* a7 */
       +2.840362437e+04, /* b1 */
       +8.404241820e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +2.619020262e+05, /* a1 */
       -7.298722030e+02, /* a2 */
       +3.317177270e+00, /* a3 */
       -4.281334360e-04, /* a4 */
       +1.036104594e-07, /* a5 */
       -9.438304330e-12, /* a6 */
       +2.725038297e-16, /* a7 */
       +3.392428060e+04, /* b1 */
       -6.679585350e-01  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +1.779004264e+08, /* a1 */
       -1.082328257e+05, /* a2 */
       +2.810778365e+01, /* a3 */
       -2.975232262e-03, /* a4 */
       +1.854997534e-07, /* a5 */
       -5.796231540e-12, /* a6 */
       +7.191720164e-17, /* a7 */
       +8.890942630e+05, /* b1 */
       -2.181728151e+02  /* b2 */        
      }
    },
  
  
  
/* species N
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +0.000000000e+00, /* a1 */
       +0.000000000e+00, /* a2 */
       +2.500000000e+00, /* a3 */
       +0.000000000e+00, /* a4 */
       +0.000000000e+00, /* a5 */
       +0.000000000e+00, /* a6 */
       +0.000000000e+00, /* a7 */
       +5.610463780e+04, /* b1 */
       +4.193905036e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +8.876501380e+04, /* a1 */
       -1.071231500e+02, /* a2 */
       +2.362188287e+00, /* a3 */
       +2.916720081e-04, /* a4 */
       -1.729515100e-07, /* a5 */
       +4.012657880e-11, /* a6 */
       -2.677227571e-15, /* a7 */
       +5.697351330e+04, /* b1 */
       +4.865231506e+00  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +5.475181050e+08, /* a1 */
       -3.107574980e+05, /* a2 */
       +6.916782740e+01, /* a3 */
       -6.847988130e-03, /* a4 */
       +3.827572400e-07, /* a5 */
       -1.098367709e-11, /* a6 */
       +1.277986024e-16, /* a7 */
       +2.550585618e+06, /* b1 */
       -5.848769753e+02  /* b2 */        
      }
    },
  
  
  
  
/* species O2+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -8.607205450e+04, /* a1 */
       +1.051875934e+03, /* a2 */
       -5.432380470e-01, /* a3 */
       +6.571166540e-03, /* a4 */
       -3.274263750e-06, /* a5 */
       +5.940645340e-11, /* a6 */
       +3.238784790e-13, /* a7 */
       +1.345544668e+05, /* b1 */
       +2.902709750e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +7.384654880e+04, /* a1 */
       -8.459559540e+02, /* a2 */
       +4.985164160e+00, /* a3 */
       -1.611010890e-04, /* a4 */
       +6.427083990e-08, /* a5 */
       -1.504939874e-11, /* a6 */
       +1.578465409e-15, /* a7 */
       +1.446321044e+05, /* b1 */
       -5.811230650e+00  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -1.562125524e+09, /* a1 */
       +1.161406778e+06, /* a2 */
       -3.302504720e+02, /* a3 */
       +4.710937520e-02, /* a4 */
       -3.354461380e-06, /* a5 */
       +1.167968599e-10, /* a6 */
       -1.589754791e-15, /* a7 */
       -8.857866270e+06, /* b1 */
       +2.852035602e+03  /* b2 */        
      }
    },
  
  
    
/* species N2+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -3.474047470e+04, /* a1 */
       +2.696222703e+02, /* a2 */
       +3.164916370e+00, /* a3 */
       -2.132239781e-03, /* a4 */
       +6.730476400e-06, /* a5 */
       -5.637304970e-09, /* a6 */
       +1.621756000e-12, /* a7 */
       +1.790004424e+05, /* b1 */
       +6.832974166e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -2.845599002e+06, /* a1 */
       +7.058893030e+03, /* a2 */
       -2.884886385e+00, /* a3 */
       +3.068677059e-03, /* a4 */
       -4.361652310e-07, /* a5 */
       +2.102514545e-11, /* a6 */
       +5.411996470e-16, /* a7 */
       +1.340388483e+05, /* b1 */
       +5.090897022e+01  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -3.712829770e+08, /* a1 */
       +3.139287234e+05, /* a2 */
       -9.603518050e+01, /* a3 */
       +1.571193286e-02, /* a4 */
       -1.175065525e-06, /* a5 */
       +4.144441230e-11, /* a6 */
       -5.621893090e-16, /* a7 */
       -2.217361867e+06, /* b1 */
       +8.436270947e+02  /* b2 */        
      }
    },

/* species O+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +0.000000000e+00, /* a1 */
       +0.000000000e+00, /* a2 */
       +2.500000000e+00, /* a3 */
       +0.000000000e+00, /* a4 */
       +0.000000000e+00, /* a5 */
       +0.000000000e+00, /* a6 */
       +0.000000000e+00, /* a7 */
       +1.879352842e+05, /* b1 */
       +4.393376760e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -2.166513208e+05, /* a1 */
       +6.665456150e+02, /* a2 */
       +1.702064364e+00, /* a3 */
       +4.714992810e-04, /* a4 */
       -1.427131823e-07, /* a5 */
       +2.016595903e-11, /* a6 */
       -9.107157762e-16, /* a7 */
       +1.837191966e+05, /* b1 */
       +1.005690382e+01  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -2.143835383e+08, /* a1 */
       +1.469518523e+05, /* a2 */
       -3.680864540e+01, /* a3 */
       +5.036164540e-03, /* a4 */
       -3.087873854e-07, /* a5 */
       +9.186834870e-12, /* a6 */
       -1.074163268e-16, /* a7 */
       -9.614208960e+05, /* b1 */
       +3.426193080e+02  /* b2 */        
      }
    },
      
  

/* 
   species N+ 
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +5.237079210e+03, /* a1 */
       +2.299958315e+00, /* a2 */
       +2.487488821e+00, /* a3 */
       +2.737490756e-05, /* a4 */
       -3.134447576e-08, /* a5 */
       +1.850111332e-11, /* a6 */
       -4.447350984e-15, /* a7 */
       +2.256284738e+05, /* b1 */
       +5.076830786e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +2.904970374e+05, /* a1 */
       -8.557908610e+02, /* a2 */
       +3.477389290e+00, /* a3 */
       -5.288267190e-04, /* a4 */
       +1.352350307e-07, /* a5 */
       -1.389834122e-11, /* a6 */
       +5.046166279e-16, /* a7 */
       +2.310809984e+05, /* b1 */
       -1.994146545e+00  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +1.646092148e+07, /* a1 */
       -1.113165218e+04, /* a2 */
       +4.976986640e+00, /* a3 */
       -2.005393583e-04, /* a4 */
       +1.022481356e-08, /* a5 */
       -2.691430863e-13, /* a6 */
       +3.539931593e-18, /* a7 */
       +3.136284696e+05, /* b1 */
       -1.706646380e+01  /* b2 */        
      }
    },
      
  
        
/* 
   species O2-
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +1.883874344e+04, /* a1 */
       +1.149551768e+02, /* a2 */
       +1.518876821e+00, /* a3 */
       +8.016111380e-03, /* a4 */
       -9.850571030e-06, /* a5 */
       +6.044196210e-09, /* a6 */
       -1.486439845e-12, /* a7 */
       -7.101538760e+03, /* b1 */
       +1.501210380e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
       -5.655208050e+04, /* a1 */
       -2.367815862e+02, /* a2 */
       +4.675833670e+00, /* a3 */
       -2.197245300e-05, /* a4 */
       +1.711509280e-08, /* a5 */
       -1.757645062e-12, /* a6 */
       +8.248172790e-17, /* a7 */
       -5.960177750e+03, /* b1 */
       -2.436885556e+00  /* b2 */        
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -5.655208050e+04, /* a1 */
       -2.367815862e+02, /* a2 */
       +4.675833670e+00, /* a3 */
       -2.197245300e-05, /* a4 */
       +1.711509280e-08, /* a5 */
       -1.757645062e-12, /* a6 */
       +8.248172790e-17, /* a7 */
       -5.960177750e+03, /* b1 */
       -2.436885556e+00  /* b2 */        
      }
    },

/* 
   species O-
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -5.695857110e+03, /* a1 */
       +1.099287334e+02, /* a2 */
       +2.184719661e+00, /* a3 */
       +5.326359800e-04, /* a4 */
       -5.298878440e-07, /* a5 */
       +2.870216236e-10, /* a6 */
       -6.524692740e-14, /* a7 */
       +1.093287498e+04, /* b1 */
       +6.729863860e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +9.769363180e+03, /* a1 */
       +7.159604780e+00, /* a2 */
       +2.494961726e+00, /* a3 */
       +1.968240938e-06, /* a4 */
       -4.304174850e-10, /* a5 */
       +4.912083080e-14, /* a6 */
       -2.271600083e-18, /* a7 */
       +1.149554438e+04, /* b1 */
       +4.837036440e+00  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +5.662391000e+02, /* a1 */
       +7.572340320e+00, /* a2 */
       +2.498352500e+00, /* a3 */
       +1.862632395e-07, /* a4 */
       -1.151227211e-11, /* a5 */
       +3.688814210e-16, /* a6 */
       -4.793297600e-21, /* a7 */
       +1.148426000e+04, /* b1 */
       +4.813406590e+00  /* b2 */        
      }
    },
    
/* 
   species O3
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -1.282314507e+04, /* a1 */
       +5.898216640e+02, /* a2 */
       -2.547496763e+00, /* a3 */
       +2.690121526e-02, /* a4 */
       -3.528258340e-05, /* a5 */
       +2.312290922e-08, /* a6 */
       -6.044893270e-12, /* a7 */
       +1.348368701e+04, /* b1 */
       +3.852218580e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
       -3.869662480e+07, /* a1 */
       +1.023344994e+05, /* a2 */
       -8.961551600e+01, /* a3 */
       +3.706144970e-02, /* a4 */
       -4.137638740e-06, /* a5 */
       -2.725018591e-10, /* a6 */
       +5.248188110e-14, /* a7 */
       -6.517918180e+05, /* b1 */
       +7.029109520e+02  /* b2 */        
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -3.869662480e+07, /* a1 */
       +1.023344994e+05, /* a2 */
       -8.961551600e+01, /* a3 */
       +3.706144970e-02, /* a4 */
       -4.137638740e-06, /* a5 */
       -2.725018591e-10, /* a6 */
       +5.248188110e-14, /* a7 */
       -6.517918180e+05, /* b1 */
       +7.029109520e+02  /* b2 */        
      }
    },

/* 
   species N2(A)  -> set to N2 for the time being
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +2.210371497e+04, /* a1 */
       -3.818461820e+02, /* a2 */
       +6.082738360e+00, /* a3 */
       -8.530914410e-03, /* a4 */
       +1.384646189e-05, /* a5 */
       -9.625793620e-09, /* a6 */
       +2.519705809e-12, /* a7 */
       +7.108460860e+02, /* b1 */
       -1.076003744e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +5.877124060e+05, /* a1 */
       -2.239249073e+03, /* a2 */
       +6.066949220e+00, /* a3 */
       -6.139685500e-04, /* a4 */
       +1.491806679e-07, /* a5 */
       -1.923105485e-11, /* a6 */
       +1.061954386e-15, /* a7 */
       +1.283210415e+04, /* b1 */
       -1.586640027e+01  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +8.310139160e+08, /* a1 */
       -6.420733540e+05, /* a2 */
       +2.020264635e+02, /* a3 */
       -3.065092046e-02, /* a4 */
       +2.486903333e-06, /* a5 */
       -9.705954110e-11, /* a6 */
       +1.437538881e-15, /* a7 */
       +4.938707040e+06, /* b1 */
       -1.672099740e+03  /* b2 */        
      }
    },
    
/* 
   species NO  
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -1.143916503e+04, /* a1 */
       +1.536467592e+02, /* a2 */
       +3.431468730e+00, /* a3 */
       -2.668592368e-03, /* a4 */
       +8.481399120e-06, /* a5 */
       -7.685111050e-09, /* a6 */
       +2.386797655e-12, /* a7 */
       +9.098214410e+03, /* b1 */
       +6.728725490e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +2.239018716e+05, /* a1 */
       -1.289651623e+03, /* a2 */
       +5.433936030e+00, /* a3 */
       -3.656034900e-04, /* a4 */
       +9.880966450e-08, /* a5 */
       -1.416076856e-11, /* a6 */
       +9.380184620e-16, /* a7 */
       +1.750317656e+04, /* b1 */
       -8.501669090e+00  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -9.575303540e+08, /* a1 */
       +5.912434480e+05, /* a2 */
       -1.384566826e+02, /* a3 */
       +1.694339403e-02, /* a4 */
       -1.007351096e-06, /* a5 */
       +2.912584076e-11, /* a6 */
       -3.295109350e-16, /* a7 */
       -4.677501240e+06, /* b1 */
       +1.242081216e+03  /* b2 */        
      }
    },
    
/* 
   species NO+  
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +1.398106635e+03, /* a1 */
       -1.590446941e+02, /* a2 */
       +5.122895400e+00, /* a3 */
       -6.394388620e-03, /* a4 */
       +1.123918342e-05, /* a5 */
       -7.988581260e-09, /* a6 */
       +2.107383677e-12, /* a7 */
       +1.187495132e+05, /* b1 */
       -4.398433810e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +6.069876900e+05, /* a1 */
       -2.278395427e+03, /* a2 */
       +6.080324670e+00, /* a3 */
       -6.066847580e-04, /* a4 */
       +1.432002611e-07, /* a5 */
       -1.747990522e-11, /* a6 */
       +8.935014060e-16, /* a7 */
       +1.322709615e+05, /* b1 */
       -1.519880037e+01  /* b2 */        
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +2.676400347e+09, /* a1 */
       -1.832948690e+06, /* a2 */
       +5.099249390e+02, /* a3 */
       -7.113819280e-02, /* a4 */
       +5.317659880e-06, /* a5 */
       -1.963208212e-10, /* a6 */
       +2.805268230e-15, /* a7 */
       +1.443308939e+07, /* b1 */
       -4.324044462e+03 /* b2 */        
      }
    },


/* 
   species H2  
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
	4.078322810E+04,-8.009185450E+02, 8.214701670E+00,-1.269714360E-02, 1.753604930E-05,
   -1.202860160E-08, 3.368093160E-12, 2.682484380E+03,-3.043788660E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
5.608123380E+05,-8.371491340E+02, 2.975363040E+00, 1.252249930E-03,-3.740718420E-07,
    5.936628250E-11,-3.606995730E-15, 5.339815850E+03,-2.202764050E+00
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
5.608123380E+05,-8.371491340E+02, 2.975363040E+00, 1.252249930E-03,-3.740718420E-07,
    5.936628250E-11,-3.606995730E-15, 5.339815850E+03,-2.202764050E+00

      }
   },


/* 
   species H2O
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-3.947960830E+04, 5.755731020E+02, 9.317826530E-01, 7.222712860E-03,-7.342557370E-06,
    4.955043490E-09,-1.336933246E-12,-3.303974310E+04, 1.724205775E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
1.034972096E+06,-2.412698562E+03, 4.646110780E+00, 2.291998307E-03,-6.836830480E-07,
    9.426468930E-11,-4.822380530E-15,-1.384286509E+04,-7.978148510E+00
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
1.034972096E+06,-2.412698562E+03, 4.646110780E+00, 2.291998307E-03,-6.836830480E-07,
    9.426468930E-11,-4.822380530E-15,-1.384286509E+04,-7.978148510E+00
      }
    },


/* 
   species H
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
0.000000000E+00, 0.000000000E+00, 2.500000000E+00, 0.000000000E+00, 0.000000000E+00,
    0.000000000E+00, 0.000000000E+00, 2.547370801E+04,-4.466828530E-01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
6.078774250E+01,-1.819354417E-01, 2.500211817E+00,-1.226512864E-07, 3.732876330E-11,
   -5.687744560E-15, 3.410210197E-19, 2.547486398E+04,-4.481917770E-01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
6.078774250E+01,-1.819354417E-01, 2.500211817E+00,-1.226512864E-07, 3.732876330E-11,
   -5.687744560E-15, 3.410210197E-19, 2.547486398E+04,-4.481917770E-01
      }
    },

/* 
   species OH
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-1.998858990E+03, 9.300136160E+01, 3.050854229E+00, 1.529529288E-03,-3.157890998E-06,
    3.315446180E-09,-1.138762683E-12, 2.991214235E+03, 4.674110790E+00
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
1.017393379E+06,-2.509957276E+03, 5.116547860E+00, 1.305299930E-04,-8.284322260E-08,
    2.006475941E-11,-1.556993656E-15, 2.019640206E+04,-1.101282337E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
1.017393379E+06,-2.509957276E+03, 5.116547860E+00, 1.305299930E-04,-8.284322260E-08,
    2.006475941E-11,-1.556993656E-15, 2.019640206E+04,-1.101282337E+01
      }
    },


/* 
   species HO2
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-7.598882540E+04, 1.329383918E+03,-4.677388240E+00, 2.508308202E-02,-3.006551588E-05,
    1.895600056E-08,-4.828567390E-12,-5.873350960E+03, 5.193602140E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
-1.810669724E+06, 4.963192030E+03,-1.039498992E+00, 4.560148530E-03,-1.061859447E-06,
    1.144567878E-10,-4.763064160E-15,-3.200817190E+04, 4.066850920E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
-1.810669724E+06, 4.963192030E+03,-1.039498992E+00, 4.560148530E-03,-1.061859447E-06,
    1.144567878E-10,-4.763064160E-15,-3.200817190E+04, 4.066850920E+01
      }
    },



/* 
   species CO
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
1.489045326E+04,-2.922285939E+02, 5.724527170E+00,-8.176235030E-03, 1.456903469E-05,
   -1.087746302E-08, 3.027941827E-12,-1.303131878E+04,-7.859241350E+00
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
4.619197250E+05,-1.944704863E+03, 5.916714180E+00,-5.664282830E-04, 1.398814540E-07,
   -1.787680361E-11, 9.620935570E-16,-2.466261084E+03,-1.387413108E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
4.619197250E+05,-1.944704863E+03, 5.916714180E+00,-5.664282830E-04, 1.398814540E-07,
   -1.787680361E-11, 9.620935570E-16,-2.466261084E+03,-1.387413108E+01
      }
    },





/* 
   species CO2
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
4.943650540E+04,-6.264116010E+02, 5.301725240E+00, 2.503813816E-03,-2.127308728E-07,
   -7.689988780E-10, 2.849677801E-13,-4.528198460E+04,-7.048279440E+00
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
1.176962419E+05,-1.788791477E+03, 8.291523190E+00,-9.223156780E-05, 4.863676880E-09,
   -1.891053312E-12, 6.330036590E-16,-3.908350590E+04,-2.652669281E+01      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
1.176962419E+05,-1.788791477E+03, 8.291523190E+00,-9.223156780E-05, 4.863676880E-09,
   -1.891053312E-12, 6.330036590E-16,-3.908350590E+04,-2.652669281E+01
      }
    },



/* 
   species CH3
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-2.876188806E+04, 5.093268660E+02, 2.002143949E-01, 1.363605829E-02,-1.433989346E-05,
    1.013556725E-08,-3.027331936E-12, 1.408271825E+04, 2.022772791E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
2.760802663E+06,-9.336531170E+03, 1.487729606E+01,-1.439429774E-03, 2.444477951E-07,
   -2.224555778E-11, 8.395065760E-16, 7.481809480E+04,-7.919682400E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
2.760802663E+06,-9.336531170E+03, 1.487729606E+01,-1.439429774E-03, 2.444477951E-07,
   -2.224555778E-11, 8.395065760E-16, 7.481809480E+04,-7.919682400E+01
      }
    },




/* 
   species CH4
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-1.766850998E+05, 2.786181020E+03,-1.202577850E+01, 3.917619290E-02,-3.619054430E-05,
    2.026853043E-08,-4.976705490E-12,-2.331314360E+04, 8.904322750E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
3.730042760E+06,-1.383501485E+04, 2.049107091E+01,-1.961974759E-03, 4.727313040E-07,
   -3.728814690E-11, 1.623737207E-15, 7.532066910E+04,-1.219124889E+02
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
3.730042760E+06,-1.383501485E+04, 2.049107091E+01,-1.961974759E-03, 4.727313040E-07,
   -3.728814690E-11, 1.623737207E-15, 7.532066910E+04,-1.219124889E+02
      }
    },




/* 
   species H2O2
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-9.279533580E+04, 1.564748385E+03,-5.976460140E+00, 3.270744520E-02,-3.932193260E-05,
    2.509255235E-08,-6.465045290E-12,-2.494004728E+04, 5.877174180E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
1.489428027E+06,-5.170821780E+03, 1.128204970E+01,-8.042397790E-05,-1.818383769E-08,
    6.947265590E-12,-4.827831900E-16, 1.418251038E+04,-4.650855660E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
1.489428027E+06,-5.170821780E+03, 1.128204970E+01,-8.042397790E-05,-1.818383769E-08,
    6.947265590E-12,-4.827831900E-16, 1.418251038E+04,-4.650855660E+01
      }
    },


/* 
   species CHO
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-1.189851887E+04, 2.151536111E+02, 2.730224028E+00, 1.806516108E-03, 4.984300570E-06,
   -5.814567920E-09, 1.869689894E-12, 2.905755640E+03, 1.136772540E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
6.949606120E+05,-3.656223380E+03, 9.604731170E+00,-1.117129278E-03, 2.875328019E-07,
   -3.626247740E-11, 1.808329595E-15, 2.543704440E+04,-3.582473720E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
6.949606120E+05,-3.656223380E+03, 9.604731170E+00,-1.117129278E-03, 2.875328019E-07,
   -3.626247740E-11, 1.808329595E-15, 2.543704440E+04,-3.582473720E+01
      }
    },




/* 
   species CH2O
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-1.173916343E+05, 1.873628846E+03,-6.890288570E+00, 2.641561665E-02,-2.186389299E-05,
    1.005693006E-08,-2.023476949E-12,-2.307351768E+04, 6.420420550E+01      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
1.700825405E+06,-7.620853840E+03, 1.472447547E+01,-1.649111754E-03, 3.292144720E-07,
   -3.495049770E-11, 1.526135000E-15, 3.146812947E+04,-7.386478500E+01      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
1.700825405E+06,-7.620853840E+03, 1.472447547E+01,-1.649111754E-03, 3.292144720E-07,
   -3.495049770E-11, 1.526135000E-15, 3.146812947E+04,-7.386478500E+01
      }
    },



/* 
   species CH3O
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
8.657117660E+04,-6.631685250E+02, 2.257455672E+00, 2.266283789E-02,-2.970566403E-05,
    2.199341353E-08,-6.588043380E-12, 4.174102130E+03, 8.174777900E+00
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
2.101188243E+06,-8.841968800E+03, 1.822645731E+01,-1.743485034E-03, 3.340434270E-07,
   -3.430673160E-11, 1.473897771E-15, 5.309582060E+04,-9.422500590E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
2.101188243E+06,-8.841968800E+03, 1.822645731E+01,-1.743485034E-03, 3.340434270E-07,
   -3.430673160E-11, 1.473897771E-15, 5.309582060E+04,-9.422500590E+01
      }
    },



/* 
   species C2H3
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-3.347896870E+04, 1.064104103E+03,-6.403857060E+00, 3.934515480E-02,-4.760046090E-05,
    3.170071350E-08,-8.633406430E-12, 3.039122649E+04, 5.809226180E+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
2.718080093E+06,-1.030956829E+04, 1.836579807E+01,-1.580131153E-03, 2.680594939E-07,
   -2.439003999E-11, 9.209096390E-16, 9.765055590E+04,-9.760086860E+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
2.718080093E+06,-1.030956829E+04, 1.836579807E+01,-1.580131153E-03, 2.680594939E-07,
   -2.439003999E-11, 9.209096390E-16, 9.765055590E+04,-9.760086860E+01
      }
    },



/* 
   species C2H4
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-1.163605836E+05, 2.554851510E+03,-1.609746428E+01, 6.625779320E-02,-7.885081860E-05,
    5.125224820E-08,-1.370340031E-11,-6.176191070E+03, 1.093338343E+02
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
3.408763670E+06,-1.374847903E+04, 2.365898074E+01,-2.423804419E-03, 4.431395660E-07,
   -4.352683390E-11, 1.775410633E-15, 8.820429380E+04,-1.371278108E+02      
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
3.408763670E+06,-1.374847903E+04, 2.365898074E+01,-2.423804419E-03, 4.431395660E-07,
   -4.352683390E-11, 1.775410633E-15, 8.820429380E+04,-1.371278108E+02
      }
    },



/* 
   species C2H5
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-1.411312551E+05, 2.714285088E+03,-1.534977725E+01, 6.451672580E-02,-7.259143960E-05,
    4.599116010E-08,-1.218367535E-11, 5.981418840E+02, 1.090966520E+02
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
4.169220400E+06,-1.662982142E+04, 2.795442134E+01,-3.051715761E-03, 5.685160040E-07,
   -5.682863600E-11, 2.355648561E-15, 1.137010087E+05,-1.639357995E+02
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
4.169220400E+06,-1.662982142E+04, 2.795442134E+01,-3.051715761E-03, 5.685160040E-07,
   -5.682863600E-11, 2.355648561E-15, 1.137010087E+05,-1.639357995E+02
      }
    },



/* 
   species C2H6
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
-1.862044161E+05, 3.406191860E+03,-1.951705092E+01, 7.565835590E-02,-8.204173220E-05,
    5.061135800E-08,-1.319281992E-11,-2.702932890E+04, 1.298140496E+02
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
5.025782130E+06,-2.033022397E+04, 3.322552930E+01,-3.836703410E-03, 7.238405860E-07,
   -7.319182500E-11, 3.065468699E-15, 1.115963950E+05,-2.039410584E+02
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,       /* Tmax [K] */
5.025782130E+06,-2.033022397E+04, 3.322552930E+01,-3.836703410E-03, 7.238405860E-07,
   -7.319182500E-11, 3.065468699E-15, 1.115963950E+05,-2.039410584E+02
      }
    },

/* species CH  -  METHYLIDYNE
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
  Props & Hf0: TPIS,v2,pt2,1979, thermo-2.inp*/
  {
  { 200.0, 1000.0,
    2.220590133E+04, -3.405411530E+02, 5.531452290E+00, -5.794964260E-03, 7.969554880E-06,
-4.465911590E-09, 9.596338320E-13, 7.240783270E+04, -9.107673050E+00},
  { 1000.0, 5900.0,
    2.060763440E+06, -5.396206660E+03, 7.856293850E+00, -7.965907450E-04, 1.764308305E-07,
-1.976386267E-11, 5.030429510E-16,                 1.062236592E+05, -3.154757439E+01},
  { 5900.0, 6000.0,
    2.060763440E+06, -5.396206660E+03, 7.856293850E+00, -7.965907450E-04, 1.764308305E-07,
-1.976386267E-11, 5.030429510E-16,                 1.062236592E+05, -3.154757439E+01}
  },


/* species NH
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
  Props: TPIS 1978 V1 pt2 p223 DelH: JPC 1989 v93 p530, thermo-2.inp*/
  {
  { 200.0, 1000.0,
    1.359651320E+04, -1.900296604E+02, 4.518496790E+00, -2.432776899E-03, 2.377587464E-06,
-2.592797084E-10, -2.659680792E-13,                 4.280972190E+04, -3.886561616E+00 },
  { 1000.0, 5900.0,
    1.958141991E+06, -5.782861300E+03, 9.335742020E+00, -2.292910311E-03, 6.076092480E-07,
-6.647942750E-11, 2.384234783E-15,                 7.898912340E+04, -4.116970400E+01
},
  { 5900.0, 6000.0,
    1.958141991E+06, -5.782861300E+03, 9.335742020E+00, -2.292910311E-03, 6.076092480E-07,
-6.647942750E-11, 2.384234783E-15,                 7.898912340E+04, -4.116970400E+01
}
  },

/* species C2H2	acetylene    
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
   Hf:TRC(10/93) w-3040. Gurvich,1991 pt1 p47 pt2 p39 from thermo.inp */
  {
  { 200.0, 1000.0,
    1.598112089E+05,-2.216644118E+03, 1.265707813E+01,-7.979651080E-03, 8.054992750E-06,
   -2.433307673E-09,-7.529233180E-14, 3.712619060E+04,-5.244338900E+01},
  { 1000.0, 5900.0,
    1.713847410E+06,-5.929106660E+03, 1.236127943E+01, 1.314186993E-04,-1.362764431E-07,
    2.712655786E-11,-1.302066204E-15, 6.266578970E+04,-5.818960590E+01},
  { 5900.0, 6000.0,
    1.713847410E+06,-5.929106660E+03, 1.236127943E+01, 1.314186993E-04,-1.362764431E-07,
    2.712655786E-11,-1.302066204E-15, 6.266578970E+04,-5.818960590E+01}
  },


/* species C3H8 
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
  Props: TRC w1350,10/85. JPCRD V2,NO2,1973,P427, thermo-2.inp*/
  {
  {200.0, 1000.0,
    -2.433144337E+05, 4.656270810E+03, -2.939466091E+01, 1.188952745E-01, -1.376308269E-04,
 8.814823910E-08, -2.342987994E-11,                -3.540335270E+04, 1.841749277E+02
},
  {1000.0, 5900.0,
     6.420731680E+06, -2.659791134E+04, 4.534356840E+01, -5.020663920E-03, 9.471216940E-07,
-9.575405230E-11, 4.009672880E-15,                 1.455582459E+05, -2.818374734E+02
},
  {5900.0, 6000.0,
     6.420731680E+06, -2.659791134E+04, 4.534356840E+01, -5.020663920E-03, 9.471216940E-07,
-9.575405230E-11, 4.009672880E-15,                 1.455582459E+05, -2.818374734E+02
}
  },

/* species C12H23	JET-A(G) from thermo.inp
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
  */  
  {
  {273.15, 1000.0,
-6.068695590E+05, 8.328259590E+03,-4.312321270E+01, 2.572390455E-01,-2.629316040E-04,
 1.644988940E-07,-4.645335140E-11,                -7.606962760E+04, 2.794305937E+02
  },
  {1000.0, 5900.0,
 1.858356102E+07,-7.677219890E+04, 1.419826133E+02,-7.437524530E-03, 5.856202550E-07,
 1.223955647E-11,-3.149201922E-15,                 4.221989520E+05, -8.986061040E+02
   },
  {5900.0, 6000.0,
 1.858356102E+07,-7.677219890E+04, 1.419826133E+02,-7.437524530E-03, 5.856202550E-07,
 1.223955647E-11,-3.149201922E-15,                 4.221989520E+05, -8.986061040E+02
   }
  },

/* species He (Helium) Ref-Elm. Moore,1971. Moore,1970a. Gordon,1999 from thermo.inp
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
  */  
  {
  { 200.0, 1000.0,
    0.000000000E+00, 0.000000000E+00, 2.500000000E+00, 0.000000000E+00, 0.000000000E+00,
    0.000000000E+00, 0.000000000E+00,-7.453750000E+02, 9.287239740E-01
  },
  { 1000.0, 5900.0,
    0.000000000E+00, 0.000000000E+00, 2.500000000E+00, 0.000000000E+00, 0.000000000E+00,
    0.000000000E+00, 0.000000000E+00,-7.453750000E+02, 9.287239740E-01
  },
  { 5900.0, 6000.0,
    0.000000000E+00, 0.000000000E+00, 2.500000000E+00, 0.000000000E+00, 0.000000000E+00,
    0.000000000E+00, 0.000000000E+00,-7.453750000E+02, 9.287239740E-01
  }
  
  },



/* species Air
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +1.009950160e+04, /* a1 */
       -1.968275610e+02, /* a2 */
       +5.009155110e+00, /* a3 */
       -5.761013730e-03, /* a4 */
       +1.066859930e-05, /* a5 */
       -7.940297970e-09, /* a6 */
       +2.185231910e-12, /* a7 */
       -1.767967310e+02, /* b1 */
       -3.921504225e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
       +2.415214430e+05, /* a1 */
       -1.257874600e+03, /* a2 */
       +5.144558670e+00, /* a3 */
       -2.138541790e-04, /* a4 */
       +7.065227840e-08, /* a5 */
       -1.071483490e-11, /* a6 */
       +6.577800150e-16, /* a7 */
       +6.462263190e+03, /* b1 */
       -8.147411905e+00  /* b2 */
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +2.415214430e+05, /* a1 */
       -1.257874600e+03, /* a2 */
       +5.144558670e+00, /* a3 */
       -2.138541790e-04, /* a4 */
       +7.065227840e-08, /* a5 */
       -1.071483490e-11, /* a6 */
       +6.577800150e-16, /* a7 */
       +6.462263190e+03, /* b1 */
       -8.147411905e+00  /* b2 */
      }
    },

/* species H2+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.150e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -3.120886060e+04, /* a1 */
        2.304622909e+02, /* a2 */
        3.335564420e+00, /* a3 */
       -2.419056763e-03, /* a4 */
        7.006022340e-06, /* a5 */
       -5.610010660e-09, /* a6 */
        1.564169746e-12, /* a7 */
        1.774104638e+05, /* b1 */
       -8.278523760e-01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +1.672225964e+06, /* a1 */
       -6.595184990e+03, /* a2 */
       +1.279321925e+01, /* a3 */
       -5.509345260e-03, /* a4 */
       +2.030669412e-06, /* a5 */
       -3.351027480e-10, /* a6 */
       +1.946089104e-14, /* a7 */
       +2.189999548e+05, /* b1 */
       -6.792710780e+01  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,        /* Tmax [K] */
       -1.822070983e+08, /* a1 */
       +1.018196269e+05, /* a2 */
       -1.245831898e+01, /* a3 */
        1.076647496e-03, /* a4 */
       -3.932290360e-08, /* a5 */
        6.285405030e-13, /* a6 */
       -2.094721880e-18, /* a7 */
       -6.513101500e+05, /* b1 */
        1.471415370e+02  /* b2 */
      }
    },
    
/* species Cs
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,            /* Tmin [K] */
       +1000.0e0,           /* Tmax [K] */
       +5.466584070e+01,    /* a1 */
       -8.279346040e-01,    /* a2 */
       +2.504942210e+00,    /* a3 */
       -1.494620690e-05,    /* a4 */
       +2.425976774e-08,    /* a5 */
       -2.013172322e-11,    /* a6 */
       +6.704271991e-15,    /* a7 */
       +8.459321390e+03,    /* b1 */
       +6.848825772e+00     /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +6.166040900e+06, /* a1 */
       -1.896175522e+04, /* a2 */
       +2.483229903e+01, /* a3 */
       -1.251977234e-02, /* a4 */
       +3.309017390e-06, /* a5 */
       -3.354012020e-10, /* a6 */
       +9.626500908e-15, /* a7 */
       +1.285111231e+05, /* b1 */
       -1.522942188e+02  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,        /* Tmax [K] */
       -9.566231720e+08, /* a1 */
       +4.321690420e+05, /* a2 */
       -6.371801020e+01, /* a3 */
       +5.246260580e-03, /* a4 */
       -2.366560159e-07, /* a5 */
       +5.848488480e-12, /* a6 */
       -6.169370441e-17, /* a7 */
       -3.585268840e+06, /* b1 */
       +6.156618174e+02  /* b2 */
      }
    },


/* species Cs+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.150e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +0.000000000e+00, /* a1 */
       +0.000000000e+00, /* a2 */
       +2.500000000e+00, /* a3 */
       +0.000000000e+00, /* a4 */
       +0.000000000e+00, /* a5 */
       +0.000000000e+00, /* a6 */
       +0.000000000e+00, /* a7 */
       +5.438737820e+04, /* b1 */
       +6.182757992e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +0.000000000e+00, /* a1 */
       +0.000000000e+00, /* a2 */
       +2.500000000e+00, /* a3 */
       +0.000000000e+00, /* a4 */
       +0.000000000e+00, /* a5 */
       +0.000000000e+00, /* a6 */
       +0.000000000e+00, /* a7 */
       +5.438737820e+04, /* b1 */
       +6.182757992e+00  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -2.479469300e+08, /* a1 */
       +1.405115456e+05, /* a2 */
       -2.805027359e+01, /* a3 */
       +3.087928133e-03, /* a4 */
       -1.273598265e-07, /* a5 */
       -3.748818380e-13, /* a6 */
       +1.214944533e-16, /* a7 */
       -1.072498017e+06, /* b1 */
       +2.756827325e+02  /* b2 */
      }
    },
    
/* species Ar
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +0.000000000e+00, /* a1 */
       +0.000000000e+00, /* a2 */
       +2.500000000e+00, /* a3 */
       +0.000000000e+00, /* a4 */
       +0.000000000e+00, /* a5 */
       +0.000000000e+00, /* a6 */
       +0.000000000e+00, /* a7 */
       -7.453750000e+02, /* b1 */
       +4.379674910e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +2.010538475e+01, /* a1 */
       -5.992661070e-02, /* a2 */
       +2.500069401e+00, /* a3 */
       -3.992141160e-08, /* a4 */
       +1.205272140e-11, /* a5 */
       -1.819015576e-15, /* a6 */
       +1.078576636e-19, /* a7 */
       -7.449939610e+02, /* b1 */
       +4.379180110e+00  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -9.951265080e+08, /* a1 */
       +6.458887260e+05, /* a2 */
       -1.675894697e+02, /* a3 */
       +2.319933363e-02, /* a4 */
       -1.721080911e-06, /* a5 */
       +6.531938460e-11, /* a6 */
       -9.740147729e-16, /* a7 */
       -5.078300340e+06, /* b1 */
       +1.465298484e+03  /* b2 */
      }
    },
   
/* species C
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +6.495031470e+02, /* a1 */
       -9.649010860e-01, /* a2 */
       +2.504675479e+00, /* a3 */
       -1.281448025e-05, /* a4 */
       +1.980133654e-08, /* a5 */
       -1.606144025e-11, /* a6 */
       +5.314483411e-15, /* a7 */
       +8.545763110e+04, /* b1 */
       +4.747924288e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -1.289136472e+05, /* a1 */
       +1.719528572e+02, /* a2 */
       +2.646044387e+00, /* a3 */
       -3.353068950e-04, /* a4 */
       +1.742092740e-07, /* a5 */
       -2.902817829e-11, /* a6 */
       +1.642182385e-15, /* a7 */
       +8.410597850e+04, /* b1 */
       +4.130047418e+00  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +4.432528010e+08, /* a1 */
       -2.886018412e+05, /* a2 */
       +7.737108320e+01, /* a3 */
       -9.715281890e-03, /* a4 */
       +6.649595330e-07, /* a5 */
       -2.230078776e-11, /* a6 */
       +2.899388702e-16, /* a7 */
       +2.355273444e+06, /* b1 */
       -6.405123160e+02  /* b2 */
      }
    },
    
/* species C2
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +5.559634510e+05, /* a1 */
       -9.980126440e+03, /* a2 */
       +6.681620370e+01, /* a3 */
       -1.743432724e-01, /* a4 */
       +2.448523051e-04, /* a5 */
       -1.703467580e-07, /* a6 */
       +4.684527730e-11, /* a7 */
       +1.445869634e+05, /* b1 */
       -3.448229700e+02  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -9.689267930e+05, /* a1 */
       +3.561092990e+03, /* a2 */
       -5.064138930e-01, /* a3 */
       +2.945154879e-03, /* a4 */
       -7.139441190e-07, /* a5 */
       +8.670657250e-11, /* a6 */
       -4.076906810e-15, /* a7 */
       +7.681796830e+04, /* b1 */
       +3.339985240e+01  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +6.315145920e+06, /* a1 */
       +1.365420661e+04, /* a2 */
       -3.996903670e+00, /* a3 */
       +1.937561376e-03, /* a4 */
       -1.584446580e-07, /* a5 */
       +5.520861660e-12, /* a6 */
       -7.253735340e-17, /* a7 */
       +9.387024990e+03, /* b1 */
       +6.614329920e+01  /* b2 */
      }
    },
    
/* species CN
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +3.949148570e+03, /* a1 */
       -1.391590572e+02, /* a2 */
       +4.930835320e+00, /* a3 */
       -6.304670510e-03, /* a4 */
       +1.256836472e-05, /* a5 */
       -9.878300500e-09, /* a6 */
       +2.843137221e-12, /* a7 */
       +5.228455380e+04, /* b1 */
       -2.763115585e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -2.228006270e+06, /* a1 */
       +5.040733390e+03, /* a2 */
       -2.121897722e-01, /* a3 */
       +1.354901134e-03, /* a4 */
       +1.325929798e-07, /* a5 */
       -6.937006370e-11, /* a6 */
       +5.494952270e-15, /* a7 */
       +1.784496132e+04, /* b1 */
       +3.282563919e+01  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -1.794798118e+08, /* a1 */
       +1.054346069e+05, /* a2 */
       -1.729624170e+01, /* a3 */
       +2.194895530e-03, /* a4 */
       -8.508938030e-08, /* a5 */
       +9.318692990e-13, /* a6 */
       +6.358139930e-18, /* a7 */
       -7.962594120e+05, /* b1 */
       +1.913139639e+02  /* b2 */
      }
    },
    
/* species NCO
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +1.136503036e+04, /* a1 */
       -2.444613367e+02, /* a2 */
       +4.671376100e+00, /* a3 */
       +2.309387548e-03, /* a4 */
       +2.798649599e-06, /* a5 */
       -4.546357380e-09, /* a6 */
       +1.692880931e-12, /* a7 */
       +1.577649188e+04, /* b1 */
       -2.171476903e-01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
       +1.089445289e+05, /* a1 */
       -1.735459316e+03, /* a2 */
       +8.655610330e+00, /* a3 */
       -4.053229260e-04, /* a4 */
       +7.599716410e-08, /* a5 */
       -7.253804150e-12, /* a6 */
       +3.244872410e-16, /* a7 */
       +2.365792776e+04, /* b1 */
       -2.619532970e+01  /* b2 */
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +1.089445289e+05, /* a1 */
       -1.735459316e+03, /* a2 */
       +8.655610330e+00, /* a3 */
       -4.053229260e-04, /* a4 */
       +7.599716410e-08, /* a5 */
       -7.253804150e-12, /* a6 */
       +3.244872410e-16, /* a7 */
       +2.365792776e+04, /* b1 */
       -2.619532970e+01  /* b2 */
      }
    },
    
/* species Ar+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -5.731209170e+04, /* a1 */
       +7.930791470e+02, /* a2 */
       -1.717121217e+00, /* a3 */
       +1.044184018e-02, /* a4 */
       -1.180207501e-05, /* a5 */
       +6.528134780e-09, /* a6 */
       -1.447558130e-12, /* a7 */
       +1.790572230e+05, /* b1 */
       +2.949150950e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -3.835965400e+05, /* a1 */
       +8.162019700e+02, /* a2 */
       +2.301342628e+00, /* a3 */
       -4.952983770e-06, /* a4 */
       +1.205108477e-08, /* a5 */
       -2.185050286e-12, /* a6 */
       +1.265493898e-16, /* a7 */
       +1.771811455e+05, /* b1 */
       +7.947507480e+00  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +1.006884827e+07, /* a1 */
       -6.624361280e+03, /* a2 */
       +4.446908200e+00, /* a3 */
       -3.017567664e-04, /* a4 */
       +2.612882069e-08, /* a5 */
       -1.201637769e-12, /* a6 */
       +2.299206903e-17, /* a7 */
       +2.349504137e+05, /* b1 */
       -1.032262257e+01  /* b2 */
      }
    },
    
/* species C+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,         /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       +2.258535929e+03, /* a1 */
       -1.574575687e+00, /* a2 */
       +2.503637730e+00, /* a3 */
       -5.202878370e-06, /* a4 */
       +4.516908390e-09, /* a5 */
       -2.181431053e-12, /* a6 */
       +4.495047033e-16, /* a7 */
       +2.168951913e+05, /* b1 */
       +4.345699505e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -1.289136472e+05, /* a1 */
       +1.719528572e+02, /* a2 */
       +2.646044387e+00, /* a3 */
       -3.353068950e-04, /* a4 */
       +1.742092740e-07, /* a5 */
       -2.902817829e-11, /* a6 */
       +1.642182385e-15, /* a7 */
       +8.410597850e+04, /* b1 */
       +4.130047418e+00  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +5.618135320e+05, /* a1 */
       -6.047058900e+03, /* a2 */
       +5.884541470e+00, /* a3 */
       -7.211894530e-04, /* a4 */
       +6.823484110e-08, /* a5 */
       -2.599878590e-12, /* a6 */
       +3.633868358e-17, /* a7 */
       +2.581370458e+05, /* b1 */
       -2.280019759e+01  /* b2 */
      }
    },
   
/* species C2+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -9.913423840e+04, /* a1 */
       +1.347170609e+03, /* a2 */
       -3.476753160e+00, /* a3 */
       +1.676429424e-02, /* a4 */
       -1.865908025e-05, /* a5 */
       +1.091134647e-08, /* a6 */
       -2.434913818e-12, /* a7 */
       +2.335454800e+05, /* b1 */
       +4.406644620e+01  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +3.836292810e+06, /* a1 */
       -6.242062450e+03, /* a2 */
       +2.779245639e+00, /* a3 */
       +6.065865860e-03, /* a4 */
       -2.452799858e-06, /* a5 */
       +3.882942500e-10, /* a6 */
       -2.190639912e-14, /* a7 */
       +2.857447553e+05, /* b1 */
       +7.297383490e-01  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       +4.992689800e+07, /* a1 */
       -2.017309121e+04, /* a2 */
       +6.342227540e+00, /* a3 */
       +6.369922600e-04, /* a4 */
       -1.036760828e-07, /* a5 */
       +4.943272920e-12, /* a6 */
       -7.988804260e-17, /* a7 */
       +4.120857610e+05, /* b1 */
       -2.112967169e+01  /* b2 */
      }
    },
    
/* species CN+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -8.302909570e+05, /* a1 */
       +8.775687500e+03, /* a2 */
       -2.977443560e+01, /* a3 */
       +4.976897060e-02, /* a4 */
       -1.302225951e-05, /* a5 */
       -2.058325353e-08, /* a6 */
       +1.126843895e-11, /* a7 */
       +1.703860539e+05, /* b1 */
       +2.039918818e+02  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -7.153463080e+06, /* a1 */
       +1.857250421e+04, /* a2 */
       -1.084534159e+01, /* a3 */
       +6.106681430e-03, /* a4 */
       -1.191208566e-06, /* a5 */
       +1.184848778e-10, /* a6 */
       -4.799838730e-15, /* a7 */
       +9.242644960e+04, /* b1 */
       +1.135340573e+02  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -2.354919695e+08, /* a1 */
       +1.433776703e+05, /* a2 */
       -2.975360271e+01, /* a3 */
       +4.280545600e-03, /* a4 */
       -2.707260413e-07, /* a5 */
       +8.178340660e-12, /* a6 */
       -9.629506200e-17, /* a7 */
       -9.229047140e+05, /* b1 */
       +2.964624987e+02  /* b2 */
      }
    },
    
/* species CO+
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +298.15e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -2.178786658e+04, /* a1 */
       +1.288857032e+02, /* a2 */
       +3.769057550e+00, /* a3 */
       -3.431730130e-03, /* a4 */
       +8.193945750e-06, /* a5 */
       -6.463814690e-09, /* a6 */
       +1.803727574e-12, /* a7 */
       +1.482345898e+05, /* b1 */
       +3.990547070e+00  /* b2 */
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       +2.316847506e+05, /* a1 */
       -1.057646148e+03, /* a2 */
       +4.554257780e+00, /* a3 */
       +4.495520320e-04, /* a4 */
       -2.489507047e-07, /* a5 */
       +5.267566420e-11, /* a6 */
       -3.289510270e-15, /* a7 */
       +1.555050724e+05, /* b1 */
       -3.873462640e+00  /* b2 */
      },
      {
       +6000.0e0,        /* Tmin [K] */ 
       +20000.0e0,       /* Tmax [K] */
       -3.035604054e+08, /* a1 */
       +2.393118392e+05, /* a2 */
       -7.034999240e+01, /* a3 */
       +1.139551440e-02, /* a4 */
       -8.315173100e-07, /* a5 */
       +2.863705515e-11, /* a6 */
       -3.803269410e-16, /* a7 */
       -1.688617704e+06, /* b1 */
       +6.291980420e+02  /* b2 */
      }
    },

/* species HNO
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.0e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -6.854764860e+04,
       +9.551627200e+02,
       -6.000720210e-01,
        7.995176750e-03,
       -6.547079160e-07,
       -3.670513400e-09,
        1.783392519e-12,
        6.435351260e+03,
        3.048166179e+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
       -5.795614980e+06,
        1.945457427e+04,
       -2.152568374e+01,
        1.797428992e-02,
       -4.976040670e-06,
        6.397924170e-10,
       -3.142619368e-14,
       -1.104192372e+05,
        1.818650338e+02
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */ 
       +6000.0e0,        /* Tmax [K] */
       -5.795614980e+06,
        1.945457427e+04,
       -2.152568374e+01,
        1.797428992e-02,
       -4.976040670e-06,
        6.397924170e-10,
       -3.142619368e-14,
       -1.104192372e+05,
        1.818650338e+02
      }
    },

/* species NO2
   pos 0: Tmin lower range limit
   pos 1: Tmax upper range limit
   pos 2-8: a1,a2,...,a7
   pos 9-10: b1,b2
*/    
    {
      {
       +200.00e0,        /* Tmin [K] */ 
       +1000.0e0,        /* Tmax [K] */
       -5.642038780e+04,
       +9.633085720e+02,
       -2.434510974e+00,
        1.927760886e-02,
       -1.874559328e-05,
        9.145497730e-09,
       -1.777647635e-12,
       -1.547925037e+03,
        4.067851210e+01
      },
      {
       +1000.0e0,        /* Tmin [K] */ 
       +(6000.0e0-dTrangemin),        /* Tmax [K] */
        7.213001570e+05,
       -3.832615200e+03,
        1.113963285e+01,
       -2.238062246e-03,
        6.547723430e-07,
       -7.611335900e-11,
        3.328361050e-15,
        2.502497403e+04,
       -4.305130040e+01
      },
      {
       +(6000.0e0-dTrangemin),        /* Tmin [K] */
       +6000.0e0,       /* Tmax [K] */
        7.213001570e+05,
       -3.832615200e+03,
        1.113963285e+01,
       -2.238062246e-03,
        6.547723430e-07,
       -7.611335900e-11,
        3.328361050e-15,
        2.502497403e+04,
       -4.305130040e+01
      }
    }
   
  };


void find_species_name(long spec, char **name){
  *name=(char *)realloc(*name,((long)strlen(speciesname[smap[spec]])+2)*sizeof(char));
  strcpy(*name,speciesname[smap[spec]]);
}

/* same as find_species_name but replaces '-' by 'minus' and '+' by 'plus' */
void find_species_variable_name(long spec, char **name){
  long cnt;
  *name=(char *)realloc(*name,((long)strlen(speciesname[smap[spec]])+10)*sizeof(char));
  strcpy(*name,speciesname[smap[spec]]);
  cnt=0;
  do {
    if ((*name)[cnt]=='-') {
      SOAP_strcut(cnt, cnt, *name);
      SOAP_strins("minus", name,  cnt);
    }
    if ((*name)[cnt]=='+') {
      SOAP_strcut(cnt, cnt, *name);
      SOAP_strins("plus", name,  cnt);
    }
    cnt++;
  } while((*name)[cnt]!=0);
}


void read_model_thermo(void){
}

double _calM(long spec){
  double tmp;
  tmp=calM[smap[spec]];
  return(tmp);
}

double _m(long spec){
  double tmp;
  tmp=calM[smap[spec]]/calA;
  return(tmp);

}

double _C(long spec){
  double tmp;
#ifdef speceminus
  if (CHEM_NEUTRAL && spec==speceminus) {
    tmp=-echarge;
  } else {
    tmp=echarge*(double)ck[smap[spec]];
  }
#else
    tmp=echarge*(double)ck[smap[spec]];  
#endif
  return(tmp);
}

long _Charge_number(long spec){
  long tmp;
#ifdef speceminus
  if (CHEM_NEUTRAL && spec==speceminus) {
    tmp=-1;
  } else {
    tmp=ck[smap[spec]];
  }
#else
    tmp=ck[smap[spec]];
#endif

  return(tmp);
}


double _s(long spec){
  double sk;
  sk=sign((double)_C(spec));
  return(sk);  
}



static void find_Pa_index(long spec, double *T, long *index, double *dTplus, long *index2, double *weight2){
  long cnt;
  *dTplus=0.0;
  *index=0;
  for (cnt=0; cnt<=2; cnt++) {
    if (*T>=Pa[smap[spec]][cnt][0] && *T<=Pa[smap[spec]][cnt][1]) {
      *index=cnt;
    } 
  }  
  if (*T<=Pa[smap[spec]][0][0]) {
    *dTplus=*T-Pa[smap[spec]][0][0];
    *T=Pa[smap[spec]][0][0];
    *index=0;
  }
  if (*T>=Pa[smap[spec]][2][1]) {
    *dTplus=*T-Pa[smap[spec]][2][1];
    *T=Pa[smap[spec]][2][1];
    *index=2;
  }
  *index2=*index;
  *weight2=0.0;
  if (*index>0 && (*T)<(Pa[smap[spec]][*index][0]+dToverlap)) {
    *index2=*index-1;
    *weight2=0.5-0.5*((*T)-Pa[smap[spec]][*index][0])/dToverlap;
  }
  if (*index<2 && *T>Pa[smap[spec]][*index][1]-dToverlap) {
    *index2=*index+1;
    *weight2=0.5-0.5*(Pa[smap[spec]][*index][1]-(*T))/dToverlap;
  }
}

double _cpk_from_T(long spec, double T){
  double cpk,dTplus;
  long index,index2;
  double T2,T3,T4,weight2;
  
  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);
  T2=T*T;
  T3=T2*T;
  T4=T3*T;
    
  cpk=calR/calM[smap[spec]]*(Pa[smap[spec]][index][2]/T2+
                             Pa[smap[spec]][index][3]/T+
                             Pa[smap[spec]][index][4]+
                             Pa[smap[spec]][index][5]*T+
			     Pa[smap[spec]][index][6]*T2+
			     Pa[smap[spec]][index][7]*T3+
			     Pa[smap[spec]][index][8]*T4);
  if (index2!=index) {
    cpk=weight2*
        calR/calM[smap[spec]]*(Pa[smap[spec]][index2][2]/T2+
                               Pa[smap[spec]][index2][3]/T+
                               Pa[smap[spec]][index2][4]+
                               Pa[smap[spec]][index2][5]*T+
	    		       Pa[smap[spec]][index2][6]*T2+
		 	       Pa[smap[spec]][index2][7]*T3+
			       Pa[smap[spec]][index2][8]*T4)
        +(1.0-weight2)*cpk;
  }
#ifdef specN2
  if (_FLUID_N2VIBMODEL && spec==specN2) {
    cpk-=_dev_dTv_from_Tv(T+dTplus);
  }
#endif
#ifdef speceminus
  if (_FLUID_EENERGY && spec==speceminus) {
    cpk=0.0;
  }
#endif
  return(cpk);
}


double _cpk_from_T_equilibrium(long spec, double T){
  double cpk,dTplus;
  long index,index2;
  double T2,T3,T4,weight2;
  
  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);
  T2=T*T;
  T3=T2*T;
  T4=T3*T;
    
  cpk=calR/calM[smap[spec]]*(Pa[smap[spec]][index][2]/T2+
                             Pa[smap[spec]][index][3]/T+
                             Pa[smap[spec]][index][4]+
                             Pa[smap[spec]][index][5]*T+
			     Pa[smap[spec]][index][6]*T2+
			     Pa[smap[spec]][index][7]*T3+
			     Pa[smap[spec]][index][8]*T4);
  if (index2!=index) {
    cpk=weight2*
        calR/calM[smap[spec]]*(Pa[smap[spec]][index2][2]/T2+
                               Pa[smap[spec]][index2][3]/T+
                               Pa[smap[spec]][index2][4]+
                               Pa[smap[spec]][index2][5]*T+
	    		       Pa[smap[spec]][index2][6]*T2+
		 	       Pa[smap[spec]][index2][7]*T3+
			       Pa[smap[spec]][index2][8]*T4)
        +(1.0-weight2)*cpk;
  }
  return(cpk);
}



double _cp_from_w_T(spec_t w, double T){
  double cpmix;
  long spec;
  cpmix=0.0e0;
  for (spec=0; spec<ns; spec++){
    cpmix=cpmix+w[spec]*_cpk_from_T(spec,T);
  }
  return(cpmix);
}


static double _Rk(long spec){
  double Rk;
  Rk=calR/calM[smap[spec]];
  return(Rk);
}


double _hk_from_T(long spec, double T){
  double tmp,dTplus;
  long index,index2;
  double T2,T3,T4,T5,weight2;
  

  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);  
  T2=T*T;
  T3=T2*T;
  T4=T3*T;
  T5=T4*T;
  
  tmp=calR/calM[smap[spec]]*(
             -Pa[smap[spec]][index][2]/T
             +Pa[smap[spec]][index][3]*log(T)
             +Pa[smap[spec]][index][4]*T
             +Pa[smap[spec]][index][5]*T2/2.0e0+
             +Pa[smap[spec]][index][6]*T3/3.0e0+
             +Pa[smap[spec]][index][7]*T4/4.0e0+
             +Pa[smap[spec]][index][8]*T5/5.0e0+
             +Pa[smap[spec]][index][9]
	     );
	 
  if (index2!=index){
    tmp=(1.0-weight2)*tmp+weight2*calR/calM[smap[spec]]*(
             -Pa[smap[spec]][index2][2]/T
             +Pa[smap[spec]][index2][3]*log(T)
             +Pa[smap[spec]][index2][4]*T
             +Pa[smap[spec]][index2][5]*T2/2.0e0+
             +Pa[smap[spec]][index2][6]*T3/3.0e0+
             +Pa[smap[spec]][index2][7]*T4/4.0e0+
             +Pa[smap[spec]][index2][8]*T5/5.0e0+
             +Pa[smap[spec]][index2][9]
	     );
  }	     
	         
#ifdef specN2
  if (_FLUID_N2VIBMODEL && spec==specN2) {
    tmp-=_ev_from_T(T+dTplus);
  }
#endif
#ifdef speceminus
  if (_FLUID_EENERGY && spec==speceminus) {
    tmp-=_he_from_Te(T+dTplus);
  }
#endif
  if (fabs(dTplus)>1e-20){
    tmp=tmp+_cpk_from_T(spec, T)*dTplus;
  }
  return(tmp);
}


double _hk_from_T_equilibrium(long spec, double T){
  double tmp,dTplus;
  long index,index2;
  double T2,T3,T4,T5,weight2;
  

  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);  
  T2=T*T;
  T3=T2*T;
  T4=T3*T;
  T5=T4*T;
  
  tmp=calR/calM[smap[spec]]*(
             -Pa[smap[spec]][index][2]/T
             +Pa[smap[spec]][index][3]*log(T)
             +Pa[smap[spec]][index][4]*T
             +Pa[smap[spec]][index][5]*T2/2.0e0+
             +Pa[smap[spec]][index][6]*T3/3.0e0+
             +Pa[smap[spec]][index][7]*T4/4.0e0+
             +Pa[smap[spec]][index][8]*T5/5.0e0+
             +Pa[smap[spec]][index][9]
	     );
	 
  if (index2!=index){
    tmp=(1.0-weight2)*tmp+weight2*calR/calM[smap[spec]]*(
             -Pa[smap[spec]][index2][2]/T
             +Pa[smap[spec]][index2][3]*log(T)
             +Pa[smap[spec]][index2][4]*T
             +Pa[smap[spec]][index2][5]*T2/2.0e0+
             +Pa[smap[spec]][index2][6]*T3/3.0e0+
             +Pa[smap[spec]][index2][7]*T4/4.0e0+
             +Pa[smap[spec]][index2][8]*T5/5.0e0+
             +Pa[smap[spec]][index2][9]
	     );
  }	     
	         
  if (fabs(dTplus)>1e-20){
    tmp=tmp+_cpk_from_T_equilibrium(spec, T)*dTplus;
  }
  return(tmp);
}



double _h_from_w_T(spec_t w, double T){
  double hmix;
  long spec;
  hmix=0.0e0;
  for (spec=0; spec<ns; spec++){
    hmix=hmix+w[spec]*_hk_from_T(spec,T);
  }
  return(hmix);
}


/* the mixture gas constant 
   - only includes contribution from species contributing to P evaluated at T
   - thus, don't include electron contribution if _FLUID_EENERGY is defined  
*/
double _R(spec_t w){
  double Rmix;
  long spec;
  Rmix=0.0e0;
  for (spec=0; spec<ns; spec++){
    Rmix=Rmix+w[spec]*calR/calM[smap[spec]];
  }
#ifdef speceminus
  if (_FLUID_EENERGY) Rmix-=w[speceminus]*calR/calM[smap[speceminus]];
#endif
  return(Rmix);
}


double _e_from_w_T(spec_t w, double T){
  double tmp,Rmix,hmix;
  Rmix=_R(w);
  hmix=_h_from_w_T(w,T);
  tmp=hmix-Rmix*T;
  return(tmp);
}


double _rho_from_w_P_T(spec_t w, double P, double T){
  double tmp,Rmix;
  Rmix=_R(w);
  assert(T!=0.0e0);
  tmp=P/Rmix/T;
  return(tmp);
}

double _T_from_rhok_P(spec_t rhok, double P){
  double sum,tmp;
  long spec;
  sum=0.0e0;
  for (spec=0; spec<ns; spec++){
    sum=sum+rhok[spec]/calM[smap[spec]];
  }
#ifdef speceminus
  if (_FLUID_EENERGY){
    sum-=rhok[speceminus]/calM[smap[speceminus]];
  }
#endif
  assert(sum!=0.0e0);
  tmp=P/calR/sum;
  return(tmp);
}


double _T_from_w_rho_P(spec_t w, double rho, double P){
  spec_t rhok;
  double T;
  long spec;
  for (spec=0; spec<ns; spec++){
    rhok[spec]=rho*w[spec];
  }
  T=_T_from_rhok_P(rhok, P);
  return(T);
}



double _a_from_w_T(spec_t w,  double T) {
  double tmp,cpmix,Rmix;

  cpmix=_cp_from_w_T(w,T);
  Rmix=_R(w);
  assert((cpmix-Rmix)!=0.0e0);
  assert(Rmix*cpmix*T/(cpmix-Rmix)>=0.0e0);
  tmp=sqrt(Rmix*cpmix*T/(cpmix-Rmix));
  return(tmp);
}


/* here, P does not include electron pressure */
double _P_from_w_rho_T(spec_t w, double rho, double T) {
  double tmp,Rmix;
  Rmix=_R(w);
  tmp=rho*Rmix*T;
  return(tmp);
}


/* here, P includes electron pressure */
double _P_from_w_rho_T_Te(spec_t w, double rho, double T, double Te) {
  double tmp,Tk;
  long spec;
  tmp=0.0;
  for (spec=0; spec<ns; spec++){
    Tk=T;
#ifdef speceminus
    if (spec==speceminus) Tk=Te; 
#endif
    tmp=tmp+rho*w[spec]*calR/calM[smap[spec]]*Tk;
  }
  return(tmp);
}



double _T_from_w_e(spec_t w, double e) {
  double Rmix,hmix,cpmix,T,dT;
  double Phi,Phip;
  long cnt;
  bool CONVERGED;

  T=1000.0e0;
  Rmix=_R(w);
  CONVERGED=FALSE;
  cnt=0;
  do {
    cpmix=_cp_from_w_T(w,T);
    hmix=_h_from_w_T(w,T);
    Phi=e+Rmix*T-hmix;
    Phip=Rmix-cpmix;
    assert(Phip!=0.0e0);
    dT=-Phi/Phip;
    T=T+dT;
    if (fabs(dT/T)<T_accuracy) CONVERGED=TRUE;
    cnt++;
#ifndef NDEBUG
    if (cnt>95){
      wfprintf(stderr,"cnt=%ld T=%E dT=%E e=%E hmix-Rmix*T=%E Rmix=%E\n",cnt,T,dT,e,hmix-Rmix*T,Rmix);
    }
#endif
    if (cnt>100) {
      fatal_error("Couldn't find root for temperature in _T_from_w_e subroutine.\n"
                  "T=%E  dT=%E e=%E Rmix=%E \n"
                      "Exiting loop after 100 steps.\n\n"
                      ,T,dT,e,Rmix);
      //if (FALSE && finite(T) && finite(e)) CONVERGED=TRUE;
      //  else fatal_error("T or e not finite.\n");
    }
  } while (!CONVERGED);
  return(T);
}

double _T_from_w_h(spec_t w, double h) {
  double hmix,cpmix,T;
  double Phi,Phip,dT;
  long cnt;
  bool CONVERGED;

  T=300.0e0;
  cnt=0;
  CONVERGED=FALSE;
  do {
    cpmix=_cp_from_w_T(w,T);
    hmix=_h_from_w_T(w,T);
    Phi=h-hmix;
    Phip=-cpmix;
    assert(Phip!=0.0e0);
    dT=-Phi/Phip;
    T=T+dT;
    if (fabs(dT/T)<T_accuracy) CONVERGED=TRUE;
    cnt++;
    if (cnt>100) {
      fatal_error("Couldn't find root for temperature in _T_from_w_h subroutine. Exiting loop after 100 steps.\n"
                  "T=%E  dT=%E \n",T,dT);
      CONVERGED=TRUE;
    }
  } while (!CONVERGED);
  return(T);
}

double _de_dP_at_constant_rho(spec_t rhok, double T){
  double rho,cp,R,tmp;
  long spec;
  spec_t w;

  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rho=rho+rhok[spec];
  }
  for (spec=0; spec<ns; spec++){
    w[spec]=rhok[spec]/rho;
  }
  cp=_cp_from_w_T(w,T);
  R=_R(w);
  tmp=(cp-R)/rho/R;
  return(tmp);
}


double _de_dT_at_constant_rho(spec_t rhok, double T){
  double rho,cp,R,tmp;
  long spec;
  spec_t w;

  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rho=rho+rhok[spec];
  }
  for (spec=0; spec<ns; spec++){
    w[spec]=rhok[spec]/rho;
  }
  cp=_cp_from_w_T(w,T);
  R=_R(w);
  tmp=cp-R;
  return(tmp);
}


void find_de_drhok_at_constant_P(spec_t rhok, double T, spec_t dedrhok){
  double rho,cp,R,h;
  long spec;
  spec_t w;

  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rho=rho+rhok[spec];
  }
  for (spec=0; spec<ns; spec++){
    w[spec]=rhok[spec]/rho;
  }
  cp=_cp_from_w_T(w,T);
  R=_R(w);
  h=_h_from_w_T(w,T);
  for (spec=0; spec<ns; spec++){
    dedrhok[spec]=_hk_from_T(spec,T)/rho-h/rho+R*T/rho-cp*T*_Rk(spec)/rho/R;
  }
#ifdef speceminus
  if (_FLUID_EENERGY) dedrhok[speceminus]=_hk_from_T(speceminus,T)/rho-h/rho+R*T/rho;
#endif
}


void find_de_drhok_at_constant_T(spec_t rhok, double T, spec_t dedrhok){
  double rho,R,h;
  long spec;
  spec_t w;

  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rho=rho+rhok[spec];
  }
  for (spec=0; spec<ns; spec++){
    w[spec]=rhok[spec]/rho;
  }
  R=_R(w);
  h=_h_from_w_T(w,T);
  for (spec=0; spec<ns; spec++){
    dedrhok[spec]=_hk_from_T(spec,T)/rho-h/rho+R*T/rho-_Rk(spec)*T/rho;
  }
#ifdef speceminus
  if (_FLUID_EENERGY) dedrhok[speceminus]=_hk_from_T(speceminus,T)/rho-h/rho+R*T/rho;
#endif
}


double _Omega11(double T, double eps){
  double tmp;
  double Tstar;
  long index;

  assert(eps!=0.0e0);
  Tstar=T/eps;
  index=2;
  if (Tstar<10.0e0) index=1;
  if (Tstar<5.0e0) index=0;
  tmp=Pd[index][0]
     +Pd[index][1]*Tstar
     +Pd[index][2]*Tstar*Tstar
     +Pd[index][3]*Tstar*Tstar*Tstar
     +Pd[index][4]*Tstar*Tstar*Tstar*Tstar;
#ifndef NDEBUG
  if (!(tmp>0.0)){
    wfprintf(stderr,"\n  T=%E\n",T);
    wfprintf(stderr,"  index=%ld\n",index);
    wfprintf(stderr,"Pd[index][0]=%E\n",Pd[index][0]); 
    wfprintf(stderr,"Pd[index][1]=%E\n",Pd[index][1]); 
    wfprintf(stderr,"Pd[index][2]=%E\n",Pd[index][2]); 
    wfprintf(stderr,"Pd[index][3]=%E\n",Pd[index][3]); 
    wfprintf(stderr,"Pd[index][4]=%E\n",Pd[index][4]);
    fatal_error("Problem in function _Omega11() part of thermo.c. Temperature may be out of polynomial bounds."); 
  }
#endif
  return(tmp);
}


double _Omega22(double T, double eps){
  double tmp;
  double Omega11,Astar;
  double Tstar;
  long index;

  Omega11=_Omega11(T,eps);
  assert(eps!=0.0e0);
  Tstar=T/eps;
  index=2;
  if (Tstar<10.0e0) index=1;
  if (Tstar<5.0e0) index=0;
  Astar=Pe[index][0]
       +Pe[index][1]*Tstar
       +Pe[index][2]*Tstar*Tstar;
#ifndef NDEBUG
  if (Astar<0.0){
    wfprintf(stderr,"\n  T=%E\n",T);
    wfprintf(stderr,"  index=%ld\n",index);
    wfprintf(stderr,"  Pe[index][0]=%E\n",Pe[index][0]); 
    wfprintf(stderr,"  Pe[index][1]=%E\n",Pe[index][1]); 
    wfprintf(stderr,"  Pe[index][2]=%E\n",Pe[index][2]);
    fatal_error("Problem in _Omega22() part of thermo.c. Temperature may be out of polynomial bounds.");
  }
#endif
  tmp=Astar*Omega11;
  return(tmp);
}


void find_nuk_eta_kappa(spec_t w, double rho, double T, double Te,
                   spec_t nuk, double *eta, double *kappa){
  long spec,k,l;
  spec_t etak,kappak,chik;
  double P,etamix,kappamix,Nn;
  double chiden,sum;
  double calD[ns][ns];
  double nuwsum,wsum;

  // first make sure the temperature is not out of polynomial bounds
  T=min(T,MAX_T_FOR_NUK_ETA_KAPPA_POLYNOMIALS);
  P=_P_from_w_rho_T(w,rho,T);
  for (spec=0; spec<ns; spec++){
    if (w[spec]<1.0E-12) w[spec]=1.0E-9;
  }
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL) { /* don't add the charged species */
      assert((Psig[smap[spec]]*Psig[smap[spec]]*_Omega22(T,Peps[smap[spec]]))!=0.0e0);
      assert((calM[smap[spec]]*T)>=0.0e0);
      etak[spec]=8.44107E-7*sqrt(calM[smap[spec]]*T)/(Psig[smap[spec]]*Psig[smap[spec]]
                        *_Omega22(T,Peps[smap[spec]]));
#ifndef NDEBUG
      if (!(etak[spec]>0.0)){
        wfprintf(stderr,"\n  _Omega22(T,Peps[smap[spec]])=%E\n",_Omega22(T,Peps[smap[spec]])); 
        wfprintf(stderr,"  sqrt(calM[smap[spec]]*T)=%E\n",sqrt(calM[smap[spec]]*T));
        fatal_error("  Problem computing etak in find_nuk_eta_kappa().");
      }
#endif
      assert(etak[spec]>0.0);
      if (numatoms[smap[spec]]>1) {
        kappak[spec]=15.0e0/4.0e0*calR/calM[smap[spec]]*etak[spec]*
                  (0.115e0+0.354e0*calM[smap[spec]]*_cpk_from_T(spec,T)/calR);
      } else {
        kappak[spec]=15.0e0/4.0e0*calR/calM[smap[spec]]*etak[spec];
      }
    }
  }
  chiden=0.0e0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL) { /* don't add the charged species */
      chiden=chiden+w[spec]/calM[smap[spec]];
    }
  }
  assert(chiden!=0.0e0);
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL) { /* don't add the charged species */
      chik[spec]=w[spec]/calM[smap[spec]]/chiden;
    }
  }
  etamix=0.0e0;
  kappamix=0.0e0;
  for (k=0; k<ns; k++){
    if (speciestype[k]==SPECIES_NEUTRAL) { /* don't add the charged species */
      sum=0.0e0;
      for (l=0; l<ns; l++){
        if (l!=k && speciestype[l]==SPECIES_NEUTRAL) { /* don't add the charged species */
          assert(etak[l]!=0.0e0);
          assert((1.0e0+calM[smap[k]]/calM[smap[l]])*8.0e0>0.0e0);
          assert((etak[k]/etak[l])>0.0e0);
          sum=sum+chik[l]/sqrt((1.0e0+calM[smap[k]]/calM[smap[l]])*8.0e0)*
                pow(1.0e0+sqrt(etak[k]/etak[l])*pow(calM[smap[l]]/calM[smap[k]],0.25e0),2.0e0);
        }
      }	
      assert((chik[k]+sum)!=0.0e0);
      etamix=etamix+chik[k]*etak[k]/(chik[k]+sum);
      assert((chik[k]+1.0654e0*sum)!=0.0e0);
      kappamix=kappamix+chik[k]*kappak[k]/(chik[k]+1.0654e0*sum);
    }
  }
  for (k=0; k<ns; k++){
    for (l=0; l<ns; l++){
      if (speciestype[k]==SPECIES_NEUTRAL && speciestype[l]==SPECIES_NEUTRAL) { /* don't add the charged species */
        assert(Peps[smap[k]]*Peps[smap[l]]>0.0e0);
        assert((P*_Omega11(T,sqrt(Peps[smap[k]]*Peps[smap[l]]))
               *(Psig[smap[k]]+Psig[smap[l]])*(Psig[smap[k]]+Psig[smap[l]]))!=0.0e0);
        assert(T*T*T*(1.0e0/calM[smap[l]]+1.0e0/calM[smap[k]])>0.0e0);
        calD[k][l]=2.381112E-5*sqrt(T*T*T*(1.0e0/calM[smap[l]]+1.0e0/calM[smap[k]]))
                 /(P*_Omega11(T,sqrt(Peps[smap[k]]*Peps[smap[l]]))
                 *(Psig[smap[k]]+Psig[smap[l]])*(Psig[smap[k]]+Psig[smap[l]]));
      }
    }
  }
  *eta=etamix;
  *kappa=kappamix;
  Nn=0.0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL) { 
      Nn+=w[spec]*rho/_calM(spec)*calA;
    }
  }
  
  for (k=0; k<ns; k++){
    switch (speciestype[k]){
      case SPECIES_NEUTRAL:
        sum=0.0e0;
        for (l=0; l<ns; l++){
          if (l!=k && speciestype[l]==SPECIES_NEUTRAL) { /* don't add the charged species */
            assert(calD[k][l]!=0.0e0);
            sum=sum+chik[l]/calD[k][l];
          }
        }
        assert((sum+1.0E-20)!=0.0e0);
        nuk[k]=rho*(1.0e0-chik[k])/(sum+1.0E-20);
      break; 
  /* now set the diffusion coefficients of the ions
   * such will only be used when solving a quasi-neutral plasma
   * such will not be used by the fluid modules when solving the drift-diffusion model
   */
      case SPECIES_IONPLUS:
        nuk[k]=_muk_from_N_Tk_Ek(Nn, T,  0.0, k)*kB*T* rho/fabs(_C(k))*(1.0+Te/T);
      break;
      case SPECIES_IONMINUS:
        nuk[k]=_muk_from_N_Tk_Ek(Nn, T,  0.0, k)*kB*T* rho/fabs(_C(k))*(1.0+Te/T);
      break;
      case SPECIES_ELECTRON:
        nuk[k]=0.0;
      break;
      default:
        fatal_error("Problem with speciestype in find_nuk_eta_kappa().");
    }
  }
  
#ifdef speceminus
  /* now set the ambipolar diffusion coefficient of the electrons from the ions diffusion coefficient
   * such will only be used when solving a quasi-neutral plasma
   * such will not be used by the fluid modules when solving the drift-diffusion model
   */
  nuwsum=0.0;
  wsum=0.0;
  for (k=0; k<ncs; k++){
    if (speciestype[k]==SPECIES_IONPLUS || speciestype[k]==SPECIES_IONMINUS){
      nuwsum+=nuk[k]*w[k]; 
      wsum+=w[k];
    }
  }
  nuk[speceminus]=nuwsum/wsum;
#endif
}


/* electron energy */

#define melectron 9.10938291e-31
#define Relectron (kB/melectron)

double _Te_from_ee(double ee){
  double Te;
  Te=ee/(1.5e0*Relectron);
  return(Te);
}

double _ee_from_Te(double Te){
  double ee;
  ee=1.5e0*Relectron*Te;
  return(ee);
}

double _dee_dTe_from_Te(double Te){
  double tmp;
  tmp=1.5e0*Relectron;
  return(tmp);
}


double _he_from_Te(double Te){
  double he;
  he=2.5e0*Relectron*Te;
  return(he);
}



/* Nitrogen vibration energy */

double _Tv_from_ev(double ev){
  double Tv;
  if (ev>evlim) {
    Tv=Thetav/log(RN2*Thetav/ev+1.0);
  } else {
    Tv=ev/evlim*Tvlim;
  } 
  return(Tv);
}

double _ev_from_T(double T){
  double ev;

  ev=T/Tvlim*evlim;
  if (ev>evlim) {
    ev=RN2*Thetav/(exp(Thetav/T)-1.0);
  }
  return(ev);
}




double _dev_dTv_from_Tv(double T){
  double _dev_dTv_from_Tv;
  if (T<Tvlim) {
    _dev_dTv_from_Tv=1.0/Tvlim*evlim;
  } else {
    _dev_dTv_from_Tv=RN2*sqr(Thetav)*exp(Thetav/T)/sqr(T)/sqr(exp(Thetav/T)-1.0);
  }
  return(_dev_dTv_from_Tv);
}




double _sk_from_T(long spec, double T){
  double tmp,dTplus,weight2;
  long index,index2;
  
  
  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);  
    
  tmp=calR/calM[smap[spec]]*(-Pa[smap[spec]][index][2]/T/T/2.0e0
                             -Pa[smap[spec]][index][3]/T
                             +Pa[smap[spec]][index][4]*log(T)
                             +Pa[smap[spec]][index][5]*T
                             +Pa[smap[spec]][index][6]*T*T/2.0e0
                             +Pa[smap[spec]][index][7]*T*T*T/3.0e0
			     +Pa[smap[spec]][index][8]*T*T*T*T/4.0e0
			     +Pa[smap[spec]][index][10]);
  if (index2!=index){
    tmp=(1.0-weight2)*tmp+weight2*calR/calM[smap[spec]]*(-Pa[smap[spec]][index2][2]/T/T/2.0e0
                             -Pa[smap[spec]][index2][3]/T
                             +Pa[smap[spec]][index2][4]*log(T)
                             +Pa[smap[spec]][index2][5]*T
                             +Pa[smap[spec]][index2][6]*T*T/2.0e0
                             +Pa[smap[spec]][index2][7]*T*T*T/3.0e0
			     +Pa[smap[spec]][index2][8]*T*T*T*T/4.0e0
			     +Pa[smap[spec]][index2][10]);

  }


  assert(T!=0.0e0);

/*
  nitrogen vibrational energy:
  ------------------------------------------------------------------------
  first notice that by definition 
  sk=\int cpk(T)/T dT;
  also, for the nitrogen vibrational energy correction: 
  cpk(T)-=_dev_dTv_from_Tv(T);
  with
  _dev_dTv_from_Tv=RN2*sqr(Thetav)*exp(Thetav/T)/sqr(T)/sqr(exp(Thetav/T)-1.0);  
  therefore, we need to take away from tmp
  tmp-=\int RN2*sqr(Thetav)*exp(Thetav/T)/sqr(exp(Thetav/T)-1.0)/T^3 dT  
  ========================================================================
*/  
#ifdef specN2
  if (_FLUID_N2VIBMODEL && spec==specN2) {
    if (T+dTplus>Tvlim){
      tmp-=RN2*(
         (1.0+1.0/(exp(Thetav/(T+dTplus))-1.0))*Thetav/(T+dTplus)-log(exp(Thetav/(T+dTplus))-1.0)   
      );
    } 
  }
#endif

#if (defined(speceminus))
  if (_FLUID_EENERGY && spec==speceminus){
    fatal_error("In _sk_from_T() part of thermo.c, spec can't be set to speceminus when the electron temperature is not equal to the gas temperature.");
  }
#endif

  tmp+=_cpk_from_T(spec, T)*log((T+dTplus)/T);  
  return(tmp);
}


double _sk_from_T_equilibrium(long spec, double T){
  double tmp,dTplus,weight2;
  long index,index2;
  
  
  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);  
    
  tmp=calR/calM[smap[spec]]*(-Pa[smap[spec]][index][2]/T/T/2.0e0
                             -Pa[smap[spec]][index][3]/T
                             +Pa[smap[spec]][index][4]*log(T)
                             +Pa[smap[spec]][index][5]*T
                             +Pa[smap[spec]][index][6]*T*T/2.0e0
                             +Pa[smap[spec]][index][7]*T*T*T/3.0e0
			     +Pa[smap[spec]][index][8]*T*T*T*T/4.0e0
			     +Pa[smap[spec]][index][10]);
  if (index2!=index){
    tmp=(1.0-weight2)*tmp+weight2*calR/calM[smap[spec]]*(-Pa[smap[spec]][index2][2]/T/T/2.0e0
                             -Pa[smap[spec]][index2][3]/T
                             +Pa[smap[spec]][index2][4]*log(T)
                             +Pa[smap[spec]][index2][5]*T
                             +Pa[smap[spec]][index2][6]*T*T/2.0e0
                             +Pa[smap[spec]][index2][7]*T*T*T/3.0e0
			     +Pa[smap[spec]][index2][8]*T*T*T*T/4.0e0
			     +Pa[smap[spec]][index2][10]);

  }


  assert(T!=0.0e0);

  tmp+=_cpk_from_T_equilibrium(spec, T)*log((T+dTplus)/T);  
  return(tmp);
}



double _s_from_w_T(spec_t w, double T){
  double smix;
  long spec;
  smix=0.0e0;
  for (spec=0; spec<ns; spec++){
    smix=smix+w[spec]*_sk_from_T(spec,T);
  }
  return(smix);
}


double _dsk_dT_from_T(long spec, double T){
  double tmp,dTplus;
  long index,index2;
  double weight2;
  
  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);  
  tmp=calR/calM[smap[spec]]*(+Pa[smap[spec]][index][2]/T/T/T
                             +Pa[smap[spec]][index][3]/T/T
                             +Pa[smap[spec]][index][4]/T
                             +Pa[smap[spec]][index][5]
                             +Pa[smap[spec]][index][6]*T
                             +Pa[smap[spec]][index][7]*T*T
			     +Pa[smap[spec]][index][8]*T*T*T);
  if (index2!=index){
    tmp=(1.0-weight2)*tmp+weight2*calR/calM[smap[spec]]*(+Pa[smap[spec]][index2][2]/T/T/T
                             +Pa[smap[spec]][index2][3]/T/T
                             +Pa[smap[spec]][index2][4]/T
                             +Pa[smap[spec]][index2][5]
                             +Pa[smap[spec]][index2][6]*T
                             +Pa[smap[spec]][index2][7]*T*T
			     +Pa[smap[spec]][index2][8]*T*T*T);

  }			       


  if (dTplus!=0.0) tmp+=_cpk_from_T(spec, T)/(T+dTplus);
#ifdef specN2
  if (_FLUID_N2VIBMODEL && spec==specN2) {
    if (T+dTplus>Tvlim) { /* if T+dTplus<Tvlim, the following won't affect tmp */
    tmp-=RN2*sqr(Thetav)*exp(Thetav/(T+dTplus))/sqr(exp(Thetav/(T+dTplus))-1.0)/(T+dTplus)/(T+dTplus)/(T+dTplus);
    } 
  }
#endif
#if (defined(speceminus))
  if (_FLUID_EENERGY && spec==speceminus){
    fatal_error("In _dsk_dT_from_T part of thermo.c, spec can't be set to speceminus when the electron temperature is not equal to the gas temperature.");
  }
#endif

  return(tmp);
}


double _dsk_dT_from_T_equilibrium(long spec, double T){
  double tmp,dTplus;
  long index,index2;
  double weight2;
  
  find_Pa_index(spec, &T,&index,&dTplus,&index2,&weight2);  
  tmp=calR/calM[smap[spec]]*(+Pa[smap[spec]][index][2]/T/T/T
                             +Pa[smap[spec]][index][3]/T/T
                             +Pa[smap[spec]][index][4]/T
                             +Pa[smap[spec]][index][5]
                             +Pa[smap[spec]][index][6]*T
                             +Pa[smap[spec]][index][7]*T*T
			     +Pa[smap[spec]][index][8]*T*T*T);
  if (index2!=index){
    tmp=(1.0-weight2)*tmp+weight2*calR/calM[smap[spec]]*(+Pa[smap[spec]][index2][2]/T/T/T
                             +Pa[smap[spec]][index2][3]/T/T
                             +Pa[smap[spec]][index2][4]/T
                             +Pa[smap[spec]][index2][5]
                             +Pa[smap[spec]][index2][6]*T
                             +Pa[smap[spec]][index2][7]*T*T
			     +Pa[smap[spec]][index2][8]*T*T*T);

  }			       


  if (dTplus!=0.0) tmp+=_cpk_from_T_equilibrium(spec, T)/(T+dTplus);

  return(tmp);
}


/* product of electron mobility and number density as a function of electron temperature */
double _mueN_from_Te(double Te){
  double mueN;
  Te=max(Te,300.0);
  mueN=3.74E19*exp(33.5/sqrt(log(Te)));
  return(mueN);
}


double _dmueN_dTe_from_Te(double Te){
  double logTe,dmueNdTe;
  if (Te>300.0){
    logTe=log(Te);
    dmueNdTe=-6.2645E20*exp(33.5/sqrt(logTe))/(Te*pow(logTe,1.5e0));
  } else {
    dmueNdTe=0.0;
  }
  return(dmueNdTe);
}



/* find the mobility of species k [m2/Vs] using the species temperature Tk [K] and electric field in the species reference frame Ek [V/m] 

  O2+, N2+, and NO+ are found from Sinnott, G., Golden, D. E., & Varney, R. N. (1968). Positive-Ion Mobilities in Dry Air. Physical Review, 170(1), 272–275. doi:10.1103/physrev.170.272 

  O2- is found from GOSHO, Y. AND HARADA, A., “A New Technique for Measuring Negative Ion Mobilities at Atmospheric Pressure,” Journal of Physics D, Vol. 16, 1983, pp. 1159–1166.
   
  H2+, Cs+, N+, O+, O- are approximated using Fig. 8 in THE MOBILITIES OF SMALL IONS THE ATMOSPHERE AND THEIR RELATIONSHIP by E. UNGETHUM, Aerosol Science, 1974, Vol. 5, pp. 25 37. 
*/ 
double _muk_from_N_Tk_Ek(double N, double Tk, double Ek, long k){
  double mu,Estar;
  mu=0.0;
  Estar=Ek/N;
#ifdef speceminus
  if (CHEM_NEUTRAL && k==speceminus){
    /* electrons */
    mu=_mueN_from_Te(Tk)/N;
  } else {
#endif
    switch (smap[k]){
      case SMAP_eminus:
        mu=_mueN_from_Te(Tk)/N;
      break;
      case SMAP_O2plus:
        mu=1.0/N*min(1.18E23/sqrt(Tk),3.61E12/sqrt(Estar));
      break;
      case SMAP_Oplus:
        mu=1.0/N*min(1.4*1.18E23/sqrt(Tk),1.4*3.61E12/sqrt(Estar)); //based on O2+ with mass adjustment
      break;
      case SMAP_N2plus:
        mu=1.0/N*min(0.75E23/sqrt(Tk),2.03E12/sqrt(Estar));
      break;
      case SMAP_Nplus:
        mu=1.0/N*min(1.4*0.75E23/sqrt(Tk),1.4*2.03E12/sqrt(Estar)); //based on N2+ with mass adjustment
      break;
      case SMAP_O2minus:
        mu=1.0/N*min(0.97E23/sqrt(Tk),3.56E19*pow(Estar,-0.1));
      break;
      case SMAP_Ominus:
        mu=1.0/N*min(1.4*0.97E23/sqrt(Tk),1.4*3.56E19*pow(Estar,-0.1)); //based on O2- with mass adjustment
      break;
      case SMAP_NOplus:
        mu=1.0/N*min(1.62E23/sqrt(Tk),4.47E12/sqrt(Estar));
      break;
      case SMAP_H2plus:
        mu=1.0/N*min(4.0*1.00E23/sqrt(Tk),4.0*2.50E12/sqrt(Estar)); //based on Air+ with mass adjustment
      break;
      case SMAP_Csplus:
        mu=1.0/N*min(3.0/8.0*1.00E23/sqrt(Tk),3.0/8.0*2.50E12/sqrt(Estar)); //based on Air+ with mass adjustment
      break;
      case SMAP_Arplus:
        mu=1.0/N*min(0.85*1.00E23/sqrt(Tk),0.85*2.50E12/sqrt(Estar)); //based on Air+ with mass adjustment
      break;
      case SMAP_Cplus:
        mu=1.0/N*min(1.55*1.00E23/sqrt(Tk),1.55*2.50E12/sqrt(Estar)); //based on Air+ with mass adjustment
      break;
      case SMAP_C2plus:
        mu=1.0/N*min(1.10*1.00E23/sqrt(Tk),1.10*2.50E12/sqrt(Estar)); //based on Air+ with mass adjustment
      break;
      case SMAP_CNplus:
        mu=1.0/N*min(1.06*1.00E23/sqrt(Tk),1.06*2.50E12/sqrt(Estar)); //based on Air+ with mass adjustment
      break;
      case SMAP_COplus:
        mu=1.0/N*min(1.02*1.00E23/sqrt(Tk),1.02*2.50E12/sqrt(Estar)); //based on Air+ with mass adjustment
      break;
      default:
        fatal_error("Mobility can't be found for species %ld",k);
    }
#ifdef speceminus
  }
#endif
  return(mu);
}


double _dmukN_dTk_from_Tk_EkoverN(double Tk, double Ekstar, long k){
  double dmuN_dT;
  dmuN_dT=0.0;
#ifdef speceminus
  if (CHEM_NEUTRAL && k==speceminus){
    /* electrons */
    dmuN_dT=_dmueN_dTe_from_Te(Tk);
  } else {
#endif
    switch (smap[k]){
      case SMAP_eminus:
        dmuN_dT=_dmueN_dTe_from_Te(Tk);
      break;
      case SMAP_O2plus:
        if (1.18E23/sqrt(Tk)<3.61E12/sqrt(Ekstar)){
          dmuN_dT=-0.5*1.18E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_Oplus:
        if (1.18E23/sqrt(Tk)<3.61E12/sqrt(Ekstar)){
          dmuN_dT=-0.5*1.4*1.18E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_N2plus:
        if (0.75E23/sqrt(Tk)<2.03E12/sqrt(Ekstar)){
          dmuN_dT=-0.5*0.75E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_Nplus:
        if (0.75E23/sqrt(Tk)<2.03E12/sqrt(Ekstar)){
          dmuN_dT=-0.5*1.4*0.75E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_O2minus:
        if (0.97E23/sqrt(Tk)<3.56E19*pow(Ekstar,-0.1)){
          dmuN_dT=-0.5*0.97E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_Ominus:
        if (0.97E23/sqrt(Tk)<3.56E19*pow(Ekstar,-0.1)){
          dmuN_dT=-0.5*1.4*0.97E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_NOplus:
        if (1.62E23/sqrt(Tk)<4.47E12/sqrt(Ekstar)){
          dmuN_dT=-0.5*1.62E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_H2plus:
        if (4.0*1.00E23/sqrt(Tk)<4.0*2.50E12/sqrt(Ekstar)){
          dmuN_dT=-0.5*4.0*1.00E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      case SMAP_Csplus:
        if (3.0/8.0*1.00E23/sqrt(Tk)<3.0/8.0*2.50E12/sqrt(Ekstar)){
          dmuN_dT=-0.5*3.0/8.0*1.00E23*pow(Tk,-1.5);
        } else {
          dmuN_dT=0.0;
        }
      break;
      default:
        fatal_error("dmukN_dTk_from_Tk_EkoverN can't be found for species %ld",k);
    }
#ifdef speceminus
  }
#endif
  return(dmuN_dT);


}



/* fraction of electron joule heating that is consumed in nitrogen vibrational energy */ 
double _zetav_from_Te(double Te){
  double etav;
  double k[11],Tepower[11];
  long cnt;
  k[0]=+1.8115947E-3;
  k[1]=+2.1238526E-5;
  k[2]=-2.2082300E-8;
  k[3]=+7.3911515E-12;
  k[4]=-8.0418868E-16;
  k[5]=+4.3999729E-20;
  k[6]=-1.4009604E-24;
  k[7]=+2.7238062E-29;
  k[8]=-3.1981279E-34;
  k[9]=+2.0887979E-39;
  k[10]=-5.8381036E-45;
  Te=max(0.0,Te);
  Te=min(60000.0,Te);
  Tepower[0]=1.0;
  Tepower[1]=Te;
  for (cnt=2; cnt<11; cnt++) Tepower[cnt]=Tepower[cnt-1]*Te;
  etav=0.0;
  for (cnt=0; cnt<11; cnt++) etav+=k[cnt]*Tepower[cnt];  
  return(etav);
}


/* electron energy loss function as a function of electron temperature */
double _zetae_from_Te(double Te){
  double k[7];
  double xi,Tepower[7];  
  long cnt;
  
  if (Te<19444.0){
    k[0]=+5.1572656E-4;
    k[1]=+3.4153708E-8;
    k[2]=-3.2100688E-11;
    k[3]=+1.0247332E-14;
    k[4]=-1.2153348E-18;
    k[5]=+7.2206246E-23;
    k[6]=-1.4498434E-27;
  } else {
    k[0]=+2.1476152E-1;
    k[1]=-4.4507259E-5;
    k[2]=+3.5155106E-9;
    k[3]=-1.3270119E-13;
    k[4]=+2.6544932E-18;
    k[5]=-2.7145800E-23;
    k[6]=+1.1197905E-28;
  }

  Te=max(0.0,Te);
  Te=min(60000.0,Te);
  Tepower[0]=1.0;
  Tepower[1]=Te;
  for (cnt=2; cnt<7; cnt++) Tepower[cnt]=Tepower[cnt-1]*Te;
  xi=0.0;
  for (cnt=0; cnt<7; cnt++) xi+=k[cnt]*Tepower[cnt];  
  return(xi);
}



double _dzetae_dTe_from_Te(double Te){
  double k[7];
  double dzetaedTe,Tepower[7];  
  long cnt;
  
  if (Te<19444.0){
    k[0]=+5.1572656E-4;
    k[1]=+3.4153708E-8;
    k[2]=-3.2100688E-11;
    k[3]=+1.0247332E-14;
    k[4]=-1.2153348E-18;
    k[5]=+7.2206246E-23;
    k[6]=-1.4498434E-27;
  } else {
    k[0]=+2.1476152E-1;
    k[1]=-4.4507259E-5;
    k[2]=+3.5155106E-9;
    k[3]=-1.3270119E-13;
    k[4]=+2.6544932E-18;
    k[5]=-2.7145800E-23;
    k[6]=+1.1197905E-28;
  }

  Te=max(0.0,Te);
  Te=min(60000.0,Te);
  Tepower[0]=0.0;
  Tepower[1]=1.0;
  for (cnt=2; cnt<7; cnt++) Tepower[cnt]=Tepower[cnt-1]*Te;
  dzetaedTe=0.0;
  for (cnt=1; cnt<7; cnt++) dzetaedTe+=k[cnt]*Tepower[cnt]*(double)cnt;  

  return(dzetaedTe);
}



void find_Te_from_EoverN(double Estar, double *Te){
  double w[9];
  double sum,Estarexp,logEstar;
  long cnt;
  Estar=max(1e-60,Estar);
  if (Estar<3e-19){ 
    w[0] = -3.69167532692495882511E+08;
    w[1] = -6.26956713747712671757E+07;
    w[2] = -4.65528490607805550098E+06;
    w[3] = -1.97394448288739687996E+05;
    w[4] = -5.22784662897089219769E+03;
    w[5] = -8.85545617874565635930E+01;
    w[6] = -9.36914737923363882821E-01;
    w[7] = -5.66073394421067171284E-03;
    w[8] = -1.49535882691330832494E-05;
    logEstar=log(Estar);
    Estarexp=logEstar;
    sum=w[0];
    for (cnt=1; cnt<9; cnt++) {
      sum+=w[cnt]*Estarexp;
      Estarexp*=logEstar;
    }
    *Te=exp(sum);
  } else {
    /* the following is a continuation of the curve assuming that the relationship between Estar and Te is linear for Estar>3e-19)*/
    *Te=59520.0+(Estar-3e-19)*1.5e23; 
  }
}


void find_EoverN_from_Te(double Te, double *Estar){
  long cnt;
  double Te_new,dTedEstar;
  cnt=0;
  *Estar=0.5e-20;
  do {
    find_Te_from_EoverN(*Estar, &Te_new);
    find_dTe_dEoverN_from_EoverN(*Estar, &dTedEstar);
    (*Estar)-=(Te_new-Te)/dTedEstar;
    cnt++;
  } while (cnt<100 && fabs((Te_new-Te)/Te)>1e-5);
  if (cnt==100) fatal_error("Problem finding Eover in find_EoverN_from_Te. Te=%E Te_new=%E",Te,Te_new);
}



/*
static double _Res_Telocal(double EoverN, double Te){
  double theta,mueN,Res,keiN2,keiO2,Qen,Qei;
  mueN=_mueN_from_Te(Te);
  Qen=3.0*kB*(Te-300.0)*fabs(_C(speceminus))*_zetae_from_Te(Te)/(2.0*emass*mueN);
  Qei=0.0;
  theta=log(EoverN);
  keiN2=exp(-0.0105809*sqr(theta)-2.40411e-75*pow(theta,46.0))*1E-6;
  keiO2=exp(-0.0102785*sqr(theta)-2.42260e-75*pow(theta,46.0))*1E-6;
  Qei+=1.947E-18*0.2*keiO2;
  Qei+=2.507E-18*0.8*keiN2;
  Res=sqr(EoverN)*mueN*fabs(_C(speceminus))-Qen-Qei; 
  return(Res);
}


void find_Te_from_EoverN_test(double EoverN, double *Te){
  double dResTelocaldTe,dTe,Tenew;
  long cnt;
  Tenew=5000.0;
  dTe=1.0;
  cnt=0;
  do {
    *Te=Tenew;
    dResTelocaldTe=(_Res_Telocal(EoverN, *Te+dTe)-_Res_Telocal(EoverN, *Te))/dTe;
    Tenew=*Te-_Res_Telocal(EoverN, *Te)/dResTelocaldTe;
    cnt++;
//   printf("%ld %E  %E\n",cnt,EoverN,*Te);
  } while(Tenew-*Te>0.5 && cnt<100);
  if (cnt>99) {
    fatal_error("Problem finding root in find_Te_from_EoverN\n");
  }
}
*/


void find_dTe_dEoverN_from_EoverN(double Estar, double *dTedEstar){
  double w[9];
  double sum,Estarexp,logEstar,dsumdlogEstar,dlogEstardEstar;
  long cnt;
  double Te;
  Estar=max(1e-40,Estar);  //make sure Estar is not zero
  if (Estar<3e-19){ 
    w[0] = -3.69167532692495882511E+08;
    w[1] = -6.26956713747712671757E+07;
    w[2] = -4.65528490607805550098E+06;
    w[3] = -1.97394448288739687996E+05;
    w[4] = -5.22784662897089219769E+03;
    w[5] = -8.85545617874565635930E+01;
    w[6] = -9.36914737923363882821E-01;
    w[7] = -5.66073394421067171284E-03;
    w[8] = -1.49535882691330832494E-05;
    logEstar=log(Estar);
    Estarexp=logEstar;
    sum=w[0];
    dsumdlogEstar=0.0;
    for (cnt=1; cnt<9; cnt++) {
      sum+=w[cnt]*Estarexp;
      dsumdlogEstar+=w[cnt]*Estarexp/logEstar*(double)cnt;
      Estarexp*=logEstar;
    }
    dlogEstardEstar=1.0/Estar;
    sum=exp(sum);
    Te=sum;
    *dTedEstar=Te*dsumdlogEstar*dlogEstardEstar;
  } else {
    *dTedEstar=1.5e23; 
  }
}


/* ionization potential in Joule of species k
 * the ionization potential is the energy needed to liberate one electron from one molecule 
 * in the following reaction: e- + M -> e- + e- + M+
 * This can be found by subtracting hmolar of products from the hmolar of the reactants
 * hmolar for all species can be found using src/test 
 * */
double _ionizationpotential(long k){
  double ionizationenergy;
  switch (smap[k]){
    case SMAP_NO:
      ionizationenergy=1.504e-18;
    break;
    case SMAP_N2:
      ionizationenergy=2.507E-18;    
    break;
    case SMAP_O2:
      ionizationenergy=1.947E-18;
    break;
    case SMAP_H2:
      ionizationenergy=2.472e-18;
    break;
    case SMAP_Cs:
      ionizationenergy=6.445e-19;
    break;
    default:
      ionizationenergy=0.0;
      fatal_error("spec %ld does not have an ionization energy specified in _ionizationpotential().",k);
  }
  return(ionizationenergy);
}

