#ifndef _THERMO_H
#define _THERMO_H


#include <src/common.h>
#include <stdio.h>

#define calR 8.314472e0  /* J/(K mol) */
#define calA 6.02214199E23  /* particules per mole */
#define Rchem 1.987192004e0  /*cal/(K mol)*/
#define echarge 1.602176462e-19 /* C */
#define emass 9.10938188E-31  /* kg */
#define epsilon0 8.854187817E-12 /* permittivity of free space */
#define kB 1.3806488E-23   /* m2*kg/(s2*K) */


const static long numatoms[SMAP_NS]=
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
   3,  /*HCO */
   4,  /*HCHO*/
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
   3,  /*NO2*/
   1,  /* H+ */
   3,  /* CH2 */
   1,  /* He+ */
   3,  /* C2H */
   2,     /* N2(A3Sigma) */
   2,     /* N2(B3Pi) */
   2,     /* N2(ap1Sigma) */
   2,     /* N2(C3Pi) */
   1,     /* O(1D) */
   1,      /* O(1S) */
   6,  /*C2H4+ */
   4,  /* HCCO */
   9,  /* C3H6p */
   5,  /* CH2CO */
   7,  /* CH3CHO */
   8,  /* C3H5 */
   10, /* C3H7n */
   4,  /* H2CC */
   3,    /*CH2(S) */
   6,    /*CH3CO */
   7,    /*C3H4p */
   6,    /*C3H3p1 */ 
   4,    /*NH3 */ 
   3,    /*NH2 */
   5,    /*NH4+ */
   4,    /*NH3+ */
   3,    /*NH2+ */ 
   5,    /*NH4 */
   3,    /*NNH */
   4,    /*N2H2 */
   5,    /*N2H3 */
   6,    /*N2H4 */
   1,    /*H- */
   3,    /*H3+ */
   3,    /*N3+ */
   4,    /*N4+ */
   3,    /*NH2- */
   2,    /*NH+ */
   3,    /*NNH+ */
   1,    /*H(2P) */
   2,    /* H2(C1Pi) */
   1,    /* N(2D) */
   3,    /*H3 */
   4,    /*N4 */
   3,    /*N3 */
   4,     /*NH3v */
   1,    /*H(3P) */
   2,    /* H2v2 */
   2,    /* H2(B1SIGMA) */
   2,    /* H2v */
   2,    /* H2v3 */
   1,    /* Ar(4S) */
   1,    /* Ar(4P) */
   2,    /* Ar2+ */
   2,    /* Ar2star */
   2    /* Ar2 */
  }; 


static const double calM[SMAP_NS]=
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
   29.01804E-3,   /*HCO */
   30.02598E-3,   /*HCHO*/
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
   1.00739E-3,     /* H+ */
   14.02658E-3,    /* CH2 */
   4.00205e-3,    /* He+ */
   25.02934e-3,    /* C2H */   
   28.0134e-3,     /* N2(A3Sigma) */
   28.0134e-3,     /* N2(B3Pi) */
   28.0134e-3,     /* N2(ap1Sigma) */
   28.0134e-3,     /* N2(C3Pi) */
   15.9994e-3,     /* O(1D) */
   15.9994e-3,     /* O(1S) */
   28.0526E-3,     /*C2H4+*/
   41.02874E-3,    /*HCCO */
   42.07974E-3,    /*C3H6p */
   42.03668E-3,    /*CH2CO */
   44.05256E-3,    /*CH3CHO */
   41.07180E-3,    /*C3H5 */
   43.08768E-3,    /*C3H7n */
   26.03728E-3,    /*H2CC */
   14.02658E-3,    /*CH2(S) */
   43.04462E-3,    /*CH3CO */
   40.06386E-3,    /*C3H4p */
   39.05592E-3,    /*C3H3p1 */
   17.0305E-3,    /*NH3 */ 
   16.0226E-3,    /*NH2 */
   18.0379E-3,    /*NH4+ */
   17.0300E-3,    /*NH3+ */
   16.0220E-3,    /*NH2+ */
   18.0385E-3,    /*NH4 */ 
   29.0213E-3,    /*NNH */
   30.0293E-3,    /*N2H2 */
   31.0372E-3,    /*N2H3 */
   32.0452E-3,    /*N2H4 */
   1.00849E-3,    /*H- */
   3.02327E-3,    /*H3+ */
   42.0196E-3,    /*N3+ */
   56.0263E-3,    /*N4+ */
   16.0231E-3,    /*NH2- */
   15.0141E-3,    /*NH+ */
   29.0208E-3,    /*NNH+ */
   1.00794E-3,    /*H(2P) */
   2.01588E-3,    /* H2(C1Pi) */
   14.0067e-3,    /* N(2D) */
   3.02327E-3,    /*H3 */
   56.0268E-3,    /*N4 */
   42.0201E-3,    /*N3 */
   17.0305E-3,    /*NH3v */
   1.00794E-3,    /*H(3P)*/
   2.01588E-3,    /*H2v2*/
   2.01588E-3,   /*H2(B1Sigma)*/   
   2.01588E-3,   /*H2v*/
   2.01588E-3,   /*H2v3*/
   39.94800E-3,   /*Ar(4S)*/
   39.94800E-3,   /*Ar(4P)*/
   79.89545E-3,   /*Ar2+*/
   79.89600E-3,   /*Ar2star*/
   79.89600E-3   /*Ar2*/
  };

/* the species number of charges */
static const long ck[SMAP_NS]=
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
   0,  /*HCO */
   0,  /*HCHO*/
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
   1,  /* H+ */
   0,  /* CH2 */
   1,  /* He+ */
   0,  /* C2H */
   0,     /* N2(A3Sigma) */
   0,     /* N2(B3Pi) */
   0,     /* N2(ap1Sigma) */
   0,     /* N2(C3Pi) */
   0,     /* O(1D) */
   0,      /* O(1S) */
   1,  /* C2H4+ */
   0,  /* HCCO */
   0,  /* C3H6p */
   0,  /* CH2CO */
   0,  /* CH3CHO */
   0,  /* C3H5 */
   0,  /* C3H7n */
   0,  /* H2CC */
   0,    /*CH2(S) */
   0,    /*CH3CO */
   0,    /*C3H4p */
   0,    /*C3H3p1 */
   0,    /*NH3 */
   0,    /*NH2 */ 
   +1,    /*NH4+ */
   +1,    /*NH3+ */ 
   +1,    /*NH2+ */
   0,    /*NH4 */
   0,    /*NNH */
   0,    /*N2H2 */
   0,    /*N2H3 */
   0,    /*N2H4 */
   -1,    /*H- */
   +1,    /*H3+ */
   +1,    /*N3+ */
   +1,    /*N4+ */
   -1,    /*NH2- */
   +1,    /*NH+ */
   +1,    /*NNH+ */
   0,    /*H(2P) */
   0,    /* H2(C1Pi) */
   0,    /* N(2D) */
   0,    /*H3 */
   0,    /*N4 */
   0,    /*N3 */
   0,    /*NH3v */
   0,    /*H(3P) */
   0,    /* H2v2 */
   0,   /* H2(B1SIGMA) */
   0,    /* H2v */
   0,    /* H2v3 */
   0,    /* Ar(4S) */
   0,    /* Ar(4P) */
   1,    /* Ar2+ */
   0,     /* Ar2star */
   0,     /* Ar2 */
  };


/* taken from - Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
              - Reviews of Modern Physics , Vol.72, No. 2, 2000
	      - P. J. Mohr AND B. N. Taylor, CODATA Recommended Values of the Fundamental Physical
	        constants, 1998 */

void read_model_thermo(void);

static inline double _calM(long spec){
  double tmp;
  tmp=calM[smap[spec]];
  return(tmp);
}

static inline double _m(long spec){
  double tmp;
  tmp=calM[smap[spec]]/calA;
  return(tmp);
}

static inline double _C(long spec){
  double tmp;
  tmp=echarge*(double)ck[smap[spec]];
  return(tmp);
}

static inline double _numatoms(long spec){
  double tmp;
  tmp=numatoms[smap[spec]];
  return(tmp);
}

static inline long _Charge_number(long spec){
  long tmp;
  tmp=ck[smap[spec]];
  return(tmp);
}

static inline double _s(long spec){
  double sk;
  sk=sign((double)_C(spec));
  return(sk);  
}

double _cp_from_w_T_equilibrium(spec_t w, double T);

double _cp_from_w_T_neutrals_equilibrium(spec_t w, double T);

double _cpk_from_T_equilibrium(long spec, double T);

double _R(spec_t w);

double _a_from_w_T_equilibrium(spec_t w,  double T);

double _e_from_w_T(spec_t w, double T);

double _h_from_w_T(spec_t w, double T);

double _hk_from_T(long spec, double T);

double _hk_from_T_equilibrium(long spec, double T);

double _P_from_w_rho_T(spec_t w, double rho, double T);

double _P_from_w_rho_T_Te(spec_t w, double rho, double T, double Te);

double _rho_from_w_P_T(spec_t w, double P, double T);

double _s_from_w_T(spec_t w, double T);

double _sk_from_T(long spec, double T);

double _sk_from_T_equilibrium(long spec, double T);

double _T_from_w_e(spec_t w, double e);

double _T_from_w_h(spec_t w, double h);

double _T_from_rhok_P(spec_t rhok, double P);

double _T_from_w_rho_P(spec_t w, double rho, double P);



		     
double _de_dP_at_constant_rho(spec_t rhok, double T);

double _de_dT_at_constant_rho(spec_t rhok, double T);

void find_de_drhok_at_constant_P(spec_t rhok, double T, spec_t dedrhok);

void find_de_drhok_at_constant_T(spec_t rhok, double T, spec_t dedrhok);

void find_species_name(long spec, char **name);

int find_neutral_spec_from_ion_spec(long specion, long *specneutral);

void find_species_variable_name(long spec, char **name);

double _dsk_dT_from_T(long spec, double T);

double _dsk_dT_from_T_equilibrium(long spec, double T);


/* electron energy */

double _Te_from_ee(double ee);

double _ee_from_Te(double Te);

double _dee_dTe_from_Te(double Te);

double _he_from_Te(double Te);


/* nitrogen vibration energy */
double _Tv_from_ev(double ev);

double _ev_from_T(double T);

double _dev_dTv_from_Tv(double T);

double _zetav_from_Te_Rodriguez2025(double Te);

double _zetav_from_Te_Parent2024(double Te);

double _zetav_from_Te_Aleksandrov(double Te);

double _zetae_from_Te(double Te);

double _dzetae_dTe_from_Te(double Te);


double _ionizationpot(long spec);

double _EoverN_from_rhok_Te(spec_t rhok, double Te);

double _Te_from_rhok_EoverN(spec_t rhok, double EoverN);

double _EoverNk_from_Te(long spec, double Te);

double _EoverNk_from_Te_chieminus(long spec, double Te, double chieminus);


#endif /* _THERMO_H */
