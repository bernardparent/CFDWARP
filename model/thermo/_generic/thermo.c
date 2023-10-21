// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Felipe Martin Rodriguez Fuentes
Copyright 1998-2018,2020 Bernard Parent
Copyright 2020 Aaron Trinh
Copyright 2001 Jason Etele
Copyright 2000 Giovanni Fusina
Copyright 2023 Ajjay Omprakas


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
#include "enthalpy.h"
#include "thermo.h"



#define T_accuracy 1.0e-12






typedef char speciesname_t[SMAP_NS];


const static speciesname_t speciesname[SMAP_NS]=
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
   "HCO",
   "HCHO",
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
   "NO2",
   "H+",
   "CH2",
   "He+",
   "C2H",
   "N2(A3Sigma)",
   "N2(B3Pi)",
   "N2(ap1Sigma)",
   "N2(C3Pi)",
   "O(1D)",
   "O(1S)",
   "C2H4+",
   "HCCO",
   "C3H6p",
   "CH2CO",
   "CH3CHO",
   "C3H5",
   "C3H7n",
   "H2CC",
   "CH2(S)",
   "CH3CO",
   "C3H4p",
   "C3H3p1",
   "NH3",
   "NH2",
   "NH4+",
   "NH3+",
   "NH2+",
   "NH4",
   "NNH",
   "N2H2",
   "N2H3",
   "N2H4",
   "H-",
   "H3+",
   "N3+",
   "N4+",
   "NH2-",
   "NH+",
   "NNH+",
   "H(2P)",
   "H2(C1Pi)",
   "N(2D)",
   "H3",
   "N4",
   "N3",
  };







/*
Molecular Weight Reference:
NASA RP-1311 volume I & II, 1994 & 1996

Units: kg/mol

*/     

const static double calM[SMAP_NS]=
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
  };


/* the species number of charges */
const static long ck[SMAP_NS]=
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
  };
  
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
  }; 
  

  

void find_species_name(long spec, char **name){
  *name=(char *)realloc(*name,((long)strlen(speciesname[smap[spec]])+2)*sizeof(char));
  strcpy(*name,speciesname[smap[spec]]);
}


int find_neutral_spec_from_ion_spec(long specion, long *specneutral){
  bool FOUND;
  long spec;
  char *ionname,*neutralname;
  ionname=(char *)malloc(sizeof(char));
  neutralname=(char *)malloc(sizeof(char));
  find_species_name(specion, &ionname);
  strrep(ionname, "+", "");
  FOUND=FALSE;
  for (spec=ncs; spec<ns; spec++){
    find_species_name(spec, &neutralname);
    if (strcmp(ionname,neutralname)==0){
       if (FOUND) fatal_error("Problem in find_neutral_spec_from_ion_spec().");
       FOUND=TRUE;
       *specneutral=spec;
    }
  }
  free(ionname);
  free(neutralname);
  return(FOUND);  
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

double _numatoms(long spec){
  double tmp;
  tmp=numatoms[smap[spec]];
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

static double _cp_from_w_T(spec_t w, double T){
  double cpmix;
  long spec;
  cpmix=0.0e0;
  for (spec=0; spec<ns; spec++){
    cpmix=cpmix+w[spec]*_cpk_from_T(spec,T);
  }
  return(cpmix);
}



double _cp_from_w_T_equilibrium(spec_t w, double T){
  double cpmix;
  long spec;
  cpmix=0.0e0;
  for (spec=0; spec<ns; spec++){
    cpmix=cpmix+w[spec]*_cpk_from_T_equilibrium(spec,T);
  }
  return(cpmix);
}


double _cp_from_w_T_neutrals_equilibrium(spec_t w, double T){
  double cpmix,wneutralsum;
  long spec;
  wneutralsum=0.0;
  for (spec=0; spec<ns; spec++) if (speciestype[spec]==SPECIES_NEUTRAL) {
    wneutralsum+=w[spec];
  }
  cpmix=0.0e0;
  for (spec=0; spec<ns; spec++) if (speciestype[spec]==SPECIES_NEUTRAL){
    cpmix+=w[spec]/wneutralsum*_cpk_from_T_equilibrium(spec,T);
  }
  return(cpmix);
}


static double _Rk(long spec){
  double Rk;
  Rk=calR/calM[smap[spec]];
  return(Rk);
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



double _a_from_w_T_equilibrium(spec_t w,  double T) {
  double tmp,cpmix,Rmix;
  long spec;

  cpmix=_cp_from_w_T_equilibrium(w,T);
  Rmix=_R(w);
  assert((cpmix-Rmix)!=0.0e0);
  if  (!(Rmix*cpmix*T/(cpmix-Rmix)>=0.0e0)){
    for (spec=0; spec<ns; spec++){
      wfprintf(stderr,"w[%ld]=%E\n",spec,w[spec]); 
    }
    wfprintf(stderr,"Rmix=%E\n",Rmix); 
    wfprintf(stderr,"cpmix=%E\n",cpmix); 
    wfprintf(stderr,"T=%E\n",T); 
  }
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
  if (ev>THERMO_EVLIM) {
    Tv=THERMO_THETAV/log(THERMO_RN2*THERMO_THETAV/ev+1.0);
  } else {
    Tv=ev/THERMO_EVLIM*THERMO_TVLIM;
  } 
  return(Tv);
}

double _ev_from_T(double T){
  double ev;

  ev=T/THERMO_TVLIM*THERMO_EVLIM;
  if (ev>THERMO_EVLIM) {
    ev=THERMO_RN2*THERMO_THETAV/(exp(THERMO_THETAV/T)-1.0);
  }
  return(ev);
}




double _dev_dTv_from_Tv(double T){
  double _dev_dTv_from_Tv;
  if (T<THERMO_TVLIM) {
    _dev_dTv_from_Tv=1.0/THERMO_TVLIM*THERMO_EVLIM;
  } else {
    _dev_dTv_from_Tv=THERMO_RN2*sqr(THERMO_THETAV)*exp(THERMO_THETAV/T)/sqr(T)/sqr(exp(THERMO_THETAV/T)-1.0);
  }
  return(_dev_dTv_from_Tv);
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





/* fraction of electron joule heating that is consumed in nitrogen vibrational energy 
 * Aleksandrov, N. L., Vysikailo, F. I., Islamov, R. S., Kochetov, I. V., Napartovich, A. P., and Pevgov, V. G., “Electron
Distribution Function in 4:1 N2-O2 Mixture,” High Temperature, Vol. 19, No. 1, 1981, pp. 17–21.
 * */ 
double _zetav_from_Te_Aleksandrov(double Te){
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
  Te=max(300.0,Te);
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

  Te=max(300.0,Te);
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

  Te=max(300.0,Te);
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


/* ionization potential in eV of species k
 * the ionization potential of species M is the energy needed to liberate one electron from one particule 
 * in the following reaction: e- + M -> e- + e- + M+
 * This can be found by subtracting hmolar of products from the hmolar of the reactants
 * hmolar and dividing by the avogadro number; hmolar for all species can be found using src/test 
 * */
double _ionizationpot(long k){
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
    case SMAP_O:
      ionizationenergy=2.181e-18;
    break;
    case SMAP_N:
      ionizationenergy=2.330e-18;
    break;
    case SMAP_NH3:
      ionizationenergy=1.634e-18;
    break;
    case SMAP_H:
      ionizationenergy=2.1789e-18;
    break;
    case SMAP_H2C1Pi:
      ionizationenergy=1.6343e-18;
    break;
    case SMAP_NH:
      ionizationenergy=2.1592e-18;
    break;
    case SMAP_H2P:
      ionizationenergy=5.4497e-19;
    break;
    default:
      ionizationenergy=0.0;
      fatal_error("spec %ld does not have an ionization energy specified in _ionizationpot().",k);
  }
  return(ionizationenergy/echarge);
}




