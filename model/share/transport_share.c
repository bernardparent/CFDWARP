// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Bernard Parent

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

#include <model/transport/_transport.h>
#include <src/common.h>
#include <model/thermo/_thermo.h>


#define LNLAMBDA_NONE 0
#define LNLAMBDA_GUPTAYOS 1
#define LNLAMBDA_RAIZER 2
#define LNLAMBDA_NRL 3
#define LNLAMBDA LNLAMBDA_RAIZER  //Use LNLAMBDA_RAIZER. Other methods are for testing purposes only.



// find the scalar thermal conductivity of the kth charged species
double _kappac_from_rhok_Tk_muk(spec_t rhok, double T, double Te, double muk, long k){
  double kappa,cp,Tk;
  if (k>=ncs) fatal_error("_kappac_from_rhok_Tk_Ek can only be used for charged species, not for species %ld.",k);
  if (smap[k]==SMAP_eminus) {
    Tk=Te;
  } else {
    Tk=T;
  }
  cp=_cpk_from_T_equilibrium(k,Tk);
  //cp=5.0/2.0*kB/_m(k);
  kappa=cp*rhok[k]*kB*Tk*muk/fabs(_C(k));
  return(kappa);
}



// find the scalar viscosity of the kth charged species
double _etac_from_rhok_Tk_muk(spec_t rhok, double T, double Te, double muk, long k){
  double Pr,eta,kappa,cp,Tk;
  if (k>=ncs) fatal_error("_etac_from_rhok_Tk_Ek can only be used for charged species, not for species %ld.",k);
  if (smap[k]==SMAP_eminus) {
    Tk=Te;
  } else {
    Tk=T;
  }
  cp=_cpk_from_T_equilibrium(k,Tk);
  kappa=_kappac_from_rhok_Tk_muk(rhok,T,Te,muk,k);
  if (speciestype[k]==SPECIES_ELECTRON){
    Pr=0.73; 
  } else {
    Pr=0.96;
  }
  eta=Pr*kappa/cp;
  return(eta);
}



void adjust_nuk_using_mobilities_given_muk(spec_t rhok, double T, double Te, chargedspec_t muk, spec_t nuk){
  long spec,k;
  spec_t w;
  double rho;
  double nuwsum,wsum;
  

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
  
  
  for (k=0; k<ncs; k++){
    switch (speciestype[k]){
      case SPECIES_IONPLUS:
        nuk[k]=muk[k]*kB*T* rho/fabs(_C(k))*(1.0+Te/T);
      break;
      case SPECIES_IONMINUS:
        nuk[k]=muk[k]*kB*T* rho/fabs(_C(k))*(1.0+Te/T);
      break;
    }
  }
  
#ifdef speceminus
  nuwsum=0.0;
  wsum=0.0;
  for (k=0; k<ncs; k++){
    if (speciestype[k]==SPECIES_IONPLUS || speciestype[k]==SPECIES_IONMINUS){
      nuwsum+=nuk[k]*max(1.0e-30,w[k]); 
      wsum+=max(1.0e-30,w[k]);
    }
  }
  nuk[speceminus]=nuwsum/(1.0e-99+wsum);
#endif

}




void adjust_muk_for_Ek_effect(long k, double Ek, double N, double *muk){
  double B,p;
  if (speciestype[k]==SPECIES_IONPLUS || speciestype[k]==SPECIES_IONMINUS){
    switch (smap[k]){
      case SMAP_O2plus:  
        B=3.61E12;
        p=-0.5;
      break;
      case SMAP_N2plus:
        B=2.03E12;
        p=-0.5;
      break;
      case SMAP_O2minus:
        B=3.56e19;
        p=-0.1;
      break;
      case SMAP_Ominus:
        B=1.4*3.56e19;
        p=-0.1;
      break;
      case SMAP_NOplus:
        B=4.47e12;
        p=-0.5;
      break;
      default:
        B=0.55/sqrt(_m(k));
        p=-0.5;
    }
    *muk=min(*muk,B*pow(Ek/N,p)/N);
  }
}


static double _lnLambda(spec_t rhok, double Te){
#ifdef speceminus  
  double Pe,rhoe,lnlambda;
  rhoe=rhok[speceminus];
  Pe = rhoe*calR*Te/_calM(speceminus)/101325.0e0; //electron pressure, atm
  /*electron pressure correction for the collision integrals of ionic species*/ /*Eq(24b)*/
  switch (LNLAMBDA){
    case LNLAMBDA_GUPTAYOS:
      lnlambda=0.5*log(2.09E-2*1E-12*Te*Te*Te*Te/Pe+1.52*pow(Te*Te*Te*Te*1E-12/Pe,2.0/3.0));
    break;
    case LNLAMBDA_RAIZER:
      // Raizer Gas Discharge Physics, page 14
      lnlambda=13.57+1.5*log(Te/11600.0)/log(10.0)-0.5*log(rhoe/_m(speceminus)/1e6)/log(10.0);
    break;
    case LNLAMBDA_NRL:
      // NRL Plasma Formulary, page 34
      lnlambda=23.0-log(sqrt(rhoe/_m(speceminus)/1e6)*pow(Te/11600.0,-1.5));
    break;
    case LNLAMBDA_NONE:
      lnlambda=1.0;
    break;
    default:
      fatal_error("LNLAMBDA set to invalid value.");
  }
#else
  lnlambda=1.0;
#endif
  lnlambda=max(1.0,lnlambda);
  return(lnlambda);
}




/* electron mobility as a function of electron temperature and mass densities*/
static double _mue_from_Nn_Ni_Te_rhok(double Nn, double Ni, double Te, spec_t rhok){
  double muen,mue,muei;
  Te=max(Te,300.0);
  muen=3.74E19*exp(33.5/sqrt(log(Te)))/Nn;
  muei=9.5e16*pow(Te,1.5)/Ni/_lnLambda(rhok,Te);
  //muei=2.5e16*pow(Te,1.5)/Ni;
  mue=1.0/(1.0/muen+1.0/muei);
  return(mue);
}


/* the weakly-ionized ion mobility corresponds to mui=min(A*Ti^n,B*Estar^p)/Nn 
 * mi is the mass of the ion in kg*/
static double _mui_from_Nn_Ni_Ti(double Nn, double Ni, double A, double Ti, double n,  double mi){
  double muin,muii,mui;
  muin=(A*pow(Ti,n))/Nn;
  muii=14.3e0/sqrt(mi)*pow(Ti,1.5)/Ni;
  mui=1.0/(1.0/muin+1.0/muii);
  return(mui);
}



/* find the mobility of species k [m2/Vs] using the species temperature Tk [K] and electric field in the species reference frame Ek [V/m] 

  O2+, N2+, and NO+ are found from Sinnott, G., Golden, D. E., & Varney, R. N. (1968). Positive-Ion Mobilities in Dry Air. Physical Review, 170(1), 272–275. doi:10.1103/physrev.170.272 

  O2- is found from GOSHO, Y. AND HARADA, A., “A New Technique for Measuring Negative Ion Mobilities at Atmospheric Pressure,” Journal of Physics D, Vol. 16, 1983, pp. 1159–1166.
   
  H2+, Cs+, N+, O+, O- are approximated using Fig. 8 in THE MOBILITIES OF SMALL IONS THE ATMOSPHERE AND THEIR RELATIONSHIP by E. UNGETHUM, Aerosol Science, 1974, Vol. 5, pp. 25 37. 
*/ 
double _muk_from_rhok_T_Te_ParentMacheret(spec_t rhok, double T, double Te, long k){
  double mu,N,Nn,Ni;
  long spec;
  mu=0.0;
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  Nn=0.0;
  for (spec=ncs; spec<ns; spec++){
    Nn+=rhok[spec]/_m(spec);
  }
  Nn=max(1e0,Nn);
  Ni=0.0;
  for (spec=0; spec<ncs; spec++){
    if (speciestype[spec]==SPECIES_IONPLUS) Ni+=rhok[spec]/_m(spec);
  }
  Ni=max(1e0,Ni);

#ifdef speceminus
  if (CHEM_NEUTRAL && k==speceminus){
    /* electrons */
    mu=_mue_from_Nn_Ni_Te_rhok(Nn,Ni,Te,rhok);
  } else {
#endif
    switch (smap[k]){
      case SMAP_eminus:
        mu=_mue_from_Nn_Ni_Te_rhok(Nn,Ni,Te,rhok);
      break;
      case SMAP_O2plus:  
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 1.18e23, T, -0.5,   _m(k));
      break;
      case SMAP_N2plus:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 0.75e23, T, -0.5,   _m(k));
      break;
      case SMAP_O2minus:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 0.97e23, T, -0.5,   _m(k));
      break;
      case SMAP_NOplus:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 1.62e23, T, -0.5,   _m(k));
      break;
      default:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 2.2e10/sqrt(_m(k)), T, -0.5,  _m(k));
    }
#ifdef speceminus
  }
#endif
  return(mu);
}




static void find_dmue_from_Nn_Ni_Te(double Nn, double Ni, double Te, double *dmuedTe, spec_t dmuedrhok){
  double logTe,dmuendTe;
  double muen,kmuei,muei,dmuendNn,dmueidTe,dmueidNi,dmuedmuen,dmuedmuei;
  spec_t dmuendrhok,dmueidrhok;
  long spec;
  if (Te>300.0){
//  muei=1e16*pow(Te,1.5)/Ni;
    muen=3.74E19*exp(33.5/sqrt(log(Te)))/Nn;
    kmuei=1.9e16;
    muei=kmuei/Ni*pow(Te,1.5);
//    mue=1.0/(1.0/muen+1.0/muei);

    
    logTe=log(Te);
    dmuendTe=-6.2645E20*exp(33.5/sqrt(logTe))/(Te*pow(logTe,1.5e0))/Nn;
    dmuendNn=-3.74E19*exp(33.5/sqrt(logTe))/sqr(Nn);
    for (spec=0; spec<ns; spec++){
      dmuendrhok[spec]=0.0;
      if (spec>=ncs) dmuendrhok[spec]=dmuendNn*1.0/_m(spec);
    }

    dmueidTe=kmuei/Ni*1.5*pow(Te,0.5);
    dmueidNi=-kmuei/sqr(Ni)*pow(Te,1.5);
    for (spec=0; spec<ns; spec++){
      dmueidrhok[spec]=0.0;
      if (speciestype[spec]==SPECIES_IONPLUS) dmueidrhok[spec]=dmueidNi*1.0/_m(spec);
    }
    
 

    dmuedmuen=1.0/sqr(muen)/sqr(1.0/muen+1.0/muei);
    dmuedmuei=1.0/sqr(muei)/sqr(1.0/muen+1.0/muei);
    *dmuedTe=dmuendTe*dmuedmuen+dmueidTe*dmuedmuei;
    for (spec=0; spec<ns; spec++) dmuedrhok[spec]=dmuendrhok[spec]*dmuedmuen+dmueidrhok[spec]*dmuedmuei;
  } else {
    *dmuedTe=0.0;
    for (spec=0; spec<ns; spec++) dmuedrhok[spec]=0.0;
  }
}



/* the mobility corresponds to mui=min(A*Ti^n,B*Estar^p)/N */
static void find_dmui_from_Nn_Ni_Ti(double Nn, double Ni, double A, double Ti, double n, double B, double Estar, double p, double mi,  double *dmuidTi, spec_t dmuidrhok){
  long spec;
  double term1,term2;
  double dmuindTi,dmuiidNi,dmuidmuin,dmuidmuii,kmuii,muii,muin,dmuiidTi;
  spec_t dmuindrhok,dmuiidrhok;
  
  for (spec=0; spec<ns; spec++) dmuindrhok[spec]=0.0;

  muin=min(A*pow(Ti,n),B*pow(Estar,p))/Nn;
  term1=A*pow(Ti,n);
  term2=B*pow(Estar,p);
  if (term1<term2){
  	dmuindTi=term1*n/Ti/Nn;
    for (spec=ncs; spec<ns; spec++) dmuindrhok[spec]=term1;
  } else {
  	dmuindTi=0.0;
    for (spec=ncs; spec<ns; spec++) dmuindrhok[spec]=term2;
  }
  for (spec=ncs; spec<ns; spec++) dmuindrhok[spec]*=-1.0/(Nn*Nn)/_m(spec);

  kmuii=14.3/sqrt(mi);
  muii=kmuii/Ni*pow(Ti,1.5);
  dmuiidTi=kmuii/Ni*1.5*pow(Ti,0.5);
  dmuiidNi=-kmuii/sqr(Ni)*pow(Ti,1.5);
  for (spec=0; spec<ns; spec++){
    dmuiidrhok[spec]=0.0;
    if (speciestype[spec]==SPECIES_IONPLUS) dmuiidrhok[spec]=dmuiidNi*1.0/_m(spec);
  }


  dmuidmuin=1.0/sqr(muin)/sqr(1.0/muin+1.0/muii);
  dmuidmuii=1.0/sqr(muii)/sqr(1.0/muin+1.0/muii);
  *dmuidTi=dmuindTi*dmuidmuin+dmuiidTi*dmuidmuii;
  for (spec=0; spec<ns; spec++) dmuidrhok[spec]=dmuindrhok[spec]*dmuidmuin+dmuiidrhok[spec]*dmuidmuii;

}


void find_dmuk_from_rhok_Tk_Ek_ParentMacheret(spec_t rhok, double Tk, double Ek, long k, double *dmukdTk, spec_t dmukdrhok){
  double N,Nn,Ni,Ekstar;
  long spec;
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  Nn=0.0;
  for (spec=ncs; spec<ns; spec++){
    Nn+=rhok[spec]/_m(spec);
  }
  Nn=max(1e0,Nn);
  Ni=0.0;
  for (spec=0; spec<ncs; spec++){
    if (speciestype[spec]==SPECIES_IONPLUS) Ni+=rhok[spec]/_m(spec);
  }
  Ni=max(1e0,Ni);
  Ekstar=Ek/N;
#ifdef speceminus
  if (CHEM_NEUTRAL && k==speceminus){
    /* electrons */
    find_dmue_from_Nn_Ni_Te(Nn, Ni, Tk, dmukdTk, dmukdrhok);
  } else {
#endif
    for (spec=0; spec<ns; spec++) dmukdrhok[spec]=0.0;
    switch (smap[k]){
      case SMAP_eminus:
        find_dmue_from_Nn_Ni_Te(Nn, Ni, Tk, dmukdTk, dmukdrhok);
      break;
      case SMAP_O2plus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 1.18e23, Tk, -0.5, 3.61E12, Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_N2plus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 0.75e23, Tk, -0.5, 2.03E12, Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_O2minus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 0.97e23, Tk, -0.5, 3.56e19, Ekstar, -0.1, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_Ominus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 1.4*0.97e23, Tk, -0.5, 1.4*3.56e19, Ekstar, -0.1, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_NOplus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 1.62e23, Tk, -0.5, 4.47e12, Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
      break;
      default:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 2.2e10/sqrt(_m(k)), Tk, -0.5, 0.55/sqrt(_m(k)), Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
    }
#ifdef speceminus
  }
#endif


}

