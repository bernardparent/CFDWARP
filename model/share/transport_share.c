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
#define LNLAMBDA LNLAMBDA_NRL  //Use LNLAMBDA_NRL. Other methods are for testing purposes only.

#define MASSDIFFCOEFF_ORANBORIS 1
#define MASSDIFFCOEFF_SCEBD 2
#define MASSDIFFCOEFF MASSDIFFCOEFF_ORANBORIS  //Use MASSDIFFCOEFF_ORANBORIS. 


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


void find_muk_from_nuk(spec_t nuk, spec_t rhok, double T, double Te, chargedspec_t muk){
  long spec,k;
  double rho;
  

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
  

  for (k=0; k<ncs; k++){
    switch (speciestype[k]){
      case SPECIES_IONPLUS:
        muk[k]=nuk[k]/(kB*T*rho)*fabs(_C(k));
      break;
      case SPECIES_IONMINUS:
        muk[k]=nuk[k]/(kB*T*rho)*fabs(_C(k));
      break;
      case SPECIES_ELECTRON:
        muk[k]=nuk[k]/(kB*Te*rho)*fabs(_C(k));
      break;
    }
  }
}


/* compute nue given nui for a quasi-neutral plasma */
void adjust_nue_given_nui(spec_t rhok, double T, double Te, spec_t nuk){
  long spec;
  double rho;
  

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];

#ifdef speceminus
  double nuwsum,wsum;
  long k;
  spec_t w;

  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
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


void adjust_nuk_using_mobilities_given_muk(spec_t rhok, double T, double Te, chargedspec_t muk, spec_t nuk){
  long spec,k;
  double rho;
  

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
    
  
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
  double nuwsum,wsum;
  spec_t w;

  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
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
  double B,p,EoverN,weight,muehi;
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
      case SMAP_NOplus:
        B=4.47e12;
        p=-0.5;
      break;
      case SMAP_O2minus:
        B=3.56e19;
        p=-0.1;
      break;
      default:
        B=0.55/sqrt(_m(k));
        p=-0.5;
    }
    *muk=min(*muk,B*pow(Ek/N,p)/N);
  }
  if (speciestype[k]==SPECIES_ELECTRON){
    EoverN=Ek/N;
    if (EoverN>1e-19){
      /* High-Power Subnanosecond Beams of Runaway Electrons Generated in Dense Gases
         Victor F Tarasenko and Sergei I Yakovlenko, Physica Scripta,
         https://iopscience.iop.org/article/10.1238/Physica.Regular.072a00041
      */
      EoverN=min(EoverN,5e-17);  // muehi is valid between 1e-19 and 5e-17 Vm2
      muehi=(4.0E19*powint(-log(EoverN)-30.0,4)+1.3e40*EoverN)/N; 
      weight=max(0.0,min(1.0,19.0+log(EoverN)/log(10.0)));
      *muk=weight*muehi+(1.0-weight)*(*muk);
    }  
  }
}


double _lnLambda(spec_t rhok, double Te){
  double lnlambda;
#ifdef speceminus  
  double Pe,rhoe;
  rhoe=rhok[speceminus];
  rhoe=max(1e-20,rhoe);
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
static double _mue_from_Nn_Ni_Te_rhok_old(double Nn, double Ni, double Te, spec_t rhok){
  double muen,mue,muei;
  muen=3.74E19*exp(33.5/sqrt(log(max(300.0,Te))))/Nn;
  muei=9.5e16*pow(Te,1.5)/Ni/_lnLambda(rhok,Te);
  //muei=2.5e16*pow(Te,1.5)/Ni;
  mue=1.0/(1.0/muen+1.0/muei);
  return(mue);
}

static double _mue_from_Te_rhok(double Te, spec_t rhok){
  double mue;
  long spec;
  spec_t Nk;
  double N,sum;
  N=0.0;
  for (spec=0; spec<ns; spec++) {
    Nk[spec]=rhok[spec]/_m(spec);
    N+=Nk[spec];
  }
  sum=0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL)
        sum+=(Nk[spec]/N)/_mueNk_from_Te_ParentMacheret(spec, Te);
#ifdef speceminus
    if (speciestype[spec]==SPECIES_IONPLUS || speciestype[spec]==SPECIES_IONMINUS)
        sum+=(Nk[spec]/N)/(1.47*powint(pi,3)*sqr(epsilon0)*_m(speceminus)/(2.0*sqr(_C(spec))*fabs(_C(speceminus))*_lnLambda(rhok,Te))*pow(8.0*kB*Te/pi/_m(speceminus),1.5));
#endif
  }
  assert(sum!=0.0);
  assert(N!=0.0);
  mue=1.0/(sum*N);
  return(mue);
}




/* the ion mobility due to collisions with neutrals is calculated as muin=A*Ti^n/Nn 
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
    mu=_mue_from_Te_rhok(Te,rhok);
  } else {
#endif
    switch (smap[k]){
      case SMAP_eminus:
        mu=_mue_from_Te_rhok(Te,rhok);
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


void find_nuk_from_Dij(spec_t rhok, spec2_t Dij, spec_t nuk){
  long k,i,j;
  double sum,rho,N;
  spec_t chi,Dm;
  
  rho=0;
  for (k=0; k<ns; k++) rho+=rhok[k];
  N=0.0;
  for (k=0; k<ns; k++) N+=rhok[k]/_m(k);
  for (k=0; k<ns; k++) chi[k]=max(1.0e-99,rhok[k]/_m(k)/N); 
  
  switch (MASSDIFFCOEFF){
    case MASSDIFFCOEFF_ORANBORIS:
      for (i=0; i<ns; i++){
        sum=0.0;
        for (j=0; j<ns; j++) {
          if (j!=i) 
            sum+=chi[j]/(Dij[i][j]); 
        }
        Dm[i]=(1.0-chi[i])/(sum+1.0e-20); 
      }
      for (i=0; i<ns; i++){
        nuk[i]=rho*Dm[i]; 
      }
    break;
    default:
      fatal_error("Problem with MASSDIFFCOEFF in find_nuk_from_Dij().");
  }
}






static double _mueN_N2(double Te){
  double mueN_N2;  
  /* 
  Data in log-log coordinates. 
  Obtained with BOLSIG+ using Morgan LXCat cross-sections 
  and experimental data from
  Grigoriev, I. S. and Meilikhov, E. Z., Handbook of Physical Quantities, CRC, Boca Raton, Florida, 1997 
  */

  /* log K */
  double Te_control[] = 
  { 
    5.70378247465620,
    6.00674317023720,
    7.52656892402778,
    9.05804529503584,
    9.28657969532927,
    9.69562262447425,
    10.0010042738634,
    10.7942349131424,
    11.3595487221077,
    13.1074773005257,
    14.4872723773202
  };
  /* log m2/Vs*/
  double mueN_control[] = 
  { 
    58.8996283915835,
    58.6143994480998,
    57.2079523809124,
    56.0504995922214,
    55.8498288967592,
    55.4608930906023,
    55.2620422318571,
    55.2281751623455,
    55.0864976593422,
    54.4855134423581,
    54.2677899585132
  };
  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  mueN_N2 = EXM_f_from_monotonespline(N, Te_control, mueN_control, Te);
  return(exp( mueN_N2 ));
}

static double _mueN_O2(double Te){
  double mueN_O2;  
  /* 
  Data in log-log coordinates. 
  Obtained with BOLSIG+ using Morgan LXCat cross-sections 
  and experimental data from
  Grigoriev, I. S. and Meilikhov, E. Z., Handbook of Physical Quantities, CRC, Boca Raton, Florida, 1997 
  */

  /* log K */
  double Te_control[] = 
  { 
    5.70378247465620,
    7.39303753148020,
    8.15517758348606,
    8.64580049988764,
    9.45446056782523,
    10.2346191249505,
    10.5829258193665,
    11.1835383764881,
    11.8792632932909,
    13.6691438195424,
    14.3029333748767
  };
  /* log m2/Vs*/
  double mueN_control[] = 
  { 
    58.8729601445013,
    57.5646273248511,
    57.2544269409688,
    56.9107008574445,
    56.3934443433482,
    56.1375109692110,
    55.7320458611028,
    55.6126991034703,
    55.1275961353042,
    54.3524788207324,
    54.2715666017936
  };
  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  mueN_O2 = EXM_f_from_monotonespline(N, Te_control, mueN_control, Te);
  return(exp( mueN_O2 ));
}

static double _mueN_NO(double Te){
  double mueN_NO;  
  /* 
  Data in log-log coordinates. 
  Obtained with BOLSIG+ using Phelps LXCat cross-sections 
  and experimental data from

  Bailey, V. and Somerville, J., “XCVIII. The behaviour of electrons in nitric oxide,” The London, Edinburgh, and Dublin
  Philosophical Magazine and Journal of Science Vol. 17, No. 116, 1934, pp. 1169, 1176.
  
  Parkes, D. and Sugden, T., “Electron attachment and detachment in nitric oxide,” Journal of the Chemical Society, Faraday
  Transactions 2: Molecular and Chemical Physics, Vol. 68, 1972, pp. 600, 614.
  */

  /* log K */
  double Te_control[] = 
  { 
    5.70378247465620,
    6.11187508837842,
    6.28632106382936,
    7.06055731604381,
    7.54084158330673,
    7.75370449660376,
    7.99226642016825,
    8.66078964015549,
    9.55059685233927,
    10.7689171050734,
    12.7979646330018,
    14.2900207133957
  };
  /* log m2/Vs*/
  double mueN_control[] = 
  { 
    58.6293380618436,
    58.5552831386806,
    58.5354062420094,
    58.0893558537861,
    57.3686124409252,
    56.9769206567492,
    56.5835313696109,
    56.3153084187583,
    55.9881405125567,
    55.4931539528205,
    54.1576130578674,
    54.0272665183761
  };
  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  mueN_NO = EXM_f_from_monotonespline(N, Te_control, mueN_control, Te);
  return(exp( mueN_NO ));
}

static double _mueN_N(double Te){
  double mueN_N;  
  /* 
  Data in log-log coordinates. 
  Obtained with BOLSIG+ using Morgan LXCat cross-sections 
  */

  /* log K */
  double Te_control[] = 
  { 
    5.70378247465620,
    6.13224239120285,
    6.63877453697547,
    7.53217816278452,
    9.51015326130485,
    9.67250020696868,
    9.90958126612667,
    10.0721001956244,
    10.2909206904909,
    11.4162185004713,
    12.8046835832505,
    14.2629425862465
  };
  /* log m2/Vs*/
  double mueN_control[] = 
  { 
    60.1673170102955,
    59.8889739096267,
    59.4735770614906,
    58.6331239896146,
    56.8895165067773,
    56.7319881654726,
    56.5247555312056,
    56.3966649580482,
    56.2403683546507,
    55.6370481324806,
    55.4009342307237,
    55.3789359833286
  };
  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  mueN_N = EXM_f_from_monotonespline(N, Te_control, mueN_control, Te);
  return(exp( mueN_N ));
}

static double _mueN_O(double Te){
  double mueN_O;  
  /* 
  Data in log-log coordinates. 
  Obtained with BOLSIG+ using Morgan LXCat cross-sections 
  */

  /* log K */
  double Te_control[] = 
  { 
    5.70378247465620,
    6.13952428150941,
    6.70544631640622,
    7.57951709833810,
    8.89727364573419,
    9.48048267293584,
    10.9925001202487,
    11.4348348805704,
    12.3839559088080,
    13.4898268792378,
    14.3823727973700
  };
  /* log m2/Vs*/
  double mueN_control[] = 
  { 
    60.2387759742777,
    59.9607027609325,
    59.4786044268034,
    58.5874379001030,
    57.2323662090014,
    56.7285331349713,
    55.6853472580936,
    55.2478418824560,
    54.1838206328367,
    54.0456594074545,
    54.1328676376063
  };
  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  mueN_O = EXM_f_from_monotonespline(N, Te_control, mueN_control, Te);
  return(exp( mueN_O ));
}

static double _mueN_NH3(double Te){
  double mueN_NH3;  
  /* 
  Data in log-log coordinates. 
  Obtained with BOLSIG+ using Morgan LXCat cross-sections 
  and experimental data from
  Lisovskiy, Valeriy, et al. "Electron drift velocity in NH3 in strong electric fields determined from rf breakdown curves." 
  Journal of Physics D: Applied Physics 38.6 (2005): 872. 
  */

  /* log K */
  double Te_control[] = 
  { 
    5.70378247465620,
    6.19053974785861,
    6.24135109789768,
    7.03617842704541,
    9.56528386809320,
    10.0472135898206,
    10.4282169546571,
    10.7073962257784,
    11.2676001481862,
    11.7445614851530,
    13.0608547382026,
    14.9141228466324
  };
  /* log m2/Vs*/
  double mueN_control[] = 
  { 
    54.0580694275312,
    54.3181707667273,
    54.3946089721228,
    54.6311710173000,
    55.2897567904578,
    55.6265008560765,
    55.7045482187792,
    55.6723602562662,
    55.2485220539277,
    54.9779646059690,
    54.6349916977022,
    54.3457514999829
  };
  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  mueN_NH3 = EXM_f_from_monotonespline(N, Te_control, mueN_control, Te);
  return(exp( mueN_NH3 ));
}


double _mueNk_from_Te_ParentMacheret(long spec, double Te){
  double mueN;
  switch (smap[spec]){
    case SMAP_N2:
      mueN=_mueN_N2(Te);
    break;
    case SMAP_O2:
      mueN=_mueN_O2(Te);
    break;
    case SMAP_NO:
      mueN=_mueN_NO(Te);
    break;
    case SMAP_N:
      mueN=_mueN_N(Te);
    break;
    case SMAP_O:
      mueN=_mueN_O(Te);
    break;
    case SMAP_NH3:
      mueN=_mueN_NH3(Te);
    break;
    case SMAP_NH3v:
      mueN=_mueN_NH3(Te);
    break;
    default:
      mueN=3.74E19*exp(33.5/sqrt(log(Te)));
  }
  return(mueN);
}
