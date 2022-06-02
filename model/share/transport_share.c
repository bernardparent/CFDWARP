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
