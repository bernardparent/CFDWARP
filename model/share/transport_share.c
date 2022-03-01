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


void find_nuk_eta_kappa(spec_t rhok, double T, double Te,
                   spec_t nuk, double *eta, double *kappa){
  
  *eta=_eta_from_rhok_T_Te(rhok,T,Te);
  *kappa=_kappa_from_rhok_T_Te(rhok, T, Te);
  find_nuk_from_rhok_T_Te(rhok, T, Te, nuk);  
}


// find the scalar thermal conductivity of the kth charged species
double _kappac_from_rhok_Tk_Ek(spec_t rhok, double T, double E, long k){
  double mu,kappa,cp;
  if (k>=ncs) fatal_error("_kappac_from_rhok_Tk_Ek can only be used for charged species, not for species %ld.",k);
  mu=_muk_from_rhok_Tk_Ek(rhok, T, E, k);
  cp=_cpk_from_T_equilibrium(k,T);
  //cp=5.0/2.0*kB/_m(k);
  kappa=cp*rhok[k]*kB*T*mu/fabs(_C(k));
  return(kappa);
}


// find the scalar viscosity of the kth charged species
double _etac_from_rhok_Tk_Ek(spec_t rhok, double T, double E, long k){
  double kappa,Pr,cp,eta;
  if (k>=ncs) fatal_error("_etac_from_rhok_Tk_Ek can only be used for charged species, not for species %ld.",k);
  kappa=_kappac_from_rhok_Tk_Ek(rhok, T, E, k);
  // Pr=cp*eta/kappa

  if (speciestype[k]==SPECIES_ELECTRON) Pr=0.73; else Pr=0.96;
  cp=_cpk_from_T_equilibrium(k,T);
  //if (fabs(cp-5.0/2.0*kB/_m(k))/min(cp,5.0/2.0*kB/_m(k))>2.2 && k==specN2plus ) printf("x");
  eta=Pr*kappa/cp;
  return(eta);
}


