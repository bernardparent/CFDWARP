// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2026 Bernard Parent


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


#include <model/share/fluid_vib.h>
#include <model/share/fluid_share.h>
#include <model/share/emfield_share.h>
#include <model/fluid/_fluid.h>
#include <model/thermo/_thermo.h>
#include <model/thermo/_generic/electronimpact.h>
#include <model/transport/_transport.h>
#include <model/emfield/_emfield.h>
#include <model/chem/_chem.h>
#include <model/metrics/_metrics.h>
#include <cycle/_cycle.h>
#include <cycle/share/cycle_share.h>
#include <cycle/share/res_share.h>
#include <model/fluid/.active/fluid.h>
#include <model/share/model_share.h>
#include <src/bdry.h>


/* ======================================================================================
   Vibrational energy equation, common to the fluid modules with _FLUID_VIBMODEL:
   - VibSpecies(), which sets the species with vibrationally excited states (the set V);
   - the reading of TAUVTMODEL, QEVMODEL and VIBDISTTYPE;
   - the vibrational energy ev (per unit mass of the vibrating species) and the
     vibrational temperature, obtained through find_vibdist() of the thermo module;
   - the V-T relaxation (TAUVTMODEL) and the e-V excitation (QEVMODEL) source terms;
   - the conduction of the vibrational energy (kappav, nuv and their turbulent counterparts).
   ====================================================================================== */

#if _FLUID_VIBMODEL

#define TAUVT_EQUILIBRIUM_CONSTANT 1e-12


/* mass fraction of the vibrating species */
double _wv(np_t np){
  long spec;
  double wv;
  wv=0.0;
  for (spec=0; spec<ns; spec++) if (_isvib(spec)) wv+=_w(np,spec);
  return(wv);
}


/* vibrational energy per unit mass of the vibrating species, J/kg; 0 if no species vibrates */
double _ev (np_t np){
  double ret,wv;
  assert_np(np,_rho(np)!=0.0e0);
  wv=_wv(np);
  ret=0.0;
  if (wv!=0.0) ret=np.bs->U[fluxev]/(_rho(np)*wv);
  return(ret);
}


/* vibrational energy per unit mass of the vibrating species at equilibrium at the gas temperature */
double _evzero (np_t np, gl_t *gl){
  spec_t w;
  find_w(np,w);
  return(_ev_from_w_Tv(gl->model.fluid.VIBDISTTYPE,w,_T(np,gl)));
}


void find_vibdist_from_np(np_t np, gl_t *gl, vibdist_t *vibdist){
  spec_t Nk;
  long spec;
  for (spec=0; spec<ns; spec++) Nk[spec]=_rhok(np,spec)/_m(spec);
  find_vibdist(gl, gl->model.fluid.VIBDISTTYPE, Nk, _T(np,gl), _ev(np), 0.0, vibdist);
}


double _Tv(np_t np, gl_t *gl){
  vibdist_t vibdist;
  find_vibdist_from_np(np, gl, &vibdist);
  return(vibdist.Tv);
}


/* vibrational temperature of species spec, from the population of its first vibrational level
   (T if spec does not vibrate) */
double _Tvk(np_t np, gl_t *gl, long spec){
  vibdist_t vibdist;
  if (spec<0 || spec>=ns) fatal_error("Species %ld is not part of the chemical solver in _Tvk().",spec);
  find_vibdist_from_np(np, gl, &vibdist);
  return(vibdist.T1k[spec]);
}


/* Tv given the mass fractions w, the gas temperature T and ev */
double _Tv_from_w_T_ev(gl_t *gl, spec_t w, double T, double ev){
  spec_t Nk;
  long spec;
  vibdist_t vibdist;
  for (spec=0; spec<ns; spec++) Nk[spec]=w[spec]/_m(spec);
  find_vibdist(gl, gl->model.fluid.VIBDISTTYPE, Nk, T, ev, 0.0, &vibdist);
  return(vibdist.Tv);
}


/* clip ev between its values at Tvmin and Tvmax */
void reformat_ev(gl_t *gl, spec_t w, double *ev, char *suffix, bool *flag){
  double evmin,evmax;
  evmin=_ev_from_w_Tv(gl->model.fluid.VIBDISTTYPE,w,gl->model.fluid.Tvmin);
  evmax=_ev_from_w_Tv(gl->model.fluid.VIBDISTTYPE,w,gl->model.fluid.Tvmax);
  if (*ev<evmin) {
    add_to_clipped_variables2(gl,"Tvmin",suffix);
    *ev=evmin;
    *flag=TRUE;
  }
  if (*ev>evmax) {
    add_to_clipped_variables2(gl,"Tvmax",suffix);
    *ev=evmax;
    *flag=TRUE;
  }
}


/* the derivatives of ev with respect to U */
void find_dev_dU(np_t np, gl_t *gl, flux_t devdU){
  long flux,spec;
  double rhov;
  for (flux=0; flux<nf; flux++) devdU[flux]=0.0;
  rhov=_rho(np)*_wv(np);
  if (rhov!=0.0){
    devdU[fluxev]=1.0/rhov;
    for (spec=0; spec<ns; spec++) 
      if (_isvib(spec)) devdU[spec]=-_ev(np)/rhov;
  }
}


/* Tv given ev and the composition and the gas temperature of node np */
double _Tv_from_np_ev(np_t np, gl_t *gl, double ev){
  spec_t Nk;
  long spec;
  vibdist_t vibdist;
  for (spec=0; spec<ns; spec++) Nk[spec]=_rhok(np,spec)/_m(spec);
  find_vibdist(gl, gl->model.fluid.VIBDISTTYPE, Nk, _T(np,gl), ev, 0.0, &vibdist);
  return(vibdist.Tv);
}


void find_dTv_dU(np_t np, gl_t *gl, flux_t dTvdU){
  long flux,spec;
  double rhov;
  vibdist_t vibdist;
  for (flux=0; flux<nf; flux++) dTvdU[flux]=0.0;
  find_vibdist_from_np(np, gl, &vibdist);
  if (!vibdist.EMPTY){
    rhov=_rho(np)*_wv(np);
    dTvdU[fluxev]=1.0/vibdist.cv/rhov;
    for (spec=0; spec<ns; spec++) 
      if (_isvib(spec)) dTvdU[spec]=-vibdist.evk[spec]/vibdist.cv/rhov;
  }
}


void read_model_fluid_actions_vib_species(char *actionname, char **argum, SOAP_codex_t *codex){
  long spec,specfound,cnt,numargum;
  bool isvibspec[ns];
  char *specname,*argumname,*list,*statename;
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (strcmp(actionname,"VibSpecies")==0) {
    if (gl->model.fluid.VIBSPECIESREAD) SOAP_fatal_error(codex,"VibSpecies() can not be called more than once.");
    specname=(char *)malloc(sizeof(char));
    argumname=(char *)malloc(sizeof(char));
    list=(char *)malloc(sizeof(char));
    statename=(char *)malloc(sizeof(char));
    for (spec=0; spec<ns; spec++) isvibspec[spec]=FALSE;
    numargum=SOAP_number_argums(*argum);
    for (cnt=0; cnt<numargum; cnt++){
      SOAP_get_argum_string(codex, &argumname, *argum, cnt);
      specfound=-1;
      for (spec=0; spec<ns; spec++){
        find_species_name(spec, &specname);
        if (strcmp(specname,argumname)==0) specfound=spec;
      }
      if (specfound<0) 
        SOAP_fatal_error(codex,"Species %s listed in VibSpecies() is not a species of the chemical solver.",argumname);
      if (isvibspec[specfound]) 
        SOAP_fatal_error(codex,"Species %s is listed more than once in VibSpecies().",argumname);
      if (!is_vib_data_available(specfound)) {
        find_vib_data_species_list(&list);
        SOAP_fatal_error(codex,"Species %s listed in VibSpecies() has no vibrational data in model/thermo/_generic/vibration.c. Species with vibrational data: %s.",argumname,list);
      }
      if (is_vib_state_in_chem(specfound,&statename)) 
        SOAP_fatal_error(codex,"Species %s can not be listed in VibSpecies() because the chemical solver includes its vibrationally excited state %s as a separate species.",argumname,statename);
      isvibspec[specfound]=TRUE;
    }
    set_vibspecies(isvibspec);
    gl->model.fluid.VIBSPECIESREAD=TRUE;
    free(specname);
    free(argumname);
    free(list);
    free(statename);
    codex->ACTIONPROCESSED=TRUE;
  }
}


/* to be called before processing the Fluid() module of the control file */
void add_vib_constants_to_codex(SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  SOAP_add_int_to_vars(codex,"TAUVTMODEL_PARENT2024",TAUVTMODEL_PARENT2024);
  SOAP_add_int_to_vars(codex,"TAUVTMODEL_MILLIKANWHITE",TAUVTMODEL_MILLIKANWHITE);
  SOAP_add_int_to_vars(codex,"TAUVTMODEL_TEQUILIBRIUM",TAUVTMODEL_TEQUILIBRIUM);
  SOAP_add_int_to_vars(codex,"QEVMODEL_NONE",QEVMODEL_NONE);
  SOAP_add_int_to_vars(codex,"QEVMODEL_ALEKSANDROV",QEVMODEL_ALEKSANDROV);
  SOAP_add_int_to_vars(codex,"QEVMODEL_PARENT2024",QEVMODEL_PARENT2024);
  SOAP_add_int_to_vars(codex,"QEVMODEL_RODRIGUEZ2025",QEVMODEL_RODRIGUEZ2025);
  SOAP_add_int_to_vars(codex,"QEVMODEL_RODRIGUEZ2026",QEVMODEL_RODRIGUEZ2026);
  SOAP_add_int_to_vars(codex,"VIBDISTTYPE_BOLTZMANN",VIBDISTTYPE_BOLTZMANN);
  /* deprecated constants: they remain defined so that the deprecation message is given */
  SOAP_add_int_to_vars(codex,"N2VIBMODEL_ALEKSANDROV",N2VIBMODEL_ALEKSANDROV);
  SOAP_add_int_to_vars(codex,"N2VIBMODEL_TEQUILIBRIUM",N2VIBMODEL_TEQUILIBRIUM);
  SOAP_add_int_to_vars(codex,"N2VIBMODEL_PARENT2024",N2VIBMODEL_PARENT2024);
  SOAP_add_int_to_vars(codex,"N2VIBMODEL_RODRIGUEZ2025",N2VIBMODEL_RODRIGUEZ2025);
  SOAP_add_int_to_vars(codex,"N2VIBMODEL_RODRIGUEZ2026",N2VIBMODEL_RODRIGUEZ2026);
  SOAP_add_int_to_vars(codex,"N2VIBMODEL_MILLIKAN",N2VIBMODEL_MILLIKAN);
  gl->model.fluid.VIBSPECIESREAD=FALSE;
}


static void fatal_error_N2VIBMODEL_deprecated(SOAP_codex_t *codex){
  int N2VIBMODEL;
  char tauvt[100],qev[100];
#if defined(_FLUID_FAVREREYNOLDS) && defined(_FLUID_PLASMA)
  bool N2VIBJOULEHEATING;
#endif
  find_int_var_from_codex(codex,"N2VIBMODEL",&N2VIBMODEL);
  switch (N2VIBMODEL){
    case N2VIBMODEL_PARENT2024:
      strcpy(tauvt,"TAUVTMODEL_PARENT2024");
      strcpy(qev,"QEVMODEL_PARENT2024");
    break;
    case N2VIBMODEL_ALEKSANDROV:
      strcpy(tauvt,"TAUVTMODEL_PARENT2024");
      strcpy(qev,"QEVMODEL_ALEKSANDROV");
    break;
    case N2VIBMODEL_RODRIGUEZ2025:
      strcpy(tauvt,"TAUVTMODEL_PARENT2024");
      strcpy(qev,"QEVMODEL_RODRIGUEZ2025");
    break;
    case N2VIBMODEL_RODRIGUEZ2026:
      strcpy(tauvt,"TAUVTMODEL_MILLIKANWHITE");
      strcpy(qev,"QEVMODEL_RODRIGUEZ2026");
    break;
    case N2VIBMODEL_MILLIKAN:
      strcpy(tauvt,"TAUVTMODEL_MILLIKANWHITE");
      strcpy(qev,"QEVMODEL_NONE");
    break;
    case N2VIBMODEL_TEQUILIBRIUM:
      strcpy(tauvt,"TAUVTMODEL_TEQUILIBRIUM");
      strcpy(qev,"QEVMODEL_NONE");
    break;
    default:
      strcpy(tauvt,"TAUVTMODEL_PARENT2024");
      strcpy(qev,"QEVMODEL_NONE");
  }
#if !_FLUID_EENERGY
  strcpy(qev,"QEVMODEL_NONE");
#endif
#if defined(_FLUID_FAVREREYNOLDS) && defined(_FLUID_PLASMA)
  /* the e-V excitation was set by N2VIBJOULEHEATING */
  strcpy(qev,"QEVMODEL_NONE");
  if (N2VIBMODEL!=N2VIBMODEL_TEQUILIBRIUM && SOAP_is_var_in_codex(codex,"N2VIBJOULEHEATING")) {
    find_bool_var_from_codex(codex,"N2VIBJOULEHEATING",&N2VIBJOULEHEATING);
    if (N2VIBJOULEHEATING) strcpy(qev,"QEVMODEL_ALEKSANDROV");
  }
#endif
  fatal_error("N2VIBMODEL is deprecated within Model(%s()). Replace it by the following:\n\n"
              "    VibSpecies(\"N2\");\n"
              "    TAUVTMODEL=%s;\n"
              "    QEVMODEL=%s;\n"
              "    VIBDISTTYPE=VIBDISTTYPE_BOLTZMANN;\n",_FLUID_ACTIONNAME,tauvt,qev);
}


/* to be called after processing the Fluid() module of the control file */
void read_model_fluid_vib(SOAP_codex_t *codex){
  bool N2ONLY;
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (SOAP_is_var_in_codex(codex, "N2VIBMODEL")) fatal_error_N2VIBMODEL_deprecated(codex);
  if (SOAP_is_var_in_codex(codex, "N2VIBJOULEHEATING")) fatal_error("The parameter N2VIBJOULEHEATING is deprecated. Replace N2VIBJOULEHEATING=YES; by QEVMODEL=QEVMODEL_ALEKSANDROV; and N2VIBJOULEHEATING=NO; by QEVMODEL=QEVMODEL_NONE;.");
  if (SOAP_is_var_in_codex(codex, "N2VIBJOULEHEATINGMODEL")) fatal_error("The parameter N2VIBJOULEHEATINGMODEL is deprecated. Use QEVMODEL instead.");
  if (!gl->model.fluid.VIBSPECIESREAD) 
    fatal_error("VibSpecies() must be called within Model(%s()). List the species with vibrationally excited states, for example VibSpecies(\"N2\");. Use VibSpecies(); if no species has vibrationally excited states.",_FLUID_ACTIONNAME);

  find_int_var_from_codex(codex,"TAUVTMODEL",&gl->model.fluid.TAUVTMODEL);
  if (gl->model.fluid.TAUVTMODEL!=TAUVTMODEL_PARENT2024 && gl->model.fluid.TAUVTMODEL!=TAUVTMODEL_MILLIKANWHITE 
   && gl->model.fluid.TAUVTMODEL!=TAUVTMODEL_TEQUILIBRIUM)
    SOAP_fatal_error(codex,"TAUVTMODEL must be set to TAUVTMODEL_PARENT2024, TAUVTMODEL_MILLIKANWHITE or TAUVTMODEL_TEQUILIBRIUM.");
  find_int_var_from_codex(codex,"QEVMODEL",&gl->model.fluid.QEVMODEL);
  if (gl->model.fluid.QEVMODEL!=QEVMODEL_NONE && gl->model.fluid.QEVMODEL!=QEVMODEL_ALEKSANDROV 
   && gl->model.fluid.QEVMODEL!=QEVMODEL_PARENT2024 && gl->model.fluid.QEVMODEL!=QEVMODEL_RODRIGUEZ2025 
   && gl->model.fluid.QEVMODEL!=QEVMODEL_RODRIGUEZ2026)
    SOAP_fatal_error(codex,"QEVMODEL must be set to QEVMODEL_NONE, QEVMODEL_ALEKSANDROV, QEVMODEL_PARENT2024, QEVMODEL_RODRIGUEZ2025 or QEVMODEL_RODRIGUEZ2026.");
  find_int_var_from_codex(codex,"VIBDISTTYPE",&gl->model.fluid.VIBDISTTYPE);
  if (gl->model.fluid.VIBDISTTYPE!=VIBDISTTYPE_BOLTZMANN)
    SOAP_fatal_error(codex,"VIBDISTTYPE must be set to VIBDISTTYPE_BOLTZMANN. Other distributions of the vibrational states are not available yet.");

  N2ONLY=FALSE;
#ifdef specN2
  N2ONLY=(_numvibspecies()==1 && _isvib(specN2));
#endif
  if (gl->model.fluid.TAUVTMODEL==TAUVTMODEL_PARENT2024 && !N2ONLY)
    SOAP_fatal_error(codex,"TAUVTMODEL_PARENT2024 is a fit for N2 and can only be used with VibSpecies(\"N2\");. With other vibrating species or with no vibrating species, use TAUVTMODEL_MILLIKANWHITE or TAUVTMODEL_TEQUILIBRIUM.");
  if (gl->model.fluid.QEVMODEL!=QEVMODEL_NONE && !N2ONLY)
    SOAP_fatal_error(codex,"The e-V excitation models are fits for N2 and require VibSpecies(\"N2\");. With other vibrating species or with no vibrating species, set QEVMODEL=QEVMODEL_NONE;.");
  if (gl->model.fluid.TAUVTMODEL==TAUVTMODEL_TEQUILIBRIUM && gl->model.fluid.QEVMODEL!=QEVMODEL_NONE)
    SOAP_fatal_error(codex,"TAUVTMODEL_TEQUILIBRIUM requires QEVMODEL=QEVMODEL_NONE;.");
}


/*
Macheret, Sergey O., et al. "Electron-beam-generated plasmas in hypersonic magnetohydrodynamic channels." 
AIAA journal 39.6 (2001): 1127-1138. 
Macheret, S. O., Martinelli, L., and Miles, R. B., "Shock Wave Propagation and Structure in Non-Uniform Gases and
Plasmas," 1999, AIAA Paper 99-0598. 
- VT relaxation time of N2 in air */
static double _tauvt_Macheret(spec_t rhok, double T){
  double N,N_O,tauvt;
  long spec;  
  tauvt=1.0e99;
  N=0.0;
  N_O=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
#ifdef specO
  N_O=rhok[specO]/_m(specO);
#else
  N_O=0.0;
#endif
  tauvt=1.0/(N*7E-16*exp(-141.0/pow(T,1.0/3.0))
       +N_O*5E-18*exp(-128.0/sqrt(T)));
  return(tauvt);
}


/*
Gnoffo, Peter A. Conservation equations and physical models for hypersonic air flows in thermal and chemical nonequilibrium. Vol. 2867. 
National Aeronautics and Space Administration, Office of Management, Scientific and Technical Information Division, 1989.
- VT relaxation time (s) of species k in collisions with the other species excluding electrons and negative ions
- Eqs (55) through (58), (61)
- the pressure is found from rhok and T so that the derivatives of 1/tau with respect to rhok and T are complete
Carter, N., Aiken, T. T., & Boyd, I. D. (2025). Sensitivity Analysis and Uncertainty Quantification of Plasma Formation in 
Nonequilibrium Hypersonic Airflows. In AIAA SCITECH 2025 Forum (p. 0950).
- Eqs (18) through (21) 
*/
static double _tauvtk_MillikanWhite(long k, spec_t rhok, double T){
  double sum, Nsum, tauvt, mu_sr, A_sr, B_sr, thetav, Nk, P;
  long spec;  
  sum=0.0;
  Nsum=0.0;
  thetav=_thetavk(k);
  P=0.0;
  for (spec=0; spec<ns; spec++) P+=rhok[spec]/_m(spec)*kB*T;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]!=SPECIES_ELECTRON && speciestype[spec]!=SPECIES_IONMINUS){
      mu_sr = 1.0e3 * _calM(k)*_calM(spec)/(_calM(k)+_calM(spec)); // reduced mass, g/mol
      A_sr = 1.16e-3 * pow(thetav,4.0/3.0) * pow(mu_sr,0.5); // Millikan-White A factor
      B_sr = 0.015*pow(mu_sr,1.0/4.0); // Millikan-White B factor
      sum += rhok[spec]/_m(spec) / ( (101325.0e0/P) * exp( A_sr * ( 1.0/pow(T,1.0/3.0) - B_sr ) - 18.42 ) ); // Gnoffo ref., eq.(55)
      Nsum += rhok[spec]/_m(spec);
    }
  }
  Nk=max(rhok[k]/_m(k),1.0e-99);
  tauvt=Nsum/sum + 1.0/(Nk*3.00E-21*pow(5.0e4/T,2.0)*pow(8*kB*T/pi/_m(k),0.5)); // apply Park high T correction
  return(tauvt);
}


/* V-T relaxation time of each vibrating species */
static void find_tauvtk(gl_t *gl, spec_t rhok, double T, spec_t tauvtk){
  long spec,k;
  double rho;
  rho=0.0;
  for (spec=0; spec<ns; spec++) {
    rho+=rhok[spec];
    tauvtk[spec]=1.0e99;
  }
  for (k=0; k<ns; k++){
    if (_isvib(k)){
      switch (gl->model.fluid.TAUVTMODEL){
        case TAUVTMODEL_TEQUILIBRIUM:
          tauvtk[k]=TAUVT_EQUILIBRIUM_CONSTANT/rho;
        break;
        case TAUVTMODEL_PARENT2024:
          tauvtk[k]=_tauvt_Macheret(rhok,T);
        break;
        case TAUVTMODEL_MILLIKANWHITE:
          tauvtk[k]=_tauvtk_MillikanWhite(k,rhok,T);
        break;
        default:
          fatal_error("Problem with TAUVTMODEL in find_tauvtk()."); 
      }
    }
  }
}


/* V-T relaxation time of each vibrating species and the derivatives of its inverse */
static void find_tauvtk_d1overtauvtk(gl_t *gl, spec_t rhok, double T, spec_t tauvtk, spec_t d1overtauvtkdT, spec2_t d1overtauvtkdrhok){
  spec_t rhok2,drhok;
  double rho, dT;
  long spec2,spec,k;
  rho=0.0;
  dT=T*1.0e-4;
  for (spec=0; spec<ns; spec++) {
    rho+=rhok[spec];
    drhok[spec]=0.0001*rhok[spec]+1.0e-30;
    tauvtk[spec]=1.0e99;
    d1overtauvtkdT[spec]=0.0;
    for (spec2=0; spec2<ns; spec2++) d1overtauvtkdrhok[spec][spec2]=0.0;
  }
  for (k=0; k<ns; k++){
    if (_isvib(k)){
      switch (gl->model.fluid.TAUVTMODEL){
        case TAUVTMODEL_TEQUILIBRIUM:
          tauvtk[k] = TAUVT_EQUILIBRIUM_CONSTANT/rho;
          for (spec=0; spec<ns; spec++) 
            d1overtauvtkdrhok[k][spec]=1.0/TAUVT_EQUILIBRIUM_CONSTANT;
        break;
        case TAUVTMODEL_PARENT2024: 
          tauvtk[k]=_tauvt_Macheret(rhok,T);
          for (spec=0; spec<ns; spec++) {
            d1overtauvtkdT[k]+=7E-16*exp(-141.0/pow(T,1.0/3.0))*rhok[spec]*calA/_calM(spec)*(141.0/3.0*pow(T,-4.0/3.0));
            d1overtauvtkdrhok[k][spec]=7E-16*exp(-141.0/pow(T,1.0/3.0))*calA/_calM(spec);
          }
#ifdef specO
          d1overtauvtkdrhok[k][specO]+=5E-18*exp(-128.0/pow(T,1.0/2.0))*calA/_calM(specO);
          d1overtauvtkdT[k]+=5E-18*exp(-128.0/pow(T,1.0/2.0))*rhok[specO]*calA/_calM(specO)*(128.0/2.0*pow(T,-3.0/2.0));
#endif
        break;
        case TAUVTMODEL_MILLIKANWHITE:
          tauvtk[k]=_tauvtk_MillikanWhite(k,rhok,T);
          d1overtauvtkdT[k]=(1.0/_tauvtk_MillikanWhite(k,rhok,T+dT)-1.0/tauvtk[k])/dT; 
          for (spec=0; spec<ns; spec++){
            for (spec2=0; spec2<ns; spec2++) rhok2[spec2]=rhok[spec2];
            rhok2[spec]+=drhok[spec]; // perturb rhok[spec] only
            d1overtauvtkdrhok[k][spec]=(1.0/_tauvtk_MillikanWhite(k,rhok2,T)-1.0/tauvtk[k])/drhok[spec];        
          }
        break;
        default:
          fatal_error("Problem with TAUVTMODEL in find_tauvtk_d1overtauvtk()."); 
      }
    }
  }
}


/* the V-T relaxation source term */
void find_Svib_VT(np_t *np, gl_t *gl, long l, flux_t S){
  spec_t rhok,tauvtk;
  double T;
  long flux,spec;
  vibdist_t vibdist;

  T=_T(np[l],gl);
  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);
  for (flux=0; flux<nf; flux++) S[flux]=0.0e0;
  find_tauvtk(gl,rhok,T,tauvtk);
  find_vibdist_from_np(np[l],gl,&vibdist);
  for (spec=0; spec<ns; spec++){
    if (_isvib(spec)){
      assert(tauvtk[spec]!=0.0);
      S[fluxev]+=rhok[spec]/tauvtk[spec]*(_evk_from_T(spec,T)-vibdist.evk[spec]);
    }
  }
}


/* the Jacobian of the V-T relaxation source term S=sum_k rhok*(evzerok-evk)/tauk; dTdU is the
   derivative of T with respect to U given by the fluid module. Only the loss part of the source term,
   -sum_k rhok*evk/tauk, is linearized: the derivatives of the gain part, sum_k rhok*evzerok/tauk, would
   make the implicit operator lose diagonal dominance, and are included only when TEST is defined so
   that the Jacobian can be verified against the numerical one. */
void find_dSvib_VT_dU(np_t *np, gl_t *gl, long l, flux_t dTdU, sqmat_t dSvibdU){
  long row,col,spec,k,flux;
  spec_t rhok,tauvtk,d1overtauvtkdT;
#ifdef TEST
  spec_t evzerok;
#endif
  spec2_t d1overtauvtkdrhok;
  double T,sumrhocv,sumfovertau;
  vibdist_t vibdist;

  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);
  for (row=0; row<nf; row++) {
    for (col=0; col<nf; col++) {
      dSvibdU[row][col]=0.0e0;
    }
  }
  find_vibdist_from_np(np[l],gl,&vibdist);
  /* when the vibrating species are absent, the source term is zero */
  if (vibdist.EMPTY) return;
  T=_T(np[l],gl);
  find_tauvtk_d1overtauvtk(gl,rhok,T,tauvtk,d1overtauvtkdT,d1overtauvtkdrhok);
  /* loss part: wrt rho*ev through Tv */
  sumrhocv=0.0;
  for (k=0; k<ns; k++) if (_isvib(k)) sumrhocv+=rhok[k]*vibdist.cvk[k];
  sumfovertau=0.0;
  for (k=0; k<ns; k++) if (_isvib(k)) sumfovertau+=rhok[k]*vibdist.cvk[k]/sumrhocv/tauvtk[k];
  dSvibdU[fluxev][fluxev]=-sumfovertau;

  /* loss part: wrt rhok at constant rho*ev, through Tv; this is zero with one vibrating species */
  for (k=0; k<ns; k++) if (_isvib(k)) 
    dSvibdU[fluxev][k]+=vibdist.evk[k]*(sumfovertau-1.0/tauvtk[k]);

  /* loss part: wrt 1/tauvtk */
  for (spec=0; spec<ns; spec++) 
    for (k=0; k<ns; k++) if (_isvib(k)) 
      dSvibdU[fluxev][spec]+=-rhok[k]*vibdist.evk[k]*d1overtauvtkdrhok[k][spec];
  for (flux=0; flux<nf; flux++)
    for (k=0; k<ns; k++) if (_isvib(k)) 
      dSvibdU[fluxev][flux]+=-rhok[k]*vibdist.evk[k]*d1overtauvtkdT[k]*dTdU[flux];

#ifdef TEST
  for (k=0; k<ns; k++) if (_isvib(k)) evzerok[k]=_evk_from_T(k,T);
  /* gain part: wrt rhok*evzerok */
  for (k=0; k<ns; k++) if (_isvib(k)) dSvibdU[fluxev][k]+=evzerok[k]/tauvtk[k];
  for (flux=0; flux<nf; flux++) 
    for (k=0; k<ns; k++) if (_isvib(k)) 
      dSvibdU[fluxev][flux]+=rhok[k]*_cvk_from_T(k,T)*dTdU[flux]/tauvtk[k];
  /* gain part: wrt 1/tauvtk */
  for (spec=0; spec<ns; spec++) 
    for (k=0; k<ns; k++) if (_isvib(k)) 
      dSvibdU[fluxev][spec]+=rhok[k]*evzerok[k]*d1overtauvtkdrhok[k][spec];
  for (flux=0; flux<nf; flux++)
    for (k=0; k<ns; k++) if (_isvib(k)) 
      dSvibdU[fluxev][flux]+=rhok[k]*evzerok[k]*d1overtauvtkdT[k]*dTdU[flux];
#endif
}


/* display the vibrational energy source term Jacobian found analytically and numerically */
void test_dSvib_dU_display(gl_t *gl, sqmat_t dSvibdU, sqmat_t dSvibdUnum){
  long row,col;
  bool ROWFILLED;
  double diff,diffmax;
  printf("\n");
  printf("Rows of the Jacobian of the V-T relaxation source term\n");
  printf("The analytical Jacobian should match the numerical one when the -DTEST compile flag is set.\n");
  printf("The column 'scaled diff' is (analytical-numerical)*Uref[col].\n\n");
  for (row=0; row<nf; row++){
    ROWFILLED=FALSE;
    for (col=0; col<nf; col++) if (dSvibdU[row][col]!=0.0 || dSvibdUnum[row][col]!=0.0) ROWFILLED=TRUE;
    if (ROWFILLED) {
      printf("row %ld:\n",row);
      printf("  col      analytical       numerical     scaled diff\n");
      diffmax=0.0;
      for (col=0; col<nf; col++){
        diff=(dSvibdU[row][col]-dSvibdUnum[row][col])*gl->cycle.fluid.Uref[col];
        diffmax=max(diffmax,fabs(diff));
        printf("  %3ld  %+14.6E  %+14.6E  %+14.6E\n",col,dSvibdU[row][col],dSvibdUnum[row][col],diff);
      }
      printf("  maximum scaled difference: %E\n\n",diffmax);
    }
  }
}


#ifdef _FLUID_PLASMA

/* e-V excitation: the part of the electron energy lost to the vibrational energy, W/m3;
   emfieldfact multiplies the Joule heating */
double _Qev(np_t *np, gl_t *gl, long l, double emfieldfact){
  double Qev;
#ifdef specN2
  long spec;
  spec_t N;
  double chiN2,Te,Tv,Tref,Ntotal,thetavib,f0,f1,f2,f3;
#endif

  Qev=0.0;
  switch (gl->model.fluid.QEVMODEL){
    case QEVMODEL_NONE:
      Qev=0.0;
    break;
#ifdef specN2
    case QEVMODEL_PARENT2024:
      chiN2=_Nk(np[l],gl,specN2)/_N(np[l],gl);
      Qev=_zetav_from_Te_Parent2024(_Te(np[l]))*chiN2*_Qjoulek(np, gl, l, speceminus)*emfieldfact;
    break;
    case QEVMODEL_ALEKSANDROV:
      chiN2=_Nk(np[l],gl,specN2)/_N(np[l],gl);
      Qev=_zetav_from_Te_Aleksandrov(_Te(np[l]))*chiN2*_Qjoulek(np, gl, l, speceminus)*emfieldfact;
    break;
    case QEVMODEL_RODRIGUEZ2025:
      Ntotal=0.0;
      for (spec=0; spec<ns; spec++){
        N[spec]=_rhok(np[l],spec)/_m(spec); 
        Ntotal+=N[spec];
      }
      Tv=_Tv(np[l],gl);
      Te=_Te(np[l]);
      Tref=300.0e0; //gas temperature at which the experiments or BOLSIG+ results of mueN and E/N were obtained
      thetavib=_thetavk(specN2); //characteristic N2 vibration temperature in K
      // electron cooling-heating due to e-V interactions
      Qev=_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverNk_from_Te_chieminus(specN2, Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(1.0-exp( thetavib/Te - thetavib/Tv ));
    break;
    case QEVMODEL_RODRIGUEZ2026:
      Ntotal=0.0;
      for (spec=0; spec<ns; spec++){
        N[spec]=_rhok(np[l],spec)/_m(spec); 
        Ntotal+=N[spec];
      }
      Tv=_Tv(np[l],gl);
      Te=_Te(np[l]);
      Tref=300.0e0; //gas temperature at which the experiments or BOLSIG+ results of mueN and E/N were obtained
      thetavib=_thetavk(specN2); //characteristic N2 vibration temperature in K
      f0 = 1.0-exp(-thetavib/Tv);
      f1 = exp(-1.0*thetavib/Tv)*(1.0-exp(-thetavib/Tv));
      f2 = exp(-2.0*thetavib/Tv)*(1.0-exp(-thetavib/Tv));
      f3 = 1.0-f0-f1-f2;
      // electron cooling-heating due to e-V interactions
      Qev=_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f0*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverNk_from_Te_chieminus(specN2, Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(1.0-exp( thetavib/Te - thetavib/Tv ))
               +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f1*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v1(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(1.0-exp( thetavib/Te - thetavib/Tv ))
               +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f2*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v2(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(1.0-exp( thetavib/Te - thetavib/Tv ))
               +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f3*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v3(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(1.0-exp( thetavib/Te - thetavib/Tv ));
    break;
#endif
    default:
      fatal_error("Problem with QEVMODEL in _Qev().");
  } 
  return(Qev);
}


/* derivatives of _Qev() with respect to U: only the de-excitation parts of the 
   RODRIGUEZ2025 and RODRIGUEZ2026 models, through Tv and Te, are included */
void find_dQev_dU(np_t *np, gl_t *gl, long l, sqmat_t dQevdU){
  long row,col;
#ifdef specN2
  spec_t N;
  long flux,spec;
  double Ntotal,Tv,Te,Tref,thetavib,f0,f1,f2,f3;
  flux_t dTvdU,dTedU;
#endif

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dQevdU[row][col]=0.0;
    }
  }	
  switch (gl->model.fluid.QEVMODEL){
#ifdef specN2
    case QEVMODEL_RODRIGUEZ2025:
      Ntotal=0.0;
      for (spec=0; spec<ns; spec++){
        N[spec]=_rhok(np[l],spec)/_m(spec); 
        Ntotal+=N[spec];
      }
      Tv=_Tv(np[l],gl);
      Te=_Te(np[l]);
      Tref=300.0e0; //gas temperature at which the experiments or BOLSIG+ results of mueN and E/N were obtained
      thetavib=_thetavk(specN2); //characteristic N2 vibration temperature in K
      // electron cooling-heating due to e-V interactions
      find_dTv_dU(np[l], gl, dTvdU);  
      find_dTe_dU(np[l], gl, dTedU);  
      for (flux=0; flux<nf; flux++){
        // wrt Tv
        dQevdU[fluxev][flux]+=_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverNk_from_Te_chieminus(specN2, Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*exp(- thetavib/Tv)*(thetavib/sqr(Tv)) )*dTvdU[flux];
        // wrt Te
        dQevdU[fluxev][flux]+=_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverNk_from_Te_chieminus(specN2, Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*(-thetavib/sqr(Te))*exp(- thetavib/Tv) )*dTedU[flux];
      }
    break;
    case QEVMODEL_RODRIGUEZ2026:
      Ntotal=0.0;
      for (spec=0; spec<ns; spec++){
        N[spec]=_rhok(np[l],spec)/_m(spec); 
        Ntotal+=N[spec];
      }
      Tv=_Tv(np[l],gl);
      Te=_Te(np[l]);
      Tref=300.0e0; //gas temperature at which the experiments or BOLSIG+ results of mueN and E/N were obtained
      thetavib=_thetavk(specN2); //characteristic N2 vibration temperature in K
      // electron cooling-heating due to e-V interactions
      find_dTv_dU(np[l], gl, dTvdU);  
      find_dTe_dU(np[l], gl, dTedU);  
      f0 = 1.0-exp(-thetavib/Tv); 
      f1 = exp(-1.0*thetavib/Tv)*(1.0-exp(-thetavib/Tv)); 
      f2 = exp(-2.0*thetavib/Tv)*(1.0-exp(-thetavib/Tv)); 
      f3 = 1.0-f0-f1-f2;
      for (flux=0; flux<nf; flux++){
        // wrt Tv
        dQevdU[fluxev][flux]+=_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f0*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverNk_from_Te_chieminus(specN2, Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*exp(- thetavib/Tv)*(thetavib/sqr(Tv)) )*dTvdU[flux]
                            +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f1*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v1(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*exp(- thetavib/Tv)*(thetavib/sqr(Tv)) )*dTvdU[flux]
                            +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f2*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v2(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*exp(- thetavib/Tv)*(thetavib/sqr(Tv)) )*dTvdU[flux]
                            +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f3*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v3(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*exp(- thetavib/Tv)*(thetavib/sqr(Tv)) )*dTvdU[flux];
        // wrt Te
        dQevdU[fluxev][flux]+=_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f0*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverNk_from_Te_chieminus(specN2, Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*(-thetavib/sqr(Te))*exp(- thetavib/Tv) )*dTedU[flux]
                            +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f1*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v1(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*(-thetavib/sqr(Te))*exp(- thetavib/Tv) )*dTedU[flux]
                            +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f2*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v2(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*(-thetavib/sqr(Te))*exp(- thetavib/Tv) )*dTedU[flux]
                            +_zetav_from_Te_Rodriguez2025(_Te(np[l]))*N[speceminus]*f3*N[specN2]*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverN_from_Te_chieminus_N2v3(Te, N[speceminus]/Ntotal))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(-exp( thetavib/Te)*(-thetavib/sqr(Te))*exp(- thetavib/Tv) )*dTedU[flux];
      }
    break;  
#endif
  }
}

#endif //_FLUID_PLASMA


/* the diffusion coefficient of the vibrational energy, rho*Dv, such that the vibrational 
   energy conduction flux is -rho*Dv*grad(ev); rho*Dv=kappav/cv at constant composition */
double _nuv (np_t *np, long l, gl_t *gl) {
  double nuv;
  double Pr,T,cp;
  spec_t w;
  T=_T(np[l],gl);
  find_w(np[l],w);
  cp=_cp_from_w_T_equilibrium(w,T);
  Pr=_eta(np,l,gl)*cp/_kappa(np,l,gl);
  nuv=_wv_from_w(w)*(_eta(np,l,gl)/Pr);  
  return(nuv);
}


double _kappav (np_t *np, long l, gl_t *gl) {
  double kappav;
  double Pr,T,cp;
  spec_t w;
  long spec;
  vibdist_t vibdist;
  find_vibdist_from_np(np[l],gl,&vibdist);
  T=_T(np[l],gl);
  find_w(np[l],w);
  cp=_cp_from_w_T_equilibrium(w,T);
  Pr=_eta(np,l,gl)*cp/_kappa(np,l,gl);
  kappav=0.0;
  for (spec=0; spec<ns; spec++) if (_isvib(spec)) kappav+=w[spec]*(_eta(np,l,gl)/Pr)*vibdist.cvk[spec];  
  return(kappav);
}


#ifdef _FLUID_FAVREREYNOLDS


/* rho*Dv including the turbulent diffusion */
double _nuvstar (np_t *np, long l, gl_t *gl) {
  double nuvstar;
  double Pr,T,cp;
  spec_t w;
  T=_T(np[l],gl);
  find_w(np[l],w);
  cp=_cp_from_w_T_equilibrium(w,T);
  Pr=_eta(np,l,gl)*cp/_kappa(np,l,gl);
  nuvstar=_wv_from_w(w)*(_eta(np,l,gl)/Pr+_etat(np,l,gl)/gl->model.fluid.Prt);  
  return(nuvstar);
}


double _kappavstar (np_t *np, long l, gl_t *gl) {
  double kappavstar;
  double Pr,T,cp;
  spec_t w;
  long spec;
  vibdist_t vibdist;
  find_vibdist_from_np(np[l],gl,&vibdist);
  T=_T(np[l],gl);
  find_w(np[l],w);
  cp=_cp_from_w_T_equilibrium(w,T);
  Pr=_eta(np,l,gl)*cp/_kappa(np,l,gl);
  kappavstar=0.0;
  for (spec=0; spec<ns; spec++) if (_isvib(spec)) kappavstar+=w[spec]*(_eta(np,l,gl)/Pr+_etat(np,l,gl)/gl->model.fluid.Prt)*vibdist.cvk[spec];  
  return(kappavstar);
}

#endif //_FLUID_FAVREREYNOLDS

#endif //_FLUID_VIBMODEL
