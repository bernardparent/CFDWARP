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

/*
Ballistic electron model (birth-potential multigroup model) in argon:
B. Parent and B. Perry, "Birth-Potential Multigroup Fluid Model for Ballistic Electrons:
Breakdown and Cathode Sheath Regimes", arXiv:2610.02139, 2026 [parent2026].
The ionization by an external electron beam (Reaction 8) follows B. Parent, M. N. Shneider and
S. O. Macheret, "Sheath governing equations in computational weakly-ionized plasmadynamics",
J. Comput. Phys. 232:234-251, 2013 [parent2013].

Species 0 is the bulk electron group; species 1..nb-1 are the volumetric ballistic groups
ordered by energy (1 is the least energetic) and species nb is the group of the electrons
emitted from the cathode. The equation numbers below refer to parent2026.

The reactions of Table 1 are not in standard form: their products carry fractional weights
(f_j, g_j, chi_bk) that depend on the potential, on the electric field and on N_Ar.
The source terms and their Jacobian are hence hard-coded here rather than built on chem_share.c.
Each reaction is written as a rate of progress omega (1/(m3 s)) and its derivatives, and every
contribution to the Jacobian goes through add_to_dW_term(), which, as in chem_share.c, keeps only
the negative contributions unless TEST is defined; the one exception is the gain of the bulk in
Reaction 6 (see find_W_dW_parent2026()). A species that is not changed by a reaction
(Ar in the elastic and excitation collisions and in the thermal runaway) is left out of both sides.
*/

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/transport/_transport.h>
#include <model/share/chem_share.h>
#include <model/share/model_share.h>
#include <model/emfield/_emfield.h>
#include <model/fluid/_fluid.h>
#include <model/share/fluid_plasma_share.h>
#include <src/control.h>
#include <exm.h>
#include "hayashi.h"
#include "chem.h"

#ifndef _FLUID_PLASMA
#error The Ar_Plasma_Ballistic_20s chemical solver can not be compiled with a fluid module without plasma support.
#endif

void write_model_chem_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    BALLISTIC=TRUE;     {FALSE: no ballistic electrons; the model then collapses to the fluid LEA model}\n"
    "    phigauge=1000.0;    {V, potential mapped to the zero of the energy axis (Eq. 3)}\n"
    "    phicathode=0.0;     {V, potential of the cathode}\n"
    "    epsemit=1.0;        {eV, energy of the electrons emitted from the cathode}\n"
    "    Tborn=1600.0;       {K, temperature at which electrons are created in the volume}\n"
    "    zeta=0.8;           {straggling parameter of the runaway partition (Eq. 12)}\n"
    "    EcritfactReaction7=5.0; {ratio of the critical field of Reaction 7, Ecrit7, to Ecrit (Eq. 10); 1e99 turns Reaction 7 off}\n"
    "    chi0=0.5;           {share of the bulk electrons that can be re-collimated (Eq. 34)}\n"
    "    kturnmax=1e12;      {1/s, cap on the turning point rate coefficient (Eq. 32)}\n"
    "    Qbeamext=0.0;       {W/m3, power deposited by an external electron beam (Reaction 8)}\n"
    "    beamdir_x=1.0;      {direction along which the ballistic groups travel}\n"
    "    beamdir_y=0.0;\n"
    "    beamdir_z=0.0;\n"
    "    UNIFORMLEVELS=TRUE; {levels uniform between 0 and phigauge-phicathode}\n"
    "    {if UNIFORMLEVELS=FALSE, give the %d levels in eV in increasing order, \n"
    "     epslevel1=0.0; epslevel2=...; up to epslevel%d}\n"
    "  );\n"
  ,_CHEM_ACTIONNAME,nb,nb);
}


static void find_group_geometry(gl_t *gl){
  long j;
  double epsmax;
  gl_model_chem_t *chem;

  chem=&(gl->model.chem);
  /* windows of the volumetric groups, Eqs. (3)-(5) */
  for (j=1; j<nb; j++){
    chem->philo[j]=chem->phigauge-chem->epslevel[j+1];
    chem->phihi[j]=chem->phigauge-chem->epslevel[j];
    chem->phimid[j]=0.5*(chem->philo[j]+chem->phihi[j]);
    chem->Deltaphi[j]=chem->phihi[j]-chem->philo[j];
  }
  /* r_j is the width of the adjacent less energetic window over the one of group j (App. A);
     group 1 has no such window and r[1] is not used */
  chem->r[1]=1.0;
  for (j=2; j<nb; j++) chem->r[j]=chem->Deltaphi[j-1]/chem->Deltaphi[j];
  /* emission group, Eq. (6) */
  chem->phiemit=chem->phicathode-chem->epsemit;
  chem->philo[nb]=chem->phiemit;
  chem->phihi[nb]=chem->phiemit;
  chem->phimid[nb]=chem->phiemit;
  chem->Deltaphi[nb]=0.0;
  chem->r[nb]=0.0;
  chem->epsborn=1.5*kB*chem->Tborn/echarge;
  /* max_eps sum_k sigma_k Deltaeps_k, so that Ecrit=N_Ar*lossmax in V/m (Eq. 10) */
  chem->lossmax=_maxloss_Hayashi(&epsmax);
}


void read_model_chem_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long numvarsinit,j,dim;
  double norm;
  char varname[100];
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl_model_chem_t *chem;

  chem=&(gl->model.chem);
  if (strcmp(actionname,_CHEM_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);

    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_CHEM_ACTIONNAME);
    gl->MODEL_CHEM_READ=TRUE;

    action_original=codex->action;
    codex->action=&read_model_chem_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_bool_var_from_codex(codex,"BALLISTIC",&chem->BALLISTIC);
    find_double_var_from_codex(codex,"phigauge",&chem->phigauge);
    find_double_var_from_codex(codex,"phicathode",&chem->phicathode);
    if (chem->phigauge<=chem->phicathode) SOAP_fatal_error(codex,"phigauge must be greater than phicathode.");
    find_double_var_from_codex(codex,"epsemit",&chem->epsemit);
    if (chem->epsemit<=0.0) SOAP_fatal_error(codex,"epsemit must be positive: the density of the emission group at the cathode is its flux divided by its speed.");
    find_double_var_from_codex(codex,"Tborn",&chem->Tborn);
    if (chem->Tborn<0.0) SOAP_fatal_error(codex,"Tborn must not be negative.");
    find_double_var_from_codex(codex,"zeta",&chem->zeta);
    if (chem->zeta<=0.0) SOAP_fatal_error(codex,"zeta must be positive.");
    find_double_var_from_codex(codex,"EcritfactReaction7",&chem->EcritfactReaction7);
    if (chem->EcritfactReaction7<1.0) SOAP_fatal_error(codex,"EcritfactReaction7 must be at least 1.");
    find_double_var_from_codex(codex,"chi0",&chem->chi0);
    if (chem->chi0<0.0 || chem->chi0>1.0) SOAP_fatal_error(codex,"chi0 must be within 0 and 1.");
    find_double_var_from_codex(codex,"kturnmax",&chem->kturnmax);
    if (chem->kturnmax<=0.0) SOAP_fatal_error(codex,"kturnmax must be positive.");
    find_double_var_from_codex(codex,"Qbeamext",&chem->Qbeamext);
    if (chem->Qbeamext<0.0) SOAP_fatal_error(codex,"Qbeamext must not be negative.");
    find_double_var_from_codex(codex,"beamdir_x",&chem->beamdir[0]);
    find_double_var_from_codex(codex,"beamdir_y",&chem->beamdir[1]);
    find_double_var_from_codex(codex,"beamdir_z",&chem->beamdir[2]);
    norm=0.0;
    for (dim=0; dim<3; dim++) norm+=sqr(chem->beamdir[dim]);
    norm=sqrt(norm);
    if (norm==0.0) SOAP_fatal_error(codex,"beamdir_x, beamdir_y and beamdir_z can not all be zero.");
    for (dim=0; dim<3; dim++) chem->beamdir[dim]/=norm;
    find_bool_var_from_codex(codex,"UNIFORMLEVELS",&chem->UNIFORMLEVELS);
    chem->epslevel[0]=0.0;
    if (chem->UNIFORMLEVELS){
      for (j=1; j<=nb; j++) chem->epslevel[j]=(double)(j-1)/(double)(nb-1)*(chem->phigauge-chem->phicathode);
    } else {
      for (j=1; j<=nb; j++){
        sprintf(varname,"epslevel%ld",j);
        find_double_var_from_codex(codex,varname,&chem->epslevel[j]);
      }
      if (chem->epslevel[1]<0.0) SOAP_fatal_error(codex,"epslevel1 must not be negative.");
      for (j=2; j<=nb; j++){
        if (chem->epslevel[j]<=chem->epslevel[j-1])
          SOAP_fatal_error(codex,"The levels epslevel1..epslevel%d must be in increasing order.",nb);
      }
    }
    find_group_geometry(gl);

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}


/* destination weights g_j (j=0..nb, nonzero only for j=1..nb-1) of an electron born at phibirth,
   and their derivatives with respect to phibirth, Eqs. (14)-(16); a birth above the window of
   group 1 goes to group 1 and a birth below the window of group nb-1 goes to group nb-1 */
static void find_g(gl_t *gl, double phibirth, double *g, double *dgdphibirth){
  long j,jA;
  double s;
  for (j=0; j<=nb; j++){
    g[j]=0.0;
    dgdphibirth[j]=0.0;
  }
  if (phibirth>=gl->model.chem.phihi[1]) {
    g[1]=1.0;
  } else if (phibirth<gl->model.chem.philo[nb-1]) {
    g[nb-1]=1.0;
  } else {
    jA=1;
    while (jA<nb-1 && phibirth<gl->model.chem.philo[jA]) jA++;
    if (jA==nb-1) {
      g[nb-1]=1.0;
    } else {
      s=(gl->model.chem.phihi[jA]-phibirth)/gl->model.chem.Deltaphi[jA];
      g[jA]=1.0-s;
      g[jA+1]=s;
      dgdphibirth[jA]=1.0/gl->model.chem.Deltaphi[jA];
      dgdphibirth[jA+1]=-1.0/gl->model.chem.Deltaphi[jA];
    }
  }
}


/* runaway fraction, Eqs. (10)-(12), with the critical field Ecritfact*Ecrit, and its derivative
   with respect to N_Ar; En is the component of E along beamdir; without ballistic electrons, no
   electron runs away, and all of them are hence created in the bulk */
static void find_frun_Ecritfact(gl_t *gl, double En, double NAr, double Ecritfact, double *frun, double *dfrundNAr){
  double Ecrit,x;
  *frun=0.0;
  *dfrundNAr=0.0;
  Ecrit=Ecritfact*gl->model.chem.lossmax*NAr;
  if (gl->model.chem.BALLISTIC && NAr>0.0 && -En>Ecrit){
    x=(-En-Ecrit)/(gl->model.chem.zeta*Ecrit);
    *frun=1.0-exp(-x);
    *dfrundNAr=-exp(-x)*(-En)/(gl->model.chem.zeta*Ecrit*NAr);
  }
}


/* runaway fraction of an electron liberated at rest, Eqs. (10)-(12) */
void find_frun(gl_t *gl, double En, double NAr, double *frun, double *dfrundNAr){
  find_frun_Ecritfact(gl, En, NAr, 1.0, frun, dfrundNAr);
}


/* exp(f(ln Te)) with f the monotone spline through the control points (x,f), held constant
   outside of them, and its derivative with respect to Te */
static double _f_from_lnspline(long N, double *x, double *f, double Te, double *dfdTe){
  double lnTe,val,lo,hi;
  lnTe=log(Te);
  if (lnTe<=x[0] || lnTe>=x[N-1]){
    *dfdTe=0.0;
    if (lnTe<=x[0]) val=exp(f[0]); else val=exp(f[N-1]);
  } else {
    val=exp(EXM_f_from_monotonespline(N, x, f, lnTe));
    lo=max(x[0],lnTe-1e-6);
    hi=min(x[N-1],lnTe+1e-6);
    *dfdTe=val*(EXM_f_from_monotonespline(N, x, f, hi)-EXM_f_from_monotonespline(N, x, f, lo))/(hi-lo)/Te;
  }
  return(val);
}


/* electron-impact ionization rate coefficient k1 in m3/s, Table 2, and its derivative with respect to Te */
double _k1(double Te, double *dk1dTe){
  double lnTe_control[] = {
    10.42099, 10.56124, 10.59131, 10.65607, 10.70699, 10.7947,
    10.92818, 11.12563, 11.45201, 13.11916, 16.1181
  };
  double lnk1_control[] = {
    -76.65280, -62.60476, -58.59835, -48.07229, -43.27593, -39.10465,
    -36.06700, -34.02993, -32.4686, -29.9596, -29.25725
  };
  return(_f_from_lnspline(sizeof(lnTe_control)/sizeof(lnTe_control[0]), lnTe_control, lnk1_control, Te, dk1dTe));
}


/* reduced mobility mu_e*N of the bulk electrons in argon in 1/(m V s), given by the transport module,
   and its derivative with respect to Te */
double _mueN(gl_t *gl, double Te, double *dmueNdTe){
  double dTe;
  dTe=1e-6*Te;
  *dmueNdTe=(_mueNk_from_Te(gl,specAr,Te+dTe)-_mueNk_from_Te(gl,specAr,Te-dTe))/(2.0*dTe);
  return(_mueNk_from_Te(gl,specAr,Te));
}


/* mobility of the argon ions in m2/(V s), given by the transport module and corrected for the
   electric field with adjust_muk_for_Ek_effect(), as in the fluid module, and, if JACOBIAN is TRUE, 
   its derivatives with respect to NAr, NAr+ and T; Emag in V/m */
static double _muArplus(gl_t *gl, spec_t rhok, double T, double Te, double Emag, bool JACOBIAN,
                        double *dmudNAr, double *dmudNArp, double *dmudT){
  spec_t rhokcopy,nuk,dmudrhok;
  chargedspec_t muk,kappac;
  double eta,kappan,N,mu,mufield;
  long k;
  *dmudNAr=0.0;
  *dmudNArp=0.0;
  *dmudT=0.0;
  N=rhok[specAr]/_m(specAr)+rhok[specArplus]/_m(specArplus);
  if (N<=0.0) return(0.0);
  /* find_nuk_eta_kappak_muk() may rewrite the densities it is given */
  for (k=0; k<ns; k++) rhokcopy[k]=rhok[k];
  find_nuk_eta_kappak_muk(gl, rhokcopy, T, Te, nuk, &eta, &kappan, kappac, muk);
  mu=muk[specArplus];
  mufield=mu;
  if (Emag>0.0) adjust_muk_for_Ek_effect(specArplus, Emag, N, &mufield);
  if (JACOBIAN){
    if (mufield<mu){
      /* field-limited mobility, proportional to N^(-1/2) */
      *dmudNAr=-0.5*mufield/N;
      *dmudNArp=-0.5*mufield/N;
    } else {
      /* the derivatives of the mobility of the transport module without the field correction (Ek=0) */
      find_dmuk_from_rhok_Tk_Ek(gl, rhok, T, 0.0, specArplus, dmudT, dmudrhok);
      *dmudNAr=dmudrhok[specAr]*_m(specAr);
      *dmudNArp=dmudrhok[specArplus]*_m(specArplus);
    }
  }
  return(mufield);
}


/* ion-impact ionization rate coefficient k2 in m3/s (Table 1, without the f_j factor) and, if JACOBIAN is TRUE,
   its derivatives with respect to NAr, NAr+ and T; Emag in V/m */
void find_k2(gl_t *gl, spec_t rhok, double T, double Te, double Emag, bool JACOBIAN,
                    double *k2, double *dk2dNAr, double *dk2dNArp, double *dk2dT){
  double N,x,F,dFdN,mu,dmudNAr,dmudNArp,dmudT;
  *k2=0.0;
  *dk2dNAr=0.0;
  *dk2dNArp=0.0;
  *dk2dT=0.0;
  N=rhok[specAr]/_m(specAr)+rhok[specArplus]/_m(specArplus);
  if (Emag<=0.0 || N<=0.0) return;
  x=Emag/N/1e-21; /* E/N in Townsends */
  F=8.2e-22*exp(-19900.0/x)+7.1e-21*exp(-46000.0/x);
  mu=_muArplus(gl, rhok, T, Te, Emag, JACOBIAN, &dmudNAr, &dmudNArp, &dmudT);
  *k2=mu*Emag*F;
  if (JACOBIAN){
    /* dF/dN=dF/dx*dx/dN with dx/dN=-x/N */
    dFdN=-(8.2e-22*19900.0*exp(-19900.0/x)+7.1e-21*46000.0*exp(-46000.0/x))/sqr(x)*x/N;
    *dk2dNAr=Emag*(dmudNAr*F+mu*dFdN);
    *dk2dNArp=Emag*(dmudNArp*F+mu*dFdN);
    *dk2dT=Emag*dmudT*F;
  }
}


/* adds one contribution to an element of the Jacobian; as in chem_share.c,
   only the negative contributions are kept unless TEST is defined */
static inline void add_to_dW_term(double *dW, double term){
#ifdef TEST
  *dW+=term;
#else
  *dW+=min(0.0,term);
#endif
}


/* the following add the contributions of a rate of progress omega (in 1/(m3 s))
   to W (kg/(m3 s)), dWdrhok (1/s) and dWdTe (kg/(m3 s K)); nu is the number of
   particles of species spec formed (nu>0) or destroyed (nu<0) per event,
   domegadN is the derivative of nu*omega with respect to the number density of species col,
   and domegadTe its derivative with respect to Te */

#ifdef TEST
/* largest rate added to W, used to check the conservation of mass and charge */
static __thread double nuomegamax;
#endif

static inline void add_W(spec_t W, long spec, double nuomega){
  W[spec]+=_m(spec)*nuomega;
#ifdef TEST
  nuomegamax=max(nuomegamax,fabs(nuomega));
#endif
}

static inline void add_dW(spec2_t dWdrhok, long spec, long col, double dnuomegadN){
  add_to_dW_term(&(dWdrhok[spec][col]),_m(spec)*dnuomegadN/_m(col));
}

static inline void add_dWdTe(spec_t dWdTe, long spec, double dnuomegadTe){
  add_to_dW_term(&(dWdTe[spec]),_m(spec)*dnuomegadTe);
}

static inline void add_dWdT(spec_t dWdT, long spec, double dnuomegadT){
  add_to_dW_term(&(dWdT[spec]),_m(spec)*dnuomegadT);
}


/* source terms of Eqs. (36), (38) and of the argon ions (Table 1), and, if JACOBIAN is TRUE,
   their Jacobian with respect to rhok, T and Te */
static void find_W_dW_parent2026(np_t np, gl_t *gl, spec_t rhok, double T, double Te, bool JACOBIAN,
                                 spec_t W, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe){
  long k,s,j,b,dim,ch;
  spec_t N;
  EXM_vec3D_t E;
  double phi,Emag,En,NAr,NArp,Nheavy,frun,dfrundNAr,frun7,dfrun7dNAr;
  double f[nb+1],dfdNAr[nb+1],gloc[nb+1],dgloc[nb+1],g[nb+1],dgdphi[nb+1];
  double omega,k1,dk1dTe,k2,dk2dNAr,dk2dNArp,dk2dT,krun,dkrundNAr,dkrundTe,mueN,dmueNdTe,phibirth;
  double eps,V,sigma,chi,c,dcdNAr,nu,dnudNAr,kturn,k8,domegadNAr,domegadNArp;

  for (k=0; k<ns; k++) {
    W[k]=0.0;
    N[k]=rhok[k]/_m(k);
  }
#ifdef TEST
  nuomegamax=0.0;
#endif
  if (JACOBIAN){
    for (k=0; k<ns; k++){
      dWdT[k]=0.0;
      dWdTe[k]=0.0;
      for (s=0; s<ns; s++) dWdrhok[k][s]=0.0;
    }
  }

  NAr=N[specAr];
  NArp=N[specArplus];
  Nheavy=NAr+NArp;
  phi=_phi(np,gl);
  find_E(np,gl,E);
  Emag=0.0;
  En=0.0;
  for (dim=0; dim<3; dim++) {
    Emag+=sqr(E[dim]);
    En+=E[dim]*gl->model.chem.beamdir[dim];
  }
  Emag=sqrt(Emag);

  /* partition coefficients f_j of an electron liberated at rest, Eqs. (13)-(17) */
  find_frun(gl, En, NAr, &frun, &dfrundNAr);
  find_g(gl, phi, gloc, dgloc);
  for (j=0; j<=nb; j++){
    f[j]=frun*gloc[j];
    dfdNAr[j]=dfrundNAr*gloc[j];
  }
  f[0]=1.0-frun;
  dfdNAr[0]=-dfrundNAr;

  /* Reaction 1: e0 + Ar -> Ar+ + e0 + e_j, rate coefficient f_j*k1 */
  if (REACTION1){
    k1=_k1(Te,&dk1dTe);
    omega=k1*N[speceminus]*NAr;
    add_W(W,speceminus,-omega);
    add_W(W,specAr,-omega);
    add_W(W,specArplus,omega);
    add_W(W,speceminus,omega);
    for (j=0; j<nb; j++) add_W(W,j,f[j]*omega);
    if (JACOBIAN){
      /* reactant side */
      for (s=0; s<2; s++){
        if (s==0) k=speceminus; else k=specAr;
        add_dW(dWdrhok,k,speceminus,-k1*NAr);
        add_dW(dWdrhok,k,specAr,-k1*N[speceminus]);
        add_dWdTe(dWdTe,k,-dk1dTe*N[speceminus]*NAr);
      }
      /* product side */
      for (s=0; s<2; s++){
        if (s==0) k=specArplus; else k=speceminus;
        add_dW(dWdrhok,k,speceminus,k1*NAr);
        add_dW(dWdrhok,k,specAr,k1*N[speceminus]);
        add_dWdTe(dWdTe,k,dk1dTe*N[speceminus]*NAr);
      }
      for (j=0; j<nb; j++){
        add_dW(dWdrhok,j,speceminus,f[j]*k1*NAr);
        add_dW(dWdrhok,j,specAr,f[j]*k1*N[speceminus]);
        add_dW(dWdrhok,j,specAr,dfdNAr[j]*omega);
        add_dWdTe(dWdTe,j,f[j]*dk1dTe*N[speceminus]*NAr);
      }
    }
  }

  /* Reaction 2: Ar+ + Ar -> Ar+ + Ar+ + e_j, rate coefficient f_j*k2 */
  if (REACTION2){
    find_k2(gl, rhok, T, Te, Emag, JACOBIAN, &k2, &dk2dNAr, &dk2dNArp, &dk2dT);
    omega=k2*NArp*NAr;
    add_W(W,specArplus,-omega);
    add_W(W,specAr,-omega);
    add_W(W,specArplus,2.0*omega);
    for (j=0; j<nb; j++) add_W(W,j,f[j]*omega);
    if (JACOBIAN && Nheavy>0.0){
      /* reactant side */
      for (s=0; s<2; s++){
        if (s==0) k=specArplus; else k=specAr;
        add_dW(dWdrhok,k,specArplus,-k2*NAr);
        add_dW(dWdrhok,k,specAr,-k2*NArp);
        add_dW(dWdrhok,k,specArplus,-dk2dNArp*NArp*NAr);
        add_dW(dWdrhok,k,specAr,-dk2dNAr*NArp*NAr);
        add_dWdT(dWdT,k,-dk2dT*NArp*NAr);
      }
      /* product side */
      add_dW(dWdrhok,specArplus,specArplus,2.0*k2*NAr);
      add_dW(dWdrhok,specArplus,specAr,2.0*k2*NArp);
      add_dW(dWdrhok,specArplus,specArplus,2.0*dk2dNArp*NArp*NAr);
      add_dW(dWdrhok,specArplus,specAr,2.0*dk2dNAr*NArp*NAr);
      add_dWdT(dWdT,specArplus,2.0*dk2dT*NArp*NAr);
      for (j=0; j<nb; j++){
        add_dW(dWdrhok,j,specArplus,f[j]*k2*NAr);
        add_dW(dWdrhok,j,specAr,f[j]*k2*NArp);
        add_dW(dWdrhok,j,specArplus,f[j]*dk2dNArp*NArp*NAr);
        add_dW(dWdrhok,j,specAr,f[j]*dk2dNAr*NArp*NAr);
        add_dW(dWdrhok,j,specAr,dfdNAr[j]*omega);
        add_dWdT(dWdT,j,f[j]*dk2dT*NArp*NAr);
      }
    }
  }

  /* Reactions 3 to 6, for every ballistic group b */
  for (b=1; b<=nb; b++){
    eps=_eps_group(gl,b,phi);
    V=sqrt(2.0*eps*echarge/emass);

    /* Reaction 3: e_b + Ar -> e0 + Ar, rate coefficient f0*sigma_el*V_b */
    if (REACTION3 && V>0.0){
      c=f[0]*_sigma_el_Hayashi(eps)*V;
      dcdNAr=dfdNAr[0]*_sigma_el_Hayashi(eps)*V;
      omega=c*N[b]*NAr;
      add_W(W,b,-omega);
      add_W(W,speceminus,omega);
      if (JACOBIAN){
        add_dW(dWdrhok,b,b,-c*NAr);
        add_dW(dWdrhok,b,specAr,-c*N[b]);
        add_dW(dWdrhok,b,specAr,-dcdNAr*N[b]*NAr);
        add_dW(dWdrhok,speceminus,b,c*NAr);
        add_dW(dWdrhok,speceminus,specAr,c*N[b]);
        add_dW(dWdrhok,speceminus,specAr,dcdNAr*N[b]*NAr);
      }
    }

    /* Reactions 4 and 5: e_b + Ar -> e_j + Ar*_k and e_b + Ar -> Ar+ + e_j + e_j',
       rate coefficient sigma_k*V_b, the primary going to group j with the weight
       chi_bk*g_j(phibirth_bk)+(1-chi_bk)*f_j and the secondary with the weight f_j */
    for (ch=0; ch<HAYASHI_NUMINEL; ch++){
      if (ch==HAYASHI_INEL_IZ && !REACTION5) continue;
      if (ch!=HAYASHI_INEL_IZ && !REACTION4) continue;
      sigma=_sigma_inel_Hayashi(ch,eps);
      if (sigma<=0.0 || V<=0.0) continue;
      chi=sqrt(max(0.0,1.0-_Deltaeps_inel_Hayashi(ch)/eps));
      phibirth=_phibar(gl,b,phi)+_Deltaeps_inel_Hayashi(ch);
      find_g(gl, phibirth, g, dgdphi);
      c=sigma*V;
      omega=c*N[b]*NAr;
      add_W(W,b,-omega);
      if (ch==HAYASHI_INEL_IZ){
        add_W(W,specAr,-omega);
        add_W(W,specArplus,omega);
      }
      for (j=0; j<nb; j++){
        nu=chi*g[j]+(1.0-chi)*f[j];
        if (ch==HAYASHI_INEL_IZ) nu+=f[j];
        add_W(W,j,nu*omega);
      }
      if (JACOBIAN){
        add_dW(dWdrhok,b,b,-c*NAr);
        add_dW(dWdrhok,b,specAr,-c*N[b]);
        if (ch==HAYASHI_INEL_IZ){
          add_dW(dWdrhok,specAr,b,-c*NAr);
          add_dW(dWdrhok,specAr,specAr,-c*N[b]);
          add_dW(dWdrhok,specArplus,b,c*NAr);
          add_dW(dWdrhok,specArplus,specAr,c*N[b]);
        }
        for (j=0; j<nb; j++){
          nu=chi*g[j]+(1.0-chi)*f[j];
          dnudNAr=(1.0-chi)*dfdNAr[j];
          if (nu!=0.0){
            add_dW(dWdrhok,j,b,nu*c*NAr);
            add_dW(dWdrhok,j,specAr,nu*c*N[b]);
          }
          if (dnudNAr!=0.0) add_dW(dWdrhok,j,specAr,dnudNAr*omega);
          if (ch==HAYASHI_INEL_IZ){
            /* the secondary electron */
            if (f[j]!=0.0){
              add_dW(dWdrhok,j,b,f[j]*c*NAr);
              add_dW(dWdrhok,j,specAr,f[j]*c*N[b]);
            }
            if (dfdNAr[j]!=0.0) add_dW(dWdrhok,j,specAr,dfdNAr[j]*omega);
          }
        }
      }
    }

    /* Reaction 6: e_b -> e0, rate coefficient kturn_b of Eq. (32), capped by kturnmax;
       the field decelerates the group where En>0 since the electron charge is negative;
       a group at rest (eps<=0, level above the local potential) has no sink otherwise and is
       returned to the bulk at kturnmax whatever the field direction; the gain of the bulk is
       kept in the Jacobian although positive, so that the transfer is implicit on both sides:
       otherwise the bulk would receive kturn*dt times what the group loses in one iteration */
    if (REACTION6 && ((eps<gl->model.chem.epsborn && En>0.0) || eps<=0.0)){
      if (V>0.0) kturn=min(gl->model.chem.kturnmax,echarge*En/(emass*V));
        else kturn=gl->model.chem.kturnmax;
      omega=kturn*N[b];
      add_W(W,b,-omega);
      add_W(W,speceminus,omega);
      if (JACOBIAN){
        add_dW(dWdrhok,b,b,-kturn);
        dWdrhok[speceminus][b]+=_m(speceminus)*kturn/_m(b);
      }
    }
  }

  /* Reaction 7: e0 + Ar -> e_j + Ar, rate coefficient chi0*g_j(phibirth_0)*frun7*|Ce|/(me*mueN),
     with phibirth_0 from Eq. (35) and frun7 the runaway fraction with Ecrit7=EcritfactReaction7*Ecrit */
  find_frun_Ecritfact(gl, En, NAr, gl->model.chem.EcritfactReaction7, &frun7, &dfrun7dNAr);
  if (frun7>0.0){
    mueN=_mueN(gl,Te,&dmueNdTe);
    krun=gl->model.chem.chi0*frun7*echarge/(emass*mueN);
    dkrundNAr=gl->model.chem.chi0*dfrun7dNAr*echarge/(emass*mueN);
    phibirth=phi-1.5*kB*Te/echarge;
    find_g(gl, phibirth, g, dgdphi);
    omega=krun*N[speceminus]*NAr;
    add_W(W,speceminus,-omega);
    for (j=1; j<nb; j++) add_W(W,j,g[j]*omega);
    if (JACOBIAN){
      dkrundTe=-krun/mueN*dmueNdTe;
      add_dW(dWdrhok,speceminus,speceminus,-krun*NAr);
      add_dW(dWdrhok,speceminus,specAr,-krun*N[speceminus]);
      add_dW(dWdrhok,speceminus,specAr,-dkrundNAr*N[speceminus]*NAr);
      add_dWdTe(dWdTe,speceminus,-dkrundTe*N[speceminus]*NAr);
      for (j=1; j<nb; j++){
        if (g[j]!=0.0){
          add_dW(dWdrhok,j,speceminus,g[j]*krun*NAr);
          add_dW(dWdrhok,j,specAr,g[j]*krun*N[speceminus]);
          add_dW(dWdrhok,j,specAr,g[j]*dkrundNAr*N[speceminus]*NAr);
          add_dWdTe(dWdTe,j,g[j]*dkrundTe*N[speceminus]*NAr);
        }
        /* the birth potential depends on Te through Eq. (35) */
        if (dgdphi[j]!=0.0) add_dWdTe(dWdTe,j,-dgdphi[j]*1.5*kB/echarge*omega);
      }
    }
  }

  /* Reaction 8: Ar -> Ar+ + e_j, rate coefficient f_j*2.36e17*Qbeamext/N with N=NAr+NAr+ [parent2013];
     the electron is liberated at rest */
  if (REACTION8 && gl->model.chem.Qbeamext>0.0 && Nheavy>0.0){
    k8=2.36e17*gl->model.chem.Qbeamext/Nheavy;
    omega=k8*NAr;
    add_W(W,specAr,-omega);
    add_W(W,specArplus,omega);
    for (j=0; j<nb; j++) add_W(W,j,f[j]*omega);
    if (JACOBIAN){
      /* omega=2.36e17*Qbeamext*NAr/(NAr+NAr+) */
      domegadNAr=k8*NArp/Nheavy;
      domegadNArp=-k8*NAr/Nheavy;
      add_dW(dWdrhok,specAr,specAr,-domegadNAr);
      add_dW(dWdrhok,specAr,specArplus,-domegadNArp);
      add_dW(dWdrhok,specArplus,specAr,domegadNAr);
      add_dW(dWdrhok,specArplus,specArplus,domegadNArp);
      for (j=0; j<nb; j++){
        add_dW(dWdrhok,j,specAr,f[j]*domegadNAr);
        add_dW(dWdrhok,j,specArplus,f[j]*domegadNArp);
        add_dW(dWdrhok,j,specAr,dfdNAr[j]*omega);
      }
    }
  }

#ifdef TEST
  /* every reaction conserves the charge, to the roundoff error of the largest rate added to W, and the mass,
     to the precision of the molar masses in thermo (the one of Ar+ differs from the one of Ar minus the one of e-
     by 4e-8 of it) */
  {
    double summass,sumcharge;
    summass=0.0;
    sumcharge=0.0;
    for (k=0; k<ns; k++){
      summass+=W[k];
      sumcharge+=W[k]/_m(k)*_C(k);
    }
    if (fabs(summass)>1e-6*nuomegamax*_m(specAr) || fabs(sumcharge)>1e-10*nuomegamax*echarge)
      fatal_error("Mass or charge not conserved in find_W_dW_parent2026(): sum W=%E kg/m3s, sum W*C/m=%E C/m3s.",summass,sumcharge);
  }
#endif
}


void find_W ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  spec2_t dWdrhok;
  spec_t dWdT,dWdTe;
  find_W_dW_parent2026(np, gl, rhok, T, Te, FALSE, W, dWdrhok, dWdT, dWdTe);
}


void find_dW_dx ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv,
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k,s;
  spec_t W;
  for (k=0; k<ns; k++){
    dWdT[k]=0.0;
    dWdTe[k]=0.0;
    dWdTv[k]=0.0;
    dWdQbeam[k]=0.0;
    for (s=0; s<ns; s++) dWdrhok[k][s]=0.0;
  }
  find_W_dW_parent2026(np, gl, rhok, T, Te, TRUE, W, dWdrhok, dWdT, dWdTe);
}


/* energy source terms in W/m3 coming from the reactions of the ballistic electron model: Q[QB_BORN], Q[QB_BEAM] and Q[QB_RUN]
   are the terms Qborn, Qbeam and Qrun of the bulk electron energy equation, Eqs. (41)-(43), and Q[QB_THRESHOLD]
   is the energy the ballistic groups lose to the thresholds of the excitations and ionizations of Reactions 4 and 5;
   if JACOBIAN is TRUE, dQdrhok, dQdT and dQdTe are their derivatives with respect to rhok, T and Te, not clipped */
void find_Q_ballistic(np_t np, gl_t *gl, spec_t rhok, double T, double Te, bool JACOBIAN,
                       double *Q, spec_t *dQdrhok, double *dQdT, double *dQdTe){
  long k,b,dim,ch,n;
  spec_t N;
  spec_t dQdN[numQballistic];
  EXM_vec3D_t E;
  double phi,Emag,En,NAr,NArp,Nheavy,frun,dfrundNAr,f0,df0dNAr,frun7,dfrun7dNAr;
  double k1,dk1dTe,k2,dk2dNAr,dk2dNArp,dk2dT,eps,V,sigma,chi,c,dcdNAr,kturn,sumiz,R;
  double mueN,dmueNdTe,krun,dkrundNAr,dkrundTe,ee,deedTe,epsborn,thres,k8;

  for (n=0; n<numQballistic; n++){
    Q[n]=0.0;
    dQdT[n]=0.0;
    dQdTe[n]=0.0;
    for (k=0; k<ns; k++) dQdN[n][k]=0.0;
  }
  for (k=0; k<ns; k++) N[k]=rhok[k]/_m(k);
  NAr=N[specAr];
  NArp=N[specArplus];
  Nheavy=NAr+NArp;
  phi=_phi(np,gl);
  find_E(np,gl,E);
  Emag=0.0;
  En=0.0;
  for (dim=0; dim<3; dim++) {
    Emag+=sqr(E[dim]);
    En+=E[dim]*gl->model.chem.beamdir[dim];
  }
  Emag=sqrt(Emag);
  find_frun(gl, En, NAr, &frun, &dfrundNAr);
  f0=1.0-frun;
  df0dNAr=-dfrundNAr;
  epsborn=gl->model.chem.epsborn*echarge;

  /* electrons liberated at rest by Reactions 1, 2, 5 and 8: the bracket of Eq. (41), R in 1/(m3 s) */
  R=0.0;
  if (REACTION1){
    k1=_k1(Te,&dk1dTe);
    R+=k1*N[speceminus]*NAr;
    dQdN[QB_BORN][speceminus]+=f0*epsborn*k1*NAr;
    dQdN[QB_BORN][specAr]+=f0*epsborn*k1*N[speceminus];
    dQdTe[QB_BORN]+=f0*epsborn*dk1dTe*N[speceminus]*NAr;
  }
  if (REACTION2 && Nheavy>0.0){
    find_k2(gl, rhok, T, Te, Emag, JACOBIAN, &k2, &dk2dNAr, &dk2dNArp, &dk2dT);
    R+=k2*NArp*NAr;
    dQdN[QB_BORN][specArplus]+=f0*epsborn*(k2*NAr+dk2dNArp*NArp*NAr);
    dQdN[QB_BORN][specAr]+=f0*epsborn*(k2*NArp+dk2dNAr*NArp*NAr);
    dQdT[QB_BORN]+=f0*epsborn*dk2dT*NArp*NAr;
  }
  if (REACTION8 && gl->model.chem.Qbeamext>0.0 && Nheavy>0.0){
    k8=2.36e17*gl->model.chem.Qbeamext/Nheavy;
    R+=k8*NAr;
    dQdN[QB_BORN][specAr]+=f0*epsborn*k8*NArp/Nheavy;
    dQdN[QB_BORN][specArplus]-=f0*epsborn*k8*NAr/Nheavy;
  }

  for (b=1; b<=nb; b++){
    eps=_eps_group(gl,b,phi);
    V=sqrt(2.0*eps*echarge/emass);

    /* Reaction 3: the share f0 of the de-beamed electrons is thermalized in the bulk (Eq. (42)) */
    if (REACTION3 && V>0.0){
      c=f0*_sigma_el_Hayashi(eps)*eps*echarge*V;
      dcdNAr=df0dNAr*_sigma_el_Hayashi(eps)*eps*echarge*V;
      Q[QB_BEAM]+=c*N[b]*NAr;
      dQdN[QB_BEAM][b]+=c*NAr;
      dQdN[QB_BEAM][specAr]+=c*N[b]+dcdNAr*N[b]*NAr;
    }

    /* Reactions 4 and 5: the primary loses the threshold, and the share 1-chi of it gives
       what it has left to the bulk (Eq. (42)); the secondary of Reaction 5 is born at rest */
    sumiz=0.0;
    for (ch=0; ch<HAYASHI_NUMINEL; ch++){
      if (ch==HAYASHI_INEL_IZ && !REACTION5) continue;
      if (ch!=HAYASHI_INEL_IZ && !REACTION4) continue;
      sigma=_sigma_inel_Hayashi(ch,eps);
      if (sigma<=0.0 || V<=0.0) continue;
      chi=sqrt(max(0.0,1.0-_Deltaeps_inel_Hayashi(ch)/eps));
      c=(1.0-chi)*sigma*max(0.0,eps-_Deltaeps_inel_Hayashi(ch))*echarge*V;
      Q[QB_BEAM]+=c*N[b]*NAr;
      dQdN[QB_BEAM][b]+=c*NAr;
      dQdN[QB_BEAM][specAr]+=c*N[b];
      thres=sigma*_Deltaeps_inel_Hayashi(ch)*echarge*V;
      Q[QB_THRESHOLD]+=thres*N[b]*NAr;
      dQdN[QB_THRESHOLD][b]+=thres*NAr;
      dQdN[QB_THRESHOLD][specAr]+=thres*N[b];
      if (ch==HAYASHI_INEL_IZ) sumiz=sigma*V;
    }
    R+=sumiz*N[b]*NAr;
    dQdN[QB_BORN][b]+=f0*epsborn*sumiz*NAr;
    dQdN[QB_BORN][specAr]+=f0*epsborn*sumiz*N[b];

    /* Reaction 6: the group gives its kinetic energy to the bulk at the turning point (Eq. (42)) */
    if (REACTION6 && ((eps<gl->model.chem.epsborn && En>0.0) || eps<=0.0)){
      if (V>0.0) kturn=min(gl->model.chem.kturnmax,echarge*En/(emass*V));
        else kturn=gl->model.chem.kturnmax;
      Q[QB_BEAM]+=kturn*N[b]*eps*echarge;
      dQdN[QB_BEAM][b]+=kturn*eps*echarge;
    }
  }
  Q[QB_BORN]=f0*epsborn*R;
  dQdN[QB_BORN][specAr]+=df0dNAr*epsborn*R;

  /* Reaction 7: the runaway electrons take their thermal energy away from the bulk, Eq. (43);
     frun7 is the runaway fraction with Ecrit7=EcritfactReaction7*Ecrit */
  find_frun_Ecritfact(gl, En, NAr, gl->model.chem.EcritfactReaction7, &frun7, &dfrun7dNAr);
  if (frun7>0.0){
    mueN=_mueN(gl,Te,&dmueNdTe);
    krun=gl->model.chem.chi0*frun7*echarge/(emass*mueN);
    dkrundNAr=gl->model.chem.chi0*dfrun7dNAr*echarge/(emass*mueN);
    dkrundTe=-krun/mueN*dmueNdTe;
    ee=emass*_ee_from_Te(Te);
    deedTe=emass*_dee_dTe_from_Te(Te);
    Q[QB_RUN]=krun*N[speceminus]*NAr*ee;
    dQdN[QB_RUN][speceminus]+=krun*NAr*ee;
    dQdN[QB_RUN][specAr]+=(krun+dkrundNAr*NAr)*N[speceminus]*ee;
    dQdTe[QB_RUN]+=N[speceminus]*NAr*(dkrundTe*ee+krun*deedTe);
  }

  if (JACOBIAN){
    for (n=0; n<numQballistic; n++) for (k=0; k<ns; k++) dQdrhok[n][k]=dQdN[n][k]/_m(k);
  }
}


/* the energy exchanged by the bulk electrons in collisions is not given through Qei with this model:
   it is given by the terms Qborn, Qbeam and Qrun of Eqs. (41)-(43) and by the relaxation term of the TEMODEL
   of the fluid module */
void find_Qei(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  *Qei=0.0;
}


void find_dQei_dx(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;
}
