// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2000-2002, 2017 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "fluid.h"
#include "fluid_conv.h"
#include "fluid_source.h"
#include "fluid_bdry.h"
#include <src/control.h>
#include <src/common.h>
#include <model/thermo/_thermo.h>
#include <model/chem/_chem.h>
#include <model/metrics/_metrics.h>
#include <model/emfield/_emfield.h>
#include <model/_model.h>
#include <cycle/_cycle.h>
#include <model/share/fluid_share.h>


static bool _FLUIDPRIMMEM(np_t np){
  return(np.FLUIDPRIMMEM);
}


void write_model_fluid_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    TURBMODEL=TURBMODEL_KOMEGA2008;\n"
    "    DILATDISSIP=DILATDISSIP_WILCOX;\n"
    "    RAPCOMP=NO;\n"
    "    TURBSOURCE=YES;\n"
    "    REACTING=YES;\n"
    "    Prt=0.9e0;\n"
    "    Sct=1.0e0;\n"
    "    ADD_ETA_TO_ETAT_WITHIN_QK=NO;  {say NO for standard kinetic energy transport model}\n"
    "    kdiv=1.0e-50;    {m2/s2}\n"
    "    psidiv=0.1e-3;   {1/s}\n"
    "    Pmin=1.0e-2;    Pmax=9.9e99;   {Pa}\n"
    "    Tmin=1.0e1;     Tmax=6.0e3;    {K}\n"
    "    Twmin=Tmin;     Twmax=Tmax;    {K}\n"
    "    kmin=1.0e-10;   kmax=9.9e99;   {m2/s2}\n"
    "    psimin=1e-10;   psimax=9.9e99; {1/s}\n"
    "    wmin=1.0e-50;                  {min mass fraction allowed in the domain}\n"
#ifdef _2D
    "    AXISYMMETRIC=NO;\n"
#endif
    "  );\n"
  ,_FLUID_ACTIONNAME);

}


void write_cycle_fluid_template(FILE **controlfile){
  long spec,dim;
  wfprintf(*controlfile,
    "  %s(\n"
    "    xiverge=1e-3;{residual convergence threshold}\n"
    "    rhoref=1.0;  {reference density in kg/m3}\n"
    "    aref=300.0;  {reference sound speed in m/s}\n" 
    "    kref=1e4;    {reference turbulence kinetic energy in m2/s2}\n" 
    "    psiref=1e8;  {reference specific dissipation rate of the TKE in 1/s if for TURBMODEL_KOMEGA*\n"
    "                  reference dissipation rate of the TKE in m2/s3 for TURBMODEL_KEPSILON}\n" 
  ,_FLUID_ACTIONNAME);
  for (spec=0; spec<ns; spec++){
    wfprintf(*controlfile,
    "    Uref[%d]=rhoref;   \n",spec+1);
  }
  for (dim=0; dim<nd; dim++){
    wfprintf(*controlfile,
    "    Uref[%d]=rhoref*aref;   \n",ns+dim+1);
  }
  wfprintf(*controlfile,
    "    Uref[%d]=rhoref*aref*aref;  \n"
    "    Uref[%d]=rhoref*kref;  \n"
    "    Uref[%d]=rhoref*psiref;  \n"
    "  );\n",nd+ns+1,nd+ns+2,nd+ns+3);  
}


void read_model_fluid_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_fluid_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_FLUID_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);

    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_FLUID_ACTIONNAME);
    SOAP_add_int_to_vars(codex,"DILATDISSIP_NONE",DILATDISSIP_NONE);
    SOAP_add_int_to_vars(codex,"DILATDISSIP_WILCOX",DILATDISSIP_WILCOX);
    SOAP_add_int_to_vars(codex,"DILATDISSIP_SARKAR",DILATDISSIP_SARKAR);
    SOAP_add_int_to_vars(codex,"TURBMODEL_KEPSILON",TURBMODEL_KEPSILON); 
    SOAP_add_int_to_vars(codex,"TURBMODEL_KOMEGA1988",TURBMODEL_KOMEGA1988); 
    SOAP_add_int_to_vars(codex,"TURBMODEL_KOMEGA2008",TURBMODEL_KOMEGA2008); 

    gl->MODEL_FLUID_READ=TRUE;
    
    action_original=codex->action;
    codex->action=&read_model_fluid_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_double_var_from_codex(codex,"Pmin",&gl->model.fluid.Pmin);
    find_double_var_from_codex(codex,"Pmax",&gl->model.fluid.Pmax);
    find_double_var_from_codex(codex,"Tmin",&gl->model.fluid.Tmin);
    find_double_var_from_codex(codex,"Tmax",&gl->model.fluid.Tmax);
    find_double_var_from_codex(codex,"Twmin",&gl->model.fluid.Twmin);
    find_double_var_from_codex(codex,"Twmax",&gl->model.fluid.Twmax);
    find_double_var_from_codex(codex,"kmin",&gl->model.fluid.kmin);
    find_double_var_from_codex(codex,"kmax",&gl->model.fluid.kmax);
    find_double_var_from_codex(codex,"psimin",&gl->model.fluid.psimin);
    find_double_var_from_codex(codex,"psimax",&gl->model.fluid.psimax);
    find_double_var_from_codex(codex,"wmin",&gl->model.fluid.wmin);
    find_double_var_from_codex(codex,"kdiv",&gl->model.fluid.kdiv);
    find_double_var_from_codex(codex,"psidiv",&gl->model.fluid.psidiv);
    find_double_var_from_codex(codex,"Sct",&gl->model.fluid.Sct);
    find_double_var_from_codex(codex,"Prt",&gl->model.fluid.Prt);
    find_bool_var_from_codex(codex,"ADD_ETA_TO_ETAT_WITHIN_QK",&gl->model.fluid.ADD_ETA_TO_ETAT_WITHIN_QK);
    find_bool_var_from_codex(codex,"TURBSOURCE",&gl->model.fluid.TURBSOURCE);
    find_int_var_from_codex(codex,"TURBMODEL",&gl->model.fluid.TURBMODEL);
    if (gl->model.fluid.TURBMODEL!=TURBMODEL_KEPSILON && gl->model.fluid.TURBMODEL!=TURBMODEL_KOMEGA1988 && gl->model.fluid.TURBMODEL!=TURBMODEL_KOMEGA2008)
      SOAP_fatal_error(codex,"TURBMODEL must be set to either TURBMODEL_KEPSILON or TURBMODEL_KOMEGA1988 or TURBMODEL_KOMEGA2008.");
    find_bool_var_from_codex(codex,"RAPCOMP",&gl->model.fluid.RAPCOMP);
    find_bool_var_from_codex(codex,"REACTING",&gl->model.fluid.REACTING);
    find_int_var_from_codex(codex,"DILATDISSIP",&gl->model.fluid.DILATDISSIP);
    if (gl->model.fluid.DILATDISSIP!=DILATDISSIP_WILCOX && gl->model.fluid.DILATDISSIP!=DILATDISSIP_SARKAR && gl->model.fluid.DILATDISSIP!=DILATDISSIP_NONE )
      SOAP_fatal_error(codex,"DILATDISSIP must be set to either DILATDISSIP_NONE, DILATDISSIP_WILCOX, DILATDISSIP_SARKAR.");

#ifdef _2D
    find_bool_var_from_codex(codex,"AXISYMMETRIC",&gl->model.fluid.AXISYMMETRIC);
#endif
    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}


void read_cycle_fluid_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_cycle_fluid_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  char tmpstr[200];
  long flux;
  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_FLUID_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_FLUID_ACTIONNAME);
    gl->CYCLE_FLUID_READ=TRUE;
    
    action_original=codex->action;
    codex->action=&read_cycle_fluid_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_double_var_from_codex(codex,"xiverge",&gl->cycle.fluid.xiverge);

    for (flux=0; flux<nf; flux++){
      sprintf(tmpstr,"Uref[%ld]",flux+1);
      find_double_var_from_codex(codex,tmpstr,&gl->cycle.fluid.Uref[flux]);
    }

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}


double _rhok (np_t np, long spec) {
  double ret;
  ret=np.bs->U[spec];
  return(ret);
}


double _rho (np_t np) {
  double ret;
  long spec;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->rhomem;
  } else {
    ret=0.0;
    for (spec=0; spec<ns; spec++){
      ret+=_rhok(np,spec);
    }
  }
  return(ret);
}


double _w(np_t np, long spec) {
  double ret;
  assert_np(np,_rho(np)!=0.0e0);
  ret=np.bs->U[spec]/_rho(np);
/*  if (_rho(np)==0.0e0) {
    if (_FLUIDPRIMMEM(np)) wfprintf(stdout,"resumed\n");
    exit(1);
  }*/
  return(ret);
}


void find_w(np_t np, spec_t w) {
  long spec;
  double sum;
  sum=_rho(np);
  assert_np(np,sum!=0.0e0);
  for (spec=0; spec<ns; spec++){
    w[spec]=np.bs->U[spec]/sum;
  }
}


double _V(np_t np, long theta) {
  double ret;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->Vmem[theta];
  } else {
    ret=np.bs->U[ns+theta]/_rho(np);
  }
  return(ret);
}


double _V_from_U(np_t np, long theta) {
  double rho,ret;
  long spec;
  rho=0.0e0;
  for (spec=0; spec<ns; spec++) rho+=np.bs->U[spec];
  ret=np.bs->U[ns+theta]/rho;
  return(ret);
}


void find_V(np_t np, dim_t V) {
  long dim;
  double rho;
  if (_FLUIDPRIMMEM(np)){
    for (dim=0; dim<nd; dim++) V[dim]=np.wk->Vmem[dim];
  } else {
    rho=_rho(np);
    for (dim=0; dim<nd; dim++) V[dim]=np.bs->U[ns+dim]/rho;
  }
}


double _q(np_t np) {
  double ret;
  long dim;
  ret=0.0;
  for (dim=0; dim<nd; dim++) ret+=sqr(_V(np,dim));
  ret=sqrt(ret);
  return(ret);
}


double _a(np_t np, gl_t *gl) {
  double sum,ret,rho,k;
  long dim,spec;
  double dPdrhoetstar,htstar,P;
  spec_t dPdrhok,w;
  dim_t dPdrhoV;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->amem;
  } else {
    rho=_rho(np);
    P=_P(np, gl);
    k=_k(np);
    find_w(np,w);
    sum=0.0e0;
    for (dim=0; dim<nd; dim++){
      sum=sum+sqr(_V(np,dim));
    }
    find_dP_dx(np, gl, &dPdrhoetstar, dPdrhok, dPdrhoV);
    htstar=np.bs->U[fluxet]/rho+2.0e0/3.0e0*k+P/rho;
//  np.wk->amem=2.0e0/3.0e0*k+dPdrhoetstar*(htstar-k-sum);
    ret=2.0e0/3.0e0*k+dPdrhoetstar*(htstar-sum-k);
    for (spec=0; spec<ns; spec++){
      ret+=dPdrhok[spec]*w[spec];
    }
    assert_np(np,ret>=0.0e0);
    ret=sqrt(max(0.0e0,ret));
  }
  return(ret);
}


double _athermo(np_t np, gl_t *gl) {
  double ret,rho,T;
  long spec;
  spec_t w;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->athermomem;
  } else {
    rho=_rho(np);
    T=_T(np,gl);
    for (spec=0; spec<ns; spec++){
      w[spec]=np.bs->U[spec]/rho;
    }
    ret=_a_from_w_T(w,T);
    assert_np(np,ret>0.0e0);
  }
  return(ret);
}


double _P (np_t np, gl_t *gl) {
  double ret,rho;
  spec_t w;
  long spec;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->Pmem;
  } else {
    rho=_rho(np);
    for (spec=0; spec<ns; spec++){
      w[spec]=np.bs->U[spec]/rho;
    }     
    ret=_P_from_w_rho_T(w,rho,_T(np,gl));
  }
  return(ret);
}


double _T (np_t np, gl_t *gl) {
  double ekin,eint,ret,rho;
  long dim;
  spec_t w;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->Tmem;
  } else {
    rho=_rho(np);
    find_w(np,w);
    assert(rho!=0.0e0);
    ekin=0.0e0;
    for (dim=0; dim<nd; dim++){
      ekin=ekin+0.5e0*sqr(np.bs->U[ns+dim]/rho);
    }
    eint=np.bs->U[fluxet]/rho /* total energy */
      -ekin  /* kinetic energy */
      -np.bs->U[fluxtke]/rho /* k */
      ; 
    ret=_T_from_w_e(w,eint);
  }
  return(ret);
}


double _k (np_t np){
  double ret;
  assert_np(np,_rho(np)!=0.0e0);
  ret=np.bs->U[fluxtke]/_rho(np);
  return(ret);
}


double _psi (np_t np){
  double ret;
  assert_np(np,_rho(np)!=0.0e0);
  ret=np.bs->U[fluxpsi]/_rho(np);
  return(ret);
}


double _etstar (np_t np){
  double ret;
  assert_np(np,_rho(np)!=0.0e0);
  ret=np.bs->U[fluxet]/_rho(np);
  return(ret);
}




double _Pstar (np_t np, gl_t *gl) {
  double ret;
  ret=_P(np,gl)+2.0e0/3.0e0*_rho(np)*_k(np);
  return(ret);
}


double _etat(np_t *np, long l, gl_t *gl) {
  double etat;
  if (_FLUIDPRIMMEM(np[l])){
    etat=np[l].wk->etat;
  } else {
    etat=_etat_mem(np,l,gl);
  }
  return(etat);
}


double _kappa(np_t np, gl_t *gl) {
  double T,cp,ret,eta,kappa;
  spec_t w,nu;

  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->kappamem;
  } else {
    T=_T(np,gl);
    find_w(np,w);
    cp=_cp_from_w_T(w,T);
    find_nuk_eta_kappa(w, _rho(np), T, nu, &eta, &kappa);
    ret=cp*_eta(np,gl)/( /* Prandtl number */(_eta(np,gl))/(kappa)*cp );

  }

  return(ret);
}


double _eta(np_t np, gl_t *gl) {
  double ret;
  spec_t w,nu;
  double eta,kappa;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->etamem;
  } else { 
    find_w(np,w); 
    find_nuk_eta_kappa(w, _rho(np), _T(np,gl), nu, &eta, &kappa);
    ret=eta;
  }
  return(ret);
}


double _nu(np_t np, gl_t *gl, long spec) {
  double ret,eta,kappa;
  spec_t nu,w;
  if (_FLUIDPRIMMEM(np)){
    ret=np.wk->numem[spec];
  } else {
    find_w(np,w); 
    find_nuk_eta_kappa(w, _rho(np), _T(np,gl), nu, &eta, &kappa);
    ret=nu[spec];  
  }
  return(ret);
}


void find_dP_dx(np_t np, gl_t *gl, double *dPdrhoetstar, spec_t dPdrhok, dim_t dPdrhoV){
  double k,T,etstar,dedP,sum,rho;
  dim_t V;
  spec_t rhok,dedrhok;
  long spec,dim;

  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=_V(np,dim);
    sum=sum+sqr(V[dim]);
  }
  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rhok[spec]=_rhok(np,spec);
    rho=rho+rhok[spec];
  }
  T=_T(np,gl);
  etstar=_etstar(np);
  k=_k(np);
  dedP=_de_dP_at_constant_rho(rhok,T);
  find_de_drhok_at_constant_P(rhok,T,dedrhok);

  /* ---- dPdrhoetstar ------*/
  assert_np(np,dedP*rho!=0.0e0);
  *dPdrhoetstar=1.0e0/(rho*dedP);

  /* ---- dPdrhok ------*/
  assert_np(np,rho!=0.0e0);
  assert_np(np,dedP!=0.0e0);
  for (spec=0; spec<ns; spec++){
      dPdrhok[spec]=-(*dPdrhoetstar)*(etstar-sum-k
                                 +dedrhok[spec]*rho);
  }

  /* ---- dPdrhou ------*/
  for (dim=0; dim<nd; dim++){
    dPdrhoV[dim]=-V[dim]*(*dPdrhoetstar);
  }
}


void find_dT_dx(np_t np, gl_t *gl, double *dTdrhoetstar, spec_t dTdrhok, dim_t dTdrhoV){
  double k,T,etstar,dedT,sum,rho;
  dim_t V;
  spec_t rhok,dedrhok;
  long spec,dim;

  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=_V(np,dim);
    sum=sum+sqr(V[dim]);
  }
  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rhok[spec]=_rhok(np,spec);
    rho=rho+rhok[spec];
  }
  T=_T(np,gl);
  etstar=_etstar(np);
  k=_k(np);
  dedT=_de_dT_at_constant_rho(rhok,T);
  find_de_drhok_at_constant_T(rhok,T,dedrhok);

  /* ---- dTdrhoetstar ------*/
  assert_np(np,dedT*rho!=0.0e0);
  *dTdrhoetstar=1.0e0/(rho*dedT);

  /* ---- dPdrhok ------*/
  assert_np(np,rho!=0.0e0);
  assert_np(np,dedT!=0.0e0);
  /*for (spec=0; spec<ns; spec++){
    dTdrhok[spec]=(-(etstar-k)/rho+sum/rho-dedrhok[spec])/dedT;
  }*/
  for (spec=0; spec<ns; spec++){
      dTdrhok[spec]=-(etstar-sum-k
                                 +dedrhok[spec]*rho)/dedT/rho;
  }

  /* ---- dPdrhou ------*/
  for (dim=0; dim<nd; dim++){
    dTdrhoV[dim]=-V[dim]*(*dTdrhoetstar);
  }

}


void find_dT_dU(np_t np, gl_t *gl, flux_t dTdU){
  long flux,dim,spec;
  double dTdrhoetstar;
  dim_t dTdrhoV;
  spec_t dTdrhok;
  /* double T1,T2; */
  for (flux=0; flux<nf; flux++) dTdU[flux]=0.0;
  find_dT_dx(np, gl, &dTdrhoetstar, dTdrhok, dTdrhoV);
  for (spec=0; spec<ns; spec++) dTdU[spec]=dTdrhok[spec];
  for (dim=0; dim<nd; dim++) dTdU[ns+dim]=dTdrhoV[dim];
  dTdU[fluxet]=dTdrhoetstar;
  dTdU[fluxtke]=-dTdrhoetstar;

  /*
 //testing!
  for (flux=0; flux<nf; flux++){
    wfprintf(stdout,"dTdU=%15.10E  ",dTdU[flux]);
    T1=_T(np);
    np.bs->U[flux]+=1.0e-6;
    find_prim_fluid(np,gl);
    T2=_T(np);
    np.bs->U[flux]-=1.0e-6;
    find_prim_fluid(np,gl);
    wfprintf(stdout,"dTdU=%15.10E\n",(T2-T1)/1.0e-6);
  }
 wfprintf(stdout,"\n\n");
*/
}


void find_Ustar(np_t np, gl_t *gl, flux_t Ustar){
  double Omega;
  long flux;

  Omega=_Omega(np,gl);
  for (flux=0; flux<nf; flux++)
    Ustar[flux]=np.bs->U[flux]*Omega;
}


void find_G(np_t np, gl_t *gl, flux_t G){
  long spec,dim;
  spec_t w;
  find_w(np,w);
  for (spec=0; spec<ns; spec++) G[spec]=w[spec];
  for (dim=0; dim<nd; dim++) G[ns+dim]=_V(np,dim);
  G[fluxet]=_T(np,gl);
  G[fluxtke]=_k(np);
  G[fluxpsi]=_psi(np);
}


void find_Kstar_interface(np_t *np, gl_t *gl, long lL, long lR, metrics_t metrics, long theta, long vartheta, sqmat_t K){
  long dim,dim2,spec,row,col;
  double Omega,alpha,sum;
  double beta[nd][nd];
  
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++) K[row][col]=0.0e0;
  }
  
    Omega=metrics.Omega;
    sum=0.0e0;
    for (dim=0; dim<nd; dim++){
      sum=sum+metrics.X2[theta][dim]*metrics.X2[vartheta][dim];
    }
    alpha=Omega*sum;
    for (row=0; row<nd; row++){
      beta[row][row]=Omega*(sum+1.0e0/3.0e0*metrics.X2[theta][row]
                                           *metrics.X2[vartheta][row]);
      for (col=0; col<nd; col++){
        if (col!=row) {
          beta[row][col]=Omega*(metrics.X2[vartheta][row]*metrics.X2[theta][col]
                -2.0e0/3.0e0*metrics.X2[theta][row]*metrics.X2[vartheta][col]);
        }
      }
    }
    for (row=0; row<nf; row++){
      for (col=0; col<nf; col++){
        K[row][col]=0.0e0;
      }
    }
    for (spec=0; spec<ns; spec++){
      K[spec][spec]=avg(_nustar(np,lL,gl,spec),_nustar(np,lR,gl,spec))*alpha;
      K[fluxet][spec]=avg(_nustar(np,lL,gl,spec),_nustar(np,lR,gl,spec))*_hk_from_T(spec,avg(_T(np[lL],gl),_T(np[lR],gl)))*alpha;
    }
    
    for (dim=0; dim<nd; dim++){
      for (dim2=0; dim2<nd; dim2++){
        K[ns+dim][ns+dim2]=avg(_etastar(np,lL,gl),_etastar(np,lR,gl))*beta[dim][dim2];
        K[fluxet][ns+dim2]=K[fluxet][ns+dim2]+
           avg(_etastar(np,lL,gl),_etastar(np,lR,gl))*beta[dim][dim2]*avg(_V(np[lL],dim),_V(np[lR],dim));
      }
    }
    
    K[fluxet][fluxet]=avg(_kappastar(np,lL,gl),_kappastar(np,lR,gl))*alpha;
    K[fluxet][fluxtke]=avg(_etakstar(np,lL,gl),_etakstar(np,lR,gl))*alpha;
    K[fluxtke][fluxtke]=avg(_etakstar(np,lL,gl),_etakstar(np,lR,gl))*alpha;
    K[fluxpsi][fluxpsi]=avg(_etapsistar(np,lL,gl),_etapsistar(np,lR,gl))*alpha; 
}


void find_dG_dUstar(np_t np, gl_t *gl, sqmat_t B){
  long dim,spec,row,col;
  spec_t w,dTdrhok;
  dim_t dTdrhoV;
  double k,psi,Omega,rho,dTdrhoetstar;
  Omega=_Omega(np,gl);
  rho=_rho(np);
  k=_k(np);
  psi=_psi(np);

  find_w(np,w);
  find_dT_dx(np, gl, &dTdrhoetstar, dTdrhok, dTdrhoV);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=0.0e0;
    }
  }
  for (row=0; row<ns; row++){
    for (col=0; col<ns; col++){
      assert_np(np,rho!=0.0e0);
      assert_np(np,Omega!=0.0e0);
      B[row][col]=-w[row]/(Omega*rho);
    }
  }
  for (spec=0; spec<ns; spec++){
    assert_np(np,rho!=0.0e0);
    assert_np(np,Omega!=0.0e0);
    B[spec][spec]=B[spec][spec]+1.0e0/(Omega*rho);
    for (dim=0; dim<nd; dim++){
      B[ns+dim][spec]=-_V(np,dim)/(Omega*rho);
    }
    B[fluxet][spec]=dTdrhok[spec]/(Omega);
    B[fluxtke][spec]=-k/(Omega*rho);
    B[fluxpsi][spec]=-psi/(Omega*rho);
  }
  for (dim=0; dim<nd; dim++){
    B[ns+dim][ns+dim]=1.0e0/(Omega*rho);
    B[fluxet][ns+dim]=dTdrhoV[dim]/(Omega);
  }
  B[fluxet][fluxet]=dTdrhoetstar/(Omega);
  B[fluxet][fluxtke]=-dTdrhoetstar/(Omega);
  B[fluxtke][fluxtke]=1.0e0/Omega/rho;
  B[fluxpsi][fluxpsi]=1.0e0/Omega/rho;
}


void reformat_T(gl_t *gl, double *T, char *suffix, bool *flag){
  if (*T<gl->model.fluid.Tmin) {
    add_to_clipped_variables2(gl,"Tmin",suffix);
    *T=gl->model.fluid.Tmin;
    *flag=TRUE;
  }
  if (*T>gl->model.fluid.Tmax) {
    add_to_clipped_variables2(gl,"Tmax",suffix);
    *T=gl->model.fluid.Tmax;
    *flag=TRUE;
  }
}


void reformat_P(gl_t *gl, double *P, char *suffix, bool *flag){
  if (*P<gl->model.fluid.Pmin) {
    add_to_clipped_variables2(gl,"Pmin",suffix);
    *P=gl->model.fluid.Pmin;
    *flag=TRUE;
  }
  if (*P>gl->model.fluid.Pmax) {
    add_to_clipped_variables2(gl,"Pmax",suffix);
    *P=gl->model.fluid.Pmax;
    *flag=TRUE;
  }
}


void reformat_w(gl_t *gl, spec_t w, char *suffix, bool *flag){
  long spec;
  double sum;
  char *speciesname;
  sum=0.0e0;
  for (spec=0; spec<ns; spec++){
    if (w[spec]<gl->model.fluid.wmin) {
      speciesname=(char *)malloc(sizeof(char));
      find_species_name(spec, &speciesname);
      SOAP_strins("wmin_", &speciesname,  0);
      add_to_clipped_variables2(gl,speciesname,suffix);
      w[spec]=gl->model.fluid.wmin;
      free(speciesname);
      *flag=TRUE;
    }
    sum=sum+w[spec];
  }
  assert(sum!=0.0e0);
  for (spec=0; spec<ns; spec++) w[spec]=w[spec]/sum;
}


void reformat_rhok(gl_t *gl, spec_t rhok, char *suffix, bool *flag){
  long spec;
  double rho;
  double rhomin,rhomax;
  spec_t w;

  rho=0.0e0;
  for (spec=0; spec<ns; spec++) rho=rho+rhok[spec];
  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
  reformat_w(gl,w,suffix,flag);
  rhomin=gl->model.fluid.Pmin/(_R(w)*gl->model.fluid.Tmax);
  rhomax=gl->model.fluid.Pmax/(_R(w)*gl->model.fluid.Tmin);
  if (rho<rhomin) {
    add_to_clipped_variables2(gl,"rhomin",suffix);
    rho=rhomin;
    *flag=TRUE;
  }
  if (rho>rhomax) {
    add_to_clipped_variables2(gl,"rhomax",suffix);
    rho=rhomax;
    *flag=TRUE;
  }
  for (spec=0; spec<ns; spec++) rhok[spec]=rho*w[spec];
}


void reformat_k_psi(gl_t *gl, double *k, double *psi, char *suffix, bool *flag){

  if (*k<gl->model.fluid.kmin){
    add_to_clipped_variables2(gl,"kmin",suffix);
    *k=gl->model.fluid.kmin;
    *flag=TRUE;
  }
  if (*k>gl->model.fluid.kmax){
    add_to_clipped_variables2(gl,"kmax",suffix);
    *k=gl->model.fluid.kmax;
    *flag=TRUE;
  }
  if (*psi<gl->model.fluid.psimin){
    add_to_clipped_variables2(gl,"psimin",suffix);
    *psi=gl->model.fluid.psimin;
    *flag=TRUE;
  }
  if (*psi>gl->model.fluid.psimax){
    add_to_clipped_variables2(gl,"psimax",suffix);
    *psi=gl->model.fluid.psimax;
    *flag=TRUE;
  }
}


void find_prim_fluid_mem(np_t *np, long l, gl_t *gl, double P, double T){
  spec_t w;
  double sum1,sum,rho,k,eta,psitilde;
  spec_t dPdrhok;
  dim_t dPdrhoV;
  double dPdrhoetstar,htstar;
  long dim,spec;

  assert_np(np[l],is_node_resumed(np[l]));
  /*reformat_T(&T,&ref_flag);
  reformat_P(&P,&ref_flag); */
  np[l].wk->Tmem=T;
  np[l].wk->Pmem=P;
  sum1=0.0e0;
  for (spec=0; spec<ns; spec++){
    sum1=sum1+np[l].bs->U[spec];
  }
  assert_np(np[l],sum1!=0.0e0);
  for (spec=0; spec<ns; spec++){
    w[spec]=np[l].bs->U[spec]/sum1;
  }
  /* reformat_w(w,&ref_flag); */
  np[l].wk->athermomem=_a_from_w_T(w,T);
  assert_np(np[l],np[l].wk->athermomem>0.0e0);
  rho=_rho_from_w_P_T(w,P,T);
  np[l].wk->rhomem=rho;
  assert_np(np[l],np[l].wk->rhomem>0.0e0);
  /* find Vmem */
  assert_np(np[l],rho!=0.0e0);
  for (dim=0; dim<nd; dim++){
    np[l].wk->Vmem[dim]=np[l].bs->U[ns+dim]/rho;
  }
  find_nuk_eta_kappa(w, rho, T,  np[l].wk->numem, &(np[l].wk->etamem), &(np[l].wk->kappamem));
  np[l].wk->kappamem=_cp_from_w_T(w,T)*np[l].wk->etamem/( /* Prandtl number */(np[l].wk->etamem)/(np[l].wk->kappamem)*(_cp_from_w_T(w,T)) );

  k=_k(np[l]);
  psitilde=_psitilde(np[l],gl);
  eta=np[l].wk->etamem;
  assert_np(np[l],psitilde!=0.0e0);

  /* find amem */
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    sum=sum+sqr(_V(np[l],dim));
  }
  find_dP_dx(np[l], gl, &dPdrhoetstar, dPdrhok, dPdrhoV);
  htstar=np[l].bs->U[fluxet]/rho+2.0e0/3.0e0*k+P/rho;

  np[l].wk->amem=2.0e0/3.0e0*k+dPdrhoetstar*(htstar-k-sum);
  for (spec=0; spec<ns; spec++){
    np[l].wk->amem+=dPdrhok[spec]*w[spec];
  }
  assert_np(np[l],np[l].wk->amem>=0.0e0);
  np[l].wk->amem=sqrt(max(0.0e0,np[l].wk->amem));
  np[l].wk->etat=_etat_from_rho_eta_k_psitilde(np, l, gl, rho, eta, k, psitilde);
  np[l].FLUIDPRIMMEM=TRUE;
}


void find_U(np_t *np, long l, gl_t *gl, spec_t rhok, dim_t V, double T, double k, double psi){
  double rho,emix,sum;
  long dim,spec;
  spec_t w;
  bool ref_flag;

  ref_flag=FALSE;
  reformat_rhok(gl,rhok,"",&ref_flag);
  reformat_T(gl,&T,"",&ref_flag);
  reformat_k_psi(gl,&k,&psi,"",&ref_flag);

  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rho=rho+rhok[spec];
    np[l].bs->U[spec]=rhok[spec];
  }
  assert_np(np[l],rho!=0.0e0);
  for (spec=0; spec<ns; spec++){
    w[spec]=rhok[spec]/rho;
  }
  emix=_e_from_w_T(w,T);
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    np[l].bs->U[ns+dim]=rho*V[dim];
    sum=sum+sqr(V[dim]);
  }
  np[l].bs->U[fluxet]=rho*(0.5e0*(sum)+k+emix);
  np[l].bs->U[fluxtke]=rho*k;
  np[l].bs->U[fluxpsi]=rho*psi;
  if (is_node_resumed(np[l])) {
    np[l].wk->Tmem=T;
    find_prim_fluid(np,l,gl);
  }
}


void find_U_2(np_t *np, long l, gl_t *gl, spec_t w, dim_t V, double P, double T, double k, double psi){
  double sum,rho,emix;
  long spec,dim;
  bool ref_flag;

  ref_flag=FALSE;
  /* reformat_w(gl,w,&ref_flag); */
  reformat_T(gl,&T,"",&ref_flag);
  reformat_P(gl,&P,"",&ref_flag);
  reformat_k_psi(gl,&k,&psi,"",&ref_flag);

  emix=_e_from_w_T(w,T);
  rho=_rho_from_w_P_T(w,P,T);
  for (spec=0; spec<ns; spec++){
    np[l].bs->U[spec]=rho*w[spec];
  }
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    np[l].bs->U[ns+dim]=rho*V[dim];
    sum=sum+sqr(V[dim]);
  }
  np[l].bs->U[fluxet]=rho*(0.5e0*(sum)+emix+k);
  np[l].bs->U[fluxtke]=rho*k;
  np[l].bs->U[fluxpsi]=rho*psi;
  if (is_node_resumed(np[l])){
    np[l].wk->Tmem=T;
    find_prim_fluid_mem(np, l, gl, P, T);
  }
}


void find_U_3(np_t *np, long l, gl_t *gl, spec_t rhok, dim_t V, double Pstar, double k, double psi){
  double T,rho,P;
  bool ref_flag;
  long spec;

  rho=0.0e0;
  for (spec=0; spec<ns; spec++){
    rho=rho+rhok[spec];
  }

  reformat_P(gl,&Pstar,"",&ref_flag);
  P=Pstar-2.0e0/3.0e0*rho*k;
  assert_np(np[l],P>0.0e0);
  ref_flag=FALSE;
  reformat_rhok(gl,rhok,"",&ref_flag);
  reformat_k_psi(gl,&k,&psi,"",&ref_flag);

  T=_T_from_rhok_P(rhok,P);
  find_U(np, l, gl, rhok, V, T, k, psi);
}


void find_prim_fluid(np_t *np, long l, gl_t *gl){
  long dim,spec;
  double rho,P,T,sum,k,psi,e;
  spec_t w;
  dim_t V;
  bool ref_flag;
  rho=0.0e0;
  for (spec=0; spec<ns; spec++) rho+=np[l].bs->U[spec];
  for (spec=0; spec<ns; spec++){
    if (np[l].bs->U[spec]/rho<gl->model.fluid.wmin){
      np[l].bs->U[spec]=rho*gl->model.fluid.wmin;
      add_to_clipped_variables(gl,"wmin");
    }
  }
  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=np[l].bs->U[spec];
  assert_np(np[l],rho!=0.0e0);
  for (spec=0; spec<ns; spec++){
    w[spec]=np[l].bs->U[spec]/rho;
  }
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=np[l].bs->U[ns+dim]/rho;
    sum=sum+sqr(V[dim]);
  }
  assert_np(np[l],rho!=0.0e0);
  k=np[l].bs->U[fluxtke]/rho;
  psi=np[l].bs->U[fluxpsi]/rho;
  e=np[l].bs->U[fluxet]/rho-0.5e0*sum-k;
  /* clip e and adjust U*/
  if (e<_e_from_w_T(w, gl->model.fluid.Tmin)){
    e=_e_from_w_T(w, gl->model.fluid.Tmin); 
    np[l].bs->U[fluxet]=rho*(e+0.5e0*sum+k);
    add_to_clipped_variables(gl,"Tmin");
  }
  if (e>_e_from_w_T(w, gl->model.fluid.Tmax)){
    e=_e_from_w_T(w, gl->model.fluid.Tmax); 
    np[l].bs->U[fluxet]=rho*(e+0.5e0*sum+k);
    add_to_clipped_variables(gl,"Tmax");
  }
  T=_T_from_w_e(w,e);
  np[l].wk->Tmem=T;
  P=_P_from_w_rho_T(w,rho,T);

  /*
  ret=2.0e0/3.0e0*k+dPdrhoetstar*(htstar-sum-k);
    for (spec=0; spec<ns; spec++){
      ret+=dPdrhok[spec]*w[spec];
    }
*/ 
  ref_flag=FALSE;
  reformat_k_psi(gl,&k,&psi,"",&ref_flag);
  reformat_P(gl,&P,"",&ref_flag);
  reformat_T(gl,&T,"",&ref_flag);
  /* not necessary I think
     ???????????????????????
     reformat_w(w,&ref_flag);
  */
  if (ref_flag) {
    find_U_2(np, l, gl, w, V, P, T, k, psi);
  }
  find_prim_fluid_mem(np, l, gl, P, T);
}



void add_dUstar_to_U(np_t *np, long l, gl_t *gl, flux_t dUstar){
  double Omega;
  long flux;
  Omega=_Omega(np[l],gl);
  for (flux=0; flux<nf; flux++){
    np[l].bs->U[flux]=np[l].bs->U[flux]+dUstar[flux]/Omega;
  }
  find_prim_fluid(np,l,gl);
}


/*  -=USED FOR TESTING PURPOSES ONLY=-
void check_domain_values(np_t *np, gl_t *gl, char *message){
  long i,j,k;
  double P;
  resume_nodes_only_in_zone_and_update_bdry_nodes(np,gl,gl->domain);
  for1DL(i,gl->domain.is,gl->domain.ie)
    for2DL(j,gl->domain.js,gl->domain.je)
      for3DL(k,gl->domain.ks,gl->domain.ke)
        P=_Pstar(np[_ai(gl,i,j,k)]);
        if (P<4000 || P>6000) {
          printf("problem here %s\n",message);
          exit(EXIT_FAILURE);
        }
      end3DL
    end2DL
  end1DL

}*/





