// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2010-2011 Bernard Parent
Copyright 2021 Prasanna Thoguluva Ranjendran


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

#include "fluid.h"
#include "fluid_conv.h"
#include "fluid_source.h"
#include "fluid_bdry.h"
#include <src/control.h>
#include <src/common.h>
#include <model/metrics/_metrics.h>
#include <model/_model.h>
#include <cycle/_cycle.h>
#include <model/share/fluid_share.h>


#if ns > 1
  #error The Navier-Stokes_perfect module can not be compiled with a chemical model that has more than 1 species. 
#endif

#define ETAMODEL_CONSTANT 1
#define ETAMODEL_SUTHERLAND 2

static bool _FLUIDPRIMMEM(np_t np){
  return(np.FLUIDPRIMMEM);
}


void write_model_fluid_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    gamma=1.4;\n"
    "    R=286.0;        {J/kgK }\n"
    "    ETAMODEL=ETAMODEL_CONSTANT;     \n"
    "    eta=1.716e-5;       {kg/ms}\n"
    "    Pr=0.71;     \n"
    "    Pmin=1.0e-2;    Pmax=9.9e99;   {Pa}\n"
    "    Tmin=1.0e1;     Tmax=26.0e3;    {K}\n"
#ifdef _2D
    "    AXISYMMETRIC=FALSE;\n"
#endif
    "    SetBodyForce(is,"if2DL("js,")if3DL("ks,")" ie,"if2DL("je,")if3DL("ke,")" 0.0{N/m3}"if2DL(",0.0{N/m3}")if3DL(",0.0{N/m3}")");\n"
    "    SetHeatDeposited(is,"if2DL("js,")if3DL("ks,")" ie,"if2DL("je,")if3DL("ke,")" 0.0 {W/m3});\n"
    "    {\n"
    "    AddHeatPoint(0.0{x,m},"if2DL("0.0{y,m},")if3DL("0.0{z,m},")" 0.1{radius,m}, 0.0{W"if2D("/m")"});\n"
    "    }\n"
    "  );\n"
  ,_FLUID_ACTIONNAME);


}


void write_cycle_fluid_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    xiverge=1e-3;  {residual convergence threshold}\n"
    "    rhoref=1.0;  {reference density in kg/m3}\n"
    "    aref=300.0;  {reference sound speed in m/s}\n" 
    "    Uref[1]=rhoref;            \n"
    "    for (dim,1,numdim,\n"
    "      Uref[1+dim]=rhoref*aref;\n"
    "    );\n"
    "    Uref[2+numdim]=rhoref*aref*aref;\n"  
    "  );\n"
  ,_FLUID_ACTIONNAME);  
}


void read_model_fluid_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
  read_model_fluid_actions_Fbody_Qadd(actionname, argum, codex);
}


void read_model_fluid_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  long numvarsinit,i,j,k,dim;

  if (strcmp(actionname,_FLUID_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);

    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_FLUID_ACTIONNAME);

    SOAP_add_int_to_vars(codex,"ETAMODEL_CONSTANT",ETAMODEL_CONSTANT);
    SOAP_add_int_to_vars(codex,"ETAMODEL_SUTHERLAND",ETAMODEL_SUTHERLAND);

    if (!gl->CONTROL_READ){
      gl->MODEL_FLUID_READ=TRUE;
      for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          (*((readcontrolarg_t *)codex->action_args)->np)[_ai(((readcontrolarg_t *)codex->action_args)->gl,i,j,k)].bs->Qadd=0.0e0;
          for (dim=0; dim<nd; dim++) {
              (*((readcontrolarg_t *)codex->action_args)->np)[_ai(((readcontrolarg_t *)codex->action_args)->gl,i,j,k)].bs->Fbody[dim]=0.0e0;
          }
      }
    }
    action_original=codex->action;
    codex->action=&read_model_fluid_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_double_var_from_codex(codex,"gamma",&gl->model.fluid.gamma);
    find_double_var_from_codex(codex,"R",&gl->model.fluid.R);
    find_int_var_from_codex(codex,"ETAMODEL",&gl->model.fluid.ETAMODEL);
    find_double_var_from_codex(codex,"eta",&gl->model.fluid.eta);
    if (gl->model.fluid.ETAMODEL==ETAMODEL_SUTHERLAND){
      find_double_var_from_codex(codex,"eta_T0",&gl->model.fluid.eta_T0);
      find_double_var_from_codex(codex,"eta_C",&gl->model.fluid.eta_C);
    } else {
      gl->model.fluid.eta_T0=0.0;
      gl->model.fluid.eta_C=0.0;
    }
    find_double_var_from_codex(codex,"Pr",&gl->model.fluid.Pr);
    find_double_var_from_codex(codex,"Pmin",&gl->model.fluid.Pmin);
    find_double_var_from_codex(codex,"Pmax",&gl->model.fluid.Pmax);
    find_double_var_from_codex(codex,"Tmin",&gl->model.fluid.Tmin);
    find_double_var_from_codex(codex,"Tmax",&gl->model.fluid.Tmax);

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


double _rho (np_t np) {
  double tmp;
  tmp=np.bs->U[0];
  return(tmp);
}


double _V(np_t np, long theta) {
  double tmp;
  if (_FLUIDPRIMMEM(np)){
    tmp=np.wk->Vmem[theta];
  } else {
    tmp=np.bs->U[1+theta]/_rho(np);
  }
  return(tmp);
}


double _V_from_U(np_t np, long theta) {
  double tmp;
  tmp=np.bs->U[1+theta]/np.bs->U[0];
  return(tmp);
}


void find_V(np_t np, dim_t V) {
  long dim;
  double rho;
  if (_FLUIDPRIMMEM(np)){
    for (dim=0; dim<nd; dim++) V[dim]=np.wk->Vmem[dim];
  } else {
    rho=_rho(np);
    for (dim=0; dim<nd; dim++) V[dim]=np.bs->U[1+dim]/rho;
  }
}


double _q(np_t np) {
  double tmp;
  long dim;
  tmp=0.0;
  for (dim=0; dim<nd; dim++) tmp+=sqr(_V(np,dim));
  tmp=sqrt(tmp);
  return(tmp);
}


double _a(np_t np, gl_t *gl) {
  double tmp,T;
  if (_FLUIDPRIMMEM(np)){
    tmp=np.wk->amem;
  } else {
    T=_T(np,gl);
    assert_np(np,T>=0.0e0);
    tmp=sqrt(max(0.0e0,gl->model.fluid.gamma*gl->model.fluid.R*T));
  }
  return(tmp);
}



double _P (np_t np, gl_t *gl) {
  double tmp,T,rho;
  if (_FLUIDPRIMMEM(np)){
    tmp=np.wk->Pmem;
  } else {
    rho=_rho(np);
    T=_T(np,gl);
    tmp=rho*gl->model.fluid.R*T;
  }
  return(tmp);
}


double _T (np_t np, gl_t *gl) {
  double ekin,eint,tmp,rho;
  long dim;
  if (_FLUIDPRIMMEM(np)){
    tmp=np.wk->Tmem;
  } else {
    rho=_rho(np);
    assert(rho!=0.0e0);
    ekin=0.0e0;
    for (dim=0; dim<nd; dim++){
      ekin=ekin+0.5e0*sqr(np.bs->U[1+dim]/rho);
    }
    eint=np.bs->U[1+nd]/rho /* total energy */
      -ekin;  /* kinetic energy */
    tmp=eint/gl->model.fluid.R*(gl->model.fluid.gamma-1.0);
  }
  return(tmp);
}





double _et (np_t np){
  double tmp;
  assert_np(np,_rho(np)!=0.0e0);
  tmp=np.bs->U[1+nd]/_rho(np);
  return(tmp);
}



double _ht (np_t np, gl_t *gl){
  double tmp;
  assert_np(np,_rho(np)!=0.0e0);
  tmp=np.bs->U[1+nd]/_rho(np)+_P(np,gl)/_rho(np);
  return(tmp);
}



double _dVi_dxj(np_t *np, long l, gl_t *gl, long i, long j){
  double sum;
  long dim;
  sum=0.0e0;
  if (is_node_valid(np[l], TYPELEVEL_FLUID_WORK)){
    for (dim=0; dim<nd; dim++){
     /* at the boundary nodes, only the derivatives of the velocity not perpendicular
        to the boundary plane are considered */
      if (is_node_valid(np[_al(gl,l,dim,+1)], TYPELEVEL_FLUID_WORK) 
       && is_node_valid(np[_al(gl,l,dim,-1)], TYPELEVEL_FLUID_WORK)) {
         sum+=_X(np[l], dim,j)*0.5*
             (_V_from_U(np[_al(gl,l,dim,+1)],i)-_V_from_U(np[_al(gl,l,dim,-1)],i));
      }
    }
  }
  return(sum);
}


void find_Ustar(np_t np, gl_t *gl, flux_t Ustar){
  double Omega;
  long flux;

  Omega=_Omega(np,gl);
  for (flux=0; flux<nf; flux++)
    Ustar[flux]=np.bs->U[flux]*Omega;
}


void find_G(np_t np, gl_t *gl, flux_t G){
  long dim;
  G[0]=0.0;
  for (dim=0; dim<nd; dim++) G[1+dim]=_V(np,dim);
  G[1+nd]=_T(np,gl);
}



void find_dT_dx(np_t np, gl_t *gl, double *dTdrhoet, double *dTdrho, dim_t dTdrhoV){
  double et,dedT,sum,rho;
  dim_t V;
  long dim;
  double dedrho;

  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=_V(np,dim);
    sum=sum+sqr(V[dim]);
  }
  rho=_rho(np);
  et=_et(np);
  dedT=gl->model.fluid.R/(gl->model.fluid.gamma-1.0); /* Cv*/
  dedrho=0.0;

  /* ---- dTdrhoet ------*/
  assert_np(np,dedT*rho!=0.0e0);
  *dTdrhoet=1.0e0/(rho*dedT);

  /* ---- dTdrho ------*/
  assert_np(np,rho!=0.0e0);
  assert_np(np,dedT!=0.0e0);

  *dTdrho=-(et-sum+dedrho*rho)/dedT/rho;

  /* ---- dTdrhou ------*/
  for (dim=0; dim<nd; dim++){
    dTdrhoV[dim]=-V[dim]*(*dTdrhoet);
  }


}


void find_dG_dUstar(np_t np, gl_t *gl, sqmat_t B){
  long dim,row,col;
  dim_t dTdrhoV;
  double dTdrho,Omega,rho,dTdrhoet;

  Omega=_Omega(np,gl);
  rho=_rho(np);

  find_dT_dx(np, gl, &dTdrhoet, &dTdrho, dTdrhoV);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=0.0e0;
    }
  }

    assert_np(np,rho!=0.0e0);
    assert_np(np,Omega!=0.0e0);
    for (dim=0; dim<nd; dim++){
      B[1+dim][0]=-_V(np,dim)/(Omega*rho);
    }
    B[1+nd][0]=dTdrho/(Omega);


  for (dim=0; dim<nd; dim++){
    B[1+dim][1+dim]=1.0e0/(Omega*rho);
    B[1+nd][1+dim]=dTdrhoV[dim]/(Omega);
  }
  B[1+nd][1+nd]=dTdrhoet/(Omega);
}







void find_Kstar_interface(np_t *np, gl_t *gl, long lL, long lR, metrics_t metrics, long theta, long vartheta, sqmat_t K, int CYCLELEVEL){
  long dim,dim2,row,col;
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
    
    for (dim=0; dim<nd; dim++){
      for (dim2=0; dim2<nd; dim2++){
        K[1+dim][1+dim2]=avg(_eta(np,lL,gl),_eta(np,lR,gl))*beta[dim][dim2];
        K[1+nd][1+dim2]=K[1+nd][1+dim2]+
           avg(_eta(np,lL,gl),_eta(np,lR,gl))*beta[dim][dim2]*avg(_V(np[lL],dim),_V(np[lR],dim));
      }
    }
    
    K[1+nd][1+nd]=avg(_kappa(np,lL,gl),_kappa(np,lR,gl))*alpha;
    

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


void reformat_rho(gl_t *gl, double *rho, char *suffix, bool *flag){
  double rhomin,rhomax;
  rhomin=gl->model.fluid.Pmin/(gl->model.fluid.R*gl->model.fluid.Tmax);
  rhomax=gl->model.fluid.Pmax/(gl->model.fluid.R*gl->model.fluid.Tmin);
  if (*rho<rhomin) {
    add_to_clipped_variables2(gl,"rhomin",suffix);
    *rho=rhomin;
    *flag=TRUE;
  }
  if (*rho>rhomax) {
    add_to_clipped_variables2(gl,"rhomax",suffix);
    *rho=rhomax;
    *flag=TRUE;
  }
}




void find_prim_fluid_mem(np_t *np, long l, gl_t *gl, double P, double T){
  double rho;
  long dim;

  assert_np(np[l],is_node_resumed(np[l]));
  /*
  reformat_T(&T,&ref_flag);
  reformat_P(&P,&ref_flag); */
  np[l].wk->Tmem=T;
  np[l].wk->Pmem=P;
  np[l].wk->amem=sqrt(gl->model.fluid.gamma*gl->model.fluid.R*T);
  rho=P/(T*gl->model.fluid.R);
  np[l].wk->rhomem=rho;
  assert_np(np[l],np[l].wk->rhomem>0.0e0);
 /* find Vmem */
  assert_np(np[l],rho!=0.0e0);
  for (dim=0; dim<nd; dim++){
    np[l].wk->Vmem[dim]=np[l].bs->U[1+dim]/rho;
  }

  np[l].FLUIDPRIMMEM=TRUE;
}


void find_U(np_t *np, long l, gl_t *gl, double rho, dim_t V, double T){
  double sum;
  long dim;
  bool ref_flag;
  
  ref_flag=FALSE;
  reformat_T(gl,&T,"",&ref_flag);
  
  np[l].bs->U[0]=rho;
  assert_np(np[l],rho!=0.0e0);

  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    np[l].bs->U[1+dim]=rho*V[dim];
    sum=sum+sqr(V[dim]);
  }
  np[l].bs->U[1+nd]=rho*(0.5e0*(sum)+gl->model.fluid.R/(gl->model.fluid.gamma-1.0)*T);
  if (is_node_resumed(np[l])) {
    np[l].wk->Tmem=T;
    find_prim_fluid(np,l,gl);
  }
}


void find_U_2(np_t *np, long l, gl_t *gl, dim_t V, double P, double T){
  double sum,rho;
  long dim;
  bool ref_flag;
 
  ref_flag=FALSE;
  reformat_T(gl,&T,"",&ref_flag);
  reformat_P(gl,&P,"",&ref_flag);

  rho=P/gl->model.fluid.R/T;
  np[l].bs->U[0]=rho;
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    np[l].bs->U[1+dim]=rho*V[dim];
    sum=sum+sqr(V[dim]);
  }
  np[l].bs->U[1+nd]=rho*(0.5e0*(sum) + gl->model.fluid.R/(gl->model.fluid.gamma-1.0)*T);

  if (is_node_resumed(np[l])){
    np[l].wk->Tmem=T;
    find_prim_fluid_mem(np, l, gl, P, T);
  }
}


void find_U_3(np_t *np, long l, gl_t *gl, double rho, dim_t V, double P){
  double T;
  bool ref_flag;

  reformat_P(gl,&P,"",&ref_flag);
  assert_np(np[l],P>0.0e0);
  T=P/rho/gl->model.fluid.R;
  find_U(np, l, gl, rho, V, T);
}



void find_prim_fluid(np_t *np, long l, gl_t *gl){
  long dim;
  double rho,P,T,sum,e,Cv;
  dim_t V;
  bool ref_flag;

  rho=np[l].bs->U[0];

  assert_np(np[l],rho!=0.0e0);

  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=np[l].bs->U[1+dim]/rho;
    sum=sum+sqr(V[dim]);
  }
  assert_np(np[l],rho!=0.0e0);

  e=np[l].bs->U[1+nd]/rho-0.5e0*sum;
  /* clip e and adjust U*/
  Cv=gl->model.fluid.R/(gl->model.fluid.gamma-1.0);
  T=e/Cv;

  if (T<gl->model.fluid.Tmin){
    T=gl->model.fluid.Tmin; 
    add_to_clipped_variables(gl,"Tmin");
  }
  if (T>gl->model.fluid.Tmax){
    T=gl->model.fluid.Tmax; 
    add_to_clipped_variables(gl,"Tmax");
  }
  np[l].bs->U[1+nd]=rho*(Cv*T+0.5e0*sum);


  np[l].wk->Tmem=T;

  P=rho*gl->model.fluid.R*T;

  ref_flag=FALSE;

  reformat_P(gl,&P,"",&ref_flag);

  reformat_T(gl,&T,"",&ref_flag);

  if (ref_flag) {
   find_U_2(np, l, gl, V, P, T);
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

double _kappa(np_t *np, long l, gl_t *gl) {
  double Pr=gl->model.fluid.Pr;
  double Cp=gl->model.fluid.gamma*gl->model.fluid.R/(gl->model.fluid.gamma-1.0);

  return(_eta_from_T(_T(np[l],gl),gl)*Cp/Pr);
}

double _eta_from_T(double T,gl_t *gl) {
  double mu_0=gl->model.fluid.eta;
  double C=gl->model.fluid.eta_C;
  double T_0=gl->model.fluid.eta_T0;

  switch (gl->model.fluid.ETAMODEL){
    case ETAMODEL_CONSTANT:
      return(mu_0);
      break;
    case ETAMODEL_SUTHERLAND:
      return(mu_0*(T_0+C)/(T+C)*sqrt(T*T*T/(T_0*T_0*T_0)));
      break;
    default:
      fatal_error("Given ETAMODEL is not a valid choice for viscosity model. Input either ETAMODEL_CONSTANT or ETAMODEL_SUTHERLAND.");
      return(-1); //to avoid compiler warning
  }
}

double _eta(np_t *np, long l, gl_t *gl) {
  return(_eta_from_T(_T(np[l],gl),gl));
}




