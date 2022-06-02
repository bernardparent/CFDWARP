// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015-2018,2020,2021 Bernard Parent
Copyright 2001 Jason Etele
Copyright 2002 Thomas E. Schwartzentruber

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


#include <model/share/fluid_share.h>
#include <model/share/emfield_share.h>
#include <model/fluid/_fluid.h>
#include <model/thermo/_thermo.h>
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

#define PECLETMAX 20.0
#define PRESERVE_MOMENTUM TRUE
#define MINRATIO_BDRY_RHO_P 0.1





void write_disc_fluid_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
#ifdef _FLUID_NEUTRALSTRANSPORT
    "    zetaA1=0.1e0;    {conditions the eigenvalue of the A jacobian}\n"
    "    zetaA2=0.2e0;     {conditions eigenvalues so they don't go out of bound; set to 0.1-0.5}\n"
#ifdef _RESCONV_ZETAADEN
    "    zetaA3=0.2e0; {conditions the eigenvalues of the inverted A jacobian}\n"
#endif
#endif //_FLUID_NEUTRALSTRANSPORT
#ifdef _FLUID_DRIFTDIFFUSION
    "    aref=300.0;       {reference speed of sound in m/s used when conditioning Dstar eigenvalues  }\n"
#endif
#ifdef _FLUID_PLASMA
    "    zetaD=0.0;        {conditions the Dstar eigenvalues for the charged species}\n"
#endif
  ,_FLUID_ACTIONNAME);
#ifdef _FLUID_PLASMA
  wfprintf(*controlfile,
    "    for (spec,1,numspec,\n"
    "      if (SPECIESTYPE[spec]==SPECIESTYPE_IONPLUS,\n"
    "        betag[spec]=1.0;\n"
    "        betaa[spec]=0.0;\n"
    "      );\n"
    "      if (SPECIESTYPE[spec]==SPECIESTYPE_IONMINUS,\n"
    "        betag[spec]=-0.01;\n"
    "        betaa[spec]=0.5;\n"
    "      );\n"
    "      if (SPECIESTYPE[spec]==SPECIESTYPE_ELECTRON,\n"
    "        betag[spec]=-0.001;\n"
    "        betaa[spec]=0.999;\n"
    "      );\n"
    "    );\n"
  );
#endif
  wfprintf(*controlfile,
    "  );\n");
}


void read_disc_fluid_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_disc_fluid_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long numvarsinit;
#ifdef _FLUID_PLASMA
  long spec;
  char tmpstr[400];
#endif
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_FLUID_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_FLUID_ACTIONNAME);


    gl->DISC_FLUID_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_fluid_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

#ifdef _FLUID_NEUTRALSTRANSPORT
    find_double_var_from_codex(codex,"zetaA1",&gl->model.fluid.zetaA1);
    find_double_var_from_codex(codex,"zetaA2",&gl->model.fluid.zetaA2);
#ifdef _RESCONV_ZETAADEN
    find_double_var_from_codex(codex,"zetaA3",&gl->model.fluid.zetaA3);
#else
    gl->model.fluid.zetaA3=1e-30;
#endif
#endif //_FLUID_NEUTRALSTRANSPORT

#ifdef _FLUID_DRIFTDIFFUSION
    find_double_var_from_codex(codex,"aref",&gl->model.fluid.aref);
#endif
#ifdef _FLUID_PLASMA
    find_double_var_from_codex(codex,"zetaD",&gl->model.fluid.zetaD);
    for (spec=0; spec<ns; spec++){
      if (spec<ncs){
        sprintf(tmpstr,"betag[%ld]",spec+1);
        find_double_var_from_codex(codex,tmpstr,&(gl->model.fluid.betag[spec]));
        sprintf(tmpstr,"betaa[%ld]",spec+1);
        find_double_var_from_codex(codex,tmpstr,&(gl->model.fluid.betaa[spec]));
      } else {
        gl->model.fluid.betag[spec]=0.0;
        gl->model.fluid.betaa[spec]=0.0;
      }    
    }
#endif

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}




void find_chi_inverse(dim2_t X, dim2_t ChiInv){
  EXM_mat_t chi,chiinv;
  long dim1,dim2;

  EXM_init_matrix(&chi,nd,nd);
  EXM_init_matrix(&chiinv,nd,nd);
  for (dim1=0; dim1<nd; dim1++){
    for (dim2=0; dim2<nd; dim2++){
      chi.cont[EXM_aim(chi.glm,dim1,dim2)]=X[dim1][dim2];
    }
  }
  EXM_invert_matrix_analytical(chi,&chiinv);
  for (dim1=0; dim1<nd; dim1++){
    for (dim2=0; dim2<nd; dim2++){
      ChiInv[dim1][dim2]=chiinv.cont[EXM_aim(chiinv.glm,dim1,dim2)];
    }
  }
  EXM_free_matrix(&chi);
  EXM_free_matrix(&chiinv);
}


#if defined(_FLUID_NEUTRALSTRANSPORT)

#if defined(_FLUID_MULTISPECIES)
/* average ratio of specific heats for the heavy species */
double _gamma(np_t np, gl_t *gl){
  spec_t w;
  double cp,R,gamma;

  find_w(np,w);
  cp=_cp_from_w_T(w, _T(np,gl));
  R=_R(w);
  gamma=cp/(cp-R);
  return(gamma);
}


double _Rgas(np_t np, gl_t *gl){
  double R;
  spec_t w;
  R=_R(w);
  return(R);
}

#else
double _gamma(np_t np, gl_t *gl){
  return(gl->model.fluid.gamma);
}


double _Rgas(np_t np, gl_t *gl){
  return(gl->model.fluid.R);
}

#endif



double _Vstar (np_t np, long theta) {
  double ret;
  long dim;
  ret=0.0e0;
  for (dim=0; dim<nd; dim++){
    ret=ret+_X(np,theta,dim)*_V(np,dim);
  }
  return(ret);
}


double _Vhat (np_t np, long theta) {
  double ret,numer,denom;
  long dim;
  numer=0.0e0;
  denom=0.0e0;
  for (dim=0; dim<nd; dim++){
    numer=numer+_X(np,theta,dim)*_V(np,dim);
    denom=denom+sqr(_X(np,theta,dim));
  }
  assert_np(np,denom>0.0e0);
  assert_np(np,denom!=0.0e0);
  ret=numer/sqrt(denom);
  return(ret);
}

double _varphi(np_t *np, gl_t *gl, long l, long dim){
  np_t np_p1,np_p0;
  double varphi,a,Vi,Xhati,Mi,gam,rho;
  long dim2;

  np_p0=np[l];
  np_p1=np[_al(gl,l,dim,+1)];
  assert_np(np_p0,is_node_resumed(np_p0));
  assert_np(np_p1,is_node_resumed(np_p1));
  Xhati=0.0e0;
  for (dim2=0; dim2<nd; dim2++) Xhati+=sqr(_X(np_p0,dim,dim2));
  Xhati=sqrt(Xhati);
  Vi=_Vstar(np_p0,dim);
  a=_a(np_p0,gl);
  Mi=Vi/a/Xhati;
  gam=_gamma(np_p0,gl);
  rho=_rho(np_p0);
  varphi=fabs( Xhati/rho/a*max(0.0e0,(1.0e0-sqr(Mi))/(1.0e0+(gam-1.0e0)*sqr(Mi)))
          *(_Pstar(np_p1,gl)-_Pstar(np_p0,gl))
         );
  return(varphi);
}


void find_P_bdry_symmetrical(np_t *np, gl_t *gl,  long lA, long lB, long lC, long lD, long theta, long thetasgn, bool BDRYDIRECFOUND, int ACCURACY, double *P){
  double sum;
  long dim;
  sum=0.0e0;
  if (PRESERVE_MOMENTUM && BDRYDIRECFOUND){
    for (dim=0; dim<nd; dim++){
      sum=sum+sqr(_X(np[lB],theta,dim));
    }
    assert(sum!=0.0);
    sum=-(double)thetasgn*_Vstar(np[lB],theta)*fabs(_Vstar(np[lB],theta))*_rho(np[lB])/sum;
  }
  *P=_f_symmetry(ACCURACY,_P(np[lB],gl),_P(np[lC],gl),_P(np[lD],gl))+sum;
  *P=max(_P(np[lB],gl)*MINRATIO_BDRY_RHO_P,max(gl->model.fluid.Pmin*1.00001e0,*P));
}


void find_Pstar_bdry_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,  bool BDRYDIRECFOUND, int ACCURACY,  double *Pstar){
  double sum;
  long dim;
  if (PRESERVE_MOMENTUM) {
    sum=0.0e0;
    if (BDRYDIRECFOUND){
      for (dim=0; dim<nd; dim++){
        sum=sum+sqr(_X(np[lB],theta,dim));
      }
      assert(sum!=0.0);
      sum=-(double)thetasgn*_Vstar(np[lB],theta)*fabs(_Vstar(np[lB],theta))*_rho(np[lB])/sum;
    }
    *Pstar=_f_symmetry(ACCURACY,_Pstar(np[lB],gl),_Pstar(np[lC],gl))
             +sum;
  } else {
    *Pstar=_f_symmetry(ACCURACY,_Pstar(np[lB],gl),_Pstar(np[lC],gl));
  }
  *Pstar=max(_Pstar(np[lB],gl)*MINRATIO_BDRY_RHO_P,max(gl->model.fluid.Pmin*1.00001e0,*Pstar));
}


void find_V_P_T_bdry_freestream(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, int ACCURACY, dim_t V, double *P, double *T, bool *OUTFLOW, bool *FREESTREAM){
  dim_t Vnew,Vstarnew,Vinf,Vextrapol;
  long dim,i,j;
  double MB,Rgas,Tnew,Pnew,gamma,Mnew2,Tinf,qnew2,Pinf,qinf2,Pstaginf;
  dim2_t chiinv;

  MB=_q(np[lB])/_a(np[lB],gl);
  gamma=_gamma(np[lB],gl);
  Rgas=_Rgas(np[lB],gl);
  for (dim=0; dim<nd; dim++) Vinf[dim]=_bdry_param(np,gl,lA,dim,TYPELEVEL_FLUID_WORK);
  Pinf=_bdry_param(np,gl,lA,nd,TYPELEVEL_FLUID_WORK);
  Tinf=_bdry_param(np,gl,lA,nd+1,TYPELEVEL_FLUID_WORK);
  qinf2=0.0;
  for (dim=0; dim<nd; dim++) qinf2+=sqr(Vinf[dim]);
  Pstaginf=Pinf*pow(1.0+(gamma-1.0)/2.0*qinf2/(gamma*Rgas*Tinf),gamma/(gamma-1.0));

  for (dim=0; dim<nd; dim++){
    Vstarnew[dim]=_f_symmetry(ACCURACY,_Vstar(np[lB],dim),_Vstar(np[lC],dim));
  }
  find_chi_inverse(np[lA].bs->X,chiinv);
  *FREESTREAM=TRUE;
  for (i=0; i<nd; i++){
    Vextrapol[i]=0.0e0;
    for (j=0; j<nd; j++){
      Vextrapol[i]+=chiinv[i][j]*Vstarnew[j];
    }
    if (sqr(Vextrapol[i]-Vinf[i])>0.3*qinf2) *FREESTREAM=FALSE;
  }

  if (*FREESTREAM){
    if (Vstarnew[theta]*thetasgn>0.0){
      /* inflow bdry */
      *OUTFLOW=FALSE;
      if (MB<1.0){
        /* subsonic inflow bdry */
        find_chi_inverse(np[lA].bs->X,chiinv);
        qnew2=0.0;
        for (dim=0; dim<nd; dim++){
          Vnew[dim]=Vextrapol[dim];
          qnew2+=sqr(Vnew[dim]);
        }
        // Hnew=Hinf
        Tnew=(Rgas*Tinf+0.5*(qinf2-qnew2)*(1.0-1.0/gamma))/Rgas;
        Mnew2=qnew2/(gamma*Rgas*Tnew);
        //Pstagnew=Pstagfreestream
        Pnew=Pstaginf/pow(1.0+(gamma-1.0)/2.0*Mnew2,gamma/(gamma-1.0));
      } else {
        /* supersonic inflow bdry */
        for (dim=0; dim<nd; dim++) Vnew[dim]=Vinf[dim];
        Pnew=Pinf;
        Tnew=Tinf;
      }

    } else {   
      /* outflow bdry */
      *OUTFLOW=TRUE;
      if (MB<1.0){
        /* subsonic outflow bdry */
        for (dim=0; dim<nd; dim++) Vnew[dim]=Vextrapol[dim];
        Tnew=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
        Pnew=Pinf;
      } else {
        /* supersonic outflow bdry */
        for (dim=0; dim<nd; dim++) Vnew[dim]=_V(np[lB],dim);
        Pnew=_P(np[lB],gl);
        Tnew=_T(np[lB],gl);
      }
    }
  } else {
    for (dim=0; dim<nd; dim++) Vnew[dim]=_V(np[lB],dim);
    if (MB>1.0){
      Pnew=_P(np[lB],gl);
    } else {
      Pnew=Pinf;
    } 
    Tnew=_T(np[lB],gl);    
    *OUTFLOW=TRUE;
  }

  *P=Pnew;
  *T=Tnew;
  for (dim=0; dim<nd; dim++) V[dim]=Vnew[dim];

}

#endif



#if defined(_FLUID_MULTISPECIES) && defined(_FLUID_NEUTRALSTRANSPORT)


void find_w_V_P_T_bdry_inflow_reservoir(np_t *np, gl_t *gl, long lA, long lB, long lC, int ACCURACY, spec_t w, dim_t V, double *P, double *T){
  long dim;
  double sum,gam,cp,R,Pstag,Tstag,M2,sum_max;

  assert_np(np[lA],is_node_resumed(np[lA]));
  find_w(np[lA],w);
  *T=_T(np[lA],gl);
  /* first, find gamma at the boundary */
  cp=_cp_from_w_T(w, *T);
  R=_R(w);
  gam=(cp)/(cp-R);

  /* find stagnation pressure and temperature at node C */
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    sum=sum+sqr(_V(np[lA],dim));
  }
  M2=sum/(gam*R*_T(np[lA],gl));
  Pstag=_P(np[lA],gl)*pow(1.0e0+(gam-1.0e0)*M2/2.0e0,gam/(gam-1.0e0));
  Tstag=_T(np[lA],gl)*(1.0e0+(gam-1.0e0)*M2/2.0e0);
  /* find the velocity at node C, extrapolated from node B*/

  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=_V(np[lB],dim);
    sum=sum+sqr(V[dim]);
  }
  sum_max=(Tstag-gl->model.fluid.Tmin)*2.0e0*gam*R/(gam-1.0e0);
  if (sum>sum_max) {
    for (dim=0; dim<nd; dim++){
      V[dim]=V[dim]*sqrt(sum_max/sum);
    }
    sum=sum_max;
  }
  *T=Tstag-(gam-1.0e0)*sum/(gam*R)/2.0e0;
  assert_np(np[lA],*T>0.0e0);
  M2=sum/(gam*R*(*T));
  *P=Pstag/pow(1.0e0+(gam-1.0e0)*M2/2.0e0,gam/(gam-1.0e0));
}


static double find_fofT(double Pstag, double Vstar, double mdot, double gamma, double R, double T){
  double B,D,F;
  double ret;
  B=Pstag*Pstag*Vstar*Vstar/(R*R*mdot*mdot);
  D=(gamma-1.0)*Vstar*Vstar/(2.0*gamma*R);
  F=-(gamma+1.0)/(gamma-1.0);
  ret=T-B/T*pow(1.0+D/T,F)+D;
  return(ret);
}


/* function by Jason Etele, to be used with Reservoir2 bdry condition */
static double find_fprimeofT(double Pstag, double Vstar, double mdot, double gamma, double R, double T){
  double B,D,F,K;
  double ret;
  B=Pstag*Pstag*Vstar*Vstar/(R*R*mdot*mdot);
  D=(gamma-1.0)*Vstar*Vstar/(2.0*gamma*R);
  F=-(gamma+1.0)/(gamma-1.0);
  K=1.0+D/T;
  ret=1.0+B*D*F/(T*T*T)*pow(K,F-1.0) + B/(T*T)*pow(K,F);
  return(ret);
}


/* function by Jason Etele, to be used with Reservoir2 bdry condition */
static double find_mdot_over_A(double Pstag, double Tstag, double M, double gamma, double R){
  double fcnM,ret;
  fcnM=M*pow(1.0+(gamma-1.0)/2.0*M*M,-(gamma+1.0)/(2.0*(gamma-1.0)));
  ret=Pstag*pow(gamma/(R*Tstag),0.5)*fcnM;
  return(ret);
}


void find_w_V_P_T_bdry_inflow_reservoir_2(np_t *np, gl_t *gl, long lA, long lB, long lC, int ACCURACY, spec_t w, dim_t V, double *P, double *T){
  int i;
  long dim;
  double sum,gam,cp,R,Pstagset,Tstag,mdotset,M2,sum_max;
  double a,Vtot,M,fofT,fprimeofT,Tnew;

  assert_np(np[lA],is_node_resumed(np[lA]));
  find_w(np[lA],w);
  *T=_T(np[lA],gl);
  /* first, find gamma at the boundary */
  cp=_cp_from_w_T(w, *T);
  R=_R(w);
  gam=(cp)/(cp-R);

  /* find stagnation pressure and massflow divided by area at node C */
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    sum=sum+sqr(_V(np[lA],dim));
  }
  M2=sum/(gam*R*(*T));
  M=sqrt(M2);
  Pstagset=_P(np[lA],gl)*pow(1.0e0+(gam-1.0e0)*M2/2.0e0,gam/(gam-1.0e0));
  Tstag=_T(np[lA],gl)*(1.0e0+(gam-1.0e0)*M2/2.0e0);
  mdotset=find_mdot_over_A(Pstagset,Tstag,M,gam,R);

  /* find the velocity at node C, extrapolated from node B*/
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=_V(np[lB],dim);
    sum=sum+sqr(V[dim]);
  }
  sum_max=(Tstag-gl->model.fluid.Tmin)*2.0e0*gam*R/(gam-1.0e0);
  if (sum>sum_max) {
    for (dim=0; dim<nd; dim++){
      V[dim]=V[dim]*sqrt(sum_max/sum);
    }
    sum=sum_max;
  }
  Vtot=sqrt(sum);

  /* use Newton's method to find temperature from ptot, mdot/A
     (which are held constant), and Vtot(which is read from the interior); */
  assert_np(np[lA],*T>0.0e0);
  for(i=0;i<7;i++){
    a=sqrt(gam*R*(*T));
    M=Vtot/a;
    fofT=find_fofT(Pstagset,Vtot,mdotset,gam,R,*T);
    fprimeofT=find_fprimeofT(Pstagset,Vtot,mdotset,gam,R,*T);
    Tnew=*T-fofT/fprimeofT;
    *T=Tnew;
  }
  M2=sum/(gam*R*(*T));  /*Note: we neglect the dependence of gamma on temp*/
  *P=Pstagset/pow(1.0e0+(gam-1.0e0)*M2/2.0e0,gam/(gam-1.0e0));
}


void find_w_V_P_T_bdry_symmetrical(np_t *np, gl_t *gl,  long lA, long lB, long lC, long lD, long theta, long thetasgn, 
                                   bool BDRYDIRECFOUND, int ACCURACY, spec_t w, dim_t V, double *P, double *T){
  dim_t Vstar;
  long dim,i,j,spec;
  double rho;
  dim2_t chiinv;

  for (dim=0; dim<nd; dim++){
    Vstar[dim]=_f_symmetry(ACCURACY,_Vstar(np[lB],dim),_Vstar(np[lC],dim),_Vstar(np[lD],dim));
  }
  Vstar[theta]=0.0e0;
  find_chi_inverse(np[lA].bs->X,chiinv);
  for (i=0; i<nd; i++){
    V[i]=0.0e0;
    for (j=0; j<nd; j++){
      V[i]=V[i]+chiinv[i][j]*Vstar[j];
    }
  }
  for (spec=0; spec<ns; spec++){
      w[spec]=_f_symmetry(ACCURACY,_w(np[lB],spec),_w(np[lC],spec),_w(np[lD],spec));
  }

  rho=_f_symmetry(ACCURACY,_rho(np[lB]),_rho(np[lC]),_rho(np[lD]));
  *T=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl),_T(np[lD],gl));

  /* here, correct the charged species so that dN/deta is zero at the symmetry surface */
  for (spec=0; spec<ncs; spec++){
    w[spec]=_w(np[lB],spec)*_rho(np[lB])/rho;
  }

  find_P_bdry_symmetrical(np, gl,  lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY, P);
}
#endif


#if _FLUID_N2VIBMODEL
double _kappav (np_t *np, long l, gl_t *gl) {
  double kappav;
  double Tv,Pr,T,cp;
  spec_t w;
  Tv=_Tv(np[l]);
  T=_T(np[l],gl);
  find_w(np[l],w);
  cp=_cp_from_w_T(w,T);
  Pr=_eta(np,l,gl)*cp/_kappa(np,l,gl);
  kappav=w[specN2]*(_eta(np,l,gl)/Pr)*_dev_dTv_from_Tv(Tv);  
  return(kappav);
}
#endif


#ifdef _FLUID_FAVREREYNOLDS


double _ktilde (np_t *np, long l, gl_t *gl){
  double ret;
  ret=_k(np[l]);
  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      ret=max(gl->model.fluid.kdiv,_k(np[l]));
    break;
    case TURBMODEL_KOMEGA1988:
      /* etat=0.09*rho*k/psi */
      ret=max(min(gl->model.fluid.kdiv,_psitilde(np[l],gl)*_eta(np,l,gl)/_rho(np[l])*0.09e0),
              _k(np[l]));   
    break;
    case TURBMODEL_KOMEGA2008:
      /* etat approx rho*k/psi */
      ret=max(min(gl->model.fluid.kdiv,_psitilde(np[l],gl)*_eta(np,l,gl)/_rho(np[l])*0.09e0),
              _k(np[l]));   
    break;
    default:
      fatal_error("Turbulence model invalid in _ktilde");
  }
  return(ret);
}


double _psitilde (np_t np, gl_t *gl){
  double ret;
  ret=max(gl->model.fluid.psidiv,_psi(np));
  return(ret);
}


double _eps (np_t np, gl_t *gl){
  double ret;
  ret=0.0;
  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      ret=max(0.0,_psi(np));
    break;
    case TURBMODEL_KOMEGA1988:
      ret=0.09*max(0.0,_k(np))*max(0.0,_psi(np));
    break;
    case TURBMODEL_KOMEGA2008:
      ret=0.09*max(0.0,_k(np))*max(0.0,_psi(np));
    break;
    default:
      fatal_error("Turbulence model invalid in _eps().");
  }
  return(ret);
}


double _omega (np_t *np, long l, gl_t *gl){
  double ret;
  ret=0.0;
  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      ret=_psi(np[l])/_ktilde(np,l,gl);
    break;
    case TURBMODEL_KOMEGA1988:
      ret=_psi(np[l]);
    break;
    case TURBMODEL_KOMEGA2008:
      ret=_psi(np[l]);
    break;
    default:
      fatal_error("Turbulence model invalid in _omega().");
  }
  return(ret);
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


double _etat_from_rho_eta_k_psitilde(np_t *np, long l, gl_t *gl, double rho, double eta, double k, double psitilde){
  double etat,Ret,sum,sum2;
  long i,j,m;
  etat=0.0;
  switch(gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      assert(eta*psitilde!=0.0e0);
      Ret=rho*sqr(max(k,0.0e0))/(eta*psitilde);
      assert((1.0e0+Ret/50.0e0)!=0.0e0);
      assert(psitilde!=0.0e0);
      etat=0.09e0*rho*sqr(max(k,0.0e0))/psitilde
                      *exp(-3.4e0/sqr(1.0e0+Ret/50.0e0));
    break;
    case TURBMODEL_KOMEGA1988:
      etat=rho*max(k,0.0e0)/psitilde;
    break;
    case TURBMODEL_KOMEGA2008:
      sum=0.0;
      for (i=0; i<nd; i++){
        for (j=0; j<nd; j++){
          sum2=_dVi_dxj(np,l,gl,i,j)+_dVi_dxj(np,l,gl,j,i);
          if (i==j) {
            for (m=0; m<nd; m++){
              sum2-=2.0/3.0*_dVi_dxj(np,l,gl,m,m);
            }
          }
          sum+=sqr(sum2);
        }
      } 
      sum=0.0;
      etat=rho*max(k,0.0e0)/max(psitilde,35.0/12.0*sqrt(0.5*sum)); 
    break;
    default:
      fatal_error("Turbulence model invalid in _etat_from_rho_eta_k_psitilde."); 
  }
  etat=max(0.0,etat);
  return(etat);
}


double _etat_mem(np_t *np, long l, gl_t *gl) {
  double etat;
  etat=_etat_from_rho_eta_k_psitilde(np, l, gl, _rho(np[l]), _eta(np,l,gl), _k(np[l]), _psitilde(np[l],gl));
  return(etat);
}


double _etastar (np_t *np, long l, gl_t *gl) {
  double etastar;
  etastar=_eta(np,l,gl)+_etat(np,l,gl);
  return(etastar); 
}


double _kappastar (np_t *np, long l, gl_t *gl) {
  double kappastar;
  double T,cp;
  spec_t w;
  T=_T(np[l],gl);
  find_w(np[l],w);
  cp=_cp_from_w_T(w,T);
  kappastar=_kappa(np,l,gl)+cp*_etat(np,l,gl)/gl->model.fluid.Prt;
  return(kappastar);
}


double _nustar (np_t *np, long l, gl_t *gl, long spec) {
  double nustar;
  nustar=_nu(np[l],gl,spec)+_etat(np,l,gl)/gl->model.fluid.Sct;
  return(nustar);
}


double _etakstar (np_t *np, long l, gl_t *gl) {
  double etakstar;
  etakstar=_eta(np,l,gl);  
  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      etakstar+=_etat(np,l,gl);
    break;
    case TURBMODEL_KOMEGA1988:
      etakstar+=_etat(np,l,gl)/2.0;
    break;
    case TURBMODEL_KOMEGA2008:
      etakstar+=3.0/5.0*_rho(np[l])*max(0.0,_k(np[l]))/_psitilde(np[l],gl);
    break;
    default:
      fatal_error("Turbulence model invalid in _etakstar().");
  }
  return(etakstar);
}


double _etapsistar (np_t *np, long l, gl_t *gl) {
  double etapsistar;
  etapsistar=_eta(np,l,gl);  
  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      etapsistar+=_etat(np,l,gl)/1.3e0;
    break;
    case TURBMODEL_KOMEGA1988:
      etapsistar+=_etat(np,l,gl)/2.0e0;
    break;
    case TURBMODEL_KOMEGA2008:
      etapsistar+=0.5e0*_rho(np[l])*max(0.0,_k(np[l]))/_psitilde(np[l],gl);
    break;
    default:
      fatal_error("Turbulence model invalid in _etapsistar().");
  }
  return(etapsistar);
}


#if _FLUID_N2VIBMODEL
double _kappavstar (np_t *np, long l, gl_t *gl) {
  double kappavstar;
  double Tv,Pr,T,cp;
  spec_t w;
  Tv=_Tv(np[l]);
  T=_T(np[l],gl);
  find_w(np[l],w);
  cp=_cp_from_w_T(w,T);
  Pr=_eta(np,l,gl)*cp/_kappa(np,l,gl);
  kappavstar=w[specN2]*(_eta(np,l,gl)/Pr+_etat(np,l,gl)/gl->model.fluid.Prt)*_dev_dTv_from_Tv(Tv);  
  return(kappavstar);
}
#endif


void find_k_psi_bdry_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, double *k, double *psi){
  double sqrdw;
  long dim;
  *k=max(gl->model.fluid.kmin,1e-99);  //don't make k=0 because this will create problems with positivity-preserving algorithms
  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      *psi=max(gl->model.fluid.psimin,1e-99);   
    break;
    case TURBMODEL_KOMEGA1988:
      sqrdw=0.0e0;
      for (dim=0; dim<nd; dim++){
        sqrdw=sqrdw+sqr(np[lB].bs->x[dim]-np[lA].bs->x[dim]);
      }
      assert_np(np[lB],_rho(np[lB])!=0.0e0);
      assert_np(np[lA],sqrdw!=0.0e0);
      /* use the Menter boundary condition here but substitute eta by etastar -> this works better when
          y+ of the near wall node is significantly higher than 1*/
      *psi=800.0e0*_etastar(np,lB,gl)/_rho(np[lB])/sqrdw;
    break;
    case TURBMODEL_KOMEGA2008:
      sqrdw=0.0e0;
      for (dim=0; dim<nd; dim++){
        sqrdw=sqrdw+sqr(np[lB].bs->x[dim]-np[lA].bs->x[dim]);
      }
      assert_np(np[lB],_rho(np[lB])!=0.0e0);
      assert_np(np[lA],sqrdw!=0.0e0);
      /* use the Menter boundary condition here but substitute eta by etastar -> this works better when
          y+ of the near wall node is significantly higher than 1*/
      *psi=800.0e0*_etastar(np,lB,gl)/_rho(np[lB])/sqrdw;
    break;
    default:
      fatal_error("Turbulence model invalid in find_k_psi_bdry_wall().");
  }

  // Other type of boundary condition for omega
 /* if (gl->model.fluid.TURBMODEL==TURBMODEL_KOMEGA2008) {
    sqrdw=0.0e0;
    sqrqB=0.0e0;
    for (dim=0; dim<nd; dim++){
      sqrdw=sqrdw+sqr(np[lB].bs->x[dim]-np[lA].bs->x[dim]);
      sqrqB=sqrqB+sqr(_V(np[lB],dim));
    }
    ksplus=1.0;
    psiC=40000.0e0/sqr(ksplus)*sqrt(sqrqB/sqrdw);
  }*/

  // Boundary condition for omega with surface roughness
    /* the following is valid as long as ksplus<5 with ksplus=ks*sqrt(tauw*rhow)/muw */
  /* ks is the surface roughness height (here fixed arbitrarily to 1 micron)*/
  /*if (gl->model.fluid.TURBMODEL==TURBMODEL_KOMEGA2008) {
    assert_np(np[lB],_rho(np[lB])!=0.0e0);
    ks=1.0e-6;    
    psiC=40000.0e0*_eta(np[lB])/_rho(np[lB])/sqr(ks);
  }*/

}


/* function returning the value of the turbulence kinetic energy production term */
double _Qk(np_t *np, long l, gl_t *gl){
  long i,j,k;
  double sum,Qk,etaeff;
#ifdef _2D
  long theta,vartheta;
#endif

  if (gl->model.fluid.ADD_ETA_TO_ETAT_WITHIN_QK) {
    etaeff=_etastar(np,l,gl);
  } else {
    etaeff=_etat(np,l,gl);
  }

  Qk=0.0;
  sum=0.0;
  for (k=0; k<nd; k++) sum+=_dVi_dxj(np,l,gl,k,k);
  for (i=0; i<nd; i++){
    for (j=0; j<nd; j++){
      Qk+=_dVi_dxj(np,l,gl,i,j)*(
             etaeff*(_dVi_dxj(np,l,gl,i,j)+_dVi_dxj(np,l,gl,j,i)-2.0/3.0*_delta(i,j)*sum)
             -2.0/3.0*_rho(np[l])*_k(np[l])*_delta(i,j)
          );
    }
  } 

#ifdef _2D
  // these additional axisymmetric terms were derived by Jason Etele and are turned off for now
  if (gl->model.fluid.AXISYMMETRIC && FALSE) {
    sum=0.0e0;
    for (theta=0; theta<nd; theta++) {
      for (vartheta=0; vartheta<nd; vartheta++) {
        sum+=_X(np[l],theta,vartheta)*0.5*(_V(np[_al(gl,l,vartheta,+1)],vartheta)
                                          -_V(np[_al(gl,l,vartheta,-1)],vartheta));
      }
    }
    Qk-=2.0e0/3.0e0*etaeff*_V(np[l],1)/_x(np[l],1)*sum;
  }
#endif
  return(Qk);
}




double _chiw(np_t *np, long l, gl_t *gl){
  double sum,sum2;
  long i,j,k,m;
  sum=0.0;
  for (i=0; i<nd; i++){
    for (j=0; j<nd; j++){
      for (k=0; k<nd; k++){
        sum2=0.0;
        if (k==i) {
          for (m=0; m<nd; m++) sum2+=_dVi_dxj(np,l,gl,m,m);
        }
        sum+=(_dVi_dxj(np,l,gl,i,j)-_dVi_dxj(np,l,gl,j,i))
            *(_dVi_dxj(np,l,gl,j,k)-_dVi_dxj(np,l,gl,k,j))
            *(_dVi_dxj(np,l,gl,k,i)+_dVi_dxj(np,l,gl,i,k)-sum2); 
      }
    }
  }
  sum=fabs(sum*50.0*50.0*50.0/9.0/9.0/9.0/_psitilde(np[l],gl)/_psitilde(np[l],gl)/_psitilde(np[l],gl));
  return(sum);
}


static double _factepsilon(np_t np, gl_t *gl, long theta){
  double Omega,sum,fact;
  long dim;

  Omega=_Omega(np, gl);
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    sum=sum+_X(np,theta,dim)*_X(np,theta,dim);
  }
  fact=Omega*sum;
  return(fact);
}


/* in St, find the source terms of the k-omega / k-epsilon turbulence models */
void find_Stnorm(np_t *np, gl_t *gl, long l, flux_t St){
  long flux,i,j,dim,theta;
  double dkdxj,domegadxj,chiw,eta,etat,Ret,Qk,sum1,sum1a,sum2,sum2a,
         Omega,rho,ktilde,psitilde,k,eps,psi;

  for (flux=0; flux<nf; flux++){
    St[flux]=0.0e0;
  }
  Omega=_Omega(np[l],gl);
  psitilde=_psitilde(np[l],gl);
  psi=_psi(np[l]);
  rho=_rho(np[l]);
  ktilde=_ktilde(np,l,gl);
  k=_k(np[l]);
  eps=_eps(np[l],gl);

  Qk=_Qk(np,l,gl);

  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      eta=_eta(np,l,gl);
      etat=_etat(np,l,gl);
      assert_np(np[l],eta!=0.0e0);
      assert_np(np[l],psitilde!=0.0e0);
      Ret=rho*sqr(max(k,0.0e0))/eta/psitilde;

      sum1=0.0e0;
      sum2=0.0e0;
      for (i=0; i<nd; i++){
       sum1a=0.0e0;
       sum2a=0.0e0;
       for (theta=0; theta<nd; theta++){
         sum1a=sum1a+0.5e0*_X(np[l],theta,i)*(
                +sqrt(fabs(_ktilde(np,_al(gl,l,theta,+1),gl)))
                -sqrt(fabs(_ktilde(np,_al(gl,l,theta,-1),gl))));
           /*the following 5 lines may cause problems. This needs further validation,
             especially in 3D. Will this work independantly of the mesh orientation
             in 2D/3D? */
         if (i!=theta) {
           sum2a=sum2a+sqr(_factepsilon(np[l],gl,theta)/Omega*(
                         +1.0e0*_Vhat(np[_al(gl,l,theta,+1)],i)
                         -2.0e0*_Vhat(np[l],i)
                         +1.0e0*_Vhat(np[_al(gl,l,theta,-1)],i)));
         }
       }
       sum1=sum1+sqr(sum1a);
       sum2=sum2+sum2a;
      }
      St[fluxtke]=(TEST_VANISH*Qk-rho*eps
                               -TEST_VANISH*2.0e0*eta*sum1);
      assert_np(np[l],ktilde!=0.0e0);
      assert_np(np[l],rho!=0.0e0);
      St[fluxpsi]=(psi/ktilde*(
                       +TEST_VANISH*1.44e0*Qk-1.92e0*rho*eps
                       +TEST_VANISH*1.92e0*rho*eps*0.3e0*exp(-sqr(Ret)))
                       +TEST_VANISH*2.0e0*eta*etat/rho*sum2);
    break;
    case TURBMODEL_KOMEGA1988:
      St[fluxtke]=(TEST_VANISH*Qk-rho*eps);
      assert_np(np[l],ktilde!=0.0e0);
      St[fluxpsi]=max(0.0,psi)/ktilde*(TEST_VANISH*5.0e0/9.0e0*Qk-5.0e0/6.0e0*rho*eps);
    break;
    case TURBMODEL_KOMEGA2008:
      St[fluxtke]=(TEST_VANISH*Qk-rho*eps);
      assert_np(np[l],ktilde!=0.0e0);
      chiw=_chiw(np,l,gl);
      sum1=0.0;
      for (j=0; j<nd; j++){
        dkdxj=0.0;
        domegadxj=0.0;
        for (dim=0; dim<nd; dim++){
          dkdxj+=_X(np[l], dim,j)*0.5*(_k(np[_al(gl,l,dim,+1)])-_k(np[_al(gl,l,dim,-1)]));
          domegadxj+=_X(np[l], dim,j)*0.5*(_psi(np[_al(gl,l,dim,+1)])-_psi(np[_al(gl,l,dim,-1)]));
        }
        sum1+=dkdxj*domegadxj;
      }
      St[fluxpsi]=(TEST_VANISH*psi/ktilde*13.0e0/25.0e0*Qk 
                       -psi/ktilde*0.7867*(1.0+85.0*chiw)/(1.0+100.0*chiw)*rho*eps
                       +TEST_VANISH*1.0/8.0*rho/psitilde*max(0.0,sum1));
    break;
    default:
      fatal_error("Turbulence model invalid in find_Stnorm().");
  }
}


static void find_xistar_and_fMt(gl_t *gl, double Mt, double *xi_star, double *fMt){
  double Mt0,fact;
  switch(gl->model.fluid.DILATDISSIP){
    case DILATDISSIP_NONE:
      *xi_star=0.0e0;
      *fMt=0.0e0;
    break;
    case DILATDISSIP_WILCOX:
      *xi_star=1.5e0;
      Mt0=0.25e0;
      if (Mt>Mt0) fact=1.0e0;
        else fact=0.0e0;
      *fMt=(sqr(Mt)-sqr(Mt0))*fact;
    break;
    case DILATDISSIP_SARKAR:
      *xi_star=1.0e0;
      *fMt=sqr(Mt);
    break;
    default:
      *xi_star=0.0e0;
      *fMt=0.0e0;
      fatal_error("Dilatational dissipation model not valid in find_xistar_and_fMt().");
  }
}


/* in St, find the turbulence source terms related to compressibility effects */
void find_Stcomp(np_t *np, gl_t *gl, long l, flux_t St){
  double sum,rho,athermo,Mt;
  long flux,theta,vartheta;
  double ktilde,xi_star,fMt,eps,psi;

  for (flux=0; flux<nf; flux++){
    St[flux]=0.0e0;
  }
  rho=_rho(np[l]);
  eps=_eps(np[l],gl);
  athermo=_athermo(np[l],gl);
#ifdef TEST
  athermo=300.0e0;
#endif
  psi=_psi(np[l]);
  ktilde=_ktilde(np,l,gl);
  if (gl->model.fluid.DILATDISSIP==DILATDISSIP_SARKAR || gl->model.fluid.DILATDISSIP==DILATDISSIP_WILCOX) {
    assert_np(np[l],athermo!=0.0e0);
    Mt=sqrt(2.0e0*max(_k(np[l]),0.0e0))/athermo;
    find_xistar_and_fMt(gl,Mt,&xi_star,&fMt);
    St[fluxtke]=-xi_star*rho*eps*fMt;
    if (gl->model.fluid.TURBMODEL==TURBMODEL_KOMEGA1988 || gl->model.fluid.TURBMODEL==TURBMODEL_KOMEGA2008) {
      St[fluxpsi]=psi/ktilde*xi_star*rho*eps*fMt;
    }
  }
  if (gl->model.fluid.RAPCOMP) {
    sum=0.0e0;
    for (theta=0; theta<nd; theta++){
      for (vartheta=0; vartheta<nd; vartheta++){
        sum=sum+_X(np[l],theta,vartheta)*0.5e0*
                   (_V(np[_al(gl,l,theta,+1)],vartheta)-
                    _V(np[_al(gl,l,theta,-1)],vartheta));
      }
    }
    St[fluxpsi]+=-rho*psi*sum;
  }
}


#ifdef _2D
/* this subroutine is by Jason Etele; the derivation of the terms
   can be found in the doc subdirectory */
void find_Saxi_FavreReynolds(np_t *np, gl_t *gl, long l, flux_t S){
  long flux,species,theta,vartheta,lp,lm;
  double x1, V1;
  double sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12;

  for (flux=0; flux<nf; flux++) S[flux]=0.0e0;

  if (gl->model.fluid.AXISYMMETRIC) {
   x1=np[l].bs->x[1];
   if (fabs(x1)<1e-15) fatal_error("No node must lie on the y=0 axis when AXISYMMETRIC is set to TRUE.");
   V1=_V(np[l],1);
   for (flux=0; flux<nf; flux++){
     S[flux]=(-1.0/x1)*V1*np[l].bs->U[flux];
   }
   // the following terms added by Jason are turned off: they need to be re-derived and their validity
   // checked when x1 is negative
   if (FALSE){
    sum1=0.0e0;
    sum2=0.0e0;
    sum3=0.0e0;
    sum4=0.0e0;
    sum5=0.0e0;
    sum6=0.0e0;
    sum7=0.0e0;
    sum8=0.0e0;
    sum9=0.0e0;
    sum10=0.0e0;
    sum11=0.0e0;
    sum12=0.0e0;
    /* species continuity */
    for (species=0; species<ns; species++) {
      for (theta=0; theta<nd; theta++) {
        lp=_al(gl, l, theta, +1);
        lm=_al(gl, l, theta, -1);
        assert_np(np[lp],np[lp].bs->x[1]>0.0e0);
        assert_np(np[lm],np[lm].bs->x[1]>0.0e0);
        sum1=sum1+_X(np[l],theta,1)*0.5e0*(+_w(np[lp],species)-_w(np[lm],species));
      }
    sum7=sum7+_nustar(np,l,gl,species)*_hk_from_T(species,_T(np[l],gl))*sum1;
    S[species]=S[species]+(1.0/x1)*(_nustar(np,l,gl,species)*sum1);
    }
    /* 2 momentum equations, 1 energy */
    for (theta=0; theta<nd; theta++) {
      lp=_al(gl, l, theta, +1);
      lm=_al(gl, l, theta, -1);
      assert_np(np[lp],_x(np[lp],1)!=0.0e0);
      assert_np(np[lm],_x(np[lm],1)!=0.0e0);
      sum2=sum2+_X(np[l],theta,1)*0.5e0*(+_V(np[lp],0)-_V(np[lm],0));
      sum3=sum3+_X(np[l],theta,0)*0.5e0*(+_V(np[lp],1)-_V(np[lm],1));
      sum4=sum4+_X(np[l],theta,0)*0.5e0*(+_etastar(np,lp,gl)*_V(np[lp],1)/_x(np[lp],1)
					 -_etastar(np,lm,gl)*_V(np[lm],1)/_x(np[lm],1));
      sum5=sum5+_X(np[l],theta,1)*0.5e0*(+_V(np[lp],1)-_V(np[lm],1));
      sum6=sum6+_X(np[l],theta,1)*0.5e0*(+_etastar(np,lp,gl)*_V(np[lp],1)/_x(np[lp],1)
					 -_etastar(np,lm,gl)*_V(np[lm],1)/_x(np[lm],1));
      sum8=sum8+_X(np[l],theta,1)*0.5e0*(+_k(np[lp])-_k(np[lm]));
      sum9=sum9+_X(np[l],theta,0)*0.5e0*(+_V(np[lp],0)-_V(np[lm],0));
      sum10=sum10+_X(np[l],theta,1)*0.5e0*(+_T(np[lp],gl)-_T(np[lm],gl));
      for (vartheta=0; vartheta<nd; vartheta++)  {
	      sum11=sum11+_X(np[l],theta,vartheta)*0.5e0*(+_etastar(np,lp,gl)*_V(np[lp],vartheta)*_V(np[lp],1)/_x(np[lp],1)
						    -_etastar(np,lm,gl)*_V(np[lm],vartheta)*_V(np[lm],1)/_x(np[lm],1));
      }
      sum12=sum12+_X(np[l],theta,1)*0.5e0*(+_psi(np[lp])-_psi(np[lm]));
    }
    S[ns]=S[ns]+(1.0/x1)*(_etastar(np,l,gl)*(sum2+sum3)-(2.0/3.0)*x1*sum4);
    S[ns+1]=S[ns+1]+(1.0/x1)*(2.0*_etastar(np,l,gl)*(sum5-V1/x1) -(2.0/3.0)*x1*sum6);
    S[ns+2]=S[ns+2]+(1.0/x1)*(sum7 + (_etakstar(np,l,gl)*sum8) + _etastar(np,l,gl)*(_V(np[l],0)*(sum2+sum3)+ V1*(4.0/3.0*sum5
				  - 2.0/3.0*sum9) + _kappastar(np,l,gl)*sum10 -(2.0/3.0)*(_etastar(np,l,gl)*V1*V1/x1 + x1*sum11)));
    /* 2 turbulence equations */
    if (gl->model.fluid.TURBSOURCE){
      S[fluxtke]=S[fluxtke]+(1.0/x1)*_etakstar(np,l,gl)*sum8;
      S[fluxpsi]=S[fluxpsi]+(1.0/x1)*_etapsistar(np,l,gl)*sum12;
    }
   }
  }
}
#else
void find_Saxi_FavreReynolds(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  for (flux=0; flux<nf; flux++){
    S[flux]=0.0e0;
  }
}
#endif


void find_dStnorm_dU(np_t *np, gl_t *gl, long l, sqmat_t dStnormdU){
  long spec,row,col;
  double chiw;
  double switch_psiclip;

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dStnormdU[row][col]=0.0e0;
    }
  }
  switch (gl->model.fluid.TURBMODEL){
    case TURBMODEL_KEPSILON:
      dStnormdU[fluxtke][fluxpsi]=-1.0e0;
      dStnormdU[fluxpsi][fluxtke]=1.92e0*sqr(_eps(np[l],gl))/sqr(_ktilde(np,l,gl));
      dStnormdU[fluxpsi][fluxpsi]=-3.84e0*_eps(np[l],gl)/_ktilde(np,l,gl);
    break;
    case TURBMODEL_KOMEGA1988:
      switch_psiclip=1.0e0;
      /*if (_psi(np[l])<=clip.psimin+1e-10) {
        switch_psiclip=0.0;
        wfprintf(stderr,"turning off psi in jacobian\n");
      } */
      for (spec=0; spec<ns; spec++){
        dStnormdU[fluxtke][spec]=_eps(np[l],gl);
        dStnormdU[fluxpsi][spec]=0.09*5.0e0/6.0e0*sqr(max(0.0,_psi(np[l])));
      }
      dStnormdU[fluxtke][fluxtke]=-0.09*max(0.0,_psi(np[l]));
      dStnormdU[fluxtke][fluxpsi]=-0.09*max(0.0,_k(np[l]))*switch_psiclip;
      dStnormdU[fluxpsi][fluxpsi]=-0.09*5.0e0/3.0e0*max(0.0,_psi(np[l]))*switch_psiclip;
    break;
    case TURBMODEL_KOMEGA2008:
      chiw=_chiw(np,l,gl);
      for (spec=0; spec<ns; spec++){
        dStnormdU[fluxtke][spec]=_eps(np[l],gl);
        dStnormdU[fluxpsi][spec]=0.0708*(1.0+85.0*chiw)/(1.0+100.0*chiw)*sqr(_psi(np[l]));
      }
      dStnormdU[fluxtke][fluxtke]=-0.09*_psi(np[l]);
      dStnormdU[fluxtke][fluxpsi]=-0.09*_k(np[l]);
      dStnormdU[fluxpsi][fluxpsi]=-max(0.0,2.0*0.09*0.7867*(1.0+85.0*chiw)/(1.0+100.0*chiw))*_psi(np[l])
         -0.09*_psi(np[l])*max(0.0,2360.1*chiw*(1.0+85.0*chiw)/sqr(1.0+100.0*chiw))
         -0.09*_psi(np[l])*max(0.0,-200.61*chiw/(1.0+100.0*chiw));
    break;
    default:
      fatal_error("Turbulence model invalid in find_dStnorm_dU().");
  }
}


void find_dStcomp_dU(np_t *np, gl_t *gl, long l, sqmat_t dStcompdU){
  double Mt,athermo,ktilde,k;
  long theta,vartheta,row,col,spec;
  double sum,xi_star,fMt;

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dStcompdU[row][col]=0.0e0;
    }
  }
  if (gl->model.fluid.RAPCOMP) {
    sum=0.0e0;
    for (theta=0; theta<nd; theta++){
      for (vartheta=0; vartheta<nd; vartheta++){
        sum=sum+_X(np[l],theta,vartheta)*0.5e0*
                   (_V(np[_al(gl,l,theta,+1)],vartheta)-
                    _V(np[_al(gl,l,theta,-1)],vartheta));
      }
    }
    dStcompdU[fluxpsi][fluxpsi]=-sum;
  }
  if (gl->model.fluid.DILATDISSIP==DILATDISSIP_SARKAR || gl->model.fluid.DILATDISSIP==DILATDISSIP_WILCOX) {
    /*  note: the following assumes athermo to be constant */
    athermo=_athermo(np[l],gl);
#ifdef TEST
    athermo=300.0e0;
#endif
    ktilde=_ktilde(np,l,gl);
    k=_k(np[l]);
    assert_np(np[l],athermo!=0.0e0);
    Mt=sqrt(2.0e0*max(k,0.0e0)/sqr(athermo));
    find_xistar_and_fMt(gl, Mt, &xi_star, &fMt);

    switch (gl->model.fluid.TURBMODEL){
      case TURBMODEL_KEPSILON:
        dStcompdU[fluxtke][fluxpsi]=-xi_star*fMt;
        dStcompdU[fluxtke][fluxtke]=-xi_star*_psi(np[l])*fMt/ktilde;
        for (spec=0; spec<ns; spec++){
          dStcompdU[fluxtke][spec]=xi_star*_psi(np[l])*fMt;
        }
      break;
      case TURBMODEL_KOMEGA1988:
        dStcompdU[fluxtke][fluxpsi]=-xi_star*fMt*k*0.09;
        dStcompdU[fluxtke][fluxtke]=-xi_star*fMt*2.0e0*_psi(np[l])*0.09;
        for (spec=0; spec<ns; spec++){
          dStcompdU[fluxtke][spec]=xi_star*2.0e0*k*_psi(np[l])*fMt*0.09;
        }
      break;
      case TURBMODEL_KOMEGA2008:
/*      St[fluxtke]=-xi_star*rho*eps*fMt;
        St[fluxtke]=-xi_star*(rho*k)*(rho*psi)/rho*fMt*0.09;
        need to check about the factor 2.0e0 below */
        dStcompdU[fluxtke][fluxpsi]=-xi_star*fMt*k*0.09;
        dStcompdU[fluxtke][fluxtke]=-xi_star*fMt*2.0e0*_psi(np[l])*0.09;
        for (spec=0; spec<ns; spec++){
          dStcompdU[fluxtke][spec]=xi_star*2.0e0*k*_psi(np[l])*fMt*0.09;
        }
      break;
      default:
        fatal_error("Turbulence model invalid in find_dStcomp_dU().");
    }
  }
}



#endif


double _Lc(metrics_t metrics){
  double Lc,Xmag;
  long dim;
  Xmag=0.0;
  for (dim=0; dim<nd; dim++) Xmag+=sqr(metrics.X[dim]);
  Xmag=sqrt(Xmag);
  Lc=1.0/Xmag;
  return(Lc);
}


double _astar_from_jacvars(jacvars_t jacvars, metrics_t metrics){
  double astar;
  astar=_a_from_jacvars(jacvars)/_Lc(metrics);
  return(astar);
}


double _Vstar_from_jacvars(jacvars_t jacvars, metrics_t metrics){
  double Vstar;
  long dim;
  Vstar=0.0;
  for (dim=0; dim<nd; dim++) Vstar+=_V_from_jacvars(jacvars,dim)*metrics.X[dim];

  return(Vstar);
}


#ifdef _FLUID_NEUTRALSTRANSPORT
void condition_Lambda_absolute_Harten(jacvars_t jacvars, metrics_t metrics, double zeta, bool SPECTRALRADIUS, sqmat_t Lambdaabs){
  long flux,dim;
  double Xmag2,alocal,a,sigma;
  // sigma is the spectral radius
  sigma=0.0;
  for (flux=0; flux<nf; flux++){
    sigma=max(sigma,Lambdaabs[flux][flux]);
  }
  a=_a_from_jacvars(jacvars);
  Xmag2=0.0;
  for (dim=0; dim<nd; dim++) {
      Xmag2+=sqr(metrics.X[dim]);
  }
  a=_a_from_jacvars(jacvars)*sqrt(Xmag2);

  if (SPECTRALRADIUS){
    alocal=sigma;
  } else {
    alocal=a;
  }
  for (flux=0; flux<nf; flux++){
    Lambdaabs[flux][flux]+=zeta*alocal;
  }
}


void condition_Lambda_plus_minus_Harten_very_old(jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  long flux;
  double sigma,zetaA1;
  sigma=0.0;
  for (flux=0; flux<nf; flux++) sigma=max(sigma,0.5*(Lambdaplus[flux][flux]-Lambdaminus[flux][flux]));
  zetaA1=jacvarsp0.zetaA1;
  //sigma+=(_a_from_jacvars(jacvarsp0)+_a_from_jacvars(jacvarsp1));
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]+=sigma*zetaA1;
    Lambdaminus[flux][flux]-=sigma*zetaA1;
  }
}


void condition_Lambda_plus_minus_Harten_old(jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  long flux;
  double zetaA1,ap0,ap1;
  flux_t sigma;
  for (flux=0; flux<nf; flux++) sigma[flux]=0.5*(Lambdaplus[flux][flux]-Lambdaminus[flux][flux]);
  zetaA1=jacvarsp0.zetaA1;
  ap0=_astar_from_jacvars(jacvarsp0,metrics);
  ap1=_astar_from_jacvars(jacvarsp1,metrics);
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]+=max(sigma[flux],ap0)*zetaA1;
    Lambdaminus[flux][flux]-=max(sigma[flux],ap1)*zetaA1;
  }
}


void condition_Lambda_plus_minus_Harten(jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  long flux;
  double zetaA1,ap0,ap1,Vp0,Vp1,aref,Vref;

  zetaA1=jacvarsp0.zetaA1;
  ap0=_astar_from_jacvars(jacvarsp0,metrics);
  ap1=_astar_from_jacvars(jacvarsp1,metrics);
  Vp0=_Vstar_from_jacvars(jacvarsp0, metrics);
  Vp1=_Vstar_from_jacvars(jacvarsp1, metrics);

  /* first, make sure the eigenvalues are within physically-realistic bounds and clip them if necessary */
  aref=max(ap0,ap1);
  Vref=max(fabs(Vp0),fabs(Vp1));
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]=min(Lambdaplus[flux][flux],Vref+aref/notzero(jacvarsp0.zetaA2,1e-99));
    Lambdaminus[flux][flux]=max(Lambdaminus[flux][flux],-Vref-aref/notzero(jacvarsp0.zetaA2,1e-99));
  }

  /* second, add entropy correction to prevent carbuncle */
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]+=(Vref+aref)*zetaA1;
    Lambdaminus[flux][flux]-=(Vref+aref)*zetaA1;
  }
}


void condition_Lambda_absolute_Peclet(jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaabs){
  long flux,dim,dim2;
  double a,Xmag2,Xmag2min,sum,Pe,zetaAc,zetaAa;
  Xmag2=0.0;
  Xmag2min=1e99;
  for (dim=0; dim<nd; dim++) {
      Xmag2+=sqr(metrics.X[dim]);
      sum=0.0;
      for (dim2=0; dim2<nd; dim2++) sum+=sqr(metrics.X2[dim][dim2]);
      Xmag2min=min(sum,Xmag2min);
  }
  a=_a_from_jacvars(jacvars)*sqrt(Xmag2);
  Pe=_Pe_from_jacvars(jacvars, metrics);
  zetaAc=jacvars.zetaA1*max(0.0,sqrt(Xmag2min/Xmag2)-1.0/max(1.0e-10,min(PECLETMAX,Pe)));
  zetaAa=jacvars.zetaA1;
  for (flux=0; flux<nf; flux++){
    if (FALSE) Lambdaabs[flux][flux]+=zetaAa*a; else Lambdaabs[flux][flux]+=zetaAc*a; 
  }
}


void condition_Lambda_plus_minus_Peclet(np_t *np, gl_t *gl, long lp0, long theta,  jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  long flux,dim,dim2;
  double sum,Xmag2min,Pe,zetaA1,Xmag2;
  double aref,Vref,ap0,ap1,Vp0,Vp1,factPe;

  zetaA1=jacvarsp0.zetaA1;
  ap0=_astar_from_jacvars(jacvarsp0,metrics);
  ap1=_astar_from_jacvars(jacvarsp1,metrics);
  Vp0=_Vstar_from_jacvars(jacvarsp0, metrics);
  Vp1=_Vstar_from_jacvars(jacvarsp1, metrics);
  
  aref=max(ap0,ap1);
  Vref=max(fabs(Vp0),fabs(Vp1));
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]=min(Lambdaplus[flux][flux],Vref+aref/notzero(jacvarsp0.zetaA2,1e-99));
    Lambdaminus[flux][flux]=max(Lambdaminus[flux][flux],-Vref-aref/notzero(jacvarsp0.zetaA2,1e-99));
  }

  Xmag2min=1e99;
  Xmag2=0.0;
  for (dim=0; dim<nd; dim++) {
    Xmag2+=sqr(metrics.X[dim]);
    sum=0.0;
    for (dim2=0; dim2<nd; dim2++) sum+=sqr(metrics.X2[dim][dim2]);
    Xmag2min=min(sum,Xmag2min);
  }
  zetaA1=jacvarsp0.zetaA1;
  Pe=0.5*(_Pe_from_jacvars(jacvarsp0, metrics)+_Pe_from_jacvars(jacvarsp1, metrics));  
  factPe=max(0.0,sqrt(Xmag2min/Xmag2)-1.0/max(1.0e-10,min(PECLETMAX,Pe)));
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]+=(aref+Vref)*zetaA1*factPe;
    Lambdaminus[flux][flux]-=(aref+Vref)*zetaA1*factPe;
  }
}




void condition_Lambda_plus_minus_Pascal(np_t *np, gl_t *gl, long lp0, long theta, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  long flux,dim,lp1;
  double zetaA1;
  double ap0,ap1,Vp0,Vp1,Vref,aref,factP,Pstarmax,Pstarmin;

  zetaA1=jacvarsp0.zetaA1;
  ap0=_astar_from_jacvars(jacvarsp0,metrics);
  ap1=_astar_from_jacvars(jacvarsp1,metrics);
  Vp0=_Vstar_from_jacvars(jacvarsp0, metrics);
  Vp1=_Vstar_from_jacvars(jacvarsp1, metrics);

  /* first, make sure the eigenvalues are within physically-realistic bounds and clip them if necessary */
  aref=max(ap0,ap1);
  Vref=max(fabs(Vp0),fabs(Vp1));
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]=min(Lambdaplus[flux][flux],Vref+aref/notzero(jacvarsp0.zetaA2,1e-99));
    Lambdaminus[flux][flux]=max(Lambdaminus[flux][flux],-Vref-aref/notzero(jacvarsp0.zetaA2,1e-99));
  }

  /* second, add entropy correction to prevent carbuncle */
  lp1=_al(gl,lp0,theta,+1);
  Pstarmax=max(_Pstar(np[lp0],gl),_Pstar(np[lp1],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta){
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp0,dim,+1)],gl));
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp0,dim,-1)],gl));
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp1,dim,+1)],gl));
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp1,dim,-1)],gl));
    }
  }  
  Pstarmin=min(_Pstar(np[lp0],gl),_Pstar(np[lp1],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta){
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp0,dim,+1)],gl));
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp0,dim,-1)],gl));
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp1,dim,+1)],gl));
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp1,dim,-1)],gl));
    }
  }  
  Pstarmin=max(1e-40,Pstarmin);
  Pstarmax=max(1e-40,Pstarmax);

  zetaA1=jacvarsp0.zetaA1;
  factP=max(0.0,pow(fabs(Pstarmax-Pstarmin)/Pstarmin,1.0));
//  factP=max(0.0,pow(fabs(Pstarmax-Pstarmin)/Pstarmin,1.0)-0.3);
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]+=(aref+Vref)*zetaA1*factP;
    Lambdaminus[flux][flux]-=(aref+Vref)*zetaA1*factP;
  }
}


void condition_Lambda_plus_minus_Parent(np_t *np, gl_t *gl, long lp0, long theta, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  long flux,spec,dim,lp1;
  double ap0,ap1,Vp0,Vp1,Vref,aref,factP,Pstarmax,Pstarmin,factcommon;
  double lambdaadd1,lambdaadd2,lambdaadd,alphamin,alphamax,alphamass;
  flux_t LUstarp0,LUstarp1,alpha,fact;
  sqmat_t Lambdap0,Lambdap1;
  bool FACTCOMMON[nf];


  /* FIRST, make sure the eigenvalues of the waves are within same order of magnitude as the spectral radius and, 
     if not, blend them with Steger-Warming Eigenvalues    */
  ap0=_astar_from_jacvars(jacvarsp0,metrics);
  ap1=_astar_from_jacvars(jacvarsp1,metrics);
  Vp0=_Vstar_from_jacvars(jacvarsp0, metrics);
  Vp1=_Vstar_from_jacvars(jacvarsp1, metrics);
  find_Lambda_from_jacvars(jacvarsp0, metrics, Lambdap0);
  find_Lambda_from_jacvars(jacvarsp1, metrics, Lambdap1);
  aref=max(ap0,ap1);
  Vref=max(fabs(Vp0),fabs(Vp1));
  for (flux=0; flux<nf; flux++) FACTCOMMON[flux]=TRUE;
#ifdef _FLUID_PLASMA
  for (flux=0; flux<nf; flux++) {
    if (flux<ncs) FACTCOMMON[flux]=FALSE;
  }
#endif
  factcommon=1.0;
  for (flux=0; flux<nf; flux++){
    fact[flux]=min(1.0,(Vref+aref/notzero(jacvarsp0.zetaA2,1e-99))/notzero(max(Lambdaplus[flux][flux],-Lambdaminus[flux][flux]),1e-99));
    if (FACTCOMMON[flux]) factcommon=min(factcommon,fact[flux]);
  }
  for (flux=0; flux<nf; flux++){
    if (FACTCOMMON[flux]) fact[flux]=factcommon;
  }
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]=(1.0-fact[flux])*max(0.0,Lambdap0[flux][flux])+fact[flux]*Lambdaplus[flux][flux];
    Lambdaminus[flux][flux]=(1.0-fact[flux])*min(0.0,Lambdap1[flux][flux])+fact[flux]*Lambdaminus[flux][flux];
  }

  /* SECOND, add eigenvalue conditioning in vicinity of pressure gradients to prevent carbuncle */
  lp1=_al(gl,lp0,theta,+1);
  Pstarmax=max(_Pstar(np[lp0],gl),_Pstar(np[lp1],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta){
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp0,dim,+1)],gl));
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp0,dim,-1)],gl));
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp1,dim,+1)],gl));
      Pstarmax=max(Pstarmax,_Pstar(np[_al(gl,lp1,dim,-1)],gl));
    }
  }  
  Pstarmin=min(_Pstar(np[lp0],gl),_Pstar(np[lp1],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta){
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp0,dim,+1)],gl));
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp0,dim,-1)],gl));
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp1,dim,+1)],gl));
      Pstarmin=min(Pstarmin,_Pstar(np[_al(gl,lp1,dim,-1)],gl));
    }
  }  
  Pstarmin=max(1e-40,Pstarmin);
  Pstarmax=max(1e-40,Pstarmax);

//  factP=max(0.0,pow(fabs(Pstarmax-Pstarmin)/Pstarmin,1.0));
  factP=max(0.0,pow(fabs(Pstarmax-Pstarmin)/Pstarmin,1.0)-0.3);
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]+=(aref+Vref)*jacvarsp0.zetaA1*factP;
    Lambdaminus[flux][flux]-=(aref+Vref)*jacvarsp0.zetaA1*factP;
  }


  /* THIRD, make sure the flux is positivity-preserving in multiple dimensions by conditioning the eigenvalues */
  find_LUstar_from_jacvars(jacvarsp1, metrics, LUstarp1);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUstarp0);

  for (flux=0; flux<nf; flux++) alpha[flux]=Lambdaplus[flux][flux]*LUstarp0[flux]/min(LUstarp0[fluxet],LUstarp0[fluxet-1]);
  alphamass=0.0;
  for (spec=0; spec<ns; spec++) alphamass+=alpha[spec];
  alphamax=alphamass;
  alphamin=alphamass;
  for (flux=fluxmom; flux<fluxet-1; flux++) {
    alphamax=max(alphamax,alpha[flux]);
    alphamin=min(alphamin,alpha[flux]);
  }
  lambdaadd1=max(0.0,alphamax-alphamin-(alpha[fluxet]+alpha[fluxet-1]));

  for (flux=0; flux<nf; flux++) alpha[flux]=-Lambdaminus[flux][flux]*LUstarp1[flux]/min(LUstarp1[fluxet],LUstarp1[fluxet-1]);
  alphamass=0.0;
  for (spec=0; spec<ns; spec++) alphamass+=alpha[spec];
  alphamax=alphamass;
  alphamin=alphamass;
  for (flux=fluxmom; flux<fluxet-1; flux++) {
    alphamax=max(alphamax,alpha[flux]);
    alphamin=min(alphamin,alpha[flux]);
  }
  lambdaadd2=max(0.0,alphamax-alphamin-(alpha[fluxet]+alpha[fluxet-1]));
  lambdaadd=0.5*max(lambdaadd1,lambdaadd2);
  for (flux=fluxet-1; flux<=fluxet; flux++){
    Lambdaplus[flux][flux]+=lambdaadd;
    Lambdaminus[flux][flux]-=lambdaadd;
  }

}


// as defined in "Computational Aerothermodynamic Simulation Issues on Unstructured Grids" by Gnoffo and White, AIAA 2004-2371
double _Re_edge_from_jacvars(jacvars_t jacvars, metrics_t metrics){
  double Re_edge,a,Xmag2,eta;
  long dim;

  a=_a_from_jacvars(jacvars);
  eta=_eta_from_jacvars(jacvars);
  Xmag2=0.0;
  for (dim=0; dim<nd; dim++) {
    Xmag2+=sqr(metrics.X[dim]);
  }
  Re_edge=jacvars.rho/max(1.0e-20,eta)*a/sqrt(Xmag2); 
  return(Re_edge);
}


/* the following is taken from "Computational Aerothermodynamic Simulation Issues on Unstructured Grids" by Gnoffo and White, AIAA 2004-2371
*/
void condition_Lambda_absolute_Gnoffo(jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaabs){
  long flux;
  double sigma,lambdaref,Re_edge;
  sigma=0.0;
  for (flux=0; flux<nf; flux++){
    sigma=max(sigma,Lambdaabs[flux][flux]);
  }
  Re_edge=_Re_edge_from_jacvars(jacvars,metrics);
  lambdaref=jacvars.zetaA1*min(1.0,pow(Re_edge/100.0,8.0))*sigma;
  for (flux=0; flux<nf; flux++){
    if (Lambdaabs[flux][flux]<2.0*lambdaref){
      Lambdaabs[flux][flux]=sqr(Lambdaabs[flux][flux])/4.0/(1.0e-20+lambdaref)+lambdaref;
    }
  }
}


void condition_Lambda_plus_minus_Gnoffo(jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  long flux;
  double ap0,ap1,Vp0,Vp1,aref,Vref,zetaA1,Re_edge;

  zetaA1=jacvarsp0.zetaA1;
  ap0=_astar_from_jacvars(jacvarsp0,metrics);
  ap1=_astar_from_jacvars(jacvarsp1,metrics);
  Vp0=_Vstar_from_jacvars(jacvarsp0, metrics);
  Vp1=_Vstar_from_jacvars(jacvarsp1, metrics);

  /* first, make sure the eigenvalues are within physically-realistic bounds and clip them if necessary */
  aref=max(ap0,ap1);
  Vref=max(fabs(Vp0),fabs(Vp1));
  for (flux=0; flux<nf; flux++){
    Lambdaplus[flux][flux]=min(Lambdaplus[flux][flux],Vref+aref/notzero(jacvarsp0.zetaA2,1e-99));
    Lambdaminus[flux][flux]=max(Lambdaminus[flux][flux],-Vref-aref/notzero(jacvarsp0.zetaA2,1e-99));
  }

  Re_edge=0.5*(_Re_edge_from_jacvars(jacvarsp0,metrics)+_Re_edge_from_jacvars(jacvarsp1,metrics));
  zetaA1*=min(1.0,pow(Re_edge/100.0,8.0));

  for (flux=0; flux<nf; flux++){  
    Lambdaplus[flux][flux]+=zetaA1*(aref+Vref);
    Lambdaminus[flux][flux]-=zetaA1*(aref+Vref);
  }
}



void find_conditioned_Lambda_absolute_from_jacvars_denom(jacvars_t jacvars, metrics_t metrics, sqmat_t Lambda){
  long flux;
  find_Lambda_from_jacvars(jacvars,metrics,Lambda);
  for (flux=0; flux<nf; flux++) Lambda[flux][flux]=max(1.0e-20,fabs(Lambda[flux][flux]));

  condition_Lambda_absolute_Harten(jacvars,metrics,jacvars.zetaA3,TRUE,Lambda);
}


void find_conditioned_Lambda_absolute_from_jacvars(jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t Lambdaabs){
  long flux;
  find_Lambda_from_jacvars(jacvars,metrics,Lambdaabs);
  if (EIGENVALCOND==EIGENVALCOND_DEFAULT) EIGENVALCOND=jacvars.EIGENVALCOND;

  for (flux=0; flux<nf; flux++) Lambdaabs[flux][flux]=fabs(Lambdaabs[flux][flux]);
  switch (EIGENVALCOND){
    case EIGENVALCOND_PECLET:
      condition_Lambda_absolute_Peclet(jacvars, metrics, Lambdaabs);
    break;
    case EIGENVALCOND_PASCAL:
      fatal_error("Eigenvalue conditioning can not be set to EIGENVALCOND_PASCAL with this flux discretization scheme.");
    break;
    case EIGENVALCOND_PARENT:
      fatal_error("Eigenvalue conditioning can not be set to EIGENVALCOND_PARENT with this flux discretization scheme.");
    break;
    case EIGENVALCOND_GNOFFO:
      condition_Lambda_absolute_Gnoffo(jacvars, metrics, Lambdaabs);
    break;
    case EIGENVALCOND_HARTEN:
      condition_Lambda_absolute_Harten(jacvars,metrics,jacvars.zetaA1, FALSE, Lambdaabs);
    break;
    case EIGENVALCOND_NONE:
    break;
    default:
      fatal_error("EIGENVALCOND=%ld is not valid in find_conditioned_Lambda_absolute_from_jacvars().",EIGENVALCOND);
  }
}


void condition_Lambda_plus_minus(np_t *np, gl_t *gl, long lp0, long theta, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  int EIGENVALCOND, sqmat_t Lambdaplus, sqmat_t Lambdaminus){
  /* make sure zetaA2 is within valid bounds */
  if (jacvarsp0.zetaA2<0.0) fatal_error("zetaA2 can not be negative.");
  if (jacvarsp0.zetaA2>1.0) fatal_error("zetaA2 can not be greater than 1.");
  if (EIGENVALCOND==EIGENVALCOND_DEFAULT) EIGENVALCOND=gl->cycle.resconv.EIGENVALCOND;

  switch (EIGENVALCOND){
    case EIGENVALCOND_PECLET:
      condition_Lambda_plus_minus_Peclet(np, gl, lp0, theta, jacvarsp0, jacvarsp1, metrics,  Lambdaplus, Lambdaminus);
    break;
    case EIGENVALCOND_PASCAL:
      condition_Lambda_plus_minus_Pascal(np, gl, lp0, theta, jacvarsp0, jacvarsp1, metrics,  Lambdaplus, Lambdaminus);
    break;
    case EIGENVALCOND_PARENT:
      condition_Lambda_plus_minus_Parent(np, gl, lp0, theta, jacvarsp0, jacvarsp1, metrics,  Lambdaplus, Lambdaminus);
    break;
    case EIGENVALCOND_GNOFFO:
      condition_Lambda_plus_minus_Gnoffo(jacvarsp0, jacvarsp1, metrics,  Lambdaplus, Lambdaminus);
    break;
    case EIGENVALCOND_HARTEN:
      condition_Lambda_plus_minus_Harten(jacvarsp0, jacvarsp1, metrics,  Lambdaplus, Lambdaminus);
    break;
    case EIGENVALCOND_NONE:
    break;
    default:
      fatal_error("EIGENVALCOND=%ld is not valid in condition_Lambda_plus_minus().",EIGENVALCOND);
  }
}
#else //not defined _FLUID_NEUTRALSTRANSPORT

void find_conditioned_Lambda_absolute_from_jacvars(jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t Lambdaabs){
  long flux;
  find_Lambda_from_jacvars(jacvars,metrics,Lambdaabs);

  for (flux=0; flux<nf; flux++) Lambdaabs[flux][flux]=fabs(Lambdaabs[flux][flux]);
}

void condition_Lambda_plus_minus(np_t *np, gl_t *gl, long lp0, long theta, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  int EIGENVALCOND, sqmat_t Lambdaplus, sqmat_t Lambdaminus){
}

#endif //_FLUID_NEUTRALSTRANSPORT

void rearrange_metrics_eigenset2(metrics_t *metrics){

#ifdef _3D
  double tol,X1,X2,X3,Xmax;
  long Xmaxindex;

  tol=5e-2;

  X1=metrics->X[0];
  X2=metrics->X[1];
  X3=metrics->X[2];
  
  Xmax=max(fabs(X1),max(fabs(X2),fabs(X3)));
  assert(Xmax!=0.0);
  Xmaxindex=0;
  if (fabs(X1)>=fabs(X2) && fabs(X1)>=fabs(X3)) Xmaxindex=0;
  if (fabs(X2)>=fabs(X1) && fabs(X2)>=fabs(X3)) Xmaxindex=1;
  if (fabs(X3)>=fabs(X1) && fabs(X3)>=fabs(X2)) Xmaxindex=2;
  if (fabs(X1+X2+X3)<Xmax*tol) {
     metrics->X[Xmaxindex]=metrics->X[Xmaxindex]*(1.0+tol-fabs(X1+X2+X3)/(Xmax));
  }
#endif
}


void find_Gamma(np_t *np, gl_t *gl, long l, sqmat_t Gamma){
  long row,col;
  double dtau_constant,dtau_local;
#ifdef _FLUID_PLASMA
  spec2_t alpha;
  double dtaufact;
#endif

  flux_t dtau_vector;
  find_constant_dtau(np, gl, l, &dtau_constant);
  find_dtau(np,gl,l,dtau_vector);
#ifdef _FLUID_PLASMA
  find_alpha(np, gl, l, alpha);
#endif
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Gamma[row][col]=0.0e0;
      if (row==col) {
        switch (gl->PRECONDITIONER){
          case PRECON_LOCALTIMESTEP:
            Gamma[row][row]=1.0/dtau_constant;
          break;
          case PRECON_LOCALTIMESTEP2:
            Gamma[row][row]=1.0/dtau_constant;
          break;
          case PRECON_CONSTANTTIMESTEP:
            Gamma[row][row]=1.0/gl->dtau;
          break;
          case PRECON_LOCALEIGENVALUE: /* as in Parent, Positivity-Preserving Dual-Time Stepping Schemes for Gas Dynamics, JCP 2018. */
            dtau_local=max(0.01*dtau_constant,min(100.0*dtau_constant,dtau_vector[row]));
            Gamma[row][row]=1.0e0/sqrt(dtau_constant*dtau_local);
          break;
          case PRECON_LOCALEIGENVALUE2: /* new way, can be more stable in some cases */
            if (row>=fluxet-nd && row<=fluxet) {
              dtau_local=dtau_constant;
            } else {
              dtau_local=max(0.01*dtau_constant,min(100.0*dtau_constant,dtau_vector[row]));
            }
            Gamma[row][row]=1.0e0/sqrt(dtau_constant*dtau_local);
          break;
          default:
            fatal_error("PRECONDITIONER %d not supported.",gl->PRECONDITIONER); 
        }
#ifdef _FLUID_PLASMA
        if (row<ncs){
          dtaufact=1.0;
          if (row!=speceminus) dtaufact=gl->cycle.fluid.zetaGammai;
            else dtaufact=min(gl->cycle.fluid.zetaGammai,gl->cycle.fluid.zetaGammae/alpha[speceminus][speceminus]);
          
          Gamma[row][row]*=1.0/dtaufact;
        }
#ifdef _FLUID_DRIFTDIFFUSION
        /* keep the electron temperature the same during the update */
        if (row==nf-1) Gamma[row][row]/=1e-20;
#endif
#if _FLUID_EENERGY
        if (gl->model.fluid.TEMODEL==TEMODEL_LOCAL && row==fluxee) Gamma[row][row]=1.0/(dtau_constant*1e-99);
#endif
#endif
      }
    }
  }
}


void find_Gamma_inverse(np_t *np, gl_t *gl, long l, sqmat_t GammaInv){
  sqmat_t Gamma;
  find_Gamma(np,gl,l,Gamma);
  invert_diagonal_matrix(Gamma, GammaInv);
}


#if _FLUID_CONVECTION
void set_jacvars_eigenconditioning_constants(gl_t *gl, jacvars_t *jacvars){
  /* the zeta prime artificial dissipation term.. */
  jacvars->zetaA1=gl->model.fluid.zetaA1;
  jacvars->zetaA3=gl->model.fluid.zetaA3;
  jacvars->zetaA2=gl->model.fluid.zetaA2;
  jacvars->EIGENVALCOND=gl->cycle.resconv.EIGENVALCOND;
}
#endif


#ifdef _FLUID_PLASMA

double _betan(long k){
  double betan;
  if (speciestype[k]==SPECIES_NEUTRAL) betan=1.0; else betan=0.0;
  return(betan);
}

double _betai(long k){
  double betai;
  if (speciestype[k]==SPECIES_IONMINUS || speciestype[k]==SPECIES_IONPLUS) betai=1.0; else betai=0.0;
  return(betai);
}

double _betac(long k){
  double betac;
  if (speciestype[k]!=SPECIES_NEUTRAL) betac=1.0; else betac=0.0;
  return(betac);
}


/* determines whether the gauss's law terms are added to the continuity equations */
double _betag(gl_t *gl, long k){  
  return(gl->model.fluid.betag[k]);
}


double _betaplus(long k){
  double betaplus;
  if (speciestype[k]==SPECIES_IONPLUS) betaplus=1.0; else betaplus=0.0;
  return(betaplus);
}


double _betaa(gl_t *gl, long k){
  return(gl->model.fluid.betaa[k]);
}


double _rhoc(np_t np, gl_t *gl){
  double rhoc;
  long spec;
  rhoc=0.0;
  for (spec=0; spec<ns; spec++){
    rhoc+=_C(spec)*_Nk(np, gl, spec);
  }
  if (CHEM_NEUTRAL) rhoc=0.0;

  return(rhoc);
}


#define Ethresholdbdry 1e-2

bool is_bdry_cathode(np_t *np, gl_t *gl, long lA, long lB, long lC,
                     long theta, long thetasgn){
  bool CATHODE;
  double dphi,dwall,fact_cathode,Nminus,Nplus;
  long dim,spec;

  CATHODE=FALSE;

  dwall=0.0;
  for (dim=0; dim<nd; dim++){
    dwall=dwall+sqr(np[lB].bs->x[dim]-np[lA].bs->x[dim]);
  }
  assert_np(np[lA],dwall>0.0e0);
  dwall=sqrt(dwall);


  Nminus=0.0;
  Nplus=0.0;
  for (spec=0; spec<ns; spec++){
    if (_Charge_number(spec)<0) Nminus+=_Nk(np[lA],gl,spec);
    if (_Charge_number(spec)>0) Nplus+=_Nk(np[lA],gl,spec);
  }
  if (Nplus>Nminus) fact_cathode=1.0; else fact_cathode=0.0;

// Approach 1
/*
  dphi=minmod(_phi(np[lA],gl)-_phi(np[lB],gl),_phi(np[lB],gl)-_phi(np[lC],gl));
*/


// Approach 2
  dphi=2.0*(_phi(np[lA],gl)-_phi(np[lB],gl));
  dphi+=_phi(np[lB],gl)-_phi(np[lC],gl);
  for (dim=0; dim<nd; dim++){
    if (dim!=theta) dphi+=1.0*(_phi(np[_al(gl,lA,dim,+1)],gl)-_phi(np[_al(gl,lB,dim,+1)],gl));
    if (dim!=theta) dphi+=1.0*(_phi(np[_al(gl,lA,dim,-1)],gl)-_phi(np[_al(gl,lB,dim,-1)],gl));
  }
  dphi/=(3.0+((double)nd-1.0)*2.0);


// Approach 3

/*  double dphi1,dphi2;
  dphi1=2.0*(_phi(np[lA],gl)-_phi(np[lB],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta) dphi1+=1.0*(_phi(np[_al(gl,lA,dim,+1)],gl)-_phi(np[_al(gl,lB,dim,+1)],gl));
    if (dim!=theta) dphi1+=1.0*(_phi(np[_al(gl,lA,dim,-1)],gl)-_phi(np[_al(gl,lB,dim,-1)],gl));
  }
  dphi1/=(2.0+((double)nd-1.0)*2.0);

  dphi2=2.0*(_phi(np[lB],gl)-_phi(np[lC],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta) dphi2+=1.0*(_phi(np[_al(gl,lB,dim,+1)],gl)-_phi(np[_al(gl,lC,dim,+1)],gl));
    if (dim!=theta) dphi2+=1.0*(_phi(np[_al(gl,lB,dim,-1)],gl)-_phi(np[_al(gl,lC,dim,-1)],gl));
  }
  dphi2/=(2.0+((double)nd-1.0)*2.0);
  dphi1=minmod(dphi1,dphi2);
  dphi2=minmod(_phi(np[lA],gl)-_phi(np[lB],gl),_phi(np[lB],gl)-_phi(np[lC],gl));
  dphi=minmod(dphi1,dphi2);
*/


  if (dphi<Ethresholdbdry*dwall*fact_cathode) CATHODE=TRUE;


  return(CATHODE);
}



/* conditioned conductivity (must not be zero because some terms are divided by sigma_positive) */ 
double _sigma_positive(np_t *np, gl_t *gl, long l){
  double sigma_positive,sigma;
  sigma=_sigma(np,gl,l);
  sigma_positive=max(sigma,gl->model.fluid.sigmadiv);
  return(sigma_positive);
}



double _alpha(np_t *np, gl_t *gl, long l, long k, long r){
  double alpha,sigma_positive;
  sigma_positive=_sigma_positive(np,gl,l);
  alpha=_m(k)/_m(r)*(_delta(r,k)+_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_Nk(np[l],gl,k)/sigma_positive);
  return(alpha);
}


void find_alpha(np_t *np, gl_t *gl, long l, spec2_t alpha){
  long k,r;
  double sigma;
  sigma=_sigma_positive(np,gl,l);
  for (k=0; k<ns; k++){
    for (r=0; r<ns; r++){
      alpha[k][r]=_calM(k)/_calM(r)*(_delta(r,k)+_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_Nk(np[l],gl,k)/sigma);
      
    }
#ifndef NDEBUG
    if (alpha[k][k]<=0.0){
      wfprintf(stderr,"alpha[%ld][%ld]=%E   \n",k,k,alpha[k][k]);
      wfprintf(stderr,"mu_k=%E N_k=%E betaa_k=%E C_k=%E \n",_mu(np,gl,l,k),_Nk(np[l],gl,k),_betaa(gl,k),_C(k));
      wfprintf(stderr," sigma=%E betaa_k*C_k*mu_k*N_k=%E\n",sigma,_betaa(gl,k)*_C(k)*_mu(np,gl,l,k)*_Nk(np[l],gl,k));
      fatal_error("alpha[k][k] is negative at node i=%ld j=%ld k=%ld.",_i(l,gl,0),_i(l,gl,1),_i(l,gl,2));
    }
#endif
  }
  
}


void find_Z(np_t *np, gl_t *gl, long l, sqmat_t Z) {
  long row,col;
  spec2_t alpha;
  find_alpha(np, gl, l, alpha);
  set_matrix_to_identity(Z);
  for (row=0; row<ns; row++){
    for (col=0; col<ns; col++){
      Z[row][col]=alpha[row][col];
    }
  }

}


void find_dZU_dU(np_t *np, gl_t *gl, long l, sqmat_t dZU_dU){
  double sigma,dsigma_drhom;
  long k,r,m;
  sigma=_sigma_positive(np,gl,l);

  set_matrix_to_identity(dZU_dU);


  for (k=0; k<ns; k++){
    for (r=0; r<ns; r++){
      /* wrt rhok[k] */
      dZU_dU[k][k]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)/sigma*_rhok(np[l],r);
      /* wrt rhok[r] */
      dZU_dU[k][r]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_rhok(np[l],k)/sigma;
      /* wrt sigma */
      for (m=0; m<ns; m++){
        dsigma_drhom=_C(m)*_mu(np,gl,l,m)/_m(m);
        dZU_dU[k][m]+=-1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_rhok(np[l],k)/sqr(sigma)*_rhok(np[l],r)*dsigma_drhom; 
      }
    }
  }

}


void find_LambdaZ(np_t *np, gl_t *gl, long l, sqmat_t LambdaZ) {
  long r,k;
  double sigma;
  sigma=_sigma_positive(np,gl,l);
  set_matrix_to_identity(LambdaZ);
  for (k=0; k<ns; k++){
    for (r=0; r<ns; r++){
      LambdaZ[k][k]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)/sigma*_rhok(np[l],r);
    }
  }

}



void find_Y1star_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, flux_t Ystar){
  double Estar;
  long spec, flux;
  metrics_t metrics;
  find_metrics_at_interface(np,gl,lL,lR,theta,&metrics);
  find_Estar_at_interface(np, gl, lL, lR, theta, &Estar);

  for (flux=0; flux<nf; flux++){
    Ystar[flux]=0.0;
  }

  for (spec=0; spec<ns; spec++){
    Ystar[spec]=metrics.Omega*Estar*_betag(gl,spec);  
  }

}



void find_Y2star_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, flux_t Ystar){
  double Jstar;
  long spec, flux;
  metrics_t metrics;
  find_metrics_at_interface(np,gl,lL,lR,theta,&metrics);
  find_Jstar_at_interface(np, gl, lL, lR, theta, &Jstar);

  for (flux=0; flux<nf; flux++){
    Ystar[flux]=0.0;
  }

  for (spec=0; spec<ns; spec++){
    Ystar[spec]=-metrics.Omega*Jstar*_betaa(gl,spec);  
  }

}

#endif


#ifndef _FLUID_PLASMA

void find_Z(np_t *np, gl_t *gl, long l, sqmat_t Z) {
  set_matrix_to_identity(Z);
}


void find_dZU_dU(np_t *np, gl_t *gl, long l, sqmat_t dZU_dU){
  set_matrix_to_identity(dZU_dU);
}


void find_LambdaZ(np_t *np, gl_t *gl, long l, sqmat_t LambdaZ) {
  set_matrix_to_identity(LambdaZ);
}


#endif



void find_init_mass_fraction_templates(char **specstr1, char **specstr2){
  char *specname;
  *specstr1=(char *)realloc(*specstr1,1000*sizeof(char));
  *specstr2=(char *)realloc(*specstr2,1000*sizeof(char));
  #if (defined(specN2) && defined(specO2))
    strcpy(*specstr1,
    "    Species(\"O2\", \"N2\", \"default\");\n"
    "    w_O2=0.235;\n"
    "    w_N2=0.765;\n"
    "    w_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",w_O2,w_N2,w_default");
  #endif
  #if (defined(specN2) && !defined(specO2))
    strcpy(*specstr1,
    "    Species(\"N2\", \"default\");\n"
    "    w_N2=1.0;\n"
    "    w_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",w_N2,w_default");
  #endif
  #if (defined(specO2) && !defined(specN2))
    strcpy(*specstr1,
    "    Species(\"O2\",  \"default\");\n"
    "    w_O2=1.0;\n"
    "    w_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",w_O2,w_default");
  #endif
  #if (defined(specCO2) && !defined(specN2) && !defined(specO2))
    strcpy(*specstr1,
    "    Species(\"CO2\", \"default\");\n"
    "    w_CO2=1.0;\n"
    "    w_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",w_CO2,w_default");
  #endif
  if (ns==1){
    specname=(char *)malloc(sizeof(char));
    find_species_variable_name(0, &specname);
    strcpy(*specstr1,"    w_");
    strcat(*specstr1,specname);
    strcat(*specstr1,"=1.0;\n");
    strcpy(*specstr2, ",w_");
    strcat(*specstr2,specname); 
    free(specname);   
  }
}


void find_init_molar_fraction_templates(char **specstr1, char **specstr2){
  char *chargedstr1,*chargedstr2,*chargedstr3;
  char *specname;
  *specstr1=(char *)realloc(*specstr1,1000*sizeof(char));
  *specstr2=(char *)realloc(*specstr2,1000*sizeof(char));
  chargedstr1=(char *)malloc(1000*sizeof(char));
  chargedstr2=(char *)malloc(1000*sizeof(char));
  chargedstr3=(char *)malloc(1000*sizeof(char));
  strcpy(chargedstr1,"");
  strcpy(chargedstr2,"");
  strcpy(chargedstr3,"");
  
  #if (defined(speceminus) && FALSE)
    long specionplus;
    char *nameionplus;
    nameionplus=(char *)malloc(1000*sizeof(char));
    specionplus=0;
    do {
      specionplus++;
    } while(_Charge_number(specionplus)<=0 && specionplus<ns);
    find_species_name(specionplus, &nameionplus);
    if (specionplus==ns) fatal_error("Problem finding a positive ion in find_init_molar_fraction_templates().");
    strcat(chargedstr3,",\"e-\",\"");
    strcat(chargedstr3,nameionplus);
    strcat(chargedstr3,"\"");
    find_species_variable_name(specionplus, &nameionplus);
    strcat(chargedstr1,
    "    chi_eminus=1e-12;\n"
    "    chi_"
    );
    strcat(chargedstr1,nameionplus);
    strcat(chargedstr1,"=1e-12;\n");
    strcat(chargedstr2,",chi_eminus,chi_");
    strcat(chargedstr2,nameionplus);
    free(nameionplus);
  #endif   
  #if (defined(speceminus))
    strcat(chargedstr3,",\"e-\"");
    strcat(chargedstr1,
    "    chi_eminus=1e-12;\n"
    );
    strcat(chargedstr2,",chi_eminus");
  #endif   
  #if (defined(specN2) && defined(specO2))
    strcpy(*specstr1,
    "    Species(\"O2\",\"N2\"");
    strcat(*specstr1,chargedstr3);
    strcat(*specstr1,",\"default\");\n"
    "    chi_O2=0.21;\n"
    "    chi_N2=0.79;\n"
    );
    strcat(*specstr1,chargedstr1);
    strcat(*specstr1,
    "    chi_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",chi_O2,chi_N2");
    strcat(*specstr2, chargedstr2);
    strcat(*specstr2, ",chi_default");
  #endif
  #if (defined(specN2) && !defined(specO2))
    strcpy(*specstr1,
    "    Species(\"N2\", \"default\");\n"
    "    chi_N2=1.0;\n"
    "    chi_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",chi_N2,chi_default");
  #endif
  #if (defined(specO2) && !defined(specN2))
    strcpy(*specstr1,
    "    Species(\"O2\",  \"default\");\n"
    "    chi_O2=1.0;\n"
    "    chi_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",chi_O2,chi_default");
  #endif
  #if (defined(specCO2) && !defined(specN2) && !defined(specO2))
    strcpy(*specstr1,
    "    Species(\"CO2\", \"default\");\n"
    "    chi_CO2=1.0;\n"
    "    chi_default=1e-30;\n"
    );  
    strcpy(*specstr2, ",chi_CO2,chi_default");
  #endif
  if (ns==1){
    specname=(char *)malloc(sizeof(char));
    find_species_variable_name(0, &specname);
    strcpy(*specstr1,"    chi_");
    strcat(*specstr1,specname);
    strcat(*specstr1,"=1.0;\n");
    strcpy(*specstr2, ",chi_");
    strcat(*specstr2,specname); 
    free(specname);   
  }

  free(chargedstr1);
  free(chargedstr2);
  free(chargedstr3);
}



void find_init_number_density_templates(char **specstr1, char **specstr2){
  char *specname;
  char *chargedstr1,*chargedstr2,*chargedstr3;
  *specstr1=(char *)realloc(*specstr1,1000*sizeof(char));
  *specstr2=(char *)realloc(*specstr2,1000*sizeof(char));
  chargedstr1=(char *)malloc(1000*sizeof(char));
  chargedstr2=(char *)malloc(1000*sizeof(char));
  chargedstr3=(char *)malloc(1000*sizeof(char));
  strcpy(chargedstr1,"");
  strcpy(chargedstr2,"");
  strcpy(chargedstr3,"");
  
  #if (defined(speceminus))
    strcat(chargedstr3,",\"e-\"");
    strcat(chargedstr1,
    "    N_eminus=1e12; {1/m3} \n"
    );
    strcat(chargedstr2,",N_eminus");
  #endif   
  #if (defined(specN2) && defined(specO2))
    strcpy(*specstr1,
    "    Species(\"O2\",\"N2\"");
    strcat(*specstr1,chargedstr3);
    strcat(*specstr1,",\"default\");\n"
    "    N_O2=0.21*1e24; {1/m3}\n"
    "    N_N2=0.79*1e24;\n"
    );
    strcat(*specstr1,chargedstr1);
    strcat(*specstr1,
    "    N_default=1e9;\n"
    );  
    strcpy(*specstr2, ",N_O2,N_N2");
    strcat(*specstr2, chargedstr2);
    strcat(*specstr2, ",N_default");
  #endif
  #if (defined(specN2) && !defined(specO2))
    strcpy(*specstr1,
    "    Species(\"N2\", \"default\");\n"
    "    N_N2=1e24; {1/m3}\n"
    "    N_default=1e9;\n"
    );  
    strcpy(*specstr2, ",N_N2,N_default");
  #endif
  #if (defined(specO2) && !defined(specN2))
    strcpy(*specstr1,
    "    Species(\"O2\",  \"default\");\n"
    "    N_O2=1e24; {1/m3}\n"
    "    N_default=1e9;\n"
    );  
    strcpy(*specstr2, ",N_O2,N_default");
  #endif
  #if (defined(specCO2) && !defined(specN2) && !defined(specO2))
    strcpy(*specstr1,
    "    Species(\"CO2\", \"default\");\n"
    "    N_CO2=1e24; {1/m3}\n"
    "    N_default=1e9;\n"
    );  
    strcpy(*specstr2, ",N_CO2,N_default");
  #endif
  #if (defined(specH2) && defined(specHe) && !defined(specN2) && !defined(specO2))
    strcpy(*specstr1,
    "    Species(\"H2\", \"He\", \"CH4\", \"default\");\n"
    "    N=1e25;         {1/m3} \n"
    "    N_H2=N*0.7994;   {1/m3}\n"
    "    N_He=N*0.187;   {1/m3}\n"
    "    N_CH4=N*0.0136; {1/m3}\n"
    "    N_default=1e9;  {1/m3}\n"
    );  
    strcpy(*specstr2, ",N_H2,N_He,N_CH4,N_default");
  #endif
  if (ns==1){
    specname=(char *)malloc(sizeof(char));
    find_species_variable_name(0, &specname);
    strcpy(*specstr1,"    N_");
    strcat(*specstr1,specname);
    strcat(*specstr1,"=1e24;\n");
    strcpy(*specstr2, ",N_");
    strcat(*specstr2,specname); 
    free(specname);   
  }
  free(chargedstr1);
  free(chargedstr2);
  free(chargedstr3);
}






#ifdef _2D

/* see page 15 in "A new parabolized Navier-Stokes code for chemically reacting flow fields"
by Dinesh Kumble Prabhu, Iowa State University
*/
void find_Saxi(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  double x1;
  dim_t Vn;
#ifdef _FLUID_PLASMA
  long spec;
  EXM_vec3D_t Vk,Ve;
#endif

  for (flux=0; flux<nf; flux++) S[flux]=0.0e0;

  if (gl->model.fluid.AXISYMMETRIC) {
    x1=np[l].bs->x[1];
    find_V(np[l],Vn);
    if (fabs(x1)<1e-15) fatal_error("No node must lie on the y=0 axis when AXISYMMETRIC is set to TRUE.");
    for (flux=0; flux<nf; flux++){
      S[flux]=(-1.0/x1)*Vn[1]*np[l].bs->U[flux];
    }
    S[fluxet]+=(-1.0/x1)*Vn[1]*_Pstar(np[l],gl);
#ifdef _FLUID_PLASMA
    for (spec=0; spec<ncs; spec++){
      find_Vk(np, gl, l, spec, Vk);
      //S[spec]=(-1.0/x1)*Vk[1]*np[l].bs->U[spec];
    }
    find_Ve_from_J(np, gl, l, Ve);
    //S[speceminus]=(-1.0/x1)*Ve[1]*np[l].bs->U[speceminus];
#endif
#if _FLUID_EENERGY
    find_Ve_from_J(np, gl, l, Ve);
    //S[fluxee]=(-1.0/x1)*Ve[1]*np[l].bs->U[fluxee];
#endif
  }
}


void find_dSaxi_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU){
  long row,col;
  long flux;
  double x1;
  dim_t Vn;
#ifdef _FLUID_PLASMA
  long spec;
  EXM_vec3D_t Vk,Ve;
#endif

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dS_dU[row][col]=0.0e0;
    }
  }

  if (gl->model.fluid.AXISYMMETRIC) {
    x1=np[l].bs->x[1];
    find_V(np[l],Vn);
    if (fabs(x1)<1e-15) fatal_error("No node must lie on the y=0 axis when AXISYMMETRIC is set to TRUE.");
    for (flux=0; flux<nf; flux++){
      dS_dU[flux][flux]=min(0.0,(-1.0/x1)*Vn[1]);
    }
#ifdef _FLUID_PLASMA
    for (spec=0; spec<ncs; spec++){
      find_Vk(np, gl, l, spec, Vk);
      //dS_dU[spec][spec]=min(0.0,(-1.0/x1)*Vk[1]);
    }
    find_Ve_from_J(np, gl, l, Ve);
    //dS_dU[speceminus][speceminus]=min(0.0,(-1.0/x1)*Ve[1]);
#endif
#if _FLUID_EENERGY
    find_Ve_from_J(np, gl, l, Ve);
    //dS_dU[fluxee][fluxee]=min(0.0,(-1.0/x1)*Ve[1]);
#endif
  }


}

#else

void find_Saxi(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  for (flux=0; flux<nf; flux++){
    S[flux]=0.0e0;
  }
}


void find_dSaxi_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dS_dU[row][col]=0.0e0;
    }
  }
}
#endif


#if defined(_FLUID_MULTISPECIES)

double _w_product_at_catalytic_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, long specR, long specP, double factprodreact){
  double wP;
  dim_t n;
  long dim,dim2;
  double VnormalR,VnormalP,coeffwA;
  
  // first find unit normal vector pointing towards fluid
  if (!find_unit_vector_normal_to_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK, n)){
    //fatal_error("couldn't find unit vector normal to boundary place in _w_product_at_catalytic_wall");
    
    return(_w(np[lB],specP)); 
  }
      
  // idea is to set the mass flux of specP coming out equal to the mass flux of specR coming in
  // first find the mass flux of spec  due to diffusion as a vector
  VnormalR=0.0;
  for (dim=0; dim<nd; dim++) {
    for (dim2=0; dim2<nd; dim2++){
      VnormalR-=n[dim2]*_X(np[lB],dim,dim2)*_nu(np[lB],gl,specR)*0.5*(_w(np[_al(gl,lB,dim,+1)],specR)-_w(np[_al(gl,lB,dim,-1)],specR));  
    }
  }

#ifdef _FLUID_PLASMA
  EXM_vec3D_t VionRB;
  if (specR<ncs){
    find_Vk(np, gl, lB, specR, VionRB);
    VnormalR=0.0;
    for (dim=0; dim<nd; dim++) VnormalR+=_rhok(np[lA],specR)*n[dim]*VionRB[dim];
  }
#endif
  
  // only consider a change to wP if the reactant species is injected in the surface  
  VnormalR=min(0.0,VnormalR);
      
  // but VnormalP must be equal to -VnormalR*factprodreact
  VnormalP=0.0;
  coeffwA=0.0;
  for (dim=0; dim<nd; dim++) {
    for (dim2=0; dim2<nd; dim2++){
      if (dim==theta){
        VnormalP-=factprodreact*thetasgn*n[dim2]*0.5*(_X(np[lB],dim,dim2)*_nu(np[lB],gl,specP)+_X(np[lA],dim,dim2)*_nu(np[lA],gl,specP))*_w(np[lB],specP);
        coeffwA-=-factprodreact*thetasgn*n[dim2]*0.5*(_X(np[lB],dim,dim2)*_nu(np[lB],gl,specP)+_X(np[lA],dim,dim2)*_nu(np[lA],gl,specP));
      } else { 
        VnormalP-=factprodreact*n[dim2]*_X(np[lB],dim,dim2)*_nu(np[lB],gl,specP)*0.5*(_w(np[_al(gl,lB,dim,+1)],specP)-_w(np[_al(gl,lB,dim,-1)],specP));   
      }
    }
  }
  // idea is that VnormalP+coeffwA*wA=-VnormalR
  // thus: wA=(-VnormalR-VnormalP)/coeffwA
  wP=max(0.0,min(1.0,(-VnormalR-VnormalP)/coeffwA));
  
  
  return(wP);
}


void update_w_at_catalytic_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, double Twall, double Tewall, long paramstart, long paramend, spec_t wwall){
  long specr,specp,numcat,cat,spec;
  double dwall,gamma,eta,kappa;
  spec_t nuk,rhok;

  assert(is_node_bdry(np[lA],TYPELEVEL_FLUID_WORK));
  dwall=_distance_between_near_bdry_node_and_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK);
  for (spec=0; spec<ns; spec++) rhok[spec]=wwall[spec]*_rho(np[lB]);
  
  find_nuk_eta_kappa(rhok,Twall,Tewall,nuk,&eta,&kappa);

  if (mod(paramend-paramstart+1,3)!=0) fatal_error("Wrong number of extra parameters to catalytic boundary condition.");
  numcat=round((double)(paramend-paramstart+1)/3.0);
  for (cat=0; cat<numcat; cat++){
    specr=round(_bdry_param(np,gl,lA,paramstart+cat*3,TYPELEVEL_FLUID_WORK))-1;
    specp=round(_bdry_param(np,gl,lA,paramstart+1+cat*3,TYPELEVEL_FLUID_WORK))-1;
    if (specr==specp) fatal_error("Wrong specification of the catalytic boundary condition. The reactant species can not be the same as the product species in a wall catalytic process.");
    if (specr<0 || specr>=ns) fatal_error("Wrong specification of the catalytic boundary condition. The reactant species number is not within bounds.");
    if (specp<0 || specp>=ns) fatal_error("Wrong specification of the catalytic boundary condition. The product species number is not within bounds.");
    gamma=_bdry_param(np,gl,lA,paramstart+2+cat*3,TYPELEVEL_FLUID_WORK);
    if (gamma==0.0) fatal_error("The catalytic recombination coefficient can not be set to zero. Set it rather to a very small positive number approaching zero.");
    if (gamma==1.0) fatal_error("The catalytic recombination coefficient can not be set to 1. Set it rather to a positive number less then and approaching 1.");
    if (gamma>1.0) fatal_error("The catalytic recombination coefficient can not be set to a value greater than 1.");
    if (gamma<0.0) fatal_error("The catalytic recombination coefficient can not be set to a value less than 0.");
    
    wwall[specr]=max(0.0,_w(np[lB],specr)-dwall*0.25*wwall[specr]*_rho(np[lB])*sqrt(8.0*calR/_calM(specr)*Twall/pi)/(nuk[specr]/gamma-nuk[specr]));
    //wwall[specp]=_w(np[lB],specp)+nuk[specr]/nuk[specp]*(_w(np[lB],specr)-wwall[specr]);
    wwall[specp]=_w_product_at_catalytic_wall(np, gl, lA, lB, lC, theta, thetasgn, specr, specp, 1.0);
    //printf("%E ",wwall[specr]);
  }
}


void update_w_V_at_injection_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, double Twall, double Tewall, 
                                  long paramstart, long paramend, spec_t wwall, dim_t Vwall){
  long dim,spec,numinj,inj;
  double sum,mdot,mdottot,dwall,eta,kappa;
  spec_t nuk,nukA,nukB,rhok;
  dim_t n;
  assert(is_node_bdry(np[lA],TYPELEVEL_FLUID_WORK));

  if (!find_unit_vector_normal_to_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK, n)){
    fatal_error("Problem finding unit vector normal to boundary plane in update_w_V_at_injection_wall().");
  }
  dwall=_distance_between_near_bdry_node_and_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK);

  for (spec=0; spec<ns; spec++) rhok[spec]=wwall[spec]*_rho(np[lB]);
  find_nuk_eta_kappa(rhok,Twall,Tewall,nukB,&eta,&kappa);
  for (spec=0; spec<ns; spec++) rhok[spec]=wwall[spec]*_rho(np[lA]);
  find_nuk_eta_kappa(rhok,Twall,Tewall,nukA,&eta,&kappa);
  for (spec=0; spec<ns; spec++) nuk[spec]=0.5*(nukB[spec]+nukA[spec]);
  
  //fatal_error("\n param=%ld",np[lA].numbdryparam);
  
  if (mod(paramend-paramstart+1,2)!=0) fatal_error("Wrong number of extra parameters to injection boundary condition.");
  numinj=round((double)(paramend-paramstart)/2.0);
  // first find the total mdot
  mdottot=0;
  for (inj=0; inj<numinj; inj++){
    mdottot+=_bdry_param(np,gl,lA,paramstart+1+inj*2,TYPELEVEL_FLUID_WORK);
  }  
  
  

  // second set the mass fractions of the non-injected species
  for (spec=0; spec<ns; spec++){
    wwall[spec]=max(0.0,_w(np[lB],spec)-dwall*(_w(np[lB],spec)*mdottot)/nuk[spec]);
  }


  // third set the mass fractions of the injected species
  for (inj=0; inj<numinj; inj++){
    spec=round(_bdry_param(np,gl,lA,paramstart+inj*2,TYPELEVEL_FLUID_WORK))-1;
    if (spec<0 || spec>=ns) fatal_error("Wrong specification of the injection boundary condition. The injection species number is not within bounds.");
    mdot=_bdry_param(np,gl,lA,paramstart+1+inj*2,TYPELEVEL_FLUID_WORK);
    if (mdot<=0.0) fatal_error("The injection mass flow rate of species %ld can not be set to a value less than or equal to 0.",spec);
    wwall[spec]=max(0.0,_w(np[lB],spec)-dwall*(_w(np[lB],spec)*mdottot-mdot)/nuk[spec]);
    //printf("%E %E\n",wwall[spec],_w(np[lB],spec));
  }
  // fourth make sure the mass fractions sum to 1
  sum=0.0;
  for (spec=0; spec<ns; spec++) sum+=wwall[spec];
  for (spec=0; spec<ns; spec++) wwall[spec]/=sum;
  
  // fifth update Vwall
  for (dim=0; dim<nd; dim++) Vwall[dim]=n[dim]*mdottot/_rho(np[lB]);
}
#endif
