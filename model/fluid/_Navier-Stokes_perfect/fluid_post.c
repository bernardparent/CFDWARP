// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2010-2011 Bernard Parent

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
#include "fluid_post.h"
#include "fluid_bdry.h"
#include "fluid_source.h"
#include <model/metrics/_metrics.h>
#include <model/_model.h>
#include <model/share/fluid_share.h>


typedef struct {
  double R,gamma,h1,q1;
} argfpot_t;

typedef struct {
  double R,gamma,h1,q1;
} argpstag_t;

static double _qcnew_funct(void *arg, double q, double P){
  double tmp,cp,h,T;
  h=-sqr(q)/2.0+((argfpot_t *)arg)->h1+sqr(((argfpot_t *)arg)->q1)/2.0;
  cp=((argfpot_t *)arg)->R*((argfpot_t *)arg)->gamma/(((argfpot_t *)arg)->gamma-1.0);
  T=h/cp;
  if (T>60000.0) {
    wfprintf(stderr,"T=%E along numerical differentiation of qc\n",T);
  }
  if (!isfinite(T))
    fatal_error("T not finite in functqc.");
  tmp=-(((argfpot_t *)arg)->R*T)/(P*q);
  return(tmp);
}




/* the thermally perfect and calorically perfect stagnation pressure */
double _Pstag(np_t np, gl_t *gl, long numsteps){
  argpstag_t arg;
  double Pstag,cp,M;

  arg.R=gl->model.fluid.R;
  arg.gamma=gl->model.fluid.gamma;
  arg.q1=_q(np);
  cp=gl->model.fluid.R*gl->model.fluid.gamma/(gl->model.fluid.gamma-1.0);
  arg.h1=cp*_T(np,gl);
  M=_q(np)/_a(np,gl);
  Pstag=_P(np,gl)*pow(1.0+(arg.gamma-1.0)/2.0*sqr(M), arg.gamma/(arg.gamma-1.0));
  return(Pstag);
}


int reversibly_expand_q_and_rho_to_given_P(np_t np, gl_t *gl, double Pc, double *qc, double *rhoc,
           long numsteps, double q_min){
  argfpot_t arg;
  double Tc,h,cp;
  long error;
  int retval;

  cp=gl->model.fluid.R*gl->model.fluid.gamma/(gl->model.fluid.gamma-1.0);

  arg.R=gl->model.fluid.R;
  arg.gamma=gl->model.fluid.gamma;
  arg.q1=_q(np);
  arg.h1=cp*_T(np,gl);
  /* if the desired back pressure is more than the flow stagnation pressure,
     it is not possible to integrate; for robustness, first check
     if the stagnation pressure of the flow is greater than the back
     pressure. If not, then neglect this streamline. */
  if (_Pstag(np,gl,numsteps)>Pc) {
    *qc=EXM_numerical_differentiation(&_qcnew_funct, &arg, EXM_NUMDIFF_MODIFIEDEULER,
                  numsteps, max(q_min,_q(np)), _P(np,gl), Pc, &error);
    h=-sqr(*qc)/2.0+arg.h1+sqr(arg.q1)/2.0;
    Tc=h/cp;
    if (Tc<1.0 || Tc>60000.0) wfprintf(stderr,"big problem here..Tc=%E\n",Tc);
    *rhoc=Pc/(arg.R*Tc);
    retval=0;
  } else {
    retval=1;
  }
  return(retval);
}





double _Tstag(np_t np, gl_t *gl){
  double cp,hstag,Tstag;
  cp=gl->model.fluid.R*gl->model.fluid.gamma/(gl->model.fluid.gamma-1.0);
  hstag=np.bs->U[1+nd]/_rho(np)+_P(np,gl)/_rho(np);
  Tstag=hstag/cp;
  return(Tstag);
}

double _Pstar(np_t np, gl_t *gl){
  double Pstar;
  Pstar=_P(np,gl);
  return(Pstar);
}

double _htstar(np_t np, gl_t *gl){
  double htstar;
  htstar=_ht(np,gl);
  return(htstar);
}


void find_post_variable_name_fluid(long varnum, char *varname){

  if (varnum<nd) sprintf(varname,"V[%ld]",varnum);
  if (varnum>=nd && varnum<nd*2) sprintf(varname,"M[%ld]",varnum-nd);

  switch (varnum) {
    case 2*nd+0:   sprintf(varname,"rho");       break;
    case 2*nd+1:   sprintf(varname,"P");         break;
    case 2*nd+2:   sprintf(varname,"T");         break;
    case 2*nd+3:   sprintf(varname,"a");         break;
  }
}



void find_post_variable_value_fluid(np_t *np, long l, gl_t *gl,
                          long varnum, double *varvalue){
  double a;

  *varvalue=0.0;
  if (is_node_valid(np[l],TYPELEVEL_FLUID)){
    a=_a(np[l],gl);
    if (varnum<nd) *varvalue=_V(np[l],varnum);
    if (varnum<2*nd && varnum>=nd) *varvalue=_V(np[l],varnum-nd)/a;
    
    
    switch (varnum) {
      case 2*nd+0:   *varvalue=_rho(np[l]);       break;
      case 2*nd+1:   *varvalue=_P(np[l],gl);       break;
      case 2*nd+2:   *varvalue=_T(np[l],gl);       break;
      case 2*nd+3:   *varvalue=_a(np[l],gl);       break;
    }
  }
}

