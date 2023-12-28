// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2005-2006 Bernard Parent

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
#include <model/thermo/_thermo.h>
#include <model/transport/_transport.h>
#include <model/metrics/_metrics.h>
#include <model/emfield/_emfield.h>
#include <model/_model.h>
#include <model/share/fluid_share.h>


typedef struct {
  double R,h1,q1,k1;
  spec_t w1;
} argfpot_t;

typedef struct {
  double R,h1,q1,k1;
  spec_t w1;
} argpstag_t;

static double _qcnew_funct(void *arg, double q, double Pstar){
  double ret,h,T;
  h=-sqr(q)/2.0+((argfpot_t *)arg)->h1+sqr(((argfpot_t *)arg)->q1)/2.0;
  T=_T_from_w_h(((argfpot_t *)arg)->w1, h);
  if (T>60000.0) {
    wfprintf(stderr,"T=%E along numerical differentiation of qc\n",T);
  }
  if (!isfinite(T))
    fatal_error("T not finite in functqc.");
  ret=-(((argfpot_t *)arg)->R*T+2.0/3.0*((argfpot_t *)arg)->k1)/(Pstar*q);
  return(ret);
}


static double _Pstag_funct(void *arg, double q){
  double ret,h,T;
  h=(-sqr(q)/2.0+((argpstag_t *)arg)->h1+sqr(((argpstag_t *)arg)->q1)/2.0);
  T=_T_from_w_h(((argpstag_t *)arg)->w1, h);
  if (T>6000.0)
    wfprintf(stderr,"T=%E (above high T polynomials limit) along integration path of _Pstag\n",T);
  ret=q/(((argpstag_t *)arg)->R*T+2.0/3.0*((argpstag_t *)arg)->k1);
  return(ret);
}



/* the thermally perfect, but calorically imperfect stagnation pressure,
   the way it should be calculated. Note that this becomes the well known
   ideal gas isentropic formula in the case of a calorically
   perfect gas; numsteps is the number of numerical integration intervals taken*/
double _Pstag(np_t np, gl_t *gl, long numsteps){
  argpstag_t arg;
  long error;
  double Pstag;

  find_w(np,arg.w1);
  arg.R=_R(arg.w1);
  arg.q1=_q(np);
  arg.k1=_k(np);
  arg.h1=_h_from_w_T(arg.w1,_T(np,gl));
  Pstag=_Pstar(np,gl)*exp(EXM_numerical_integration(&_Pstag_funct, &arg, EXM_NUMINTEG_POLY2,
                       numsteps, 0.0e0, _q(np), &error));
  return(Pstag);
}


int reversibly_expand_q_and_rho_to_given_P(np_t np, gl_t *gl, double Pstarc, double *qc, double *rhoc,
           long numsteps, double q_min){
  argfpot_t arg;
  double Tc;
  long error;
  int retval;

  find_w(np,arg.w1);
  arg.R=_R(arg.w1);
  arg.q1=_q(np);
  arg.k1=_k(np);
  arg.h1=_h_from_w_T(arg.w1,_T(np,gl));
  /* if the desired back pressure is more than the flow stagnation pressure,
     it is not possible to integrate; for robustness, first check
     if the stagnation pressure of the flow is greater than the back
     pressure. If not, then neglect this streamline. */
  if (_Pstag(np,gl,numsteps)>Pstarc) {
    *qc=EXM_numerical_differentiation(&_qcnew_funct, &arg, EXM_NUMDIFF_MODIFIEDEULER,
                  numsteps, max(q_min,_q(np)), _Pstar(np,gl), Pstarc, &error);
    Tc=_T_from_w_h(arg.w1, -sqr(*qc)/2.0+arg.h1+sqr(arg.q1)/2.0);
    if (Tc<1.0 || Tc>60000.0) wfprintf(stderr,"big problem here..Tc=%E\n",Tc);
    *rhoc=Pstarc/(arg.R*Tc+2.0/3.0*arg.k1);
    retval=0;
  } else {
    retval=1;
  }
  return(retval);
}




double _Tstag(np_t np, gl_t *gl){
  double hstag,Tstag;
  spec_t w;
  find_w(np,w);
  hstag=np.bs->U[fluxet]/_rho(np)-_k(np)-_w(np,specN2)*_ev(np)+_P(np,gl)/_rho(np);
  Tstag=_T_from_w_h(w, hstag);
  return(Tstag);
}


void find_post_variable_name_fluid(long varnum, char *varname){
  char *speciesname;
  speciesname=(char *)malloc(sizeof(char));
  if (varnum<ns) {
    find_species_name(varnum,&speciesname);
    sprintf(varname,"chi_%s",speciesname);
  }
  free(speciesname);

  if (varnum>=ns && varnum<ns+nd) sprintf(varname,"V[%ld]",varnum-ns);

  switch (varnum) {
    case nd+ns+0:   sprintf(varname,"rho");       break;
    case nd+ns+1:   sprintf(varname,"P");         break;
    case nd+ns+2:   sprintf(varname,"Pstar");     break;
    case nd+ns+3:   sprintf(varname,"T");         break;
    case nd+ns+4:   sprintf(varname,"Tv");         break;
    case nd+ns+5:   sprintf(varname,"a");         break;
    case nd+ns+6:   sprintf(varname,"eta");        break;
    case nd+ns+7:   sprintf(varname,"etastar");    break;
    case nd+ns+8:   sprintf(varname,"kappastar"); break;
    case nd+ns+9:   sprintf(varname,"k");         break;
    case nd+ns+10:  sprintf(varname,"epsilon");       break;
    case nd+ns+11:  sprintf(varname,"omega");     break;
    case nd+ns+12:  sprintf(varname,"gamma");     break;
    case nd+ns+13:  sprintf(varname,"N");     break;
  }
}



void find_post_variable_value_fluid(np_t *np, long l, gl_t *gl,
                          long varnum, double *varvalue){
  spec_t rhok,nu;
  double eta,kappa,N;
  long spec;

  *varvalue=0.0;
  if (is_node_valid(np[l],TYPELEVEL_FLUID)){
    N=0.0;
    for (spec=0; spec<ns; spec++) N+=_rhok(np[l],spec)/_m(spec);
    if (varnum<ns) *varvalue=_rhok(np[l],varnum)/_m(varnum)/N;
    if (varnum>=ns && varnum<ns+nd) *varvalue=_V(np[l],varnum-ns);
    //assert(is_node_resumed(np[l]));
    find_rhok(np[l],rhok);
    find_nuk_eta_kappa(gl, rhok, _T(np[l],gl), _Te(np[l],gl), nu, &eta, &kappa);
    
    switch (varnum) {
      case nd+ns+0:   *varvalue=_rho(np[l]);       break;
      case nd+ns+1:   *varvalue=_P(np[l],gl);         break;
      case nd+ns+2:   *varvalue=_Pstar(np[l],gl);     break;
      case nd+ns+3:   *varvalue=_T(np[l],gl);         break;
      case nd+ns+4:   *varvalue=_Tv(np[l]);         break;
      case nd+ns+5:   *varvalue=_a(np[l],gl);         break;
      case nd+ns+6:   *varvalue=_eta(np,l,gl);        break;
      case nd+ns+7:   *varvalue=_etastar(np,l,gl);    break;
      case nd+ns+8:   *varvalue=_kappastar(np,l,gl); break;
      case nd+ns+9:   *varvalue=_k(np[l]);         break;
      case nd+ns+10:  *varvalue=_eps(np[l],gl);       break;
      case nd+ns+11:  *varvalue=_omega(np,l,gl);     break;
      case nd+ns+12:  *varvalue=_gamma(np[l],gl);     break;
      case nd+ns+13:  *varvalue=N;     break;
    }
  }
}

