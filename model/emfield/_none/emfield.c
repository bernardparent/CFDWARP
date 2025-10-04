// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2005-2018 Bernard Parent

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

#include <src/common.h>
#include <model/emfield/_emfield.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/chem/_chem.h>
#include <model/metrics/_metrics.h>
#include <model/fluid/_fluid.h>
#include <cycle/_cycle.h>
#include <src/control.h>
#include <src/bdry.h>


void write_model_emfield_template(FILE **controlfile){

}

void write_cycle_emfield_template(FILE **controlfile){

}

void write_disc_emfield_template(FILE **controlfile){

}


void read_disc_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->DISC_EMFIELD_READ=TRUE;  
}


void read_model_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->MODEL_EMFIELD_READ=TRUE;  
}


void read_cycle_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->CYCLE_EMFIELD_READ=TRUE;  
}


double _phi(np_t np, gl_t *gl){
  return(0.0);
}



void find_emfield_force(np_t *np, gl_t *gl, long l, dim_t Femfield){
  long dim;
  for (dim=0; dim<nd; dim++) Femfield[dim]=0.0;
}


double _E_dot_J(np_t *np, gl_t *gl, long l){
  double EdotJ;
  EdotJ=0.0;
  return(EdotJ);
}

double _E_dot_Je(np_t *np, gl_t *gl, long l){
  double EdotJe;
  EdotJe=0.0;
  return(EdotJe);
}



void find_sigma(np_t *np, gl_t *gl, long l, double *sigma){
  *sigma=0.0;
}


void find_B(np_t np, gl_t *gl, EXM_vec3D_t B){
  long dim;
  for (dim=0; dim<3; dim++) B[dim]=0.0;
}

void find_E(np_t np, gl_t *gl, EXM_vec3D_t E){
  long dim;
  for (dim=0; dim<3; dim++) E[dim]=0.0;
}


void find_Ek(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Ek){
  long dim;
  for (dim=0; dim<3; dim++) Ek[dim]=0.0;
}


void find_E_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, EXM_vec3D_t E){
  long n;
  for (n=0; n<3; n++) E[n]=0.0e0;
}

double _mu(np_t *np, gl_t *gl, long l, long spec){
  double muk;
#ifdef _FLUID_PLASMA
  double N;
  N=_N(np[l],gl);
  muk=_muk_from_N_Tk_Ek(N, _Tk(np,gl,l,spec), 0.0, spec);
#else
  muk=0.0;
#endif
  return(muk);
}

double _sigma(np_t *np, gl_t *gl, long l){
  double sigma;
#ifdef _FLUID_PLASMA
  long spec;
  sigma=0.0;
  for (spec=0; spec<ns; spec++){
    if (_Charge_number(spec)!=0) {
      sigma+=fabs(_C(spec))*max(0.0,_Nk(np[l],gl,spec))*_mu(np,gl,l,spec);
    }
  }
#else
  sigma=0.0;
#endif
  return(sigma);
}



void find_Ee_for_Townsend_ionization(np_t *np, gl_t *gl, long l, EXM_vec3D_t E){
  long i;
  for (i=0; i<3; i++) E[i]=0.0;
}

void find_Jstar_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, double *Jstar){
  *Jstar=0.0;
}


void find_Estar_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, double *Estar){
  *Estar=0.0;
}


void find_DeltaVk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t DeltaVk){
  long i;
  for (i=0; i<3; i++) DeltaVk[i]=0.0;  
}

void find_Vk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Vk){
  dim_t Vn;
  long dim;
  find_V(np[l],Vn);    
  for (dim=0; dim<3; dim++) Vk[dim]=0.0;
  for (dim=0; dim<nd; dim++) Vk[dim]=Vn[dim];
}
