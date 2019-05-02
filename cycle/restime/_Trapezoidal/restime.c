// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015-2018 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cycle/restime/_restime.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>


void write_disc_restime_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    weightm1_convection=0.5;\n"
    "    weightm1_default=0.5;\n"
    "  );\n"
  ,_RESTIME_ACTIONNAME);
}


void read_disc_restime_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_disc_restime_actions(char *actionname, char **argum, SOAP_codex_t *codex){

  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_RESTIME_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_RESTIME_ACTIONNAME);

    gl->DISC_RESTIME_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_restime_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_double_var_from_codex(codex,"weightm1_convection",&gl->cycle.restime.weightm1_trapezoidal_convection);
    find_double_var_from_codex(codex,"weightm1_default",&gl->cycle.restime.weightm1_trapezoidal_default);

    if (gl->cycle.restime.weightm1_trapezoidal_convection<0.0 || gl->cycle.restime.weightm1_trapezoidal_convection>1.0)
      SOAP_fatal_error(codex,"In the Trapezoidal() module, weightm1_convection must lie within 0 and 1.");
    if (gl->cycle.restime.weightm1_trapezoidal_default<0.0 || gl->cycle.restime.weightm1_trapezoidal_default>1.0)
      SOAP_fatal_error(codex,"In the Trapezoidal() module, weightm1_default must lie within 0 and 1.");

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }


}


void add_Z_dUstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  long l;
  long flux;
  flux_t tmp,dRes;
  sqmat_t Z;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    for (flux=0; flux<nf; flux++){
      tmp[flux]=_Omega(np[l],gl)*(np[l].bs->U[flux]-np[l].bs->Um1[flux])/gl->dt;
    }      
    find_Z(np,gl,l,Z);
    multiply_matrix_and_vector(Z, tmp, dRes);
    for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=dRes[flux];
  }
}


