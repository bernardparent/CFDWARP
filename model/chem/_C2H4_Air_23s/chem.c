// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2023 Bernard Parent

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

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/share/chem_share.h>
#include <src/control.h>
#include "zettervall2017.h"
#include "zettervall2017b.h"
#include "konnovzettervall.h"

#define CHEMMODEL_NONE 1
#define CHEMMODEL_ZETTERVALL2017 2
#define CHEMMODEL_ZETTERVALL2017B 3
#define CHEMMODEL_KONNOVZETTERVALL 4



void write_model_chem_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    CHEMMODEL=CHEMMODEL_ZETTERVALL2017;\n"
    "    TMIN_LINDEMANN=500.0;\n"
    "  );\n"
  ,_CHEM_ACTIONNAME);
}



void read_model_chem_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_CHEM_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);

    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_CHEM_ACTIONNAME);
    SOAP_add_int_to_vars(codex,"CHEMMODEL_NONE",CHEMMODEL_NONE); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_ZETTERVALL2017",CHEMMODEL_ZETTERVALL2017); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_ZETTERVALL2017B",CHEMMODEL_ZETTERVALL2017B); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_KONNOVZETTERVALL",CHEMMODEL_KONNOVZETTERVALL); 
    gl->MODEL_CHEM_READ=TRUE;

    action_original=codex->action;
    codex->action=&read_model_chem_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"CHEMMODEL",&gl->model.chem.CHEMMODEL);
    if (gl->model.chem.CHEMMODEL!=CHEMMODEL_ZETTERVALL2017 && gl->model.chem.CHEMMODEL!=CHEMMODEL_ZETTERVALL2017B 
     && gl->model.chem.CHEMMODEL!=CHEMMODEL_KONNOVZETTERVALL && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)
      SOAP_fatal_error(codex,"CHEMMODEL must be set to either CHEMMODEL_KONNOVZETTERVALL or CHEMMODEL_ZETTERVALL2017 or CHEMMODEL_NONE or CHEMMODEL_ZETTERVALL2017B.");

    find_double_var_from_codex(codex,"TMIN_LINDEMANN",&gl->model.chem.TMIN_LINDEMANN);

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}






void find_W_None ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;

  for ( k = 0; k < ns; k++ ) {
    W[k] = 0.0;
  }
}


void find_dW_dx_None ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */

  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }
}


void find_W ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  switch (gl->model.chem.CHEMMODEL){
    case CHEMMODEL_ZETTERVALL2017: 
      find_W_Zettervall2017 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_ZETTERVALL2017B: 
      find_W_Zettervall2017b ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_KONNOVZETTERVALL: 
      find_W_KonnovZettervall ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_NONE: 
      find_W_None ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );    
    break;
    default:
      fatal_error("Problem with CHEMMODEL in find_W() within chem.c");
  }
}


void find_dW_dx ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  switch (gl->model.chem.CHEMMODEL){
    case CHEMMODEL_ZETTERVALL2017: 
      find_dW_dx_Zettervall2017 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_ZETTERVALL2017B: 
      find_dW_dx_Zettervall2017b ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_KONNOVZETTERVALL: 
      find_dW_dx_KonnovZettervall ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_NONE: 
      find_dW_dx_None ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    default:
      fatal_error("Problem with CHEMMODEL in find_W() within chem.c");
  }

}





void find_Qei(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  *Qei=0.0;  
}


void find_dQei_dx(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
  
}


