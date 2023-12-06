// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2018,2021 Bernard Parent
Copyright 2021 Prasanna Thoguluva Rajendran

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
#include "macheret2007old.h"
#include "macheret2007.h"
#include "rajendran2022.h"

#define CHEMMODEL_NONE 1
#define CHEMMODEL_MACHERET2007 2
#define CHEMMODEL_MACHERET2007OLD 3
#define CHEMMODEL_RAJENDRAN2022 4


#define Estarmin 1e-40


void write_model_chem_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    CHEMMODEL=CHEMMODEL_MACHERET2007;\n"
    "    QEISOURCETERMS=TRUE; {include electron energy cooling due to electron impact}\n"
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
    SOAP_add_int_to_vars(codex,"CHEMMODEL_MACHERET2007",CHEMMODEL_MACHERET2007); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_MACHERET2007OLD",CHEMMODEL_MACHERET2007OLD); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_RAJENDRAN2022",CHEMMODEL_RAJENDRAN2022); 
    gl->MODEL_CHEM_READ=TRUE;

    action_original=codex->action;
    codex->action=&read_model_chem_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"CHEMMODEL",&gl->model.chem.CHEMMODEL);
    if (gl->model.chem.CHEMMODEL!=CHEMMODEL_MACHERET2007 && gl->model.chem.CHEMMODEL!=CHEMMODEL_MACHERET2007OLD && gl->model.chem.CHEMMODEL!=CHEMMODEL_RAJENDRAN2022 
        && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)
      SOAP_fatal_error(codex,"CHEMMODEL must be set to either CHEMMODEL_MACHERET2007 or CHEMMODEL_MACHERET2007OLD or CHEMMODEL_NONE or CHEMMODEL_RAJENDRAN2022.");
    find_bool_var_from_codex(codex,"QEISOURCETERMS",&gl->model.chem.QEISOURCETERMS);

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
    case CHEMMODEL_MACHERET2007: 
      find_W_Macheret2007 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_MACHERET2007OLD: 
      find_W_Macheret2007old ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_RAJENDRAN2022: 
      find_W_Rajendran2022 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
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
    case CHEMMODEL_MACHERET2007: 
      find_dW_dx_Macheret2007 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_MACHERET2007OLD: 
      find_dW_dx_Macheret2007old ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_RAJENDRAN2022: 
      find_dW_dx_Rajendran2022 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_NONE: 
      find_dW_dx_None ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    default:
      fatal_error("Problem with CHEMMODEL in find_W() within chem.c");
  }
}


void find_Qei(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  
  *Qei=0.0;  
  if (gl->model.chem.QEISOURCETERMS){
    switch (gl->model.chem.CHEMMODEL){
      case CHEMMODEL_MACHERET2007: 
        find_Qei_Macheret2007 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_MACHERET2007OLD: 
        find_Qei_Macheret2007old ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_RAJENDRAN2022: 
        find_Qei_Rajendran2022 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_NONE: 
        *Qei=0.0;
      break;
      default:
        fatal_error("Problem with CHEMMODEL in find_Qei() within chem.c");
    }
  }
}


void find_dQei_dx(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
  
  if (gl->model.chem.QEISOURCETERMS){
    switch (gl->model.chem.CHEMMODEL){
      case CHEMMODEL_MACHERET2007: 
        find_dQei_dx_Macheret2007 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_MACHERET2007OLD: 
        find_dQei_dx_Macheret2007old ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_RAJENDRAN2022: 
        find_dQei_dx_Rajendran2022 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_NONE: 
        *dQeidTe=0.0;
      break;
      default:
        fatal_error("Problem with CHEMMODEL in find_Qei() within chem.c");
    }
  }
}
