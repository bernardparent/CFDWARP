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

#include <src/common.h>
#include <model/emfield/_emfield.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/chem/_chem.h>
#include <model/metrics/_metrics.h>
#include <src/control.h>
#include <src/bdry.h>

#ifndef _FLUID_PLASMA
#error The fixed ebeam module can not be compiled with a fluid module without plasma support.
#endif


void write_model_beam_template(FILE **controlfile){
  wfprintf(*controlfile,
    " \n"
    "  EbeamFixed(\n"
    "    SetQbeam(is,"
                     if2DL("js,")
                     if3DL("ks,")
                          " ie,"
                     if2DL("je,")
                     if3DL("ke,")
    "  0.0{W/m3});\n"
    "  );\n"
  );
}



void read_and_init_beam_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
  long i,j,k;
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (strcmp(actionname,"SetQbeam")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=nd*2+1)
      SOAP_fatal_error(codex,"Number of arguments not equal to %ld in SetQbeam(); action.",nd*2+1);

/*
*((readcontrolarg_t *)codex->action_args)->np
*/
    for_1DL(i,SOAP_get_argum_long(codex,*argum,0), SOAP_get_argum_long(codex,*argum,nd)){
      for_2DL(j,SOAP_get_argum_long(codex,*argum,1), SOAP_get_argum_long(codex,*argum,nd+1)){
        for_3DL(k,SOAP_get_argum_long(codex,*argum,2), SOAP_get_argum_long(codex,*argum,nd+2)){
          if (is_node_in_zone(i, j, k, gl->domain_lim)){
            (*((readcontrolarg_t *)codex->action_args)->np)[_ai(((readcontrolarg_t *)codex->action_args)->gl,i,j,k)].bs->Qbeam=SOAP_get_argum_double(codex,*argum,2*nd);
          }
        }
      }
    }
    codex->ACTIONPROCESSED=TRUE;
  }
}


void read_model_beam_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long i,j,k;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_BEAM_ACTIONNAME)==0) {
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_BEAM_ACTIONNAME);
    if (!gl->CONTROL_READ){
      gl->MODEL_BEAM_READ=TRUE;   
      for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
	      (*((readcontrolarg_t *)codex->action_args)->np)[_ai(((readcontrolarg_t *)codex->action_args)->gl,i,j,k)].bs->Qbeam=0.0;
      }
    }
    action_original=codex->action;
    codex->action=&read_and_init_beam_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;
    codex->ACTIONPROCESSED=TRUE;
  }
}




double _Qbeam(np_t np, gl_t *gl){
  double Qbeam;
    Qbeam=np.bs->Qbeam;
  return(Qbeam);
}


double _dQbeam_dN(np_t np, gl_t *gl){
  return(0.0);
}
