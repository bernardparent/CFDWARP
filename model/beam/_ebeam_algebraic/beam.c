// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2005-2006 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#error The algebraic ebeam module can not be compiled with a fluid module without plasma support.
#endif


void write_model_beam_template(FILE **controlfile){
  wfprintf(*controlfile,
    " \n"
    "  EbeamAlgebraic(\n"
    "    SetQbeamoverN(is, "
                     if2D("js, ")
                     if3D("ks, ")
                          "ie, "
                     if2D("je, ")
                     if3D("ke, ")
    "  {QbeamoverN, W}0.0);\n"
    "  );\n"
  );
}



void read_and_init_beam_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
  long i,j,k;

  if (strcmp(actionname,"SetQbeamoverN")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=nd*2+1)
      SOAP_fatal_error(codex,"Number of arguments not equal to %ld in SetQbeamoverN(); action.",nd*2+1);

/*
*((readcontrolarg_t *)codex->action_args)->np
*/
    for1DL(i,SOAP_get_argum_long(codex,*argum,0), SOAP_get_argum_long(codex,*argum,nd))
      for2DL(j,SOAP_get_argum_long(codex,*argum,1), SOAP_get_argum_long(codex,*argum,nd+1))
        for3DL(k,SOAP_get_argum_long(codex,*argum,2), SOAP_get_argum_long(codex,*argum,nd+2))
          (*((readcontrolarg_t *)codex->action_args)->np)[_ai(((readcontrolarg_t *)codex->action_args)->gl,i,j,k)].bs->QbeamoverN=SOAP_get_argum_double(codex,*argum,2*nd);
        end3DL
      end2DL
    end1DL
    codex->ACTIONPROCESSED=TRUE;
  }
}


void read_model_beam_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long i,j,k;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_BEAM_ACTIONNAME)==0) {
    wfprintf(stdout,"%s..",_BEAM_ACTIONNAME);
    gl->MODEL_BEAM_READ=TRUE;  
  
    for1DL(i,gl->domain_lim.is,gl->domain_lim.ie)
      for2DL(j,gl->domain_lim.js,gl->domain_lim.je)
        for3DL(k,gl->domain_lim.ks,gl->domain_lim.ke)
	      (*((readcontrolarg_t *)codex->action_args)->np)[_ai(((readcontrolarg_t *)codex->action_args)->gl,i,j,k)].bs->QbeamoverN=0.0;
	end3DL
      end2DL
    end1DL

    action_original=codex->action;
    codex->action=&read_and_init_beam_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;
    codex->ACTIONPROCESSED=TRUE;
  }
}




double _Qbeam(np_t np, gl_t *gl){
  double Qbeam;
  Qbeam=np.bs->QbeamoverN*_N(np,gl);
  return(Qbeam);
}


double _dQbeam_dN(np_t np, gl_t *gl){
  double dQbeamdN;
  dQbeamdN=np.bs->QbeamoverN;
  return(dQbeamdN);
}

