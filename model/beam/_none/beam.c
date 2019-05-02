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
#include <src/control.h>


void write_model_beam_template(FILE **controlfile){
}



void read_and_init_beam_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_beam_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->MODEL_BEAM_READ=TRUE;  
}




double _Qbeam(np_t np, gl_t *gl){
  double Qbeam;
    Qbeam=0.0e0;
  return(Qbeam);
}


double _dQbeam_dN(np_t np, gl_t *gl){
  return(0.0);
}
