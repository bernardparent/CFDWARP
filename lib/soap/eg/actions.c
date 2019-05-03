// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000 Bernard Parent

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

#include <stdlib.h>
#include <stdio.h>
#include <soap.h>
#include <string.h>

typedef struct {
  double x;
} action_args_t;

void actions(char *actionname, char **argum, SOAP_codex_t *codex){
  if (strcmp(actionname,"assign_to_x")==0) {
    SOAP_substitute_all_argums(argum,codex);
    sscanf(*argum,"%lg",&(((action_args_t *)codex->action_args)->x));
  }
}

int main(){
  char code[]=""
              " assign_to_x(10.0e0);"
              "";
  SOAP_codex_t codex;
  action_args_t action_args;

  SOAP_init_codex(&codex,"?");
  codex.ACTION=TRUE;
  codex.action=&actions;
  codex.action_args=(void *)&action_args;
  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
  SOAP_free_codex(&codex);

  printf("x=%E\n",action_args.x);
  return(0);
}

