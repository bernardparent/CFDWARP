// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <soap.h>
#include <string.h>

typedef struct {
  double *x;
} function_args_t;


void functions(char *functionname, char **argum, char **returnstr, SOAP_codex_t *codex){
  long cnt;

  if (strcmp(functionname,"x")==0) {
    SOAP_substitute_all_argums(argum,codex);
    sscanf(*argum,"%ld",&cnt);
    *returnstr=(char *)realloc(*returnstr,30*sizeof(char));
    sprintf(*returnstr,"%E",((function_args_t *)codex->function_args)->x[cnt]);
  }

}

int main(){
  char code[]=""
              " sum=0;"
              " for (cnt,0,30,"
              "   sum=sum+x(cnt);"
              " );"
              " writeln(sum);"
              "";
  SOAP_codex_t codex;
  long cnt;
  function_args_t function_args;
  double x[31];

  for (cnt=0; cnt<=30; cnt++) x[cnt]=(double)cnt;
  function_args.x=x;

  SOAP_init_codex(&codex,"?");
  codex.FUNCTION=TRUE;
  codex.function=&functions;
  codex.function_args=(void *)&function_args;
  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
  SOAP_free_codex(&codex);

  return(0);
}
