// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <soap.h>
#include <string.h>

void functions(char *functionname, char **argum, char **returnstr, SOAP_codex_t *codex){
  double base,expon;

  if (strcmp(functionname,"pow")==0) {
    SOAP_substitute_all_argums(argum,codex);
    sscanf(*argum,"%lg,%lg",&base,&expon);
    *returnstr=(char *)realloc(*returnstr,30*sizeof(char));
    sprintf(*returnstr,"%E",pow(base,expon));
  }

}

int main(){
  char code[]=""
              " k=0;"
              " while (k<=10,"
              "   writeln(k,pow(k,2),k^2);"
              "   k=k+1;"
              " );"
              "";
  SOAP_codex_t codex;

  SOAP_init_codex(&codex,"?");
  codex.FUNCTION=TRUE;
  codex.function=&functions;
  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
  SOAP_free_codex(&codex);

  return(0);
}
