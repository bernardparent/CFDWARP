// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <soap.h>


int main(){
  char code[]=""
              " mike=10;"
              " tmp[mike]=10;"
              " array[10][0]=tulip;"
              " array[11][0]=tulip;"
              " array[12][0]=tulip;"
              " array[13][0]=tulip;"
              " k=10;"
              " bern=0;"
              " while (k<=13,"
              "   writeln(array[(k+1)/10*10-1][bern]);"
              "   k=k+1;"
              " );"
              " writeln(10^array[tmp[mike]][0]);"
              "";
  SOAP_codex_t codex;

  SOAP_init_codex(&codex,"?");
  SOAP_add_to_vars(&codex,"rose","10.4");
  SOAP_add_to_vars(&codex,"tulip","5");
  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
  SOAP_free_codex(&codex);

  return(0);
}
