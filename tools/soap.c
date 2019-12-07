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

#include <soap.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>




int main( int argc, char **argv ){
 SOAP_codex_t codex;
 char *code,*filename;

 if (argc==2) {
   filename=(char *)malloc(4000*sizeof(char));
   strcpy(filename,argv[1]);
   code=(char *)malloc(sizeof(char));
   SOAP_store_file_as_string(filename,&code);
   /* delete first line if it starts by #*/
   if (code[0]=='#') {
     do SOAP_strcut(0,0,code); while (code[0]!='\n');
   }
   SOAP_init_codex(&codex,filename);
   SOAP_insert_line_numbers_in_code(&code,1);
   SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
   SOAP_free_codex(&codex);
   free(filename);
   free(code);
   return(1);
 } else {
   fprintf(stderr,"soap needs one argument: the filename of the file containing the code. \n"
                  "example: soap test.soap  \n");
 }
  return(0);
}
