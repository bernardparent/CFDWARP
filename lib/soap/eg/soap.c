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
