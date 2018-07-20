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
