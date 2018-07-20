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
