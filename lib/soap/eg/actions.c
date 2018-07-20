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

