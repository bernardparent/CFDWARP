
#include <soap.h>


int main(){
  char code[]=""
              " k=0;"
              " while (k<=10,"
              "   writeln(k);"
              "   k=k+1;"
              " );"
              "";
  SOAP_codex_t codex;

  SOAP_init_codex(&codex,"?");
  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
  SOAP_free_codex(&codex);

  return(0);
}
