
#include <soap.h>


int main(){
  char code[]=""
              " k=tulip;"
              " while (k<=10,"
              "   writeln(k*rose);"
              "   k=k+1;"
              " );"
              " writeln(x+10);"
              "";
  SOAP_codex_t codex;

  SOAP_init_codex(&codex,"?");
  SOAP_add_to_vars(&codex,"rose","10.4");
  SOAP_add_to_vars(&codex,"tulip","5");
  SOAP_add_to_vars(&codex,"x","10.0e0");

  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
  SOAP_free_codex(&codex);

  return(0);
}
