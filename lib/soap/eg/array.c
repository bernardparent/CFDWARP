
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
