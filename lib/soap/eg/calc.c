#include <soap.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

int main ( int argc, char **argv ) {
  SOAP_codex_t codex;
  char *code, *format;

  if ( argc >= 2 ) {
    code = ( char * ) malloc ( 4000 * sizeof ( char ) );
    strcpy ( code, argv[argc - 1] );
    SOAP_init_codex ( &codex, "?" );

    if ( argc > 2 ) {
      format = ( char * ) malloc ( 4000 * sizeof ( char ) );
      strcpy ( format, argv[1] );
      SOAP_strins ( "\",", &code, 0 );
      SOAP_strins ( format, &code, 0 );
      SOAP_strins ( "printf(\"", &code, 0 );
    } else {
      SOAP_strins ( "writeln(", &code, 0 );
    }
    SOAP_strins ( ");", &code, ( long ) strlen ( code ) );
    SOAP_process_code ( code, &codex, SOAP_VARS_KEEP_ALL );
    SOAP_free_codex ( &codex );
    return ( EXIT_SUCCESS );
  } else {
    fprintf ( stderr, "calc needs at least one argument: the string containing the arithmetic. \n"
              "examples:\n"
              "  calc '10+10' will return 20. \n"
              "  calc '%%15.15E\\n' '10+10' will return 20 formatted according to %%15.15E.\n" );
  }
  return ( EXIT_FAILURE );
}
