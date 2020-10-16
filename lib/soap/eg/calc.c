// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2000 Bernard Parent

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
    // add some physical constants in SI units
    SOAP_add_to_vars(&codex, "kB", "1.380649E-23");
    SOAP_add_to_vars(&codex, "me", "9.10938356E-31");
    SOAP_add_to_vars(&codex, "e", "1.60217662E-19");

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
