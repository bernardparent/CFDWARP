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

#include <gridg.h>
#include <string.h>
#include <stdlib.h>

#define INSERT_LINE_NUMBERS TRUE

void PrintGrid1D ( GRIDG_gl1d_t gl1d, GRIDG_xgrid_t * xgrid ) {
  long i;
  
  printf( "\nNode\tX\t\tAcs\n" );
  printf( "------------------------------------\n" );
  
  for ( i = gl1d.is; i <= gl1d.ie; i++ ) {
      printf ( "%ld\t%.5E\t%.5E\n", 
                i,
                xgrid[GRIDG_ai1 ( gl1d, i )].x,
                xgrid[GRIDG_ai1 ( gl1d, i )].Acs );
  }
  printf( "\n" );
}

int main ( int argc, char **argv ) {
  GRIDG_xgrid_t *xgrid;
  GRIDG_gl1d_t gl1d;
  SOAP_codex_t codex;
  char *code;

  const char *GRID_DEFINITION = 
    "  is=1;\n"
    "  ie=40;\n"
    "  Size ( is, ie );\n"
    "  Point( is, 0.0, 1.0 );\n"
    "  Point( ie, 0.3, 1.0 );\n"
    "  JoinCorners( is, ie, FE, 0.5, 0.001, 1.0 );\n";

  if ( argc == 1 ) {
    printf ( "Generating grid internally..\n" );
    
    code = ( char * ) malloc ( (strlen(GRID_DEFINITION) + 1) * sizeof ( char ) );
    strcpy ( code, GRID_DEFINITION );
    
    SOAP_init_codex ( &codex, "internal_grid_string" );
    
    if ( INSERT_LINE_NUMBERS ) {
      SOAP_insert_line_numbers_in_code ( &code, 1 );
    }
    
    // Directly process the un-wrapped string using the library function
    GRIDG_read_grid_1D_from_argum ( code, &codex, &gl1d, &xgrid );
    
    SOAP_free_codex ( &codex );
    free ( code );
    
    PrintGrid1D ( gl1d, xgrid );
    
  } else {
    fprintf ( stderr, "printgrid1D does not require any arguments. \n"
              "example: ./printgrid1D \n" );
    exit ( 1 );
  }
  return ( 0 );
}
