// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2026 Bernard Parent

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

typedef struct {
  GRIDG_xgrid_t **xgrid;
  GRIDG_gl1d_t *gl1d;
} action_args_t;

void WritePost1D ( char *filename, GRIDG_gl1d_t gl1d, GRIDG_xgrid_t * xgrid ) {
  long i;
  FILE *postfile;

  postfile = fopen ( filename, "w" );
  // GRIDG_xgrid_t contains x and Acs (cross-sectional area)
  fprintf ( postfile, "VARIABLES=\"X\",\"Acs\"\n" );
  fprintf ( postfile, "ZONE I=%ld F=POINT\n", gl1d.ie - gl1d.is + 1 );
  
  for ( i = gl1d.is; i <= gl1d.ie; i++ ) {
      fprintf ( postfile, "%.5E  %.5E\n", 
                xgrid[GRIDG_ai1 ( gl1d, i )].x,
                xgrid[GRIDG_ai1 ( gl1d, i )].Acs );
  }
  fclose ( postfile );
}

void actions ( char *actionname, char **argum, SOAP_codex_t * codex ) {
  GRIDG_xgrid_t **xgrid;
  GRIDG_gl1d_t *gl1d;
  
  xgrid = ( ( action_args_t * ) codex->action_args )->xgrid;
  gl1d = ( ( action_args_t * ) codex->action_args )->gl1d;

  if ( strcmp ( actionname, "Grid" ) == 0 ) {
    printf ( "  Grid.." );
    GRIDG_read_grid_1D_from_argum ( *argum, codex, gl1d, xgrid );
    codex->ACTIONPROCESSED = TRUE;
    printf ( "[done]\n" );
  }
}

int main ( int argc, char **argv ) {
  GRIDG_xgrid_t *xgrid;
  GRIDG_gl1d_t gl1d;
  bool PROBLEM;
  char infilename[100], outfilename[100];
  SOAP_codex_t codex;
  char *code;
  action_args_t action_args;

  if ( argc == 3 ) {
    strcpy ( infilename, argv[1] );
    strcpy ( outfilename, argv[2] );
    printf ( "Reading grid file %s..\n", infilename );
    
    if ( INSERT_LINE_NUMBERS ) {
      code = ( char * ) malloc ( sizeof ( char ) );
      SOAP_store_file_as_string ( infilename, &code );
      SOAP_init_codex ( &codex, infilename );
      SOAP_insert_line_numbers_in_code ( &code, 1 );
      
      codex.ACTION = TRUE;
      codex.action = &actions;
      action_args.xgrid = &xgrid;
      action_args.gl1d = &gl1d;
      codex.action_args = ( void * ) &action_args;
      
      SOAP_process_code ( code, &codex, SOAP_VARS_KEEP_ALL );
      SOAP_free_codex ( &codex );
    } else {
      GRIDG_read_grid_1D_from_file ( infilename, &gl1d, &xgrid, TRUE, &PROBLEM );
    }
    
    WritePost1D ( outfilename, gl1d, xgrid );
    printf ( "[done]\n" );
  } else {
    fprintf ( stderr, "tecgrid1D needs two arguments: the file containing the grid and the postfile. \n"
              "example: tecgrid1D gridfile postfile \n" );
    exit ( 1 );
  }
  return ( 0 );
}
