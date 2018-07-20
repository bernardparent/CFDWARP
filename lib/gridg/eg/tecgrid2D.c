#include <gridg.h>
#include <string.h>
#include <stdlib.h>

#define INSERT_LINE_NUMBERS TRUE

typedef struct {
  GRIDG_xygrid_t **xygrid;
  GRIDG_gl2d_t *gl2d;
} action_args_t;

void WritePost2D ( char *filename, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t * xygrid ) {
  long i, j;
  FILE *postfile;

  postfile = fopen ( filename, "w" );
  fprintf ( postfile, "VARIABLES=\"X\",\"Y\"\n" );
  fprintf ( postfile, "ZONE I=%ld, J=%ld F=POINT\n", gl2d.ie - gl2d.is + 1, gl2d.je - gl2d.js + 1 );
  for ( j = gl2d.js; j <= gl2d.je; j++ ) {
    for ( i = gl2d.is; i <= gl2d.ie; i++ )
      fprintf ( postfile, "%.5E  %.5E\n", xygrid[GRIDG_ai2 ( gl2d, i, j )].x,
                xygrid[GRIDG_ai2 ( gl2d, i, j )].y );
  }
  fclose ( postfile );
}

void actions ( char *actionname, char **argum, SOAP_codex_t * codex ) {
  GRIDG_xygrid_t **xygrid;
  GRIDG_gl2d_t *gl2d;
  xygrid = ( ( action_args_t * ) codex->action_args )->xygrid;
  gl2d = ( ( action_args_t * ) codex->action_args )->gl2d;

  if ( strcmp ( actionname, "Grid" ) == 0 ) {
    printf ( "  Grid.." );
    GRIDG_read_grid_2D_from_argum ( *argum, codex, gl2d, xygrid );
    codex->ACTIONPROCESSED = TRUE;
    printf ( "[done]\n" );
  }
}

int main ( int argc, char **argv ) {
  GRIDG_xygrid_t *xygrid;
  GRIDG_gl2d_t gl2d;
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
      action_args.xygrid = &xygrid;
      action_args.gl2d = &gl2d;
      codex.action_args = ( void * ) &action_args;
      SOAP_process_code ( code, &codex, SOAP_VARS_KEEP_ALL );
      SOAP_free_codex ( &codex );
    } else {
      GRIDG_read_grid_2D_from_file ( infilename, &gl2d, &xygrid, TRUE, &PROBLEM );
    }
    WritePost2D ( outfilename, gl2d, xygrid );
    printf ( "[done]\n" );
  } else {
    fprintf ( stderr, "tecgrid needs two arguments: the file containing the grid and the postfile. \n"
              "example: tecgrid2D gridfile postfile \n" );
    exit ( 1 );
  }
  return ( 0 );
}
