#include <gridg.h>
#include <string.h>
#include <stdlib.h>

#define INSERT_LINE_NUMBERS TRUE

typedef struct {
  GRIDG_xyzgrid_t **xyzgrid;
  GRIDG_gl3d_t *gl3d;
} action_args_t;

void WritePost3D ( char *filename, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t * xyzgrid ) {
  long i, j, k;
  FILE *postfile;

  postfile = fopen ( filename, "w" );
  fprintf ( postfile, "VARIABLES=\"X\",\"Y\",\"Z\"\n" );
  fprintf ( postfile, "ZONE I=%ld, J=%ld, K=%ld F=POINT\n",
            gl3d.ie - gl3d.is + 1, gl3d.je - gl3d.js + 1, gl3d.ke - gl3d.ks + 1 );
  for ( k = gl3d.ks; k <= gl3d.ke; k++ ) {
    for ( j = gl3d.js; j <= gl3d.je; j++ ) {
      for ( i = gl3d.is; i <= gl3d.ie; i++ )
        fprintf ( postfile, "%.5E  %.5E  %.5E\n", xyzgrid[GRIDG_ai3 ( gl3d, i, j, k )].x,
                  xyzgrid[GRIDG_ai3 ( gl3d, i, j, k )].y, xyzgrid[GRIDG_ai3 ( gl3d, i, j, k )].z );
    }
  }
  fclose ( postfile );
}

void actions ( char *actionname, char **argum, SOAP_codex_t * codex ) {
  GRIDG_xyzgrid_t **xyzgrid;
  GRIDG_gl3d_t *gl3d;
  xyzgrid = ( ( action_args_t * ) codex->action_args )->xyzgrid;
  gl3d = ( ( action_args_t * ) codex->action_args )->gl3d;

  if ( strcmp ( actionname, "Grid" ) == 0 ) {
    printf ( "  Grid.." );
    GRIDG_read_grid_3D_from_argum ( *argum, codex, gl3d, xyzgrid );
    printf ( "[done]\n" );
    codex->ACTIONPROCESSED = TRUE;
  }
}

int main ( int argc, char **argv ) {
  GRIDG_xyzgrid_t *xyzgrid;
  GRIDG_gl3d_t gl3d;
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
      action_args.xyzgrid = &xyzgrid;
      action_args.gl3d = &gl3d;
      codex.action_args = ( void * ) &action_args;
      SOAP_process_code ( code, &codex, SOAP_VARS_KEEP_ALL );
      SOAP_free_codex ( &codex );
    } else {
      GRIDG_read_grid_3D_from_file ( infilename, &gl3d, &xyzgrid, FALSE, &PROBLEM );
    }
    WritePost3D ( outfilename, gl3d, xyzgrid );
    printf ( "[done]\n" );
  } else {
    fprintf ( stderr, "tecgrid needs two arguments: the file containing the grid and the postfile. \n"
              "example: tecgrid3D gridfile postfile \n" );
    exit ( 1 );
  }
  return ( 0 );
}
