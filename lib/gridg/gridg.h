#ifndef _GRIDG_H
#define _GRIDG_H

#include <memdebug.h>
#include <soap.h>
#include <exm.h>
#include <stdio.h>

typedef struct {
  double x,Acs;
  bool INIT;
} GRIDG_xgrid_t;

typedef struct {
  double x,y;
  bool INIT;
} GRIDG_xygrid_t;

typedef struct {
  double x,y,z;
  bool INIT;
} GRIDG_xyzgrid_t;

typedef struct {
  long is,ie,js,je,ks,ke;
} GRIDG_gl3d_t;

typedef struct {
  long is,ie,js,je;
} GRIDG_gl2d_t;

typedef struct {
  long is,ie;
} GRIDG_gl1d_t;

long GRIDG_ai3(GRIDG_gl3d_t gl, long i, long j, long k);

long GRIDG_ai2(GRIDG_gl2d_t gl, long i, long j);

long GRIDG_ai1(GRIDG_gl1d_t gl, long i);

void GRIDG_read_grid_1D_from_file(char *filename, GRIDG_gl1d_t *gl1d, GRIDG_xgrid_t **xgrid,
                bool VERBOSE, bool *Problem);

void GRIDG_read_grid_1D_from_argum(char *argum, SOAP_codex_t *codex, GRIDG_gl1d_t *gl1d, GRIDG_xgrid_t **xgrid);

void GRIDG_write_grid_1D_to_file(FILE **controlfile);


void GRIDG_read_grid_2D_from_file(char *filename, GRIDG_gl2d_t *gl2d, GRIDG_xygrid_t **xygrid,
                bool VERBOSE, bool *Problem);

void GRIDG_read_grid_2D_from_argum(char *argum, SOAP_codex_t *codex, GRIDG_gl2d_t *gl2d, GRIDG_xygrid_t **xygrid);


void GRIDG_write_grid_2D_to_file(FILE **controlfile);


void GRIDG_read_grid_3D_from_file(char *filename, GRIDG_gl3d_t *gl3d, GRIDG_xyzgrid_t **xyzgrid,
                bool VERBOSE, bool *Problem);

void GRIDG_read_grid_3D_from_argum(char *argum, SOAP_codex_t *codex, GRIDG_gl3d_t *gl3d, GRIDG_xyzgrid_t **xyzgrid);

void GRIDG_write_grid_3D_to_file(FILE **controlfile);


#endif /* _GRIDG_H */



