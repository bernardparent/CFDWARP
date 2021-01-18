#ifndef _SHARE_H
#define _SHARE_H

#include "gridg.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define EOS 0

int find_segment_NO(double *seg, long N, double garb0,
                   double L, double garb1, double garb2, SOAP_codex_t *codex);

/* the following return 0 on success and 1 on failure */
int find_segment_1(double *seg, long N, double NLoverN,
                 double L, double epsL, double epsH, SOAP_codex_t *codex);

int find_segment_2(double *seg, long N, double NLoverN,
                 double L, double dzeta1, double epsL, SOAP_codex_t *codex);

int find_segment_3(double *seg, long N, double NLoverN,
                 double L, double dzeta1, double dzeta2, SOAP_codex_t *codex);

int find_segment_4(double *seg, long N, double NLoverN,
                 double L, double dzeta1, double epsH, SOAP_codex_t *codex);

int find_segment_5(double *seg, long N, double NLoverN,
                 double L, double epsL, double dzeta2, SOAP_codex_t *codex);

int find_x_along_i(GRIDG_xgrid_t *xgrid, double segval1, double segval2,
                 double NLoverN,
                 GRIDG_gl1d_t gl, long i1, long i2, SOAP_codex_t *codex, 
                 int(*find_segment)(double *, long , double,
                                    double , double , double, SOAP_codex_t * ));

int find_xy_along_ij(GRIDG_xygrid_t *xygrid, double segval1, double segval2,
                   double NLoverN,
                   GRIDG_gl2d_t gl, long i1, long j1, long i2, long j2, SOAP_codex_t *codex, 
                   int(*find_segment)(double *, long , double,
                                      double , double , double, SOAP_codex_t * ));

int find_xyz_along_ijk(GRIDG_xyzgrid_t *xyzgrid, double segval1, double segval2,
                     double NLoverN,
                     GRIDG_gl3d_t gl, long i1, long j1, long k1, long i2, long j2, long k2, SOAP_codex_t *codex, 
                     int(*find_segment)(double *, long , double,
                         double , double , double , SOAP_codex_t *));

void erode_bdry(GRIDG_xygrid_t *xygrid, long erofact,
               GRIDG_gl2d_t gl, long i1, long j1, long i2, long j2);

double read_double_after_X_chars(FILE *infile, long X, char *varname, bool VERBOSE);

long read_long_after_X_chars(FILE *infile, long X, char *varname, bool VERBOSE);

void grab_code(FILE *infile, char **code, long *codelength);

long sgndif(long i1, long i2);

/* returns 0 on success, 1 on failure */
int find_segment_type(int (**FindSegment)(double *, long , double, double , double , double, SOAP_codex_t * ),
                     char modeL, char modeH, SOAP_codex_t *codex);

void replace_equal_sign(char **expr);

/* will open file named *filename, look for an action
   named *actionname using the interpret.c library and will
   associate the function *RXaction with optional arguments
   *RXaction_args to *action. Should *actionname not be present
   in *filename, a module will be appended using the
   function *WriteRXaction, and then read.*/
void read_RX_grid(char *filename, char *actionname,
                void (*RXaction)(char **, void *, SOAP_codex_t *),
                void (*WriteRXaction)(FILE **),
                void *RXaction_args, bool VERBOSE);

void GRIDG_fatal_error(const char *formatstr, ...);


#endif
