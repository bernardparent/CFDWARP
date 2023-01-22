#ifndef _CONTROL_H
#define _CONTROL_H


#include <src/common.h>
#include <model/_model.h>
#include <cycle/_cycle.h>
#include <src/data.h>

typedef struct {
  np_t **np;
  //np_t npextranode;
  //double Aextranode;
  gl_t *gl;
  np_t *np_post;
  gl_t gl_post;
  zone_t domain_post;
  input_t *input;
  double CFLinit;
  bool VERBOSE,POSTMODULE,CYCLEMODULE,RESETITERCOUNT,GRIDONLY;
  //bool EXTRANODE;
  long module_level;
  long TYPELEVEL;
} readcontrolarg_t;

void read_control_functions(char *functionname, char **argum,
                           char **returnstr, SOAP_codex_t *codex);

void write_control(char *filename);

void read_control(char *control_filename, input_t input, bool CYCLEMODULE, bool POSTMODULE, bool GRIDONLY,
                 bool RESETITERCOUNT, double CFLinit, np_t **np, gl_t *gl);

double _x_DISTMPI_compatible(np_t *np, gl_t *gl, long i, long j, long k, long dim);

void write_modules(FILE *outputfile);

void write_license ( FILE * outputfile );

void  write_hline(FILE *outputfile, long indent, long width);

void find_metrics_on_all_nodes(np_t *np, gl_t *gl, zone_t zone);

#endif /* _CONTROL_H */
