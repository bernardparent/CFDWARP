#ifndef _POST_H
#define _POST_H


#include <src/common.h>
#include <model/_model.h>

/* from the domain defined by *np and *gl, and the array *xcut specifying
   the x-positions of the cuts, with numcut the number of cuts, find
   a new domain described by npcut and glcut which will have numcut
   i stations.
*/
void create_domain_of_cuts_along_x(np_t *np, gl_t *gl, np_t **npcut, gl_t *glcut,
                              double *xcut, long numcut);


#ifdef _2DL
void create_domain_of_cuts_along_y(np_t *np, gl_t *gl, np_t **npcut, gl_t *glcut,
                              double *ycut, long numcut);
#endif


#ifdef _3D
void create_domain_of_cuts_along_z(np_t *np, gl_t *gl, np_t **npcut, gl_t *glcut,
                              double *zcut, long numcut);
#endif

/* given the domain specified by np and gl, output the nodes part of zone
   in a post-processor file compatible with the post-processor specified
   in the string *postprocessor (currently, only nodplot and tecplot
   are possible). The post-processor file is named *filename
*/
void write_post_file(np_t *np, gl_t *gl, zone_t zone, char *filename,
                   char *postprocessor, bool GRIDONLY);


void read_post(char *argum, SOAP_codex_t *codex);

void write_post_template(FILE **controlfile);


#endif /* _POST_H */
