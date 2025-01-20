#ifndef _METRICS_H
#define _METRICS_H

#include <src/common.h>

#define BDRYMETRICS_NORMAL 0
#define BDRYMETRICS_CENTERED 1


void write_metrics_template(FILE **controlfile);

void read_metrics(char *argum, SOAP_codex_t *codex);

double _X(np_t np, long dim1, long dim2);

double _x(np_t np, long dim);

void update_metrics_at_node(np_t *np, gl_t *gl, long l);

void update_metrics_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta);

void find_Omega_and_X_at_node(np_t *np, gl_t *gl, long l, double *Omega, dim2_t X);

void find_metrics_at_node(np_t *np, gl_t *gl, long l,
                               long theta, metrics_t *metrics);


void find_Omega_and_X_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta,
                          double *Omega, dim2_t X);

void find_metrics_at_interface(np_t *np, gl_t *gl, long lL, long lR,
                               long theta, metrics_t *metrics);

void find_metrics_at_interface_2(np_t *np, gl_t *gl, long lL, long lR,
                               long theta, long theta2, metrics_t *metrics);

void update_metrics_at_bdry_node(np_t *np, gl_t *gl, int TYPELEVEL, long l_A, long l_B, long l_C,
                       int BDRYMETRICS);

void find_metrics_from_base_level(np_t *np, gl_t *gl, long l, metrics_t *metrics);

double _Omega(np_t np, gl_t *gl);

bool find_unit_vector_normal_to_boundary_plane(np_t *np, gl_t *gl, long lA, long lB, long lC, int TYPELEVEL, dim_t n);

bool find_distance_between_near_bdry_node_and_boundary_plane(np_t *np, gl_t *gl, long lA, long lB, long lC, int TYPELEVEL, double *dwall);
		   
#ifdef _2D
void find_side_projected_area_of_axisymmetric_cell(np_t *np, gl_t *gl, long l, int METRICSRMIN, dim_t projarea);
#endif
           
#endif
