#ifndef _CYCLE_H
#define _CYCLE_H

#include <model/_model.h>
#include <src/common.h>
#include <src/control.h>



void find_dtau(np_t *np, gl_t *gl, long l, flux_t dtau);

void find_constant_dtau(np_t *np, gl_t *gl, long l, double *dtaumin);

void read_cycle(char *argum, SOAP_codex_t *codex);

void read_disc(char *argum, SOAP_codex_t *codex);

void init_cycle(char *argum, SOAP_codex_t *codex);

void write_cycle_template(FILE **controlfile);

void write_disc_template(FILE **controlfile);

void perform_one_iteration(np_t *np, gl_t *gl);

void resume_nodes_only_in_zone_and_update_bdry_nodes(np_t *np, gl_t *gl, zone_t zone);

void resume_nodes_only_in_zone(np_t *np, gl_t *gl, zone_t zone);

void resume_nodes_in_zone(np_t *np, gl_t *gl, zone_t zone);

void resume_nodes_specified_in_function(np_t *np, gl_t *gl,
                                bool(*FUNCT)(gl_t *, long, long, long));

void check_residual(np_t *np, gl_t *gl, zone_t zone);

void write_runtime_template(FILE **controlfile);

void runtime_actions_cycle_specific(char *actionname, char **argum, SOAP_codex_t *codex);


#endif /* _CYCLE_H */ 
