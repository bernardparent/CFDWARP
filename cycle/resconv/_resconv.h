#ifndef _RESCONV_H
#define _RESCONV_H

#include <src/common.h>
#include <model/_model.h>

void init_Fstar_interfaces(np_t *np, gl_t *gl, zone_t zone);

void find_Fstar_interfaces(np_t *np, gl_t *gl, long theta, long ls, long le);

void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl);

void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda);

void write_disc_resconv_template(FILE **controlfile);

void read_disc_resconv_actions(char *actionname, char **argum, SOAP_codex_t *codex);

#endif /* _RESCONV_H */ 
