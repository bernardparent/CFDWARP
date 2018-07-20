#ifndef _RESTIME_H
#define _RESTIME_H

#include <src/common.h>
#include <model/_model.h>


#ifdef _RESTIME_CDF

void find_Fstarxt_interface(np_t *np, gl_t *gl, long theta, long l, jacvars_t jacvars, metrics_t metrics, flux_t dFstar);

void find_dFstarxt_dUstar_interface(np_t *np, gl_t *gl, long theta, long l, jacvars_t jacvarsp1h, metrics_t metricsp1h, sqmat_t B, sqmat_t C);

void find_Lambdaxt(np_t *np, gl_t *gl, long l, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxt);

void find_Lambdaxt_interface(np_t *np, gl_t *gl, long lp0, long lp1, long theta, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxt);

void find_Deltaxt_interface(np_t *np, gl_t *gl, long lp0, long lp1, long theta, 
                            jacvars_t jacvarsp0, jacvars_t jacvarsp1,  metrics_t metrics, flux_t Deltaxt);

/* jacvars and metrics must be evaluated between node lp0 and node lp1 */
void find_Lambdaxt_plus_dtau_from_jacvars(np_t *np, gl_t *gl, long lp0, long lp1, long theta, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxtplus);

/* jacvars and metrics must be evaluated between node lp0 and node lp1 */
void find_Lambdaxt_minus_dtau_from_jacvars(np_t *np, gl_t *gl, long lp0, long lp1, long theta, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxtminus);

void find_Lambdaxt_plus_dtau(np_t *np, gl_t *gl, long l, long theta, sqmat_t Lambdaxtplus);

void find_Lambdaxt_minus_dtau(np_t *np, gl_t *gl, long l, long theta, sqmat_t Lambdaxtminus);

#endif

void write_disc_restime_template(FILE **controlfile);

void read_disc_restime_actions(char *actionname, char **argum, SOAP_codex_t *codex);

void add_Z_dUstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal);

void update_U_predictor_corrector(np_t *np, gl_t *gl, zone_t zone);


#endif /* _RESTIME_H */ 
