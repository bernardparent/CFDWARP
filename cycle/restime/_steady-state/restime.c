#include <cycle/restime/_restime.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>


void write_disc_restime_template(FILE **controlfile){
}


void read_disc_restime_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_disc_restime_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->DISC_RESTIME_READ=TRUE;
}


void find_Fstarxt_interface(np_t *np, gl_t *gl, long theta, long l, jacvars_t jacvars, metrics_t metrics, flux_t dFstar){
  long flux; 
  for (flux=0; flux<nf; flux++){
    dFstar[flux]=0.0;
  }
}


void add_Z_dUstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
}

