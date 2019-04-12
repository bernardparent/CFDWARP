#include <cycle/resconv/_resconv.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <cycle/share/res_share.h>
#include <src/control.h>

#if (_FLUID_CONVECTION)
  #error The fluid module has convection terms: choose a scheme other than "none" for the convection terms discretization
#endif


void write_disc_resconv_template(FILE **controlfile){

}


void read_disc_resconv_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){

}


void read_disc_resconv_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->DISC_RESCONV_READ=TRUE;
  gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_NONE;
}


void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
}



void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  sqmat_t Lambda_minus,Lambda_plus;
  long flux;
  find_Lambda_minus_dtau_FVS(np,gl,l,theta,EIGENVALCOND_NONE,Lambda_minus);
  find_Lambda_plus_dtau_FVS(np,gl,l,theta,EIGENVALCOND_NONE,Lambda_plus);
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=Lambda_plus[flux][flux]-Lambda_minus[flux][flux];
}

