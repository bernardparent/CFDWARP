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


void add_Z_dUstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  long l;
  long flux;
  flux_t tmp,dRes;
  sqmat_t Z;

    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      for (flux=0; flux<nf; flux++){
        tmp[flux]=_Omega(np[l],gl)*(np[l].bs->U[flux]
                         -np[l].bs->Um1[flux]
                         )/gl->dt;
      }
      find_Z(np,gl,l,Z);
      multiply_matrix_and_vector(Z, tmp, dRes);
      for (flux=0; flux<nf; flux++) {
        np[l].wk->Res[flux]+=fact*dRes[flux];
      }
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
      for (flux=0; flux<nf; flux++) {
        np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*dRes[flux];
      }
#endif
    }
}


