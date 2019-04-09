#include <cycle/restime/_restime.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>

#ifndef _RESCONV_MUSCL
#error The Trapezoidal+ time discretization requires the MUSCL spatial discretization 
#endif

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
  flux_t tmp,dRes,Up1h,Um1h;
  sqmat_t Z;

    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      for (flux=0; flux<nf; flux++) {
        Up1h[flux]= 
              +1.0*np[l].bs->U[flux]; 
        Um1h[flux]= 
              +1.0*np[l].bs->Um1[flux]; 
      }
        
      for (flux=0; flux<nf; flux++){
        tmp[flux]=_Omega(np[l],gl)*(Up1h[flux]-Um1h[flux])/gl->dt;
      }
      
      find_Z(np,gl,l,Z);
      multiply_matrix_and_vector(Z, tmp, dRes);
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=fact*dRes[flux];
    }
}



