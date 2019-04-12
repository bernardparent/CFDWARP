#include <cycle/restime/_restime.h>
#include <cycle/share/cycle_share.h>
#include <src/control.h>

#ifndef _TS_NONE
  #error RK4 can not be compiled with a fluid relaxation scheme other than none 
#endif

#ifndef _CYCLE_PREDICTOR_CORRECTOR
  #error RK4 can not be compiled with a cycle other than predictor_corrector 
#endif


void write_disc_restime_template(FILE **controlfile){
}


void read_disc_restime_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_disc_restime_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->DISC_RESTIME_READ=TRUE;
  gl->numsubiter_pc=4;
}


void update_U_local(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  flux_t dUstar;
  double Omega;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    Omega=_Omega(np[l],gl);
    for (flux=0; flux<nf; flux++) {
      switch (gl->subiter_pc){
        case 0:
          np[l].bs->dUstar1[flux]=-np[l].wk->Res[flux]*gl->dt;
          dUstar[flux]=Omega*np[l].bs->Um1[flux]+0.5*np[l].bs->dUstar1[flux]-Omega*np[l].bs->U[flux];
        break;
        case 1:
          np[l].bs->dUstar2[flux]=-np[l].wk->Res[flux]*gl->dt;
          dUstar[flux]=Omega*np[l].bs->Um1[flux]+0.5*np[l].bs->dUstar2[flux]-Omega*np[l].bs->U[flux];
        break;
        case 2:
          np[l].bs->dUstar3[flux]=-np[l].wk->Res[flux]*gl->dt;
          dUstar[flux]=Omega*np[l].bs->Um1[flux]+np[l].bs->dUstar3[flux]-Omega*np[l].bs->U[flux];
        break;
        case 3:
          np[l].bs->dUstar4[flux]=-np[l].wk->Res[flux]*gl->dt;
          dUstar[flux]=Omega*np[l].bs->Um1[flux]+(np[l].bs->dUstar1[flux]+2.0*np[l].bs->dUstar2[flux]+2.0*np[l].bs->dUstar3[flux]+np[l].bs->dUstar4[flux])/6.0-Omega*np[l].bs->U[flux];
        break;
        default:
          fatal_error("gl->subiter can not be set to %ld for RK4.",gl->subiter_pc);
      }
    }
    add_dUstar_to_U(np,l,gl,dUstar);
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    thread_lock_global_set(gl,THREADTYPE_ALL);
    gl->effiter_U+=1.0/(double)(gl->nn);
    thread_lock_global_unset(gl,THREADTYPE_ALL);
  }
}



void update_U_predictor_corrector(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np,gl,zone,&update_U_local,SWEEPTYPE_I,TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
}




void add_Z_dUstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){

}


