#include <cycle/ts/_ts.h>
#include <cycle/share/cycle_share.h>

#ifdef _CYCLE_PREDICTOR_CORRECTOR
#error the explicit relaxation is not compatible with the predictor-corrector cycle
#endif


static void find_dUtil(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  flux_t dUstar;
  double dtau;
  sqmat_t GammaInv;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) 
      find_Gamma_inverse(np, gl, l, GammaInv);
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
      multiply_matrix_and_vector(GammaInv,np[l].wk->Res,dUstar);
      for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]=-dUstar[flux];
    } else {
      find_constant_dtau(np,gl,l,&dtau);
      for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]=-np[l].wk->Res[flux]*dtau;
    }
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }
}




void update_Util(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    add_dUstar_to_U(np,l,gl,np[l].wk->dUstar);
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    thread_lock_global_set(gl,THREADTYPE_ALL);
    gl->effiter_U+=1.0/(double)(gl->nn);
    thread_lock_global_unset(gl,THREADTYPE_ALL);
  }
}





void update_U(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np,gl,zone,&find_dUtil,SWEEPTYPE_I,TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_Util,SWEEPTYPE_I,TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
}


void find_dU(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np,gl,zone,&find_dUtil,SWEEPTYPE_I,TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
}

