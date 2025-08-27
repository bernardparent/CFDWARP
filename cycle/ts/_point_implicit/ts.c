// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2025 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
   of conditions and the following disclaimer in the documentation and/or other
   materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cycle/share/cycle_share.h>
#include <cycle/ts/_ts.h>
#include <cycle/share/ts_share.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>

#ifdef _CYCLE_PREDICTOR_CORRECTOR
#error the block ADI relaxation is not compatible with the predictor-corrector cycle
#endif

#define alpha 1.0e0

#define fluxmass (fluxmom-1)
#define nfr (nf-fluxmass)



void update_dUtilde(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux,row,col,is,ie,js,je,ks,ke;
  sqmat_t dS_dU;

  find_ijk_from_l(gl, ls, &is, &js, &ks);
  find_ijk_from_l(gl, le, &ie, &je, &ke);

  EXM_mat_t matchem,matcheminv,matdU1,matdU2;

  if (TRUE){
    EXM_init_matrix(&matchem, nf, nf);
    EXM_init_matrix(&matcheminv, nf, nf);
    EXM_init_matrix(&matdU1, nf, 1);
    EXM_init_matrix(&matdU2, nf, 1);
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      thread_lock_node_set(np,l,THREADTYPE_ZONE);
      find_dSstar_dUstar(np,gl,l,dS_dU);
//      find_dSchem_dU(np, gl, l, dS_dU); 
      if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
        fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in update_dUtilde() within ts.c.");     
      }
      for (row=0; row<nf; row++){
        for (col=0; col<nf; col++){
          matchem.cont[EXM_aim(matchem.glm,row,col)]=-dS_dU[row][col]*alpha*np[l].wk->dtau; 
        }
      }
      for (row=0; row<nf; row++) matchem.cont[EXM_aim(matchem.glm,row,row)]+=1.0e0;         

      EXM_invert_matrix_partial_pivoting(matchem, &matcheminv);
      for (flux=0; flux<nf; flux++) matdU1.cont[EXM_aim(matdU1.glm,flux,0)]=np[l].wk->dUstar[flux];
      EXM_multiply_matrices(matcheminv, matdU1, &matdU2);
      for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]=matdU2.cont[EXM_aim(matdU2.glm,flux,0)];
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
    EXM_free_matrix(&matchem);
    EXM_free_matrix(&matcheminv);
    EXM_free_matrix(&matdU1);
    EXM_free_matrix(&matdU2);
  }
}





void init_dUtilde(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  double dtau;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
      fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in init_dUtilde() in ts.c");
    } else {
      find_constant_dtau(np,gl,l,&dtau);
      np[l].wk->dtau=dtau;
      for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]=-np[l].wk->Res[flux]*dtau;
    }
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }
}


void update_dU(np_t *np, gl_t *gl, long theta, long ls, long le){
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

    sweep_with_1D_segments(np,gl,zone,&init_dUtilde,SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
    sweep_with_1D_segments(np,gl,zone,&update_dUtilde,SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
    sweep_with_1D_segments(np,gl,zone,&update_dU,SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);

}



void find_dU(np_t *np, gl_t *gl, zone_t zone){

  sweep_with_1D_segments(np,gl,zone,&init_dUtilde,SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_dUtilde,SWEEPTYPE_IJK, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
}

