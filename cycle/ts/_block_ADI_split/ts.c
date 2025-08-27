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


static double _w_local(np_t np, long spec){
  double ret;
/*  switch (spec){
    case specO2: 
      ret=0.5;
    break;  
    case specCO2:
      ret=0.5;
    break;
    default:
      ret=1e-99;
  }*/
  ret=_w(np,spec);
  return(ret);
}

void set_matrices_in_block_TDMA_line_split(np_t *np, gl_t *gl, long theta, long l, 
                    sqmat_t M1, sqmat_t M2, sqmat_t M3, sqmat_t M4, long L,
                    double *AA, double *BB, double *CC, double *RHS){
  long row,col,spec;

  for (row=0; row<nfr; row++){
    for (col=0; col<nfr; col++){
      AA[EXM_mi(nfr,L,row,col)]=M1[row+fluxmass][col+fluxmass];
      BB[EXM_mi(nfr,L,row,col)]=M2[row+fluxmass][col+fluxmass];
      CC[EXM_mi(nfr,L,row,col)]=M3[row+fluxmass][col+fluxmass];
      RHS[EXM_mi(nfr,L,row,col)]=0.0;
    }
    RHS[EXM_mi(nfr,L,row,0)]=M4[row+fluxmass][0];
  }
  RHS[EXM_mi(nfr,L,0,0)]=0.0;
  for (spec=0; spec<ns; spec++) RHS[EXM_mi(nfr,L,0,0)]+=M4[spec][0];
  
  for (col=0; col<nfr; col++){
    AA[EXM_mi(nfr,L,0,col)]=0.0;
    BB[EXM_mi(nfr,L,0,col)]=0.0;
    CC[EXM_mi(nfr,L,0,col)]=0.0;
    for (row=0; row<ns; row++){
      AA[EXM_mi(nfr,L,0,col)]+=M1[row][col+fluxmass];
      BB[EXM_mi(nfr,L,0,col)]+=M2[row][col+fluxmass];
      CC[EXM_mi(nfr,L,0,col)]+=M3[row][col+fluxmass];
    }
  }  
  
  for (row=1; row<nfr; row++){
    AA[EXM_mi(nfr,L,row,0)]=0.0;
    BB[EXM_mi(nfr,L,row,0)]=0.0;
    CC[EXM_mi(nfr,L,row,0)]=0.0;
    for (col=0; col<ns; col++){
      AA[EXM_mi(nfr,L,row,0)]+=M1[row+fluxmass][col]*_w_local(np[l],col);
      BB[EXM_mi(nfr,L,row,0)]+=M2[row+fluxmass][col]*_w_local(np[l],col);
      CC[EXM_mi(nfr,L,row,0)]+=M3[row+fluxmass][col]*_w_local(np[l],col);
    }
  }
  
}


void set_matrices_in_block_TDMA_line(sqmat_t M1, sqmat_t M2, sqmat_t M3, sqmat_t M4, long L,
                    double *AA, double *BB, double *CC, double *RHS){
  long row,col;

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      AA[EXM_mi(nf,L,row,col)]=M1[row][col];
      BB[EXM_mi(nf,L,row,col)]=M2[row][col];
      CC[EXM_mi(nf,L,row,col)]=M3[row][col];
      RHS[EXM_mi(nf,L,row,col)]=M4[row][col];
    }
  }
}


void update_dUtilde_split(np_t *np, gl_t *gl, long theta, long ls, long le){
  long is,ie,js,je,ks,ke;
  long jj,l,flux,maxsize,row,col,spec;
  sqmat_t A1,B1,C1,D1,A2,B2,C2,garb,Mtmp,Am1h,Bm1h,Bp1h,Cp1h;
  double *AA, *BB, *CC, *RHS;
  double dtau,ressum;
  long maxnumb;
  sqmat_t dS_dU;

  set_matrix_to_zero(garb);
  find_ijk_from_l(gl, ls, &is, &js, &ks);
  find_ijk_from_l(gl, le, &ie, &je, &ke);
  if (theta==0){
    maxsize=ie-is+2;
  } else {
    if (theta==1){
      maxsize=je-js+2;
    } else {
      maxsize=ke-ks+2;
    }
  }
  maxnumb=((maxsize+1)*nfr*nfr);
  AA = (double *) malloc(maxnumb*sizeof(double));
  BB = (double *) malloc(maxnumb*sizeof(double));
  CC = (double *) malloc(maxnumb*sizeof(double));
  RHS = (double *) malloc(maxnumb*sizeof(double));

  /*-------Here: Boundary node at left */
  jj=0;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_minus_one(ls,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_plus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, C1);
  set_matrices_in_block_TDMA_line_split(np, gl, theta, l, A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  for(l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
   set_matrix_to_zero(A1);
   set_matrix_to_zero(B1);
   set_matrix_to_zero(C1);

   add_TDMA_jacobians_non_conservative(np, gl, theta, l, A1, B1, C1, FALSE);
   if (l==ls) {
     find_TDMA_jacobians_conservative(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);
   } else {
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         Am1h[row][col]=Bp1h[row][col];
         Bm1h[row][col]=Cp1h[row][col];
       }
     }
   }
   find_TDMA_jacobians_conservative(np, gl, theta, l, Bp1h, Cp1h);
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       A1[row][col]-=Am1h[row][col];
       B1[row][col]-=Bm1h[row][col];
       B1[row][col]+=Bp1h[row][col];
       C1[row][col]+=Cp1h[row][col];
     }
   }


   if (gl->PRECONDITIONER==PRECON_LOCALTIMESTEP) {
     dtau=np[l].wk->dtau;
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         Mtmp[row][col]=0.0e0;
       }
       Mtmp[row][row]=alpha*dtau;
     }
     multiply_diagonal_matrix_and_matrix(Mtmp,A1,A2);
     multiply_diagonal_matrix_and_matrix(Mtmp,B1,B2);
     multiply_diagonal_matrix_and_matrix(Mtmp,C1,C2);
   } else {
     fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in update_dUtilde() within ts.c.");     
   }
   for (flux=0; flux<nf; flux++) B2[flux][flux]+=1.0e0;
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       A1[row][col]=A2[row][col];
       B1[row][col]=B2[row][col];
       C1[row][col]=C2[row][col];
     }
   }
   set_matrix_to_zero(D1);
   for (flux=0; flux<nf; flux++) D1[flux][0]=np[l].wk->dUstar[flux];


   jj++;
   set_matrices_in_block_TDMA_line_split(np, gl, theta, l, A1, B1, C1, D1, jj, AA, BB, CC, RHS);
  }

  /*-------Here: Boundary node at right */
  jj++;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_plus_one(le,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_minus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, A1);
  set_matrices_in_block_TDMA_line_split(np, gl, theta, l, A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  /* solve the TDMA */
  EXM_solve_block_TDMA(AA, BB, CC, RHS, jj, nfr);

  /*--------Here: add RHS of TDMA to Ustar.*/
  jj=0;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    jj++;
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    //for (flux=0; flux<=fluxmass; flux++) np[l].wk->dUstar[flux]=0.0;
    for (flux=0; flux<nfr; flux++) np[l].wk->dUstar[flux+fluxmass]=RHS[EXM_mi(nfr,jj,flux,0)];
    //ressum=0.0;
    //for (spec=0; spec<ns; spec++) ressum+=np[l].wk->Res[flux];
    for (spec=0; spec<ns; spec++) np[l].wk->dUstar[spec]=_w_local(np[l],spec)*RHS[EXM_mi(nfr,jj,0,0)];
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }

  free(AA);
  free(BB);
  free(CC);
  free(RHS);
  

  EXM_mat_t matchem,matcheminv,matdU1,matdU2;
  double dUstarsum;
  spec_t deltadUstar;

  // here, invert the chemical jacobian to find the new mass fractions
  if (theta==nd-1 && FALSE){
    EXM_init_matrix(&matchem, ns, ns);
    EXM_init_matrix(&matcheminv, ns, ns);
    EXM_init_matrix(&matdU1, ns, 1);
    EXM_init_matrix(&matdU2, ns, 1);
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      thread_lock_node_set(np,l,THREADTYPE_ZONE);
      find_dSchem_dU(np, gl, l, dS_dU); 
      if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
        fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in update_dUtilde() within ts.c.");     
      }
      for (row=0; row<ns; row++){
        for (col=0; col<ns; col++){
          matchem.cont[EXM_aim(matchem.glm,row,col)]=-dS_dU[row][col]*alpha*np[l].wk->dtau; 
        }
        matchem.cont[EXM_aim(matchem.glm,row,row)]+=1.0e0;         
      }
      EXM_invert_matrix_partial_pivoting(matchem, &matcheminv);
      for (spec=0; spec<ns; spec++) matdU1.cont[EXM_aim(matdU1.glm,spec,0)]=np[l].wk->dUstar[spec];
      EXM_multiply_matrices(matcheminv, matdU1, &matdU2);

      for (spec=0; spec<ns; spec++) deltadUstar[spec]=matdU2.cont[EXM_aim(matdU2.glm,spec,0)]-np[l].wk->dUstar[spec];
      for (spec=0; spec<ns; spec++) np[l].wk->dUstar[spec]=matdU2.cont[EXM_aim(matdU2.glm,spec,0)];
      dUstarsum=0.0;
      for (spec=0; spec<ns; spec++) dUstarsum+=np[l].wk->dUstar[spec];
      for (spec=0; spec<ns; spec++) np[l].wk->dUstar[fluxet]+=deltadUstar[spec]*_hk_from_T(spec,298.0);
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
    EXM_free_matrix(&matchem);
    EXM_free_matrix(&matcheminv);
    EXM_free_matrix(&matdU1);
    EXM_free_matrix(&matdU2);
    
    
    
  }
  if (theta==nd-1 && FALSE){
    EXM_init_matrix(&matchem, ns+1, ns+1);
    EXM_init_matrix(&matcheminv, ns+1, ns+1);
    EXM_init_matrix(&matdU1, ns+1, 1);
    EXM_init_matrix(&matdU2, ns+1, 1);
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      thread_lock_node_set(np,l,THREADTYPE_ZONE);
      find_dSchem_dU(np, gl, l, dS_dU); 
      if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
        fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in update_dUtilde() within ts.c.");     
      }
      for (row=0; row<ns; row++){
        for (col=0; col<ns; col++){
          matchem.cont[EXM_aim(matchem.glm,row,col)]=-dS_dU[row][col]*alpha*np[l].wk->dtau; 
        }
        matchem.cont[EXM_aim(matchem.glm,row,ns)]=-dS_dU[row][fluxet]*alpha*np[l].wk->dtau; 
      }
      for (col=0; col<ns; col++){
        matchem.cont[EXM_aim(matchem.glm,ns,col)]=-dS_dU[fluxet][col]*alpha*np[l].wk->dtau; 
      }
      matchem.cont[EXM_aim(matchem.glm,ns,ns)]=-dS_dU[fluxet][fluxet]*alpha*np[l].wk->dtau;
      for (row=0; row<ns+1; row++) matchem.cont[EXM_aim(matchem.glm,row,row)]+=1.0e0;         

      EXM_invert_matrix_partial_pivoting(matchem, &matcheminv);
      for (spec=0; spec<ns; spec++) matdU1.cont[EXM_aim(matdU1.glm,spec,0)]=np[l].wk->dUstar[spec];
      matdU1.cont[EXM_aim(matdU1.glm,ns,0)]=np[l].wk->dUstar[fluxet];
      EXM_multiply_matrices(matcheminv, matdU1, &matdU2);
      for (spec=0; spec<ns; spec++) np[l].wk->dUstar[spec]=matdU2.cont[EXM_aim(matdU2.glm,spec,0)];
      np[l].wk->dUstar[fluxet]=matdU2.cont[EXM_aim(matdU2.glm,ns,0)];
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
    EXM_free_matrix(&matchem);
    EXM_free_matrix(&matcheminv);
    EXM_free_matrix(&matdU1);
    EXM_free_matrix(&matdU2);
  }
  if (theta==nd-1 && TRUE){
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
      for (flux=0; flux<nf; flux++) {
        if (flux<nf) matdU1.cont[EXM_aim(matdU1.glm,flux,0)]= np[l].wk->dUstar[flux]; 
        else matdU1.cont[EXM_aim(matdU1.glm,flux,0)]=np[l].wk->dUstar[flux];
      }
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




void update_dUtilde_chem(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux,row,col,spec,is,ie,js,je,ks,ke;
  sqmat_t dS_dU;

  find_ijk_from_l(gl, ls, &is, &js, &ks);
  find_ijk_from_l(gl, le, &ie, &je, &ke);

  EXM_mat_t matchem,matcheminv,matdU1,matdU2;
  double dUstarsum;
  spec_t deltadUstar;

  // here, invert the chemical jacobian to find the new mass fractions
  if (FALSE){
    EXM_init_matrix(&matchem, ns, ns);
    EXM_init_matrix(&matcheminv, ns, ns);
    EXM_init_matrix(&matdU1, ns, 1);
    EXM_init_matrix(&matdU2, ns, 1);
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      thread_lock_node_set(np,l,THREADTYPE_ZONE);
      find_dSchem_dU(np, gl, l, dS_dU); 
      if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
        fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in update_dUtilde() within ts.c.");     
      }
      for (row=0; row<ns; row++){
        for (col=0; col<ns; col++){
          matchem.cont[EXM_aim(matchem.glm,row,col)]=-dS_dU[row][col]*alpha*np[l].wk->dtau; 
        }
        matchem.cont[EXM_aim(matchem.glm,row,row)]+=1.0e0;         
      }
      EXM_invert_matrix_partial_pivoting(matchem, &matcheminv);
      for (spec=0; spec<ns; spec++) matdU1.cont[EXM_aim(matdU1.glm,spec,0)]=np[l].wk->dUstar[spec];
      EXM_multiply_matrices(matcheminv, matdU1, &matdU2);

      for (spec=0; spec<ns; spec++) deltadUstar[spec]=matdU2.cont[EXM_aim(matdU2.glm,spec,0)];
      for (spec=0; spec<ns; spec++) np[l].wk->dUstar[spec]=matdU2.cont[EXM_aim(matdU2.glm,spec,0)];
      dUstarsum=0.0;
      for (spec=0; spec<ns; spec++) dUstarsum+=np[l].wk->dUstar[spec];
      for (spec=0; spec<ns; spec++) np[l].wk->dUstar[fluxet]+=deltadUstar[spec]*_hk_from_T(spec,298.0);
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
    EXM_free_matrix(&matchem);
    EXM_free_matrix(&matcheminv);
    EXM_free_matrix(&matdU1);
    EXM_free_matrix(&matdU2);
    
    
    
  }
      double Omega;
  
  if (FALSE){
    EXM_init_matrix(&matchem, ns+1, ns+1);
    EXM_init_matrix(&matcheminv, ns+1, ns+1);
    EXM_init_matrix(&matdU1, ns+1, 1);
    EXM_init_matrix(&matdU2, ns+1, 1);
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      thread_lock_node_set(np,l,THREADTYPE_ZONE);
      if (l==ls || l==le){
        
      } else {
      find_dSchem_dU(np, gl, l, dS_dU); 
      if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
        fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in update_dUtilde() within ts.c.");     
      }
      for (row=0; row<ns; row++){
        for (col=0; col<ns; col++){
          matchem.cont[EXM_aim(matchem.glm,row,col)]=-dS_dU[row][col]*alpha*np[l].wk->dtau; 
        }
        matchem.cont[EXM_aim(matchem.glm,row,ns)]=-dS_dU[row][fluxet]*alpha*np[l].wk->dtau; 
      }
      for (col=0; col<ns; col++){
        matchem.cont[EXM_aim(matchem.glm,ns,col)]=-dS_dU[fluxet][col]*alpha*np[l].wk->dtau; 
      }
      matchem.cont[EXM_aim(matchem.glm,ns,ns)]=-dS_dU[fluxet][fluxet]*alpha*np[l].wk->dtau;
      for (row=0; row<ns+1; row++) matchem.cont[EXM_aim(matchem.glm,row,row)]+=1.0e0;         

      EXM_invert_matrix_partial_pivoting(matchem, &matcheminv);
      for (spec=0; spec<ns; spec++) matdU1.cont[EXM_aim(matdU1.glm,spec,0)]=np[l].wk->dUstar[spec];
      matdU1.cont[EXM_aim(matdU1.glm,ns,0)]=np[l].wk->dUstar[fluxet];
      EXM_multiply_matrices(matcheminv, matdU1, &matdU2);
      Omega=_Omega(np[l],gl);
      for (spec=0; spec<ns; spec++) {
        //np[l].bs->U[spec]+=matdU2.cont[EXM_aim(matdU2.glm,spec,0)]/Omega;
        //np[l].wk->dUstar[spec]=0.0;
        np[l].wk->dUstar[spec]=matdU2.cont[EXM_aim(matdU2.glm,spec,0)];
      }
//      np[l].bs->U[fluxet]+=matdU2.cont[EXM_aim(matdU2.glm,ns,0)]/Omega;
//      np[l].wk->dUstar[fluxet]=0.0;
        np[l].wk->dUstar[fluxet]=matdU2.cont[EXM_aim(matdU2.glm,ns,0)];
      }
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
    EXM_free_matrix(&matchem);
    EXM_free_matrix(&matcheminv);
    EXM_free_matrix(&matdU1);
    EXM_free_matrix(&matdU2);
  }
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
      Omega=_Omega(np[l],gl);
      for (flux=0; flux<nf; flux++) {
        np[l].bs->U[flux]+=matdU2.cont[EXM_aim(matdU2.glm,flux,0)]/Omega;
        np[l].wk->dUstar[flux]=0.0;
        //np[l].wk->dUstar[flux]=matdU2.cont[EXM_aim(matdU2.glm,flux,0)];
      }
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
    EXM_free_matrix(&matchem);
    EXM_free_matrix(&matcheminv);
    EXM_free_matrix(&matdU1);
    EXM_free_matrix(&matdU2);
  }
}



void update_dUtilde_nochem(np_t *np, gl_t *gl, long theta, long ls, long le){
  long is,ie,js,je,ks,ke;
  long jj,l,flux,maxsize,row,col,spec;
  sqmat_t A1,B1,C1,D1,A2,B2,C2,garb,Mtmp,Am1h,Bm1h,Bp1h,Cp1h;
  double *AA, *BB, *CC, *RHS;
  double dtau;
  long maxnumb;
  sqmat_t dS_dU;

  set_matrix_to_zero(garb);
  find_ijk_from_l(gl, ls, &is, &js, &ks);
  find_ijk_from_l(gl, le, &ie, &je, &ke);
  if (theta==0){
    maxsize=ie-is+2;
  } else {
    if (theta==1){
      maxsize=je-js+2;
    } else {
      maxsize=ke-ks+2;
    }
  }
  maxnumb=((maxsize+1)*nfr*nfr);
  AA = (double *) malloc(maxnumb*sizeof(double));
  BB = (double *) malloc(maxnumb*sizeof(double));
  CC = (double *) malloc(maxnumb*sizeof(double));
  RHS = (double *) malloc(maxnumb*sizeof(double));

  /*-------Here: Boundary node at left */
  jj=0;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_minus_one(ls,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_plus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, C1);
  set_matrices_in_block_TDMA_line_split(np, gl, theta, l, A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  for(l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
   set_matrix_to_zero(A1);
   set_matrix_to_zero(B1);
   set_matrix_to_zero(C1);


   add_TDMA_jacobians_non_conservative(np, gl, theta, l, A1, B1, C1, FALSE);
   if (l==ls) {
     find_TDMA_jacobians_conservative(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);
   } else {
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         Am1h[row][col]=Bp1h[row][col];
         Bm1h[row][col]=Cp1h[row][col];
       }
     }
   }
   find_TDMA_jacobians_conservative(np, gl, theta, l, Bp1h, Cp1h);
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       A1[row][col]-=Am1h[row][col];
       B1[row][col]-=Bm1h[row][col];
       B1[row][col]+=Bp1h[row][col];
       C1[row][col]+=Cp1h[row][col];
     }
   }


   if (gl->PRECONDITIONER==PRECON_LOCALTIMESTEP) {
     dtau=np[l].wk->dtau;
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         Mtmp[row][col]=0.0e0;
       }
       Mtmp[row][row]=alpha*dtau;
     }
     multiply_diagonal_matrix_and_matrix(Mtmp,A1,A2);
     multiply_diagonal_matrix_and_matrix(Mtmp,B1,B2);
     multiply_diagonal_matrix_and_matrix(Mtmp,C1,C2);
   } else {
     fatal_error("PRECONDITIONER must be set to PRECON_LOCALTIMESTEP in update_dUtilde() within ts.c.");     
   }
   for (flux=0; flux<nf; flux++) B2[flux][flux]+=1.0e0;
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       A1[row][col]=A2[row][col];
       B1[row][col]=B2[row][col];
       C1[row][col]=C2[row][col];
     }
   }
   set_matrix_to_zero(D1);
   for (flux=0; flux<nf; flux++) D1[flux][0]=np[l].wk->dUstar[flux];


   jj++;
   set_matrices_in_block_TDMA_line_split(np, gl, theta, l, A1, B1, C1, D1, jj, AA, BB, CC, RHS);
  }

  /*-------Here: Boundary node at right */
  jj++;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_plus_one(le,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_minus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, A1);
  set_matrices_in_block_TDMA_line_split(np, gl, theta, l, A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  /* solve the TDMA */
  EXM_solve_block_TDMA(AA, BB, CC, RHS, jj, nfr);

  /*--------Here: add RHS of TDMA to Ustar.*/
  jj=0;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    jj++;
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    //for (flux=0; flux<=fluxmass; flux++) np[l].wk->dUstar[flux]=0.0;
    for (flux=0; flux<nfr; flux++) np[l].wk->dUstar[flux+fluxmass]=RHS[EXM_mi(nfr,jj,flux,0)];
    for (spec=0; spec<ns; spec++) np[l].wk->dUstar[spec]=_w_local(np[l],spec)*RHS[EXM_mi(nfr,jj,0,0)];
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }

  free(AA);
  free(BB);
  free(CC);
  free(RHS);
  

}




void update_dUtilde_standard(np_t *np, gl_t *gl, long theta, long ls, long le){
  long is,ie,js,je,ks,ke;
  long jj,l,flux,maxsize,row,col;
  sqmat_t A1,B1,C1,D1,A2,B2,C2,GammaInv,garb,Mtmp,Am1h,Bm1h,Bp1h,Cp1h;
  double *AA, *BB, *CC, *RHS;
  double dtau;
  long maxnumb;

  set_matrix_to_zero(garb);
  find_ijk_from_l(gl, ls, &is, &js, &ks);
  find_ijk_from_l(gl, le, &ie, &je, &ke);
  if (theta==0){
    maxsize=ie-is+2;
  } else {
    if (theta==1){
      maxsize=je-js+2;
    } else {
      maxsize=ke-ks+2;
    }
  }
  maxnumb=((maxsize+1)*nf*nf);
  AA = (double *) malloc(maxnumb*sizeof(double));
  BB = (double *) malloc(maxnumb*sizeof(double));
  CC = (double *) malloc(maxnumb*sizeof(double));
  RHS = (double *) malloc(maxnumb*sizeof(double));

  /*-------Here: Boundary node at left */
  jj=0;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_minus_one(ls,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_plus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, C1);
  set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  for(l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
   set_matrix_to_zero(A1);
   set_matrix_to_zero(B1);
   set_matrix_to_zero(C1);


   add_TDMA_jacobians_non_conservative(np, gl, theta, l, A1, B1, C1, (theta==0));
   if (l==ls) {
     find_TDMA_jacobians_conservative(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);
   } else {
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         Am1h[row][col]=Bp1h[row][col];
         Bm1h[row][col]=Cp1h[row][col];
       }
     }
   }
   find_TDMA_jacobians_conservative(np, gl, theta, l, Bp1h, Cp1h);
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       A1[row][col]-=Am1h[row][col];
       B1[row][col]-=Bm1h[row][col];
       B1[row][col]+=Bp1h[row][col];
       C1[row][col]+=Cp1h[row][col];
     }
   }


   if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
     
     find_Gamma_inverse(np,gl,l,GammaInv);
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         GammaInv[row][col]*=alpha;
       }
     }
     multiply_matrix_and_matrix(GammaInv,A1,A2);
     multiply_matrix_and_matrix(GammaInv,B1,B2);
     multiply_matrix_and_matrix(GammaInv,C1,C2);
   } else {
     //find_constant_dtau(np,gl,l,&dtau);
     dtau=np[l].wk->dtau;
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         Mtmp[row][col]=0.0e0;
       }
       Mtmp[row][row]=alpha*dtau;
     }
     multiply_diagonal_matrix_and_matrix(Mtmp,A1,A2);
     multiply_diagonal_matrix_and_matrix(Mtmp,B1,B2);
     multiply_diagonal_matrix_and_matrix(Mtmp,C1,C2);
   }
   for (flux=0; flux<nf; flux++) B2[flux][flux]+=1.0e0;
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       A1[row][col]=A2[row][col];
       B1[row][col]=B2[row][col];
       C1[row][col]=C2[row][col];
     }
   }
   set_matrix_to_zero(D1);
#ifdef FLUX_NORM
   if (theta==0) for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]/=_flux_norm(flux);
#endif
   for (flux=0; flux<nf; flux++) D1[flux][0]=np[l].wk->dUstar[flux];


#ifdef FLUX_NORM
   recondition_matrix(A1);
   recondition_matrix(B1);
   recondition_matrix(C1);
#endif

//   recondition_matrix_negative(A1);
//   recondition_matrix_positive(B1);
//   recondition_matrix_negative(C1);



   jj++;
   set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);
  }

  /*-------Here: Boundary node at right */
  jj++;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_plus_one(le,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_minus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, A1);
  set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  /* solve the TDMA */
  EXM_solve_block_TDMA(AA, BB, CC, RHS, jj, nf);

  /*--------Here: add RHS of TDMA to Ustar.*/
  jj=0;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    jj++;
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]=RHS[EXM_mi(nf,jj,flux,0)];

#ifdef FLUX_NORM
    if (theta==nd-1) {
      for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]*=_flux_norm(flux);
    }
#endif

    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }

  free(AA);
  free(BB);
  free(CC);
  free(RHS);
}



void update_dUtilde_chem2(np_t *np, gl_t *gl, long theta, long ls, long le){
  long is,ie,js,je,ks,ke;
  long jj,l,flux,maxsize,row,col;
  sqmat_t dSstar_dUstar,A1,B1,C1,D1,A2,B2,C2,garb,Mtmp;
  double *AA, *BB, *CC, *RHS;
  double dtau;
  long maxnumb;

  set_matrix_to_zero(garb);
  find_ijk_from_l(gl, ls, &is, &js, &ks);
  find_ijk_from_l(gl, le, &ie, &je, &ke);
  if (theta==0){
    maxsize=ie-is+2;
  } else {
    if (theta==1){
      maxsize=je-js+2;
    } else {
      maxsize=ke-ks+2;
    }
  }
  maxnumb=((maxsize+1)*nf*nf);
  AA = (double *) malloc(maxnumb*sizeof(double));
  BB = (double *) malloc(maxnumb*sizeof(double));
  CC = (double *) malloc(maxnumb*sizeof(double));
  RHS = (double *) malloc(maxnumb*sizeof(double));

  /*-------Here: Boundary node at left */
  jj=0;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_minus_one(ls,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_plus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, C1);
  set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  for(l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
   set_matrix_to_zero(A1);
   set_matrix_to_zero(B1);
   set_matrix_to_zero(C1);


   find_dSstar_dUstar(np,gl,l,dSstar_dUstar);
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       B1[row][col]+=-dSstar_dUstar[row][col];
     }
   }

   if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
     fatal_error("Problem here..");     
   } else {
     dtau=np[l].wk->dtau;
     for (row=0; row<nf; row++){
       for (col=0; col<nf; col++){
         Mtmp[row][col]=0.0e0;
       }
       Mtmp[row][row]=alpha*dtau;
     }
     multiply_diagonal_matrix_and_matrix(Mtmp,A1,A2);
     multiply_diagonal_matrix_and_matrix(Mtmp,B1,B2);
     multiply_diagonal_matrix_and_matrix(Mtmp,C1,C2);
   }
   for (flux=0; flux<nf; flux++) B2[flux][flux]+=1.0e0;
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       A1[row][col]=A2[row][col];
       B1[row][col]=B2[row][col];
       C1[row][col]=C2[row][col];
     }
   }
   set_matrix_to_zero(D1);
   for (flux=0; flux<nf; flux++) D1[flux][0]=np[l].wk->dUstar[flux];




   jj++;
   set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);
  }

  /*-------Here: Boundary node at right */
  jj++;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_plus_one(le,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_minus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, A1);
  set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  /* solve the TDMA */
  EXM_solve_block_TDMA(AA, BB, CC, RHS, jj, nf);

  /*--------Here: add RHS of TDMA to Ustar.*/
  jj=0;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    jj++;
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]=RHS[EXM_mi(nf,jj,flux,0)];


    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }

  free(AA);
  free(BB);
  free(CC);
  free(RHS);
}


void update_dUtilde_chem3(np_t *np, gl_t *gl, long theta, long ls, long le){
  long is,ie,js,je,ks,ke;
  long jj,l,flux,maxsize,row,col;
  sqmat_t dSstar_dUstar,A1,B1,C1,D1,garb,Mtmp;
  double *AA, *BB, *CC, *RHS;
  double dtau;
  long maxnumb;

  set_matrix_to_zero(garb);
  find_ijk_from_l(gl, ls, &is, &js, &ks);
  find_ijk_from_l(gl, le, &ie, &je, &ke);
  if (theta==0){
    maxsize=ie-is+2;
  } else {
    if (theta==1){
      maxsize=je-js+2;
    } else {
      maxsize=ke-ks+2;
    }
  }
  maxnumb=((maxsize+1)*nf*nf);
  AA = (double *) malloc(maxnumb*sizeof(double));
  BB = (double *) malloc(maxnumb*sizeof(double));
  CC = (double *) malloc(maxnumb*sizeof(double));
  RHS = (double *) malloc(maxnumb*sizeof(double));

  /*-------Here: Boundary node at left */
  jj=0;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_minus_one(ls,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_plus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, C1);
  set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  for(l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
   set_matrix_to_zero(A1);
   set_matrix_to_zero(B1);
   set_matrix_to_zero(C1);


   find_dSstar_dUstar(np,gl,l,dSstar_dUstar);
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
       B1[row][col]+=-dSstar_dUstar[row][col]*np[l].wk->dtau*alpha;
     }
   }

   for (flux=0; flux<nf; flux++) B1[flux][flux]+=1.0e0;
   set_matrix_to_zero(D1);
   for (flux=0; flux<nf; flux++) D1[flux][0]=np[l].wk->dUstar[flux];




   jj++;
   set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);
  }

  /*-------Here: Boundary node at right */
  jj++;
  set_matrix_to_zero(A1);
  set_matrix_to_identity(B1);
  set_matrix_to_zero(C1);
  set_matrix_to_zero(D1);
  l=_l_plus_one(le,gl,theta);
  find_bdry_jacobian(np, gl, _node_type(np[l],TYPELEVEL_FLUID_WORK), l, _l_minus_one(l,gl,theta), TYPELEVEL_FLUID_WORK, B1, A1);
  set_matrices_in_block_TDMA_line(A1, B1, C1, D1, jj, AA, BB, CC, RHS);

  /* solve the TDMA */
  EXM_solve_block_TDMA(AA, BB, CC, RHS, jj, nf);

  /*--------Here: add RHS of TDMA to Ustar.*/
  jj=0;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    jj++;
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    for (flux=0; flux<nf; flux++) np[l].wk->dUstar[flux]=RHS[EXM_mi(nf,jj,flux,0)];


    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }

  free(AA);
  free(BB);
  free(CC);
  free(RHS);
}



void update_dUtilde(np_t *np, gl_t *gl, long theta, long ls, long le){
//  if (theta==nd-1 || TRUE) update_dUtilde_split(np,gl,theta,ls,le);
//    else update_dUtilde_split(np,gl,theta,ls,le);
//  if (theta==0) update_dUtilde_chem(np,gl,theta,ls,le);
  if (theta==0) update_dUtilde_chem(np,gl,theta,ls,le);
  update_dUtilde_nochem(np,gl,theta,ls,le);
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
    sweep_with_1D_segments(np,gl,zone,&update_dUtilde,SWEEPTYPE_IJK, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
    sweep_with_1D_segments(np,gl,zone,&update_dU,SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);

}



void find_dU(np_t *np, gl_t *gl, zone_t zone){

  sweep_with_1D_segments(np,gl,zone,&init_dUtilde,SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_dUtilde,SWEEPTYPE_IJK, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
}

