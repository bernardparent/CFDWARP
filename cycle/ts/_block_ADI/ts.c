#include <cycle/share/cycle_share.h>
#include <cycle/ts/_ts.h>
#include <cycle/share/ts_share.h>
#include <model/_model.h>

#ifdef _CYCLE_PREDICTOR_CORRECTOR
#error the block ADI relaxation is not compatible with the predictor-corrector cycle
#endif

#define alpha 1.0e0

//#define FLUX_NORM

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


#ifdef FLUX_NORM
static void recondition_matrix(sqmat_t M){
  long row,col,flux;
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
         M[row][col]*=_flux_norm(col)/_flux_norm(row);
     /*   if (row<ncs){
          for (flux=ncs; flux<nf; flux++) M[row][flux]=0.0;
          for (flux=ncs; flux<nf; flux++) M[flux][row]=0.0;
        } */
     }
   }
}
#endif

/*
static void recondition_matrix_positive(sqmat_t M){
  long row,col,flux;
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
        M[row][col]=max(M[row][col],0.0);
     }
   }
}

static void recondition_matrix_negative(sqmat_t M){
  long row,col,flux;
   for (row=0; row<nf; row++){
     for (col=0; col<nf; col++){
        M[row][col]=min(M[row][col],0.0);
     }
   }
}

*/


void update_dUtilde(np_t *np, gl_t *gl, long theta, long ls, long le){
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


void init_dUtilde(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  flux_t dUtilde;
  double dtau;
  sqmat_t GammaInv;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) 
      find_Gamma_inverse(np,gl,l,GammaInv);
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    if (gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP) {
      for (flux=0; flux<nf; flux++) dUtilde[flux]=-np[l].wk->Res[flux];
      multiply_matrix_and_vector(GammaInv,dUtilde,np[l].wk->dUstar);
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

