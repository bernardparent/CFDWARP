#include <cycle/res/_res.h>
#include <cycle/resconv/_resconv.h>
#include <cycle/restime/_restime.h>
#include <cycle/ressource/_ressource.h>
#include <cycle/share/res_share.h>
#include <cycle/share/cycle_share.h>

#define DUSTAR_SECONDORDER TRUE
#define YSTARDH_SECONDORDER TRUE
			

/* set TVDLIMITER to LIMITER_MINMOD for Roe (1986) minmod limiter
                     LIMITER_MINMOD2 for most compressive minmod type limiter
                     LIMITER_VANLEER for the Van Leer limiter
                     LIMITER_SUPERBEE for most compressive limiter possible
*/
#define TVDLIMITER LIMITER_VANLEER




      
void add_dKstar_dG_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  long l,vartheta,flux;
  flux_t Gp0[nd],Gm1[nd],dGm1h,tmpm1h;
  flux_t Gp0m1[nd],Gp0p1[nd],Gm1p1[nd],Gm1m1[nd];
  sqmat_t Km1h;
  metrics_t metricsm1h;

  for (l=ls; l!=_l_plus_one(_l_plus_one(le,gl,theta),gl,theta); l=_l_plus_one(l,gl,theta)){
    find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0),
                             theta, &metricsm1h);
    for (vartheta=0; vartheta<nd; vartheta++){
      find_Kstar_interface(np,gl,_al(gl,l,theta,-1),_al(gl,l,theta,+0),metricsm1h,theta,vartheta,Km1h);

      if (theta==vartheta) {
        if (l==ls) find_G(np[_al(gl,l,theta,-1)],gl,Gm1[vartheta]);
          else copy_vector(Gp0[vartheta],Gm1[vartheta]);
        find_G(np[_al(gl,l,theta,+0)],gl,Gp0[vartheta]);
        for (flux=0; flux<nf; flux++){
          dGm1h[flux]=Gp0[vartheta][flux]-Gm1[vartheta][flux];
        }
      } else {
        if (l==ls) {
          find_G(np[_all(gl,l,theta,-1,vartheta,+1)],gl,Gm1p1[vartheta]);
          find_G(np[_all(gl,l,theta,-1,vartheta,-1)],gl,Gm1m1[vartheta]);
        } else {
          copy_vector(Gp0p1[vartheta],Gm1p1[vartheta]);
          copy_vector(Gp0m1[vartheta],Gm1m1[vartheta]);
        }
        find_G(np[_all(gl,l,theta,+0,vartheta,+1)],gl,Gp0p1[vartheta]);
        find_G(np[_all(gl,l,theta,+0,vartheta,-1)],gl,Gp0m1[vartheta]);
        for (flux=0; flux<nf; flux++){
          dGm1h[flux]=0.25e0*(Gp0p1[vartheta][flux]+Gm1p1[vartheta][flux]
                             -Gp0m1[vartheta][flux]-Gm1m1[vartheta][flux]);
        }
      }
      multiply_matrix_and_vector(Km1h,dGm1h,tmpm1h);

#ifdef _RESCONV_INCLUDES_DIFFUSION
      for (flux=0; flux<nf; flux++) np[_al(gl,l,theta,-1)].wk->Fp1h_diffusion[theta][flux]=-tmpm1h[flux];
#else
      for (flux=0; flux<nf; flux++){
        if (l!=_l_plus_one(le,gl,theta)) np[l].wk->Res[flux]+=fact*tmpm1h[flux];
        if (l!=ls) np[_al(gl,l,theta,-1)].wk->Res[flux]-=fact*tmpm1h[flux];
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
        if (l!=_l_plus_one(le,gl,theta)) np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*tmpm1h[flux];
        if (l!=ls) np[_al(gl,l,theta,-1)].bs->Res_trapezoidal[flux]-=fact_trapezoidal*tmpm1h[flux];
#endif
      }
#endif
    }
  }

}





#ifdef _FLUID_PLASMA

void add_dDstarU_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  long l,flux;
  flux_t tmpm1h,DUstarplusm1,DUstarminusp0,DUstarplusp0,DUstarplusm2,DUstarminusm1,DUstarminusp1;
  flux_t Ustarm1,Ustarp0,Ustarm2,Ustarp1;
  sqmat_t Dstarplusm1,Dstarminusp0,Dstarplusp0,Dstarplusm2,Dstarminusm1,Dstarminusp1;

    /*  FVS-TVD discretization of d(Dstar*Ustar)/dX */
  for (l=ls; l!=_l_plus_one(_l_plus_one(le,gl,theta),gl,theta); l=_l_plus_one(l,gl,theta)){
    find_Dstarplus(np, gl, _al(gl,l,theta,-1), theta, Dstarplusm1);
    find_Dstarminus(np, gl, l, theta, Dstarminusp0);
    find_Ustar(np[_al(gl,l,theta,-1)], gl, Ustarm1);
    find_Ustar(np[_al(gl,l,theta,+0)], gl, Ustarp0);
    multiply_matrix_and_vector(Dstarplusm1,Ustarm1,DUstarplusm1);
    multiply_matrix_and_vector(Dstarminusp0,Ustarp0,DUstarminusp0);

    for (flux=0; flux<nf; flux++){
      tmpm1h[flux]=DUstarplusm1[flux]+DUstarminusp0[flux]; 
    }
    if (DUSTAR_SECONDORDER && !is_node_bdry(np[_al(gl,l,theta,-1)],TYPELEVEL_FLUID_WORK) && !is_node_bdry(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK)) {
      find_Dstarplus(np, gl, _al(gl,l,theta,+0), theta, Dstarplusp0);
      find_Dstarplus(np, gl, _al(gl,l,theta,-2), theta, Dstarplusm2);
      find_Dstarminus(np, gl, _al(gl,l,theta,-1), theta, Dstarminusm1);
      find_Dstarminus(np, gl, _al(gl,l,theta,+1), theta, Dstarminusp1);
      find_Ustar(np[_al(gl,l,theta,-2)], gl, Ustarm2);
      find_Ustar(np[_al(gl,l,theta,+1)], gl, Ustarp1);
      multiply_matrix_and_vector(Dstarplusp0,Ustarp0,DUstarplusp0);
      multiply_matrix_and_vector(Dstarplusm2,Ustarm2,DUstarplusm2);
      multiply_matrix_and_vector(Dstarminusm1,Ustarm1,DUstarminusm1);
      multiply_matrix_and_vector(Dstarminusp1,Ustarp1,DUstarminusp1);

      for (flux=0; flux<nf; flux++) {
        tmpm1h[flux]+=0.5*_limiter_TVD((DUstarplusp0[flux]-DUstarplusm1[flux])/notzero(DUstarplusm1[flux]-DUstarplusm2[flux],1e-99),TVDLIMITER)*(DUstarplusm1[flux]-DUstarplusm2[flux])
                     +0.5*_limiter_TVD((DUstarminusm1[flux]-DUstarminusp0[flux])/notzero(DUstarminusp0[flux]-DUstarminusp1[flux],1e-99),TVDLIMITER)*(DUstarminusp0[flux]-DUstarminusp1[flux]);
      }            

    }

    for (flux=0; flux<nf; flux++){
      if (l!=_l_plus_one(le,gl,theta)) np[l].wk->Res[flux]-=fact*tmpm1h[flux];
      if (l!=ls) np[_al(gl,l,theta,-1)].wk->Res[flux]+=fact*tmpm1h[flux];
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
      if (l!=_l_plus_one(le,gl,theta)) np[l].bs->Res_trapezoidal[flux]-=fact_trapezoidal*tmpm1h[flux];
      if (l!=ls) np[_al(gl,l,theta,-1)].bs->Res_trapezoidal[flux]+=fact_trapezoidal*tmpm1h[flux];
#endif
    }

  }


}



void add_Ystar_dH_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  long l,flux;
  flux_t tmpm1h,tmpp1h;
  flux_t Ystarm1h,Ystarp1h,Hp0,Hm1,Hp1,Hm2,Hp2;
  int cnt;


    /* first-order discretization of Ystar*dH/dX 
       note: the stencil is not conservative hence why tmpp1h is not equal to tmpm1h
    */ 
  for (cnt=1; cnt<=2; cnt++){
   for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    switch (cnt){
      case 1:
        find_Y1star_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, Ystarm1h);
        find_Y1star_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, Ystarp1h);
        find_H1(np, gl, _al(gl,l,theta,-1), Hm1);
        find_H1(np, gl, _al(gl,l,theta,+0), Hp0);
        find_H1(np, gl, _al(gl,l,theta,+1), Hp1);
      break;
      case 2:
        find_Y2star_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, Ystarm1h);
        find_Y2star_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, Ystarp1h);
        find_H2(np, gl, _al(gl,l,theta,-1), Hm1);
        find_H2(np, gl, _al(gl,l,theta,+0), Hp0);
        find_H2(np, gl, _al(gl,l,theta,+1), Hp1);
      break;
      default:
        fatal_error("cnt can not be set to %d in add_Ystar_dH_residual().",cnt);
    }
    for (flux=0; flux<nf; flux++){
      tmpm1h[flux]=0.5*(Ystarm1h[flux]+fabs(Ystarm1h[flux]))*(Hp0[flux]-Hm1[flux]);
      tmpp1h[flux]=0.5*(Ystarp1h[flux]-fabs(Ystarp1h[flux]))*(Hp1[flux]-Hp0[flux]);
    }
    if (YSTARDH_SECONDORDER && !is_node_bdry(np[_al(gl,l,theta,-1)],TYPELEVEL_FLUID_WORK) && !is_node_bdry(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK)) {
      switch (cnt){
        case 1: 
          find_H1(np, gl, _al(gl,l,theta,-2), Hm2);
        break;
        case 2: 
          find_H2(np, gl, _al(gl,l,theta,-2), Hm2);
        break;
        default:
          fatal_error("cnt can not be set to %d in add_Ystar_dH_residual().",cnt);
      }

      for (flux=0; flux<nf; flux++){
        tmpm1h[flux]=0.5*(Ystarm1h[flux]+fabs(Ystarm1h[flux]))*(_f_TVD2(Hm1[flux],Hp0[flux],Hp1[flux],TVDLIMITER)-_f_TVD2(Hm2[flux],Hm1[flux],Hp0[flux],TVDLIMITER));
      }
    }
    if (YSTARDH_SECONDORDER && !is_node_bdry(np[_al(gl,l,theta,+1)],TYPELEVEL_FLUID_WORK) && !is_node_bdry(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK)) {
      switch (cnt){
        case 1: 
          find_H1(np, gl, _al(gl,l,theta,+2), Hp2);
        break;
        case 2: 
          find_H2(np, gl, _al(gl,l,theta,+2), Hp2);
        break;
        default:
          fatal_error("cnt can not be set to %d in add_Ystar_dH_residual().",cnt);
      }


      for (flux=0; flux<nf; flux++){
        tmpp1h[flux]=0.5*(Ystarp1h[flux]-fabs(Ystarp1h[flux]))*(_f_TVD2(Hp2[flux],Hp1[flux],Hp0[flux],TVDLIMITER)-_f_TVD2(Hp1[flux],Hp0[flux],Hm1[flux],TVDLIMITER));

      }
    }

    for (flux=0; flux<nf; flux++){
      np[l].wk->Res[flux]+=fact*(+tmpm1h[flux]+tmpp1h[flux]);
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
      np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*(+tmpm1h[flux]+tmpp1h[flux]);
#endif
      
    }

   }
  }
}





#endif






void update_residual(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;
#if (defined(_RESTIME_CDF_TRAPEZOIDAL) || defined(_RESTIME_TRAPEZOIDAL))
  long flux;
#endif
  double fact,fact_trapezoidal;
//   wfprintf(stderr,"residual updated from i=%ld j=%ld to i=%ld j=%ld\n",np[ls].i,np[ls].j,np[le].i,np[le].j); 

#if (defined(_RESTIME_CDF_TRAPEZOIDAL) || defined(_RESTIME_TRAPEZOIDAL))
  fact=0.5;
  fact_trapezoidal=0.5;
#else
  fact=1.0;
  fact_trapezoidal=0.0;
#endif

  if (_FLUID_SOURCE) add_Sstar_residual(theta,ls,le,np,gl,fact,fact_trapezoidal);

  if (_FLUID_DIFFUSION) add_dKstar_dG_residual(theta,ls,le,np,gl,fact,fact_trapezoidal);

#ifdef _FLUID_PLASMA
  add_Ystar_dH_residual(theta,ls,le,np,gl,fact,fact_trapezoidal);
  add_dDstarU_residual(theta,ls,le,np,gl,fact,fact_trapezoidal);
#endif

#ifdef _RESTIME_CDF
  fact=1.0;
  fact_trapezoidal=0.0;
#endif
  if (_FLUID_CONVECTION) add_dFstar_residual(theta,ls,le,np,gl,fact,fact_trapezoidal);


#ifdef UNSTEADY 
  if (theta==nd-1) {
    fact=1.0;
    fact_trapezoidal=0.0;
    add_Z_dUstar_residual(theta, ls, le, np, gl,fact,fact_trapezoidal);
  }
#endif


#if (defined(_RESTIME_CDF_TRAPEZOIDAL) || defined(_RESTIME_TRAPEZOIDAL))
  if (theta==nd-1) {
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
      for (flux=0; flux<nf; flux++){
        np[l].wk->Res[flux]+=np[l].bs->Res_trapezoidal_m1[flux];
      }
    }
  }
#endif


  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
    np[l].wk->xi=_xi(np[l],gl,np[l].wk->Res);
    if (theta==0){
      thread_lock_global_set(gl, THREADTYPE_LOOP);
      gl->effiter_R+=1.0/(double)(gl->nn);
      thread_lock_global_unset(gl, THREADTYPE_LOOP);
    }
  }
}


void init_residual(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
#ifdef _RESCONV_DELTA_LAMBDA_STORAGE    
  long dim;
#endif
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]=0.0e0;
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
    for (flux=0; flux<nf; flux++) np[l].bs->Res_trapezoidal[flux]=0.0e0;
#endif
#ifdef _RESCONV_DELTA_LAMBDA_STORAGE    
    for (dim=0; dim<nd; dim++){
      for (flux=0; flux<nf; flux++) np[l].bs->Delta_Lambda[dim][flux]=0.0;
    }
#endif
  }
}


void unlock_residual(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
  }
}



void find_residual(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np, gl, zone, &init_residual, SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#ifdef _RESCONV_STORAGE_FSTAR
  zone_t zoneext;
  zoneext=_zone_intersection(gl->domain_all,_zone_expansion(zone,+0));  //can't expand zone when using MPI and periodic BCs
  init_Fstar_interfaces(np, gl, zoneext);
  sweep_with_1D_segments(np, gl, zoneext, &find_Fstar_interfaces, SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np, gl, zone, &update_residual, SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#ifdef _2DL
  init_Fstar_interfaces(np, gl, zoneext);
  sweep_with_1D_segments(np, gl, zoneext, &find_Fstar_interfaces, SWEEPTYPE_J, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np, gl, zone, &update_residual, SWEEPTYPE_J, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#endif
#ifdef _3DL
  init_Fstar_interfaces(np, gl, zoneext);
  sweep_with_1D_segments(np, gl, zoneext, &find_Fstar_interfaces, SWEEPTYPE_K, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np, gl, zone, &update_residual, SWEEPTYPE_K, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#endif
#else
  sweep_with_1D_segments(np, gl, zone, &update_residual, SWEEPTYPE_IJK, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
#endif
#ifdef ZONETHREADS
  sweep_with_1D_segments(np, gl, zone, &unlock_residual, SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#endif
}


void find_residual_per_dimension(np_t *np, gl_t *gl, zone_t zone, int sweeptype){
  sweep_with_1D_segments(np, gl, zone, &init_residual, SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np, gl, zone, &update_residual, sweeptype, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
#ifdef ZONETHREADS
  sweep_with_1D_segments(np, gl, zone, &unlock_residual, SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#endif
}

