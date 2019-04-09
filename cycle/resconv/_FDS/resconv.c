#include <cycle/resconv/_resconv.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <cycle/share/res_share.h>
#include <src/control.h>

#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif

#define limiter 2

#define ACCURACY_FIRSTORDER 1
#define ACCURACY_SECONDORDER 2

void write_disc_resconv_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    AVERAGING=AVERAGING_ROE;\n"
    "    ACCURACY=ACCURACY_SECONDORDER;\n"
    "    EIGENVALCOND=EIGENVALCOND_GNOFFO;\n"
    "  );\n"
  ,_RESCONV_ACTIONNAME,numiter_FDSplus_default);

}


void read_disc_resconv_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_disc_resconv_actions(char *actionname, char **argum, SOAP_codex_t *codex){

  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_RESCONV_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_RESCONV_ACTIONNAME);

    SOAP_add_int_to_vars(codex,"EIGENVALCOND_HARTEN",EIGENVALCOND_HARTEN); 
    SOAP_add_int_to_vars(codex,"EIGENVALCOND_PECLET",EIGENVALCOND_PECLET);
    SOAP_add_int_to_vars(codex,"EIGENVALCOND_GNOFFO",EIGENVALCOND_GNOFFO);
    SOAP_add_int_to_vars(codex,"AVERAGING_ROE",AVERAGING_ROE);
    SOAP_add_int_to_vars(codex,"AVERAGING_ARITH",AVERAGING_ARITH);
    SOAP_add_int_to_vars(codex,"ACCURACY_FIRSTORDER",ACCURACY_FIRSTORDER);
    SOAP_add_int_to_vars(codex,"ACCURACY_SECONDORDER",ACCURACY_SECONDORDER);

    gl->DISC_RESCONV_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_resconv_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"EIGENVALCOND",&gl->cycle.resconv.EIGENVALCOND);
    if (gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PECLET && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_HARTEN && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_GNOFFO)
      SOAP_fatal_error(codex,"EIGENVALCOND must be set to either EIGENVALCOND_PECLET, EIGENVALCOND_HARTEN, EIGENVALCOND_GNOFFO.");

    find_int_var_from_codex(codex,"AVERAGING",&gl->cycle.resconv.AVERAGING);
    if (gl->cycle.resconv.AVERAGING!=AVERAGING_ROE  && gl->cycle.resconv.AVERAGING!=AVERAGING_ARITH)
      SOAP_fatal_error(codex,"AVERAGING must be set to either AVERAGING_ROE or AVERAGING_ARITH.");

    find_int_var_from_codex(codex,"ACCURACY",&gl->cycle.resconv.ACCURACY);
    if (gl->cycle.resconv.ACCURACY!=ACCURACY_FIRSTORDER && gl->cycle.resconv.ACCURACY!=ACCURACY_SECONDORDER)
      SOAP_fatal_error(codex,"ACCURACY must be set to either ACCURACY_FIRSTORDER or ACCURACY_SECONDORDER.");

    gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FDS;

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}




static void find_jacvars_at_interface_local(np_t *np, gl_t *gl, metrics_t metrics, long lL, long lR, long theta, jacvars_t *jacvars){
  jacvars_t jacvarsL,jacvarsR;

  find_jacvars(np[lL], gl, metrics, theta, &jacvarsL);
  find_jacvars(np[lR], gl, metrics, theta, &jacvarsR);

  find_jacvars_at_interface_from_jacvars(jacvarsL, jacvarsR, gl, theta, metrics, gl->cycle.resconv.AVERAGING, jacvars);

}


static void find_alpha(jacvars_t jacvars, metrics_t metrics, np_t np_m, np_t np_p,
                       gl_t *gl, flux_t alphak){

  sqmat_t L;
  flux_t mattmp;
  long flux;
  flux_t Up,Um;

  find_L_from_jacvars(jacvars, metrics, L);
  find_Ustar_given_metrics(np_m, gl, metrics, Um);
  find_Ustar_given_metrics(np_p, gl, metrics, Up);
  for (flux=0; flux<nf; flux++) mattmp[flux]=Up[flux]-Um[flux];
  multiply_matrix_and_vector(L,mattmp,alphak);
}


static void find_g(flux_t alpham1,flux_t alphap0,flux_t alphap1,flux_t g){
  long flux;
  for (flux=0; flux<nf; flux++){
    g[flux]=0.0e0;
    if (limiter==1) {
      g[flux]=minmod(alpham1[flux],alphap0[flux])
             +minmod(alphap0[flux],alphap1[flux])
             -alphap0[flux];
    }
    if (limiter==2) {
      g[flux]=minmod3(alpham1[flux],alphap0[flux],alphap1[flux]);
    }
  }
}


static void find_Fstar_interface_2o(np_t *np, gl_t *gl, metrics_t metrics, long lm3h, long lm1h, long lp1h, long lp3h,
                     long theta, jacvars_t jacvarsm1, jacvars_t jacvarsp0, jacvars_t jacvarsp1,
                     flux_t Fint){
  flux_t Fm1h,Fp1h,alpham1,alphap0,alphap1,g,fluxtmp;
  sqmat_t R,lambdap;
  long flux;
  jacvars_t jacvarsm1h,jacvarsp1h;
#ifdef _RESTIME_CDF
  flux_t dFstarxt;
#endif


  find_alpha(jacvarsm1, metrics, np[lm3h], np[lm1h], gl, alpham1);
  find_alpha(jacvarsp0, metrics, np[lm1h], np[lp1h], gl, alphap0);
  find_alpha(jacvarsp1, metrics, np[lp1h], np[lp3h], gl, alphap1);
  find_g(alpham1,alphap0,alphap1,g);
  find_conditioned_Lambda_absolute_from_jacvars(jacvarsp0, metrics, gl->cycle.resconv.EIGENVALCOND, lambdap);
  find_Linv_from_jacvars(jacvarsp0, metrics, R);
  for (flux=0; flux<nf; flux++)
    fluxtmp[flux]=-lambdap[flux][flux]*alphap0[flux]
                  +lambdap[flux][flux]*g[flux];
  multiply_matrix_and_vector(R,fluxtmp,Fint);


  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_Fstar_from_jacvars(jacvarsm1h, metrics, Fm1h);
  find_Fstar_from_jacvars(jacvarsp1h, metrics, Fp1h);

  for (flux=0; flux<nf; flux++) Fint[flux]=0.5e0*(Fint[flux]+Fm1h[flux]+Fp1h[flux]);
#ifdef _RESTIME_CDF
  /* now add the spacetime cross-diffusion terms */
    find_Fstarxt_interface(np, gl, theta, lm1h, jacvarsp0, metrics, dFstarxt);
    for (flux=0; flux<nf; flux++) Fint[flux]+=dFstarxt[flux];
#endif

}


static void find_Fstar_interface_1o(np_t *np, gl_t *gl, metrics_t metrics, long lm1h, long lp1h,
                     long theta, jacvars_t jacvarsp0, flux_t Fint){
  flux_t Fm1h,Fp1h,alphap0,fluxtmp;
  sqmat_t R,lambdap;
  long flux;
  jacvars_t jacvarsm1h,jacvarsp1h;
#ifdef _RESTIME_CDF
  flux_t dFstarxt;
#endif


  find_alpha(jacvarsp0, metrics, np[lm1h], np[lp1h], gl, alphap0);
  find_conditioned_Lambda_absolute_from_jacvars(jacvarsp0, metrics, gl->cycle.resconv.EIGENVALCOND, lambdap);
  find_Linv_from_jacvars(jacvarsp0, metrics, R);
  for (flux=0; flux<nf; flux++)
    fluxtmp[flux]=-lambdap[flux][flux]*alphap0[flux];
  multiply_matrix_and_vector(R,fluxtmp,Fint);

  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_Fstar_from_jacvars(jacvarsm1h, metrics, Fm1h);
  find_Fstar_from_jacvars(jacvarsp1h, metrics, Fp1h);


  for (flux=0; flux<nf; flux++) Fint[flux]=0.5e0*(Fint[flux]+Fm1h[flux]+Fp1h[flux]);

#ifdef _RESTIME_CDF
  /* now add the spacetime cross-diffusion terms */
    find_Fstarxt_interface(np, gl, theta, lm1h, jacvarsp0, metrics, dFstarxt);
    for (flux=0; flux<nf; flux++) Fint[flux]+=dFstarxt[flux];
#endif
}



static void find_Delta_Lambda_for_dtau_local(np_t *np, gl_t *gl, long l, long theta, jacvars_t jacvarsm1h, jacvars_t jacvarsp1h, metrics_t metricsm1h, metrics_t metricsp1h, flux_t Delta_Lambda){
  sqmat_t Lambda_minus,Lambda_plus;
  long flux;
  find_Lambda_minus_dtau_FDS_from_jacvars(jacvarsm1h, metricsm1h, gl->cycle.resconv.EIGENVALCOND, Lambda_minus);
  find_Lambda_plus_dtau_FDS_from_jacvars(jacvarsp1h, metricsp1h, gl->cycle.resconv.EIGENVALCOND, Lambda_plus);

#ifdef _RESTIME_CDF
  sqmat_t Lambdaxt_minus,Lambdaxt_plus;
  long row,col;
  find_Lambdaxt_minus_dtau_from_jacvars(np, gl, _al(gl,l,theta,-1), l, theta, jacvarsm1h, metricsm1h, Lambdaxt_minus);
//  find_Lambdaxt_minus_dtau(np, gl, l, theta, Lambdaxt_minus);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Lambda_minus[row][col]-=Lambdaxt_minus[row][col];
    }
  }
  find_Lambdaxt_plus_dtau_from_jacvars(np, gl, l, _al(gl,l,theta,+1), theta, jacvarsp1h, metricsp1h, Lambdaxt_plus);
//  find_Lambdaxt_plus_dtau(np, gl, l, theta, Lambdaxt_plus);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Lambda_plus[row][col]-=Lambdaxt_plus[row][col];
    }
  }
#endif
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=Lambda_plus[flux][flux]-Lambda_minus[flux][flux];
}





void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  jacvars_t *jacvarstmp,*jacvarsm3h,*jacvarsm1h,*jacvarsp1h,*jacvarsp3h;
  flux_t Fm1h,Fp1h;
  long flux,l;
  metrics_t metricsm1h,metricsp1h;

  jacvarsm3h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsm1h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsp1h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsp3h= (jacvars_t *) malloc(sizeof(jacvars_t));

  l=_l_minus_one(ls,gl,theta);
  find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metricsp1h);
  if (is_node_bdry(np[l],TYPELEVEL_FLUID_WORK))
    find_jacvars_at_interface_local(np,gl,metricsp1h,l,l,theta,jacvarsm1h);
    else find_jacvars_at_interface_local(np,gl,metricsp1h,_al(gl,l,theta,-1),_al(gl,l,theta,+0),theta,jacvarsm1h);
  find_jacvars_at_interface_local(np,gl,metricsp1h,_al(gl,l,theta,+0),_al(gl,l,theta,+1),theta,jacvarsp1h);
  find_jacvars_at_interface_local(np,gl,metricsp1h,_al(gl,l,theta,+1),_al(gl,l,theta,+2),theta,jacvarsp3h);


    if (!is_node_bdry(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK) && gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER) {
      find_Fstar_interface_2o(np, gl,metricsp1h,_al(gl,l,theta,-1),_al(gl,l,theta,+0),_al(gl,l,theta,+1),_al(gl,l,theta,+2),
               theta, *jacvarsm1h,*jacvarsp1h,*jacvarsp3h,Fp1h);
    } else {
      find_Fstar_interface_1o(np, gl, metricsp1h,_al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, *jacvarsp1h, Fp1h);
    }


  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    jacvarstmp=jacvarsm3h;
    jacvarsm3h=jacvarsm1h;
    jacvarsm1h=jacvarsp1h;
    jacvarsp1h=jacvarsp3h;
    jacvarsp3h=jacvarstmp;
    metricsm1h=metricsp1h;
    find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metricsp1h);

    for (flux=0; flux<nf; flux++) Fm1h[flux]=Fp1h[flux];

    if (!is_node_bdry(np[_al(gl,l,theta,+1)],TYPELEVEL_FLUID_WORK))
    find_jacvars_at_interface_local(np,gl,metricsp1h,_al(gl,l,theta,+1),_al(gl,l,theta,+2),theta,jacvarsp3h);
    else find_jacvars_at_interface_local(np,gl,metricsp1h,_al(gl,l,theta,+1),_al(gl,l,theta,+1),theta,jacvarsp3h);

    if (!is_node_bdry(np[_al(gl,l,theta,+1)],TYPELEVEL_FLUID_WORK)  && gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER) {
      find_Fstar_interface_2o(np, gl, metricsp1h,_al(gl,l,theta,-1),_al(gl,l,theta,+0),_al(gl,l,theta,+1),_al(gl,l,theta,+2),
               theta, *jacvarsm1h,*jacvarsp1h,*jacvarsp3h,Fp1h);
    } else {
      find_Fstar_interface_1o(np, gl, metricsp1h,_al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, *jacvarsp1h, Fp1h);
    }

    for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=fact*(Fp1h[flux]-Fm1h[flux]);
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
    for (flux=0; flux<nf; flux++) np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*(Fp1h[flux]-Fm1h[flux]);
#endif
#ifdef _RESCONV_DELTA_LAMBDA_STORAGE
    find_Delta_Lambda_for_dtau_local(np, gl, l, theta, *jacvarsm1h, *jacvarsp1h, metricsm1h, metricsp1h, np[l].bs->Delta_Lambda[theta]); 
#endif

  }

  free(jacvarsm3h);
  free(jacvarsm1h);
  free(jacvarsp1h);
  free(jacvarsp3h);
}


#ifdef _RESCONV_DELTA_LAMBDA_STORAGE

void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  long flux;
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=np[l].bs->Delta_Lambda[theta][flux];
}

#else

void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  find_Delta_Lambda_for_dtau_local(np, gl, l, theta, Delta_Lambda);
}

#endif



