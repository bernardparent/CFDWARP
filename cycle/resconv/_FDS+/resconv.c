#include <cycle/resconv/_resconv.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <cycle/share/res_share.h>
#include <src/control.h>

#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif


#define ACCURACY_FIRSTORDER 1
#define ACCURACY_SECONDORDER 2

/* set the following to true to condition the Roe scheme eigenvalues (prior to the positivity-preserving modifications) */
#define ROE_EIGENVALUE_CONDITIONING TRUE

/* set the following to true to condition the eigenvalues in FVS form after the positivity-preserving modifications */
#define EIGENVALUE_CONDITIONING FALSE

/* set the following to true to get the Roe scheme eigenvalues */
#define ROE_EIGENVALUES FALSE




void write_disc_resconv_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    AVERAGING=AVERAGING_ROE;\n"
    "    ACCURACY=ACCURACY_SECONDORDER;\n"
    "    numiter=%d;\n"
    "    xi=1.99;\n"
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
      SOAP_fatal_error(codex,"AVERAGING must be set to either AVERAGING_ROE, or AVERAGING_ARITH.");

    find_int_var_from_codex(codex,"ACCURACY",&gl->cycle.resconv.ACCURACY);
    if (gl->cycle.resconv.ACCURACY!=ACCURACY_FIRSTORDER && gl->cycle.resconv.ACCURACY!=ACCURACY_SECONDORDER)
      SOAP_fatal_error(codex,"ACCURACY must be set to either ACCURACY_FIRSTORDER or ACCURACY_SECONDORDER.");

    find_double_var_from_codex(codex,"xi",&gl->cycle.resconv.xi);
    if (gl->cycle.resconv.xi<0.0 || gl->cycle.resconv.xi>2.0)
      SOAP_fatal_error(codex,"xi must be set to a value between 0 and 2.");

    find_int_var_from_codex(codex,"numiter",&gl->cycle.resconv.numiter);
    if (gl->cycle.resconv.numiter<0)
      SOAP_fatal_error(codex,"numiter must be set to a positive integer.");

    gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FDS;

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}


static void find_Lambda_minus_plus_Roeplus(np_t *np, gl_t *gl, long lp0, long lp1, jacvars_t jacvarsp1h, metrics_t metrics, sqmat_t lambdaminus, sqmat_t lambdaplus, long theta){
  sqmat_t Yminus,Yplus,Zplus,Zminus,lambdap0,lambdap1,Lp0,Linvp0,Linvp1,Lp1,Lp1h,Linvp1h,lambdapp1h;
  long row,col,flux,cnt;
  flux_t Up0,fluxtmp,fluxtmp2,fluxtmp3,LUp0,LUp1,Up1;
  jacvars_t jacvarsp0,jacvarsp1;
  np_t npp0,npp1;
#ifdef _RESTIME_CDF
  sqmat_t lambdaxtp1h,lambdaxtp0,lambdaxtp1;
  flux_t Deltaxt_p1h;
#endif
  npp0=np[lp0];
  npp1=np[lp1];

  find_jacvars(npp0,gl,metrics,theta,&jacvarsp0);
  find_jacvars(npp1,gl,metrics,theta,&jacvarsp1);


  find_Lambda_from_jacvars(jacvarsp0, metrics, lambdap0);
  find_Lambda_from_jacvars(jacvarsp1, metrics, lambdap1);
#ifdef _RESTIME_CDF
  find_Lambdaxt_interface(np, gl, lp0, lp1, theta, jacvarsp1h, metrics, lambdaxtp1h);
  find_Lambdaxt(np, gl, lp0, jacvarsp0, metrics, lambdaxtp0);
  find_Lambdaxt(np, gl, lp1, jacvarsp1, metrics, lambdaxtp1);
  find_Deltaxt_interface(np, gl, lp0, lp1, theta, jacvarsp0, jacvarsp1, metrics, Deltaxt_p1h);
  for (flux=0; flux<nf; flux++) lambdap0[flux][flux]+=-lambdaxtp0[flux][flux]*Deltaxt_p1h[flux];
  for (flux=0; flux<nf; flux++) lambdap1[flux][flux]+=-lambdaxtp1[flux][flux]*Deltaxt_p1h[flux];
#endif


  find_L_from_jacvars(jacvarsp1h, metrics, Lp1h);
  find_Linv_from_jacvars(jacvarsp1h, metrics, Linvp1h);
  if (ROE_EIGENVALUE_CONDITIONING){
    find_conditioned_Lambda_absolute_from_jacvars(jacvarsp1h, metrics, gl->cycle.resconv.EIGENVALCOND, lambdapp1h);
  } else {
    find_Lambda_from_jacvars(jacvarsp1h, metrics, lambdapp1h);
    for (flux=0; flux<nf; flux++) lambdapp1h[flux][flux]=fabs(lambdapp1h[flux][flux]);
  }

#ifdef _RESTIME_CDF
  for (flux=0; flux<nf; flux++) lambdapp1h[flux][flux]+=-fabs(lambdaxtp1h[flux][flux])*Deltaxt_p1h[flux];
#endif

  find_L_from_jacvars(jacvarsp0, metrics, Lp0);
  find_Ustar_given_metrics(npp0,gl,metrics,Up0);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0);

  find_L_from_jacvars(jacvarsp1, metrics, Lp1);
  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_Linv_from_jacvars(jacvarsp1, metrics, Linvp1);
  find_Ustar_given_metrics(npp1,gl,metrics,Up1);
  find_LUstar_from_jacvars(jacvarsp1, metrics, LUp1);


  multiply_matrix_and_vector(Lp1h,Up0,fluxtmp);
  multiply_matrix_and_vector(lambdapp1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(Linvp1h,fluxtmp2,fluxtmp);
  multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);


  multiply_matrix_and_vector(Lp1h,Up1,fluxtmp);
  multiply_matrix_and_vector(lambdapp1h,fluxtmp,fluxtmp3);
  multiply_matrix_and_vector(Linvp1h,fluxtmp3,fluxtmp);
  multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaminus[row][col]=0.0;
      lambdaplus[row][col]=0.0;
      Yminus[row][col]=0.0;
      Yplus[row][col]=0.0;
      Zminus[row][col]=0.0;
      Zplus[row][col]=0.0;
    }
  }

    for (flux=0; flux<nf; flux++) {
      Yminus[flux][flux]=min(0.0,
        0.5*lambdap0[flux][flux]+0.5*fluxtmp2[flux]/notzero(LUp0[flux],1.0e-30) );
      Yplus[flux][flux]=max(0.0,
        0.5*lambdap0[flux][flux]+0.5*fluxtmp2[flux]/notzero(LUp0[flux],1.0e-30) );
      Zplus[flux][flux]=max(0.0,
        0.5*lambdap1[flux][flux]-0.5*fluxtmp3[flux]/notzero(LUp1[flux],1.0e-30) );
      Zminus[flux][flux]=min(0.0,
        0.5*lambdap1[flux][flux]-0.5*fluxtmp3[flux]/notzero(LUp1[flux],1.0e-30) );
    }


  
  for (cnt=0; cnt<gl->cycle.resconv.numiter; cnt++){

    multiply_matrix_and_vector(Zplus,LUp1,fluxtmp2);
    multiply_matrix_and_vector(Linvp1,fluxtmp2,fluxtmp);
    multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);

    multiply_matrix_and_vector(Yminus,LUp0,fluxtmp3);
    multiply_matrix_and_vector(Linvp0,fluxtmp3,fluxtmp);
    multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);
    for (flux=0; flux<nf; flux++) {
      Yminus[flux][flux]=min(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1.0e-30));
      Zplus[flux][flux]=max(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1.0e-30));
    }
    for (flux=0; flux<nf; flux++) {
      Yplus[flux][flux]=max(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1.0e-30));
      Zminus[flux][flux]=min(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1.0e-30));
    }
    
  }

  
  for (flux=0; flux<nf; flux++) {
    lambdaplus[flux][flux]=Yplus[flux][flux]+Zplus[flux][flux];
    lambdaminus[flux][flux]=Yminus[flux][flux]+Zminus[flux][flux];
  }
  if (EIGENVALUE_CONDITIONING) condition_Lambda_plus_minus(np, gl, lp0, theta, jacvarsp0, jacvarsp1, metrics, gl->cycle.resconv.EIGENVALCOND, lambdaplus, lambdaminus);
  
}


/* the following gives exactly the Roe scheme but with an arithmetic average at the interface */
static void find_Lambda_minus_plus_Roe(np_t *np, gl_t *gl, long lp0, long lp1, jacvars_t jacvarsp1h, metrics_t metrics, sqmat_t lambdaminus, sqmat_t lambdaplus, long theta){
  sqmat_t lambdap0,lambdap1,Lp0,Linvp0,Linvp1,Lp1,Lp1h,Linvp1h,lambdapp1h;
  long row,col,flux;
  flux_t Up0,fluxtmp,fluxtmp2,fluxtmp3,LUp0,LQp1,Up1;
  jacvars_t jacvarsp0,jacvarsp1;
  np_t npp0,npp1;

  npp0=np[lp0];
  npp1=np[lp1];

  find_jacvars(npp0,gl,metrics,theta,&jacvarsp0);
  find_jacvars(npp1,gl,metrics,theta,&jacvarsp1);


  find_Lambda_from_jacvars(jacvarsp0, metrics, lambdap0);
  find_Lambda_from_jacvars(jacvarsp1, metrics, lambdap1);


  find_L_from_jacvars(jacvarsp1h, metrics, Lp1h);
  find_Linv_from_jacvars(jacvarsp1h, metrics, Linvp1h);
  if (ROE_EIGENVALUE_CONDITIONING){
    find_conditioned_Lambda_absolute_from_jacvars(jacvarsp1h, metrics, gl->cycle.resconv.EIGENVALCOND, lambdapp1h);
  } else {
    find_Lambda_from_jacvars(jacvarsp1h, metrics, lambdapp1h);
    for (flux=0; flux<nf; flux++) lambdapp1h[flux][flux]=fabs(lambdapp1h[flux][flux]);
  }


  find_L_from_jacvars(jacvarsp0, metrics, Lp0);
  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_Ustar_given_metrics(npp0,gl,metrics,Up0);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0);

  find_L_from_jacvars(jacvarsp1, metrics, Lp1);
  find_Linv_from_jacvars(jacvarsp1, metrics, Linvp1);
  find_Ustar_given_metrics(npp1,gl,metrics,Up1);
  find_LUstar_from_jacvars(jacvarsp1, metrics, LQp1);


  multiply_matrix_and_vector(Lp1h,Up0,fluxtmp);
  multiply_matrix_and_vector(lambdapp1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(Linvp1h,fluxtmp2,fluxtmp);
  multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);


  multiply_matrix_and_vector(Lp1h,Up1,fluxtmp);
  multiply_matrix_and_vector(lambdapp1h,fluxtmp,fluxtmp3);
  multiply_matrix_and_vector(Linvp1h,fluxtmp3,fluxtmp);
  multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaminus[row][col]=0.0;
      lambdaplus[row][col]=0.0;
    }
  }


  for (flux=0; flux<nf; flux++) {
    lambdaplus[flux][flux]=0.5*lambdap0[flux][flux]+0.5*fluxtmp2[flux]/notzero(LUp0[flux],1.0e-30);
    lambdaminus[flux][flux]=0.5*lambdap1[flux][flux]-0.5*fluxtmp3[flux]/notzero(LQp1[flux],1.0e-30);
  }


}


static void find_Lambda_minus_plus(np_t *np, gl_t *gl, long lp0, long lp1, jacvars_t jacvarsp1h, metrics_t metrics, sqmat_t lambdaminus, sqmat_t lambdaplus, long theta){
   if (ROE_EIGENVALUES){
     find_Lambda_minus_plus_Roe(np,gl,lp0,lp1,jacvarsp1h, metrics, lambdaminus, lambdaplus, theta);
   } else {
     find_Lambda_minus_plus_Roeplus(np,gl,lp0,lp1,jacvarsp1h, metrics, lambdaminus, lambdaplus, theta);
   }
}



static void find_Fstar_interface_1o(np_t *np, gl_t *gl, long lm1h, long lp1h,
                     long theta, jacvars_t jacvarsp0, flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h){
  flux_t Fm1h,Fp1h,fluxtmp,fluxtmp2;
  /* flux_t Ustar; */
  sqmat_t R;
  /* sqmat_t L; */
  long flux;
  metrics_t metrics;
  jacvars_t jacvarsm1h,jacvarsp1h;


  find_metrics_at_interface(np, gl, lm1h, lp1h, theta, &metrics);
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);

  find_Lambda_minus_plus(np,gl,lm1h,lp1h,jacvarsp0, metrics, lambdaminusp1h, lambdaplusm1h, theta);

  for (flux=0; flux<nf; flux++){
    np[lm1h].bs->Delta_Lambda[theta][flux]+=lambdaplusm1h[flux][flux];
    np[lp1h].bs->Delta_Lambda[theta][flux]-=lambdaminusp1h[flux][flux];
  }

  find_Linv_from_jacvars(jacvarsm1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsm1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaplusm1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fm1h);


  find_Linv_from_jacvars(jacvarsp1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsp1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaminusp1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fp1h);


  for (flux=0; flux<nf; flux++) Fint[flux]=Fm1h[flux]+Fp1h[flux];

}



static void find_M(jacvars_t jacvars_p0, metrics_t metrics, np_t np_m1h, np_t np_p1h,
                       gl_t *gl, flux_t Mp0){

  sqmat_t L;
  flux_t mattmp;
  long flux;
  flux_t Up1h,Um1h;

  find_L_from_jacvars(jacvars_p0, metrics, L);
  find_Ustar_given_metrics(np_m1h, gl, metrics, Um1h);
  find_Ustar_given_metrics(np_p1h, gl, metrics, Up1h);
  for (flux=0; flux<nf; flux++) mattmp[flux]=Up1h[flux]-Um1h[flux];
  multiply_matrix_and_vector(L,mattmp,Mp0);
}




static void find_Gminus(np_t *np, long lm1, long lp0, jacvars_t jacvarsp0, jacvars_t jacvarsm1h, gl_t *gl, metrics_t metrics, long theta, flux_t Gminusp0){
  sqmat_t lambdaminusp0,mattmp;
  flux_t LUp0;


  find_Lambda_minus_plus(np, gl, lm1, lp0, jacvarsm1h, metrics, lambdaminusp0, mattmp, theta);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0); 

  multiply_diagonal_matrix_and_vector(lambdaminusp0,LUp0,Gminusp0);
  
}


static void find_Gplus(np_t *np, long lp0, long lp1, jacvars_t jacvarsp0, jacvars_t jacvarsp1h, gl_t *gl, metrics_t metrics, long theta, flux_t Gplusp0){
  sqmat_t lambdaplusp0,mattmp;
  flux_t LUp0;


  find_Lambda_minus_plus(np, gl, lp0, lp1, jacvarsp1h, metrics, mattmp, lambdaplusp0, theta);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0); 
  multiply_diagonal_matrix_and_vector(lambdaplusp0,LUp0,Gplusp0);
  
}




static void find_Fstar_interface_2o(np_t *np, gl_t *gl, long lm3h, long lm1h, long lp1h, long lp3h,
                     long theta, jacvars_t jacvarsm1, jacvars_t jacvarsp0, jacvars_t jacvarsp1,
                     flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h){
  flux_t Mm1,Mp0,Mp1,fluxtmp,fluxtmp1,fluxtmp2,fluxtmp3,fluxtmp4,
         Gminusm1h,Gminusp1h,Gplusm1h,Gplusp1h,Psiminusp0,Psiplusp0;
  sqmat_t Linvp0,lambdapp0,lambdap0,Linvm1h,Linvp1h,Lp1h,Lm1h;
  jacvars_t jacvarsp1h,jacvarsm1h;
  long flux;
  metrics_t metrics;



  find_metrics_at_interface(np, gl, lm1h, lp1h, theta, &metrics);

  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);

  find_Lambda_minus_plus(np,gl,lm1h,lp1h,jacvarsp0, metrics, lambdaminusp1h, lambdaplusm1h, theta);



  find_M(jacvarsm1, metrics, np[lm3h], np[lm1h], gl, Mm1);
  find_M(jacvarsp0, metrics, np[lm1h], np[lp1h], gl, Mp0);
  find_M(jacvarsp1, metrics, np[lp1h], np[lp3h], gl, Mp1);
  find_Lambda_from_jacvars(jacvarsp0, metrics, lambdap0);
  if (ROE_EIGENVALUE_CONDITIONING){
    find_conditioned_Lambda_absolute_from_jacvars(jacvarsp0, metrics, gl->cycle.resconv.EIGENVALCOND, lambdapp0);
  } else {
    for (flux=0; flux<nf; flux++) lambdapp0[flux][flux]=fabs(lambdap0[flux][flux]);
  }
  find_Gminus(np,lm3h,lm1h, jacvarsm1h, jacvarsm1, gl, metrics, theta, Gminusm1h);
  find_Gminus(np,lm1h,lp1h, jacvarsp1h, jacvarsp0, gl, metrics, theta, Gminusp1h);
  find_Gplus(np,lm1h,lp1h, jacvarsm1h, jacvarsp0, gl, metrics, theta, Gplusm1h);
  find_Gplus(np,lp1h,lp3h, jacvarsp1h, jacvarsp1, gl, metrics, theta, Gplusp1h);


  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_L_from_jacvars(jacvarsp1h, metrics, Lp1h);
  find_L_from_jacvars(jacvarsm1h, metrics, Lm1h);
  find_Linv_from_jacvars(jacvarsp1h, metrics, Linvp1h);
  find_Linv_from_jacvars(jacvarsm1h, metrics, Linvm1h);


  /* here find Psiminusp0 */
  for (flux=0; flux<nf; flux++)
    fluxtmp[flux]=(lambdap0[flux][flux]-lambdapp0[flux][flux])*minmod3(Mm1[flux],Mp0[flux],Mp1[flux]);
  multiply_matrix_and_vector(Linvp0,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(Lp1h,fluxtmp2,fluxtmp);

  for (flux=0; flux<nf; flux++)
    Psiminusp0[flux]=max(
                      -fabs(gl->cycle.resconv.xi*Gminusp1h[flux]/notzero(Gminusp1h[flux]-Gminusm1h[flux],1.0e-30)),
                       min(fluxtmp[flux]/2.0/notzero(Gminusp1h[flux]-Gminusm1h[flux],1.0e-30),
                          fabs(gl->cycle.resconv.xi*Gminusp1h[flux]/notzero(Gminusp1h[flux]-Gminusm1h[flux],1.0e-30)))  
                        );


  /*here find Psiplusp0 */
  for (flux=0; flux<nf; flux++)
    fluxtmp[flux]=(lambdap0[flux][flux]+lambdapp0[flux][flux])*minmod3(Mm1[flux],Mp0[flux],Mp1[flux]);
  multiply_matrix_and_vector(Linvp0,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(Lm1h,fluxtmp2,fluxtmp);

  for (flux=0; flux<nf; flux++)
    Psiplusp0[flux]=max(
                      -fabs(gl->cycle.resconv.xi*Gplusm1h[flux]/notzero(Gplusp1h[flux]-Gplusm1h[flux],1.0e-30)),
                       min(fluxtmp[flux]/2.0/notzero(Gplusp1h[flux]-Gplusm1h[flux],1.0e-30),
                          fabs(gl->cycle.resconv.xi*Gplusm1h[flux]/notzero(Gplusp1h[flux]-Gplusm1h[flux],1.0e-30)))  
                        );

  multiply_matrix_and_vector(Linvm1h,Gplusm1h,fluxtmp1);
  multiply_matrix_and_vector(Linvp1h,Gminusp1h,fluxtmp2);
  for (flux=0; flux<nf; flux++)
     fluxtmp[flux]=0.5*Psiplusp0[flux]*(Gplusp1h[flux]-Gplusm1h[flux]);
  multiply_matrix_and_vector(Linvm1h,fluxtmp,fluxtmp3);
  for (flux=0; flux<nf; flux++)
     fluxtmp[flux]=-0.5*Psiminusp0[flux]*(Gminusp1h[flux]-Gminusm1h[flux]);
  multiply_matrix_and_vector(Linvp1h,fluxtmp,fluxtmp4);


  for (flux=0; flux<nf; flux++) Fint[flux]=fluxtmp1[flux]+fluxtmp2[flux]+fluxtmp3[flux]+fluxtmp4[flux];

}




void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  jacvars_t *jacvarstmp,*jacvarsm3h,*jacvarsm1h,*jacvarsp1h,*jacvarsp3h;
  flux_t Fm1h,Fp1h;
  long flux,l;
  sqmat_t Lambdaminus_p0,Lambdaplus_p0,Lambdaminus_p1;

  jacvarsm3h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsm1h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsp1h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsp3h= (jacvars_t *) malloc(sizeof(jacvars_t));


  l=_l_minus_one(ls,gl,theta);
  if (is_node_bdry(np[l],TYPELEVEL_FLUID_WORK))
    find_jacvars_at_interface(np,gl,l,l,theta,gl->cycle.resconv.AVERAGING,jacvarsm1h);
    else find_jacvars_at_interface(np,gl,_al(gl,l,theta,-1),_al(gl,l,theta,+0),theta,gl->cycle.resconv.AVERAGING,jacvarsm1h);
  find_jacvars_at_interface(np,gl,_al(gl,l,theta,+0),_al(gl,l,theta,+1),theta,gl->cycle.resconv.AVERAGING,jacvarsp1h);
  find_jacvars_at_interface(np,gl,_al(gl,l,theta,+1),_al(gl,l,theta,+2),theta,gl->cycle.resconv.AVERAGING,jacvarsp3h);


    if (!is_node_bdry(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK)  && gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER) {
      find_Fstar_interface_2o(np, gl,_al(gl,l,theta,-1),_al(gl,l,theta,+0),_al(gl,l,theta,+1),_al(gl,l,theta,+2),
               theta, *jacvarsm1h,*jacvarsp1h,*jacvarsp3h,Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    } else {
      find_Fstar_interface_1o(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, *jacvarsp1h, Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    }




  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    jacvarstmp=jacvarsm3h;
    jacvarsm3h=jacvarsm1h;
    jacvarsm1h=jacvarsp1h;
    jacvarsp1h=jacvarsp3h;
    jacvarsp3h=jacvarstmp;

    for (flux=0; flux<nf; flux++) Fm1h[flux]=Fp1h[flux];
    for (flux=0; flux<nf; flux++) Lambdaminus_p0[flux][flux]=Lambdaminus_p1[flux][flux];

    if (!is_node_bdry(np[_al(gl,l,theta,+1)],TYPELEVEL_FLUID_WORK))
    find_jacvars_at_interface(np,gl,_al(gl,l,theta,+1),_al(gl,l,theta,+2),theta,gl->cycle.resconv.AVERAGING,jacvarsp3h);
    else find_jacvars_at_interface(np,gl,_al(gl,l,theta,+1),_al(gl,l,theta,+1),theta,gl->cycle.resconv.AVERAGING,jacvarsp3h);

    if (!is_node_bdry(np[_al(gl,l,theta,+1)],TYPELEVEL_FLUID_WORK) && gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER) {
      find_Fstar_interface_2o(np, gl, _al(gl,l,theta,-1),_al(gl,l,theta,+0),_al(gl,l,theta,+1),_al(gl,l,theta,+2),
               theta, *jacvarsm1h,*jacvarsp1h,*jacvarsp3h,Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    } else {
      find_Fstar_interface_1o(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, *jacvarsp1h, Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    }
    for (flux=0; flux<nf; flux++){
      np[l].bs->Delta_Lambda[theta][flux]=-Lambdaminus_p0[flux][flux]    +Lambdaplus_p0[flux][flux];
    }



#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
    for (flux=0; flux<nf; flux++) {
      np[l].bs->trapezoidalm1_next[flux]+=gl->cycle.restime.weightm1_trapezoidal_convection*(Fp1h[flux]-Fm1h[flux]);
      np[l].wk->Res[flux]+=(1.0-gl->cycle.restime.weightm1_trapezoidal_convection)*(Fp1h[flux]-Fm1h[flux]);
    }
#else
    for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=(Fp1h[flux]-Fm1h[flux]);
#endif



  }

  free(jacvarsm3h);
  free(jacvarsm1h);
  free(jacvarsp1h);
  free(jacvarsp3h);
}




void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  long flux;
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=np[l].bs->Delta_Lambda[theta][flux];
}





/* same as for the Yee-Roe non-positivity-preserving scheme */
void find_Delta_Lambda_for_dtau_non_positive(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  sqmat_t Lambda_minus,Lambda_plus;
  long flux;
  jacvars_t jacvarsm1h,jacvarsp1h;
  metrics_t metricsm1h, metricsp1h;


  find_metrics_at_interface(np,gl,l,_al(gl,l,theta,+1),theta,&metricsp1h);
  find_metrics_at_interface(np,gl,_al(gl,l,theta,-1),l,theta,&metricsm1h);


  find_jacvars_at_interface(np, gl, l, _al(gl,l,theta,+1), theta, gl->cycle.resconv.AVERAGING,&jacvarsp1h);
  find_jacvars_at_interface(np, gl, _al(gl,l,theta,-1), l, theta, gl->cycle.resconv.AVERAGING,&jacvarsm1h);

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



