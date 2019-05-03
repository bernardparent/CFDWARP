// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2011-2012 Bernard Parent

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


#include <cycle/resconv/_resconv.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <cycle/share/res_share.h>
#include <src/control.h>

#define ACCURACY_FIRSTORDER 1
#define ACCURACY_SECONDORDER 2

#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif



void write_disc_resconv_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    LIMITER=LIMITER_VANLEER;\n"
    "    ACCURACY=ACCURACY_SECONDORDER;\n"
    "    xi=1.99;\n"
    "    EIGENVALCOND=EIGENVALCOND_GNOFFO;\n"
    "  );\n"
  ,_RESCONV_ACTIONNAME);

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
    SOAP_add_int_to_vars(codex,"LIMITER_MINMOD",LIMITER_MINMOD);
    SOAP_add_int_to_vars(codex,"LIMITER_VANLEER",LIMITER_VANLEER);
    SOAP_add_int_to_vars(codex,"LIMITER_SUPERBEE",LIMITER_SUPERBEE);
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
    find_int_var_from_codex(codex,"LIMITER",&gl->cycle.resconv.LIMITER);
    if (gl->cycle.resconv.LIMITER!=LIMITER_MINMOD && gl->cycle.resconv.LIMITER!=LIMITER_VANLEER && gl->cycle.resconv.LIMITER!=LIMITER_SUPERBEE)
      SOAP_fatal_error(codex,"LIMITER must be set to either LIMITER_MINMOD, LIMITER_VANLEER, or LIMITER_SUPERBEE.");

    find_int_var_from_codex(codex,"ACCURACY",&gl->cycle.resconv.ACCURACY);
    if (gl->cycle.resconv.ACCURACY!=ACCURACY_FIRSTORDER && gl->cycle.resconv.ACCURACY!=ACCURACY_SECONDORDER)
      SOAP_fatal_error(codex,"ACCURACY must be set to either ACCURACY_FIRSTORDER or ACCURACY_SECONDORDER.");

    find_double_var_from_codex(codex,"xi",&gl->cycle.resconv.xi);
    if (gl->cycle.resconv.xi<0.0 || gl->cycle.resconv.xi>2.0)
      SOAP_fatal_error(codex,"xi must be set to a value between 0 and 2.");

    gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FVSPLUS;
    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}




static void find_Fplus_Fminus(np_t *np, gl_t *gl, long l, long theta, metrics_t metrics, flux_t Deltaxt_p1h, flux_t Fplus, flux_t Fminus, flux_t Gplus, flux_t Gminus, sqmat_t Linv, sqmat_t L, sqmat_t Lambdaminus, sqmat_t Lambdaplus ){
  long row,col;
  sqmat_t lambdap,lambda;
  flux_t fluxtmp;
  jacvars_t jacvars;
#ifdef _RESTIME_CDF
  sqmat_t lambdaxt;
  long flux;
#endif


  find_jacvars(np[l],gl,metrics,theta,&jacvars);
  find_conditioned_Lambda_absolute_from_jacvars(jacvars, metrics, gl->cycle.resconv.EIGENVALCOND, lambdap);
  find_Lambda_from_jacvars(jacvars, metrics, lambda);
  //set_matrix_to_zero(lambdap);
  //for (flux=0; flux<nf; flux++) lambdap[flux][flux]=fabs(lambda[flux][flux]);
#ifdef _RESTIME_CDF
  find_Lambdaxt(np, gl, l, jacvars, metrics, lambdaxt);
  for (flux=0; flux<nf; flux++) {
    lambda[flux][flux]+=-lambdaxt[flux][flux]*Deltaxt_p1h[flux];
    lambdap[flux][flux]+=-fabs(lambdaxt[flux][flux])*Deltaxt_p1h[flux];
  }
#endif


  for (row=0; row<nf; row++){  
    for (col=0; col<nf; col++){
      Lambdaplus[row][col]=0.5*max(0.0,lambda[row][col]+lambdap[row][col])+0.5*max(0.0,lambda[row][col]-lambdap[row][col]);
      Lambdaminus[row][col]=0.5*min(0.0,lambda[row][col]+lambdap[row][col])+0.5*min(0.0,lambda[row][col]-lambdap[row][col]);
      
    }
  }

//  condition_Lambda_plus_minus(jacvars, jacvars, metrics, Lambdaplus, Lambdaminus);
  find_Linv_from_jacvars(jacvars, metrics, Linv);
  find_L_from_jacvars(jacvars, metrics, L); 
  find_LUstar_from_jacvars(jacvars, metrics, fluxtmp);

  multiply_diagonal_matrix_and_vector(Lambdaplus,fluxtmp,Gplus);

  multiply_matrix_and_vector(Linv,Gplus,Fplus);
  
  multiply_diagonal_matrix_and_vector(Lambdaminus,fluxtmp,Gminus);

  multiply_matrix_and_vector(Linv,Gminus,Fminus);

}



static void find_Fstar_interface_2o(np_t *np, gl_t *gl, long lm3h, long lm1h, long lp1h, long lp3h,
                     long theta, flux_t Fint, sqmat_t Lambdaminus_p1h, sqmat_t Lambdaplus_m1h){
  flux_t Fplusm1h,Fplusm3h,Fplusp1h,Fplusp3h;
  flux_t Fminusm1h,Fminusm3h,Fminusp1h,Fminusp3h;
  flux_t Gplusm1h,Gplusm3h,Gplusp1h,Gplusp3h;
  flux_t Gminusm1h,Gminusm3h,Gminusp1h,Gminusp3h;
  sqmat_t Linvm1h, Linvm3h, Linvp1h, Linvp3h;
  sqmat_t Lm1h, Lm3h, Lp1h, Lp3h, mattmp1, mattmp2;
  flux_t psiplus,psiminus,fluxtmp,fluxtmpplus,fluxtmpminus,phiplus,phiminus;
  metrics_t metrics;
  flux_t Deltaxt_p0;
  long flux;
#ifdef _RESTIME_CDF
  jacvars_t jacvarsm1h,jacvarsp1h;
#endif

  find_metrics_at_interface(np, gl, lm1h, lp1h, theta, &metrics);

#ifdef _RESTIME_CDF
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_Deltaxt_interface(np, gl, lm1h, lp1h, theta, jacvarsm1h, jacvarsp1h, metrics, Deltaxt_p0);
#else
  for (flux=0; flux<nf; flux++) Deltaxt_p0[flux]=0.0;
#endif


  find_Fplus_Fminus(np, gl, lm1h, theta, metrics, Deltaxt_p0, Fplusm1h, Fminusm1h, Gplusm1h, Gminusm1h, Linvm1h, Lm1h, mattmp1, Lambdaplus_m1h);
  find_Fplus_Fminus(np, gl, lm3h, theta, metrics, Deltaxt_p0, Fplusm3h, Fminusm3h, Gplusm3h, Gminusm3h, Linvm3h, Lm3h, mattmp1, mattmp2);
  find_Fplus_Fminus(np, gl, lp1h, theta, metrics, Deltaxt_p0, Fplusp1h, Fminusp1h, Gplusp1h, Gminusp1h, Linvp1h, Lp1h, Lambdaminus_p1h, mattmp2);
  find_Fplus_Fminus(np, gl, lp3h, theta, metrics, Deltaxt_p0, Fplusp3h, Fminusp3h, Gplusp3h, Gminusp3h, Linvp3h, Lp3h, mattmp1, mattmp2);


  for (flux=0; flux<nf; flux++){
    phiplus[flux]=_limiter_TVD((Fplusp1h[flux]-Fplusm1h[flux])/notzero(Fplusm1h[flux]-Fplusm3h[flux],1e-99),gl->cycle.resconv.LIMITER);
    phiminus[flux]=_limiter_TVD((Fminusm1h[flux]-Fminusp1h[flux])/notzero(Fminusp1h[flux]-Fminusp3h[flux],1e-99),gl->cycle.resconv.LIMITER);
  }
  
  for (flux=0; flux<nf; flux++){
    fluxtmp[flux]=phiminus[flux]*(Fminusp1h[flux]-Fminusp3h[flux]);
  }
  multiply_matrix_and_vector(Lp1h,fluxtmp,fluxtmpminus);
  
  for (flux=0; flux<nf; flux++){
    fluxtmp[flux]=phiplus[flux]*(Fplusm1h[flux]-Fplusm3h[flux]);
  }
  multiply_matrix_and_vector(Lm1h,fluxtmp,fluxtmpplus);

  for (flux=0; flux<nf; flux++) {
    if (TRUE){
      /* includes positivity-preserving limiting process */
      psiplus[flux]=max(-fabs(gl->cycle.resconv.xi*Gplusm1h[flux]/notzero(Gplusm1h[flux]-Gplusm3h[flux],1.0e-40)),
                    min(fluxtmpplus[flux]/notzero(Gplusm1h[flux]-Gplusm3h[flux],1.0e-40),
                        fabs(gl->cycle.resconv.xi*Gplusm1h[flux]/notzero(Gplusm1h[flux]-Gplusm3h[flux],1.0e-40))
                   ) );
      psiminus[flux]=max(-fabs(gl->cycle.resconv.xi*Gminusp1h[flux]/notzero(Gminusp1h[flux]-Gminusp3h[flux],1.0e-40)), 
                    min(fluxtmpminus[flux]/notzero(Gminusp1h[flux]-Gminusp3h[flux],1.0e-40),
                        fabs(gl->cycle.resconv.xi*Gminusp1h[flux]/notzero(Gminusp1h[flux]-Gminusp3h[flux],1.0e-40)) 
                   ) );
    } else {
      /* excludes positivity-preserving limiting process */
      psiplus[flux]=fluxtmpplus[flux]/notzero(Gplusm1h[flux]-Gplusm3h[flux],1.0e-40);
      psiminus[flux]=fluxtmpminus[flux]/notzero(Gminusp1h[flux]-Gminusp3h[flux],1.0e-40); 
    }              

  }


  for (flux=0; flux<nf; flux++){
    fluxtmp[flux]=psiplus[flux]*notzero(Gplusm1h[flux]-Gplusm3h[flux],1.0e-40);
  }
  multiply_matrix_and_vector(Linvm1h,fluxtmp,fluxtmpplus);

  for (flux=0; flux<nf; flux++){
    fluxtmp[flux]=psiminus[flux]*notzero(Gminusp1h[flux]-Gminusp3h[flux],1.0e-40);
  }
  multiply_matrix_and_vector(Linvp1h,fluxtmp,fluxtmpminus);

  for (flux=0; flux<nf; flux++) {
     Fint[flux]=Fplusm1h[flux]+0.5*fluxtmpplus[flux]
               +Fminusp1h[flux]+0.5*fluxtmpminus[flux];
  }

}



static void find_Fstar_interface_1o(np_t *np, gl_t *gl, long lm1h, long lp1h, 
                     long theta, flux_t Fint, sqmat_t Lambdaminus_p1h, sqmat_t Lambdaplus_m1h){
  flux_t Fplusm1h,Fplusp1h,Gplusm1h,Gplusp1h;
  flux_t Fminusm1h,Fminusp1h,Gminusm1h,Gminusp1h;
  sqmat_t Linvm1h,Linvp1h,Lm1h,Lp1h,mattmp1, mattmp2;
  metrics_t metrics;
  long flux;
  flux_t Deltaxt_p0;
#ifdef _RESTIME_CDF
  jacvars_t jacvarsm1h,jacvarsp1h;
#endif

  find_metrics_at_interface(np, gl, lm1h, lp1h, theta, &metrics);

#ifdef _RESTIME_CDF
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_Deltaxt_interface(np, gl, lm1h, lp1h, theta, jacvarsm1h, jacvarsp1h, metrics, Deltaxt_p0);
#else
  for (flux=0; flux<nf; flux++) Deltaxt_p0[flux]=0.0;
#endif

  find_Fplus_Fminus(np, gl, lm1h, theta, metrics, Deltaxt_p0, Fplusm1h, Fminusm1h, Gplusm1h, Gminusm1h, Linvm1h, Lm1h, mattmp1, Lambdaplus_m1h);
  find_Fplus_Fminus(np, gl, lp1h, theta, metrics, Deltaxt_p0, Fplusp1h, Fminusp1h, Gplusp1h, Gminusp1h, Linvp1h, Lp1h, Lambdaminus_p1h, mattmp2);

  for (flux=0; flux<nf; flux++) {
     Fint[flux]=Fplusm1h[flux]+Fminusp1h[flux];
  }

}




void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  flux_t Fm1h,Fp1h;
  long flux,l;
  sqmat_t Lambdaminus_p0,Lambdaplus_p0,Lambdaminus_p1;

  l=_l_minus_one(ls,gl,theta);

    if (!is_node_bdry(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK) && gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER) {
      find_Fstar_interface_2o(np, gl,_al(gl,l,theta,-1),_al(gl,l,theta,+0),_al(gl,l,theta,+1),_al(gl,l,theta,+2),
               theta, Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    } else {
      find_Fstar_interface_1o(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta,  Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    }


  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){

    for (flux=0; flux<nf; flux++) Fm1h[flux]=Fp1h[flux];
    for (flux=0; flux<nf; flux++) Lambdaminus_p0[flux][flux]=Lambdaminus_p1[flux][flux];


    if (!is_node_bdry(np[_al(gl,l,theta,+1)],TYPELEVEL_FLUID_WORK) && gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER) {
      find_Fstar_interface_2o(np, gl, _al(gl,l,theta,-1),_al(gl,l,theta,+0),_al(gl,l,theta,+1),_al(gl,l,theta,+2),
               theta, Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    } else {
      find_Fstar_interface_1o(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta,  Fp1h,Lambdaminus_p1,Lambdaplus_p0);
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

}



void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  long flux;
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=np[l].bs->Delta_Lambda[theta][flux];
}





