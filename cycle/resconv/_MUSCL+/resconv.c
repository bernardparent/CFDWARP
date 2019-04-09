#include <cycle/resconv/_resconv.h>
#include <cycle/share/res_share.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>

#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif


#define FLUX_FDSplus 1
#define FLUX_FVSplus 2
#define FLUX_FDSplusFilter 4
#define FLUX_FVSplusFilter 5
#define INTERPOL_AOWENO5 1
#define INTERPOL_WENO5 2
#define INTERPOL_WENO3 3
#define INTERPOL_TVD2_MINMOD 4
#define INTERPOL_TVD2_VANLEER 5
#define INTERPOL_TVD2_SUPERBEE 6
#define INTERPOL_TVD5 7
#define INTERPOL_CWENO3 8
#define INTERPOL_AOWENO7 9


void write_disc_resconv_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    FLUX=FLUX_FDSplus;\n"
    "    numiter=%ld;\n"
    "    AVERAGING=AVERAGING_ARITH;\n"
    "    AOWENO_TYPE=AOWENO_TYPE_COMPRESSIVE;\n"
    "    AOWENO_gammalo=0.95;\n"
    "    AOWENO_gammahi=0.9999;\n"
    "    INTERPOL=INTERPOL_AOWENO7;\n"
    "    EIGENVALCOND=EIGENVALCOND_PASCAL;\n"
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
    SOAP_add_int_to_vars(codex,"EIGENVALCOND_PASCAL",EIGENVALCOND_PASCAL);
    SOAP_add_int_to_vars(codex,"EIGENVALCOND_GNOFFO",EIGENVALCOND_GNOFFO);
    SOAP_add_int_to_vars(codex,"AVERAGING_ROE",AVERAGING_ROE);
    SOAP_add_int_to_vars(codex,"AVERAGING_ARITH",AVERAGING_ARITH);
    SOAP_add_int_to_vars(codex,"FLUX_FDSplus",FLUX_FDSplus); 
    SOAP_add_int_to_vars(codex,"FLUX_FDSplusFilter",FLUX_FDSplusFilter); 
    SOAP_add_int_to_vars(codex,"FLUX_FVSplusFilter",FLUX_FVSplusFilter); 
    SOAP_add_int_to_vars(codex,"FLUX_FVSplus",FLUX_FVSplus); 
    SOAP_add_int_to_vars(codex,"INTERPOL_AOWENO7",INTERPOL_AOWENO7); 
    SOAP_add_int_to_vars(codex,"INTERPOL_AOWENO5",INTERPOL_AOWENO5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO5",INTERPOL_WENO5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO3",INTERPOL_WENO3); 
    SOAP_add_int_to_vars(codex,"INTERPOL_CWENO3",INTERPOL_CWENO3); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_MINMOD",INTERPOL_TVD2_MINMOD); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_VANLEER",INTERPOL_TVD2_VANLEER); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_SUPERBEE",INTERPOL_TVD2_SUPERBEE); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD5",INTERPOL_TVD5); 
    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_COMPRESSIVE",AOWENO_TYPE_COMPRESSIVE); 
    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_DIFFUSIVE",AOWENO_TYPE_DIFFUSIVE); 

    gl->DISC_RESCONV_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_resconv_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"EIGENVALCOND",&gl->cycle.resconv.EIGENVALCOND);
    if (gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PECLET && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PASCAL && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_HARTEN && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_GNOFFO)
      SOAP_fatal_error(codex,"EIGENVALCOND must be set to either EIGENVALCOND_PECLET, EIGENVALCOND_PASCAL, EIGENVALCOND_HARTEN, EIGENVALCOND_GNOFFO.");
    find_int_var_from_codex(codex,"FLUX",&gl->cycle.resconv.FLUX);
    if (gl->cycle.resconv.FLUX!=FLUX_FDSplus && gl->cycle.resconv.FLUX!=FLUX_FVSplus && gl->cycle.resconv.FLUX!=FLUX_FDSplusFilter && gl->cycle.resconv.FLUX!=FLUX_FVSplusFilter)
      SOAP_fatal_error(codex,"FLUX must be set to either FLUX_FDSplus or FLUX_FVSplus or FLUX_FDSplusFilter or FLUX_FVSplusFilter.");
    find_int_var_from_codex(codex,"numiter",&gl->cycle.resconv.numiter);

    find_int_var_from_codex(codex,"AVERAGING",&gl->cycle.resconv.AVERAGING);
    if (gl->cycle.resconv.AVERAGING!=AVERAGING_ROE  && gl->cycle.resconv.AVERAGING!=AVERAGING_ARITH )
      SOAP_fatal_error(codex,"AVERAGING must be set to either AVERAGING_ROE or AVERAGING_ARITH.");

    find_int_var_from_codex(codex,"INTERPOL",&gl->cycle.resconv.INTERPOL);
    if (gl->cycle.resconv.INTERPOL!=INTERPOL_AOWENO5 && gl->cycle.resconv.INTERPOL!=INTERPOL_AOWENO7 && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO3
     && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO5 && gl->cycle.resconv.INTERPOL!=INTERPOL_CWENO3
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_MINMOD && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_VANLEER 
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_SUPERBEE && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD5)
      SOAP_fatal_error(codex,"INTERPOL must be set to either INTERPOL_AOWENO7, INTERPOL_AOWENO5, INTERPOL_WENO5, INTERPOL_WENO3, INTERPOL_CWENO3, INTERPOL_TVD2_MINMOD, INTERPOL_TVD2_VANLEER, INTERPOL_TVD2_SUPERBEE, or INTERPOL_TVD5.");
    find_int_var_from_codex(codex,"AOWENO_TYPE",&gl->cycle.resconv.AOWENO_TYPE);
    if (gl->cycle.resconv.AOWENO_TYPE!=AOWENO_TYPE_COMPRESSIVE && gl->cycle.resconv.AOWENO_TYPE!=AOWENO_TYPE_DIFFUSIVE)
      SOAP_fatal_error(codex,"AOWENO_TYPE must be set to either AOWENO_TYPE_COMPRESSIVE or AOWENO_TYPE_DIFFUSIVE.");

    find_double_var_from_codex(codex,"AOWENO_gammalo",&gl->cycle.resconv.AOWENO_gammalo);
    find_double_var_from_codex(codex,"AOWENO_gammahi",&gl->cycle.resconv.AOWENO_gammahi);

    gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FDSPLUS;
    if (gl->cycle.resconv.FLUX==FLUX_FVSplus || gl->cycle.resconv.FLUX==FLUX_FVSplusFilter) gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FVSPLUS;


    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}


static void find_Fstar_interface(np_t *np, gl_t *gl, long l, long theta, musclvarscycle_t musclvarsm7h, musclvarscycle_t musclvarsm5h, musclvarscycle_t musclvarsm3h, musclvarscycle_t musclvarsm1h, musclvarscycle_t musclvarsp1h, musclvarscycle_t musclvarsp3h, musclvarscycle_t musclvarsp5h, musclvarscycle_t musclvarsp7h, flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h){
  metrics_t metrics;
  long flux,nodes_from_bdry;
  flux_t musclvarsL,musclvarsR,Finttmp;
  int AOWENO_TYPE;
  double gammalo,gammahi;

  find_metrics_at_interface(np, gl, l, _al(gl,l,theta,+1), theta, &metrics);
  nodes_from_bdry=_nodes_from_bdry_through_links_limited(np,gl,l,theta,TYPELEVEL_FLUID_WORK,hbw_resconv_fluid);
  for (flux=0; flux<nf; flux++) {
    AOWENO_TYPE=gl->cycle.resconv.AOWENO_TYPE;
    gammalo=gl->cycle.resconv.AOWENO_gammalo;
    gammahi=gl->cycle.resconv.AOWENO_gammahi;
    switch (gl->cycle.resconv.INTERPOL){
      case INTERPOL_AOWENO7:
        if (nodes_from_bdry>=4){
          musclvarsL[flux]=_f_AOWENO7(musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],0.5,gammalo,gammahi,AOWENO_TYPE);
          musclvarsR[flux]=_f_AOWENO7(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux],-0.5,gammalo,gammahi,AOWENO_TYPE); 
        } else {
          if (nodes_from_bdry>=3){
            musclvarsL[flux]=_f_AOWENO5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],0.5,0.95,0.95,AOWENO_TYPE);
            musclvarsR[flux]=_f_AOWENO5(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],-0.5,0.95,0.95,AOWENO_TYPE); 
          } else {
//            musclvarsL[flux]=_f_TVD5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
//            musclvarsR[flux]=_f_TVD5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
            musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_VANLEER);
            musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_VANLEER);
          } 
        }
      break;
      case INTERPOL_AOWENO5:
        if (nodes_from_bdry>=3){
          musclvarsL[flux]=_f_AOWENO5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],0.5,gammalo,gammahi,AOWENO_TYPE);
          musclvarsR[flux]=_f_AOWENO5(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],-0.5,gammalo,gammahi,AOWENO_TYPE); 
        } else {
          musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_VANLEER);
          musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_VANLEER);
//          musclvarsL[flux]=_f_TVD5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
//          musclvarsR[flux]=_f_TVD5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
        }
      break;
      case INTERPOL_TVD5:
        musclvarsL[flux]=_f_TVD5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
        musclvarsR[flux]=_f_TVD5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
      break;
      case INTERPOL_TVD2_MINMOD:
        musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_MINMOD);
        musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_MINMOD);
      break;
      case INTERPOL_TVD2_VANLEER:
        musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_VANLEER);
        musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_VANLEER);
      break;
      case INTERPOL_TVD2_SUPERBEE:
        musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_SUPERBEE);
        musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_SUPERBEE);
      break;
      case INTERPOL_WENO5:
        musclvarsL[flux]=_f_WENO5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
        musclvarsR[flux]=_f_WENO5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
      break;
      case INTERPOL_WENO3:
        musclvarsL[flux]=_f_WENO3(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux]);
        musclvarsR[flux]=_f_WENO3(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux]);
      break;
      case INTERPOL_CWENO3:
        if (nodes_from_bdry>=3 ){
          musclvarsL[flux]=_f_CWENO3(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],0.5);
          musclvarsR[flux]=_f_CWENO3(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],-0.5);
        } else {
          musclvarsL[flux]=_f_TVD5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
          musclvarsR[flux]=_f_TVD5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
        } 
      break;
      default:
        fatal_error("gl->cycle.resconv.INTERPOL invalid in find_Fstar_interface().");
    } 
  }

  switch (gl->cycle.resconv.FLUX){
    case FLUX_FDSplus:
      find_Fstar_interface_FDSplus_muscl(np, gl,  l, _al(gl,l,theta,+1), theta, musclvarsL, musclvarsR, metrics,   gl->cycle.resconv.numiter, gl->cycle.resconv.EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint, lambdaminusp1h,  lambdaplusm1h);
    break;
    case FLUX_FVSplus:
      find_Fstar_interface_FVSplus_muscl(np, gl,  l, _al(gl,l,theta,+1), theta, musclvarsL, musclvarsR, metrics,   gl->cycle.resconv.numiter, gl->cycle.resconv.EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint, lambdaminusp1h,  lambdaplusm1h);
    break;
    case FLUX_FDSplusFilter:
      find_Fstar_interface_FDS_muscl(np, gl, l, _al(gl,l,theta,+1), theta, musclvarsL, musclvarsR, metrics, EIGENVALCOND_NONE, gl->cycle.resconv.AVERAGING, Finttmp);
      filter_Fstar_interface_positivity_preserving(np, gl, l, theta, metrics, gl->cycle.resconv.numiter, gl->cycle.resconv.EIGENVALCOND, Finttmp, Fint, lambdaminusp1h,  lambdaplusm1h);
    break;
    case FLUX_FVSplusFilter:
      find_Fstar_interface_FVS_muscl(np, gl, l, _al(gl,l,theta,+1), theta, musclvarsL, musclvarsR, metrics, EIGENVALCOND_NONE, gl->cycle.resconv.AVERAGING, Finttmp);
      filter_Fstar_interface_positivity_preserving(np, gl, l, theta, metrics, gl->cycle.resconv.numiter, gl->cycle.resconv.EIGENVALCOND, Finttmp, Fint, lambdaminusp1h,  lambdaplusm1h);
    break;
    default:
      fatal_error("gl->cycle.resconv.FLUX must be set to either FLUX_FDSplus, FLUX_FVSplus, or FLUX_FDSplusFilter, or FLUX_FVSplusFilter in find_Fstar_interface().");
  }


}


void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  flux_t Fm1h,Fp1h;
  long flux,l;
  sqmat_t Lambdaminus_p0,Lambdaplus_p0,Lambdaminus_p1;
  musclvarscycle_t musclvarsm4,musclvarsm3,musclvarsm2,musclvarsm1,musclvarsp0,musclvarsp1,
         musclvarsp2,musclvarsp3,musclvarsp4;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){

    if(l==ls){
      if (hbw_resconv_fluid>=4) find_musclvarscycle_offset(np,gl,l,theta,-4,musclvarsm4);
      find_musclvarscycle_offset(np,gl,l,theta,-3,musclvarsm3);
      find_musclvarscycle_offset(np,gl,l,theta,-2,musclvarsm2);
      find_musclvarscycle_offset(np,gl,l,theta,-1,musclvarsm1);
      find_musclvarscycle_offset(np,gl,l,theta,+0,musclvarsp0);
      find_musclvarscycle_offset(np,gl,l,theta,+1,musclvarsp1);
      find_musclvarscycle_offset(np,gl,l,theta,+2,musclvarsp2);
      find_musclvarscycle_offset(np,gl,l,theta,+3,musclvarsp3);
      if (hbw_resconv_fluid>=4) find_musclvarscycle_offset(np,gl,l,theta,+4,musclvarsp4);


      find_Fstar_interface(np, gl, _al(gl,l,theta,-1), theta, musclvarsm4, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, Fm1h,Lambdaminus_p1,Lambdaplus_p0);
    } else {
      for (flux=0; flux<nf; flux++) {
        Fm1h[flux]=Fp1h[flux];
        if (hbw_resconv_fluid>=4) musclvarsm3[flux]=musclvarsm2[flux];
        musclvarsm2[flux]=musclvarsm1[flux];
        musclvarsm1[flux]=musclvarsp0[flux];
        musclvarsp0[flux]=musclvarsp1[flux];
        musclvarsp1[flux]=musclvarsp2[flux];
        musclvarsp2[flux]=musclvarsp3[flux];
        if (hbw_resconv_fluid>=4) musclvarsp3[flux]=musclvarsp4[flux];
      }
    }


    for (flux=0; flux<nf; flux++) {
      Lambdaminus_p0[flux][flux]=Lambdaminus_p1[flux][flux];
    }

    find_musclvarscycle_offset(np,gl,l,theta,+3,musclvarsp3);


    find_Fstar_interface(np, gl, l, theta, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, musclvarsp4, Fp1h,Lambdaminus_p1,Lambdaplus_p0);
    for (flux=0; flux<nf; flux++){
      np[l].wk->Res[flux]+=fact*(Fp1h[flux]-Fm1h[flux]);
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
      np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*(Fp1h[flux]-Fm1h[flux]);
#endif
      np[l].bs->Delta_Lambda[theta][flux]=-Lambdaminus_p0[flux][flux]    +Lambdaplus_p0[flux][flux];
    }
  }

}


void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  long flux;
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=np[l].bs->Delta_Lambda[theta][flux];
}


