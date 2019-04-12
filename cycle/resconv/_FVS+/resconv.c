#include <cycle/resconv/_resconv.h>
#include <cycle/share/res_share.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>

#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif


#define INTERPOL_AOWENO5 1
#define INTERPOL_AOWENO7 2
#define INTERPOL_WENO3 4
#define INTERPOL_WENO5 5
#define INTERPOL_WENO7 6
#define INTERPOL_TVD2_MINMOD 8
#define INTERPOL_TVD2_VANLEER 9
#define INTERPOL_TVD2_SUPERBEE 10
#define INTERPOL_TVD5 11
#define INTERPOL_CWENO3 12
#define INTERPOL_FIRSTORDER 14
#define INTERPOL_BDF2 15



void write_disc_resconv_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    numiter=4;\n"
    "    AOWENO_TYPE=AOWENO_TYPE_DIFFUSIVE;\n"
    "    AOWENO_gammalo=0.85;\n"
    "    AOWENO_gammahi=0.85;\n"
    "    INTERPOL=INTERPOL_AOWENO5;\n"
    "    EIGENVALCOND=EIGENVALCOND_PASCAL;\n"
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
    SOAP_add_int_to_vars(codex,"EIGENVALCOND_PASCAL",EIGENVALCOND_PASCAL);
    SOAP_add_int_to_vars(codex,"INTERPOL_FIRSTORDER",INTERPOL_FIRSTORDER);
    SOAP_add_int_to_vars(codex,"INTERPOL_BDF2",INTERPOL_BDF2);
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_MINMOD",INTERPOL_TVD2_MINMOD); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_VANLEER",INTERPOL_TVD2_VANLEER); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_SUPERBEE",INTERPOL_TVD2_SUPERBEE); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD5",INTERPOL_TVD5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_AOWENO5",INTERPOL_AOWENO5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_AOWENO7",INTERPOL_AOWENO7); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO3",INTERPOL_WENO3); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO5",INTERPOL_WENO5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO7",INTERPOL_WENO7); 
    SOAP_add_int_to_vars(codex,"INTERPOL_CWENO3",INTERPOL_CWENO3); 

    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_COMPRESSIVE",AOWENO_TYPE_COMPRESSIVE); 
    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_DIFFUSIVE",AOWENO_TYPE_DIFFUSIVE); 

    gl->DISC_RESCONV_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_resconv_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"EIGENVALCOND",&gl->cycle.resconv.EIGENVALCOND);
    if (gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PECLET && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_HARTEN && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_GNOFFO && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PASCAL)
      SOAP_fatal_error(codex,"EIGENVALCOND must be set to either EIGENVALCOND_PECLET, EIGENVALCOND_HARTEN, EIGENVALCOND_GNOFFO, EIGENVALCOND_PASCAL.");

    find_int_var_from_codex(codex,"numiter",&gl->cycle.resconv.numiter);
    find_int_var_from_codex(codex,"INTERPOL",&gl->cycle.resconv.INTERPOL);
    if (
        gl->cycle.resconv.INTERPOL!=INTERPOL_FIRSTORDER  
#if (hbw_resconv_fluid>=2)
     && gl->cycle.resconv.INTERPOL!=INTERPOL_BDF2  
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_MINMOD && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_VANLEER 
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_SUPERBEE && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO3
#endif
#if (hbw_resconv_fluid>=3)
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD5 && gl->cycle.resconv.INTERPOL!=INTERPOL_AOWENO5 
     && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO5 && gl->cycle.resconv.INTERPOL!=INTERPOL_CWENO3
#endif
#if (hbw_resconv_fluid>=4)
     && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO7 && gl->cycle.resconv.INTERPOL!=INTERPOL_AOWENO7
#endif
     )
      SOAP_fatal_error(codex,"INTERPOL must be set to either "
        "INTERPOL_FIRSTORDER, "
#if (hbw_resconv_fluid>=2)      
        "INTERPOL_BDF2, INTERPOL_TVD2_MINMOD, INTERPOL_TVD2_VANLEER, INTERPOL_TVD2_SUPERBEE, INTERPOL_WENO3, "
#endif
#if (hbw_resconv_fluid>=3)      
        "INTERPOL_AOWENO5, INTERPOL_WENO5, INTERPOL_CWENO3, INTERPOL_TVD5,"
#endif
#if (hbw_resconv_fluid>=4)      
        "INTERPOL_AOWENO7, INTERPOL_WENO7,"
#endif
      );
    find_int_var_from_codex(codex,"AOWENO_TYPE",&gl->cycle.resconv.AOWENO_TYPE);
    if (gl->cycle.resconv.AOWENO_TYPE!=AOWENO_TYPE_COMPRESSIVE && gl->cycle.resconv.AOWENO_TYPE!=AOWENO_TYPE_DIFFUSIVE)
      SOAP_fatal_error(codex,"AOWENO_TYPE must be set to either AOWENO_TYPE_COMPRESSIVE or AOWENO_TYPE_DIFFUSIVE.");

    find_double_var_from_codex(codex,"AOWENO_gammalo",&gl->cycle.resconv.AOWENO_gammalo);
    find_double_var_from_codex(codex,"AOWENO_gammahi",&gl->cycle.resconv.AOWENO_gammahi);

    gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FVSPLUS;

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}


static void find_Fplus_Fminus(np_t *np, gl_t *gl, long l, long theta, metrics_t metrics,  flux_t Deltaxt_p1h, flux_t Fplus, flux_t Fminus ){
  long row,col;
  sqmat_t lambdap,lambda;
  flux_t fluxtmp;
  jacvars_t jacvars;
  flux_t Gplus, Gminus;
  sqmat_t Linv, L;
   sqmat_t Lambdaminus,  Lambdaplus;
#ifdef _RESTIME_CDF
  sqmat_t lambdaxt;
  long flux;
#endif

  assert(is_node_resumed(np[l]));
  assert(is_node_in_domain_lim(l,gl));
  assert(is_node_valid(np[l],TYPELEVEL_FLUID));

  find_jacvars(np[l],gl,metrics,theta,&jacvars);
  find_Lambda_from_jacvars(jacvars, metrics, lambda);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdap[row][col]=fabs(lambda[row][col]);
    }
  }
#ifdef _RESTIME_CDF
  find_Lambdaxt(np, gl, l, jacvars, metrics, lambdaxt);
  for (flux=0; flux<nf; flux++) {
    lambda[flux][flux]+=-lambdaxt[flux][flux]*Deltaxt_p1h[flux];
    lambdap[flux][flux]+=-fabs(lambdaxt[flux][flux])*Deltaxt_p1h[flux];
  }
#endif


  for (row=0; row<nf; row++){  
    for (col=0; col<nf; col++){
      Lambdaplus[row][col]=0.5*(lambda[row][col]+lambdap[row][col]);
      Lambdaminus[row][col]=0.5*(lambda[row][col]-lambdap[row][col]);
      
    }
  }


  find_Linv_from_jacvars(jacvars, metrics, Linv);
  find_L_from_jacvars(jacvars, metrics, L); 
  find_LUstar_from_jacvars(jacvars, metrics, fluxtmp);

  multiply_diagonal_matrix_and_vector(Lambdaplus,fluxtmp,Gplus);

  multiply_matrix_and_vector(Linv,Gplus,Fplus);
  
  multiply_diagonal_matrix_and_vector(Lambdaminus,fluxtmp,Gminus);

  multiply_matrix_and_vector(Linv,Gminus,Fminus);

}


/* find Fplusminus at node _al(gl,l,theta,offset) if it's a valid node
   if not, then find Fplusminus at valid node in between l and _al(gl,l,theta,offset) that is closest
   to _al(gl,l,theta,offset) */
void find_Fplus_Fminus_offset(np_t *np, gl_t *gl, long l, long theta, long offset, metrics_t metrics,  flux_t Deltaxt_p1h, flux_t Fplus, flux_t Fminus){
  long cnt;

  if (!is_node_valid(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK)) fatal_error("must start with a valid node in find_Fplus_Fminus_offset");
  cnt=0;
  do {
    if (offset>0) cnt++; else cnt--;
    assert(is_node_in_domain_lim(l, gl));
  } while(is_node_valid(np[_al(gl,l,theta,cnt)],TYPELEVEL_FLUID_WORK) && abs(cnt)<=abs(offset));
  if (offset>0) cnt=min(offset,cnt-1); else cnt=max(cnt+1,offset);
  if (!is_node_valid(np[_al(gl,l,theta,cnt)],TYPELEVEL_FLUID_WORK)) fatal_error("problem finding Fplus/Fminus in find_Fplus_Fminus_offset");
  find_Fplus_Fminus(np, gl, _al(gl,l,theta,cnt), theta, metrics,  Deltaxt_p1h, Fplus, Fminus);
}



static void find_Fstar_interface(np_t *np, gl_t *gl, long l, long theta, flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h){
  metrics_t metrics;
  long flux,nodes_from_bdry;
  flux_t Fplusm1h,Fminusm1h,Fplusp1h,Fminusp1h;
  flux_t Fplusm3h,Fminusm3h,Fplusp3h,Fminusp3h;
  flux_t Fplusm5h,Fminusm5h,Fplusp5h,Fminusp5h;
  flux_t Fminusm7h,Fminusp7h,Fplusm7h,Fplusp7h;
  int AOWENO_TYPE;
  double gammalo,gammahi;
  flux_t Deltaxt,FplusL,FminusR,Finttmp;
  long hbw_interpol;
#ifdef _RESTIME_CDF  
  jacvars_t jacvarsm1h,jacvarsp1h;
#endif

  find_metrics_at_interface(np, gl, l, _al(gl,l,theta,+1), theta, &metrics);
#ifdef _RESTIME_CDF
  find_jacvars(np[l],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[_al(gl,l,theta,+1)],gl,metrics,theta,&jacvarsp1h);
  find_Deltaxt_interface(np, gl, l, _al(gl,l,theta,+1), theta, jacvarsm1h, jacvarsp1h, metrics, Deltaxt);
#else
  for (flux=0; flux<nf; flux++) Deltaxt[flux]=0.0;
#endif

  hbw_interpol=0;
  if (gl->cycle.resconv.INTERPOL==INTERPOL_FIRSTORDER) hbw_interpol=1;
  if (gl->cycle.resconv.INTERPOL==INTERPOL_BDF2 || gl->cycle.resconv.INTERPOL==INTERPOL_WENO3 || gl->cycle.resconv.INTERPOL==INTERPOL_TVD2_VANLEER || gl->cycle.resconv.INTERPOL==INTERPOL_TVD2_MINMOD || gl->cycle.resconv.INTERPOL==INTERPOL_TVD2_SUPERBEE) hbw_interpol=2;
  if (gl->cycle.resconv.INTERPOL==INTERPOL_TVD5 || gl->cycle.resconv.INTERPOL==INTERPOL_CWENO3 || gl->cycle.resconv.INTERPOL==INTERPOL_WENO5 || gl->cycle.resconv.INTERPOL==INTERPOL_AOWENO5) hbw_interpol=3;
  if (gl->cycle.resconv.INTERPOL==INTERPOL_WENO7 || gl->cycle.resconv.INTERPOL==INTERPOL_AOWENO7) hbw_interpol=4;
  if (hbw_interpol==0) fatal_error("Problem finding hbw_interpol in find_Fstar_interface");

  find_Fplus_Fminus_offset(np, gl, l, theta, +0, metrics,  Deltaxt, Fplusm1h, Fminusm1h);
  find_Fplus_Fminus_offset(np, gl, l, theta, +1, metrics,  Deltaxt, Fplusp1h, Fminusp1h);
  if (hbw_interpol>=2) {
    find_Fplus_Fminus_offset(np, gl, l, theta, -1, metrics,  Deltaxt, Fplusm3h, Fminusm3h);
    find_Fplus_Fminus_offset(np, gl, l, theta, +2, metrics,  Deltaxt, Fplusp3h, Fminusp3h);
  }
  if (hbw_interpol>=3) {
    find_Fplus_Fminus_offset(np, gl, l, theta, -2, metrics,  Deltaxt, Fplusm5h, Fminusm5h);
    find_Fplus_Fminus_offset(np, gl, l, theta, +3, metrics,  Deltaxt, Fplusp5h, Fminusp5h);
  }
  if (hbw_interpol>=4) {
    find_Fplus_Fminus_offset(np, gl, l, theta, -3, metrics,  Deltaxt, Fplusm7h, Fminusm7h);
    find_Fplus_Fminus_offset(np, gl, l, theta, +4, metrics,  Deltaxt, Fplusp7h, Fminusp7h);
  }

  nodes_from_bdry=_nodes_from_bdry_limited(np,gl,l,theta,TYPELEVEL_FLUID_WORK,hbw_resconv_fluid);
  AOWENO_TYPE=gl->cycle.resconv.AOWENO_TYPE;
  gammalo=gl->cycle.resconv.AOWENO_gammalo;
  gammahi=gl->cycle.resconv.AOWENO_gammahi;


  for (flux=0; flux<nf; flux++) {
    switch (gl->cycle.resconv.INTERPOL){
      case INTERPOL_FIRSTORDER:
        FplusL[flux]=Fplusm1h[flux];
        FminusR[flux]=Fminusp1h[flux];
      break;
#if (hbw_resconv_fluid>=2)
      case INTERPOL_BDF2:
        FplusL[flux]=-0.5*Fplusm3h[flux]+1.5*Fplusm1h[flux];
        FminusR[flux]=-0.5*Fminusp3h[flux]+1.5*Fminusp1h[flux];
      break;
      case INTERPOL_TVD2_MINMOD:
        FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_MINMOD);
        FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_MINMOD);
      break;
      case INTERPOL_TVD2_VANLEER:
        FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_VANLEER);
        FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_VANLEER);
      break;
      case INTERPOL_TVD2_SUPERBEE:
        FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_SUPERBEE);
        FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_SUPERBEE);
      break;
      case INTERPOL_WENO3:
        FplusL[flux]=_f_WENO3(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux]);
        FminusR[flux]=_f_WENO3(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux]);
      break;
#endif
#if (hbw_resconv_fluid>=3)
      case INTERPOL_TVD5:
        if (nodes_from_bdry>=3){
          FplusL[flux]=_f_TVD5(Fplusm5h[flux],Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],Fplusp3h[flux]);
          FminusR[flux]=_f_TVD5(Fminusp5h[flux],Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],Fminusm3h[flux]);
        } else {
          FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_VANLEER);
          FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_VANLEER);
        }
      break;
      case INTERPOL_CWENO3:
        if (nodes_from_bdry>=3){
          FplusL[flux]=_f_CWENO3(Fplusm5h[flux],Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],Fplusp3h[flux],0.5);
          FminusR[flux]=_f_CWENO3(Fminusp5h[flux],Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],Fminusm3h[flux],-0.5);
        } else {
          FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_VANLEER);
          FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_VANLEER);
        }
      break;
      case INTERPOL_WENO5:
        if (nodes_from_bdry>=3){
          FplusL[flux]=_f_WENO5(Fplusm5h[flux],Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],Fplusp3h[flux]);
          FminusR[flux]=_f_WENO5(Fminusp5h[flux],Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],Fminusm3h[flux]);
        } else {
          FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_VANLEER);
          FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_VANLEER);
        }
      break;
      case INTERPOL_AOWENO5:
        if (nodes_from_bdry>=3){
          FplusL[flux]=_f_AOWENO5(Fplusm5h[flux],Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],Fplusp3h[flux],0.5,gammalo,gammahi,AOWENO_TYPE);
          FminusR[flux]=_f_AOWENO5(Fminusp5h[flux],Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],Fminusm3h[flux],0.5,gammalo,gammahi,AOWENO_TYPE);
        } else {
          FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_VANLEER);
          FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_VANLEER);
        }
      break;
#endif
#if (hbw_resconv_fluid>=4)
      case INTERPOL_WENO7:
        if (nodes_from_bdry>=4){
          FplusL[flux]=_f_WENO7(Fplusm7h[flux],Fplusm5h[flux],Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],Fplusp3h[flux],Fplusp5h[flux]);
          FminusR[flux]=_f_WENO7(Fminusp7h[flux],Fminusp5h[flux],Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],Fminusm3h[flux],Fminusm5h[flux]);
        } else {
          FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_VANLEER);
          FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_VANLEER);
        }
      break;
      case INTERPOL_AOWENO7:
        if (nodes_from_bdry>=4){
          FplusL[flux]=_f_AOWENO7(Fplusm7h[flux],Fplusm5h[flux],Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],Fplusp3h[flux],Fplusp5h[flux],0.5,gammalo,gammahi,AOWENO_TYPE);
          FminusR[flux]=_f_AOWENO7(Fminusp7h[flux],Fminusp5h[flux],Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],Fminusm3h[flux],Fminusm5h[flux],0.5,gammalo,gammahi,AOWENO_TYPE);
        } else {
          FplusL[flux]=_f_TVD2(Fplusm3h[flux],Fplusm1h[flux],Fplusp1h[flux],LIMITER_VANLEER);
          FminusR[flux]=_f_TVD2(Fminusp3h[flux],Fminusp1h[flux],Fminusm1h[flux],LIMITER_VANLEER);
        }
      break;
#endif
      default:
        fatal_error("gl->cycle.resconv.INTERPOL invalid in find_Fstar_interface().");
    }
    Finttmp[flux]=FplusL[flux]+FminusR[flux];
  } 
  
  filter_Fstar_interface_positivity_preserving(np, gl, l, theta, metrics, gl->cycle.resconv.numiter, gl->cycle.resconv.EIGENVALCOND, Finttmp, Fint, lambdaminusp1h,  lambdaplusm1h);

}




void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  flux_t Fm1h,Fp1h;
  long flux,l;
  sqmat_t Lambdaminus_p0,Lambdaplus_p0,Lambdaminus_p1;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){

    if(l==ls){
      find_Fstar_interface(np, gl, _al(gl,l,theta,-1), theta,  Fm1h,Lambdaminus_p1,Lambdaplus_p0);
    } else {
      for (flux=0; flux<nf; flux++) {
        Fm1h[flux]=Fp1h[flux];
      }
    }

    for (flux=0; flux<nf; flux++) {
      Lambdaminus_p0[flux][flux]=Lambdaminus_p1[flux][flux];
    }

    find_Fstar_interface(np, gl, l, theta,  Fp1h,Lambdaminus_p1,Lambdaplus_p0);

#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
    for (flux=0; flux<nf; flux++) {
      np[l].bs->Res_trapezoidal[flux]+=gl->cycle.restime.weightm1_trapezoidal_convection*(Fp1h[flux]-Fm1h[flux]);
      np[l].wk->Res[flux]+=(1.0-gl->cycle.restime.weightm1_trapezoidal_convection)*(Fp1h[flux]-Fm1h[flux]);
    }
#else
    for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=(Fp1h[flux]-Fm1h[flux]);
#endif

    for (flux=0; flux<nf; flux++){
      np[l].bs->Delta_Lambda[theta][flux]=-Lambdaminus_p0[flux][flux]    +Lambdaplus_p0[flux][flux];
    }
  }

}


void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  long flux;
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=np[l].bs->Delta_Lambda[theta][flux];
}

