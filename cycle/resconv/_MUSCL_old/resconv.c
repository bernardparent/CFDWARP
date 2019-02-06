#include <cycle/resconv/_resconv.h>
#include <cycle/share/res_share.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>

#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif


#define FLUX_FDS 1
#define FLUX_FVS 2
#define INTERPOL_AOWENO5 1
#define INTERPOL_AOWENO7 2
#define INTERPOL_AOWENO9 3
#define INTERPOL_WENO3 4
#define INTERPOL_WENO5 5
#define INTERPOL_WENO7 6
#define INTERPOL_WENO9 7
#define INTERPOL_TVD2_MINMOD 8
#define INTERPOL_TVD2_VANLEER 9
#define INTERPOL_TVD2_SUPERBEE 10
#define INTERPOL_TVD2_SMART 16
#define INTERPOL_TVD5 11
#define INTERPOL_CWENO3 12
#define INTERPOL_CWENO5 13
#define INTERPOL_FIRSTORDER 14
#define INTERPOL_SECONDORDER 15
#define FACEINTEG_CENTRAL1 21
#define FACEINTEG_CENTRAL3 22
#define FACEINTEG_CENTRAL5 23
#define FACEINTEG_AOWENO5 24



void write_disc_resconv_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    FLUX=FLUX_FDS;\n"
    "    AVERAGING=AVERAGING_ROE;\n"
    "    AOWENO_TYPE=AOWENO_TYPE_DIFFUSIVE;\n"
    "    AOWENO_gammalo=0.95;\n"
    "    AOWENO_gammahi=0.999;\n"
    "    INTERPOL=INTERPOL_AOWENO7;\n"
    "    FACEINTEG=FACEINTEG_CENTRAL3;\n"
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
    SOAP_add_int_to_vars(codex,"AVERAGING_ROE",AVERAGING_ROE);
    SOAP_add_int_to_vars(codex,"AVERAGING_ARITH",AVERAGING_ARITH);
    SOAP_add_int_to_vars(codex,"FLUX_FDS",FLUX_FDS); 
    SOAP_add_int_to_vars(codex,"FLUX_FVS",FLUX_FVS); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_MINMOD",INTERPOL_TVD2_MINMOD); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_VANLEER",INTERPOL_TVD2_VANLEER); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_SUPERBEE",INTERPOL_TVD2_SUPERBEE); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_SMART",INTERPOL_TVD2_SMART); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD5",INTERPOL_TVD5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_AOWENO5",INTERPOL_AOWENO5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_AOWENO7",INTERPOL_AOWENO7); 
    SOAP_add_int_to_vars(codex,"INTERPOL_AOWENO9",INTERPOL_AOWENO9); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO3",INTERPOL_WENO3); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO5",INTERPOL_WENO5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO7",INTERPOL_WENO7); 
    SOAP_add_int_to_vars(codex,"INTERPOL_WENO9",INTERPOL_WENO9); 
    SOAP_add_int_to_vars(codex,"INTERPOL_CWENO3",INTERPOL_CWENO3); 
    SOAP_add_int_to_vars(codex,"INTERPOL_CWENO5",INTERPOL_CWENO5); 
    SOAP_add_int_to_vars(codex,"INTERPOL_FIRSTORDER",INTERPOL_FIRSTORDER); 
    SOAP_add_int_to_vars(codex,"INTERPOL_SECONDORDER",INTERPOL_SECONDORDER); 
    SOAP_add_int_to_vars(codex,"FACEINTEG_CENTRAL1",FACEINTEG_CENTRAL1); 
    SOAP_add_int_to_vars(codex,"FACEINTEG_CENTRAL3",FACEINTEG_CENTRAL3); 
    SOAP_add_int_to_vars(codex,"FACEINTEG_CENTRAL5",FACEINTEG_CENTRAL5); 
    SOAP_add_int_to_vars(codex,"FACEINTEG_AOWENO5",FACEINTEG_AOWENO5); 

    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_COMPRESSIVE",AOWENO_TYPE_COMPRESSIVE); 
    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_DIFFUSIVE",AOWENO_TYPE_DIFFUSIVE); 

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

    find_int_var_from_codex(codex,"FACEINTEG",&gl->cycle.resconv.FACEINTEG);
    if (gl->cycle.resconv.FACEINTEG!=FACEINTEG_CENTRAL1  && gl->cycle.resconv.FACEINTEG!=FACEINTEG_CENTRAL3 && gl->cycle.resconv.FACEINTEG!=FACEINTEG_CENTRAL5 && gl->cycle.resconv.FACEINTEG!=FACEINTEG_AOWENO5)
      SOAP_fatal_error(codex,"FACEINTEG must be set to either FACEINTEG_CENTRAL1 or FACEINTEG_CENTRAL3 or FACEINTEG_CENTRAL5 or FACEINTEG_AOWENO5.");

    find_int_var_from_codex(codex,"FLUX",&gl->cycle.resconv.FLUX);
    if (gl->cycle.resconv.FLUX!=FLUX_FDS && gl->cycle.resconv.FLUX!=FLUX_FVS)
      SOAP_fatal_error(codex,"FLUX must be set to either FLUX_FDS or FLUX_FVS.");
    find_int_var_from_codex(codex,"INTERPOL",&gl->cycle.resconv.INTERPOL);
    if (
#if (hbw_resconv_fluid>=2)
        gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_MINMOD && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_VANLEER 
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_SUPERBEE && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_SMART && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO3
     && gl->cycle.resconv.INTERPOL!=INTERPOL_FIRSTORDER && gl->cycle.resconv.INTERPOL!=INTERPOL_SECONDORDER
#endif
#if (hbw_resconv_fluid>=3)
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD5 && gl->cycle.resconv.INTERPOL!=INTERPOL_AOWENO5 
     && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO5 && gl->cycle.resconv.INTERPOL!=INTERPOL_CWENO3
#endif
#if (hbw_resconv_fluid>=4)      
     && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO7 && gl->cycle.resconv.INTERPOL!=INTERPOL_AOWENO7 
#endif
#if (hbw_resconv_fluid>=5)      
     && gl->cycle.resconv.INTERPOL!=INTERPOL_AOWENO9 && gl->cycle.resconv.INTERPOL!=INTERPOL_WENO9
     && gl->cycle.resconv.INTERPOL!=INTERPOL_CWENO5
#endif
     )
      SOAP_fatal_error(codex,"INTERPOL must be set to either "
#if (hbw_resconv_fluid>=2)      
        "INTERPOL_FIRSTORDER, INTERPOL_SECONDORDER, INTERPOL_TVD2_MINMOD, INTERPOL_TVD2_VANLEER, INTERPOL_TVD2_SUPERBEE, INTERPOL_TVD2_SMART, INTERPOL_WENO3, "
#endif
#if (hbw_resconv_fluid>=3)      
        "INTERPOL_AOWENO5, INTERPOL_WENO5, INTERPOL_CWENO3, INTERPOL_TVD5,"
#endif
#if (hbw_resconv_fluid>=4)      
        " INTERPOL_AOWENO7, INTERPOL_WENO7,"
#endif
#if (hbw_resconv_fluid>=5)      
        " INTERPOL_WENO9, INTERPOL_AOWENO9, INTERPOL_CWENO5."
#endif
      );
    find_int_var_from_codex(codex,"AOWENO_TYPE",&gl->cycle.resconv.AOWENO_TYPE);
    if (gl->cycle.resconv.AOWENO_TYPE!=AOWENO_TYPE_COMPRESSIVE && gl->cycle.resconv.AOWENO_TYPE!=AOWENO_TYPE_DIFFUSIVE)
      SOAP_fatal_error(codex,"AOWENO_TYPE must be set to either AOWENO_TYPE_COMPRESSIVE or AOWENO_TYPE_DIFFUSIVE.");

    find_double_var_from_codex(codex,"AOWENO_gammalo",&gl->cycle.resconv.AOWENO_gammalo);
    find_double_var_from_codex(codex,"AOWENO_gammahi",&gl->cycle.resconv.AOWENO_gammahi);

    gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FDS;
    if (gl->cycle.resconv.FLUX==FLUX_FVS) gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FVS;


    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}


static void find_Fstar_interface(np_t *np, gl_t *gl, long l, long theta, flux_t musclvarsm9h, flux_t musclvarsm7h, flux_t musclvarsm5h, flux_t musclvarsm3h, flux_t musclvarsm1h, flux_t musclvarsp1h, flux_t musclvarsp3h, flux_t musclvarsp5h, flux_t musclvarsp7h, flux_t musclvarsp9h, flux_t Fint){
  metrics_t metrics;
  long flux,nodes_from_bdry;
  flux_t musclvarsL,musclvarsR;
  int AOWENO_TYPE;
  double gammalo,gammahi;
#ifdef _RESTIME_CDF  
  sqmat_t lambdaminusp1h, lambdaplusm1h;
#endif

  find_metrics_at_interface(np, gl, l, _al(gl,l,theta,+1), theta, &metrics);
  nodes_from_bdry=_nodes_from_bdry_through_links_limited(np,gl,l,theta,TYPELEVEL_FLUID_WORK,hbw_resconv_fluid);
  for (flux=0; flux<nf; flux++) {
    AOWENO_TYPE=gl->cycle.resconv.AOWENO_TYPE;
    gammalo=gl->cycle.resconv.AOWENO_gammalo;
    gammahi=gl->cycle.resconv.AOWENO_gammahi;

    switch (gl->cycle.resconv.INTERPOL){
#if (hbw_resconv_fluid>=2)
      case INTERPOL_FIRSTORDER:
        musclvarsL[flux]=musclvarsm1h[flux];
        musclvarsR[flux]=musclvarsp1h[flux];
      break;
      case INTERPOL_SECONDORDER:
        musclvarsL[flux]=1.5*musclvarsm1h[flux]-0.5*musclvarsm3h[flux];
        musclvarsR[flux]=1.5*musclvarsp1h[flux]-0.5*musclvarsp3h[flux];
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
      case INTERPOL_TVD2_SMART:
        musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_SMART);
        musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_SMART);
      break;
      case INTERPOL_WENO3:
        musclvarsL[flux]=_f_WENO3(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux]);
        musclvarsR[flux]=_f_WENO3(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux]);
      break;
#endif

#if (hbw_resconv_fluid>=3)
      case INTERPOL_TVD5:
        if (nodes_from_bdry>=3){
          musclvarsL[flux]=_f_TVD5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
          musclvarsR[flux]=_f_TVD5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
        } else {
          musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_VANLEER);
          musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_VANLEER);
        } 
      break;
      case INTERPOL_WENO5:
        musclvarsL[flux]=_f_WENO5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
        musclvarsR[flux]=_f_WENO5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
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
      case INTERPOL_CWENO3:
        if (nodes_from_bdry>=3 ){
          musclvarsL[flux]=_f_CWENO3(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],0.5);
          musclvarsR[flux]=_f_CWENO3(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],-0.5);
        } else {
          musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_VANLEER);
          musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_VANLEER);
//          musclvarsL[flux]=_f_TVD5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
//          musclvarsR[flux]=_f_TVD5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
        } 
      break;
#endif

#if (hbw_resconv_fluid>=4)
      case INTERPOL_WENO7:
        musclvarsL[flux]=_f_WENO7(musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux]);
        musclvarsR[flux]=_f_WENO7(musclvarsp7h[flux],musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux],musclvarsm5h[flux]);
      break;
      case INTERPOL_AOWENO7:
        if (nodes_from_bdry>=4 ){
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
#endif

#if (hbw_resconv_fluid>=5)
      case INTERPOL_WENO9:
        musclvarsL[flux]=_f_WENO9(musclvarsm9h[flux],musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux]);
        musclvarsR[flux]=_f_WENO9(musclvarsp9h[flux],musclvarsp7h[flux],musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux],musclvarsm5h[flux],musclvarsm7h[flux]);
      break;
      case INTERPOL_CWENO5:
        if (nodes_from_bdry>=5 ){
          musclvarsL[flux]=_f_CWENO5(musclvarsm9h[flux],musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux],0.5);
          musclvarsR[flux]=_f_CWENO5(musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux],musclvarsp9h[flux],-0.5);
        } else {
          musclvarsL[flux]=_f_WENO9(musclvarsm9h[flux],musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux]);
          musclvarsR[flux]=_f_WENO9(musclvarsp9h[flux],musclvarsp7h[flux],musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux],musclvarsm5h[flux],musclvarsm7h[flux]);
        } 
      break;
      case INTERPOL_AOWENO9:
        if (nodes_from_bdry>=5 ){
          musclvarsL[flux]=_f_AOWENO9(musclvarsm9h[flux],musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux],0.5,gammalo,gammahi,AOWENO_TYPE);
          musclvarsR[flux]=_f_AOWENO9(musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux],musclvarsp9h[flux],-0.5,gammalo,gammahi,AOWENO_TYPE); 
        } else {
          if (nodes_from_bdry>=4 ){
            musclvarsL[flux]=_f_AOWENO7(musclvarsm7h[flux],musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],0.5,0.95,0.95,AOWENO_TYPE);
            musclvarsR[flux]=_f_AOWENO7(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],musclvarsp7h[flux],-0.5,0.95,0.95,AOWENO_TYPE); 
          } else {
            if (nodes_from_bdry>=3){
              musclvarsL[flux]=_f_AOWENO5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],0.5,0.95,0.95,AOWENO_TYPE);
              musclvarsR[flux]=_f_AOWENO5(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux],musclvarsp5h[flux],-0.5,0.95,0.95,AOWENO_TYPE); 
            } else {
              musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_VANLEER);
              musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_VANLEER);
//            musclvarsL[flux]=_f_TVD5(musclvarsm5h[flux],musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],musclvarsp3h[flux]);
//            musclvarsR[flux]=_f_TVD5(musclvarsp5h[flux],musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],musclvarsm3h[flux]);
            }
          } 
        } 
      break;
#endif

      default:
        fatal_error("gl->cycle.resconv.INTERPOL invalid in find_Fstar_interface().");
    } 

  }


  switch (gl->cycle.resconv.FLUX){
    case FLUX_FDS:
      find_Fstar_interface_FDS_muscl(np, gl, l, _al(gl,l,theta,+1), theta, musclvarsL, musclvarsR, metrics, gl->cycle.resconv.EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint);
    break;
    case FLUX_FVS:
      find_Fstar_interface_FVS_muscl(np, gl, l, _al(gl,l,theta,+1), theta, musclvarsL, musclvarsR, metrics, gl->cycle.resconv.EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint);
    break;
    default:
      fatal_error("gl->cycle.resconv.FLUX must be set to either FLUX_FDS or FLUX_FVS in find_Fstar_interface().");
  }


}




void find_Fstar_interfaces(np_t *np, gl_t *gl, long theta, long ls, long le){
  flux_t Fm1h,Fp1h;
  long flux,l;
  flux_t musclvarsm5,musclvarsm4,musclvarsm3,musclvarsm2,musclvarsm1,musclvarsp0,musclvarsp1,
         musclvarsp2,musclvarsp3,musclvarsp4,musclvarsp5;

  if (gl->cycle.resconv.FACEINTEG!=FACEINTEG_CENTRAL1){
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){

      if(l==ls){

        if (hbw_resconv_fluid>=5) find_musclvars_offset(np,gl,l,theta,-5,musclvarsm5);
        if (hbw_resconv_fluid>=4) find_musclvars_offset(np,gl,l,theta,-4,musclvarsm4);
        if (hbw_resconv_fluid>=3) find_musclvars_offset(np,gl,l,theta,-3,musclvarsm3);
        if (hbw_resconv_fluid>=2) find_musclvars_offset(np,gl,l,theta,-2,musclvarsm2);
        find_musclvars_offset(np,gl,l,theta,-1,musclvarsm1);
        find_musclvars_offset(np,gl,l,theta,+0,musclvarsp0);
        find_musclvars_offset(np,gl,l,theta,+1,musclvarsp1);
        if (hbw_resconv_fluid>=2) find_musclvars_offset(np,gl,l,theta,+2,musclvarsp2);
        if (hbw_resconv_fluid>=3) find_musclvars_offset(np,gl,l,theta,+3,musclvarsp3);
        if (hbw_resconv_fluid>=4) find_musclvars_offset(np,gl,l,theta,+4,musclvarsp4);
        if (hbw_resconv_fluid>=5) find_musclvars_offset(np,gl,l,theta,+5,musclvarsp5);
        find_Fstar_interface(np, gl, _al(gl,l,theta,-1), theta, musclvarsm5, musclvarsm4, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, musclvarsp4, Fm1h);

      } else {

        for (flux=0; flux<nf; flux++) {
          Fm1h[flux]=Fp1h[flux];
          if (hbw_resconv_fluid>=5) musclvarsm4[flux]=musclvarsm3[flux];
          if (hbw_resconv_fluid>=4) musclvarsm3[flux]=musclvarsm2[flux];
          if (hbw_resconv_fluid>=3) musclvarsm2[flux]=musclvarsm1[flux];
          if (hbw_resconv_fluid>=2) musclvarsm1[flux]=musclvarsp0[flux];
          musclvarsp0[flux]=musclvarsp1[flux];
          if (hbw_resconv_fluid>=2) musclvarsp1[flux]=musclvarsp2[flux];
          if (hbw_resconv_fluid>=3) musclvarsp2[flux]=musclvarsp3[flux];
          if (hbw_resconv_fluid>=4) musclvarsp3[flux]=musclvarsp4[flux];
          if (hbw_resconv_fluid>=5) musclvarsp4[flux]=musclvarsp5[flux];
        }

      }

      switch (hbw_resconv_fluid){
        case 1:
          find_musclvars_offset(np,gl,l,theta,+1,musclvarsp1);
        break;
        case 2:
          find_musclvars_offset(np,gl,l,theta,+2,musclvarsp2);
        break;
        case 3:
          find_musclvars_offset(np,gl,l,theta,+3,musclvarsp3);
        break;
        case 4:
          find_musclvars_offset(np,gl,l,theta,+4,musclvarsp4);
        break;
        case 5:
          find_musclvars_offset(np,gl,l,theta,+5,musclvarsp5);
        break;
        default:
          fatal_error("hbw_resconv_fluid can not be set to %ld",hbw_resconv_fluid);
      }
      find_Fstar_interface(np, gl, l, theta, musclvarsm4, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, musclvarsp4, musclvarsp5, Fp1h);
      np[l].wk->Fp1h=(double *)realloc(np[l].wk->Fp1h,nf*sizeof(double));
      np[_al(gl,l,theta,-1)].wk->Fp1h=(double *)realloc(np[_al(gl,l,theta,-1)].wk->Fp1h,nf*sizeof(double));
      for (flux=0; flux<nf; flux++) {
        assert(is_node_valid(np[l],TYPELEVEL_FLUID));
        np[l].wk->Fp1h[flux]=Fp1h[flux];
        assert(is_node_valid(np[_al(gl,l,theta,-1)],TYPELEVEL_FLUID));
        np[_al(gl,l,theta,-1)].wk->Fp1h[flux]=Fm1h[flux];

      }
    }
  }

}


void init_Fstar_interfaces(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;
  zone_t newzone;
  newzone=_zone_intersection(gl->domain_all,_zone_expansion(zone,+1));
  for1DL(i,newzone.is,newzone.ie)
    for2DL(j,newzone.js,newzone.je)
      for3DL(k,newzone.ks,newzone.ke)
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)){
          free(np[_ai(gl,i,j,k)].wk->Fp1h);
          np[_ai(gl,i,j,k)].wk->Fp1h=NULL;
        }
      end3DL
    end2DL
  end1DL
}


int find_face_integrated_flux_CENTRAL5(np_t *np, gl_t *gl, long l, long dim, flux_t F){
  int error;
  long flux;
  long lm2,lm1,lp0,lp1,lp2;
  lm2=_al(gl,l,dim,-2);
  lm1=_al(gl,l,dim,-1);
  lp0=_al(gl,l,dim,+0);
  lp1=_al(gl,l,dim,+1);
  lp2=_al(gl,l,dim,+2);
  error=4;
  if (is_node_valid(np[lp0],TYPELEVEL_FLUID_WORK) && np[lp0].wk->Fp1h!=NULL){
    if (is_node_valid(np[lm2],TYPELEVEL_FLUID_WORK) && is_node_valid(np[lp2],TYPELEVEL_FLUID_WORK) && np[lm2].wk->Fp1h!=NULL && np[lm1].wk->Fp1h!=NULL && np[lp0].wk->Fp1h!=NULL && np[lp1].wk->Fp1h!=NULL && np[lp2].wk->Fp1h!=NULL){
      error=0;
      for (flux=0; flux<nf; flux++) F[flux]=(23.0*np[lm2].wk->Fp1h[flux]-332.0*np[lm1].wk->Fp1h[flux]+6378.0*np[lp0].wk->Fp1h[flux]-332.0*np[lp1].wk->Fp1h[flux]+23.0*np[lp2].wk->Fp1h[flux])/5760.0;
    } else {
      if (is_node_valid(np[lm2],TYPELEVEL_FLUID_WORK) && is_node_valid(np[lp2],TYPELEVEL_FLUID_WORK) && np[lm1].wk->Fp1h!=NULL && np[lp0].wk->Fp1h!=NULL && np[lp1].wk->Fp1h!=NULL){
        error=1;
        for (flux=0; flux<nf; flux++) F[flux]=(np[lm1].wk->Fp1h[flux]+22.0*np[lp0].wk->Fp1h[flux]+np[lp1].wk->Fp1h[flux])/24.0;
      } else {
        if (np[lp0].wk->Fp1h!=NULL){
          error=2;
          for (flux=0; flux<nf; flux++) F[flux]=np[lp0].wk->Fp1h[flux];
        } else {
          error=3;
        }
      }
    }
  }
  return(error);
}


double integrated_flux_AOWENO5(double Fm2, double Fm1, double Fp0, double Fp1, double Fp2){
  double ret;
  double F3,F5,epsilon;
  double gammahi,w3,w5,wtil3,wtil5,c31,c32,beta3,beta5,c51,c52,c53,c54;
  gammahi=0.9;
  epsilon=1e-12;
  F3=(Fm1+22.0*Fp0+Fp1)/24.0;
  F5=(23.0*Fm2-332.0*Fm1+6378.0*Fp0-332.0*Fp1+23.0*Fp2)/5760.0;  
  c31=(Fp1+Fm1)/2.0-Fp0;
  c32=(Fp1+Fm1)/2.0-Fm1;
  beta3=13.0/3.0*sqr(c31)+sqr(c32);
  c51=(Fm2-4.0*Fm1+6.0*Fp0-4.0*Fp1+Fp2)/24.0;
  c52=(Fp2-Fm2-2.0*Fp1+2.0*Fm1)/12.0;
  c53=(Fm2-16.0*Fm1+30.0*Fp0-16.0*Fp1+Fp2)/24.0;
  c54=(-8.0*Fp1+8.0*Fm1+Fp2-Fm2)/(-12.0);
  beta5=87617.0/140.0*sqr(c51)+3129.0/80.0*sqr(c52)+21.0/5.0*c51*c53+13.0/3.0*sqr(c53)+0.5*c52*c54+sqr(c54);
  w3=(1.0-gammahi)/sqr(beta3+epsilon);
  w5=gammahi/sqr(beta5+epsilon);
  wtil3=w3/(w3+w5);
  wtil5=w5/(w3+w5);
  ret=wtil5/gammahi*(F5-(1.0-gammahi)*F3)+wtil3*F3;
  
  return(ret);
}


int find_face_integrated_flux_AOWENO5(np_t *np, gl_t *gl, long l, long dim, flux_t F){
  int error;
  long flux;
  long lm2,lm1,lp0,lp1,lp2;
  double Fm2,Fm1,Fp0,Fp1,Fp2;
  lm2=_al(gl,l,dim,-2);
  lm1=_al(gl,l,dim,-1);
  lp0=_al(gl,l,dim,+0);
  lp1=_al(gl,l,dim,+1);
  lp2=_al(gl,l,dim,+2);
  error=4;
  if (is_node_valid(np[lp0],TYPELEVEL_FLUID_WORK) && np[lp0].wk->Fp1h!=NULL){
    if (is_node_valid(np[lm2],TYPELEVEL_FLUID_WORK) && is_node_valid(np[lp2],TYPELEVEL_FLUID_WORK) && np[lm2].wk->Fp1h!=NULL && np[lm1].wk->Fp1h!=NULL && np[lp0].wk->Fp1h!=NULL && np[lp1].wk->Fp1h!=NULL && np[lp2].wk->Fp1h!=NULL){
      error=0;
      for (flux=0; flux<nf; flux++) {
        Fm2=np[lm2].wk->Fp1h[flux];
        Fm1=np[lm1].wk->Fp1h[flux];
        Fp0=np[lp0].wk->Fp1h[flux];
        Fp1=np[lp1].wk->Fp1h[flux];
        Fp2=np[lp2].wk->Fp1h[flux];
        F[flux]=integrated_flux_AOWENO5(Fm2, Fm1, Fp0, Fp1, Fp2);
      }
    } else {
      if (is_node_valid(np[lm2],TYPELEVEL_FLUID_WORK) && is_node_valid(np[lp2],TYPELEVEL_FLUID_WORK) && np[lm1].wk->Fp1h!=NULL && np[lp0].wk->Fp1h!=NULL && np[lp1].wk->Fp1h!=NULL){
        error=1;
        for (flux=0; flux<nf; flux++) F[flux]=(np[lm1].wk->Fp1h[flux]+22.0*np[lp0].wk->Fp1h[flux]+np[lp1].wk->Fp1h[flux])/24.0;
      } else {
        if (np[lp0].wk->Fp1h!=NULL){
          error=2;
          for (flux=0; flux<nf; flux++) F[flux]=np[lp0].wk->Fp1h[flux];
        } else {
          error=3;
        }
      }
    }
  }
  return(error);
}




static void integrate_Fstar_interface(np_t *np, gl_t *gl, long l, long theta, flux_t Fint){
  long dim,flux;
  flux_t Fp0,Fm1,Fp1;
#ifdef _3D
  long i,j;
  flux_t Fm2,Fp2;
  int error_im1,error_im2,error_ip0,error_ip1,error_ip2;
  long dim2;
  long lp1p0,lm1p0,lp0p1,lp0m1,lp1p1,lp1m1,lm1p1,lm1m1;
#endif
  assert(is_node_valid(np[l],TYPELEVEL_FLUID));
  switch (gl->cycle.resconv.FACEINTEG){
    case FACEINTEG_CENTRAL1:
      for (flux=0; flux<nf; flux++) Fint[flux]=np[l].wk->Fp1h[flux];
    break;
    case FACEINTEG_CENTRAL3:
#ifdef _2D
      dim=mod(theta+1,nd);
      if (np[_al(gl,l,dim,+1)].wk->Fp1h!=NULL && np[_al(gl,l,dim,-1)].wk->Fp1h!=NULL){
        for (flux=0; flux<nf; flux++) {
          Fp0[flux]=np[l].wk->Fp1h[flux];
          Fp1[flux]=np[_al(gl,l,dim,+1)].wk->Fp1h[flux];
          Fm1[flux]=np[_al(gl,l,dim,-1)].wk->Fp1h[flux];
          Fint[flux]=(22.0*Fp0[flux]+Fp1[flux]+Fm1[flux])/24.0;
        }
      } else {
        for (flux=0; flux<nf; flux++) Fint[flux]=np[l].wk->Fp1h[flux];
      }
#endif
#ifdef _3D
      dim=mod(theta+1,nd);
      dim2=mod(theta+2,nd);
      assert(np[l].wk->Fp1h!=NULL);
      lp1p0=_al(gl,l,dim,+1);
      lm1p0=_al(gl,l,dim,-1);
      lp0p1=_al(gl,l,dim2,+1);
      lp0m1=_al(gl,l,dim2,-1);
      lp1p1=_all(gl,l,dim2,+1,dim,+1);
      lp1m1=_all(gl,l,dim2,-1,dim,+1);
      lm1p1=_all(gl,l,dim2,+1,dim,-1);
      lm1m1=_all(gl,l,dim2,-1,dim,-1);
      if (np[lp1p0].wk->Fp1h!=NULL && np[lm1p0].wk->Fp1h!=NULL){
        for (flux=0; flux<nf; flux++) {
        if (np[lp0p1].wk->Fp1h!=NULL && np[lp0m1].wk->Fp1h!=NULL) {
          Fp0[flux]=(22.0*np[l].wk->Fp1h[flux]+np[lp0p1].wk->Fp1h[flux]+np[lp0m1].wk->Fp1h[flux])/24.0;
        } else {
          Fp0[flux]=np[l].wk->Fp1h[flux];
        }
        if (np[lp1p1].wk->Fp1h!=NULL && np[lp1m1].wk->Fp1h!=NULL) {
          Fp1[flux]=(22.0*np[lp1p0].wk->Fp1h[flux]+np[lp1p1].wk->Fp1h[flux]+np[lp1m1].wk->Fp1h[flux])/24.0;
        } else {
          Fp1[flux]=np[lp1p0].wk->Fp1h[flux];
        }
        if (np[lm1p1].wk->Fp1h!=NULL && np[lm1m1].wk->Fp1h!=NULL) {
          Fm1[flux]=(22.0*np[lm1p0].wk->Fp1h[flux]+np[lm1p1].wk->Fp1h[flux]+np[lm1m1].wk->Fp1h[flux])/24.0;
        } else {
          Fm1[flux]=np[lm1p0].wk->Fp1h[flux];
        }
        Fint[flux]=(22.0*Fp0[flux]+Fp1[flux]+Fm1[flux])/24.0;
        }
      } else {
        for (flux=0; flux<nf; flux++) Fint[flux]=np[l].wk->Fp1h[flux];
      }
#endif
    break;
    case FACEINTEG_CENTRAL5:
#ifdef _2D
      dim=mod(theta+1,nd);
      find_face_integrated_flux_CENTRAL5(np,gl,l,dim,Fint);
#endif
#ifdef _3D
      i=mod(theta+1,nd);
      j=mod(theta+2,nd);
      assert(np[l].wk->Fp1h!=NULL);
       
      error_im2=find_face_integrated_flux_CENTRAL5(np,gl,_al(gl,l,i,-2),j,Fm2);
      error_ip2=find_face_integrated_flux_CENTRAL5(np,gl,_al(gl,l,i,+2),j,Fp2);  
      error_im1=find_face_integrated_flux_CENTRAL5(np,gl,_al(gl,l,i,-1),j,Fm1);
      error_ip1=find_face_integrated_flux_CENTRAL5(np,gl,_al(gl,l,i,+1),j,Fp1);
      error_ip0=find_face_integrated_flux_CENTRAL5(np,gl,_al(gl,l,i,+0),j,Fp0);
      
      if (error_im2==0 && error_ip2==0 && error_im1==0 && error_ip1==0 && error_ip0==0){
        for (flux=0; flux<nf; flux++) Fint[flux]=(23.0*Fm2[flux]-332.0*Fm1[flux]+6378.0*Fp0[flux]-332.0*Fp1[flux]+23.0*Fp2[flux])/5760.0;
      } else {
        if (error_im1==0 && error_ip1==0 && error_ip0==0){
          for (flux=0; flux<nf; flux++) Fint[flux]=(22.0*Fp0[flux]+Fp1[flux]+Fm1[flux])/24.0;
        } else {
          for (flux=0; flux<nf; flux++) Fint[flux]=np[l].wk->Fp1h[flux];
        }
      }

#endif

    break;
    case FACEINTEG_AOWENO5:
#ifdef _2D
      dim=mod(theta+1,nd);
      find_face_integrated_flux_AOWENO5(np,gl,l,dim,Fint);
#endif
#ifdef _3D
      i=mod(theta+1,nd);
      j=mod(theta+2,nd);
      assert(np[l].wk->Fp1h!=NULL);
       
      error_im2=find_face_integrated_flux_AOWENO5(np,gl,_al(gl,l,i,-2),j,Fm2);
      error_ip2=find_face_integrated_flux_AOWENO5(np,gl,_al(gl,l,i,+2),j,Fp2);  
      error_im1=find_face_integrated_flux_AOWENO5(np,gl,_al(gl,l,i,-1),j,Fm1);
      error_ip1=find_face_integrated_flux_AOWENO5(np,gl,_al(gl,l,i,+1),j,Fp1);
      error_ip0=find_face_integrated_flux_AOWENO5(np,gl,_al(gl,l,i,+0),j,Fp0);
      
      if (error_im2==0 && error_ip2==0 && error_im1==0 && error_ip1==0 && error_ip0==0){
        for (flux=0; flux<nf; flux++) Fint[flux]=integrated_flux_AOWENO5(Fm2[flux], Fm1[flux], Fp0[flux], Fp1[flux], Fp2[flux]);
      } else {
        if (error_im1==0 && error_ip1==0 && error_ip0==0){
          for (flux=0; flux<nf; flux++) Fint[flux]=(22.0*Fp0[flux]+Fp1[flux]+Fm1[flux])/24.0;
        } else {
          for (flux=0; flux<nf; flux++) Fint[flux]=np[l].wk->Fp1h[flux];
        }
      }

#endif

    break;
    default:
      fatal_error("FACEINTEG can not be set to %d.",gl->cycle.resconv.FACEINTEG);
  }
}


void add_dFstar_residual_new(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  flux_t Fm1h,Fp1h;
  long flux,l;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){

    if(l==ls){
      integrate_Fstar_interface(np, gl, _al(gl,l,theta,-1), theta, Fm1h);
    } else {
      for (flux=0; flux<nf; flux++) {
        Fm1h[flux]=Fp1h[flux];
      }
    }
    integrate_Fstar_interface(np, gl, l, theta, Fp1h);

    for (flux=0; flux<nf; flux++) {
      np[l].wk->Res[flux]+=fact*(Fp1h[flux]-Fm1h[flux]);
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
      np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*(Fp1h[flux]-Fm1h[flux]); 
#endif
    }
  }

}



void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  flux_t Fm1h,Fp1h;
  long flux,l;
  flux_t musclvarsm5,musclvarsm4,musclvarsm3,musclvarsm2,musclvarsm1,musclvarsp0,musclvarsp1,
         musclvarsp2,musclvarsp3,musclvarsp4,musclvarsp5;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){

    switch (gl->cycle.resconv.FACEINTEG){


      case FACEINTEG_CENTRAL1:
        if(l==ls){
          if (hbw_resconv_fluid>=5) find_musclvars_offset(np,gl,l,theta,-5,musclvarsm5);
          if (hbw_resconv_fluid>=4) find_musclvars_offset(np,gl,l,theta,-4,musclvarsm4);
          if (hbw_resconv_fluid>=3) find_musclvars_offset(np,gl,l,theta,-3,musclvarsm3);
          if (hbw_resconv_fluid>=2) find_musclvars_offset(np,gl,l,theta,-2,musclvarsm2);
          find_musclvars_offset(np,gl,l,theta,-1,musclvarsm1);
          find_musclvars_offset(np,gl,l,theta,+0,musclvarsp0);
          find_musclvars_offset(np,gl,l,theta,+1,musclvarsp1);
          if (hbw_resconv_fluid>=2) find_musclvars_offset(np,gl,l,theta,+2,musclvarsp2);
          if (hbw_resconv_fluid>=3) find_musclvars_offset(np,gl,l,theta,+3,musclvarsp3);
          if (hbw_resconv_fluid>=4) find_musclvars_offset(np,gl,l,theta,+4,musclvarsp4);
          if (hbw_resconv_fluid>=5) find_musclvars_offset(np,gl,l,theta,+5,musclvarsp5);
          find_Fstar_interface(np, gl, _al(gl,l,theta,-1), theta, musclvarsm5, musclvarsm4, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, musclvarsp4, Fm1h);

        } else {

          for (flux=0; flux<nf; flux++) {
            Fm1h[flux]=Fp1h[flux];
            if (hbw_resconv_fluid>=5) musclvarsm4[flux]=musclvarsm3[flux];
            if (hbw_resconv_fluid>=4) musclvarsm3[flux]=musclvarsm2[flux];
            if (hbw_resconv_fluid>=3) musclvarsm2[flux]=musclvarsm1[flux];
            if (hbw_resconv_fluid>=2) musclvarsm1[flux]=musclvarsp0[flux];
            musclvarsp0[flux]=musclvarsp1[flux];
            if (hbw_resconv_fluid>=2) musclvarsp1[flux]=musclvarsp2[flux];
            if (hbw_resconv_fluid>=3) musclvarsp2[flux]=musclvarsp3[flux];
            if (hbw_resconv_fluid>=4) musclvarsp3[flux]=musclvarsp4[flux];
            if (hbw_resconv_fluid>=5) musclvarsp4[flux]=musclvarsp5[flux];
          }

        }

        switch (hbw_resconv_fluid){
          case 1:
            find_musclvars_offset(np,gl,l,theta,+1,musclvarsp1);
          break;
          case 2:
            find_musclvars_offset(np,gl,l,theta,+2,musclvarsp2);
          break;
          case 3:
            find_musclvars_offset(np,gl,l,theta,+3,musclvarsp3);
          break;
          case 4:
            find_musclvars_offset(np,gl,l,theta,+4,musclvarsp4);
          break;
          case 5:
            find_musclvars_offset(np,gl,l,theta,+5,musclvarsp5);
          break;
          default:
            fatal_error("hbw_resconv_fluid can not be set to %ld",hbw_resconv_fluid);
        }
        find_Fstar_interface(np, gl, l, theta, musclvarsm4, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, musclvarsp4, musclvarsp5, Fp1h);
      break;

      default:
        if(l==ls){
          integrate_Fstar_interface(np, gl, _al(gl,l,theta,-1), theta, Fm1h);
        } else {
          for (flux=0; flux<nf; flux++) {
            Fm1h[flux]=Fp1h[flux];
          }
        }
        integrate_Fstar_interface(np, gl, l, theta, Fp1h);

    }

    for (flux=0; flux<nf; flux++) {
      np[l].wk->Res[flux]+=fact*(Fp1h[flux]-Fm1h[flux]);
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
      np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*(Fp1h[flux]-Fm1h[flux]); 
#endif
    }
  }

}







void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  sqmat_t Lambda_minus,Lambda_plus;
  long flux;
  find_Lambda_minus_dtau_FDS(np,gl,l,theta,gl->cycle.resconv.EIGENVALCOND,gl->cycle.resconv.AVERAGING,Lambda_minus);
  find_Lambda_plus_dtau_FDS(np,gl,l,theta,gl->cycle.resconv.EIGENVALCOND,gl->cycle.resconv.AVERAGING,Lambda_plus);
  for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=Lambda_plus[flux][flux]-Lambda_minus[flux][flux];
}

