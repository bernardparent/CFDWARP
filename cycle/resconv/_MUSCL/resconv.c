// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015-2018, 2024 Bernard Parent

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
#include <cycle/share/res_share.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>

#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif


#define FLUX_FDS 1
#define FLUX_FVS 2
#define FLUX_FDSR 3
#define FLUX_KNP 4
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
#define INTERPOL_TVD2_MINMOD2 17
#define INTERPOL_TVD5 11
#define INTERPOL_CWENO3 12
#define INTERPOL_CWENO5 13
#define INTERPOL_FIRSTORDER 14
#define INTERPOL_SECONDORDER 15
#define POSFILTER_NONE 1
#define POSFILTER_PARENT 2
#define POSFILTER_TEST 3



void write_disc_resconv_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    FLUX=FLUX_FDS;\n"
    "    AVERAGING=AVERAGING_ROE;\n"
    "    AOWENO_TYPE=AOWENO_TYPE_DIFFUSIVE;\n"
    "    AOWENO_gammalo=0.85;\n"
    "    AOWENO_gammahi=0.85;\n"
    "    INTERPOL=INTERPOL_TVD2_VANLEER;\n"
    "    EIGENVALCOND=EIGENVALCOND_PECLET;\n"
    "    POSFILTER=POSFILTER_NONE;\n"
    "    POSFILTER_numiter=4;\n"
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
    SOAP_add_int_to_vars(codex,"EIGENVALCOND_PARENT",EIGENVALCOND_PARENT);
    SOAP_add_int_to_vars(codex,"AVERAGING_ROE",AVERAGING_ROE);
    SOAP_add_int_to_vars(codex,"AVERAGING_ARITH",AVERAGING_ARITH);
    SOAP_add_int_to_vars(codex,"FLUX_FDS",FLUX_FDS); 
    SOAP_add_int_to_vars(codex,"FLUX_FDSR",FLUX_FDSR); 
    SOAP_add_int_to_vars(codex,"FLUX_FVS",FLUX_FVS); 
    SOAP_add_int_to_vars(codex,"FLUX_KNP",FLUX_KNP); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_MINMOD",INTERPOL_TVD2_MINMOD); 
    SOAP_add_int_to_vars(codex,"INTERPOL_TVD2_MINMOD2",INTERPOL_TVD2_MINMOD2); 
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
    SOAP_add_int_to_vars(codex,"POSFILTER_PARENT",POSFILTER_PARENT); 
    SOAP_add_int_to_vars(codex,"POSFILTER_TEST",POSFILTER_TEST); 
    SOAP_add_int_to_vars(codex,"POSFILTER_NONE",POSFILTER_NONE); 

    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_COMPRESSIVE",AOWENO_TYPE_COMPRESSIVE); 
    SOAP_add_int_to_vars(codex,"AOWENO_TYPE_DIFFUSIVE",AOWENO_TYPE_DIFFUSIVE); 

    gl->DISC_RESCONV_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_resconv_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"EIGENVALCOND",&gl->cycle.resconv.EIGENVALCOND);
    if (gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PECLET && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_HARTEN && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_GNOFFO && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PASCAL && gl->cycle.resconv.EIGENVALCOND!=EIGENVALCOND_PARENT)
      SOAP_fatal_error(codex,"EIGENVALCOND must be set to either EIGENVALCOND_PECLET, EIGENVALCOND_HARTEN, EIGENVALCOND_GNOFFO, EIGENVALCOND_PASCAL, EIGENVALCOND_PARENT.");

    find_int_var_from_codex(codex,"AVERAGING",&gl->cycle.resconv.AVERAGING);
    if (gl->cycle.resconv.AVERAGING!=AVERAGING_ROE  && gl->cycle.resconv.AVERAGING!=AVERAGING_ARITH)
      SOAP_fatal_error(codex,"AVERAGING must be set to either AVERAGING_ROE or AVERAGING_ARITH.");


    find_int_var_from_codex(codex,"POSFILTER",&gl->cycle.resconv.POSFILTER);
    if (gl->cycle.resconv.POSFILTER!=POSFILTER_NONE  && gl->cycle.resconv.POSFILTER!=POSFILTER_PARENT && gl->cycle.resconv.POSFILTER!=POSFILTER_TEST)
      SOAP_fatal_error(codex,"POSFILTER must be set to either POSFILTER_NONE or POSFILTER_PARENT or POSFILTER_TEST.");
    find_int_var_from_codex(codex,"POSFILTER_numiter",&gl->cycle.resconv.POSFILTER_numiter);
    if (gl->cycle.resconv.POSFILTER_numiter<1) SOAP_fatal_error(codex,"POSFILTER_numiter must be set to an integer greater or equal to 1.");


    find_int_var_from_codex(codex,"FLUX",&gl->cycle.resconv.FLUX);
    if (gl->cycle.resconv.FLUX!=FLUX_FDS && gl->cycle.resconv.FLUX!=FLUX_FVS && gl->cycle.resconv.FLUX!=FLUX_FDSR && gl->cycle.resconv.FLUX!=FLUX_KNP)
      SOAP_fatal_error(codex,"FLUX must be set to either FLUX_FDS or FLUX_FDSR or FLUX_FVS or FLUX_KNP.");
    find_int_var_from_codex(codex,"INTERPOL",&gl->cycle.resconv.INTERPOL);
    if (
#if (hbw_resconv_fluid>=2)
        gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_MINMOD && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_MINMOD2 
     && gl->cycle.resconv.INTERPOL!=INTERPOL_TVD2_VANLEER 
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
        "INTERPOL_FIRSTORDER, INTERPOL_SECONDORDER, INTERPOL_TVD2_MINMOD, INTERPOL_TVD2_MINMOD2, INTERPOL_TVD2_VANLEER, INTERPOL_TVD2_SUPERBEE, INTERPOL_TVD2_SMART, INTERPOL_WENO3, "
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
    if (gl->cycle.resconv.FLUX==FLUX_FDSR) gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FDSR;
    if (gl->cycle.resconv.FLUX==FLUX_FDS && gl->cycle.resconv.POSFILTER!=POSFILTER_NONE) gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FDSPLUS;
    if (gl->cycle.resconv.FLUX==FLUX_FVS && gl->cycle.resconv.POSFILTER!=POSFILTER_NONE) gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FVSPLUS;
    if (gl->cycle.resconv.FLUX==FLUX_KNP) gl->cycle.resconv.CONVJACOBIAN=CONVJACOBIAN_FVS;

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}


static void find_Fstar_interface(np_t *np, gl_t *gl, long l, long theta, flux_t musclvarsm9h, flux_t musclvarsm7h, flux_t musclvarsm5h, flux_t musclvarsm3h, flux_t musclvarsm1h, flux_t musclvarsp1h, flux_t musclvarsp3h, flux_t musclvarsp5h, flux_t musclvarsp7h, flux_t musclvarsp9h, flux_t Fint){
  metrics_t metrics;
  long flux,nodes_from_bdry;
  flux_t musclvarsL,musclvarsR;
  int AOWENO_TYPE,EIGENVALCOND;
  double gammalo,gammahi;

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
      case INTERPOL_TVD2_MINMOD2:
        musclvarsL[flux]=_f_TVD2(musclvarsm3h[flux],musclvarsm1h[flux],musclvarsp1h[flux],LIMITER_MINMOD2);
        musclvarsR[flux]=_f_TVD2(musclvarsp3h[flux],musclvarsp1h[flux],musclvarsm1h[flux],LIMITER_MINMOD2);
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


  if (gl->cycle.resconv.POSFILTER==POSFILTER_NONE){
    EIGENVALCOND=gl->cycle.resconv.EIGENVALCOND;
  } else {
    EIGENVALCOND=EIGENVALCOND_NONE;
  }
  switch (gl->cycle.resconv.FLUX){
    case FLUX_FDS:
      find_Fstar_interface_FDS_muscl(np, gl, l, _al(gl,l,theta,+1),  theta, musclvarsL, musclvarsR,
                     metrics, EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint);
    break;
    case FLUX_FDSR:
      find_Fstar_interface_FDSR_muscl(np, gl, l, _al(gl,l,theta,+1),  theta, musclvarsL, musclvarsR,
                     metrics, EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint);
    break;
    case FLUX_FVS:
      find_Fstar_interface_FVS_muscl(np, gl, l, _al(gl,l,theta,+1),  theta, musclvarsL,  musclvarsR,
                     metrics, EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint);
    break;
    case FLUX_KNP:
      find_Fstar_interface_KNP_muscl(gl, theta, musclvarsL, musclvarsR,
                     metrics, EIGENVALCOND, gl->cycle.resconv.AVERAGING, Fint);
    break;
    default:
      fatal_error("gl->cycle.resconv.FLUX must be set to either FLUX_FDS or FLUX_FDSR or FLUX_FVS or FLUX_KNP in find_Fstar_interface().");
  }


}


static void find_Fstar_interface_trapezoidal(np_t *np, gl_t *gl, long l, long theta, musclvarscycle_t musclvarsm9h, musclvarscycle_t musclvarsm7h, musclvarscycle_t musclvarsm5h, musclvarscycle_t musclvarsm3h, musclvarscycle_t musclvarsm1h, musclvarscycle_t musclvarsp1h, musclvarscycle_t musclvarsp3h, musclvarscycle_t musclvarsp5h, musclvarscycle_t musclvarsp7h, musclvarscycle_t musclvarsp9h, flux_t Fint){
  flux_t Fintp0,Fintm1;
  long flux;
  find_Fstar_interface(np, gl, l, theta, musclvarsm9h, musclvarsm7h, musclvarsm5h, musclvarsm3h, musclvarsm1h, musclvarsp1h, musclvarsp3h, musclvarsp5h, musclvarsp7h, musclvarsp9h, Fintp0);
#ifdef _RESTIME_TRAPEZOIDAL_MUSCL
  find_Fstar_interface(np, gl, l, theta, &(musclvarsm9h[nf]), &(musclvarsm7h[nf]), &(musclvarsm5h[nf]), &(musclvarsm3h[nf]), &(musclvarsm1h[nf]), &(musclvarsp1h[nf]), &(musclvarsp3h[nf]), &(musclvarsp5h[nf]), &(musclvarsp7h[nf]), &(musclvarsp9h[nf]), Fintm1);
#else
  for (flux=0; flux<nf; flux++) Fintm1[flux]=Fintp0[flux];
#endif
  for (flux=0; flux<nf; flux++) Fint[flux]=0.5*(Fintp0[flux]+Fintm1[flux]);
}





void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  flux_t Fm1h,Fp1h,Ftmp;
  long flux,l;
  metrics_t metrics;
  sqmat_t Lambdaminus_p0,Lambdaplus_p0,Lambdaminus_p1;
  musclvarscycle_t musclvarsm5,musclvarsm4,musclvarsm3,musclvarsm2,musclvarsm1,musclvarsp0,musclvarsp1,
         musclvarsp2,musclvarsp3,musclvarsp4,musclvarsp5;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){

        if(l==ls){
          if (hbw_resconv_fluid>=5) find_musclvarscycle_offset(np,gl,l,theta,-5,musclvarsm5);
          if (hbw_resconv_fluid>=4) find_musclvarscycle_offset(np,gl,l,theta,-4,musclvarsm4);
          if (hbw_resconv_fluid>=3) find_musclvarscycle_offset(np,gl,l,theta,-3,musclvarsm3);
          if (hbw_resconv_fluid>=2) find_musclvarscycle_offset(np,gl,l,theta,-2,musclvarsm2);
          find_musclvarscycle_offset(np,gl,l,theta,-1,musclvarsm1);
          find_musclvarscycle_offset(np,gl,l,theta,+0,musclvarsp0);
          find_musclvarscycle_offset(np,gl,l,theta,+1,musclvarsp1);
          if (hbw_resconv_fluid>=2) find_musclvarscycle_offset(np,gl,l,theta,+2,musclvarsp2);
          if (hbw_resconv_fluid>=3) find_musclvarscycle_offset(np,gl,l,theta,+3,musclvarsp3);
          if (hbw_resconv_fluid>=4) find_musclvarscycle_offset(np,gl,l,theta,+4,musclvarsp4);
          if (hbw_resconv_fluid>=5) find_musclvarscycle_offset(np,gl,l,theta,+5,musclvarsp5);
          find_Fstar_interface_trapezoidal(np, gl, _al(gl,l,theta,-1), theta, musclvarsm5, musclvarsm4, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, musclvarsp4, Fm1h);
          if (gl->cycle.resconv.POSFILTER==POSFILTER_PARENT || gl->cycle.resconv.POSFILTER==POSFILTER_TEST){
            // apply positivity-preserving filter
            find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), l, theta, &metrics);
            for (flux=0; flux<nf; flux++){
              Ftmp[flux]=Fm1h[flux];
            }
            switch (gl->cycle.resconv.POSFILTER) {
              case POSFILTER_PARENT:
                filter_Fstar_interface_positivity_preserving_PARENT(np, gl, _al(gl,l,theta,-1), theta, metrics, gl->cycle.resconv.POSFILTER_numiter, gl->cycle.resconv.EIGENVALCOND, Ftmp, Fm1h, Lambdaminus_p1, Lambdaplus_p0);
              break;
              case POSFILTER_TEST:
                filter_Fstar_interface_positivity_preserving_TEST(np, gl, _al(gl,l,theta,-1), theta, metrics, gl->cycle.resconv.POSFILTER_numiter, gl->cycle.resconv.EIGENVALCOND, Ftmp, Fm1h, Lambdaminus_p1, Lambdaplus_p0);
              break;
              default:
                fatal_error("POSFILTER must be set to POSFILTER_PARENT or POSFILTER_TEST in add_dFstar_residual().");
            }
          }

        } else {

          for (flux=0; flux<nf; flux++) {
            Fm1h[flux]=Fp1h[flux];
          }
          for (flux=0; flux<nmc; flux++) {
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

        for (flux=0; flux<nf; flux++) {
          Lambdaminus_p0[flux][flux]=Lambdaminus_p1[flux][flux];
        }
        switch (hbw_resconv_fluid){
          case 1:
            find_musclvarscycle_offset(np,gl,l,theta,+1,musclvarsp1);
          break;
          case 2:
            find_musclvarscycle_offset(np,gl,l,theta,+2,musclvarsp2);
          break;
          case 3:
            find_musclvarscycle_offset(np,gl,l,theta,+3,musclvarsp3);
          break;
          case 4:
            find_musclvarscycle_offset(np,gl,l,theta,+4,musclvarsp4);
          break;
          case 5:
            find_musclvarscycle_offset(np,gl,l,theta,+5,musclvarsp5);
          break;
          default:
            fatal_error("hbw_resconv_fluid can not be set to %ld",hbw_resconv_fluid);
        }
        find_Fstar_interface_trapezoidal(np, gl, l, theta, musclvarsm4, musclvarsm3, musclvarsm2, musclvarsm1, musclvarsp0, musclvarsp1, musclvarsp2, musclvarsp3, musclvarsp4, musclvarsp5, Fp1h);
        if (gl->cycle.resconv.POSFILTER==POSFILTER_PARENT || gl->cycle.resconv.POSFILTER==POSFILTER_TEST){
          for (flux=0; flux<nf; flux++) {
            Lambdaminus_p0[flux][flux]=Lambdaminus_p1[flux][flux];
          }
          // apply positivity-preserving filter
          find_metrics_at_interface(np, gl, l, _al(gl,l,theta,+1), theta, &metrics);
          for (flux=0; flux<nf; flux++){
            Ftmp[flux]=Fp1h[flux];
          }
          switch (gl->cycle.resconv.POSFILTER){
            case POSFILTER_PARENT:
              filter_Fstar_interface_positivity_preserving_PARENT(np, gl, l, theta, metrics, gl->cycle.resconv.POSFILTER_numiter, gl->cycle.resconv.EIGENVALCOND, Ftmp, Fp1h, Lambdaminus_p1, Lambdaplus_p0);
            break;
            case POSFILTER_TEST:
              filter_Fstar_interface_positivity_preserving_TEST(np, gl, l, theta, metrics, gl->cycle.resconv.POSFILTER_numiter, gl->cycle.resconv.EIGENVALCOND, Ftmp, Fp1h, Lambdaminus_p1, Lambdaplus_p0);
            break;
            default:
              fatal_error("POSFILTER must be set to POSFILTER_PARENT or POSFILTER_TEST in add_dFstar_residual().");
          }
        }

#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
    for (flux=0; flux<nf; flux++) {
      np[l].bs->trapezoidalm1_next[flux]+=gl->cycle.restime.weightm1_trapezoidal_convection*(Fp1h[flux]-Fm1h[flux]);
      np[l].wk->Res[flux]+=(1.0-gl->cycle.restime.weightm1_trapezoidal_convection)*(Fp1h[flux]-Fm1h[flux]);
    }
#else
    for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=(Fp1h[flux]-Fm1h[flux]);
#endif

    for (flux=0; flux<nf; flux++) {
      np[l].bs->Delta_Lambda[theta][flux]=-Lambdaminus_p0[flux][flux]    +Lambdaplus_p0[flux][flux];
    }
  }

}







void find_Delta_Lambda_for_dtau(np_t *np, gl_t *gl, long l, long theta, flux_t Delta_Lambda){
  sqmat_t Lambda_minus,Lambda_plus;
  long flux;
  switch (gl->cycle.resconv.POSFILTER){
    case POSFILTER_NONE:
      find_Lambda_minus_dtau_FDS(np,gl,l,theta,gl->cycle.resconv.EIGENVALCOND,gl->cycle.resconv.AVERAGING,Lambda_minus);
      find_Lambda_plus_dtau_FDS(np,gl,l,theta,gl->cycle.resconv.EIGENVALCOND,gl->cycle.resconv.AVERAGING,Lambda_plus);
      for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=Lambda_plus[flux][flux]-Lambda_minus[flux][flux];
    break;
    case POSFILTER_PARENT:
      for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=np[l].bs->Delta_Lambda[theta][flux];
    break;
    case POSFILTER_TEST:
      for (flux=0; flux<nf; flux++) Delta_Lambda[flux]=np[l].bs->Delta_Lambda[theta][flux];
    break;
    default:
      fatal_error("POSFILTER must be set to either POSFILTER_NONE or POSFILTER_PARENT or POSFILTER_TEST, not to %d",gl->cycle.resconv.POSFILTER);
  }
}



