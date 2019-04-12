#include <cycle/resconv/_resconv.h>
#include <cycle/restime/_restime.h>
#include <cycle/share/cycle_share.h>
#include <cycle/share/res_share.h>
#include <model/_model.h>
#include <src/control.h>


#if (!_FLUID_CONVECTION)
  #error The fluid module specifies no convection terms: choose "none" for the convection terms discretization
#endif

/* in 2D turn off the multidimensional terms if they are a function of CNTBDRYMAX or more boundary nodes */
#ifdef _2D
#define CNTBDRYMAX 2
#endif

/* in 3D turn off the multidimensional terms if they are a function of CNTBDRYMAX or more boundary nodes */
#ifdef _3D
#define CNTBDRYMAX 2
#endif

#define CROSS_DIFFUSION_FACT 1.0
#define FLUXDISCTYPE_ORIGINAL 0
#define FLUXDISCTYPE_EFFICIENT 1
#define FLUXDISCTYPE FLUXDISCTYPE_EFFICIENT


#define AVERAGING_ARITH 0
#define AVERAGING_ROE 1
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



static double _f_cross(double Um1hp1, double Up1hp1, double Um1hp0, double Up1hp0, double Um1hm1, double Up1hm1,
                         np_t *np, gl_t *gl, long lm1hp1, long lp1hp1, long lm1hp0, 
                                   long lp1hp0, long lm1hm1, long lp1hm1){
  double ret,ret1,fact;
  long cntbdry;


  cntbdry=0;
  if (is_node_bdry_no_cross(np[lm1hp1],TYPELEVEL_FLUID_WORK)) cntbdry++;
  if (is_node_bdry_no_cross(np[lp1hp1],TYPELEVEL_FLUID_WORK)) cntbdry++;
  if (is_node_bdry_no_cross(np[lm1hp0],TYPELEVEL_FLUID_WORK)) cntbdry++;
  if (is_node_bdry_no_cross(np[lp1hp0],TYPELEVEL_FLUID_WORK)) cntbdry++;
  if (is_node_bdry_no_cross(np[lm1hm1],TYPELEVEL_FLUID_WORK)) cntbdry++;
  if (is_node_bdry_no_cross(np[lp1hm1],TYPELEVEL_FLUID_WORK)) cntbdry++;

  ret1=0.25*(Um1hp1+Up1hp1-Um1hm1-Up1hm1);
  fact=(minmod(minmod(Um1hp1-Um1hp0,Um1hp0-Um1hm1),minmod(Up1hp1-Up1hp0,Up1hp0-Up1hm1)))/(1.0e-20+ret1);
  if (gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER){
    ret=ret1*(1.0-fact);
  } else {
    ret=ret1;
  }
  if (cntbdry>=CNTBDRYMAX) ret=0.0;
  return(ret);
}



static void find_jacvars_at_interface_normal(np_t *np, gl_t *gl, metrics_t metrics, long lL, long lR, long theta,  jacvars_t *jacvars){
  find_jacvars_at_interface(np, gl, lL, lR, theta, gl->cycle.resconv.AVERAGING, jacvars);
}


static void find_jacvars_at_interface_cross(np_t *np, gl_t *gl, metrics_t metrics, long lL, long lR, long theta, jacvars_t *jacvars){
  find_jacvars_at_interface(np, gl, lL, lR, theta, gl->cycle.resconv.AVERAGING, jacvars);
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



static void find_Aabs_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t Aabs){
  sqmat_t mattmp,R,L,lambdap;

  find_conditioned_Lambda_absolute_from_jacvars_denom(jacvars, metrics, lambdap);
  find_Linv_from_jacvars(jacvars, metrics, R);
  find_L_from_jacvars(jacvars, metrics, L);
  multiply_diagonal_matrix_and_matrix(lambdap,L,mattmp);
  multiply_matrix_and_matrix(R,mattmp,Aabs);
}


static void find_g(flux_t alpham1,flux_t alphap0,flux_t alphap1,flux_t g){
  long flux;
  for (flux=0; flux<nf; flux++){
    g[flux]=minmod3(alpham1[flux],alphap0[flux],alphap1[flux]);
  }
}


static void find_Fstar_interface_2o_Yee(np_t *np, gl_t *gl, metrics_t metrics, long lm1h, long lp1h, 
                     long theta, jacvars_t jacvarsm1, jacvars_t jacvarsp0, jacvars_t jacvarsp1, flux_t Fint){
  flux_t alpham1,alphap0,alphap1,g,fluxtmp;
  sqmat_t R,lambdap;
  long flux,lm3h,lp3h;

  lp3h=_al(gl,lp1h,theta,+1);
  lm3h=_al(gl,lm1h,theta,-1);

  find_alpha(jacvarsm1, metrics, np[lm3h], np[lm1h], gl, alpham1);
  find_alpha(jacvarsp0, metrics, np[lm1h], np[lp1h], gl, alphap0);
  find_alpha(jacvarsp1, metrics, np[lp1h], np[lp3h], gl, alphap1);
  find_g(alpham1,alphap0,alphap1,g);
  find_conditioned_Lambda_absolute_from_jacvars(jacvarsp0, metrics, gl->cycle.resconv.EIGENVALCOND, lambdap);
  find_Linv_from_jacvars(jacvarsp0, metrics, R);
  for (flux=0; flux<nf; flux++)
    fluxtmp[flux]=
                  +lambdap[flux][flux]*g[flux];
  multiply_matrix_and_vector(R,fluxtmp,Fint);



  for (flux=0; flux<nf; flux++) Fint[flux]=0.5e0*(Fint[flux]);
}






static double _Lc(np_t np, long theta){
  double Lc,Xmag;
  long dim;
  Xmag=0.0;
  for (dim=0; dim<nd; dim++) Xmag+=sqr(_X(np,theta,dim));
  Xmag=sqrt(Xmag);
  Lc=1.0/Xmag;
  return(Lc);
}


static double _astar(np_t np, gl_t *gl, long theta){
  double astar;
  astar=_a(np,gl)/_Lc(np,theta);
  return(astar);
}


static void find_Theta(np_t *np, gl_t *gl, long lm1h, long lp1h, long theta1, long theta2, double *Theta){
  double Vstar1,Vstarabs1,sigma1;
  double Vstarabs2,sigma2;  

  Vstar1=minmod(_Vstar(np[lm1h], theta1),_Vstar(np[lp1h], theta1));

  Vstarabs1=0.5*fabs(_Vstar(np[lm1h], theta1)+_Vstar(np[lp1h], theta1));
  sigma1=Vstarabs1+0.5*(_astar(np[lm1h],gl,theta1)+_astar(np[lp1h],gl,theta1));
  Vstarabs1+=gl->model.fluid.zetaA3*sigma1;  

  Vstarabs2=0.5*fabs(_Vstar(np[lm1h], theta2)+_Vstar(np[lp1h], theta2));
  sigma2=Vstarabs2+0.5*(_astar(np[lm1h],gl,theta2)+_astar(np[lp1h],gl,theta2));
  Vstarabs2+=gl->model.fluid.zetaA3*sigma2;  

  *Theta=Vstar1/(Vstarabs1+Vstarabs2);
}



static void find_Fstar_interface_efficient(np_t *np, gl_t *gl, metrics_t metrics, long lm1h, long lp1h,
                     long theta, jacvars_t jacvarsm1, jacvars_t jacvarsp0, jacvars_t jacvarsp1, flux_t Fint){
  long flux;
#ifdef _RESTIME_CDF
  flux_t dFstarxt;
#endif
  flux_t F_normal_2o,Fm1h,Fp1h,alphap0,fluxtmp,fluxtmp2,Um1hp1,Um1hm1,Up1hp1,Up1hm1,Um1hp0,Up1hp0,fluxadd1;
  sqmat_t dUstardUprime,Ax,Ay,R,L,lambdap;
  double Theta;
  metrics_t metricsx,metricsy;
  jacvars_t jacvarsx,jacvarsy,jacvarsm1h,jacvarsp1h;
#ifdef _3D
  flux_t fluxadd2;
  sqmat_t Az;
  jacvars_t jacvarsz;
  metrics_t metricsz;
#endif



  find_alpha(jacvarsp0, metrics, np[lm1h], np[lp1h], gl, alphap0);
  find_conditioned_Lambda_absolute_from_jacvars(jacvarsp0, metrics, gl->cycle.resconv.EIGENVALCOND, lambdap);
  find_Linv_from_jacvars(jacvarsp0, metrics, R);
  find_L_from_jacvars(jacvarsp0, metrics, L);
  for (flux=0; flux<nf; flux++)
    fluxtmp[flux]=-lambdap[flux][flux]*alphap0[flux];
  multiply_matrix_and_vector(R,fluxtmp,Fint);

  /* check that the same solution is obtained whether the following is set to TRUE or FALSE */
  if (TRUE){
    find_Fstar_given_metrics(np[lm1h], gl, metrics, theta, Fm1h);
    find_Fstar_given_metrics(np[lp1h], gl, metrics, theta, Fp1h);
  } else {
    find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
    find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
    find_Fstar_from_jacvars(jacvarsm1h, metrics, Fm1h);
    find_Fstar_from_jacvars(jacvarsp1h, metrics, Fp1h);
  }

  for (flux=0; flux<nf; flux++) Fint[flux]=0.5e0*(Fint[flux]+Fm1h[flux]+Fp1h[flux]);

  /* add second-order yee terms */
  if (!is_node_bdry(np[lm1h],TYPELEVEL_FLUID_WORK) && !is_node_bdry(np[lp1h],TYPELEVEL_FLUID_WORK) && gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER){

    find_Fstar_interface_2o_Yee(np, gl, metrics, lm1h, lp1h, 
                     theta, jacvarsm1, jacvarsp0, jacvarsp1, F_normal_2o);
    for (flux=0; flux<nf; flux++) Fint[flux]+=F_normal_2o[flux];
  }


#ifdef _RESTIME_CDF
  /* now add the spacetime cross-diffusion terms */
    find_Fstarxt_interface(np, gl, theta, lm1h, jacvarsp0, metrics, dFstarxt);
    for (flux=0; flux<nf; flux++) Fint[flux]+=dFstarxt[flux];
#endif

  /* now add the spatial cross-diffusion terms */
#ifdef _2D

    if (theta==0) {
      metricsx=metrics;
      jacvarsx=jacvarsp0; 
      find_Theta(np,gl,lm1h,lp1h,1,0,&Theta);
      find_A_from_jacvars(jacvarsx, metricsx, Ax);
      find_dUstar_dUprime_from_jacvars(jacvarsx, metricsx, dUstardUprime);
      find_Uprime(np[_al(gl,lm1h,1,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,1,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lm1h,1,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lp1h,1,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,1,-1)], gl, Up1hm1);
      find_Uprime(np[_al(gl,lp1h,1,+0)], gl, Up1hp0);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=
          _f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],
             np,gl,_al(gl,lm1h,1,+1), _al(gl,lp1h,1,+1), _al(gl,lm1h,1,+0),
                _al(gl,lp1h,1,+0), _al(gl,lm1h,1,-1), _al(gl,lp1h,1,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Ax,fluxtmp2,fluxadd1); 
    }

    if (theta==1) {
      metricsy=metrics; 
      jacvarsy=jacvarsp0;
      find_Theta(np,gl,lm1h,lp1h,0,1,&Theta);
      find_A_from_jacvars(jacvarsy, metricsy, Ay);
      find_dUstar_dUprime_from_jacvars(jacvarsy, metricsy, dUstardUprime);
      find_Uprime(np[_al(gl,lm1h,0,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,0,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lm1h,0,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lp1h,0,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,0,-1)], gl, Up1hm1);
      find_Uprime(np[_al(gl,lp1h,0,+0)], gl, Up1hp0);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],
             np,gl,_al(gl,lm1h,0,+1),_al(gl,lp1h,0,+1),_al(gl,lm1h,0,+0),
                _al(gl,lp1h,0,+0),_al(gl,lm1h,0,-1),_al(gl,lp1h,0,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Ay,fluxtmp2,fluxadd1);
    }


    for (flux=0; flux<nf; flux++) {
      Fint[flux]-=CROSS_DIFFUSION_FACT*fluxadd1[flux];
      if (isnan(Fint[flux])){
        fatal_error("Problem computing Fint[%ld] in resconv.c; fluxadd1[%ld]=%E. "
                   "Make sure zetaA3 is not zero.",flux,flux,fluxadd1[flux]);
      }
    }
#endif

#ifdef _3D

    if (theta==0) {
      metricsx=metrics; 
      jacvarsx=jacvarsp0;
      find_A_from_jacvars(jacvarsx, metricsx, Ax);
      find_dUstar_dUprime_from_jacvars(jacvarsx, metricsx, dUstardUprime);

      find_Theta(np,gl,lm1h,lp1h,1,0,&Theta);
      find_Uprime(np[_al(gl,lm1h,1,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,1,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,1,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,1,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,1,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,1,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,1,+1), _al(gl,lp1h,1,+1), _al(gl,lm1h,1,+0),
                _al(gl,lp1h,1,+0), _al(gl,lm1h,1,-1), _al(gl,lp1h,1,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Ax,fluxtmp2,fluxadd1);

      find_Theta(np,gl,lm1h,lp1h,2,0,&Theta);
      find_Uprime(np[_al(gl,lm1h,2,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,2,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,2,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,2,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,2,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,2,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,2,+1), _al(gl,lp1h,2,+1), _al(gl,lm1h,2,+0),
                _al(gl,lp1h,2,+0), _al(gl,lm1h,2,-1), _al(gl,lp1h,2,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Ax,fluxtmp2,fluxadd2);
    }


    if (theta==1) {
      metricsy=metrics; 
      jacvarsy=jacvarsp0;
      find_A_from_jacvars(jacvarsy, metricsy, Ay);
      find_dUstar_dUprime_from_jacvars(jacvarsy, metricsy, dUstardUprime);

      find_Theta(np,gl,lm1h,lp1h,0,1,&Theta);
      find_Uprime(np[_al(gl,lm1h,0,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,0,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,0,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,0,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,0,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,0,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,0,+1), _al(gl,lp1h,0,+1), _al(gl,lm1h,0,+0),
                _al(gl,lp1h,0,+0), _al(gl,lm1h,0,-1), _al(gl,lp1h,0,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Ay,fluxtmp2,fluxadd1);

      find_Theta(np,gl,lm1h,lp1h,2,1,&Theta);
      find_Uprime(np[_al(gl,lm1h,2,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,2,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,2,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,2,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,2,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,2,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,2,+1), _al(gl,lp1h,2,+1), _al(gl,lm1h,2,+0),
                _al(gl,lp1h,2,+0), _al(gl,lm1h,2,-1), _al(gl,lp1h,2,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Ay,fluxtmp2,fluxadd2);
    }

    if (theta==2) {
      metricsz=metrics; 
      jacvarsz=jacvarsp0;
      find_A_from_jacvars(jacvarsz, metricsz, Az);
      find_dUstar_dUprime_from_jacvars(jacvarsz, metricsz, dUstardUprime);

      find_Theta(np,gl,lm1h,lp1h,0,2,&Theta);
      find_Uprime(np[_al(gl,lm1h,0,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,0,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,0,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,0,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,0,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,0,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,0,+1), _al(gl,lp1h,0,+1), _al(gl,lm1h,0,+0),
                _al(gl,lp1h,0,+0), _al(gl,lm1h,0,-1), _al(gl,lp1h,0,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Az,fluxtmp2,fluxadd1);

      find_Theta(np,gl,lm1h,lp1h,1,2,&Theta);
      find_Uprime(np[_al(gl,lm1h,1,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,1,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,1,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,1,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,1,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,1,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,1,+1), _al(gl,lp1h,1,+1), _al(gl,lm1h,1,+0),
                _al(gl,lp1h,1,+0), _al(gl,lm1h,1,-1), _al(gl,lp1h,1,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=fluxtmp[flux]*Theta;
      multiply_matrix_and_vector(Az,fluxtmp2,fluxadd2);
    }

    for (flux=0; flux<nf; flux++) {
      Fint[flux]-=CROSS_DIFFUSION_FACT*(fluxadd1[flux]+fluxadd2[flux]);
      if (isnan(Fint[flux])){
        fatal_error("Problem computing Fint[%ld] in resconv.c; fluxadd1[%ld]=%E, fluxadd2[%ld]=%E. "
                   "Make sure zetaA3 is not zero.",flux,flux,fluxadd1[flux],fluxadd2[flux]);
      }
    }
#endif


}



static void find_Fstar_interface_original(np_t *np, gl_t *gl, metrics_t metrics, long lm1h, long lp1h,
                     long theta, jacvars_t jacvarsm1, jacvars_t jacvarsp0, jacvars_t jacvarsp1, flux_t Fint){
  long flux;
#ifdef _RESTIME_CDF
  flux_t dFstarxt;
#endif
  flux_t F_normal_2o,Fm1h,Fp1h,alphap0,fluxtmp,fluxtmp2,Um1hp1,Um1hm1,Up1hp1,Up1hm1,Um1hp0,Up1hp0,fluxadd1;
  sqmat_t dUstardUprime,Ax,Ay,Aabsx,Aabsy,R,L,lambdap,Aden,Adeninv;
  jacvars_t jacvarsx,jacvarsy;
  metrics_t metricsx,metricsy;
  jacvars_t jacvarsm1h,jacvarsp1h;
#ifdef _3D
  flux_t fluxadd2;
  sqmat_t Az,Aabsz;
  jacvars_t jacvarsz;
  metrics_t metricsz;
#endif



  find_alpha(jacvarsp0, metrics, np[lm1h], np[lp1h], gl, alphap0);
  find_conditioned_Lambda_absolute_from_jacvars(jacvarsp0, metrics, gl->cycle.resconv.EIGENVALCOND, lambdap);
  find_Linv_from_jacvars(jacvarsp0, metrics, R);
  find_L_from_jacvars(jacvarsp0, metrics, L);
  for (flux=0; flux<nf; flux++)
    fluxtmp[flux]=-lambdap[flux][flux]*alphap0[flux];
  multiply_matrix_and_vector(R,fluxtmp,Fint);

  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_Fstar_from_jacvars(jacvarsm1h, metrics, Fm1h);
  find_Fstar_from_jacvars(jacvarsp1h, metrics, Fp1h);

  for (flux=0; flux<nf; flux++) Fint[flux]=0.5e0*(Fint[flux]+Fm1h[flux]+Fp1h[flux]);

  /* add second-order yee terms */
  if (!is_node_bdry(np[lm1h],TYPELEVEL_FLUID_WORK) && !is_node_bdry(np[lp1h],TYPELEVEL_FLUID_WORK) &&  gl->cycle.resconv.ACCURACY==ACCURACY_SECONDORDER){

    find_Fstar_interface_2o_Yee(np, gl, metrics, lm1h, lp1h, 
                     theta, jacvarsm1, jacvarsp0, jacvarsp1, F_normal_2o);
    for (flux=0; flux<nf; flux++) Fint[flux]+=F_normal_2o[flux];
  }


#ifdef _RESTIME_CDF
  /* now add the spacetime cross-diffusion terms */
    find_Fstarxt_interface(np, gl, theta, lm1h, jacvarsp0, metrics, dFstarxt);
    for (flux=0; flux<nf; flux++) Fint[flux]+=dFstarxt[flux];
#endif

  /* now add the spatial cross-diffusion terms */
#ifdef _2D
    find_metrics_at_interface_2(np, gl, lm1h, lp1h, theta, 0, &metricsx);
    find_metrics_at_interface_2(np, gl, lm1h, lp1h, theta, 1, &metricsy);

    find_jacvars_at_interface_cross(np,gl,metricsx, lm1h,lp1h,0,&jacvarsx);
    find_jacvars_at_interface_cross(np,gl,metricsy, lm1h,lp1h,1,&jacvarsy);



    find_A_from_jacvars(jacvarsx, metricsx, Ax);
    find_A_from_jacvars(jacvarsy, metricsy, Ay);
    find_Aabs_from_jacvars(jacvarsx, metricsx, Aabsx);
    find_Aabs_from_jacvars(jacvarsy, metricsy, Aabsy);

    add_two_matrices(Aabsx,Aabsy,Aden);
    invert_matrix(Aden,Adeninv);

    if (theta==0) {
      find_dUstar_dUprime_from_jacvars(jacvarsx, metricsx, dUstardUprime);
      find_Uprime(np[_al(gl,lm1h,1,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,1,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lm1h,1,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lp1h,1,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,1,-1)], gl, Up1hm1);
      find_Uprime(np[_al(gl,lp1h,1,+0)], gl, Up1hp0);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=
          _f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],
             np,gl,_al(gl,lm1h,1,+1), _al(gl,lp1h,1,+1), _al(gl,lm1h,1,+0),
                _al(gl,lp1h,1,+0), _al(gl,lm1h,1,-1), _al(gl,lp1h,1,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ay,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ax,fluxtmp,fluxadd1);
 
    }

    if (theta==1) {
      find_dUstar_dUprime_from_jacvars(jacvarsy, metricsy, dUstardUprime);
      find_Uprime(np[_al(gl,lm1h,0,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,0,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lm1h,0,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lp1h,0,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,0,-1)], gl, Up1hm1);
      find_Uprime(np[_al(gl,lp1h,0,+0)], gl, Up1hp0);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],
             np,gl,_al(gl,lm1h,0,+1),_al(gl,lp1h,0,+1),_al(gl,lm1h,0,+0),
                _al(gl,lp1h,0,+0),_al(gl,lm1h,0,-1),_al(gl,lp1h,0,-1));
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ax,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ay,fluxtmp,fluxadd1);
    }


    for (flux=0; flux<nf; flux++) {
      Fint[flux]-=CROSS_DIFFUSION_FACT*fluxadd1[flux];
      if (isnan(Fint[flux])){
        fatal_error("Problem computing Fint[%ld] in resconv.c; fluxadd1[%ld]=%E. "
                   "Make sure zetaA3 is not zero.",flux,flux,fluxadd1[flux]);
      }

    }
#endif

#ifdef _3D

    find_metrics_at_interface_2(np, gl, lm1h, lp1h, theta, 0, &metricsx);
    find_metrics_at_interface_2(np, gl, lm1h, lp1h, theta, 1, &metricsy);
    find_metrics_at_interface_2(np, gl, lm1h, lp1h, theta, 2, &metricsz);

    find_jacvars_at_interface_cross(np,gl,metricsx,lm1h,lp1h,0,&jacvarsx);
    find_jacvars_at_interface_cross(np,gl,metricsy,lm1h,lp1h,1,&jacvarsy);
    find_jacvars_at_interface_cross(np,gl,metricsz,lm1h,lp1h,2,&jacvarsz);

    find_A_from_jacvars(jacvarsx, metricsx, Ax);
    find_A_from_jacvars(jacvarsy, metricsy, Ay);
    find_A_from_jacvars(jacvarsz, metricsz, Az);
    find_Aabs_from_jacvars(jacvarsx, metricsx, Aabsx);
    find_Aabs_from_jacvars(jacvarsy, metricsy, Aabsy);
    find_Aabs_from_jacvars(jacvarsz, metricsz, Aabsz);


    if (theta==0) {
      find_dUstar_dUprime_from_jacvars(jacvarsx, metricsx, dUstardUprime);
      find_Uprime(np[_al(gl,lm1h,1,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,1,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,1,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,1,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,1,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,1,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,1,+1), _al(gl,lp1h,1,+1), _al(gl,lm1h,1,+0),
                _al(gl,lp1h,1,+0), _al(gl,lm1h,1,-1), _al(gl,lp1h,1,-1));
      add_two_matrices(Aabsx,Aabsy,Aden);
      invert_matrix(Aden,Adeninv);
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ay,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ax,fluxtmp,fluxadd1);





      find_Uprime(np[_al(gl,lm1h,2,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,2,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,2,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,2,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,2,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,2,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,2,+1), _al(gl,lp1h,2,+1), _al(gl,lm1h,2,+0),
                _al(gl,lp1h,2,+0), _al(gl,lm1h,2,-1), _al(gl,lp1h,2,-1));
      add_two_matrices(Aabsx,Aabsz,Aden);
      invert_matrix(Aden,Adeninv);
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Az,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ax,fluxtmp,fluxadd2);
    }

    if (theta==1) {
      find_dUstar_dUprime_from_jacvars(jacvarsy, metricsy, dUstardUprime);
      find_Uprime(np[_al(gl,lm1h,0,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,0,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,0,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,0,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,0,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,0,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,0,+1), _al(gl,lp1h,0,+1), _al(gl,lm1h,0,+0),
                _al(gl,lp1h,0,+0), _al(gl,lm1h,0,-1), _al(gl,lp1h,0,-1));
      add_two_matrices(Aabsx,Aabsy,Aden);
      invert_matrix(Aden,Adeninv);
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ax,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ay,fluxtmp,fluxadd1);

      find_Uprime(np[_al(gl,lm1h,2,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,2,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,2,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,2,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,2,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,2,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,2,+1), _al(gl,lp1h,2,+1), _al(gl,lm1h,2,+0),
                _al(gl,lp1h,2,+0), _al(gl,lm1h,2,-1), _al(gl,lp1h,2,-1));
      add_two_matrices(Aabsz,Aabsy,Aden);
      invert_matrix(Aden,Adeninv);
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Az,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ay,fluxtmp,fluxadd2);
    }

    if (theta==2) {
      find_dUstar_dUprime_from_jacvars(jacvarsz, metricsz, dUstardUprime);
      find_Uprime(np[_al(gl,lm1h,0,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,0,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,0,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,0,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,0,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,0,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,0,+1), _al(gl,lp1h,0,+1), _al(gl,lm1h,0,+0),
                _al(gl,lp1h,0,+0), _al(gl,lm1h,0,-1), _al(gl,lp1h,0,-1));
      add_two_matrices(Aabsx,Aabsz,Aden);
      invert_matrix(Aden,Adeninv);
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ax,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Az,fluxtmp,fluxadd1);

      find_Uprime(np[_al(gl,lm1h,1,+1)], gl, Um1hp1);
      find_Uprime(np[_al(gl,lm1h,1,+0)], gl, Um1hp0);
      find_Uprime(np[_al(gl,lm1h,1,-1)], gl, Um1hm1);
      find_Uprime(np[_al(gl,lp1h,1,+1)], gl, Up1hp1);
      find_Uprime(np[_al(gl,lp1h,1,+0)], gl, Up1hp0);
      find_Uprime(np[_al(gl,lp1h,1,-1)], gl, Up1hm1);
      for (flux=0; flux<nf; flux++) fluxtmp2[flux]=_f_cross(Um1hp1[flux],Up1hp1[flux],Um1hp0[flux],Up1hp0[flux],Um1hm1[flux],Up1hm1[flux],np,gl,_al(gl,lm1h,1,+1), _al(gl,lp1h,1,+1), _al(gl,lm1h,1,+0),
                _al(gl,lp1h,1,+0), _al(gl,lm1h,1,-1), _al(gl,lp1h,1,-1));
      add_two_matrices(Aabsz,Aabsy,Aden);
      invert_matrix(Aden,Adeninv);
      multiply_matrix_and_vector(dUstardUprime,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Ay,fluxtmp,fluxtmp2);
      multiply_matrix_and_vector(Adeninv,fluxtmp2,fluxtmp);
      multiply_matrix_and_vector(Az,fluxtmp,fluxadd2);
    }

    for (flux=0; flux<nf; flux++) {
      Fint[flux]-=CROSS_DIFFUSION_FACT*(fluxadd1[flux]+fluxadd2[flux]);
      if (isnan(Fint[flux])){
        fatal_error("Problem computing Fint[%ld] in resconv.c; fluxadd1[%ld]=%E, fluxadd2[%ld]=%E. "
                   "Make sure zetaA3 is not zero.",flux,flux,fluxadd1[flux],fluxadd2[flux]);
      }
    }
#endif


}

static void find_Fstar_interface(np_t *np, gl_t *gl, metrics_t metrics, long lm1h, long lp1h,
                     long theta, jacvars_t jacvarsm1, jacvars_t jacvarsp0, jacvars_t jacvarsp1, flux_t Fint){
  switch (FLUXDISCTYPE){
   case FLUXDISCTYPE_ORIGINAL:
    find_Fstar_interface_original(np, gl, metrics, lm1h, lp1h, theta, jacvarsm1, jacvarsp0, jacvarsp1, Fint);
   break;
   case FLUXDISCTYPE_EFFICIENT: 
    find_Fstar_interface_efficient(np, gl, metrics, lm1h, lp1h, theta, jacvarsm1, jacvarsp0, jacvarsp1, Fint);
   break;
  }

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




void add_dFstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  jacvars_t *jacvarstmp,*jacvarsm3h,*jacvarsm1h,*jacvarsp1h,*jacvarsp3h;
  flux_t Fm1h,Fp1h;
  long flux,l;
  metrics_t metricsp1h,metricsm1h;

  jacvarsm3h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsm1h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsp1h= (jacvars_t *) malloc(sizeof(jacvars_t));
  jacvarsp3h= (jacvars_t *) malloc(sizeof(jacvars_t));

  l=_l_minus_one(ls,gl,theta);
  find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metricsp1h);
  if (is_node_bdry(np[l],TYPELEVEL_FLUID_WORK))
    find_jacvars_at_interface_normal(np,gl,metricsp1h,l,l,theta,jacvarsm1h);
    else find_jacvars_at_interface_normal(np,gl,metricsp1h,_al(gl,l,theta,-1),_al(gl,l,theta,+0),theta,jacvarsm1h);
  find_jacvars_at_interface_normal(np,gl,metricsp1h,_al(gl,l,theta,+0),_al(gl,l,theta,+1),theta,jacvarsp1h);
  find_jacvars_at_interface_normal(np,gl,metricsp1h,_al(gl,l,theta,+1),_al(gl,l,theta,+2),theta,jacvarsp3h);

  find_Fstar_interface(np, gl,metricsp1h,_al(gl,l,theta,+0),_al(gl,l,theta,+1),
               theta, *jacvarsm1h,*jacvarsp1h,*jacvarsp3h,Fp1h);


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
    find_jacvars_at_interface_normal(np,gl,metricsp1h,_al(gl,l,theta,+1),_al(gl,l,theta,+2),theta,jacvarsp3h);
    else find_jacvars_at_interface_normal(np,gl,metricsp1h,_al(gl,l,theta,+1),_al(gl,l,theta,+1),theta,jacvarsp3h);

    find_Fstar_interface(np, gl, metricsp1h,_al(gl,l,theta,+0),_al(gl,l,theta,+1),
               theta, *jacvarsm1h,*jacvarsp1h,*jacvarsp3h,Fp1h);

#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
    for (flux=0; flux<nf; flux++) {
      np[l].bs->Res_trapezoidal[flux]+=gl->cycle.restime.weightm1_trapezoidal_convection*(Fp1h[flux]-Fm1h[flux]);
      np[l].wk->Res[flux]+=(1.0-gl->cycle.restime.weightm1_trapezoidal_convection)*(Fp1h[flux]-Fm1h[flux]);
    }
#else
    for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=(Fp1h[flux]-Fm1h[flux]);
#endif


#ifdef _RESCONV_DELTA_LAMBDA_STORAGE
    find_Delta_Lambda_for_dtau_local(np, gl, l, theta, *jacvarsm1h, *jacvarsp1h, metricsm1h, metricsp1h,  np[l].bs->Delta_Lambda[theta]); 
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
