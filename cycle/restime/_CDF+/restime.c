#include <cycle/restime/_restime.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <cycle/share/res_share.h>
#include <src/control.h>



void write_disc_restime_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    xi1=0.5; {except for momentum and energy fluxes}\n"
    "    xi2=0.5; {momentum fluxes}\n"
    "    xi3=0.5; {total energy flux}\n"
    "  );\n"
  ,_RESTIME_ACTIONNAME);

}

void read_disc_restime_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_disc_restime_actions(char *actionname, char **argum, SOAP_codex_t *codex){

  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_RESTIME_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_RESTIME_ACTIONNAME);


    gl->DISC_RESTIME_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_restime_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_double_var_from_codex(codex,"xi1",&gl->cycle.restime.xi1);
    find_double_var_from_codex(codex,"xi2",&gl->cycle.restime.xi2);
    find_double_var_from_codex(codex,"xi3",&gl->cycle.restime.xi3);
    if (gl->cycle.restime.xi1<0.0)
      SOAP_fatal_error(codex,"xi1 must be set to a value equal or greater to 0.");
    if (gl->cycle.restime.xi2<0.0)
      SOAP_fatal_error(codex,"xi2 must be set to a value equal or greater to 0.");
    if (gl->cycle.restime.xi3<0.0)
      SOAP_fatal_error(codex,"xi3 must be set to a value equal or greater to 0.");

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}



void find_Lambdaxt(np_t *np, gl_t *gl, long l, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxt){
  long flux;
  sqmat_t Lambda,LambdaZ;

  find_LambdaZ(np,gl,l,LambdaZ);
  set_matrix_to_zero(Lambdaxt);
  find_Lambda_from_jacvars(jacvars, metrics, Lambda);

  for (flux=0; flux<nf; flux++) {
    assert(LambdaZ[flux][flux]>=0.0);
    Lambdaxt[flux][flux]=Lambda[flux][flux]/max(1.0,gl->dt*fabs(Lambda[flux][flux])/max(1e-99,LambdaZ[flux][flux]));
  }
}


void find_Lambdaxt_interface(np_t *np, gl_t *gl, long lp0, long lp1, long theta, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxt){
  long row,col,flux;
  sqmat_t Lambda,LambdaZp0,LambdaZp1,LambdaZ;
  

  find_LambdaZ(np,gl,lp0,LambdaZp0);
  find_LambdaZ(np,gl,lp1,LambdaZp1);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      LambdaZ[row][col]=0.5*(LambdaZp0[row][col]+LambdaZp1[row][col]);
    }
  }
  find_Lambda_from_jacvars(jacvars, metrics, Lambda);

  set_matrix_to_zero(Lambdaxt);
  for (flux=0; flux<nf; flux++) {
    Lambdaxt[flux][flux]=Lambda[flux][flux]/max(1.0,gl->dt*fabs(Lambda[flux][flux])/max(1e-99,LambdaZ[flux][flux]));
  }
}



static double _psi(gl_t *gl, double val0, double val1, long flux){
  double xi,ret;
  xi=gl->cycle.restime.xi1;
  if (flux>=fluxmom && flux<fluxmom+nd){
    xi=gl->cycle.restime.xi2;
  }
  if (flux==fluxet){
    xi=gl->cycle.restime.xi3;
  }
  ret=minmod(val0,val1)*min(1.0,xi*max(fabs(val1),fabs(val0))/notzero(min(fabs(val0),fabs(val1)),1e-40));
  return(ret);
}



static void find_Deltaxt(np_t *np, gl_t *gl, long l, long theta, jacvars_t jacvars, metrics_t metrics, flux_t Deltaxt){
  flux_t vectmp,dU,LUstarp0;
  sqmat_t L;
  long flux; 


  find_LUstar_from_jacvars(jacvars,metrics,LUstarp0);
  find_L_from_jacvars(jacvars,metrics,L);
  for (flux=0; flux<nf; flux++) dU[flux]=np[l].bs->U[flux]-np[l].bs->Um1[flux];
  multiply_matrix_and_vector(L,dU,vectmp);

  /* find Deltaxt using phi */
  for (flux=0; flux<nf; flux++) Deltaxt[flux]=metrics.Omega*vectmp[flux]/notzero(LUstarp0[flux],1e-99);

}


void find_Deltaxt_interface(np_t *np, gl_t *gl, long lp0, long lp1, long theta, 
                            jacvars_t jacvarsp0, jacvars_t jacvarsp1,  metrics_t metrics, flux_t Deltaxt){
  long flux;
  flux_t Deltaxt_p0,Deltaxt_p1;

  find_Deltaxt(np,gl,lp0,theta,jacvarsp0,metrics,Deltaxt_p0);
  find_Deltaxt(np,gl,lp1,theta,jacvarsp1,metrics,Deltaxt_p1);  

  for (flux=0; flux<nf; flux++){
    Deltaxt[flux]=_psi(gl,Deltaxt_p0[flux],Deltaxt_p1[flux],flux); 
  }
}



/* jacvars and metrics must be evaluated between node lp0 and node lp1 */
void find_Lambdaxt_plus_dtau_from_jacvars(np_t *np, gl_t *gl, long lp0, long lp1, long theta, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxtplus){
  long flux;
  sqmat_t Lambdaxt;
  find_Lambdaxt_interface(np, gl, lp0, lp1, theta, jacvars, metrics, Lambdaxt);
  set_matrix_to_zero(Lambdaxtplus);
  for (flux=0; flux<nf; flux++) Lambdaxtplus[flux][flux]=0.5*(Lambdaxt[flux][flux]+fabs(Lambdaxt[flux][flux]));
}


void find_Lambdaxt_plus_dtau(np_t *np, gl_t *gl, long l, long theta, sqmat_t Lambdaxtplus){
  jacvars_t jacvarsL,jacvarsR,jacvars;
  metrics_t metrics;
  long lp0,lp1;
  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);
  find_metrics_at_interface(np,gl,lp0,lp1,theta,&metrics);
  find_jacvars(np[lp0], gl, metrics, theta, &jacvarsL);
  find_jacvars(np[lp1], gl, metrics, theta, &jacvarsR);
  find_jacvars_at_interface_Roe_average(jacvarsL,jacvarsR,gl,theta,&jacvars);
  find_Lambdaxt_plus_dtau_from_jacvars(np, gl, lp0, lp1, theta, jacvars, metrics, Lambdaxtplus);
}


/* jacvars and metrics must be evaluated between node lp0 and node lp1 */
void find_Lambdaxt_minus_dtau_from_jacvars(np_t *np, gl_t *gl, long lp0, long lp1, long theta, jacvars_t jacvars, metrics_t metrics, sqmat_t Lambdaxtminus){
  long flux;
  sqmat_t Lambdaxt;
  find_Lambdaxt_interface(np, gl, lp0, lp1, theta, jacvars, metrics, Lambdaxt);
  set_matrix_to_zero(Lambdaxtminus);
  for (flux=0; flux<nf; flux++) Lambdaxtminus[flux][flux]=0.5*(Lambdaxt[flux][flux]-fabs(Lambdaxt[flux][flux]));
}


void find_Lambdaxt_minus_dtau(np_t *np, gl_t *gl, long l, long theta, sqmat_t Lambdaxtminus){
  long lp0,lm1;
  metrics_t metrics;
  jacvars_t jacvarsL,jacvarsR,jacvars;
  lp0=_al(gl,l,theta,+0);
  lm1=_al(gl,l,theta,-1);
  find_metrics_at_interface(np,gl,lm1,lp0,theta,&metrics);
  find_jacvars(np[lm1], gl, metrics, theta, &jacvarsL);
  find_jacvars(np[lp0], gl, metrics, theta, &jacvarsR);
  find_jacvars_at_interface_Roe_average(jacvarsL,jacvarsR,gl,theta,&jacvars);

  find_Lambdaxt_minus_dtau_from_jacvars(np, gl, lm1, lp0, theta, jacvars, metrics, Lambdaxtminus);

}


void add_Z_dUstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl, double fact, double fact_trapezoidal){
  long l;
  long flux;
  flux_t tmp,dRes;
  sqmat_t Z;

    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      for (flux=0; flux<nf; flux++){
        tmp[flux]=_Omega(np[l],gl)*(np[l].bs->U[flux]-np[l].bs->Um1[flux])/gl->dt;
      }      
      find_Z(np,gl,l,Z);
      multiply_matrix_and_vector(Z, tmp, dRes);
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=fact*dRes[flux];
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
      for (flux=0; flux<nf; flux++) np[l].bs->Res_trapezoidal[flux]+=fact_trapezoidal*dRes[flux];
#endif
    }
}



/* find space-time flux contribution -> this can only be used with pure FDS schemes without positivity-preservation and without Reconstruction-Evolution */ 
void find_Fstarxt_interface(np_t *np, gl_t *gl, long theta, long l, jacvars_t jacvars, metrics_t metrics, flux_t dFstar){
  sqmat_t Linv,Lambdaxt;
  long flux,lp1,flux2;
  flux_t LUstar,Deltaxt_p0,Deltaxt_p1,tmpflux;

  lp1=_al(gl,l,theta,+1);

  find_Linv_from_jacvars(jacvars, metrics, Linv);
  find_LUstar_from_jacvars(jacvars, metrics, LUstar);
  find_Lambdaxt_interface(np, gl, l, lp1, theta, jacvars, metrics, Lambdaxt);
  find_Deltaxt(np, gl, l, theta, jacvars, metrics, Deltaxt_p0);
  find_Deltaxt(np, gl, lp1, theta, jacvars, metrics, Deltaxt_p1);
  for (flux=0; flux<nf; flux++) tmpflux[flux]=-0.5*(Lambdaxt[flux][flux]-fabs(Lambdaxt[flux][flux]))*Deltaxt_p1[flux]*LUstar[flux]-0.5*(Lambdaxt[flux][flux]+fabs(Lambdaxt[flux][flux]))*Deltaxt_p0[flux]*LUstar[flux];

  multiply_matrix_and_vector(Linv,tmpflux,dFstar);

  for (flux=0; flux<nf; flux++){
    if (isnan(dFstar[flux])){
      for (flux2=0; flux2<nf; flux2++) wfprintf(stderr,"LUstar[%ld]=%E\n",flux2,LUstar[flux2]);
      wfprintf(stderr,"\n");
      for (flux2=0; flux2<nf; flux2++) wfprintf(stderr,"Deltaxt_p0[%ld]=%E\n",flux2,Deltaxt_p0[flux2]);
      wfprintf(stderr,"\n");
      fatal_error("Problem computing dFstar[%ld] in restime.c;   dFstar[%ld]=%E."
                  ,flux,flux,dFstar[flux]);
      
    }
  }
}



void find_dFstarxt_dUstar_interface(np_t *np, gl_t *gl, long theta, long l, jacvars_t jacvarsp1h, metrics_t metricsp1h, sqmat_t B, sqmat_t C){
  sqmat_t tmpmat,Aplus,Aminus,L,Linv,Lambdaxtplus,Lambdaxtminus,Lambdaxt;
  double Omegap0,Omegap1,Omegap1h;
  long lp0,lp1,row,col,flux;

  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);


  find_Linv_from_jacvars(jacvarsp1h, metricsp1h, Linv);
  find_L_from_jacvars(jacvarsp1h, metricsp1h, L);
  find_Lambdaxt_interface(np, gl, l, lp1, theta, jacvarsp1h, metricsp1h, Lambdaxt);
  set_matrix_to_zero(Lambdaxtminus);
  set_matrix_to_zero(Lambdaxtplus);
  for (flux=0; flux<nf; flux++){
    Lambdaxtplus[flux][flux]=0.5*(Lambdaxt[flux][flux]+fabs(Lambdaxt[flux][flux]));
    Lambdaxtminus[flux][flux]=0.5*(Lambdaxt[flux][flux]-fabs(Lambdaxt[flux][flux]));
  }

  multiply_diagonal_matrix_and_matrix(Lambdaxtplus,L,tmpmat);
  multiply_matrix_and_matrix(Linv,tmpmat,Aplus);

  multiply_diagonal_matrix_and_matrix(Lambdaxtminus,L,tmpmat);
  multiply_matrix_and_matrix(Linv,tmpmat,Aminus);

  Omegap1h=metricsp1h.Omega;
  Omegap0=_Omega(np[lp0],gl);
  Omegap1=_Omega(np[lp1],gl);  

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=Omegap1h/Omegap0*(-Aplus[row][col]);
      C[row][col]=Omegap1h/Omegap1*(-Aminus[row][col]);
    }
  }


}


