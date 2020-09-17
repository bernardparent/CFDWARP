// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2015 Bernard Parent

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

#include <cycle/share/ts_share.h>
#include <cycle/share/res_share.h>
#include <cycle/restime/_restime.h>
#include <model/_model.h>



void find_Kstar_dG_dUstar_interface(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C){
  sqmat_t Kp1h;
  sqmat_t Bp0,Bp1;
  long lp0,lp1,row,col;
  metrics_t metricsp1h;

  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);

  find_metrics_at_interface(np, gl, lp0, lp1, theta, &metricsp1h);

  find_dG_dUstar(np[lp0],gl,Bp0);
  find_dG_dUstar(np[lp1],gl,Bp1);

  find_Kstar_interface(np,gl,lp0,lp1,metricsp1h,theta,theta,Kp1h,CYCLELEVEL_TS);
 
  multiply_matrix_and_matrix(Kp1h,Bp0,B);
  multiply_matrix_and_matrix(Kp1h,Bp1,C);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      C[row][col]=-C[row][col];
    }
  }
}


void add_Kstar_dG_dUstar_to_TDMA_check(np_t *np, gl_t *gl, long theta, long l, sqmat_t A, sqmat_t B, sqmat_t C){
  sqmat_t Bp1h,Cp1h,Am1h,Bm1h;
  long row,col;

  find_Kstar_dG_dUstar_interface(np, gl, theta, l, Bp1h, Cp1h);
  find_Kstar_dG_dUstar_interface(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      A[row][col]+=-Am1h[row][col];
      B[row][col]+=+Bp1h[row][col]-Bm1h[row][col];
      C[row][col]+=+Cp1h[row][col];
    }
  }
}



void add_Kstar_dG_dUstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C){
  sqmat_t Km1h,Kp1h;
  sqmat_t Bm1,Bp0,Bp1;
  sqmat_t A1,B1,B2,C1;
  long lm1,lp0,lp1,row,col;
  metrics_t metricsm1h,metricsp1h;

  lm1=_al(gl,l,theta,-1);
  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);

  find_metrics_at_interface(np, gl, lm1, lp0, theta, &metricsm1h);
  find_metrics_at_interface(np, gl, lp0, lp1, theta, &metricsp1h);

  find_dG_dUstar(np[lp1],gl,Bp1);
  find_dG_dUstar(np[lp0],gl,Bp0);
  find_dG_dUstar(np[lm1],gl,Bm1);

  find_Kstar_interface(np,gl,lm1,lp0,metricsm1h,theta,theta,Km1h,CYCLELEVEL_TS);
  find_Kstar_interface(np,gl,lp0,lp1,metricsp1h,theta,theta,Kp1h,CYCLELEVEL_TS);
  //make_matrix_positive(Km1h);
  //make_matrix_positive(Kp1h);
  multiply_matrix_and_matrix(Kp1h,Bp0,B1);
  multiply_matrix_and_matrix(Km1h,Bp0,B2);
  multiply_matrix_and_matrix(Km1h,Bm1,A1);
  multiply_matrix_and_matrix(Kp1h,Bp1,C1);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      A[row][col]+=-fact*A1[row][col];
      B[row][col]+=fact*(B1[row][col]+B2[row][col]);
      C[row][col]+=-fact*C1[row][col];
    }
  }
}




#ifdef _FLUID_PLASMA


void find_Dstar_interface(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C){
  sqmat_t Dstarplusp0, Dstarminusp1;
  long lp0,lp1,row,col;
  

  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);
  find_Dstarplus(np, gl, lp0, theta, Dstarplusp0);
  find_Dstarminus(np, gl, lp1, theta, Dstarminusp1);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=+Dstarplusp0[row][col];
      C[row][col]=+Dstarminusp1[row][col];
    }
  }

}


void add_Dstar_to_TDMA_check(np_t *np, gl_t *gl, long theta, long l, sqmat_t A, sqmat_t B, sqmat_t C){
  sqmat_t Bp1h,Cp1h,Am1h,Bm1h;
  long row,col;

  find_Dstar_interface(np, gl, theta, l, Bp1h, Cp1h);
  find_Dstar_interface(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      A[row][col]+=-Am1h[row][col];
      B[row][col]+=+Bp1h[row][col]-Bm1h[row][col];
      C[row][col]+=+Cp1h[row][col];
    }
  }


}



void add_Dstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C){
  sqmat_t Dstarplusm1, Dstarplusp0, Dstarminusp0, Dstarminusp1;
  long lm1,lp0,lp1,row,col;

  lm1=_al(gl,l,theta,-1);
  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);

  find_Dstarplus(np, gl, lm1, theta, Dstarplusm1);
  find_Dstarplus(np, gl, lp0, theta, Dstarplusp0);
  find_Dstarminus(np, gl, lp0, theta, Dstarminusp0);
  find_Dstarminus(np, gl, lp1, theta, Dstarminusp1);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      A[row][col]+=-fact*Dstarplusm1[row][col];
      B[row][col]+=+fact*(Dstarplusp0[row][col]-Dstarminusp0[row][col]);
      C[row][col]+=+fact*Dstarminusp1[row][col];
    }
  }
}



void add_Ystar_dH_dUstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C){
  flux_t Ystarm1h,Ystarp1h;
  long row,col;
  sqmat_t dHdUstarm1,dHdUstarp0,dHdUstarp1;
  sqmat_t matm1h,matp1h;
  sqmat_t matm1h_times_dHdUstarm1,matm1h_times_dHdUstarp0,matp1h_times_dHdUstarp0,matp1h_times_dHdUstarp1;
  metrics_t metricsp0;
  int cnt;

  find_metrics_at_node(np,gl, l, theta, &metricsp0);

  for (cnt=1; cnt<=2; cnt++){
    switch (cnt){
      case 1:
        find_Y1star_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, Ystarm1h);
        find_Y1star_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, Ystarp1h);
        find_dH1_dUstar(np, gl, _al(gl,l,theta,-1), dHdUstarm1);  
        find_dH1_dUstar(np, gl, _al(gl,l,theta,+0), dHdUstarp0);  
        find_dH1_dUstar(np, gl, _al(gl,l,theta,+1), dHdUstarp1);  
      break;
      case 2:
        find_Y2star_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, Ystarm1h);
        find_Y2star_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, Ystarp1h);
        find_dH2_dUstar(np, gl, _al(gl,l,theta,-1), dHdUstarm1);  
        find_dH2_dUstar(np, gl, _al(gl,l,theta,+0), dHdUstarp0);  
        find_dH2_dUstar(np, gl, _al(gl,l,theta,+1), dHdUstarp1);  
      break;
      default:
        fatal_error("cnt can not be set to %d in add_Ystar_dH_dUstar_to_TDMA().",cnt);
    }
    for (row=0; row<nf; row++){
      for (col=0; col<nf; col++){
        matm1h[row][col]=0.0;
        matp1h[row][col]=0.0;
      }
    }
    for (row=0; row<nf; row++){
      matm1h[row][row]=0.5*(Ystarm1h[row]+fabs(Ystarm1h[row]));
      matp1h[row][row]=0.5*(Ystarp1h[row]-fabs(Ystarp1h[row]));
    }
    multiply_diagonal_matrix_and_matrix(matm1h, dHdUstarm1, matm1h_times_dHdUstarm1);
    multiply_diagonal_matrix_and_matrix(matm1h, dHdUstarp0, matm1h_times_dHdUstarp0);
    multiply_diagonal_matrix_and_matrix(matp1h, dHdUstarp0, matp1h_times_dHdUstarp0);
    multiply_diagonal_matrix_and_matrix(matp1h, dHdUstarp1, matp1h_times_dHdUstarp1);


    for (row=0; row<nf; row++){
      for (col=0; col<nf; col++){
        A[row][col]+=-fact*matm1h_times_dHdUstarm1[row][col];
        B[row][col]+=fact*(matm1h_times_dHdUstarp0[row][col]-matp1h_times_dHdUstarp0[row][col]);
        C[row][col]+=fact*matp1h_times_dHdUstarp1[row][col];
      }
    }

  }

}


#endif


void add_dSstar_dUstar_to_TDMA(np_t *np, gl_t *gl, long l,  double fact, sqmat_t B){
  sqmat_t dSstar_dUstar;
  long row,col;

#ifdef _RESSOURCE_LAMBDAMINUS_STORAGE
  sqmat_t Lambdaminus,Linv,L,mattmp;
  long theta;
  jacvars_t jacvars;
  metrics_t metrics;
  theta=0;
  find_metrics_at_node(np, gl, l, theta, &metrics);
  find_jacvars(np[l], gl, metrics, theta, &jacvars);
  find_Linv_from_jacvars(jacvars, metrics, Linv);
  find_L_from_jacvars(jacvars, metrics, L);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Lambdaminus[row][col]=0.0;
    }
    Lambdaminus[row][row]=np[l].bs->Lambda_S_minus[row];
  }
  multiply_diagonal_matrix_and_matrix(Lambdaminus,L,mattmp);
  multiply_matrix_and_matrix(Linv,mattmp,dSstar_dUstar);
#else
  find_dSstar_dUstar(np,gl,l,dSstar_dUstar);
#endif
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]+=-fact*dSstar_dUstar[row][col];
    }
  }
}


#ifdef UNSTEADY
void add_Z_dUstardt_dUstar_to_TDMA(np_t *np, gl_t *gl, long l,  sqmat_t B){
  long row,col;
  sqmat_t dZU_dU;

  find_dZU_dU(np, gl, l, dZU_dU);
  /* make sure the jacobian has no negative diagonal element */
  for (row=0; row<nf; row++) dZU_dU[row][row]=max(0.0,dZU_dU[row][row]);
  /* the linearization is done assuming frozen Z */
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]+=1.0e0/gl->dt*dZU_dU[row][col];
    }
  }
}
#endif


void find_bdry_jacobian(np_t *np, gl_t *gl, long bdrytype, long lA, long lB,
                 int TYPELEVEL, sqmat_t DA, sqmat_t DB){
  set_matrix_to_identity(DA);
  set_matrix_to_zero(DB);
}


static void find_Aabs_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t Aabs){
  sqmat_t mattmp,R,L,lambdap;
  long flux;

  find_Lambda_from_jacvars(jacvars, metrics, lambdap);
  for (flux=0; flux<nf; flux++) lambdap[flux][flux]=fabs(lambdap[flux][flux]);
  find_conditioned_Lambda_absolute_from_jacvars(jacvars, metrics, EIGENVALCOND_DEFAULT, lambdap);
  find_Linv_from_jacvars(jacvars, metrics, R);
  find_L_from_jacvars(jacvars, metrics, L);
  multiply_diagonal_matrix_and_matrix(lambdap,L,mattmp);
  multiply_matrix_and_matrix(R,mattmp,Aabs);
}


void find_dFstar_dUstar_FDS_interface(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C){
  sqmat_t Ap1,Ap0;
  sqmat_t Aabs;
  jacvars_t jacvarsL,jacvarsR,jacvarsp1h;
  double Omegap0,Omegap1,Omegap1h;
  long lp0,lp1;
  metrics_t metricsp1h;
  long row,col;

  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);

  Omegap1h=(_Omega(np[lp0],gl)+_Omega(np[lp1],gl))*0.5e0;
  Omegap0=_Omega(np[lp0],gl);
  Omegap1=_Omega(np[lp1],gl);

  find_metrics_at_interface(np,gl,lp0,lp1,theta,&metricsp1h);
  find_jacvars(np[lp0], gl, metricsp1h, theta, &jacvarsL);
  find_jacvars(np[lp1], gl, metricsp1h, theta, &jacvarsR);
  find_jacvars_at_interface_arith_average(jacvarsL,jacvarsR,gl,theta,&jacvarsp1h);
  //find_jacvars_at_interface_Roe_average(jacvarsL,jacvarsR,gl,theta,&jacvarsp1h);

  find_Aabs_from_jacvars(jacvarsp1h, metricsp1h, Aabs);
  find_dFstar_dUstar(np[l], gl, theta, Ap0);
  find_dFstar_dUstar(np[lp1], gl, theta, Ap1);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=0.5*(Ap0[row][col]+Omegap1h/Omegap0*Aabs[row][col]);
      C[row][col]=0.5*(Ap1[row][col]-Omegap1h/Omegap1*Aabs[row][col]);
    }
  }

#if (defined(_RESTIME_CDF))
  sqmat_t B1,C1;
  find_dFstarxt_dUstar_interface(np, gl, theta, l, jacvarsp1h, metricsp1h, B1, C1);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]+=B1[row][col];
      C[row][col]+=C1[row][col];
    }
  }
#endif  
}


static void find_Aplus_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t Aplus){
  sqmat_t Lambdap,Lambda,mattmp,Linv,L;
  long flux;

  find_conditioned_Lambda_absolute_from_jacvars(jacvars, metrics, EIGENVALCOND_DEFAULT, Lambdap);
  find_Lambda_from_jacvars(jacvars, metrics, Lambda);
  for (flux=0; flux<nf; flux++) Lambdap[flux][flux]=0.5*(Lambda[flux][flux]+Lambdap[flux][flux]);
  find_Linv_from_jacvars(jacvars, metrics, Linv);
  find_L_from_jacvars(jacvars, metrics, L);
  multiply_matrix_and_matrix(Lambdap,L,mattmp);
  multiply_matrix_and_matrix(Linv,mattmp,Aplus);
}


static void find_Aminus_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t Aminus){
  sqmat_t Lambdap,Lambda,mattmp,Linv,L;
  long flux;

  find_conditioned_Lambda_absolute_from_jacvars(jacvars, metrics, EIGENVALCOND_DEFAULT, Lambdap);
  find_Lambda_from_jacvars(jacvars, metrics, Lambda);
  for (flux=0; flux<nf; flux++) Lambdap[flux][flux]=0.5*(Lambda[flux][flux]-Lambdap[flux][flux]);
  find_Linv_from_jacvars(jacvars, metrics, Linv);
  find_L_from_jacvars(jacvars, metrics, L);
  multiply_matrix_and_matrix(Lambdap,L,mattmp);
  multiply_matrix_and_matrix(Linv,mattmp,Aminus);
}


void find_dFstar_dUstar_FVS_interface(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C){
  jacvars_t jacvarsp1,jacvarsp0;
  sqmat_t Aminusp1,Aplusp0;
  long lp0,lp1,row,col;
  metrics_t metricsp1h;

  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);

  find_metrics_at_interface(np,gl,lp0,lp1,theta,&metricsp1h);
  find_jacvars(np[lp1],gl,metricsp1h,theta,&jacvarsp1);
  find_jacvars(np[lp0],gl,metricsp1h,theta,&jacvarsp0);

  find_Aplus_from_jacvars(jacvarsp0,metricsp1h,Aplusp0);
  find_Aminus_from_jacvars(jacvarsp1,metricsp1h,Aminusp1);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=metricsp1h.Omega/_Omega(np[lp0],gl)*Aplusp0[row][col];
      C[row][col]=metricsp1h.Omega/_Omega(np[lp1],gl)*Aminusp1[row][col];
    }
  }

#if (defined(_RESTIME_CDF))
  sqmat_t B1,C1;
  jacvars_t jacvarsp1h;
  find_jacvars_at_interface_arith_average(jacvarsp0,jacvarsp1,gl,theta,&jacvarsp1h);
  find_dFstarxt_dUstar_interface(np, gl, theta, l, jacvarsp1h, metricsp1h, B1, C1);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]+=B1[row][col];
      C[row][col]+=C1[row][col];
    }
  }
#endif  

}


static void find_A_plus_minus_from_jacvars_FVSplus(np_t *np, gl_t *gl, long theta, long l, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics, sqmat_t Aplusp0, sqmat_t Aminusp1){
  sqmat_t Lambdap0,mattmp,Linvp0,Lp0;
  sqmat_t Lambdap1,Linvp1,Lp1,Lambdaplusp0,Lambdaminusp1;
  long flux;
  flux_t tmpplusp0,tmpminusp1;

  find_Lambda_from_jacvars(jacvarsp0, metrics, Lambdap0);
  find_Lambda_from_jacvars(jacvarsp1, metrics, Lambdap1);
  set_matrix_to_zero(Lambdaplusp0);
  set_matrix_to_zero(Lambdaminusp1);
  for (flux=0; flux<nf; flux++) Lambdaplusp0[flux][flux]=0.5*(Lambdap0[flux][flux]+fabs(Lambdap0[flux][flux]));
  for (flux=0; flux<nf; flux++) Lambdaminusp1[flux][flux]=0.5*(Lambdap1[flux][flux]-fabs(Lambdap1[flux][flux]));

#if (defined(_RESTIME_CDF))
  sqmat_t Lambdaxtp0,Lambdaxtp1;
  long lp0,lp1;
  lp0=l;
  lp1=_al(gl,l,theta,+1);
  find_Lambdaxt(np, gl, lp0, jacvarsp0, metrics, Lambdaxtp0);
  find_Lambdaxt(np, gl, lp1, jacvarsp1, metrics, Lambdaxtp1);
  for (flux=0; flux<nf; flux++) {
    //don't make this function of Deltaxt 
    Lambdaplusp0[flux][flux]-=0.5*(Lambdaxtp0[flux][flux]+fabs(Lambdaxtp0[flux][flux]));
    Lambdaminusp1[flux][flux]-=0.5*(Lambdaxtp1[flux][flux]-fabs(Lambdaxtp1[flux][flux]));
  }
#endif
  for (flux=0; flux<nf; flux++){
    tmpplusp0[flux]=Lambdaplusp0[flux][flux];
    tmpminusp1[flux]=Lambdaminusp1[flux][flux];
  }
  for (flux=0; flux<nf;  flux++) {
    Lambdaplusp0[flux][flux]=max(0.0,tmpplusp0[flux])+max(0.0,tmpminusp1[flux]);
    Lambdaminusp1[flux][flux]=min(0.0,tmpminusp1[flux])+min(0.0,tmpplusp0[flux]);
  }

  condition_Lambda_plus_minus(np, gl, l, theta, jacvarsp0, jacvarsp1, metrics, EIGENVALCOND_DEFAULT, Lambdaplusp0,Lambdaminusp1);

  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_L_from_jacvars(jacvarsp0, metrics, Lp0);
  multiply_matrix_and_matrix(Lambdaplusp0,Lp0,mattmp);
  multiply_matrix_and_matrix(Linvp0,mattmp,Aplusp0);

  find_Linv_from_jacvars(jacvarsp1, metrics, Linvp1);
  find_L_from_jacvars(jacvarsp1, metrics, Lp1);
  multiply_matrix_and_matrix(Lambdaminusp1,Lp1,mattmp);
  multiply_matrix_and_matrix(Linvp1,mattmp,Aminusp1);
}


void find_dFstar_dUstar_FVSplus_interface(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C){
  jacvars_t jacvarsp1,jacvarsp0;
  sqmat_t Aminusp1,Aplusp0;
  long lp0,lp1,row,col;
  metrics_t metricsp1h;

  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);

  find_metrics_at_interface(np,gl,lp0,lp1,theta,&metricsp1h);
  find_jacvars(np[lp1],gl,metricsp1h,theta,&jacvarsp1);
  find_jacvars(np[lp0],gl,metricsp1h,theta,&jacvarsp0);

  find_A_plus_minus_from_jacvars_FVSplus(np,gl,theta,l,jacvarsp0, jacvarsp1, metricsp1h, Aplusp0, Aminusp1);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=metricsp1h.Omega/_Omega(np[lp0],gl)*Aplusp0[row][col];
      C[row][col]=metricsp1h.Omega/_Omega(np[lp1],gl)*Aminusp1[row][col];
    }
  }
}


static void find_A_plus_minus_from_jacvars_FDSplus(np_t *np, gl_t *gl, long theta, long l, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics, sqmat_t Aplusp0, sqmat_t Aminusp1){
  sqmat_t Lambdap0,Lambdap1h,mattmp,Linvp0,Lp0;
  sqmat_t Lambdap1,Linvp1,Lp1,Lambdaplusp0,Lambdaminusp1;
  long flux;
  flux_t tmpplusp0,tmpminusp1;
  jacvars_t jacvarsp1h;

  find_jacvars_at_interface_arith_average(jacvarsp0,jacvarsp1,gl,theta,&jacvarsp1h);

  find_Lambda_from_jacvars(jacvarsp0, metrics, Lambdap0);
  find_Lambda_from_jacvars(jacvarsp1, metrics, Lambdap1);
  find_Lambda_from_jacvars(jacvarsp1h, metrics, Lambdap1h);
  set_matrix_to_zero(Lambdaplusp0);
  set_matrix_to_zero(Lambdaminusp1);
  for (flux=0; flux<nf; flux++) Lambdaplusp0[flux][flux]=0.5*(Lambdap0[flux][flux]+fabs(Lambdap1h[flux][flux]));
  for (flux=0; flux<nf; flux++) Lambdaminusp1[flux][flux]=0.5*(Lambdap1[flux][flux]-fabs(Lambdap1h[flux][flux]));

#if (defined(_RESTIME_CDF))
  sqmat_t Lambdaxt;
  long lp0,lp1;
  lp0=l;
  lp1=_al(gl,l,theta,+1);
  find_Lambdaxt_interface(np, gl, lp0, lp1, theta, jacvarsp1h, metrics, Lambdaxt);
  for (flux=0; flux<nf; flux++) {
    //don't make this function of Deltaxt 
    Lambdaplusp0[flux][flux]-=0.5*(Lambdaxt[flux][flux]+fabs(Lambdaxt[flux][flux]));
    Lambdaminusp1[flux][flux]-=0.5*(Lambdaxt[flux][flux]-fabs(Lambdaxt[flux][flux]));
  }
#endif
  for (flux=0; flux<nf; flux++){
    tmpplusp0[flux]=Lambdaplusp0[flux][flux];
    tmpminusp1[flux]=Lambdaminusp1[flux][flux];
  }
  for (flux=0; flux<nf;  flux++) {
    Lambdaplusp0[flux][flux]=max(0.0,tmpplusp0[flux])+max(0.0,tmpminusp1[flux]);
    Lambdaminusp1[flux][flux]=min(0.0,tmpminusp1[flux])+min(0.0,tmpplusp0[flux]);
  }

  condition_Lambda_plus_minus(np, gl, l, theta, jacvarsp0, jacvarsp1, metrics, EIGENVALCOND_DEFAULT, Lambdaplusp0,Lambdaminusp1);

  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_L_from_jacvars(jacvarsp0, metrics, Lp0);
  multiply_matrix_and_matrix(Lambdaplusp0,Lp0,mattmp);
  multiply_matrix_and_matrix(Linvp0,mattmp,Aplusp0);

  find_Linv_from_jacvars(jacvarsp1, metrics, Linvp1);
  find_L_from_jacvars(jacvarsp1, metrics, Lp1);
  multiply_matrix_and_matrix(Lambdaminusp1,Lp1,mattmp);
  multiply_matrix_and_matrix(Linvp1,mattmp,Aminusp1);
}


void find_dFstar_dUstar_FDSplus_interface(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C){
  jacvars_t jacvarsp1,jacvarsp0;
  sqmat_t Aminusp1,Aplusp0;
  long lp0,lp1,row,col;
  metrics_t metricsp1h;

  lp0=_al(gl,l,theta,+0);
  lp1=_al(gl,l,theta,+1);

  find_metrics_at_interface(np,gl,lp0,lp1,theta,&metricsp1h);
  find_jacvars(np[lp1],gl,metricsp1h,theta,&jacvarsp1);
  find_jacvars(np[lp0],gl,metricsp1h,theta,&jacvarsp0);

  find_A_plus_minus_from_jacvars_FDSplus(np,gl,theta,l,jacvarsp0, jacvarsp1, metricsp1h, Aplusp0, Aminusp1);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]=metricsp1h.Omega/_Omega(np[lp0],gl)*Aplusp0[row][col];
      C[row][col]=metricsp1h.Omega/_Omega(np[lp1],gl)*Aminusp1[row][col];
    }
  }
}


static void condition_dFstar_dUstar_interface(sqmat_t A, sqmat_t B){
/*  long row;
  // this makes the integration process more robust in some cases for some (not yet fully understood) reason; 
  // note: need to look into this more closely later
  for (row=0; row<ns; row++) {
    B[fluxet][row]=max(0.0,B[fluxet][row]);
    A[fluxet][row]=min(0.0,A[fluxet][row]);
  }
*/ 
}


void add_dFstar_dUstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C){
  sqmat_t Bp1h,Cp1h,Am1h,Bm1h;
  long row,col;

  set_matrix_to_zero(Bp1h);
  set_matrix_to_zero(Cp1h);
  set_matrix_to_zero(Am1h);
  set_matrix_to_zero(Bm1h);

  switch (gl->cycle.resconv.CONVJACOBIAN){ 
    case CONVJACOBIAN_FVS:
      find_dFstar_dUstar_FVS_interface(np, gl, theta, l, Bp1h, Cp1h);
      find_dFstar_dUstar_FVS_interface(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);
    break;
    case CONVJACOBIAN_FVSPLUS:
      find_dFstar_dUstar_FVSplus_interface(np, gl, theta, l, Bp1h, Cp1h);
      find_dFstar_dUstar_FVSplus_interface(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);
    break;
    case CONVJACOBIAN_FDS:
      find_dFstar_dUstar_FDS_interface(np, gl, theta, l, Bp1h, Cp1h);
      find_dFstar_dUstar_FDS_interface(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);
    break;
    case CONVJACOBIAN_FDSPLUS:
      find_dFstar_dUstar_FDSplus_interface(np, gl, theta, l, Bp1h, Cp1h);
      find_dFstar_dUstar_FDSplus_interface(np, gl, theta, _al(gl,l,theta,-1), Am1h, Bm1h);
    break;
    default:
      fatal_error("gl->cycle.resconv.CONVJACOBIAN can not be set to %d.",gl->cycle.resconv.CONVJACOBIAN);
  }
  condition_dFstar_dUstar_interface(Bp1h,Cp1h);
  condition_dFstar_dUstar_interface(Am1h,Bm1h);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      A[row][col]+=-fact*Am1h[row][col];
      B[row][col]+=fact*(Bp1h[row][col]-Bm1h[row][col]);
      C[row][col]+=fact*Cp1h[row][col];
    }
  }

}


void find_TDMA_jacobians_conservative(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C){
  sqmat_t B1,C1;
  long row,col;
  double weightp0_default,weightp0_convection;
  weightp0_default=1.0;
  weightp0_convection=1.0;
#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
  weightp0_default=1.0-gl->cycle.restime.weightm1_trapezoidal_default;
  weightp0_convection=1.0-gl->cycle.restime.weightm1_trapezoidal_convection;
#endif
#ifdef _RESTIME_TRAPEZOIDAL_MUSCL
  weightp0_default=1.0;
  weightp0_convection=0.5;
#endif

  set_matrix_to_zero(B);
  set_matrix_to_zero(C);

  if (_FLUID_CONVECTION) {
    
    set_matrix_to_zero(B1);
    set_matrix_to_zero(C1);

    switch (gl->cycle.resconv.CONVJACOBIAN){ 
      case CONVJACOBIAN_FVS:
        find_dFstar_dUstar_FVS_interface(np, gl, theta, l, B1, C1);
      break;
      case CONVJACOBIAN_FVSPLUS:
        find_dFstar_dUstar_FVSplus_interface(np, gl, theta, l, B1, C1);
      break;
      case CONVJACOBIAN_FDS:
        find_dFstar_dUstar_FDS_interface(np, gl, theta, l, B1, C1);
      break;
      case CONVJACOBIAN_FDSPLUS:
        find_dFstar_dUstar_FDSplus_interface(np, gl, theta, l, B1, C1);
      break;
      default:
        fatal_error("gl->cycle.resconv.CONVJACOBIAN can not be set to %d.",gl->cycle.resconv.CONVJACOBIAN);
    }

    condition_dFstar_dUstar_interface(B1,C1);

    for (row=0; row<nf; row++){
      for (col=0; col<nf; col++){
        B[row][col]=weightp0_convection*B1[row][col];
        C[row][col]=weightp0_convection*C1[row][col];
      }
    }

  }



  if (_FLUID_DIFFUSION){
    find_Kstar_dG_dUstar_interface(np, gl, theta, l, B1, C1);
    for (row=0; row<nf; row++){
      for (col=0; col<nf; col++){
        B[row][col]+=weightp0_default*B1[row][col];
        C[row][col]+=weightp0_default*C1[row][col];
      }
    }
  }

#ifdef _FLUID_PLASMA
  find_Dstar_interface(np, gl, theta, l, B1, C1);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      B[row][col]+=weightp0_default*B1[row][col];
      C[row][col]+=weightp0_default*C1[row][col];
    }
  }
#endif

}


void add_TDMA_jacobians_non_conservative(np_t *np, gl_t *gl, long theta, long l, sqmat_t A, sqmat_t B, sqmat_t C, bool SOURCETERM){
  double weightp0_default;
  weightp0_default=1.0;
#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
  weightp0_default=1.0-gl->cycle.restime.weightm1_trapezoidal_default;
#endif
#ifdef _RESTIME_TRAPEZOIDAL_MUSCL
  weightp0_default=1.0;
#endif

#ifdef UNSTEADY
  if (SOURCETERM) add_Z_dUstardt_dUstar_to_TDMA(np, gl, l, B);
#endif   

  if (SOURCETERM && _FLUID_SOURCE) add_dSstar_dUstar_to_TDMA(np, gl, l, weightp0_default, B);

#ifdef _FLUID_PLASMA
  add_Ystar_dH_dUstar_to_TDMA(np, gl, theta, l, weightp0_default, A, B, C);
#endif

// need to make sure that the diagonal elements of B are positive for stability purposes
// ?????
//  for (flux=0; flux<nf; flux++) B[flux][flux]=max(0.0,B[flux][flux]);
}


void add_TDMA_jacobian_diagonally_dominant(np_t *np, gl_t *gl, long theta, long l, sqmat_t B){
  long dim;
  sqmat_t garb;
  double weightp0_default,weightp0_convection;
  weightp0_default=1.0;
  weightp0_convection=1.0;
#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
  weightp0_default=1.0-gl->cycle.restime.weightm1_trapezoidal_default;
  weightp0_convection=1.0-gl->cycle.restime.weightm1_trapezoidal_convection;
#endif
#ifdef _RESTIME_TRAPEZOIDAL_MUSCL
  weightp0_default=1.0;
  weightp0_convection=0.5;
#endif
  set_matrix_to_zero(garb);

  for (dim=0; dim<nd; dim++){
    if (dim!=theta){
      if (_FLUID_CONVECTION) add_dFstar_dUstar_to_TDMA(np, gl, dim, l, weightp0_convection, garb, B, garb);
      if (_FLUID_DIFFUSION) add_Kstar_dG_dUstar_to_TDMA(np, gl, dim, l, weightp0_default, garb, B, garb);
#ifdef _FLUID_PLASMA
      add_Ystar_dH_dUstar_to_TDMA(np, gl, dim, l, weightp0_default, garb, B, garb);
      add_Dstar_to_TDMA(np, gl, dim, l, weightp0_default, garb, B, garb);
#endif
    }
  }
}

