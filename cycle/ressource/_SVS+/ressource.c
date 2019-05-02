// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2016-2018 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <cycle/ressource/_ressource.h>
#include <cycle/share/cycle_share.h>


void add_Sstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  long l;
  long flux,row,col;
  metrics_t metrics;
  jacvars_t jacvarsp0,jacvarsm1;
  flux_t S,LUp0,LSm1,LSp0,LUm1,fluxtmp,Splus,Sminus;
  sqmat_t Lambdaminus,Lambdaplus,Linvm1,Linvp0,Lp0,Lm1;

  if (theta==0) {
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      find_Sstar(np,gl,l,S);

      find_metrics_at_node(np, gl, l, theta, &metrics);
      find_jacvars(np[l], gl, metrics, theta, &jacvarsp0);
#ifdef UNSTEADY
      find_jacvars_from_U(np[l].bs->Um1, metrics, gl, theta, &jacvarsm1);
#else
      jacvarsm1=jacvarsp0;
#endif


      /* need to find LSp0,  LUp0, LUm1, Linvm1, Linvp0, Lp0, Lm1 */

      find_L_from_jacvars(jacvarsp0, metrics, Lp0);
      find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
      find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0);

      find_L_from_jacvars(jacvarsm1, metrics, Lm1);
      find_Linv_from_jacvars(jacvarsm1, metrics, Linvm1);
      find_LUstar_from_jacvars(jacvarsm1, metrics, LUm1);
      
      multiply_matrix_and_vector(Lp0,S,LSp0);
      multiply_matrix_and_vector(Lm1,S,LSm1);

      for (row=0; row<nf; row++){
        for (col=0; col<nf; col++){
          Lambdaminus[row][col]=0.0;
          Lambdaplus[row][col]=0.0;
        }
        Lambdaminus[row][row]=0.5*(LSp0[row]/notzero(LUp0[row],1e-99)-fabs(LSp0[row]/notzero(LUp0[row],1e-99)));
        Lambdaplus[row][row]=0.5*(LSm1[row]/notzero(LUm1[row],1e-99)+fabs(LSm1[row]/notzero(LUm1[row],1e-99)));
      }



      multiply_diagonal_matrix_and_vector(Lambdaminus,LUp0,fluxtmp);
      multiply_matrix_and_vector(Linvp0,fluxtmp,Sminus);

      multiply_diagonal_matrix_and_vector(Lambdaplus,LUm1,fluxtmp);
      multiply_matrix_and_vector(Linvm1,fluxtmp,Splus);

#ifdef _RESSOURCE_LAMBDAMINUS_STORAGE
      for (flux=0; flux<nf; flux++) np[l].bs->Lambda_S_minus[flux]=Lambdaminus[flux][flux];
#endif
#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]-=(1.0-gl->cycle.restime.weightm1_trapezoidal_default)*(Splus[flux]+Sminus[flux]);
      for (flux=0; flux<nf; flux++) np[l].bs->trapezoidalm1_next[flux]-=gl->cycle.restime.weightm1_trapezoidal_default*(Splus[flux]+Sminus[flux]);
#else
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]-=(Splus[flux]+Sminus[flux]);
#endif

    }
  }
}


