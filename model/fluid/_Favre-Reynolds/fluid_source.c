// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2000-2002, 2017 Bernard Parent

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

#include "fluid.h"
#include "fluid_source.h"
#include <model/metrics/_metrics.h>
#include <model/chem/_chem.h>
#include <model/thermo/_thermo.h>
#include <model/_model.h>
#include <model/emfield/_emfield.h>
#include <model/share/fluid_share.h>
#include <model/share/chem_share.h>


void find_Schem(np_t *np, gl_t *gl, long l, flux_t S){
  spec_t W,rhok;
  long flux,spec;
  double Estar;
  for (flux=0; flux<nf; flux++){
    S[flux]=0.0;
  }
  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);
  Estar=0.0;
  find_W(np[l],gl, rhok, _T(np[l],gl), _T(np[l],gl), _T(np[l],gl), Estar, _Qbeam(np[l],gl), W);
  for (spec=0; spec<ns; spec++){
    S[spec]=W[spec];
  }
}


void find_Sstar(np_t *np, gl_t *gl, long l, flux_t S){
  flux_t Schem,St_norm,St_comp,Saxisymmetric,Sheatforces;
  long flux;

  if (gl->model.fluid.REACTING) find_Schem(np,gl,l,Schem);
    else set_vector_to_zero(Schem);
  if (gl->model.fluid.TURBSOURCE) find_Stnorm(np, gl, l, St_norm);
    else set_vector_to_zero(St_norm);
  if (gl->model.fluid.TURBSOURCE) find_Stcomp(np, gl, l, St_comp);
    else set_vector_to_zero(St_comp);
  find_Saxi(np, gl, l, Saxisymmetric);
  find_Sheatforces(np, gl, l, Sheatforces);

  for (flux=0; flux<nf; flux++){
    S[flux]=_Omega(np[l],gl)*(
            Saxisymmetric[flux]+St_norm[flux]+
            St_comp[flux]+Schem[flux]+Sheatforces[flux]
            );
  }
    
}


static void find_dSchem_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU){
  long k,s;
  spec_t rhok,dWdT,dWdTe,dWdTv,dWdQbeam;
  spec2_t dWdrhok;
  flux_t dTdU;
  double Estar;

  Estar=0.0;  
  for (s=0; s<ns; s++) rhok[s]=_rhok(np[l],s);

  find_dW_dx(np[l],gl, rhok, _T(np[l],gl), _T(np[l],gl), _T(np[l],gl), Estar,_Qbeam(np[l],gl),
           dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam);
  find_dT_dU(np[l], gl, dTdU);
  for ( k = 0; k < nf; k++ ){
    for ( s = 0; s < nf; s++ ){
      dS_dU[k][s]=0.0e0;
    }
  }
	
  for ( k = 0; k < ns; k++ ){
    for (s=0; s<ns; s++)  dS_dU[k][s]+=dWdrhok[k][s];
    for (s=0; s<nf; s++)  dS_dU[k][s]+=dWdT[k]*dTdU[s];
    for (s=0; s<nf; s++)  dS_dU[k][s]+=dWdTv[k]*dTdU[s];
    for (s=0; s<nf; s++)  dS_dU[k][s]+=dWdTe[k]*dTdU[s];
  }
}


void test_dSchem_dU(np_t *np, gl_t *gl, long l){
  long spec;
  double Estar;
  spec_t rhok;
  Estar=0.0;
  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);
  test_dW_dx(np[l],gl, gl->cycle.fluid.Uref, rhok, _T(np[l],gl), _T(np[l],gl), _T(np[l],gl), Estar, _Qbeam(np[l],gl));
}


void find_dSstar_dUstar(np_t *np, gl_t *gl, long l, sqmat_t dSstar_dUstar){
  long col,row;
  sqmat_t dSchemdU,dStnormdU,dStcompdU,dSaxidU;

  if (gl->model.fluid.REACTING) find_dSchem_dU(np,gl,l,dSchemdU);
    else set_matrix_to_zero(dSchemdU);
  if (gl->model.fluid.TURBSOURCE) find_dStnorm_dU(np,gl,l,dStnormdU);
    else set_matrix_to_zero(dStnormdU);
  if (gl->model.fluid.TURBSOURCE) find_dStcomp_dU(np,gl,l,dStcompdU);
    else set_matrix_to_zero(dStcompdU);
  find_dSaxi_dU(np,gl,l,dSaxidU);
  
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dSstar_dUstar[row][col]=dSchemdU[row][col]+dSaxidU[row][col]
                 +dStnormdU[row][col]+dStcompdU[row][col];
    }
  }

#ifndef TEST
  for (row=0; row<nf; row++) dSstar_dUstar[row][row]=min(0.0e0,dSstar_dUstar[row][row]);
#endif
}
