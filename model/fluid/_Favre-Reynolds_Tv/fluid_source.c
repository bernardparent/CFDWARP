// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2005-2006 Bernard Parent

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
#include <model/emfield/_emfield.h>
#include <model/chem/_chem.h>
#include <model/thermo/_thermo.h>
#include <model/beam/_beam.h>
#include <model/_model.h>
#include <model/share/fluid_share.h>
#include <model/share/fluid_vib.h>
#include <model/share/chem_share.h>


void find_Schem(np_t *np, gl_t *gl, long l, flux_t S){
  spec_t W,rhok;
  vibdist_t vibdist;
  long flux,spec;
  double EoverN,Te;

  for (flux=0; flux<nf; flux++){
    S[flux]=0.0;
  }
  EoverN=0.0;
  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);
  switch (gl->model.fluid.TEMODEL){ 
    case TEMODEL_TEQUILIBRIUM:
      Te=_T(np[l],gl);
    break;
    case TEMODEL_TVEQUILIBRIUM:
      Te=_Tv(np[l],gl);
    break;
    default:
      Te=0.0;  // needed to avoid warning about unitialized Te
      fatal_error("Problem with TEMODEL in find_Schem().");
  }  
  
  find_W(np[l],gl, rhok, _T(np[l],gl), Te, _Tv(np[l],gl), EoverN, _Qbeam(np[l],gl), W);
  for (spec=0; spec<ns; spec++){
    S[spec]=W[spec];
  }
  find_vibdist_from_np(np[l],gl,&vibdist);
  for (spec=0; spec<ns; spec++){
    if (_isvib(spec)) S[fluxev]+=W[spec]*vibdist.evk[spec];
  }
}


void find_Sstar(np_t *np, gl_t *gl, long l, flux_t S){
  flux_t Schem,St_norm,St_comp,Saxisymmetric,Svib,Sheatforces;
  long flux;

  if (gl->model.fluid.REACTING) find_Schem(np,gl,l,Schem);
    else set_vector_to_zero(Schem);
  if (gl->model.fluid.TURBSOURCE) find_Stnorm(np, gl, l, St_norm);
    else set_vector_to_zero(St_norm);
  if (gl->model.fluid.TURBSOURCE) find_Stcomp(np, gl, l, St_comp);
    else set_vector_to_zero(St_comp);
  find_Svib_VT(np, gl, l, Svib);
  find_Saxi(np, gl, l, Saxisymmetric);
  find_Sheatforces(np, gl, l, Sheatforces);

  for (flux=0; flux<nf; flux++){
    S[flux]=_Omega(np[l],gl)*(
            Saxisymmetric[flux]+St_norm[flux]+
            St_comp[flux]+Schem[flux]+Svib[flux]+Sheatforces[flux]
            );
  }
}


static void find_dSvib_dU(np_t *np, gl_t *gl, long l, sqmat_t dSvibdU){
  flux_t dTdU;

  find_dT_dU(np[l], gl, dTdU);
  find_dSvib_VT_dU(np, gl, l, dTdU, dSvibdU);
}


static void find_dSchem_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU){
  long k,s;
  spec_t rhok,dWdT,dWdTe,dWdTv,dWdQbeam;
  spec2_t dWdrhok;
  flux_t dTdU,dTvdU,dTedU;
  double EoverN;
  flux_t W;
  vibdist_t vibdist;
    
  EoverN=0.0;
  for (s=0; s<ns; s++) rhok[s]=_rhok(np[l],s);

  find_dW_dx(np[l],gl, rhok, _T(np[l],gl), _T(np[l],gl), _Tv(np[l],gl), EoverN, _Qbeam(np[l],gl),
           dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam);
  find_dT_dU(np[l], gl, dTdU);
  find_dTv_dU(np[l], gl, dTvdU);
  switch (gl->model.fluid.TEMODEL){ 
    case TEMODEL_TEQUILIBRIUM:
      find_dT_dU(np[l], gl, dTedU);
    break;
    case TEMODEL_TVEQUILIBRIUM:
      find_dTv_dU(np[l], gl, dTedU);
    break;
    default:
      fatal_error("Problem with TEMODEL in find_dSchem_dU().");
  }  
  for ( k = 0; k < nf; k++ ){
    for ( s = 0; s < nf; s++ ){
      dS_dU[k][s]=0.0e0;
    }
  }
	
  for ( k = 0; k < ns; k++ ){
    for (s=0; s<ns; s++)  dS_dU[k][s]+=dWdrhok[k][s];
    for (s=0; s<nf; s++)  dS_dU[k][s]+=dWdT[k]*dTdU[s];
    for (s=0; s<nf; s++)  dS_dU[k][s]+=dWdTv[k]*dTvdU[s];
    for (s=0; s<nf; s++)  dS_dU[k][s]+=dWdTe[k]*dTedU[s];
  }

  find_W(np[l],gl,rhok, _T(np[l],gl), _T(np[l],gl), _Tv(np[l],gl), EoverN, _Qbeam(np[l],gl), W);

//  S[fluxev]=sum_k W[k]*evk;
// here add the contribution to the jacobian coming from the chemical source term in the vib energy
// equation

  find_vibdist_from_np(np[l],gl,&vibdist);
  for (k=0; k<ns; k++){
    if (_isvib(k)){
      for (s=0; s<ns; s++)  dS_dU[fluxev][s]+=dWdrhok[k][s]*vibdist.evk[k];
      for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=dWdT[k]*dTdU[s]*vibdist.evk[k];
      for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=dWdTv[k]*dTvdU[s]*vibdist.evk[k];
      for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=dWdTe[k]*dTedU[s]*vibdist.evk[k];
#ifdef TEST
      for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=vibdist.cvk[k]*dTvdU[s]*W[k];
#else
      /* when species k is formed (W[k]>0) this is a gain term, left out of the Jacobian */
      if (W[k]<0.0) for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=vibdist.cvk[k]*dTvdU[s]*W[k];
#endif
    }
  }

#ifndef TEST
  for (k=0; k<nf; k++) dS_dU[k][k]=min(0.0e0,dS_dU[k][k]);
#endif
}


void test_dSchem_dU(np_t *np, gl_t *gl, long l){
  long spec;
  double Estar;
  spec_t rhok;
  Estar=0.0;
  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);

  test_dW_dx(np[l],gl, gl->cycle.fluid.Uref, rhok,_T(np[l],gl), _T(np[l],gl), _Tv(np[l],gl), Estar, _Qbeam(np[l],gl));
}


/* the vibrational energy source term tested by test_dSvib_dU(): V-T relaxation */
static void find_Svibstar(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  find_Svib_VT(np, gl, l, S);
  for (flux=0; flux<nf; flux++) S[flux]*=_Omega(np[l],gl);
}


void test_dSvib_dU(np_t *np, gl_t *gl, long l){
  sqmat_t dSvibdU,dSvibdUnum;
  find_dSvib_dU(np,gl,l,dSvibdU);
  find_numerical_jacobian_3(np, l, gl, &find_Svibstar, dSvibdUnum);
  test_dSvib_dU_display(gl, dSvibdU, dSvibdUnum);
}


void find_dSstar_dUstar(np_t *np, gl_t *gl, long l, sqmat_t dSstar_dUstar){
  long col,row;
  sqmat_t dSchemdU,dStnormdU,dStcompdU,dSaxidU,dSvibdU;

  if (gl->model.fluid.REACTING) find_dSchem_dU(np,gl,l,dSchemdU);
    else set_matrix_to_zero(dSchemdU);
  if (gl->model.fluid.TURBSOURCE) find_dStnorm_dU(np,gl,l,dStnormdU);
    else set_matrix_to_zero(dStnormdU);
  if (gl->model.fluid.TURBSOURCE) find_dStcomp_dU(np,gl,l,dStcompdU);
    else set_matrix_to_zero(dStcompdU);
  find_dSvib_dU(np,gl,l,dSvibdU);
  find_dSaxi_dU(np,gl,l,dSaxidU);
  
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dSstar_dUstar[row][col]=dSchemdU[row][col]+dSaxidU[row][col]
                 +dStnormdU[row][col]+dStcompdU[row][col]
                 +dSvibdU[row][col];
    }
  }

#ifndef TEST
  for (row=0; row<nf; row++) dSstar_dUstar[row][row]=min(0.0e0,dSstar_dUstar[row][row]);
#endif

}

