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


double _tauvt(np_t np, gl_t *gl){
  double N_tot,N_O,T,tauvt;
  long spec;
  
  T=_T(np,gl);
  tauvt=1.0e99;
  switch (gl->model.fluid.N2VIBMODEL){ 
    case N2VIBMODEL_MACHERET:
      N_O=0.0;
#ifdef specO
      N_O=_rhok(np,specO)*calA/_calM(specO);
#endif
      N_tot=0.0;
      for (spec=0; spec<ns; spec++) N_tot+=_rhok(np,spec)*calA/_calM(spec);
      tauvt=1.0/(N_tot*7E-16*exp(-141.0/pow(T,1.0/3.0))
            +N_O*5E-18*exp(-128.0/sqrt(T)));
    break;
    case N2VIBMODEL_MILLIKAN:
      fatal_error("Millikan N2 vibrational energy model not implemented in this version.");
    break;
    default:
      fatal_error("N2 vibrational model must be set to Macheret.");
  }
  return(tauvt);
}


static void find_Svib(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  
  for (flux=0; flux<nf; flux++) S[flux]=0.0e0;
  assert_np(np[l],_tauvt(np[l],gl)!=0.0);
  S[fluxev]=(_rhok(np[l],specN2)/_tauvt(np[l],gl)*(_evzero(np[l],gl)-_ev(np[l])));
}


void find_Schem(np_t *np, gl_t *gl, long l, flux_t S){
  spec_t W,rhok;
  long flux,spec;
  double EoverN;

  for (flux=0; flux<nf; flux++){
    S[flux]=0.0;
  }
  EoverN=0.0;
  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);
  find_W(rhok, _T(np[l],gl), _T(np[l],gl), _Tv(np[l]), EoverN, _Qbeam(np[l],gl), W);
  for (spec=0; spec<ns; spec++){
    S[spec]=W[spec];
  }
  S[fluxev]=W[specN2]*_ev(np[l]);
}


void find_Sstar(np_t *np, gl_t *gl, long l, flux_t S){
  flux_t Schem,St_norm,St_comp,Saxisymmetric,Svib;
  long flux;

  if (gl->model.fluid.REACTING) find_Schem(np,gl,l,Schem);
    else set_vector_to_zero(Schem);
  if (gl->model.fluid.TURBSOURCE) find_Stnorm(np, gl, l, St_norm);
    else set_vector_to_zero(St_norm);
  if (gl->model.fluid.TURBSOURCE) find_Stcomp(np, gl, l, St_comp);
    else set_vector_to_zero(St_comp);
  find_Svib(np, gl, l, Svib);
  find_Saxi(np, gl, l, Saxisymmetric);

  for (flux=0; flux<nf; flux++){
    S[flux]=_Omega(np[l],gl)*(
            Saxisymmetric[flux]+St_norm[flux]+
            St_comp[flux]+Schem[flux]+Svib[flux]
            );
  }
}


static void find_dSvib_dU(np_t *np, gl_t *gl, long l, sqmat_t dSvibdU){
  long row,col,spec,dim,flux;
  spec_t dTdrhok;
  double dTdrhoetstar;
  dim_t dTdrhoV;
  flux_t dTdU;
  double fact_evzero;

  /* for good linearization, fix fact_evzero to 1.0 ; if not stable, try fact_evzero=0.0; */
  fact_evzero=1.0;

  for (row=0; row<nf; row++) {
    for (col=0; col<nf; col++) {
      dSvibdU[row][col]=0.0e0;
    }
  }
  find_dT_dx(np[l], gl, &dTdrhoetstar, dTdrhok, dTdrhoV);
  for (flux=0; flux<nf; flux++) dTdU[flux]=0.0;
  for (spec=0; spec<ns; spec++) dTdU[spec]=dTdrhok[spec];
  for (dim=0; dim<nd; dim++) dTdU[ns+dim]=dTdrhoV[dim];
  dTdU[fluxet]=dTdrhoetstar;
  
  /* wrt rhoN2*ev */
  dSvibdU[fluxev][fluxev]=-1.0/_tauvt(np[l],gl);

  /* wrt rhoN2*evzero */
  dSvibdU[fluxev][specN2]=fact_evzero*1.0/_tauvt(np[l],gl)*_evzero(np[l],gl);
  for (flux=0; flux<nf; flux++)
    dSvibdU[fluxev][flux]+=fact_evzero*1.0/_tauvt(np[l],gl)*_rhok(np[l],specN2)*_dev_dTv_from_Tv(_T(np[l],gl))*dTdU[flux];

  /* wrt to (1/tauvt) */
  if (gl->model.fluid.N2VIBMODEL==N2VIBMODEL_MACHERET) {
 
    for (spec=0; spec<ns; spec++)
      dSvibdU[fluxev][spec]+=
        -_rhok(np[l],specN2)*(_ev(np[l])-fact_evzero*_evzero(np[l],gl))
          *7E-16*exp(-141.0/pow(_T(np[l],gl),1.0/3.0))*calA/_calM(spec);
#ifdef specO
    dSvibdU[fluxev][specO]+=
        -_rhok(np[l],specN2)*(_ev(np[l])-fact_evzero*_evzero(np[l],gl))*5E-18
        *exp(-128.0/pow(_T(np[l],gl),1.0/2.0))*calA/_calM(specO);
#endif

    for (flux=0; flux<nf; flux++){
      for (spec=0; spec<ns; spec++)
        dSvibdU[fluxev][flux]+=
        -_rhok(np[l],specN2)*(_ev(np[l])-fact_evzero*_evzero(np[l],gl))*7E-16*exp(-141.0/pow(_T(np[l],gl),1.0/3.0))
        *_rhok(np[l],spec)*calA/_calM(spec)
        *(141.0/3.0*pow(_T(np[l],gl),-4.0/3.0))*dTdU[flux];
#ifdef specO
      dSvibdU[fluxev][flux]+=
        -_rhok(np[l],specN2)*(_ev(np[l])-fact_evzero*_evzero(np[l],gl))*5E-18*exp(-128.0/pow(_T(np[l],gl),1.0/2.0))
        *_rhok(np[l],specO)*calA/_calM(specO)
        *(128.0/2.0*pow(_T(np[l],gl),-3.0/2.0))*dTdU[flux];
#endif
    }
    dSvibdU[fluxev][fluxtke]=-dSvibdU[fluxev][fluxet];
    dSvibdU[fluxev][fluxev]+=-dSvibdU[fluxev][fluxet];
  }

  /* wrt to (1/tauvt) */
  if (gl->model.fluid.N2VIBMODEL==N2VIBMODEL_MILLIKAN) {
    fatal_error("The Millikan model for N2 Vibration is not implemented in this version.");
  }
}


static void find_dSchem_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU){
  long k,s;
  spec_t mu,rhok,dWdT,dWdTe,dWdTv,dWdQbeam;
  spec2_t dWdrhok;
  flux_t dTdU,dTvdU,dTedU;
  double ev,devdTv,EoverN;
  flux_t W;
    
  EoverN=0.0;
  for (s=0; s<ns; s++) rhok[s]=_rhok(np[l],s);
  for (s=0; s<ns; s++) mu[s]=_mu(np,gl,l,s);

  find_dW_dx(rhok, mu, _T(np[l],gl), _T(np[l],gl), _Tv(np[l]), EoverN, _Qbeam(np[l],gl),
           dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam);
  find_dT_dU(np[l], gl, dTdU);
  find_dTv_dU(np[l], gl, dTvdU);
  find_dT_dU(np[l], gl, dTedU);
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

  find_W(rhok, _T(np[l],gl), _T(np[l],gl), _Tv(np[l]), EoverN, _Qbeam(np[l],gl), W);

//  S[fluxev]=W[specN2]*_ev(np[l]);
// here add the contribution to the jacobian coming from the chemical source term in the nitrogen vib energy
// equation

  ev=_ev(np[l]);
  devdTv=_dev_dTv_from_Tv(_Tv(np[l]));

  for (s=0; s<ns; s++)  dS_dU[fluxev][s]+=dWdrhok[specN2][s]*ev;
  for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=dWdT[specN2]*dTdU[s]*ev;
  for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=dWdTv[specN2]*dTvdU[s]*ev;
  for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=dWdTe[specN2]*dTedU[s]*ev;
  for (s=0; s<nf; s++)  dS_dU[fluxev][s]+=devdTv*dTvdU[s]*W[specN2];

}


void test_dSchem_dU(np_t *np, gl_t *gl, long l){
  long spec;
  double Estar;
  spec_t rhok,mu;
  Estar=0.0;
  for (spec=0; spec<ns; spec++) rhok[spec]=_rhok(np[l],spec);
  for (spec=0; spec<ns; spec++) mu[spec]=_mu(np,gl,l,spec);

  test_dW_dx(gl->cycle.fluid.Uref, rhok, mu,_T(np[l],gl), _T(np[l],gl), _Tv(np[l]), Estar, _Qbeam(np[l],gl));
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

