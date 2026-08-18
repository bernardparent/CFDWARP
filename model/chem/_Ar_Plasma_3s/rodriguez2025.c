// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2025 Felipe Martin Rodriguez Fuentes

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

/*
This is the 3 species version of the argon mechanism: Ar(4s) and every reaction it takes
part in have been removed, so that the only remaining process is the ionization of ground
state argon by electron impact.  It is reaction 2 of the 4 species mechanism, renamed
reaction 1 here since it is the only one left; its spline control points are unchanged.
The electron energy lost to the inelastic processes other than ionization is accounted for
by the source term of the electron temperature transport model, following the -ne*N*k*L
term of Boeuf and Pitchford, Phys. Rev. E 51:1376, 1995.
*/

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/share/chem_share.h>
#include <model/share/model_share.h>
#include <exm.h>
#include "chem.h"

#define TeV 11604.52500617 /* Kelvin per eV */

const static bool REACTION[2]=
  {
   FALSE,/* reaction 0 */
   TRUE  /* reaction 1 */
  };


/* Reaction 1: e- + Ar -> Ar+ + e- + e-, ground ionization (15.8 eV)
   (reaction 2 of the 4 species and 8 species mechanisms) */
static double _kf1(np_t np, gl_t *gl, double Te){
  double kf1;
  double Te_control[] =
  {
    10.82646,
    10.96672,
    10.99798,
    11.06005,
    11.11007,
    11.19590,
    11.32563,
    11.52376,
    11.88015,
    13.52579,
    16.11810
  };
  double kf_control[] =
  {
    -62.83592,
    -48.78568,
    -44.60907,
    -34.45233,
    -29.62413,
    -25.43377,
    -22.37717,
    -20.26565,
    -18.58446,
    -16.15063,
    -15.47303
  };
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf1 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf1 ));
}



static double _averaged_rate_local(np_t np, gl_t *gl, double Te, spec_t rhok,  long react1, long react2, char* id, double (*_kf)(np_t, gl_t *, double)){
  double kfave,rate;
  averagedrates_id_type id2;
  strcpy(id2,id);
  strins("Ne",id2,0);
  rate=_kf(np,gl,Te);
  kfave=min(_averaged_rate(np,gl,id,rate),_averaged_rate(np,gl,id2,rhok[react1]*rhok[react2]*rate)/(rhok[react1]*rhok[react2]));
  return(kfave);
}


void find_W_rodriguez2025 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;
  spec_t N;

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);

  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
    W[k] = 0.0;
  }

  /* Ground ionization */
  if (REACTION[1]){
     add_to_W_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), N, W );
  }

}

void find_dW_dx_rodriguez2025 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv,
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */
  spec_t N;
  double dkfdTe,dkfdT,dkfdTv,dTe;

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);

  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }

  /* find properties needed by add_to_dW* functions in proper units */
  dTe=0.001*Te;

  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA; /* particules/cm^3 */
  }

  /* Ground ionization */
  if (REACTION[1]){
     dkfdTe = (_kf1(np,gl,Te+dTe)-_kf1(np,gl,Te))/dTe;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     add_to_dW_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
  }

}


void find_Qei_rodriguez2025(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){

  *Qei=0.0;

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);

  if (REACTION[1]) {
    add_to_Qei(gl,Te, specAr, _ionizationpot(specAr), _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), rhok, Qei);
  }

}


void find_dQei_dx_rodriguez2025(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){

  long spec;

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);

  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;

  if (REACTION[1]) {
    add_to_dQei(gl,Te, specAr, _ionizationpot(specAr),  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), 0.0,  rhok, dQeidrhok, dQeidTe);
  }

}




