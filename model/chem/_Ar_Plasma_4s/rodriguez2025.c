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

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/share/chem_share.h>
#include <model/share/model_share.h>
#include <exm.h>
#include "chem.h"

#define TeV 11604.52500617 /* Kelvin per eV */

const static bool REACTION[6]=
  {
   FALSE,/* reaction 0 */
   TRUE, /* reaction 1 */  
   TRUE, /* reaction 2 */  
   TRUE, /* reaction 3 */  
   TRUE, /* reaction 4 */    
   TRUE  /* reaction 5 */      
  };


/* Reaction 1: e- + Ar -> Ar(4S) + e- */
static double _kf1(np_t np, gl_t *gl, double Te){
  double kf1;
  double Te_control[] = 
  { 
      10.426842881798962,
      10.448523459349468,
      10.470214825626545,
      10.491909642440874,
      10.513597807717549,
      10.535276034154350,
      10.556989418911494,
      10.664874402038480,
      10.772731554624395,
      10.875529583951693,
      10.953649266438417,
      11.000116036607594,
      11.027047756209608,
      11.051076630287564,
      11.094221332741899,
      11.147475633462962,
      11.213571618121964,
      11.296357056900984,
      11.408899180737061,
      11.579718158153540,
      11.868814688682924,
      12.414796442316845,
      13.288117035527231
  };
  double kf_control[] = 
  {
      -56.869684956605269,
      -54.404672039572588,
      -52.103565751682126,
      -49.959327497916512,
      -47.964703967294902,
      -46.112343131835964,
      -44.383921834455208,
      -37.537754530580209,
      -33.018627724559451,
      -30.138261231283884,
      -28.349398518564950,
      -27.204099650650324,
      -26.399869134541088,
      -25.612349870737688,
      -24.462370524449220,
      -23.521465989760557,
      -22.696982097503433,
      -21.964006711281961,
      -21.300984395752920,
      -20.665875899910027,
      -20.001696422123771,
      -19.306608627984254,
      -18.751721637059706
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf1 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf1 ));
}

/* Reaction 1b: e- + Ar(4S) -> Ar + e- */
static double _kf1b(np_t np, gl_t *gl, double Te){
  double kf1b;
  double Te_control[] = 
  { 
      7.056565294775152,
      7.144008245290196,
      7.231442544679870,
      7.318875711639040,
      7.756082145842774,
      8.193285498178808,
      8.630481723462170,
      9.067682491454400,
      9.504885757995108,
      9.942080835870144,
      10.379285686831711,
      10.816478536289583,
      11.253681524477567,
      11.690886432854930,
      12.128082432505435,
      12.565281340689033,
      13.002483161927154,
      13.439679960957569,
      13.876880741357384
  };
  double kf_control[] = 
  {
      -22.824698350602993,
      -22.776030005269941,
      -22.725920550195134,
      -22.674436716075363,
      -22.400056003131319,
      -22.115611403066836,
      -21.824214853638153,
      -21.475621767143039,
      -21.174992446488098,
      -20.910727983962762,
      -20.633699698363252,
      -20.445577788802343,
      -20.389540584118510,
      -20.321349351090326,
      -20.234468770020630,
      -20.156639821034659,
      -20.096556313074224,
      -20.056970974462708,
      -20.041991021497978
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf1b = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf1b ));
}

/* Reaction 2: e- + Ar -> Ar+ + e- + e-, ground ionization (15.8 eV) */
static double _kf2(np_t np, gl_t *gl, double Te){
  double kf2;
  double Te_control[] = 
  { 
      10.875088281069146,
      10.912487679643895,
      10.949884469009515,
      11.011412129714213,
      11.045534526414871,
      11.072643764652126,
      11.102249934346027,
      11.136706126584885,
      11.176447640785092,
      11.305060536824511,
      11.403405587266139,
      11.530026084437790,
      11.701365725431661,
      11.902593205561555,
      12.041336990945240,
      12.651815562607892,
      12.934777153688277,
      13.240606963469691,
      13.567439517032776,
      13.917929978318670
  };
  double kf_control[] = 
  {
      -58.825682075352518,
      -55.332610629548185,
      -51.640435700618660,
      -43.897630696415902,
      -37.810011871997588,
      -33.538810177283501,
      -30.382531162623028,
      -27.974231402306742,
      -26.071651984746421,
      -23.320367638733110,
      -21.940639364555771,
      -20.590289473321992,
      -19.451855530092949,
      -18.586209020525686,
      -18.119979292620695,
      -17.219108819433202,
      -16.960114138457886,
      -16.745179847964220,
      -16.570789168410105,
      -16.435650161496891
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf2 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf2 ));
}

/* Reaction 3: e- + Ar(4S) -> Ar+ + e- + e-, stepwise ionization (4.2 eV) */
static double _kf3(np_t np, gl_t *gl, double Te){
  double kf3;
  double Te_control[] = 
  { 
      8.980244441846285,
      9.067683829844169,
      9.155122953647219,
      9.242563761130038,
      9.329994454365105,
      9.417438163639705,
      9.854645091237762,
      10.291841888744406,
      10.729043480407910,
      11.166240079684437,
      11.603439832281747,
      12.040644407581713,
      12.477843380721994,
      12.915041159107833,
      13.352240806449975,
      13.789439666545883
   };
  double kf_control[] = 
  {
      -47.456767300525378,
      -42.772291770219525,
      -38.271519966111072,
      -34.209659715623900,
      -30.870622744088788,
      -28.280463421464386,
      -21.639804099453745,
      -18.795494927954628,
      -17.320932434134033,
      -16.582841128545077,
      -16.177700230694612,
      -15.905277220722361,
      -15.698234662566422,
      -15.540146101096299,
      -15.425613691608950,
      -15.352604553275206
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf3 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf3 ));
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
  spec_t X, N, Nlimit;
  double R;

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
    W[k] = 0.0;
    Nlimit[k]=N[k];
    if (k==speceminus) Nlimit[k]=min(1e8,Nlimit[k]);
  }

  R=Rchem;

  /* Metastable excitation and inverse reaction */    
  
  if (REACTION[1]){
     if (gl->model.chem.RATE_AVERAGING) {
       if (gl->model.chem.EXCITED_STATES) {
          add_to_W_2r2p ( speceminus, specAr, specAr4S, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), N, W );
          if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4S, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf1b", &_kf1b), N, W );
       }
     } else {
       if (gl->model.chem.EXCITED_STATES) {
          add_to_W_2r2p ( speceminus, specAr, specAr4S, speceminus, _kf1(np,gl,Te), N, W);
          if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4S, specAr, speceminus, _kf1b(np,gl,Te), N, W);
       }
     } 
  }

  /* Ground ionization */

  if (REACTION[2]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), N, W );
     } else {
       add_to_W_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _kf2(np,gl,Te), N, W );
     }
  }

  /* Stepwise ionization */

  if (REACTION[3]){
     if (gl->model.chem.RATE_AVERAGING){
       if (gl->model.chem.EXCITED_STATES) add_to_W_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), N, W );
     } else {
       if (gl->model.chem.EXCITED_STATES) add_to_W_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _kf3(np,gl,Te), N, W );
     }
  }

  /* Penning ionization */  

  if (REACTION[4]){
     add_to_W_fw_2r3p ( specAr4S, specAr4S, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, W );
  }

  /* Ar(4s) quenching */

  if (REACTION[5]){
     add_to_W_fw_2r2p ( specAr4S, specAr, specAr, specAr, 2.3000E-15 * calA, 0.0, 0.0 * R, T, X, W );
  }

}

void find_dW_dx_rodriguez2025 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */
  spec_t X, N, Nlimit;
  double R;
  double dkfdTe,dkfdT,dkfdTv;
  
  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);  

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA;
    Nlimit[k]=N[k];
    if (k==speceminus) Nlimit[k]=min(1e8,Nlimit[k]);
  }

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
  R=Rchem;
  
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA; /* particules/cm^3 */
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }

  /* Ar electron impact reactions */    
  
  if (REACTION[1]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
        if (gl->model.chem.EXCITED_STATES) {
          add_to_dW_2r2p ( speceminus, specAr, specAr4S, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
          if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4S, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf1b", &_kf1b), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
        }
     } else {
       if (gl->model.chem.EXCITED_STATES) {
          add_to_dW_2r2p ( speceminus, specAr, specAr4S, speceminus, _kf1(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
          if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4S, specAr, speceminus, _kf1b(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
       }
     } 
  }

  if (REACTION[2]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       add_to_dW_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _kf2(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     }
  }

  if (REACTION[3]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       if (gl->model.chem.EXCITED_STATES) add_to_dW_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       if (gl->model.chem.EXCITED_STATES) add_to_dW_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _kf3(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     }
  }

  /* Penning ionization */

  if (REACTION[4]){
     add_to_dW_fw_2r3p ( specAr4S, specAr4S, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

  /* Ar(4s) quenching */

  if (REACTION[5]){
     add_to_dW_fw_2r2p ( specAr4S, specAr, specAr, specAr, 2.3000E-15 * calA, 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

}


void find_Qei_rodriguez2025(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
 
  *Qei=0.0;  

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);  

  if (REACTION[1]) {
     if (gl->model.chem.RATE_AVERAGING) {
        if (gl->model.chem.EXCITED_STATES) add_to_Qei(gl,Te, specAr, 11.6, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), rhok, Qei);
     } else {
        if (gl->model.chem.EXCITED_STATES) add_to_Qei(gl,Te, specAr, 11.6, _kf1(np,gl,Te), rhok, Qei); // Ar->Ar(4s) (11.6 eV)
     }
  }

  if (REACTION[2]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr, _ionizationpot(specAr), _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr, _ionizationpot(specAr), _kf2(np,gl,Te), rhok, Qei); // ground ionization 15.8 eV
     } 
  }

  if (REACTION[3]) {
     if (gl->model.chem.RATE_AVERAGING) {
       if (gl->model.chem.EXCITED_STATES) add_to_Qei(gl,Te, specAr4S, 4.2, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), rhok, Qei);
     } else {
       if (gl->model.chem.EXCITED_STATES) add_to_Qei(gl,Te, specAr4S, 4.2, _kf3(np,gl,Te), rhok, Qei); // Ar(4s)->Ar+ (15.8 - 11.6 = 4.2 eV)
     } 
  }

}


void find_dQei_dx_rodriguez2025(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){

  long spec;

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);  

  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  

  if (REACTION[1]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr, 11.6,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr, 11.6, _kf1(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // Ar->Ar(4s) (11.6 eV)
     } 
  }

  if (REACTION[2]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr, _ionizationpot(specAr),  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr, _ionizationpot(specAr), _kf2(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // ground ionization 15.8 eV
     } 
  }

  if (REACTION[3]) {
     if (gl->model.chem.RATE_AVERAGING) {
       if (gl->model.chem.EXCITED_STATES) add_to_dQei(gl,Te, specAr4S, 4.2,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       if (gl->model.chem.EXCITED_STATES) add_to_dQei(gl,Te, specAr4S, 4.2, _kf3(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // Ar(4s)->Ar+ (15.8 - 11.6 = 4.2 eV)
     } 
  }

}




