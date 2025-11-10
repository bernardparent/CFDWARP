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

const static bool REACTION[19]=
  {
   FALSE,/* reaction 0 */
   TRUE, /* reaction 1 */  
   TRUE, /* reaction 2 */  
   TRUE, /* reaction 3 */  
   TRUE, /* reaction 4 */      
   TRUE, /* reaction 5 */     
   TRUE, /* reaction 6 */     
   TRUE, /* reaction 7 */   
   TRUE, /* reaction 8 */  
   TRUE, /* reaction 9 */   
   TRUE, /* reaction 10 */     
   TRUE, /* reaction 11 */  
   TRUE, /* reaction 12 */    
   TRUE, /* reaction 13 */    
   TRUE, /* reaction 14 */    
   TRUE, /* reaction 15 */  
   TRUE, /* reaction 16 */    
   TRUE, /* reaction 17 */     
   TRUE  /* reaction 18 */ 
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

/* Reaction 2: e- + Ar -> Ar(4P) + e- */
static double _kf2(np_t np, gl_t *gl, double Te){
  double kf2;
  double Te_control[] = 
  { 
      10.545828497305541,
      10.551699185915687,
      10.563440141861250,
      10.581053479420989,
      10.617364187739453,
      10.639167146168779,
      10.660956549643897,
      10.682757845781653,
      10.704560252415627,
      10.726358764742853,
      10.756864296943597,
      10.848224687720862,
      10.939596057436621,
      11.030961979761358,
      11.122328610126058,
      11.213687453616906,
      11.305057679677573,
      11.684953836439020,
      12.119116023023208,
      12.553286698947188,
      12.987448812866520,
      13.421613769783100,
      13.855780564417739
  };
  double kf_control[] = 
  {
      -59.826320025802893,
      -57.162501118008493,
      -55.464446095053788,
      -53.664094074581968,
      -50.397975324848979,
      -48.613577246834659,
      -46.970503240798443,
      -45.440545950879297,
      -44.027964642682733,
      -42.725069990925633,
      -41.062433215526546,
      -36.985351957914851,
      -33.589553914801556,
      -29.310216477779111,
      -25.030192420311302,
      -23.114958684613242,
      -21.990086148535692,
      -20.009130457648123,
      -19.168231313391697,
      -18.741025136244875,
      -18.480135675123705,
      -18.309747783246063,
      -18.206699034430542
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf2 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf2 ));
}

/* Reaction 2b: e- + Ar(4P) -> Ar + e- */
static double _kf2b(np_t np, gl_t *gl, double Te){
  double kf2b;
  double Te_control[] = 
  { 
      7.056565294775152,
      7.144008245290196,
      7.231442544679870,
      7.318875711639040,
      7.406327400007066,
      8.105843771524167,
      8.805361589904287,
      9.504885757995108,
      10.204404733029481,
      10.903921770873028,
      11.603441952336190,
      12.302952321960307,
      13.002483161927154,
      13.702002108241009
  };
  double kf_control[] = 
  {
      -23.262089852566909,
      -23.211268875602414,
      -23.158409082031511,
      -23.103446052063450,
      -23.046265904683072,
      -22.473058367353161,
      -21.674318835632857,
      -20.934714293265319,
      -20.262984314598913,
      -19.894334178865051,
      -19.705260025951397,
      -19.616414251366280,
      -19.594353920192376,
      -19.622080195525843
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf2b = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf2b ));
}

/* Reaction 3: e- + Ar(4S) -> Ar(4P) + e- */
static double _kf3(np_t np, gl_t *gl, double Te){
  double kf3;
  double Te_control[] = 
  { 
      8.641538368676400,
      8.666003207209252,
      8.728448288630180,
      8.848329763990707,
      8.882750345018193,
      8.969248049096372,
      9.200349244613665,
      9.306947205875320,
      9.360099936804787,
      9.835546044410583,
      10.470741420315500,
      11.105947357962146,
      11.741137983549356,
      12.376330992529057,
      13.011526632968218,
      13.646722101800350,
  };
  double kf_control[] = 
  {
      -55.121541432069527,
      -53.665373499501214,
      -50.046166698584322,
      -43.751677041594846,
      -42.155168197348900,
      -38.621012560539896,
      -32.111896182645033,
      -24.565678828026108,
      -19.703606572095222,
      -16.078749945616888,
      -14.966479864467910,
      -14.392116362296845,
      -14.053992391608563,
      -13.815797599156657,
      -13.638846642147865,
      -13.529670180391864
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf3 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf3 ));
}

/* Reaction 3b: e- + Ar(4P) -> Ar(4S) + e- */
static double _kf3b(np_t np, gl_t *gl, double Te){
  double kf3b;
  double Te_control[] = 
  { 
      7.056565294775152,
      7.144008245290196,
      7.231442544679870,
      7.318883404331368,
      7.756082145842774,
      8.193285498178808,
      8.630481723462170,
      9.067683829844169,
      9.504885757995108,
      9.942080835870144,
      10.379282081363659,
      10.816483193430843,
      11.253680020588133,
      11.690876720109380,
      12.128082432505435,
      12.565281340689033,
      13.002480545424085
  };
  double kf_control[] = 
  {
      -13.073206613774024,
      -13.096642025854697,
      -13.122268381916543,
      -13.150346272040810,
      -13.333689530040569,
      -13.566496551403407,
      -13.744195052547679,
      -13.880438488245078,
      -13.958950305180123,
      -13.991905312997696,
      -13.967449186182682,
      -13.896953829302708,
      -13.807006817473054,
      -13.721546738642060,
      -13.642582986812188,
      -13.570480278092608,
      -13.509107325143860
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf3b = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf3b ));
}

/* Reaction 4: e- + Ar -> Ar+ + e- + e-, ground ionization (15.8 eV) */
static double _kf4(np_t np, gl_t *gl, double Te){
  double kf4;
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
  kf4 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf4 ));
}

/* Reaction 5: e- + Ar(4S) -> Ar+ + e- + e-, stepwise ionization (4.2 eV) */
static double _kf5(np_t np, gl_t *gl, double Te){
  double kf5;
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
  kf5 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf5 ));
}

/* Reaction 6: e- + Ar(4P) -> Ar+ + e- + e-, stepwise ionization */
static double _kf6(np_t np, gl_t *gl, double Te){
  double kf6;
  double Te_control[] = 
  { 
      7.056565294775152,
      7.144008245290196,
      7.231442544679870,
      7.318883404331368,
      7.406327400007066,
      7.493768059470271,
      7.930962625044073,
      8.368163048991413,
      8.805363329735192,
      9.242563761130038,
      9.679761789441532,
      10.116961081654866,
      10.554163170723623,
      10.991362977169448,
      11.428561795411676,
      11.865761450084150,
      12.302962854928497,
      12.740162068771941,
      13.177361490415228,
      13.614560240133805
  };
  double kf_control[] = 
  {
      -59.678668592908387,
      -55.926524213203280,
      -52.489341015944966,
      -49.340568375320565,
      -46.455904265533476,
      -43.813077321745801,
      -33.569912188429463,
      -26.946243865712802,
      -22.653804731377736,
      -19.864383903196909,
      -18.045963544792983,
      -16.858166738229990,
      -16.084970400526490,
      -15.589964633035976,
      -15.285808460513575,
      -15.114988894507908,
      -15.039097807293787,
      -15.031701028272783,
      -15.073028312560323,
      -15.147665279383535
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf6 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf6 ));
}

// Reaction 7: e- + Ar2star -> Ar2+ + e- + e-
static double _kf7(np_t np, gl_t *gl, double Te){
  double kf7;
  kf7 = 9.00E-8 * pow( Te / TeV, 0.7) * exp( -3.66 * TeV / Te );
  return(kf7);
}

// Reaction 8: e- + Ar2star -> Ar + Ar + e-
static double _kf8(np_t np, gl_t *gl, double Te){
  double kf8;
  kf8 = 1.00E-9;
  return(kf8);
}

// Reaction 9: e- + Ar2+ -> Ar(4p) + Ar 
static double _kf9(np_t np, gl_t *gl, double Te){
  double kf9;
  kf9 = 5.40E-8 * pow( Te / TeV, -0.66);
  return(kf9);
}

// Reaction 10: e- + Ar+ -> Ar(4p)  
static double _kf10(np_t np, gl_t *gl, double Te){
  double kf10;
  kf10 = 4.00E-13 * pow( Te / TeV, -0.50);
  return(kf10);
}

// Reaction 11: e- + e- + Ar+ -> Ar(4p) + e- 
static double _kf11(np_t np, gl_t *gl, double Te){
  double kf11;
  kf11 = 5.00E-27 * pow( Te / TeV, -4.5); // cm6/s  
  return(kf11);
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


static double _averaged_rate_local_3r(np_t np, gl_t *gl, double Te, spec_t rhok,  long react1, long react2, long react3, char* id, double (*_kf)(np_t, gl_t *, double)){
  double kfave,rate;
  averagedrates_id_type id2;
  strcpy(id2,id);
  strins("Ne",id2,0);
  rate=_kf(np,gl,Te);
  kfave=min(_averaged_rate(np,gl,id,rate),_averaged_rate(np,gl,id2,rhok[react1]*rhok[react2]*rhok[react3]*rate)/(rhok[react1]*rhok[react2]*rhok[react3]));
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

  /* Ar electron impact reactions */    
  
  if (REACTION[1]){
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_W_2r2p ( speceminus, specAr, specAr4S, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), N, W );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4S, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf1b", &_kf1b), N, W );
     } else {
       add_to_W_2r2p ( speceminus, specAr, specAr4S, speceminus, _kf1(np,gl,Te), N, W);
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4S, specAr, speceminus, _kf1b(np,gl,Te), N, W);
     } 
  }

  if (REACTION[2]){
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_W_2r2p ( speceminus, specAr, specAr4P, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), N, W );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4P, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf2b", &_kf2b), N, W );
     } else {
       add_to_W_2r2p ( speceminus, specAr, specAr4P, speceminus, _kf2(np,gl,Te), N, W );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4P, specAr, speceminus, _kf2b(np,gl,Te), N, W);
     } 
  }

  if (REACTION[3]){
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_W_2r2p ( speceminus, specAr4S, specAr4P, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), N, W );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4P, specAr4S, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf3b", &_kf3b), N, W );
     } else {
       add_to_W_2r2p ( speceminus, specAr4S, specAr4P, speceminus, _kf3(np,gl,Te), N, W );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_W_2r2p ( speceminus, specAr4P, specAr4S, speceminus, _kf3b(np,gl,Te), N, W);
     } 
  }

  if (REACTION[4]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf4", &_kf4), N, W );
     } else {
       add_to_W_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _kf4(np,gl,Te), N, W );
     }
  }

  if (REACTION[5]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf5", &_kf5), N, W );
     } else {
       add_to_W_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _kf5(np,gl,Te), N, W );
     }
  }

  if (REACTION[6]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r3p ( speceminus, specAr4P, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf6", &_kf6), N, W );
     } else {
       add_to_W_2r3p ( speceminus, specAr4P, specArplus, speceminus, speceminus, _kf6(np,gl,Te), N, W );
     }
  }

  if (REACTION[7]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r3p ( speceminus, specAr2star, specAr2plus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf7", &_kf7), N, W );
     } else {
       add_to_W_fw_2r3p ( speceminus, specAr2star, specAr2plus, speceminus, speceminus, 9.0000E-8 * calA * pow (TeV, -0.7), 0.7, 3.66 * R * TeV, Te, X, W );
     }
  }

  if (REACTION[8]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r3p ( speceminus, specAr2star, specAr, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf8", &_kf8), N, W );
     } else {
       add_to_W_fw_2r3p ( speceminus, specAr2star, specAr, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, Te, X, W );
     }
  }

  if (REACTION[9]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r2p ( speceminus, specAr2plus, specAr4P, specAr, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2plus, "kf9", &_kf9), N, W );
     } else {
       add_to_W_fw_2r2p ( speceminus, specAr2plus, specAr4P, specAr, 5.4000E-8 * calA * pow (TeV, 0.66), -0.66, 0.0 * R, Te, X, W );
     }
  }

  if (REACTION[10]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_2r1p ( speceminus, specArplus, specAr4P, _averaged_rate_local(np, gl, Te, rhok, speceminus, specArplus, "kf10", &_kf10), N, W );
     } else {
       add_to_W_fw_2r1p ( speceminus, specArplus, specAr4P, 4.0000E-13 * calA * pow (TeV, 0.5), -0.5, 0.0 * R, Te, X, W );
     }
  }

  if (REACTION[11]){
     if (gl->model.chem.RATE_AVERAGING){
       add_to_W_3r2p ( speceminus, speceminus, specArplus, specAr4P, speceminus, _averaged_rate_local_3r(np, gl, Te, rhok, speceminus, speceminus, specArplus, "kf11", &_kf11), N, W);
     } else {
       add_to_W_fw_3r2p ( speceminus, speceminus, specArplus, specAr4P, speceminus, 5.0000E-27 * sqr(calA) * pow (TeV, 4.5), -4.5, 0.0 * R, Te, X, W );
     }
  }

  /* Penning ionization */  

  if (REACTION[12]){
     add_to_W_fw_2r3p ( specAr4S, specAr4S, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, W );
  }

  if (REACTION[13]){
     add_to_W_fw_2r3p ( specAr4S, specAr4P, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, W );
  }

  if (REACTION[14]){
     add_to_W_fw_2r3p ( specAr4P, specAr4P, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, W );
  }

  if (REACTION[15]){
     add_to_W_fw_2r4p ( specAr2star, specAr2star, specAr2plus, specAr, specAr, speceminus, 5.0000E-10 * calA, 0.0, 0.0 * R, T, X, W );
  }

  /* Three-body neutral reactions */  

  if (REACTION[16]){
     add_to_W_fw_3r2p ( specAr4S, specAr, specAr, specAr2star, specAr, 1.0000E-32 * sqr(calA), 0.0, 0.0 * R, T, X, W );
  }

  if (REACTION[17]){
     add_to_W_fw_3r2p ( specAr4P, specAr, specAr, specAr2star, specAr, 1.0000E-32 * sqr(calA), 0.0, 0.0 * R, T, X, W );
  }

  /* Charge transfer */

  if (REACTION[18]){
     add_to_W_fw_3r2p ( specArplus, specAr, specAr, specAr2plus, specAr, 2.5000E-31 * sqr(calA), 0.0, 0.0 * R, T, X, W );
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
       add_to_dW_2r2p ( speceminus, specAr, specAr4S, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4S, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf1b", &_kf1b), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       add_to_dW_2r2p ( speceminus, specAr, specAr4S, speceminus, _kf1(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4S, specAr, speceminus, _kf1b(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } 
  }

  if (REACTION[2]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r2p ( speceminus, specAr, specAr4P, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4P, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf2b", &_kf2b), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       add_to_dW_2r2p ( speceminus, specAr, specAr4P, speceminus, _kf2(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4P, specAr, speceminus, _kf2b(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } 
  }  

  if (REACTION[3]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r2p ( speceminus, specAr4S, specAr4P, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4P, specAr4S, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf3b", &_kf3b), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       add_to_dW_2r2p ( speceminus, specAr4S, specAr4P, speceminus, _kf3(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
       if (gl->model.chem.SUPERELASTIC_COLLISIONS) add_to_dW_2r2p ( speceminus, specAr4P, specAr4S, speceminus, _kf3b(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } 
  }

  if (REACTION[4]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf4", &_kf4), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       add_to_dW_2r3p ( speceminus, specAr, specArplus, speceminus, speceminus, _kf4(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     }
  }

  if (REACTION[5]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf5", &_kf5), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       add_to_dW_2r3p ( speceminus, specAr4S, specArplus, speceminus, speceminus, _kf5(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     }
  }

  if (REACTION[6]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r3p ( speceminus, specAr4P, specArplus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf6", &_kf6), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     } else {
       add_to_dW_2r3p ( speceminus, specAr4P, specArplus, speceminus, speceminus, _kf6(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe  );
     }
  }

  if (REACTION[7]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r3p ( speceminus, specAr2star, specAr2plus, speceminus, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf7", &_kf7), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
     } else {
       add_to_dW_fw_2r3p ( speceminus, specAr2star, specAr2plus, speceminus, speceminus, 9.0000E-8 * calA * pow (TeV, -0.7), 0.7, 3.66 * R * TeV, Te, X, dWdTe, dWdrhok );
     }
  }

  if (REACTION[8]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r3p ( speceminus, specAr2star, specAr, specAr, speceminus, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf8", &_kf8), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
     } else {
       add_to_dW_fw_2r3p ( speceminus, specAr2star, specAr, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, Te, X, dWdTe, dWdrhok );
     }
  }

  if (REACTION[9]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r2p ( speceminus, specAr2plus, specAr4P, specAr, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2plus, "kf9", &_kf9), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
     } else {
       add_to_dW_fw_2r2p ( speceminus, specAr2plus, specAr4P, specAr, 5.4000E-8 * calA * pow (TeV, 0.66), -0.66, 0.0 * R, Te, X, dWdTe, dWdrhok );
     }
  }

  if (REACTION[10]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_2r1p ( speceminus, specArplus, specAr4P, _averaged_rate_local(np, gl, Te, rhok, speceminus, specArplus, "kf10", &_kf10), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
     } else {
       add_to_dW_fw_2r1p ( speceminus, specArplus, specAr4P, 4.0000E-13 * calA * pow (TeV, 0.5), -0.5, 0.0 * R, Te, X, dWdTe, dWdrhok );
     }
  }

  if (REACTION[11]){
     dkfdTe = 0.0;
     dkfdT = 0.0;
     dkfdTv = 0.0;
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dW_3r2p ( speceminus, speceminus, specArplus, specAr4P, speceminus, _averaged_rate_local_3r(np, gl, Te, rhok, speceminus, speceminus, specArplus, "kf11", &_kf11), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
     } else {
      add_to_dW_fw_3r2p ( speceminus, speceminus, specArplus, specAr4P, speceminus, 5.0000E-27 * sqr(calA) * pow (TeV, 4.5), -4.5, 0.0 * R, Te, X, dWdTe, dWdrhok );
     }
  }

  /* Penning ionization */

  if (REACTION[12]){
     add_to_dW_fw_2r3p ( specAr4S, specAr4S, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

  if (REACTION[13]){
     add_to_dW_fw_2r3p ( specAr4S, specAr4P, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

  if (REACTION[14]){
     add_to_dW_fw_2r3p ( specAr4P, specAr4P, specArplus, specAr, speceminus, 1.0000E-9 * calA, 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

  if (REACTION[15]){
     add_to_dW_fw_2r4p ( specAr2star, specAr2star, specAr2plus, specAr, specAr, speceminus, 5.0000E-10 * calA, 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

  /* Three-body neutral reactions */

  if (REACTION[16]){
     add_to_dW_fw_3r2p ( specAr4S, specAr, specAr, specAr2star, specAr, 1.0000E-32 * sqr(calA), 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

  if (REACTION[17]){
     add_to_dW_fw_3r2p ( specAr4P, specAr, specAr, specAr2star, specAr, 1.0000E-32 * sqr(calA), 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }

  /* Charge transfer */

  if (REACTION[18]){
     add_to_dW_fw_3r2p ( specArplus, specAr, specAr, specAr2plus, specAr, 2.5000E-31 * sqr(calA), 0.0, 0.0 * R, T, X, dWdT, dWdrhok );
  }
  
}


void find_Qei_rodriguez2025(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
 
  *Qei=0.0;  

  if (gl->model.chem.TE_FROM_TOWNSEND) Te=_Te_from_rhok_EoverN(rhok, Estar);  

  if (REACTION[1]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr, 11.6, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf1", &_kf1), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr, 11.6, _kf1(np,gl,Te), rhok, Qei); // Ar->Ar(4s) (11.6 eV)
     }
  }

  if (REACTION[2]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr, 13.1, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr, 13.1, _kf2(np,gl,Te), rhok, Qei); // Ar->Ar(4p) (13.1 eV)
     }
  }

  if (REACTION[3]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr4S, 1.5, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr4S, 1.5, _kf3(np,gl,Te), rhok, Qei); // Ar(4s)->Ar(4p) (1.5 eV)
     }
  }

  if (REACTION[4]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr, _ionizationpot(specAr), _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf4", &_kf4), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr, _ionizationpot(specAr), _kf4(np,gl,Te), rhok, Qei); // ground ionization 15.8 eV
     } 
  }

  if (REACTION[5]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr4S, 4.2, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf5", &_kf5), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr4S, 4.2, _kf5(np,gl,Te), rhok, Qei); // Ar(4s)->Ar+ (15.8 - 11.6 = 4.2 eV)
     } 
  }

  if (REACTION[6]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr4P, 2.7, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf6", &_kf6), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr4P, 2.7, _kf6(np,gl,Te), rhok, Qei); // Ar(4p)->Ar+ (15.8 - 13.1 = 2.7 eV)
     } 
  }

  if (REACTION[7]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr2star, 4.2, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf7", &_kf7), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr2star, 4.2, 9.0000E-8 * pow (TeV, -0.7) * pow (Te, 0.7) * exp ( -3.66 * TeV / Te ), rhok, Qei); // (15.8 - 11.6 = 4.2 eV) estimate assuming Ar2* = Ar-Ar(4s)
     } 
  }

  if (REACTION[8]) {     
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_Qei(gl,Te, specAr2star, 4.2, _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf8", &_kf8), rhok, Qei);
     } else {
       add_to_Qei(gl,Te, specAr2star, 4.2, 1.0000E-9, rhok, Qei); // (15.8 - 11.6 = 4.2 eV) estimate assuming Ar2* = Ar-Ar(4s)
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
       add_to_dQei(gl,Te, specAr, 13.1,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf2", &_kf2), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr, 13.1, _kf2(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // Ar->Ar(4p) (13.1 eV)
     } 
  }

  if (REACTION[3]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr4S, 1.5,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf3", &_kf3), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr4S, 1.5, _kf3(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // Ar(4s)->Ar(4p) (1.5 eV)
     } 
  }

  if (REACTION[4]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr, _ionizationpot(specAr),  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr, "kf4", &_kf4), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr, _ionizationpot(specAr), _kf4(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // ground ionization 15.8 eV
     } 
  }

  if (REACTION[5]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr4S, 4.2,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4S, "kf5", &_kf5), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr4S, 4.2, _kf5(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // Ar(4s)->Ar+ (15.8 - 11.6 = 4.2 eV)
     } 
  }

  if (REACTION[6]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr4P, 2.7,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr4P, "kf6", &_kf6), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr4P, 2.7, _kf6(np,gl,Te), 0.0,  rhok, dQeidrhok, dQeidTe); // Ar(4p)->Ar+ (15.8 - 13.1 = 2.7 eV)
     } 
  }

  if (REACTION[7]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr2star, 4.2,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf7", &_kf7), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr2star, 4.2, 9.0000E-8 * pow (TeV, -0.7) * pow (Te, 0.7) * exp ( -3.66 * TeV / Te ), 0.0,  rhok, dQeidrhok, dQeidTe); // (15.8 - 11.6 = 4.2 eV) estimate assuming Ar2* = Ar-Ar(4s)
     } 
  }

  if (REACTION[8]) {
     if (gl->model.chem.RATE_AVERAGING) {
       add_to_dQei(gl,Te, specAr2star, 4.2,  _averaged_rate_local(np, gl, Te, rhok, speceminus, specAr2star, "kf8", &_kf8), 0.0,  rhok, dQeidrhok, dQeidTe);
     } else {
       add_to_dQei(gl,Te, specAr2star, 4.2, 1.0000E-9, 0.0,  rhok, dQeidrhok, dQeidTe); // (15.8 - 11.6 = 4.2 eV) estimate assuming Ar2* = Ar-Ar(4s)
     } 
  }

}




