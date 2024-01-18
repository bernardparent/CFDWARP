// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2024 Bernard Parent


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

#include <model/thermo/_thermo.h>
#include "thermo.h"
#include "electronimpact.h"



static double _EoverN_from_Te_NH3_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.10929975860069,
    6.11802895395280,
    6.12794891320948,
    6.14177343902385,
    6.15368424981004,
    6.16496717549049,
    6.18165640067889,
    6.21839234216407,
    6.31275619941210,
    7.17523682449177,
    9.47870962282596,
    10.3608846014599,
    10.6417495125494,
    10.8518297395922,
    11.0729483155267,
    11.3203712980607,
    11.6090728281791,
    11.9699552044317,
    12.4387641453030,
    13.0019859033808,
    13.4976713123232,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -50.8303924467099,
    -48.9790281007709,
    -48.5163354048294,
    -48.0534418938969,
    -47.5906813086665,
    -47.1278679597192,
    -46.6649294172650,
    -46.2021759732358,
    -45.7390833021391,
    -45.2765139689848,
    -44.8136175259262,
    -44.3507792573246,
    -43.8879191697673,
    -43.4248617142043,
    -42.9624792432396,
    -42.4995016046822,
    -42.0365808964492,
    -41.5737970249387,
    -41.1107739782095,
    -40.6481239826754,
    -40.2777716253226,
    -39.1439465808988
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}


static double _Te_from_EoverN_NH3_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -52.9594571388631,
    -52.4038493723252,
    -51.8485869527994,
    -51.2932608926103,
    -50.7378652050785,
    -50.1825029591136,
    -49.6268955495910,
    -49.0717268260039,
    -48.5163354048294,
    -47.9608944260011,
    -47.4053354387613,
    -46.8499873585642,
    -46.2946655386444,
    -45.7390833021391,
    -45.1837613889478,
    -44.6285935256383,
    -44.0730479720138,
    -43.5180050459235,
    -42.9624792432396,
    -42.4070355604433,
    -41.8514969069593,
    -41.2963890154631,
    -40.7404682680663,
    -40.2777716253226,
    -39.1439465808988
  };
  /* log K */
  double Te_data[] = 
  { 
    6.10904186046902,
    6.10904186046902,
    6.10904186046902,
    6.10904186046902,
    6.10929975860069,
    6.10981535541541,
    6.11161785365031,
    6.11649403423471,
    6.12794891320948,
    6.14401754917895,
    6.15786887130107,
    6.17346753025588,
    6.21073302654868,
    6.31275619941210,
    7.64657171641471,
    9.98615776582372,
    10.5488214073834,
    10.8098270078175,
    11.0729483155267,
    11.3739200660876,
    11.7423933837798,
    12.2377868894461,
    12.8837446475744,
    13.4976713123232,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}



static double _EoverN_from_Te_NH3_Snoeckx2023(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Snoeckx theoretical cross-sections */
  /* Snoeckx, Ramses, Jonathan Tennyson, and Min Suk Cha. 
    "Theoretical cross sections for electron collisions relevant for ammonia discharges part 1: NH3, NH2, and NH." 
     Plasma Sources Science and Technology 32.11 (2023): 115020. */  
  /* log K */
  double Te_data[] = 
  { 
    6.10955759023806,
    6.11161785365031,
    6.11546944398115,
    6.12592205667572,
    6.15565566493842,
    6.23490269197890,
    6.86782317017740,
    9.10206218040117,
    9.88328702559398,
    10.0900749326623,
    10.2211961047901,
    10.3520317218042,
    10.5011413997384,
    10.6771659081349,
    10.8936492135379,
    11.1641550837464,
    11.4921326963762,
    11.8952253241307,
    12.3737044155629,
    12.8951699952379,
    13.4588141114680,
    14.1023418716538,
    14.9385027840463,
    16.3712656820747
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -51.1111598186747,
    -49.9404967391899,
    -49.4479117000320,
    -48.9550372841689,
    -48.4620949055093,
    -47.9690245520843,
    -47.4763213353225,
    -46.9836139726580,
    -46.4906720465465,
    -45.9981610929529,
    -45.5047371894991,
    -45.0121321118466,
    -44.5191449917828,
    -44.0263205950510,
    -43.5331992606644,
    -43.0405884842886,
    -42.5476470927791,
    -42.0549701265258,
    -41.5620547651484,
    -41.0694660403064,
    -40.5764943265361,
    -40.0835059397123,
    -39.5905462323655,
    -39.1593648336272
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_NH3_Snoeckx2023(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Snoeckx theoretical cross-sections */
  /* Snoeckx, Ramses, Jonathan Tennyson, and Min Suk Cha. 
    "Theoretical cross sections for electron collisions relevant for ammonia discharges part 1: NH3, NH2, and NH." 
     Plasma Sources Science and Technology 32.11 (2023): 115020. */  
  /* log Vm2 */
  double EoverN_data[] = 
  { 
    -51.1111598186747,
    -49.9404967391899,
    -49.4479117000320,
    -48.9550372841689,
    -48.4620949055093,
    -47.9690245520843,
    -47.4763213353225,
    -46.9836139726580,
    -46.4906720465465,
    -45.9981610929529,
    -45.5047371894991,
    -45.0121321118466,
    -44.5191449917828,
    -44.0263205950510,
    -43.5331992606644,
    -43.0405884842886,
    -42.5476470927791,
    -42.0549701265258,
    -41.5620547651484,
    -41.0694660403064,
    -40.5764943265361,
    -40.0835059397123,
    -39.5905462323655,
    -39.1593648336272
  };
  /* log K */
  double Te_data[] = 
  { 
    6.10955759023806,
    6.11161785365031,
    6.11546944398115,
    6.12592205667572,
    6.15565566493842,
    6.23490269197890,
    6.86782317017740,
    9.10206218040117,
    9.88328702559398,
    10.0900749326623,
    10.2211961047901,
    10.3520317218042,
    10.5011413997384,
    10.6771659081349,
    10.8936492135379,
    11.1641550837464,
    11.4921326963762,
    11.8952253241307,
    12.3737044155629,
    12.8951699952379,
    13.4588141114680,
    14.1023418716538,
    14.9385027840463,
    16.3712656820747
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}


static double _EoverN_from_Te_N2_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.11777329755360,
    6.12744258408474,
    6.14600809474788,
    6.17901439121605,
    6.23308161248632,
    6.31422787002356,
    6.42639291337572,
    6.57026966875156,
    6.74253924698289,
    6.93496247301502,
    7.12888595635390,
    7.29794161434954,
    7.44590102095255,
    7.60352996515614,
    7.80088076190862,
    8.05462032845107,
    8.34937281280070,
    8.63306670696211,
    8.85827509485550,
    9.02311805437549,
    9.14756478528784,
    9.24362729651278,
    9.31530298798804,
    9.36711855741750,
    9.40603397366717,
    9.44165160928006,
    9.51272947569662,
    9.80575743920766,
    10.3855503423225,
    10.8663340661637,
    11.1894893389638,
    11.3959453547267,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.8453075314907,
    -54.4291331089220,
    -54.0124270696525,
    -53.5958459856044,
    -53.1793559600124,
    -52.7630683248577,
    -52.3464360027824,
    -51.9294806425851,
    -51.5130088549234,
    -51.0966180879540,
    -50.6800383239008,
    -50.2634795189951,
    -49.8468311138376,
    -49.4304530527132,
    -49.0138059524549,
    -48.5972506316385,
    -48.1803336457515,
    -47.7642805312345,
    -47.3477899531147,
    -46.9311786186324,
    -46.5145312853787,
    -46.0979552443984,
    -45.6815185659174,
    -45.2650643362336,
    -44.8484293006694,
    -44.4317197335711,
    -44.0152983078927,
    -43.5989741084572,
    -43.1820998634711,
    -42.7655412882071,
    -42.4532943089311,
    -38.0453342922307
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_N2_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.5330453984035,
    -53.8041251270642,
    -53.0752047474074,
    -52.3464360027824,
    -51.6173759324195,
    -50.8882998167709,
    -50.1591316616517,
    -49.4304530527132,
    -48.7014359841603,
    -47.9724317104059,
    -47.2434167668114,
    -46.5145312853787,
    -45.7854988191063,
    -45.0565999475825,
    -44.3276860427394,
    -43.5989741084572,
    -42.8699050432541,
    -42.1406793547864,
    -41.3076396750262,
    -40.6829260296843,
    -39.1439465808988
  };
  /* log K */
  double Te_data[] = 
  { 
    6.11777329755360,
    6.14027456290012,
    6.20277976825152,
    6.33933031868524,
    6.57026966875156,
    6.88577697379399,
    7.21698201618018,
    7.48313936609267,
    7.80088076190862,
    8.27444100426920,
    8.75401408453109,
    9.05736939093704,
    9.24362729651278,
    9.35544352583699,
    9.42306371351197,
    9.55307108040562,
    10.3855503423225,
    11.1173530734675,
    11.6034705726304,
    12.2293194383470,
    12.7963582069535,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_O2_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.15491683963699,
    6.18954080420304,
    6.24147644340909,
    6.31086086347565,
    6.39511056094682,
    6.49251468754874,
    6.61074683231689,
    6.75937146194990,
    6.92929994372840,
    7.09866647079291,
    7.24304486171689,
    7.34474724226770,
    7.41003510776406,
    7.45802238148490,
    7.51018544476948,
    7.58542728102636,
    7.69155356116069,
    7.82621044630694,
    7.99383443303362,
    8.25592748654253,
    8.76185869285190,
    9.38188987473781,
    9.84773040258699,
    10.1380164437308,
    10.3249928984717,
    10.4680429858440,
    10.5960742922339,
    10.7310922426638,
    10.8945110968061,
    11.1022446850151,
    11.3535785551212,
    11.5584838509519,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  {
  -55.2620422318571,
  -54.8453075314907,
  -54.4291331089220,
  -54.0124270696525,
  -53.5958459856044,
  -53.1793559600124,
  -52.7630683248577,
  -52.3464360027824,
  -51.9294806425851,
  -51.5130088549234,
  -51.0966180879540,
  -50.6800383239008,
  -50.2634795189951,
  -49.8468311138376,
  -49.4304530527132,
  -49.0138059524549,
  -48.5972506316385,
  -48.1803336457515,
  -47.7642805312345,
  -47.3477899531147,
  -46.9311786186324,
  -46.5145312853787,
  -46.0979552443984,
  -45.6815185659174,
  -45.2650643362336,
  -44.8484293006694,
  -44.4317197335711,
  -44.0152983078927,
  -43.5989741084572,
  -43.1820998634711,
  -42.7655412882071,
  -42.4532943089311,
  -38.0453342922307
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_O2_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.8453075314907,
    -54.4291331089220,
    -54.0124270696525,
    -53.5958459856044,
    -53.1793559600124,
    -52.7630683248577,
    -52.3464360027824,
    -51.9294806425851,
    -51.5130088549234,
    -51.0966180879540,
    -50.6800383239008,
    -50.2634795189951,
    -49.8468311138376,
    -49.4304530527132,
    -49.0138059524549,
    -48.5972506316385,
    -48.1803336457515,
    -47.7642805312345,
    -47.3477899531147,
    -46.9311786186324,
    -46.5145312853787,
    -46.0979552443984,
    -45.6815185659174,
    -45.2650643362336,
    -44.8484293006694,
    -44.4317197335711,
    -44.0152983078927,
    -43.5989741084572,
    -43.1820998634711,
    -42.7655412882071,
    -42.4532943089311,
    -39.1439465808988
  };
  /* log K */
  double Te_data[] = 
  { 
    6.15491683963699,
    6.18954080420304,
    6.24147644340909,
    6.31086086347565,
    6.39511056094682,
    6.49251468754874,
    6.61074683231689,
    6.75937146194990,
    6.92929994372840,
    7.09866647079291,
    7.24304486171689,
    7.34474724226770,
    7.41003510776406,
    7.45802238148490,
    7.51018544476948,
    7.58542728102636,
    7.69155356116069,
    7.82621044630694,
    7.99383443303362,
    8.25592748654253,
    8.76185869285190,
    9.38188987473781,
    9.84773040258699,
    10.1380164437308,
    10.3249928984717,
    10.4680429858440,
    10.5960742922339,
    10.7310922426638,
    10.8945110968061,
    11.1022446850151,
    11.3535785551212,
    11.5584838509519,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}


double _EoverN_from_rhok_Te(spec_t rhok, double Te){
  double Estar,Nsum;
  long spec;
  Estar=0.0;
  Nsum=0.0;
  for (spec=0; spec<ns; spec++){
    switch (smap[spec]){
      case SMAP_NH3:
        Estar+=rhok[spec]/_m(spec)*_EoverN_from_Te_NH3_Morgan(Te);
        Nsum+=rhok[spec]/_m(spec);
      case SMAP_N2:
        Estar+=rhok[spec]/_m(spec)*_EoverN_from_Te_N2_Morgan(Te);
        Nsum+=rhok[spec]/_m(spec);
      case SMAP_O2:
        Estar+=rhok[spec]/_m(spec)*_EoverN_from_Te_O2_Morgan(Te);
        Nsum+=rhok[spec]/_m(spec);
      break;
    }
  }
  if (Nsum==0.0) fatal_error("Couldn't find any valid species when determining Estar in _EoverN_from_Nk_Te().");
  Estar/=Nsum;
  return(Estar);
}


double _Te_from_rhok_EoverN(spec_t rhok, double EoverN){
  double Te,Nsum;
  long spec;
  Te=0.0;
  Nsum=0.0;
  for (spec=0; spec<ns; spec++){
    switch (smap[spec]){
      case SMAP_NH3:
        Te+=rhok[spec]/_m(spec)*_Te_from_EoverN_NH3_Morgan(EoverN);
        Nsum+=rhok[spec]/_m(spec);
      case SMAP_N2:
        Te+=rhok[spec]/_m(spec)*_Te_from_EoverN_N2_Morgan(EoverN);
        Nsum+=rhok[spec]/_m(spec);
      case SMAP_O2:
        Te+=rhok[spec]/_m(spec)*_Te_from_EoverN_O2_Morgan(EoverN);
        Nsum+=rhok[spec]/_m(spec);
      break;
    }
  }
  if (Nsum==0.0) fatal_error("Couldn't find any valid species when determining Te in _Te_from_Nk_EoverN().");
  Te/=Nsum;
  return(Te);
}

