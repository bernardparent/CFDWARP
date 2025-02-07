// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2024 Bernard Parent
Copyright 2024,2025 Felipe Martin Rodriguez Fuentes


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
#include <model/transport/_transport.h>
#include "thermo.h"
#include "electronimpact.h"


static double _EoverN_from_Te_H_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections
    Two excitation cross-section sets
    Ionization cross section in the Morgan database is sourced from
    Fite, Wade L., and R. T. Brackmann. "Collisions of electrons with hydrogen atoms. I. Ionization." 
    Physical Review 112.4 (1958): 1141. 
  */
  /* log K */
  double Te_data[] = 
  { 
    6.11033068652705,
    6.11393058756052,
    6.14077443794177,
    6.23216982760437,
    6.68882940347612,
    7.85542292975633,
    8.88272619071966,
    9.86053155270631,
    10.5666170813485,
    11.1170083053207,
    13.0941975227273,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -52.1379167435347,
    -51.3570438413659,
    -50.3807566097887,
    -49.5997764374038,
    -48.4282564730751,
    -47.0614794348485,
    -46.2804862910852,
    -45.6947269609332,
    -44.9135483965923,
    -43.5469925827967,
    -41.5943358686612,
    -40.4532799008825
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_H_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections
     log Vm^2*/
  double EoverN_data[] = 
  { 
    -52.1379167435347,
    -51.3570438413659,
    -50.3807566097887,
    -49.5997764374038,
    -48.4282564730751,
    -47.0614794348485,
    -46.2804862910852,
    -45.6947269609332,
    -44.9135483965923,
    -43.5469925827967,
    -41.5943358686612,
    -40.4532799008825
  };
  /* log K */
  double Te_data[] = 
  { 
    6.11033068652705,
    6.11393058756052,
    6.14077443794177,
    6.23216982760437,
    6.68882940347612,
    7.85542292975633,
    8.88272619071966,
    9.86053155270631,
    10.5666170813485,
    11.1170083053207,
    13.0941975227273,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_H2_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections 
    Added the Phelps' data for the 13.86 eV excitation process to the Morgan database. 
    Morgan's original data didn't include the 13.86 eV excitation process.
  */

  /* log K */
  double Te_data[] = 
  { 
    6.11033068652705,
    6.11393058756052,
    6.13676842357083,
    6.20442238272387,
    6.51132074413782,
    7.37429148857443,
    8.37645272369656,
    8.84949075064230,
    9.34484860305746,
    9.93632539207584,
    10.7526685939775,
    11.4948541084563,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -53.3095821519381,
    -52.5283242840063,
    -51.5521354105511,
    -50.7712733061790,
    -49.5997764374038,
    -48.2329546677074,
    -46.8661131279339,
    -46.0851552361401,
    -45.3040664940369,
    -44.7182809883681,
    -43.9373759163476,
    -42.7659153297301,
    -40.0355447001826
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_H2_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections
     log Vm^2*/
  double EoverN_data[] = 
  { 
    -53.3095821519381,
    -52.5283242840063,
    -51.5521354105511,
    -50.7712733061790,
    -49.5997764374038,
    -48.2329546677074,
    -46.8661131279339,
    -46.0851552361401,
    -45.3040664940369,
    -44.7182809883681,
    -43.9373759163476,
    -42.7659153297301,
    -40.0355447001826
  };
  /* log K */
  double Te_data[] = 
  { 
    6.11033068652705,
    6.11393058756052,
    6.13676842357083,
    6.20442238272387,
    6.51132074413782,
    7.37429148857443,
    8.37645272369656,
    8.84949075064230,
    9.34484860305746,
    9.93632539207584,
    10.7526685939775,
    11.4948541084563,
    15.4249484703984L
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}


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



  /* Data in log-log coordinates. Obtained with BOLSIG+ using Snoeckx theoretical cross-sections */
  /* Snoeckx, Ramses, Jonathan Tennyson, and Min Suk Cha. 
    "Theoretical cross sections for electron collisions relevant for ammonia discharges part 1: NH3, NH2, and NH." 
     Plasma Sources Science and Technology 32.11 (2023): 115020. */  
  /* log K */
/*static double _EoverN_from_Te_NH3_Snoeckx2023(double Te){
  double EoverN;
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
  // log Vm^2
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
*/



  /* Data in log-log coordinates. Obtained with BOLSIG+ using Snoeckx theoretical cross-sections */
  /* Snoeckx, Ramses, Jonathan Tennyson, and Min Suk Cha. 
    "Theoretical cross sections for electron collisions relevant for ammonia discharges part 1: NH3, NH2, and NH." 
     Plasma Sources Science and Technology 32.11 (2023): 115020. */  
  /* log Vm2 */
/*static double _Te_from_EoverN_NH3_Snoeckx2023(double EoverN){
  double Te;
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
  // log K 
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
*/

static double _EoverN_from_Te_NH3v_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.13324994796315,
    6.21026693472015,
    6.40551929023296,
    6.46434883279277,
    6.54022449032017,
    6.61246330305119,
    7.03247749924518,
    9.36513245944587,
    10.5660189952063,
    10.9749679071665,
    12.5993962383727,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -49.7948263288472,
    -49.0139993573487,
    -48.0380174235713,
    -47.2566751645404,
    -46.4757383535745,
    -46.0851552361401,
    -45.3040664940369,
    -44.7182809883681,
    -43.9373759163476,
    -43.1563425604737,
    -41.2035854952824,
    -39.4586573257385
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_NH3v_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -49.7948263288472,
    -49.0139993573487,
    -48.0380174235713,
    -47.2566751645404,
    -46.4757383535745,
    -46.0851552361401,
    -45.3040664940369,
    -44.7182809883681,
    -43.9373759163476,
    -43.1563425604737,
    -41.2035854952824,
    -39.4586573257385
  };
  /* log K */
  double Te_data[] = 
  { 
    6.13324994796315,
    6.21026693472015,
    6.40551929023296,
    6.46434883279277,
    6.54022449032017,
    6.61246330305119,
    7.03247749924518,
    9.36513245944587,
    10.5660189952063,
    10.9749679071665,
    12.5993962383727,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}



static double _EoverN_from_Te_N2(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    5.78359961896139,
    6.00674317027560,
    7.05656529477427,
    9.05804529498440,
    9.69562262438953,
    10.0010042739407,
    11.4313155204774,
    11.9302347337974,
    14.1582412663937,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -51.8608448501949,
    -51.0135469898077,
    -49.5582597572009,
    -46.7448490404409,
    -44.4422639474468,
    -43.7491167668869,
    -42.3992270876608,
    -41.6847888630171,
    -39.5409971904059,
    -38.7384814727906
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_N2(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -51.8608448501949,
    -51.0135469898077,
    -49.5582597572009,
    -46.7448490404409,
    -44.4422639474468,
    -43.7491167668869,
    -42.3992270876608,
    -41.6847888630171,
    -39.5409971904059,
    -38.7384814727906
  };
  /* log K */
  double Te_data[] = 
  { 
    5.78359961896139,
    6.00674317027560,
    7.05656529477427,
    9.05804529498440,
    9.69562262438953,
    10.0010042739407,
    11.4313155204774,
    11.9302347337974,
    14.1582412663937,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_O2(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.15485809401642,
    7.39303753139549,
    8.64580049989085,
    9.45446056757264,
    10.2346191251222,
    10.9062872740050,
    11.1835383764881,
    11.8792632932909,
    12.6549872537726,
    13.6691438195424,
    14.3029333748767,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  {
    -50.6568720458690,
    -48.3542869528750,
    -46.7448490404409,
    -46.0517018598809,
    -44.9530895712128,
    -43.5718077519165,
    -43.0405884842886,
    -41.9779001988901,
    -40.9153153604791,
    -39.8524104655706,
    -39.3210418187212,
    -38.7384814727906
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_O2(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -50.6568720458690,
    -48.3542869528750,
    -46.7448490404409,
    -46.0517018598809,
    -44.9530895712128,
    -43.5718077519165,
    -43.0405884842886,
    -41.9779001988901,
    -40.9153153604791,
    -39.8524104655706,
    -39.3210418187212,
    -38.7384814727906
  };
  /* log K */
  double Te_data[] = 
  { 
    6.15485809401642,
    7.39303753139549,
    8.64580049989085,
    9.45446056757264,
    10.2346191251222,
    10.9062872740050,
    11.1835383764881,
    11.8792632932909,
    12.6549872537726,
    13.6691438195424,
    14.3029333748767,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_NO(double Te){
  double EoverN;
  /* Data in log-log coordinates */
  /* log K */
  double Te_data[] = 
  { 
    4.60517018598809,
    5.03874736289747,
    5.73837712326189,
    6.04379878382206,
    6.90775527898214,
    8.11634666097751,
    8.36842970883516,
    9.40042667908016,
    9.92446419681838,
    10.2589051510024,
    11.3016260743750,
    12.0194099250342,
    13.1593475780142,
    14.9141228466324
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  {
    -53.3604828483550,
    -52.4569498135940,
    -51.1678011430320,
    -50.6571555835696,
    -49.3445762009764,
    -45.9733528934699,
    -45.3366953449681,
    -44.6078652696987,
    -44.2445612886146,
    -43.7195579646453,
    -42.4083429210688,
    -41.3594369670419,
    -40.0475679096723,
    -38.7384814727906
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_NO(double EoverN){
  double Te;
  /* Data in log-log coordinates */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -53.3604828483550,
    -52.4569498135940,
    -51.1678011430320,
    -50.6571555835696,
    -49.3445762009764,
    -45.9733528934699,
    -45.3366953449681,
    -44.6078652696987,
    -44.2445612886146,
    -43.7195579646453,
    -42.4083429210688,
    -41.3594369670419,
    -40.0475679096723,
    -38.7384814727906
  };
  /* log K */
  double Te_data[] = 
  { 
    4.60517018598809,
    5.03874736289747,
    5.73837712326189,
    6.04379878382206,
    6.90775527898214,
    8.11634666097751,
    8.36842970883516,
    9.40042667908016,
    9.92446419681838,
    10.2589051510024,
    11.3016260743750,
    12.0194099250342,
    13.1593475780142,
    14.9141228466324
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_N_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.13224239120285,
    6.21236263929437,
    6.44739154799051,
    7.58719354583644,
    8.72559199862003,
    9.27739788589606,
    9.47782191748582,
    9.59343168349298,
    10.0915182814815,
    11.1390059524936,
    11.7497463580851,
    13.5255054859667,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.4161739642795,
    -53.5702872769338,
    -51.4558242853259,
    -49.7640546701805,
    -48.9181617977308,
    -48.0721200611113,
    -47.2264395389385,
    -45.1116946013894,
    -43.4205327031147,
    -42.5744697919428,
    -40.8823548746298,
    -39.6300795920744
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}


static double _Te_from_EoverN_N_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates, Morgan cross-section database*/
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.4161739642795,
    -53.5702872769338,
    -51.4558242853259,
    -49.7640546701805,
    -48.9181617977308,
    -48.0721200611113,
    -47.2264395389385,
    -45.1116946013894,
    -43.4205327031147,
    -42.5744697919428,
    -40.8823548746298,
    -39.6300795920744
  };
  /* log K */
  double Te_data[] = 
  { 
    6.13224239120285,
    6.21236263929437,
    6.44739154799051,
    7.58719354583644,
    8.72559199862003,
    9.27739788589606,
    9.47782191748582,
    9.59343168349298,
    10.0915182814815,
    11.1390059524936,
    11.7497463580851,
    13.5255054859667,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_O_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.13952428150941,
    6.23921451075502,
    6.50509127253047,
    7.63149184394679,
    8.67714095161883,
    9.13625680520932,
    9.26681577656552,
    9.38481813451690,
    9.57989105446622,
    10.1629440884153,
    10.8830304118408,
    11.2996171771576,
    11.6030464844137,
    12.4065258944741,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.4161739642795,
    -53.5702872769338,
    -51.4558242853259,
    -49.7640546701805,
    -48.9181617977308,
    -48.0721200611113,
    -47.2264395389385,
    -46.3806226803493,
    -45.1116946013894,
    -44.2659604629784,
    -43.4205327031147,
    -42.5744697919428,
    -40.8823548746298,
    -39.0041846385236
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_O_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates, Morgan cross-section database*/
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.4161739642795,
    -53.5702872769338,
    -51.4558242853259,
    -49.7640546701805,
    -48.9181617977308,
    -48.0721200611113,
    -47.2264395389385,
    -46.3806226803493,
    -45.1116946013894,
    -44.2659604629784,
    -43.4205327031147,
    -42.5744697919428,
    -40.8823548746298,
    -39.0041846385236
  };
  /* log K */
  double Te_data[] = 
  { 
    6.13952428150941,
    6.23921451075502,
    6.50509127253047,
    7.63149184394679,
    8.67714095161883,
    9.13625680520932,
    9.26681577656552,
    9.38481813451690,
    9.57989105446622,
    10.1629440884153,
    10.8830304118408,
    11.2996171771576,
    11.6030464844137,
    12.4065258944741,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_H2O_Triniti(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.11007305416700,
    6.11598187033107,
    6.15170894052578,
    6.42036539465608,
    6.77922951437664,
    7.53899144401857,
    9.76394860688878,
    10.9426549181759,
    11.3090525580932,
    12.7398048681991,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -49.3036175388273,
    -48.2400658087849,
    -47.1756319565333,
    -45.6863645428635,
    -45.0477666181547,
    -44.6221085743659,
    -44.1963418092433,
    -43.5576703023159,
    -42.9196919680344,
    -40.7916057061531,
    -38.8685901581376
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_H2O_Triniti(double EoverN){
  double Te;
  /* Data in log-log coordinates, Triniti cross-section database*/
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -49.3036175388273,
    -48.2400658087849,
    -47.1756319565333,
    -45.6863645428635,
    -45.0477666181547,
    -44.6221085743659,
    -44.1963418092433,
    -43.5576703023159,
    -42.9196919680344,
    -40.7916057061531,
    -38.8685901581376
  };
  /* log K */
  double Te_data[] = 
  { 
    6.11007305416700,
    6.11598187033107,
    6.15170894052578,
    6.42036539465608,
    6.77922951437664,
    7.53899144401857,
    9.76394860688878,
    10.9426549181759,
    11.3090525580932,
    12.7398048681991,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_C2H4_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.10981535541541,
    6.11213225695384,
    6.11546944398115,
    6.11879553122299,
    6.14873864131225,
    6.27655218119267,
    6.49022696107802,
    7.51145528620996,
    9.93238727165571,
    11.1644839770235,
    12.6660370899592,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.4110369666816,
    -53.9852837578055,
    -53.7725146160956,
    -52.9211984267460,
    -50.7932646707746,
    -48.6651694184972,
    -46.5371846759751,
    -44.4092160934006,
    -42.2812424187746,
    -40.1531753084322,
    -38.4507994003388
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_C2H4_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates, Morgan cross-section database*/
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -54.4110369666816,
    -53.9852837578055,
    -53.7725146160956,
    -52.9211984267460,
    -50.7932646707746,
    -48.6651694184972,
    -46.5371846759751,
    -44.4092160934006,
    -42.2812424187746,
    -40.1531753084322,
    -38.4507994003388
  };
  /* log K */
  double Te_data[] = 
  { 
    6.10981535541541,
    6.11213225695384,
    6.11546944398115,
    6.11879553122299,
    6.14873864131225,
    6.27655218119267,
    6.49022696107802,
    7.51145528620996,
    9.93238727165571,
    11.1644839770235,
    12.6660370899592,
    15.4249484703984
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

static double _EoverN_from_Te_CO_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.11136055273548,
    6.11238935941058,
    6.11393058756052,
    6.12032692509042,
    6.12617563852062,
    6.19120529156866,
    6.47549601761968,
    6.90457893777639,
    7.14823248330010,
    7.26601551895648,
    7.49546517896868,
    8.04126205931249,
    8.71326839967549,
    9.05180967765335,
    9.15089544894786,
    9.42774317923393,
    10.0823416896766,
    11.0675282480573,
    12.9578315739641,
    17.1609507896773
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -55.0493531384468,
    -54.8361211151816,
    -54.4110369666816,
    -54.1980213879586,
    -53.3468858844851,
    -52.2829475994561,
    -51.2188155408136,
    -50.1548853707702,
    -49.0907592581733,
    -48.0271438115423,
    -46.9627563249598,
    -45.8989807728633,
    -44.8347141184775,
    -44.4092160934006,
    -43.9835740781014,
    -43.5576703023159,
    -42.7067221582170,
    -40.5790111063881,
    -38.4507994003388
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_CO_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates, Morgan cross-section database*/
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -55.2620422318571,
    -55.0493531384468,
    -54.8361211151816,
    -54.4110369666816,
    -54.1980213879586,
    -53.3468858844851,
    -52.2829475994561,
    -51.2188155408136,
    -50.1548853707702,
    -49.0907592581733,
    -48.0271438115423,
    -46.9627563249598,
    -45.8989807728633,
    -44.8347141184775,
    -44.4092160934006,
    -43.9835740781014,
    -43.5576703023159,
    -42.7067221582170,
    -40.5790111063881,
    -38.4507994003388
  };
  /* log K */
  double Te_data[] = 
  { 
    6.11136055273548,
    6.11238935941058,
    6.11393058756052,
    6.12032692509042,
    6.12617563852062,
    6.19120529156866,
    6.47549601761968,
    6.90457893777639,
    7.14823248330010,
    7.26601551895648,
    7.49546517896868,
    8.04126205931249,
    8.71326839967549,
    9.05180967765335,
    9.15089544894786,
    9.42774317923393,
    10.0823416896766,
    11.0675282480573,
    12.9578315739641,
    17.1609507896773
  };

  int N = sizeof(EoverN_data)/sizeof(EoverN_data[0]);
  EoverN = ( min( EoverN_data[N-1], max( log( EoverN) , EoverN_data[0] ) ) );
  Te = exp( EXM_f_from_monotonespline(N, EoverN_data, Te_data, EoverN) );
  return(Te);
}

double _EoverNk_from_Te(long spec, double Te){
  double Estar;
  switch (smap[spec]){
    case SMAP_NH3:
      Estar=_EoverN_from_Te_NH3_Morgan(Te);
    break;
    case SMAP_NH3v:
      Estar=_EoverN_from_Te_NH3v_Morgan(Te);
    break;
    case SMAP_N2:
      Estar=_EoverN_from_Te_N2(Te);
    break;
    case SMAP_O2:
      Estar=_EoverN_from_Te_O2(Te);
    break;
    case SMAP_NO:
      Estar=_EoverN_from_Te_NO(Te);
    break;
    case SMAP_H:
      Estar=_EoverN_from_Te_H_Morgan(Te);
    break;
    case SMAP_H2:
      Estar=_EoverN_from_Te_H2_Morgan(Te);
    break;
    case SMAP_N:
      Estar=_EoverN_from_Te_N_Morgan(Te);
    break;
    case SMAP_O:
      Estar=_EoverN_from_Te_O_Morgan(Te);
    break;
    case SMAP_H2O:
      Estar=_EoverN_from_Te_H2O_Triniti(Te);
    break;
    case SMAP_C2H4:
      Estar=_EoverN_from_Te_C2H4_Morgan(Te);
    break;
    case SMAP_CO:
      Estar=_EoverN_from_Te_CO_Morgan(Te);
    break;
    default:
      Estar=0.0;
  }
  return(Estar);
}


double _Tek_from_EoverN(long spec, double EoverN){
  double Te;
  switch (smap[spec]){
    case SMAP_NH3:
      Te=_Te_from_EoverN_NH3_Morgan(EoverN);
    break;
    case SMAP_NH3v:
      Te=_Te_from_EoverN_NH3v_Morgan(EoverN);
    break;
    case SMAP_N2:
      Te=_Te_from_EoverN_N2(EoverN);
    break;
    case SMAP_O2:
      Te=_Te_from_EoverN_O2(EoverN);
    break;
    case SMAP_NO:
      Te=_Te_from_EoverN_NO(EoverN);
    break;
    case SMAP_H:
      Te=_Te_from_EoverN_H_Morgan(EoverN);
    break;
    case SMAP_H2:
      Te=_Te_from_EoverN_H2_Morgan(EoverN);
    break;
    case SMAP_N:
      Te=_Te_from_EoverN_N_Morgan(EoverN);
    break;
    case SMAP_O:
      Te=_Te_from_EoverN_O_Morgan(EoverN);
    break;
    case SMAP_H2O:
      Te=_Te_from_EoverN_H2O_Triniti(EoverN);
    break;
    case SMAP_C2H4:
      Te=_Te_from_EoverN_C2H4_Morgan(EoverN);
    break;
    case SMAP_CO:
      Te=_Te_from_EoverN_CO_Morgan(EoverN);
    break;
    default:
      Te=0.0;
  }
  return(Te);
}


double _EoverN_from_rhok_Te(spec_t rhok, double Te){
  double Estar,Nsum,EoverNk;
  long spec;
  Estar=0.0;
  Nsum=0.0;
  for (spec=0; spec<ns; spec++){
    EoverNk=_EoverNk_from_Te(spec, Te);
    if (EoverNk>1e-30){
      Estar+=rhok[spec]/_m(spec)*EoverNk;
      Nsum+=rhok[spec]/_m(spec);
    }
  }
  if (Nsum==0.0) fatal_error("Couldn't find any valid species when determining Estar in _EoverN_from_rhok_Te().");
  Estar/=Nsum;
  return(Estar);
}


double _Te_from_rhok_EoverN(spec_t rhok, double EoverN){
  double Te,Nsum,Tek;
  long spec;
  Te=0.0;
  Nsum=0.0;
  for (spec=0; spec<ns; spec++){
    Tek=_Tek_from_EoverN(spec, EoverN);
    if (Tek>1e-30){
      Te+=rhok[spec]/_m(spec)*Tek;
      Nsum+=rhok[spec]/_m(spec);
    }
  }
  if (Nsum==0.0) fatal_error("Couldn't find any valid species when determining Te in _Te_from_rhok_EoverN().");
  Te/=Nsum;
  return(Te);
}

/* find the product of inelastic forward (k_kl) and backward (kb_lk) rate coefficients and the activation energy
 * for the kth N2vib species process
 * sum_l (k_kl*e_kl-kb_lk*e_lk) (cooling-heating from N2vib)
 * k_kl, kb_lk in m^3/s 
 * e_kl, e_lk in Joule */
static double _ke_kl_N2vib_from_Te(long k, double Te) {
  double ke_kl, Te_data[14], ke_kl_data[14]; 
  switch (k){
    case 1: 
      /* Te */
      Te_data[0] =  618.908000328523;
      Te_data[1] =  2262.88237620116;
      Te_data[2] =  3371.50133178963;
      Te_data[3] =  4438.34399735592;
      Te_data[4] =  6282.68983833492;
      Te_data[5] =  8943.22060474715;
      Te_data[6] =  19696.7471104552;
      Te_data[7] =  28887.5309153338;
      Te_data[8] =  37250.5252697730;
      Te_data[9] =  45319.5383240561;
      Te_data[10] = 62679.9077332711;
      Te_data[11] = 156506.360583075;
      Te_data[12] = 373510.978198263;
      Te_data[13] = 1059879.95056260;
      /* sum_l(k_kl*e_kl) */
      ke_kl_data[0] =   8.91699416535870e-42;
      ke_kl_data[1] =  -1.70432917313655e-36;
      ke_kl_data[2] =   1.37357807186088e-36;
      ke_kl_data[3] =   1.43680797272169e-35;
      ke_kl_data[4] =   6.71910422618799e-35;
      ke_kl_data[5] =   1.72383551343784e-34;
      ke_kl_data[6] =   3.49340836324648e-34;
      ke_kl_data[7] =   3.53621371637705e-34;
      ke_kl_data[8] =   3.21752476210811e-34;
      ke_kl_data[9] =   2.81368172191871e-34;
      ke_kl_data[10] =  2.00034236128683e-34;
      ke_kl_data[11] =  1.57846121546353e-35;
      ke_kl_data[12] = -5.40441736086305e-35;
      ke_kl_data[13] = -6.84031145203270e-35;
    break; 
    case 2: 
      Te_data[0] =  618.908000328523;
      Te_data[1] =  2262.88237620116;
      Te_data[2] =  3371.50133178963;
      Te_data[3] =  4438.34399735592;
      Te_data[4] =  6282.68983833492;
      Te_data[5] =  8943.22060474715;
      Te_data[6] =  19696.7471104552;
      Te_data[7] =  28887.5309153338;
      Te_data[8] =  37250.5252697730;
      Te_data[9] =  45319.5383240561;
      Te_data[10] = 62679.9077332711;
      Te_data[11] = 156506.360583075;
      Te_data[12] = 373510.978198263;
      Te_data[13] = 1059879.95056260;

      ke_kl_data[0] =  -1.01257563268800e-43;
      ke_kl_data[1] =  -1.29299691848918e-35;
      ke_kl_data[2] =  -7.70935673183447e-35;
      ke_kl_data[3] =  -1.80171427548097e-34;
      ke_kl_data[4] =  -3.62237717357694e-34;
      ke_kl_data[5] =  -4.92662265377810e-34;
      ke_kl_data[6] =  -4.06897910377454e-34;
      ke_kl_data[7] =  -2.82188326410815e-34;
      ke_kl_data[8] =  -1.95617115257576e-34;
      ke_kl_data[9] =  -1.36113556812124e-34;
      ke_kl_data[10] = -6.43003150699854e-35;
      ke_kl_data[11] = -6.22254963289554e-36;
      ke_kl_data[12] = -8.63653314557700e-37;
      ke_kl_data[13] =  2.28995901944352e-37;
    break;
    case 3: 
      Te_data[0] =  618.908000328523;
      Te_data[1] =  2262.88237620116;
      Te_data[2] =  3371.50133178963;
      Te_data[3] =  4438.34399735592;
      Te_data[4] =  6282.68983833492;
      Te_data[5] =  12269.8511065130;
      Te_data[6] =  19696.7471104552;
      Te_data[7] =  28887.5309153338;
      Te_data[8] =  37250.5252697730;
      Te_data[9] =  45319.5383240561;
      Te_data[10] = 62679.9077332711;
      Te_data[11] = 156506.360583075;
      Te_data[12] = 373510.978198263;
      Te_data[13] = 1059879.95056260;

      ke_kl_data[0] =  -9.02458032633180e-42;
      ke_kl_data[1] =  -5.06243404007705e-35;
      ke_kl_data[2] =  -2.30178356365543e-34;
      ke_kl_data[3] =  -4.69510780972978e-34;
      ke_kl_data[4] =  -8.52275296973685e-34;
      ke_kl_data[5] =  -1.19595436150910e-33;
      ke_kl_data[6] =  -1.05864013352746e-33;
      ke_kl_data[7] =  -8.25496196277683e-34;
      ke_kl_data[8] =  -6.40319504837904e-34;
      ke_kl_data[9] =  -4.99215968899187e-34;
      ke_kl_data[10] = -3.04855600993321e-34;
      ke_kl_data[11] = -7.33050444278219e-35;
      ke_kl_data[12] = -2.24841457932390e-35;
      ke_kl_data[13] = -5.95457117126933e-36;
    break;
    case 4: 
      Te_data[0] =  618.908000328523;
      Te_data[1] =  2262.88237620116;
      Te_data[2] =  3371.50133178963;
      Te_data[3] =  4438.34399735592;
      Te_data[4] =  6282.68983833492;
      Te_data[5] =  12269.8511065130;
      Te_data[6] =  19696.7471104552;
      Te_data[7] =  28887.5309153338;
      Te_data[8] =  37250.5252697730;
      Te_data[9] =  45319.5383240561;
      Te_data[10] = 62679.9077332711;
      Te_data[11] = 156506.360583075;
      Te_data[12] = 373510.978198263;
      Te_data[13] = 1059879.95056260;

      ke_kl_data[0] =  -1.78213311353088e-41;
      ke_kl_data[1] =  -1.29077325753885e-34;
      ke_kl_data[2] =  -4.53261377376819e-34;
      ke_kl_data[3] =  -8.15407930884039e-34;
      ke_kl_data[4] =  -1.33123734995416e-33;
      ke_kl_data[5] =  -1.75979076036975e-33;
      ke_kl_data[6] =  -1.57887617921174e-33;
      ke_kl_data[7] =  -1.26480934431422e-33;
      ke_kl_data[8] =  -1.00712743104408e-33;
      ke_kl_data[9] =  -8.05337609868695e-34;
      ke_kl_data[10] = -5.17348875324510e-34;
      ke_kl_data[11] = -1.39626681561028e-34;
      ke_kl_data[12] = -4.58876413673634e-35;
      ke_kl_data[13] = -1.41882305935205e-35;
    break;
    default:
      fatal_error("kth N2vib process level in _ke_kl_from_Te() set to invalid value.");
  }
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( Te , Te_data[0] ) ) );
  ke_kl = EXM_f_from_monotonespline(N, Te_data, ke_kl_data, Te); 
  return(ke_kl);
}

/* find the fraction fk of vibrationally excited N2 level k wrt total N2 (ground state + vib. excited states)
 * assuming a Boltzmann distribution in Tv
 * energy e_k in eV for the transition from ground state to vib. level k */
double _fk_N2vib_from_Tv(long k, double Tv) {
  long l,nvib;
  double Zvib,fk;
  static const double e_k[] = {0.288, 0.572, 0.855, 1.133,1.4082,
  1.6794, 1.9470, 2.2111, 2.4717, 2.7287,
  2.9821, 3.2320, 3.4782, 3.7208, 3.9598,
  4.1951, 4.4268, 4.6547, 4.8789, 5.0993}; // eV
  nvib = 4;
  if (k < 1 || k > nvib) fatal_error("Invalid vibrational energy level in _fk_from_Tv().");
  // find partition function to ensure fractions add up to unity
  Zvib=0.0;
  for (l=1; l<=nvib; l++) Zvib+=exp(-fabs(_C(speceminus))*e_k[l-1]/(kB*Tv));
  Zvib+=1.0;
  // find fraction of vib. excited N2 state k
  fk=exp(-fabs(_C(speceminus))*e_k[k-1]/(kB*Tv))/Zvib;
  return(fk);
}

/* find Qe from excitation and de-excitation of vibrational N2 levels by electron impact
 * and include Qe from ground state N2VO = N2*(1-sum_fk)*/
double _QeN2vib_from_rhok_T_Tv_Te(gl_t *gl, spec_t rhok, double T, double Tv, double Te){
  long nvib,viblevel,sum;
  double Qe,Tref;
  nvib=4;
  Tref=300.0e0; //gas temperature at which the experiments or BOLSIG+ results of mueN and E/N were obtained
  Qe=0.0;
  sum=0.0;
  for (viblevel=1; viblevel<=nvib; viblevel++) sum+=_fk_N2vib_from_Tv(viblevel,Tv);
  assert(sum<1.0e0);
  // electron cooling-heating source term from N2 vibrational levels
  for (viblevel=1; viblevel<=nvib; viblevel++) Qe+=rhok[speceminus]/_m(speceminus)*rhok[specN2]/_m(specN2)*_fk_N2vib_from_Tv(viblevel,Tv)*_ke_kl_N2vib_from_Te(viblevel,Te);
  // electron cooling-heating source term from N2 in ground state
  Qe+=rhok[speceminus]/_m(speceminus)*rhok[specN2]/_m(specN2)*fabs(_C(speceminus))*_mueNk_from_Te(gl,specN2,Te)*max(sqr(_EoverNk_from_Te(specN2, Te))-3.0*kB*max(0.0,Te-Tref)/(_m(specN2)*sqr(_mueNk_from_Te(gl,specN2,Te))),0.0)*(1.0-sum);
  Qe+=3.0*kB*fabs(_C(speceminus))*(Te-T)*rhok[speceminus]*rhok[specN2]/(_m(speceminus)*sqr(_m(specN2))*_mueNk_from_Te(gl,specN2,Te))*(1.0-sum);
  return(Qe);
}

/* find derivatives of Qe from excitation and de-excitation of vibrational N2 levels by electron impact
 * and include the derivatives of Qe from ground state N2VO = N2*(1-sum_fk)*/
void find_dQeN2vib_from_rhok_T_Tv_Te(gl_t *gl, spec_t rhok, double T, double Tv, double Te, spec_t dQeN2vibdrhok, double *dQeN2vibdT, double *dQeN2vibdTv, double *dQeN2vibdTe){
  long spec,viblevel,nvib;
  double dTv,dTe,sum;
  dTe=10.0;
  dTv=dTe;
  nvib=4;
  sum=0.0;
  for (viblevel=1; viblevel<=nvib; viblevel++) sum+=_fk_N2vib_from_Tv(viblevel,Tv);
  for (spec=0; spec<ns; spec++){
    if (spec==specN2 || spec==speceminus) dQeN2vibdrhok[spec] = _QeN2vib_from_rhok_T_Tv_Te(gl,rhok,T,Tv,Te)/rhok[spec];
    else dQeN2vibdrhok[spec]=0.0;
  }
  *dQeN2vibdT=3.0*kB*fabs(_C(speceminus))*rhok[speceminus]*rhok[specN2]/(_m(speceminus)*sqr(_m(specN2))*_mueNk_from_Te(gl,specN2,Te))*(1.0-sum);
  *dQeN2vibdTv=(_QeN2vib_from_rhok_T_Tv_Te(gl,rhok,T,Tv+dTv,Te)-_QeN2vib_from_rhok_T_Tv_Te(gl,rhok,T,Tv,Te))/dTv;
  *dQeN2vibdTe=(_QeN2vib_from_rhok_T_Tv_Te(gl,rhok,T,Tv,Te+dTe)-_QeN2vib_from_rhok_T_Tv_Te(gl,rhok,T,Tv,Te))/dTe;
}