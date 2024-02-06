// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2024 Bernard Parent
Copyright 2024 Felipe Martin Rodriguez Fuentes


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

static double _EoverN_from_Te_H_Morgan(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    6.10955759023806,
    6.11058825252976,
    6.11264639578264,
    6.11905092644849,
    6.13551324579841,
    6.17732947085113,
    6.26932719081473,
    6.44204443466952,
    6.71252448865311,
    7.08223304152285,
    7.55126153661038,
    8.13156771780325,
    8.84331222217878,
    9.70203062068486,
    10.3866246486031,
    10.6277867054200,
    10.8077150746143,
    11.0309995603238,
    11.3463395697431,
    11.7988853189096,
    12.2742147582549,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -52.9710237752345,
    -52.0336121941841,
    -51.5130088549234,
    -50.9924846520783,
    -50.4720536088765,
    -49.9508086400425,
    -49.4304530527132,
    -48.9097613283098,
    -48.3889819226718,
    -47.8681639417493,
    -47.3477899531147,
    -46.8269271514765,
    -46.3062070876381,
    -45.7854988191063,
    -45.2650643362336,
    -44.7441801799264,
    -44.2234497542938,
    -43.7031878349985,
    -43.1820998634711,
    -42.6615654255782,
    -42.2448171725761,
    -40.3479193852247
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}

static double _Te_from_EoverN_H_Morgan(double EoverN){
  double Te;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -52.9710237752345,
    -52.0336121941841,
    -51.5130088549234,
    -50.9924846520783,
    -50.4720536088765,
    -49.9508086400425,
    -49.4304530527132,
    -48.9097613283098,
    -48.3889819226718,
    -47.8681639417493,
    -47.3477899531147,
    -46.8269271514765,
    -46.3062070876381,
    -45.7854988191063,
    -45.2650643362336,
    -44.7441801799264,
    -44.2234497542938,
    -43.7031878349985,
    -43.1820998634711,
    -42.6615654255782,
    -42.2448171725761,
    -39.1439465808988
  };
  /* log K */
  double Te_data[] = 
  { 
    6.10955759023806,
    6.11058825252976,
    6.11264639578264,
    6.11905092644849,
    6.13551324579841,
    6.17732947085113,
    6.26932719081473,
    6.44204443466952,
    6.71252448865311,
    7.08223304152285,
    7.55126153661038,
    8.13156771780325,
    8.84331222217878,
    9.70203062068486,
    10.3866246486031,
    10.6277867054200,
    10.8077150746143,
    11.0309995603238,
    11.3463395697431,
    11.7988853189096,
    12.2742147582549,
    15.4249484703984
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


static double _EoverN_from_Te_N2(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    5.78359961896139,
    5.88538231327133,
    6.00674317027560,
    6.18906472706955,
    7.05656529477427,
    7.52656892402001,
    7.74971247533422,
    8.08618471195543,
    8.86485406595354,
    9.05804529498440,
    9.18479700062354,
    9.28657969493348,
    9.54147194456227,
    9.69562262438953,
    9.76461549587648,
    10.0010042739407,
    10.6516833690737,
    12.7963582069535,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -51.8608448501949,
    -51.3500192264290,
    -51.0135469898077,
    -50.6568720458690,
    -49.5582597572009,
    -49.0474341334349,
    -48.7109618968137,
    -48.3542869528750,
    -47.2556746642069,
    -46.7448490404409,
    -46.4083768038197,
    -46.0517018598809,
    -44.9530895712128,
    -44.4422639474468,
    -44.1057917108256,
    -43.7491167668869,
    -43.3907432663126,
    -40.6829260296843,
    -39.1439465808988
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
    -51.3500192264290,
    -51.0135469898077,
    -50.6568720458690,
    -49.5582597572009,
    -49.0474341334349,
    -48.7109618968137,
    -48.3542869528750,
    -47.2556746642069,
    -46.7448490404409,
    -46.4083768038197,
    -46.0517018598809,
    -44.9530895712128,
    -44.4422639474468,
    -44.1057917108256,
    -43.7491167668869,
    -43.7031878349985,
    -43.4944745485133,
    -42.7655412882071,
    -40.6829260296843,
    -32.9293384824766,
  };
  /* log K */
  double Te_data[] = 
  { 
    5.78359961896139,
    5.88538231327133,
    6.00674317027560,
    6.18906472706955,
    7.05656529477427,
    7.52656892402001,
    7.74971247533422,
    8.08618471195543,
    8.86485406595354,
    9.05804529498440,
    9.18479700062354,
    9.28657969493348,
    9.54147194456227,
    9.69562262438953,
    9.76461549587648,
    10.0010042739407,
    10.2371160053208,
    10.5247980777726,
    11.1894893389638,
    12.7963582069535,
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
    5.65766872569018,
    5.88610403145016,
    6.47697236288968,
    7.39303753139549,
    8.15517758344238,
    8.64580049989085,
    9.04443964292862,
    9.45446056757264,
    10.2346191251222,
    10.4238611247607,
    11.1022446850151,
    11.5584838509519,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  {
    -50.6568720458690,
    -50.0690853809669,
    -49.4041090773736,
    -48.3542869528750,
    -47.2556746642069,
    -46.7448490404409,
    -46.4083768038197,
    -46.0517018598809,
    -44.9530895712128,
    -44.4422639474468,
    -43.1820998634711,
    -42.4532943089311,
    -39.1439465808988
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
    -50.0690853809669,
    -49.4041090773736,
    -48.3542869528750,
    -47.2556746642069,
    -46.7448490404409,
    -46.4083768038197,
    -46.0517018598809,
    -44.9530895712128,
    -44.4422639474468,
    -42.6615654255782,
    -42.4532943089311,
    -32.9293384824766
  };
  /* log K */
  double Te_data[] = 
  { 
    5.65766872569018,
    5.88610403145016,
    6.47697236288968,
    7.39303753139549,
    8.15517758344238,
    8.64580049989085,
    9.04443964292862,
    9.45446056757264,
    10.2346191251222,
    10.4238611247607,
    11.4209369942094,
    11.5584838509519,
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
    2.30258509299405,
    3.15709869850392,
    4.15325682806993,
    5.60445309962002,
    6.40079140504183,
    6.90775527898214,
    8.07432899557906,
    8.67059265954218,
    9.54396882476086,
    10.3367243480594,
    11.1840222084466,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  {
    -60.9779329059717,
    -57.3277129484975,
    -54.4138087601033,
    -51.4008648205280,
    -50.0932061213430,
    -49.3445762009764,
    -45.5773327731255,
    -45.1609340852378,
    -44.6399834464209,
    -43.5989741084572,
    -42.5573177634926,
    -38.0453342922307
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}



double _EoverNk_from_Te(long spec, double Te){
  double Estar;
  switch (smap[spec]){
    case SMAP_NH3:
      Estar=_EoverN_from_Te_NH3_Morgan(Te);
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
    default:
      Estar=0.0;
  }
  return(Estar);
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
      break;
      case SMAP_N2:
        Estar+=rhok[spec]/_m(spec)*_EoverN_from_Te_N2(Te);
        Nsum+=rhok[spec]/_m(spec);
      break;
      case SMAP_O2:
        Estar+=rhok[spec]/_m(spec)*_EoverN_from_Te_O2(Te);
        Nsum+=rhok[spec]/_m(spec);
      break;
      case SMAP_NO:
        Estar+=rhok[spec]/_m(spec)*_EoverN_from_Te_NO(Te);
        Nsum+=rhok[spec]/_m(spec);
      break;
      case SMAP_H:
        Estar+=rhok[spec]/_m(spec)*_EoverN_from_Te_H_Morgan(Te);
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
      break;
      case SMAP_N2:
        Te+=rhok[spec]/_m(spec)*_Te_from_EoverN_N2(EoverN);
        Nsum+=rhok[spec]/_m(spec);
      break;
      case SMAP_O2:
        Te+=rhok[spec]/_m(spec)*_Te_from_EoverN_O2(EoverN);
        Nsum+=rhok[spec]/_m(spec);
      break;
      case SMAP_H:
        Te+=rhok[spec]/_m(spec)*_Te_from_EoverN_H_Morgan(EoverN);
        Nsum+=rhok[spec]/_m(spec);
      break;
    }
  }
  if (Nsum==0.0) fatal_error("Couldn't find any valid species when determining Te in _Te_from_Nk_EoverN().");
  Te/=Nsum;
  return(Te);
}


static double _mueN_N2(double Te){
  double mueN_N2;  
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_control[] = 
  { 
    2.30258509299405,
    4.30406509320417,
    6.11777329755360,
    6.16032225820390,
    6.31422787002356,
    6.65309818932878,
    7.12888595635390,
    7.52092804412992,
    8.05462032845107,
    8.75401408453109,
    9.14756478528784,
    9.34312262669112,
    9.44165160928006,
    10.0847647584658,
    11.1894893389638,
    11.9794616745625,
    12.7963582069535
  };
  /* log m2/Vs*/
  double mueN_control[] = 
  { 
    67.4168215829997,
    61.6589718870732,
    58.8604497828070,
    58.7654953142236,
    58.5220600000565,
    58.1127487333608,
    57.6332201163168,
    57.3072804521954,
    56.9306920039230,
    56.4577812312559,
    56.0953860429103,
    55.7885443350081,
    55.5658436861888,
    55.3408534122814,
    55.1394225183092,
    54.8599215249993,
    54.5710926348406
  };
  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  mueN_N2 = EXM_f_from_monotonespline(N, Te_control, mueN_control, Te);
  return(exp( mueN_N2 ));
}



double _mueNk_from_Te(long spec, double Te){
  double mueN;
  switch (smap[spec]){
    case SMAP_N2:
      mueN=_mueN_N2(Te);
    break;
    default:
      mueN=3.74E19*exp(33.5/sqrt(log(Te)));
  }
  return(mueN);
}