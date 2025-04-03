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
    5.92986212593964,
    5.93754619632850,
    5.96025036515485,
    6.01951836176432,
    6.10955759023806,
    6.28826344061357,
    7.72182727184468,
    9.69918769055403,
    10.3722047880748,
    10.6161899130309,
    10.7970879820400,
    10.9836739376205,
    11.1922117363535,
    11.4268988308571,
    11.6922646857972,
    12.0017727835481,
    12.3707560104786,
    12.7924922227570,
    13.2499497567448,
    13.7387991976562
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -48.6234744426906,
    -48.2329546677074,
    -46.6708553605325,
    -46.2804862910852,
    -45.8895830104045,
    -45.4995423726219,
    -45.1089639579919,
    -44.7182809883681,
    -44.3278644071383,
    -43.9373759163476,
    -43.5469925827967,
    -43.1563425604737,
    -42.7659153297301,
    -42.3754011879738,
    -41.9848998006652,
    -41.5943358686612,
    -41.2035854952824,
    -40.8131344977387,
    -40.4230021913491,
    -40.0323789846445
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
    -48.6234744426906,
    -48.2329546677074,
    -46.6708553605325,
    -46.2804862910852,
    -45.8895830104045,
    -45.4995423726219,
    -45.1089639579919,
    -44.7182809883681,
    -44.3278644071383,
    -43.9373759163476,
    -43.5469925827967,
    -43.1563425604737,
    -42.7659153297301,
    -42.3754011879738,
    -41.9848998006652,
    -41.5943358686612,
    -41.2035854952824,
    -40.8131344977387,
    -40.4230021913491,
    -40.0323789846445
  };
  /* log K */
  double Te_data[] = 
  { 
    5.92986212593964,
    5.93754619632850,
    5.96025036515485,
    6.01951836176432,
    6.10955759023806,
    6.28826344061357,
    7.72182727184468,
    9.69918769055403,
    10.3722047880748,
    10.6161899130309,
    10.7970879820400,
    10.9836739376205,
    11.1922117363535,
    11.4268988308571,
    11.6922646857972,
    12.0017727835481,
    12.3707560104786,
    12.7924922227570,
    13.2499497567448,
    13.7387991976562
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
    6.22805647628365,
    6.29620936452088,
    6.34931918983483,
    6.39414128685717,
    6.45781878963581,
    6.55799756840379,
    6.76085105062523,
    7.68196882324770,
    9.27337461074692,
    10.0397187861214,
    10.3329550258431,
    10.5024173169769,
    10.6503095537732,
    10.8079499542396,
    10.9872137643256,
    11.1928512247845,
    11.4262662797575,
    11.6912942830039,
    12.0003482818827,
    12.3687855664913,
    12.7915233868989,
    13.2489277842491,
    13.7382978816914
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -48.6234744426906,
    -48.2329546677074,
    -47.8426616489384,
    -47.4520951453803,
    -47.0614794348485,
    -46.6708553605325,
    -46.2804862910852,
    -45.8895830104045,
    -45.4995423726219,
    -45.1089639579919,
    -44.7182809883681,
    -44.3278644071383,
    -43.9373759163476,
    -43.5469925827967,
    -43.1563425604737,
    -42.7659153297301,
    -42.3754011879738,
    -41.9848998006652,
    -41.5943358686612,
    -41.2035854952824,
    -40.8131344977387,
    -40.4230021913491,
    -40.0323789846445
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
    -48.6234744426906,
    -48.2329546677074,
    -47.8426616489384,
    -47.4520951453803,
    -47.0614794348485,
    -46.6708553605325,
    -46.2804862910852,
    -45.8895830104045,
    -45.4995423726219,
    -45.1089639579919,
    -44.7182809883681,
    -44.3278644071383,
    -43.9373759163476,
    -43.5469925827967,
    -43.1563425604737,
    -42.7659153297301,
    -42.3754011879738,
    -41.9848998006652,
    -41.5943358686612,
    -41.2035854952824,
    -40.8131344977387,
    -40.4230021913491,
    -40.0323789846445
  };
  /* log K */
  double Te_data[] = 
  { 
    6.22805647628365,
    6.29620936452088,
    6.34931918983483,
    6.39414128685717,
    6.45781878963581,
    6.55799756840379,
    6.76085105062523,
    7.68196882324770,
    9.27337461074692,
    10.0397187861214,
    10.3329550258431,
    10.5024173169769,
    10.6503095537732,
    10.8079499542396,
    10.9872137643256,
    11.1928512247845,
    11.4262662797575,
    11.6912942830039,
    12.0003482818827,
    12.3687855664913,
    12.7915233868989,
    13.2489277842491,
    13.7382978816914
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


/* 
Compute the influence of the ionization degree on E/N vs Te for N2
E/N vs Te splines are obtained for electron molar fractions of 1e-6 to 1e-2
Given the electron molar fraction chieminus, interpolate at a given Te to find E/N 
*/
static double _EoverN_from_Te_chieminus_N2(double Te, double chieminus){
  long k;
  double chie_data[] = {1.000E-6,1.000E-5,1.000E-4,1.000E-3,1.000E-2};
  int N_chie = sizeof(chie_data)/sizeof(chie_data[0]);
  double EoverN_chie_data[N_chie],Tek[N_chie],EoverN,Te_data[12],EoverN_data[12];
  // spline data at chieminus = 1e-6,1e-5,1e-4,1e-3,1e-2
  for (k=0; k<N_chie; k++){  
    switch (k){
      case 0:
        Te_data[0]=5.7683209957938; Te_data[1]=5.9895807688236; Te_data[2]=6.5808712966722; Te_data[3]=7.1628887632684;
        Te_data[4]=8.0619851978272; Te_data[5]=8.4079571483189; Te_data[6]=9.0635173434530; Te_data[7]=9.8359370093448;
        Te_data[8]=10.2153848302326; Te_data[9]=11.3844718843909; Te_data[10]=12.2033581841564; Te_data[11]=13.8155105579643;
        
        EoverN_data[0]=-51.8608448501949; EoverN_data[1]=-51.0557104756815; EoverN_data[2]=-50.1104688751939; EoverN_data[3]=-49.3864420173586;
        EoverN_data[4]=-48.3392573343656; EoverN_data[5]=-47.8592084179916; EoverN_data[6]=-46.6399631066307; EoverN_data[7]=-44.0774427408722;
        EoverN_data[8]=-43.4726005050452; EoverN_data[9]=-42.4528700839551; EoverN_data[10]=-41.3655827304183; EoverN_data[11]=-39.8080228556486;
      break;
      case 1:
        Te_data[0]=5.7683209957938; Te_data[1]=5.9895807688236; Te_data[2]=6.5808712966722; Te_data[3]=7.1628887632684;
        Te_data[4]=8.0619851978272; Te_data[5]=8.4079571483189; Te_data[6]=9.0635173434530; Te_data[7]=9.6333495008421;
        Te_data[8]=10.2153848302326; Te_data[9]=11.3844718843909; Te_data[10]=12.2033581841564; Te_data[11]=13.8155105579643;
        
        EoverN_data[0]=-51.8608448501949; EoverN_data[1]=-51.0557104756815; EoverN_data[2]=-50.1085328420543; EoverN_data[3]=-49.3801397401679;
        EoverN_data[4]=-48.3083794171069; EoverN_data[5]=-47.7656757851313; EoverN_data[6]=-46.3349910379067; EoverN_data[7]=-44.5138801470076;
        EoverN_data[8]=-43.4726005050452; EoverN_data[9]=-42.4528700839551; EoverN_data[10]=-41.3655827304183; EoverN_data[11]=-39.8080228556486;
      break;
      case 2:
        Te_data[0]=5.7683209957938; Te_data[1]=5.9895807688236; Te_data[2]=6.5808712966722; Te_data[3]=7.1628887632684;
        Te_data[4]=8.0619851978272; Te_data[5]=8.4079571483189; Te_data[6]=9.0635173434530; Te_data[7]=9.5975225839140;
        Te_data[8]=10.2153848302326; Te_data[9]=11.3844718843909; Te_data[10]=12.2033581841564; Te_data[11]=13.8155105579643;
        
        EoverN_data[0]=-51.8608448501949; EoverN_data[1]=-51.0557104756815; EoverN_data[2]=-50.1069099506727; EoverN_data[3]=-49.3727751273235;
        EoverN_data[4]=-48.2419739062285; EoverN_data[5]=-47.4771310463929; EoverN_data[6]=-45.6569866292733; EoverN_data[7]=-44.3535890924510;
        EoverN_data[8]=-43.4775996021468; EoverN_data[9]=-42.4528700839551; EoverN_data[10]=-41.3655827304183; EoverN_data[11]=-39.8080228556486;
      break;
      case 3:
        Te_data[0]=5.7683209957938; Te_data[1]=5.9895807688236; Te_data[2]=6.5808712966722; Te_data[3]=7.1628887632684;
        Te_data[4]=8.0619851978272; Te_data[5]=8.4079571483189; Te_data[6]=9.0635173434530; Te_data[7]=9.8359370093448;
        Te_data[8]=10.2153848302326; Te_data[9]=11.3844718843909; Te_data[10]=12.2033581841564; Te_data[11]=13.8155105579643;
        
        EoverN_data[0]=-51.8608448501949; EoverN_data[1]=-51.0557104756815; EoverN_data[2]=-50.1064620441829; EoverN_data[3]=-49.3691609881841;
        EoverN_data[4]=-48.1138496265006; EoverN_data[5]=-47.0154858747258; EoverN_data[6]=-45.0964427747303; EoverN_data[7]=-43.9094349253216;
        EoverN_data[8]=-43.4724699527492; EoverN_data[9]=-42.4528700839551; EoverN_data[10]=-41.3655827304183; EoverN_data[11]=-39.8080228556486;
      break;
      case 4:
        Te_data[0]=5.7683209957938; Te_data[1]=5.9895807688236; Te_data[2]=6.5808712966722; Te_data[3]=7.1628887632684;
        Te_data[4]=8.0619851978272; Te_data[5]=8.4079571483189; Te_data[6]=9.0635173434530; Te_data[7]=9.8359370093448;
        Te_data[8]=10.2153848302326; Te_data[9]=11.3844718843909; Te_data[10]=12.2033581841564; Te_data[11]=13.8155105579643;

        EoverN_data[0]=-51.8608448501949; EoverN_data[1]=-51.0557104756815; EoverN_data[2]=-50.1064682435513; EoverN_data[3]=-49.3685992138692;
        EoverN_data[4]=-48.0548406010350; EoverN_data[5]=-46.8277415123705; EoverN_data[6]=-44.8901986329883; EoverN_data[7]=-43.8229867640502;
        EoverN_data[8]=-43.4632414442052; EoverN_data[9]=-42.4528700839551; EoverN_data[10]=-41.3655827304183; EoverN_data[11]=-39.8080228556486;
      break;
      default:
        fatal_error("Couldn't find EoverN from electron molar fraction data in _EoverN_from_Te_chieminus_N2().");
    }
    int N = sizeof(Te_data)/sizeof(Te_data[0]);
    Tek[k] = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
    // inerpolate over Te to find E/N at each ionization fraction data
    EoverN_chie_data[k] = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Tek[k]) );
  }
  // interpolate over chieminus to find final E/N
  EoverN = EXM_f_from_monotonespline(N_chie, chie_data, EoverN_chie_data, max(chie_data[0],min( chieminus, chie_data[N_chie-1])));
  return(EoverN);
};

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


double _EoverNk_from_Te_chieminus(long spec, double Te, double chieminus){
  double Estar;
  switch (smap[spec]){
    case SMAP_N2:
      Estar=_EoverN_from_Te_chieminus_N2(Te, chieminus);
    break;
    default:
      Estar=_EoverNk_from_Te(spec, Te);
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