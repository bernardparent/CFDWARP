// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Ajjay Omprakas, Prasanna Thoguluva Rajendran
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
#include <model/share/chem_share.h>
#include <model/share/model_share.h>


/* set all reactions to true except for testing purposes */
const static bool REACTION[36]=
  {
   FALSE, /* reaction 0 */
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
   TRUE, /* reaction 13 */ /* modified from original kim2021.c */
   TRUE, /* reaction 14 */ /* modified from original kim2021.c */
   TRUE, /* reaction 15 */ /* modified from original kim2021.c */
   TRUE, /* reaction 16 */ /* modified from original kim2021.c */
   TRUE, /* reaction 17 */
   TRUE, /* reaction 18 */
   TRUE, /* reaction 19 */ 
   TRUE, /* reaction 20 */
   TRUE, /* reaction 21 */
   TRUE, /* reaction 22 */
   TRUE, /* reaction 23 */
   TRUE, /* reaction 24 */
   TRUE, /* reaction 25 */
   TRUE, /* reaction 26 */ /* modified from original kim2021.c */
   TRUE, /* reaction 27 */
   TRUE, /* reaction 28 */
   TRUE, /* reaction 29 */ /* modified from original kim2021.c */
   TRUE, /* reaction 30 */ /* N2 ionization to N2+ added for completeness, also included in Qei */
   TRUE, /* reaction 31 */ /* O2 ionization to O2+ added for completeness, also included in Qei */ 
   TRUE, /* reaction 32 */ /* NO ionization to NO+ added for completeness, also included in Qei */
   TRUE, /* reaction 33 */ /* N2+ EI recombination to N2 */
   TRUE, /* reaction 34 */ /* O2+ EI recombination to O2 */
   TRUE, /* reaction 35 */ /* NO+ EI recombination to NO */
  };

#define specEND -1


const static long specM1[]=
  {
   specN2, specN2plus, specEND
  };

const static long specM2[]=
  {
   specO2, specNO, specO2plus, specNOplus, specEND
  };

const static long specM3[]=
  {
   specN, specNplus, specEND
  };

const static long specM4[]=
  {
   specO, specOplus, specEND
  };
  
const static long specM5[]=
  {
   specN2, specNO, specN2plus, specNOplus, specEND
  };
  
const static long specM6[]=
  {
   specO2, specO2plus, specEND
  };
  
const static long specM7[]=
  {
   specN2, specO2, specN2plus, specO2plus, specEND
  };
  
const static long specM8[]=
  {
   specNO, specN, specO, specNOplus, specNplus, specOplus, specEND
  };

double _kb_polynomial_rodriguez2024(long numreactant,double A1, double A2, double A3, double A4, double A5, double A, double n, double E, double T)
{
  double ke,kf,kb,Z, R;
  R=Rchem;
  Z = 10000.0/T;
  ke = exp(A1/Z + A2 + A3*log(Z) + A4*Z + A5*Z*Z);
  switch (numreactant){
    case 2:
      kf=A/calA*pow(T,n)*exp(-E/(R*T));
    break;
    case 3:
      kf=A/sqr(calA)*pow(T,n)*exp(-E/(R*T));    
    break;
    default:
      fatal_error("numreactant can not be set to %ld in _kb_polynomial_rodriguez2024 within rodriguez2024.c",numreactant);
      kf=0.0;
  }
  assert(ke!=0.0);
  kb = kf/ke;
  
  return(kb);
  
}

double _dkbdT_polynomial_rodriguez2024(long numreactant,double A1, double A2, double A3, double A4, double A5, double A, double n, double E, double T)
{
  double dkbdT,R;
  R=Rchem;
  
  switch (numreactant){
    case 2:
      dkbdT = (A*pow(T,n)*exp(-E/(R*T))*exp(- A2 - (A1*T)/10000.0 - A3*log(10000.0/T) - (10000.0*A4)/T - (100000000.0*A5)/pow(T,2.0))*(A3/T - A1/10000.0 + (10000.0*A4)/pow(T,2.0) + (200000000.0*A5)/pow(T,3.0)))/calA + (A*n*pow(T,(n - 1.0))*exp(-E/(R*T))*exp(- A2 - (A1*T)/10000.0 - A3*log(10000.0/T) - (10000.0*A4)/T - (100000000.0*A5)/pow(T,2.0)))/calA + (A*E*pow(T,n)*exp(-E/(R*T))*exp(- A2 - (A1*T)/10000.0 - A3*log(10000.0/T) - (10000.0*A4)/T - (100000000.0*A5)/pow(T,2.0)))/(R*calA*pow(T,2.0));
    break;
    case 3:
      dkbdT =  (A*pow(T,n)*exp(-E/(R*T))*exp(- A2 - (A1*T)/10000.0 - A3*log(10000.0/T) - (10000.0*A4)/T - (100000000.0*A5)/pow(T,2))*(A3/T - A1/10000.0 + (10000.0*A4)/pow(T,2) + (200000000.0*A5)/pow(T,3.0)))/sqr(calA) + (A*n*pow(T,(n - 1))*exp(-E/(R*T))*exp(- A2 - (A1*T)/10000.0 - A3*log(10000.0/T) - (10000.0*A4)/T - (100000000.0*A5)/pow(T,2.0)))/sqr(calA) + (A*E*pow(T,n)*exp(-E/(R*T))*exp(- A2 - (A1*T)/10000.0 - A3*log(10000.0/T) - (10000.0*A4)/T - (100000000.0*A5)/pow(T,2.0)))/(R*sqr(calA)*pow(T,2.0));
    break;
    default:
      fatal_error("numreactant can not be set to %ld in _dkbdT_polynomial_rodriguez2024 within rodriguez2024.c",numreactant);
      dkbdT=0.0;
  }
 
  return(dkbdT);
  
}

static double _kf30(np_t np, gl_t *gl, double Te){
  /* N2 ionization */
  double kf30;
  double Te_control[] = 
  { 
    9.46979690785538,
    9.51956710917422,
    9.84281034346154,
    10.1366387167620,
    10.7022762534868,
    11.6533000019413,
    11.8041027220493,
    12.6902842369606,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -41.1651192144546,
    -37.0450701265803,
    -30.9051454150272,
    -28.6391530514683,
    -25.2959688327970,
    -19.9307259125084,
    -19.3648566793161,
    -17.2892786324613,
    -15.7126305428502
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf30 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf30 ));
}

static double _kf31(np_t np, gl_t *gl, double Te){
  /* O2 ionization */
  double kf31;
  double Te_control[] = 
  { 
    10.2400211744112,
    10.4029544399414,
    10.5381130290901,
    10.6701822643877,
    10.7424402398679,
    10.9126526678721,
    11.2567702476959,
    11.5315128937124,
    13.1215130159495,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -44.0086650000949,
    -36.3676148723802,
    -31.2399809468666,
    -27.6637509405515,
    -26.2590792610331,
    -23.9600510823582,
    -21.3796954345712,
    -20.0957244134845,
    -16.5643827535867,
    -15.2808481264246
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf31 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf31 ));
}

static double _kf32(np_t np, gl_t *gl, double Te){
  /* NO ionization */
  double kf32;
  double Te_control[] = 
  { 
    8.48735489578072,
    8.88417441689603,
    9.26307991901759,
    9.56535121835222,
    9.79028324262506,
    10.0982265002167,
    10.5172886557505,
    11.6903229376143,
    11.8306340172242,
    13.1133488685331,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -43.8225480449728,
    -39.5269655001771,
    -36.0264392504189,
    -33.1178390903473,
    -30.6688261659045,
    -26.8094807206002,
    -22.8654342085346,
    -18.6548852426260,
    -18.3163207286281,
    -16.4545478877795,
    -15.5993018575432
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf32 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf32 ));
}


static double _kf15(np_t np, gl_t *gl, double Te){
  /* N ionization */
  double kf15;
  double Te_control[] = 
  {  
    9.97325436098781,
    10.1148027736797,
    10.2972026114240,
    10.5460787574850,
    11.0330643193868,
    12.1594735791518,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -49.6819351185621,
    -39.5250605981431,
    -32.5921522151217,
    -27.7250021793994,
    -23.2627131533956,
    -18.8262458570609,
    -16.4747705948971
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf15 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf15 ));
}

static double _kf16(np_t np, gl_t *gl, double Te){
  /* O ionization */
  double kf16;
  double Te_control[] = 
  {  
    9.74441278855896,
    9.94024806451964,
    10.2101556529438,
    10.5776176382447,
    11.0606203856489,
    12.0830741580182,
    13.3653926898139,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -45.7081121554908,
    -36.5368222983865,
    -30.2191572203654,
    -26.0516301248467,
    -22.8011086572626,
    -18.4691359642218,
    -16.8692957629669,
    -16.3412392022725
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf16 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf16 ));
}


static double _kf26b(np_t np, gl_t *gl, double Te){
  double kf26b;
  double Te_control[] = 
  { 
    2.96134635039148,
    5.71525385948086,
    9.17412056693579,
    9.86376375536634,
    10.2874072392099,
    11.4501196047481,
    12.1488911688257,
    12.3695988451580,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -15.1096633905206,
    -15.6250837536637,
    -16.9584559651989,
    -17.4850816245454,
    -18.2694199501118,
    -18.4823072424797,
    -18.7176080892025,
    -18.7960420527258,
    -19.8069751050723
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf26b = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf26b ));
}

static double _kf29(np_t np, gl_t *gl, double Te){
  double kf29;
  double Te_control[] = 
  { 
    9.56616455715265,
    9.64207114311839,
    9.84834371140720,
    9.98080156662319,
    10.2742403990660,
    10.8311636242907,
    10.9936710805699,
    11.1448917846708,
    11.4428329237054,
    11.6893506477953,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -45.6975300461603,
    -41.8099435457700,
    -35.8671792776148,
    -33.6512967448907,
    -30.3852199730567,
    -26.0004089682630,
    -24.6913302332582,
    -23.4618783436313,
    -21.1865719447444,
    -19.8271778123898,
    -13.9208710736221
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf29 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf29 ));
}


void find_W_Rodriguez2024 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                       double Qbeam, spec_t W ) {
  double N[ns];
  double R,kb,TTv,Tblim,Teblim;
  long k;
  spec_t X;

  /* find properties needed by add_to_W* functions */
  R=1.987192004e0;
  TTv=sqrt(Tv*T);
  Tblim=min(max(300.0,T),32000.0);
  Teblim=min(max(300.0,Te),32000.0);
  
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
  }


  if (REACTION[1]) { 
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 1.216e20, -1.2140, 113200.0*R, TTv, X, W );
      

      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 1.216e20, -1.2140, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 1.216e20, -1.2140, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM1[k], specN, specN,   specN2, specM1[k], kb , N, W);      
    }
  }
      

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, TTv, X, W );
      
      
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 7.0e21, -1.6, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 7.0e21, -1.6, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM2[k], specN, specN,   specN2, specM2[k], kb , N, W); 
    }
  }

  if (REACTION[3]){ 
    for (k=0; specM3[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM3[k],   specN, specN, specM3[k], 3.591e20, -1.226, 113200.0*R, TTv, X, W );
      
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.591e20, -1.226, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.591e20, -1.226, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM3[k], specN, specN,   specN2, specM3[k], kb , N, W);
    }
  }

  if (REACTION[4]){ 
    for (k=0; specM4[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM4[k],   specN, specN, specM4[k], 3.0e22, -1.6, 113200.0*R, TTv, X, W );
      
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.0e22, -1.6, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.0e22, -1.6, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM4[k], specN, specN,   specN2, specM4[k], kb , N, W);
    }
  }

  if (REACTION[5]){
    for (k=0; specM5[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM5[k],   specO, specO, specM5[k], 3.354e15, -0.2726, 59500.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.354e15, -0.2726, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.354e15, -0.2726, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM5[k], specO, specO,   specO2, specM5[k], kb , N, W);
    }
  }

  if (REACTION[6]){
    for (k=0; specM6[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM6[k],   specO, specO, specM6[k], 1.117e25, -2.585, 59500.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.117e25, -2.585, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.117e25, -2.585, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM6[k], specO, specO,   specO2, specM6[k], kb , N, W);
    }
  }

  if (REACTION[7]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM3[k],   specO, specO, specM3[k], 1.0e22, -1.500, 59500.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.0e22, -1.500, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.0e22, -1.500, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM3[k], specO, specO,   specO2, specM3[k], kb , N, W);
    }
  }

  if (REACTION[8]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM4[k],   specO, specO, specM4[k], 3.00e21, -1.50, 59500.0*R, TTv, X, W );
     
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.00e21, -1.50, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.00e21, -1.50, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM4[k], specO, specO,   specO2, specM4[k], kb , N, W);
    }
  }

  if (REACTION[9]){
    for (k=0; specM7[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM7[k],   specN, specO, specM7[k], 1.450e15, 0.0, 75200.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 1.450e15, 0.0, 75200.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 1.450e15, 0.0, 75200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM7[k], specN, specO,   specNO, specM7[k], kb , N, W);
    }
  }

  if (REACTION[10]){
    for (k=0; specM8[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM8[k],   specN, specO, specM8[k], 9.640e14, 0.0, 75200.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 9.640e14, 0.0, 75200.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 9.640e14, 0.0, 75200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM8[k], specN, specO,   specNO, specM8[k], kb , N, W);
    }
  }

  if (REACTION[11]){
    add_to_W_fw_2r2p ( specN2, specO,   specNO, specN, 5.700e12, 0.42, 42938.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -1.789e-1, 1.728, -2.172e-1, -3.733, -2.285e-4, 5.700e12, 0.42, 42938.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 7.441e-2, 2.852, 1.054, -5.303, 1.909e-1, 5.700e12, 0.42, 42938.0*R, Tblim); 
        
      add_to_W_2r2p ( specNO, specN,  specN2, specO, kb , N, W);
  }

  if (REACTION[12]){
    add_to_W_fw_2r2p ( specNO, specO,   specO2, specN, 8.4000e12, 0.0, 19400.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -1.673e-1, -1.390, -1.656e-1, -1.551, -1.102e-4, 8.4000e12, 0.0, 19400.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 1.744e-2, -1.563, 2.354e-1, -1.233, -3.250e-1, 8.4000e12, 0.0, 19400.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specN,  specO, specNO, kb , N, W);
  }

  if (REACTION[13]){
    add_to_W_fw_2r2p ( specN, specO, specNOplus, speceminus, 8.80e08, 1.0, 31900.0*R, T, X, W );
    add_to_W_fw_2r2p ( speceminus, specNOplus, specN, specO, 3.00E-7 * calA * pow ( 300, 0.56 ), -0.56, 0.0 * R, Te, X, W );
  }

  if (REACTION[14]){
    add_to_W_fw_2r2p ( specO, specO, specO2plus, speceminus, 7.10e02, 2.70, 80600.0*R, T, X, W );
    add_to_W_fw_2r2p ( speceminus, specO2plus, specO, specO, 2.40E-7 * calA * pow ( 300, 0.70 ), -0.70, 0.0 * R, Te, X, W );
  }

  if (REACTION[15]){
    add_to_W_2r3p ( speceminus, specN, specNplus, speceminus, speceminus, _kf15(np,gl,Te), N, W );
    add_to_W_fw_3r2p ( specNplus, speceminus, speceminus, speceminus, specN, 2.20E40, -4.5, 0.0 * R, Te, X, W );
  }

  if (REACTION[16]){
    add_to_W_2r3p ( speceminus, specO, specOplus, speceminus, speceminus, _kf16(np,gl,Te), N, W );
    add_to_W_fw_3r2p ( specOplus, speceminus, speceminus, speceminus, specO, 2.20E40, -4.5, 0.0 * R, Te, X, W );
  }

  if (REACTION[17]){
    add_to_W_fw_2r2p ( specN2, specO2plus,   specN2plus, specO2, 9.90e12, 0.0, 40700.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -1.970e-1, 1.031, -2.049e-1, -4.005, -9.866e-5, 9.90e12, 0.0, 40700.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 4.428e-2, 1.462, 6.450e-1, -4.638, -3.600e-2, 9.90e12, 0.0, 40700.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specN2plus,  specO2plus, specN2, kb , N, W);
  }

  if (REACTION[18]){
    add_to_W_fw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.40e13, -1.08, 12800.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -6.300e-2, -5.419e-1, -4.449e-2, -1.266, -1.089e-4, 3.40e13, -1.08, 12800.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 3.518e-2, -7.715e-2, 4.619e-1, -1.923, 9.394e-2, 3.40e13, -1.08, 12800.0*R, Tblim); 
        
      add_to_W_2r2p ( specN2, specOplus,  specNOplus, specN, kb , N, W);
  }

  if (REACTION[19]){
    add_to_W_fw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.00e12, 0.5, 77200.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -4.092e-1, -7.798e-1, -4.273e-1, -7.618, -4.475e-4, 1.00e12, 0.5, 77200.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 1.270e-1, 6.372e-1, 1.751, -9.528, -4.012e-2, 1.00e12, 0.5, 77200.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specNplus,  specNOplus, specO, kb , N, W);
  }

  if (REACTION[20]){
    add_to_W_fw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.40e13, 0.41, 32600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -3.056e-2, 1.762, -8.580e-2, -3.265, -2.168e-4, 2.40e13, 0.41, 32600.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 9.436e-2, 3.786, 1.322, -6.192, 7.723e-1, 2.40e13, 0.41, 32600.0*R, Tblim); 
        
      add_to_W_2r2p ( specNO, specO2plus,  specO2, specNOplus, kb , N, W);
  }

  if (REACTION[21]){
    add_to_W_fw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.20e13, 0.0, 35500.0*R, T, X, W );
   
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -4.862e-2, 1.066, -7.350e-2, -3.537, -8.701e-5, 7.20e13, 0.0, 35500.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 6.423e-2, 2.396, 9.127e-1, -5.527, 5.453e-1, 7.20e13, 0.0, 35500.0*R, Tblim); 
        
      add_to_W_2r2p ( specO, specN2plus,  specN, specNOplus, kb , N, W);
  }

  if (REACTION[22]){
    add_to_W_fw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.70e13, 0.14, 28600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -2.113e-1, -1.152, -1.758e-1, -2.803, -1.205e-4, 8.70e13, 0.14, 28600.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 1.523e-2, -1.586, 1.942e-1, -2.103, -4.874e-1, 8.70e13, 0.14, 28600.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specNplus,  specN, specO2plus, kb , N, W);
  }

  if (REACTION[23]){
    add_to_W_fw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.40e05, 1.90, 26600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -1.673e-1, -1.965, -1.656e-1, -2.619, -1.102e-4, 1.40e05, 1.90, 26600.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 1.744e-2, -2.138, 2.354e-1, -2.302, -3.250e-1, 1.40e05, 1.90, 26600.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specNplus,  specNO, specOplus, kb , N, W);
  }

  if (REACTION[24]){
    add_to_W_fw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.20e12, 0.29, 48600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, -1.978e-1, 3.723e-1, -2.514e-1, -4.815, -3.270e-4, 7.20e12, 0.29, 48600.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 1.118e-1, 2.223, 1.557, -7.425, 4.473e-1, 7.20e12, 0.29, 48600.0*R, Tblim); 
        
      add_to_W_2r2p ( specN, specO2plus,  specO, specNOplus, kb , N, W);
  }

  if (REACTION[25]){
    add_to_W_fw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.10e11, 0.36, 22800.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial_rodriguez2024(2, 1.438e-2, 1.608, -2.901e-2, -2.271, 2.188e-5, 9.10e11, 0.36, 22800.0*R, Tblim);
      else
        kb=_kb_polynomial_rodriguez2024(2, 2.905e-2, 2.473, 4.508e-1, -3.603, 4.514e-1, 9.10e11, 0.36, 22800.0*R, Tblim); 
        
      add_to_W_2r2p ( specO, specN2plus,  specN2, specOplus, kb , N, W);
  }
  
  if (REACTION[26]){
    add_to_W_fw_2r2p ( specN, specN,   specN2plus, speceminus, 4.40e07, 1.5, 67500.0*R, T, X, W );
    add_to_W_2r2p ( speceminus, specN2plus, specN, specN, _kf26b(np,gl,Te), N, W );
  } 
  
  if (REACTION[27]){
    add_to_W_fw_2r2p ( specNplus, specN2,   specN2plus, specN, 1.00e12, 0.5, 12200.0*R, T, X, W );
    add_to_W_fw_2r2p ( specN2plus, specN,   specNplus, specN2,  1.00e12, 0.5, 12200.0*R, T, X, W );
  }  
  
  if (REACTION[28]){
    add_to_W_fw_2r2p ( specO2plus, specO,   specOplus, specO2, 4.00e12, -0.09, 18000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specOplus, specO2,   specO2plus, specO,  4.00e12, -0.09, 18000.0*R, T, X, W );
  }   
  
  if (REACTION[29]){
    add_to_W_2r3p ( speceminus, specN2, specN, specN, speceminus, _kf29(np,gl,Te), N, W );
    add_to_W_fw_3r2p ( specN, specN, speceminus,   speceminus, specN2,   1.20e25, -1.6, 113200.0*R, Te, X, W );
  }   

  if (REACTION[30]){
    add_to_W_2r3p ( speceminus, specN2, specN2plus, speceminus, speceminus, _kf30(np,gl,Te), N, W ); 
  }   

  if (REACTION[31]){
    add_to_W_2r3p ( speceminus, specO2, specO2plus, speceminus, speceminus, _kf31(np,gl,Te), N, W );
  }  

  if (REACTION[32]){
    add_to_W_2r3p ( speceminus, specNO, specNOplus, speceminus, speceminus, _kf32(np,gl,Te), N, W );
  }  

  if (REACTION[33]){
    add_to_W_fw_3r2p ( specN2plus, speceminus, speceminus, speceminus, specN2, 2.20E40, -4.5, 0.0 * R, Te, X, W );
  }

  if (REACTION[34]){
    add_to_W_fw_3r2p ( specO2plus, speceminus, speceminus, speceminus, specO2, 2.20E40, -4.5, 0.0 * R, Te, X, W );
  }

  if (REACTION[35]){
    add_to_W_fw_3r2p ( specNOplus, speceminus, speceminus, speceminus, specNO, 2.20E40, -4.5, 0.0 * R, Te, X, W );
  }

}



/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx_Rodriguez2024 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s, spec;                    
  spec_t N;
  double TTv,TTe,TvTe,Tblim,Teblim,R,kb,dkbdT, dkbdTv, dkbdTe, dkfdTe, dkfdT, dkfdTv;
  spec_t dWdTTv,dWdTTe,dWdTvTe;
  spec_t X;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  
  R=Rchem;
  TTv=sqrt(T*Tv);
  TTe=sqrt(T*Te);
  TvTe=sqrt(Tv*Te);
  Tblim=min(max(300.0,T),32000.0);
  Teblim=min(max(300.0,Te),32000.0);
  /* initialize all derivatives to zero */
  for ( s = 0; s < ns; s++ ) {
    dWdTTv[s] = 0.0;
    dWdTTe[s] = 0.0;
    dWdTvTe[s] = 0.0;
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }

  /* find properties needed by add_to_dW* functions in proper units */
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA;
  }


  if (REACTION[1]) {
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 1.216e20, -1.2140, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 1.216e20, -1.2140, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 1.216e20, -1.2140, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 1.216e20, -1.2140, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 1.216e20, -1.2140, 113200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM1[k], specN, specN,   specN2, specM1[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }
      

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 7.0e21, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 7.0e21, -1.6, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 7.0e21, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 7.0e21, -1.6, 113200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM2[k], specN, specN,   specN2, specM2[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe); 
    }
  }

  if (REACTION[3]){
    for (k=0; specM3[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specN2, specM3[k],   specN, specN, specM3[k], 3.591e20, -1.226, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.591e20, -1.226, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.591e20, -1.226, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.591e20, -1.226, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.591e20, -1.226, 113200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM3[k], specN, specN,   specN2, specM3[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[4]){
    for (k=0; specM4[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specN2, specM4[k],   specN, specN, specM4[k], 3.0e22, -1.6, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.0e22, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.0e22, -1.6, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.0e22, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.0e22, -1.6, 113200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM4[k], specN, specN,   specN2, specM4[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[5]){
    for (k=0; specM5[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM5[k],   specO, specO, specM5[k], 3.354e15, -0.2726, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.354e15, -0.2726, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.354e15, -0.2726, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.354e15, -0.2726, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.354e15, -0.2726, 59500.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM5[k], specO, specO,   specO2, specM5[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[6]){
    for (k=0; specM6[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM6[k],   specO, specO, specM6[k], 1.117e25, -2.585, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.117e25, -2.585, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.117e25, -2.585, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.117e25, -2.585, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.117e25, -2.585, 59500.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM6[k], specO, specO,   specO2, specM6[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[7]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM3[k],   specO, specO, specM3[k], 1.0e22, -1.500, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.0e22, -1.500, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.0e22, -1.500, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.0e22, -1.500, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.0e22, -1.500, 59500.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM3[k], specO, specO,   specO2, specM3[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }



  if (REACTION[8]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM4[k],   specO, specO, specM4[k], 3.00e21, -1.50, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.00e21, -1.50, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.00e21, -1.50, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.00e21, -1.50, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.00e21, -1.50, 59500.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM4[k], specO, specO,   specO2, specM4[k],   kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[9]){
    for (k=0; specM7[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM7[k],   specN, specO, specM7[k], 1.450e15, 0.0, 75200.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 1.450e15, 0.0, 75200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 1.450e15, 0.0, 75200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 1.450e15, 0.0, 75200.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 1.450e15, 0.0, 75200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM7[k], specN, specO,   specNO, specM7[k],  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[10]){
    for (k=0; specM8[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM8[k],   specN, specO, specM8[k], 9.640e14, 0.0, 75200.0*R, TTv, X, dWdTTv, dWdrhok );
      
      if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 9.640e14, 0.0, 75200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 9.640e14, 0.0, 75200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 9.640e14, 0.0, 75200.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 9.640e14, 0.0, 75200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( specM8[k], specN, specO,   specNO, specM8[k],  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    }
  }

  if (REACTION[11]){
    add_to_dW_fw_2r2p ( specN2, specO,   specNO, specN, 5.700e12, 0.42, 42938.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -1.789e-1, 1.728, -2.172e-1, -3.733, -2.285e-4, 5.700e12, 0.42, 42938.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -1.789e-1, 1.728, -2.172e-1, -3.733, -2.285e-4, 5.700e12, 0.42, 42938.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 7.441e-2, 2.852, 1.054, -5.303, 1.909e-1, 5.700e12, 0.42, 42938.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 7.441e-2, 2.852, 1.054, -5.303, 1.909e-1, 5.700e12, 0.42, 42938.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specNO, specN,  specN2, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[12]){
    add_to_dW_fw_2r2p ( specNO, specO,   specO2, specN, 8.4000e12, 0.0, 19400.0*R, T, X, dWdT, dWdrhok);
   
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -1.673e-1, -1.390, -1.656e-1, -1.551, -1.102e-4, 8.4000e12, 0.0, 19400.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -1.673e-1, -1.390, -1.656e-1, -1.551, -1.102e-4, 8.4000e12, 0.0, 19400.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 1.744e-2, -1.563, 2.354e-1, -1.233, -3.250e-1, 8.4000e12, 0.0, 19400.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 1.744e-2, -1.563, 2.354e-1, -1.233, -3.250e-1, 8.4000e12, 0.0, 19400.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specN,  specNO, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }


  if (REACTION[13]){
    add_to_dW_fw_2r2p ( specN, specO, specNOplus, speceminus, 8.80e08, 1.0, 31900.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_fw_2r2p ( specNOplus, speceminus, specN, specO, 3.00E-7 * calA * pow ( 300, 0.56 ), -0.56, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[14]){
    add_to_dW_fw_2r2p ( specO, specO, specO2plus, speceminus, 7.10e02, 2.70, 80600.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_fw_2r2p ( specO2plus, speceminus, specO, specO, 2.40E-7 * calA * pow ( 300, 0.70 ), -0.70, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[15]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( speceminus, specN, specNplus, speceminus, speceminus, _kf15(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specNplus, speceminus, speceminus, speceminus, specN, 2.20E40, -4.5, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[16]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( speceminus, specO, specOplus, speceminus, speceminus, _kf16(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specOplus, speceminus, speceminus, speceminus, specO, 2.20E40, -4.5, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[17]){
    add_to_dW_fw_2r2p ( specN2, specO2plus,   specN2plus, specO2, 9.90e12, 0.0, 40700.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -1.970e-1, 1.031, -2.049e-1, -4.005, -9.866e-5, 9.90e12, 0.0, 40700.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -1.970e-1, 1.031, -2.049e-1, -4.005, -9.866e-5, 9.90e12, 0.0, 40700.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 4.428e-2, 1.462, 6.450e-1, -4.638, -3.600e-2, 9.90e12, 0.0, 40700.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 4.428e-2, 1.462, 6.450e-1, -4.638, -3.600e-2, 9.90e12, 0.0, 40700.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specN2plus, specO2plus,  specN2,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[18]){
    add_to_dW_fw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.40e13, -1.08, 12800.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -6.300e-2, -5.419e-1, -4.449e-2, -1.266, -1.089e-4, 3.40e13, -1.08, 12800.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -6.300e-2, -5.419e-1, -4.449e-2, -1.266, -1.089e-4, 3.40e13, -1.08, 12800.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 3.518e-2, -7.715e-2, 4.619e-1, -1.923, 9.394e-2, 3.40e13, -1.08, 12800.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 3.518e-2, -7.715e-2, 4.619e-1, -1.923, 9.394e-2, 3.40e13, -1.08, 12800.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specN2, specOplus,  specNOplus, specN,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[19]){
    add_to_dW_fw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.00e12, 0.5, 77200.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -4.092e-1, -7.798e-1, -4.273e-1, -7.618, -4.475e-4, 1.00e12, 0.5, 77200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -4.092e-1, -7.798e-1, -4.273e-1, -7.618, -4.475e-4, 1.00e12, 0.5, 77200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 1.270e-1, 6.372e-1, 1.751, -9.528, -4.012e-2, 1.00e12, 0.5, 77200.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 1.270e-1, 6.372e-1, 1.751, -9.528, -4.012e-2, 1.00e12, 0.5, 77200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specNplus,  specNOplus, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[20]){
    add_to_dW_fw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.40e13, 0.41, 32600.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -3.056e-2, 1.762, -8.580e-2, -3.265, -2.168e-4, 2.40e13, 0.41, 32600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -3.056e-2, 1.762, -8.580e-2, -3.265, -2.168e-4, 2.40e13, 0.41, 32600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 9.436e-2, 3.786, 1.322, -6.192, 7.723e-1, 2.40e13, 0.41, 32600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 9.436e-2, 3.786, 1.322, -6.192, 7.723e-1, 2.40e13, 0.41, 32600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specNO, specO2plus,  specO2, specNOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[21]){
    add_to_dW_fw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.20e13, 0.0, 35500.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -4.862e-2, 1.066, -7.350e-2, -3.537, -8.701e-5, 7.20e13, 0.0, 35500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -4.862e-2, 1.066, -7.350e-2, -3.537, -8.701e-5, 7.20e13, 0.0, 35500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 6.423e-2, 2.396, 9.127e-1, -5.527, 5.453e-1, 7.20e13, 0.0, 35500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 6.423e-2, 2.396, 9.127e-1, -5.527, 5.453e-1, 7.20e13, 0.0, 35500.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO, specN2plus,  specN, specNOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[22]){
    add_to_dW_fw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.70e13, 0.14, 28600.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -2.113e-1, -1.152, -1.758e-1, -2.803, -1.205e-4, 8.70e13, 0.14, 28600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -2.113e-1, -1.152, -1.758e-1, -2.803, -1.205e-4, 8.70e13, 0.14, 28600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 1.523e-2, -1.586, 1.942e-1, -2.103, -4.874e-1, 8.70e13, 0.14, 28600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 1.523e-2, -1.586, 1.942e-1, -2.103, -4.874e-1, 8.70e13, 0.14, 28600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specNplus,  specN, specO2plus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[23]){
    add_to_dW_fw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.40e05, 1.90, 26600.0*R, T, X, dWdT, dWdrhok);
   
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -1.673e-1, -1.965, -1.656e-1, -2.619, -1.102e-4, 1.40e05, 1.90, 26600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -1.673e-1, -1.965, -1.656e-1, -2.619, -1.102e-4, 1.40e05, 1.90, 26600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 1.744e-2, -2.138, 2.354e-1, -2.302, -3.250e-1, 1.40e05, 1.90, 26600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 1.744e-2, -2.138, 2.354e-1, -2.302, -3.250e-1, 1.40e05, 1.90, 26600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specNplus,  specNO, specOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[24]){
    add_to_dW_fw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.20e12, 0.29, 48600.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, -1.978e-1, 3.723e-1, -2.514e-1, -4.815, -3.270e-4, 7.20e12, 0.29, 48600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, -1.978e-1, 3.723e-1, -2.514e-1, -4.815, -3.270e-4, 7.20e12, 0.29, 48600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 1.118e-1, 2.223, 1.557, -7.425, 4.473e-1, 7.20e12, 0.29, 48600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 1.118e-1, 2.223, 1.557, -7.425, 4.473e-1, 7.20e12, 0.29, 48600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specN, specO2plus,  specO, specNOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[25]){
    add_to_dW_fw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.10e11, 0.36, 22800.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial_rodriguez2024(2, 1.438e-2, 1.608, -2.901e-2, -2.271, 2.188e-5, 9.10e11, 0.36, 22800.0*R, Tblim);
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 1.438e-2, 1.608, -2.901e-2, -2.271, 2.188e-5, 9.10e11, 0.36, 22800.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial_rodriguez2024(2, 2.905e-2, 2.473, 4.508e-1, -3.603, 4.514e-1, 9.10e11, 0.36, 22800.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial_rodriguez2024(2, 2.905e-2, 2.473, 4.508e-1, -3.603, 4.514e-1, 9.10e11, 0.36, 22800.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO, specN2plus,  specN2, specOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (REACTION[26]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_fw_2r2p ( specN, specN,   specN2plus, speceminus, 4.40e07, 1.5, 67500.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_2r2p ( specN2plus, speceminus, specN, specN, _kf26b(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }  
  
  if (REACTION[27]){
    add_to_dW_fw_2r2p ( specNplus, specN2,   specN2plus, specN, 1.00e12, 0.5, 12200.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_fw_2r2p ( specN2plus, specN,   specNplus, specN2, 1.00e12, 0.5, 12200.0*R, T, X, dWdT, dWdrhok);
  }  

  if (REACTION[28]){
    add_to_dW_fw_2r2p ( specO2plus, specO,   specOplus, specO2, 4.00e12, -0.09, 18000.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_fw_2r2p ( specOplus, specO2,   specO2plus, specO, 4.00e12, -0.09, 18000.0*R, T, X, dWdT, dWdrhok);
  } 

  if (REACTION[29]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( speceminus, specN2, specN, specN, speceminus, _kf29(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specN, specN, speceminus, speceminus, specN2, 1.20e25, -1.6, 113200.0*R, Te, X, dWdTe, dWdrhok);
  }   

  if (REACTION[30]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( speceminus, specN2, specN2plus, speceminus, speceminus, _kf30(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }

  if (REACTION[31]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( speceminus, specO2, specO2plus, speceminus, speceminus, _kf31(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }

  if (REACTION[32]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( speceminus, specNO, specNOplus, speceminus, speceminus, _kf32(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }

  if (REACTION[33]){
    add_to_dW_fw_3r2p ( specN2plus, speceminus, speceminus, speceminus, specN2, 2.20E40, -4.5, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[34]){
    add_to_dW_fw_3r2p ( specO2plus, speceminus, speceminus, speceminus, specO2, 2.20E40, -4.5, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[35]){
    add_to_dW_fw_3r2p ( specNOplus, speceminus, speceminus, speceminus, specNO, 2.20E40, -4.5, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }  


  for (spec=0; spec<ns; spec++){
    dWdT[spec]+=dWdTTv[spec]*0.5/TTv*Tv;
    dWdTv[spec]+=dWdTTv[spec]*0.5/TTv*T;

    dWdT[spec]+=dWdTTe[spec]*0.5/TTe*Te;
    dWdTe[spec]+=dWdTTe[spec]*0.5/TTe*T;

    dWdTv[spec]+=dWdTvTe[spec]*0.5/TvTe*Te;
    dWdTe[spec]+=dWdTvTe[spec]*0.5/TvTe*Tv;
  }
  
}





void find_Qei_Rodriguez2024( np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){

    if (REACTION[15]) 
      add_to_Qei(gl,Te,specN,_ionizationpot(specN), _kf15(np,gl,Te), rhok, Qei);
    if (REACTION[16]) 
      add_to_Qei(gl,Te,specO,_ionizationpot(specO), _kf16(np,gl,Te), rhok, Qei);
    if (REACTION[30]) 
      add_to_Qei(gl,Te,specN2,_ionizationpot(specN2), _kf30(np,gl,Te), rhok, Qei);
    if (REACTION[31]) 
      add_to_Qei(gl,Te,specO2,_ionizationpot(specO2), _kf31(np,gl,Te), rhok, Qei);
    if (REACTION[32]) 
      add_to_Qei(gl,Te,specNO,_ionizationpot(specNO), _kf32(np,gl,Te), rhok, Qei);

 
}



void find_dQei_dx_Rodriguez2024( np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){

    if (REACTION[15]) 
      add_to_dQei(gl,Te,specN,_ionizationpot(specN), _kf15(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[16]) 
      add_to_dQei(gl,Te,specO,_ionizationpot(specO), _kf16(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[30])           
      add_to_dQei(gl,Te,specN2,_ionizationpot(specN2), _kf30(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[31]) 
      add_to_dQei(gl,Te,specO2,_ionizationpot(specO2), _kf31(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[32]) 
      add_to_dQei(gl,Te,specNO,_ionizationpot(specNO), _kf32(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);

}
