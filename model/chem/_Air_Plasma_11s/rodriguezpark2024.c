// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2021 Bernard Parent
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

#define EXPON_TTV 1.0
#define EXPON_TTE 0.0
#define EXPON_TVTE 0.0

/* set all reactions to true except for testing purposes */
const static bool REACTION[31]=
  {
   TRUE, /* reaction 0 */
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
   TRUE, /* reaction 18 */
   TRUE, /* reaction 19 */ 
   TRUE, /* reaction 20 */
   TRUE, /* reaction 21 */
   TRUE, /* reaction 22 */
   TRUE, /* reaction 23 */
   TRUE, /* reaction 24 */
   TRUE, /* reaction 25 */
   TRUE, /* reaction 26 */
   TRUE, /* reaction 27 */
   TRUE, /* reaction 28 */
   TRUE, /* reaction 29 */
   TRUE, /* reaction 30 */
  };

#define specEND -1

const static long specM1[]=
  {
   specN, specO, specNplus, specOplus, specEND
  };

const static long specM2[]=
  {
   specN2, specO2, specNO, specN2plus, specO2plus, specNOplus, specEND
  };

const static long specM3[]=
  {
   specN, specO, specNO, specNplus, specOplus,  specEND
  };

const static long specM4[]=
  {
   specN2, specO2, specN2plus, specO2plus, specNOplus,  specEND
  };

static double _kf3(np_t np, gl_t *gl, double Te){
  double kf3;
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
  kf3 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf3 ));
}  

static double _kf12b(np_t np, gl_t *gl, double Te){
  double kf12b;
  double Te_control[] = 
  { 
    2.96134635039148,
    5.71525385948086,
    9.17412056693579,
    12.3695988451580,
    14.9141228466324
  };
  double kf_control[] = 
  {
    -15.1096633905206,
    -15.6250837536637,
    -16.9584559651989,
    -19.0191856040400,
    -19.9348084765821
  };  
  int N = sizeof(Te_control)/sizeof(Te_control[0]);
  Te=min(Te_control[N-1],max(log( Te ),Te_control[0]));
  kf12b = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf12b ));
}

static double _kf26(np_t np, gl_t *gl, double Te){
  /* O ionization */
  double kf26;
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
  kf26 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf26 ));
}

static double _kf27(np_t np, gl_t *gl, double Te){
  /* N ionization */
  double kf27;
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
  kf27 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf27 ));
}

static double _kf28(np_t np, gl_t *gl, double Te){
  /* O2 ionization */
  double kf28;
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
  kf28 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf28 ));
}

static double _kf29(np_t np, gl_t *gl, double Te){
  /* N2 ionization */
  double kf29;
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
  kf29 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf29 ));
}

static double _kf30(np_t np, gl_t *gl, double Te){
  /* NO ionization */
  double kf30;
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
  kf30 = EXM_f_from_monotonespline(N, Te_control, kf_control, Te);
  return(exp( kf30 ));
}


void find_W_RodriguezPark2024 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                       double Qbeam, spec_t W ) {
  double N[ns];
  double R,kf,TTv,TvTe,TTe;
  long k;
  spec_t X;

  /* find properties needed by add_to_W* functions */
  R=Rchem;
  TTv=pow(T,EXPON_TTV)*pow(Tv,1.0-EXPON_TTV);
  TTe=pow(T,EXPON_TTE)*pow(Te,1.0-EXPON_TTE);
  TvTe=pow(Tv,EXPON_TVTE)*pow(Te,1.0-EXPON_TVTE);

  
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k] = 0.0;
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
  }


  if (REACTION[1]) {
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, T, X, W );
    }
  }
      

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, T, X, W );
    }
  }

  if (REACTION[3]){
    add_to_W_2r3p ( specN2, speceminus, specN, specN, speceminus, _kf3(np,gl,Te), N, W );
    add_to_W_bw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 3.0e24, -1.6, 113200.0*R, TTe, X, W );
  }

  if (REACTION[4]){
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, T, X, W );
    }
  }

  if (REACTION[5]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, T, X, W );
    }
  }

  if (REACTION[6]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, T, X, W );
    }
  }

  if (REACTION[7]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, TTv, X, W );
      add_to_W_bw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, T, X, W );
    }
  }

  if (REACTION[8]){
    add_to_W_fwbw_2r2p ( specNO, specO,   specN, specO2, 8.4e12, 0.0, 19400.0*R, T, X, W );
  }

  if (REACTION[9]){
    add_to_W_fwbw_2r2p ( specN2, specO,   specNO, specN, 5.7e12, 0.42, 42938.0*R, T, X, W );
  }

  if (REACTION[10]){
    add_to_W_fw_2r2p ( specN, specO,   specNOplus, speceminus, 5.3e12, 0.0, 32000.0*R, T, X, W );
    add_to_W_fw_2r2p ( specNOplus, speceminus, specN, specO, 3.00E-7 * calA * pow ( 300, 0.56 ), -0.56, 0.0 * R, Te, X, W );
  }

  if (REACTION[11]){
    add_to_W_fw_2r2p ( specO, specO,   specO2plus, speceminus, 1.1e13, 0.0, 81200.0*R, T, X, W );
    add_to_W_fw_2r2p ( specO2plus, speceminus, specO, specO, 2.40E-7 * calA * pow ( 300, 0.70 ), -0.70, 0.0 * R, Te, X, W );
  }

  if (REACTION[12]){
    add_to_W_fw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, T, X, W );
    add_to_W_2r2p ( specN2plus, speceminus, specN, specN, _kf12b(np,gl,Te), N, W );
  }

  if (REACTION[13]){
    add_to_W_fwbw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.0e12, 0.5, 77200.0*R, T, X, W );
  }

  if (REACTION[14]){
    add_to_W_fwbw_2r2p ( specNplus, specN2,   specN2plus, specN, 1.0e12, 0.5, 12200.0*R, T, X, W );
  }

  if (REACTION[15]){
    add_to_W_fwbw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.7e13, 0.14, 28600.0*R, T, X, W );
  }

  if (REACTION[16]){
    add_to_W_fwbw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.4e5, 1.90, 26600.0*R, T, X, W );
  }

  if (REACTION[17]){
    add_to_W_fwbw_2r2p ( specO2plus, specN2,   specN2plus, specO2, 9.9e12, 0.00, 40700.0*R, T, X, W );
  }

  if (REACTION[18]){
    add_to_W_fwbw_2r2p ( specO2plus, specO,   specOplus, specO2, 4.0e12, -0.09, 18000.0*R, T, X, W );
  }

  if (REACTION[19]){
    add_to_W_fwbw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.4e13, -1.08, 12800.0*R, T, X, W );
  }

  if (REACTION[20]){
    add_to_W_fwbw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.4e13, 0.41, 32600.0*R, T, X, W );
  }

  if (REACTION[21]){
    add_to_W_fwbw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.2e12, 0.29, 48600.0*R, T, X, W );
  }

  if (REACTION[22]){
    add_to_W_fwbw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.1e11, 0.36, 22800.0*R, T, X, W );
  }

  if (REACTION[23]){
    add_to_W_fwbw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.2e13, 0.00, 35500.0*R, T, X, W );
  }

  if (REACTION[24]){
    kf=_kf_Arrhenius(2, 1.07e11, -0.52, 0.0, Te);
    add_to_W_2r1p ( specOplus, speceminus,   specO,   kf , N, W);
  }

  if (REACTION[25]){
    kf=_kf_Arrhenius(2, 1.52e11, -0.48, 0.0, Te);
    add_to_W_2r1p ( specNplus, speceminus,   specN,   kf , N, W);
  }

  if (REACTION[26]) {
    add_to_W_2r3p ( specO, speceminus, specOplus, speceminus, speceminus, _kf26(np,gl,Te), N, W );
    add_to_W_fw_3r2p ( specOplus, speceminus, speceminus, specO, speceminus, 2.2E40, -4.5, 0.0, Te, X, W );
  }

  if (REACTION[27]){
    add_to_W_2r3p ( specN, speceminus, specNplus, speceminus, speceminus, _kf27(np,gl,Te), N, W );
    add_to_W_fw_3r2p ( specNplus, speceminus, speceminus, specN, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, W );
  }

  if (REACTION[28]){
    add_to_W_2r3p ( specO2, speceminus, specO2plus, speceminus, speceminus, _kf28(np,gl,Te), N, W );
    add_to_W_fw_3r2p ( specO2plus, speceminus, speceminus, specO2, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, W );
  }

  if (REACTION[29]){
    add_to_W_2r3p ( specN2, speceminus, specN2plus, speceminus, speceminus, _kf29(np,gl,Te), N, W ); 
    add_to_W_fw_3r2p ( specN2plus, speceminus, speceminus, specN2, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, W );
  }

  if (REACTION[30]){
    add_to_W_2r3p ( specNO, speceminus, specNOplus, speceminus, speceminus, _kf30(np,gl,Te), N, W );
    add_to_W_fw_3r2p ( specNOplus, speceminus, speceminus, specNO, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, W );
  }


}



/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx_RodriguezPark2024 ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s, spec;                    
  spec_t N;
  double TTv,TTe,TvTe,R,kf,dkfdTe,dkfdT,dkfdTv;
  spec_t dWdTTv,dWdTTe,dWdTvTe;
  spec_t X;

  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  
  R=Rchem;
  TTv=pow(T,EXPON_TTV)*pow(Tv,1.0-EXPON_TTV);
  TTe=pow(T,EXPON_TTE)*pow(Te,1.0-EXPON_TTE);
  TvTe=pow(Tv,EXPON_TVTE)*pow(Te,1.0-EXPON_TVTE);
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
      add_to_dW_fw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specN2, specM1[k],   specN, specN, specM1[k], 3.0e22, -1.6, 113200.0*R, T, X, dWdT, dWdrhok );
    }
  }
      

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[3]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specN2, speceminus, specN, specN, speceminus, _kf3(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_bw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 3.0e24, -1.6, 113200.0*R, TTe, X, dWdTTe, dWdrhok );
  }

  if (REACTION[4]){
    for (k=0; specM1[k]!=specEND; k++) {
      add_to_dW_fw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specO2, specM1[k],   specO, specO, specM1[k], 1.0e22, -1.5, 59500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[5]){
    for (k=0; specM2[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specO2, specM2[k],   specO, specO, specM2[k], 2.0e21, -1.5, 59500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[6]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specNO, specM3[k],   specN, specO, specM3[k], 1.1e17, 0.0, 75500.0*R, T, X, dWdT, dWdrhok );
    }
  }

  if (REACTION[7]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_dW_fw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, TTv, X, dWdTTv, dWdrhok );
      add_to_dW_bw_2r3p ( specNO, specM4[k],   specN, specO, specM4[k], 5.0e15, 0.0, 75500.0*R, T, X, dWdT, dWdrhok );
    }
  }



  if (REACTION[8]){
    add_to_dW_fwbw_2r2p ( specNO, specO,   specN, specO2, 8.4e12, 0.0, 19400.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[9]){
    add_to_dW_fwbw_2r2p ( specN2, specO,   specNO, specN, 5.7e12, 0.42, 42938.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[10]){
    add_to_dW_fw_2r2p ( specN, specO,   specNOplus, speceminus, 5.3e12, 0.0, 32000.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_fw_2r2p ( specNOplus, speceminus, specN, specO, 3.00E-7 * calA * pow ( 300, 0.56 ), -0.56, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[11]){
    add_to_dW_fw_2r2p ( specO, specO,   specO2plus, speceminus, 1.1e13, 0.0, 81200.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_fw_2r2p ( specO2plus, speceminus, specO, specO, 2.40E-7 * calA * pow ( 300, 0.70 ), -0.70, 0.0 * R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[12]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_fw_2r2p ( specN, specN,   specN2plus, speceminus, 2.0e13, 0.0, 67700.0*R, T, X, dWdT, dWdrhok);
    add_to_dW_2r2p ( specN2plus, speceminus, specN, specN, _kf12b(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
  }


  if (REACTION[13]){
    add_to_dW_fwbw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.0e12, 0.5, 77200.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[14]){
    add_to_dW_fwbw_2r2p ( specNplus, specN2,   specN2plus, specN, 1.0e12, 0.5, 12200.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[15]){
    add_to_dW_fwbw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.7e13, 0.14, 28600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[16]){
    add_to_dW_fwbw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.4e5, 1.90, 26600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[17]){
    add_to_dW_fwbw_2r2p ( specO2plus, specN2,   specN2plus, specO2, 9.9e12, 0.00, 40700.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[18]){
    add_to_dW_fwbw_2r2p ( specO2plus, specO,   specOplus, specO2, 4.0e12, -0.09, 18000.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[19]){
    add_to_dW_fwbw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.4e13, -1.08, 12800.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[20]){
    add_to_dW_fwbw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.4e13, 0.41, 32600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[21]){
    add_to_dW_fwbw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.2e12, 0.29, 48600.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[22]){
    add_to_dW_fwbw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.1e11, 0.36, 22800.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[23]){
    add_to_dW_fwbw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.2e13, 0.00, 35500.0*R, T, X, dWdT, dWdrhok);
  }

  if (REACTION[24]){
    kf=_kf_Arrhenius(2, 1.07e11, -0.52, 0.0, Te);
    dkfdT=0.0;
    dkfdTv=0.0;
    dkfdTe=_dkfdT_Arrhenius(2, 1.07e11, -0.52, 0.0, Te);
    add_to_dW_2r1p ( specOplus, speceminus,   specO,   kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[25]){
    kf=_kf_Arrhenius(2, 1.52e11, -0.48, 0.0, Te);
    dkfdT=0.0;
    dkfdTv=0.0;
    dkfdTe=_dkfdT_Arrhenius(2, 1.52e11, -0.48, 0.0, Te);
    add_to_dW_2r1p ( specNplus, speceminus,   specN,   kf , N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[26]) {
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specO, speceminus, specOplus, speceminus, speceminus, _kf26(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specOplus, speceminus, speceminus, specO, speceminus, 2.2E40, -4.5, 0.0, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[27]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specN, speceminus, specNplus, speceminus, speceminus, _kf27(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specNplus, speceminus, speceminus, specN, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[28]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specO2, speceminus, specO2plus, speceminus, speceminus, _kf28(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specO2plus, speceminus, speceminus, specO2, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[29]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specN2, speceminus, specN2plus, speceminus, speceminus, _kf29(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specN2plus, speceminus, speceminus, specN2, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (REACTION[30]){
    dkfdTe = 0.0;
    dkfdT = 0.0;
    dkfdTv = 0.0;
    add_to_dW_2r3p ( specNO, speceminus, specNOplus, speceminus, speceminus, _kf30(np,gl,Te), N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe );
    add_to_dW_fw_3r2p ( specNOplus, speceminus, speceminus, specNO, speceminus, 2.2E40, -4.50, 0.0*R, Te, X, dWdTe, dWdrhok );
  }



  for (spec=0; spec<ns; spec++){
    dWdT[spec]+=dWdTTv[spec]*EXPON_TTV*TTv/T;
    dWdTv[spec]+=dWdTTv[spec]*(1.0-EXPON_TTV)*TTv/Tv;

    dWdT[spec]+=dWdTTe[spec]*EXPON_TTE*TTe/T;
    dWdTe[spec]+=dWdTTe[spec]*(1.0-EXPON_TTE)*TTe/Te;

    dWdTv[spec]+=dWdTvTe[spec]*EXPON_TVTE*TvTe/Tv;
    dWdTe[spec]+=dWdTvTe[spec]*(1.0-EXPON_TVTE)*TvTe/Te;

  }
  
}



void find_Qei_RodriguezPark2024( np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){

    if (REACTION[26]) 
      add_to_Qei(gl,Te,specO,_ionizationpot(specO), _kf26(np,gl,Te), rhok, Qei);
    if (REACTION[27]) 
      add_to_Qei(gl,Te,specN,_ionizationpot(specN), _kf27(np,gl,Te), rhok, Qei);
    if (REACTION[28]) 
      add_to_Qei(gl,Te,specO2,_ionizationpot(specO2), _kf28(np,gl,Te), rhok, Qei);
    if (REACTION[29]) 
      add_to_Qei(gl,Te,specN2,_ionizationpot(specN2), _kf29(np,gl,Te), rhok, Qei);
    if (REACTION[30]) 
      add_to_Qei(gl,Te,specNO,_ionizationpot(specNO), _kf30(np,gl,Te), rhok, Qei);
 
}



void find_dQei_dx_RodriguezPark2024( np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){

    if (REACTION[26]) 
      add_to_dQei(gl,Te,specO,_ionizationpot(specO), _kf26(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[27]) 
      add_to_dQei(gl,Te,specN,_ionizationpot(specN), _kf27(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[28]) 
      add_to_dQei(gl,Te,specO2,_ionizationpot(specO2), _kf28(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[29]) 
      add_to_dQei(gl,Te,specN2,_ionizationpot(specN2), _kf29(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[30]) 
      add_to_dQei(gl,Te,specNO,_ionizationpot(specNO), _kf30(np,gl,Te), 0.0, rhok, dQeidrhok, dQeidTe);

}
