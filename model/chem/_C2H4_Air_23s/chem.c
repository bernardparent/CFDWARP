// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Felipe Martin Rodriguez Fuentes

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

/* taken from  Zettervall, N., Fureby, C., and Nilsson, E. J., 
  “Small skeletal kinetic reaction mechanism for ethylene–air combustion,” 
   Energy and Fuels, Vol. 31, No. 12, 2017, pp. 14138–14149. */
   
const static bool REACTION[67]=
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
   TRUE,  /* reaction 25 */
   TRUE,  /* reaction 26 */
   TRUE, /* reaction 27 */
   TRUE, /* reaction 28 */
   TRUE, /* reaction 29 */
   TRUE, /* reaction 30 */
   TRUE, /* reaction 31 */
   TRUE, /* reaction 32 */
   TRUE, /* reaction 33 */
   TRUE, /* reaction 34 */
   TRUE, /* reaction 35 */
   TRUE, /* reaction 36 */
   TRUE, /* reaction 37 */
   TRUE, /* reaction 38 */
   TRUE, /* reaction 39 */
   TRUE, /* reaction 40 */
   TRUE, /* reaction 41 */
   TRUE, /* reaction 42 */
   TRUE, /* reaction 43 */
   TRUE, /* reaction 44 */
   TRUE, /* reaction 45 */
   TRUE, /* reaction 46 */
   TRUE, /* reaction 47 */
   TRUE, /* reaction 48 */
   TRUE, /* reaction 49 */
   TRUE, /* reaction 50 */
   TRUE, /* reaction 51 */
   TRUE, /* reaction 52 */
   TRUE, /* reaction 53 */
   TRUE, /* reaction 54 */
   TRUE, /* reaction 55 */
   TRUE, /* reaction 56 */
   TRUE, /* reaction 57 */
   TRUE, /* reaction 58 */
   TRUE, /* reaction 59 */
   TRUE, /* reaction 60 */
   TRUE, /* reaction 61 */
   TRUE, /* reaction 62 */
   TRUE, /* reaction 63 */
   TRUE, /* reaction 64 */
   TRUE, /* reaction 65 */
   TRUE, /* reaction 66 */
  };   
   
   
#define specEND -1

   const static long specM[]=
  {
   specH, specH2, specO, specO2, specOH, specH2O, specN2, specH2O2,
   specC2H, specC2H2, specC2H3, specC2H4, specC2H5, specCH, specCH2, 
   specCH3, specCH4, specHO2, specCHO, specCH2O, specCH3O, specCO,  
   specCO2, specEND
  };
  
  
  /* Third-body efficiencies */
  
  /* Smooke, M. D.; Giovangigli, V. Formulation of the Premixed
     and Nonpremixed Test Problems. In Reduced Chemical Mechansims
     and Asymptotic Approximations for Methane-Air Flames; Springer-
     Verlag: New York, 1991; p 384. */
	 
  const static double eff1[]=                       
  {
   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
   1.0
  };
  
  const static double eff2[]=                       
  {
   1.0, 1.0, 1.0, 0.4, 1.0, 6.5, 0.4, 1.0,
   1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 
   1.0, 6.5, 1.0, 1.0, 1.0, 1.0, 0.75,  
   1.5
  };
  
  const static double eff3[]=                       
  {
   1.0, 1.0, 1.0, 0.4, 1.0, 6.5, 0.4, 1.0,
   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
   1.0, 6.5, 1.0, 1.0, 1.0, 1.0, 0.75,  
   1.5
  };  
  

void write_model_chem_template(FILE **controlfile){
}

void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
}

/*  Third-body recombinations/dissociations: Reaction no. (2), (15), (21), (25), (26), (35), (45), (56), (61), (62), (65), (66) */
/* "Fall-off" reactions (+M): Reaction no. (25), (26) Lindemann approx*/
/*  Units of s, mole, cm3, cal, and K */

void find_W ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
 
  long k;
  spec_t X;

  for ( k = 0; k < ns; k++ ) {
    X[k]   = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
    W[k]   = 0.0;
  }
  

  if (REACTION[1]){
  // reaction (1): C2H5 + H -> CH3 + CH3
  add_to_W_fw_2r2p ( specC2H5, specH, specCH3, specCH3, 1.50e+13, 0.0, 0.0, T, X, W );
  }
  
  if (REACTION[2]){
	 for (k=0; specM[k]!=specEND; k++){
  // reaction (2): C2H4 + H + M-> C2H5 + M
  add_to_W_fw_3r2p ( specC2H4, specH, specM[k], specC2H5, specM[k], eff1[k]*4.17e+10, 0.0, 11030.0, T, X, W );
	 }
  } 
  
  if (REACTION[3]){
  // reaction (3): C2H4 + H -> C2H5
  add_to_W_fw_2r1p ( specC2H4, specH, specC2H5,       1.80e+12, 1.0, 9400.0, T, X, W ); 
  } 
  
  if (REACTION[4]){
  // reaction (4): C2H4 -> C2H3 + H
  add_to_W_fw_1r2p ( specC2H4, specC2H3, specH,       1.00e+12, 0.0, 48200.0, T, X, W );
  } 
  
  if (REACTION[5]){
  // reaction (5): C2H4 + H -> C2H3 + H2
  add_to_W_fw_2r2p ( specC2H4, specH, specC2H3, specH2,  5.00e+12, 0.0, 5000.0, T, X, W );
  } 
  
  if (REACTION[6]){
  // reaction (6): C2H4 + O -> CHO + CH3
  add_to_W_fw_2r2p ( specC2H4, specO, specCHO, specCH3,  3.31e+12, 0.0, 1130.0, T, X, W );
  } 
  
  if (REACTION[7]){
  // reaction (7): CHO + CH3 -> C2H4 + O
  add_to_W_fw_2r2p ( specCHO, specCH3, specC2H4, specO,  1.58e+11, 0.0, 31180.0, T, X, W );
  } 
  
  if (REACTION[8]){
  // reaction (8): C2H4 + OH -> C2H3 + H2O
  add_to_W_fw_2r2p ( specC2H4, specOH, specC2H3, specH2O,  4.79e+12, 0.0, 1230.0, T, X, W );
  } 
  
  if (REACTION[9]){
  // reaction (9): C2H3 + H2O -> C2H4 + OH
  add_to_W_fw_2r2p ( specC2H3, specH2O, specC2H4, specOH,  1.20e+12, 0.0, 14000.0, T, X, W );
  } 
  
  if (REACTION[10]){  
  // reaction (10): C2H4 + CH3 -> C2H3 + CH4
  add_to_W_fw_2r2p ( specC2H4, specCH3, specC2H3, specCH4,  1.00e+13, 0.0, 13000.0, T, X, W );
  } 
  
  if (REACTION[11]){
  // reaction (11): C2H3 + CH4 -> C2H4 + CH3
  add_to_W_fw_2r2p ( specC2H3, specCH4, specC2H4, specCH3,  3.02e+13, 0.0, 12580.0, T, X, W );
  } 
  ;
  if (REACTION[12]){
  // reaction (12): C2H3 -> C2H + H2
  add_to_W_fw_1r2p ( specC2H3, specC2H, specH2,  1.80e+12, 0.0, 44800.0, T, X, W );
  } 
  
  if (REACTION[13]){
  // reaction (13): C2H3 + H -> C2H2 + H2
  add_to_W_fw_2r2p ( specC2H3, specH, specC2H2, specH2,  2.00e+13, 0.0, 2500.0, T, X, W );
  } 
  
  if (REACTION[14]){
  // reaction (14): C2H3 + OH -> C2H2 + H2O
  add_to_W_fw_2r2p ( specC2H3, specOH, specC2H2, specH2O,  7.00e+14, 0.0, 1100.0, T, X, W );
  } 
  
  if (REACTION[15]){
	 for (k=0; specM[k]!=specEND; k++){	  
  // reaction (15): C2H2 + H + M -> C2H3 + M
  add_to_W_fw_3r2p ( specC2H2, specH, specM[k], specC2H3, specM[k], eff1[k]*1.23e+11, 1.0, 10360.0, T, X, W );
     }
  } 
  
  if (REACTION[16]){
  // reaction (16): C2H2 + H -> C2H + H2
  add_to_W_fw_2r2p ( specC2H2, specH, specC2H, specH2,  2.00e+14, 0.0, 19000.0, T, X, W );
  } 
  
  if (REACTION[17]){  
  // reaction (17): C2H2 + OH -> C2H + H2O
  add_to_W_fw_2r2p ( specC2H2, specOH, specC2H, specH2O,  8.00e+12, 0.0, 5000.0, T, X, W );  
  } 
  
  if (REACTION[18]){  
  // reaction (18): C2H + H2O -> C2H2 + OH
  add_to_W_fw_2r2p ( specC2H, specH2O, specC2H2, specOH,  5.37e+12, 0.0, 16360.0, T, X, W );    
  } 
  
  if (REACTION[19]){  
  // reaction (19): C2H2 + O -> C2H + OH
  add_to_W_fw_2r2p ( specC2H2, specO, specC2H, specOH,  3.24e+15, 0.6, 12000.0, T, X, W );   
  } 
   
  if (REACTION[20]){ 
  // reaction (20): C2H + OH -> C2H2 + O
  add_to_W_fw_2r2p ( specC2H, specOH, specC2H2, specO,  2.95e+14, 0.6, 910.0, T, X, W );   
  } 
  
  if (REACTION[21]){
	 for (k=0; specM[k]!=specEND; k++){		  
  // reaction (21): C2H + H + M -> C2H2 + M
    add_to_W_fw_3r2p ( specC2H, specH, specM[k], specC2H2, specM[k], eff1[k]*1.10e+9, 1.0, 770.0, T, X, W );   
   }  
  } 
  
  if (REACTION[22]){
  // reaction (22): C2H + O2  -> CHO + CO
  add_to_W_fw_2r2p ( specC2H, specO2, specCHO, specCO,  1.00e+13, 0.0, 7000.0, T, X, W );  
  } 
  
  if (REACTION[23]){
  // reaction (23): CHO + CO  -> C2H + O2
  add_to_W_fw_2r2p ( specCHO, specCO, specC2H, specO2,  8.51e+12, 0.0, 138400.0, T, X, W ); 
  } 
  
  if (REACTION[24]){
  // reaction (24): C2H + OH  -> CHO + CH
  add_to_W_fw_2r2p ( specC2H, specOH, specCHO, specCH,  8.00e+12, 0.0, 5000.0, T, X, W ); 
  } 
  
  if (REACTION[25]){
  // reaction (25): CH4(+M) -> CH3 + H(+M) 	  	  
	 for (k=0; specM[k]!=specEND; k++){	 
    add_to_W_fw_2r3p_Lindemann ( specCH4, specM[k], specCH3, specH, specM[k], eff2[k]*6.30e+14, 0.0, 104000.0, eff2[k]*1.00e+17, 0.0, 86000.0, T, X, W );
	 }
  } 
  
  if (REACTION[26]){
  // reaction (26): CH3 + H(+M)  -> CH4(+M) 
	 for (k=0; specM[k]!=specEND; k++){
    add_to_W_fw_3r2p_Lindemann ( specCH3, specH, specM[k], specCH4, specM[k], eff2[k]*5.20e+12, 0.0, -1.310e+3, eff2[k]*8.25e+14, 0.0, -1.9310e+4,  T, X, W );
	 }  
  } 
  
  if (REACTION[27]){
  // reaction (27): CH4 + H  -> CH3 + H2
  add_to_W_fw_2r2p ( specCH4, specH, specCH3, specH2,  2.20e+4, 3.0, 8750.0, T, X, W );   
  } 
  
  if (REACTION[28]){
  // reaction (28): CH3 + H2  -> CH4 + H
  add_to_W_fw_2r2p ( specCH3, specH2, specCH4, specH,  9.57e+2, 3.0, 8750.0, T, X, W );   
  } 
  
  if (REACTION[29]){  
  // reaction (29): CH4 + OH  -> CH3 + H2O
  add_to_W_fw_2r2p ( specCH4, specOH, specCH3, specH2O,  1.60e+6, 2.1, 2460.0, T, X, W ); 
  } 
  
  if (REACTION[30]){  
  // reaction (30): CH3 + H2O  -> CH4 + OH
  add_to_W_fw_2r2p ( specCH3, specH2O, specCH4, specOH,  3.02e+5, 2.1, 17422.0, T, X, W );  
  } 
  
  if (REACTION[31]){
  // reaction (31) CH3 + O  -> CH2O + H
  add_to_W_fw_2r2p ( specCH3, specO, specCH2O, specH,  6.80e+13, 0.0, 0.0, T, X, W ); 
  } 
  
  if (REACTION[32]){
  // reaction (32) CH3 + O2  -> CH3O + O
  add_to_W_fw_2r2p ( specCH3, specO2, specCH3O, specO,  3.00e+13, 0.0, 25652.0, T, X, W );     
  } 
  
  if (REACTION[33]){  
  // reaction (33) CH3 + OH  -> CH2 + H2O
  add_to_W_fw_2r2p ( specCH3, specOH, specCH2, specH2O,  7.60e+6, 2.0, 5000.0, T, X, W );  
  } 
    
  if (REACTION[34]){  
  // reaction (34) CH3O + H  -> CH2O + H2
  add_to_W_fw_2r2p ( specCH3O, specH, specCH2O, specH2,  2.00e+13, 0.0, 0.0, T, X, W );   
  } 
  
  if (REACTION[35]){ 
	 for (k=0; specM[k]!=specEND; k++){  
  // reaction (35) CH3O + M  -> CH2O + H + M
    add_to_W_fw_2r3p ( specCH3O, specM[k], specCH2O, specH, specM[k] ,  eff3[k]*2.40e+13, 0.0, 28812.0, T, X, W );
   }  
  } 
  
  if (REACTION[36]){  
  // reaction (36) CH2 + O  -> CO + H2
  add_to_W_fw_2r2p ( specCH2, specO, specCO, specH2,  3.00e+13, 0.0, 0.0, T, X, W );   
  } 
  
  if (REACTION[37]){  
  // reaction (37) CH2 + OH  -> CH + H2O
  add_to_W_fw_2r2p ( specCH2, specOH, specCH, specH2O,  1.13e+7, 2.0, 3000.0, T, X, W );   
  } 
  
  if (REACTION[38]){  
  // reaction (38) CH2O + H  -> CHO + H2
  add_to_W_fw_2r2p ( specCH2O, specH, specCHO, specH2,  5.00e+13, 0.0, 3991.0, T, X, W );   
  } 
  
  if (REACTION[39]){  
  // reaction (39) CH2O + OH  -> CHO + H2O
  add_to_W_fw_2r2p ( specCH2O, specOH, specCHO, specH2O,  1.40e+14, 0.0, 1100.0, T, X, W );  
  } 
  
  if (REACTION[40]){  
  // reaction (40) CH + O  -> CO + H
  add_to_W_fw_2r2p ( specCH, specO, specCO, specH,  5.70e+13, 0.0, 0.0, T, X, W );  
  } 
  
  if (REACTION[41]){  
  // reaction (41) CH + OH  -> CHO + H
  add_to_W_fw_2r2p ( specCH, specOH, specCHO, specH,  3.00e+13, 0.0, 0.0, T, X, W );    
  } 
    
  if (REACTION[42]){  
  // reaction (42) CH + O2  -> CHO + O
  add_to_W_fw_2r2p ( specCH, specO2, specCHO, specO,  3.30e+13, 0.0, 0.0, T, X, W );  
  } 
  
  if (REACTION[43]){  
  // reaction (43) CH + CO2  -> CHO + CO
  add_to_W_fw_2r2p ( specCH, specCO2, specCHO, specCO,  8.40e+13, 0.0, 200.0, T, X, W );   
  } 
  
  if (REACTION[44]){  
  // reaction (44) CHO + H  -> CO + H2
  add_to_W_fw_2r2p ( specCHO, specH, specCO, specH2,  4.00e+13, 0.0, 0.0, T, X, W );   
  } 
  
  if (REACTION[45]){ 
	 for (k=0; specM[k]!=specEND; k++){  
  // reaction (45) CHO + M  -> CO + H + M
    add_to_W_fw_2r3p ( specCHO, specM[k], specCO, specH, specM[k], eff3[k]*1.60e+14, 0.0, 14700.0, T, X, W );
   }  
  } 
  
  if (REACTION[46]){  
  // reaction (46) CO + OH  -> CO2 + H
  add_to_W_fw_2r2p ( specCO, specOH, specCO2, specH,  1.51e+7, 1.3, -758.0, T, X, W );  
  } 
  
  if (REACTION[47]){  
  // reaction (47) CO2 + H  -> CO + OH
  add_to_W_fw_2r2p ( specCO2, specH, specCO, specOH,  1.57e+9, 1.3, 19800.0, T, X, W );  
  } 
  
  if (REACTION[48]){  
  // reaction (48) H + O2  -> OH + O
  add_to_W_fw_2r2p ( specH, specO2, specOH, specO,  5.00e+14, 0.0, 16800.0, T, X, W );  
  } 
  
  if (REACTION[49]){  
  // reaction (49) OH + O  -> H + O2
  add_to_W_fw_2r2p ( specOH, specO, specH, specO2,  1.20e+13, 0.0, 690.0, T, X, W );   
  } 
  
  if (REACTION[50]){  
  // reaction (50) O + H2 -> OH + H
  add_to_W_fw_2r2p ( specO, specH2, specOH, specH,  1.80e+10, 1.0, 8826.0, T, X, W );   
  } 
  
  if (REACTION[51]){  
  // reaction (51) OH + H -> O + H2
  add_to_W_fw_2r2p ( specOH, specH, specO, specH2,  8.00e+9, 1.0, 6760.0, T, X, W );   
  } 
  
  if (REACTION[52]){  
  // reaction (52) H2 + OH -> H2O + H
  add_to_W_fw_2r2p ( specH2, specOH, specH2O, specH,  1.17e+9, 1.3, 3626.0, T, X, W );  
  } 
  
  if (REACTION[53]){  
  // reaction (53) H2O + H -> H2 + OH
  add_to_W_fw_2r2p ( specH2O, specH, specH2, specOH,  6.00e+9, 1.3, 18588.0, T, X, W );   
  } 
  
  if (REACTION[54]){  
  // reaction (54) OH + OH -> O + H2O
  add_to_W_fw_2r2p ( specOH, specOH, specO, specH2O,  6.00e+8, 1.3, 0.0, T, X, W );   
  } 
  
  if (REACTION[55]){  
  // reaction (55) O + H2O -> OH + OH
  add_to_W_fw_2r2p ( specO, specH2O, specOH, specOH,  4.00e+9, 1.3, 17029.0, T, X, W ); 
  } 
  
  if (REACTION[56]){
	 for (k=0; specM[k]!=specEND; k++){	  
  // reaction (56) H + O2 + M -> HO2 + M
    add_to_W_fw_3r2p ( specH, specO2, specM[k], specHO2, specM[k] ,  eff3[k]*2.50e+18, -0.8, 0.0, T, X, W ); 
   }  
  } 
  
  if (REACTION[57]){  
  // reaction (57) H + HO2 -> OH + OH
  add_to_W_fw_2r2p ( specH, specHO2, specOH, specOH,  1.50e+14, 0.0, 1004.0, T, X, W );  
  } 
  
  if (REACTION[58]){  
  // reaction (58) H + HO2 -> H2 + O2
  add_to_W_fw_2r2p ( specH, specHO2, specH2, specO2,  2.50e+13, 0.0, 700.0, T, X, W );   
  } 
  
  if (REACTION[59]){  
  // reaction (59) OH + HO2 -> H2O + O2
  add_to_W_fw_2r2p ( specOH, specHO2, specH2O, specO2,  2.00e+13, 0.0, 1000.0, T, X, W ); 
  } 
  
  if (REACTION[60]){  
  // reaction (60) HO2 + HO2 -> H2O2 + O2
  add_to_W_fw_2r2p ( specHO2, specHO2, specH2O2, specO2,  2.00e+14, 0.0, 0.0, T, X, W );   
  } 
  
  if (REACTION[61]){
	 for (k=0; specM[k]!=specEND; k++){	 	  
  // reaction (61) H2O2 + M -> OH + OH + M
    add_to_W_fw_2r3p ( specH2O2, specM[k], specOH, specOH, specM[k],  eff3[k]*1.30e+17, 0.0, 45500.0, T, X, W );    
   }
  } 
  
  if (REACTION[62]){
	 for (k=0; specM[k]!=specEND; k++){		  
  // reaction (62) OH + OH + M -> H2O2 + M 
  add_to_W_fw_3r2p ( specOH, specOH, specM[k], specH2O2, specM[k],  eff3[k]*9.86e+14, 0.0, -5070.0, T, X, W );
     }  
  } 
  
  if (REACTION[63]){  
  // reaction (63) H2O2 + OH -> H2O + HO2
  add_to_W_fw_2r2p ( specH2O2, specOH, specH2O, specHO2,  1.00e+13, 0.0, 1800.0, T, X, W );   
  } 
  
  if (REACTION[64]){  
  // reaction (64) H2O + HO2 -> H2O2 + OH
  add_to_W_fw_2r2p ( specH2O, specHO2, specH2O2, specOH,  2.86e+13, 0.0, 32790.0, T, X, W );  
  } 
  
  if (REACTION[65]){  
	 for (k=0; specM[k]!=specEND; k++){	
  // reaction (65) OH + H + M -> H2O + M 
    add_to_W_fw_3r2p ( specOH, specH, specM[k], specH2O, specM[k],  eff3[k]*2.20e+22, -2.0, 0.0, T, X, W );
   }  
  } 
  
  if (REACTION[66]){  
	 for (k=0; specM[k]!=specEND; k++){	  
  // reaction (66) H + H + M -> H2 + M
    add_to_W_fw_3r2p ( specH, specH, specM[k], specH2, specM[k],  eff3[k]*1.80e+18, -1.0, 0.0, T, X, W );  
   }
  } 




}

void find_dW_dx ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
                    
   long k, s;  
   spec_t X;  
   
   for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
   }

   for ( s = 0; s < ns; s++ ) {
    dWdT[s]     = 0.0;
    dWdTe[s]    = 0.0;
    dWdTv[s]    = 0.0;
    dWdQbeam[s] = 0.0;
     for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
     }
   }


  if (REACTION[1]){
  // reaction (1): C2H5 + H -> CH3 + CH3
  add_to_dW_fw_2r2p ( specC2H5, specH, specCH3, specCH3, 1.50e+13, 0.0, 0.0, T, X, dWdT, dWdrhok );
  }
  
  if (REACTION[2]){
	 for (k=0; specM[k]!=specEND; k++){
  // reaction (2): C2H4 + H + M-> C2H5 + M
    add_to_dW_fw_3r2p ( specC2H4, specH, specM[k], specC2H5, specM[k], eff1[k]*4.17e+10, 0.0, 11030.0, T, X, dWdT, dWdrhok );
	 }
  } 
  
  if (REACTION[3]){
  // reaction (3): C2H4 + H -> C2H5
  add_to_dW_fw_2r1p ( specC2H4, specH, specC2H5,       1.80e+12, 1.0, 9400.0, T, X, dWdT, dWdrhok ); 
  } 
  
  if (REACTION[4]){
  // reaction (4): C2H4 -> C2H3 + H
  add_to_dW_fw_1r2p ( specC2H4, specC2H3, specH,       1.00e+12, 0.0, 48200.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[5]){
  // reaction (5): C2H4 + H -> C2H3 + H2
  add_to_dW_fw_2r2p ( specC2H4, specH, specC2H3, specH2,  5.00e+12, 0.0, 5000.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[6]){
  // reaction (6): C2H4 + O -> CHO + CH3
  add_to_dW_fw_2r2p ( specC2H4, specO, specCHO, specCH3,  3.31e+12, 0.0, 1130.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[7]){
  // reaction (7): CHO + CH3 -> C2H4 + O
  add_to_dW_fw_2r2p ( specCHO, specCH3, specC2H4, specO,  1.58e+11, 0.0, 31180.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[8]){
  // reaction (8): C2H4 + OH -> C2H3 + H2O
  add_to_dW_fw_2r2p ( specC2H4, specOH, specC2H3, specH2O,  4.79e+12, 0.0, 1230.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[9]){
  // reaction (9): C2H3 + H2O -> C2H4 + OH
  add_to_dW_fw_2r2p ( specC2H3, specH2O, specC2H4, specOH,  1.20e+12, 0.0, 14000.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[10]){  
  // reaction (10): C2H4 + CH3 -> C2H3 + CH4
  add_to_dW_fw_2r2p ( specC2H4, specCH3, specC2H3, specCH4,  1.00e+13, 0.0, 13000.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[11]){
  // reaction (11): C2H3 + CH4 -> C2H4 + CH3
  add_to_dW_fw_2r2p ( specC2H3, specCH4, specC2H4, specCH3,  3.02e+13, 0.0, 12580.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[12]){
  // reaction (12): C2H3 -> C2H + H2
  add_to_dW_fw_1r2p ( specC2H3, specC2H, specH2,  1.80e+12, 0.0, 44800, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[13]){
  // reaction (13): C2H3 + H -> C2H2 + H2
  add_to_dW_fw_2r2p ( specC2H3, specH, specC2H2, specH2,  2.00e+13, 0.0, 2500.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[14]){
  // reaction (14): C2H3 + OH -> C2H2 + H2O
  add_to_dW_fw_2r2p ( specC2H3, specOH, specC2H2, specH2O,  7.00e+14, 0.0, 1100.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[15]){
	 for (k=0; specM[k]!=specEND; k++){	  
  // reaction (15): C2H2 + H + M -> C2H3 + M
    add_to_dW_fw_3r2p ( specC2H2, specH, specM[k], specC2H3, specM[k], eff1[k]*1.23e+11, 1.0, 10360.0, T, X, dWdT, dWdrhok );
   }
  } 
  
  if (REACTION[16]){
  // reaction (16): C2H2 + H -> C2H + H2
  add_to_dW_fw_2r2p ( specC2H2, specH, specC2H, specH2,  2.00e+14, 0.0, 19000.0, T, X, dWdT, dWdrhok );
  } 
  
  if (REACTION[17]){  
  // reaction (17): C2H2 + OH -> C2H + H2O
  add_to_dW_fw_2r2p ( specC2H2, specOH, specC2H, specH2O,  8.00e+12, 0.0, 5000.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[18]){  
  // reaction (18): C2H + H2O -> C2H2 + OH
  add_to_dW_fw_2r2p ( specC2H, specH2O, specC2H2, specOH,  5.37e+12, 0.0, 16360.0, T, X, dWdT, dWdrhok );    
  } 
  
  if (REACTION[19]){  
  // reaction (19): C2H2 + O -> C2H + OH
  add_to_dW_fw_2r2p ( specC2H2, specO, specC2H, specOH,  3.24e+15, 0.6, 12000.0, T, X, dWdT, dWdrhok );   
  } 
   
  if (REACTION[20]){ 
  // reaction (20): C2H + OH -> C2H2 + O
  add_to_dW_fw_2r2p ( specC2H, specOH, specC2H2, specO,  2.95e+14, 0.6, 910.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[21]){
	 for (k=0; specM[k]!=specEND; k++){		  
  // reaction (21): C2H + H + M -> C2H2 + M
    add_to_dW_fw_3r2p ( specC2H, specH, specM[k], specC2H2, specM[k], eff1[k]*1.10e+9, 1.0, 770.0, T, X, dWdT, dWdrhok );   
   }  
  } 
  
  if (REACTION[22]){
  // reaction (22): C2H + O2  -> CHO + CO
  add_to_dW_fw_2r2p ( specC2H, specO2, specCHO, specCO,  1.00e+13, 0.0, 7000.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[23]){
  // reaction (23): CHO + CO  -> C2H + O2
  add_to_dW_fw_2r2p ( specCHO, specCO, specC2H, specO2,  8.51e+12, 0.0, 138400.0, T, X, dWdT, dWdrhok ); 
  } 
  
  if (REACTION[24]){
  // reaction (24): C2H + OH  -> CHO + CH
  add_to_dW_fw_2r2p ( specC2H, specOH, specCHO, specCH,  8.00e+12, 0.0, 5000.0, T, X, dWdT, dWdrhok ); 
  } 
  
  if (REACTION[25]){
  // reaction (25): CH4(+M) -> CH3 + H(+M)  
   for (k=0; specM[k]!=specEND; k++){	
    add_to_dW_fw_2r3p_Lindemann ( specCH4, specM[k], specCH3, specH, specM[k], eff2[k]*6.30e+14, 0.0, 104000.0, eff2[k]*1.00e+17, 0.0, 86000.0,  T, X, dWdT, dWdrhok );
	 }
  } 
  
  if (REACTION[26]){
  // reaction (26): CH3 + H(+M)  -> CH4(+M)
   for (k=0; specM[k]!=specEND; k++){	
    add_to_dW_fw_3r2p_Lindemann ( specCH3, specH, specM[k], specCH4, specM[k], eff2[k]*5.20e+12, 0.0, -1.310e+3, eff2[k]*8.25e+14, 0.0, -1.9310e+4,  T,X, dWdT, dWdrhok );
	 }    
  } 

  
  if (REACTION[27]){
  // reaction (27): CH4 + H  -> CH3 + H2
  add_to_dW_fw_2r2p ( specCH4, specH, specCH3, specH2,  2.20e+4, 3.0, 8750.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[28]){
  // reaction (28): CH3 + H2  -> CH4 + H
  add_to_dW_fw_2r2p ( specCH3, specH2, specCH4, specH,  9.57e+2, 3.0, 8750.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[29]){  
  // reaction (29): CH4 + OH  -> CH3 + H2O
  add_to_dW_fw_2r2p ( specCH4, specOH, specCH3, specH2O,  1.60e+6, 2.1, 2460.0, T, X, dWdT, dWdrhok ); 
  } 
  
  if (REACTION[30]){  
  // reaction (30): CH3 + H2O  -> CH4 + OH
  add_to_dW_fw_2r2p ( specCH3, specH2O, specCH4, specOH,  3.02e+5, 2.1, 17422.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[31]){
  // reaction (31) CH3 + O  -> CH2O + H
  add_to_dW_fw_2r2p ( specCH3, specO, specCH2O, specH,  6.80e+13, 0.0, 0.0, T, X, dWdT, dWdrhok ); 
  } 
  
  if (REACTION[32]){
  // reaction (32) CH3 + O2  -> CH3O + O
  add_to_dW_fw_2r2p ( specCH3, specO2, specCH3O, specO,  3.00e+13, 0.0, 25652.0, T, X, dWdT, dWdrhok );     
  } 
  
  if (REACTION[33]){  
  // reaction (33) CH3 + OH  -> CH2 + H2O
  add_to_dW_fw_2r2p ( specCH3, specOH, specCH2, specH2O,  7.60e+6, 2.0, 5000.0, T, X, dWdT, dWdrhok );  
  } 
    
  if (REACTION[34]){  
  // reaction (34) CH3O + H  -> CH2O + H2
  add_to_dW_fw_2r2p ( specCH3O, specH, specCH2O, specH2,  2.00e+13, 0.0, 0.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[35]){ 
	 for (k=0; specM[k]!=specEND; k++){  
  // reaction (35) CH3O + M  -> CH2O + H + M
    add_to_dW_fw_2r3p ( specCH3O, specM[k], specCH2O, specH, specM[k] ,  eff3[k]*2.40e+13, 0.0, 28812.0, T, X, dWdT, dWdrhok );
   }  
  } 
  
  if (REACTION[36]){  
  // reaction (36) CH2 + O  -> CO + H2
  add_to_dW_fw_2r2p ( specCH2, specO, specCO, specH2,  3.00e+13, 0.0, 0.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[37]){  
  // reaction (37) CH2 + OH  -> CH + H2O
  add_to_dW_fw_2r2p ( specCH2, specOH, specCH, specH2O,  1.13e+7, 2.0, 3000.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[38]){  
  // reaction (38) CH2O + H  -> CHO + H2
  add_to_dW_fw_2r2p ( specCH2O, specH, specCHO, specH2,  5.00e+13, 0.0, 3991.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[39]){  
  // reaction (39) CH2O + OH  -> CHO + H2O
  add_to_dW_fw_2r2p ( specCH2O, specOH, specCHO, specH2O,  1.40e+14, 0.0, 1100.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[40]){  
  // reaction (40) CH + O  -> CO + H
  add_to_dW_fw_2r2p ( specCH, specO, specCO, specH,  5.70e+13, 0.0, 0.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[41]){  
  // reaction (41) CH + OH  -> CHO + H
  add_to_dW_fw_2r2p ( specCH, specOH, specCHO, specH,  3.00e+13, 0.0, 0.0, T, X, dWdT, dWdrhok );    
  } 
    
  if (REACTION[42]){  
  // reaction (42) CH + O2  -> CHO + O
  add_to_dW_fw_2r2p ( specCH, specO2, specCHO, specO,  3.30e+13, 0.0, 0.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[43]){  
  // reaction (43) CH + CO2  -> CHO + CO
  add_to_dW_fw_2r2p ( specCH, specCO2, specCHO, specCO,  8.40e+13, 0.0, 200.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[44]){  
  // reaction (44) CHO + H  -> CO + H2
  add_to_dW_fw_2r2p ( specCHO, specH, specCO, specH2,  4.00e+13, 0.0, 0.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[45]){ 
	 for (k=0; specM[k]!=specEND; k++){  
  // reaction (45) CHO + M  -> CO + H + M
    add_to_dW_fw_2r3p ( specCHO, specM[k], specCO, specH, specM[k] ,  eff3[k]*1.60e+14, 0.0, 14700.0, T, X, dWdT, dWdrhok );
   }  
  } 
  
  if (REACTION[46]){  
  // reaction (46) CO + OH  -> CO2 + H
  add_to_dW_fw_2r2p ( specCO, specOH, specCO2, specH,  1.51e+7, 1.3, -758.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[47]){  
  // reaction (47) CO2 + H  -> CO + OH
  add_to_dW_fw_2r2p ( specCO2, specH, specCO, specOH,  1.57e+9, 1.3, 19800.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[48]){  
  // reaction (48) H + O2  -> OH + O
  add_to_dW_fw_2r2p ( specH, specO2, specOH, specO,  5.00e+14, 0.0, 16800.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[49]){  
  // reaction (49) OH + O  -> H + O2
  add_to_dW_fw_2r2p ( specOH, specO, specH, specO2,  1.20e+13, 0.0, 690.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[50]){  
  // reaction (50) O + H2 -> OH + H
  add_to_dW_fw_2r2p ( specO, specH2, specOH, specH,  1.80e+10, 1.0, 8826.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[51]){  
  // reaction (51) OH + H -> O + H2
  add_to_dW_fw_2r2p ( specOH, specH, specO, specH2,  8.00e+9, 1.0, 6760.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[52]){  
  // reaction (52) H2 + OH -> H2O + H
  add_to_dW_fw_2r2p ( specH2, specOH, specH2O, specH,  1.17e+9, 1.3, 3626.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[53]){  
  // reaction (53) H2O + H -> H2 + OH
  add_to_dW_fw_2r2p ( specH2O, specH, specH2, specOH,  6.00e+9, 1.3, 18588.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[54]){  
  // reaction (54) OH + OH -> O + H2O
  add_to_dW_fw_2r2p ( specOH, specOH, specO, specH2O,  6.00e+8, 1.3, 0.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[55]){  
  // reaction (55) O + H2O -> OH + OH
  add_to_dW_fw_2r2p ( specO, specH2O, specOH, specOH,  4.00e+9, 1.3, 17029.0, T, X, dWdT, dWdrhok ); 
  } 
  
  if (REACTION[56]){
	 for (k=0; specM[k]!=specEND; k++){	  
  // reaction (56) H + O2 + M -> HO2 + M
    add_to_dW_fw_3r2p ( specH, specO2, specM[k], specHO2, specM[k] ,  eff3[k]*2.50e+18, -0.8, 0, T, X, dWdT, dWdrhok ); 
   }  
  } 
  
  if (REACTION[57]){  
  // reaction (57) H + HO2 -> OH + OH
  add_to_dW_fw_2r2p ( specH, specHO2, specOH, specOH,  1.50e+14, 0.0, 1004.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[58]){  
  // reaction (58) H + HO2 -> H2 + O2
  add_to_dW_fw_2r2p ( specH, specHO2, specH2, specO2,  2.50e+13, 0.0, 700.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[59]){  
  // reaction (59) OH + HO2 -> H2O + O2
  add_to_dW_fw_2r2p ( specOH, specHO2, specH2O, specO2,  2.00e+13, 0.0, 1000.0, T, X, dWdT, dWdrhok ); 
  } 
  
  if (REACTION[60]){  
  // reaction (60) HO2 + HO2 -> H2O2 + O2
  add_to_dW_fw_2r2p ( specHO2, specHO2, specH2O2, specO2,  2.00e+14, 0.0, 0.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[61]){
	 for (k=0; specM[k]!=specEND; k++){	 	  
  // reaction (61) H2O2 + M -> OH + OH + M
    add_to_dW_fw_2r3p ( specH2O2, specM[k], specOH, specOH, specM[k],  eff3[k]*1.30e+17, 0.0, 45500.0, T, X, dWdT, dWdrhok );    
   }
  } 
  
  if (REACTION[62]){
	 for (k=0; specM[k]!=specEND; k++){		  
  // reaction (62) OH + OH + M -> H2O2 + M 
    add_to_dW_fw_3r2p ( specOH, specOH, specM[k], specH2O2, specM[k],  eff3[k]*9.86e+14, 0.0, -5070.0, T, X, dWdT, dWdrhok );
   }  
  } 
  
  if (REACTION[63]){  
  // reaction (63) H2O2 + OH -> H2O + HO2
  add_to_dW_fw_2r2p ( specH2O2, specOH, specH2O, specHO2,  1.00e+13, 0.0, 1800.0, T, X, dWdT, dWdrhok );   
  } 
  
  if (REACTION[64]){  
  // reaction (64) H2O + HO2 -> H2O2 + OH
  add_to_dW_fw_2r2p ( specH2O, specHO2, specH2O2, specOH,  2.86e+13, 0.0, 32790.0, T, X, dWdT, dWdrhok );  
  } 
  
  if (REACTION[65]){  
	 for (k=0; specM[k]!=specEND; k++){	
  // reaction (65) OH + H + M -> H2O + M 
    add_to_dW_fw_3r2p ( specOH, specH, specM[k], specH2O, specM[k],  eff3[k]*2.20e+22, -2.0, 0.0, T, X, dWdT, dWdrhok );
   }  
  } 
  
  if (REACTION[66]){  
	 for (k=0; specM[k]!=specEND; k++){	  
  // reaction (66) H + H + M -> H2 + M
    add_to_dW_fw_3r2p ( specH, specH, specM[k], specH2, specM[k],  eff3[k]*1.80e+18, -1.0, 0.0, T, X, dWdT, dWdrhok );  
   }
  } 


}

void find_Qei(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){  
  *Qei=0.0;
}


void find_dQei_dx(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
}

