// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Ajjay Omprakas
Copyright 2022 Prasanna Thoguluva Rajendran

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


/* set all reactions to true except for testing purposes */
const static bool REACTION[32]=
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
  
      

double _kb_polynomial(long numreactant,double A1, double A2, double A3, double A4, double A5, double A, double n, double E, double T)
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
      fatal_error("numreactant can not be set to %ld in _kb_polynomial within kim2021.c",numreactant);
      kf=0.0;
  }
  assert(ke!=0.0);
  kb = kf/ke;
  
  return(kb);
  
}

double _dkbdT_polynomial(long numreactant,double A1, double A2, double A3, double A4, double A5, double A, double n, double E, double T)
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
      fatal_error("numreactant can not be set to %ld in _dkbdT_polynomial within kim2021.c",numreactant);
      dkbdT=0.0;
  }
 
  return(dkbdT);
  
}


void find_W_Kim2021 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
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
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 1.216e20, -1.2140, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 1.216e20, -1.2140, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM1[k], specN, specN,   specN2, specM1[k], kb , N, W);      
    }
  }
      

  if (REACTION[2]){
    for (k=0; specM2[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM2[k],   specN, specN, specM2[k], 7.0e21, -1.6, 113200.0*R, TTv, X, W );
      
      
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 7.0e21, -1.6, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 7.0e21, -1.6, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM2[k], specN, specN,   specN2, specM2[k], kb , N, W); 
    }
  }

  if (REACTION[3]){ 
    for (k=0; specM3[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM3[k],   specN, specN, specM3[k], 3.591e20, -1.226, 113200.0*R, TTv, X, W );
      
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.591e20, -1.226, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.591e20, -1.226, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM3[k], specN, specN,   specN2, specM3[k], kb , N, W);
    }
  }

  if (REACTION[4]){ 
    for (k=0; specM4[k]!=specEND; k++) {
      add_to_W_fw_2r3p ( specN2, specM4[k],   specN, specN, specM4[k], 3.0e22, -1.6, 113200.0*R, TTv, X, W );
      
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.0e22, -1.6, 113200.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.0e22, -1.6, 113200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM4[k], specN, specN,   specN2, specM4[k], kb , N, W);
    }
  }

  if (REACTION[5]){
    for (k=0; specM5[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM5[k],   specO, specO, specM5[k], 3.354e15, -0.2726, 59500.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.354e15, -0.2726, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.354e15, -0.2726, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM5[k], specO, specO,   specO2, specM5[k], kb , N, W);
    }
  }

  if (REACTION[6]){
    for (k=0; specM6[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM6[k],   specO, specO, specM6[k], 1.117e25, -2.585, 59500.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.117e25, -2.585, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.117e25, -2.585, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM6[k], specO, specO,   specO2, specM6[k], kb , N, W);
    }
  }

  if (REACTION[7]){
    for (k=0; specM3[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM3[k],   specO, specO, specM3[k], 1.0e22, -1.500, 59500.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.0e22, -1.500, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.0e22, -1.500, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM3[k], specO, specO,   specO2, specM3[k], kb , N, W);
    }
  }

  if (REACTION[8]){
    for (k=0; specM4[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specO2, specM4[k],   specO, specO, specM4[k], 3.00e21, -1.50, 59500.0*R, TTv, X, W );
     
      if (T < 20000.0)
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.00e21, -1.50, 59500.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.00e21, -1.50, 59500.0*R, Tblim); 
        
      add_to_W_3r2p ( specM4[k], specO, specO,   specO2, specM4[k], kb , N, W);
    }
  }

  if (REACTION[9]){
    for (k=0; specM7[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM7[k],   specN, specO, specM7[k], 1.450e15, 0.0, 75200.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 1.450e15, 0.0, 75200.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 1.450e15, 0.0, 75200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM7[k], specN, specO,   specNO, specM7[k], kb , N, W);
    }
  }

  if (REACTION[10]){
    for (k=0; specM8[k]!=specEND; k++){
      add_to_W_fw_2r3p ( specNO, specM8[k],   specN, specO, specM8[k], 9.640e14, 0.0, 75200.0*R, TTv, X, W );
      
      if (T < 20000.0)
        kb=_kb_polynomial(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 9.640e14, 0.0, 75200.0*R, Tblim);
      else
        kb=_kb_polynomial(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 9.640e14, 0.0, 75200.0*R, Tblim); 
        
      add_to_W_3r2p ( specM8[k], specN, specO,   specNO, specM8[k], kb , N, W);
    }
  }

  if (REACTION[11]){
    add_to_W_fw_2r2p ( specN2, specO,   specNO, specN, 5.700e12, 0.42, 42938.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -1.789e-1, 1.728, -2.172e-1, -3.733, -2.285e-4, 5.700e12, 0.42, 42938.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 7.441e-2, 2.852, 1.054, -5.303, 1.909e-1, 5.700e12, 0.42, 42938.0*R, Tblim); 
        
      add_to_W_2r2p ( specNO, specN,  specN2, specO, kb , N, W);
  }

  if (REACTION[12]){
    add_to_W_fw_2r2p ( specNO, specO,   specO2, specN, 8.4000e12, 0.0, 19400.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -1.673e-1, -1.390, -1.656e-1, -1.551, -1.102e-4, 8.4000e12, 0.0, 19400.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 1.744e-2, -1.563, 2.354e-1, -1.233, -3.250e-1, 8.4000e12, 0.0, 19400.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specN,  specO, specNO, kb , N, W);
  }

  if (REACTION[13]){
    add_to_W_fw_2r2p ( specN, specO,   specNOplus, speceminus, 8.80e08, 1.0, 31900.0*R, T, X, W );
    
    if (Te < 20000.0)
        kb=_kb_polynomial(2, 1.907e-1, -7.976, -1.848, -3.255, -1.662e-3, 8.80e08, 1.0, 31900.0*R, Teblim);
      else
        kb=_kb_polynomial(2, -1.361e-1, -6.297, -1.866, -5.991, 1.384, 8.80e08, 1.0, 31900.0*R, Teblim); 
        
      add_to_W_2r2p ( specNOplus, speceminus,  specO, specN, kb , N, W);
  }

  if (REACTION[14]){
    add_to_W_fw_2r2p ( specO, specO,   specO2plus, speceminus, 7.10e02, 2.70, 80600.0*R, T, X, W );
   
    if (Te < 20000.0)
        kb=_kb_polynomial(2, -7.183e-3, -7.603, -2.099, -8.070, -1.989e-3, 7.10e02, 2.70, 80600.0*R, Teblim);
      else
        kb=_kb_polynomial(2, -2.428e-2, -4.074, -3.091e-1, -1.342e1, 1.831, 7.10e02, 2.70, 80600.0*R, Teblim); 
        
      add_to_W_2r2p ( specO2plus, speceminus,  specO, specO, kb , N, W);
  }

  if (REACTION[15]){
    add_to_W_fw_2r3p ( specN, speceminus,   specNplus, speceminus, speceminus, 2.50e034, -3.820, 168600.0*R, Te, X, W );
    
    if (Te < 20000.0)
        kb=_kb_polynomial(3, -1.217, -3.006, -2.354, -1.675e1, -2.195e-3, 2.50e034, -3.820, 168600.0*R, Teblim);
      else
        kb=_kb_polynomial(3, 1.103, 2.700, 7.541, -2.188e1, -2.910, 2.50e034, -3.820, 168600.0*R, Teblim); 
        
      add_to_W_3r2p ( speceminus, speceminus, specNplus,  speceminus, specN, kb , N, W);
      
  }

  if (REACTION[16]){
    add_to_W_fw_2r3p ( specO, speceminus,   specOplus, speceminus, speceminus, 3.90e033, -3.78, 158500.0*R, Te, X, W );
    
    if (Te < 20000.0)
        kb=_kb_polynomial(3, -2.660e-1, -5.390, -1.747, -1.575e1, -7.662e-4, 3.90e033, -3.78, 158500.0*R, Teblim);
      else
        kb=_kb_polynomial(3, 1.789, 1.711e1, 1.527e1, -4.975e1, 9.411, 3.90e033, -3.78, 158500.0*R, Teblim); 
        
      add_to_W_3r2p ( speceminus, speceminus, specOplus,  speceminus, specO, kb , N, W);
      
  }

  if (REACTION[17]){
    add_to_W_fw_2r2p ( specN2, specO2plus,   specN2plus, specO2, 9.90e12, 0.0, 40700.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -1.970e-1, 1.031, -2.049e-1, -4.005, -9.866e-5, 9.90e12, 0.0, 40700.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 4.428e-2, 1.462, 6.450e-1, -4.638, -3.600e-2, 9.90e12, 0.0, 40700.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specN2plus,  specO2plus, specN2, kb , N, W);
  }

  if (REACTION[18]){
    add_to_W_fw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.40e13, -1.08, 12800.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -6.300e-2, -5.419e-1, -4.449e-2, -1.266, -1.089e-4, 3.40e13, -1.08, 12800.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 3.518e-2, -7.715e-2, 4.619e-1, -1.923, 9.394e-2, 3.40e13, -1.08, 12800.0*R, Tblim); 
        
      add_to_W_2r2p ( specN2, specOplus,  specNOplus, specN, kb , N, W);
  }

  if (REACTION[19]){
    add_to_W_fw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.00e12, 0.5, 77200.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -4.092e-1, -7.798e-1, -4.273e-1, -7.618, -4.475e-4, 1.00e12, 0.5, 77200.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 1.270e-1, 6.372e-1, 1.751, -9.528, -4.012e-2, 1.00e12, 0.5, 77200.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specNplus,  specNOplus, specO, kb , N, W);
  }

  if (REACTION[20]){
    add_to_W_fw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.40e13, 0.41, 32600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -3.056e-2, 1.762, -8.580e-2, -3.265, -2.168e-4, 2.40e13, 0.41, 32600.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 9.436e-2, 3.786, 1.322, -6.192, 7.723e-1, 2.40e13, 0.41, 32600.0*R, Tblim); 
        
      add_to_W_2r2p ( specNO, specO2plus,  specO2, specNOplus, kb , N, W);
  }

  if (REACTION[21]){
    add_to_W_fw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.20e13, 0.0, 35500.0*R, T, X, W );
   
    if (T < 20000.0)
        kb=_kb_polynomial(2, -4.862e-2, 1.066, -7.350e-2, -3.537, -8.701e-5, 7.20e13, 0.0, 35500.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 6.423e-2, 2.396, 9.127e-1, -5.527, 5.453e-1, 7.20e13, 0.0, 35500.0*R, Tblim); 
        
      add_to_W_2r2p ( specO, specN2plus,  specN, specNOplus, kb , N, W);
  }

  if (REACTION[22]){
    add_to_W_fw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.70e13, 0.14, 28600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -2.113e-1, -1.152, -1.758e-1, -2.803, -1.205e-4, 8.70e13, 0.14, 28600.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 1.523e-2, -1.586, 1.942e-1, -2.103, -4.874e-1, 8.70e13, 0.14, 28600.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specNplus,  specN, specO2plus, kb , N, W);
  }

  if (REACTION[23]){
    add_to_W_fw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.40e05, 1.90, 26600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -1.673e-1, -1.965, -1.656e-1, -2.619, -1.102e-4, 1.40e05, 1.90, 26600.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 1.744e-2, -2.138, 2.354e-1, -2.302, -3.250e-1, 1.40e05, 1.90, 26600.0*R, Tblim); 
        
      add_to_W_2r2p ( specO2, specNplus,  specNO, specOplus, kb , N, W);
  }

  if (REACTION[24]){
    add_to_W_fw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.20e12, 0.29, 48600.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, -1.978e-1, 3.723e-1, -2.514e-1, -4.815, -3.270e-4, 7.20e12, 0.29, 48600.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 1.118e-1, 2.223, 1.557, -7.425, 4.473e-1, 7.20e12, 0.29, 48600.0*R, Tblim); 
        
      add_to_W_2r2p ( specN, specO2plus,  specO, specNOplus, kb , N, W);
  }

  if (REACTION[25]){
    add_to_W_fw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.10e11, 0.36, 22800.0*R, T, X, W );
    
    if (T < 20000.0)
        kb=_kb_polynomial(2, 1.438e-2, 1.608, -2.901e-2, -2.271, 2.188e-5, 9.10e11, 0.36, 22800.0*R, Tblim);
      else
        kb=_kb_polynomial(2, 2.905e-2, 2.473, 4.508e-1, -3.603, 4.514e-1, 9.10e11, 0.36, 22800.0*R, Tblim); 
        
      add_to_W_2r2p ( specO, specN2plus,  specN2, specOplus, kb , N, W);
  }
  
  if (REACTION[26]){
    add_to_W_fw_2r2p ( specN, specN,   specN2plus, speceminus, 4.40e07, 1.5, 67500.0*R, T, X, W );
    
    if (Te < 20000.0)
        kb=_kb_polynomial(2, 1.420e-1, -6.909, -1.922, -6.792, -1.749e-3, 4.40e07, 1.5, 67500.0*R, Teblim);
      else
        kb=_kb_polynomial(2, -7.185e-2, -3.901, -9.534e-1, -1.152e1, 1.929, 4.40e07, 1.5, 67500.0*R, Teblim); 
        
      add_to_W_2r2p ( speceminus, specN2plus,  specN, specN, kb , N, W);
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
    add_to_W_fw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 1.20e25, -1.6, 113200.0*R, Te, X, W );
    add_to_W_fw_3r2p ( specN, specN, speceminus,   specN2, speceminus,   1.20e25, -1.6, 113200.0*R, Te, X, W );
  }   


}



/* Verify the validity of the dW terms at node i=10,j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU 
 * Make sure to verify the dW terms over a wide range of temperatures and mass fractions 
 * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives
 * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()
 */ 

void find_dW_dx_Kim2021 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s, spec;                    
  spec_t N;
  double TTv,TTe,TvTe,Tblim,Teblim,R,kb,dkbdT, dkbdTv, dkbdTe;
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
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 1.216e20, -1.2140, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 1.216e20, -1.2140, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 1.216e20, -1.2140, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 1.216e20, -1.2140, 113200.0*R, Tblim);
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
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 7.0e21, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 7.0e21, -1.6, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 7.0e21, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 7.0e21, -1.6, 113200.0*R, Tblim);
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
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.591e20, -1.226, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.591e20, -1.226, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.591e20, -1.226, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.591e20, -1.226, 113200.0*R, Tblim);
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
        kb=_kb_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.0e22, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 2.491, 7.155e-1, 2.091, -1.169e1, 5.921e-3, 3.0e22, -1.6, 113200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.0e22, -1.6, 113200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, -8.575e-1, -5.071e-1, -9.211, -1.734e1, 1.024e1, 3.0e22, -1.6, 113200.0*R, Tblim);
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
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.354e15, -0.2726, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.354e15, -0.2726, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.354e15, -0.2726, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.354e15, -0.2726, 59500.0*R, Tblim);
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
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.117e25, -2.585, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.117e25, -2.585, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.117e25, -2.585, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.117e25, -2.585, 59500.0*R, Tblim);
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
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.0e22, -1.500, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 1.0e22, -1.500, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.0e22, -1.500, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 1.0e22, -1.500, 59500.0*R, Tblim);
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
        kb=_kb_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.00e21, -1.50, 59500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 1.567, 1.217, 1.909, -6.281, 5.237e-3, 3.00e21, -1.50, 59500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.00e21, -1.50, 59500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(3, -2.087, -3.315e1, -2.595e1, 4.370e1, -1.195e1, 3.00e21, -1.50, 59500.0*R, Tblim);
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
        kb=_kb_polynomial(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 1.450e15, 0.0, 75200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 1.450e15, 0.0, 75200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 1.450e15, 0.0, 75200.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 1.450e15, 0.0, 75200.0*R, Tblim);
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
        kb=_kb_polynomial(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 9.640e14, 0.0, 75200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(3, 2.093, -6.229e-1, 2.028, -7.872, 5.586e-3, 9.640e14, 0.0, 75200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 9.640e14, 0.0, 75200.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(3, -1.640, -2.142e1, -1.964e1, 1.910e1, -2.422, 9.640e14, 0.0, 75200.0*R, Tblim);
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
        kb=_kb_polynomial(2, -1.789e-1, 1.728, -2.172e-1, -3.733, -2.285e-4, 5.700e12, 0.42, 42938.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -1.789e-1, 1.728, -2.172e-1, -3.733, -2.285e-4, 5.700e12, 0.42, 42938.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 7.441e-2, 2.852, 1.054, -5.303, 1.909e-1, 5.700e12, 0.42, 42938.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 7.441e-2, 2.852, 1.054, -5.303, 1.909e-1, 5.700e12, 0.42, 42938.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specNO, specN,  specN2, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[12]){
    add_to_dW_fw_2r2p ( specNO, specO,   specO2, specN, 8.4000e12, 0.0, 19400.0*R, T, X, dWdT, dWdrhok);
   
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -1.673e-1, -1.390, -1.656e-1, -1.551, -1.102e-4, 8.4000e12, 0.0, 19400.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -1.673e-1, -1.390, -1.656e-1, -1.551, -1.102e-4, 8.4000e12, 0.0, 19400.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 1.744e-2, -1.563, 2.354e-1, -1.233, -3.250e-1, 8.4000e12, 0.0, 19400.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 1.744e-2, -1.563, 2.354e-1, -1.233, -3.250e-1, 8.4000e12, 0.0, 19400.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specN,  specNO, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }


  if (REACTION[13]){
    add_to_dW_fw_2r2p ( specN, specO,   specNOplus, speceminus, 8.80e08, 1.0, 31900.0*R, T, X, dWdT, dWdrhok);
    
    if (Te < 20000.0)
      {
        kb=_kb_polynomial(2, 1.907e-1, -7.976, -1.848, -3.255, -1.662e-3, 8.80e08, 1.0, 31900.0*R, Teblim);
        dkbdTe=_dkbdT_polynomial(2, 1.907e-1, -7.976, -1.848, -3.255, -1.662e-3, 8.80e08, 1.0, 31900.0*R, Teblim);
      }
      else
      {
        kb=_kb_polynomial(2, -1.361e-1, -6.297, -1.866, -5.991, 1.384, 8.80e08, 1.0, 31900.0*R, Teblim); 
        dkbdTe=_dkbdT_polynomial(2, -1.361e-1, -6.297, -1.866, -5.991, 1.384, 8.80e08, 1.0, 31900.0*R, Teblim);
      }
      
      dkbdT=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specNOplus, speceminus,  specO, specN,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[14]){
    add_to_dW_fw_2r2p ( specO, specO,   specO2plus, speceminus, 7.10e02, 2.70, 80600.0*R, T, X, dWdT, dWdrhok);
    
    if (Te < 20000.0)
      {
        kb=_kb_polynomial(2, -7.183e-3, -7.603, -2.099, -8.070, -1.989e-3, 7.10e02, 2.70, 80600.0*R, Teblim);
        dkbdTe=_dkbdT_polynomial(2, -7.183e-3, -7.603, -2.099, -8.070, -1.989e-3, 7.10e02, 2.70, 80600.0*R, Teblim);
      }
      else
      {
        kb=_kb_polynomial(2, -2.428e-2, -4.074, -3.091e-1, -1.342e1, 1.831, 7.10e02, 2.70, 80600.0*R, Teblim); 
        dkbdTe=_dkbdT_polynomial(2, -2.428e-2, -4.074, -3.091e-1, -1.342e1, 1.831, 7.10e02, 2.70, 80600.0*R, Teblim);
      }
      
      dkbdT=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2plus, speceminus,  specO, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
    
  }

  if (REACTION[15]){
    add_to_dW_fw_2r3p ( specN, speceminus,   specNplus, speceminus, speceminus, 2.50e034, -3.820, 168600.0*R, Te, X, dWdTe, dWdrhok);
    
    if (Te < 20000.0)
      {
        kb=_kb_polynomial(3, -1.217, -3.006, -2.354, -1.675e1, -2.195e-3, 2.50e034, -3.820, 168600.0*R, Teblim);
        dkbdTe=_dkbdT_polynomial(3, -1.217, -3.006, -2.354, -1.675e1, -2.195e-3, 2.50e034, -3.820, 168600.0*R, Teblim);
      }
      else
      {
        kb=_kb_polynomial(3, 1.103, 2.700, 7.541, -2.188e1, -2.910, 2.50e034, -3.820, 168600.0*R, Teblim); 
        dkbdTe=_dkbdT_polynomial(3, 1.103, 2.700, 7.541, -2.188e1, -2.910, 2.50e034, -3.820, 168600.0*R, Teblim);
      }
      
      dkbdT=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( speceminus, speceminus, specNplus,  speceminus, specN,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
     
  }

  if (REACTION[16]){
    add_to_dW_fw_2r3p ( specO, speceminus,   specOplus, speceminus, speceminus, 3.90e033, -3.78, 158500.0*R, Te, X, dWdTe, dWdrhok);
    
    if (Te < 20000.0)
      {
        kb=_kb_polynomial(3, -2.660e-1, -5.390, -1.747, -1.575e1, -7.662e-4, 3.90e033, -3.78, 158500.0*R, Teblim);
        dkbdTe=_dkbdT_polynomial(3, -2.660e-1, -5.390, -1.747, -1.575e1, -7.662e-4, 3.90e033, -3.78, 158500.0*R, Teblim);
      }
      else
      {
        kb=_kb_polynomial(3, 1.789, 1.711e1, 1.527e1, -4.975e1, 9.411, 3.90e033, -3.78, 158500.0*R, Teblim); 
        dkbdTe=_dkbdT_polynomial(3, 1.789, 1.711e1, 1.527e1, -4.975e1, 9.411, 3.90e033, -3.78, 158500.0*R, Teblim);
      }
      
      dkbdT=0.0;
      dkbdTv=0.0;
      add_to_dW_3r2p ( speceminus, speceminus, specOplus,  speceminus, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
     
  }

  if (REACTION[17]){
    add_to_dW_fw_2r2p ( specN2, specO2plus,   specN2plus, specO2, 9.90e12, 0.0, 40700.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -1.970e-1, 1.031, -2.049e-1, -4.005, -9.866e-5, 9.90e12, 0.0, 40700.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -1.970e-1, 1.031, -2.049e-1, -4.005, -9.866e-5, 9.90e12, 0.0, 40700.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 4.428e-2, 1.462, 6.450e-1, -4.638, -3.600e-2, 9.90e12, 0.0, 40700.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 4.428e-2, 1.462, 6.450e-1, -4.638, -3.600e-2, 9.90e12, 0.0, 40700.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specN2plus, specO2plus,  specN2,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[18]){
    add_to_dW_fw_2r2p ( specNOplus, specN,   specOplus, specN2, 3.40e13, -1.08, 12800.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -6.300e-2, -5.419e-1, -4.449e-2, -1.266, -1.089e-4, 3.40e13, -1.08, 12800.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -6.300e-2, -5.419e-1, -4.449e-2, -1.266, -1.089e-4, 3.40e13, -1.08, 12800.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 3.518e-2, -7.715e-2, 4.619e-1, -1.923, 9.394e-2, 3.40e13, -1.08, 12800.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 3.518e-2, -7.715e-2, 4.619e-1, -1.923, 9.394e-2, 3.40e13, -1.08, 12800.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specN2, specOplus,  specNOplus, specN,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[19]){
    add_to_dW_fw_2r2p ( specNOplus, specO,   specNplus, specO2, 1.00e12, 0.5, 77200.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -4.092e-1, -7.798e-1, -4.273e-1, -7.618, -4.475e-4, 1.00e12, 0.5, 77200.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -4.092e-1, -7.798e-1, -4.273e-1, -7.618, -4.475e-4, 1.00e12, 0.5, 77200.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 1.270e-1, 6.372e-1, 1.751, -9.528, -4.012e-2, 1.00e12, 0.5, 77200.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 1.270e-1, 6.372e-1, 1.751, -9.528, -4.012e-2, 1.00e12, 0.5, 77200.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specNplus,  specNOplus, specO,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[20]){
    add_to_dW_fw_2r2p ( specNOplus, specO2,   specO2plus, specNO, 2.40e13, 0.41, 32600.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -3.056e-2, 1.762, -8.580e-2, -3.265, -2.168e-4, 2.40e13, 0.41, 32600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -3.056e-2, 1.762, -8.580e-2, -3.265, -2.168e-4, 2.40e13, 0.41, 32600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 9.436e-2, 3.786, 1.322, -6.192, 7.723e-1, 2.40e13, 0.41, 32600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 9.436e-2, 3.786, 1.322, -6.192, 7.723e-1, 2.40e13, 0.41, 32600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specNO, specO2plus,  specO2, specNOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[21]){
    add_to_dW_fw_2r2p ( specNOplus, specN,   specN2plus, specO, 7.20e13, 0.0, 35500.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -4.862e-2, 1.066, -7.350e-2, -3.537, -8.701e-5, 7.20e13, 0.0, 35500.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -4.862e-2, 1.066, -7.350e-2, -3.537, -8.701e-5, 7.20e13, 0.0, 35500.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 6.423e-2, 2.396, 9.127e-1, -5.527, 5.453e-1, 7.20e13, 0.0, 35500.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 6.423e-2, 2.396, 9.127e-1, -5.527, 5.453e-1, 7.20e13, 0.0, 35500.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO, specN2plus,  specN, specNOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[22]){
    add_to_dW_fw_2r2p ( specO2plus, specN,   specNplus, specO2, 8.70e13, 0.14, 28600.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -2.113e-1, -1.152, -1.758e-1, -2.803, -1.205e-4, 8.70e13, 0.14, 28600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -2.113e-1, -1.152, -1.758e-1, -2.803, -1.205e-4, 8.70e13, 0.14, 28600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 1.523e-2, -1.586, 1.942e-1, -2.103, -4.874e-1, 8.70e13, 0.14, 28600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 1.523e-2, -1.586, 1.942e-1, -2.103, -4.874e-1, 8.70e13, 0.14, 28600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specNplus,  specN, specO2plus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[23]){
    add_to_dW_fw_2r2p ( specOplus, specNO,   specNplus, specO2, 1.40e05, 1.90, 26600.0*R, T, X, dWdT, dWdrhok);
   
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -1.673e-1, -1.965, -1.656e-1, -2.619, -1.102e-4, 1.40e05, 1.90, 26600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -1.673e-1, -1.965, -1.656e-1, -2.619, -1.102e-4, 1.40e05, 1.90, 26600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 1.744e-2, -2.138, 2.354e-1, -2.302, -3.250e-1, 1.40e05, 1.90, 26600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 1.744e-2, -2.138, 2.354e-1, -2.302, -3.250e-1, 1.40e05, 1.90, 26600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO2, specNplus,  specNO, specOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[24]){
    add_to_dW_fw_2r2p ( specNOplus, specO,   specO2plus, specN, 7.20e12, 0.29, 48600.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, -1.978e-1, 3.723e-1, -2.514e-1, -4.815, -3.270e-4, 7.20e12, 0.29, 48600.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, -1.978e-1, 3.723e-1, -2.514e-1, -4.815, -3.270e-4, 7.20e12, 0.29, 48600.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 1.118e-1, 2.223, 1.557, -7.425, 4.473e-1, 7.20e12, 0.29, 48600.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 1.118e-1, 2.223, 1.557, -7.425, 4.473e-1, 7.20e12, 0.29, 48600.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specN, specO2plus,  specO, specNOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (REACTION[25]){
    add_to_dW_fw_2r2p ( specOplus, specN2,   specN2plus, specO, 9.10e11, 0.36, 22800.0*R, T, X, dWdT, dWdrhok);
    
    if (T < 20000.0)
      {
        kb=_kb_polynomial(2, 1.438e-2, 1.608, -2.901e-2, -2.271, 2.188e-5, 9.10e11, 0.36, 22800.0*R, Tblim);
        dkbdT=_dkbdT_polynomial(2, 1.438e-2, 1.608, -2.901e-2, -2.271, 2.188e-5, 9.10e11, 0.36, 22800.0*R, Tblim);
      }
      else
      {
        kb=_kb_polynomial(2, 2.905e-2, 2.473, 4.508e-1, -3.603, 4.514e-1, 9.10e11, 0.36, 22800.0*R, Tblim); 
        dkbdT=_dkbdT_polynomial(2, 2.905e-2, 2.473, 4.508e-1, -3.603, 4.514e-1, 9.10e11, 0.36, 22800.0*R, Tblim);
      }
      
      dkbdTe=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( specO, specN2plus,  specN2, specOplus,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (REACTION[26]){
    add_to_dW_fw_2r2p ( specN, specN,   specN2plus, speceminus, 4.40e07, 1.5, 67500.0*R, T, X, dWdT, dWdrhok);
    
    if (Te < 20000.0)
      {
        kb=_kb_polynomial(2, 1.420e-1, -6.909, -1.922, -6.792, -1.749e-3, 4.40e07, 1.5, 67500.0*R, Teblim);
        dkbdTe=_dkbdT_polynomial(2, 1.420e-1, -6.909, -1.922, -6.792, -1.749e-3, 4.40e07, 1.5, 67500.0*R, Teblim);
      }
      else
      {
        kb=_kb_polynomial(2, -7.185e-2, -3.901, -9.534e-1, -1.152e1, 1.929, 4.40e07, 1.5, 67500.0*R, Teblim); 
        dkbdTe=_dkbdT_polynomial(2, -7.185e-2, -3.901, -9.534e-1, -1.152e1, 1.929, 4.40e07, 1.5, 67500.0*R, Teblim);
      }
      
      dkbdT=0.0;
      dkbdTv=0.0;
      add_to_dW_2r2p ( speceminus, specN2plus,  specN, specN,  kb , N, dkbdT, dkbdTv, dkbdTe, dWdrhok, dWdT, dWdTv, dWdTe);
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
    add_to_dW_fw_2r3p ( specN2, speceminus,   specN, specN, speceminus, 1.20e25, -1.6, 113200.0*R, Te, X, dWdTe, dWdrhok);
    add_to_dW_fw_3r2p ( specN, specN, speceminus,   specN2, speceminus, 1.20e25, -1.6, 113200.0*R, Te, X, dWdTe, dWdrhok);
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





void find_Qei_Kim2021(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){

    if (REACTION[16]) 
      add_to_Qei(gl,Te,specO,_ionizationpot(specO), 3.9e33/calA*pow(Te,-3.78)*exp(-158500.0/Te), rhok, Qei);
    if (REACTION[15]) 
      add_to_Qei(gl,Te,specN,_ionizationpot(specN), 2.5e34/calA*pow(Te,-3.82)*exp(-168600.0/Te), rhok, Qei);
 
}



void find_dQei_dx_Kim2021(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){

    if (REACTION[16]) 
      add_to_dQei(gl,Te,specO,_ionizationpot(specO), 3.9e33/calA*pow(Te,-3.78)*exp(-158500.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);
    if (REACTION[15]) 
      add_to_dQei(gl,Te,specN,_ionizationpot(specN), 2.5e34/calA*pow(Te,-3.82)*exp(-168600.0/Te), 0.0, rhok, dQeidrhok, dQeidTe);

}
