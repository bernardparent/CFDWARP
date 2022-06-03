// SPDX-License-Identifier: BSD-2-Clause
/*
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


/* 
 * Refs.

  Capitelli, M. and Bruno, D. and Larrichiuta, A. "Fundamental Aspects of Plasma Chemical Physics, Transport", Springer, 2015

 * */

#include <model/thermo/_thermo.h>
#include <model/transport/_transport.h>
#include <assert.h>

#define Runiv 1.987 //cal/(g-mol K)
#define kBol 1.38066E-16 // erg/K

#define EPC_NONE 0
#define EPC_GUPTAYOS 1
#define EPC_RAIZER 2
#define EPC_NRL 3
#define EPC EPC_GUPTAYOS  //Use EPC_GUPTAYOS. Other methods are for testing purposes only.

#define METHOD1 1  // Eq(27), first Chapman-Enskog approximation
#define METHOD2 2  // Eq(30), approximation to method 1, yields results close to method 1 for air mixture
#define METHOD3 3  // Eq(40a), approximation to method 1 ??? or method 2??
#define METHOD4 4   // Eq(49) : 5-temperature model, thermal non-equilibirum
#define METHOD METHOD2  //Use METHOD2. Other methods are for testing purposes only


double collision_integral_curvefit11(long s, long r, double T){
  double A,B,C,D;
  
  T=min(max(1000.0,T),30000.0);
  
  switch (smap[s]){
    case SMAP_eminus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O: if(T >= 1000.0 && T< 9000.0) {A=0.0164; B=-0.2431; C=1.1231; D=-1.5561;}
          else if(T >= 9000.0 && T<= 30000.0) {A=-0.2027; B=5.6428; C=-51.5646; D=155.6091;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_N: A=0.0; B=0.0; C=0.0; D=1.6094; break;    
        case SMAP_NO: if(T >= 1000.0 && T< 8000.0) {A=-0.2202; B=5.2265; C=-40.5659; D=104.7126;}
          else if(T >= 8000.0 && T<= 30000.0) {A=-0.2871; B=8.3757; C=-81.3787; D=265.6292;}
          else fatal_error("Temperature of %.15E K is out of range of validitity (1000-30000 K) for species pair %ld,%ld.",s,r);
          break;    
        case SMAP_O2: if(T >= 1000.0 && T< 9000.0) {A=0.0241; B=-0.3467; C=1.3887; D=-0.0110;}
          else if(T >= 9000.0 && T<= 30000.0) {A=0.0025; B=-0.0742; C=0.7235; D=-0.2116;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;    
        case SMAP_N2: A=0.1147; B=-2.8945; C=24.5080; D=-67.3691; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_Oplus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O: A=0.0; B=-0.0034; C=-0.0572; D=4.9901; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_Nplus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N: A=0.0; B=-0.0033; C=-0.0572; D=5.0452; break;    
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_O2plus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_N2plus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_NOplus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=23.8237; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_NO: A=0.0; B=-0.0047; C=-0.0551; D=4.8737; break;    
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_O:
      switch (smap[r]){
        case SMAP_eminus: if(T >= 1000.0 && T< 9000.0) {A=0.0164; B=-0.2431; C=1.1231; D=-1.5561;}
          else if(T >= 9000.0 && T<= 30000.0) {A=-0.2027; B=5.6428; C=-51.5646; D=155.6091;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_Oplus: A=0.0; B=-0.0034; C=-0.0572; D=4.9901; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O: A=0.0; B=-0.0034; C=-0.0572; D=4.9901; break;
        case SMAP_N: A=0.0; B=0.0048; C=-0.4195; D=5.7774; break;    
        case SMAP_NO: A=0.0; B=-0.0179; C=0.0152; D=3.9996; break;    
        case SMAP_O2: A=0.0; B=-0.0226; C=0.1300; D=3.3363; break;    
        case SMAP_N2: A=0.0; B=-0.0139; C=-0.0825; D=4.5785; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_N:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=0.0; D=1.6094; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_Nplus: A=0.0; B=-0.0033; C=-0.0572; D=5.0452; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O: A=0.0; B=0.0048; C=-0.4195; D=5.7774; break;
        case SMAP_N: A=0.0; B=-0.0033; C=-0.0572; D=5.0452; break;    
        case SMAP_NO: A=0.0; B=-0.0185; C=0.0118; D=4.0590; break;    
        case SMAP_O2: A=0.0; B=-0.0179; C=0.0152; D=3.9996; break;    
        case SMAP_N2: A=0.0; B=-0.0194; C=0.0119; D=4.1055; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_NO:
      switch (smap[r]){
        case SMAP_eminus: if(T >= 1000.0 && T< 8000.0) {A=-0.2202; B=5.2265; C=-40.5659; D=104.7126;}
          else if(T >= 8000.0 && T<= 30000.0) {A=-0.2871; B=8.3757; C=-81.3787; D=265.6292;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_NOplus: A=0.0; B=-0.0047; C=-0.0551; D=4.8737; break;
        case SMAP_O: A=0.0; B=-0.0179; C=0.0152; D=3.9996; break;
        case SMAP_N: A=0.0; B=-0.0185; C=0.0118; D=4.0590; break;    
        case SMAP_NO: A=0.0; B=-0.0364; C=0.3825; D=2.4718; break;    
        case SMAP_O2: A=0.0; B=-0.0438; C=0.5352; D=1.7252; break;    
        case SMAP_N2: A=0.0; B=-0.0291; C=0.2324; D=3.2082; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_O2:
      switch (smap[r]){
        case SMAP_eminus: if(T >= 1000.0 && T< 9000.0) {A=0.0241; B=-0.3467; C=1.3887; D=-0.0110;}
          else if(T >= 9000.0 && T<= 30000.0) {A=0.0025; B=-0.0742; C=0.7235; D=-0.2116;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O: A=0.0; B=-0.0226; C=0.1300; D=3.3363; break;
        case SMAP_N: A=0.0; B=-0.0179; C=0.0152; D=3.9996; break;    
        case SMAP_NO: A=0.0; B=-0.0438; C=0.5352; D=1.7252; break;    
        case SMAP_O2: A=0.0; B=-0.0410; C=0.4977; D=1.8302; break;    
        case SMAP_N2: A=0.0; B=-0.0465; C=0.5729; D=1.6185; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_N2:
      switch (smap[r]){
        case SMAP_eminus: A=0.1147; B=-2.8945; C=24.5080; D=-67.3691; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.8543; break;
        case SMAP_O: A=0.0; B=-0.0139; C=-0.0825; D=4.5785; break;
        case SMAP_N: A=0.0; B=-0.0194; C=0.0119; D=4.1055; break;    
        case SMAP_NO: A=0.0; B=-0.0291; C=0.2324; D=3.2082; break;    
        case SMAP_O2: A=0.0; B=-0.0465; C=0.5729; D=1.6185; break;    
        case SMAP_N2: A=0.0; B=-0.0112; C=-0.1182; D=4.8464; break;    
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    
    default:
      fatal_error("spec %ld does not have collision integral data.",s);
  }
  
  return(exp(D)*pow(T,A*log(T)*log(T)+B*log(T)+C));
}

double collision_integral_curvefit22(long s, long r, double T){
  double A,B,C,D;
  
  T=min(max(1000.0,T),30000.0);
  switch (smap[s]){
    case SMAP_eminus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_O: if(T >= 1000.0 && T< 9000.0) {A=0.0164; B=-0.2431; C=1.1231; D=-1.5561;}
          else if(T >= 9000.0 && T<= 30000.0) {A=-0.2027; B=5.6428; C=-51.5646; D=155.6091;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_N: A=0.0; B=0.0; C=0.0; D=1.6094; break;
        case SMAP_NO: if(T >= 1000.0 && T< 8000.0) {A=-0.2202; B=5.2265; C=-40.5659; D=104.7126;}
          else if(T >= 8000.0 && T<= 30000.0) {A=-0.2871; B=8.3757; C=-81.3787; D=265.6292;}
          else fatal_error("Temperature of %.15E K is out of range of validitity (1000-30000 K) for species pair %ld,%ld.",s,r);
          break;
        case SMAP_O2: if(T >= 1000.0 && T< 9000.0) {A=0.0241; B=-0.3467; C=1.3887; D=-0.0110;}
          else if(T >= 9000.0 && T<= 30000.0) {A=0.0025; B=-0.0742; C=0.7235; D=-0.2116;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_N2: A=0.1147; B=-2.8945; C=24.5080; D=-67.3691; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_Oplus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4235; D=6.7787; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_Nplus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4146; D=6.9078; break;
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_O2plus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_N2plus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_NOplus:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=-2.0; D=24.3061; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-2.0; D=24.3602; break;
        case SMAP_O: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NO: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_O:
      switch (smap[r]){
        case SMAP_eminus: if(T >= 1000.0 && T< 9000.0) {A=0.0164; B=-0.2431; C=1.1231; D=-1.5561;}
          else if(T >= 9000.0 && T<= 30000.0) {A=-0.2027; B=5.6428; C=-51.5646; D=155.6091;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4235; D=6.7787; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O: A=0.0; B=-0.0207; C=0.0780; D=3.5658; break;
        case SMAP_N: A=0.0; B=0.0065; C=-0.4467; D=6.0426; break;
        case SMAP_NO: A=0.0; B=-0.0203; C=0.0730; D=3.8818; break;
        case SMAP_O2: A=0.0; B=-0.0247; C=0.1783; D=3.2517; break;
        case SMAP_N2: A=0.0; B=-0.0169; C=-0.0143; D=4.4195; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_N:
      switch (smap[r]){
        case SMAP_eminus: A=0.0; B=0.0; C=0.0; D=1.6094; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4146; D=6.9078; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O: A=0.0; B=0.0065; C=-0.4467; D=6.0426; break;
        case SMAP_N: A=0.0; B=-0.0118; C=-0.0960; D=4.3252; break;
        case SMAP_NO: A=0.0; B=-0.0196; C=0.0478; D=4.0321; break;
        case SMAP_O2: A=0.0; B=-0.0203; C=0.0730; D=3.8818; break;
        case SMAP_N2: A=0.0; B=-0.0190; C=0.0239; D=4.1782; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_NO:
      switch (smap[r]){
        case SMAP_eminus: if(T >= 1000.0 && T< 8000.0) {A=-0.2202; B=5.2265; C=-40.5659; D=104.7126;}
          else if(T >= 8000.0 && T<= 30000.0) {A=-0.2871; B=8.3757; C=-81.3787; D=265.6292;}
          else fatal_error("Temperature of %.15E K is out of range of validitity (1000-30000 K) for species pair %ld,%ld.",s,r);
          break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O: A=0.0; B=-0.0203; C=0.0730; D=3.8818; break;
        case SMAP_N: A=0.0; B=-0.0196; C=0.0478; D=4.0321; break;
        case SMAP_NO: A=0.0; B=-0.0453; C=0.5624; D=1.7669; break;
        case SMAP_O2: A=0.0; B=-0.0522; C=0.7045; D=1.0738; break;
        case SMAP_N2: A=0.0; B=-0.0385; C=0.4226; D=2.4507; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_O2:
      switch (smap[r]){
        case SMAP_eminus: if(T >= 1000.0 && T< 9000.0) {A=0.0241; B=-0.3467; C=1.3887; D=-0.0110;}
          else if(T >= 9000.0 && T<= 30000.0) {A=0.0025; B=-0.0742; C=0.7235; D=-0.2116;}
          else fatal_error("Temperature is out of range of validitity for species pair %ld,%ld.",s,r);
          break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O: A=0.0; B=-0.0247; C=0.1783; D=3.2517; break;
        case SMAP_N: A=0.0; B=-0.0203; C=0.0730; D=3.8818; break;
        case SMAP_NO: A=0.0; B=-0.0522; C=0.7045; D=1.0738; break;
        case SMAP_O2: A=0.0; B=-0.0485; C=0.6475; D=1.2607; break;
        case SMAP_N2: A=0.0; B=-0.0558; C=0.7590; D=0.8995; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_N2:
      switch (smap[r]){
        case SMAP_eminus: A=0.1147; B=-2.8945; C=24.5080; D=-67.3691; break;
        case SMAP_Oplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_Nplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_N2plus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_NOplus: A=0.0; B=0.0; C=-0.4000; D=6.7760; break;
        case SMAP_O: A=0.0; B=-0.0169; C=-0.0143; D=4.4195; break;
        case SMAP_N: A=0.0; B=-0.0190; C=0.0239; D=4.1782; break;
        case SMAP_NO: A=0.0; B=-0.0385; C=0.4226; D=2.4507; break;
        case SMAP_O2: A=0.0; B=-0.0558; C=0.7590; D=0.8995; break;
        case SMAP_N2: A=0.0; B=-0.0203; C=0.0683; D=4.0900; break;
        default: fatal_error("No collision integral data for species pair %ld,%ld.",s,r);
      }
    break;    
    
    default:
      fatal_error("spec %ld does not have collision integral data.",s);
  }
  
  return(exp(D)*pow(T,A*log(T)*log(T)+B*log(T)+C));
}

double collision_crosssection_ratio(long s, long r, double T){
  double A,B,C;
  
  switch (smap[s]){
    case SMAP_eminus:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_N:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_NO:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_O2:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_N2:  A=0.0; B=0.0; C=0.0; break;
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_Oplus:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O:  A=0.0002; B=0.0; C=0.0549; break;
        case SMAP_N:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NO:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2:  A=0.0; B=0.0; C=0.1933; break;
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_Nplus:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N:  A=0.0002; B=0.0002; C=0.0537; break;
        case SMAP_NO:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2:  A=0.0; B=0.0; C=0.1933; break;
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_O2plus:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NO:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2:  A=0.0; B=0.0; C=0.1933; break;
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_N2plus:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NO:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2:  A=0.0; B=0.0; C=0.1933; break;
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_NOplus:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.4463; break;
        case SMAP_O:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NO:  A=0.0003; B=-0.0006; C=0.0632; break;
        case SMAP_O2:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2:  A=0.0; B=0.0; C=0.1933; break;
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_O:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_Oplus:  A=0.0002; B=0.0; C=0.0549; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O:  A=0.0002; B=0.0; C=0.0549; break;
        case SMAP_N:  A=0.0147; B=-0.2628; C=1.2943; break;
        case SMAP_NO:  A=0.0033; B=-0.0366; C=0.2332; break;
        case SMAP_O2:  A=0.0024; B=-0.0245; C=0.1808; break;
        case SMAP_N2:  A=0.0042; B=-0.0471; C=0.2747; break;   
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;
    case SMAP_N:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_Nplus:  A=0.0002; B=0.0002; C=0.0537; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O:  A=0.0147; B=-0.2628; C=1.2943; break;
        case SMAP_N:  A=0.0002; B=0.0002; C=0.0537; break;    
        case SMAP_NO:  A=0.0038; B=-0.0425; C=0.2574; break;    
        case SMAP_O2:  A=0.0033; B=-0.0366; C=0.2332; break;    
        case SMAP_N2:  A=0.0043; B=-0.0494; C=0.2850; break;    
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_NO:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NOplus:  A=0.0003; B=-0.0006; C=0.0632; break;
        case SMAP_O:  A=0.0033; B=-0.0366; C=0.2332; break;
        case SMAP_N:  A=0.0038; B=-0.0425; C=0.2574; break;    
        case SMAP_NO:  A=-0.0027; B=0.0700; C=-0.2553; break;    
        case SMAP_O2:  A=-0.0010; B=0.0410; C=-0.1312; break;    
        case SMAP_N2:  A=-0.0045; B=0.1010; C=-0.3872; break;    
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_O2:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O:  A=0.0024; B=-0.0245; C=0.1808; break;
        case SMAP_N:  A=0.0033; B=-0.0366; C=0.2332; break;    
        case SMAP_NO:  A=-0.0010; B=0.0410; C=-0.1312; break;    
        case SMAP_O2:  A=0.0001; B=0.0181; C=-0.0306; break;    
        case SMAP_N2:  A=-0.0019; B=0.0602; C=-0.2175; break;   
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;    
    case SMAP_N2:
      switch (smap[r]){
        case SMAP_eminus:  A=0.0; B=0.0; C=0.0; break;
        case SMAP_Oplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_Nplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_N2plus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_NOplus:  A=0.0; B=0.0; C=0.1933; break;
        case SMAP_O:  A=0.0042; B=-0.0471; C=0.2747; break;
        case SMAP_N:  A=0.0043; B=-0.0494; C=0.2850; break;    
        case SMAP_NO:  A=-0.0045; B=0.1010; C=-0.3872; break;    
        case SMAP_O2:  A=-0.0019; B=0.0602; C=-0.2175; break;    
        case SMAP_N2:  A=-0.0073; B=0.1444; C=-0.5625; break;    
        default: fatal_error("No collision cross section ratio data for species pair %ld,%ld.",s,r);
      }
    break;    
    
    default:
      fatal_error("spec %ld does not have collision cross section ratio data.",s);
  }
  return(exp(C)*pow(T,A*log(T)+B));
}

static double _cpk_from_T_GUPTAYOS(long k, double T){
  double a1,a2,a3,a4,a5;

  T=min(max(300.0,T),30000.0);
  a1=0.0; a2=0.0; a3=0.0; a4=0.0; a5=0.0; 

  switch (smap[k]){
    case SMAP_eminus: 
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.25E1; a2=0.0E0; a3=0.0E0; a4=0.0E0; a5=0.0E0;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.25E1; a2=0.0E0; a3=0.0E0; a4=0.0E0; a5=0.0E0;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.25E1; a2=0.0E0; a3=0.0E0; a4=0.0E0; a5=0.0E0;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.25E1; a2=0.0E0; a3=0.0E0; a4=0.0E0; a5=0.0E0;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.25E1; a2=0.0E0; a3=0.0E0; a4=0.0E0; a5=0.0E0;}
    break;
    case SMAP_Oplus:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.24985E1; a2=0.11411E-4; a3=-0.29761E-7; a4=0.32247E-10; a5=-0.12376E-13;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.25060E1; a2=-0.14464E-4; a3=0.12446E-7; a4=-0.46858E-11; a5=0.65549E-15;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.29440E1; a2=-0.41080E-3; a3=0.91560E-7; a4=-0.58480E-11; a5=0.11900E-15;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.12784E1; a2=0.40866E-3; a3=-0.21731E-7; a4=0.33252E-12; a5=0.63160E-18;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.12889E1; a2=0.43343E-3; a3=-0.26758E-7; a4=0.62159E-12; a5=-0.45131E-17;}
    break;
    case SMAP_Nplus:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.27270E1; a2=-0.28200E-3; a3=0.11050E-6; a4=-0.15510E-10; a5=0.78470E-15;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.27270E1; a2=-0.28200E-3; a3=0.11050E-6; a4=-0.15510E-10; a5=0.78470E-15;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.24990E1; a2=-0.37250E-5; a3=0.11470E-7; a4=-0.11020E-11; a5=0.30780E-16;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.23856E1; a2=0.83495E-4; a3=-0.58815E-8; a4=0.18850E-12; a5=-0.16120E-17;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.22286E1; a2=0.12458E-3; a3=-0.87636E-8; a4=0.26204E-12; a5=-0.21674E-17;}
    break;
    case SMAP_O2plus:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.32430E1; a2=0.11740E-2; a3=-0.39000E-6; a4=0.54370E-10; a5=-0.23920E-14;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.32430E1; a2=0.11740E-2; a3=-0.39000E-6; a4=0.54370E-10; a5=-0.23920E-14;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.51690E1; a2=-0.86200E-3; a3=0.20410E-6; a4=-0.13000E-10; a5=0.24940E-15;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=-0.28017E0; a2=0.16674E-2; a3=-0.12107E-6; a4=0.32113E-11; a5=-0.28349E-16;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.20445E1; a2=0.10313E-2; a3=-0.74046E-7; a4=0.19257E-11; a5=-0.17461E-16;}
    break;
    case SMAP_N2plus:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.35498E1; a2=-0.60810E-3; a3=0.14690E-5; a4=-0.65091E-10; a5=-0.35649E-12;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.33970E1; a2=0.45250E-3; a3=0.12720E-6; a4=-0.38790E-10; a5=0.24590E-14;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.33780E1; a2=0.86290E-3; a3=-0.12760E-6; a4=0.80870E-11; a5=-0.18800E-15;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.43942E1; a2=0.18868E-3; a3=-0.71272E-8; a4=-0.17511E-12; a5=0.67176E-17;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.39493E1; a2=0.36795E-3; a3=-0.26910E-7; a4=0.67110E-12; a5=-0.58244E-17;}
    break;
    case SMAP_NOplus:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.35294E1; a2=-0.30342E-3; a3=0.38544E-6; a4=0.10519E-8; a5=-0.72777E-12;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.32152E1; a2=0.99742E-3; a3=-0.29030E-6; a4=0.36925E-10; a5=-0.15994E-14;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.26896E1; a2=0.13796E-2; a3=-0.33985E-6; a4=0.33776E-10; a5=-0.10427E-14;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.59346E1; a2=-0.13178E-2; a3=0.23297E-6; a4=-0.11733E-10; a5=0.18402E-15;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=-0.51595E1; a2=0.26290E-2; a3=-0.16254E-6; a4=0.39381E-11; a5=-0.34311E-16;}
    break;
    case SMAP_O:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.28236E1; a2=-0.89478E-3; a3=0.83060E-6; a4=-0.16837E-9; a5=-0.73205E-13;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.25421E1; a2=-0.27551E-4; a3=-0.31028E-8; a4=0.45511E-11; a5=-0.43681E-15;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.25460E1; a2=-0.59520E-4; a3=0.27010E-7; a4=-0.27980E-11; a5=0.93800E-16;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=-0.97871E-2; a2=0.12450E-2; a3=-0.16154E-6; a4=0.80380E-11; a5=-0.12624E-15;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.16428E2; a2=-0.39313E-2; a3=0.29840E-6; a4=-0.81613E-11; a5=0.75004E-16;}
    break;
    case SMAP_N:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.25031E1; a2=-0.21800E-4; a3=0.54205E-7; a4=-0.56476E-10; a5=0.20999E-13;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.24820E1; a2=0.69258E-4; a3=-0.63065E-7; a4=0.18387E-10; a5=-0.11747E-14;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.27480E1; a2=-0.39090E-3; a3=0.13380E-6; a4=-0.11910E-10; a5=0.33690E-15;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=-0.12280E1; a2=0.19268E-2; a3=-0.24370E-6; a4=0.12193E-10; a5=-0.19918E-15;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.15520E2; a2=-0.38858E-2; a3=0.32288E-6; a4=-0.96053E-11; a5=0.95472E-16;}
    break;    
    case SMAP_NO:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.35887E1; a2=-0.12479E-2; a3=0.39786E-5; a4=-0.28651E-8; a5=0.63015E-12;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.32047E1; a2=0.12705E-2; a3=-0.46603E-6; a4=0.75007E-10; a5=-0.42314E-14;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.38543E1; a2=0.23409E-3; a3=-0.21354E-7; a4=0.16689E-11; a5=-0.49070E-16;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.43309E1; a2=-0.58086E-4; a3=0.28059E-7; a4=-0.15694E-11; a5=0.24104E-16;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.23507E1; a2=0.58643E-3; a3=-0.31316E-7; a4=0.60495E-12; a5=-0.40557E-17;}
    break;    
    case SMAP_O2:
      if(T>=300.0e0 && T<1000.0e0)          {a1=0.36146E1; a2=-0.18598E-2; a3=0.70814E-5;  a4=-0.68070E-8;  a5=0.21628E-11;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.35949E1; a2=0.75213E-3;  a3=-0.18732E-6; a4=0.27913E-10;  a5=-0.15774E-14;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.38599E1; a2=0.32510E-3;  a3=-0.92131E-8; a4=-0.78684E-12; a5=0.29426E-16;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.34867E1; a2=0.52384E-3;  a3=-0.39123E-7; a4=0.10094E-11;  a5=-0.88718E-17;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=0.39620E1; a2=0.39446E-3;  a3=-0.29506E-7; a4=0.73975E-12;  a5=-0.64209E-17;}
    break;    
    case SMAP_N2:
      if(T>=300.0e0 && T<1000.0e0)          {a1=3.6748E+00; a2=-1.2081E-03; a3=2.3240E-06;  a4=-6.3218E-10;  a5=-2.2577E-13;}
      else if(T>=1000.0e0 && T<6000.0e0)    {a1=0.32125E1;  a2=0.10137E-2;  a3=-0.30467E-6; a4=0.41091E-10;  a5=-0.2017E-14;}
      else if(T>=6000.0e0 && T<15000.0e0)   {a1=0.31811E1;  a2=0.89745E-3;  a3=-0.20216E-6; a4=0.18266E-10;  a5=-0.50334E-15;}
      else if(T>=15000.0e0 && T<25000.0e0)  {a1=0.96377E1;  a2=-0.25728E-2; a3=0.33020E-6;  a4=-0.14315E-10; a5=0.20333E-15;}
      else if(T>=25000.0e0 && T<=30000.0e0) {a1=-0.51681E1; a2=0.23337E-2;  a3=-0.12953E-6; a4=0.27872E-11;  a5=-0.21360E-16;}
    break;    
    default:
      fatal_error("No specific heat data for species pair %ld.",k);
  }
 
  
  return(Runiv*(a1+a2*T+a3*T*T+a4*T*T*T+a5*T*T*T*T));
}



double EXM_matrix_determinant(EXM_mat_t mat){//improve this, triangular matrix is enough
  long row,col,row2;
  double fact,det;
  EXM_mat_t mattmp;
  EXM_init_matrix(&mattmp,mat.glm.numrow,mat.glm.numcol);
  if (mat.glm.numrow==mat.glm.numcol) {
    for (row=0; row<mat.glm.numrow; row++){
      for (col=0; col<mat.glm.numcol; col++){
        mattmp.cont[EXM_aim(mat.glm,row,col)]=mat.cont[EXM_aim(mat.glm,row,col)];
      }
    }
// make the non-diagonal elements zero for mattmp 
    for (row=0; row<mat.glm.numrow; row++){
      for (row2=0; row2<mat.glm.numrow; row2++){
        if (row2!=row) {
          if (mattmp.cont[EXM_aim(mat.glm,row,row)]==0.0){
            EXM_free_matrix(&mattmp);
            return(0.0e0);
          }     
          fact=-mattmp.cont[EXM_aim(mat.glm,row2,row)]/mattmp.cont[EXM_aim(mat.glm,row,row)];
          for (col=0; col<mat.glm.numcol; col++){
            mattmp.cont[EXM_aim(mat.glm,row2,col)]+=fact*mattmp.cont[EXM_aim(mat.glm,row,col)];
          }
        }
      }    
    }
   //EXM_display_matrix(mattmp);
  } else {
    printf("determinant cannot be calculated\n");
    printf("number of rows not equal to number of columns \n");
  }
  det=1.0e0;
  for (row=0; row<mat.glm.numrow; row++){
    det*=mattmp.cont[EXM_aim(mat.glm,row,row)];
  }
  EXM_free_matrix(&mattmp);
  return(det);
}

static double _epc(double rhoe, double Te){
  double Pe;
  double epc;
  
  Pe = rhoe*calR*Te/_calM(speceminus)/101325.0e0; //electron pressure, atm
  /*electron pressure correction for the collision integrals of ionic species*/ /*Eq(24b)*/
  switch (EPC){
    case EPC_GUPTAYOS:
      epc=0.5*log(2.09E-2*1E-12*Te*Te*Te*Te/Pe+1.52*pow(Te*Te*Te*Te*1E-12/Pe,2.0/3.0));
    break;
    case EPC_RAIZER:
      // Raizer Gas Discharge Physics, page 14
      epc=13.57+1.5*log(Te/11600.0)/log(10.0)-0.5*log(rhoe/_m(speceminus)/1e6)/log(10.0);
    break;
    case EPC_NRL:
      // NRL Plasma Formulary, page 34
      epc=23.0-log(sqrt(rhoe/_m(speceminus)/1e6)*pow(Te/11600.0,-1.5));
    break;
    case EPC_NONE:
      epc=1.0;
    break;
    default:
      fatal_error("EPC set to invalid value.");
  }
  return(epc);
}


static void find_aij_Ai_Delta1_Delta2(spec_t rhok, double T, double epc, spec2_t aij, spec_t Ai, spec2_t Delta1, spec2_t Delta2){
  double N,Bstar;
  long spec,i,j,l;
  spec2_t Bij;
  spec_t chik,M;
  
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  for (spec=0; spec<ns; spec++) chik[spec]=rhok[spec]/_m(spec)/N;
  for (spec=0; spec<ns; spec++) M[spec]=_calM(spec)*1000.0; // g/g-mole

  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      Delta1[i][j]=collision_integral_curvefit11(i,j,T)*(8.0e0/3.0e0)*sqrt((2.0e0*M[i]*M[j])/(pi*Runiv*T*(M[i] + M[j])))*1.546e-20; // 1e-20 m^2 = 1 Angstrom^2
      Delta2[i][j]=collision_integral_curvefit22(i,j,T)*(0.2*16.0e0)*sqrt((2.0e0*M[i]*M[j])/(pi*Runiv*T*(M[i] + M[j])))*1.546e-20;
      if(_Charge_number(i)!=0 && _Charge_number(j)!=0) {Delta1[i][j]*=epc; Delta2[i][j]*=epc;}
      
      Bstar=collision_crosssection_ratio(i,j,T);
      aij[i][j]=4.184E7*2.0*M[i]*M[j]*((0.5*33.0-0.2*18.0*Bstar)*Delta1[i][j]-4.0*Delta2[i][j])/(15.0*kBol*(M[i]+M[j])*(M[i]+M[j]));
      Bij[i][j]=4.184E7*2.0*(8.0*M[i]*M[j]*Delta2[i][j]+(M[i]-M[j])*(9.0*M[i]-7.5*M[j]+3.6*Bstar*M[j])*Delta1[i][j])/(15.0*kBol*(M[i]+M[j])*(M[i]+M[j]));
    }
  }
  for (i=0; i<ns; i++){
    Ai[i]=0.0;
    for (l=0; l<ns; l++){
      Ai[i]+=chik[l]*Bij[i][l];
    }
  }

}


// Eq(30), approximation to method 1, yields results close to method 1 for air mixture
double _kappa_from_rhok_T_Te_METHOD2(spec_t rhok, double T, double Te){
  long spec,i,j;
  spec_t chik,  Ai;
  spec2_t Delta1,Delta2;
  double N, den, kappa_tr, kappa_int;
  double num,aav;
  spec2_t aij;

  find_aij_Ai_Delta1_Delta2(rhok, T, _epc(rhok[speceminus], Te), aij, Ai, Delta1, Delta2);
  
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  for (spec=0; spec<ns; spec++) chik[spec]=rhok[spec]/_m(spec)/N;

  kappa_tr     = 0.0e0;
  kappa_int    = 0.0e0;
  num=0.0;
  den=0.0;
  
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      num+=chik[i]*chik[j]*(1.0/Ai[i]-1.0/Ai[j])*(1.0/Ai[i]-1.0/Ai[j])*aij[i][j];
      den+=num/aij[i][j];
    }
  }
  aav=num/den;
  num=0.0;
  for (i=0; i<ns; i++){
    num+=chik[i]/(Ai[i]+aav);
  }
  kappa_tr=num/(1.0-aav*num);
  
  /*kappa_int=kappa_rot+kappa_vib+kappa_el*/
  den=0.0;
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      den+=chik[j]*Delta1[i][j];
    }
    kappa_int+=(_cpk_from_T_GUPTAYOS(i,T)/Runiv-2.5)*chik[i]/den;
    den=0.0;
  }
  kappa_int*=2.3901E-8*kBol;
  
  return((kappa_tr+kappa_int)*418.4); // cal/cm-s-K to W/m-K
}





// Eq(27), first Chapman-Enskog approximation, bug: off by factor of -1 : needs fixing
double _kappa_from_rhok_T_Te_METHOD1(spec_t rhok, double T, double Te){
  long spec,i,j,l;
  spec_t chik,  Ai;
  spec2_t Delta1,Delta2;
  double N, den, kappa_tr, kappa_int;
  
  double Aij[ns][ns];
  EXM_mat_t tmp1,tmp2;
  double tmp3=0.0;
  spec2_t aij;

  find_aij_Ai_Delta1_Delta2(rhok, T, _epc(rhok[speceminus], Te), aij, Ai, Delta1, Delta2);
  
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  for (spec=0; spec<ns; spec++) chik[spec]=rhok[spec]/_m(spec)/N;

  kappa_tr     = 0.0e0;
  kappa_int    = 0.0e0;

  EXM_init_matrix(&tmp1,ns+1,ns+1);
  EXM_init_matrix(&tmp2,ns,ns);
  EXM_mat_t *numtr =&tmp1;
  EXM_mat_t *dentr =&tmp2;
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      for (l=0; l<ns; l++){
        tmp3+=chik[i]*chik[l]*aij[i][l];
      }
      Aij[i][j]=-chik[i]*chik[j]*aij[i][j]+krodelta(i,j)*(chik[i]*Ai[i]+tmp3);
      tmp3=0.0;
    }
  }
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      dentr->cont[EXM_aim(dentr->glm,i,j)]=Aij[i][j];
      numtr->cont[EXM_aim(numtr->glm,i,j)]=Aij[i][j];
    }
  }
  for (i=0; i<ns; i++){ numtr->cont[EXM_aim(numtr->glm,i,ns)]=chik[i]; }
  for (j=0; j<ns; j++){ numtr->cont[EXM_aim(numtr->glm,ns,j)]=chik[j]; }
  kappa_tr=EXM_matrix_determinant(tmp1)/EXM_matrix_determinant(tmp2);
  kappa_tr=fabs(kappa_tr); //bug: off by factor of -1 : needs fixing
  //EXM_display_matrix(tmp1);

  /*kappa_int=kappa_rot+kappa_vib+kappa_el*/
  den=0.0;
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      den+=chik[j]*Delta1[i][j];
    }
    kappa_int+=(_cpk_from_T_GUPTAYOS(i,T)/Runiv-2.5)*chik[i]/den;
    den=0.0;
  }
  kappa_int*=2.3901E-8*kBol;
  
  return((kappa_tr+kappa_int)*418.4); // cal/cm-s-K to W/m-K
}




// Eq(40a), approximation to method 1 ??? or method 2??
double _kappa_from_rhok_T_Te_METHOD3(spec_t rhok, double T, double Te){
  long spec,i,j;
  spec_t chik,  Ai;
  spec2_t Delta1,Delta2;
  double N, den, kappa_tr, kappa_int;
  double asr;
  spec2_t aij;

  find_aij_Ai_Delta1_Delta2(rhok, T, _epc(rhok[speceminus], Te), aij, Ai, Delta1, Delta2);
  
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  for (spec=0; spec<ns; spec++) chik[spec]=rhok[spec]/_m(spec)/N;

  kappa_tr     = 0.0e0;
  kappa_int    = 0.0e0;

  den=0.0;
  kappa_tr=0.0;
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      asr  = 1.0e0 + (1.0e0-(_m(i)/_m(j))*(0.45e0 - 2.54e0*(_m(i)/_m(j))))/((1.0e0+_m(i)/_m(j))*(1.0e0+_m(i)/_m(j)));
      den+=asr*chik[j]*Delta2[i][j];
    }
    kappa_tr+=chik[i]/den;
    den=0.0;
  }
  kappa_tr*=2.3901E-8*0.25*15.0*kBol;

  /*kappa_int=kappa_rot+kappa_vib+kappa_el*/
  den=0.0;
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      den+=chik[j]*Delta1[i][j];
    }
    kappa_int+=(_cpk_from_T_GUPTAYOS(i,T)/Runiv-2.5)*chik[i]/den;
    den=0.0;
  }
  kappa_int*=2.3901E-8*kBol;
    
  return((kappa_tr+kappa_int)*418.4); // cal/cm-s-K to W/m-K
}




// Eq(49) : 5-temperature model, thermal non-equilibirum
double _kappa_from_rhok_T_Te_METHOD4(spec_t rhok, double T, double Te){
  long spec,i,j;
  spec_t chik,  Ai, M;
  spec2_t Delta1,Delta2;
  double N, den, kappa_tr, kappa_int;
  double asr,kappa_e,kappa_vib,epc;
  double Delta1ie[ns],Delta2ie[ns];
  spec2_t aij;
  epc=_epc(rhok[speceminus], Te);
  find_aij_Ai_Delta1_Delta2(rhok, T, epc, aij, Ai, Delta1, Delta2);
  for (spec=0; spec<ns; spec++) M[spec]=_calM(spec)*1000.0; // g/g-mole
  
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  for (spec=0; spec<ns; spec++) chik[spec]=rhok[spec]/_m(spec)/N;

  kappa_tr     = 0.0e0;
  kappa_int    = 0.0e0;

  /*kappa_int=kappa_rot+kappa_vib+kappa_el*/
  den=0.0;
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      den+=chik[j]*Delta1[i][j];
    }
    kappa_int+=(_cpk_from_T_GUPTAYOS(i,T)/Runiv-2.5)*chik[i]/den;
    den=0.0;
  }
  kappa_int*=2.3901E-8*kBol;
  
  for (i=0; i<ns; i++){
    j=speceminus;
    Delta1ie[i]=collision_integral_curvefit11(i,j,Te)*(8.0e0/3.0e0)*sqrt((2.0e0*M[i]*M[j])/(pi*Runiv*Te*(M[i] + M[j])))*1.546e-20;
    Delta2ie[i]=collision_integral_curvefit22(i,j,Te)*(0.2*16.0e0)*sqrt((2.0e0*M[i]*M[j])/(pi*Runiv*Te*(M[i] + M[j])))*1.546e-20;
    if(_Charge_number(i)!=0) {Delta1ie[i]*=epc;Delta2ie[i]*=epc;}
  }
  den=0.0;
  kappa_tr=0.0;
  kappa_e=0.0;
  kappa_vib=0.0;
  
  for (i=0; i<ns; i++){
    if (speciestype[i]!=SPECIES_ELECTRON) { 
      for (j=0; j<ns; j++){
        if (speciestype[j]!=SPECIES_ELECTRON) { 
          asr  = 1.0e0 + (1.0e0-(_m(i)/_m(j))*(0.45e0 - 2.54e0*(_m(i)/_m(j))))/((1.0e0+_m(i)/_m(j))*(1.0e0+_m(i)/_m(j)));
          den+=asr*chik[j]*Delta2[i][j]+3.54*chik[speceminus]*Delta2ie[i];
        }
      }
      kappa_tr+=chik[i]/den;
      den=0.0;
    }
  }
  kappa_tr*=2.3901E-8*0.25*15.0*kBol;
  
  den=0.0;
  for (i=0; i<ns; i++){
    if (_numatoms(i)>1) { //molecular species 
      for (j=0; j<ns; j++){
        if (speciestype[j]!=SPECIES_ELECTRON) { 
          asr  = 1.0e0 + (1.0e0-(_m(i)/_m(j))*(0.45e0 - 2.54e0*(_m(i)/_m(j))))/((1.0e0+_m(i)/_m(j))*(1.0e0+_m(i)/_m(j)));
          den+=chik[j]*Delta1[i][j]+chik[speceminus]*Delta1ie[i];
        }
      }
      kappa_vib+=chik[i]/den;
      den=0.0;
    }
  }
  kappa_vib*=2.3901E-8*kBol;
  
  den=0.0;
  for (j=0; j<ns; j++){
    if (speciestype[j]!=SPECIES_ELECTRON) { 
      asr  = 1.0e0 + (1.0e0-(_m(i)/_m(j))*(0.45e0 - 2.54e0*(_m(i)/_m(j))))/((1.0e0+_m(i)/_m(j))*(1.0e0+_m(i)/_m(j)));
      den+=1.45*chik[j]*Delta2[speceminus][j]+chik[speceminus]*Delta2ie[speceminus];
    }
  }
  kappa_e=2.3901E-8*kBol*chik[speceminus]/den;
  kappa_int=kappa_vib+kappa_e; //just for calculation purposes
  
  return((kappa_tr+kappa_int)*418.4); // cal/cm-s-K to W/m-K
}



double _kappa_from_rhok_T_Te(spec_t rhok, double T, double Te){
  double kappa;
  switch (METHOD){
    case (METHOD1):
      kappa=_kappa_from_rhok_T_Te_METHOD1(rhok, T, Te);
    break;
    case (METHOD2):
      kappa=_kappa_from_rhok_T_Te_METHOD2(rhok, T, Te);
    break;
    case (METHOD3):
      kappa=_kappa_from_rhok_T_Te_METHOD3(rhok, T, Te);
    break;
    case (METHOD4):
      kappa=_kappa_from_rhok_T_Te_METHOD4(rhok, T, Te);
    break;
    default:
      fatal_error("METHOD can not be set to %ld",METHOD);
  }
  return(kappa);
}





double _kappan_from_rhok_T_Te(spec_t rhok, double T, double Te){
  return(1.0);

}

double _eta_from_rhok_T_Te(spec_t rhok, double T, double Te){
  return(1.0);

}

double _etan_from_rhok_T_Te(spec_t rhok, double T, double Te){
  return(1.0);
}

void find_nuk_from_rhok_T_Te(spec_t rhok, double T, double Te, spec_t nuk){

}

double _muk_from_rhok_T_Te_Ek(spec_t rhok, double T, double Te, double Ek, long k){
  long spec,i,j;
  spec_t chik,  Ai;
  spec2_t Delta1,Delta2;
  double N, den, kappak,Tk;
  double num,aav,muk,cpk;
  spec2_t aij;

  if (k==speceminus) Tk=Te; else Tk=T;
  find_aij_Ai_Delta1_Delta2(rhok, Tk, _epc(rhok[speceminus], Te), aij, Ai, Delta1, Delta2);
  
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  for (spec=0; spec<ns; spec++) chik[spec]=rhok[spec]/_m(spec)/N; 
  num=0.0;
  den=0.0;
  for (i=0; i<ns; i++){
    for (j=0; j<ns; j++){
      num+=chik[i]*chik[j]*(1.0/Ai[i]-1.0/Ai[j])*(1.0/Ai[i]-1.0/Ai[j])*aij[i][j];
      den+=num/aij[i][j];
    }
  }
  aav=num/den;
  num=0.0;
  for (i=0; i<ns; i++){
    num+=chik[i]/(Ai[i]+aav);
  }
  
  kappak=chik[k]/(Ai[k]+aav)/(1.0-aav*num)*418.4;

  cpk=_cpk_from_T_equilibrium(k,Tk);
  muk=kappak/(cpk*rhok[k]*kB*Tk)*fabs(_C(k));
  
  return(muk);
}


void find_dmuk_from_rhok_Tk_Ek(spec_t rhok, double Tk, double Ek, long k, double *dmukdTk, spec_t dmukdrhok){
  long spec;
  *dmukdTk=0.0;
  for (spec=0; spec<ns; spec++) dmukdrhok[spec]=0.0;
}



