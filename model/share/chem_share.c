// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2018-2020 Bernard Parent
Copyright 2020 Ajjay Omprakas

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

#include <model/share/chem_share.h>
#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>

#define Kcmin 1.0e-40



/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   W in kg m^(-3) s^(-1)
*/
void add_to_W_2r2p ( int specR1, int specR2,
                     int specP1, int specP2,
                     double kf, spec_t N, spec_t W){
  double Wp;
  Wp=kf*N[specR1]*N[specR2];
  W[specR1]-=Wp/ calA * _calM(specR1) * 1.0e6;
  W[specR2]-=Wp/ calA * _calM(specR2) * 1.0e6;
  W[specP1]+=Wp/ calA * _calM(specP1) * 1.0e6;
  W[specP2]+=Wp/ calA * _calM(specP2) * 1.0e6;                            
}

/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   W in kg m^(-3) s^(-1)
*/
void add_to_W_2r1p ( int specR1, int specR2,
                     int specP1, 
                     double kf, spec_t N, spec_t W){
  double Wp;
  Wp=kf*N[specR1]*N[specR2];
  W[specR1]-=Wp/ calA * _calM(specR1) * 1.0e6;
  W[specR2]-=Wp/ calA * _calM(specR2) * 1.0e6;
  W[specP1]+=Wp/ calA * _calM(specP1) * 1.0e6;
}


/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   W in kg m^(-3) s^(-1)
*/
void add_to_W_2r3p ( int specR1, int specR2,
                     int specP1, int specP2, int specP3,
                     double kf, spec_t N, spec_t W){
  double Wp;
  Wp=kf*N[specR1]*N[specR2];
  W[specR1]-=Wp/ calA * _calM(specR1) * 1.0e6;
  W[specR2]-=Wp/ calA * _calM(specR2) * 1.0e6;
  W[specP1]+=Wp/ calA * _calM(specP1) * 1.0e6;
  W[specP2]+=Wp/ calA * _calM(specP2) * 1.0e6;                            
  W[specP3]+=Wp/ calA * _calM(specP3) * 1.0e6;                            
}



/* kf in cm^6 s^(-1) 
   N in cm^(-3)
   W in kg m^(-3) s^(-1)
*/
void add_to_W_3r2p ( int specR1, int specR2, int specR3,
                     int specP1, int specP2,
                     double kf, spec_t N, spec_t W){
  double Wp;
  Wp=kf*N[specR1]*N[specR2]*N[specR3];
  W[specR1]-=Wp/ calA * _calM(specR1) * 1.0e6;
  W[specR2]-=Wp/ calA * _calM(specR2) * 1.0e6;
  W[specR3]-=Wp/ calA * _calM(specR3) * 1.0e6;
  W[specP1]+=Wp/ calA * _calM(specP1) * 1.0e6;
  W[specP2]+=Wp/ calA * _calM(specP2) * 1.0e6;                            
}


/* kf in cm^6 s^(-1) 
   N in cm^(-3)
   W in kg m^(-3) s^(-1)
*/
void add_to_W_3r3p ( int specR1, int specR2, int specR3,
                     int specP1, int specP2, int specP3,
                     double kf, spec_t N, spec_t W){
  double Wp;
  Wp=kf*N[specR1]*N[specR2]*N[specR3];
  W[specR1]-=Wp/ calA * _calM(specR1) * 1.0e6;
  W[specR2]-=Wp/ calA * _calM(specR2) * 1.0e6;
  W[specR3]-=Wp/ calA * _calM(specR3) * 1.0e6;
  W[specP1]+=Wp/ calA * _calM(specP1) * 1.0e6;
  W[specP2]+=Wp/ calA * _calM(specP2) * 1.0e6;                            
  W[specP3]+=Wp/ calA * _calM(specP3) * 1.0e6;                            
}



/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   dkfdT in cm^3 s^(-1) K^(-1)
   dkfdTv in cm^3 s^(-1) K^(-1)
   dkfdTe in cm^3 s^(-1) K^(-1)
   dWdrhok in s^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdTv in kg m^(-3) s^(-1) K^(-1)
   dWdTe in kg m^(-3) s^(-1) K^(-1)
*/
void add_to_dW_2r1p ( int specR1, int specR2, int specP1, double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe){
  
  dWdTe[specR1]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTe[specR2]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTe[specP1]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  
  dWdTv[specR1]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTv[specR2]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTv[specP1]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  
  dWdT[specR1]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdT[specR2]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdT[specP1]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  
  
  dWdrhok[specR1][specR1]-=kf*N[specR2];
  dWdrhok[specR1][specR2]-=kf*N[specR1]* _calM(specR1) / _calM(specR2);

  dWdrhok[specR2][specR1]-=kf*N[specR2]* _calM(specR2) / _calM(specR1);
  dWdrhok[specR2][specR2]-=kf*N[specR1];
  
  dWdrhok[specP1][specR1]+=kf*N[specR2]* _calM(specP1) / _calM(specR1);
  dWdrhok[specP1][specR2]+=kf*N[specR1]* _calM(specP1) / _calM(specR2);

}


/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   dkfdT in cm^3 s^(-1) K^(-1)
   dkfdTv in cm^3 s^(-1) K^(-1)
   dkfdTe in cm^3 s^(-1) K^(-1)
   dWdrhok in s^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdTv in kg m^(-3) s^(-1) K^(-1)
   dWdTe in kg m^(-3) s^(-1) K^(-1)
*/
void add_to_dW_2r2p ( int specR1, int specR2, int specP1, int specP2, double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe){
  
  dWdTe[specR1]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTe[specR2]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTe[specP1]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTe[specP2]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  
  dWdTv[specR1]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTv[specR2]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTv[specP1]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTv[specP2]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  
  dWdT[specR1]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdT[specR2]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdT[specP1]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdT[specP2]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  
  
  dWdrhok[specR1][specR1]-=kf*N[specR2];
  dWdrhok[specR1][specR2]-=kf*N[specR1]* _calM(specR1) / _calM(specR2);

  dWdrhok[specR2][specR1]-=kf*N[specR2]* _calM(specR2) / _calM(specR1);
  dWdrhok[specR2][specR2]-=kf*N[specR1];
  
  dWdrhok[specP1][specR1]+=kf*N[specR2]* _calM(specP1) / _calM(specR1);
  dWdrhok[specP1][specR2]+=kf*N[specR1]* _calM(specP1) / _calM(specR2);

  dWdrhok[specP2][specR1]+=kf*N[specR2]* _calM(specP2) / _calM(specR1);
  dWdrhok[specP2][specR2]+=kf*N[specR1]* _calM(specP2) / _calM(specR2);
}


/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   dkfdT in cm^3 s^(-1) K^(-1)
   dkfdTv in cm^3 s^(-1) K^(-1)
   dkfdTe in cm^3 s^(-1) K^(-1)
   dWdrhok in s^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdTv in kg m^(-3) s^(-1) K^(-1)
   dWdTe in kg m^(-3) s^(-1) K^(-1)
*/
void add_to_dW_2r3p ( int specR1, int specR2, int specP1, int specP2, int specP3, double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe){
  
  dWdTe[specR1]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTe[specR2]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTe[specP1]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTe[specP2]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  dWdTe[specP3]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP3) * 1.0e6;
  
  dWdTv[specR1]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTv[specR2]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTv[specP1]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTv[specP2]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  dWdTv[specP3]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP3) * 1.0e6;
  
  dWdT[specR1]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdT[specR2]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdT[specP1]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdT[specP2]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  dWdT[specP3]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP3) * 1.0e6;
  
  
  dWdrhok[specR1][specR1]-=kf*N[specR2];
  dWdrhok[specR1][specR2]-=kf*N[specR1]* _calM(specR1) / _calM(specR2);

  dWdrhok[specR2][specR1]-=kf*N[specR2]* _calM(specR2) / _calM(specR1);
  dWdrhok[specR2][specR2]-=kf*N[specR1];
  
  dWdrhok[specP1][specR1]+=kf*N[specR2]* _calM(specP1) / _calM(specR1);
  dWdrhok[specP1][specR2]+=kf*N[specR1]* _calM(specP1) / _calM(specR2);

  dWdrhok[specP2][specR1]+=kf*N[specR2]* _calM(specP2) / _calM(specR1);
  dWdrhok[specP2][specR2]+=kf*N[specR1]* _calM(specP2) / _calM(specR2);

  dWdrhok[specP3][specR1]+=kf*N[specR2]* _calM(specP3) / _calM(specR1);
  dWdrhok[specP3][specR2]+=kf*N[specR1]* _calM(specP3) / _calM(specR2);
}



/* kf in cm^6 s^(-1) 
   N in cm^(-3)
   dkfdT in cm^6 s^(-1) K^(-1)
   dkfdTv in cm^6 s^(-1) K^(-1)
   dkfdTe in cm^6 s^(-1) K^(-1)
   dWdrhok in s^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdTv in kg m^(-3) s^(-1) K^(-1)
   dWdTe in kg m^(-3) s^(-1) K^(-1)
*/
void add_to_dW_3r3p ( int specR1, int specR2, int specR3, int specP1, int specP2, int specP3, 
                      double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe){
  
  dWdTe[specR1]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTe[specR2]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTe[specR3]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR3) * 1.0e6;
  dWdTe[specP1]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTe[specP2]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  dWdTe[specP3]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP3) * 1.0e6;
  
  dWdTv[specR1]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTv[specR2]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTv[specR3]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR3) * 1.0e6;
  dWdTv[specP1]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTv[specP2]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  dWdTv[specP3]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP3) * 1.0e6;
  
  dWdT[specR1]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdT[specR2]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdT[specR3]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR3) * 1.0e6;
  dWdT[specP1]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdT[specP2]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  dWdT[specP3]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP3) * 1.0e6;
  
  
  dWdrhok[specR1][specR1]-=kf*N[specR2]*N[specR3];
  dWdrhok[specR1][specR2]-=kf*N[specR1]*N[specR3]* _calM(specR1) / _calM(specR2);
  dWdrhok[specR1][specR3]-=kf*N[specR1]*N[specR2]* _calM(specR1) / _calM(specR3);

  dWdrhok[specR2][specR1]-=kf*N[specR2]*N[specR3]* _calM(specR2) / _calM(specR1);
  dWdrhok[specR2][specR2]-=kf*N[specR1]*N[specR3];
  dWdrhok[specR2][specR3]-=kf*N[specR1]*N[specR2]* _calM(specR2) / _calM(specR3);
  

  dWdrhok[specR3][specR1]-=kf*N[specR2]*N[specR3]* _calM(specR3) / _calM(specR1);
  dWdrhok[specR3][specR2]-=kf*N[specR1]*N[specR3]* _calM(specR3) / _calM(specR2);
  dWdrhok[specR3][specR3]-=kf*N[specR1]*N[specR2];
  
  
  dWdrhok[specP1][specR1]+=kf*N[specR2]*N[specR3]* _calM(specP1) / _calM(specR1);
  dWdrhok[specP1][specR2]+=kf*N[specR1]*N[specR3]* _calM(specP1) / _calM(specR2);
  dWdrhok[specP1][specR3]+=kf*N[specR1]*N[specR2]* _calM(specP1) / _calM(specR3);

  dWdrhok[specP2][specR1]+=kf*N[specR2]*N[specR3]* _calM(specP2) / _calM(specR1);
  dWdrhok[specP2][specR2]+=kf*N[specR1]*N[specR3]* _calM(specP2) / _calM(specR2);
  dWdrhok[specP2][specR3]+=kf*N[specR1]*N[specR2]* _calM(specP2) / _calM(specR3);


  dWdrhok[specP3][specR1]+=kf*N[specR2]*N[specR3]* _calM(specP3) / _calM(specR1);
  dWdrhok[specP3][specR2]+=kf*N[specR1]*N[specR3]* _calM(specP3) / _calM(specR2);
  dWdrhok[specP3][specR3]+=kf*N[specR1]*N[specR2]* _calM(specP3) / _calM(specR3);
}


/* kf in cm^6 s^(-1) 
   N in cm^(-3)
   dkfdT in cm^6 s^(-1) K^(-1)
   dkfdTv in cm^6 s^(-1) K^(-1)
   dkfdTe in cm^6 s^(-1) K^(-1)
   dWdrhok in s^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdTv in kg m^(-3) s^(-1) K^(-1)
   dWdTe in kg m^(-3) s^(-1) K^(-1)
*/
void add_to_dW_3r2p ( int specR1, int specR2, int specR3, int specP1, int specP2,  
                      double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe){
  
  dWdTe[specR1]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTe[specR2]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTe[specR3]-=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specR3) * 1.0e6;
  dWdTe[specP1]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTe[specP2]+=dkfdTe*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  
  dWdTv[specR1]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdTv[specR2]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdTv[specR3]-=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specR3) * 1.0e6;
  dWdTv[specP1]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdTv[specP2]+=dkfdTv*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  
  dWdT[specR1]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR1) * 1.0e6;
  dWdT[specR2]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR2) * 1.0e6;
  dWdT[specR3]-=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specR3) * 1.0e6;
  dWdT[specP1]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP1) * 1.0e6;
  dWdT[specP2]+=dkfdT*N[specR1]*N[specR2]/ calA * _calM(specP2) * 1.0e6;
  
  
  dWdrhok[specR1][specR1]-=kf*N[specR2]*N[specR3];
  dWdrhok[specR1][specR2]-=kf*N[specR1]*N[specR3]* _calM(specR1) / _calM(specR2);
  dWdrhok[specR1][specR3]-=kf*N[specR1]*N[specR2]* _calM(specR1) / _calM(specR3);

  dWdrhok[specR2][specR1]-=kf*N[specR2]*N[specR3]* _calM(specR2) / _calM(specR1);
  dWdrhok[specR2][specR2]-=kf*N[specR1]*N[specR3];
  dWdrhok[specR2][specR3]-=kf*N[specR1]*N[specR2]* _calM(specR2) / _calM(specR3);
  

  dWdrhok[specR3][specR1]-=kf*N[specR2]*N[specR3]* _calM(specR3) / _calM(specR1);
  dWdrhok[specR3][specR2]-=kf*N[specR1]*N[specR3]* _calM(specR3) / _calM(specR2);
  dWdrhok[specR3][specR3]-=kf*N[specR1]*N[specR2];
  
  
  dWdrhok[specP1][specR1]+=kf*N[specR2]*N[specR3]* _calM(specP1) / _calM(specR1);
  dWdrhok[specP1][specR2]+=kf*N[specR1]*N[specR3]* _calM(specP1) / _calM(specR2);
  dWdrhok[specP1][specR3]+=kf*N[specR1]*N[specR2]* _calM(specP1) / _calM(specR3);

  dWdrhok[specP2][specR1]+=kf*N[specR2]*N[specR3]* _calM(specP2) / _calM(specR1);
  dWdrhok[specP2][specR2]+=kf*N[specR1]*N[specR3]* _calM(specP2) / _calM(specR2);
  dWdrhok[specP2][specR3]+=kf*N[specR1]*N[specR2]* _calM(specP2) / _calM(specR3);

}


/* find Gs [J/mole] for species spec at a temperature T [K] */
static double _Gs(long spec, double T){
  double Gs;
  Gs=( _hk_from_T_equilibrium ( spec, T ) - T * _sk_from_T_equilibrium ( spec, T ) ) * _calM ( spec );
  return(Gs);
}


/* find dGsdT [J/moleK] for species spec at a temperature T [K] */ 
static double _dGsdT(long spec, double T){
  double dGsdT; 
  dGsdT = ( _cpk_from_T_equilibrium ( spec, T ) - _sk_from_T_equilibrium ( spec, T ) - T * _dsk_dT_from_T_equilibrium ( spec, T )
       ) * _calM ( spec );
  return(dGsdT);
}


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fwbw_2r2p(int specR1, int specR2,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t W){
  add_to_W_fw_2r2p(specR1, specR2, specP1, specP2, A, n, E, T, X, W);
  add_to_W_bw_2r2p(specR1, specR2, specP1, specP2, A, n, E, T, X, W);
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_2r2p(int specR1, int specR2,
                      int specP1, int specP2,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  sum=kf*X[specR1]*X[specR2];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_bw_2r2p(int specR1, int specR2,
                      int specP1, int specP2,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  kb=kf/Kc; /* cm^3 (mole s)^(-1) */
  sum=-kb*X[specP1]*X[specP2]; /* mole cm^(-3) (s)^(-1) */
  if (exp(-dG0/(calR*T))>Kcmin){
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6;
    W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
  }
} 




/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fwbw_3r2p(int specR1, int specR2, int specR3,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t W){
  add_to_W_fw_3r2p(specR1, specR2, specR3, specP1, specP2, A, n, E, T, X, W);
  add_to_W_bw_3r2p(specR1, specR2, specR3, specP1, specP2, A, n, E, T, X, W);
} 



/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_3r2p(int specR1, int specR2, int specR3,
                      int specP1, int specP2,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  sum=kf*X[specR1]*X[specR2]*X[specR3]; /* mole cm^(-3) (s)^(-1) */
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specR3]-=_calM(specR3)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_bw_3r2p(int specR1, int specR2, int specR3,
                      int specP1, int specP2,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T)
     -_Gs(specR3,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* mole/cm^3 */
  kb=kf/Kc; /* cm^3 (mole s)^(-1) */
  sum=-kb*X[specP1]*X[specP2]; /* mole cm^(-3) (s)^(-1) */
  if (exp(-dG0/(calR*T))>Kcmin){
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specR3]-=_calM(specR3)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6;
    W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
  }
} 



/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fwbw_2r1p(int specR1, int specR2,
                            int specP1,
                            double A, double n, double E, double T, spec_t X, spec_t W){
  add_to_W_fw_2r1p(specR1, specR2, specP1, A, n, E, T, X, W);
  add_to_W_bw_2r1p(specR1, specR2, specP1, A, n, E, T, X, W);
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_2r1p(int specR1, int specR2,
                            int specP1,
                            double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  sum=kf*X[specR1]*X[specR2]; /* mole cm^(-3) (s)^(-1) */
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6; /* kg m^(-3) s^(-1) */
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_bw_2r1p(int specR1, int specR2,
                      int specP1,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* mole/cm^3 */
  kb=kf/Kc; /* cm^3 (mole s)^(-1) */
  sum=-kb*X[specP1]; /* mole cm^(-3) (s)^(-1) */
  if (exp(-dG0/(calR*T))>Kcmin){
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6; /* kg m^(-3) s^(-1) */
  }
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fwbw_2r3p(int specR1, int specR2,
                            int specP1, int specP2, int specP3,
                            double A, double n, double E, double T, spec_t X, spec_t W){
  add_to_W_fw_2r3p(specR1, specR2, specP1, specP2, specP3, A, n, E, T, X, W);
  add_to_W_bw_2r3p(specR1, specR2, specP1, specP2, specP3, A, n, E, T, X, W);
}



/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_2r3p(int specR1, int specR2,
                      int specP1, int specP2, int specP3,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  sum=kf*X[specR1]*X[specR2];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; 
  W[specP3]+=_calM(specP3)*sum*1e6; 
  /* kg m^(-3) s^(-1) */
}


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_bw_2r3p(int specR1, int specR2,
                      int specP1, int specP2, int specP3,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     +_Gs(specP3,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  if (exp(-dG0/(calR*T))>Kcmin){
    Kc *= 1E-6 * THERMO_P_REF/(calR*T); /* mole/cm^3 */
    kb=kf/Kc; /*  cm^6 (mole)^(-2) s^(-1) */
    sum=-kb*X[specP1]*X[specP2]*X[specP3]; /* mole cm^(-3) (s)^(-1) */
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6;
    W[specP2]+=_calM(specP2)*sum*1e6; 
    W[specP3]+=_calM(specP3)*sum*1e6; 
    /* kg m^(-3) s^(-1) */
  }
}



/* A in cm^3 (mole s)^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fwbw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t W)
{
 
  add_to_W_fw_1r2p(specR1, specP1, specP2, A, n, E, T, X, W);
  add_to_W_bw_1r2p(specR1, specP1, specP2, A, n, E, T, X, W);
 
}

/* A in cm^3 (mole s)^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t W)
{
 
  double kf,sum;
 
  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0)); /* cm^3 (mol s)^(-1) */
  sum=kf*X[specR1];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6;        /* kg m^(-3) s^(-1) */     
 
}

/* A in cm^3 (mole s)^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_bw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t W)
{
 
  double kf,dG0,Kc,kb,sum;
 
  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T);                      /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  if (exp(-dG0/(calR*T))>Kcmin)
  {
   Kc *= 1e-6 * THERMO_P_REF/(calR*T);    /* mole/cm^3 */
    kb=kf/Kc;                             /* cm^3 (mole s)^(-1) */
    sum=-kb*X[specP1]*X[specP2];
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6;
    W[specP2]+=_calM(specP2)*sum*1e6;     /* kg m^(-3) s^(-1) */
   
  }
 
}


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fwbw_2r2p(int specR1, int specR2,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, 
                             spec_t dWdT, spec2_t dWdrhok){
  add_to_dW_fw_2r2p(specR1, specR2, specP1, specP2, A, n, E, T, X, dWdT, dWdrhok);
  add_to_dW_bw_2r2p(specR1, specR2, specP1, specP2, A, n, E, T, X, dWdT, dWdrhok);
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_2r2p(int specR1, int specR2,
                       int specP1, int specP2,
                       double A, double n, double E, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,sum,dkfdT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  sum=X[specR1]*X[specR2]; /* mole cm^(-3) (s)^(-1) */
/*  W[specR1]-=_calM(specR1)*kf*1e6*sum;
  W[specR2]-=_calM(specR2)*kf*1e6*sum;
  W[specP1]+=_calM(specP1)*kf*1e6*sum;
  W[specP2]+=_calM(specP2)*kf*1e6*sum;*/ /* kg m^(-3) s^(-1) */
  /* dWdT */
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum ; 
 
  /* dW[specR1]dX */
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR1][specR2]-=_calM(specR1)/_calM(specR2)*kf*X[specR1];

  /* dW[specR2]dX */
  dWdrhok[specR2][specR1]-=_calM(specR2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR2][specR2]-=_calM(specR2)/_calM(specR2)*kf*X[specR1];

  /* dW[specP1]dX */
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP1][specR2]+=_calM(specP1)/_calM(specR2)*kf*X[specR1];

  /* dW[specP2]dX */
  dWdrhok[specP2][specR1]+=_calM(specP2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP2][specR2]+=_calM(specP2)/_calM(specR2)*kf*X[specR1];

} 



/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_bw_2r2p(int specR1, int specR2,
                       int specP1, int specP2,
                       double A, double n, double E, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  sum=-X[specP1]*X[specP2]/Kc; /* mole cm^(-3) (s)^(-1) */
  /* dWdT */
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T);
  dsumdT=X[specP1]*X[specP2]/Kc*dKcdToverKc;
  if (exp(-dG0/(calR*T))<Kcmin){
    dKcdToverKc=0.0;
    dsumdT=0.0;
  }
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  if (exp(-dG0/(calR*T))>Kcmin){
    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
    dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum + _calM(specP2)*kf*1e6*dsumdT; 
 
    /* dW[specR1]dX */
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR1][specP2]-=-_calM(specR1)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specR2]dX */
    dWdrhok[specR2][specP1]-=-_calM(specR2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR2][specP2]-=-_calM(specR2)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specP1]dX */
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP1][specP2]+=-_calM(specP1)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specP2]dX */
    dWdrhok[specP2][specP1]+=-_calM(specP2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP2][specP2]+=-_calM(specP2)/_calM(specP2)*kf*X[specP1]/Kc;
  }
} 



/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fwbw_3r2p(int specR1, int specR2, int specR3,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, 
                             spec_t dWdT, spec2_t dWdrhok){
  add_to_dW_fw_3r2p(specR1, specR2, specR3, specP1, specP2, A, n, E, T, X, dWdT, dWdrhok);
  add_to_dW_bw_3r2p(specR1, specR2, specR3, specP1, specP2, A, n, E, T, X, dWdT, dWdrhok);
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_3r2p(int specR1, int specR2, int specR3,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, 
                             spec_t dWdT, spec2_t dWdrhok){
  double kf,sum,dkfdT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */

  sum=X[specR1]*X[specR2]*X[specR3]; /* mole cm^(-3) (s)^(-1) */

  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specR3] -= _calM(specR3)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum ; 

  /* dW[specR1]dX */
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specR1][specR2]-=_calM(specR1)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specR1][specR3]-=_calM(specR1)/_calM(specR3)*kf*X[specR1]*X[specR2];

  /* dW[specR2]dX */
  dWdrhok[specR2][specR1]-=_calM(specR2)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specR2][specR2]-=_calM(specR2)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specR2][specR3]-=_calM(specR2)/_calM(specR3)*kf*X[specR1]*X[specR2];

  /* dW[specR3]dX */
  dWdrhok[specR3][specR1]-=_calM(specR3)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specR3][specR2]-=_calM(specR3)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specR3][specR3]-=_calM(specR3)/_calM(specR3)*kf*X[specR1]*X[specR2];

  /* dW[specP1]dX */
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specP1][specR2]+=_calM(specP1)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specP1][specR3]+=_calM(specP1)/_calM(specR3)*kf*X[specR1]*X[specR2];

  /* dW[specP2]dX */
  dWdrhok[specP2][specR1]+=_calM(specP2)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specP2][specR2]+=_calM(specP2)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specP2][specR3]+=_calM(specP2)/_calM(specR3)*kf*X[specR1]*X[specR2];
  
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_bw_3r2p(int specR1, int specR2, int specR3,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, 
                             spec_t dWdT, spec2_t dWdrhok){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T)
     -_Gs(specR3,T); /* J/mole */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* mole/cm^3 */

  sum=-X[specP1]*X[specP2]/Kc; /* mole cm^(-3) (s)^(-1) */

  /* dWdT */
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T)-_dGsdT(specR3,T);
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T)+1.0/T;

  dsumdT=X[specP1]*X[specP2]/Kc*dKcdToverKc;
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  if (exp(-dG0/(calR*T))>Kcmin){

    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specR3] -= _calM(specR3)*dkfdT*1e6*sum + _calM(specR3)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
    dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum + _calM(specP2)*kf*1e6*dsumdT; 

    /* dW[specR1]dX */
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR1][specP2]-=-_calM(specR1)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specR2]dX */
    dWdrhok[specR2][specP1]-=-_calM(specR2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR2][specP2]-=-_calM(specR2)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specR3]dX */
    dWdrhok[specR3][specP1]-=-_calM(specR3)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR3][specP2]-=-_calM(specR3)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specP1]dX */
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP1][specP2]+=-_calM(specP1)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specP2]dX */
    dWdrhok[specP2][specP1]+=-_calM(specP2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP2][specP2]+=-_calM(specP2)/_calM(specP2)*kf*X[specP1]/Kc;
  }
} 



/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fwbw_2r1p(int specR1, int specR2,
                             int specP1, 
                             double A, double n, double E, double T, spec_t X, 
                             spec_t dWdT, spec2_t dWdrhok){
  add_to_dW_fw_2r1p(specR1, specR2, specP1, A, n, E, T, X, dWdT, dWdrhok);
  add_to_dW_bw_2r1p(specR1, specR2, specP1, A, n, E, T, X, dWdT, dWdrhok);
} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_2r1p(int specR1, int specR2,
                       int specP1, 
                       double A, double n, double E, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,sum,dkfdT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */

  sum=X[specR1]*X[specR2]; /* mole cm^(-3) (s)^(-1) */

  /* dWdT */
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  

  /* dW[specR1]dX */
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR1][specR2]-=_calM(specR1)/_calM(specR2)*kf*X[specR1];

  /* dW[specR2]dX */
  dWdrhok[specR2][specR1]-=_calM(specR2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR2][specR2]-=_calM(specR2)/_calM(specR2)*kf*X[specR1];

  /* dW[specP1]dX */
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP1][specR2]+=_calM(specP1)/_calM(specR2)*kf*X[specR1];

} 


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_bw_2r1p(int specR1, int specR2,
                       int specP1, 
                       double A, double n, double E, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* mole/cm^3 */

  sum=-X[specP1]/Kc; /* mole cm^(-3) (s)^(-1) */

  /* dWdT */
  dG0dT=_dGsdT(specP1,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T)+1.0/T;

  dsumdT=X[specP1]/Kc*dKcdToverKc;
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  if (exp(-dG0/(calR*T))>Kcmin){
    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
  

    /* dW[specR1]dX */
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf/Kc;

    /* dW[specR2]dX */
    dWdrhok[specR2][specP1]-=-_calM(specR2)/_calM(specP1)*kf/Kc;

    /* dW[specP1]dX */
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf/Kc;
  }

} 




/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fwbw_2r3p(int specR1, int specR2,
                             int specP1, int specP2, int specP3,
                             double A, double n, double E, double T, spec_t X, 
                             spec_t dWdT, spec2_t dWdrhok){
  add_to_dW_fw_2r3p(specR1, specR2, specP1, specP2, specP3, A, n, E, T, X, dWdT, dWdrhok);
  add_to_dW_bw_2r3p(specR1, specR2, specP1, specP2, specP3, A, n, E, T, X, dWdT, dWdrhok);
}



/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_2r3p(int specR1, int specR2,
                       int specP1, int specP2, int specP3,
                       double A, double n, double E, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,sum,dkfdT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  sum=X[specR1]*X[specR2]; /* mole cm^(-3) (s)^(-1) */

  // dWdT 
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum ; 
  dWdT[specP3] += _calM(specP3)*dkfdT*1e6*sum ; 

  // dW[specR1]dX 
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR1][specR2]-=_calM(specR1)/_calM(specR2)*kf*X[specR1];

  // dW[specR2]dX 
  dWdrhok[specR2][specR1]-=_calM(specR2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR2][specR2]-=_calM(specR2)/_calM(specR2)*kf*X[specR1];

  // dW[specP1]dX 
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP1][specR2]+=_calM(specP1)/_calM(specR2)*kf*X[specR1];

  // dW[specP2]dX 
  dWdrhok[specP2][specR1]+=_calM(specP2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP2][specR2]+=_calM(specP2)/_calM(specR2)*kf*X[specR1];

  // dW[specP3]dX 
  dWdrhok[specP3][specR1]+=_calM(specP3)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP3][specR2]+=_calM(specP3)/_calM(specR2)*kf*X[specR1];
  
}


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_bw_2r3p(int specR1, int specR2,
                       int specP1, int specP2, int specP3,
                       double A, double n, double E, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     +_Gs(specP3,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=1E-6 * THERMO_P_REF/(calR*T) * max(Kcmin, exp(-dG0/(calR*T))); /* mole/cm^3 */

  sum=-X[specP1]*X[specP2]*X[specP3]/Kc; /* mole cm^(-3) (s)^(-1) */

  // dWdT 
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)+_dGsdT(specP3,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  dKcdToverKc=( dG0/(calR*sqr(T)) - dG0dT/(calR*T) - 1.0/T);
  if (exp(-dG0/(calR*T))<Kcmin) dKcdToverKc=-1.0/T;


  dsumdT=X[specP1]*X[specP2]*X[specP3]/Kc*dKcdToverKc;

  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  if (exp(-dG0/(calR*T))>Kcmin ){

    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
    dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum + _calM(specP2)*kf*1e6*dsumdT; 
    dWdT[specP3] += _calM(specP3)*dkfdT*1e6*sum + _calM(specP3)*kf*1e6*dsumdT; 

    // dW[specR1]dX 
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specR1][specP2]-=-_calM(specR1)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specR1][specP3]-=-_calM(specR1)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc;

    // dW[specR2]dX 
    dWdrhok[specR2][specP1]-=-_calM(specR2)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specR2][specP2]-=-_calM(specR2)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specR2][specP3]-=-_calM(specR2)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc;

    // dW[specP1]dX 
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specP1][specP2]+=-_calM(specP1)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specP1][specP3]+=-_calM(specP1)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc;

    // dW[specP2]dX 
    dWdrhok[specP2][specP1]+=-_calM(specP2)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specP2][specP2]+=-_calM(specP2)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specP2][specP3]+=-_calM(specP2)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc;

    // dW[specP3]dX 
    dWdrhok[specP3][specP1]+=-_calM(specP3)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specP3][specP2]+=-_calM(specP3)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specP3][specP3]+=-_calM(specP3)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc; 
  }
}


/* A in cm^3 (mole s)^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fwbw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t dWdT, spec2_t dWdrhok)
{
 
  add_to_dW_fw_1r2p(specR1, specP1, specP2, A, n, E, T, X, dWdT, dWdrhok);
  add_to_dW_bw_1r2p(specR1, specP1, specP2, A, n, E, T, X, dWdT, dWdrhok);
 
}


/* A in cm^3 (mole s)^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t dWdT, spec2_t dWdrhok)
{
 
  double kf,sum,dkfdT;
 
  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0)); /* cm^3 (mole s)^(-1) */
  sum=X[specR1];
  dkfdT=kf*n/T + kf*(E/(sqr(T)*1.987192004e0));
 
  dWdT[specR1]-=_calM(specR1)*dkfdT*1e6*sum;
  dWdT[specP1]+=_calM(specP1)*dkfdT*1e6*sum;
  dWdT[specP2]+=_calM(specP2)*dkfdT*1e6*sum;
 
  /* dW[specR1]dX */
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf;
 
  /* dW[specP1]dX */
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf;
 
  /* dW[specP2]dX */
  dWdrhok[specP2][specR1]+=_calM(specP2)/_calM(specR1)*kf;   

}


/* A in cm^3 (mole s)^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_bw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t dWdT, spec2_t dWdrhok)
{
 
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;
 
  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T);
  Kc=1e-6 * THERMO_P_REF/(calR*T)*max(Kcmin,exp(-dG0/(calR*T)));
  sum=-X[specP1]*X[specP2]/Kc;
 
  /* dWdT */
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)-_dGsdT(specR1,T);
  dKcdToverKc=((dG0/(calR*sqr(T)))-dG0dT/(calR*T)-1.0/T);
 
    if(exp(-dG0/(calR*T))<Kcmin)
  {
    dKcdToverKc=-1.0/T;
  }
 
  dsumdT=X[specP1]*X[specP2]/Kc*dKcdToverKc;
 

 
  dkfdT=kf*n/T + kf*(E/(sqr(T)*1.987192004e0));
 
  if (exp(-dG0/(calR*T))>Kcmin)
  {
   
    dWdT[specR1]-=_calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specP1]+=_calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
    dWdT[specP2]+=_calM(specP2)*dkfdT*1e6*sum + _calM(specP2)*kf*1e6*dsumdT;

    /* dW[specR1]dX */
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR1][specP2]-=-_calM(specR1)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specP1]dX */
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP1][specP2]+=-_calM(specP1)/_calM(specP2)*kf*X[specP1]/Kc;

    /* dW[specP2]dX */
    dWdrhok[specP2][specP1]+=-_calM(specP2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP2][specP2]+=-_calM(specP2)/_calM(specP2)*kf*X[specP1]/Kc;       
       
  }      
}


void test_dW_dx(gl_t *gl, spec_t rhokref, spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam){
  long r,s; 
  spec_t dWdT,dWdTe,dWdTv,dWdQbeam;
  spec2_t dWrhok,dWrhok2;   
  spec_t rhok2,dWdT2,W,W2;
  double drhok,dQbeam,dT,N,N2,Estar2;
 
  find_W(gl, rhok, T, Te, Tv, Estar, Qbeam, W);      
  find_dW_dx(gl, rhok, mu, T, Te, Tv, Estar, Qbeam,
              dWrhok, dWdT, dWdTe, dWdTv, dWdQbeam);

  N=0.0;
  for (s=0; s<ns; s++) N+=rhok[s]/_m(s);   

  for (s=0; s<ns; s++) wfprintf(stdout,"rho[%ld]=%E kg/m3\n",s,rhok[s]);    
  wfprintf(stdout,"T=%E K\nTe=%E K\nTv=%E K\nE/N=%E Vm2\nQbeam=%E W/m3\n",T,Te,Tv,Estar,Qbeam);

  // find dWdQbeam numerically    
   wfprintf(stdout,"\n\ndWdQbeam:\n");
   dQbeam=Qbeam/10000.0+1.0e2;
   find_W(gl, rhok, T, Te,Tv,Estar, Qbeam+dQbeam, W2);
   for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dQbeam;
   for (s=0; s<ns; s++){  
     wfprintf(stdout,"%15.15E  %15.15E\n",dWdQbeam[s],dWdT2[s]);
   } 
      
  /* find dWrhok numerically */
  for (s=0; s<ns; s++){
    for (r=0; r<ns; r++) rhok2[r]=rhok[r];
    drhok=max(0.0,rhok[s]/1e3)+1.0e-8*rhokref[s];
    rhok2[s]+=drhok;
    N2=N+(rhok2[s]-rhok[s])/_m(s);
    Estar2=Estar*N/N2;
    find_W(gl, rhok2, T, Te,Tv,Estar2,Qbeam, W2);
    for (r=0; r<ns; r++) dWrhok2[r][s]=(W2[r]-W[r])/(drhok);
  }
  for (r=0; r<ns; r++){
    wfprintf(stdout,"\n\ndW[%ld]drhok:\n",r);
    for (s=0; s<ns; s++){  
      wfprintf(stdout,"%15.15E  %15.15E\n",dWrhok[s][r],dWrhok2[s][r]);
    }
  }
  for (s=0; s<ns; s++){
    wfprintf(stdout,"\n\ndWdrho[%ld]:\n",s);
    for (r=0; r<ns; r++){  
      wfprintf(stdout,"%15.15E  %15.15E\n",dWrhok[s][r],dWrhok2[s][r]);
    }
  }

  // find dWdT numerically 
   wfprintf(stdout,"\n\ndWdT:\n");
   dT=T/10000.0;
  find_W(gl, rhok, T+dT, Te,Tv,Estar,Qbeam, W2);
  for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dT;
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdT[s],dWdT2[s]);
  } 


  // find dWdTv numerically 
   wfprintf(stdout,"\n\ndWdTv:\n");
  dT=Tv/10000.0;
  find_W(gl, rhok, T, Te,Tv+dT,Estar,Qbeam, W2);
  for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dT;
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdTv[s],dWdT2[s]);
  }

  
  // find dWdTe numerically 
   wfprintf(stdout,"\n\ndWdTe:\n");
  dT=Te/10000.0;
  find_W(gl, rhok, T, Te+dT,Tv,Estar,Qbeam, W2);
  for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dT;
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdTe[s],dWdT2[s]);
  }
    
  exit(1);  
}



/* find the forward reaction rate cofficient kf = A * (calA)^(1-numreactant) * T^n * exp(-E/(R*T))
 * 
 * numreactant : the number of reactants
 * A           : pre-exponential factor in cm^3 * (mole s)^(-1) K^(-n)
 * n           : temperature exponent
 * E           : activation energy of the reaction in cal/mole
 * R           : universal gas constant = 1.9872	cal/(K mol)
 * kf          : reaction rate coefficient in cm^3/s (2 reactants) or cm^6/s (3 reactants)
 * T           : temperature in K
 */
double _kf_Arrhenius(long numreactant, double A, double n, double E, double T){
  double kf,R;
  R=1.9872;
  switch (numreactant){
    case 2:
      kf=A/calA*pow(T,n)*exp(-E/(R*T));
    break;
    case 3:
      kf=A/sqr(calA)*pow(T,n)*exp(-E/(R*T));    
    break;
    default:
      fatal_error("numreactant can not be set to %ld in _kf_Arrhenius",numreactant);
      kf=0.0;
  }
  return(kf);
}


/* find derivative of kf with respect to T with kf=A * (calA)^(1-numreactant) * T^n * exp(-E/(R*T))
 * 
 * numreactant : the number of reactants
 * A           : pre-exponential factor in cm^3 * (mole s)^(-1) K^(-n)
 * n           : temperature exponent
 * E           : activation energy of the reaction in cal/mole
 * R           : universal gas constant = 1.9872	cal/(K mol)
 * dkfdT       : cm^3/sK (2 reactants) or cm^6/sK (3 reactants)
 * T           : temperature in K
 */
double _dkfdT_Arrhenius(long numreactant, double A, double n, double E, double T){
  double dkfdT,R;
  R=1.9872;
  switch (numreactant){
    case 2:
      dkfdT=(A/calA)*n*pow(T,n-1.0)*exp(-E/(R*T))+E/(R*T*T)*(A/calA)*pow(T,n)*exp(-E/(R*T));
    break;
    case 3:
      dkfdT=(A/sqr(calA))*n*pow(T,n-1.0)*exp(-E/(R*T))+E/(R*T*T)*(A/sqr(calA))*pow(T,n)*exp(-E/(R*T));
    break;
    default:
      fatal_error("numreactant can not be set to %ld in _kf_Arrhenius",numreactant);
      dkfdT=0.0;
  }
  return(dkfdT);
}


/* add contributions to Qei coming from one electron impact chemical reaction
 * spec is the neutral molecule the electrons impact
 * kf is the reaction rate in cm3/s
 * rhok is the partial densities in kg/m3
 * Qei is the heat removed from the electrons in J/m3
 */
void add_to_Qei(long spec, double kf, spec_t rhok, double *Qei){
#ifdef TEST
  kf=1e-17;
#endif
  (*Qei) += kf * 1E-6* _ionizationpotential(spec) * rhok[speceminus] / _calM ( speceminus ) * rhok[spec] / _calM ( spec ) * sqr(calA);
}


/* add contributions to dQei_dx coming from one electron impact chemical reaction
 * spec is the neutral molecule the electrons impact
 * kf is the reaction rate in cm3/s
 * dkfdTe is the derivative of kf with respect to Te in cm3/(s K)
 * rhok is the partial densities in kg/m3
 * dQeidrhok is in J/kg 
 * dQeidTe is in J/(m3 K)
 */
void add_to_dQei(long spec, double kf, double dkfdTe, spec_t rhok, spec_t dQeidrhok, double *dQeidTe){
#ifdef TEST
  kf=1e-17;
#endif
  dQeidrhok[spec] += kf *1e-6* _ionizationpotential(spec) * rhok[speceminus] / _calM ( speceminus )  / _calM ( spec ) * sqr(calA);
  dQeidrhok[speceminus] += kf *1e-6* _ionizationpotential(spec) / _calM ( speceminus ) * rhok[spec] / _calM ( spec ) * sqr(calA);
}
