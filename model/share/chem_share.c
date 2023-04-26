// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2018-2020 Bernard Parent
Copyright 2020 Ajjay Omprakas
Copyright 2021-2022 Prasanna Thoguluva Rajendran

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

#define TREF_EXCI 298.0

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
  
  dWdTe[specR1]-=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR1) * 1.0e6;
  dWdTe[specR2]-=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR2) * 1.0e6;
  dWdTe[specR3]-=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR3) * 1.0e6;
  dWdTe[specP1]+=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP1) * 1.0e6;
  dWdTe[specP2]+=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP2) * 1.0e6;
  dWdTe[specP3]+=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP3) * 1.0e6;
  
  dWdTv[specR1]-=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR1) * 1.0e6;
  dWdTv[specR2]-=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR2) * 1.0e6;
  dWdTv[specR3]-=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR3) * 1.0e6;
  dWdTv[specP1]+=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP1) * 1.0e6;
  dWdTv[specP2]+=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP2) * 1.0e6;
  dWdTv[specP3]+=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP3) * 1.0e6;
  
  dWdT[specR1]-=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR1) * 1.0e6;
  dWdT[specR2]-=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR2) * 1.0e6;
  dWdT[specR3]-=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR3) * 1.0e6;
  dWdT[specP1]+=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP1) * 1.0e6;
  dWdT[specP2]+=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP2) * 1.0e6;
  dWdT[specP3]+=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP3) * 1.0e6;
  
  
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
  
  dWdTe[specR1]-=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR1) * 1.0e6;
  dWdTe[specR2]-=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR2) * 1.0e6;
  dWdTe[specR3]-=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR3) * 1.0e6;
  dWdTe[specP1]+=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP1) * 1.0e6;
  dWdTe[specP2]+=dkfdTe*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP2) * 1.0e6;
  
  dWdTv[specR1]-=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR1) * 1.0e6;
  dWdTv[specR2]-=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR2) * 1.0e6;
  dWdTv[specR3]-=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR3) * 1.0e6;
  dWdTv[specP1]+=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP1) * 1.0e6;
  dWdTv[specP2]+=dkfdTv*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP2) * 1.0e6;
  
  dWdT[specR1]-=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR1) * 1.0e6;
  dWdT[specR2]-=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR2) * 1.0e6;
  dWdT[specR3]-=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specR3) * 1.0e6;
  dWdT[specP1]+=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP1) * 1.0e6;
  dWdT[specP2]+=dkfdT*N[specR1]*N[specR2]*N[specR3]/ calA * _calM(specP2) * 1.0e6;
  
  
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
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



/* A in cm^6 mole^(-2) s^(-1) K^(-n)
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



/* A in cm^6 mole^(-2) s^(-1) K^(-n)
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_3r2p(int specR1, int specR2, int specR3,
                      int specP1, int specP2,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^6 mole^(-2) s^(-1) */
  sum=kf*X[specR1]*X[specR2]*X[specR3]; /* mole cm^(-3) (s)^(-1) */
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specR3]-=_calM(specR3)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
} 


/* A in cm^6 mole^(-2) s^(-1) K^(-n)
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_bw_3r2p(int specR1, int specR2, int specR3,
                      int specP1, int specP2,
                      double A, double n, double E, double T, spec_t X, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^6 mole^(-2) s^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T)
     -_Gs(specR3,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* cm^3/mole */
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* cm^3/mole */
  kb=kf/Kc; /*  s^(-1) */
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
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



/* A in s^(-1) K^(-n)
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

/* A in s^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t W)
{
 
  double kf,sum;
 
  kf=A*pow(T,n)*exp(-E/(T*Rchem)); /* s^(-1) */
  sum=kf*X[specR1];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6;        /* kg m^(-3) s^(-1) */     
 
}

/* A in s^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_bw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t W)
{
 
  double kf,dG0,Kc,kb,sum;
 
  kf=A*pow(T,n)*exp(-E/(T*Rchem));/* s^(-1) */
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



/* A in s^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   W in kg m^(-3) s^(-1)
*/
void add_to_W_fw_1r1p(int specR1, int specP1, double A, double n, double E, double T, spec_t X, spec_t W)
{
 
  double kf,sum;
 
  kf=A*pow(T,n)*exp(-E/(T*Rchem)); /* s^(-1) */
  sum=kf*X[specR1];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6; /* kg m^(-3) s^(-1) */     
 
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
  sum=X[specR1]*X[specR2];
/*  W[specR1]-=_calM(specR1)*kf*1e6*sum;
  W[specR2]-=_calM(specR2)*kf*1e6*sum;
  W[specP1]+=_calM(specP1)*kf*1e6*sum;
  W[specP2]+=_calM(specP2)*kf*1e6*sum;*/ /* kg m^(-3) s^(-1) */
  /* dWdT */
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  sum=-X[specP1]*X[specP2]/Kc; /* mole^2 cm^(-6) */
  /* dWdT */
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T);
  dsumdT=X[specP1]*X[specP2]/Kc*dKcdToverKc;
  if (exp(-dG0/(calR*T))<Kcmin){
    dKcdToverKc=0.0;
    dsumdT=0.0;
  }
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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



/* A in cm^6 mole^(-2) s^(-1) K^(-n)
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


/* A in cm^6 mole^(-2) s^(-1) K^(-n)
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */

  sum=X[specR1]*X[specR2]*X[specR3]; /* mole^3 cm^(-9) */

  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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


/* A in cm^6 mole^(-2) s^(-1) K^(-n)
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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^6 mole^(-2) s^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T)
     -_Gs(specR3,T); /* J/mole */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* cm^3/mole */

  sum=-X[specP1]*X[specP2]/Kc; /* mole^3 cm^(-9) */

  /* dWdT */
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T)-_dGsdT(specR3,T);
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T)+1.0/T;

  dsumdT=X[specP1]*X[specP2]/Kc*dKcdToverKc;
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */

  sum=X[specR1]*X[specR2];

  /* dWdT */
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* cm^3/mole */

  sum=-X[specP1]/Kc; /* mole^2 cm^(-6) */

  /* dWdT */
  dG0dT=_dGsdT(specP1,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T)+1.0/T;

  dsumdT=X[specP1]/Kc*dKcdToverKc;
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
  sum=X[specR1]*X[specR2];

  // dWdT 
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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

  kf=A*pow(T,n)*exp(-E/(T*Rchem));  /* cm^3 (mole s)^(-1) */
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     +_Gs(specP3,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=1E-6 * THERMO_P_REF/(calR*T) * max(Kcmin, exp(-dG0/(calR*T))); /* mole/cm^3 */

  sum=-X[specP1]*X[specP2]*X[specP3]/Kc; /* mole^2 cm^(-6)*/

  // dWdT 
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)+_dGsdT(specP3,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  dKcdToverKc=( dG0/(calR*sqr(T)) - dG0dT/(calR*T) - 1.0/T);
  if (exp(-dG0/(calR*T))<Kcmin) dKcdToverKc=-1.0/T;


  dsumdT=X[specP1]*X[specP2]*X[specP3]/Kc*dKcdToverKc;

  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*Rchem));

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



/* A in s^(-1) K^(-n)
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


/* A in s^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t dWdT, spec2_t dWdrhok)
{
 
  double kf,sum,dkfdT;
 
  kf=A*pow(T,n)*exp(-E/(T*Rchem)); /* s^(-1) */
  sum=X[specR1];
  dkfdT=kf*n/T + kf*(E/(sqr(T)*Rchem));
 
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


/* A in s^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_bw_1r2p(int specR1, int specP1, int specP2, double A, double n, double E, double T, spec_t X, spec_t dWdT, spec2_t dWdrhok)
{
 
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;
 
  kf=A*pow(T,n)*exp(-E/(T*Rchem));
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
 

 
  dkfdT=kf*n/T + kf*(E/(sqr(T)*Rchem));
 
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

/* A in s^(-1) K^(-n)
   E in cal mole^(-1)
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_1r1p(int specR1, int specP1, double A, double n, double E, double T, spec_t X, spec_t dWdT, spec2_t dWdrhok)
{
 
  double kf,sum,dkfdT;
 
  kf=A*pow(T,n)*exp(-E/(T*Rchem)); /* s^(-1) */
  sum=X[specR1];
  dkfdT=kf*n/T + kf*(E/(sqr(T)*Rchem));
 
  dWdT[specR1]-=_calM(specR1)*dkfdT*1e6*sum;
  dWdT[specP1]+=_calM(specP1)*dkfdT*1e6*sum;
 
  /* dW[specR1]dX */
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf;
 
  /* dW[specP1]dX */
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf;
 

}


/* 
   k = A*T^n*exp(E1/(R*T)+E2/(R*T)^2+E3/(R*T)^3+E4/(R*T)^4)
   for A in cm^3 (mole s)^(-1) K^(-n), kf will be in cm^3 (mole s)^(-1)
   for A in cm^3 (s)^(-1) K^(-n), kf will be in cm^3 (s)^(-1)     
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
*/
double _kf_fit4(double A, double n, double E1, double E2, double E3, double E4, double T){
  double kf;
  kf=A*pow(T,n)*exp(E1/(Rchem*T)+E2/(Rchem*Rchem*T*T)+E3/(Rchem*Rchem*Rchem*T*T*T)
    +E4/(Rchem*Rchem*Rchem*Rchem*T*T*T*T)); 
  return(kf);
}


/* 
   kf = A*T^n*exp(E1/(R*T)+E2/(R*T)^2+E3/(R*T)^3+E4/(R*T)^4)
   for A in cm^3 (mole s)^(-1) K^(-n), dkfdT will be in  cm^3 (K mole s)^(-1)
   for A in cm^3 (s)^(-1) K^(-n), dkfdT will be in  cm^3 (K s)^(-1)
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
*/
double _dkfdT_fit4(double A, double n, double E1, double E2, double E3, double E4, double T){
  double dkfdT,kf;
  kf=_kf_fit4(A,n,E1,E2,E3,E4,T);
  dkfdT= kf*n/T 
       + kf*(-E1/(Rchem*T*T)-2.0*E2/(Rchem*Rchem*T*T*T)-3.0*E3/(Rchem*Rchem*Rchem*T*T*T*T)
       - 4.0*E4/(Rchem*Rchem*Rchem*Rchem*T*T*T*T*T));
  return(dkfdT);
}


/* 
   k = A*T^n*exp(E1/(R*T)+E2/(R*T)^2+E3/(R*T)^3+E4/(R*T)^4)
   A in cm^3 (mole s)^(-1) K^(-n) 
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_W_fw_2r2p_fit4(int specR1, int specR2,
                           int specP1, int specP2,
                           double A, double n, double E1, double E2, double E3, double E4, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=_kf_fit4(A,n,E1,E2,E3,E4,T);  /* cm^3 (mole s)^(-1) */
  sum=kf*X[specR1]*X[specR2];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
                             
}


/* 
   k = A*T^n*exp(E1/(R*T)+E2/(R*T)^2+E3/(R*T)^3+E4/(R*T)^4)
   A in cm^3 (mole s)^(-1) K^(-n) 
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_W_fw_2r3p_fit4(int specR1, int specR2,
                      int specP1, int specP2, int specP3,
                      double A, double n, double E1, double E2, double E3, double E4, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=_kf_fit4(A,n,E1,E2,E3,E4,T);  /* cm^3 (mole s)^(-1) */
  sum=kf*X[specR1]*X[specR2];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; 
  W[specP3]+=_calM(specP3)*sum*1e6; 
  /* kg m^(-3) s^(-1) */
}



/* 
   k = A*T^n*exp(E1/(R*T)+E2/(R*T)^2+E3/(R*T)^3+E4/(R*T)^4)
   A in cm^3 (mole s)^(-1) K^(-n) 
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_W_fw_2r4p_fit4(int specR1, int specR2,
                      int specP1, int specP2, int specP3, int specP4,
                      double A, double n, double E1, double E2, double E3, double E4, double T, spec_t X, spec_t W){
  double kf,sum;

  kf=_kf_fit4(A,n,E1,E2,E3,E4,T);  /* cm^3 (mole s)^(-1) */
  sum=kf*X[specR1]*X[specR2];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; 
  W[specP3]+=_calM(specP3)*sum*1e6; 
  W[specP4]+=_calM(specP4)*sum*1e6; 
  /* kg m^(-3) s^(-1) */
}


/* A in cm^3 (mole s)^(-1) K^(-n) 
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_2r2p_fit4(int specR1, int specR2,
                            int specP1, int specP2,
                            double A, double n, double E1, double E2, double E3, double E4,
                            double T, spec_t X, 
                            spec_t dWdT, spec2_t dWdrhok){
  double kf,sum,dkfdT;

  kf=_kf_fit4(A,n,E1,E2,E3,E4,T);  /* cm^3 (mole s)^(-1) */
  sum=X[specR1]*X[specR2]; /* mole^2 cm^(-6) */
/*  W[specR1]-=_calM(specR1)*kf*1e6*sum;
  W[specR2]-=_calM(specR2)*kf*1e6*sum;
  W[specP1]+=_calM(specP1)*kf*1e6*sum;
  W[specP2]+=_calM(specP2)*kf*1e6*sum;*/ /* kg m^(-3) s^(-1) */
  /* dWdT */
  dkfdT=_dkfdT_fit4(A, n, E1, E2, E3, E4, T);

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
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_2r3p_fit4(int specR1, int specR2,
                       int specP1, int specP2, int specP3,
                       double A, double n, double E1, double E2, double E3, double E4, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,sum,dkfdT;

  kf=_kf_fit4(A,n,E1,E2,E3,E4,T);  /* cm^3 (mole s)^(-1) */

  sum=X[specR1]*X[specR2]; /* mole^2 cm^(-6) */

  // dWdT 
  dkfdT=_dkfdT_fit4(A, n, E1, E2, E3, E4, T);

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
   E1 in cal/mole
   E2 in (cal/mole)^2
   E3 in (cal/mole)^3
   E4 in (cal/mole)^4
   T in Kelvin
   X in mole/cm3
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/
void add_to_dW_fw_2r4p_fit4(int specR1, int specR2,
                       int specP1, int specP2, int specP3, int specP4,
                       double A, double n, double E1, double E2, double E3, double E4, double T, spec_t X, 
                       spec_t dWdT, spec2_t dWdrhok){
  double kf,sum,dkfdT;

  kf=_kf_fit4(A,n,E1,E2,E3,E4,T);  /* cm^3 (mole s)^(-1) */

  sum=X[specR1]*X[specR2]; /* mole^2 cm^(-6) */

  // dWdT 
  dkfdT=_dkfdT_fit4(A, n, E1, E2, E3, E4, T);

  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum ; 
  dWdT[specP3] += _calM(specP3)*dkfdT*1e6*sum ; 
  dWdT[specP4] += _calM(specP4)*dkfdT*1e6*sum ; 

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

  // dW[specP4]dX 
  dWdrhok[specP4][specR1]+=_calM(specP4)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP4][specR2]+=_calM(specP4)/_calM(specR2)*kf*X[specR1];
  
}


/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0) 
 Ainf in s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_fw_2r2p_Lindemann(int specR1, int specR2,
                                int specP1, int specP2,
                                double Ainf, double ninf, double Einf,								
                                double A0, double n0, double E0,
                                double T, spec_t X, spec_t W){
                                  
  double kf,sum, k0, kinf;
  
  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  } else {
    kf = 0.0; 
  }
  
  sum=kf*X[specR1]*X[specR2];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */


}




/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0) 
 Ainf in s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_fw_2r3p_Lindemann(int specR1, int specR2,
                                int specP1, int specP2, int specP3,
                                double Ainf, double ninf, double Einf,								
                                double A0, double n0, double E0,
                                double T, spec_t X, spec_t W){
  double kf,sum, k0, kinf;
  
  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  } else {
    kf = 0.0; 
  }

  sum=kf*X[specR1]*X[specR2];
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; 
  W[specP3]+=_calM(specP3)*sum*1e6; 
}

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0) 
 Ainf in s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 

void add_to_W_bw_2r2p_Lindemann(int specR1, int specR2,
                                int specP1, int specP2,
                                double Ainf, double ninf, double Einf,								
                                double A0, double n0, double E0,
                                double T, spec_t X, spec_t W){
  double k0, kinf;
  double kf,dG0,Kc,kb,sum;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  } else {
    kf = 0.0; 
  }
  
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

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0) 
 Ainf in s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_bw_2r3p_Lindemann(int specR1, int specR2,
                                int specP1, int specP2, int specP3,
                                double Ainf, double ninf, double Einf,								
                                double A0, double n0, double E0,
                                double T, spec_t X, spec_t W){
  double k0, kinf;
  double kf,dG0,Kc,kb,sum;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  } else {
    kf = 0.0; 
  }

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

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0) 
 Ainf in s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_fwbw_2r2p_Lindemann(int specR1, int specR2,
                                  int specP1, int specP2, 
                                  double Ainf, double ninf, double Einf,								
                                  double A0, double n0, double E0,
                                  double T, spec_t X, spec_t W){
  add_to_W_fw_2r2p_Lindemann(specR1, specR2, specP1, specP2,
                             Ainf, ninf, Einf, A0, n0, E0, T, X, W);
  add_to_W_bw_2r2p_Lindemann(specR1, specR2, specP1, specP2, 
                             Ainf, ninf, Einf, A0, n0, E0, T, X, W);
}

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0) 
 Ainf in s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_fwbw_2r3p_Lindemann(int specR1, int specR2,
                                  int specP1, int specP2, int specP3,
                                  double Ainf, double ninf, double Einf,								
                                  double A0, double n0, double E0,
                                  double T, spec_t X, spec_t W){
  add_to_W_fw_2r3p_Lindemann(specR1, specR2, specP1, specP2, specP3,
                             Ainf, ninf, Einf, A0, n0, E0, T, X, W);
  add_to_W_bw_2r3p_Lindemann(specR1, specR2, specP1, specP2, specP3,
                             Ainf, ninf, Einf, A0, n0, E0, T, X, W);
}

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0)
 Ainf in  s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_fw_2r2p_Lindemann(int specR1, int specR2,
                                 int specP1, int specP2,
                                 double Ainf, double ninf, double Einf, 								 
                                 double A0, double n0, double E0,
                                 double T, spec_t X, 
                                 spec_t dWdT, spec2_t dWdrhok){
									 
  double kf,sum, k0, kinf, dkfdT, dk0dT, dkinfdT, dkfdX2;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  
  // dWdT 
    
    dk0dT  = k0*n0/T + k0*(E0/(sqr(T)*Rchem));
    dkinfdT = kinf*ninf/T + kinf*(Einf/(sqr(T)*Rchem));
  
    dkfdX2 = -(k0 * k0 * kinf) / sqr(kinf + k0*X[specR2]); /*  cm^6 (mole^2 s)^(-1)   */ 
 
    dkfdT = kinf*dk0dT / ( kinf + k0*X[specR2] ) - k0 * ( X[specR2]*kinf*dk0dT - X[specR2]*k0*dkinfdT ) / sqr(kinf + k0*X[specR2]);

  } else {
    kf = 0.0; 
    dkfdT = 0.0; 
    dkfdX2 = 0.0;
  }
  
  sum=X[specR1]*X[specR2]; 

  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum ; 

  // dW[specR1]dX 
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR1][specR2]-=_calM(specR1)/_calM(specR2)*kf*X[specR1] + _calM(specR1)/_calM(specR2) * sum * dkfdX2;

  // dW[specR2]dX 
  dWdrhok[specR2][specR1]-=_calM(specR2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR2][specR2]-=_calM(specR2)/_calM(specR2)*kf*X[specR1] + _calM(specR2)/_calM(specR2) * sum * dkfdX2;

  // dW[specP1]dX 
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP1][specR2]+=_calM(specP1)/_calM(specR2)*kf*X[specR1] + _calM(specP1)/_calM(specR2) * sum * dkfdX2;

  // dW[specP2]dX 
  dWdrhok[specP2][specR1]+=_calM(specP2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP2][specR2]+=_calM(specP2)/_calM(specR2)*kf*X[specR1] + _calM(specP2)/_calM(specR2) * sum * dkfdX2;
  
}


/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0)
 Ainf in  s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_fw_2r3p_Lindemann(int specR1, int specR2,
                                 int specP1, int specP2, int specP3,
                                 double Ainf, double ninf, double Einf, 								 
                                 double A0, double n0, double E0,
                                 double T, spec_t X, 
                                 spec_t dWdT, spec2_t dWdrhok){
									 
  double kf,sum, k0, kinf, dkfdT, dk0dT, dkinfdT, dkfdX2;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  
  // dWdT 
    
    dk0dT  = k0*n0/T + k0*(E0/(sqr(T)*Rchem));
    dkinfdT = kinf*ninf/T + kinf*(Einf/(sqr(T)*Rchem));
  
    dkfdX2 = -(k0 * k0 * kinf) / sqr(kinf + k0*X[specR2]); /*  cm^6 (mole^2 s)^(-1)   */ 
 
    dkfdT = kinf*dk0dT / ( kinf + k0*X[specR2] ) - k0 * ( X[specR2]*kinf*dk0dT - X[specR2]*k0*dkinfdT ) / sqr(kinf + k0*X[specR2]);

  } else {
    kf = 0.0; 
    dkfdT = 0.0; 
    dkfdX2 = 0.0;
  }
  
  sum=X[specR1]*X[specR2]; 

  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum ; 
  dWdT[specP3] += _calM(specP3)*dkfdT*1e6*sum ; 

  // dW[specR1]dX 
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR1][specR2]-=_calM(specR1)/_calM(specR2)*kf*X[specR1] + _calM(specR1)/_calM(specR2) * sum * dkfdX2;

  // dW[specR2]dX 
  dWdrhok[specR2][specR1]-=_calM(specR2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specR2][specR2]-=_calM(specR2)/_calM(specR2)*kf*X[specR1] + _calM(specR2)/_calM(specR2) * sum * dkfdX2;

  // dW[specP1]dX 
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP1][specR2]+=_calM(specP1)/_calM(specR2)*kf*X[specR1] + _calM(specP1)/_calM(specR2) * sum * dkfdX2;

  // dW[specP2]dX 
  dWdrhok[specP2][specR1]+=_calM(specP2)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP2][specR2]+=_calM(specP2)/_calM(specR2)*kf*X[specR1] + _calM(specP2)/_calM(specR2) * sum * dkfdX2;

  // dW[specP3]dX 
  dWdrhok[specP3][specR1]+=_calM(specP3)/_calM(specR1)*kf*X[specR2];
  dWdrhok[specP3][specR2]+=_calM(specP3)/_calM(specR2)*kf*X[specR1] + _calM(specP3)/_calM(specR2) * sum * dkfdX2;
  
}

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0)
 Ainf in  s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_bw_2r2p_Lindemann(int specR1, int specR2,
                                 int specP1, int specP2,
                                 double Ainf, double ninf, double Einf, 								 
                                 double A0, double n0, double E0,
                                 double T, spec_t X, 
                                 spec_t dWdT, spec2_t dWdrhok){
  double kinf,k0;
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT,dk0dT,dkinfdT,dkfdX2;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  
    dk0dT  = k0*n0/T + k0*(E0/(sqr(T)*Rchem));
    dkinfdT = kinf*ninf/T + kinf*(Einf/(sqr(T)*Rchem));
  
    dkfdX2 = -(k0 * k0 * kinf) / sqr(kinf + k0*X[specR2]); /*  cm^6 (mole^2 s)^(-1)   */ 
 
    dkfdT = kinf*dk0dT / ( kinf + k0*X[specR2] ) - k0 * ( X[specR2]*kinf*dk0dT - X[specR2]*k0*dkinfdT ) / sqr(kinf + k0*X[specR2]);
    
  } else {
    kf = 0.0; 
    dkfdT = 0.0; 
    dkfdX2 = 0.0;
  }
  
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  sum=-X[specP1]*X[specP2]/Kc; /* mole^2 cm^(-6) */
  /* dWdT */
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T);
  dsumdT=X[specP1]*X[specP2]/Kc*dKcdToverKc;
  if (exp(-dG0/(calR*T))<Kcmin){
    dKcdToverKc=0.0;
    dsumdT=0.0;
  }
  
  if (exp(-dG0/(calR*T))>Kcmin){
      
    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
    dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum + _calM(specP2)*kf*1e6*dsumdT; 
 
    /* dW[specR1]dX */
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR1][specP2]-=-_calM(specR1)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specR1)/_calM(specP2)*dkfdX2*X[specP1]*X[specP2]/Kc;

    /* dW[specR2]dX */
    dWdrhok[specR2][specP1]-=-_calM(specR2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR2][specP2]-=-_calM(specR2)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specR2)/_calM(specP2)*dkfdX2*X[specP1]*X[specP2]/Kc;

    /* dW[specP1]dX */
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP1][specP2]+=-_calM(specP1)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specP1)/_calM(specP2)*dkfdX2*X[specP1]*X[specP2]/Kc;

    /* dW[specP2]dX */
    dWdrhok[specP2][specP1]+=-_calM(specP2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP2][specP2]+=-_calM(specP2)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specP2)/_calM(specP2)*dkfdX2*X[specP1]*X[specP2]/Kc;

  }
  
}

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0)
 Ainf in  s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_bw_2r3p_Lindemann(int specR1, int specR2,
                                 int specP1, int specP2, int specP3,
                                 double Ainf, double ninf, double Einf, 								 
                                 double A0, double n0, double E0,
                                 double T, spec_t X, 
                                 spec_t dWdT, spec2_t dWdrhok){
  double kinf,k0;
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT,dk0dT,dkinfdT,dkfdX2;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /*  cm^3 (mole s)^(-1) */
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem)); /* s^(-1)  */ 
    kf = k0*kinf/(kinf + k0*X[specR2]);  /*  cm^3 (mole s)^(-1) */
  
    dk0dT  = k0*n0/T + k0*(E0/(sqr(T)*Rchem));
    dkinfdT = kinf*ninf/T + kinf*(Einf/(sqr(T)*Rchem));
  
    dkfdX2 = -(k0 * k0 * kinf) / sqr(kinf + k0*X[specR2]); /*  cm^6 (mole^2 s)^(-1)   */ 
 
    dkfdT = kinf*dk0dT / ( kinf + k0*X[specR2] ) - k0 * ( X[specR2]*kinf*dk0dT - X[specR2]*k0*dkinfdT ) / sqr(kinf + k0*X[specR2]);

  } else {
    kf = 0.0; 
    dkfdT = 0.0; 
    dkfdX2 = 0.0;
  }
  
  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     +_Gs(specP3,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T); /* J/mole */
  Kc=1E-6 * THERMO_P_REF/(calR*T) * max(Kcmin, exp(-dG0/(calR*T))); /* mole/cm^3 */

  sum=-X[specP1]*X[specP2]*X[specP3]/Kc; /* mole^2 cm^(-6)*/

  // dWdT 
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)+_dGsdT(specP3,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T);
  dKcdToverKc=( dG0/(calR*sqr(T)) - dG0dT/(calR*T) - 1.0/T);
  if (exp(-dG0/(calR*T))<Kcmin) dKcdToverKc=-1.0/T;

  dsumdT=X[specP1]*X[specP2]*X[specP3]/Kc*dKcdToverKc;

  if (exp(-dG0/(calR*T))>Kcmin ){

    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
    dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum + _calM(specP2)*kf*1e6*dsumdT; 
    dWdT[specP3] += _calM(specP3)*dkfdT*1e6*sum + _calM(specP3)*kf*1e6*dsumdT; 

    // dW[specR1]dX 
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specR1][specP2]-=-_calM(specR1)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specR1][specP3]-=-_calM(specR1)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc - _calM(specR1)/_calM(specP3)*dkfdX2*X[specP1]*X[specP2]*X[specP3]/Kc ;

    // dW[specR2]dX 
    dWdrhok[specR2][specP1]-=-_calM(specR2)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specR2][specP2]-=-_calM(specR2)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specR2][specP3]-=-_calM(specR2)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc - _calM(specR2)/_calM(specP3)*dkfdX2*X[specP1]*X[specP2]*X[specP3]/Kc;

    // dW[specP1]dX 
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specP1][specP2]+=-_calM(specP1)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specP1][specP3]+=-_calM(specP1)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc - _calM(specP1)/_calM(specP3)*dkfdX2*X[specP1]*X[specP2]*X[specP3]/Kc;

    // dW[specP2]dX 
    dWdrhok[specP2][specP1]+=-_calM(specP2)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specP2][specP2]+=-_calM(specP2)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specP2][specP3]+=-_calM(specP2)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc - _calM(specP2)/_calM(specP3)*dkfdX2*X[specP1]*X[specP2]*X[specP3]/Kc;

    // dW[specP3]dX 
    dWdrhok[specP3][specP1]+=-_calM(specP3)/_calM(specP1)*kf*X[specP2]*X[specP3]/Kc;
    dWdrhok[specP3][specP2]+=-_calM(specP3)/_calM(specP2)*kf*X[specP1]*X[specP3]/Kc;
    dWdrhok[specP3][specP3]+=-_calM(specP3)/_calM(specP3)*kf*X[specP1]*X[specP2]/Kc - _calM(specP3)/_calM(specP3)*dkfdX2*X[specP1]*X[specP2]*X[specP3]/Kc ;
  }
  
}

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0)
 Ainf in  s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_fwbw_2r2p_Lindemann(int specR1, int specR2,
                                   int specP1, int specP2, 
                                   double Ainf, double ninf, double Einf, 								 
                                   double A0, double n0, double E0,
                                   double T, spec_t X, 
                                   spec_t dWdT, spec2_t dWdrhok){
  add_to_dW_fw_2r2p_Lindemann(specR1, specR2, specP1, specP2,
                              Ainf, ninf, Einf, A0, n0, E0, T, X, dWdT, dWdrhok);
  add_to_dW_bw_2r2p_Lindemann(specR1, specR2, specP1, specP2,
                              Ainf, ninf, Einf, A0, n0, E0, T, X, dWdT, dWdrhok);
}

/*
 Note : kf  = k0/(1+k0*X[specR2]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^3 (mole s)^(-1) K^(-n0)
 Ainf in  s^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_fwbw_2r3p_Lindemann(int specR1, int specR2,
                                   int specP1, int specP2, int specP3,
                                   double Ainf, double ninf, double Einf, 								 
                                   double A0, double n0, double E0,
                                   double T, spec_t X, 
                                   spec_t dWdT, spec2_t dWdrhok){
  add_to_dW_fw_2r3p_Lindemann(specR1, specR2, specP1, specP2, specP3,
                              Ainf, ninf, Einf, A0, n0, E0, T, X, dWdT, dWdrhok);
  add_to_dW_bw_2r3p_Lindemann(specR1, specR2, specP1, specP2, specP3,
                              Ainf, ninf, Einf, A0, n0, E0, T, X, dWdT, dWdrhok);
}




/*
 Note : kf  = k0/(1+k0*X[specR3]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^6 mole^(-2) s^(-1) K^(-n0) 
 Ainf in cm^3 (mole s)^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_fw_3r2p_Lindemann(int specR1, int specR2, int specR3,
                                int specP1, int specP2,
                                double Ainf, double ninf, double Einf, 
                                double A0,  double n0 , double E0, 								
							                	double T,  spec_t X, spec_t W){
  double kf,sum, k0, kinf;
  
  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /* cm^6 mole^(-2) s^(-1) */ 
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem));  /*  cm^3 (mole s)^(-1)   */ 
    kf = k0*kinf/(kinf + k0*X[specR3]);  /* cm^6 mole^(-2) s^(-1) */
  } else {
    kf = 0.0; 
  }
  
  sum=kf*X[specR1]*X[specR2]*X[specR3]; 
  W[specR1]-=_calM(specR1)*sum*1e6;
  W[specR2]-=_calM(specR2)*sum*1e6;
  W[specR3]-=_calM(specR3)*sum*1e6;
  W[specP1]+=_calM(specP1)*sum*1e6;
  W[specP2]+=_calM(specP2)*sum*1e6; 
} 


/*
 Note : kf  = k0/(1+k0*X[specR3]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^6 mole^(-2) s^(-1) K^(-n0) 
 Ainf in cm^3 (mole s)^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_bw_3r2p_Lindemann(int specR1, int specR2, int specR3,
                                int specP1, int specP2,
                                double Ainf, double ninf, double Einf, 
                                double A0,  double n0 , double E0, 								
                                double T, spec_t X, spec_t W){
  double kf,dG0,Kc,kb,sum,k0,kinf;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /* cm^6 mole^(-2) s^(-1) */ 
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem));  /*  cm^3 (mole s)^(-1)   */ 
    kf = k0*kinf/(kinf + k0*X[specR3]);  /* cm^6 mole^(-2) s^(-1) */
  } else {
    kf = 0.0; 
  }

  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T)
     -_Gs(specR3,T); /* J/mole */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* cm^3/mole */
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


/*
 Note : kf  = k0/(1+k0*X[specR3]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 A0 in cm^6 mole^(-2) s^(-1) K^(-n0) 
 Ainf in cm^3 (mole s)^(-1) K^(-ninf)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 W in kg m^(-3) s^(-1)
*/ 
void add_to_W_fwbw_3r2p_Lindemann(int specR1, int specR2, int specR3,
                                  int specP1, int specP2,
                                  double Ainf, double ninf, double Einf, 
                                  double A0,  double n0 , double E0, 								
                                  double T, spec_t X, spec_t W){
  add_to_W_fw_3r2p_Lindemann(specR1, specR2, specR3, specP1, specP2, Ainf, ninf, Einf, A0, n0, E0, T, X, W);
  add_to_W_bw_3r2p_Lindemann(specR1, specR2, specR3, specP1, specP2, Ainf, ninf, Einf, A0, n0, E0, T, X, W);
} 




/*
 Note : kf  = k0/(1+k0*X[specR3]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 Ainf in cm^3 (mole s)^(-1) K^(-ninf) 
 A0 in cm^6 mole^(-2) s^(-1) K^(-n0)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_fw_3r2p_Lindemann(int specR1, int specR2, int specR3,
                                 int specP1, int specP2,
                                 double Ainf, double ninf, double Einf, 
                                 double A0, double n0, double E0,
                                 double T, spec_t X, 
                                 spec_t dWdT, spec2_t dWdrhok){
									 
  double kf, sum, k0, kinf, dkfdT, dk0dT, dkinfdT, dkfdX3;

  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /* cm^6 mole^(-2) s^(-1) */ 
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem));  /*  cm^3 (mole s)^(-1)   */ 
    kf = k0*kinf/(kinf + k0*X[specR3]);  /* cm^6 mole^(-2) s^(-1) */

    dk0dT  = k0*n0/T + k0*(E0/(sqr(T)*Rchem));
    dkinfdT = kinf*ninf/T + kinf*(Einf/(sqr(T)*Rchem));
  
    dkfdX3 = -(k0 * k0 * kinf) / sqr(kinf + k0*X[specR3]); /*  cm^6 (mole^2 s)^(-1)   */ 
 
    dkfdT = kinf*dk0dT / ( kinf + k0*X[specR3] ) - k0 * ( X[specR3]*kinf*dk0dT - X[specR3]*k0*dkinfdT ) / sqr(kinf + k0*X[specR3]);
    
    
  } else {
    kf = 0.0; 
    dkfdT = 0.0; 
    dkfdX3 = 0.0; 
  }
  
  // dWdT 
  
  sum=X[specR1]*X[specR2]*X[specR3]; 
  
  dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum ;
  dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum ;
  dWdT[specR3] -= _calM(specR3)*dkfdT*1e6*sum ;
  dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum ;
  dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum ; 

  /* dW[specR1]dX */
  dWdrhok[specR1][specR1]-=_calM(specR1)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specR1][specR2]-=_calM(specR1)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specR1][specR3]-=_calM(specR1)/_calM(specR3)*kf*X[specR1]*X[specR2] + _calM(specR1)/_calM(specR3) * sum * dkfdX3;

  /* dW[specR2]dX */
  dWdrhok[specR2][specR1]-=_calM(specR2)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specR2][specR2]-=_calM(specR2)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specR2][specR3]-=_calM(specR2)/_calM(specR3)*kf*X[specR1]*X[specR2] + _calM(specR2)/_calM(specR3) * sum * dkfdX3;

  /* dW[specR3]dX */
  dWdrhok[specR3][specR1]-=_calM(specR3)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specR3][specR2]-=_calM(specR3)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specR3][specR3]-=_calM(specR3)/_calM(specR3)*kf*X[specR1]*X[specR2] + _calM(specR3)/_calM(specR3) * sum * dkfdX3;
  
  /* dW[specP1]dX */
  dWdrhok[specP1][specR1]+=_calM(specP1)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specP1][specR2]+=_calM(specP1)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specP1][specR3]+=_calM(specP1)/_calM(specR3)*kf*X[specR1]*X[specR2] + _calM(specP1)/_calM(specR3) * sum * dkfdX3;

  /* dW[specP2]dX */
  dWdrhok[specP2][specR1]+=_calM(specP2)/_calM(specR1)*kf*X[specR2]*X[specR3];
  dWdrhok[specP2][specR2]+=_calM(specP2)/_calM(specR2)*kf*X[specR1]*X[specR3];
  dWdrhok[specP2][specR3]+=_calM(specP2)/_calM(specR3)*kf*X[specR1]*X[specR2] + _calM(specP2)/_calM(specR3) * sum * dkfdX3;
  
} 


/*
 Note : kf  = k0/(1+k0*X[specR3]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 Ainf in cm^3 (mole s)^(-1) K^(-ninf) 
 A0 in cm^6 mole^(-2) s^(-1) K^(-n0)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_bw_3r2p_Lindemann(int specR1, int specR2, int specR3,
                                 int specP1, int specP2,
                                 double Ainf, double ninf, double Einf, 
                                 double A0, double n0, double E0,
                                 double T, spec_t X, 
                                 spec_t dWdT, spec2_t dWdrhok){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT,k0,kinf,dk0dT,dkinfdT,dkfdX3;


  if (A0 > 0.0 || Ainf > 0.0) {
    k0  = A0*pow(T,n0)*exp(-E0/(T*Rchem));    /* cm^6 mole^(-2) s^(-1) */ 
    kinf = Ainf*pow(T,ninf)*exp(-Einf/(T*Rchem));  /*  cm^3 (mole s)^(-1)   */ 
    kf = k0*kinf/(kinf + k0*X[specR3]);  /* cm^6 mole^(-2) s^(-1) */
  
    dk0dT  = k0*n0/T + k0*(E0/(sqr(T)*Rchem));
    dkinfdT = kinf*ninf/T + kinf*(Einf/(sqr(T)*Rchem));
  
    dkfdX3 = -(k0 * k0 * kinf) / sqr(kinf + k0*X[specR3]); /*  cm^6 (mole^2 s)^(-1)   */ 
 
    dkfdT = kinf*dk0dT / ( kinf + k0*X[specR3] ) - k0 * ( X[specR3]*kinf*dk0dT - X[specR3]*k0*dkinfdT ) / sqr(kinf + k0*X[specR3]);
    
  } else {
    kf = 0.0; 
    dkfdT = 0.0; 
    dkfdX3 = 0.0; 
  }

  dG0=_Gs(specP1,T)
     +_Gs(specP2,T)
     -_Gs(specR1,T)
     -_Gs(specR2,T)
     -_Gs(specR3,T); /* J/mole */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* cm^3/mole */

  sum=-X[specP1]*X[specP2]/Kc; /* mole^3 cm^(-9) */

  /* dWdT */
  
  dG0dT=_dGsdT(specP1,T)+_dGsdT(specP2,T)-_dGsdT(specR1,T)-_dGsdT(specR2,T)-_dGsdT(specR3,T);
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T)+1.0/T;

  dsumdT=X[specP1]*X[specP2]/Kc*dKcdToverKc;
  
  if (exp(-dG0/(calR*T))>Kcmin){

    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specR3] -= _calM(specR3)*dkfdT*1e6*sum + _calM(specR3)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
    dWdT[specP2] += _calM(specP2)*dkfdT*1e6*sum + _calM(specP2)*kf*1e6*dsumdT; 

    /* dW[specR1]dX */
    dWdrhok[specR1][specP1]-=-_calM(specR1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR1][specP2]-=-_calM(specR1)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specR1)/_calM(specP2)*dkfdX3*X[specP1]*X[specP2]/Kc;

    /* dW[specR2]dX */
    dWdrhok[specR2][specP1]-=-_calM(specR2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR2][specP2]-=-_calM(specR2)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specR2)/_calM(specP2)*dkfdX3*X[specP1]*X[specP2]/Kc;

    /* dW[specR3]dX */
    dWdrhok[specR3][specP1]-=-_calM(specR3)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specR3][specP2]-=-_calM(specR3)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specR3)/_calM(specP2)*dkfdX3*X[specP1]*X[specP2]/Kc;

    /* dW[specP1]dX */
    dWdrhok[specP1][specP1]+=-_calM(specP1)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP1][specP2]+=-_calM(specP1)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specP1)/_calM(specP2)*dkfdX3*X[specP1]*X[specP2]/Kc;

    /* dW[specP2]dX */
    dWdrhok[specP2][specP1]+=-_calM(specP2)/_calM(specP1)*kf*X[specP2]/Kc;
    dWdrhok[specP2][specP2]+=-_calM(specP2)/_calM(specP2)*kf*X[specP1]/Kc - _calM(specP2)/_calM(specP2)*dkfdX3*X[specP1]*X[specP2]/Kc;
  }
} 


/*
 Note : kf  = k0/(1+k0*X[specR3]/kinf)
        with kinf determined from Ainf, ninf, Einf
             k0 determined from A0, n0, E0
 Units:
 Ainf in cm^3 (mole s)^(-1) K^(-ninf) 
 A0 in cm^6 mole^(-2) s^(-1) K^(-n0)
 Einf,E0 in cal/mole
 T in Kelvin
 X in mole/cm3
 dWdT in kg m^(-3) s^(-1) K^(-1)
 dWdrhok in s^(-1)
*/ 
void add_to_dW_fwbw_3r2p_Lindemann(int specR1, int specR2, int specR3,
                                   int specP1, int specP2,
                                   double Ainf, double ninf, double Einf, 
                                   double A0, double n0, double E0,
                                   double T, spec_t X, 
                                   spec_t dWdT, spec2_t dWdrhok){
  add_to_dW_fw_3r2p_Lindemann(specR1, specR2, specR3, specP1, specP2,
                              Ainf, ninf, Einf, A0, n0, E0, T, X, dWdT, dWdrhok);
  add_to_dW_bw_3r2p_Lindemann(specR1, specR2, specR3, specP1, specP2,
                              Ainf, ninf, Einf, A0, n0, E0, T, X, dWdT, dWdrhok);
}




void test_dW_dx(gl_t *gl, spec_t rhokref, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam){
  long r,s; 
  spec_t dWdT,dWdTe,dWdTv,dWdQbeam;
  spec2_t dWrhok,dWrhok2;   
  spec_t rhok2,dWdT2,W,W2, X;
  double drhok,dQbeam,dT,N,N2,Estar2;
 
  find_W(gl, rhok, T, Te, Tv, Estar, Qbeam, W);      
  find_dW_dx(gl, rhok, T, Te, Tv, Estar, Qbeam,
              dWrhok, dWdT, dWdTe, dWdTv, dWdQbeam);

  N=0.0;
  for (s=0; s<ns; s++) N+=rhok[s]/_m(s);   
  for (s=0; s<ns; s++) X[s]= rhok[s] / _calM ( s ) * 1.0e-06;   
  for (s=0; s<ns; s++) wfprintf(stdout,"X[%ld]=%E mol/cm3\n",s,X[s]);  
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
  R=Rchem;
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
  R=Rchem;
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


#ifdef speceminus
/* add contributions to Qei coming from one electron impact chemical reaction
 * spec is the neutral molecule the electrons impact
 * exci is the excitation energy in eV
 * kf is the reaction rate in cm3/s
 * rhok is the partial densities in kg/m3
 * Qei is the heat removed from the electrons in J/m3
 */
void add_to_Qei(long spec, double exci, double kf, spec_t rhok, double *Qei){
#ifdef TEST
  kf=1e-17;
#endif
  (*Qei) += kf * 1E-6* exci * echarge * rhok[speceminus] / _calM ( speceminus ) * rhok[spec] / _calM ( spec ) * sqr(calA);
}


/* add contributions to dQei_dx coming from one electron impact chemical reaction
 * spec is the neutral molecule the electrons impact
 * exci is the excitation energy in eV
 * kf is the reaction rate in cm3/s
 * dkfdTe is the derivative of kf with respect to Te in cm3/(s K)
 * rhok is the partial densities in kg/m3
 * dQeidrhok is in J/kg 
 * dQeidTe is in J/(m3 K)
 */
void add_to_dQei(long spec, double exci, double kf, double dkfdTe, spec_t rhok, spec_t dQeidrhok, double *dQeidTe){
#ifdef TEST
  kf=1e-17;
#endif
  dQeidrhok[spec] += kf *1e-6* exci * echarge * rhok[speceminus] / _calM ( speceminus )  / _calM ( spec ) * sqr(calA);
  dQeidrhok[speceminus] += kf *1e-6* exci * echarge / _calM ( speceminus ) * rhok[spec] / _calM ( spec ) * sqr(calA);
}


/* find excitation energy in eV using the difference in enthalpy of formation */
double _exci_2r2p(long react1, long react2, long prod1, long prod2){
  double exci;
  // find excitation energy needed per one particule in J 
  exci=_hk_from_T_equilibrium(prod1, TREF_EXCI)*_m(prod1)+_hk_from_T_equilibrium(prod2, TREF_EXCI)*_m(prod2)
      -_hk_from_T_equilibrium(react1, TREF_EXCI)*_m(react1)-_hk_from_T_equilibrium(react2, TREF_EXCI)*_m(react2);
  // convert J to eV
  exci=exci*6.242e18;
  return(exci);
}


/* find excitation energy in eV using the difference in enthalpy of formation */
double _exci_2r3p(long react1, long react2, long prod1, long prod2, long prod3){
  double exci;
  // find excitation energy needed per one particule in J 
  exci=_hk_from_T_equilibrium(prod1, TREF_EXCI)*_m(prod1)+_hk_from_T_equilibrium(prod2, TREF_EXCI)*_m(prod2)
      +_hk_from_T_equilibrium(prod3, TREF_EXCI)*_m(prod3)
      -_hk_from_T_equilibrium(react1, TREF_EXCI)*_m(react1)-_hk_from_T_equilibrium(react2, TREF_EXCI)*_m(react2);
  // convert J to eV
  exci=exci*6.242e18;
  return(exci);
}

/* find excitation energy in eV using the difference in enthalpy of formation */
double _exci_2r4p(long react1, long react2, long prod1, long prod2, long prod3, long prod4){
  double exci;
  // find excitation energy needed per one particule in J 
  exci=_hk_from_T_equilibrium(prod1, TREF_EXCI)*_m(prod1)+_hk_from_T_equilibrium(prod2, TREF_EXCI)*_m(prod2)
      +_hk_from_T_equilibrium(prod3, TREF_EXCI)*_m(prod3)+_hk_from_T_equilibrium(prod4, TREF_EXCI)*_m(prod4)
      -_hk_from_T_equilibrium(react1, TREF_EXCI)*_m(react1)-_hk_from_T_equilibrium(react2, TREF_EXCI)*_m(react2);
  // convert J to eV
  exci=exci*6.242e18;
  return(exci);
}

#endif
