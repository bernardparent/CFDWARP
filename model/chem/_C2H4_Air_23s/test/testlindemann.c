#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define calR 8.314472e0  /* J/(K mol) */
#define calA 6.02214199E23  /* particules per mole */
#define Rchem 1.987192004e0  /*cal/(K mol)*/
#define echarge 1.602176462e-19 /* C */
#define emass 9.10938188E-31  /* kg */
#define epsilon0 8.854187817E-12 /* permittivity of free space */
#define kB 1.3806488E-23   /* m2*kg/(s2*K) */


#define specCH3 0
#define specH 1
#define specM 2
#define specCH4 3

#define ns 4


typedef double spec_t[ns];


const static double calM[ns]=
  {
   15.03452E-3,   /* CH3*/
   1.00794E-3,    /*H */
   28.0e-3,    /*M */
   16.04246E-3,   /* CH4 */
};

double _calM(long spec){
  double tmp;
  tmp=calM[spec];
  return(tmp);
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
  W[specR1]-=_calM(specR1)*sum*1.0e6;
  W[specR2]-=_calM(specR2)*sum*1.0e6;
  W[specR3]-=_calM(specR3)*sum*1.0e6;
  W[specP1]+=_calM(specP1)*sum*1.0e6;
  W[specP2]+=_calM(specP2)*sum*1.0e6; 
} 


void main(){
  spec_t X,W;
  double eff,T;
  long spec;
  
  /* test the reaction
CH3 + H (+ M) => CH4 (+ M)     5.20000E+12      0      -1310     	
			  LOW/ 8.25000E+14      0      -19310 /    
C2H4/ 3.0/ CH4/ 6.5/ CO/ 0.75/ CO2/ 1.5/ H2/ 1.0/ H2O/ 6.5/ N2/ 0.4/ O2/ 0.4
*/

  T=1000.0; /* K */
  eff=0.5;
  for (spec=0; spec<ns; spec++) W[spec]=0.0;
  X[0]=1e-6; /* moles per cm3 */
  X[1]=1e-6;
  X[2]=1e-6;
  X[3]=1e-6;
  add_to_W_fw_3r2p_Lindemann(specCH3, specH, specM, specCH4, specM, eff*5.20000E+12, 0.0, -1310.0, eff*8.25000E+14, 0.0, -19310.0, T, X, W);
  for (spec=0; spec<ns; spec++) printf("W[ld]=%E\n",W[spec]);
}
