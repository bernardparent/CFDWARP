#include <model/share/chem_share.h>
#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>

#define Kcmin 1.0e-40


/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_W_Arrhenius2r2p(int specR1, int specR2,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     +Gs[specP2]
     -Gs[specR1]
     -Gs[specR2]; /* J/gmol */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  kb=kf/Kc; /* cm^3 (gmol s)^(-1) */
  sum=kf*X[specR1]*X[specR2]
     -kb*X[specP1]*X[specP2]; /* gmol cm^(-3) (s)^(-1) */
  if (exp(-dG0/(calR*T))>Kcmin){
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6;
    W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
  }
} 


/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_W_Arrhenius3r2p(int specR1, int specR2, int specR3,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     +Gs[specP2]
     -Gs[specR1]
     -Gs[specR2]
     -Gs[specR3]; /* J/gmol */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* gmol/cm^3 */
  kb=kf/Kc; /* cm^3 (gmol s)^(-1) */
  sum=kf*X[specR1]*X[specR2]*X[specR3]
     -kb*X[specP1]*X[specP2]; /* gmol cm^(-3) (s)^(-1) */
  if (exp(-dG0/(calR*T))>Kcmin){
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specR3]-=_calM(specR3)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6;
    W[specP2]+=_calM(specP2)*sum*1e6; /* kg m^(-3) s^(-1) */
  }
} 



/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_W_Arrhenius2r1p(int specR1, int specR2,
                            int specP1,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     -Gs[specR1]
     -Gs[specR2]; /* J/gmol */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* gmol/cm^3 */
  kb=kf/Kc; /* cm^3 (gmol s)^(-1) */
  sum=kf*X[specR1]*X[specR2]
     -kb*X[specP1]; /* gmol cm^(-3) (s)^(-1) */
  if (exp(-dG0/(calR*T))>Kcmin){
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6; /* kg m^(-3) s^(-1) */
  }
} 


/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_W_Arrhenius2r3p(int specR1, int specR2,
                            int specP1, int specP2, int specP3,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W){
  double kf,dG0,Kc,kb,sum;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     +Gs[specP2]
     +Gs[specP3]
     -Gs[specR1]
     -Gs[specR2]; /* J/gmol */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  if (exp(-dG0/(calR*T))>Kcmin){
    Kc *= 1E-6 * THERMO_P_REF/(calR*T); /* gmol/cm^3 */
    kb=kf/Kc; /*  cm^6 (gmol)^(-2) s^(-1) */
    sum=kf*X[specR1]*X[specR2]
     -kb*X[specP1]*X[specP2]*X[specP3]; /* gmol cm^(-3) (s)^(-1) */
    W[specR1]-=_calM(specR1)*sum*1e6;
    W[specR2]-=_calM(specR2)*sum*1e6;
    W[specP1]+=_calM(specP1)*sum*1e6;
    W[specP2]+=_calM(specP2)*sum*1e6; 
    W[specP3]+=_calM(specP3)*sum*1e6; 
    /* kg m^(-3) s^(-1) */
  }
}



/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_dW_Arrhenius2r2p(int specR1, int specR2,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     +Gs[specP2]
     -Gs[specR1]
     -Gs[specR2]; /* J/gmol */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  sum=X[specR1]*X[specR2]
     -X[specP1]*X[specP2]/Kc; /* gmol cm^(-3) (s)^(-1) */
/*  W[specR1]-=_calM(specR1)*kf*1e6*sum;
  W[specR2]-=_calM(specR2)*kf*1e6*sum;
  W[specP1]+=_calM(specP1)*kf*1e6*sum;
  W[specP2]+=_calM(specP2)*kf*1e6*sum;*/ /* kg m^(-3) s^(-1) */
  /* dWdT */
  dG0dT=dGsdT[specP1]+dGsdT[specP2]-dGsdT[specR1]-dGsdT[specR2];
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
    dWdX[specR1][specR1]-=_calM(specR1)*kf*1e6*X[specR2];
    dWdX[specR1][specR2]-=_calM(specR1)*kf*1e6*X[specR1];
    dWdX[specR1][specP1]-=-_calM(specR1)*kf*1e6*X[specP2]/Kc;
    dWdX[specR1][specP2]-=-_calM(specR1)*kf*1e6*X[specP1]/Kc;

    /* dW[specR2]dX */
    dWdX[specR2][specR1]-=_calM(specR2)*kf*1e6*X[specR2];
    dWdX[specR2][specR2]-=_calM(specR2)*kf*1e6*X[specR1];
    dWdX[specR2][specP1]-=-_calM(specR2)*kf*1e6*X[specP2]/Kc;
    dWdX[specR2][specP2]-=-_calM(specR2)*kf*1e6*X[specP1]/Kc;

    /* dW[specP1]dX */
    dWdX[specP1][specR1]+=_calM(specP1)*kf*1e6*X[specR2];
    dWdX[specP1][specR2]+=_calM(specP1)*kf*1e6*X[specR1];
    dWdX[specP1][specP1]+=-_calM(specP1)*kf*1e6*X[specP2]/Kc;
    dWdX[specP1][specP2]+=-_calM(specP1)*kf*1e6*X[specP1]/Kc;

    /* dW[specP2]dX */
    dWdX[specP2][specR1]+=_calM(specP2)*kf*1e6*X[specR2];
    dWdX[specP2][specR2]+=_calM(specP2)*kf*1e6*X[specR1];
    dWdX[specP2][specP1]+=-_calM(specP2)*kf*1e6*X[specP2]/Kc;
    dWdX[specP2][specP2]+=-_calM(specP2)*kf*1e6*X[specP1]/Kc;
  }
} 




/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_dW_Arrhenius3r2p(int specR1, int specR2, int specR3,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     +Gs[specP2]
     -Gs[specR1]
     -Gs[specR2]
     -Gs[specR3]; /* J/gmol */
  Kc=max(Kcmin,exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* gmol/cm^3 */

  sum=X[specR1]*X[specR2]*X[specR3]
     -X[specP1]*X[specP2]/Kc; /* gmol cm^(-3) (s)^(-1) */
/*  W[specR1]-=_calM(specR1)*kf*1e6*sum;
  W[specR2]-=_calM(specR2)*kf*1e6*sum;
  W[specP1]+=_calM(specP1)*kf*1e6*sum;
  W[specP2]+=_calM(specP2)*kf*1e6*sum;*/ /* kg m^(-3) s^(-1) */

  /* dWdT */
  dG0dT=dGsdT[specP1]+dGsdT[specP2]-dGsdT[specR1]-dGsdT[specR2]-dGsdT[specR3];
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
    dWdX[specR1][specR1]-=_calM(specR1)*kf*1e6*X[specR2]*X[specR3];
    dWdX[specR1][specR2]-=_calM(specR1)*kf*1e6*X[specR1]*X[specR3];
    dWdX[specR1][specR3]-=_calM(specR1)*kf*1e6*X[specR1]*X[specR2];
    dWdX[specR1][specP1]-=-_calM(specR1)*kf*1e6*X[specP2]/Kc;
    dWdX[specR1][specP2]-=-_calM(specR1)*kf*1e6*X[specP1]/Kc;

    /* dW[specR2]dX */
    dWdX[specR2][specR1]-=_calM(specR2)*kf*1e6*X[specR2]*X[specR3];
    dWdX[specR2][specR2]-=_calM(specR2)*kf*1e6*X[specR1]*X[specR3];
    dWdX[specR2][specR3]-=_calM(specR2)*kf*1e6*X[specR1]*X[specR2];
    dWdX[specR2][specP1]-=-_calM(specR2)*kf*1e6*X[specP2]/Kc;
    dWdX[specR2][specP2]-=-_calM(specR2)*kf*1e6*X[specP1]/Kc;

    /* dW[specR3]dX */
    dWdX[specR3][specR1]-=_calM(specR3)*kf*1e6*X[specR2]*X[specR3];
    dWdX[specR3][specR2]-=_calM(specR3)*kf*1e6*X[specR1]*X[specR3];
    dWdX[specR3][specR3]-=_calM(specR3)*kf*1e6*X[specR1]*X[specR2];
    dWdX[specR3][specP1]-=-_calM(specR3)*kf*1e6*X[specP2]/Kc;
    dWdX[specR3][specP2]-=-_calM(specR3)*kf*1e6*X[specP1]/Kc;

    /* dW[specP1]dX */
    dWdX[specP1][specR1]+=_calM(specP1)*kf*1e6*X[specR2]*X[specR3];
    dWdX[specP1][specR2]+=_calM(specP1)*kf*1e6*X[specR1]*X[specR3];
    dWdX[specP1][specR3]+=_calM(specP1)*kf*1e6*X[specR1]*X[specR2];
    dWdX[specP1][specP1]+=-_calM(specP1)*kf*1e6*X[specP2]/Kc;
    dWdX[specP1][specP2]+=-_calM(specP1)*kf*1e6*X[specP1]/Kc;

    /* dW[specP2]dX */
    dWdX[specP2][specR1]+=_calM(specP2)*kf*1e6*X[specR2]*X[specR3];
    dWdX[specP2][specR2]+=_calM(specP2)*kf*1e6*X[specR1]*X[specR3];
    dWdX[specP2][specR3]+=_calM(specP2)*kf*1e6*X[specR1]*X[specR2];
    dWdX[specP2][specP1]+=-_calM(specP2)*kf*1e6*X[specP2]/Kc;
    dWdX[specP2][specP2]+=-_calM(specP2)*kf*1e6*X[specP1]/Kc;
  }
} 




/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_dW_Arrhenius2r1p(int specR1, int specR2,
                             int specP1, 
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     -Gs[specR1]
     -Gs[specR2]; /* J/gmol */
  Kc=max(Kcmin, exp(-dG0/(calR*T)));
  Kc /= 1E-6 * THERMO_P_REF/(calR*T); /* gmol/cm^3 */

  sum=X[specR1]*X[specR2]
     -X[specP1]/Kc; /* gmol cm^(-3) (s)^(-1) */

  /* dWdT */
  dG0dT=dGsdT[specP1]-dGsdT[specR1]-dGsdT[specR2];
  
  dKcdToverKc=(dG0/(calR*sqr(T)))-dG0dT/(calR*T)+1.0/T;

  dsumdT=X[specP1]/Kc*dKcdToverKc;
  dkfdT= kf*n/T 
       + kf*(E/(sqr(T)*1.987192004e0));

  if (exp(-dG0/(calR*T))>Kcmin){
    dWdT[specR1] -= _calM(specR1)*dkfdT*1e6*sum + _calM(specR1)*kf*1e6*dsumdT;
    dWdT[specR2] -= _calM(specR2)*dkfdT*1e6*sum + _calM(specR2)*kf*1e6*dsumdT;
    dWdT[specP1] += _calM(specP1)*dkfdT*1e6*sum + _calM(specP1)*kf*1e6*dsumdT;
  

    /* dW[specR1]dX */
    dWdX[specR1][specR1]-=_calM(specR1)*kf*1e6*X[specR2];
    dWdX[specR1][specR2]-=_calM(specR1)*kf*1e6*X[specR1];
    dWdX[specR1][specP1]-=-_calM(specR1)*kf*1e6/Kc;

    /* dW[specR2]dX */
    dWdX[specR2][specR1]-=_calM(specR2)*kf*1e6*X[specR2];
    dWdX[specR2][specR2]-=_calM(specR2)*kf*1e6*X[specR1];
    dWdX[specR2][specP1]-=-_calM(specR2)*kf*1e6/Kc;

    /* dW[specP1]dX */
    dWdX[specP1][specR1]+=_calM(specP1)*kf*1e6*X[specR2];
    dWdX[specP1][specR2]+=_calM(specP1)*kf*1e6*X[specR1];
    dWdX[specP1][specP1]+=-_calM(specP1)*kf*1e6/Kc;
  }

} 

/* A in cm^3 (gmol s)^(-1) K^(-n) 
   E in cal gmol^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/
void add_to_dW_Arrhenius2r3p(int specR1, int specR2,
                             int specP1, int specP2, int specP3,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX){
  double kf,dG0,Kc,sum,dkfdT,dsumdT,dKcdToverKc,dG0dT;

  kf=A*pow(T,n)*exp(-E/(T*1.987192004e0));  /* cm^3 (gmol s)^(-1) */
  dG0=Gs[specP1]
     +Gs[specP2]
     +Gs[specP3]
     -Gs[specR1]
     -Gs[specR2]; /* J/gmol */
  Kc=1E-6 * THERMO_P_REF/(calR*T) * max(Kcmin, exp(-dG0/(calR*T))); /* gmol/cm^3 */
  
  sum=X[specR1]*X[specR2]
     -X[specP1]*X[specP2]*X[specP3]/Kc; /* gmol cm^(-3) (s)^(-1) */
/*  W[specR1]-=_calM(specR1)*kf*1e6*sum;
  W[specR2]-=_calM(specR2)*kf*1e6*sum;
  W[specP1]+=_calM(specP1)*kf*1e6*sum;
  W[specP2]+=_calM(specP2)*kf*1e6*sum;
  W[specP3]+=_calM(specP3)*kf*1e6*sum;*/
 /* kg m^(-3) s^(-1) */

  // dWdT 
  dG0dT=dGsdT[specP1]+dGsdT[specP2]+dGsdT[specP3]-dGsdT[specR1]-dGsdT[specR2];
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
    dWdX[specR1][specR1]-=_calM(specR1)*kf*1e6*X[specR2];
    dWdX[specR1][specR2]-=_calM(specR1)*kf*1e6*X[specR1];
    dWdX[specR1][specP1]-=-_calM(specR1)*kf*1e6*X[specP2]*X[specP3]/Kc;
    dWdX[specR1][specP2]-=-_calM(specR1)*kf*1e6*X[specP1]*X[specP3]/Kc;
    dWdX[specR1][specP3]-=-_calM(specR1)*kf*1e6*X[specP1]*X[specP2]/Kc;

    // dW[specR2]dX 
    dWdX[specR2][specR1]-=_calM(specR2)*kf*1e6*X[specR2];
    dWdX[specR2][specR2]-=_calM(specR2)*kf*1e6*X[specR1];
    dWdX[specR2][specP1]-=-_calM(specR2)*kf*1e6*X[specP2]*X[specP3]/Kc;
    dWdX[specR2][specP2]-=-_calM(specR2)*kf*1e6*X[specP1]*X[specP3]/Kc;
    dWdX[specR2][specP3]-=-_calM(specR2)*kf*1e6*X[specP1]*X[specP2]/Kc;

    // dW[specP1]dX 
    dWdX[specP1][specR1]+=_calM(specP1)*kf*1e6*X[specR2];
    dWdX[specP1][specR2]+=_calM(specP1)*kf*1e6*X[specR1];
    dWdX[specP1][specP1]+=-_calM(specP1)*kf*1e6*X[specP2]*X[specP3]/Kc;
    dWdX[specP1][specP2]+=-_calM(specP1)*kf*1e6*X[specP1]*X[specP3]/Kc;
    dWdX[specP1][specP3]+=-_calM(specP1)*kf*1e6*X[specP1]*X[specP2]/Kc;

    // dW[specP2]dX 
    dWdX[specP2][specR1]+=_calM(specP2)*kf*1e6*X[specR2];
    dWdX[specP2][specR2]+=_calM(specP2)*kf*1e6*X[specR1];
    dWdX[specP2][specP1]+=-_calM(specP2)*kf*1e6*X[specP2]*X[specP3]/Kc;
    dWdX[specP2][specP2]+=-_calM(specP2)*kf*1e6*X[specP1]*X[specP3]/Kc;
    dWdX[specP2][specP3]+=-_calM(specP2)*kf*1e6*X[specP1]*X[specP2]/Kc;

    // dW[specP3]dX 
    dWdX[specP3][specR1]+=_calM(specP3)*kf*1e6*X[specR2];
    dWdX[specP3][specR2]+=_calM(specP3)*kf*1e6*X[specR1];
    dWdX[specP3][specP1]+=-_calM(specP3)*kf*1e6*X[specP2]*X[specP3]/Kc;
    dWdX[specP3][specP2]+=-_calM(specP3)*kf*1e6*X[specP1]*X[specP3]/Kc;
    dWdX[specP3][specP3]+=-_calM(specP3)*kf*1e6*X[specP1]*X[specP2]/Kc; 
  }
}


void test_dW_dx(spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam){
  long r,s; 
  spec_t dWdT,dWdTe,dWdTv,dWdQbeam;
  spec2_t dWrhok,dWrhok2;   
  spec_t rhok2,dWdT2,W,W2;
  double dQbeam,dT,N,N2,Estar2;
 
  find_W(rhok, T, Te, Tv, Estar, Qbeam, W);      
  find_dW_dx(rhok, mu, T, Te, Tv, Estar, Qbeam,
              dWrhok, dWdT, dWdTe, dWdTv, dWdQbeam);

  N=0.0;
  for (s=0; s<ns; s++) N+=rhok[s]/_m(s);   

  for (s=0; s<ns; s++) wfprintf(stdout,"rho[%ld]=%E kg/m3\n",s,rhok[s]);    
  wfprintf(stdout,"T=%E K\nTe=%E K\nTv=%E K\nE/N=%E Vm2\nQbeam=%E W/m3\n",T,Te,Tv,Estar,Qbeam);

  // find dWdQbeam numerically    
   wfprintf(stdout,"\n\ndWdQbeam:\n");
   dQbeam=Qbeam/10000.0+1.0e2;
   find_W(rhok, T, Te,Tv,Estar, Qbeam+dQbeam, W2);
   for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dQbeam;
   for (s=0; s<ns; s++){  
     wfprintf(stdout,"%15.15E  %15.15E\n",dWdQbeam[s],dWdT2[s]);
   } 
      
  /* find dWrhok numerically */
  for (s=0; s<ns; s++){
    for (r=0; r<ns; r++) rhok2[r]=rhok[r];
    rhok2[s]+=rhok[s]/1e3+1.0e-10;
    N2=N+(rhok2[s]-rhok[s])/_m(s);
    Estar2=Estar*N/N2;
    find_W(rhok2, T, Te,Tv,Estar2,Qbeam, W2);
    for (r=0; r<ns; r++) dWrhok2[r][s]=(W2[r]-W[r])/(rhok[s]/1e3+1.0e-10);
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
   dT=T/100000.0;
  find_W(rhok, T+dT, Te,Tv,Estar,Qbeam, W2);
  for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dT;
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdT[s],dWdT2[s]);
  } 
  
  // find dWdTe numerically 
   wfprintf(stdout,"\n\ndWdTe:\n");
  dT=Te/100000.0;
  find_W(rhok, T, Te+dT,Tv,Estar,Qbeam, W2);
  for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dT;
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdTe[s],dWdT2[s]);
  }
    
  exit(1);  
}

