#include <exm.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#define echarge 1.60217646E-19 /* elementary charge in Coulomb */   
#define kB 1.3806503e-23     /* Boltzman constant */
#define epsilon0 8.854188E-12 /* permittivity of free space */


/* drift velocity of nitrogen in m/s
   from Ch. 21 of Handbook of Physical Quantities by Grigoriev, I. S. and Meilikhov, E. Z. 1997
   valid range: 3e-20<Estar<100E-20 Vm2
   can be used until Estar=240e-20 Vm2 with perhaps not so bad accuracy*/
double _vdr_N2(double Estar){
  double vdr;

vdr= 9.39370742160896043060E+03
+  1.01176190755016803982E+24*Estar
- 5.14736796346349010349E+41*sqr(Estar)
+ 1.01087741046674631584E+59*Estar*sqr(Estar);
  return(vdr);
}


/* drift velocity of oxygen in m/s
   from Ch. 21 of Handbook of Physical Quantities by Grigoriev, I. S. and Meilikhov, E. Z. 1997
   valid range: 3e-20<Estar<100E-20 Vm2
   can be used until Estar=240e-20 Vm2 with perhaps not so bad accuracy */
double _vdr_O2(double Estar){
  double vdr;


vdr=  2.75467534779719135258E+04
+  1.15609508977947801813E+24*Estar
 -5.74009130605601932491E+41*sqr(Estar)
+ 1.03418565951837087954E+59*Estar*sqr(Estar);

  return(vdr);
}


/* townsend ionization coefficient for N2 in cm3/s 
   from A. K. Mnatsakanyan, G. V. Naidis, Processes of Formation and Decay of Charged Particles in Nitrogen-Oxygen Plasmas, in: B. M. Smirnov (Ed.), Khimiia Plazmy [Plasma Chemistry], vol. 14, Energoatomizdat, Moscow, Russia, 227–255, in russian, 1987.
   valid range: 3e-20<E/N<30e-20 Vm2   
*/
double _kN2_1(double T, double Estar){
  double kN2;
  kN2=pow(10.0,-8.3-36.5e-20/Estar);  
  return(kN2);
}



/* townsend ionization coefficient for N2 in cm3/s
   from Raizer, Gas discharge physics, p. 56
   valid range: 27<E/P<200 V/(cm*Torr)   at T=300K
   valid range: 20.3<E/P<150.4 V/(m*Pa)   at T=300K
   valid range: 8.4e-20<E/N<62e-20 Vm2   
   
*/
double _kN2_2(double T, double Estar){
  double kN2,vdr;
  double EoverP,alphaoverP;

  EoverP=Estar/kB/T;
  /* change the units from V/(m*Pa) to V/(cm*Torr) */
  EoverP=EoverP*133.0/100.0;
  alphaoverP=8.8*exp(-275.0/EoverP);
  /* change the units from 1/(cm*Torr) to 1/(m*Pa) */
  alphaoverP=alphaoverP/133.0*100.0;
  vdr=_vdr_N2(Estar);
  kN2=alphaoverP*kB*T*vdr*1e6;
  return(kN2);
}


/* townsend ionization coefficient for N2 in cm3/s 
   from Raizer, Gas discharge physics, p. 56
   valid range: 100<E/P<600 V/(cm*Torr)   at T=300K
   valid range: 75.2<E/P<451 V/(m*Pa)   at T=300K
   valid range: 31.1e-20<E/N<187e-20 Vm2   
*/
double _kN2_3(double T, double Estar){
  double kN2,vdr;
  double EoverP,alphaoverP;


  EoverP=Estar/kB/T;
  /* change the units from V/(m*Pa) to V/(cm*Torr) */
  EoverP=EoverP*133.0/100.0;
  alphaoverP=12.0*exp(-342.0/EoverP);
  /* change the units from 1/(cm*Torr) to 1/(m*Pa) */
  alphaoverP=alphaoverP/133.0*100.0;
  vdr=_vdr_N2(Estar);
  kN2=alphaoverP*kB*T*vdr*1e6;  
  return(kN2);
}


/* townsend ionization coefficient for N2 in cm3/s 
   valid range: 3.0e-20<E/N<187e-20 Vm2   
*/
double _kN2_123(double T, double Estar){
  double kN2,kN2_1,kN2_2,kN2_3;
  double fact1,fact2,fact3;

  kN2=0.0; /* to avoid compiler warning */

    kN2_1=_kN2_1(T,Estar);
    kN2_2=_kN2_2(T,Estar);
    kN2_3=_kN2_3(T,Estar);
    if (Estar<8.4e-20) kN2=kN2_1;
    if (Estar>=8.4e-20 && Estar<30e-20) {
      fact2=(Estar-8.4e-20)/(30e-20-8.4e-20);
      fact1=1.0-fact2;
      kN2=fact1*kN2_1+fact2*kN2_2;
    }
    if (Estar>=30e-20 && Estar<62e-20) {
      fact3=(Estar-30e-20)/(62e-20-30e-20);
      fact2=1.0-fact3;
      kN2=fact2*kN2_2+fact3*kN2_3;
    }
    if (Estar>=62e-20) kN2=kN2_3;

  return(kN2);
}




/* townsend ionization coefficient for O2 in cm3/s 
   from A. K. Mnatsakanyan, G. V. Naidis, Processes of Formation and Decay of Charged Particles in Nitrogen-Oxygen Plasmas, in: B. M. Smirnov (Ed.), Khimiia Plazmy [Plasma Chemistry], vol. 14, Energoatomizdat, Moscow, Russia, 227–255, in russian, 1987.
   valid range: 3e-20<E/N<30e-20 Vm2   
*/
double _kO2_1(double T, double Estar){
  double kO2;
  kO2=pow(10.0,-8.8-28.1e-20/Estar);  
  return(kO2);
}


/* townsend ionization coefficient for Air in cm3/s 
   from A. K. Mnatsakanyan, G. V. Naidis, Processes of Formation and Decay of Charged Particles in Nitrogen-Oxygen Plasmas, in: B. M. Smirnov (Ed.), Khimiia Plazmy [Plasma Chemistry], vol. 14, Energoatomizdat, Moscow, Russia, 227–255, in russian, 1987.
   valid range: 3e-20<E/N<30e-20 Vm2   
*/
double _kAir_1(double T, double Estar){
  double kAir;
  kAir=pow(10.0,-8.25-36.5e-20/Estar);  
  return(kAir);
}

/* townsend ionization coefficient for Air in Vm4 
   from A. K. Mnatsakanyan, G. V. Naidis, Processes of Formation and Decay of Charged Particles in Nitrogen-Oxygen Plasmas, in: B. M. Smirnov (Ed.), Khimiia Plazmy [Plasma Chemistry], vol. 14, Energoatomizdat, Moscow, Russia, 227–255, in russian, 1987.
   valid range: 3e-20<E/N<30e-20 Vm2   
*/
double _kovermueN_Air_1(double T, double Estar){
  double kAir,vdr,kovermueN_Air;
  kAir=pow(10.0,-8.25-36.5e-20/Estar);
  vdr=0.8*_vdr_N2(Estar)+0.2*_vdr_O2(Estar);
  kovermueN_Air=kAir/vdr*Estar*1e-6;  
  return(kovermueN_Air);
}


/* townsend ionization coefficient for Air in cm3/s 
   from Raizer, Gas discharge physics, p. 56
   valid range: 44<E/P<176 V/(cm*Torr)   at T=300K
   valid range: 33.1<E/P<132.3 V/(m*Pa)   at T=300K
   valid range: 13.7e-20<E/N<54.8e-20 Vm2   
*/
double _kAir_2(double T, double Estar){
  double kAir,vdr;
  double EoverP,alphaoverP;


  EoverP=Estar/kB/T;
  /* change the units from V/(m*Pa) to V/(cm*Torr) */
  EoverP=EoverP*133.0/100.0;
  alphaoverP=1.17e-4*sqr(EoverP-32.2);
  /* change the units from 1/(cm*Torr) to 1/(m*Pa) */
  alphaoverP=alphaoverP/133.0*100.0;
  vdr=0.8*_vdr_N2(Estar)+0.2*_vdr_O2(Estar);
  kAir=alphaoverP*kB*T*vdr*1e6;  
  return(kAir);
}



/* townsend ionization coefficient for Air in Vm4 
   from Raizer, Gas discharge physics, p. 56
   valid range: 44<E/P<176 V/(cm*Torr)   at T=300K
   valid range: 33.1<E/P<132.3 V/(m*Pa)   at T=300K
   valid range: 13.7e-20<E/N<54.8e-20 Vm2   
*/
double _kovermueN_Air_2(double T, double Estar){
  double kovermueN_Air;
  double EoverP,alphaoverP;


  EoverP=Estar/kB/T;
  /* change the units from V/(m*Pa) to V/(cm*Torr) */
  EoverP=EoverP*133.0/100.0;
  alphaoverP=1.17e-4*sqr(EoverP-32.2);
  /* change the units from 1/(cm*Torr) to 1/(m*Pa) */
  alphaoverP=alphaoverP/133.0*100.0;
  kovermueN_Air=alphaoverP*kB*T*Estar;  
  return(kovermueN_Air);
}


/* townsend ionization coefficient for Air in cm3/s
   from Raizer, Gas discharge physics, p. 56
   valid range: 100<E/P<800 V/(cm*Torr)   at T=300K
   valid range: 75.2<E/P<601.5 V/(m*Pa)   at T=300K
   valid range: 31.1e-20<E/N<249e-20 Vm2   
   
*/
double _kAir_3(double T, double Estar){
  double kAir,vdr;
  double EoverP,alphaoverP;

  EoverP=Estar/kB/T;
  /* change the units from V/(m*Pa) to V/(cm*Torr) */
  EoverP=EoverP*133.0/100.0;
  alphaoverP=15.0*exp(-365.0/EoverP);
  /* change the units from 1/(cm*Torr) to 1/(m*Pa) */
  alphaoverP=alphaoverP/133.0*100.0;
  vdr=0.8*_vdr_N2(Estar)+0.2*_vdr_O2(Estar);
  kAir=alphaoverP*kB*T*vdr*1e6;
  return(kAir);
}


/* townsend ionization coefficient for Air in Vm4
   from Raizer, Gas discharge physics, p. 56
   valid range: 100<E/P<800 V/(cm*Torr)   at T=300K
   valid range: 75.2<E/P<601.5 V/(m*Pa)   at T=300K
   valid range: 31.1e-20<E/N<249e-20 Vm2   
   
*/
double _kovermueN_Air_3(double T, double Estar){
  double kovermueN_Air;
  double EoverP,alphaoverP;

  EoverP=Estar/kB/T;
  /* change the units from V/(m*Pa) to V/(cm*Torr) */
  EoverP=EoverP*133.0/100.0;
  alphaoverP=15.0*exp(-365.0/EoverP);
  /* change the units from 1/(cm*Torr) to 1/(m*Pa) */
  alphaoverP=alphaoverP/133.0*100.0;
  kovermueN_Air=alphaoverP*kB*T*Estar;
  return(kovermueN_Air);
}


/* townsend ionization coefficient for Air in cm3/s 
   valid range: 3.0e-20<E/N<249e-20 Vm2   
*/
double _kAir_123(double T, double Estar){
  double kAir,kAir_1,kAir_2,kAir_3;
  double fact1,fact2,fact3;

  kAir=0.0;  /* to avoid compiler warning */
    kAir_1=_kAir_1(T,Estar);
    kAir_2=_kAir_2(T,Estar);
    kAir_3=_kAir_3(T,Estar);
    if (Estar<13.7e-20) kAir=kAir_1;
    if (Estar>=13.7e-20 && Estar<30e-20) {
      fact2=(Estar-13.7e-20)/(30e-20-13.7e-20);
      fact1=1.0-fact2;
      kAir=fact1*kAir_1+fact2*kAir_2;
    }
    if (Estar>=30e-20 && Estar<54.8e-20) {
      fact3=(Estar-30e-20)/(54.8e-20-30e-20);
      fact2=1.0-fact3;
      kAir=fact2*kAir_2+fact3*kAir_3;
    }
    if (Estar>=54.8e-20) kAir=kAir_3;

  return(kAir);
}


/* townsend ionization coefficient for Air in Vm4 
   valid range: 3.0e-20<E/N<249e-20 Vm2   
*/
double _kovermueN_Air_123(double T, double Estar){
  double kovermueN_Air,kovermueN_Air_1,kovermueN_Air_2,kovermueN_Air_3;
  double fact1,fact2,fact3;

  kovermueN_Air=0.0; /* to avoid compiler warning */
    kovermueN_Air_1=_kovermueN_Air_1(T,Estar);
    kovermueN_Air_2=_kovermueN_Air_2(T,Estar);
    kovermueN_Air_3=_kovermueN_Air_3(T,Estar);
    if (Estar<13.7e-20) kovermueN_Air=kovermueN_Air_1;
    if (Estar>=13.7e-20 && Estar<30e-20) {
      fact2=(Estar-13.7e-20)/(30e-20-13.7e-20);
      fact1=1.0-fact2;
      kovermueN_Air=fact1*kovermueN_Air_1+fact2*kovermueN_Air_2;
    }
    if (Estar>=30e-20 && Estar<54.8e-20) {
      fact3=(Estar-30e-20)/(54.8e-20-30e-20);
      fact2=1.0-fact3;
      kovermueN_Air=fact2*kovermueN_Air_2+fact3*kovermueN_Air_3;
    }
    if (Estar>=54.8e-20) kovermueN_Air=kovermueN_Air_3;
  return(kovermueN_Air);
}




/* townsend ionization coefficient for N2 in cm3/s 
   valid range: 3.0e-20<E/N<187e-20 Vm2   
   Can also be used outside of this range (will give meaningful results)
*/
double _kN2_new(double T, double Estar){
  double kN2;
//log(kN2)=-0.0105809*pow(log(Estar),2.0) - 2.40411e-75*pow(log(Estar),46.0)

  kN2=T/300.0*exp(-0.0105809*pow(log(Estar),2.0) - 2.40411E-75*pow(log(Estar),46.0));
  return(kN2);
}


/* townsend ionization coefficient for Air in cm3/s 
   valid range: 3.0e-20<E/N<240e-20 Vm2   
   Can also be used outside of this range (will give meaningful results)
*/
double _kAir_new(double T, double Estar){
  double kAir;
  kAir=T/300.0*exp(-0.0105031*pow(log(Estar),2.0) - 2.40983E-75*pow(log(Estar),46.0));
  return(kAir);
}


/* townsend ionization coefficient for Air in Vm4
   valid range: 3.0e-20<E/N<240e-20 Vm2   
   Can also be used outside of this range (will give meaningful results)
*/
double _kovermueN_Air_new(double T, double Estar){
  double kovermueN_Air;
  kovermueN_Air=Estar*(T/300.0)*exp( - 14.803*pow(fabs(log(Estar)),0.2)-4.84797E-72*pow(log(Estar),44.0))*1e-6;
  return(kovermueN_Air);
}


int main( int argc, char **argv ){
  double Estar,T;

  T=300.0e0; /* cannot change this when using _kN2_1, _kN2_2, _kN2_3 */

  /* N2 */
  if (FALSE){
    for (Estar=3e-20; Estar<=187e-20; Estar*=1.07){
//      printf("\{%5.3f,%5.3f\}, ",log(Estar),log(_kN2_123(T,Estar)));
      printf("%E  %E  %E\n",Estar,_kN2_123(T,Estar),_kN2_new(T,Estar));
    }
  }

  /* Air */
  if (TRUE){  
    for (Estar=3e-20; Estar<=240e-20; Estar*=1.0813){
//    printf("\{%5.3f,%5.3f\}, ",log(Estar),log(_kAir_123(T,Estar)));
//      printf("%E  %E  %E\n",Estar,_kAir_123(T,Estar),_kAir_new(T,Estar));
//      printf("{%10.9f, %10.9f}, ",log(Estar),log(_kovermueN_Air_123(T,Estar)/Estar*1e6));
      printf("%E  %E  %E    %5.0f%% \n",Estar,_kovermueN_Air_123(T,Estar),_kovermueN_Air_new(T,Estar),(_kovermueN_Air_new(T,Estar)-_kovermueN_Air_123(T,Estar))/max(_kovermueN_Air_new(T,Estar),_kovermueN_Air_123(T,Estar))*100.0 );
    }
  }



  return(EXIT_SUCCESS);
}

