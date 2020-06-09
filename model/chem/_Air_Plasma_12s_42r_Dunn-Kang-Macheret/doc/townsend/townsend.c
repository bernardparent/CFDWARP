#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sqr(a)   ((a)*(a))
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))


void find_Te_from_EoverN(double Estar, double *Te){
  double w[9];
  double sum,Estarexp,logEstar;
  long cnt;
  Estar=max(0.0,Estar);
  if (Estar<3e-19){ 
    w[0] = -3.69167532692495882511E+08;
    w[1] = -6.26956713747712671757E+07;
    w[2] = -4.65528490607805550098E+06;
    w[3] = -1.97394448288739687996E+05;
    w[4] = -5.22784662897089219769E+03;
    w[5] = -8.85545617874565635930E+01;
    w[6] = -9.36914737923363882821E-01;
    w[7] = -5.66073394421067171284E-03;
    w[8] = -1.49535882691330832494E-05;
    logEstar=log(Estar);
    Estarexp=logEstar;
    sum=w[0];
    for (cnt=1; cnt<9; cnt++) {
      sum+=w[cnt]*Estarexp;
      Estarexp*=logEstar;
    }
    *Te=exp(sum);
  } else {
    /* the following is a continuation of the curve assuming that the relationship between Estar and Te is linear for Estar>3e-19)*/
    *Te=59520.0+(Estar-3e-19)*1.5e23; 
  }
}


int main(){
  double kO2,kN2,kH2,EoverN,Te,Te_eV;

  EoverN=1e-21;
  do{
    find_Te_from_EoverN(EoverN,&Te);
    kO2=exp(-0.0102785*pow(log(EoverN),2.0)-2.42260E-75*pow(log(EoverN),46.0));
    kN2=exp(-0.0105809*pow(log(EoverN),2.0)-2.40411E-75*pow(log(EoverN),46.0));
    Te_eV=Te/1.16045E4;
    kH2=3.1228E-8*pow(Te_eV,0.17156)*exp(-20.07734/Te_eV);
//    printf("EoverN=%E  Te=%E  kN2=%E cm3/s  kO2=%E cm3/s  kH2=%E cm3/s\n",EoverN,Te,kN2,kO2,kH2);
    printf("%E  %E  %E  %E  %E  %E\n",EoverN,Te,kN2,kO2,kH2,kH2/kN2);
    EoverN*=1.2;
  } while(EoverN<187e-20);
//  } while(EoverN<3e-19);
  return(0);
}

