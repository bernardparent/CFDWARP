#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))

void find_Te_from_EoverN(double Estar, double *Te){
  double w[9];
  double sum,Estarexp,logEstar;
  long cnt;
  Estar=max(1e-60,Estar);
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



/* product of electron mobility and number density as a function of electron temperature */
double _mueN_from_Te(double Te){
  double mueN;
  mueN=3.74E19*exp(33.5/sqrt(log(Te)));
  return(mueN);
}


int main(){
  double EoverN,mueN,Te,Ve,Testar;
  EoverN=9e-20;
  find_Te_from_EoverN(EoverN, &Te);
  mueN=_mueN_from_Te(Te);
  Ve=mueN*EoverN;
  printf("EoverN=%E\nTe=%E\nmueN=%E\nVe=%E\nVe2=%E\n",EoverN,Te,mueN,Ve,pow(Te,1.0)*3.0);
  Testar=3000.0+9.0/1.5e7*Te*Te;
  printf("Testar=%E\n",Testar); 
}
