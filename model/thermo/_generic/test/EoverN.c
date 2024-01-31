#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))
#define echarge 1.602176462e-19 /* C */
#define emass 9.10938188E-31  /* kg */
#define epsilon0 8.854187817E-12 /* permittivity of free space */
#define kB 1.3806488E-23   /* m2*kg/(s2*K) */


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



/* returns f at thisx given x,f,b at each data point */
double EXM_f_from_monotonespline(long N, double *x, double *f, double thisx){
  long i;
  double thisf,dx,t,deltam1,deltap0,deltap1,mp0,mp1,alpha,beta,tau;
  mp0=mp1=0.0; //to avoid compiler warning
  /* first find i that is such that x[i]<=thisx<=x[i+1] */
  i=-1;
  do {
    i++;
  } while(i<(N-1) && !(x[i]<=thisx && x[i+1]>=thisx));
  if (i>=N-1){
    fprintf(stderr,"Couldn't find an interval for x=%E in EXM_f_from_spline.",thisx);
  }
  assert(i<N-1);
  dx=x[i+1]-x[i];
  t=(thisx-x[i])/dx;
  deltap0=(f[i+1]-f[i])/dx;
  if(i>=1 && i<=N-3) {
    deltam1=(f[i]-f[i-1])/(x[i]-x[i-1]); 
    mp0=0.5*(deltam1+deltap0); 
    if(deltam1*deltap0<0.0) mp0=0.0; 
    deltap1=(f[i+2]-f[i+1])/(x[i+2]-x[i+1]); 
    mp1=0.5*(deltap0+deltap1); 
    if(deltap0*deltap1<0.0) mp1=0.0;
  }
  else if (i==0) {
    deltap1=(f[i+2]-f[i+1])/(x[i+2]-x[i+1]);
    mp0=deltap0;
    mp1=0.5*(deltap0+deltap1); 
    if(deltap0*deltap1<0.0) mp1=0.0; 
  }
  else if (i==N-2) {
    deltam1=(f[i]-f[i-1])/(x[i]-x[i-1]); 
    mp0=0.5*(deltam1+deltap0); 
    if(deltam1*deltap0<0.0) mp0=0.0; 
    mp1=deltap0;
  }
  else fprintf(stderr,"Input to EXM_f_from_monotonespline() out of range.");
  alpha=mp0/deltap0;
  beta=mp1/deltap0;
  if(alpha*alpha+beta*beta>9.0) {
    tau=3.0/sqrt(alpha*alpha+beta*beta);
    mp0=tau*alpha*deltap0;
    mp1=tau*beta*deltap0;
  }
  thisf=f[i]*(2.0*t*t*t-3.0*t*t+1)+dx*mp0*(t*t*t-2.0*t*t+t)+f[i+1]*(-2.0*t*t*t+3.0*t*t)+dx*mp1*(t*t*t-t*t);
  return(thisf);
}


static double _EoverN_from_Te_N2(double Te){
  double EoverN;
  /* Data in log-log coordinates. Obtained with BOLSIG+ using Morgan LXCat cross-sections */
  /* log K */
  double Te_data[] = 
  { 
    5.78359961896139,
    5.88538231327133,
    6.00674317027560,
    6.18906472706955,
    7.05656529477427,
    7.52656892402001,
    7.74971247533422,
    8.08618471195543,
    8.86485406595354,
    9.05804529498440,
    9.18479700062354,
    9.28657969493348,
    9.54147194456227,
    9.69562262438953,
    9.76461549587648,
    10.0010042739407,
    10.6516833690737,
    12.7963582069535,
    15.4249484703984
  };
  /* log Vm^2*/
  double EoverN_data[] = 
  { 
    -51.8608448501949,
    -51.3500192264290,
    -51.0135469898077,
    -50.6568720458690,
    -49.5582597572009,
    -49.0474341334349,
    -48.7109618968137,
    -48.3542869528750,
    -47.2556746642069,
    -46.7448490404409,
    -46.4083768038197,
    -46.0517018598809,
    -44.9530895712128,
    -44.4422639474468,
    -44.1057917108256,
    -43.7491167668869,
    -43.3907432663126,
    -40.6829260296843,
    -39.1439465808988
  };
  
  int N = sizeof(Te_data)/sizeof(Te_data[0]);
  Te = ( min( Te_data[N-1], max( log( Te ), Te_data[0] ) ) );
  EoverN = exp( EXM_f_from_monotonespline(N, Te_data, EoverN_data, Te) );
  return(EoverN);
}


/* electron energy loss function as a function of electron temperature */
double _zetae_from_Te(double Te){
  double k[7];
  double xi,Tepower[7];  
  long cnt;
  
  if (Te<19444.0){
    k[0]=+5.1572656E-4;
    k[1]=+3.4153708E-8;
    k[2]=-3.2100688E-11;
    k[3]=+1.0247332E-14;
    k[4]=-1.2153348E-18;
    k[5]=+7.2206246E-23;
    k[6]=-1.4498434E-27;
  } else {
    k[0]=+2.1476152E-1;
    k[1]=-4.4507259E-5;
    k[2]=+3.5155106E-9;
    k[3]=-1.3270119E-13;
    k[4]=+2.6544932E-18;
    k[5]=-2.7145800E-23;
    k[6]=+1.1197905E-28;
  }

  Te=max(300.0,Te);
  Te=min(60000.0,Te);
  Tepower[0]=1.0;
  Tepower[1]=Te;
  for (cnt=2; cnt<7; cnt++) Tepower[cnt]=Tepower[cnt-1]*Te;
  xi=0.0;
  for (cnt=0; cnt<7; cnt++) xi+=k[cnt]*Tepower[cnt];  
  return(xi);
}




void main(){
  double Te,mueN,Qen_over_NeN,Qe_over_NeN,EoverN;
  for (Te=300.0; Te<100000.0; Te*=2.0){
    mueN=_mueN_from_Te(Te);
    Qen_over_NeN=3.0*kB*Te*echarge*_zetae_from_Te(Te)/(2.0*emass*mueN);
    EoverN=_EoverN_from_Te_N2(Te);
    Qe_over_NeN=echarge*mueN*EoverN*EoverN;
    printf("Te=%EK  EoverN=%E Td  Qe_over_NeN=%E   Qen_over_NeN=%E\n",Te,EoverN/1e-21, Qe_over_NeN,Qen_over_NeN);
  } 


  for (Te=300.0; Te<100000.0; Te*=2.0){
    mueN=_mueN_from_Te(Te);
    EoverN=_EoverN_from_Te_N2(Te);
    printf("Te=%EK  EoverN=%E Td  zetae_new=%E   zetae_old=%E\n",Te,EoverN/1e-21, echarge*mueN*EoverN*EoverN/(3.0*kB*Te*echarge)*(2.0*emass*mueN), _zetae_from_Te(Te));
  } 

}
