

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define sqr(a) (a*a)

/* n   is the number of steps
   x1  is lower limit
   x2  is upper limit
   error is 0 if no error is found
*/

double Integral(double(*FUNCT)(void *, double), void *arg_ptr,
                long n, double x1, double x2, long *error){
  double f,x,sum,dx;

  *error=0;
  dx=(x2-x1)/(double)n;
  sum=0.0e0;
  for (x=x1+dx*0.5e0; x<x2; x+=dx){
    f=(*FUNCT)(arg_ptr,x);
    sum+=f*dx;
  }
  return(sum);
}

double Integral2(double(*FUNCT)(void *, double), void *arg_ptr,
                 long n, double x1, double x2, long *error){
  double f1,f2,f3,x,sum,dx;

  *error=0;
  dx=(x2-x1)/(double)n;
  f1=(*FUNCT)(arg_ptr,x1);
  f2=(*FUNCT)(arg_ptr,x1+dx);
  f3=(*FUNCT)(arg_ptr,x1+2.0*dx);
  sum=dx*(27.0/24.0*f2+9.0/24.0*f1);
  for (x=x1+2.0*dx; x<x2-1.5*dx; x+=dx){
    f1=f2;
    f2=f3;
    f3=(*FUNCT)(arg_ptr,x+dx);
    sum+=dx*(f1/24.0+11.0/12.0*f2+f3/24.0);
  }
  f1=f2;
  f2=f3;
  f3=(*FUNCT)(arg_ptr,x2);
  sum+=dx*(27.0/24.0*f2+9.0/24.0*f3);
  return(sum);
}



double EXM_numerical_integration(double(*FUNCT)(void *, double), void *arg_ptr,
                    long METHOD, long n, double x1, double x2, long *error){
  double f1,f2,f3,x,sum,dx;

  *error=0;
  dx=(x2-x1)/(double)n;
  if (METHOD==EXM_NUMINTEG_RECTANGLES) {
    sum=0.0e0;
    for (x=x1+dx*0.5e0; x<x2; x+=dx){
      f=(*FUNCT)(arg_ptr,x);
      sum+=f*dx;
    }
  }
  if (METHOD==EXM_NUMINTEG_POLY2) { /* integrate using second degree polynomials */
    f1=(*FUNCT)(arg_ptr,x1);
    f2=(*FUNCT)(arg_ptr,x1+dx);
    f3=(*FUNCT)(arg_ptr,x1+2.0*dx);
    sum=dx*(27.0/24.0*f2+9.0/24.0*f1);
    for (x=x1+2.0*dx; x<x2-1.5*dx; x+=dx){
      f1=f2;
      f2=f3;
      f3=(*FUNCT)(arg_ptr,x+dx);
      sum+=dx*(f1/24.0+11.0/12.0*f2+f3/24.0);
    }
    f1=f2;
    f2=f3;
    f3=(*FUNCT)(arg_ptr,x2);
    sum+=dx*(27.0/24.0*f2+9.0/24.0*f3);
  }
  return(sum);
}





typedef struct {
  double cp,R,h1,V1,k1;
} arg1_t;

double funct1(void *arg1, double V){
  double tmp,T;
  T=(-V*V/2.0e0+((arg1_t *)arg1)->h1+sqr(((arg1_t *)arg1)->V1)/2.0e0)
    /((arg1_t *)arg1)->cp;
  tmp=V/(((arg1_t *)arg1)->R*T+2.0e0/3.0e0*((arg1_t *)arg1)->k1);
  return(tmp);
}



void main(){
  double R,cp,Cv,gamma;
  double dV,T1,P1,V1star,k1;
  double Pstag1,Pstag2;
  long error;
  arg1_t arg1;
  gamma=1.4e0;
  cp=1000.0e0;
  Cv=cp/gamma;
  R=cp-Cv;
  T1=300.0e0;
  P1=100000.0e0;
  V1star=3500.0e0;
  dV=2.0e-4;
  k1=1000.0e0;
  Pstag1=P1*pow(1.0e0+(gamma-1.0e0)/2.0e0*V1star*V1star/(2.0e0/3.0e0*gamma*k1+(gamma-1.0e0)*cp*T1),gamma/(gamma-1.0e0));
  arg1.cp=cp;  arg1.R=R;  arg1.V1=V1star;  arg1.h1=cp*T1; arg1.k1=k1;
  Pstag2=P1*exp(Integral2(&funct1, &arg1, 100, 0.0e0, V1star, &error));
  printf("Pstag1=%E  Pstag2=%E\n",Pstag1,Pstag2);

}
