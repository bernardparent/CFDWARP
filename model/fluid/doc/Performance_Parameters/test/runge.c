

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define sqr(a) (a*a)

#define EXM_FORWARDEULER 0
#define EXM_MODIFIEDEULER 1
#define EXM_IMPROVEDEULER 2
#define EXM_RUNGEKUTTA 3


/* numerically differentiate FUNCT (which returns dx/dt)
   from t1 to t2, starting from x1, and returning x2 */
double EXP_NumDiff(double(*FUNCT)(void *, double, double), void *arg_ptr,
                 long METHOD, long n, double x1, double t1, double t2, long *error){
  double x,dt,t;
  double x_1,x_2,x_3;
  long cnt;

  *error=0;
  dt=(t2-t1)/(double)n;
  x=x1;
  t=t1;
  if (METHOD==EXM_FORWARDEULER)
    for (cnt=0; cnt<n; cnt++) {
      x+=dt*(*FUNCT)(arg_ptr,x,t);
      t+=dt;
    }
  if (METHOD==EXM_RUNGEKUTTA)
    for (cnt=0; cnt<n; cnt++) {
      x_1=x+0.5*dt*(*FUNCT)(arg_ptr,x,t);
      x_2=x+0.5*dt*(*FUNCT)(arg_ptr,x_1,t+0.5*dt);
      x_3=x+dt*(*FUNCT)(arg_ptr,x_2,t+0.5*dt);
      x+=(x_1-x)/3.0
        +(x_2-x)*2.0/3.0
        +(x_3-x)/3.0
        +dt/6.0*(*FUNCT)(arg_ptr,x_2,t+dt);
      t+=dt;
    }
  if (METHOD==EXM_IMPROVEDEULER)
    for (cnt=0; cnt<n; cnt++) {
      x+=0.5*dt*(*FUNCT)(arg_ptr,x,t)
        +0.5*dt*(*FUNCT)(arg_ptr,x+dt*(*FUNCT)(arg_ptr,x,t),t+dt);
      t+=dt;
    }
  if (METHOD==EXM_MODIFIEDEULER)
    for (cnt=0; cnt<n; cnt++) {
      x+=dt*(*FUNCT)(arg_ptr,x+0.5e0*dt*(*FUNCT)(arg_ptr,x,t),t+dt*0.5);
      t+=dt;
    }
  return(x);
}


typedef struct {
  double cp,R,h2,v2,k2;
} arg1_t;


/* gives dv/dPstar */
double funct1(void *arg1, double V, double Pstar){
  double tmp,T,v2,h2,R,cp,k2;
  v2=((arg1_t *)arg1)->v2;
  R=((arg1_t *)arg1)->R;
  cp=((arg1_t *)arg1)->cp;
  h2=((arg1_t *)arg1)->h2;
  k2=((arg1_t *)arg1)->k2;
  T=(sqr(v2)/2.0-sqr(V)/2.0+h2)/cp;
  tmp=-(k2*2.0/3.0+R*T)/(Pstar*V);
  return(tmp);
}


void main(){
  long error;
  double Pstar2,Pstar3,v3,v3exact,gamma;
  arg1_t arg1;
  arg1.cp=1000.0e0;
  arg1.R=286.0e0;
  arg1.h2=arg1.cp*450.0e0;
  arg1.v2=300.0e0;
  arg1.k2=1000.0e0;
  Pstar2=6000.0e0;
  Pstar3=100.0e0;
  gamma=arg1.cp/(arg1.cp-arg1.R);
  v3=DiffEquation(&funct1, &arg1, EXM_MODIFIEDEULER, 1000, arg1.v2, Pstar2, Pstar3, &error);

  v3exact=sqrt(2.0*arg1.h2+4.0/3.0*arg1.k2*arg1.cp/arg1.R+sqr(arg1.v2)-(2.0*arg1.h2+4.0/3.0*arg1.k2*arg1.cp/arg1.R)*pow(Pstar3/Pstar2,(gamma-1.0)/gamma));
  printf("v3=%E  v3exact=%E\n",v3,v3exact);
}
