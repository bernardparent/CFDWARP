
#include <exm.h>
#include <extstdio.h>
#include <stdlib.h>

#define nf 3
#define Cp 1000.0e0
#define R 300.0e0

typedef double flux_t[nf];
typedef double sqmat_t[nf][nf];


void multiply_matrix_and_vector(sqmat_t sqmat, flux_t mat1, flux_t mat2){
  long row,col;
  double tmp;
  
  for (row=0; row<nf; row++){
     tmp=0.0e0;
     for (col=0; col<nf; col++){
     /* printf("row=%ld  col=%ld  sqmat=%E   mat1=%E\n",row,col,sqmat[row][col],mat1[col]); */
       tmp=tmp+sqmat[row][col]*mat1[col];     
     }
     
     mat2[row]=tmp;
  }
}

void find_dF_dU(double H, double u, double X, sqmat_t A){
  double PE,PD;
  
  PE=R/(Cp-R);
  PD=PE-2.0e0;

  A[0][0]=0.0e0;
  A[0][1]=X;  
  A[0][2]=0.0e0;

  A[1][0]=0.5e0*X*u*u*PD;
  A[1][1]=-X*u*PD;
  A[1][2]=X*PE;

  A[2][0]=X*u*(0.5e0*PE*u*u-H);
  A[2][1]=X*u*(H/u-PE*u);
  A[2][2]=X*u*(1.0e0+PE);

}

void find_F(double rhoOmega, double H, double u, double X, flux_t F){
  double POmega;
  POmega=rhoOmega*R*(H-u*u/2.0e0)/Cp;
  F[0]=rhoOmega*u*X;
  F[1]=rhoOmega*u*u*X+X*POmega;
  F[2]=rhoOmega*u*X*H;
}

void find_Q(double rhoOmega, double H, double u, flux_t U){
  
  U[0]=rhoOmega;
  U[1]=rhoOmega*u;
  U[2]=rhoOmega*((H-u*u*0.5e0)/Cp*(Cp-R)+u*u*0.5e0);  

}


void find_average(double rhoOmega_L, double H_L, double u_L, double X_L,
                  double rhoOmega_R, double H_R, double u_R, double X_R,
                  flux_t Q_L, flux_t Q_R, flux_t F_L, flux_t F_R,
                  double *H, double *u, double *X){
  double fact1L, fact1R, fact2R, fact2L, tmpL, tmpR, den;
  double PD,PE,a,b,c,c2,c1;
  flux_t dUstar,dF;
  long flux;

  tmpL=rhoOmega_L*u_L;
  tmpR=-rhoOmega_R*u_R;
  fact1L=tmpL/(tmpL+tmpR);
  fact1R=tmpR/(tmpL+tmpR);
  *X=fact1L*X_L+fact1R*X_R;


  tmpL=sqrt(rhoOmega_L);
  tmpR=sqrt(rhoOmega_R);
  fact2L=tmpL/(tmpL+tmpR);
  fact2R=tmpR/(tmpL+tmpR);
  *u=fact2L*u_L+fact2R*u_R;

  PE=R/(Cp-R);
  PD=PE-2.0e0;
  for (flux=0; flux<nf; flux++){
    dUstar[flux]=Q_R[flux]-Q_L[flux];
    dF[flux]=F_R[flux]-F_L[flux];
  }
  a=0.5e0*PD*dUstar[0];
  b=-PD*dUstar[1];
 
  c=PE*dUstar[2]-dF[1]/(*X);
  *u=(-b-sqrt(b*b-4.0e0*a*c))/(2.0e0*a);

  *H=(-0.5e0*(*u)*(*u)*(*u)*PE*dUstar[0]+(*u)*(*u)*PE*dUstar[1]
      -(*u)*(1.0e0+PE)*dUstar[2]+dF[2]/(*X))/(dUstar[1]-(*u)*dUstar[0]);
}

int main(){
  flux_t Q_R,Q_L,F_R,F_L,dUstar,dF,dF2;
  sqmat_t A;
  double rhoOmega_R,rhoOmega_L,H_R,H_L,u_R,u_L,X_R,X_L,H_A,u_A,X_A;
  long flux;

  seed_random();
  rhoOmega_R=random_double(1.0e0,10.0e0);
  rhoOmega_L=random_double(1.0e0,10.0e0);
  u_R=random_double(-1000.0e0,1000.0e0);
  u_L=random_double(-1000.0e0,1000.0e0);
  H_R=random_double(0.5e0*u_R*u_R+Cp*100.0e0,0.5e0*u_R*u_R+Cp*10000.0e0);
  H_L=random_double(0.5e0*u_L*u_L+Cp*100.0e0,0.5e0*u_L*u_L+Cp*10000.0e0);
  X_R=random_double(-1.0e0,1.0e0);
  X_L=random_double(-1.0e0,1.0e0);
/*  X_R=X_L; */

  find_Q(rhoOmega_R,H_R,u_R,Q_R);
  find_Q(rhoOmega_L,H_L,u_L,Q_L);
  find_F(rhoOmega_R,H_R,u_R,X_R,F_R);
  find_F(rhoOmega_L,H_L,u_L,X_L,F_L);
  find_average(rhoOmega_L, H_L, u_L, X_L, 
               rhoOmega_R, H_R, u_R, X_R, 
               Q_L,Q_R,F_L,F_R,
               &H_A, &u_A, &X_A);
  find_dF_dU(H_A,u_A,X_A,A);

  for (flux=0; flux<nf; flux++){
    dUstar[flux]=Q_R[flux]-Q_L[flux];
    dF[flux]=F_R[flux]-F_L[flux];
  }  

  multiply_matrix_and_vector(A, dUstar, dF2);  

  for (flux=0; flux<nf; flux++){
    printf("%E    %E \n",dF[flux],dF2[flux]);
  }  
  return(1);
}