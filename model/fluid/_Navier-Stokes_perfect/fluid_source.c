#include "fluid.h"
#include "fluid_source.h"
#include <model/metrics/_metrics.h>
#include <model/_model.h>



#ifdef _2D
static void find_Saxi(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  double x1, V1;

  for (flux=0; flux<nf; flux++) S[flux]=0.0e0;

  if (gl->model.fluid.AXISYMMETRIC) {
    x1=np[l].bs->x[1];
    V1=_V(np[l],1);
    if (x1<1e-15) fatal_error("No node must lie on or below the y=0 axis when AXISYMMETRIC is set to TRUE.");
    for (flux=0; flux<nf; flux++){
      S[flux]=(-1.0/x1)*V1*np[l].bs->U[flux];
    }
  }
}
#else
static void find_Saxi(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  for (flux=0; flux<nf; flux++){
    S[flux]=0.0e0;
  }
}
#endif



void find_Sstar(np_t *np, gl_t *gl, long l, flux_t S){
  long flux;
  long dim;
  flux_t Saxi;

  for (flux=0; flux<nf; flux++){
    S[flux]=0.0;
  }
  S[nd+1]=np[l].bs->Qadd*_Omega(np[l],gl);
  for (dim=0; dim<nd; dim++) {
    S[1+dim]=np[l].bs->Fbody[dim]*_Omega(np[l],gl);
    S[1+nd]+=np[l].bs->Fbody[dim]*_V(np[l],dim)*_Omega(np[l],gl);
  }

  find_Saxi(np,gl,l,Saxi);
  for (flux=0; flux<nf; flux++) S[flux]+=_Omega(np[l],gl)*(Saxi[flux]);
}


void test_dSchem_dU(np_t *np, gl_t *gl, long l){
  wfprintf(stdout,"\n\nThe NavierStokesPerfect() module can not be used with chemical solvers.\n");
}


#ifdef _2D
static void find_dSaxi_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU){
  long row,col,flux;
  double x1,V1,rho;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dS_dU[row][col]=0.0e0;
    }
  }
  if (gl->model.fluid.AXISYMMETRIC) {
    x1=np[l].bs->x[1];
    V1=_V(np[l],1);
    rho=_rho(np[l]);
    assert(x1!=0.0);
    /* recall S[flux]=(-1.0/x1)*(rhov)/(rho)*U[flux] */
    /* wrt U */
    for (flux=0; flux<nf; flux++){
      dS_dU[flux][flux]+=(-1.0/x1)*V1;
    }
    /* wrt rho */
    for (row=0; row<nf; row++){
      dS_dU[row][0]+=(1.0/x1)*V1/rho*np[l].bs->U[row];
    }
    /* wrt rhoV1 */
    for (row=0; row<nf; row++){
      dS_dU[row][2]+=-(1.0/x1)/rho*np[l].bs->U[row];
    }

  }

}
#else
static void find_dSaxi_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dS_dU[row][col]=0.0e0;
    }
  }
}
#endif


void find_dSstar_dUstar(np_t *np, gl_t *gl, long l, sqmat_t dSstar_dUstar){
  long col,row;
  sqmat_t dSaxidU;
 
  find_dSaxi_dU(np,gl,l,dSaxidU);
  
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dSstar_dUstar[row][col]=dSaxidU[row][col];
    }
  }

}

