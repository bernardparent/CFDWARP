#include <model/metrics/_metrics.h>



double _Omega(np_t np, gl_t *gl){
  double tmp;
  tmp=np.bs->Omega;
  return(tmp);
}


double _X(np_t np, long dim1, long dim2){
  double tmp;
  tmp=np.bs->X[dim1][dim2];
  return(tmp);
}


double _x(np_t np, long dim){
  double tmp;
  tmp=np.bs->x[dim];
  return(tmp);
}


void find_Omega_and_X_at_node(np_t *np, gl_t *gl, long l, double *Omega, dim2_t X){
#ifdef _2D
  long i,j;
  *Omega=0.25e0*(np[_al(gl,l,0,+1)].bs->x[0]-np[_al(gl,l,0,-1)].bs->x[0])
                   *(np[_al(gl,l,1,+1)].bs->x[1]-np[_al(gl,l,1,-1)].bs->x[1])
            -0.25e0*(np[_al(gl,l,1,+1)].bs->x[0]-np[_al(gl,l,1,-1)].bs->x[0])
                   *(np[_al(gl,l,0,+1)].bs->x[1]-np[_al(gl,l,0,-1)].bs->x[1]);
  for (i=0; i<nd; i++){
    for (j=0; j<nd; j++){
      X[i][j]=(krodelta(i,j)-0.5e0)/(*Omega)*(
          (+np[_al(gl,l,mod(i+1,nd),+1)].bs->x[mod(j+1,nd)]
           -np[_al(gl,l,mod(i+1,nd),-1)].bs->x[mod(j+1,nd)])
       );
     }
  }
#endif
#ifdef _3D
  long i,j;
  double tmp;
  tmp=0.0e0;
  for (i=0; i<nd; i++){
    tmp=tmp+0.125e0*(np[_al(gl,l,mod(i+0,nd),+1)].bs->x[0]-np[_al(gl,l,mod(i+0,nd),-1)].bs->x[0])
                   *(np[_al(gl,l,mod(i+1,nd),+1)].bs->x[1]-np[_al(gl,l,mod(i+1,nd),-1)].bs->x[1])
                   *(np[_al(gl,l,mod(i+2,nd),+1)].bs->x[2]-np[_al(gl,l,mod(i+2,nd),-1)].bs->x[2])
           -0.125e0*(np[_al(gl,l,mod(i+0,nd),+1)].bs->x[0]-np[_al(gl,l,mod(i+0,nd),-1)].bs->x[0])
                   *(np[_al(gl,l,mod(i+2,nd),+1)].bs->x[1]-np[_al(gl,l,mod(i+2,nd),-1)].bs->x[1])
                   *(np[_al(gl,l,mod(i+1,nd),+1)].bs->x[2]-np[_al(gl,l,mod(i+1,nd),-1)].bs->x[2]);
  }
  *Omega=tmp;
  for (i=0; i<nd; i++){
    for (j=0; j<nd; j++){
      X[i][j]=0.25e0/(*Omega)*(
          (+np[_al(gl,l,mod(i+1,nd),+1)].bs->x[mod(j+1,nd)]
           -np[_al(gl,l,mod(i+1,nd),-1)].bs->x[mod(j+1,nd)])
         *(+np[_al(gl,l,mod(i+2,nd),+1)].bs->x[mod(j+2,nd)]
           -np[_al(gl,l,mod(i+2,nd),-1)].bs->x[mod(j+2,nd)])
         -(+np[_al(gl,l,mod(i+1,nd),+1)].bs->x[mod(j+2,nd)]
           -np[_al(gl,l,mod(i+1,nd),-1)].bs->x[mod(j+2,nd)])
         *(+np[_al(gl,l,mod(i+2,nd),+1)].bs->x[mod(j+1,nd)]
           -np[_al(gl,l,mod(i+2,nd),-1)].bs->x[mod(j+1,nd)])
       );
    }
  }
#endif
}


void update_metrics_at_node(np_t *np, gl_t *gl, long l) {
  find_Omega_and_X_at_node(np, gl, l, &(np[l].bs->Omega), np[l].bs->X);
}


static double _dx_dX_int(np_t *np, gl_t *gl, long lL, long lR, long theta,
                       long dimx, long dimX){
  double tmp;
  if (dimX==theta) {
    tmp=np[lR].bs->x[dimx]-np[lL].bs->x[dimx];
  } else {
    tmp=0.25e0*(np[_al(gl,lL,dimX,+1)].bs->x[dimx]-np[_al(gl,lL,dimX,-1)].bs->x[dimx])
       +0.25e0*(np[_al(gl,lR,dimX,+1)].bs->x[dimx]-np[_al(gl,lR,dimX,-1)].bs->x[dimx]);
  }
  return(tmp);
}


void find_Omega_and_X_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta,
                          double *Omega, dim2_t X) {
#ifdef _2D
  long i,j;
  *Omega=_dx_dX_int(np,gl,lL,lR,theta,0,0)
               *_dx_dX_int(np,gl,lL,lR,theta,1,1)
               -_dx_dX_int(np,gl,lL,lR,theta,0,1)
               *_dx_dX_int(np,gl,lL,lR,theta,1,0);
  for (i=0; i<nd; i++){
    for (j=0; j<nd; j++){
      X[i][j]=(2.0e0*krodelta(i,j)-1.0e0)/notzero(*Omega,1e-99)*
         _dx_dX_int(np,gl,lL,lR,theta,mod(j+1,nd),mod(i+1,nd));
    }
  }
#endif
#ifdef _3D
  long i,j;

  *Omega=0.0e0;
  for (i=0; i<nd; i++){
    *Omega+=_dx_dX_int(np,gl,lL,lR,theta,0,i)
                  *_dx_dX_int(np,gl,lL,lR,theta,1,mod(i+1,nd))
                  *_dx_dX_int(np,gl,lL,lR,theta,2,mod(i+2,nd))
                  -_dx_dX_int(np,gl,lL,lR,theta,0,i)
                  *_dx_dX_int(np,gl,lL,lR,theta,1,mod(i+2,nd))
                  *_dx_dX_int(np,gl,lL,lR,theta,2,mod(i+1,nd));
  }
  for (i=0; i<nd; i++){
    for (j=0; j<nd; j++){
      X[i][j]=1.0e0/notzero(*Omega,1e-99)*(
         _dx_dX_int(np,gl,lL,lR,theta,mod(j+1,nd),mod(i+1,nd))*
                        _dx_dX_int(np,gl,lL,lR,theta,mod(j+2,nd),mod(i+2,nd))
        -_dx_dX_int(np,gl,lL,lR,theta,mod(j+2,nd),mod(i+1,nd))*
                        _dx_dX_int(np,gl,lL,lR,theta,mod(j+1,nd),mod(i+2,nd))
       );
    }
  }
#endif
}


void update_metrics_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta) {
  find_Omega_and_X_at_interface(np, gl, lL, lR, theta,
                          &(np[lL].bs->Omega_int[theta]), np[lL].bs->X_int[theta]);
}



/* lL is the node on the left of the interface
   lR is the node on the right of the interface
   theta is the dimension perpendicular to the interface
*/
void find_metrics_at_interface(np_t *np, gl_t *gl, long lL, long lR,
                               long theta, metrics_t *metrics){
  long dim,dim2,dim1;
  dim2_t X;
  double Omega;
   

  if (FALSE){
    /* for testing purposes only */
    find_Omega_and_X_at_interface(np, gl, lL, lR, theta, &Omega, X);
    metrics->Omega=Omega;
    for (dim=0; dim<nd; dim++) metrics->X[dim]=X[theta][dim];
    for (dim1=0; dim1<nd; dim1++) {
      for (dim2=0; dim2<nd; dim2++) {
        metrics->X2[dim1][dim2]=X[dim1][dim2];
      }
    }
  } else {
    metrics->Omega=np[lL].bs->Omega_int[theta];
    for (dim=0; dim<nd; dim++) metrics->X[dim]=np[lL].bs->X_int[theta][theta][dim];
    for (dim1=0; dim1<nd; dim1++) {
      for (dim2=0; dim2<nd; dim2++) {
        metrics->X2[dim1][dim2]=np[lL].bs->X_int[theta][dim1][dim2];
      }
    }
  }
}


/* lL is the node on the left of the interface
   lR is the node on the right of the interface
   theta is the dimension perpendicular to the interface
   theta2 is such that metrics->X[dim] corresponds to X[theta2][dim]
*/
void find_metrics_at_interface_2(np_t *np, gl_t *gl, long lL, long lR,
                               long theta, long theta2, metrics_t *metrics){
  long dim,dim2,dim1;
  dim2_t X;
  double Omega;
  find_Omega_and_X_at_interface(np, gl, lL, lR, theta, &Omega, X);
  metrics->Omega=Omega;
  for (dim=0; dim<nd; dim++) metrics->X[dim]=X[theta2][dim];
  for (dim1=0; dim1<nd; dim1++) {
    for (dim2=0; dim2<nd; dim2++) {
      metrics->X2[dim1][dim2]=X[dim1][dim2];
    }
  }
}


static double _fBCm1(double phi1, double phi2){
  double tmp;
  tmp=phi2+0.0e0*phi1;
  return(tmp);
}


static void update_bdry_metrics_transmissive(np_t npC, np_t npB, np_t npA){
  long dim1,dim2;
  npA.bs->Omega=_fBCm1(npC.bs->Omega,npB.bs->Omega);
  for (dim1=0; dim1<nd; dim1++){
    for (dim2=0; dim2<nd; dim2++){
      npA.bs->X[dim1][dim2]=_fBCm1(npC.bs->X[dim1][dim2],npB.bs->X[dim1][dim2]);
    }
  }
}


void update_metrics_at_bdry_node(np_t *np, gl_t *gl, int TYPELEVEL,
                       long l_A, long l_B, long l_C, int BDRYMETRICS){
  if (_node_type(np[l_A],TYPELEVEL)>=0) {
    switch (BDRYMETRICS) {
      case BDRYMETRICS_NORMAL:
        update_bdry_metrics_transmissive(np[l_C],np[l_B], np[l_A]);
      break;
      case BDRYMETRICS_CENTERED:
        update_metrics_at_node(np, gl, l_A);
      break;
    }
  } else {
    fatal_error("Called update_metrics_at_bdry_node for non-bdry node (?).");
  }
}


void find_metrics_at_node(np_t *np, gl_t *gl, long l,
                          long theta, metrics_t *metrics){
  long dim,dim2,dim1;
  metrics->Omega=_Omega(np[l],gl);
  for (dim=0; dim<nd; dim++) metrics->X[dim]=_X(np[l],theta,dim);
  for (dim1=0; dim1<nd; dim1++) {
    for (dim2=0; dim2<nd; dim2++) {
      metrics->X2[dim1][dim2]=_X(np[l],dim1,dim2);
    }
  }
}

