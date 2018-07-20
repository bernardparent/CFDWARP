#include <model/metrics/_metrics.h>


double _Jinv(np_t np){
  double tmp;
  tmp=np.mt->Jinv;
  return(tmp);
}

double _Omega(np_t np, gl_t *gl){
  double tmp;
  tmp=_Jinv(np);
#ifdef _1D
  tmp=_Jinv(np)*np.bs->Acs;
#endif
  return(tmp);
}

double _X(np_t np, long dim1, long dim2){
  double tmp;
  tmp=np.mt->X[dim1][dim2];
  return(tmp);
}


double _x(np_t np, long dim){
  double tmp;
  tmp=np.bs->x[dim];
  return(tmp);
}


static double dxdX_centered(np_t *np, gl_t *gl, long l, long dimx, long dimX){
  double tmp;
  tmp=0.5e0*(np[_al(gl,l,dimX,+1)].bs->x[dimx]-np[_al(gl,l,dimX,-1)].bs->x[dimx]);
  return(tmp);
}


/* second order accurate backward derivative or second order accurate central
   derivative when possible */
static double dxdX_2ndorder(np_t *np, gl_t *gl, long l, long dimx, long dimX){
  double tmp;
  if (!is_node_valid(np[_al(gl,l,dimX,+1)],TYPELEVEL_FLUID)) {
    assert(is_node_valid(np[_al(gl,l,dimX,-2)],TYPELEVEL_FLUID));
    assert(is_node_valid(np[_al(gl,l,dimX,-1)],TYPELEVEL_FLUID));
    assert(is_node_valid(np[_al(gl,l,dimX,+0)],TYPELEVEL_FLUID));
    tmp=1.5e0*np[_al(gl,l,dimX,+0)].bs->x[dimx]-2.0e0*np[_al(gl,l,dimX,-1)].bs->x[dimx]
       +0.5e0*np[_al(gl,l,dimX,-2)].bs->x[dimx];
  } else {
    if (!is_node_valid(np[_al(gl,l,dimX,-1)],TYPELEVEL_FLUID)) {
      assert(is_node_valid(np[_al(gl,l,dimX,+2)],TYPELEVEL_FLUID));
      assert(is_node_valid(np[_al(gl,l,dimX,+1)],TYPELEVEL_FLUID));
      assert(is_node_valid(np[_al(gl,l,dimX,+0)],TYPELEVEL_FLUID));
      tmp=-1.5e0*np[_al(gl,l,dimX,+0)].bs->x[dimx]+2.0e0*np[_al(gl,l,dimX,+1)].bs->x[dimx]
          -0.5e0*np[_al(gl,l,dimX,+2)].bs->x[dimx];
    } else {
      assert(is_node_valid(np[_al(gl,l,dimX,+1)],TYPELEVEL_FLUID));
      assert(is_node_valid(np[_al(gl,l,dimX,+0)],TYPELEVEL_FLUID));
      assert(is_node_valid(np[_al(gl,l,dimX,-1)],TYPELEVEL_FLUID));
      tmp=0.5e0*(np[_al(gl,l,dimX,+1)].bs->x[dimx]-np[_al(gl,l,dimX,-1)].bs->x[dimx]);
    }
  }
  return(tmp);
}


/* first order accurate backward derivative, or second order central when possible */
static double dxdX_1storder(np_t *np, gl_t *gl, long l, long dimx, long dimX){
  double tmp;
  if (!is_node_valid(np[_al(gl,l,dimX,+1)],TYPELEVEL_FLUID)) {
    assert(is_node_valid(np[_al(gl,l,dimX,-1)],TYPELEVEL_FLUID));
    assert(is_node_valid(np[_al(gl,l,dimX,+0)],TYPELEVEL_FLUID));
    tmp=1.0e0*np[_al(gl,l,dimX,+0)].bs->x[dimx]-1.0e0*np[_al(gl,l,dimX,-1)].bs->x[dimx];
  } else {
    if (!is_node_valid(np[_al(gl,l,dimX,-1)],TYPELEVEL_FLUID)) {
      assert(is_node_valid(np[_al(gl,l,dimX,+1)],TYPELEVEL_FLUID));
      assert(is_node_valid(np[_al(gl,l,dimX,+0)],TYPELEVEL_FLUID));
      tmp=-1.0e0*np[_al(gl,l,dimX,+0)].bs->x[dimx]+1.0e0*np[_al(gl,l,dimX,+1)].bs->x[dimx];
    } else {
      assert(is_node_valid(np[_al(gl,l,dimX,+1)],TYPELEVEL_FLUID));
      assert(is_node_valid(np[_al(gl,l,dimX,+0)],TYPELEVEL_FLUID));
      assert(is_node_valid(np[_al(gl,l,dimX,-1)],TYPELEVEL_FLUID));
      tmp=0.5e0*(np[_al(gl,l,dimX,+1)].bs->x[dimx]-np[_al(gl,l,dimX,-1)].bs->x[dimx]);
    }
  }
  return(tmp);
}


/* returns 0 in case of no error, and 1 in case of negative Jinv found */
static int FindMetricsLocal(np_t *np, gl_t *gl, long l,
                            double(dxdX(np_t *, gl_t *, long , long , long ))  ) {
  #ifdef _1D
    np[l].mt->Jinv=dxdX(np,gl,l,0,0);
    np[l].mt->X[0][0]=1.0e0/np[l].mt->Jinv;
  #endif
  #ifdef _2D
    long i,j;
    np[l].mt->Jinv=dxdX(np,gl,l,0,0)*dxdX(np,gl,l,1,1)
                  -dxdX(np,gl,l,0,1)*dxdX(np,gl,l,1,0);
    for (i=0; i<nd; i++){
      for (j=0; j<nd; j++){
        np[l].mt->X[i][j]=(2.0e0*krodelta(i,j)-1.0e0)/np[l].mt->Jinv*
           dxdX(np,gl,l,mod(j+1,nd),mod(i+1,nd));
      }
    }
  #endif
  #ifdef _3D
    long i,j;
    np[l].mt->Jinv=0.0e0;
    for (i=0; i<nd; i++){
      np[l].mt->Jinv+=dxdX(np,gl,l,0,i)*dxdX(np,gl,l,1,mod(i+1,nd))*dxdX(np,gl,l,2,mod(i+2,nd))
                     -dxdX(np,gl,l,0,i)*dxdX(np,gl,l,1,mod(i+2,nd))*dxdX(np,gl,l,2,mod(i+1,nd));
    }
    for (i=0; i<nd; i++){
      for (j=0; j<nd; j++){
        np[l].mt->X[i][j]=1.0e0/np[l].mt->Jinv*(
             dxdX(np,gl,l,mod(j+1,nd),mod(i+1,nd))*dxdX(np,gl,l,mod(j+2,nd),mod(i+2,nd))
            -dxdX(np,gl,l,mod(j+2,nd),mod(i+1,nd))*dxdX(np,gl,l,mod(j+1,nd),mod(i+2,nd))
         );
      }
    }
  #endif
  if (is_node_valid(np[l],TYPELEVEL_FLUID) && np[l].mt->Jinv<0.0e0) return(1); else return(0);

}

/* this subroutine is called for inner nodes */
void update_metrics_at_node(np_t *np, gl_t *gl, long l) {
  FindMetricsLocal(np,gl,l,&dxdX_centered);
}


/* subroutine called for boundary nodes */
void update_metrics_at_bdry_node(np_t *np, gl_t *gl, int TYPELEVEL,
                       long l_A, long l_B, long l_C, int BDRYMETRICS){
  if (_node_type(np[l_A],TYPELEVEL)>=0) {
    switch (BDRYMETRICS) {
      case BDRYMETRICS_NORMAL:
        if (FindMetricsLocal(np, gl, l_A, &dxdX_2ndorder)!=0){
          if (FindMetricsLocal(np, gl, l_A, &dxdX_1storder)!=0){
             fatal_error("Jinv is negative or node invalid at i=%ld"
               #ifdef _2DL
                 " j=%ld"
               #endif
               #ifdef _3DL
                 " k=%ld"
               #endif
                 ,_i(l_A,gl,0)
               #ifdef _2DL
                 ,_i(l_A,gl,1)
               #endif
               #ifdef _3DL
                 ,_i(l_A,gl,2)
               #endif
             );
          }
        }
      break;
      case BDRYMETRICS_CENTERED: /* for symmetry boundary, this one forces a centered discretization */
        FindMetricsLocal(np, gl, l_A, &dxdX_centered);
      break;
    }
  } else {
    fatal_error("Called update_metrics_at_bdry_node for a non-boundary node.");
  }
}


void find_Omega_and_X_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta,
                          double *Jinv, dim2_t X) {
  long i,j;
  *Jinv=0.0e0;
  for (i=0; i<nd; i++){
    *Jinv=0.5*(_Jinv(np[lL])+_Jinv(np[lR]));
  }
  assert_np(np[lL],*Jinv>0.0);
  for (i=0; i<nd; i++){
    for (j=0; j<nd; j++){
      X[i][j]=0.5*(_X(np[lL],i,j)+_X(np[lR],i,j));
    }
  }
}



void find_metrics_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, metrics_t *metrics){
  long dim,dim2;
  dim2_t X;
  double Jinv;
  find_Omega_and_X_at_interface(np, gl, lL, lR, theta, &Jinv, X);
#ifdef _1D
  metrics->Omega=0.5*(_Omega(np[lL],gl)+_Omega(np[lR],gl));
#else
  metrics->Omega=Jinv;
#endif
  for (dim=0; dim<nd; dim++) metrics->X[dim]=X[theta][dim];
  for (dim2=0; dim2<nd; dim2++){
    for (dim=0; dim<nd; dim++) metrics->X2[dim2][dim]=X[dim2][dim];
  }
}


