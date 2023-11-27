// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1999-2004, 2019, 2021 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
   of conditions and the following disclaimer in the documentation and/or other
   materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <model/metrics/_metrics.h>
#include <lib/exm/exm.h>
#include <src/bdry.h>


double _Omega(np_t np, gl_t *gl){
  double tmp;
  tmp=np.bs->Omega;
  return(tmp);
}


double _X(np_t np, long dim1, long dim2){
  double tmp;
/*  assert(dim1>=0);
  assert(dim1<nd);
  assert(dim2>=0);
  assert(dim2<nd); */
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


#ifdef _3D
static void find_point_in_between_8_nodes(np_t *np, gl_t *gl, long l, long dim1, long dim1sgn, long dim2, long dim2sgn, long dim3, long dim3sgn, EXM_vec3D_t a){
  long dim; 

  for (dim=0; dim<nd; dim++){
    a[dim]=0.125*(
             _x(np[l],dim)
            +_x(np[_al(gl,l,dim1,dim1sgn)],dim)
            +_x(np[_al(gl,l,dim2,dim2sgn)],dim)
            +_x(np[_all(gl,l,dim1,dim1sgn,dim2,dim2sgn)],dim)
            +_x(np[_al(gl,l,dim3,dim3sgn)],dim)
            +_x(np[_all(gl,l,dim1,dim1sgn,dim3,dim3sgn)],dim)
            +_x(np[_all(gl,l,dim2,dim2sgn,dim3,dim3sgn)],dim)
            +_x(np[_al(gl,_all(gl,l,dim1,dim1sgn,dim2,dim2sgn),dim3,dim3sgn)],dim)
           );
  }
}
#endif


// these metrics at the interface preserve freestream
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
  EXM_vec3D_t a,b,c,d,ba,ca,da,abc,adc;

  *Omega=0.0e0;
  for (i=0; i<nd; i++){
    *Omega+=_dx_dX_int(np,gl,lL,lR,theta,0,i)
                  *_dx_dX_int(np,gl,lL,lR,theta,1,mod(i+1,nd))
                  *_dx_dX_int(np,gl,lL,lR,theta,2,mod(i+2,nd))
                  -_dx_dX_int(np,gl,lL,lR,theta,0,i)
                  *_dx_dX_int(np,gl,lL,lR,theta,1,mod(i+2,nd))
                  *_dx_dX_int(np,gl,lL,lR,theta,2,mod(i+1,nd));
  }
  find_point_in_between_8_nodes(np,gl,lL,theta,+1,mod(theta+1,nd),+1,mod(theta+2,nd),-1,a);
  find_point_in_between_8_nodes(np,gl,lL,theta,+1,mod(theta+1,nd),+1,mod(theta+2,nd),+1,b);
  find_point_in_between_8_nodes(np,gl,lL,theta,+1,mod(theta+1,nd),-1,mod(theta+2,nd),+1,c);
  find_point_in_between_8_nodes(np,gl,lL,theta,+1,mod(theta+1,nd),-1,mod(theta+2,nd),-1,d);

  for (i=0; i<nd; i++){
    ba[i]=b[i]-a[i];
    ca[i]=c[i]-a[i];
    da[i]=d[i]-a[i];
  }

  abc[0]=0.5*(ba[1]*ca[2]-ba[2]*ca[1]);
  abc[1]=0.5*(ba[2]*ca[0]-ba[0]*ca[2]);
  abc[2]=0.5*(ba[0]*ca[1]-ba[1]*ca[0]);

  adc[0]=-0.5*(da[1]*ca[2]-da[2]*ca[1]);
  adc[1]=-0.5*(da[2]*ca[0]-da[0]*ca[2]);
  adc[2]=-0.5*(da[0]*ca[1]-da[1]*ca[0]);

  for (i=0; i<nd; i++){
    if (i==theta){
      for (j=0; j<nd; j++){  
        X[i][j]=(abc[j]+adc[j])/notzero(*Omega,1e-99);
      }
//      if (i==theta) printf("%E  %E  %E\n",X[i][0],X[i][1],X[i][2]);
    } else {
      for (j=0; j<nd; j++){
        X[i][j]=1.0e0/notzero(*Omega,1e-99)*(
           _dx_dX_int(np,gl,lL,lR,theta,mod(j+1,nd),mod(i+1,nd))*
                        _dx_dX_int(np,gl,lL,lR,theta,mod(j+2,nd),mod(i+2,nd))
          -_dx_dX_int(np,gl,lL,lR,theta,mod(j+2,nd),mod(i+1,nd))*
                        _dx_dX_int(np,gl,lL,lR,theta,mod(j+1,nd),mod(i+2,nd))
        );
      }
    }
     // if (i==theta) printf("%E  %E  %E\n\n",X[i][0],X[i][1],X[i][2]);

  }
#endif
}


// these metrics at the interface do not preserve freestream in 3D
void find_Omega_and_X_at_interface_old(np_t *np, gl_t *gl, long lL, long lR, long theta,
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


// finds the unit normal vector perpendicular to the boundary surface and pointing towards the fluid
bool find_unit_vector_normal_to_boundary_plane(np_t *np, gl_t *gl, long lA, long lB, long lC, int TYPELEVEL, dim_t n){
  long cnt,nodefound,lA2,lA3,dim,iA,jA,kA,iB,jB,kB;
  EXM_vec3D_t vecA2,vecA3,vecB,vecnormal;
  double vecmag;
  nodefound=1;
  assert(is_node_bdry(np[lA], TYPELEVEL));
  
  // first check if lB is perpendicular to the surface along the generalized coordinates (can't deal with a corner node here)
  find_ijk_from_l(gl, lA, &iA, &jA, &kA);
  find_ijk_from_l(gl, lB, &iB, &jB, &kB);
  
  cnt=0;
  if (iA!=iB) cnt++;
  if (jA!=jB) cnt++;
#ifdef _3D
  if (kA!=kB) cnt++;
#endif  
  if (cnt!=1){
    //fatal_error("Inner node B is misaligned with boundary node A in find_unit_vector_normal_to_boundary_plane(): iA=%ld jA=%ld kA=%ld iB=%ld jB=%ld kB=%ld \n",iA,jA,kA,iB,jB,kB); 
    return(FALSE);
  }

  lA2=0; // to prevent compiler warning
  lA3=0; // to prevent compiler warning
  for (dim=0; dim<nd; dim++){
    if (is_node_bdry_with_single_direc(np,gl,_al(gl,lA,dim,+1), TYPELEVEL)){
      switch (nodefound){
        case 1: 
          lA2=_al(gl,lA,dim,+1); 
          nodefound++; 
        break;
        case 2: 
          lA3=_al(gl,lA,dim,+1); 
          nodefound++; 
        break;
      }
      
    } else {
      if (is_node_bdry_with_single_direc(np,gl,_al(gl,lA,dim,-1), TYPELEVEL)){
        switch (nodefound){
          case 1: 
            lA2=_al(gl,lA,dim,-1); 
            nodefound++; 
          break;
          case 2: 
            lA3=_al(gl,lA,dim,-1); 
            nodefound++; 
          break;
        }
        
      }
    }
  }
//  fprintf(stderr,"nodefound=%ld\n",nodefound);
  if (nodefound!=nd) return(FALSE);
  for (dim=0; dim<3; dim++) {
    if (dim<nd){
      vecA2[dim]=_x(np[lA2],dim)-_x(np[lA],dim);
      if (nd==3) vecA3[dim]=_x(np[lA3],dim)-_x(np[lA],dim); 
      vecB[dim]=_x(np[lB],dim)-_x(np[lA],dim); 
    } else { 
      vecA2[dim]=0.0;
      vecA3[dim]=0.0;
      vecB[dim]=0.0; 
    }
  }
  if (nd<3) {
     vecA3[0]=0.0;
     vecA3[1]=0.0;
     vecA3[2]=1.0;
  }
  EXM_cross_product(vecA2, vecA3, vecnormal);
  // now normalize the vector
  vecmag=EXM_vector_magnitude(vecnormal);
  for (dim=0; dim<3; dim++) vecnormal[dim]/=vecmag;
  // now make sure the vector points towards node lB
  if (EXM_dot_product(vecnormal,vecB)<0.0){
    //reverse the sign of vecnormal
    for (dim=0; dim<3; dim++) vecnormal[dim]=-vecnormal[dim];
  }
  // set n to vecnormal
  for (dim=0; dim<nd; dim++) n[dim]=vecnormal[dim];
  return(TRUE);
}


double _distance_between_near_bdry_node_and_boundary_plane(np_t *np, gl_t *gl, long lA, long lB, long lC, int TYPELEVEL){
  double dwall;
  dim_t n,vecBA;
  long i,j,k;
  long dim;
  if (!find_unit_vector_normal_to_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL, n)){
    find_ijk_from_l_all(gl, lA, &i, &j, &k);
    fatal_error("Couldn't find unit vector normal to boundary plane in _distance_between_near_bdry_node_and_boundary_plane at the node (%ld,%ld,%ld). This generally occurs when the boundary node is a corner node.",i,j,k);
    dwall=0.0;
  } else {
    for (dim=0; dim<nd; dim++){
      vecBA[dim]=_x(np[lB],dim)-_x(np[lA],dim); 
    }
    // multiply BA by n
    dwall=0.0;
    for (dim=0; dim<nd; dim++) dwall+=vecBA[dim]*n[dim];
    if (dwall<=0.0) {
      find_ijk_from_l_all(gl, lA, &i, &j, &k);
      fatal_error("Problem in _distance_between_near_bdry_node_and_boundary_plane at the node (%ld,%ld,%ld): dwall can not be equal to %E.",i,j,k,dwall);
    }
  }
  return(dwall);
}
