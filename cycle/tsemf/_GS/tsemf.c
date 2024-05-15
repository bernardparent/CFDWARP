// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015-2018 Bernard Parent

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

#include <cycle/tsemf/_tsemf.h>
#include <cycle/share/cycle_share.h>
#include <cycle/share/tsemf_share.h>
#include <src/bdry.h>

#ifndef EMFIELD
  #error the SOR2 tsemf.c must be compiled with EMFIELD
#endif
#define SOR_SWEEP_FORWARD 1
#define SOR_SWEEP_BACKWARD 2






static void add_to_tsemf_coefficients(np_t *np, long l, double coeff){
  long cnt;
  bool FOUND;
  FOUND=FALSE;
  cnt=0;
  do {
    if (np->bs->tsemfnode[cnt]==l){
      FOUND=TRUE;
      np->bs->tsemfcoeff2[cnt]+=coeff;
    }
    cnt++;
  } while(cnt<np->bs->tsemfnodenum);
  if (!FOUND){
    np->bs->tsemfcoeff2=(double *)realloc(np->bs->tsemfcoeff2,(1+np->bs->tsemfnodenum)*sizeof(double));
    np->bs->tsemfnode=(long *)realloc(np->bs->tsemfnode,(1+np->bs->tsemfnodenum)*sizeof(long));
    np->bs->tsemfnode[np->bs->tsemfnodenum]=l;
    np->bs->tsemfcoeff2[np->bs->tsemfnodenum]=coeff;
    np->bs->tsemfnodenum+=1;
  }
}


static void init_dUstar_emfield_SOR2(np_t *np, gl_t *gl, long flux, zone_t zone){
  long l,l2,dim,dim2;
  double fact;
  long i,j,k;

  for_ijk(zone,is,js,ks,ie,je,ke){
        l=_ai(gl,i,j,k);
        if (is_node_valid(np[l],TYPELEVEL_EMFIELD)){
          assert_np(np[l],is_node_valid(np[l],TYPELEVEL_EMFIELD));
          np[l].bs->dUstaremfield[flux]=0.0;
	      if (is_node_inner(np[l],TYPELEVEL_EMFIELD)) {
            np[l].bs->tsemfcoeff2=(double *)realloc(np[l].bs->tsemfcoeff2,sizeof(double));
            np[l].bs->tsemfnode=(long *)realloc(np[l].bs->tsemfnode,sizeof(long));
            np[l].bs->tsemfnode[0]=l;
            np[l].bs->tsemfcoeff2[0]=np[l].bs->coeffp0sum[flux];
            np[l].bs->tsemfnodenum=1;
            np[l].bs->tsemf_rhs=-np[l].bs->Resemfield[flux];

            for (dim=0; dim<nd; dim++) {
              /* do node on the right */
              l2=_al(gl,l,dim,+1);
              if (is_node_inner(np[l2],TYPELEVEL_EMFIELD)){
                fact=-np[l].bs->coeffp1[dim][flux]/np[l2].bs->coeffp0sum[flux];
                for (dim2=0; dim2<nd; dim2++){
                  add_to_tsemf_coefficients(&(np[l]),_al(gl,l2,dim2,-1),fact*np[l2].bs->coeffm1[dim2][flux]);
                  add_to_tsemf_coefficients(&(np[l]),_al(gl,l2,dim2,+1),fact*np[l2].bs->coeffp1[dim2][flux]);
                }
                np[l].bs->tsemf_rhs+=-fact*np[l2].bs->Resemfield[flux];
              } else {
                add_to_tsemf_coefficients(&(np[l]),_al(gl,l,dim,+1),np[l].bs->coeffp1[dim][flux]);
              }
              /* do node on the left */
              l2=_al(gl,l,dim,-1);
              if (is_node_inner(np[l2],TYPELEVEL_EMFIELD)){
                fact=-np[l].bs->coeffm1[dim][flux]/np[l2].bs->coeffp0sum[flux];
                for (dim2=0; dim2<nd; dim2++){
                  add_to_tsemf_coefficients(&(np[l]),_al(gl,l2,dim2,-1),fact*np[l2].bs->coeffm1[dim2][flux]);
                  add_to_tsemf_coefficients(&(np[l]),_al(gl,l2,dim2,+1),fact*np[l2].bs->coeffp1[dim2][flux]);
                }
                np[l].bs->tsemf_rhs+=-fact*np[l2].bs->Resemfield[flux];
              } else {
                add_to_tsemf_coefficients(&(np[l]),_al(gl,l,dim,-1),np[l].bs->coeffm1[dim][flux]);
              }
            }
          }
        }
  }
}




#ifdef DISTMPI


void update_dUstar_emfield_SOR2_node(np_t *np, gl_t *gl, long plane, long planetype, long l, long flux, int SOR_SWEEP){
  long dim,theta,thetasgn;
  double sum,RHS,Cp0,Cp1,dtau;
  long cnt,l2;

  if (is_node_valid(np[l],TYPELEVEL_EMFIELD)) { 
	if (is_node_inner(np[l],TYPELEVEL_EMFIELD)) {


      /* for inner node */
      if (planetype==0){
        sum=-np[l].bs->Resemfield[flux];
        for (dim=0; dim<nd; dim++){
          sum-=np[l].bs->coeffp1[dim][flux]*np[_al(gl,l,dim,+1)].bs->dUstaremfield[flux]
              +np[l].bs->coeffm1[dim][flux]*np[_al(gl,l,dim,-1)].bs->dUstaremfield[flux];
        }
        dtau=np[l].bs->dtauemfield[flux];
        RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/(np[l].bs->coeffp0sum[flux]+1.0/dtau)*sum;
      } else {
        sum=np[l].bs->tsemf_rhs;
        for (cnt=1; cnt<np[l].bs->tsemfnodenum; cnt++) {
          l2=np[l].bs->tsemfnode[cnt];
          assert(is_node_valid(np[l2],TYPELEVEL_EMFIELD));
          sum-=np[l].bs->tsemfcoeff2[cnt]*np[l2].bs->dUstaremfield[flux];
        }
        dtau=np[l].bs->dtauemfield[flux];
        RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/(np[l].bs->tsemfcoeff2[0]+1.0/dtau)*sum;

      }

      np[l].bs->dUstaremfield[flux]=RHS; 


#ifndef NDEBUG
      long i,j,k;
      find_ijk_from_l(gl,l,&i,&j,&k);
      switch (SOR_SWEEP){
        case SOR_SWEEP_FORWARD:
          if (!np[_ai(gl,i-1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i-1,j,k);
          if (np[_ai(gl,i+1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i+1,j,k);

          if (!np[_ai(gl,i,j-1,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j-1,k);
          if (np[_ai(gl,i,j+1,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j+1,k);

#ifdef _3DL
          if (!np[_ai(gl,i,j,k-1)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j,k-1);
          if (np[_ai(gl,i,j,k+1)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j,k+1);
#endif
        break;
        case SOR_SWEEP_BACKWARD:
          if (!np[_ai(gl,i+1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i+1,j,k);
          if (np[_ai(gl,i-1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i-1,j,k);

          if (!np[_ai(gl,i,j+1,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j+1,k);
          if (np[_ai(gl,i,j-1,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j-1,k);

#ifdef _3DL
          if (!np[_ai(gl,i,j,k+1)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j,k+1);
          if (np[_ai(gl,i,j,k-1)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR2_node()\n",i,j,k-1);
#endif
        break;
        default:
          fatal_error("SOR_SWEEP must be set to either SOR_SWEEP_BACKWARD or SOR_SWEEP_FORWARD");

      }
#endif
    } else {
      // for bdry node  
      if (find_bdry_direc(np, gl, l, TYPELEVEL_EMFIELD, &theta, &thetasgn)){
        find_linearization_coefficients_bdry_node_emfield(np, gl, l, theta, thetasgn, flux, 
                                     _node_type(np[l], TYPELEVEL_EMFIELD), &Cp0, &Cp1,&sum);
        sum-=Cp1*np[_al(gl,l,theta,thetasgn)].bs->dUstaremfield[flux];
        RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/Cp0*sum;
        np[l].bs->dUstaremfield[flux]=RHS;   

#ifndef NDEBUG
        switch (SOR_SWEEP){
          case SOR_SWEEP_FORWARD:
            if (thetasgn>0){
              if (np[_al(gl,l,theta,+1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node wrongly updated in update_dUstar_emfield_SOR2_node()");
            } else {
              if (!np[_al(gl,l,theta,-1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node not updated in update_dUstar_emfield_SOR2_node()");
            }
          break;
          case SOR_SWEEP_BACKWARD:
            if (thetasgn>0){
              if (!np[_al(gl,l,theta,+1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node not updated in update_dUstar_emfield_SOR2_node()");
            } else {
              if (np[_al(gl,l,theta,-1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node wrongly updated in update_dUstar_emfield_SOR2_node()");
            }
          break;
          default:
            fatal_error("SOR_SWEEP must be set to either SOR_SWEEP_BACKWARD or SOR_SWEEP_FORWARD");
        }
#endif

      }
    }
#ifndef NDEBUG
    np[l].bs->TSEMF_UPDATED=TRUE;
#endif
  }

}




static void find_mpivars_in_zone(np_t *np, gl_t *gl, long is, long js, long ks, long ie, long je, long ke, long flux, int *cntvars, double **mpivars){
  long i,j,k;
  *cntvars=0;
  for_1DL(i,is,ie){
    for_2DL(j,js,je){
      for_3DL(k,ks,ke){
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
          (*mpivars)=(double *)realloc(*mpivars,sizeof(double)*(*cntvars+1));
          (*mpivars)[*cntvars]=np[_ai(gl,i,j,k)].bs->dUstaremfield[flux];
          (*cntvars)++;
        }    
      }
    }
  }
}


static void copy_mpivars_in_zone(np_t *np, gl_t *gl, long is, long js, long ks, long ie, long je, long ke, long flux, int numvars, double *mpivars){
  long i,j,k;
  int cntvars;
  cntvars=0;
  for_1DL(i,is,ie){
    for_2DL(j,js,je){
      for_3DL(k,ks,ke){
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
          np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]=mpivars[cntvars];
          cntvars++;
#ifndef NDEBUG
          np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED=TRUE;
#endif
        }    
      }
    }
  }
#ifndef NDEBUG
  if (cntvars!=numvars) printf("cntvars=%d  numvars=%d\n",cntvars,numvars);
  assert(cntvars==numvars);
#endif
}



static void update_dUstar_emfield_SOR2_forward(np_t *np, gl_t *gl, long flux, zone_t zone, long numiter){
  long i,iter,cnt;
  long j,k,plane,planestart,planeend,jplusk;
#ifdef _3D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+(zone.ke-zone.ks+1)+1];
  long planetype[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+(zone.ke-zone.ks+1)+1];
#endif
#ifdef _2D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+1];
  long planetype[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+1];
#endif
  int numvars,rank,thisrank;
  double *mpivars;

  int packsize,buffersize,bbuffersize;
  double *buffer,*bbuffer;
  MPI_Status MPI_Status1;

  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &packsize );
  
  buffersize = min(INT_MAX,2*(zone.ie-zone.is)*(zone.je-zone.js)if3DL(*(zone.ke-zone.ks)) * (MPI_BSEND_OVERHEAD + packsize));
  buffer = (double *)malloc( buffersize );
  mpivars=(double *)malloc(sizeof(double));


  /* find lplane */
  planestart=1;
  planeend=((zone.ie-zone.is+1)+(zone.je-zone.js+1));
#ifdef _3DL
  planeend+=(zone.ke-zone.ks+1);
#endif
  for (plane=planestart; plane<=planeend; plane++){
    lplane[plane]=(long *)malloc(sizeof(long));
    cnt=0;
    lplane[plane][0]=0; 
    if (mod(plane,2)==0) planetype[plane]=0; else planetype[plane]=1;
    for (i=zone.is; i<=zone.ie; i++){
      jplusk=plane-(i-zone.is+1);
#ifdef _3DL
      for (k=zone.ks; k<=zone.ke; k++){
        j=jplusk-(k-zone.ks)+zone.js;
#else
        k=0;
        j=jplusk+zone.js;
#endif
        if (j>=zone.js && j<=zone.je){
          cnt++;
          lplane[plane]=(long *)realloc(lplane[plane],sizeof(long)*(cnt+1));
          lplane[plane][0]=cnt; 
          lplane[plane][cnt]=_ai(gl,i,j,k);
          if ((i==zone.is || i==zone.ie) && (j==zone.js || j==zone.je) if3D(&& (k==zone.ks || k==zone.ke))) planetype[plane]=0;

        }
#ifdef _3DL
      }
#endif
    }
  }


  for (iter=0; iter<numiter; iter++)  {
    MPI_Buffer_attach( buffer, buffersize );
#ifndef NDEBUG
    for_ijk(zone,is-1,js-1,ks-1,ie+1,je+1,ke+1){
          np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED=FALSE;
    }
#endif

    /* receive data from other processes before the threading starts */

    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is-2,zone.js,zone.ks,zone.is-1,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js-2,zone.ks,zone.ie,zone.js-1,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks-2,zone.ie,zone.je,zone.ks-1, flux, numvars,mpivars);
    }
#endif




    /* the threaded loop */
    for (plane=planestart; plane<=planeend; plane++){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt) schedule(static)  
#endif
      for (cnt=1; cnt<=lplane[plane][0]; cnt++){
        update_dUstar_emfield_SOR2_node(np, gl, plane, planetype[plane], lplane[plane][cnt], flux, SOR_SWEEP_FORWARD);
      }
    }


    /* exchange data with other processes after the threading */

    find_mpivars_in_zone(np,gl,zone.ie-1,zone.js,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.is+1,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.je-1,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.js+1,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }
#ifdef _3DL
    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke-1,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.je,zone.ks+1, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }
#endif




    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.ie+1,zone.js,zone.ks,zone.ie+2,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.je+1,zone.ks,zone.ie,zone.je+2,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26214,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke+1,zone.ie,zone.je,zone.ke+2, flux, numvars,mpivars);
    }
#endif




    /* the corner fluxes */
/*    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

#ifdef _3D
    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }


#endif


    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.js-1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.je+1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.js-1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.je+1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

#ifdef _3D

    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.js-1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.je+1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.js-1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.je+1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

#endif
*/

#ifndef NDEBUG
    for_ijk(zone,is,js,ks,ie,je,ke){
          if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD) && !np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED) 
            fatal_error("Node not updated correctly at i=%ld j=%ld k=%ld.",i,j,k);
    }
#endif

    MPI_Buffer_detach( &bbuffer, &bbuffersize );
  }

  for (plane=planestart; plane<=planeend; plane++){
    free(lplane[plane]);
  }
  free(buffer);
  free(mpivars);

}




static void update_dUstar_emfield_SOR2_backward(np_t *np, gl_t *gl, long flux, zone_t zone, long numiter){
  long i,iter,cnt;
  long j,k,plane,planestart,planeend,jplusk;
#ifdef _3D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+(zone.ke-zone.ks+1)+1];
  long planetype[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+(zone.ke-zone.ks+1)+1];
#endif
#ifdef _2D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+1];
  long planetype[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+1];
#endif
  int numvars,rank,thisrank;
  double *mpivars;

  int packsize,buffersize,bbuffersize;
  double *buffer,*bbuffer;
  MPI_Status MPI_Status1;

  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &packsize );
  
  buffersize = 2*(zone.ie-zone.is)*(zone.je-zone.js)if3DL(*(zone.ke-zone.ks)) * (MPI_BSEND_OVERHEAD + packsize);
  buffer = (double *)malloc( buffersize );
  mpivars=(double *)malloc(sizeof(double));


  planestart=1;
  planeend=((zone.ie-zone.is+1)+(zone.je-zone.js+1));
#ifdef _3DL
  planeend+=(zone.ke-zone.ks+1);
#endif
  for (plane=planeend; plane>=planestart; plane--){
    lplane[plane]=(long *)malloc(sizeof(long));
    cnt=0;
    lplane[plane][0]=0; 
    if (mod(plane,2)==0) planetype[plane]=0; else planetype[plane]=1;
    for (i=zone.ie; i>=zone.is; i--){
      jplusk=plane-(i-zone.is+1);
#ifdef _3DL
      for (k=zone.ke; k>=zone.ks; k--){
        j=jplusk-(k-zone.ks)+zone.js;
#else
        k=0;
        j=jplusk+zone.js;
#endif
        if (j>=zone.js && j<=zone.je){
          cnt++;
          lplane[plane]=(long *)realloc(lplane[plane],sizeof(long)*(cnt+1));
          lplane[plane][0]=cnt; 
          lplane[plane][cnt]=_ai(gl,i,j,k);
          if ((i==zone.is || i==zone.ie) && (j==zone.js || j==zone.je) if3D(&& (k==zone.ks || k==zone.ke))) planetype[plane]=0;
        }
#ifdef _3DL
      }
#endif
    }
  }



  for (iter=0; iter<numiter; iter++)  {
    MPI_Buffer_attach( buffer, buffersize );
#ifndef NDEBUG
    for_ijk(zone,is-1,js-1,ks-1,ie+1,je+1,ke+1){
          np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED=FALSE;
    }
#endif

    /* receive data from other processes before the threading starts */

    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.ie+1,zone.js,zone.ks,zone.ie+2,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.je+1,zone.ks,zone.ie,zone.je+2,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke+1,zone.ie,zone.je,zone.ke+2, flux, numvars,mpivars);
    }
#endif


    /* the threaded loop */
    for (plane=planeend; plane>=planestart; plane--){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt) schedule(static)  
#endif
      for (cnt=1; cnt<=lplane[plane][0]; cnt++){
        update_dUstar_emfield_SOR2_node(np, gl, plane, planetype[plane], lplane[plane][cnt], flux, SOR_SWEEP_BACKWARD);
      }
    }


    /* exchange data with other processes after the threading */

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.is+1,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.ie-1,zone.js,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.js+1,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.je-1,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }
#ifdef _3DL
    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.je,zone.ks+1, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke-1,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }
#endif


    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is-2,zone.js,zone.ks,zone.is-1,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js-2,zone.ks,zone.ie,zone.js-1,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,26234,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks-2,zone.ie,zone.je,zone.ks-1, flux, numvars,mpivars);
    }
#endif



    /* the corner fluxes */

/*
    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

#ifdef _3D
    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.is,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.js,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&(np[_ai(gl,zone.ie,zone.je,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR2");
    }


#endif


    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.js-1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.je+1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.js-1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.je+1,zone.ks-1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

#ifdef _3D

    rank=_node_rank(gl,zone.is-1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.js-1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.is-1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.is-1,zone.je+1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.js-1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.js-1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

    rank=_node_rank(gl,zone.ie+1,zone.je+1,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&(np[_ai(gl,zone.ie+1,zone.je+1,zone.ke+1)].bs->dUstaremfield[flux]),1,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR2");
    }

#endif
*/


#ifndef NDEBUG
    for_ijk(zone,is,js,ks,ie,je,ke){
          if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD) && !np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED) 
            fatal_error("Node not updated correctly at i=%ld j=%ld k=%ld.",i,j,k);
    }
#endif

    MPI_Buffer_detach( &bbuffer, &bbuffersize );
  }

  for (plane=planestart; plane<=planeend; plane++){
    free(lplane[plane]);
  }

  free(buffer);
  free(mpivars);

}



void update_dUstar_emfield_SOR2(np_t *np, gl_t *gl, long flux, zone_t zone){
  long cycle,numcycle,numsubiter;
  find_tsemf_SOR_numsubiter_numcycle(gl->numsubiter_tsemf, TSEMF_SOR_NUMSUBITEROPT, &numsubiter, &numcycle);
  for (cycle=0; cycle<numcycle; cycle++){
    update_dUstar_emfield_SOR2_forward(np, gl, flux, zone, numsubiter);
    update_dUstar_emfield_SOR2_backward(np, gl, flux, zone, numsubiter);
  }
}


#else //not DISTMPI



static void update_dUstar_emfield_SOR2_istation(np_t *np, gl_t *gl, long flux, long i, zone_t zone, int SOR_SWEEP){
  long j,k,l,l2,dim,theta,thetasgn,cnt;
  double sum,RHS,Cp0,Cp1,dtau;

  for_2DL(j,zone.js,zone.je){
    for_3DL(k,zone.ks,zone.ke){
      switch (SOR_SWEEP) {
        case SOR_SWEEP_FORWARD:
          l=_ai(gl,i,j,k);
        break;
        case SOR_SWEEP_BACKWARD: 
          l=_ai(gl,i,zone.je-(j-zone.js),zone.ke-(k-zone.ks));
        break;
      }
      if (is_node_valid(np[l],TYPELEVEL_EMFIELD)) { 
	    if (is_node_inner(np[l],TYPELEVEL_EMFIELD)) {
          /* for inner node */
          if (mod(i+j+k,2)==0){
            sum=-np[l].bs->Resemfield[flux];
            for (dim=0; dim<nd; dim++){
              sum-=np[l].bs->coeffp1[dim][flux]*np[_al(gl,l,dim,+1)].bs->dUstaremfield[flux]
                  +np[l].bs->coeffm1[dim][flux]*np[_al(gl,l,dim,-1)].bs->dUstaremfield[flux];
            }
            dtau=np[l].bs->dtauemfield[flux];
            RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/(np[l].bs->coeffp0sum[flux]+1.0/dtau)*sum;
          } else {
            sum=np[l].bs->tsemf_rhs;
            for (cnt=1; cnt<np[l].bs->tsemfnodenum; cnt++) {
              l2=np[l].bs->tsemfnode[cnt];
              assert(is_node_valid(np[l2],TYPELEVEL_EMFIELD));
              sum-=np[l].bs->tsemfcoeff2[cnt]*np[l2].bs->dUstaremfield[flux];
            }
            dtau=np[l].bs->dtauemfield[flux];
            RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/(np[l].bs->tsemfcoeff2[0]+1.0/dtau)*sum;

          }
          np[l].bs->dUstaremfield[flux]=RHS;
        } else {
          /* for bdry node */ 
          if (find_bdry_direc(np, gl, l, TYPELEVEL_EMFIELD, &theta, &thetasgn)){
            find_linearization_coefficients_bdry_node_emfield(np, gl, l, theta, thetasgn, flux, 
                                     _node_type(np[l], TYPELEVEL_EMFIELD), &Cp0, &Cp1,&sum);
            sum-=Cp1*np[_al(gl,l,theta,thetasgn)].bs->dUstaremfield[flux];
            RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/Cp0*sum;
            np[l].bs->dUstaremfield[flux]=RHS;   
          }
        }
      }
    }
  }
}



void update_dUstar_emfield_SOR2(np_t *np, gl_t *gl, long flux, zone_t zone){
  long i,cnt,ilocal,numsubiter,numcycle,cycle;

  find_tsemf_SOR_numsubiter_numcycle(gl->numsubiter_tsemf, TSEMF_SOR_NUMSUBITEROPT, &numsubiter, &numcycle);

  for (cycle=0; cycle<numcycle; cycle++){

    for (i=zone.is; i<=zone.ie+(numsubiter-1)*3; i++){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt,ilocal) schedule(dynamic)  
#endif
      for (cnt=0; cnt<numsubiter; cnt++)  {
        ilocal=i-cnt*3;
        if (ilocal>=zone.is && ilocal<=zone.ie) {
          update_dUstar_emfield_SOR2_istation(np, gl, flux, ilocal, zone, SOR_SWEEP_FORWARD);
        }
      }
    }

    for (i=zone.ie; i>=zone.is-(numsubiter-1)*3; i--){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt,ilocal) schedule(dynamic)  
#endif
      for (cnt=0; cnt<numsubiter; cnt++)  {
        ilocal=i+cnt*3;
        if (ilocal>=zone.is && ilocal<=zone.ie) {
          update_dUstar_emfield_SOR2_istation(np, gl, flux, ilocal, zone, SOR_SWEEP_BACKWARD);
        }
      }
    }
  }
}

#endif //DISTMPI

void update_U_emfield_SOR2(np_t *np, gl_t *gl, zone_t zone){
  long flux;
  for (flux=0; flux<nfe; flux++){
    init_dUstar_emfield_SOR2(np, gl, flux, zone);
    update_dUstar_emfield_SOR2(np, gl, flux, zone);
  }
  sweep_with_1D_segments(np,gl,zone,&update_U_from_dUstar_emfield_without_relaxation,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
}




void update_U_emfield(np_t *np, gl_t *gl, zone_t zone){
 
  switch (gl->tsemfmethod) {
    case TSEMF_DEFAULT:  /* default when no scheme is specified within control file */
      update_U_emfield_SOR2(np, gl, zone);
    break;
    case TSEMF_NEWTON: 
      update_U_emfield_Newton(np, gl, zone);
    break;
    case TSEMF_ADI: 
      update_U_emfield_ADI(np, gl, zone);
    break;
    case TSEMF_SOR: 
      update_U_emfield_SOR(np, gl, zone);
    break;
    case TSEMF_SOR2: 
      update_U_emfield_SOR2(np, gl, zone);
    break;
    default:
      fatal_error("The EMF pseudotime stepping method can not be set to %ld. It must be set to TSEMF_ADI, TSEMF_SOR, or TSEMF_SOR2.",gl->tsemfmethod);
  }
}

