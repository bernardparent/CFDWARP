#include <cycle/tsemf/_tsemf.h>
#include <cycle/share/cycle_share.h>
#include <cycle/share/tsemf_share.h>
#include <src/bdry.h>

#ifdef EMFIELD



void find_tsemf_SOR_numsubiter_numcycle(long numsubitertot,long numsubiteropt, long *numsubiter, long *numcycle){
  *numsubiter=min(numsubiteropt,numsubitertot/2);
  if (*numsubiter==numsubiteropt) *numcycle=max(1,round((double)numsubitertot*0.5/(double)numsubiteropt)); else *numcycle=1;
  *numsubiter=numsubitertot/(2*(*numcycle));
}


void update_U_from_dUstar_emfield(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    for (flux=0; flux<nfe; flux++) np[l].bs->dUstaremfield[flux]*=gl->relaxEMF;   
    add_dUstar_to_U_emfield(np, gl, l, np[l].bs->dUstaremfield);
    thread_lock_global_set(gl,THREADTYPE_ALL);
    gl->effiter_U_emfield+=1.0e0/(double)gl->nn;
    thread_lock_global_unset(gl,THREADTYPE_ALL);
  }
}


void update_U_from_dUstar_emfield_without_relaxation(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    for (flux=0; flux<nfe; flux++) np[l].bs->dUstaremfield[flux]*=gl->relaxEMF;   
    add_dUstar_to_U_emfield(np, gl, l, np[l].bs->dUstaremfield);
    thread_lock_global_set(gl,THREADTYPE_ALL);
    gl->effiter_U_emfield+=1.0e0/(double)gl->nn;
    thread_lock_global_unset(gl,THREADTYPE_ALL);
  }
}


void init_dUstar_emfield_ADI(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  double dtau;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    assert_np(np[l],is_node_inner(np[l],TYPELEVEL_EMFIELD));
    for (flux=0; flux<nfe; flux++) {
      find_dtau_emfield(np,gl,l,flux,&dtau);  
      np[l].bs->dUstaremfield[flux]=-np[l].bs->Resemfield[flux]*dtau;
    }
  }
}




void update_dUstar_emfield_ADI(np_t *np, gl_t *gl, long theta, long ls, long le){
  long jj,l,flux,cnt;
  EXM_tdmaline_t *tdma;
  double dtau;

  cnt=0;
  for(l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) cnt++;
  tdma=(EXM_tdmaline_t *)malloc((cnt+5)*sizeof(EXM_tdmaline_t));

  for (flux=0; flux<nfe; flux++){
    jj=0;
    for(l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)) {
      jj++;

      find_dtau_emfield(np,gl,l,flux,&dtau);
      find_linearization_coefficients_inner_node_emfield(np, gl, l, theta, flux, &(tdma[jj].val[0]), &(tdma[jj].val[1]), &(tdma[jj].val[2]));

      tdma[jj].val[1]+=1.0/dtau;
      tdma[jj].val[3]=np[l].bs->dUstaremfield[flux]/dtau;

    }


    tdma[0].val[0]=0.0;
    find_linearization_coefficients_bdry_node_emfield(np, gl, _l_minus_one(ls,gl,theta), theta, +1,
                       flux, _node_type(np[_l_minus_one(ls,gl,theta)], TYPELEVEL_EMFIELD),
	  	     &(tdma[0].val[1]), &(tdma[0].val[2]), &(tdma[0].val[3]));

    tdma[jj+1].val[2]=0.0;
    find_linearization_coefficients_bdry_node_emfield(np, gl, _l_plus_one(le,gl,theta), theta, -1,
                     flux, _node_type(np[_l_plus_one(le,gl,theta)], TYPELEVEL_EMFIELD),
		     &(tdma[jj+1].val[1]), &(tdma[jj+1].val[0]), &(tdma[jj+1].val[3]));

    /* solve the TDMA */
    solve_TDMA_emfield(np, gl, theta, _l_minus_one(ls,gl,theta), _l_plus_one(le,gl,theta), TYPELEVEL_EMFIELD, tdma, jj+2);


    /*--------Here: add RHS of TDMA to Ustar. */
    jj=-1;
    for (l=_l_minus_one(ls,gl,theta); l!=_l_plus_one(_l_plus_one(le,gl,theta),gl,theta); l=_l_plus_one(l,gl,theta)){
      jj++;
      np[l].bs->dUstaremfield[flux]=tdma[jj].val[3]/tdma[jj].val[1];
    }

  }
  free(tdma);
}



void update_U_emfield_ADI(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np,gl,zone,&init_dUstar_emfield_ADI,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_dUstar_emfield_ADI,SWEEPTYPE_IJK, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_U_from_dUstar_emfield,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
}


#ifdef _TSEMF_STORE_COEFFICIENTS
static void init_dUstar_emfield_SOR(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    assert_np(np[l],is_node_valid(np[l],TYPELEVEL_EMFIELD));
    for (flux=0; flux<nfe; flux++) {
      np[l].bs->dUstaremfield[flux]=0.0;
    }
  }
}

#define SOR_SWEEP_FORWARD 1
#define SOR_SWEEP_BACKWARD 2


#ifdef DISTMPI


void update_dUstar_emfield_SOR_node(np_t *np, gl_t *gl, long l, long flux, int SOR_SWEEP){
  long dim,theta,thetasgn;
  double sum,RHS,Cp0,Cp1,dtau;
#ifndef NDEBUG
  long i,j,k;
#endif

  if (is_node_valid(np[l],TYPELEVEL_EMFIELD)) { 
	if (is_node_inner(np[l],TYPELEVEL_EMFIELD)) {
      // for inner node 
      sum=-np[l].bs->Resemfield[flux];
      for (dim=0; dim<nd; dim++){
        sum-=np[l].bs->coeffp1[dim][flux]*np[_al(gl,l,dim,+1)].bs->dUstaremfield[flux]
            +np[l].bs->coeffm1[dim][flux]*np[_al(gl,l,dim,-1)].bs->dUstaremfield[flux];
      }
#ifndef NDEBUG
      find_ijk_from_l(gl,l,&i,&j,&k);
      switch (SOR_SWEEP){
        case SOR_SWEEP_FORWARD:
          if (!np[_ai(gl,i-1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i-1,j,k);
          if (np[_ai(gl,i+1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i+1,j,k);

          if (!np[_ai(gl,i,j-1,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j-1,k);
          if (np[_ai(gl,i,j+1,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j+1,k);

#ifdef _3DL
          if (!np[_ai(gl,i,j,k-1)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j,k-1);
          if (np[_ai(gl,i,j,k+1)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j,k+1);
#endif
        break;
        case SOR_SWEEP_BACKWARD:
          if (!np[_ai(gl,i+1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i+1,j,k);
          if (np[_ai(gl,i-1,j,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i-1,j,k);

          if (!np[_ai(gl,i,j+1,k)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j+1,k);
          if (np[_ai(gl,i,j-1,k)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j-1,k);

#ifdef _3DL
          if (!np[_ai(gl,i,j,k+1)].bs->TSEMF_UPDATED) fatal_error("Node not updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j,k+1);
          if (np[_ai(gl,i,j,k-1)].bs->TSEMF_UPDATED) fatal_error("Node updated at i=%ld j=%ld k=%ld in update_dUstar_emfield_SOR_node()\n",i,j,k-1);
#endif
        break;
        default:
          fatal_error("SOR_SWEEP must be set to either SOR_SWEEP_BACKWARD or SOR_SWEEP_FORWARD");

      }
#endif
      dtau=np[l].bs->dtauemfield[flux];
      RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/(np[l].bs->coeffp0sum[flux]+1.0/dtau)*sum;
      np[l].bs->dUstaremfield[flux]=RHS; 
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
              if (np[_al(gl,l,theta,+1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node wrongly updated in update_dUstar_emfield_SOR_node()");
            } else {
              if (!np[_al(gl,l,theta,-1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node not updated in update_dUstar_emfield_SOR_node()");
            }
          break;
          case SOR_SWEEP_BACKWARD:
            if (thetasgn>0){
              if (!np[_al(gl,l,theta,+1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node not updated in update_dUstar_emfield_SOR_node()");
            } else {
              if (np[_al(gl,l,theta,-1)].bs->TSEMF_UPDATED) fatal_error("Near-bdry node wrongly updated in update_dUstar_emfield_SOR_node()");
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
  for1DL(i,is,ie)
    for2DL(j,js,je)
      for3DL(k,ks,ke)
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
          (*mpivars)=(double *)realloc(*mpivars,sizeof(double)*(*cntvars+1));
          (*mpivars)[*cntvars]=np[_ai(gl,i,j,k)].bs->dUstaremfield[flux];
          (*cntvars)++;
        }    
      end3DL
    end2DL
  end1DL
}


static void copy_mpivars_in_zone(np_t *np, gl_t *gl, long is, long js, long ks, long ie, long je, long ke, long flux, int numvars, double *mpivars){
  long i,j,k;
  int cntvars;
  cntvars=0;
  for1DL(i,is,ie)
    for2DL(j,js,je)
      for3DL(k,ks,ke)
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
          np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]=mpivars[cntvars];
          cntvars++;
#ifndef NDEBUG
          np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED=TRUE;
#endif
        }    
      end3DL
    end2DL
  end1DL
#ifndef NDEBUG
  if (cntvars!=numvars) printf("cntvars=%d  numvars=%d\n",cntvars,numvars);
  assert(cntvars==numvars);
#endif
}



static void update_dUstar_emfield_SOR_forward(np_t *np, gl_t *gl, long flux, zone_t zone, long numiter){
  long i,iter,cnt;
  long j,k,plane,planestart,planeend,jplusk;
#ifdef _3D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+(zone.ke-zone.ks+1)+1];
#endif
#ifdef _2D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+1];
#endif
  int numvars,rank,thisrank;
  double *mpivars;

  int packsize,buffersize,bbuffersize;
  double *buffer,*bbuffer;
  MPI_Status MPI_Status1;

  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &packsize );
  
  buffersize = (zone.ie-zone.is)*(zone.je-zone.js)if3DL(*(zone.ke-zone.ks)) * (MPI_BSEND_OVERHEAD + packsize);
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
        }
#ifdef _3DL
      }
#endif
    }
  }


  for (iter=0; iter<numiter; iter++)  {
    MPI_Buffer_attach( buffer, buffersize );
#ifndef NDEBUG
    for1DL(i,zone.is-1,zone.ie+1)
      for2DL(j,zone.js-1,zone.je+1)
        for3DL(k,zone.ks-1,zone.ke+1)
          np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED=FALSE;
        end3DL
      end2DL
    end1DL
#endif

    /* receive data from other processes before the threading starts */

    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is-1,zone.js,zone.ks,zone.is-1,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js-1,zone.ks,zone.ie,zone.js-1,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks-1,zone.ie,zone.je,zone.ks-1, flux, numvars,mpivars);
    }
#endif


    /* the threaded loop */
    for (plane=planestart; plane<=planeend; plane++){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt) schedule(static)  
#endif
      for (cnt=1; cnt<=lplane[plane][0]; cnt++){
        update_dUstar_emfield_SOR_node(np, gl, lplane[plane][cnt], flux, SOR_SWEEP_FORWARD);
      }
    }


    /* exchange data with other processes after the threading */

    find_mpivars_in_zone(np,gl,zone.ie,zone.js,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.is,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.je,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.js,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }
#ifdef _3DL
    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.je,zone.ks, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }
#endif


    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.ie+1,zone.js,zone.ks,zone.ie+1,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.je+1,zone.ks,zone.ie,zone.je+1,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke+1,zone.ie,zone.je,zone.ke+1, flux, numvars,mpivars);
    }
#endif


#ifndef NDEBUG
    for1DL(i,zone.is,zone.ie)
      for2DL(j,zone.js,zone.je)
        for3DL(k,zone.ks,zone.ke)
          if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD) && !np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED) 
            fatal_error("Node not updated correctly at i=%ld j=%ld k=%ld.",i,j,k);
        end3DL
      end2DL
    end1DL
#endif

    MPI_Buffer_detach( &bbuffer, &bbuffersize );
  }

  for (plane=planestart; plane<=planeend; plane++){
    free(lplane[plane]);
  }

  free(buffer);
  free(mpivars);

}




static void update_dUstar_emfield_SOR_backward(np_t *np, gl_t *gl, long flux, zone_t zone, long numiter){
  long i,iter,cnt;
  long j,k,plane,planestart,planeend,jplusk;
#ifdef _3D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+(zone.ke-zone.ks+1)+1];
#endif
#ifdef _2D
  long *lplane[((zone.ie-zone.is+1)+(zone.je-zone.js+1))+1];
#endif
  int numvars,rank,thisrank;
  double *mpivars;

  int packsize,buffersize,bbuffersize;
  double *buffer,*bbuffer;
  MPI_Status MPI_Status1;

  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &packsize );
  
  buffersize = (zone.ie-zone.is)*(zone.je-zone.js)if3DL(*(zone.ke-zone.ks)) * (MPI_BSEND_OVERHEAD + packsize);
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
        }
#ifdef _3DL
      }
#endif
    }
  }



  for (iter=0; iter<numiter; iter++)  {
    MPI_Buffer_attach( buffer, buffersize );
#ifndef NDEBUG
    for1DL(i,zone.is-1,zone.ie+1)
      for2DL(j,zone.js-1,zone.je+1)
        for3DL(k,zone.ks-1,zone.ke+1)
          np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED=FALSE;
        end3DL
      end2DL
    end1DL
#endif

    /* receive data from other processes before the threading starts */

    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.ie+1,zone.js,zone.ks,zone.ie+1,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.je+1,zone.ks,zone.ie,zone.je+1,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke+1,zone.ie,zone.je,zone.ke+1, flux, numvars,mpivars);
    }
#endif


    /* the threaded loop */
    for (plane=planeend; plane>=planestart; plane--){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt) schedule(static)  
#endif
      for (cnt=1; cnt<=lplane[plane][0]; cnt++){
        update_dUstar_emfield_SOR_node(np, gl, lplane[plane][cnt], flux, SOR_SWEEP_BACKWARD);
      }
    }


    /* exchange data with other processes after the threading */

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.is,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.ie,zone.js,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.ie+1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.js,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.je,zone.ks,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.je+1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }
#ifdef _3DL
    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks,zone.ie,zone.je,zone.ks, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }

    find_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ke,zone.ie,zone.je,zone.ke, flux, &numvars,&mpivars);
    rank=_node_rank(gl,zone.is,zone.js,zone.ke+1);
    if (rank!=thisrank){
      if (MPI_Bsend(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
      if (MPI_Bsend(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in update_dUstar_emfield_SOR");
    }
#endif


    rank=_node_rank(gl,zone.is-1,zone.js,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is-1,zone.js,zone.ks,zone.is-1,zone.je,zone.ke, flux, numvars,mpivars);
    }

    rank=_node_rank(gl,zone.is,zone.js-1,zone.ks);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js-1,zone.ks,zone.ie,zone.js-1,zone.ke, flux, numvars,mpivars);
    }

#ifdef _3DL
    rank=_node_rank(gl,zone.is,zone.js,zone.ks-1);
    if (rank!=thisrank){
      if (MPI_Recv(&numvars,1,MPI_INT,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      mpivars=(double *)realloc(mpivars,sizeof(double)*numvars);
      if (MPI_Recv(mpivars,numvars,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in update_dUstar_emfield_SOR");
      copy_mpivars_in_zone(np,gl,zone.is,zone.js,zone.ks-1,zone.ie,zone.je,zone.ks-1, flux, numvars,mpivars);
    }
#endif


#ifndef NDEBUG
    for1DL(i,zone.is,zone.ie)
      for2DL(j,zone.js,zone.je)
        for3DL(k,zone.ks,zone.ke)
          if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD) && !np[_ai(gl,i,j,k)].bs->TSEMF_UPDATED) 
            fatal_error("Node not updated correctly at i=%ld j=%ld k=%ld.",i,j,k);
        end3DL
      end2DL
    end1DL
#endif

    MPI_Buffer_detach( &bbuffer, &bbuffersize );
  }

  for (plane=planestart; plane<=planeend; plane++){
    free(lplane[plane]);
  }

  free(buffer);
  free(mpivars);

}


void update_dUstar_emfield_SOR(np_t *np, gl_t *gl, long flux, zone_t zone){
  long cycle,numcycle,numsubiter;
  find_tsemf_SOR_numsubiter_numcycle(gl->numsubiter_tsemf, TSEMF_SOR_NUMSUBITEROPT, &numsubiter, &numcycle);
  for (cycle=0; cycle<numcycle; cycle++){
    update_dUstar_emfield_SOR_forward(np, gl, flux, zone, numsubiter);
    update_dUstar_emfield_SOR_backward(np, gl, flux, zone, numsubiter);
  }
}


#else //not DISTMPI

void update_dUstar_emfield_SOR_istation(np_t *np, gl_t *gl, long flux, long i, zone_t zone, int SOR_SWEEP, long iter){
  long j,k,l,dim,theta,thetasgn;
  double sum,RHS,Cp0,Cp1,dtau;

  for2DL(j,zone.js,zone.je)
    for3DL(k,zone.ks,zone.ke)
      l=0;  // to avoid compiler warning
      switch (SOR_SWEEP) {
        case SOR_SWEEP_FORWARD:
          l=_ai(gl,i,j,k);
        break;
        case SOR_SWEEP_BACKWARD: 
          l=_ai(gl,i,zone.je-(j-zone.js),zone.ke-(k-zone.ks));
        break;
        default:
          fatal_error("SOR_SWEEP must be set to either SOR_SWEEP_FORWARD or SOR_SWEEP_BACKWARD in update_dUstar_emfield_SOR_istation.");
      }
      if (is_node_valid(np[l],TYPELEVEL_EMFIELD)) { 
	    if (is_node_inner(np[l],TYPELEVEL_EMFIELD)) {
          /* for inner node */
          sum=-np[l].bs->Resemfield[flux];
          for (dim=0; dim<nd; dim++){
            sum-=np[l].bs->coeffp1[dim][flux]*np[_al(gl,l,dim,+1)].bs->dUstaremfield[flux]
                +np[l].bs->coeffm1[dim][flux]*np[_al(gl,l,dim,-1)].bs->dUstaremfield[flux];
          }
          dtau=np[l].bs->dtauemfield[flux];
          RHS=(1.0-gl->relaxEMF)*np[l].bs->dUstaremfield[flux]+gl->relaxEMF/(np[l].bs->coeffp0sum[flux]+1.0/dtau)*sum;
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
    end3DL
  end2DL
}




void update_dUstar_emfield_SOR(np_t *np, gl_t *gl, long flux, zone_t zone){
  long i,cnt,ilocal,numsubiter,numcycle,cycle;

  find_tsemf_SOR_numsubiter_numcycle(gl->numsubiter_tsemf, TSEMF_SOR_NUMSUBITEROPT, &numsubiter, &numcycle);

  for (cycle=0; cycle<numcycle; cycle++){

    for (i=zone.is; i<=zone.ie+(numsubiter-1)*2; i++){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt,ilocal) schedule(dynamic)  
#endif
      for (cnt=0; cnt<numsubiter; cnt++)  {
        ilocal=i-cnt*2;
        if (ilocal>=zone.is && ilocal<=zone.ie) {
          update_dUstar_emfield_SOR_istation(np, gl, flux, ilocal, zone, SOR_SWEEP_FORWARD,cnt);
        }
      }
    }

    for (i=zone.ie; i>=zone.is-(numsubiter-1)*2; i--){
#ifdef OPENMPTHREADS
      #pragma omp parallel for private(cnt,ilocal) schedule(dynamic)  
#endif
      for (cnt=0; cnt<numsubiter; cnt++)  {
        ilocal=i+cnt*2;
        if (ilocal>=zone.is && ilocal<=zone.ie) {
          update_dUstar_emfield_SOR_istation(np, gl, flux, ilocal, zone, SOR_SWEEP_BACKWARD,cnt);
        }
      }
    }
  }
}



#endif //DISTMPI


void update_U_emfield_SOR(np_t *np, gl_t *gl, zone_t zone){
  long flux;
  sweep_with_1D_segments(np,gl,zone,&init_dUstar_emfield_SOR,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_valid,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  for (flux=0; flux<nfe; flux++){
    update_dUstar_emfield_SOR(np, gl, flux, zone);
  }
  sweep_with_1D_segments(np,gl,zone,&update_U_from_dUstar_emfield_without_relaxation,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);

}

#endif



/* exact inversion with the bandwidth corresponding to ie-is */
void update_dUstar_emfield_Newton_ij_1(np_t *np, gl_t *gl, long k, long flux, zone_t zone){
  long line,i,j,hbw,l;
  double C_ip1,C_im1,C_ip0,C_jp0,C_jp1,C_jm1,dtau;
  double *xdma;
  EXM_gl2D_t xdmagl;
  bool DIREC_FOUND;
  long theta,thetasgn;

  xdmagl.is=0;
  xdmagl.ie=(zone.ie-zone.is+1)*2+1+2;
  xdmagl.js=0;
  xdmagl.je=(zone.je-zone.js+1)
           *(zone.ie-zone.is+1)-1;
  hbw=(xdmagl.ie-xdmagl.is+1)/2-1;
  xdma=(double *)malloc((xdmagl.je-xdmagl.js+1)*(xdmagl.ie-xdmagl.is+1)*sizeof(double));

  /* set it up */
  /* init xdma */
  for (i=xdmagl.is; i<=xdmagl.ie; i++){
    for (j=xdmagl.js; j<=xdmagl.je; j++){
      xdma[EXM_ai2(xdmagl,i,j)]=0.0;
    }
  }
  line=0;
  for2DL(j,zone.js,zone.je)
    for1DL(i,zone.is,zone.ie)
      /* for inner node */
      l=_ai(gl,i,j,k);
      if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
          find_linearization_coefficients_inner_node_emfield(np, gl, l, 0, flux, &C_im1, &C_ip0, &C_ip1);
          find_linearization_coefficients_inner_node_emfield(np, gl, l, 1, flux, &C_jm1, &C_jp0, &C_jp1);

          find_dtau_emfield(np,gl,l,flux,&dtau); 
          /* i+0,j+0 */
          xdma[EXM_ai2(xdmagl,hbw,line)]=C_ip0+C_jp0+1.0/dtau;
          /* i-1,j+0 */ 
          xdma[EXM_ai2(xdmagl,hbw-1,line)]=C_im1;
          /* i+1,j+0 */ 
          xdma[EXM_ai2(xdmagl,hbw+1,line)]=C_ip1;
          /* i+0,j-1 */ 
          xdma[EXM_ai2(xdmagl,hbw-hbw+1,line)]=C_jm1;
          /* i+0,j+1 */ 
          xdma[EXM_ai2(xdmagl,hbw+hbw-1,line)]=C_jp1;
          xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]=np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]/dtau;           
        } else {
          xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
          DIREC_FOUND=find_bdry_direc(np, gl, l, TYPELEVEL_EMFIELD, &theta, &thetasgn);              
          if (DIREC_FOUND && theta==0) {
            find_linearization_coefficients_bdry_node_emfield(np, gl, l, theta, thetasgn, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+thetasgn,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
          } else {
            if (DIREC_FOUND && theta==1) {
              find_linearization_coefficients_bdry_node_emfield(np, gl, l, theta, thetasgn, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+thetasgn*(hbw-1),line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
            } else {
              if (is_node_inner(np[_ai(gl,i+1,j+1,k)], TYPELEVEL_EMFIELD)) {
                find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+hbw,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
              } else {
                if (is_node_inner(np[_ai(gl,i-1,j+1,k)], TYPELEVEL_EMFIELD)) {
                  find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+hbw-2,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
                } else {
                  if (is_node_inner(np[_ai(gl,i+1,j-1,k)], TYPELEVEL_EMFIELD)) {
                    find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw-hbw+2,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
                  } else {
                    if (is_node_inner(np[_ai(gl,i-1,j-1,k)], TYPELEVEL_EMFIELD)) {
                      find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw-hbw,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
                    } else {
                      //printf("Problem in update_dUstar_emfield_XDMA at node (%ld,%ld,%ld)\n",i,j,k);
                      xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
      } 
      line++;
    end1DL
  end2DL
  EXM_solve_XDMA(xdma, xdmagl);
  line=0;
  /* update dUstar */
  for2DL(j,zone.js,zone.je)
    for1DL(i,zone.is,zone.ie)
      if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
        np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]=xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]
                               /xdma[EXM_ai2(xdmagl,hbw,line)];
      }
      line++;
    end1DL
  end2DL

  free(xdma);
}


/* exact inversion with the bandwidth corresponding to je-js */
void update_dUstar_emfield_Newton_ij_2(np_t *np, gl_t *gl, long k, long flux, zone_t zone){
  long line,i,j,hbw,l;
  double *xdma;
  EXM_gl2D_t xdmagl;
  double C_ip1,C_im1,C_ip0,C_jp0,C_jp1,C_jm1,dtau;
  bool DIREC_FOUND;
  long theta,thetasgn;

  xdmagl.is=0;
  xdmagl.ie=(zone.je-zone.js+1)*2+1+2;
  xdmagl.js=0;
  xdmagl.je=(zone.je-zone.js+1)
           *(zone.ie-zone.is+1)-1;
  hbw=(xdmagl.ie-xdmagl.is+1)/2-1;
  xdma=(double *)malloc((xdmagl.je-xdmagl.js+1)*(xdmagl.ie-xdmagl.is+1)*sizeof(double));

  /* set it up */
  /* init xdma */
  for (i=xdmagl.is; i<=xdmagl.ie; i++){
    for (j=xdmagl.js; j<=xdmagl.je; j++){
      xdma[EXM_ai2(xdmagl,i,j)]=0.0;
    }
  }
  line=0;
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      l=_ai(gl,i,j,k);
      if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) { 
        /* for inner node */
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {

          find_linearization_coefficients_inner_node_emfield(np, gl, l, 0, flux, &C_im1, &C_ip0, &C_ip1);
          find_linearization_coefficients_inner_node_emfield(np, gl, l, 1, flux, &C_jm1, &C_jp0, &C_jp1);

          find_dtau_emfield(np,gl,l,flux,&dtau); 

          /* j+0,i+0 */
          xdma[EXM_ai2(xdmagl,hbw,line)]=(C_ip0)+(C_jp0)+1.0/dtau;    
          /* j-1,i+0 */ 
          xdma[EXM_ai2(xdmagl,hbw-1,line)]=C_jm1;
          /* j+1,i+0 */ 
          xdma[EXM_ai2(xdmagl,hbw+1,line)]=C_jp1;
          /* j+0,i-1 */ 
          xdma[EXM_ai2(xdmagl,hbw-hbw+1,line)]=C_im1;
          /* j+0,i+1 */ 
          xdma[EXM_ai2(xdmagl,hbw+hbw-1,line)]=C_ip1;
          xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]=np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]/dtau;
        } else {
          xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
          DIREC_FOUND=find_bdry_direc(np, gl, l, TYPELEVEL_EMFIELD, &theta, &thetasgn);

          if (DIREC_FOUND && theta==0) {
            find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), theta, thetasgn, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+thetasgn*(hbw-1),line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)])
            );
          } else {
            if (DIREC_FOUND && theta==1) {
              find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), theta, thetasgn, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+thetasgn,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)])
				   );
            } else {
              if (is_node_inner(np[_ai(gl,i+1,j+1,k)], TYPELEVEL_EMFIELD)) {
                find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+hbw,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)])
				);
              } else {
                if (is_node_inner(np[_ai(gl,i-1,j+1,k)], TYPELEVEL_EMFIELD)) {
                  find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw-hbw+2,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)])
				        );
                } else {
                  if (is_node_inner(np[_ai(gl,i+1,j-1,k)], TYPELEVEL_EMFIELD)) {
                    find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+hbw-2,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)])
				    );
                  } else {
                    if (is_node_inner(np[_ai(gl,i-1,j-1,k)], TYPELEVEL_EMFIELD)) {
                      find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 0, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw-hbw,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)])
                      );
                    } else {
                      //printf("Problem in EEF_UpdateQtildeXDMAj at node (%ld,%ld,%ld)\n",i,j,k);
                      xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
	  } 
      line++;
    end2DL
  end1DL
  EXM_solve_XDMA(xdma, xdmagl);
  line=0;
  /* update dU */
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
        np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]=xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]
                              /xdma[EXM_ai2(xdmagl,hbw,line)];
      }
      line++;
    end2DL
  end1DL

  free(xdma);
}


void update_dUstar_emfield_Newton_ij(np_t *np, gl_t *gl, long k, long flux, zone_t zone){
  /* first check if zone is valid */
  if (zone.js!=gl->domain_all.js || zone.je!=gl->domain_all.je
   || zone.is!=gl->domain_all.is || zone.ie!=gl->domain_all.ie) {
    fatal_error("The tsemf time stepping Newton_ij can not be used if the domain is split along i or j with MPI.");
  }
  if (zone.je-zone.js<zone.ie-zone.is) {
    update_dUstar_emfield_Newton_ij_2(np, gl, k, flux, zone);
  } else {
    update_dUstar_emfield_Newton_ij_1(np, gl, k, flux, zone);
  }
}


void update_dUstar_emfield_Newton_jk(np_t *np, gl_t *gl, long i, long flux, zone_t zone){
  long line,j,k,hbw,l;
  double C_kp1,C_km1,C_kp0,C_jp0,C_jp1,C_jm1,dtau;
  double *xdma;
  EXM_gl2D_t xdmagl;
  bool DIREC_FOUND;
  long theta,thetasgn;


  /* first check if zone is valid */
  if (zone.js!=gl->domain_all.js || zone.je!=gl->domain_all.je
   || zone.ks!=gl->domain_all.ks || zone.ke!=gl->domain_all.ke) {
    fatal_error("The tsemf time stepping Newton_jk can not be used if the domain is split along j or k with MPI.");
  }

  xdmagl.is=0;
  xdmagl.ie=(zone.ke-zone.ks+1)*2+1+2;
  xdmagl.js=0;
  xdmagl.je=(zone.je-zone.js+1)
           *(zone.ke-zone.ks+1)-1;
  hbw=(xdmagl.ie-xdmagl.is+1)/2-1;
  xdma=(double *)malloc((xdmagl.je-xdmagl.js+1)*(xdmagl.ie-xdmagl.is+1)*sizeof(double));

  /* set it up */
  /* init xdma */
  for (k=xdmagl.is; k<=xdmagl.ie; k++){
    for (j=xdmagl.js; j<=xdmagl.je; j++){
      xdma[EXM_ai2(xdmagl,k,j)]=0.0;
    }
  }
  line=0;
  for2DL(j,zone.js,zone.je)
    for3DL(k,zone.ks,zone.ke)
      /* for inner node */
      l=_ai(gl,i,j,k);
	  if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
          find_linearization_coefficients_inner_node_emfield(np, gl, l, 2, flux, &C_km1, &C_kp0, &C_kp1);
          find_linearization_coefficients_inner_node_emfield(np, gl, l, 1, flux, &C_jm1, &C_jp0, &C_jp1);
          find_dtau_emfield(np,gl,l,flux,&dtau); 
          /* k+0,j+0 */
          xdma[EXM_ai2(xdmagl,hbw,line)]=C_kp0+C_jp0+1.0/dtau;
          /* k-1,j+0 */ 
          xdma[EXM_ai2(xdmagl,hbw-1,line)]=C_km1;
          /* k+1,j+0 */ 
          xdma[EXM_ai2(xdmagl,hbw+1,line)]=C_kp1;
          /* k+0,j-1 */ 
          xdma[EXM_ai2(xdmagl,hbw-hbw+1,line)]=C_jm1;
          /* k+0,j+1 */ 
          xdma[EXM_ai2(xdmagl,hbw+hbw-1,line)]=C_jp1;
          assert(dtau!=0.0);
          xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]=np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]/dtau;            
	    } else {
          xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
          DIREC_FOUND=find_bdry_direc(np, gl, l, TYPELEVEL_EMFIELD, &theta, &thetasgn);            
          if (DIREC_FOUND && theta==2) {
            find_linearization_coefficients_bdry_node_emfield(np, gl, l, theta, thetasgn, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+thetasgn,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
          } else {
            if (DIREC_FOUND && theta==0) {
              xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
            } else {
              if (DIREC_FOUND && theta==1) {
                find_linearization_coefficients_bdry_node_emfield(np, gl, l, theta, thetasgn, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+thetasgn*(hbw-1),line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
              } else {
                if (is_node_inner(np[_ai(gl,i,j+1,k+1)], TYPELEVEL_EMFIELD)) {
                  find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 2, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+hbw,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
                } else {
                  if (is_node_inner(np[_ai(gl,i,j+1,k-1)], TYPELEVEL_EMFIELD)) {
                    find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 2, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw+hbw-2,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
                  } else {
                    if (is_node_inner(np[_ai(gl,i,j-1,k+1)], TYPELEVEL_EMFIELD)) {
                      find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 2, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw-hbw+2,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
                    } else {
                      if (is_node_inner(np[_ai(gl,i,j-1,k-1)], TYPELEVEL_EMFIELD)) {
                        find_linearization_coefficients_bdry_node_emfield(np, gl, _ai(gl,i,j,k), 2, +0, flux, 
                                   _node_type(np[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD), 
                                   &(xdma[EXM_ai2(xdmagl,hbw,line)]), 
                                   &(xdma[EXM_ai2(xdmagl,hbw-hbw,line)]),
                                   &(xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]));
                      } else {
                        //printf("Problem in update_dUstar_emfield_XDMA at node (%ld,%ld,%ld)\n",i,j,k);
                        xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
                      }
                    }
                  }
                }
              }
            }
          }
        }
	  } else {
	    xdma[EXM_ai2(xdmagl,hbw,line)]=1.0;
	  } 
      line++;
    end3DL
  end2DL
  EXM_solve_XDMA(xdma, xdmagl);
  line=0;
  /* update dUstar */
  for2DL(j,zone.js,zone.je)
    for3DL(k,zone.ks,zone.ke)
      if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
        assert(xdma[EXM_ai2(xdmagl,hbw,line)]!=0.0);
        np[_ai(gl,i,j,k)].bs->dUstaremfield[flux]=xdma[EXM_ai2(xdmagl,xdmagl.ie,line)]
                               /xdma[EXM_ai2(xdmagl,hbw,line)];
      }
      line++;
    end3DL
  end2DL
  free(xdma);
}


void update_U_emfield_ADIi(np_t *np, gl_t *gl, zone_t zone){
  long flux,i;
  sweep_with_1D_segments(np,gl,zone,&init_dUstar_emfield_ADI,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);

#ifndef _3D
  fatal_error("The tsemf pseudotime method TSEMF_ADIi can not be used in 2D or 1D.\n");
#endif

#ifdef OPENMPTHREADS
  #pragma omp parallel for private(i,flux) schedule(dynamic) 
#endif
  for (i=zone.is; i<=zone.ie; i++) {
    for (flux=0; flux<nfe; flux++) {
      update_dUstar_emfield_Newton_jk(np, gl, i, flux, zone);
    }
  }
  sweep_with_1D_segments(np,gl,zone,&update_dUstar_emfield_ADI,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_U_from_dUstar_emfield,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
}


void update_U_emfield_ADIk(np_t *np, gl_t *gl, zone_t zone){
  long flux;
  sweep_with_1D_segments(np,gl,zone,&init_dUstar_emfield_ADI,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);


#ifdef _3D
  long k;
#ifdef OPENMPTHREADS
  #pragma omp parallel for private(k,flux) schedule(dynamic) 
#endif
  for (k=zone.ks; k<=zone.ke; k++) {
    for (flux=0; flux<nfe; flux++) {
      update_dUstar_emfield_Newton_ij(np, gl, k, flux, zone);
    }
  }
#else
  for (flux=0; flux<nfe; flux++) {
    update_dUstar_emfield_Newton_ij(np, gl, 1, flux, zone);
  }
#endif
  sweep_with_1D_segments(np,gl,zone,&update_dUstar_emfield_ADI,SWEEPTYPE_K, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_U_from_dUstar_emfield,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
}


void init_dUstar_emfield_Newton(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  double dtau;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    assert_np(np[l],is_node_inner(np[l],TYPELEVEL_EMFIELD));
    for (flux=0; flux<nfe; flux++) {
      find_dtau_emfield(np,gl,l,flux,&dtau);  
      np[l].bs->dUstaremfield[flux]=-np[l].bs->Resemfield[flux]*dtau;
    }
  }
}


void update_U_emfield_Newton(np_t *np, gl_t *gl, zone_t zone){
  long flux;
#ifdef _3D
  fatal_error("Newton method can not be used in 3D.");
#endif
  sweep_with_1D_segments(np,gl,zone,&init_dUstar_emfield_Newton,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  for (flux=0; flux<nfe; flux++){
    update_dUstar_emfield_Newton_ij(np, gl, 1, flux, zone);
  }
  sweep_with_1D_segments(np,gl,zone,&update_U_from_dUstar_emfield,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);

}


#endif
