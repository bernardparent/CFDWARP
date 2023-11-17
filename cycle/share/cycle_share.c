// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2018,2021 Bernard Parent
Copyright 2020 Minindu Weerakoon
Copyright 2001 Giovanni Fusina
Copyright 2002 Thomas E. Schwartzentruber
Copyright 2021 Prasanna Thoguluva Rajendran

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

#include <cycle/share/cycle_share.h>
#include <src/data.h>
#include <src/common.h>
#include <src/bdry.h>
#include <src/init.h>
#include <src/post.h>
#include <cycle/ts/_ts.h>
#include <cycle/tsemf/_tsemf.h>
#include <cycle/_cycle.h>
#include <cycle/res/_res.h>
#include <cycle/resconv/_resconv.h>
#include <cycle/restime/_restime.h>
#include <model/fluid/_fluid.h>
#include <model/emfield/_emfield.h>
#include <model/metrics/_metrics.h>
#include <model/fluid/_fluid.h>
#include <model/chem/_chem.h>


#ifdef OPENMPTHREADS
  #define  maxloopthread LONG_MAX
  #define  maxzonethread LONG_MAX
#else
  #define  maxloopthread 256
  #define  maxzonethread 256
#endif

#define MAXRATIO_DTAUMAX_DTAUMIN 100.0

typedef struct {
  np_t *np;
  gl_t *gl;
  long theta,ls,le;
} segment_t;


typedef struct {
  np_t *np;
  gl_t *gl;
  long theta,ls,le;
  void (*funct)(np_t *, gl_t *, long, long, long);
} segmentarg_t;


typedef struct {
  np_t *np;
  gl_t *gl;
  zone_t zone;
  void (*funct)(np_t *, gl_t *, zone_t);
} threadzone_t;


void *segmentfunct(void *segmentarg){
  (((segmentarg_t *) segmentarg)->funct)(
     ((segmentarg_t *) segmentarg)->np,
     ((segmentarg_t *) segmentarg)->gl,
     ((segmentarg_t *) segmentarg)->theta,
     ((segmentarg_t *) segmentarg)->ls,
     ((segmentarg_t *) segmentarg)->le);
  return(NULL);
}


void find_musclvarscycle(np_t np, gl_t *gl, musclvarscycle_t musclvars){
  find_musclvars(np,gl,musclvars);
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_MUSCLVARS
  long flux;
//  for (flux=0; flux<nf; flux++) musclvars[nf+flux]=musclvars[flux];
  for (flux=0; flux<nf; flux++) musclvars[nf+flux]=np.bs->trapezoidalm1[flux];
#endif
}


static void execute_function_on_all_segments(segmentarg_t *segmentarg, long numsegment,  int SEGMENTWORK){
  if (
#if !defined(POSIXTHREADS) && !defined(OPENMPTHREADS)
     TRUE
#else
     (SEGMENTWORK==SEGMENTWORK_LIGHT && segmentarg[0].gl->NOSHORTTHREADS)
#endif
     ){
    long cnt;
    for (cnt=0; cnt<numsegment; cnt++){
      segmentarg[cnt].funct(segmentarg[cnt].np,segmentarg[cnt].gl,segmentarg[cnt].theta,segmentarg[cnt].ls,segmentarg[cnt].le);
    }
  } else {
#ifdef POSIXTHREADS
    long cnt;
    void *retval;
    pthread_t *pthread;
    pthread=(pthread_t *)malloc((numsegment+3)*sizeof(pthread_t));
    for (cnt=0; cnt<numsegment; cnt++){
      if (pthread_create(&((pthread)[cnt]), NULL, segmentfunct, (void *)(&(segmentarg[cnt]))))
        fatal_error("Cannot create thread.");
    }
    for (cnt=0; cnt<numsegment; cnt++){
      if (pthread_join(pthread[cnt],&retval))
        fatal_error("Cannot join thread %ld.",cnt);
    }
    free(pthread);
#endif
#ifdef OPENMPTHREADS
    long cnt;
#pragma omp parallel for private(cnt) schedule(dynamic) 
    for (cnt=0; cnt<numsegment; cnt++){
      segmentarg[cnt].funct(segmentarg[cnt].np,segmentarg[cnt].gl,segmentarg[cnt].theta,segmentarg[cnt].ls,segmentarg[cnt].le);
    }
#endif
  }  
}


static void create_segments(np_t *np, gl_t *gl, long theta, long ls, long le,
            void funct(np_t *, gl_t *, long, long, long),
            segmentarg_t *segmentarg, long *cntsegment, bool COUNTFLAG, int TYPELEVEL, 
            bool is_node_valid_local(np_t, int)){
  long l,lm1,ls_local,le_local;
  bool INSIDE;
  l=ls;
  ls_local=ls; /* only needed to avoid compiler warning */
  INSIDE=FALSE;
  do {
    lm1=l;
    l=_l_plus_one(l,gl,theta);
    if ((!INSIDE) && (is_node_valid_local(np[l],TYPELEVEL))) {
      ls_local=lm1;
      INSIDE=TRUE;
    }
    if ((INSIDE) && ((!is_node_valid_local(np[l],TYPELEVEL)) || (l==le))){
      le_local=l;
      if (!COUNTFLAG) {
        segmentarg[*cntsegment].np=np;
        segmentarg[*cntsegment].gl=gl;
        segmentarg[*cntsegment].theta=theta;
        segmentarg[*cntsegment].ls=_l_plus_one(ls_local,gl,theta);
        segmentarg[*cntsegment].le=_l_minus_one(le_local,gl,theta);
        segmentarg[*cntsegment].funct=funct;
      }
      (*cntsegment)++;
      INSIDE=FALSE;
    }
  } while (l!=le);
  if (INSIDE) fatal_error("Problem setting up segments.");
}


void sweep_with_1D_segments(np_t *np, gl_t *gl, zone_t zone,
                         void funct(np_t *, gl_t *, long, long, long),
                         int sweeptype, int TYPELEVEL, bool is_node_valid_local(np_t, int),
                         int SEGMENTWORK, int GRIDLEVEL){
  long j,k,cntsegment,numthread;
  ifn1D( long i; )
  segmentarg_t *segmentarg;
  int cnt;
  bool COUNTFLAG;
  numthread=0;
  assert(is_zone_in_zone(zone,gl->domain_all));
  segmentarg=(segmentarg_t *)malloc(sizeof(segmentarg_t));
  /* do this loop twice: the first time just to count.. */
  for (cnt=0; cnt<2; cnt++){
    if (cnt==0) COUNTFLAG=TRUE; else COUNTFLAG=FALSE;
    if (!COUNTFLAG) segmentarg=(segmentarg_t *)realloc(segmentarg,numthread*sizeof(segmentarg_t));
    /* the first dimension loop */
    if (sweeptype==SWEEPTYPE_IJK || sweeptype==SWEEPTYPE_I) {
      cntsegment=0;
      for_2DL(j,zone.js,zone.je){
       if (mod(j-gl->domain_all.js,GRIDLEVEL)==0){
        for_3DL(k,zone.ks,zone.ke){
         if (mod(k-gl->domain_all.ks,GRIDLEVEL)==0){
          create_segments(np,gl,0,_ai(gl,zone.is-1,j,k),_ai(gl,zone.ie+1,j,k),
                        funct, segmentarg,&cntsegment, (bool)COUNTFLAG, TYPELEVEL,is_node_valid_local);
          if (cntsegment>=maxloopthread) {
            numthread=max(numthread,cntsegment);
            if (!COUNTFLAG) execute_function_on_all_segments(segmentarg,cntsegment,SEGMENTWORK);
            cntsegment=0;
          }
         }
        }
       }
      }
      if (cntsegment>0 && !COUNTFLAG) execute_function_on_all_segments(segmentarg,cntsegment,SEGMENTWORK);
      numthread=max(numthread,cntsegment);
    }
    /* the second dimension loop */
#ifdef _2DL
    if (sweeptype==SWEEPTYPE_IJK || sweeptype==SWEEPTYPE_J) {
      cntsegment=0;
      for_1DL(i,zone.is,zone.ie){
       if (mod(i-gl->domain_all.is,GRIDLEVEL)==0){
        for_3DL(k,zone.ks,zone.ke){
         if (mod(k-gl->domain_all.ks,GRIDLEVEL)==0){
          create_segments(np,gl,1,_ai(gl,i,zone.js-1,k),_ai(gl,i,zone.je+1,k),
                          funct, segmentarg,&cntsegment,(bool)COUNTFLAG, TYPELEVEL,is_node_valid_local);
          if (cntsegment>=maxloopthread) {
            numthread=max(numthread,cntsegment);
            if (!COUNTFLAG) execute_function_on_all_segments(segmentarg,cntsegment,SEGMENTWORK);
            cntsegment=0;
          }
         }
        }
       }
      }
      if (cntsegment>0 && !COUNTFLAG) execute_function_on_all_segments(segmentarg,cntsegment,SEGMENTWORK);
      numthread=max(numthread,cntsegment);
    }
#endif
    /* the third dimension loop */
#ifdef _3DL
    if (sweeptype==SWEEPTYPE_IJK || sweeptype==SWEEPTYPE_K) {
      cntsegment=0;
      for_1DL(i,zone.is,zone.ie){
       if (mod(i-gl->domain_all.is,GRIDLEVEL)==0){
        for_2DL(j,zone.js,zone.je){
         if (mod(j-gl->domain_all.js,GRIDLEVEL)==0){
          create_segments(np,gl,2,_ai(gl,i,j,zone.ks-1),_ai(gl,i,j,zone.ke+1),
                          funct, segmentarg, &cntsegment,(bool)COUNTFLAG, TYPELEVEL,is_node_valid_local);
          if (cntsegment>=maxloopthread) {
            numthread=max(numthread,cntsegment);
            if (!COUNTFLAG) execute_function_on_all_segments(segmentarg,cntsegment,SEGMENTWORK);
            cntsegment=0;
          }
         }
        }
       }
      }
      if (cntsegment>0 && !COUNTFLAG) execute_function_on_all_segments(segmentarg,cntsegment,SEGMENTWORK);
      numthread=max(numthread,cntsegment);
    }
#endif
  }
  free(segmentarg);
}


/* the following first sets the offset to 0, then 1, then -1 */
static long _node_offset_from_cnt(long cnt){
  long offset;
  offset=0;
  if (cnt==0) offset=0;
  if (cnt==1) offset=1;
  if (cnt==2) offset=-1;
  return(offset);
}


void update_bdry_node(np_t *np, gl_t *gl, long l){
  long dim,dimsgn,l_C,l_B,l_A,l_D;
  bool BDRYDIRECFOUND;
#ifdef _2DL
  long offset1,offset2,cnt1,cnt2;
#endif
#ifdef _3D
  long offset3,cnt3;
#endif
  bool UPDATED;
  assert(is_node_bdry(np[l],TYPELEVEL_FLUID_WORK));
  UPDATED=FALSE;
  BDRYDIRECFOUND=find_bdry_direc(np, gl, l, TYPELEVEL_FLUID_WORK, &dim, &dimsgn);

  if (is_node_link(np[l],TYPELEVEL_FLUID_WORK)) {
    // in case the boundary node is a link, U has already been updated: simply update the prim variables 
    find_prim_fluid(np, l, gl);
    UPDATED=TRUE;
  }

  if (BDRYDIRECFOUND && !UPDATED){
    l_A=l;
    l_B=_al(gl,l,dim,dimsgn);
    l_C=_al(gl,l,dim,dimsgn*2);
    if (is_node_inner(np[_al(gl,l,dim,dimsgn*3)],TYPELEVEL_FLUID_WORK)) l_D=_al(gl,l,dim,dimsgn*3);
      else l_D=l_C;
    assert(is_node_inner(np[l_C],TYPELEVEL_FLUID_WORK));
    assert(is_node_inner(np[l_B],TYPELEVEL_FLUID_WORK));
    update_bdry_fluid(np,gl,l_A,l_B,l_C,l_D,dim,dimsgn,BDRYDIRECFOUND,TYPELEVEL_FLUID_WORK);
    UPDATED=TRUE;
  }

  /* now, do the corners */
  if (!UPDATED) {
#ifdef _2D
    for (cnt1=0; cnt1<=2; cnt1++){
      for (cnt2=0; cnt2<=2; cnt2++){
        offset1=_node_offset_from_cnt(cnt1);
        offset2=_node_offset_from_cnt(cnt2);
        l_C=_all(gl,l,0,offset1*2,1,offset2*2);
        l_B=_all(gl,l,0,offset1,1,offset2);
        l_A=l;
        l_D=l_C;
        if  (  is_node_inner(np[l_B],TYPELEVEL_FLUID_WORK)
            && is_node_inner(np[l_C],TYPELEVEL_FLUID_WORK) && !UPDATED){
          update_bdry_fluid(np,gl,l_A,l_B,l_C,l_D,dim,dimsgn,BDRYDIRECFOUND,TYPELEVEL_FLUID_WORK);
          UPDATED=TRUE;
        }
      }
    }
#endif
#ifdef _3D
    for (cnt1=0; cnt1<=2; cnt1++){
      for (cnt2=0; cnt2<=2; cnt2++){
        for (cnt3=0; cnt3<=2; cnt3++){
          offset1=_node_offset_from_cnt(cnt1);
          offset2=_node_offset_from_cnt(cnt2);
          offset3=_node_offset_from_cnt(cnt3);
          l_C=_al(gl,
                      _al(gl,
                         _al(gl,l,0,offset1*2),
                       1,offset2*2),
                    2,offset3*2);
          l_B=_al(gl,
                      _al(gl,
                         _al(gl,l,0,offset1),
                       1,offset2),
                    2,offset3);
          l_A=l;
          l_D=l_C;
          if  (  is_node_inner(np[l_B],TYPELEVEL_FLUID_WORK)
              && is_node_inner(np[l_C],TYPELEVEL_FLUID_WORK) && !UPDATED){
            update_bdry_fluid(np,gl,l_A,l_B,l_C,l_D,dim,dimsgn,BDRYDIRECFOUND,TYPELEVEL_FLUID_WORK);
            UPDATED=TRUE;
          }
        }
      }
    }
#endif
  }
  if (!UPDATED) {
    fatal_error("Problem updating boundary node in update_bdry_node() function.");
  }
}


void update_bdry_nodes_on_segment(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (is_node_bdry(np[l],TYPELEVEL_FLUID_WORK)){
      thread_lock_node_set(np,l,THREADTYPE_ZONE);
      update_bdry_node(np, gl, l);
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
  }
}


void update_bdry_nodes(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np, gl, zone, &update_bdry_nodes_on_segment, SWEEPTYPE_I, TYPELEVEL_FLUID_WORK,&is_node_valid,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
}


#ifdef DISTMPI

#define numfluidvars (nf+1+max(0,hbw_resconv_fluid-1)*nmc)
#define numlinkvars ((hbw_resconv_fluid-1)*nmc)
#define DOUBLE_INT_MAX 100000000000000

typedef double sendvars_t[max(nfe,numfluidvars)];
typedef struct {
  sendvars_t vars;
  int proc;
  long l;
  bool SENT;
} sendnode_t;


void update_linked_nodes_2(np_t *np, gl_t *gl, int TYPELEVEL){

  int rankrecv,numproc,ranksend,thisrank;
  long i,j,k;
  zone_t zonesend,zonerecv,zone;
  MPI_Status MPI_Status1;

  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

  /* here we need to mpi the linkmusclvars */
  for (ranksend=0; ranksend<numproc; ranksend++){
    zonesend=_domain_from_rank(ranksend,gl);
    for (rankrecv=0; rankrecv<numproc; rankrecv++){
      if (rankrecv!=ranksend && (ranksend==thisrank || rankrecv==thisrank)){
        zonerecv=_domain_lim_from_rank(rankrecv,gl);
        if (is_zone_intersecting_zone(zonesend,zonerecv)){
          zone=_zone_intersection(zonesend,zonerecv);
          for_ijk(zone,is,js,ks,ie,je,ke){
                if (ranksend==thisrank) {
//                  if (is_node_link(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)) printf("x");
                  if (is_node_link(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK) && is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)){
                    assert(np[_ai(gl,i,j,k)].numlinkmusclvars!=0);
                    assert(np[_ai(gl,i,j,k)].linkmusclvars!=NULL);
                    MPI_Send(&np[_ai(gl,i,j,k)].numlinkmusclvars,1,MPI_INT,rankrecv,0,MPI_COMM_WORLD);
                    MPI_Send(np[_ai(gl,i,j,k)].linkmusclvars,numlinkvars,MPI_DOUBLE,rankrecv,0,MPI_COMM_WORLD);
                    
                  }
                }
                if (rankrecv==thisrank) {
                  if (is_node_link(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK) && is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)){
                    MPI_Recv(&np[_ai(gl,i,j,k)].numlinkmusclvars,1,MPI_INT,ranksend,0,MPI_COMM_WORLD,&MPI_Status1);
                    assert(np[_ai(gl,i,j,k)].linkmusclvars!=NULL);
                    MPI_Recv(np[_ai(gl,i,j,k)].linkmusclvars,numlinkvars,MPI_DOUBLE,ranksend,0,MPI_COMM_WORLD,&MPI_Status1);
                  }
                }
          }
        }
      }
    }
  } 
   
  MPI_Barrier(MPI_COMM_WORLD);
}


void update_linked_nodes(np_t *np, gl_t *gl, int TYPELEVEL){
  long i,j,k,l1,l2,flux,offset,l,cntlink;
  MPI_Status MPI_Status1;
  musclvarscycle_t musclvars;
  sendvars_t mpivars;
  int thisrank,numproc,rank2,rank1,thisproc;
  int packsize,buffersize,bbuffersize;
  double *buffer,*bbuffer;
  sendnode_t *sendnode;
  long numsendvars,numvars,numsend,cntsend,cnt;
  double *sendvars;
  int *recvproc;
  int cntproc;
  bool FOUND;
  zone_t zone;

  zone=gl->domain;
  switch (TYPELEVEL){
    case TYPELEVEL_FLUID:
      numvars=numfluidvars;
    break;
    case TYPELEVEL_FLUID_WORK:
      numvars=numfluidvars;
    break;
#ifdef EMFIELD
    case TYPELEVEL_EMFIELD:
      numvars=nfe;
    break;
#endif
    default:
      fatal_error("TYPELEVEL can not be set to %d.\n",TYPELEVEL);
      numvars=0;
  }

  sendnode=(sendnode_t *)malloc(sizeof(sendnode_t));
  sendvars=(double *)malloc(sizeof(double));
  cntsend=0;

  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &packsize );
  
  recvproc=(int *)malloc((numproc+2)*sizeof(int));
  buffersize = min(INT_MAX,nmc*(zone.ie-zone.is)*(zone.je-zone.js)if3DL(*(zone.ke-zone.ks)) * (MPI_BSEND_OVERHEAD + packsize));
  buffer = (double *)malloc( buffersize );

  MPI_Buffer_attach( buffer, buffersize );


  for_ijk(zone,is,js,ks,ie,je,ke){
        np[_ai(gl,i,j,k)].numlinkmusclvars=0;
  }  

  /* first send the packets */
  cntsend=0;
  for_ijk(zone,is,js,ks,ie,je,ke){
        if (is_node_link(np[_ai(gl,i,j,k)],TYPELEVEL)){
#ifdef _CYCLE_MULTIZONE
          fatal_error("Linked nodes can not be used with Multizone cycle yet. Need to update update_linked_nodes() function.");
#endif
#ifdef _CYCLE_MULTIZONE_MARCHING
          fatal_error("Linked nodes can not be used with MultizoneMarching cycle yet. Need to update update_linked_nodes() function.");
#endif
          if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL)){
            for (cntlink=0; cntlink<_num_node_link(np[_ai(gl,i,j,k)],TYPELEVEL); cntlink++){
              l1=_ai_all(gl,i,j,k);
              l2=_node_link(np[_ai(gl,i,j,k)],cntlink,TYPELEVEL);
              rank1=_node_rank(gl, i, j, k);
              rank2=_node_rank(gl, _i_all(l2,gl,0), _i_all(l2,gl,1), _i_all(l2,gl,2));
           
              if (rank1==thisrank) {
                if (TYPELEVEL==TYPELEVEL_FLUID_WORK || TYPELEVEL==TYPELEVEL_FLUID){
                  for (flux=0; flux<nf; flux++) mpivars[flux]=np[_l_from_l_all(gl,l1)].bs->U[flux];
                  mpivars[nf]=(double)_nodes_between_link_and_bdry_limited(np, gl, _l_from_l_all(gl,l1), l2, TYPELEVEL, max(0,hbw_resconv_fluid-1));
                  for (offset=1; offset<hbw_resconv_fluid; offset++) {
//                  find_prim_fluid(np, _al_link(np, gl, _l_from_l_all(gl,l1), offset, TYPELEVEL), gl);
                    find_musclvarscycle(np[_al_link(np, gl, _l_from_l_all(gl,l1), l2, offset, TYPELEVEL)], gl, musclvars);
                    for (flux=0; flux<nmc; flux++) mpivars[1+flux+nf+(offset-1)*nmc]=musclvars[flux];
                  }
                  if (rank1!=rank2){
                    for (flux=0; flux<numvars; flux++) sendnode[cntsend].vars[flux]=mpivars[flux];
                    sendnode[cntsend].proc=(int)rank2;
                    sendnode[cntsend].l=l2;
                    sendnode[cntsend].SENT=FALSE;
                    cntsend++;
                    sendnode=(sendnode_t *)realloc(sendnode,(cntsend+1)*sizeof(sendnode_t));
                  } else {
                    /* no need to send with MPI*/
                    //printf("\n --(%ld,%ld,%ld) %d",i,j,k,thisrank);

                    l=_l_from_l_all(gl,l2);
                    for (flux=0; flux<nf; flux++) np[l].bs->U[flux]=mpivars[flux];
                    assert(np[l].linkmusclvars!=NULL);
                    assert(is_node_bdry(np[l],TYPELEVEL));
                    assert(is_node_link(np[l],TYPELEVEL));
                    np[l].numlinkmusclvars=(short)round(mpivars[nf]);
                    for (offset=1; offset<hbw_resconv_fluid; offset++) {
                      for (flux=0; flux<nmc; flux++) np[l].linkmusclvars[flux+(offset-1)*nmc]=mpivars[1+flux+nf+(offset-1)*nmc];
                    }
                   
                  } 
                }
#ifdef EMFIELD
                if (TYPELEVEL==TYPELEVEL_EMFIELD){
                  for (flux=0; flux<numvars; flux++) mpivars[flux]=np[_l_from_l_all(gl,l1)].bs->Uemfield[flux];
                  if (rank1!=rank2) {
                    for (flux=0; flux<numvars; flux++) sendnode[cntsend].vars[flux]=mpivars[flux];
                    sendnode[cntsend].proc=(int)rank2;
                    sendnode[cntsend].l=l2;
                    sendnode[cntsend].SENT=FALSE;
                    cntsend++;
                    sendnode=(sendnode_t *)realloc(sendnode,(cntsend+1)*sizeof(sendnode_t));
                  } else {
                    /* no need to send with MPI */
                    for (flux=0; flux<nfe; flux++) np[_l_from_l_all(gl,l2)].bs->Uemfield[flux]=mpivars[flux];
                  }
                }
#endif

              }
            }
          }
        }
  }

  numsend=cntsend;


  /* send nodes in block one proc at a time */
  
  do {
    thisproc=-1;
    numsendvars=0;
    for (cntsend=0; cntsend<numsend; cntsend++){
      if (thisproc==-1 && !sendnode[cntsend].SENT) thisproc=sendnode[cntsend].proc;
      if (sendnode[cntsend].proc==thisproc){
        assert(!sendnode[cntsend].SENT);
        sendvars=(double *)realloc(sendvars,(numsendvars+2*numvars)*sizeof(double));
        for (flux=0; flux<numvars; flux++) sendvars[numsendvars+flux]=sendnode[cntsend].vars[flux];
        numsendvars+=numvars;
#ifndef NDEBUG
        sendvars[numsendvars]=(double)mod(sendnode[cntsend].l,DOUBLE_INT_MAX);
        numsendvars++;
#endif
        sendnode[cntsend].SENT=TRUE;
      }
    }
    if (thisproc!=-1){
      if (MPI_Bsend(&numsendvars,1,MPI_LONG,thisproc,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("Problem with MPI_Bsend in update_linked_nodes().");    
      if (MPI_Bsend(sendvars,numsendvars,MPI_DOUBLE,thisproc,0,MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("Problem with MPI_Bsend in update_linked_nodes().");    
    }
  } while (thisproc!=-1);


  for (cnt=0; cnt<(numproc+2); cnt++){
    recvproc[cnt]=-1;
  }


  for_ijk(zone,is,js,ks,ie,je,ke){
        if (is_node_link(np[_ai(gl,i,j,k)],TYPELEVEL) && is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL)){
          l1=_node_link(np[_ai(gl,i,j,k)],0,TYPELEVEL);
          rank2=_node_rank(gl, i, j, k);
          rank1=_node_rank(gl, _i_all(l1,gl,0), _i_all(l1,gl,1), _i_all(l1,gl,2));
          if (rank1!=rank2 && rank2==thisrank){
            /* rank1 is one process that we will need to get data from; store it in recvproc */
            cntproc=0;
            while(recvproc[cntproc]!=-1 && recvproc[cntproc]!=rank1 ) {
              cntproc++;
            } 
            assert(cntproc<numproc);
            recvproc[cntproc]=rank1;
          }
        }
  }


  cntproc=0;
  while (recvproc[cntproc]!=-1) {
    thisproc=recvproc[cntproc];
    MPI_Recv(&numsendvars,1,MPI_LONG,thisproc,0,MPI_COMM_WORLD,&MPI_Status1);
    sendvars=(double *)realloc(sendvars,numsendvars*sizeof(double));
    MPI_Recv(sendvars,numsendvars,MPI_DOUBLE,thisproc,0,MPI_COMM_WORLD,&MPI_Status1);
    cntsend=0;
    for_ijk(zone,is,js,ks,ie,je,ke){
          if (is_node_link(np[_ai(gl,i,j,k)],TYPELEVEL) && is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL)){
            l2=_ai_all(gl,i,j,k);
            assert(is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL));
            l1=_node_link(np[_ai(gl,i,j,k)],0,TYPELEVEL);
            rank2=_node_rank(gl, i, j, k);
            rank1=_node_rank(gl, _i_all(l1,gl,0), _i_all(l1,gl,1), _i_all(l1,gl,2));
            if (rank1!=rank2 && rank2==thisrank){
              if (thisproc==rank1){
                cntsend=0;
                FOUND=FALSE;
                do {
                  if (mod(l2,DOUBLE_INT_MAX)==(long)sendvars[cntsend+numvars]) FOUND=TRUE;
                   else cntsend+=numvars+1;
                } while (!FOUND && cntsend<numsendvars);
                if (!FOUND) fatal_error("Problem finding node number passed by MPI. "); 
                for (flux=0; flux<numvars; flux++) mpivars[flux]=sendvars[cntsend+flux];
                
#ifndef NDEBUG
                assert(mod(l2,DOUBLE_INT_MAX)==(long)sendvars[cntsend+numvars]);                
#endif
                l=_l_from_l_all(gl,l2);
                assert(is_node_bdry(np[l],TYPELEVEL));
                assert(is_node_link(np[l],TYPELEVEL));
                if (TYPELEVEL==TYPELEVEL_FLUID_WORK || TYPELEVEL==TYPELEVEL_FLUID){
                  for (flux=0; flux<nf; flux++) np[l].bs->U[flux]=mpivars[flux];
                  assert(np[l].linkmusclvars!=NULL);
                  np[l].numlinkmusclvars=(short)round(mpivars[nf]);
                  for (offset=1; offset<hbw_resconv_fluid; offset++) {
                    for (flux=0; flux<nmc; flux++)
                      np[l].linkmusclvars[flux+(offset-1)*nmc]=mpivars[1+flux+nf+(offset-1)*nmc];
                  }
                }
#ifdef EMFIELD
                if (TYPELEVEL==TYPELEVEL_EMFIELD){
                  for (flux=0; flux<nfe; flux++) np[l].bs->Uemfield[flux]=mpivars[flux];
                }
#endif
              }
            }
          }
    }
    cntproc++;
  }

  MPI_Buffer_detach( &bbuffer, &bbuffersize );
  free(buffer);
  MPI_Barrier(MPI_COMM_WORLD);
  free(sendnode);
  free(recvproc);
  free(sendvars);

  update_linked_nodes_2(np, gl, TYPELEVEL);
}


#else//DISTMPI


void update_linked_nodes(np_t *np, gl_t *gl, int TYPELEVEL){
  long i,j,k,l1,l2,flux;
  for_ijk(gl->domain,is,js,ks,ie,je,ke){
        l1=_ai(gl,i,j,k);
        if (is_node_bdry(np[l1],TYPELEVEL) && is_node_link(np[l1],TYPELEVEL)){
#ifdef _CYCLE_MULTIZONE
          fatal_error("Linked nodes can not be used with Multizone cycle yet. Need to update update_linked_nodes() function.");
#endif
#ifdef _CYCLE_MULTIZONE_MARCHING
          fatal_error("Linked nodes can not be used with MultizoneMarching cycle yet. Need to update update_linked_nodes() function.");
#endif
          assert(is_node_bdry(np[l1],TYPELEVEL));
          l2=_node_link(np[l1],0,TYPELEVEL);
          if (TYPELEVEL==TYPELEVEL_FLUID_WORK || TYPELEVEL==TYPELEVEL_FLUID){
            for (flux=0; flux<nf; flux++) np[l1].bs->U[flux]=np[l2].bs->U[flux];
          }
#ifdef EMFIELD
          if (TYPELEVEL==TYPELEVEL_EMFIELD){
            for (flux=0; flux<nfe; flux++) np[l1].bs->Uemfield[flux]=np[l2].bs->Uemfield[flux];
          }
#endif
        }
  }
}


#endif//DISTMPI


static bool is_node_in_region(bool(*FUNCT)(gl_t *, long, long, long),
                         gl_t *gl, long i, long j, long k){
  bool tmp;
  tmp=FUNCT(gl,i,j,k);
  return(tmp);
}


static bool is_node_in_region_extended_by_bb(bool(*FUNCT)(gl_t *, long, long, long),
                         gl_t *gl, long i, long j, long k){
  bool tmp;
  long cnti,cntj,cntk;
  tmp=FALSE;
  for_1DL(cnti,i-hbw_bdry_fluid,i+hbw_bdry_fluid){
    for_2DL(cntj,j-hbw_bdry_fluid,j+hbw_bdry_fluid){
      for_3DL(cntk,k-hbw_bdry_fluid,k+hbw_bdry_fluid){
        if (FUNCT(gl,cnti,cntj,cntk)) tmp=TRUE;
      }
    }
  }
  return(tmp);
}


void resume_nodes_specified_in_function(np_t *np, gl_t *gl,
                                bool(*FUNCT)(gl_t *, long, long, long)){
  long i,j,k;
  long *noderes;
  long *bdryres;
  long numnoderes,numbdryres,cnt;
  copy_base_to_work_node_type(np,gl,gl->domain_lim);
  noderes=(long *)malloc((gl->domain.ie-gl->domain.is+4)if2DL(*(gl->domain.je-gl->domain.js+4))
                         if3DL(*(gl->domain.ke-gl->domain.ks+4))*sizeof(long));
  bdryres=(long *)malloc((gl->domain.ie-gl->domain.is+4)if2DL(*(gl->domain.je-gl->domain.js+4))
                         if3DL(*(gl->domain.ke-gl->domain.ks+4))*sizeof(long));
  numnoderes=0;
  numbdryres=0;
  for_ijk(gl->domain,is-1,js-1,ks-1,ie+1,je+1,ke+1){
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)
           && is_node_in_region_extended_by_bb(FUNCT,gl,i,j,k)) {
          if (resume_node(&(np[_ai(gl,i,j,k)]))  ) {
            if (is_node_in_region(FUNCT,gl,i,j,k) ){
              bdryres[numbdryres]=_ai(gl,i,j,k);
              numbdryres++;
            }
            noderes[numnoderes]=_ai(gl,i,j,k);
            numnoderes++;
          }
        } else {
          suspend_node(&(np[_ai(gl,i,j,k)]));
        }
  }
  /* rebuild the working variables of the inner nodes of the nodes resumed*/
  for (cnt=0; cnt<numnoderes; cnt++){
    if (is_node_resumed(np[noderes[cnt]]) && is_node_inner(np[noderes[cnt]],TYPELEVEL_FLUID)){
       find_prim_fluid(np,noderes[cnt],gl);
    }
  }
  /* rebuild the working variables of the boundary nodes of the nodes resumed*/
  for (cnt=0; cnt<numbdryres; cnt++){
    if (is_node_resumed(np[bdryres[cnt]]) && is_node_bdry(np[bdryres[cnt]],TYPELEVEL_FLUID)) {
       update_bdry_node(np,gl,bdryres[cnt]);
    }
  }
  /* suspend all nodes needed only to compute the boundary nodes.
     this is necessary to ensure that all non-suspended nodes are properly updated.*/
  for_ijk(gl->domain,is-1,js-1,ks-1,ie+1,je+1,ke+1){
        if (!(is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID) && is_node_in_region(FUNCT,gl,i,j,k)))
          suspend_node(&(np[_ai(gl,i,j,k)]));
  }
  free(noderes);
  free(bdryres);
}


void resume_nodes_only_in_zone_and_update_bdry_nodes(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;
  long *noderes;
  long *bdryres;
  long numnoderes,numbdryres,cnt;
  copy_base_to_work_node_type(np,gl,gl->domain_lim);
  noderes=(long *)malloc((gl->domain_lim.ie-gl->domain_lim.is+1)
                         if2DL(*(gl->domain_lim.je-gl->domain_lim.js+1))
                         if3DL(*(gl->domain_lim.ke-gl->domain_lim.ks+1))*sizeof(long));
  bdryres=(long *)malloc((gl->domain_lim.ie-gl->domain_lim.is+1)
                         if2DL(*(gl->domain_lim.je-gl->domain_lim.js+1))
                         if3DL(*(gl->domain_lim.ke-gl->domain_lim.ks+1))*sizeof(long));
  numnoderes=0;
  numbdryres=0;
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)
                 && (i>=zone.is-hbw_bdry_fluid) && (i<=zone.ie+hbw_bdry_fluid)
           if2DL(&& (j>=zone.js-hbw_bdry_fluid) && (j<=zone.je+hbw_bdry_fluid))
           if3DL(&& (k>=zone.ks-hbw_bdry_fluid) && (k<=zone.ke+hbw_bdry_fluid))) {
          if (resume_node(&(np[_ai(gl,i,j,k)]))  ) {
            if (is_node_in_zone(i,j,k,zone)){
              bdryres[numbdryres]=_ai(gl,i,j,k);
              numbdryres++;
            }
            noderes[numnoderes]=_ai(gl,i,j,k);
            numnoderes++;
          }
        } else {
          suspend_node(&(np[_ai(gl,i,j,k)]));
        }
  }
  /* rebuild the working variables of the inner nodes of the nodes resumed*/
#ifdef OPENMPTHREADS
#pragma omp parallel for private(cnt) schedule(dynamic)
#endif
  for (cnt=0; cnt<numnoderes; cnt++){
    if (is_node_resumed(np[noderes[cnt]]) && is_node_inner(np[noderes[cnt]],TYPELEVEL_FLUID)){
       find_prim_fluid(np,noderes[cnt],gl);
    }
  }
  free(noderes);
  /* rebuild the working variables of the boundary nodes of the nodes resumed*/
#ifdef OPENMPTHREADS
#pragma omp parallel for private(cnt) schedule(dynamic)
#endif
  for (cnt=0; cnt<numbdryres; cnt++){
    if (is_node_resumed(np[bdryres[cnt]]) && is_node_bdry(np[bdryres[cnt]],TYPELEVEL_FLUID)) {
      find_ijk_from_l(gl, bdryres[cnt], &i, &j, &k);
      if (is_node_in_zone(i, j, k, gl->domain)){
        update_bdry_node(np,gl,bdryres[cnt]);
      } else {
        find_prim_fluid(np,bdryres[cnt],gl);
      }
    }
  }
  free(bdryres);
  /* suspend all nodes needed only to compute the boundary nodes.
     this is necessary to ensure that all non-suspended nodes are properly updated.*/
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        if (!(is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)
           && is_node_in_zone(i,j,k,zone)))
          suspend_node(&(np[_ai(gl,i,j,k)]));
  }
}


void resume_nodes_in_zone(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;
  long *noderes;
  long numnoderes,cnt;
  zone_t zoneint;
  copy_base_to_work_node_type(np,gl,gl->domain_lim);
  noderes=(long *)malloc((gl->domain_lim.ie-gl->domain_lim.is+1)
                         if2DL(*(gl->domain_lim.je-gl->domain_lim.js+1))
                         if3DL(*(gl->domain_lim.ke-gl->domain_lim.ks+1))*sizeof(long));
  numnoderes=0;
  zoneint=_zone_intersection(gl->domain_lim,zone);
  for_ijk(zoneint,is,js,ks,ie,je,ke){
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)) {
          if (resume_node(&(np[_ai(gl,i,j,k)]))  ) {
            noderes[numnoderes]=_ai(gl,i,j,k);
            numnoderes++;
          }
        } 
  }
  /* rebuild the working variables of the inner nodes of the nodes resumed*/
#ifdef OPENMPTHREADS
#pragma omp parallel for private(cnt) schedule(dynamic)
#endif
  for (cnt=0; cnt<numnoderes; cnt++){
    if (is_node_resumed(np[noderes[cnt]]) && is_node_valid(np[noderes[cnt]],TYPELEVEL_FLUID)){
       find_prim_fluid(np,noderes[cnt],gl);
    }
  }
  free(noderes);
}


void resume_nodes_only_in_zone(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;
  long *noderes;
  long numnoderes,cnt;
  copy_base_to_work_node_type(np,gl,gl->domain_lim);
  noderes=(long *)malloc((gl->domain_lim.ie-gl->domain_lim.is+1)
                         if2DL(*(gl->domain_lim.je-gl->domain_lim.js+1))
                         if3DL(*(gl->domain_lim.ke-gl->domain_lim.ks+1))*sizeof(long));
  numnoderes=0;
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)
                 && (i>=zone.is) && (i<=zone.ie)
           if2DL(&& (j>=zone.js) && (j<=zone.je))
           if3DL(&& (k>=zone.ks) && (k<=zone.ke))) {
          if (resume_node(&(np[_ai(gl,i,j,k)]))  ) {
            noderes[numnoderes]=_ai(gl,i,j,k);
            numnoderes++;
          }
        } else {
          suspend_node(&(np[_ai(gl,i,j,k)]));
        }
  }

  /* rebuild the working variables of the inner nodes of the nodes resumed*/
#ifdef OPENMPTHREADS
#pragma omp parallel for private(cnt) schedule(dynamic)
#endif
  for (cnt=0; cnt<numnoderes; cnt++){
    if (is_node_resumed(np[noderes[cnt]]) && is_node_valid(np[noderes[cnt]],TYPELEVEL_FLUID)){
       find_prim_fluid(np,noderes[cnt],gl);
    }
  }
  free(noderes);
  /* suspend all nodes needed only to compute the boundary nodes.
     this is necessary to ensure that all non-suspended nodes are properly updated.*/
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        if (!(is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)
           && is_node_in_zone(i,j,k,zone)))
          suspend_node(&(np[_ai(gl,i,j,k)]));
  }
}


#ifdef UNSTEADY



void increase_time_level(np_t *np, gl_t *gl){
  long i,j,k,flux,l;
  gl->time+=gl->dt;
  gl->iter=0;
  add_double_to_codex(&(gl->cycle.codex),"time",gl->time);
  gl->timelevelcnt++;  
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        l=_ai(gl,i,j,k);
        if ((is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID))){
          for (flux=0; flux<nf; flux++){
#if _RESTIME_BW > 3
            np[l].bs->Um3[flux]=np[l].bs->Um2[flux];
#endif
#if _RESTIME_BW > 2
            np[l].bs->Um2[flux]=np[l].bs->Um1[flux];
#endif
            np[l].bs->Um1[flux]=np[l].bs->U[flux];
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
            np[l].bs->trapezoidalm1[flux]=np[l].bs->trapezoidalm1_next[flux];
#endif
          }
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_MUSCLVARS
          find_musclvars(np[l],gl,np[l].bs->trapezoidalm1);
#endif
        }
#ifdef EMFIELD
        if ((is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD))){
          for (flux=0; flux<nfe; flux++){
            np[l].bs->Uemfieldm1[flux]=np[l].bs->Uemfield[flux];
          }
        }
#endif
  }
}
#endif//UNSTEADY



void runtime_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  char *oldfilename;
  
  oldfilename=(char *)malloc(sizeof(char)*(5+strlen((((readcontrolarg_t *)codex->action_args)->gl->output_filename))));
  strcpy(oldfilename,(((readcontrolarg_t *)codex->action_args)->gl->output_filename));
  if (strcmp(actionname,"WriteDataFile")==0) {
    if (SOAP_number_argums(*argum)==1){
      SOAP_substitute_all_argums(argum, codex);         
      SOAP_get_argum_string(codex,&(((readcontrolarg_t *)codex->action_args)->gl->output_filename),*argum,0);
    }
    if (SOAP_number_argums(*argum)>1){
      SOAP_fatal_error(codex,"Action WriteDataFile() can not be called with more than 1 argument. Either it is called with one argument (a string containing the data file name) or with no argument. If no argument is given, the default data file name as specified on the command line will be used.");
    }
    write_data_file(*((readcontrolarg_t *)codex->action_args)->np,
                    ((readcontrolarg_t *)codex->action_args)->gl);
    codex->ACTIONPROCESSED=TRUE;
  }
  strcpy((((readcontrolarg_t *)codex->action_args)->gl->output_filename),oldfilename);
  free(oldfilename);
  if (strcmp(actionname,"Init")==0) {
    read_init(*argum, codex);
    codex->action=&runtime_actions;
    ((readcontrolarg_t *)codex->action_args)->gl->RESIDUAL_ALTERED=TRUE;
#ifdef EMFIELD
    ((readcontrolarg_t *)codex->action_args)->gl->RESIDUAL_ALTERED_EMFIELD=TRUE;
#endif
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Model")==0) {
    read_model(*argum, codex);
    codex->action=&runtime_actions;
    ((readcontrolarg_t *)codex->action_args)->gl->RESIDUAL_ALTERED=TRUE;
#ifdef EMFIELD
    ((readcontrolarg_t *)codex->action_args)->gl->RESIDUAL_ALTERED_EMFIELD=TRUE;
#endif
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Disc")==0) {
    read_disc(*argum, codex);
    codex->action=&runtime_actions;
    ((readcontrolarg_t *)codex->action_args)->gl->RESIDUAL_ALTERED=TRUE;
#ifdef EMFIELD
    ((readcontrolarg_t *)codex->action_args)->gl->RESIDUAL_ALTERED_EMFIELD=TRUE;
#endif
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Cycle")==0) {
    read_cycle(*argum, codex);
    codex->action=&runtime_actions;
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Post")==0) {
    read_post(*argum, codex);
    codex->action=&runtime_actions;
    codex->ACTIONPROCESSED=TRUE;
  }
  runtime_actions_cycle_specific(actionname,argum,codex);
}


void write_cycle_template(FILE **controlfile){
  wfprintf(*controlfile,
  "\n\n"
  "Cycle(\n"
  );
  write_cycle_fluid_template(controlfile);
#ifdef EMFIELD
  write_cycle_emfield_template(controlfile);
#endif
  write_runtime_template(controlfile);
  wfprintf(*controlfile,
  ");\n"
  );
}


void read_cycle_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl;
  gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_CYCLE_ACTIONNAME)==0 && !gl->CONTROL_READ) {
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_CYCLE_ACTIONNAME);
    ((readcontrolarg_t *)codex->action_args)->
       gl->cycle.code_runtime=(char *)malloc((strlen(*argum)+2)*sizeof(char));
    strcpy(((readcontrolarg_t *)codex->action_args)->gl->cycle.code_runtime,*argum);
    ((readcontrolarg_t *)codex->action_args)->gl->cycle.RUNTIMEMODULEFOUND=TRUE;
    codex->ACTIONPROCESSED=TRUE;
  }
  read_cycle_fluid_actions(actionname, argum, codex);
  read_cycle_emfield_actions(actionname, argum, codex);
}


void read_cycle(char *argum, SOAP_codex_t *codexcontrol){
  gl_t *gl;
  gl=((readcontrolarg_t *)codexcontrol->action_args)->gl;
  if (!gl->CONTROL_READ){
    gl->cycle.RUNTIMEMODULEFOUND=FALSE;
  }
  codexcontrol->action=&read_cycle_actions;
  SOAP_process_code(argum, codexcontrol, SOAP_VARS_KEEP_ALL);
  if (!gl->CONTROL_READ){
    if (!gl->CYCLE_FLUID_READ) 
      fatal_error("The fluid module %s() was not found within Cycle().",_FLUID_ACTIONNAME);
    if (!gl->CYCLE_EMFIELD_READ) 
      fatal_error("The emfield module %s() was not found within Cycle().",_EMFIELD_ACTIONNAME);
    if (!gl->cycle.RUNTIMEMODULEFOUND)
      fatal_error("The module %s() was not found within Cycle().",_CYCLE_ACTIONNAME);
    init_cycle(argum,codexcontrol);
  }
}


void write_disc_template(FILE **controlfile){
  wfprintf(*controlfile,
  "\n\n"
  "Disc(\n"
  );
  write_disc_fluid_template(controlfile);
#ifdef EMFIELD
  write_disc_emfield_template(controlfile);
#endif
  write_disc_resconv_template(controlfile);
  write_disc_restime_template(controlfile);
  wfprintf(*controlfile,
  ");\n"
  );
}


void read_disc_actions(char *actionname, char **argum, SOAP_codex_t *codex){
//  gl_t *gl;
//  gl=((readcontrolarg_t *)codex->action_args)->gl;
  read_disc_fluid_actions(actionname, argum, codex);
  read_disc_emfield_actions(actionname, argum, codex);
  read_disc_resconv_actions(actionname, argum, codex);
  read_disc_restime_actions(actionname, argum, codex);
}


void read_disc(char *argum, SOAP_codex_t *codexcontrol){
  gl_t *gl;
  gl=((readcontrolarg_t *)codexcontrol->action_args)->gl;
  codexcontrol->action=&read_disc_actions;
  gl->DISC_FLUID_READ=FALSE;
  gl->DISC_EMFIELD_READ=FALSE;
  gl->DISC_RESCONV_READ=FALSE;
  gl->DISC_RESTIME_READ=FALSE;
  SOAP_process_code(argum, codexcontrol, SOAP_VARS_KEEP_ALL);
  if (!gl->CONTROL_READ){
    if (!gl->DISC_FLUID_READ) 
      fatal_error("The fluid module %s() was not found within Disc().",_FLUID_ACTIONNAME);
    if (!gl->DISC_EMFIELD_READ) 
      fatal_error("The emfield module %s() was not found within Disc().",_EMFIELD_ACTIONNAME);
    if (!gl->DISC_RESCONV_READ) 
      fatal_error("The residual convection module %s() was not found within Disc().",_RESCONV_ACTIONNAME);
    if (!gl->DISC_RESTIME_READ) 
      fatal_error("The residual time module %s() was not found within Disc().",_RESTIME_ACTIONNAME);
  }
}


#ifdef DISTMPI
/* not used anymore */
void  MPI_Allreduce_Sum_Cliplist(char **cliplist_str){
  int rank,numproc,proc,thiscliplist_len;
  char *cliplistmem_str,*thiscliplist_str;
  cliplistmem_str=(char *)malloc((strlen(*cliplist_str)+10)*sizeof(char));
  strcpy(cliplistmem_str,*cliplist_str);
  thiscliplist_str=(char *)malloc(sizeof(char));
  strcpy(*cliplist_str,"");
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  for (proc=0; proc<numproc; proc++){
    if (proc==rank) {
      thiscliplist_len=strlen(cliplistmem_str);
      thiscliplist_str=(char *)realloc(thiscliplist_str,sizeof(char)*(thiscliplist_len+1));
      strcpy(thiscliplist_str,cliplistmem_str);
    }
    MPI_Bcast(&thiscliplist_len,1,MPI_INT,proc,MPI_COMM_WORLD);
    thiscliplist_str=(char *)realloc(thiscliplist_str,sizeof(char)*(thiscliplist_len+1));
    MPI_Bcast(thiscliplist_str,thiscliplist_len+1,MPI_CHAR,proc,MPI_COMM_WORLD);
    *cliplist_str=(char *)realloc(*cliplist_str,sizeof(char)*(strlen(*cliplist_str)+thiscliplist_len+1));
    strcat(*cliplist_str,thiscliplist_str);
  }
  free(cliplistmem_str);
  free(thiscliplist_str);
}


void find_clipped_variables_all(gl_t *gl){
  int rank,numproc,proc,cnt;
  int thisclipnamenum,thisclipname_len;
  char *thisclipname;
  long thisclipnum;
  reset_clipped_variables_all(gl);
  thisclipname=(char *)malloc(sizeof(char));
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  for (proc=0; proc<numproc; proc++){
    if (proc==rank) {
      thisclipnamenum=gl->model.clipnamenum;
    }
    MPI_Bcast(&thisclipnamenum,1,MPI_INT,proc,MPI_COMM_WORLD);
    for (cnt=0; cnt<thisclipnamenum; cnt++){
      if (proc==rank) {
        thisclipname_len=strlen(gl->model.clipname[cnt]);
      }
      MPI_Bcast(&thisclipname_len,1,MPI_INT,proc,MPI_COMM_WORLD);
      thisclipname=(char *)realloc(thisclipname,sizeof(char)*(thisclipname_len+1));
      if (proc==rank) {
        strcpy(thisclipname,gl->model.clipname[cnt]);
        thisclipnum=gl->model.clipnum[cnt];
      }
      MPI_Bcast(thisclipname,thisclipname_len+1,MPI_CHAR,proc,MPI_COMM_WORLD);
      MPI_Bcast(&thisclipnum,1,MPI_LONG,proc,MPI_COMM_WORLD);
      add_to_clipped_variables_all(gl, thisclipname, thisclipnum);
//      if (rank==0) printf("\n_%s(%ld)%d_",thisclipname,thisclipnum,proc);
    }
  }
  free(thisclipname);
}
#endif


void update_runtime_codex_xi_from_gl(gl_t *gl, SOAP_codex_t *codex){
  char *cliplist_str;
#ifdef DISTMPI
  int rank,proc;
  long ijk_ximax;
  struct {
    double ximax;
    int rank;
  } ximaxrank,ximaxrank_max;
#ifdef EMFIELD
  long ijk_ximax_emfield;
  struct {
    double ximax;
    int rank;
  } ximaxrank_emfield,ximaxrank_max_emfield;
#endif
#endif//DISTMPI

#ifdef DISTMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  if (rank!=0) codex->SCREENOUTPUT=FALSE;
#endif
  cliplist_str=(char *)malloc(sizeof(char));
#ifdef DISTMPI
  find_clipped_variables_all(gl);
  find_clipped_variables_list_all(gl,&cliplist_str);
  add_string_to_codex(codex,"clipinfo",cliplist_str);
  find_clipped_muscl_variables_list_all(gl,&cliplist_str);
  add_string_to_codex(codex,"clipinfo_muscl",cliplist_str);
  find_clipped_bdry_variables_list_all(gl,&cliplist_str);
  add_string_to_codex(codex,"clipinfo_bdry",cliplist_str);
#else
  find_clipped_variables_list(gl,&cliplist_str);
  add_string_to_codex(codex,"clipinfo",cliplist_str);
  find_clipped_muscl_variables_list(gl,&cliplist_str);
  add_string_to_codex(codex,"clipinfo_muscl",cliplist_str);
  find_clipped_bdry_variables_list(gl,&cliplist_str);
  add_string_to_codex(codex,"clipinfo_bdry",cliplist_str);
  //MPI_Allreduce_Sum_Cliplist(&cliplist_str);
#endif
  free(cliplist_str);
#ifdef DISTMPI
  ximaxrank.ximax=gl->ximax;
  ximaxrank.rank=rank;
  MPI_Allreduce(&ximaxrank, &ximaxrank_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  add_double_to_codex(codex,"ximax",ximaxrank_max.ximax);
  ijk_ximax=gl->flux_ximax+1;
  MPI_Bcast(&ijk_ximax,1,MPI_LONG,ximaxrank_max.rank,MPI_COMM_WORLD);
  add_int_to_codex(codex,"flux_ximax",ijk_ximax);
  ijk_ximax=gl->i_ximax;
  MPI_Bcast(&ijk_ximax,1,MPI_LONG,ximaxrank_max.rank,MPI_COMM_WORLD);
  add_int_to_codex(codex,"i_ximax",ijk_ximax);
#ifdef EMFIELD
  ximaxrank_emfield.ximax=gl->ximax_emfield;
  ximaxrank_emfield.rank=rank;
  MPI_Allreduce(&ximaxrank_emfield, &ximaxrank_max_emfield, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  add_double_to_codex(codex,"ximax_emfield",ximaxrank_max_emfield.ximax);
  ijk_ximax_emfield=gl->i_ximax_emfield;
  MPI_Bcast(&ijk_ximax_emfield,1,MPI_LONG,ximaxrank_max_emfield.rank,MPI_COMM_WORLD);
  add_int_to_codex(codex,"i_ximax_emfield",ijk_ximax_emfield);
#endif
#ifdef _2DL
  ijk_ximax=gl->j_ximax;
  MPI_Bcast(&ijk_ximax,1,MPI_LONG,ximaxrank_max.rank,MPI_COMM_WORLD);
  add_int_to_codex(codex,"j_ximax",ijk_ximax);
#ifdef EMFIELD
  ijk_ximax_emfield=gl->j_ximax_emfield;
  MPI_Bcast(&ijk_ximax_emfield,1,MPI_LONG,ximaxrank_max_emfield.rank,MPI_COMM_WORLD);
  add_int_to_codex(codex,"j_ximax_emfield",ijk_ximax_emfield);
#endif
#endif//_2DL
#ifdef _3DL
  ijk_ximax=gl->k_ximax;
  MPI_Bcast(&ijk_ximax,1,MPI_LONG,ximaxrank_max.rank,MPI_COMM_WORLD);
  add_int_to_codex(codex,"k_ximax",ijk_ximax);
#ifdef EMFIELD
  ijk_ximax_emfield=gl->k_ximax_emfield;
  MPI_Bcast(&ijk_ximax_emfield,1,MPI_LONG,ximaxrank_max_emfield.rank,MPI_COMM_WORLD);
  add_int_to_codex(codex,"k_ximax_emfield",ijk_ximax_emfield);
#endif
#endif//_3DL
#else//DISTMPI
  add_double_to_codex(codex,"ximax",gl->ximax);
  add_int_to_codex(codex,"flux_ximax",gl->flux_ximax+1);
  add_int_to_codex(codex,"i_ximax",gl->i_ximax);
#ifdef EMFIELD
  add_double_to_codex(codex,"ximax_emfield",gl->ximax_emfield);
  add_int_to_codex(codex,"i_ximax_emfield",gl->i_ximax_emfield);
#endif
#ifdef _2DL
  add_int_to_codex(codex,"j_ximax",gl->j_ximax);
#ifdef EMFIELD
  add_int_to_codex(codex,"j_ximax_emfield",gl->j_ximax_emfield);
#endif
#endif//_2DL
#ifdef _3DL
  add_int_to_codex(codex,"k_ximax",gl->k_ximax);
#ifdef EMFIELD
  add_int_to_codex(codex,"k_ximax_emfield",gl->k_ximax_emfield);
#endif
#endif//_3DL
#endif//DISTMPI
}


void update_runtime_codex_vars_except_xi_from_gl(gl_t *gl, SOAP_codex_t *codex){
#ifdef DISTMPI
  double effiter_U_sum,effiter_R_sum;
  int rank,proc;
#ifdef EMFIELD
  double effiter_U_sum_emfield,effiter_R_sum_emfield;
#endif
#endif//DISTMPI
#ifdef DISTMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  if (rank!=0) codex->SCREENOUTPUT=FALSE;
#endif
  add_int_to_codex(codex,"iter",   gl->iter);
  add_double_to_codex(codex,"xiverge",gl->cycle.fluid.xiverge);
  add_string_to_codex(codex,"outputfilename",   gl->output_filename);
#ifdef EMFIELD
  add_double_to_codex(codex,"xiverge_emfield",gl->cycle.emfield.xiverge);
#endif
#if defined(UNSTEADY)
  add_double_to_codex(codex,"time",gl->time);  
#endif
  add_double_to_codex(codex,"CFL",gl->CFL);
  add_int_to_codex(codex,"timelevelcnt",gl->timelevelcnt);
#ifdef UNSTEADY
  add_double_to_codex(codex,"dt",gl->dt);
#endif
#ifdef _CYCLE_MULTIZONE_MARCHING
  add_double_to_codex(codex,"window.is",gl->window.is);
  add_double_to_codex(codex,"window.ie",gl->window.ie);
  add_int_to_codex(&(gl->cycle.codex), "numzones_updated", 0);
  add_int_to_codex(&(gl->cycle.codex), "numzones_total", 0);
#endif
#ifdef _CYCLE_MULTIZONE
  add_int_to_codex(&(gl->cycle.codex), "numzones_updated", 0);
  add_int_to_codex(&(gl->cycle.codex), "numzones_total", 0);
#endif
#ifdef DISTMPI
  MPI_Allreduce(&gl->effiter_U, &effiter_U_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  add_double_to_codex(codex,"effiter_U",effiter_U_sum);
  MPI_Allreduce(&gl->effiter_R, &effiter_R_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  add_double_to_codex(codex,"effiter_R",effiter_R_sum);
#ifdef EMFIELD
  MPI_Allreduce(&gl->effiter_U_emfield, &effiter_U_sum_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  add_double_to_codex(codex,"effiter_U_emfield",effiter_U_sum_emfield);
  MPI_Allreduce(&gl->effiter_R_emfield, &effiter_R_sum_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  add_double_to_codex(codex,"effiter_R_emfield",effiter_R_sum_emfield);
#endif
#else//DISTMPI
  add_double_to_codex(codex,"effiter_U",gl->effiter_U);
  add_double_to_codex(codex,"effiter_R",gl->effiter_R);
#ifdef EMFIELD
  add_double_to_codex(codex,"Lc",gl->Lc);
//  add_double_to_codex(codex,"relaxEMF",gl->relaxEMF);
  add_double_to_codex(codex,"effiter_U_emfield",gl->effiter_U_emfield);
  add_double_to_codex(codex,"effiter_R_emfield",gl->effiter_R_emfield);
#endif
#endif//DISTMPI
#ifdef _AVERAGEDRATES
  add_int_to_codex(&(gl->cycle.codex), "AVERAGEDRATES_ON", AVERAGEDRATES_ON);
  add_int_to_codex(&(gl->cycle.codex), "AVERAGEDRATES_OFF", AVERAGEDRATES_OFF);
#endif
}


void add_constants_to_codex(gl_t *gl, SOAP_codex_t *codex){
  char str[100];
  sprintf(str, "%d", TSEMF_DEFAULT);
  SOAP_add_to_vars(codex,"TSEMF_DEFAULT",str);
  sprintf(str, "%d", TSEMF_ADI);
  SOAP_add_to_vars(codex,"TSEMF_ADI",str);
  sprintf(str, "%d", TSEMF_DDADI);
  SOAP_add_to_vars(codex,"TSEMF_DDADI",str);
  sprintf(str, "%d", TSEMF_IMAF);
  SOAP_add_to_vars(codex,"TSEMF_IMAF",str);
  sprintf(str, "%d", TSEMF_ADIIMAF);
  SOAP_add_to_vars(codex,"TSEMF_ADIIMAF",str);
  sprintf(str, "%d", TSEMF_NEWTON);
  SOAP_add_to_vars(codex,"TSEMF_NEWTON",str);
  sprintf(str, "%d", TSEMF_ADIi);
  SOAP_add_to_vars(codex,"TSEMF_ADIi",str);
  sprintf(str, "%d", TSEMF_ADIk);
  SOAP_add_to_vars(codex,"TSEMF_ADIk",str);
  sprintf(str, "%d", TSEMF_IMAFk);
  SOAP_add_to_vars(codex,"TSEMF_IMAFk",str);
  sprintf(str, "%d", TSEMF_IMAFi);
  SOAP_add_to_vars(codex,"TSEMF_IMAFi",str);
  sprintf(str, "%d", TSEMF_SOR);
  SOAP_add_to_vars(codex,"TSEMF_SOR",str);
  sprintf(str, "%d", TSEMF_SOR2);
  SOAP_add_to_vars(codex,"TSEMF_SOR2",str);

  sprintf(str, "%d", PRECON_CONSTANTTIMESTEP);
  SOAP_add_to_vars(codex,"PRECON_CONSTANTTIMESTEP",str);
  sprintf(str, "%d", PRECON_LOCALTIMESTEP);
  SOAP_add_to_vars(codex,"PRECON_LOCALTIMESTEP",str);
  sprintf(str, "%d", PRECON_LOCALTIMESTEP2);
  SOAP_add_to_vars(codex,"PRECON_LOCALTIMESTEP2",str);
  sprintf(str, "%d", PRECON_LOCALEIGENVALUE);
  SOAP_add_to_vars(codex,"PRECON_LOCALEIGENVALUE",str);
  sprintf(str, "%d", PRECON_LOCALEIGENVALUE2);
  SOAP_add_to_vars(codex,"PRECON_LOCALEIGENVALUE2",str);

}


void process_code_runtime(np_t *np, gl_t *gl, char *code_runtime, SOAP_codex_t *codex){
  char *code;
  SOAP_vars_t *varsmem;
  readcontrolarg_t Runtimearg;

  varsmem=(SOAP_vars_t *)malloc(sizeof(SOAP_vars_t));
  SOAP_copy_all_vars(codex->vars, &varsmem);
  Runtimearg.np=&np;
  Runtimearg.gl=gl;
  Runtimearg.input=(input_t *)malloc(sizeof(input_t));
  Runtimearg.input->READDATAFILE=FALSE;
  Runtimearg.TYPELEVEL=TYPELEVEL_FLUID;
  Runtimearg.module_level=0;
  Runtimearg.POSTMODULE=FALSE;
  Runtimearg.CYCLEMODULE=FALSE;
  Runtimearg.RESETITERCOUNT=FALSE;
  Runtimearg.CFLinit=-1.0; //setting CFLinit to an invalid value will not change the CFL
  Runtimearg.VERBOSE=FALSE;
  Runtimearg.gl_post=*gl;
  Runtimearg.domain_post=gl->domain;
  Runtimearg.np_post=np;
  if (!gl->cycle.RUNTIMEMODULEFOUND)
    fatal_error("The %s() module was not found within Cycle().",_CYCLE_ACTIONNAME);
  code=(char *)malloc((strlen(code_runtime)+2)*sizeof(char));
  strcpy(code,code_runtime);
  codex->ACTION=TRUE;
  codex->action=&runtime_actions;
  codex->action_args=(void *)&Runtimearg;
  ((readcontrolarg_t *)codex->action_args)->np=&np;
  ((readcontrolarg_t *)codex->action_args)->gl=gl;
/*  if (codex->action_being_processed==NULL){
    codex->action_being_processed=(char *)malloc((strlen(_CYCLE_ACTIONNAME)+2)*sizeof(char));
    strcpy(codex->action_being_processed,_CYCLE_ACTIONNAME);
  }*/
  codex->VERBOSE=FALSE;
  codex->SCREENOUTPUT=TRUE;
  add_constants_to_codex(gl, codex);
  update_runtime_codex_xi_from_gl(gl, codex);
  update_runtime_codex_vars_except_xi_from_gl(gl,codex);
  SOAP_process_code(code, codex, SOAP_VARS_KEEP_ALL);
  gl->CFL=SOAP_var_value(codex,"CFL");
#ifdef UNSTEADY
  gl->dt=SOAP_var_value(codex,"dt");
#endif
  gl->ximax=SOAP_var_value(codex,"ximax");
  assert(gl->CFL>=0.0e0);
  /* here, make sure that all changes to vars within runtime module are erased, because those
     will not be written to datafile 
     -> CFL and ximax and dt are exception to this, and this is why they are probed 
        through SOAP_var_value above */  
  if (gl->RESETRUNTIMEVARS){
    SOAP_free_all_vars(codex->vars);
    SOAP_copy_all_vars(varsmem,&(codex->vars));
  }
  free(Runtimearg.input);
  SOAP_free_all_vars(varsmem);
  free(varsmem);
  free(code);
  reset_clipped_variables(gl);
}


void find_ximax(np_t *np, gl_t *gl, zone_t zone, int IJK_UPDATE){
  long i,j,k;
  double xi;
  gl->ximax=0.0e0;
  for_ijk(zone,is,js,ks,ie,je,ke){
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)) {
          assert(is_node_resumed(np[_ai(gl,i,j,k)]));
          xi=np[_ai(gl,i,j,k)].wk->xi;
          if (xi<-1.0e99 || isnan(xi)) {
            fatal_error("problem with xi (xi=%E) at i=%ld, j=%ld, k=%ld.",xi,i,j,k);
          }
          if (xi>=gl->ximax) {
            gl->ximax=xi;
            if (IJK_UPDATE==IJK_UPDATE_YES) {
              gl->flux_ximax=np[_ai(gl,i,j,k)].wk->flux_xi;
              gl->i_ximax=i;
              gl->j_ximax=j;
              gl->k_ximax=k;
            }
          }
        }
  }
}



#if defined(UNSTEADY) && defined(_AVERAGEDRATES)
void set_averaged_rates_to_zero(np_t *np, gl_t *gl){
  long i,j,k,cnt;
  gl->averagedrates_time=0.0;
  for_ijk(gl->domain,is,js,ks,ie,je,ke){
    if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)) {
      assert(is_node_resumed(np[_ai(gl,i,j,k)]));
      for (cnt=0; cnt<numaveragedrates; cnt++) np[_ai(gl,i,j,k)].bs->averagedrates[cnt]=0.0;
    }
  }
}


void add_to_averaged_rates(np_t *np, gl_t *gl){
  long i,j,k;
  flux_t Sstar;
  gl->AVERAGEDRATES=AVERAGEDRATES_ADD;
  for_ijk(gl->domain,is,js,ks,ie,je,ke){
    if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)) {
      assert(is_node_resumed(np[_ai(gl,i,j,k)]));
      //find_Schem(np,gl,_ai(gl,i,j,k),Sstar);
      find_Sstar(np,gl,_ai(gl,i,j,k),Sstar);
    }
  }
  gl->averagedrates_time+=gl->dt;
}

#endif

/*
static void PrintZones(zone_t *zones, long numzone){
  long cnt;

  for (cnt=0; cnt<numzone; cnt++){
    printf("%ld   is=%ld js=%ld    ie=%ld je=%ld\n",cnt,zones[cnt].is,zones[cnt].js,
    zones[cnt].ie,zones[cnt].je);
  }
  printf("\n");
}
*/

static void rearrange_overlapping_zones(zone_t *zones, long numzone){
  long cnt1,cnt2;

  /* PrintZone(zones,numzones); */
  for (cnt1=0; cnt1<numzone; cnt1++){
    for (cnt2=0; cnt2<numzone; cnt2++){
      if (cnt2!=cnt1){
        /* do overlap along i :
           make ie of zones[cnt1] smaller and is of zones[cnt2] bigger */
        if (if3DL( zones[cnt1].ks==zones[cnt2].ks && )
            if2DL( zones[cnt1].js==zones[cnt2].js && )
            if3DL( zones[cnt1].ke==zones[cnt2].ke && )
            if2DL( zones[cnt1].je==zones[cnt2].je && )
                   zones[cnt1].ie< zones[cnt2].ie &&
                   zones[cnt1].ie>=zones[cnt2].is) {
          zones[cnt1].ie=(zones[cnt1].ie+zones[cnt2].is)/2;
          zones[cnt2].is=zones[cnt1].ie+1;
          if (    zones[cnt1].is>zones[cnt1].ie
               || zones[cnt2].is>zones[cnt2].ie )
             fatal_error("Problem modifying zones along i.");
        }
      }
    }
  }
#ifdef _2DL
  for (cnt1=0; cnt1<numzone; cnt1++){
    for (cnt2=0; cnt2<numzone; cnt2++){
      if (cnt2!=cnt1){
        /* do overlap along j :
           make je of zones[cnt1] smaller and js of zones[cnt2] bigger*/
        if (if3DL( zones[cnt1].ks==zones[cnt2].ks && )
                   zones[cnt1].is==zones[cnt2].is &&
            if3DL( zones[cnt1].ke==zones[cnt2].ke && )
                   zones[cnt1].ie==zones[cnt2].ie &&
                   zones[cnt1].je< zones[cnt2].je &&
                   zones[cnt1].je>=zones[cnt2].js) {
          zones[cnt1].je=(zones[cnt1].je+zones[cnt2].js)/2;
          zones[cnt2].js=zones[cnt1].je+1;
          if (    zones[cnt1].js>zones[cnt1].je
               || zones[cnt2].js>zones[cnt2].je )
             fatal_error("Problem modifying zones along j.");
        }
      }
    }
  }
#endif
#ifdef _3DL
  for (cnt1=0; cnt1<numzone; cnt1++){
    for (cnt2=0; cnt2<numzone; cnt2++){
      if (cnt2!=cnt1){
        /* do overlap along k :
           make je of zones[cnt1] smaller and js of zones[cnt2] bigger*/
        if (zones[cnt1].is==zones[cnt2].is &&
            zones[cnt1].js==zones[cnt2].js &&
            zones[cnt1].ie==zones[cnt2].ie &&
            zones[cnt1].je==zones[cnt2].je &&
            zones[cnt1].ke< zones[cnt2].ke &&
            zones[cnt1].ke>=zones[cnt2].ks) {
          zones[cnt1].ke=(zones[cnt1].ke+zones[cnt2].ks)/2;
          zones[cnt2].ks=zones[cnt1].ke+1;
          if (    zones[cnt1].ks>zones[cnt1].ke
               || zones[cnt2].ks>zones[cnt2].ke )
             fatal_error("Problem modifying zones along k.");
        }
      }
    }
  }
#endif
  /* PrintZone(zone,numzone); */
}


/* setup multizone situated inside zone */
void setup_multizone(np_t *np, gl_t *gl, zone_t zone, zone_t lim, double xiverge,
                    long zonelength, bool UPDATE_ALL_ZONES, multizone_t *multizone){
  long cnt;
  long numsubzones;
  zone_t *subzones;
  double ximax;
  long i,j,k;
  /* find the zones for the ts process */
  subzones=(zone_t *)malloc(sizeof(zone_t));
  find_subzones_in_zone_given_zonelength(zonelength, zone, &numsubzones, &subzones);
  /* find out which zones need to be updated */
  multizone->numzones_ts=0;
  multizone->ts=(zone_t *)malloc(numsubzones*sizeof(zone_t));
  for (cnt=0; cnt<numsubzones; cnt++){
    ximax=0.0e0;
    for_ijk(subzones[cnt],is,js,ks,ie,je,ke){
          if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)) {
            ximax=max(ximax,np[_ai(gl,i,j,k)].wk->xi);
          }
    }
    if (ximax>xiverge || UPDATE_ALL_ZONES) {
      multizone->ts[multizone->numzones_ts]=subzones[cnt];
      (multizone->numzones_ts)++;
    }
  }
  /* setup res and bdry, limited by lim_is,lim_js, etc*/
  multizone->bdry=(zone_t *)malloc(multizone->numzones_ts*sizeof(zone_t));
  multizone->res=(zone_t *)malloc(multizone->numzones_ts*sizeof(zone_t));
  for (cnt=0; cnt<multizone->numzones_ts; cnt++){
    multizone->bdry[cnt].is=max(lim.is,multizone->ts[cnt].is-hbw_bdry_fluid);
    multizone->bdry[cnt].ie=min(lim.ie,multizone->ts[cnt].ie+hbw_bdry_fluid);
#ifdef _2DL
    multizone->bdry[cnt].js=max(lim.js,multizone->ts[cnt].js-hbw_bdry_fluid);
    multizone->bdry[cnt].je=min(lim.je,multizone->ts[cnt].je+hbw_bdry_fluid);
#else
    multizone->bdry[cnt].js=0;
    multizone->bdry[cnt].je=0;
#endif
#ifdef _3DL
    multizone->bdry[cnt].ks=max(lim.ks,multizone->ts[cnt].ks-hbw_bdry_fluid);
    multizone->bdry[cnt].ke=min(lim.ke,multizone->ts[cnt].ke+hbw_bdry_fluid);
#else
    multizone->bdry[cnt].ks=0;
    multizone->bdry[cnt].ke=0;
#endif
    multizone->res[cnt].is=max(lim.is,multizone->ts[cnt].is-hbw_bdry_fluid-hbw_res_fluid);
    multizone->res[cnt].ie=min(lim.ie,multizone->ts[cnt].ie+hbw_bdry_fluid+hbw_res_fluid);
#ifdef _2DL
    multizone->res[cnt].js=max(lim.js,multizone->ts[cnt].js-hbw_bdry_fluid-hbw_res_fluid);
    multizone->res[cnt].je=min(lim.je,multizone->ts[cnt].je+hbw_bdry_fluid+hbw_res_fluid);
#else
    multizone->res[cnt].js=0;
    multizone->res[cnt].je=0;
#endif
#ifdef _3DL
    multizone->res[cnt].ks=max(lim.ks,multizone->ts[cnt].ks-hbw_bdry_fluid-hbw_res_fluid);
    multizone->res[cnt].ke=min(lim.ke,multizone->ts[cnt].ke+hbw_bdry_fluid+hbw_res_fluid);
#else
    multizone->res[cnt].ks=0;
    multizone->res[cnt].ke=0;
#endif
  }
  multizone->numzones_total=numsubzones;
  multizone->numzones_res=multizone->numzones_ts;
  multizone->numzones_bdry=multizone->numzones_ts;
  free(subzones);
  rearrange_overlapping_zones(multizone->res,multizone->numzones_res);
}


void *thread_zone(void *threadzone){
  np_t * np = ((threadzone_t *) threadzone)->np;
  gl_t * gl = ((threadzone_t *) threadzone)->gl;
  zone_t zone = ((threadzone_t *) threadzone)->zone;
  ((threadzone_t *) threadzone)->funct(np,gl,zone);
  return(NULL);
}


void create_thread_zone(np_t *np, gl_t * gl, zone_t zone, void (*funct)(np_t *, gl_t *, zone_t zone),
                      pthread_t *pthread, threadzone_t *threadzone){
  threadzone->np=np;
  threadzone->gl=gl;
  threadzone->zone=zone;
  threadzone->funct=funct;
#ifdef ZONETHREADS
  if (pthread_create(pthread, NULL, thread_zone, threadzone))
    fatal_error("Cannot create thread.");
#else
  (*thread_zone)(threadzone);
#endif
}


void join_all_threads_zone(long numthread, pthread_t *pthread, bool COUNTFLAG){
#ifdef ZONETHREADS
  long thread;
  void *retval;
  if (!COUNTFLAG) {
    for (thread=0; thread<numthread; thread++){
      if (pthread_join(pthread[thread],&retval))
        fatal_error("Cannot join thread %ld.",thread);
    }
  }
#endif
}


static void update_U_from_dUstar_1(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    thread_lock_node_set(np,l,THREADTYPE_ZONE);
    add_dUstar_to_U(np,l,gl,np[l].wk->dUstar);
    thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    /* - if not using SMALLTHREADS, only need to lock for the loop threads, 
         since gl is local for the zone thread 
       - if using SMALLTHREADS, then need to lock for both the loop and zone threads
       For now, lock for both the loop and zone threads */
    thread_lock_global_set(gl,THREADTYPE_ALL);
    gl->effiter_U+=1.0/(double)(gl->nn);
    thread_lock_global_unset(gl,THREADTYPE_ALL);
  }
}


static void update_U_from_dUstar(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np,gl,zone,&update_U_from_dUstar_1,SWEEPTYPE_I,TYPELEVEL_FLUID_WORK,
                         &is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);
}


long _numthread_optimized(long numzone){
  long l,cnt,lmax,numthread;
  numthread=numzone;
  if (numzone>maxzonethread) {
    lmax=0;
    for (cnt=1; cnt<=maxzonethread; cnt++){
      l=mod(numzone,cnt);
      if (l==0) l=cnt;
      if (l>lmax) {
        numthread=cnt;
        lmax=l;
      }
    }
  }
  return(numthread);
}


#ifdef DISTMPI


void exchange_U(np_t *np, gl_t *gl){
  int bl,rankrecv,numproc,ranksend,thisrank,pack_size_Ulocal,pack_size_cnt;
  long i,j,k,flux,iterator,cnt,prevcnt,total;
  long primcnt=0;
  long bufsize=0;
  long *recvcnt,*sendcnt,*processcnt;
  long *primnodenums=NULL;
  long *processnodenums;
  double *buf,*bufptr;
  zone_t zonesend,zonerecv,zone;
  flux_t *recvUlocal;
  flux_t *sendUlocal=NULL;
  MPI_Status MPI_Status1;
 
  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  recvcnt=(long *)malloc(numproc*sizeof(long));
  sendcnt=(long *)malloc(numproc*sizeof(long));
  processcnt=(long *)malloc(numproc*sizeof(long));
  processnodenums=(long *)malloc(numproc*sizeof(long));
  for (i=0; i<numproc; i++){
    sendcnt[i]=0;
    processcnt[i]=0;
  }
  for (ranksend=0; ranksend<numproc; ranksend++){
    zonesend=_domain_from_rank(ranksend,gl);
    for (rankrecv=0; rankrecv<numproc; rankrecv++){
      if (rankrecv!=ranksend && (ranksend==thisrank || rankrecv==thisrank)){
        zonerecv=_domain_lim_from_rank(rankrecv,gl);
        if (is_zone_intersecting_zone(zonesend,zonerecv)){
          zone=_zone_intersection(zonesend,zonerecv);
          for_ijk(zone,is,js,ks,ie,je,ke){
                if (ranksend==thisrank) {
                  for (total=0,iterator=0; iterator<numproc; iterator++) total+=sendcnt[iterator];
                  for (prevcnt=0,iterator=0; iterator<=rankrecv; iterator++) prevcnt+=sendcnt[iterator];
                  sendUlocal=(flux_t *)realloc(sendUlocal,(total+1)*sizeof(flux_t));
                  for (iterator=prevcnt+1;iterator<total+1;iterator++){
                    for (flux=0; flux<nf; flux++) *(*(sendUlocal + iterator) + flux)=*(*(sendUlocal + iterator-1) + flux);
                  }
                  for (flux=0; flux<nf; flux++) *(*(sendUlocal + prevcnt) + flux)=np[_ai(gl,i,j,k)].bs->U[flux];
                  sendcnt[rankrecv]++;
                }
                if (rankrecv==thisrank){
                  for (prevcnt=0,iterator=0; iterator<=ranksend; iterator++) prevcnt+=processcnt[iterator];
                  processnodenums=(long *)realloc(processnodenums,(prevcnt+1)*sizeof(long));
                  processnodenums[prevcnt]=_ai(gl,i,j,k);
                  processcnt[ranksend]++;
                  if (is_node_resumed(np[_ai(gl,i,j,k)])){
                    primnodenums=(long *)realloc(primnodenums,(primcnt+1)*sizeof(long));
                    primnodenums[primcnt]=_ai(gl,i,j,k);
                    primcnt++;
                  }
                }
          }
        }
      }
    }
  }
 
  if(numproc != 1){
    for (rankrecv=0; rankrecv<numproc; rankrecv++){
      if (thisrank!=rankrecv){
        MPI_Pack_size(nf*sendcnt[rankrecv],MPI_DOUBLE,MPI_COMM_WORLD,&pack_size_Ulocal);
        MPI_Pack_size(1,MPI_LONG,MPI_COMM_WORLD,&pack_size_cnt);
        bufsize+=(2*MPI_BSEND_OVERHEAD)+pack_size_Ulocal+pack_size_cnt;
      }
    }
    buf=(double *)malloc(bufsize);
    MPI_Buffer_attach(buf, bufsize);
    for (rankrecv=0; rankrecv<numproc; rankrecv++){
      if (thisrank!=rankrecv){
        for (prevcnt=0,iterator=0; iterator<rankrecv; iterator++) prevcnt+=sendcnt[iterator];
        MPI_Bsend(&sendcnt[rankrecv],1,MPI_LONG,rankrecv,1,MPI_COMM_WORLD);
        MPI_Bsend(&sendUlocal[prevcnt],nf*sendcnt[rankrecv],MPI_DOUBLE,rankrecv,0,MPI_COMM_WORLD);
      }
    }
    free(sendUlocal);
    for (ranksend=0; ranksend<numproc; ranksend++){
      if (thisrank!=ranksend){
        MPI_Recv(&recvcnt[ranksend],1,MPI_LONG,ranksend,1,MPI_COMM_WORLD,&MPI_Status1);
        recvUlocal=(flux_t *)malloc(recvcnt[ranksend]*sizeof(flux_t));
        MPI_Recv(recvUlocal,recvcnt[ranksend]*nf,MPI_DOUBLE,ranksend,0,MPI_COMM_WORLD,&MPI_Status1);
        for (cnt=0; cnt<recvcnt[ranksend]; cnt++){
          for (prevcnt=0,iterator=0; iterator<ranksend; iterator++) prevcnt+=processcnt[iterator];
          for (flux=0; flux<nf; flux++) np[processnodenums[prevcnt+cnt]].bs->U[flux]=*(*(recvUlocal + cnt) + flux);
        }
        free(recvUlocal);
      }
    }

#ifdef OPENMPTHREADS
#pragma omp parallel for private(cnt) schedule(dynamic)
#endif
    for (cnt=0; cnt<primcnt; cnt++) find_prim_fluid(np,primnodenums[cnt],gl);
 
    MPI_Buffer_detach(&bufptr,&bl);
    free(buf);
  }
  
  free(processnodenums);
  free(primnodenums);
  free(processcnt);
  free(recvcnt);
  free(sendcnt);
  MPI_Barrier(MPI_COMM_WORLD);
}


void exchange_U_old(np_t *np, gl_t *gl){  //same as above but without the MPI_Buffer
  int rankrecv,numproc,ranksend,thisrank;
  long i,j,k,flux;
  long cnt = 0;
  long *nodenums = NULL;
  zone_t zonesend,zonerecv,zone;
  flux_t Ulocal;
  MPI_Status MPI_Status1;
 
  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  for (ranksend=0; ranksend<numproc; ranksend++){
    zonesend=_domain_from_rank(ranksend,gl);
    for (rankrecv=0; rankrecv<numproc; rankrecv++){
      if (rankrecv!=ranksend && (ranksend==thisrank || rankrecv==thisrank)){
        zonerecv=_domain_lim_from_rank(rankrecv,gl);
        if (is_zone_intersecting_zone(zonesend,zonerecv)){
          zone=_zone_intersection(zonesend,zonerecv);
          for_ijk(zone,is,js,ks,ie,je,ke){
                if (ranksend==thisrank) {
                  for (flux=0; flux<nf; flux++) Ulocal[flux]=np[_ai(gl,i,j,k)].bs->U[flux];
                  MPI_Send(Ulocal,nf,MPI_DOUBLE,rankrecv,0,MPI_COMM_WORLD);
                }
                if (rankrecv==thisrank) {
                  MPI_Recv(Ulocal,nf,MPI_DOUBLE,ranksend,0,MPI_COMM_WORLD,&MPI_Status1);
                  for (flux=0; flux<nf; flux++) np[_ai(gl,i,j,k)].bs->U[flux]=Ulocal[flux];
                  if (is_node_resumed(np[_ai(gl,i,j,k)])){
                    nodenums=(long *)realloc(nodenums,(cnt+1)*sizeof(long));
                    nodenums[cnt]=_ai(gl,i,j,k);
                    cnt++;
                  }
                }
          }
        }
      }
    }
  }
 
#ifdef OPENMPTHREADS
#pragma omp parallel for private(i) schedule(dynamic)
#endif
  for (i=0; i<cnt; i++) find_prim_fluid(np,nodenums[i],gl);
  free(nodenums);
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif


void update_U_with_multizone(np_t *np, gl_t *gl, multizone_t multizone){
  long cnt,numzonethread,cntthread;
  pthread_t *pthread;
  threadzone_t *threadzone;

  /* Find dUstar for inner nodes*/
  numzonethread=_numthread_optimized(multizone.numzones_ts);
  pthread=(pthread_t *)malloc(numzonethread*sizeof(pthread_t));
  threadzone=(threadzone_t *)malloc(numzonethread*sizeof(threadzone_t));
  cntthread=0;
  for (cnt=0; cnt<multizone.numzones_ts; cnt++) {
    create_thread_zone(np, gl, multizone.ts[cnt], &find_dU, &(pthread[cntthread]), &(threadzone[cntthread]));
    cntthread++;
    if (cntthread==numzonethread) {
      join_all_threads_zone(cntthread, pthread, FALSE);
      cntthread=0;
    }
  }
  if (cntthread>0) join_all_threads_zone(cntthread, pthread, FALSE);
  for (cnt=0; cnt<multizone.numzones_ts; cnt++) update_U_from_dUstar(np, gl, multizone.ts[cnt]);
  free(pthread);
  free(threadzone);
}


void update_bdry_nodes_with_multizone(np_t *np, gl_t *gl, multizone_t multizone){
  long cnt;
  for (cnt=0; cnt<multizone.numzones_bdry; cnt++) update_bdry_nodes(np, gl, multizone.bdry[cnt]);
}


void find_residual_with_multizone(np_t *np, gl_t *gl, multizone_t multizone){
  long cnt,numzonethread,cntthread;
  pthread_t *pthread;
  threadzone_t *threadzone;
  numzonethread=_numthread_optimized(multizone.numzones_res);
  pthread=(pthread_t *)malloc(numzonethread*sizeof(pthread_t));
  threadzone=(threadzone_t *)malloc(numzonethread*sizeof(threadzone_t));
  cntthread=0;
  for (cnt=0; cnt<multizone.numzones_res; cnt++) {
    create_thread_zone(np, gl, multizone.res[cnt], &find_residual, &(pthread[cntthread]), &(threadzone[cntthread]));
    cntthread++;
    if (cntthread==numzonethread) {
      join_all_threads_zone(cntthread, pthread, FALSE);
      cntthread=0;
    }
  }
  if (cntthread>0) {
    join_all_threads_zone(cntthread, pthread, FALSE);
  }
  free(pthread);
  free(threadzone);
}


void solve_multizone(np_t *np, gl_t *gl, multizone_t multizone){
  update_U_with_multizone(np,gl,multizone);
  update_bdry_nodes_with_multizone(np,gl,multizone);
  find_residual_with_multizone(np,gl,multizone);
}


void free_multizone(multizone_t *multizone){
  free(multizone->res);
  free(multizone->bdry);
  free(multizone->ts);
}


void check_residual(np_t *np, gl_t *gl, zone_t zone){
  resume_nodes_in_zone(np, gl, zone);
#ifdef EMFIELD
  update_prim_emfield_mem_in_zone(np, gl, zone);
#endif
  find_residual(np, gl, zone);
  find_ximax(np,gl,zone,IJK_UPDATE_YES);
#ifdef EMFIELD
  find_residual_emfield(np,gl,zone);
  find_ximax_emfield(np, gl, zone);
#endif
#ifdef DISTMPI
  int rank,proc;
  struct {
    double ximax;
    int rank;
  } ximaxrank,ximaxrank_max;
#ifdef EMFIELD
  struct {
    double ximax;
    int rank;
  } ximaxrank_emfield,ximaxrank_max_emfield;
#endif
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  ximaxrank.ximax=gl->ximax;
  ximaxrank.rank=rank;
  MPI_Allreduce(&ximaxrank, &ximaxrank_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  gl->ximax=ximaxrank_max.ximax;
  MPI_Bcast(&(gl->i_ximax),1,MPI_LONG,ximaxrank_max.rank,MPI_COMM_WORLD);
#ifdef EMFIELD
  ximaxrank_emfield.ximax=gl->ximax_emfield;
  ximaxrank_emfield.rank=rank;
  MPI_Allreduce(&ximaxrank_emfield, &ximaxrank_max_emfield, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  gl->ximax_emfield=ximaxrank_max_emfield.ximax;
  MPI_Bcast(&(gl->i_ximax_emfield),1,MPI_LONG,ximaxrank_max_emfield.rank,MPI_COMM_WORLD);
#endif
#ifdef _2DL
  MPI_Bcast(&(gl->j_ximax),1,MPI_LONG,ximaxrank_max.rank,MPI_COMM_WORLD);
#ifdef EMFIELD
  MPI_Bcast(&(gl->j_ximax_emfield),1,MPI_LONG,ximaxrank_max_emfield.rank,MPI_COMM_WORLD);
#endif
#endif //_2DL
#ifdef _3DL
  MPI_Bcast(&(gl->k_ximax),1,MPI_LONG,ximaxrank_max.rank,MPI_COMM_WORLD);
#ifdef EMFIELD
  MPI_Bcast(&(gl->k_ximax_emfield),1,MPI_LONG,ximaxrank_max_emfield.rank,MPI_COMM_WORLD);
#endif
#endif //_3DL
#endif //DISTMPI
}


void find_xi(np_t np, gl_t *gl, flux_t Res, double *xi, long *flux_xi){
  long flux;
  double xitmp;
  assert_np(np,is_node_resumed(np));
  *xi=0.0;
  *flux_xi=-1;
  for (flux=0; flux<nf; flux++) {
    xitmp=fabs(Res[flux]/_Omega(np,gl)/gl->cycle.fluid.Uref[flux]);
    if (xitmp>*xi) {
      *xi=max(*xi,xitmp);
      *flux_xi=flux;
    }
    if (isnan(xitmp)){
      fatal_error("problem computing xitmp in function find_xi() in cycle_share.c;\n   xitmp=%E\n   Res[%ld]=%E\n   Omega=%E\n   Uref[%ld]=%E\n",xitmp,flux,Res[flux],_Omega(np,gl),flux,gl->cycle.fluid.Uref[flux]);
    } 
  }
}


static void find_Delta_Lambda_for_dtau_local(np_t *np, gl_t *gl, long l, long dim, flux_t Delta_Lambda){
  long offset,maxoffset,flux,dim2;
  flux_t Delta_Lambda_tmp;
  
  find_Delta_Lambda_for_dtau(np, gl, l, dim, Delta_Lambda);
  
  if (gl->PRECONDITIONER==PRECON_LOCALTIMESTEP2){
    maxoffset=1;
    for (dim2=dim; dim2<=dim; dim2++){
      for (offset=1; offset<=maxoffset; offset++){
        if (is_node_inner(np[_al(gl,l,dim2,-offset)],TYPELEVEL_FLUID_WORK)){
          find_Delta_Lambda_for_dtau(np, gl, _al(gl,l,dim2,-offset), dim, Delta_Lambda_tmp);
          for (flux=0; flux<nf; flux++)
            Delta_Lambda[flux]=max(Delta_Lambda[flux],Delta_Lambda_tmp[flux]);
        }
        if (is_node_inner(np[_al(gl,l,dim2,+offset)],TYPELEVEL_FLUID_WORK)){
          find_Delta_Lambda_for_dtau(np, gl, _al(gl,l,dim2,+offset), dim, Delta_Lambda_tmp);
          for (flux=0; flux<nf; flux++)
            Delta_Lambda[flux]=max(Delta_Lambda[flux],Delta_Lambda_tmp[flux]);
        }
      }
    }
  }
}


void find_dtau(np_t *np, gl_t *gl, long l, flux_t dtau){
  double dtaumin,dtaumax;
  long dim,flux;
  double dtaulocal[nf][nd];
#ifdef UNSTEADY
  sqmat_t LambdaZ;
#endif
  flux_t Delta_Lambda;
  assert_np(np[l],is_node_inner(np[l],TYPELEVEL_FLUID_WORK));
  if (gl->PRECONDITIONER!=PRECON_CONSTANTTIMESTEP){
#ifdef UNSTEADY
    find_LambdaZ(np,gl,l,LambdaZ);
    set_matrix_to_identity(LambdaZ); //turn off effect of LambdaZ -> seems to be detrimental not beneficial
    for (dim=0; dim<nd; dim++){
      find_Delta_Lambda_for_dtau_local(np, gl, l, dim, Delta_Lambda);
      for (flux=0; flux<nf; flux++){
        assert(LambdaZ[flux][flux]>0.0);
        dtaulocal[flux][dim]=gl->dt/LambdaZ[flux][flux]/notzero(Delta_Lambda[flux]*gl->dt/LambdaZ[flux][flux]+1.0,1e-39);
      }
    }
#else
    for (dim=0; dim<nd; dim++){
      find_Delta_Lambda_for_dtau_local(np, gl, l, dim, Delta_Lambda);
      for (flux=0; flux<nf; flux++){
        dtaulocal[flux][dim]=1.0/notzero(Delta_Lambda[flux],1e-39);
      }
    }
#endif
   /* find optimal dtaus for each flux */
   for (flux=0; flux<nf; flux++){
    dtaumin=1.0e99;
    dtaumax=0.0e0;
    for (dim=0; dim<nd; dim++){
      dtaumin=min(dtaulocal[flux][dim],dtaumin);
      dtaumax=max(dtaulocal[flux][dim],dtaumax);
    }
    dtaumax=min(dtaumin*MAXRATIO_DTAUMAX_DTAUMIN,dtaumax);
    dtau[flux]=gl->CFL*pow(dtaumin,1.0e0-gl->sigma1)*pow(dtaumax,gl->sigma1);
   }
  } else {
    for (flux=0; flux<nf; flux++){    
      dtau[flux]=gl->dtau;
    }
  }
}


void find_constant_dtau(np_t *np, gl_t *gl, long l, double *dtau){
  long flux;
  flux_t dtau_vector;
  double dtaumin,dtaumax;

  find_dtau(np,gl,l,dtau_vector);
  /* average min and max dtau */
  dtaumin=1.0e99;
  dtaumax=-1.0e99;
  for (flux=0; flux<nf; flux++) dtaumin=min(dtaumin,dtau_vector[flux]);
  for (flux=0; flux<nf; flux++) dtaumax=max(dtaumax,dtau_vector[flux]);
  dtaumax=min(dtaumin*MAXRATIO_DTAUMAX_DTAUMIN,dtaumax);
  *dtau=pow(dtaumin,1.0-gl->sigma2)*pow(dtaumax,gl->sigma2);
}


#ifdef EMFIELD
#ifdef DISTMPI

void exchange_U_emfield(np_t *np, gl_t *gl){
  int rankrecv,numproc,ranksend,thisrank;
  long i,j,k,flux;
  zone_t zonesend,zonerecv,zone;
  fluxemfield_t Ulocal;
  MPI_Status MPI_Status1;

  MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  
  for (ranksend=0; ranksend<numproc; ranksend++){
    zonesend=_domain_from_rank(ranksend,gl);
    for (rankrecv=0; rankrecv<numproc; rankrecv++){
      if (rankrecv!=ranksend && (ranksend==thisrank || rankrecv==thisrank)){
        zonerecv=_domain_lim_from_rank(rankrecv,gl);
        if (is_zone_intersecting_zone(zonesend,zonerecv)){
          zone=_zone_intersection(zonesend,zonerecv);
          for_ijk(zone,is,js,ks,ie,je,ke){
                if (ranksend==thisrank) {
                  for (flux=0; flux<nfe; flux++) Ulocal[flux]=np[_ai(gl,i,j,k)].bs->Uemfield[flux];
                  MPI_Send(Ulocal,nfe,MPI_DOUBLE,rankrecv,0,MPI_COMM_WORLD);
                }
                if (rankrecv==thisrank) {
                  MPI_Recv(Ulocal,nfe,MPI_DOUBLE,ranksend,0,MPI_COMM_WORLD,&MPI_Status1);
                  for (flux=0; flux<nfe; flux++) np[_ai(gl,i,j,k)].bs->Uemfield[flux]=Ulocal[flux];
                }
          }
        }
      }
    }
  } 
   
  MPI_Barrier(MPI_COMM_WORLD);
}


void exchange_U_emfield_old(np_t *np, gl_t *gl){
  int rank;
  long i,j,k,flux;
  fluxemfield_t Ulocal;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for_ijk (gl->domain_all,is,js,ks,ie,je,ke){
        if (rank==_node_rank(gl, i, j, k) && is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
          for (flux=0; flux<nfe; flux++) {
            Ulocal[flux]=np[_ai(gl,i,j,k)].bs->Uemfield[flux];
          }
        }
        MPI_Bcast_Node(Ulocal,nfe,MPI_DOUBLE,_node_rank(gl,i,j,k),MPI_COMM_WORLD,i,j,k,gl);
        if (is_node_in_zone(i,j,k,gl->domain_lim) && is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
          for (flux=0; flux<nfe; flux++) {
            np[_ai(gl,i,j,k)].bs->Uemfield[flux]=Ulocal[flux];
          }
        }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif


void update_prim_emfield_mem_in_zone_1(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;
  //printf("(%ld,%ld) to (%ld,%ld)\n",_i(ls,gl,0),_i(ls,gl,1),_i(le,gl,0),_i(le,gl,1));
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (is_node_valid(np[l],TYPELEVEL_EMFIELD)){
      find_prim_emfield_mem_1(np, gl, l);
    }
  }
}


void update_prim_emfield_mem_in_zone_2(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (is_node_valid(np[l],TYPELEVEL_EMFIELD)){
      find_prim_emfield_mem_2(np, gl, l);
    }
  }
}


void update_prim_emfield_mem_in_zone_3(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (is_node_valid(np[l],TYPELEVEL_EMFIELD)){
      find_prim_emfield_mem_3(np, gl, l);
    }
  }
}


#ifdef _TSEMF_STORE_COEFFICIENTS
void update_prim_emfield_mem_in_zone_4(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,dim,flux;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (is_node_inner(np[l],TYPELEVEL_EMFIELD)){ 
      for (flux=0; flux<nfe; flux++){
        find_dtau_emfield(np,gl,l,flux,&(np[l].bs->dtauemfield[flux]));  
        np[l].bs->coeffp0sum[flux]=0.0;
        for (dim=0; dim<nd; dim++){
          find_linearization_coefficients_inner_node_emfield_interface(np, gl, l, dim, flux, &(np[l].bs->coeffm1[dim][flux]), &(np[l].bs->coeffp0[dim][flux]), &(np[l].bs->coeffp1[dim][flux]));
          np[l].bs->coeffp0sum[flux]+=np[l].bs->coeffp0[dim][flux];
        }
        find_linearization_coefficients_inner_node_emfield(np, gl, l, &(np[l].bs->tsemfcoeff[0][flux]), &(gl->tsemfcoeffzone));
      }
    }
  }
}
#endif



void update_prim_emfield_mem_in_zone(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np,gl,zone,&update_prim_emfield_mem_in_zone_1,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_valid,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_prim_emfield_mem_in_zone_2,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_valid,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_prim_emfield_mem_in_zone_3,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_valid,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#ifdef _TSEMF_STORE_COEFFICIENTS
  sweep_with_1D_segments(np,gl,zone,&update_prim_emfield_mem_in_zone_4,SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_valid,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
#endif
}


void add_convection_residual_emfield(long theta, long ls, long le, np_t *np, gl_t *gl){
  long l,flux;
  fluxemfield_t Fm1h;
  for (l=ls; l!=_l_plus_one(_l_plus_one(le,gl,theta),gl,theta); l=_l_plus_one(l,gl,theta)){
    find_Fstar_interface_emfield(np,gl,_al(gl,l,theta,-1),_al(gl,l,theta,+0),theta,Fm1h);
    for (flux=0; flux<nfe; flux++){
      if (l!=_l_plus_one(le,gl,theta)) np[l].bs->Resemfield[flux]-=Fm1h[flux];
      if (l!=ls) np[_al(gl,l,theta,-1)].bs->Resemfield[flux]+=Fm1h[flux];
    }
  }
}


void add_source_residual_emfield(long theta, long ls, long le, np_t *np, gl_t *gl){
  long l;
  long flux;
  fluxemfield_t S;
  if (theta==0) {
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      find_Sstar_emfield(np,gl,l,S);
      for (flux=0; flux<nfe; flux++) np[l].bs->Resemfield[flux]-=S[flux];
    }
  }
}


void update_residual_emfield(np_t *np, gl_t *gl, long theta, long ls, long le){
  add_convection_residual_emfield(theta,ls,le,np,gl);
  add_source_residual_emfield(theta,ls,le,np,gl);
}


void initialize_residual_emfield(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l,flux;
  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    for (flux=0; flux<nfe; flux++) np[l].bs->Resemfield[flux]=0.0e0;
    gl->effiter_R_emfield+=1.0e0/(double)gl->nn;
  }
}


void update_bdry_node_emfield(np_t *np, gl_t *gl, long l){
  long dim,l_C,l_B,l_A;
  long dimsgn; 
  bool BDRYDIRECFOUND;
#ifdef _2DL
  long dim1; long dim2;
#endif
#ifdef _3D
  long dim3;
#endif
  bool UPDATED;
  assert(is_node_bdry(np[l],TYPELEVEL_EMFIELD));
  UPDATED=FALSE;
  BDRYDIRECFOUND=find_bdry_direc(np, gl, l, TYPELEVEL_EMFIELD, &dim, &dimsgn);

  if (is_node_link(np[l],TYPELEVEL_EMFIELD)) {
    // in case the boundary node is a link, Uemf has already been updated 
    UPDATED=TRUE;
  }

  if (!UPDATED && BDRYDIRECFOUND){
    l_A=l;
    l_B=_al(gl,l,dim,dimsgn);
    l_C=_al(gl,l,dim,dimsgn*2);
    assert(is_node_inner(np[l_C],TYPELEVEL_EMFIELD));
    assert(is_node_inner(np[l_B],TYPELEVEL_EMFIELD));
    update_bdry_emfield(np,gl,l_A,l_B,l_C,dim,dimsgn,BDRYDIRECFOUND,TYPELEVEL_EMFIELD);
    UPDATED=TRUE;
  }
  /* now, do the corners */
  if (!UPDATED) {
#ifdef _2D
    for (dim1=-1; dim1<=1; dim1++){
      for (dim2=-1; dim2<=1; dim2++){
        l_C=_all(gl,l,0,dim1*2,1,dim2*2);
        l_B=_all(gl,l,0,dim1,1,dim2);
        l_A=l;
        if  (  is_node_inner(np[l_B],TYPELEVEL_EMFIELD)
            && is_node_inner(np[l_C],TYPELEVEL_EMFIELD) && !UPDATED){
          update_bdry_emfield(np,gl,l_A,l_B,l_C,dim,dimsgn,BDRYDIRECFOUND,TYPELEVEL_EMFIELD);
          UPDATED=TRUE;
        }
      }
    }
#endif
#ifdef _3D
    for (dim1=-1; dim1<=1; dim1++){
      for (dim2=-1; dim2<=1; dim2++){
        for (dim3=-1; dim3<=1; dim3++){
          l_C=_al(gl,
                      _al(gl,
                         _al(gl,l,0,dim1*2),
                       1,dim2*2),
                    2,dim3*2);
          l_B=_al(gl,
                      _al(gl,
                         _al(gl,l,0,dim1),
                       1,dim2),
                    2,dim3);
          l_A=l;
          if  (  is_node_inner(np[l_B],TYPELEVEL_EMFIELD)
              && is_node_inner(np[l_C],TYPELEVEL_EMFIELD) && !UPDATED){
            update_bdry_emfield(np,gl,l_A,l_B,l_C,dim,dimsgn,BDRYDIRECFOUND,TYPELEVEL_EMFIELD);
            UPDATED=TRUE;
          }
        }
      }
    }
#endif
  }
}


void update_bdry_nodes_on_segment_emfield(np_t *np, gl_t *gl, long theta, long ls, long le){
  long l;

  for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
    if (is_node_bdry(np[l],TYPELEVEL_EMFIELD)){
      thread_lock_node_set(np,l,THREADTYPE_ZONE);
      update_bdry_node_emfield(np, gl, l);
      thread_lock_node_unset(np,l,THREADTYPE_ZONE);
    }
  }
}


void update_bdry_nodes_emfield(np_t *np, gl_t *gl, zone_t zone){
  sweep_with_1D_segments(np, gl, zone, &update_bdry_nodes_on_segment_emfield, SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_valid,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
}


void find_residual_emfield(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;

  /* now, let's find the residual and store it in bs->dUstaremfield*/
  sweep_with_1D_segments(np,gl,zone,&initialize_residual_emfield, SWEEPTYPE_I, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_LIGHT,GRIDLEVEL_ONE);
  sweep_with_1D_segments(np,gl,zone,&update_residual_emfield, SWEEPTYPE_IJK, TYPELEVEL_EMFIELD,&is_node_inner,SEGMENTWORK_HEAVY,GRIDLEVEL_ONE);

  /* let's find max residual, and put it in gl*/
  for_ijk(zone,is,js,ks,ie,je,ke){
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
          np[_ai(gl,i,j,k)].bs->_xi_emfield=_xi_emfield(np[_ai(gl,i,j,k)],gl,np[_ai(gl,i,j,k)].bs->Resemfield);
        }
  }
}


void find_ximax_emfield(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;
  gl->ximax_emfield=0.0e0;
  for_ijk(zone,is,js,ks,ie,je,ke){
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD) && np[_ai(gl,i,j,k)].bs->_xi_emfield>=gl->ximax_emfield) {
	  gl->ximax_emfield=np[_ai(gl,i,j,k)].bs->_xi_emfield;
	  gl->i_ximax_emfield=i;
	  gl->j_ximax_emfield=j;
	  gl->k_ximax_emfield=k;
	}
  }
}


void read_UpdateEMField_arguments(char **argum, SOAP_codex_t *codex, gl_t *gl){
  SOAP_substitute_all_argums(argum, codex);
  gl->Lc=SOAP_get_argum_double(codex,*argum,0);
  gl->relaxEMF=SOAP_get_argum_double(codex,*argum,1);
  gl->numsubiter_tsemf=4; /* make the default number of subiterations equal to 4 */
  gl->tsemfmethod=TSEMF_DEFAULT;
  if (gl->Lc<=0.0) fatal_error("The length scale Lc must be positive when calling UpdateEMField().");
  if (gl->relaxEMF<=0.0) fatal_error("The relaxation factor relaxEMF must be positive when calling UpdateEMField().");
  if (gl->relaxEMF>2.0) fatal_error("The relaxation factor relaxEMF must be less than 2 when calling UpdateEMField().");
  if (gl->numsubiter_tsemf<=0.0) fatal_error("The number of subiterations subiter_tsemf must be positive when calling UpdateEMField().");
#ifdef UNSTEADY
  gl->dt=SOAP_get_argum_double(codex,*argum,2);
  if (gl->dt<=0.0) fatal_error("The time step dt must be positive when calling UpdateEMField().");
  if (SOAP_number_argums(*argum)>3) gl->tsemfmethod=SOAP_get_argum_long(codex,*argum,3);
  if (SOAP_number_argums(*argum)>4){
    if (gl->tsemfmethod==TSEMF_SOR || gl->tsemfmethod==TSEMF_SOR2 || gl->tsemfmethod==TSEMF_ADIIMAF || gl->tsemfmethod==TSEMF_IMAF || gl->tsemfmethod==TSEMF_IMAFk  || gl->tsemfmethod==TSEMF_IMAFi) 
      gl->numsubiter_tsemf=SOAP_get_argum_long(codex,*argum,4);
    else fatal_error("UpdateEMField accepts the number of subiterations as a 5th argument only if TSEMF_SOR, TSEMF_SOR2, TSEMF_ADIIMAF, TSEMF_IMAF, TSMEF_IMAFk, TSMEF_IMAFi is specified.");
  }
#else 
  if (SOAP_number_argums(*argum)>2) gl->tsemfmethod=SOAP_get_argum_long(codex,*argum,2);
  if (SOAP_number_argums(*argum)>3) {
    if (gl->tsemfmethod==TSEMF_SOR || gl->tsemfmethod==TSEMF_SOR2 || gl->tsemfmethod==TSEMF_ADIIMAF || gl->tsemfmethod==TSEMF_IMAF || gl->tsemfmethod==TSEMF_IMAFk  || gl->tsemfmethod==TSEMF_IMAFi) 
      gl->numsubiter_tsemf=SOAP_get_argum_long(codex,*argum,3);
    else fatal_error("UpdateEMField accepts the number of subiterations as a 4th argument only if TSEMF_SOR, TSEMF_SOR2, TSEMF_ADIIMAF, TSEMF_IMAF, TSEMF_IMAFk, TSMEF_IMAFi is specified.");
  }
#endif
}


void solve_TDMA_emfield(np_t *np, gl_t *gl, long theta, long ls, long le, int TYPELEVEL, EXM_tdmaline_t *tdma, long numlines){
#ifdef DISTMPI
  long line,cnt,i,j,k,i_s,j_s,k_s;
  double tmp;
  MPI_Status MPI_Status1;

  if (gl->EM_MPIBDRY_EXPLICIT){
    EXM_solve_TDMA(tdma, numlines);
  } else {
    /* if ls node is inner node, need to obtain the tdma[0] from another process that owns ls */
    if (is_node_inner(np[ls],TYPELEVEL)){
      find_ijk_from_l(gl, ls, &i, &j, &k);
      assert(_ai_all(gl,i,j,k)<LONG_MAX);
      if (MPI_Recv(tdma[0].val,4,MPI_DOUBLE,_node_rank(gl,i,j,k),_ai_all(gl,i,j,k),MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in solve_TDMA_emfield");
      assert(tdma[0].val[0]==0.0);
    } 
    for (line=0; line<numlines-1; line++){
      assert(tdma[line].val[1]!=0.0);
      tmp = -(tdma[line+1].val[0] / tdma[line].val[1]);
      for (cnt = 1; cnt <= 2; cnt++)
        tdma[line+1].val[cnt - 1] += tdma[line].val[cnt] * tmp;
      tdma[line+1].val[3] += tdma[line].val[3] * tmp;
      tdma[line+1].val[0] = 0.0;
    }
    /* if le node is inner node, need to send the tdma[numlines-2] to another process that owns le */
    if (is_node_inner(np[le],TYPELEVEL)){
      find_ijk_from_l(gl, le, &i, &j, &k);
      find_ijk_from_l(gl, _l_minus_one(le,gl,theta), &i_s, &j_s, &k_s);
      assert(_ai_all(gl,i,j,k)<LONG_MAX);
      if (MPI_Send(tdma[numlines-2].val,4,MPI_DOUBLE,_node_rank(gl,i,j,k),_ai_all(gl,i_s,j_s,k_s),MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in solve_TDMA_emfield");
    } 
    /* if le node is inner node, need to obtain the tdma[numlines-1] from another process that owns le */
    if (is_node_inner(np[le],TYPELEVEL)){
      find_ijk_from_l(gl, le, &i, &j, &k);
      assert(_ai_all(gl,i,j,k)<LONG_MAX);
      if (MPI_Recv(tdma[numlines-1].val,4,MPI_DOUBLE,_node_rank(gl,i,j,k),_ai_all(gl,i,j,k),MPI_COMM_WORLD,&MPI_Status1)!=MPI_SUCCESS) fatal_error("MPI_Recv problem in solve_TDMA_emfield");
      assert(tdma[numlines-1].val[2]==0.0);
    } 
    for (line=numlines-1; line>0; line--){
      assert(tdma[line].val[1]!=0.0);
      tdma[line].val[3] /= tdma[line].val[1];
      tdma[line].val[1] = 1.0;
      tdma[line-1].val[3] -= tdma[line].val[3] * tdma[line-1].val[2];
      tdma[line-1].val[2] = 0.0;
    }
    assert(tdma[0].val[1]!=0.0);
    tdma[0].val[3] /= tdma[0].val[1];
    tdma[0].val[1] = 1.0;
    /* if ls node is inner node, need to send the tdma[1] to another process that owns ls */
    if (is_node_inner(np[ls],TYPELEVEL)){
      find_ijk_from_l(gl, ls, &i, &j, &k);
      find_ijk_from_l(gl, _l_plus_one(ls,gl,theta), &i_s, &j_s, &k_s);
      assert(_ai_all(gl,i,j,k)<LONG_MAX);
      if (MPI_Send(tdma[1].val,4,MPI_DOUBLE,_node_rank(gl,i,j,k),_ai_all(gl,i_s,j_s,k_s),MPI_COMM_WORLD)!=MPI_SUCCESS) fatal_error("MPI_Send problem in solve_TDMA_emfield");
    } 
  }
#else
  EXM_solve_TDMA(tdma, numlines);
#endif
}
#endif//EMFIELD



