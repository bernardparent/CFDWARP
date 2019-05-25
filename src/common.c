// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1999-2006 Bernard Parent

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

#include <src/common.h>
#include <model/_model.h>
#include <stdarg.h>
#include <cycle/_cycle.h>
#include <src/bdry.h>


#ifdef MLOCK
  #include <sys/mman.h>
#endif


size_t wfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream){
#ifdef DISTMPI
  int rank,proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  if (rank==0) {
    return(fwrite(ptr,size,nmemb,stream));
  } else {
    return(0);
  }
#else
  return(fwrite(ptr,size,nmemb,stream));
#endif
}

FILE *wfopen(const char *path, const char *mode){
#ifdef DISTMPI
  int rank,proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  if (rank==0) {
    return(fopen(path,mode));
  } else {
    return(NULL);
  }
#else
  return(fopen(path,mode));
#endif
}


int wfclose(FILE *fp){
#ifdef DISTMPI
  int rank,proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  if (rank==0) {
    return(fclose(fp));
  } else {
    return(0);
  }
#else
  return(fclose(fp));
#endif

}


/* this function calls vfprintf with the same arguments
   as it is given. If the call is to stdout, flushes it.
   Needed to deal with MPI, so that only process 0 outputs things
   to stdout. The output to stdout from other processes is ignored*/
int wfprintf(FILE *stream, const char *format, ...){
  va_list ap;
  int retval;
  FILE *streamlocal;
  streamlocal=stream;
  // make the stream equal to stdout if it is set to stderr
  if (streamlocal==stderr) streamlocal=stdout;
  retval=0;
#ifdef DISTMPI
  int rank,proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  if (rank==0 || stream==stderr) {
    if (stream==stderr && proc>1) fprintf(stdout,"[MPI rank=%d]",rank);
    va_start(ap, format);
    retval=vfprintf(streamlocal, format, ap);
    va_end(ap);
//      fprintf(stdout,"[rank=%d]",rank);
    if (streamlocal==stdout) fflush(stdout);
  }
#else
  va_start(ap, format);
  retval=vfprintf(streamlocal, format, ap);
  va_end(ap);
  if (streamlocal==stdout) fflush(stdout);
#endif
  return(retval);
}


/* this function calls wfprintf with the same arguments
   as it is given and exits.*/
void fatal_error(const char *formatstr, ...){
  va_list ap;
  char *newstr;
  int term_width,term_height;
  newstr=(char *)malloc(10000*sizeof(char));
  wfprintf(stderr,"\n\n");
  va_start(ap, formatstr);
  vsprintf(newstr,formatstr, ap);
  va_end(ap);
  find_terminal_window_size(&term_width,&term_height);
  wfprintf(stderr,"%s",strwrp(newstr,min(term_width-1,70)));
  free(newstr);
  wfprintf(stderr,"\n\nCFDWARP fatal error. Exiting.\n\n");
#ifdef DISTMPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#endif
  exit(EXIT_FAILURE);
}



long _ai(gl_t *gl, long i, long j, long k) {
  long ii;
  ii=(i-gl->domain_lim.is);
  if2DL(
    ii=ii*(gl->domain_lim.je-gl->domain_lim.js+1)+(j-gl->domain_lim.js);
    if3DL(ii=ii*(gl->domain_lim.ke-gl->domain_lim.ks+1)+(k-gl->domain_lim.ks);)
  )
  return(ii);
}


long _ai_all(gl_t *gl, long i, long j, long k) {
  long ii;
  ii=(i-gl->domain_lim_all.is);
  if2DL(
    ii=ii*(gl->domain_lim_all.je-gl->domain_lim_all.js+1)+(j-gl->domain_lim_all.js);
    if3DL(ii=ii*(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1)+(k-gl->domain_lim_all.ks);)
  )
  return(ii);
}


void find_ijk_from_l(gl_t *gl, long l, long *i, long *j, long *k){
#ifdef _2D
    *j=mod(l,gl->domain_lim.je-gl->domain_lim.js+1)+gl->domain_lim.js;
    *i=l/(gl->domain_lim.je-gl->domain_lim.js+1)+gl->domain_lim.is;
    *k=0;
#endif
#ifdef _3D
    *k=mod(l,gl->domain_lim.ke-gl->domain_lim.ks+1)+gl->domain_lim.ks;
    *j=mod(l/(gl->domain_lim.ke-gl->domain_lim.ks+1),(gl->domain_lim.je-gl->domain_lim.js+1))
          +gl->domain_lim.js;
    *i=(l/(gl->domain_lim.ke-gl->domain_lim.ks+1))/(gl->domain_lim.je-gl->domain_lim.js+1)
          +gl->domain_lim.is;
#endif
}


void find_ijk_from_l_all(gl_t *gl, long l, long *i, long *j, long *k){
#ifdef _2D
    *j=mod(l,gl->domain_lim_all.je-gl->domain_lim_all.js+1)+gl->domain_lim_all.js;
    *i=l/(gl->domain_lim_all.je-gl->domain_lim_all.js+1)+gl->domain_lim_all.is;
    *k=0;
#endif
#ifdef _3D
    *k=mod(l,gl->domain_lim_all.ke-gl->domain_lim_all.ks+1)+gl->domain_lim_all.ks;
    *j=mod(l/(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1),(gl->domain_lim_all.je-gl->domain_lim_all.js+1))
          +gl->domain_lim_all.js;
    *i=(l/(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1))/(gl->domain_lim_all.je-gl->domain_lim_all.js+1)
          +gl->domain_lim_all.is;
#endif
}


long _l_from_l_all(gl_t *gl, long l_all){
  long i,j,k;
  long l;
  find_ijk_from_l_all(gl, l_all, &i, &j, &k);
  l=_ai(gl,i,j,k);
  assert(is_node_in_domain_lim(l,gl));
  return(l);
}
 

long _i(long l, gl_t *gl, long dim){
  long i[3];
  find_ijk_from_l(gl,l,&i[0],&i[1],&i[2]);
  return(i[dim]);
}


long _i_all(long l_all, gl_t *gl, long dim){
  long i[3];
  find_ijk_from_l_all(gl,l_all,&i[0],&i[1],&i[2]);
  return(i[dim]);
}


long _al(gl_t *gl, long l, long theta, long offset) {
  long ii;
  ii=l+offset;
#ifdef _2D
  if (theta==0) ii=l+offset*(gl->domain_lim.je-gl->domain_lim.js+1);
#endif
#ifdef _3D
  if (theta==0) {
    ii=l+offset*(gl->domain_lim.je-gl->domain_lim.js+1)
               *(gl->domain_lim.ke-gl->domain_lim.ks+1);
  }
  if (theta==1) {
    ii=l+offset*(gl->domain_lim.ke-gl->domain_lim.ks+1);
  }
#endif
  return(ii);
}


long _al_check(gl_t *gl, long l, long theta, long offset) {
  long ii;
  ii=_al(gl,l,theta,offset);
  assert(is_node_in_domain_lim(ii,gl));
  return(ii);
}


long _all(gl_t *gl, long l, long theta1, long offset1, long theta2, long offset2) {
  long ii;
  double ltmp;
  ltmp=_al(gl, l, theta1, offset1);
  ii=_al(gl, ltmp, theta2, offset2);
  return(ii);
}


/* find the node number of a node distanced by offset node(s) from llink perpendicular to link surface  
   note: llink must be an inner node and a link node; lbdry is the bdry node that is linked to llink */
long _al_link(np_t *np, gl_t *gl, long llink, long lbdry, long offset, int TYPELEVEL){
  long cnt,theta,thetasgn,lret; 
  assert(is_node_link(np[llink],TYPELEVEL));
  assert(is_node_inner(np[llink],TYPELEVEL));
  cnt=0;
  find_link_direc(np, gl, llink, lbdry, TYPELEVEL, &theta, &thetasgn);
  do {
    cnt++;
  } while(
       cnt<offset
    && is_node_inner(np[_al(gl,llink,theta,cnt*thetasgn)],TYPELEVEL) 
#ifdef DISTMPI
    && is_node_in_domain_lim(_al(gl,llink,theta,cnt*thetasgn), gl) 
#endif
    );

  lret=_al(gl,llink,theta,thetasgn*cnt);
  if (!is_node_valid(np[lret],TYPELEVEL)) fatal_error("Problem in _al_link.");
  return(lret);
}


/* find the node number of the nearest inner node from node l */
bool find_l_of_nearest_inner_node(np_t *np, gl_t *gl, long l, int nodetype, long *linner){
  bool FOUND;
  long dim,dim2;

  FOUND=FALSE;
  /* first check if l is an inner node */
  if (is_node_inner(np[l],nodetype)) {
    *linner=l;
    FOUND=TRUE;
  }
  /* second check if the nodes opposite the faces are inner */
  if (!FOUND){
    for (dim=0; dim<nd; dim++){
      if (is_node_inner(np[_al(gl,l,dim,+1)],nodetype) && !FOUND){
        *linner=_al(gl,l,dim,+1);
        FOUND=TRUE;
      }
      if (is_node_inner(np[_al(gl,l,dim,-1)],nodetype) && !FOUND){
        *linner=_al(gl,l,dim,-1);
        FOUND=TRUE;
      }
    }
  }
  /* third check if the corner nodes are inner */
  if (!FOUND){
    for (dim=0; dim<nd; dim++){
      for (dim2=0; dim2<nd; dim2++){
        if (dim2!=dim){
          if (is_node_inner(np[_all(gl,l,dim,+1,dim2,+1)],nodetype) && !FOUND){
            *linner=_all(gl,l,dim,+1,dim2,+1);
            FOUND=TRUE;
          }
          if (is_node_inner(np[_all(gl,l,dim,+1,dim2,-1)],nodetype) && !FOUND){
            *linner=_all(gl,l,dim,+1,dim2,-1);
            FOUND=TRUE;
          }
          if (is_node_inner(np[_all(gl,l,dim,-1,dim2,+1)],nodetype) && !FOUND){
            *linner=_all(gl,l,dim,-1,dim2,+1);
            FOUND=TRUE;
          }
          if (is_node_inner(np[_all(gl,l,dim,-1,dim2,-1)],nodetype) && !FOUND){
            *linner=_all(gl,l,dim,-1,dim2,-1);
            FOUND=TRUE;
          }
        }
      }
    }
  }
  return(FOUND);
}


long _l_plus_one(long l, gl_t *gl, long theta){
  long tmp;
  tmp=_al(gl,l,theta,+1);
  return(tmp);
}

long _l_minus_one(long l, gl_t *gl, long theta){
  long tmp;
  tmp=_al(gl,l,theta,-1);
  return(tmp);
}


long _l_plus_two(long l, gl_t *gl, long theta){
  long tmp;
  tmp=_al(gl,l,theta,+2);
  return(tmp);
}

long _l_minus_two(long l, gl_t *gl, long theta){
  long tmp;
  tmp=_al(gl,l,theta,-2);
  return(tmp);
}


long _nodes_from_bdry(np_t *np, gl_t *gl, long l, long theta, int TYPELEVEL){
  long cnt;
  cnt=0; 
  if (!is_node_valid(np[l],TYPELEVEL)) fatal_error("_nodes_from_bdry() must be called with l a valid node.");
  if (!is_node_bdry(np[l],TYPELEVEL)){
    do {
      cnt++; 
    } while(
            is_node_inner(np[_al(gl,l,theta,cnt)],TYPELEVEL) 
         && is_node_inner(np[_al(gl,l,theta,-cnt)],TYPELEVEL)
#ifdef DISTMPI
//         && is_node_in_domain_lim(_al(gl,l,theta,-cnt), gl)  
//         && is_node_in_domain_lim(_al(gl,l,theta,cnt), gl) 
#endif
      );
  }
#ifdef DISTMPI
  if (!is_node_in_domain_lim(_al(gl,l,theta,cnt), gl) || !is_node_in_domain_lim(_al(gl,l,theta,-cnt), gl))
    fatal_error("MPI problem in _nodes_from_bdry() part of common.c.");
#endif
  return(cnt);
}


long _nodes_from_bdry_limited(np_t *np, gl_t *gl, long l, long theta, int TYPELEVEL, long limit){
  long cnt;
  cnt=0; 
  if (!is_node_valid(np[l],TYPELEVEL)) fatal_error("_nodes_from_bdry() must be called with l a valid node.");
  if (!is_node_bdry(np[l],TYPELEVEL)){
    do {
      cnt++; 
    } while(
            cnt<limit
         && is_node_inner(np[_al(gl,l,theta,cnt)],TYPELEVEL) 
         && is_node_inner(np[_al(gl,l,theta,-cnt)],TYPELEVEL)
#ifdef DISTMPI
//         && is_node_in_domain_lim(_al(gl,l,theta,-cnt), gl)  
//         && is_node_in_domain_lim(_al(gl,l,theta,cnt), gl) 
#endif
      );
  }
#ifdef DISTMPI
  if (!is_node_in_domain_lim(_al(gl,l,theta,cnt), gl) || !is_node_in_domain_lim(_al(gl,l,theta,-cnt), gl))
    fatal_error("MPI problem in _nodes_from_bdry_limited part of common.c.");
#endif
  return(cnt);
}


long _nodes_from_bdry_through_links_limited(np_t *np, gl_t *gl, long l, long theta, int TYPELEVEL, long limit){
  long cnt1,cnt2,cnt,cnt1link,cnt2link;
#ifndef DISTMPI
  long llink,lbdry,thetalink,thetalinksgn;
#endif
  cnt1=0;
  cnt1link=0;
  cnt2=0; 
  cnt2link=0;
  if (!is_node_valid(np[l],TYPELEVEL)) fatal_error("_nodes_from_bdry() must be called with l a valid node.");
  if (!is_node_bdry(np[l],TYPELEVEL)){
    do {
      cnt1++; 
    } while(
            cnt1<limit
         && is_node_inner(np[_al(gl,l,theta,cnt1)],TYPELEVEL) 
#ifdef DISTMPI
         && is_node_in_domain_lim(_al(gl,l,theta,cnt1), gl) 
#endif
      );
    do {
      cnt2++; 
    } while(
            cnt2<limit
         && is_node_inner(np[_al(gl,l,theta,-cnt2)],TYPELEVEL)
#ifdef DISTMPI
         && is_node_in_domain_lim(_al(gl,l,theta,-cnt2), gl)  
#endif
      );
  } 

  if (is_node_link(np[_al(gl,l,theta,cnt1)],TYPELEVEL) && cnt1<limit){
#ifdef DISTMPI
    cnt1link=np[_al(gl,l,theta,cnt1)].numlinkmusclvars+1;
#else
    lbdry=_al(gl,l,theta,cnt1);
    assert(is_node_bdry(np[lbdry],TYPELEVEL));
    llink=_node_link(np[lbdry],0,TYPELEVEL);
    find_link_direc(np, gl, llink, lbdry, TYPELEVEL, &thetalink, &thetalinksgn);
    do {
      cnt1link++;
    } while(
          cnt1+cnt1link<limit
       && is_node_inner(np[_al(gl,llink,thetalink,cnt1link*thetalinksgn)],TYPELEVEL) 
      );
#endif
  }

  if (is_node_link(np[_al(gl,l,theta,-cnt2)],TYPELEVEL) && cnt2<limit){
#ifdef DISTMPI
    cnt2link=np[_al(gl,l,theta,-cnt2)].numlinkmusclvars+1;
#else
    lbdry=_al(gl,l,theta,-cnt2);
    assert(is_node_bdry(np[lbdry],TYPELEVEL));
    llink=_node_link(np[lbdry],0,TYPELEVEL);
    find_link_direc(np, gl, llink, lbdry, TYPELEVEL, &thetalink, &thetalinksgn);
    do {
      cnt2link++;
    } while(
          cnt2+cnt2link<limit
       && is_node_inner(np[_al(gl,llink,thetalink,cnt2link*thetalinksgn)],TYPELEVEL) 
      );
#endif
  }

  cnt=min(limit,min(cnt1+cnt1link,cnt2+cnt2link));
/*  if (theta==0 && cnt1+cnt1link<5) {
    printf("%ld%ld%ld%ld-",cnt1+cnt1link,cnt2+cnt2link,cnt1,cnt2);
  }*/
  return(cnt);
}


long _nodes_between_link_and_bdry_limited(np_t *np, gl_t *gl, long llink, long lbdry, int TYPELEVEL, long limit){
  long cnt,theta,thetasgn; 
  assert(is_node_link(np[llink],TYPELEVEL));
  assert(is_node_inner(np[llink],TYPELEVEL));
  cnt=0;
  find_link_direc(np, gl, llink, lbdry, TYPELEVEL, &theta, &thetasgn);
  do {
    cnt++;
  } while(
       cnt<=limit
    && is_node_valid(np[_al(gl,llink,theta,cnt*thetasgn)],TYPELEVEL) 
#ifdef DISTMPI
    && is_node_in_domain_lim(_al(gl,llink,theta,cnt*thetasgn), gl) 
#endif
    );
  cnt--;
  return(cnt);
}




void multiply_matrix_by_constant(double multfact, sqmat_t A){
  long row,col;

  for (row=0; row<nf; row++) {
    for (col=0; col<nf; col++) {
       A[row][col]=A[row][col]*multfact;
    }
  }
}

void display_matrix(sqmat_t A){
  long col,row;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      wfprintf(stdout,"%+12.5E ",A[row][col]);
    }
    wfprintf(stdout,"\n");
  }
  wfprintf(stdout,"\n");
}

void multiply_matrix_and_matrix(sqmat_t A, sqmat_t B, sqmat_t C){
  long row,col,flux;
  double tmp;

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      tmp=0.0e0;
      for (flux=0; flux<nf; flux++){
        tmp=tmp+A[row][flux]*B[flux][col];
      }
      C[row][col]=tmp;
    }
  }
}

void arithmetic_mean_of_two_matrices(sqmat_t A, sqmat_t B, sqmat_t C){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      C[row][col]=0.5e0*(A[row][col]+B[row][col]);
    }
  }
}



void harmonic_mean_of_two_matrices(sqmat_t A, sqmat_t B, sqmat_t C){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      C[row][col]=2.0e0/(1.0/notzero(A[row][col],1.0e-39)+1.0/notzero(B[row][col],1.0e-39));
    }
  }
}


void make_matrix_positive(sqmat_t A){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      A[row][col]=max(0.0,A[row][col]);
    }
  }
}


void make_matrix_diagonal(sqmat_t A){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      if (row!=col) A[row][col]=0.0;
    }
  }
}


void multiply_matrix_and_diagonal_matrix(sqmat_t A, sqmat_t B, sqmat_t C){
  long row,col;

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      C[row][col]=A[row][col]*B[col][col];
    }
  }
}

void multiply_diagonal_matrix_and_matrix(sqmat_t A, sqmat_t B, sqmat_t C){
  long row,col;

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      C[row][col]=A[row][row]*B[row][col];
    }
  }
}

void set_matrix_to_identity(sqmat_t A){
  long row,col;

  for (row=0; row<nf; row++) {
    for (col=0; col<nf; col++) {
       A[row][col]=0.0e0;
    }
  }
  for (row=0; row<nf; row++) {
    A[row][row]=1.0e0;
  }
}

void set_matrix_to_zero(sqmat_t A){
  long row,col;

  for (row=0; row<nf; row++) {
    for (col=0; col<nf; col++) {
       A[row][col]=0.0e0;
    }
  }
}

void set_matrix_to_matrix(sqmat_t A, sqmat_t B){
  long row,col;

  for (row=0; row<nf; row++) {
    for (col=0; col<nf; col++) {
       A[row][col]=B[row][col];
    }
  }
}


void set_vector_to_zero(flux_t S){
  long flux;
  for (flux=0; flux<nf; flux++) S[flux]=0.0e0;
}




/* C=A+B */
void add_two_matrices(sqmat_t A, sqmat_t B, sqmat_t C){
  long row,col;

  for (row=0; row<nf; row++) {
    for (col=0; col<nf; col++) {
       C[row][col]=A[row][col]+B[row][col];
    }
  }
}

void copy_matrix(sqmat_t AA, sqmat_t BB){
  long row,col;

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      BB[row][col]=AA[row][col];
    }
  }
}


void copy_vector(flux_t AA, flux_t BB){
  long row;
  for (row=0; row<nf; row++){
    BB[row]=AA[row];
  }
}

void multiply_matrix_and_vector(sqmat_t sqmat, flux_t mat1, flux_t mat2){
  long row,col;
  double tmp;

  for (row=0; row<nf; row++){
     tmp=0.0e0;
     for (col=0; col<nf; col++){
     /* wfprintf(stdout,"row=%ld  col=%ld  sqmat=%E   mat1=%E\n",row,col,sqmat[row][col],mat1[col]); */
       tmp=tmp+sqmat[row][col]*mat1[col];
     }

     mat2[row]=tmp;
  }
}

double _subdet(sqmat_t mat, long numrow){
  double tmp;
  sqmat_t matnew;
  long col2,row2,col,row;
  double tmp2;
  if (numrow==2) {
    tmp=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
  } else {
    for (row=0; row<nf; row++){
      for (col=0; col<nf; col++){
        matnew[row][col]=0.0; //needed to avoid compiler warning
      }
    }
    tmp=0.0e0;
    for (col=0; col<numrow; col++){
      for (col2=0; col2<nf; col2++){
        for (row2=1; row2<nf; row2++){
           matnew[row2-1][col2]=mat[row2][col2];
        }
      }
      for (col2=col; col2<(numrow-1); col2++){
        for (row2=0; row2<(numrow-1); row2++){
          matnew[row2][col2]=matnew[row2][col2+1];
        }
      }
      if (mod(col,2)==0) {
        tmp2=1.0e0;
      } else {
        tmp2=-1.0e0;
      }
      tmp=tmp+tmp2*mat[0][col]*_subdet(matnew,numrow-1);
    }
  }
  return(tmp);
}

void find_determinant(sqmat_t mat1, double *det){
  *det=_subdet(mat1,nf);
}

void det_mat_old(sqmat_t mat1, double *det){
  long row,rowc,colc;
  double fact;

  for (row=0; row<nf; row++){
    for (rowc=0; rowc<nf; rowc++){
      if (rowc!=row) {
        fact=-mat1[rowc][row]/mat1[row][row];
        for (colc=0; colc<nf; colc++) {
          mat1[rowc][colc]=mat1[rowc][colc]+
            fact*mat1[row][colc];
        }
      }
    }
  }
  *det=1.0e0;
  for (row=0; row<nf; row++){
    *det=*det*mat1[row][row];
  }
}

void multiply_diagonal_matrix_and_vector(sqmat_t sqmat, flux_t mat1, flux_t mat2){
  long row;

  for (row=0; row<nf; row++){
     mat2[row]=sqmat[row][row]*mat1[row];
  }
}


void invert_matrix_gaussian_elimination(sqmat_t mattmp, sqmat_t mat2){
  long row,row2,col;
  double multfact;
  sqmat_t mat1;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      mat1[row][col]=mattmp[row][col];
    }
  }
  set_matrix_to_identity(mat2);
    for (row=0; row<nf; row++){
      for (row2=0; row2<nf; row2++){
        if (row2 != row) {
          assert(mat1[row][row]!=0.0);
          multfact=-mat1[row2][row]/(mat1[row][row]);

          for (col=0; col<nf; col++){
            mat1[row2][col]=mat1[row2][col]+mat1[row][col]*multfact;
            mat2[row2][col]=mat2[row2][col]+mat2[row][col]*multfact;
	  }
        }
      }
    }
    for (row=0; row<nf; row++){
      for (col=0; col<nf; col++){
        assert(mat1[row][row]!=0.0);
        mat2[row][col]=mat2[row][col]/(mat1[row][row]);
      }
      mat1[row][row]=1.0e0;
    }
}


// function by Jaehyuk Lee 
void invert_matrix_partial_pivoting(sqmat_t mattmp, sqmat_t mat2){
  long p, row,row2,col;
  double temp,m,multfact;
  sqmat_t mat1;
  for(row=0;row<nf;row++){
    for(col=0;col<nf;col++){
      mat1[row][col]=mattmp[row][col];
    }
  }
  set_matrix_to_identity(mat2);
  for(row=0;row<nf-1;row++){
    p=row;
    m=mat1[row][row];
    for(row2=row+1;row2<nf;row2++){
      if(fabs(m)<fabs(mat1[row2][row])){
        p=row2;
        m=mat1[row2][row];
      }
    }
    if(p!=row){
      for(col=0;col<nf;col++){
      temp=mat1[p][col];
      mat1[p][col]=mat1[row][col];
      mat1[row][col]=temp;
      temp=mat2[p][col];
      mat2[p][col]=mat2[row][col];
      mat2[row][col]=temp;
      }
    }
    for(row2=row+1;row2<nf;row2++){
      assert(mat1[row][row]!=0.0e0);
      multfact=-(mat1[row2][row])/(mat1[row][row]);
      for(col=0;col<nf;col++){
        mat2[row2][col]=mat2[row2][col]+mat2[row][col]*multfact;
      }
      for(col=row;col<nf;col++){
        mat1[row2][col]=mat1[row2][col]+mat1[row][col]*multfact;
      }mat1[row2][row]=0.0e0;
    }
  }
  for(row=nf-1;row>=1;row--){
    for(row2=row-1;row2>=0;row2--){
      assert(mat1[row][row]!=0.0e0);
      multfact=-mat1[row2][row]/mat1[row][row];
      for(col=0;col<nf;col++){
        mat2[row2][col]=mat2[row2][col]+multfact*mat2[row][col];
      }
    }
  }
  for(row=0;row<nf;row++){
    for(col=0;col<nf;col++){
      assert(mat1[row][row]!=0.0e0);
      mat2[row][col]=mat2[row][col]/mat1[row][row];
    }
  }
}


void invert_matrix(sqmat_t mattmp, sqmat_t mat2){
  invert_matrix_gaussian_elimination(mattmp,mat2);
  //invert_matrix_partial_pivoting(mattmp,mat2);
}


void invert_diagonal_matrix(sqmat_t mat1, sqmat_t mat2){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      mat2[row][col]=0.0;
    }
  }
  for (row=0; row<nf; row++){
    assert(mat1[row][row]!=0.0);
    mat2[row][row]=1.0/mat1[row][row];
  }
}



void create_node(np_t *np, long i, long j, long k){
  np->bs=(npbs_t *) malloc(sizeof(npbs_t));
  np->bdryparam=NULL;
  np->linkarray=NULL;
  if (np->bs==NULL) {
    perror("malloc");
    fatal_error("Problem malloc'ing np->bs in create_node subroutine.");
  }
  np->wk=(npwk_t *) malloc(sizeof(npwk_t));
  if (np->wk==NULL) {
    perror("malloc");
    fatal_error("Problem malloc'ing np->wk in create_node subroutine.");
  }
  np->status='A';
  np->FLUIDPRIMMEM=FALSE;
#ifdef ZONETHREADS
  thread_lock_init(&(np->wk->lock), THREADTYPE_ZONE);
#endif
  np->type=-2;
  np->type_wk=-2;
#ifdef EMFIELD
  np->type_emf=-2;
#endif
#ifdef _TSEMF_SOR2
  np->bs->tsemfcoeff=NULL;
  np->bs->tsemfnode=NULL;
  np->bs->tsemfnodenum=0;
#endif
#ifdef MLOCK
  if (!(mlock(np->bs,sizeof(npbs_t)) &&
          mlock(np->wk,sizeof(npwk_t)) &&
          mlock(np,sizeof(np_t))
       )) wfprintf(stderr,"problem locking memory.\n");
#endif
}

bool is_node_resumed(np_t np){
  if (np.status=='A') {
    return(TRUE);
  } else {
    return(FALSE);
  }
}

bool is_node_suspended(np_t np){
  if (np.status=='S') {
    return(TRUE);
  } else {
    return(FALSE);
  }
}

void suspend_node(np_t *np){
  if (is_node_resumed(*np)) {
#ifdef MLOCK
    if (!(munlock(np->bs,sizeof(npbs_t)) &&
            munlock(np->wk,sizeof(npwk_t)) &&
            munlock(np,sizeof(np_t))
         )) wfprintf(stderr,"problem unlocking memory.\n");
#endif
#ifdef ZONETHREADS
    thread_lock_destroy(&(np->wk->lock),THREADTYPE_ZONE);
#endif
    free(np->wk);
    np->status='S';
    np->FLUIDPRIMMEM=FALSE;
  }
}


bool resume_node(np_t *np){
  if (is_node_suspended(*np)) {
    np->wk=(npwk_t *) malloc(sizeof(npwk_t));
    if (np->wk==NULL) {
      perror("malloc");
      fatal_error("Problem malloc'ing np->wk in resume_node subroutine.");
    }
    np->type_wk=np->type;
    np->wk->xi=-10.0e99; /* this line should not be necessary, but when
                         xi is read later on in cycle_share.c in find_ximax,
                         xi is sometimes not a numerical quantity. Need to
                         understand why this is so.. */
    np->status='A';
#ifdef ZONETHREADS
    thread_lock_init(&(np->wk->lock), THREADTYPE_ZONE);
#endif
#ifdef MLOCK
    if (!(mlock(np->bs,sizeof(npbs_t)) &&
            mlock(np->wk,sizeof(npwk_t)) &&
            mlock(np,sizeof(np_t))
         )) wfprintf(stderr,"problem locking memory.\n");
#endif
#ifdef _RESCONV_STORAGE_FSTAR
    np->wk->Fp1h=NULL;
#endif

    return(TRUE);
  } else {
    return(FALSE);
  }
}


void dispose_node(np_t *np){
  if (np->status!='W') {
    free(np->linkarray);
    free(np->bdryparam);
    free(np->bs);
  }
  if (is_node_resumed(*np)) {
#ifdef ZONETHREADS
    thread_lock_destroy(&(np->wk->lock),THREADTYPE_ZONE);
#endif
    free(np->wk);
  }
  np->status='D';
  np->FLUIDPRIMMEM=FALSE;
}


void find_numerical_jacobian(np_t *np, long l, gl_t *gl, void(*FUNCT)(np_t, gl_t *, long, flux_t), long theta, sqmat_t Ak){
 long row,col,flux;
 flux_t Fstar,Fstarnew,Ustar,dUstar,dUstar2;

 for (flux=0; flux<nf; flux++){
   Ustar[flux]=np[l].bs->U[flux]*_Omega(np[l],gl);
//   dUstar[flux]=Ustar[flux]/1000000.0e0+1.0e-10;
//   dUstar[flux]=Ustar[flux]/1000000.0e0+1.0e-15;
   dUstar[flux]=gl->cycle.fluid.Uref[flux]/1.0e6*_Omega(np[l],gl);
 }
 (*FUNCT)(np[l],gl,theta,Fstar);
 for (col=0; col<nf; col++){
   for (flux=0; flux<nf; flux++){
     np[l].bs->U[flux]=Ustar[flux]/_Omega(np[l],gl);
     dUstar2[flux]=0.0e0;
   }
   dUstar2[col]=dUstar[col];
   add_dUstar_to_U(np,l,gl,dUstar2);
   (*FUNCT)(np[l], gl, theta, Fstarnew);
   for (row=0; row<nf; row++){
     assert(dUstar[col]!=0.0e0);
     Ak[row][col]=(Fstarnew[row]-Fstar[row])/dUstar[col];
   }
 }
 for (flux=0; flux<nf; flux++){
   np[l].bs->U[flux]=Ustar[flux]/_Omega(np[l],gl);
 }
 find_prim_fluid(np,l,gl);
}


void find_numerical_jacobian_2(np_t *np, long l, gl_t *gl, void(*FUNCT)(np_t, gl_t *, flux_t), sqmat_t Ak){
 long row,col,flux;
 flux_t Fstar,Fstarnew,U,dU,dUstar;

 for (flux=0; flux<nf; flux++){
   U[flux]=np[l].bs->U[flux];
   dU[flux]=gl->cycle.fluid.Uref[flux]/1e6;

 }
 (*FUNCT)(np[l], gl, Fstar);
 for (col=0; col<nf; col++){
   for (flux=0; flux<nf; flux++){
     np[l].bs->U[flux]=U[flux];
     dUstar[flux]=0.0e0;
   }
   dUstar[col]=dU[col]*_Omega(np[l],gl);
   add_dUstar_to_U(np,l,gl,dUstar);
   (*FUNCT)(np[l], gl, Fstarnew);
   for (row=0; row<nf; row++){
     assert(dU[col]!=0.0e0);
     Ak[row][col]=(Fstarnew[row]-Fstar[row])/(dU[col]*_Omega(np[l],gl));
   }
 }
 for (flux=0; flux<nf; flux++){
   np[l].bs->U[flux]=U[flux];
 }
 find_prim_fluid(np,l,gl);
}





void find_numerical_jacobian_3(np_t *np, long l, gl_t *gl, void(*FUNCT)(np_t *, gl_t *, long, flux_t), sqmat_t Ak){
 long row,col,flux;
 flux_t Fstar,Fstarnew,Ustar,dUstar,dUstar2;

 for (flux=0; flux<nf; flux++){
   Ustar[flux]=np[l].bs->U[flux]*_Omega(np[l],gl);
//   dUstar[flux]=Ustar[flux]/10000000.0e0+1.0e-10;
   //dUstar[flux]=Ustar[flux]/10000000.0e0+1.0e-24;
   dUstar[flux]=gl->cycle.fluid.Uref[flux]/1e6*_Omega(np[l],gl);

   if (dUstar[flux]==0.0e0) {
     dUstar[flux]=1.0e-30;
   }
 }
 (*FUNCT)(np,gl,l,Fstar);
 for (col=0; col<nf; col++){
   for (flux=0; flux<nf; flux++){
     np[l].bs->U[flux]=Ustar[flux]/_Omega(np[l],gl);
     dUstar2[flux]=0.0e0;
   }
   dUstar2[col]=dUstar[col];
   add_dUstar_to_U(np,l,gl,dUstar2);
   (*FUNCT)(np, gl, l, Fstarnew);
   for (row=0; row<nf; row++){
     assert(dUstar[col]!=0.0e0);
     Ak[row][col]=(Fstarnew[row]-Fstar[row])/dUstar[col];
   }
 }
 for (flux=0; flux<nf; flux++){
   np[l].bs->U[flux]=Ustar[flux]/_Omega(np[l],gl);
 }
 find_prim_fluid(np,l,gl);
}


/* find the Jacobian A of the flux F such that A*U=F */
void find_homogeneous_jacobian(np_t *np, gl_t *gl, long l, flux_t F, sqmat_t A){
 long row,col;
 flux_t LU,LF;
 sqmat_t L,mattmp,Lambda,Linv;
  long theta;
  jacvars_t jacvars;
  metrics_t metrics;

 /* the idea here is
    A*U=F
    Linv*Lambda*LU=F
    Lambda*LU=L*F
    Lambda[k]=(L*F)[k]/LU[k]

    A=Linv*Lambda*L
 */

  theta=0;
  find_metrics_at_node(np, gl, l, theta, &metrics);
  find_jacvars(np[l], gl, metrics, theta, &jacvars);
  find_Linv_from_jacvars(jacvars, metrics, Linv);
  find_L_from_jacvars(jacvars, metrics, L);
  find_LUstar_from_jacvars(jacvars, metrics, LU);
  multiply_matrix_and_vector(L,F,LF);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Lambda[row][col]=0.0;
    }
    if (LU[row]<=0.0) fatal_error("Characteristic variable can not be less or equal to zero in function find_homogeneous_jacobian().");
    Lambda[row][row]=LF[row]/LU[row];
  }
  multiply_diagonal_matrix_and_matrix(Lambda,L,mattmp);
  multiply_matrix_and_matrix(Linv,mattmp,A);

  /* check */
/* 
  flux_t fluxtmp;
  multiply_matrix_and_vector(A,np[l].bs->U,fluxtmp);
  for (row=0; row<nf; row++){
    printf("%E  %E\n",fluxtmp[row]*_Omega(np[l],gl),F[row]);
  } */
}







bool is_node_inner(np_t np, int nodetype){
  bool tmp;
  long type;
  if (is_node_valid(np,nodetype)) {
    /* need to initialize type to avoid compiler warning */
    type=np.type;
    if (nodetype==TYPELEVEL_FLUID) {
      type=np.type;
    }
    if (nodetype==TYPELEVEL_FLUID_WORK) {
      type=np.type_wk;
    }
#ifdef EMFIELD
    if (nodetype==TYPELEVEL_EMFIELD) {
      type=np.type_emf;
    }
#endif
    if (type==NODETYPE_INNER) {
      tmp=TRUE;
    } else {
      tmp=FALSE;
    }
  } else {
    tmp=FALSE;
  }
  return(tmp);
}


bool is_node_bdry(np_t np, int TYPELEVEL){
  bool tmp;
  long type;
  if (is_node_valid(np,TYPELEVEL)) {
    /* need to initialize type to avoid compiler warning */
    type=np.type;
    if (TYPELEVEL==TYPELEVEL_FLUID) {
      type=np.type;
    }
    if (TYPELEVEL==TYPELEVEL_FLUID_WORK) {
      type=np.type_wk;
    }
#ifdef EMFIELD
    if (TYPELEVEL==TYPELEVEL_EMFIELD) {
      type=np.type_emf;
    }
#endif
    if (type>=NODETYPE_BDRY) {
      tmp=TRUE;
    } else {
      tmp=FALSE;
    }
  } else {
    tmp=FALSE;
  }
  return(tmp);
}


bool is_node_link(np_t np, int TYPELEVEL){
  bool tmp;
  tmp=FALSE;
  if (np.numlink>0 && (TYPELEVEL==TYPELEVEL_FLUID || TYPELEVEL==TYPELEVEL_FLUID_WORK)) tmp=TRUE; 
#ifdef EMFIELD
  if (np.numlink_emf>0 && TYPELEVEL==TYPELEVEL_EMFIELD) tmp=TRUE; 
#endif
  return(tmp);
}


bool is_node_inner_near_bdry(np_t *np, gl_t *gl, long l, int TYPELEVEL){
  bool NEARBDRY;
  long dim;
  NEARBDRY=FALSE;
  if (is_node_inner(np[l],TYPELEVEL)){
    for (dim=0; dim<nd; dim++){
      if (is_node_bdry(np[_al(gl,l,dim,+1)],TYPELEVEL)) NEARBDRY=TRUE;
      if (is_node_bdry(np[_al(gl,l,dim,-1)],TYPELEVEL)) NEARBDRY=TRUE;
    }
  }
  return(NEARBDRY);
}


/* n can be set to 0,1,2,3 */
bool is_node_inner_within_n_nodes_of_bdry(np_t *np, gl_t *gl, long l, long n, int TYPELEVEL){
  bool NEARBDRY;
  long dim,cnt;
  NEARBDRY=FALSE;
  assert(n>=0);
  assert(n<=3);
  if (is_node_inner(np[l],TYPELEVEL)){
    for (dim=0; dim<nd; dim++){
      for (cnt=1; cnt<=n; cnt++){
        if (is_node_bdry(np[_al(gl,l,dim,+cnt)],TYPELEVEL)) NEARBDRY=TRUE;
        if (is_node_bdry(np[_al(gl,l,dim,-cnt)],TYPELEVEL)) NEARBDRY=TRUE;
      }
    }
  }
  return(NEARBDRY);
}


bool is_node_in_zone(long i, long j, long k, zone_t zone){
  bool tmp;
  if (   i>=zone.is && i<=zone.ie
#ifdef _2DL
        && j>=zone.js && j<=zone.je
#endif
#ifdef _3DL
        && k>=zone.ks && k<=zone.ke
#endif
     ) tmp=TRUE; else tmp=FALSE;
  return(tmp);
}


bool is_node_in_zone_2(long l, gl_t *gl, zone_t zone){
  bool tmp;
  long i,j,k;
  find_ijk_from_l(gl, l, &i, &j, &k);

  if (   i>=zone.is && i<=zone.ie
#ifdef _2DL
        && j>=zone.js && j<=zone.je
#endif
#ifdef _3DL
        && k>=zone.ks && k<=zone.ke
#endif
     ) tmp=TRUE; else tmp=FALSE;
  return(tmp);
}


/* returns TRUE if zone1 is enclosed within zone2, FALSE otherwise */
bool is_zone_in_zone(zone_t zone1, zone_t zone2){
  bool RET;
  RET=TRUE;
  if (!is_node_in_zone(zone1.is,zone1.js,zone1.ks,zone2)) RET=FALSE;
  if (!is_node_in_zone(zone1.ie,zone1.je,zone1.ke,zone2)) RET=FALSE;
  return(RET);
}


/* returns TRUE if zone1 intersects zone2, FALSE otherwise */
bool is_zone_intersecting_zone(zone_t zone1, zone_t zone2){
  bool FOUND;
  FOUND=FALSE;
  if (!FOUND && is_node_in_zone(zone1.is,zone1.js,zone1.ks,zone2)) FOUND=TRUE;
  if (!FOUND && is_node_in_zone(zone1.is,zone1.je,zone1.ks,zone2)) FOUND=TRUE;
  if (!FOUND && is_node_in_zone(zone1.ie,zone1.js,zone1.ks,zone2)) FOUND=TRUE;
  if (!FOUND && is_node_in_zone(zone1.ie,zone1.je,zone1.ks,zone2)) FOUND=TRUE;
#ifdef _3DL
  if (!FOUND && is_node_in_zone(zone1.is,zone1.js,zone1.ke,zone2)) FOUND=TRUE;
  if (!FOUND && is_node_in_zone(zone1.is,zone1.je,zone1.ke,zone2)) FOUND=TRUE;
  if (!FOUND && is_node_in_zone(zone1.ie,zone1.js,zone1.ke,zone2)) FOUND=TRUE;
  if (!FOUND && is_node_in_zone(zone1.ie,zone1.je,zone1.ke,zone2)) FOUND=TRUE;
#endif
  return(FOUND);
}


void find_subzones_in_zone_given_zonelength(long zonelength, zone_t zone, long *numzone, zone_t **subzones){
  long cntzone,zonei,numzonei,numzonej,numzonek;
#ifdef _2DL
  long zonej;
#endif
#ifdef _3DL
  long zonek;
#endif
  numzonei=max(round((zone.ie-zone.is)/(double)zonelength),1);
#ifdef _2DL
  numzonej=max(round((zone.je-zone.js)/(double)zonelength),1);
#else
  numzonej=1;
#endif
#ifdef _3DL
  numzonek=max(round((zone.ke-zone.ks)/(double)zonelength),1);
#else
  numzonek=1;
#endif
  *numzone=numzonei*numzonej*numzonek;
  *subzones=(zone_t *)realloc(*subzones,(*numzone)*sizeof(zone_t));
  for (cntzone=0; cntzone<(*numzone); cntzone++){
    zonei=mod(cntzone,numzonei)+1;
    (*subzones)[cntzone].is=zone.is+(zone.ie-zone.is+1)*(zonei-1)/numzonei+min(zonei-1,1);
    (*subzones)[cntzone].ie=min(zone.ie,zone.is+(zone.ie-zone.is+1)*(zonei)/numzonei);
#ifdef _2DL
    zonej=mod(cntzone/(numzonei),numzonej)+1;
    (*subzones)[cntzone].js=zone.js+(zone.je-zone.js+1)*(zonej-1)/numzonej+min(zonej-1,1);
    (*subzones)[cntzone].je=min(zone.je,zone.js+(zone.je-zone.js+1)*(zonej)/numzonej);
#endif
#ifdef _3DL
    zonek=mod(cntzone/(numzonei*numzonej),numzonek)+1;
    (*subzones)[cntzone].ks=zone.ks+(zone.ke-zone.ks+1)*(zonek-1)/numzonek+min(zonek-1,1);
    (*subzones)[cntzone].ke=min(zone.ke,zone.ks+(zone.ke-zone.ks+1)*(zonek)/numzonek);
#endif
  }
}


void find_subzones_in_zone_given_numsubzone(zone_t zone, long numsubzonedesired, long *numsubzone, zone_t **subzones){
  long zonelength,zonelengthmax,numsubzonelocal,zonelengthlocal;
  zonelengthmax=max(zone.ie-zone.is+1,zone.je-zone.js+1);
  bool VALID;
#ifdef _3DL
  zonelengthmax=max(zonelengthmax,zone.ke-zone.ks+1);
#endif
  *numsubzone=0;
  VALID=FALSE;
  zonelength=1;
  for (zonelengthlocal=3; zonelengthlocal<zonelengthmax; zonelengthlocal++){
    find_subzones_in_zone_given_zonelength(zonelengthlocal, zone, &numsubzonelocal, subzones);
    if (fabs(numsubzonedesired-numsubzonelocal)<fabs(numsubzonedesired-*numsubzone)){
      zonelength=zonelengthlocal;
      *numsubzone=numsubzonelocal;
    }
    if (numsubzonelocal==1){
      VALID=TRUE;
    }
  }
  find_subzones_in_zone_given_zonelength(zonelength, zone, numsubzone, subzones);
  if (!VALID) fatal_error("Problem with algorithm in find_subzones_in_zone_given_numsubzone().");
}



bool is_node_in_domain_lim(long l, gl_t *gl){
  bool tmp;
  long i,j,k;
  find_ijk_from_l(gl, l, &i, &j, &k);
  if (   i>=gl->domain_lim.is && i<=gl->domain_lim.ie
#ifdef _2DL
        && j>=gl->domain_lim.js && j<=gl->domain_lim.je
#endif
#ifdef _3DL
        && k>=gl->domain_lim.ks && k<=gl->domain_lim.ke
#endif
     ) tmp=TRUE; else tmp=FALSE;
  return(tmp);
}


bool is_node_valid(np_t np, int TYPELEVEL){
  bool tmp;
  long type;
  /* need to initialize type to avoid compiler warning */
  type=np.type;
  if (TYPELEVEL==TYPELEVEL_FLUID) {
    type=np.type;
  }
  if (TYPELEVEL==TYPELEVEL_FLUID_WORK) {
    type=np.type_wk;
  }
#ifdef EMFIELD
  if (TYPELEVEL==TYPELEVEL_EMFIELD) {
    type=np.type_emf;
  }
#endif
  if (type>=NODETYPE_INNER) {
    tmp=TRUE;
  } else {
    tmp=FALSE;
  }
  return(tmp);
}


long _node_link(np_t np, long cntlink, int TYPELEVEL){
  long tmp;
  /* by default, set to no link */
  tmp=LINK_NONE;
  if (TYPELEVEL==TYPELEVEL_FLUID || TYPELEVEL==TYPELEVEL_FLUID_WORK) {
    assert(cntlink<np.numlink);
    tmp=np.linkarray[cntlink*2];
  }
#ifdef EMFIELD
  if (TYPELEVEL==TYPELEVEL_EMFIELD) {
    assert(cntlink<np.numlink_emf);
    tmp=np.linkarray_emf[cntlink*2];
  }
#endif
  return(tmp);
}

/*
long _node_link_old(np_t np, long cntlink, int TYPELEVEL){
  long tmp;
//  assert(is_node_bdry(np,TYPELEVEL));
  tmp=LINK_NONE;
  if (TYPELEVEL==TYPELEVEL_FLUID || TYPELEVEL==TYPELEVEL_FLUID_WORK) {
    tmp=np.link;
  }
#ifdef EMFIELD
  if (TYPELEVEL==TYPELEVEL_EMFIELD) {
    tmp=np.link_emf;
  }
#endif
  return(tmp);
}
*/

long _num_node_link(np_t np, int TYPELEVEL){
  long tmp;
  tmp=0;
  if (TYPELEVEL==TYPELEVEL_FLUID || TYPELEVEL==TYPELEVEL_FLUID_WORK) {
    tmp=np.numlink;
  }
#ifdef EMFIELD
  if (TYPELEVEL==TYPELEVEL_EMFIELD) {
    tmp=np.numlink_emf;
  }
#endif
  return(tmp);
}


long _node_type(np_t np, int TYPELEVEL){
  long tmp;
  /* need to initialize tmp to avoid compiler warning */
  tmp=np.type;
  if (TYPELEVEL==TYPELEVEL_FLUID) {
    tmp=np.type;
  }
#ifdef EMFIELD
  if (TYPELEVEL==TYPELEVEL_EMFIELD) {
    tmp=np.type_emf;
  }
#endif
  if (TYPELEVEL==TYPELEVEL_FLUID_WORK) {
    if (is_node_resumed(np)) {
      tmp=np.type_wk;
    } else {
      tmp=np.type;
    }
  }
  return(tmp);
}


int find_zone_intersection(zone_t zone1, zone_t zone2, zone_t *zoneint){
  int retval;
  retval=0;
  zoneint->is=max(zone2.is,zone1.is);
  zoneint->ie=min(zone2.ie,zone1.ie);
  if (zoneint->is>zoneint->ie) retval=1;
#ifdef _2DL
  zoneint->js=max(zone2.js,zone1.js);
  zoneint->je=min(zone2.je,zone1.je);
  if (zoneint->js>zoneint->je) retval=2;
#endif
#ifdef _3DL
  zoneint->ks=max(zone2.ks,zone1.ks);
  zoneint->ke=min(zone2.ke,zone1.ke);
  if (zoneint->ks>zoneint->ke) retval=3;
#endif
  return(retval);
}


zone_t _zone_from_point(long i, long j, long k){
  zone_t zone;
  zone.is=i;   zone.ie=i;
  zone.js=j;   zone.je=j;
  zone.ks=k;   zone.ke=k;
  return(zone);
}


zone_t _zone_expansion(zone_t zone, long expnum){
  zone.is-=expnum;
  zone.ie+=expnum;
#ifdef _2DL
  zone.js-=expnum;
  zone.je+=expnum;
#endif
#ifdef _3DL
  zone.ks-=expnum;
  zone.ke+=expnum;
#endif
  return(zone);
}


zone_t _zone_intersection(zone_t zone1, zone_t zone2){
  zone_t zone;
  int PROBLEM;
  PROBLEM=find_zone_intersection(zone1,zone2,&zone);
#ifdef _2D
  if (PROBLEM)   fatal_error("In _zone_intersection(), can't find intersection between zone1 [%ld,%ld <-> %ld,%ld] and zone 2 [%ld,%ld <-> %ld,%ld]. This is probably due to some problem within Grid() or within Bdry().",zone1.is,zone1.js,zone1.ie,zone1.je,zone2.is,zone2.js,zone2.ie,zone2.je);
#endif

#ifdef _3D
  if (PROBLEM)   fatal_error("In _zone_intersection(), can't find intersection between zone1 [%ld,%ld,%ld <-> %ld,%ld,%ld] and zone 2 [%ld,%ld,%ld <-> %ld,%ld,%ld]. This is probably due to some problem within Grid() or within Bdry().",zone1.is,zone1.js,zone1.ks,zone1.ie,zone1.je,zone1.ke,zone2.is,zone2.js,zone2.ks,zone2.ie,zone2.je,zone2.ke);
#endif

  return(zone);
}





/* was used in src/post.c only */
zone_t _zone_intersection_experimental(zone_t zone1, zone_t zone2){
  zone1.is=min(max(zone2.is,zone1.is),min(zone2.ie,zone1.ie));
  zone1.ie=max(max(zone2.is,zone1.is),min(zone2.ie,zone1.ie));
#ifdef _2DL
  zone1.js=min(max(zone2.js,zone1.js),min(zone2.je,zone1.je));
  zone1.je=max(max(zone2.js,zone1.js),min(zone2.je,zone1.je));
#endif
#ifdef _3DL
  zone1.ks=min(max(zone2.ks,zone1.ks),min(zone2.ke,zone1.ke));
  zone1.ke=max(max(zone2.ks,zone1.ks),min(zone2.ke,zone1.ke));
#endif
  return(zone1);
}



void validate_data_structure(np_t *np,gl_t *gl){
  long i,j,k,i2,j2,k2;
  for_1DL(i,gl->domain_lim.is,gl->domain_lim.ie){
    for_2DL(j,gl->domain_lim.js,gl->domain_lim.je){
      for_3DL(k,gl->domain_lim.ks,gl->domain_lim.ke){
         find_ijk_from_l(gl, _ai(gl,i,j,k), &i2, &j2, &k2);
#ifndef NDEBUG
         np[_ai(gl,i,j,k)].i=i;
         np[_ai(gl,i,j,k)].j=j;
         np[_ai(gl,i,j,k)].k=k;
#endif
         if (i!=i2 if2DL( || j!=j2) if3DL( || k!=k2) )
#ifdef _2D
           fatal_error("Problem in the data structure at i=%ld j=%ld "
                      "(i2=%ld j2=%ld).",i,j,i2,j2);
#endif
#ifdef _3D
           fatal_error("Problem in the data structure at i=%ld j=%ld k=%ld "
                      "(i2=%ld j2=%ld k2=%ld).",i,j,k,i2,j2,k2);
#endif

      }
    }
  }
}


void init_data_structure(np_t **np, gl_t *gl, zone_t domain, zone_t domain_all){
  gl->domain=domain;
  gl->domain_all=domain_all;
  gl->window=gl->domain;
  gl->domain_lim=_domain_lim_from_domain(gl->domain,gl);
  gl->domain_lim_all=_domain_lim_from_domain(gl->domain_all,gl);
  *np = (np_t *) malloc((gl->domain_lim.ie-gl->domain_lim.is+1)
#ifdef _2DL
                         *(gl->domain_lim.je-gl->domain_lim.js+1)
#endif
#ifdef _3DL
                         *(gl->domain_lim.ke-gl->domain_lim.ks+1)
#endif
                 *sizeof(np_t));
  if (*np==NULL) {
    perror("malloc");
    fatal_error("Problem malloc'ing *np in create_domain_of_cuts_along_x subroutine.");
  }
  validate_data_structure(*np,gl);
}


void init_data_structure_and_create_nodes(np_t **np, gl_t *gl, zone_t domain, zone_t domain_all){
  long i,j,k;
  init_data_structure(np,gl,domain,domain_all);
  for_1DL (i,gl->domain_lim.is,gl->domain_lim.ie){
    for_2DL (j,gl->domain_lim.js,gl->domain_lim.je){
      for_3DL (k,gl->domain_lim.ks,gl->domain_lim.ke){
        create_node(&((*np)[_ai(gl,i,j,k)]), i, j, k);
        suspend_node(&((*np)[_ai(gl,i,j,k)]));
      }
    }
  }
}



/* numxnode is only used here, in this function */
//#define numxnode max(2,hbw_res_fluid-1)
#define numxnode 1


/* extend domain by extra nodes necessary for computation */ 
zone_t _domain_lim_from_domain(zone_t domain, gl_t *gl){
  zone_t zone;
  long tmp;

  tmp=max(numxnode,numxnode+nodeoverlap);
  zone.is=max(gl->domain_all.is-numxnode, domain.is-tmp);
  zone.ie=min(gl->domain_all.ie+numxnode, domain.ie+tmp);
  /* need to initialize zone to avoid compiler warning */
  zone.js=0;
  zone.je=0;
  zone.ks=0;
  zone.ke=0;
  #ifdef _2DL
    zone.js=max(gl->domain_all.js-numxnode, domain.js-tmp);
    zone.je=min(gl->domain_all.je+numxnode, domain.je+tmp);
  #endif
  #ifdef _3DL
    zone.ks=max(gl->domain_all.ks-numxnode, domain.ks-tmp);
    zone.ke=min(gl->domain_all.ke+numxnode, domain.ke+tmp);
  #endif
  return(zone);
}


#ifdef DISTMPI

int _node_rank(gl_t *gl, long i, long j, long k){
  int rank;
  zone_t domainlim;


  if (!gl->DISTDOMAIN) return(0);
  domainlim=_domain_lim_from_domain(gl->domain_all, gl);

#ifdef _2D
  rank=(i-domainlim.is)*gl->numdomain_i/(domainlim.ie-domainlim.is+2)
      +gl->numdomain_i*((j-domainlim.js)*gl->numdomain_j/(domainlim.je-domainlim.js+2));
#endif
#ifdef _3D
  rank=(i-domainlim.is)*gl->numdomain_i/(domainlim.ie-domainlim.is+2)
      +gl->numdomain_i*((j-domainlim.js)*gl->numdomain_j/(domainlim.je-domainlim.js+2))
      +gl->numdomain_i*gl->numdomain_j*((k-domainlim.ks)*gl->numdomain_k/(domainlim.ke-domainlim.ks+2));
#endif

  return(rank);
}


zone_t _domain_from_rank_mem(int rank, gl_t *gl){
  zone_t zone;
  long i,j,k;

  if (!gl->DISTDOMAIN) {
    fatal_error("gl->DISTDOMAIN must be set to TRUE when calling _domain_from_rank_mem().");
  }
  zone.is=gl->domain_all.ie;
  zone.ie=gl->domain_all.is;
  zone.js=gl->domain_all.je;
  zone.je=gl->domain_all.js;
  zone.ks=gl->domain_all.ke;
  zone.ke=gl->domain_all.ks;
  for_1DL(i,gl->domain_all.is,gl->domain_all.ie){
    for_2DL(j,gl->domain_all.js,gl->domain_all.je){
      for_3DL(k,gl->domain_all.ks,gl->domain_all.ke){
        if (_node_rank(gl, i, j, k)==rank){
          zone.is=min(zone.is,i);
          zone.ie=max(zone.ie,i);
          zone.js=min(zone.js,j);
          zone.je=max(zone.je,j);
          zone.ks=min(zone.ks,k);
          zone.ke=max(zone.ke,k);
        }
      }
    }
  }
  if (zone.is>zone.ie || zone.js>zone.je || zone.ks>zone.ke){
    fatal_error("Problem in _domain_from_rank_mem() numdomain_i=%ld numdomain_j=%ld numdomain_k=%ld rank=%d zone.is=%ld .ie=%ld .js=%ld .je=%ld .ks=%ld .ke=%ld .",gl->numdomain_i,gl->numdomain_j,gl->numdomain_k,rank,zone.is,zone.ie,zone.js,zone.je,zone.ks,zone.ke);
  }
  return(zone);
  
}



zone_t _domain_from_rank(int rank, gl_t *gl){
  return(gl->domain_from_rank[rank]);
}


zone_t _domain_lim_from_rank(int rank, gl_t *gl){
  return(gl->domain_lim_from_rank[rank]);
}


zone_t _domain_lim_from_rank_mem(int rank, gl_t *gl){
  zone_t zone,domain;
  domain=_domain_from_rank_mem(rank,gl);
  zone=_domain_lim_from_domain(domain,gl);
  return(zone);
}

static long MPI_IBsend_cnt;

int MPI_IBsend(void* message, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm){
  int retval;
  MPI_IBsend_cnt++;
  if (mod(MPI_IBsend_cnt,500)==0 && FALSE)
    retval=MPI_Ssend(message,count,datatype,dest,tag,comm);
  else retval=MPI_Send(message,count,datatype,dest,tag,comm);
  return(retval);
}

/* returns: MPI_NO_DATA_EXCHANGED if no exchange takes place
            MPI_DATA_SENT  if process sends packet to another process
            MPI_DATA_RECEIVED if process receives packet from another process */
int MPI_Bcast_Node(void *message, int count, MPI_Datatype datatype, int root,
                   MPI_Comm comm, long i, long j, long k, gl_t *gl){
  int ret;
  int rank,numproc,proc;
  MPI_Status MPI_Status1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  ret=MPI_NO_DATA_EXCHANGED;
  if (rank==root) {
    for (proc=0; proc<numproc; proc++){
      if (proc!=root){
        if (is_node_in_zone(i,j,k,_domain_lim_from_rank(proc,gl))){
          MPI_IBsend(message,count,datatype,proc,0,MPI_COMM_WORLD);
          ret=MPI_DATA_SENT;
        }
      }
    }
  } else {
    if (is_node_in_zone(i,j,k,gl->domain_lim)) {
      MPI_Recv(message,count,datatype,root,0,MPI_COMM_WORLD,&MPI_Status1);
      ret=MPI_DATA_RECEIVED;
    }
  }

  /* ret=MPI_Bcast(message,count,datatype,root,comm); */
  return(ret);
}

#endif



/* Kronecker delta */
double _delta(long r, long k){
  double delta;
  if (r==k) delta=1.0; else delta=0.0;
  return(delta);
}


void add_int_to_codex(SOAP_codex_t *codex, char *varname, long varvalue){
  char *tmpstr;
  tmpstr=(char *)malloc(50*sizeof(char));
  sprintf(tmpstr,"%ld",varvalue);
  SOAP_add_to_vars(codex,varname,tmpstr);
  free(tmpstr);
}

void add_string_to_codex(SOAP_codex_t *codex, char *varname, char *varvalue){
  char *str;
  str=(char *)malloc((strlen(varvalue)+2)*sizeof(char));
  strcpy(str,varvalue);
  SOAP_strins("\"", &str,  0);
  SOAP_strins("\"", &str,  strlen(str));
  SOAP_add_to_vars(codex,varname,str);
  free(str);
}


void add_double_to_codex(SOAP_codex_t *codex, char *varname, double varvalue){
  char *tmpstr;
  tmpstr=(char *)malloc(50*sizeof(char));
  sprintf(tmpstr,"%20.15E",varvalue);
  SOAP_add_to_vars(codex,varname,tmpstr);
  free(tmpstr);
}


void find_double_var_from_codex(SOAP_codex_t *codex, char *name, double *var){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (!gl->CONTROL_READ){
    *var=SOAP_var_value(codex, name);
  } else {
    if (SOAP_is_var_in_codex(codex, name)){
      *var=SOAP_var_value(codex, name);
    }
  }
}


void find_bool_var_from_codex(SOAP_codex_t *codex, char *name, bool *var){
  double vardouble;
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (!gl->CONTROL_READ){
    vardouble=SOAP_var_value(codex, name);
    if (!SOAP_is_double_a_long(vardouble) || !(round(vardouble)==0 || round(vardouble)==1)) 
      SOAP_fatal_error(codex,"Value given to variable \"%s\" must be of boolean type.",name);
    *var=round(vardouble);
  } else {
    if (SOAP_is_var_in_codex(codex, name)){
      vardouble=SOAP_var_value(codex, name);
      if (!SOAP_is_double_a_long(vardouble) || !(round(vardouble)==0 || round(vardouble)==1)) 
        SOAP_fatal_error(codex,"Value given to variable \"%s\" must be of boolean type.",name);
      *var=round(vardouble);
    }
  }
}


void find_int_var_from_codex(SOAP_codex_t *codex, char *name, int *var){
  double vardouble;
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (!gl->CONTROL_READ){
    vardouble=SOAP_var_value(codex, name);
    if (!SOAP_is_double_a_long(vardouble)) 
      SOAP_fatal_error(codex,"Value given to variable \"%s\" must be of integer type.",name);
   *var=round(vardouble);
  } else {
    if (SOAP_is_var_in_codex(codex, name)){
      vardouble=SOAP_var_value(codex, name);
      if (!SOAP_is_double_a_long(vardouble)) 
        SOAP_fatal_error(codex,"Value given to variable \"%s\" must be of integer type.",name);
      *var=round(vardouble);
    }
  }
}



void find_zone_from_argum(char *argum, long startpos, gl_t *gl, SOAP_codex_t *codex, zone_t *zone){
  zone->is=SOAP_get_argum_long(codex,argum,startpos);
  zone->ie=SOAP_get_argum_long(codex,argum,startpos+nd);
  if (zone->ie<zone->is) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be less than is (%ld).",zone->ie,zone->is);
  if (zone->ie>gl->domain_all.ie)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be greater than domain_all.ie (%ld).",zone->ie,gl->domain_all.ie);
  if (zone->is<gl->domain_all.is)
    SOAP_fatal_error(codex,"Zone boundaries invalid: is (%ld) can not be less than domain_all.is (%ld).",zone->is,gl->domain_all.is);
  
#ifdef _2DL
  zone->js=SOAP_get_argum_long(codex,argum,startpos+1);
  zone->je=SOAP_get_argum_long(codex,argum,startpos+1+nd);
  if (zone->je<zone->js) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: je (%ld) can not be less than js (%ld).",zone->je,zone->js);
  if (zone->je>gl->domain_all.je)
    SOAP_fatal_error(codex,"Zone boundaries invalid: je (%ld) can not be greater than domain_all.je (%ld).",zone->je,gl->domain_all.je);
  if (zone->js<gl->domain_all.js)
    SOAP_fatal_error(codex,"Zone boundaries invalid: js (%ld) can not be less than domain_all.js (%ld).",zone->js,gl->domain_all.js);
#endif
#ifdef _3DL
  zone->ks=SOAP_get_argum_long(codex,argum,startpos+2);
  zone->ke=SOAP_get_argum_long(codex,argum,startpos+2+nd);
  if (zone->ke<zone->ks) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: ke (%ld) can not be less than ks (%ld).",zone->ke,zone->ks);
  if (zone->ke>gl->domain_all.ke)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ke (%ld) can not be greater than domain_all.ke (%ld).",zone->ke,gl->domain_all.ke);
  if (zone->ks<gl->domain_all.ks)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ks (%ld) can not be less than domain_all.ks (%ld).",zone->ks,gl->domain_all.ks);
#endif
}


void thread_lock_destroy(thread_lock_t *lockvar, int threadtype){
  switch (threadtype){
    case THREADTYPE_ALL:
      #ifdef PTHREADS
        pthread_mutex_destroy(lockvar);
      #endif
      #ifdef OPENMPTHREADS
        omp_destroy_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_destroy(lockvar);
      #endif
      #ifdef OPENMPTHREADS
        omp_destroy_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_destroy(lockvar);
      #endif
    break;   
    case THREADTYPE_OPENMP:
      #ifdef OPENMPTHREADS
        omp_destroy_lock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX:
      #ifdef PTHREADS
        pthread_mutex_destroy(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_destroy(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_destroy(lockvar);
      #endif
    break;
  }
}


void thread_lock_init(thread_lock_t *lockvar, int threadtype){
  switch (threadtype){
    case THREADTYPE_ALL:
      #ifdef PTHREADS
        pthread_mutex_init(lockvar, NULL);
      #endif
      #ifdef OPENMPTHREADS
        omp_init_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_init(lockvar, NULL);
      #endif
      #ifdef OPENMPTHREADS
        omp_init_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_init(lockvar, NULL);
      #endif
    break;   
    case THREADTYPE_OPENMP:
      #ifdef OPENMPTHREADS
        omp_init_lock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX:
      #ifdef PTHREADS
        pthread_mutex_init(lockvar, NULL);
      #endif
    break;
    case THREADTYPE_POSIX_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_init(lockvar, NULL);
      #endif
    break;
    case THREADTYPE_POSIX_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_init(lockvar, NULL);
      #endif
    break;
  }
}


void thread_lock_set(thread_lock_t *lockvar, int threadtype){
  switch (threadtype){
    case THREADTYPE_ALL:
      #ifdef PTHREADS
        pthread_mutex_lock(lockvar);
      #endif
      #ifdef OPENMPTHREADS
        omp_set_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_lock(lockvar);
      #endif
      #ifdef OPENMPTHREADS
        omp_set_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_OPENMP:
      #ifdef OPENMPTHREADS
        omp_set_lock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX:
      #ifdef PTHREADS
        pthread_mutex_lock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_lock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_lock(lockvar);
      #endif
    break;
  }
}


void thread_lock_unset(thread_lock_t *lockvar, int threadtype){
  switch (threadtype){
    case THREADTYPE_ALL:
      #ifdef PTHREADS
        pthread_mutex_unlock(lockvar);
      #endif
      #ifdef OPENMPTHREADS
        omp_unset_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_unlock(lockvar);
      #endif
      #ifdef OPENMPTHREADS
        omp_unset_lock(lockvar);
      #endif
    break;   
    case THREADTYPE_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_unlock(lockvar);
      #endif
    break;   
    case THREADTYPE_OPENMP:
      #ifdef OPENMPTHREADS
        omp_unset_lock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX:
      #ifdef PTHREADS
        pthread_mutex_unlock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX_ZONE:
      #ifdef ZONETHREADS
        pthread_mutex_unlock(lockvar);
      #endif
    break;
    case THREADTYPE_POSIX_LOOP:
      #ifdef POSIXTHREADS
        pthread_mutex_unlock(lockvar);
      #endif
    break;
  }
}



void thread_lock_node_set(np_t *np, long l, int threadtype){
  if (threadtype!=THREADTYPE_ZONE && threadtype!=THREADTYPE_POSIX_ZONE) 
     fatal_error("Can't unlock node for threadtype %ld in thread_lock_node_unset.",threadtype);
#ifdef ZONETHREADS
  thread_lock_set(&(np[l].wk->lock),threadtype);
#endif
}


void thread_lock_node_unset(np_t *np, long l, int threadtype){
  if (threadtype!=THREADTYPE_ZONE && threadtype!=THREADTYPE_POSIX_ZONE) 
     fatal_error("Can't unlock node for threadtype %ld in thread_lock_node_unset.",threadtype);
#ifdef ZONETHREADS
  thread_lock_unset(&(np[l].wk->lock),threadtype);
#endif
}


void thread_lock_global_set(gl_t *gl, int threadtype){
  thread_lock_set(&(gl->lock),threadtype);
}


void thread_lock_global_unset(gl_t *gl, int threadtype){
  thread_lock_unset(&(gl->lock),threadtype);
}


double distance2_between_nodes(np_t *np, gl_t *gl, long l1, long l2){
  double dist2;
  long dim;
  dist2=0.0;
  for (dim=0; dim<nd; dim++){
    dist2+=sqr(_x(np[l1],dim)-_x(np[l2],dim));
  }
  return(dist2);
}


void write_options_row ( FILE * outputfile, char *col1, char *col2, char *col3, int linewidth, int lengthcol1, int lengthcol2 ) {
  char *linestr;

  linestr=(char *)malloc(sizeof(char)*(10000+strlen(col1)+strlen(col2)+strlen(col3)+linewidth));
  strcpy ( linestr, "  " );

  strcat ( linestr, col1 );
  while ( strlen ( linestr ) < ( lengthcol1 ) ) {
    strcat ( linestr, " " );
  }

  strcat ( linestr, col2 );
  while ( strlen ( linestr ) < ( lengthcol1 + lengthcol2 ) ) {
    strcat ( linestr, " " );
  }

  strcat ( linestr, col3 );
  wfprintf ( outputfile, strwrpind ( linestr, linewidth, -( lengthcol1 + lengthcol2 ) ) );
  wfprintf ( outputfile, "\n" );
  free(linestr);
}



