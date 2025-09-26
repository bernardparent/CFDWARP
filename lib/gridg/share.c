// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000,2022 Bernard Parent

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

#include "share.h"
#include <stdio.h>
#include <stdarg.h>

#define EOS 0

typedef struct {
  bool FOUND;
  char *actionname;
  void (*action)(char **, void *, SOAP_codex_t *);
  void *action_args;
} RXaction_args2_t;


long GRIDG_ai3(GRIDG_gl3d_t gl, long i, long j, long k) {
  long ii;
  ii=(i-gl.is+1);
  ii=ii*(gl.je-gl.js+3)+(j-gl.js+1);
  ii=ii*(gl.ke-gl.ks+3)+(k-gl.ks+1);
  return(ii);
}


long GRIDG_al3(GRIDG_gl3d_t gl, long l, long dim, long offset) {
  long ii;
  
  ii=l+offset;
  if (dim==0) {
    ii=l+offset*(gl.je-gl.js+3)
               *(gl.ke-gl.ks+3);
  }
  if (dim==1) {
    ii=l+offset*(gl.ke-gl.ks+3);
  }
  
  return(ii);
}


void GRIDG_find_ijk_from_l(GRIDG_gl3d_t gl, long l, long *i, long *j, long *k){
    *k=mod(l,gl.ke-gl.ks+3)+gl.ks-1;
    *j=mod(l/(gl.ke-gl.ks+3),(gl.je-gl.js+3))
          +gl.js-1;
    *i=(l/(gl.ke-gl.ks+3))/(gl.je-gl.js+3)
          +gl.is-1;
}



long GRIDG_ai2(GRIDG_gl2d_t gl, long i, long j) {
  long ii;
  ii=(i-gl.is+1)*(gl.je-gl.js+3)+(j-gl.js+1);
  return(ii);
}

long GRIDG_ai1(GRIDG_gl1d_t gl, long i) {
  long ii;
  ii=(i-gl.is+1);
  return(ii);
}


void read_x_chars_from_file(FILE *infile, long cntmax){
  long cnt;
  char tmp;
  for (cnt=0; cnt<cntmax; cnt++){
    if (fscanf(infile, "%c", &tmp)!=1) GRIDG_fatal_error("The function fscanf could not read properly in read_x_chars_from_file.");
  }

}


int find_root_eps(double *eps, double(*FUNCT)(void *, double), void *arg_ptr, char *calledfunction, SOAP_codex_t *codex){
  long IFLAG;

  *eps=EXM_find_root_zero_in(FUNCT,arg_ptr,0.2,5.0,1.0e-12,1.0e-12,&IFLAG);
  if (IFLAG==4) {
    *eps=EXM_find_root_Newton_Raphson(FUNCT, arg_ptr,  1.0, 0.01, 1.0e-12, 1.0e-12, &IFLAG);
    if (IFLAG!=1 || !isfinite(*eps)) {
#ifdef NDEBUG
      SOAP_fatal_error(codex,"Problem finding root in grid segment.");
#else
      SOAP_fatal_error(codex,"Problem finding root in subroutine %s.",calledfunction);
#endif
    }
  }
  return(0);
}



void find_segment_0(double *seg, long NL, long NH,
                  double dzeta1, double epsL, double epsH, SOAP_codex_t *codex){
  double dzeta,zeta;
  long cntseg,cnt;
  dzeta=dzeta1;
  zeta=0.0e0;
  cntseg=0;
  for (cnt=1; cnt<=NL-1; cnt++){
    seg[cntseg]=zeta;
    cntseg=cntseg+1;
    zeta=zeta+dzeta;
    dzeta=dzeta*epsL;
  }
  dzeta=dzeta/epsL;
  for (cnt=1; cnt<=NH; cnt++){
    seg[cntseg]=zeta;
    cntseg=cntseg+1;
    zeta=zeta+dzeta;
    dzeta=dzeta*epsH;
  }
}

typedef struct {
  double epsL,epsH,L;
  long NH,NL;
} arg1_t;


double segment_1_funct(void *arg1, double dzeta1){
  double error;
  long cnt;
  double sum1,sum2;

  sum1=0.0e0;
  for (cnt=1; cnt<=((arg1_t *)arg1)->NL-1; cnt++){
    sum1=sum1+dzeta1*pow(((arg1_t *)arg1)->epsL,(double)(cnt-1));
  }
  sum2=0.0e0;
  for (cnt=1; cnt<=((arg1_t *)arg1)->NH-1; cnt++){
    sum2=sum2+dzeta1*pow(((arg1_t *)arg1)->epsL,(double)(((arg1_t *)arg1)->NL-2))*
         pow(((arg1_t *)arg1)->epsH,(double)(cnt-1));
  }
  error=-((arg1_t *)arg1)->L+sum1+sum2;
  return(error);
}


/* N is the total number of nodes */
int find_segment_1(double *seg, long N, double NLoverN,
                  double L, double epsL, double epsH, SOAP_codex_t *codex){
  double dzeta1;
  long NL,NH;
  long IFLAG;
  arg1_t arg1;
  int err;

  NL=NLoverN*N;        NH=N-NL+1;
  arg1.NL=NL;          arg1.NH=NH;
  arg1.epsH=epsH;      arg1.epsL=epsL;
  arg1.L=L;

  dzeta1=EXM_find_root_zero_in(&segment_1_funct, &arg1, 0.0e0, L, 1.0e-12, 1.0e-12, &IFLAG);
  if (IFLAG==4) {
#ifdef NDEBUG
    SOAP_fatal_error(codex,"Problem finding root within grid segment.");
#else
    SOAP_fatal_error(codex,"Problem finding root in find_segment_1.");
#endif
    err=1;
  } else {
    err=0;
  }
  find_segment_0(seg, NL, NH, dzeta1, epsL, epsH, codex);
  return(err);
}



typedef struct {
  double dzeta1,epsL,L;
  long N,NH,NL;
} arg2_t;


double segment_2_funct(void *arg2, double epsH){
  double error;
  long cnt;
  double sum1,sum2;

  sum1=0.0e0;
  for (cnt=1; cnt<=((arg2_t *)arg2)->NL-1; cnt++){
    sum1=sum1+((arg2_t *)arg2)->dzeta1*pow(((arg2_t *)arg2)->epsL,(double)(cnt-1));
  }
  sum2=0.0e0;
  for (cnt=1; cnt<=((arg2_t *)arg2)->NH-1; cnt++){
    sum2=sum2+((arg2_t *)arg2)->dzeta1*pow(((arg2_t *)arg2)->epsL,(double)(((arg2_t *)arg2)->NL-2))*
         pow(epsH,(double)(cnt-1));
  }
  error=-((arg2_t *)arg2)->L+sum1+sum2;
  return(error);
}

/* N is the total number of nodes */
int find_segment_2(double *seg, long N, double NLoverN,
                 double L, double dzeta1, double epsL, SOAP_codex_t *codex){
  double epsH;
  long NL,NH;
  arg2_t arg2;
  int err;

  NL=N*NLoverN;        NH=N-NL+1;
  arg2.NL=NL;          arg2.NH=NH;
  arg2.dzeta1=dzeta1;  arg2.epsL=epsL;
  arg2.L=L;            arg2.N=N;

  err=find_root_eps(&epsH,&segment_2_funct, &arg2, "find_segment_2",codex);
  find_segment_0(seg, NL, NH, dzeta1, epsL, epsH, codex);
  return(err);
}


typedef struct {
  double dzeta1,dzeta2,L;
  long N,NH,NL;
} arg3_t;

double segment_3_funct(void *arg3, double epsL){
  double error;
  long cnt;
  double sum1,sum2;

  sum1=0.0e0;
  for (cnt=1; cnt<=((arg3_t *)arg3)->NL-1; cnt++){
    sum1=sum1+((arg3_t *)arg3)->dzeta1*pow(epsL,(double)(cnt-1));
  }
  sum2=0.0e0;
  for (cnt=1; cnt<=((arg3_t *)arg3)->NH-1; cnt++){
    sum2=sum2+((arg3_t *)arg3)->dzeta2
       *pow(((arg3_t *)arg3)->dzeta1*pow(epsL,(double)(((arg3_t *)arg3)->NL-2))/((arg3_t *)arg3)->dzeta2
                        ,(double)(cnt-((arg3_t *)arg3)->NH+1)/(double)(2-((arg3_t *)arg3)->NH));
  }
  error=-((arg3_t *)arg3)->L+sum1+sum2;
  return(error);
}


/* N is the total number of nodes */
int find_segment_3(double *seg, long N, double NLoverN,
                  double L, double dzeta1, double dzeta2, SOAP_codex_t *codex){
  double epsL,epsH;
  long NL,NH;
  int err;
  arg3_t arg3;

  NL=N*NLoverN;        NH=N-NL+1;
  arg3.NL=NL;          arg3.NH=NH;
  arg3.dzeta1=dzeta1;  arg3.dzeta2=dzeta2;
  arg3.L=L;            arg3.N=N;

  /* find epsL */
  err=find_root_eps(&epsL,&segment_3_funct, &arg3, "find_segment_3",codex);
  epsH=pow(dzeta1/dzeta2*pow(epsL,(double)(NL-2)),1.0e0/(double)(2-NH));

  find_segment_0(seg, NL, NH, dzeta1, epsL, epsH, codex);
  return(err);
}


typedef struct {
  double dzeta1,epsH,L;
  long N,NH,NL;
} arg4_t;

double segment_4_funct(void *arg4, double epsL){
  double error;
  long cnt;
  double sum1,sum2;

  sum1=0.0e0;
  for (cnt=1; cnt<=((arg4_t *)arg4)->NL-1; cnt++){
    sum1=sum1+((arg4_t *)arg4)->dzeta1*pow(epsL,(double)(cnt-1));
  }
  sum2=0.0e0;
  for (cnt=1; cnt<=((arg4_t *)arg4)->NH-1; cnt++){
    sum2=sum2+((arg4_t *)arg4)->dzeta1*pow(epsL,(double)(((arg4_t *)arg4)->NL-2))*
         pow(((arg4_t *)arg4)->epsH,(double)(cnt-1));
  }
  error=-((arg4_t *)arg4)->L+sum1+sum2;
  return(error);
}



/* N is the total number of nodes */
int find_segment_4(double *seg, long N, double NLoverN,
                 double L, double dzeta1, double epsH, SOAP_codex_t *codex){
  double epsL;
  long NL,NH;
  int err;
  arg4_t arg4;

  NL=N*NLoverN;        NH=N-NL+1;
  arg4.NL=NL;          arg4.NH=NH;
  arg4.dzeta1=dzeta1;  arg4.epsH=epsH;
  arg4.L=L;            arg4.N=N;

  err=find_root_eps(&epsL,&segment_4_funct, &arg4, "find_segment_4", codex);
  find_segment_0(seg, NL, NH, dzeta1, epsL, epsH, codex);
  return(err);
}


/* N is the total number of nodes */
int find_segment_5(double *seg, long N, double NLoverN,
                 double L, double epsL, double dzeta2, SOAP_codex_t *codex){
   double *seg2;
   long cnt;
   int err;
   seg2= (double *) malloc ((N+4)*sizeof(double));
   err=find_segment_4(seg2, N, 1.0e0-NLoverN, L, dzeta2, 1.0e0/epsL, codex);
   for (cnt=0; cnt<N; cnt++){
     seg[cnt]=L-seg2[N-cnt-1];
   }
   seg[0]=0.0e0;
   free(seg2);
   return(err);
}


/* N is the total number of nodes */
int find_segment_6(double *seg, long N, double NMoverN,
                 double L, double dzeta1, double dzeta2, SOAP_codex_t *codex){
   double *segL,*segH,*segM;
   long cnt,NL,NH,NM,i2,i1;
   int err;
   NL=N/2+1;
   NH=N-NL+1;
   segL= (double *) malloc ((NL+4)*sizeof(double));
   segH= (double *) malloc ((NH+4)*sizeof(double));
   err=find_segment_4(segL, NL, 1.0e0-NMoverN, 0.5e0*L, dzeta1, 1.0e0, codex);
   err+=find_segment_5(segH, NH, NMoverN, 0.5e0*L, 1.0e0, dzeta2, codex);
   for (cnt=0; cnt<NL; cnt++) seg[cnt]=segL[cnt];
   for (cnt=1; cnt<NH; cnt++) seg[cnt+NL-1]=segH[cnt]+0.5e0*L;
   free(segL);
   free(segH);

   i1=(long)(((double)N-(double)N*NMoverN)/2.0e0)+1;
   i2=N-i1;
   NM=i2-i1+1;
   segM= (double *) malloc ((NM+4)*sizeof(double));
   err+=find_segment_3(segM, NM, 0.5e0, seg[i2]-seg[i1],
                seg[i1+1]-seg[i1], seg[i2]-seg[i2-1],codex);
   for (cnt=1; cnt<NM; cnt++) seg[i1+cnt]=segM[cnt]+seg[i1];
   free(segM);
   return(err);
}


/* N is the total number of nodes */
int find_segment_NO(double *seg, long N, double garb0,
                   double L, double garb1, double garb2, SOAP_codex_t *codex){
   if (N==1) {
     seg[0]=0.0e0;
   } else {
     SOAP_fatal_error(codex,"The function find_segment_NO can only be used on a segment with N=1.");
   }
   return(0);
}


void find_x_along_segment(GRIDG_xgrid_t *xgrid, GRIDG_gl1d_t gl, double *seg,
                       long i1, long i2){
  double L;
  long i,N,direc;
  double x1,x2,A1,A2;

  x1=xgrid[GRIDG_ai1(gl,i1)].x;
  A1=xgrid[GRIDG_ai1(gl,i1)].Acs;
  x2=xgrid[GRIDG_ai1(gl,i2)].x;
  A2=xgrid[GRIDG_ai1(gl,i2)].Acs;
  direc=+1;
  if (i2<i1) direc=-1;
  N=labs(i2-i1)+1;
  L=seg[N-1];
  assert(L!=0.0e0);
  for (i=0; i<N; i++){
    xgrid[GRIDG_ai1(gl,i1+direc*i)].x=(double)seg[i]/L*(x2-x1)+x1;
    xgrid[GRIDG_ai1(gl,i1+direc*i)].Acs=(double)seg[i]/L*(A2-A1)+A1;
    xgrid[GRIDG_ai1(gl,i1+direc*i)].INIT=TRUE;
  }
}


void find_xy_along_segment(GRIDG_xygrid_t *xygrid, GRIDG_gl2d_t gl, double *seg,
                        long i1, long j1, long i2, long j2){
  double L;
  long i,j,N,direc;
  double x1,x2,y1,y2;

  x1=xygrid[GRIDG_ai2(gl,i1,j1)].x;
  y1=xygrid[GRIDG_ai2(gl,i1,j1)].y;
  x2=xygrid[GRIDG_ai2(gl,i2,j2)].x;
  y2=xygrid[GRIDG_ai2(gl,i2,j2)].y;

  direc=+1;
  if (i1!=i2) {
    if (i2<i1) direc=-1;
    N=labs(i2-i1)+1;
    L=seg[N-1];
    assert(L!=0.0e0);
    for (i=0; i<N; i++){
      xygrid[GRIDG_ai2(gl,i1+direc*i,j1)].x=(double)seg[i]/L*(x2-x1)+x1;
      xygrid[GRIDG_ai2(gl,i1+direc*i,j1)].y=(double)seg[i]/L*(y2-y1)+y1;
      xygrid[GRIDG_ai2(gl,i1+direc*i,j1)].INIT=TRUE;
    }
  }
  if (j1!=j2) {
    if (j2<j1) direc=-1;
    N=labs(j2-j1)+1;
    L=seg[N-1];
    assert(L!=0.0e0);
    for (j=0; j<N; j++){
      xygrid[GRIDG_ai2(gl,i1,j1+direc*j)].x=seg[j]/L*(x2-x1)+x1;
      xygrid[GRIDG_ai2(gl,i1,j1+direc*j)].y=seg[j]/L*(y2-y1)+y1;
      xygrid[GRIDG_ai2(gl,i1,j1+direc*j)].INIT=TRUE;
    }
  }
}

void find_xyz_along_segment(GRIDG_xyzgrid_t *xyzgrid, GRIDG_gl3d_t gl, double *seg,
                            long i1, long j1, long k1, long i2, long j2, long k2){
  double L;
  long i,j,k,N,direc;
  double x1,x2,y1,y2,z1,z2;

  x1=xyzgrid[GRIDG_ai3(gl,i1,j1,k1)].x;
  y1=xyzgrid[GRIDG_ai3(gl,i1,j1,k1)].y;
  z1=xyzgrid[GRIDG_ai3(gl,i1,j1,k1)].z;
  x2=xyzgrid[GRIDG_ai3(gl,i2,j2,k2)].x;
  y2=xyzgrid[GRIDG_ai3(gl,i2,j2,k2)].y;
  z2=xyzgrid[GRIDG_ai3(gl,i2,j2,k2)].z;

  direc=+1;
  if (i1!=i2) {
    if (i2<i1) direc=-1;
    N=labs(i2-i1)+1;
    L=seg[N-1];
    assert(L!=0.0e0);
    for (i=0; i<N; i++){
      xyzgrid[GRIDG_ai3(gl,i1+direc*i,j1,k1)].x=(double)seg[i]/L*(x2-x1)+x1;
      xyzgrid[GRIDG_ai3(gl,i1+direc*i,j1,k1)].y=(double)seg[i]/L*(y2-y1)+y1;
      xyzgrid[GRIDG_ai3(gl,i1+direc*i,j1,k1)].z=(double)seg[i]/L*(z2-z1)+z1;
      xyzgrid[GRIDG_ai3(gl,i1+direc*i,j1,k1)].INIT=TRUE;
    }
  }
  if (j1!=j2) {
    if (j2<j1) direc=-1;
    N=labs(j2-j1)+1;
    L=seg[N-1];
    assert(L!=0.0e0);
    for (j=0; j<N; j++){
      xyzgrid[GRIDG_ai3(gl,i1,j1+direc*j,k1)].x=(double)seg[j]/L*(x2-x1)+x1;
      xyzgrid[GRIDG_ai3(gl,i1,j1+direc*j,k1)].y=(double)seg[j]/L*(y2-y1)+y1;
      xyzgrid[GRIDG_ai3(gl,i1,j1+direc*j,k1)].z=(double)seg[j]/L*(z2-z1)+z1;
      xyzgrid[GRIDG_ai3(gl,i1,j1+direc*j,k1)].INIT=TRUE;
    }
  }
  if (k1!=k2) {
    if (k2<k1) direc=-1;
    N=labs(k2-k1)+1;
    L=seg[N-1];
    assert(L!=0.0e0);
    for (k=0; k<N; k++){
      xyzgrid[GRIDG_ai3(gl,i1,j1,k1+direc*k)].x=(double)seg[k]/L*(x2-x1)+x1;
      xyzgrid[GRIDG_ai3(gl,i1,j1,k1+direc*k)].y=(double)seg[k]/L*(y2-y1)+y1;
      xyzgrid[GRIDG_ai3(gl,i1,j1,k1+direc*k)].z=(double)seg[k]/L*(z2-z1)+z1;
      xyzgrid[GRIDG_ai3(gl,i1,j1,k1+direc*k)].INIT=TRUE;
    }
  }
}

int find_x_along_i(GRIDG_xgrid_t *xgrid, double segval1, double segval2,
                 double NLoverN,
                 GRIDG_gl1d_t gl, long i1, long i2, SOAP_codex_t *codex,
                 int(*find_segment)(double *, long , double,
                                    double , double , double, SOAP_codex_t * )){
  double L;
  long N;
  double *seg;
  int err;

  seg= (double *) malloc ((gl.ie-gl.is+3)*sizeof(double));
  N=labs(i2-i1)+1;
  L=fabs(xgrid[GRIDG_ai1(gl,i2)].x-xgrid[GRIDG_ai1(gl,i1)].x);
  err=(*find_segment)(seg, N, NLoverN, L, segval1, segval2, codex);
  find_x_along_segment(xgrid, gl, seg, i1, i2);
  free(seg);
  return(err);
}



int find_xy_along_ij(GRIDG_xygrid_t *xygrid, double segval1, double segval2,
                   double NLoverN,
                   GRIDG_gl2d_t gl, long i1, long j1, long i2, long j2, SOAP_codex_t *codex,
                   int(*find_segment)(double *, long , double,
                                      double , double , double, SOAP_codex_t * )){
  double L;
  long N;
  double *seg;
  int err;

  seg= (double *) malloc (max(gl.je-gl.js+3,gl.ie-gl.is+3)*sizeof(double));
  N=max(labs(j2-j1)+1,labs(i2-i1)+1);
  L=sqrt(sqr(xygrid[GRIDG_ai2(gl,i2,j2)].x-xygrid[GRIDG_ai2(gl,i1,j1)].x) +
         sqr(xygrid[GRIDG_ai2(gl,i2,j2)].y-xygrid[GRIDG_ai2(gl,i1,j1)].y));
  err=(*find_segment)(seg, N, NLoverN, L, segval1, segval2, codex);
  find_xy_along_segment(xygrid, gl, seg, i1, j1, i2, j2);
  free(seg);
  return(err);
}


int find_xyz_along_ijk(GRIDG_xyzgrid_t *xyzgrid, double segval1, double segval2,
                        double NLoverN,
                        GRIDG_gl3d_t gl, long i1, long j1, long k1, long i2, long j2, long k2, SOAP_codex_t *codex,
                        int(*find_segment)(double *, long , double,
                                           double , double , double, SOAP_codex_t * )){
  double L;
  long N;
  double *seg;
  int err;

  seg= (double *) malloc (max(gl.ke-gl.ks+3,max(gl.je-gl.js+3,gl.ie-gl.is+3))*sizeof(double));
  N=max(labs(k2-k1)+1,max(labs(j2-j1)+1,labs(i2-i1)+1));
  L=sqrt(sqr(xyzgrid[GRIDG_ai3(gl,i2,j2,k2)].x-xyzgrid[GRIDG_ai3(gl,i1,j1,k1)].x) +
         sqr(xyzgrid[GRIDG_ai3(gl,i2,j2,k2)].y-xyzgrid[GRIDG_ai3(gl,i1,j1,k1)].y) +
         sqr(xyzgrid[GRIDG_ai3(gl,i2,j2,k2)].z-xyzgrid[GRIDG_ai3(gl,i1,j1,k1)].z));
  err=(*find_segment)(seg, N, NLoverN, L, segval1, segval2, codex);
  find_xyz_along_segment(xyzgrid, gl, seg, i1, j1, k1, i2, j2, k2);
  free(seg);
  return(err);
}


void erode_bdry(GRIDG_xygrid_t *xygrid, long erofact,
                GRIDG_gl2d_t gl, long i1, long j1, long i2, long j2){
  long i,j,ero;
  GRIDG_xygrid_t *xygrid2;
  xygrid2 = (GRIDG_xygrid_t *) malloc((gl.ie-gl.is+1+4)*(gl.je-gl.js+1+4)
             *sizeof(GRIDG_xygrid_t));
  for (ero=1; ero<=erofact; ero++){
    if (i1==i2) {
      xygrid2[GRIDG_ai2(gl,i1,j1)].x=xygrid[GRIDG_ai2(gl,i1,j1)].x;
      xygrid2[GRIDG_ai2(gl,i1,j1)].y=xygrid[GRIDG_ai2(gl,i1,j1)].y;
      xygrid2[GRIDG_ai2(gl,i1,j2)].x=xygrid[GRIDG_ai2(gl,i1,j2)].x;
      xygrid2[GRIDG_ai2(gl,i1,j2)].y=xygrid[GRIDG_ai2(gl,i1,j2)].y;
      for (j=j1+1; j<j2; j++){
        xygrid2[GRIDG_ai2(gl,i1,j)].x=0.4e0*xygrid[GRIDG_ai2(gl,i1,j-1)].x
                               +0.2e0*xygrid[GRIDG_ai2(gl,i1,j+0)].x
                               +0.4e0*xygrid[GRIDG_ai2(gl,i1,j+1)].x;
        xygrid2[GRIDG_ai2(gl,i1,j)].y=0.4e0*xygrid[GRIDG_ai2(gl,i1,j-1)].y
                               +0.2e0*xygrid[GRIDG_ai2(gl,i1,j+0)].y
                               +0.4e0*xygrid[GRIDG_ai2(gl,i1,j+1)].y;
      }
      for (j=j1+1; j<j2; j++){
        xygrid[GRIDG_ai2(gl,i1,j)].x=0.4e0*xygrid2[GRIDG_ai2(gl,i1,j-1)].x
                              +0.2e0*xygrid2[GRIDG_ai2(gl,i1,j+0)].x
                              +0.4e0*xygrid2[GRIDG_ai2(gl,i1,j+1)].x;
        xygrid[GRIDG_ai2(gl,i1,j)].y=0.4e0*xygrid2[GRIDG_ai2(gl,i1,j-1)].y
                              +0.2e0*xygrid2[GRIDG_ai2(gl,i1,j+0)].y
                              +0.4e0*xygrid2[GRIDG_ai2(gl,i1,j+1)].y;
        xygrid[GRIDG_ai2(gl,i1,j)].INIT=TRUE;
      }

    } else {
      xygrid2[GRIDG_ai2(gl,i1,j1)].x=xygrid[GRIDG_ai2(gl,i1,j1)].x;
      xygrid2[GRIDG_ai2(gl,i1,j1)].y=xygrid[GRIDG_ai2(gl,i1,j1)].y;
      xygrid2[GRIDG_ai2(gl,i2,j1)].x=xygrid[GRIDG_ai2(gl,i2,j1)].x;
      xygrid2[GRIDG_ai2(gl,i2,j1)].y=xygrid[GRIDG_ai2(gl,i2,j1)].y;
      for (i=i1+1; i<i2; i++){
        xygrid2[GRIDG_ai2(gl,i,j1)].x=0.4e0*xygrid[GRIDG_ai2(gl,i-1,j1)].x
                               +0.2e0*xygrid[GRIDG_ai2(gl,i+0,j1)].x
                               +0.4e0*xygrid[GRIDG_ai2(gl,i+1,j1)].x;
        xygrid2[GRIDG_ai2(gl,i,j1)].y=0.4e0*xygrid[GRIDG_ai2(gl,i-1,j1)].y
                               +0.2e0*xygrid[GRIDG_ai2(gl,i+0,j1)].y
                               +0.4e0*xygrid[GRIDG_ai2(gl,i+1,j1)].y;
      }
      for (i=i1+1; i<i2; i++){
        xygrid[GRIDG_ai2(gl,i,j1)].x=0.4e0*xygrid2[GRIDG_ai2(gl,i-1,j1)].x
                              +0.2e0*xygrid2[GRIDG_ai2(gl,i+0,j1)].x
                              +0.4e0*xygrid2[GRIDG_ai2(gl,i+1,j1)].x;
        xygrid[GRIDG_ai2(gl,i,j1)].y=0.4e0*xygrid2[GRIDG_ai2(gl,i-1,j1)].y
                              +0.2e0*xygrid2[GRIDG_ai2(gl,i+0,j1)].y
                              +0.4e0*xygrid2[GRIDG_ai2(gl,i+1,j1)].y;
        xygrid[GRIDG_ai2(gl,i,j1)].INIT=TRUE;
      }
    }
  }
  free(xygrid2);
}



double read_double_after_X_chars(FILE *infile, long X, char *varname,
                           bool VERBOSE){
  double tmp;
  read_x_chars_from_file(infile,X);
  if (fscanf(infile, "%lg%*[^\n]", &tmp)!=1) GRIDG_fatal_error("The function fscanf could not read properly in read_double_after_X_chars().");
  if (VERBOSE) printf("%s=%E\n",varname,tmp);
  return(tmp);
}


long read_long_after_X_chars(FILE *infile, long X, char *varname,
                       bool VERBOSE){
  long tmp;
  read_x_chars_from_file(infile,X);
  if (fscanf(infile, "%ld%*[^\n]", &tmp)!=1) GRIDG_fatal_error("The function fscanf could not read properly in read_long_after_X_chars().");
  if (VERBOSE) printf("%s=%ld\n",varname,tmp);
  return(tmp);
}


void grab_code(FILE *infile, char **code, long *codelength){
  long cnt,bracketcnt;
  char tmp;
  *code=(char *)malloc(sizeof(char));
  do {
    tmp=fgetc(infile);
  } while (tmp!='(');
  bracketcnt=1;
  cnt=0;
  do {
    *code=(char *)realloc(*code,(cnt+1)*sizeof(char));
    (*code)[cnt]=fgetc(infile);
    if ((*code)[cnt]=='(') bracketcnt++;
    if ((*code)[cnt]==')') bracketcnt--;
    cnt++;
  } while (bracketcnt!=0);
  *codelength=cnt-1;
}


long sgndif(long i1, long i2){
  long tmp;
  tmp=0;
  if (i2>i1) tmp=-1;
  if (i2<i1) tmp=+1;
  return(tmp);
}

/* returns 0 on success, 1 on failure */


int find_segment_type(int (**FindSegment)(double *, long , double, double , double , double, SOAP_codex_t * ),
                     char modeL, char modeH, SOAP_codex_t *codex){
  bool FOUND;
  FOUND=FALSE;
  if (modeL=='E' && modeH=='E') {
    *FindSegment=&find_segment_1;
    FOUND=TRUE;
  }
  if ((modeL=='F' || modeL=='G') && modeH=='E') {
    *FindSegment=&find_segment_4;
    FOUND=TRUE;
  }
  if (modeL=='E' && (modeH=='F' || modeH=='G')) {
    *FindSegment=&find_segment_5;
    FOUND=TRUE;
  }
  if ((modeL=='F' || modeL=='G') && (modeH=='F' || modeH=='G')) {
    *FindSegment=&find_segment_3;
    FOUND=TRUE;
  }
  if ((modeL=='f' || modeL=='g') && (modeH=='f' || modeH=='g')) {
    *FindSegment=&find_segment_6;
    FOUND=TRUE;
  }
  if (modeL=='N' && modeH=='O') {
    *FindSegment=&find_segment_NO;
    FOUND=TRUE;
  }
  if (!FOUND) {
    SOAP_fatal_error(codex,"Could not find a segment corresponding to the mode >%c%c<.",modeL,modeH);
  } 
  return(0);
}



static void RXaction2(char *actionname, char **argum, SOAP_codex_t *codex){
  bool *FOUND;
  char *actionname2;
  SOAP_codex_t codex2;
  FOUND=&(((RXaction_args2_t *)(codex->action_args))->FOUND);
  actionname2=((RXaction_args2_t *)(codex->action_args))->actionname;
  if (strcmp(actionname,actionname2)==0) {
    codex2.ACTION=FALSE;
    codex2.FUNCTION=FALSE;
    codex2.VERBOSE=codex->VERBOSE;
    codex2.vars=codex->vars;
    codex2.filename=codex->filename;
    *FOUND=TRUE;
    ((RXaction_args2_t *)(codex->action_args))->action
        (argum,((RXaction_args2_t *)(codex->action_args))->action_args,&codex2);
    codex->vars=codex2.vars;
  }
}



void read_RX_grid(char *filename, char *actionname,
                void (*RXaction)(char **, void *, SOAP_codex_t *),
                void (*WriteRXaction)(FILE **),
                void *RXaction_args, bool VERBOSE){

  char *code;
  FILE *file_a;
  SOAP_codex_t codex;
  RXaction_args2_t RXaction_args2;

  SOAP_init_codex(&codex,filename);

  codex.VERBOSE=VERBOSE;
  codex.ACTION=TRUE;
  codex.action_args=&RXaction_args2;
  codex.action=&RXaction2;

  RXaction_args2.FOUND=FALSE;
  RXaction_args2.actionname=actionname;
  RXaction_args2.action=RXaction;
  RXaction_args2.action_args=RXaction_args;

  code=(char *)malloc(sizeof(char));
  SOAP_store_file_as_string(filename, &code);
  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);

  if (!RXaction_args2.FOUND) {
    if (VERBOSE) printf("%s was not found in file %s. Appending it.\n",actionname,filename);
    fprintf(stdout,"Module missing from control file %s. "
                   "Appending defaults.\n",filename);
    file_a = fopen(filename, "a");
    if (file_a == NULL) GRIDG_fatal_error("Could not open file %s for appending.", filename);
    fprintf(file_a,"\n{ the following has been automatically appended }");
    (*WriteRXaction)(&file_a);
    fclose(file_a);
    SOAP_store_file_as_string(filename, &code);
    SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
    if (!RXaction_args2.FOUND) {
      GRIDG_fatal_error("Even after appending, the module could not be found.");
    }
  }

  free(code);
  SOAP_free_codex(&codex);
}


void replace_equal_sign(char **expr){
  long cnt;
  *expr=(char *)realloc(*expr,((long)strlen(*expr)+4)*sizeof(char));
  cnt=0;
  while (((*expr)[cnt]!='=') && (*expr)[cnt]!=EOS) cnt++;
  if ((*expr)[cnt]=='=') {
    SOAP_strcut(cnt,cnt,*expr);
    SOAP_strins("-(",expr,cnt);
    SOAP_strins(")",expr,(long)strlen(*expr));
  }
}


void GRIDG_fatal_error(const char *formatstr, ...){
  va_list ap;
  char *newstr;
  int term_width,term_height;
  newstr=(char *)malloc(10000*sizeof(char));
  fprintf(stderr,"\n\n");
  va_start(ap, formatstr);
  vsprintf(newstr,formatstr, ap);
   va_end(ap);
  find_terminal_window_size(&term_width,&term_height);
  fprintf(stderr,"%s",strwrp(newstr,min(term_width-1,70)));
  free(newstr);

  fprintf(stderr,"\n\nGRIDG fatal error. Exiting.\n\n");
  exit(EXIT_FAILURE);
}


