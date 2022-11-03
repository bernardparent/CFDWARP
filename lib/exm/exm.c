// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2002 Bernard Parent
Copyright 2019 Jaehyuk Lee

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


/* 

EXM: External Module 
Functions that can be used with CFDWARP or any other code 

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <stdarg.h>
#include "exm.h"



void EXM_fatal_error(const char *formatstr, ...){
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

  fprintf(stderr,"\n\nEXM fatal error. Exiting.\n\n");
  exit(EXIT_FAILURE);
}


long EXM_ai3(EXM_gl3D_t gl, long i, long j, long k) {
  long ii;
  long isr,jsr,jer,ksr,ker;
  isr=gl.is;
  jsr=gl.js;
  jer=gl.je;
  ksr=gl.ks;
  ker=gl.ke;
  ii=(i-isr)*(jer-jsr+1)+(j-jsr);
  ii=ii*(ker-ksr+1)+(k-ksr);
  return(ii);
}

long EXM_ai2(EXM_gl2D_t gl, long i, long j) {
  long ii;
  long isr,jsr,jer;
  isr=gl.is;
  jsr=gl.js;
  jer=gl.je;
  ii=(i-isr)*(jer-jsr+1)+(j-jsr);
  return(ii);
}

long EXM_ai1(EXM_gl1D_t gl, long i) {
  long ii;
  long isr;
  isr=gl.is;
  ii=(i-isr);
  return(ii);
}

long EXM_aim(EXM_glm_t gl, long row, long col) {
  long ii;
  ii=row*gl.numcol+col;
  return(ii);
}

/*===========================================================================
find root solver
============================================================================*/

double SIGN(NUM1, NUM2)
double NUM1, NUM2;
{
  double NUM3;
  double tmp;

  if (NUM2<0.0e0) tmp=-1.0e0; else tmp=1.0e0;
  NUM3 = NUM1 * tmp;
  return NUM3;
}




/*    BASED ON A METHOD BY T J DEKKER
      WRITTEN BY L F SHAMPINE AND H A WATTS
      MODIFIED FOR THE MATH LIBRARY BY C B BAILEY
      TRANSLATED FROM FORTRAN TO PASCAL TO C BY B PARENT

      ABSTRACT
         EXM_find_root_zero_in SEARCHES FOR A ZERO OF A FUNCTION F(X) BETWEEN
         THE GIVEN VALUES B AND C UNTIL THE WIDTH OF THE INTERVAL
         (B,C) HAS COLLAPSED TO WITHIN A TOLERANCE SPECIFIED BY
         THE STOPPING CRITERION, ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
         THE METHOD USED IS AN EFFICIENT COMBINATION OF BISECTION AND
         THE SECANT RULE.  IN ORDER TO INSURE THAT EXM_find_root_zero_in WILL CONVERGE
         TO A ZERO, THE USER SHOULD PICK VALUES FOR B AND C AT WHICH
         THE FUNCTION DIFFERS IN SIGN.

      DESCRIPTION OF ARGUMENTS
      F,B,C,RE AND AE ARE INPUT PARAMETERS
      B,C AND IFLAG ARE OUTPUT PARAMETERS
         F     - NAME OF THE REAL VALUED EXTERNAL FUNCTION.  THIS NAME
                 MUST BE IN AN EXTERNAL STATEMENT IN THE CALLING
                 PROGRAM.  F MUST BE A FUNCTION OF ONE REAL ARGUMENT.
         minval- ONE END OF THE INTERVAL (B,C).  THE VALUE RETURNED FOR
                 B USUALLY IS THE BETTER APPROXIMATION TO A ZERO OF F.
         maxval- THE OTHER END OF THE INTERVAL (B,C)
         relerr- RELATIVE ERROR USED FOR RW IN THE STOPPING CRITERION.
                 IF THE REQUESTED RE IS LESS THAN MACHINE PRECISION,
                 THEN RW IS SET TO APPROXIMATELY MACHINE PRECISION.
         abserr- ABSOLUTE ERROR USED IN THE STOPPING CRITERION.  IF THE
                 GIVEN INTERVAL (B,C) CONTAINS THE ORIGIN, THEN A
                 NONZERO VALUE SHOULD BE CHOSEN FOR AE.
         IFLAG - A STATUS CODE.  USER MUST CHECK IFLAG AFTER EACH CALL.
                 CONTROL RETURNS TO THE USER FROM EXM_find_root_zero_in IN ALL CASES.
                 XERROR DOES NOT PROCESS DIAGNOSTICS IN THESE CASES.
                  1 B IS WITHIN THE REQUESTED TOLERANCE OF A ZERO.
                    THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
                    TOLERANCE, THE FUNCTION CHANGES SIGN IN (B,C), AND
                    F(X) DECREASED IN MAGNITUDE AS (B,C) COLLAPSED.
                  2 F(B) = 0.  HOWEVER, THE INTERVAL (B,C) MAY NOT HAVE
                    COLLAPSED TO THE REQUESTED TOLERANCE.
                  3 B MAY BE NEAR A SINGULAR POINT OF F(X).
                    THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
                    TOLERANCE AND THE FUNCTION CHANGES SIGN IN (B,C) BUT
                    F(X) INCREASED IN MAGNITUDE AS (B,C) COLLAPSED,I.E.
                      ABS(F(B OUT)) .GT. MAX(ABS(F(B IN)),ABS(F(C IN)))
                  4 NO CHANGE IN SIGN OF F(X) WAS FOUND ALTHOUGH THE
                    INTERVAL (B,C) COLLAPSED TO THE REQUESTED TOLERANCE.
                    THE USER MUST EXAMINE THIS CASE AND DECIDE WHETHER
                    B IS NEAR A LOCAL MINIMUM OF F(X), OR B IS NEAR A
                    ZERO OF EVEN MULTIPLICITY, OR NEITHER OF THESE.
                  5 TOO MANY (.GT. 500) FUNCTION EVALUATIONS USED.

      REFERENCES
        1.  L F SHAMPINE AND H A WATTS, EXM_find_root_zero_in, A ROOT-SOLVING CODE,
            SC-TM-70-631, SEPT 1970.
        2.  T J DEKKER, FINDING A ZERO BY MEANS OF SUCCESSIVE LINEAR
            INTERPOLATION, *CONSTRUCTIVE ASPECTS OF THE FUNDAMENTAL
            THEOREM OF ALGEBRA*, EDITED BY B DEJON AND P HENRICI, 1969.


      ER IS TWO TIMES THE COMPUTER UNIT ROUNDOFF VALUE WHICH IS
      DEFINED HERE TO BE THE VALUE FOR THE IBM PC DOUBLE PRECISION*/




double EXM_find_root_zero_in(double(*FUNCT)(void *, double), void *arg_ptr,
            double minval, double maxval,
            double relerr, double abserr, long *IFLAG){
  double T, A, P, U, TOL, ACMB, CMB, FX, FC, FB, FA, ACBS, ER, RW, AW;
  long IC, KOUNT;

  *IFLAG = 0;
  ER = 2.0e-13;
  RW = max(relerr, ER);
  AW = max(abserr, 0.0e0);
  IC = 0;
  ACBS = fabs(minval - maxval);
  A = maxval;
  T = A;
  FA = (*FUNCT)(arg_ptr,T);
  T = minval;
  FB = (*FUNCT)(arg_ptr,T);
  FC = FA;
  KOUNT = 2;
  FX = max(fabs(FB), fabs(FC));
_L1:
  if (fabs(FC) > fabs(FB))
    goto _L2;
  A = minval;
  FA = FB;
  minval = maxval;
  FB = FC;
  maxval = A;
  FC = FA;
_L2:
  if (FB == 0.0e0)
    *IFLAG = 2;
  CMB = 0.5e0 * (maxval - minval);
  ACMB = fabs(CMB);
  TOL = RW * fabs(minval) + AW;
  if (ACMB <= TOL) {
    *IFLAG = 1;
    if (SIGN(1.0e0, FB) == SIGN(1.0e0, FC))
      *IFLAG = 4;
    if (fabs(FB) > FX)
      *IFLAG = 3;
  }
  P = (minval - A) * FB;
  U = FA - FB;
  if (P >= 0.0e0)
    goto _L3;
  P = -P;
  U = -U;
_L3:
  A = minval;
  FA = FB;
  IC++;
  if (IC < 4)
    goto _L4;
  if (8.0e0 * ACMB >= ACBS)
    goto _L6;
  IC = 0;
  ACBS = ACMB;
_L4:
  if (P > fabs(U) * TOL)
    goto _L5;
  minval += SIGN(TOL, CMB);
  goto _L7;
_L5:
  if (P >= CMB * U)
    goto _L6;
  assert(U!=0.0e0);
  minval += P / U;
  goto _L7;
_L6:
  minval = 0.5e0 * (maxval + minval);
_L7:
  T = minval;
  FB = (*FUNCT)(arg_ptr,T);
  if (FB == 0.0e0)
    *IFLAG = 2;
  if (SIGN(1.0e0, FB) != SIGN(1.0e0, FC))
    goto _L8;
  maxval = A;
  FC = FA;
_L8:
  KOUNT++;
  if (KOUNT > 500)
    *IFLAG = 5;
  if (*IFLAG == 0)
    goto _L1;

  return(minval);
}


/*
   root_guess: enter a value as a first guess of the root
   droot_init: a small value, usually set to about 10^8 times smaller than
               the maximum value of a root you would expect to
               obtain
   relerr: the maximum admissible relative error on the residual to
           obtain convergence
   abserr: the maximum admissible absolute error on the residual
           to obtain convergence
   *IFLAG= 1: convergence has been obtained correctly;
           2: problem in convergence, too many iterations
           3: droot_init can not be zero
           4: problem finding adequate dx from droot_init (droot_init may be too small or too large) 
           5: the function provided returned NaN or a non finite number
*/

double EXM_find_root_Newton_Raphson(
              double(*FUNCT)(void *, double), void *arg_ptr,
              double root_guess, double droot_init,
              double relerr, double abserr, long *IFLAG){
  double x1,x2,dx;
  double res1,res2,resref;
  long cnt;
  bool PROBLEM,NEEDEDADJUSTMENT;
  x1=0.0;
  *IFLAG=0;
  if (droot_init==0.0) *IFLAG=3;
  /* check if droot_init is proper and adjust it if required */
  PROBLEM=FALSE;
  NEEDEDADJUSTMENT=FALSE;
  res1=(*FUNCT)(arg_ptr,x1);
  cnt=0;
  do {
    res2=(*FUNCT)(arg_ptr,x1+droot_init);
    if (res2-res1==0.0) {
      PROBLEM=TRUE;
      NEEDEDADJUSTMENT=TRUE;
      droot_init*=10.0;
    } else {
      PROBLEM=FALSE;
    }
    cnt++;
  } while (PROBLEM && cnt<=100);
  if (NEEDEDADJUSTMENT) droot_init*=100.0;
  if (cnt>=100) *IFLAG=4;
  resref=max(fabs(res2),max(fabs(res1),fabs((*FUNCT)(arg_ptr,x1+droot_init))));
  if (isnan(resref) || !isfinite(resref)) *IFLAG=5;
  if (*IFLAG==0){
    PROBLEM=FALSE;
    x1=root_guess;
    res1=(*FUNCT)(arg_ptr,x1);
    dx=droot_init;
    cnt=0;
    do {
      cnt++;
      x2=x1+dx;
      res2=(*FUNCT)(arg_ptr,x2);
      if (res2-res1==0.0 || x2-x1==0.0){
        droot_init*=2.0;
        dx=droot_init;
        PROBLEM=TRUE;
      } else {
        PROBLEM=FALSE;
        dx=-res2/(res2-res1)*(x2-x1);
      }
      res1=res2;
      x1=x2;
    } while ((fabs(res1)>=abserr || PROBLEM) && (fabs(res1/resref)>=relerr || PROBLEM)
           && (cnt<=300));
    *IFLAG=1;
    if (cnt>=300) *IFLAG=2;
    if (isnan(res1) || !isfinite(res1)) *IFLAG=5;
  }
  return(x1);
}






/*============================================================================
XDMA solver
=============================================================================*/

/*
#define maxXDMAthread 5


typedef struct {
  EXM_gl2D_t *gl;
  long hbw,line,line2min,line2max;
  double *xdma;
} XDMAthreadptr_t;


void *XDMAthreadfunct(void *XDMAthreadptr){
  EXM_gl2D_t *gl;
  long col,hbw,line,line2min,line2max,line2,dc;
  double fact;
  double *xdma;
  
  xdma=((XDMAthreadptr_t *)XDMAthreadptr)->xdma;
  gl=((XDMAthreadptr_t *)XDMAthreadptr)->gl;
  hbw=((XDMAthreadptr_t *)XDMAthreadptr)->hbw;
  line=((XDMAthreadptr_t *)XDMAthreadptr)->line;
  line2min=((XDMAthreadptr_t *)XDMAthreadptr)->line2min;
  line2max=((XDMAthreadptr_t *)XDMAthreadptr)->line2max;
  for (line2=line2min; line2<=line2max; line2++){
    dc=(line2-line);
    // here the idea is to add line*fact to line2, with dc added to the
    // column index of line 
    fact=-xdma[EXM_ai2(*gl,hbw-dc,line2)]/xdma[EXM_ai2(*gl,hbw,line)];
    for (col=hbw; col<gl->ie; col++){
      xdma[EXM_ai2(*gl,col-dc,line2)]+=fact*xdma[EXM_ai2(*gl,col,line)];
    }
    xdma[EXM_ai2(*gl,gl->ie,line2)]+=fact*xdma[EXM_ai2(*gl,gl->ie,line)];
  }
  return(NULL);
}


// solves a x-diagonal matrix as defined by the pointer variable xdma with the size
// in gl
void EXM_SolveXDMAthread(double *xdma, EXM_gl2D_t gl){
  long numXDMAthread,cntXDMAthread,hbw,line,line2,dc,linemax,col,line2min,line2max,numline;
  double fact;
  void *retval;
  pthread_t XDMAthread[maxXDMAthread];
  XDMAthreadptr_t XDMAthreadptr[maxXDMAthread];

// NOTE: gl.js must be equal to zero
//        gl.is must be equal to zero 

  hbw=(gl.ie-gl.is+1)/2-1;
  linemax=gl.je;

  for (line=0; line<linemax; line++){
    line2min=line+1;
    line2max=min(line+hbw,linemax);
    cntXDMAthread=0;
    numline=max(1,round((line2max-line2min)/maxXDMAthread));
    do {
      XDMAthreadptr[cntXDMAthread].gl=&gl;
      XDMAthreadptr[cntXDMAthread].hbw=hbw;
      XDMAthreadptr[cntXDMAthread].line=line;
      XDMAthreadptr[cntXDMAthread].xdma=xdma;    
      XDMAthreadptr[cntXDMAthread].line2min=line2min+cntXDMAthread*numline;
      XDMAthreadptr[cntXDMAthread].line2max=min(line2max,line2min+(cntXDMAthread+1)*numline-1);
      if (cntXDMAthread==maxXDMAthread-1) XDMAthreadptr[cntXDMAthread].line2max=line2max;
      if (line==0 && FALSE) fprintf(stderr,"%ld  %ld   %ld %ld  %ld  %ld\n",line2min,line2max,XDMAthreadptr[cntXDMAthread].line2min,XDMAthreadptr[cntXDMAthread].line2max,hbw,numline);
      if (pthread_create(&((XDMAthread)[cntXDMAthread]), NULL, &XDMAthreadfunct,
                         (void *)(&(XDMAthreadptr[cntXDMAthread]))))
        fprintf(stderr,"Cannot create XDMA thread.\n");
      cntXDMAthread++;
    } while (cntXDMAthread<maxXDMAthread || XDMAthreadptr[cntXDMAthread-1].line2max!=line2max);
    numXDMAthread=cntXDMAthread;
    for (cntXDMAthread=0; cntXDMAthread<numXDMAthread; cntXDMAthread++){
      if (pthread_join(XDMAthread[cntXDMAthread],&retval))
        fprintf(stderr,"Cannot join XDMA thread %ld.\n",cntXDMAthread);
    }
  }


  for (line=linemax; line>0; line--){
    for (line2=line-1; line2>=max(line-hbw,0); line2--){
      dc=(line2-line);
      // here the idea is to add line*fact to line2, with dc added to the
      //   column index of line 
      fact=-xdma[EXM_ai2(gl,hbw-dc,line2)]/xdma[EXM_ai2(gl,hbw,line)];
      col=hbw;
      xdma[EXM_ai2(gl,col-dc,line2)]+=fact*xdma[EXM_ai2(gl,col,line)];

      xdma[EXM_ai2(gl,gl.ie,line2)]+=fact*xdma[EXM_ai2(gl,gl.ie,line)];
    }
  }

  for (line=0; line<=linemax; line++){
     xdma[EXM_ai2(gl,gl.ie,line)]/=xdma[EXM_ai2(gl,hbw,line)];
     xdma[EXM_ai2(gl,hbw,line)]=1.0e0;
  }

}

*/



/*solves a x-diagonal matrix as defined by the pointer variable xdma with the size
 in gl*/
void EXM_solve_XDMA(double *xdma, EXM_gl2D_t gl){
  long hbw,line,line2,dc,linemax,col;
  double fact,sum,aux;

/* NOTE: gl.js must be equal to zero
         gl.is must be equal to zero */

  hbw=(gl.ie-gl.is+1)/2-1;
  linemax=gl.je;

  for (line=0; line<linemax; line++){

#if defined(OPENMPTHREADS) 
#pragma omp parallel for private(line2,dc,fact,col) schedule(static) 
#endif
    for (line2=line+1; line2<=min(line+hbw,linemax); line2++){
      dc=(line2-line);
      // here the idea is to add line*fact to line2, with dc added to the
      //   column index of line 
      assert(xdma[EXM_ai2(gl,hbw,line)]!=0.0);
      fact=-xdma[EXM_ai2(gl,hbw-dc,line2)]/xdma[EXM_ai2(gl,hbw,line)];
      
      for (col=hbw; col<gl.ie; col++){
        xdma[EXM_ai2(gl,col-dc,line2)]+=fact*xdma[EXM_ai2(gl,col,line)];
      }
      xdma[EXM_ai2(gl,gl.ie,line2)]+=fact*xdma[EXM_ai2(gl,gl.ie,line)];
    }
  }


  xdma[EXM_ai2(gl,gl.ie,linemax)]/=xdma[EXM_ai2(gl,hbw,linemax)];
  xdma[EXM_ai2(gl,hbw,linemax)]=1.0;
  for (line=linemax-1; line>=0; line--){
    sum=0.0;
#if defined(OPENMPTHREADS) 
#pragma omp parallel for reduction(+:sum) private(line2,aux) schedule(static)  
#endif
    for (line2=line+1; line2<=min(line+hbw,linemax); line2++){
      aux=xdma[EXM_ai2(gl,hbw+line2-line,line)]*xdma[EXM_ai2(gl,gl.ie,line2)];
      sum=sum+aux;
    }
    xdma[EXM_ai2(gl,gl.ie,line)]-=sum;
    assert(xdma[EXM_ai2(gl,hbw,line)]!=0.0);
    xdma[EXM_ai2(gl,gl.ie,line)]/=xdma[EXM_ai2(gl,hbw,line)];
    xdma[EXM_ai2(gl,hbw,line)]=1.0;
  }

}



void EXM_solve_XDMA_old(double *xdma, EXM_gl2D_t gl){
  long hbw,line,line2,dc,linemax,col;
  double fact;

/* NOTE: gl.js must be equal to zero
         gl.is must be equal to zero */

  hbw=(gl.ie-gl.is+1)/2-1;
  linemax=gl.je;

  for (line=0; line<linemax; line++){
    for (line2=line+1; line2<=min(line+hbw,linemax); line2++){
      dc=(line2-line);
      /* here the idea is to add line*fact to line2, with dc added to the
         column index of line */
      fact=-xdma[EXM_ai2(gl,hbw-dc,line2)]/xdma[EXM_ai2(gl,hbw,line)];
      for (col=hbw; col<gl.ie; col++){
        xdma[EXM_ai2(gl,col-dc,line2)]+=fact*xdma[EXM_ai2(gl,col,line)];
      }
      xdma[EXM_ai2(gl,gl.ie,line2)]+=fact*xdma[EXM_ai2(gl,gl.ie,line)];
    }
  }


  for (line=linemax; line>0; line--){
    for (line2=line-1; line2>=max(line-hbw,0); line2--){
      dc=(line2-line);
      /* here the idea is to add line*fact to line2, with dc added to the
         column index of line */
      fact=-xdma[EXM_ai2(gl,hbw-dc,line2)]/xdma[EXM_ai2(gl,hbw,line)];
      col=hbw;
      xdma[EXM_ai2(gl,col-dc,line2)]+=fact*xdma[EXM_ai2(gl,col,line)];

      xdma[EXM_ai2(gl,gl.ie,line2)]+=fact*xdma[EXM_ai2(gl,gl.ie,line)];
    }
  }

  for (line=0; line<=linemax; line++){
     xdma[EXM_ai2(gl,gl.ie,line)]/=xdma[EXM_ai2(gl,hbw,line)];
     xdma[EXM_ai2(gl,hbw,line)]=1.0e0;
  }

}



/*============================================================================
TDMA solver
=============================================================================*/

/*solves a tri-diagonal matrix as defined by the pointer variables tdm lines*/
void EXM_solve_TDMA(EXM_tdmaline_t *tdma, long numlines)
{
  long line,cnt;
  double tmp;

  /*tdma[0].val[0] must be equal to zero*/
  /*tdma[numlines-1].val[2] must be equal to zero*/
  for (line=0; line<numlines-1; line++){
    tmp = -(tdma[line+1].val[0] / tdma[line].val[1]);
    for (cnt = 1; cnt <= 2; cnt++)
      tdma[line+1].val[cnt - 1] += tdma[line].val[cnt] * tmp;
    tdma[line+1].val[3] += tdma[line].val[3] * tmp;
    tdma[line+1].val[0] = 0.0;
  }
  for (line=numlines-1; line>0; line--){
    tdma[line].val[3] /= tdma[line].val[1];
    tdma[line].val[1] = 1.0;
    tdma[line-1].val[3] -= tdma[line].val[3] * tdma[line-1].val[2];
    tdma[line-1].val[2] = 0.0;
  }
  tdma[0].val[3] /= tdma[0].val[1];
  tdma[0].val[1] = 1.0;
}


/*============================================================================
PDMA solver
=============================================================================*/

/*solves a penta-diagonal matrix as defined by the pointer variables pdm lines*/
void EXM_solve_PDMA(EXM_pdmaline_t *pdma, long numlines)
{
  long line;
  double tmp;

  /* first sweep downwards */
  for (line=1; line<numlines-1; line++){
    tmp = -(pdma[line+1].val[0] / pdma[line].val[1]);
    pdma[line+1].val[0] = 0.0;
    pdma[line+1].val[1] += pdma[line].val[2] * tmp;
    pdma[line+1].val[2] += pdma[line].val[3] * tmp;
    pdma[line+1].val[3] += pdma[line].val[4] * tmp;
    pdma[line+1].val[5] += pdma[line].val[5] * tmp;
  }

  /* second sweep downwards */
  for (line=0; line<numlines-1; line++){
    tmp = -(pdma[line+1].val[1] / pdma[line].val[2]);
    pdma[line+1].val[1] = 0.0;
    pdma[line+1].val[2] += pdma[line].val[3] * tmp;
    pdma[line+1].val[3] += pdma[line].val[4] * tmp;
    pdma[line+1].val[5] += pdma[line].val[5] * tmp;
  }

  /* first sweep upwards */ 
  for (line=numlines-2; line>0; line--){
    tmp = -(pdma[line-1].val[4] / pdma[line].val[3]);
    pdma[line-1].val[4] =0.0;
    pdma[line-1].val[3] += pdma[line].val[2]*tmp ;
    pdma[line-1].val[2] += pdma[line].val[1]*tmp ;
    pdma[line-1].val[1] += pdma[line].val[0]*tmp ;
    pdma[line-1].val[5] += pdma[line].val[5]*tmp;

  }


  /* second sweep upwards */
  for (line=numlines-1; line>0; line--){
    pdma[line].val[5] /= pdma[line].val[2];
    pdma[line].val[2] = 1.0;
    pdma[line-1].val[5] -= pdma[line].val[5] * pdma[line-1].val[3];
    pdma[line-1].val[3] = 0.0;
  }


  pdma[0].val[5] /= pdma[0].val[2];
  pdma[0].val[2] = 1.0;
}











/*======================================================================
            Block TDMA solver
  ====================================================================== */


long EXM_mi(long k, long line, long row, long col){
  long tmp;
  tmp=line*k*k+col*k+row;
  return(tmp);
}


static void mat_add_tdma(double *mat1, long line1, double *mat2, long line2,
                  double *mat3, long line3, long k){
  long row,col;
  for (row=0; row<k; row++){
    for (col=0; col<k; col++){
       mat3[EXM_mi(k,line3,row,col)]=
           mat1[EXM_mi(k,line1,row,col)]
          +mat2[EXM_mi(k,line2,row,col)];
    }
  }
}



static void mat_mult_tdma(double *mat1, long line1, double *mat2,long line2,
                 double *mat3, long line3, long k){
  long row,col,cnt;
  for (row=0; row<k; row++){
    for (col=0; col<k; col++){
       mat3[EXM_mi(k,line3,row,col)]=0.0e0;
       for (cnt=0; cnt<k; cnt++){
           mat3[EXM_mi(k,line3,row,col)]=mat3[EXM_mi(k,line3,row,col)]+
             mat1[EXM_mi(k,line1,row,cnt)]*mat2[EXM_mi(k,line2,cnt,col)];
       }
    }
  }
}
	




static void mat_init_i_tdma(double *mat1, long line1, long k){
  long m,n;
  for (m=0; m<k; m++){
    for (n=0; n<k; n++){
      mat1[EXM_mi(k,line1,m,n)]=0.0e0;
    }
  }
  for (m=0; m<k; m++){
    mat1[EXM_mi(k,line1,m,m)]=1.0e0;
  }
} 



static void mat_inv_tdma(double *mat1,long line1,double *mat2,long line2,long k){
  long row,row2,col;
  double multfact;
/*
     Idea: init mat2 as identity; gaussian elimination on mat1/mat2 
*/
  mat_init_i_tdma(mat2,line2,k);
  for (row=0; row<k; row++){
    for (row2=0; row2<k; row2++){
      if (row2 != row) {
        multfact=-mat1[EXM_mi(k,line1,row2,row)]
                 /(mat1[EXM_mi(k,line1,row,row)]+0.0e-30);
/*
           Add line row*multfact to line row2
*/ 
        for (col=0; col<k; col++){
          mat1[EXM_mi(k,line1,row2,col)]=mat1[EXM_mi(k,line1,row2,col)]+
                      mat1[EXM_mi(k,line1,row,col)]*multfact;
          mat2[EXM_mi(k,line2,row2,col)]=mat2[EXM_mi(k,line2,row2,col)]+
                     mat2[EXM_mi(k,line2,row,col)]*multfact;
        }
      }
    }
  }

  for (row=0; row<k; row++){
    for (col=0; col<k; col++){
      assert(mat1[EXM_mi(k,line1,row,row)]!=0.0e0);
      mat2[EXM_mi(k,line2,row,col)]=mat2[EXM_mi(k,line2,row,col)]/(mat1[EXM_mi(k,line1,row,row)]);
    }
    mat1[EXM_mi(k,line1,row,row)]=1.e0;
  }
}



static void mat_equal_tdma(double *mat1,long line1,double *mat2,long line2, long k){
  long row,col;
  for (row=0; row<k; row++){
    for (col=0; col<k; col++){
      mat2[EXM_mi(k,line2,row,col)]=mat1[EXM_mi(k,line1,row,col)];
    }
  }
}


static void find_multfact(double *mat1,long line1,double *mat2,long line2,
                   double *mat3,long line3, long k){
  long row,col;
/*
  Idea: mat3*mat1=-mat2
  or: mat3=-mat2*mat1_inv
*/
  mat_equal_tdma(mat1,line1,mat3,5,k);
  mat_inv_tdma(mat3,5,mat3,4,k);
  mat_mult_tdma(mat2,line2,mat3,4,mat3,line3,k);
  for (row=0; row<k; row++){
    for (col=0; col<k; col++){
      mat3[EXM_mi(k,line3,row,col)]=-mat3[EXM_mi(k,line3,row,col)];
    }
  }
}

/* solve block TDMA with first line at line=0 and last line at
   line=linemax */
void EXM_solve_block_TDMA(double *AA, double *BB, double *CC, double *RHS,
                    long linemax, long k){
  long line;
  double *TMP;

  TMP=(double *) malloc(6*k*k*sizeof(double));

/*  --------------------------------------------------------------
    Sweep Downward
    -------------------------------------------------------------- */
  for (line=0; line<linemax; line++){
    find_multfact(BB,line,AA,line+1,TMP,1,k);
    mat_mult_tdma(TMP,1,CC,line,TMP,2,k);
    mat_add_tdma(TMP,2,BB,line+1,TMP,3,k);
    mat_equal_tdma(TMP,3,BB,line+1,k);
    mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
    mat_add_tdma(TMP,2,RHS,line+1,TMP,3,k);
    mat_equal_tdma(TMP,3,RHS,line+1,k);
  }

/*   --------------------------------------------------------------
     Sweep Upward
     -------------------------------------------------------------- */
  for (line=linemax; line>0; line--){
    find_multfact(BB,line,CC,line-1,TMP,1,k);
    mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
    mat_add_tdma(TMP,2,RHS,line-1,TMP,3,k);
    mat_equal_tdma(TMP,3,RHS,line-1,k);
  }
/*   --------------------------------------------------------------
     Make BB identity
     -------------------------------------------------------------- */
  for (line=0; line<=linemax; line++){
    mat_inv_tdma(BB,line,TMP,1,k);
    mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
    mat_equal_tdma(TMP,2,RHS,line,k);
  }
  free(TMP);
}



void EXM_solve_block_TDMA_and_check(double *AA, double *BB, double *CC, double *DD,
                            long linemax, long k){
  long line,cnt;
  double *A,*B,*C,*D;
  double *TMP;
  bool PROBLEM;

  TMP=(double *) malloc(20*k*k*sizeof(double));

  A=(double *) malloc((linemax+2)*k*k*sizeof(double));
  B=(double *) malloc((linemax+2)*k*k*sizeof(double));
  C=(double *) malloc((linemax+2)*k*k*sizeof(double));
  D=(double *) malloc((linemax+2)*k*k*sizeof(double));

  for (line=0; line<=linemax; line++){
    mat_equal_tdma(AA,line,A,line,k);
    mat_equal_tdma(BB,line,B,line,k);
    mat_equal_tdma(CC,line,C,line,k);
    mat_equal_tdma(DD,line,D,line,k);
  }

  EXM_solve_block_TDMA(AA, BB, CC, DD, linemax, k);
  /* now, compare  A[line]*DD[line-1]+B[line]*DD[line]+C[line]*DD[line+1] and D[line] */
  for (line=0; line<=linemax; line++){
    mat_mult_tdma(B,line,DD,line,TMP,10,k);

    if (line!=linemax) {
      mat_mult_tdma(C,line,DD,line+1,TMP,11,k);
      mat_add_tdma(TMP,10,TMP,11,TMP,12,k);
      mat_equal_tdma(TMP,12,TMP,10,k);
    }

    if (line!=0) {
      mat_mult_tdma(A,line,DD,line-1,TMP,11,k);
      mat_add_tdma(TMP,10,TMP,11,TMP,12,k);
      mat_equal_tdma(TMP,12,TMP,10,k);
    }
    PROBLEM=FALSE;
    for (cnt=0; cnt<k; cnt++) if (fabs(TMP[EXM_mi(k,10,cnt,0)]-D[EXM_mi(k,line,cnt,0)])>1.0e-10)
      PROBLEM=TRUE;

    if (PROBLEM) {
      for (cnt=0; cnt<k; cnt++) printf("%E  ",TMP[EXM_mi(k,10,cnt,0)]);
      printf("\n");
      for (cnt=0; cnt<k; cnt++) printf("%E  ",D[EXM_mi(k,line,cnt,0)]);
      printf("\n\n");
    }
  }
  free(A);
  free(B);
  free(C);
  free(D);
  free(TMP);

}




/* solve block PDMA with first line at line=0 and last line at
   line=linemax */
void EXM_solve_block_PDMA(double *AA, double *BB, double *CC, double *DD,
                    double *EE, double *RHS, long linemax, long k){
  long line;
  double *TMP;

  TMP=(double *) malloc(6*k*k*sizeof(double));

/*     --------------------------------------------------------------
       Sweep Downward
       -------------------------------------------------------------- */
  for (line=1; line<linemax; line++){
/*       ------------First line Downwards------------------------ */
    find_multfact(CC,line,BB,line+1,TMP,1,k);
    mat_mult_tdma(TMP,1,DD,line,TMP,2,k);
    mat_add_tdma(TMP,2,CC,line+1,TMP,3,k);
    mat_equal_tdma(TMP,3,CC,line+1,k);
    mat_mult_tdma(TMP,1,EE,line,TMP,2,k);
    mat_add_tdma(TMP,2,DD,line+1,TMP,3,k);
    mat_equal_tdma(TMP,3,DD,line+1,k);
    mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
    mat_add_tdma(TMP,2,RHS,line+1,TMP,3,k);
    mat_equal_tdma(TMP,3,RHS,line+1,k);
/*       ------------Second Line Downwards----------------------- */
    if (line != linemax-1) {
      find_multfact(CC,line,AA,line+2,TMP,1,k);
      mat_mult_tdma(TMP,1,DD,line,TMP,2,k);
      mat_add_tdma(TMP,2,BB,line+2,TMP,3,k);
      mat_equal_tdma(TMP,3,BB,line+2,k);
      mat_mult_tdma(TMP,1,EE,line,TMP,2,k);
      mat_add_tdma(TMP,2,CC,line+2,TMP,3,k);
      mat_equal_tdma(TMP,3,CC,line+2,k);
      mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
      mat_add_tdma(TMP,2,RHS,line+2,TMP,3,k);
      mat_equal_tdma(TMP,3,RHS,line+2,k);
    }
  }



/*     --------------------------------------------------------------
       Sweep Upward
       -------------------------------------------------------------- */
  for (line=linemax; line>0; line--){
/*       -------First Line upwards------------------------------ */
    find_multfact(CC,line,DD,line-1,TMP,1,k);
    mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
    mat_add_tdma(TMP,2,RHS,line-1,TMP,3,k);
    mat_equal_tdma(TMP,3,RHS,line-1,k);
/*       -------Second Line Upwards----------------------------- */
    if (line != 1) {
      find_multfact(CC,line,EE,line-2,TMP,1,k);
      mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
      mat_add_tdma(TMP,2,RHS,line-2,TMP,3,k);
      mat_equal_tdma(TMP,3,RHS,line-2,k);
    }
  }


/*     --------------------------------------------------------------
       Make CC identity
       -------------------------------------------------------------- */
  for (line=0; line<=linemax; line++){
    mat_inv_tdma(CC,line,TMP,1,k);
    mat_mult_tdma(TMP,1,RHS,line,TMP,2,k);
    mat_equal_tdma(TMP,2,RHS,line,k);
  }
  free(TMP);
}


/* Other Math subroutines of interes*/

double sign_old(double x){
  double tmp;
  if (x<0) {
    tmp=-1.0e0;
  } else {
    tmp=1.0e0;
  }
  return(tmp);
}






long mod_old(long numb, long div){
  long tmp;
  long tmp2;
  long tmp3;
  tmp=numb/div;
  tmp2=tmp*div;
  tmp3=numb-tmp2;
  return (tmp3);
}

long mod_old2 (long a, long b)
{
   if(b < 0) //you can check for b == 0 separately and do what you want
     return mod(-a, -b);   
   long ret = a % b;
   if(ret < 0)
     ret+=b;
   return ret;
}


double powint(double x, long y){
  double sum;
  long cnt;
  sum=1.0e0;
  for (cnt=1; cnt<=labs(y); cnt++){
    sum=sum*x;
  }
  return(sum);
}

double krodelta(long i, long j){
  double tmp;
  if (i==j) {
    tmp=1.0e0;
  } else {
    tmp=0.0e0;
  }
  return(tmp);
}

void EXM_init_matrix(EXM_mat_t *mat, long numrow, long numcol){
  long row,col;
  mat->cont=(double *)malloc(numrow*numcol*sizeof(double));
  mat->glm.numcol=numcol;
  mat->glm.numrow=numrow;
  for (row=0; row<numrow; row++){
    for (col=0; col<numcol; col++){
      mat->cont[EXM_aim(mat->glm,row,col)]=0.0e0;
    }
  }
}


void EXM_reinit_matrix(EXM_mat_t *mat, long numrow, long numcol){
  long row,col;
  mat->cont=(double *)realloc(mat->cont,numrow*numcol*sizeof(double));
  mat->glm.numcol=numcol;
  mat->glm.numrow=numrow;
  for (row=0; row<numrow; row++){
    for (col=0; col<numcol; col++){
      mat->cont[EXM_aim(mat->glm,row,col)]=0.0e0;
    }
  }
}


void EXM_free_matrix(EXM_mat_t *mat){
  free(mat->cont);
}

void EXM_init_identity_matrix(EXM_mat_t *mat, long numrow, long numcol){
  long row,col;
  mat->cont=(double *)malloc(numrow*numcol*sizeof(double));
  mat->glm.numcol=numcol;
  mat->glm.numrow=numrow;
  for (row=0; row<numrow; row++){
    for (col=0; col<numcol; col++){
      if (row==col) {
        mat->cont[EXM_aim(mat->glm,row,col)]=1.0e0;
      } else {
        mat->cont[EXM_aim(mat->glm,row,col)]=0.0e0;
      }
    }
  }
}


void EXM_display_matrix(EXM_mat_t mat){
  long row,col;

  printf("numrow=%ld  numcol=%ld\n",mat.glm.numrow,mat.glm.numcol);
  for (row=0; row<mat.glm.numrow; row++){
    for (col=0; col<mat.glm.numcol; col++){
      printf("%E  ",mat.cont[EXM_aim(mat.glm,row,col)]);
    }
    printf("\n");
  }
  printf("\n");
}

void EXM_multiply_matrices(EXM_mat_t mat1, EXM_mat_t mat2, EXM_mat_t *matr){
  long cnt,numrow,numcol,row,col;
  numrow=mat1.glm.numrow;
  numcol=mat2.glm.numcol;
  EXM_reinit_matrix(matr,numrow,numcol);
  if (mat2.glm.numrow==mat1.glm.numcol) {
    for (row=0; row<numrow; row++){
      for (col=0; col<numcol; col++){
         matr->cont[EXM_aim(matr->glm,row,col)]=0.0e0;
         for (cnt=0; cnt<mat1.glm.numcol; cnt++){
           matr->cont[EXM_aim(matr->glm,row,col)]=matr->cont[EXM_aim(matr->glm,row,col)]
                               +mat1.cont[EXM_aim(mat1.glm,row,cnt)]
                               *mat2.cont[EXM_aim(mat2.glm,cnt,col)];
         }
      }
    }
  } else {
    printf("cannot multiply the two matrices\n");
  } 


}



void EXM_invert_matrix_gaussian_elimination(EXM_mat_t mat, EXM_mat_t *matinv){
  long row,col,row2;
  double fact;
  EXM_mat_t mattmp;
  EXM_init_matrix(&mattmp,mat.glm.numrow,mat.glm.numcol);
  if (mat.glm.numrow==mat.glm.numcol) {
    EXM_reinit_matrix(matinv,mat.glm.numrow,mat.glm.numcol);
    for (row=0; row<mat.glm.numrow; row++){
      for (col=0; col<mat.glm.numcol; col++){
        mattmp.cont[EXM_aim(mat.glm,row,col)]=mat.cont[EXM_aim(mat.glm,row,col)];
        matinv->cont[EXM_aim(mat.glm,row,col)]=0.0;
      }
      matinv->cont[EXM_aim(mat.glm,row,row)]=1.0;
    }
// make the non-diagonal elements zero for mattmp 
    for (row=0; row<mat.glm.numrow; row++){
      for (row2=0; row2<mat.glm.numrow; row2++){
        if (row2!=row) {
					if (mattmp.cont[EXM_aim(mat.glm,row,row)]==0.0){
            printf("matrix cannot be inverted\n");
          }          
          fact=-mattmp.cont[EXM_aim(mat.glm,row2,row)]/mattmp.cont[EXM_aim(mat.glm,row,row)];
          for (col=0; col<mat.glm.numcol; col++){
            mattmp.cont[EXM_aim(mat.glm,row2,col)]+=fact*mattmp.cont[EXM_aim(mat.glm,row,col)];
            matinv->cont[EXM_aim(mat.glm,row2,col)]+=fact*matinv->cont[EXM_aim(mat.glm,row,col)];
          }
        }
      }    
    }
//   EXM_display_matrix(mattmp);
// make the diagonal elements equal to 1 for mattmp
    for (row=0; row<mat.glm.numrow; row++){
		  if (mattmp.cont[EXM_aim(mat.glm,row,row)]==0.0){
        printf("matrix cannot be inverted\n");
      }          
      fact=1.0/mattmp.cont[EXM_aim(mat.glm,row,row)];
      for (col=0; col<mat.glm.numcol; col++){
        mattmp.cont[EXM_aim(mat.glm,row,col)]*=fact;
        matinv->cont[EXM_aim(mat.glm,row,col)]*=fact;
      }
    }
  } else {
    printf("matrix cannot be inverted\n");
    printf("number of rows not equal to number of columns \n");
  }
  EXM_free_matrix(&mattmp);
}


// function by Jaehyuk Lee 
void EXM_invert_matrix_partial_pivoting(EXM_mat_t mat, EXM_mat_t *matinv){
  long pivot,row,row2,col;
  double temp,pivotval,fact;
  EXM_mat_t mattmp;
  EXM_init_matrix(&mattmp,mat.glm.numrow,mat.glm.numcol);
  if (mat.glm.numrow==mat.glm.numcol) {
    EXM_reinit_matrix(matinv,mat.glm.numrow,mat.glm.numcol);
    for (row=0; row<mat.glm.numrow; row++){
      for (col=0; col<mat.glm.numcol; col++){
        mattmp.cont[EXM_aim(mat.glm,row,col)]=mat.cont[EXM_aim(mat.glm,row,col)];
        matinv->cont[EXM_aim(mat.glm,row,col)]=0.0;
      }
      matinv->cont[EXM_aim(mat.glm,row,row)]=1.0;
    }
    //process starts from here
    for(row=0;row<mat.glm.numrow-1;row++){
      //partial pivoting start
      pivot=row;
      pivotval=mattmp.cont[EXM_aim(mat.glm,row,row)];
      for(row2=row+1;row2<mat.glm.numrow;row2++){
        if(fabs(pivotval)<fabs(mattmp.cont[EXM_aim(mat.glm,row2,row)])){
          pivot=row2;
          pivotval=mattmp.cont[EXM_aim(mat.glm,row2,row)];
        }
      }
      if(pivot!=row){
        for(col=0;col<mat.glm.numrow;col++){
          //partial pivoting of mattmp
          temp=mattmp.cont[EXM_aim(mat.glm,pivot,col)];
          mattmp.cont[EXM_aim(mat.glm,pivot,col)]=mattmp.cont[EXM_aim(mat.glm,row,col)];
          mattmp.cont[EXM_aim(mat.glm,row,col)]=temp;
          //partial pivoting of matinv
          temp=matinv->cont[EXM_aim(mat.glm,pivot,col)];
          matinv->cont[EXM_aim(mat.glm,pivot,col)]=matinv->cont[EXM_aim(mat.glm,row,col)];
          matinv->cont[EXM_aim(mat.glm,row,col)]=temp;
        }
      }//partial pivoting complete
      //forward substitution start
      for(row2=row+1;row2<mat.glm.numrow;row2++){
        assert(mattmp.cont[EXM_aim(mat.glm,row,row)]!=0.0e0);
        fact=-(mattmp.cont[EXM_aim(mat.glm,row2,row)])/(mattmp.cont[EXM_aim(mat.glm,row,row)]);
        for(col=0;col<mat.glm.numrow;col++){
          matinv->cont[EXM_aim(mat.glm,row2,col)]=matinv->cont[EXM_aim(mat.glm,row2,col)]+matinv->cont[EXM_aim(mat.glm,row,col)]*fact;
        }
        for(col=row;col<mat.glm.numrow;col++){
          mattmp.cont[EXM_aim(mat.glm,row2,col)]=mattmp.cont[EXM_aim(mat.glm,row2,col)]+mattmp.cont[EXM_aim(mat.glm,row,col)]*fact;
        }
        mattmp.cont[EXM_aim(mat.glm,row2,row)]=0.0e0;
      }//forward substitution complete
    }
    //backward substitution start
    for(row=mat.glm.numrow-1;row>=1;row--){
      //multiply fact over matinv
      for(row2=row-1;row2>=0;row2--){
        assert(mattmp.cont[EXM_aim(mat.glm,row,row)]!=0.0e0);
        fact=-mattmp.cont[EXM_aim(mat.glm,row2,row)]/mattmp.cont[EXM_aim(mat.glm,row,row)];
        for(col=0;col<mat.glm.numrow;col++){
          matinv->cont[EXM_aim(mat.glm,row2,col)]=matinv->cont[EXM_aim(mat.glm,row2,col)]+matinv->cont[EXM_aim(mat.glm,row,col)]*fact;
        }
      }
    }
    //devide matinv by diagonal of mattmp
    for(row=0;row<mat.glm.numrow;row++){
      for(col=0;col<mat.glm.numrow;col++){
        assert(mattmp.cont[EXM_aim(mat.glm,row,row)]!=0.0e0);
        matinv->cont[EXM_aim(mat.glm,row,col)]=matinv->cont[EXM_aim(mat.glm,row,col)]/mattmp.cont[EXM_aim(mat.glm,row,row)];
      }
    }//backward substitution complete
    //EXM_display_matrix(mattmp);
  }
  else{
    printf("matrix cannot be inverted\n");
    printf("number of rows not equal to number of columns \n");
  }
  EXM_free_matrix(&mattmp);
}


void EXM_invert_matrix(EXM_mat_t mat, EXM_mat_t *matinv){
  //EXM_invert_matrix_gaussian_elimination(mat,matinv);
  EXM_invert_matrix_partial_pivoting(mat,matinv);
}


void EXM_invert_matrix_analytical(EXM_mat_t mat, EXM_mat_t *matinv){
  double den;
  
  if ((mat.glm.numrow==mat.glm.numcol) && (mat.glm.numrow<4)) {
    EXM_reinit_matrix(matinv,mat.glm.numrow,mat.glm.numcol);
    if (mat.glm.numrow==1){
      den=mat.cont[EXM_aim(mat.glm,0,0)];
      assert(den!=0.0);
      matinv->cont[EXM_aim(matinv->glm,0,0)]=1.0e0/den;
    }
    if (mat.glm.numrow==2){
      den=mat.cont[EXM_aim(mat.glm,0,0)]*mat.cont[EXM_aim(mat.glm,1,1)]
         -mat.cont[EXM_aim(mat.glm,0,1)]*mat.cont[EXM_aim(mat.glm,1,0)];
      assert(den!=0.0);
      matinv->cont[EXM_aim(matinv->glm,0,0)]=mat.cont[EXM_aim(mat.glm,1,1)]/den;
      matinv->cont[EXM_aim(matinv->glm,0,1)]=-mat.cont[EXM_aim(mat.glm,0,1)]/den;
      matinv->cont[EXM_aim(matinv->glm,1,0)]=-mat.cont[EXM_aim(mat.glm,1,0)]/den;
      matinv->cont[EXM_aim(matinv->glm,1,1)]=mat.cont[EXM_aim(mat.glm,0,0)]/den;
    }
    if (mat.glm.numrow==3){
      den=mat.cont[EXM_aim(mat.glm,0,0)]*mat.cont[EXM_aim(mat.glm,1,1)]*mat.cont[EXM_aim(mat.glm,2,2)]
         -mat.cont[EXM_aim(mat.glm,0,0)]*mat.cont[EXM_aim(mat.glm,1,2)]*mat.cont[EXM_aim(mat.glm,2,1)]
         -mat.cont[EXM_aim(mat.glm,1,0)]*mat.cont[EXM_aim(mat.glm,0,1)]*mat.cont[EXM_aim(mat.glm,2,2)]
         +mat.cont[EXM_aim(mat.glm,1,0)]*mat.cont[EXM_aim(mat.glm,0,2)]*mat.cont[EXM_aim(mat.glm,2,1)]
         +mat.cont[EXM_aim(mat.glm,2,0)]*mat.cont[EXM_aim(mat.glm,0,1)]*mat.cont[EXM_aim(mat.glm,1,2)]
         -mat.cont[EXM_aim(mat.glm,2,0)]*mat.cont[EXM_aim(mat.glm,0,2)]*mat.cont[EXM_aim(mat.glm,1,1)];
      assert(den!=0.0);
      matinv->cont[EXM_aim(mat.glm,0,0)]=(mat.cont[EXM_aim(mat.glm,1,1)]*mat.cont[EXM_aim(mat.glm,2,2)]-mat.cont[EXM_aim(mat.glm,1,2)]*mat.cont[EXM_aim(mat.glm,2,1)])/den;
      matinv->cont[EXM_aim(mat.glm,0,1)]=(mat.cont[EXM_aim(mat.glm,2,1)]*mat.cont[EXM_aim(mat.glm,0,2)]-mat.cont[EXM_aim(mat.glm,2,2)]*mat.cont[EXM_aim(mat.glm,0,1)])/den;
      matinv->cont[EXM_aim(mat.glm,0,2)]=(mat.cont[EXM_aim(mat.glm,0,1)]*mat.cont[EXM_aim(mat.glm,1,2)]-mat.cont[EXM_aim(mat.glm,0,2)]*mat.cont[EXM_aim(mat.glm,1,1)])/den;

      matinv->cont[EXM_aim(mat.glm,1,0)]=(mat.cont[EXM_aim(mat.glm,1,2)]*mat.cont[EXM_aim(mat.glm,2,0)]-mat.cont[EXM_aim(mat.glm,1,0)]*mat.cont[EXM_aim(mat.glm,2,2)])/den;
      matinv->cont[EXM_aim(mat.glm,1,1)]=(mat.cont[EXM_aim(mat.glm,2,2)]*mat.cont[EXM_aim(mat.glm,0,0)]-mat.cont[EXM_aim(mat.glm,2,0)]*mat.cont[EXM_aim(mat.glm,0,2)])/den;
      matinv->cont[EXM_aim(mat.glm,1,2)]=(mat.cont[EXM_aim(mat.glm,0,2)]*mat.cont[EXM_aim(mat.glm,1,0)]-mat.cont[EXM_aim(mat.glm,0,0)]*mat.cont[EXM_aim(mat.glm,1,2)])/den;

      matinv->cont[EXM_aim(mat.glm,2,0)]=(mat.cont[EXM_aim(mat.glm,1,0)]*mat.cont[EXM_aim(mat.glm,2,1)]-mat.cont[EXM_aim(mat.glm,1,1)]*mat.cont[EXM_aim(mat.glm,2,0)])/den;
      matinv->cont[EXM_aim(mat.glm,2,1)]=(mat.cont[EXM_aim(mat.glm,2,0)]*mat.cont[EXM_aim(mat.glm,0,1)]-mat.cont[EXM_aim(mat.glm,2,1)]*mat.cont[EXM_aim(mat.glm,0,0)])/den;
      matinv->cont[EXM_aim(mat.glm,2,2)]=(mat.cont[EXM_aim(mat.glm,0,0)]*mat.cont[EXM_aim(mat.glm,1,1)]-mat.cont[EXM_aim(mat.glm,0,1)]*mat.cont[EXM_aim(mat.glm,1,0)])/den;
    }    
  } else {
    printf("mat.contrix cannot be inverted\n");
    printf("number of rows not equal to number of columns \n");
  }
}
  

double rad(double angle){
 double tmp;
 tmp=angle*pi/180.0e0;
 return(tmp);
}

double deg(double angle){
 double tmp;
 tmp=angle/pi*180.0e0;
 return(tmp);
}


/* find orthogonal vector to the three points pa,
   pb and pc in the xyz frame of reference */
void EXM_find_orthogonal_vector(EXM_vec3D_t pa, EXM_vec3D_t pb, EXM_vec3D_t pc,
                    EXM_vec3D_t orthovect){
  EXM_vec3D_t dp1,dp2;
  long cnt;

  for (cnt=0; cnt<3; cnt++){
    dp1[cnt]=pa[cnt]-pc[cnt];
    dp2[cnt]=pb[cnt]-pc[cnt];
  }
  orthovect[0]=dp1[1]*dp2[2]-dp1[2]*dp2[1];
  orthovect[1]=dp1[2]*dp2[0]-dp1[0]*dp2[2];
  orthovect[2]=dp1[0]*dp2[1]-dp1[1]*dp2[0];
}


/* find plane A*x+B*y+C*z=D from orthogonal vector
   and one point on plane p0 */
void EXM_find_plane(EXM_vec3D_t orthovect, EXM_vec3D_t p0,
               double *A, double *B, double *C, double *D){
  *A=orthovect[0];
  *B=orthovect[1];
  *C=orthovect[2];
  *D=(*A)*p0[0]+(*B)*p0[1]+(*C)*p0[2];
}

/* find pp, a point on the plane A-B-C-D which also
   lies on the line composed of vect and p0 */
void EXM_find_point_in_plane_on_vector(EXM_vec3D_t vect, EXM_vec3D_t p0,
               double A, double B, double C, double D,
               EXM_vec3D_t pp){
  double t;
  long cnt;
  assert(A*vect[0]+B*vect[1]+C*vect[2]!=0.0e0);
  t=(D-A*p0[0]-B*p0[1]-C*p0[2])/(A*vect[0]+B*vect[1]+C*vect[2]);
  for (cnt=0; cnt<3; cnt++)
    pp[cnt]=p0[cnt]+vect[cnt]*t;
}


/* the plane is defined by points pa,pb and pc while
   the point to be mirrored is pp_o. The
   mirrored point is pp_m; pp_p is the point on the plane
   nearest to pp_o, midway between pp_o and pp_m */
void EXM_mirror_point_wrt_plane(EXM_vec3D_t pa, EXM_vec3D_t pb,
                         EXM_vec3D_t pc, EXM_vec3D_t pp_o,
                         EXM_vec3D_t pp_m){
  EXM_vec3D_t orthovect;
  double A,B,C,D;
  EXM_vec3D_t pp_p;
  long cnt;

  EXM_find_orthogonal_vector(pa,pb,pc,orthovect);
  EXM_find_plane(orthovect, pa, &A, &B, &C, &D);
  EXM_find_point_in_plane_on_vector(orthovect, pp_o, A, B, C, D, pp_p);
  for (cnt=0; cnt<3; cnt++)
    pp_m[cnt]=2.0e0*pp_p[cnt]-pp_o[cnt];
}


void EXM_display_vector(EXM_vec3D_t pp){
  printf("x=%E  y=%E  z=%E \n",pp[0],pp[1],pp[2]);
}


double EXM_dot_product(EXM_vec3D_t vec1, EXM_vec3D_t vec2){
  double sum;
  long dim;

  sum=0.0e0;
  for (dim=0; dim<3; dim++){
    sum=sum+vec1[dim]*vec2[dim];
  }
  return(sum);
}


void EXM_cross_product(EXM_vec3D_t vec1, EXM_vec3D_t vec2, EXM_vec3D_t prod){

  prod[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
  prod[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
  prod[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
}


double EXM_vector_magnitude(EXM_vec3D_t vec){
  long dim;
  double sum;
  sum=0.0e0;
  for (dim=0; dim<3; dim++){
    sum=sum+vec[dim]*vec[dim];
  }
  sum=sqrt(sum);
  return(sum);
}


double EXM_angle_between_vectors(EXM_vec3D_t vec1, EXM_vec3D_t vec2){
  double theta;
  EXM_vec3D_t vec1norm,vec2norm;
  EXM_normalize_vector(vec1, vec1norm);
  EXM_normalize_vector(vec2, vec2norm);
  theta=acos(EXM_dot_product(vec1norm,vec2norm));
  return(theta);
}


void EXM_normalize_vector(EXM_vec3D_t vec, EXM_vec3D_t vecnorm){
  long dim;
  double vecmag;
  vecmag=EXM_vector_magnitude(vec);
  assert(vecmag!=0.0);
  for (dim=0; dim<3; dim++) vecnorm[dim]=vec[dim]/vecmag;
}


/* rotation matrix that rotates vec1 into vec2 such that vec2norm=R * vec1norm */
void EXM_find_rotation_matrix(EXM_vec3D_t vec1, EXM_vec3D_t vec2, EXM_mat3x3_t R){
  double theta,costheta,sintheta;
  EXM_vec3D_t s,snorm;
  theta=EXM_angle_between_vectors(vec1, vec2);
  EXM_cross_product(vec1, vec2, s);
  if (EXM_vector_magnitude(s)!=0.0) {
    EXM_normalize_vector(s,snorm);
  } else {
    snorm[0]=0.0; 
    snorm[1]=0.0; 
    snorm[2]=0.0;
  }
  costheta=cos(theta);
  sintheta=sin(theta);
  
  // from https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  // https://math.stackexchange.com/questions/4025666/how-to-rotate-a-vector-in-3d-space-around-arbitrary-axis
  R[0][0]=costheta+sqr(snorm[0])*(1.0-costheta);
  R[0][1]=snorm[0]*snorm[1]*(1.0-costheta)-snorm[2]*sintheta;
  R[0][2]=snorm[0]*snorm[2]*(1.0-costheta)+snorm[1]*sintheta;
  
  R[1][0]=snorm[1]*snorm[0]*(1.0-costheta)+snorm[2]*sintheta;
  R[1][1]=costheta+sqr(snorm[1])*(1.0-costheta);
  R[1][2]=snorm[1]*snorm[2]*(1.0-costheta)-snorm[0]*sintheta;
  
  R[2][0]=snorm[2]*snorm[0]*(1.0-costheta)-snorm[1]*sintheta;
  R[2][1]=snorm[2]*snorm[1]*(1.0-costheta)+snorm[0]*sintheta;
  R[2][2]=costheta+sqr(snorm[2])*(1.0-costheta);
}


void EXM_multiply_matrix_vector(EXM_mat3x3_t mat, EXM_vec3D_t vec, EXM_vec3D_t res){
  long row,col;
  for (row=0; row<3; row++){
    res[row]=0.0;
    for (col=0; col<3; col++){
      res[row]+=mat[row][col]*vec[col];
    }
  }
}


double EXM_area_quadrilateral(EXM_vec3D_t A, EXM_vec3D_t B, EXM_vec3D_t C, EXM_vec3D_t D){
  EXM_vec3D_t x,y,z;
  double theta,phi,area;
  long dim;
  for (dim=0; dim<3; dim++){
    x[dim]=A[dim]-B[dim];
    y[dim]=C[dim]-B[dim];
    z[dim]=D[dim]-B[dim];
  }
  theta=EXM_angle_between_vectors(y,z);
  phi=EXM_angle_between_vectors(x,y);
  area=0.5e0*EXM_vector_magnitude(y)*(EXM_vector_magnitude(z)*fabs(sin(theta))+EXM_vector_magnitude(x)*fabs(sin(phi)));
  return(area);
}


/* numerically differentiate FUNCT (which returns dx/dt)
   from t1 to t2, starting from x1, and returning x2 */
double EXM_numerical_differentiation(double(*FUNCT)(void *, double, double), void *arg_ptr,
                 long METHOD, long n, double x1, double t1, double t2, long *error){
  double x,dt,t;
  double x_1,x_2,x_3;
  long cnt;

  *error=0;
  dt=(t2-t1)/(double)n;
  x=x1;
  t=t1;
  if (METHOD==EXM_NUMDIFF_FORWARDEULER)
    for (cnt=0; cnt<n; cnt++) {
      x+=dt*(*FUNCT)(arg_ptr,x,t);
      t+=dt;
    }
  if (METHOD==EXM_NUMDIFF_RUNGEKUTTA)
    for (cnt=0; cnt<n; cnt++) {
      x_1=x+0.5*dt*(*FUNCT)(arg_ptr,x,t);
      x_2=x+0.5*dt*(*FUNCT)(arg_ptr,x_1,t+0.5*dt);
      x_3=x+dt*(*FUNCT)(arg_ptr,x_2,t+0.5*dt);
      x+=(x_1-x)/3.0
        +(x_2-x)*2.0/3.0
        +(x_3-x)/3.0
        +dt/6.0*(*FUNCT)(arg_ptr,x_2,t+dt);
      t+=dt;
    }
  if (METHOD==EXM_NUMDIFF_IMPROVEDEULER)
    for (cnt=0; cnt<n; cnt++) {
      x+=0.5*dt*(*FUNCT)(arg_ptr,x,t)
        +0.5*dt*(*FUNCT)(arg_ptr,x+dt*(*FUNCT)(arg_ptr,x,t),t+dt);
      t+=dt;
    }
  if (METHOD==EXM_NUMDIFF_MODIFIEDEULER)
    for (cnt=0; cnt<n; cnt++) {
      x+=dt*(*FUNCT)(arg_ptr,x+0.5e0*dt*(*FUNCT)(arg_ptr,x,t),t+dt*0.5);
      t+=dt;
    }
  return(x);
}


double EXM_numerical_integration(double(*FUNCT)(void *, double), void *arg_ptr,
                    long METHOD, long n, double x1, double x2, long *error){
  double f,f1,f2,f3,x,sum,dx;

  sum=0.0;
  *error=0;
  dx=(x2-x1)/(double)n;
  if (METHOD==EXM_NUMINTEG_RECTANGLES) {
    sum=0.0e0;
    for (x=x1+dx*0.5e0; x<x2; x+=dx){
      f=(*FUNCT)(arg_ptr,x);
      sum+=f*dx;
    }
  }
  if (METHOD==EXM_NUMINTEG_POLY2) { /* integrate using second degree polynomials */
    f1=(*FUNCT)(arg_ptr,x1);
    f2=(*FUNCT)(arg_ptr,x1+dx);
    f3=(*FUNCT)(arg_ptr,x1+2.0*dx);
    sum=dx*(27.0/24.0*f2+9.0/24.0*f1);
    for (x=x1+2.0*dx; x<x2-1.5*dx; x+=dx){
      f1=f2;
      f2=f3;
      f3=(*FUNCT)(arg_ptr,x+dx);
      sum+=dx*(f1/24.0+11.0/12.0*f2+f3/24.0);
    }
    f1=f2;
    f2=f3;
    f3=(*FUNCT)(arg_ptr,x2);
    sum+=dx*(27.0/24.0*f2+9.0/24.0*f3);
  }
  return(sum);
}


/* numerically integrate the function FUNCT(x) between
   x=x1 and x=x2, returning the value of the integral */
void EXM_numerical_integration_vector(void(*FUNCT)(void *, double, EXM_vec3D_t vector), void *arg_ptr,
                        long METHOD, long n, double x1, double x2, long *error,
                        EXM_vec3D_t sum){
  double x,dx;
  long dim;
  EXM_vec3D_t f,f1,f2,f3;

  for (dim=0; dim<3; dim++) sum[dim]=0.0;
  *error=0;
  dx=(x2-x1)/(double)n;
  if (METHOD==EXM_NUMINTEG_RECTANGLES) {
    for (dim=0; dim<3; dim++) sum[dim]=0.0;
    for (x=x1+dx*0.5e0; x<x2; x+=dx){
      (*FUNCT)(arg_ptr,x,f);
      for (dim=0; dim<3; dim++) sum[dim]+=f[dim]*dx;
    }
  }
  if (METHOD==EXM_NUMINTEG_POLY2) { /* integrate using second degree polynomials */
    (*FUNCT)(arg_ptr,x1,f1);
    (*FUNCT)(arg_ptr,x1+dx,f2);
    (*FUNCT)(arg_ptr,x1+2.0*dx,f3);
    for (dim=0; dim<3; dim++) sum[dim]=dx*(27.0/24.0*f2[dim]+9.0/24.0*f1[dim]);
    for (x=x1+2.0*dx; x<x2-1.5*dx; x+=dx){
      for (dim=0; dim<3; dim++) {
        f1[dim]=f2[dim];
        f2[dim]=f3[dim];
      }
      (*FUNCT)(arg_ptr,x+dx,f3);
      for (dim=0; dim<3; dim++) sum[dim]+=dx*(f1[dim]/24.0+11.0/12.0*f2[dim]+f3[dim]/24.0);
    }
    for (dim=0; dim<3; dim++){
      f1[dim]=f2[dim];
      f2[dim]=f3[dim];
    }
    (*FUNCT)(arg_ptr,x2,f3);
    for (dim=0; dim<3; dim++) sum[dim]+=dx*(27.0/24.0*f2[dim]+9.0/24.0*f3[dim]);
  }
}


/* returns f at thisx given x,f at each data point */
double EXM_f_from_line(long N, double *x, double *f, double thisx){
  long i;
  double thisf,dxi;
  if (N<2) EXM_fatal_error("Number of data points supplied for linear interpolation must be at least 2.");
  /* first find i that is such that x[i]<=thisx<=x[i+1] */
  i=-1;
  do {
    i++;
  } while(i<(N-1) && !(x[i]<=thisx && x[i+1]>=thisx));
  if (i>=N-1){
    EXM_fatal_error("Couldn't find an interval for x=%E in EXM_f_from_spline.",thisx);
  }
  assert(i<N-1);
  dxi=x[i+1]-x[i];
  thisf=f[i]+((f[i+1]-f[i])/dxi)*(thisx-x[i]);
  return(thisf);
}


/* find b at each node given f and x on all nodes from i=0 to i=N-1 
   note: *b must have been malloced and given enough memory space prior to calling this function*/
void EXM_find_spline(long N, double *x, double *f, double *b){
  EXM_pdmaline_t *pdma;
  double dxim1,dxip0,dfim1,dfip0;
  long n;
  pdma=(EXM_pdmaline_t *)malloc(N*sizeof(EXM_pdmaline_t));
  /* do inner nodes first */
  for (n=1; n<N-1; n++){
    dxim1=x[n]-x[n-1];
    dxip0=x[n+1]-x[n];
    dfim1=f[n]-f[n-1];
    dfip0=f[n+1]-f[n];
    pdma[n].val[0]=0.0;
    pdma[n].val[1]=dxim1;
    pdma[n].val[2]=2.0*(dxim1+dxip0);
    pdma[n].val[3]=dxip0;
    pdma[n].val[4]=0.0;
    assert(dxip0!=0.0);
    assert(dxim1!=0.0);
    pdma[n].val[5]=3.0*(dfip0/dxip0-dfim1/dxim1);
  }
  /* do left bdry node */
  pdma[0].val[0]=0.0;
  pdma[0].val[1]=0.0;
  pdma[0].val[2]=-(x[2]-x[1]);
  pdma[0].val[3]=(x[1]-x[0])+(x[2]-x[1]);
  pdma[0].val[4]=-(x[1]-x[0]);
  pdma[0].val[5]=0.0;
  /* do right bdry node */
  pdma[N-1].val[0]=-(x[N-1]-x[N-2]);
  pdma[N-1].val[1]=(x[N-2]-x[N-3])+(x[N-1]-x[N-2]);
  pdma[N-1].val[2]=-(x[N-2]-x[N-3]);
  pdma[N-1].val[3]=0.0;
  pdma[N-1].val[4]=0.0;
  pdma[N-1].val[5]=0.0;

  EXM_solve_PDMA(pdma, N);
  for (n=0; n<N; n++) {
    assert(pdma[n].val[2]!=0.0);
    b[n]=pdma[n].val[5]/pdma[n].val[2];
  }
  free(pdma);
}


/* returns f at thisx given x,f,b at each data point */
double EXM_f_from_spline(long N, double *x, double *f, double *b, double thisx){
  long i;
  double thisf,di,ci,bi,ai,dxi;
  /* first find i that is such that x[i]<=thisx<=x[i+1] */
  i=-1;
  do {
    i++;
  } while(i<(N-1) && !(x[i]<=thisx && x[i+1]>=thisx));
  if (i>=N-1){
    EXM_fatal_error("Couldn't find an interval for x=%E in EXM_f_from_spline.",thisx);
  }
  assert(i<N-1);
  dxi=x[i+1]-x[i];
  ai=(b[i+1]-b[i])/(3.0*dxi);
  bi=b[i];
  ci=(f[i+1]-f[i])/dxi-b[i]*dxi-(b[i+1]-b[i])/3.0*dxi;
  di=f[i];
  thisf=ai*(thisx-x[i])*(thisx-x[i])*(thisx-x[i])+bi*(thisx-x[i])*(thisx-x[i])+ci*(thisx-x[i])+di;
  return(thisf);
}

/* returns f at thisx given x,f,b at each data point */
double EXM_f_from_monotonespline(long N, double *x, double *f, double thisx){
  long i;
  double thisf,dx,t,deltam1,deltap0,deltap1,mp0,mp1,alpha,beta,tau;
  mp0=mp1=0.0; //to avoid compiler warning
  /* first find i that is such that x[i]<=thisx<=x[i+1] */
  i=-1;
  do {
    i++;
  } while(i<(N-1) && !(x[i]<=thisx && x[i+1]>=thisx));
  if (i>=N-1){
    EXM_fatal_error("Couldn't find an interval for x=%E in EXM_f_from_spline.",thisx);
  }
  assert(i<N-1);
  dx=x[i+1]-x[i];
  t=(thisx-x[i])/dx;
  deltap0=(f[i+1]-f[i])/dx;
  if(i>=1 && i<=N-3) {
    deltam1=(f[i]-f[i-1])/(x[i]-x[i-1]); 
    mp0=0.5*(deltam1+deltap0); 
    if(deltam1*deltap0<0.0) mp0=0.0; 
    deltap1=(f[i+2]-f[i+1])/(x[i+2]-x[i+1]); 
    mp1=0.5*(deltap0+deltap1); 
    if(deltap0*deltap1<0.0) mp1=0.0;
  }
  else if (i==0) {
    deltap1=(f[i+2]-f[i+1])/(x[i+2]-x[i+1]);
    mp0=deltap0;
    mp1=0.5*(deltap0+deltap1); 
    if(deltap0*deltap1<0.0) mp1=0.0; 
  }
  else if (i==N-2) {
    deltam1=(f[i]-f[i-1])/(x[i]-x[i-1]); 
    mp0=0.5*(deltam1+deltap0); 
    if(deltam1*deltap0<0.0) mp0=0.0; 
    mp1=deltap0;
  }
  else EXM_fatal_error("Input to EXM_f_from_monotonespline() out of range.");
  alpha=mp0/deltap0;
  beta=mp1/deltap0;
  if(alpha*alpha+beta*beta>9.0) {
    tau=3.0/sqrt(alpha*alpha+beta*beta);
    mp0=tau*alpha*deltap0;
    mp1=tau*beta*deltap0;
  }
  thisf=f[i]*(2.0*t*t*t-3.0*t*t+1)+dx*mp0*(t*t*t-2.0*t*t+t)+f[i+1]*(-2.0*t*t*t+3.0*t*t)+dx*mp1*(t*t*t-t*t);
  return(thisf);
}





#define EOS 0

/* insert str1 into str2 at the position pos; 
   make sure *str2 has enough memory allocated */
char *strins(char *str1, char *str2,  long pos){
  long len1,len2,i;
  len1=(long)strlen(str1);
  len2=(long)strlen(str2);

  for (i=len2; i>=pos; i--) (str2)[i+len1]=(str2)[i];
  for (i=0; i<len1; i++) (str2)[pos+i]=str1[i];
  (str2)[len2+len1]=EOS;
  return str2;
}


/* add line breaks without breaking words with width the maximum number of characters per line */
char *strwrp(char *str, int width){
  long cnt,cntbreak;
  bool CONTINUE;

  CONTINUE=TRUE;
  cntbreak=0;
  cnt=0;
  do {
    cnt++;
    if (str[cnt]=='\n') cntbreak=cnt;
    if (cnt-cntbreak>width){
      cntbreak=cnt;
      do {
        cntbreak--;
      } while(str[cntbreak]!=' ' && str[cntbreak]!='-' && cntbreak>0);
      if (cntbreak>0){
        if (str[cntbreak]=='-') {
          strins("\n",str,cntbreak+1);
          cntbreak++;
        }
        str[cntbreak]='\n';
        cnt+=1;
      } else {
        // problem breaking line..
        CONTINUE=FALSE;
      }
    }        
  } while (CONTINUE && cnt<strlen(str)-1); 
  return str;
}


char *strrep(char *str, char *orig, char *repl)
{
    char *p;
    char *strtmp;
    strtmp=(char *)malloc(sizeof(char)*(strlen(str)+strlen(repl)+2));
    if(!(p = strstr(str, orig)))
        return str;
    strncpy(strtmp, str, p-str);
    strtmp[p-str] = EOS;
    sprintf(strtmp+(p-str), "%s%s", repl, p+strlen(orig));
    strcpy(str,strtmp);
    free(strtmp);
    return str;
}



/* add indent spaces to second and subsequent lines; make sure str has enough memory allocated */
char *strind(char *str, int indent){
  long cnt,cnt2;
  static char whitespace[2];
  strcpy(whitespace," ");
  if (indent<0){
    /* hang indent */
    cnt=0;
    do {
      cnt++;
      if (str[cnt]=='\n' && str[cnt+1]!=EOS) {
        for (cnt2=0; cnt2<(-indent); cnt2++) strins(whitespace,str,cnt+1);
        cnt-=indent;
      }
    } while (str[cnt]!=EOS);
  } else {
    /* indent */
    for (cnt=0; cnt<indent; cnt++) strins(whitespace,str,0);
  }
  return str;
}


/* add line breaks without breaking words with width the maximum number of characters per line 
   and indent the number of indented characters (either negative or positive)  */
char *strwrpind(char *str, int width, int indent){
  long cnt,cnt2,cntbreak;
  bool CONTINUE;
  static char whitespace[2];
  strcpy(whitespace," ");

  if (indent>0){
    for (cnt=0; cnt<indent; cnt++) strins(whitespace,str,0);
  }

  CONTINUE=TRUE;
  cntbreak=0;
  cnt=0;
  do {
    cnt++;
    if (str[cnt]=='\n') cntbreak=cnt;
    if (cnt-cntbreak>width-1){
      cntbreak=cnt;
      do {
        cntbreak--;
      } while(str[cntbreak]!=' ' && str[cntbreak]!='-' && cntbreak>0);
      if (cntbreak>0){
        if (str[cntbreak]=='-') {
          strins("\n",str,cntbreak+1);
          cntbreak++;
        }
        str[cntbreak]='\n';
        cnt++;
        if (indent<0){
          // do the hang indent here
          for (cnt2=0; cnt2<-indent; cnt2++){
            strins(whitespace,str,cntbreak+1);
          }
        }
      } else {
        // problem breaking line..
        CONTINUE=FALSE;
      }
    }        
  } while (CONTINUE && cnt<strlen(str)-1); 
  return str;
}


void find_terminal_window_size(int *width, int *height){
  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  *width=w.ws_col;
  *height=w.ws_row;
}


double avg_harmonic(double arg1, double arg2){
  double ret;
  if (arg1==0.0) EXM_fatal_error("Problem: arg1 is zero in avg_harmonic().");
  if (arg2==0.0) EXM_fatal_error("Problem: arg2 is zero in avg_harmonic()."); 
  if (arg1<0.0) EXM_fatal_error("Problem: arg1 is negative in avg_harmonic().");
  if (arg2<0.0) EXM_fatal_error("Problem: arg2 is negative in avg_harmonic().");
  ret=2.0/(1.0/arg1+1.0/arg2);
  return(ret);
}




static int argrank(int argc, char **argv, char *arg){
  int cnt,tmp;
  bool FOUND;
  tmp=0;
  FOUND=FALSE;
  for (cnt=1; cnt<argc; cnt++){
   if (strcmp(argv[cnt],arg) == 0) {
     tmp=cnt;
     if (!FOUND) FOUND=TRUE; else EXM_fatal_error("The flag %s can not be called twice.",arg);
   }
  }
  return(tmp);
}



int process_flag_string(int argc, char **argv, char *flag, char **arg){
  int RET;
  int flagrank; 
  flagrank=argrank(argc,argv,flag);
  if (flagrank!=0) {
    if (argc<=flagrank+1) EXM_fatal_error("Missing argument after flag %s.",flag);
    *arg=(char *)realloc(*arg,(2+(long)strlen(argv[flagrank+1]))*sizeof(char));
    strcpy(*arg,argv[flagrank+1]);
    strcpy(argv[flagrank+0],"\0");
    strcpy(argv[flagrank+1],"\0");
    RET=2;
  } else {
    RET=0;
  }
  return(RET);
}


int process_flag_int(int argc, char **argv, char *flag, int *arg){
  int RET;
  int flagrank; 
  int eos=EOS;
  flagrank=argrank(argc,argv,flag);
  if (flagrank!=0) {
    if (argc<=flagrank+1) EXM_fatal_error("Missing argument after flag %s.",flag);
    if (sscanf(argv[flagrank+1], "%d%n", arg,&eos)!=1 || argv[flagrank+1][eos]!=EOS) {
      EXM_fatal_error("Expecting integer argument after %s flag.",flag);
    }
    strcpy(argv[flagrank+0],"\0");
    strcpy(argv[flagrank+1],"\0");
    RET=2;
  } else {
    RET=0;
  }
  return(RET);
}


int process_flag_long(int argc, char **argv, char *flag, long *arg){
  int RET;
  int flagrank; 
  int eos=EOS;
  flagrank=argrank(argc,argv,flag);
  if (flagrank!=0) {
    if (argc<=flagrank+1) EXM_fatal_error("Missing argument after flag %s.",flag);
    if (sscanf(argv[flagrank+1], "%ld%n", arg,&eos)!=1 || argv[flagrank+1][eos]!=EOS) {
      EXM_fatal_error("Expecting long integer argument after %s flag.",flag);
    }
    strcpy(argv[flagrank+0],"\0");
    strcpy(argv[flagrank+1],"\0");
    RET=2;
  } else {
    RET=0;
  }
  return(RET);
}


int process_flag_double(int argc, char **argv, char *flag, double *arg){
  int RET;
  int flagrank; 
  int eos=EOS;
  flagrank=argrank(argc,argv,flag);
  if (flagrank!=0) {
    if (argc<=flagrank+1) EXM_fatal_error("Missing argument after flag %s.",flag);
    if (sscanf(argv[flagrank+1], "%lf%n", arg,&eos)!=1 || argv[flagrank+1][eos]!=EOS) {
      EXM_fatal_error("Expecting integer argument after %s flag.",flag);
    }
    strcpy(argv[flagrank+0],"\0");
    strcpy(argv[flagrank+1],"\0");
    RET=2;
  } else {
    RET=0;
  }
  return(RET);
}


int process_flag_int_multiple(int argc, char **argv, char *flag, int **arg){
  int RET,cnt;
  int flagrank; 
  int eos=EOS;
  bool CONTINUE;
  flagrank=argrank(argc,argv,flag);
  if (flagrank!=0) {
    if (argc<=flagrank+1) EXM_fatal_error("Missing argument after flag %s.",flag);
    cnt=0;
    CONTINUE=TRUE;
    do {
      cnt++;
      *arg=(int *)realloc(*arg,(cnt+2)*sizeof(int));
      if ((flagrank+cnt)>=argc || sscanf(argv[flagrank+cnt], "%d%n", &((*arg)[cnt-1]),&eos)!=1 
        || argv[flagrank+cnt][eos]!=EOS) CONTINUE=FALSE; else CONTINUE=TRUE;
    } while(CONTINUE);
    RET=cnt;
    for (cnt=0; cnt<RET; cnt++) strcpy(argv[flagrank+cnt],"\0");
  } else {
    RET=0;
  }
  return(RET);
}


int process_flag_double_multiple(int argc, char **argv, char *flag, double **arg){
  int RET,cnt;
  int flagrank; 
  bool CONTINUE;
  flagrank=argrank(argc,argv,flag);
  if (flagrank!=0) {
    if (argc<=flagrank+1) EXM_fatal_error("Missing argument after flag %s.",flag);
    cnt=0;
    CONTINUE=TRUE;
    do {
      cnt++;
      *arg=(double *)realloc(*arg,(cnt+2)*sizeof(double));
      if ((flagrank+cnt)>=argc || sscanf(argv[flagrank+cnt], "%lg", &((*arg)[cnt-1]))!=1) 
        CONTINUE=FALSE; else CONTINUE=TRUE;
    } while(CONTINUE);
    RET=cnt;
    for (cnt=0; cnt<RET; cnt++) strcpy(argv[flagrank+cnt],"\0");
  } else {
    RET=0;
  }
  return(RET);
}


int process_flag(int argc, char **argv, char *flag){
  int RET;
  int flagrank; 
  flagrank=argrank(argc,argv,flag);
  if (flagrank!=0) {
    strcpy(argv[flagrank],"\0");
    RET=1;
  } else {
    RET=0;
  }
  return(RET);
}


int find_remaining_options(int argc, char **argv, char **options){
  long cnt;
  int RET;
  RET=0;
  *options=(char *)realloc(*options,sizeof(char));
  (*options)[0]=EOS;
  /* start cnt at 1 to exclude the executable */
  for (cnt=1; cnt<argc; cnt++){
    if ((argv[cnt])[0]!='\0'){
      RET++;
      *options=(char *)realloc(*options,(2+strlen(*options)+strlen(argv[cnt]))*sizeof(char));
      strcat(*options,argv[cnt]);
      strcat(*options," ");
    }
  }
  return(RET);
}




double min3(double val1, double val2, double val3){
  return(min(val1,min(val2,val3)));
}


double max3(double val1, double val2, double val3){
  return(max(val1,max(val2,val3)));
}


double notzero(double val, double sub){
  if (val==0.0) val=sub;
  return(val);
}


double minmod(double x, double y){
  double tmp;
  tmp=sign(x)*max(0.0e0,min(fabs(x),sign(x)*y));
  return(tmp);
}


double minmod_old(double val1, double val2){
  double ret;
  if (fabs(val1)>fabs(val2)) ret=val2; else ret=val1;
  if (val1*val2<0.0) ret=0.0;
  return(ret);
}


double minmod3(double x, double y, double z){
  double tmp;
  tmp=sign(x)*max(0.0e0,min(sign(x)*z,min(fabs(x),sign(x)*y)));
  return(tmp);
}


double maxmag(double x, double y){
  double ret;
  if (fabs(x)>fabs(y)) ret=x; else ret=y;
  return(ret);
}



#ifdef linux



#include <execinfo.h>

/* perform a backtrace on linux systems;
 * make sure to compile and link your code with debugging symbols (-g flag)
 * make sure to link your code with the -rdynamic flag
 * leave this function commented because it is not ANSI C compliant
*/
void output_backtrace(void) {
  void* callstack[128];
  int i, frames = backtrace(callstack, 128);
  char** strs = backtrace_symbols(callstack, frames);
  for (i = 0; i < frames; ++i) {
    printf("%s\n", strs[i]);
  }
  free(strs);
}

#else

void output_backtrace(void) {
  printf("No backtrace available. Use a linux OS to find backtrace.\n");
}

#endif

