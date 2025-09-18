// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000 Bernard Parent.
Copyright 2002 Derrick C. Alexander.

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
#include "2Dbase.h"
#include "2Dextra.h"

typedef struct {
  GRIDG_gl2d_t *gl2d;
  GRIDG_xygrid_t **xygrid;
} actionsarg2D_t;


static void verify_zone_validity(long is, long js, long ie, long je, GRIDG_gl2d_t gl2d, SOAP_codex_t *codex){
  if (ie<is) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be less than is (%ld).",ie,is);
  if (ie>gl2d.ie)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be greater than domain.ie (%ld) as specified in Size().",ie,gl2d.ie);
  if (is<gl2d.is)
    SOAP_fatal_error(codex,"Zone boundaries invalid: is (%ld) can not be less than domain.is (%ld) as specified in Size().",is,gl2d.is);
  if (je<js) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: je (%ld) can not be less than js (%ld).",je,js);
  if (je>gl2d.je)
    SOAP_fatal_error(codex,"Zone boundaries invalid: je (%ld) can not be greater than domain.je (%ld) as specified in Size().",je,gl2d.je);
  if (js<gl2d.js)
    SOAP_fatal_error(codex,"Zone boundaries invalid: js (%ld) can not be less than domain.js (%ld) as specified in Size().",js,gl2d.js);
}



static void verify_node_position_validity(long i, long j, GRIDG_gl2d_t gl2d, SOAP_codex_t *codex){
  if (i>gl2d.ie)
    SOAP_fatal_error(codex,"Node location invalid: i (%ld) can not be greater than domain.ie (%ld) as specified in Size().",i,gl2d.ie);
  if (i<gl2d.is)
    SOAP_fatal_error(codex,"Node location invalid: i (%ld) can not be less than domain.is (%ld) as specified in Size().",i,gl2d.is);
  if (j>gl2d.je)
    SOAP_fatal_error(codex,"Node location invalid: j (%ld) can not be greater than domain.je (%ld) as specified in Size().",j,gl2d.je);
  if (j<gl2d.js)
    SOAP_fatal_error(codex,"Node location invalid: j (%ld) can not be less than domain.js (%ld) as specified in Size().",j,gl2d.js);
}



void GRIDG_write_grid_2D_to_file(FILE **controlfile){
  fprintf(*controlfile," \n");
  fprintf(*controlfile,
    "\n"
    "is=1;\n"
    "ie=30;\n"
    "js=1;\n"
    "je=30;\n"
    "\n"
    "Grid(\n"
    "  Size(is,js, ie,je);\n"
    "  Point(is,js, 0.0e0,0.0e0);\n"
    "  Point(ie,js, 1.0e0,0.0e0);\n"
    "  Point(ie,je, 1.0e0,1.0e0);\n"
    "  Point(is,je, 0.0e0,1.0e0);\n"
    "  JoinCorners(is,js, ie,je,  EE,0.5e0,1.0e0,1.0e0,  EE,0.5e0,1.0e0,1.0e0);\n"
    ");\n"
  );
}


static double find_distance(GRIDG_xygrid_t *xygrid, GRIDG_gl2d_t gl2d,
                    long i1, long j1,
                    long i2, long j2){
  double dist;
  dist=sqrt(sqr(xygrid[GRIDG_ai2(gl2d,i2,j2)].x-xygrid[GRIDG_ai2(gl2d,i1,j1)].x)
           +sqr(xygrid[GRIDG_ai2(gl2d,i2,j2)].y-xygrid[GRIDG_ai2(gl2d,i1,j1)].y)
           );
  return(dist);
}

void Size2D(GRIDG_gl2d_t *gl, GRIDG_xygrid_t **xygrid,
            long is, long js, long ie, long je){
  long i,j;
  gl->is=is;  gl->js=js;
  gl->ie=ie;  gl->je=je;

  *xygrid = (GRIDG_xygrid_t *) malloc((gl->ie-gl->is+1+4)*(gl->je-gl->js+1+4)
                                     *sizeof(GRIDG_xygrid_t));
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      (*xygrid)[GRIDG_ai2(*gl,i,j)].x=(double)(i-is)/(double)(ie-is);
      (*xygrid)[GRIDG_ai2(*gl,i,j)].y=(double)(j-js)/(double)(je-js);
      (*xygrid)[GRIDG_ai2(*gl,i,j)].INIT=FALSE;
    }
  }
}


static void Size_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t *gl2d, GRIDG_xygrid_t **xygrid){
  long is,js,ie,je;
  int eos = EOS;

  SOAP_substitute_all_argums(argum, codex);

  if (sscanf(*argum,"%ld,%ld,%ld,%ld%n",&is,&js,&ie,&je,&eos)!=4 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Size().");
  }
  if (is>=ie) SOAP_fatal_error(codex,"In Size(), is (%ld) must be less than ie (%ld).",is,ie);
  if (js>=je) SOAP_fatal_error(codex,"In Size(), js (%ld) must be less than je (%ld).",js,je);
  Size2D(gl2d, xygrid, is, js, ie, je);
}


void Point2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid,
             long i, long j, double x, double y){
  xygrid[GRIDG_ai2(gl2d,i,j)].x=x;
  xygrid[GRIDG_ai2(gl2d,i,j)].y=y;
  xygrid[GRIDG_ai2(gl2d,i,j)].INIT=TRUE;
}


static void Point_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long i,j;
  double x,y;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%lg,%lg%n",&i,&j,&x,&y,&eos)!=4 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Point().");
  }
  verify_node_position_validity(i, j, gl2d, codex);
  Point2D(gl2d, xygrid, i, j, x, y);
}


void Corners2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid,
                   long i1, long j1, long i2, long j2,
                   double x1, double y1, double x2, double y2){
  xygrid[GRIDG_ai2(gl2d,i1,j1)].x=x1;  xygrid[GRIDG_ai2(gl2d,i1,j1)].y=y1;  xygrid[GRIDG_ai2(gl2d,i1,j1)].INIT=TRUE;
  xygrid[GRIDG_ai2(gl2d,i2,j1)].x=x2;  xygrid[GRIDG_ai2(gl2d,i2,j1)].y=y1;  xygrid[GRIDG_ai2(gl2d,i2,j1)].INIT=TRUE;
  xygrid[GRIDG_ai2(gl2d,i2,j2)].x=x2;  xygrid[GRIDG_ai2(gl2d,i2,j2)].y=y2;  xygrid[GRIDG_ai2(gl2d,i2,j2)].INIT=TRUE;
  xygrid[GRIDG_ai2(gl2d,i1,j2)].x=x1;  xygrid[GRIDG_ai2(gl2d,i1,j2)].y=y2;  xygrid[GRIDG_ai2(gl2d,i1,j2)].INIT=TRUE;
}





static void Corners_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long i1,j1,i2,j2;
  double x1,y1,x2,y2;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%lg,%lg,%lg,%lg%n",
         &i1,&j1,&i2,&j2,&x1,&y1,&x2,&y2,&eos)!=8 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Corners().");
  }
  verify_zone_validity(i1, j1, i2, j2, gl2d, codex);
  Corners2D(gl2d, xygrid,i1,j1,i2,j2,x1,y1,x2,y2);
}

static int find_xy_along_ij_local(GRIDG_xygrid_t *xygrid, double segval1, double segval2,
                     double NLoverN, char modeL, char modeH,
                     GRIDG_gl2d_t gl2d, SOAP_codex_t *codex, long i1, long j1, long i2, long j2,
                     int(*find_segment)(double *, long , double,
                           double , double , double, SOAP_codex_t * )){
  int err;
  if (!xygrid[GRIDG_ai2(gl2d,i1,j1)].INIT){
    SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Grid segment could not be formed.",i1,j1);
  }
  if (!xygrid[GRIDG_ai2(gl2d,i2,j2)].INIT){
    SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Grid segment could not be formed.",i2,j2);
  }
  if (modeL=='G' || modeL=='g') {
    if (!xygrid[GRIDG_ai2(gl2d,i1+sgndif(i1,i2),j1+sgndif(j1,j2))].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Grid segment could not be formed.",i1+sgndif(i1,i2),j1+sgndif(j1,j2));
    }
    segval1=find_distance(xygrid, gl2d, i1, j1, i1+sgndif(i1,i2), j1+sgndif(j1,j2));
  }
  if (modeH=='G' || modeH=='g') {
    if (!xygrid[GRIDG_ai2(gl2d,i2-sgndif(i1,i2),j2-sgndif(j1,j2))].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Grid segment could not be formed.",i2-sgndif(i1,i2),j2-sgndif(j1,j2));
    }
    segval2=find_distance(xygrid, gl2d, i2, j2, i2-sgndif(i1,i2), j2-sgndif(j1,j2));
  }
  err=find_xy_along_ij(xygrid, segval1, segval2, NLoverN,
                    gl2d, i1, j1, i2, j2, codex, find_segment);
  return(err);
}


int JoinCorners2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
             long i1, long j1, long i2, long j2,
             char modeL_i, char modeH_i, double NLoverN_i, double valL_i, double valH_i,
             char modeL_j, char modeH_j, double NLoverN_j, double valL_j, double valH_j){
  long i;
  int err;
  int (*FindSegment_i)(double *, long , double, double , double , double, SOAP_codex_t * );
  int (*FindSegment_j)(double *, long , double, double , double , double, SOAP_codex_t * );

  find_segment_type(&FindSegment_i, modeL_i, modeH_i, codex);
  find_segment_type(&FindSegment_j, modeL_j, modeH_j, codex);

  err=find_xy_along_ij_local(xygrid, valL_i, valH_i, NLoverN_i, modeL_i, modeH_i, gl2d, codex,
                  i1, j1, i2, j1, FindSegment_i);
  err+=find_xy_along_ij_local(xygrid, valL_i, valH_i, NLoverN_i, modeL_i, modeH_i, gl2d, codex,
                  i1, j2, i2, j2, FindSegment_i);
  for (i=i1; i<=i2 && err==0; i++){
    err+=find_xy_along_ij_local(xygrid, valL_j, valH_j, NLoverN_j, modeL_j, modeH_j, gl2d, codex,
                  i, j1, i, j2, FindSegment_j);
  }
  return(err);
}


static void JoinCorners_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long i1,j1,i2,j2,cnt;
  int err;
  char modeL_i,modeH_i,modeL_j,modeH_j;
  double NLoverN_i,valL_i,valH_i,NLoverN_j,valL_j,valH_j;
  int eos=EOS;
  modeL_i='0';
  modeL_j='0';
  modeH_i='0';
  modeH_j='0';
  for (cnt=0; cnt<SOAP_number_argums(*argum); cnt++)
    if (cnt!=4 && cnt!=8) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,"
               "%c%c,%lg,%lg,%lg,"
               "%c%c,%lg,%lg,%lg%n",
         &i1,&j1,&i2,&j2,
         &modeL_i,&modeH_i,&NLoverN_i,&valL_i,&valH_i,
         &modeL_j,&modeH_j,&NLoverN_j,&valL_j,&valH_j,&eos)!=14 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in JoinCorners().");
  }
  if ((i1==i2) && !(modeL_i=='N' && modeH_i=='O')) {
     modeL_i='N';
     modeH_i='O';
  }
  if ((j1==j2) && !(modeL_j=='N' && modeH_j=='O')) {
     modeL_j='N';
     modeH_j='O';
  }
  verify_zone_validity(i1, j1, i2, j2, gl2d, codex);

  err=JoinCorners2D(gl2d, xygrid, codex, i1, j1, i2, j2,
              modeL_i, modeH_i, NLoverN_i, valL_i, valH_i,
              modeL_j, modeH_j, NLoverN_j, valL_j, valH_j);
  if (err!=0) SOAP_fatal_error(codex,"Error occurred while executing JoinCorners() command.");
}




/* Algorithm by Derrick Alexander */
void JoinFaces2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
               long i1, long j1, long i2, long j2){
  long i,j;
  double fy_start,fy_end,fx_top,fx_bot,y1,y2;
  bool VALID;

  VALID=FALSE;
  for (i=i1+1; i<i2; i++){
    for (j=j1+1; j<j2; j++){
      VALID=TRUE;
      if (!xygrid[GRIDG_ai2(gl2d,i1,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. JoinFaces command failed.",i1,j);
      }
      if (!xygrid[GRIDG_ai2(gl2d,i2,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. JoinFaces command failed.",i2,j);
      }
      if (!xygrid[GRIDG_ai2(gl2d,i,j1)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. JoinFaces command failed.",i,j1);
      }
      if (!xygrid[GRIDG_ai2(gl2d,i,j2)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. JoinFaces command failed.",i,j2);
      }

      fy_start=(xygrid[GRIDG_ai2(gl2d,i1,j)].y-xygrid[GRIDG_ai2(gl2d,i1,j1)].y)
              /(xygrid[GRIDG_ai2(gl2d,i1,j2)].y-xygrid[GRIDG_ai2(gl2d,i1,j1)].y);
      fy_end=  (xygrid[GRIDG_ai2(gl2d,i2,j)].y-xygrid[GRIDG_ai2(gl2d,i2,j1)].y)
              /(xygrid[GRIDG_ai2(gl2d,i2,j2)].y-xygrid[GRIDG_ai2(gl2d,i2,j1)].y);
      fx_top=  (xygrid[GRIDG_ai2(gl2d,i,j2)].x-xygrid[GRIDG_ai2(gl2d,i1,j2)].x)
              /(xygrid[GRIDG_ai2(gl2d,i2,j2)].x-xygrid[GRIDG_ai2(gl2d,i1,j2)].x);
      fx_bot=  (xygrid[GRIDG_ai2(gl2d,i,j1)].x-xygrid[GRIDG_ai2(gl2d,i1,j1)].x)
              /(xygrid[GRIDG_ai2(gl2d,i2,j1)].x-xygrid[GRIDG_ai2(gl2d,i1,j1)].x);

      y1=-(fy_end-fy_start)*(fx_top-fx_bot)*xygrid[GRIDG_ai2(gl2d,i,j1)].y
         +xygrid[GRIDG_ai2(gl2d,i,j1)].y;
      y2=((fy_end-fy_start)*fx_bot+fy_start)*
         (xygrid[GRIDG_ai2(gl2d,i,j2)].y-xygrid[GRIDG_ai2(gl2d,i,j1)].y);
      xygrid[GRIDG_ai2(gl2d,i,j)].y=(y1+y2)/(1.0e0-(fy_end-fy_start)*(fx_top-fx_bot));
      xygrid[GRIDG_ai2(gl2d,i,j)].x=((fx_top-fx_bot)
                              *(xygrid[GRIDG_ai2(gl2d,i,j)].y-xygrid[GRIDG_ai2(gl2d,i,j1)].y)
                              /(xygrid[GRIDG_ai2(gl2d,i,j2)].y-xygrid[GRIDG_ai2(gl2d,i,j1)].y)+fx_bot)
                              *(xygrid[GRIDG_ai2(gl2d,i2,j)].x-xygrid[GRIDG_ai2(gl2d,i1,j)].x)
                             +xygrid[GRIDG_ai2(gl2d,i1,j)].x;
      xygrid[GRIDG_ai2(gl2d,i,j)].INIT=TRUE;
    }
  }

  if (!VALID){
    SOAP_fatal_error(codex,"When using the JoinFaces() command, set i2>i1+1 and j2>j1+1.");
  }
}


int Join2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long i1, long j1, long i2, long j2,
           char index, char modeL, char modeH,
           double NLoverN, double valL, double valH){
  long i,j;
  int err;
  int (*FindSegment)(double *, long , double, double , double , double, SOAP_codex_t * );

  find_segment_type(&FindSegment, modeL, modeH, codex);
  err=0;
  switch (index) {
    case 'i':
      for (j=j1; j<=j2 && err==0; j++){
        err+=find_xy_along_ij_local(xygrid, valL, valH,
                          NLoverN, modeL, modeH, gl2d, codex,
                          i1, j, i2, j, FindSegment);
      }
    break;
    case 'j':
      for (i=i1; i<=i2 && err==0; i++){
        err+=find_xy_along_ij_local(xygrid, valL, valH,
                          NLoverN, modeL, modeH, gl2d, codex,
                          i, j1, i, j2, FindSegment);
      }
    break;
  }
  return(err);
}


static void Join_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long i1,j1,i2,j2,cnt;
  char index,modeL,modeH;
  double NLoverN,valL,valH;
  int err;
  int eos=EOS;

  for (cnt=0; cnt<9; cnt++) if (cnt!=4 && cnt!=5) {
    SOAP_substitute_argum(argum, cnt, codex);
  }
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%c,%c%c,%lg,%lg,%lg%n",
         &i1,&j1,&i2,&j2,&index,&modeL,&modeH,&NLoverN,
         &valL,&valH,&eos)!=10 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Join().");
  }
  verify_zone_validity(i1, j1, i2, j2, gl2d, codex);
  if (index!='i' && index!='j') SOAP_fatal_error(codex,"In Join(), the index name must be either 'i' or 'j'.");
  err=Join2D(gl2d, xygrid, codex, i1, j1, i2, j2, index, modeL, modeH, NLoverN, valL, valH);
  if (err!=0) SOAP_fatal_error(codex,"Error occurred while executing Join() command.");
}


static void JoinFaces_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long i1,j1,i2,j2;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld%n",&i1,&j1,&i2,&j2,&eos)!=4 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in JoinFaces().");
  }
  verify_zone_validity(i1, j1, i2, j2, gl2d, codex);
  JoinFaces2D(gl2d, xygrid, codex, i1, j1, i2, j2);
}


void Translate2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           double dx, double dy){
  long i,j;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Translate command failed.",i,j);
      }
      xygrid[GRIDG_ai2(gl2d,i,j)].x=xygrid[GRIDG_ai2(gl2d,i,j)].x+dx;
      xygrid[GRIDG_ai2(gl2d,i,j)].y=xygrid[GRIDG_ai2(gl2d,i,j)].y+dy;
    }
  }
}


static void translate_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je;
  double dx,dy;
  int eos=EOS;
  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%lg,%lg%n",
               &is,&js,&ie,&je,&dx,&dy,&eos)!=6 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Translate().");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);
  Translate2D(gl2d, xygrid, codex, is, js, ie, je, dx, dy);
}


void Scale2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           double x_orig, double y_orig, double scalefactx, double scalefacty){
  long i,j;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Scale command failed.",i,j);
      }
      xygrid[GRIDG_ai2(gl2d,i,j)].x=x_orig+(xygrid[GRIDG_ai2(gl2d,i,j)].x-x_orig)*scalefactx;
      xygrid[GRIDG_ai2(gl2d,i,j)].y=y_orig+(xygrid[GRIDG_ai2(gl2d,i,j)].y-y_orig)*scalefacty;
    }
  }
}

static void Scale_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je;
  double x_orig,y_orig,scalefactx,scalefacty;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%lg,%lg,%lg,%lg%n",
               &is,&js,&ie,&je,&x_orig,&y_orig,&scalefactx,&scalefacty,&eos)!=8 
    || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Scale().");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);
  if (scalefactx<=0.0 || scalefacty<=0.0) SOAP_fatal_error(codex,"In Scale(), the scale factors must all be positive.");
  Scale2D(gl2d, xygrid, codex, is, js, ie, je, x_orig, y_orig, scalefactx, scalefacty);
}


void ReverseIndex2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je, char index){
  long i,j;
  GRIDG_xygrid_t tmp;

  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (index=='i' && i<ie-(i-is)) {
        if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. ReverseIndex command failed.",i,j);
        }
        tmp=xygrid[GRIDG_ai2(gl2d,i,j)];
        xygrid[GRIDG_ai2(gl2d,i,j)]=xygrid[GRIDG_ai2(gl2d,ie-(i-is),j)];
        xygrid[GRIDG_ai2(gl2d,ie-(i-is),j)]=tmp;
      }
      if (index=='j' && j<je-(j-js)) {
        if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. ReverseIndex command failed.",i,j);
        }
        tmp=xygrid[GRIDG_ai2(gl2d,i,j)];
        xygrid[GRIDG_ai2(gl2d,i,j)]=xygrid[GRIDG_ai2(gl2d,i,je-(j-js))];
        xygrid[GRIDG_ai2(gl2d,i,je-(j-js))]=tmp;
      }
    }
  }
}


static void ReverseIndex_argum(char **argum, SOAP_codex_t *codex,
                         GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je,cnt;
  char index;
  int eos=EOS;

  for (cnt=0; cnt<4; cnt++) SOAP_substitute_argum(argum, cnt, codex);

  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%c%n",
         &is,&js,&ie,&je,&index,&eos)!=5 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in ReverseIndex().");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);
  if (index!='i' && index!='j') SOAP_fatal_error(codex,"In ReverseIndex(), the index name must be either 'i' or 'j'.");
  ReverseIndex2D(gl2d, xygrid, codex, is, js, ie, je, index);
}


void Rotate2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           double pivot_x, double pivot_y,
           double angle){
  long i,j;
  double r,phi,dx,dy;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Rotate command failed.",i,j);
      }
      dx=xygrid[GRIDG_ai2(gl2d,i,j)].x-pivot_x;
      dy=xygrid[GRIDG_ai2(gl2d,i,j)].y-pivot_y;
      r=sqrt(sqr(dx)+sqr(dy));
      if (dx==0.0e0) phi=0.5e0*pi*sign(dy); 
        else phi=atan(dy/dx);
      if (dx<0.0e0) phi=phi+pi;
      xygrid[GRIDG_ai2(gl2d,i,j)].x=pivot_x+r*cos(phi+angle);
      xygrid[GRIDG_ai2(gl2d,i,j)].y=pivot_y+r*sin(phi+angle);
    }
  }
}


static void Rotate_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je;
  double angle,pivot_x,pivot_y;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%lg,%lg,%lg%n",
         &is,&js,&ie,&je,&pivot_x,&pivot_y,&angle,&eos)!=7 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Rotate().");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);
  Rotate2D(gl2d, xygrid, codex, is, js, ie, je, pivot_x, pivot_y, angle);
}


void Mirror2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           char axis_name, double axis_pos){
  long i,j;

  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Mirror command failed.",i,j);
      }

      if (axis_name=='x') {
        xygrid[GRIDG_ai2(gl2d,i,j)].x=2.0e0*axis_pos-xygrid[GRIDG_ai2(gl2d,i,j)].x;
      }
      if (axis_name=='y') {
        xygrid[GRIDG_ai2(gl2d,i,j)].y=2.0e0*axis_pos-xygrid[GRIDG_ai2(gl2d,i,j)].y;
      }
    }
  }
}


static void Mirror_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je;
  char axis_name;
  double axis_pos;
  long cnt;
  int eos=EOS;

  for (cnt=0; cnt<6; cnt++) if (cnt!=4) SOAP_substitute_argum(argum, cnt, codex);

  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%c,%lg%n",
         &is,&js,&ie,&je,&axis_name,&axis_pos,&eos)!=6 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Mirror().");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);
  if (axis_name!='x' && axis_name!='y') SOAP_fatal_error(codex,"In Mirror(), the axis name must be either 'x' or 'y'.");
  Mirror2D(gl2d, xygrid, codex, is, js, ie, je, axis_name, axis_pos);
}



typedef struct {
  double val1;
  char *expr;
  char varchange;
  SOAP_codex_t codex;
} argequation2D_t;


static double Equation2D(void *argequation2D, double val2){
  double res,y,x;
  char ystr[40],xstr[40];
  char varchange;
  char *expr;
  SOAP_codex_t *codex;
  int eos=EOS;

  codex=&(((argequation2D_t *)argequation2D)->codex);
  varchange=((argequation2D_t *)argequation2D)->varchange;
  x=0.0e0; y=0.0e0;
  if (varchange=='x') {
    y=((argequation2D_t *)argequation2D)->val1;
    x=val2;
  }
  if (varchange=='y') {
    x=((argequation2D_t *)argequation2D)->val1;
    y=val2;
  }
  expr=(char *)malloc((3+(long)strlen(((argequation2D_t *)argequation2D)->expr))*sizeof(char));
  strcpy(expr,((argequation2D_t *)argequation2D)->expr);
  /* here, setup codex and solve expression */
  sprintf(ystr,"%15.15E",y);
  sprintf(xstr,"%15.15E",x);
  SOAP_add_to_vars(codex, "x", xstr);
  SOAP_add_to_vars(codex, "y", ystr);
 // printf("res=%E x=%E  y=%E  xstr=%s ystr=%s\n",res,x,y,xstr,ystr);
  SOAP_substitute_expression(&expr, codex);
  /* done */
  if (sscanf(expr,"%lg%n",&res,&eos)!=1 || expr[eos]!=EOS){
    GRIDG_fatal_error("Problem evaluating expression within Equation() command.");
  }
  
  free(expr);
  return(res);
}


static void Equation_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je,i,j;
  double drootinit;
  char varchange;
  argequation2D_t argequation2D;
  char *expr;
  long IFLAG,cnt;

  for (cnt=0; cnt<4; cnt++) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%c,",
         &is,&js,&ie,&je,&varchange)!=5){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Equation();");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);
  if (varchange!='x' && varchange!='y') SOAP_fatal_error(codex,"In Equation(), the variable to be changed must be either 'x' or 'y'.");

  expr=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&expr,*argum,5);
  replace_equal_sign(&expr);
  argequation2D.expr=(char *)malloc(((long)strlen(expr)+30)*sizeof(char));
  strcpy(argequation2D.expr,expr);
  SOAP_copy_codex(codex,&(argequation2D.codex));
  argequation2D.codex.vars=codex->vars;

  drootinit=0.0;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      drootinit=max(drootinit,fabs(xygrid[GRIDG_ai2(gl2d,i,j)].y));
      drootinit=max(drootinit,fabs(xygrid[GRIDG_ai2(gl2d,i,j)].x));
    }
  }
  drootinit/=1e7;
  drootinit=max(1.0e-20,drootinit);
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Equation command failed.",i,j);
      }

      argequation2D.varchange=varchange;
      if (varchange=='x') {
        argequation2D.val1=xygrid[GRIDG_ai2(gl2d,i,j)].y;
        xygrid[GRIDG_ai2(gl2d,i,j)].x=
          EXM_find_root_Newton_Raphson(&Equation2D, &argequation2D, xygrid[GRIDG_ai2(gl2d,i,j)].x,
                         drootinit, 1.0e-10, 1.0e-20, &IFLAG);
      }
      if (varchange=='y') {
        argequation2D.val1=xygrid[GRIDG_ai2(gl2d,i,j)].x;
        xygrid[GRIDG_ai2(gl2d,i,j)].y=
          EXM_find_root_Newton_Raphson(&Equation2D, &argequation2D, xygrid[GRIDG_ai2(gl2d,i,j)].y,
                         drootinit, 1.0e-10, 1.0e-20, &IFLAG);
      }
      if (IFLAG!=1) SOAP_fatal_error(codex,"In Equation(), root could not be obtained with Newton-Raphson solver for node %ld,%ld.",i,j);
    }
  }

  free(expr);
  codex->vars=argequation2D.codex.vars;
  free( argequation2D.expr );
  SOAP_free_codex_copy(&(argequation2D.codex));
}


static void Spline_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je,i,j;
  long cnt;
  long N,n;
  double *f,*b,*x;
  double thisx;
  char *expr;

  for (cnt=0; cnt<4; cnt++) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,",
         &is,&js,&ie,&je)!=4){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Spline();");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);

  expr=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&expr,*argum,4);

  if (strcmp(expr,"y(x)")!=0 && strcmp(expr,"x(y)")!=0) SOAP_fatal_error(codex,"In Spline(), the 5th argument must be either y(x) or x(y).");

  N=SOAP_number_argums(*argum)-5;
  if (mod(N,2)!=0) SOAP_fatal_error(codex,"Number of arguments within Spline must be 5 + some even number.");
  N=N/2;
  if (N<4) SOAP_fatal_error(codex,"Number of data points supplied within Spline must be at least 4.");

  x=(double *)malloc(N*sizeof(double));
  f=(double *)malloc(N*sizeof(double));
  b=(double *)malloc(N*sizeof(double));
  
  for (n=0; n<N; n++) {
    SOAP_substitute_argum(argum, 5+n*2, codex);    
    x[n]=SOAP_get_argum_double(codex, *argum, 5+n*2);
    SOAP_substitute_argum(argum, 5+n*2+1, codex);    
    f[n]=SOAP_get_argum_double(codex, *argum, 5+n*2+1);
  }
  /* check if data points are valid (x[n+1]>x[n]) */
  for (n=0; n<N-1; n++){
    if (x[n+1]<=x[n]) SOAP_fatal_error(codex, "Data points supplied to spline must be such that x[i+1]>x[i] if y(x) or such that y[i+1]>y[i] if x(y)."); 
  }

  EXM_find_spline(N, x, f, b);
  
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT)
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Spline command failed.",i,j);
      switch (expr[2]){
        case 'x':
          thisx=xygrid[GRIDG_ai2(gl2d,i,j)].x;
        break;
        case 'y':
          thisx=xygrid[GRIDG_ai2(gl2d,i,j)].y;
        break;
        default:
          thisx=0.0; //to avoid compiler warning
          GRIDG_fatal_error("expr[2] can not be set to %c in Spline_argum().",expr[2]);
      }
      if (thisx>=x[0] && thisx<=x[N-1]) { 
        switch (expr[0]){
          case 'x':
            xygrid[GRIDG_ai2(gl2d,i,j)].x=EXM_f_from_spline(N, x, f, b, thisx);
          break;
          case 'y':
            xygrid[GRIDG_ai2(gl2d,i,j)].y=EXM_f_from_spline(N, x, f, b, thisx);
          break;
          default:
            GRIDG_fatal_error("expr[0] can not be set to %c in Spline_argum().",expr[0]);
        }
      } else {
        SOAP_fatal_error(codex,"Point %E is outside of the supplied range %E to %E.",thisx,x[0],x[N-1]);
      }
    }
  }


  free(expr);
  free(x);
  free(b);
  free(f);
}


void Copy2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           long id, long jd){
  long i,j;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
        SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Copy command failed.",i,j);
      }

      xygrid[GRIDG_ai2(gl2d,id-is+i,jd-js+j)].x=xygrid[GRIDG_ai2(gl2d,i,j)].x;
      xygrid[GRIDG_ai2(gl2d,id-is+i,jd-js+j)].y=xygrid[GRIDG_ai2(gl2d,i,j)].y;
      xygrid[GRIDG_ai2(gl2d,id-is+i,jd-js+j)].INIT=TRUE;
    }
  }
}


static void Copy_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long is,js,ie,je,id,jd;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld%n",&is,&js,&ie,&je,&id,&jd,&eos)!=6 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Copy().");
  }
  verify_zone_validity(is, js, ie, je, gl2d, codex);
  verify_node_position_validity(id, jd,  gl2d, codex);
  Copy2D(gl2d, xygrid, codex, is, js, ie, je, id, jd);
}


static void actions(char *actionname, char **argum, SOAP_codex_t *codex){
  GRIDG_gl2d_t *gl2d;
  GRIDG_xygrid_t **xygrid;

  gl2d=((actionsarg2D_t *)codex->action_args)->gl2d;
  xygrid=((actionsarg2D_t *)codex->action_args)->xygrid;


  if (strcmp(actionname,"Size")==0) {
    Size_argum(argum,codex,gl2d,xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Point")==0) {
    Point_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Corners")==0) {
    Corners_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"JoinCorners")==0) {
    JoinCorners_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Mirror")==0) {
    Mirror_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"ReverseIndex")==0) {
    ReverseIndex_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Join")==0) {
    Join_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"JoinFaces")==0) {
    JoinFaces_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Translate")==0) {
    translate_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Scale")==0) {
    Scale_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Rotate")==0) {
    Rotate_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Equation")==0) {
    Equation_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Spline")==0) {
    Spline_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Copy")==0) {
    Copy_argum(argum,codex,*gl2d,*xygrid);
    codex->ACTIONPROCESSED=TRUE;
  }
}

static void functions(char *functionname, char **argum,
                      char **returnstr, SOAP_codex_t *codex){
  GRIDG_gl2d_t *gl2d;
  GRIDG_xygrid_t **xygrid;
  long i,j;
  int eos=EOS;

  gl2d=((actionsarg2D_t *)codex->action_args)->gl2d;
  xygrid=((actionsarg2D_t *)codex->action_args)->xygrid;

  if (strcmp(functionname,"_x")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%ld%n",&i,&j,&eos)!=2 || (*argum)[eos]!=EOS)
       SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in _x().");
    verify_node_position_validity(i, j, *gl2d, codex);
    if (!(*xygrid)[GRIDG_ai2(*gl2d,i,j)].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized, _x command failed.",i,j);
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",(*xygrid)[GRIDG_ai2(*gl2d,i,j)].x);
  }
  if (strcmp(functionname,"_y")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%ld%n",&i,&j,&eos)!=2 || (*argum)[eos]!=EOS)
       SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in _y().");
    verify_node_position_validity(i, j, *gl2d, codex);
    if (!(*xygrid)[GRIDG_ai2(*gl2d,i,j)].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized, _y command failed.",i,j);
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",(*xygrid)[GRIDG_ai2(*gl2d,i,j)].y);
  }
}


static long add_one(long i, long imax){
  if (i==imax) i--; else i++;
  return(i);
}

static double _x(GRIDG_xygrid_t grid, long dim){
  if (dim==0) return(grid.x);
  return(grid.y);
}

/* this function finds the x,y coordinates of "ghost nodes" located outside the grid
   limits by "mirroring" the inner within the grid with respect to the boundary plane */ 
void find_xy_on_ghost_nodes_mirror(GRIDG_xygrid_t *grid, GRIDG_gl2d_t gl){
  long dim,i,j;
  EXM_vec3D_t pa,pb,pc,pp_o,pp_m;

  pa[0]=1.0e0;   pa[1]=0.0e0;   pa[2]=0.0e0;
  pb[0]=0.0e0;   pb[1]=1.0e0;   pb[2]=0.0e0;
  pc[0]=0.0e0;   pc[1]=0.0e0;   pc[2]=1.0e0;
  pp_o[0]=0.0e0; pp_o[1]=0.0e0; pp_o[2]=0.0e0;

  /* here the idea is to find a plane in 3D represented
     by 3 points which lie on the boundary surface:
     pa,pb and pc. */

  /* do j=ajs a3 j=aje */
  for(i=gl.is; i<=gl.ie; i++){
    /* do j=aje */
    for (dim=0; dim<2; dim++){
      pa[dim]=_x(grid[GRIDG_ai2(gl,add_one(i,gl.ie),gl.je)],dim);
      pb[dim]=_x(grid[GRIDG_ai2(gl,i,gl.je)],dim);
      pc[dim]=_x(grid[GRIDG_ai2(gl,i,gl.je)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,i,gl.je-1)],dim);
    }
    EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
    grid[GRIDG_ai2(gl,i,gl.je+1)].x=pp_m[0];
    grid[GRIDG_ai2(gl,i,gl.je+1)].y=pp_m[1];
    /* do j=ajs */
    for (dim=0; dim<2; dim++){
      pa[dim]=_x(grid[GRIDG_ai2(gl,add_one(i,gl.ie),gl.js)],dim);
      pb[dim]=_x(grid[GRIDG_ai2(gl,i,gl.js)],dim);
      pc[dim]=_x(grid[GRIDG_ai2(gl,i,gl.js)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,i,gl.js+1)],dim);
    }
    EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
    grid[GRIDG_ai2(gl,i,gl.js-1)].x=pp_m[0];
    grid[GRIDG_ai2(gl,i,gl.js-1)].y=pp_m[1];
  }

  /* do i=GRIDG_ai3s and i=GRIDG_ai3e */
  /* copy x,y,z of is plane to is-1 plane and of ie plane to ie+1 plane*/
  for(j=gl.js-1; j<=gl.je+1; j++){
    /* do i=GRIDG_ai3e */
    for (dim=0; dim<2; dim++){
      pa[dim]=_x(grid[GRIDG_ai2(gl,gl.ie,j)],dim);
      pb[dim]=_x(grid[GRIDG_ai2(gl,gl.ie,add_one(j,gl.je+1))],dim);
      pc[dim]=_x(grid[GRIDG_ai2(gl,gl.ie,j)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,gl.ie-1,j)],dim);
    }
    EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
    grid[GRIDG_ai2(gl,gl.ie+1,j)].x=pp_m[0];
    grid[GRIDG_ai2(gl,gl.ie+1,j)].y=pp_m[1];
    /* do i=GRIDG_ai3s */
    for (dim=0; dim<2; dim++){
      pa[dim]=_x(grid[GRIDG_ai2(gl,gl.is,j)],dim);
      pb[dim]=_x(grid[GRIDG_ai2(gl,gl.is,add_one(j,gl.je+1))],dim);
      pc[dim]=_x(grid[GRIDG_ai2(gl,gl.is,j)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,gl.is+1,j)],dim);
    }
    EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
    grid[GRIDG_ai2(gl,gl.is-1,j)].x=pp_m[0];
    grid[GRIDG_ai2(gl,gl.is-1,j)].y=pp_m[1];
  }
}



/* this function finds the x,y coordinates of "ghost nodes" located outside the grid
   limits by "mirroring" the inner within the grid with respect to the boundary plane */ 
void find_xy_on_ghost_nodes_extrapolate(GRIDG_xygrid_t *grid, GRIDG_gl2d_t gl){
  long dim,i,j;
  EXM_vec3D_t pc,pp_o;

  pc[0]=0.0e0;   pc[1]=0.0e0;   pc[2]=1.0e0;
  pp_o[0]=0.0e0; pp_o[1]=0.0e0; pp_o[2]=0.0e0;

  /* here the idea is to find a plane in 3D represented
     by 3 points which lie on the boundary surface:
     pa,pb and pc. */

  /* do j=ajs a3 j=aje */
  for(i=gl.is; i<=gl.ie; i++){
    /* do j=aje */
    for (dim=0; dim<2; dim++){
      pc[dim]=_x(grid[GRIDG_ai2(gl,i,gl.je)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,i,gl.je-1)],dim);
    }
    grid[GRIDG_ai2(gl,i,gl.je+1)].x=2.0*pc[0]-pp_o[0];
    grid[GRIDG_ai2(gl,i,gl.je+1)].y=2.0*pc[1]-pp_o[1];
    /* do j=ajs */
    for (dim=0; dim<2; dim++){
      pc[dim]=_x(grid[GRIDG_ai2(gl,i,gl.js)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,i,gl.js+1)],dim);
    }
    grid[GRIDG_ai2(gl,i,gl.js-1)].x=2.0*pc[0]-pp_o[0];
    grid[GRIDG_ai2(gl,i,gl.js-1)].y=2.0*pc[1]-pp_o[1];
  }

  /* do i=GRIDG_ai3s and i=GRIDG_ai3e */
  /* copy x,y,z of is plane to is-1 plane and of ie plane to ie+1 plane*/
  for(j=gl.js-1; j<=gl.je+1; j++){
    /* do i=GRIDG_ai3e */
    for (dim=0; dim<2; dim++){
      pc[dim]=_x(grid[GRIDG_ai2(gl,gl.ie,j)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,gl.ie-1,j)],dim);
    }
    grid[GRIDG_ai2(gl,gl.ie+1,j)].x=2.0*pc[0]-pp_o[0];
    grid[GRIDG_ai2(gl,gl.ie+1,j)].y=2.0*pc[1]-pp_o[1];
    /* do i=GRIDG_ai3s */
    for (dim=0; dim<2; dim++){
      pc[dim]=_x(grid[GRIDG_ai2(gl,gl.is,j)],dim);
      pp_o[dim]=_x(grid[GRIDG_ai2(gl,gl.is+1,j)],dim);
    }
    grid[GRIDG_ai2(gl,gl.is-1,j)].x=2.0*pc[0]-pp_o[0];
    grid[GRIDG_ai2(gl,gl.is-1,j)].y=2.0*pc[1]-pp_o[1];
  }
}



void RXGrid2D(char **argum, void *action_args, SOAP_codex_t *codex){
  SOAP_codex_t codex2;

  codex2.vars=codex->vars;
  codex2.ACTION=TRUE;
  codex2.action=&actions;
  codex2.action_args=action_args;
  codex2.FUNCTION=TRUE;
  codex2.function=&functions;
  codex2.function_args=action_args;
  codex2.VERBOSE=codex->VERBOSE;
  codex2.filename=codex->filename;

  SOAP_process_code(*argum, &codex2, SOAP_VARS_KEEP_ALL);
  codex->vars=codex2.vars;
}


static void verify_all_nodes_initialized(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid){
  long i,j;

  for(i=gl2d.is; i<=gl2d.ie; i++){
    for(j=gl2d.js; j<=gl2d.je; j++){
      if (!xygrid[GRIDG_ai2(gl2d,i,j)].INIT){
        GRIDG_fatal_error("The x,y values were not provided to node (%ld,%ld). Grid generation failed.",i,j);
      }
    }
  }
}


void GRIDG_read_grid_2D_from_file(char *filename, GRIDG_gl2d_t *gl2d, GRIDG_xygrid_t **xygrid,
                bool VERBOSE, bool *Problem){
  actionsarg2D_t actionsarg2D;

  actionsarg2D.gl2d=gl2d;
  actionsarg2D.xygrid=xygrid;

  read_RX_grid(filename, "Grid", &RXGrid2D, &GRIDG_write_grid_2D_to_file,
             &actionsarg2D, VERBOSE);
  verify_all_nodes_initialized(*gl2d, *xygrid);  
  find_xy_on_ghost_nodes_extrapolate(*xygrid, *gl2d);
}


void GRIDG_read_grid_2D_from_argum(char *argum, SOAP_codex_t *codex, GRIDG_gl2d_t *gl2d, GRIDG_xygrid_t **xygrid){
  actionsarg2D_t actionsarg2D;

  actionsarg2D.gl2d=gl2d;
  actionsarg2D.xygrid=xygrid;

  codex->ACTION=TRUE;
  codex->action=&actions;
  codex->action_args=&actionsarg2D;
  codex->FUNCTION=TRUE;
  codex->function=&functions;
  codex->function_args=&actionsarg2D;
  SOAP_process_code(argum, codex, SOAP_VARS_KEEP_ALL);
  verify_all_nodes_initialized(*gl2d, *xygrid);
  find_xy_on_ghost_nodes_extrapolate(*xygrid, *gl2d);
}


