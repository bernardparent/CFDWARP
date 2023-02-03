// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000,2022 Bernard Parent.
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
#include "3Dbase.h"
#include "3Dextra.h"


typedef struct {
  GRIDG_gl3d_t *gl3d;
  GRIDG_xyzgrid_t **xyzgrid;
} actionsarg3D_t;


static void verify_zone_validity(long is, long js, long ks, long ie, long je, long ke, GRIDG_gl3d_t gl3d, SOAP_codex_t *codex){
  if (ie<is) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be less than is (%ld).",ie,is);
  if (ie>gl3d.ie)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be greater than domain.ie (%ld) as specified in Size().",ie,gl3d.ie);
  if (is<gl3d.is)
    SOAP_fatal_error(codex,"Zone boundaries invalid: is (%ld) can not be less than domain.is (%ld) as specified in Size().",is,gl3d.is);

  if (je<js) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: je (%ld) can not be less than js (%ld).",je,js);
  if (je>gl3d.je)
    SOAP_fatal_error(codex,"Zone boundaries invalid: je (%ld) can not be greater than domain.je (%ld) as specified in Size().",je,gl3d.je);
  if (js<gl3d.js)
    SOAP_fatal_error(codex,"Zone boundaries invalid: js (%ld) can not be less than domain.js (%ld) as specified in Size().",js,gl3d.js);

  if (ke<ks) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: ke (%ld) can not be less than ks (%ld).",ke,ks);
  if (ke>gl3d.ke)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ke (%ld) can not be greater than domain.ke (%ld) as specified in Size().",ke,gl3d.ke);
  if (ks<gl3d.ks)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ks (%ld) can not be less than domain.ks (%ld) as specified in Size().",ks,gl3d.ks);

}


static void verify_node_validity(long i, long j, long k, GRIDG_gl3d_t gl3d, SOAP_codex_t *codex){
  if (i>gl3d.ie)
    SOAP_fatal_error(codex,"Node location invalid: i (%ld) can not be greater than domain.ie (%ld) as specified in Size().",i,gl3d.ie);
  if (i<gl3d.is)
    SOAP_fatal_error(codex,"Node location invalid: i (%ld) can not be less than domain.is (%ld) as specified in Size().",i,gl3d.is);

  if (j>gl3d.je)
    SOAP_fatal_error(codex,"Node location invalid: j (%ld) can not be greater than domain.je (%ld) as specified in Size().",j,gl3d.je);
  if (j<gl3d.js)
    SOAP_fatal_error(codex,"Node location invalid: j (%ld) can not be less than domain.js (%ld) as specified in Size().",j,gl3d.js);

  if (k>gl3d.ke)
    SOAP_fatal_error(codex,"Node location invalid: k (%ld) can not be greater than domain.ke (%ld) as specified in Size().",k,gl3d.ke);
  if (k<gl3d.ks)
    SOAP_fatal_error(codex,"Node location invalid: k (%ld) can not be less than domain.ks (%ld) as specified in Size().",k,gl3d.ks);
}



void GRIDG_write_grid_3D_to_file(FILE **controlfile){
  fprintf(*controlfile," \n");
  fprintf(*controlfile,
    "\n"
    "is=1;\n"
    "ie=30;\n"
    "js=1;\n"
    "je=30;\n"
    "ks=1;\n"
    "ke=30;\n"
    "\n"
    "Grid(\n"
    "  Size(is,js,ks, ie,je,ke);\n"
    "  Point(is,js,ks, 0.0e0,0.0e0,0.0e0);\n"
    "  Point(ie,js,ks, 1.0e0,0.0e0,0.0e0);\n"
    "  Point(ie,je,ks, 1.0e0,1.0e0,0.0e0);\n"
    "  Point(is,je,ks, 0.0e0,1.0e0,0.0e0);\n"
    "  Point(is,js,ke, 0.0e0,0.0e0,1.0e0);\n"
    "  Point(ie,js,ke, 1.0e0,0.0e0,1.0e0);\n"
    "  Point(ie,je,ke, 1.0e0,1.0e0,1.0e0);\n"
    "  Point(is,je,ke, 0.0e0,1.0e0,1.0e0);\n"
    "  JoinCorners(is,js,ks, ie,je,ke,  EE,0.5e0,1.0e0,1.0e0,  EE,0.5e0,1.0e0,1.0e0,  EE,0.5e0,1.0e0,1.0e0);\n"
    ");\n"
  );
}



static double find_distance(GRIDG_xyzgrid_t *xyzgrid, GRIDG_gl3d_t gl3d,
                    long i1, long j1, long k1,
                    long i2, long j2, long k2){
  double dist;
  dist=sqrt(sqr(xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].x-xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].x)
           +sqr(xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].y-xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].y)
           +sqr(xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].z-xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].z));
  return(dist);
}


void Size3D(GRIDG_gl3d_t *gl, GRIDG_xyzgrid_t **xyzgrid,
            long is, long js, long ks, long ie, long je, long ke){
  long i,j,k;

  gl->is=is;  gl->js=js;  gl->ks=ks;
  gl->ie=ie;  gl->je=je;  gl->ke=ke;

  *xyzgrid = (GRIDG_xyzgrid_t *) malloc(
          (gl->ie-gl->is+1+4)*(gl->je-gl->js+1+4)*(gl->ke-gl->ks+1+4)
           *sizeof(GRIDG_xyzgrid_t));
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        (*xyzgrid)[GRIDG_ai3(*gl,i,j,k)].INIT=FALSE;
      }
    }
  }

}


static void Size_argum(char **argum, SOAP_codex_t *codex,  GRIDG_gl3d_t *gl3d, GRIDG_xyzgrid_t **xyzgrid){
  long is,js,ks,ie,je,ke;
  int eos=EOS;
  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld%n",&is,&js,&ks,&ie,&je,&ke,&eos)!=6 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Size().");
  }
  if (is>=ie) SOAP_fatal_error(codex,"In Size(), is (%ld) must be less than ie (%ld).",is,ie);
  if (js>=je) SOAP_fatal_error(codex,"In Size(), js (%ld) must be less than je (%ld).",js,je);
  if (ks>=ke) SOAP_fatal_error(codex,"In Size(), ks (%ld) must be less than ke (%ld).",ks,ke);
  Size3D(gl3d, xyzgrid, is, js, ks, ie, je, ke);
}


void Point3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid,
             long i, long j, long k, double x, double y, double z){
  xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=x;
  xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=y;
  xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=z;
  xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT=TRUE;
}


static void Point_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long i,j,k;
  double x,y,z;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%lg,%lg,%lg%n",&i,&j,&k,&x,&y,&z,&eos)!=6 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Point().");
  }
  verify_node_validity(i, j, k, gl3d, codex);
  Point3D(gl3d, xyzgrid, i, j, k, x, y, z);
}


void Corners3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid,
                   long i1, long j1, long k1, long i2, long j2, long k2,
                   double x1, double y1, double z1, double x2, double y2, double z2){
  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].x=x1;  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].y=y1;  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].z=z1;  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].INIT=TRUE;
  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k1)].x=x2;  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k1)].y=y1;  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k1)].z=z1;  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k1)].INIT=TRUE;
  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k1)].x=x2;  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k1)].y=y2;  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k1)].z=z1;  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k1)].INIT=TRUE;
  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k1)].x=x1;  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k1)].y=y2;  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k1)].z=z1;  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k1)].INIT=TRUE;
  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k2)].x=x1;  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k2)].y=y1;  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k2)].z=z2;  xyzgrid[GRIDG_ai3(gl3d,i1,j1,k2)].INIT=TRUE;
  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k2)].x=x2;  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k2)].y=y1;  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k2)].z=z2;  xyzgrid[GRIDG_ai3(gl3d,i2,j1,k2)].INIT=TRUE;
  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].x=x2;  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].y=y2;  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].z=z2;  xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].INIT=TRUE;
  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k2)].x=x1;  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k2)].y=y2;  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k2)].z=z2;  xyzgrid[GRIDG_ai3(gl3d,i1,j2,k2)].INIT=TRUE;
}


static void Corners_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long i1,j1,k1,i2,j2,k2;
  double x1,y1,z1,x2,y2,z2;
  int eos=EOS;
  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%lg,%lg,%lg,%lg,%lg,%lg%n",
         &i1,&j1,&k1,&i2,&j2,&k2,&x1,&y1,&z1,&x2,&y2,&z2,&eos)!=12 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Corners().");
  }
  verify_zone_validity(i1, j1, k1, i2, j2, k2, gl3d, codex);
  Corners3D(gl3d, xyzgrid,i1,j1,k1,i2,j2,k2,x1,y1,z1,x2,y2,z2);
}


static int find_xyz_along_ijk_local(GRIDG_xyzgrid_t *xyzgrid, double segval1, double segval2,
                     double NLoverN, char modeL, char modeH,
                     GRIDG_gl3d_t gl3d, SOAP_codex_t *codex, long i1, long j1, long k1, long i2, long j2, long k2,
                     int(*find_segment)(double *, long , double,
                     double , double , double, SOAP_codex_t * )){
  int err;
  if (!xyzgrid[GRIDG_ai3(gl3d,i1,j1,k1)].INIT){
    SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Grid segment could not be formed.",i1,j1,k1);
  }
  if (!xyzgrid[GRIDG_ai3(gl3d,i2,j2,k2)].INIT){
    SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Grid segment could not be formed.",i2,j2,k2);
  }

  if (modeL=='G' || modeL=='g') {
    if (!xyzgrid[GRIDG_ai3(gl3d,i1+sgndif(i1,i2),j1+sgndif(j1,j2),k1+sgndif(k1,k2))].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Grid segment could not be formed.",i1+sgndif(i1,i2),j1+sgndif(j1,j2),k1+sgndif(k1,k2));
    }
    segval1=find_distance(xyzgrid, gl3d, i1, j1, k1, i1+sgndif(i1,i2), j1+sgndif(j1,j2), k1+sgndif(k1,k2));
  }
  if (modeH=='G' || modeH=='g') {
    if (!xyzgrid[GRIDG_ai3(gl3d,i2-sgndif(i1,i2),j2-sgndif(j1,j2),k2-sgndif(k1,k2))].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Grid segment could not be formed.",i2-sgndif(i1,i2),j2-sgndif(j1,j2),k2-sgndif(k1,k2));
    }
    segval2=find_distance(xyzgrid, gl3d, i2, j2, k2, i2-sgndif(i1,i2), j2-sgndif(j1,j2), k2-sgndif(k1,k2));
  }

  err=find_xyz_along_ijk(xyzgrid, segval1, segval2, NLoverN, gl3d, i1, j1, k1, i2, j2, k2, codex, find_segment);
  return(err);
}


int JoinCorners3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
                 long i1, long j1, long k1, long i2, long j2, long k2,
                 char modeL_i, char modeH_i, double NLoverN_i, double valL_i, double valH_i,
                 char modeL_j, char modeH_j, double NLoverN_j, double valL_j, double valH_j,
                 char modeL_k, char modeH_k, double NLoverN_k, double valL_k, double valH_k){
  long i,j;
  int err;
  int (*FindSegment_i)(double *, long , double, double , double , double , SOAP_codex_t *);
  int (*FindSegment_j)(double *, long , double, double , double , double , SOAP_codex_t *);
  int (*FindSegment_k)(double *, long , double, double , double , double , SOAP_codex_t *);

  find_segment_type(&FindSegment_i, modeL_i, modeH_i, codex);
  find_segment_type(&FindSegment_j, modeL_j, modeH_j, codex);
  find_segment_type(&FindSegment_k, modeL_k, modeH_k, codex);

  err=find_xyz_along_ijk_local(xyzgrid, valL_i, valH_i,
                  NLoverN_i, modeL_i, modeH_i, gl3d, codex, 
                  i1, j1, k1, i2, j1, k1, FindSegment_i);
  err+=find_xyz_along_ijk_local(xyzgrid, valL_i, valH_i,
                  NLoverN_i, modeL_i, modeH_i, gl3d, codex, 
                  i1, j2, k1, i2, j2, k1, FindSegment_i);
  for (i=i1; i<=i2 && err==0; i++){
    err+=find_xyz_along_ijk_local(xyzgrid, valL_j, valH_j,
                    NLoverN_j, modeL_j, modeH_j, gl3d, codex, 
                    i, j1, k1, i, j2, k1, FindSegment_j);
  }
  err+=find_xyz_along_ijk_local(xyzgrid, valL_i, valH_i,
                  NLoverN_i, modeL_i, modeH_i, gl3d, codex, 
                  i1, j1, k2, i2, j1, k2, FindSegment_i);
  err+=find_xyz_along_ijk_local(xyzgrid, valL_i, valH_i,
                  NLoverN_i, modeL_i, modeH_i, gl3d, codex, 
                  i1, j2, k2, i2, j2, k2, FindSegment_i);
  for (i=i1; i<=i2 && err==0; i++){
    err+=find_xyz_along_ijk_local(xyzgrid, valL_j, valH_j,
                    NLoverN_j, modeL_j, modeH_j, gl3d, codex, 
                    i, j1, k2, i, j2, k2, FindSegment_j);
  }

  for (i=i1; i<=i2 && err==0; i++){
    for (j=j1; j<=j2 && err==0; j++){
      err+=find_xyz_along_ijk_local(xyzgrid, valL_k, valH_k,
                      NLoverN_k, modeL_k, modeH_k, gl3d, codex, 
                      i, j, k1, i, j, k2, FindSegment_k);
    }
  }
  return(err);
}


static void JoinCorners_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long i1,j1,k1,i2,j2,k2,cnt;
  int err;
  char modeL_i,modeH_i,modeL_j,modeH_j,modeL_k,modeH_k;
  double NLoverN_i,valL_i,valH_i,NLoverN_j,valL_j,valH_j,
         NLoverN_k,valL_k,valH_k;
  int eos=EOS;

  modeL_i='0';
  modeL_j='0';
  modeL_k='0';
  modeH_i='0';
  modeH_j='0';
  modeH_k='0';
  for (cnt=0; cnt<SOAP_number_argums(*argum); cnt++)
    if (cnt!=6 && cnt!=10 && cnt!=14) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,"
               "%c%c,%lg,%lg,%lg,"
               "%c%c,%lg,%lg,%lg,"
               "%c%c,%lg,%lg,%lg%n",
         &i1,&j1,&k1,&i2,&j2,&k2,
         &modeL_i,&modeH_i,&NLoverN_i,&valL_i,&valH_i,
         &modeL_j,&modeH_j,&NLoverN_j,&valL_j,&valH_j,
         &modeL_k,&modeH_k,&NLoverN_k,&valL_k,&valH_k,&eos)!=21 || (*argum)[eos]!=EOS) {
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
  if ((k1==k2) && !(modeL_k=='N' && modeH_k=='O')) {
     modeL_k='N';
     modeH_k='O';
  }
  verify_zone_validity(i1, j1, k1, i2, j2, k2, gl3d, codex);
  err=JoinCorners3D(gl3d, xyzgrid, codex, i1, j1, k1, i2, j2, k2,
              modeL_i, modeH_i, NLoverN_i, valL_i, valH_i,
              modeL_j, modeH_j, NLoverN_j, valL_j, valH_j,
              modeL_k, modeH_k, NLoverN_k, valL_k, valH_k);
  if (err!=0) SOAP_fatal_error(codex,"Error while executing JoinCorners() command.");
}


int Join3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long i1, long j1, long k1, long i2, long j2,
           long k2, char index, char modeL, char modeH,
           double NLoverN, double valL, double valH){
  long i,j,k;
  int (*FindSegment)(double *, long , double, double , double , double, SOAP_codex_t * );
  int err;
  find_segment_type(&FindSegment, modeL, modeH, codex);

  err=0;
  switch (index) {

    case 'i':
      for (j=j1; j<=j2 && err==0; j++){
        for (k=k1; k<=k2 && err==0; k++){
          err+=find_xyz_along_ijk_local(xyzgrid, valL, valH,
                          NLoverN, modeL, modeH, gl3d, codex,
                          i1, j, k, i2, j, k, FindSegment);
        }
      }
    break;
    case 'j':
      for (i=i1; i<=i2 && err==0; i++){
        for (k=k1; k<=k2 && err==0; k++){
          err+=find_xyz_along_ijk_local(xyzgrid, valL, valH,
                          NLoverN, modeL, modeH, gl3d, codex,
                          i, j1, k, i, j2, k, FindSegment);
        }
      }
    break;

    case 'k':
      for (i=i1; i<=i2 && err==0; i++){
        for (j=j1; j<=j2 && err==0; j++){
          err+=find_xyz_along_ijk_local(xyzgrid, valL, valH,
                          NLoverN, modeL, modeH, gl3d, codex,
                          i, j, k1, i, j, k2, FindSegment);
        }
      }
    break;
  }
  return(err);
}


static void Join_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long i1,j1,k1,i2,j2,k2,cnt;
  char index,modeL,modeH;
  double NLoverN,valL,valH;
  int err;
  int eos=EOS;
  for (cnt=0; cnt<11; cnt++) if (cnt!=6 && cnt!=7) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%c,%c%c,%lg,%lg,%lg%n",
         &i1,&j1,&k1,&i2,&j2,&k2,&index,&modeL,&modeH,&NLoverN,
         &valL,&valH,&eos)!=12 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Join().");
  }
  verify_zone_validity(i1, j1, k1, i2, j2, k2, gl3d, codex);
  if (index!='i' && index!='j' && index!='k') SOAP_fatal_error(codex,"In Join(), the index name must be either 'i', 'j', or 'k'.");

  err=Join3D(gl3d, xyzgrid, codex,
            i1, j1, k1, i2, j2, k2, index, modeL, modeH,
            NLoverN, valL, valH);
  if (err!=0) SOAP_fatal_error(codex,"Error occurred while executing Join() command.");

}


#define x(i,j,k)   xyz[GRIDG_ai3(gl3d,i,j,k)].x
#define y(i,j,k)   xyz[GRIDG_ai3(gl3d,i,j,k)].y
#define z(i,j,k)   xyz[GRIDG_ai3(gl3d,i,j,k)].z

/* algorithm by Derrick Alexander */
void JoinFaces3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyz, SOAP_codex_t *codex,
               long i1, long j1, long k1, long i2, long j2, long k2){
  long i,j,k,cnt;
  bool VALID;
  double ax_y,ax_z,ay_x,ay_z,az_x,az_y;
  double fx_top,fx_bottom,fx_left,fx_right,fy_start,fy_end,fy_bottom,fy_top,
         fz_left,fz_right,fz_start,fz_end;
  double x_denom,y_denom1,y_denom2,y_denom,y1,y2,y3,y4,y5,y_other,x1,x2,x3,x4,yx,z1,z2;
  double li1,li2,lj1,lj2,lk1,lk2,idiff,jdiff,kdiff;

  VALID=FALSE;
  for (cnt=1; cnt<4; cnt++){
    for (i=i1+1; i<i2; i++){
      for (j=j1+1; j<j2; j++){
        for (k=k1+1; k<k2; k++){
          VALID=TRUE;
          if (!xyz[GRIDG_ai3(gl3d,i1,j,k)].INIT){
            SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. JoinFaces command failed.",i1,j,k);
          }
          if (!xyz[GRIDG_ai3(gl3d,i2,j,k)].INIT){
            SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. JoinFaces command failed.",i2,j,k);
          }
          if (!xyz[GRIDG_ai3(gl3d,i,j1,k)].INIT){
            SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. JoinFaces command failed.",i,j1,k);
          }
          if (!xyz[GRIDG_ai3(gl3d,i,j2,k)].INIT){
            SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. JoinFaces command failed.",i,j2,k);
          }
          if (!xyz[GRIDG_ai3(gl3d,i,j,k1)].INIT){
            SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. JoinFaces command failed.",i,j,k1);
          }
          if (!xyz[GRIDG_ai3(gl3d,i,j,k2)].INIT){
            SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. JoinFaces command failed.",i,j,k2);
          }

          if (cnt==1) {
            ax_y=0.5;
            ax_z=0.5;
            ay_x=0.5;
            ay_z=0.5;
            az_x=0.5;
            az_y=0.5;
          } else {
            li1=sqrt(sqr(x(i,j,k)-x(i1,j,k))+sqr(y(i,j,k)-y(i1,j,k))+sqr(z(i,j,k)-z(i1,j,k)));
            li2=sqrt(sqr(x(i,j,k)-x(i2,j,k))+sqr(y(i,j,k)-y(i2,j,k))+sqr(z(i,j,k)-z(i2,j,k)));

            lj1=sqrt(sqr(x(i,j,k)-x(i,j1,k))+sqr(y(i,j,k)-y(i,j1,k))+sqr(z(i,j,k)-z(i,j1,k)));
            lj2=sqrt(sqr(x(i,j,k)-x(i,j2,k))+sqr(y(i,j,k)-y(i,j2,k))+sqr(z(i,j,k)-z(i,j2,k)));

            lk1=sqrt(sqr(x(i,j,k)-x(i,j,k1))+sqr(y(i,j,k)-y(i,j,k1))+sqr(z(i,j,k)-z(i,j,k1)));
            lk2=sqrt(sqr(x(i,j,k)-x(i,j,k2))+sqr(y(i,j,k)-y(i,j,k2))+sqr(z(i,j,k)-z(i,j,k2)));

            if (li1<=li2) idiff=li1; else idiff=li2;
            if (lj1<=lj2) jdiff=lj1; else jdiff=lj2;
            if (lk1<=lk2) kdiff=lk1; else kdiff=lk2;

            /* x factors */
            ax_y=kdiff/(jdiff+kdiff);
            ax_z=1.0e0-ax_y;
            /* y factors */
            ay_x=kdiff/(idiff+kdiff);
            ay_z=1.0e0-ay_x;
            /* z factors */
            az_x=jdiff/(idiff+jdiff);
            az_y=1.0e0-az_x;
          }

          /* x factors */
          fx_top=   (x(i ,j ,k2)-x(i1,j ,k2))/(x(i2,j ,k2)-x(i1,j ,k2));
          fx_bottom=(x(i ,j ,k1)-x(i1,j ,k1))/(x(i2,j ,k1)-x(i1,j ,k1));
          fx_left=  (x(i ,j1,k )-x(i1,j1,k ))/(x(i2,j1,k )-x(i1,j1,k ));
          fx_right= (x(i ,j2,k )-x(i1,j2,k ))/(x(i2,j2,k )-x(i1,j2,k ));
          /* y factors */
          fy_start= (y(i1,j ,k )-y(i1,j1,k ))/(y(i1,j2,k )-y(i1,j1,k ));
          fy_end=   (y(i2,j ,k )-y(i2,j1,k ))/(y(i2,j2,k )-y(i2,j1,k ));
          fy_bottom=(y(i ,j ,k1)-y(i ,j1,k1))/(y(i ,j2,k1)-y(i ,j1,k1));
          fy_top=   (y(i ,j ,k2)-y(i ,j1,k2))/(y(i ,j2,k2)-y(i ,j1,k2));
          /* z factors */
          fz_left=  (z(i ,j1,k )-z(i ,j1,k1))/(z(i ,j1,k2)-z(i ,j1,k1));
          fz_right= (z(i ,j2,k )-z(i ,j2,k1))/(z(i ,j2,k2)-z(i ,j2,k1));
          fz_start= (z(i1,j ,k )-z(i1,j ,k1))/(z(i1,j ,k2)-z(i1,j ,k1));
          fz_end=   (z(i2,j ,k )-z(i2,j ,k1))/(z(i2,j ,k2)-z(i2,j ,k1));

          x_denom=1.0e0-(ax_z*(fx_top-fx_bottom)*az_x*(fz_end-fz_start));

          y_denom1=ay_x*(fy_end-fy_start)+ay_z*(fy_top-fy_bottom)*az_x*(fz_end-fz_start);
          y_denom2=ax_z*(fx_top-fx_bottom)*az_y*(fz_right-fz_left)+ax_y*(fx_right-fx_left);
          y_denom=1.0e0-(y_denom1*y_denom2)/x_denom-ay_z*(fy_top-fy_bottom)*az_y*(fz_right-fz_left);

          y1=ay_x*(fy_end-fy_start)+ay_z*(fy_top-fy_bottom)*az_x*(fz_end-fz_start);
          y2=(y(i,j2,k)-y(i,j1,k))*(-x(i1,j,k))/(x(i2,j,k)-x(i1,j,k));
          y3=ay_x*fy_start;
          y4=ay_z*(fy_top-fy_bottom)*az_y*(fz_right-fz_left)*(-y(i,j1,k));
          y5=ay_z*(fy_top-fy_bottom)*(az_y*fz_left+az_x*fz_start)+ay_z*fy_bottom;
          y_other=y1*y2+y4+(y3+y5)*(y(i,j2,k)-y(i,j1,k))+y(i,j1,k);

          x1=(ax_z*(fx_top-fx_bottom)*az_y*(fz_right-fz_left)+ax_y*(fx_right-fx_left));
          x2=(x(i2,j,k)-x(i1,j,k))/(y(i,j2,k)-y(i,j1,k));
          x3=ax_z*(fx_top-fx_bottom)*az_x*(fz_end-fz_start)*(-x(i1,j,k));
          x4=ax_z*(fx_top-fx_bottom)*(az_y*fz_left+az_x*fz_start)+ax_z*fx_bottom+ax_y*fx_left;
          yx=(y_denom1/x2)*(x1*x2*(-y(i,j1,k))+x3+x4*(x(i2,j,k)-x(i1,j,k))+x(i1,j,k))/x_denom;

          y(i,j,k)=(y_other+yx)/y_denom;

          x(i,j,k)=((y(i,j,k)-y(i,j1,k))*x1*x2+x3+x4*(x(i2,j,k)-x(i1,j,k))+x(i1,j,k))/x_denom;

          z1=az_y*((fz_right-fz_left)*(y(i,j,k)-y(i,j1,k))/(y(i,j2,k)-y(i,j1,k))+fz_left);
          z2=az_x*((fz_end-fz_start)*(x(i,j,k)-x(i1,j,k))/(x(i2,j,k)-x(i1,j,k))+fz_start);
          z(i,j,k)=(z1+z2)*(z(i,j,k2)-z(i,j,k1))+z(i,j,k1);
          xyz[GRIDG_ai3(gl3d,i,j,k)].INIT=TRUE;
        }
      }
    }
  }

  if (!VALID){
    SOAP_fatal_error(codex,"When using the JoinFaces() command, set i2>i1+1 and j2>j1+1 and k2>k1+1.");
  }
}

static void JoinFaces_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long i1,j1,k1,i2,j2,k2;
  int eos=EOS;
  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld%n",&i1,&j1,&k1,&i2,&j2,&k2,&eos)!=6 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in JoinFaces().");
  }
  verify_zone_validity(i1, j1, k1, i2, j2, k2, gl3d, codex);
  JoinFaces3D(gl3d, xyzgrid, codex, i1, j1, k1, i2, j2, k2);
}


static double _distance(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyz, long l1, long l2){
  double dist;
  dist=sqrt(sqr(xyz[l2].x-xyz[l1].x)
           +sqr(xyz[l2].y-xyz[l1].y)
           +sqr(xyz[l2].z-xyz[l1].z));
  return(dist);
}


void Plane(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyz, SOAP_codex_t *codex,
           long i1, long j1, long k1, long i2, long j2, long k2){
  bool VALID;
  long dim1,dim2,l1,l2,l11,l12,l22,l21,cnt1,cnt2,l;
  double fact;
  EXM_vec3D_t vec21start,vec21end,vec21,vecstart,vecend,vecstartrot,vecendrot;
  EXM_mat3x3_t Rstart,Rend;
  double scalestart,scaleend;

  VALID=FALSE;
  l11=0; //avoid compiler warning 
  l12=0; //avoid compiler warning
  l21=0; //avoid compiler warning
  l22=0; //avoid compiler warning
  if (k2==k1){
    VALID=TRUE;
    dim1=0;
    dim2=1;  
    l11=GRIDG_ai3(gl3d, i1,  j1,  k1);
    l12=GRIDG_ai3(gl3d, i1,  j2,  k1);
    l21=GRIDG_ai3(gl3d, i2,  j1,  k1);
    l22=GRIDG_ai3(gl3d, i2,  j2,  k1);
  }
  if (j2==j1){
    VALID=TRUE;
    dim1=0;
    dim2=2;  
    l11=GRIDG_ai3(gl3d, i1,  j1,  k1);
    l12=GRIDG_ai3(gl3d, i1,  j1,  k2);
    l21=GRIDG_ai3(gl3d, i2,  j1,  k1);
    l22=GRIDG_ai3(gl3d, i2,  j1,  k2);
  }
  if (i2==i1){
    VALID=TRUE;
    dim1=1;
    dim2=2;  
    l11=GRIDG_ai3(gl3d, i1,  j1,  k1);
    l12=GRIDG_ai3(gl3d, i1,  j1,  k2);
    l21=GRIDG_ai3(gl3d, i1,  j2,  k1);
    l22=GRIDG_ai3(gl3d, i1,  j2,  k2);
  }
  if (!VALID){
    SOAP_fatal_error(codex,"When using the Plane() command, either set i1=i2 or j1=j2 or k1=k2.");
  }
  
  // test node numbering
  /*
  long i,j,k;
  for (i=gl3d.is; i<gl3d.ie; i++){
    for (j=gl3d.js; j<gl3d.je; j++){
     for (k=gl3d.ks; k<gl3d.ke; k++){
       l=GRIDG_ai3(gl3d,i,j,k);
       GRIDG_find_ijk_from_l(gl3d, l, &i2, &j2, &k2); 
       if (i!=i2 || j!=j2 || k!=k2) SOAP_fatal_error(codex,"Problem i=%ld j=%ld k=%ld i2=%ld j2=%ld k2=%ld\n",i,j,k,i2,j2,k2);
     }
    }
  }
  */
  
  // check if nodes have been initialized around the plane 
  for (cnt1=0; GRIDG_al3(gl3d,l11,dim1,cnt1+1)!=l21; cnt1++) {    
    l=GRIDG_al3(gl3d,l11,dim1,cnt1);
    GRIDG_find_ijk_from_l(gl3d, l, &i2, &j2, &k2);
    if (!xyz[l].INIT) 
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) needs to be initialized prior to calling Plane().\n",i2,j2,k2);
    l=GRIDG_al3(gl3d,l12,dim1,cnt1);
    GRIDG_find_ijk_from_l(gl3d, l, &i2, &j2, &k2);
    if (!xyz[l].INIT) 
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) needs to be initialized prior to calling Plane().\n",i2,j2,k2);
  }
  for (cnt2=0; GRIDG_al3(gl3d,l11,dim2,cnt2+1)!=l12; cnt2++) {    
    l=GRIDG_al3(gl3d,l11,dim2,cnt2);
    GRIDG_find_ijk_from_l(gl3d, l, &i2, &j2, &k2);
    if (!xyz[l].INIT) 
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) needs to be initialized prior to calling Plane().\n",i2,j2,k2);
    l=GRIDG_al3(gl3d,l21,dim2,cnt2);
    GRIDG_find_ijk_from_l(gl3d, l, &i2, &j2, &k2);
    if (!xyz[l].INIT) 
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) needs to be initialized prior to calling Plane().\n",i2,j2,k2);
  }


  vec21start[0]=xyz[l12].x-xyz[l11].x;
  vec21start[1]=xyz[l12].y-xyz[l11].y;
  vec21start[2]=xyz[l12].z-xyz[l11].z;
  vec21end[0]=xyz[l22].x-xyz[l21].x;
  vec21end[1]=xyz[l22].y-xyz[l21].y;
  vec21end[2]=xyz[l22].z-xyz[l21].z;

  for (cnt1=1; GRIDG_al3(gl3d,l11,dim1,cnt1)!=l21; cnt1++){
    l1=GRIDG_al3(gl3d,l11,dim1,cnt1);
    l2=GRIDG_al3(gl3d,l12,dim1,cnt1);
    fact=_distance(gl3d,xyz,l1,l11)/(1e-20+_distance(gl3d,xyz,l21,l11));
    vec21[0]=xyz[l2].x-xyz[l1].x;
    vec21[1]=xyz[l2].y-xyz[l1].y;
    vec21[2]=xyz[l2].z-xyz[l1].z;
    assert(EXM_vector_magnitude(vec21start)!=0.0);
    assert(EXM_vector_magnitude(vec21end)!=0.0);
    assert(EXM_vector_magnitude(vec21)!=0.0);
    scalestart=EXM_vector_magnitude(vec21)/EXM_vector_magnitude(vec21start);
    scaleend=EXM_vector_magnitude(vec21)/EXM_vector_magnitude(vec21end);
    EXM_find_rotation_matrix(vec21start, vec21, Rstart);
    EXM_find_rotation_matrix(vec21end, vec21, Rend);

    for (cnt2=1; GRIDG_al3(gl3d,l1,dim2,cnt2)!=l2; cnt2++){
      vecstart[0]=xyz[GRIDG_al3(gl3d,l11,dim2,cnt2)].x-xyz[l11].x;
      vecstart[1]=xyz[GRIDG_al3(gl3d,l11,dim2,cnt2)].y-xyz[l11].y;
      vecstart[2]=xyz[GRIDG_al3(gl3d,l11,dim2,cnt2)].z-xyz[l11].z;
      EXM_multiply_matrix_vector(Rstart,vecstart,vecstartrot);  
      vecend[0]=xyz[GRIDG_al3(gl3d,l21,dim2,cnt2)].x-xyz[l21].x;
      vecend[1]=xyz[GRIDG_al3(gl3d,l21,dim2,cnt2)].y-xyz[l21].y;
      vecend[2]=xyz[GRIDG_al3(gl3d,l21,dim2,cnt2)].z-xyz[l21].z;
      EXM_multiply_matrix_vector(Rend,vecend,vecendrot);  

      xyz[GRIDG_al3(gl3d,l1,dim2,cnt2)].x=xyz[l1].x+(1.0-fact)*vecstartrot[0]*scalestart+fact*vecendrot[0]*scaleend;
      xyz[GRIDG_al3(gl3d,l1,dim2,cnt2)].y=xyz[l1].y+(1.0-fact)*vecstartrot[1]*scalestart+fact*vecendrot[1]*scaleend;
      xyz[GRIDG_al3(gl3d,l1,dim2,cnt2)].z=xyz[l1].z+(1.0-fact)*vecstartrot[2]*scalestart+fact*vecendrot[2]*scaleend;

      xyz[GRIDG_al3(gl3d,l1,dim2,cnt2)].INIT=TRUE;    
    }
  }  
}


static void Plane_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long i1,j1,k1,i2,j2,k2;
  int eos=EOS;
  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld%n",&i1,&j1,&k1,&i2,&j2,&k2,&eos)!=6 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Plane().");
  }
  verify_zone_validity(i1, j1, k1, i2, j2, k2, gl3d, codex);
  Plane(gl3d, xyzgrid, codex, i1, j1, k1, i2, j2, k2);
}


void Translate3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           double dx, double dy, double dz){
  long i,j,k;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Translate command failed.",i,j,k);
        }
        xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x+dx;
        xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y+dy;
        xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z+dz;
      }
    }
  }
}


static void translate_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke;
  double dx,dy,dz;
  int eos=EOS;
  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%lg,%lg,%lg%n",
               &is,&js,&ks,&ie,&je,&ke,&dx,&dy,&dz,&eos)!=9 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in translate_argum().");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);
  Translate3D(gl3d, xyzgrid, codex, is, js, ks, ie, je, ke, dx, dy, dz);
}


void Rotate3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           double pivot_x, double pivot_y, double pivot_z,
           double angle, char axis){
  long i,j,k;
  double r,phi,dx,dy,dz;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Rotate command failed.",i,j,k);
        }
        if (axis=='z'){
          dx=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x-pivot_x;
          dy=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y-pivot_y;
          r=sqrt(+sqr(dx)+sqr(dy));
          if (dx==0.0e0) phi=0.5e0*pi*sign(dy); else phi=atan(dy/dx);
          if (dx<0.0e0) phi=phi+pi;
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=pivot_x+r*cos(phi+angle);
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=pivot_y+r*sin(phi+angle);
        }
        if (axis=='x'){
          dy=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y-pivot_y;
          dz=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z-pivot_z;
          r=sqrt(+sqr(dy)+sqr(dz));
          if (dy==0.0e0) phi=0.5e0*pi*sign(dz); else phi=atan(dz/dy);
          if (dy<0.0e0) phi=phi+pi;
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=pivot_y+r*cos(phi+angle);
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=pivot_z+r*sin(phi+angle);
        }
        if (axis=='y'){
          dz=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z-pivot_z;
          dx=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x-pivot_x;
          r=sqrt(+sqr(dz)+sqr(dx));
          if (dz==0.0e0) phi=0.5e0*pi*sign(dx); else phi=atan(dx/dz);
          if (dz<0.0e0) phi=phi+pi;
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=pivot_z+r*cos(phi+angle);
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=pivot_x+r*sin(phi+angle);
        }

      }
    }
  }
}



static void Rotate_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke,cnt;
  double angle,pivot_x,pivot_y,pivot_z;
  char axis;
  int eos=EOS;
  for (cnt=0; cnt<10; cnt++) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%lg,%lg,%lg,%lg,%c%n",
         &is,&js,&ks,&ie,&je,&ke,&pivot_x,&pivot_y,&pivot_z,&angle,&axis,&eos)!=11 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Rotate().");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);
  if (axis!='x' && axis!='y' && axis!='z') SOAP_fatal_error(codex,"In Rotate(), the axis name must be either 'x', 'y', or 'z'.");
  Rotate3D(gl3d, xyzgrid, codex, is, js, ks, ie, je, ke,
                pivot_x, pivot_y, pivot_z, angle, axis);
}


void Mirror3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           char axis_name, double axis_pos){
  long i,j,k;

  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Mirror command failed.",i,j,k);
        }
        if (axis_name=='x') {
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=2.0e0*axis_pos-xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x;
        }
        if (axis_name=='y') {
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=2.0e0*axis_pos-xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y;
        }
        if (axis_name=='z') {
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=2.0e0*axis_pos-xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z;
        }
      }
    }
  }
}


static void Mirror_argum(char **argum, SOAP_codex_t *codex,
                         GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke,cnt;
  char axis_name;
  double axis_pos;
  int eos=EOS;

  for (cnt=0; cnt<8; cnt++) if (cnt!=6) SOAP_substitute_argum(argum, cnt, codex);

  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%c,%lg%n",
         &is,&js,&ks,&ie,&je,&ke,&axis_name,&axis_pos,&eos)!=8 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Mirror().");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);
  if (axis_name!='x' && axis_name!='y' && axis_name!='z') SOAP_fatal_error(codex,"In Mirror(), the axis name must be either 'x', 'y', or 'z'.");
  Mirror3D(gl3d, xyzgrid, codex, is, js, ks, ie, je, ke, axis_name, axis_pos);
}




void ReverseIndex3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           char index){
  long i,j,k;
  GRIDG_xyzgrid_t tmp;

  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. ReverseIndex command failed.",i,j,k);
        }
        if (index=='i' && i<ie-(i-is)) {
          tmp=xyzgrid[GRIDG_ai3(gl3d,i,j,k)];
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)]=xyzgrid[GRIDG_ai3(gl3d,ie-(i-is),j,k)];
          xyzgrid[GRIDG_ai3(gl3d,ie-(i-is),j,k)]=tmp;
        }
        if (index=='j' && j<je-(j-js)) {
          tmp=xyzgrid[GRIDG_ai3(gl3d,i,j,k)];
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)]=xyzgrid[GRIDG_ai3(gl3d,i,je-(j-js),k)];
          xyzgrid[GRIDG_ai3(gl3d,i,je-(j-js),k)]=tmp;
        }
        if (index=='k' && k<ke-(k-ks)) {
          tmp=xyzgrid[GRIDG_ai3(gl3d,i,j,k)];
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)]=xyzgrid[GRIDG_ai3(gl3d,i,j,ke-(k-ks))];
          xyzgrid[GRIDG_ai3(gl3d,i,j,ke-(k-ks))]=tmp;
        }
      }
    }
  }
}


static void ReverseIndex_argum(char **argum, SOAP_codex_t *codex,
                         GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke,cnt;
  char index;
  int eos=EOS;

  for (cnt=0; cnt<6; cnt++) SOAP_substitute_argum(argum, cnt, codex);

  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%c%n",
         &is,&js,&ks,&ie,&je,&ke,&index,&eos)!=7 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in ReverseIndex().");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);
  if (index!='i' && index!='j' && index!='k') SOAP_fatal_error(codex,"In ReverseIndex(), the index name must be either 'i', 'j', or 'k'.");
  ReverseIndex3D(gl3d, xyzgrid, codex, is, js, ks, ie, je, ke, index);
}


void Scale3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           double x_orig, double y_orig, double z_orig, double scalefactx,
	   double scalefacty, double scalefactz){
  long i,j,k;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Scale command failed.",i,j,k);
        }
        xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=x_orig+(xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x-x_orig)*scalefactx;
        xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=y_orig+(xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y-y_orig)*scalefacty;
        xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=z_orig+(xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z-z_orig)*scalefactz;
      }
    }
  }
}


static void Scale_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke;
  double x_orig,y_orig,z_orig,scalefactx,scalefacty,scalefactz;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%lg,%lg,%lg,%lg,%lg,%lg%n",
             &is,&js,&ks,&ie,&je,&ke,&x_orig,&y_orig,&z_orig,&scalefactx,
	     &scalefacty,&scalefactz,&eos)!=12 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Scale().");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);
  if (scalefactx<=0.0 || scalefacty<=0.0 || scalefactz<=0.0) SOAP_fatal_error(codex,"In Scale(), the scale factors must all be positive.");
  Scale3D(gl3d, xyzgrid, codex, is, js, ks, ie, je, ke, x_orig, y_orig, z_orig, scalefactx,
          scalefacty,scalefactz);
}




typedef struct {
  double val1,val2;
  char *expr;
  char varchange;
  SOAP_codex_t codex;
} argequation3D_t;


static double Equation3D(void *argequation3D, double val3){
  double res,z,y,x;
  char ystr[40],xstr[40],zstr[40];
  char varchange;
  char *expr;
  SOAP_codex_t *codex;
  int eos=EOS;

  codex=&(((argequation3D_t *)argequation3D)->codex);
  varchange=((argequation3D_t *)argequation3D)->varchange;
  x=0.0e0; y=0.0e0; z=0.0e0;
  if (varchange=='x') {
    y=((argequation3D_t *)argequation3D)->val1;
    z=((argequation3D_t *)argequation3D)->val2;
    x=val3;
  }
  if (varchange=='y') {
    x=((argequation3D_t *)argequation3D)->val1;
    z=((argequation3D_t *)argequation3D)->val2;
    y=val3;
  }
  if (varchange=='z') {
    x=((argequation3D_t *)argequation3D)->val1;
    y=((argequation3D_t *)argequation3D)->val2;
    z=val3;
  }
  expr=(char *)malloc((3+(long)strlen(((argequation3D_t *)argequation3D)->expr))*sizeof(char));
  strcpy(expr,((argequation3D_t *)argequation3D)->expr);
  /* here, setup codex and solve expression */
  sprintf(zstr,"%15.15E",z);
  sprintf(ystr,"%15.15E",y);
  sprintf(xstr,"%15.15E",x);
  SOAP_add_to_vars(codex, "x", xstr);
  SOAP_add_to_vars(codex, "y", ystr);
  SOAP_add_to_vars(codex, "z", zstr);
  SOAP_substitute_expression(&expr, codex);
  /* done */
  if (sscanf(expr,"%lg%n",&res,&eos)!=1 || expr[eos]!=EOS){
    GRIDG_fatal_error("Problem evaluating expression within Equation() command.");
  }
  free(expr);
  return(res);
}


static void Equation_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke,i,j,k;
  char varchange;
  double drootinit;
  argequation3D_t argequation3D;
  char *expr;
  long IFLAG,cnt;

  for (cnt=0; cnt<6; cnt++) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%c,",
         &is,&js,&ks,&ie,&je,&ke,&varchange)!=7){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Equation().");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);
  if (varchange!='x' && varchange!='y' && varchange!='z') SOAP_fatal_error(codex,"In Equation(), the variable to be changed must be either 'x', 'y', or 'z'.");

  expr=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&expr,*argum,7);
  replace_equal_sign(&expr);
  argequation3D.expr=(char *)malloc(((long)strlen(expr)+30)*sizeof(char));
  strcpy(argequation3D.expr,expr);
  SOAP_copy_codex(codex,&(argequation3D.codex));

  drootinit=0.0;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        drootinit=max(drootinit,fabs(xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z));
        drootinit=max(drootinit,fabs(xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y));
        drootinit=max(drootinit,fabs(xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x));
      }
    }
  }
  drootinit/=1e7;
  drootinit=max(1.0e-20,drootinit);


  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Equation command failed.",i,j,k);
        }
        argequation3D.varchange=varchange;
        if (varchange=='x') {
          argequation3D.val1=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y;
          argequation3D.val2=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z;
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=
            EXM_find_root_Newton_Raphson(&Equation3D, &argequation3D, xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x,
                           drootinit, 1.0e-10, 1.0e-50, &IFLAG);
        }
        if (varchange=='y') {
          argequation3D.val1=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x;
          argequation3D.val2=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z;
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=
            EXM_find_root_Newton_Raphson(&Equation3D, &argequation3D, xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y,
                           drootinit, 1.0e-10, 1.0e-50, &IFLAG);
        }
        if (varchange=='z') {
          argequation3D.val1=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x;
          argequation3D.val2=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y;
          xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=
            EXM_find_root_Newton_Raphson(&Equation3D, &argequation3D, xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z,
                           drootinit, 1.0e-10, 1.0e-50, &IFLAG);
        }
        if (IFLAG!=1) SOAP_fatal_error(codex,"In Equation(), root could not be obtained with Newton-Raphson solver for node %ld,%ld,%ld.",i,j,k);
      }
    }
  }
  free(expr);
  codex->vars=argequation3D.codex.vars;
  free( argequation3D.expr ); 
  SOAP_free_codex_copy(&(argequation3D.codex));
}


static void Spline_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke,i,j,k;
  long cnt;
  long N,n;
  double *f,*b,*x;
  double thisx;
  char *expr;

  for (cnt=0; cnt<6; cnt++) SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,",
         &is,&js,&ks,&ie,&je,&ke)!=6){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Spline();");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);

  expr=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&expr,*argum,6);

  if (strcmp(expr,"y(x)")!=0 && strcmp(expr,"y(z)")!=0 &&
      strcmp(expr,"z(x)")!=0 && strcmp(expr,"z(y)")!=0 &&
      strcmp(expr,"x(y)")!=0 && strcmp(expr,"x(z)")!=0) SOAP_fatal_error(codex,"In Spline(), the 7th argument must be either x(y), x(z), y(x), y(z), z(x), z(y).");

  N=SOAP_number_argums(*argum)-7;
  if (mod(N,2)!=0) SOAP_fatal_error(codex,"Number of arguments within Spline must be 7 + some even number.");
  N=N/2;
  if (N<4) SOAP_fatal_error(codex,"Number of data points supplied within Spline must be at least 4.");

  x=(double *)malloc(N*sizeof(double));
  f=(double *)malloc(N*sizeof(double));
  b=(double *)malloc(N*sizeof(double));
  
  for (n=0; n<N; n++) {
    SOAP_substitute_argum(argum, 7+n*2, codex);    
    x[n]=SOAP_get_argum_double(codex, *argum, 7+n*2);
    SOAP_substitute_argum(argum, 7+n*2+1, codex);    
    f[n]=SOAP_get_argum_double(codex, *argum, 7+n*2+1);
  }
  /* check if data points are valid (x[n+1]>x[n]) */
  for (n=0; n<N-1; n++){
    if (x[n+1]<=x[n]) SOAP_fatal_error(codex, "Data points supplied to spline must be such that x[i+1]>x[i] if y(x) or such that y[i+1]>y[i] if x(y)."); 
  }

  EXM_find_spline(N, x, f, b);
  
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT)
          SOAP_fatal_error(codex,"Node (%ld,%ld) not yet initialized. Spline command failed.",i,j);
        switch (expr[2]){
          case 'x':
            thisx=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x;
          break;
          case 'y':
            thisx=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y;
          break;
          case 'z':
            thisx=xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z;
          break;
          default:
            thisx=0.0; //to avoid compiler warning
            GRIDG_fatal_error("expr[2] can not be set to %c in Spline_argum().",expr[2]);
        }
        if (thisx>=x[0] && thisx<=x[N-1]) { 
          switch (expr[0]){
            case 'x':
              xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x=EXM_f_from_spline(N, x, f, b, thisx);
            break;
            case 'y':
              xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y=EXM_f_from_spline(N, x, f, b, thisx);
            break;
            case 'z':
              xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z=EXM_f_from_spline(N, x, f, b, thisx);
            break;
            default:
              GRIDG_fatal_error("expr[0] can not be set to %c in Spline_argum().",expr[0]);
          }
        } else {
          SOAP_fatal_error(codex,"Point %E is outside of the supplied range %E to %E.",thisx,x[0],x[N-1]);
        }
      }
    }
  }
  free(expr);
  free(x);
  free(b);
  free(f);
}


void Copy3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           long id, long jd, long kd){
  long i,j,k;
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized. Copy command failed.",i,j,k);
        }
        xyzgrid[GRIDG_ai3(gl3d,id-is+i,jd-js+j,kd-ks+k)].x
           =xyzgrid[GRIDG_ai3(gl3d,i,j,k)].x;
        xyzgrid[GRIDG_ai3(gl3d,id-is+i,jd-js+j,kd-ks+k)].y
           =xyzgrid[GRIDG_ai3(gl3d,i,j,k)].y;
        xyzgrid[GRIDG_ai3(gl3d,id-is+i,jd-js+j,kd-ks+k)].z
           =xyzgrid[GRIDG_ai3(gl3d,i,j,k)].z;
        xyzgrid[GRIDG_ai3(gl3d,id-is+i,jd-js+j,kd-ks+k)].INIT=TRUE;
      }
    }
  }
}


static void Copy_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long is,js,ks,ie,je,ke,id,jd,kd;
  int eos=EOS;
  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld%n",
       &is,&js,&ks,&ie,&je,&ke,&id,&jd,&kd,&eos)!=9 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Copy().");
  }
  verify_zone_validity(is, js, ks, ie, je, ke, gl3d, codex);
  verify_node_validity(id, jd, kd, gl3d, codex);
  Copy3D(gl3d, xyzgrid, codex, is, js, ks, ie, je, ke, id, jd, kd);
}


static void actions(char *action, char **argum, SOAP_codex_t *codex){
  GRIDG_gl3d_t *gl3d;
  GRIDG_xyzgrid_t **xyzgrid;

  gl3d=((actionsarg3D_t *)codex->action_args)->gl3d;
  xyzgrid=((actionsarg3D_t *)codex->action_args)->xyzgrid;

  if (strcmp(action,"Size")==0) {
    Size_argum(argum,codex,gl3d,xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Point")==0) {
    Point_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Corners")==0) {
    Corners_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"JoinCorners")==0) {
    JoinCorners_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Join")==0) {
    Join_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"JoinFaces")==0) {
    JoinFaces_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Plane")==0) {
    Plane_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Translate")==0) {
    translate_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Mirror")==0) {
    Mirror_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"ReverseIndex")==0) {
    ReverseIndex_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Scale")==0) {
    Scale_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Rotate")==0) {
    Rotate_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Equation")==0) {
    Equation_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Copy")==0) {
    Copy_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(action,"Spline")==0) {
    Spline_argum(argum,codex,*gl3d,*xyzgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
}

static void functions(char *functionname, char **argum,
                      char **returnstr, SOAP_codex_t *codex){
  GRIDG_gl3d_t *gl3d;
  GRIDG_xyzgrid_t **xyzgrid;
  long i,j,k;
  int eos=EOS;

  gl3d=((actionsarg3D_t *)codex->action_args)->gl3d;
  xyzgrid=((actionsarg3D_t *)codex->action_args)->xyzgrid;

  if (strcmp(functionname,"_x")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%ld,%ld%n",&i,&j,&k,&eos)!=3 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in _x().");
    }
    verify_node_validity(i, j, k, *gl3d, codex);
    if (!(*xyzgrid)[GRIDG_ai3(*gl3d,i,j,k)].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized, _x command failed.",i,j,k);
    }

    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",(*xyzgrid)[GRIDG_ai3(*gl3d,i,j,k)].x);
  }
  if (strcmp(functionname,"_y")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%ld,%ld%n",&i,&j,&k,&eos)!=3 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in _y().");
    }
    verify_node_validity(i, j, k, *gl3d, codex);
    if (!(*xyzgrid)[GRIDG_ai3(*gl3d,i,j,k)].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized, _y command failed.",i,j,k);
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",(*xyzgrid)[GRIDG_ai3(*gl3d,i,j,k)].y);
  }
  if (strcmp(functionname,"_z")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%ld,%ld%n",&i,&j,&k,&eos)!=3 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in _z().");
    }
    verify_node_validity(i, j, k, *gl3d, codex);
    if (!(*xyzgrid)[GRIDG_ai3(*gl3d,i,j,k)].INIT){
      SOAP_fatal_error(codex,"Node (%ld,%ld,%ld) not yet initialized, _z command failed.",i,j,k);
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",(*xyzgrid)[GRIDG_ai3(*gl3d,i,j,k)].z);
  }
}


static long add_one(long i, long imax){
  if (i==imax) i--; else i++;
  return(i);
}

static double _x(GRIDG_xyzgrid_t grid, long dim){
  if (dim==0) return(grid.x);
  if (dim==1) return(grid.y);
  return(grid.z);
}

/* this function finds the x,y,z coordinates of "ghost nodes" located outside the grid
   limits by "mirroring" the inner within the grid with respect to the boundary plane */ 
void find_xyz_on_ghost_nodes(GRIDG_xyzgrid_t *grid, GRIDG_gl3d_t gl){
  long dim,i,j,k;
  EXM_vec3D_t pa,pb,pc,pp_o,pp_m;

  pa[0]=1.0e0;   pa[1]=0.0e0;   pa[2]=0.0e0;
  pb[0]=0.0e0;   pb[1]=1.0e0;   pb[2]=0.0e0;
  pc[0]=0.0e0;   pc[1]=0.0e0;   pc[2]=1.0e0;
  pp_o[0]=0.0e0; pp_o[1]=0.0e0; pp_o[2]=0.0e0;

  /* here the idea is to find a plane in 3D represented
     by 3 points which lie on the boundary surface:
     pa,pb and pc. */
  /* do k=aks and k=ake */  /* grid[GRIDG_ai3(gl,i,j,k)].x */
  for (i=gl.is; i<=gl.ie; i++){
    for (j=gl.js; j<=gl.je; j++){
      /* do k=ake */
      for (dim=0; dim<3; dim++){
        pa[dim]=_x(grid[GRIDG_ai3(gl,add_one(i,gl.ie),j,gl.ke)],dim);
        pb[dim]=_x(grid[GRIDG_ai3(gl,i,add_one(j,gl.je),gl.ke)],dim);
        pc[dim]=_x(grid[GRIDG_ai3(gl,i,j,gl.ke)],dim);
        pp_o[dim]=_x(grid[GRIDG_ai3(gl,i,j,gl.ke-1)],dim);
      }
      EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
      grid[GRIDG_ai3(gl,i,j,gl.ke+1)].x=pp_m[0];
      grid[GRIDG_ai3(gl,i,j,gl.ke+1)].y=pp_m[1];
      grid[GRIDG_ai3(gl,i,j,gl.ke+1)].z=pp_m[2];
      /* do k=aks */
      for (dim=0; dim<3; dim++){
        pa[dim]=_x(grid[GRIDG_ai3(gl,add_one(i,gl.ie),j,gl.ks)],dim);
        pb[dim]=_x(grid[GRIDG_ai3(gl,i,add_one(j,gl.je),gl.ks)],dim);
        pc[dim]=_x(grid[GRIDG_ai3(gl,i,j,gl.ks)],dim);
        pp_o[dim]=_x(grid[GRIDG_ai3(gl,i,j,gl.ks+1)],dim);
      }
      EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
      grid[GRIDG_ai3(gl,i,j,gl.ks-1)].x=pp_m[0];
      grid[GRIDG_ai3(gl,i,j,gl.ks-1)].y=pp_m[1];
      grid[GRIDG_ai3(gl,i,j,gl.ks-1)].z=pp_m[2];
    }
  }

  /* do j=ajs a3 j=aje */
  for(i=gl.is; i<=gl.ie; i++){
    for(k=gl.ks-1; k<=gl.ke+1; k++){
      /* do j=aje */
      for (dim=0; dim<3; dim++){
        pa[dim]=_x(grid[GRIDG_ai3(gl,add_one(i,gl.ie),gl.je,k)],dim);
        pb[dim]=_x(grid[GRIDG_ai3(gl,i,gl.je,k)],dim);
        pc[dim]=_x(grid[GRIDG_ai3(gl,i,gl.je,add_one(k,gl.ke+1))],dim);
        pp_o[dim]=_x(grid[GRIDG_ai3(gl,i,gl.je-1,k)],dim);
      }
      EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
      grid[GRIDG_ai3(gl,i,gl.je+1,k)].x=pp_m[0];
      grid[GRIDG_ai3(gl,i,gl.je+1,k)].y=pp_m[1];
      grid[GRIDG_ai3(gl,i,gl.je+1,k)].z=pp_m[2];
      /* do j=ajs */
      for (dim=0; dim<3; dim++){
        pa[dim]=_x(grid[GRIDG_ai3(gl,add_one(i,gl.ie),gl.js,k)],dim);
        pb[dim]=_x(grid[GRIDG_ai3(gl,i,gl.js,k)],dim);
        pc[dim]=_x(grid[GRIDG_ai3(gl,i,gl.js,add_one(k,gl.ke+1))],dim);
        pp_o[dim]=_x(grid[GRIDG_ai3(gl,i,gl.js+1,k)],dim);
      }
      EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
      grid[GRIDG_ai3(gl,i,gl.js-1,k)].x=pp_m[0];
      grid[GRIDG_ai3(gl,i,gl.js-1,k)].y=pp_m[1];
      grid[GRIDG_ai3(gl,i,gl.js-1,k)].z=pp_m[2];
    }
  }

  /* do i=GRIDG_ai3s and i=GRIDG_ai3e */
  /* copy x,y,z of is plane to is-1 plane and of ie plane to ie+1 plane*/
  for(j=gl.js-1; j<=gl.je+1; j++){
    for(k=gl.ks-1; k<=gl.ke+1; k++){
      /* do i=GRIDG_ai3e */
      for (dim=0; dim<3; dim++){
        pa[dim]=_x(grid[GRIDG_ai3(gl,gl.ie,j,k)],dim);
        pb[dim]=_x(grid[GRIDG_ai3(gl,gl.ie,add_one(j,gl.je+1),k)],dim);
        pc[dim]=_x(grid[GRIDG_ai3(gl,gl.ie,j,add_one(k,gl.ke+1))],dim);
        pp_o[dim]=_x(grid[GRIDG_ai3(gl,gl.ie-1,j,k)],dim);
      }
      EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
      grid[GRIDG_ai3(gl,gl.ie+1,j,k)].x=pp_m[0];
      grid[GRIDG_ai3(gl,gl.ie+1,j,k)].y=pp_m[1];
      grid[GRIDG_ai3(gl,gl.ie+1,j,k)].z=pp_m[2];
      /* do i=GRIDG_ai3s */
      for (dim=0; dim<3; dim++){
        pa[dim]=_x(grid[GRIDG_ai3(gl,gl.is,j,k)],dim);
        pb[dim]=_x(grid[GRIDG_ai3(gl,gl.is,add_one(j,gl.je+1),k)],dim);
        pc[dim]=_x(grid[GRIDG_ai3(gl,gl.is,j,add_one(k,gl.ke+1))],dim);
        pp_o[dim]=_x(grid[GRIDG_ai3(gl,gl.is+1,j,k)],dim);
      }
      EXM_mirror_point_wrt_plane(pa, pb, pc, pp_o, pp_m);
      grid[GRIDG_ai3(gl,gl.is-1,j,k)].x=pp_m[0];
      grid[GRIDG_ai3(gl,gl.is-1,j,k)].y=pp_m[1];
      grid[GRIDG_ai3(gl,gl.is-1,j,k)].z=pp_m[2];
    }
  }
}



void RXGrid3D(char **argum, void *action_args, SOAP_codex_t *codex){
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


static void verify_all_nodes_initialized(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid){
  long i,j,k;

  for(i=gl3d.is; i<=gl3d.ie; i++){
    for(j=gl3d.js; j<=gl3d.je; j++){
      for(k=gl3d.ks; k<=gl3d.ke; k++){
        if (!xyzgrid[GRIDG_ai3(gl3d,i,j,k)].INIT){
          GRIDG_fatal_error("The x,y,z values were not provided to node (%ld,%ld,%ld). Grid generation failed.",i,j,k);
        }
      }
    }
  }
}


void GRIDG_read_grid_3D_from_file(char *filename, GRIDG_gl3d_t *gl3d, GRIDG_xyzgrid_t **xyzgrid,
                bool VERBOSE, bool *Problem){
  actionsarg3D_t actionsarg3D;

  actionsarg3D.gl3d=gl3d;
  actionsarg3D.xyzgrid=xyzgrid;

  read_RX_grid(filename, "Grid", &RXGrid3D, &GRIDG_write_grid_3D_to_file,
             &actionsarg3D, VERBOSE);
  verify_all_nodes_initialized(*gl3d, *xyzgrid);  
  find_xyz_on_ghost_nodes(*xyzgrid, *gl3d);
}


void GRIDG_read_grid_3D_from_argum(char *argum, SOAP_codex_t *codex, GRIDG_gl3d_t *gl3d, GRIDG_xyzgrid_t **xyzgrid){
  actionsarg3D_t actionsarg3D;

  actionsarg3D.gl3d=gl3d;
  actionsarg3D.xyzgrid=xyzgrid;

  codex->ACTION=TRUE;
  codex->action=&actions;
  codex->action_args=&actionsarg3D;
  codex->FUNCTION=TRUE;
  codex->function=&functions;
  codex->function_args=&actionsarg3D;
  SOAP_process_code(argum, codex, SOAP_VARS_KEEP_ALL);
  verify_all_nodes_initialized(*gl3d, *xyzgrid);  
  find_xyz_on_ghost_nodes(*xyzgrid, *gl3d);
}
