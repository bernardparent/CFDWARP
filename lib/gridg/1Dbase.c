// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2000 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "share.h"
#include "1Dbase.h"

typedef struct {
  GRIDG_gl1d_t *gl1d;
  GRIDG_xgrid_t **xgrid;
} actionsarg1D_t;



static void verify_zone_validity(long is, long ie, GRIDG_gl1d_t gl1d, SOAP_codex_t *codex){
  if (ie<is) 
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be less than is (%ld).",ie,is);
  if (ie>gl1d.ie)
    SOAP_fatal_error(codex,"Zone boundaries invalid: ie (%ld) can not be greater than domain.ie (%ld) as specified in Size().",ie,gl1d.ie);
  if (is<gl1d.is)
    SOAP_fatal_error(codex,"Zone boundaries invalid: is (%ld) can not be less than domain.is (%ld) as specified in Size().",is,gl1d.is);
}


static void verify_node_validity(long i,  GRIDG_gl1d_t gl1d, SOAP_codex_t *codex){
  if (i>gl1d.ie)
    SOAP_fatal_error(codex,"Zone boundaries invalid: i (%ld) can not be greater than domain.ie (%ld) as specified in Size().",i,gl1d.ie);
  if (i<gl1d.is)
    SOAP_fatal_error(codex,"Zone boundaries invalid: i (%ld) can not be less than domain.is (%ld) as specified in Size().",i,gl1d.is);
}


void GRIDG_write_grid_1D_to_file(FILE **controlfile){
  fprintf(*controlfile," \n");
  fprintf(*controlfile," \n");
  fprintf(*controlfile,"Grid(\n");
  fprintf(*controlfile,"  Size(1, 30);\n");
  fprintf(*controlfile,"  Point(1,  0.0e0,1.0e0);\n");
  fprintf(*controlfile,"  Point(10, 1.0e0,1.0e0);\n");
  fprintf(*controlfile,"  Point(20, 2.0e0,1.2e0);\n");
  fprintf(*controlfile,"  Point(30, 2.5e0,2.0e0);\n");
  fprintf(*controlfile,"  JoinCorners(1, 10, EE,0.5e0,1.0e0,1.0e0);\n");
  fprintf(*controlfile,"  JoinCorners(20, 30, EE,0.5e0,1.0e0,1.0e0);\n");
  fprintf(*controlfile,"  Join(10, 20, EE, 0.5e0, 1.0e0, 1.0e0);\n");
  fprintf(*controlfile,");\n");
}



static double find_distance(GRIDG_xgrid_t *xgrid, GRIDG_gl1d_t gl1d,
                    long i1, long i2){
  double dist;
  dist=fabs(xgrid[GRIDG_ai1(gl1d,i2)].x
           -xgrid[GRIDG_ai1(gl1d,i1)].x);
  return(dist);
}

void Size1D(GRIDG_gl1d_t *gl, GRIDG_xgrid_t **xgrid,
            long is, long ie){
  long i;

  gl->is=is;
  gl->ie=ie;

  *xgrid = (GRIDG_xgrid_t *) malloc((gl->ie-gl->is+1+4)*sizeof(GRIDG_xgrid_t));

  for (i=is; i<=ie; i++){
    (*xgrid)[GRIDG_ai1(*gl,i)].x=(double)(i-is)/(double)(ie-is);
  }
}


static void Size_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl1d_t *gl1d, GRIDG_xgrid_t **xgrid){
  long is,ie;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);

  if (sscanf(*argum,"%ld,%ld%n",&is,&ie,&eos)!=2 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Size().");
  }
  if (is>=ie) SOAP_fatal_error(codex,"In Size(), is (%ld) must be less than ie (%ld).",is,ie);
  Size1D(gl1d, xgrid, is, ie);
}


void Point1D(GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid,
             long i, double x, double Acs){
  xgrid[GRIDG_ai1(gl1d,i)].x=x;
  xgrid[GRIDG_ai1(gl1d,i)].Acs=Acs;
}


static void Point_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid){
  long i;
  double x,Acs;
  int eos=EOS;

  SOAP_substitute_all_argums(argum, codex);
  if (sscanf(*argum,"%ld,%lg,%lg%n",&i,&x,&Acs,&eos)!=3 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Point().");
  }
  verify_node_validity(i,  gl1d, codex);
  Point1D(gl1d, xgrid, i, x, Acs);
}


static int find_x_along_i_local(GRIDG_xgrid_t *xgrid, double segval1, double segval2,
                     double NLoverN, char modeL, char modeH,
                     GRIDG_gl1d_t gl, long i1, long i2, SOAP_codex_t *codex,
                     int(*find_segment)(double *, long , double,
                           double , double , double, SOAP_codex_t * )){
  int err;
  if (modeL=='G' || modeL=='g') segval1=find_distance(xgrid, gl, i1,
                    i1+sgndif(i1,i2));
  if (modeH=='G' || modeH=='g') segval2=find_distance(xgrid, gl, i2,
                    i2-sgndif(i1,i2));

  err=find_x_along_i(xgrid, segval1, segval2, NLoverN,
                  gl, i1, i2, codex, find_segment);
  return(err);
}


int JoinCorners1D(GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid, SOAP_codex_t *codex,
            long i1, long i2,
            char modeL_i, char modeH_i, double NLoverN_i, double valL_i, double valH_i){
  int (*FindSegment_i)(double *, long , double, double , double , double, SOAP_codex_t * );
  int err;

  find_segment_type(&FindSegment_i, modeL_i, modeH_i, codex);

  err=find_x_along_i_local(xgrid, valL_i, valH_i,
                        NLoverN_i, modeL_i, modeH_i, gl1d,
                        i1, i2, codex, FindSegment_i);
  return(err);
}


static void JoinCorners_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid){
  long i1,i2,cnt;
  char modeL_i,modeH_i;
  double NLoverN_i,valL_i,valH_i;
  int err;
  int eos=EOS;

  for (cnt=0; cnt<6; cnt++) if (cnt!=2)
    SOAP_substitute_argum(argum, cnt, codex);
  if (sscanf(*argum,"%ld,%ld,%c%c,%lg,%lg,%lg%n",
       &i1,&i2,&modeL_i,&modeH_i,&NLoverN_i,&valL_i,&valH_i,&eos)!=7 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in JoinCorners().");
  }
  verify_zone_validity(i1, i2, gl1d, codex);
  err=JoinCorners1D(gl1d, xgrid, codex, i1, i2, modeL_i, modeH_i, NLoverN_i, valL_i, valH_i);
  if (err!=0) SOAP_fatal_error(codex,"Error occurred while executing JoinCorners() command.");
}


int Join1D(GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid, SOAP_codex_t *codex,
                 long i1, long i2, char modeL, char modeH,
                 double NLoverN, double valL, double valH){
  int (*FindSegment)(double *, long , double, double , double , double, SOAP_codex_t * );
  int err;
  find_segment_type(&FindSegment, modeL, modeH, codex);

  err=find_x_along_i_local(xgrid, valL, valH,
                        NLoverN, modeL, modeH, gl1d,
                        i1, i2, codex, FindSegment);
  return(err);
}


static void Join_argum(char **argum, SOAP_codex_t *codex, GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid){
  long i1,i2,cnt;
  char modeL,modeH;
  double NLoverN,valL,valH;
  int err;
  int eos=EOS;

  for (cnt=0; cnt<6; cnt++) if (cnt!=2) {
    SOAP_substitute_argum(argum, cnt, codex);
  }
  if (sscanf(*argum,"%ld,%ld,%c%c,%lg,%lg,%lg%n",
         &i1,&i2,&modeL,&modeH,&NLoverN,&valL,&valH,&eos)!=7 || (*argum)[eos]!=EOS){
    SOAP_fatal_error(codex,"One or more argument(s) could not be read properly in Join().");
  }
  verify_zone_validity(i1, i2, gl1d, codex);
  err=Join1D(gl1d, xgrid, codex, i1, i2, modeL, modeH, NLoverN, valL, valH);
  if (err!=0) SOAP_fatal_error(codex,"Error occurred while executing Join() command.");
}




static void actions(char *actionname, char **argum, SOAP_codex_t *codex){
  GRIDG_gl1d_t *gl1d;
  GRIDG_xgrid_t **xgrid;

  gl1d=((actionsarg1D_t *)codex->action_args)->gl1d;
  xgrid=((actionsarg1D_t *)codex->action_args)->xgrid;
  if (strcmp(actionname,"Size")==0) {
    Size_argum(argum,codex,gl1d,xgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Point")==0) {
    Point_argum(argum,codex,*gl1d,*xgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"JoinCorners")==0) {
    JoinCorners_argum(argum,codex,*gl1d,*xgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Join")==0) {
    Join_argum(argum,codex,*gl1d,*xgrid);
    codex->ACTIONPROCESSED=TRUE;
  }
}

static void functions(char *functionname, char **argum,
                      char **returnstr, SOAP_codex_t *codex){
  GRIDG_gl1d_t *gl1d;
  GRIDG_xgrid_t **xgrid;
  long i;
  int eos=EOS;

  gl1d=((actionsarg1D_t *)codex->action_args)->gl1d;
  xgrid=((actionsarg1D_t *)codex->action_args)->xgrid;

  if (strcmp(functionname,"_x")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld%n",&i,&eos)!=1 || (*argum)[eos]!=EOS)
       SOAP_fatal_error(codex,"One argument could not be read properly in _x().");
    verify_node_validity(i,  *gl1d, codex);
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",(*xgrid)[GRIDG_ai1(*gl1d,i)].x);
  }
}


void RXGrid1D(char **argum, void *action_args, SOAP_codex_t *codex){

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


static void find_x_bdry(GRIDG_xgrid_t *grid, GRIDG_gl1d_t gl){

  grid[GRIDG_ai1(gl,gl.is-1)].x=grid[GRIDG_ai1(gl,gl.is)].x-
     (grid[GRIDG_ai1(gl,gl.is+1)].x-grid[GRIDG_ai1(gl,gl.is)].x);
  grid[GRIDG_ai1(gl,gl.is-1)].Acs=grid[GRIDG_ai1(gl,gl.is)].Acs-
     (grid[GRIDG_ai1(gl,gl.is+1)].Acs-grid[GRIDG_ai1(gl,gl.is)].Acs);
  grid[GRIDG_ai1(gl,gl.ie+1)].x=grid[GRIDG_ai1(gl,gl.ie)].x+
     (grid[GRIDG_ai1(gl,gl.ie)].x-grid[GRIDG_ai1(gl,gl.ie-1)].x);
  grid[GRIDG_ai1(gl,gl.ie+1)].Acs=grid[GRIDG_ai1(gl,gl.ie)].Acs+
     (grid[GRIDG_ai1(gl,gl.ie)].Acs-grid[GRIDG_ai1(gl,gl.ie-1)].Acs);

}




void GRIDG_read_grid_1D_from_file(char *filename, GRIDG_gl1d_t *gl1d, GRIDG_xgrid_t **xgrid,
                bool VERBOSE, bool *Problem){
  actionsarg1D_t actionsarg1D;
  long cnt;

  actionsarg1D.gl1d=gl1d;
  actionsarg1D.xgrid=xgrid;

  read_RX_grid(filename, "Grid", &RXGrid1D, &GRIDG_write_grid_1D_to_file,
             &actionsarg1D, VERBOSE);
  find_x_bdry(*xgrid, *gl1d);
  if (VERBOSE) {
    for (cnt=gl1d->is; cnt<=gl1d->ie; cnt++){
      printf("cnt=%ld  x=%E  Acs=%E\n",cnt,(*xgrid)[GRIDG_ai1(*gl1d,cnt)].x,
              (*xgrid)[GRIDG_ai1(*gl1d,cnt)].Acs);
    }
  }
}


void GRIDG_read_grid_1D_from_argum(char *argum, SOAP_codex_t *codex, GRIDG_gl1d_t *gl1d, GRIDG_xgrid_t **xgrid){
  actionsarg1D_t actionsarg1D;

  actionsarg1D.gl1d=gl1d;
  actionsarg1D.xgrid=xgrid;

  codex->ACTION=TRUE;
  codex->action=&actions;
  codex->action_args=&actionsarg1D;
  codex->FUNCTION=TRUE;
  codex->function=&functions;
  codex->function_args=&actionsarg1D;
  SOAP_process_code(argum, codex, SOAP_VARS_KEEP_ALL);
  find_x_bdry(*xgrid, *gl1d);
}
