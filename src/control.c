// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1999-2015 Bernard Parent

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <src/control.h>
#include <src/init.h>
#include <src/bdry.h>
#include <src/data.h>
#include <src/post.h>
#include <model/metrics/_metrics.h>
#include <model/thermo/_thermo.h>
#include <model/_model.h>
#include <src/version.hh>

#define EOS 0
#define EOL '\n'
#define lengthcol1 35

#ifdef _2D
  typedef GRIDG_gl2d_t glgrid_t;
#endif
#ifdef _3D
  typedef GRIDG_gl3d_t glgrid_t;
#endif

#define numbib 40
#define maxstrlen 4096

typedef char longstr_t[maxstrlen];

const static longstr_t bibname[numbib]=
  {
   "parent2002a",
   "parent2011a",
   "parent2012a",
   "parent2013a",
   "parent2014a",
   "parent2015a",
   "parent2015b",
   "parent2016a",
   "mcbride2002a",
   "maccormack2001a",
   "bardina1994a", 
   "briley1980a", 
   "maccormack1997a",
   "beam1976a", 
   "steger1981a",
   "anderson1986a",
   "vanleer1974a",
   "roe1981a",
   "yee1990a", 
   "parent2016b",
   "parent2018a",
   "jachimowsky1988a",
   "yungster1994a",
   "kundu1999a",
   "jiang1996a",   
   "liu1994a",
   "bardina1987a",
   "peaceman1955a",
   "balsara2000a",
   "dunn1973a",
   "balsara2009a",
   "balsara2016a",
   "levy2000a",
   "dumbser2007a",
   "gottlieb1998a",
   "young1950a",
   "frankel1950a",
   "gnoffo2004a",
   "parent2017a",
   "parent2019a"
  };

const static longstr_t bib[numbib]=
  {//parent2002a
   "B Parent and JP Sislian. \"The Use of Domain Decomposition in Accelerating the Convergence of Quasihyperbolic Systems\", Journal of Computational Physics 179:140-169, 2002.",

   //parent2011a
   "B Parent, MN Shneider, and SO Macheret. \"Generalized Ohm's Law and Potential Equation in Computational Weakly-Ionized Plasmadynamics\", Journal of Computational Physics 230:1439-1453, 2011.",

   //parent2012a
   "B Parent. \"Positivity-Preserving High-Resolution Schemes for Systems of Conservation Laws\", Journal of Computational Physics 231:173-189, 2012.",

   //parent2013a
   "B Parent. \"Positivity-Preserving Flux Difference Splitting Schemes\", Journal of Computational Physics 243:194-209, 2013.",

   //parent2014a
   "B Parent, SO Macheret, and MN Shneider. \"Electron and Ion Transport Equations in Computational Weakly-Ionized Plasmadynamics\", Journal of Computational Physics 259:51-69, 2014.",

   //parent2015a
   "B Parent. \"Multidimensional Flux Difference Splitting Schemes\", AIAA Journal 53:1936-1948, 2015.",

   //parent2015b
   "B Parent, SO Macheret, and MN Shneider. \"Modeling Weakly-Ionized Plasmas in Magnetic Field: A New Computationally-Efficient Approach\", Journal of Computational Physics 300:779-799, 2015.",

   //parent2016a
   "B Parent, MN Shneider, and SO Macheret. \"Detailed Modeling of Plasmas for Computational Aerodynamics\", AIAA Journal, 54:898-911, 2016.",

   //mcbride2002a
   "BJ McBride, MJ Zehe, and S Gordon. \"NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species\", NASA TP 2002-211556, 2002.",

   //maccormack2001a
   "RW MacCormack. \"Iterative Modified Approximate Factorization\", Computers and Fluids 30:917-925, 2001.",

   //bardina1994a
   "JE Bardina, \"Three-Dimensional Navier-Stokes Method with Two-Equation Turbulence Models for Efficient Numerical Simulation of Hypersonic Flows\", AIAA-94-2950, 1994.",

   //briley1980a
   "WR Briley and H McDonald. \"On the Structure and Use of Linearized Block Implicit Schemes\", Journal of Computational Physics 34:54-73, 1980.",

   //maccormack1997a
   "RW MacCormack. \"A New Implicit Algorithm for Fluid Flow\", 13th Computational Fluid Dynamics Conference, AIAA Paper 97-2100, 1997.",

   //beam1976a
   "RM Beam and RF Warming. \"An Implicit Finite-Difference Algorithm for Hyperbolic Systems in Conservation-Law-Form\", Journal of Computational Physics 22:87-110, 1976.",

   //steger1981a
   "JL Steger and RF Warming. \"Flux Vector Splitting of the Inviscid Gasdynamic Equations with Application to Finite-Difference Methods\", Journal of Computational Physics 40:263-293, 1981.",

   //anderson1986a
   "WK Anderson, JL Thomas, and B Van Leer. \"Comparison of Finite Volume Flux Vector Splittings for the Euler Equations\", AIAA Journal 24:1453-1460, 1986.",

   //vanleer1974a
   "B Van Leer. \"Towards the Ultimate Conservation Scheme II. Monotonicity and Conservation Combined in a Second-Order Scheme\", Journal of Computational Physics 14:361-370, 1974.",

   //roe1981a
   "PL Roe. \"Approximate Riemann Solvers, Parameter Vectors, and Difference Schemes\", Journal of Computational Physics 43:357-372, 1981.",

   //yee1990a
   "HC Yee, GH Klopfer, and JL Montagne. \"High-Resolution Shock-Capturing Schemes for Inviscid and Viscous Hypersonic Flows\", Journal of Computational Physics 88:31-61, 1990.",

   //parent2016b
   "B Parent. \"Multidimensional High Resolution Schemes for Viscous Hypersonic Flows\", AIAA Journal, 55:141-152, 2017.",

   //parent2018a
   "B Parent. \"Positivity-Preserving Dual Time Stepping Schemes for Gas Dynamics\", Journal of Computational Physics 361:391-411, 2018.",

   //jachimowsky1988a
   "CJ Jachimowsky. \"An Analytical Study of the Hydrogen-Air Reaction Mechanism With Application To Scramjet Combustion\", NASA TP-2791, 1988.",

   //yungster1994a
   "S Yungster and MJ Rabinowitz. \"Computation of Shock-Induced Combustion Using a Detailed Methane-Air Mechanism\", Journal of Propulsion and Power 10:609-617, 1994.",

   //kundu1999a
   "KP Kundu, PF Penko, and TJ VanOverbeke. \"A Practical Mechanism for Computing Combustion in Gas Turbine Engines\", 35th Joint Propulsion Conference and Exhibit, AIAA Paper 99-2218, 1999.",

   //jiang1996a
   "G Jiang and CW Shu. \"Efficient Implementation of Weighted ENO Schemes\", Journal of Computational Physics 126:202-228, 1996.",

   //liu1994a
   "XD Liu, S Osher, and T Chan. \"Weighted Essentially Non-Oscillatory Schemes\", Journal of Computational Physics, 115:200-212, 1994.",

   //bardina1987a
   "J Bardina and CK Lombard. \"Three Dimensional Hypersonic Flow Simulations with the CSCM Implicit Upwind Navier-Stokes Method\", Proceedings of the 8th Computational Fluid Dynamics Conference, AIAA Paper 87-1114, 1987.",

   //peaceman155a
   "DW Peaceman and HH Rachford. \"The Numerical Solution of Parabolic and Elliptic Differential Equations\", J. Soc. Ind. Appl. Math. 3:28-41, 1955.",

   //balsara2000a
   "DS Balsara and CW Shu. \"Monotonicity Preserving Weighted Essentially Non-Oscillatory Schemes with Increasingly High Order of Accuracy\", Journal of Computational Physics 160:405-452, 2000.",

   //dunn1973a
   "MG Dunn and SW Kang. \"Theoretical and Experimental Studies of Reentry Plasmas\", NASA CR-2232, 1973.",

   //balsara2009a
   "DS Balsara, T Rumpf, M Dumbser, CD Munz. \"Efficient, High Accuracy ADER-WENO Schemes for Hydrodynamics and Divergence-free Magnetohydrodynamics\", Journal of Computational Physics 228:2480-2516, 2009.",

   //balsara2016a
   "DS Balsara, S Garain, and CW Shu. \"An Efficient Class of WENO schemes with Adaptive Order\", Journal of Computational Physics 326:780-804, 2016.",

   //levy2000a
   "D Levy, G Puppo, and G Russo. \"Compact Central WENO Schemes for Multidimensional Conservation Laws\", SIAM Journal of Scientific Computing 22:656-672, 2000.",

   //dumbser2007a
   "M Dumbser and M Kaser. \"Arbitrary High Order Non-oscillatory Finite Volume Schemes on Unstructured Meshes for Linear Hyperbolic Systems\", Journal of Computational Physics, 221:693-723, 2007.",

   //gottlieb1998a
   "S Gottlieb and CW Shu. \"Total Variation Diminishing Runge Kutta Schemes\", Mathematics of Computation, 67:73-85, 1998.",

   //young1950a
   "DM Young. \"Iterative Methods for Solving Partial Difference Equations of Elliptic Type\", PhD thesis, Harvard University, May 1950.",

   //frankel1950a
   "S Frankel, \"Convergence Rates of Iterative Treatments of Partial Differential Equations\", Mathematics of Computation 4:56-75, 1950.",

   //gnoffo2004a
   "P Gnoffo and  JA White. \"Computational Aerothermodynamic Simulation Issues on Unstructured Grids\", 37th AIAA Thermophysics Conference, AIAA Paper 2004-2371, 2004.",

   //parent2017a
   "B Parent. \"Multidimensional High-Resolution Schemes for Viscous Hypersonic Flows\", AIAA Journal, 55:141-152, 2017.",

   //parent2019a
   "B Parent. \"Making a Flux Positivity-Preserving: A General Purpose Filter for the Euler Equations\", AIAA Paper 2019-0906, AIAA Scitech, San Diego CA, 2019."

  };





char *replace_cites(char *str, int *cite){
  static char  newstr[maxstrlen];
  static char  repl[maxstrlen];
  static char  orig[maxstrlen];
  static char  temp[maxstrlen];
  long cntbib,cntcite,numcite;

  strcpy(newstr,str);

  for (cntbib=0; cntbib<numbib; cntbib++){
    if (strstr(newstr,bibname[cntbib])) {
      strcpy(orig,bibname[cntbib]);
      if (cite[cntbib]!=0) {
        numcite=cite[cntbib];
      } else {
        numcite=0;
        for (cntcite=0; cntcite<numbib; cntcite++){
          numcite=max(numcite,cite[cntcite]);
        }
        numcite++;
        cite[cntbib]=numcite;
      }
      sprintf(repl,"%ld",numcite);
      strcpy(temp,strrep(newstr,orig,repl)); 
      strcpy(newstr,temp);
    }
  }
  return(newstr);
}



char *center_string(char *str, int width, int indent){
  long cnt,spaces;
  static char  newstr[maxstrlen];
  strcpy(newstr,str);
  spaces=round(((width-indent)-strlen(str))/2.0)+indent;
  for (cnt=0; cnt<spaces; cnt++) strins(" ",newstr,0);
  return newstr;
}


void  write_hline(FILE *outputfile, long width, long indent ){
  long cnt;
  for (cnt=0; cnt<indent; cnt++) wfprintf(outputfile," ");
  for (cnt=indent; cnt<width; cnt++)
    wfprintf(outputfile,"_");
  wfprintf(outputfile,"\n\n");
}


void write_citations(FILE *outputfile, int *cite, int linewidth){
  long cnt,cnt2;
  static char  str[maxstrlen];

  for (cnt=1; cnt<=numbib; cnt++){
    for (cnt2=0; cnt2<numbib; cnt2++){
      if (cite[cnt2]==cnt) {
        if (cnt<10) {
          sprintf(str,"  [%ld]  %s\n",cnt,bib[cnt2]);
        } else {
          sprintf(str,"  [%ld] %s\n",cnt,bib[cnt2]);
        }
        wfprintf(outputfile,"%s",strwrpind(str, linewidth,-7));
      }
    }
  }
}


void write_modules_row(FILE *outputfile, char *col1, char *col2, int linewidth){
  static char  linestr[maxstrlen];


  strcpy(linestr,"  ");
  strcat(linestr,col1);
  while (strlen(linestr)<lengthcol1) {
    strcat(linestr," ");
  }
  strcat(linestr,col2);  
  wfprintf(outputfile,strwrpind(linestr,linewidth,-lengthcol1));
  wfprintf(outputfile,"\n");
}


void write_license ( FILE * outputfile ) {
  char linestr[10000];
  int linewidth,term_width,term_height;
  find_terminal_window_size(&term_width,&term_height);
  linewidth=min(maxlinewidth,max(40,term_width-2));
  linewidth=min(linewidth,60);
  if (outputfile!=stdout) linewidth=maxlinewidth;

  wfprintf(outputfile,"\n");
  sprintf(linestr,"  Copyright (c) 1998-2019, Bernard Parent\n  All rights reserved.\n\n  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:\n");
  wfprintf(outputfile,"%s\n",strwrpind(linestr, linewidth, -2));

  sprintf(linestr,"  1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.");
  wfprintf(outputfile,"%s\n\n",strwrpind(linestr, linewidth, -5));

  sprintf(linestr,"  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.");
  wfprintf(outputfile,"%s\n\n",strwrpind(linestr, linewidth, -5));

  sprintf(linestr,"  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.");
  wfprintf(outputfile,"%s\n\n",strwrpind(linestr, linewidth, -2));

  sprintf(linestr,"  The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of the CFDWARP project.");
  wfprintf(outputfile,"%s\n\n",strwrpind(linestr, linewidth, -2));

}


void write_modules(FILE *outputfile){
#ifdef _FLUID_MULTISPECIES
  long spec;
  char *speciesname=(char *)malloc(sizeof(char));
#endif
  int cite[numbib];
  long cnt;
  int linewidth,term_width,term_height;
  static char  linestr[maxstrlen];
#ifdef _USERNAME
  static char  str[maxstrlen];
#endif
  find_terminal_window_size(&term_width,&term_height);
  linewidth=min(maxlinewidth,max(50,term_width-2));
  if (outputfile!=stdout) linewidth=maxlinewidth;
  for (cnt=0; cnt<numbib; cnt++) cite[cnt]=0;

  write_hline(outputfile, linewidth, 2);
  wfprintf(outputfile,"%s\n",center_string("C F D W A R P",linewidth,1));
/*  wfprintf(outputfile,"\n");
  sprintf(linestr,"v. %s",VERSION);
  wfprintf(outputfile,"%s\n",center_string(linestr,linewidth,1));
*/
  wfprintf(outputfile,"\n");
  wfprintf(outputfile,"%s\n",center_string("A CFD Code for Plasma & Reactive Flow",linewidth,1));
  write_hline(outputfile, linewidth, 2);
  write_modules_row(outputfile,"Version",VERSION,linewidth);


  linestr[0]=EOS;
#ifdef __GNUC__
  #ifdef __INTEL_COMPILER
    strcat(linestr,"icc ");
    strcat(linestr,__VERSION__);
  #else
    strcat(linestr,"gcc ");
    strcat(linestr,__VERSION__);
  #endif
#else
  if (__STDC__) {
    strcat(linestr," using ANSI C compiler");
  }
#endif
#ifdef POSIXTHREADS
  strcat(linestr," with POSIX loop threads");
#endif
#ifdef ZONETHREADS
  strcat(linestr," with POSIX zone threads");
#endif
#ifdef OPENMPTHREADS
  strcat(linestr," with OPENMP threads");
#endif
#ifdef DISTMPI
  strcat(linestr," with distributed memory MPI");
#endif
#ifndef NDEBUG
  strcat(linestr," with assertions");
#endif
#ifdef M32
  strcat(linestr," with 32-bit-architecture compatibility");
#endif 

  write_modules_row(outputfile,"Compiler",linestr,linewidth);

  sprintf(linestr,"%s at %s",__DATE__,__TIME__);
#ifdef _USERNAME
  sprintf(str,  " by %s on %s",_USERNAME,_HOSTNAME);
  strcat(linestr,str);
#endif
  write_modules_row(outputfile,"Compiled on",linestr,linewidth);

  write_modules_row(outputfile,"Cycle Strategy",replace_cites(_CYCLE_METHOD,cite),linewidth);
  write_modules_row(outputfile,"Fluid Relaxation",replace_cites(_TS_METHOD,cite),linewidth);
  write_modules_row(outputfile,"EM Field Relaxation",replace_cites(_TSEMF_METHOD,cite),linewidth);
  write_modules_row(outputfile,"Convection Discretization",replace_cites(_RESCONV_METHOD,cite),linewidth);
  write_modules_row(outputfile,"Temporal Discretization",replace_cites(_RESTIME_METHOD,cite),linewidth);
  write_modules_row(outputfile,"Fluid Model",replace_cites(_FLUID_METHOD,cite),linewidth);
#ifdef _FLUID_MULTISPECIES
  write_modules_row(outputfile,"Thermodynamic Model",replace_cites(_THERMO_METHOD,cite),linewidth);
#else
  write_modules_row(outputfile,"Thermodynamic Model","Perfect Gas",linewidth);
#endif
  write_modules_row(outputfile,"Chemical Model",replace_cites(_CHEM_METHOD,cite),linewidth);
  write_modules_row(outputfile,"EM Field Model",replace_cites(_EMFIELD_METHOD,cite),linewidth);
  write_modules_row(outputfile,"Beam Model",replace_cites(_BEAM_METHOD,cite),linewidth);
  sprintf(linestr,"%d",nd);
  write_modules_row(outputfile,"Number of Dimensions",linestr,linewidth);
  sprintf(linestr,"%d",nf);
  write_modules_row(outputfile,"Number of Flux Components",linestr,linewidth);
#ifdef _FLUID_MULTISPECIES
  sprintf(linestr,"%d",ns);
  write_modules_row(outputfile,"Number of Species",linestr,linewidth);
  sprintf(linestr,"%d",ncs);
  write_modules_row(outputfile,"Number of Charged Species",linestr,linewidth);

  strcpy(linestr,"");  
  for (spec=0; spec<ns; spec++) {
    find_species_name(spec, &speciesname);
    if (spec>0) strcat(linestr,", ");
    strcat(linestr,speciesname);
  }
  free(speciesname);
  write_modules_row(outputfile,"List of Species",linestr,linewidth);
#endif
  write_hline(outputfile, linewidth, 2);
  write_citations(outputfile,cite,linewidth);

  wfprintf(outputfile,"\n");

  if (outputfile==stdout) {
    sprintf(linestr,"  CFDWARP is a Copyright (c) 1998-2019 of Bernard Parent. Use 'warp -l' for licensing terms.");
    wfprintf(outputfile,"%s\n",strwrpind(linestr, linewidth, -2));
    wfprintf(outputfile,"\n");
    wfprintf(outputfile,"%s\n","  Use 'warp -h' to list command line options.");
    wfprintf(outputfile,"\n");
  } else {

  }
}


void write_control(char *filename){
  FILE *controlfile;
#ifdef DISTMPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0) {
#endif
    controlfile = fopen(filename, "w");
    wfprintf(stdout,"Writing to control file %s...",filename);
    wfprintf(controlfile,"{\n");
    write_modules(controlfile);
    wfprintf(controlfile,"\n}\n\n");
#ifdef _3D
    GRIDG_write_grid_3D_to_file(&controlfile);
#endif
#ifdef _2D
    GRIDG_write_grid_2D_to_file(&controlfile);
#endif
    write_bdry_template(&controlfile);
    write_model_template(&controlfile);
    write_init_template(&controlfile);
    write_disc_template(&controlfile);
    write_cycle_template(&controlfile);
    write_post_template(&controlfile);
    fclose(controlfile);
    wfprintf(stdout,"done;\n");
#ifdef DISTMPI
  }
#endif
}


double _x_DISTMPI_global(np_t *np, gl_t *gl, SOAP_codex_t *codex, long i, long j, long k, long dim){
  double x;
#ifdef DISTMPI
  int rank;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==_node_rank(gl,i,j,k)) {
#endif
    if (i<gl->domain.is || i>gl->domain.ie
#ifdef _2DL
        || j<gl->domain.js || j>gl->domain.je
#endif
#ifdef _3DL
        || k<gl->domain.ks || k>gl->domain.ke
#endif
       ) SOAP_fatal_error(codex,"Function _%c() has been called with i=%ld,j=%ld,k=%ld that is out of bounds.",(int)'x'+dim,i,j,k);
    x=_x(np[_ai(gl,i,j,k)],dim);
#ifdef DISTMPI
  }
  MPI_Bcast(&x,1,MPI_DOUBLE,_node_rank(gl,i,j,k),MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  return(x);
}


double _x_DISTMPI_local(np_t *np, gl_t *gl, SOAP_codex_t *codex, long i, long j, long k, long dim){
  double x;
  zone_t domain_local;
  find_zone_intersection(gl->domain_lim,gl->domain_all,&domain_local);

  if (i<domain_local.is || i>domain_local.ie
#ifdef _2DL
      || j<domain_local.js || j>domain_local.je
#endif
#ifdef _3DL
      || k<domain_local.ks || k>domain_local.ke
#endif
     ) SOAP_fatal_error(codex,"Function _%c() has been called with i=%ld,j=%ld,k=%ld that is out of bounds. "
                   "Make sure that domain.is<=i<=domain.ie, domain.js<=j<=domain.je, domain.ks<=k<=domain.ke "
                   "where 'domain' refers to the zone accessible by one MPI process. "
                   "Although seldom necessary, to access x,y,z outside of the domain, use the "
                   "_x_global(), _y_global, _z_global() functions [DON'T USE THE *_global() FUNCTIONS UNLESS ABSOLUTELY "
                   "NECESSARY AS THIS WILL SLOW DOWN THE CODE WHEN USING MPI]." ,(int)'x'+dim,i,j,k);
  x=_x(np[_ai(gl,i,j,k)],dim);
  
  return(x);
}


double _Omega_DISTMPI_global(np_t *np, gl_t *gl, long i, long j, long k){
  double Omega;
#ifdef DISTMPI
  int rank;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==_node_rank(gl,i,j,k)) {
#endif
    if (i<gl->domain.is || i>gl->domain.ie
#ifdef _2DL
        || j<gl->domain.js || j>gl->domain.je
#endif
#ifdef _3DL
        || k<gl->domain.ks || k>gl->domain.ke
#endif
       ) fatal_error("Function _Omega() has been called with i,j,k that is out of bounds.");
    Omega=_Omega(np[_ai(gl,i,j,k)],gl);
#ifdef DISTMPI
  }
  MPI_Bcast(&Omega,1,MPI_DOUBLE,_node_rank(gl,i,j,k),MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  return(Omega);
}


void read_control_functions(char *functionname, char **argum,
                            char **returnstr, SOAP_codex_t *codex){
  np_t **np;
  gl_t *gl;
  long i,j,k;
  int eos=EOS;
  np=((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;
#ifdef _2D
  k=0;
  if (strcmp(functionname,"_x")==0 || strcmp(functionname,"_y")==0 
   || strcmp(functionname,"_x_global")==0 || strcmp(functionname,"_y_global")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%ld%n",&i,&j,&eos)!=2 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"Problem reading arguments in function %s().",functionname);
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
  }
  if (strcmp(functionname,"_x")==0) sprintf(*returnstr,"%E",_x_DISTMPI_local(*np,gl,codex,i,j,k,0));
  if (strcmp(functionname,"_y")==0) sprintf(*returnstr,"%E",_x_DISTMPI_local(*np,gl,codex,i,j,k,1));
  if (strcmp(functionname,"_x_global")==0) sprintf(*returnstr,"%E",_x_DISTMPI_global(*np,gl,codex,i,j,k,0));
  if (strcmp(functionname,"_y_global")==0) sprintf(*returnstr,"%E",_x_DISTMPI_global(*np,gl,codex,i,j,k,1));
#endif
#ifdef _3D
  if (strcmp(functionname,"_x")==0 || strcmp(functionname,"_y")==0 || strcmp(functionname,"_z")==0
   || strcmp(functionname,"_x_global")==0 || strcmp(functionname,"_y_global")==0 || strcmp(functionname,"_z_global")==0) {
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%ld,%ld%n",&i,&j,&k,&eos)!=3 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"Problem reading arguments in function %s().",functionname);
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
  }
  if (strcmp(functionname,"_x")==0) sprintf(*returnstr,"%E",_x_DISTMPI_local(*np,gl,codex,i,j,k,0));
  if (strcmp(functionname,"_y")==0) sprintf(*returnstr,"%E",_x_DISTMPI_local(*np,gl,codex,i,j,k,1));
  if (strcmp(functionname,"_z")==0) sprintf(*returnstr,"%E",_x_DISTMPI_local(*np,gl,codex,i,j,k,2));
  if (strcmp(functionname,"_x_global")==0) sprintf(*returnstr,"%E",_x_DISTMPI_global(*np,gl,codex,i,j,k,0));
  if (strcmp(functionname,"_y_global")==0) sprintf(*returnstr,"%E",_x_DISTMPI_global(*np,gl,codex,i,j,k,1));
  if (strcmp(functionname,"_z_global")==0) sprintf(*returnstr,"%E",_x_DISTMPI_global(*np,gl,codex,i,j,k,2));
#endif

  if (strcmp(functionname,"_Omega")==0) {
    if (!gl->METRICS_INITIALIZED){
      SOAP_fatal_error(codex,"Function _Omega() can not be called prior to the initialization of the metrics at the end of the Bdry() module.");
    }
    SOAP_substitute_all_argums(argum,codex);
#ifdef _2D
    if (sscanf(*argum, "%ld,%ld%n",&i,&j,&eos)!=2 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"Problem reading arguments in function _Omega().");
    }
#endif
#ifdef _3D
    if (sscanf(*argum, "%ld,%ld,%ld%n",&i,&j,&k,&eos)!=3 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"Problem reading arguments in function _Omega().");
    }
#endif
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_Omega_DISTMPI_global(*np,gl,i,j,k));
  }
}


#ifdef DISTMPI
static void broadcast_soap_vars(SOAP_codex_t *codex){
  long cnt,numvar;
  long value_length,name_length;
  char *name;
  char *value;
  int rank;

  name=(char *)malloc(sizeof(char));
  value=(char *)malloc(sizeof(char));
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (numvar=0; codex->vars[numvar].name!=NULL; numvar++);
  MPI_Bcast(&numvar,1,MPI_LONG,0,MPI_COMM_WORLD);
  for (cnt=0; cnt<numvar; cnt++) {
    if (rank==0) {
      value_length=strlen(codex->vars[cnt].value)+1;
      name_length=strlen(codex->vars[cnt].name)+1;
    }
    MPI_Bcast(&value_length,1,MPI_LONG,0,MPI_COMM_WORLD);
    MPI_Bcast(&name_length,1,MPI_LONG,0,MPI_COMM_WORLD);

    name=(char *)realloc(name,(name_length+3)*sizeof(char));
    value=(char *)realloc(value,(value_length+3)*sizeof(char));
    if (rank==0) {
      strcpy(name,codex->vars[cnt].name);
      strcpy(value,codex->vars[cnt].value);
    }
    MPI_Bcast(value,value_length,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(name,name_length,MPI_CHAR,0,MPI_COMM_WORLD);

    SOAP_add_to_vars(codex, name, value);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

/* the following first sets the offset to 0, then 1, then -1 */
static long _node_offset_from_cnt(long cnt){
  long offset;
  offset=0;
  if (cnt==0) offset=0;
  if (cnt==1) offset=1;
  if (cnt==2) offset=-1;
  return(offset);
}


void find_metrics_on_all_nodes(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k,TYPELEVEL,BDRYMETRICS;
  long l_C,l_B,l_A,l,theta,thetasgn,dim;
#ifdef _2DL
  long offset1,offset2,cnt1,cnt2;
#endif
#ifdef _3D
  long offset3,cnt3;
#endif
  bool UPDATED;

#ifdef EMFIELD
  TYPELEVEL=TYPELEVEL_EMFIELD;
#else
  TYPELEVEL=TYPELEVEL_FLUID;
#endif 

  /* first do the inner nodes */
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef EMFIELD
        // here  check if EMFIELD metrics is compatible with fluid metrics 
        if (is_node_valid(np[l],TYPELEVEL_FLUID) && !is_node_valid(np[l],TYPELEVEL_EMFIELD)){
          fatal_error("Emfield node must be valid where fluid node is valid at i=%ld j=%ld k=%ld.",_i(l,gl,0),_i(l,gl,1),_i(l,gl,2));
        }
#endif
        if (is_node_valid(np[l],TYPELEVEL) && !is_node_bdry(np[l],TYPELEVEL)){
          update_metrics_at_node(np,gl,l);
        }
        for (dim=0; dim<nd; dim++){
#ifdef EMFIELD
          // here  check if EMFIELD metrics is compatible with fluid metrics 
          if (is_node_valid(np[_al(gl,l,dim,+1)],TYPELEVEL_FLUID) && !is_node_valid(np[_al(gl,l,dim,+1)],TYPELEVEL_EMFIELD)){
            fatal_error("Emfield node must be valid where fluid node is valid at i=%ld j=%ld k=%ld.",_i(_al(gl,l,dim,+1),gl,0),_i(_al(gl,l,dim,+1),gl,1),_i(_al(gl,l,dim,+1),gl,2));
          }
#endif
        }
      end3DL
    end2DL
  end1DL


  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        for (dim=0; dim<nd; dim++){
          l=_ai(gl,i,j,k);
          if (is_node_valid(np[l],TYPELEVEL) && is_node_valid(np[_al(gl,l,dim,+1)],TYPELEVEL)){
            update_metrics_at_interface(np,gl,l,_al(gl,l,dim,+1),dim);
          }
        }
      end3DL
    end2DL
  end1DL


  /* then do the boundary nodes */
  /* note: we should use BDRYMETRICS_CENTERED instead of BDRYMETRICS_NORMAL on symmetry planes
     for best results: need to improve this */
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);

        if (is_node_bdry(np[l],TYPELEVEL)) {
          UPDATED=FALSE;
          if (find_bdry_direc(np, gl, l, TYPELEVEL, &theta, &thetasgn)){
            l_A=l;
            l_B=_al(gl,l,theta,thetasgn);
            l_C=_al(gl,l,theta,2*thetasgn);
            if (is_node_bdry_symmetry_plane(np[l_A])){
              BDRYMETRICS=BDRYMETRICS_CENTERED;
            } else {
              BDRYMETRICS=BDRYMETRICS_NORMAL;
            }
            update_metrics_at_bdry_node(np, gl, TYPELEVEL, l_A, l_B, l_C, BDRYMETRICS);
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
                if  (  is_node_inner(np[l_B],TYPELEVEL)
                  && is_node_inner(np[l_C],TYPELEVEL) && !UPDATED){
                  if (is_node_bdry_symmetry_plane(np[l_A])){
                    BDRYMETRICS=BDRYMETRICS_CENTERED;
                  } else {
                    BDRYMETRICS=BDRYMETRICS_NORMAL;
                  }
                  update_metrics_at_bdry_node(np, gl, TYPELEVEL, l_A, l_B, l_C,
                       BDRYMETRICS);
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
                  if  (  is_node_inner(np[l_B],TYPELEVEL)
                       &&is_node_inner(np[l_C],TYPELEVEL) && !UPDATED){
                    if (is_node_bdry_symmetry_plane(np[l_A])){
                      BDRYMETRICS=BDRYMETRICS_CENTERED;
                    } else {
                      BDRYMETRICS=BDRYMETRICS_NORMAL;
                    }
                    update_metrics_at_bdry_node(np, gl, TYPELEVEL, l_A, l_B, l_C,
                         BDRYMETRICS);
                    UPDATED=TRUE;
                  }
                }
              }
            }
            #endif
          }
          if (!UPDATED) {
            
            display_node_type_window_local_process(stdout, np, gl, TYPELEVEL, _i(l,gl,0), _i(l,gl,1), _i(l,gl,2), 40);
            fatal_error("Problem updating boundary node in find_metrics_on_all_nodes function at i=%ld j=%ld k=%ld.",_i(l,gl,0),_i(l,gl,1),_i(l,gl,2));
          }
        }
      end3DL
    end2DL
  end1DL
  gl->METRICS_INITIALIZED=TRUE;
  //display_node_type_window(stdout, np, gl, TYPELEVEL, 1, 15, 1, 20);
}




void readcontrol_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  np_t **np;
  gl_t *gl;
  bool GRIDGENERATED;
  long cnt,cnt2,cnt_mem;
  zone_t domain_local;
  SOAP_codex_t codex_mem;
  input_t *input;
  long i,j,k;
#ifdef _3D
  GRIDG_xyzgrid_t *xgrid;
#endif
#ifdef _2D
  GRIDG_xygrid_t *xgrid;
#endif
  glgrid_t glgrid;
  double val[3];
#ifdef DISTMPI
  int rank,proc;
#endif



  cnt_mem=-1;
  do {
   cnt_mem++;
  } while ((codex->vars)[cnt_mem].name!=NULL);

  SOAP_copy_codex(codex,&codex_mem);

  np=((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;


  GRIDGENERATED=FALSE;
#ifdef DISTMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
#endif
#ifdef _2D
    if (strcmp(actionname,"Grid")==0) {
      ((readcontrolarg_t *)codex->action_args)->module_level++;
#ifdef DISTMPI
      if (rank==0) {
#endif
        wfprintf(stdout,"Grid..");
        GRIDG_read_grid_2D_from_argum(*argum, codex, &glgrid, &xgrid);
#ifdef DISTMPI
      }
//      wfprintf(stdout,"[broadcasting]");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.is,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.ie,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.js,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.je,1,MPI_LONG,0,MPI_COMM_WORLD);
#endif
      GRIDGENERATED=TRUE;
      codex->ACTIONPROCESSED=TRUE;
    }
#endif
#ifdef _3D
    if (strcmp(actionname,"Grid")==0) {
      ((readcontrolarg_t *)codex->action_args)->module_level++;
#ifdef DISTMPI
      if (rank==0) {
#endif
        wfprintf(stdout,"Grid..");
        GRIDG_read_grid_3D_from_argum(*argum, codex, &glgrid, &xgrid);
#ifdef DISTMPI
      }
//      wfprintf(stdout,"[broadcasting]");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.is,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.ie,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.js,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.je,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.ks,1,MPI_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(&glgrid.ke,1,MPI_LONG,0,MPI_COMM_WORLD);
#endif
      GRIDGENERATED=TRUE;
      codex->ACTIONPROCESSED=TRUE;
    }
#endif
  if (GRIDGENERATED){
    
#ifdef DISTMPI
    gl->DISTDOMAIN=TRUE;
    broadcast_soap_vars(codex);
#endif
    if (glgrid.is!=1) fatal_error("Size() within Grid() must be called with is=1.");
    gl->domain_all.is=glgrid.is;
    gl->domain_all.ie=glgrid.ie;
#ifdef _2DL
    if (glgrid.js!=1) fatal_error("Size() within Grid() must be called with js=1.");
    gl->domain_all.js=glgrid.js;
    gl->domain_all.je=glgrid.je;
#else
    gl->domain_all.js=1;
    gl->domain_all.je=1;
#endif
#ifdef _3DL
    if (glgrid.ks!=1) fatal_error("Size() within Grid() must be called with ks=1.");
    gl->domain_all.ks=glgrid.ks;
    gl->domain_all.ke=glgrid.ke;
#else
    gl->domain_all.ks=1;
    gl->domain_all.ke=1;    
#endif

    gl->domain=gl->domain_all;

#ifdef DISTMPI
    zone_t domainlim;
    int ndomi,ndomj,ndomk;
    double ndomi_opt,ndomj_opt,ndomk_opt,err,errthis;
    domainlim=_domain_lim_from_domain(gl->domain_all, gl);

    switch (gl->DISTDOMAIN_METHOD){
      case DISTDOMAIN_METHOD_AUTO:
#ifdef _2D
        //  proci*proj=proc;
        //  proci/procj=(domainlim.ie-domainlim.is+1)/(domainlim.je-domainlim.js+1);
        //  sqr(proci)=proc*(domainlim.ie-domainlim.is+1)/(domainlim.je-domainlim.js+1);
        ndomi_opt=min((double)proc,max(1.0,(sqrt((double)proc*(double)(domainlim.ie-domainlim.is+1)/(double)(domainlim.je-domainlim.js+1)))));
        ndomj_opt=(double)proc/ndomi_opt;
        ndomk_opt=1.0;
#endif
#ifdef _3D
        //  proci*proj*prock=proc;
        //  procj/proci=(domainlim.je-domainlim.js+1)/(domainlim.ie-domainlim.is+1);
        //  prock/proci=(domainlim.ke-domainlim.ks+1)/(domainlim.ie-domainlim.is+1);
        //  proci=pow(proc/((domainlim.je-domainlim.js+1)/(domainlim.ie-domainlim.is+1)*(domainlim.ke-domainlim.ks+1)/(domainlim.ie-domainlim.is+1)),1/3);
        ndomi_opt=pow((double)proc/((double)(domainlim.je-domainlim.js+1)/(double)(domainlim.ie-domainlim.is+1)*(double)(domainlim.ke-domainlim.ks+1)/(double)(domainlim.ie-domainlim.is+1)),1.0/3.0);
        ndomj_opt=ndomi_opt*(double)(domainlim.je-domainlim.js+1)/(double)(domainlim.ie-domainlim.is+1);
        ndomk_opt=(double)proc/ndomj_opt/ndomi_opt;
#endif

        gl->numdomain_i=1;
        gl->numdomain_j=1;
        gl->numdomain_k=1;
    
        //  wfprintf(stdout,"[%dx%d opt domains]",ndomi_opt,ndomj_opt);
        err=1e99;
        for1DL(ndomi,1,proc)
          for2DL(ndomj,1,proc)
            for3DL(ndomk,1,proc)
              errthis=fabs((double)ndomi-ndomi_opt)/ndomi_opt+fabs((double)ndomj-ndomj_opt)/ndomj_opt +fabs((double)ndomk-ndomk_opt)/ndomk_opt;
              if (errthis<err && ndomi*ndomj*ndomk==proc){
                gl->numdomain_i=ndomi;
                gl->numdomain_j=ndomj;
                gl->numdomain_k=ndomk;
                err=errthis;
              }
            end3DL
          end2DL
        end1DL
      break;
      case DISTDOMAIN_METHOD_USERSPEC:
#ifdef _2D
        if (gl->numdomain_i*gl->numdomain_j!=proc){
          fatal_error("The number of MPI domains specified is %ldx%ld=%ld. Such must be equal to the number of MPI processes proc=%ld.",gl->numdomain_i,gl->numdomain_j,gl->numdomain_i*gl->numdomain_j,proc);
        }
#endif
#ifdef _3D
        if (gl->numdomain_i*gl->numdomain_j*gl->numdomain_k!=proc){
          fatal_error("The number of MPI domains specified is %ldx%ldx%ld=%ld. Such must be equal to the number of MPI processes proc=%ld.",gl->numdomain_i,gl->numdomain_j,gl->numdomain_k,gl->numdomain_i*gl->numdomain_j*gl->numdomain_k,proc);
        }
#endif
      break;
      default:
        fatal_error("DISTDOMAIN_METHOD can not be set to %ld.",gl->DISTDOMAIN_METHOD);
    }


    wfprintf(stdout,"%dx%dx%d domains..",gl->numdomain_i,gl->numdomain_j,gl->numdomain_k);
    if (gl->numdomain_i*gl->numdomain_j*gl->numdomain_k!=proc){
      fatal_error("Number of MPI processes doesn't match grid size: numdomain_i=%ld numdomain_j=%ld numdomain_k=%ld. ",gl->numdomain_i,gl->numdomain_j,gl->numdomain_k);
    }
    
    gl->domain_from_rank=(zone_t *)malloc(sizeof(zone_t)*(2+proc));
    gl->domain_lim_from_rank=(zone_t *)malloc(sizeof(zone_t)*(2+proc));
    for (cnt=0; cnt<proc; cnt++) {
      gl->domain_from_rank[cnt]=_domain_from_rank_mem(cnt, gl);
      gl->domain_lim_from_rank[cnt]=_domain_lim_from_rank_mem(cnt, gl);
    }
    gl->domain=_domain_from_rank(rank, gl);

   
#endif



    init_data_structure_and_create_nodes(np, gl, gl->domain, gl->domain_all);
#ifdef DISTMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for1DL (i,gl->domain_all.is-1,gl->domain_all.ie+1)
      for2DL (j,gl->domain_all.js-1,gl->domain_all.je+1)
        for3DL (k,gl->domain_all.ks-1,gl->domain_all.ke+1)
#ifdef DISTMPI
          if (rank==0) {
#endif
#ifdef _3D
            val[0]=xgrid[GRIDG_ai3(glgrid,i+glgrid.is-gl->domain_all.is,j+glgrid.js-gl->domain_all.js,k+glgrid.ks-gl->domain_all.ks)].x;
            val[1]=xgrid[GRIDG_ai3(glgrid,i+glgrid.is-gl->domain_all.is,j+glgrid.js-gl->domain_all.js,k+glgrid.ks-gl->domain_all.ks)].y;
            val[2]=xgrid[GRIDG_ai3(glgrid,i+glgrid.is-gl->domain_all.is,j+glgrid.js-gl->domain_all.js,k+glgrid.ks-gl->domain_all.ks)].z;
#endif
#ifdef _2D
            val[0]=xgrid[GRIDG_ai2(glgrid,i+glgrid.is-gl->domain_all.is,j+glgrid.js-gl->domain_all.js)].x;
            val[1]=xgrid[GRIDG_ai2(glgrid,i+glgrid.is-gl->domain_all.is,j+glgrid.js-gl->domain_all.js)].y;
#endif
#ifdef DISTMPI
          }
          MPI_Bcast_Node(val,3,MPI_DOUBLE,0,MPI_COMM_WORLD,i,j,k,gl);
#endif

          if (is_node_in_zone(i,j,k,gl->domain_lim)) {
#ifdef _3D
          (*np)[_ai(gl,i,j,k)].bs->x[0]=val[0];
          (*np)[_ai(gl,i,j,k)].bs->x[1]=val[1];
          (*np)[_ai(gl,i,j,k)].bs->x[2]=val[2];
#endif
#ifdef _2D
          (*np)[_ai(gl,i,j,k)].bs->x[0]=val[0];
          (*np)[_ai(gl,i,j,k)].bs->x[1]=val[1];
#endif
          }
        end3DL
      end2DL
    end1DL
#ifdef DISTMPI
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) free(xgrid);
#else
    free(xgrid);
#endif

    wfprintf(stdout,"done;\n");
    gl->iter=1;
    gl->sigma1=0.5e0;
    gl->sigma2=0.5e0;
    gl->CFL=1.0e0;


    gl->effiter_U=0.0e0;
    gl->effiter_R=0.0e0;
#ifdef EMFIELD
    gl->Lc=1.0e0;
    gl->effiter_U_emfield=0.0e0;
    gl->effiter_R_emfield=0.0e0;
#endif

    find_zone_intersection(gl->domain_lim,gl->domain_all,&domain_local);
    add_int_to_codex(codex,"domain.is",   domain_local.is);
    add_int_to_codex(codex,"domain.ie",   domain_local.ie);
#ifdef _2DL
    add_int_to_codex(codex,"domain.js",   domain_local.js);
    add_int_to_codex(codex,"domain.je",   domain_local.je);
#endif
#ifdef _3D
    add_int_to_codex(codex,"domain.ks",   domain_local.ks);
    add_int_to_codex(codex,"domain.ke",   domain_local.ke);
#endif

  }

  if (strcmp(actionname,"Bdry")==0 && !((readcontrolarg_t *)codex->action_args)->GRIDONLY) {
    if (((readcontrolarg_t *)codex->action_args)->module_level!=1)
      fatal_error("The Grid%ldd() module was not found.",nd);
    wfprintf(stdout,"Bdry..");
    read_bdry(*argum, codex);
#ifdef DISTMPI
    find_metrics_on_all_nodes((*np), gl, _zone_intersection(gl->domain_all,_zone_expansion(gl->domain,+1)));
#else
    find_metrics_on_all_nodes((*np), gl, gl->domain);
#endif
#ifdef DISTMPI
    double metrics[nd*nd+1];
    long dim,dim2;
    int TYPELEVEL;
//    wfprintf(stdout,"[broadcasting]");
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef EMFIELD
    TYPELEVEL=TYPELEVEL_EMFIELD;
#else
    TYPELEVEL=TYPELEVEL_FLUID;
#endif
    for1DL (i,gl->domain_all.is,gl->domain_all.ie)
      for2DL (j,gl->domain_all.js,gl->domain_all.je)
        for3DL (k,gl->domain_all.ks,gl->domain_all.ke)
          if (rank==_node_rank(gl, i, j, k) && is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL)) {
            metrics[0]=(*np)[_ai(gl,i,j,k)].bs->Omega;
            for (dim=0; dim<nd; dim++) {
              for (dim2=0; dim2<nd; dim2++) {
                metrics[1+dim*nd+dim2]=(*np)[_ai(gl,i,j,k)].bs->X[dim][dim2];
              }
            }
          }
          MPI_Bcast_Node(metrics,nd*nd+1,MPI_DOUBLE,_node_rank(gl,i,j,k),MPI_COMM_WORLD,i,j,k,gl);
          if (is_node_in_zone(i,j,k,gl->domain_lim) && is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL)) {
            (*np)[_ai(gl,i,j,k)].bs->Omega=metrics[0];
            for (dim=0; dim<nd; dim++) {
              for (dim2=0; dim2<nd; dim2++) {
                (*np)[_ai(gl,i,j,k)].bs->X[dim][dim2]=metrics[1+dim*nd+dim2];
              }
            }              
          }
        end3DL
      end2DL
    end1DL
    MPI_Barrier(MPI_COMM_WORLD);

/*  if (is_node_in_zone(5,23,0,gl->domain_lim)){
    //printf("V=%E\n",_V((*np)[_ai(gl,5,24,0)],0));
    //printf("Vstar=%E\n",_Vstar((*np)[_ai(gl,5,24,0)],1));
    //printf("X=%E\n",_X((*np)[_ai(gl,5,23,0)],0,0)); 
  }*/

#endif

    wfprintf(stdout,"done;\n");
    ((readcontrolarg_t *)codex->action_args)->module_level++;
    input=((readcontrolarg_t *)codex->action_args)->input;
    for1DL(i,gl->domain_lim.is,gl->domain_lim.ie)
     for2DL(j,gl->domain_lim.js,gl->domain_lim.je)
      for3DL(k,gl->domain_lim.ks,gl->domain_lim.ke)
        (*np)[_ai(gl,i,j,k)].INIT_FLUID=FALSE;        
        #ifdef EMFIELD
          (*np)[_ai(gl,i,j,k)].INIT_EMFIELD=FALSE;        
        #endif
      end3DL
     end2DL
    end1DL

    if (!input->INTERPOLATION){
      read_data_file(*input, *np, gl);
#ifdef UNSTEADY
      add_double_to_codex(codex,"time",  gl->time);
#endif
    }
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,ControlActionName_Model)==0 && !((readcontrolarg_t *)codex->action_args)->GRIDONLY) {
    if (((readcontrolarg_t *)codex->action_args)->module_level!=2)
      fatal_error("The Bdry() module was not found.");
    wfprintf(stdout,"%s..",ControlActionName_Model);
    read_model(*argum, codex);
    wfprintf(stdout,"done;\n");
    ((readcontrolarg_t *)codex->action_args)->module_level++;
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Init")==0  && !((readcontrolarg_t *)codex->action_args)->GRIDONLY) {
    if (((readcontrolarg_t *)codex->action_args)->module_level!=3)
      fatal_error("The %s() module was not found.",ControlActionName_Model);
    input=((readcontrolarg_t *)codex->action_args)->input;
    if (input->INTERPOLATION){
      read_data_file(*input, *np, gl);
#ifdef UNSTEADY
      add_double_to_codex(codex,"time",  gl->time);
#endif
    }
    if (!((readcontrolarg_t *)codex->action_args)->input->READDATAFILE){
      wfprintf(stdout,"Init..");
      read_init(*argum, codex);
      wfprintf(stdout,"done;\n");
    }
    ((readcontrolarg_t *)codex->action_args)->module_level++;
    if (((readcontrolarg_t *)codex->action_args)->RESETITERCOUNT) {
       ((readcontrolarg_t *)codex->action_args)->gl->iter=1;
       ((readcontrolarg_t *)codex->action_args)->gl->effiter_U=0.0;
       ((readcontrolarg_t *)codex->action_args)->gl->effiter_R=0.0;
#ifdef EMFIELD
       ((readcontrolarg_t *)codex->action_args)->gl->effiter_U_emfield=0.0;
       ((readcontrolarg_t *)codex->action_args)->gl->effiter_R_emfield=0.0;
#endif
    }
    for1DL(i,gl->domain_lim.is,gl->domain_lim.ie)
     for2DL(j,gl->domain_lim.js,gl->domain_lim.je)
      for3DL(k,gl->domain_lim.ks,gl->domain_lim.ke)
        if (is_node_valid((*np)[_ai(gl,i,j,k)], TYPELEVEL_FLUID) &&  !(*np)[_ai(gl,i,j,k)].INIT_FLUID){
          fatal_error("The fluid properties at node %ld"if2DL(",%ld")if3DL(",%ld")" were not initialized properly.",i
#ifdef _2DL
  ,j
#endif
#ifdef _3DL
  ,k
#endif
          );
        }        
        #ifdef EMFIELD
        if (is_node_valid((*np)[_ai(gl,i,j,k)], TYPELEVEL_EMFIELD) &&  !(*np)[_ai(gl,i,j,k)].INIT_EMFIELD){
          fatal_error("The EM properties at node %ld"if2DL(",%ld")if3DL(",%ld")" were not initialized properly.",i
#ifdef _2DL
  ,j
#endif
#ifdef _3DL
  ,k
#endif
          );
        }
        #endif
      end3DL
     end2DL
    end1DL
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Disc")==0  && !((readcontrolarg_t *)codex->action_args)->GRIDONLY) {
    if (((readcontrolarg_t *)codex->action_args)->module_level!=4)
      fatal_error("The Init() module was not found.");
    wfprintf(stdout,"Disc..");
    read_disc(*argum, codex);
    wfprintf(stdout,"done;\n");
    ((readcontrolarg_t *)codex->action_args)->module_level++;
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Cycle")==0 && ((readcontrolarg_t *)codex->action_args)->CYCLEMODULE) {
    if (((readcontrolarg_t *)codex->action_args)->module_level!=5)
      fatal_error("The Disc() module was not found.");
    wfprintf(stdout,"Cycle..");
    read_cycle(*argum, codex);
    wfprintf(stdout,"done;\n");
    ((readcontrolarg_t *)codex->action_args)->module_level++;
    codex->ACTIONPROCESSED=TRUE;
  }
  if (strcmp(actionname,"Post")==0 && ((readcontrolarg_t *)codex->action_args)->POSTMODULE) {
    if (((readcontrolarg_t *)codex->action_args)->module_level!=6)
      fatal_error("The Cycle() module was not found.");
    wfprintf(stdout,"Post..");
    read_post(*argum, codex);
//    wfprintf(stdout,"done;\n");
    ((readcontrolarg_t *)codex->action_args)->module_level++;
    codex->ACTIONPROCESSED=TRUE;
  }
  codex->ACTION=codex_mem.ACTION;
  codex->action=codex_mem.action;
  codex->action_args=codex_mem.action_args;
  codex->FUNCTION=codex_mem.FUNCTION;
  codex->function=codex_mem.function;
  codex->function_args=codex_mem.function_args;

  cnt=-1;
  do {
   cnt++;
  } while ((codex->vars)[cnt].name!=NULL);

  for (cnt2=cnt_mem; cnt2<cnt; cnt2++) {
    free((codex->vars)[cnt2].name);
    free((codex->vars)[cnt2].value);
  }
  (codex->vars)[cnt_mem].name=NULL;
  SOAP_free_codex_copy(&codex_mem);
}


void read_control(char *control_filename, input_t input, bool CYCLEMODULE, bool POSTMODULE, bool GRIDONLY,
                 bool RESETITERCOUNT, np_t **np, gl_t *gl){
  SOAP_codex_t codex;
  long i,j,k;
  readcontrolarg_t readcontrolarg;
  char *code;
#ifdef DISTMPI
  long nn_sum;
  int rank;
#endif
  gl->CONTROL_READ=FALSE;
  gl->INIT_FLUID_READ=FALSE;
  gl->INIT_EMFIELD_READ=FALSE;
  gl->CYCLE_FLUID_READ=FALSE;
  gl->CYCLE_EMFIELD_READ=FALSE;
  gl->METRICS_INITIALIZED=FALSE;
  wfprintf(stdout,"Reading CFDWARP control file %s..\n",control_filename);
  SOAP_init_codex(&codex,control_filename);
#ifdef DISTMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank!=0) {
    codex.SCREENOUTPUT=FALSE;
    codex.FILEOUTPUT=FALSE;
    codex.SYSTEMCALL=FALSE;
  }
#endif
  readcontrolarg.np=np;
  readcontrolarg.gl=gl;
  readcontrolarg.VERBOSE=TRUE;
  readcontrolarg.POSTMODULE=POSTMODULE;
  readcontrolarg.GRIDONLY=GRIDONLY;
  readcontrolarg.CYCLEMODULE=CYCLEMODULE;
  readcontrolarg.RESETITERCOUNT=RESETITERCOUNT;
  readcontrolarg.input=&input;
  readcontrolarg.module_level=0;
  code=(char *)malloc(sizeof(char));
  SOAP_store_file_as_string(control_filename, &code);
  codex.ACTION=TRUE;
  codex.action=&readcontrol_actions;
  codex.action_args=(void *)&readcontrolarg;
  codex.FUNCTION=TRUE;
  codex.function=&read_control_functions;
  codex.function_args=(void *)&readcontrolarg;
  SOAP_insert_line_numbers_in_code(&code,1);
#ifdef UNSTEADY
  add_int_to_codex(&codex,"UNSTEADY",   TRUE);
  add_double_to_codex(&codex,"time",  gl->time);
#else
  add_int_to_codex(&codex,"UNSTEADY",   FALSE);
#endif
  add_string_to_codex(&codex,"domain.is",   "UNDEFINED");
  add_string_to_codex(&codex,"domain.ie",   "UNDEFINED");
#ifdef _2DL
  add_string_to_codex(&codex,"domain.js",   "UNDEFINED");
  add_string_to_codex(&codex,"domain.je",   "UNDEFINED");
#endif
#ifdef _3D
  add_string_to_codex(&codex,"domain.ks",   "UNDEFINED");
  add_string_to_codex(&codex,"domain.ke",   "UNDEFINED");
#endif

  thread_lock_init(&(gl->lock),THREADTYPE_ALL);
  gl->nn=LONG_MAX;
    /* we need to set gl->nn to a value prior to ProcessCode or gl->effiter_R will not be updated correctly in the function  update_residual() part of the file res.c when called by the cycle module PRIOR to the first iteration*/

  SOAP_process_code(code, &codex, SOAP_VARS_KEEP_ALL);
  if (readcontrolarg.module_level!=7 && POSTMODULE)
    fatal_error("The Post() module was not found.");
  if (readcontrolarg.module_level!=6 && CYCLEMODULE && !POSTMODULE)
    fatal_error("The Cycle() module was not found.");
  if (readcontrolarg.module_level!=5 && !CYCLEMODULE && !POSTMODULE && !GRIDONLY)
    fatal_error("The Disc() module was not found.");
//  if (readcontrolarg.module_level!=4 && !CYCLEMODULE  && !POSTMODULE)
//    fatal_error("The Init() module was not found.");
  if (!gl->CYCLE_FLUID_READ && CYCLEMODULE) 
    fatal_error("The fluid module %s() was not found within Cycle().",_FLUID_ACTIONNAME);
  if (!gl->CYCLE_EMFIELD_READ && CYCLEMODULE) 
    fatal_error("The emfield module %s() was not found within Cycle().",_EMFIELD_ACTIONNAME);
  free(code);
  wfprintf(stdout,"done.\n");

  /* now compute gl->nn correctly */
  gl->nn=0;
  for1DL(i,gl->domain.is,gl->domain.ie)
    for2DL(j,gl->domain.js,gl->domain.je)
      for3DL(k,gl->domain.ks,gl->domain.ke)
        if (is_node_inner((*np)[_ai(gl,i,j,k)],TYPELEVEL_FLUID)) gl->nn++;
      end3DL
    end2DL
  end1DL
#ifdef DISTMPI
  MPI_Allreduce(&gl->nn, &nn_sum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  gl->nn=nn_sum;
#endif
  gl->CONTROL_READ=TRUE;

  if ( POSTMODULE || GRIDONLY )
    SOAP_free_codex ( &codex );
  else
    SOAP_free_codex_copy ( &codex );
}

