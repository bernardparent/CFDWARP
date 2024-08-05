// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2010-2011 Bernard Parent

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

#include "fluid.h"
#include "fluid_init.h"
#include <model/_model.h>
#include <src/init.h>

#define INIT_TYPE1 1
#define INIT_TYPE2 2
#define INIT_TYPE3 3

typedef struct {
  double P,T;
  dim_t V;
} Init1_t;


void write_init_fluid_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n"
  ,_FLUID_ACTIONNAME);
  wfprintf(*controlfile,
  "    {\n"
  "    _______________________________________________________________________________________\n"
  "\n"
  "    Initial Condition Type       Parameters\n"
  "    _______________________________________________________________________________________\n"
  "\n"
  "    INIT_TYPE1                   V[1]..V[nd],  P,  T     \n"
  "    INIT_TYPE2                   M[1]..M[nd],  P,  T     \n"
  "    INIT_TYPE3                   V[1]..V[nd],  P,  rho   \n"
  "    _______________________________________________________________________________________\n"
  "    }\n"
  "    Mx=2;\n"
  if2DL(
  "    My=0;\n"
  )
  if3DL(
  "    Mz=0;\n"
  )
  "    P=10000; {Pa}\n"
  "    T=300; {K}\n"
  "    All(INIT_TYPE2,Mx"if2DL(",My")if3DL(",Mz")",P,T);\n"
  "    {\n"
  "    Bdry(BDRY_WALLTFIXED1, INIT_TYPE2, Mx"if2DL(",My")if3DL(",Mz")",P,T);\n"
  "    Region(is,"if2DL("js,")if3DL("ks,")"  ie,"if2DL("je,")if3DL("ke,")"  INIT_TYPE2, Mx"if2DL(",My")if3DL(",Mz")",P,T);\n"
  "    }\n"
  "  );\n"
  );

}


void add_init_types_fluid_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"INIT_TYPE1",   INIT_TYPE1);
  add_int_to_codex(codex,"INIT_TYPE2",   INIT_TYPE2);
  add_int_to_codex(codex,"INIT_TYPE3",   INIT_TYPE3);
}



void find_U_init_2(np_t *np, long l, gl_t *gl, Init1_t Init1){
  find_U_2(np, l, gl, Init1.V, Init1.P, Init1.T);
}


/* V[dim], P,T */
void init_node_1(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim;
  Init1_t Init1;
  
  verify_positivity_of_determinative_property(initvar,nd,nd+1);

  for (dim=0; dim<nd; dim++) Init1.V[dim]=initvar[dim];
  Init1.P=initvar[nd];
  Init1.T=initvar[nd+1];
  //printf("-- P=%E T=%E\n",Init1.P,Init1.T);
  find_U_init_2(np, l, gl, Init1);
}


/* M[dim], P,T */
void init_node_2(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim;
  double a;
  Init1_t Init1;

  verify_positivity_of_determinative_property(initvar,nd,nd+1);

  Init1.P=initvar[nd];
  Init1.T=initvar[nd+1];
  a=sqrt(gl->model.fluid.gamma*gl->model.fluid.R*Init1.T);
  for (dim=0; dim<nd; dim++) Init1.V[dim]=initvar[dim]*a;
  find_U_init_2(np, l, gl, Init1);
}


/* V[dim], P,rho */
void init_node_3(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim;
  double T,P,rho;
  Init1_t Init1;

  verify_positivity_of_determinative_property(initvar,nd,nd+1);

  for (dim=0; dim<nd; dim++) Init1.V[dim]=initvar[dim];
  P=initvar[nd];
  rho=initvar[nd+1];
  T=P/rho/gl->model.fluid.R;
  Init1.P=P;
  Init1.T=T;
  //printf("-- P=%E T=%E\n",Init1.P,Init1.T);
  find_U_init_2(np, l, gl, Init1);
}


void init_node_fluid(np_t *np, long l, gl_t *gl, long inittype, initvar_t initvar){
  switch (inittype) {
    case INIT_TYPE1: 
      init_node_1(np, l, gl, initvar); 
    break;
    case INIT_TYPE2: 
      init_node_2(np, l, gl, initvar); 
    break;
    case INIT_TYPE3: 
      init_node_3(np, l, gl, initvar); 
    break;
    default:
      fatal_error("Initial condition type %ld invalid.",inittype);
  }
}


void find_default_initvar_name(initvarname_t *initvar_name){
  long dim;
  initvarname_t string;

  for (dim=0; dim<nd; dim++){
    sprintf(string,"V[%ld]",dim);
    strcpy(initvar_name[dim],string);
  }
  strcpy(initvar_name[nd],"P");
  strcpy(initvar_name[nd+1],"T");
} 

/* find initvar (the type of initvar must match defaultinitvartypefluid) */
void find_default_initvar(np_t *np, gl_t *gl, long l, initvar_t initvar){
  long dim;
  for (dim=0; dim<nd; dim++) initvar[dim]=_V(np[l],dim);
  initvar[nd]=_P(np[l],gl);
  initvar[nd+1]=_T(np[l],gl);
}


