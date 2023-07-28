// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2005-2006 Bernard Parent

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
#include <model/thermo/_thermo.h>
#include <model/transport/_transport.h>
#include <model/_model.h>
#include <src/init.h>
#include <model/share/fluid_share.h>

#define INIT_TYPE1 1
#define INIT_TYPE2 2
#define INIT_TYPE3 3
#define INIT_TYPE4 4
#define INIT_TYPE5 5
#define INIT_TYPE6 6


typedef struct {
  double T,Tv,k,psi;
  dim_t V;
  spec_t rho;
} Init1_t;


typedef struct {
  double T,Tv,P,k,psi;
  dim_t V;
  spec_t w;
} Init2_t;



void write_init_fluid_template(FILE **controlfile){
  char *specname,*initstr,*wdefault,*specstr1,*specstr2;
  wdefault=(char *)malloc(sizeof(char));
  specname=(char *)malloc(sizeof(char));
  initstr=(char *)malloc(sizeof(char)*1000);
  specstr1=(char *)malloc(sizeof(char));
  specstr2=(char *)malloc(sizeof(char));

  find_init_mass_fraction_templates(&specstr1,&specstr2);

  
  strcpy(initstr,"INIT_TYPE2,Mx"if2DL(",My")if3DL(",Mz")",P,T");
  strcat(initstr,specstr2);
  strcat(initstr,",k,psi,Tv");
    
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
  "    INIT_TYPE1                   V[1]..V[nd],  T,  rho, w[1]..w[ns],     k, psi, Tv\n"
  "    INIT_TYPE2                   M[1]..M[nd],  P,  T,   w[1]..w[ns],     k, psi, Tv\n"
  "    INIT_TYPE3                   M[1]..M[nd],  Re, T,   w[1]..w[ns],     k, psi, Tv\n"
  "    INIT_TYPE4                   Mmag, angles, P,  T,   w[1]..w[ns],     k, psi, Tv\n"
  "    INIT_TYPE5                   V[1]..V[nd],  P,  T,   w[1]..w[ns],     k, psi, Tv\n"
  "    INIT_TYPE6                   V[1]..V[nd],  P,  T,   chi[1]..chi[ns], k, psi, Tv\n"
  "    _______________________________________________________________________________________\n"
  "\n"
  "    (a) In the freestream set psi to 110*q_infty for the k-omega models.\n"
  "    (b) In the freestream make sure that 1E-5*q_infty^2/ReL<k<0.1*q_infty^2/ReL.\n"
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
  "    k=1e-6; {J/kg}\n"
  "    psi=110*sqrt(Mx^2"if2DL("+My^2")if3DL("+Mz^2")")*sqrt(1.4*287*T); {1/s for TURBMODEL_KOMEGA*}\n"
  "    Tv=300; {K}\n"
  "%s"
  "    All(%s);\n"
  "    {\n"
  "    Bdry(BDRY_WALLTFIXED1, %s);\n"
  "    Region(is,"if2DL("js,")if3DL("ks,")"  ie,"if2DL("je,")if3DL("ke,")" %s);\n"
  "    }\n"
  "  );\n",specstr1,initstr,initstr,initstr
  );
  free(specname);
  free(wdefault);
  free(initstr);
  free(specstr1);
  free(specstr2);
}


void add_init_types_fluid_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"INIT_TYPE1",   INIT_TYPE1);
  add_int_to_codex(codex,"INIT_TYPE2",   INIT_TYPE2);
  add_int_to_codex(codex,"INIT_TYPE3",   INIT_TYPE3);
  add_int_to_codex(codex,"INIT_TYPE4",   INIT_TYPE4);
  add_int_to_codex(codex,"INIT_TYPE5",   INIT_TYPE5);
  add_int_to_codex(codex,"INIT_TYPE6",   INIT_TYPE6);
}


void find_U_init_1(np_t *np, long l, gl_t *gl, Init1_t Init1){
  find_U(np, l, gl, Init1.rho, Init1.V, Init1.T, Init1.k, Init1.psi, Init1.Tv);
}


void find_U_init_2(np_t *np, long l, gl_t *gl, Init2_t Init2){
  bool CLIPPED;
  reformat_w(gl, Init2.w, "",&CLIPPED);
  find_U_2(np, l, gl, Init2.w, Init2.V, Init2.P, Init2.T, Init2.k, Init2.psi, Init2.Tv);
}



/* V[dim], T,rho,w[spec],k,psi,Tv */
void init_node_1(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim,spec;
  Init1_t Init1;
  initvar_t values;
  reformat_initvar_species_fractions(gl, initvar, values,nd+2);
  verify_positivity_of_determinative_property(values,nd,nd+ns+4);

  for (dim=0; dim<nd; dim++) Init1.V[dim]=values[dim];
  Init1.T=values[nd];
  for (spec=0; spec<ns; spec++) Init1.rho[spec]=values[nd+2+spec]*values[nd+1];
  Init1.k=values[nd+2+ns];
  Init1.psi=values[nd+3+ns];
  Init1.Tv=values[nd+4+ns];
  find_U_init_1(np, l, gl, Init1);
}

/* M[dim], P,T,w[spec],k,psi,Tv */
void init_node_2(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim,spec;
  double a;
  Init2_t Init2;
  initvar_t values;
  reformat_initvar_species_fractions(gl, initvar, values,nd+2);
  verify_positivity_of_determinative_property(values,nd,nd+ns+4);

  Init2.P=values[nd];
  Init2.T=values[nd+1];
  for (spec=0; spec<ns; spec++) Init2.w[spec]=values[nd+2+spec];
  a=_a_from_w_T_equilibrium(Init2.w,Init2.T);
  for (dim=0; dim<nd; dim++) Init2.V[dim]=values[dim]*a;
  Init2.k=values[nd+2+ns];
  Init2.psi=values[nd+3+ns];
  Init2.Tv=values[nd+4+ns];
  find_U_init_2(np, l, gl, Init2);
}


/* M[dim], Re,T,w[spec],k,psi,Tv */
void init_node_3(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim,spec;
  double a,Re,eta,kappa,sum,rho;
  spec_t nuk,rhok;
  Init2_t Init2;
  initvar_t values;
  reformat_initvar_species_fractions(gl, initvar, values,nd+2);
  verify_positivity_of_determinative_property(values,nd,nd+ns+4);

  Init2.T=values[nd+1];
  Re=values[nd];
  for (spec=0; spec<ns; spec++) Init2.w[spec]=values[nd+2+spec];
  a=_a_from_w_T_equilibrium(Init2.w,Init2.T);
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    Init2.V[dim]=values[dim]*a;
    sum=sum+sqr(Init2.V[dim]);
  }
  rho=1.0e0;
  Init2.Tv=values[nd+4+ns];
  for (spec=0; spec<ns; spec++) rhok[spec]=Init2.w[spec]*rho;
  find_nuk_eta_kappa(rhok, Init2.T, _Te_from_T_Tv(gl,Init2.T,Init2.Tv), nuk, &eta, &kappa);
  rho=Re*eta/sqrt(sum);
  Init2.P=_P_from_w_rho_T(Init2.w,rho,Init2.T);
  Init2.k=values[nd+2+ns];
  Init2.psi=values[nd+3+ns];
  find_U_init_2(np, l, gl, Init2);
}


/* Mflight and angles in degrees (if necessary), P,T,w[spec],k,psi,Tv */
void init_node_4(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim,spec;
  double a,theta,M_flight,phi;
  Init2_t Init2;
  initvar_t values;
  dim_t Mk;
  reformat_initvar_species_fractions(gl, initvar, values,nd+2);
  verify_positivity_of_determinative_property(values,nd,nd+ns+4);

  Init2.P=values[nd];
  Init2.T=values[nd+1];
  theta=0.0e0;
  phi=0.0e0;  
#ifdef _2DL
  theta=values[1];
#endif
#ifdef _3DL
  phi=values[2];
#endif
  M_flight=values[0];
  Mk[0]=M_flight*cos(theta)*cos(phi);
#ifdef _2DL
  Mk[1]=M_flight*sin(theta)*cos(phi);
#endif
#ifdef _3DL
  Mk[2]=M_flight*sin(phi);
#endif
  for (spec=0; spec<ns; spec++) Init2.w[spec]=values[nd+2+spec];
  a=_a_from_w_T_equilibrium(Init2.w,Init2.T);
  for (dim=0; dim<nd; dim++) Init2.V[dim]=Mk[dim]*a;
  Init2.k=values[nd+2+ns];
  Init2.psi=values[nd+3+ns];
  Init2.Tv=values[nd+4+ns];
  find_U_init_2(np, l, gl, Init2);
}


/* V[dim], P,T,w[spec],k,psi,Tv */
void init_node_5(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim,spec;
  Init2_t Init2;
  initvar_t values;
  reformat_initvar_species_fractions(gl, initvar, values,nd+2);
  verify_positivity_of_determinative_property(values,nd,nd+ns+4);

  Init2.P=values[nd];
  Init2.T=values[nd+1];
  for (spec=0; spec<ns; spec++) Init2.w[spec]=values[nd+2+spec];
  for (dim=0; dim<nd; dim++) Init2.V[dim]=values[dim];
  Init2.k=values[nd+2+ns];
  Init2.psi=values[nd+3+ns];
  Init2.Tv=values[nd+4+ns];
  find_U_init_2(np, l, gl, Init2);
}

/* V[dim], P,T,chi[spec],k,psi,Tv */
void init_node_6(np_t *np, long l, gl_t *gl, initvar_t initvar){
  long dim,spec;
  double rho;
  Init2_t Init2;
  initvar_t values;
  reformat_initvar_species_fractions(gl, initvar, values,nd+2);
  verify_positivity_of_determinative_property(values,nd,nd+ns+4);

  Init2.P=values[nd];
  Init2.T=values[nd+1];
  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=values[nd+2+spec]*_calM(spec)/calA;
  for (spec=0; spec<ns; spec++) Init2.w[spec]=values[nd+2+spec]*_calM(spec)/calA/rho;
  for (dim=0; dim<nd; dim++) Init2.V[dim]=values[dim];
  Init2.k=values[nd+2+ns];
  Init2.psi=values[nd+3+ns];
  Init2.Tv=values[nd+4+ns];
  find_U_init_2(np, l, gl, Init2);
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
    case INIT_TYPE4: 
      init_node_4(np, l, gl, initvar); 
    break;
    case INIT_TYPE5: 
      init_node_5(np, l, gl, initvar); 
    break;
    case INIT_TYPE6: 
      init_node_6(np, l, gl, initvar); 
    break;
    default:
      fatal_error("Initial condition type %ld invalid.",inittype);
  }
}

void find_default_initvar_name(initvarname_t *initvar_name){
  char *speciesname;
  long dim,spec;
  initvarname_t string;
  speciesname=(char *)malloc(sizeof(char));

  for (dim=0; dim<nd; dim++){
    sprintf(string,"V[%ld]",dim);
    strcpy(initvar_name[dim],string);
  }
  strcpy(initvar_name[nd],"P");
  strcpy(initvar_name[nd+1],"T");
  for (spec=0; spec<ns; spec++){
    find_species_name(spec,&speciesname);
    sprintf(string,"w_%s",speciesname); 
    strcpy(initvar_name[nd+2+spec],string);
  }
  strcpy(initvar_name[nd+2+ns],"k");
  strcpy(initvar_name[nd+3+ns],"psi");
  strcpy(initvar_name[nd+4+ns],"Tv");
  free(speciesname);
} 

/* find initvar (the type of initvar must match defaultinitvartypefluid) */
void find_default_initvar(np_t *np, gl_t *gl, long l, initvar_t initvar){
  long dim,spec;
  for (dim=0; dim<nd; dim++) initvar[dim]=_V(np[l],dim);
  initvar[nd]=_P(np[l],gl);
  initvar[nd+1]=_T(np[l],gl);
  for (spec=0; spec<ns; spec++) initvar[nd+2+spec]=_w(np[l],spec);
  initvar[nd+2+ns]=_k(np[l]);
  initvar[nd+3+ns]=_psi(np[l]);
  initvar[nd+4+ns]=_Tv(np[l]);
}


