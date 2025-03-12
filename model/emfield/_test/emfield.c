// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2023 Bernard Parent

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
#include <model/emfield/_emfield.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/transport/_transport.h>
#include <model/chem/_chem.h>
#include <model/metrics/_metrics.h>
#include <model/fluid/_fluid.h>
#include <cycle/_cycle.h>
#include <src/control.h>
#include <src/bdry.h>

#define BDRYEMF_ELECTRODE 0

#define INITEMF_TYPE1 1


void write_bdry_emfield_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n",_EMFIELD_ACTIONNAME
  );
  wfprintf(*controlfile,
  "    {\n"
  "    _________________________________________________________________________________________\n"
  "\n"
  "    EMField Bdry Condition Type   ID   Description\n"
  "    _________________________________________________________________________________________\n"
  "\n"
  "    BDRYEMF_ELECTRODE             %c    Electrode, phi fixed\n"
  "    _________________________________________________________________________________________\n"
  "\n"
  "    }\n"
  "    All(BDRYEMF_ELECTRODE);\n"
  "    {\n"
  "    Region(is" if2DL(",js") if3DL(",ks") ",  ie" if2DL(",je") if3DL(",ke") ",  BDRYEMF_ELECTRODE);\n"
  "    }\n"
  "  );\n",_bdry_ID(BDRYEMF_ELECTRODE));
}



void write_init_emfield_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n"
  ,_EMFIELD_ACTIONNAME);
  wfprintf(*controlfile,
  "    {\n"
  "    _______________________________________________________________________________________\n"
  "\n"
  "    EMField Initial Condition Type         Parameters\n"
  "    _______________________________________________________________________________________\n"
  "\n"
  "    INITEMF_TYPE1                          phi\n"
  "    _______________________________________________________________________________________\n"
  "\n"
  "    }\n"
  "    phi=0.0; {V}\n"
  "    All(INITEMF_TYPE1,phi);\n"
  "    {\n"
  "    Bdry(BDRYEMF_ELECTRODE, INITEMF_TYPE1,phi);\n"
  "    Region(is,"if2DL("js,")if3DL("ks,")"  ie,"if2DL("je,")if3DL("ke,")"INITEMF_TYPE1,phi);\n"
  "    }\n"
  "  );\n"
  );
}



void find_prim_emfield_mem_1(np_t *np, gl_t *gl, long l){
  long m,n,linner;
  EXM_vec3D_t E;
  bool INNER_NODE_FOUND;
  
  INNER_NODE_FOUND=find_l_of_nearest_inner_node(np, gl, l, TYPELEVEL_EMFIELD, &linner);
  
  if (INNER_NODE_FOUND){
      assert_np(np[linner],is_node_inner(np[linner],TYPELEVEL_EMFIELD));
      
      for (n=0; n<3; n++) E[n]=0.0e0;
      for (n=0; n<nd; n++) {
        for (m=0; m<nd; m++){
              E[n]-=_X(np[linner],m,n)*0.5e0*
	       (_phi(np[_al(gl,linner,m,+1)],gl)-_phi(np[_al(gl,linner,m,-1)],gl));
        }
      }
      for (n=0; n<3; n++) {
        np[l].bs->E[n]=E[n];
      }      

  } else {
    for (n=0; n<3; n++) {
      np[l].bs->E[n]=0.0;
    }
  }     


}


void find_prim_emfield_mem_2(np_t *np, gl_t *gl, long l){

}


void find_prim_emfield_mem_3(np_t *np, gl_t *gl, long l){

}



void add_init_types_emfield_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"INITEMF_TYPE1",   INITEMF_TYPE1);
}



void add_bdry_types_emfield_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"BDRYEMF_ELECTRODE",   BDRYEMF_ELECTRODE);
}



void write_model_emfield_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    sigmasolid=1000; {S/m}\n"
    "    sigmafluid=1000; {S/m}\n"
    "    epsilonr=1.0; \n"
    "  );\n",_EMFIELD_ACTIONNAME
  );
}

void write_cycle_emfield_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    xiverge_emfield=1e-2; {residual convergence threshold for the potential equation}\n"
    "    sigmaref=3e-6; {reference value for the conductivity in S/m -> used to determine the local time step of the potential equation}\n" 
    "    Uref_emfield[1]=100.0e0; {reference value for the electric field potential in Volts -> used to determine xi}\n" 
    "  );\n"
  ,_EMFIELD_ACTIONNAME);
}


void write_disc_emfield_template(FILE **controlfile){

}

void init_node_emfield(np_t np, gl_t *gl, long inittype, initvar_emfield_t initvar){
  switch (inittype) {
    case INITEMF_TYPE1: 
      np.bs->Uemfield[0]=initvar[0]; 
    break;
    default:
      fatal_error("EMF Initial condition type %ld invalid.",inittype);
  }
}

/* find initvar_emfield (the type of initvar must match defaultinitvartypeemfield) */
void find_default_initvar_emfield(np_t *np, gl_t *gl, long l, initvar_emfield_t initvar){
  initvar[0]=np[l].bs->Uemfield[0];
}


void find_default_initvar_name_emfield(initvarname_t *initvar_name_emfield){
  strcpy(initvar_name_emfield[0],"phi");
}

void read_disc_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->DISC_EMFIELD_READ=TRUE;  
}



void read_model_emfield_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_EMFIELD_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_EMFIELD_ACTIONNAME);
    if (!gl->CONTROL_READ){
      gl->MODEL_EMFIELD_READ=TRUE;  
    }

    action_original=codex->action;
    codex->action=&read_model_emfield_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_double_var_from_codex(codex,"epsilonr",&gl->model.emfield.epsilonr);
    find_double_var_from_codex(codex,"sigmasolid",&gl->model.emfield.sigmasolid);
    find_double_var_from_codex(codex,"sigmafluid",&gl->model.emfield.sigmafluid);
    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}



void read_cycle_emfield_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){

}


void read_cycle_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  char tmpstr[200];
  long flux;
  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_EMFIELD_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE)  wfprintf(stdout,"%s..",_EMFIELD_ACTIONNAME);
    gl->CYCLE_EMFIELD_READ=TRUE;  
    action_original=codex->action;
    codex->action=&read_cycle_emfield_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;
    find_double_var_from_codex(codex,"xiverge_emfield",&gl->cycle.emfield.xiverge);
    find_double_var_from_codex(codex,"sigmaref",&gl->cycle.emfield.sigmaref);
    for (flux=0; flux<nfe; flux++){
      sprintf(tmpstr,"Uref_emfield[%ld]",flux+1);
      find_double_var_from_codex(codex,tmpstr,&gl->cycle.emfield.Uref[flux]);
    }
    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}



double _phi(np_t np, gl_t *gl){
  return(np.bs->Uemfield[0]);
}


#ifdef UNSTEADY
double _phim1(np_t np, gl_t *gl){
  return(np.bs->Uemfieldm1[0]);
}
#endif



void find_emfield_force(np_t *np, gl_t *gl, long l, dim_t Femfield){
  long dim;
  for (dim=0; dim<nd; dim++) Femfield[dim]=0.0;
}


double _E_dot_J(np_t *np, gl_t *gl, long l){
  double EdotJ;
  EdotJ=0.0;
  return(EdotJ);
}

double _E_dot_Je(np_t *np, gl_t *gl, long l){
  double EdotJe;
  EdotJe=0.0;
  return(EdotJe);
}



void find_sigma(np_t *np, gl_t *gl, long l, double *sigma){
  if (is_node_valid(np[l],TYPELEVEL_FLUID)) {
    *sigma=gl->model.emfield.sigmafluid;
  } else {
    *sigma=gl->model.emfield.sigmasolid;
  }
}


void find_B(np_t np, gl_t *gl, EXM_vec3D_t B){
  long dim;
  for (dim=0; dim<3; dim++) B[dim]=0.0;
}

void find_E(np_t np, gl_t *gl, EXM_vec3D_t E){
  long dim;
  assert_np(np,is_node_valid(np,TYPELEVEL_EMFIELD));
  for (dim=0; dim<3; dim++) E[dim]=np.bs->E[dim];
}


void find_Ek(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Ek){
  find_E(np[l], gl, Ek);
}


/* find dphi/dXj at interface between node l and node _al(gl,l,i,+1) */
double _dphidXj(np_t *np, gl_t *gl, long l, long i, long j){
  double dphidXj;
  assert(is_node_valid(np[l],TYPELEVEL_EMFIELD));
  assert(is_node_valid(np[_al(gl,l,i,+1)],TYPELEVEL_EMFIELD));

  if (i!=j) {
    assert_np(np[_al(gl,l,j,+1)],is_node_valid(np[_al(gl,l,j,+1)],TYPELEVEL_EMFIELD));
    assert_np(np[_al(gl,l,j,-1)],is_node_valid(np[_al(gl,l,j,-1)],TYPELEVEL_EMFIELD));
    assert_np(np[_all(gl,l,j,+1,i,+1)],is_node_valid(np[_all(gl,l,j,+1,i,+1)],TYPELEVEL_EMFIELD));
    assert_np(np[_all(gl,l,j,-1,i,+1)],is_node_valid(np[_all(gl,l,j,-1,i,+1)],TYPELEVEL_EMFIELD));
    dphidXj=0.25*(_phi(np[_al(gl,l,j,+1)],gl)      -_phi(np[_al(gl,l,j,-1)],gl)
                  +_phi(np[_all(gl,l,j,+1,i,+1)],gl)-_phi(np[_all(gl,l,j,-1,i,+1)],gl));
  } else {
    assert_np(np[_al(gl,l,j,+1)],is_node_valid(np[_al(gl,l,j,+1)],TYPELEVEL_EMFIELD));
    assert_np(np[_al(gl,l,j,+0)],is_node_valid(np[_al(gl,l,j,+0)],TYPELEVEL_EMFIELD));
    dphidXj=(_phi(np[_al(gl,l,j,+1)],gl)-_phi(np[_al(gl,l,j,+0)],gl));
  }
  return(dphidXj);
}


void find_E_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, EXM_vec3D_t E){
  long n,m;
  metrics_t metrics;
  find_metrics_at_interface(np, gl, lL, lR, theta, &metrics);
  for (n=0; n<3; n++) E[n]=0.0e0;
  for (n=0; n<nd; n++) {
    for (m=0; m<nd; m++){
      if (m==theta){
        assert_np(np[lR],is_node_valid(np[lR],TYPELEVEL_EMFIELD));
        assert_np(np[lL],is_node_valid(np[lL],TYPELEVEL_EMFIELD));
        assert_np(np[lL],is_node_in_domain_lim(lL,gl));
        assert_np(np[lR],is_node_in_domain_lim(lR,gl));
        E[n]-=metrics.X2[m][n]*
	           (_phi(np[lR],gl)-_phi(np[lL],gl));
      } else {
        if (is_node_valid(np[_al(gl,lR,m,+1)],TYPELEVEL_EMFIELD) 
            && is_node_valid(np[_al(gl,lL,m,+1)],TYPELEVEL_EMFIELD)
            && is_node_valid(np[_al(gl,lR,m,-1)],TYPELEVEL_EMFIELD)
            && is_node_valid(np[_al(gl,lL,m,-1)],TYPELEVEL_EMFIELD)
            ){
          assert_np(np[_al(gl,lR,m,+1)],is_node_in_domain_lim(_al(gl,lR,m,+1),gl));
          assert_np(np[_al(gl,lR,m,-1)],is_node_in_domain_lim(_al(gl,lR,m,-1),gl));
          assert_np(np[_al(gl,lL,m,+1)],is_node_in_domain_lim(_al(gl,lL,m,+1),gl));
          assert_np(np[_al(gl,lL,m,-1)],is_node_in_domain_lim(_al(gl,lL,m,-1),gl));
          E[n]-=metrics.X2[m][n]*0.25*
	           (_phi(np[_al(gl,lR,m,+1)],gl)-_phi(np[_al(gl,lR,m,-1)],gl)
             +_phi(np[_al(gl,lL,m,+1)],gl)-_phi(np[_al(gl,lL,m,-1)],gl) );
        } else {
          if (is_node_valid(np[_al(gl,lR,m,+1)],TYPELEVEL_EMFIELD) 
              && is_node_valid(np[_al(gl,lL,m,+1)],TYPELEVEL_EMFIELD)
              && is_node_valid(np[_al(gl,lR,m,+0)],TYPELEVEL_EMFIELD)
              && is_node_valid(np[_al(gl,lL,m,+0)],TYPELEVEL_EMFIELD)
              ){
            assert_np(np[_al(gl,lR,m,+1)],is_node_in_domain_lim(_al(gl,lR,m,+1),gl));
            assert_np(np[_al(gl,lR,m,+0)],is_node_in_domain_lim(_al(gl,lR,m,+0),gl));
            assert_np(np[_al(gl,lL,m,+1)],is_node_in_domain_lim(_al(gl,lL,m,+1),gl));
            assert_np(np[_al(gl,lL,m,+0)],is_node_in_domain_lim(_al(gl,lL,m,+0),gl));
            E[n]-=metrics.X2[m][n]*0.5*
	             (_phi(np[_al(gl,lR,m,+1)],gl)-_phi(np[_al(gl,lR,m,+0)],gl)
               +_phi(np[_al(gl,lL,m,+1)],gl)-_phi(np[_al(gl,lL,m,+0)],gl) );
          } else {
            assert_np(np[_al(gl,lR,m,+0)],is_node_valid(np[_al(gl,lR,m,+0)],TYPELEVEL_EMFIELD));
            assert_np(np[_al(gl,lR,m,-1)],is_node_valid(np[_al(gl,lR,m,-1)],TYPELEVEL_EMFIELD));
            assert_np(np[_al(gl,lL,m,+0)],is_node_valid(np[_al(gl,lL,m,+0)],TYPELEVEL_EMFIELD));
            assert_np(np[_al(gl,lL,m,-1)],is_node_valid(np[_al(gl,lL,m,-1)],TYPELEVEL_EMFIELD));
            assert_np(np[_al(gl,lR,m,+0)],is_node_in_domain_lim(_al(gl,lR,m,+0),gl));
            assert_np(np[_al(gl,lR,m,-1)],is_node_in_domain_lim(_al(gl,lR,m,-1),gl));
            assert_np(np[_al(gl,lL,m,+0)],is_node_in_domain_lim(_al(gl,lL,m,+0),gl));
            assert_np(np[_al(gl,lL,m,-1)],is_node_in_domain_lim(_al(gl,lL,m,-1),gl));
            E[n]-=metrics.X2[m][n]*0.5*
	             (_phi(np[_al(gl,lR,m,+0)],gl)-_phi(np[_al(gl,lR,m,-1)],gl)
               +_phi(np[_al(gl,lL,m,+0)],gl)-_phi(np[_al(gl,lL,m,-1)],gl) );
          }
        }
      }
    }
  }
}




void find_Ee_for_Townsend_ionization(np_t *np, gl_t *gl, long l, EXM_vec3D_t E){
  long i;
  for (i=0; i<3; i++) E[i]=0.0;
}



void find_Estar_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, double *Estar){
  long j,n;
  double dphidXj;
  metrics_t metrics;

  find_metrics_at_interface(np,gl,lL,lR,theta,&metrics);
  *Estar=0.0;
  for (j=0; j<nd; j++){
    dphidXj=_dphidXj(np, gl, lL, theta, j);
    for (n=0; n<nd; n++){
      *Estar+=-metrics.X2[j][n]*metrics.X2[theta][n]*dphidXj;
    }
  }
}



void find_DeltaVk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t DeltaVk){
  long i;
  for (i=0; i<3; i++) DeltaVk[i]=0.0;  
}

void find_Vk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Vk){
  dim_t Vn;
  long dim;
  find_V(np[l],Vn);    
  for (dim=0; dim<3; dim++) Vk[dim]=0.0;
  for (dim=0; dim<nd; dim++) Vk[dim]=Vn[dim];
}


void find_linearization_coefficients_bdry_node_emfield(np_t *np, gl_t *gl, long lA, long theta, long thetasgn,
                        long flux, long bdrytype, double *valA, double *valB, double *valRHS){
  
  *valA=1.0e0;
  *valB=0.0e0;
  *valRHS=0.0e0;
}




double _sigma(np_t *np, gl_t *gl, long l){
  double sigma;
  if (is_node_valid(np[l],TYPELEVEL_FLUID)) {
    sigma=gl->model.emfield.sigmafluid;
  } else {
    sigma=gl->model.emfield.sigmasolid;
  }
  return(sigma);
}


/* find sigma between nodes np[l] and np[al(gl,l,dim,+1)] */
double _sigma_interface(np_t *np, gl_t *gl, long l, long dim){
  double sigmaint;
  sigmaint=0.5*(_sigma(np,gl,l)+_sigma(np,gl,_al(gl,l,dim,+1)));
  return(sigmaint);
}


void find_linearization_coefficients_inner_node_emfield(np_t *np, gl_t *gl, long l, long flux, double *coeff, EXM_gl3D_t *coeffgl){
  long n,coeffl,inter,coefflL,coefflR,theta,j,lL,lR,i,k;
  double sigmaint;
  metrics_t metrics;
#ifdef UNSTEADY
  double epsilonr;
#endif
  
  coeffgl->is=-2;
  coeffgl->ie=+2;
  coeffgl->js=-2;
  coeffgl->je=+2;
#ifdef _3D
  coeffgl->ks=-2;
  coeffgl->ke=+2;
#else
  coeffgl->ks=0;
  coeffgl->ke=0;
#endif
  for (i=coeffgl->is; i<=coeffgl->ie; i++){
    for (j=coeffgl->js; j<=coeffgl->je; j++){
      for (k=coeffgl->ks; k<=coeffgl->ke; k++){
        coeff[EXM_ai3(*coeffgl,i,j,k)]=0.0;
      }
    }
  }
  coeffl=EXM_ai3(*coeffgl,0,0,0);
  
  for (theta=0; theta<nd; theta++){
    for (inter=-1; inter<=+1; inter+=2){
      if (inter==+1){
        lL=l;
        lR=_al(gl,l,theta,+1);
        coefflL=coeffl;
        coefflR=EXM_al3(*coeffgl,coeffl,theta,+1);
      } else {
        lL=_al(gl,l,theta,-1);
        lR=l;        
        coefflL=EXM_al3(*coeffgl,coeffl,theta,-1);
        coefflR=coeffl;
      }
     
      find_metrics_at_interface(np, gl, lL, lR, theta, &metrics);
      
      
      sigmaint=_sigma_interface(np, gl, lL, theta);
      for (j=0; j<nd; j++){
        for (n=0; n<nd; n++){
          if (theta!=j) {
//            *Jstar+=-sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*0.25*(_phi(np[_al(gl,lL,j,+1)],gl)      -_phi(np[_al(gl,lL,j,-1)],gl)
//                  +_phi(np[_all(gl,lL,j,+1,theta,+1)],gl)-_phi(np[_all(gl,lL,j,-1,theta,+1)],gl));
            coeff[EXM_al3(*coeffgl,coefflL,j,+1)]-=(double)inter*sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*0.25*metrics.Omega/_Omega(np[_al(gl,lL,j,+1)],gl);
            coeff[EXM_al3(*coeffgl,coefflL,j,-1)]+=(double)inter*sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*0.25*metrics.Omega/_Omega(np[_al(gl,lL,j,-1)],gl);
            coeff[EXM_all3(*coeffgl,coefflL,j,+1,theta,+1)]-=(double)inter*sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*0.25*metrics.Omega/_Omega(np[_all(gl,lL,j,+1,theta,+1)],gl);
            coeff[EXM_all3(*coeffgl,coefflL,j,-1,theta,+1)]+=(double)inter*sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*0.25*metrics.Omega/_Omega(np[_all(gl,lL,j,-1,theta,+1)],gl);
          } else {
//            *Jstar+=-sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*(_phi(np[_al(gl,lL,j,+1)],gl)-_phi(np[_al(gl,lL,j,+0)],gl));
            coeff[EXM_al3(*coeffgl,coefflL,j,+1)]-=(double)inter*sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*metrics.Omega/_Omega(np[_al(gl,lL,j,+1)],gl);
            coeff[EXM_al3(*coeffgl,coefflL,j,+0)]+=(double)inter*sigmaint*metrics.X2[j][n]*metrics.X2[theta][n]*metrics.Omega/_Omega(np[_al(gl,lL,j,+0)],gl);            
          }
        }
      }
      
      
#ifdef UNSTEADY
      epsilonr=gl->model.emfield.epsilonr;
      for (j=0; j<nd; j++){
        if (theta!=j) {
          for (n=0; n<nd; n++){
//          dphidXj=_phi(np[_al(gl,lL,j,+1)],gl)      -_phi(np[_al(gl,lL,j,-1)],gl)
//                  +_phi(np[_al(gl,lR,j,+1)],gl)-_phi(np[_al(gl,lR,j,-1)],gl);
         //   Jstar+=-0.25*epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]*dphidXj/notzero(gl->dt,1e-99);
            coeff[EXM_al3(*coeffgl,coefflL,j,+1)]-= +(double)inter*0.25*epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]/notzero(gl->dt,1e-99)*metrics.Omega/_Omega(np[_al(gl,lL,j,+1)],gl);
            coeff[EXM_al3(*coeffgl,coefflL,j,-1)]-= -(double)inter*0.25*epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]/notzero(gl->dt,1e-99)*metrics.Omega/_Omega(np[_al(gl,lL,j,-1)],gl);
            coeff[EXM_al3(*coeffgl,coefflR,j,+1)]-= +(double)inter*0.25*epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]/notzero(gl->dt,1e-99)*metrics.Omega/_Omega(np[_al(gl,lR,j,+1)],gl);
            coeff[EXM_al3(*coeffgl,coefflR,j,-1)]-= -(double)inter*0.25*epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]/notzero(gl->dt,1e-99)*metrics.Omega/_Omega(np[_al(gl,lR,j,-1)],gl);            
          }
        } else {
          for (n=0; n<nd; n++){
//          dphidXj=(_phi(np[lR],gl)-_phi(np[lL],gl));
        //    Jstar+=-epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]*dphidXj/notzero(gl->dt,1e-99);
            coeff[EXM_al3(*coeffgl,coefflR,j,+0)]+= -(double)inter*epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]/notzero(gl->dt,1e-99)*metrics.Omega/_Omega(np[_al(gl,lR,j,+0)],gl);
            coeff[EXM_al3(*coeffgl,coefflL,j,+0)]+= +(double)inter*epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[theta][n]/notzero(gl->dt,1e-99)*metrics.Omega/_Omega(np[_al(gl,lL,j,+0)],gl);
          }
        }
        
      }
#endif
    }
  }
}


void find_dtau_emfield(np_t *np, gl_t *gl, long l, long flux, double *dtau){
  long i,j;
  double sum,dx,dxmin,dxmax;
  long dim;
  double sigmap1h,sigmam1h,sigmamax;

  assert_np(np[l],is_node_inner(np[l],TYPELEVEL_EMFIELD));
  dxmin=1e99;
  dxmax=0.0;
  for (i=0; i<nd; i++) {
    sum=0.0;
    for (j=0; j<nd; j++) sum+=sqr(_X(np[l],i,j));
    sum=sqrt(sum);
    dxmin=min(dxmin,fabs(1.0/sum));
    dxmax=max(dxmax,fabs(1.0/sum));
  }
  dx=sqrt(dxmin*dxmax);

  sigmamax=0.0;
  for (dim=0; dim<nd; dim++){
    sigmap1h=_sigma_interface(np, gl, l, dim); 
    sigmam1h=_sigma_interface(np, gl, _al(gl,l,dim,-1), dim);
#ifdef UNSTEADY
    sigmap1h+=epsilon0*gl->model.emfield.epsilonr/gl->dt; 
    sigmam1h+=epsilon0*gl->model.emfield.epsilonr/gl->dt;
#endif
    sigmamax=max(max(sigmamax,sigmap1h),sigmam1h);
  }
  *dtau=gl->Lc*dx/(sigmamax+gl->cycle.emfield.sigmaref); 

  assert_np(np[l],*dtau>0.0);
}



/* find the TDMA coefficients for the phi equation along theta dimension */
void find_linearization_coefficients_inner_node_emfield_interface(np_t *np, gl_t *gl, long l, long theta, long flux, double *C_1, double *C_2, double *C_3){
  metrics_t metricsp1h,metricsm1h;
  double sigmam1h,sigmap1h,Xsquarem1h,Xsquarep1h;
  long n;
#ifdef UNSTEADY
  double epsilonr_m1h,epsilonr_p1h;
#endif

  assert(flux==0);
  find_metrics_at_interface(np, gl, l, _al(gl,l,theta,+1), theta, &metricsp1h);
  find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), l, theta, &metricsm1h);

  Xsquarem1h=0.0;
  Xsquarep1h=0.0;
  for (n=0; n<nd; n++){
    Xsquarem1h+=sqr(metricsm1h.X2[theta][n]);
    Xsquarep1h+=sqr(metricsp1h.X2[theta][n]);
  }

  sigmam1h=_sigma_interface(np, gl, _al(gl,l,theta,-1), theta);
  sigmap1h=_sigma_interface(np, gl, l, theta);

  *C_1=-sigmam1h*metricsm1h.Omega*Xsquarem1h/_Omega(np[_al(gl,l,theta,-1)],gl);
  *C_2=sigmam1h*metricsm1h.Omega*Xsquarem1h/_Omega(np[l],gl)
      +sigmap1h*metricsp1h.Omega*Xsquarep1h/_Omega(np[l],gl);
  *C_3=-sigmap1h*metricsp1h.Omega*Xsquarep1h/_Omega(np[_al(gl,l,theta,+1)],gl);

#ifdef UNSTEADY
  assert(_Omega(np[_al(gl,l,theta,-1)],gl)!=0.0);
  assert(_Omega(np[_al(gl,l,theta,+1)],gl)!=0.0);
  assert(_Omega(np[l],gl)!=0.0);
  assert(gl->dt!=0.0);
  epsilonr_m1h=gl->model.emfield.epsilonr;
  epsilonr_p1h=gl->model.emfield.epsilonr;
  *C_1+=-epsilon0*epsilonr_m1h*metricsm1h.Omega*Xsquarem1h/_Omega(np[_al(gl,l,theta,-1)],gl)/gl->dt;
  *C_2+=epsilon0*epsilonr_m1h*metricsm1h.Omega*Xsquarem1h/_Omega(np[l],gl)/gl->dt
      +epsilon0*epsilonr_p1h*metricsp1h.Omega*Xsquarep1h/_Omega(np[l],gl)/gl->dt;
  *C_3+=-epsilon0*epsilonr_p1h*metricsp1h.Omega*Xsquarep1h/_Omega(np[_al(gl,l,theta,+1)],gl)/gl->dt;
  
#endif



}

#ifdef UNSTEADY
/* find dphim1/dXj at interface between node l and node _al(gl,l,i,+1) */
double _dphim1dXj(np_t *np, gl_t *gl, long l, long i, long j){
  double dphidXj;
  if (i!=j) {
    assert_np(np[_al(gl,l,j,+1)],is_node_valid(np[_al(gl,l,j,+1)],TYPELEVEL_EMFIELD));
    assert_np(np[_al(gl,l,j,-1)],is_node_valid(np[_al(gl,l,j,-1)],TYPELEVEL_EMFIELD));
    assert_np(np[_all(gl,l,j,+1,i,+1)],is_node_valid(np[_all(gl,l,j,+1,i,+1)],TYPELEVEL_EMFIELD));
    assert_np(np[_all(gl,l,j,-1,i,+1)],is_node_valid(np[_all(gl,l,j,-1,i,+1)],TYPELEVEL_EMFIELD));
    dphidXj=0.25*(_phim1(np[_al(gl,l,j,+1)],gl)      -_phim1(np[_al(gl,l,j,-1)],gl)
                  +_phim1(np[_all(gl,l,j,+1,i,+1)],gl)-_phim1(np[_all(gl,l,j,-1,i,+1)],gl));
  } else {
    assert_np(np[_al(gl,l,j,+1)],is_node_valid(np[_al(gl,l,j,+1)],TYPELEVEL_EMFIELD));
    assert_np(np[_al(gl,l,j,+0)],is_node_valid(np[_al(gl,l,j,+0)],TYPELEVEL_EMFIELD));
    dphidXj=(_phim1(np[_al(gl,l,j,+1)],gl)-_phim1(np[_al(gl,l,j,+0)],gl));
  }
  return(dphidXj);
}
#endif


#ifdef UNSTEADY
/* find ith component of the displacement current between node l and node _al(gl,l,i,+1) */
double _Jdispstar_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long i){
  long j,n;
  double Jstar,epsilonr;
  double dphidXj,dphim1dXj;

  Jstar=0.0;
  epsilonr=gl->model.emfield.epsilonr;
  for (j=0; j<nd; j++){
    dphidXj=_dphidXj(np, gl, l, i, j);
    dphim1dXj=_dphim1dXj(np, gl, l, i, j);
    for (n=0; n<nd; n++){
      Jstar+=-epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[i][n]*dphidXj/notzero(gl->dt,1e-99);
      Jstar-=-epsilon0*epsilonr*metrics.X2[j][n]*metrics.X2[i][n]*dphim1dXj/notzero(gl->dt,1e-99);
    }
  }
  return(Jstar);
}
#endif



/* find ith component of the current between node l and node _al(gl,l,i,+1) */
static double _Jstar_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long i){
  long lL,lR;
  double Estar,Jstar;
  lL=l;
  lR=_al(gl,l,i,+1);
  find_Estar_at_interface(np, gl, l, lR, i, &Estar);
  Jstar=Estar*_sigma_interface(np, gl, l, i);
  return(Jstar);
}


static double _Jtstar_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long i){
  double Jtstar;
  Jtstar=_Jstar_interface(np, gl, metrics, l, i);
#ifdef UNSTEADY
  Jtstar+=_Jdispstar_interface(np, gl, metrics, l, i);
#endif
  return(Jtstar);
}




void find_Fstar_interface_emfield(np_t *np, gl_t *gl, long lL, long lR, long i, fluxemfield_t Fstar){
  metrics_t metrics;
  find_metrics_at_interface(np, gl, lL, lR, i, &metrics);
  Fstar[0]=_Jtstar_interface(np,gl,metrics,lL,i)*metrics.Omega;
}



void find_Sstar_emfield(np_t *np, gl_t *gl, long l, fluxemfield_t S){
  S[0]=0.0;
}



void update_bdry_emfield(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, bool BDRYDIRECFOUND, int TYPELEVEL){

 /* NOTE: No relaxation of the properties is allowed at the emfield boundary nodes */
  if (_node_type(np[lA],TYPELEVEL)==BDRYEMF_ELECTRODE) {
    /* phi is specified, so don't do anything */
  }
}


double _xi_emfield(np_t np, gl_t *gl, fluxemfield_t Res){
  double xiphi;
  assert_np(np,is_node_inner(np,TYPELEVEL_EMFIELD));
  xiphi=fabs(Res[0]/_Omega(np,gl)/gl->cycle.emfield.Uref[0]);
  return(xiphi);
}

void find_post_variable_name_emfield(long varnum, char *varname){
  switch (varnum) {
    case 0:  
      sprintf(varname,"phi");
    break;
    case 1:
      sprintf(varname,"sigma");
    break;
  }
  if (varnum>=2 && varnum<=1+nd) sprintf(varname,"E[%ld]",varnum-2);
}


void find_post_variable_value_emfield(np_t *np, long l, gl_t *gl, long varnum, double *varvalue){
  EXM_vec3D_t E;
  double sigma;
  long n;

  if (is_node_valid(np[l],TYPELEVEL_EMFIELD)) {
    find_E(np[l], gl, E);
  } else {
    for (n=0; n<3; n++) {
      E[n]=0.0e0;
    }
  }
  sigma=_sigma(np,gl,l);

  *varvalue=0.0;
  if (is_node_valid(np[l],TYPELEVEL_EMFIELD)){
    switch (varnum) {
      case 0:
        *varvalue=_phi(np[l],gl);
      break;
      case 1:
        *varvalue=sigma;
      break;
    }
    if (varnum>=2 && varnum<=1+nd) *varvalue=E[varnum-2];
  }
}


void add_dUstar_to_U_emfield(np_t *np, gl_t *gl, long l, fluxemfield_t dUstar){
  assert_np(np[l],is_node_inner(np[l],TYPELEVEL_EMFIELD));
  np[l].bs->Uemfield[0]=np[l].bs->Uemfield[0]+dUstar[0]/_Omega(np[l],gl);
}
