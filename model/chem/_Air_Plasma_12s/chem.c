// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2021 Bernard Parent

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

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/share/chem_share.h>
#include <src/control.h>
#include "parentdunn2021.h"

#define CHEMMODEL_NONE 1
#define CHEMMODEL_PARENTDUNN2021 2



/* set all reactions to true except for testing purposes */
const static bool NEGATIVEIONREACTION[7]=
  {
   TRUE, /* reaction 0 */
   TRUE, /* reaction 1 */
   TRUE, /* reaction 2 */
   TRUE, /* reaction 3 */
   TRUE, /* reaction 4 */
   TRUE, /* reaction 5 */
   TRUE, /* reaction 6 */
  };


void write_model_chem_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    CHEMMODEL=CHEMMODEL_PARENTDUNN2021;\n"
    "    NEGATIVEIONREACTIONS=FALSE; {include reactions function of EoverN}\n"
    "    QEISOURCETERMS=TRUE; {include electron energy cooling due to electron impact}\n"
    "  );\n"
  ,_CHEM_ACTIONNAME);
}



void read_model_chem_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_CHEM_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);

    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_CHEM_ACTIONNAME);
    SOAP_add_int_to_vars(codex,"CHEMMODEL_NONE",CHEMMODEL_NONE); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_PARENTDUNN2021",CHEMMODEL_PARENTDUNN2021);
    gl->MODEL_CHEM_READ=TRUE;

    action_original=codex->action;
    codex->action=&read_model_chem_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"CHEMMODEL",&gl->model.chem.CHEMMODEL);
    if (gl->model.chem.CHEMMODEL!=CHEMMODEL_PARENTDUNN2021
        && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)
      SOAP_fatal_error(codex,"CHEMMODEL must be set to CHEMMODEL_PARENTDUNN2021.");
    find_bool_var_from_codex(codex,"NEGATIVEIONREACTIONS",&gl->model.chem.NEGATIVEIONREACTIONS);
    find_bool_var_from_codex(codex,"QEISOURCETERMS",&gl->model.chem.QEISOURCETERMS);

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}



void add_W_NegativeIon ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  double N[ns];
  double R;
  long k;
  spec_t X;
  
  R=1.9872;
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  
    
  /* find properties needed by add_to_W* functions */
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
  }


    
  if (NEGATIVEIONREACTION[1]){
  }

  
}



void add_dW_dx_NegativeIon ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, 
                  spec_t dWdTv, spec_t dWdQbeam ) {
  long k;  
  spec_t N,X;
  double R;
  
  R=1.9872;
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  

  /* find properties needed by add_to_dW* functions in proper units */
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA;
  }

  if (NEGATIVEIONREACTION[1]) {
  }
  


}





void find_W_None ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  long k;

  for ( k = 0; k < ns; k++ ) {
    W[k] = 0.0;
  }
}


void find_dW_dx_None ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  long k, s;                    /* counters */

  for ( s = 0; s < ns; s++ ) {
    dWdT[s] = 0.0;
    dWdTe[s] = 0.0;
    dWdTv[s] = 0.0;
    dWdQbeam[s] = 0.0;
    for ( k = 0; k < ns; k++ ) {
      dWdrhok[s][k] = 0.0;
    }
  }
}


void find_W ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  switch (gl->model.chem.CHEMMODEL){
    case CHEMMODEL_PARENTDUNN2021: 
      find_W_ParentDunn2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_NONE: 
      find_W_None ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );    
    break;
    default:
      fatal_error("Problem with CHEMMODEL in find_W() within chem.c");
  }
  if (gl->model.chem.NEGATIVEIONREACTIONS && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  
    add_W_NegativeIon ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
}


void find_dW_dx ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  switch (gl->model.chem.CHEMMODEL){
    case CHEMMODEL_PARENTDUNN2021: 
      find_dW_dx_ParentDunn2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_NONE: 
      find_dW_dx_None ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    default:
      fatal_error("Problem with CHEMMODEL in find_W() within chem.c");
  }
  if (gl->model.chem.NEGATIVEIONREACTIONS  && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  
    add_dW_dx_NegativeIon ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );

}





void find_Qei(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  double theta;
  
  *Qei=0.0;  
  if (gl->model.chem.QEISOURCETERMS){
    switch (gl->model.chem.CHEMMODEL){
      case CHEMMODEL_PARENTDUNN2021: 
        find_Qei_ParentDunn2021 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_NONE: 
        *Qei=0.0;
      break;
      default:
        fatal_error("Problem with CHEMMODEL in find_Qei() within chem.c");
    }

    theta=log(Estar);

    if (gl->model.chem.NEGATIVEIONREACTIONS && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  {
      if (NEGATIVEIONREACTION[1]) 
        add_to_Qei(specN2, exp(-0.0105809*sqr(theta)-2.40411e-75*pow(theta,46.0)), rhok, Qei);
      if (NEGATIVEIONREACTION[2]) 
        add_to_Qei(specO2, exp(-0.0102785*sqr(theta)-2.42260e-75*pow(theta,46.0)), rhok, Qei);
      if (NEGATIVEIONREACTION[3]) 
        add_to_Qei(specNO, exp ( -5.9890E-6 * pow ( theta , 4.0 ) + 2.5988E-84 * pow ( theta, 51.0 ) ), rhok, Qei);
      if (NEGATIVEIONREACTION[5]) 
        add_to_Qei(specN, exp(-9.3740E-3*sqr(theta)-3.3250e-23*pow(theta,14.0)), rhok, Qei);
      if (NEGATIVEIONREACTION[6]) 
        add_to_Qei(specO, exp(-1.0729E-2*sqr(theta)+1.6762E-87*pow(theta,53.0)), rhok, Qei);
    }
  } 
}


void find_dQei_dx(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  double theta;
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
  
  if (gl->model.chem.QEISOURCETERMS){
    switch (gl->model.chem.CHEMMODEL){
      case CHEMMODEL_PARENTDUNN2021: 
        find_dQei_dx_ParentDunn2021 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_NONE: 
        *dQeidTe=0.0;
      break;
      default:
        fatal_error("Problem with CHEMMODEL in find_Qei() within chem.c");
    }
  
  
  
    theta=log(Estar);

    if (gl->model.chem.NEGATIVEIONREACTIONS && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  {
      if (NEGATIVEIONREACTION[1]) 
        add_to_dQei(specN2, exp(-0.0105809*sqr(theta)-2.40411e-75*pow(theta,46.0)), 0.0, rhok, dQeidrhok, dQeidTe);
      if (NEGATIVEIONREACTION[2]) 
        add_to_dQei(specO2, exp(-0.0102785*sqr(theta)-2.42260e-75*pow(theta,46.0)), 0.0, rhok, dQeidrhok, dQeidTe);
      if (NEGATIVEIONREACTION[3]) 
        add_to_dQei(specNO, exp ( -5.9890E-6 * pow ( theta , 4.0 ) + 2.5988E-84 * pow ( theta, 51.0 ) ), 0.0, rhok, dQeidrhok, dQeidTe);
      if (NEGATIVEIONREACTION[5]) 
      add_to_dQei(specN, exp(-9.3740E-3*sqr(theta)-3.3250e-23*pow(theta,14.0)), 0.0, rhok, dQeidrhok, dQeidTe);
      if (NEGATIVEIONREACTION[6]) 
      add_to_dQei(specO, exp(-1.0729E-2*sqr(theta)+1.6762E-87*pow(theta,53.0)), 0.0, rhok, dQeidrhok, dQeidTe);
    }
  }
}


void find_We ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, 
                          double *We_create, double *We_destroy ){

  switch (gl->model.chem.CHEMMODEL){
    default:
      fatal_error("CHEMMODEL must be set to CHEMMODEL_LENARD1964 for two-temperature model of We*ee.");
  }
}
