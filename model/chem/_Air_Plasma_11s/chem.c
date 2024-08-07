// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2018,2021 Bernard Parent
Copyright 2021,2022 Prasanna Thoguluva Rajendran

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
#include "dunnkang1973.h"
#include "park1993.h"
#include "boyd2007.h"
#include "parentdunn2021.h"
#include "parentpark2021.h"
#include "lenard1964.h"
#include "farbar2013.h"
#include "kim2021.h"
#include "parent2023.h"
#include "parent2023b.h"
#include "thoguluva2023.h"
#include "rodriguez2024.h"
#include "rodriguezpark2024.h"

#define CHEMMODEL_NONE 1
#define CHEMMODEL_DUNNKANG1973 2
#define CHEMMODEL_PARK1993 3
#define CHEMMODEL_BOYD2007 4
#define CHEMMODEL_LENARD1964 5
#define CHEMMODEL_FARBAR2013 6
#define CHEMMODEL_PARENTDUNN2021 7
#define CHEMMODEL_PARENTPARK2021 8
#define CHEMMODEL_KIM2021 9
#define CHEMMODEL_PARENT2023 10
#define CHEMMODEL_PARENT2023B 11
#define CHEMMODEL_THOGULUVA2023 12
#define CHEMMODEL_RODRIGUEZ2024 13
#define CHEMMODEL_RODRIGUEZPARK2024 14

#define Estarmin 1e-40



void write_model_chem_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    CHEMMODEL=CHEMMODEL_DUNNKANG1973;\n"
    "    for (react,1,6, ADDITIONALREACTION[react]=FALSE;); {include reactions function of EoverN}\n"
    "    TOWNSENDIONIZATIONIMPLICIT=FALSE; {keep this to FALSE generally}\n"
    "    QEISOURCETERMS=TRUE; {include electron energy cooling due to electron impact}\n"
    "  );\n"
  ,_CHEM_ACTIONNAME);
}



void read_model_chem_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  long numvarsinit,cnt;
  char addreactstr[100];
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_CHEM_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);

    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_CHEM_ACTIONNAME);
    SOAP_add_int_to_vars(codex,"CHEMMODEL_NONE",CHEMMODEL_NONE); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_DUNNKANG1973",CHEMMODEL_DUNNKANG1973); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_PARK1993",CHEMMODEL_PARK1993); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_BOYD2007",CHEMMODEL_BOYD2007); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_LENARD1964",CHEMMODEL_LENARD1964);
    SOAP_add_int_to_vars(codex,"CHEMMODEL_FARBAR2013",CHEMMODEL_FARBAR2013);
    SOAP_add_int_to_vars(codex,"CHEMMODEL_PARENTDUNN2021",CHEMMODEL_PARENTDUNN2021);
    SOAP_add_int_to_vars(codex,"CHEMMODEL_PARENTPARK2021",CHEMMODEL_PARENTPARK2021);
    SOAP_add_int_to_vars(codex,"CHEMMODEL_KIM2021",CHEMMODEL_KIM2021); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_PARENT2023",CHEMMODEL_PARENT2023); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_PARENT2023B",CHEMMODEL_PARENT2023B); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_THOGULUVA2023",CHEMMODEL_THOGULUVA2023); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_RODRIGUEZ2024",CHEMMODEL_RODRIGUEZ2024); 
    SOAP_add_int_to_vars(codex,"CHEMMODEL_RODRIGUEZPARK2024",CHEMMODEL_RODRIGUEZPARK2024); 
    gl->MODEL_CHEM_READ=TRUE;

    action_original=codex->action;
    codex->action=&read_model_chem_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"CHEMMODEL",&gl->model.chem.CHEMMODEL);
    if (gl->model.chem.CHEMMODEL!=CHEMMODEL_DUNNKANG1973 && gl->model.chem.CHEMMODEL!=CHEMMODEL_PARK1993 
        && gl->model.chem.CHEMMODEL!=CHEMMODEL_BOYD2007 && gl->model.chem.CHEMMODEL!=CHEMMODEL_LENARD1964
        && gl->model.chem.CHEMMODEL!=CHEMMODEL_FARBAR2013 && gl->model.chem.CHEMMODEL!=CHEMMODEL_PARENTDUNN2021
        && gl->model.chem.CHEMMODEL!=CHEMMODEL_PARENTPARK2021 && gl->model.chem.CHEMMODEL!=CHEMMODEL_KIM2021
        && gl->model.chem.CHEMMODEL!=CHEMMODEL_PARENT2023 && gl->model.chem.CHEMMODEL!=CHEMMODEL_PARENT2023B 
        &&  gl->model.chem.CHEMMODEL!=CHEMMODEL_THOGULUVA2023 && gl->model.chem.CHEMMODEL!=CHEMMODEL_RODRIGUEZ2024 && gl->model.chem.CHEMMODEL!=CHEMMODEL_RODRIGUEZPARK2024 && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)
      SOAP_fatal_error(codex,"CHEMMODEL must be set to either CHEMMODEL_DUNNKANG1973 or CHEMMODEL_NONE or CHEMMODEL_BOYD2007 or CHEMMODEL_PARK1993 or CHEMMODEL_LENARD1964 or CHEMMODEL_FARBAR2013 or CHEMMODEL_PARENTDUNN2021 or CHEMMODEL_PARENTPARK2021 or CHEMMODEL_KIM2021 or CHEMMODEL_PARENT2023 or CHEMMODEL_PARENT2023B or CHEMMODEL_THOGULUVA2023 or CHEMMODEL_RODRIGUEZ2024 or CHEMMODEL_RODRIGUEZPARK2024.");
    for (cnt=1; cnt<=6; cnt++) {
      sprintf(addreactstr,"ADDITIONALREACTION[%ld]",cnt);
      find_bool_var_from_codex(codex,addreactstr,&gl->model.chem.ADDITIONALREACTION[cnt]);
    }
    find_bool_var_from_codex(codex,"TOWNSENDIONIZATIONIMPLICIT",&gl->model.chem.TOWNSENDIONIZATIONIMPLICIT);
    find_bool_var_from_codex(codex,"QEISOURCETERMS",&gl->model.chem.QEISOURCETERMS);

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }
}



void add_W_Additional ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  double N[ns];
  double theta,R;
  long k;
  spec_t X;
  double Estar_from_Te,Te_from_Estar;
  
  R=Rchem;
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  
  Estar_from_Te=_EoverN_from_rhok_Te(rhok, Te);
  //if (Estar_from_Te>Estar) Estar=Estar_from_Te;   //??? needs to be verified
  Te_from_Estar=_Te_from_rhok_EoverN(rhok, Estar);
  
  
  Te_from_Estar=max(300.0,Te_from_Estar); 
    
  /* find properties needed by add_to_W* functions */
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM (k ) * 1e-6 * calA;  /* particules/cm^3 */
  }
  Estar = max ( Estarmin, Estar );
  theta = log ( Estar );


    
  if (gl->model.chem.ADDITIONALREACTION[1])
      add_to_W_2r3p ( specN2, speceminus,   specN2plus, speceminus, speceminus,   exp ( -0.0105809 * sqr ( theta ) - 2.40411e-75 * pow ( theta, 46.0 ) ), N, W);

  if (gl->model.chem.ADDITIONALREACTION[2])
      add_to_W_2r3p ( specO2, speceminus,   specO2plus, speceminus, speceminus,   exp ( -0.0102785 * sqr ( theta ) - 2.42260e-75 * pow ( theta, 46.0 ) ), N, W);

  if (gl->model.chem.ADDITIONALREACTION[3])
      add_to_W_2r3p ( specNO, speceminus,   specNOplus, speceminus, speceminus,   exp ( -5.9890E-6 * pow ( theta , 4.0 ) + 2.5988E-84 * pow ( theta, 51.0 ) ), N, W);

  if (gl->model.chem.ADDITIONALREACTION[4]) 
      add_to_W_fw_3r2p ( specNOplus, speceminus, speceminus,   specNO, speceminus, 2.2e40, -4.5, 0.0*R, Te, X, W );

  if (gl->model.chem.ADDITIONALREACTION[5])
      add_to_W_2r3p ( specN, speceminus,   specNplus, speceminus, speceminus,   exp(-9.3740E-3*sqr(theta)-3.3250e-23*pow(theta,14.0)), N, W);

  if (gl->model.chem.ADDITIONALREACTION[6])
      add_to_W_2r3p ( specO, speceminus,   specOplus, speceminus, speceminus,   exp(-1.0729E-2*sqr(theta)+1.6762E-87*pow(theta,53.0)), N, W);
  
}



void add_dW_dx_Additional ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, 
                  spec_t dWdTv, spec_t dWdQbeam ) {
  long k;  
  spec_t N,X;
  double R,theta,kf,dkfdTe,dkfdT,dkfdTv,Te_from_Estar;
  
  R=Rchem;
  for ( k = 0; k < ns; k++ ) {
    X[k] = rhok[k] / _calM ( k ) * 1.0e-06;     /* mole/cm3 */
  }
  

  /* find properties needed by add_to_dW* functions in proper units */
  for ( k = 0; k < ns; k++ ) {
    N[k] = rhok[k] / _calM ( k ) * 1e-6 * calA;
  }
  Estar = max ( Estarmin, Estar );
  theta = log ( Estar );

  Te_from_Estar=_Te_from_rhok_EoverN(rhok, Estar);
  Te_from_Estar=max(300.0,Te_from_Estar); 
  

  if (gl->model.chem.ADDITIONALREACTION[1] && gl->model.chem.TOWNSENDIONIZATIONIMPLICIT) {
    kf=exp ( -0.0105809 * sqr ( theta ) - 2.40411e-75 * pow ( theta, 46.0 ) );
    dkfdTe = 0.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specN2, speceminus,   specN2plus, speceminus, speceminus,   kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (gl->model.chem.ADDITIONALREACTION[2] && gl->model.chem.TOWNSENDIONIZATIONIMPLICIT) {
    kf=exp ( -0.0102785 * sqr ( theta ) - 2.42260e-75 * pow ( theta, 46.0 ) );
    dkfdTe = 0.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specO2, speceminus,   specO2plus, speceminus, speceminus,   kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (gl->model.chem.ADDITIONALREACTION[3] && gl->model.chem.TOWNSENDIONIZATIONIMPLICIT){
    kf=exp ( -5.9890E-6 * pow ( theta , 4.0 ) + 2.5988E-84 * pow ( theta, 51.0 ) );
    dkfdTe = 0.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specNO, speceminus,   specNOplus, speceminus, speceminus,   kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }

  if (gl->model.chem.ADDITIONALREACTION[4]){
    add_to_dW_fw_3r2p ( specNOplus, speceminus, speceminus,   specNO, speceminus, 2.2e40, -4.5, 0.0*R, Te, X, dWdTe, dWdrhok );
  }

  if (gl->model.chem.ADDITIONALREACTION[5] && gl->model.chem.TOWNSENDIONIZATIONIMPLICIT) {
    kf=exp(-9.3740E-3*sqr(theta)-3.3250e-23*pow(theta,14.0));
    dkfdTe = 0.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specN, speceminus,   specNplus, speceminus, speceminus,   kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
  }
  
  if (gl->model.chem.ADDITIONALREACTION[6] && gl->model.chem.TOWNSENDIONIZATIONIMPLICIT) {
    kf=exp(-1.0729E-2*sqr(theta)+1.6762E-87*pow(theta,53.0));
    dkfdTe = 0.0;
    dkfdT=0.0;
    dkfdTv=0.0;
    add_to_dW_2r3p ( specO, speceminus,   specOplus, speceminus, speceminus,   kf, N, dkfdT, dkfdTv, dkfdTe, dWdrhok, dWdT, dWdTv, dWdTe);
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


void find_W ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {
  switch (gl->model.chem.CHEMMODEL){
    case CHEMMODEL_DUNNKANG1973: 
      find_W_DunnKang1973 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_PARK1993: 
      find_W_Park1993 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_BOYD2007: 
      find_W_Boyd2007 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_LENARD1964: 
      find_W_Lenard1964 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_FARBAR2013: 
      find_W_Farbar2013 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_PARENTDUNN2021: 
      find_W_ParentDunn2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_PARENTPARK2021: 
      find_W_ParentPark2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_KIM2021: 
      find_W_Kim2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_PARENT2023: 
      find_W_Parent2023 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_PARENT2023B: 
      find_W_Parent2023b ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_THOGULUVA2023: 
      find_W_Thoguluva2023 ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_RODRIGUEZ2024: 
      find_W_Rodriguez2024 (np, gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_RODRIGUEZPARK2024: 
      find_W_RodriguezPark2024 (np, gl, rhok, T, Te, Tv, Estar, Qbeam, W );
    break;
    case CHEMMODEL_NONE: 
      find_W_None ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );    
    break;
    default:
      fatal_error("Problem with CHEMMODEL in find_W() within chem.c");
  }
  if (gl->model.chem.ADDITIONALREACTION && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  
    add_W_Additional ( gl, rhok, T, Te, Tv, Estar, Qbeam, W );
}


void find_dW_dx ( np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam,
                  spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {
  switch (gl->model.chem.CHEMMODEL){
    case CHEMMODEL_DUNNKANG1973: 
      find_dW_dx_DunnKang1973 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_PARK1993: 
      find_dW_dx_Park1993 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_BOYD2007: 
      find_dW_dx_Boyd2007 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_LENARD1964: 
      find_dW_dx_Lenard1964 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_FARBAR2013: 
      find_dW_dx_Farbar2013 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_PARENTDUNN2021: 
      find_dW_dx_ParentDunn2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_PARENTPARK2021: 
      find_dW_dx_ParentPark2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_KIM2021: 
      find_dW_dx_Kim2021 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_PARENT2023: 
      find_dW_dx_Parent2023 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_PARENT2023B: 
      find_dW_dx_Parent2023b ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_THOGULUVA2023: 
      find_dW_dx_Thoguluva2023 ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_RODRIGUEZ2024: 
      find_dW_dx_Rodriguez2024 (np, gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_RODRIGUEZPARK2024: 
      find_dW_dx_RodriguezPark2024 (np, gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    case CHEMMODEL_NONE: 
      find_dW_dx_None ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );
    break;
    default:
      fatal_error("Problem with CHEMMODEL in find_W() within chem.c");
  }
  if (gl->model.chem.ADDITIONALREACTION  && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  
    add_dW_dx_Additional ( gl, rhok, T, Te, Tv, Estar, Qbeam, dWdrhok, dWdT, dWdTe, dWdTv, dWdQbeam );

}





void find_Qei(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei){
  double theta;
  
  *Qei=0.0;  
  if (gl->model.chem.QEISOURCETERMS){
    switch (gl->model.chem.CHEMMODEL){
      case CHEMMODEL_DUNNKANG1973: 
        find_Qei_DunnKang1973 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_PARK1993: 
        find_Qei_Park1993 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_BOYD2007: 
        find_Qei_Boyd2007 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_LENARD1964: 
        find_Qei_Lenard1964 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_FARBAR2013: 
        find_Qei_Farbar2013 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_PARENTDUNN2021: 
        find_Qei_ParentDunn2021 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_PARENTPARK2021: 
        find_Qei_ParentPark2021 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_KIM2021: 
        find_Qei_Kim2021 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_PARENT2023: 
        find_Qei_Parent2023 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_PARENT2023B: 
        find_Qei_Parent2023b ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_THOGULUVA2023: 
        find_Qei_Thoguluva2023 ( gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_RODRIGUEZ2024: 
        find_Qei_Rodriguez2024 (np, gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_RODRIGUEZPARK2024: 
        find_Qei_RodriguezPark2024 (np, gl, rhok, Estar, Te, Qei );
      break;
      case CHEMMODEL_NONE: 
        *Qei=0.0;
      break;
      default:
        fatal_error("Problem with CHEMMODEL in find_Qei() within chem.c");
    }

    theta=log(Estar);

    if (gl->model.chem.ADDITIONALREACTION && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  {
      if (gl->model.chem.ADDITIONALREACTION[1]) 
        add_to_Qei(gl,Te,specN2,_ionizationpot(specN2), exp(-0.0105809*sqr(theta)-2.40411e-75*pow(theta,46.0)), rhok, Qei);
      if (gl->model.chem.ADDITIONALREACTION[2]) 
        add_to_Qei(gl,Te,specO2,_ionizationpot(specO2), exp(-0.0102785*sqr(theta)-2.42260e-75*pow(theta,46.0)), rhok, Qei);
      if (gl->model.chem.ADDITIONALREACTION[3]) 
        add_to_Qei(gl,Te,specNO,_ionizationpot(specNO), exp ( -5.9890E-6 * pow ( theta , 4.0 ) + 2.5988E-84 * pow ( theta, 51.0 ) ), rhok, Qei);
      if (gl->model.chem.ADDITIONALREACTION[5]) 
        add_to_Qei(gl,Te,specN,_ionizationpot(specN), exp(-9.3740E-3*sqr(theta)-3.3250e-23*pow(theta,14.0)), rhok, Qei);
      if (gl->model.chem.ADDITIONALREACTION[6]) 
        add_to_Qei(gl,Te,specO,_ionizationpot(specO), exp(-1.0729E-2*sqr(theta)+1.6762E-87*pow(theta,53.0)), rhok, Qei);
    }
  } 
}


void find_dQei_dx(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe){
  double theta;
  long spec;
  
  for (spec=0; spec<ns; spec++) dQeidrhok[spec]=0.0;
  *dQeidTe=0.0;  
  
  if (gl->model.chem.QEISOURCETERMS){
    switch (gl->model.chem.CHEMMODEL){
      case CHEMMODEL_DUNNKANG1973: 
        find_dQei_dx_DunnKang1973 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_PARK1993: 
        find_dQei_dx_Park1993 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_BOYD2007: 
        find_dQei_dx_Boyd2007 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_LENARD1964: 
        find_dQei_dx_Lenard1964 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_FARBAR2013: 
        find_dQei_dx_Farbar2013 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_PARENTDUNN2021: 
        find_dQei_dx_ParentDunn2021 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_PARENTPARK2021: 
        find_dQei_dx_ParentPark2021 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_KIM2021: 
        find_dQei_dx_Kim2021 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_PARENT2023: 
        find_dQei_dx_Parent2023 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_PARENT2023B: 
        find_dQei_dx_Parent2023b ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_THOGULUVA2023: 
        find_dQei_dx_Thoguluva2023 ( gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_RODRIGUEZ2024: 
        find_dQei_dx_Rodriguez2024 (np, gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_RODRIGUEZPARK2024: 
        find_dQei_dx_RodriguezPark2024 (np, gl, rhok, Estar, Te, dQeidrhok, dQeidTe );
      break;
      case CHEMMODEL_NONE: 
        *dQeidTe=0.0;
      break;
      default:
        fatal_error("Problem with CHEMMODEL in find_Qei() within chem.c");
    }
  
  
  
    theta=log(Estar);

    if (gl->model.chem.ADDITIONALREACTION && gl->model.chem.CHEMMODEL!=CHEMMODEL_NONE)  {
      if (gl->model.chem.ADDITIONALREACTION[1]) 
        add_to_dQei(gl,Te,specN2,_ionizationpot(specN2), exp(-0.0105809*sqr(theta)-2.40411e-75*pow(theta,46.0)), 0.0, rhok, dQeidrhok, dQeidTe);
      if (gl->model.chem.ADDITIONALREACTION[2]) 
        add_to_dQei(gl,Te,specO2,_ionizationpot(specO2), exp(-0.0102785*sqr(theta)-2.42260e-75*pow(theta,46.0)), 0.0, rhok, dQeidrhok, dQeidTe);
      if (gl->model.chem.ADDITIONALREACTION[3]) 
        add_to_dQei(gl,Te,specNO,_ionizationpot(specNO), exp ( -5.9890E-6 * pow ( theta , 4.0 ) + 2.5988E-84 * pow ( theta, 51.0 ) ), 0.0, rhok, dQeidrhok, dQeidTe);
      if (gl->model.chem.ADDITIONALREACTION[5]) 
      add_to_dQei(gl,Te,specN,_ionizationpot(specN), exp(-9.3740E-3*sqr(theta)-3.3250e-23*pow(theta,14.0)), 0.0, rhok, dQeidrhok, dQeidTe);
      if (gl->model.chem.ADDITIONALREACTION[6]) 
      add_to_dQei(gl,Te,specO,_ionizationpot(specO), exp(-1.0729E-2*sqr(theta)+1.6762E-87*pow(theta,53.0)), 0.0, rhok, dQeidrhok, dQeidTe);
    }
  }
}


void find_We ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, 
                          double *We_create, double *We_destroy ){

  switch (gl->model.chem.CHEMMODEL){
    case CHEMMODEL_LENARD1964: 
      find_We_Lenard1964 ( gl, rhok, T, Te, Tv, Estar, Qbeam, We_create, We_destroy );
    break;
    default:
      fatal_error("CHEMMODEL must be set to CHEMMODEL_LENARD1964 for two-temperature model of We*ee.");
  }
}
