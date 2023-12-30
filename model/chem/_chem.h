#ifndef _CHEM_H
#define _CHEM_H

#include <src/common.h>


void write_model_chem_template(FILE **controlfile);

void read_model_chem_actions(char *actionname, char **argum, SOAP_codex_t *codex);

void find_W(np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W);

void find_dW_dx(np_t np, gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam,
              spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam);

void find_Qei(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei);

void find_dQei_dx(np_t np, gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe);

void find_We ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, 
                          double *We_create, double *We_destroy );


#ifdef _AVERAGEDRATES_CHEM
void find_Qchemave (np_t np, gl_t *gl, spec_t rhok, double Estar, double T, double Te, double *Q );
#endif

#endif /* _CHEM_H */ 
