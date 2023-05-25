#ifndef PARENT2023B_H

#define PARENT2023B_H
#include <src/common.h>
  
void find_W_Parent2023b ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                  double Qbeam, spec_t W );

void find_dW_dx_Parent2023b ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, 
                  spec_t dWdTv, spec_t dWdQbeam );

void find_dQei_dx_Parent2023b(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe);

void find_Qei_Parent2023b(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei);

#endif
