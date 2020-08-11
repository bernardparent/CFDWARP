#ifndef LENARD1964_H

#define LENARD1964_H
#include <src/common.h>
  
void find_W_Lenard1964 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                  double Qbeam, spec_t W );

void find_dW_dx_Lenard1964 ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, 
                  double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, 
                  spec_t dWdTv, spec_t dWdQbeam );

void find_dQei_dx_Lenard1964(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe);

void find_Qei_Lenard1964(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei);

void find_We_Lenard1964 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, 
                          double *We_create, double *We_destroy );

#endif
