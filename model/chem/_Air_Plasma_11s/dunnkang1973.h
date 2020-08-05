#ifndef DUNNKANG1973_H

#define DUNNKANG1973_H
#include <src/common.h>
  
void find_W_DunnKang1973 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                  double Qbeam, spec_t W );

void find_dW_dx_DunnKang1973 ( gl_t *gl, spec_t rhok, spec_t mu, double T, double Te, double Tv, 
                  double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, 
                  spec_t dWdTv, spec_t dWdQbeam );

#endif
