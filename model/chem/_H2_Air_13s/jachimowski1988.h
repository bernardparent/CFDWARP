#ifndef JACHIMOWSKI1988_H

#define JACHIMOWSKI1988_H
#include <src/common.h>
  
void find_W_Jachimowski1988 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, 
                  double Qbeam, spec_t W );

void find_dW_dx_Jachimowski1988 ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, 
                  double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, 
                  spec_t dWdTv, spec_t dWdQbeam );

#endif
