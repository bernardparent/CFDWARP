#ifndef _CHEM_H
#define _CHEM_H

#include <src/common.h>


void find_W(spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W);

void find_dW_dx(spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam,
              spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, 
	      spec_t dWdQbeam);

void test_dW_dx(spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam);


#endif /* _CHEM_H */ 
