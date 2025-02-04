#ifndef ELECTRONIMPACT_H

#define ELECTRONIMPACT_H
#include <src/common.h>
  
double _EoverN_from_Nk_Te(spec_t N, double Te);

double _Te_from_Nk_EoverN(spec_t N, double EoverN);

double _EoverNk_from_Te(long spec, double Te);

double _Tek_from_EoverN(long spec, double EoverN);

double _fk_from_Tv(long k, double Tv);

double _QeN2vib_from_rhok_T_Tv_Te(gl_t *gl, spec_t rhok, double T, double Tv, double Te);

#endif
