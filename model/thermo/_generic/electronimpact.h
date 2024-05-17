#ifndef ELECTRONIMPACT_H

#define ELECTRONIMPACT_H
#include <src/common.h>
  
double _EoverN_from_Nk_Te(spec_t N, double Te);

double _Te_from_Nk_EoverN(spec_t N, double EoverN);

double _EoverNk_from_Te(long spec, double Te);

double _Tek_from_EoverN(long spec, double EoverN);


#endif
