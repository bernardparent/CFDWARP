#ifndef THERMO_H

#define THERMO_H
#include <src/common.h>
  

#define THERMO_RN2 296.8E0
#define THERMO_THETAV 3353.0E0
#define THERMO_EVLIM 0.1
#define THERMO_TVLIM (THERMO_THETAV/log(THERMO_RN2*THERMO_THETAV/THERMO_EVLIM+1.0))

double _numatoms(long spec);


#endif
