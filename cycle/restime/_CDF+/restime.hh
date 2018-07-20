#define _RESTIME_METHOD "Positivity-Preserving Cross-Difference Formula [parent2018a]"
#define _RESTIME_CDFPLUS
#define _RESTIME_BW 2
#define _RESTIME_CDF
#define UNSTEADY
//#define _RESTIME_STORAGE_TRAPEZOIDAL
//#define _RESTIME_CDF_TRAPEZOIDAL

#define _RESTIME_ACTIONNAME "CDFplus"

typedef struct {
  double xi1,xi2,xi3;
} gl_cycle_restime_t;

