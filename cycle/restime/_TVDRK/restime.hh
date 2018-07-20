#define _RESTIME_METHOD "Total Variation Diminishing Runge Kutta [gottlieb1998a]"
#define _RESTIME_TVDRK
#define _RESTIME_BW 2
#define _RESTIME_PREDICTOR_CORRECTOR
#define UNSTEADY


#define _RESTIME_ACTIONNAME "TVDRK"

typedef struct {
  int METHOD;
} gl_cycle_restime_t;

