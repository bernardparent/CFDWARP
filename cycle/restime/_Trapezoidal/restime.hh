#define _RESTIME_METHOD "Trapezoidal 2nd-Order"
#define _RESTIME_BW 2
#define _RESTIME_TRAPEZOIDAL
#define _RESTIME_TRAPEZOIDAL_RESIDUAL
#define _RESTIME_STORAGE_TRAPEZOIDAL
#define _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
#define UNSTEADY
#define _RESTIME_ACTIONNAME "Trapezoidal"

typedef struct {
  double weightm1_trapezoidal_default,weightm1_trapezoidal_convection;
} gl_cycle_restime_t;

