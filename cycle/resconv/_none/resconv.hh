#define _RESCONV_METHOD "None"
#define _RESCONV_NONE
#define hbw_resconv_fluid 0
#define hbw_resconvplasma_fluid hbw_resconv_fluid  //half bandwidth of stencil that depends on emfield mem props
#define _RESCONV_ACTIONNAME "None"

typedef struct {
  int CONVJACOBIAN;
} gl_cycle_resconv_t;

