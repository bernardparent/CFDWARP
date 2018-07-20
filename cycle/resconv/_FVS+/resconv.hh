#define _RESCONV_METHOD "Positivity-Preserving FVS Van Leer Second-Order [parent2012a]"
#define _RESCONV_FVSPLUS
#define _RESCONV_DELTA_LAMBDA_STORAGE
#define _RESCONV_POSITIVITY_PRESERVING
#define _RESCONV_LAMBDA_ABSOLUTE_CONDITIONING
#define hbw_resconv_fluid 2 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid hbw_resconv_fluid  //half bandwidth of stencil that depends on emfield mem props
#define CONVJACOBIAN CONVJACOBIAN_FVSPLUS
#define _RESCONV_ACTIONNAME "FVSplus"

typedef struct {
  int LIMITER,ACCURACY;
  double xi;
  int EIGENVALCOND;
} gl_cycle_resconv_t;

