#define _RESCONV_METHOD "Flux Difference Splitting [roe1981a] with Symmetric Total Variation Diminishing [yee1990a], and eigenvalue conditioning HARTEN, GNOFFO [gnoffo2004a], PECLET [parent2017a]"
#define _RESCONV_FDS
#define _RESCONV_LAMBDA_ABSOLUTE_CONDITIONING
#define hbw_resconv_fluid 2 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid hbw_resconv_fluid  //half bandwidth of stencil that depends on emfield mem props
#define _RESCONV_DELTA_LAMBDA_STORAGE
#define CONVJACOBIAN CONVJACOBIAN_FDS
#define _RESCONV_ACTIONNAME "FDS"

typedef struct {
  int AVERAGING,ACCURACY;
  int EIGENVALCOND;
} gl_cycle_resconv_t;

