#define _RESCONV_METHOD "Positivity-Preserving FDS Symmetric Total Variation Diminishing 2nd Order [parent2013a], and eigenvalue conditioning HARTEN, GNOFFO [gnoffo2004a], PECLET [parent2017a]"
#define _RESCONV_FDSPLUS
#define _RESCONV_DELTA_LAMBDA_STORAGE
#define _RESCONV_LAMBDA_ABSOLUTE_CONDITIONING
#define _RESCONV_POSITIVITY_PRESERVING
#define hbw_resconv_fluid 2   //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid hbw_resconv_fluid  //half bandwidth of stencil that depends on emfield mem props
#define CONVJACOBIAN CONVJACOBIAN_FDS
#define _RESCONV_ACTIONNAME "FDSplus"

typedef struct {
  int AVERAGING,numiter,ACCURACY;
  double xi;
  int EIGENVALCOND;
} gl_cycle_resconv_t;

