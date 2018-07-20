#define _RESCONV_METHOD "FDS Cauchy Second-Order [parent2016b], and eigenvalue conditioning HARTEN, GNOFFO [gnoffo2004a], PECLET [parent2017a]"
#define _RESCONV_FDS_CAUCHY2
#define _RESCONV_ZETAADEN
#define _RESCONV_LAMBDA_ABSOLUTE_CONDITIONING
#define hbw_resconv_fluid 2 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid hbw_resconv_fluid  //half bandwidth of stencil that depends on emfield mem props
#define _RESCONV_DELTA_LAMBDA_STORAGE
#define CONVJACOBIAN CONVJACOBIAN_FDS
#define _RESCONV_ACTIONNAME "FDSCauchy"

typedef struct {
  int AVERAGING,ACCURACY;
  int EIGENVALCOND;
} gl_cycle_resconv_t;

