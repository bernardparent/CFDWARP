#define _RESCONV_METHOD "Positivity-Preserving FVS [steger1981a] and flux interpolation with TVD [anderson1986a], WENO [jiang1996a], CWENO [dumbser2007a], AOWENO [balsara2016a]"
#define _RESCONV_FVSplus
#define hbw_resconv_fluid 4 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid hbw_resconv_fluid  //half bandwidth of stencil that depends on emfield mem props 
#define _RESCONV_POSITIVITY_PRESERVING
#define _RESCONV_DELTA_LAMBDA_STORAGE
#define CONVJACOBIAN CONVJACOBIAN_FVSPLUS
#define _RESCONV_ACTIONNAME "FVSplus"

typedef struct {
  int INTERPOL;
  int AOWENO_TYPE;
  double AOWENO_gammalo, AOWENO_gammahi;
  int EIGENVALCOND;
  int numiter;
} gl_cycle_resconv_t;



