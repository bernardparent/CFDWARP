#define _RESCONV_METHOD "FVS [steger1981a] and flux interpolation with TVD [anderson1986a], WENO [jiang1996a], CWENO [dumbser2007a], AOWENO [balsara2016a]"
#define _RESCONV_FVS
#define hbw_resconv_fluid 4 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid hbw_resconv_fluid  //half bandwidth of stencil that depends on emfield mem props 
#define CONVJACOBIAN CONVJACOBIAN_FVS
#define _RESCONV_ACTIONNAME "FVS"

typedef struct {
  int INTERPOL;
  int AOWENO_TYPE;
  double AOWENO_gammalo, AOWENO_gammahi;
  int EIGENVALCOND;
} gl_cycle_resconv_t;



