#define _RESCONV_METHOD "MUSCL with flux FVS [steger1981a], FDS [roe1981a], and primitive interpolation TVD [anderson1986a], WENO [jiang1996a], CWENO [dumbser2007a], AOWENO [balsara2016a], and eigenvalue conditioning HARTEN, GNOFFO [gnoffo2004a], PECLET [parent2017a]"
#define _RESCONV_MUSCL
#define hbw_resconv_fluid 5 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid 1  //half bandwidth of stencil that depends on emfield mem props
#define _RESCONV_LAMBDA_ABSOLUTE_CONDITIONING
#define CONVJACOBIAN CONVJACOBIAN_FDS
#define _RESCONV_ACTIONNAME "MUSCL"

typedef struct {
  int FLUX,INTERPOL;
  int AOWENO_TYPE;
  double AOWENO_gammalo, AOWENO_gammahi;
  int EIGENVALCOND,AVERAGING;
} gl_cycle_resconv_t;



