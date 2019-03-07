#define _RESCONV_METHOD "Reconstruction-Evolution MUSCL with flux FVS [steger1981a], FDS [roe1981a], and primitive interpolation TVD [anderson1986a], WENO [jiang1996a], CWENO [dumbser2007a], AOWENO [balsara2016a], and eigenvalue conditioning HARTEN, GNOFFO [gnoffo2004a], PECLET [parent2017a], PASCAL [parent2018a], and PARENT positivity-preserving filter [parent2019a]"
#define _RESCONV_MUSCL
#define _RESCONV_MUSCLPLUS
#define _RESCONV_DELTA_LAMBDA_STORAGE
#define _RESCONV_POSITIVITY_PRESERVING
//#define _RESCONV_INCLUDES_DIFFUSION
#define hbw_resconv_fluid 5 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid 1  //half bandwidth of stencil that depends on emfield mem props
#define _RESCONV_ACTIONNAME "MUSCL"

typedef struct {
  int FLUX,INTERPOL,CONVJACOBIAN;
  int AOWENO_TYPE;
  double AOWENO_gammalo, AOWENO_gammahi;
  int EIGENVALCOND,AVERAGING,FACEINTEG,POSFILTER_numiter,POSFILTER;
} gl_cycle_resconv_t;



