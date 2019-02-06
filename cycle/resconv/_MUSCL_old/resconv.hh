#define _RESCONV_METHOD "Multidimensional MUSCL with flux FVS [steger1981a], FDS [roe1981a], and primitive interpolation TVD [anderson1986a], WENO [jiang1996a], CWENO [dumbser2007a], AOWENO [balsara2016a], and eigenvalue conditioning HARTEN, GNOFFO [gnoffo2004a], PECLET [parent2017a]"
#define _RESCONV_MUSCL
#define hbw_resconv_fluid 5 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid 1  //half bandwidth of stencil that depends on emfield mem props
#define CONVJACOBIAN CONVJACOBIAN_FDS
#define _RESCONV_ACTIONNAME "MUSCL"
#define _RESCONV_STORAGE_FSTAR

typedef struct {
  int FLUX,INTERPOL;
  int AOWENO_TYPE;
  double AOWENO_gammalo, AOWENO_gammahi;
  int EIGENVALCOND,AVERAGING,FACEINTEG;
} gl_cycle_resconv_t;



