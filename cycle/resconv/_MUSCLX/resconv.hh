#define _RESCONV_METHOD "MUSCLX with flux FVS [steger1981a], FDS [roe1981a], and primitive interpolation TVD [anderson1986a], WENO [jiang1996a], CWENO [dumbser2007a], AOWENO [balsara2016a], and eigenvalue conditioning HARTEN, GNOFFO [gnoffo2004a], PECLET [parent2017a]"
#define _RESCONV_MUSCLX
#define hbw_resconv_fluid 5 //maximum half bandwidth of stencil (excluding center node)
#define hbw_resconvplasma_fluid 1  //half bandwidth of stencil that depends on emfield mem props
#define _RESCONV_ACTIONNAME "MUSCLX"
#define _RESCONV_STORAGE_FSTAR

typedef struct {
  int FLUX,INTERPOL;
  int AOWENO_TYPE;
  double AOWENO_gammalo, AOWENO_gammahi;
  int CONVJACOBIAN,EIGENVALCOND,AVERAGING,FACEINTEG;
} gl_cycle_resconv_t;



