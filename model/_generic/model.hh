#define _MODEL_METHOD "GENERIC"
#define _MODEL_GENERIC

#include <model/thermo/.active/thermo.hh>
#include <model/chem/.active/chem.hh>
#include <model/fluid/.active/fluid.hh>
#include <model/metrics/.active/metrics.hh>
#include <model/emfield/.active/emfield.hh>
#include <model/beam/.active/beam.hh>


#define ControlActionName_Model "Model"
#define numpostvar (totalpostvarfluid+totalpostvaremfield)
#define numinitvar (totalinitvarfluid)
#ifdef EMFIELD
  #define numinitvar_emfield (totalinitvaremfield)
#endif



typedef struct {
  dim_t X;
  dim2_t X2;
  double Omega;
} metrics_t;


typedef char* clipname_t;
typedef long clipnum_t;

typedef struct {
  long clipnumtot;
  long clipnamenum;
  clipname_t *clipname;
  clipnum_t *clipnum;
  thread_lock_t clip_lock;
#ifdef DISTMPI
  long clipnumtot_all;
  long clipnamenum_all;
  clipname_t *clipname_all;
  clipnum_t *clipnum_all;
#endif
  
#ifdef EMFIELD
  gl_model_emfield_t emfield;
#endif
  gl_model_fluid_t fluid;
  gl_model_chem_t chem;
} gl_model_t;
