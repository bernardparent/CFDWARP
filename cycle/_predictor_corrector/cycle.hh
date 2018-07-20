#define _CYCLE_METHOD "Predictor-Corrector"
#define _CYCLE_PREDICTOR_CORRECTOR

#define  _CYCLE_ACTIONNAME "PredictorCorrector"
#include <cycle/res/.active/res.hh>
#include <cycle/ts/.active/ts.hh>
#include <cycle/tsemf/.active/tsemf.hh>
#include <model/fluid/.active/cycle_fluid.hh>
#include <model/emfield/.active/cycle_emfield.hh>


typedef struct {
  bool RUNTIMEMODULEFOUND;
  char *code_runtime;
  SOAP_codex_t codex;
  gl_cycle_fluid_t fluid;
  gl_cycle_resconv_t resconv;
  gl_cycle_restime_t restime;

#ifdef EMFIELD
  gl_cycle_emfield_t emfield;
#endif

} gl_cycle_t;
