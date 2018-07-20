#define _CYCLE_METHOD "Multizone Marching [parent2002a]"
#define _CYCLE_MULTIZONE_MARCHING

#define _CYCLE_ACTIONNAME "MultizoneMarching"
#include <cycle/res/.active/res.hh>
#include <cycle/ts/.active/ts.hh>
#include <cycle/tsemf/.active/tsemf.hh>
#include <model/fluid/.active/cycle_fluid.hh>
#include <model/emfield/.active/cycle_emfield.hh>

typedef struct {
  long entrance,phi1,phi2,phi3;
  double varphiverge;
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
