#include <src/common.h>
#include <model/emfield/_emfield.h>
#include <model/_model.h>
#include <src/control.h>


void write_model_beam_template(FILE **controlfile){
}



void read_and_init_beam_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_model_beam_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  gl->MODEL_BEAM_READ=TRUE;  
}




double _Qbeam(np_t np, gl_t *gl){
  double Qbeam;
    Qbeam=0.0e0;
  return(Qbeam);
}


double _dQbeam_dN(np_t np, gl_t *gl){
  return(0.0);
}
