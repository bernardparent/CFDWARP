#ifndef _BEAM_H
#define _BEAM_H

#include <src/common.h>

  void write_model_beam_template(FILE **controlfile);
  
  void read_model_beam_actions(char *actionname, char **argum, SOAP_codex_t *codex);
  
  double _Qbeam(np_t np, gl_t *gl);

  double _dQbeam_dN(np_t np, gl_t *gl);
  
  
#endif 
