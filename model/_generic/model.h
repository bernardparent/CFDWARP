#ifndef MODEL_H

  #define MODEL_H
  #include <src/common.h>

  void add_to_clipped_variables(gl_t *gl, char *str);

  void find_clipped_variables_list(gl_t *gl, char **cliplist);

  void find_clipped_muscl_variables_list(gl_t *gl, char **cliplist);

  void find_clipped_bdry_variables_list(gl_t *gl, char **cliplist);
 
  void write_model_template(FILE **controlfile);
 
  void read_model(char *argum, SOAP_codex_t *codex);
  
#endif
