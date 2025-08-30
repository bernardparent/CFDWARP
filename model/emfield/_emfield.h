#ifndef _EMFIELD_H
#define _EMFIELD_H

#include <src/common.h>


  void write_disc_emfield_template(FILE **controlfile);

  void write_cycle_emfield_template(FILE **controlfile);

  void write_model_emfield_template(FILE **controlfile);

  void write_init_emfield_template(FILE **controlfile);

  void read_model_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex);
  
  void read_cycle_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex);

  void read_disc_emfield_actions(char *actionname, char **argum, SOAP_codex_t *codex);

  void find_B(np_t np, gl_t *gl, EXM_vec3D_t B);

  void find_E(np_t np, gl_t *gl, EXM_vec3D_t E);
    
  void find_Ek(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Ek);

  void find_Ek_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, long spec, EXM_vec3D_t Ek);

  void find_J(np_t *np, gl_t *gl, long l, EXM_vec3D_t J);

  void find_post_variable_name_emfield(long varnum, char *varname);

  void find_post_variable_value_emfield(np_t *np, long l, gl_t *gl, long varnum, double *varvalue);

  void find_emfield_force(np_t *np, gl_t *gl, long l, dim_t Femfield);

  double _sigma(np_t *np, gl_t *gl, long l);

  void find_DeltaVk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t DeltaVk);

  void find_DeltaVk_thermal(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t DeltaVk);

  void find_DeltaVk_drift_from_mu(np_t *np, gl_t *gl, long l, long spec, double mu, EXM_vec3D_t DeltaVk);

  double _phi(np_t np, gl_t *gl);

  double _Te_EMF(np_t np, gl_t *gl);

  void find_Ve_from_J(np_t *np, gl_t *gl, long l, EXM_vec3D_t Ve);

  void find_Vk_from_Vk_at_interfaces(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Vk);

#ifdef UNSTEADY
  double _phim1(np_t np, gl_t *gl);
#endif

  double _Ekmag(np_t *np, gl_t *gl, long l, long spec);

  double _Eemag_smoothed(np_t *np, gl_t *gl, long l);

  void find_E_smoothed(np_t *np, gl_t *gl, long l, dim_t E);

  bool is_node_bdry_dielectric(np_t np);

  void find_DeltaVk_from_mu(np_t *np, gl_t *gl, long l, long spec, double mu, EXM_vec3D_t DeltaVk);

#endif 
