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

  double _phi(np_t np, gl_t *gl);
    
  void find_B(np_t np, gl_t *gl, EXM_vec3D_t B);

  void find_E(np_t np, gl_t *gl, EXM_vec3D_t E);
    
  void find_J(np_t *np, gl_t *gl, long l, EXM_vec3D_t J);

  void find_post_variable_name_emfield(long varnum, char *varname);

  void find_post_variable_value_emfield(np_t *np, long l, gl_t *gl, long varnum, double *varvalue);

  void find_Ek(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Ek);

  void find_emfield_force(np_t *np, gl_t *gl, long l, dim_t Femfield);

  double _Qjoulek(np_t *np, gl_t *gl, long l, long k);

  double _E_dot_J(np_t *np, gl_t *gl, long l);

  double _E_dot_J_recast(np_t *np, gl_t *gl, long l);

  double _E_dot_Je_over_mueNe(np_t *np, gl_t *gl, long l);

  void find_Vk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Vk);

  double _mu(np_t *np, gl_t *gl, long l, long spec);

  double _sigma(np_t *np, gl_t *gl, long l);

  void find_Ee_for_Townsend_ionization(np_t *np, gl_t *gl, long l, EXM_vec3D_t E);

  void find_J_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, EXM_vec3D_t J);

  void find_Jstar_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, double *Jstar);

  void find_Estar_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, double *Estar);

  void find_E_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, EXM_vec3D_t E);
 
  void find_Ek_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, long spec, EXM_vec3D_t Ek);


/*  
  void find_dTe_dNk(np_t *np, gl_t *gl, long l, spec_t dTedNk);
*/

  long _time_marching_method_emfield(long fluxe, gl_t *gl);
  
  double _dTedXj_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long i, long j);

  void find_DeltaVk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t DeltaVk);

  double _Ekmag(np_t *np, gl_t *gl, long l, long spec);

  double _Eemag_smoothed(np_t *np, gl_t *gl, long l);

  double _epsilonr(np_t np, gl_t *gl);
    
  double _epsilonr_interface(np_t *np, gl_t *gl, long l, long theta);

#endif 
