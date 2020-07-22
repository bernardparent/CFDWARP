#ifndef _EMFIELD_SHARE_H
#define _EMFIELD_SHARE_H

#include <model/_model.h>
#include <src/common.h>
#define EMFIELD_PRESSUREGRADIENT TRUE

  double _sigma_interface(np_t *np, gl_t *gl, long l, long dim);

  double _Jstar_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long i);

  double _Jtstar_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long i);

  void find_J_at_node_from_Jstar_at_interfaces(np_t *np, gl_t *gl, long l, double(*_Jstar_funct)(np_t *, gl_t *, metrics_t, long, long), EXM_vec3D_t J);

  void find_Ee_for_Townsend_ionization(np_t *np, gl_t *gl, long l, EXM_vec3D_t E);

  void find_Estar_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, double *Estar);

  void find_E_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, EXM_vec3D_t E);
 
  void find_Jstar_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, double *Jstar);

  double _Qjoulek(np_t *np, gl_t *gl, long l, long k);

  double _E_dot_J(np_t *np, gl_t *gl, long l);

  double _E_dot_J_recast(np_t *np, gl_t *gl, long l);

  double _E_dot_Je_over_mueNe(np_t *np, gl_t *gl, long l);

  double _E_dot_Je(np_t *np, gl_t *gl, long l);

  void find_Vk(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Vk);

  double _dTedXj_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long i, long j);

  double _epsilonr(np_t np, gl_t *gl);
    
  double _epsilonr_interface(np_t *np, gl_t *gl, long l, long theta);

  double _Nk_max(np_t *np, gl_t *gl, long l, long spec);

  void find_Vk_from_mu(np_t *np, gl_t *gl, long l, double mu, long spec, EXM_vec3D_t Vk);

  void find_Vk_at_interface(np_t *np, gl_t *gl, long l, long theta, long spec, EXM_vec3D_t Vk);

  void find_Ve_from_J(np_t *np, gl_t *gl, long l, EXM_vec3D_t Ve);

  void find_Ve_from_Je(np_t *np, gl_t *gl, long l, EXM_vec3D_t Ve);

  void find_Vk_from_Vk_at_interfaces(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Vk);

  void find_Jk_from_Jk_at_interfaces(np_t *np, gl_t *gl, long l, long spec, EXM_vec3D_t Jk);

#endif /* _EMFIELD_SHARE_H */
