#ifndef FLUID_H

  #define FLUID_H
  #include <src/common.h>
  #include <model/fluid/_fluid.h>


  double _P (np_t np, gl_t *gl);

  double _q (np_t np);

  double _a (np_t np, gl_t *gl);

  double _rho (np_t np);

  double _athermo(np_t np, gl_t *gl);

  void find_dP_dx(np_t np, gl_t *gl, double *dPdrhoetstar, spec_t dPdrhok, dim_t dPdrhoV);

  void find_dT_dx(np_t np, gl_t *gl, double *dTdrhoetstar, spec_t dTdrhok, dim_t dTdrhoV);

  void find_dT_dU(np_t np, gl_t *gl, flux_t dTdU);
    
  double _rhok (np_t np, long spec);

  double _V(np_t np, long theta);

  void find_V(np_t np, dim_t V);

  double _etstar (np_t np);

  void find_w(np_t np, spec_t w);

  void reformat_w(gl_t *gl, spec_t w, char *suffix, bool *flag);

  void reformat_rhok(gl_t *gl, spec_t rhok, char *suffix, bool *flag);

  void reformat_T(gl_t *gl, double *T, char *suffix, bool *flag);

  void reformat_P(gl_t *gl, double *P, char *suffix,  bool *flag);

  double _T (np_t np, gl_t *gl);
  
  double _w(np_t np, long spec);

  double _eta (np_t *np, long l, gl_t *gl);

  double _alpha(np_t np, gl_t *gl, long theta, long vartheta);

  double _nu(np_t np, gl_t *gl, long spec);

  double _kappa(np_t *np, long l, gl_t *gl);

  void find_U(np_t *np, long l, gl_t *gl, spec_t rhok, dim_t V, double T);

  void find_U_2(np_t *np, long l, gl_t *gl, spec_t w, dim_t V, double P, double T);

  void find_U_3(np_t *np, long l, gl_t *gl, spec_t rhok, dim_t V, double Pstar);

  void add_dUstar_to_U(np_t *np, long l, gl_t *gl, flux_t dUstar);
  
  void find_prim_fluid(np_t *np, long l, gl_t *gl);

#endif
