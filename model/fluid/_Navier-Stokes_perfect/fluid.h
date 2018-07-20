#ifndef FLUID_H

  #define FLUID_H
  #include <src/common.h>
  #include <model/fluid/_fluid.h>


  double _P (np_t np, gl_t *gl);

  double _q (np_t np);

  double _a (np_t np, gl_t *gl);

  double _rho (np_t np);  
    
  double _V(np_t np, long theta);

  void find_V(np_t np, dim_t V);

  double _T (np_t np, gl_t *gl);
  
  double _et (np_t np);

  double _ht (np_t np, gl_t *gl);

  void find_U(np_t *np, long l, gl_t *gl, double rho, dim_t V, double T);

  void find_U_2(np_t *np, long l, gl_t *gl, dim_t V, double P, double T);

  void find_U_3(np_t *np, long l, gl_t *gl, double rho, dim_t V, double P);

  void add_dUstar_to_U(np_t *np, long l, gl_t *gl, flux_t dUstar);
  
  void find_prim_fluid(np_t *np, long l, gl_t *gl);

  double _dVi_dxj(np_t *np, long l, gl_t *gl, long i, long j);

  void reformat_rho(gl_t *gl, double *rho, char *suffix,  bool *flag);

  void reformat_T(gl_t *gl, double *T, char *suffix, bool *flag);

  void reformat_P(gl_t *gl, double *P, char *suffix, bool *flag);


#endif
