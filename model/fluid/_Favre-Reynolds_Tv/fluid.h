#ifndef FLUID_H

  #define FLUID_H
  #include <src/common.h>
  #include <model/fluid/_fluid.h>

  #define N2VIBMODEL_MACHERET 1
  #define N2VIBMODEL_MILLIKAN 2
  
  bool _FLUIDPRIMMEM(np_t np);

  double _P (np_t np, gl_t *gl);

  double _q (np_t np);

  double _a (np_t np, gl_t *gl);

  double _rho (np_t np);

  double _k (np_t np);

  double _psi (np_t np);

  double _Tv(np_t np);
  
  double _omega (np_t np, gl_t *gl);

  double _athermo(np_t np, gl_t *gl);

  void find_dP_dx(np_t np, gl_t *gl, double *dPdrhoetstar, spec_t dPdrhok, dim_t dPdrhoV);

  void find_dT_dx(np_t np, gl_t *gl, double *dTdrhoetstar, spec_t dTdrhok, dim_t dTdrhoV);

  void find_dT_dU(np_t np, gl_t *gl, flux_t dTdU);
  
  void find_dTv_dU(np_t np, gl_t *gl, flux_t dTvdU);
  
  double _rhok (np_t np, long spec);

  double _V(np_t np, long theta);

  void find_V(np_t np, dim_t V);

  double _etstar (np_t np);

  double _Pstar (np_t np, gl_t *gl);

  void find_w(np_t np, spec_t w);

  void reformat_w(gl_t *gl, spec_t w, char *suffix,  bool *flag);

  void reformat_rhok(gl_t *gl, spec_t rhok, char *suffix,  bool *flag);

  void reformat_k_psi(gl_t *gl, double *k, double *psi, char *suffix,  bool *flag);

  void reformat_T(gl_t *gl, double *T, char *suffix,  bool *flag);

  void reformat_Tv(gl_t *gl, double T, double *Tv, char *suffix,  bool *flag);

  void reformat_P(gl_t *gl, double *P, char *suffix,  bool *flag);

  double _T (np_t np, gl_t *gl);
  
  double _w(np_t np, long spec);

  double _etastar (np_t *np, long l, gl_t *gl);

  double _psitilde (np_t np, gl_t *gl);

  double _ktilde (np_t np, gl_t *gl);

  double _ev (np_t np);

  double _evzero (np_t np, gl_t *gl);

  double _Tv_from_ev(double ev);

  double _eps (np_t np, gl_t *gl);

  double _eta (np_t np, gl_t *gl);

  double _etat (np_t *np, long l, gl_t *gl);

  double _nustar (np_t *np, long l, gl_t *gl, long spec);

  double _etakstar (np_t *np, long l, gl_t *gl);

  double _etapsistar (np_t *np, long l, gl_t *gl);

  double _kappastar (np_t *np, long l, gl_t *gl);

  void find_U(np_t *np, long l, gl_t *gl, spec_t rhok, dim_t V, double T, double k, double psi, double Tv);

  void find_U_2(np_t *np, long l, gl_t *gl, spec_t w, dim_t V, double P, double T, double k, double psi, double Tv);

  void find_U_3(np_t *np, long l, gl_t *gl, spec_t rhok, dim_t V, double Pstar, double k, double psi, double Tv);

  void add_dUstar_to_U(np_t *np, long l, gl_t *gl, flux_t dUstar);
  
  void find_prim_fluid(np_t *np, long l, gl_t *gl);

  double _dVi_dxj(np_t *np, long l, gl_t *gl, long i, long j);

#endif
