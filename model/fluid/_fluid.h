#ifndef _FLUID_H

  #define _FLUID_H
  

  //void add_to_clipped_variables(gl_t *gl, char *str);

  /* needed by emfield module */
  //void find_metrics_from_base_level(np_t *np, gl_t *gl, long l, metrics_t *metrics);

 
  /* total number density */
  double _N(np_t np, gl_t *gl);
  
  double _Nk(np_t np, gl_t *gl, long spec);
  
  double _rhok (np_t np, long spec);

  double _Nk_from_U (flux_t U, long spec);

  double _Tk(np_t *np, gl_t *gl, long l, long spec);
  
  double _Pk(np_t *np, gl_t *gl, long l, long spec);

  double _rho (np_t np);
  
  double _T (np_t np, gl_t *gl);
  
  double _P (np_t np, gl_t *gl);

  void find_w(np_t np, spec_t w);
  
  double _rhoc(np_t np, gl_t *gl);
  
  /* neutral number density */
  double _Nn(np_t np, gl_t *gl);
 
  double _V(np_t np, long theta);
    
  void find_V(np_t np, dim_t V);
    
  void read_model_fluid_actions(char *actionname, char **argum, SOAP_codex_t *codex);
  
  void read_cycle_fluid_actions(char *actionname, char **argum, SOAP_codex_t *codex);

  void read_disc_fluid_actions(char *actionname, char **argum, SOAP_codex_t *codex);

  void write_model_fluid_template(FILE **controlfile);

  void write_cycle_fluid_template(FILE **controlfile);

  void write_disc_fluid_template(FILE **controlfile);

  void write_init_fluid_template(FILE **controlfile);

  void find_post_variable_name_fluid(long varnum, char *varname);

  void find_post_variable_value_fluid(np_t *np, long l, gl_t *gl, long varnum, double *varvalue);
  
  bool is_node_bdry_symmetry_plane_fluid(np_t np);

  double _Pe_from_jacvars(jacvars_t jacvars, metrics_t metrics);

  double _a(np_t np, gl_t *gl);

  double _a_from_jacvars(jacvars_t jacvars);

  double _eta_from_jacvars(jacvars_t jacvars);

  double _V_from_jacvars(jacvars_t jacvars, long dim);

  double _Pstar (np_t np, gl_t *gl);

  double _Re_edge_from_jacvars(jacvars_t jacvars, metrics_t metrics);

  double _Te_m1(np_t np, gl_t *gl);

  double _eta(np_t np, gl_t *gl);

  double _V_from_U(np_t np, long theta);

  double _kappa(np_t np, gl_t *gl);

  double _eta(np_t np, gl_t *gl);

  double _nu(np_t np, gl_t *gl, long spec);


#ifdef _FLUID_FAVREREYNOLDS

  double _dVi_dxj(np_t *np, long l, gl_t *gl, long i, long j);

  double _etastar (np_t *np, long l, gl_t *gl);

  double _psitilde (np_t np, gl_t *gl);

  double _k (np_t np);

  double _etat (np_t *np, long l, gl_t *gl);

  double _psi (np_t np);

  double _ktilde (np_t np, gl_t *gl);

  double _eps (np_t np, gl_t *gl);

  double _Vhat (np_t np, long theta);

  double _athermo(np_t np, gl_t *gl);

  double _nustar (np_t *np, long l, gl_t *gl, long spec);

  double _etakstar (np_t *np, long l, gl_t *gl);

  double _kappastar (np_t *np, long l, gl_t *gl);

  double _etapsistar (np_t *np, long l, gl_t *gl);

  
#endif

#if _FLUID_N2VIBMODEL
  double _Tv(np_t np);
#endif

#ifdef _FLUID_PLASMA
  void find_alpha(np_t *np, gl_t *gl, long l, spec2_t alpha);

#endif

#endif
