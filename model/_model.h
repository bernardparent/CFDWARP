#ifndef _MODEL_H
#define _MODEL_H

#include <src/common.h>

#define FLUIDBDRYINTERFACE_IN 0
#define FLUIDBDRYINTERFACE_OUT 1
#define FLUIDBDRYINTERFACE_IN_OUT 2


long max_order_accuracy_near_discontinuities(long musclvarflux);

void find_clipped_variables_list(gl_t *gl, char **cliplist);

void find_clipped_muscl_variables_list(gl_t *gl, char **cliplist);

void find_clipped_bdry_variables_list(gl_t *gl, char **cliplist);
 
bool is_node_bdry_symmetry_plane(np_t np);

long _time_marching_method_emfield(long fluxe, gl_t *gl);

void write_init_fluid_template(FILE **controlfile);

void write_bdry_fluid_template(FILE **controlfile);

void write_cycle_fluid_template(FILE **controlfile);

void write_model_template(FILE **controlfile);

void read_model(char *argum, SOAP_codex_t *codex);

void add_bdry_types_fluid_to_codex(SOAP_codex_t *codex);

void add_bdry_types_emfield_to_codex(SOAP_codex_t *codex);

void add_init_types_fluid_to_codex(SOAP_codex_t *codex);

void add_init_types_emfield_to_codex(SOAP_codex_t *codex);

void reset_clipped_variables(gl_t *gl);

void free_clipped_variables(gl_t *gl);

double _Omega(np_t np, gl_t *gl);

double _X(np_t np, long dim1, long dim2);

double _varphi(np_t *np, gl_t *gl, long l, long dim);

void find_dFstar_dUstar(np_t np, gl_t *gl, long theta, sqmat_t A);

void find_dG_dUstar(np_t np, gl_t *gl, sqmat_t B);

void find_dSstar_dUstar(np_t *np, gl_t *gl, long l, sqmat_t C);

void find_Fstar(np_t np, gl_t *gl, long theta, flux_t F);

void find_Fstar_given_metrics(np_t np, gl_t *gl, metrics_t metrics, long theta, flux_t Fstar);

void find_G(np_t np, gl_t *gl, flux_t G);

void find_H1(np_t *np, gl_t *gl, long l, flux_t H);

void find_H2(np_t *np, gl_t *gl, long l, flux_t H);

void find_dH1_dUstar(np_t *np, gl_t *gl, long l, sqmat_t dHdUstar);

void find_dH2_dUstar(np_t *np, gl_t *gl, long l, sqmat_t dHdUstar);

void find_Y1star_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, flux_t Ystar);

void find_Y2star_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, flux_t Ystar);

void find_metrics_at_interface(np_t *np, gl_t *gl, long lL, long lR,
                              long theta, metrics_t *metrics);

void find_Kstar_interface(np_t *np, gl_t *gl, long lL, long lR, metrics_t metrics, long theta, long vartheta, sqmat_t K, int CYCLELEVEL);

void find_Ustar(np_t np, gl_t *gl, flux_t Ustar);

void find_Sstar(np_t *np, gl_t *gl, long l, flux_t S);

void find_Dstar(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstar);

void find_Dstarplus(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstarplus);

void find_Dstarminus(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstarminus);

void find_D2star_interface(np_t *np, gl_t *gl, long l, long theta, flux_t Dstar);

void find_Dstar_test(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstar);

void find_Z(np_t *np, gl_t *gl, long l, sqmat_t Z);

void find_dZU_dU(np_t *np, gl_t *gl, long l, sqmat_t dZU_dU);

void find_LambdaZ(np_t *np, gl_t *gl, long l, sqmat_t LambdaZ);

void find_J(np_t *np, gl_t *gl, long l, EXM_vec3D_t J);

void add_dUstar_to_U(np_t *np, long l, gl_t *gl, flux_t dUstar);

void find_prim_fluid(np_t *np, long l, gl_t *gl);

void update_metrics_at_node(np_t *np, gl_t *gl, long l);

void find_metrics_from_base_level(np_t *np, gl_t *gl, long l, metrics_t *metrics);

double _x(np_t np, long dim);

void init_node_fluid(np_t *np, long l, gl_t *gl, long inittype, initvar_t values);

void update_bdry_fluid(np_t *np, gl_t *gl, long lA, long lB, long lC, long lD, long theta, long thetasgn, bool BDRYDIRECFOUND, int TYPELEVEL);

void find_metrics_at_interface(np_t *np, gl_t *gl, long lL, long lR,
                               long theta, metrics_t *metrics);

void find_metrics_at_interface_2(np_t *np, gl_t *gl, long lL, long lR,
                               long theta, long theta2, metrics_t *metrics);

void find_metrics_at_node(np_t *np, gl_t *gl, long l,
                               long theta, metrics_t *metrics);

void find_Fstar_from_jacvars(jacvars_t jacvars, metrics_t metrics, flux_t F);

void find_jacvars_from_U(flux_t U, metrics_t metrics, gl_t *gl, long theta, jacvars_t *jacvars);

void find_A_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t A);

void find_Linv_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t R);

void find_L_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t L);

void find_LUstar_from_jacvars(jacvars_t jacvars,  metrics_t metrics, flux_t LUstar);

void find_Ustar_from_jacvars(jacvars_t jacvars,  metrics_t metrics, flux_t Ustar);

void find_Lambda_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t lambda);

void find_jacvars_at_interface_Roe_average(jacvars_t jacvarsL, jacvars_t jacvarsR, gl_t *gl, long theta, jacvars_t *jacvars);

void find_jacvars_at_interface_arith_average(jacvars_t jacvarsL, jacvars_t jacvarsR, gl_t *gl, long theta, jacvars_t *jacvars);
 
void find_jacvars(np_t np, gl_t *gl, metrics_t metrics, long theta, jacvars_t *jacvars);

void find_musclvars(np_t np, gl_t *gl, flux_t musclvars);

bool is_musclvar_a_charged_mass_fraction(long flux);

bool is_musclvar_a_mass_fraction(long flux);

void find_jacvars_from_musclvars(flux_t musclvars, metrics_t metrics, gl_t *gl, long theta, jacvars_t *jacvars);

void find_musclvars_from_jacvars(jacvars_t jacvars, gl_t *gl, flux_t musclvars);

void find_Ustar_from_musclvars(flux_t musclvars, metrics_t metrics, gl_t *gl, flux_t Ustar);


void find_Ustar_given_metrics(np_t np, gl_t *gl, metrics_t metrics, flux_t Ustar);


int reversibly_expand_q_and_rho_to_given_P(np_t np, gl_t *gl, double Pc, double *qc, double *rhoc,
           long numsteps, double q_min);

double _rho (np_t np);

double _Pstar (np_t np, gl_t *gl);

double _T (np_t np, gl_t *gl);

double _htstar (np_t np, gl_t *gl);

double _V(np_t np, long theta);

double _w(np_t np, long spec);

double _Pstag(np_t np, gl_t *gl, long numsteps);

double _Tstag(np_t np, gl_t *gl);

void find_post_proc_var_name(long varnum, char *varname);

void find_post_proc_var_value(np_t *np, long l, gl_t *gl, long varnum, double *varvalue);

double _a(np_t np, gl_t *gl);


double _T(np_t np, gl_t *gl);

void find_default_initvar(np_t *np, gl_t *gl, long l, initvar_t initvar);


#ifdef EMFIELD

void write_init_emfield_template(FILE **controlfile);

void write_bdry_emfield_template(FILE **controlfile);

void write_cycle_emfield_template(FILE **controlfile);

void update_bdry_emfield(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, bool BDRYDIRECFOUND, int TYPELEVEL);

void update_Te_local(np_t *np, gl_t *gl, long l);

void add_dUstar_to_U_emfield(np_t *np, gl_t *gl, long l, fluxemfield_t dUstar);

void find_prim_emfield_mem_1(np_t *np, gl_t *gl, long l);

void find_prim_emfield_mem_2(np_t *np, gl_t *gl, long l);

void find_prim_emfield_mem_3(np_t *np, gl_t *gl, long l);

void find_Sstar_emfield(np_t *np, gl_t *gl, long l, fluxemfield_t S);

void find_Fstar_interface_emfield(np_t *np, gl_t *gl, long lL, long lR, long theta, fluxemfield_t F);

double _U_emfield(np_t np, gl_t *gl, long flux);

double _xi_emfield(np_t np, gl_t *gl, fluxemfield_t Res);

void find_dSstar_dUstar_emfield(np_t *np, gl_t *gl, long l, long flux, double *C);

void init_node_emfield(np_t np, gl_t *gl, long inittype, initvar_emfield_t values);

void find_default_initvar_emfield(np_t *np, gl_t *gl, long l, initvar_emfield_t initvar);

void find_dFstar_dUstar_interface_emfield(np_t *np, gl_t *gl, long lL, long lR, long i, long flux, 
                                       double *dFstardUstarLL, double *dFstardUstarL, 
                                       double *dFstardUstarR, double *dFstardUstarRR);

void find_dtau_emfield(np_t *np, gl_t *gl, long l, long flux, double *dtau);

void find_linearization_coefficients_bdry_node_emfield(np_t *np, gl_t *gl, long lA, long theta, long thetasgn,
                        long flux, long bdrytype, double *valA, double *valB, double *valRHS);

void find_linearization_coefficients_inner_node_emfield(np_t *np, gl_t *gl, long l, long theta, long flux, double *valC, double *valB, double *valA);

#endif

void add_to_clipped_variables(gl_t *gl, char *str);

void add_to_clipped_variables2(gl_t *gl, char *str, char *suffix);

#ifdef DISTMPI
void find_clipped_variables_list_all(gl_t *gl, char **cliplist);

void find_clipped_muscl_variables_list_all(gl_t *gl, char **cliplist);

void find_clipped_bdry_variables_list_all(gl_t *gl, char **cliplist);

void add_to_clipped_variables_all(gl_t *gl, char *str, long clipnum);

void reset_clipped_variables_all(gl_t *gl);
 
#endif

void find_species_name(long spec, char **name);

double _Vstar (np_t np, long theta);

void test_dSchem_dU(np_t *np, gl_t *gl, long l);

void find_emfield_force(np_t *np, gl_t *gl, long l, dim_t Femfield);

double _E_dot_J(np_t *np, gl_t *gl, long l);

double _E_dot_J_recast(np_t *np, gl_t *gl, long l);

double _Qbeam(np_t np, gl_t *gl);

void find_Gamma(np_t *np, gl_t *gl, long l, sqmat_t Gamma);

void find_Gamma_inverse(np_t *np, gl_t *gl, long l, sqmat_t GammaInv);

double _flux_norm(long flux);

void condition_Lambda_plus_minus(np_t *np, gl_t *gl, long lp0, long theta, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  int EIGENVALCOND, sqmat_t Lambdaplus, sqmat_t Lambdaminus);

void find_conditioned_Lambda_absolute_from_jacvars(jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t Lambda);

void find_conditioned_Lambda_absolute_from_jacvars_denom(jacvars_t jacvars, metrics_t metrics, sqmat_t Lambda);

void find_Uprime(np_t np, gl_t *gl, flux_t Uprime);

void find_Uprime_from_jacvars(jacvars_t jacvars, flux_t Uprime);

void find_Uprime_from_U(flux_t U, flux_t Uprime);

void find_dUstar_dUprime_from_jacvars(jacvars_t jacvars, metrics_t metrics, sqmat_t dUstardUprime);

bool is_node_bdry_no_cross(np_t np, int TYPELEVEL);


#endif /* _MODEL_H */
