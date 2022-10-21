#ifndef _FLUID_SHARE_H
#define _FLUID_SHARE_H

#include <model/_model.h>
#include <src/common.h>


#ifdef TEST
  #define TEST_VANISH 0.0e0
#else
  #define TEST_VANISH 1.0e0
#endif

void find_Sheatforces(np_t *np, gl_t *gl, long l, flux_t S);

void read_model_fluid_actions_Fbody_Qadd(char *actionname, char **argum, SOAP_codex_t *codex);

void find_chi_inverse(dim2_t X, dim2_t ChiInv);

#if defined(_FLUID_NEUTRALSTRANSPORT)

double _Vstar (np_t np, long theta);

double _Vhat (np_t np, long theta);

double _gamma(np_t np, gl_t *gl);

double _Rgas(np_t np, gl_t *gl);

void find_P_bdry_symmetrical(np_t *np, gl_t *gl,  long lA, long lB, long lC, long lD, long theta, long thetasgn, bool BDRYDIRECFOUND, int ACCURACY, double *P);

void find_Pstar_bdry_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, bool BDRYDIRECFOUND, int ACCURACY, double *Pstar);

void find_V_P_T_bdry_freestream(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, int ACCURACY, dim_t V, double *P, double *T, bool *OUTFLOW, bool *FREESTREAM);
#endif

#if defined(_FLUID_MULTISPECIES) && defined(_FLUID_NEUTRALSTRANSPORT)
void find_w_V_P_T_bdry_inflow_reservoir(np_t *np, gl_t *gl, long lA, long lB, long lC, int ACCURACY, spec_t w, dim_t V, double *P, double *T);

void find_w_V_P_T_bdry_inflow_reservoir_2(np_t *np, gl_t *gl, long lA, long lB, long lC, int ACCURACY,  spec_t w, dim_t V, double *P, double *T);

void find_w_V_P_T_bdry_symmetrical(np_t *np, gl_t *gl,  long lA, long lB, long lC, long lD, long theta, long thetasgn, bool BDRYDIRECFOUND, int ACCURACY, spec_t w, dim_t V, double *P, double *T);



#endif






#if _FLUID_N2VIBMODEL
double _kappav (np_t *np, long l, gl_t *gl);
#endif


#ifdef _FLUID_FAVREREYNOLDS

#define TURBMODEL_KEPSILON 1
#define TURBMODEL_KOMEGA1988 2
#define TURBMODEL_KOMEGA2008 3
#define DILATDISSIP_NONE 0
#define DILATDISSIP_WILCOX 1
#define DILATDISSIP_SARKAR 2

double _ktilde (np_t *np, long l, gl_t *gl);

double _psitilde (np_t np, gl_t *gl);

double _eps (np_t np, gl_t *gl);

double _omega (np_t *np, long l, gl_t *gl);

double _dVi_dxj(np_t *np, long l, gl_t *gl, long i, long j);

double _etat_from_rho_eta_k_psitilde(np_t *np, long l, gl_t *gl, double rho, double eta, double k, double psitilde);

double _etat_mem(np_t *np, long l, gl_t *gl);

double _kappastar (np_t *np, long l, gl_t *gl);

#if _FLUID_N2VIBMODEL
double _kappavstar (np_t *np, long l, gl_t *gl);
#endif

double _nustar (np_t *np, long l, gl_t *gl, long spec);

double _etakstar (np_t *np, long l, gl_t *gl);

double _etapsistar (np_t *np, long l, gl_t *gl);

void find_k_psi_bdry_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, double *k, double *psi);

double _Qk(np_t *np, long l, gl_t *gl);

double _chiw(np_t *np, long l, gl_t *gl);

void find_Stnorm(np_t *np, gl_t *gl, long l, flux_t St);

void find_Stcomp(np_t *np, gl_t *gl, long l, flux_t St);

void find_dStnorm_dU(np_t *np, gl_t *gl, long l, sqmat_t dStnormdU);

void find_dStcomp_dU(np_t *np, gl_t *gl, long l, sqmat_t dStcompdU);


#endif


void find_conditioned_Lambda_absolute_from_jacvars_denom(jacvars_t jacvars, metrics_t metrics, sqmat_t Lambda);

void find_conditioned_Lambda_absolute_from_jacvars(jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t Lambdaabs);

void condition_Lambda_plus_minus(np_t *np, gl_t *gl, long lp0, long theta, jacvars_t jacvarsp0, jacvars_t jacvarsp1, metrics_t metrics,  int EIGENVALCOND, sqmat_t Lambdaplus, sqmat_t Lambdaminus);

void rearrange_metrics_eigenset2(metrics_t *metrics);

#if _FLUID_CONVECTION
void set_jacvars_eigenconditioning_constants(gl_t *gl, jacvars_t *jacvars);
#endif

#ifdef _FLUID_PLASMA
double _betan(long k);

double _betac(long k);

double _betag(gl_t *gl, long k);

double _betaa(gl_t *gl, long k);

double _betai(long k);

double _betaplus(long k);

double _betaminus(long k);

double _rhoc(np_t np, gl_t *gl);

bool is_bdry_cathode(np_t *np, gl_t *gl, long lA, long lB, long lC,
                     long theta, long thetasgn);

double _sigma_positive(np_t *np, gl_t *gl, long l);

double _alpha(np_t *np, gl_t *gl, long l, long k, long r);

void find_alpha(np_t *np, gl_t *gl, long l, spec2_t alpha);

void find_Z(np_t *np, gl_t *gl, long l, sqmat_t Z);

void find_dZU_dU(np_t *np, gl_t *gl, long l, sqmat_t dZU_dU);

void find_LambdaZ(np_t *np, gl_t *gl, long l, sqmat_t LambdaZ);


#endif

void find_init_mass_fraction_templates(char **specstr1, char **specstr2);

void find_init_molar_fraction_templates(char **specstr1, char **specstr2);

void find_init_number_density_templates(char **specstr1, char **specstr2);



void find_Saxi(np_t *np, gl_t *gl, long l, flux_t S);

void find_dSaxi_dU(np_t *np, gl_t *gl, long l, sqmat_t dS_dU);

double _w_product_at_catalytic_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, long specR, long specP, double factprodreact);


void update_w_at_catalytic_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, double Twall, double Tewall, long paramstart, long paramend, spec_t wwall);

void update_w_V_at_injection_wall(np_t *np, gl_t *gl, long lA, long lB, long lC, double Twall, double Tewall, 
                                  long paramstart, long paramend, spec_t wwall, dim_t Vwall);


#endif /* _FLUID_SHARE_H */
