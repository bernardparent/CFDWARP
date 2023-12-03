#ifndef _RES_SHARE_H
#define _RES_SHARE_H

#include <model/_model.h>
#include <src/common.h>
#include <cycle/share/cycle_share.h>

#define LIMITER_FIRSTORDER 0
#define LIMITER_MINMOD 1
#define LIMITER_MINMOD2 2
#define LIMITER_VANLEER 3
#define LIMITER_SUPERBEE 4
#define LIMITER_VENKATAKRISHNAN 5
#define LIMITER_KOREN 6
#define LIMITER_VANALBADA 7
#define LIMITER_OSPRE 8
#define LIMITER_OSHER 9
#define LIMITER_SWEBY 10
#define LIMITER_SMART 11

#define numiter_FVSplus_default 2
#define numiter_FDSplus_default 2

#define AOWENO_TYPE_DIFFUSIVE 0
#define AOWENO_TYPE_COMPRESSIVE 1

#ifdef _FLUID_NEUTRALSTRANSPORT
  #define EIGENVALCOND_HARTEN 0
  #define EIGENVALCOND_GNOFFO 1
  #define EIGENVALCOND_PECLET 2
  #define EIGENVALCOND_PASCAL 3
  #define EIGENVALCOND_PARENT 4
#endif

#define EIGENVALCOND_NONE -1
#define EIGENVALCOND_DEFAULT -2

#define AVERAGING_ARITH 0
#define AVERAGING_ROE 1



void find_musclvarscycle_offset(np_t *np, gl_t *gl, long l, long theta, long offset, musclvarscycle_t musclvars);

double _f_AOWENO5(double um2, double um1, double up0, double up1, 
                   double up2, double x, double gammalo, double gammahi, int AOWENO_TYPE);

double _f_AOWENO7(double um3, double um2, double um1, double up0, double up1, 
                   double up2, double up3, double x, double gammalo, double gammahi, int AOWENO_TYPE);

double _f_AOWENO9(double um4, double um3, double um2, double um1, double up0, 
                   double up1, double up2, double up3, double up4, double x, double gammalo, double gammahi, int AOWENO_TYPE);

double _f_CWENO3(double um2, double um1, double up0, double up1, double up2, double x);

double _f_CWENO5(double um4, double um3, double um2, double um1, double up0, double up1, 
                     double up2, double up3, double up4, double x);

double _f_CWENO6(double um5, double um4, double um3, double um2, double um1, double up0, double up1, 
                     double up2, double up3, double up4, double up5, double x);

double _f_WENO9(double um4, double um3, double um2, double um1, double up0, double up1, double up2, double up3, double up4);

double _f_WENO7(double um3, double um2, double um1, double up0, double up1, double up2, double up3);

double _f_WENO5(double um2, double um1, double up0, double up1, double up2);

double _f_WENO3(double um1, double up0, double up1);

double _limiter_TVD(double r, int LIMITER);

double _f_TVD5(double um2, double um1, double up0, double up1, double up2);

double _f_TVD4(double um2, double um1, double up0, double up1);

double _f_TVD3(double um1, double up0, double up1);

double _f_TVD2(double um1, double up0, double up1, int LIMITER);

double _f_TVD_AO5(double um2, double um1, double up0, double up1, double up2);

double _f_WENO_TVD5(double um2, double um1, double up0, double up1, double up2);

void find_Lambda_minus_dtau_FVS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, sqmat_t lambdaminus);

void find_Lambda_plus_dtau_FVS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, sqmat_t lambdaplus);

void find_Lambda_minus_dtau_FDS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, int AVERAGING, sqmat_t lambdaminus);

void find_Lambda_plus_dtau_FDS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, int AVERAGING, sqmat_t lambdaplus);

void find_Lambda_minus_dtau_FDS_from_jacvars(jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t lambdaminus);

void find_Lambda_plus_dtau_FDS_from_jacvars(jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t lambdaplus);

void find_Lambda_plus_minus_dtau_FDS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, int AVERAGING, sqmat_t lambdaplus, sqmat_t lambdaminus);

void find_Fstar_interface_FDSplus_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h, long theta, 
                     flux_t musclvarsm1h, flux_t musclvarsp1h, metrics_t metrics,  long numiter, 
                     int EIGENVALCOND, int AVERAGING, flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h);

void find_Fstar_interface_FDS_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h,  long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, int AVERAGING, flux_t Fint);

void find_Fstar_interface_FDSR_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h,  long theta, flux_t musclvarsm1h, 
                     flux_t musclvarsp1h, metrics_t metrics, int EIGENVALCOND, int AVERAGING, flux_t Fint);

void find_Fstar_interface_FVSplus_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h, long theta, 
                     flux_t musclvarsm1h, flux_t musclvarsp1h, metrics_t metrics,  long numiter, 
                     int EIGENVALCOND, int AVERAGING, flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h);

void find_Fstar_interface_FVS_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h,  long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, int AVERAGING, flux_t Fint);

void find_jacvars_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, int AVERAGING, jacvars_t *jacvars);

void find_jacvars_at_interface_from_jacvars(jacvars_t jacvarsL, jacvars_t jacvarsR, gl_t *gl, long theta, metrics_t metrics, int AVERAGING, jacvars_t *jacvars);

void filter_Fstar_interface_positivity_preserving_PARENT(np_t *np, gl_t *gl, long lp0, long theta, metrics_t metrics, long numiter, int EIGENVALCOND, flux_t Fint, flux_t Fpositive, sqmat_t lambdaminus, sqmat_t lambdaplus);

void filter_Fstar_interface_positivity_preserving_TEST(np_t *np, gl_t *gl, long lp0, long theta, metrics_t metrics, long numiter, int EIGENVALCOND, flux_t Fint, flux_t Fpositive, sqmat_t lambdaminus, sqmat_t lambdaplus);



#endif /* _RES_SHARE_H */
