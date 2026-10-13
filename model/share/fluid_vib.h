#ifndef _FLUID_VIB_H
#define _FLUID_VIB_H

#include <model/_model.h>
#include <src/common.h>
#include <model/thermo/_thermo.h>


#if _FLUID_VIBMODEL
/* model of the V-T relaxation time */
#define TAUVTMODEL_PARENT2024 1
#define TAUVTMODEL_MILLIKANWHITE 2
#define TAUVTMODEL_TEQUILIBRIUM 3

/* model of the e-V excitation */
#define QEVMODEL_NONE 0
#define QEVMODEL_ALEKSANDROV 1
#define QEVMODEL_PARENT2024 2
#define QEVMODEL_RODRIGUEZ2025 3
#define QEVMODEL_RODRIGUEZ2026 4

/* deprecated; defined only so that a deprecation message can be given */
#define N2VIBMODEL_ALEKSANDROV 101
#define N2VIBMODEL_TEQUILIBRIUM 102
#define N2VIBMODEL_PARENT2024 103
#define N2VIBMODEL_RODRIGUEZ2025 104
#define N2VIBMODEL_RODRIGUEZ2026 105
#define N2VIBMODEL_MILLIKAN 106

double _kappav (np_t *np, long l, gl_t *gl);

double _nuv (np_t *np, long l, gl_t *gl);

void reformat_ev(gl_t *gl, spec_t w, double *ev, char *suffix, bool *flag);

void find_dev_dU(np_t np, gl_t *gl, flux_t devdU);

double _wv(np_t np);

double _ev (np_t np);

double _evzero (np_t np, gl_t *gl);

void find_vibdist_from_np(np_t np, gl_t *gl, vibdist_t *vibdist);

double _Tv(np_t np, gl_t *gl);

double _Tvk(np_t np, gl_t *gl, long spec);

double _Tv_from_w_T_ev(gl_t *gl, spec_t w, double T, double ev);

double _Tv_from_np_ev(np_t np, gl_t *gl, double ev);

void find_dTv_dU(np_t np, gl_t *gl, flux_t dTvdU);

void read_model_fluid_actions_vib_species(char *actionname, char **argum, SOAP_codex_t *codex);

void add_vib_constants_to_codex(SOAP_codex_t *codex);

void read_model_fluid_vib(SOAP_codex_t *codex);

void find_Svib_VT(np_t *np, gl_t *gl, long l, flux_t S);

void find_dSvib_VT_dU(np_t *np, gl_t *gl, long l, flux_t dTdU, sqmat_t dSvibdU);

void test_dSvib_dU_display(gl_t *gl, sqmat_t dSvibdU, sqmat_t dSvibdUnum);

#ifdef _FLUID_PLASMA
double _Qev(np_t *np, gl_t *gl, long l, double emfieldfact);

void find_dQev_dU(np_t *np, gl_t *gl, long l, sqmat_t dQevdU);
#endif

#ifdef _FLUID_FAVREREYNOLDS
double _kappavstar (np_t *np, long l, gl_t *gl);

double _nuvstar (np_t *np, long l, gl_t *gl);
#endif

#endif

#endif /* _FLUID_VIB_H */
