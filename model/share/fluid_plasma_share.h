#ifndef _FLUID_PLASMA_SHARE_H
#define _FLUID_PLASMA_SHARE_H

#include <model/_model.h>
#include <src/common.h>

#ifdef _FLUID_PLASMA
double _betan(long k);

static inline double _betac(long k){
  double betac;
  if (speciestype[k]!=SPECIES_NEUTRAL) betac=1.0; else betac=0.0;
  return(betac);
}

double _betag(gl_t *gl, long k);

double _betaa(gl_t *gl, long k);

double _betai(long k);

double _betaplus(long k);

bool is_bdry_cathode(np_t *np, gl_t *gl, long lA, long lB, long lC,
                     long theta, long thetasgn);

double _sigma_positive(np_t *np, gl_t *gl, long l);

double _alpha(np_t *np, gl_t *gl, long l, long k, long r);

double _Vk_drift_perp_to_bdry(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, long spec);

#endif

#if defined(_FLUID_PLASMA) && defined(_FLUID_NEUTRALSTRANSPORT)
double _Te_local(np_t *np, gl_t *gl, long l, double Temin);

double _Pkmem(np_t *np, gl_t *gl, long l, long spec);

double _Nnmem(np_t np, gl_t *gl);

double _Nmem(np_t np, gl_t *gl);

void find_dQbeam_dU(np_t np, gl_t *gl, flux_t dQbeamdU);

#endif

/* Ballistic electron model of parent2026, shared by the plasma fluid modules.
   It is compiled when the chem module defines _CHEM_BALLISTIC and the fluid module defines _FLUID_PLASMA.

   The chem module must provide (see model/chem/_chem.h):
     - nb ballistic groups, with the species ordered as the bulk electrons (speceminus=0), the groups 1 to nb,
       and then the heavy species; specegroup(j) is the species of group j
     - in gl->model.chem: BALLISTIC, beamdir, epsborn, kturnmax, phiemit and the potential windows phimid, Deltaphi and r
     - numQballistic, QB_BORN, QB_BEAM, QB_RUN, QB_THRESHOLD and find_Q_ballistic()

   The fluid module must provide:
     - an electron energy equation (fluxee), and in gl->model.fluid: chibmin, INCLUDE_E_EFFECTS_ON_ELECTRON_MOBILITY and
       the TEMODELs TEMODEL_CONSTANT, TEMODEL_TEQUILIBRIUM and TEMODEL_LOCAL (and TEMODEL_TVEQUILIBRIUM and
       TEMODEL_TRANSPORTEMF if defined); any other TEMODEL transports the electron energy
     - the numpostvarballistic post variables of the groups, placed at the end of those of the fluid module

   The transport module must give the ballistic groups no mobility, no thermal conductivity and no viscosity
   in find_nuk_eta_kappak_muk(), and no mobility in find_dmuk_from_rhok_Tk_Ek(). */

#if defined(_CHEM_BALLISTIC) && defined(_FLUID_PLASMA)

  /* potentials, energies and speeds of the groups */

  double _phibar(gl_t *gl, long j, double phi);

  double _eps_group(gl_t *gl, long j, double phi);

  bool _is_ballistic_group(long spec);

  /* current of the groups */
  double _Jstar_ballistic_groups_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long theta);

  /* control file */
  void read_model_fluid_actions_ballistic(SOAP_codex_t *codex, gl_t *gl);

  /* source terms */
  void add_Semfield_ballistic_groups(np_t *np, gl_t *gl, long l, flux_t S);

  void add_Schem_ballistic(np_t *np, gl_t *gl, long l, spec_t rhok, double T, double Te, flux_t S);

  void add_dSchem_dU_ballistic(np_t *np, gl_t *gl, long l, spec_t rhok, double T, double Te, flux_t dTdU, flux_t dTedU, sqmat_t dS_dU);

  /* boundary conditions */
  void update_Nk_electrons_at_cathode_ballistic(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, double gamma, spec_t Nk);

  void update_Nk_electrons_at_anode_ballistic(np_t *np, gl_t *gl, long lA, long lB, spec_t Nk);

  /* post variables */
  void find_post_variable_name_ballistic(long varnum, char *varname);

  void find_post_variable_value_ballistic(np_t *np, long l, gl_t *gl, long varnum, double *varvalue);

#endif

#endif /* _FLUID_PLASMA_SHARE_H */
