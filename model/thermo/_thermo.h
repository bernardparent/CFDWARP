#ifndef _THERMO_H
#define _THERMO_H


#include <src/common.h>
#include <stdio.h>

#define calR 8.314472e0  /* J/(K mol) */
#define calA 6.02214199E23  /* particules per mole */
#define Rchem 1.987192004e0  /*cal/(K mol)*/
#define echarge 1.602176462e-19 /* C */
#define emass 9.10938188E-31  /* kg */
#define epsilon0 8.854187817E-12 /* permittivity of free space */
#define kB 1.3806488E-23   /* m2*kg/(s2*K) */


/* taken from - Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
              - Reviews of Modern Physics , Vol.72, No. 2, 2000
	      - P. J. Mohr AND B. N. Taylor, CODATA Recommended Values of the Fundamental Physical
	        constants, 1998 */

void read_model_thermo(void);

double _calM(long spec);

double _m(long spec);

double _C(long spec);

long _Charge_number(long spec);

double _s(long spec);

double _cp_from_w_T_equilibrium(spec_t w, double T);

double _cp_from_w_T_neutrals_equilibrium(spec_t w, double T);

double _cpk_from_T_equilibrium(long spec, double T);

double _R(spec_t w);

double _a_from_w_T_equilibrium(spec_t w,  double T);

double _e_from_w_T(spec_t w, double T);

double _h_from_w_T(spec_t w, double T);

double _hk_from_T(long spec, double T);

double _hk_from_T_equilibrium(long spec, double T);

double _P_from_w_rho_T(spec_t w, double rho, double T);

double _P_from_w_rho_T_Te(spec_t w, double rho, double T, double Te);

double _rho_from_w_P_T(spec_t w, double P, double T);

double _s_from_w_T(spec_t w, double T);

double _sk_from_T(long spec, double T);

double _sk_from_T_equilibrium(long spec, double T);

double _T_from_w_e(spec_t w, double e);

double _T_from_w_h(spec_t w, double h);

double _T_from_rhok_P(spec_t rhok, double P);

double _T_from_w_rho_P(spec_t w, double rho, double P);



		     
double _de_dP_at_constant_rho(spec_t rhok, double T);

double _de_dT_at_constant_rho(spec_t rhok, double T);

void find_de_drhok_at_constant_P(spec_t rhok, double T, spec_t dedrhok);

void find_de_drhok_at_constant_T(spec_t rhok, double T, spec_t dedrhok);

void find_species_name(long spec, char **name);

int find_neutral_spec_from_ion_spec(long specion, long *specneutral);

void find_species_variable_name(long spec, char **name);

double _dsk_dT_from_T(long spec, double T);

double _dsk_dT_from_T_equilibrium(long spec, double T);


/* electron energy */

double _Te_from_ee(double ee);

double _ee_from_Te(double Te);

double _dee_dTe_from_Te(double Te);

double _he_from_Te(double Te);


/* nitrogen vibration energy */
double _Tv_from_ev(double ev);

double _ev_from_T(double T);

double _dev_dTv_from_Tv(double T);

double _zetav_from_Te_Parent2024(double Te);

double _zetav_from_Te_Aleksandrov(double Te);

double _zetae_from_Te(double Te);

double _dzetae_dTe_from_Te(double Te);

void find_Te_from_EoverN(double Estar, double *Te);

void find_dTe_dEoverN_from_EoverN(double Estar, double *dTedEstar);

void find_EoverN_from_Te(double Te, double *EoverN);

double _ionizationpot(long spec);

double _numatoms(long spec);

double _EoverN_from_rhok_Te(spec_t rhok, double Te);

double _Te_from_rhok_EoverN(spec_t rhok, double EoverN);

double _EoverNk_from_Te(long spec, double Te);


#endif /* _THERMO_H */
