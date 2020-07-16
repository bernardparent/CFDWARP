#ifndef _THERMO_H
#define _THERMO_H


#include <src/common.h>
#include <stdio.h>

#define calR 8.314472e0
#define calA 6.02214199E23  
#define echarge 1.602176462e-19
#define emass 9.10938188E-31
#define epsilon0 8.854187817E-12 /* permittivity of free space */
#define kB 1.3806488E-23


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

double _cp_from_w_T(spec_t w, double T);

double _cpk_from_T(long spec, double T);

double _cpk_from_T_equilibrium(long spec, double T);

double _R(spec_t w);

double _a_from_w_T(spec_t w,  double T);

double _e_from_w_T(spec_t w, double T);

double _h_from_w_T(spec_t w, double T);

double _hk_from_T(long spec, double T);

double _hk_from_T_equilibrium(long spec, double T);

double _P_from_w_rho_T(spec_t w, double rho, double T);

double _P_from_w_rho_T_Te(spec_t w, double rho, double T, double Te);

double _rho_from_w_P_T(spec_t w, double P, double T);

double _s_from_w_T(spec_t w, double T);

double _dTskdT_from_T(long spec, double T);

double _sk_from_T(long spec, double T);

double _sk_from_T_equilibrium(long spec, double T);

double _T_from_w_e(spec_t w, double e);

double _T_from_w_h(spec_t w, double h);

double _T_from_rhok_P(spec_t rhok, double P);

double _T_from_w_rho_P(spec_t w, double rho, double P);

void find_nuk_eta_kappa(spec_t w, double rho, double T, 
                   spec_t nuk, double *eta, double *kappa);
		     
double _de_dP_at_constant_rho(spec_t rhok, double T);

double _de_dT_at_constant_rho(spec_t rhok, double T);

void find_de_drhok_at_constant_P(spec_t rhok, double T, spec_t dedrhok);

void find_de_drhok_at_constant_T(spec_t rhok, double T, spec_t dedrhok);

void find_species_name(long spec, char **name);

void find_species_variable_name(long spec, char **name);

void find_default_init_mass_fraction_string(long spec, char **wdefault);

void find_default_init_number_density_string(long spec, char **Ndefault);

void find_default_init_number_density_fraction_string(long spec, char **chidefault);

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

double _dmueN_dTe_from_Te(double Te);

double _mueN_from_Te(double Te);

double _dmukN_dTk_from_Tk_EkoverN(double Tk, double Ekstar, long k);


/* find the mobility of species k using the species temperature Tk and electric field in the species reference frame Ek (N is the total number density of the plasma, Tk the species temperature and Ek the electric field in the species frame)*/ 
double _muk_from_N_Tk_Ek(double N, double Tk, double Ek, long k);

double _zetav_from_Te(double Te);

double _zetae_from_Te(double Te);

double _dzetae_dTe_from_Te(double Te);

void find_Te_from_EoverN(double Estar, double *Te);

void find_dTe_dEoverN_from_EoverN(double Estar, double *dTedEstar);

void find_EoverN_from_Te(double Te, double *EoverN);


#endif /* _THERMO_H */
