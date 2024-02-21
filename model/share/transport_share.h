#ifndef _TRANSPORT_SHARE_H
#define _TRANSPORT_SHARE_H

#include <model/_model.h>
#include <src/common.h>


void find_nuk_from_Dij(spec_t rhok, spec2_t Dij, spec_t nuk);

void adjust_nue_given_nui(spec_t rhok, double T, double Te, spec_t nuk);

double _kappac_from_rhok_Tk_muk(spec_t rhok, double T, double Te, double muk, long k);

double _etac_from_rhok_Tk_muk(spec_t rhok, double T, double Te, double muk, long k);

void adjust_nuk_using_mobilities_given_muk(spec_t rhok, double T, double Te, chargedspec_t muk, spec_t nuk);

void find_muk_from_nuk(spec_t nuk, spec_t rhok, double T, double Te, chargedspec_t muk);


void adjust_muk_for_Ek_effect(long k, double Ek, double N, double *muk);

void find_dmuk_from_rhok_Tk_Ek_ParentMacheret(spec_t rhok, double Tk, double Ek, long k, double *dmukdTk, spec_t dmukdrhok);

double _muk_from_rhok_T_Te_ParentMacheret(spec_t rhok, double T, double Te, long k);

double _lnLambda(spec_t rhok, double Te);

double _mueNk_from_Te_ParentMacheret(long spec, double Te);

#endif /* _TRANSPORT_SHARE_H */
