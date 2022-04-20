#ifndef TRANSPORT_H

#define TRANSPORT_H
#include <src/common.h>
#include <model/share/transport_share.h>  

double _kappa_from_rhok_T_Te(spec_t rhok, double T, double Te);

double _kappan_from_rhok_T_Te(spec_t rhok, double T, double Te);

double _eta_from_rhok_T_Te(spec_t rhok, double T, double Te);

double _etan_from_rhok_T_Te(spec_t rhok, double T, double Te);

void find_nuk_from_rhok_T_Te(spec_t rhok, double T, double Te, spec_t nuk);

double _muk_from_rhok_T_Te_Ek(spec_t rhok, double T, double Te, double Ek, long k);

void find_dmuk_from_rhok_Tk_Ek(spec_t rhok, double Tk, double Ek, long k, double *dmukdTk, spec_t dmukdrhok);


#endif
