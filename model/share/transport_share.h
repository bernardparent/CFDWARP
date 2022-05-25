#ifndef _TRANSPORT_SHARE_H
#define _TRANSPORT_SHARE_H

#include <model/_model.h>
#include <src/common.h>


double _kappac_from_rhok_Tk_Ek(spec_t rhok, double T, double Te, double E, long k);

double _etac_from_rhok_Tk_Ek(spec_t rhok, double T, double Te, double E, long k);

void adjust_nuk_using_mobilities(spec_t rhok, double T, double Te, spec_t nuk);

void adjust_muk_for_Ek_effect(long k, double Ek, double N, double *muk);

#endif /* _TRANSPORT_SHARE_H */
