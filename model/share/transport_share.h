#ifndef _TRANSPORT_SHARE_H
#define _TRANSPORT_SHARE_H

#include <model/_model.h>
#include <src/common.h>

void find_nuk_eta_kappa(spec_t rhok, double T, double Te,
                   spec_t nuk, double *eta, double *kappa);

double _kappac_from_rhok_Tk_Ek(spec_t rhok, double T, double Te, double E, long k);

double _etac_from_rhok_Tk_Ek(spec_t rhok, double T, double Te, double E, long k);

#endif /* _TRANSPORT_SHARE_H */
