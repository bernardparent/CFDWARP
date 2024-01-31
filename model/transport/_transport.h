#ifndef TRANSPORT_H

#define TRANSPORT_H
#include <src/common.h>
#include <model/share/transport_share.h>  


void write_model_transport_template(FILE **controlfile);

void read_model_transport_actions(char *actionname, char **argum, SOAP_codex_t *codex);

void find_dmuk_from_rhok_Tk_Ek(gl_t *gl, spec_t rhok, double Tk, double Ek, long k, double *dmukdTk, spec_t dmukdrhok);

void find_nuk_eta_kappa(gl_t *gl, spec_t rhok, double T, double Te, spec_t nuk, double *eta, double *kappa);

void find_nuk_eta_kappak_muk(gl_t *gl, spec_t rhok, double T, double Te,
                             spec_t nuk, double *eta, double *kappan, chargedspec_t kappac, chargedspec_t muk);

#endif
