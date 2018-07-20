#ifndef _CHEM_SHARE_H
#define _CHEM_SHARE_H

#include <model/_model.h>
#include <src/common.h>

void add_to_W_Arrhenius2r2p(int specR1, int specR2,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

void add_to_W_Arrhenius2r3p(int specR1, int specR2,
                            int specP1, int specP2, int specP3,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

void add_to_dW_Arrhenius2r2p(int specR1, int specR2,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX);

void add_to_dW_Arrhenius2r3p(int specR1, int specR2,
                             int specP1, int specP2, int specP3,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX);





void add_to_W_Arrhenius3r2p(int specR1, int specR2, int specR3,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

void add_to_W_Arrhenius2r1p(int specR1, int specR2,
                            int specP1,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

void add_to_dW_Arrhenius2r1p(int specR1, int specR2,
                             int specP1, 
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX);

void add_to_dW_Arrhenius3r2p(int specR1, int specR2, int specR3,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdX);




void test_dW_dx(spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam);

#endif /* _CHEM_SHARE_H */
