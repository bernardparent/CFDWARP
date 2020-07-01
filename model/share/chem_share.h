#ifndef _CHEM_SHARE_H
#define _CHEM_SHARE_H

#include <model/_model.h>
#include <src/common.h>


/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   W in kg m^(-3) s^(-1)
*/

void add_to_W_2r1p ( int specR1, int specR2,
                     int specP1, 
                     double kf, spec_t N, spec_t W);

void add_to_W_2r2p ( int specR1, int specR2,
                     int specP1, int specP2,
                     double kf, spec_t N, spec_t W);

void add_to_W_2r3p ( int specR1, int specR2,
                     int specP1, int specP2, int specP3,
                     double kf, spec_t N, spec_t W);


/* kf in cm^6 s^(-1) 
   N in cm^(-3)
   W in kg m^(-3) s^(-1)
*/

void add_to_W_3r2p ( int specR1, int specR2, int specR3,
                     int specP1, int specP2,
                     double kf, spec_t N, spec_t W);

void add_to_W_3r3p ( int specR1, int specR2, int specR3,
                     int specP1, int specP2, int specP3,
                     double kf, spec_t N, spec_t W);


/* kf in cm^3 s^(-1) 
   N in cm^(-3)
   dkfdT in cm^3 s^(-1) K^(-1)
   dkfdTv in cm^3 s^(-1) K^(-1)
   dkfdTe in cm^3 s^(-1) K^(-1)
   dWdrhok in s^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdTv in kg m^(-3) s^(-1) K^(-1)
   dWdTe in kg m^(-3) s^(-1) K^(-1)
*/
void add_to_dW_2r1p ( int specR1, int specR2, int specP1, double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe);

void add_to_dW_2r2p ( int specR1, int specR2, int specP1, int specP2, double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe);


void add_to_dW_2r3p ( int specR1, int specR2, int specP1, int specP2, int specP3, double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe);

/* kf in cm^6 s^(-1) 
   N in cm^(-3)
   dkfdT in cm^6 s^(-1) K^(-1)
   dkfdTv in cm^6 s^(-1) K^(-1)
   dkfdTe in cm^6 s^(-1) K^(-1)
   dWdrhok in s^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdTv in kg m^(-3) s^(-1) K^(-1)
   dWdTe in kg m^(-3) s^(-1) K^(-1)
*/

void add_to_dW_3r2p ( int specR1, int specR2, int specR3, int specP1, int specP2,  
                      double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe);
                      
void add_to_dW_3r3p ( int specR1, int specR2, int specR3, int specP1, int specP2, int specP3, 
                      double kf, spec_t N, 
                      double dkfdT, double dkfdTv, double dkfdTe, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTv, spec_t dWdTe);

/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   W in kg m^(-3) s^(-1)
*/

void add_to_W_fwbw_2r1p(int specR1, int specR2,
                            int specP1,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

void add_to_W_fwbw_2r2p(int specR1, int specR2,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

void add_to_W_fwbw_2r3p(int specR1, int specR2,
                            int specP1, int specP2, int specP3,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

void add_to_W_fwbw_3r2p(int specR1, int specR2, int specR3,
                            int specP1, int specP2,
                            double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t W);

/* A in cm^3 (mole s)^(-1) K^(-n) 
   E in cal mole^(-1) 
   T in Kelvin
   X in mole/cm3
   Gs in J/mole
   dGsdT in J mole^(-1) K^(-1)
   dWdT in kg m^(-3) s^(-1) K^(-1)
   dWdrhok in s^(-1)
*/

void add_to_dW_fwbw_2r1p(int specR1, int specR2,
                             int specP1, 
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdrhok);

void add_to_dW_fwbw_2r2p(int specR1, int specR2,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdrhok);

void add_to_dW_fwbw_2r3p(int specR1, int specR2,
                             int specP1, int specP2, int specP3,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdrhok);

void add_to_dW_fwbw_3r2p(int specR1, int specR2, int specR3,
                             int specP1, int specP2,
                             double A, double n, double E, double T, spec_t X, spec_t Gs, spec_t dGsdT,
                             spec_t dWdT, spec2_t dWdrhok);


/* find the forward reaction rate cofficient kf = A * (calA)^(1-numreactant) * T^n * exp(-E/(R*T))
 * 
 * numreactant : the number of reactants
 * A           : pre-exponential factor in cm^3 * (mole s)^(-1) K^(-n)
 * n           : temperature exponent
 * E           : activation energy of the reaction in cal/mole
 * kf          : reaction rate coefficient in cm^3/s (2 reactants) or cm^6/s (3 reactants)
 * T           : temperature in K
 */
double _kf_Arrhenius(long numreactant, double A, double n, double E, double T);


/* find derivative of kf with respect to T with kf=A * (calA)^(1-numreactant) * T^n * exp(-E/(R*T))
 * 
 * numreactant : the number of reactants
 * A           : pre-exponential factor in cm^3 * (mole s)^(-1) K^(-n)
 * n           : temperature exponent
 * E           : activation energy of the reaction in cal/mole
 * dkfdT       : cm^3/sK (2 reactants) or cm^6/sK (3 reactants)
 * T           : temperature in K
 */
double _dkfdT_Arrhenius(long numreactant, double A, double n, double E, double T);






void test_dW_dx(spec_t rhokref, spec_t rhok, spec_t mu, double T, double Te, double Tv, double Estar, double Qbeam);

#endif /* _CHEM_SHARE_H */
