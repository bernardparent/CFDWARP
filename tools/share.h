#ifndef _SHARE_H
#define _SHARE_H

#include <exm.h>
#include <memdebug.h>
#include <string.h>

void FindConvMachNumber(double gamma1, double R1, double P1, double T1, double U1,  
                        double gamma2, double R2, double P2, double T2, double U2,  
                        double *MC1, double *MC2, double *UC, bool *noroot);

void FindShearLayerGrowth(double gamma1, double R1, double P1, double T1, double U1,  
                    double gamma2, double R2, double P2, double T2, double U2,  
                    double UC, double MC1, double MC2, double *dshear, 
                    double *betac, double *betad);

/* from wedge angle del (in radians) find the props after the shock 
   from a phi limit of phi1<phi<phi2*/
void FindWedgeShockProps(double gamma, double del, double M1,  
                     double *M2, double *T2overT1, double *P2overP1,
                     double *phiout, bool *valid);

/* adiabatic wall only!! and M<6*/
double VanDriestIICorrel(double M);

double SkinFrictionCoeff(double M, double Re, double x, bool TURB, bool COMP);

double BdryLayerThickness(double M, double Re, double x, bool TURB, bool COMP);

/* from M1 (Mach number) and gam (gamma), find the shock 
   properties */
void FindNormalShockProps(double M1, double gam, double *M2, double *P2overP1, 
                      double *T2overT1, double *U2overU1, 
                      double *Pstag2overPstag1);

void FindTandHfromPstdatm(double P, double *H, double *T);

void FindExtCompInletProps(double M1, double L, double Rgas,
               double T3, double Pdyn, double gam,
               double *T1, double *T2, double *P1,
               double *P2, double *P3, double *Ms,
               double *M2, double *M3, double *W1,
               double *mdot, double *altitude);



#endif /* _SHARE_H */ 
