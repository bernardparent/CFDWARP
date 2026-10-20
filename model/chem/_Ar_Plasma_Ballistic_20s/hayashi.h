#ifndef HAYASHI_H

#define HAYASHI_H
#include <src/common.h>


double _sigma_el_Hayashi(double eps);

double _sigma_inel_Hayashi(long k, double eps);

double _Deltaeps_inel_Hayashi(long k);

double _maxloss_Hayashi(double *epsmax);

#endif
