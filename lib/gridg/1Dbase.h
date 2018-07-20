#ifndef _1DBASE_H
#define _1DBASE_H

#include "gridg.h"
#include <assert.h>

void Size1D(GRIDG_gl1d_t *gl, GRIDG_xgrid_t **xgrid,
            long is, long ie);

void Point1D(GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid,
             long i, double x, double Acs);


int JoinCorners1D(GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid, SOAP_codex_t *codex,
            long i1, long i2,
            char modeL_i, char modeH_i, double NLoverN_i, double valL_i, double valH_i);

int Join1D(GRIDG_gl1d_t gl1d, GRIDG_xgrid_t *xgrid, SOAP_codex_t *codex,
           long i1, long i2,
           char modeL, char modeH,
           double NLoverN, double valL, double valH);

#endif
