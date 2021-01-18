#ifndef _2DBASE_H
#define _2DBASE_H

#include "gridg.h"

void Size2D(GRIDG_gl2d_t *gl, GRIDG_xygrid_t **xygrid,
            long is, long js, long ie, long je);

void Point2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid,
             long i, long j, double x, double y);

void Corners2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid,
                   long i1, long j1, long i2, long j2,
                   double x1, double y1, double x2, double y2);

int JoinCorners2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
            long i1, long j1, long i2, long j2,
            char modeL_i, char modeH_i, double NLoverN_i, double valL_i, double valH_i,
            char modeL_j, char modeH_j, double NLoverN_j, double valL_j, double valH_j);

int Join2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long i1, long j1, long i2, long j2,
           char index, char modeL, char modeH,
           double NLoverN, double valL, double valH);

void Translate2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           double dx, double dy);

void Rotate2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           double pivot_x, double pivot_y,
           double angle);

void Copy2D(GRIDG_gl2d_t gl2d, GRIDG_xygrid_t *xygrid, SOAP_codex_t *codex,
           long is, long js, long ie, long je,
           long id, long jd);


#endif
