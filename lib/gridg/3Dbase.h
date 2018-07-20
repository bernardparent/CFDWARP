#ifndef _3DBASE_H
#define _3DBASE_H

#include "gridg.h"
#include <assert.h>

void Size3D(GRIDG_gl3d_t *gl, GRIDG_xyzgrid_t **xyzgrid,
            long is, long js, long ks, long ie, long je, long ke);

void Point3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid,
             long i, long j, long k, double x, double y, double z);

void Corners3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid,
                   long i1, long j1, long k1, long i2, long j2, long k2,
                   double x1, double y1, double z1, double x2, double y2, double z2);

int JoinCorners3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
            long i1, long j1, long k1, long i2, long j2, long k2,
            char modeL_i, char modeH_i, double NLoverN_i, double valL_i, double valH_i,
            char modeL_j, char modeH_j, double NLoverN_j, double valL_j, double valH_j,
            char modeL_k, char modeH_k, double NLoverN_k, double valL_k, double valH_k);

int Join3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long i1, long j1, long k1, long i2, long j2,
           long k2, char index, char modeL, char modeH,
           double NLoverN, double valL, double valH);

void JoinAll3D(GRIDG_gl3d_t gl, GRIDG_xyzgrid_t *xyz, SOAP_codex_t *codex,
               long i1, long j1, long k1, long i2, long j2, long k2);

void Plane(GRIDG_gl3d_t gl, GRIDG_xyzgrid_t *xyz, SOAP_codex_t *codex,
           long i1, long j1, long k1, long i2, long j2, long k2);

void Translate3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           double dx, double dy, double dz);

void Rotate3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           double pivot_x, double pivot_y, double pivot_z,
           double angle, char axis);

void Copy3D(GRIDG_gl3d_t gl3d, GRIDG_xyzgrid_t *xyzgrid, SOAP_codex_t *codex,
           long is, long js, long ks, long ie, long je, long ke,
           long id, long jd, long kd);


#endif
