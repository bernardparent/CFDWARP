#ifndef _TS_SHARE_H
#define _TS_SHARE_H

#include <model/_model.h>
#include <src/common.h>


void find_bdry_jacobian(np_t *np, gl_t *gl, long bdrytype, long lA, long lB,
                 int TYPELEVEL, sqmat_t DA, sqmat_t DB);


void find_TDMA_jacobians_conservative(np_t *np, gl_t *gl, long theta, long l, sqmat_t B, sqmat_t C);

void add_TDMA_jacobians_non_conservative(np_t *np, gl_t *gl, long theta, long l, sqmat_t A, sqmat_t B, sqmat_t C, bool SOURCETERM);

void add_TDMA_jacobian_diagonally_dominant(np_t *np, gl_t *gl, long theta, long l, sqmat_t B);

void add_Kstar_dG_dUstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C);

void add_Dstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C);

void add_Ystar_dH_dUstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C);

void add_dSstar_dUstar_to_TDMA(np_t *np, gl_t *gl, long l,  double fact, sqmat_t B);

void add_Z_dUstardt_dUstar_to_TDMA(np_t *np, gl_t *gl, long l,  sqmat_t B);

void add_dFstar_dUstar_to_TDMA(np_t *np, gl_t *gl, long theta, long l, double fact, sqmat_t A, sqmat_t B, sqmat_t C);


#endif /* _TS_SHARE_H */
