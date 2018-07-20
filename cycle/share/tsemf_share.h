#ifndef _TSEMF_SHARE_H
#define _TSEMF_SHARE_H

#define TSEMF_SOR_NUMSUBITEROPT 200

#include <model/_model.h>
#include <src/common.h>


void find_tsemf_SOR_numsubiter_numcycle(long numsubitertot,long numsubiteropt, long *numsubiter, long *numcycle);

void update_U_from_dUstar_emfield(np_t *np, gl_t *gl, long theta, long ls, long le);

void update_U_from_dUstar_emfield_without_relaxation(np_t *np, gl_t *gl, long theta, long ls, long le);

void update_U_and_Res_from_dUstar_emfield(np_t *np, gl_t *gl, long theta, long ls, long le);

void update_U_emfield_ADI(np_t *np, gl_t *gl, zone_t zone);

void update_U_emfield_SOR(np_t *np, gl_t *gl, zone_t zone);

void init_dUstar_emfield_ADI(np_t *np, gl_t *gl, long theta, long ls, long le);

void update_dUstar_emfield_ADI(np_t *np, gl_t *gl, long theta, long ls, long le);

void update_U_emfield_ADIi(np_t *np, gl_t *gl, zone_t zone);

void update_U_emfield_ADIk(np_t *np, gl_t *gl, zone_t zone);

void update_U_emfield_Newton(np_t *np, gl_t *gl, zone_t zone);

void update_dUstar_emfield_SOR_istation(np_t *np, gl_t *gl, long flux, long i, zone_t zone, int SOR_SWEEP, long iter);


#endif /* _TSEMF_SHARE_H */
