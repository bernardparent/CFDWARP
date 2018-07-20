#ifndef _RES_H
#define _RES_H

#include <src/common.h>

/* find the residual of all nodes comprised in zone */

void find_residual(np_t *np, gl_t *gl, zone_t zone);

void find_residual_per_dimension(np_t *np, gl_t *gl, zone_t zone, int sweeptype);


#endif /* _RES_H */ 
