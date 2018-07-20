#ifndef _TS_H
#define _TS_H

#include <src/common.h>
#include <model/_model.h>

void update_U(np_t *np, gl_t *gl, zone_t zone);

void find_dU(np_t *np, gl_t *gl, zone_t zone);

#endif /* _TS_H */ 
