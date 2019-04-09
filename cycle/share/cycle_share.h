#ifndef _CYCLE_SHARE_H
#define _CYCLE_SHARE_H

#include <model/_model.h>
#include <src/common.h>

#define TSEMF_DEFAULT 0
#define TSEMF_ADI    1
#define TSEMF_DDADI  2
#define TSEMF_IMAF   3
#define TSEMF_NEWTON 4
#define TSEMF_ADIIMAF   5
#define TSEMF_ADIi   6
#define TSEMF_ADIk   7
#define TSEMF_IMAFk   8
#define TSEMF_IMAFi   9
#define TSEMF_SOR    12
#define TSEMF_SOR2    13

#define SWEEPTYPE_IJK 0
#define SWEEPTYPE_I 1
#define SWEEPTYPE_J 2
#define SWEEPTYPE_K 3

#define IJK_UPDATE_YES TRUE
#define IJK_UPDATE_NO  FALSE

#define SEGMENTWORK_HEAVY 1
#define SEGMENTWORK_LIGHT 2

#define GRIDLEVEL_ONE 1
#define GRIDLEVEL_TWO 2

#define PRECON_CONSTANTTIMESTEP 10
#define PRECON_LOCALTIMESTEP 11
#define PRECON_LOCALEIGENVALUE 12
#define PRECON_LOCALEIGENVALUE2 13

#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_MUSCLVARS
#define nmc (2*nf)
#else
#define nmc nf
#endif

typedef double musclvarscycle_t[nmc];


typedef struct{
  zone_t *res;
  zone_t *ts;
  zone_t *bdry;
  long numzones_ts,numzones_total,numzones_res,numzones_bdry;
#ifdef DISTMPI
  zone_t domain;
#endif
} multizone_t;

void find_musclvarscycle(np_t np, gl_t *gl, musclvarscycle_t musclvars);

void sweep_with_1D_segments(np_t *np, gl_t *gl, zone_t zone,
                   void funct(np_t *, gl_t *, long, long, long), int sweeptype,
		   int TYPELEVEL, bool is_node_valid_local(np_t, int), int SEGMENTWORK, int GRIDLEVEL);

/* resume nodes from gl->is to gl->ie inclusively, and
   suspend all the others */
void resume_nodes_only_in_zone_and_update_bdry_nodes(np_t *np, gl_t *gl, zone_t zone);

void resume_nodes_specified_in_function(np_t *np, gl_t *gl,
                                bool(*FUNCT)(gl_t *, long, long, long));

void resume_nodes_only_in_zone(np_t *np, gl_t *gl, zone_t zone);

void update_bdry_nodes(np_t *np, gl_t *gl, zone_t zone);

void update_linked_nodes(np_t *np, gl_t *gl, int TYPELEVEL);

void add_int_to_codex(SOAP_codex_t *codex, char *varname, long varvalue);

void add_double_to_codex(SOAP_codex_t *codex, char *varname, double varvalue);

void process_code_runtime(np_t *np, gl_t *gl, char *code_runtime, SOAP_codex_t *codex);

/* finds the maximum _xi in zone, and puts it in gl->ximax.
   if IJK_UPDATE is set to IJK_UPDATE_YES, the position of ximax will be
   stored in gl->i_ximax, gl->j_ximax, gl->j_ximax */
void find_ximax(np_t *np, gl_t *gl, zone_t zone, int IJK_UPDATE);

void update_runtime_codex_vars_except_xi_from_gl(gl_t *gl, SOAP_codex_t *codex);

void update_runtime_codex_xi_from_gl(gl_t *gl, SOAP_codex_t *codex);

void check_residual(np_t *np, gl_t *gl, zone_t zone);

/* setup a multizone variable, with U updated in zone,
   the boundary nodes updated in zone+hbw_bdry_fluid and the residual updated
   in zone+hbw_bdry_fluid+hbw_res_fluid; lim defines the limits of the domain that _can_ be
   updated*/
void setup_multizone(np_t *np, gl_t *gl, zone_t zone, zone_t lim, double xiverge,
                    long zonelength, bool UPDATE_ALL_ZONES, multizone_t *multizone);

void update_bdry_nodes_with_multizone(np_t *np, gl_t *gl, multizone_t multizone);

void find_residual_with_multizone(np_t *np, gl_t *gl, multizone_t multizone);

void update_U_with_multizone(np_t *np, gl_t *gl, multizone_t multizone);

void solve_multizone(np_t *np, gl_t *gl, multizone_t multizone);

void solve_multizone_reverse(np_t *np, gl_t *gl, multizone_t multizone);

void free_multizone(multizone_t *multizone);

void exchange_U(np_t *np, gl_t *gl);

double _xi(np_t np, gl_t *gl, flux_t Res);

void find_constant_dtau(np_t *np, gl_t *gl, long l, double *dtau);

void find_dtau(np_t *np, gl_t *gl, long l, flux_t dtau);

void read_cycle_actions(char *actionname, char **argum, SOAP_codex_t *codex);

#ifdef UNSTEADY
void increase_time_level(np_t *np, gl_t *gl);
#endif

#ifdef EMFIELD
void exchange_U_emfield(np_t *np, gl_t *gl);

void setup_multizone_emfield(np_t *np, gl_t *gl, zone_t zone, zone_t lim, double xiverge,
                             long zonelength, multizone_t *multizone);

void solve_multizone_emfield(np_t *np, gl_t *gl, multizone_t multizone);

void update_bdry_nodes_emfield(np_t *np, gl_t *gl, zone_t zone);

void find_residual_emfield(np_t *np, gl_t *gl, zone_t zone);

void find_ximax_emfield(np_t *np, gl_t *gl, zone_t zone);

void reupdate_nearby_bdry_nodes_emfield(np_t *np, gl_t *gl, long l);

void update_prim_emfield_mem_in_zone(np_t *np, gl_t *gl, zone_t zone);

void update_Te_local_in_zone(np_t *np, gl_t *gl, zone_t zone);

void read_UpdateEMField_arguments(char **argum, SOAP_codex_t *codex, gl_t *gl);

void solve_TDMA_emfield(np_t *np, gl_t *gl, long theta, long ls, long le, 
                        int TYPELEVEL, EXM_tdmaline_t *tdma, long numlines);
#endif

#endif /* _CYCLE_SHARE_H */
