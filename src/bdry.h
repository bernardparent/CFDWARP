#ifndef _BDRY_H
#define _BDRY_H


#include <src/common.h>
#include <model/_model.h>


double _bdry_param(np_t *np, gl_t *gl, long l, short param, int TYPELEVEL);

char _bdry_ID(long nodetype);

/* given the argument to the Bdry() module in the control file,
   and the codex (code extra variables) this subroutine processes the argument
   and initializes the boundary nodes
*/
void read_bdry(char *argum, SOAP_codex_t *codex);

/* given the pointer to the controlfile file, this appends to the
   file a template Bdry(); to specify the boundary nodes
*/
void write_bdry_template(FILE **controlfile);

/* in the zone defined by zone, the nodes' types are updated
   according to the type specified. If TYPELEVEL is set to TYPELEVEL_FLUID_WORK,
   the node types (the node types are the boundary conditions) used
   by the iterative process will be modified (temporary node types).
   If TYPELEVEL is set to TYPELEVEL_FLUID, then the permanent node type
   is altered
*/
void update_node_type(np_t *np, gl_t *gl, int TYPELEVEL, long type, zone_t zone);


/* same as above, but only updates the nodes that were previously indicated
   as inner nodes
*/
void update_inner_to_bdry_node_type(np_t *np, gl_t *gl, int TYPELEVEL, long type,
                               zone_t zone);


/* same as above, but only update the nodes previously identified as
   inner and bdry nodes (not the unused nodes)
*/
void update_inner_and_bdry_to_node_type(np_t *np, gl_t *gl, int TYPELEVEL, long type,
                                  zone_t zone);


/* For the bdry node number l_A, finds the boundary direction, with
   the variable theta and thetasgn (thetasgn can be either +1 or -1)
   returns TRUE if found, FALSE if not found. The boundary direction
   here means in which direction are the inner nodes on which the
   boundary node depends on. For example, if the boundary node depends
   on two inner nodes in the positive j direction, then theta=1 and
   thetasgn=1
*/
bool find_bdry_direc(np_t *np, gl_t *gl, long l_A, int TYPELEVEL, long *theta,
                     long *thetasgn);

void find_link_direc(np_t *np, gl_t *gl, long llink, int TYPELEVEL,
                     long *theta, long *thetasgn);

bool is_multiple_bdry_direc(np_t *np, gl_t *gl, long l_A, int TYPELEVEL);

/* for all nodes in zone, adjust the nodes' type: i.e. make sure that
   the boundary nodes are valid and modify them accordingly.
   TYPELEVEL can be set to TYPELEVEL_FLUID_WORK or TYPELEVEL_FLUID to modify
   the node types on the temporary level or the permanent level (BASE)
*/
void adjust_node_type(np_t *np, gl_t *gl, zone_t zone, int TYPELEVEL);

void display_zone(FILE *outfile, zone_t zone);

void display_node_type_window(FILE *outfile, np_t *np, gl_t *gl, int TYPELEVEL,
                        long im, long jm, long km, long bw);

void display_node_type_window_local_process(FILE *outfile, np_t *np, gl_t *gl, int TYPELEVEL,
                        long im, long jm, long km, long bw);


/* show, to stdout, the node types on level TYPELEVEL (either
   TYPELEVEL_FLUID or TYPELEVEL_FLUID_WORK) in zone
*/
void display_node_type(np_t *np, gl_t *gl, zone_t zone, int TYPELEVEL);


/* For all nodes part of zone, copy the node type from base to
   work (i.e. copy the permanent node type to the temporary node type
*/
void copy_base_to_work_node_type(np_t *np, gl_t *gl, zone_t zone);

#endif /* _BDRY_H */
