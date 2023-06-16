#ifndef _COMMON_H
#define _COMMON_H




#if (defined(POSIXTHREADS) || defined(ZONETHREADS))
  #define PTHREADS
#endif

#if (defined(PTHREADS) || defined(OPENMPTHREADS))
  #define THREADS
#endif

#define THREADTYPE_ALL 0
#define THREADTYPE_POSIX 1
#define THREADTYPE_POSIX_LOOP 2
#define THREADTYPE_POSIX_ZONE 3
#define THREADTYPE_OPENMP 4
#define THREADTYPE_LOOP 5
#define THREADTYPE_ZONE 6

#define CYCLELEVEL_RES 1
#define CYCLELEVEL_TS 2

#ifdef DISTMPI
  #define DISTDOMAIN_METHOD_AUTO 1
  #define DISTDOMAIN_METHOD_USERSPEC   2
#endif

#define nodeoverlap (max(hbw_res_emfield,max(hbw_res_fluid,hbw_bdry_fluid)+hbw_mem_fluid))

#ifdef PTHREADS
  #define thread_lock_t pthread_mutex_t
#endif
#ifdef OPENMPTHREADS
  #define thread_lock_t omp_lock_t
#endif
#if (!defined(PTHREADS) && !defined(OPENMPTHREADS))
  #define thread_lock_t int
#endif


#ifdef DISTMPI
  #include "mpi.h"
  #define MPI_NO_DATA_EXCHANGED 0
  #define MPI_DATA_SENT 1
  #define MPI_DATA_RECEIVED 2
#endif

#ifdef _1D
  #error CFDWARP can not be compiled in 1 dimension
#endif

#ifdef _2D
  #define nd 2
  #define _2DL
#endif

#ifdef _3D
  #define nd 3
  #define _2DL
  #define _3DL
#endif


  #define for_1DL(i,is,ie)   for (i=is; i<=ie; i++)
  #define if1DL(a) a

#ifdef _2DL
  #define for_2DL(j,js,je)   for (j=js; j<=je; j++)
  #define if2DL(a) a
  #define ifn2DL(a)
#else
  #define for_2DL(j,js,je)   for (j=js; j<=js; j++)
  #define if2DL(a)
  #define ifn2DL(a) a
#endif

#ifdef _3DL
  #define for_3DL(k,ks,ke)   for (k=ks; k<=ke; k++)
  #define if3DL(a) a
  #define ifn3DL(a)
#else
  #define for_3DL(k,ks,ke)   for (k=ks; k<=ks; k++)
  #define if3DL(a)
  #define ifn3DL(a) a
#endif

  #define if1D(a)
  #define ifn1D(a) a

#ifdef _2D
  #define if2D(a) a
  #define ifn2D(a)
#else
  #define if2D(a)
  #define ifn2D(a) a
#endif

#ifdef _3D
  #define if3D(a) a
  #define ifn3D(a)
#else
  #define if3D(a)
  #define ifn3D(a) a
#endif

#ifdef _1D
  #define for_ijk(zone,is,js,ks,ie,je,ke)  for (i=zone.is; i<=zone.ie; i++) for (j=zone.js; j<=zone.js; j++) for (k=zone.ks; k<=zone.ks; k++)
  #define for_jik(zone,is,js,ks,ie,je,ke)  for (j=zone.js; j<=zone.js; j++) for (i=zone.is; i<=zone.ie; i++) for (k=zone.ks; k<=zone.ks; k++)
  #define for_kij(zone,is,js,ks,ie,je,ke)  for (k=zone.ks; k<=zone.ks; k++) for (i=zone.is; i<=zone.ie; i++) for (j=zone.js; j<=zone.js; j++) 
  #define for_kji(zone,is,js,ks,ie,je,ke)  for (k=zone.ks; k<=zone.ks; k++) for (j=zone.js; j<=zone.js; j++) for (i=zone.is; i<=zone.ie; i++)  
#endif

#ifdef _2D
  #define for_ijk(zone,is,js,ks,ie,je,ke)  for (i=zone.is; i<=zone.ie; i++) for (j=zone.js; j<=zone.je; j++) for (k=zone.ks; k<=zone.ks; k++)
  #define for_jik(zone,is,js,ks,ie,je,ke)  for (j=zone.js; j<=zone.je; j++) for (i=zone.is; i<=zone.ie; i++) for (k=zone.ks; k<=zone.ks; k++)
  #define for_kij(zone,is,js,ks,ie,je,ke)  for (k=zone.ks; k<=zone.ks; k++) for (i=zone.is; i<=zone.ie; i++) for (j=zone.js; j<=zone.je; j++) 
  #define for_kji(zone,is,js,ks,ie,je,ke)  for (k=zone.ks; k<=zone.ks; k++) for (j=zone.js; j<=zone.je; j++) for (i=zone.is; i<=zone.ie; i++)  
#endif

#ifdef _3D
  #define for_ijk(zone,is,js,ks,ie,je,ke)  for (i=zone.is; i<=zone.ie; i++) for (j=zone.js; j<=zone.je; j++) for (k=zone.ks; k<=zone.ke; k++)
  #define for_jik(zone,is,js,ks,ie,je,ke)  for (j=zone.js; j<=zone.je; j++) for (i=zone.is; i<=zone.ie; i++) for (k=zone.ks; k<=zone.ke; k++)
  #define for_kij(zone,is,js,ks,ie,je,ke)  for (k=zone.ks; k<=zone.ke; k++) for (i=zone.is; i<=zone.ie; i++) for (j=zone.js; j<=zone.je; j++) 
  #define for_kji(zone,is,js,ks,ie,je,ke)  for (k=zone.ks; k<=zone.ke; k++) for (j=zone.js; j<=zone.je; j++) for (i=zone.is; i<=zone.ie; i++)  
#endif

#ifdef NDEBUG
  #define ifndebug(a) a
  #define ifdebug(a)
  #define assert_np(np,expr)           ((void) 0)
#else
  #define ifdebug(a) a
  #define ifndebug(a)
  #define assert_np(np,expr)                                                   \
      ((expr) ? 0 : wfprintf(stderr,"ASSERTION %s FAILED at node i=%ld  j=%ld"\
                                   "  k=%ld, in function %s of file=%s at line=%d.\n",      \
                   __STRING(expr),(np).i,(np).j,(np).k,__FUNCTION__,__FILE__, __LINE__));
#endif

#define CHARACTER_BLOCK '*'
#define TYPELEVEL_FLUID 1
#define TYPELEVEL_FLUID_WORK 0
#define NODETYPE_UNUSED -3
#define NODETYPE_INNER -1
#define NODETYPE_BDRY 0 /* and more */
#define LINK_NONE -3
#define LINKDIM_NONE -3
#define LINKDIMSGN_NONE -3
#define LINKDIMSGN_BASE 10
#define EOS 0



#define MAX_LINE_WIDTH 97

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <exm.h>
#include <gridg.h>
#include <memdebug.h>
#include <string.h>
#include <pthread.h>
#ifdef OPENMPTHREADS
  #include <omp.h>
#endif
#include <limits.h>

typedef struct{
 long is,js,ks,ie,je,ke;
} zone_t;


typedef struct {
  long garbage; /* necessary to avoid compiler hang on AIX */
} voidarg_t;

typedef double dim_t[nd];
typedef double vec2_t[3][3];
typedef long longdim_t[nd];
typedef double dim2_t[nd][nd];

#include <model/.active/model.hh>
#include <cycle/.active/cycle.hh>
#ifdef EMFIELD
  #define TYPELEVEL_EMFIELD 2
#endif

typedef double sqmat_t[nf][nf];
typedef double initvar_t[numinitvar];
#ifdef EMFIELD
  typedef double initvar_emfield_t[numinitvar_emfield];
#endif

typedef struct {
  zone_t domain,window,domain_all,domain_lim,domain_lim_all;
  long iter,nn;
#ifdef _CYCLE_PREDICTOR_CORRECTOR
  long numsubiter_pc,subiter_pc;
#endif
  char *output_filename, *post_filename;
  bool CONVERGED,OUTPUTBINARYMPI,OUTPUTASCII,OUTPUTINTERPOLATION;
  int PRECONDITIONER;
  double effiter_U,effiter_R,ximax,xiave,CFL,sigma1,sigma2,dtau;
#ifdef EMFIELD
  double ximax_emfield,effiter_U_emfield,effiter_R_emfield;
  long i_ximax_emfield,j_ximax_emfield,k_ximax_emfield;
  double Lc,relaxEMF;
  long tsemfmethod,numsubiter_tsemf;
#endif
  long flux_ximax,i_ximax,j_ximax,k_ximax;
#if defined(UNSTEADY)
  double time,dt;
#endif
  thread_lock_t lock;
#ifdef _TS_BLOCK_IMAF
  long subiter_ts;
#endif
#ifdef DISTMPI
  bool DISTDOMAIN;
  int DISTDOMAIN_METHOD;
#endif
#if (defined(_TSEMF_IMAF) || defined(_TSEMF_IMAF_ADI))
  long subiter_tsemf;
#endif
  gl_cycle_t cycle;
  gl_model_t model;
  bool MODEL_FLUID_READ,MODEL_CHEM_READ,CYCLE_FLUID_READ,MODEL_EMFIELD_READ,CYCLE_EMFIELD_READ,MODEL_BEAM_READ,DISC_FLUID_READ,DISC_EMFIELD_READ,DISC_RESCONV_READ,DISC_RESTIME_READ;
  bool INIT_FLUID_READ,INIT_EMFIELD_READ,BDRY_FLUID_READ,BDRY_EMFIELD_READ,METRICS_INITIALIZED;
  bool CONTROL_READ,RESETRUNTIMEVARS,RESIDUAL_ALTERED;
  long initspecies[ns];
  long nsinit;
#ifdef EMFIELD
  bool RESIDUAL_ALTERED_EMFIELD;
#endif
#ifdef DISTMPI
  long numdomain_i,numdomain_j,numdomain_k;
  zone_t *domain_from_rank,*domain_lim_from_rank;
#endif
  bool REPORTCLIPPING;
#if defined(EMFIELD) && defined(DISTMPI)
  bool EM_MPIBDRY_EXPLICIT;
#endif
#if defined(POSIXTHREADS) || defined(OPENMPTHREADS)
  bool NOSHORTTHREADS;
#endif
#ifdef _TSEMF_STORE_COEFFICIENTS
  EXM_gl3D_t tsemfcoeffzone;
#endif
} gl_t;


/* the np working variables */
typedef struct {
   double xi,dtau;
   long flux_xi;
   flux_t Res;

   /* zone threading lock */
#ifdef ZONETHREADS
   thread_lock_t lock;
#endif

#ifdef _FLUID_NAVIERSTOKES
   spec_t numem;
   double Pmem,amem,athermomem;
   double etamem,kappamem,rhomem;
   double Tmem;
   dim_t Vmem;   
#endif


#ifdef _FLUID_FAVREREYNOLDS
   spec_t numem;
   double Pmem,amem,athermomem;
   double etamem,kappamem,etat,rhomem;
   double Tmem;
   dim_t Vmem;   
#endif

#ifdef _FLUID_PLASMA
   double Temem;
   chargedspec_t mumem;
   chargedspec_t Nkmem;
   chargedspec_t Pkmem;
   chargedspec_t kappacmem;
   double Nmem,Nnmem,kappanmem;
#endif

#ifdef _FLUID_NAVIERSTOKESPERFECT
   double Pmem,amem;
   double rhomem;
   double Tmem;
   dim_t Vmem;
#endif



   /* --Approximate Factorization specific */
   flux_t dUstar;
#ifdef _TS_BLOCK_IMAF
   flux_t dUinc;
#endif
#ifdef _TS_STORAGE_JACOBIANS
   sqmat_t TDMA_A[nd],TDMA_B[nd],TDMA_C[nd];
#endif

#ifdef _RESCONV_INCLUDES_DIFFUSION
   flux_t Fp1h_diffusion[nd];
#endif
#ifdef _RESCONV_STORAGE_FSTAR
   double *Fp1h;
#endif


} npwk_t;


/* the np base variables from which the working and metrics variables
   can be reconstructed*/
typedef struct {
   dim_t x;
   flux_t U;
   double Omega,Omega_int[nd];
   dim2_t X,X_int[nd];

#ifdef _FLUID_FBODY_QADD
   dim_t Fbody;
   double Qadd;
#endif
#ifdef EMFIELD
   fluxemfield_t Uemfield;
   fluxemfield_t Resemfield,dUstaremfield;
   double _xi_emfield;
#endif

#ifdef _TSEMF_STORE_COEFFICIENTS
   fluxemfield_t coeffm1[nd],coeffp0[nd],coeffp1[nd];
   fluxemfield_t coeffp0sum;   
   fluxemfield_t dtauemfield;
   fluxemfield_t tsemfcoeff[27];
#endif

#if (defined(_TSEMF_IMAF) || defined(_TSEMF_IMAF_ADI))
   fluxemfield_t dUincemfield;
#endif

#ifdef _TSEMF_SPLIT
   fluxemfield_t dUstaremfield_tmp;
#endif

#ifdef _EMFIELD_EPOTENTIAL_BFIXED
   EXM_vec3D_t B,E,J;
   dim_t EMFoverL;
   double EMFmaxpower;
   double sigmasolid,epsilonr;
   double sigmamem,Eemagmem,Eemagsmoothedmem;
   dim_t Esmoothedmem;
#endif
#ifdef _EMFIELD_EPOTENTIAL
   EXM_vec3D_t E;
   dim_t EMFoverL;
   double EMFmaxpower;
   double sigmasolid,epsilonr;
   double sigmamem,Eemagsmoothedmem;
   dim_t Esmoothedmem;
#endif
#ifdef _BEAM_EBEAM_FIXED
   double Qbeam;
#endif
#ifdef _BEAM_EBEAM_ALGEBRAIC
   double QbeamoverN;
#endif
#ifdef UNSTEADY
   flux_t Um1;
#ifdef EMFIELD
   flux_t Uemfieldm1;
#endif
#if _RESTIME_BW > 2
   flux_t Um2;
#endif
#if _RESTIME_BW > 3
   flux_t Um3;
#endif
#endif//UNSTEADY
#ifdef _RESTIME_RK4
   flux_t dUstar1,dUstar2,dUstar3,dUstar4;
#endif
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL 
   flux_t trapezoidalm1;
#endif
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL 
   flux_t trapezoidalm1_next;
#endif
#ifdef _RESCONV_DELTA_LAMBDA_STORAGE
   flux_t Delta_Lambda[nd];
#endif
#ifdef _RESSOURCE_LAMBDAMINUS_STORAGE
   flux_t Lambda_S_minus;
#endif
#ifdef _TSEMF_SOR2
  double *tsemfcoeff;
  long *tsemfnode;
  long tsemfnodenum;
  double tsemf_rhs;
#endif
#ifndef NDEBUG
  bool TSEMF_UPDATED;
#endif
#ifdef _FLUID_FBODY_QADD
   bool FBODY,QADD;   
#endif

} npbs_t;


typedef struct {
   npbs_t *bs;
   npwk_t *wk;
   char status;
   bool INIT_FLUID;
   bool FLUIDPRIMMEM;
   long type,type_wk;
   short numbdryparam;
   double *bdryparam;
   //long link;
   long *linkarray;
   short numlink;
#ifdef EMFIELD
   short numbdryparam_emf;
   double *bdryparam_emf;
   long type_emf;
   //long link_emf;
   bool INIT_EMFIELD;
   long *linkarray_emf;
   short numlink_emf;
#endif
#ifdef DISTMPI
   short  numlinkmusclvars;
   double *linkmusclvars;
#endif
#ifndef NDEBUG
     long i,j,k;
#endif
} np_t;

size_t wfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);

FILE *wfopen(const char *path, const char *mode);

int wfclose(FILE *fp);

int wfprintf(FILE *stream, const char *format, ...);

void fatal_error(const char *formatstr, ...);

long _i(long l, gl_t *gl, long dim);

long _i_all(long l, gl_t *gl, long dim);

void find_ijk_from_l(gl_t *gl, long l, long *i, long *j, long *k);

void find_ijk_from_l_all(gl_t *gl, long l, long *i, long *j, long *k);

long _l_from_l_all(gl_t *gl, long l_all);

void multiply_matrix_by_constant(double multfact, sqmat_t A);

void display_matrix(sqmat_t A);

void multiply_matrix_and_matrix(sqmat_t A, sqmat_t B, sqmat_t C);

void make_matrix_positive(sqmat_t A);

void make_matrix_diagonal(sqmat_t A);

void arithmetic_mean_of_two_matrices(sqmat_t A, sqmat_t B, sqmat_t C);

void harmonic_mean_of_two_matrices(sqmat_t A, sqmat_t B, sqmat_t C);

void multiply_diagonal_matrix_and_matrix(sqmat_t A, sqmat_t B, sqmat_t C);

void multiply_matrix_and_diagonal_matrix(sqmat_t A, sqmat_t B, sqmat_t C);

void set_matrix_to_identity(sqmat_t A);

void set_matrix_to_zero(sqmat_t A);

void set_matrix_to_matrix(sqmat_t A, sqmat_t B);

void set_vector_to_zero(flux_t S);

void add_two_matrices(sqmat_t A, sqmat_t B, sqmat_t C);

void subtract_two_matrices(sqmat_t A, sqmat_t B, sqmat_t C);

void copy_matrix(sqmat_t AA, sqmat_t BB);

void copy_vector(flux_t AA, flux_t BB);

void multiply_matrix_and_vector(sqmat_t sqmat, flux_t mat1, flux_t mat2);

void multiply_diagonal_matrix_and_vector(sqmat_t sqmat, flux_t mat1, flux_t mat2);

void invert_matrix(sqmat_t mat1, sqmat_t mat2);

void invert_diagonal_matrix(sqmat_t mat1, sqmat_t mat2);

void find_determinant(sqmat_t mat1, double *det);

long _ai(gl_t *gl, long i, long j, long k);

long _ai_all(gl_t *gl, long i, long j, long k);

long _al(gl_t *gl, long l, long theta, long offset);

long _al_check(gl_t *gl, long l, long theta, long offset);

long _all(gl_t *gl, long l, long theta1, long offset1, long theta2, long offset2);

long _alll(gl_t *gl, long l, long theta1, long offset1, long theta2, long offset2, long theta3, long offset3);

long _al_link(np_t *np, gl_t *gl, long llink, long lbdry, long offset, int TYPELEVEL);

long _l_plus_one(long l, gl_t *gl, long theta);

long _l_minus_one(long l, gl_t *gl, long theta);

long _l_plus_two(long l, gl_t *gl, long theta);

long _l_minus_two(long l, gl_t *gl, long theta);

bool find_l_of_nearest_inner_node(np_t *np, gl_t *gl, long l, int TYPELEVEL, long *linner);

long _nodes_from_bdry(np_t *np, gl_t *gl, long l, long theta, int TYPELEVEL);

long _nodes_from_bdry_limited(np_t *np, gl_t *gl, long l, long theta, int TYPELEVEL, long limit);

long _nodes_from_bdry_through_links_limited(np_t *np, gl_t *gl, long l, long theta, int TYPELEVEL, long limit);

long _nodes_between_link_and_bdry_limited(np_t *np, gl_t *gl, long llink, long lbdry, int TYPELEVEL, long limit);

void find_numerical_jacobian(np_t *np, long l, gl_t *gl, void(*FUNCT)(np_t, gl_t *, long, flux_t), long theta, sqmat_t Ak);

void find_numerical_jacobian_2(np_t *np, long l, gl_t *gl, void(*FUNCT)(np_t, gl_t *, flux_t), sqmat_t Ak);

void find_numerical_jacobian_3(np_t *np, long l, gl_t *gl, void(*FUNCT)(np_t *, gl_t *, long, flux_t),
                 sqmat_t Ak);

void find_homogeneous_jacobian(np_t *np, gl_t *gl, long l, flux_t F, sqmat_t A);

/* the arguments i,j,k are here only used for debugging purposes
   if ever the code jams on this node in particular */
void create_node(np_t *np, long i, long j, long k);

void suspend_node(np_t *np);

bool resume_node(np_t *np);

void dispose_node(np_t *np);

bool is_node_valid(np_t np, int TYPELEVEL);

bool is_node_resumed(np_t np);

bool is_node_suspended(np_t np);

bool is_node_inner(np_t np, int TYPELEVEL);

bool is_node_bdry(np_t np, int TYPELEVEL);

bool is_node_link(np_t np, int TYPELEVEL);

bool is_node_inner_near_bdry(np_t *np, gl_t *gl, long l, int TYPELEVEL);

bool is_node_inner_within_n_nodes_of_bdry(np_t *np, gl_t *gl, long l, long n, int TYPELEVEL);

long _node_type(np_t np, int TYPELEVEL);

long _node_link(np_t np, long cntlink, int TYPELEVEL);

long _num_node_link(np_t np, int TYPELEVEL);

bool is_node_in_zone(long i, long j, long k, zone_t zone);

bool is_node_in_zone_2(long l, gl_t *gl, zone_t zone);

bool is_zone_in_zone(zone_t zone1, zone_t zone2);

bool is_zone_intersecting_zone(zone_t zone1, zone_t zone2);

void find_subzones_in_zone_given_zonelength(long zonelength, zone_t zone, long *numzone, zone_t **subzones);

void find_subzones_in_zone_given_numsubzone(zone_t zone, long numsubzonedesired, long *numsubzone, zone_t **subzones);

bool is_node_in_domain_lim(long l, gl_t *gl);

zone_t _zone_from_point(long i, long j, long k);

zone_t _zone_expansion(zone_t zone, long expnum);

int find_zone_intersection(zone_t zone1, zone_t zone2, zone_t *zoneint);

zone_t _zone_intersection_old(zone_t zone1, zone_t zone2);

zone_t _zone_intersection(zone_t zone1, zone_t zone2);

zone_t _domain_lim_from_domain(zone_t domain, gl_t *gl);

void init_data_structure(np_t **np, gl_t *gl, zone_t domain, zone_t domain_all);

void init_data_structure_and_create_nodes(np_t **np, gl_t *gl, zone_t domain, zone_t domain_all);

void validate_data_structure(np_t *np,gl_t *gl);

#ifdef DISTMPI

int _node_rank(gl_t *gl, long i, long j, long k);

int MPI_Bcast_Node(void *message, int count, MPI_Datatype datatype, int root,
                   MPI_Comm comm,long i, long j, long k, gl_t *gl);

/* the intermittenly blocking send */
int MPI_IBsend(void* message, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm);

zone_t _domain_from_rank(int rank, gl_t *gl);

zone_t _domain_from_rank_mem(int rank, gl_t *gl);

zone_t _domain_lim_from_rank(int rank, gl_t *gl);

zone_t _domain_lim_from_rank_mem(int rank, gl_t *gl);

#endif//DISTMPI

double _delta(long r, long k);

void add_int_to_codex(SOAP_codex_t *codex, char *varname, long varvalue);

void add_string_to_codex(SOAP_codex_t *codex, char *varname, char *varvalue);

void add_double_to_codex(SOAP_codex_t *codex, char *varname, double varvalue);

void find_double_var_from_codex(SOAP_codex_t *codex, char *name, double *var);

void find_bool_var_from_codex(SOAP_codex_t *codex, char *name, bool *var);

void find_int_var_from_codex(SOAP_codex_t *codex, char *name, int *var);

void find_zone_from_argum(char *argum, long startpos, gl_t *gl, SOAP_codex_t *codex, zone_t *zone);

void thread_lock_node_set(np_t *np, long l, int threadtype);

void thread_lock_node_unset(np_t *np, long l, int threadtype);

void thread_lock_global_set(gl_t *gl, int threadtype);

void thread_lock_global_unset(gl_t *gl, int threadtype);

void thread_lock_set(thread_lock_t *lockvar, int threadtype);

void thread_lock_unset(thread_lock_t *lockvar, int threadtype);

void thread_lock_destroy(thread_lock_t *lockvar, int threadtype);

void thread_lock_init(thread_lock_t *lockvar, int threadtype);

double distance2_between_nodes(np_t *np, gl_t *gl, long l1, long l2);

void write_options_row ( FILE * outputfile, char *col1, char *col2, char *col3, int linewidth, int lengthcol1, int lengthcol2 );

double _smooth(np_t *np, gl_t *gl, long l, int TYPELEVEL, double(*FUNCT)(np_t));

double _smooth2(np_t *np, gl_t *gl, long l, int TYPELEVEL, double(*FUNCT)(np_t *, gl_t *, long));

double _smooth3(np_t *np, gl_t *gl, long l, long spec, long dim, int TYPELEVEL, double(*FUNCT)(np_t *, gl_t *, long, long, long));


#endif /* _COMMON_H */
