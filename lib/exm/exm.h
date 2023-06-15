#ifndef _EXM_H
#define _EXM_H


#include <memdebug.h>
#include <math.h>

#define EXM_NUMDIFF_FORWARDEULER 0
#define EXM_NUMDIFF_MODIFIEDEULER 1
#define EXM_NUMDIFF_IMPROVEDEULER 2
#define EXM_NUMDIFF_RUNGEKUTTA 3
#define EXM_NUMINTEG_RECTANGLES 0
#define EXM_NUMINTEG_POLY2 1
//#define mod(a,b) ((a)%(b))
//the following is the correct mod function but a bit more expensive to compute
#define mod(a,b) ((((a)%(b))+(b))%(b))
#define avg(a,b) (0.5e0*(a+b))
//#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )
#define sign(a) ( ( (a) < 0 )  ?  (-1)   : (1) )
#define TRUE    1
#define FALSE   0
typedef unsigned char bool;

typedef double EXM_vec3D_t[3];
typedef double EXM_mat3x3_t[3][3];

#ifndef tan
  #define tan(a)   (sin(a)/cos(a))
#endif


#ifndef round
//  #define round(a) (floor(a+0.5e0))
  #define round(a) (a<0?ceil((a)-0.5):floor((a)+0.5))
#endif


#define sqr(a)   ((a)*(a))
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))
#define pi  3.14159265358979323846


/* to be used in conjunction with EXM_ai3 */
typedef struct {
  long is,ie,js,je,ks,ke;
} EXM_gl3D_t;

/* to be used in conjunction with EXM_ai2 */
typedef struct {
  long is,ie,js,je;
} EXM_gl2D_t;

/* to be used in conjunction with EXM_aim */
typedef struct {
  long numcol,numrow;
} EXM_glm_t;

/* to be used in conjunction with EXM_ai1 */
typedef struct {
  long is,ie;
} EXM_gl1D_t;

typedef struct {
  double *cont;
  EXM_glm_t glm;
} EXM_mat_t;

typedef struct {
  double val[4];   /*0-2:LHS    3:RHS*/
} EXM_tdmaline_t;   /*tdm line*/

typedef struct {
  double val[6];   /*0-4:LHS    5:RHS*/
} EXM_pdmaline_t;   /*pdm line*/



double powint(double x, long y);

double krodelta(long i, long j);

double rad(double angle);

double deg(double angle);

long EXM_aim(EXM_glm_t gl, long row, long col);

long EXM_ai3(EXM_gl3D_t gl3d, long i, long j, long k);

long EXM_ai2(EXM_gl2D_t gl2d, long i, long j);

long EXM_ai1(EXM_gl1D_t gl1d, long i);

long EXM_al3(EXM_gl3D_t gl, long l, long theta, long offset);

long EXM_all3(EXM_gl3D_t gl, long l, long theta1, long offset1, long theta2, long offset2);

double EXM_find_root_zero_in(double(*FUNCT)(void *, double), void *arg_ptr,
              double minval, double maxval,
              double relerr, double abserr, long *IFLAG);

double EXM_find_root_Newton_Raphson(
              double(*FUNCT)(void *, double), void *arg_ptr,
              double root_guess, double droot_init,
              double relerr, double abserr, long *IFLAG);

void EXM_solve_XDMA(double *xdma, EXM_gl2D_t gl);

void EXM_solve_TDMA(EXM_tdmaline_t *tdma, long numlines);

void EXM_solve_PDMA(EXM_pdmaline_t *pdma, long numlines);

long EXM_mi(long k, long line, long row, long col);

void EXM_solve_block_PDMA(double *AA, double *BB, double *CC, double *DD,
                    double *EE, double *RHS, long linemax, long k);

void EXM_solve_block_TDMA(double *AA, double *BB, double *CC, double *RHS,
                    long linemax,long k);

void EXM_solve_block_TDMA_and_check(double *AA, double *BB, double *CC, double *RHS,
                    long linemax,long k);

void EXM_init_matrix(EXM_mat_t *mat, long numrow, long numcol);

void EXM_free_matrix(EXM_mat_t *mat);

void EXM_init_random_matrix(EXM_mat_t *mat, long numrow, long numcol);

void EXM_display_matrix(EXM_mat_t mat);

void EXM_multiply_matrices(EXM_mat_t mat1, EXM_mat_t mat2, EXM_mat_t *matr);

void EXM_init_identity_matrix(EXM_mat_t *mat, long numrow, long numcol);

void EXM_invert_matrix(EXM_mat_t mat, EXM_mat_t *matinv);

void EXM_invert_matrix_partial_pivoting(EXM_mat_t mat, EXM_mat_t *matinv);

void EXM_invert_matrix_gaussian_elimination(EXM_mat_t mat, EXM_mat_t *matinv);

void EXM_invert_matrix_analytical(EXM_mat_t mat, EXM_mat_t *matinv);

/* find orthogonal vector to the three points pa,
   pb and pc in the xyz frame of reference */
void EXM_find_orthogonal_vector(EXM_vec3D_t pa, EXM_vec3D_t pb, EXM_vec3D_t pc,
                    EXM_vec3D_t orthovect);

/* find plane A*x+B*y+C*z=D from orthogonal vector
   and one point on plane p0 */
void EXM_find_plane(EXM_vec3D_t orthovect, EXM_vec3D_t p0,
               double *A, double *B, double *C, double *D);

/* find pp, a point on the plane A-B-C-D which also
   lies on the line composed of vect and p0 */
void EXM_find_point_in_plane_on_vector(EXM_vec3D_t vect, EXM_vec3D_t p0,
               double A, double B, double C, double D,
               EXM_vec3D_t pp);

/* the plane is defined by points pa,pb and pc while
   the point to be mirrored is pp_o. The
   mirrored point is pp_m; pp_p is the point on the plane
   nearest to pp_o, midway between pp_o and pp_m */
void EXM_mirror_point_wrt_plane(EXM_vec3D_t pa, EXM_vec3D_t pb,
                         EXM_vec3D_t pc, EXM_vec3D_t pp_o,
                         EXM_vec3D_t pp_m);

void EXM_display_vector(EXM_vec3D_t pp);

double EXM_dot_product(EXM_vec3D_t vec1, EXM_vec3D_t vec2);

void EXM_cross_product(EXM_vec3D_t vec1, EXM_vec3D_t vec2, EXM_vec3D_t prod);

double EXM_vector_magnitude(EXM_vec3D_t vec);

double EXM_angle_between_vectors(EXM_vec3D_t vec1, EXM_vec3D_t vec2);

void EXM_normalize_vector(EXM_vec3D_t vec, EXM_vec3D_t vecnorm);

void EXM_find_rotation_matrix(EXM_vec3D_t vec1, EXM_vec3D_t vec2, EXM_mat3x3_t R);

void EXM_multiply_matrix_vector(EXM_mat3x3_t mat, EXM_vec3D_t vec, EXM_vec3D_t res);

double EXM_area_quadrilateral(EXM_vec3D_t A, EXM_vec3D_t B, EXM_vec3D_t C, EXM_vec3D_t D);

/* numerically differentiate FUNCT (which returns dx/dt)
   from t1 to t2, starting from x1, and returning x2 */
double EXM_numerical_differentiation(double(*FUNCT)(void *, double, double), void *arg_ptr,
                   long METHOD, long n, double x1, double t1, double t2, long *error);

/* numerically integrate the function FUNCT(x) between
   x=x1 and x=x2, returning the value of the integral */
double EXM_numerical_integration(double(*FUNCT)(void *, double), void *arg_ptr,
                    long METHOD, long n, double x1, double x2, long *error);

/* numerically integrate the function FUNCT(x) between
   x=x1 and x=x2, returning the value of the integral */
void EXM_numerical_integration_vector(void(*FUNCT)(void *, double, EXM_vec3D_t vector), void *arg_ptr,
                        long METHOD, long n, double x1, double x2, long *error,
                        EXM_vec3D_t sumvector);

double EXM_f_from_line(long N, double *x, double *f, double thisx);

void EXM_find_spline(long N, double *x, double *f, double *b);

double EXM_f_from_spline(long N, double *x, double *f, double *b, double thisx);

double EXM_f_from_monotonespline(long N, double *x, double *f, double thisx);

/* insert str1 into str2 at the position pos; returns *str2
   make sure *str2 has enough memory allocated */
char *strins(char *str1, char *str2,  long pos);

/* add line breaks without breaking words with width the maximum number of characters per line */
char *strwrp(char *str, int width);

/* replace the first occurrence of orig by repl within str; returns *str */
char *strrep(char *str, char *orig, char *repl);

/* add indent spaces to second and subsequent lines; make sure str has enough memory allocated */
char *strind(char *str, int indent);

/* add line breaks without breaking words with width the maximum number of characters per line 
   and indent the number of indented characters (either negative or positive)  */
char *strwrpind(char *str, int width, int indent);

/* find the current terminal window size in characters */
void find_terminal_window_size(int *width, int *height);

double avg_harmonic(double arg1, double arg2);

int process_flag_string(int argc, char **argv, char *flag, char **arg);

int process_flag_int(int argc, char **argv, char *flag, int *arg);

int process_flag_long(int argc, char **argv, char *flag, long *arg);

int process_flag_double(int argc, char **argv, char *flag, double *arg);

int process_flag_int_multiple(int argc, char **argv, char *flag, int **arg);

int process_flag_double_multiple(int argc, char **argv, char *flag, double **arg);

int process_flag(int argc, char **argv, char *flag);

int find_remaining_options(int argc, char **argv, char **options);

double min3(double val1, double val2, double val3);

double max3(double val1, double val2, double val3);

double notzero(double val, double sub);

double minmod(double x, double y);

double minmod3(double x, double y, double z);

double maxmag(double x, double y);

void output_backtrace(void);


#ifndef NDEBUG
  #define assert_str(x) #x
  #define assert(x) if (!(x)) { printf("Assertion failed: (%s), function %s, file %s, line %d.\n", assert_str(x), __PRETTY_FUNCTION__, __FILE__, __LINE__); output_backtrace(); abort(); }
#else
  #define __ASSERT_VOID_CAST (void)
  #define assert(x)		(__ASSERT_VOID_CAST (0))
#endif


#endif /* _EXM_H */
