#define _FLUID_METHOD "Navier-Stokes Perfect Gas"
#define _FLUID_NAVIERSTOKESPERFECT
#define _FLUID_NEUTRALSTRANSPORT
#define _FLUID_CONVECTION TRUE
#define _FLUID_DIFFUSION TRUE
#define _FLUID_SOURCE TRUE
#define _FLUID_ACTIONNAME "NavierStokesPerfect"
#define _FLUID_N2VIBMODEL FALSE
#define _FLUID_EENERGY FALSE
#define fluxmom 1
#define fluxet (1+nd)
#define defaultinitvartypefluid 1
#define totalinitvartypefluid 3

#define NODETYPE_BDRYOUTFLOW 1
#define nf (2+nd)
#define totalpostvarfluid (4+nd*2)
#define totalinitvarfluid (nd+2)
#define hbw_bdry_fluid 3
// hbw_mem_fluid is the neighbor-nodes bandwidth   on which the node fluid props depend when being resumed
#define hbw_mem_fluid 0
// hbw_mem_fluid_metrics is the neighbor-nodes-metrics bandwidth on which the node fluid props depend
#define hbw_mem_fluid_metrics 1

#define FLUX_LAMBDA_CONVECTIVE 0



typedef double spec_t[ns];
typedef double chargedspec_t[ncs];
typedef double spec2_t[ns][ns];
typedef double flux_t[nf];

typedef struct {
  dim_t V;
  double rho,a2,gamma,eta;
  double zetaA1,zetaA3,zetaA2;
  int EIGENVALCOND;
} jacvars_t;





typedef struct {
  double zetaA1,zetaA3,zetaA2;
  double Pmin,Pmax,Tmin,Tmax;
  double gamma,R,eta,kappa;
  bool AXISYMMETRIC;
} gl_model_fluid_t;

