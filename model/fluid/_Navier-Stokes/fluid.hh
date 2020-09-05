#define _FLUID_METHOD "Navier-Stokes Multispecies"
#define _FLUID_NAVIERSTOKES
#define _FLUID_MULTISPECIES
#define _FLUID_NEUTRALSTRANSPORT
#define _FLUID_DIFFUSION TRUE
#define _FLUID_CONVECTION TRUE
#define _FLUID_SOURCE TRUE
#define _FLUID_N2VIBMODEL FALSE
#define _FLUID_EENERGY FALSE
#define _FLUID_ACTIONNAME "NavierStokes"
#define fluxmom ns
#define fluxet (ns+nd)
#define defaultinitvartypefluid 5
#define totalinitvartypefluid 6
#define FLUX_LAMBDA_CONVECTIVE 0

#define NODETYPE_BDRYOUTFLOW 1
#define nf (ns+nd+1)
#define totalpostvarfluid (8+nd+ns)
#define totalinitvarfluid (nd+ns+2)
#define hbw_bdry_fluid 2
// hbw_mem_fluid is the neighbor-nodes bandwidth   on which the node fluid props depend when being resumed
#define hbw_mem_fluid 0
// hbw_mem_fluid_metrics is the neighbor-nodes-metrics bandwidth on which the node fluid props depend
#define hbw_mem_fluid_metrics 1

typedef double spec_t[ns];
typedef double chargedspec_t[ncs];
typedef double spec2_t[ns][ns];
typedef double flux_t[nf];

typedef struct {
  double dPdrhoetstar,htstar;
  spec_t w,dPdrhok;
  dim_t V,dPdrhou;
  double rho,P,T;
  double zetaA1,zetaA3,zetaA2;
  int EIGENVALCOND;
} jacvars_t;




typedef struct {
#ifdef _2D
  bool AXISYMMETRIC;
#endif
  bool DIFFUSION,REACTING;
  double zetaA1,zetaA3,zetaA2;
  double Pmin,Pmax,Tmin,Tmax,wmin,aref,Eref,
         Twmin,Twmax;
} gl_model_fluid_t;

