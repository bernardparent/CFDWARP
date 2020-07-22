#define _FLUID_METHOD "Favre-Reynolds Multispecies N2 Vibrational Energy Transport"
#define _FLUID_FAVREREYNOLDS
#define _FLUID_MULTISPECIES
#define _FLUID_NEUTRALSTRANSPORT
#define _FLUID_DIFFUSION TRUE
#define _FLUID_CONVECTION TRUE
#define _FLUID_SOURCE TRUE
#define _FLUID_N2VIBMODEL TRUE
#define _FLUID_EENERGY FALSE
#define _FLUID_ACTIONNAME "FavreReynoldsTv"
#define defaultinitvartypefluid 5
#define totalinitvartypefluid 6
#define FLUX_LAMBDA_CONVECTIVE 0

#define NODETYPE_BDRYOUTFLOW 1
#define nf (ns+nd+4)
#define fluxmom ns
#define fluxet (ns+nd)
#define fluxtke (ns+nd+1)
#define fluxpsi (ns+nd+2)
#define fluxev (ns+nd+3)

#define totalpostvarfluid (13+nd+ns)
#define totalinitvarfluid (nd+ns+5)
#define hbw_bdry_fluid 2
// hbw_mem_fluid is the neighbor-nodes bandwidth   on which the node fluid props depend when being resumed
#define hbw_mem_fluid 1

typedef double spec_t[ns];
typedef double chargedspec_t[ncs];
typedef double spec2_t[ns][ns];
typedef double flux_t[nf];

typedef struct {
  double k,psi,dPdrhoetstar,htstar,ev;
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
  bool RAPCOMP,TURBSOURCE,
       REACTING,ADD_ETA_TO_ETAT_WITHIN_QK;
  double kdiv,psidiv,zetaA1,zetaA3,zetaA2,Prt,Sct;
  double Pmin,Pmax,Tmin,Tmax,wmin,kmin,kmax,
         psimin,psimax,aref,Eref,evref,kref,psiref,
         Tvmin,Tvmax,Twmin,Twmax;
  int DILATDISSIP,TURBMODEL,N2VIBMODEL,TEMODEL;
} gl_model_fluid_t;

