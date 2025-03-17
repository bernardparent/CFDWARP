#define _EMFIELD_TEST
#define EMFIELD
#define _EMFIELD_METHOD "Test"
#define _EMFIELD_ACTIONNAME "EMFieldTest"

#define hbw_res_emfield 1
// hbw_mem_emfield is the neighbor-nodes bandwidth   on which the node emfield props depend when being resumed
#define hbw_mem_emfield 1
#define nfe (1)
#define totalinitvaremfield (1)
#define totalpostvaremfield (2+nd)
#define totalinitvartypeemfield 1
#define defaultinitvartypeemfield 1
#define tsemfcoeffhbw 2 

typedef double fluxemfield_t[nfe];



typedef struct {
  double sigmasolid,sigmafluid,epsilonr;
} gl_model_emfield_t;
