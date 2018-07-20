#define _RES_METHOD "Generic"
#define _RES_GENERIC

#include <cycle/resconv/.active/resconv.hh>
#include <cycle/restime/.active/restime.hh>
#include <cycle/ressource/.active/ressource.hh>

/*
hbw_resplasma_fluid is the half band width of the stencil (excluding center node) that depends on emfield mem properties
hbw_resconvplasma_fluid is defined in cycle/resconv and varies depending on the flux discretization scheme
2 is the half bandwidth of a TVD stencil used to discretize dDstarU term
*/
#ifdef _FLUID_PLASMA
#define hbw_resplasma_fluid max(hbw_resconvplasma_fluid,2)  
#else
#define hbw_resplasma_fluid 0
#endif

#define hbw_res_fluid max(hbw_resconv_fluid,hbw_resplasma_fluid)

