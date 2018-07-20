#define _CHEM_AIR10S19R
#define _CHEM_METHOD "air low temperature 10 species 19 reactions"

#define ns 10

#define specN2 2
#define specO 3
#define speceminus 0
#define specO2plus 5
#define specN2plus 6

const static long smap[ns] = {
  SMAP_eminus,
  SMAP_O2,
  SMAP_N2,
  SMAP_O,
  SMAP_N,
  SMAP_O2plus,
  SMAP_N2plus,
  SMAP_Oplus,
  SMAP_Nplus,
  SMAP_O2minus
};



