#define _CHEM_AIRLOWTEMP
#define _CHEM_METHOD "air low temperature"

#define ns 15

#define specN2 10
#define specO 1
#define speceminus 0

const static long smap[ns] = {
  SMAP_eminus,
  SMAP_O,
  SMAP_Ominus,
  SMAP_Oplus,
  SMAP_O2,
  SMAP_O2minus,
  SMAP_O2plus,
  SMAP_O3,
  SMAP_N,
  SMAP_Nplus,
  SMAP_N2,
  SMAP_N2plus,
  SMAP_N2A,
  SMAP_NO,
  SMAP_NOplus
};


