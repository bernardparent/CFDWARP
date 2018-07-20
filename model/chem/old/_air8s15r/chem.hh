#define _CHEM_AIR8S15R
#define _CHEM_METHOD "air low temperature 8 species 15 reactions"

#define ns 8

#define speceminus 0
#define specO2 1
#define specN2 2
#define specO 3
#define specO2plus 5
#define specN2plus 6

#define SPECIES_ELECTRON 0
#define SPECIES_NEUTRAL  1
#define SPECIES_IONPLUS  2
#define SPECIES_IONMINUS 3


const static long smap[ns] = {
  SMAP_eminus,
  SMAP_O2,
  SMAP_N2,
  SMAP_O,
  SMAP_N,
  SMAP_O2plus,
  SMAP_N2plus,
  SMAP_O2minus
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONMINUS
};


