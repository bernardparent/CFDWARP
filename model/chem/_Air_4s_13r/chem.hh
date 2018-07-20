#define _CHEM_AIR4S13R
#define _CHEM_METHOD "Air Plasma 4 species 13 reactions"

#define CHEM_NEUTRAL TRUE

#define ns 4
#define ncs 0

#define specO2 0
#define specN2 1
#define specO 2
#define specN 3



#define SPECIES_ELECTRON 0
#define SPECIES_NEUTRAL  1
#define SPECIES_IONPLUS  2
#define SPECIES_IONMINUS 3

/*
Species ordering:
1. Electrons
2. Negative ions
3. Positive ions
4. Neutrals

If there are is no electron species, give speceminus a rank of -1
*/


const static long smap[ns] = {
  SMAP_O2,
  SMAP_N2,
  SMAP_O,
  SMAP_N
};

const static long speciestype[ns] = {
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


