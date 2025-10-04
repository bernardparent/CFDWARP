#define _CHEM_H2AIR3S0R
#define _CHEM_METHOD "H2 Air 3 species 0 reaction"

#define ns 3
#define ncs 0
 

#define specO2 0
#define specN2 1


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
  SMAP_H2
};


#define SPECIES_ELECTRON 0
#define SPECIES_NEUTRAL  1
#define SPECIES_IONPLUS  2
#define SPECIES_IONMINUS 3


const static long speciestype[ns] = {
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


typedef struct {
} gl_model_chem_t;
