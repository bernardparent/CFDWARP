#define _CHEM_N2PLASMA3S3R
#define _CHEM_METHOD "N2 Plasma 3 species 3 reactions"
#define _CHEM_PLASMA

#define CHEM_NEUTRAL FALSE

#define ns 3
#define ncs 2



#define speceminus 0
#define specN2plus 1
#define specN2 2



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
  SMAP_eminus,
  SMAP_N2plus,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_IONPLUS,
  SPECIES_NEUTRAL
};

typedef struct {
} gl_model_chem_t;

