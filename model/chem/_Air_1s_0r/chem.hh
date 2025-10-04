#define _CHEM_NONE
#define _CHEM_METHOD "Air 1 species 0 reaction"

#define ns 1
#define ncs 0
 


const static long smap[ns] = {
  SMAP_Air
};


#define SPECIES_ELECTRON 0
#define SPECIES_NEUTRAL  1
#define SPECIES_IONPLUS  2
#define SPECIES_IONMINUS 3


const static long speciestype[ns] = {
  SPECIES_NEUTRAL
};

typedef struct {
} gl_model_chem_t;
