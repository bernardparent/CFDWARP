#define _CHEM_AIR5S
#define _CHEM_METHOD "Air 5s Lenard [lenard1964a], Dunn-Kang [dunn1973a], Park [park1993a], Boyd [boyd2007a]"
#define _CHEM_ACTIONNAME "Air5s"

#define ns 5
#define ncs 0
 

#define specO 0
#define specN 1
#define specNO 2
#define specO2 3
#define specN2 4

const static long smap[ns] = {
  SMAP_O,
  SMAP_N,
  SMAP_NO,
  SMAP_O2,
  SMAP_N2
};


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


const static long speciestype[ns] = {
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


typedef struct {
  int CHEMMODEL;
} gl_model_chem_t;

