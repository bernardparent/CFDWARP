#define _CHEM_AIRPLASMA12S42RDUNNKANGMACHERET
#define _CHEM_METHOD "Air Plasma 12 species 42 reactions Dung-Kang-Macheret [dunn1973a,parent2016a]"
#define _CHEM_PLASMA

#define CHEM_NEUTRAL FALSE

#define ns 12
#define ncs 7

#define speceminus 0
#define specO2minus 1
#define specOplus 2
#define specNplus 3
#define specO2plus 4
#define specN2plus 5
#define specNOplus 6
#define specO 7
#define specN 8
#define specNO 9
#define specO2 10
#define specN2 11




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
  SMAP_O2minus,
  SMAP_Oplus,
  SMAP_Nplus,
  SMAP_O2plus,
  SMAP_N2plus,
  SMAP_NOplus,
  SMAP_O,
  SMAP_N,
  SMAP_NO,
  SMAP_O2,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_IONMINUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


typedef struct {
} gl_model_chem_t;
