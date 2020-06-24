#define _CHEM_CSAIRPLASMA
#define _CHEM_METHOD "Cesium Air Plasma 13 species 45 reactions Lenard"
#define _CHEM_PLASMA

#define CHEM_NEUTRAL FALSE

#define ns 13    /* number of species */
#define ncs 7    /* number of charged species (electrons, ions) */

#define speceminus 0
#define specOminus 1
#define specO2minus 2
#define specO2plus 3
#define specN2plus 4
#define specNOplus 5
#define specCsplus 6
#define specCs 7
#define specO 8
#define specN 9
#define specNO 10
#define specO2 11
#define specN2 12




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
  SMAP_Ominus,
  SMAP_O2minus,
  SMAP_O2plus,
  SMAP_N2plus,
  SMAP_NOplus,
  SMAP_Csplus,
  SMAP_Cs,
  SMAP_O,
  SMAP_N,
  SMAP_NO,
  SMAP_O2,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_IONMINUS,
  SPECIES_IONMINUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


