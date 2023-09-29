#define _CHEM_NH3PLASMA10S
#define _CHEM_METHOD "NH3 Plasma 10s Bavafa [bavafa2008]"
#define _CHEM_PLASMA
#define _CHEM_ACTIONNAME "NH3Plasma10s"
#define _CHEM_WEEE_TWOTEMPERATURE
#define _AVERAGEDRATES

#define numaveragedrates 2

const static long averagedrates_react[numaveragedrates] = {
  1,  
  2
};

#define ns 10
#define ncs 4
#define CHEM_NEUTRAL FALSE

#define speceminus 0
#define specNH2plus 1
#define specNH3plus 2
#define specNH4plus 3
#define specNH 4
#define specNH2 5
#define specNH3 6
#define specH 7
#define specNH4 8
#define specN2 9

const static long smap[ns] = {
  SMAP_eminus,
  SMAP_NH2plus,
  SMAP_NH3plus,
  SMAP_NH4plus,
  SMAP_NH,
  SMAP_NH2,
  SMAP_NH3,
  SMAP_H,
  SMAP_NH4,
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
  SPECIES_ELECTRON,
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


typedef struct {
  int CHEMMODEL;
  bool QEISOURCETERMS;
} gl_model_chem_t;


