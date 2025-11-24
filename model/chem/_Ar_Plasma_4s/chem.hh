#define _CHEM_ARPLASMA4S
#define _CHEM_METHOD "Ar Plasma 4 species [rodriguez2025]"
#define _CHEM_PLASMA
#define _CHEM_ACTIONNAME "ArPlasma4s"
#define _AVERAGEDRATES
#define _AVERAGEDRATES_CHEM
 
#define numaveragedrates_chem 8

const static averagedrates_id_type averagedrates_chem_id[numaveragedrates_chem] = {
  "kf1",
  "kf1b",
  "kf2",
  "kf3",
  "Nekf1",
  "Nekf1b",
  "Nekf2",
  "Nekf3"
};

#define ns 5
#define ncs 2

#define speceminus 0
#define specArplus 1
#define specAr4S 2
#define specAr 3
#define specN2 4


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
  SMAP_Arplus,
  SMAP_Ar4S,
  SMAP_Ar,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_IONPLUS,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


typedef struct {
  int CHEMMODEL;
  bool QEISOURCETERMS,TE_FROM_TOWNSEND, EXCITED_STATES, SUPERELASTIC_COLLISIONS, RATE_AVERAGING;
  double EoverNmin;
} gl_model_chem_t;

