#define _CHEM_ARPLASMA3S
#define _CHEM_METHOD "Ar Plasma 3 species [rodriguez2025]"
#define _CHEM_PLASMA
#define _CHEM_ACTIONNAME "ArPlasma3s"
#define _AVERAGEDRATES
#define _AVERAGEDRATES_CHEM

#define numaveragedrates_chem 2

const static averagedrates_id_type averagedrates_chem_id[numaveragedrates_chem] = {
  "kf1",
  "Nekf1"
};

#define ns 4
#define ncs 2

#define speceminus 0
#define specArplus 1
#define specAr 2
#define specN2 3


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
  SMAP_Ar,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_IONPLUS,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


typedef struct {
  int CHEMMODEL;
  bool QEISOURCETERMS,TE_FROM_TOWNSEND;
  double EoverNmin;
} gl_model_chem_t;

