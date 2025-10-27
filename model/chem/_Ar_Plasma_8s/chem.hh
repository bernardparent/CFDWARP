#define _CHEM_ARPLASMA8S
#define _CHEM_METHOD "Ar Plasma 8 species [rodriguez2025]"
#define _CHEM_PLASMA
#define _CHEM_ACTIONNAME "ArPlasma8s"

 

#define ns 9
#define ncs 3

#define speceminus 0
#define specArplus 1
#define specAr2plus 2
#define specAr4S 3
#define specAr4P 4
#define specAr2star 5
#define specAr2 6
#define specAr 7
#define specN2 8


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
  SMAP_Ar2plus,
  SMAP_Ar4S,
  SMAP_Ar4P,
  SMAP_Ar2star,
  SMAP_Ar2,
  SMAP_Ar,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
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
  bool QEISOURCETERMS,TE_FROM_TOWNSEND,SUPERELASTIC_COLLISIONS;
  double EoverNmin;
} gl_model_chem_t;

