#define _CHEM_H2_AIR_13S
#define _CHEM_METHOD "H2-Air 13 species Jachimowski [jachimowski1988a]"
#define _CHEM_ACTIONNAME "H2Air13s"

#define CHEM_NEUTRAL TRUE
#define ns 13
#define ncs 0

#define specH2 0
#define specO2 1
#define specH 2
#define specO 3
#define specOH 4
#define specH2O 5
#define specHO2 6
#define specH2O2 7
#define specN 8
#define specNO 9
#define specNO2 10
#define specHNO 11
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
  SMAP_H2,
  SMAP_O2,
  SMAP_H,
  SMAP_O,
  SMAP_OH,
  SMAP_H2O,
  SMAP_HO2,
  SMAP_H2O2,
  SMAP_N,
  SMAP_NO,
  SMAP_NO2,
  SMAP_HNO,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};




typedef struct {
  int CHEMMODEL;
} gl_model_chem_t;
