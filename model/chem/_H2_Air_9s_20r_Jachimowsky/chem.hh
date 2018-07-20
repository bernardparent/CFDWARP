#define _CHEM_H2_AIR_9S_20R_JACHIMOWSKY
#define _CHEM_METHOD "H2-Air 9 species 20 reactions Jachimowsky [jachimowsky1988a]"

#define CHEM_NEUTRAL TRUE
#define ns 9
#define ncs 0

#define specH2 0
#define specO2 1
#define specH 2
#define specO 3
#define specOH 4
#define specH2O 5
#define specHO2 6
#define specH2O2 7
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
  SMAP_H2,
  SMAP_O2,
  SMAP_H,
  SMAP_O,
  SMAP_OH,
  SMAP_H2O,
  SMAP_HO2,
  SMAP_H2O2,
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
  SPECIES_NEUTRAL
};


