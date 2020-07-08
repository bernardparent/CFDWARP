#define _CHEM_H2_AIR_17S_LENARDSMITH
#define _CHEM_METHOD "H2 Air Plasma 17 species Lenard-Smith"

#define CHEM_NEUTRAL FALSE
#define ns 17
#define ncs 6

#define speceminus 0
#define specO2minus 1
#define specO2plus 2
#define specN2plus 3
#define specNOplus 4
#define specH2plus 5
#define specH 6
#define specO 7
#define specN 8
#define specNO 9
#define specOH 10
#define specH2O 11
#define specHO2 12
#define specH2O2 13
#define specH2 14
#define specO2 15
#define specN2 16



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


#define speceminus 0
#define specO2minus 1
#define specO2plus 2
#define specN2plus 3
#define specNOplus 4
#define specH2plus 5
#define specH 6
#define specO 7
#define specN 8
#define specNO 9
#define specOH 10
#define specH2O 11
#define specHO2 12
#define specH2O2 13
#define specH2 14
#define specO2 15
#define specN2 16

const static long smap[ns] = {
  SMAP_eminus,
  SMAP_O2minus,
  SMAP_O2plus,
  SMAP_N2plus,
  SMAP_NOplus,
  SMAP_H2plus,
  SMAP_H,
  SMAP_O,
  SMAP_N,
  SMAP_NO,
  SMAP_OH,
  SMAP_H2O,
  SMAP_HO2,
  SMAP_H2O2,
  SMAP_H2,
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




