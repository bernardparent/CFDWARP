#define _CHEM_C3H8_AIR_16S_23R_KUNDU
#define _CHEM_METHOD "C3H8-Air 16 species 23 reactions Kundu [kundu1999a]"

#define CHEM_NEUTRAL TRUE
#define ns 16
#define ncs 0

#define   specH  0
#define   specO  1
#define   specN  2
#define   specH2 3
#define   specO2 4
#define   specN2 5
#define   specOH 6
#define   specCH 7
#define   specCO 8
#define   specNH 9
#define   specNO 10
#define   specCO2 11
#define   specHO2 12
#define   specH2O 13
#define   specC2H2 14
#define   specC3H8 15
//#define   specC12H23 15



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
  SMAP_H,
  SMAP_O,
  SMAP_N,
  SMAP_H2,
  SMAP_O2,
  SMAP_N2,
  SMAP_OH,
  SMAP_CH,
  SMAP_CO,
  SMAP_NH,
  SMAP_NO,
  SMAP_CO2,
  SMAP_HO2,
  SMAP_H2O,
  SMAP_C2H2,
  SMAP_C3H8
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
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};

