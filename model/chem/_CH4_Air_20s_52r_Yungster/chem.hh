#define _CHEM_CH4_AIR_20S_52R_YUNGSTER
#define _CHEM_METHOD "CH4-Air 20 species 52 reaction Yungster [yungster1994a]"

#define CHEM_NEUTRAL TRUE
#define ns 20
#define ncs 0

#define specO2 1
#define specN2 19
#define specO 0

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
  SMAP_O,
  SMAP_O2,
  SMAP_H2,
  SMAP_H2O,  
  SMAP_H,
  SMAP_OH,  
  SMAP_HO2,  
  SMAP_CO,
  SMAP_CO2,
  SMAP_CH3,
  SMAP_CH4,  
  SMAP_H2O2,
  SMAP_CHO,
  SMAP_CH2O,
  SMAP_CH3O,
  SMAP_C2H3,
  SMAP_C2H4,
  SMAP_C2H5,
  SMAP_C2H6,
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
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


