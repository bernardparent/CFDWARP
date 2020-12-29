#define _CHEM_CO2PLASMA23S
#define _CHEM_METHOD "CO2 Plasma 23s Park [park1994a]"
#define _CHEM_PLASMA

#define CHEM_NEUTRAL FALSE

#define ns 23    /* number of species */
#define ncs 11    /* number of charged species (electrons, ions) */


#define speceminus   0
#define specArplus   1
#define specCplus    2
#define specNplus    3
#define specOplus    4
#define specC2plus   5
#define specN2plus   6
#define specO2plus   7
#define specCNplus   8
#define specCOplus   9
#define specNOplus   10
#define specC2       11
#define specAr       12
#define specC        13
#define specN        14
#define specO        15
#define specN2       16
#define specO2       17
#define specCN       18
#define specCO       19
#define specNO       20
#define specCO2      21
#define specNCO      22

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
  SMAP_Cplus,
  SMAP_Nplus,
  SMAP_Oplus,
  SMAP_C2plus,
  SMAP_N2plus,
  SMAP_O2plus,
  SMAP_CNplus,
  SMAP_COplus,
  SMAP_NOplus,  
  SMAP_C2,
  SMAP_Ar,
  SMAP_C,
  SMAP_N,
  SMAP_O,
  SMAP_N2,
  SMAP_O2,
  SMAP_CN,
  SMAP_CO,
  SMAP_NO,
  SMAP_CO2,
  SMAP_NCO
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
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
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


typedef struct {
} gl_model_chem_t;
