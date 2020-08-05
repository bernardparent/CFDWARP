#define _CHEM_AIRPLASMA11S
#define _CHEM_METHOD "Air Plasma 11s Dunn-Kang [dunn1973a]"
#define _CHEM_PLASMA
#define _CHEM_ACTIONNAME "AirPlasma11s"

#define ns 11
#define ncs 6
#define CHEM_NEUTRAL FALSE

#define speceminus 0
#define specOplus 1
#define specNplus 2
#define specO2plus 3
#define specN2plus 4
#define specNOplus 5
#define specO 6
#define specN 7
#define specNO 8
#define specO2 9
#define specN2 10

const static long smap[ns] = {
  SMAP_eminus,
  SMAP_Oplus,
  SMAP_Nplus,
  SMAP_O2plus,
  SMAP_N2plus,
  SMAP_NOplus,
  SMAP_O,
  SMAP_N,
  SMAP_NO,
  SMAP_O2,
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
  SPECIES_IONPLUS,
  SPECIES_IONPLUS,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL
};


typedef struct {
  int CHEMMODEL;
} gl_model_chem_t;

