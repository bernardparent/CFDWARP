#define _CHEM_AIRPLASMA8S
#define _CHEM_METHOD "Air Plasma 9 species Macheret2007 [parent2007a], Rajendran2022, Bowersox2024"
#define _CHEM_PLASMA
#define _CHEM_ACTIONNAME "AirPlasma9s"

#define CHEM_NEUTRAL FALSE

#define ns 9
#define ncs 4

#define speceminus 0
#define specO2minus 1
#define specO2plus 2
#define specN2plus 3
#define specO 4
#define specN 5
#define specNO 6
#define specO2 7
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
  SMAP_O2minus,
  SMAP_O2plus,
  SMAP_N2plus,
  SMAP_O,
  SMAP_N,
  SMAP_NO,
  SMAP_O2,
  SMAP_N2
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_IONMINUS,
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
  bool QEISOURCETERMS;
  double Tminchem,Teminchem;
} gl_model_chem_t;

