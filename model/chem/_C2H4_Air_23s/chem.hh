#define _CHEM_C2H4_AIR_23S_66R_ZETTERVALL
#define _CHEM_METHOD "C2H4-Air 23 species 66 reactions Zettervall [zettervall2017a] and Zettervall-Konnov"
#define _CHEM_ACTIONNAME "C2H4Air23s"

#define CHEM_NEUTRAL TRUE
#define ns 25
#define ncs 0

#define specH    0
#define specH2   1
#define specO    2
#define specO2   3
#define specOH   4
#define specH2O  5
#define specN2   6
#define specH2O2 7
#define specC2H  8
#define specC2H2 9
#define specC2H3 10
#define specC2H4 11
#define specC2H5 12
#define specCH   13
#define specCH2  14
#define specCH3  15
#define specCH4  16
#define specHO2  17
#define specCHO  18
#define specCH2O 19
#define specCH3O 20
#define specCO   21
#define specCO2  22
#define specNO  23
#define specN  24

#define SPECIES_ELECTRON 0
#define SPECIES_NEUTRAL  1
#define SPECIES_IONPLUS  2
#define SPECIES_IONMINUS 3


#define SMAP_CHO SMAP_HCO
#define SMAP_CH2O SMAP_HCHO
#define specHCO  specCHO

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
  SMAP_H2,  
  SMAP_O,  
  SMAP_O2,   
  SMAP_OH,  
  SMAP_H2O, 
  SMAP_N2,  
  SMAP_H2O2,
  SMAP_C2H, 
  SMAP_C2H2,
  SMAP_C2H3, 
  SMAP_C2H4,
  SMAP_C2H5,
  SMAP_CH, 
  SMAP_CH2, 
  SMAP_CH3, 
  SMAP_CH4, 
  SMAP_HO2, 
  SMAP_CHO, 
  SMAP_CH2O,
  SMAP_CH3O,
  SMAP_CO, 
  SMAP_CO2,
  SMAP_NO,
  SMAP_N
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
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL,
  SPECIES_NEUTRAL  
};

typedef struct {
  int CHEMMODEL;
  double TMIN_LINDEMANN;
  bool LINDEMANNREACTIONS;
} gl_model_chem_t;

