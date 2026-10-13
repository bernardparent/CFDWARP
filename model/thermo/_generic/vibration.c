// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2026 Bernard Parent


Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
   of conditions and the following disclaimer in the documentation and/or other
   materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
  Vibrational thermodynamics of the molecules:
  - vibdata[]: the vibrational constants of each molecule;
  - the set V of species with vibrationally excited states, given by VibSpecies() in the
    control file;
  - e_vk(T), c_vk(T) and the vibrational entropy of each species at equilibrium (Boltzmann
    distribution at T);
  - find_vibdist(), the closure that gives the distribution of the vibrational states of
    the species in V from the vibrational energy e_v (J per kg of vibrating species). All
    other files obtain the vibrational state through find_vibdist() and never use the
    Boltzmann formula directly, so that another distribution (Treanor-Gordiets) can later
    be added here only.
*/

#include <model/thermo/_thermo.h>
#include <soap.h>
#include "thermo.h"


/* hc/kB in K cm */
#define VIB_C2 1.4387769
/* below e_vk=VIB_EVLIM (J/kg), e_vk is extended linearly in T so that the inversion
   e_v -> T_v remains defined as e_v -> 0 */
#define VIB_EVLIM 0.1
/* the vibrating species are considered absent when their mass fraction is below this */
#define VIB_WVMIN 1.0e-10
#define VIBMODES_MAX 4
#define VIB_NEWTON_ITERMAX 100
#define VIB_NEWTON_TOL 1.0e-13


typedef struct {
  long g;                 /* degeneracy */
  double we, wexe, weye;  /* Dunham constants, cm-1 */
} vibmode_t;


typedef struct {
  long smap;          /* SMAP_N2, SMAP_CO, ... */
  char name[8];
  double D0;          /* dissociation energy from v=0, cm-1 */
  double muosc;       /* oscillator reduced mass, amu (diatomics) */
  long nmode;
  vibmode_t mode[VIBMODES_MAX];
  double xij[VIBMODES_MAX][VIBMODES_MAX]; /* anharmonicity of polyatomics, cm-1 */
  /* if not zero, the following replace the spacing of the first level of the first mode
     (in K) and the gas constant of the species (in J/kgK) obtained from the constants above;
     this is done for N2 only, so that the results of earlier versions are unchanged */
  double theta1, R;
} vibdata_t;


/* vibdata[]: one entry per molecule, giving in this order:
     smap    the species (SMAP_N2, SMAP_O2, ...);
     name    the name of the species, as written in VibSpecies();
     D0      the dissociation energy measured from the level v=0, in cm-1 (1 cm-1 corresponds
             to VIB_C2=1.4387769 K); for a polyatomic molecule, the energy needed to break its
             weakest bond (OC-O, H-OH, H-NH2, O-NO, O-O2, H-CH3);
     muosc   the reduced mass of the two atoms, in amu (diatomics only, 0 for polyatomics);
     nmode   the number of vibrational modes, a degenerate mode being counted once;
     mode[]  for each mode, {g, we, wexe, weye}:
               g     the degeneracy of the mode;
               we, wexe, weye
                     the Dunham constants omega_e, omega_e x_e and omega_e y_e, in cm-1, of the
                     vibrational term G(v)=we(v+1/2)-wexe(v+1/2)^2+weye(v+1/2)^3; for the
                     polyatomics, we is the observed fundamental frequency G(1)-G(0), in cm-1,
                     and wexe and weye are 0;
     xij     the anharmonic coupling constants x_ij between the modes of the polyatomics, in
             cm-1; all 0 for now;
     theta1  if not 0, replaces the characteristic temperature of the first mode, in K;
     R       if not 0, replaces the gas constant of the species, in J/(kg K).
   init_vibdata() obtains the characteristic temperature of each mode, in K, as
   theta=VIB_C2*(G(1)-G(0))=VIB_C2*(we-2*wexe+13/4*weye). D0, muosc and xij are not used by the
   Boltzmann distribution; they are kept for the distributions with anharmonic levels.
   For N2, theta1=3353 K and R=296.8 J/(kg K) are the values used by the N2 vibrational model
   of earlier versions (THERMO_THETAV and THERMO_RN2 in thermo.h of 26oct06a), so that their
   results are unchanged; the Dunham constants and calR/_calM() would give 3352.2 K and
   296.80 J/(kg K).
   Sources:
   - diatomics (N2, O2, NO, CO, H2, OH), ground electronic state: we, wexe, weye and D0 from
     K.P. Huber and G. Herzberg, Molecular Spectra and Molecular Structure IV: Constants of
     Diatomic Molecules, Van Nostrand Reinhold, 1979, as given in the NIST Chemistry WebBook
     (NIST Standard Reference Database 69); D0 given in eV there is converted to cm-1; muosc
     is obtained from the masses of the atoms 14N, 16O, 12C and 1H;
   - polyatomics (CO2, H2O, NH3, NO2, O3, CH4): g and we are the experimental fundamental
     frequencies listed in the NIST Computational Chemistry Comparison and Benchmark Database
     (CCCBDB, NIST Standard Reference Database 101); D0 is the energy of the weakest bond at 0 K,
     with the products in their ground state, from:
       CO2 -> CO+O    43976.12+-0.15 cm-1: S. Gong, P. Wang and Y. Mo, J. Phys. Chem. Lett. 15,
                      10842 (2024);
       H2O -> H+OH    41145.92+-0.12 cm-1: O.V. Boyarkin, M.A. Koshelev, O. Aseev,
                      P. Maksyutenko, T.R. Rizzo, N.F. Zobov, L. Lodi, J. Tennyson and
                      O.L. Polyansky, Chem. Phys. Lett. 568-569, 14 (2013);
       NH3 -> H+NH2   37115+-20 cm-1: D.H. Mordaunt, M.N.R. Ashfold and R.N. Dixon, J. Chem.
                      Phys. 104, 6460 (1996);
       NO2 -> O+NO    25128.57+-0.05 cm-1: R. Jost, J. Nygard, A. Pasinski and A. Delon,
                      J. Chem. Phys. 105, 1287 (1996);
       O3  -> O+O2    8563.5+-3.5 cm-1: B. Ruscic, Active Thermochemical Tables, ver. 1.112 of
                      the Thermochemical Network (2013), as given by S. Ndengue, R. Dawes,
                      X.-G. Wang, T. Carrington, Z. Sun and H. Guo, J. Chem. Phys. 144, 074302
                      (2016);
       CH4 -> H+CH3   432.373 kJ/mol=36143.6 cm-1: B. Ruscic, J. Phys. Chem. A 119, 7810
                      (2015). */
static const vibdata_t vibdata[]={
  /* each entry has the form
       {smap, name, D0, muosc, nmode,
          { {g, we, wexe, weye}, ...one per mode... },
          xij, theta1, R },
     with D0, we, wexe, weye and xij in cm-1, muosc in amu, theta1 in K and R in J/(kg K) */

  {SMAP_N2,    "N2",    78714.0,     7.00154,   1,
     { {1,  2358.57,    14.324,  -0.00226} },
     {{0.0}},   3353.0,   296.8 },

  {SMAP_O2,    "O2",    41260.0,     7.99746,   1,
     { {1,  1580.19,    11.98,    0.0474} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_NO,    "NO",    52400.0,     7.46844,   1,
     { {1,  1904.20,    14.075,   0.0} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_CO,    "CO",    89460.0,     6.85621,   1,
     { {1,  2169.81,    13.288,   0.0105} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_H2,    "H2",    36118.0,     0.50391,   1,
     { {1,  4401.21,   121.34,    0.8129} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_OH,    "OH",    35420.0,     0.94809,   1,
     { {1,  3737.761,   84.8813,   0.5409} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_CO2,   "CO2",   43976.1,     0.0,       3,
     { {1,  1333.0,      0.0,     0.0},
       {2,   667.4,      0.0,     0.0},
       {1,  2349.1,      0.0,     0.0} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_H2O,   "H2O",   41145.9,     0.0,       3,
     { {1,  3657.1,      0.0,     0.0},
       {1,  1594.7,      0.0,     0.0},
       {1,  3755.9,      0.0,     0.0} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_NH3,   "NH3",   37115.0,     0.0,       4,
     { {1,  3337.0,      0.0,     0.0},
       {1,   950.0,      0.0,     0.0},
       {2,  3444.0,      0.0,     0.0},
       {2,  1627.0,      0.0,     0.0} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_NO2,   "NO2",   25128.6,     0.0,       3,
     { {1,  1318.0,      0.0,     0.0},
       {1,   750.0,      0.0,     0.0},
       {1,  1618.0,      0.0,     0.0} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_O3,    "O3",     8563.5,     0.0,       3,
     { {1,  1103.0,      0.0,     0.0},
       {1,   701.0,      0.0,     0.0},
       {1,  1042.0,      0.0,     0.0} },
     {{0.0}},      0.0,     0.0 },

  {SMAP_CH4,   "CH4",   36143.6,     0.0,       4,
     { {1,  2917.0,      0.0,     0.0},
       {2,  1534.0,      0.0,     0.0},
       {3,  3019.0,      0.0,     0.0},
       {3,  1306.0,      0.0,     0.0} },
     {{0.0}},      0.0,     0.0 },
};

#define NUMVIBDATA ((long)(sizeof(vibdata)/sizeof(vibdata[0])))


/* species whose vibrational states are solved for as separate species by some chemical solvers */
typedef struct {
  long smapmolecule;
  long smapstate;
} vibstate_t;

static const vibstate_t vibstate[]={
  {SMAP_NH3, SMAP_NH3v},
  {SMAP_H2,  SMAP_H2v},
  {SMAP_H2,  SMAP_H2v2},
  {SMAP_H2,  SMAP_H2v3},
  {SMAP_N2,  SMAP_N2v1},
  {SMAP_N2,  SMAP_N2v2},
  {SMAP_N2,  SMAP_N2v3},
  {SMAP_N2,  SMAP_N2v4},
  {SMAP_N2,  SMAP_N2v5},
  {SMAP_N2,  SMAP_N2v6},
  {SMAP_N2,  SMAP_N2v7},
  {SMAP_N2,  SMAP_N2v8},
};

#define NUMVIBSTATE ((long)(sizeof(vibstate)/sizeof(vibstate[0])))


/* data derived from vibdata[] for the species of the chemical solver, set once by
   set_vibspecies() */
static bool VIBINIT=FALSE;
static bool VIBSPECIESSET=FALSE;
static bool isvib[ns];
static long numvib=0;
static long nmodek[ns];
static long gk[ns][VIBMODES_MAX];
static double thetak[ns][VIBMODES_MAX];
static double Rk[ns];
static double Tvlimk[ns];


static long _vibdata_index(long spec){
  long cnt,ret;
  ret=-1;
  for (cnt=0; cnt<NUMVIBDATA; cnt++){
    if (vibdata[cnt].smap==smap[spec]) ret=cnt;
  }
  return(ret);
}


bool is_vib_data_available(long spec){
  return(_vibdata_index(spec)>=0);
}


/* the list of species with vibrational data, e.g. "N2, O2, NO" */
void find_vib_data_species_list(char **list){
  long cnt;
  *list=(char *)realloc(*list,(NUMVIBDATA*10+2)*sizeof(char));
  strcpy(*list,"");
  for (cnt=0; cnt<NUMVIBDATA; cnt++){
    if (cnt>0) strcat(*list,", ");
    strcat(*list,vibdata[cnt].name);
  }
}


/* returns TRUE if the chemical solver includes some vibrational states of species spec as
   separate species; the name of the first one found is then given in statename */
bool is_vib_state_in_chem(long spec, char **statename){
  long cnt,spec2;
  bool FOUND;
  FOUND=FALSE;
  for (cnt=0; cnt<NUMVIBSTATE; cnt++){
    if (vibstate[cnt].smapmolecule==smap[spec]){
      for (spec2=0; spec2<ns; spec2++){
        if (smap[spec2]==vibstate[cnt].smapstate && !FOUND) {
          FOUND=TRUE;
          find_species_name(spec2, statename);
        }
      }
    }
  }
  return(FOUND);
}


static double _evk_mode(double gR, double theta, double T){
  return(gR*theta/(exp(theta/T)-1.0));
}


/* e_vk(T) without the linear extension */
static double _evk_Boltzmann(long spec, double T){
  long m;
  double ev;
  ev=0.0;
  for (m=0; m<nmodek[spec]; m++) ev+=_evk_mode((double)gk[spec][m]*Rk[spec],thetak[spec][m],T);
  return(ev);
}


static void init_vibdata(void){
  long spec,m,index,cnt;
  double Tlo,Thi,Tmid;
  for (spec=0; spec<ns; spec++){
    index=_vibdata_index(spec);
    nmodek[spec]=0;
    Rk[spec]=calR/_calM(spec);
    Tvlimk[spec]=0.0;
    if (index>=0){
      nmodek[spec]=vibdata[index].nmode;
      if (vibdata[index].R!=0.0) Rk[spec]=vibdata[index].R;
      for (m=0; m<nmodek[spec]; m++){
        gk[spec][m]=vibdata[index].mode[m].g;
        thetak[spec][m]=VIB_C2*(vibdata[index].mode[m].we-2.0*vibdata[index].mode[m].wexe
                              +13.0/4.0*vibdata[index].mode[m].weye);
      }
      if (vibdata[index].theta1!=0.0) thetak[spec][0]=vibdata[index].theta1;
      /* temperature at which e_vk=VIB_EVLIM */
      if (nmodek[spec]==1 && gk[spec][0]==1) {
        Tvlimk[spec]=thetak[spec][0]/log(Rk[spec]*thetak[spec][0]/VIB_EVLIM+1.0);
      } else {
        Tlo=1.0e-3;
        Thi=1.0e4;
        for (cnt=0; cnt<200; cnt++){
          Tmid=0.5*(Tlo+Thi);
          if (_evk_Boltzmann(spec,Tmid)>VIB_EVLIM) Thi=Tmid; else Tlo=Tmid;
        }
        Tvlimk[spec]=0.5*(Tlo+Thi);
      }
    }
  }
  VIBINIT=TRUE;
}


/* set the species with vibrationally excited states; this can be done only once, before
   any vibrational or enthalpy function is called with the species in V */
void set_vibspecies(bool *isvibspec){
  long spec;
  if (VIBSPECIESSET) fatal_error("set_vibspecies() can not be called more than once.");
  if (!VIBINIT) init_vibdata();
  numvib=0;
  for (spec=0; spec<ns; spec++) {
    isvib[spec]=isvibspec[spec];
    if (isvib[spec]) {
      if (!is_vib_data_available(spec)) fatal_error("Species %ld has no vibrational data in set_vibspecies().",spec);
      numvib++;
    }
  }
  VIBSPECIESSET=TRUE;
}


bool _isvib(long spec){
  if (!VIBSPECIESSET) fatal_error("The species with vibrationally excited states have not been set yet; VibSpecies() must be called in the Fluid() module of the control file before the vibrational energy of any species is needed.");
  return(isvib[spec]);
}


long _numvibspecies(void){
  if (!VIBSPECIESSET) fatal_error("The species with vibrationally excited states have not been set yet in _numvibspecies().");
  return(numvib);
}


static void check_spec(long spec){
  if (!VIBINIT) init_vibdata();
  if (nmodek[spec]==0) fatal_error("Species %ld has no vibrational data in vibration.c.",spec);
}


/* the characteristic vibrational temperature of the mode of species spec with the lowest
   energy, K */
double _thetavk(long spec){
  long m;
  double theta;
  check_spec(spec);
  theta=thetak[spec][0];
  for (m=1; m<nmodek[spec]; m++) theta=min(theta,thetak[spec][m]);
  return(theta);
}


/* vibrational energy of species spec at equilibrium at T, J/kg */
double _evk_from_T(long spec, double T){
  long m;
  double ev;
  check_spec(spec);
  ev=T/Tvlimk[spec]*VIB_EVLIM;
  if (ev>VIB_EVLIM) {
    ev=0.0;
    for (m=0; m<nmodek[spec]; m++) ev+=_evk_mode((double)gk[spec][m]*Rk[spec],thetak[spec][m],T);
  }
  return(ev);
}


/* d e_vk/dT at equilibrium at T, J/kgK */
double _cvk_from_T(long spec, double T){
  long m;
  double cv,theta,gR;
  check_spec(spec);
  if (T<Tvlimk[spec]) {
    cv=1.0/Tvlimk[spec]*VIB_EVLIM;
  } else {
    cv=0.0;
    for (m=0; m<nmodek[spec]; m++) {
      gR=(double)gk[spec][m]*Rk[spec];
      theta=thetak[spec][m];
      cv+=gR*sqr(theta)*exp(theta/T)/sqr(T)/sqr(exp(theta/T)-1.0);
    }
  }
  return(cv);
}


/* vibrational entropy of species spec at equilibrium at T, J/kgK; zero below the
   temperature where e_vk is extended linearly */
double _svk_from_T(long spec, double T){
  long m;
  double sv,theta,gR;
  check_spec(spec);
  sv=0.0;
  if (T>Tvlimk[spec]) {
    for (m=0; m<nmodek[spec]; m++) {
      gR=(double)gk[spec][m]*Rk[spec];
      theta=thetak[spec][m];
      sv+=gR*((1.0+1.0/(exp(theta/T)-1.0))*theta/T-log(exp(theta/T)-1.0));
    }
  }
  return(sv);
}


/* d s_vk/dT at equilibrium at T, J/kgK2 */
double _dsvk_dT_from_T(long spec, double T){
  long m;
  double dsvdT,theta,gR;
  check_spec(spec);
  dsvdT=0.0;
  if (T>Tvlimk[spec]) {
    for (m=0; m<nmodek[spec]; m++) {
      gR=(double)gk[spec][m]*Rk[spec];
      theta=thetak[spec][m];
      dsvdT+=gR*sqr(theta)*exp(theta/T)/sqr(exp(theta/T)-1.0)/T/T/T;
    }
  }
  return(dsvdT);
}


/* mass fraction of the vibrating species */
double _wv_from_w(spec_t w){
  long spec;
  double wv;
  wv=0.0;
  for (spec=0; spec<ns; spec++) if (_isvib(spec)) wv+=w[spec];
  return(wv);
}


/* vibrational energy per unit mass of vibrating species of a distribution with the
   characteristic temperature Tv; for VIBDISTTYPE_BOLTZMANN, this is the Boltzmann
   distribution at Tv; with Tv=T this gives the equilibrium vibrational energy whatever the
   distribution type; returns 0 if no species vibrates */
double _ev_from_w_Tv(int VIBDISTTYPE, spec_t w, double Tv){
  long spec;
  double wv,ev;
  ev=0.0;
  switch (VIBDISTTYPE){
    case VIBDISTTYPE_BOLTZMANN:
      wv=_wv_from_w(w);
      if (wv>0.0) {
        for (spec=0; spec<ns; spec++) if (_isvib(spec)) ev+=(w[spec]/wv)*_evk_from_T(spec,Tv);
      }
    break;
    default:
      fatal_error("VIBDISTTYPE=%d is not valid in _ev_from_w_Tv(): only VIBDISTTYPE_BOLTZMANN "
                  "is available.",VIBDISTTYPE);
  }
  return(ev);
}


/* Boltzmann distribution at a common Tv for all species in V: solve
   sum_k (rho_k/rhoV) e_vk(Tv) = ev for Tv */
static void find_vibdist_Boltzmann(spec_t rhok, double rhoV, double T, double ev, vibdist_t *vibdist){
  long spec,specV,cnt;
  double Tv,f,dfdTv,dTv,Tvlo,Tvhi,evlo;
  bool CONVERGED;

  specV=-1;
  for (spec=0; spec<ns; spec++) if (isvib[spec]) specV=spec;
  if (numvib==1 && nmodek[specV]==1) {
    /* one vibrating species with one mode: analytical inversion */
    if (ev>VIB_EVLIM) {
      Tv=thetak[specV][0]/log((double)gk[specV][0]*Rk[specV]*thetak[specV][0]/ev+1.0);
    } else {
      Tv=ev/VIB_EVLIM*Tvlimk[specV];
    }
    vibdist->Tv=Tv;
    vibdist->evk[specV]=ev;
    vibdist->cvk[specV]=_cvk_from_T(specV,Tv);
    vibdist->cv=vibdist->cvk[specV];
  } else {
    /* Newton iterations kept within a bracket; e_v(Tv) increases monotonically */
    Tvlo=0.0;
    Tvhi=1.0e7;
    evlo=0.0;
    for (spec=0; spec<ns; spec++) if (isvib[spec]) evlo+=rhok[spec]/rhoV*_evk_from_T(spec,Tvlo);
    if (ev<=evlo) {
      /* below the linear extension: e_v is linear in Tv */
      dfdTv=0.0;
      for (spec=0; spec<ns; spec++) if (isvib[spec]) dfdTv+=rhok[spec]/rhoV*_cvk_from_T(spec,0.0);
      Tv=ev/dfdTv;
    } else {
      Tv=max(1.0,T);
      cnt=0;
      do {
        f=-ev;
        dfdTv=0.0;
        for (spec=0; spec<ns; spec++) {
          if (isvib[spec]) {
            f+=rhok[spec]/rhoV*_evk_from_T(spec,Tv);
            dfdTv+=rhok[spec]/rhoV*_cvk_from_T(spec,Tv);
          }
        }
        if (f>0.0) Tvhi=min(Tvhi,Tv); else Tvlo=max(Tvlo,Tv);
        dTv=-f/dfdTv;
        if (Tv+dTv<=Tvlo || Tv+dTv>=Tvhi) dTv=0.5*(Tvlo+Tvhi)-Tv;
        Tv+=dTv;
        cnt++;
        CONVERGED=(fabs(dTv)<VIB_NEWTON_TOL*Tv);
      } while (!CONVERGED && cnt<VIB_NEWTON_ITERMAX);
      if (!CONVERGED) fatal_error("Couldn't find Tv from ev=%E J/kg in find_vibdist().",ev);
    }
    vibdist->Tv=Tv;
    vibdist->cv=0.0;
    for (spec=0; spec<ns; spec++) {
      if (isvib[spec]) {
        vibdist->evk[spec]=_evk_from_T(spec,Tv);
        vibdist->cvk[spec]=_cvk_from_T(spec,Tv);
        vibdist->cv+=rhok[spec]/rhoV*vibdist->cvk[spec];
      }
    }
  }
}


/* The distribution of the vibrational states of the species in V given
   Nk:    the number densities of all species, 1/m3
   T:     the translational temperature, K
   ev:    the vibrational energy per unit mass of the vibrating species, J/kg
   devdt: the time rate of change of ev, J/kgs (not used by the Boltzmann distribution)
   In vibdist, the derivatives cvk and cv are taken at constant composition. When the
   vibrating species are absent (V empty or their mass fraction below VIB_WVMIN), Tv=T and
   all derivatives are zero. */
void find_vibdist(gl_t *gl, int VIBDISTTYPE, spec_t Nk, double T, double ev, double devdt, vibdist_t *vibdist){
  long spec;
  double rho,rhoV;
  spec_t rhok;

  if (!VIBSPECIESSET) fatal_error("The species with vibrationally excited states must be set before calling find_vibdist().");
  if (VIBDISTTYPE!=VIBDISTTYPE_BOLTZMANN)
    fatal_error("VIBDISTTYPE=%d is not valid in find_vibdist(): only VIBDISTTYPE_BOLTZMANN "
                "is available.",VIBDISTTYPE);
  vibdist->VIBDISTTYPE=VIBDISTTYPE;
  vibdist->ev=ev;
  rho=0.0;
  rhoV=0.0;
  for (spec=0; spec<ns; spec++){
    rhok[spec]=Nk[spec]*_m(spec);
    rho+=rhok[spec];
    if (isvib[spec]) rhoV+=rhok[spec];
    vibdist->evk[spec]=0.0;
    vibdist->cvk[spec]=0.0;
  }
  vibdist->EMPTY=(numvib==0 || !(rhoV>VIB_WVMIN*rho));
  if (vibdist->EMPTY) {
    vibdist->Tv=T;
    vibdist->cv=0.0;
    for (spec=0; spec<ns; spec++) if (isvib[spec]) vibdist->evk[spec]=_evk_from_T(spec,T);
  } else {
    find_vibdist_Boltzmann(rhok, rhoV, T, ev, vibdist);
  }
  for (spec=0; spec<ns; spec++){
    if (isvib[spec]) vibdist->T1k[spec]=vibdist->Tv; else vibdist->T1k[spec]=T;
  }
}
