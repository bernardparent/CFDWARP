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

#define _CHEM_ARPLASMABALLISTIC20S
#define _CHEM_BALLISTIC
#define _CHEM_METHOD "Ar Plasma Ballistic Electron 20 species [parent2026]"
#define _CHEM_PLASMA
#define _CHEM_ACTIONNAME "ArPlasmaBallistic20s"

/* species 0 is the bulk electron group; species 1..nb are the ballistic groups,
   with 1..nb-1 the volumetric groups and nb the emission group */
#define nb 17
#define ns (nb+3)
#define ncs (nb+2)

#define speceminus 0
#define specegroup(j) (j)
#define specArplus (nb+1)
#define specAr (nb+2)

/* reactions of the ballistic electron model, used by find_W() and by find_Q_ballistic() */
#define REACTION1 TRUE  /* e0 + Ar -> Ar+ + e0 + e_j  (electron-impact ionization by the bulk) */
#define REACTION2 TRUE  /* Ar+ + Ar -> Ar+ + Ar+ + e_j  (ion-impact ionization) */
#define REACTION3 TRUE  /* e_b + Ar -> e0 + Ar  (elastic de-beaming) */
#define REACTION4 TRUE  /* e_b + Ar -> e_j + Ar*_k  (excitation cascade and scattering) */
#define REACTION5 TRUE  /* e_b + Ar -> Ar+ + e_j + e_j'  (ionization cascade and scattering) */
#define REACTION6 TRUE  /* e_b -> e0  (turning point) */
#define REACTION8 TRUE  /* Ar -> Ar+ + e_j  (ionization by an external electron beam) */

#define HAYASHI_NUMINEL 26  /* 25 excitation channels and 1 ionization channel */
#define HAYASHI_INEL_IZ 25  /* the ionization channel */

/* energy source terms given by find_Q_ballistic() to model/share/fluid_plasma_share.c */
#define numQballistic 4
#define QB_BORN 0
#define QB_BEAM 1
#define QB_RUN 2
#define QB_THRESHOLD 3

#define SPECIES_ELECTRON 0
#define SPECIES_NEUTRAL  1
#define SPECIES_IONPLUS  2
#define SPECIES_IONMINUS 3

/*
Species ordering:
1. Electrons (bulk group first, then the ballistic groups 1..nb)
2. Negative ions
3. Positive ions
4. Neutrals
*/

const static long smap[ns] = {
  SMAP_eminus,
  SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus,
  SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus,
  SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus, SMAP_eminus,
  SMAP_Arplus,
  SMAP_Ar
};

const static long speciestype[ns] = {
  SPECIES_ELECTRON,
  SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON,
  SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON,
  SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON, SPECIES_ELECTRON,
  SPECIES_IONPLUS,
  SPECIES_NEUTRAL
};

/* names used in place of speciesname[smap[spec]] (see find_species_name() in thermo.c) */
const static char * const speciesname_chem[ns] = {
  "e-",
  "e1-",  "e2-",  "e3-",  "e4-",  "e5-",  "e6-",
  "e7-",  "e8-",  "e9-",  "e10-", "e11-", "e12-",
  "e13-", "e14-", "e15-", "e16-", "e17-",
  "Ar+",
  "Ar"
};

typedef struct {
  double phigauge,phicathode,epsemit,Tborn,zeta,EcritfactReaction7,chi0,kturnmax,Qbeamext;
  bool UNIFORMLEVELS;
  bool BALLISTIC;          /* FALSE: no ballistic electrons, the model collapses to the fluid LEA model */
  double epslevel[nb+1];   /* levels in eV, epslevel[1..nb] (Eq. 2 of parent2026) */
  EXM_vec3D_t beamdir;     /* unit vector along which the ballistic groups travel */
  /* the following are computed from the above once the control file is read */
  double philo[nb+1],phihi[nb+1],phimid[nb+1],Deltaphi[nb+1],r[nb+1];
  double phiemit,epsborn,lossmax;
} gl_model_chem_t;
