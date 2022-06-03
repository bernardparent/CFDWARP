// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2018,2020 Bernard Parent
Copyright 2020 Aaron Trinh
Copyright 2001 Jason Etele
Copyright 2000 Giovanni Fusina


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
 * Refs.
 * Dixon-Lewis, G. "Computer Modeling of Combustion Reactions in Flowing Systems with Transport" in "Combustion Chemistry" edited by Gardiner, W.C. Springer-Verlag, NY, 1984.
 * */

#include <model/thermo/_thermo.h>
#include <model/transport/_transport.h>

// maximum reduced temperature (T/Peps) used to determine the viscosity, thermal conductivity, and mass diffusion coefficients with Leonard Jones potentials
#define TSTAR_MAX 100.0

#define INCLUDE_ELECTRONS_IN_ETA TRUE
#define INCLUDE_ELECTRONS_IN_NU FALSE
#define INCLUDE_ELECTRONS_IN_KAPPA TRUE




const static double Pd[3][5]=
 {
   {2.3527333E+0, -1.3589968E+0, 5.2202460E-1, -9.4262883E-2, 6.4354629E-3},
   {1.2660308E+0, -1.6441443E-1, 2.2945928E-2, -1.6324168E-3, 4.5833672E-5},
   {8.5263337E-1, -1.3552911E-2, 2.6162080E-4, -2.4647654E-6, 8.6538568E-9}
 };

const static double Pe[3][3]=
 {
   {1.1077725E+0, -9.4802344E-3, +1.6918277E-3},
   {1.0871429E+0, +3.1964282E-3, -8.9285689E-5},
   {1.1059000E+0, +6.5136364E-4, -3.4090910E-6}
 };


/*
transport properties: Lennard-Jones Potential Parameters: epsilon & sigma

Primary Reference: all species except those not mentioned in the secondary reference
Svehla, R.A,"Estimated Viscosities and Thermal Conductivites of Gases at high temperature,"
	NASA TR R-132, 1962
CFDWARP/model/thermo/_generic/ref/NASA-TR-R-132.pdf

Secondary Reference: HO2, HCO, C2H3, C2H5, C4H8, CH2O, CH3, CH3O, HNO, NO2
Sandia National Laboratories, "A Fortran Computer Code Package for the Evaluation of Gas Phase
	Multicomponent Transport Properties," SAND86-8246, 1986
Sandia National Laboratories: SAND86-8246 Update/Revision, 1998
CFDWARP/model/transport/_dixonlewis/ref/SAND86-8246.pdf

the transport properties of CHO, H2CO, C4H8O, C6H12O, C8H16, and C12H24
are calculated as           HCO, CH2O, C4H8, C2H5OC2H5, C6H12, and C6H12 respectively

*/

/*The epsilon parameter values below are sourced from one of two tables above and include
 * Boltzmann's constant inside already. Thus the units are in Kelvin. (epsilon/k_b)
 * 
 * Units of sigma: A/10 (dA, deciangstrom)
 * 
 * */

const static double Peps[SMAP_NS]=
  {
   106.7e0,     /* e- */ /* !! unknown value: fixed to the one of O */
   106.7e0,     /* O2 */
   71.4e0,      /* N2 */
   106.7e0,     /* O */
   71.4e0,      /* N */
   106.7e0,     /* O2+ */ /* !! unknown value: fixed to the one of O2 */
   71.4e0,      /* N2+ */ /* !! unknown value: fixed to the one of N2 */
   106.7e0,     /* O+ */  /* !! unknown value: fixed to the one of O */
   71.4e0,      /* N+ */  /* !! unknown value: fixed to the one of N */
   106.7e0,     /* O2- */ /* !! unknown value: fixed to the one of O2 */
   106.7e0,     /* O- */  /* !! unknown value: fixed to the one of O */
   106.7,       /* O3 */  /* !! unknown value: fixed to the one of O2 */
   71.4,        /* N2(A) */  /* !! unknown value: fixed to the one of N2 */
   116.7e0,     /* NO */
   116.7e0,     /* NO+ */  /* !! unknown value: fixed to the one of NO */
   59.7e0,      /* H2 */
   809.1e0,     /*H2O */
   37.0e0,      /*H */
   79.8e0,      /*OH*/
   107.4e0,     /*HO2*/
   91.7e0,      /*CO*/
   195.2e0,     /*CO2*/
   144.0e0,     /* CH3*/
   148.6e0,     /* CH4 */
   289.3E+0,    /*H2O2*/
   498.0E+0,    /*CHO */
   498.0E+0,    /*CH2O*/
   417.0E+0,    /*CH3O*/
   209.0E+0,    /*C2H3*/
   224.7E+0,    /*C2H4*/
   252.3E+0,    /*C2H5*/
   215.7E+0,    /*C2H6*/
   68.6E+0,     //CH from NASA TR R-132.pdf
   65.3E+0,     //NH from NASA TR R-132.pdf
   231.8E+0,    //C2H2
   237.1E+0,    //C3H8 from NASA TR R-132.pdf
   297.1E+0,    //C12H23 from NASA TR R-132.pdf
   10.22E+0,    //He
   78.46e0,     /* Air -> obtained from N2 and O2 assuming 4:1 ratio */
   59.7e0,      /* H2+ (unknown!, fixed to the one of H2) */
   71.4,        /* Cs (unknown!, fixed to one of N2) */
   71.4,        /* Cs+ (unknown!, fixed to one of N2*/
   93.3,        /*Ar*/
   30.6,        /*C*/ 
   78.8,        /*C2*/ 
   75.0,        /*CN*/
   232.4,       /*NCO*/  
   93.3,        /*Ar+*/  /* !! unknown value: fixed to the one of Ar */
   30.6,        /*C+*/  /* !! unknown value: fixed to the one of C */
   78.8,        /*C2+*/  /* !! unknown value: fixed to the one of C2 */
   75.0,        /*CN+*/  /* !! unknown value: fixed to the one of CN */
   91.7,        /*CO+*/  /* !! unknown value: fixed to the one of CO */
   116.7,        /*HNO*/
   200.0,        /*NO2*/
   37.0,        /* H+ */ /* !! unknown value: fixed to the one of H */
   144.0,         /* CH2 */
   10.22,         /* He+  !! unknown value: fixed to the one of He*/
  };

const static double Psig[SMAP_NS]=
  {
   0.34e0,      /* e- */   /* !! unknown value: fixed to the one of O */
   0.3467e0,    /* O2 */
   0.3798e0,    /* N2 */
   0.34e0,      /* O */
   0.3298e0,    /* N */
   0.84e0,      /* O2+ */  /* !! determined from the mobility at 300K */
   1.22e0,      /* N2+ */  /* !! determined from the mobility at 300K */
   0.34e0,      /* O+ */   /* !! unknown value: fixed to the one of O */
   0.3298e0,    /* N+ */   /* !! unknown value: fixed to the one of N */
   0.99e0,      /* O2- */  /* !! determined from the mobility at 300K */
   0.34e0,      /* O- */  /* !! unknown value: fixed to the one of O */
   0.3467e0,    /* O3 */   /* !! unknown value: fixed to the one of O2 */
   0.3798e0,    /* N2(A) */  /* !! unknown value: fixed to the one of N2 */
   0.3492e0,    /* NO */   
   0.3492e0,    /* NO+ */  /* !! unknown value: fixed to the one of NO */
   0.2827E+0,   /* H2 */
   0.2641E+0,   /*H2O */
   0.2708E+0,   /*H */
   0.3147E+0,   /*OH*/
   0.3458E+0,   /*HO2*/
   0.3690E+0,   /*CO*/
   0.3941E+0,   /*CO2*/
   0.3800E+0,   /* CH3*/
   0.3758E+0,   /* CH4 */
   0.4196E+0,   /*H2O2*/
   0.3590E+0,   /*CHO */
   0.3590E+0,   /*CH2O*/
   0.3690E+0,   /*CH3O*/
   0.4100E+0,   /*C2H3*/
   0.4163E+0,   /*C2H4*/
   0.4302E+0,   /*C2H5*/
   0.4443E+0,   /*C2H6*/
   0.3370E+0,   //CH from NASA TR R-132.pdf
   0.3312E+0,   //NH from NASA TR R-132.pdf
   0.4033E+0,   //C2H2
   0.5118E+0,   //C3H8 from NASA TR R-132.pdf
   0.6182E+0,   //C12H23 from NASA TR R-132.pdf     
   0.2551E+0,   //He
   0.3732e0,    /* Air -> obtained from N2 and O2 assuming a 4:1 ratio */
   0.2827E+0,   /* H2+ !! fixed to the one of H2 */
   0.3567E+0,   /* Cs, unknown!, fixed to the one of Na */
   0.3567E+0,   /* Cs+, unknown!, fixed to the one of Na */
   0.3542E+0,   /*Ar*/
   0.3385E+0,   /*C*/
   0.3913E+0,   /*C2*/
   0.3856E+0,   /*CN*/
   0.3828E+0,   /*NCO*/
   0.3543E+0,   /* Ar+ !! fixed to the one of Ar */   
   0.3385E+0,   /*C+ !! fixed to the one of C */
   0.3913E+0,   /*C2+* !! fixed to the one of C2 */
   0.3856E+0,   /*CN+* !! fixed to the one of CN */
   0.3690E+0,   /*CO+* !! fixed to the one of CO */  
   0.3492E+0,   /*HNO*/
   0.2500E+0,   /*NO2*/
   0.2708E+0,   /* H+  !! fixed to the one of H */ 
   0.3800E+0,   /* CH2 */
   0.2551E+0,   /* He+ !! fixed to the one of He */
  };

static double _Omega11(double T, double eps){
  double tmp;
  double Tstar;
  long index;

  assert(eps!=0.0e0);
  Tstar=T/eps;
  Tstar=min(TSTAR_MAX,Tstar);
  index=2;
  if (Tstar<10.0e0) index=1;
  if (Tstar<5.0e0) index=0;
  tmp=Pd[index][0]
     +Pd[index][1]*Tstar
     +Pd[index][2]*Tstar*Tstar
     +Pd[index][3]*Tstar*Tstar*Tstar
     +Pd[index][4]*Tstar*Tstar*Tstar*Tstar;
#ifndef NDEBUG
  if (!(tmp>0.0)){
    wfprintf(stderr,"\n  T=%E\n",T);
    wfprintf(stderr,"  index=%ld\n",index);
    wfprintf(stderr,"Pd[index][0]=%E\n",Pd[index][0]); 
    wfprintf(stderr,"Pd[index][1]=%E\n",Pd[index][1]); 
    wfprintf(stderr,"Pd[index][2]=%E\n",Pd[index][2]); 
    wfprintf(stderr,"Pd[index][3]=%E\n",Pd[index][3]); 
    wfprintf(stderr,"Pd[index][4]=%E\n",Pd[index][4]);
    fatal_error("Problem in function _Omega11() part of thermo.c. Temperature may be out of polynomial bounds."); 
  }
#endif
  return(tmp);
}


static double _Omega22(double T, double eps){
  double tmp;
  double Omega11,Astar;
  double Tstar;
  long index;

  Omega11=_Omega11(T,eps);
  assert(eps!=0.0e0);
  Tstar=T/eps;
  Tstar=min(TSTAR_MAX,Tstar);
  index=2;
  if (Tstar<10.0e0) index=1;
  if (Tstar<5.0e0) index=0;
  Astar=Pe[index][0]
       +Pe[index][1]*Tstar
       +Pe[index][2]*Tstar*Tstar;
#ifndef NDEBUG
  if (Astar<0.0){
    wfprintf(stderr,"\n  T=%E\n",T);
    wfprintf(stderr,"  index=%ld\n",index);
    wfprintf(stderr,"  Pe[index][0]=%E\n",Pe[index][0]); 
    wfprintf(stderr,"  Pe[index][1]=%E\n",Pe[index][1]); 
    wfprintf(stderr,"  Pe[index][2]=%E\n",Pe[index][2]);
    fatal_error("Problem in _Omega22() part of thermo.c. Temperature may be out of polynomial bounds.");
  }
#endif
  tmp=Astar*Omega11;
  return(tmp);
}




static double _etak_from_T(long spec, double T){
  double etak;
  etak=0.0;
  switch (speciestype[spec]) {
    case SPECIES_IONPLUS:
      etak=0.00118*sqrt(_m(spec))*pow(T,2.5e0);
    break;
    case SPECIES_ELECTRON:
      etak=6.35e-4*sqrt(_m(spec))*pow(T,2.5e0);
    break;
    default: 
      assert((Psig[smap[spec]]*Psig[smap[spec]]*_Omega22(T,Peps[smap[spec]]))!=0.0e0);
      assert((_calM(spec)*T)>=0.0e0);
      etak=8.44107E-7*sqrt(_calM(spec)*T)/(Psig[smap[spec]]*Psig[smap[spec]]
            *_Omega22(T,Peps[smap[spec]]));
#ifndef NDEBUG
      if (!(etak>0.0)){
        wfprintf(stderr,"\n  _Omega22(T,Peps[smap[spec]])=%E\n",_Omega22(T,Peps[smap[spec]])); 
        wfprintf(stderr,"  sqrt(_calM(spec)*T)=%E\n",sqrt(_calM(spec)*T));
        fatal_error("  Problem computing etak in _etak_from_T().");
      }
#endif
    break;
  }
  assert(etak>0.0);
  return(etak);
}





static double _kappak_from_T(long spec, double T){
  double etak,kappak;
  etak=_etak_from_T(spec,T);
  if (_numatoms(spec)>1) {
    kappak=15.0e0/4.0e0*calR/_calM(spec)*etak*
                  (0.115e0+0.354e0*_calM(spec)*_cpk_from_T(spec,T)/calR);
  } else {
    kappak=15.0e0/4.0e0*calR/_calM(spec)*etak;
  }  
  return(kappak);
}



double _kappa_from_rhok_T_Te(spec_t rhok, double T, double Te){
  long spec,k,l;
  spec_t kappak,chik,etak,w;
  double rho,chisum,kappamix,sum;

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
  
  for (spec=0; spec<ns; spec++){
    kappak[spec]=_kappak_from_T(spec,T); 
    etak[spec]=_etak_from_T(spec,T);
  }
  chisum=0.0e0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_KAPPA) 
      chisum+=w[spec]/_calM(spec);
  }
  assert(chisum!=0.0e0);
  for (spec=0; spec<ns; spec++){
    chik[spec]=w[spec]/_calM(spec)/chisum;
  }
  kappamix=0.0e0;
  for (k=0; k<ns; k++){
    if (speciestype[k]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_KAPPA) { 
      sum=0.0e0;
      for (l=0; l<ns; l++){
        if (l!=k && (speciestype[l]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_KAPPA)) { 
          assert(etak[l]!=0.0e0);
          assert((1.0e0+_calM(k)/_calM(l))*8.0e0>0.0e0);
          assert((etak[k]/etak[l])>0.0e0);
          sum=sum+chik[l]/sqrt((1.0e0+_calM(k)/_calM(l))*8.0e0)*
                pow(1.0e0+sqrt(etak[k]/etak[l])*pow(_calM(l)/_calM(k),0.25e0),2.0e0);
        }
      }	
      assert((chik[k]+1.0654e0*sum)!=0.0e0);
      kappamix=kappamix+chik[k]*kappak[k]/(chik[k]+1.0654e0*sum);
    }
  }

    
  return(kappamix);                   
}


/* find kappa for the neutral species mixture */
double _kappan_from_rhok_T_Te(spec_t rhok, double T, double Te){
  long spec,k,l;
  spec_t etak,kappak,chik,w;
  double chisum,kappamix,sum,sum2,rho;
  
  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
  
  for (spec=0; spec<ns; spec++){
    kappak[spec]=_kappak_from_T(spec,T); 
    etak[spec]=_etak_from_T(spec,T);
  }
  chisum=0.0e0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL) 
      chisum+=w[spec]/_calM(spec);
  }
  assert(chisum!=0.0e0);
  for (spec=0; spec<ns; spec++){
    chik[spec]=w[spec]/_calM(spec)/chisum;
  }
  kappamix=0.0e0;
  for (k=0; k<ns; k++){
    if (speciestype[k]==SPECIES_NEUTRAL) { 
      sum=0.0e0;
      for (l=0; l<ns; l++){
        if (l!=k && (speciestype[l]==SPECIES_NEUTRAL)) { 
          assert(etak[l]!=0.0e0);
          assert((1.0e0+_calM(k)/_calM(l))*8.0e0>0.0e0);
          assert((etak[k]/etak[l])>0.0e0);
          sum=sum+chik[l]/sqrt((1.0e0+_calM(k)/_calM(l))*8.0e0)*
                pow(1.0e0+sqrt(etak[k]/etak[l])*pow(_calM(l)/_calM(k),0.25e0),2.0e0);
        }
      }	
      assert((chik[k]+1.0654e0*sum)!=0.0e0);
      kappamix=kappamix+chik[k]*kappak[k]/(chik[k]+1.0654e0*sum);
    }
  }


  /* make an adjustment to *kappa when charged species are not included */
  sum=0.0;
  sum2=0.0;
  for (k=0; k<ns; k++){
    if (speciestype[k]==SPECIES_NEUTRAL) { 
      sum+=w[k]/_calM(k);
    }    
    sum2+=w[k]/_calM(k);
  }
  return(kappamix*sum/sum2);
}




double _eta_from_rhok_T_Te(spec_t rhok, double T, double Te){
  long spec,k,l;
  spec_t etak,chik,w;
  double chisum,etamix,sum,rho;

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
  
  for (spec=0; spec<ns; spec++){
    etak[spec]=_etak_from_T(spec,T);
  }
  chisum=0.0e0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_ETA) 
      chisum+=w[spec]/_calM(spec);
  }
  assert(chisum!=0.0e0);
  for (spec=0; spec<ns; spec++){
    chik[spec]=w[spec]/_calM(spec)/chisum;
  }
  etamix=0.0e0;
  for (k=0; k<ns; k++){
    if (speciestype[k]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_ETA) { 
      sum=0.0e0;
      for (l=0; l<ns; l++){
        if (l!=k && (speciestype[l]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_ETA)) { 
          assert(etak[l]!=0.0e0);
          assert((1.0e0+_calM(k)/_calM(l))*8.0e0>0.0e0);
          assert((etak[k]/etak[l])>0.0e0);
          sum=sum+chik[l]/sqrt((1.0e0+_calM(k)/_calM(l))*8.0e0)*
                pow(1.0e0+sqrt(etak[k]/etak[l])*pow(_calM(l)/_calM(k),0.25e0),2.0e0);
        }
      }	
      assert((chik[k]+sum)!=0.0e0);
      etamix+=chik[k]*etak[k]/(chik[k]+sum);
    }
  }
  return(etamix);

}



/* find eta for the neutral species mixture */
double _etan_from_rhok_T_Te(spec_t rhok, double T, double Te){
  long spec,k,l;
  spec_t etak,chik,w;
  double chisum,etamix,sum,sum2,rho;

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
  
  for (spec=0; spec<ns; spec++){
    etak[spec]=_etak_from_T(spec,T);
  }
  chisum=0.0e0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL) 
      chisum+=w[spec]/_calM(spec);
  }
  assert(chisum!=0.0e0);
  for (spec=0; spec<ns; spec++){
    chik[spec]=w[spec]/_calM(spec)/chisum;
  }
  etamix=0.0e0;
  for (k=0; k<ns; k++){
    if (speciestype[k]==SPECIES_NEUTRAL) { 
      sum=0.0e0;
      for (l=0; l<ns; l++){
        if (l!=k && (speciestype[l]==SPECIES_NEUTRAL)) { 
          assert(etak[l]!=0.0e0);
          assert((1.0e0+_calM(k)/_calM(l))*8.0e0>0.0e0);
          assert((etak[k]/etak[l])>0.0e0);
          sum=sum+chik[l]/sqrt((1.0e0+_calM(k)/_calM(l))*8.0e0)*
                pow(1.0e0+sqrt(etak[k]/etak[l])*pow(_calM(l)/_calM(k),0.25e0),2.0e0);
        }
      }	
      assert((chik[k]+sum)!=0.0e0);
      etamix+=chik[k]*etak[k]/(chik[k]+sum);
    }
  }
  /* make an adjustment to *eta when charged species are not included */
  sum=0.0;
  sum2=0.0;
  for (k=0; k<ns; k++){
    if (speciestype[k]==SPECIES_NEUTRAL) { 
      sum+=w[k]/_calM(k);
    }    
    sum2+=w[k]/_calM(k);
  }
  return(etamix*sum/sum2);

}


void find_nuk_from_rhok_T_Te_muk(spec_t rhok, double T, double Te, chargedspec_t muk, spec_t nuk){
  long spec,k,l;
  spec_t chik,w;
  double P,rho;
  double chisum,sum;
  double calD[ns][ns];

  rho=0.0;
  for (spec=0; spec<ns; spec++) rho+=rhok[spec];
  for (spec=0; spec<ns; spec++) w[spec]=rhok[spec]/rho;
  
  // first make sure the temperature is not out of polynomial bounds
  P=_P_from_w_rho_T(w,rho,T);
  assert(P!=0.0);
  for (spec=0; spec<ns; spec++) rhok[spec]=rho*w[spec];

  chisum=0.0e0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_NU) 
      chisum+=w[spec]/_calM(spec);
  }
  assert(chisum!=0.0e0);
  for (spec=0; spec<ns; spec++){
    chik[spec]=w[spec]/_calM(spec)/chisum;
  }
  
  for (k=0; k<ns; k++){
    for (l=0; l<ns; l++){
      if ((speciestype[k]!=SPECIES_ELECTRON && speciestype[l]!=SPECIES_ELECTRON) || INCLUDE_ELECTRONS_IN_NU) { 
        assert(Peps[smap[k]]*Peps[smap[l]]>0.0e0);
        assert(P!=0.0e0);
        assert(_Omega11(T,sqrt(Peps[smap[k]]*Peps[smap[l]])));
        assert((P*_Omega11(T,sqrt(Peps[smap[k]]*Peps[smap[l]]))
               *(Psig[smap[k]]+Psig[smap[l]])*(Psig[smap[k]]+Psig[smap[l]]))!=0.0e0);
        assert(T*T*T*(1.0e0/_calM(l)+1.0e0/_calM(k))>0.0e0);
        calD[k][l]=2.381112E-5*sqrt(T*T*T*(1.0e0/_calM(l)+1.0e0/_calM(k)))
                 /(P*_Omega11(T,sqrt(Peps[smap[k]]*Peps[smap[l]]))
                 *(Psig[smap[k]]+Psig[smap[l]])*(Psig[smap[k]]+Psig[smap[l]]));
      }
    }
  }
  
  
  for (k=0; k<ns; k++){
    sum=0.0e0;
    for (l=0; l<ns; l++){
      if (l!=k && (speciestype[l]!=SPECIES_ELECTRON || INCLUDE_ELECTRONS_IN_NU)) { 
        assert(calD[k][l]!=0.0e0);
        sum=sum+chik[l]/calD[k][l];
      }
    }
    assert((sum+1.0E-20)!=0.0e0);
    nuk[k]=rho*(1.0e0-chik[k])/(sum+1.0E-20);
  }
  
  adjust_nuk_using_mobilities_given_muk(rhok, T, Te, muk, nuk);

}




/* electron mobility as a function of electron temperature and mass densities*/
static double _mue_from_Nn_Ni_Te(double Nn, double Ni, double Te){
  double muen,mue,muei;
  Te=max(Te,300.0);
  muen=3.74E19*exp(33.5/sqrt(log(Te)))/Nn;
  muei=1.9e16*pow(Te,1.5)/Ni;
  mue=1.0/(1.0/muen+1.0/muei);
  return(mue);
}


/* the weakly-ionized ion mobility corresponds to mui=min(A*Ti^n,B*Estar^p)/Nn 
 * mi is the mass of the ion in kg*/
static double _mui_from_Nn_Ni_Ti(double Nn, double Ni, double A, double Ti, double n, double B, double Estar, double p, double mi){
  double muin,muii,mui;
  muin=min(A*pow(Ti,n),B*pow(Estar,p))/Nn;
  muii=14.3e0/sqrt(mi)*pow(Ti,1.5)/Ni;
  mui=1.0/(1.0/muin+1.0/muii);
  return(mui);
}



/* find the mobility of species k [m2/Vs] using the species temperature Tk [K] and electric field in the species reference frame Ek [V/m] 

  O2+, N2+, and NO+ are found from Sinnott, G., Golden, D. E., & Varney, R. N. (1968). Positive-Ion Mobilities in Dry Air. Physical Review, 170(1), 272–275. doi:10.1103/physrev.170.272 

  O2- is found from GOSHO, Y. AND HARADA, A., “A New Technique for Measuring Negative Ion Mobilities at Atmospheric Pressure,” Journal of Physics D, Vol. 16, 1983, pp. 1159–1166.
   
  H2+, Cs+, N+, O+, O- are approximated using Fig. 8 in THE MOBILITIES OF SMALL IONS THE ATMOSPHERE AND THEIR RELATIONSHIP by E. UNGETHUM, Aerosol Science, 1974, Vol. 5, pp. 25 37. 
*/ 
double _muk_from_rhok_T_Te_Ek(spec_t rhok, double T, double Te, double Ek, long k){
  double mu,Estar,N,Nn,Ni;
  long spec;
  mu=0.0;
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  Nn=0.0;
  for (spec=ncs; spec<ns; spec++){
    Nn+=rhok[spec]/_m(spec);
  }
  Nn=max(1e0,Nn);
  Ni=0.0;
  for (spec=0; spec<ncs; spec++){
    if (speciestype[spec]==SPECIES_IONPLUS) Ni+=rhok[spec]/_m(spec);
  }
  Ni=max(1e0,Ni);

  Estar=Ek/N;
#ifdef speceminus
  if (CHEM_NEUTRAL && k==speceminus){
    /* electrons */
    mu=_mue_from_Nn_Ni_Te(Nn,Ni,Te);
  } else {
#endif
    switch (smap[k]){
      case SMAP_eminus:
        mu=_mue_from_Nn_Ni_Te(Nn,Ni,Te);
      break;
      case SMAP_O2plus:  
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 1.18e23, T, -0.5, 3.61E12, Estar, -0.5,  _m(k));
      break;
      case SMAP_N2plus:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 0.75e23, T, -0.5, 2.03E12, Estar, -0.5,  _m(k));
      break;
      case SMAP_O2minus:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 0.97e23, T, -0.5, 3.56e19, Estar, -0.1,  _m(k));
      break;
      case SMAP_Ominus:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 1.4*0.97e23, T, -0.5, 1.4*3.56e19, Estar, -0.1,  _m(k));
      break;
      case SMAP_NOplus:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 1.62e23, T, -0.5, 4.47e12, Estar, -0.5,  _m(k));
      break;
      default:
        mu=_mui_from_Nn_Ni_Ti(Nn,Ni, 2.2e10/sqrt(_m(k)), T, -0.5, 0.55/sqrt(_m(k)), Estar, -0.5,  _m(k));
    }
#ifdef speceminus
  }
#endif
  return(mu);
}




static void find_dmue_from_Nn_Ni_Te(double Nn, double Ni, double Te, double *dmuedTe, spec_t dmuedrhok){
  double logTe,dmuendTe;
  double muen,kmuei,muei,dmuendNn,dmueidTe,dmueidNi,dmuedmuen,dmuedmuei;
  spec_t dmuendrhok,dmueidrhok;
  long spec;
  if (Te>300.0){
//  muei=1e16*pow(Te,1.5)/Ni;
    muen=3.74E19*exp(33.5/sqrt(log(Te)))/Nn;
    kmuei=1.9e16;
    muei=kmuei/Ni*pow(Te,1.5);
//    mue=1.0/(1.0/muen+1.0/muei);

    
    logTe=log(Te);
    dmuendTe=-6.2645E20*exp(33.5/sqrt(logTe))/(Te*pow(logTe,1.5e0))/Nn;
    dmuendNn=-3.74E19*exp(33.5/sqrt(logTe))/sqr(Nn);
    for (spec=0; spec<ns; spec++){
      dmuendrhok[spec]=0.0;
      if (spec>=ncs) dmuendrhok[spec]=dmuendNn*1.0/_m(spec);
    }

    dmueidTe=kmuei/Ni*1.5*pow(Te,0.5);
    dmueidNi=-kmuei/sqr(Ni)*pow(Te,1.5);
    for (spec=0; spec<ns; spec++){
      dmueidrhok[spec]=0.0;
      if (speciestype[spec]==SPECIES_IONPLUS) dmueidrhok[spec]=dmueidNi*1.0/_m(spec);
    }
    
 

    dmuedmuen=1.0/sqr(muen)/sqr(1.0/muen+1.0/muei);
    dmuedmuei=1.0/sqr(muei)/sqr(1.0/muen+1.0/muei);
    *dmuedTe=dmuendTe*dmuedmuen+dmueidTe*dmuedmuei;
    for (spec=0; spec<ns; spec++) dmuedrhok[spec]=dmuendrhok[spec]*dmuedmuen+dmueidrhok[spec]*dmuedmuei;
  } else {
    *dmuedTe=0.0;
    for (spec=0; spec<ns; spec++) dmuedrhok[spec]=0.0;
  }
}



/* the mobility corresponds to mui=min(A*Ti^n,B*Estar^p)/N */
static void find_dmui_from_Nn_Ni_Ti(double Nn, double Ni, double A, double Ti, double n, double B, double Estar, double p, double mi,  double *dmuidTi, spec_t dmuidrhok){
  long spec;
  double term1,term2;
  double dmuindTi,dmuiidNi,dmuidmuin,dmuidmuii,kmuii,muii,muin,dmuiidTi;
  spec_t dmuindrhok,dmuiidrhok;
  
  for (spec=0; spec<ns; spec++) dmuindrhok[spec]=0.0;

  muin=min(A*pow(Ti,n),B*pow(Estar,p))/Nn;
  term1=A*pow(Ti,n);
  term2=B*pow(Estar,p);
  if (term1<term2){
  	dmuindTi=term1*n/Ti/Nn;
    for (spec=ncs; spec<ns; spec++) dmuindrhok[spec]=term1;
  } else {
  	dmuindTi=0.0;
    for (spec=ncs; spec<ns; spec++) dmuindrhok[spec]=term2;
  }
  for (spec=ncs; spec<ns; spec++) dmuindrhok[spec]*=-1.0/(Nn*Nn)/_m(spec);

  kmuii=14.3/sqrt(mi);
  muii=kmuii/Ni*pow(Ti,1.5);
  dmuiidTi=kmuii/Ni*1.5*pow(Ti,0.5);
  dmuiidNi=-kmuii/sqr(Ni)*pow(Ti,1.5);
  for (spec=0; spec<ns; spec++){
    dmuiidrhok[spec]=0.0;
    if (speciestype[spec]==SPECIES_IONPLUS) dmuiidrhok[spec]=dmuiidNi*1.0/_m(spec);
  }


  dmuidmuin=1.0/sqr(muin)/sqr(1.0/muin+1.0/muii);
  dmuidmuii=1.0/sqr(muii)/sqr(1.0/muin+1.0/muii);
  *dmuidTi=dmuindTi*dmuidmuin+dmuiidTi*dmuidmuii;
  for (spec=0; spec<ns; spec++) dmuidrhok[spec]=dmuindrhok[spec]*dmuidmuin+dmuiidrhok[spec]*dmuidmuii;

}


void find_dmuk_from_rhok_Tk_Ek(spec_t rhok, double Tk, double Ek, long k, double *dmukdTk, spec_t dmukdrhok){
  double N,Nn,Ni,Ekstar;
  long spec;
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=rhok[spec]/_m(spec);
  Nn=0.0;
  for (spec=ncs; spec<ns; spec++){
    Nn+=rhok[spec]/_m(spec);
  }
  Nn=max(1e0,Nn);
  Ni=0.0;
  for (spec=0; spec<ncs; spec++){
    if (speciestype[spec]==SPECIES_IONPLUS) Ni+=rhok[spec]/_m(spec);
  }
  Ni=max(1e0,Ni);
  Ekstar=Ek/N;
#ifdef speceminus
  if (CHEM_NEUTRAL && k==speceminus){
    /* electrons */
    find_dmue_from_Nn_Ni_Te(Nn, Ni, Tk, dmukdTk, dmukdrhok);
  } else {
#endif
    for (spec=0; spec<ns; spec++) dmukdrhok[spec]=0.0;
    switch (smap[k]){
      case SMAP_eminus:
        find_dmue_from_Nn_Ni_Te(Nn, Ni, Tk, dmukdTk, dmukdrhok);
      break;
      case SMAP_O2plus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 1.18e23, Tk, -0.5, 3.61E12, Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_N2plus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 0.75e23, Tk, -0.5, 2.03E12, Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_O2minus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 0.97e23, Tk, -0.5, 3.56e19, Ekstar, -0.1, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_Ominus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 1.4*0.97e23, Tk, -0.5, 1.4*3.56e19, Ekstar, -0.1, _m(k),  dmukdTk, dmukdrhok);
      break;
      case SMAP_NOplus:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 1.62e23, Tk, -0.5, 4.47e12, Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
      break;
      default:
        find_dmui_from_Nn_Ni_Ti(Nn, Ni, 2.2e10/sqrt(_m(k)), Tk, -0.5, 0.55/sqrt(_m(k)), Ekstar, -0.5, _m(k),  dmukdTk, dmukdrhok);
    }
#ifdef speceminus
  }
#endif


}


void find_nuk_eta_kappak_muk(spec_t rhok, double T, double Te,
                   spec_t nuk, double *eta, double *kappan, chargedspec_t kappac, chargedspec_t muk){
  double Ek;
  long spec;
  
  *eta=_eta_from_rhok_T_Te(rhok,T,Te);

  Ek=0.0;
  for (spec=0; spec<ncs; spec++) {
    muk[spec]=_muk_from_rhok_T_Te_Ek(rhok, T, Te, Ek, spec);
    kappac[spec]=_kappac_from_rhok_Tk_muk(rhok, T, Te, muk[spec], spec);
  }
  *kappan=_kappan_from_rhok_T_Te(rhok, T, Te);
  find_nuk_from_rhok_T_Te_muk(rhok, T, Te, muk, nuk);
}


void find_nuk_eta_kappa(spec_t rhok, double T, double Te,
                   spec_t nuk, double *eta, double *kappa){
  chargedspec_t muk,kappac;
  double kappan;
  long spec;
  find_nuk_eta_kappak_muk(rhok, T, Te, nuk, eta, &kappan, kappac, muk);
  *kappa=kappan;
  for (spec=0; spec<ncs; spec++) *kappa+=kappac[spec];
}


