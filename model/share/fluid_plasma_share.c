// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015-2018,2020,2021,2026 Bernard Parent

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

#include <model/share/fluid_plasma_share.h>
#include <model/share/fluid_share.h>
#include <model/share/emfield_share.h>
#include <model/fluid/_fluid.h>
#include <model/thermo/_thermo.h>
#include <model/transport/_transport.h>
#include <model/emfield/_emfield.h>
#include <model/chem/_chem.h>
#include <model/metrics/_metrics.h>
#include <model/beam/_beam.h>
#include <model/fluid/.active/fluid.h>
#include <model/share/model_share.h>
#include <src/control.h>
#include <src/bdry.h>


#ifdef _FLUID_PLASMA

double _betan(long k){
  double betan;
  if (speciestype[k]==SPECIES_NEUTRAL) betan=1.0; else betan=0.0;
  return(betan);
}

double _betai(long k){
  double betai;
  if (speciestype[k]==SPECIES_IONMINUS || speciestype[k]==SPECIES_IONPLUS) betai=1.0; else betai=0.0;
  return(betai);
}



/* determines whether the gauss's law terms are added to the continuity equations */
double _betag(gl_t *gl, long k){  
  return(gl->model.fluid.betag[k]);
}


double _betaplus(long k){
  double betaplus;
  if (speciestype[k]==SPECIES_IONPLUS) betaplus=1.0; else betaplus=0.0;
  return(betaplus);
}


double _betaa(gl_t *gl, long k){
  return(gl->model.fluid.betaa[k]);
}


double _rhoc(np_t np, gl_t *gl){
  double rhoc;
  long spec;
  rhoc=0.0;
  for (spec=0; spec<ns; spec++){
    rhoc+=_C(spec)*_Nk(np, gl, spec);
  }

  return(rhoc);
}


#define Ethresholdbdry 1e-2

bool is_bdry_cathode(np_t *np, gl_t *gl, long lA, long lB, long lC,
                     long theta, long thetasgn){
  bool CATHODE;
  double dphi,dwall,fact_cathode,Nminus,Nplus;
  long dim,spec;

  CATHODE=FALSE;

  dwall=0.0;
  for (dim=0; dim<nd; dim++){
    dwall=dwall+sqr(np[lB].bs->x[dim]-np[lA].bs->x[dim]);
  }
  assert_np(np[lA],dwall>0.0e0);
  dwall=sqrt(dwall);


  Nminus=0.0;
  Nplus=0.0;
  for (spec=0; spec<ns; spec++){
    if (_Charge_number(spec)<0) Nminus+=_Nk(np[lA],gl,spec);
    if (_Charge_number(spec)>0) Nplus+=_Nk(np[lA],gl,spec);
  }
  if (Nplus>Nminus) fact_cathode=1.0; else fact_cathode=0.0;

// Approach 1
/*
  dphi=minmod(_phi(np[lA],gl)-_phi(np[lB],gl),_phi(np[lB],gl)-_phi(np[lC],gl));
*/


// Approach 2
  dphi=2.0*(_phi(np[lA],gl)-_phi(np[lB],gl));
  dphi+=_phi(np[lB],gl)-_phi(np[lC],gl);
  for (dim=0; dim<nd; dim++){
    if (dim!=theta) dphi+=1.0*(_phi(np[_al(gl,lA,dim,+1)],gl)-_phi(np[_al(gl,lB,dim,+1)],gl));
    if (dim!=theta) dphi+=1.0*(_phi(np[_al(gl,lA,dim,-1)],gl)-_phi(np[_al(gl,lB,dim,-1)],gl));
  }
  dphi/=(3.0+((double)nd-1.0)*2.0);


// Approach 3

/*  double dphi1,dphi2;
  dphi1=2.0*(_phi(np[lA],gl)-_phi(np[lB],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta) dphi1+=1.0*(_phi(np[_al(gl,lA,dim,+1)],gl)-_phi(np[_al(gl,lB,dim,+1)],gl));
    if (dim!=theta) dphi1+=1.0*(_phi(np[_al(gl,lA,dim,-1)],gl)-_phi(np[_al(gl,lB,dim,-1)],gl));
  }
  dphi1/=(2.0+((double)nd-1.0)*2.0);

  dphi2=2.0*(_phi(np[lB],gl)-_phi(np[lC],gl));
  for (dim=0; dim<nd; dim++){
    if (dim!=theta) dphi2+=1.0*(_phi(np[_al(gl,lB,dim,+1)],gl)-_phi(np[_al(gl,lC,dim,+1)],gl));
    if (dim!=theta) dphi2+=1.0*(_phi(np[_al(gl,lB,dim,-1)],gl)-_phi(np[_al(gl,lC,dim,-1)],gl));
  }
  dphi2/=(2.0+((double)nd-1.0)*2.0);
  dphi1=minmod(dphi1,dphi2);
  dphi2=minmod(_phi(np[lA],gl)-_phi(np[lB],gl),_phi(np[lB],gl)-_phi(np[lC],gl));
  dphi=minmod(dphi1,dphi2);
*/


  if (dphi<Ethresholdbdry*dwall*fact_cathode) CATHODE=TRUE;


  return(CATHODE);
}



/* conditioned conductivity (must not be zero because some terms are divided by sigma_positive) */ 
double _sigma_positive(np_t *np, gl_t *gl, long l){
  double sigma_positive,sigma;
  sigma=_sigma(np,gl,l);
  sigma_positive=max(sigma,gl->model.fluid.sigmadiv);
  return(sigma_positive);
}



double _alpha(np_t *np, gl_t *gl, long l, long k, long r){
  double alpha,sigma_positive;
  sigma_positive=_sigma_positive(np,gl,l);
  alpha=_m(k)/_m(r)*(_delta(r,k)+_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_Nk(np[l],gl,k)/sigma_positive);
  return(alpha);
}


void find_alpha(np_t *np, gl_t *gl, long l, spec2_t alpha){
  long k,r;
  double sigma;
  sigma=_sigma_positive(np,gl,l);
  for (k=0; k<ns; k++){
    for (r=0; r<ns; r++){
      alpha[k][r]=_calM(k)/_calM(r)*(_delta(r,k)+_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_Nk(np[l],gl,k)/sigma);
      
    }
#ifndef NDEBUG
    if (alpha[k][k]<=0.0){
      wfprintf(stderr,"alpha[%ld][%ld]=%E   \n",k,k,alpha[k][k]);
      wfprintf(stderr,"mu_k=%E N_k=%E betaa_k=%E C_k=%E \n",_mu(np,gl,l,k),_Nk(np[l],gl,k),_betaa(gl,k),_C(k));
      wfprintf(stderr," sigma=%E betaa_k*C_k*mu_k*N_k=%E\n",sigma,_betaa(gl,k)*_C(k)*_mu(np,gl,l,k)*_Nk(np[l],gl,k));
      fatal_error("alpha[k][k] is negative at node i=%ld j=%ld k=%ld.",_i(l,gl,0),_i(l,gl,1),_i(l,gl,2));
    }
#endif
  }
  
}


void find_Z(np_t *np, gl_t *gl, long l, sqmat_t Z) {
  long row,col;
  spec2_t alpha;
  find_alpha(np, gl, l, alpha);
  set_matrix_to_identity(Z);
  for (row=0; row<ns; row++){
    for (col=0; col<ns; col++){
      Z[row][col]=alpha[row][col];
    }
  }

}


void find_dZU_dU(np_t *np, gl_t *gl, long l, sqmat_t dZU_dU){
  double sigma,dsigma_drhom;
  long k,r,m;
  sigma=_sigma_positive(np,gl,l);

  set_matrix_to_identity(dZU_dU);


  for (k=0; k<ns; k++){
    for (r=0; r<ns; r++){
      /* wrt rhok[k] */
      dZU_dU[k][k]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)/sigma*_rhok(np[l],r);
      /* wrt rhok[r] */
      dZU_dU[k][r]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_rhok(np[l],k)/sigma;
      /* wrt sigma */
      for (m=0; m<ns; m++){
        dsigma_drhom=_C(m)*_mu(np,gl,l,m)/_m(m);
        dZU_dU[k][m]+=-1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_rhok(np[l],k)/sqr(sigma)*_rhok(np[l],r)*dsigma_drhom; 
      }
    }
  }

}


void find_LambdaZ(np_t *np, gl_t *gl, long l, sqmat_t LambdaZ) {
  long r,k;
  double sigma;
  sigma=_sigma_positive(np,gl,l);
  set_matrix_to_identity(LambdaZ);
  for (k=0; k<ns; k++){
    for (r=0; r<ns; r++){
      LambdaZ[k][k]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)/sigma*_rhok(np[l],r);
    }
  }

}



void find_Y1star_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, flux_t Ystar){
  double Estar;
  long spec, flux;
  metrics_t metrics;
  find_metrics_at_interface(np,gl,lL,lR,theta,&metrics);
  find_Estar_at_interface(np, gl, lL, lR, theta, &Estar);

  for (flux=0; flux<nf; flux++){
    Ystar[flux]=0.0;
  }

  for (spec=0; spec<ns; spec++){
    Ystar[spec]=metrics.Omega*Estar*_betag(gl,spec);  
  }

}



void find_Y2star_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, flux_t Ystar){
  double Jstar;
  long spec, flux;
  metrics_t metrics;
  find_metrics_at_interface(np,gl,lL,lR,theta,&metrics);
  find_Jstar_at_interface(np, gl, lL, lR, theta, &Jstar);

  for (flux=0; flux<nf; flux++){
    Ystar[flux]=0.0;
  }

  for (spec=0; spec<ns; spec++){
    Ystar[spec]=-metrics.Omega*Jstar*_betaa(gl,spec);  
  }

}


/* drift velocity of species spec towards the boundary plane at lA, extrapolated from lB and lC; positive towards the wall */
double _Vk_drift_perp_to_bdry(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, long spec){
  EXM_vec3D_t Vk,Vk_m1,Vk_m2;
  double Vkperp;
  long dim;
  dim_t n;
  assert(spec<ncs);
  find_Vk_drift(np, gl, lC, spec, Vk_m2);
  find_Vk_drift(np, gl, lB, spec, Vk_m1);
  for (dim=0; dim<3; dim++) Vk[dim]=-Vk_m2[dim]+2.0*Vk_m1[dim];
  if (!find_unit_vector_normal_to_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK, n)) fatal_error("Problem finding n in _Vk_drift_perp_to_bdry() in fluid_plasma_share.c");
  Vkperp=0.0;
  for (dim=0; dim<nd; dim++) Vkperp-=Vk[dim]*n[dim];
  return(Vkperp);
}


#endif



/* Ballistic electron model of Parent and Perry, arXiv:2610.02139, 2026 [parent2026]: 
   the parts of the model that do not depend on the gas, shared by the plasma fluid modules; the reactions and their
   energy source terms are in the chem module; see fluid_plasma_share.h for what the chem and fluid modules must provide */

#if defined(_CHEM_BALLISTIC) && defined(_FLUID_PLASMA)

/* unit vector along which the ballistic groups travel */
static void find_beamdir(gl_t *gl, EXM_vec3D_t beamdir){
  long dim;
  for (dim=0; dim<3; dim++) beamdir[dim]=gl->model.chem.beamdir[dim];
}


/* dimensionless offset xi_j of the mean group potential, Eqs. (A.4) and (A.5) */
static double _xi(gl_t *gl, long j, double phi){
  double y,z,r,xi;
  y=(phi-gl->model.chem.phimid[j])/gl->model.chem.Deltaphi[j];
  /* group 1 has no less energetic window: phi should not exceed phigauge; if it does, xi is held */
  if (j==1) y=min(y,0.5);
  r=gl->model.chem.r[j];
  z=y-0.5;
  if (j<nb-1){
    if (y<=-0.5) xi=-0.5;
      else if (y<=0.5) xi=(4.0*y-1.0)/6.0;
      else if (z<=r) xi=(r+6.0*r*z+6.0*r*z*z-3.0*z*z-4.0*z*z*z)/(6.0*(r+2.0*r*z-z*z));
      else xi=(1.0+2.0*r)/6.0;
  } else {
    if (y<=-0.5) xi=-0.5;
      else if (y<=0.5) xi=(2.0*y-1.0)/4.0;
      else if (z<=r) xi=(6.0*r*z+6.0*r*z*z-3.0*z*z-4.0*z*z*z)/(6.0*(2.0*r+2.0*r*z-z*z));
      else xi=r*(3.0+2.0*r)/(6.0*(2.0+r));
  }
  return(xi);
}


/* mean group potential in V, Eq. (18), and phi_nb=phi_emit for the emission group */
double _phibar(gl_t *gl, long j, double phi){
  double phibar;
  if (j==nb) phibar=gl->model.chem.phiemit;
    else phibar=gl->model.chem.phimid[j]+gl->model.chem.Deltaphi[j]*_xi(gl,j,phi);
  return(phibar);
}


/* kinetic energy of group j=1..nb in eV at the potential phi, Eq. (7) */
double _eps_group(gl_t *gl, long j, double phi){
  return(max(0.0,phi-_phibar(gl,j,phi)));
}


/* speed of group j=1..nb in m/s at the potential phi, Eq. (8) */
static double _Ve_group(gl_t *gl, long j, double phi){
  return(sqrt(2.0*_eps_group(gl,j,phi)*echarge/emass));
}


/* TRUE if spec is one of the ballistic electron groups 1..nb (parent2026) */
bool _is_ballistic_group(long spec){
  return(spec>=0 && spec<ns && speciestype[spec]==SPECIES_ELECTRON && spec!=speceminus);
}


/* contravariant speed V_eb n.X of the ballistic group spec along theta, given the metrics X[dim]=X_theta,dim
   at the potential phi; the group travels at the speed of Eq. (8) of parent2026 along the unit vector n
   set by beamdir_x, beamdir_y and beamdir_z in the chem action */
static double _Vstar_ballistic_group(gl_t *gl, long spec, double phi, double *X){
  EXM_vec3D_t n;
  double nstar;
  long dim;
  find_beamdir(gl, n);
  nstar=0.0;
  for (dim=0; dim<nd; dim++) nstar+=n[dim]*X[dim];
  return(_Ve_group(gl,spec,phi)*nstar);
}


/* contravariant current of the ballistic groups (beam convection only) between node l and node _al(gl,l,theta,+1),
   upwinded as in find_Dstar_interface(); added to the total current in _Jtstar_interface() so that the Ohm-based
   potential equation sees the beams; it is not part of _Jstar_interface(), so E.J does not include the work done
   on the groups, whose energy reaches the total energy through Qbeam and the thresholds instead */
double _Jstar_ballistic_groups_interface(np_t *np, gl_t *gl, metrics_t metrics, long l, long theta){
  long k,lR,lup;
  double Vbstar,Jstar;

  lR=_al(gl,l,theta,+1);
  Jstar=0.0;
  if (is_node_valid(np[l],TYPELEVEL_FLUID) && is_node_valid(np[lR],TYPELEVEL_FLUID)){
    for (k=0; k<ncs; k++){
      if (_is_ballistic_group(k)) {
        lup=l;
        Vbstar=_Vstar_ballistic_group(gl,k,_phi(np[l],gl),metrics.X2[theta]);
        if (Vbstar<0.0) {
          lup=lR;
          Vbstar=_Vstar_ballistic_group(gl,k,_phi(np[lR],gl),metrics.X2[theta]);
        }
        Jstar+=_C(k)*_Nk(np[lup],gl,k)*Vbstar;
      }
    }
  }
  return(Jstar);
}


/* divergence of the current of the ballistic groups at node l, in A/m3, discretized as in the Ohm-based
   potential equation (interface currents of _Jtstar_interface(), zero at the dielectric boundaries) */
static double _divJb(np_t *np, gl_t *gl, long l){
  long theta,lL;
  double divJb;
  metrics_t metrics;

  divJb=0.0;
  for (theta=0; theta<nd; theta++){
    lL=_al(gl,l,theta,-1);
    find_metrics_at_interface(np,gl,l,_al(gl,l,theta,+1),theta,&metrics);
    if (!is_node_bdry_dielectric(np[l]) && !is_node_bdry_dielectric(np[_al(gl,l,theta,+1)]))
      divJb+=_Jstar_ballistic_groups_interface(np,gl,metrics,l,theta)*metrics.Omega;
    find_metrics_at_interface(np,gl,lL,l,theta,&metrics);
    if (!is_node_bdry_dielectric(np[lL]) && !is_node_bdry_dielectric(np[l]))
      divJb-=_Jstar_ballistic_groups_interface(np,gl,metrics,lL,theta)*metrics.Omega;
  }
  return(divJb/_Omega(np[l],gl));
}


/* reads the parameters of the fluid module specific to the ballistic electron model and checks the others; called once
   the other parameters of the fluid module are read */
void read_model_fluid_actions_ballistic(SOAP_codex_t *codex, gl_t *gl){
  find_double_var_from_codex(codex,"chibmin",&gl->model.fluid.chibmin);
  if (gl->model.fluid.INCLUDE_E_EFFECTS_ON_ELECTRON_MOBILITY)
    SOAP_fatal_error(codex,"INCLUDE_E_EFFECTS_ON_ELECTRON_MOBILITY must be set to NO with the ballistic electron model: the mobility of the bulk electrons is a function of Te, and the ballistic groups have none.");
}


/* beam convection of the ballistic groups, added to the matrix Dstar of find_Dstarmat() */
static void add_Dstarmat_ballistic_groups(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstar){
  double X[nd];
  long dim,k;
  for (dim=0; dim<nd; dim++) X[dim]=_X(np[l],theta,dim);
  for (k=0; k<ncs; k++){
    if (_is_ballistic_group(k)) Dstar[k][k]+=_Vstar_ballistic_group(gl,k,_phi(np[l],gl),X);
  }
}


/* at a boundary node, only the beam convection of the ballistic groups is kept in Dstar: it carries the emission
   group from the cathode into the domain */
static void set_Dstarmat_to_zero_except_ballistic_groups(sqmat_t Dstar){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      if (row!=col || !_is_ballistic_group(row)) Dstar[row][col]=0.0;
    }
  }
}


/* beam convection of the ballistic groups at the interface between lL and lR, with the speed of the upwind node,
   added to the vector Dstar of find_Dstar_interface() */
static void add_Dstar_interface_ballistic_groups(np_t *np, gl_t *gl, long lL, long lR, metrics_t metrics, long theta, flux_t Dstar){
  double Vbstar;
  long k;
  for (k=0; k<ncs; k++){
    if (_is_ballistic_group(k)) {
      Vbstar=_Vstar_ballistic_group(gl,k,_phi(np[lL],gl),metrics.X2[theta]);
      if (Vbstar<0.0) Vbstar=_Vstar_ballistic_group(gl,k,_phi(np[lR],gl),metrics.X2[theta]);
      Dstar[k]+=Vbstar;
    }
  }
}


/* beam convection of the ballistic groups at node l, added to the vector Dstar of find_Dstar() */
static void add_Dstar_ballistic_groups(np_t *np, gl_t *gl, long l, metrics_t metricsint, long theta, flux_t Dstar){
  long k;
  for (k=0; k<ncs; k++){
    if (_is_ballistic_group(k)) Dstar[k]+=_Vstar_ballistic_group(gl,k,_phi(np[l],gl),metricsint.X2[theta]);
  }
}

#endif



/* functions formerly in the fluid_plasma.c of _Navier-Stokes_plasma and _Favre-Reynolds_plasma */

#if defined(_FLUID_PLASMA) && defined(_FLUID_NEUTRALSTRANSPORT)

double _Te_local(np_t *np, gl_t *gl, long l, double Temin){
  double Te,EoverN,Emag;
  spec_t rhok;

  Emag=_Eemag_smoothed(np,gl,l);
  EoverN=Emag/_Nn(np[l],gl);
  find_rhok(np[l],rhok);
  Te=_Te_from_rhok_EoverN(rhok, EoverN);
  Te=max(Te,Temin);

  return(Te);
}


double _Tk(np_t *np, gl_t *gl, long l, long spec){
  double Tk;
  if (spec==speceminus) Tk=_Te(np[l]); else Tk=_T(np[l],gl);
  return(Tk);
}


double _Pkmem(np_t *np, gl_t *gl, long l, long spec){
  double Nk,Tk;    
  Nk=_Nk(np[l],gl,spec);
  Tk=_Tk(np,gl,l,spec);
  return(kB*Tk*Nk); 
}


double _Pk(np_t *np, gl_t *gl, long l, long spec){
  double Pk;
  if (_FLUIDPRIMMEM(np[l]) && spec<ncs){
    Pk=np[l].wk->Pkmem[spec];
  } else {
    Pk=_Pkmem(np,gl,l,spec);
  }
  return(Pk);
}


double _Nnmem(np_t np, gl_t *gl){
  double Nn;
  long spec; 
  Nn=0.0;
  for (spec=0; spec<ns; spec++){
    if (speciestype[spec]==SPECIES_NEUTRAL) { 
      Nn+=np.bs->U[spec]/_calM(spec)*calA;
    }
  }
  return(Nn);
}


double _Nn(np_t np, gl_t *gl){
  double Nn;
  if (_FLUIDPRIMMEM(np)){
    Nn=np.wk->Nnmem;
  } else {
    Nn=_Nnmem(np,gl);
  }
  return(Nn);
}


double _Nmem(np_t np, gl_t *gl){
  double N;
  long spec;
  
  N=0.0;
  for (spec=0; spec<ns; spec++){
    N+=np.bs->U[spec]/_calM(spec)*calA;
  }
  return(N);
}


double _N(np_t np, gl_t *gl){
  double N;
  if (_FLUIDPRIMMEM(np)){
    N=np.wk->Nmem;
  } else {
    N=_Nmem(np,gl);
  }
  return(N);
}


void find_dQbeam_dU(np_t np, gl_t *gl, flux_t dQbeamdU){
  long flux,spec;
  double dQbeamdN;
  flux_t dNdU;
  dQbeamdN=_dQbeam_dN(np,gl);
  for (flux=0; flux<nf; flux++) dNdU[flux]=0.0;
  for (spec=0; spec<ns; spec++) dNdU[spec]=1.0/_calM(spec)*calA;
  for (flux=0; flux<nf; flux++) dQbeamdU[flux]=dNdU[flux]*dQbeamdN;
}


static double _et_denom(double et, gl_t *gl){
  double ret;
  if (fabs(et)<gl->model.fluid.etdiv) {
    ret=gl->model.fluid.etdiv; 
  } else {
    ret=et;
  }
  return(ret);
}


/* TRUE if TEMODEL solves a transport equation for the energy of the bulk electrons */
static bool is_TEMODEL_transport(gl_t *gl){
  long TEMODEL;
  bool TRANSPORT;
  TEMODEL=gl->model.fluid.TEMODEL;
  TRANSPORT=(TEMODEL!=TEMODEL_CONSTANT && TEMODEL!=TEMODEL_TEQUILIBRIUM && TEMODEL!=TEMODEL_LOCAL);
#ifdef TEMODEL_TVEQUILIBRIUM
  if (TEMODEL==TEMODEL_TVEQUILIBRIUM) TRANSPORT=FALSE;
#endif
#ifdef TEMODEL_TRANSPORTEMF
  if (TEMODEL==TEMODEL_TRANSPORTEMF) TRANSPORT=FALSE;
#endif
  return(TRANSPORT);
}


void find_Dstar(np_t *np, gl_t *gl, long l, long theta, metrics_t metricsint, flux_t Dstar) {
  double Vstar;
  double DeltaVrstar,sigma;
  long flux,r,row,k,i,m;
  long dim,spec;
  dim_t dTedx,V,dTkdx;
  EXM_vec3D_t Vk,Ve,E,DeltaVr;
  double T,Te,etstar,Vkstar,Vestar,hk,Estar;
  
  //find_metrics_at_node(np, gl, l, theta, &metricsint);
  for (row=0; row<nf; row++){
    Dstar[row]=0.0;
  }  
  for (dim=0; dim<nd; dim++){
    V[dim]=_V(np[l],dim);
  }
  find_E(np[l], gl, E);
  Estar=0.0;
  Vstar=0.0;
  for (dim=0; dim<nd; dim++) {
    Estar+=E[dim]*metricsint.X2[theta][dim];
    Vstar+=V[dim]*metricsint.X2[theta][dim];
  }
  T=_T(np[l],gl);
  Te=_Te(np[l]);
  etstar=_etstar(np[l]);
  sigma=_sigma_positive(np,gl,l);

  for (r=0; r<ns; r++){
    find_DeltaVk(np, gl, l, r, DeltaVr);
    DeltaVrstar=0.0;
    for (dim=0; dim<nd; dim++) {
      DeltaVrstar+=DeltaVr[dim]*metricsint.X2[theta][dim];
    }
    Dstar[r]+=_betac(r)*DeltaVrstar;
    
    for (k=0; k<ns; k++){
      Dstar[k]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_betac(k)*_betac(r)*(DeltaVrstar+Vstar)/sigma*_rhok(np[l],r);
    }
  }

  for (k=0; k<ncs; k++){
    Dstar[k]+=(_s(k)-_betag(gl,k)+_betaa(gl,k))*_mu(np,gl,l,k)*Estar;
  }
#ifdef _CHEM_BALLISTIC
  add_Dstar_ballistic_groups(np, gl, l, metricsint, theta, Dstar);
#endif
  
  for (spec=0; spec<ns; spec++){
#ifdef _CHEM_BALLISTIC
    /* the ballistic groups carry directed kinetic energy, not enthalpy: they are left out of the enthalpy flux */
    if (_is_ballistic_group(spec)) continue;
#endif
    if (_Charge_number(spec)!=0){
      //find_Vk(np, gl, l, spec, Vk);
      if (spec==speceminus) {
        find_Ve_from_J(np, gl, l, Vk);
        //find_Ve_from_Je(np, gl, l, Vk);
      } else {
        find_Vk_from_Vk_at_interfaces(np, gl, l, spec, Vk);
      }
      
      for (i=0; i<nd; i++){
        dTkdx[i]=0.0;
        for (m=0; m<nd; m++) {
          dTkdx[i]+=_dTkdXj(np,gl,l,m,spec)*_X(np[l],m,i)/2.0;
        }
      }
      // get rid of the charged species thermal diffusion here, which we will add in Kstar matrix for stability
      //for (dim=0; dim<nd; dim++) Vk[dim]+=_mu(np,gl,l,spec)*kB/fabs(_C(spec))*dTkdx[dim];
    
      Vkstar=0.0e0;
      for (dim=0; dim<nd; dim++){
        Vkstar+=metricsint.X2[theta][dim]*Vk[dim];
      }
      if (spec==speceminus) hk=_hk_from_T_equilibrium(speceminus,Te); 
        else hk=_hk_from_T_equilibrium(spec,T);
      
      Dstar[fluxet]+=1.0/_et_denom(etstar,gl)*_betac(spec)*_w(np[l],spec)*hk*(Vkstar-Vstar);
    }
  }  

  if (is_TEMODEL_transport(gl)){
    //find_Vk(np, gl, l, speceminus, Ve);
    //find_Vk_from_Vk_at_interfaces(np, gl, l, speceminus, Ve);
    find_Ve_from_J(np, gl, l, Ve);
    //find_Ve_from_Je(np, gl, l, Ve);

    for (i=0; i<nd; i++){
      dTedx[i]=0.0;
      for (m=0; m<nd; m++) {
        dTedx[i]+=_dTedXj(np,gl,l,m)*_X(np[l],m,i)/2.0;
      }
    }
    // get rid of the electron thermal diffusion due to pressure gradients within electron velocity (we will add this in Kstar matrix for stability)
    //for (dim=0; dim<nd; dim++) Ve[dim]+=_mu(np,gl,l,speceminus)*kB/fabs(_C(speceminus))*dTedx[dim];


    Vestar=0.0e0;
    for (dim=0; dim<nd; dim++){
      Vestar+=metricsint.X2[theta][dim]*Ve[dim];
    }
    Dstar[fluxee]=5.0/3.0*(Vestar-Vstar); // check this later
  }

  if (is_node_bdry(np[l],TYPELEVEL_FLUID)){
#ifdef _CHEM_BALLISTIC
    /* at a boundary node, only the beam convection of the ballistic groups is kept: it carries the emission
       group from the cathode into the domain */
    for (flux=0; flux<nf; flux++) if (!_is_ballistic_group(flux)) Dstar[flux]=0.0;
    for (flux=0; flux<nf; flux++) Dstar[flux]*=metricsint.Omega/_Omega(np[l],gl);    
#else
    for (flux=0; flux<nf; flux++) Dstar[flux]=0.0;
#endif
  } else {
    for (flux=0; flux<nf; flux++) Dstar[flux]*=metricsint.Omega/_Omega(np[l],gl);    
  }
}


void find_Dstarplus(np_t *np, gl_t *gl, long l, long theta, metrics_t metricsint, flux_t Dstarplus) {
  flux_t Dstar;
  long dim, flux;
  double Xsum,astar,Vestar;

  //find_metrics_at_node(np, gl, l, theta, &metricsint);

  /*
  dim_t Ve;
  find_Ve_from_J(np, gl, l, Ve);
  Vestar=0.0e0;
  for (dim=0; dim<nd; dim++){
    Vestar+=metricsint.X2[theta][dim]*Ve[dim];
  }*/
  Vestar=0.0;

  /* find astar (for eigenvalue conditioning) */
  Xsum=0.0;
  for (dim=0; dim<nd; dim++) Xsum+=sqr(metricsint.X2[theta][dim]);
  Xsum=sqrt(Xsum);
  astar=_a(np[l],gl)*Xsum;

  /* first set Dstarplus to Dstar */
  find_Dstar(np, gl, l, theta, metricsint, Dstar);
  for (flux=0; flux<nf; flux++){
    Dstarplus[flux]=Dstar[flux];
  } 
 
  /* make all diagonal elements positive and condition the eigenvalues of the charged species*/
  for (flux=0; flux<nf; flux++){
    Dstarplus[flux]=0.5*(Dstar[flux]+fabs(Dstar[flux]));
    if (flux<ncs) Dstarplus[flux]=0.5*(Dstar[flux]+max(astar*gl->model.fluid.zetaD,fabs(Dstar[flux])));
    if (flux==fluxee) Dstarplus[flux]=0.5*(Dstar[flux]+max((astar+fabs(Vestar))*gl->model.fluid.zetaD,fabs(Dstar[flux])));
  }
  
}


void find_Dstarminus(np_t *np, gl_t *gl, long l, long theta, metrics_t metricsint, flux_t Dstarminus) {
  flux_t Dstar;
  long flux, dim;
  double Xsum,astar,Vestar;

  //find_metrics_at_node(np, gl, l, theta, &metricsint);

  /*
  dim_t Ve;
  find_Ve_from_J(np, gl, l, Ve);
  Vestar=0.0e0;
  for (dim=0; dim<nd; dim++){
    Vestar+=metricsint.X2[theta][dim]*Ve[dim];
  }*/
  Vestar=0.0;

 
  /* find astar (for eigenvalue conditioning) */
  Xsum=0.0;
  for (dim=0; dim<nd; dim++) Xsum+=sqr(metricsint.X2[theta][dim]);
  Xsum=sqrt(Xsum);
  astar=_a(np[l],gl)*Xsum;


  /* first set Dstarplus to Dstar */
  find_Dstar(np, gl, l, theta, metricsint, Dstar);
  for (flux=0; flux<nf; flux++){
    Dstarminus[flux]=Dstar[flux];
  } 

  /* make all diagonal elements negative */
  for (flux=0; flux<nf; flux++){
    Dstarminus[flux]=0.5*(Dstar[flux]-fabs(Dstar[flux]));
    if (flux<ncs) Dstarminus[flux]=0.5*(Dstar[flux]-max(astar*gl->model.fluid.zetaD,fabs(Dstar[flux])));
    if (flux==fluxee) Dstarminus[flux]=0.5*(Dstar[flux]-max((fabs(Vestar)+astar)*gl->model.fluid.zetaD,fabs(Dstar[flux])));
  }
}


double _DeltaVkstar(np_t *np, gl_t *gl, long l, long k, long theta){
  double ret;
  long dim;
  EXM_vec3D_t DeltaVk;
  find_DeltaVk(np, gl, l, k, DeltaVk);

  ret=0.0e0;
  for (dim=0; dim<nd; dim++){
    ret=ret+_X(np[l],theta,dim)*DeltaVk[dim];
  }
  return(ret);
}


double _DeltaVkstar_interface(np_t *np, gl_t *gl, long lL, long lR, long k, long theta, metrics_t metrics){
  double ret;
  long dim;
  EXM_vec3D_t DeltaVkL,DeltaVkR;
  find_DeltaVk(np, gl, lL, k, DeltaVkL);
  find_DeltaVk(np, gl, lR, k, DeltaVkR);

  ret=0.0e0;
  for (dim=0; dim<nd; dim++){
    ret=ret+metrics.X2[theta][dim]*minmod(DeltaVkL[dim],DeltaVkR[dim]);
  }
  return(ret);
}


void find_Dstarmat(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstar) {
  double Vstar;
  double DeltaVrstar,sigma;
  long r,row,col,k,i,m;
  long dim,spec;
  dim_t dTedx,dTdx;
  EXM_vec3D_t Vk,Ve,E;
  double T,Te,etstar,Vkstar,Vestar,hk,Estar;
  
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Dstar[row][col]=0.0;
    }
  }  

  Vstar=_Vstar(np[l],theta);
  find_E(np[l], gl, E);
  Estar=0.0;
  for (dim=0; dim<nd; dim++) Estar+=E[dim]*_X(np[l],theta,dim);
  T=_T(np[l],gl);
  Te=_Te(np[l]);
  etstar=_etstar(np[l]);
  sigma=_sigma_positive(np,gl,l);

  for (r=0; r<ns; r++){
    DeltaVrstar=_DeltaVkstar(np,gl,l,r,theta);
    Dstar[r][r]+=_betac(r)*DeltaVrstar;
    
    for (k=0; k<ns; k++){
      Dstar[k][k]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*_mu(np,gl,l,k)*_betac(k)*_betac(r)*(DeltaVrstar+Vstar)/sigma*_rhok(np[l],r);
    }
  }

  for (k=0; k<ncs; k++){
    Dstar[k][k]+=(_s(k)-_betag(gl,k)+_betaa(gl,k))*_mu(np,gl,l,k)*Estar;
  }
#ifdef _CHEM_BALLISTIC
  add_Dstarmat_ballistic_groups(np, gl, l, theta, Dstar);
#endif
  
  for (i=0; i<nd; i++){
    dTedx[i]=0.0;
    for (m=0; m<nd; m++) {
      dTedx[i]+=_dTedXj(np,gl,l,m)*_X(np[l],m,i)/2.0;
    }
  }

  for (i=0; i<nd; i++){
    dTdx[i]=0.0;
    for (m=0; m<nd; m++) {
      dTdx[i]+=_dTkdXj(np,gl,l,m,ns-1)*_X(np[l],m,i)/2.0;
    }
  }
  
  
  for (spec=0; spec<ns; spec++){
#ifdef _CHEM_BALLISTIC
    /* the ballistic groups carry directed kinetic energy, not enthalpy: they are left out of the enthalpy flux */
    if (_is_ballistic_group(spec)) continue;
#endif
    if (_Charge_number(spec)!=0){
      //find_Vk(np, gl, l, spec, Vk);
      if (spec==speceminus) {
        find_Ve_from_J(np, gl, l, Vk);
        //find_Ve_from_Je(np, gl, l, Vk);
      } else {
        find_Vk_from_Vk_at_interfaces(np, gl, l, spec, Vk);
      }
      // get rid of the charged species thermal diffusion here, which we will add in Kstar matrix for stability
      for (dim=0; dim<nd; dim++) {
        if (spec==speceminus) Vk[dim]+=_mu(np,gl,l,spec)*kB/fabs(_C(spec))*dTedx[dim];
          else Vk[dim]+=_mu(np,gl,l,spec)*kB/fabs(_C(spec))*dTdx[dim];
      }
      Vkstar=0.0e0;
      for (dim=0; dim<nd; dim++){
        Vkstar+=_X(np[l],theta,dim)*Vk[dim];
      }
      if (spec==speceminus) hk=_hk_from_T_equilibrium(speceminus,Te); 
        else hk=_hk_from_T_equilibrium(spec,T);
      
      Dstar[fluxet][fluxet]+=1.0/_et_denom(etstar,gl)*_betac(spec)*_w(np[l],spec)*hk*(Vkstar-Vstar);


        
    }
  }  

  if (is_TEMODEL_transport(gl)){
    //find_Vk(np, gl, l, speceminus, Ve);
    //find_Vk_from_Vk_at_interfaces(np, gl, l, speceminus, Ve);
    find_Ve_from_J(np, gl, l, Ve);
    //find_Ve_from_Je(np, gl, l, Ve);

    // get rid of the electron thermal diffusion here, which we will add in Kstar matrix for stability
    for (dim=0; dim<nd; dim++) Ve[dim]+=_mu(np,gl,l,speceminus)*kB/fabs(_C(speceminus))*dTedx[dim];


    Vestar=0.0e0;
    for (dim=0; dim<nd; dim++){
      Vestar+=_X(np[l],theta,dim)*Ve[dim];
    }
    Dstar[fluxee][fluxee]=5.0/3.0*(Vestar-Vstar); // check this later
  }

  if (is_node_bdry(np[l],TYPELEVEL_FLUID)){
#ifdef _CHEM_BALLISTIC
    set_Dstarmat_to_zero_except_ballistic_groups(Dstar);
#else
    set_matrix_to_zero(Dstar);
#endif
  }
}


void find_Dstar_interface(np_t *np, gl_t *gl, long lL, long theta, flux_t Dstar) {
  double Vstar;
  double DeltaVrstar,sigma;
  long r,row,k,lR,i,m;
  long dim,spec,flux;
  EXM_vec3D_t Vk,Ve,VeR,VeL;
  dim_t dTedx;
  double T,Te,etstar,Vkstar,Vestar,hk,Estar;
  metrics_t metrics;
  
  for (row=0; row<nf; row++){
    Dstar[row]=0.0;
  }  

  lR=_al(gl,lL,theta,+1);
  find_metrics_at_interface(np, gl, lL, lR, theta, &metrics);
  Vstar=0.5*(_Vstar(np[lL],theta)+_Vstar(np[lR],theta));
  find_Estar_at_interface(np, gl, lL, lR, theta, &Estar);
  T=0.5*(_T(np[lL],gl)+_T(np[lR],gl));
  Te=0.5*(_Te(np[lL])+_Te(np[lR]));
  etstar=0.5*(_etstar(np[lL])+_etstar(np[lR]));
  sigma=0.5*(_sigma_positive(np,gl,lL)+_sigma_positive(np,gl,lR));

  for (r=0; r<ns; r++){
    DeltaVrstar=_DeltaVkstar_interface(np, gl, lL, lR, r, theta, metrics);
    Dstar[r]+=_betac(r)*DeltaVrstar;
    
    for (k=0; k<ns; k++){
      Dstar[k]+=1.0/_m(r)*_betaa(gl,k)*_C(r)*0.5*(_mu(np,gl,lL,k)+_mu(np,gl,lR,k))*_betac(k)*_betac(r)*(DeltaVrstar+Vstar)/sigma*0.5*(_rhok(np[lL],r)+_rhok(np[lR],r));
    }
  }

  for (k=0; k<ncs; k++){
    Dstar[k]+=(_s(k)-_betag(gl,k)+_betaa(gl,k))*0.5*(_mu(np,gl,lL,k)+_mu(np,gl,lR,k))*Estar;
  }
#ifdef _CHEM_BALLISTIC
  add_Dstar_interface_ballistic_groups(np, gl, lL, lR, metrics, theta, Dstar);
#endif
  
  for (spec=0; spec<ns; spec++){
#ifdef _CHEM_BALLISTIC
    /* the ballistic groups carry directed kinetic energy, not enthalpy: they are left out of the enthalpy flux */
    if (_is_ballistic_group(spec)) continue;
#endif
    if (_Charge_number(spec)!=0){
      find_Vk_at_interface(np, gl, lL, theta, spec, Vk);
      Vkstar=0.0e0;
      for (dim=0; dim<nd; dim++){
        Vkstar+=metrics.X2[theta][dim]*Vk[dim];
      }
      if (spec==speceminus) hk=_hk_from_T_equilibrium(speceminus,Te); 
        else hk=_hk_from_T_equilibrium(spec,T);
      Dstar[fluxet]+=1.0/_et_denom(etstar,gl)*_betac(spec)*0.5*(_w(np[lL],spec)+_w(np[lR],spec))*hk*(Vkstar-Vstar);
    }
  }  

  if (is_TEMODEL_transport(gl)){
    //find_Vk_at_interface(np, gl, lL, theta, speceminus, Ve);
    find_Ve_from_J(np, gl, lL, VeL);
    find_Ve_from_J(np, gl, lR, VeR);
    for (dim=0; dim<nd; dim++) Ve[dim]=minmod(VeL[dim],VeR[dim]); 
    for (i=0; i<nd; i++){
      dTedx[i]=0.0;
      for (m=0; m<nd; m++) {
        dTedx[i]+=_dTedXj_interface_upwind(np, gl, metrics, lL, theta, m)*metrics.X2[m][i];
      }
    }
    // get rid of the electron thermal diffusion here, which we will add in Kstar matrix for stability
    for (dim=0; dim<nd; dim++) Ve[dim]+=min(_mu(np,gl,lL,speceminus),_mu(np,gl,lR,speceminus))*kB/fabs(_C(speceminus))*dTedx[dim];

    Vestar=0.0e0;
    for (dim=0; dim<nd; dim++){
      Vestar+=metrics.X2[theta][dim]*Ve[dim];
    }
    Dstar[fluxee]=5.0/3.0*(Vestar-Vstar); // check this later
  }
  for (flux=0; flux<nf; flux++) Dstar[flux]*=metrics.Omega;

  if (is_node_bdry(np[lL],TYPELEVEL_FLUID) || is_node_bdry(np[lR],TYPELEVEL_FLUID)){
   // for (flux=0; flux<nf; flux++) Dstar[flux]=0.0;
  }
}



void find_Dstarmatplus(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstarplus) {
  sqmat_t Dstar;
  long row, col, dim;
  double Xsum,astar,Vestar;


  /* find astar (for eigenvalue conditioning) */
  Xsum=0.0;
  for (dim=0; dim<nd; dim++) Xsum+=sqr(_X(np[l],theta,dim));
  Xsum=sqrt(Xsum);
  astar=_a(np[l],gl)*Xsum;
/*  find_Ve_from_J(np, gl, l, Ve);
  Vestar=0.0e0;
  for (dim=0; dim<nd; dim++){
    Vestar+=_X(np[l],theta,dim)*Ve[dim];
  }
*/
  Vestar=0.0;



  /* first set Dstarplus to Dstar */
  find_Dstarmat(np, gl, l, theta, Dstar);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Dstarplus[row][col]=Dstar[row][col];
    }
  } 
 
  /* make all diagonal elements positive and condition the eigenvalues of the charged species*/
  for (row=0; row<nf; row++){
    Dstarplus[row][row]=0.5*(Dstar[row][row]+fabs(Dstar[row][row]));
    if (row<ncs) Dstarplus[row][row]=0.5*(Dstar[row][row]+max(astar*gl->model.fluid.zetaD,fabs(Dstar[row][row])));
    if (row==fluxee) Dstarplus[row][row]=0.5*(Dstar[row][row]+max((astar+fabs(Vestar))*gl->model.fluid.zetaD,fabs(Dstar[row][row])));
    
  }
  
  /* fix negative species rows 
     Notes: - this would yield exactly the Steger-Warming scheme if only one type of negative species is present
            - when more than one type of negative species are present, this will still 
              result in Dstarplus+Dstarminus=Dstar*/
  for (row=0; row<ncs; row++){
    for (col=0; col<ncs; col++){
      if (col!=row)
        Dstarplus[row][col]=Dstar[row][col]/notzero(Dstar[col][col]-Dstar[row][row],1e-99)
                           *(Dstarplus[col][col]-Dstarplus[row][row]);
    }
  }


}


void find_Dstarmatminus(np_t *np, gl_t *gl, long l, long theta, sqmat_t Dstarminus) {
  sqmat_t Dstar;
  long row, col, dim;
  double Xsum,astar,Vestar;
  /*
  find_Ve_from_J(np, gl, l, Ve);
  Vestar=0.0e0;
  for (dim=0; dim<nd; dim++){
    Vestar+=_X(np[l],theta,dim)*Ve[dim];
  }*/
  Vestar=0.0;

 
  /* find astar (for eigenvalue conditioning) */
  Xsum=0.0;
  for (dim=0; dim<nd; dim++) Xsum+=sqr(_X(np[l],theta,dim));
  Xsum=sqrt(Xsum);
  astar=_a(np[l],gl)*Xsum;


  /* first set Dstarplus to Dstar */
  find_Dstarmat(np, gl, l, theta, Dstar);
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      Dstarminus[row][col]=Dstar[row][col];
    }
  } 

  /* make all diagonal elements negative */
  for (row=0; row<nf; row++){
    Dstarminus[row][row]=0.5*(Dstar[row][row]-fabs(Dstar[row][row]));
    if (row<ncs) Dstarminus[row][row]=0.5*(Dstar[row][row]-max(astar*gl->model.fluid.zetaD,fabs(Dstar[row][row])));
    if (row==fluxee) Dstarminus[row][row]=0.5*(Dstar[row][row]-max((fabs(Vestar)+astar)*gl->model.fluid.zetaD,fabs(Dstar[row][row])));
  }
  
  /* fix negative species rows 
     Notes: - this would yield exactly the Steger-Warming scheme if only one type of negative species is present
            - when more than one type of negative species are present, this will still 
              result in Dstarplus+Dstarminus=Dstar*/
  for (row=0; row<ncs; row++){
    for (col=0; col<ncs; col++){
      if (col!=row)
        Dstarminus[row][col]=Dstar[row][col]/notzero(Dstar[col][col]-Dstar[row][row],1e-99)
                           *(Dstarminus[col][col]-Dstarminus[row][row]);
    }
  }


}


void find_H1(np_t *np, gl_t *gl, long l, flux_t H){
  long spec,flux;
  
  for (flux=0; flux<nf; flux++) {
    H[flux]=0.0;
  }

  for (spec=0; spec<ncs; spec++){
    
    H[spec]=_rhok(np[l],spec)*_mu(np,gl,l,spec);
  }
}


void find_H2(np_t *np, gl_t *gl, long l, flux_t H){
  double sigma_positive,mu;
  long spec,flux;
  for (flux=0; flux<nf; flux++) {
    H[flux]=0.0;
  }
  sigma_positive=_sigma_positive(np, gl, l);
  for (spec=0; spec<ncs; spec++){
    mu=_mu(np,gl,l,spec);
    H[spec]=_rhok(np[l],spec)*mu/sigma_positive;
  }
}


void find_dH1_dUstar(np_t *np, gl_t *gl, long l, sqmat_t dHdUstar){
  long row,col,spec;
  double Omega;
  spec_t mu;
#ifdef TEST
  spec_t rhok,dmukdrhok;
  flux_t dTkdU;
  long flux;
  double dmukdTk;
#endif

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dHdUstar[row][col]=0.0;
    }
  }	
  for (spec=0; spec<ns; spec++) {
    mu[spec]=_mu(np,gl,l,spec);
  }

  Omega=_Omega(np[l],gl);

//    H[spec]=rhok[spec]*mu[spec];

  /* wrt to rhok */
  for (spec=0; spec<ns; spec++){
      dHdUstar[spec][spec]+=1.0/Omega*mu[spec];
  }



#ifdef TEST
  /* wrt to muk */
  for (spec=0; spec<ns; spec++) {
    rhok[spec]=_rhok(np[l],spec);
  }
  for (spec=0; spec<ncs; spec++){
      if (spec==speceminus) find_dTe_dU(np[l], gl,  dTkdU); else find_dT_dU(np[l], gl, dTkdU);
      find_dmuk_from_rhok_Tk_Ek(gl, rhok, _Tk(np,gl,l,spec), _Ekmag(np,gl,l,spec), spec, &dmukdTk, dmukdrhok);
      for (flux=0; flux<nf; flux++) dHdUstar[spec][flux]+=1.0/Omega*rhok[spec]*(dmukdTk*dTkdU[flux]);
      for (flux=0; flux<ns; flux++) dHdUstar[spec][flux]+=1.0/Omega*rhok[spec]*(dmukdrhok[flux]);
  }
#endif
}


void find_dH2_dUstar(np_t *np, gl_t *gl, long l, sqmat_t dHdUstar){
  long row,col,spec;
  double Omega,sigma_positive;
  spec_t mu;
#ifdef TEST
  spec_t rhok,dmukdrhok;
  flux_t dTkdU,dsigmadU;
  flux_t dmukdU[ncs];
  long flux;
  double dmukdTk;
#endif

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      dHdUstar[row][col]=0.0;
    }
  }	
  for (spec=0; spec<ns; spec++) {
    mu[spec]=_mu(np,gl,l,spec);
  }

  Omega=_Omega(np[l],gl);
  sigma_positive=_sigma_positive(np,gl,l);

//  H[spec]=rhok[spec]*mu[spec]/sigma_positive;
//  H[spec]=rhok[spec]*mu[spec]/notzero(sigma,1e-99);


  /* wrt to rhok */
  for (spec=0; spec<ns; spec++){
      dHdUstar[spec][spec]+=1.0/Omega*mu[spec]/sigma_positive;
  }



#ifdef TEST
  /* wrt to muk and to sigma=sum_c |C_c| N_c mu_c */
  for (spec=0; spec<ns; spec++) {
    rhok[spec]=_rhok(np[l],spec);
  }
  for (spec=0; spec<ncs; spec++){
    for (flux=0; flux<nf; flux++) dmukdU[spec][flux]=0.0;
    if (spec==speceminus) find_dTe_dU(np[l], gl,  dTkdU); else find_dT_dU(np[l], gl, dTkdU);
    find_dmuk_from_rhok_Tk_Ek(gl, rhok, _Tk(np,gl,l,spec), _Ekmag(np,gl,l,spec), spec, &dmukdTk, dmukdrhok);
    for (flux=0; flux<nf; flux++) dmukdU[spec][flux]+=dmukdTk*dTkdU[flux];
    for (flux=0; flux<ns; flux++) dmukdU[spec][flux]+=dmukdrhok[flux];
  }
  for (flux=0; flux<nf; flux++) dsigmadU[flux]=0.0;
  if (_sigma(np,gl,l)>gl->model.fluid.sigmadiv){
    for (spec=0; spec<ncs; spec++){
      dsigmadU[spec]+=fabs(_C(spec))*mu[spec]/_m(spec);
      for (flux=0; flux<nf; flux++) dsigmadU[flux]+=fabs(_C(spec))*rhok[spec]/_m(spec)*dmukdU[spec][flux];
    }
  }
  for (spec=0; spec<ncs; spec++){
    for (flux=0; flux<nf; flux++) 
      dHdUstar[spec][flux]+=1.0/Omega*rhok[spec]*(dmukdU[spec][flux]/sigma_positive-mu[spec]*dsigmadU[flux]/sqr(sigma_positive));
  }
#endif
}

#endif



/* ballistic electron model (parent2026), continued: source terms, boundary conditions and post variables */

#if defined(_CHEM_BALLISTIC) && defined(_FLUID_PLASMA)

/* the ambipolar form (betaa) of the drift flux mu_k rho_k J/sigma drops (mu_k rho_k/sigma) div(J), taking
   div(J)=-d(rhoc)/dt; but J excludes the ballistic groups, so that div(J)=-d(rhoc)/dt-div(Jb): the missing
   part is added to the source terms S of find_Semfield() */
void add_Semfield_ballistic_groups(np_t *np, gl_t *gl, long l, flux_t S){
  double divJb,sigma;
  long spec;
  divJb=_divJb(np,gl,l);
  sigma=_sigma_positive(np,gl,l);
  for (spec=0; spec<ns; spec++){
    S[spec]-=_betaa(gl,spec)*_mu(np,gl,l,spec)*_rhok(np[l],spec)/sigma*divJb;
  }
}


/* adds one contribution to an element of the source Jacobian; as in the chem module,
   only the negative contributions are kept unless TEST is defined */
static inline void add_to_dS_term(double *dS, double term){
#ifdef TEST
  *dS+=term;
#else
  *dS+=min(0.0,term);
#endif
}


/* energy source terms of the reactions of the ballistic electron model (parent2026), added to the source terms S of find_Schem(): the bulk
   electrons receive Qborn and Qbeam and lose Qrun (Eqs. (41)-(43)); the kinetic energy of the ballistic groups is not
   part of the total energy, which hence receives what the groups lose in collisions (Qbeam and the thresholds of the
   excitations and ionizations) and loses what the bulk gives to the groups through runaway (Qrun) */
void add_Schem_ballistic(np_t *np, gl_t *gl, long l, spec_t rhok, double T, double Te, flux_t S){
  double Q[numQballistic],dQdT[numQballistic],dQdTe[numQballistic];
  spec_t dQdrhok[numQballistic];
  find_Q_ballistic(np[l], gl, rhok, T, Te, FALSE, Q, dQdrhok, dQdT, dQdTe);
  if (is_TEMODEL_transport(gl)) {
    S[fluxee]+=Q[QB_BORN]+Q[QB_BEAM]-Q[QB_RUN];
  }
  S[fluxet]+=Q[QB_BEAM]+Q[QB_THRESHOLD]-Q[QB_RUN];
}


/* Jacobian of the energy source terms of add_Schem_ballistic(), added to dS_dU of find_dSchem_dU(); Qborn, Qbeam
   and the thresholds are gains and Qrun is a loss: each contribution goes through add_to_dS_term() */
void add_dSchem_dU_ballistic(np_t *np, gl_t *gl, long l, spec_t rhok, double T, double Te, flux_t dTdU, flux_t dTedU, sqmat_t dS_dU){
  long s,cnt,row,n;
  double sgn;
  double Q[numQballistic],dQdT[numQballistic],dQdTe[numQballistic];
  spec_t dQdrhok[numQballistic];
  find_Q_ballistic(np[l], gl, rhok, T, Te, TRUE, Q, dQdrhok, dQdT, dQdTe);
  for (cnt=0; cnt<6; cnt++){
    switch (cnt){
      case 0: row=fluxee; n=QB_BORN;      sgn=+1.0; break;
      case 1: row=fluxee; n=QB_BEAM;      sgn=+1.0; break;
      case 2: row=fluxee; n=QB_RUN;       sgn=-1.0; break;
      case 3: row=fluxet; n=QB_BEAM;      sgn=+1.0; break;
      case 4: row=fluxet; n=QB_THRESHOLD; sgn=+1.0; break;
      default: row=fluxet; n=QB_RUN;      sgn=-1.0; break;
    }
    if (row==fluxee && !is_TEMODEL_transport(gl)) continue;
    for (s=0; s<ns; s++) add_to_dS_term(&(dS_dU[row][s]),sgn*dQdrhok[n][s]);
    for (s=0; s<nf; s++) add_to_dS_term(&(dS_dU[row][s]),sgn*dQdT[n]*dTdU[s]);
    for (s=0; s<nf; s++) add_to_dS_term(&(dS_dU[row][s]),sgn*dQdTe[n]*dTedU[s]);
  }
}


/* electron densities Nk at the cathode wall node lA, once those of the ions are set: with ballistic electrons, the
   secondary electrons of an electrode that the beam direction enters all go to the emission group nb, whose flux normal
   to the wall is gamma times the ion flux to the wall (Eq. (50) of parent2026); the bulk and the volumetric groups receive
   nothing from the cathode (Eq. (51)); those of a dielectric wall or of a wall that the beam direction does not enter go
   to the bulk, as without ballistic electrons; gamma is the secondary emission coefficient of the wall */
void update_Nk_electrons_at_cathode_ballistic(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn, double gamma, spec_t Nk){
  double Gammai,Vemitn,beamdirn,Newall,relaxNe_bdry;
  long dim,spec;
  dim_t n;
  EXM_vec3D_t beamdir;

  relaxNe_bdry=gl->cycle.fluid.relaxNe_bdry;
  beamdirn=0.0;
  if (gl->model.chem.BALLISTIC && !is_node_bdry_dielectric(np[lA])){
    if (!find_unit_vector_normal_to_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK, n))
      fatal_error("Problem finding the unit vector normal to the cathode in update_Nk_electrons_at_cathode_ballistic().");
    find_beamdir(gl, beamdir);
    for (dim=0; dim<nd; dim++) beamdirn+=beamdir[dim]*n[dim];
  }
  if (beamdirn<=0.0){
    /* without ballistic electrons, or at a dielectric wall or a wall that the beam direction does not enter, the
       secondary electrons go to the bulk, whose drift flux is then gamma times the ion drift flux, as in the
       drift-diffusion model; the groups receive nothing */
    Newall=0.0;
    for (spec=0; spec<ncs; spec++){
      if (speciestype[spec]==SPECIES_IONPLUS) Newall+=_Nk(np[lB],gl,spec)*_mu(np,gl,lB,spec);
    }
    Newall*=gamma/_mu(np,gl,lB,speceminus);
    for (spec=0; spec<ns; spec++){
      if (speciestype[spec]==SPECIES_ELECTRON){
        if (spec==speceminus) Nk[spec]=_Nk(np[lA],gl,spec)+relaxNe_bdry*(Newall-_Nk(np[lA],gl,spec));
          else Nk[spec]=_Nk(np[lA],gl,spec)+relaxNe_bdry*(0.0-_Nk(np[lA],gl,spec));
      }
    }
  } else {
    Gammai=0.0;
    for (spec=0; spec<ncs; spec++){
      if (speciestype[spec]==SPECIES_IONPLUS) 
        Gammai+=Nk[spec]*max(0.0,_Vk_drift_perp_to_bdry(np, gl, lA, lB, lC, theta, thetasgn, spec));
    }
    Vemitn=_Ve_group(gl,specegroup(nb),_phi(np[lA],gl))*beamdirn;
    if (Vemitn<=0.0) fatal_error("The emission group has no speed at the cathode: the potential there must be above phicathode-epsemit.");
    for (spec=0; spec<ns; spec++){
      if (speciestype[spec]==SPECIES_ELECTRON){
        if (spec==specegroup(nb)) Newall=gamma*Gammai/Vemitn; else Newall=0.0;
        Nk[spec]=_Nk(np[lA],gl,spec)+relaxNe_bdry*(Newall-_Nk(np[lA],gl,spec));
      }
    }
  }
}


/* bulk electron density Nk[speceminus] at the anode wall node lA of a wall that does not absorb the electrons: zero
   gradient (Eq. (52) of parent2026); the ballistic groups flow out, with Nk set to the near-wall value by the fluid
   module */
void update_Nk_electrons_at_anode_ballistic(np_t *np, gl_t *gl, long lA, long lB, spec_t Nk){
  Nk[speceminus]=_Nk(np[lA],gl,speceminus)+gl->cycle.fluid.relaxNe_bdry*(_Nk(np[lB],gl,speceminus)-_Nk(np[lA],gl,speceminus));
}


/* names of the numpostvarballistic post variables of the ballistic groups, varnum=0..numpostvarballistic-1: the kinetic
   energy in eV and the speed in m/s of each group, the total electron number density in 1/m3, the mean energy of all
   the electrons in eV, and the current density of the ballistic groups along beamdir in A/m2 */
void find_post_variable_name_ballistic(long varnum, char *varname){
  if (varnum>=0 && varnum<nb) sprintf(varname,"eps_e%ld",varnum+1);
  if (varnum>=nb && varnum<2*nb) sprintf(varname,"V_e%ld",varnum-nb+1);
  if (varnum==2*nb) sprintf(varname,"Ne_total");
  if (varnum==2*nb+1) sprintf(varname,"eps_mean");
  if (varnum==2*nb+2) sprintf(varname,"J_beam");
}


/* values of the post variables of find_post_variable_name_ballistic() at node l */
void find_post_variable_value_ballistic(np_t *np, long l, gl_t *gl, long varnum, double *varvalue){
  double phi,Ne,epssum,Jbeam;
  long spec;
  phi=_phi(np[l],gl);
  if (varnum>=0 && varnum<nb) *varvalue=_eps_group(gl,varnum+1,phi);
  if (varnum>=nb && varnum<2*nb) *varvalue=_Ve_group(gl,varnum-nb+1,phi);
  if (varnum>=2*nb){
    Ne=_Nk(np[l],gl,speceminus);
    epssum=_Nk(np[l],gl,speceminus)*1.5*kB*_Te(np[l])/echarge;
    Jbeam=0.0;
    for (spec=0; spec<ns; spec++){
      if (_is_ballistic_group(spec)){
        Ne+=_Nk(np[l],gl,spec);
        epssum+=_Nk(np[l],gl,spec)*_eps_group(gl,spec,phi);
        Jbeam+=_C(spec)*_Nk(np[l],gl,spec)*_Ve_group(gl,spec,phi);
      }
    }
    if (varnum==2*nb) *varvalue=Ne;
    if (varnum==2*nb+1) *varvalue=epssum/notzero(Ne,1e-99);
    if (varnum==2*nb+2) *varvalue=Jbeam;
  }
}

#endif
