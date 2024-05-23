// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2000-2002,2017,2020 Bernard Parent

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

#include "fluid.h"
#include "fluid_bdry.h"
#include "fluid_source.h"
#include <src/bdry.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
#include <model/emfield/_emfield.h>
#include <model/share/fluid_share.h>
#include <model/share/model_share.h>
#include <model/_model.h>



#define BDRY_INFLOWSUPERSONIC  0 
#define BDRY_INFLOWSUBSONIC1 7
#define BDRY_INFLOWSUBSONICMASSFLOWFIXED1 16
#define BDRY_OUTFLOWSUPERSONIC1 1
#define BDRY_OUTFLOWSUBSONIC1 4
#define BDRY_OUTFLOWSUBSONICMFIXED1 15
#define BDRY_SYMMETRICAL1 9
#define BDRY_SYMMETRICAL2 14
#define BDRY_WALLTFIXED1 3
#define BDRY_WALLADIABATIC1 6
#define BDRY_SLIPWALL1 12
#define BDRY_FREESTREAM1 2
#define BDRY_WALLTFIXEDINJECTION1 17
#define BDRY_WALLABLATION1 20
#define BDRY_INFLOWINJECTION1 18
#define BDRY_INFLOWINJECTIONCHOKED1 19


void write_bdry_fluid_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n",_FLUID_ACTIONNAME);
  wfprintf(*controlfile,
    "    {\n"
    "    _________________________________________________________________________________________\n"
    "\n" 
    "    Boundary Condition Type           ID  Description\n"
    "    _________________________________________________________________________________________\n"
    "\n" 
    "    BDRY_INFLOWSUPERSONIC             %c   Inflow, supersonic\n"
    "    BDRY_INFLOWSUBSONIC1              %c   Inflow, subsonic, Tstag, Pstag fixed\n"
    "    BDRY_INFLOWSUBSONICMASSFLOWFIXED1 %c   Inflow, subsonic, Pstag, Massflow/Area fixed\n"
    "    BDRY_INFLOWINJECTION1             %c   Inflow, param Tstag, Pstag, specCs \n"
    "    BDRY_INFLOWINJECTIONCHOKED1       %c   Inflow, param Tstag, Pstag, specCs \n"
    "    BDRY_OUTFLOWSUPERSONIC1           %c   Outflow, supersonic\n"
    "    BDRY_OUTFLOWSUBSONIC1             %c   Outflow, subsonic, P fixed, param P\n"
    "    BDRY_OUTFLOWSUBSONICMFIXED1       %c   Outflow, subsonic, M fixed, param M\n"
    "    BDRY_SYMMETRICAL2                 %c   Symmetrical, 2nd order\n"
    "    BDRY_SYMMETRICAL1                 %c   Symmetrical, 1st order\n"
    "    BDRY_WALLTFIXED1                  %c   Wall, T specified, param Twall\n"
    "    BDRY_WALLTFIXEDINJECTION1         %c   Wall, T specified, param Twall, specCs,mdotCs[kg/m2s], ...\n"
    "    BDRY_WALLABLATION1                %c   Wall, param specC,chiC,specCs,chiCs,...\n"
    "    BDRY_WALLADIABATIC1               %c   Wall, Adiabatic\n"
    "    BDRY_SLIPWALL1                    %c   Slip wall, 1st order, Adiabatic\n"
    "    BDRY_FREESTREAM1                  %c   Freestream, 1o, params Vx,Vy,"if3DL("Vz,")" P, T\n"
    "    _________________________________________________________________________________________\n"
    "    }\n"
    "    Twall=300.0; {K}\n"
    "    All(BDRY_WALLTFIXED1,Twall);\n"
    "    Plane(\"i\",is,BDRY_INFLOWSUPERSONIC);\n"
    "    Plane(\"i\",ie,BDRY_OUTFLOWSUPERSONIC1);\n"
    "    Plane(\"j\",js,BDRY_WALLTFIXED1,Twall);\n"
    "    Plane(\"j\",je,BDRY_WALLTFIXED1,Twall);\n"
#ifdef _3D
    "    Plane(\"k\",ks,BDRY_SYMMETRICAL2);\n"
    "    Plane(\"k\",ke,BDRY_SYMMETRICAL2);\n"
#endif
    "    {\n"
    "    Region(is" if2DL(",js") if3DL(",ks") ",  ie" if2DL(",je") if3DL(",ke") ",  BDRY_INFLOWSUPERSONIC);\n"
    "    }\n"
    "  );\n",_bdry_ID(BDRY_INFLOWSUPERSONIC),_bdry_ID(BDRY_INFLOWSUBSONIC1),
             _bdry_ID(BDRY_INFLOWSUBSONICMASSFLOWFIXED1),_bdry_ID(BDRY_INFLOWINJECTION1),_bdry_ID(BDRY_INFLOWINJECTIONCHOKED1),_bdry_ID(BDRY_OUTFLOWSUPERSONIC1),
             _bdry_ID(BDRY_OUTFLOWSUBSONIC1),_bdry_ID(BDRY_OUTFLOWSUBSONICMFIXED1),
             _bdry_ID(BDRY_SYMMETRICAL2),_bdry_ID(BDRY_SYMMETRICAL1),_bdry_ID(BDRY_WALLTFIXED1),_bdry_ID(BDRY_WALLTFIXEDINJECTION1),
             _bdry_ID(BDRY_WALLABLATION1),_bdry_ID(BDRY_WALLADIABATIC1),_bdry_ID(BDRY_SLIPWALL1),_bdry_ID(BDRY_FREESTREAM1)
  );
}


void add_bdry_types_fluid_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"BDRY_INFLOWSUPERSONIC",   BDRY_INFLOWSUPERSONIC);
  add_int_to_codex(codex,"BDRY_INFLOWSUBSONIC1",   BDRY_INFLOWSUBSONIC1);
  add_int_to_codex(codex,"BDRY_INFLOWSUBSONICMASSFLOWFIXED1",  BDRY_INFLOWSUBSONICMASSFLOWFIXED1 );
  add_int_to_codex(codex,"BDRY_INFLOWINJECTION1",  BDRY_INFLOWINJECTION1 );
  add_int_to_codex(codex,"BDRY_INFLOWINJECTIONCHOKED1",  BDRY_INFLOWINJECTIONCHOKED1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUPERSONIC1",  BDRY_OUTFLOWSUPERSONIC1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUBSONIC1",  BDRY_OUTFLOWSUBSONIC1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUBSONICMFIXED1",  BDRY_OUTFLOWSUBSONICMFIXED1 );
  add_int_to_codex(codex,"BDRY_SYMMETRICAL1", BDRY_SYMMETRICAL1  );
  add_int_to_codex(codex,"BDRY_SYMMETRICAL2", BDRY_SYMMETRICAL2  );
  add_int_to_codex(codex,"BDRY_WALLTFIXED1", BDRY_WALLTFIXED1  );
  add_int_to_codex(codex,"BDRY_WALLTFIXEDINJECTION1", BDRY_WALLTFIXEDINJECTION1  );
  add_int_to_codex(codex,"BDRY_WALLABLATION1", BDRY_WALLABLATION1  );
  add_int_to_codex(codex,"BDRY_INFLOWINJECTION1", BDRY_INFLOWINJECTION1  );
  add_int_to_codex(codex,"BDRY_WALLADIABATIC1",  BDRY_WALLADIABATIC1 );
  add_int_to_codex(codex,"BDRY_SLIPWALL1", BDRY_SLIPWALL1  );
  add_int_to_codex(codex,"BDRY_FREESTREAM1",   BDRY_FREESTREAM1);
}



/* turn off cross-derivative terms within FDS at bdry nodes except for slip walls and symmetry surfaces */
bool is_node_bdry_no_cross(np_t np, int TYPELEVEL){
  bool RET;
  RET=FALSE;
  if (is_node_bdry(np,TYPELEVEL)
      && !(_node_type(np,TYPELEVEL)==BDRY_SLIPWALL1
       || _node_type(np,TYPELEVEL)==BDRY_SYMMETRICAL1
       || _node_type(np,TYPELEVEL)==BDRY_SYMMETRICAL2
      )
   ) RET=TRUE;
  return(RET);
}



static void update_bdry_back_pressure(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, int ACCURACY){
  long spec,dim;
  spec_t rhok;
  double Pstar,k,psi;
  dim_t V;
  bool ref_flag;

  assert_np(np[lA],is_node_resumed(np[lA]));
  Pstar=_bdry_param(np,gl,lA,0,TYPELEVEL_FLUID_WORK);
  for (spec=0; spec<ns; spec++){
    rhok[spec]=_f_extrapol(ACCURACY,_rhok(np[lB],spec),_rhok(np[lC],spec));
  }
  reformat_rhok(gl,rhok,"_bdry",&ref_flag);

  for (dim=0; dim<nd; dim++){
    V[dim]=_f_extrapol(ACCURACY,_V(np[lB],dim),_V(np[lC],dim));
  }
  k=_f_extrapol(ACCURACY,_k(np[lB]),_k(np[lC]));
  psi=_f_extrapol(ACCURACY,_psi(np[lB]),_psi(np[lC]));
  find_U_3(np, lA, gl, rhok, V, Pstar, k, psi);

}


static void update_bdry_outflow_Mach(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta,  int ACCURACY){
  long spec,dim;
  spec_t rhok;
  double P, Mdes, Mextr;
  double V_temp, Vmag,k,psi;
  dim_t V;
  bool ref_flag;

  assert_np(np[lA],is_node_resumed(np[lA]));

  /* Store desired Mach number (Mdes) at C */
  Mdes=_bdry_param(np,gl,lA,0,TYPELEVEL_FLUID_WORK);


  /* Extrapolate Mach number (Mextr) at C from A and B */
  P=_f_extrapol(ACCURACY,_P(np[lB],gl), _P(np[lC],gl));

  for (spec=0; spec<ns; spec++){
    rhok[spec]=_f_extrapol(ACCURACY,_rhok(np[lB],spec),_rhok(np[lC],spec));
  }
  reformat_rhok(gl,rhok,"_bdry",&ref_flag);

  for (dim=0; dim<nd; dim++){
    V[dim]=_f_extrapol(ACCURACY,_V(np[lB],dim),_V(np[lC],dim));
  }

  k=_f_extrapol(ACCURACY,_k(np[lB]),_k(np[lC]));
  psi=_f_extrapol(ACCURACY,_psi(np[lB]),_psi(np[lC]));
  find_U_3(np, lA, gl, rhok, V, P, k,psi);

  V_temp=0.0;
  for (dim=0; dim<nd; dim++){
    V_temp=V_temp+sqr(_V(np[lA],dim));
  }
  Vmag = sqrt(V_temp);
  Mextr=Vmag/_a(np[lA],gl);


  /* Modify velocity components to give Mdes */
  for (dim=0; dim<nd; dim++){
    V[dim]=( _V(np[lA],dim)*Mdes/Mextr );
  }
  find_U_3(np, lA, gl, rhok, V, P, k, psi);
}




static void update_bdry_inflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                               bool BDRYDIRECFOUND){
  double P,T,k,psi;
  spec_t w;
  dim_t V;
  long spec,dim;

  assert_np(np[lA],is_node_resumed(np[lA]));
  P=_P(np[lA],gl);
  T=_T(np[lA],gl);
  for (spec=0; spec<ns; spec++) w[spec]=_w(np[lA],spec);
  for (dim=0; dim<nd; dim++) V[dim]=_V(np[lA],dim);
  k=_k(np[lA]);
  psi=_psi(np[lA]);
  find_U_2(np, lA, gl, w, V, P, T, k, psi);
}



static void update_bdry_inflow_reservoir(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                         bool BDRYDIRECFOUND, int ACCURACY){
  spec_t w;
  double k,psi,P,T;
  dim_t V;
  find_w_V_P_T_bdry_inflow_reservoir(np, gl, lA, lB, lC, ACCURACY, w, V, &P, &T);
  k=_k(np[lA]);
  psi=_psi(np[lA]);
  find_U_2(np, lA, gl, w, V, P, T, k, psi);
}


static void update_bdry_inflow_reservoir_2(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                           bool BDRYDIRECFOUND, int ACCURACY){
  spec_t w;
  double P,T,k,psi;
  dim_t V;
  find_w_V_P_T_bdry_inflow_reservoir_2(np, gl, lA, lB, lC, ACCURACY, w, V, &P, &T);
  k=_k(np[lA]);
  psi=_psi(np[lA]);
  find_U_2(np, lA, gl, w, V, P, T, k, psi);
}


static void update_bdry_outflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                bool BDRYDIRECFOUND, int ACCURACY){
  double P,T,k,psi;
  spec_t w;
  dim_t V;
  long spec,dim;
  bool ref_flag;

  assert_np(np[lA],is_node_resumed(np[lA]));
  P=_f_extrapol(ACCURACY,_P(np[lB],gl),_P(np[lC],gl));
  T=_f_extrapol(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
  for (spec=0; spec<ns; spec++){
    w[spec]=_f_extrapol(ACCURACY,_w(np[lB],spec),_w(np[lC],spec));
  }
  reformat_w(gl,w,"_bdry",&ref_flag);
  for (dim=0; dim<nd; dim++){
    V[dim]=_f_extrapol(ACCURACY,_V(np[lB],dim),_V(np[lC],dim));
  }
  k=_f_extrapol(ACCURACY,_k(np[lB]),_k(np[lC]));
  psi=_f_extrapol(ACCURACY,_psi(np[lB]),_psi(np[lC]));
  find_U_2(np, lA, gl, w, V, P, T, k, psi);
}


static void update_bdry_symmetrical(np_t *np, gl_t *gl, long lA, long lB, long lC, long lD, long theta, long thetasgn,
                                    bool BDRYDIRECFOUND, int ACCURACY){
  spec_t w;
  double k,psi,P,T;
  dim_t V;
  bool ref_flag;
  find_w_V_P_T_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND,  ACCURACY, w, V, &P, &T);
  reformat_w(gl,w,"_bdry",&ref_flag);

  k=_f_symmetry(ACCURACY,_k(np[lB]),_k(np[lC]));
  psi=_f_symmetry(ACCURACY,_psi(np[lB]),_psi(np[lC]));

  find_U_2(np, lA, gl, w, V, P, T, k, psi);
}


static void update_bdry_freestream(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                   bool BDRYDIRECFOUND, int ACCURACY){
  dim_t V;
  long spec;
  double k,psi,T,P;
  spec_t w;
  bool OUTFLOW, FREESTREAM;

  find_V_P_T_bdry_freestream(np, gl, lA, lB, lC, theta, thetasgn, ACCURACY, V, &P, &T, &OUTFLOW, &FREESTREAM);
  if (OUTFLOW){
    for (spec=0; spec<ns; spec++) w[spec]=_w(np[lB],spec);
    k=_k(np[lB]);
    psi=_psi(np[lB]);
  } else {
    for (spec=0; spec<ns; spec++) w[spec]=_w(np[lA],spec);
    k=_k(np[lA]);
    psi=_psi(np[lA]);
  }
  find_U_2(np, lA, gl, w, V, P, T, k, psi);

}


static void update_bdry_wall(np_t *np, gl_t *gl, long lA, long lB, long lC,
                            long theta, long thetasgn,
                            bool ADIABATIC, bool INJECTION, bool BDRYDIRECFOUND, int ACCURACY){
  spec_t wwall;
  double kwall,psiwall,Twall,Pwall;
  long dim,spec;
  dim_t Vwall;
  spec_t nukA,nukB;
  bool ref_flag;

  if (ADIABATIC) {
    Twall=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
  } else {
    Twall=_bdry_param(np,gl,lA,0,TYPELEVEL_FLUID_WORK);
  }
  /* clip Twall so that it remains within the user-specified bounds */
  Twall=max(gl->model.fluid.Twmin,min(gl->model.fluid.Twmax,Twall));
  
  for (dim=0; dim<nd; dim++){
    Vwall[dim]=0.0e0;
  }
  for (spec=0; spec<ns; spec++){
    wwall[spec]=_f_symmetry(ACCURACY,_w(np[lB],spec),_w(np[lC],spec));
  }

  if (INJECTION) {
    for (spec=0; spec<ns; spec++){
      nukA[spec]=_nustar(np,lA,gl,spec);
      nukB[spec]=_nustar(np,lB,gl,spec); 
    }
    update_w_V_at_injection_wall(np, gl, lA, lB, lC, nukA, nukB, 1, np[lA].numbdryparam-1, wwall, Vwall);
  }

  for (spec=0; spec<ncs; spec++) 
    wwall[spec]=0.0;

  reformat_w(gl,wwall,"_bdry",&ref_flag);

  find_k_psi_bdry_wall(np, gl, lA, lB, lC, &kwall, &psiwall);

  find_Pstar_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY, &Pwall);

  find_U_2(np, lA,gl,wwall,Vwall,Pwall,Twall,kwall,psiwall);
}



/*
R.N. Gupta, Aerothermodynamic Analysis of Stardust Sample Return
Capsule with Coupled Radiation and Ablation
 JOURNAL OF SPACECRAFT AND ROCKETS
Vol. 37, No. 4, July–August 2000

Phenolic Impregnated Ceramic Ablator(PICA)

0.001 atm < p < 1.00 atm

sublimation temperature and heat of ablation for PICA (with 92% carbon,
4.9% oxygen, 2.2% hydrogen, and 0.9% nitrogen by mass)

 
 */
static void update_bdry_wall_ablation(np_t *np, gl_t *gl, long lA, long lB, long lC,
                            long theta, long thetasgn,
                            bool BDRYDIRECFOUND, int ACCURACY){
  spec_t wwall;
  double kwall,psiwall,Twall,Pwall;
  long dim,spec,j;
  double alpha1[6],alpha2[6],alpha3[6];
  dim_t Vwall;
  spec_t nukA,nukB;
  bool ref_flag;
  double wablationwall;
  long paramstart,paramend,numinj,inj;


  find_Pstar_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY, &Pwall);

  
  wablationwall=0.0;
  paramstart=1;
  paramend=np[lA].numbdryparam-1;
  if (mod(paramend-paramstart+1,2)!=0) fatal_error("Wrong number of extra parameters to ablation boundary condition.");
  numinj=round((double)(paramend-paramstart)/2.0);
  for (inj=0; inj<numinj; inj++){
    spec=round(_bdry_param(np,gl,lA,paramstart+inj*2,TYPELEVEL_FLUID_WORK))-1;
    if (spec<0 || spec>=ns) fatal_error("Wrong specification of the ablation boundary condition. The ablation species number is not within bounds.");
    wablationwall+=_w(np[lA],spec);
  }
  assert(wablationwall>=0.0);
  assert(wablationwall<=1.0001);
  
  
  alpha1[1]=3790.0;
  alpha1[2]=86.795;
  alpha1[3]=-2980.0;
  alpha1[4]=-8250.2;
  alpha1[5]=-7631.7;
  alpha2[1]=329.94;
  alpha2[2]=-66.703;
  alpha2[3]=-1524.6;
  alpha2[4]=-4340.9;
  alpha2[5]=-3885.7;
  alpha3[1]=20.386;
  alpha3[2]=-17.654;
  alpha3[3]=-268.62;
  alpha3[4]=-771.00;
  alpha3[5]=-684.89;

  Twall=0.0;
  for (j=1; j<=5; j++){
    Twall+=alpha1[j]*powint(log(wablationwall),j-1)
          +log(Pwall/101300.0)*alpha2[j]*powint(log(wablationwall),j-1)
          +sqr(log(Pwall/101300.0))*alpha3[j]*powint(log(wablationwall),j-1);
  }
  /* clip Twall so that it remains within the user-specified bounds */
  Twall=max(gl->model.fluid.Twmin,min(gl->model.fluid.Twmax,Twall));
  
  for (dim=0; dim<nd; dim++){
    Vwall[dim]=0.0e0;
  }
  for (spec=0; spec<ns; spec++){
    wwall[spec]=_f_symmetry(ACCURACY,_w(np[lB],spec),_w(np[lC],spec));
  }

  for (spec=0; spec<ns; spec++){
    nukA[spec]=_nustar(np,lA,gl,spec);
    nukB[spec]=_nustar(np,lB,gl,spec); 
  }
  // here need to change np[lA].bdryparam before sending it to update_w_V_at_injection_wall
  update_w_V_at_injection_wall(np, gl, lA, lB, lC, nukA, nukB, 1, np[lA].numbdryparam-1, wwall, Vwall);

  for (spec=0; spec<ncs; spec++) 
    wwall[spec]=0.0;

  reformat_w(gl,wwall,"_bdry",&ref_flag);

  find_k_psi_bdry_wall(np, gl, lA, lB, lC, &kwall, &psiwall);


  find_U_2(np, lA,gl,wwall,Vwall,Pwall,Twall,kwall,psiwall);
}



static void update_bdry_inflow_injection_old(np_t *np, gl_t *gl, long lA, long lB, long lC,
                                         long theta, long thetasgn,
                                         bool BDRYDIRECFOUND, int ACCURACY){
  spec_t wwall;
  double Rk,cpk,rhowall,Tstag,kwall,psiwall,Twall,Pwall,mdot;
  long dim,spec,specinj;
  dim_t Vwall,n;
  bool ref_flag;

//  Twall=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
  Tstag=_bdry_param(np,gl,lA,0,TYPELEVEL_FLUID_WORK);
  specinj=round(_bdry_param(np,gl,lA,1,TYPELEVEL_FLUID_WORK))-1;
  mdot=_bdry_param(np,gl,lA,2,TYPELEVEL_FLUID_WORK);

  
  if (!find_unit_vector_normal_to_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK, n)) 
    fatal_error("Problem finding unit normal vector in update_bdry_inflow_injection()");
  
  
  Pwall=_Pstar(np[lB],gl);

/*
  Twall=Tstag-Vwallmag2/(2.0*cp);
  (mdot/rhowall)=Vwallmag;
  Twall=Tstag-sqr(mdot/rhowall)/(2.0*cp);
  Pwall/(R*Twall)=rhowall;
  Twall=Tstag-sqr(mdot*(R*Twall)/Pwall)/(2.0*cp);
  sqr(Twall)*sqr(mdot*R/Pwall)/(2.0*cp)+Twall-Tstag=0;
  a*T^2+b*T+c=0;
  a=sqr(mdot*R/Pwall)/(2.0*cp);
  b=1
  c=-Tstag
  Twall=(-b+sqrt(b*b-4*a*c))/(2*a);
  Twall=(-1.0+sqrt(1.0+4*sqr(mdot*R/Pwall)/(2.0*cp)*Tstag))/(2*sqr(mdot*R/Pwall)/(2.0*cp));
*/
  Rk=kB/_m(specinj);
  cpk=_cpk_from_T_equilibrium(specinj, _T(np[lA],gl));
  Twall=(-1.0+sqrt(1.0+4.0*sqr(mdot*Rk/Pwall)/(2.0*cpk)*Tstag))/(2.0*sqr(mdot*Rk/Pwall)/(2.0*cpk));
  
  rhowall=Pwall/Rk/Twall;
  for (dim=0; dim<nd; dim++){
    Vwall[dim]=n[dim]*mdot/rhowall;
  }

  for (spec=0; spec<ns; spec++){
    if (spec==specinj) wwall[spec]=1.0; 
      else wwall[spec]=1e-99;
  }



  reformat_w(gl,wwall,"_bdry",&ref_flag);
  kwall=_k(np[lB]);
  psiwall=_psi(np[lB]);


  find_U_2(np, lA,gl,wwall,Vwall,Pwall,Twall,kwall,psiwall);
}




static void update_bdry_inflow_injection(np_t *np, gl_t *gl, long lA, long lB, long lC,
                                         long theta, long thetasgn,
                                         bool FORCECHOKING, bool BDRYDIRECFOUND, int ACCURACY){
  spec_t wwall;
  double Rk,cpk,gammak,Mwall,Tstag,Pstag,kwall,psiwall,Twall,Pwall;
  long dim,spec,specinj;
  dim_t Vwall,n;
  bool ref_flag;

//  Twall=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
  Tstag=_bdry_param(np,gl,lA,0,TYPELEVEL_FLUID_WORK);
  Pstag=_bdry_param(np,gl,lA,1,TYPELEVEL_FLUID_WORK);
  specinj=round(_bdry_param(np,gl,lA,2,TYPELEVEL_FLUID_WORK))-1;

  
  if (!find_unit_vector_normal_to_boundary_plane(np, gl, lA, lB, lC, TYPELEVEL_FLUID_WORK, n)) 
    fatal_error("Problem finding unit normal vector in update_bdry_inflow_injection()");
  
  
  Pwall=_Pstar(np[lB],gl);

  Rk=kB/_m(specinj);
  cpk=_cpk_from_T_equilibrium(specinj, _T(np[lA],gl));
  gammak=cpk/(cpk-Rk);
  
  // first assume choking at the injection point (M=1)
  Mwall=1.0;
  Pwall=Pstag/pow(1.0+(gammak-1.0)/2.0,gammak/(gammak-1.0));
  // check if choking is possible; if not, set the pressure at the boundary equal to the pressure nearby
  if (!FORCECHOKING && Pwall<_Pstar(np[lC],gl)){
    Pwall=min(Pstag,_Pstar(np[lC],gl));
    // find Mach number at the wall knowing the pressure and stagnation pressure
    if ((pow(Pstag/Pwall,(gammak-1.0)/gammak)-1.0)*2.0/(gammak-1.0)<0.0) {
      wfprintf(stderr,"\n\n gammak=%E  Pstag=%E  Pwall=%E \n\n",gammak,Pstag,Pwall);
    }
    assert_np(np[lA],(pow(Pstag/Pwall,(gammak-1.0)/gammak)-1.0)*2.0/(gammak-1.0)>=0.0);
    Mwall=sqrt((pow(Pstag/Pwall,(gammak-1.0)/gammak)-1.0)*2.0/(gammak-1.0));
    assert_np(np[lA],Mwall<1.0);
  }
  

  Twall=Tstag/(1.0+(gammak-1.0)/2.0*sqr(Mwall));
  for (dim=0; dim<nd; dim++){
    Vwall[dim]=n[dim]*Mwall*sqrt(gammak*Rk*Twall);
  }

  for (spec=0; spec<ns; spec++){
    if (spec==specinj) wwall[spec]=1.0; 
      else wwall[spec]=1e-99;
  }



  reformat_w(gl,wwall,"_bdry",&ref_flag);
  kwall=1.0e-10;
  psiwall=100.0*Mwall*sqrt(gammak*Rk*Twall); 


  find_U_2(np, lA,gl,wwall,Vwall,Pwall,Twall,kwall,psiwall);
}



bool is_node_bdry_symmetry_plane_fluid(np_t np){
  bool RET;
  if (is_node_bdry(np,TYPELEVEL_FLUID) && (_node_type(np,TYPELEVEL_FLUID)==BDRY_SYMMETRICAL1 
   || _node_type(np,TYPELEVEL_FLUID)==BDRY_SYMMETRICAL2)) RET=TRUE; else RET=FALSE;
  return(RET);
}


bool is_node_bdry_wall_fluid(np_t np, gl_t *gl){
  bool RET;
  if (is_node_bdry(np,TYPELEVEL_FLUID) && (_node_type(np,TYPELEVEL_FLUID)==BDRY_WALLADIABATIC1 
   || _node_type(np,TYPELEVEL_FLUID)==BDRY_WALLTFIXED1 || _node_type(np,TYPELEVEL_FLUID)==BDRY_WALLTFIXEDINJECTION1 
   || _node_type(np,TYPELEVEL_FLUID)==BDRY_WALLABLATION1)) RET=TRUE; else RET=FALSE;
  return(RET);    
}

void update_bdry_fluid(np_t *np, gl_t *gl, long lA, long lB, long lC, long lD, long theta, long thetasgn, bool BDRYDIRECFOUND, int TYPELEVEL){

  switch (_node_type(np[lA],TYPELEVEL)) {

    case BDRY_INFLOWSUPERSONIC:
      update_bdry_inflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND);
    break;

    case BDRY_OUTFLOWSUPERSONIC1:
      update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLTFIXED1:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLTFIXEDINJECTION1:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, TRUE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLABLATION1:
      update_bdry_wall_ablation(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_OUTFLOWSUBSONIC1:
      update_bdry_back_pressure(np, gl, lA, lB, lC, theta, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLADIABATIC1: 
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, TRUE, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;
    
    case BDRY_INFLOWSUBSONIC1:
      update_bdry_inflow_reservoir(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_SYMMETRICAL1 :
      if (BDRYDIRECFOUND)
         update_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_SLIPWALL1:
      if (BDRYDIRECFOUND)
        update_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_SYMMETRICAL2:
      if (BDRYDIRECFOUND)
        update_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_SECONDORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_OUTFLOWSUBSONICMFIXED1:
      update_bdry_outflow_Mach(np, gl, lA, lB, lC, theta, ACCURACY_FIRSTORDER);
    break;

    case BDRY_INFLOWSUBSONICMASSFLOWFIXED1:
      update_bdry_inflow_reservoir_2(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_FREESTREAM1:
      if (BDRYDIRECFOUND)
        update_bdry_freestream(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_INFLOWINJECTION1:
      if (BDRYDIRECFOUND)
        update_bdry_inflow_injection(np, gl, lA, lB, lC, theta, thetasgn, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_INFLOWINJECTIONCHOKED1:
      if (BDRYDIRECFOUND)
        update_bdry_inflow_injection(np, gl, lA, lB, lC, theta, thetasgn, TRUE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    default:
      fatal_error("Boundary condition %ld not valid.",_node_type(np[lA],TYPELEVEL));
  }
}








