// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2010-2011,2020-2021 Bernard Parent

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
#include <model/share/fluid_share.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
#include <model/emfield/_emfield.h>
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
#define BDRY_WALLTFIXEDCATALYTIC1 17
#define BDRY_WALLTFIXEDINJECTION1 18
#define BDRY_INFLOWINJECTION1 19


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
    "    BDRY_INFLOWSUBSONIC1              %c   Inflow, subsonic, Tstag, Pstag fixed, 1o\n"
    "    BDRY_INFLOWSUBSONICMASSFLOWFIXED1 %c   Inflow, subsonic, Pstag, Massflow/Area fixed, 1o\n"
    "    BDRY_INFLOWINJECTION1             %c   Inflow, param Tstag, Pstag, specCs \n"
    "    BDRY_OUTFLOWSUPERSONIC1           %c   Outflow, supersonic, 1o\n"
    "    BDRY_OUTFLOWSUBSONIC1             %c   Outflow, subsonic, P fixed, param P, 1o\n"
    "    BDRY_OUTFLOWSUBSONICMFIXED1       %c   Outflow, subsonic, M fixed, param M, 1o\n"
    "    BDRY_SYMMETRICAL2                 %c   Symmetrical, 2o\n"
    "    BDRY_SYMMETRICAL1                 %c   Symmetrical, 1o\n"
    "    BDRY_WALLTFIXED1                  %c   Wall, T specified, param Twall\n"
    "    BDRY_WALLTFIXEDCATALYTIC1         %c   Wall, T specified, param Twall, specO,specO2,gammaO, ...\n"
    "    BDRY_WALLTFIXEDINJECTION1         %c   Wall, T specified, param Twall, specCs,mdotCs[kg/m2s], ...\n"
    "    BDRY_WALLADIABATIC1               %c   Wall, Adiabatic\n"
    "    BDRY_SLIPWALL1                    %c   Slip wall, Adiabatic, 1o\n"
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
             _bdry_ID(BDRY_INFLOWSUBSONICMASSFLOWFIXED1),_bdry_ID(BDRY_INFLOWINJECTION1),_bdry_ID(BDRY_OUTFLOWSUPERSONIC1),
             _bdry_ID(BDRY_OUTFLOWSUBSONIC1),_bdry_ID(BDRY_OUTFLOWSUBSONICMFIXED1),
             _bdry_ID(BDRY_SYMMETRICAL2),_bdry_ID(BDRY_SYMMETRICAL1),_bdry_ID(BDRY_WALLTFIXED1),_bdry_ID(BDRY_WALLTFIXEDCATALYTIC1),
             _bdry_ID(BDRY_WALLTFIXEDINJECTION1),
             _bdry_ID(BDRY_WALLADIABATIC1),_bdry_ID(BDRY_SLIPWALL1),_bdry_ID(BDRY_FREESTREAM1)
  );
}


void add_bdry_types_fluid_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"BDRY_INFLOWSUPERSONIC",   BDRY_INFLOWSUPERSONIC);
  add_int_to_codex(codex,"BDRY_INFLOWSUBSONIC1",   BDRY_INFLOWSUBSONIC1);
  add_int_to_codex(codex,"BDRY_INFLOWSUBSONICMASSFLOWFIXED1",  BDRY_INFLOWSUBSONICMASSFLOWFIXED1 );
  add_int_to_codex(codex,"BDRY_INFLOWINJECTION1",  BDRY_INFLOWINJECTION1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUPERSONIC1",  BDRY_OUTFLOWSUPERSONIC1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUBSONIC1",  BDRY_OUTFLOWSUBSONIC1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUBSONICMFIXED1",  BDRY_OUTFLOWSUBSONICMFIXED1 );
  add_int_to_codex(codex,"BDRY_SYMMETRICAL1", BDRY_SYMMETRICAL1  );
  add_int_to_codex(codex,"BDRY_SYMMETRICAL2", BDRY_SYMMETRICAL2  );
  add_int_to_codex(codex,"BDRY_WALLTFIXED1", BDRY_WALLTFIXED1  );
  add_int_to_codex(codex,"BDRY_WALLTFIXEDCATALYTIC1", BDRY_WALLTFIXEDCATALYTIC1  );
  add_int_to_codex(codex,"BDRY_WALLTFIXEDINJECTION1", BDRY_WALLTFIXEDINJECTION1  );
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
  double P;
  dim_t V;
  bool ref_flag;

  assert_np(np[lA],is_node_resumed(np[lA]));
  P=_bdry_param(np,gl,lA,0,TYPELEVEL_FLUID_WORK);
  for (spec=0; spec<ns; spec++){
    rhok[spec]=_f_extrapol(ACCURACY,_rhok(np[lB],spec),_rhok(np[lC],spec));
  }

  for (dim=0; dim<nd; dim++){
    V[dim]=_f_extrapol(ACCURACY,_V(np[lB],dim),_V(np[lC],dim));
  }
  reformat_rhok(gl,rhok,"_bdry",&ref_flag);

  find_U_3(np, lA, gl, rhok, V, P);
}


static void update_bdry_outflow_Mach(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, int ACCURACY){
  long spec,dim;
  spec_t rhok;
  double P, Mdes, Mextr;
  double V_temp, Vmag;
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

  find_U_3(np, lA, gl, rhok, V, P);

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
  find_U_3(np, lA, gl, rhok, V, P);
}


static void update_bdry_inflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                               bool BDRYDIRECFOUND){
  double P,T;
  spec_t w;
  dim_t V;
  long spec,dim;

  assert_np(np[lA],is_node_resumed(np[lA]));
  P=_P(np[lA],gl);
  T=_T(np[lA],gl);
  for (spec=0; spec<ns; spec++) w[spec]=_w(np[lA],spec);
  for (dim=0; dim<nd; dim++) V[dim]=_V(np[lA],dim);
  find_U_2(np, lA, gl, w, V, P, T);
}


static void update_bdry_inflow_reservoir(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                         bool BDRYDIRECFOUND, int ACCURACY){
  spec_t w;
  double P,T;
  dim_t V;
  find_w_V_P_T_bdry_inflow_reservoir(np, gl, lA, lB, lC, ACCURACY, w, V, &P, &T);
  find_U_2(np, lA, gl, w, V, P, T);
}


static void update_bdry_inflow_reservoir_2(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                           bool BDRYDIRECFOUND, int ACCURACY){
  spec_t w;
  double P,T;
  dim_t V;
  find_w_V_P_T_bdry_inflow_reservoir_2(np, gl, lA, lB, lC, ACCURACY, w, V, &P, &T);
  find_U_2(np, lA, gl, w, V, P, T);
}


static void update_bdry_outflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                bool BDRYDIRECFOUND, int ACCURACY){
  double P,T;
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
  for (dim=0; dim<nd; dim++){
    V[dim]=_f_extrapol(ACCURACY,_V(np[lB],dim),_V(np[lC],dim));
  }

  reformat_w(gl,w,"_bdry",&ref_flag);
  find_U_2(np, lA, gl, w, V, P, T);
}


static void update_bdry_symmetrical(np_t *np, gl_t *gl, long lA, long lB, long lC, long lD, long theta, long thetasgn,
                                    bool BDRYDIRECFOUND, int ACCURACY){
  spec_t w;
  double P,T;
  dim_t V;
  bool ref_flag;
  find_w_V_P_T_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND,  ACCURACY, w, V, &P, &T);
  reformat_w(gl,w,"_bdry",&ref_flag);
//  assert_np(np[lA],is_node_resumed(np[lA]));
  find_U_2(np, lA, gl, w, V, P, T);
}


static void update_bdry_freestream(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                   bool BDRYDIRECFOUND, int ACCURACY){
  dim_t V;
  long spec;
  double T,P;
  spec_t w;
  bool OUTFLOW, FREESTREAM;

  find_V_P_T_bdry_freestream(np, gl, lA, lB, lC, theta, thetasgn, ACCURACY, V, &P, &T, &OUTFLOW, &FREESTREAM);
  if (OUTFLOW){
    for (spec=0; spec<ns; spec++) w[spec]=_w(np[lB],spec);
  } else {
    for (spec=0; spec<ns; spec++) w[spec]=_w(np[lA],spec);
  }
  find_U_2(np, lA, gl, w, V, P, T);

}


static void update_bdry_wall(np_t *np, gl_t *gl, long lA, long lB, long lC,
                            long theta, long thetasgn,
                            bool ADIABATIC, bool CATALYTIC, bool INJECTION, bool BDRYDIRECFOUND, int ACCURACY){
  spec_t wwall;
  double Twall,Pwall;
  long dim,spec,specneutral;
  spec_t nukA,nukB;
  bool ref_flag;
  dim_t Vwall;

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
//    wwall[spec]=_f_symmetry(ACCURACY,_w(np[lB],spec),_w(np[lC],spec));
    wwall[spec]=_w(np[lB],spec);
  }

  if (CATALYTIC) update_w_at_catalytic_wall(np, gl, lA, lB, lC, theta, thetasgn, Twall, Twall, 1, np[lA].numbdryparam-1, wwall);

  if (INJECTION) {
    for (spec=0; spec<ns; spec++){
      nukA[spec]=_nu(np[lA],gl,spec);
      nukB[spec]=_nu(np[lB],gl,spec); 
    }
    update_w_V_at_injection_wall(np, gl, lA, lB, lC, nukA, nukB, 1, np[lA].numbdryparam-1, wwall, Vwall);
  }
  
  for (spec=0; spec<ncs; spec++){
    // set ion and electron densities to zero at the surface
    wwall[spec]=0.0;
    // make sure that no net mass flow (due to diffusion) goes through the boundary
    if (speciestype[spec]==SPECIES_IONPLUS) {
      if (!find_neutral_spec_from_ion_spec(spec,&specneutral))
        fatal_error("Couldn't find a neutral species associated with positive ion species %ld.",spec); 
      wwall[specneutral]+=_w_product_at_catalytic_wall(np, gl, lA, lB, lC, theta, thetasgn, spec, specneutral, 1.0)-_w(np[lB],specneutral);
    }
  }

  reformat_w(gl,wwall,"_bdry",&ref_flag);

  find_Pstar_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY, &Pwall);

  find_U_2(np, lA,gl,wwall,Vwall,Pwall,Twall);

}

static void update_bdry_inflow_injection(np_t *np, gl_t *gl, long lA, long lB, long lC,
                                         long theta, long thetasgn,
                                         bool BDRYDIRECFOUND, int ACCURACY){
  spec_t wwall;
  double Rk,cpk,gammak,Mwall,Tstag,Pstag,Twall,Pwall;
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
  if (Pwall<_Pstar(np[lC],gl)){
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


  find_U_2(np, lA,gl,wwall,Vwall,Pwall,Twall);
}


bool is_node_bdry_symmetry_plane_fluid(np_t np){
  bool RET;
  if (is_node_bdry(np,TYPELEVEL_FLUID) && (_node_type(np,TYPELEVEL_FLUID)==BDRY_SYMMETRICAL1 
   || _node_type(np,TYPELEVEL_FLUID)==BDRY_SYMMETRICAL2)) RET=TRUE; else RET=FALSE;
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
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, FALSE, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLTFIXEDCATALYTIC1:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, TRUE, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLTFIXEDINJECTION1:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, FALSE, TRUE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLADIABATIC1: 
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, TRUE, FALSE, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_OUTFLOWSUBSONIC1:
      update_bdry_back_pressure(np, gl, lA, lB, lC, theta, ACCURACY_FIRSTORDER);
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
      update_bdry_outflow_Mach(np, gl, lA, lB, lC,  theta, ACCURACY_FIRSTORDER);
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
        update_bdry_inflow_injection(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    default:
      fatal_error("Boundary condition %ld not valid.",_node_type(np[lA],TYPELEVEL));
  }
}


