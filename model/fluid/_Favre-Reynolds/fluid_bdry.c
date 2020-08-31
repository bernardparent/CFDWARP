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
    "    BDRY_OUTFLOWSUPERSONIC1           %c   Outflow, supersonic\n"
    "    BDRY_OUTFLOWSUBSONIC1             %c   Outflow, subsonic, P fixed\n"
    "    BDRY_OUTFLOWSUBSONICMFIXED1       %c   Outflow, subsonic, M fixed\n"
    "    BDRY_SYMMETRICAL2                 %c   Symmetrical, 2nd order\n"
    "    BDRY_SYMMETRICAL1                 %c   Symmetrical, 1st order\n"
    "    BDRY_WALLTFIXED1                  %c   Wall, T specified, param Twall\n"
    "    BDRY_WALLADIABATIC1               %c   Wall, Adiabatic\n"
    "    BDRY_SLIPWALL1                    %c   Slip wall, 1st order, Adiabatic\n"
    "    BDRY_FREESTREAM1                  %c   Freestream, 1o, params Vx,Vy,"if3DL("Vz,")" P, T\n"
    "    _________________________________________________________________________________________\n"
    "    }\n"
    "    All(BDRY_WALLTFIXED1);\n"
    "    Plane(\"i\",is,BDRY_INFLOWSUPERSONIC);\n"
    "    Plane(\"i\",ie,BDRY_OUTFLOWSUPERSONIC1);\n"
    "    Twall=300.0; {K}\n"
    "    Plane(\"j\",js,BDRY_WALLTFIXED1,Twall);\n"
    "    Plane(\"j\",je,BDRY_WALLTFIXED1,Twall);\n"
#ifdef _3D
    "    Plane(\"k\",ks,BDRY_SYMMETRICAL2);\n"
    "    Plane(\"k\",ke,BDRY_SYMMETRICAL2);\n"
#endif
    "    {\n"
    "    Cut(is" if2DL(",js") if3DL(",ks") ",  ie" if2DL(",je") if3DL(",ke") ");\n"
    "    Region(is" if2DL(",js") if3DL(",ks") ",  ie" if2DL(",je") if3DL(",ke") ",  BDRY_INFLOWSUPERSONIC);\n"
    "    Link(i1" if2DL(",j1") if3DL(",k1") ",  i2" if2DL(",j2") if3DL(",k2") ");\n"
    "    }\n"
    "  );\n",_bdry_ID(BDRY_INFLOWSUPERSONIC),_bdry_ID(BDRY_INFLOWSUBSONIC1),
             _bdry_ID(BDRY_INFLOWSUBSONICMASSFLOWFIXED1),_bdry_ID(BDRY_OUTFLOWSUPERSONIC1),
             _bdry_ID(BDRY_OUTFLOWSUBSONIC1),_bdry_ID(BDRY_OUTFLOWSUBSONICMFIXED1),
             _bdry_ID(BDRY_SYMMETRICAL2),_bdry_ID(BDRY_SYMMETRICAL1),_bdry_ID(BDRY_WALLTFIXED1),
             _bdry_ID(BDRY_WALLADIABATIC1),_bdry_ID(BDRY_SLIPWALL1),_bdry_ID(BDRY_FREESTREAM1)
  );
}


void add_bdry_types_fluid_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"BDRY_INFLOWSUPERSONIC",   BDRY_INFLOWSUPERSONIC);
  add_int_to_codex(codex,"BDRY_INFLOWSUBSONIC1",   BDRY_INFLOWSUBSONIC1);
  add_int_to_codex(codex,"BDRY_INFLOWSUBSONICMASSFLOWFIXED1",  BDRY_INFLOWSUBSONICMASSFLOWFIXED1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUPERSONIC1",  BDRY_OUTFLOWSUPERSONIC1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUBSONIC1",  BDRY_OUTFLOWSUBSONIC1 );
  add_int_to_codex(codex,"BDRY_OUTFLOWSUBSONICMFIXED1",  BDRY_OUTFLOWSUBSONICMFIXED1 );
  add_int_to_codex(codex,"BDRY_SYMMETRICAL1", BDRY_SYMMETRICAL1  );
  add_int_to_codex(codex,"BDRY_SYMMETRICAL2", BDRY_SYMMETRICAL2  );
  add_int_to_codex(codex,"BDRY_WALLTFIXED1", BDRY_WALLTFIXED1  );
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
  Pstar=_Pstar(np[lA],gl);
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
  V_temp=0.0;
  for (dim=0; dim<nd; dim++){
    V_temp=V_temp+sqr(_V(np[lA],dim));
  }
  Vmag=sqrt(V_temp);
  Mdes=Vmag/_a(np[lA],gl);


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
                            bool ADIABATIC, bool BDRYDIRECFOUND, int ACCURACY){
  spec_t wwall;
  double kwall,psiwall,Twall,Pwall;
  long dim,spec;
  dim_t Vwall;
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
  for (spec=0; spec<ncs; spec++) 
    wwall[spec]=0.0;

  reformat_w(gl,wwall,"_bdry",&ref_flag);

  find_k_psi_bdry_wall(np, gl, lA, lB, lC, &kwall, &psiwall);

  find_Pstar_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY, &Pwall);

  find_U_2(np, lA,gl,wwall,Vwall,Pwall,Twall,kwall,psiwall);
}


bool is_node_bdry_symmetry_plane_fluid(np_t np){
  double RET;
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
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_OUTFLOWSUBSONIC1:
      update_bdry_back_pressure(np, gl, lA, lB, lC, theta, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLADIABATIC1: 
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, TRUE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
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

    default:
      fatal_error("Boundary condition %ld not valid.",_node_type(np[lA],TYPELEVEL));
  }
}








