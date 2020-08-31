// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2010-2011,2020 Bernard Parent

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
#include <model/metrics/_metrics.h>
#include <model/share/fluid_share.h>
#include <model/_model.h>
#include <model/share/model_share.h>
#include <soap.h>


#define BDRY_INFLOWSUPERSONIC 0
#define BDRY_OUTFLOWSUPERSONIC1 1
#define BDRY_SYMMETRICAL1 2
#define BDRY_SYMMETRICAL2 3
#define BDRY_WALLTFIXED1 4
#define BDRY_WALLADIABATIC1 5
#define BDRY_INFLOWSUBSONIC1 6
#define BDRY_OUTFLOWSUBSONIC1 7
#define BDRY_SYMMETRICAL3 8
#define BDRY_FREESTREAM1 9
#define BDRY_WALLTFIXED2 10
#define BDRY_WALLADIABATIC2 11

#define SYMMETRY_TOTAL_ENTHALPY_CONSERVED FALSE

void write_bdry_fluid_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n",_FLUID_ACTIONNAME);
  wfprintf(*controlfile,
    "    {\n"
    "    _________________________________________________________________________________________\n"
    "\n" 
    "    Boundary Condition Type    ID    Description\n"
    "    _________________________________________________________________________________________\n"
    "\n" 
    "    BDRY_INFLOWSUPERSONIC      %c     Inflow, supersonic\n"
    "    BDRY_OUTFLOWSUPERSONIC1    %c     Outflow, supersonic, 1o\n"
    "    BDRY_SYMMETRICAL1          %c     Symmetrical, 1o\n"
    "    BDRY_SYMMETRICAL2          %c     Symmetrical, 2o\n"
    "    BDRY_SYMMETRICAL3          %c     Symmetrical, 3o\n"
    "    BDRY_WALLTFIXED1           %c     Wall, T specified, 1o, param Twall\n"
    "    BDRY_WALLTFIXED2           %c     Wall, T specified, 2o, param Twall\n"
    "    BDRY_WALLADIABATIC1        %c     Wall, adiabatic, 1o\n"
    "    BDRY_WALLADIABATIC2        %c     Wall, adiabatic, 2o\n"
    "    BDRY_INFLOWSUBSONIC1       %c     Subsonic Inflow 1o (Constant Tstag, Pstag at inflow)\n"
    "    BDRY_OUTFLOWSUBSONIC1      %c     Subsonic Outflow 1o (Constant P at outflow)\n"
    "    BDRY_FREESTREAM1           %c     Freestream, 1o, params Vx,Vy,"if3DL("Vz,")" P, T\n"
    "    _________________________________________________________________________________________\n"
    "    }\n"
    "    All(BDRY_OUTFLOWSUPERSONIC1);\n"
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
    "  );\n",_bdry_ID(BDRY_INFLOWSUPERSONIC),_bdry_ID(BDRY_OUTFLOWSUPERSONIC1),
             _bdry_ID(BDRY_SYMMETRICAL1),_bdry_ID(BDRY_SYMMETRICAL2),_bdry_ID(BDRY_SYMMETRICAL3),
             _bdry_ID(BDRY_WALLTFIXED1),_bdry_ID(BDRY_WALLTFIXED2),_bdry_ID(BDRY_WALLADIABATIC1),
             _bdry_ID(BDRY_WALLADIABATIC2),_bdry_ID(BDRY_INFLOWSUBSONIC1),
             _bdry_ID(BDRY_OUTFLOWSUBSONIC1),_bdry_ID(BDRY_FREESTREAM1)
  );
}


void add_bdry_types_fluid_to_codex(SOAP_codex_t *codex){
  add_int_to_codex(codex,"BDRY_INFLOWSUPERSONIC",   BDRY_INFLOWSUPERSONIC);
  add_int_to_codex(codex,"BDRY_OUTFLOWSUPERSONIC1",   BDRY_OUTFLOWSUPERSONIC1);
  add_int_to_codex(codex,"BDRY_SYMMETRICAL1",   BDRY_SYMMETRICAL1 );
  add_int_to_codex(codex,"BDRY_SYMMETRICAL2",   BDRY_SYMMETRICAL2);
  add_int_to_codex(codex,"BDRY_SYMMETRICAL3",   BDRY_SYMMETRICAL3);
  add_int_to_codex(codex,"BDRY_WALLTFIXED1",   BDRY_WALLTFIXED1);
  add_int_to_codex(codex,"BDRY_WALLTFIXED2",   BDRY_WALLTFIXED2);
  add_int_to_codex(codex,"BDRY_WALLADIABATIC1",   BDRY_WALLADIABATIC1);
  add_int_to_codex(codex,"BDRY_WALLADIABATIC2",   BDRY_WALLADIABATIC2);
  add_int_to_codex(codex,"BDRY_INFLOWSUBSONIC1",   BDRY_INFLOWSUBSONIC1);
  add_int_to_codex(codex,"BDRY_OUTFLOWSUBSONIC1",   BDRY_OUTFLOWSUBSONIC1);
  add_int_to_codex(codex,"BDRY_FREESTREAM1",   BDRY_FREESTREAM1);
}


/* turn off cross-derivative terms within FDS at bdry nodes except for slip walls */
bool is_node_bdry_no_cross(np_t np, int TYPELEVEL){
  bool RET;
  RET=FALSE;
  if (is_node_bdry(np,TYPELEVEL)
      && !(_node_type(np,TYPELEVEL)==BDRY_SYMMETRICAL1 
       || _node_type(np,TYPELEVEL)==BDRY_SYMMETRICAL2
       || _node_type(np,TYPELEVEL)==BDRY_SYMMETRICAL3
      )
   ) RET=TRUE;
  return(RET);
}


static void update_bdry_inflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                               bool BDRYDIRECFOUND){
  double P,T;
  dim_t V;
  long dim;

  assert_np(np[lA],is_node_resumed(np[lA]));
  P=_P(np[lA],gl);
  T=_T(np[lA],gl);
  for (dim=0; dim<nd; dim++) V[dim]=_V(np[lA],dim);
  //printf("---%E  %E  %E  %E\n",P,T,V[0],V[1]);
  find_U_2(np, lA, gl, V, P, T);
}


static void update_bdry_outflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                bool BDRYDIRECFOUND, int ACCURACY){
  double P,T;
  dim_t V;
  long dim;

  assert_np(np[lA],is_node_resumed(np[lA]));
  P=_f_extrapol(ACCURACY,_P(np[lB],gl),_P(np[lC],gl));
  T=_f_extrapol(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
  for (dim=0; dim<nd; dim++){
    V[dim]=_f_extrapol(ACCURACY,_V(np[lB],dim),_V(np[lC],dim));
  }
  find_U_2(np, lA, gl, V, P, T);
}


static void update_bdry_symmetrical(np_t *np, gl_t *gl, long lA, long lB, long lC, long lD, long theta, long thetasgn,
                                    bool BDRYDIRECFOUND, int ACCURACY){
  dim_t V,Vstar;
  long dim,i,j;
  double h,T,P;
  dim2_t chiinv;

  T=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl),_T(np[lD],gl));
  for (dim=0; dim<nd; dim++) Vstar[dim]=_f_symmetry(ACCURACY,_Vstar(np[lB],dim),_Vstar(np[lC],dim),_Vstar(np[lD],dim));

  Vstar[theta]=0.0e0;
  find_chi_inverse(np[lA].bs->X,chiinv);
  for (i=0; i<nd; i++){
    V[i]=0.0e0;
    for (j=0; j<nd; j++){
      V[i]=V[i]+chiinv[i][j]*Vstar[j];
    }
  }
  if (!BDRYDIRECFOUND) for (dim=0; dim<nd; dim++) V[dim]=0.0;
  find_P_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY, &P);
  if (SYMMETRY_TOTAL_ENTHALPY_CONSERVED){  
    /* impose conservation of total enthalpy */
    h=_ht(np[lB],gl);
    for (i=0; i<nd; i++) h-=0.5*sqr(V[i]);
    /* find TC from hC */
    T=max(gl->model.fluid.Tmin,h/((gl->model.fluid.gamma*gl->model.fluid.R)/(gl->model.fluid.gamma-1.0)));
  }
  
  find_U_2(np, lA, gl, V, P, T);
}


static void update_bdry_subsonic_inflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                        bool BDRYDIRECFOUND, int ACCURACY){
  dim_t Vnew,Vstarnew;
  long dim,i,j;
  double Rgas,Tnew,Pnew,gamma,Mnew,Mold,qnew,Tstag,Pstag;
  dim2_t chiinv;

  for (dim=0; dim<nd; dim++){
    Vstarnew[dim]=_f_symmetry(ACCURACY,_Vstar(np[lB],dim),_Vstar(np[lC],dim));
  }
  if (Vstarnew[theta]*thetasgn>0.0){
    find_chi_inverse(np[lA].bs->X,chiinv);
    qnew=0.0;
    for (i=0; i<nd; i++){
      Vnew[i]=0.0e0;
      for (j=0; j<nd; j++){
        Vnew[i]=Vnew[i]+chiinv[i][j]*Vstarnew[j];
      }
      qnew+=sqr(Vnew[i]);
    }
    qnew=sqrt(qnew);
  } else {
    /* set the new flow speed at the inflow to zero if the fluid is flowing towards the boundary */
    qnew=0.0;
    for (i=0; i<nd; i++){
      Vnew[i]=0.0;
    }
  }
  gamma=gl->model.fluid.gamma;
  Rgas=gl->model.fluid.R;
    
  Mold=_q(np[lA])/_a(np[lA],gl);
  Tstag=_T(np[lA],gl)*(1.0+(gamma-1.0)/2.0*sqr(Mold));
  Pstag=_P(np[lA],gl)*pow(1.0+(gamma-1.0)/2.0*sqr(Mold),gamma/(gamma-1.0));

  Tnew=Tstag-(gamma-1.0)/2.0*sqr(qnew)/(gamma*Rgas);
  Mnew=qnew/sqrt(gamma*Rgas*Tnew);
  if (Mnew<1.0) {
    /* only update the inflow boundary if the new Mach number is expected to be less than 1 */
    Pnew=Pstag/pow(1.0+(gamma-1.0)/2.0*sqr(Mnew),gamma/(gamma-1.0));
    assert_np(np[lA],is_node_resumed(np[lA]));
    find_U_2(np, lA, gl, Vnew, Pnew, Tnew);
  }
  
}


/* subsonic outflow BC with constant P */
static void update_bdry_subsonic_outflow(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                         bool BDRYDIRECFOUND, int ACCURACY){
  dim_t Vnew,Vstarnew;
  long dim,i,j;
  double Tnew,Pnew,qnew;
  dim2_t chiinv;

  for (dim=0; dim<nd; dim++){
    Vstarnew[dim]=_f_symmetry(ACCURACY,_Vstar(np[lB],dim),_Vstar(np[lC],dim));
  }
  Tnew=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
  if (Vstarnew[theta]*thetasgn<0.0){
    find_chi_inverse(np[lA].bs->X,chiinv);
    qnew=0.0;
    for (i=0; i<nd; i++){
      Vnew[i]=0.0e0;
      for (j=0; j<nd; j++){
        Vnew[i]=Vnew[i]+chiinv[i][j]*Vstarnew[j];
      }
      qnew+=sqr(Vnew[i]);
    }
    qnew=sqrt(qnew);
  } else {
    /* set the new flow speed at the outflow to zero if the fluid is flowing away from the boundary */
    qnew=0.0;
    for (i=0; i<nd; i++){
      Vnew[i]=0.0;
    }
  }
    
  Pnew=_P(np[lA],gl);

  assert_np(np[lA],is_node_resumed(np[lA]));
  find_U_2(np, lA, gl, Vnew, Pnew, Tnew);
}


static void update_bdry_freestream(np_t *np, gl_t *gl, long lA, long lB, long lC, long theta, long thetasgn,
                                   bool BDRYDIRECFOUND, int ACCURACY){
  dim_t V;
  double T,P;
  bool OUTFLOW, FREESTREAM;

  find_V_P_T_bdry_freestream(np, gl, lA, lB, lC, theta, thetasgn, ACCURACY, V, &P, &T, &OUTFLOW, &FREESTREAM);
  find_U_2(np, lA, gl, V, P, T);
}


static void update_bdry_wall(np_t *np, gl_t *gl, long lA, long lB, long lC,
                            long theta, long thetasgn,
                            bool ADIABATIC, bool BDRYDIRECFOUND, int ACCURACY){
  double Twall,Pwall;
  long dim;
  dim_t Vwall;

  if (ADIABATIC) {
    Twall=_f_symmetry(ACCURACY,_T(np[lB],gl),_T(np[lC],gl));
  } else {
    Twall=_bdry_param(np,gl,lA,0,TYPELEVEL_FLUID_WORK);
  }
  /* clip Twall so that it remains within the user-specified bounds */
  Twall=max(gl->model.fluid.Tmin,min(gl->model.fluid.Tmax,Twall));
  
  for (dim=0; dim<nd; dim++) Vwall[dim]=0.0e0;

  find_Pstar_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY, &Pwall);

  find_U_2(np, lA,gl,Vwall,Pwall,Twall);
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

    case BDRY_SYMMETRICAL1 :
      if (BDRYDIRECFOUND)
         update_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_SYMMETRICAL2:
      if (BDRYDIRECFOUND)
        update_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_SECONDORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_SYMMETRICAL3:
      if (BDRYDIRECFOUND)
        update_bdry_symmetrical(np, gl, lA, lB, lC, lD, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_THIRDORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLTFIXED1:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLADIABATIC1:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, TRUE, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_WALLTFIXED2:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, FALSE, BDRYDIRECFOUND, ACCURACY_SECONDORDER);
    break;

    case BDRY_WALLADIABATIC2:
      update_bdry_wall(np, gl, lA, lB, lC, theta, thetasgn, TRUE, BDRYDIRECFOUND, ACCURACY_SECONDORDER);
    break;

    case BDRY_INFLOWSUBSONIC1:
      if (BDRYDIRECFOUND)
        update_bdry_subsonic_inflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
    break;

    case BDRY_OUTFLOWSUBSONIC1:
      if (BDRYDIRECFOUND)
        update_bdry_subsonic_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
      else update_bdry_outflow(np, gl, lA, lB, lC, theta, thetasgn, BDRYDIRECFOUND, ACCURACY_FIRSTORDER);
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


