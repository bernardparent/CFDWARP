// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2001 Giovanni Fusina

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

#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/.active/model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>

#define maxreact 30



void write_chem_template(FILE **controlfile){
}

void read_and_init_chem_actions(char *action, char **argum, SOAP_codex_t *codex){

}

void find_Schem(np_t np, gl_t *gl, flux_t S){
  long flux;
  long k,r,s;
  double conc[ns];
  double T,w[ns],rho,frr[maxreact],Kc[maxreact];
  double Gs[ns],dGs[maxreact];
  double third;

/*  FindChemData; */


T=_T(np);

for (s=0; s<ns; s++){
  w[s]=_w(np,s);
}

rho=_rho(np);



/*Find concentrations: mol/cm3 and Gibbs: J/mol */

  for (k=0; k<ns; k++){
        Gs[k]=0.0e0;
  }

  for (k=0; k<ns; k++){
    conc[k] = w[k] * rho / _calM(k) * 1.0e-06;
    Gs[k]=(_hk_from_T(k,T)-T*_sk_from_T(k,T))*_calM(k);
  }

        for (r=0; r<maxreact; r++){
          frr[r]=0.0e0;
          dGs[r]=0.0e0;
        }



/*Difference in Gibbs for each reaction */



        dGs[1]=2.0e0*Gs[4]-(Gs[0]+Gs[1]);
        dGs[2]=(Gs[4]+Gs[3])-(Gs[2]+Gs[1]);
        dGs[3]=(Gs[4]+Gs[2])-(Gs[3]+Gs[0]);
        dGs[4]=(Gs[5]+Gs[2])-(Gs[4]+Gs[0]);
        dGs[5]=(Gs[5]+Gs[3])-(2.0e0*Gs[4]);
        dGs[6]=Gs[5]-(Gs[2]+Gs[4]);
        dGs[7]=Gs[0]-2.0e0*Gs[2];


/*      Rgas = 8.3144126d0*/

/*determine Kp (storing results in Kc since this is ultimately 
         what is sought here).

       [Kp] = (atm) ** (n-m)*/


        for (r=0; r<maxreact; r++){
          Kc[r]=0.0e0;
        }

        for (r=1; r<maxreact; r++){
          Kc[r]=exp(-dGs[r]/(8.3144126e0*T));
        }



/*
*       determine Kc by multiplying Kp by (RT/101325*100**3) ** (m-n)
*
*       [Kc] = (cm**3/mol) ** (m-n)
*/


        Kc[6] = Kc[6] * (82.056872e0 * T);
        Kc[7] = Kc[7] * (82.056872e0 * T);




/*       CALCULATE FORWARD REACTION RATE
*       R = 1.987192004 cal/mol-K
*       [Kf] = (cm**3/mol) ** (m-1) / sec
*/




        frr[1]=1.7e+13*exp(-48150.0e0/(1.987192004e0*T));
        frr[2]=1.42e+14*exp(-16400.0e0/(1.987192004e0*T));
        frr[3]=2.07e+14*exp(-13750.0e0/(1.987192004e0*T));
        frr[4]=3.16e+07*pow(T,1.8e0)*exp(-3030.0e0/(1.987192004e0*T));
        frr[5]=5.50e+13*exp(-7000.0e0/(1.987192004e0*T));
        frr[6]=2.21e+22*pow(T,(-2.0e0));
        frr[7]=6.53e+17/T;



/*
*       calculate the reaction rates
*
*
*       calculate third body efficiencies
*/

        third = conc[0];
        for (k=1; k<ns; k++){ 
          third = third + conc[k];
        }

/*
*       find the species production rates [kg/s]
*/


  for (flux=0; flux<nf; flux++){
    S[flux]=0.0e0;
  }


        S[0]=(-frr[1]*(conc[0]*conc[1]-1.0e0/Kc[1]*(conc[4]*conc[4]))
           -frr[3]*(conc[3]*conc[0]-1.0e0/Kc[3]*(conc[4]*conc[2]))
           -frr[4]*(conc[4]*conc[0]-1.0e0/Kc[4]*(conc[5]*conc[2]))
           +frr[7]*(conc[2]*conc[2]-1.0e0/Kc[7]*(conc[0]))*third
           )*1.0e+06*_calM(0)*_Omega(np,gl);

       S[1]=(-frr[1]*(conc[0]*conc[1]-1.0e0/Kc[1]*(conc[4]*conc[4]))
           -frr[2]*(conc[2]*conc[1]-1.0e0/Kc[2]*(conc[4]*conc[3]))
           )*1.0e+06*_calM(1)*_Omega(np,gl);

       S[2]=(-frr[2]*(conc[2]*conc[1]-1.0e0/Kc[2]*(conc[4]*conc[3]))
           +frr[3]*(conc[3]*conc[0]-1.0e0/Kc[3]*(conc[4]*conc[2]))
           +frr[4]*(conc[4]*conc[0]-1.0e0/Kc[4]*(conc[5]*conc[2]))
           -frr[6]*(conc[2]*conc[4]-1.0e0/Kc[6]*(conc[5]))*third
         -2.0e0*frr[7]*(conc[2]*conc[2]-1.0e0/Kc[7]*(conc[0]))*third
           )*1.0e+06*_calM(2)*_Omega(np,gl);

       S[3]=(+frr[2]*(conc[2]*conc[1]-1.0e0/Kc[2]*(conc[4]*conc[3]))
           -frr[3]*(conc[3]*conc[0]-1.0e0/Kc[3]*(conc[4]*conc[2]))
           +frr[5]*(conc[4]*conc[4]-1.0e0/Kc[5]*(conc[5]*conc[3]))
           )*1.0e+06*_calM(3)*_Omega(np,gl);

       S[4]=(+2.0e0*frr[1]*(conc[0]*conc[1]-1.0e0/Kc[1]*(conc[4]*conc[4]))
             +frr[2]*(conc[2]*conc[1]-1.0e0/Kc[2]*(conc[4]*conc[3]))
             +frr[3]*(conc[3]*conc[0]-1.0e0/Kc[3]*(conc[4]*conc[2]))
             -frr[4]*(conc[4]*conc[0]-1.0e0/Kc[4]*(conc[5]*conc[2]))
           -2.0e0*frr[5]*(conc[4]*conc[4]-1.0e0/Kc[5]*(conc[5]*conc[3]))
             -frr[6]*(conc[2]*conc[4]-1.0e0/Kc[6]*(conc[5]))*third
             )*1.0e+06*_calM(4)*_Omega(np,gl);

       S[5]=(+frr[4]*(conc[4]*conc[0]-1.0e0/Kc[4]*(conc[5]*conc[2]))
        +frr[5]*(conc[4]*conc[4]-1.0e0/Kc[5]*(conc[5]*conc[3]))
        +frr[6]*(conc[2]*conc[4]-1.0e0/Kc[6]*(conc[5]))*third
        )*1.0e+06*_calM(5)*_Omega(np,gl);



}


void find_prim_fluid_mem_chem(double P, double T, np_t np){
  spec_t w;
  double sum1,rho;
  long dim,spec;
  bool ref_flag;

  assert(is_node_resumed(np));
  ref_flag=FALSE;
  np.wk->Tmem=T;
  /* np.wk->Pmem=P; */
  sum1=0.0e0;
  for (spec=0; spec<ns; spec++) sum1=sum1+np.bs->U[spec];
  assert(sum1!=0.0e0);
  for (spec=0; spec<ns; spec++) w[spec]=np.bs->U[spec]/sum1;
  rho=_rho_from_w_P_T(w,P,T);
  np.wk->rhomem=rho;
  assert(np.wk->rhomem>0.0e0);
  /* find Vmem */
  assert(rho!=0.0e0);
  for (dim=0; dim<nd; dim++) np.wk->Vmem[dim]=np.bs->U[ns+dim]/rho;
}

void find_prim_chem(np_t np){
  long dim,spec;
  double rho,T,P,e,sum,sum1;
  spec_t w;
  dim_t V;

  sum1=0.0e0;
  for (spec=0; spec<ns; spec++) sum1=sum1+np.bs->U[spec];
  assert(sum1!=0.0e0);
  for (spec=0; spec<ns; spec++) w[spec]=np.bs->U[spec]/sum1;
  rho=sum1;
  sum=0.0e0;
  for (dim=0; dim<nd; dim++){
    V[dim]=np.bs->U[ns+dim]/sum1;
    sum=sum+sqr(V[dim]);
  }
  assert(rho!=0.0e0);
  e=np.bs->U[ns+nd]/sum1-0.5e0*sum
   -np.bs->U[ns+nd+1]/sum1;
  T=_T_from_w_e(w,e);
  np.wk->Tmem=T;
  P=_P_from_w_rho_T(w,rho,T);

  find_prim_fluid_mem_chem(P, T, np);
}


void add_dU_to_Utilde_chem(flux_t dUstar, np_t np, gl_t *gl){
  double Omega;
  long flux;

  Omega=_Omega(np,gl);
  assert(Omega>0.0e0);
  for (flux=0; flux<nf; flux++){
    np.bs->U[flux]=np.bs->U[flux]+dUstar[flux]/Omega;
  }
  find_prim_chem(np);
}




void find_numerical_jacobian_chem(void(*FUNCT)(np_t, gl_t *, flux_t),
                       np_t np, gl_t *gl, sqmat_t Ak){
 long row,col,flux;
 flux_t F,Fnew,Ustar,dUstar,dUstar2;

 for (flux=0; flux<nf; flux++){
   Ustar[flux]=np.bs->U[flux]*_Omega(np,gl);
   dUstar[flux]=Ustar[flux]/100000000.0e0+1.0e-15;

/*   if(flux<ns){
     dUstar[flux]=1.0e-11*_Omega(np)*_rho(np);
   } else {
   if(flux<ns+nd){
     dUstar[flux]=1.0e-11*_Omega(np)*_rho(np)*_a(np);
   }  else {
        dUstar[flux]=Ustar[flux]/100000000.0e0+1.0e-12*_Omega(np);
      }
   }
*/
}
 (*FUNCT)(np,gl,F);
 for (col=0; col<nf; col++){
   for (flux=0; flux<nf; flux++){
     np.bs->U[flux]=Ustar[flux]/_Omega(np,gl);
     dUstar2[flux]=0.0e0;
   }
   dUstar2[col]=dUstar[col];
   add_dU_to_Utilde_chem(dUstar2,np,gl);
   find_prim_chem(np);
   (*FUNCT)(np, gl,Fnew);
   for (row=0; row<nf; row++){
     assert(dUstar[col]!=0.0e0);
     Ak[row][col]=(Fnew[row]-F[row])/dUstar[col];
   }
 }
 for (flux=0; flux<nf; flux++){
   np.bs->U[flux]=Ustar[flux]/_Omega(np,gl);
 }
 find_prim_chem(np);
}

/*---------------------------------------------------*/

void find_dSchem_dU(np_t np, gl_t *gl, sqmat_t C){

 find_numerical_jacobian_chem(&(find_Schem), np, gl, C);
 /* display_matrix(C);*/
}



