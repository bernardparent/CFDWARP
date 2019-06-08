// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2002 Timothy Hui


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

/*******************************************************************************
   Define the number of chemical reactions for the Lutz combustion model.
   There are 21 reactions here because C counts arrays from 0, and the 0 value
   is not used.  This is done for convenience and clarity so the reaction 1
   corresponds to reaction 1, instead of 0 to 1, 1 to 2, etc.
*******************************************************************************/

#define maxreact 21

/* Declare functions to access various properties: temperature, velocity,
   composition, and density */


void write_chem_template(FILE **controlfile){
}

void read_and_init_chem_actions(char *action, char **argum, SOAP_codex_t *codex){

}

void find_Schem(np_t np, gl_t *gl, flux_t S)
{
	long 	flux;
  	long 	k,r,s;
  	double 	X[ns];
  	double 	T,w[ns],rho,frr[maxreact],Kc[maxreact];
  	double 	Gs[ns],dGs[maxreact];
  	double 	third,third5,third13,third14;
	
	/* Get values for temperature, species composition, and density */

	T   = _T(np);
	rho = _rho(np);

	for ( s = 0; s < ns; s++ )
	{
	    w[s] = _w(np,s);
	}

	/* Calculate properties for each species */

	/* Initialize matrices to store values for Gibbs free energy [J/mol] */

  	for ( k = 0; k < ns; k++ )
	{
	    Gs[k] = 0.0e0;
  	}

	/* Initialize matrices to store values for forward reaction rate,
	   difference in Gibbs free energy */

	for ( k = 0; k < ns; k++ )
	{
	    X[k]  = w[k]*rho/_calM(k)*1.0e-06;
	    Gs[k] = (_hk_from_T(k,T) - T*_sk_from_T(k,T))*_calM(k);
	}

  	for ( r = 0; r < maxreact; r++ )
	{
	    frr[r] = 0.0e0;
	    dGs[r] = 0.0e0;
	}

	/* Calculate difference in Gibbs for each reaction */

	dGs[1]  = 2.0e0*Gs[4]       -	(Gs[0] + Gs[1]);
	dGs[2]  = (Gs[5] + Gs[2])   - 	(Gs[4] + Gs[0]);
	dGs[3]  = (Gs[1] + Gs[2])   - 	(Gs[3] + Gs[4]);
	dGs[4]  = (Gs[4] + Gs[2])   - 	(Gs[3] + Gs[0]);
	dGs[5]  = Gs[6]             - 	(Gs[2] + Gs[1]);
	dGs[6]  = (Gs[5] + Gs[1])   - 	(Gs[4] + Gs[6]);
	dGs[7]  = 2.0e0*Gs[4]       - 	(Gs[2] + Gs[6]);
	dGs[8]  = (Gs[1] + Gs[4])   - 	(Gs[3] + Gs[6]);
	dGs[9]  = (Gs[3] + Gs[5])   - 	2.0e0*Gs[4];
	dGs[10] = Gs[0]             - 	2.0e0*Gs[2];
	dGs[11] = 2.0e0*Gs[0]       - 	(2.0e0*Gs[2] + Gs[0]);
	dGs[12] = (Gs[0] + Gs[5])   - 	(2.0e0*Gs[2] + Gs[5]);
	dGs[13] = Gs[5]             - 	(Gs[2] + Gs[4]);
	dGs[14] = Gs[4]             - 	(Gs[2] + Gs[3]);
	dGs[15] = Gs[1]             - 	2.0e0*Gs[3];
	dGs[16] = (Gs[0] + Gs[1])   - 	(Gs[2] + Gs[6]);
	dGs[17] = (Gs[7] + Gs[1])   - 	2.0e0*Gs[6];
	dGs[18] = 2.0e0*Gs[4]       - 	Gs[7];
	dGs[19] = (Gs[6] + Gs[0])   - 	(Gs[7] + Gs[2]);
	dGs[20] = (Gs[5] + Gs[6])   - 	(Gs[7] + Gs[4]);

	/* Determine Kp (storing results in Kc since this is ultimately
	   what is sought here). [(atm) ** (n-m)] */
	
	for ( r = 0; r < maxreact; r++ )
	{	
	    Kc[r] = 0.0e0;
	}

	for ( r = 1; r < maxreact; r++ )
	{
	    Kc[r] = exp(-dGs[r]/(8.3144126e0*T));
	}

	/* Determine Kc by multiplying Kp by (RT/101325*100**3)**(m-n)
	  [(cm**3/mol) ** (m-n)] */

	Kc[5]  = Kc[5]  * (82.056872e0 * T);
   	Kc[10] = Kc[10] * (82.056872e0 * T);
   	Kc[11] = Kc[11] * (82.056872e0 * T);
   	Kc[12] = Kc[12] * (82.056872e0 * T);
   	Kc[13] = Kc[13] * (82.056872e0 * T);
   	Kc[14] = Kc[14] * (82.056872e0 * T);
	Kc[15] = Kc[15] * (82.056872e0 * T);
	Kc[18] = Kc[18] / (82.056872e0 * T);

	/* Calculate forward reaction rate, here a different R is used
	   where R = 1.987192004 cal/g-K. [(cm**3/mol)**(m-1)/s] */

	frr[1]  = 1.70e+13 * exp(-47780.0e0/(1.987192004e0*T));
	frr[2]  = 1.17e+09 * pow(T, (1.300e0))*exp(-3626.0e0/(1.987192004e0*T));
	frr[3]  = 4.00e+14 * pow(T, (-0.500e0));
	frr[4]  = 5.06e+04 * pow(T, (2.670e0))*exp(-6290.0e0/(1.987192004e0*T));
	frr[5]  = 3.61e+17 * pow(T, (-0.720e0));
	frr[6]  = 7.50e+12;
	frr[7]  = 1.40e+14 * exp(-1073.0e0/(1.987192004e0*T));
	frr[8]  = 1.40e+13 * exp(-1073.0e0/(1.987192004e0*T));
	frr[9]  = 6.00e+08 * pow(T, (1.300e0));
	frr[10] = 1.00e+18 * pow(T, (-1.000e0));
	frr[11] = 9.20e+16 * pow(T, (-0.600e0));
	frr[12] = 6.00e+19 * pow(T, (-1.250e0));
	frr[13] = 1.60e+22 * pow(T, (-2.000e0));
	frr[14] = 6.20e+16 * pow(T, (-0.600e0));
	frr[15] = 1.89e+13 * exp(1788.0e0/(1.987192004e0*T));
	frr[16] = 1.25e+13;
	frr[17] = 2.00e+12;
	frr[18] = 1.30e+17 * exp(-45500.0e0/(1.987192004e0*T));
	frr[19] = 1.60e+12 * exp(-3800.0e0/(1.987192004e0*T));
	frr[20] = 1.00e+13 * exp(-1800.0e0/(1.987192004e0*T));

	/* Calculate third body efficiencies */

	third = 0.0e0;
	for ( k = 0; k < ns; k++ )
	{
	    third = third + X[k];
	}

	third5	= third + 1.86e0 * X[0] + 17.60e0 * X[5] + 0.26e0 * X[8];
	third13	= third                 +  4.00e0 * X[5];
	third14 = third                 +  4.00e0 * X[5]; 	

	/* Calculate the species production rates [kg/s] */

	for ( flux = 0; flux < nf; flux++ )
	{
	    S[flux]=0.0e0;
  	}

	S[0] = (-frr[1]  * (X[0]*X[1]      - 1.0e0/Kc[1]*(X[4]*X[4]))
	        -frr[2]  * (X[4]*X[0]      - 1.0e0/Kc[2]*(X[5]*X[2]))
	        -frr[4]  * (X[3]*X[0]      - 1.0e0/Kc[4]*(X[4]*X[2]))
	        +frr[10] * (X[2]*X[2] 	   - 1.0e0/Kc[10]*(X[0]))*third
	        +frr[11] * (X[2]*X[2]*X[0] - 1.0e0/Kc[11]*(X[0]*X[0]))
	        +frr[12] * (X[2]*X[2]*X[5] - 1.0e0/Kc[12]*(X[0]*X[5]))
	        +frr[16] * (X[2]*X[6]      - 1.0e0/Kc[16]*(X[0]*X[1]))
	        +frr[19] * (X[7]*X[2]      - 1.0e0/Kc[19]*(X[6]*X[0]))
	       )*1.0e+06*_calM(0)*_Omega(np,gl);

  	S[1] = (-frr[1]  * (X[0]*X[1] - 1.0e0/Kc[1]*(X[4]*X[4]))
	        +frr[3]  * (X[3]*X[4] - 1.0e0/Kc[3]*(X[1]*X[2]))
	        -frr[5]  * (X[2]*X[1] - 1.0e0/Kc[5]*(X[6]))*third5
	        +frr[6]  * (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	        +frr[8]  * (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	        +frr[15] * (X[3]*X[3] - 1.0e0/Kc[15]*(X[1]))*third
	        +frr[16] * (X[2]*X[6] - 1.0e0/Kc[16]*(X[0]*X[1]))
	        +frr[17] * (X[6]*X[6] - 1.0e0/Kc[17]*(X[7]*X[1]))
	       )*1.0e+06*_calM(1)*_Omega(np,gl);

 	S[2] = (+frr[2]         * (X[4]*X[0] 	  - 1.0e0/Kc[2]*(X[5]*X[2]))
	        +frr[3]         * (X[3]*X[4] 	  - 1.0e0/Kc[3]*(X[1]*X[2]))
	        +frr[4]         * (X[3]*X[0] 	  - 1.0e0/Kc[4]*(X[4]*X[2]))
	        -frr[5]         * (X[2]*X[1] 	  - 1.0e0/Kc[5]*(X[6]))*third5
	        -frr[7]         * (X[2]*X[6] 	  - 1.0e0/Kc[7]*(X[4]*X[4]))
	        -2.0e0*frr[10]  * (X[2]*X[2] 	  - 1.0e0/Kc[10]*(X[0]))*third
	        -2.0e0*frr[11]  * (X[2]*X[2]*X[0] - 1.0e0/Kc[11]*(X[0]*X[0]))
	        -2.0e0*frr[12]  * (X[2]*X[2]*X[5] - 1.0e0/Kc[12]*(X[0]*X[5]))
	        -frr[13]        * (X[2]*X[4]      - 1.0e0/Kc[13]*(X[5]))*third13
	        -frr[14]        * (X[2]*X[3]      - 1.0e0/Kc[14]*(X[4]))*third14
	        -frr[16]        * (X[2]*X[6]      - 1.0e0/Kc[16]*(X[0]*X[1]))
	        -frr[19]        * (X[7]*X[2]      - 1.0e0/Kc[19]*(X[6]*X[0]))
	       )*1.0e+06*_calM(2)*_Omega(np,gl);

  	S[3] = (-frr[3]         * (X[3]*X[4] - 1.0e0/Kc[3]*(X[1]*X[2]))
	        -frr[4]         * (X[3]*X[0] - 1.0e0/Kc[4]*(X[4]*X[2]))
	        -frr[8]         * (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	        +frr[9]         * (X[4]*X[4] - 1.0e0/Kc[9]*(X[3]*X[5]))
	        -frr[14]        * (X[2]*X[3] - 1.0e0/Kc[14]*(X[4]))*third14
	        -2.0e0*frr[15]  * (X[3]*X[3] - 1.0e0/Kc[15]*(X[1]))*third
	       )*1.0e+06*_calM(3)*_Omega(np,gl);

  	S[4] = (+2.0e0*frr[1]	* (X[0]*X[1] - 1.0e0/Kc[1]*(X[4]*X[4]))
	        -frr[2]       	* (X[4]*X[0] - 1.0e0/Kc[2]*(X[5]*X[2]))
	        -frr[3]       	* (X[3]*X[4] - 1.0e0/Kc[3]*(X[1]*X[2]))
	        +frr[4]         * (X[3]*X[0] - 1.0e0/Kc[4]*(X[4]*X[2]))
	        -frr[6]         * (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	        +2.0e0*frr[7]  	* (X[2]*X[6] - 1.0e0/Kc[7]*(X[4]*X[4]))
	        +frr[8]       	* (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	        -2.0e0*frr[9]  	* (X[4]*X[4] - 1.0e0/Kc[9]*(X[3]*X[5]))
	        -frr[13]      	* (X[2]*X[4] - 1.0e0/Kc[13]*(X[5]))*third13
	        +frr[14]      	* (X[2]*X[3] - 1.0e0/Kc[14]*(X[4]))*third14
	        +2.0e0*frr[18] 	* (X[7]      - 1.0e0/Kc[18]*(X[4]*X[4]))*third
	        -frr[20]      	* (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	       )*1.0e+06*_calM(4)*_Omega(np,gl);

  	S[5] = (+frr[2]		* (X[4]*X[0] - 1.0e0/Kc[2]*(X[5]*X[2]))
	        +frr[6]    	* (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	        +frr[9]    	* (X[4]*X[4] - 1.0e0/Kc[9]*(X[3]*X[5]))
	        +frr[13]   	* (X[2]*X[4] - 1.0e0/Kc[13]*X[5])*third13
	        +frr[20] 	* (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	       )*1.0e+06*_calM(5)*_Omega(np,gl);

  	S[6] = (+frr[5]       	* (X[2]*X[1] - 1.0e0/Kc[5]*(X[6]))*third5
	        -frr[6]       	* (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	        -frr[7]       	* (X[2]*X[6] - 1.0e0/Kc[7]*(X[4]*X[4]))
	        -frr[8]       	* (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	        -frr[16] 		* (X[2]*X[6] - 1.0e0/Kc[16]*(X[0]*X[1]))
	        -2.0e0*frr[17]	* (X[6]*X[6] - 1.0e0/Kc[17]*(X[7]*X[1]))
	        +frr[19]       	* (X[7]*X[2] - 1.0e0/Kc[19]*(X[6]*X[0]))
	        +frr[20]      	* (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	       )*1.0e+06*_calM(6)*_Omega(np,gl);
	
  	S[7] = (+frr[17]       * (X[6]*X[6] - 1.0e0/Kc[17]*(X[7]*X[1]))
	        -frr[18]       * (X[7] - 1.0e0/Kc[18]*(X[4]*X[4]))*third
	        -frr[19]       * (X[7]*X[2] - 1.0e0/Kc[19]*(X[6]*X[0]))
	        -frr[20]       * (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	       )*1.0e+06*_calM(7)*_Omega(np,gl);

}

void find_prim_fluid_mem_chem(double P, double T, np_t np){
  spec_t w;
  double sum1,rho;
  long dim,spec;
  bool ref_flag;

  assert(_ResumedNode(np));
  ref_flag=FALSE;
  np.wk->Tmem=T;
  /* np.wk->Pmem=P;*/
  sum1=0.0e0;
  for (spec=0; spec<ns; spec++) sum1=sum1+np.bs->U[spec];
  assert(sum1!=0.0e0);
  for (spec=0; spec<ns; spec++) w[spec]=np.bs->U[spec]/sum1;
  rho=_rho_from_w_P_T(w,P,T);
  np.wk->rhomem=rho;
  assert(np.wk->rhomem>0.0e0);
   /*find Vmem*/
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
   dUstar[flux]=Ustar[flux]/100000.0e0+1.0e-15;
 }
 (*FUNCT)(np,gl,F);
 for (col=0; col<nf; col++){
   for (flux=0; flux<nf; flux++){
     np.bs->U[flux]=Ustar[flux]/_Omega(np,gl);
     dUstar2[flux]=0.0e0;
   }
   dUstar2[col]=dUstar[col];
   add_dUstar_to_U(np,gl,dUstar2);
   find_prim_fluid(np,gl);
   (*FUNCT)(np, gl,Fnew);
   for (row=0; row<nf; row++){
     assert(dUstar[col]!=0.0e0);
     Ak[row][col]=(Fnew[row]-F[row])/dUstar[col];
   }
 }
 for (flux=0; flux<nf; flux++){
   np.bs->U[flux]=Ustar[flux]/_Omega(np,gl);
 }
 find_prim_fluid(np,gl);
}

/***********************************************************************
 This function calculates the analytical Jacobian for the chemical
 source term.
************************************************************************/

void find_dSchem_dU(np_t np, gl_t *gl, sqmat_t C)
{
	long 	k,r,s;  			/* counters */
  	long 	row,col;
  	double 	T,rho,frr[maxreact],Kc[maxreact];
  	double 	dGs[maxreact];
  	double 	third,third5,third13,third14;
  	double 	Cv,Rgas;
	double	deriv5,deriv10,deriv13,deriv14,deriv15,deriv18;
  	double 	dGsdT[ns],dXdrho[ns];
  	double 	dWsdX[ns][ns];
 	double 	dWsdT[ns],dWsdT1[ns];
  	double 	dfrrdT[maxreact];
  	double 	dKcdT[maxreact];
  	double 	dTdrhoE;
  	spec_t 	Cvs,dTdrhok,X,w,Gs;
  	dim_t 	dTdrhoV;
  	double 	hmix,CPmix,Rmix;

	Rgas = 8.3144126e0;

	/* Get values for temperature, species composition, and density */
	
	T = _T(np);
	for ( s = 0; s < ns; s++ )
	{
	    w[s] = _w(np,s);
	}
	rho = _rho(np);

	/* Calculate properties for each species */

	/* Initialize matrices to store values for Gibbs free energy [J/mol]
	   and specific heat at constant volume [J/g-K] and dWi/dT */

  	for ( k = 0; k < ns; k++ )
	{
	    Gs[k]     = 0.0e0;
	    Cvs[k]    = 0.0e0;
	    dWsdT1[k] = 0.0e0;
	    dWsdT[k]  = 0.0e0;
	    dGsdT[k]  = 0.0e0;
	}

	/* Initialize matrices to store values for forward reaction rate
	   difference in Gibbs free energy, rate of change of forward
	   reaction rate wrt temperature, and rate of change of equilibrium
	   constant wrt temperature. */

	for ( r = 0; r < maxreact; r++ )
	{
	    frr[r]    = 0.0e0;
	    dGs[r]    = 0.0e0;
	    dfrrdT[r] = 0.0e0;
	    dKcdT[r]  = 0.0e0;
	}

  	/* Calculate concentrations [mol/cm**3], Gibbs free energy [J/mol],
	   and specific heat at constant volume [J/g-K] */

	for ( k = 0; k < ns; k++ )

	{
	    X[k]   = w[k]*rho/_calM(k)*1.0e-06;
	    Gs[k]  = (_hk_from_T(k,T) - T * _sk_from_T(k,T))*_calM(k);
	    Cvs[k] = _cpk_from_T(k,T) - Rgas/_calM(k);
  	}

	/* Calculate properties for the entire mixture */	

	/* Calculate specific heat at constant volume [J/g-K] */
	
	Cv = 0.0e0;
	for ( k = 0; k < ns; k++ )
	{
	    Cv = Cv + w[k]*Cvs[k];
	}

	/* Calculate enthalpy [J/mol], specific heat at constant pressure
	   [J/g-K], and gas constant [J/g-K] */

	hmix	= 0.0e0;
	CPmix 	= 0.0e0;
	Rmix	= 0.0e0;
	for( k = 0; k < ns; k++ )
	{
	    hmix  = hmix  + w[k]*_hk_from_T(k,T);
	    CPmix = CPmix + w[k]*_cpk_from_T(k,T);
	    Rmix  = Rmix  + w[k]*(Rgas/_calM(k));
	}

	/* Calculate rate of change of Gibbs free energy wrt temperature */

	for ( k = 0; k < ns; k++ )
	{
	    dGsdT[k] = (Cvs[k] - _sk_from_T(k,T) + T*_dsk_dT_from_T(k,T)
	                - Rgas/_calM(k)
	               )*_calM(k);
  	}

	/* Calculate difference in Gibbs for each reaction */

	dGs[1]  = 2.0e0*Gs[4]       -	(Gs[0] + Gs[1]);
	dGs[2]  = (Gs[5] + Gs[2])   - 	(Gs[4] + Gs[0]);
	dGs[3]  = (Gs[1] + Gs[2])   - 	(Gs[3] + Gs[4]);
	dGs[4]  = (Gs[4] + Gs[2])   - 	(Gs[3] + Gs[0]);
	dGs[5]  = Gs[6]             - 	(Gs[2] + Gs[1]);
	dGs[6]  = (Gs[5] + Gs[1])   - 	(Gs[4] + Gs[6]);
	dGs[7]  = 2.0e0*Gs[4]       - 	(Gs[2] + Gs[6]);
	dGs[8]  = (Gs[1] + Gs[4])   - 	(Gs[3] + Gs[6]);
	dGs[9]  = (Gs[3] + Gs[5])   - 	2.0e0*Gs[4];
	dGs[10] = Gs[0]             - 	2.0e0*Gs[2];
	dGs[11] = 2.0e0*Gs[0]       - 	(2.0e0*Gs[2] + Gs[0]);
	dGs[12] = (Gs[0] + Gs[5])   - 	(2.0e0*Gs[2] + Gs[5]);
	dGs[13] = Gs[5]             - 	(Gs[2] + Gs[4]);
	dGs[14] = Gs[4]             - 	(Gs[2] + Gs[3]);
	dGs[15] = Gs[1]             - 	2.0e0*Gs[3];
	dGs[16] = (Gs[0] + Gs[1])   - 	(Gs[2] + Gs[6]);
	dGs[17] = (Gs[7] + Gs[1])   - 	2.0e0*Gs[6];
	dGs[18] = 2.0e0*Gs[4]       - 	Gs[7];
	dGs[19] = (Gs[6] + Gs[0])   - 	(Gs[7] + Gs[2]);
	dGs[20] = (Gs[5] + Gs[6])   - 	(Gs[7] + Gs[4]);

	/* Determine Kp (storing results in Kc since this is ultimately
	   what is sought here). [(atm) ** (n-m)] */

	for ( r = 0; r < maxreact; r++ )
	{
	    Kc[r] = 0.0e0;
    }

   	for ( r = 1; r < maxreact; r++ )
	{
	    Kc[r] = exp(-dGs[r]/(8.3144126e0*T));
   	}

	/* Determine Kc by multiplying Kp by (RT/101325*100**3)**(m-n)
	  [(cm**3/mol) ** (m-n)] */

	Kc[5]  = Kc[5]  * (82.056872e0 * T);
   	Kc[10] = Kc[10] * (82.056872e0 * T);
   	Kc[11] = Kc[11] * (82.056872e0 * T);
   	Kc[12] = Kc[12] * (82.056872e0 * T);
   	Kc[13] = Kc[13] * (82.056872e0 * T);
   	Kc[14] = Kc[14] * (82.056872e0 * T);
	Kc[15] = Kc[15] * (82.056872e0 * T);
	Kc[18] = Kc[18] / (82.056872e0 * T);

	
	/* Calculate forward reaction rate, here a different R is used
	   where R = 1.987192004 cal/g-K. [(cm**3/mol)**(m-1)/s] */

	frr[1]  = 1.70e+13 * exp(-47780.0e0/(1.987192004e0*T));
	frr[2]  = 1.17e+09 * pow(T, (1.300e0))*exp(-3626.0e0/(1.987192004e0*T));
	frr[3]  = 4.00e+14 * pow(T, (-0.500e0));
	frr[4]  = 5.06e+04 * pow(T, (2.670e0))*exp(-6290.0e0/(1.987192004e0*T));
	frr[5]  = 3.61e+17 * pow(T, (-0.720e0));
	frr[6]  = 7.50e+12;
	frr[7]  = 1.40e+14 * exp(-1073.0e0/(1.987192004e0*T));
	frr[8]  = 1.40e+13 * exp(-1073.0e0/(1.987192004e0*T));
	frr[9]  = 6.00e+08 * pow(T, (1.300e0));
	frr[10] = 1.00e+18 * pow(T, (-1.000e0));
	frr[11] = 9.20e+16 * pow(T, (-0.600e0));
	frr[12] = 6.00e+19 * pow(T, (-1.250e0));
	frr[13] = 1.60e+22 * pow(T, (-2.000e0));
	frr[14] = 6.20e+16 * pow(T, (-0.600e0));
	frr[15] = 1.89e+13 * exp(1788.0e0/(1.987192004e0*T));
	frr[16] = 1.25e+13;
	frr[17] = 2.00e+12;
	frr[18] = 1.30e+17 * exp(-45500.0e0/(1.987192004e0*T));
	frr[19] = 1.60e+12 * exp(-3800.0e0/(1.987192004e0*T));
	frr[20] = 1.00e+13 * exp(-1800.0e0/(1.987192004e0*T));

	/* Calculate third body efficiencies */

  	third = 0.0e0;
   	for ( k = 0; k < ns; k++ )
	{
	    third = third + X[k];
   	}
	third5	= third + 1.86e0 * X[0] + 17.6e0 * X[5] + 0.26e0 * X[8];
	third13	= third 				+  4.0e0 * X[5];
	third14 = third 				+  4.0e0 * X[5]; 	

/******************************************************************************/
	/* The following portion of this function will calculate the
	   chemical Jacobian C = dW/dUstar.  This derivative cannot be
	   obtained directly; the chain rule will be used and several
	   other derivatives will have to be determined. */

	/* For all the terms in the Jacobian, the derivative dWi/dT is
	   required; this has to be further broken down into more
	   derivatives.

	   dWi/dT = sum j=1->maxreact (  dWi/dfrrj*dfrrj/dT | Kcj=const
	                              + dWi/dKcj*dKcj/dT | frrj=const  )

	*/   	
/******************************************************************************/

/*-----------------------------dfrrj/dT------------------------------------*/
	/* First calculate dfrrj/dT */

	dfrrdT[1]  = frr[1]*(47780.0e0/(1.987192004e0*T*T));
   	dfrrdT[2]  = frr[2]*(3626.0e0/(1.987192004e0*T*T) + 1.300e0/T);
	dfrrdT[3]  = frr[3]*(-0.500e0/T);
	dfrrdT[4]  = frr[4]*(6290.0e0/(1.987192004e0*T*T) + 2.670e0/T);
   	dfrrdT[5]  = frr[5]*(-0.720e0/T);
   	dfrrdT[6]  = 0.0e0;
   	dfrrdT[7]  = frr[7]*(1073.0e0/(1.987192004e0*T*T));
   	dfrrdT[8]  = frr[8]*(1073.0e0/(1.987192004e0*T*T));
   	dfrrdT[9]  = frr[9]*(1.300e0/T);
   	dfrrdT[10] = frr[10]*(-1.000e0/T);
   	dfrrdT[11] = frr[11]*(-0.600e0/T);
   	dfrrdT[12] = frr[12]*(-1.250e0/T);
   	dfrrdT[13] = frr[13]*(-2.000e0/T);
   	dfrrdT[14] = frr[14]*(-0.600e0/T);
   	dfrrdT[15] = frr[15]*(-1788.0e0/(1.987192004e0*T*T));
   	dfrrdT[16] = 0.0e0;
   	dfrrdT[17] = 0.0e0;
   	dfrrdT[18] = frr[18]*(45500.0e0/(1.987192004e0*T*T));
   	dfrrdT[19] = frr[19]*(3800.0e0/(1.987192004e0*T*T));
   	dfrrdT[20] = frr[20]*(1800.0e0/(1.987192004e0*T*T));

/*-----------------------------dWi/dT------------------------------------*/	
	
	/* Calculate dWi/dfrrj*dfrrj/dT and
	   add dWi/dKcj*dKcj/dT to the result to get
	   dWi/dT*/
	
	dWsdT1[0] = (-dfrrdT[1]  * (X[0]*X[1]      - 1.0e0/Kc[1]*(X[4]*X[4]))
	             -dfrrdT[2]  * (X[4]*X[0]      - 1.0e0/Kc[2]*(X[5]*X[2]))
	             -dfrrdT[4]  * (X[3]*X[0]      - 1.0e0/Kc[4]*(X[4]*X[2]))
	             +dfrrdT[10] * (X[2]*X[2]      - 1.0e0/Kc[10]*(X[0]))*third
	             +dfrrdT[11] * (X[2]*X[2]*X[0] - 1.0e0/Kc[11]*(X[0]*X[0]))
	             +dfrrdT[12] * (X[2]*X[2]*X[5] - 1.0e0/Kc[12]*(X[0]*X[5]))
	             +dfrrdT[16] * (X[2]*X[6]      - 1.0e0/Kc[16]*(X[0]*X[1]))
	             +dfrrdT[19] * (X[7]*X[2]      - 1.0e0/Kc[19]*(X[6]*X[0]))
	            )*1.0e+6*_calM(0);

	dWsdT1[1] = (-dfrrdT[1]  * (X[0]*X[1] - 1.0e0/Kc[1]*(X[4]*X[4]))
	             +dfrrdT[3]  * (X[3]*X[4] - 1.0e0/Kc[3]*(X[1]*X[2]))
	             -dfrrdT[5]  * (X[2]*X[1] - 1.0e0/Kc[5]*(X[6]))*third5
	             +dfrrdT[6]  * (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	             +dfrrdT[8]  * (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	             +dfrrdT[15] * (X[3]*X[3] - 1.0e0/Kc[15]*(X[1]))*third
	             +dfrrdT[16] * (X[2]*X[6] - 1.0e0/Kc[16]*(X[0]*X[1]))
	             +dfrrdT[17] * (X[6]*X[6] - 1.0e0/Kc[17]*(X[7]*X[1]))
	            )*1.0e+6*_calM(1);

	dWsdT1[2] = (+dfrrdT[2]        * (X[4]*X[0]      - 1.0e0/Kc[2]*(X[5]*X[2]))
	             +dfrrdT[3]        * (X[3]*X[4]      - 1.0e0/Kc[3]*(X[1]*X[2]))
	             +dfrrdT[4]        * (X[3]*X[0]      - 1.0e0/Kc[4]*(X[4]*X[2]))
	             -dfrrdT[5]        * (X[2]*X[1]      - 1.0e0/Kc[5]*(X[6]))*third5
	             -dfrrdT[7]        * (X[2]*X[6]      - 1.0e0/Kc[7]*(X[4]*X[4]))
	             -2.0e0*dfrrdT[10] * (X[2]*X[2]      - 1.0e0/Kc[10]*(X[0]))*third
	             -2.0e0*dfrrdT[11] * (X[2]*X[2]*X[0] - 1.0e0/Kc[11]*(X[0]*X[0]))
	             -2.0e0*dfrrdT[12] * (X[2]*X[2]*X[5] - 1.0e0/Kc[12]*(X[0]*X[5]))
	             -dfrrdT[13]       * (X[2]*X[4]      - 1.0e0/Kc[13]*(X[5]))*third13
	             -dfrrdT[14]       * (X[2]*X[3]      - 1.0e0/Kc[14]*(X[4]))*third14
	             -dfrrdT[16]       * (X[2]*X[6]      - 1.0e0/Kc[16]*(X[0]*X[1]))
	             -dfrrdT[19]       * (X[7]*X[2]      - 1.0e0/Kc[19]*(X[6]*X[0]))
	            )*1.0e+6*_calM(2);

	dWsdT1[3] = (-dfrrdT[3]        * (X[3]*X[4] - 1.0e0/Kc[3]*(X[1]*X[2]))
	             -dfrrdT[4]        * (X[3]*X[0] - 1.0e0/Kc[4]*(X[4]*X[2]))
	             -dfrrdT[8]        * (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	             +dfrrdT[9]        * (X[4]*X[4] - 1.0e0/Kc[9]*(X[3]*X[5]))
	             -dfrrdT[14]       * (X[2]*X[3] - 1.0e0/Kc[14]*(X[4]))*third14
	             -2.0e0*dfrrdT[15] * (X[3]*X[3] - 1.0e0/Kc[15]*(X[1]))*third
	            )*1.0e+6*_calM(3);

	dWsdT1[4] = (+2.0e0*dfrrdT[1]  * (X[0]*X[1] - 1.0e0/Kc[1]*(X[4]*X[4]))
	             -dfrrdT[2]        * (X[4]*X[0] - 1.0e0/Kc[2]*(X[5]*X[2]))
	             -dfrrdT[3]        * (X[3]*X[4] - 1.0e0/Kc[3]*(X[1]*X[2]))
	             +dfrrdT[4]        * (X[3]*X[0] - 1.0e0/Kc[4]*(X[4]*X[2]))
	             -dfrrdT[6]        * (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	             +2.0e0*dfrrdT[7]  * (X[2]*X[6] - 1.0e0/Kc[7]*(X[4]*X[4]))
	             +dfrrdT[8]        * (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	             -2.0e0*dfrrdT[9]  * (X[4]*X[4] - 1.0e0/Kc[9]*(X[3]*X[5]))
	             -dfrrdT[13]       * (X[2]*X[4] - 1.0e0/Kc[13]*(X[5]))*third13
	             +dfrrdT[14]       * (X[2]*X[3] - 1.0e0/Kc[14]*(X[4]))*third14
	             +2.0e0*dfrrdT[18] * (X[7]      - 1.0e0/Kc[18]*(X[4]*X[4]))*third
	             -dfrrdT[20]       * (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	            )*1.0e+6*_calM(4);

	dWsdT1[5] = (+dfrrdT[2]  * (X[4]*X[0] - 1.0e0/Kc[2]*(X[5]*X[2]))
	             +dfrrdT[6]  * (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	             +dfrrdT[9]  * (X[4]*X[4] - 1.0e0/Kc[9]*(X[3]*X[5]))
	             +dfrrdT[13] * (X[2]*X[4] - 1.0e0/Kc[13]*X[5])*third13
	             +dfrrdT[20] * (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	            )*1.0e+6*_calM(5);

	dWsdT1[6] = (+dfrrdT[5]        * (X[2]*X[1] - 1.0e0/Kc[5]*(X[6]))*third5
	             -dfrrdT[6]        * (X[4]*X[6] - 1.0e0/Kc[6]*(X[5]*X[1]))
	             -dfrrdT[7]        * (X[2]*X[6] - 1.0e0/Kc[7]*(X[4]*X[4]))
	             -dfrrdT[8]        * (X[3]*X[6] - 1.0e0/Kc[8]*(X[1]*X[4]))
	             -dfrrdT[16]       * (X[2]*X[6] - 1.0e0/Kc[16]*(X[0]*X[1]))
	             -2.0e0*dfrrdT[17] * (X[6]*X[6] - 1.0e0/Kc[17]*(X[7]*X[1]))
	             +dfrrdT[19]       * (X[7]*X[2] - 1.0e0/Kc[19]*(X[6]*X[0]))
	             +dfrrdT[20]       * (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	            )*1.0e+6*_calM(6);

	dWsdT1[7] = (+dfrrdT[17] * (X[6]*X[6] - 1.0e0/Kc[17]*(X[7]*X[1]))
	             -dfrrdT[18] * (X[7]      - 1.0e0/Kc[18]*(X[4]*X[4]))*third
	             -dfrrdT[19] * (X[7]*X[2] - 1.0e0/Kc[19]*(X[6]*X[0]))
	             -dfrrdT[20] * (X[7]*X[4] - 1.0e0/Kc[20]*(X[5]*X[6]))
	            )*1.0e+6*_calM(7);

/*-----------------------------dKcj/dT------------------------------------*/
	/* First determine dKcj/dT, rate of change of the equilibrium
	   constant wrt temperature */
	
	dKcdT[1]  =  (dGs[1]/(Rgas*T*T))
	            -(2.0e0*dGsdT[4]        - (dGsdT[0] + dGsdT[1]))/(Rgas*T);

   	dKcdT[2]  =  (dGs[2]/(Rgas*T*T))
	            -((dGsdT[5] + dGsdT[2]) - (dGsdT[4] + dGsdT[0]))/(Rgas*T);

   	dKcdT[3]  =  (dGs[3]/(Rgas*T*T))
	            -((dGsdT[1] + dGsdT[2]) - (dGsdT[3] + dGsdT[4]))/(Rgas*T);

   	dKcdT[4]  =  (dGs[4]/(Rgas*T*T))
	            -((dGsdT[4] + dGsdT[2]) - (dGsdT[3] + dGsdT[0]))/(Rgas*T);

   	dKcdT[5]  =  (dGs[5]/(Rgas*T*T))
	            -(dGsdT[6]              - (dGsdT[2] + dGsdT[1]))/(Rgas*T)
	            +(1.0e0/T);

   	dKcdT[6]  =  (dGs[6]/(Rgas*T*T))
	            -((dGsdT[5] + dGsdT[1]) - (dGsdT[4] + dGsdT[6]))/(Rgas*T);

   	dKcdT[7]  =  (dGs[7]/(Rgas*T*T))
	            -(2.0e0*dGsdT[4]        - (dGsdT[2] + dGsdT[6]))/(Rgas*T);                	

   	dKcdT[8]  =  (dGs[8]/(Rgas*T*T))
	            -((dGsdT[1] + dGsdT[4]) - (dGsdT[3] + dGsdT[6]))/(Rgas*T);

   	dKcdT[9]  =  (dGs[9]/(Rgas*T*T))
	            -((dGsdT[3] + dGsdT[5]) - 2.0e0*dGsdT[4])/(Rgas*T);
   	
	dKcdT[10] =  (dGs[10]/(Rgas*T*T))
	            -(dGsdT[0]              - 2.0e0*dGsdT[2])/(Rgas*T)
	            +(1.0e0/T);

   	dKcdT[11] =  (dGs[11]/(Rgas*T*T))
	            -(2.0e0*dGsdT[0]        - (2.0e0*dGsdT[2] + dGsdT[0]))/(Rgas*T)
	            +(1.0e0/T);

   	dKcdT[12] =  (dGs[12]/(Rgas*T*T))
	            -((dGsdT[0] + dGsdT[5]) - (2.0e0*dGsdT[2] + dGsdT[5]))/(Rgas*T)
	            +(1.0e0/T);

   	dKcdT[13] =  (dGs[13]/(Rgas*T*T))
	            -(dGsdT[5]              - (dGsdT[2] + dGsdT[4]))/(Rgas*T)
	            +(1.0e0/T);

   	dKcdT[14] =  (dGs[14]/(Rgas*T*T))
	            -(dGsdT[4]              - (dGsdT[2] + dGsdT[3]))/(Rgas*T)
	            +(1.0e0/T);

   	dKcdT[15] =  (dGs[15]/(Rgas*T*T))
	            -(dGsdT[1]              - 2.0e0*dGsdT[3])/(Rgas*T)
	            +(1.0e0/T);

   	dKcdT[16] =  (dGs[16]/(Rgas*T*T))
	            -((dGsdT[0] + dGsdT[1]) - (dGsdT[2] + dGsdT[6]))/(Rgas*T);

  	dKcdT[17] =  (dGs[17]/(Rgas*T*T))
	            -((dGsdT[7] + dGsdT[1]) - 2.0e0*dGsdT[6])/(Rgas*T);

   	dKcdT[18] =  (dGs[18]/(Rgas*T*T))
	            -(2.0e0*dGsdT[4]        - dGsdT[7])/(Rgas*T)
	            -(1.0e0/T);

   	dKcdT[19] =  (dGs[19]/(Rgas*T*T))
	            -((dGsdT[6] + dGsdT[0]) - (dGsdT[7] + dGsdT[2]))/(Rgas*T);

   	dKcdT[20] =  (dGs[20]/(Rgas*T*T))
	            -((dGsdT[5] + dGsdT[6]) - (dGsdT[7] + dGsdT[4]))/(Rgas*T);

/*--------------------------dWi/dKcj*dKcj/dT---------------------------*/
	/* Calculate dWi/dKcj*dKcj/dT */

	dWsdT[0] = (-dKcdT[1]  * (frr[1]/Kc[1]*(X[4]*X[4]))
	            -dKcdT[2]  * (frr[2]/Kc[2]*(X[5]*X[2]))
	            -dKcdT[4]  * (frr[4]/Kc[4]*(X[4]*X[2]))
	            +dKcdT[10] * (frr[10]/Kc[10]*(X[0]))*third
	            +dKcdT[11] * (frr[11]/Kc[11]*(X[0]*X[0]))
	            +dKcdT[12] * (frr[12]/Kc[12]*(X[0]*X[5]))
	            +dKcdT[16] * (frr[16]/Kc[16]*(X[0]*X[1]))
	            +dKcdT[19] * (frr[19]/Kc[19]*(X[6]*X[0]))
	           )*1.0e+6*_calM(0) + dWsdT1[0];

 	dWsdT[1] = (-dKcdT[1]  * (frr[1]/Kc[1]*(X[4]*X[4]))
	            +dKcdT[3]  * (frr[3]/Kc[3]*(X[1]*X[2]))
	            -dKcdT[5]  * (frr[5]/Kc[5]*(X[6]))*third5
	            +dKcdT[6]  * (frr[6]/Kc[6]*(X[5]*X[1]))
	            +dKcdT[8]  * (frr[8]/Kc[8]*(X[1]*X[4]))
	            +dKcdT[15] * (frr[15]/Kc[15]*(X[1]))*third
	            +dKcdT[16] * (frr[16]/Kc[16]*(X[0]*X[1]))
	            +dKcdT[17] * (frr[17]/Kc[17]*(X[7]*X[1]))
	           )*1.0e+6*_calM(1) + dWsdT1[1];

  	dWsdT[2] = (+dKcdT[2]        * (frr[2]/Kc[2]*(X[5]*X[2]))
	            +dKcdT[3]        * (frr[3]/Kc[3]*(X[1]*X[2]))
	            +dKcdT[4]        * (frr[4]/Kc[4]*(X[4]*X[2]))
	            -dKcdT[5]        * (frr[5]/Kc[5]*(X[6]))*third5
	            -dKcdT[7]        * (frr[7]/Kc[7]*(X[4]*X[4]))
	            -2.0e0*dKcdT[10] * (frr[10]/Kc[10]*(X[0]))*third
	            -2.0e0*dKcdT[11] * (frr[11]/Kc[11]*(X[0]*X[0]))
	            -2.0e0*dKcdT[12] * (frr[12]/Kc[12]*(X[0]*X[5]))
	            -dKcdT[13]       * (frr[13]/Kc[13]*(X[5]))*third13
	            -dKcdT[14]       * (frr[14]/Kc[14]*(X[4]))*third14
	            -dKcdT[16]       * (frr[16]/Kc[16]*(X[0]*X[1]))
	            -dKcdT[19]       * (frr[19]/Kc[19]*(X[6]*X[0]))
	           )*1.0e+6*_calM(2) + dWsdT1[2];
	
	dWsdT[3] = (-dKcdT[3]        * (frr[3]/Kc[3]*(X[1]*X[2]))
	            -dKcdT[4]        * (frr[4]/Kc[4]*(X[4]*X[2]))
	            -dKcdT[8]        * (frr[8]/Kc[8]*(X[1]*X[4]))
	            +dKcdT[9]        * (frr[9]/Kc[9]*(X[3]*X[5]))
	            -dKcdT[14]       * (frr[14]/Kc[14]*(X[4]))*third14
	            -2.0e0*dKcdT[15] * (frr[15]/Kc[15]*(X[1]))*third
	           )*1.0e+6*_calM(3) + dWsdT1[3];

  	dWsdT[4] = (+2.0e0*dKcdT[1]	* (frr[1]/Kc[1]*(X[4]*X[4]))
	            -dKcdT[2]        * (frr[2]/Kc[2]*(X[5]*X[2]))
	            -dKcdT[3]        * (frr[3]/Kc[3]*(X[1]*X[2]))
	            +dKcdT[4]        * (frr[4]/Kc[4]*(X[4]*X[2]))
	            -dKcdT[6]        * (frr[6]/Kc[6]*(X[5]*X[1]))
	            +2.0e0*dKcdT[7]  * (frr[7]/Kc[7]*(X[4]*X[4]))
	            +dKcdT[8]        * (frr[8]/Kc[8]*(X[1]*X[4]))
	            -2.0e0*dKcdT[9]  * (frr[9]/Kc[9]*(X[3]*X[5]))
	            -dKcdT[13]       * (frr[13]/Kc[13]*(X[5]))*third13
	            +dKcdT[14]       * (frr[14]/Kc[14]*(X[4]))*third14
	            +2.0e0*dKcdT[18] * (frr[18]/Kc[18]*(X[4]*X[4]))*third
	            -dKcdT[20]       * (frr[20]/Kc[20]*(X[5]*X[6]))
	           )*1.0e+6*_calM(4) + dWsdT1[4];

  	dWsdT[5] = (+dKcdT[2]  * (frr[2]/Kc[2]*(X[5]*X[2]))
	            +dKcdT[6]  * (frr[6]/Kc[6]*(X[5]*X[1]))
	            +dKcdT[9]  * (frr[9]/Kc[9]*(X[3]*X[5]))
	            +dKcdT[13] * (frr[13]/Kc[13]*X[5])*third13
	            +dKcdT[20] * (frr[20]/Kc[20]*(X[5]*X[6]))
	           )*1.0e+6*_calM(5) + dWsdT1[5];

  	dWsdT[6] = (+dKcdT[5]        * (frr[5]/Kc[5]*(X[6]))*third5
	            -dKcdT[6]        * (frr[6]/Kc[6]*(X[5]*X[1]))
	            -dKcdT[7]        * (frr[7]/Kc[7]*(X[4]*X[4]))
	            -dKcdT[8]        * (frr[8]/Kc[8]*(X[1]*X[4]))
	            -dKcdT[16]       * (frr[16]/Kc[16]*(X[0]*X[1]))
	            -2.0e0*dKcdT[17] * (frr[17]/Kc[17]*(X[7]*X[1]))
	            +dKcdT[19]       * (frr[19]/Kc[19]*(X[6]*X[0]))
	            +dKcdT[20]       * (frr[20]/Kc[20]*(X[5]*X[6]))
	           )*1.0e+6*_calM(6) + dWsdT1[6];

   	dWsdT[7] = (+dKcdT[17] * (frr[17]/Kc[17]*(X[7]*X[1]))
	            -dKcdT[18] * (frr[18]/Kc[18]*(X[4]*X[4]))*third
	            -dKcdT[19] * (frr[19]/Kc[19]*(X[6]*X[0]))
	            -dKcdT[20] * (frr[20]/Kc[20]*(X[5]*X[6]))
	           )*1.0e+6*_calM(7) + dWsdT1[7];

/*------------------------------------------------------------------------*/
	/* Now that dWi/dT is found, let's determine the derivatives required
	   in order to calculate the Jacobian for the terms dWi/drhoj.  This
	   term is composed of the following partial derivatives:

	   dWi/drhoj = dWi/dXj*dXj/drhoj + dWi/dT*dT/drhoj

	   Let's first determine dWi/dXj */

/*-----------------------------dWi/dXj------------------------------------*/
	/* Initialize the matrix the store values for dWsdX */

	for ( row = 0; row < ns; row++ )
	{
	    for ( col = 0; col < ns; col++ )
	    {
	        dWsdX[row][col]=0.0e0;
	    }
	}
	
	deriv5	= frr[5]  * (X[2]*X[1] - 1.0e0/Kc[5]*X[6]);
	deriv10 = frr[10] * (X[2]*X[2] - 1.0e0/Kc[10]*X[0]);
	deriv13 = frr[13] * (X[2]*X[4] - 1.0e0/Kc[13]*X[5]);
	deriv14 = frr[14] * (X[2]*X[3] - 1.0e0/Kc[14]*X[4]);
	deriv15 = frr[15] * (X[3]*X[3] - 1.0e0/Kc[15]*X[1]);
	deriv18 = frr[18] * (X[7]      - 1.0e0/Kc[18]*X[4]*X[4]);

	dWsdX[0][0] = (- frr[1]*X[1] - frr[2]*X[4]   - frr[4]*X[3]
	               + deriv10     - frr[10]/Kc[10]*third + frr[11]*X[2]*X[2]
	               - 2.0e0*frr[11]/Kc[11]*X[0] - frr[12]/Kc[12]*X[5]
	               - frr[16]/Kc[16]*X[1] - frr[19]/Kc[19]*X[6]
	              )*1.0e+6*_calM(0);

   	dWsdX[0][1] = (- frr[1]*X[0] + deriv10 - frr[16]/Kc[16]*X[0]
	              )*1.0e+6*_calM(0);

  	dWsdX[0][2] = (+ frr[2]/Kc[2]*X[5] + frr[4]/Kc[4]*X[4] + deriv10
	               + 2.0e0*frr[10]*X[2]*third   + 2.0e0*frr[11]*X[2]*X[0]
	               + frr[16]*X[6] + frr[19]*X[7]
	              )*1.0e+6*_calM(0);

   	dWsdX[0][3] = (- frr[4]*X[0] + deriv10
	              )*1.0e+6*_calM(0);

   	dWsdX[0][4] = (2.0e0*frr[1]/Kc[1]*X[4] - frr[2]*X[0] + frr[4]/Kc[4]*X[2]
	               + deriv10
	              )*1.0e+6*_calM(0);

  	dWsdX[0][5] = (+ frr[2]/Kc[2]*X[2] + deriv10 + frr[12]*X[2]*X[2]
	               - frr[12]/Kc[12]*X[0]
	              )*1.0e+6*_calM(0);

   	dWsdX[0][6] = (+ deriv10 + frr[16]*X[2] - frr[19]/Kc[19]*X[0]
	              )*1.0e+6*_calM(0);

  	dWsdX[0][7] = (+ deriv10 + frr[19]*X[2]
	              )*1.0e+6*_calM(0);

   	dWsdX[0][8] = (+ deriv10
	              )*1.0e+6*_calM(0);

  	dWsdX[1][0] = (- frr[1]*X[1] - 2.86e0*deriv5 + deriv15
	               - frr[16]/Kc[16]*X[1]
	              )*1.0e+6*_calM(1);

  	dWsdX[1][1] = (- frr[1]*X[0]  - frr[3]/Kc[3]*X[2] - deriv5
	               - frr[5]*X[2]*third5 - frr[6]/Kc[6]*X[5]
	               - frr[8]/Kc[8]*X[4] + deriv15 - frr[15]/Kc[15]*third
	               - frr[16]/Kc[16]*X[0] - frr[17]/Kc[17]*X[7]
	              )*1.0e+6*_calM(1);

  	dWsdX[1][2] = (- frr[3]/Kc[3]*X[1] - deriv5 - frr[5]*X[1]*third5
	               + deriv15 + frr[16]*X[6]
	              )*1.0e+6*_calM(1);

   	dWsdX[1][3] = (+ frr[3]*X[4] - deriv5 + frr[8]*X[6] + deriv15
	               + 2.0e0*frr[15]*X[3]*third
	              )*1.0e+6*_calM(1);

 	dWsdX[1][4] = (+ 2.0e0*frr[1]/Kc[1]*X[4] + frr[3]*X[3] - deriv5
	               + frr[6]*X[6] - frr[8]/Kc[8]*X[1] + deriv15
	              )*1.0e+6*_calM(1);

   	dWsdX[1][5] = (- 18.6e0*deriv5 - frr[6]/Kc[6]*X[1] + deriv15
	              )*1.0e+6*_calM(1);

   	dWsdX[1][6] = (- deriv5  + frr[5]/Kc[5]*third5 + frr[6]*X[4]
	               + frr[8]*X[3] + deriv15 + frr[16]*X[2]
	               + 2.0e0*frr[17]*X[6]
	              )*1.0e+6*_calM(1);

   	dWsdX[1][7] = (- deriv5 + deriv15 - frr[17]/Kc[17]*X[1]
	              )*1.0e+6*_calM(1);

	dWsdX[1][8] = (- 1.26e0*deriv5 + deriv15
	              )*1.0e+6*_calM(1);

  	dWsdX[2][0] = (+ frr[2]*X[4] + frr[4]*X[3] - 2.86e0*deriv5
	               - 2.0e0*deriv10 + 2.0e0*frr[10]/Kc[10]*third
	               + 4.0e0*frr[11]/Kc[11]*X[0] + 2.0e0*frr[12]/Kc[12]*X[5]
	               - deriv13 - deriv14 + frr[16]/Kc[16]*X[1]
	               + frr[19]/Kc[19]*X[6]
	              )*1.0e+6*_calM(2);

	dWsdX[2][1] = (- frr[3]/Kc[3]*X[2] - deriv5 - frr[5]*X[2]*third5
	               - 2.0e0*deriv10 - deriv13 - deriv14 + frr[16]/Kc[16]*X[0]
	              )*1.0e+6*_calM(2);

   	dWsdX[2][2] = (- frr[2]/Kc[2]*X[5] - frr[3]/Kc[3]*X[1]
	               - frr[4]/Kc[4]*X[4] - deriv5 - frr[5]*X[1]*third5
	               - frr[7]*X[6] - 2.0e0*deriv10 - 4.0e0*frr[10]*X[2]*third
	               - 4.0e0*frr[11]*X[2]*X[0] - 4.0e0*frr[12]*X[2]*X[5]
	               - deriv13 - frr[13]*X[4]*third13 - deriv14
	               - frr[14]*X[3]*third14 - frr[16]*X[6] - frr[19]*X[7]
	              )*1.0e+6*_calM(2);

  	dWsdX[2][3] = (+ frr[3]*X[4] + frr[4]*X[0] - deriv5 - 2.0e0*deriv10
	               - deriv13 - deriv14 - frr[14]*X[2]*third14
	              )*1.0e+6*_calM(2);
   	
	dWsdX[2][4] = (+ frr[2]*X[0] + frr[3]*X[3] - frr[4]/Kc[4]*X[2] - deriv5
	               + 2.0e0*frr[7]/Kc[7]*X[4] - 2.0e0*deriv10 - deriv13
	               - frr[13]*X[2]*third13 - deriv14 + frr[14]/Kc[14]*third14
	              )*1.0e+6*_calM(2);


   	dWsdX[2][5] = (- frr[2]/Kc[2]*X[2] - 18.6e0*deriv5 - 2.0e0*deriv10
	               - 2.0e0*frr[12]*X[2]*X[2] + 2.0e0*frr[12]/Kc[12]*X[0]
	               - 5.0e0*deriv13 + frr[13]/Kc[13]*third13 - 5.0e0*deriv14
	              )*1.0e+6*_calM(2);

   	dWsdX[2][6] = (- deriv5 + frr[5]/Kc[5]*third5 - frr[7]*X[2]
	               - 2.0e0*deriv10 - deriv13 - deriv14 - frr[16]*X[2]
	               + frr[19]/Kc[19]*X[0]
	              )*1.0e+6*_calM(2);

   	dWsdX[2][7] = (- deriv5 - 2.0e0*deriv10 -deriv13 -deriv14 - frr[19]*X[2]
	              )*1.0e+6*_calM(2);

	dWsdX[2][8] = (- 1.26e0*deriv5 - 2.0e0*deriv10 -deriv13 - deriv14
	              )*1.0e+6*_calM(2);

  	dWsdX[3][0] = (- frr[4]*X[3] - deriv14 - 2.0e0*deriv15
	              )*1.0e+6*_calM(3);

   	dWsdX[3][1] = (+ frr[3]/Kc[3]*X[2] + frr[8]/Kc[8]*X[4] - deriv14
	               - 2.0e0*deriv15 + 2.0e0*frr[15]/Kc[15]*third
	              )*1.0e+6*_calM(3);

   	dWsdX[3][2] = (+ frr[3]/Kc[3]*X[1] + frr[4]/Kc[4]*X[4] - deriv14
	               - frr[14]*X[3]*third14 - 2.0e0*deriv15
	              )*1.0e+6*_calM(3);

   	dWsdX[3][3] = (- frr[3]*X[4] - frr[4]*X[0] - frr[8]*X[6]
	               - frr[9]/Kc[9]*X[5] - deriv14 - frr[14]*X[2]*third14
	               - 2.0e0*deriv15 - 4.0e0*frr[15]*X[3]*third
	              )*1.0e+6*_calM(3);

 	dWsdX[3][4] = (- frr[3]*X[3] + frr[4]/Kc[4]*X[2] + frr[8]/Kc[8]*X[1]
	               + 2.0e0*frr[9]*X[4] - deriv14 + frr[14]/Kc[14]*third14
	               - 2.0e0*deriv15
	              )*1.0e+6*_calM(3);

   	dWsdX[3][5] = (- frr[9]/Kc[9]*X[3] - 5.0e0*deriv14 - 2.0e0*deriv15
	              )*1.0e+6*_calM(3);

  	dWsdX[3][6] = (- frr[8]*X[3] - deriv14 - 2.0e0*deriv15
	             )*1.0e+6*_calM(3);

  	dWsdX[3][7] = (- deriv14 - 2.0e0*deriv15
	              )*1.0e+6*_calM(3);

   	dWsdX[3][8] = (- deriv14 - 2.0e0*deriv15
	              )*1.0e+6*_calM(3);

   	dWsdX[4][0] =(+ 2.0e0*frr[1]*X[1] - frr[2]*X[4] + frr[4]*X[3] - deriv13
	              + deriv14 + 2.0e0*deriv18
	             )*1.0e+6*_calM(4);

   	dWsdX[4][1] = (+ 2.0e0*frr[1]*X[0] + frr[3]/Kc[3]*X[2]
	               + frr[6]/Kc[6]*X[5] - frr[8]/Kc[8]*X[4] - deriv13
	               + deriv14 + 2.0e0*deriv18
	              )*1.0e+6*_calM(4);

  	dWsdX[4][2] = (+ frr[2]/Kc[2]*X[5] + frr[3]/Kc[3]*X[1]
	               - frr[4]/Kc[4]*X[4] + 2.0e0*frr[7]*X[6] - deriv13
	               - frr[13]*X[4]*third13 + deriv14 + frr[14]*X[3]*third14
	               + 2.0e0*deriv18
	              )*1.0e+6*_calM(4);

   	dWsdX[4][3] = (- frr[3]*X[4] + frr[4]*X[0] + frr[8]*X[6]
	               + 2.0e0*frr[9]/Kc[9]*X[5] - deriv13 + deriv14
	               + frr[14]*X[2]*third14 + 2.0e0*deriv18
	              )*1.0e+6*_calM(4);

  	dWsdX[4][4] = (- 4.0e0*frr[1]/Kc[1]*X[4] - frr[2]*X[0] - frr[3]*X[3]
	               - frr[4]/Kc[4]*X[2] - frr[6]*X[6]
	               - 4.0e0*frr[7]/Kc[7]*X[4] - frr[8]/Kc[8]*X[1]
	               - 4.0e0*frr[9]*X[4] - deriv13 - frr[13]*X[2]*third13
	               + deriv14 - frr[14]/Kc[14]*third14 + 2.0e0*deriv18
	               - 4.0e0*frr[18]/Kc[18]*X[4]*third - frr[20]*X[7]
	              )*1.0e+6*_calM(4);

  	dWsdX[4][5] = (+ frr[2]/Kc[2]*X[2] + frr[6]/Kc[6]*X[1]
	               + 2.0e0*frr[9]/Kc[9]*X[3] - 5.0e0*deriv13
	               + frr[13]/Kc[13]*third13 + 5.0e0*deriv14 + 2.0e0*deriv18
	               + frr[20]/Kc[20]*X[6]
	              )*1.0e+6*_calM(4);

  	dWsdX[4][6] = (- frr[6]*X[4] + 2.0e0*frr[7]*X[2] + frr[8]*X[3] - deriv13
	               + deriv14 + 2.0e0*deriv18 + frr[20]/Kc[20]*X[5]
	              )*1.0e+6*_calM(4);

   	dWsdX[4][7] = (- deriv13 + deriv14 + 2.0e0*deriv18 + 2.0e0*frr[18]*third
	               - frr[20]*X[4]
	              )*1.0e+6*_calM(4);

  	dWsdX[4][8] = (- deriv13 + deriv14 + 2.0e0*deriv18
	              )*1.0e+6*_calM(4);

  	dWsdX[5][0] = (+ frr[2]*X[4] + deriv13
	              )*1.0e+6*_calM(5);

   	dWsdX[5][1] = (- frr[6]/Kc[6]*X[5] + deriv13
	              )*1.0e+6*_calM(5);

   	dWsdX[5][2] = (- frr[2]/Kc[2]*X[5] + deriv13 + frr[13]*X[4]*third13
	              )*1.0e+6*_calM(5);

   	dWsdX[5][3] = (- frr[9]/Kc[9]*X[5] + deriv13
	              )*1.0e+6*_calM(5);

  	dWsdX[5][4] = (+ frr[2]*X[0] + frr[6]*X[6] + 2.0e0*frr[9]*X[4] + deriv13
	               + frr[13]*X[2]*third13 + frr[20]*X[7]
	              )*1.0e+6*_calM(5);

  	dWsdX[5][5] = (- frr[2]/Kc[2]*X[2] - frr[6]/Kc[6]*X[1]
	               - frr[9]/Kc[9]*X[3] + deriv13 - frr[13]/Kc[13]*third13
	               - frr[20]/Kc[20]*X[6]
	              )*1.0e+6*_calM(5);

  	dWsdX[5][6] = (+ frr[6]*X[4] + deriv13 - frr[20]/Kc[20]*X[5]
	              )*1.0e+6*_calM(5);

  	dWsdX[5][7] = (+ deriv13 + frr[20]*X[4]
	              )*1.0e+6*_calM(5);

   	dWsdX[5][8] = (+ deriv13
	              )*1.0e+6*_calM(5);

  	dWsdX[6][0] = (+ 2.86e0*deriv5 + frr[16]/Kc[16]*X[1]
	               - frr[19]/Kc[19]*X[6]
	              )*1.0e+6*_calM(6);

   	dWsdX[6][1] = (+ deriv5 + frr[5]*X[2]*third5 + frr[6]/Kc[6]*X[5]
	               + frr[8]/Kc[8]*X[4] + frr[16]/Kc[16]*X[0]
	               + 2.0e0*frr[17]/Kc[17]*X[7]
	              )*1.0e+6*_calM(6);

   	dWsdX[6][2] = (+ deriv5 + frr[5]*X[1]*third5 - frr[7]*X[6]
	               - frr[16]*X[6] + frr[19]*X[7]
	              )*1.0e+6*_calM(6);

   	dWsdX[6][3] = (+ deriv5 - frr[8]*X[6]
	              )*1.0e+6*_calM(6);

  	dWsdX[6][4] = (+ deriv5 - frr[6]*X[6] + 2.0e0*frr[7]/Kc[7]*X[4]
	               + frr[8]/Kc[8]*X[1] + frr[20]*X[7]
	              )*1.0e+6*_calM(6);

  	dWsdX[6][5] = (+ 18.6e0*deriv5 + frr[6]/Kc[6]*X[1] - frr[20]/Kc[20]*X[6]
	              )*1.0e+6*_calM(6);

   	dWsdX[6][6] = (+ deriv5 - frr[5]/Kc[5]*third5 - frr[6]*X[4]
	               - frr[7]*X[2] - frr[8]*X[3] - frr[16]*X[2]
	               - 4.0e0*frr[17]*X[6] - frr[19]/Kc[19]*X[0]
	               - frr[20]/Kc[20]*X[5]
	              )*1.0e+6*_calM(6);

   	dWsdX[6][7] = (+ deriv5 + 2.0e0*frr[17]/Kc[17]*X[1] + frr[19]*X[2]
	               + frr[20]*X[4]
	              )*1.0e+6*_calM(6);

  	dWsdX[6][8] = (+ 1.26e0*deriv5
	              )*1.0e+6*_calM(6);

   	dWsdX[7][0] = (- deriv18 + frr[19]/Kc[19]*X[6]
	              )*1.0e+6*_calM(7);

   	dWsdX[7][1] = (- frr[17]/Kc[17]*X[7] - deriv18
	              )*1.0e+6*_calM(7);

   	dWsdX[7][2] = (- deriv18 - frr[19]*X[7]
	              )*1.0e+6*_calM(7);

   	dWsdX[7][3] = (- deriv18
	              )*1.0e+6*_calM(7);

   	dWsdX[7][4] = (- deriv18 + 2.0e0*frr[18]/Kc[18]*X[4]*third
	               - frr[20]*X[7]
	              )*1.0e+6*_calM(7);

   	dWsdX[7][5] = (- deriv18 + frr[20]/Kc[20]*X[6]
	              )*1.0e+6*_calM(7);

   	dWsdX[7][6] = (+ 2.0e0*frr[17]*X[6] - deriv18 + frr[19]/Kc[19]*X[0]
	               + frr[20]/Kc[20]*X[5]
	              )*1.0e+6*_calM(7);

   	dWsdX[7][7] = (- frr[17]/Kc[17]*X[1] - deriv18 - frr[18]*third
	               - frr[19]*X[2] - frr[20]*X[4]
	              )*1.0e+6*_calM(7);

	dWsdX[7][8] = (- deriv18
	              )*1.0e+6*_calM(7);

/*-----------------------------dXj/drhoj--------------------------------*/

	/* Now calculate dXj/drhoj */

	for ( k = 0 ; k < ns; k++ )
	{
	    dXdrho[k] = (1.0e0 / _calM(k)) * 1.0e-6;
   	}

/*-----------dT/drhoj,dT/drho_v,dT/drho_E,dT/drho_K---------------------*/

	/* Initialize variables */

	for ( k = 0; k < ns; k++ )
	{
	    dTdrhok[k] = 0.0e0;
   	}
	
	for ( k = 0; k < nd; k++ )
	{
	    dTdrhoV[k] = 0.0e0;
	}

   	dTdrhoE = 0.0e0;

	find_dT_dx(np, gl, &dTdrhoE, dTdrhok, dTdrhoV);

	/* Initialize Jacobian matrix C */

	for ( k = 0; k < nf; k++ )
	{
	    for ( s = 0; s < nf; s++ )
	    {
	        C[k][s] = 0.0e0;
	    }
	}
	
	for ( k = 0; k < ns; k++ )
	{
	    for ( s = 0; s < ns; s++ )
	    {
	        C[k][s]=(dWsdX[k][s]*dXdrho[s] + dWsdT[k]*dTdrhok[s]);
	    }

	    /* Calculate the terms involving the various components of energy,
	       their location in the C matrix will depend on the number of
	       dimensions     */

	    /* The first term after the reaction species will always be the
	       x-velocity component */

	    C[k][ns] = (dTdrhoV[0] * dWsdT[k]);

	    /* For 2D and 3D, the second term will be the y-velocity component*/	

	    if ( nd != 1 )
	    {
	        C[k][ns+1] = (dTdrhoV[1] * dWsdT[k]);
	
	        if ( nd == 2 )
	        {
	            /* For 2D, the third term will be the internal energy, and
	               the fourth term will be the turbulent kinetic energy   */
	            C[k][ns+2] = (dTdrhoE * dWsdT[k]);
	            C[k][ns+3] = (-dTdrhoE * dWsdT[k]);
	        }
	        else if ( nd == 3 )
	        {
	            /* For 3D, the third term will be the z-velocity component,
	               and the fourth term will be the internal energy, and
	               the fifth term will be the turbulent kinetic energy */
	            C[k][ns+2] = (dTdrhoV[2] * dWsdT[k]);
	            C[k][ns+3] = (dTdrhoE * dWsdT[k]);
	            C[k][ns+4] = (-dTdrhoE * dWsdT[k]);
	        }
	    }
	    /* For 1D, the second term will be the internal energy, and the
	       third term will be the turbulent kinetic energy */
	    else
	    {
	        C[k][ns+1] = (dTdrhoE * dWsdT[k]);
	        C[k][ns+2] = (-dTdrhoE * dWsdT[k]);
	    }
	}

}
