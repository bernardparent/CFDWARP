#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/_kpsi/model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
#include <src/common.h> /* only needed for display_matrix() function */

#define nr 24 /* Note: ns is specified in chem.hh */
#define Rcal 1.986E0  /*[cal/(gmol K)], this is different than 
	calR = 8.314 J/gmol K in _thermo.h*/

static double Cf[nr]=
  {
    7.5e6 /*4*5.36e5*/, 3.0e10, 1.0e11, 2.0e13, 1.0e12, 1.0e14, 3.0e2, 5.0e13, 6.0e11, 7.0e8, 1.0e13,
    2.0e18, 1.0e17, 7.0e19, 2.24e14, 1.74e13, 8.41e13, 5.75e13, 5.6e11, 2.4e15, 6.0e13,
	     /*1.36e14*/1.4e14, 4.0e13, 6.4e9
  };
static double n_f[nr]=
  {
    1.5e0, 1.0e0, 1.5e0, 0.0e0, 0.0e0, 0.0e0, 2.5e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, -1.0e0,
    -1.0e0, -1.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 1.0e0
  };
static double E[nr]=
  {
    7.9e3, 4.5e3, 1.0e4, 0.0e0, 0.0e0, 4.0e3, 0.0e0, 3.0e3, 5.5e3, 2.325e4, 3.5e3, 0.0e0,
    0.0e0, 0.0e0, 8.4e3, 4.73e3, 1.005e4, 9.0e3, 5.4e2, 2.295e4, 0.0e0, /*3.775e4*/3.79e4, 0.0e0, 3.14e3
  };

/* Species Map for Kerosene 19 species, 24 reaction model 
   NOTE:SPECIES 3 AND 10-19 FROM thermo.hh ARE NOT USED IN THIS 
	COMBUSTION MODEL */

/* 0:"O"     1:"O2"     2:"H2"        3:"N2"    4:"H2O"   5:"H"      6:"OH"
   7:"HO2"   8:"C2H2"   9:"C2H4"     10:"CO"   11:"CO2"  12:"CH3"   13:"CH4"   
  14: "H2CO" 15:"*C6H12O"   16:"*C12H24"  17:"N"    18:"NO"   20:"THIRD BODY"*/

static double stnumR[nr][ns+1]=
  {                 
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0},
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
    {1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
    {0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
    {0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
    {1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0},
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}
  };

static double stnumP[nr][ns+1]=
  {                
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,5,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,0,0,0,0,1,0,2,1,0,1,0,0,0,0,0,0,0},
    {0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
    {0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0},
    {0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
    {0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0},
    {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
    {0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0},
    {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}
  };

void write_chem_template(FILE **controlfile){

}

void read_and_init_chem_actions(char *action, char **argum, SOAP_codex_t *codex){
	
}

double _Xk(long k, np_t np){
  /* species concentration [mol/cm^3] */
  double tmp;
  tmp=_w(np,k)*_rho(np)/_calM(k)*1.0e-6;
  return(tmp);
}

double _Gk(long k, np_t np){
  /* gibbs energy of species [J/mol] */
  double tmp,T;
  T=_T(np);
  tmp=(_hk_from_T(k,T)-T*_sk_from_T(k,T))*_calM(k);
  return(tmp);
}

double _deltaGk(long r, np_t np){
  /* find the difference in gibbs energy for a given reaction */
  long k;
  double tmp;
  tmp=0.0e0;
  for (k=0; k<ns; k++){
    if(stnumP[r][k]!=0 || stnumR[r][k]!=0){
      tmp=tmp+(stnumP[r][k]*_Gk(k,np)-stnumR[r][k]*_Gk(k,np));
    }
  }
  return(tmp);
}

double _Kp(long r, np_t np){
  /* find the equilibrium constant (based on pressures) for each reaction */
  /* UNITS: deltaGk[N m/mol], calR[N m/mol K], T[K] therefore Kp [no units] */
  double tmp,T;
  T=_T(np);
  tmp=exp(-_deltaGk(r,np)/(calR*T));
  return(tmp);
  }

double _Kc(long r, np_t np){
  /* convert Kp into Kc (based on concentrations) for each reaction 
     UNITS:  calR[N m/mol K], T[K], 1.0e6[cm^3/m^3], 101325[N/m^2 atm], 
     therefore Kc multiples of [cm^3/mol] */
  long k;
  double tmp,sum,T;
  sum=0.0e0;
  T=_T(np);
  for (k=0; k<ns; k++){
    if(stnumP[r][k]!=0 || stnumR[r][k]!=0){
      sum=sum+(stnumP[r][k]-stnumR[r][k]);
    }
  }
  tmp=_Kp(r,np)*pow(1/(calR*T*1.0e6/101325),sum);
  return(tmp);
}

double _kf(long r, np_t np){
  /* forward reaction rate based on the modified Arrhenius equation 
     [{(cm^3/mol)^x}/s] where x is a whole number multiple depending on the 
     reaction under consideration */
  double tmp,T;
  T=_T(np);
  tmp=Cf[r]*pow(T,n_f[r])*exp(-E[r]/T);
  return(tmp);
}

double _dGkdT(long k, np_t np){
  /* derivative of gibbs energy of species k wrt temperature [J/(mol K)] */
  double tmp,T;
  T=_T(np);
  tmp=(_cpk_from_T(k,T)-_dTskdT_from_T(k,T))*_calM(k);
  return(tmp);
}

double _ddeltaGkdT(long r, np_t np){
  /* find the difference in the derivatives of gibbs energies wrt 
     temperature for a given reaction */
  long k;
  double tmp;
  tmp=0.0e0;
  for (k=0; k<ns; k++){
    if(stnumP[r][k]!=0 || stnumR[r][k]!=0){
      tmp=tmp+(stnumP[r][k]*_dGkdT(k,np)-stnumR[r][k]*_dGkdT(k,np));
    }
  }
  return(tmp);
}

double _dkfdT(long r, np_t np){
  /* derivative of the forward reaction rate wrt the density of species k */
  double tmp,T;
  T=_T(np);
  tmp=_kf(r,np)*(n_f[r]/T + E[r]/(T*T));
  return(tmp);
}

double _dKcdT(long r, np_t np){
  /* derivative of the equilibrium constant (based on concentrations) wrt 
     the density of species k */
  long k;
  double tmp,sum,T;
  T=_T(np);
  sum=0.0e0;
  for (k=0; k<ns; k++){
    if(stnumP[r][k]!=0 || stnumR[r][k]!=0){
      sum=sum+(stnumP[r][k]-stnumR[r][k]);
    }
  }
  tmp=_Kc(r,np)*(_deltaGk(r,np)/(calR*T*T) -_ddeltaGkdT(r,np)/(calR*T) 
		 - sum/T);
  return(tmp);
}

double _dXdrhok(long k){
  /* derivative of the concentration input in [mol/cm^3] of species k 
     wrt the density [kg/m^3] of species k */
  double tmp;
  tmp=1.0e-6/_calM(k);
  return(tmp);
}

void find_Schem(np_t np, flux_t S){
  long flux;
  long r, k, p;
  double T, sum;
  double third;
  double sum2, sum3;
  double dXdt[ns];

  T=_T(np);
  for (flux=0; flux<nf; flux++){
    S[flux]=0.0e0;
  }

  /* calculate the sum of the concentrations for use in third body 
     calculations */
  third=0.0e0;
  for (k=0; k<ns; k++){
    third=third+_Xk(k,np);
  }

  /* find the time rate of change of concentrations d[X]/dt [mol/cm^3 s] 
     for each species */
  for (k=0; k<ns; k++){
    sum=0.0e0;
    for (r=0; r<nr; r++){
      sum2=1.0e0;
      sum3=1.0e0;
      for (p=0; p<ns; p++){
	sum2=sum2*pow(_Xk(p,np),stnumR[r][p]);
	sum3=sum3*pow(_Xk(p,np),stnumP[r][p]);
      }
      if (stnumR[r][ns]==1){
	sum2=sum2*third;
	sum3=sum3*third;
      }
      /* note that for the first three reactions (0,1,2) Kc is not used as 
	 there is no backwards reaction (Kerosene model) */
      if (r<3){
	sum=sum+(stnumP[r][k]-stnumR[r][k])*(_kf(r,np)*sum2);
      }
      else{
	sum=sum+(stnumP[r][k]-stnumR[r][k])*(_kf(r,np)*sum2 
					     - _kf(r,np)/_Kc(r,np)*sum3);
      }
    }
    dXdt[k]=sum;
  }
  /* convert dXdt [mol_k/cm^3 s] into dwdt [kg_k/m^3 s] and insert into 
     chemical source term */
  for (k=0; k<ns; k++){
    S[k]=(1.0e6)*_calM(k)*dXdt[k]*_Omega(np);
  }
					   
} /* end of find_Schem */

void find_dSchem_dU(np_t np, sqmat_t C){
  
  long r, k, p, q;
  long trip;
  double sumR[nr],sumP[nr];
  double sum2,sum3,sum4[ns];
  double sumA, sumB;
  double stcoeff[nr][ns];
  double third;
  double tbderiv[ns];
  double dwkdT[ns];
  double dwkdX[ns][ns];
  double dTdrhoE;
  spec_t dTdrhok;
  dim_t dTdrhoV;
  sqmat_t C2;

  for (k=0; k<nf; k++){
    for (p=0; p<nf; p++){
      C[k][p]=0.0e0;
    }
  }

  for (r=0; r<nr; r++){
    sumR[r]=1.0e0;
    sumP[r]=1.0e0;
    for (k=0; k<ns; k++){
      sumR[r]=sumR[r]*pow(_Xk(k,np),stnumR[r][k]);
      sumP[r]=sumP[r]*pow(_Xk(k,np),stnumP[r][k]);
    }
  }

  third=0.0e0;
  for (k=0; k<ns; k++){
    third=third+_Xk(k,np);

    for (p=0; p<ns; p++){
      dwkdX[k][p]=0.0e0;
    }

    tbderiv[k]=0.0e0;
    for (r=0; r<nr; r++){
      stcoeff[r][k]=stnumP[r][k]-stnumR[r][k];
      
      if (stcoeff[r][k]!=0 && stnumR[r][ns]!=0){
	tbderiv[k]=tbderiv[k]+stcoeff[r][k]*_kf(r,np)*
	  (sumR[r]-sumP[r]/_Kc(r,np));
      } 
    }
  }

  for (k=0; k<ns; k++){

    trip=0;

    /* find the derivative of dw/dt wrt temperature */
    sumA=0.0e0;
    sumB=0.0e0;
    for (r=0; r<nr; r++){
      sum2=sumR[r];
      sum3=sumP[r];
      if (stcoeff[r][k]!=0){
	if (stnumR[r][ns]!=0){
	  sum2=sumR[r]*third;
	  sum3=sumP[r]*third;
	}  
	if (r<3){
	  sumA=sumA+stcoeff[r][k]*_dkfdT(r,np)*sum2;
	  sumB=0;
	}
	else{
	  sumA=sumA+stcoeff[r][k]*(_dkfdT(r,np)*(sum2 - sum3/_Kc(r,np)));
	  sumB=sumB+stcoeff[r][k]*(_kf(r,np)/_Kc(r,np))*
	    sum3*_dKcdT(r,np)/_Kc(r,np);
	}
      }
    }
    dwkdT[k]=(sumB+sumA)*1.0e6*_calM(k);

    /* find the derivative of dw/dt wrt concentration */
    for (r=0; r<nr; r++){ 
      if (stcoeff[r][k]!=0){
	for (p=0; p<ns; p++){
	  sum4[p]=0.0e0;
	  if (stnumR[r][p]!=0){
	    sum2=1.0e0;
	    for (q=0; q<ns; q++){
	      if (q==p){
		sum2=sum2*stnumR[r][q]*pow(_Xk(q,np),stnumR[r][q]-1);
	      }
	      else{
		sum2=sum2*pow(_Xk(q,np),stnumR[r][q]);
	      }
	    }
	    if (stnumR[r][ns]==1){
	      sum2=sum2*third;
	    }	
	    sum4[p]=sum4[p]+stcoeff[r][k]*(_kf(r,np)*sum2);
	  }
	}
	for (p=0; p<ns; p++){
	  dwkdX[k][p]=dwkdX[k][p]+1.0e6*_calM(k)*sum4[p];
	}
	
	for (p=0; p<ns; p++){
	  sum4[p]=0.0e0;
	  if (stnumP[r][p]!=0 && r>2){
	    sum3=1.0e0;
	    for (q=0; q<ns; q++){
	      if (q==p){
		sum3=sum3*stnumP[r][q]*pow(_Xk(q,np),stnumP[r][q]-1);
	      }
	      else{
		sum3=sum3*pow(_Xk(q,np),stnumP[r][q]);
	      }
	    }
	    if (stnumP[r][ns]==1){
	      sum3=sum3*third;
	    }
	    sum4[p]=sum4[p]+stcoeff[r][k]*(-_kf(r,np)*sum3/_Kc(r,np));
	  }
	}
	for (p=0; p<ns; p++){
	  dwkdX[k][p]=dwkdX[k][p]+1.0e6*_calM(k)*sum4[p];
	}
	
	if (stnumR[r][ns]==1 && trip==0){
	  trip=1;
	  for (p=0; p<ns; p++){
	    dwkdX[k][p]=dwkdX[k][p]+1.0e6*_calM(k)*tbderiv[k];
	  }
	}
      }
    }
  } /* end of loop though the species k=0->(ns-1) */

  /* assemble the terms for the chemical Jacobian */
  find_dT_dx(np, &dTdrhoE, dTdrhok, dTdrhoV);
  for (k=0; k<ns; k++){
    /* derivative of S[k] wrt to rhok[ns] */
    for (p=0; p<ns; p++){
      C[k][p]=(dwkdX[k][p]*_dXdrhok(p) + dwkdT[k]*dTdrhok[p]);
    }
    /* derivative of S[k] wrt to rhov[nd] */
    for (p=ns; p<ns+nd; p++){
      C[k][p]=(dwkdT[k]*dTdrhoV[p-ns]);
    } 
    /* derivative of S[k] wrt to rhoE */
    C[k][ns+nd]=(dwkdT[k]*dTdrhoE);
    /* derivative of S[k] wrt to rhok */
    C[k][ns+nd+1]=(-dTdrhoE*dwkdT[k]);
  }

  /* This section checks the analytical with the numerical Jacobian */
  /*
    find_numerical_jacobian_2(&(find_Schem), np, C2);
    for (k=0; k<ns+nd+1; k++){
    for (p=0; p<ns+nd+1; p++){
    printf("\nAna C[%ld][%ld] = %.7e\tNum C2[%ld][%ld] = %.7e\tDiff = %.8e",k,p,C[k][p],k,p,C2[k][p],
    (C[k][p]-C2[k][p])/C2[k][p]);
    }
    }
  */

} /* end of find_dSchem_dU */
