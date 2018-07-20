#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
#include <src/common.h> /* only needed for display_matrix() function */

#define nr 24 /* Note: ns is specified in chem.hh */
#define Rcal 1.986E0  /*[cal/(gmol K)], this is different than calR =
			8.314 J/gmol K in _thermo.h*/

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

void write_chem_template(FILE **controlfile){

}

void read_and_init_chem_actions(char *action, char **argum, SOAP_codex_t *codex){

}

double _Xk(long k, double rhok){
  /* species concentration [mol/cm^3] */
  double tmp;
  tmp=rhok/_calM(k)*1.0e-6;
  return(tmp);
}

double _Gk(long k, double T){
  /* gibbs energy of species [J/mol] */
  double tmp;
  tmp=(_hk_from_T(k,T)-T*_sk_from_T(k,T))*_calM(k);
  return(tmp);
}

double _kf(long r, double T){
  /* forward reaction rate based on the modified Arrhenius equation 
     [{(cm^3/mol)^x}/s] where x is a whole number multiple depending on the 
     reaction under consideration */
  double tmp;
  tmp=Cf[r]*pow(T,n_f[r])*exp(-E[r]/T);
  return(tmp);
}

double _dGkdT(long k, double T){
  /* derivative of gibbs energy of species k wrt temperature [J/(mol K)] */
  double tmp;
  tmp=(_cpk_from_T(k,T)-_dTskdT_from_T(k,T))*_calM(k);
  return(tmp);
}

double _dkfdT(long r, double T){
  /* derivative of the forward reaction rate wrt the density of species k */
  double tmp;
  tmp=_kf(r,T)*(n_f[r]/T + E[r]/(T*T));
  return(tmp);
}


void find_W(spec_t rhok, double T, double Te, double Tv, double Efieldstar, double Qbeam, spec_t W){
  long r, k;
  double X[ns], Gs[ns];
  double dGs[nr], Kc[nr], kf[nr], react[nr];
  double third;

  for (r=0; r<ns; r++){
    W[r]=0.0e0;
  }

  third=0.0e0;
  for(k=0; k<ns; k++){
    X[k]=_Xk(k,rhok[k]);
    Gs[k]=_Gk(k,T);
    third=third + X[k];
  }


  /* Calculate the difference between the gibbs free energy between the 
     products and reactants of each reaction */

  dGs[0]= 2.0 * Gs[15] - (Gs[16] + Gs[1]) ;
  dGs[1]= Gs[14] + Gs[12] + 5.0 * Gs[9] - (Gs[16] + Gs[6]);
  dGs[2]= Gs[7] + Gs[10] + Gs[12] + 2.0 * Gs[9] - (Gs[15] + Gs[1]);
  dGs[3]= Gs[14] + Gs[5] - (Gs[12] + Gs[0]);
  dGs[4]= Gs[14] + Gs[6] - (Gs[12] + Gs[1]);
  dGs[5]= Gs[4] + Gs[10] + Gs[5] - (Gs[14] + Gs[6]);
  dGs[6]= 2.0 * Gs[14] - (Gs[9] + Gs[1]);
  dGs[7]= Gs[12] + Gs[14] - (Gs[9] + Gs[6]);
  dGs[8]= Gs[13] + Gs[5] - (Gs[12] + Gs[2]);
  dGs[9]= Gs[8] + Gs[2] - Gs[9];
  dGs[10]= Gs[12] + Gs[10] - (Gs[8] + Gs[6]);
  dGs[11]= Gs[2] - 2.0 * Gs[5];
  dGs[12]= Gs[1] - 2.0 * Gs[0];
  dGs[13]= Gs[4] - (Gs[6] + Gs[5]);
  dGs[14]= Gs[6] + Gs[0] - (Gs[5] + Gs[1]);
  dGs[15]= Gs[6] + Gs[5] - (Gs[0] + Gs[2]);
  dGs[16]= Gs[2] + Gs[6] - (Gs[5] + Gs[4]);
  dGs[17]= 2.0 * Gs[6] - (Gs[0] + Gs[4]);
  dGs[18]= Gs[11] + Gs[5] - (Gs[10] + Gs[6]);
  dGs[19]= Gs[5] + Gs[1] - Gs[7];
  dGs[20]= 2.0 * Gs[6] - (Gs[7] + Gs[5]);
  dGs[21]= Gs[17] + Gs[18] - (Gs[3] + Gs[0]);
  dGs[22]= Gs[5] + Gs[18] - (Gs[17] + Gs[6]);
  dGs[23]= Gs[18] + Gs[0] - (Gs[17] + Gs[1]);

  /* determine Kp (storing results in Kc since this is ultimately what is 
     sought here */

  for(r=0; r<nr; r++){
    Kc[r]=exp(-dGs[r]/(calR*T));
    kf[r]=_kf(r,T);
  }

  /* determine Kc by multiplying Kp by (RT/101325*100**3) ** (prod-react) 
     [Kc] = (cm**3/mol) ** (prod-react) */

  Kc[5] = Kc[5] / (82.056872e0 * T);
  Kc[9] = Kc[9] / (82.056872e0 * T);
  Kc[11] = Kc[11] * (82.056872e0 * T);
  Kc[12] = Kc[12] * (82.056872e0 * T);
  Kc[13] = Kc[13] * (82.056872e0 * T);
  Kc[19] = Kc[19] / (82.056872e0 * T);


  /* calculate the individual reactions */

  react[0] = kf[0] * X[16] * X[1];
  react[1] = kf[1] * X[16] * X[6];
  react[2] = kf[2] * X[15] * X[1];
  react[3] = kf[3] * (X[12] * X[0] - X[14] * X[5] / Kc[3]);
  react[4] = kf[4] * (X[12] * X[1] - X[14] * X[6] / Kc[4]);
  react[5] = kf[5] * (X[14] * X[6] - X[4] * X[10] * X[5] / Kc[5]);
  react[6] = kf[6] * (X[9] * X[1] - X[14] * X[14] / Kc[6]);
  react[7] = kf[7] * (X[9] * X[6] - X[12] * X[14] / Kc[7]);
  react[8] = kf[8] * (X[12] * X[2] - X[13] * X[5] / Kc[8]);
  react[9] = kf[9] * (X[9] - X[8] * X[2] / Kc[9]);
  react[10] = kf[10] * (X[8] * X[6] - X[12] * X[10] / Kc[10]);
  react[11] = kf[11] * (X[5] * X[5] * third - X[2] * third / Kc[11]);
  react[12] = kf[12] * (X[0] * X[0] * third - X[1] * third / Kc[12]);
  react[13] = kf[13] * (X[6] * X[5] * third - X[14] * third / Kc[13]);
  react[14] = kf[14] * (X[5] * X[1] - X[6] * X[0] / Kc[14]);
  react[15] = kf[15] * (X[0] * X[2] - X[6] * X[5] / Kc[15]);
  react[16] = kf[16] * (X[5] * X[4] - X[2] * X[6] / Kc[16]);
  react[17] = kf[17] * (X[0] * X[4] - X[6] * X[6] / Kc[17]);
  react[18] = kf[18] * (X[10] * X[6] - X[11] * X[5] / Kc[18]);
  react[19] = kf[19] * (X[7] * third - X[5] * X[1] * third / Kc[19]);
  react[20] = kf[20] * (X[7] * X[5] - X[6] * X[6] / Kc[20]);
  react[21] = kf[21] * (X[0] * X[3] - X[17] * X[18] / Kc[21]);
  react[22] = kf[22] * (X[17] * X[6] - X[5] * X[18]/ Kc[22]);
  react[23] = kf[23] * (X[17] * X[1] - X[18] * X[0] / Kc[23]);
  
  /* calculate the chemical source term [kg_k/(m^3 s)] */

  W[0] = (-1.0 * react[3] - 2.0 * react[12] + react[14] - react[15] - react[17] - react[21] + react[23]) 
	    * 1.0e6 * _calM(0);
  W[1] = (-1.0 * react[0] - react[2] - react[4] - react[6] +
	    react[12] - react[14] + react[19] - react[23]) 
    * 1.0e6 * _calM(1);
  W[2] = (-1.0 * react[8] + react[9] + react[11] - react[15] + react[16])
    * 1.0e6 * _calM(2);
  W[3] = (-1.0 * react[21]) * 1.0e6 * _calM(3);
  W[4] = (react[13] - react[16] - react[17] + react[5]) * 1.0e6 * _calM(4);
  W[5] = (react[3] + react[5] + react[8] - 2.0 * react[11] -
	   react[13] - react[14] + react[15] - react[16] +
	   react[18] + react[19] - react[20] + react[22])
    * 1.0e6 * _calM(5);
  W[6] = (-1.0 * react[1] + react[4] - react[5] - react[7] -
	    react[10] - react[13] + react[14] + react[15] +
	    react[16] + 2.0 * react[17] - react[18] - react[22] + 2.0 * react[20])
    * 1.0e6 * _calM(6);
  W[7] = (react[2] - react[19] - react[20]) * 1.0e6 * _calM(7);
  W[8] = (react[9] - react[10]) * 1.0e6 * _calM(8);
  W[9] = (5.0 * react[1] + 2.0 * react[2] - react[6] - react[7] - react[9])
    * 1.0e6 * _calM(9);
  W[10] = (react[2] + react[5] + react[10] - react[18]) * 1.0e6 * _calM(10);
  W[11] = react[18] * 1.0e6 * _calM(11);
  W[12] = (react[1] + react[2] - react[3] - react[4] + react[7] - react[8] + react[10]) 
    * 1.0e6 * _calM(12);
  W[13] = react[8] * 1.0e6 * _calM(13);
  W[14] = (react[1] + react[3] + react[4] - react[5] + 2.0 * react[6] + react[7])
    * 1.0e6 * _calM(14);
  W[15] = (2.0 * react[0] - react[2]) * 1.0e6 * _calM(15);
  W[16] = -1.0 * (react[0] + react[1]) * 1.0e6 * _calM(16);
  W[17] = (react[21] - react[22] - react[23]) * 1.0e6 * _calM(17);
  W[18] = (react[21] + react[22] + react[23]) * 1.0e6 * _calM(18);


} 

  
void find_dW_dx(spec_t rhok, double T, double Te, double Tv, double Efieldstar, double qbeam,
              spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdEfieldstar, 
	      spec_t dWdQbeam){
  long r, k, p;
  double third;
  double Kc[nr], kf[nr], kb[nr], dGs[nr], react[nr], dKcdT[nr];
  double X[ns], Gs[ns], dGsdT[ns];
  double dWsdX[ns][ns];


  /* initialize derivatives to zero */
  for (r=0; r<ns; r++){
    dWdT[r]=0.0;
    dWdTe[r]=0.0;
    dWdTv[r]=0.0;
    dWdEfieldstar[r]=0.0;
    dWdQbeam[r]=0.0;
    for (k=0; k<ns; k++){
      dWdrhok[r][k]=0.0;
    }
  }
  
  third=0.0;
  for(k=0; k<ns; k++){
    X[k]=_Xk(k,rhok[k]);
    Gs[k]=_Gk(k,T);
    dGsdT[k]=_dGkdT(k,T);
    third=third + X[k];
  }

  /* Calculate the difference between the gibbs free energy between the products 
   and reactants of each reaction */

  dGs[0]= 2.0 * Gs[15] - (Gs[16] + Gs[1]) ;
  dGs[1]= Gs[14] + Gs[12] + 5.0 * Gs[9] - (Gs[16] + Gs[6]);
  dGs[2]= Gs[7] + Gs[10] + Gs[12] + 2.0 * Gs[9] - (Gs[15] + Gs[1]);
  dGs[3]= Gs[14] + Gs[5] - (Gs[12] + Gs[0]);
  dGs[4]= Gs[14] + Gs[6] - (Gs[12] + Gs[1]);
  dGs[5]= Gs[4] + Gs[10] + Gs[5] - (Gs[14] + Gs[6]);
  dGs[6]= 2.0 * Gs[14] - (Gs[9] + Gs[1]);
  dGs[7]= Gs[12] + Gs[14] - (Gs[9] + Gs[6]);
  dGs[8]= Gs[13] + Gs[5] - (Gs[12] + Gs[2]);
  dGs[9]= Gs[8] + Gs[2] - Gs[9];
  dGs[10]= Gs[12] + Gs[10] - (Gs[8] + Gs[6]);
  dGs[11]= Gs[2] - 2.0 * Gs[5];
  dGs[12]= Gs[1] - 2.0 * Gs[0];
  dGs[13]= Gs[4] - (Gs[6] + Gs[5]);
  dGs[14]= Gs[6] + Gs[0] - (Gs[5] + Gs[1]);
  dGs[15]= Gs[6] + Gs[5] - (Gs[0] + Gs[2]);
  dGs[16]= Gs[2] + Gs[6] - (Gs[5] + Gs[4]);
  dGs[17]= 2.0 * Gs[6] - (Gs[0] + Gs[4]);
  dGs[18]= Gs[11] + Gs[5] - (Gs[10] + Gs[6]);
  dGs[19]= Gs[5] + Gs[1] - Gs[7];
  dGs[20]= 2.0 * Gs[6] - (Gs[7] + Gs[5]);
  dGs[21]= Gs[17] + Gs[18] - (Gs[3] + Gs[0]);
  dGs[22]= Gs[5] + Gs[18] - (Gs[17] + Gs[6]);
  dGs[23]= Gs[0] + Gs[18] - (Gs[17] + Gs[1]);


  /* determine Kc by multiplying Kp by (RT/101325*100**3) ** (prod-react) 
     [Kc] = (cm**3/mol) ** (prod-react) */

  for(r=0; r<nr; r++){
    Kc[r]=exp(-dGs[r]/(calR*T)); 
  }

  Kc[5] = Kc[5] / (82.056872e0 * T);
  Kc[9] = Kc[9] / (82.056872e0 * T);
  Kc[11] = Kc[11] * (82.056872e0 * T);
  Kc[12] = Kc[12] * (82.056872e0 * T);
  Kc[13] = Kc[13] * (82.056872e0 * T);
  Kc[19] = Kc[19] / (82.056872e0 * T);
   
  for(r=0; r<nr; r++){
    kf[r]=_kf(r,T);
    kb[r]=kf[r]/Kc[r];
  }

  /* find the derivative of each source term wrt temperature */
 
  react[0] = _dkfdT(0,T) * X[16] * X[1];
  react[1] = _dkfdT(1,T) * X[16] * X[6];
  react[2] = _dkfdT(2,T) * X[15] * X[1];
  react[3] = _dkfdT(3,T) * (X[12] * X[0] - X[14] * X[5] / Kc[3]);
  react[4] = _dkfdT(4,T) * (X[12] * X[1] - X[14] * X[6] / Kc[4]);
  react[5] = _dkfdT(5,T) * (X[14] *X[6] - X[4] * X[10] * X[5] / Kc[5]);
  react[6] = _dkfdT(6,T) * (X[9] * X[1] - X[14] * X[14] / Kc[6]);
  react[7] = _dkfdT(7,T) * (X[9] * X[6] - X[12] * X[14] / Kc[7]);
  react[8] = _dkfdT(8,T) * (X[12] * X[2] - X[13] * X[5] / Kc[8]);
  react[9] = _dkfdT(9,T) * (X[9] - X[8] * X[2] / Kc[9]);
  react[10] = _dkfdT(10,T) * (X[8] * X[6] - X[12] * X[10] / Kc[10]);
  react[11] = _dkfdT(11,T) * (X[5] * X[5] * third - X[2] * third / Kc[11]);
  react[12] = _dkfdT(12,T) * (X[0] * X[0] * third - X[1] * third / Kc[12]);
  react[13] = _dkfdT(13,T) * (X[6] * X[5] * third - X[4] * third / Kc[13]);
  react[14] = _dkfdT(14,T) * (X[5] * X[1] - X[6] * X[0] / Kc[14]);
  react[15] = _dkfdT(15,T) * (X[0] * X[2] - X[6] * X[5] / Kc[15]);
  react[16] = _dkfdT(16,T) * (X[5] * X[4] - X[2] * X[6] / Kc[16]);
  react[17] = _dkfdT(17,T) * (X[0] * X[4] - X[6] * X[6] / Kc[17]);
  react[18] = _dkfdT(18,T) * (X[10] * X[6] - X[11] * X[5] / Kc[18]);
  react[19] = _dkfdT(19,T) * (X[7] * third - X[5] * X[1] * third / Kc[19]);
  react[20] = _dkfdT(20,T) * (X[7] * X[5] - X[6] * X[6] / Kc[20]);
  react[21] = _dkfdT(21,T) * (X[0] * X[3] - X[17] * X[18] / Kc[21]);
  react[22] = _dkfdT(22,T) * (X[17] * X[6] - X[5] * X[18] / Kc[22]);
  react[23] = _dkfdT(23,T) * (X[17] * X[1] - X[18] * X[0] / Kc[23]);

  /**************/
  dWdT[0] = (-1.0 * react[3] - 2.0 * react[12] + react[14] - react[15] - react[17] 
	      - react[21] + react[23]) 
    * 1.0e6 * _calM(0);
  dWdT[1] = (-1.0 * react[0] - react[2] - react[4] - react[6] +
	      react[12] - react[14] + react[19] - react[23]) 
    * 1.0e6 * _calM(1);
  dWdT[2] = (-1.0 * react[8] + react[9] + react[11] - react[15] + react[16])
    * 1.0e6 * _calM(2);
  dWdT[3] = (-1.0 * react[21]) * 1.0e6 * _calM(3);
  dWdT[4] = (react[13] - react[16] - react[17] + react[5]) * 1.0e6 * _calM(4);
  dWdT[5] = (react[3] + react[5] + react[8] - 2.0 * react[11] -
	      react[13] - react[14] + react[15] - react[16] +
	      react[18] + react[19] - react[20] + react[22])
    * 1.0e6 * _calM(5);
  dWdT[6] = (-1.0 * react[1] + react[4] - react[5] - react[7] -
	      react[10] - react[13] + react[14] + react[15] +
	      react[16] + 2.0 * react[17] - react[18] - react[22] + 2.0 * react[20])
    * 1.0e6 * _calM(6);
  dWdT[7] = (react[2] - react[19] - react[20]) * 1.0e6 * _calM(7);
  dWdT[8] = (react[9] - react[10]) * 1.0e6 * _calM(8);
  dWdT[9] = (5.0 * react[1] + 2.0 * react[2] - react[6] - react[7] - react[9])
    * 1.0e6 * _calM(9);
  dWdT[10] = (react[2] + react[5] + react[10] - react[18]) * 1.0e6 * _calM(10);
  dWdT[11] = react[18] * 1.0e6 * _calM(11);
  dWdT[12] = (react[1] + react[2] - react[3] - react[4] + react[7] - react[8] + react[10]) 
    * 1.0e6 * _calM(12);
  dWdT[13] = react[8] * 1.0e6 * _calM(13);
  dWdT[14] = (react[1] + react[3] + react[4] - react[5] + 2.0 * react[6] + react[7])
    * 1.0e6 * _calM(14);
  dWdT[15] = (2.0 * react[0] - react[2]) * 1.0e6 * _calM(15);
  dWdT[16] = -1.0 * (react[0] + react[1]) * 1.0e6 * _calM(16);
  dWdT[17] = (react[21] - react[22] - react[23]) * 1.0e6 * _calM(17);
  dWdT[18] = (react[21] + react[22] + react[23]) * 1.0e6 * _calM(18);

  /*****************/
  dKcdT[0]= 0.0e0;
  dKcdT[1]= 0.0e0;
  dKcdT[2]= 0.0e0;
  dKcdT[3]= dGs[3]/(calR*T*T) - (dGsdT[14] + dGsdT[5] - (dGsdT[12] + dGsdT[0])) / calR / T;
  dKcdT[4]= dGs[4]/(calR*T*T) - (dGsdT[14] + dGsdT[6] - (dGsdT[12] + dGsdT[1])) / calR / T;
  dKcdT[5]= -1.0e0 / T + dGs[5]/(calR*T*T) - (dGsdT[4] + dGsdT[10] + dGsdT[5] - 
					      (dGsdT[14] + dGsdT[6])) / calR / T;
  dKcdT[6]= dGs[6]/(calR*T*T) - (2.0e0 * dGsdT[14] - (dGsdT[9] + dGsdT[1])) / calR / T;
  dKcdT[7]= dGs[7]/(calR*T*T) - (dGsdT[12] + dGsdT[14] - (dGsdT[9] + dGsdT[6])) / calR / T;
  dKcdT[8]= dGs[8]/(calR*T*T) - (dGsdT[13] + dGsdT[5] - (dGsdT[12] + dGsdT[2])) / calR / T;
  dKcdT[9]= -1.0e0 / T + dGs[9]/(calR*T*T) - (dGsdT[8] + dGsdT[2] - dGsdT[9]) / calR / T;
  dKcdT[10]= dGs[10]/(calR*T*T) - (dGsdT[12] + dGsdT[10] - (dGsdT[8] + dGsdT[6])) / calR / T;
  dKcdT[11]= 1.0e0 / T + dGs[11]/(calR*T*T) - (dGsdT[2] - 2.0e0 * dGsdT[5]) / calR / T;
  dKcdT[12]= 1.0e0 / T + dGs[12]/(calR*T*T) - (dGsdT[1] - 2.0e0 * dGsdT[0]) / calR /T;
  dKcdT[13]= 1.0e0 / T + dGs[13]/(calR*T*T) - (dGsdT[4] - (dGsdT[6] + dGsdT[5])) / calR / T;
  dKcdT[14]= dGs[14]/(calR*T*T) - (dGsdT[6] + dGsdT[0] - (dGsdT[5] + dGsdT[1])) / calR / T;
  dKcdT[15]= dGs[15]/(calR*T*T) - (dGsdT[6] + dGsdT[5] - (dGsdT[0] + dGsdT[2])) / calR / T;
  dKcdT[16]= dGs[16]/(calR*T*T) - (dGsdT[2] + dGsdT[6] - (dGsdT[5] + dGsdT[4])) / calR / T;
  dKcdT[17]= dGs[17]/(calR*T*T) - (2.0e0 * dGsdT[6] - (dGsdT[0] + dGsdT[4])) / calR / T;
  dKcdT[18]= dGs[18]/(calR*T*T) - (dGsdT[11] + dGsdT[5] - (dGsdT[10] + dGsdT[6])) / calR / T;
  dKcdT[19]= -1.0e0 / T + dGs[19]/(calR*T*T) - (dGsdT[5] + dGsdT[1] - dGsdT[7]) / calR / T;
  dKcdT[20]= dGs[20]/(calR*T*T) - (2.0e0 * dGsdT[6] - (dGsdT[7] + dGsdT[5])) / calR / T;
  dKcdT[21]= dGs[21]/(calR*T*T) - (dGsdT[17] + dGsdT[18] - (dGsdT[0] + dGsdT[3])) / calR / T;
  dKcdT[22]= dGs[22]/(calR*T*T) - (dGsdT[5] + dGsdT[18] - (dGsdT[17] + dGsdT[6])) / calR / T;
  dKcdT[23]= dGs[23]/(calR*T*T) - (dGsdT[18] + dGsdT[0] - (dGsdT[17] + dGsdT[1])) / calR / T;

  /******************/
  react[0] = 0.0e0;
  react[1] = 0.0e0;
  react[2] = 0.0e0;
  react[3] = kf[3] / Kc[3] * X[14] * X[5] * dKcdT[3];
  react[4] = kf[4] / Kc[4] * X[14] * X[6] * dKcdT[4];
  react[5] = kf[5] / Kc[5] * X[4] * X[10] * X[5] * dKcdT[5];
  react[6] = kf[6] / Kc[6] * X[14] * X[14] * dKcdT[6];
  react[7] = kf[7] / Kc[7] * X[12] * X[14] * dKcdT[7];
  react[8] = kf[8] / Kc[8] * X[13] * X[5] * dKcdT[8];
  react[9] = kf[9] / Kc[9] * X[8] * X[2] * dKcdT[9];
  react[10] = kf[10] / Kc[10] * X[12] * X[10] * dKcdT[10];
  react[11] = kf[11] / Kc[11] * X[2] * third * dKcdT[11];
  react[12] = kf[12] / Kc[12] * X[1] * third * dKcdT[12];
  react[13] = kf[13] / Kc[13] * X[4] * third * dKcdT[13];
  react[14] = kf[14] / Kc[14] * X[6] * X[0] * dKcdT[14];
  react[15] = kf[15] / Kc[15] * X[6] * X[5] * dKcdT[15];
  react[16] = kf[16] / Kc[16] * X[2] * X[6] * dKcdT[16];
  react[17] = kf[17] / Kc[17] * X[6] * X[6] * dKcdT[17];
  react[18] = kf[18] / Kc[18] * X[8] * X[5] * dKcdT[18];
  react[19] = kf[19] / Kc[19] * X[5] * X[1] * third * dKcdT[19];
  react[20] = kf[20] / Kc[20] * X[6] * X[6] * dKcdT[20];
  react[21] = kf[21] * X[17] * X[18] / Kc[21] * dKcdT[21];
  react[22] = kf[22] * X[5] * X[18] / Kc[22] * dKcdT[22];
  react[23] = kf[23] * X[18] * X[0] / Kc[23] * dKcdT[23];

  /* get the final answer */
  dWdT[0] = dWdT[0] + (-1.0 * react[3] - 2.0 * react[12] + react[14] - react[15] - react[17] 
			 - react[21] + react[23]) * 1.0e6 * _calM(0);
  dWdT[1] = dWdT[1] + (-1.0 * react[0] - react[2] - react[4] - react[6] +
			 react[12] - react[14] + react[19] - react[23]) 
                         * 1.0e6 * _calM(1);
  dWdT[2] = dWdT[2] + (-1.0 * react[8] + react[9] + react[11] - react[15] + react[16])
                         * 1.0e6 * _calM(2);
  dWdT[3] = dWdT[3] + (-1.0 * react[21]) * 1.0e6 * _calM(3);
  dWdT[4] = dWdT[4] + (react[13] - react[16] - react[17] + react[5]) * 1.0e6 * _calM(4);
  dWdT[5] = dWdT[5] + (react[3] + react[5] + react[8] - 2.0 * react[11] -
			 react[13] - react[14] + react[15] - react[16] +
			 react[18] + react[19] - react[20] + react[22])
                         * 1.0e6 * _calM(5);
  dWdT[6] = dWdT[6] + (-1.0 * react[1] + react[4] - react[5] - react[7] -
			 react[10] - react[13] + react[14] + react[15] +
			 react[16] + 2.0 * react[17] - react[18] - react[22] + 2.0 * react[20])
                         * 1.0e6 * _calM(6);
  dWdT[7] = dWdT[7] + (react[2] - react[19] - react[20]) * 1.0e6 * _calM(7);
  dWdT[8] = dWdT[8] + (react[9] - react[10]) * 1.0e6 * _calM(8);
  dWdT[9] = dWdT[9] + (5.0 * react[1] + 2.0 * react[2] - react[6] - react[7] - react[9])
                         * 1.0e6 * _calM(9);
  dWdT[10] = dWdT[10] + (react[2] + react[5] + react[10] - react[18]) 
                           * 1.0e6 * _calM(10);
  dWdT[11] = dWdT[11] + react[18] * 1.0e6 * _calM(11);
  dWdT[12] = dWdT[12] + (react[1] + react[2] - react[3] - react[4] + react[7] 
			   - react[8] + react[10]) * 1.0e6 * _calM(12);
  dWdT[13] = dWdT[13] + react[8] * 1.0e6 * _calM(13);
  dWdT[14] = dWdT[14] + (react[1] + react[3] + react[4] - react[5] + 2.0 * react[6] + react[7])
                           * 1.0e6 * _calM(14);
  dWdT[15] = dWdT[15] + (2.0 * react[0] - react[2]) * 1.0e6 * _calM(15);
  dWdT[16] = dWdT[16] + -1.0 * (react[0] + react[1]) * 1.0e6 * _calM(16);
  dWdT[17] = dWdT[17] + (react[21] - react[22] - react[23]) * 1.0e6 * _calM(17);
  dWdT[18] = dWdT[18] + (react[21] + react[22] + react[23]) * 1.0e6 * _calM(18);

  /* find the derivative of each source term wrt concentration */

  for (k=0; k<ns; k++){
    for (p=0; p<ns; p++){
      dWsdX[k][p]=0.0e0;
    }
  }

    dWsdX[0][0] = 1.0e6*_calM(0) * (-1.0e0*kf[3]*X[12] 
					      - 4.0e0*kf[12]*X[0]*third - kb[14]*X[6]
					      - kf[15]*X[2] - kf[17]*X[4] - kf[21]*X[3]
					      - kb[23]*X[18]);
    dWsdX[0][1] = 1.0e6*_calM(0) * (2.0e0*kb[12]*third + kf[14]*X[5] + kf[23]*X[17]);
    dWsdX[0][2] = 1.0e6*_calM(0) * (-1.0e0*kf[15]*X[0]);
    dWsdX[0][3] = 1.0e6*_calM(0) * (-1.0e0*kf[21]*X[0]);
    dWsdX[0][4] = 1.0e6*_calM(0) * (-1.0e0*kf[17]*X[0]);
    dWsdX[0][5] = 1.0e6*_calM(0) * (kb[3]*X[14] + kf[14]*X[1] + kb[15]*X[6]);
    dWsdX[0][6] = 1.0e6*_calM(0) * (-1.0e0*kb[14]*X[0] + kb[15]*X[5] + kb[17]*2.0e0*X[6]);
    dWsdX[0][12] = 1.0e6*_calM(0) * (-1.0e0*kf[3]*X[0]);
    dWsdX[0][14] = 1.0e6*_calM(0) * (kb[3]*X[5]);
    dWsdX[0][17] = 1.0e6*_calM(0) * (kb[21]*X[18] + kf[23]*X[1]);
    dWsdX[0][18] = 1.0e6*_calM(0) * (kb[21]*X[17] - kb[23]*X[0]);
    for(k=0; k<ns; k++){
      dWsdX[0][k] = dWsdX[0][k] + 1.0e6*_calM(0)*(-2.0e0*kf[12]*X[0]*X[0] 
							     + 2.0e0*kb[12]*X[1]);
    }

    dWsdX[1][0] = 1.0e6*_calM(1) * (2.0e0*kf[12]*X[0]*third + kb[14]*X[6] + kb[23]*X[18]);
    dWsdX[1][1] = 1.0e6*_calM(1) * (-1.0e0*kf[0]*X[16] - kf[2]*X[15]
					      - kf[4]*X[12] - kf[6]*X[9] - kb[12]*third
					      - kf[14]*X[5] - kb[19]*X[5]*third - kf[23]*X[17]);
    dWsdX[1][5] = 1.0e6*_calM(1) * (-1.0e0*kf[14]*X[1] - kb[19]*X[1]*third);
    dWsdX[1][6] = 1.0e6*_calM(1) * (kb[4]*X[14] + kb[14]*X[0]);
    dWsdX[1][7] = 1.0e6*_calM(1) * (kf[19]*third);
    dWsdX[1][9] = 1.0e6*_calM(1) * (-1.0e0*kf[6]*X[1]);
    dWsdX[1][12] = 1.0e6*_calM(1) * (-1.0e0*kf[4]*X[1]);
    dWsdX[1][14] = 1.0e6*_calM(1) * (kb[4]*X[6] + 2.0e0*kb[6]*X[14]);
    dWsdX[1][15] = 1.0e6*_calM(1) * (-1.0e0*kf[2]*X[1]);
    dWsdX[1][16] = 1.0e6*_calM(1) * (-1.0e0*kf[0]*X[1]);
    dWsdX[1][17] = 1.0e6*_calM(1) * (-kf[23]*X[1]);
    dWsdX[1][18] = 1.0e6*_calM(1) * (kb[23]*X[0]);
    for(k=0; k<ns; k++){
      dWsdX[1][k] = dWsdX[1][k] + 1.0e6*_calM(1) * (kf[12]*X[0]*X[0]
							     - kb[12]*X[1] + kf[19]*X[7]
							     - kb[19]*X[5]*X[1]);
    }

    dWsdX[2][0] = 1.0e6*_calM(2) * (-1.0e0*kf[15]*X[2]);
    dWsdX[2][2] = 1.0e6*_calM(2) * (-1.0e0*kf[8]*X[12] - kb[9]*X[8]
					      -kb[11]*third - kf[15]*X[0] - kb[16]*X[6]);
    dWsdX[2][4] = 1.0e6*_calM(2) * (kf[16]*X[5]);
    dWsdX[2][5] = 1.0e6*_calM(2) * (kb[8]*X[13] + 2.0e0*kf[11]*X[5]*third
					      + kb[15]*X[6] + kf[15]*X[4]);
    dWsdX[2][6] = 1.0e6*_calM(2) * (kb[15]*X[5] - kb[16]*X[2]);
    dWsdX[2][8] = 1.0e6*_calM(2) * (-1.0e0*kb[9]*X[2]);
    dWsdX[2][9] = 1.0e6*_calM(2) * (kf[9]);
    dWsdX[2][12] = 1.0e6*_calM(2) * (-1.0e0*kf[8]*X[2]);
    dWsdX[2][13] = 1.0e6*_calM(2) * (kb[8]*X[5]);
    for(k=0; k<ns; k++){
      dWsdX[2][k] = dWsdX[2][k] + 1.0e6*_calM(2) * (kf[11]*X[5]*X[5] - kb[11]*X[2]);
    }
  
    dWsdX[3][0] = 1.0e6*_calM(3) * (-1.0e0*kf[21]*X[3]);
    dWsdX[3][3] = 1.0e6*_calM(3) * (-1.0e0*kf[21]*X[0]);
    dWsdX[3][17] = 1.0e6*_calM(3) * (kb[21]*X[18]);
    dWsdX[3][18] = 1.0e6*_calM(3) * (kb[21]*X[17]);

    dWsdX[4][0] = 1.0e6*_calM(4) * (-1.0e0*kf[17]*X[4]);
    dWsdX[4][2] = 1.0e6*_calM(4) * (kb[16]*X[6]);
    dWsdX[4][4] = 1.0e6*_calM(4) * (-1.0e0*kb[13]*third - kf[16]*X[5]
					      - kf[17]*X[0] - kb[5]*X[10]*X[5]);
    dWsdX[4][5] = 1.0e6*_calM(4) * (-1.0e0*kf[15]*X[4]+kf[13]*X[6]*third- kb[5]*X[4]*X[10]);
    dWsdX[4][6] = 1.0e6*_calM(4) * (kf[13]*X[5]*third + kb[16]*X[2]
					      + 2.0e0*kb[17]*X[6] + kf[5]*X[14]);
    dWsdX[4][10] = 1.0e6*_calM(4) * (-kb[5]*X[4]*X[5]);
    dWsdX[4][14] = 1.0e6*_calM(4) * (kf[5]*X[6]);
    for(k=0; k<ns; k++){
      dWsdX[4][k] = dWsdX[4][k] + 1.0e6*_calM(4) * (kf[13]*X[6]*X[5] - kb[13]*X[4]);
    }
      
    dWsdX[5][0] = 1.0e6*_calM(5) * (kf[3]*X[12] + kb[14]*X[6] + kf[15]*X[2]);
    dWsdX[5][1] = 1.0e6*_calM(5) * (-1.0e0*kf[14]*X[5]-kb[19]*X[5]*third);
    dWsdX[5][2] = 1.0e6*_calM(5) * (kf[8]*X[12] + 2.0e0*kb[11]*third
					      + kf[15]*X[0] + kb[16]*X[6]);
    dWsdX[5][4] = 1.0e6*_calM(5) * (-1.0e0*kb[5]*X[10]*X[5] + kb[13]*third - kf[16]*X[5]);
    dWsdX[5][5] = 1.0e6*_calM(5) * (-1.0e0*kb[3]*X[14] - kb[5]*X[4]*X[10]
					      - kb[8]*X[13] - 4.0e0*kf[11]*X[5]*third
					      - kf[13]*X[6]*third - kf[14]*X[1]
					      - kb[15]*X[6] - kf[15]*X[4] - kb[18]*X[11]
					      - kb[19]*X[1]*third - kf[20]*X[7]
					      - kb[22]*X[18]);
    dWsdX[5][6] = 1.0e6*_calM(5) * (kf[5]*X[14] - kf[13]*X[5]*third
					      + kb[14]*X[0] - kb[15]*X[5] + kb[16]*X[2]
					      + kf[18]*X[10] + 2.0e0*kb[20]*X[6]
					      + kf[22]*X[17]);
    dWsdX[5][7] = 1.0e6*_calM(5) * (kf[19]*third - kf[20]*X[5]);
    dWsdX[5][10] = 1.0e6*_calM(5) * (-1.0e0*kb[5]*X[4]*X[5] + kf[18]*X[6]);
    dWsdX[5][11] = 1.0e6*_calM(5) * (-1.0e0*kb[18]*X[5]);
    dWsdX[5][12] = 1.0e6*_calM(5) * (kf[3]*X[0] + kf[8]*X[2]);
    dWsdX[5][13] = 1.0e6*_calM(5) * (-1.0e0*kb[8]*X[5]);
    dWsdX[5][14] = 1.0e6*_calM(5) * (-1.0e0*kb[3]*X[5] + kf[5]*X[6]);
    dWsdX[5][17] = 1.0e6*_calM(5) * (kf[22]*X[6]);
    dWsdX[5][18] = 1.0e6*_calM(5) * (-kb[22]*X[5]);
    for(k=0; k<ns; k++){
      dWsdX[5][k] = dWsdX[5][k] + 1.0e6*_calM(5) * 
	(-2.0e0*kf[11]*X[5]*X[5] + 2.0e0*kb[11]*X[2] - kf[13]*X[6]*X[5]
	 + kb[13]*X[4] + kf[19]*X[7]- kb[19]*X[5]*X[1]);
    }
 
    dWsdX[6][0] = 1.0e6*_calM(6) * (-1.0e0*kb[14]*X[6] + kf[15]*X[2] + 2.0e0*kf[17]*X[4]);
    dWsdX[6][1] = 1.0e6*_calM(6) * (kf[4]*X[12] + kf[14]*X[5]);
    dWsdX[6][2] = 1.0e6*_calM(6) * (kf[15]*X[0] - kb[16]*X[6]);
    dWsdX[6][4] = 1.0e6*_calM(6) * (kb[5]*X[10]*X[5]+kb[13]*third+
					      kf[16]*X[5]+2.0e0*kf[17]*X[0]);
    dWsdX[6][5] = 1.0e6*_calM(6) * (kb[5]*X[4]*X[10] 
					      - kf[13]*X[6]*third + kf[14]*X[1] - kb[15]*X[6]
					      + kf[15]*X[4] + kb[18]*X[11] + 2.0e0*kf[20]*X[7]
					      + kb[22]*X[18]);
    dWsdX[6][6] = 1.0e6*_calM(6) * (-1.0e0*kf[1]*X[16] - kb[4]*X[14]
					      - kf[5]*X[14] - kf[7]*X[9] - kf[10]*X[8]
					      - kf[13]*X[5]*third - kb[14]*X[0]
					      - kb[15]*X[5] - kb[16]*X[2] - 4.0e0*kb[17]*X[6]
					      - kf[18]*X[10] - 4.0e0*kb[20]*X[6] - kf[22]*X[17]);
    dWsdX[6][7] = 1.0e6*_calM(6) * (2.0e0*kf[20]*X[5]);
    dWsdX[6][8] = 1.0e6*_calM(6) * (-1.0e0*kf[10]*X[6]);
    dWsdX[6][9] = 1.0e6*_calM(6) * (-1.0e0*kf[7]*X[6]);
    dWsdX[6][10] = 1.0e6*_calM(6) * (kb[5]*X[4]*X[5] + kb[10]*X[12] - kf[18]*X[6]);
    dWsdX[6][11] = 1.0e6*_calM(6) * (kb[18]*X[5]);
    dWsdX[6][12] = 1.0e6*_calM(6) * (kf[4]*X[1] + kb[7]*X[14] + kb[10]*X[10]);
    dWsdX[6][14] = 1.0e6*_calM(6) * (-1.0e0*kb[4]*X[6] - kf[5]*X[6] + kb[7]*X[12]);
    dWsdX[6][16] = 1.0e6*_calM(6) * (-1.0e0*kf[1]*X[6]);
    dWsdX[6][17] = 1.0e6*_calM(6) * (-1.0e0*kf[22]*X[6]);
    dWsdX[6][18] = 1.0e6*_calM(6) * (kb[22]*X[5]);
    for(k=0; k<ns; k++){
      dWsdX[6][k] = dWsdX[6][k] + 1.0e6*_calM(6) *
	(-1.0e0*kf[13]*X[6]*X[5] + kb[13]*X[4]);
    }

    dWsdX[7][1] = 1.0e6*_calM(7) * (kb[19]*X[5]*third + kf[2]*X[15]);
    dWsdX[7][5] = 1.0e6*_calM(7) * (kb[19]*X[1]*third - kf[20]*X[7]);
    dWsdX[7][6] = 1.0e6*_calM(7) * (2.0e0*kb[20]*X[6]);
    dWsdX[7][7] = 1.0e6*_calM(7) * (-1.0e0*kf[19]*third - kf[20]*X[5]);
    dWsdX[7][15] = 1.0e6*_calM(7) * (kf[2]*X[1]);
    for(k=0; k<ns; k++){
      dWsdX[7][k] = dWsdX[7][k] + 1.0e6*_calM(7) * (-1.0e0*kf[19]*X[7] + kb[19]*X[5]*X[1]);
    }
  
    dWsdX[8][2] = 1.0e6*_calM(8) * (-1.0e0*kb[9]*X[8]);
    dWsdX[8][6] = 1.0e6*_calM(8) * (-1.0e0*kf[10]*X[8]);
    dWsdX[8][8] = 1.0e6*_calM(8) * (-1.0e0*kb[9]*X[2] - kf[10]*X[6]);
    dWsdX[8][9] = 1.0e6*_calM(8) * (kf[9]);
    dWsdX[8][10] = 1.0e6*_calM(8) * (kb[10]*X[12]);
    dWsdX[8][12] = 1.0e6*_calM(8) * (kb[10]*X[10]);

    dWsdX[9][1] = 1.0e6*_calM(9) * (-1.0e0*kf[6]*X[9] + 2.0e0*kf[2]*X[15]);
    dWsdX[9][2] = 1.0e6*_calM(9) * (kb[9]*X[8]);
    dWsdX[9][6] = 1.0e6*_calM(9) * (-1.0e0*kf[7]*X[9] + 5.0e0*kf[1]*X[16]);
    dWsdX[9][8] = 1.0e6*_calM(9) * (kb[9]*X[2]);
    dWsdX[9][9] = 1.0e6*_calM(9) * (-1.0e0*kf[6]*X[1] - kf[7]*X[6] - kf[9]);
    dWsdX[9][12] = 1.0e6*_calM(9) * (kb[7]*X[14]);
    dWsdX[9][14] = 1.0e6*_calM(9) * (2.0e0*kb[6]*X[14] + kb[7]*X[12]);
    dWsdX[9][15] = 1.0e6*_calM(9) * (2.0e0*kf[2]*X[1]);
    dWsdX[9][16] = 1.0e6*_calM(9) * (5.0e0*kf[1]*X[6]);

    dWsdX[10][1] = 1.0e6*_calM(10) * (kf[2]*X[15]);
    dWsdX[10][4] = 1.0e6*_calM(10) * (-1.0e0*kb[5]*X[10]*X[5]);
    dWsdX[10][5] = 1.0e6*_calM(10) * (-1.0e0*kb[5]*X[4]*X[10] + kb[18]*X[11]);
    dWsdX[10][6] = 1.0e6*_calM(10) * (kf[5]*X[14] + kf[10]*X[8] - kf[18]*X[10]);
    dWsdX[10][8] = 1.0e6*_calM(10) * (kf[10]*X[6]);
    dWsdX[10][10] = 1.0e6*_calM(10) * (-1.0e0*kb[5]*X[4]*X[5] - kb[10]*X[12] - kf[18]*X[6]);
    dWsdX[10][11] = 1.0e6*_calM(10) * (kb[18]*X[5]);
    dWsdX[10][12] = 1.0e6*_calM(10) * (-1.0e0*kb[10]*X[10]);
    dWsdX[10][14] = 1.0e6*_calM(10) * (kf[5]*X[6]);
    dWsdX[10][15] = 1.0e6*_calM(10) * (kf[2]*X[1]);

    dWsdX[11][5] = 1.0e6*_calM(11) * (-1.0e0*kb[18]*X[11]);
    dWsdX[11][6] = 1.0e6*_calM(11) * (kf[18]*X[10]);
    dWsdX[11][10] = 1.0e6*_calM(11) * (kf[18]*X[6]);
    dWsdX[11][11] = 1.0e6*_calM(11) * (-1.0e0*kb[18]*X[5]);

    dWsdX[12][0] = 1.0e6*_calM(12) * (-1.0e0*kf[3]*X[12]);
    dWsdX[12][1] = 1.0e6*_calM(12) * (kf[2]*X[15]-kf[4]*X[12]);
    dWsdX[12][2] = 1.0e6*_calM(12) * (-1.0e0*kf[8]*X[12]);
    dWsdX[12][5] = 1.0e6*_calM(12) * (kb[3]*X[14]+kb[8]*X[13]);
    dWsdX[12][6] = 1.0e6*_calM(12) * (kf[1]*X[16] + kb[4]*X[14] + 
						kf[7]*X[9]+kf[10]*X[8]);
    dWsdX[12][8] = 1.0e6*_calM(12) * (kf[10]*X[6]);
    dWsdX[12][9] = 1.0e6*_calM(12) * (kf[7]*X[6]);
    dWsdX[12][10] = 1.0e6*_calM(12) * (-1.0e0*kb[10]*X[12]);
    dWsdX[12][12] = 1.0e6*_calM(12) * (-1.0e0*kf[3]*X[0] - kf[4]*X[1] -
						 kb[7]*X[14] - kf[8]*X[2] - kb[10]*X[10]);
    dWsdX[12][13] = 1.0e6*_calM(12) * (kb[8]*X[5]);
    dWsdX[12][14] = 1.0e6*_calM(12) * (kb[3]*X[5] + kb[4]*X[6] -
						 kb[7]*X[12]);
    dWsdX[12][15] = 1.0e6*_calM(12) * (kf[2]*X[1]);
    dWsdX[12][16] = 1.0e6*_calM(12) * (kf[1]*X[6]);
  
    dWsdX[13][2] = 1.0e6*_calM(13) * (kf[8]*X[12]);
    dWsdX[13][5] = 1.0e6*_calM(13) * (-1.0e0*kb[8]*X[13]);
    dWsdX[13][12] = 1.0e6*_calM(13) * (kf[8]*X[2]);	
    dWsdX[13][13] = 1.0e6*_calM(13) * (-1.0e0*kb[8]*X[5]);
 
    dWsdX[14][0] = 1.0e6*_calM(14) * (kf[3]*X[12]);
    dWsdX[14][1] = 1.0e6*_calM(14) * (kf[4]*X[12] + 2.0e0*kf[6]*X[9]);
    dWsdX[14][4] = 1.0e6*_calM(14) * (kb[5]*X[10]*X[5]);
    dWsdX[14][5] = 1.0e6*_calM(14) * (-1.0e0*kb[3]*X[14]+kb[5]*X[4]*X[10]);
    dWsdX[14][6] = 1.0e6*_calM(14) * (-1.0e0*kb[4]*X[14] - kf[5]*X[14]
						+ kf[7]*X[9] + kf[1]*X[16]);
    dWsdX[14][9] = 1.0e6*_calM(14) * (2.0e0*kf[6]*X[1] + kf[7]*X[6]);
    dWsdX[14][10] = 1.0e6*_calM(14) * (kb[5]*X[4]*X[5]);
    dWsdX[14][12] = 1.0e6*_calM(14) * (kf[3]*X[0] + kf[4]*X[1] - kb[7]*X[14]);
    dWsdX[14][14] = 1.0e6*_calM(14) * (-1.0e0*kb[3]*X[5] - kb[4]*X[6]
						 - kf[5]*X[6] - 4.0e0*kb[6]*X[14] - kb[7]*X[12]);
    dWsdX[14][16] = 1.0e6*_calM(14) * (kf[1]*X[6]);

    dWsdX[15][1] = 1.0e6*_calM(15) * (2.0e0*kf[0]*X[16] - kf[2]*X[15]);
    dWsdX[15][15] = 1.0e6*_calM(15) * (-1.0e0*kf[2]*X[1]);
    dWsdX[15][16] = 1.0e6*_calM(15) * (2.0e0*kf[0]*X[1]);

    dWsdX[16][1] = 1.0e6*_calM(16) * (-1.0e0*kf[0]*X[16]);
    dWsdX[16][6] = 1.0e6*_calM(16) * (-1.0e0*kf[1]*X[16]);
    dWsdX[16][16] = 1.0e6*_calM(16) * (-1.0e0*kf[0]*X[1] - kf[1]*X[6]);

    dWsdX[17][0] = 1.0e6*_calM(17) * (kf[21]*X[3] + kb[23]*X[18]);
    dWsdX[17][1] = 1.0e6*_calM(17) * (-kf[23]*X[17]);
    dWsdX[17][3] = 1.0e6*_calM(17) * (kf[21]*X[0]);
    dWsdX[17][5] = 1.0e6*_calM(17) * (kb[22]*X[18]);
    dWsdX[17][6] = 1.0e6*_calM(17) * (-1.0e0*kf[22]*X[17]);
    dWsdX[17][17] = 1.0e6*_calM(17) * (-1.0e0*kb[21]*X[18] - kf[23]*X[1] - kf[22]*X[6]);
    dWsdX[17][18] = 1.0e6*_calM(17) * (-1.0e0*kb[21]*X[17] + kb[23]*X[0] + kb[22]*X[5]);

    dWsdX[18][0] = 1.0e6*_calM(18) * (kf[21]*X[3] - kb[23]*X[18]);
    dWsdX[18][1] = 1.0e6*_calM(18) * (kf[23]*X[17]);
    dWsdX[18][3] = 1.0e6*_calM(18) * (kf[21]*X[0]);
    dWsdX[18][5] = 1.0e6*_calM(18) * (-1.0e0*kb[22]*X[18]);
    dWsdX[18][6] = 1.0e6*_calM(18) * (kf[22]*X[17]);
    dWsdX[18][17] = 1.0e6*_calM(18) * (-1.0e0*kb[21]*X[18] + kf[23]*X[1] + kf[22]*X[6]);
    dWsdX[18][18] = 1.0e6*_calM(18) * (-1.0e0*kb[21]*X[17] - kb[23]*X[0] - kb[22]*X[5]);

  for (k=0; k<ns; k++){
    for (p=0; p<ns; p++){
      dWdrhok[k][p]=dWsdX[k][p]*1.0e-6/_calM(p);
    }
  }
    
} 
