

void FindWOriginal(spec_t rhok, double T, double Te, double Tv, double Efieldstar, double Qbeam, spec_t W){
  long spec,reaction,spec2;
  double kf[nr];
  spec_t nk,dnkdt;
  double vartheta,qi,N;
  
  for (spec=0; spec<ns; spec++) W[spec]=0.0;
  /* find the gas temperature and the electron temperature */
  vartheta=Efieldstar/1.0e-20;
    
  /* set the reaction rates for each reaction */
  kf[0]=pow(10.0,-8.3-36.5/vartheta);  /* 1a */
  kf[1]=pow(10.0,-8.8-28.1/vartheta);  /* 1n */
  kf[2]=2.0e-7*pow(300.0/Te,0.7);  /* 2a */
  kf[3]=2.8e-7*pow(300.0/Te,0.5);  /* 2b */
  kf[4]=2.0e-7*pow(300.0/T,0.5);   /* 3a */
  kf[5]=kf[4];                     /* 3b */
  kf[6]=2.0e-25*pow(300.0/T,2.5);  /* 4a */
  kf[7]=kf[6];                     /* 4b */
  kf[8]=kf[6];                     /* 4c */
  kf[9]=kf[6];                     /* 4d */
  kf[10]=1.4e-29*300.0/Te*exp(-600.0/T)*exp(700.0*(Te-T)/Te/T); /* 5a */
  kf[11]=1.07e-31*sqr(300.0/Te)*exp(-70.0/T)*exp(1500.0*(Te-T)/Te/T); /* 5b */
  kf[12]=8.6e-10*exp(-6030.0/T)*(1.0-exp(-1570.0/T)); /* 6 */
  kf[13]=0.0;
  kf[14]=0.0;
  /* find the mole fraction Xk in moles/cm^3 */
  for (spec=0; spec<ns; spec++){
    nk[spec]=rhok[spec]/_calM(spec)*1e-6*calA;
  } 
  
  for (reaction=0; reaction<nr; reaction++){
    /* find dXkdt for reaction n */
    for (spec=0; spec<ns; spec++){
      dnkdt[spec]=((double)(mP[reaction][spec])-(double)(mR[reaction][spec]))*kf[reaction];
      for (spec2=0; spec2<ns; spec2++) {
        if (mR[reaction][spec2]==1) {
          dnkdt[spec]=dnkdt[spec]*nk[spec2];
        }
        if (mR[reaction][spec2]==2) {
          dnkdt[spec]=dnkdt[spec]*nk[spec2]*nk[spec2];
        }
        if (mR[reaction][spec2]==3) {
          dnkdt[spec]=dnkdt[spec]*nk[spec2]*nk[spec2]*nk[spec2];
        }
      }	
    }
    /* add dXkdt to W*/ 
    for (spec=0; spec<ns; spec++){
      W[spec]+=1.0e6*_calM(spec)*dnkdt[spec]/calA;
    }
  }

  /* ionization rates due to the electron beam */
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=nk[spec];
  qi=Qbeam/34.0/1.6022e-19;
  W[speceminus]+=(nk[specO2]/N+nk[specN2]/N)*qi*_calM(speceminus)/calA; /* eminus */
  W[specO2plus]+=nk[specO2]/N*qi*_calM(specO2plus)/calA;       /* O2plus */
  W[specN2plus]+=nk[specN2]/N*qi*_calM(specN2plus)/calA; /* N2plus */
  W[specO2]-=nk[specO2]/N*qi*_calM(specO2)/calA;
  W[specN2]-=nk[specN2]/N*qi*_calM(specN2)/calA;
}
