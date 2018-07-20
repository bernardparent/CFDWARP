#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/.active/model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>

#define nr 19


const static long mR[nr][ns]=
 {/*e-  O2  N2  O   N   O2+ N2+ O+  N+  O2-   */
   {1,  0,  0,  0,  0,  1,  0,  0,  0,  0 },  /* 1a */
   {1,  0,  0,  0,  0,  0,  1,  0,  0,  0 },  /* 1b */
   {0,  0,  0,  0,  0,  0,  1,  0,  0,  1 },  /* 2a */
   {0,  0,  0,  0,  0,  1,  0,  0,  0,  1 },  /* 2b */
   {0,  0,  0,  0,  0,  0,  0,  0,  1,  1 },  /* 2c */
   {0,  0,  0,  0,  0,  0,  0,  1,  0,  1 },  /* 2d */
   {0,  0,  1,  0,  0,  0,  1,  0,  0,  1 },  /* 3a */
   {0,  0,  1,  0,  0,  1,  0,  0,  0,  1 },  /* 3b */
   {0,  0,  1,  0,  0,  0,  0,  0,  1,  1 },  /* 3c */
   {0,  0,  1,  0,  0,  0,  0,  1,  0,  1 },  /* 3d */
   {0,  1,  0,  0,  0,  0,  1,  0,  0,  1 },  /* 3e */
   {0,  1,  0,  0,  0,  1,  0,  0,  0,  1 },  /* 3f */
   {0,  1,  0,  0,  0,  0,  0,  0,  1,  1 },  /* 3g */
   {0,  1,  0,  0,  0,  0,  0,  1,  0,  1 },  /* 3h */
   {1,  2,  0,  0,  0,  0,  0,  0,  0,  0 },  /* 4a */
   {1,  1,  1,  0,  0,  0,  0,  0,  0,  0 },  /* 4b */
   {0,  1,  0,  0,  0,  0,  0,  0,  0,  1 },  /* 5  */
   {1,  0,  1,  0,  0,  0,  0,  0,  0,  0 },  /* 6a  */
   {1,  1,  0,  0,  0,  0,  0,  0,  0,  0 },  /* 6b  */
 };

 
const static long mP[nr][ns]=
 {/*e-  O2  N2  O   N   O2+ N2+ O+  N+  O2-   */
   {0,  0,  0,  2,  0,  0,  0,  0,  0,  0 },  /* 1a */
   {0,  0,  0,  0,  2,  0,  0,  0,  0,  0 },  /* 1b */
   {0,  1,  1,  0,  0,  0,  0,  0,  0,  0 },  /* 2a */
   {0,  2,  0,  0,  0,  0,  0,  0,  0,  0 },  /* 2b */
   {0,  1,  0,  0,  1,  0,  0,  0,  0,  0 },  /* 2c */
   {0,  1,  0,  1,  0,  0,  0,  0,  0,  0 },  /* 2d */
   {0,  1,  2,  0,  0,  0,  0,  0,  0,  0 },  /* 3a */
   {0,  2,  1,  0,  0,  0,  0,  0,  0,  0 },  /* 3b */
   {0,  1,  1,  0,  1,  0,  0,  0,  0,  0 },  /* 3c */
   {0,  1,  1,  1,  0,  0,  0,  0,  0,  0 },  /* 3d */
   {0,  2,  1,  0,  0,  0,  0,  0,  0,  0 },  /* 3e */
   {0,  3,  0,  0,  0,  0,  0,  0,  0,  0 },  /* 3f */
   {0,  2,  0,  0,  1,  0,  0,  0,  0,  0 },  /* 3g */
   {0,  2,  0,  1,  0,  0,  0,  0,  0,  0 },  /* 3h */
   {0,  1,  0,  0,  0,  0,  0,  0,  0,  1 },  /* 4a */
   {0,  0,  1,  0,  0,  0,  0,  0,  0,  1 },  /* 4b */
   {1,  2,  0,  0,  0,  0,  0,  0,  0,  0 },  /* 5  */
   {2,  0,  0,  0,  0,  0,  1,  0,  0,  0 },  /* 6a  */
   {2,  0,  0,  0,  0,  1,  0,  0,  0,  0 },  /* 6b  */
 };
  


void write_chem_template(FILE **controlfile){
}


void read_and_init_chem_actions(char *action, char **argum, SOAP_codex_t *codex){

}


void find_Schem(np_t np, gl_t *gl, flux_t S){
  long flux,spec,reaction,spec2;
  double kf[nr];
  spec_t nk,dnkdt;
  double Te,T,Tv,kv;
  double vartheta;
  
  for (flux=0; flux<nf; flux++) S[flux]=0.0;
  
  /* find the gas temperature and the electron temperature */
  T=_T(np);
  Te=_Te(np,gl);
  Tv=_Tv(np);
  vartheta=_Efieldstar(np,gl)/1.0e-20;
  kv=pow(10.0,-28.3*exp(-3353.0/Tv)/sqr(vartheta));
      
  /* set the reaction rates for each reaction */
  kf[0]=2.0e-7*pow(300.0/Te,0.7);  /* 1a */
  kf[1]=2.8e-7*pow(300.0/Te,0.5);  /* 1b */
  kf[2]=2.0e-7*pow(300.0/T,0.5);   /* 2a */
  kf[3]=kf[2];                     /* 2b */
  kf[4]=kf[2];                     /* 2c */
  kf[5]=kf[2];                     /* 2d */
  kf[6]=2.0e-25*pow(300.0/T,2.5);  /* 3a */
  kf[7]=kf[6];                     /* 3b */
  kf[8]=kf[6];                     /* 3c */
  kf[9]=kf[6];                     /* 3d */
  kf[10]=kf[6];                    /* 3e */
  kf[11]=kf[6];                    /* 3f */
  kf[12]=kf[6];                    /* 3g */
  kf[13]=kf[6];                    /* 3h */
  kf[14]=1.4e-29*300.0/Te*exp(-600.0/T)*exp(700.0*(Te-T)/Te/T); /* 4a */
  kf[15]=1.07e-31*sqr(300.0/Te)*exp(-70.0/T)*exp(1500.0*(Te-T)/Te/T); /* 4b */
  kf[16]=8.6e-10*exp(-6030.0/T)*(1.0-exp(-1570.0/T)); /* 5 */
  kf[17]=kv*pow(10.0,-8.3-36.5/vartheta);  /* 6a */
  kf[18]=kv*pow(10.0,-8.8-28.1/vartheta);  /* 6b */

  /* find the mole fraction Xk in moles/cm^3 */
  for (spec=0; spec<ns; spec++){
    nk[spec]=_rhok(np,spec)/_calM(spec)*1e-6*calA;
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
    /* add dXkdt to S */ 
    for (spec=0; spec<ns; spec++){
      S[spec]+=_Omega(np,gl)*1.0e6*_calM(spec)*dnkdt[spec]/calA;
    }
  }
}


void find_numerical_jacobian_chem(void(*FUNCT)(np_t, gl_t *, flux_t), np_t np, gl_t *gl, sqmat_t Ak){
 long row,col,flux;
 flux_t F,Fnew,Ustar,dUstar,dUstar2;

 for (flux=0; flux<nf; flux++){
   Ustar[flux]=np.bs->U[flux]*_Omega(np,gl);
   dUstar[flux]=Ustar[flux]/1000000.0e0+1.0e-15*_Omega(np,gl);
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
   (*FUNCT)(np, gl, Fnew);
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



void find_dSchem_dU(np_t np, gl_t *gl, sqmat_t C){
  long row,col;
  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      C[row][col]=0.0;
    }
  }	
  find_numerical_jacobian_chem(&find_Schem, np,gl, C);
}
