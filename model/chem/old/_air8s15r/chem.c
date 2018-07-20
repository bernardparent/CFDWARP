#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/metrics/_metrics.h>
//#include <model/.active/model_eef.h>

#define nr 15

#define ATTACHMENT TRUE
#define TOWNSEND TRUE

const static long mR[nr][ns]=
 {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
   {1,  0,  1,  0,  0,  0,  0,  0 },  /* 1a */
   {1,  1,  0,  0,  0,  0,  0,  0 },  /* 1b */
   {1,  0,  0,  0,  0,  1,  0,  0 },  /* 2a */
   {1,  0,  0,  0,  0,  0,  1,  0 },  /* 2b */
   {0,  0,  0,  0,  0,  0,  1,  1 },  /* 3a */
   {0,  0,  0,  0,  0,  1,  0,  1 },  /* 3b */
   {0,  0,  1,  0,  0,  0,  1,  1 },  /* 4a */
   {0,  0,  1,  0,  0,  1,  0,  1 },  /* 4b */
   {0,  1,  0,  0,  0,  0,  1,  1 },  /* 4c */
   {0,  1,  0,  0,  0,  1,  0,  1 },  /* 4d */
   {1,  2,  0,  0,  0,  0,  0,  0 },  /* 5a */
   {1,  1,  1,  0,  0,  0,  0,  0 },  /* 5b */
   {0,  1,  0,  0,  0,  0,  0,  1 },  /* 6  */
   {0,  1,  0,  0,  0,  0,  0,  0 },  /* 7a */
   {0,  0,  1,  0,  0,  0,  0,  0 },  /* 7b */
 };

 
const static long mP[nr][ns]=
 {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
   {2,  0,  0,  0,  0,  0,  1,  0 },  /* 1a */
   {2,  0,  0,  0,  0,  1,  0,  0 },  /* 1b */
   {0,  0,  0,  2,  0,  0,  0,  0 },  /* 2a */
   {0,  0,  0,  0,  2,  0,  0,  0 },  /* 2b */
   {0,  1,  1,  0,  0,  0,  0,  0 },  /* 3a */
   {0,  2,  0,  0,  0,  0,  0,  0 },  /* 3b */
   {0,  1,  2,  0,  0,  0,  0,  0 },  /* 4a */
   {0,  2,  1,  0,  0,  0,  0,  0 },  /* 4b */
   {0,  2,  1,  0,  0,  0,  0,  0 },  /* 4c */
   {0,  3,  0,  0,  0,  0,  0,  0 },  /* 4d */
   {0,  1,  0,  0,  0,  0,  0,  1 },  /* 5a */
   {0,  0,  1,  0,  0,  0,  0,  1 },  /* 5b */
   {1,  2,  0,  0,  0,  0,  0,  0 },  /* 6  */
   {1,  0,  0,  0,  0,  1,  0,  0 },  /* 7a */
   {1,  0,  0,  0,  0,  0,  1,  0 },  /* 7b */
 };
  




  /* Remember: 
     by definition
       W[0]=qi*nk[1]/N  
     where 
     then
     Wp[13]=kf[13]*nk[1]=qi*nk[1]/N
     or 
     kf[13]=qi/N
     
  */

double _k1a(double EoverN){
  double kdata[62],EoverNdata[62];
  double k,fact;
  long cnt;
  
  EoverN=max(1.01E-20,EoverN); 
  EoverN=min(3.3e-18,EoverN); 
  k=0.0; /* to avoid compiler warning */

  EoverNdata[0]=1.000000E-20;  kdata[0]=1.584893E-45; 
  EoverNdata[1]=1.100000E-20;  kdata[1]=3.297477E-42; 
  EoverNdata[2]=1.210000E-20;  kdata[2]=3.425396E-39; 
  EoverNdata[3]=1.331000E-20;  kdata[3]=1.892386E-36; 
  EoverNdata[4]=1.464100E-20;  kdata[4]=5.888557E-34; 
  EoverNdata[5]=1.610510E-20;  kdata[5]=1.087356E-31; 
  EoverNdata[6]=1.771561E-20;  kdata[6]=1.249400E-29; 
  EoverNdata[7]=1.948717E-20;  kdata[7]=9.326715E-28; 
  EoverNdata[8]=2.143589E-20;  kdata[8]=4.704144E-26; 
  EoverNdata[9]=2.357948E-20;  kdata[9]=1.661257E-24; 
  EoverNdata[10]=2.593742E-20;  kdata[10]=4.242970E-23; 
  EoverNdata[11]=2.853117E-20;  kdata[11]=8.071842E-22; 
  EoverNdata[12]=3.138428E-20;  kdata[12]=1.174830E-20; 
  EoverNdata[13]=3.452271E-20;  kdata[13]=1.340449E-19; 
  EoverNdata[14]=3.797498E-20;  kdata[14]=1.225771E-18; 
  EoverNdata[15]=4.177248E-20;  kdata[15]=9.166218E-18; 
  EoverNdata[16]=4.594973E-20;  kdata[16]=5.708691E-17; 
  EoverNdata[17]=5.054470E-20;  kdata[17]=3.010715E-16; 
  EoverNdata[18]=5.559917E-20;  kdata[18]=1.365068E-15; 
  EoverNdata[19]=6.115909E-20;  kdata[19]=5.394589E-15; 
  EoverNdata[20]=6.727500E-20;  kdata[20]=1.881515E-14; 
  EoverNdata[21]=7.400250E-20;  kdata[21]=5.857795E-14; 
  EoverNdata[22]=8.140275E-20;  kdata[22]=1.644834E-13; 
  EoverNdata[23]=8.954302E-20;  kdata[23]=4.204818E-13; 
  EoverNdata[24]=9.849733E-20;  kdata[24]=9.869966E-13; 
  EoverNdata[25]=1.083471E-19;  kdata[25]=2.143858E-12; 
  EoverNdata[26]=1.191818E-19;  kdata[26]=4.339613E-12; 
  EoverNdata[27]=1.310999E-19;  kdata[27]=8.238811E-12; 
  EoverNdata[28]=1.442099E-19;  kdata[28]=1.475597E-11; 
  EoverNdata[29]=1.586309E-19;  kdata[29]=2.506465E-11; 
  EoverNdata[30]=1.744940E-19;  kdata[30]=4.057306E-11; 
  EoverNdata[31]=1.919434E-19;  kdata[31]=6.286342E-11; 
  EoverNdata[32]=2.111378E-19;  kdata[32]=9.359894E-11; 
  EoverNdata[33]=2.322515E-19;  kdata[33]=1.344089E-10; 
  EoverNdata[34]=2.554767E-19;  kdata[34]=1.867662E-10; 
  EoverNdata[35]=2.810244E-19;  kdata[35]=2.518722E-10; 
  EoverNdata[36]=3.091268E-19;  kdata[36]=3.543204E-10; 
  EoverNdata[37]=3.400395E-19;  kdata[37]=5.277020E-10; 
  EoverNdata[38]=3.740434E-19;  kdata[38]=7.639723E-10; 
  EoverNdata[39]=4.114478E-19;  kdata[39]=1.078014E-09; 
  EoverNdata[40]=4.525926E-19;  kdata[40]=1.486210E-09; 
  EoverNdata[41]=4.978518E-19;  kdata[41]=2.006305E-09; 
  EoverNdata[42]=5.476370E-19;  kdata[42]=2.577740E-09; 
  EoverNdata[43]=6.024007E-19;  kdata[43]=3.252933E-09; 
  EoverNdata[44]=6.626408E-19;  kdata[44]=4.044664E-09; 
  EoverNdata[45]=7.289048E-19;  kdata[45]=4.962910E-09; 
  EoverNdata[46]=8.017953E-19;  kdata[46]=6.016556E-09; 
  EoverNdata[47]=8.819749E-19;  kdata[47]=7.149926E-09; 
  EoverNdata[48]=9.701723E-19;  kdata[48]=8.417356E-09; 
  EoverNdata[49]=1.067190E-18;  kdata[49]=9.502134E-09; 
  EoverNdata[50]=1.173909E-18;  kdata[50]=1.046326E-08; 
  EoverNdata[51]=1.291299E-18;  kdata[51]=1.142111E-08; 
  EoverNdata[52]=1.420429E-18;  kdata[52]=1.236778E-08; 
  EoverNdata[53]=1.562472E-18;  kdata[53]=1.329631E-08; 
  EoverNdata[54]=1.718719E-18;  kdata[54]=1.420078E-08; 
  EoverNdata[55]=1.890591E-18;  kdata[55]=1.507631E-08; 
  EoverNdata[56]=2.079651E-18;  kdata[56]=1.591900E-08; 
  EoverNdata[57]=2.287616E-18;  kdata[57]=1.672589E-08; 
  EoverNdata[58]=2.516377E-18;  kdata[58]=1.749487E-08; 
  EoverNdata[59]=2.768015E-18;  kdata[59]=1.822457E-08; 
  EoverNdata[60]=3.044816E-18;  kdata[60]=1.891432E-08; 
  EoverNdata[61]=3.349298E-18;  kdata[61]=1.956399E-08; 

  for (cnt=0; cnt<61; cnt++){
    if (EoverN>=EoverNdata[cnt] && EoverN<=EoverNdata[cnt+1]) {
      fact=(EoverN-EoverNdata[cnt])/(EoverNdata[cnt+1]-EoverNdata[cnt]);
      k=kdata[cnt]+fact*(kdata[cnt+1]-kdata[cnt]);
    }
  }
  return(k);
}


double _k1b(double EoverN){
  double kdata[62],EoverNdata[62];
  double k,fact;
  long cnt;
  
  EoverN=max(1.01E-20,EoverN); 
  EoverN=min(3.3e-18,EoverN); 
  k=0.0; /* to avoid compiler warning */
  EoverNdata[0]=1.000000E-20;  kdata[0]=1.258925E-37; 
  EoverNdata[1]=1.100000E-20;  kdata[1]=4.513833E-35; 
  EoverNdata[2]=1.210000E-20;  kdata[2]=9.481117E-33; 
  EoverNdata[3]=1.331000E-20;  kdata[3]=1.224769E-30; 
  EoverNdata[4]=1.464100E-20;  kdata[4]=1.017002E-28; 
  EoverNdata[5]=1.610510E-20;  kdata[5]=5.650812E-27; 
  EoverNdata[6]=1.771561E-20;  kdata[6]=2.179127E-25; 
  EoverNdata[7]=1.948717E-20;  kdata[7]=6.029161E-24; 
  EoverNdata[8]=2.143589E-20;  kdata[8]=1.233510E-22; 
  EoverNdata[9]=2.357948E-20;  kdata[9]=1.918037E-21; 
  EoverNdata[10]=2.593742E-20;  kdata[10]=2.323986E-20; 
  EoverNdata[11]=2.853117E-20;  kdata[11]=2.244509E-19; 
  EoverNdata[12]=3.138428E-20;  kdata[12]=1.763900E-18; 
  EoverNdata[13]=3.452271E-20;  kdata[13]=1.149293E-17; 
  EoverNdata[14]=3.797498E-20;  kdata[14]=6.315268E-17; 
  EoverNdata[15]=4.177248E-20;  kdata[15]=2.972237E-16; 
  EoverNdata[16]=4.594973E-20;  kdata[16]=1.215126E-15; 
  EoverNdata[17]=5.054470E-20;  kdata[17]=4.370839E-15; 
  EoverNdata[18]=5.559917E-20;  kdata[18]=1.399484E-14; 
  EoverNdata[19]=6.115909E-20;  kdata[19]=4.031118E-14; 
  EoverNdata[20]=6.727500E-20;  kdata[20]=1.054664E-13; 
  EoverNdata[21]=7.400250E-20;  kdata[21]=2.528314E-13; 
  EoverNdata[22]=8.140275E-20;  kdata[22]=5.597939E-13; 
  EoverNdata[23]=8.954302E-20;  kdata[23]=1.153038E-12; 
  EoverNdata[24]=9.849733E-20;  kdata[24]=2.223979E-12; 
  EoverNdata[25]=1.083471E-19;  kdata[25]=4.040940E-12; 
  EoverNdata[26]=1.191818E-19;  kdata[26]=6.954354E-12; 
  EoverNdata[27]=1.310999E-19;  kdata[27]=1.139192E-11; 
  EoverNdata[28]=1.442099E-19;  kdata[28]=1.784234E-11; 
  EoverNdata[29]=1.586309E-19;  kdata[29]=2.682827E-11; 
  EoverNdata[30]=1.744940E-19;  kdata[30]=3.887135E-11; 
  EoverNdata[31]=1.919434E-19;  kdata[31]=5.445365E-11; 
  EoverNdata[32]=2.111378E-19;  kdata[32]=7.398019E-11; 
  EoverNdata[33]=2.322515E-19;  kdata[33]=9.774733E-11; 
  EoverNdata[34]=2.554767E-19;  kdata[34]=1.259202E-10; 
  EoverNdata[35]=2.810244E-19;  kdata[35]=1.585210E-10; 
  EoverNdata[36]=3.091268E-19;  kdata[36]=4.016661E-10; 
  EoverNdata[37]=3.400395E-19;  kdata[37]=5.911514E-10; 
  EoverNdata[38]=3.740434E-19;  kdata[38]=8.459842E-10; 
  EoverNdata[39]=4.114478E-19;  kdata[39]=1.180412E-09; 
  EoverNdata[40]=4.525926E-19;  kdata[40]=1.609825E-09; 
  EoverNdata[41]=4.978518E-19;  kdata[41]=2.150605E-09; 
  EoverNdata[42]=5.476370E-19;  kdata[42]=2.829366E-09; 
  EoverNdata[43]=6.024007E-19;  kdata[43]=3.660356E-09; 
  EoverNdata[44]=6.626408E-19;  kdata[44]=4.663386E-09; 
  EoverNdata[45]=7.289048E-19;  kdata[45]=5.859519E-09; 
  EoverNdata[46]=8.017953E-19;  kdata[46]=7.265025E-09; 
  EoverNdata[47]=8.819749E-19;  kdata[47]=8.620912E-09; 
  EoverNdata[48]=9.701723E-19;  kdata[48]=1.013412E-08; 
  EoverNdata[49]=1.067190E-18;  kdata[49]=1.143479E-08; 
  EoverNdata[50]=1.173909E-18;  kdata[50]=1.259140E-08; 
  EoverNdata[51]=1.291299E-18;  kdata[51]=1.374407E-08; 
  EoverNdata[52]=1.420429E-18;  kdata[52]=1.488328E-08; 
  EoverNdata[53]=1.562472E-18;  kdata[53]=1.600066E-08; 
  EoverNdata[54]=1.718719E-18;  kdata[54]=1.708910E-08; 
  EoverNdata[55]=1.890591E-18;  kdata[55]=1.814271E-08; 
  EoverNdata[56]=2.079651E-18;  kdata[56]=1.915680E-08; 
  EoverNdata[57]=2.287616E-18;  kdata[57]=2.012780E-08; 
  EoverNdata[58]=2.516377E-18;  kdata[58]=2.105318E-08; 
  EoverNdata[59]=2.768015E-18;  kdata[59]=2.193130E-08; 
  EoverNdata[60]=3.044816E-18;  kdata[60]=2.276133E-08; 
  EoverNdata[61]=3.349298E-18;  kdata[61]=2.354314E-08; 

  for (cnt=0; cnt<61; cnt++){
    if (EoverN>=EoverNdata[cnt] && EoverN<=EoverNdata[cnt+1]) {
      fact=(EoverN-EoverNdata[cnt])/(EoverNdata[cnt+1]-EoverNdata[cnt]);
      k=kdata[cnt]+fact*(kdata[cnt+1]-kdata[cnt]);
    }
  }
  return(k);
}




void find_W(spec_t rhok, double T, double Te, double Tv, double Efieldstar, double Qbeam, spec_t W){
  long spec;
  double kf[nr];
  double Wp[nr];
  spec_t nk,calM;
  double vartheta,N;
  
  for (spec=0; spec<ns; spec++) W[spec]=0.0;
  /* find the gas temperature and the electron temperature */
  assert(Efieldstar!=0.0);     
  vartheta=Efieldstar/1.0e-20;
  /* find the mole fraction Xk in moles/cm^3 */
  for (spec=0; spec<ns; spec++){
    calM[spec]=_calM(spec);
    nk[spec]=rhok[spec]/calM[spec]*1e-6*calA;
  } 
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=nk[spec]; /* N is in 1/cm3 */
    
  /* set the reaction rates for each reaction */
  if (TOWNSEND) {
    kf[0]=pow(10.0,-8.3-36.5/vartheta);  /* 1a */
    kf[0]=_k1a(Efieldstar);
    kf[1]=pow(10.0,-8.8-28.1/vartheta);  /* 1b */
    kf[1]=_k1b(Efieldstar);
  } else {
    kf[0]=0.0;
    kf[1]=0.0;
  }
  kf[2]=2.0e-7*pow(300.0/Te,0.7);  /* 2a */
  kf[3]=2.8e-7*pow(300.0/Te,0.5);  /* 2b */
  kf[4]=2.0e-7*pow(300.0/T,0.5);   /* 3a */
  kf[5]=kf[4];                     /* 3b */
  kf[6]=2.0e-25*pow(300.0/T,2.5);  /* 4a */
  kf[7]=kf[6];                     /* 4b */
  kf[8]=kf[6];                     /* 4c */
  kf[9]=kf[6];                     /* 4d */
  if (ATTACHMENT) {
    kf[10]=1.4e-29*300.0/Te*exp(-600.0/T)*exp(700.0*(Te-T)/Te/T); /* 5a */
    kf[11]=1.07e-31*sqr(300.0/Te)*exp(-70.0/T)*exp(1500.0*(Te-T)/Te/T); /* 5b */
  } else {
    kf[10]=0.0;
    kf[11]=0.0;
  }
  kf[12]=8.6e-10*exp(-6030.0/T)*(1.0-exp(-1570.0/T)); /* 6 */
  kf[13]=1.8e11*Qbeam/N;  /* 7a */  /* Qbeam is in watts per m3, N in 1/cm3 */
  kf[14]=kf[13];            /* 7b */
  
  
  /* reaction 7b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  1,  0,  0,  0,  0,  0 },  /* 7b */
// MP=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 7b */
  Wp[14]=kf[14]*nk[2];
  W[0]+=+Wp[14];
  W[2]+=-Wp[14];
  W[6]+=+Wp[14];  
 
  /* reaction 7a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  0,  0,  0 },  /* 7a */
// MP=  {1,  0,  0,  0,  0,  1,  0,  0 },  /* 7a */
  Wp[13]=kf[13]*nk[1];
  W[0]+=+Wp[13];
  W[1]+=-Wp[13];
  W[5]+=+Wp[13];  
  
  /* reaction 6 */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  0,  0,  1 },  /* 6  */
// MP=  {1,  2,  0,  0,  0,  0,  0,  0 },  /* 6  */
  Wp[12]=kf[12]*nk[1]*nk[7];
  W[0]+=+Wp[12];
  W[1]+=+Wp[12];
  W[7]+=-Wp[12];

  /* reaction 5b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  1,  1,  0,  0,  0,  0,  0 },  /* 5b */
// MP=  {0,  0,  1,  0,  0,  0,  0,  1 },  /* 5b */
  Wp[11]=kf[11]*nk[0]*nk[1]*nk[2];
  W[0]+=-Wp[11];
  W[1]+=-Wp[11];
  W[7]+=+Wp[11];
  
  /* reaction 5a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  2,  0,  0,  0,  0,  0,  0 },  /* 5a */
// MP=  {0,  1,  0,  0,  0,  0,  0,  1 },  /* 5a */
  Wp[10]=kf[10]*nk[0]*sqr(nk[1]);
  W[0]+=-Wp[10];
  W[1]+=-Wp[10];
  W[7]+=+Wp[10];

  /* reaction 4d */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  1,  0,  1 },  /* 4d */
// MP=  {0,  3,  0,  0,  0,  0,  0,  0 },  /* 4d */
  Wp[9]=kf[9]*nk[1]*nk[5]*nk[7];
  W[1]+=+2.0*Wp[9];
  W[5]+=-Wp[9];
  W[7]+=-Wp[9];

  /* reaction 4c */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  0,  1,  1 },  /* 4c */
// MP=  {0,  2,  1,  0,  0,  0,  0,  0 },  /* 4c */
  Wp[8]=kf[8]*nk[1]*nk[6]*nk[7];
  W[1]+=+Wp[8];
  W[2]+=+Wp[8];
  W[6]+=-Wp[8];
  W[7]+=-Wp[8];

  /* reaction 4b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  1,  0,  0,  1,  0,  1 },  /* 4b */
// MP=  {0,  2,  1,  0,  0,  0,  0,  0 },  /* 4b */
  Wp[7]=kf[7]*nk[2]*nk[5]*nk[7];
  W[1]+=+2.0*Wp[7];
  W[5]+=-Wp[7];
  W[7]+=-Wp[7];

  /* reaction 4a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  1,  0,  0,  0,  1,  1 },  /* 4a */
// MP=  {0,  1,  2,  0,  0,  0,  0,  0 },  /* 4a */
  Wp[6]=kf[6]*nk[2]*nk[6]*nk[7];
  W[1]+=+Wp[6];
  W[2]+=+Wp[6];
  W[6]+=-Wp[6];
  W[7]+=-Wp[6];

  /* reaction 3b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 3b */
// MP=  {0,  2,  0,  0,  0,  0,  0,  0 },  /* 3b */
  Wp[5]=kf[5]*nk[5]*nk[7];
  W[1]+=+2.0*Wp[5];
  W[5]+=-Wp[5];
  W[7]+=-Wp[5];

  /* reaction 3a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  0,  0,  0,  0,  1,  1 },  /* 3a */
// MP=  {0,  1,  1,  0,  0,  0,  0,  0 },  /* 3a */
  Wp[4]=kf[4]*nk[6]*nk[7];
  W[1]+=+Wp[4];
  W[2]+=+Wp[4];
  W[6]+=-Wp[4];
  W[7]+=-Wp[4];

  /* reaction 2b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 2b */
// MP=  {0,  0,  0,  0,  2,  0,  0,  0 },  /* 2b */
  Wp[3]=kf[3]*nk[0]*nk[6];
  W[0]+=-Wp[3];
  W[4]+=+2.0*Wp[3];
  W[6]+=-Wp[3];
      
  /* reaction 2a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  0,  0,  0,  0,  1,  0,  0 },  /* 2a */
// MP=  {0,  0,  0,  2,  0,  0,  0,  0 },  /* 2a */
  Wp[2]=kf[2]*nk[0]*nk[5];
  W[0]+=-Wp[2];
  W[3]+=+2.0*Wp[2];
  W[5]+=-Wp[2];
  
  /* reaction 1b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  1,  0,  0,  0,  0,  0,  0 },  /* 1b */
// MP=  {2,  0,  0,  0,  0,  1,  0,  0 },  /* 1b */
  Wp[1]=kf[1]*nk[0]*nk[1];
  W[0]+=+Wp[1];
  W[1]+=-Wp[1];
  W[5]+=+Wp[1];
  
  /* reaction 1a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  0,  1,  0,  0,  0,  0,  0 },  /* 1a */
// MP=  {2,  0,  0,  0,  0,  0,  1,  0 },  /* 1a */
  Wp[0]=kf[0]*nk[0]*nk[2];
  W[0]+=+Wp[0];
  W[2]+=-Wp[0];
  W[6]+=+Wp[0];

  for (spec=0; spec<ns; spec++) W[spec]=W[spec]/calA*calM[spec]*1.0e6;
}



void find_dW_dx(spec_t rhok, double T, double Te, double Tv, double Efieldstar, double Qbeam,
              spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdEfieldstar, 
	      spec_t dWdQbeam){
  long  k,r,s,spec;  			/* counters */
  spec_t calM,nk;
  double dkdTe,dkdT,dkdE,dkdQb;
  double vartheta,N;
  double kf[nr];
  double dEfieldstar;
  
  /* first, initialize all derivatives to zero */
  for (s=0; s<ns; s++){
    dWdT[s]=0.0;
    dWdTe[s]=0.0;
    dWdTv[s]=0.0;
    dWdEfieldstar[s]=0.0;
    dWdQbeam[s]=0.0;
    for (k=0; k<ns; k++){
      dWdrhok[s][k]=0.0;
    }
  }
      
  assert(Efieldstar!=0.0);     
  vartheta=Efieldstar/1.0e-20;
  /* find the mole fraction Xk in moles/cm^3 */
  for (spec=0; spec<ns; spec++){
    calM[spec]=_calM(spec);
    nk[spec]=rhok[spec]/calM[spec]*1e-6*calA;
  } 
  N=0.0;
  for (spec=0; spec<ns; spec++) N+=nk[spec]; /* N is in 1/cm3 */
    
  /* set the reaction rates for each reaction */

  
  
  /* reaction 1a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  0,  1,  0,  0,  0,  0,  0 },  /* 1a */
// MP=  {2,  0,  0,  0,  0,  0,  1,  0 },  /* 1a */
  if (TOWNSEND) {
    kf[0]=pow(10.0,-8.3-36.5/vartheta);  /* 1a */
    kf[0]=_k1a(Efieldstar);
  } else {
    kf[0]=0.0;
  }
//  Wp[0]=kf[0]*nk[0]*nk[2];
//  W[0]+=+Wp[0];
//  W[2]+=-Wp[0];
//  W[6]+=+Wp[0];
  dWdrhok[0][0]+=+kf[0]*nk[2];
  dWdrhok[0][2]+=+kf[0]*nk[0];
  dWdrhok[2][0]+=-kf[0]*nk[2];
  dWdrhok[2][2]+=-kf[0]*nk[0];
  dWdrhok[6][0]+=+kf[0]*nk[2];
  dWdrhok[6][2]+=+kf[0]*nk[0];
  dkdE=kf[0]*log(10.0)*36.5/sqr(vartheta)*nk[0]*nk[2];
      
      /* to be used in conjunction with the derivatives of k1a and k1b wrt Efieldstar */
  dEfieldstar=max(1.0e-23,fabs(Efieldstar)/1000.0);

  dkdE=(_k1a(Efieldstar+dEfieldstar)-_k1a(Efieldstar))/dEfieldstar*nk[0]*nk[2]*1.0e-20;
  if (TOWNSEND){
    dWdEfieldstar[0]+=+dkdE;  
    dWdEfieldstar[2]+=-dkdE;  
    dWdEfieldstar[6]+=+dkdE;  
  }
  /* reaction 1b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  1,  0,  0,  0,  0,  0,  0 },  /* 1b */
// MP=  {2,  0,  0,  0,  0,  1,  0,  0 },  /* 1b */
  if (TOWNSEND) {
    kf[1]=pow(10.0,-8.8-28.1/vartheta);  /* 1b */
    kf[1]=_k1b(Efieldstar);
  } else {
    kf[1]=0.0;
  }
//  if (vartheta>1e-19) printf("vartheta=%E  %E  %E\n",vartheta,pow(10.0,-8.8-28.1/vartheta),_k1b(Efieldstar));
  
//  Wp[1]=kf[1]*nk[0]*nk[1];
//  W[0]+=+Wp[1];
//  W[1]+=-Wp[1];
//  W[5]+=+Wp[1];  
  dWdrhok[0][0]+=+kf[1]*nk[1];
  dWdrhok[0][1]+=+kf[1]*nk[0];
  dWdrhok[1][0]+=-kf[1]*nk[1];
  dWdrhok[1][1]+=-kf[1]*nk[0];
  dWdrhok[5][0]+=+kf[1]*nk[1];
  dWdrhok[5][1]+=+kf[1]*nk[0];
  dkdE=kf[1]*log(10.0)*28.1/sqr(vartheta)*nk[0]*nk[1];
  dkdE=(_k1b(Efieldstar+dEfieldstar)-_k1b(Efieldstar))/dEfieldstar*nk[0]*nk[1]*1e-20;
  //if (vartheta>1e-19) printf("vartheta=%E  %E  %E\n",vartheta,kf[1]*log(10.0)*28.1/sqr(vartheta)*nk[0]*nk[1],(_k1b(Efieldstar+dEfieldstar)-_k1b(Efieldstar))/dEfieldstar*nk[0]*nk[1]*1e-20);
  if (TOWNSEND) { 
    dWdEfieldstar[0]+=+dkdE;      
    dWdEfieldstar[1]+=-dkdE;      
    dWdEfieldstar[5]+=+dkdE;      
  }
  /* reaction 2a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  0,  0,  0,  0,  1,  0,  0 },  /* 2a */
// MP=  {0,  0,  0,  2,  0,  0,  0,  0 },  /* 2a */
  kf[2]=2.0e-7*pow(300.0/Te,0.7);  /* 2a */
//  Wp[2]=kf[2]*nk[0]*nk[5];
//  W[0]+=-Wp[2];
//  W[3]+=+2.0*Wp[2];
//  W[5]+=-Wp[2];
  dWdrhok[0][0]+=-kf[2]*nk[5];
  dWdrhok[0][5]+=-kf[2]*nk[0];
  dWdrhok[3][0]+=+2.0*kf[2]*nk[5];
  dWdrhok[3][5]+=+2.0*kf[2]*nk[0];
  dWdrhok[5][0]+=-kf[2]*nk[5];
  dWdrhok[5][5]+=-kf[2]*nk[0];
  dkdTe=-0.7*2.0e-7*pow(Te/300.0,-1.7)/300.0*nk[0]*nk[5];
  dWdTe[0]+=-dkdTe;
  dWdTe[3]+=+2.0*dkdTe;
  dWdTe[5]+=-dkdTe;
    
  /* reaction 2b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 2b */
// MP=  {0,  0,  0,  0,  2,  0,  0,  0 },  /* 2b */
  kf[3]=2.8e-7*pow(300.0/Te,0.5);  /* 2b */
//  Wp[3]=kf[3]*nk[0]*nk[6];
//  W[0]+=-Wp[3];
//  W[4]+=+2.0*Wp[3];
//  W[6]+=-Wp[3];
  dWdrhok[0][0]+=-kf[3]*nk[6];
  dWdrhok[0][6]+=-kf[3]*nk[0];
  dWdrhok[4][0]+=+2.0*kf[3]*nk[6];
  dWdrhok[4][6]+=+2.0*kf[3]*nk[0];
  dWdrhok[6][0]+=-kf[3]*nk[6];
  dWdrhok[6][6]+=-kf[3]*nk[0];
  dkdTe=-0.5*2.8e-7*pow(Te/300.0,-1.5)/300.0*nk[0]*nk[6];
  dWdTe[0]+=-dkdTe;
  dWdTe[4]+=+2.0*dkdTe;
  dWdTe[6]+=-dkdTe;
  
  /* reaction 3a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  0,  0,  0,  0,  1,  1 },  /* 3a */
// MP=  {0,  1,  1,  0,  0,  0,  0,  0 },  /* 3a */
  kf[4]=2.0e-7*pow(300.0/T,0.5);   /* 3a */
//  Wp[4]=kf[4]*nk[6]*nk[7];
//  W[1]+=+Wp[4];
//  W[2]+=+Wp[4];
//  W[6]+=-Wp[4];
//  W[7]+=-Wp[4];
  dWdrhok[1][6]+=+kf[4]*nk[7];
  dWdrhok[1][7]+=+kf[4]*nk[6];
  dWdrhok[2][6]+=+kf[4]*nk[7];
  dWdrhok[2][7]+=+kf[4]*nk[6];
  dWdrhok[6][6]+=-kf[4]*nk[7];
  dWdrhok[6][7]+=-kf[4]*nk[6];
  dWdrhok[7][6]+=-kf[4]*nk[7];
  dWdrhok[7][7]+=-kf[4]*nk[6];
  dkdT=-0.5*2.0e-7*pow(T/300.0,-1.5)/300.0*nk[6]*nk[7];
  dWdT[1]+=+dkdT;
  dWdT[2]+=+dkdT;
  dWdT[6]+=-dkdT;
  dWdT[7]+=-dkdT;
  
  
  /* reaction 3b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  0,  0,  0,  1,  0,  1 },  /* 3b */
// MP=  {0,  2,  0,  0,  0,  0,  0,  0 },  /* 3b */
  kf[5]=kf[4];                     /* 3b */
//  Wp[5]=kf[5]*nk[5]*nk[7];
//  W[1]+=+2.0*Wp[5];
//  W[5]+=-Wp[5];
//  W[7]+=-Wp[5];
  dWdrhok[1][5]+=+2.0*kf[5]*nk[7];
  dWdrhok[1][7]+=+2.0*kf[5]*nk[5];
  dWdrhok[5][5]+=-kf[5]*nk[7];
  dWdrhok[5][7]+=-kf[5]*nk[5];
  dWdrhok[7][5]+=-kf[5]*nk[7];
  dWdrhok[7][7]+=-kf[5]*nk[5];
  dkdT=-0.5*2.0e-7*pow(T/300.0,-1.5)/300.0*nk[5]*nk[7];
  dWdT[1]+=+2.0*dkdT;
  dWdT[5]+=-dkdT;
  dWdT[7]+=-dkdT;
  
  /* reaction 4a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  1,  0,  0,  0,  1,  1 },  /* 4a */
// MP=  {0,  1,  2,  0,  0,  0,  0,  0 },  /* 4a */
  kf[6]=2.0e-25*pow(300.0/T,2.5);  /* 4a */
//  Wp[6]=kf[6]*nk[2]*nk[6]*nk[7];
//  W[1]+=+Wp[6];
//  W[2]+=+Wp[6];
//  W[6]+=-Wp[6];
//  W[7]+=-Wp[6];
  dWdrhok[1][2]+=+kf[6]*nk[6]*nk[7];
  dWdrhok[1][6]+=+kf[6]*nk[2]*nk[7];
  dWdrhok[1][7]+=+kf[6]*nk[2]*nk[6];
  dWdrhok[2][2]+=+kf[6]*nk[6]*nk[7];
  dWdrhok[2][6]+=+kf[6]*nk[2]*nk[7];
  dWdrhok[2][7]+=+kf[6]*nk[2]*nk[6];
  dWdrhok[6][2]+=-kf[6]*nk[6]*nk[7];
  dWdrhok[6][6]+=-kf[6]*nk[2]*nk[7];
  dWdrhok[6][7]+=-kf[6]*nk[2]*nk[6];
  dWdrhok[7][2]+=-kf[6]*nk[6]*nk[7];
  dWdrhok[7][6]+=-kf[6]*nk[2]*nk[7];
  dWdrhok[7][7]+=-kf[6]*nk[2]*nk[6];
  dkdT=-2.5*kf[6]/T*nk[2]*nk[6]*nk[7];
  dWdT[1]+=+dkdT;
  dWdT[2]+=+dkdT;
  dWdT[6]+=-dkdT;
  dWdT[7]+=-dkdT;
  
  /* reaction 4b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  1,  0,  0,  1,  0,  1 },  /* 4b */
// MP=  {0,  2,  1,  0,  0,  0,  0,  0 },  /* 4b */
  kf[7]=kf[6];                     /* 4b */
//  Wp[7]=kf[7]*nk[2]*nk[5]*nk[7];
//  W[1]+=+2.0*Wp[7];
//  W[5]+=-Wp[7];
//  W[7]+=-Wp[7];
  dWdrhok[1][2]+=+2.0*kf[7]*nk[5]*nk[7];
  dWdrhok[1][5]+=+2.0*kf[7]*nk[2]*nk[7];
  dWdrhok[1][7]+=+2.0*kf[7]*nk[2]*nk[5];
  dWdrhok[5][2]+=-kf[7]*nk[5]*nk[7];
  dWdrhok[5][5]+=-kf[7]*nk[2]*nk[7];
  dWdrhok[5][7]+=-kf[7]*nk[2]*nk[5];
  dWdrhok[7][2]+=-kf[7]*nk[5]*nk[7];
  dWdrhok[7][5]+=-kf[7]*nk[2]*nk[7];
  dWdrhok[7][7]+=-kf[7]*nk[2]*nk[5];
  dkdT=-2.5*kf[7]/T*nk[2]*nk[5]*nk[7];
  dWdT[1]+=+2.0*dkdT;   /* 2.0e-25*pow(300.0/T,2.5) */
  dWdT[5]+=-dkdT;
  dWdT[7]+=-dkdT;
  
  /* reaction 4c */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  0,  1,  1 },  /* 4c */
// MP=  {0,  2,  1,  0,  0,  0,  0,  0 },  /* 4c */
  kf[8]=kf[6];                     /* 4c */
//  Wp[8]=kf[8]*nk[1]*nk[6]*nk[7];
//  W[1]+=+Wp[8];
//  W[2]+=+Wp[8];
//  W[6]+=-Wp[8];
//  W[7]+=-Wp[8];
  dWdrhok[1][1]+=+kf[8]*nk[6]*nk[7];
  dWdrhok[1][6]+=+kf[8]*nk[1]*nk[7];
  dWdrhok[1][7]+=+kf[8]*nk[1]*nk[6];
  dWdrhok[2][1]+=+kf[8]*nk[6]*nk[7];
  dWdrhok[2][6]+=+kf[8]*nk[1]*nk[7];
  dWdrhok[2][7]+=+kf[8]*nk[1]*nk[6];
  dWdrhok[6][1]+=-kf[8]*nk[6]*nk[7];
  dWdrhok[6][6]+=-kf[8]*nk[1]*nk[7];
  dWdrhok[6][7]+=-kf[8]*nk[1]*nk[6];
  dWdrhok[7][1]+=-kf[8]*nk[6]*nk[7];
  dWdrhok[7][6]+=-kf[8]*nk[1]*nk[7];
  dWdrhok[7][7]+=-kf[8]*nk[1]*nk[6];
  dkdT=-2.5*kf[8]/T*nk[1]*nk[6]*nk[7];
  dWdT[1]+=+dkdT;  
  dWdT[2]+=+dkdT;
  dWdT[6]+=-dkdT;
  dWdT[7]+=-dkdT;
  
  
  /* reaction 4d */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  1,  0,  1 },  /* 4d */
// MP=  {0,  3,  0,  0,  0,  0,  0,  0 },  /* 4d */
  kf[9]=kf[6];                     /* 4d */
//  Wp[9]=kf[9]*nk[1]*nk[5]*nk[7];
//  W[1]+=+2.0*Wp[9];
//  W[5]+=-Wp[9];
//  W[7]+=-Wp[9];
  dWdrhok[1][1]+=+2.0*kf[9]*nk[5]*nk[7];
  dWdrhok[1][5]+=+2.0*kf[9]*nk[1]*nk[7];
  dWdrhok[1][7]+=+2.0*kf[9]*nk[1]*nk[5];
  dWdrhok[5][1]+=-kf[9]*nk[5]*nk[7];
  dWdrhok[5][5]+=-kf[9]*nk[1]*nk[7];
  dWdrhok[5][7]+=-kf[9]*nk[1]*nk[5];
  dWdrhok[7][1]+=-kf[9]*nk[5]*nk[7];
  dWdrhok[7][5]+=-kf[9]*nk[1]*nk[7];
  dWdrhok[7][7]+=-kf[9]*nk[1]*nk[5];
  dkdT=-2.5*kf[9]/T*nk[1]*nk[5]*nk[7];
  dWdT[1]+=+2.0*dkdT;  /* 2.0e-25*pow(300.0/T,2.5) */
  dWdT[5]+=-dkdT;
  dWdT[7]+=-dkdT;


  /* reaction 5a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  2,  0,  0,  0,  0,  0,  0 },  /* 5a */
// MP=  {0,  1,  0,  0,  0,  0,  0,  1 },  /* 5a */
  if (ATTACHMENT){
    kf[10]=1.4e-29*300.0/Te*exp(-600.0/T)*exp(700.0*(Te-T)/Te/T); /* 5a */
  } else {
    kf[10]=0.0;
  }
//  Wp[10]=kf[10]*nk[0]*sqr(nk[1]);
//  W[0]+=-Wp[10];
//  W[1]+=-Wp[10];
//  W[7]+=+Wp[10];
  dWdrhok[0][0]+=-kf[10]*sqr(nk[1]);
  dWdrhok[0][1]+=-kf[10]*nk[0]*2.0*nk[1];
  dWdrhok[1][0]+=-kf[10]*sqr(nk[1]);
  dWdrhok[1][1]+=-kf[10]*nk[0]*2.0*nk[1];
  dWdrhok[7][0]+=+kf[10]*sqr(nk[1]);
  dWdrhok[7][1]+=+kf[10]*nk[0]*2.0*nk[1];
  dkdTe=(-1.0/Te+700.0/sqr(Te))*kf[10]*nk[0]*sqr(nk[1]);
  dkdT= (600.0/sqr(T)-700.0/Te/T-700.0*(Te-T)/Te/sqr(T))*kf[10]*nk[0]*sqr(nk[1]);
  dWdTe[0]+=-dkdTe;
  dWdTe[1]+=-dkdTe;
  dWdTe[7]+=+dkdTe;
  dWdT[0]+=-dkdT;
  dWdT[1]+=-dkdT;
  dWdT[7]+=+dkdT;
  
  /* reaction 5b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {1,  1,  1,  0,  0,  0,  0,  0 },  /* 5b */
// MP=  {0,  0,  1,  0,  0,  0,  0,  1 },  /* 5b */
  if (ATTACHMENT) {
    kf[11]=1.07e-31*sqr(300.0/Te)*exp(-70.0/T)*exp(1500.0*(Te-T)/Te/T); /* 5b */
  } else {
    kf[11]=0.0;
  }
//  Wp[11]=kf[11]*nk[0]*nk[1]*nk[2];
//  W[0]+=-Wp[11];
//  W[1]+=-Wp[11];
//  W[7]+=+Wp[11];
  dWdrhok[0][0]+=-kf[11]*nk[1]*nk[2];
  dWdrhok[0][1]+=-kf[11]*nk[0]*nk[2];
  dWdrhok[0][2]+=-kf[11]*nk[0]*nk[1];
  dWdrhok[1][0]+=-kf[11]*nk[1]*nk[2];
  dWdrhok[1][1]+=-kf[11]*nk[0]*nk[2];
  dWdrhok[1][2]+=-kf[11]*nk[0]*nk[1];
  dWdrhok[7][0]+=+kf[11]*nk[1]*nk[2];
  dWdrhok[7][1]+=+kf[11]*nk[0]*nk[2];
  dWdrhok[7][2]+=+kf[11]*nk[0]*nk[1];
  dkdTe=(-2.0/Te+1500.0/sqr(Te))*kf[11]*nk[0]*nk[1]*nk[2];
  dkdT=(70.0/sqr(T)-1500.0/sqr(T))*kf[11]*nk[0]*nk[1]*nk[2];
  dWdTe[0]+=-dkdTe;
  dWdTe[1]+=-dkdTe;
  dWdTe[7]+=+dkdTe;
  dWdT[0]+=-dkdT;
  dWdT[1]+=-dkdT;
  dWdT[7]+=+dkdT;
      
  /* reaction 6 */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  0,  0,  1 },  /* 6  */
// MP=  {1,  2,  0,  0,  0,  0,  0,  0 },  /* 6  */
  kf[12]=8.6e-10*exp(-6030.0/T)*(1.0-exp(-1570.0/T)); /* 6 */
//  Wp[12]=kf[12]*nk[1]*nk[7];
//  W[0]+=+Wp[12];
//  W[1]+=+Wp[12];
//  W[7]+=-Wp[12];
  dWdrhok[0][1]+=+kf[12]*nk[7];
  dWdrhok[0][7]+=+kf[12]*nk[1];
  dWdrhok[1][1]+=+kf[12]*nk[7];
  dWdrhok[1][7]+=+kf[12]*nk[1];
  dWdrhok[7][1]+=-kf[12]*nk[7];
  dWdrhok[7][7]+=-kf[12]*nk[1];
      
  dkdT=(
        kf[12]*6030.0/sqr(T)-8.6e-10*exp(-6030.0/T)*exp(-1570.0/T)*1570.0/sqr(T)
       )*nk[1]*nk[7];
      
  dWdT[0]+=+dkdT;
  dWdT[1]+=+dkdT;
  dWdT[7]+=-dkdT;
  
  /* reaction 7a */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  1,  0,  0,  0,  0,  0,  0 },  /* 7a */
// MP=  {1,  0,  0,  0,  0,  1,  0,  0 },  /* 7a */
  kf[13]=1.8e11*Qbeam/N;  /* 7a */  /* Qbeam is in watts per m3, N in 1/cm3 */
//  Wp[13]=kf[13]*nk[1];
//  W[0]+=+Wp[13];
//  W[1]+=-Wp[13];
//  W[5]+=+Wp[13];  
  dWdrhok[0][1]+=+kf[13];
  dWdrhok[1][1]+=-kf[13];
  dWdrhok[5][1]+=+kf[13];
  for (s=0; s<ns; s++){
    dWdrhok[0][s]+=-kf[13]*nk[1]/N;
    dWdrhok[1][s]+=+kf[13]*nk[1]/N;
    dWdrhok[5][s]+=-kf[13]*nk[1]/N;
  }
  dkdQb=1.8e11/N*nk[1];
  dWdQbeam[0]+=+dkdQb; 
  dWdQbeam[1]+=-dkdQb; 
  dWdQbeam[5]+=+dkdQb; 
  
  /* reaction 7b */
//    {/*e-  O2  N2  O   N   O2+ N2+ O2-   */
// MR=  {0,  0,  1,  0,  0,  0,  0,  0 },  /* 7b */
// MP=  {1,  0,  0,  0,  0,  0,  1,  0 },  /* 7b */
  kf[14]=kf[13];            /* 7b */
//  Wp[14]=kf[14]*nk[2];
//  W[0]+=+Wp[14];
//  W[2]+=-Wp[14];
//  W[6]+=+Wp[14];  
  dWdrhok[0][2]+=+kf[14];
  dWdrhok[2][2]+=-kf[14];
  dWdrhok[6][2]+=+kf[14];
  for (s=0; s<ns; s++){
    dWdrhok[0][s]+=-kf[14]*nk[2]/N;
    dWdrhok[2][s]+=+kf[14]*nk[2]/N;
    dWdrhok[6][s]+=-kf[14]*nk[2]/N;
  }
  dkdQb=1.8e11/N*nk[2];
  dWdQbeam[0]+=+dkdQb; 
  dWdQbeam[2]+=-dkdQb; 
  dWdQbeam[6]+=+dkdQb; 
  
  for (s=0; s<ns; s++) {
    for (r=0; r<ns; r++) {
      dWdrhok[s][r]=dWdrhok[s][r]*calM[s]/calM[r];
    }
    dWdT[s]=dWdT[s]/calA*calM[s]*1.0e6;
    dWdTe[s]=dWdTe[s]/calA*calM[s]*1.0e6;
    dWdEfieldstar[s]=dWdEfieldstar[s]/calA*calM[s]*1.0e6*1.0e20;
    dWdQbeam[s]=dWdQbeam[s]/calA*calM[s]*1.0e6;
  }
              
}




void test_derivative(spec_t rhok, double T, double Te, double Tv, double Efieldstar, double Qbeam){
  long s; 
  spec_t dWdT,dWdTe,dWdTv,dWdEfieldstar,dWdQbeam;
  spec2_t dWdrhok;   
  spec_t dWdT2,W,W2;
  double dQbeam;
 
  find_W(rhok, T, Te, Tv, Efieldstar, Qbeam, W);      
  find_dW_dx(rhok, T, Te, Tv, Efieldstar, Qbeam,
              dWdrhok, dWdT, dWdTe, dWdTv, dWdEfieldstar, dWdQbeam);
   
   /* find dWdQbeam numerically */
   dQbeam=Qbeam/1000000.0+1.0e-2;
   find_W(rhok, T, Te,Tv,Efieldstar,Qbeam+dQbeam, W2);
   for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dQbeam;
   for (s=0; s<ns; s++){  
     wfprintf(stdout,"%15.15E  %15.15E\n",dWdQbeam[s],dWdT2[s]);
   }
      
  /* find dWdrhok numerically */
  /*for (s=0; s<ns; s++){
    for (r=0; r<ns; r++) rhok2[r]=rhok[r];
    rhok2[s]+=rhok[s]/1000000.0+1.0e-10;
    find_W(rhok2, T, Te,Tv,Efieldstar,Qbeam, W2);
    for (r=0; r<ns; r++) dWdrhok2[r][s]=(W2[r]-W[r])/(rhok[s]/1000000.0+1.0e-10);
  }
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdrhok[s][7],dWdrhok2[s][7]);
  }*/

  /* find dWdT numerically */
  /* dT=T/1000000000.0;
  find_W(rhok, T+dT, Te,Tv,Efieldstar,Qbeam, W2);
  for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dT;
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdT[s],dWdT2[s]);
  } */
  
  /* find dWdTe numerically */
  /*dT=Te/10000000.0;
  find_W(rhok, T, Te+dT,Tv,Efieldstar,Qbeam, W2);
  for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dT;
  for (s=0; s<ns; s++){  
    wfprintf(stdout,"%15.15E  %15.15E\n",dWdTe[s],dWdT2[s]);
  }*/

   /* find dWdEfieldstar numerically */
/*   dEfieldstar=Efieldstar/10000000.0+1.0e-30;
   find_W(rhok, T, Te,Tv,Efieldstar+dEfieldstar,Qbeam, W2);
   for (s=0; s<ns; s++) dWdT2[s]=(W2[s]-W[s])/dEfieldstar;
   for (s=0; s<ns; s++){  
     wfprintf(stdout,"%15.15E  %15.15E\n",dWdEfieldstar[s],dWdT2[s]);
   }*/
    
  exit(1);  
}

