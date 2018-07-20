#include <stdio.h>
#include <stdlib.h>
#include <exm.h>

#include <math.h>
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))
#define calR           8.314472e0
#define calA 6.02214199E23  
#define echarge 1.602176462e-19
#define emass 9.10938188E-31
#define imass 5e-26
#define epsilon0 8.854187817E-12
#define pi  3.14159265358979323846
#define kB 1.3807E-23
#define sqr(a)   ((a)*(a))
#define Ipot 1.602e-18
//#define Ipot 0.0



double _fq_k1a(double EoverN){
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


double _fq_k1b(double EoverN){
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



double _fq(double EoverN){
  double fq;
  fq=0.765*_fq_k1a(EoverN)+0.235*_fq_k1b(EoverN);
  fq=fq*1e-6; /* change from cm3/2 to m3/s */
  return(fq);
}


double _Teprime(double EeffoverN){
  double Te,fact;
  double EeffoverNdata[18];
  double Tedata[18];
  long cnt;
  
  Te=0.0;
  EeffoverNdata[0]=0.003E-20;   Tedata[0]=286.8;
  EeffoverNdata[1]=0.005E-20;   Tedata[1]=317.5;
  EeffoverNdata[2]=0.007E-20;   Tedata[2]=358.5;
  EeffoverNdata[3]=0.01E-20;    Tedata[3]=430.2;
  EeffoverNdata[4]=0.03E-20;    Tedata[4]=1024.1;
  EeffoverNdata[5]=0.05E-20;    Tedata[5]=1638.6;
  EeffoverNdata[6]=0.07E-20;    Tedata[6]=2048.3;
  EeffoverNdata[7]=0.1E-20;     Tedata[7]=2866.4;
  EeffoverNdata[8]=0.3E-20;     Tedata[8]=6233.0;
  EeffoverNdata[9]=0.5E-20;     Tedata[9]=7906.5;
  EeffoverNdata[10]=0.7E-20;    Tedata[10]=9448.8;
  EeffoverNdata[11]=1.0E-20;    Tedata[11]=11256.9;
  EeffoverNdata[12]=3.0E-20;    Tedata[12]=17198.6;
  EeffoverNdata[13]=5.0E-20;    Tedata[13]=20337.8;
  EeffoverNdata[14]=7.0E-20;    Tedata[14]=21498.3;
  EeffoverNdata[15]=10.0E-20;   Tedata[15]=26140.3;
  EeffoverNdata[16]=20.0E-20;   Tedata[16]=46874.9;
  EeffoverNdata[17]=30.0E-20;   Tedata[17]=57904.3;
  if (EeffoverN>=EeffoverNdata[17]) Te=Tedata[17];
  if (EeffoverN<=EeffoverNdata[0]) Te=Tedata[0];
  for (cnt=0; cnt<17; cnt++){
    if (EeffoverN>=EeffoverNdata[cnt] && EeffoverN<=EeffoverNdata[cnt+1]) {
      fact=(EeffoverN-EeffoverNdata[cnt])/(EeffoverNdata[cnt+1]-EeffoverNdata[cnt]);
      Te=Tedata[cnt]+fact*(Tedata[cnt+1]-Tedata[cnt]);
    }
  }
  return(Te);
}

double _mueNnprime(double EeffoverNn){
  double vdrO2data[23],vdrN2data[23];
  double EeffoverNdata[23];
  double cO2,cN2,muenNn,fact,vdrO2,vdrN2,vdr,EeffoverN;
  long cnt;
 
  EeffoverN=max(0.00301E-20,EeffoverNn); 
  EeffoverN=min(9999.99E-20,EeffoverN); 
  cO2=0.235;
  cN2=0.765;
  vdrO2=0.0; /* to avoid compiler warning */
  vdrN2=0.0; /* idem */
  EeffoverNdata[0]=0.003E-20;  vdrN2data[0]=1.1E3;   vdrO2data[0]=2.6E3;
  EeffoverNdata[1]=0.005E-20;  vdrN2data[1]=1.6E3;   vdrO2data[1]=3.7E3;
  EeffoverNdata[2]=0.007E-20;  vdrN2data[2]=2.0E3;   vdrO2data[2]=4.1E3;
  EeffoverNdata[3]=0.01E-20;   vdrN2data[3]=2.4E3;   vdrO2data[3]=4.3E3;
  EeffoverNdata[4]=0.03E-20;   vdrN2data[4]=3.1E3;   vdrO2data[4]=5.4E3;
  EeffoverNdata[5]=0.05E-20;   vdrN2data[5]=3.5E3;   vdrO2data[5]=7.3E3;
  EeffoverNdata[6]=0.07E-20;   vdrN2data[6]=3.9E3;   vdrO2data[6]=9.2E3;
  EeffoverNdata[7]=0.10E-20;   vdrN2data[7]=4.4E3;   vdrO2data[7]=10.0E3;
  EeffoverNdata[8]=0.30E-20;   vdrN2data[8]=7.5E3;   vdrO2data[8]=22.0E3;
  EeffoverNdata[9]=0.50E-20;   vdrN2data[9]=11.0E3;  vdrO2data[9]=26.0E3;
  EeffoverNdata[10]=0.70E-20;  vdrN2data[10]=14.0E3; vdrO2data[10]=28.0E3;
  EeffoverNdata[11]=1.0E-20;   vdrN2data[11]=18.0E3; vdrO2data[11]=31.0E3;
  EeffoverNdata[12]=3.0E-20;   vdrN2data[12]=42.0E3; vdrO2data[12]=72.0E3;
  EeffoverNdata[13]=5.0E-20;   vdrN2data[13]=61.0E3; vdrO2data[13]=110.0E3;
  EeffoverNdata[14]=7.0E-20;   vdrN2data[14]=79.0E3; vdrO2data[14]=120.0E3;
  EeffoverNdata[15]=10.0E-20;  vdrN2data[15]=100.0E3; vdrO2data[15]=160.0E3;
  EeffoverNdata[16]=20.0E-20;  vdrN2data[16]=200.0E3; vdrO2data[16]=260.0E3;
  EeffoverNdata[17]=30.0E-20;  vdrN2data[17]=290.0E3; vdrO2data[17]=330.0E3;
  EeffoverNdata[18]=50.0E-20;  vdrN2data[18]=420.0E3; vdrO2data[18]=450.0E3;
  EeffoverNdata[19]=80.0E-20;  vdrN2data[19]=530.0E3; vdrO2data[19]=640.0E3;
  EeffoverNdata[20]=100.0E-20; vdrN2data[20]=590.0E3; vdrO2data[20]=710.0E3;
  EeffoverNdata[21]=1000.0E-20; vdrN2data[21]=2500.0E3; vdrO2data[21]=3500.0E3; //wrong
  EeffoverNdata[22]=10000.0E-20; vdrN2data[22]=12500.0E3; vdrO2data[22]=17500.0E3; //wrong

  for (cnt=0; cnt<22; cnt++){
    if (EeffoverN>=EeffoverNdata[cnt] && EeffoverN<=EeffoverNdata[cnt+1]) {
      fact=(EeffoverN-EeffoverNdata[cnt])/(EeffoverNdata[cnt+1]-EeffoverNdata[cnt]);
      vdrN2=vdrN2data[cnt]+fact*(vdrN2data[cnt+1]-vdrN2data[cnt]);
      vdrO2=vdrO2data[cnt]+fact*(vdrO2data[cnt+1]-vdrO2data[cnt]);
    }
  }
  vdr=(cO2*vdrO2+cN2*vdrN2)/(cO2+cN2);
  muenNn=vdr/EeffoverN;
  
  return(muenNn);
}

double _delta(double EoverN){  
  double delta; 
  delta=2.0/3.0*emass/kB/_Teprime(EoverN)*(
    sqr(_mueNnprime(EoverN))*sqr(EoverN)-Ipot/echarge*_mueNnprime(EoverN)*_fq(EoverN)
  );
  delta=max(delta,2.0*emass/imass);
  return(delta);
}

double _deltafromTe(double Te){
  long cnt;
  double delta,fact;
  double Tedata[51],deltadata[51];
  
  Tedata[0]=2.868000E+02;  deltadata[0]= 3.198217E-04;
  Tedata[1]=2.957518E+02;  deltadata[1]= 3.997107E-04;
  Tedata[2]=3.067522E+02;  deltadata[2]= 5.009238E-04;
  Tedata[3]=3.207755E+02;  deltadata[3]= 6.194431E-04;
  Tedata[4]=3.419306E+02;  deltadata[4]= 6.994563E-04;
  Tedata[5]=3.687790E+02;  deltadata[5]= 7.719741E-04;
  Tedata[6]=4.042948E+02;  deltadata[6]= 8.042780E-04;
  Tedata[7]=4.509663E+02;  deltadata[7]= 8.057628E-04;
  Tedata[8]=5.145096E+02;  deltadata[8]= 7.486142E-04;
  Tedata[9]=5.907615E+02;  deltadata[9]= 6.976825E-04;
  Tedata[10]=6.822638E+02; deltadata[10]=6.533609E-04;
  Tedata[11]=7.920666E+02; deltadata[11]=6.158861E-04;
  Tedata[12]=9.238299E+02; deltadata[12]=5.853838E-04;
  Tedata[13]=1.083952E+03; deltadata[13]=5.596592E-04;
  Tedata[14]=1.280273E+03; deltadata[14]=5.371722E-04;
  Tedata[15]=1.515857E+03; deltadata[15]=5.223037E-04;
  Tedata[16]=1.745248E+03; deltadata[16]=5.307070E-04;
  Tedata[17]=1.971427E+03; deltadata[17]=5.587345E-04;
  Tedata[18]=2.307279E+03; deltadata[18]=5.407750E-04;
  Tedata[19]=2.740855E+03; deltadata[19]=5.083820E-04;
  Tedata[20]=3.110065E+03; deltadata[20]=5.248225E-04;
  Tedata[21]=3.495459E+03; deltadata[21]=5.625134E-04;
  Tedata[22]=3.957930E+03; deltadata[22]=6.084156E-04;
  Tedata[23]=4.512896E+03; deltadata[23]=6.641662E-04;
  Tedata[24]=5.178855E+03; deltadata[24]=7.317160E-04;
  Tedata[25]=5.978007E+03; deltadata[25]=8.133948E-04;
  Tedata[26]=6.582945E+03; deltadata[26]=9.090121E-04;
  Tedata[27]=7.154984E+03; deltadata[27]=1.023063E-03;
  Tedata[28]=7.841430E+03; deltadata[28]=1.160601E-03;
  Tedata[29]=8.605688E+03; deltadata[29]=1.272458E-03;
  Tedata[30]=9.501849E+03; deltadata[30]=1.401560E-03;
  Tedata[31]=1.035624E+04; deltadata[31]=1.562332E-03;
  Tedata[32]=1.131832E+04; deltadata[32]=1.770451E-03;
  Tedata[33]=1.192478E+04; deltadata[33]=2.160441E-03;
  Tedata[34]=1.265252E+04; deltadata[34]=2.654011E-03;
  Tedata[35]=1.352582E+04; deltadata[35]=3.277121E-03;
  Tedata[36]=1.457377E+04; deltadata[36]=4.061173E-03;
  Tedata[37]=1.583131E+04; deltadata[37]=5.043836E-03;
  Tedata[38]=1.727350E+04; deltadata[38]=6.266917E-03;
  Tedata[39]=1.823024E+04; deltadata[39]=7.773400E-03;
  Tedata[40]=1.937833E+04; deltadata[40]=9.691279E-03;
  Tedata[41]=2.049241E+04; deltadata[41]=1.196490E-02;
  Tedata[42]=2.110358E+04; deltadata[42]=1.441023E-02;
  Tedata[43]=2.240148E+04; deltadata[43]=1.719452E-02;
  Tedata[44]=2.474838E+04; deltadata[44]=2.012997E-02;
  Tedata[45]=2.804897E+04; deltadata[45]=2.379294E-02;
  Tedata[46]=3.257763E+04; deltadata[46]=2.831001E-02;
  Tedata[47]=3.801201E+04; deltadata[47]=3.366693E-02;
  Tedata[48]=4.453328E+04; deltadata[48]=4.001813E-02;
  Tedata[49]=4.979196E+04; deltadata[49]=4.840760E-02;
  Tedata[50]=5.478713E+04; deltadata[50]=5.914255E-02;
                          
  delta=deltadata[0];     
                          
  for (cnt=0; cnt<50; cnt++){
    if (Te>=Tedata[cnt] && Te<=Tedata[cnt+1]) {
      fact=(Te-Tedata[cnt])/(Tedata[cnt+1]-Tedata[cnt]);
      delta=deltadata[cnt]+fact*(deltadata[cnt+1]-deltadata[cnt]);
    }
  }
  if (Te>=Tedata[50]) delta=deltadata[50];
  
  return(delta);
}


double _mueNnfromTe(double Te){
  long cnt;
  double mueNn,fact;
  double Tedata[51],mueNndata[51];
  
  Tedata[0]=2.868000E+02;  mueNndata[0]= 4.836229E+25;
  Tedata[1]=2.957518E+02;  mueNndata[1]= 4.575291E+25;
  Tedata[2]=3.067522E+02;  mueNndata[2]= 4.346909E+25;
  Tedata[3]=3.207755E+02;  mueNndata[3]= 4.119276E+25;
  Tedata[4]=3.419306E+02;  mueNndata[4]= 3.766064E+25;
  Tedata[5]=3.687790E+02;  mueNndata[5]= 3.424062E+25;
  Tedata[6]=4.042948E+02;  mueNndata[6]= 3.049496E+25;
  Tedata[7]=4.509663E+02;  mueNndata[7]= 2.686398E+25;
  Tedata[8]=5.145096E+02;  mueNndata[8]= 2.304831E+25;
  Tedata[9]=5.907615E+02;  mueNndata[9]= 1.986860E+25;
  Tedata[10]=6.822638E+02; mueNndata[10]=1.721883E+25;
  Tedata[11]=7.920666E+02; mueNndata[11]=1.501069E+25;
  Tedata[12]=9.238299E+02; mueNndata[12]=1.317058E+25;
  Tedata[13]=1.083952E+03; mueNndata[13]=1.162449E+25;
  Tedata[14]=1.280273E+03; mueNndata[14]=1.031416E+25;
  Tedata[15]=1.515857E+03; mueNndata[15]=9.222218E+24;
  Tedata[16]=1.745248E+03; mueNndata[16]=8.312265E+24;
  Tedata[17]=1.971427E+03; mueNndata[17]=7.553971E+24;
  Tedata[18]=2.307279E+03; mueNndata[18]=6.699761E+24;
  Tedata[19]=2.740855E+03; mueNndata[19]=5.900078E+24;
  Tedata[20]=3.110065E+03; mueNndata[20]=5.321443E+24;
  Tedata[21]=3.495459E+03; mueNndata[21]=4.867161E+24;
  Tedata[22]=3.957930E+03; mueNndata[22]=4.488593E+24;
  Tedata[23]=4.512896E+03; mueNndata[23]=4.173119E+24;
  Tedata[24]=5.178855E+03; mueNndata[24]=3.910224E+24;
  Tedata[25]=5.978007E+03; mueNndata[25]=3.691145E+24;
  Tedata[26]=6.582945E+03; mueNndata[26]=3.412290E+24;
  Tedata[27]=7.154984E+03; mueNndata[27]=3.145033E+24;
  Tedata[28]=7.841430E+03; mueNndata[28]=2.922319E+24;
  Tedata[29]=8.605688E+03; mueNndata[29]=2.671295E+24;
  Tedata[30]=9.501849E+03; mueNndata[30]=2.454912E+24;
  Tedata[31]=1.035624E+04; mueNndata[31]=2.254927E+24;
  Tedata[32]=1.131832E+04; mueNndata[32]=2.091204E+24;
  Tedata[33]=1.192478E+04; mueNndata[33]=1.975962E+24;
  Tedata[34]=1.265252E+04; mueNndata[34]=1.879927E+24;
  Tedata[35]=1.352582E+04; mueNndata[35]=1.799897E+24;
  Tedata[36]=1.457377E+04; mueNndata[36]=1.733206E+24;
  Tedata[37]=1.583131E+04; mueNndata[37]=1.677630E+24;
  Tedata[38]=1.727350E+04; mueNndata[38]=1.627770E+24;
  Tedata[39]=1.823024E+04; mueNndata[39]=1.552017E+24;
  Tedata[40]=1.937833E+04; mueNndata[40]=1.488889E+24;
  Tedata[41]=2.049241E+04; mueNndata[41]=1.417702E+24;
  Tedata[42]=2.110358E+04; mueNndata[42]=1.315751E+24;
  Tedata[43]=2.240148E+04; mueNndata[43]=1.234090E+24;
  Tedata[44]=2.474838E+04; mueNndata[44]=1.169880E+24;
  Tedata[45]=2.804897E+04; mueNndata[45]=1.129115E+24;
  Tedata[46]=3.257763E+04; mueNndata[46]=1.107596E+24;
  Tedata[47]=3.801201E+04; mueNndata[47]=1.089663E+24;
  Tedata[48]=4.453328E+04; mueNndata[48]=1.074719E+24;
  Tedata[49]=4.979196E+04; mueNndata[49]=1.045097E+24;
  Tedata[50]=5.478713E+04; mueNndata[50]=1.013081E+24;
                          
  mueNn=mueNndata[0];     
                          
  for (cnt=0; cnt<50; cnt++){
    if (Te>=Tedata[cnt] && Te<=Tedata[cnt+1]) {
      fact=(Te-Tedata[cnt])/(Tedata[cnt+1]-Tedata[cnt]);
      mueNn=mueNndata[cnt]+fact*(mueNndata[cnt+1]-mueNndata[cnt]);
    }
  }
  if (Te>=Tedata[50]) mueNn=mueNndata[50]*sqrt(Tedata[50]/Te);
  if (Te<=Tedata[0]) mueNn=mueNndata[0]*sqrt(Tedata[0]/Te);
  
  return(mueNn);
}



int main(){
  double EoverN;
  for (EoverN=1.0e-23; EoverN<3.0e-19; EoverN*=1.2){
  printf("%E %E %E %E %E\n",EoverN,_Teprime(EoverN),_mueNnprime(EoverN),_delta(EoverN),_fq(EoverN));
  }
  return(0);
}
