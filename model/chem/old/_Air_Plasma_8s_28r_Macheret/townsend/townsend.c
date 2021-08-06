#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sqr(a)   ((a)*(a))
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))
#define TRUE 0
#define TRUE 1



double _vdrN2(double EeffoverN){
  double vdrN2data[21],EeffoverNdata[21];
  double vdrN2,fact;
  long cnt;
  
  EeffoverN=max(0.00301E-20,EeffoverN); 
  EeffoverN=min(99.99E-20,EeffoverN); 
  vdrN2=0.0; /* to avoid compiler warning */
  EeffoverNdata[0]=0.003E-20;  vdrN2data[0]=1.1E3;   
  EeffoverNdata[1]=0.005E-20;  vdrN2data[1]=1.6E3;   
  EeffoverNdata[2]=0.007E-20;  vdrN2data[2]=2.0E3;   
  EeffoverNdata[3]=0.01E-20;   vdrN2data[3]=2.4E3;   
  EeffoverNdata[4]=0.03E-20;   vdrN2data[4]=3.1E3;   
  EeffoverNdata[5]=0.05E-20;   vdrN2data[5]=3.5E3;   
  EeffoverNdata[6]=0.07E-20;   vdrN2data[6]=3.9E3;   
  EeffoverNdata[7]=0.10E-20;   vdrN2data[7]=4.4E3;   
  EeffoverNdata[8]=0.30E-20;   vdrN2data[8]=7.5E3;   
  EeffoverNdata[9]=0.50E-20;   vdrN2data[9]=11.0E3;  
  EeffoverNdata[10]=0.70E-20;  vdrN2data[10]=14.0E3; 
  EeffoverNdata[11]=1.0E-20;   vdrN2data[11]=18.0E3; 
  EeffoverNdata[12]=3.0E-20;   vdrN2data[12]=42.0E3; 
  EeffoverNdata[13]=5.0E-20;   vdrN2data[13]=61.0E3; 
  EeffoverNdata[14]=7.0E-20;   vdrN2data[14]=79.0E3; 
  EeffoverNdata[15]=10.0E-20;  vdrN2data[15]=100.0E3; 
  EeffoverNdata[16]=20.0E-20;  vdrN2data[16]=200.0E3; 
  EeffoverNdata[17]=30.0E-20;  vdrN2data[17]=290.0E3; 
  EeffoverNdata[18]=50.0E-20;  vdrN2data[18]=420.0E3; 
  EeffoverNdata[19]=80.0E-20;  vdrN2data[19]=530.0E3; 
  EeffoverNdata[20]=100.0E-20; vdrN2data[20]=590.0E3; 

  for (cnt=0; cnt<20; cnt++){
    if (EeffoverN>=EeffoverNdata[cnt] && EeffoverN<=EeffoverNdata[cnt+1]) {
      fact=(EeffoverN-EeffoverNdata[cnt])/(EeffoverNdata[cnt+1]-EeffoverNdata[cnt]);
      vdrN2=vdrN2data[cnt]+fact*(vdrN2data[cnt+1]-vdrN2data[cnt]);
    }
  }
  return(vdrN2);
}

double _vdrO2(double EeffoverN){
  double vdrO2data[21],EeffoverNdata[21];
  double vdrO2,fact;
  long cnt;
  
  EeffoverN=max(0.00301E-20,EeffoverN); 
  EeffoverN=min(99.99E-20,EeffoverN); 
  vdrO2=0.0; /* to avoid compiler warning */
  EeffoverNdata[0]=0.003E-20;  vdrO2data[0]=2.6E3;   
  EeffoverNdata[1]=0.005E-20;  vdrO2data[1]=3.7E3;   
  EeffoverNdata[2]=0.007E-20;  vdrO2data[2]=4.1E3;   
  EeffoverNdata[3]=0.01E-20;   vdrO2data[3]=4.3E3;   
  EeffoverNdata[4]=0.03E-20;   vdrO2data[4]=5.4E3;   
  EeffoverNdata[5]=0.05E-20;   vdrO2data[5]=7.3E3;   
  EeffoverNdata[6]=0.07E-20;   vdrO2data[6]=9.2E3;   
  EeffoverNdata[7]=0.10E-20;   vdrO2data[7]=10.0E3;  
  EeffoverNdata[8]=0.30E-20;   vdrO2data[8]=22.0E3;  
  EeffoverNdata[9]=0.50E-20;   vdrO2data[9]=26.0E3;  
  EeffoverNdata[10]=0.70E-20;  vdrO2data[10]=28.0E3; 
  EeffoverNdata[11]=1.0E-20;   vdrO2data[11]=31.0E3; 
  EeffoverNdata[12]=3.0E-20;   vdrO2data[12]=72.0E3; 
  EeffoverNdata[13]=5.0E-20;   vdrO2data[13]=110.0E3; 
  EeffoverNdata[14]=7.0E-20;   vdrO2data[14]=120.0E3; 
  EeffoverNdata[15]=10.0E-20;   vdrO2data[15]=160.0E3; 
  EeffoverNdata[16]=20.0E-20;   vdrO2data[16]=260.0E3; 
  EeffoverNdata[17]=30.0E-20;   vdrO2data[17]=330.0E3; 
  EeffoverNdata[18]=50.0E-20;   vdrO2data[18]=450.0E3; 
  EeffoverNdata[19]=80.0E-20;   vdrO2data[19]=640.0E3; 
  EeffoverNdata[20]=100.0E-20;  vdrO2data[20]=710.0E3; 

  for (cnt=0; cnt<20; cnt++){
    if (EeffoverN>=EeffoverNdata[cnt] && EeffoverN<=EeffoverNdata[cnt+1]) {
      fact=(EeffoverN-EeffoverNdata[cnt])/(EeffoverNdata[cnt+1]-EeffoverNdata[cnt]);
      vdrO2=vdrO2data[cnt]+fact*(vdrO2data[cnt+1]-vdrO2data[cnt]);
    }
  }
  return(vdrO2);
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


void FindkN2vsEoverN(){
  double EoverN,vartheta,k,knew,k1,k2,k3;
  long cnt;
  
  cnt=0;
  for (EoverN=1e-20; EoverN<=350.0e-20; EoverN*=1.1){
    
    vartheta=EoverN/1.0e-20;
    if ((vartheta>=3.0 && vartheta <= 30.0) || TRUE) {
      k1=pow(10.0,-8.3-36.5/vartheta);
    } else {
      k1=0.0;
    }
    if ((vartheta>=13.6 && vartheta <= 54.5) || TRUE) {
      k2=_vdrN2(EoverN)*1E6*pow(10.0,-24.4406)*sqr(3.2271E20*EoverN-32.2);
    } else {
      k2=0.0;
    }
    if ((vartheta>=31.0 && vartheta <= 247.9) || TRUE) {
      k3=_vdrN2(EoverN)*1E6*4.6483E-20*exp(-1.1311E-18/EoverN);
    } else {
      k3=0.0;
    }
    if (vartheta<=30.0) k=k1; else k=k3;
    printf("EoverNdata[%ld]=%E;  kdata[%ld]=%E; \n",cnt,EoverN,cnt,k);
    cnt++;
  }

}


void FindkN2vsTe(){
  double EoverN,vartheta,k,knew,k1,k2,k3;
  long cnt;
  
  cnt=0;
  for (EoverN=1e-20; EoverN<=30.0e-20; EoverN*=1.1){
    
    vartheta=EoverN/1.0e-20;
    if ((vartheta>=3.0 && vartheta <= 30.0) || TRUE) {
      k1=pow(10.0,-8.3-36.5/vartheta);
    } else {
      k1=0.0;
    }
    if ((vartheta>=13.6 && vartheta <= 54.5) || TRUE) {
      k2=_vdrN2(EoverN)*1E6*pow(10.0,-24.4406)*sqr(3.2271E20*EoverN-32.2);
    } else {
      k2=0.0;
    }
    if ((vartheta>=31.0 && vartheta <= 247.9) || TRUE) {
      k3=_vdrN2(EoverN)*1E6*4.6483E-20*exp(-1.1311E-18/EoverN);
    } else {
      k3=0.0;
    }
    if (vartheta<=30.0) k=k1; else k=k3;
    printf("Tedata[%ld]=%E;  kdata[%ld]=%E; \n",cnt,_Teprime(EoverN),cnt,k);
    cnt++;
  }

}


void FindkO2vsEoverN(){
  double EoverN,vartheta,k,knew,k1,k2,k3;
  long cnt;

  cnt=0;
  for (EoverN=1e-20; EoverN<=350.0e-20; EoverN*=1.1){
    
    vartheta=EoverN/1.0e-20;
    if ((vartheta>=3.0 && vartheta <= 30.0) || TRUE) {
      k1=pow(10.0,-8.8-28.1/vartheta);
    } else {
      k1=0.0;
    }
    if ((vartheta>=13.6 && vartheta <= 54.5) || TRUE) {
      k2=_vdrO2(EoverN)*1E6*pow(10.0,-24.4406)*sqr(3.2271E20*EoverN-32.2);
    } else {
      k2=0.0;
    }
    if ((vartheta>=31.0 && vartheta <= 247.9) || TRUE) {
      k3=_vdrO2(EoverN)*1E6*4.6483E-20*exp(-1.1311E-18/EoverN);
    } else {
      k3=0.0;
    }
    if (vartheta<=30.0) k=k1; else k=k3;
    printf("EoverNdata[%ld]=%E;  kdata[%ld]=%E; \n",cnt,EoverN,cnt,k);
    cnt++;
  }

}


void FindkO2vsTe(){
  double EoverN,vartheta,k,knew,k1,k2,k3;
  long cnt;

  cnt=0;
  for (EoverN=1e-20; EoverN<=30.0e-20; EoverN*=1.1){
    
    vartheta=EoverN/1.0e-20;
    if ((vartheta>=3.0 && vartheta <= 30.0) || TRUE) {
      k1=pow(10.0,-8.8-28.1/vartheta);
    } else {
      k1=0.0;
    }
    if ((vartheta>=13.6 && vartheta <= 54.5) || TRUE) {
      k2=_vdrO2(EoverN)*1E6*pow(10.0,-24.4406)*sqr(3.2271E20*EoverN-32.2);
    } else {
      k2=0.0;
    }
    if ((vartheta>=31.0 && vartheta <= 247.9) || TRUE) {
      k3=_vdrO2(EoverN)*1E6*4.6483E-20*exp(-1.1311E-18/EoverN);
    } else {
      k3=0.0;
    }
    if (vartheta<=30.0) k=k1; else k=k3;
    printf("Tedata[%ld]=%E;  kdata[%ld]=%E; \n",cnt,_Teprime(EoverN),cnt,k);
    cnt++;
  }

}


int main(){

  FindkO2vsTe();
  return(0);
}

