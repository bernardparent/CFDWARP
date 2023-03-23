// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2022 Prasanna Thoguluva Rajendran

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

/*
The molar fractions have been obtained using curve fits proposed by A. D'Angola et al.[1][2].
[1]  AD' Angola, G Colonna, C Gorse, and M. Capitelli. "Thermodynamic and transport properties 
     in equilibrium air plasmas in a wide pressure and temperature range", The European Physical
     Journal D 46:129-150, 2008.
[2]  AD' Angola, G Colonna, C Gorse, and M. Capitelli. "Thermodynamic properties of high temperature 
     air in local thermodynamic equilibrium: II accurate analytical expression or electron molar
     fractions", The European Physical Journal D 65:453-457, 2011.
*/

#include "share.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>


int chkarg (int argc, char **argv, char *arg) {
  int cnt, tmp;

  tmp = 0;
  for (cnt = 1; cnt < argc; cnt++) {
    if (strcmp (argv[cnt], arg) == 0) {
      tmp = cnt;
    }
  }
  return (tmp);
}

double _chi_N2(double T, double P){
  double sum,logP,a1,a2,c1,c2,Delta1,Delta2,sigma1,sigma2;
  
  logP=log(P);
  a2=exp(-4.037422e-1-7.147343e-4*logP+4.492235e-4*pow(logP,2.0)+9.648313e-5*pow(logP,3.0)-1.284083e-8*pow(logP,4.0));
  a1=0.8-a2;
  c1=exp(8.110148+4.553237e-2*logP-8.193725e-4*pow(logP,2.0)-2.156896e-4*pow(logP,3.0));
  c2=exp(8.812799+5.665883e-2*logP+1.293767e-3*pow(logP,2.0));
  Delta1=exp(6.561427+1.422222e-1*logP-7.476314e-4*pow(logP,2.0)-8.715905e-4*pow(logP,3.0));
  Delta2=exp(7.016774+1.058804e-1*logP+3.292541e-3*pow(logP,2.0)+2.267238e-4*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sum=a1*sigma1+a2*sigma2;
  return(0.8-sum);
}

double _chi_N2plus(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3;
  
  logP=log(P);
  a1=exp(-9.746298+4.199007e-1*logP+3.143417e-3*pow(logP,2.0)+3.882378e-4*pow(logP,3.0));
  a3=-exp(-9.992503+4.689265e-1*logP+1.182887e-3*pow(logP,2.0)-1.176687e-4*pow(logP,3.0));
  a2=-a3-a1;
  c1=exp(8.884490+5.573065e-2*logP+1.616982e-3*pow(logP,2.0)+6.738352e-5*pow(logP,3.0));
  c2=exp(9.203463+7.494796e-2*logP+2.541069e-3*pow(logP,2.0)+7.257196e-5*pow(logP,3.0)+6.051419e-6*pow(logP,4.0));
  c3=exp(9.449201+6.238298e-2*logP+1.564727e-3*pow(logP,2.0)+5.575073e-5*pow(logP,3.0));
  Delta1=exp(6.552069+1.058201e-1*logP+3.989777e-3*pow(logP,2.0)+1.801416e-4*pow(logP,3.0));
  Delta2=exp(7.294752-1.099569e-3*logP+4.040325e-3*pow(logP,2.0)+2.717526e-3*pow(logP,3.0)-5.081078e-5*pow(logP,4.0)-3.474609e-5*pow(logP,5.0));
  Delta3=exp(7.762006+1.260807e-1*logP+2.223845e-3*pow(logP,2.0)-1.231135e-4*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_N(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=8.188731e-1+2.581889e-3*logP+1.395103e-4*pow(logP,2.0);
  a3=-exp(-3.552214+4.085111e-1*logP-2.961084e-2*pow(logP,2.0));
  a2=-a3-a1;
  c1=exp(8.812279+5.474146e-2*logP+1.019131e-3*pow(logP,2.0));
  cm=exp(8.405373+4.371184e-2*logP+1.893389e-3*pow(logP,2.0)+1.927737e-4*pow(logP,3.0));
  c2=exp(9.516473+6.520807e-2*logP+1.270979e-3*pow(logP,2.0)-3.857140e-5*pow(logP,3.0)-5.540006e-6*pow(logP,4.0));
  c3=exp(9.864059+3.659617e-2*logP+7.907898e-3*pow(logP,2.0));
  Delta1=exp(7.051737+1.128378e-1*logP+2.407727e-3*pow(logP,2.0)-1.247502e-5*pow(logP,3.0));
  Deltam=exp(6.923056+1.987863e-1*logP+3.645361e-3*pow(logP,2.0)-5.777817e-4*pow(logP,3.0));
  Delta2=exp(7.949412+1.206502e-1*logP+1.785666e-3*pow(logP,2.0)-3.344976e-5*pow(logP,3.0));
  Delta3=exp(8.814892+5.421480e-2*logP+1.537056e-3*pow(logP,2.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_Nplus(double T, double P){
  double sum,logP,a1,a2,a3,a4,c1,c2,c3,c4,Delta1,Delta2,Delta3,Delta4,sigma1,sigma2,sigma3,sigma4,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=exp(-1.211184+2.634222e-4*logP+2.560470e-3*pow(logP,2.0));
  a2=exp(-2.230927+2.047906e-2*logP-2.220684e-3*pow(logP,2.0));
  a4=-exp(-1.200529-3.074481e-2*logP+4.780888e-3*pow(logP,2.0)+8.341989e-4*pow(logP,3.0)+6.160353e-6*pow(logP,4.0)-2.708386e-6*pow(logP,5.0));
  a3=-a4-a2-a1;
  c1=exp(9.494309+5.588021e-2*logP+2.479295e-3*pow(logP,2.0)+5.228102e-4*pow(logP,3.0)+5.047984e-5*pow(logP,4.0)-1.606423e-6*pow(logP,5.0)-8.671283e-7*pow(logP,6.0)-5.919943e-8*pow(logP,7.0));
  cm=exp(9.276675+8.451791e-2*logP-7.509912e-3*pow(logP,2.0)+1.762683e-3*pow(logP,3.0)-2.856325e-4*pow(logP,4.0)+3.392616e-5*pow(logP,5.0)-5.010435e-6*pow(logP,6.0)+3.875277e-7*pow(logP,7.0));
  c2=exp(9.511003+8.651691e-2*logP-5.130145e-5*pow(logP,2.0)-2.847046e-4*pow(logP,3.0));
  c3=exp(1.037880e+1+6.497169e-2*logP+3.027193e-3*pow(logP,2.0)+1.559114e-4*pow(logP,3.0)-2.230902e-7*pow(logP,4.0)+3.440378e-6*pow(logP,5.0));
  c4=exp(1.025494e+1+6.494578e-2*logP+1.277401e-3*pow(logP,2.0));
  Delta1=exp(8.228344+2.288911e-1*logP-7.989931e-4*pow(logP,2.0)-1.145501e-3*pow(logP,3.0));
  Deltam=exp(7.931270-4.388112e-2*logP+2.643605e-2*pow(logP,2.0)-1.501361e-3*pow(logP,3.0)-2.178943e-4*pow(logP,4.0)+2.476492e-5*pow(logP,5.0));
  Delta2=exp(7.645166+8.574186e-2*logP+2.708947e-4*pow(logP,2.0)+6.210369e-4*pow(logP,3.0));
  Delta3=exp(8.810742+1.305064e-1*logP-1.083168e-3*pow(logP,2.0)+4.025862e-5*pow(logP,3.0)+1.348428e-4*pow(logP,4.0)-2.273123e-5*pow(logP,5.0));
  Delta4=exp(8.187912+1.182600e-1*logP+6.307194e-3*pow(logP,2.0)+2.948945e-4*pow(logP,3.0)+1.136590e-6*pow(logP,4.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sigma4=1.0/(1.0+exp(-2.0*(T-c4)/Delta4));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3+a4*sigma4;
  return(fmax(0.0,sum));
}

double _chi_Nplusplus(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=exp(-1.320561+4.613513e-3*logP+1.563146e-3*pow(logP,2.0)+9.805924e-5*pow(logP,3.0));
  a3=-exp(-2.441955+1.600937e-2*logP-1.796504e-2*pow(logP,2.0)+4.445771e-5*pow(logP,3.0));
  a2=-a3-a1;
  c1=exp(1.018105e+1+6.182886e-2*logP+4.542717e-4*pow(logP,2.0)+1.665348e-4*pow(logP,3.0)-1.688929e-5*pow(logP,4.0));
  cm=exp(1.020635e+1+6.787015e-2*logP+2.930559e-3*pow(logP,2.0)-2.387278e-5*pow(logP,3.0)-1.580874e-5*pow(logP,4.0));
  c2=exp(1.071778e+1+6.267958e-2*logP+1.384143e-3*pow(logP,2.0)+1.803319e-5*pow(logP,3.0));
  c3=exp(1.081164e+1+6.929471e-2*logP+3.005312e-3*pow(logP,2.0)+5.422861e-5*pow(logP,3.0));
  Delta1=exp(8.328213+7.134338e-2*logP+8.440573e-3*pow(logP,2.0)-1.913632e-4*pow(logP,3.0));
  Deltam=exp(8.637046+1.730873e-1*logP-2.312739e-3*pow(logP,2.0)-1.253255e-4*pow(logP,3.0)+6.870714e-5*pow(logP,4.0));
  Delta2=exp(8.558877+1.280075e-1*logP+7.408166e-3*pow(logP,2.0)-6.068102e-5*pow(logP,3.0)-3.499092e-5*pow(logP,4.0));
  Delta3=exp(9.008121+1.059058e-1*logP+3.835047e-3*pow(logP,2.0)-5.778232e-4*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_Nplusplusplus(double T, double P){
  double sum,logP,a1,c1,c2,cm,Delta1,Delta2,Deltam,sigma1,sigma2,sigmam;
  
  logP=log(P);
  a1=exp(-1.339800+1.954622e-2*logP-3.939015e-3*pow(logP,2.0)-4.170049e-4*pow(logP,3.0));
  c1=exp(1.070665e+1+6.722548e-2*logP+6.769799e-5*pow(logP,2.0)+4.111595e-5*pow(logP,3.0));
  cm=exp(1.066404e+1+5.711793e-2*logP+1.063676e-3*pow(logP,2.0)-1.137507e-6*pow(logP,3.0));
  c2=exp(1.105085e+1+5.890335e-2*logP+1.918852e-3*pow(logP,2.0)+9.521033e-5*pow(logP,3.0));
  Delta1=exp(9.340050+5.929963e-2*logP+1.505109e-3*pow(logP,2.0)+2.034159e-4*pow(logP,3.0));
  Deltam=exp(8.726521+1.521811e-1*logP+2.430293e-3*pow(logP,2.0)-4.716643e-4*pow(logP,3.0));
  Delta2=exp(9.258763+1.273121e-1*logP-6.021997e-4*pow(logP,2.0)-2.540618e-4*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sum=a1*sigmam*sigma1-a1*sigma2;
  return(fmax(0.0,sum));
}

double _chi_Nplusplusplusplus(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=exp(-1.849635-4.491118e-3*logP-3.702617e-4*pow(logP,2.0));
  a3=-exp(-6.074622e-1+6.073274e-1*logP+9.963043e-2*pow(logP,2.0)+5.415504e-3*pow(logP,3.0));
  a2=-a3-a1;
  c1=exp(1.100960e+1+7.368509e-2*logP+1.075589e-3*pow(logP,2.0));
  cm=exp(1.100986e+1+4.882927e-2*logP+3.853047e-4*pow(logP,2.0)-1.475967e-6*pow(logP,3.0));
  c2=exp(1.206372e+1-1.734608e-3*logP-1.447988e-2*pow(logP,2.0)+1.590266e-3*pow(logP,3.0));
  c3=exp(1.280436e+1-1.896326e-1*logP+2.801196e-2*pow(logP,2.0));
  Delta1=exp(9.329126+7.704932e-2*logP+2.666225e-3*pow(logP,2.0));
  Deltam=exp(9.006971+1.074664e-1*logP-1.472426e-3*pow(logP,2.0)-2.722012e-4*pow(logP,3.0));
  Delta2=exp(1.019997e+1-1.423777e-1*logP-4.095877e-2*pow(logP,2.0)+2.180861e-3*pow(logP,3.0)+2.368183e-4*pow(logP,4.0));
  Delta3=exp(1.103058e+1-2.553162e-1*logP+2.330651e-2*pow(logP,2.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_O2(double T, double P){
  double sum,logP,a1,a2,c1,c2,Delta1,Delta2,sigma1,sigma2;
  
  logP=log(P);
  a2=exp(-1.685730-3.728595e-2*logP-5.172240e-3*pow(logP,2.0)+2.021941e-4*pow(logP,3.0)+6.195083e-5*pow(logP,4.0)-5.999263e-6*pow(logP,5.0));
  a1=0.2-a2;
  c1=exp(7.851965-4.971670e-2*logP-1.438515e-2*pow(logP,2.0)-8.639710e-4*pow(logP,3.0));
  c2=exp(8.148167+4.575379e-2*logP+1.841872e-4*pow(logP,2.0));
  Delta1=exp(6.500431+7.318423e-2*logP-2.704126e-3*pow(logP,2.0)-2.824658e-4*pow(logP,3.0));
  Delta2=exp(6.459154+1.486515e-1*logP+5.919587e-3*pow(logP,2.0)-3.159509e-5*pow(logP,3.0)-4.048213e-5*pow(logP,4.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sum=a1*sigma1+a2*sigma2;
  return(0.2-sum);
}

double _chi_O2plus(double T, double P){
  double sum,logP,a1,a2,a3,a4,c1,c2,c3,c4,Delta1,Delta2,Delta3,Delta4,sigma1,sigma2,sigma3,sigma4;
  
  logP=log(P);
  a1=exp(-1.373444e+1+6.627381e-1*logP-1.950471e-2*pow(logP,2.0)+7.469315e-4*pow(logP,3.0)+1.358278e-4*pow(logP,4.0));
  a2=exp(-1.419853e+1+4.889623e-1*logP-6.123742e-3*pow(logP,2.0)+5.940913e-4*pow(logP,3.0)+9.783232e-5*pow(logP,4.0));
  a4=-exp(-1.342851e+1+6.025406e-1*logP-1.482459e-2*pow(logP,2.0)+1.461126e-4*pow(logP,3.0)+1.408990e-4*pow(logP,4.0));
  a3=-a4-a2-a1;
  c1=exp(8.794853+4.659480e-2*logP+5.610403e-4*pow(logP,2.0)+1.044006e-4*pow(logP,3.0)-1.835079e-5*pow(logP,4.0));
  c2=exp(8.991604+5.142449e-2*logP+1.298498e-3*pow(logP,2.0)+4.051458e-4*pow(logP,3.0)+1.170299e-5*pow(logP,4.0));
  c3=exp(9.563817+7.340431e-2*logP+7.915772e-4*pow(logP,2.0)-1.592330e-4*pow(logP,3.0)-1.027704e-5*pow(logP,4.0));
  c4=exp(8.900254+3.563862e-2*logP+1.399785e-3*pow(logP,2.0)+1.003372e-4*pow(logP,3.0)-3.618984e-5*pow(logP,4.0));
  Delta1=exp(7.268996+9.440745e-2*logP-2.146537e-3*pow(logP,2.0)+4.167152e-5*pow(logP,3.0)+3.077941e-5*pow(logP,4.0));
  Delta2=exp(7.456563+1.277214e-1*logP+8.479742e-3*pow(logP,2.0)+8.341173e-4*pow(logP,3.0)-1.597360e-4*pow(logP,4.0));
  Delta3=exp(7.834428+1.245447e-1*logP+4.949361e-3*pow(logP,2.0)+3.875066e-5*pow(logP,3.0)-2.966365e-5*pow(logP,4.0));
  Delta4=exp(7.450971+9.288765e-2*logP-1.491663e-3*pow(logP,2.0)+7.510663e-4*pow(logP,3.0)-9.458429e-5*pow(logP,4.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sigma4=1.0/(1.0+exp(-2.0*(T-c4)/Delta4));
  sum=a1*sigma1+a2*sigma2+a3*sigma3+a4*sigma4;
  return(fmax(0.0,sum));
}

double _chi_O2minus(double T, double P){
  double sum,logP,a1,a2,a3,a4,c1,c2,c3,c4,Delta1,Delta2,Delta3,Delta4,sigma1,sigma2,sigma3,sigma4;
  
  logP=log(P);
  a1=exp(-2.009128e+1+1.218472*logP-1.023713e-2*pow(logP,2.0)-3.693487e-4*pow(logP,3.0));
  a3=-exp(-2.169571e+1+1.231117*logP+1.792651e-3*pow(logP,2.0)+2.558252e-4*pow(logP,3.0)-1.732401e-4*pow(logP,4.0)+8.498995e-6*pow(logP,5.0)+1.264359e-6*pow(logP,6.0));
  a4=-exp(-2.472870e+1+1.526884*logP+1.203852e-3*pow(logP,2.0)+1.430794e-4*pow(logP,3.0));
  a2=-a4-a3-a1;
  c1=exp(8.151744+5.269669e-2*logP+1.328087e-3*pow(logP,2.0)+9.918314e-5*pow(logP,3.0)+6.931618e-6*pow(logP,4.0));
  c2=exp(8.327753+6.884887e-2*logP+2.843931e-3*pow(logP,2.0)+1.083879e-4*pow(logP,3.0));
  c3=exp(8.601320+7.342289e-2*logP-9.411900e-4*pow(logP,2.0)-1.339663e-4*pow(logP,3.0)+3.126379e-5*pow(logP,4.0));
  c4=exp(9.428115+5.014640e-2*logP-3.340382e-4*pow(logP,2.0)+7.998702e-6*pow(logP,3.0));
  Delta1=exp(6.140093+1.051897e-1*logP+2.939827e-3*pow(logP,2.0)+1.422812e-4*pow(logP,3.0));
  Delta2=exp(6.644117+1.374513e-1*logP+4.095263e-3*pow(logP,2.0)+8.402722e-5*pow(logP,3.0)-1.242256e-5*pow(logP,4.0)-7.990825e-6*pow(logP,5.0)+8.075101e-7*pow(logP,6.0)+2.001120e-7*pow(logP,7.0));
  Delta3=exp(6.985630+8.256947e-2*logP+9.999196e-3*pow(logP,2.0)-2.953396e-5*pow(logP,3.0)-1.526330e-4*pow(logP,4.0));
  Delta4=exp(7.530896+1.558330e-1*logP+4.905502e-3*pow(logP,2.0)-8.411242e-4*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sigma4=1.0/(1.0+exp(-2.0*(T-c4)/Delta4));
  sum=a1*sigma1+a2*sigma2+a3*sigma3+a4*sigma4;
  return(fmax(0.0,sum));
}

double _chi_O(double T, double P){
  double sum,logP,a1,a2,a3,a4,c1,c2,c3,c4,Delta1,Delta2,Delta3,Delta4,sigma1,sigma2,sigma3,sigma4,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=exp(-1.139281-1.050647e-2*logP-1.022007e-3*pow(logP,2.0)-4.830320e-5*pow(logP,3.0)-3.531305e-6*pow(logP,4.0)-2.296630e-7*pow(logP,5.0));
  a2=-exp(-4.979500+2.665257e-1*logP+1.458327e-2*pow(logP,2.0)-2.533456e-3*pow(logP,3.0)-3.704428e-4*pow(logP,4.0)+2.339924e-5*pow(logP,5.0));
  a4=-exp(-1.615594-5.157778e-3*logP-1.550658e-3*pow(logP,2.0)-1.264223e-4*pow(logP,3.0)+2.343404e-5*pow(logP,4.0)+3.184705e-6*pow(logP,5.0));
  a3=-a4-a2-a1;
  c1=exp(8.145639+5.431612e-2*logP+2.023998e-3*pow(logP,2.0)+1.003745e-4*pow(logP,3.0));
  cm=exp(7.940783+6.741169e-2*logP-2.087042e-3*pow(logP,2.0)+3.972481e-4*pow(logP,3.0)-3.481686e-5*pow(logP,4.0)+1.485858e-6*pow(logP,5.0));
  c2=exp(9.811744+7.436247e-2*logP-1.239267e-4*pow(logP,2.0)+6.132060e-4*pow(logP,3.0));
  c3=exp(8.819734+5.805213e-2*logP+1.501067e-3*pow(logP,2.0)+2.511693e-5*pow(logP,3.0));
  c4=exp(9.554277+6.746571e-2*logP+8.910292e-4*pow(logP,2.0)-4.496226e-5*pow(logP,3.0));
  Delta1=exp(6.576786+1.491330e-1*logP+3.724868e-3*pow(logP,2.0)-1.382563e-4*pow(logP,3.0)+1.947915e-6*pow(logP,4.0)+1.082756e-6*pow(logP,5.0));
  Deltam=exp(6.664764+4.575484e-2*logP+4.557480e-3*pow(logP,2.0));
  Delta2=exp(9.044853+5.997097e-3*logP+4.532508e-4*pow(logP,2.0)+6.756744e-4*pow(logP,3.0));
  Delta3=exp(6.918048+9.326905e-2*logP+2.506390e-3*pow(logP,2.0)+1.395474e-4*pow(logP,3.0));
  Delta4=exp(8.033301+1.233674e-1*logP+1.651217e-3*pow(logP,2.0)-3.811131e-5*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sigma4=1.0/(1.0+exp(-2.0*(T-c4)/Delta4));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3+a4*sigma4;
  return(fmax(0.0,sum));
}

double _chi_Ominus(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3;
  
  logP=log(P);
  a1=exp(-1.492297e+1+9.064321e-1*logP-8.724265e-3*pow(logP,2.0)-2.165125e-4*pow(logP,3.0)+1.166368e-4*pow(logP,4.0));
  a2=exp(-1.175041e+1+7.618857e-1*logP+1.501595e-3*pow(logP,2.0)+1.781504e-4*pow(logP,3.0)+9.215991e-6*pow(logP,4.0));
  a3=-a1-a2;
  c1=exp(8.415326+5.157258e-2*logP+2.024706e-3*pow(logP,2.0)+1.312425e-4*pow(logP,3.0)+1.315036e-5*pow(logP,4.0));
  c2=exp(9.270258+5.316281e-2*logP+1.482070e-3*pow(logP,2.0)-5.476676e-5*pow(logP,3.0)-9.733849e-6*pow(logP,4.0));
  c3=exp(9.598507+6.569448e-2*logP+5.303147e-4*pow(logP,2.0)-9.613381e-5*pow(logP,3.0)-7.330576e-6*pow(logP,4.0));
  Delta1=exp(6.462668+6.272626e-2*logP-6.193918e-3*pow(logP,2.0)+6.376014e-4*pow(logP,3.0)+2.245471e-4*pow(logP,4.0));
  Delta2=exp(7.724023+9.838486e-2*logP+4.215920e-3*pow(logP,2.0)-6.990084e-5*pow(logP,3.0)-3.230965e-5*pow(logP,4.0));
  Delta3=exp(7.809041+1.423877e-1*logP+4.366188e-3*pow(logP,2.0)-8.184536e-5*pow(logP,3.0)-1.524608e-5*pow(logP,4.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_Oplus(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=exp(-2.319093-7.610174e-3*logP-1.953269e-3*pow(logP,2.0)-3.002482e-4*pow(logP,3.0)-1.751192e-5*pow(logP,4.0));
  a3=-exp(-1.436906+2.872384e-1*logP+2.978317e-2*pow(logP,2.0)+1.769679e-4*pow(logP,3.0)-9.414001e-5*pow(logP,4.0));
  a2=-a3-a1;
  c1=exp(9.588569+6.997026e-2*logP+9.769379e-4*pow(logP,2.0)-6.246775e-5*pow(logP,3.0)-4.877947e-6*pow(logP,4.0));
  cm=exp(9.115221+6.168847e-2*logP+2.270858e-3*pow(logP,2.0)+1.412631e-4*pow(logP,3.0));
  c2=exp(1.020364e+1+6.299762e-2*logP-1.091887e-3*pow(logP,2.0)+3.702998e-5*pow(logP,3.0));
  c3=exp(1.027215e+1+4.672465e-2*logP+1.597850e-4*pow(logP,2.0)+9.311678e-6*pow(logP,3.0));
  Delta1=exp(8.044970+1.175891e-1*logP+1.645336e-3*pow(logP,2.0)-9.489377e-5*pow(logP,3.0)-9.694619e-6*pow(logP,4.0));
  Deltam=exp(7.651684+1.477558e-1*logP-1.967294e-3*pow(logP,2.0)-9.075769e-4*pow(logP,3.0));
  Delta2=exp(8.680331+1.325526e-1*logP+2.754338e-3*pow(logP,2.0)-7.964755e-5*pow(logP,3.0));
  Delta3=exp(8.696369+1.339624e-1*logP+1.995427e-3*pow(logP,2.0)-3.323281e-5*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_Oplusplus(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=7.063013e-2-5.187789e-4*logP-9.288238e-6*pow(logP,2.0);
  a3=-exp(-2.991458-5.757422e-2*logP-3.835760e-3*pow(logP,2.0));
  a2=-a3-a1;
  c1=exp(1.029003e+1+4.517420e-2*logP-1.618224e-5*pow(logP,2.0)+2.245678e-4*pow(logP,3.0)+3.130833e-5*pow(logP,4.0)-2.423868e-6*pow(logP,5.0)-3.903368e-7*pow(logP,6.0));
  cm=exp(1.029386e+1+8.048612e-2*logP-4.497818e-4*pow(logP,2.0)+3.852087e-5*pow(logP,3.0));
  c2=exp(1.082680e+1+7.388982e-2*logP+9.267668e-4*pow(logP,2.0));
  c3=exp(1.078471e+1+5.999115e-2*logP+1.044468e-3*pow(logP,2.0));
  Delta1=exp(8.449025+1.233942e-1*logP-3.128794e-3*pow(logP,2.0)-5.456369e-4*pow(logP,3.0)+5.445584e-5*pow(logP,4.0)+6.520078e-6*pow(logP,5.0));
  Deltam=exp(8.843594+4.195145e-2*logP+1.187095e-2*pow(logP,2.0)+1.964457e-4*pow(logP,3.0)-4.989937e-5*pow(logP,4.0)-9.711143e-7*pow(logP,5.0));
  Delta2=exp(9.267200+5.532633e-2*logP+2.362320e-3*pow(logP,2.0)+6.299569e-4*pow(logP,3.0)+1.122230e-5*pow(logP,4.0)-2.869166e-6*pow(logP,5.0)-4.451869e-7*pow(logP,6.0));
  Delta3=exp(8.785646+9.165132e-2*logP+9.925663e-4*pow(logP,2.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_Oplusplusplus(double T, double P){
  double sum,logP,a1,a2,c1,c2,cm,Delta1,Delta2,Deltam,sigma1,sigma2,sigmam;
  
  logP=log(P);
  a1=exp(-2.760009+3.495500e-3*logP-5.357864e-3*pow(logP,2.0)-2.144466e-4*pow(logP,3.0)+9.251860e-6*pow(logP,4.0)-9.005345e-7*pow(logP,5.0));
  a2=-a1;
  c1=exp(1.074207e+1+5.260080e-2*logP+4.936255e-4*pow(logP,2.0)-4.405321e-5*pow(logP,3.0)-3.025027e-6*pow(logP,4.0)-5.425422e-7*pow(logP,5.0));
  cm=exp(1.077506e+1+6.587529e-2*logP+2.491665e-4*pow(logP,2.0)+1.077355e-4*pow(logP,3.0));
  c2=exp(1.111558e+1+5.973321e-2*logP+2.038965e-3*pow(logP,2.0)+9.054082e-5*pow(logP,3.0));
  Delta1=exp(8.835975+1.411710e-1*logP+2.773994e-3*pow(logP,2.0)-6.211959e-4*pow(logP,3.0)+3.813517e-6*pow(logP,4.0)+1.323357e-5*pow(logP,5.0)-1.119305e-6*pow(logP,6.0)-3.062376e-7*pow(logP,7.0));
  Deltam=exp(9.367809+3.868631e-2*logP-7.976461e-4*pow(logP,2.0)+6.108727e-4*pow(logP,3.0));
  Delta2=exp(9.317898+1.146590e-1*logP-4.219919e-4*pow(logP,2.0)-1.986513e-4*pow(logP,3.0)-9.501572e-6*pow(logP,4.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sum=a1*sigmam*sigma1+a2*sigma2;
  return(fmax(0.0,sum));
}

double _chi_Oplusplusplusplus(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=exp(-3.273424-9.222532e-3*logP-2.546540e-3*pow(logP,2.0)-6.142466e-4*pow(logP,3.0)-6.803461e-5*pow(logP,4.0)-1.480622e-6*pow(logP,5.0));
  a3=-exp(-3.227410-4.108171e-3*logP-6.841752e-4*pow(logP,2.0)-3.928651e-5*pow(logP,3.0));
  a2=-a3-a1;
  c1=exp(1.114079e+1+6.128099e-2*logP+1.305781e-3*pow(logP,2.0)-4.745385e-5*pow(logP,3.0)-1.294845e-5*pow(logP,4.0)-6.416314e-7*pow(logP,5.0));
  cm=exp(1.097387e+1+5.385207e-2*logP+3.454294e-5*pow(logP,2.0)-9.334055e-5*pow(logP,3.0));
  c2=exp(1.133963e+1+5.445065e-2*logP-3.976441e-4*pow(logP,2.0)+1.251159e-4*pow(logP,3.0));
  c3=exp(1.473199e+1-3.158041e-1*logP-3.070674e-2*pow(logP,2.0)+3.070674e-2*pow(logP,3.0));
  Delta1=exp(9.124558+1.015232e-1*logP-1.452067e-3*pow(logP,2.0)-4.363441e-4*pow(logP,3.0)-9.737843e-6*pow(logP,4.0)+1.643326e-6*pow(logP,5.0));
  Deltam=exp(9.008289+5.266326e-2*logP-2.558320e-4*pow(logP,2.0)+3.532844e-5*pow(logP,3.0));
  Delta2=exp(9.165912+3.362575e-2*logP+1.118630e-3*pow(logP,2.0)-3.084012e-4*pow(logP,3.0)-7.665827e-5*pow(logP,4.0));
  Delta3=exp(1.306288e+1-3.228563e-12*logP-3.275522e-2*pow(logP,2.0)+6.750116e-3*pow(logP,3.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_NO(double T, double P){
  double sum,logP,a1,a2,a3,c1,c2,c3,Delta1,Delta2,Delta3,sigma1,sigma2,sigma3;
  
  logP=log(P);
  a1=exp(-2.397641+9.644207e-2*logP);
  a3=-exp(-2.923272+1.671984e-1*logP);
  a2=-a3-a1;
  c1=exp(7.942600+2.917164e-2*logP+6.775381e-4*pow(logP,2.0)+2.209082e-5*pow(logP,3.0));
  c2=exp(8.274503+6.655621e-2*logP+2.214534e-3*pow(logP,2.0)+3.856329e-5*pow(logP,3.0));
  c3=exp(8.364477+7.365241e-2*logP+2.771836e-3*pow(logP,2.0)-5.013391e-6*pow(logP,3.0)-5.293616e-6*pow(logP,4.0));
  Delta1=exp(6.780323+6.029139e-2*logP+4.276063e-4*pow(logP,2.0));
  Delta2=exp(6.495225+7.930874e-2*logP-1.952605e-3*pow(logP,2.0)-7.384374e-4*pow(logP,3.0)-5.231985e-5*pow(logP,4.0));
  Delta3=exp(7.549495+9.399569e-2*logP);
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sum=a1*sigma1+a2*sigma2+a3*sigma3;
  return(fmax(0.0,sum));
}

double _chi_NOplus(double T, double P){
  double sum,logP,a1,a2,a3,a4,c1,c2,c3,c4,Delta1,Delta2,Delta3,Delta4,sigma1,sigma2,sigma3,sigma4;
  
  logP=log(P);
  a1=exp(-7.135266+4.617651e-2*logP-7.097386e-4*pow(logP,2.0));
  a2=exp(-8.753925+1.392942e-1*logP+1.317873e-2*pow(logP,2.0));
  a4=-exp(-8.772596+1.382471e-1*logP-1.513718e-3*pow(logP,2.0)+1.822779e-3*pow(logP,3.0)+5.774867e-5*pow(logP,4.0));
  a3=-a4-a2-a1;
  c1=exp(8.740893+4.144123e-2*logP+3.456197e-4*pow(logP,2.0));
  c2=exp(8.817743+4.865084e-2*logP+7.358462e-4*pow(logP,2.0));
  c3=exp(8.899564+6.228872e-2*logP+1.910295e-3*pow(logP,2.0)+5.292903e-5*pow(logP,3.0));
  c4=exp(9.221935+8.005371e-2*logP+3.728793e-3*pow(logP,2.0)-1.235847e-4*pow(logP,3.0)-6.058282e-6*pow(logP,4.0));
  Delta1=exp(6.996599+6.789593e-2*logP+1.320085e-3*pow(logP,2.0)+2.143434e-5*pow(logP,3.0)+6.597691e-6*pow(logP,4.0));
  Delta2=exp(6.260938+9.417073e-2*logP+7.841151e-3*pow(logP,2.0));
  Delta3=exp(7.246371+1.012940e-1*logP+4.389279e-3*pow(logP,2.0)-2.344414e-5*pow(logP,3.0)-1.533963e-5*pow(logP,4.0));
  Delta4=exp(7.940040+1.021609e-1*logP+5.411563e-3*pow(logP,2.0)+1.592304e-4*pow(logP,3.0)-7.583651e-5*pow(logP,4.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sigma4=1.0/(1.0+exp(-2.0*(T-c4)/Delta4));
  sum=a1*sigma1+a2*sigma2+a3*sigma3+a4*sigma4;
  return(fmax(0.0,sum));
}

double _chi_eminus_improved(double T, double P){ // improved : 10.1140/epjd/e2011-20424-5
  double sum,logP,a1,a2,a3,a3b,a4,c0,c1,c2,c3,c4,Delta0,Delta1,Delta2,Delta3,Delta4,sigma0,sigma1,sigma2,sigma3,cm,Deltam,sigmam,gamma4;
  double delta0,delta1,delta2,delta3,delta4,delta5,delta6;
  
  P=P*1.0132515; //atm to bar
  logP=log(P);
  delta0=-1.453520e+002-2.449447e+000*logP-1.115521e-001*pow(logP,2.0)+4.815904e-002*pow(logP,3.0)+4.816233e-003*pow(logP,4.0)-4.844801e-004*pow(logP,5.0)-6.055194e-005*pow(logP,6.0);
  delta1=1.369547e-001-2.858663e-003*logP+3.205470e-005*pow(logP,2.0)-8.338958e-005*pow(logP,3.0)-5.701341e-006*pow(logP,4.0)+8.585860e-007*pow(logP,5.0)+8.583472e-008*pow(logP,6.0);
  delta2=-6.345758e-005-1.679107e-006*logP+4.450866e-008*pow(logP,2.0)+4.948604e-008*pow(logP,3.0)+2.313164e-009*pow(logP,4.0)-5.251591e-010*pow(logP,5.0)-4.461556e-011*pow(logP,6.0);
  delta3=1.617445e-008+4.919459e-010*logP-2.922150e-011*pow(logP,2.0)-1.375302e-011*pow(logP,3.0)-4.284662e-013*pow(logP,4.0)+1.504247e-013*pow(logP,5.0)+1.156213e-014*pow(logP,6.0);
  delta4=-2.299085e-012-7.370074e-014*logP+6.536403e-015*pow(logP,2.0)+1.931609e-015*pow(logP,3.0)+3.938597e-017*pow(logP,4.0)-2.191122e-017*pow(logP,5.0)-1.631548e-018*pow(logP,6.0);
  delta5=1.702203e-016+5.392941e-018*logP-6.228118e-019*pow(logP,2.0)-1.322524e-019*pow(logP,3.0)-1.891956e-021*pow(logP,4.0)+1.569686e-021*pow(logP,5.0)+1.216511e-022*pow(logP,6.0);
  delta6=-5.100309e-021-1.529637e-022*logP+2.147344e-023*pow(logP,2.0)+3.497704e-024*pow(logP,3.0)+4.747413e-026*pow(logP,4.0)-4.387925e-026*pow(logP,5.0)-3.786008e-027*pow(logP,6.0);
  a1=exp(-3.897330e-001+4.088303e-003*logP+2.004021e-003*pow(logP,2.0)+6.989005e-004*pow(logP,3.0)+6.796657e-005*pow(logP,4.0)-1.460238e-006*pow(logP,5.0));
  a2=exp(-1.599518e+000-3.681243e-002*logP-1.499672e-002*pow(logP,2.0)-4.875898e-003*pow(logP,3.0)-9.278204e-005*pow(logP,4.0)+8.792226e-005*pow(logP,5.0)+1.273088e-005*pow(logP,6.0));
  a3b=exp(-3.465436e-001-2.831472e-003*logP-1.021467e-003*pow(logP,2.0)-7.753035e-005*pow(logP,3.0));
  a3=a3b-a2-a1;
  a4=exp(-5.123507e+000+3.850431e-001*logP+6.380724e-003*pow(logP,2.0)-1.088865e-002*pow(logP,3.0)+1.189765e-003*pow(logP,4.0));
  c0=8.047945e+003+2.117097e+002*logP-6.978618e+000*logP*logP;
  c1=exp(9.510842e+000+6.332838e-002*logP+5.805702e-004*pow(logP,2.0)-6.640588e-005*pow(logP,3.0));
  cm=exp(6.343867e+000+1.473478e+000*logP-2.628976e-001*pow(logP,2.0)+2.653667e-002*pow(logP,3.0)-1.170989e-003*pow(logP,4.0));
  c2=exp(1.025313e+001+6.613035e-002*logP+2.106960e-003*pow(logP,2.0)+1.249059e-004*pow(logP,3.0)-3.254728e-006*pow(logP,4.0)-1.073094e-006*pow(logP,5.0)-4.149968e-007*pow(logP,6.0)-4.918145e-008*pow(logP,7.0));
  c3=exp(1.009244e+001+5.691765e-002*logP+2.642057e-003*pow(logP,2.0)+3.297719e-005*pow(logP,3.0));
  c4=exp(9.247533e+000+1.007851e-001*logP-1.058914e-002*pow(logP,2.0)+9.089247e-004*pow(logP,3.0));
  Delta0=5.000000e+001+6.863114e-002*logP;
  Delta1=exp(7.970326e+000+1.336309e-001*logP+3.159430e-003*pow(logP,2.0)+3.346003e-004*pow(logP,3.0)+5.512825e-005*pow(logP,4.0));
  Deltam=exp(1.029159e+001+3.502366e-002*logP-1.043994e-002*pow(logP,2.0)-7.498040e-004*pow(logP,3.0)+1.464646e-004*pow(logP,4.0)+1.031691e-005*pow(logP,5.0)-3.878009e-007*pow(logP,6.0));
  Delta2=exp(8.461864e+000+1.033435e-001*logP-6.800325e-003*pow(logP,2.0)-2.171111e-003*pow(logP,3.0)+8.042855e-005*pow(logP,4.0)+3.126866e-005*pow(logP,5.0)+3.548083e-006*pow(logP,6.0)+1.732832e-007*pow(logP,7.0));
  Delta3=exp(9.041428e+000+9.809302e-002*logP+1.899235e-003 *pow(logP,2.0)-1.329754e-004*pow(logP,3.0)-2.357106e-005*pow(logP,4.0));
  Delta4=exp(8.132477e+000+9.307459e-002 *logP-5.524780e-003*pow(logP,2.0)+2.540886e-004*pow(logP,3.0)+1.232437e-004*pow(logP,4.0));
  sigma0=1.0/(1.0+exp(-2.0*(T-c0)/Delta0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  gamma4=exp(-(T-c4)*(T-c4)/(Delta4*Delta4));
  sum=(1.0-sigma0)*exp(delta0+delta1*T+delta2*T*T+delta3*T*T*T+delta4*T*T*T*T+delta5*T*T*T*T*T+delta6*T*T*T*T*T*T)+a1*sigma0*sigma1*sigmam+sigma0*(a2*sigma2+a3*sigma3-a4*gamma4);
  return(sum);
}

double _chi_eminus_default(double T, double P){
  double sum,logP,a1,a2,a3,a4,a5,a6,c1,c2,c3,c4,c5,c6,Delta1,Delta2,Delta3,Delta4,Delta5,Delta6,sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,cm,Deltam,sigmam;
  
  logP=log(P);
  a1=exp(-3.932487e-1+7.116035e-4*logP+4.083493e-4*pow(logP,2.0)+3.307562e-4*pow(logP,3.0)+2.215248e-5*pow(logP,4.0)-4.020145e-6*pow(logP,5.0));
  a2=exp(-1.599518-3.681243e-2*logP-1.499672e-2*pow(logP,2.0)-4.875898e-3*pow(logP,3.0)-9.278204e-5*pow(logP,4.0)+8.792226e-5*pow(logP,5.0)+1.273088e-5*pow(logP,6.0));
  a3=exp(-3.031217-1.236964e-2*logP+4.999807e-3*pow(logP,2.0)+4.130827e-4*pow(logP,3.0)-5.879976e-5*pow(logP,4.0)-5.643378e-6*pow(logP,5.0)-2.118479e-7*pow(logP,6.0)-8.835667e-8*pow(logP,7.0));
  a4=exp(-3.096101+5.690833e-2*logP+1.063005e-2*pow(logP,2.0)+8.066239e-4*pow(logP,3.0));
  a6=-exp(-3.465436e-1-2.831472e-3*logP-1.021467e-3*pow(logP,2.0)-7.753035e-5*pow(logP,3.0));
  a5=-a6-a2-a1;
  c1=exp(9.514823+6.426626e-2*logP+3.538392e-4*pow(logP,2.0)-1.093881e-4*pow(logP,3.0));
  cm=exp(6.343867+1.473478*logP-2.628976e-1*pow(logP,2.0)+2.653667e-2*pow(logP,3.0)-1.170989e-3*pow(logP,4.0));
  c2=exp(1.025313e+1+6.613035e-2*logP+2.106960e-3*pow(logP,2.0)+1.249059e-4*pow(logP,3.0)-3.254728e-6*pow(logP,4.0)-1.073094e-6*pow(logP,5.0)-4.149968e-7*pow(logP,6.0)-4.918145e-8*pow(logP,7.0));
  c3=exp(1.074247e+1+6.026184e-2*logP+6.834881e-4*pow(logP,2.0)+1.412968e-6*pow(logP,3.0));
  c4=exp(1.106632e+1+5.734452e-2*logP+1.326880e-3*pow(logP,2.0)+4.870977e-5*pow(logP,3.0));
  c5=exp(1.009244e+1+5.691765e-2*logP+2.642057e-3*pow(logP,2.0)+3.297719e-5*pow(logP,3.0));
  c6=exp(1.219498e+2-3.565001*logP+7.046916e-1*pow(logP,2.0)+3.062083e-1*pow(logP,3.0)-2.940975e-2*pow(logP,4.0));
  Delta1=exp(7.931006+1.174176e-1*logP-5.369256e-4*pow(logP,2.0)-1.640676e-4*pow(logP,3.0)+2.876393e-5*pow(logP,4.0));
  Deltam=exp(1.029159e+1+3.502366e-2*logP-1.043994e-2*pow(logP,2.0)-7.498040e-4*pow(logP,3.0)+1.464646e-4*pow(logP,4.0)+1.031691e-5*pow(logP,5.0)-3.878009e-7*pow(logP,6.0));
  Delta2=exp(8.461864+1.033435e-1*logP-6.800325e-3*pow(logP,2.0)-2.171111e-3*pow(logP,3.0)+8.042855e-5*pow(logP,4.0)+3.126866e-5*pow(logP,5.0)+3.548083e-6*pow(logP,6.0)+1.732832e-7*pow(logP,7.0));
  Delta3=exp(8.457103+1.570495e-1*logP+2.577271e-2*pow(logP,2.0)-4.699755e-4*pow(logP,3.0)-7.340190e-4*pow(logP,4.0)-1.521958e-6*pow(logP,5.0)+7.337098e-6*pow(logP,6.0)-1.937258e-8*pow(logP,7.0));
  Delta4=exp(9.134358+1.817063e-1*logP+8.463508e-3*pow(logP,2.0));
  Delta5=exp(9.041428+9.809302e-2*logP+1.899235e-3*pow(logP,2.0)-1.329754e-4*pow(logP,3.0)-2.357106e-5*pow(logP,4.0));
  Delta6=exp(1.163952e+2-3.232407*logP+6.981116e-1*pow(logP,2.0)+2.997466e-1*pow(logP,3.0)-2.749064e-2*pow(logP,4.0));
  sigma1=1.0/(1.0+exp(-2.0*(T-c1)/Delta1));
  sigmam=1.0/(1.0+exp(-2.0*(T-cm)/Deltam));
  sigma2=1.0/(1.0+exp(-2.0*(T-c2)/Delta2));
  sigma3=1.0/(1.0+exp(-2.0*(T-c3)/Delta3));
  sigma4=1.0/(1.0+exp(-2.0*(T-c4)/Delta4));
  sigma5=1.0/(1.0+exp(-2.0*(T-c5)/Delta5));
  sigma6=1.0/(1.0+exp(-2.0*(T-c6)/Delta6));
  sum=a1*sigmam*sigma1+a2*sigma2+a3*sigma3+a4*sigma4+a5*sigma5+a6*sigma6;
  return(fmax(0.0,sum));
}

double _chi_eminus(double T, double P){
  if(T<100.0 || T>20000.0) return(_chi_eminus_default(T,P));
  else return(_chi_eminus_improved(T,P));
}


int main (int argc, char **argv) {
  double T,P,Tbegin,Tend,Tincrem;
  bool validOptions = TRUE;
  
  if (chkarg (argc, argv, "-Tbegin") == 0)
    validOptions = FALSE;       //check if argument flags exist
  if (chkarg (argc, argv, "-Tend") == 0)
    validOptions = FALSE;
  if (chkarg (argc, argv, "-P") == 0)
    validOptions = FALSE;

  if (chkarg (argc, argv, "-Tbegin") != 0) {
    if (sscanf (argv[chkarg (argc, argv, "-Tbegin") + 1], "%lg", &Tbegin) != 1)
      validOptions = FALSE;
  }
  if (chkarg (argc, argv, "-Tend") != 0) {
    if (sscanf (argv[chkarg (argc, argv, "-Tend") + 1], "%lg", &Tend) != 1)
      validOptions = FALSE;
  }
  if (chkarg (argc, argv, "-Tincrem") != 0) {
    if (sscanf (argv[chkarg (argc, argv, "-Tincrem") + 1], "%lg", &Tincrem) != 1)
      validOptions = FALSE;
  }
  if (chkarg (argc, argv, "-P") != 0) {
    if (sscanf (argv[chkarg (argc, argv, "-P") + 1], "%lg", &P) != 1)
      validOptions = FALSE;
  }

  
  if (!validOptions) {
    fprintf (stderr, "\nFlags:\n\n"
             "Flag               \tArg                                     \tArg Type  \tRequired?\n"
             "-------------------------------------------------------------------------------------------------\n"
             "-Tbegin            \tbeginning of temperature range, Kelvin  \tdouble       \tY\n"
             "-Tend              \tend of temperature range, Kelvin        \tdouble       \tY\n"
             "-Tincrem           \tincrement in temperature, Kelvin        \tdouble       \tY\n"
             "-P                 \tpressure, atm                           \tdouble       \tY\n"
             "Eg: \n"
             "equilairplasma -Tbegin 1000.0 -Tend 60000.0 -Tincrem 1000.0 -P 1.0 \nwill output the equilibrium air-plasma molar fractions over the temperature range 1000-60000 K\n"
             "for a pressure of 1 atm to stdout .\n");
    exit (1);
  }

  sscanf (argv[chkarg (argc, argv, "-Tbegin") + 1], "%lg", &Tbegin);
  sscanf (argv[chkarg (argc, argv, "-Tend") + 1], "%lg", &Tend);
  sscanf (argv[chkarg (argc, argv, "-Tincrem") + 1], "%lg", &Tincrem);
  sscanf (argv[chkarg (argc, argv, "-P") + 1], "%lg", &P);

  printf("Pressure = %g atm. Temperature range is %g to %g K.\n",P,Tbegin,Tend);
  if(P<=1E-2 || P>=1E2) fprintf (stderr, "Ensure that the pressure lies between 0.01 and 100 atm. Aborting.\n");
  for(T=Tbegin;T<=Tend;T=T+Tincrem){
    if(T<=49.9 || T>=60000.1) fprintf (stderr, "Ensure that the temperatures lie between 50 and 60000 K. Aborting.\n");
  }
  if (Tbegin==Tend){
    fprintf(stdout,"T          = %E K\n",T);
    fprintf(stdout,"chi_N2     = %E \n",_chi_N2(T,P));
    fprintf(stdout,"chi_N2+    = %E \n",_chi_N2plus(T,P));
    fprintf(stdout,"chi_N      = %E \n",_chi_N(T,P));
    fprintf(stdout,"chi_N+     = %E \n",_chi_Nplus(T,P));
    fprintf(stdout,"chi_N++    = %E \n",_chi_Nplusplus(T,P));
    fprintf(stdout,"chi_N+++   = %E \n",_chi_Nplusplusplus(T,P));
    fprintf(stdout,"chi_N++++  = %E \n",_chi_Nplusplusplusplus(T,P));
    fprintf(stdout,"chi_O2     = %E \n",_chi_O2(T,P));
    fprintf(stdout,"chi_O2+    = %E \n",_chi_O2plus(T,P));
    fprintf(stdout,"chi_O2-    = %E \n",_chi_O2minus(T,P));
    fprintf(stdout,"chi_O      = %E \n",_chi_O(T,P));
    fprintf(stdout,"chi_O-     = %E \n",_chi_Ominus(T,P));
    fprintf(stdout,"chi_O+     = %E \n",_chi_Oplus(T,P));
    fprintf(stdout,"chi_O++    = %E \n",_chi_Oplusplus(T,P));
    fprintf(stdout,"chi_O+++   = %E \n",_chi_Oplusplusplus(T,P));
    fprintf(stdout,"chi_O++++  = %E \n",_chi_Oplusplusplusplus(T,P));
    fprintf(stdout,"chi_NO     = %E \n",_chi_NO(T,P));
    fprintf(stdout,"chi_NO+    = %E \n",_chi_NOplus(T,P));
    fprintf(stdout,"chi_e-     = %E \n",_chi_eminus(T,P));
  } else {
    fprintf(stdout,"T,K\tN2\tN2+\tN\tN+\tN++\tN+++\tN++++\tO2\tO2+\tO2-\tO\tO-\tO+\tO++\tO+++\tO++++\tNO\tNO+\te-\n");
    for(T=Tbegin;T<=Tend;T=T+Tincrem){
      fprintf (stdout,"%.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E %.14E\n",T,_chi_N2(T,P),_chi_N2plus(T,P),_chi_N(T,P),_chi_Nplus(T,P),_chi_Nplusplus(T,P),_chi_Nplusplusplus(T,P),_chi_Nplusplusplusplus(T,P),_chi_O2(T,P),_chi_O2plus(T,P),_chi_O2minus(T,P),_chi_O(T,P),_chi_Ominus(T,P),_chi_Oplus(T,P),_chi_Oplusplus(T,P),_chi_Oplusplusplus(T,P),_chi_Oplusplusplusplus(T,P),_chi_NO(T,P),_chi_NOplus(T,P),_chi_eminus(T,P));
    }
  }
  
  
  return (0);
}
