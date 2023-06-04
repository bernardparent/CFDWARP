// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2023 Bernard Parent

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

#include <exm.h>
#include "share.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Computer Physics Communications 175 (2006) 545-552   (for non-magnetized plasma only)



int main ( int argc, char **argv ) {
  bool VALIDOPTIONS=TRUE;
  char *options;
  int RET;
  double e,ame,eps0,ccc,x,Ne,fcoll,anu,fom,omega,omegap,rab1,rab2,rab3,rab4,rab5;
  double b,absorp1,absorp2,absorp3,att1,att2,att3;
  double muen,mue,muei,lnlambda,Nn,Ni,Te;


  options = NULL;
  
  if (process_flag_double(argc, argv, "-t", &x)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-f", &fom)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-Ne", &Ne)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-Nn", &Nn)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-Te", &Te)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-Ni", &Ni)!=2) Ni=Ne;

  if ( !VALIDOPTIONS ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag   \tArg                              \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-t     \t<slab thickness [m]>             \tdouble   \tY\n"
              "-Te    \t<electron temperature [K]>       \tdouble   \tY\n"
              "-f     \t<wave frequency [Hz]>            \tdouble   \tY\n"
              "-Ne    \t<electron number density [1/m3]> \tdouble   \tY\n"
              "-Nn    \t<neutrals number density [1/m3]> \tdouble   \tY\n"
              "-Ni     \t<ion number density [1/m3]>     \tdouble   \tN\n"
               "\n\n"
              "Eg.: ./atten -t 0.01 -Ne 1e18 -Nn 1e23 -Te 10000 -f 1e9\n"
              "\n\n");
    exit (EXIT_FAILURE);
  }
  
  RET = find_remaining_options ( argc, argv, &options );
  if ( RET >= 1 ) {
    fprintf ( stderr, "\n\nThe following command line options could not be processed:\n%s\n\n", options );
    exit (EXIT_FAILURE);
  }
    
	e=1.6e-19; // elementary charge in Coulomb
	ame=9.1e-31; // electron mass in kg
	eps0=8.85e-12;  // permittivity of free space
	ccc=3.0e8;


  // Collision frequency (linear, GHz converted to 1/s)
  fcoll=100.0;
  anu= fcoll * 1.0e9;
    
  // find collision frequency anu from mobility in 1/s */
  lnlambda=23.0-log(sqrt(Ne/1e6)*pow(Te/11600.0,-1.5));
  muen=3.74e19*exp(33.5/sqrt(log(max(300.0,Te))))/Nn;
  muei=9.5e16*(pow(Te,1.5))/Ni/lnlambda;
  mue=1.0/(1.0/muen+1.0/muei);
  anu=e/(ame*mue);
    
  omega= 2.*pi* fom ;
  // Plasma frequency (circular, 1/s)
  omegap=sqrt(Ne*sqr(e)/ame/eps0);

  rab1=sqr(omega*omega) + anu*anu*sqr(omega);
  rab2=sqr(omegap) * sqr(omega) / rab1;
  rab3=anu*omega*sqr(omegap) / rab1;
  rab4=sqrt( sqr(1.0-rab2) + sqr(rab3) );
  rab5=1.0-rab2;

  b=sqr(omega/ccc) * (rab4-rab5)/2.0;

  absorp1=sqrt(b);
  absorp2=sqr(omegap)/ccc/anu;
  absorp3=sqr(omegap)*anu/ccc/sqr(omega);

  att1=exp(-absorp1*x);
  att2=exp(-absorp2*x);
  att3=exp(-absorp3*x);


  printf("Ne                     = %E 1/m3\n",Ne);
  printf("Ni                     = %E 1/m3\n",Ni);
  printf("Nn                     = %E 1/m3\n",Nn);
  printf("Te                     = %E K\n",Te);
  printf("plasma layer thickness = %E m\n",x);
  printf("muei                   = %E m2/Vs\n",muei);
  printf("muen                   = %E m2/Vs\n",muen);
  printf("mue                    = %E m2/Vs\n",mue);
  printf("e- collision frequency = %E Hz\n",anu);
  printf("EM wave frequency      = %E Hz\n",fom);
  printf("I/I0, att1^2           = %E\n",sqr(att1));
  printf("att2^2                 = %E\n",sqr(att2));
  printf("att3^2                 = %E\n",sqr(att3));


  return ( EXIT_SUCCESS );
}

