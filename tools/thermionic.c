// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2019 Bernard Parent

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

#include "share.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define A0 1.20173E6  // A/(m2K2)
#define kB 8.61733E-5 // Boltzmann constant in eV/K
#define echarge 1.6e-19 // coulombs



int main (int argc, char **argv) {
  double T,W,lambdaR,tmp;
  double Ag;
  bool VALIDOPTIONS = TRUE;
  char *options;
  int RET;
  options = NULL;
  
  if (process_flag_double(argc, argv, "-T", &T)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-W", &W)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-lambdaR", &tmp)==2) lambdaR=tmp; else lambdaR=0.5;

  Ag=A0*lambdaR;


  if (!VALIDOPTIONS) {
    fprintf (stderr, "\nPROBLEM WITH FLAGS\n\n\nRequired and Optional Flags:\n\n"
             "Flag      Arg                     Arg Type     Required?\n"
             "----------------------------------------------------------\n"
             "-W        work function in eV     double       Y\n"
             "-T        temperature in K        double       Y\n"
             "-lambdaR  material constant       double       N\n"
             "          default: 0.5                          \n"
             "Eg: \n"
             "./thermionic -W 3.0 -T 2800.0\nwill give the thermionic current of a material with a work function of 3 eV and a temperature of 2800 K.\n");
    exit (EXIT_FAILURE);
  }
  
  RET = find_remaining_options ( argc, argv, &options );
  if ( RET >= 1 ) {
    fprintf ( stderr, "\n\nThe following command line options could not be processed:\n%s\n\n", options );
    exit (EXIT_FAILURE);
  } else {
    printf("T             = %E K\n",T);
    printf("Work function = %E eV\n",W);
    printf("J_thermionic  = %E A/m2\n",Ag*T*T*exp(-W/(kB*T)));
  }
  free ( options );
  return(EXIT_SUCCESS);
}
