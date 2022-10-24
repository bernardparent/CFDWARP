// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2001,2022 Bernard Parent

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



typedef struct {
  double TyoverTx;
  double gamma;
} arg2_t;


static double _errfunct2 ( void *arg, double Ms ) {
  double tmp, TyoverTx, gamma;
  gamma=( ( arg2_t * ) arg )->gamma;
  TyoverTx=(2.0*gamma/(gamma+1.0)*Ms*Ms-(gamma-1.0)/(gamma+1.0))*((gamma-1.0)/(gamma+1.0)+2.0/(gamma+1.0)/Ms/Ms);
  tmp = TyoverTx - ( ( arg2_t * ) arg )->TyoverTx;
  return ( tmp );
}


void find_Ms ( double TyoverTx, double gamma, double *Ms ) {
  long IFLAG;
  arg2_t arg;
  arg.TyoverTx = TyoverTx;
  arg.gamma=gamma;
  *Ms = EXM_find_root_zero_in ( &( _errfunct2 ), &arg, 1.0e0, 20.0e0, 1.0e-12, 1.0e-12, &IFLAG );
  if ( IFLAG == 4 ) {
    fprintf ( stderr, "couldn't find root in find_Ms, exiting..\n" );
    exit ( 1 );
  }
}


int main ( int argc, char **argv ) {
  bool VALIDOPTIONS=TRUE;
  double tmp, Ms, U0, P0, M0, T0, H1, W0, gamma, Rgas, Pdyn, U, P, T, H, W, M,
         mdot, altitude, PyoverPx, TyoverTx, Te, cp, Tstag, Pstag, Tstag0, Pstag0, Wold, Mold, phi, delta;
  long numshock,cnt;
  char *options;
  int RET;
  options = NULL;
  
  if (process_flag_double(argc, argv, "-M", &M0)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-W", &W0)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-Pdyn", &Pdyn)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_double(argc, argv, "-Te", &Te)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_long(argc, argv, "-shocks", &numshock)!=2) VALIDOPTIONS=FALSE;  
  if (process_flag_double(argc, argv, "-gamma", &tmp)==2) gamma=tmp; else gamma=1.4;
  if (process_flag_double(argc, argv, "-Rgas", &tmp)==2) Rgas=tmp; else Rgas=287.0;

  if ( !VALIDOPTIONS ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag   \tArg                              \tArg Type \tRequired? \n"
              "----------------------------------------------------------------------------\n"
              "-M     \t<flight Mach number [m/s]>       \tdouble   \tY\n"
              "-Pdyn  \t<flight dynamic pressure [Pa]>   \tdouble   \tY\n"
              "-Te    \t<temperature at inlet exit [K]>  \tdouble   \tY\n"
              "-shocks\t<number of shocks>               \tdouble   \tY\n"
              "-W     \t<height of the inlet [m]>        \tdouble   \tY\n"
              "-gamma \t<ratio of specific heats>        \tdouble   \tN\n"
              "-Rgas  \t<gas constant [J/kgK]>           \tdouble   \tN\n"
               "\n\n"
              "Eg.: ./inlet -M 7.0 -Pdyn 66000 -Te 1000 -shocks 2 -W 1.0\n"
              "\n\n");
    exit (EXIT_FAILURE);
  }
  
  RET = find_remaining_options ( argc, argv, &options );
  if ( RET >= 1 ) {
    fprintf ( stderr, "\n\nThe following command line options could not be processed:\n%s\n\n", options );
    exit (EXIT_FAILURE);
  }
    
  cp=gamma/(gamma-1.0)*Rgas;
  P0 = 2.0e0 * Pdyn / gamma / pow ( M0, 2.0 );
  FindTandHfromPstdatm ( P0, &altitude, &T0 );
  U0 = M0 * sqrt ( gamma * Rgas * T0 );
  H1=0.5*U0*U0+cp*T0;    
  TyoverTx=pow(Te/T0,1.0/numshock); 
  find_Ms ( TyoverTx, gamma, &Ms) ;
  PyoverPx=2.0*gamma/(gamma+1.0)*Ms*Ms-(gamma-1.0)/(gamma+1.0);
  M0=U0/sqrt(gamma*Rgas*T0);
  Pstag0=P0*pow(1.0+(gamma-1.0)/2.0*M0*M0,gamma/(gamma-1.0));
  Tstag0=T0*(1.0+(gamma-1.0)/2.0*M0*M0);

  mdot=W0*U0*P0/T0/Rgas;
  fprintf ( stdout,
             "Altitude: %E m\n"
             "Ms      : %E\n"
             "mdot    : %E kg/sm\n"
             "gamma   : %E\n"
             "Rgas    : %E J/kgK\n",
             altitude, Ms, mdot, gamma, Rgas );
    
  fprintf (stdout, "\n"
              "W0      : %E m\n"
              "P0      : %E Pa\n"
              "T0      : %E K\n"
              "U0      : %E m/s\n"
              "M0      : %E\n"
              "Tstag0  : %E\n"
              "Pstag0  : %E Pa\n",
             W0, P0, T0, U0, M0, Tstag0,Pstag0 );
  M=M0;
  W=W0;
  for (cnt=1; cnt<=numshock; cnt++){
    Mold=M;
    Wold=W;
    // total enthalpy is conserved through the shocks
    H=H1;
    // all shocks have equal strength so T2/T0=T3/T2=T4/T3 etc
    T=pow(TyoverTx,(double)cnt)*T0;
    P=pow(PyoverPx,(double)cnt)*P0;
    U=sqrt((H-cp*T)*2.0);
    M=U/sqrt(gamma*Rgas*T);
    Pstag=P*pow(1.0+(gamma-1.0)/2.0*M*M,gamma/(gamma-1.0));
    Tstag=T*(1.0+(gamma-1.0)/2.0*M*M);
    // mass conservation yields flow height W
    W = W0 * ( P0 / T0 * U0 ) / ( P / T * U );
    // find shock angle
    phi=asin(Ms/Mold);
    // find turning angle
    delta=phi-asin(W/Wold*sin(phi));
    fprintf (stdout, "\n"
              "W%ld      : %E m\n"
              "P%ld      : %E Pa\n"
              "T%ld      : %E K\n"
              "U%ld      : %E m/s\n"
              "M%ld      : %E\n"
              "Pstag%ld  : %E Pa\n"
              "Tstag%ld  : %E K\n"
              "phi%ld%ld   : %E deg\n"
              "delta%ld%ld : %E deg\n",
             cnt,W, cnt,P, cnt,T, cnt,U, cnt,M, cnt,Pstag, cnt,Tstag, cnt-1,cnt,phi/pi*180.0, cnt-1,cnt,delta/pi*180.0 );
  }

  return ( EXIT_SUCCESS );
}
