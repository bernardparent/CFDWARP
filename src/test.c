// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2001-2002 Bernard Parent

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

#include <src/common.h>
#include <model/_model.h>
#include <model/.active/model.h>
#include <model/thermo/_thermo.h>
#include <src/control.h>
#include <cycle/share/cycle_share.h>
#include <cycle/share/res_share.h>
#include <stdarg.h>


#define VERBOSE FALSE
#define lengthcol1 17
#define lengthcol2 14

void seed_random ( void ) {
  time_t timer;
  time ( &timer );
  srand ( timer );
}

void printf_ijk_position_of_node_l ( gl_t * gl, long l ) {
  long i, j, k;
  char strtmp[100], str[100];
  find_ijk_from_l ( gl, l, &i, &j, &k );
  sprintf ( str, "i=%ld", i - gl->domain_all.is + 1 );

#ifdef _2DL
  sprintf ( strtmp, " j=%ld", j - gl->domain_all.js + 1 );
  strcat ( str, strtmp );
#endif

#ifdef _3DL
  sprintf ( strtmp, " k=%ld", k - gl->domain_all.ks + 1 );
  strcat ( str, strtmp );
#endif
  printf ( "%s", str );
}

void test_A ( np_t * np, long l, gl_t * gl, long theta ) {
  sqmat_t A, Anum;

  printf ( "\n" );
  printf ( "The A jacobian should be equal to the numerically determined one.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );
  find_dFstar_dUstar ( np[l], gl, theta, A );
  find_numerical_jacobian ( np, l, gl, &find_Fstar, theta, Anum );

  printf ( "A jacobian (analytical):\n" );
  display_matrix ( A );
  printf ( "A jacobian (numerically determined):\n" );
  display_matrix ( Anum );
}

void test_A_from_jacvars ( np_t * np, long l, gl_t * gl, long theta ) {
  sqmat_t A2, Anum;
  jacvars_t jacvars;
  metrics_t metrics;

  printf ( "\n" );
  printf ( "Testing A jacobian from jacvars compared to numerically determined one.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  find_metrics_at_node ( np, gl, l, theta, &metrics );
  find_jacvars ( np[l], gl, metrics, theta, &jacvars );
  find_A_from_jacvars ( jacvars, metrics, A2 );
  find_numerical_jacobian ( np, l, gl, &find_Fstar, theta, Anum );

  printf ( "A jacobian (analytical):\n" );
  display_matrix ( A2 );
  printf ( "A jacobian (numerically determined):\n" );
  display_matrix ( Anum );
}

void test_dG_dUstar ( np_t * np, long l, gl_t * gl ) {
  sqmat_t Bk, Bknum;

  printf ( "\n" );
  printf ( "Testing dG/dUstar jacobian compared to numerically determined one.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  find_dG_dUstar ( np[l], gl, Bk );
  find_numerical_jacobian_2 ( np, l, gl, &find_G, Bknum );

  printf ( "dG/dUstar jacobian (analytical):\n" );
  display_matrix ( Bk );
  printf ( "dG/dUstar jacobian (numerical):\n" );
  display_matrix ( Bknum );
}

void test_dUstar_dUprime ( np_t * np, long l, gl_t * gl ) {
  sqmat_t dUstardUprime, dUprimedUstarnum, dUstardUprimenum;
  metrics_t metrics;
  jacvars_t jacvars;
  long theta;

  printf ( "\n" );
  printf ( "Testing dUstar/dUprime jacobian compared to numerically determined one.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  theta = 0;                    //the value given to theta shouldn't change the jacobian
  find_metrics_at_node ( np, gl, l, theta, &metrics );
  find_jacvars ( np[l], gl, metrics, theta, &jacvars );

  find_dUstar_dUprime_from_jacvars ( jacvars, metrics, dUstardUprime );
  find_numerical_jacobian_2 ( np, l, gl, &find_Uprime, dUprimedUstarnum );
  invert_matrix ( dUprimedUstarnum, dUstardUprimenum );

  printf ( "dUstar/dUprime jacobian (analytical):\n" );
  display_matrix ( dUstardUprime );
  printf ( "dUstar/dUprime jacobian (numerical):\n" );
  display_matrix ( dUstardUprimenum );
}


#ifdef _FLUID_PLASMA
void test_dH1_dUstar ( np_t * np, gl_t * gl, long l ) {
  sqmat_t C, Cnum;

  printf ( "\n" );
  printf ( "Testing dH/dUstar jacobian compared to numerically determined one.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  find_dH1_dUstar ( np, gl, l, C );
  find_numerical_jacobian_3 ( np, l, gl, &find_H1, Cnum );

  printf ( "dH/dUstar jacobian (analytical):\n" );
  display_matrix ( C );
  printf ( "dH/dUstar jacobian (numerical):\n" );
  display_matrix ( Cnum );
}


void test_dH2_dUstar ( np_t * np, gl_t * gl, long l ) {
  sqmat_t C, Cnum;

  printf ( "\n" );
  printf ( "Testing dH/dUstar jacobian compared to numerically determined one.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  find_dH2_dUstar ( np, gl, l, C );
  find_numerical_jacobian_3 ( np, l, gl, &find_H2, Cnum );

  printf ( "dH/dUstar jacobian (analytical):\n" );
  display_matrix ( C );
  printf ( "dH/dUstar jacobian (numerical):\n" );
  display_matrix ( Cnum );
}


void test_Dstar ( np_t * np, gl_t * gl, long l ) {
  sqmat_t Dstar, Dstarplus, Dstarminus, Dstarsum, Dstartest;
  flux_t Ustar, DstartesttimesUstar, DstartimesUstar;
  long theta, flux;

  printf ( "\n" );
  printf ( "Testing  Dstar plus and minus matrices.\n" );
  printf ( "Dstar should be equal to Dstarplus+Dstarminus.\n" );
  printf ( "The matrices are determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  theta = 0;
  find_Dstarplus ( np, gl, l, theta, Dstarplus );
  find_Dstarminus ( np, gl, l, theta, Dstarminus );
  find_Dstar_test ( np, gl, l, theta, Dstartest );
  find_Dstar ( np, gl, l, theta, Dstar );
  find_Ustar ( np[l], gl, Ustar );

  add_two_matrices ( Dstarminus, Dstarplus, Dstarsum );

  printf ( "Dstarplus:\n" );
  display_matrix ( Dstarplus );

  printf ( "Dstarminus:\n" );
  display_matrix ( Dstarminus );

  printf ( "Dstar:\n" );
  display_matrix ( Dstar );
  printf ( "Dstarplus+Dstarminus:\n" );
  display_matrix ( Dstarsum );

  multiply_matrix_and_vector ( Dstartest, Ustar, DstartesttimesUstar );
  multiply_matrix_and_vector ( Dstar, Ustar, DstartimesUstar );

  printf ( "\n\nBoth columns should be equal:\n" );
  for ( flux = 0; flux < nf; flux++ ) {
    printf ( "DstartesttimesUstar[%3.3ld]=%12.12E  DstartimesUstar[%3.3ld]=%12.12E \n", flux,
             DstartesttimesUstar[flux], flux, DstartimesUstar[flux] );
  }
}


#endif


void find_ZUstar(np_t *np, gl_t *gl, long l, flux_t ZUstar){
  sqmat_t Z;
  long flux;
  flux_t Ustar;
  find_Z(np, gl, l, Z);
  for (flux=0; flux<nf; flux++) Ustar[flux]=np[l].bs->U[flux]*_Omega(np[l],gl);
  multiply_matrix_and_vector(Z,Ustar,ZUstar);
}


void test_dZU_dU ( np_t * np, gl_t * gl, long l ) {
  sqmat_t C, Cnum;

  printf ( "\n" );
  printf
    ( "dZU/dU jacobian should be equal to numerically determined one when -DTEST compile flag is set.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  find_dZU_dU ( np, gl, l, C );
  find_numerical_jacobian_3 ( np, l, gl, &find_ZUstar, Cnum );

  printf ( "dZU/dU jacobian (analytical):\n" );
  display_matrix ( C );
  printf ( "dZU/dU jacobian (numerical):\n" );
  display_matrix ( Cnum );
}


void test_LambdaZ ( np_t * np, gl_t * gl, long l ) {
  sqmat_t LambdaZ,Z;
  flux_t Ustar, LambdaZtimesUstar, ZtimesUstar;
  long flux;

  printf ( "\n" );
  printf ( "LambdaZ is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );


  for (flux=0; flux<nf; flux++) Ustar[flux]=np[l].bs->U[flux]*_Omega(np[l],gl);
  find_Z(np, gl, l, Z);
  find_LambdaZ(np, gl, l, LambdaZ);


  printf ( "\n\nZ diagonal elements:\n" );
  for ( flux = 0; flux < nf; flux++ ) {
    printf ( "Z[%3.3ld]=%12.12E \n", flux,
             Z[flux][flux] );
  }

  printf ( "\n\nLambdaZ diagonal elements:\n" );
  for ( flux = 0; flux < nf; flux++ ) {
    printf ( "LambdaZ[%3.3ld]=%12.12E \n", flux,
             LambdaZ[flux][flux] );
  }


  multiply_matrix_and_vector ( Z, Ustar, ZtimesUstar );
  multiply_matrix_and_vector ( LambdaZ, Ustar, LambdaZtimesUstar );

  printf ( "\n\nBoth columns should be equal:\n" );
  for ( flux = 0; flux < nf; flux++ ) {
    printf ( "ZtimesUstar[%3.3ld]=%12.12E  LambdaZtimesUstar[%3.3ld]=%12.12E \n", flux,
             ZtimesUstar[flux], flux, LambdaZtimesUstar[flux] );
  }
}



void test_dS_dU ( np_t * np, gl_t * gl, long l ) {
  sqmat_t C, Cnum;

  printf ( "\n" );
  printf
    ( "dSstar/dUstar jacobian should be equal to numerically determined one when -DTEST compile flag is set.\n" );
  printf ( "The numerical jacobian is function of the Uref vector in the control file.\n" );
  printf ( "The jacobian is determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  find_dSstar_dUstar ( np, gl, l, C );
  find_numerical_jacobian_3 ( np, l, gl, &find_Sstar, Cnum );

  printf ( "dSstar/dUstar jacobian (analytical):\n" );
  display_matrix ( C );
  printf ( "dSstar/dUstar jacobian (numerical):\n" );
  display_matrix ( Cnum );
}

void test_dSchem_dU_local ( np_t * np, gl_t * gl, long l ) {

  printf ( "\n" );
  printf ( "Testing terms within dSchem/dUstar Jacobian.\n" );
  printf ( "The terms are determined at node " );
  printf_ijk_position_of_node_l ( gl, l );
  printf ( ".\n\n" );

  test_dSchem_dU ( np, gl, l );

}



void test_Lambda ( np_t * np, gl_t * gl, long lL, long lR, long theta ) {
  sqmat_t A, A2, lambdak;
  jacvars_t jacvars, jacvarsL, jacvarsR;
  long cnt, flux;
  double det, det2;
  metrics_t metrics;

  printf ( "\n" );
  printf
    ( "Testing eigenvalues by checking  if det(A-Lambda*I) is orders of magnitude less than det(A-2.2*Lambda*I).\n" );
  printf ( "The terms are determined at the interface between node " );
  printf_ijk_position_of_node_l ( gl, lL );
  printf ( " and node " );
  printf_ijk_position_of_node_l ( gl, lR );
  printf ( ".\n\n" );

  find_metrics_at_interface ( np, gl, lL, lR, theta, &metrics );
  find_jacvars ( np[lL], gl, metrics, theta, &jacvarsL );
  find_jacvars ( np[lR], gl, metrics, theta, &jacvarsR );
  find_jacvars_at_interface_from_jacvars ( jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING_ROE, &jacvars );

  for ( cnt = 0; cnt < nf; cnt++ ) {
    find_A_from_jacvars ( jacvars, metrics, A );
    find_A_from_jacvars ( jacvars, metrics, A2 );
    find_Lambda_from_jacvars ( jacvars, metrics, lambdak );
    for ( flux = 0; flux < nf; flux++ ) {
      A[flux][flux] = A[flux][flux] - lambdak[cnt][cnt];
    }
    for ( flux = 0; flux < nf; flux++ ) {
      A2[flux][flux] = A2[flux][flux] - 2.2e0 * lambdak[cnt][cnt];
    }
    find_determinant ( A, &det );
    find_determinant ( A2, &det2 );
    printf ( "det(A-Lambda[%ld]*I)=%+12.5E   det(A-2.2*Lambda[%ld]*I)=%+12.5E\n", cnt, det, cnt, det2 );
  }
  printf ( "\n\n" );
}

void test_Lambda_2 ( np_t * np, gl_t * gl, long lL, long lR, long theta ) {
  sqmat_t A, lambda;
  jacvars_t jacvars, jacvarsL, jacvarsR;
  long cnt, flux;
  double det, lambda_local, lambda_max, lambda_min, tmp;
  metrics_t metrics;

  printf ( "\n" );
  printf ( "Finding eigenvalues numerically by trying out det(A-lambda*I)\n" );
  printf ( "The terms are determined at the interface between node " );
  printf_ijk_position_of_node_l ( gl, lL );
  printf ( " and node " );
  printf_ijk_position_of_node_l ( gl, lR );
  printf ( ".\n\n" );

  find_metrics_at_interface ( np, gl, lL, lR, theta, &metrics );
  find_jacvars ( np[lL], gl, metrics, theta, &jacvarsL );
  find_jacvars ( np[lR], gl, metrics, theta, &jacvarsR );
  find_jacvars_at_interface_from_jacvars ( jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING_ROE, &jacvars );

  find_Lambda_from_jacvars ( jacvars, metrics, lambda );
  lambda_max = 0.0;
  for ( flux = 0; flux < nf; flux++ )
    lambda_max = max ( lambda[flux][flux], lambda_max );
  lambda_min = 9.9e99;
  for ( flux = 0; flux < nf; flux++ )
    lambda_min = min ( lambda[flux][flux], lambda_min );
  find_A_from_jacvars ( jacvars, metrics, A );
  printf ( "\nLambda_max=%+12.5E  Lambda_min=%+12.5E\n", lambda_max, lambda_min );

  tmp = lambda_max - lambda_min;
  lambda_min = lambda_min - tmp * 3.0;
  lambda_max = lambda_max + tmp * 3.0;

  for ( cnt = 0; cnt < 1000; cnt++ ) {
    lambda_local = ( lambda_max - lambda_min ) * ( double ) cnt / 1000.0e0 + lambda_min;
    for ( flux = 0; flux < nf; flux++ )
      A[flux][flux] = A[flux][flux] - lambda_local;
    find_determinant ( A, &det );
    for ( flux = 0; flux < nf; flux++ )
      A[flux][flux] = A[flux][flux] + lambda_local;
    printf ( "Lambda=%+12.5E det(A-Lambda*I)=%+12.5E\n", lambda_local, det );
  }
}

void test_LUstar ( np_t * np, gl_t * gl, long lL, long lR, long theta ) {
  sqmat_t L;
  flux_t Ustar, LtimesUstar, LUstar;
  metrics_t metrics;
  jacvars_t jacvarsL, jacvarsR, jacvars;
  long flux;

  printf ( "\n" );
  printf ( "L*Ustar should be equal to the LUstar vector.\n" );
  printf ( "The terms are determined at the interface between node " );
  printf_ijk_position_of_node_l ( gl, lL );
  printf ( " and node " );
  printf_ijk_position_of_node_l ( gl, lR );
  printf ( ".\n\n" );

  find_metrics_at_interface ( np, gl, lL, lR, theta, &metrics );
  find_jacvars ( np[lL], gl, metrics, theta, &jacvarsL );
  find_jacvars ( np[lR], gl, metrics, theta, &jacvarsR );
  find_jacvars_at_interface_from_jacvars ( jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING_ROE, &jacvars );

  find_L_from_jacvars ( jacvars, metrics, L );
  find_Ustar_from_jacvars ( jacvars, metrics, Ustar );
  multiply_matrix_and_vector ( L, Ustar, LtimesUstar );
  find_LUstar_from_jacvars ( jacvars, metrics, LUstar );

  for ( flux = 0; flux < nf; flux++ ) {
    printf ( "L_times_Ustar[%ld]=%+14.7E    LUstar[%ld]=%+14.7E\n", flux, LtimesUstar[flux], flux,
             LUstar[flux] );

  }
}



void test_Linv ( np_t * np, gl_t * gl, long lL, long lR, long theta ) {
  sqmat_t R, A, lambdak, Z1, Z2;
  metrics_t metrics;
  jacvars_t jacvarsL, jacvarsR, jacvars;

  printf ( "\n" );
  printf ( "Linv*Lambda should be equal to A*Linv.\n" );
  printf ( "The terms are determined at the interface between node " );
  printf_ijk_position_of_node_l ( gl, lL );
  printf ( " and node " );
  printf_ijk_position_of_node_l ( gl, lR );
  printf ( ".\n\n" );

  find_metrics_at_interface ( np, gl, lL, lR, theta, &metrics );
  find_jacvars ( np[lL], gl, metrics, theta, &jacvarsL );
  find_jacvars ( np[lR], gl, metrics, theta, &jacvarsR );
  find_jacvars_at_interface_from_jacvars ( jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING_ROE, &jacvars );

  find_Linv_from_jacvars ( jacvars, metrics, R );
  find_A_from_jacvars ( jacvars, metrics, A );
  find_Lambda_from_jacvars ( jacvars, metrics, lambdak );
  multiply_matrix_and_matrix ( R, lambdak, Z1 );
  multiply_matrix_and_matrix ( A, R, Z2 );

  printf ( "Linv * Lambda \n" );
  display_matrix ( Z1 );
  printf ( "A * Linv \n" );
  display_matrix ( Z2 );
}


void test_L ( np_t * np, gl_t * gl, long lL, long lR, long theta ) {
  sqmat_t Linv, L, Mult;
// sqmar_t Lnum;
  jacvars_t jacvarsL, jacvarsR, jacvars;
  metrics_t metrics;

  printf ( "\n" );
  printf ( "Linv*L should give the identity matrix.\n" );
  printf ( "The terms are determined at the interface between node " );
  printf_ijk_position_of_node_l ( gl, lL );
  printf ( " and node " );
  printf_ijk_position_of_node_l ( gl, lR );
  printf ( ".\n\n" );

  find_metrics_at_interface ( np, gl, lL, lR, theta, &metrics );
  find_jacvars ( np[lL], gl, metrics, theta, &jacvarsL );
  find_jacvars ( np[lR], gl, metrics, theta, &jacvarsR );
  find_jacvars_at_interface_from_jacvars ( jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING_ROE, &jacvars );

  find_Linv_from_jacvars ( jacvars, metrics, Linv );
  find_L_from_jacvars ( jacvars, metrics, L );
  multiply_matrix_and_matrix ( Linv, L, Mult );

  printf ( "Linv*L Analytical \n" );
  display_matrix ( Mult );

/*
   invert_matrix(Linv,Lnum);
   printf("L Analytical (found from the L subroutine)  \n");
   display_matrix(L);
   printf("L Numerical (found from the Linv subroutine - inverted) \n");
   display_matrix(Lnum);
*/
}


void test_LinvLambdaL ( np_t * np, gl_t * gl, long lL, long lR, long theta ) {
  sqmat_t R, L, lambdak, Z1, Z2, A;
  jacvars_t jacvarsL, jacvarsR, jacvars;
  metrics_t metrics;

  printf ( "\n" );
  printf ( "Linv*Lambda*L should be equal to A jacobian.\n" );
  printf ( "The terms are determined at the interface between node " );
  printf_ijk_position_of_node_l ( gl, lL );
  printf ( " and node " );
  printf_ijk_position_of_node_l ( gl, lR );
  printf ( ".\n\n" );

  find_metrics_at_interface ( np, gl, lL, lR, theta, &metrics );
  find_jacvars ( np[lL], gl, metrics, theta, &jacvarsL );
  find_jacvars ( np[lR], gl, metrics, theta, &jacvarsR );
  find_jacvars_at_interface_from_jacvars ( jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING_ROE, &jacvars );

  find_Linv_from_jacvars ( jacvars, metrics, R );
  find_L_from_jacvars ( jacvars, metrics, L );
  find_Lambda_from_jacvars ( jacvars, metrics, lambdak );
  multiply_matrix_and_matrix ( lambdak, L, Z1 );
  multiply_matrix_and_matrix ( R, Z1, Z2 );
  find_A_from_jacvars ( jacvars, metrics, A );

  printf ( "Linv * Lambda * L:\n" );
  display_matrix ( Z2 );
  printf ( "A jacobian:\n" );
  display_matrix ( A );
}

void test_Roe_average ( np_t * np, gl_t * gl, long lL, long lR, long theta ) {
  sqmat_t A;
  flux_t dUstar, FL, FR, F, F2;
  long flux;
  metrics_t metrics;
  jacvars_t jacvarsL, jacvarsR, jacvars;

  printf ( "\n" );
  printf ( "Roe Averaging implied that Delta Fstar=A*(Delta Ustar).\n" );
  printf ( "The terms are determined at the interface between node " );
  printf_ijk_position_of_node_l ( gl, lL );
  printf ( " and node " );
  printf_ijk_position_of_node_l ( gl, lR );
  printf ( ".\n" );
  printf ( "NOTE: flow conditions must not be constant between the two nodes.\n" );
  printf ( "\n" );

  find_metrics_at_interface ( np, gl, lL, lR, theta, &metrics );
  find_jacvars ( np[lL], gl, metrics, theta, &jacvarsL );
  find_jacvars ( np[lR], gl, metrics, theta, &jacvarsR );
  find_jacvars_at_interface_from_jacvars ( jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING_ROE, &jacvars );

  find_A_from_jacvars ( jacvars, metrics, A );

  find_Fstar ( np[lL], gl, theta, FL );
  find_Fstar ( np[lR], gl, theta, FR );
  for ( flux = 0; flux < nf; flux++ ) {
    F[flux] = FR[flux] - FL[flux];
    dUstar[flux] = np[lR].bs->U[flux] * _Omega ( np[lR], gl ) - np[lL].bs->U[flux] * _Omega ( np[lL], gl );
    printf ( "dUstar[%ld]=%+12.5E\n", flux, dUstar[flux] );
  }
  multiply_matrix_and_vector ( A, dUstar, F2 );

  printf ( "\n" );
  printf ( "Delta Fstar     A*(Delta Ustar)\n" );
  for ( flux = 0; flux < nf; flux++ ) {
    printf ( "%+14.7E  %+14.7E \n", F[flux], F2[flux] );
  }
}


void test_jacvars ( np_t * np, gl_t * gl, long l, long theta ) {
  metrics_t metrics;
  jacvars_t jacvars_from_U, jacvars, jacvars_from_musclvars;
  sqmat_t A;
  flux_t musclvars;

  find_metrics_at_node ( np, gl, l, theta, &metrics );
  find_musclvars ( np[l], gl, musclvars );
  find_jacvars_from_U ( np[l].bs->U, metrics, gl, theta, &jacvars_from_U );
  find_jacvars ( np[l], gl, metrics, theta, &jacvars );
  find_jacvars_from_musclvars ( musclvars, metrics, gl, theta, &jacvars_from_musclvars );

  find_A_from_jacvars ( jacvars, metrics, A );
  printf ( "A Jacobian from jacvars\n" );
  display_matrix ( A );

  find_A_from_jacvars ( jacvars_from_musclvars, metrics, A );
  printf ( "A Jacobian from jacvars from musclvars\n" );
  display_matrix ( A );

  find_A_from_jacvars ( jacvars_from_U, metrics, A );
  printf ( "A Jacobian from jacvars from U\n" );
  display_matrix ( A );

}

void print_column_species_names ( long digits ) {
  long spec, cnt;
  char *speciesname = ( char * ) malloc ( sizeof ( char ) );

  printf ( " T" );
  for ( cnt = 0; cnt < digits; cnt++ )
    printf ( " " );
  for ( spec = 0; spec < ns; spec++ ) {
    find_species_name ( spec, &speciesname );
    printf ( "%s", speciesname );
    for ( cnt = 0; cnt < 1 + digits - strlen ( speciesname ); cnt++ )
      printf ( " " );
  }
  printf ( "\n" );
  free ( speciesname );

}

void test_h ( double Tmin, double Tmax, double dT ) {
  double T, h;
  spec_t w;
  long spec, spec2;
  printf ( "\n" );
  printf ( "Species specific enthalpies [J/kg] as function of temperature [K].\n" );
  printf ( "Tmin=%EK Tmax=%EK dT=%EK.\n", Tmin, Tmax, dT );
  printf ( "\n" );
  print_column_species_names ( 12 );
  T = Tmin;
  do {
    wfprintf ( stdout, "%12.5E ", T );
    for ( spec = 0; spec < ns; spec++ ) {
      for ( spec2 = 0; spec2 < ns; spec2++ )
        w[spec2] = 0.0;
      w[spec] = 1.0;
      h = _h_from_w_T ( w, T );
      wfprintf ( stdout, "%+12.5E ", h );
    }
    wfprintf ( stdout, "\n" );
    T += dT;
  } while ( T < Tmax );

}

void test_cp ( double Tmin, double Tmax, double dT ) {
  double T, Cp;
  spec_t w;
  long spec, spec2;
  printf ( "\n" );
  printf ( "Species specific heats [J/kgK] as function of temperature [K].\n" );
  printf ( "Tmin=%EK Tmax=%EK dT=%EK.\n", Tmin, Tmax, dT );
  printf ( "\n" );
  print_column_species_names ( 12 );
  T = Tmin;
  do {
    wfprintf ( stdout, "%12.5E ", T );
    for ( spec = 0; spec < ns; spec++ ) {
      for ( spec2 = 0; spec2 < ns; spec2++ )
        w[spec2] = 0.0;
      w[spec] = 1.0;
      Cp = _cp_from_w_T ( w, T );
      wfprintf ( stdout, "%+12.5E ", Cp );
    }
    wfprintf ( stdout, "\n" );
    T += dT;
  } while ( T < Tmax );

}

void test_s ( double Tmin, double Tmax, double dT ) {
  double T, s;
  spec_t w;
  long spec, spec2;
  printf ( "\n" );
  printf ( "Species specific entropies [J/kgK] as function of temperature [K].\n" );
  printf ( "Tmin=%EK Tmax=%EK dT=%EK.\n", Tmin, Tmax, dT );
  printf ( "\n" );
  print_column_species_names ( 12 );
  T = Tmin;
  do {
    wfprintf ( stdout, "%12.5E ", T );
    for ( spec = 0; spec < ns; spec++ ) {
      for ( spec2 = 0; spec2 < ns; spec2++ )
        w[spec2] = 0.0;
      w[spec] = 1.0;
      s = _s_from_w_T ( w, T );
      wfprintf ( stdout, "%+12.5E ", s );
    }
    wfprintf ( stdout, "\n" );
    T += dT;
  } while ( T < Tmax );

}

void test_dsdT ( double Tmin, double Tmax, double dT ) {
  double T, dsdT;
  long spec;
  printf ( "\n" );
  printf ( "Species dsdT [J/(kg K^2)] as function of temperature [K]. \n" );
  printf ( "Tmin=%EK Tmax=%EK dT=%EK.\n", Tmin, Tmax, dT );
  printf ( "\n" );
  print_column_species_names ( 12 );
  T = Tmin;
  do {
    wfprintf ( stdout, "%12.5E ", T );
    for ( spec = 0; spec < ns; spec++ ) {
      dsdT = _dsk_dT_from_T ( spec, T );

      wfprintf ( stdout, "%+12.5E ", dsdT );
    }
    wfprintf ( stdout, "\n" );
    T += dT;
  } while ( T < Tmax );

}

int chkarg ( int argc, char **argv, char *arg ) {
  int cnt, tmp;
  tmp = 0;
  for ( cnt = 1; cnt < argc; cnt++ ) {
    if ( strcmp ( argv[cnt], arg ) == 0 ) {
      tmp = cnt;
    }
  }
  return ( tmp );
}

int main ( int argc, char **argv ) {
  np_t npL, npR;
  long theta, cnt, lL, lR, node_i, node_j, node_k;
  char *control_filename;
  np_t *npArray;
  gl_t gl;
  bool CONTROL;
  int RET;
  input_t input;
  int *argint, argoneint;
  char tmpstr[1000];

#ifdef _FLUID_MULTISPECIES
  double Tmin, Tmax, dT;
#endif
  int term_width, term_height, linewidth;

  find_terminal_window_size ( &term_width, &term_height );
  linewidth = min ( MAX_LINE_WIDTH, max ( 43, term_width - 2 ) );

  argint = NULL;
  input.name = NULL;
  input.READDATAFILE = FALSE;
#ifdef UNSTEADY
  input.M1 = FALSE;
  input.name_m1 = NULL;
#if _RESTIME_BW > 2
  input.M2 = FALSE;
  input.name_m2 = NULL;
#endif
#if _RESTIME_BW > 3
  input.M3 = FALSE;
  input.name_m3 = NULL;
#endif
#endif //UNSTEADY
  input.ASCII = FALSE;
  input.INTERPOLATION = FALSE;

  seed_random (  );
  create_node ( &npL, 0, 0, 0 );
  create_node ( &npR, 0, 0, 0 );

  CONTROL = FALSE;
  control_filename = NULL;
  input.READDATAFILE = FALSE;
#ifdef UNSTEADY
  input.M1 = FALSE;
#endif
  if ( process_flag_string ( argc, argv, "-r", &control_filename ) )
    CONTROL = TRUE;

  if ( CONTROL && argc >= 5 + nd ) {
    if ( process_flag_string ( argc, argv, "-i", &input.name ) )
      input.READDATAFILE = TRUE;

    read_control ( control_filename, input, TRUE, FALSE, FALSE, FALSE, &npArray, &gl );
    resume_nodes_only_in_zone_and_update_bdry_nodes ( npArray, &gl, gl.domain );
    update_bdry_nodes ( npArray, &gl, gl.domain );


    theta = 0;
    RET = process_flag_int ( argc, argv, "-dim", &argoneint );
    if ( RET ) {
      if ( RET != 2 )
        fatal_error ( "After -dim flag, you must supply 1 integer argument." );
      if ( argoneint<0 || argoneint>=nd ) 
        fatal_error ( "After -dim flag, the dimension must be 0"if2DL(", 1")if3DL(", 2")"\n");
      theta=(long)argoneint;
    }

    node_i = 1;
    node_j = 1;
    node_k = 1;
    RET = process_flag_int_multiple ( argc, argv, "-node", &argint );

    if ( RET ) {
      if ( RET != nd + 1 )
        fatal_error ( "After -node flag, you must supply %d integer arguments.", nd );
      node_i = argint[0];
#ifdef _2DL
      node_j = argint[1];
#endif
#ifdef _3DL
      node_k = argint[2];
#endif
      printf ( "Evaluating properties at node (%ld,%ld,%ld).\n", node_i, node_j, node_k );
      lL = _ai ( &gl, node_i, node_j, node_k );
      if ( !is_node_in_zone_2 ( lL, &gl, gl.domain ) )
        fatal_error ( "Node specified after -node flag is not within domain." );
      if ( !is_node_in_zone_2 ( _al(&gl, lL, theta, +1), &gl, gl.domain ) )
        fatal_error ( "Node specified after -node flag leads to node on its right not being within domain." );
    } else {
      fatal_error ( "You must specify -node flag." );
    }



    lL = _ai ( &gl, node_i, node_j, node_k );
    lR = _al(&gl, lL, theta, +1);
    dispose_node ( &npL );
    dispose_node ( &npR );
    npL = npArray[lL];
    npR = npArray[lR];

    for ( cnt = 3; cnt < argc; cnt++ ) {
      if ( strcmp ( "A", argv[cnt] ) == 0 )
        test_A ( npArray, lL, &gl, theta );
      if ( strcmp ( "Ajacvars", argv[cnt] ) == 0 )
        test_A_from_jacvars ( npArray, lL, &gl, theta );
      if ( strcmp ( "dGdUstar", argv[cnt] ) == 0 )
        test_dG_dUstar ( npArray, lL, &gl );
#ifdef _FLUID_PLASMA
      if ( strcmp ( "dH1dUstar", argv[cnt] ) == 0 )
        test_dH1_dUstar ( npArray, &gl, _ai ( &gl, node_i, node_j, node_k ) );
      if ( strcmp ( "dH2dUstar", argv[cnt] ) == 0 )
        test_dH2_dUstar ( npArray, &gl, _ai ( &gl, node_i, node_j, node_k ) );
      if ( strcmp ( "Dstar", argv[cnt] ) == 0 )
        test_Dstar ( npArray, &gl, _ai ( &gl, node_i, node_j, node_k ) );
#endif
      if ( strcmp ( "dZdU", argv[cnt] ) == 0 )
        test_dZU_dU ( npArray, &gl, _ai ( &gl, node_i, node_j, node_k ) );
      if ( strcmp ( "dSdU", argv[cnt] ) == 0 )
        test_dS_dU ( npArray, &gl, _ai ( &gl, node_i, node_j, node_k ) );
      if ( strcmp ( "dSchemdU", argv[cnt] ) == 0 )
        test_dSchem_dU_local ( npArray, &gl, _ai ( &gl, node_i, node_j, node_k ) );
      if ( strcmp ( "Roe", argv[cnt] ) == 0 )
        test_Roe_average ( npArray, &gl, lL, lR, theta );
      if ( strcmp ( "Lambda", argv[cnt] ) == 0 ) {
        test_Lambda ( npArray, &gl, lL, lR, theta );
        test_LinvLambdaL ( npArray, &gl, lL, lR, theta );
      }
      if ( strcmp ( "Lambda2", argv[cnt] ) == 0 )
        test_Lambda_2 ( npArray, &gl, lL, lR, theta );
      if ( strcmp ( "LambdaZ", argv[cnt] ) == 0 )
        test_LambdaZ ( npArray, &gl, _ai ( &gl, node_i, node_j, node_k ) );
      if ( strcmp ( "Linv", argv[cnt] ) == 0 )
        test_Linv ( npArray, &gl, lL, lR, theta );
      if ( strcmp ( "L", argv[cnt] ) == 0 )
        test_L ( npArray, &gl, lL, lR, theta );
      if ( strcmp ( "LUstar", argv[cnt] ) == 0 )
        test_LUstar ( npArray, &gl, lL, lR, theta );
      if ( strcmp ( "dUstardUprime", argv[cnt] ) == 0 )
        test_dUstar_dUprime ( npArray, lL, &gl );
      if ( strcmp ( "jacvars", argv[cnt] ) == 0 )
        test_jacvars ( npArray, &gl, lL, theta );
    }

  } else {
#ifdef _FLUID_MULTISPECIES
    if ( argc == 5 ) {
      sscanf ( argv[2], "%lg", &Tmin );
      sscanf ( argv[3], "%lg", &Tmax );
      sscanf ( argv[4], "%lg", &dT );
      if ( strcmp ( "h", argv[1] ) == 0 )
        test_h ( Tmin, Tmax, dT );
      if ( strcmp ( "cp", argv[1] ) == 0 )
        test_cp ( Tmin, Tmax, dT );
      if ( strcmp ( "s", argv[1] ) == 0 )
        test_s ( Tmin, Tmax, dT );
      if ( strcmp ( "dsdT", argv[1] ) == 0 )
        test_dsdT ( Tmin, Tmax, dT );
    } else {
#endif
      write_hline ( stderr, linewidth, 2 );
      write_options_row ( stderr, "Flag", "Argument(s)", "Description", linewidth, lengthcol1,
                          lengthcol2 );
      write_hline ( stderr, linewidth, 2 );

      write_options_row ( stderr, "-r", "string", "control file name", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "-i", "string", "Input binary data file", linewidth, lengthcol1,
                          lengthcol2 );
      sprintf(tmpstr,"%d int",nd); 
      write_options_row ( stderr, "-node", tmpstr, "i,j"if3D(",k")" indices", linewidth, lengthcol1, lengthcol2 );
      write_options_row ( stderr, "-dim", "1 int", "dimension", linewidth, lengthcol1, lengthcol2 );
      write_options_row ( stderr, "A", "none", 
                          "dFstar/dUstar Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" A", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "Ajacvars", "none", 
                          "dFstar/dUstar Jacobian from jacvars  ./test -r test.wrp -node 10 10"if3D(" 10")" Ajacvars",
                          linewidth, lengthcol1, lengthcol2 );
      write_options_row ( stderr, "dGdUstar", "none", 
                          "dG/dUstar Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" dGdUstar", linewidth,
                          lengthcol1, lengthcol2 );
#ifdef _FLUID_PLASMA
      write_options_row ( stderr, "dH1dUstar", "none", 
                          "dH1/dUstar Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" dH1dUstar", linewidth,
                          lengthcol1, lengthcol2 );
      write_options_row ( stderr, "dH2dUstar", "none", 
                          "dH2/dUstar Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" dH2dUstar", linewidth,
                          lengthcol1, lengthcol2 );
      write_options_row ( stderr, "Dstar", "none", 
                          "Dstar plus/minus matrices  ./test -r test.wrp -node 10 10"if3D(" 10")" Dstar", linewidth,
                          lengthcol1, lengthcol2 );
#endif
      write_options_row ( stderr, "dZdU", "none", 
                          "dZ/dU Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" dZdU", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "dSdU", "none", 
                          "dS/dU Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" dSdU", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "dSchemdU", "none", 
                          "dSchem/dU Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" dSchemdU", linewidth,
                          lengthcol1, lengthcol2 );
      write_options_row ( stderr, "Linv", "none", 
                          "Right eigenvectors  ./test -r test.wrp -node 10 10"if3D(" 10")" Linv", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "L", "none", 
                          "Left eigenvectors  ./test -r test.wrp -node 10 10"if3D(" 10")" L ", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "Lambda", "none", 
                          "eigenvalues  ./test -r test.wrp -node 10 10"if3D(" 10")" Lambda ", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "Lambda2", "none", 
                          "test eigenvalues  ./test -r test.wrp -node 10 10"if3D(" 10")" Lambda2 ", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "LambdaZ", "none", 
                          "test Z matrix eigenvalues  ./test -r test.wrp -node 10 10"if3D(" 10")" LambdaZ", linewidth,
                          lengthcol1, lengthcol2 );
      write_options_row ( stderr, "LUstar", "none", 
                          "characteristic variables  ./test -r test.wrp -node 10 10"if3D(" 10")" LUstar", linewidth,
                          lengthcol1, lengthcol2 );
      write_options_row ( stderr, "dUstardUprime", "none", 
                          "dUstar/dUprime Jacobian  ./test -r test.wrp -node 10 10"if3D(" 10")" dUstardUprime", linewidth,
                          lengthcol1, lengthcol2 );
      write_options_row ( stderr, "Roe", "none", 
                          "Roe Averaging  ./test -r test.wrp -node 10 10"if3D(" 10")" Roe", linewidth, lengthcol1,
                          lengthcol2 );
      write_options_row ( stderr, "jacvars", "none", 
                          "jacvars from U  ./test -r test.wrp -node 10 10"if3D(" 10")" jacvars", linewidth, lengthcol1,
                          lengthcol2 );
#ifdef _FLUID_MULTISPECIES
      write_options_row ( stderr, "h", "none", "species specific enthalpy  ./test h 500 2000 2",
                          linewidth, lengthcol1, lengthcol2 );
      write_options_row ( stderr, "cp", "none", "species specific heat  ./test cp 500 2000 2",
                          linewidth, lengthcol1, lengthcol2 );
      write_options_row ( stderr, "s", "none", "species specific entropy  ./test s 500 2000 2",
                          linewidth, lengthcol1, lengthcol2 );
      write_options_row ( stderr, "dsdT", "none", 
                          "species dsdT derivative  ./test dsdT 500 2000 2", linewidth, lengthcol1,
                          lengthcol2 );
#endif
      write_hline ( stderr, linewidth, 2 );


      exit ( EXIT_FAILURE );
#ifdef _FLUID_MULTISPECIES
    }
#endif
  }
  //dispose_node ( &npL );
  //dispose_node ( &npR );
/* 
  long i,j,k;
  for_ijk ( gl.domain_lim,is,js,ks,ie,je,ke ){
        dispose_node ( &( npArray[_ai ( &gl, i, j, k )] ) );
  }
  free_clipped_variables( &gl );
  free( gl.cycle.code_runtime );
  SOAP_free_codex ( &gl.cycle.codex );
  free ( npArray ); 
  free ( control_filename );
  free ( input.name );
  free ( argint );*/
  return ( EXIT_SUCCESS );
}
