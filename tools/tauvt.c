/* Vibrational-Translational (V-T) relaxation time of N2 for N2-O2 mixture's:
   Millikan and White gave an empirical correlation for obtaining the
   V-T relaxation based on molecular collision's[1]. It can be applied to 
   mixture by considering the individual collisions between the species. And 
   Macheret expressed N2 vibrational relaxation time in terms of total
   number density and translational temperature (T) [2]. This program,
   gives the (V-T) relaxation time of N2 using both Millikan - White correlation
   and Macheret correlation. And from the experiments(compression & expansion) 
   of Sebacher, correlation is derived for N2 relaxation in air. They are also implemented
   in this code.     
[1]. MILLIKAN, R. AND WHITE, D., “Systematics of Vibrational Relaxation,”The Journal of Chemical Physics, Vol. 39,
No. 12, 1963, pp. 3209–3213.
[2]. PARENT, B., MACHERET, S. O., SHNEIDER, M. N., AND HARADA, N., “Numerical Study of an Electron-
Beam-confined Faraday accelerator,” Journal of Propulsion and Power, Vol. 23, No. 5, 2007, pp. 1023–1031.
[3]. ECKSTROM, D.J., “Vibrational relaxation of shock-heated N2 by atomic oxygen using the ir tracer
method,” The Journal of Chemical Physics, Vol.59, No. 6, 1973, pp. 2787-2795.  
[4]. SEBACHER, D.I., GUY, R. W., “Vibrational relaxation in expanding N2 and Air,” NASA TM X-71988, August 1974.  
*/
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#define Theta_N2 3395.0e0       /*characteristic vibrational temperature of N2 in Kelvin */
#define Theta_O2 2239.0e0       /*characteristic vibrational temperature of O2 in Kelvin */
#define M_O2     32.0e0         /* molecular weight of O2 in g/mol */
#define M_O      16.00e0        /* molecular weight of O in g/mol */
#define M_N2     28.014e0       /* molecular weight of N2 in g/mol */
#define P_conv   101300.0e0     /* converting the pressure in pa to atm */
#define R_N2     296.8e0        /* individual gas constant of N2 in J/Kg.K */
#define R_O2     259.812e0      /* individual gas constant of O2 in J/Kg.K */
#define R_O      519.625e0      /* individual gas constant of O in J/Kg.K  */
#define Ava_num  6.0221413e23   /* Avagadro's constant 1/mol */
#define min(x, y) (((x) < (y)) ? (x) : (y))     /* minimum of the two numbers  */

/*  "check_arg" function: This function will check 
     whether the flag's which are needed for calculation 
     exists or not   
   > argc : Number of arguments passed into the program from
     the command line
   > argv: Array of arguments passed into the program from
     the command line 
   > arg: Stores the flag argument like "-P","-T","-w_O2"   
   > This function compares the argument in "argv" and "arg"   
     which is passed into the function. If argument matches, it will 
     return the location of argument or else it will return 0.  
*/
int check_arg ( int argc, char **argv, char *arg ) {
  int i, temp;
  temp = 0;
  for ( i = 1; i < argc; i++ ) {
    if ( strcmp ( argv[i], arg ) == 0 ) {
      temp = i;
    }
  }
  return temp;
}

/* tau_millikan_2species function: This function calculates the 
   relaxation rate between two species using Millikan correlation
   > M1: Molecular weight of relaxing species with vibrational energy
   > M2: Molecular weight of the colliding species 
   > theta: Characteristic vibrational temperature of the relaxing species with vibrational energy
   > T: Translational temperature in Kelvin
   > P: Pressure in atm  
   > For N2-O2 mixture we know the collision probabilities, molecular weight's of the relaxing and colliding species are passed 
     into the function for corresponding translational temperature and pressure. This function will return
     the V-T relaxation rate for the considered collision.
*/
double tau_millikan_2species ( double M1, double M2, double theta, double T, double P ) {
  double mu, A, B;
  double tau;
          /************************************************************/
  /* the following formula's are taken from this paper [1]
     > see p.3211 for tau relation. 
     > Relaxation coefficients A and B is given in p.3213 eqn(3.a)
   */
         /************************************************************/
  mu = ( ( M1 * M2 ) / ( M1 + M2 ) );   /* calculating the equivalent molecular weight between collisions */
  A = 0.00116 * pow ( mu, ( 1.0 / 2.0 ) ) * pow ( theta, ( 4.0 / 3.0 ) );       /* Millikan's relaxation coefficient A */
  B = 0.015 * pow ( mu, ( 1.0 / 4.0 ) );        /*Millikan's relaxation coefficient B */
  tau = ( exp ( ( A * ( pow ( T, -( 1.0 / 3.0 ) ) - B ) - 18.42 ) ) ) / P;      /*Relaxation time calculation, here P is specified in atm */
  return tau;
}

/* tau_millikan function: This function calculates the total
   V-T relaxation rate of the mixture
   > w_O2: Mass fraction of oxygen
   > w_O:  Mass fraction of atomic oxygen
   > T: Translational temperature in Kelvin
   > P: Pressure in Pa  
   > Individual relaxation rate calculated for different collisions in tau_millikan_2species function
     is used in this function to determine the N2 V-T relaxation rate of the mixture. 
*/
double tau_millikan ( double w_O2, double w_O, double P, double T ) {
  double w_N2, chi_O2, chi_O, chi_N2, P_atm;
  double tau_N2_N2, tau_N2_O2, tau_N2_O, tau;
  w_N2 = 1.0 - ( w_O2 + w_O );
  P_atm = P / P_conv;
  chi_O2 = ( w_O2 / M_O2 ) / ( ( w_O2 / M_O2 ) + ( w_N2 / M_N2 ) + ( w_O / M_O ) );
  chi_O = ( w_O / M_O ) / ( ( w_O2 / M_O2 ) + ( w_N2 / M_N2 ) + ( w_O / M_O ) );
  chi_N2 = 1.0 - ( chi_O2 + chi_O );
  tau_N2_N2 = tau_millikan_2species ( M_N2, M_N2, Theta_N2, T, P_atm );
  tau_N2_O2 = tau_millikan_2species ( M_N2, M_O2, Theta_N2, T, P_atm );
  tau_N2_O = tau_millikan_2species ( M_N2, M_O, Theta_N2, T, P_atm );
 /**********************************************************/
  /*the following formula is taken from this paper[3]        
     see p.2791 eqn.23, T_A is the V-T relaxation parameter of N2, and the general form of this equation is 
     expressed in p.2790 eqn. 16 */
/***********************************************************/
  tau = 1.0 / ( ( chi_N2 / tau_N2_N2 ) + ( chi_O2 / tau_N2_O2 ) + ( chi_O / tau_N2_O ) );
  return tau;
}

/* tau_macheret function: This function calculates the total
   V-T relaxation rate of the mixture using Macheret correlation
   > w_O2: Mass fraction of oxygen
   > w_O:  Mass fraction of atomic oxygen
   > P: Pressure in Pa  
   > T: Translational temperature in Kelvin
   > Total relaxation rate of the mixture is determined from total number density and translational temperature.
     Macheret's relaxation rate includes the correction term which is the coupling between vibration and 
     dissociation(inclusion of oxygen atom mass fraction in the equation). 
*/
double tau_macheret ( double w_O2, double w_O, double P, double T ) {
  double w_N2, N_N2, N_O2, N_O, N, tau;
  w_N2 = 1.0 - ( w_O2 + w_O );
  N_N2 = ( ( ( P * w_N2 ) / ( R_N2 * T ) ) * Ava_num ) / ( M_N2 / 1000.0 );     /* Number density of N2 in 1/m^3 */
  N_O2 = ( ( ( P * w_O2 ) / ( R_O2 * T ) ) * Ava_num ) / ( M_O2 / 1000.0 );     /* Number density of O2 in 1/m^3  */
  N_O = ( ( ( P * w_O ) / ( R_O * T ) ) * Ava_num ) / ( M_O / 1000.0 ); /* Number density of O in 1/m^3 */
  N = N_N2 + N_O2 + N_O;        /* total number density of the mixture in 1/m^3  */
         /************************************************************/
  /* the following formula is taken from this paper [2]
     see p.1025 for tau relation. */
         /************************************************************/
  tau =
    1.0 / ( N *
            ( ( 7.0 * pow ( 10.0, -16.0 ) * exp ( -141.0 * pow ( T, -( 1.0 / 3.0 ) ) ) ) +
              ( w_O * 5.0 * pow ( 10.0, -18.0 ) * exp ( -128.0 * pow ( T, -( 1.0 / 2.0 ) ) ) ) ) );

  return tau;
}

/* tau_sebacher_compression function: This function calculates N2 V-T relaxation
   in air through a compression process. 
   > w_O2: Mass fraction of oxygen
   > w_O:  Mass fraction of atomic oxygen
   > P: Pressure in Pa  
   > T: Translational temperature in Kelvin
   > For compression experiments(Ex.shock tube) N2 vibrational relaxation rate in air is determined and relaxation rate variation with temperature is  shown by Sebacher et al [4]. From that data correlation is derived and implemented here.        
*/

double tau_sebacher_compression ( double w_O2, double w_O, double P, double T ) {
  double w_N2, N_N2, N_O2, N_O, N, tau;
  w_N2 = 1.0 - ( w_O2 + w_O );
  N_N2 = ( ( ( P * w_N2 ) / ( R_N2 * T ) ) * Ava_num ) / ( M_N2 / 1000.0 );     /* Number density of N2 in 1/m^3 */
  N_O2 = ( ( ( P * w_O2 ) / ( R_O2 * T ) ) * Ava_num ) / ( M_O2 / 1000.0 );     /* Number density of O2 in 1/m^3  */
  N_O = ( ( ( P * w_O ) / ( R_O * T ) ) * Ava_num ) / ( M_O / 1000.0 ); /* Number density of O in 1/m^3 */
  N = N_N2 + N_O2 + N_O;        /* total number density of the mixture in 1/m^3  */
        /************************************************************/
  /* the following correlation is fitted using Levenberg–Marquardt least square fitting method 
     for the data taken from this report [4] see Figure-4 (p.18) for N2 V-T relaxation in air 
     data for compression experiments */
         /************************************************************/
  tau = 1.0 / ( N * 1.0109 * pow ( 10.0, -16.0 ) * exp ( -135.581 * pow ( T, -( 1.0 / 3.0 ) ) ) );
  return tau;
}

/* tau_sebacher_expansion function: This function calculates N2 V-T relaxation
   in air through a expansion process.  
   > w_O2: Mass fraction of oxygen
   > w_O:  Mass fraction of atomic oxygen
   > P: Pressure in Pa  
   > T: Translational temperature in Kelvin
> For expansion experiments(Ex.Nozzle studies) N2 vibrational relaxation rate in air is determined and relaxation rate variation with temperature is shown by Sebacher et al [4]. From that data correlation is derived and implemented here.  
*/

double tau_sebacher_expansion ( double w_O2, double w_O, double P, double T ) {
  double w_N2, N_N2, N_O2, N_O, N, tau;
  w_N2 = 1.0 - ( w_O2 + w_O );
  N_N2 = ( ( ( P * w_N2 ) / ( R_N2 * T ) ) * Ava_num ) / ( M_N2 / 1000.0 );     /* Number density of N2 in 1/m^3 */
  N_O2 = ( ( ( P * w_O2 ) / ( R_O2 * T ) ) * Ava_num ) / ( M_O2 / 1000.0 );     /* Number density of O2 in 1/m^3  */
  N_O = ( ( ( P * w_O ) / ( R_O * T ) ) * Ava_num ) / ( M_O / 1000.0 ); /* Number density of O in 1/m^3 */
  N = N_N2 + N_O2 + N_O;        /* total number density of the mixture in 1/m^3  */
        /************************************************************/
  /* the following correlation is fitted using Levenberg–Marquardt least square fitting method 
     for the data taken from this report [4] see Figure-4 (p.18) for N2 V-T relaxation in air 
     data for expansion experiments */
         /************************************************************/
  tau = 1.0 / ( N * 3.2834 * pow ( 10.0, -14.0 ) * exp ( -160.1098 * pow ( T, -( 1.0 / 3.0 ) ) ) );
  return tau;
}

/* main program*/
int main ( int argc, char **argv ) {
  double P, T, w_O2, w_O;
  double P1, T1, w_O2_1, w_O_1;
  double P2, T2, w_O2_2, w_O_2;
  double tau_vt_macheret, tau_vt_millikan, tau_sebacher_comp, tau_sebacher_expan;
  int steps, cnt;
  bool check = true;
  if ( check_arg ( argc, argv, "-P" ) == 0 )
    check = false;
  else {
    sscanf ( argv[check_arg ( argc, argv, "-P" ) + 1], "%lg", &P1 );
    if ( sscanf ( argv[min ( check_arg ( argc, argv, "-P" ) + 2, argc - 1 )], "%lg", &P2 ) != 1 )
      P2 = P1;
  }
  if ( check_arg ( argc, argv, "-T" ) == 0 )
    check = false;
  else {
    sscanf ( argv[check_arg ( argc, argv, "-T" ) + 1], "%lg", &T1 );
    if ( sscanf ( argv[min ( check_arg ( argc, argv, "-T" ) + 2, argc - 1 )], "%lg", &T2 ) != 1 )
      T2 = T1;
  }
  if ( check_arg ( argc, argv, "-wO2" ) == 0 )
    check = false;
  else {
    sscanf ( argv[check_arg ( argc, argv, "-wO2" ) + 1], "%lg", &w_O2_1 );
    if ( sscanf ( argv[min ( check_arg ( argc, argv, "-wO2" ) + 2, argc - 1 )], "%lg", &w_O2_2 ) != 1 )
      w_O2_2 = w_O2_1;
  }
  if ( check_arg ( argc, argv, "-wO" ) == 0 )
    check = false;
  else {
    sscanf ( argv[check_arg ( argc, argv, "-wO" ) + 1], "%lg", &w_O_1 );
    if ( sscanf ( argv[min ( check_arg ( argc, argv, "-wO" ) + 2, argc - 1 )], "%lg", &w_O_2 ) != 1 )
      w_O_2 = w_O_1;
  }
  if ( check_arg ( argc, argv, "-steps" ) == 0 )
    steps = 1;
  else
    sscanf ( argv[check_arg ( argc, argv, "-steps" ) + 1], "%d", &steps );

  if ( check == false ) {
    fprintf ( stderr, "\nFlags Missing!!\n\n"
              "Required Flag's:\n\n"
              "Flag      \tArg                 \tType    \tRequired? \n"
              "---------------------------------------------------------------\n"
              "-P    \t\t<Pressure(Pa)>        \tdouble   \tY\n"
              "-T    \t\t<Temperature(K)>       \tdouble   \tY\n"
              "-wO2  \t\t<Mass fraction of O2> \tdouble  \tY\n"
              "-wO  \t\t<Mass fraction of O>   \tdouble  \tY\n"
              "-steps    \t<number of steps>      \tint      \tN\n" "\n" );
  }

  else {
    printf ( "\n" );
    printf
      ( "P [Pa]        T [K] \t    wO2 \t  wO \t        Sebacher-     Macheret[s]   Sebacher-     Millikan[s]\n" );
    printf ( "\t\t\t\t\t\t\texpansion[s]                compression[s]\n" );
    printf ( "\n" );
    if ( steps == 1 ) {
      P = P2;
      T = T2;
      w_O2 = w_O2_2;
      w_O = w_O_2;
      tau_vt_millikan = tau_millikan ( w_O2, w_O, P, T );
      tau_vt_macheret = tau_macheret ( w_O2, w_O, P, T );
      tau_sebacher_comp = tau_sebacher_compression ( w_O2, w_O, P, T );
      tau_sebacher_expan = tau_sebacher_expansion ( w_O2, w_O, P, T );
      printf ( "%E  %E  %E  %E  %E  %E  %E  %E\n", P, T, w_O2, w_O, tau_sebacher_expan, tau_vt_macheret,
               tau_sebacher_comp, tau_vt_millikan );
    } else {
      for ( cnt = 0; cnt < steps; cnt++ ) {
        T = ( ( ( double ) cnt * ( ( T2 - T1 ) / ( steps - 1 ) ) ) + T1 );
        P = ( ( ( double ) cnt * ( ( P2 - P1 ) / ( steps - 1 ) ) ) + P1 );
        w_O2 = ( ( ( double ) cnt * ( ( w_O2_2 - w_O2_1 ) / ( steps - 1 ) ) ) + w_O2_1 );
        w_O = ( ( ( double ) cnt * ( ( w_O_2 - w_O_1 ) / ( steps - 1 ) ) ) + w_O_1 );
        tau_vt_millikan = tau_millikan ( w_O2, w_O, P, T );
        tau_vt_macheret = tau_macheret ( w_O2, w_O, P, T );
        tau_sebacher_comp = tau_sebacher_compression ( w_O2, w_O, P, T );
        tau_sebacher_expan = tau_sebacher_expansion ( w_O2, w_O, P, T );
        printf ( "%E  %E  %E  %E  %E  %E  %E  %E\n", P, T, w_O2, w_O, tau_sebacher_expan, tau_vt_macheret,
                 tau_sebacher_comp, tau_vt_millikan );
      }
    }
  }
  return 0;
}
