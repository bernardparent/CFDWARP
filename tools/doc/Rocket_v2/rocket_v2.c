/*
  This program called rocket_v2.c designs an axisymmetric rocket nozzle.  There
are two main computational components.  
        
  The first is the calculation of the rocket chamber conditons by assuming 
equilibrium combustion of one of four possible mixtures:

(A) Kerosene[C_12H_24] + Air 
(B) Kerosene + O_2
(C) H_2 + Air
(D) H_2 + O_2.   

  (More combustion mixutures can be added by creating an additional section in 
the "initialize" function without modifying any other subroutines.)  
NOTE:  the SPECIES_PRESENT and ATOMS_PRESENT values can be left at 19 and 4, 
but this will lead to slower convergence times than if the values are changed  
to reflect the lesser number of species and atoms required for any but the 
Kerosene + Air reaction. The post combustion composition is calculated using 
the Gibbs Minimization Technique which solves the system [Gibbs][X]=[B] where 
[Gibbs] is a square symmetric matrix of dimension "distinct atoms + 2". 
 
  The resutls of this analysis are sent to a Tecplot file "eqchemfieldinal.plt" 
which contains the final massfractions and kmol amounts of each species.
 
  The second main component of the code is the design of an axisymmetric 
nozzle using the assumption that the mixture is a perfect gas which
undergoes frozen expansion (i.e. constant ratio of specific heats, gamma, 
throughout the nozzle).  The design is accomplished by assuming radial
flow within a section of the nozzle and uses Kliegel's analysis for the 
throat region.  An upper contour is found which produces a specified
Mach number distribution along the axisymmetric axis.  This contour is then
corrected for the presence of a boundary layer using the method detailed
by Edenfield.  NOTE:  if the nozzle contour fails at low exit Mach numbers,
try increasing the in number of points along the centreline distribution
(i.e. IATOB, IBTOC, etc.) to ensure the characteristic mesh is tight enough
(as the Mach number approaches unity, the characteristics tend to vertical
hence making their intersection difficult to find unless the starting points
are sufficiently close to one another)  Also, if the contour portion fails
try increasing "jlimit" in the "findcontour" function. 

  The results of the nozzle design are sent to a Tecplot file "nozzlecont.plt"
which contains the inviscid as well as viscous profiles and the boundary
layer profile.  The file "centreveldist.plt" is a Tecplot file containing the
prescribed Mach number distribution along the nozzle axis which is output
to ensure a smooth transition between nozzle flow zones.  A summary of the
important design and calculated variables is output to a text file
"rocketsum.out".

  There is also a subroutine that can be used to design a contoured 
subsonic portion of the nozzle upstream of the throat.  This contour
is taken from Prof. Sislian's class notes.  The upper wall points of 
this section of the nozzle are written to the file "subnozzlecont.plt".

  Please see the document "Axisymmetric Rocket Nozzle Design" for 
further details and references.
*/

#include <math.h> /* rememeber to compile using gcc filename.c -lm 
		     (Note: log=ln) */
#include <stdio.h>

/******************************************************************************
 INPUT SECTION
******************************************************************************/

/*Equilibrium Combustion Input Parameters*/
/*---------------------------------------------------------*/

#define PI 3.14159265359
#define RUNIVERSAL 8314.3        /* [J/kgmol K] */
#define SPECIES 19   /*total number of possible species and atoms*/  
#define ATOMS 4      /*that can be considered */

#define SPECIES_PRESENT 19       /* These should be 19 and 4 for OXIDIZER=0, */
#define ATOMS_PRESENT 4          /* or 16 and 3 for OXIDIZER=1, or 6 and 2 for
				    OXIDIZER=2*/
#define INITIAL_TEMPERATURE 300  /* Fuel injection temperature, [K] */
#define INITIAL_PRESSURE 1.5E6   /* Combustor pressure, which is const [Pa] */
#define EQUIVALENCE_RATIO 1.0   
#define OXIDIZER 1               /* 0 = ker+air, 1 = ker+O2, 2 = H2+O2, 
				    3 = H2+air */
#define KGFUELUSED 1             /* this is the actual amount of fuel being 
				    injected [kg/s]
                                    for use in calculating actual mass flow */

/*Nozzle Contour Design Parameters*/
/*----------------------------------------------------------*/
/* Physical Variables */

#define MEXIT 5           /*design Mach# at exit, must set higher than 2.4*/
#define THROAT 1E-2       /*throat radius [m]*/
#define MAXANGLE 10E0     /*maximum expansion angle of nozzle (pt C)[degrees]*/
#define VISCOUS 1         /*set to 1 for bndry layer correction*/
#define TWALL 800         /*nozzle wall temp [K] (needed for bndry layer)*/
#define SUBSONIC 1        /*set to 1 to calculate a subsonic contour*/
#define CHAMBERM 0.1      /*mach number at end of combustion chamber*/
#define SUBLENGTH .25     /*length of subsonic portion of nozzle as percentage
                            of total supersonic length*/


/* Numerical Variables */

#define IBTOA 200       /*num nodes between points A and B*/
#define ICTOB 100       /*num nodes between points B and C*/
#define IDTOC 100       /*num nodes between points C and D*/
#define IFTOD 100       /*num nodes between points D and F*/
#define JTOT 100        /*num nodes between wall and axis of sym at exit*/
#define MAXJARRAY 700   /*num nodes the characteristic mesh is allowed to 
			  grow to in the j direction*/
#define ISUB 100        /*number of points along the subsonic contour*/

/******************************************************************************
 DECLARE STRUCT TYPES
******************************************************************************/

struct eqmodel  
  {  
    double kmol[SPECIES];/* number of elements in array = number of species */
    double kmolsum;      /* total number of kmols,=sum of kmol[k], k=1->ns */
    double temp;
    double pres;         /* the value stored here is in atm */
    double h[SPECIES];
    double cp[SPECIES];
    double s[SPECIES];
    double g[SPECIES];
    double pi[ATOMS];
    double dlnkmol[SPECIES];
    double dlnkmolsum;
    double dlntemp;
  };

struct eqlist
  {
    double totatoms[ATOMS];
    double ftoaair;          /* [kmolf/kmol], by mol */
    double ftoao2;           /* [kmolf/kmolo2], by mol */
    double fh2toao2;         /* [kmolh2/kmolo2], by mol */
    double fh2toair;         /* [kmolh2/kmol], by mol */
    double atof;             /* [kga/kgf] OR [kgo2/kgf] */
    double phi;
    double Molwair;          /* [kgair/kmolair] */
    double hstart;           /* [J/kgmix] */
    double Rmix;             /* Rmix and Cpmix pertain to whole mixture */   
  };

struct eqout
  {
    double mass[SPECIES];
    double massfrac[SPECIES];
    double gammamix;
    double Cpmix;
    double rhotot;
    double masstot;
    double Molwmix;
    double Rgasmix;
    double Tstar;
    double pstar;
    double fuelin;
    double massflowout;
    double achamber;
  };

struct nozmodel
  {
    double w[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double rorcr[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double R, K, epsbar, hexitbar, hthroatbar;
    double wa, wb, wc, wd, wf;
    double ra, rb, rc, rd;
    double xdforcr, xaborcr, xincorcr;
    double thetab, thetac, thetad;/*these are integrated angles, not physical*/
    double psibarc;
  };

struct nozlist
  {
    double gamma, Rmix, ttot, ptot, mexit;
    double hthroat, thetawall;
    double rhotot, Cpmix, astar, rhostar, rhostarfrac, tstar, pstar;
    double texit, pexit, aexit, achamber;
  };

struct nozout
  {
    double rcr, hexit, hchamb, Rave;
    double xworcr[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double yworcr[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double thetaw[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double qw[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double xw[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double yw[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double ybound[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double Reyinv[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double Tinv[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double Reyref[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double Tref[IBTOA+ICTOB+IDTOC+IFTOD+1];
    double xsub[ISUB];
    double ysub[ISUB];
  };

struct moc
  {
    double x4, y4, u4, v4, a4, m4;
  }; 

struct slope
  {
    double bprime, bdprime;
    double aprime, adprime, atprime;
  };

struct poly
  {
    double D[5];
    double C[5];
  };


/******************************************************************************
   PROTOTYPING OF FUNCTIONS USED IN THE PROGRAM (AS THEY DO NOT RETURN THE 
   STANDARD INTEGER VALUES)
******************************************************************************/

/*Equilibrium Chemistry functions-------------------------------------*/
void eqinitialize(int choice, int ns, struct eqmodel work[1], struct eqlist 
		constant[1], double Molw[SPECIES], double a[SPECIES][2][7], 
		double Rgas);
void makeGibbs(int ns,int na,struct eqmodel work[1], struct eqlist constant[1], 
	       double Gibbs[ATOMS+2][ATOMS+2], 
	       double eta[ATOMS][SPECIES], double Molw[SPECIES], 
	       double a[SPECIES][2][7], double Rgas);
void makeB(int ns, int na, struct eqmodel work[1], struct eqlist constant[1], 
	   double B[ATOMS+2], double eta[ATOMS][SPECIES],
	   double Molw[SPECIES], double a[SPECIES][2][7], double Rgas);
void gauss(int dim, double A[ATOMS+2][ATOMS+2], double B[ATOMS+2]); 
void extract(int na, int ns, struct eqmodel work[1], double B[na+2], 
	     double eta[ATOMS][SPECIES], double Molw[SPECIES], 
	     double a[SPECIES][2][7], double Rgas);
void update(int ns, struct eqmodel work[1]);
void properties(int ns, double Molw[SPECIES], struct eqmodel work[1], 
		struct eqlist constant[1], struct eqout eqfinal[1], 
		double a[SPECIES][2][7], double Rgas);
void printscreen1(struct eqmodel work[1], struct eqout eqfinal[1], 
		  struct eqlist constant[1], int ns, FILE *summary);
void printoutput1(struct eqmodel work[1], struct eqout eqfinal[1], 
		  FILE *output, int numspec);

double enthalpyfind(double Rgas, double Molw[SPECIES], double a[SPECIES][2][7],
		    double temp, int spec);
double cpfind(double Rgas, double Molw[SPECIES], double a[SPECIES][2][7], 
	      double temp, int spec);
double entropyfind(double Rgas, double Molw[SPECIES], double a[SPECIES][2][7],
		   double temp, int spec);
double converge(int ns, struct eqmodel work[1]);

/*Nozzle functions-----------------------------------------------*/
void nozinitialize(struct nozlist *set, struct eqmodel work[1], 
		   struct eqout eqeqfinal[1]);
void finddfpoly(struct nozlist set, struct nozmodel *node, struct slope *deriv,
	       struct poly *coeff);
void radialbounds(struct nozlist set, struct nozmodel *node);
void findabpoly(struct nozlist set, struct nozmodel *node, struct slope *deriv,
	       struct poly *coeff);
void radialveldist(struct nozlist set, struct nozmodel *node);
void dfveldist(struct nozmodel *node, struct poly coeff);
void abveldist(struct nozmodel *node, struct poly coeff);
void findcontour(struct nozlist set,struct nozmodel node, struct nozout *final);
void findmocvalues(double u1, double V1, double a1, double m1, double y1, 
		   double u2, double v2, double a2, double m2, double y2,
		   double x1, double x2, int j, struct nozlist set,
		   struct moc *mocvars);
void boundarylayer(struct nozlist set, struct nozout *final,
		   struct eqmodel work[1], double Molw[SPECIES], 
		   double a[SPECIES][2][7], double epsilonk[SPECIES], 
		   double sigmak[SPECIES], double omega11polyd[3][5], 
		   double astarpolye[3][3], struct eqout eqfinal[1], 
		   double Rgas, int ns);
void findsubcontour(struct nozlist set, struct nozout *final);
double findw(double gamma, double mach);
double finda(struct nozlist set, double w);
double findm(struct nozlist set, double w);
double findrhofrac(struct nozlist set, double w); 
double findrbar(struct nozlist set, double w);
double findangle(struct nozlist set, double w);
double findpsibar(struct nozlist set, double w, double theta);
double findmassflow(double gamma, double S);
double findhexitbar(struct nozlist set, struct nozmodel node);
double findrcr(struct nozlist set, struct nozmodel node);
double findwprime(double w, double rbar, double gamma);
double findwdbleprime(double w, double rbar, double gamma);
double findsonicprime(struct nozlist set, double S, double R1);
double findsonicdbleprime(struct nozlist set, double S, double R1);
double findsonictpleprime(struct nozlist set, double S, double R1);
double findsonicpt(struct nozlist set, double S, double R1);
double findxab(struct slope deriv, struct nozmodel node);
double findxaberror(struct slope deriv, struct nozmodel node);
double findomega11(double omega11polyd[3][5], double epsilonk[SPECIES], 
		   double T, int species);
double findmu(double epsilonk[SPECIES], double sigmak[SPECIES],
	      double omega11polyd[3][5], double astarpolye[3][3],
	      double Molw[SPECIES], struct eqout eqfinal[1], double T, int ns);

void screenoutput(struct nozlist set, struct nozout final, FILE *summary);
void screenoutput2(struct nozlist set,struct nozmodel node, struct nozout final,
		   struct poly coeff, FILE *summary);
void finaloutput(FILE *output, struct nozout final);
void subsonicoutput(FILE *output, struct nozout final);
void velocityoutput(FILE *output, struct nozlist set, struct nozmodel node,
		    struct nozout final);




/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
main()
{

  /* REQUIRED POLYNOMIALS FOR FINDING VARIOUS QUANTITIES OF INTEREST */
   
  /* The species used in this program are numbered as follows */
  /*   0. H2       5. H2O       10. CO2         15. C12H24    */
  /*   1. O2       6. HO2       11. CH3         16. N         */
  /*   2. H        7. C2H2      12. CH4         17. NO        */
  /*   3. O        8. C2H4      13. H2CO        18. N2        */
  /*   4. OH       9. CO        14. C6H120                    */
     
  /* Molecular weight of species k, [kgk/kmol] */
  double Molw[SPECIES] = {2.01588, 31.9988, 1.00794, 15.9994, 17.00734,
			  18.01528, 33.00674, 26.03788, 28.05376, 28.0104, 
			  44.0098, 15.03482, 16.04276, 30.02628, 100.16068, 
			  168.32256, 14.0067, 30.00061, 28.0134};

  /* The next four polynomials are used for finding viscosity, taken
     from Gardiner */
  double epsilonk[SPECIES] = 
   {59.7e0, 106.7e0, 37.0e0, 106.7e0, 79.8e0, 260.0e0, 168.0e0, 231.8e0, 
    224.7e0, 91.7e0, 195.2e0, 312.0e0, 148.6e0, 99, 99, 99, 71.4e0, 119.0e0, 
    71.4e0};

  double sigmak[SPECIES]=
   {0.2827e0, 0.3467e0, 0.2070e0, 0.3050e0, 0.3147e0, 0.2800e0, 0.3068E+0,
    0.4033e0, 0.4163e0, 0.3690e0, 0.3941e0, 0.3644e0, 0.3758e0, 99, 99, 99, 
    0.3798e0, 0.3470e0, 0.3798e0};

  double omega11polyd[3][5]=
  {
   {2.3527333E+0, -1.3589968E+0, 5.2202460E-1, -9.4262883E-2, 6.4354629E-3},
   {1.2660308E+0, -1.6441443E-1, 2.2945928E-2, -1.6324168E-3, 4.5833672E-5},
   {8.5263337E-1, -1.3552911E-2, 2.6162080E-4, -2.4647654E-6, 8.6538568E-9}
  };

  double astarpolye[3][3]=
  {
   {1.1077725E+0, -9.4802344E-3, +1.6918277E-3},
   {1.0871429E+0, +3.1964282E-3, -8.9285689E-5},
   {1.1059000E+0, +6.5136364E-4, -3.4090910E-6}
  };


  /* Number of kg-atoms of atom i per kmol of species k [katoms/kmol] 
     (Note: katoms = kg-atoms) */
  double eta[ATOMS][SPECIES] = 
  {
    /* atom 0, H */
    {2.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 2.0, 4.0, 0.0, 0.0, 
     3.0, 4.0, 2.0, 12.0, 24.0, 0.0, 0.0, 0.0},
    /* atom 1, O */
    {0.0, 2.0, 0.0, 1.0, 1.0, 1.0, 2.0, 0.0, 0.0, 1.0, 2.0, 
     0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0},
    /* atom 2, C */
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 
     1.0, 1.0, 1.0, 6.0, 12.0, 0.0, 0.0, 0.0},
    /* atom 3, N */
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0}
  };


  /* NASA polynomial co-efficient matrix "a" taken from NASA TM 4513 (McBride,
     Gordon, and Reno)
     Note: the underscored species (i.e. _C12H24) are taken from Jurgen's 
     19 species kerosene code*/
  double a[SPECIES][2][7]=
  {
    /* species 0, H2 */
    {
      {+2.34433112E+0, +7.98052075E-3, -1.94781510E-5, +2.01572094E-8, 
       -7.37611761E-12, -9.17935173E+2, +6.83010238E-1},
      {+2.93286579E+0, +8.26607967E-4, -1.46402335E-7, +1.54100359E-11, 
       -6.88804432E-16, -8.13065597E+2, -1.02432887E+0}
    },
    /* species 1, O2 */
    {
      {+3.78245636E+0, -2.99673415E-3, +9.84730200E-6, -9.68129508E-9, 
       +3.24372836E-12, -1.06394356E+3, +3.65767573E+0},
      {+3.66096083E+0, +6.56365523E-4, -1.41149485E-7, +2.05797658E-11, 
       -1.29913248E-15, -1.21597725E+3, +3.41536184E+0}
    },
    /* species 2, H */  
    {
      {+2.50000000E+0, +0.00000000E+0, +0.00000000E+0, +0.00000000E+0, 
       +0.00000000E+0, +2.54736599E+4, -4.46682853E-1},
      {+2.50000286E+0, -5.65334214E-9, +3.63251723E-12, -9.19949720E-16, 
       +7.95260746E-20, +2.54736589E+4, -4.46698494E-1}
    },
    /* species 3, O */
    {
      {+3.16826710E+0, -3.27931884E-3, +6.64306396E-6, -6.12806624E-9, 
       +2.11265971E-12, +2.91222592E+4, +2.05193346E+0},
      {+2.54363697E+0, -2.73162486E-5, -4.19029520E-9, +4.95481845E-12, 
       -4.79553694E-16, +2.92260120E+4, +4.92229457E+0}
    },
    /* species 4, OH */
    {
      {+3.99201543E+0, -2.40131752E-3, +4.61793841E-6, -3.88113333E-9, 
       +1.36411470E-12, +3.61508056E+3, -1.03925458E-1},
      {+2.83864607E+0, +1.10725586E-3, -2.93914978E-7, +4.20524247E-11, 
       -2.42169092E-15, +3.94395852E+3, +5.84452662E+0}
    },
    /* species 5, H2O */
    {
      {+4.19864056E+0, -2.03643410E-3, +6.52040211E-6, -5.48797062E-9, 
       +1.77197817E-12, -3.02937267E+4, -8.49032208E-1},
      {+2.67703787E+0, +2.97318329E-3, -7.73769690E-7, +9.44336689E-11, 
       -4.26900959E-15, -2.98858938E+4, +6.88255571E+0}
    },
    /* species 6, HO2 */
    {
      {+4.30179801E+0, -4.74912051E-3, +2.11582891E-5, -2.42763894E-8,
      +2.29225124E-12, +2.94808040E+2, +3.71666245E+0},
      {+4.17228728E+0, +1.88117647E-3, -3.46277408E-7, +1.94657853E-11,
      +1.76254294E-16, +6.18102964E+1, +2.95767746E+0}
    },
    /* species 7, C2H2 (acetylene) */
    {
      {+8.08681094E-1, +2.33615629E-2, -3.55171815E-5, +2.80152437E-8,
       -8.50072974E-12, +2.64289807E+4, +1.39397051E+1},
      {+4.65878504E+0, +4.88396547E-3, -1.60828775E-6, +2.46974226E-10, 
       -1.38605680E-14, +2.57594044E+4, -3.99834772E+0}
    },
    /* species 8, C2H4 */
    {
      {+3.95920148E+0, -7.57052247E-3, +5.70990292E-5, -6.91588753E-8, 
       +2.69884373E-11, +5.08977593E+3, +4.09733096E+0},
      {+3.99182761E+0, +1.04833910E-2, -3.71721385E-6, +5.94628514E-10, 
       -3.53630526E-14, +4.26865819E+3, -2.69052151E-1}
    },
    /* species 9, CO */
    {
      {+3.57953347E+0, -6.10353680E-4, +1.01681433E-6, +9.07005884E-10, 
       -9.04424499E-13, -1.43440860E+4, +3.50840928E+0},
      {+3.04848583E+0, +1.35172818E-3, -4.85794075E-7, +7.88536486E-11,
       -4.69807489E-15, -1.42661171E+4, +6.01709790E+0}
    },
    /* species 10 _CO2 */
    { 
      {0.24007797E+01, 0.87350957E-02,-0.66070878E-05, 0.20021861E-08, 
       0.63274039E-15,-0.48377527E+05, 0.96951457E+01},
      {0.44608041E+01, 0.30981719E-02,-0.12392571E-05, 0.22741325E-09,
       -0.15525954E-13,-0.48961442E+05,-0.98635982E+00}
    },
    /* species 11, CH3 */
    {
      {+3.67359040E+0, +2.01095175E-3, +5.73021856E-6, -6.87117425E-9,
       +2.54385734E-12, +1.64449988E+4, +1.60456433E+0},
      {+2.96866033E+0, +5.80717546E-3, -1.97778534E-6, +3.07278752E-10,
       -1.78853897E-14, +1.65388869E+4, +4.77944503E+0}
    },
    /* species 12, CH4 */
    {
      {+5.14987613E+0, -1.36709788E-2, +4.91800599E-5, -4.84743026E-8,
       +1.66693956E-11, -1.02466476E+4, -4.64130376E+0},
      {+1.63552643E+0, +1.00842795E-2, -3.36916254E-6, +5.34958667E-10,
       -3.15518833E-14, -1.00056455E+4, -9.99313326E+0}
    },
    /* species 13, H2CO (formaldehyde HCHO) */
    {
      {+4.79372315E+0, -9.90833369E-3, +3.73220008E-5, -3.79285261E-8, 
       +1.31772652E-11, -1.43089567E+4, +6.02812900E-1},
      {+3.16952654E+0, +6.19320583E-3, -2.25056377E-6, +3.65975680E-10,
       -2.20149470E-14, -1.44784444E+4, +6.04209449E+0}
    },
    /* species 14, _C6H12O */
    {
      {+4.6835000E-1, +6.8027000E-2, -3.9907000E-5, +9.8998000E-10, 
       +0.0000000E+0, -3.6987000E+4, +2.7960000E+1},
      {+2.9578000E+1, +0.0000000E-0, +0.0000000E-0, +0.0000000E-0, 
       +0.0000000E+0, +4.5138300E+4, -1.2478900E+2}
    },
    /* species 15, _C12H24 (kerosene) */
    {
      {+4.9915568E-3, +1.0064852E-4, -1.5911406E-8, -2.8208800E-11, 
       +1.1964748E-14, -1.9543438E+1, +1.6904105E-2},
      {+7.3484623E-2, +0.0000000E+0, +0.0000000E+0, +0.0000000E+0,
       +0.0000000E+0, -4.7675298E+1, -3.6994817E-1}
    },
    /* species 16, N */
    {
      {+2.50000000E+0, +0.00000000E+0, +0.00000000E+0, +0.00000000E+0,
       +0.00000000E+0, +5.61046378E+4, +4.19390932E+0},
      {+2.41594293E+0, +1.74890600E-4, -1.19023667E-7, +3.02262387E-11,
       -2.03609790E-15, +5.61337748E+4, +4.64960986E+0}
    },
    /* species 17, NO */  
    {
      {+4.21859896E+0, -4.63988124E-3, +1.10443049E-5, -9.34055507E-9,
       +2.80554874E-12, +9.84509964E+3, +2.28061001E+0},
      {+3.26071234E+0, +1.19101135E-3, -4.29122646E-7, +6.94481463E-11,
       -4.03295681E-15, +9.92143132E+3, +6.36900518E+0}
    },
    /* species 18, N2 */ 
    {
      {+3.53100528E+0, -1.23660987E-4, -5.02999437E-7, +2.43530612E-9, 
       -1.40881235E-12, -1.04697628E+3, +2.96747468E+0},
      {+2.95257626E+0, +1.39690057E-3, -4.92631691E-7, +7.86010367E-11, 
       -4.60755321E-15, -9.23948645E+2, +5.87189252E+0}
       }
  };

  /***************************************************************************
   PROGRAM BEGINS HERE
  ****************************************************************************/
 
  FILE *summary;
  FILE *eqchemfieldinal;
  FILE *centreline;
  FILE *wallcontour;
  FILE *subcontour;
  int choice = OXIDIZER;      /* Choice = 0 for air, = 1 for pure O2 */
  int count = 0;              /* counts the iterations performed */
  int ns=SPECIES_PRESENT;     /* these values are 19 and 4 for combust. w/air*/
  int na=ATOMS_PRESENT;       /* or 16 and 3 for combustion w/ pure oxygen */
  double tolerance = 1E-7;    /* this is the convergence criteria */
  double test = tolerance+1;  /* convergence tester */
  double Rgas = RUNIVERSAL;   /* [J/kgmol K], note kgmol = kmol */

  int k, p;
  struct eqmodel work[1];
  struct eqlist constant[1];
  struct eqout eqfinal[1];
  double Gibbs[na+2][na+2];     /* the coefficient matrix to be solved */
  double B[na+2];               /* Gibbs*X=B */
  double X[na+2];
  struct nozlist set;
  struct nozmodel node;
  struct poly coeff;
  struct slope deriv;
  struct nozout final; 

  summary = fopen("rocketsum.out", "w");
  eqchemfieldinal = fopen("eqchemfieldinal.plt", "w");
  centreline = fopen("centveldist.plt", "w");
  wallcontour = fopen("nozzlecont.plt", "w");
  subcontour = fopen("subnozzlecont.plt", "w");

  constant[0].phi = EQUIVALENCE_RATIO;     
  work[0].temp = INITIAL_TEMPERATURE;       
  work[0].pres = INITIAL_PRESSURE/101325;   /* calculated as [atm] */
  eqfinal[0].fuelin = KGFUELUSED;

  eqinitialize(choice, ns, work, constant, Molw, a, Rgas);
  work[0].temp=3800;                 /*Initial equilibrium temperature guess */
  printf("\n%d\t%.5f\t%.5f",count,test,work[0].temp);
  while(test>tolerance) 
    /*for(k=0; k<100; k++) for use if convergence is not being obtained */ 
    {
      makeGibbs(ns, na, work, constant, Gibbs, eta, Molw, a, Rgas);
      makeB(ns, na, work, constant, B, eta, Molw, a, Rgas);
      gauss(na+2, Gibbs, B);
      extract(na, ns, work, B, eta, Molw, a, Rgas);
      update(ns, work);
      test=converge(ns, work);
      count = count + 1;
      printf("\n%d\t%.7f\t%.5f",count,test,work[0].temp);
    }
  properties(ns, Molw, work, constant, eqfinal, a, Rgas);
  printscreen1(work, eqfinal, constant, ns, summary);
  printoutput1(work, eqfinal, eqchemfieldinal, ns);

  /*At this point the nozzle chamber conditions have been determined,
    with the total temperature being held in work[0].temp and the 
    total pressure being INITIAL_PRESSURE (since constant pressure 
    combustion is assumed.  The mixture properties are held in 
    the array eqfinal[0]*/

  nozinitialize(&set, work, eqfinal);
  finddfpoly(set, &node, &deriv, &coeff);
  radialbounds(set, &node);
  findabpoly(set, &node, &deriv, &coeff);
  dfveldist(&node, coeff);
  radialveldist(set, &node);
  abveldist(&node, coeff);

  node.hexitbar=findhexitbar(set, node);
  final.rcr=findrcr(set, node);
  final.hexit=node.hexitbar*final.rcr;
  node.hthroatbar=set.hthroat/final.rcr;
  
  findcontour(set, node, &final);
  if(VISCOUS==1)
    {boundarylayer(set, &final, work, Molw, a, epsilonk, sigmak, 
		   omega11polyd, astarpolye, eqfinal, Rgas, ns);}
  if(SUBSONIC==1)
    {
      findsubcontour(set, &final);
      subsonicoutput(subcontour, final);
    }
  screenoutput(set, final, summary);
  screenoutput2(set, node, final, coeff, summary);
  velocityoutput(centreline, set, node, final); 
  finaloutput(wallcontour, final);

} 
/******************************************************************************
  PROGRAM ENDS HERE
******************************************************************************/

double enthalpyfind(double Rgas, double Molw[SPECIES], double a[SPECIES][2][7],
		    double T, int spec) /*[J/kmolk]*/
  {     
    int index = 1;     
    double enth;
    if (T < 1000.00)
      {index = 0;}
    enth = Rgas*T*( a[spec][index][0] + a[spec][index][1]*T/2 + 
		    a[spec][index][2]*T*T/3 + 
		    a[spec][index][3]*T*T*T/4 +
		    a[spec][index][4]*T*T*T*T/5 + 
		    a[spec][index][5]/T );
    return(enth); /* Note: enth/Molw[spec] yields units of [J/kgk] */
  }

double cpfind(double Rgas, double Molw[SPECIES], double a[SPECIES][2][7], 
	      double T, int spec) /* [J/kmolk K] */
  {
    int index = 1;
    double cp;
    if (T < 1000.00)
      {index = 0;}
    cp = Rgas*(a[spec][index][0] + a[spec][index][1]*T + 
	       a[spec][index][2]*T*T + 
	       a[spec][index][3]*T*T*T + 
	       a[spec][index][4]*T*T*T*T);
    return(cp); /* Note: cp/Molw[spec] yields units of [J/kgk K] */
  }

double entropyfind(double Rgas, double Molw[SPECIES], double a[SPECIES][2][7],
		   double T, int spec)/*[J/kmolk K]*/
  {
    int index = 1;
    double entr;
    if (T < 1000.00)
      {index = 0;}
    entr = Rgas*(a[spec][index][0]*log(T) + a[spec][index][1]*T + 
		 a[spec][index][2]*T*T/2 +
		 a[spec][index][3]*T*T*T/3 + 
		 a[spec][index][4]*T*T*T*T/4 + 
		 a[spec][index][6]);
    return(entr); /* Note: entr/Molw[spec] yields units of [J/kgk K] */
  }

void eqinitialize(int choice, int ns, struct eqmodel work[1], struct eqlist 
		constant[1], double Molw[SPECIES], double a[SPECIES][2][7], 
		double Rgas)
  {
    int k;
    double hold, hold2, hold3;
    /* stoichiometric fuel to air ratio using air, by mol */   
    constant[0].ftoaair = 1/(4.76*18.0);   
    /* stoichiometric fuel to air ratio using O2, by mol */
    constant[0].ftoao2 = 1/18.0;    
    /* stoichiometric fuel to air ratio using H2 and O2, by mol */
    constant[0].fh2toao2 = 2;
    /* stoichiometric fuel to air ratio using H2 and air, by mol */
    constant[0].fh2toair = 2/4.76;
    /* (21/79:O2/N2) [kgair/kmolair] */
    constant[0].Molwair = 4.76*(0.21*Molw[1] + 0.79*Molw[18]);  
    /* NOTE: 1 kmol of air has 4.76 kmols in it (ie. 4.76 [kmol/kmolair]
             *[kgair/kmol] = [kgair/kmolair]) */
    work[0].kmolsum = 0.1; /* Initial solution for final kmols of species */ 
    for(k=0; k<ns; k++)
      {
	work[0].kmol[k]=work[0].kmolsum/ns;
      }
/* -------------------------------------------------------------------------
   Calculate the total number of atoms of each type per kg of mixture 
   (fuel plus either air or O2)
   -----------------------------------------------------------------------*/

    if(choice==0) /*phi*C12H24 + 18(O2 + 3.76N2) --> 12H2O + 12CO2 + 67.68N2*/
      {
	constant[0].atof=1/(constant[0].phi*constant[0].ftoaair*Molw[15]/
			   (constant[0].Molwair/4.76));           
	/*[kgair/kgf]*/ 
	constant[0].totatoms[0]=(24/Molw[15])/(1+constant[0].atof);
	/*[katomsH/kgmix]*/
	constant[0].totatoms[1]=(2/constant[0].Molwair)*constant[0].atof/
	  (1+constant[0].atof);                                 
	/*[katomsO/kgmix]*/
	constant[0].totatoms[2]=(12/Molw[15])/(1+constant[0].atof);
	/*[katomsC/kgmix]*/
	constant[0].totatoms[3]=(2*3.76)/constant[0].Molwair*constant[0].atof/
	  (1+constant[0].atof);               
	/*[katomsN/kgmix]*/

	/* Calculate the total enthalpy of the mixture */
	hold=0; 
	hold2=0; 
	hold3=0;
	constant[0].hstart=0;
	hold=enthalpyfind(Rgas, Molw, a, work[0].temp, 1);  /*[J/kmolO2]*/ 
	hold2=enthalpyfind(Rgas, Molw, a, work[0].temp, 18); /*[J/kmolN2]*/
	hold3=enthalpyfind(Rgas, Molw, a, work[0].temp, 15);/*[J/kmolC12H24]*/
	constant[0].hstart=((1/constant[0].Molwair)*constant[0].atof*hold + 
			 (3.76/constant[0].Molwair)*constant[0].atof*hold2 + 
			 (1/Molw[15])*hold3)/(1+constant[0].atof);/*[J/kgmix]*/
      }
    if(choice==1) /* phi*C12H24 + 18O2 --> 12H2O + 12CO2 */
      {
	constant[0].atof=1/(constant[0].phi*constant[0].ftoao2*Molw[15]/
			    Molw[1]);   
	/* [kgo2/kgf] */ 
	constant[0].totatoms[0]=(24/Molw[15])/(1+constant[0].atof);      
	/* [katomsH/kgmix] */
	constant[0].totatoms[1]=(2/Molw[1])*constant[0].atof/
	  (1+constant[0].atof);  
	/* [katomsO/kgmix] */
	constant[0].totatoms[2]=(12/Molw[15])/(1+constant[0].atof);       
	/* [katomsC/kgmix] */
	constant[0].totatoms[3]=0;               
	/*[katomsN/kgmix]*/

	/* Calculate the total enthalpy of the mixture */
	
	hold=0; 
	hold3=0;
	constant[0].hstart=0;
	hold=enthalpyfind(Rgas, Molw, a, work[0].temp, 1);      
	/* [J/kmolO2] */ 
	hold3=enthalpyfind(Rgas, Molw, a, work[0].temp, 15);     
	/* [J/kmolC12H24] */
	constant[0].hstart=((1/Molw[1])*constant[0].atof*hold + 
			 (1/Molw[15])*hold3)/(1+constant[0].atof); 
	/* [J/kgmix] */
      }

    if(choice==2) /* phi*2H2 + O2 --> 2H2O */
      {
	constant[0].atof=1/(constant[0].phi*constant[0].fh2toao2*Molw[0]/
			   Molw[1]);           
	/*[kgair/kgf]*/ 
	constant[0].totatoms[0]=(2/Molw[0])/(1+constant[0].atof);
	/*[katomsH/kgmix]*/
	constant[0].totatoms[1]=(2/Molw[1])*constant[0].atof/
	  (1+constant[0].atof);                                 
	/*[katomsO/kgmix]*/
	constant[0].totatoms[2]=0;
	/*[katomsC/kgmix]*/
	constant[0].totatoms[3]=0;               
	/*[katomsN/kgmix]*/

	/* Calculate the total enthalpy of the mixture */
	hold=0; 
	hold2=0; 
	constant[0].hstart=0;
	hold=enthalpyfind(Rgas, Molw, a, work[0].temp, 1);       
	/*[J/kmolO2]*/ 
	hold2=enthalpyfind(Rgas, Molw, a, work[0].temp, 0);      
	/*[J/kmolH2]*/
	constant[0].hstart=((1/Molw[1])*constant[0].atof*hold + 
			 (1/Molw[0])*hold2)/(1+constant[0].atof); 
	/*[J/kgmix]*/
      }

    if(choice==3) /* phi*2H2 + (O2 + 3.76N2) --> 2H2O + 3.76N2 */
      {
	constant[0].atof=1/(constant[0].phi*constant[0].fh2toair*Molw[0]/
			   (constant[0].Molwair/4.76));           
	/*[kgair/kgf]*/ 
	constant[0].totatoms[0]=(2/Molw[0])/(1+constant[0].atof);
	/*[katomsH/kgmix]*/
	constant[0].totatoms[1]=(2/constant[0].Molwair)*constant[0].atof/
	  (1+constant[0].atof);                                 
	/*[katomsO/kgmix]*/
	constant[0].totatoms[2]=0;
	/*[katomsC/kgmix]*/
	constant[0].totatoms[3]=(2*3.76)/constant[0].Molwair*constant[0].atof/
	  (1+constant[0].atof);               
	/*[katomsN/kgmix]*/

	/* Calculate the total enthalpy of the mixture */
	hold=0; 
	hold2=0; 
	hold3=0;
	constant[0].hstart=0;
	hold=enthalpyfind(Rgas, Molw, a, work[0].temp, 1);  /*[J/kmolO2]*/ 
	hold2=enthalpyfind(Rgas, Molw, a, work[0].temp, 18); /*[J/kmolN2]*/
	hold3=enthalpyfind(Rgas, Molw, a, work[0].temp, 0);/*[J/kmolH2]*/
	constant[0].hstart=((1/constant[0].Molwair)*constant[0].atof*hold + 
			 (3.76/constant[0].Molwair)*constant[0].atof*hold2 + 
			 (1/Molw[0])*hold3)/(1+constant[0].atof);/*[J/kgmix]*/
      }
  }

void makeGibbs(int ns, int na, struct eqmodel work[1],struct eqlist constant[1],
	       double Gibbs[na+2][na+2], double eta[ATOMS][SPECIES],
	       double Molw[SPECIES], double a[SPECIES][2][7], double Rgas)
  {
    int i, j, k;        
    double hold, hold2, hold3;

    /* calculate the various thermodynamic variables required from the
       NASA polynomials */
    for(k=0; k<ns; k++)
      {
	work[0].h[k]=0;
	work[0].cp[k]=0;
	work[0].s[k]=0;
	work[0].g[k]=0;
      }
    for(k=0; k<ns; k++)
      {
	work[0].h[k]=enthalpyfind(Rgas, Molw, a, work[0].temp, k);
	work[0].cp[k]=cpfind(Rgas, Molw, a, work[0].temp, k);
	work[0].s[k]=entropyfind(Rgas, Molw, a, work[0].temp, k);
	work[0].g[k]=work[0].h[k]-work[0].temp*work[0].s[k];
      }
    /* first 0-->(na-1) rows of the matrix */
    for(i=0; i<na; i++)
      {
	for(j=0; j<na; j++)
	  {
	    hold=0;
	    for(k=0; k<ns; k++)
	      {
		hold=hold+eta[i][k]*eta[j][k]*work[0].kmol[k];}
	    Gibbs[i][j]=hold;
	  }
	hold2=0;
	hold3=0;
	for(k=0; k<ns; k++)
	  {
	    hold2=hold2+eta[i][k]*work[0].kmol[k];
	    hold3=hold3+eta[i][k]*work[0].kmol[k]*work[0].h[k]/
	      (Rgas*work[0].temp);
	  }
	Gibbs[i][na]=hold2;
	Gibbs[i][na+1]=hold3;
      }
    /* row na of the matrix */
    for(j=0; j<na; j++)
      {
	hold=0;
	for(k=0; k<ns; k++)
	  {
	    hold=hold+eta[j][k]*work[0].kmol[k];
	  }
	Gibbs[na][j]=hold;
      }
    hold2=0;
    hold3=0;
    for(k=0; k<ns; k++)
      {
	hold2=hold2+work[0].kmol[k];
	hold3=hold3+work[0].kmol[k]*work[0].h[k]/(Rgas*work[0].temp);
      }
    Gibbs[na][na]=hold2-work[0].kmolsum;
    Gibbs[na][na+1]=hold3;
    /* row (na+1) of the matrix */
    for(j=0; j<na; j++)
      {
	hold=0;
	for(k=0; k<ns; k++)
	  {
	    hold=hold+eta[j][k]*work[0].kmol[k]*work[0].h[k]/
	      (Rgas*work[0].temp);
	  }
	Gibbs[na+1][j]=hold;
      }
    hold2=0;
    hold3=0;
    for(k=0; k<ns; k++)
      {
	hold2=hold2+work[0].kmol[k]*work[0].h[k]/(Rgas*work[0].temp);
	hold3=hold3+work[0].kmol[k]*( work[0].cp[k]/Rgas + (work[0].h[k]/
				   (Rgas*work[0].temp))*(work[0].h[k]/
				   (Rgas*work[0].temp)) );
      }
    Gibbs[na+1][na]=hold2;
    Gibbs[na+1][na+1]=hold3;
  }

void makeB(int ns, int na, struct eqmodel work[1], struct eqlist constant[1], 
	   double B[na+2], double eta[ATOMS][SPECIES],
	   double Molw[SPECIES], double a[SPECIES][2][7], double Rgas)
  {
    int i, k;
    double gkovrt[ns], hold2, hold3;
    for(i=0; i<na; i++)
      {
	hold2=0;
	for(k=0; k<ns; k++)
	  {
	    gkovrt[k]=0;
	    if(work[0].kmol[k]>1E-50)
	      {gkovrt[k]=work[0].g[k]/(Rgas*work[0].temp) + 
		 log(work[0].kmol[k]) - log(work[0].kmolsum) + 
		 log(work[0].pres);}
	    hold2=hold2 + eta[i][k]*work[0].kmol[k]*(gkovrt[k]-1);
	  }
	B[i]=constant[0].totatoms[i]+hold2;
      }
    hold2=0;
    hold3=0;
    for(k=0; k<ns; k++)
      {
	hold2=hold2 + work[0].kmol[k]*(gkovrt[k]-1);
	hold3=hold3 + (work[0].kmol[k]/(Rgas*work[0].temp))*(gkovrt[k]-1)*
	  enthalpyfind(Rgas, Molw, a, work[0].temp, k); 
      }
    B[na]=work[0].kmolsum + hold2;
    B[na+1]=constant[0].hstart/(Rgas*work[0].temp) + hold3;
  }
	
void gauss(int dim, double A[dim][dim], double B[dim]) 
  {
    /* This subroutine performs gaussian elimination with partial pivoting
       on the system Ax=B */

    int i, j, k, rowpos, colpos;
    double biggest, hold, pivotinv;
    double pivot[dim];

    for(i=0; i<dim; i++)
      {pivot[i]=0;}   

    for(i=0; i<dim; i++)
      {
	/* find the abolute largest value element in the entire matrix and 
	   its position */
	biggest=0;
	for(j=0; j<dim; j++)
	  {
	    /* prevent rows that have already been reduced to having 
	       a 1 on the diagonal from being checked again */
	    if(pivot[j]!=1) 
	      {
		for(k=0; k<dim; k++)
		  {
		    if(pivot[k]==0)
		      {
			if(fabs(A[j][k])>biggest)
			   {
			     biggest=fabs(A[j][k]);
			     rowpos=j;
			     colpos=k;
			   }
		      }
		    else
		      {
			if(pivot[k]>1)    
			  {printf("SINGULAR MATRIX type 1\n");}
		      } 
		  }
	      }
	  }
	/*-----------*/
	pivot[colpos]=pivot[colpos]+1;
	if(rowpos!=colpos)
	  {
	    for(k=0; k<dim; k++)
	      {
		hold=0;
		hold=A[rowpos][k];         /* swaps rows so that the largest 
					      element lies on the diagonal */
		A[rowpos][k]=A[colpos][k];
		A[colpos][k]=hold;
	      }
	    hold=0;
	    hold=B[rowpos];
	    B[rowpos]=B[colpos];
	    B[colpos]=hold;
	  }
	if(A[rowpos][colpos]==0)
	  {
	    /*assert(A[colpos][colpos]!=0);*/
	    printf("SINGULAR MATRIX type 2\n");
	  }
	/*-----------*/
	/* divides through by largest making diag pos = 1 */
	pivotinv=1/A[colpos][colpos];
	for(k=0; k<dim; k++)
	  {A[colpos][k]=A[colpos][k]*pivotinv;}   
	B[colpos]=B[colpos]*pivotinv;

	/*-----------*/
	/* creates zero elements above and below the pivoted row diag */
	for(j=0; j<dim; j++)
	  {
	    if(j!=colpos)
	      {
		hold=0;
		hold=A[j][colpos];
		for(k=0; k<dim; k++)
		  {A[j][k]= A[j][k] - A[colpos][k]*hold;}   
		B[j]=B[j] - B[colpos]*hold;
	      }
	  }
      }
  }

void extract(int na, int ns, struct eqmodel work[1], double B[na+2], 
	     double eta[ATOMS][SPECIES], double Molw[SPECIES], 
	     double a[SPECIES][2][7], double Rgas)
     /*This function calculates the desired variables from the 
       results of the gaussian elimination which yields the 
       B matrix */
  {
     int i, k;
     double hold2;
     double gkovrt[ns];
     work[0].dlnkmolsum=B[na];
     work[0].dlntemp=B[na+1];
     for(k=0; k<ns; k++)
       {
	 hold2=0;
	 gkovrt[k]=0;
	 work[0].dlnkmol[k]=0;
	 if(work[0].kmol[k]>1E-50)
	   {gkovrt[k]=work[0].g[k]/(Rgas*work[0].temp) + 
	      log(work[0].kmol[k]) - log(work[0].kmolsum) + 
	      log(work[0].pres);}
	 hold2=work[0].dlnkmolsum + (work[0].h[k]/(Rgas*work[0].temp))*
	   work[0].dlntemp - gkovrt[k];
	 for(i=0; i<na; i++)
	   {
	     work[0].dlnkmol[k]=work[0].dlnkmol[k] + B[i]*eta[i][k];
	   }
	 work[0].dlnkmol[k]=work[0].dlnkmol[k] + hold2;
       }
  }
	  
void update(int ns, struct eqmodel work[1])
  {
    /* calculate the relaxation parameter for numerical stability */
    int k;
    double minor, max, hold; 
    double newlnT, newlnkmolsum, newlnkmol[ns];
    double relax = 1;
    double relax1 = 1;
    double relax2 = 1;
    minor=10E-8*work[0].kmolsum;  /* minor species upper limit */
    max=fabs(work[0].dlnkmolsum);
    if(fabs(work[0].dlntemp)>max)
      {max=fabs(work[0].dlntemp);}
    for(k=0; k<ns; k++)
      {
	if(work[0].kmol[k]>minor) 
	  {
	    if(work[0].dlnkmol[k]>max)
	      {max=work[0].dlnkmol[k];}
	  }
	if(work[0].kmol[k]<minor) 
	  {
	    if(work[0].dlnkmol[k]>0)
	      {
		relax2=(log(1.0E-4)-log(work[0].kmol[k])+log(work[0].kmolsum))/
		  (work[0].dlnkmol[k] - work[0].dlnkmolsum);
	      }
	    if(fabs(relax2)<relax)
	      {relax=fabs(relax2);}
	  }
      }
    relax1=2/max;
    if(fabs(relax1)<relax)
      {relax=relax1;}
    
    /* update the new ln-ed variables */
    newlnT=log(work[0].temp) + relax*work[0].dlntemp;
    newlnkmolsum=log(work[0].kmolsum) + relax*work[0].dlnkmolsum;
    for(k=0; k<ns; k++)
      {newlnkmol[k]=log(work[0].kmol[k]) + relax*work[0].dlnkmol[k];}

    /* update the working matrix */
    work[0].temp=exp(newlnT);
    work[0].kmolsum=exp(newlnkmolsum);
    for(k=0; k<ns; k++)
      {work[0].kmol[k]=exp(newlnkmol[k]);}
  }
  
double converge(int ns, struct eqmodel work[1])
  {
    /*this calculates the convergence criteria for the equilibrium
      chemistry routine*/
    int k;
    double converge;
    double sum=0;
    double hold;
    for(k=0; k<ns; k++)
      {sum=sum+work[0].kmol[k];}
    converge=fabs(work[0].kmolsum-sum)/work[0].kmolsum;
    /*    for(k=0; k<ns; k++)
      {
	hold = fabs(work[0].kmol[k]);
	if(hold > converge)
	  {converge=hold;}
	  }*/
    return(converge);
  }

void properties(int ns, double Molw[SPECIES], struct eqmodel work[1], 
		struct eqlist constant[1], struct eqout eqfinal[1], 
		double a[SPECIES][2][7], double Rgas)
  {
    /*this calculates _all the quantities of interest and places them
      into the array eqfinal*/
    int k;
    double ratio, exponent;
    double cpsum = 0;
    double kg = 0;
    for(k=0; k<ns; k++)
      {
	work[0].cp[k]=cpfind(Rgas, Molw, a, work[0].temp, k);
	cpsum=cpsum+work[0].cp[k]*work[0].kmol[k];
	eqfinal[0].mass[k]=work[0].kmol[k]*Molw[k];
	kg=kg+eqfinal[0].mass[k];
      }
    eqfinal[0].masstot=kg;
    for(k=0; k<ns; k++)
      {eqfinal[0].massfrac[k]=eqfinal[0].mass[k]/eqfinal[0].masstot;}
    eqfinal[0].Cpmix=cpsum;
    eqfinal[0].rhotot=work[0].pres*101325/
      (Rgas/eqfinal[0].Molwmix*work[0].temp);
    eqfinal[0].Molwmix=eqfinal[0].masstot/work[0].kmolsum;
    eqfinal[0].Rgasmix=Rgas/eqfinal[0].Molwmix;
    eqfinal[0].gammamix=cpsum/(cpsum-Rgas/eqfinal[0].Molwmix);
    eqfinal[0].Tstar=work[0].temp/(1+(eqfinal[0].gammamix-1)/2);
    ratio=work[0].temp/eqfinal[0].Tstar;
    exponent=eqfinal[0].gammamix/(eqfinal[0].gammamix-1);
    eqfinal[0].pstar=(work[0].pres*101325)/pow(ratio,exponent);
    eqfinal[0].massflowout = eqfinal[0].fuelin + 
      constant[0].atof*eqfinal[0].fuelin;
    eqfinal[0].achamber=sqrt(eqfinal[0].gammamix*eqfinal[0].Rgasmix
			     *work[0].temp);
  }

void printscreen1(struct eqmodel work[1], struct eqout eqfinal[1], 
		  struct eqlist constant[1], int ns, FILE *summary)
  {
    int k;
    fprintf(summary,"The following is a summary of the pertinent variables "
	    "calculated for the\nrocket design process.  Presented first are "
	    "combustion chamber and subsequent\nisentropic expansion "
	    "variables followed by the initial design parameters for\nthe "
	    "nozzle contour which in turn are folllowed finally by the "
	    "pertinent\nvariables from the nozzle contour calculations.\n"
	    "The species map is as follows...\n\n"
	    "/*   0. H2       5. H2O       10. CO2         15. C12H24   */\n"
	    "/*   1. O2       6. HO2       11. CH3         16. N        */\n"
	    "/*   2. H        7. C2H2      12. CH4         17. NO       */\n"
	    "/*   3. O        8. C2H4      13. H2CO        18. N2       */\n"
	    "/*   4. OH       9. CO        14. C6H120                   */\n");
    fprintf(summary,"\n\nThis number should be equal to 1 ....%.8f\n\n",
	    eqfinal[0].masstot);
    if(OXIDIZER==0)
      {
	fprintf(summary,"The combustion equation is: phi*C12H24 + 18(O2 + "
		"3.76N2) --> 12H2O + 12CO2 + 67.68N2\n\n");
      }
    if(OXIDIZER==1)
      {
	fprintf(summary,"The combustion equation is: phi*C12H24 + 18O2 "
		"--> 12H2O + 12CO2\n\n");
      }
    if(OXIDIZER==2)
      {
	fprintf(summary,"The combustion equation is: phi*2H2 + O2 --> "
		"2H2O\n\n");
      }
    if(OXIDIZER==3)
      {
	fprintf(summary,"The combustion equation is: phi*2H2 + "
		"(O2 + 3.76N2) --> 2H2O + 3.76N2\n\n");
      }
    fprintf(summary, "The relevant species are (those with massfractions " 
	    "greater than 0.0001):\n");
    for(k=0;k<ns;k++)
      {
	if(eqfinal[0].massfrac[k]>1.0e-4)
	  {fprintf(summary, "%d = %.3e\n",k,eqfinal[0].massfrac[k]);}
      }
    fprintf(summary,"\n\n");
    fprintf(summary,"Equivalence Ratio = %.3f\nEquilbrium Temp = %.2f K\n"
	    "Mixture Cp = %.5f J/kmol K"
	    "\nMixture gamma = %.7f\nMixture Gas Constant = %.5f J/kmol K\n"
	    "Mixture Molecular Weight = %.3f kg/kmol\nThroat Temp = %.2f K \n"
	    "Throat Pres = %.2e Pa\nTotal Fuel Input = %.2f kg/s\n"
	    "Total Mass Flow Out = %.2f kg/s\n",
	    constant[0].phi,work[0].temp,eqfinal[0].Cpmix,eqfinal[0].gammamix,
	    eqfinal[0].Rgasmix, eqfinal[0].Molwmix, eqfinal[0].Tstar, 
	    eqfinal[0].pstar, eqfinal[0].fuelin, eqfinal[0].massflowout);
  }

void printoutput1(struct eqmodel work[1],struct eqout eqfinal[1], FILE *output, 
		 int numspec)
  {
    int p, ns;

    fprintf(output,"title=\"Equilibrium Chemistry Combustion\"\n");
    fprintf(output,"variables=\"Species\", \"frac\"\n");
    fprintf(output,"zone T=\"Massfraction\"\n");
    fprintf(output, "I=%d, J=1, K=1, F=POINT\n",numspec);
    for(p=0; p<numspec; p++)
      {fprintf(output,"%d\t%.5e\n",p, eqfinal[0].massfrac[p]);}
    fprintf(output,"\nzone T=\"Molfraction\"\n");
    fprintf(output, "I=%d, J=1, K=1, F=POINT\n",numspec);
    for(p=0; p<numspec; p++)
      {fprintf(output,"%d\t%.5e\n",p, work[0].kmol[p]/work[0].kmolsum);}      
    printf("\n");
  }  

/***************************************************************************
 **************************************************************************
 NOZZLE DESIGN FUNCTIONS
 **************************************************************************
 ***************************************************************************/
void nozinitialize(struct nozlist *set, struct eqmodel work[1], 
		   struct eqout eqfinal[1])
  {
    set->gamma=eqfinal[0].gammamix;
    set->ttot=work[0].temp;
    set->ptot=INITIAL_PRESSURE;
    set->mexit=MEXIT;
    set->hthroat=THROAT;
    set->thetawall=(PI/180)*MAXANGLE; /*convert input degrees to radians*/
    set->Rmix=eqfinal[0].Rgasmix;
    set->Cpmix=eqfinal[0].Cpmix;      /*set->Rmix*set->gamma/(set->gamma -1);*/
    set->astar=sqrt(2*set->gamma*set->Rmix*set->ttot/(set->gamma +1));
    set->rhotot=set->ptot/(set->Rmix*set->ttot);
    set->rhostar=set->rhotot*pow(2/(set->gamma+1),1/(set->gamma-1));
    set->rhostarfrac=pow(2/(set->gamma+1),1/(set->gamma-1));
    set->tstar=set->ttot*2/(set->gamma+1);
    set->pstar=set->ptot*pow(2/(set->gamma+1),set->gamma/(set->gamma-1));
    set->texit=set->ttot/(1+((set->gamma-1)/2)*MEXIT*MEXIT);
    set->pexit=set->ptot/pow((1+((set->gamma-1)/2)*MEXIT*MEXIT), 
			     set->gamma/(set->gamma-1));
    set->achamber=eqfinal[0].achamber;
    set->aexit=sqrt(set->gamma*set->Rmix*set->texit);
  }

void finddfpoly(struct nozlist set, struct nozmodel *node, struct slope *deriv,
	       struct poly *coeff)
  {
    /*this function determines the coefficients for the non-dimensional
      velocity distribution between pts F and D (which is not a radial
      flow region)*/

    double wf, wd, wdiff;
    double rd, xdf, d0, d1, d2, d3, d4;
    double wprime, wdblprime, test;
    int i;
    wprime=0;
    wdblprime=1;
    test=-0.75*wprime*wprime/wdblprime;
    wf=findw(set.gamma, set.mexit);  /*get Wf*/
    wdiff=0.015*wf;                  /*initial guess at Wf-Wd, 1.5% of Wf*/
    while(wdiff>test)
      {
	wdiff=wdiff-(5.0e-4);
	wd=wf-wdiff;
	rd=findrbar(set, wd);   
	wprime=findwprime(wd, rd, set.gamma);
	wdblprime=findwdbleprime(wd, rd, set.gamma);
	test=-0.75*wprime*wprime/wdblprime;
      }
    xdf=(-3*wprime/wdblprime)*(1-sqrt(1+4*wdiff*wdblprime/(3*wprime*wprime)));
    d0=wd;
    d1=xdf*wprime;
    d2=xdf*xdf*wdblprime/2;
    d3=10*wdiff-6*d1-3*d2;
    d4=-15*wdiff+8*d1+3*d2;

    /*place the values found in the appropriate structures*/
    node->wf=wf;
    node->wd=wd;
    node->rd=rd;
    node->xdforcr=xdf;
    coeff->D[0]=d0;
    coeff->D[1]=d1;
    coeff->D[2]=d2;
    coeff->D[3]=d3;
    coeff->D[4]=d4;
  }

void radialbounds(struct nozlist set, struct nozmodel *node)
  {
    /*this function finds the critical points within the radial flow region,
      i.e. the non-dimensional radius and velocity at points D, C, and B.
      Also, the value of Psi at point C is determined, which becomes the 
      limiting value with which one defines the nozzle contour*/
    
    double angdiff, theta, thetab, thetac, thetad;
    double w, wlo, whi, wb, wc, wd;
    double rbarc, rbarb;
    double psibarc;
    double test;

    angdiff=set.thetawall;
    wd=node->wd;
    thetad=findangle(set,wd);
    if(thetad<2*angdiff)
      {printf("\nCannot integrate the turning angle from the sonic line, "
	      "try adjusting Wdiff\n");}

    /*find the non-dim vel for thetac in the radial flow region*/
    wlo=1;
    whi=wd;
    test=100;
    while(test>1.0e-7)
      {
	w=(whi+wlo)/2;
	theta=findangle(set,w);
	if((thetad-theta)<angdiff)
	  {whi=w;}
	if((thetad-theta)>angdiff)
	  {wlo=w;}
	test=(whi-wlo)/whi;
      }
    wc=w;
    thetac=findangle(set,wc);
    rbarc=findrbar(set,wc);
    psibarc=findpsibar(set,wc,set.thetawall);
 
   /*find the non-dim vel for thetab in the radial flow region*/
    whi=wc;
    wlo=1;
    angdiff=2*angdiff;
    test=100;
    while(test>1.0e-7)
      {
	w=(whi+wlo)/2;
	theta=findangle(set,w);
	if((thetad-theta)<angdiff)
	  {whi=w;}
	if((thetad-theta)>angdiff)
	  {wlo=w;}
	test=(whi-wlo)/whi;
      }
    wb=w;
    thetab=findangle(set,wb);
    rbarb=findrbar(set,wb);

    /*insert values into the proper structures*/
    node->wc=wc;
    node->wb=wb;
    node->rc=rbarc;
    node->rb=rbarb;
    node->thetac=thetac;
    node->thetab=thetab;
    node->thetad=thetad;
    node->psibarc=psibarc;
  }

void findabpoly(struct nozlist set, struct nozmodel *node, struct slope *deriv,
	       struct poly *coeff)
  {
    /*this function determines the coefficients for the non-dimensional
      velocity distribution between pts B and A (which is not a radial
      flow region).  Also determined here is the non-dimensional distance
      of the sonic point from the throat (eps) as well as the distance
      between the origin of the radial co-ordinates and the throat (xinc).
      Other values found here are the ratio of massflows (K) between this 
      nozzle (with a curved sonic line) and a theoretical nozzle with 
      a straight sonic line at the throat and the ratio of the throat
      arc radius to the throat half height (R)*/

    double K;      /*ratio of massflows, NOT equal to k*/
    double R1;     /*ratio of Rcr/hthroat*/
    double R;      /*ratio of throat arc radius to throat radius*/
    double epsbar; /*distance from the geometric throat to sonic point*/
    double Rprev, Rstep, Rlimit, Rlowerror, Rmin, Rmax, S;
    double rb, wb, xab, chi, chiprev, errorst, errormin, dir, test, test2;
    double c0, c1, c2, c3, c4;
    
    wb=node->wb;
    rb=node->rb;	
    deriv->bprime=findwprime(wb, rb, set.gamma);
    deriv->bdprime=findwdbleprime(wb, rb, set.gamma);

    /*initial values*/
    K=1;          /*initial guess*/
    Rlimit=300;   /*maximum value*/
    Rstep=5.0e-1; /*incremental value of R*/
    Rprev=0;      /*initializing the prev R value to ensure loop is entered*/
    R=15;         /*starting value*/
    Rlowerror=R;  /*initial R value yielding the lowest error*/
    test2=fabs(R-Rprev)/(2*(R+Rprev));
    
    while(test2>1.0e-7)
      {
	R1=K/(2*sin(set.thetawall/2));  /*updated by new K each iteration*/
	S=R+1;
	deriv->aprime=findsonicprime(set, S, R1);
	deriv->adprime=findsonicdbleprime(set, S, R1);
	deriv->atprime=findsonictpleprime(set, S, R1);
	node->xaborcr=findxab(*deriv,*node);
	chi=findxaberror(*deriv,*node);
	chiprev=chi;
	errorst=chi;
	errormin=fabs(errorst);
	
	/*continuous coarse loop (crossing the root)*/
	while(chi/chiprev>0 && R<Rlimit)
	  {
	    chiprev=chi;
	    R=R+Rstep;
	    S=R+1;
	    deriv->aprime=findsonicprime(set, S, R1);
	    deriv->adprime=findsonicdbleprime(set, S, R1);
	    deriv->atprime=findsonictpleprime(set, S, R1);
	    node->xaborcr=findxab(*deriv,*node);
	    chi=findxaberror(*deriv,*node);
	    if(fabs(chi)<errormin) /*incase chi never reaches zero*/
	      {
		Rlowerror=R;
		errormin=fabs(chi);
	      }
	  }
	if(chi/chiprev<=0 && R<Rlimit)
	  {
	    dir=chiprev/fabs(chi);
	    Rmin=R-Rstep;
	    Rmax=R;
	    test=(Rmax-Rmin)/R;
	    /*continuous fine loop (homing in on the R value to yield chi=0)*/
	    while(test>1.0e-7)
	      {
		R=(Rmin+Rmax)/2;
		S=R+1;
		deriv->aprime=findsonicprime(set, S, R1);
		deriv->adprime=findsonicdbleprime(set, S, R1);
		deriv->atprime=findsonictpleprime(set, S, R1);
		node->xaborcr=findxab(*deriv,*node);
		chi=findxaberror(*deriv,*node);
		if(dir*chi>0)
		  {Rmin=R;}
		if(dir*chi<=0)
		  {Rmax=R;}
		test=(Rmax-Rmin)/R;
	      }
	  }
	else
	  {
	    R=Rlowerror;
	    S=R+1;
	    deriv->aprime=findsonicprime(set, S, R1);
	    deriv->adprime=findsonicdbleprime(set, S, R1);
	    deriv->atprime=findsonictpleprime(set, S, R1);
	    node->xaborcr=findxab(*deriv,*node);
	    chi=findxaberror(*deriv,*node);
	  }
	test2=fabs(R-Rprev)/(2*(R+Rprev));
	K=findmassflow(set.gamma,S);
	Rprev=R;
      }
    xab=node->xaborcr;
    c0=1;
    c1=xab*deriv->aprime;
    c2=(xab*xab*deriv->adprime)/2;
    c3=4*(wb-1)-xab*deriv->bprime-3*c1-2*c2;
    c4=-3*(wb-1)+xab*deriv->bprime+2*c1+c2;
    R1=K/(2*sin(set.thetawall/2));
    epsbar=findsonicpt(set,S,R1);

    /*place the values found in the appropriate structures*/
    node->R=R;
    node->K=K;
    node->epsbar=epsbar;
    node->xincorcr=node->rb-xab-epsbar; /*this value lies upstream of throat*/
    coeff->C[0]=c0;
    coeff->C[1]=c1;
    coeff->C[2]=c2;
    coeff->C[3]=c3;
    coeff->C[4]=c4;
  }

void dfveldist(struct nozmodel *node, struct poly coeff)
  {
    /*this function initializes the non-dimensional velocity and position 
      arrays along the centreline between points F and D by using the
      appropriate polynomial and the set number of points IFTOD*/

    int i, st, end;
    double pos, zeta, d0, d1, d2, d3, d4;
    
    st=0;
    end=IFTOD;
    d0=coeff.D[0];
    d1=coeff.D[1];
    d2=coeff.D[2];
    d3=coeff.D[3];
    d4=coeff.D[4];
    for(i=st; i<end; i++) /*this produces pt F, but not D*/
      {
	pos=i;
	zeta=1-pos/IFTOD;
	node->w[i]=d0 + d1*zeta + d2*zeta*zeta + d3*zeta*zeta*zeta
	  + d4*zeta*zeta*zeta*zeta;
	node->rorcr[i]=zeta*node->xdforcr + node->rd; 
      }
  }

void radialveldist(struct nozlist set, struct nozmodel *node)
  {
    /*this function initializes the non-dimensional velocity and position 
      arrays along the centreline between points D and B by using the 
      properties of a radially flowing perfect gas and the desired number
      of points within this region, IDTOC and ICTOB*/

    int i, st, end;
    double wlo, whi;
    
    /*for the region between points D and C going towards the throat*/
    st=IFTOD;
    end=IFTOD+IDTOC;
    wlo=node->wc;
    whi=node->wd;
    for(i=(st+1); i<end; i++) /*this does not produce the boundary values*/
      {
	node->w[i]=whi-(whi-wlo)*(i-st)/IDTOC;
	node->rorcr[i]=findrbar(set,node->w[i]);
      }
    /*boundary values*/
    node->w[st]=node->wd;
    node->w[end]=node->wc;
    node->rorcr[st]=node->rd;
    node->rorcr[end]=node->rc;

    /*for the region between points C and B going towards the throat*/
    st=IFTOD+IDTOC;
    end=IFTOD+IDTOC+ICTOB;
    wlo=node->wb;
    whi=node->wc;
    for(i=(st+1); i<end; i++) /*this does not produce the boundary values*/
      {
	node->w[i]=whi-(whi-wlo)*(i-st)/ICTOB;
	node->rorcr[i]=findrbar(set,node->w[i]);
      }
    /*boundary values*/
    node->w[end]=node->wb;
    node->rorcr[end]=node->rb;
  }

void abveldist(struct nozmodel *node, struct poly coeff)
  {
    /*this function initializes the non-dimensional velocity and position 
      arrays along the centreline between points B and A by using the
      appropriate polynomial and the set number of points IBTOA*/

    int i, st, end;
    double pos, eta, c0, c1, c2, c3, c4;
    
    st=IFTOD+IDTOC+ICTOB;
    end=IFTOD+IDTOC+ICTOB+IBTOA;
    c0=coeff.C[0];
    c1=coeff.C[1];
    c2=coeff.C[2];
    c3=coeff.C[3];
    c4=coeff.C[4];
    for(i=(st+1); i<(end+1); i++) /*this produces pt A, but not B*/
      {
	pos=i;
	eta=(1-(pos-st)/IBTOA)*(1-(pos-st)/IBTOA);
	node->w[i]=c0 + c1*eta + c2*eta*eta + c3*eta*eta*eta
	  + c4*eta*eta*eta*eta;
	node->rorcr[i]=node->xincorcr + node->epsbar + (eta)*node->xaborcr;
      }
  }

void findcontour(struct nozlist set, struct nozmodel node, struct nozout *final)
  {
    int i, j;
    int ist, iend, jend, jlimit, grow, hit;
    double k, gamma, xinc, mf, af, rf;
    double psibar[MAXJARRAY], delpsibar[MAXJARRAY];
    double wa, wb, rhofraca, rhofracb, anga, angb, tempa, tempb;
    double xb, ub, vb, qa, qb, thetaa, thetab, percent, arcrad, Save;
    double x1[MAXJARRAY], y1[MAXJARRAY], u1[MAXJARRAY], V1[MAXJARRAY];
    double a1[MAXJARRAY], m1[MAXJARRAY];
    double x2[MAXJARRAY], y2[MAXJARRAY], u2[MAXJARRAY], v2[MAXJARRAY];
    double a2[MAXJARRAY], m2[MAXJARRAY];
    double dy, yb, dpsib;

    struct moc mocvars;

    ist=0;
    arcrad=0;
    iend=IFTOD+IDTOC+ICTOB+IBTOA;
    jend=JTOT+6; /*starting value only, gets updated as mesh is built, but
		  is greater than JTOT to start incase psiwall is not reached 
		  until slightly past hexitbar due to the approximation used 
		  in finding the incremental stream function*/
    gamma=set.gamma;
    xinc=node.xincorcr;
    k=(gamma+1)/(gamma-1);
    dy=node.hexitbar/JTOT;

    /*initialize final characteristic line eminating from pt F*/
    mf=findm(set,node.wf);
    af=finda(set,node.wf);
    rf=node.rorcr[0];
    for(j=0; j<jend; j++)
      {
	y1[j]=0+j*dy;
	x1[j]=(rf + y1[j]*sqrt(mf*mf-1)) - xinc;
	u1[j]=node.wf*set.astar;
	V1[j]=0;
	a1[j]=af;
	m1[j]=mf;
      }
    /*record the nozzle exit parameters*/
    final->yworcr[0]=node.hexitbar;
    final->qw[0]=node.wf*set.astar;
    final->thetaw[0]=0;
    final->xworcr[0]=(rf + node.hexitbar*sqrt(mf*mf-1)) - xinc;
   
    for(i=ist; i<iend; i++)
      {
	hit=0;
	for(j=0; j<jend; j++)
	  {
	    if(j==0)
	      {
		x2[j]=node.rorcr[i+1] - xinc;
		y2[j]=0;
		u2[j]=node.w[i+1]*set.astar;
		v2[j]=0;
		a2[j]=finda(set,node.w[i+1]);
		m2[j]=findm(set,node.w[i+1]);
		psibar[j]=0;
		yb=0;
		dpsib=0;
	      }
	    else
	      {
		x2[j]=mocvars.x4;
		y2[j]=mocvars.y4;
		u2[j]=mocvars.u4;
		v2[j]=mocvars.v4;
		a2[j]=mocvars.a4;
		m2[j]=mocvars.m4;

		/*calculate stream function to see if upper wall is reached*/
		qa=sqrt(u2[j-1]*u2[j-1]+v2[j-1]*v2[j-1]);
		qb=sqrt(u2[j]*u2[j]+v2[j]*v2[j]);
		wa=qa/set.astar;
		wb=qb/set.astar;
		rhofraca=findrhofrac(set,wa);
		rhofracb=findrhofrac(set,wb);
		thetaa=atan(v2[j-1]/u2[j-1]);
		thetab=atan(v2[j]/u2[j]);
		anga=asin(1/m2[j-1])+thetaa;
		angb=asin(1/m2[j])+thetab;
		tempa=rhofraca*wa/(4*m2[j-1]*sin(anga));
		tempb=rhofracb*wb/(4*m2[j]*sin(angb));
		delpsibar[j]=(tempa+tempb)*(y2[j]*y2[j]-y2[j-1]*y2[j-1]);
		psibar[j]=psibar[j-1]+delpsibar[j];
	      }

	    if(psibar[j]>=node.psibarc && hit==0)
	      {
		percent=(node.psibarc-psibar[j-1])/(psibar[j]-psibar[j-1]);
		final->yworcr[i+1]=sqrt(y2[j-1]*y2[j-1] + 
				       percent*(y2[j]*y2[j]-y2[j-1]*y2[j-1]));
		final->qw[i+1]=qa+percent*(qb-qa);
		final->thetaw[i+1]=thetaa+percent*(thetab-thetaa);
		final->xworcr[i+1]=x2[j-1] + (x2[j]-x2[j-1])*
		  (final->yworcr[i+1]-y2[j-1])/(y2[j]-y2[j-1]);
		jlimit=j;
		hit=1;
	      }
	    findmocvalues(u1[j], V1[j], a1[j], m1[j], y1[j], x1[j], 
			  u2[j], v2[j], a2[j], m2[j], y2[j], x2[j], 
			  j, set, &mocvars);	
	  }/*end of j loop*/

	printf("i = %d\tY = %.9f\tX = %.5f\n",i, 
	       final->yworcr[i], final->xworcr[i]);

	for(j=0; j<(jend+1); j++)
	  {
	    y1[j]=y2[j];
	    x1[j]=x2[j];
	    u1[j]=u2[j];
	    V1[j]=v2[j];
	    a1[j]=a2[j];
	    m1[j]=m2[j];
	  }
	jend=jlimit+10; /*allows for the characteristic mesh to vary in size
			 in the y direction*/

	if(i>=(iend-6) && i<(iend-1)) /*calculate the ave. R value based on
					the wall co-ordinates near the
					throat*/
	  {
	    arcrad=arcrad+ (pow(final->xworcr[i+1],2) +
			    pow(final->yworcr[i+1]-node.hthroatbar,2))/
	      (2*(final->yworcr[i+1]-node.hthroatbar));
	  }
      }/*end of i loop*/

    /*record the nozzle input parameters*/
    final->yworcr[IFTOD+IDTOC+ICTOB+IBTOA]=node.hthroatbar;
    final->Rave=arcrad/(5*node.hthroatbar);
    Save=final->Rave+1;
    final->qw[IFTOD+IDTOC+ICTOB+IBTOA]=(1+1/(4*Save)-(14*gamma+15)/
	     (288*Save*Save)+(2364*gamma*gamma - 3915*gamma + 14337)/
             (82944*Save*Save*Save))*set.astar;
    final->thetaw[IFTOD+IDTOC+ICTOB+IBTOA]=0;
    final->xworcr[IFTOD+IDTOC+ICTOB+IBTOA]=0;

  }

void findmocvalues(double u1in, double v1in, double a1in, double m1in, 
		   double y1in, double x1in,
		   double u2in, double v2in, double a2in, double m2in, 
		   double y2in, double x2in, int j, struct nozlist set, 
		   struct moc *mocvars)
  {
    int iter; 
    double alphap, alpham, thetam, thetap, lambdam, lambdap;
    double Um, Up, Rm, Rp, Sm, Sp, Tp, Tm;
    double w1, w2; 
    double x1, y1, u1, V1, a1, m1, x2, y2, u2, v2, a2, m2;
    double x4, y4, u4, v4, w4;
    
    iter=0;
    x1=x1in; x2=x2in; y1=y1in; y2=y2in; u1=u1in; u2=u2in;
    V1=v1in; v2=v2in; a1=a1in; a2=a2in; m1=m1in; m2=m2in;

    while(iter<5) /*sets the number of times through the pred/corr process*/
      {

	/*minus*/
	alpham=asin(1/m1);
	thetam=atan(V1/u1);
	lambdam=tan(thetam - alpham);
	Um=u1*u1 - a1*a1;
	Rm=2*u1*V1 - Um*lambdam;
	if(j==0)
	  {Sm=0;}
	else
	  {Sm=a1*a1*V1/y1;}
	
	/*plus*/
	alphap=asin(1/m2);
	thetap=atan(v2/u2);
	lambdap=tan(thetap + alphap);
	Up=u2*u2 - a2*a2;
	Rp=2*u2*v2 - Up*lambdap;
	if(j==0)
	  {Sp=0;}
	else
	  {Sp=a2*a2*v2/y2;}

	/*solve for pt 4*/
	x4=(y1in - y2in - lambdam*x1in + lambdap*x2in)/(lambdap - lambdam);
	y4=y1in+lambdam*(x4-x1in);
	Tm=Sm*(x4-x1)+Um*u1+Rm*V1;
	Tp=Sp*(x4-x2)+Up*u2+Rp*v2;
	u4=(Tm*Rp - Tp*Rm)/(Um*Rp - Up*Rm);
	v4=(Tp*Um - Tm*Up)/(Um*Rp - Up*Rm);

	/*calculate the corrector values*/
	y1=(y1in+y4)/2;
	u1=(u1in+u4)/2;
	V1=(v1in+v4)/2;
	w1=sqrt(u1*u1+V1*V1)/set.astar;
	a1=finda(set,w1);
	m1=findm(set,w1);
	
	y2=(y2in+y4)/2;
	u2=(u2in+u4)/2;
	v2=(v2in+v4)/2;
	w2=sqrt(u2*u2+v2*v2)/set.astar;
	a2=finda(set,w2);
	m2=findm(set,w2);
	iter=iter+1;
      }

    w4=sqrt(u4*u4+v4*v4)/set.astar;
    mocvars->x4=x4;
    mocvars->y4=y4;
    mocvars->u4=u4;
    mocvars->v4=v4;
    mocvars->a4=finda(set,w4);
    mocvars->m4=findm(set,w4);
  }

void boundarylayer(struct nozlist set, struct nozout *final,
		   struct eqmodel work[1], double Molw[SPECIES], 
		   double a[SPECIES][2][7], double epsilonk[SPECIES], 
		   double sigmak[SPECIES], double omega11polyd[3][5], 
		   double astarpolye[3][3], struct eqout eqfinal[1], 
		   double Rgas, int ns)
  {
    /*this function calculates the required adjustment to the inviscid
      nozzle contour to account for the development of a turbulent
      boundary layer (based on Edenfield's analysis)*/
    
    int i, k, st, end;
    double gamma, cp, cpwall, cptemp, R, Twall, Ttot, ptot;
    double thetawall, x, scale, cphold;
    double muinv, hinv, rhoinv, pinv, Tinv, qinv, Minv, Reyinv;
    double muref, href, rhoref, Tref, Reyref;
    double T, Thigh, Tlow, h;
    double htot, hwall, haw, deltaphys, deltablayer;

    st=0;
    end=IFTOD+IDTOC+ICTOB+IBTOA;
    Twall=TWALL;
    scale=1;/*set.hthroat/(final->yworcr[end-1]*final->rcr);*/
    gamma=set.gamma;
    cp=set.Cpmix;
    R=set.Rmix;
    Ttot=set.ttot;
    ptot=set.ptot;
    htot=cp*Ttot;
    printf("Correcting inviscid contour for boundary layer");
    for(i=st;i<end;i++) /*will not calculate a displacement for the throat*/
      {
	thetawall=final->thetaw[i];
	x=final->xworcr[i]*final->rcr;
	/*calculate the conditions at the edge of the boundary layer
	  (which are equal to the conditions along the inviscid contour)*/
	qinv=final->qw[i];
	hinv=htot-0.5*qinv*qinv;
	Tinv=Ttot - qinv*qinv/(2*cp);
	Minv=qinv/sqrt(gamma*R*Tinv);
	pinv=ptot*pow((Tinv/Ttot),(gamma/(gamma-1)));
	rhoinv=pinv/(R*Tinv);
	muinv=findmu(epsilonk, sigmak, omega11polyd, astarpolye, Molw, 
		      eqfinal, Tinv, ns);
	Reyinv=(rhoinv*qinv*x/muinv)*scale;
	final->Tinv[i]=Tinv;
	final->Reyinv[i]=Reyinv;
	/*calculate the reference conditions NOTE: do not use
	 the NASA polynomials for enthalpy as these are values
	 with respect to a specific reference condition, which
	 is not the case for the value of href, while cp is the 
	 CHANGE in enthalpy and hence the reference condition is
	 irrelevant*/
	cpwall=0;
	for(k=0; k<ns; k++)
	  {
	    cptemp=cpfind(Rgas, Molw, a, Twall, k);
	    cpwall=cpwall+work[0].kmol[k]*cptemp;
	  }
	hwall=cpwall*Twall;
	haw=hinv + 0.9*(htot-hinv);
	href=0.5*(hwall+hinv) + 0.22*(haw-hinv);
	Thigh=6000;
	Tlow=300;
	while(((Thigh-Tlow)/Thigh)>1.0e-6)
	  {
	    cphold=0;
	    for(k=0; k<ns; k++)
	      {
		T=0.5e0*(Thigh+Tlow);
		cptemp=cpfind(Rgas, Molw, a, T, k);
		cphold=cphold+work[0].kmol[k]*cptemp;
	      }
	    h=cphold*T;
	    if(h>href)
	      {Thigh=T;}
	    if(h<href)
	      {Tlow=T;}
	  }
	Tref=T;
	rhoref=pinv/(R*Tref);
	muref=findmu(epsilonk, sigmak, omega11polyd, astarpolye, Molw, 
		     eqfinal, Tref, ns);
	Reyref=(rhoref*qinv*x/muref)*scale;
	final->Tref[i]=Tref;
	final->Reyref[i]=Reyref;
	/*calculate the new nozzle contour and boundary layer position*/
	deltaphys=(0.42*x/pow(Reyref,0.2775))*scale;
	deltablayer=(1.1*x*(0.195*pow(Minv,0.375)/pow(Reyinv,0.166)))*scale;
	final->xw[i]=x*scale;
	final->yw[i]=final->yworcr[i]*final->rcr + deltaphys/cos(thetawall);
	final->ybound[i]=final->yw[i] - deltablayer/cos(thetawall);
	printf(".");
      }
    final->xw[end]=0;
    final->yw[end]=set.hthroat;
    final->ybound[end]=final->yw[end];
    final->Tinv[end]=set.tstar;
    final->Reyinv[end]=0;
    printf("done\n");
  }

void findsubcontour(struct nozlist set, struct nozout *final)
  {
    /*calculates a subsonic contoured section for the nozzle*/
    int i;
    double points, length, M, step, x, y;
    double a, b, temp, exponent, h, hchamb, gamma;

    points=ISUB;
    length=SUBLENGTH*final->xworcr[0]*final->rcr;
    M=CHAMBERM;
    step=length/(points-1);
    h=set.hthroat;
    gamma=set.gamma;
    exponent=(gamma+1)/(2*(gamma-1));
    temp=1/M*pow(( 2/(gamma+1)*(1 + ((gamma-1)/2)*M*M) ),exponent);
    hchamb=sqrt(h*h*temp);
    final->hchamb=hchamb;
    for(i=0;i<points;i++)
      {
	x=length-step*i;
	a=(x*x)/(length*length);
	b=(h/hchamb)*(h/hchamb);
	y=h/sqrt( 1-(1-b)*((1-a)*(1-a)/(1+(a/3)*(a/3)*(a/3))) );
	final->xsub[i]=-(i*step);
	final->ysub[i]=y;
      }
    printf("Subsonic contour complete.\n");
  }

double finda(struct nozlist set, double w)
  {
    /*finds the speed of sound for a perfect gas*/
    double gamma, R, Cp, To, T, a, q, M;
    gamma=set.gamma;
    R=set.Rmix;
    Cp=set.Cpmix;
    To=set.ttot;
    q=w*set.astar;
    T=To-q*q/(2*Cp);
    a=sqrt(gamma*R*T);
    return(a);
  }

double findm(struct nozlist set, double w)
  {
    /*finds the Mach number for a perfect gas*/
    double a, q, M;
    q=w*set.astar;
    a=finda(set, w);
    M=q/a;
    return(M);
  }

double findrhofrac(struct nozlist set, double w) 
  {
    /*finds the ratio of density to total density for a perfect gas*/
    double rhofrac, gamma, m;
    m=findm(set, w);
    gamma=set.gamma;
    rhofrac=pow(1+((gamma-1)/2)*m*m,(-1/(gamma-1)));
    return(rhofrac);
  }

double findpsibar(struct nozlist set, double w, double theta)
  {
    /*finds the non-dimensional stream function in a radial flow region
      integrated from the centreline (theta=0) to theta*/
    double rhofrac, rbar, psibar;
    rhofrac=findrhofrac(set,w);
    rbar=findrbar(set,w);
    psibar=rhofrac*w*rbar*rbar*(1-cos(theta));
    return(psibar);
  }

double findmassflow(double gamma, double S)
  {
    /*finds the ratio of massflows between a nozzle throat with a 
      curved sonic line and a straight one (Kliegel)*/
    double a, b, c, K;
    a=(gamma+1)/(96*S*S);
    b=(8*gamma-27)/(24*S);
    c=(754*gamma*gamma-757*gamma+3615)/(2880*S*S);
    K=sqrt(1-a*(1-b+c));
    return(K);
  }

double findhexitbar(struct nozlist set, struct nozmodel node)
  {
    /*finds the non-dimensional exit radius*/
    double wf, rhofracf, psiwall, hexitbar;
    wf=node.wf;
    psiwall=node.psibarc;
    rhofracf=findrhofrac(set,wf);
    hexitbar=sqrt((2*psiwall)/(rhofracf*wf));
    return(hexitbar);
  }

double findrcr(struct nozlist set, struct nozmodel node)
  {
    /*finds the radius at which point the flow is sonic in the 
      radial flow region*/
    double wf, rhofracf, psiwall, K, hexitbar, hbar, rcr;
    wf=node.wf;
    psiwall=node.psibarc;
    K=node.K;
    hexitbar=node.hexitbar;
    rhofracf=findrhofrac(set,wf);
    hbar=sqrt((rhofracf/set.rhostarfrac)*wf)*hexitbar/K;
    rcr=set.hthroat/hbar;
    return(rcr);
  }

double findw(double gamma, double mach)
  { 
    /*finds the non-dimensional velocity in the radial flow region
      based on the Mach number*/
    double k, w;
    k=(gamma+1)/(gamma-1);
    w=sqrt(k*mach*mach/(k-1+mach*mach));
    return(w);
  }

double findwprime(double w, double rbar, double gamma)
  {
    /*valid in the radial flow region only*/
    double k, wprime;
    k=(gamma+1)/(gamma-1);
    wprime=2*w*(k-w*w)/(k*(w*w-1)*rbar);
    return(wprime);
  }

double findwdbleprime(double w, double rbar, double gamma)
  {
    /*valid in the radial flow region only*/
    double k, wprime, wdblprime;
    k=(gamma+1)/(gamma-1);
    wprime=2*w*(k-w*w)/(k*(w*w-1)*rbar);
    wdblprime=-wprime*wprime*(3*k-6*w*w+(k+2)*w*w*w*w)/(2*w*(k-w*w)*(w*w-1));
    return(wdblprime);
  }

double findsonicprime(struct nozlist set, double S, double R1)
  {
    /*valid at the sonic point only (Kliegel)*/
    double gamma, lambda, wprime;
    double a, b;
    gamma=set.gamma;
    lambda=sqrt(2/((gamma+1)*S));
    a=(4*gamma-3)/(24*S);
    b=(652*gamma*gamma+15*gamma+333)/(6912*S*S);
    wprime=lambda*R1*(1-a+b);
    return(wprime);
  }

double findsonicdbleprime(struct nozlist set, double S, double R1)
  {
    /*valid at the sonic point only (Kliegel)*/
    double gamma, lambda, wdbleprime;
    double a, b;
    gamma=set.gamma;
    lambda=sqrt(2/((gamma+1)*S));
    a=2*gamma/3;
    b=(4*gamma*gamma-69*gamma+15)/(96*S);
    wdbleprime=(lambda*R1)*(lambda*R1)*(1-a+b);
    return(wdbleprime);
  }

double findsonictpleprime(struct nozlist set, double S, double R1)
  {
    /*valid at the sonic point only (Kliegel)*/
    double gamma, lambda, wtpleprime;
    double a, b;
    gamma=set.gamma;
    lambda=sqrt(2/((gamma+1)*S));
    a=(4*gamma*gamma-57*gamma+27)/24;
    wtpleprime=(lambda*R1)*(lambda*R1)*(lambda*R1)*a;
    return(wtpleprime);
  }

double findsonicpt(struct nozlist set, double S, double R1)
  {
    /*the non-dimensional distance between the throat and the point on 
      the centreline at which the flow is sonic (Kliegel)*/
    double gamma, lambda, eps, a, b, c;
    gamma=set.gamma;
    lambda=sqrt(2/((gamma+1)*S));
    a=1/(4*lambda*S);
    b=(4*gamma-15)/(72*S);
    c=(412*gamma*gamma+270*gamma+909)/(10368*S*S);
    eps=(a/R1)*(1-b+c);
    return(eps);
  }

double findrbar(struct nozlist set, double w)
  {
    /*valid in the radial flow region only*/
    double k, q, gamma, astar, m, rbar;
    gamma=set.gamma;
    astar=set.astar;
    k=(gamma+1)/(gamma-1);
    m=findm(set, w);
    rbar=sqrt(pow((k-1+m*m)/k,k/2)/m);
    return(rbar);
  }

double findangle(struct nozlist set, double w)
  {
    /*valid in the radial flow region only*/
    double k, gamma, m, theta;
    gamma=set.gamma;
    k=(gamma+1)/(gamma-1);
    m=findm(set, w); 
    theta=sqrt(k)*atan(sqrt((m*m-1)/k)) - atan(sqrt(m*m-1));
    theta=theta/2;
    return(theta);
  }

double findxab(struct slope deriv, struct nozmodel node)
  {
    /*finds the non-dimensional distance between the sonic point
      and the beginning of the radial flow region, along the centreline*/
    double a, b, c, xab;
    a=deriv.adprime-deriv.bdprime;
    b=6*(deriv.aprime+deriv.bprime);
    c=-12*(node.wb-1);
    xab=(-b+sqrt(b*b-4*a*c))/(2*a);
    return(xab);
  }

double findxaberror(struct slope deriv, struct nozmodel node)
  {
    /*finds the error given a value of Xab, with the correct value
      of Xab this should tend to zero*/
    double a, b, c, d, xab, error;
    xab=node.xaborcr;
    a=deriv.atprime*xab*xab*xab;
    b=6*deriv.adprime*xab*xab;
    c=6*(3*deriv.aprime+deriv.bprime)*xab;
    d=-24*(node.wb-1);
    error=a+b+c+d;
    return(error);
  }

double findomega11(double omega11polyd[3][5], double epsilonk[SPECIES], 
		   double T, int species)
  {
    /*finds the co-efficient Omega_11 for use in finding mu*/
    int index;
    double Tstar, eps, omega11;
    
    eps=epsilonk[species];
    Tstar=T/eps;
    index=2;
    if (Tstar<10.0e0) index=1;
    if (Tstar<5.0e0) index=0;
    omega11=omega11polyd[index][0]
      +omega11polyd[index][1]*Tstar
      +omega11polyd[index][2]*Tstar*Tstar
      +omega11polyd[index][3]*Tstar*Tstar*Tstar
      +omega11polyd[index][4]*Tstar*Tstar*Tstar*Tstar;
    return(omega11);
  }

double findastarformu(double astarpolye[3][3], double epsilonk[SPECIES],
		      double T, int species)
  {
    /*finds the co-efficient A^* for use in finding mu*/
    int index;
    double Tstar, eps, astar;
   
    eps=epsilonk[species];
    Tstar=T/eps;
    index=2;
    if (Tstar<10.0e0) index=1;
    if (Tstar<5.0e0) index=0;
    astar=astarpolye[index][0]
      +astarpolye[index][1]*Tstar
      +astarpolye[index][2]*Tstar*Tstar;
    return(astar);
  }

double findmu(double epsilonk[SPECIES], double sigmak[SPECIES],
	      double omega11polyd[3][5], double astarpolye[3][3],
	      double Molw[SPECIES], struct eqout eqfinal[1], double T, int ns)
  { 
    /*finds the coefficient of viscosity at the given temperature
     NOTE: Molw must be in kg/gmol not kg/kgmol thus the 1.0e-3 factor in 
     calculating mu[k]; star values are NOT sonic conditions*/
    int k, l;
    double minor, tstark, astar, omega11, omega22;
    double chi[SPECIES], mu[SPECIES];
    double sum, num, denom, phikl, sum2, sum3;

    minor=1.0e-4;    /*sets the cutoff for the minimum species massfraction*/
    sum=0;
    for(k=0;k<ns;k++)
      {sum=sum+eqfinal[0].massfrac[k]/Molw[k];} 

    sum3=0;
    for(k=0;k<ns;k++)
      {
	if(eqfinal[0].massfrac[k]>minor)
	  {
	    astar=findastarformu(astarpolye, epsilonk, T, k);
	    omega11=findomega11(omega11polyd, epsilonk, T, k);
	    omega22=astar*omega11;
	    mu[k]=8.44107e-7*(sqrt(Molw[k]*(1.0e-3)*T))/(sigmak[k]*
			      sigmak[k]*omega22);
	    chi[k]=eqfinal[0].massfrac[k]/Molw[k]/sum;
	  }
      }

    for(k=0;k<ns;k++)
      {
	if(eqfinal[0].massfrac[k]>minor)
	  {
	    sum2=0;
	    for(l=0;l<ns;l++)
	      {
		if(eqfinal[0].massfrac[l]>minor)
		  {
		    if(k!=l)
		      {
			num=1+sqrt(mu[k]/mu[l])*pow(Molw[l]/Molw[k],0.25);
			denom=sqrt(8*(1+Molw[k]/Molw[l]));
			phikl=num*num/denom;
			sum2=sum2+chi[l]*phikl;
		      }
		  }
	      }
	    sum3=sum3+mu[k]*chi[k]/(chi[k]+sum2);
	  }
      }
    return(sum3);
  }

void screenoutput(struct nozlist set, struct nozout final, FILE *summary)
  {
    fprintf(summary,"\nThe generated nozzle has the following "
	    "characteristics:\n\n"
	    "Exit Mach Number: %.2f\n"
	    "Maximum Expansion Angle: %.8f rad, %.2f degrees\n"
	    "Chamber Temperature: %.4f K\n"
	    "Throat Temperature: %.4f K\n"
	    "Exit Temperature: %.4f K\n\n"
	    "Chamber Pressure: %.4e Pa\n"
	    "Throat Pressure: %.4e Pa\n"
	    "Exit Pressure: %.4e Pa\n\n"
	    "Chamber Sound Speed: %.6f m/s\n"
	    "Throat Sound Speed: %.6f m/s\n"
	    "Exit Sound Speed: %.6f m/s\n\n"
	    "Mixture Specific Heat Ratio: %.6f\n"
	    "Mixture Gas Constant: %.6f J/kmol K\n"
	    "Sonic Density Ratio (at throat): %.6f\n\n"
	    "Chamber Exit Height: %.6f m\n"
	    "Nozzle Throat Height: %.6f m\n"
	    "Nozzle Exit Height: %.6f m\n\n"
	    "The pertinent calculated nozzle variables are...",
	    set.mexit,set.thetawall,MAXANGLE,set.ttot,set.tstar,set.texit,
	    set.ptot,set.pstar,set.pexit,set.achamber,set.astar,set.aexit,
	    set.gamma,set.Rmix,set.rhostarfrac, final.hchamb,
	    set.hthroat,final.hexit);
  }

void screenoutput2(struct nozlist set, struct nozmodel node, 
		   struct nozout final, struct poly coeff, FILE *summary)
  {
    fprintf(summary,"\n\nTheta_D = %.7f\nTheta_C = %.7f\nTheta_B = %.7f\n"
	    "W_F = %.5f\nW_D = %.5f\nW_C = %.5f\nW_B = %.5f\n"
	    "rbar_F = %.5f\nrbar_D = %.5f\nrbar_C = %.5f\nrbar_B = %.5f\n"
	    "Twice Max Expansion Angle = %.7f [rad] = %.2f [degrees]\n"
	    "Psiwallbar = %.7f\n"
	    "D_0 = %.4f\tD_1 = %.4f\tD_2 = %.4f\tD_3 = %.4f\tD_4 = %.4f\n"
	    "C_0 = %.4f\tC_1 = %.4f\tC_2 = %.4f\tC_3 = %.4f\tC_4 = %.4f\n"
	    "X_df = %.5f\nX_ab = %.5f\nX_inc = %.5f\nEpsilonbar = %.7f\n"
	    "Hexitbar = %.6f\nHthroatbar = %.6f\nR_cr = %.6f [m]\n"
	    "Hexit = %.6f [m]\nHthroat = %.6f [m]\nR = %.5f\t"
	    "R_ave = %.5f (These values should be close)\n",
	    node.thetad,node.thetac,node.thetab,node.wf,node.wd,node.wc,
	    node.wb,node.rorcr[0],node.rd,node.rc,node.rb,2*set.thetawall,
	    2*set.thetawall*(180/PI),node.psibarc,
	    coeff.D[0],coeff.D[1],coeff.D[2],coeff.D[3],coeff.D[4],
	    coeff.C[0],coeff.C[1],coeff.C[2],coeff.C[3],coeff.C[4],
	    node.xdforcr,node.xaborcr,node.xincorcr,node.epsbar,
	    node.hexitbar,node.hthroatbar,final.rcr,
	    final.hexit,set.hthroat,node.R,final.Rave);
  }

void finaloutput(FILE *output, struct nozout final)
  {
    int i;
    double itot;
    itot=IFTOD+IDTOC+ICTOB+IBTOA;
    fprintf(output,"title=\"Nozzle Upper Surface Co-ordinates\""
	    "\nvariables=\"X[m]\",\"Y[m]\",\"Xinv[m]\",\"Yinv[m]\","
	    "\"Ybound[m]\",\"Xbar\",\"Ybar\",\"Qwall[m/s]\","
	    "\"Theta_wall[rad]\",\"Tinv\",\"Rey#inv\",\"Tref\",\"Rey#ref\""
	    "\nZONE T=\"r/r_cr = %.7f\"\n",final.rcr);
    fprintf(output, "I=%.0f, J=1, K=1, F=POINT\n",itot+1);
    for(i=0; i<(itot+1); i++)
      {fprintf(output,"%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t"
	       "%.7f\t%.7f\t%.7f\t%.7f\n",
	       final.xw[i], final.yw[i],final.xworcr[i]*final.rcr,
	       final.yworcr[i]*final.rcr,final.ybound[i],final.xworcr[i],
	       final.yworcr[i],final.qw[i],final.thetaw[i],
	       final.Tinv[i], final.Reyinv[i],final.Tref[i], final.Reyref[i]);}
  }

void subsonicoutput(FILE *output, struct nozout final)
  {
    int i;
    double itot;
    itot=ISUB;
    fprintf(output,"title=\"Nozzle Upper Surface Co-ordinates (Subsonic)\""
	    "\nvariables=\"Xsub[m]\",\"Ysub[m]\"\n"
	    "\nZONE T=\"subsonic\"\n");
    fprintf(output, "I=%.0f, J=1, K=1, F=POINT\n",itot);
    for(i=0; i<ISUB; i++)
      {fprintf(output,"%.7f\t%.7f\n",final.xsub[i], final.ysub[i]);}
  }

void velocityoutput(FILE *output, struct nozlist set, struct nozmodel node,
		    struct nozout final)
  {
    int i;
    double itot;
    double a, mach;
    itot=IFTOD+IDTOC+ICTOB+IBTOA;
    fprintf(output,"title=\"Centreline Velocity Distribution\""
	    "\nvariables=\"R/R_cr\",\"R[m]\",\"W\",\"MachTheo\"\n"
	    "ZONE T=\"Centreline\"\n");
    fprintf(output, "I=%.0f, J=1, K=1, F=POINT\n",itot+1);
    for(i=0; i<(itot+1); i++)
      {
	a=finda(set,node.w[i]);
	mach=node.w[i]*set.astar/a;
	fprintf(output,"%.7f\t%.7f\t%.7f\t%.7f\n",node.rorcr[i],
		node.rorcr[i]*final.rcr,node.w[i],mach);
      }
  }
