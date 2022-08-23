// SPDX-License-Identifier: BSD-2-Clause
/*

Copyright 2022 Spencer LaFoley

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
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "share.h"
#include <unistd.h>

int reaction_class(char Line[], int reaction);
// This function checks the reaction class for the incoming reactions such as REACTIONS KELVIN MOLECULES
/* Input parameters:
 *   1. Line[] - the current line being checked for a reaction class
 *   2. reaction - this integer represents the previous reaction class given
 * Output parameters:
 *   1. (returned int) reaction - this integer represents the new reaction class if it changed from the previous
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified
 *      = reaction when there is no change from the previous reaction class      */

int check_process(char Line[], int tot, int loc, FILE* output);
// This function checks a reaction to see if there are any unrecognizable processes it cannot compute such as "LOW" and "TROE"
/* Input parameters:
 *   1. Line[] - the current line being checked for unknown processes
 *   2. tot - the total number of reactions thus far
 *   3. loc - the current line location of the reaction under examination
 * Output parameters:
 *   1. output file - the function outputs an error to the output file in the form of a comment if there are unrecognizable processes */

void build_current_function(int run, FILE* input, FILE* output);
// This function builds the appropriate function based on the current cycle run and the given information from the input file then prints it to the output file
/* Input parameters:
 *   1. run - the indicator for which function is currently being built
 *    = 1 when building add_W_Ionization_CHEMKIN()
 *    = 2 when building add_dW_dx_Ionization_CHEMKIN()
 *    = 3 when building find_Qei()
 *    = 4 when building find_dQei_dx
 *   2. input file - the file to retrieve new lines of data from */
 /* Output parameters:
  *   1. output file - the function built is printed to the output file */


void print_function_header(int run, FILE* output);
// This function prints the appropriate function header based on the current run, prior to data retrieval
/* Input parameters:
 *   1. run - the indicator for the current function being built
 *    = 1 when building add_W_Ionization_CHEMKIN()
 *    = 2 when building add_dW_dx_Ionization_CHEMKIN()
 *    = 3 when building find_Qei()
 *    = 4 when building find_dQei_dx */
 /* Output parameters:
  *   1. output file - the function headers are printed in the output file */


void check_all_indicators_formfit_electrontemperature_Q(char oldLine[], char prevLine[], char lastLine[], char newLine[], int* ind, int* e, int* Q);
// This function checks multiple indicators such as form fit, electron temp, and Qei for a given reaction
/* Input parameters:
 *   1. oldLine[], prevLine[], lastLine[], newLine[] - the current reaction line being searched and the following 3 lines of information from the input file */
 /* Output parameters:
  *   1. ind - the indicator for form fit or standard fit
  *      = 0 when reaction is standard fit
  *      = 1 when reaction is form fit
  *   2. e - the indicator for an electron/gas temp reaction
  *      = 0 when reaction uses gas temperature
  *      = 1 when reaction uses electron temperature
  *   3. Q - the indicator for a Qei reaction
  *      = 0 when reaction does not have "EXCI" value
  *      = 1 when reaction does have "EXCI" value   */

void print_output_line(char oldLine[], char prevLine[], char lastLine[], char newLine[], int ind, int e, int Q, int way, int tot, int run, int numR, int numP, FILE* output, int M, int reaction);
// This function prints the necessary output lines according to the reaction information provided and the current function being built
/* Input parameters:
 *   1. oldLine[], prevLine[], lastLine[], newLine[] - the current reaction line being searched and the following 3 lines of information from the input file
 *   2. ind - the indicator for form fit or standard fit
 *      = 0 when reaction is standard fit
 *      = 1 when reaction is form fit
 *   3. e - the indicator for an electron temp reaction
 *      = 0 when reaction uses gas temperature
 *      = 1 when reaction uses electron temperature
 *   4. Q - the indicator for a Qei reaction
 *      = 0 when reaction does not have "EXCI" value
 *      = 1 when reaction does have "EXCI" value
 *   5. way - the indicator for a reaction that runs forwards or backwards
 *      = 1 when reaction is forward
 *      = 2 when reaction is forward and backward
 *   6. run - the indicator for what function is currently being built
 *    = 1 when building add_W_Ionization_CHEMKIN()
 *    = 2 when building add_dW_dx_Ionization_CHEMKIN()
 *    = 3 when building find_Qei()
 *    = 4 when building find_dQei_dx
 *   7. tot - the indicator for total number of reactions detected thus far, +1 for every new reaction
 *   8. numR - the number of reactants in the reaction
 *   9. numP - the number of products  in the reaction
 *   10. reaction - the reaction class of the output being built
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified  */
 /* Output parameters:
  *   1. output file - the appropriate output lines for each function are printed to the output file */


void print_elements_add_ionization(char oldLine[], char prevLine[], FILE* output, int numR, int numP, int step, int M);
// This function prints the elements on the add Ionization functions
/* Input parameters:
 *   1. oldLine[] - the current reaction line being searched
 *   2. numR - the number of reactants in the reaction
 *   3. numP - the number of products in the reaction */
 /* Output parameters:
  *   1. output file - the elements for a reaction added in the add_Ionization functions is printed to the output file  */


void print_numbers_add_ionization(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, int step, int M, int reaction);
// This function prints the data numbers for the add_Ionization functions
/* Input parameters:
 *   1. oldLine[], prevLine[] - the current reaction line being searched and the following line of information from the input file
 *   2. numR - the number of reactants in the reaction
 *   3. numP - the number of products in the reaction
 *   4. ind - the indicator for form fit or standard fit
 *      = 0 when reaction is standard fit
 *      = 1 when reaction is form fit 
 *   5. reaction - the reaction class of the output being built
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified  */
 /* Output parameters:
  *   1. output file - the data numbers for a reaction in the add_Ionization functions are printed to the output file  */

void print_info_find_Qei(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, char newLine[], char lastLine[], int Q, int e, int run, int reaction);
// This function prints information on the find Qei functions such as primary element, numbers, fit, etc
/* Input parameters:
 *   1. oldLine[], prevLine[], lastLine[], newLine[] - the current reaction line being searched and the following three lines of information from
 *      the input file
 *   2. numR - the number of reactants in the reaction
 *   3. numP - the number of products in the reaction
 *   4. ind - the indicator for form fit or standard fit
 *      = 0 when reaction is standard fit
 *      = 1 when reaction is form fit
 *   5. Q - the indicator for a Qei reaction
 *      = 0 when reaction does not have "EXCI" value
 *      = 1 when reaction does have "EXCI" value
 *   6. e - the indicator for an electron temp reaction
 *      = 0 when reaction uses gas temperature
 *      = 1 when reaction uses electron temperature
 *   7. run - the indicator for what function is currently being built
 *    = 1 when building add_W_Ionization_CHEMKIN()
 *    = 2 when building add_dW_dx_Ionization_CHEMKIN()
 *    = 3 when building find_Qei()
 *    = 4 when building find_dQei_dx 
 *   8. reaction - the reaction class of the output being built
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified  */
 /* Output parameters:
  *   1. output file - the information for a reaction in the find_Qei functions is printed to the output file  */


void print_numbers_find_Qei(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, int run, int reaction, int e);
// This function prints the data numbers excluding the EXCI variable for a reaction added to the find_Qei functions
/* Input parameters:
 *   1. oldLine[], prevLine[] - the current reaction line being searched and the following line of information from the input file
 *   2. numR - the number of reactants in the reaction
 *   3. numP - the number of products in the reaction
 *   4. ind - the indicator for form fit or standard fit
 *      = 0 when reaction is standard fit
 *      = 1 when reaction is form fit
 *   5. run - the indicator for what function is currently being built
 *    = 1 when building add_W_Ionization_CHEMKIN()
 *    = 2 when building add_dW_dx_Ionization_CHEMKIN()
 *    = 3 when building find_Qei()
 *    = 4 when building find_dQei_dx 
 *   6. reaction - the reaction class of the output being built
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified
 *   7. e - the indicator for Te or T
 *      = 0 for gas temperature
 *      = 1 for electron temperature   */
 /* Output parameters:
  *   1. output file - the data numbers for a reaction in the find_Qei functions is printed to the output file  */


void print_Qei_value(char line[], FILE* output);
// This function prints the Qei function value
/* Input parameters:
 *   1. line[] - the current line containing the EXCI value  */
 /* Output parameters:
  *   1. output file - the EXCI value is printed in the ouput file  */

int check_Q(char line[]);
// This function indicates whether or not there is the EXCI value for the Qei function
/* Input parameters:
 *   1. line[] - the current line being searched for "EXCI"   */
 /* Output parameters:
  *   1. (returned int) ind - the indicator for a reaction to be included in find_Qei function
  *      = 0 when reaction does not have "EXCI" value
  *      = 1 when reaction does have "EXCI" value  */


int check_electron_temperature(char line[]);
// This function checks if a reaction should use Te or T
/* Input parameters:
 *   1. line[] - the current line being searched for "TDEP/E/"     */
 /* Output parameters:
  *   1. (returned int) result - the indicator for electron temp or gas temp
  *      = 0 when reaction uses gas temperature
  *      = 1 when reaction uses electron temperature  */


int check_fitted_form(char Line_str[]);
// This function classifies a reaction between fitted and standard form
/* Input parameters:
 *   1. Line_str[] - the current line being searched for "FIT"   */
 /* Output parameters:
  *   1. (returned int) result - the indicator for fitted form or standard form
  *      = 0 when reaction is standard fit
  *      = 1 when reaction is form fit      */

int check_reaction_current_line(char Line_str[], int* way, int* tot, FILE* output, int loc);
// This functions checks a line for a reaction and also the direction of the reaction
/* Input parameters:
 *   1. Line_str[] - the current line of input being searched for a reaction
 *   2. loc - the indicator for the location of the line within the input file being searched */
 /* Output parameters:
  *   1. way - the indicator for a forward or backward reaction
  *      = 1 when reaction is forward
  *      = 2 when reaction is forward and backward
  *   2. tot - the total number of reactions detected thus far, +1 for each new reaction
  *   3. output file - if a reaction is ignored or unreadable an error is printed to the output file
  *   4. (returned int) num - the indicator if there is a reaction on the current line
  *      = 0 when there is no reaction on the current line
  *      = 1 when there is a reaction on the current line */


int number_of_reactants(char Line_str[]);
// This function determines the reactant count from the provided reaction line
/* Input parameters:
 *   1. Line_str[] - the current line containing the reaction   */
 /* Output parameters:
  *   1. (returned int) num - the number of reactants present in the reaction */


int number_of_products(char Line_str[]);
// This function determines the product count from the provided reaction line
/* Input parameters:
 *   1. Line_str[] - the current line containing the reaction   */
 /* Output parameters:
  *   1. (returned int) num - the number of products present in the reaction */

int check_M_reaction(char oldLine[], char prevLine[]);
// This function determines whether or not a reaction includes a third body and then determines how many steps are necessary to build each reaction substitution
/* Input parameters:
 *   1. oldLine[], prevLine[] - the current reaction line being searched and the following line of information from the input file   */
 /* Output parameters:
  *   1. (returned int) M - the indicator for a reaction with a third body
  *       = -2 when there is no third body or the third body does not have specified element substitutions
  *       = int > 0 when there is a third body with specified elements/factors, the int represents how many different elements are present  */


int main(int argc, char **argv) {
  char *inpfile, *outfile;
  bool VALIDOPTIONS = TRUE;
  char *options;
  int RET, run;
  options = NULL;
  
  
  inpfile =(char *)malloc(400*sizeof(char));
  outfile =(char *)malloc(400*sizeof(char));
  
  if (process_flag_string(argc, argv, "-i", &inpfile)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_string(argc, argv, "-o", &outfile)!=2) VALIDOPTIONS=FALSE;           // Confirm correct flag usage from the terminal
                                                                                        // Display error if incorrect flags
  if (!VALIDOPTIONS) {
    fprintf (stderr, "\nPROBLEM WITH FLAGS\n\n\nRequired and Optional Flags:\n\n"
             "Flag      Arg                     Arg Type     Required?\n"
             "----------------------------------------------------------\n"
             "-i        Input file              string       Y\n"
             "-o        Output file             string       Y\n\n"
             "Eg: \n"
             "./chemkin2cfdwarp -i input.txt -o output.c \nWill convert input.txt in Chemkin format to output.c in CFDWARP format\n\n");
    exit (EXIT_FAILURE);
  }
  RET = find_remaining_options ( argc, argv, &options );
  if ( RET >= 1 ) {
    fprintf ( stderr, "\n\nThe following command line options could not be processed:\n%s\n\n", options );
    exit (EXIT_FAILURE);
  }
  else {
    
	FILE* input = fopen(inpfile, "r");
  if (input == NULL) {                                      // Confirm input file exists in directory
    fprintf ( stderr, "\n\nThe following input file could not be found:\n%s\n\n", inpfile );
    exit (EXIT_FAILURE);
  }
	FILE* output = fopen(outfile, "w");
  
  
	for (run = 0; run != 4; run++) {              // Begin cycle for creating functions in output file

		print_function_header(run, output);
		build_current_function(run, input, output);
		fprintf(output, "}\n\n");
    }
    
  fclose(input);
	fclose(output);
  }
                                                    // Close the files and free the pointers after the output file has been built
	free ( options );
  free ( inpfile );
  free ( outfile );
  return(EXIT_SUCCESS);
}

int number_of_reactants(char Line_str[]) {
    int i = 0, j = 0, num = 1, m = 0;
    for (j = 0; Line_str[i] != '='; j++) {                           // Examine the elements of the current reaction up to a "=" sign
         if (Line_str[i] == '\t' || Line_str[i] == '\r' || Line_str[i] == '<' || Line_str[i] == '>' || Line_str[i] == ' ')
            j--, i++;
        else if (Line_str[i] == '+') {
            if (Line_str[i - 1] == ' ')
                num++;
            else if (Line_str[i + 1] == ' ' || Line_str[i+1] == '\t') {
                m = i + 1;
                while (Line_str[m] == ' ' || Line_str[m] == '\t')
                  m++;
                if (Line_str[m] != '+' && Line_str[m] != '<' && Line_str[m] != '=')
                  num++;
              }
            else if (Line_str[i + 1] != ' ' && Line_str[i + 1] != '+' && Line_str[i + 1] != '<' && Line_str[i + 1] != '=')
                num++;
            
            i++;
        }
        else if (Line_str[i] == '(' || Line_str[i] == ')') {
            i++, j--;
        }
        else if (Line_str[i] == '<' || Line_str[i] == '>' || Line_str[i] == '=')
            break;
        else
            i = i + 1;
    }
    return num; // returns the number of reactants
}

int number_of_products(char Line_str[]) {
    int i = 0, num = 1, m = 0;

    while (Line_str[i] != '=')
        i++;
    i++;
    if (Line_str[i] == '>')
        i++;
    while (Line_str[i + 2] != '.' && Line_str[i + 2] != 'e') {          // Examine the elements of the current reaction after a "=" sign
        if (Line_str[i] == '\t' || Line_str[i] == '\r')
            i++;
        else if (Line_str[i] == ' ')
            i++;
        else if (Line_str[i] == '+') {
            if (Line_str[i - 1] == ' ')
                num++;
            else if (Line_str[i + 1] == ' ' || Line_str[i+1] == '\t') {
                m = i + 1;
                while (Line_str[m] == ' ' || Line_str[m] == '\t')
                  m++;
                if (Line_str[m] != '+' && Line_str[m] != '<' && Line_str[m] != '=')
                  num++;
              }
            else if (Line_str[i + 1] != ' ' && Line_str[i + 1] != '+' && Line_str[i + 1] != '<' && Line_str[i + 1] != '=')
                num++;
            
            i++;                          // The presence of a '+' sign after the equal sign adds to the product count
        }
        else if (Line_str[i] == '(' || Line_str[i] == ')') {
            i++;
        }
        else {
          i++;
        }
    }
    return num; // returns the number of products
}

int check_electron_temperature(char line[]) {
    int i = 0, result = 0;                            // searches line for "TDEP/E/"
    for (i = 0; line[i] != '\n' && line[i] != '\r'; i++) {
        if (line[i] == 'T' && line[i + 1] == 'D' && line[i + 2] == 'E' && line[i + 3] == 'P' && line[i + 4] == '/' && line[i + 5] == 'E' && line[i + 6] == '/')
            result = 1;
    }
    return result; // 1 is returned if Te, 0 is returned if T
}

int check_fitted_form(char Line_str[]) {
    int i, num = 0;                             // searches line for "FIT"
    for (i = 0; Line_str[i] != '\n' && Line_str[i] != '\r'; i++) {
        if (Line_str[i] == 'F' && Line_str[i + 1] == 'I' && Line_str[i + 2] == 'T') {
            num = 1;
        }
    }
    return num; // returns 0 if the reaction is a standard form, 1 if it is fitted form
}
void print_function_header(int run, FILE* output) {       // print function header based on current run  
    if (run == 0) {
        fprintf(output, "void add_W_Ionization_CHEMKIN ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {\n  double N[ns];\n  long k, specM;\n  spec_t X;\n\n  for ( k = 0; k < ns; k++ ) {\n  \tX[k] = rhok[k] / _calM (k) * 1.0e-6;\t\t/* mole/cm^3 */\n  }\n\n");
    }
    if (run == 1) {
        fprintf(output, "/* Verify the validity of the dW terms at node i=10, j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU\n * Make sure to verify the dW terms over a wide range of temperatures and mass fractions\n * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives\n * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()\n */\n\nvoid add_dW_dx_Ionization_CHEMKIN ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {\n  long k, s, specM;\n  double kf, dkfdTe, dkfdT, dkfdTv;\n  spec_t X;\n\n  for ( k = 0; k < ns; k++ ) {\n  \tX[k] = rhok[k] / _calM (k) * 1.0e-6;\t\t/* mole/cm^3 */\n  }\n\n\n  /* find properties needed by add_to_dW* functions in proper units */\n\n");
    }
    if (run == 2) {
        fprintf(output, "void find_Qei(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei) {\n  double theta;\n  long specM;\n\n  *Qei = 0.0;\n\n  switch (gl->model.chem.IONIZATIONMODEL) {\n  \tcase IONIZATIONMODEL_CHEMKIN:\n");
    }
    if (run == 3) {
        fprintf(output, "void find_dQei_dx(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe) {\n  double theta;\n  long spec, specM;\n\n  for ( spec = 0; spec < ns; spec++ )\n  \tdQeidrhok[spec] = 0.0;\n\n  *dQeidTe = 0.0;\n\n  switch (gl->model.chem.IONIZATIONMODEL) {\n  \tcase IONIZATIONMODEL_CHEMKIN:\n");
    }
}

int check_Q(char line[]) {
    int i = 0, ind = 0;                 // searches line for "EXCI"

    for (i = 0; line[i] != '\0'; i++) {
        if (line[i] == 'E' && line[i + 1] == 'X' && line[i + 2] == 'C' && line[i + 3] == 'I') {
            ind = 1;
            break;
        }
    }
    return ind;             // returns 1 if "EXCI" is present, returns 0 if not
}

void print_Qei_value(char line[], FILE* output) {
    int i = 0, k = 0, j = 0;
    char Qei[10]; // initialize Qei[] for EXCI value
    while (line[i] != '/')
        i++;
    if (line[i] == '/')
        i++;
    while (line[i] != '/') {
        if (line[i] == ' ' || line[i] == '\t')
            i++;
        else {
            Qei[k] = line[i];
            k++, i++;
        }
    }
    Qei[k] = '\0';
    for (j = 0; Qei[j] != '\0'; j++) { // print Qei[] EXCI value
        fprintf(output, "%c", Qei[j]);
    }
    fprintf(output, ", ");
}

void print_info_find_Qei(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, char newLine[], char lastLine[], int Q, int e, int run, int reaction) {
    int i = 0, j = 0, k = 0, m = 0;
    char reactants[10][12]; // create a string array to hold reactant names
    for (j = 0; oldLine[i] != '='; j++) {
        k = 0;
        if (oldLine[i] == '\t' || oldLine[i] == ' ' || oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '\n' || oldLine[i] == '\r' || oldLine[i] == '(' || oldLine[i] == ')' || oldLine[i] == '+')
            j--, i++;
        else {
            while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r' && oldLine[i] != '+') {
                if (oldLine[i] == '(' || oldLine[i] == ')')
                    i++;
                else if (oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '=')
                    break;
                else {
                    reactants[j][k] = oldLine[i];
                    i = i + 1, k = k + 1;
                }
            }
            if (oldLine[i] == '+' && (oldLine[i + 1] == ' ' || oldLine[i + 1] == '\t')) {
                m = i + 1;
                while (oldLine[m] == ' ' || oldLine[m] == '\t') {
                    m++;
                }
                if (oldLine[m] == '+' || oldLine[m] == '<' || oldLine[m] == '=') {
                        reactants[j][k] = oldLine[i];
                        k++, i++;
                }
            }
            else if (oldLine[i] == '+' && oldLine[i + 1] != ' ' && oldLine[i + 1] != '\t') {
                if (oldLine[i + 1] == '+') {
                    reactants[j][k] = oldLine[i];
                    k++, i++;
                }
            }
            reactants[j][k] = '\0';
            if (reactants[j][0] == 'E' && reactants[j][1] == '\0')
                reactants[j][0] = 'e', reactants[j][1] = 'm', reactants[j][2] = 'i', reactants[j][3] = 'n', reactants[j][4] = 'u', reactants[j][5] = 's', reactants[j][6] = '\0';
            if (reactants[j][k - 1] == '+' && reactants[j][k] == '\0')
                reactants[j][k - 1] = 'p', reactants[j][k] = 'l', reactants[j][k + 1] = 'u', reactants[j][k + 2] = 's', reactants[j][k + 3] = '\0';
        }
    }
    for (j = 0; j < numR; j++) {
        if (reactants[j][0] == 'e' && reactants[j][1] == 'm') // print primary element from reactant
            j++;
        if (j < numR)
            for (i = 0; reactants[j][i] != '\0'; i++)
                fprintf(output, "%c", reactants[j][i]);
    }
    fprintf(output, ", ");
    if (Q == 1)                           // print EXCI value for reactions included in find_Qei functions
        print_Qei_value(newLine, output);
    if (Q == 2)
        print_Qei_value(prevLine, output);
    if (Q == 3)
        print_Qei_value(lastLine, output);
    if (ind == 1) { // printing fitted form numbers
        fprintf(output, "_kf_fit4(");
        print_numbers_find_Qei(oldLine, prevLine, numR, numP, ind, output, run, reaction, e);
    }
    else if (ind == 0) { // printing standard form numbers
        fprintf(output, "_kf(");
        print_numbers_find_Qei(oldLine, prevLine, numR, numP, ind, output, run, reaction, e);
    }
}

void print_numbers_add_ionization(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, int step, int M, int reaction) {
    char numsF[6][12], numsS[3][12], factor[15];
    int count = 0, n = 0, pow = 0;
    if (ind == 1) { // builds an array to store the numbers of a fitted form reaction
        int i = 0, h = 0, j = 0, k = 0;
        if (M >= 0) {
            while (count != step - 1) {
                if (prevLine[n] == '/')
                    count++;
                n++;
            }
            if (prevLine[n] == '/')
                n++;
            while (prevLine[n] != '/') {
                while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\n')
                    n++;
                if (prevLine[n] == '/')
                    break;
                factor[k] = prevLine[n];
                k++;
                n++;
            }
            factor[k] = '\0';
        }
        while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e')
            i++;
        while (prevLine[h + 2] != '.' && prevLine[h + 2] != 'e')
            h++;
        for (j = 0; j < 6; j++) {
            k = 0;
            if (j < 2) {
                while (oldLine[i] == ' ' || oldLine[i] == '\t' || oldLine[i] == '/')
                    i++;
                while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r' && oldLine[i] != '\t' && oldLine[i] != '/') {
                    numsF[j][k] = oldLine[i];
                    if (oldLine[i] == '-' && oldLine[i + 1] == ' ')
                        i++;
                    i++, k++;
                }
                if (k != 0)
                    numsF[j][k] = '\0';
                if (numsF[j][k - 1] == '\n' || numsF[j][k - 1] == '\r')
                    numsF[j][k - 1] = '\0';
                if ((numsF[j][0] == '0' && numsF[j][1] == '\0') || (numsF[j][0] == '.' && numsF[j][1] == '0' && numsF[j][2] == '0'))
                    numsF[j][0] = '0', numsF[j][1] = '.', numsF[j][2] = '0', numsF[j][3] = '\0';
                int v = strlen(numsF[j]), dec = 0, g = 0;
                for (g = 0; g != v; g++)
                    if (numsF[j][g] == '.')
                        dec = 1;
                g = 0;
                if (dec == 0) {
                    numsF[j][v] = '.', numsF[j][v + 1] = '0', numsF[j][v + 2] = '\0';
                }
            }
            else {
                while (prevLine[h] == ' ' || prevLine[h] == '\t' || prevLine[h] == '/')
                    h++;
                while (prevLine[h] != ' ' && prevLine[h] != '\n' && prevLine[h] != '\r' && prevLine[h] != '\t' && prevLine[h] != '/') {
                    numsF[j][k] = prevLine[h];
                    if (prevLine[h] == '-' && prevLine[h + 1] == ' ')
                        h++;
                    h++, k++;
                }
                if (k != 0)
                    numsF[j][k] = '\0';
                if (numsF[j][k - 1] == '\n' || numsF[j][k - 1] == '\r')
                    numsF[j][k - 1] = '\0';
                if ((numsF[j][0] == '0' && numsF[j][1] == '\0') || (numsF[j][0] == '.' && numsF[j][1] == '0' && numsF[j][2] == '0')) {
                    numsF[j][0] = '0', numsF[j][1] = '.', numsF[j][2] = '0', numsF[j][3] = '\0';
                }
                if (numsF[j][1] == 'e' || numsF[j][1] == 'E') {
                    for (int v = strlen(numsF[j]); v != 0; v--)
                        numsF[j][v + 2] = numsF[j][v];
                    numsF[j][2] = '0', numsF[j][1] = '.';
                }
                int v = strlen(numsF[j]), dec = 0;
                int g = 0;
                for (g = 0; g != v; g++)
                    if (numsF[j][g] == '.')
                        dec = 1;
                g = 0;
                if (dec == 0)
                    numsF[j][v] = '.', numsF[j][v + 1] = '0', numsF[j][v + 2] = '\0';
            }
        }
        j = 0, pow = -2;
        if (M >= 0) {
            while (factor[j] != '\0') {
                fprintf(output, "%c", factor[j]);
                j++;
            }
            fprintf(output, "*");
        }
        for (j = 0; j < 6; j++) {
            pow++;
            for (i = 0; numsF[j][i] != '\0'; i++) {   // print fitted form numbers
                fprintf(output, "%c", numsF[j][i]);
            }
            if (j == 0 && (reaction == 2 || reaction == 3))
                fprintf(output, "*calA, ");
            else if (j == 2 && (reaction == 1 || reaction == 3))
                fprintf(output, "*Rchem, ");
            else if (j > 2 && (reaction == 1 || reaction == 3)) {
                fprintf(output, "*powint(Rchem, %d), ", pow);
            }
            else
                fprintf(output, ", ");
        }
    }
    else { // builds an array to store the numbers of a standard form reaction
        int i = 0, j = 0, k = 0, n = 0, pow = 0;
        count = 0;
        if (M >= 0) {
            while (count != step - 1) {
                if (prevLine[n] == '/')
                    count++;
                n++;
            }
            if (prevLine[n] == '/')
                n++;
            while (prevLine[n] != '/') {
                while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\n')
                    n++;
                if (prevLine[n] == '/')
                    break;
                factor[k] = prevLine[n];
                k++;
                n++;
            }
            factor[k] = '\0';
        }
        while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e')
            i++;
        for (j = 0; j < 3; j++) {
            k = 0;
            while (oldLine[i] == ' ' || oldLine[i] == '\t' || oldLine[i] == '/')
                i++;
            while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r' && oldLine[i] != '\t' && oldLine[i] != '/') {
                numsS[j][k] = oldLine[i];
                if (oldLine[i] == '-' && oldLine[i + 1] == ' ')
                    i++;
                i++, k++;
            }
            if (k != 0)
                numsS[j][k] = '\0';
            if (numsS[j][k - 1] == '\n' || numsS[j][k - 1] == '\r')
                numsS[j][k - 1] = '\0';
            if ((numsS[j][0] == '0' && numsS[j][1] == '\0') || (numsS[j][0] == '.' && numsS[j][1] == '0' && numsS[j][2] == '0'))
                numsS[j][0] = '0', numsS[j][1] = '.', numsS[j][2] = '0', numsS[j][3] = '\0';
            if (numsS[j][1] == 'e' || numsS[j][1] == 'E') {
                for (int v = strlen(numsS[j]); v != 0; v--)
                    numsS[j][v + 2] = numsS[j][v];
                numsS[j][2] = '0', numsS[j][1] = '.';
            }
            int v = strlen(numsS[j]), dec = 0;
            for (int g = 0; g != v; g++)
                if (numsS[j][g] == '.')
                    dec = 1;
            if (dec == 0)
                numsS[j][v] = '.', numsS[j][v + 1] = '0', numsS[j][v + 2] = '\0';
        }
        j = 0, pow = 0;
        if (M >= 0) {
            while (factor[j] != '\0') {
                fprintf(output, "%c", factor[j]);
                j++;
            }
            fprintf(output, "*");
        }
        for (j = 0; j < 3; j++) {
            pow++;
            for (i = 0; numsS[j][i] != '\0'; i++)
                fprintf(output, "%c", numsS[j][i]);           // print standard form numbers
            if (j == 0 && (reaction == 2 || reaction == 3))
                fprintf(output, "*calA, ");
            else if (j == 2 && (reaction == 1 || reaction == 3))
                fprintf(output, "*Rchem, ");
            else if (j > 2 && (reaction == 1 || reaction == 3)) {
                fprintf(output, "*powint(Rchem, %d), ", pow);
            }
            else
                fprintf(output, ", ");
        }
    }
}

void print_numbers_find_Qei(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, int run, int reaction, int e) {
    char numsF[6][12], numsS[3][12];
    int pow = 0;
    if (ind == 1) { // builds an array to store the numbers of a fitted form reaction
        int i = 0, h = 0, j = 0, k = 0;
        while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e')
            i++;
        while (prevLine[h + 2] != '.' && prevLine[h + 2] != 'e')
            h++;
        for (j = 0; j < 6; j++) {
            k = 0;
            if (j < 2) {
                while (oldLine[i] == ' ' || oldLine[i] == '\t' || oldLine[i] == '/')
                    i++;
                while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r' && oldLine[i] != '\t' && oldLine[i] != '/') {
                    numsF[j][k] = oldLine[i];
                    if (oldLine[i] == '-' && oldLine[i + 1] == ' ')
                        i++;
                    i++, k++;
                }
                if (k != 0)
                    numsF[j][k] = '\0';
                if (numsF[j][k - 1] == '\n' || numsF[j][k - 1] == '\r')
                    numsF[j][k - 1] = '\0';
                if ((numsF[j][0] == '0' && numsF[j][1] == '\0') || (numsF[j][0] == '.' && numsF[j][1] == '0' && numsF[j][2] == '0'))
                    numsF[j][0] = '0', numsF[j][1] = '.', numsF[j][2] = '0', numsF[j][3] = '\0';
                int v = strlen(numsF[j]), dec = 0, g = 0;
                for (g = 0; g != v; g++)
                    if (numsF[j][g] == '.')
                        dec = 1;
                g = 0;
                if (dec == 0) {
                    numsF[j][v] = '.', numsF[j][v + 1] = '0', numsF[j][v + 2] = '\0';
                }
            }
            else {
                while (prevLine[h] == ' ' || prevLine[h] == '\t' || prevLine[h] == '/')
                    h++;
                while (prevLine[h] != ' ' && prevLine[h] != '\n' && prevLine[h] != '\r' && prevLine[h] != '\t' && prevLine[h] != '/') {
                    numsF[j][k] = prevLine[h];
                    if (prevLine[h] == '-' && prevLine[h + 1] == ' ')
                        h++;
                    h++, k++;
                }
                if (k != 0)
                    numsF[j][k] = '\0';
                if ((numsF[j][0] == '0' && numsF[j][1] == '\0') || (numsF[j][0] == '.' && numsF[j][1] == '0' && numsF[j][2] == '0')) {
                    numsF[j][0] = '0', numsF[j][1] = '.', numsF[j][2] = '0', numsF[j][3] = '\0';
                }
                if (numsF[j][1] == 'e' || numsF[j][1] == 'E') {
                    for (int v = strlen(numsF[j]); v != 0; v--)
                        numsF[j][v + 2] = numsF[j][v];
                    numsF[j][2] = '0', numsF[j][1] = '.';
                }
                int v = strlen(numsF[j]), dec = 0;
                int g = 0;
                for (g = 0; g != v; g++)
                    if (numsF[j][g] == '.')
                        dec = 1;
                g = 0;
                if (dec == 0)
                    numsF[j][v] = '.', numsF[j][v + 1] = '0', numsF[j][v + 2] = '\0';
            }
        }
        pow = -2;
        for (j = 0; j < 6; j++) {
            pow++;
            for (i = 0; numsF[j][i] != '\0'; i++) {     // print the numbers of a fit form reaction
                fprintf(output, "%c", numsF[j][i]);
            }
            if (j == 0 && (reaction == 2 || reaction == 3))
                fprintf(output, "*calA, ");
            else if (j == 2 && (reaction == 1 || reaction == 3))
                fprintf(output, "*Rchem, ");
            else if (j > 2 && (reaction == 1 || reaction == 3)) {
              if (j != 5)
                fprintf(output, "*powint(Rchem, %d), ", pow);
              else
                fprintf(output, "*powint(Rchem, %d)", pow);
            }
            else if (j != 5)
                fprintf(output, ", ");
            
        }
        if (e == 1)
              fprintf(output, ", Te)/calA");
            else
              fprintf(output, ", T)/calA");
            if (run == 2 && j == 6)
                fprintf(output, ", rhok, Qei);\n");
            else if (run == 3 && j == 6)
                fprintf(output, ", 0.0, rhok, dQeidrhok, dQeidTe);\n");
    }
    else { // builds an array to store the numbers of a standard form reaction
        int i = 0, j = 0, k = 0;
        while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e')
            i++;
        for (j = 0; j < 3; j++) {
            k = 0;
            while (oldLine[i] == ' ' || oldLine[i] == '\t' || oldLine[i] == '/')
                i++;
            while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r' && oldLine[i] != '\t' && oldLine[i] != '/') {
                numsS[j][k] = oldLine[i];
                if (oldLine[i] == '-' && oldLine[i + 1] == ' ')
                    i++;
                i++, k++;
            }
            if (k != 0)
                numsS[j][k] = '\0';
            if ((numsS[j][0] == '0' && numsS[j][1] == '\0') || (numsS[j][0] == '.' && numsS[j][1] == '0' && numsS[j][2] == '0'))
                numsS[j][0] = '0', numsS[j][1] = '.', numsS[j][2] = '0', numsS[j][3] = '\0';
            if (numsS[j][1] == 'e' || numsS[j][1] == 'E') {
                for (int v = strlen(numsS[j]); v != 0; v--)
                    numsS[j][v + 2] = numsS[j][v];
                numsS[j][2] = '0', numsS[j][1] = '.';
            }
            int v = strlen(numsS[j]), dec = 0;
            for (int g = 0; g != v; g++)
                if (numsS[j][g] == '.')
                    dec = 1;
            if (dec == 0)
                numsS[j][v] = '.', numsS[j][v + 1] = '0', numsS[j][v + 2] = '\0';
        }
        pow = 0;
        for (j = 0; j < 3; j++) {
            for (i = 0; numsS[j][i] != '\0'; i++)
                fprintf(output, "%c", numsS[j][i]);           // print numbers of standard reaction
            if (j == 0 && (reaction == 2 || reaction == 3))
                fprintf(output, "*calA, ");
            else if (j == 2 && (reaction == 1 || reaction == 3))
                fprintf(output, "*Rchem, ");
            else if (j > 2 && (reaction == 1 || reaction == 3)) {
              if (j != 2)
                fprintf(output, "*powint(Rchem, %d), ", pow);
              else
                fprintf(output, "*powint(Rchem, %d)", pow);
            }
            else if (j != 5)
                fprintf(output, ", ");
            
        }
        if (e == 1)
              fprintf(output, ", Te)/calA");
        else
              fprintf(output, ", T)/calA");
        if (run == 2 && j == 3)
                fprintf(output, ", rhok, Qei);\n");
        else if (run == 3 && j == 3)
                fprintf(output, ", 0.0, rhok, dQeidrhok, dQeidTe);\n");
    }
    fprintf(output, "\n");
}

void print_elements_add_ionization(char oldLine[], char prevLine[], FILE* output, int numR, int numP, int step, int M) {
    int i = 0, j = 0, k = 0, n = 0, count = 0, m = 0;
    char reactants[10][12]; // create a string array to hold reactant names
    for (j = 0; oldLine[i] != '='; j++) {
        k = 0;
        if (oldLine[i] == '\t' || oldLine[i] == ' ' || oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '\n' || oldLine[i] == '\r' || oldLine[i] == '+' || oldLine[i] == '(' || oldLine[i] == ')') {
            j--, i++;
        }
        else {
            while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r' && oldLine[i] != '+') {
                if (oldLine[i] == '(' || oldLine[i] == ')') {
                    i++;
                }
                else if (oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '=')
                    break;
                else {
                    reactants[j][k] = oldLine[i];
                    i = i + 1, k = k + 1;
                }
            }
            if (oldLine[i] == '+' && (oldLine[i + 1] == ' ' || oldLine[i + 1] == '\t')) {
                m = i + 1;
                while (oldLine[m] == ' ' || oldLine[m] == '\t') {
                    m++;
                }
                if (oldLine[m] == '+' || oldLine[m] == '<' || oldLine[m] == '=') {
                        reactants[j][k] = oldLine[i];
                        k++, i++;
                    }
            }
            else if (oldLine[i] == '+' && oldLine[i + 1] != ' ' && oldLine[i + 1] != '\t') {
                if (oldLine[i + 1] == '+' || oldLine[i + 1] == '<' || oldLine[i + 1] == '=') {
                    reactants[j][k] = oldLine[i];
                    k++, i++;
                }
            }
            reactants[j][k] = '\0';
            if (reactants[j][0] == 'E' && reactants[j][1] == '\0')
                reactants[j][0] = 'e', reactants[j][1] = 'm', reactants[j][2] = 'i', reactants[j][3] = 'n', reactants[j][4] = 'u', reactants[j][5] = 's', reactants[j][6] = '\0';
            if (reactants[j][k - 1] == '+' && reactants[j][k] == '\0')
                reactants[j][k - 1] = 'p', reactants[j][k] = 'l', reactants[j][k + 1] = 'u', reactants[j][k + 2] = 's', reactants[j][k + 3] = '\0';
            k = 0, n = 0;
            if (reactants[j][0] == 'M' && reactants[j][1] == '\0' && M >= 0) {
                while (count != step - 1) {
                    n++;
                    if (prevLine[n] == '/')
                        count++;
                }
                if (prevLine[n] == '/')
                    n++;
                while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\r' || prevLine[n] == '\n')
                    n++;
                while (prevLine[n] != '/') {
                    while (prevLine[n] == ' ' || prevLine[n] == '\t')
                        n++;
                    if (prevLine[n] == '/')
                        break;
                    reactants[j][k] = prevLine[n];
                    k++;
                    n++;
                }
                count++;

                reactants[j][k] = '\0';
            }
        }
    }
    i = 0, j = 0, k = 0;
    char products[10][12]; // create a string array to hold product names
    while (oldLine[i] != '=')
        i++;
    i++;
    if (oldLine[i] == '>')
        i++;
    while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e') {
        k = 0;
        if (oldLine[i] == '\t')
            i++;
        else if (oldLine[i] == ' ')
            i++;
        else if (oldLine[i] == '+')
            i++;
        else if (oldLine[i] == '(' || oldLine[i] == ')') {
            oldLine[i] = ' ';
            i++;
        }
        else {
            while (oldLine[i] != ' ' && oldLine[i] != '\t' && oldLine[i] != '+') {
                if (oldLine[i] == '(' || oldLine[i] == ')')
                    i++;
                else {
                    products[j][k] = oldLine[i];
                    i++, k++;
                }
            }
            if (oldLine[i] == '+' && (oldLine[i + 1] == ' ' || oldLine[i + 1] == '\t')) {
                m = i + 1;
                    while (oldLine[m] == ' ' || oldLine[m] == '\t') {
                        m++;
                      }
                    if (oldLine[m] == '+' || oldLine[m] == '<' || oldLine[m] == '=') {
                            products[j][k] = oldLine[i];
                            k++, i++;
                      }
            }
            else if (oldLine[i] == '+' && oldLine[i + 1] != ' ' && oldLine[i + 1] != '\t') {
                if (oldLine[i + 1] == '+' || oldLine[i + 1] == '<' || oldLine[i + 1] == '=') {
                    products[j][k] = oldLine[i];
                    k++, i++;
                }
            }
            products[j][k] = '\0';
            if (products[j][0] == 'E' && products[j][1] == '\0')
                products[j][0] = 'e', products[j][1] = 'm', products[j][2] = 'i', products[j][3] = 'n', products[j][4] = 'u', products[j][5] = 's', products[j][6] = '\0';
            if (products[j][k - 1] == '+' && products[j][k] == '\0')
                products[j][k - 1] = 'p', products[j][k] = 'l', products[j][k + 1] = 'u', products[j][k + 2] = 's', products[j][k + 3] = '\0';
            k = 0, n = 0, count = 0;
            if (products[j][0] == 'M' && products[j][1] == '\0' && M >= 0) {
                while (count != step - 1) {
                    n++;
                    if (prevLine[n] == '/')
                        count++;
                }
                if (prevLine[n] == '/')
                    n++;
                while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\r' || prevLine[n] == '\n')
                    n++;
                while (prevLine[n] != '/') {
                    while (prevLine[n] == ' ' || prevLine[n] == '\t')
                        n++;
                    if (prevLine[n] == '/')
                        break;
                    products[j][k] = prevLine[n];
                    k++;
                    n++;
                    count++;
                    products[j][k] = '\0';
                }
            }
            j++;
        }
    }
    for (j = 0; j < numR; j++) { // printing reaction names
        if (j == 0) {
            fprintf(output, "spec");
            for (i = 0; reactants[j][i] != '\0'; i++)
                fprintf(output, "%c", reactants[j][i]);
        }
        else
            for (i = 0; reactants[j][i] != '\0'; i++)
                fprintf(output, "%c", reactants[j][i]);
        fprintf(output, ", spec");
    }
    for (j = 0; j < numP; j++) { // printing product names
        if (j != (numP - 1)) {
            for (i = 0; products[j][i] != '\0'; i++)
                fprintf(output, "%c", products[j][i]);
            fprintf(output, ", spec");
        }
        else {
            for (i = 0; products[j][i] != '\0'; i++)
                fprintf(output, "%c", products[j][i]);
            fprintf(output, ", ");
        }
    }
}


void print_output_line(char oldLine[], char prevLine[], char lastLine[], char newLine[], int ind, int e, int Q, int way, int tot, int run, int numR, int numP, FILE* output, int M, int reaction) {

    int step = 0;
    while (M > 0 || M == -2 || M == -4) {
        step++;
        M = M - 1;
        if (M == -3 || M == -5)
            step = 0;
        if (run == 0) { // Begin printing first function
            if (step == 0 || step == 1) {
                fprintf(output, "  if (IONIZATIONREACTION[%d])", tot); // begin printing information in output file
            }
            if (step == 1 || M == -5)
                fprintf(output, "\t{\n");
            else
                fprintf(output, "\n");
            if (M == -5)
                fprintf(output, "\t  for (specM = 0; specM < ns; specM++)\n\t");
            if (way == 2)
                fprintf(output, "\t  add_to_W_fwbw_%dr%dp", numR, numP);
            else
                fprintf(output, "\t  add_to_W_fw_%dr%dp", numR, numP);
            if (ind == 1)
                fprintf(output, "_fit4("); // printing fit information
            else
                fprintf(output, "(");
            print_elements_add_ionization(oldLine, prevLine, output, numR, numP, step, M);
            step++;
            if (M == -3 || M == -5)
                step = 0;
            print_numbers_add_ionization(oldLine, prevLine, numR, numP, ind, output, step, M, reaction);
            if (e == 1) // printing temperature Te vs T information
                fprintf(output, "Te, X, W);");
            else
                fprintf(output, "T, X, W);");
            if (M == -3)
                fprintf(output, "\n\n");
        }
        else if (run == 1) {
            if (step == 0 || step == 1) {
                fprintf(output, "  if (IONIZATIONREACTION[%d])", tot); // begin printing information in output file
            }
            if (step == 1 || M == -5)
                fprintf(output, "\t{\n");
            else
                fprintf(output, "\n");
            if (M == -5)
                fprintf(output, "\t  for (specM = 0; specM < ns; specM++)\n\t");
            if (way == 2)
                fprintf(output, "\t  add_to_dW_fwbw_%dr%dp", numR, numP);
            else
                fprintf(output, "\t  add_to_dW_fw_%dr%dp", numR, numP);
            if (ind == 1)
                fprintf(output, "_fit4("); // printing fit information
            else
                fprintf(output, "(");
            print_elements_add_ionization(oldLine, prevLine, output, numR, numP, step, M);
            step++;
            if (M == -3 || M == -5)
                step = 0;
            print_numbers_add_ionization(oldLine, prevLine, numR, numP, ind, output, step, M, reaction);
            if (e == 1) // printing temperature Te vs T information
                fprintf(output, "Te, X, dWdTe, dWdrhok);");
            else
                fprintf(output, "T, X, dWdT, dWdrhok);");
            if (M == -3)
                fprintf(output, "\n\n");
        }
        else if (run == 2 && Q != 0) {
            fprintf(output, "  \t\tif (IONIZATIONREACTION[%d])\n", tot);
            fprintf(output, "  \t\t\tadd_to_Qei(spec");                  // printing find_Qei functions
            print_info_find_Qei(oldLine, prevLine, numR, numP, ind, output, newLine, lastLine, Q, e, run, reaction);

        }
        else if (run == 3 && Q != 0) {
            fprintf(output, "  \t\tif (IONIZATIONREACTION[%d])\n", tot);
            fprintf(output, "  \t\t\tadd_to_dQei(spec");                // printing find_Qei functions
            print_info_find_Qei(oldLine, prevLine, numR, numP, ind, output, newLine, lastLine, Q, e, run, reaction);

        }

    }
    if (run == 0 || run == 1) {
        if (M != -3)
            fprintf(output, "\n\t}\n\n");
    }
}


void check_all_indicators_formfit_electrontemperature_Q(char oldLine[], char prevLine[], char lastLine[], char newLine[], int* ind, int* e, int* Q) {
    *ind = check_fitted_form(prevLine);                           // classify a reaction using the indicators: ind, e, Q
    if (*ind == 1) {
        *e = check_electron_temperature(lastLine);
        if (check_Q(newLine) == 1)
            *Q = 1;
    }
    else {
        *e = check_electron_temperature(prevLine);
        if (check_Q(prevLine) == 1)
            *Q = 2;
        if (*Q == 0)
            if (check_Q(lastLine) == 1)
                *Q = 3;
    }
}

void build_current_function(int run, FILE* input, FILE* output) {
    char newLine[1000], oldLine[1000], prevLine[1000], lastLine[1000];
    int loc = 0, numR = 0, numP = 0, tot = 0, ind = 0, check = 0, e = 0, way = 10, Q = 0, M = 0, flag1 = 0, flag2 = 0, flag3 = 0, flag4 = 0, reaction = 0;

    rewind(input);
    fgets(oldLine, 1000, input);
    fgets(prevLine, 1000, input);
    fgets(lastLine, 1000, input);

    while (fgets(newLine, 1000, input) != NULL) {  // searches each line until the end of the document
        loc++;
        if (loc == 41)
            loc = loc;
        reaction = reaction_class(oldLine, reaction);
        Q = 0, way = 10, e = 0, check = 0, ind = 0;
        check = check_reaction_current_line(oldLine, &way, &tot, output, loc);
        if (check == 1 || check == 3) {
            flag1 = check_process(oldLine, tot, loc, output);                     // flag a reaction as unreadable if any of the information lines
            flag2 = check_process(prevLine, tot, loc, output);                    // contain unrecognizable processes
            if (flag2 != 3 && flag2 != 1) {
                flag3 = check_process(lastLine, tot, loc, output);
                if (flag3 != 3 && flag3 != 1)
                    flag4 = check_process(newLine, tot, loc, output);
            }
            if (flag1 == 1 || flag2 == 1 || flag3 == 1 || flag4 == 1)
                check = 0;
        }
        if (check == 1 || check == 3) {
            check_all_indicators_formfit_electrontemperature_Q(oldLine, prevLine, lastLine, newLine, &ind, &e, &Q);
            M = check_M_reaction(oldLine, prevLine);
            numR = number_of_reactants(oldLine);
            numP = number_of_products(oldLine);
            print_output_line(oldLine, prevLine, lastLine, newLine, ind, e, Q, way, tot, run, numR, numP, output, M, reaction);
        }
        strcpy(oldLine, prevLine); // shifts line storage to next line
        strcpy(prevLine, lastLine);
        strcpy(lastLine, newLine);

    }
}


int check_reaction_current_line(char Line_str[], int* way, int* tot, FILE* output, int loc) {
    int i, num = 0;
    for (i = 0; Line_str[i] != '\n' && Line_str[i] != '\r'; i++) {              // determines the presence of a reaction by "=>" or "<=>"
        if (Line_str[i] == '=' && Line_str[i + 1] == '>' && Line_str[i - 1] == '<') {
            *tot = *tot + 1;
            num = 1;
            *way = 2;                                                 // way indicates whether the reaction is forwards or backwards
            break;
        }
        if (Line_str[i] == '=' && Line_str[i+1] == '>') {
            *tot = *tot + 1;
            num = 1;
            *way = 1;
            break;
        }
    }
    if (Line_str[0] == '!' && num == 1) {
        num = 0;
        fprintf(output, "\n  \t\t\t// ************ERROR: Reaction #%d ignored on line %d due to commenting out.\n\n", *tot, loc);
        return num;
        }
    if (num == 1)
      return num;
    else {
      
      for (i = 0; Line_str[i] != '\n' && Line_str[i] != '\r'; i++) {    
        if (Line_str[i] == '=') {
          num = 1;
          *tot = *tot + 1;
          *way = 2;
        }
      }
      if (Line_str[0] == '!' && num == 1) {
        num = 0;
        return num;
      }
    }
    if (num == 1)
        fprintf(output, "\n  \t// **********WARNING: The following reaction, Reaction #%d on line %d, interpreted '=' as '<=>'.\n", *tot, loc);
    


    return num; // returns 0 if there is no reaction and 1 if there is a reaction
}

int check_M_reaction(char oldLine[], char prevLine[]) {
    int i = 0, j = 0, k = 0, M = 0, m = 0;
    char reactants[10][12]; // create a string array to hold reactant names
    for (j = 0; oldLine[i] != '='; j++) {
        k = 0;
        if (oldLine[i] == '\t' || oldLine[i] == ' ' || oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '\n' || oldLine[i] == '\r' || oldLine[i] == '+' || oldLine[i] == '(' || oldLine[i] == ')') {
            j--, i++;
        }
        else {
            while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r' && oldLine[i] != '+') {
                if (oldLine[i] == '(' || oldLine[i] == ')') {
                    i++;
                }
                else if (oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '=') {
                    break;
                }
                else {
                    reactants[j][k] = oldLine[i];
                    i = i + 1, k = k + 1;
                }
            }
            if (oldLine[i] == '+' && (oldLine[i + 1] == ' ' || oldLine[i + 1] == '\t')) {
                m = i + 1;
                while (oldLine[m] == ' ' || oldLine[m] == '\t') {
                    m++;
                }
                if (oldLine[m] == '+' || oldLine[m] == '<' || oldLine[m] == '=') {
                        reactants[j][k] = oldLine[i];
                        k++, i++;
                    }
            }
            else if (oldLine[i] == '+' && oldLine[i + 1] != ' ' && oldLine[i + 1] != '\t') {
                if (oldLine[i + 1] == '+') {
                    reactants[j][k] = oldLine[i];
                    k++, i++;
                }
            }
            reactants[j][k] = '\0';
            if (reactants[j][0] == 'E' && reactants[j][1] == '\0')
                reactants[j][0] = 'e', reactants[j][1] = 'm', reactants[j][2] = 'i', reactants[j][3] = 'n', reactants[j][4] = 'u', reactants[j][5] = 's', reactants[j][6] = '\0';
            if (reactants[j][0] == 'M' && reactants[j][1] == '\0')              // determine if any include an M third body
                M = 1;
        }
    }
    if (M == 1) {
        int i = 0, count = 0;
        for (i = 0; prevLine[i] != '\n' && prevLine[i] != '\r'; i++) {
            if (prevLine[i] == '=' && prevLine[i + 1] == '>') {
                M = -4;                                                 // if there is a third body but not enough info provided for substitutions
                break;
            }
        }
        if (M == 1) {
            for (i = 0; prevLine[i] != '\n' && prevLine[i] != '\r'; i++) {
                if (prevLine[i] == '/')
                    count = count + 1;
            }
            M = count / 2;                                    // count for how many substitutions are necessary if third body info is present
            if (count == 0)
                M = -4;
        }
    }
    else
        M = -2;
    return M;                   // returns number of substitutions necessary if third body is present, returns -2 if no third body, returns -4 if third body
                                // but not enough information

}

int check_process(char Line[], int tot, int loc, FILE* output) {
    char word[12];
    char wordbank[44][15] = { "E", "N2", "N2(A3Sigma)", "N2(B3Pi)", "N2(ap1Sigma)", "N2(C3Pi)", "N2+", "O2", "O2+", "N", "O", "O(1D)", "O(1S)", "C2H4+", "AR", "HE", "H", "H2", "OH", "HO2", "H2O", "H2O2", "CO", "CO2", "HCO", "CH3", "CH4", "C2H6", "CH2O", "C2H5", "CH2", "CH3O", "CH2OH", "CH", "C2H2", "C2H4", "C2H3", "CH3OH", "CH3HCO", "C2H", "CH2CO", "HCCO", "NO", "END" };
    int i, k = 0;

    for (i = 0; Line[i] != '\n' && Line[i] != '\r'; i++) {              // determines the presence of a reaction by "=>" or "<=>"
        if (Line[i] == '=')
            return 3;
    }
    i = 0;
    if (Line[i] == '!')
        return 0;
    while (Line[i] != '\n' && Line[i] != '\r' && Line[i] != '/') {
        while (Line[i] == ' ')
            i++;
        if (Line[i] == '/')
            break;
        word[k] = Line[i];
        k++, i++;
    }
    word[k] = '\0';
    i = 0;
    if (word[i] == 'F' && word[i + 1] == 'I' && word[i + 2] == 'T')
        return 0;
    else if (word[i] == 'T' && word[i + 1] == 'D' && word[i + 2] == 'E' && word[i + 3] == 'P')
        return 0;
    else if (word[i] == 'E' && word[i + 1] == 'X' && word[i + 2] == 'C' && word[i + 3] == 'I')
        return 0;
    else if (word[i] == 'D' && word[i + 1] == 'U' && word[i + 2] == 'P')
        return 0;
    else if (word[i] == '\0')
        return 0;
    for (i = 0; i < 44; i++)
        if (strcmp(word, wordbank[i]) == 0)
            return 0;
    fprintf(output, "\n  \t\t\t// ************ERROR: Reaction #%d ignored on line %d due to unknown process: \"%s\".\n\n", tot, loc, word);

    return 1;

}

int reaction_class(char Line[], int reaction) {
    int i = 0, ind = 0, j = 0, k = 0, K = 0, M = 0;
    char word[3][50];
    if (Line[i] == '!')
        return reaction;
    while (Line[i] == ' ' || Line[i] == '\t')
        i++;
    while (Line[i] != ' ' && Line[i] != '\t' && Line[i] != '\n' && Line[i] != '\r') {
        word[j][k] = Line[i];
        i++, k++;
    }
    word[j][k] = '\0';
    if (strcmp(word[0], "REACTIONS") == 0)
        ind = 1;
    else
        return reaction;
    j++;
    if (ind == 1) {
        while (Line[i] != '\n' && Line[i] != '\r') {
            k = 0;
            while (Line[i] == ' ' || Line[i] == '\t')
                i++;
            while (Line[i] != ' ' && Line[i] != '\t' && Line[i] != '\n' && Line[i] != '\r') {
                word[j][k] = Line[i];
                i++, k++;
            }
            word[j][k] = '\0';
            j++;
        }
        if (strcmp(word[1], "KELVIN") == 0)
            K = 1;
        if (strcmp(word[1], "MOLECULES") == 0)
            M = 1;
        if (strcmp(word[2], "KELVIN") == 0)
            K = 1;
        if (strcmp(word[2], "MOLECULES") == 0)
            M = 1;
        if (M == 1 && K == 1)
            return 3;
        if (M == 1 && K == 0)
            return 2;
        if (M == 0 && K == 1)
            return 1;
        if (M == 0 && K == 0)
            return 0;
    }
    return reaction;                                        // returns an integer for the reaction class:  returns 0 if there is no class specified,
                                                            // returns 1 if only Kelvin is specified, returns 2 if only molecules is specified, returns 3 if both are specified
}

