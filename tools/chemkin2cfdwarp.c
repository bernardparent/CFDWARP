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
 *      = reaction  --  when there is no change from the previous reaction class      */

int check_process(char Line[], int tot, int loc, FILE* output, int num, int run, char* species);
// This function checks a reaction to see if there are any unrecognizable processes it cannot compute such as "LOW" and "TROE". This function also classifies much of a reaction by determing how many information lines are given after a reaction is listed in the input file. 
/* Input parameters:
 *   1. Line[] - the current input line being checked for unknown processes
 *   2. tot - the total number of reactions thus far 
 *   3. loc - the current line location of the reaction under examination
 *   4. run - the indicator for which function is currently being built (here it is used to avoid printing an error during run = -1)
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx
 *   5. num - the indicator for whether Line[] is being checked for processes or within another function such as print_elements_chemkin or check_M_reaction
 *      = 0  when being used within print_elements_chemkin or check_M_reaction
 *      = 1  when processes are being checked in build_current_function
 *   6. species - the list of species to be compared with known processes
 * Output parameters:
 *   1. output file - the function outputs an error to the output file in the form of a comment if there are unrecognizable processes 
 *   2. Return integer - the returned integer specifies the type of line in the input
 *      = -2 when there is a "LOW" process specified
 *      = 0 when line is commented out or a known process is read
 *      = 1 when there is an unrecognizeable process specified
 *      = 2 when there is a "HIGH" process specified
 *      = 3 when there is a reaction is read
 *      = 10 when there is a species read     */


void write_species(FILE* input, char** species, int* loc);
// This function builds a species list when prompted by "SPECIES" in the input file. The species list is built in a memory allocated character array  char* p and the address of this character array is stored in char** species which is then able to be returned from the function and allow for species list recovery.
/* Input parameters:
 *   1. input file - the file to retrieve new lines of data from
 *   2. loc - the current line location of the input under examination (updated value returned) */
/* Output parameters:
 *   1. species - returns the pointer to the character array of the species list */

int build_current_function(int run, FILE* input, FILE* output);
// This function builds the appropriate function based on the current cycle run and the given information from the input file then prints it to the output file
/* Input parameters:
 *   1. run - the indicator for which function is currently being built
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx
 *   2. input file - the file to retrieve new lines of data from */
 /* Output parameters:
  *   1. output file - the function built is printed to the output file */


void print_function_header(int run, FILE* output);
// This function prints the appropriate function header based on the current run, prior to data retrieval
/* Input parameters:
 *   1. run - the indicator for which function is currently being built
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx */
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
  *      = -2 when reaction does have "EXCI" value that is not specified
  *      = 0 when reaction does not have "EXCI" value
  *      = 1 when reaction does have "EXCI" value on the fourth line of information
  *      = 2 when reaction does have "EXCI" value on the second line of information
  *      = 3 when reaction does have "EXCI" value on the third line of information   */

void print_output_line(char oldLine[], char prevLine[], char lastLine[], char newLine[], int ind, int e, int Q, int way, int tot, int run, int numR, int numP, FILE* output, int M, int reaction, int h, int loc, int ply, char* species);
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
 *       = -2 when reaction does have "EXCI" value that is not specified
  *      = 0 when reaction does not have "EXCI" value
  *      = 1 when reaction does have "EXCI" value on the fourth line of information
  *      = 2 when reaction does have "EXCI" value on the second line of information
  *      = 3 when reaction does have "EXCI" value on the third line of information  
 *   5. way - the indicator for a reaction that runs forwards or backwards
 *      = 1 when reaction is forward
 *      = 2 when reaction is forward and backward
 *   6. run - the indicator for which function is currently being built
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx
 *   7. tot - the indicator for total number of reactions detected thus far, +1 for every new reaction
 *   8. numR - the number of reactants in the reaction
 *   9. numP - the number of products  in the reaction
 *   10. M - the indicator for the presence of a third body reaction or the number of necessary third body substitutions if yes
 *   11. reaction - the reaction class of the output being built
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified 
 *   12. h - the indicator for HIGH/LOW process and what line it is located on
 *   13. loc - the current line location of the input under examination
 *   14. ply - the indicator for a third body presence in the reaction, regardless of substitutions
 *   15. species - the species list to be used in third body reactions and process comparisons
 *  */
 /* Output parameters:
  *   1. output file - the appropriate output lines for each function are printed to the output file */


void print_elements_chemkin(char oldLine[], char prevLine[], char lastLine[], char newLine[], FILE* output, int numR, int numP, int step, int M, int ind, int run, int way, int tot, int loc, int w, char* species);
// This function prints the elements on the find_W_CHEMKIN functions
/* Input parameters:
 *   1. oldLine[], prevLine[], lastLine[], newLine[] - the current reaction line being searched and the following 3 lines of information from the input file
 *   2. numR - the number of reactants in the reaction
 *   3. numP - the number of products in the reaction 
 *   4. step - the current iteration of substitution steps made for a third body reaction
 *   5. M - the number of subsitutions necessary or an indicator of a third body reaction
 *   6. ind - the indicator for form fit or standard fit
 *      = 0 when reaction is standard fit
 *      = 1 when reaction is form fit
 *   7. run - the indicator for which function is currently being built
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx 
 *   8. way - the indicator for a reaction that runs forwards or backwards
 *      = 1 when reaction is forward
 *      = 2 when reaction is forward and backward
 *   9. tot - the indicator for total number of reactions detected thus far, +1 for every new reaction
 *   10. loc - the current line location of the input under examination
 *   11. w - this is an indicator for high or low processes and which order to print variables 
 *   12. species - the species list to be used in third body reactions and process comparisons */
 /* Output parameters:
  *   1. output file - the elements for a reaction added in the _chemkin functions is printed to the output file  */


void print_numbers_chemkin(char oldLine[], char prevLine[], char lastLine[], char newLine[], int numR, int numP, int ind, FILE* output, int step, int M, int reaction, int w);
// This function prints the data numbers for the find_W_CHEMKIN functions
/* Input parameters:
 *   1. oldLine[], prevLine[], lastLine[], newLine[] - the current reaction line being searched and the following 3 lines of information from the input file
 *   2. numR - the number of reactants in the reaction
 *   3. numP - the number of products in the reaction
 *   4. step - the current iteration of substitution steps made for a third body reaction
 *   5. ind - the indicator for form fit or standard fit
 *      = 0 when reaction is standard fit
 *      = 1 when reaction is form fit
 *   6. reaction - the reaction class of the output being built
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified  
 *   7. M - the number of subsitutions necessary or an indicator of a third body reaction
 *   8. reaction - the reaction class of the output being built
 *   9. w - this is an indicator for high or low processes and which order to print variables   */
 /* Output parameters:
  *   1. output file - the data numbers for a reaction in the add_chemkin functions are printed to the output file  */

void print_info_find_Qei(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, char newLine[], char lastLine[], int Q, int e, int run, int reaction, int tot);
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
 *       = -2 when reaction does have "EXCI" value that is not specified
  *      = 0 when reaction does not have "EXCI" value
  *      = 1 when reaction does have "EXCI" value on the fourth line of information
  *      = 2 when reaction does have "EXCI" value on the second line of information
  *      = 3 when reaction does have "EXCI" value on the third line of information  
 *   6. e - the indicator for an electron temp reaction
 *      = 0 when reaction uses gas temperature
 *      = 1 when reaction uses electron temperature
 *   7. run - the indicator for which function is currently being built
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx
 *   8. reaction - the reaction class of the output being built
 *      = 0 when no reaction class is specified
 *      = 1 when only Kelvin is specified
 *      = 2 when only Molecules is specified
 *      = 3 when both Kelvin and Molecules are specified 
 *   9. tot - total number of reacitons */
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
 *   5. run - the indicator for which function is currently being built
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx
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
  *      = 1 when reaction does have "EXCI" value 
  *      = 2 when reaction does have EXCI value but it is not provided */


int check_electron_temperature(char line[]);
// This function checks if a reaction should use Te or T
/* Input parameters:
 *   1. line[] - the current line being searched for "TDEP/E/"     */
 /* Output parameters:
  *   1. (returned int) result - the indicator for electron temp or gas temp
  *      = 0 when reaction uses gas temperature
  *      = 1 when reaction uses electron temperature  */


int check_fitted_form(char Line_str[]);
// This function classifies a reaction between fitted and standard form and returns an integer indicating the classification
/* Input parameters:
 *   1. Line_str[] - the current line being searched for "FIT"   */
 /* Output parameters:
  *   1. (returned int) result - the indicator for fitted form or standard form
  *      = 0 when reaction is standard fit
  *      = 1 when reaction is form fit      */

int check_input(FILE* input);
// This function is to check whether or not the input file is a valid input file
/* Input parameters:
 *   1. input - the input file being examined   */
 /* Output parameters:
  *   1. (returned int) result - the indicator for fitted form or standard form
  *      = 0 when it is not a valid input file
  *      = 1 when it is a valid input file     */

int check_reaction_current_line(char Line_str[], int* way, int* tot, FILE* output, int loc, int run);
// This functions checks a line for a reaction and also the direction of the reaction
/* Input parameters:
 *   1. Line_str[] - the current line of input being searched for a reaction
 *   2. loc - the indicator for the location of the line within the input file being searched
 *   3. run - the indicator for which function is currently being built (here it is used to avoid printing an error during run = -1)
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx */
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

int check_M_reaction(char oldLine[], char prevLine[], char lastLine[], char newLine[], int tot, int loc, FILE* output, int* ply, int run, char* species);
// This function determines whether or not a reaction includes a third body and then determines how many steps are necessary to build each reaction substitution
/* Input parameters:
 *   1. oldLine[], prevLine[], lastLine[], newLine[] - the current reaction line being searched and the following lines of information from the input file   
 *   2. tot - the total number of reactions detected thus far
 *   3. loc - the indicator for the location of the line within the input file being searched
 *   4. ply - this hold the M indicator during the function and only indicates the presence of a third body reaction
 *   5. run - the indicator for which function is currently being built
 *    = -1 when building REACTION[]
 *    = 0 when building find_W_CHEMKIN()
 *    = 1 when building find_dW_dx_CHEMKIN ()
 *    = 2 when building find_Qei()
 *    = 3 when building find_dQei_dx 
 *   6. species - the list of species to be compared with known processes   */
 /* Output parameters:
  *   1. output file - if a reaction is ignored or unreadable an error is printed to the output file
  *   2. (returned int) M - the indicator for a reaction with a third body
  *       = -2 when there is no third body
  *       = -4 for the third body does not have specified element substitutions
  *       = int > 0 when there is a third body with specified elements/factors, the int represents how many different elements are present  */





int main(int argc, char **argv) {
  char *inpfile, *outfile;
  bool VALIDOPTIONS = TRUE;
  char *options;
  int RET, run, i;
  options = NULL;
  i = 1;
  
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
	
  if (check_input(input) == 0) {
    fprintf ( stderr, "\n\nThe following input file was not a valid input file. If the input file is valid, check for a missing SPECIES list \n%s\n\n", inpfile );
    exit (EXIT_FAILURE);
  }
  FILE* output = fopen(outfile, "w");
	for (run = -1; run != 4; run++) {              // Begin cycle for creating functions in output file

		print_function_header(run, output);
		i = build_current_function(run, input, output);
    if (i == 0) {
      fprintf(stderr, "\n\nERROR\nThe following input file does not provide a species list:\n%s\n\n", inpfile);
      exit (EXIT_FAILURE);
    }
    if (run != -1)
      fprintf(output, "\n}\n\n");
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
            else if (Line_str[i + 1] == ' ' || Line_str[i + 1] == '\t') {
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
            else if (Line_str[i + 1] == ' ' || Line_str[i + 1] == '\t') {
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
        fprintf(output, "void find_W_CHEMKIN ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {\n  double eff;\n  long k, specM;\n  spec_t X;\n\n  for ( k = 0; k < ns; k++ ) {\n    X[k] = rhok[k] / _calM (k) * 1.0e-6;    /* mole/cm^3 */\n    W[k] = 0.0;\n  }\n");
    }
    if (run == 1) {
        fprintf(output, "void find_dW_dx_CHEMKIN ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdT, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {\n  long s, k, specM;\n  double eff;\n  spec_t X;\n\n  for ( k = 0; k < ns; k++ ) {\n    X[k] = rhok[k] / _calM (k) * 1.0e-6;    /* mole/cm^3 */\n  }\n\n  for ( s = 0; s < ns; s++ ) {\n    dWdT[s] = 0.0;\n    dWdTe[s] = 0.0;\n    dWdTv[s] = 0.0;\n    dWdQbeam[s] = 0.0;\n    for ( k = 0; k < ns; k++) {\n      dWdrhok[s][k] = 0.0;\n    }\n  }\n");
    }
    if (run == 2) {
        fprintf(output, "void find_Qei_CHEMKIN(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei) {\n  long specM;\n\n  *Qei = 0.0;\n\n");
    }
    if (run == 3) {
        fprintf(output, "void find_dQei_dx_CHEMKIN(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe) {\n  long spec, specM;\n\n  for ( spec = 0; spec < ns; spec++ )\n    dQeidrhok[spec] = 0.0;\n\n  *dQeidTe = 0.0;\n\n");
    }
}

int check_Q(char line[]) {
    int i = 0, ind = 0;                 // searches line for "EXCI"

    for (i = 0; line[i] != '\0'; i++) {
        if (i != 0 && line[i-1] == '!' && line[i] == 'E' && line[i + 1] == 'X' && line[i + 2] == 'C' && line[i + 3] == 'I') {
            ind = 2;
            break;
          }
        if (line[i] == '!' && line[i+1] != 'E') {
            ind = 0;
            break;
          }
        if (line[i] == 'E' && line[i + 1] == 'X' && line[i + 2] == 'C' && line[i + 3] == 'I') {
            ind = 1;
            break;
          }
    }
    return ind;             // returns 1 if "EXCI" is present, returns 0 if not, returns 2 for not provided EXCI value
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

void print_info_find_Qei(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, char newLine[], char lastLine[], int Q, int e, int run, int reaction, int tot) {
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
    char products[10][12]; // create a string array to hold product names
    i = 0, j = 0;
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
            k = 0;
            j++;
        }
    }
    
    if (run == 2) {
      fprintf(output, "  if (REACTION[%d])\n", tot);
      fprintf(output, "    add_to_Qei(spec");                  // printing find_Qei functions
    }
    else if (run == 3) {
      fprintf(output, "  if (REACTION[%d])\n", tot);
      fprintf(output, "    add_to_dQei(spec");                  // printing find_Qei functions
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
    if (Q == -2) {
      fprintf(output, "_exci_%dr%dp(spec", numR, numP);
      for (j = 0; j < numR; j++) {
        i = 0;
        while (reactants[j][i] != '\0') {
          fprintf(output, "%c", reactants[j][i]);
          i++;
        }
        fprintf(output, ", spec");
      }
      for (j = 0; j < numP; j++) {
        i = 0;
        while (products[j][i] != '\0') {
          fprintf(output, "%c", products[j][i]);
          i++;
        }
        if (j != numP - 1)
          fprintf(output, ", spec");
        else if (j == numP - 1)
          fprintf(output, "), ");
        }
    }
    if (ind == 1) { // printing fitted form numbers
        fprintf(output, "_kf_fit4(");
        print_numbers_find_Qei(oldLine, prevLine, numR, numP, ind, output, run, reaction, e);
    }
    else if (ind == 0) { // printing standard form numbers
        fprintf(output, "_kf(");
        print_numbers_find_Qei(oldLine, prevLine, numR, numP, ind, output, run, reaction, e);
    }
}

void print_numbers_chemkin(char oldLine[], char prevLine[], char lastLine[], char newLine[], int numR, int numP, int ind, FILE* output, int step, int M, int reaction, int w) {
    char numsF[6][12], numsS[3][12], extra[3][12];
    int pow = 0;
    if (w != 0) {
        int i = 0, j = 0, k = 0;
        if (w == -1 || w == 1) {
            while (prevLine[i] != '/')
                i++;
            for (j = 0; j < 3; j++) {
                k = 0;
                if (prevLine[i] == '/')
                    i++;
                while (prevLine[i] == ' ' || prevLine[i] == '\t' || prevLine[i] == '/')
                    i++;
                while (prevLine[i] != ' ' && prevLine[i] != '\t' && prevLine[i] != '/') {
                    extra[j][k] = prevLine[i];
                    i++, k++;
                }
                extra[j][k] = '\0';
                if (extra[j][k - 1] == '\n' || extra[j][k - 1] == '\r')
                    extra[j][k - 1] = '\0';
                if ((extra[j][0] == '0' && extra[j][1] == '\0') || (extra[j][0] == '.' && extra[j][1] == '0' && extra[j][2] == '0'))
                    extra[j][0] = '0', extra[j][1] = '.', extra[j][2] = '0', extra[j][3] = '\0';
                if (extra[j][1] == 'e' || extra[j][1] == 'E') {
                    for (int v = strlen(extra[j]); v != 0; v--)
                        extra[j][v + 2] = extra[j][v];
                    extra[j][2] = '0', extra[j][1] = '.';
                }
                int v = strlen(extra[j]), dec = 0;
                for (int g = 0; g != v; g++)
                    if (extra[j][g] == '.')
                        dec = 1;
                if (dec == 0)
                    extra[j][v] = '.', extra[j][v + 1] = '0', extra[j][v + 2] = '\0';
            }
        }
        if (w == -2 || w == 2) {
            while (lastLine[i] != '/')
                i++;
            for (j = 0; j < 3; j++) {
                k = 0;
                if (lastLine[i] == '/')
                    i++;
                while (lastLine[i] == ' ' || lastLine[i] == '\t' || lastLine[i] == '/')
                    i++;
                while (lastLine[i] != ' ' && lastLine[i] != '\t' && lastLine[i] != '/') {
                    extra[j][k] = lastLine[i];
                    i++, k++;
                }
                extra[j][k] = '\0';
                if (extra[j][k - 1] == '\n' || extra[j][k - 1] == '\r')
                    extra[j][k - 1] = '\0';
                if ((extra[j][0] == '0' && extra[j][1] == '\0') || (extra[j][0] == '.' && extra[j][1] == '0' && extra[j][2] == '0'))
                    extra[j][0] = '0', extra[j][1] = '.', extra[j][2] = '0', extra[j][3] = '\0';
                if (extra[j][1] == 'e' || extra[j][1] == 'E') {
                    for (int v = strlen(extra[j]); v != 0; v--)
                        extra[j][v + 2] = extra[j][v];
                    extra[j][2] = '0', extra[j][1] = '.';
                }
                int v = strlen(extra[j]), dec = 0;
                for (int g = 0; g != v; g++)
                    if (extra[j][g] == '.')
                        dec = 1;
                if (dec == 0)
                    extra[j][v] = '.', extra[j][v + 1] = '0', extra[j][v + 2] = '\0';
            }
        }
        if (w == -3 || w == 3) {
            while (newLine[i] != '/')
                i++;
            for (j = 0; j < 3; j++) {
                k = 0;
                if (newLine[i] == '/')
                    i++;
                while (newLine[i] == ' ' || newLine[i] == '\t' || newLine[i] == '/')
                    i++;
                while (newLine[i] != ' ' && newLine[i] != '\t' && newLine[i] != '/') {
                    extra[j][k] = newLine[i];
                    i++, k++;
                }
                extra[j][k] = '\0';
                if (extra[j][k - 1] == '\n' || extra[j][k - 1] == '\r')
                    extra[j][k - 1] = '\0';
                if ((extra[j][0] == '0' && extra[j][1] == '\0') || (extra[j][0] == '.' && extra[j][1] == '0' && extra[j][2] == '0'))
                    extra[j][0] = '0', extra[j][1] = '.', extra[j][2] = '0', extra[j][3] = '\0';
                if (extra[j][1] == 'e' || extra[j][1] == 'E') {
                    for (int v = strlen(extra[j]); v != 0; v--)
                        extra[j][v + 2] = extra[j][v];
                    extra[j][2] = '0', extra[j][1] = '.';
                }
                int v = strlen(extra[j]), dec = 0;
                for (int g = 0; g != v; g++)
                    if (extra[j][g] == '.')
                        dec = 1;
                if (dec == 0)
                    extra[j][v] = '.', extra[j][v + 1] = '0', extra[j][v + 2] = '\0';
            }
        }

    }
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
        if (M >= 0)
            fprintf(output, "eff*");
        for (j = 0; j < 6; j++) {
            pow++;
            for (i = 0; numsF[j][i] != '\0'; i++) {   // print fitted form numbers
                fprintf(output, "%c", numsF[j][i]);
            }
            if (j == 0 && (reaction == 2 || reaction == 3))
                fprintf(output, "*powint(calA, %d), ", numR - 1);
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
        int i = 0, j = 0, k = 0, pow = 0;
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
        if (w > 0) {
            j = 0, pow = 0;
            if (M >= 0)
                fprintf(output, "eff*");
            for (j = 0; j < 3; j++) {
                pow++;
                for (i = 0; extra[j][i] != '\0'; i++)
                    fprintf(output, "%c", extra[j][i]);           // print standard form numbers
                if (j == 0 && (reaction == 2 || reaction == 3))
                    fprintf(output, "*powint(calA, %d), ", numR - 1);
                else if (j == 2 && (reaction == 1 || reaction == 3))
                    fprintf(output, "*Rchem, ");
                else if (j > 2 && (reaction == 1 || reaction == 3)) {
                    fprintf(output, "*powint(Rchem, %d), ", pow);
                }
                else
                    fprintf(output, ", ");
            }
        }
        j = 0, pow = 0;
        if (M >= 0)
            fprintf(output, "eff*");
        for (j = 0; j < 3; j++) {
            pow++;
            for (i = 0; numsS[j][i] != '\0'; i++)
                fprintf(output, "%c", numsS[j][i]);           // print standard form numbers
            if (j == 0 && (reaction == 2 || reaction == 3))
                fprintf(output, "*powint(calA, %d), ", numR - 1);
            else if (j == 2 && (reaction == 1 || reaction == 3))
                fprintf(output, "*Rchem, ");
            else if (j > 2 && (reaction == 1 || reaction == 3)) {
                fprintf(output, "*powint(Rchem, %d), ", pow);
            }
            else
                fprintf(output, ", ");
        }
        if (w < 0) {
            j = 0, pow = 0;
            if (M >= 0)
                fprintf(output, "eff*");
            for (j = 0; j < 3; j++) {
                pow++;
                for (i = 0; extra[j][i] != '\0'; i++)
                    fprintf(output, "%c", extra[j][i]);           // print standard form numbers
                if (j == 0 && (reaction == 2 || reaction == 3))
                    fprintf(output, "*powint(calA, %d), ", numR - 1);
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
                fprintf(output, "*powint(calA, %d), ", numR - 1);
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
            fprintf(output, ", Te)/powint(calA, %d)", numR - 1);
        else
            fprintf(output, ", T)/powint(calA, %d)", numR - 1);
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
                fprintf(output, "*powint(calA, %d), ", numR - 1);
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
            fprintf(output, ", Te)/powint(calA, %d)", numR - 1);
        else
            fprintf(output, ", T)/powint(calA, %d)", numR - 1);
        if (run == 2 && j == 3)
            fprintf(output, ", rhok, Qei);\n");
        else if (run == 3 && j == 3)
            fprintf(output, ", 0.0, rhok, dQeidrhok, dQeidTe);\n");
    }
    fprintf(output, "\n");
}

void print_elements_chemkin(char oldLine[], char prevLine[], char lastLine[], char newLine[], FILE* output, int numR, int numP, int step, int M, int ind, int run, int way, int tot, int loc, int w, char* species) {
    int i = 0, j = 0, k = 0, n = 0, m = 0, b = 0;
    char reactants[10][12], factor[10][15], body[10][15]; // create a string array to hold reactant names and factor info for third bodies
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
        }
    }
    i = 0, j = 0, k = 0, b = 0;
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
            k = 0, n = 0;
            j++;
        }
    }
    if (ind == 1 || ind == 0) { // builds an array to store the numbers of a third body reaction
        int k = 0, capp = M;
        n = 0;
        if (M >= 0 || M == -5) {
            while (capp >= 0) {
                if (check_process(prevLine, tot, loc, output, 0, run, species) == 10) {

                    while (prevLine[n] != '/' && M >= 0) {
                        while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\n')
                            n++;
                        if (prevLine[n] == '/')
                            break;
                        body[b][k] = prevLine[n];
                        k++;
                        n++;
                        while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\n')
                            n++;
                        if (prevLine[n] == '/')
                            break;
                    }
                    n++;
                    body[b][k] = '\0';
                    k = 0;
                    while (prevLine[n] != '/' && M >= 0) {
                        while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\n')
                            n++;
                        if (prevLine[n] == '/')
                            break;
                        factor[b][k] = prevLine[n];
                        k++;
                        n++;
                        while (prevLine[n] == ' ' || prevLine[n] == '\t' || prevLine[n] == '\n')
                            n++;
                        if (prevLine[n] == '/')
                            break;
                    }
                    n++;
                    factor[b][k] = '\0';
                    capp--, k = 0;
                    b++;
                }
                else if (check_process(lastLine, tot, loc, output, 0, run, species) == 10) {

                    while (lastLine[n] != '/' && M >= 0) {
                        while (lastLine[n] == ' ' || lastLine[n] == '\t' || lastLine[n] == '\n')
                            n++;
                        if (lastLine[n] == '/')
                            break;
                        body[b][k] = lastLine[n];
                        k++;
                        n++;
                        while (lastLine[n] == ' ' || lastLine[n] == '\t' || lastLine[n] == '\n')
                            n++;
                        if (lastLine[n] == '/')
                            break;
                    }
                    n++;
                    body[b][k] = '\0';
                    k = 0;
                    while (lastLine[n] != '/' && M >= 0) {
                        while (lastLine[n] == ' ' || lastLine[n] == '\t' || lastLine[n] == '\n')
                            n++;
                        if (lastLine[n] == '/')
                            break;
                        factor[b][k] = lastLine[n];
                        k++;
                        n++;
                        while (lastLine[n] == ' ' || lastLine[n] == '\t' || lastLine[n] == '\n')
                            n++;
                        if (lastLine[n] == '/')
                            break;
                    }
                    n++;
                    factor[b][k] = '\0';
                    capp--, k = 0;
                    b++;
                }
                else if (check_process(newLine, tot, loc, output, 0, run, species) == 10) {

                    while (newLine[n] != '/' && M >= 0) {
                        while (newLine[n] == ' ' || newLine[n] == '\t' || newLine[n] == '\n')
                            n++;
                        if (newLine[n] == '/')
                            break;
                        body[b][k] = newLine[n];
                        k++;
                        n++;
                        while (newLine[n] == ' ' || newLine[n] == '\t' || newLine[n] == '\n')
                            n++;
                        if (newLine[n] == '/')
                            break;
                    }
                    n++;
                    body[b][k] = '\0';
                    k = 0;
                    while (newLine[n] != '/' && M >= 0) {
                        while (newLine[n] == ' ' || newLine[n] == '\t' || newLine[n] == '\n')
                            n++;
                        if (newLine[n] == '/')
                            break;
                        factor[b][k] = newLine[n];
                        k++;
                        n++;
                        while (newLine[n] == ' ' || newLine[n] == '\t' || newLine[n] == '\n')
                            n++;
                        if (newLine[n] == '/')
                            break;
                    }
                    n++;
                    factor[b][k] = '\0';
                    capp--, k = 0;
                    b++;
                }
            }
        }
    }
    int cap = b - 1;
    b = 0, k = 0;
    if (M >= 0) {
        while (cap >= 0) {
            k = 0;
            fprintf(output, "\n        case spec");
            while (body[b][k] != '\0') {
                fprintf(output, "%c", body[b][k]);
                k++;
            }
            k = 0;
            fprintf(output, ": eff = ");
            while (factor[b][k] != '\0') {
                fprintf(output, "%c", factor[b][k]);
                k++;
            }

            fprintf(output, "; break;");
            cap--, b++;;
        }
        fprintf(output, "\n        default: eff = 1.0;\n      }\n");
    }
    for (j = 0; j < numR; j++) { // printing reaction names
        if (j == 0) {
            if (M == -5 || M >= 0) {
                if (run == 0) {
                    if (way == 2)
                        fprintf(output, "      add_to_W_fwbw_%dr%dp", numR, numP);
                    else if (way == 1)
                        fprintf(output, "      add_to_W_fw_%dr%dp", numR, numP);
                    if (ind == 1 && w != 0)
                        fprintf(output, "_fit4_Lindemann("); // printing fit information
                    else if (w != 0)
                        fprintf(output, "_Lindemann(");
                    else
                        fprintf(output, "(");
                }
                if (run == 1) {
                    if (way == 2)
                        fprintf(output, "      add_to_dW_fwbw_%dr%dp", numR, numP);
                    else if (way == 1)
                        fprintf(output, "      add_to_dW_fw_%dr%dp", numR, numP);
                    if (ind == 1 && w != 0)
                        fprintf(output, "_fit4_Lindemann("); // printing fit information
                    else if (w != 0)
                        fprintf(output, "_Lindemann(");
                    else
                        fprintf(output, "(");
                }
            }
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


void print_output_line(char oldLine[], char prevLine[], char lastLine[], char newLine[], int ind, int e, int Q, int way, int tot, int run, int numR, int numP, FILE* output, int M, int reaction, int h, int loc, int ply, char* species) {

    int step = 0;
    if (M > 0 || M == -2 || M == -4) {
        step++;
        M = M - 1;
        if (M == -3 || M == -5)
            step = 0;
        if (run == 0) { // Begin printing first function
            if (step == 0 || step == 1) {
                fprintf(output, "\n  if (REACTION[%d])", tot); // begin printing information in output file
            }
            if (step == 1 || M == -5)
                fprintf(output, "  {\n");
            else
                fprintf(output, "\n");
            if (M >= 0 && ply == 1)
                fprintf(output, "    for (specM = 0; specM < ns; specM++) {\n      switch (specM) {");
            else if (M == -5 && ply == 1)
                fprintf(output, "    for (specM = 0; specM < ns; specM++) {\n");
            if (M >= 0 && ply == -1)
                fprintf(output, "    for (spec = 0; spec < ns; spec++) {\n      switch (specM) {");
            else if (M == -5 && ply == -1)
                fprintf(output, "    for (spec = 0; spec < ns; spec++) {\n");
            if (way == 2 && M == -3)
                fprintf(output, "    add_to_W_fwbw_%dr%dp", numR, numP);
            else if (way == 1 && M == -3)
                fprintf(output, "    add_to_W_fw_%dr%dp", numR, numP);
            if (ind == 1 && M == -3 && h == 0)
                fprintf(output, "_fit4("); // printing fit information
            else if (M == -3 && h == 0)
                fprintf(output, "(");
            else if (ind == 1 && M == -3 && h != 0)
                fprintf(output, "_fit4_Lindemann(");
            else if (M == -3 && h != 0)
                fprintf(output, "_Lindemann(");
            print_elements_chemkin(oldLine, prevLine, lastLine, newLine, output, numR, numP, step, M, ind, run, way, tot, loc, h, species);
            step++;
            if (M == -3 || M == -5)
                step = 0;
            print_numbers_chemkin(oldLine, prevLine, lastLine, newLine, numR, numP, ind, output, step, M, reaction, h);
            if (e == 1) // printing temperature Te vs T information
                fprintf(output, "Te, X, W);");
            else
                fprintf(output, "T, X, W);");
            if (M == -3)
                fprintf(output, "\n");
        }
        else if (run == 1) {
            if (step == 0 || step == 1) {
                fprintf(output, "\n  if (REACTION[%d])", tot); // begin printing information in output file
            }
            if (step == 1 || M == -5)
                fprintf(output, "  {\n");
            else
                fprintf(output, "\n");
            if (M >= 0 && ply == 1)
                fprintf(output, "    for (specM = 0; specM < ns; specM++) {\n      switch (specM) {");
            else if (M == -5 && ply == 1)
                fprintf(output, "    for (specM = 0; specM < ns; specM++) {\n");
            if (M >= 0 && ply == -1)
                fprintf(output, "    for (spec = 0; spec < ns; spec++) {\n      switch (specM) {");
            else if (M == -5 && ply == -1)
                fprintf(output, "    for (spec = 0; spec < ns; spec++) {\n");
            if (way == 2 && M == -3)
                fprintf(output, "    add_to_dW_fwbw_%dr%dp", numR, numP);
            else if (way == 1 && M == -3)
                fprintf(output, "    add_to_dW_fw_%dr%dp", numR, numP);
            if (ind == 1 && M == -3 && h == 0)
                fprintf(output, "_fit4("); // printing fit information
            else if (M == -3 && h == 0)
                fprintf(output, "(");
            else if (ind == 1 && M == -3 && h != 0)
                fprintf(output, "_fit4_Lindemann(");
            else if (M == -3 && h != 0)
                fprintf(output, "_Lindemann(");
            print_elements_chemkin(oldLine, prevLine, lastLine, newLine, output, numR, numP, step, M, ind, run, way, tot, loc, h, species);
            step++;
            if (M == -3 || M == -5)
                step = 0;
            print_numbers_chemkin(oldLine, prevLine, lastLine, newLine, numR, numP, ind, output, step, M, reaction, h);
            if (e == 1) // printing temperature Te vs T information
                fprintf(output, "Te, X, dWdTe, dWdrhok);");
            else
                fprintf(output, "T, X, dWdT, dWdrhok);");
            if (M == -3)
                fprintf(output, "\n");
        }
        else if (run == 2 && Q != 0) {
            print_info_find_Qei(oldLine, prevLine, numR, numP, ind, output, newLine, lastLine, Q, e, run, reaction, tot);

        }
        else if (run == 3 && Q != 0) {
            print_info_find_Qei(oldLine, prevLine, numR, numP, ind, output, newLine, lastLine, Q, e, run, reaction, tot);

        }

    }
    if (run == 0 || run == 1) {
        if (M != -3)
            fprintf(output, "\n    }\n  }\n");
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
    if (*ind == 1 && *Q == 0)
      *Q = -2;                              // For fit reactions there should be a EXCI value, thus *Q = -2 if there is none provided or commented out
}

int build_current_function(int run, FILE* input, FILE* output) {
        
    char *newLine, *oldLine, *prevLine, *lastLine, *scoutLine,  *word, *species;
    
    species = (char*)malloc(sizeof(char)*1);
    newLine = (char *)malloc(1000*sizeof(char));
    oldLine = (char *)malloc(1000*sizeof(char));
    prevLine = (char *)malloc(1000*sizeof(char));
    lastLine = (char *)malloc(1000*sizeof(char));
    scoutLine = (char *)malloc(1000*sizeof(char));
    word = (char *)malloc(200*sizeof(char));
    
    
    int loc = 0, numR = 0, numP = 0, tot = 0, ind = 0, check = 0, e = 0, way = 10, Q = 0, M = 0, flag1 = 0, flag2 = 0, flag3 = 0, flag4 = 0, reaction = 0, h = 0, ply = 0, i = 0, marker = 0;

    rewind(input);
    fgets(oldLine, 1000, input);
    fgets(prevLine, 1000, input);
    fgets(lastLine, 1000, input);
    fgets(newLine, 1000, input);
    
    oldLine = (char *)realloc(oldLine, (strlen(oldLine) + 1) * sizeof(char));
    prevLine = (char *)realloc(prevLine, (strlen(prevLine) + 1) * sizeof(char));    
    lastLine = (char *)realloc(lastLine, (strlen(lastLine) + 1) * sizeof(char));
    newLine = (char *)realloc(newLine, (strlen(newLine) + 1) * sizeof(char));    

    species[0] = 'L';
    
    while (marker <= 3) {
        loc++;
        flag1 = 0, flag2 = 0, flag3 = 0, flag4 = 0;
        if (species[0] == 'L') {
          for (i = 0; newLine[i] != '\r' && newLine[i] != '\n' && newLine[i] != ' '; i++) {
            word[i] = newLine[i];
          }
          word[i] = '\0';
          if (strcmp(word, "SPECIES") == 0) {
            free (species);
            write_species(input, &species, &loc);
          }
        }
        reaction = reaction_class(oldLine, reaction);
        Q = 0, way = 10, e = 0, check = 0, ind = 0, h = 0;
        check = check_reaction_current_line(oldLine, &way, &tot, output, loc, run);
        if (check == 1) {
            flag1 = check_process(oldLine, tot, loc, output, 1, run, species);
            flag2 = check_process(prevLine, tot, loc, output, 1, run, species);
            if (flag2 != 3 && flag2 != 1) { 
                flag3 = check_process(lastLine, tot, loc, output, 1, run, species); 
                if (flag3 != 3 && flag3 != 1) { 
                    flag4 = check_process(newLine, tot, loc, output, 1, run, species);
                  }
            }
            if (flag1 == 1 || flag2 == 1 || flag3 == 1 || flag4 == 1) 
                check = 0;
            if (flag2 == 2 || flag2 == -2)
                h = flag2 / 2;
            if (flag3 == 2 || flag3 == -2)
                h = flag3;
            if (flag4 == 2 || flag4 == -2)
                h = flag4 / 2 * 3;
        }
        if (check == 1) {
            check_all_indicators_formfit_electrontemperature_Q(oldLine, prevLine, lastLine, newLine, &ind, &e, &Q);
            M = check_M_reaction(oldLine, prevLine, lastLine, newLine, tot, loc, output, &ply, run, species);
            numR = number_of_reactants(oldLine);
            numP = number_of_products(oldLine);
            if (run != -1)
              print_output_line(oldLine, prevLine, lastLine, newLine, ind, e, Q, way, tot, run, numR, numP, output, M, reaction, h, loc, ply, species);
        }
        if (marker < 4) {
          oldLine = (char *)realloc(oldLine, (strlen(prevLine) + 1) * sizeof(char));
          strcpy(oldLine, prevLine);
          if (marker == 3)
            marker++;
        }
        
        if (marker < 3) {
          prevLine = (char *)realloc(prevLine, (strlen(lastLine) + 1) * sizeof(char));
          strcpy(prevLine, lastLine);
          if (marker == 2)
            marker++;
        }
        
        if (marker < 2) {   
          lastLine = (char *)realloc(lastLine, (strlen(newLine) + 1) * sizeof(char));    
          strcpy(lastLine, newLine);
          if (marker == 1)
            marker++;
        }
        
        if (marker < 1) {
          newLine = (char *)realloc(newLine, (strlen(scoutLine) + 1) * sizeof(char));    
          strcpy(newLine, scoutLine);
        }
        
        
        if (marker == 0) {
          if (fgets(scoutLine, 1000, input) != NULL) {
            newLine = (char *)realloc(newLine, (strlen(scoutLine) + 1) * sizeof(char));    
            strcpy(newLine, scoutLine);
          }
          else {
            marker = 1;
          }
        }
        
        
    }
    if (species[0] == 'L') {
      free ( lastLine );
      free ( oldLine );
      free ( prevLine );
      free ( newLine );
      free ( scoutLine );
      free ( word );
      free ( species );
      return 0;
    }
    if (run == -1) {
      fprintf(output, "const static bool REACTION[%d] = {\n", tot+1);
      for (int num = 0; num < tot+1; num++) {
        fprintf(output,"  TRUE, /* reaction[%d] */\n", num);
      }
      fprintf(output, "};\n\n");
    }
    free ( lastLine );
    free ( oldLine );
    free ( prevLine );
    free ( newLine );
    free ( scoutLine );
    free ( word );
    free ( species );
    return 1;
}


int check_reaction_current_line(char Line_str[], int* way, int* tot, FILE* output, int loc, int run) {
    int i, num = 0;
    for (i = 0; Line_str[i] != '\n' && Line_str[i] != '\r'; i++) {              // determines the presence of a reaction by "=>" or "<=>"
        if (Line_str[i] == '=' && Line_str[i + 1] == '>' && Line_str[i - 1] == '<') {
            *tot = *tot + 1;
            num = 1;
            *way = 2;                                                 // way indicates whether the reaction is forwards or backwards
            break;
        }
        if (Line_str[i] == '=' && Line_str[i + 1] == '>') {
            *tot = *tot + 1;
            num = 1;
            *way = 1;
            break;
        }
    }
    if (Line_str[0] == '!' && num == 1) {
        num = 0;
        if (run != -1)
          fprintf(output, "\n        // ************ERROR: Reaction #%d ignored on line %d due to commenting out.\n\n", *tot, loc);
        return num;
    }
    if (num == 1)
        return num;
    else {

        for (i = 0; Line_str[i] != '\n' && Line_str[i] != '\r' && Line_str[i] != '!'; i++) {
            if (Line_str[i] == '=') {
                num = 1;
                *tot = *tot + 1;
                *way = 2;
            }
        }
        if (Line_str[0] == '!' && num == 1) {
            num = 0;
            *tot = *tot - 1;
            return num;
        }
    }


    return num; // returns 0 if there is no reaction and 1 if there is a reaction
}

int check_M_reaction(char oldLine[], char prevLine[], char lastLine[], char newLine[], int tot, int loc, FILE* output, int* ply, int run, char* species) {
    int i = 0, j = 0, k = 0, M = 0, m = 0;
    char reactants[10][12]; // create a string array to hold reactant names
    for (j = 0; oldLine[i] != '='; j++) {
        k = 0;
        if (oldLine[i] == '\t' || oldLine[i] == ' ' || oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '\n' || oldLine[i] == '\r') {
            j--, i++;
        }
        else {
            while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\r') {
                if (oldLine[i] == '(' || oldLine[i] == ')')
                    i++;
                if (k != 0 && oldLine[i] == '(')
                    break;
                if (k != 0 && oldLine[i] == '+' && oldLine[i - 1] != '(' && oldLine[i - 1] != ')')
                    break;
                if (oldLine[i] == '+' && oldLine[i - 1] != '(' && oldLine[i - 1] != ')' && oldLine[i + 1] != '(' && oldLine[i + 1] != ')' && oldLine[i] != '=')
                    break;
                else if (oldLine[i] == '<' || oldLine[i] == '>' || oldLine[i] == '=') {
                    break;
                }
                else {
                    reactants[j][k] = oldLine[i];
                    i = i + 1, k = k + 1;
                }
            }
            if (oldLine[i] == '+' && (oldLine[i + 1] == ' ' || oldLine[i + 1] == '\t' || oldLine[i + 1] == '(')) {
                m = i + 1;
                while (oldLine[m] == ' ' || oldLine[m] == '\t' || oldLine[m] == '(' || oldLine[m] == ')') {
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
            if (k != 0)
                reactants[j][k] = '\0';
            else
                j--, i++;
            if (reactants[j][0] == 'E' && reactants[j][1] == '\0')
                reactants[j][0] = 'e', reactants[j][1] = 'm', reactants[j][2] = 'i', reactants[j][3] = 'n', reactants[j][4] = 'u', reactants[j][5] = 's', reactants[j][6] = '\0';
            if (reactants[j][0] == '(' && reactants[j][1] == '+' && reactants[j][2] == 'M' && reactants[j][3] == ')')
                M = 1;
            if (reactants[j][0] == 'M' && (reactants[j][1] == '\0' || reactants[j][1] == ' '))              // determine if any include an M third body
                M = 1;
        }
    }
    if (M == 1 || M == -1) {
        int i = 0, count = 0;
        for (i = 0; prevLine[i] != '\n' && prevLine[i] != '\r'; i++) {
            if (prevLine[i] == '=') {
                *ply = M;
                M = -4;                                                 // if there is a third body but not enough info provided for substitutions
                break;
            }
        }
        if (M == 1 || M == -1) {
            if (check_process(prevLine, tot, loc, output, 0, run, species) == 10) {
                for (i = 0; prevLine[i] != '\n' && prevLine[i] != '\r'; i++) {
                    if (prevLine[i] == '/')
                        count = count + 1;
                }
                *ply = M;
                M = count / 2;                                    // count for how many substitutions are necessary if third body info is present
                    if (count == 0)
                        M = -4;
            }
            else if (check_process(lastLine, tot, loc, output, 0, run, species) == 10) {
                for (i = 0; lastLine[i] != '\n' && lastLine[i] != '\r'; i++) {
                    if (lastLine[i] == '/')
                        count = count + 1;
                }
                *ply = M;
                M = count / 2;                                       // count for how many substitutions are necessary if third body info is present
                    if (count == 0)
                        M = -4;
            }
            else if (check_process(newLine, tot, loc, output, 0, run, species) == 10) {
                for (i = 0; newLine[i] != '\n' && newLine[i] != '\r'; i++) {
                    if (newLine[i] == '/')
                        count = count + 1;
                }
                *ply = M;
                M = count / 2;                                    // count for how many substitutions are necessary if third body info is present
                if (count == 0)
                    M = -4;
            }
            else
                M = -4;
        }
    }
    else
        M = -2;
    return M;                   // returns number of substitutions necessary if third body is present, returns -2 if no third body, returns -4 if third body
                                // but not enough information

}

int check_process(char Line[], int tot, int loc, FILE* output, int num, int run, char* species) {
    char word[12];
    int i, k = 0;

    for (i = 0; Line[i] != '\n' && Line[i] != '\r'; i++) {              // determines the presence of a reaction by "=>" or "<=>"
        if (Line[i] == '=')
            return 3;
        if (Line[i] == '!')
            return 0;
    }
    i = 0;
    if (Line[i] == '!')
        return 0;
    while (Line[i] != '\n' && Line[i] != '\r' && Line[i] != '/') {
        while (Line[i] == ' ' || Line[i] == '\t')
            i++;
        if (Line[i] == '/' || Line[i] == '\n' || Line[i] == '\r')
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
    else if (word[i] == 'H' && word[i + 1] == 'I' && word[i + 2] == 'G' && word[i+3] == 'H')
        return 2;
    else if (word[i] == 'L' && word[i + 1] == 'O' && word[i + 2] == 'W')
        return -2;
    else if (word[i] == 'T' && word[i + 1] == 'R' && word[i + 2] == 'O' && word[i + 3] == 'E') {
        if (num == 1 && run != -1 && run < 2) {
            fprintf(output, "\n\n        // **********WARNING: Reaction #%d on line %d specifies the \"%s\" process but we here use the \"Lindemann\" process instead:\n", tot, loc, word);
        }
        return 0;
    }
    else if (word[i] == '\0')
        return 0;
    
    if (run != -1) {
      i = 0;
      int j = 0;
      char wordbank[30];
      while (species[j] != '\0') {
        i = 0;
        while (species[j] != ' ' && species[j] != '\0') {
          wordbank[i] = species[j];
          i++, j++;
        }
        if (species[j] == ' ') {
          wordbank[i] = '\0';
          j++;
        }
        if ((strcmp(word, wordbank) == 0) || (strcmp(word, "END") == 0))
          return 10;
      }
    }
    if (run != -1) {
      fprintf(output, "\n        // ************ERROR: Reaction #%d ignored on line %d due to unknown process or element: \"%s\".", tot, loc, word);
      fprintf(output, "\n        // ************ERROR: If the unknown operator above is an element, make sure that the species list is provided and correct.\n");
    }

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

void write_species(FILE* input, char** species, int* loc) {
  char line[500]; 
  int i = 0, j = 0;
  fgets(line, 500, input);
  *loc = *loc + 1;
  (*species) = (char*)malloc(1*sizeof(char));
  while (line[0] != 'E' || line[1] != 'N' || line[2] != 'D') {
    if (line[0] != '!' && line[0] != '\n' && line[0] != '\r') {
      while (line[i] != '\n' && line[i] != '\r') {
        while (line[i] != ' ' && line[i] != '\n' && line[i] != '\t' && line[i] != '\r') {
          (*species) = (char*) realloc(*species, (j+1)*sizeof(char));
          (*species)[j] = line[i];
          i++, j++;
        }
        if (i != 0) {
          if (line[i-1] != ' ' && line[i-1] != '\t') {
            (*species) = (char *) realloc(*species, (j+1)*sizeof(char));
            (*species)[j] = ' ';
            j++;
          }
        }
        if (line[i] == ' ' || line[i] == '\t') {
          while (line[i] == ' ' || line[i] == '\t')
            i++;
        }
      }
    }
    i = 0;
    fgets(line, 500, input);
    *loc = *loc + 1;
  }
  (*species) = (char *) realloc (*species, (j+1)*sizeof(char));
  (*species)[j] = '\0';                                       // Set the species pointer to the local variable and the returned address will reflect the list changes
}

int check_input(FILE* input) {
    char *newLine, *word, *species;
    int i = 0;
    int loc = 0;
    
    species = (char*)malloc(sizeof(char)*1);
    newLine = (char *)malloc(1000*sizeof(char));
    word = (char *)malloc(200*sizeof(char));
    
    rewind(input);
    
    
    species[0] = 'L';
    
    if (input == NULL)
      return 0;
    while (fgets(newLine, 1000, input) != NULL) {
        loc++;
        newLine = (char *)realloc(newLine, (strlen(newLine) + 1) * sizeof(char));
        if (species[0] == 'L') {
          for (i = 0; newLine[i] != '\r' && newLine[i] != '\n' && newLine[i] != ' '; i++) {
            word[i] = newLine[i];
          }
          word[i] = '\0';
          if (strcmp(word, "SPECIES") == 0) {
            free (species);
            write_species(input, &species, &loc);
          }
        }
        
        
        newLine = (char *)realloc(newLine, 1000 * sizeof(char));
    }
    if (species[0] == 'L') {
        free ( newLine );
        free ( word );
        free ( species );
        return 0;
      }
    free ( newLine );
    free ( word );
    free ( species );
    
    return 1;
}

