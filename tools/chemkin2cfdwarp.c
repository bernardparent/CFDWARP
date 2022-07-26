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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "share.h"

void make_function(int run, FILE* input, FILE* output);
void print_head(int run, FILE* output);
void check_above(char oldLine[], char prevLine[], char lastLine[], char newLine[], int* ind, int* e, int* Q);
void print_begin(char oldLine[], char prevLine[], char lastLine[], char newLine[], int ind, int e, int Q, int way, int tot, int run, int numR, int numP, FILE* output);
void print_elem0(char oldLine[], FILE* output, int numR, int numP);
void print_num0(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output);
void print_num1(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, int run);
void print_Q(char line[], FILE* output);
void print_Q2(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, char newLine[], char lastLine[], int Q, int e, int run);
int check_Q(char line[]);
int check_elec(char line[]);
int check_fit(char Line_str[]);
int search_current_line(char Line_str[], int* way, int* tot, FILE* output, int loc);
int number_of_reactants(char Line_str[]);
int number_of_products(char Line_str[]);
int process_flag_string(int argc, char **argv, char *flag, char** arg);
int find_remaining_options(int argc, char **argv, char **options);

int main(int argc, char **argv) {

  char *inpfile, *outfile;
  bool VALIDOPTIONS = TRUE;
  char *options;
  int RET, run;
  options = NULL;
  
  
   inpfile =(char *)malloc(400*sizeof(char));
   outfile =(char *)malloc(400*sizeof(char));
  
  if (process_flag_string(argc, argv, "-i", &inpfile)!=2) VALIDOPTIONS=FALSE;
  if (process_flag_string(argc, argv, "-o", &outfile)!=2) VALIDOPTIONS=FALSE;
  
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
  if (input == NULL) {
    fprintf ( stderr, "\n\nThe following input file could not be found:\n%s\n\n", inpfile );
    exit (EXIT_FAILURE);
  }
	FILE* output = fopen(outfile, "w");
  
  
	for (run = 0; run != 4; run++) {

		print_head(run, output);
		make_function(run, input, output);
		fprintf(output, "}\n\n");
    }
    
  fclose(input);
	fclose(output);
  }

	free ( options );
  return(EXIT_SUCCESS);
}





int number_of_reactants(char Line_str[]) {									// This function determines the reactant count
	int i = 0, j = 0, num = 1;
	for (j = 0; Line_str[i] != '='; j++) {
		if (Line_str[i] == '\t')
			j--, i++;
		else if (Line_str[i] == '<' || Line_str[i] == '>')
			j--, i++;
		else if (Line_str[i] == ' ')
			j--, i++;
		else if (Line_str[i] == '+' && Line_str[i - 1] == ' ')
			num++, i++, j--;
		else if (Line_str[i] == '(' || Line_str[i] == ')') {
			Line_str[i] = ' ';
			i++, j--;
		}
		else {
			while (Line_str[i] != ' ' && Line_str[i] != '\n') {
				if (Line_str[i] == '(' || Line_str[i] == ')')
					i++;
				else if (Line_str[i] == '<' || Line_str[i] == '>')
					break;
				else {
					i = i + 1;
				}
			}
		}
	}
	return num;																// returns the number of reactants
}

int number_of_products(char Line_str[]) {									// This function determines the product count
	int i = 0, j = 0, num = 1;

	while (Line_str[i] != '>')
		i++;
	i++;
	while (Line_str[i + 2] != '.' && Line_str[i + 2] != 'e') {
		if (Line_str[i] == '\t')
			i++;
		else if (Line_str[i] == ' ')
			i++;
		else if (Line_str[i] == '+' && Line_str[i - 1] == ' ')
			num++, i++;
		else if (Line_str[i] == '(' || Line_str[i] == ')') {
			Line_str[i] = ' ';
			i++;
		}
		else {
			while (Line_str[i] != ' ' && Line_str[i] != '\t') {
				if (Line_str[i] == '(' || Line_str[i] == ')')
					i++;
				else {
					i++;
				}
			}
			j++;
		}
	}
	return num;																// returns the number of products
}

int check_elec(char line[]) {												// This function checks if a reaction should use Te or T
	int i = 0, result = 0;
	for (i = 0; line[i] != '\n'; i++) {
		if (line[i] == 'T' && line[i + 1] == 'D' && line[i + 2] == 'E' && line[i + 3] == 'P' && line[i + 4] == '/' && line[i + 5] == 'E' && line[i + 6] == '/')
			result = 1;
	}
	return result;															// 1 is returned if Te, 0 is returned if T
}

int check_fit(char Line_str[]) {											// This function classifies a reaction between fitted and standard
	int i, num = 0;
	for (i = 0; Line_str[i] != '\n'; i++) {
		if (Line_str[i] == 'F' && Line_str[i + 1] == 'I' && Line_str[i + 2] == 'T') {
			num = 1;
		}
	}
	return num;																// returns 0 if the reaction is a standard form, 1 if it is fitted form
}
void print_head(int run, FILE* output) {
	if (run == 0) {
		fprintf(output, "void add_W_Ionization_Adamovich ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec_t W ) {\n  double N[ns];\n  double R;\n  long k;\n  spec_t X;\n\n  for ( k = 0; k < ns; k++ ) {\n  \tN[k] = rhok[k] / _calM (k) * 1.0e-6 * calA;\t\t/* particules/cm^3 */\n  \tX[k] = rhok[k] / _calM (k) * 1.0e-6;\t\t/* mole/cm^3 */\n  }\n  R = 1.9872;\n\n");
	}
	if (run == 1) {
		fprintf(output, "/* Verify the validity of the dW terms at node i=10, j=10 using the command ./test -r control.wrp -node 10 10 dSchemdU\n * Make sure to verify the dW terms over a wide range of temperatures and mass fractions\n * Note that the verification using ./test is done by comparing the analytical expressions to numerical derivatives\n * The numerical derivatives depend strongly on the values given to Uref[] within Cycle()\n */\n\nvoid add_dW_dx_Ionization_Adamovich ( gl_t *gl, spec_t rhok, double T, double Te, double Tv, double Estar, double Qbeam, spec2_t dWdrhok, spec_t dWdTe, spec_t dWdTv, spec_t dWdQbeam ) {\n  long k, s;\n  spec_t N;\n  double R, kf, dkfdTe, dkfdT, dkfdTv;\n  spec_t X;\n\n  for ( k = 0; k < ns; k++ ) {\n  \tN[k] = rhok[k] / _calM (k) * 1.0e-6 * calA;\t\t/* particules/cm^3 */\n  \tX[k] = rhok[k] / _calM (k) * 1.0e-6;\t\t/* mole/cm^3 */\n  }\n\n\n  /* find properties needed by add_to_dW* functions in proper units */\n  for ( k = 0; k < ns; k++ ) {\n  \tN[k] = rhok[k] / _calM (k) * 1.0e-6 * calA;\n  }\n\n  R = 1.9872;\n\n\n");
	}
	if (run == 2) {
		fprintf(output, "void find_Qei(gl_t *gl, spec_t rhok, double Estar, double Te, double *Qei) {\n  double theta;\n\n  *Qei = 0.0;\n\n  switch (gl->model.chem.IONIZATIONMODEL) {\n  \tcase IONIZATIONMODEL_ADAMOVICH:\n");
	}
	if (run == 3) {
		fprintf(output, "void find_dQei_dx(gl_t *gl, spec_t rhok, double Estar, double Te, spec_t dQeidrhok, double *dQeidTe) {\n  double theta;\n  long spec;\n\n  for ( spec = 0; spec < ns; spec++ )\n  \tdQeidrhok[spec] = 0.0;\n\n  *dQeidTe = 0.0;\n\n  switch (gl->model.chem.IONIZATIONMODEL) {\n  \tcase IONIZATIONMODEL_ADAMOVICH:\n");
	}
}

int check_Q(char line[]) {
	int i = 0, ind = 0;

	for (i = 0; line[i] != '\0'; i++) {
		if (line[i] == 'E' && line[i + 1] == 'X' && line[i + 2] == 'C' && line[i + 3] == 'I') {
			ind = 1;
			break;
		}
	}
	return ind;
}

void print_Q(char line[], FILE* output) {
	int i = 0, k = 0, j = 0;
	char Qei[10];									// initialize Qei[]
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
	for (j = 0; Qei[j] != '\0'; j++) {					// print Qei[]
		fprintf(output, "%c", Qei[j]);
	}
	fprintf(output, ", ");
}

void print_Q2(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, char newLine[], char lastLine[], int Q, int e, int run) {
	int i = 0, j = 0, k = 0;
	char reactants[10][12];										// create a string array to hold reactant names
	for (j = 0; oldLine[i] != '='; j++) {
		k = 0;
		if (oldLine[i] == '\t' || oldLine[i] == ' ' || oldLine[i] == '<' || oldLine[i] == '>')
			j--, i++;
		else if (oldLine[i] == '+' && oldLine[i - 1] == ' ')
			i++, j--;
		else if (oldLine[i] == '(' || oldLine[i] == ')') {
			oldLine[i] = ' ';
			i++, j--;
		}
		else {
			while (oldLine[i] != ' ' && oldLine[i] != '\n') {
				if (oldLine[i] == '(' || oldLine[i] == ')')
					i++;
				else if (oldLine[i] == '<' || oldLine[i] == '>')
					break;
				else {
					reactants[j][k] = oldLine[i];
					i = i + 1, k = k + 1;
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
		if (reactants[j][0] == 'e' && reactants[j][1] == 'm')						// print primary element
			j++;
		if (j < numR)
			for (i = 0; reactants[j][i] != '\0'; i++)
				fprintf(output, "%c", reactants[j][i]);
	}
	fprintf(output, ", ");
	if (Q == 1)
		print_Q(newLine, output);
	if (Q == 2)
		print_Q(prevLine, output);
	if (Q == 3)
		print_Q(lastLine, output);
	if (ind == 1) {												// printing fitted form numbers
		fprintf(output, "_kfit4(");
		if (e == 1)
			fprintf(output, "Te, ");
		else
			fprintf(output, "T, ");
		print_num1(oldLine, prevLine, numR, numP, ind, output, run);
	}
	else if (ind == 0) {											// printing standard form numbers
		fprintf(output, "_k(");
		if (e == 1)
			fprintf(output, "Te, ");
		else
			fprintf(output, "T, ");
		print_num1(oldLine, prevLine, numR, numP, ind, output, run);
	}
}

void print_num0(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output) {
	char numsF[6][12], numsS[3][12];
	if (ind == 1) {										// builds an array to store the numbers of a fitted form reaction
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
				while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\t' && oldLine[i] != '/') {
					numsF[j][k] = oldLine[i];
					if (oldLine[i] == '-' && oldLine[i + 1] == ' ')
						i++;
					i++, k++;
				}
				if (k != 0)
					numsF[j][k] = '\0';
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
				while (prevLine[h] != ' ' && prevLine[h] != '\n' && prevLine[h] != '\t' && prevLine[h] != '/') {
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
		for (j = 0; j < 6; j++) {
			for (i = 0; numsF[j][i] != '\0'; i++) {
				fprintf(output, "%c", numsF[j][i]);
			}
			if (j == 0)
				fprintf(output, "*calA, ");
			else
				fprintf(output, ", ");
		}
	}
	else {												// builds an array to store the numbers of a standard form reaction
		int	i = 0, j = 0, k = 0;
		while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e')
			i++;
		for (j = 0; j < 3; j++) {
			k = 0;
			while (oldLine[i] == ' ' || oldLine[i] == '\t' || oldLine[i] == '/')
				i++;
			while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\t' && oldLine[i] != '/') {
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
		for (j = 0; j < 3; j++) {
			for (i = 0; numsS[j][i] != '\0'; i++)
				fprintf(output, "%c", numsS[j][i]);
			fprintf(output, ", ");
		}
	}
}

void print_num1(char oldLine[], char prevLine[], int numR, int numP, int ind, FILE* output, int run) {
	char numsF[6][12], numsS[3][12];
	if (ind == 1) {										// builds an array to store the numbers of a fitted form reaction
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
				while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\t' && oldLine[i] != '/') {
					numsF[j][k] = oldLine[i];
					if (oldLine[i] == '-' && oldLine[i + 1] == ' ')
						i++;
					i++, k++;
				}
				if (k != 0)
					numsF[j][k] = '\0';
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
				while (prevLine[h] != ' ' && prevLine[h] != '\n' && prevLine[h] != '\t' && prevLine[h] != '/') {
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
		for (j = 0; j < 6; j++) {
			for (i = 0; numsF[j][i] != '\0'; i++) {
				fprintf(output, "%c", numsF[j][i]);
			}
			if (j != 5)
				fprintf(output, ", ");
			else if (run == 2)
				fprintf(output, "), rhok, Qei);\n");
			else if (run == 3)
				fprintf(output, "), 0.0, rhok, dQeidrhok, dQeidTe);\n");
		}
	}
	else {												// builds an array to store the numbers of a standard form reaction
		int	i = 0, j = 0, k = 0;
		while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e')
			i++;
		for (j = 0; j < 3; j++) {
			k = 0;
			while (oldLine[i] == ' ' || oldLine[i] == '\t' || oldLine[i] == '/')
				i++;
			while (oldLine[i] != ' ' && oldLine[i] != '\n' && oldLine[i] != '\t' && oldLine[i] != '/') {
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
		for (j = 0; j < 3; j++) {
			for (i = 0; numsS[j][i] != '\0'; i++)
				fprintf(output, "%c", numsS[j][i]);
			if (j != 5)
				fprintf(output, ", ");
			else if (run == 2)
				fprintf(output, "), rhok, Qei);\n");
			else if (run == 3)
				fprintf(output, "), 0.0, rhok, dQeidrhok, dQeidTe);\n");
		}
	}
	fprintf(output, "\n");
}

void print_elem0(char oldLine[], FILE* output, int numR, int numP) {
	int i = 0, j = 0, k = 0;
	char reactants[10][12];										// create a string array to hold reactant names
	for (j = 0; oldLine[i] != '='; j++) {
		k = 0;
		if (oldLine[i] == '\t' || oldLine[i] == ' ' || oldLine[i] == '<' || oldLine[i] == '>') {
			j--, i++;
		}
		else if (oldLine[i] == '+' && oldLine[i - 1] == ' ') {
			i++, j--;
		}
		else if (oldLine[i] == '(' || oldLine[i] == ')') {
			oldLine[i] = ' ';
			i++, j--;
		}
		else {
			while (oldLine[i] != ' ' && oldLine[i] != '\n') {
				if (oldLine[i] == '(' || oldLine[i] == ')') {
					i++;
				}
				else if (oldLine[i] == '<' || oldLine[i] == '>')
					break;
				else {
					reactants[j][k] = oldLine[i];
					i = i + 1, k = k + 1;
				}
			}
			reactants[j][k] = '\0';
			if (reactants[j][0] == 'E' && reactants[j][1] == '\0')
				reactants[j][0] = 'e', reactants[j][1] = 'm', reactants[j][2] = 'i', reactants[j][3] = 'n', reactants[j][4] = 'u', reactants[j][5] = 's', reactants[j][6] = '\0';
			if (reactants[j][k - 1] == '+' && reactants[j][k] == '\0')
				reactants[j][k - 1] = 'p', reactants[j][k] = 'l', reactants[j][k + 1] = 'u', reactants[j][k + 2] = 's', reactants[j][k + 3] = '\0';
		}
	}
	i = 0, j = 0, k = 0;
	char products[10][12];							// create a string array to hold product names
	while (oldLine[i] != '>')
		i++;
	i++;
	while (oldLine[i + 2] != '.' && oldLine[i + 2] != 'e') {
		k = 0;
		if (oldLine[i] == '\t')
			i++;
		else if (oldLine[i] == ' ')
			i++;
		else if (oldLine[i] == '+' && oldLine[i - 1] == ' ')
			i++;
		else if (oldLine[i] == '(' || oldLine[i] == ')') {
			oldLine[i] = ' ';
			i++;
		}
		else {
			while (oldLine[i] != ' ' && oldLine[i] != '\t') {
				if (oldLine[i] == '(' || oldLine[i] == ')')
					i++;
				else {
					products[j][k] = oldLine[i];
					i++, k++;
				}
			}
			products[j][k] = '\0';
			if (products[j][0] == 'E' && products[j][1] == '\0')
				products[j][0] = 'e', products[j][1] = 'm', products[j][2] = 'i', products[j][3] = 'n', products[j][4] = 'u', products[j][5] = 's', products[j][6] = '\0';
			if (products[j][k - 1] == '+' && products[j][k] == '\0')
				products[j][k - 1] = 'p', products[j][k] = 'l', products[j][k + 1] = 'u', products[j][k + 2] = 's', products[j][k + 3] = '\0';
			j++;
		}
	}
	for (j = 0; j < numR; j++) {									// printing reaction names
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
	for (j = 0; j < numP; j++) {									// printing product names
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


void print_begin(char oldLine[], char prevLine[], char lastLine[], char newLine[], int ind, int e, int Q, int way, int tot, int run, int numR, int numP, FILE* output) {

	if (run == 0) {																		// Begin printing first function
		fprintf(output, "  if (IONIZATIONREACTION[%d])\n", tot);						// begin printing information in output file
		if (way == 2)
			fprintf(output, "\t  add_to_W_fwbw_%dr%dp", numR, numP);
		else
			fprintf(output, "\t  add_to_W_fw_%dr%dp", numR, numP);
		if (ind == 1)
			fprintf(output, "_fit4(");									// printing fit information
		else
			fprintf(output, "(");
		print_elem0(oldLine, output, numR, numP);
		print_num0(oldLine, prevLine, numR, numP, ind, output);
		if (e == 1)														// printing temperature Te vs T information
			fprintf(output, "Te, X, W);");
		else
			fprintf(output, "T, X, W);");
		fprintf(output, "\n\n");
	}
	else if (run == 1) {
		fprintf(output, "  if (IONIZATIONREACTION[%d])\n", tot);						// begin printing information in output file
		if (way == 2)
			fprintf(output, "\t  add_to_dW_fwbw_%dr%dp", numR, numP);
		else
			fprintf(output, "\t  add_to_dW_fw_%dr%dp", numR, numP);
		if (ind == 1)
			fprintf(output, "_fit4(");									// printing fit information
		else
			fprintf(output, "(");
		print_elem0(oldLine, output, numR, numP);
		print_num0(oldLine, prevLine, numR, numP, ind, output);
		if (e == 1)														// printing temperature Te vs T information
			fprintf(output, "Te, X, dWdTe, dWdrhok);");
		else
			fprintf(output, "T, X, dWdT, dWdrhok);");
		fprintf(output, "\n\n");
	}
	else if (run == 2 && Q != 0) {
		fprintf(output, "  \t\tif (IONIZATIONREACTION[%d])\n", tot);
		fprintf(output, "  \t\t\tadd_to_Qei(spec");
		print_Q2(oldLine, prevLine, numR, numP, ind, output, newLine, lastLine, Q, e, run);					// func HERE__________printing R/P/#/ending of run 2

	}
	else if (run == 3 && Q == 1) {
		fprintf(output, "  \t\tif (IONIZATIONREACTION[%d])\n", tot);
		fprintf(output, "  \t\t\tadd_to_dQei(spec");
		print_Q2(oldLine, prevLine, numR, numP, ind, output, newLine, lastLine, Q, e, run);			// func HERE__________printing R/P/#/ending of run 3

	}
}


void check_above(char oldLine[], char prevLine[], char lastLine[], char newLine[], int* ind, int* e, int* Q) {

	*ind = check_fit(prevLine);
	if (*ind == 1) {
		*e = check_elec(lastLine);								// call function to check Te or T in reaction
		if (check_Q(newLine) == 1)
			*Q = 1;
	}
	else {
		*e = check_elec(prevLine);
		if (check_Q(prevLine) == 1)
			*Q = 2;
		if (*Q == 0)
			if (check_Q(lastLine) == 1)
				*Q = 3;
	}
}

void make_function(int run, FILE* input, FILE* output) {
	char newLine[300], oldLine[300], prevLine[300], lastLine[300];
	int loc = 0, numR = 0, numP = 0, tot = 0, ind = 0, check = 0, e = 0, way = 10, Q = 0;

	rewind(input);
	fgets(oldLine, 300, input);
	fgets(prevLine, 300, input);
	fgets(lastLine, 300, input);

	while (fgets(newLine, 300, input) != NULL) {					  // searches each line until the end of the document
		loc++;
		Q = 0, way = 10, e = 0, check = 0, ind = 0;
		if (loc == 41 && run == 3)
			loc = loc;
		check = search_current_line(oldLine, &way, &tot, output, loc);
		if (tot == 18)
			tot = tot;
		if (check == 1 || check == 3) {
			ind = check_fit(prevLine);
			check_above(oldLine, prevLine, lastLine, newLine, &ind, &e, &Q);
			numR = number_of_reactants(oldLine);
			numP = number_of_products(oldLine);
			print_begin(oldLine, prevLine, lastLine, newLine, ind, e, Q, way, tot, run, numR, numP, output);
		}
		strcpy(oldLine, prevLine);									// shifts line storage to next line
		strcpy(prevLine, lastLine);
		strcpy(lastLine, newLine);

	}
}


int search_current_line(char Line_str[], int* way, int* tot, FILE* output, int loc) {						// This functions checks a line for a reaction
	int i, num = 0;
	for (i = 0; Line_str[i] != '\n'; i++) {
		if (Line_str[i] == '=' && Line_str[i + 1] == '>' && (Line_str[i - 2] != ' ' || Line_str[i + 2] != ' ') && Line_str[i - 1] == '<') {
			*tot = *tot + 1;
			fprintf(output, "\n  \t\t\t// ************ERROR: Unable to read reaction #%d on line %d.\n\n", *tot, loc);
			break;
		}
		if (Line_str[i] == '=' && Line_str[i + 1] == '>' && ((Line_str[i - 1] != ' ' && Line_str[i - 1] != '<') || Line_str[i + 2] != ' ')) {
			*tot = *tot + 1;
			fprintf(output, "\n  \t\t\t// ************ERROR: Unable to read reaction #%d on line %d.\n\n", *tot, loc);
			break;
		}
		if (Line_str[i] == '=' && Line_str[i + 1] == '>' && Line_str[i - 1] == '<') {
			*tot = *tot + 1;
			num = 1;
			*way = 2;
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
		fprintf(output, "\n  \t\t\t// ************ERROR: Reaction #%d ignored on line %d due to commenting out.\n\n", *tot, loc);
	}

	return num;																 // returns 0 if there is no reaction and 1 if there is a reaction
}
