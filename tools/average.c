#include "share.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

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

void detectRowCols ( FILE * filepointer, int *numrows, int *numcols );
void CleanFile ( char *fileName );

int main ( int argc, char **argv ) {
  char output[255], inputfileprefix[255];
  char fileName[255], ichar[255], firstfileName[255], newfileName[255];
  int start, end, increm, expected_numfiles, num_files, i, j, numrows, numcols, numrows_infile,
    numcols_infile, i_element, j_element;
  bool validOptions = TRUE;
  int scanned;
  FILE *filepointer;
  expected_numfiles = 0;
  num_files = 0;
  numrows=numcols=numrows_infile=numcols_infile=0;
  if ( chkarg ( argc, argv, "-inputfileprefix" ) == 0 )
    validOptions = FALSE;       //check if argument flags exist
  if ( chkarg ( argc, argv, "-output" ) == 0 )
    validOptions = FALSE;
  if ( chkarg ( argc, argv, "-start" ) == 0 )
    validOptions = FALSE;
  if ( chkarg ( argc, argv, "-end" ) == 0 )
    validOptions = FALSE;

  if ( chkarg ( argc, argv, "-start" ) != 0 ) {
    if ( sscanf ( argv[chkarg ( argc, argv, "-start" ) + 1], "%d", &start ) != 1 )
      validOptions = FALSE;
  }
  if ( chkarg ( argc, argv, "-end" ) != 0 ) {
    if ( sscanf ( argv[chkarg ( argc, argv, "-end" ) + 1], "%d", &end ) != 1 )
      validOptions = FALSE;
  }
  if ( chkarg ( argc, argv, "-increm" ) != 0 ) {
    if ( sscanf ( argv[chkarg ( argc, argv, "-increm" ) + 1], "%d", &increm ) != 1 )
      validOptions = FALSE;
  }

  if ( !validOptions ) {
    fprintf ( stderr, "\nFlags:\n\n"
              "Flag               \tArg                     \tArg Type     \tRequired?\n"
              "-----------------------------------------------------------------------------------\n"
              "-inputfileprefix   \tinput file prefix       \tstring       \tY\n"
              "-output            \toutput text file        \tstring       \tY\n"
              "-start             \tcounter start           \tinteger      \tY\n"
              "-end               \tcounter end             \tinteger      \tY\n"
              "-increm            \tcounter increments      \tinteger      \tN\n"
              "                   \t(i=start, start+increm, ..)\n"
              "Eg: \n"
              "average -inputfileprefix Cf_x. -start 10 -end 30 -output avg_out.txt\nwill average the values in all columns in the specified files (Cf_x.10,Cf_x.11...Cf_x.30) and output the average to the file avg_out.txt .\n" );
    exit ( 1 );
  }

  sscanf ( argv[chkarg ( argc, argv, "-inputfileprefix" ) + 1], "%s", inputfileprefix );        //get the values corresponding to flags
  sscanf ( argv[chkarg ( argc, argv, "-output" ) + 1], "%s", output );
  sscanf ( argv[chkarg ( argc, argv, "-start" ) + 1], "%d", &start );
  sscanf ( argv[chkarg ( argc, argv, "-end" ) + 1], "%d", &end );
  if ( start == end ) {
    printf ( "\nTrying to average using 1 file. Please input more than 1 file.\n" );
    exit ( 1 );
  }

  int suffix_max_length, prefix_max_length;
  char suffix[256];
  sprintf ( suffix, "%d", end );
  suffix_max_length = strlen ( suffix );
  prefix_max_length = strlen ( inputfileprefix );
  if ( suffix_max_length + prefix_max_length > 250 ) {
    printf
      ( "Filenames are too long. Please ensure that the length of the filenames does not exceed 250 characters. Aborting.\n" );
    exit ( 1 );
  }

  if ( strlen ( output ) > 254 ) {
    printf ( "Output file name is too long. Please choose a shorter name. Aborting\n" );
    exit ( 1 );
  }

  if ( chkarg ( argc, argv, "-increm" ) != 0 )
    sscanf ( argv[chkarg ( argc, argv, "-increm" ) + 1], "%d", &increm );
  else
    increm = 1;
  //printf("inputfileprefix=%s\noutput=%s\nstart=%d\nend=%d\nincrem=%d\n",inputfileprefix,output,start,end,increm);
  strcpy ( fileName, inputfileprefix ); //check how many files exist and ensure files with ____ prefix do not exist already
  for ( i = start; i <= end; i = i + increm ) {
    expected_numfiles = expected_numfiles + 1;
    sprintf ( ichar, "%d", i );
    strcat ( fileName, ichar );
    strcpy ( newfileName, "____" );
    strcat ( newfileName, fileName );
    if ( access ( newfileName, R_OK ) != -1 ) {
      printf
        ( "Warning: Ensure files such as %s with four underscores preceding the name of the files to be averaged over do not exist. This program will create files with ____ as prefix to the existing file names and delete them. Please remove files such as %s ... before proceeding. Aborting.\n\n",
          newfileName, newfileName );
      exit ( 1 );
    }
    if ( access ( fileName, R_OK ) != -1 ) {
      num_files = num_files + 1;
    } else
      printf ( "File %s not found. Skipping.\n", fileName );
    strcpy ( fileName, inputfileprefix );
  }
  printf ( "Found %d/%d files.\n", num_files, expected_numfiles );

  strcpy ( firstfileName, inputfileprefix );    //check if first file exists. If it does, find the number of rows and cols from it.
  sprintf ( ichar, "%d", start );
  strcat ( firstfileName, ichar );
  if ( !( access ( firstfileName, R_OK ) != -1 ) ) {
    printf
      ( "\n\nError. Make sure that atleast the first file (%s) exists. The number of rows and columns are detected from the first file.\n\n",
        fileName );
    exit ( 1 );
  } else {
    CleanFile ( firstfileName );
    strcpy ( fileName, "____" );
    strcat ( fileName, firstfileName );
    filepointer = fopen ( fileName, "r" );
    detectRowCols ( filepointer, &numrows, &numcols );
    //printf("\nNumber of rows=%d,Number of columns=%d\n",numrows,numcols);
    if ( filepointer != NULL )
      fclose ( filepointer );
  }

  for ( i = start + increm; i <= end; i = i + increm ) {        //check if number of rows and cols are same in all files which exist.
    sprintf ( ichar, "%d", i );
    strcpy ( fileName, inputfileprefix );
    strcat ( fileName, ichar );
    if ( access ( fileName, R_OK ) != -1 ) {
      CleanFile ( fileName );
      strcpy ( newfileName, "____" );
      strcat ( newfileName, fileName );
      filepointer = fopen ( newfileName, "r" );
      detectRowCols ( filepointer, &numrows_infile, &numcols_infile );
      if ( filepointer != NULL )
        fclose ( filepointer );
      if ( ( numrows_infile != numrows ) || ( numcols_infile != numcols ) ) {
        printf
          ( "\n\nDifferent number of rows and columns in file %s compared to the first file %s. Aborting.\n\n",
            fileName, firstfileName );
        printf ( "numrowsinfile=%d,numcolsinfile=%d,numrows=%d,numcols=%d\n", numrows_infile, numcols_infile,
                 numrows, numcols );
        exit ( 1 );
      } else {
        numrows_infile = 0;
        numcols_infile = 0;
      }
    }
  }

  double avg_matrix[numrows][numcols], current_matrix[numrows][numcols];

  for ( i = 0; i < numrows; i++ ) {     //initialize elements of matrix  that will store the elements to 0.
    for ( j = 0; j < numcols; j++ ) {
      avg_matrix[i][j] = 0.0;
    }
  }

  for ( i = start; i <= end; i = i + increm ) { //sum up elements in the files
    strcpy ( fileName, inputfileprefix );
    sprintf ( ichar, "%d", i );
    strcat ( fileName, ichar );
    if ( access ( fileName, R_OK ) != -1 ) {
      CleanFile ( fileName );
      strcpy ( newfileName, "____" );
      strcat ( newfileName, fileName );
      filepointer = fopen ( newfileName, "r" );
      for ( i_element = 0; i_element < numrows; i_element++ ) {
        for ( j_element = 0; j_element < numcols; j_element++ ) {
          scanned = fscanf ( filepointer, "%lf%*c", &current_matrix[i_element][j_element] );
          if ( isnan ( current_matrix[i_element][j_element] ) == 1
               || isinf ( current_matrix[i_element][j_element] ) == 1 ) {
            printf
              ( "Number read is Nan or Inf. Please check if element in row %d and column %d of file %s is not NaN or Inf. Aborting\n\n",
                i_element + 1, j_element + 1, fileName );
            for ( i = start; i <= end; i = i + increm ) {       //remove generated files before aborting
              sprintf ( ichar, "%d", i );
              strcpy ( fileName, "____" );
              strcat ( fileName, inputfileprefix );
              strcat ( fileName, ichar );
              remove ( fileName );
            }
            if ( filepointer != NULL )
              fclose ( filepointer );
            exit ( 1 );
          }
          if ( scanned != 1 ) {
            printf ( "\n\nElement in row %d and column %d of file %s is not a number. Aborting.\n\n",
                     i_element + 1, j_element + 1, fileName );
            for ( i = start; i <= end; i = i + increm ) {       //remove generated files before aborting
              sprintf ( ichar, "%d", i );
              strcpy ( fileName, "____" );
              strcat ( fileName, inputfileprefix );
              strcat ( fileName, ichar );
              remove ( fileName );
            }
            if ( filepointer != NULL )
              fclose ( filepointer );
            exit ( 1 );
          }
          avg_matrix[i_element][j_element] =
            avg_matrix[i_element][j_element] + current_matrix[i_element][j_element];
        }
      }
      if ( filepointer != NULL )
        fclose ( filepointer );
    }
  }

  for ( i = start; i <= end; i = i + increm ) { //remove generated files after having stored them in the matrix.
    sprintf ( ichar, "%d", i );
    strcpy ( fileName, "____" );
    strcat ( fileName, inputfileprefix );
    strcat ( fileName, ichar );
    remove ( fileName );
  }

  for ( i = 0; i < numrows; i++ ) {     //average
    for ( j = 0; j < numcols; j++ ) {
      avg_matrix[i][j] = avg_matrix[i][j] / num_files;
    }
  }

  filepointer = fopen ( output, "w" );
  if ( access ( output, W_OK ) != -1 ) {        //print to file
    for ( i = 0; i < numrows; i++ ) {
      for ( j = 0; j < numcols; j++ ) {
        fprintf ( filepointer, "%.12E", avg_matrix[i][j] );
        if ( j != ( numcols - 1 ) )
          fprintf ( filepointer, "\t" );
      }
      fprintf ( filepointer, "\n" );
    }
  }
  if ( filepointer != NULL )
    fclose ( filepointer );
  printf ( "Averaged over %d files. Results outputted to %s. Done.\n", num_files, output );
  return ( 0 );
}

void CleanFile ( char *fileName ) {
  char newFileName[255] = "____";
  char ch;
  int lastchardelim = 0;
  int new_line = 1;
  int first_char_encountered = 0;
  FILE *filepointer;
  FILE *tempfilepointer;
  filepointer = fopen ( fileName, "r" );
  if ( access ( fileName, R_OK ) != -1 ) {
    strcat ( newFileName, fileName );
    tempfilepointer = fopen ( newFileName, "w" );
    if ( access ( newFileName, W_OK ) != -1 ) {
      while ( ( ch = fgetc ( filepointer ) ) != EOF ) {
        if ( ch == ' ' || ch == '\t' || ch == ',' || ch == '\n' ) {
          if ( ch == '\n' ) {
            if ( new_line != 1 ) {
              fputc ( '\n', tempfilepointer );
              new_line = 1;
              lastchardelim = 1;
              first_char_encountered = 0;
            } else if ( new_line == 1 )
              lastchardelim = 1;
          } else {
            lastchardelim = 1;
            if ( first_char_encountered == 1 )
              new_line = 0;
            else
              new_line = 1;
          }
        } else {
          if ( lastchardelim == 1 ) {
            if ( new_line == 1 ) {
              fputc ( ch, tempfilepointer );
              lastchardelim = 0;
              new_line = 0;
              first_char_encountered = 1;
            } else {
              fputc ( ' ', tempfilepointer );
              fputc ( ch, tempfilepointer );
              lastchardelim = 0;
              new_line = 0;
              first_char_encountered = 1;
            }
          } else {
            fputc ( ch, tempfilepointer );
            lastchardelim = 0;
            new_line = 0;
            first_char_encountered = 1;
          }
        }
      }
    } else {
      printf ( "Aborting. Could not write file %s. Check write permission in directory.\n\n", newFileName );
      if ( tempfilepointer != NULL )
        fclose ( tempfilepointer );
      if ( filepointer != NULL )
        fclose ( filepointer );
      exit ( 1 );
    }
  } else {
    printf ( "Aborting. Could not read file %s\n\n", fileName );
    if ( filepointer != NULL )
      fclose ( filepointer );
    exit ( 1 );
  }
  if ( filepointer != NULL )
    fclose ( filepointer );
  if ( tempfilepointer != NULL )
    fclose ( tempfilepointer );
}

void detectRowCols ( FILE * filepointer, int *numrows, int *numcols ) {
  char ch;
  int oneColFile = 0;
  int numspaces_perline = 0, elements_perline = 1;

  while ( ( ch = fgetc ( filepointer ) ) != EOF ) {
    if ( ch == ' ' ) {
      elements_perline++;
    }
    if ( ch == '\n' ) {
      break;
    }
  }

  if ( elements_perline == 1 )
    oneColFile = 1;
  fseek ( filepointer, 0, 0 );
  if ( oneColFile == 1 ) {
    while ( ( ch = fgetc ( filepointer ) ) != EOF ) {
      if ( ch == '\n' ) {
        *numrows = *numrows + 1;
      }
      if ( ch == ' ' ) {
        printf ( "\n\nDifferent number of elements per line. Aborting.\n\n" );
        exit ( 1 );
      }
    }
    *numcols = 1;
    return;
  }

  while ( ( ch = fgetc ( filepointer ) ) != EOF ) {
    if ( ch == ' ' ) {
      numspaces_perline++;
    }
    if ( ch == '\n' ) {
      if ( ( numspaces_perline + 1 ) != elements_perline ) {
        printf
          ( "\n\nEncountered different number of elements in a row compared to the first row. Aborting.\n\n" );
      } else {
        numspaces_perline = 0;
        *numrows = *numrows + 1;
      }
    }
  }
  *numcols = elements_perline;
  return;
}
