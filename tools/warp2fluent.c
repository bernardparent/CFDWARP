// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2025 Felipe Martin Rodriguez Fuentes

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

/* 
This tool interpolates solution variables from a CFDWARP interpolation file to the face centers of a FLUENT mesh.

inputs: CFDWARP interpolation data file and a FLUENT mesh file containing only face center coordinates.
outputs: A .dat file of interpolated variables to each FLUENT mesh face; and code UDF.c to load the interpolated profiles in the FLUENT GUI

The FLUENT region should be contained or overlap in the CFDWARP domain. Normally, this CFWWARP domain will be an outflow boundary,
which would be the FLUENT inflow to continue a simulation there.

This impementation does not feed back CFDWARP data outside a boundary to the FLUENT inflow.

To run this tool at each time step for a transient CFDWARP simulation, use WriteInterpolationFile() within the control file Cycle() with a system() call.
*/

#include "share.h"
#include <exm.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define EOS 0

/*
This function writes the UDF code that is interpreted inside the FLUENT GUI to load the interpolated profiles
*/
void write_udf_file(const char *udf_filename, char **varnames, long numvars, long numdim, const char *interpolated_filename) {
  FILE *udf = fopen(udf_filename, "w");
  if (!udf) {
    fprintf(stderr, "Could not open %s for writing\n", udf_filename);
    return;
  }

  fprintf(udf, "#include \"udf.h\"\n");
  fprintf(udf, "#include <stdio.h>\n");
  fprintf(udf, "#include <stdlib.h>\n");
  fprintf(udf, "#include <math.h>\n\n");

  fprintf(udf, "#define MAX_NODES 200000\n");
  fprintf(udf, "#define TOLERANCE 1.0e-6\n\n");
  fprintf(udf, "#define FILENAME \"C:\\\\Users\\\\mrffe\\\\Downloads\\\\warp_output.dat\"\n\n");

  fprintf(udf, "int data_read = 0;\n");
  fprintf(udf, "int num_nodes = 0;\n");
  if (numdim == 2) {
    fprintf(udf, "float x_data[MAX_NODES], y_data[MAX_NODES];\n");
  } else {
    fprintf(udf, "float x_data[MAX_NODES], y_data[MAX_NODES], z_data[MAX_NODES];\n");
  }
  for (long i = 0; i < numvars; i++) {
    fprintf(udf, "float %s_data[MAX_NODES];\n", varnames[i]);
  }
  fprintf(udf, "\n");
  // read the data interpolated to the FLUENT mesh
  fprintf(udf,
    "void read_warp_data(void)\n"
    "{\n"
    "    FILE *fp;\n"
    "    int temp_nodes = 0;\n"
    "    fp = fopen(FILENAME, \"r\");\n"
    "    if (fp == NULL) {\n"
    "        printf(\"Error: Could not open file.\");\n"
    "        return;\n"
    "    }\n\n"
  );

  fprintf(udf, "    while (1) {\n");
  if (numdim == 2) {
    fprintf(udf,
        "        if (fscanf(fp, \"%%e %%e\", &x_data[temp_nodes], &y_data[temp_nodes]) != 2) break;\n");
  } else {
    fprintf(udf,
        "        if (fscanf(fp, \"%%e %%e %%e\", &x_data[temp_nodes], &y_data[temp_nodes], &z_data[temp_nodes]) != 3) break;\n");
  }

  for (long i = 0; i < numvars; i++) {
    fprintf(udf,
      "        if (fscanf(fp, \"%%e\", &%s_data[temp_nodes]) != 1) break;\n",
      varnames[i]);
    }

  fprintf(udf,
      "        temp_nodes++;\n"
      "        if (temp_nodes >= MAX_NODES) {\n"
      "            printf(\"Warning: MAX_NODES exceeded, truncating input.\\\\n\");\n"
      "            break;\n"
      "        }\n"
      "    }\n"
      "    fclose(fp);\n"
      "    num_nodes = temp_nodes;\n"
      "    data_read = 1;\n"
      "    printf(\"Successfully read %%d nodes from %s\\\\n\", num_nodes);\n"
      "}\n\n",
      interpolated_filename
  );

  // write profiles for each interpolated variable
  for (long i = 0; i < numvars; i++) {
    fprintf(udf,
      "DEFINE_PROFILE(%s_profile, thread, position)\n"
      "{\n"
      "    face_t f;\n"
      "    real xyz[ND_ND];\n"
      "    float x_f, y_f%s;\n"
      "    float dx, dy%s, dist, val;\n"
      "    int i, matched;\n"
      "    if (!data_read) read_warp_data();\n"
      "    begin_f_loop(f, thread)\n"
      "    {\n"
      "        F_CENTROID(xyz, f, thread);\n"
      "        x_f = xyz[0];\n"
      "        y_f = (ND_ND > 1) ? xyz[1] : 0.0;\n",
            varnames[i],
            (numdim == 3 ? ", z_f" : ""),
            (numdim == 3 ? ", dz" : "")
        );

        if (numdim == 3) {
            fprintf(udf, "        z_f = (ND_ND > 2) ? xyz[2] : 0.0;\n");
        }

        fprintf(udf,
            "        matched = 0;\n"
            "        for (i = 0; i < num_nodes; i++) {\n"
            "            dx = x_f - x_data[i];\n"
            "            dy = y_f - y_data[i];\n");
        if (numdim == 3) {
            fprintf(udf, "            dz = z_f - z_data[i];\n");
            fprintf(udf, "            dist = sqrt(dx*dx + dy*dy + dz*dz);\n");
        } else {
            fprintf(udf, "            dist = sqrt(dx*dx + dy*dy);\n");
        }
        fprintf(udf,
            "            if (dist < TOLERANCE) {\n"
            "                val = %s_data[i];\n"
            "                F_PROFILE(f, thread, position) = val;\n"
            "                matched = 1;\n"
            "                break;\n"
            "            }\n"
            "        }\n"
            "        if (!matched) F_PROFILE(f, thread, position) = 0.0;\n"
            "    }\n"
            "    end_f_loop(f, thread)\n"
            "}\n\n",
            varnames[i]
        );
    }

    fclose(udf);
}


void find_interpolation_weight(double *x_mesh, double *x_file, double *dx1_file, double *dx2_file, double *dx3_file,
                               int nd, double radiusmax2, double *thisweight)
{
  double distance;
  EXM_mat_t mat1,mat2,mat3,mat1inv;
  long dim;
  distance=0.0;
  for (dim=0; dim<nd; dim++) distance+=sqr(x_mesh[dim]-x_file[dim]);
  *thisweight=0.0;
  if (distance<radiusmax2) {
    EXM_init_matrix(&mat1, nd, nd);
    for (dim=0; dim<nd; dim++){
      mat1.cont[EXM_aim(mat1.glm,dim,0)]=dx1_file[dim];
      mat1.cont[EXM_aim(mat1.glm,dim,1)]=dx2_file[dim];
      if (nd>2) mat1.cont[EXM_aim(mat1.glm,dim,2)]=dx3_file[dim];
    }
    EXM_init_matrix(&mat1inv, nd, nd);
  
    EXM_invert_matrix_analytical(mat1, &mat1inv);
    EXM_init_matrix(&mat2, nd, 1);
    for (dim=0; dim<nd; dim++){
      mat2.cont[EXM_aim(mat2.glm,dim,0)]=x_mesh[dim]-x_file[dim];
    }
    EXM_init_matrix(&mat3, nd, 1);
    EXM_multiply_matrices(mat1inv, mat2, &mat3);
    *thisweight=0.0;
    for (dim=0; dim<nd; dim++) *thisweight+=fabs(pow(fabs(mat3.cont[EXM_aim(mat3.glm,dim,0)]),3.0));
    *thisweight=fabs(pow(*thisweight,1.0/3.0));
    *thisweight=max(1e-16,max(0.0003-(*thisweight)*0.00001,1.0-(*thisweight)));
    EXM_free_matrix(&mat1);
    EXM_free_matrix(&mat1inv);
    EXM_free_matrix(&mat2);
    EXM_free_matrix(&mat3);
  }
}


int main(int argc, char **argv) {
  FILE *datafile, *meshdata, *fluent_output;
  bool VALIDOPTIONS = TRUE, IS_3D = FALSE, ISREAD_EMFIELD = TRUE;
  char *options, *input, *mesh, *output;
  char **initvar_fluid_names_file = NULL;
  char **initvar_emfield_names_file = NULL;
  int RET;
  options = NULL;
  bool FORMAT001;
  long cnt,cnterror,dim,l_file,l_mesh;
  long numflux_read,numspec_read,numdim_read,numnodes,numfaces,windowis,windowie,iter,numvars_fluid_file;
  long numflux_emfield_read,numdim_emfield_read,numnodes_emfield,numvars_emfield_file;
  double effiter_U,effiter_R,effiter_U_emfield,effiter_R_emfield;
  double tmp_dt,tmp_time,CFL,Lc;
  double *radiusmax2_file,*radiusmax2_emfield_file,thisweight,*weight,thisweight_emfield,*weight_emfield;
  char data_format_str[100], initvar_fluid_str_file[500], initvar_emfield_str_file[500];
  numfaces=0;
  
  input =(char *)malloc(400*sizeof(char));
  mesh =(char *)malloc(400*sizeof(char));
  output =(char *)malloc(400*sizeof(char));

  if (process_flag_string(argc, argv, "-ii", &input)!=2) VALIDOPTIONS=FALSE; 
  if (process_flag_string(argc, argv, "-m", &mesh)!=2) VALIDOPTIONS=FALSE; 
  if (process_flag_string(argc, argv, "-o", &output)!=2) VALIDOPTIONS=FALSE;    

  if (!VALIDOPTIONS) {
    fprintf (stderr, "\nPROBLEM WITH FLAGS\n\n\nRequired and Optional Flags:\n\n"
             "Flag      Arg                     Arg Type     Required?\n"
             "----------------------------------------------------------\n"
             "-ii       CFDWARP datafile        string       YES\n"
             "-m        FLUENT mesh             string       YES\n"
             "-o        Output file             string       YES\n"
             "----------------------------------------------------------\n"
             "Eg: \n"
             "./warp2fluent -ii datafile -m fluentmesh.dat -o warp_output.dat\nWill read input CFDWARP interpolation datafile (2D or 3D) obtained with the -oi flag, the FLUENT mesh and assign the CFDWARP fluid and/or emfield variables to the FUENT face centers.\n\n");
    exit (EXIT_FAILURE);
  }
  RET = find_remaining_options ( argc, argv, &options );
  if ( RET >= 1 ) {
    fprintf ( stderr, "\n\nThe following command line options could not be processed: %s\n\n", options );
    exit (EXIT_FAILURE);
  }
  else {
    datafile = fopen(input, "r");
    if (datafile==NULL) {
      fprintf ( stderr, "\n\nThe following input file could not be found: %s\n\n", input );
      exit (EXIT_FAILURE);
    }
    meshdata = fopen(mesh, "r");
    if (meshdata==NULL) {
      fprintf ( stderr, "\n\nThe following input file could not be found: %s\n\n", mesh );
      exit (EXIT_FAILURE);
    }
  }

  // Do the fluid properties first. This would include Tv and Te last if plasma

  for (cnt=0; cnt<16; cnt++){
    if (fscanf(datafile,"%c",&(data_format_str[cnt]))!=1) {
      fprintf ( stderr,"Problem with fscanf when reading CFDWARP datafile.");
    }
  }

  data_format_str[16]=EOS;
  fprintf ( stdout, "Reading CFDWARP data file...");
  FORMAT001=FALSE;  
  if (strcmp("WARPINTFORMAT001",data_format_str)==0) {
    FORMAT001=TRUE;
  }
  fprintf ( stdout, "done.\n");

  if (FORMAT001) {
    if (fscanf(datafile," numnodes=%ld nf=%ld nd=%ld ns=%ld windowis=%ld windowie=%ld iter=%ld effiter_U=%lg effiter_R=%lg CFL=%lg time=%lg dt=%lg vars_fluid=\"%[^\"]\" vars_emfield=\"%[^\"]\"%*[^\n]",
        &numnodes,&numflux_read,&numdim_read,&numspec_read,&windowis,&windowie,
        &iter,&effiter_U,&effiter_R,&CFL,&tmp_time,&tmp_dt,initvar_fluid_str_file,initvar_emfield_str_file)!=14) fprintf(stderr,"Problem with fscanf in CFDWARP input file preamble.");
  
    fgetc(datafile);
  } else {
    fprintf ( stderr,"Interpolation file format unknown.");
  }
  find_words_from_string(initvar_fluid_str_file, " ", &initvar_fluid_names_file, &numvars_fluid_file);
  find_words_from_string(initvar_emfield_str_file, " ", &initvar_emfield_names_file, &numvars_emfield_file);
  for (cnt = 0; cnt < numvars_emfield_file; cnt++){
    if (strcmp(initvar_emfield_names_file[cnt],"NONE")==0){
      numvars_emfield_file=0;
      ISREAD_EMFIELD=FALSE;
    }
  }

  double (*initvar_fluid_file)[numvars_fluid_file] = malloc(numnodes*(numvars_fluid_file)*sizeof(double));
  if (numdim_read>2) IS_3D = TRUE;

  // Rename the velocity components here from "V[i]" to make it easier to parse in the FLUENT UDF
  strcpy(initvar_fluid_names_file[0], "V_X");
  strcpy(initvar_fluid_names_file[1], "V_Y");
  if (IS_3D) strcpy(initvar_fluid_names_file[2], "V_Z");

  typedef double dim_t[numdim_read];
  dim_t *dx1_file,*dx2_file,*dx3_file,*x_file,*x_mesh;

  x_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
  dx1_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
  dx2_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
  dx3_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
  radiusmax2_file=(double *)malloc(numnodes*sizeof(double));
  
  for (l_file=0; l_file<numnodes; l_file++){
    cnterror=0;
    if (fread(initvar_fluid_file[l_file], (numvars_fluid_file)*sizeof(double), 1, datafile)!=1) cnterror++;
    if (fread(x_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (fread(dx1_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (fread(dx2_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (IS_3D) {
        if (fread(dx3_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    } else {
      for (dim=0; dim<numdim_read; dim++) dx3_file[l_file][dim]=0.0;
    }
    if (cnterror>0) fprintf ( stderr,"Could not read all data properly.");
    radiusmax2_file[l_file]=0.0e0;
    for (dim=0; dim<numdim_read; dim++){
      radiusmax2_file[l_file]+=sqr(fabs(dx1_file[l_file][dim])
                                  +fabs(dx2_file[l_file][dim])
                                  +fabs(dx3_file[l_file][dim]));
    }
    radiusmax2_file[l_file]*=1.1;
  }

  // Next, do the EMFIELD properties if any (phi)

  // Note that the EMFIELD part might have a different number of nodes if some fluid regions have been cut,
  // but fluid and EMFIELD properties (phi) will be interpolated to a single FLUENT bdry mesh

  if (ISREAD_EMFIELD){
    for (cnt=0; cnt<16; cnt++){
      if (fscanf(datafile,"%c",&(data_format_str[cnt]))!=1){
        fprintf( stderr,"Problem with fscanf in emfield part of interpolation file.");
      }
    }
    data_format_str[16]=EOS;
    FORMAT001=FALSE;  
    if (strcmp("WARPINTFORMAT001",data_format_str)==0) {
      FORMAT001=TRUE;
    }

    if (FORMAT001) {
      if (fscanf(datafile," numnodes_emfield=%ld nfe=%ld nd=%ld Lc=%lg effiter_U_emfield=%lg effiter_R_emfield=%lg%*[^\n]",
              &numnodes_emfield,&numflux_emfield_read,&numdim_emfield_read,&Lc,&effiter_U_emfield,&effiter_R_emfield)!=6){
              fprintf ( stderr,"Problem reading emfield preambule in interpolation file.");
      }
      fgetc(datafile);
    } else {
      fprintf( stderr,"Interpolation file format unknown for EMfield variables.");
    }
  }

  double (*initvar_emfield_file)[numvars_emfield_file] = malloc(numnodes_emfield*numvars_emfield_file*sizeof(double));

  dim_t *dx1_emfield_file,*dx2_emfield_file,*dx3_emfield_file,*x_emfield_file;

  x_emfield_file=(dim_t *)malloc(numnodes_emfield*sizeof(dim_t));
  dx1_emfield_file=(dim_t *)malloc(numnodes_emfield*sizeof(dim_t));
  dx2_emfield_file=(dim_t *)malloc(numnodes_emfield*sizeof(dim_t));
  dx3_emfield_file=(dim_t *)malloc(numnodes_emfield*sizeof(dim_t));
  radiusmax2_emfield_file=(double *)malloc(numnodes_emfield*sizeof(double));

  if (ISREAD_EMFIELD) {
    for (l_file=0; l_file<numnodes_emfield; l_file++){
      cnterror=0;
      if (fread(initvar_emfield_file[l_file], numvars_emfield_file*sizeof(double), 1, datafile)!=1) cnterror++;
    if (fread(x_emfield_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (fread(dx1_emfield_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (fread(dx2_emfield_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (IS_3D) {
        if (fread(dx3_emfield_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    } else {
      for (dim=0; dim<numdim_emfield_read; dim++) dx3_emfield_file[l_file][dim]=0.0;
    }
    if (cnterror>0) fprintf ( stderr,"Could not read all data properly.");
    radiusmax2_emfield_file[l_file]=0.0e0;
    for (dim=0; dim<numdim_emfield_read; dim++){
      radiusmax2_emfield_file[l_file]+=sqr(fabs(dx1_emfield_file[l_file][dim])
                                  +fabs(dx2_emfield_file[l_file][dim])
                                  +fabs(dx3_emfield_file[l_file][dim]));
    }
    radiusmax2_emfield_file[l_file]*=1.1;
    }
  }

  // Print some output to screen to verify

  printf("Dimensions = %ld\n",numdim_read);
  printf("Number of fluid variables = %ld\n",numvars_fluid_file);
  printf("Number of emfield variables = %ld\n",numvars_emfield_file);
  printf("Fluid variables = ");
  for (cnt = 0; cnt < numvars_fluid_file; cnt++) printf("%s ",initvar_fluid_names_file[cnt]);
  printf("\n");
  printf("EMField variables = ");
  for (cnt = 0; cnt < numvars_emfield_file; cnt++) printf("%s ",initvar_emfield_names_file[cnt]);
  printf("\n");

  // Parse FLUENT inputs

  fluent_output = fopen(output, "w");
  if (!fluent_output) {
    fprintf(stderr, "Could not open FLUENT output file %s for writing\n", output);
    return 1;
  }

  fprintf ( stdout, "Reading FLUENT boundary mesh face data...");
  meshdata = fopen(mesh, "r");
  if (!meshdata) {
    fprintf(stderr, "Could not open FLUENT mesh file %s\n", mesh);
    return 1;
  }

  double tmp[numdim_read];
  if (IS_3D) {
    while (fscanf(meshdata, "%lf %lf %lf", &tmp[0], &tmp[1], &tmp[2]) == 3) {
      numfaces++;
    }
  } else {
    while (fscanf(meshdata, "%lf %lf", &tmp[0], &tmp[1]) == 2) {
      numfaces++;
    }
  }
  printf("%ld FLUENT face centers read...",numfaces);

  x_mesh=(dim_t *)malloc(numfaces*sizeof(dim_t));
  rewind(meshdata);
  for (l_mesh = 0; l_mesh < numfaces; l_mesh++) {
    if (IS_3D) {
      if (fscanf(meshdata, "%lf %lf %lf", &x_mesh[l_mesh][0], &x_mesh[l_mesh][1], &x_mesh[l_mesh][2]) != 3) 
      return 1;
    } else {
      if (fscanf(meshdata, "%lf %lf", &x_mesh[l_mesh][0], &x_mesh[l_mesh][1]) != 2) 
      return 1;      
    }
  }
  fprintf ( stdout, "done.\n");

  // Interpolate to FLUENT mesh

  // FLUID
  double (*interpolated_vars)[numvars_fluid_file] = malloc(numfaces * sizeof(*interpolated_vars));
  weight=(double *)malloc(sizeof(double)*numfaces);
  for (l_mesh = 0; l_mesh < numfaces; l_mesh++) weight[l_mesh]=0.0;
  fprintf ( stdout, "Interpolating values from CFDWARP to FLUENT face centers...");
  for (l_file = 0; l_file < numnodes; l_file++) {  
    for (l_mesh = 0; l_mesh < numfaces; l_mesh++) {
      find_interpolation_weight(x_mesh[l_mesh], x_file[l_file], dx1_file[l_file], dx2_file[l_file], dx3_file[l_file], numdim_read, radiusmax2_file[l_file], &thisweight);
      if (thisweight>1e-99) {
        weight[l_mesh]+=thisweight;
        for (cnt=0; cnt<numvars_fluid_file; cnt++) {
          interpolated_vars[l_mesh][cnt]+=thisweight*initvar_fluid_file[l_file][cnt];
        }
      }
    } 
  }

  // EMFIELD 
  double (*interpolated_emfield_vars)[numvars_emfield_file] = malloc(numfaces * sizeof(*interpolated_emfield_vars));
  weight_emfield=(double *)malloc(sizeof(double)*numfaces);
  for (l_mesh = 0; l_mesh < numfaces; l_mesh++) weight_emfield[l_mesh]=0.0;
  if (ISREAD_EMFIELD) {
    for (l_file = 0; l_file < numnodes_emfield; l_file++) {  
      for (l_mesh = 0; l_mesh < numfaces; l_mesh++) {
        find_interpolation_weight(x_mesh[l_mesh], x_emfield_file[l_file], dx1_emfield_file[l_file], dx2_emfield_file[l_file], dx3_emfield_file[l_file], numdim_emfield_read, radiusmax2_emfield_file[l_file], &thisweight_emfield);
        if (thisweight_emfield>1e-99) {
          weight_emfield[l_mesh]+=thisweight_emfield;
          for (cnt=0; cnt<numvars_emfield_file; cnt++) {
            interpolated_emfield_vars[l_mesh][cnt]+=thisweight_emfield*initvar_emfield_file[l_file][cnt];
          }
        }
      } 
    }
  }

  // Write to output
 for (l_mesh = 0; l_mesh < numfaces; l_mesh++) {
    if (IS_3D) {
      fprintf(fluent_output,"%E %E %E ",x_mesh[l_mesh][0],x_mesh[l_mesh][1],x_mesh[l_mesh][2]);
    } else {
      fprintf(fluent_output,"%E %E ",x_mesh[l_mesh][0],x_mesh[l_mesh][1]);
    }
    if (weight[l_mesh]>1e-99){
      for (cnt=0; cnt<numvars_fluid_file; cnt++) {
      interpolated_vars[l_mesh][cnt]=interpolated_vars[l_mesh][cnt]/weight[l_mesh];
        fprintf(fluent_output,"%E ",interpolated_vars[l_mesh][cnt]);
      }
    } 
    if (weight_emfield[l_mesh]>1e-99){
      for (cnt=0; cnt<numvars_emfield_file; cnt++) {
      interpolated_emfield_vars[l_mesh][cnt]=interpolated_emfield_vars[l_mesh][cnt]/weight_emfield[l_mesh];
        fprintf(fluent_output,"%E ",interpolated_emfield_vars[l_mesh][cnt]);
      }
    } 
    fprintf(fluent_output,"\n");
  }
  fprintf ( stdout, "done.\n");

  // Write Fluent UDF
  fprintf ( stdout, "Writing FLUENT UDF...");
  write_udf_file("UDF.c",initvar_fluid_names_file,numvars_fluid_file,numdim_read,output);
  fprintf ( stdout, "done.\n");

  fprintf ( stdout, "Find interpolated data file and UDF.c in the directory where the warp2fluent executable is located.\n");

  fclose(datafile);
  fclose(meshdata);
  fclose(fluent_output);

  free(input);
  free(mesh);
  free(output);
	free(options);
  free(initvar_fluid_file);
  free(initvar_emfield_file);
  free(interpolated_vars);
  free(interpolated_emfield_vars);
  free(x_file);
  free(x_emfield_file);
  free(x_mesh);
  free(dx1_file);
  free(dx1_emfield_file);
  free(dx2_file);
  free(dx2_emfield_file);
  free(dx3_file);
  free(dx3_emfield_file);
  free(weight);
  free(weight_emfield);
  for (int i = 0; i < numvars_fluid_file; i++) free(initvar_fluid_names_file[i]);
  free(initvar_fluid_names_file);
  for (int i = 0; i < numvars_emfield_file; i++) free(initvar_emfield_names_file[i]);
  free(initvar_emfield_names_file);

  return(EXIT_SUCCESS);
}



