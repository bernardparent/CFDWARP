// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1999-2016 Bernard Parent

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

#include <src/data.h>
#include <model/_model.h>
#include <cycle/_cycle.h>
#include <unistd.h>

#define dt_steady 1.0e99

#define DATATYPE_BINARY 1
#define DATATYPE_ASCII 2

#define INITVARTYPE_FLUID 1
#define INITVARTYPE_EMFIELD 2

#ifdef _3DL
  #define SUBZONE_DESIRED_WIDTH 40
#else
  #define SUBZONE_DESIRED_WIDTH 40
#endif

#define MIN_NUMSUBZONE_PER_THREAD 3

typedef struct {
  long newindex;
  double newvalue;
} mapvar_t;  

long _ai_mpidatafile(gl_t *gl, long i, long j, long k) {
  long ii;
  ii=(i-gl->domain_all.is);
  if2DL(
    ii=ii*(gl->domain_all.je-gl->domain_all.js+1)+(j-gl->domain_all.js);
    if3DL(ii=ii*(gl->domain_all.ke-gl->domain_all.ks+1)+(k-gl->domain_all.ks);)
  )
  return(ii);
}


void find_NODEVALID_on_domain_all(np_t *np, gl_t *gl, int TYPELEVEL, bool *NODEVALID){
  long i,j,k;
#ifdef DISTMPI
  int rank,thisrank;
  int THISNODEVALID;
#endif
  
  for_ijk(gl->domain_lim_all,is,js,ks,ie,je,ke){
        NODEVALID[_ai_all(gl,i,j,k)]=FALSE;
  }

#ifdef DISTMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        if (j==gl->domain_all.js && k==gl->domain_all.ks) MPI_Barrier(MPI_COMM_WORLD);
        thisrank=_node_rank(gl, i, j, k);
        if (thisrank==rank) THISNODEVALID=(int)(is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL));
        MPI_Bcast(&THISNODEVALID,1,MPI_INT,thisrank,MPI_COMM_WORLD);           
        assert(THISNODEVALID==TRUE || THISNODEVALID==FALSE);          
        NODEVALID[_ai_all(gl,i,j,k)]=(bool)THISNODEVALID;
  }
  MPI_Barrier(MPI_COMM_WORLD);
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        NODEVALID[_ai_all(gl,i,j,k)]=is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL);
  }
#endif

	
}


void read_data_file_binary_ascii(char *filename, np_t *np, gl_t *gl, long level, int DATATYPE){
  FILE *datafile;
  char data_format_str[100];
  long i,j,k,flux,cnt,tmp1,tmp2,tmp3;
  double CFLmem;
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
  flux_t Res;
  bool NORES=FALSE;
  long NOREScount=0;
#endif
#ifndef UNSTEADY
  double tmp_double;
  long tmp_long;
#endif
  bool FORMAT010;
  bool *NODEVALID;
  flux_t U;
#ifdef EMFIELD
  double Lcmem;
  fluxemfield_t Uemfield;
#endif
#ifdef DISTMPI
  int rank,numproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  NODEVALID=(bool *)malloc((gl->domain_lim_all.ie-gl->domain_lim_all.is+1)
#ifdef _2DL
                          *(gl->domain_lim_all.je-gl->domain_lim_all.js+1)
#endif
#ifdef _3DL
                          *(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1)
#endif
                          *sizeof(bool));


  CFLmem=gl->CFL;
#ifdef EMFIELD
  Lcmem=gl->Lc;
#endif

  datafile = fopen(filename, "r");
  if (datafile==NULL)
    fatal_error("Having problems opening datafile %s.",filename);
  for (cnt=0; cnt<16; cnt++) {
    if (fscanf(datafile,"%c",&(data_format_str[cnt]))!=1) fatal_error("Problem with fscanf in read_data_file_binary().");
  }
  data_format_str[16]=EOS;
  wfprintf(stdout,"Reading data file %s ",filename);
  if (level!=0) wfprintf(stdout,"to time level minus %ld ",level);
  FORMAT010=FALSE;
  switch (DATATYPE){
    case DATATYPE_BINARY:
      if (strcmp("WARPBINFORMAT010",data_format_str)==0) {
        wfprintf(stdout,"in CFDWARP binary format 010..");
        FORMAT010=TRUE;
      }
    break;
    case DATATYPE_ASCII:
      if (strcmp("WARPASCFORMAT010",data_format_str)==0) {
        wfprintf(stdout,"in CFDWARP ASCII format 010..");
        FORMAT010=TRUE;
      }
    break;
  }
  if (FORMAT010) {
    if (level==0) {
      if (fscanf(datafile," windowis=%ld windowie=%ld iter=%ld effiter_U=%lg effiter_R=%lg CFL=%lg",
             &(gl->window.is),&(gl->window.ie),
             &(gl->iter),&(gl->effiter_U),&(gl->effiter_R),&(gl->CFL))!=6) fatal_error("Problem with fscanf in read_data_file_binary().");
      if (fscanf(datafile," nd=%ld ns=%ld nf=%ld",
             &tmp1,&tmp2,
             &tmp3)!=3) fatal_error("Problem with fscanf in read_data_file_binary().");
      if (tmp1!=nd) fatal_error("Data file has %ld dimensions but CFDWARP is compiled with %ld dimensions.",tmp1,nd);
      if (tmp2!=ns) fatal_error("Data file has %ld species but CFDWARP is compiled with %ld species.",tmp2,ns);
      if (tmp3!=nf) fatal_error("Data file has %ld fluxes but CFDWARP is compiled with %ld fluxes.",tmp3,nf);
      if (fscanf(datafile," is=%ld ie=%ld",
             &tmp1,&tmp2)!=2) fatal_error("Problem with fscanf in read_data_file_binary().");
      if ((tmp2-tmp1)!=(gl->domain_all.ie-gl->domain_all.is)) fatal_error("Data file has %ld grid lines along i but the control file specifies %ld grid lines.",tmp2-tmp1,(gl->domain_all.ie-gl->domain_all.is));
#ifdef _2DL
      if (fscanf(datafile," js=%ld je=%ld",
             &tmp1,&tmp2)!=2) fatal_error("Problem with fscanf in read_data_file_binary().");
      if ((tmp2-tmp1)!=(gl->domain_all.je-gl->domain_all.js)) fatal_error("Data file has %ld grid lines along j but the control file specifies %ld grid lines.",tmp2-tmp1,(gl->domain_all.je-gl->domain_all.js));
#endif
#ifdef _3DL
      if (fscanf(datafile," ks=%ld ke=%ld",
             &tmp1,&tmp2)!=2) fatal_error("Problem with fscanf in read_data_file_binary().");
      if ((tmp2-tmp1)!=(gl->domain_all.ke-gl->domain_all.js)) fatal_error("Data file has %ld grid lines along k but the control file specifies %ld grid lines.",tmp2-tmp1,(gl->domain_all.ke-gl->domain_all.ks));
#endif

#if defined(UNSTEADY)
      if (fscanf(datafile," time=%lg",&(gl->time))!=1) fatal_error("Problem reading time variable within fscanf in read_data_file_binary().");
#else
      if (fscanf(datafile," time=%lg",&tmp_double)!=1) fatal_error("Problem reading time variable within fscanf in read_data_file_binary().");
#endif


#ifdef UNSTEADY
      if (fscanf(datafile," dt=%lg",&(gl->dt))!=1) fatal_error("Problem reading dt variable within fscanf in read_data_file_binary().");
#else
      if (fscanf(datafile," dt=%lg",&tmp_double)!=1) fatal_error("Problem reading dt variable within fscanf in read_data_file_binary().");
#endif

#ifdef UNSTEADY
      if (fscanf(datafile," timelevel=%ld",&(gl->timelevel))!=1) fatal_error("Problem reading timelevel variable within fscanf in read_data_file_binary().");
#else
      if (fscanf(datafile," timelevel=%ld",&tmp_long)!=1) fatal_error("Problem reading timelevel variable within fscanf in read_data_file_binary().");
#endif
      
#ifdef EMFIELD
      if (fscanf(datafile," Lc=%lg effiter_U_emfield=%lg effiter_R_emfield=%lg",&(gl->Lc),&(gl->effiter_U_emfield),&(gl->effiter_R_emfield))!=3) fatal_error("Problem reading EMFIELD variables within fscanf in read_data_file_binary().");
#endif


#if defined(UNSTEADY) && defined(_AVERAGEDRATES)
      if (fscanf(datafile," averagedrates_time=%lg AVERAGEDRATES=%ld",&(gl->averagedrates_time),&(gl->AVERAGEDRATES))!=2) fatal_error("Problem reading AVERAGEDRATES variables within fscanf in read_data_file_binary().");
#endif

      if (fscanf(datafile,"%*[^\n]")!=0) fatal_error("Problem with fscanf in read_data_file_binary().");


    } else {
      if (fscanf(datafile," %*[^\n]")!=0) fatal_error("Problem with fscanf in read_data_file_binary().");
    }
    fgetc(datafile);
  }

  if (!FORMAT010){
    fatal_error("Data file format invalid.");
  }

  wfprintf(stdout,"fluid.");
  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_FLUID, NODEVALID);
  wfprintf(stdout,".");

  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){

#ifdef DISTMPI
        if (rank==0) {
#endif
          if (NODEVALID[_ai_all(gl,i,j,k)]) {
            switch (DATATYPE){
              case DATATYPE_BINARY:
                if (fread(U, sizeof(flux_t), 1, datafile)!=1)
                  fatal_error("Could not read all data properly.");
              break;
              case DATATYPE_ASCII:
                for (flux=0; flux<nf; flux++){
                  if (fscanf(datafile,"%lg%*[^\n]",&(U[flux]))!=1)
                    fatal_error("Could not read all data properly.");
                }
              break;
              default:
                fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
            }
          }
#ifdef DISTMPI
        }
        MPI_Bcast_Node(&U, nf, MPI_DOUBLE, 0, MPI_COMM_WORLD, i, j, k, gl);
        if (j==gl->domain_all.js && k==gl->domain_all.ks) MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (is_node_in_zone(i,j,k,gl->domain_lim)) {
          for (flux=0; flux<nf; flux++){
            if (level==0) np[_ai(gl,i,j,k)].bs->U[flux]=U[flux];
#ifdef UNSTEADY
            if (level==1) np[_ai(gl,i,j,k)].bs->Um1[flux]=U[flux];
            #if _RESTIME_BW > 2 
              if (level==2) np[_ai(gl,i,j,k)].bs->Um2[flux]=U[flux];
            #endif
            #if _RESTIME_BW > 3 
              if (level==3) np[_ai(gl,i,j,k)].bs->Um3[flux]=U[flux];
            #endif
#endif
          }
          np[_ai(gl,i,j,k)].INIT_FLUID=TRUE;
        }
  }


#ifdef EMFIELD
  wfprintf(stdout,"emfield.");
  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_EMFIELD, NODEVALID);
  wfprintf(stdout,".");
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
#ifdef DISTMPI
        if (rank==0) {
#endif
          if (NODEVALID[_ai_all(gl,i,j,k)]) {
            switch (DATATYPE){
              case DATATYPE_BINARY:
                if (fread(Uemfield, sizeof(fluxemfield_t), 1, datafile)!=1)
                  fatal_error("Could not read all data properly.");
              break;
              case DATATYPE_ASCII:
                for (flux=0; flux<nfe; flux++){
                  if (fscanf(datafile,"%lg%*[^\n]",&(Uemfield[flux]))!=1)
                    fatal_error("Could not read all data properly.");
                }
              break;
              default:
                fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
            }
          }
#ifdef DISTMPI
        }
        MPI_Bcast_Node(&Uemfield, nfe, MPI_DOUBLE, 0, MPI_COMM_WORLD, i, j, k, gl);
        if (j==gl->domain_all.js && k==gl->domain_all.ks) MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (is_node_in_zone(i,j,k,gl->domain_lim)) {
          for (flux=0; flux<nfe; flux++) {
            if (level==0) np[_ai(gl,i,j,k)].bs->Uemfield[flux]=Uemfield[flux];
#ifdef UNSTEADY
            if (level==1) np[_ai(gl,i,j,k)].bs->Uemfieldm1[flux]=Uemfield[flux];
#endif
          }
          np[_ai(gl,i,j,k)].INIT_EMFIELD=TRUE;
        }
  }
#endif



#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
  wfprintf(stdout,"trap.");
  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_FLUID, NODEVALID);
  wfprintf(stdout,".");
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){

#ifdef DISTMPI
        if (rank==0) {
#endif
          if (NODEVALID[_ai_all(gl,i,j,k)]) {
            switch (DATATYPE){
              case DATATYPE_BINARY:
                if (fread(Res, sizeof(flux_t), 1, datafile)!=1){
                  NORES=TRUE;
                  NOREScount++;
                }
              break;
              case DATATYPE_ASCII:
                for (flux=0; flux<nf; flux++){
                  if (fscanf(datafile,"%lg%*[^\n]",&(Res[flux]))!=1){
                    NORES=TRUE;
                    NOREScount++;
                  }  
                }
              break;
              default:
                fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
            }

          }
#ifdef DISTMPI
        }
		MPI_Bcast(&NORES, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        if (!NORES) MPI_Bcast_Node(&Res, nf, MPI_DOUBLE, 0, MPI_COMM_WORLD, i, j, k, gl);
        if (j==gl->domain_all.js && k==gl->domain_all.ks) MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (is_node_in_zone(i,j,k,gl->domain_lim)) {
          for (flux=0; flux<nf; flux++){
            if (!NORES) np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux]=Res[flux];
            else if (NORES) np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux]=0.0;
          }
        }
		NORES=FALSE;
  }
  if(NOREScount>0) wfprintf(stdout,"WARNING: The residual at the previous time step could not be found within the data file %s. The residual has been set to zero..",filename);
#endif


#if defined(UNSTEADY) && defined(_AVERAGEDRATES)
  double averagedrates[numaveragedrates];
  bool NOAVERAGEDRATES=FALSE;
  long NOAVERAGEDRATEScount=0;
  
  wfprintf(stdout,"averagedrates.");
  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_FLUID, NODEVALID);
  wfprintf(stdout,".");
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){

#ifdef DISTMPI
        if (rank==0) {
#endif
          if (NODEVALID[_ai_all(gl,i,j,k)]) {
            switch (DATATYPE){
              case DATATYPE_BINARY:
                if (fread(averagedrates, numaveragedrates*sizeof(double), 1, datafile)!=1){
                  NOAVERAGEDRATES=TRUE;
                  NOAVERAGEDRATEScount++;
                }
              break;
              case DATATYPE_ASCII:
                for (flux=0; flux<numaveragedrates; flux++){
                  if (fscanf(datafile,"%lg%*[^\n]",&(averagedrates[flux]))!=1){
                    NOAVERAGEDRATES=TRUE;
                    NOAVERAGEDRATEScount++;
                  }  
                }
              break;
              default:
                fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
            }

          }
#ifdef DISTMPI
        }
		MPI_Bcast(&NOAVERAGEDRATES, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        if (!NOAVERAGEDRATES) MPI_Bcast_Node(&averagedrates, numaveragedrates, MPI_DOUBLE, 0, MPI_COMM_WORLD, i, j, k, gl);
        if (j==gl->domain_all.js && k==gl->domain_all.ks) MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (is_node_in_zone(i,j,k,gl->domain_lim)) {
          for (flux=0; flux<numaveragedrates; flux++){
            if (!NOAVERAGEDRATES) np[_ai(gl,i,j,k)].bs->averagedrates[flux]=averagedrates[flux];
            else if (NOAVERAGEDRATES) np[_ai(gl,i,j,k)].bs->averagedrates[flux]=0.0;
          }
        }
		NOAVERAGEDRATES=FALSE;
  }
  if(NOAVERAGEDRATEScount>0) wfprintf(stdout,"WARNING: The averaged rates could not be found within the data file %s. The averaged rates have been set to zero..",filename);
#endif


  fclose(datafile);

  wfprintf(stdout,"done;\n");
  if (level!=0) gl->CFL=CFLmem;
#ifdef EMFIELD
  if (level!=0) gl->Lc=Lcmem;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank!=0) {
    gl->effiter_U=0.0;
    gl->effiter_R=0.0;
    #ifdef EMFIELD
    gl->effiter_U_emfield=0.0;
    gl->effiter_R_emfield=0.0;   
    #endif
  }
#endif

  free(NODEVALID);
}


void read_data_file_mpi(char *filename, np_t *np, gl_t *gl, long level){
  FILE *datafile;
  char data_format_str[100];
  long i,j,k,flux,cnt,tmp1,tmp2,tmp3;
  double CFLmem;
#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
  flux_t Res;
  bool NORES=FALSE;
  long NOREScount=0;
#endif
#ifndef UNSTEADY
  double tmp_double;
#endif
  bool FORMAT010;
  flux_t U;
#ifdef EMFIELD
  double Lcmem;
  fluxemfield_t Uemfield;
#endif
#ifdef DISTMPI
  int rank,numproc,fd;
  long offset=0;
  zone_t domain;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CFLmem=gl->CFL;
#ifdef EMFIELD
  Lcmem=gl->Lc;
#endif

  datafile = fopen(filename, "r");
  if (datafile==NULL)
    fatal_error("Having problems opening datafile %s.",filename);
  for (cnt=0; cnt<16; cnt++) {
    if (fscanf(datafile,"%c",&(data_format_str[cnt]))!=1) fatal_error("Problem with fscanf in read_data_file_mpi().");
  }
  data_format_str[16]=EOS;
  wfprintf(stdout,"Reading data file %s ",filename);
  if (level!=0) wfprintf(stdout,"to time level minus %ld ",level);
  FORMAT010=FALSE;
  if (strcmp("WARPMPBFORMAT010",data_format_str)==0) {
    wfprintf(stdout,"in CFDWARP binary (MPI) format 010..");
    FORMAT010=TRUE;
  }

  if (FORMAT010) {
    if (level==0) {
      if (fscanf(datafile," windowis=%ld windowie=%ld iter=%ld effiter_U=%lg effiter_R=%lg CFL=%lg",
             &(gl->window.is),&(gl->window.ie),
             &(gl->iter),&(gl->effiter_U),&(gl->effiter_R),&(gl->CFL))!=6) fatal_error("Problem with fscanf in read_data_file_mpi().");
      if (fscanf(datafile," nd=%ld ns=%ld nf=%ld",
             &tmp1,&tmp2,
             &tmp3)!=3) fatal_error("Problem with fscanf in read_data_file_mpi().");
      if (tmp1!=nd) fatal_error("Data file has %ld dimensions but CFDWARP is compiled with %ld dimensions.",tmp1,nd);
      if (tmp2!=ns) fatal_error("Data file has %ld species but CFDWARP is compiled with %ld species.",tmp2,ns);
      if (tmp3!=nf) fatal_error("Data file has %ld fluxes but CFDWARP is compiled with %ld fluxes.",tmp3,nf);
      if (fscanf(datafile," is=%ld ie=%ld",
             &tmp1,&tmp2)!=2) fatal_error("Problem with fscanf in read_data_file_mpi().");
      if ((tmp2-tmp1)!=(gl->domain_all.ie-gl->domain_all.is)) fatal_error("Data file has %ld grid lines along i but the control file specifies %ld grid lines.",tmp2-tmp1,(gl->domain_all.ie-gl->domain_all.is));
#ifdef _2DL
      if (fscanf(datafile," js=%ld je=%ld",
             &tmp1,&tmp2)!=2) fatal_error("Problem with fscanf in read_data_file_mpi().");
      if ((tmp2-tmp1)!=(gl->domain_all.je-gl->domain_all.js)) fatal_error("Data file has %ld grid lines along j but the control file specifies %ld grid lines.",tmp2-tmp1,(gl->domain_all.je-gl->domain_all.js));
#endif
#ifdef _3DL
      if (fscanf(datafile," ks=%ld ke=%ld",
             &tmp1,&tmp2)!=2) fatal_error("Problem with fscanf in read_data_file_mpi().");
      if ((tmp2-tmp1)!=(gl->domain_all.ke-gl->domain_all.js)) fatal_error("Data file has %ld grid lines along k but the control file specifies %ld grid lines.",tmp2-tmp1,(gl->domain_all.ke-gl->domain_all.ks));
#endif

#if defined(UNSTEADY)
      if (fscanf(datafile," time=%lg",&(gl->time))!=1) fatal_error("Problem reading time variable within fscanf in read_data_file_mpi().");
#else
      if (fscanf(datafile," time=%lg",&tmp_double)!=1) fatal_error("Problem reading time variable within fscanf in read_data_file_mpi().");
#endif


#ifdef UNSTEADY
      if (fscanf(datafile," dt=%lg",&(gl->dt))!=1) fatal_error("Problem reading dt variable within fscanf in read_data_file_mpi().");
#else
      if (fscanf(datafile," dt=%lg",&tmp_double)!=1) fatal_error("Problem reading dt variable within fscanf in read_data_file_mpi().");
#endif


#ifdef EMFIELD
      if (fscanf(datafile," Lc=%lg effiter_U_emfield=%lg effiter_R_emfield=%lg",&(gl->Lc),&(gl->effiter_U_emfield),&(gl->effiter_R_emfield))!=3) fatal_error("Problem reading EMFIELD variables within fscanf in read_data_file_mpi().");
#endif


      if (fscanf(datafile,"%*[^\n]")!=0) fatal_error("Problem with fscanf in read_data_file_mpi().");


    } else {
      if (fscanf(datafile," %*[^\n]")!=0) fatal_error("Problem with fscanf in read_data_file_mpi().");
    }
    fgetc(datafile);
  }

  if (!FORMAT010){
    fatal_error("Data file format invalid.");
  }
#ifdef DISTMPI
  fd=fileno(datafile);
  if(level==0) offset=ftell(datafile);
  domain=_zone_intersection(gl->domain_all,gl->domain_lim);
#endif

  wfprintf(stdout,"fluid.");
  wfprintf(stdout,".");

#ifdef DISTMPI
  for_ijk(domain,is,js,ks,ie,je,ke){
    if(pread(fd,U,sizeof(flux_t),_ai_mpidatafile(gl,i,j,k)*sizeof(flux_t)+offset)==-1)
      fatal_error("Could not read all data properly.");
    for (flux=0; flux<nf; flux++){
      if (level==0) np[_ai(gl,i,j,k)].bs->U[flux]=U[flux];
#ifdef UNSTEADY
      if (level==1) np[_ai(gl,i,j,k)].bs->Um1[flux]=U[flux];
      #if _RESTIME_BW > 2
        if (level==2) np[_ai(gl,i,j,k)].bs->Um2[flux]=U[flux];
      #endif
      #if _RESTIME_BW > 3
        if (level==3) np[_ai(gl,i,j,k)].bs->Um3[flux]=U[flux];
      #endif
#endif
    }
    np[_ai(gl,i,j,k)].INIT_FLUID=TRUE;
  }
  offset+=((gl->domain_all.ie-gl->domain_all.is+1)*(gl->domain_all.je-gl->domain_all.js+1)*(gl->domain_all.ke-gl->domain_all.ks+1)*sizeof(flux_t));
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    if (fread(U, sizeof(flux_t), 1, datafile)!=1)
      fatal_error("Could not read all data properly.");
    for (flux=0; flux<nf; flux++){
      if (level==0) np[_ai(gl,i,j,k)].bs->U[flux]=U[flux];
#ifdef UNSTEADY
      if (level==1) np[_ai(gl,i,j,k)].bs->Um1[flux]=U[flux];
      #if _RESTIME_BW > 2
        if (level==2) np[_ai(gl,i,j,k)].bs->Um2[flux]=U[flux];
      #endif
      #if _RESTIME_BW > 3
        if (level==3) np[_ai(gl,i,j,k)].bs->Um3[flux]=U[flux];
      #endif
#endif
    }
    np[_ai(gl,i,j,k)].INIT_FLUID=TRUE;
  }
#endif


#ifdef EMFIELD
  wfprintf(stdout,"emfield.");
  wfprintf(stdout,".");
#ifdef DISTMPI
  for_ijk(domain,is,js,ks,ie,je,ke){
    if(pread(fd,Uemfield,sizeof(fluxemfield_t),_ai_mpidatafile(gl,i,j,k)*sizeof(fluxemfield_t)+offset)==-1)
      fatal_error("Could not read all data properly.");
    for (flux=0; flux<nfe; flux++) {
      if (level==0) np[_ai(gl,i,j,k)].bs->Uemfield[flux]=Uemfield[flux];
#ifdef UNSTEADY
      if (level==1) np[_ai(gl,i,j,k)].bs->Uemfieldm1[flux]=Uemfield[flux];
#endif
    }
    np[_ai(gl,i,j,k)].INIT_EMFIELD=TRUE;
  }
  offset+=((gl->domain_all.ie-gl->domain_all.is+1)*(gl->domain_all.je-gl->domain_all.js+1)*(gl->domain_all.ke-gl->domain_all.ks+1)*sizeof(fluxemfield_t));
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    if (fread(Uemfield, sizeof(fluxemfield_t), 1, datafile)!=1)
      fatal_error("Could not read all data properly.");
    for (flux=0; flux<nfe; flux++) {
      if (level==0) np[_ai(gl,i,j,k)].bs->Uemfield[flux]=Uemfield[flux];
#ifdef UNSTEADY
      if (level==1) np[_ai(gl,i,j,k)].bs->Uemfieldm1[flux]=Uemfield[flux];
#endif
    }
    np[_ai(gl,i,j,k)].INIT_EMFIELD=TRUE;
  }
#endif
#endif


#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
  wfprintf(stdout,"trap.");
  wfprintf(stdout,".");
#ifdef DISTMPI
  for_ijk(domain,is,js,ks,ie,je,ke){
    if(pread(fd,Res,sizeof(flux_t),_ai_mpidatafile(gl,i,j,k)*sizeof(flux_t)+offset)==-1){
      printf("failed to read value at %ld,%ld,%ld;",i,j,k);
      NORES=TRUE;
      NOREScount++;
    }
    for (flux=0; flux<nf; flux++) {
      if (!NORES) np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux]=Res[flux];
      else if (NORES) np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux]=0.0;
    }
    NORES=FALSE;
  }
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    if (fread(Res, sizeof(flux_t), 1, datafile)!=1){
      NORES=TRUE;
      NOREScount++;
    }
    for (flux=0; flux<nf; flux++) {
      if (!NORES) np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux]=Res[flux];
      else if (NORES) np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux]=0.0;
    }
    NORES=FALSE;
  }
#endif
  if(NOREScount>0) wfprintf(stdout,"WARNING: The residual at the previous time step could not be found within the data file %s. The residual has been set to zero..",filename);
#endif

  fclose(datafile);

  wfprintf(stdout,"done;\n");
  if (level!=0) gl->CFL=CFLmem;
#ifdef EMFIELD
  if (level!=0) gl->Lc=Lcmem;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank!=0) {
    gl->effiter_U=0.0;
    gl->effiter_R=0.0;
    #ifdef EMFIELD
    gl->effiter_U_emfield=0.0;
    gl->effiter_R_emfield=0.0;
    #endif
  }
#endif

}


void write_data_file_binary_ascii(char *filename, np_t *np, gl_t *gl, int DATATYPE){
  FILE *datafile;
  long i,j,k,sizeoftmpflux;
  bool *NODEVALID;
#ifdef EMFIELD
  double effiter_U_emfield,effiter_R_emfield;
#endif
  long flux;
  double effiter_U,effiter_R;
#ifdef DISTMPI
  zone_t domain;
  flux_t U;
  int rank,proc,numproc;
  MPI_Status MPI_Status1;
#endif
  long ntmpflux=nf;
#ifdef EMFIELD
  ntmpflux=max(ntmpflux,nfe);
#endif
#if defined(UNSTEADY) && defined(_AVERAGEDRATES)
  ntmpflux=max(ntmpflux,numaveragedrates);
  double averagedrates[numaveragedrates];
#endif  
  sizeoftmpflux=ntmpflux*sizeof(double);
  typedef double tmpflux_t[ntmpflux];
  tmpflux_t *fluxtmp;


  fluxtmp=(tmpflux_t *)malloc((gl->domain_lim_all.ie-gl->domain_lim_all.is+1)
#ifdef _2DL
                          *(gl->domain_lim_all.je-gl->domain_lim_all.js+1)
#endif
#ifdef _3DL
                          *(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1)
#endif
                          *sizeoftmpflux);
  NODEVALID=(bool *)malloc((gl->domain_lim_all.ie-gl->domain_lim_all.is+1)
#ifdef _2DL
                          *(gl->domain_lim_all.je-gl->domain_lim_all.js+1)
#endif
#ifdef _3DL
                          *(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1)
#endif
                          *sizeof(bool));



  effiter_U=gl->effiter_U;
  effiter_R=gl->effiter_R;
#ifdef EMFIELD
  effiter_U_emfield=gl->effiter_U_emfield;
  effiter_R_emfield=gl->effiter_R_emfield;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Allreduce(&gl->effiter_U, &effiter_U, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&gl->effiter_R, &effiter_R, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef EMFIELD
  MPI_Allreduce(&gl->effiter_U_emfield, &effiter_U_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&gl->effiter_R_emfield, &effiter_R_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif


  datafile = wfopen(filename, "w");
  wfprintf(stdout,"Writing to CFDWARP ");
  switch (DATATYPE){
    case DATATYPE_BINARY:
      wfprintf(stdout,"binary");
      wfprintf(datafile,"WARPBINFORMAT010");
    break;
    case DATATYPE_ASCII:
      wfprintf(stdout,"ASCII");
      wfprintf(datafile,"WARPASCFORMAT010");
    break;
    default:
      fatal_error("DATATYPE must be either DATATYPE_BINARY or DATATYPE_ASCII.");
  }
  wfprintf(stdout," data file %s..",filename);
  wfprintf(datafile," windowis=%ld windowie=%ld iter=%ld effiter_U=%E effiter_R=%E CFL=%E",
           gl->window.is,gl->window.ie,gl->iter,effiter_U,effiter_R,gl->CFL);
  wfprintf(datafile," nd=%ld ns=%ld nf=%ld",nd,ns,nf);
  wfprintf(datafile," is=%ld ie=%ld",gl->domain_all.is,gl->domain_all.ie);
#ifdef _2DL
  wfprintf(datafile," js=%ld je=%ld",gl->domain_all.js,gl->domain_all.je);
#endif
#ifdef _3DL
  wfprintf(datafile," ks=%ld ke=%ld",gl->domain_all.ks,gl->domain_all.ke);
#endif

#if defined(UNSTEADY)
  wfprintf(datafile," time=%E",gl->time);
#else
  wfprintf(datafile," time=%E",0.0);
#endif

#ifdef UNSTEADY
  wfprintf(datafile," dt=%E",gl->dt);
#else
  wfprintf(datafile," dt=%E",dt_steady);
#endif

#if defined(UNSTEADY)
  wfprintf(datafile," timelevel=%ld",gl->timelevel);
#else
  wfprintf(datafile," timelevel=%ld",0);
#endif

#ifdef EMFIELD
  wfprintf(datafile," Lc=%E effiter_U_emfield=%E effiter_R_emfield=%E",gl->Lc,effiter_U_emfield,effiter_R_emfield);
#endif

#if defined(_AVERAGEDRATES)
  #if defined(UNSTEADY) 
    wfprintf(datafile," averagedrates_time=%E AVERAGEDRATES=%ld",gl->averagedrates_time,gl->AVERAGEDRATES);
  #else
    wfprintf(datafile," averagedrates_time=0.0E0 AVERAGEDRATES=0");
  #endif
#endif

  wfprintf(datafile,"\n");

  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_FLUID, NODEVALID);
#ifdef DISTMPI
  for (proc=0; proc<numproc; proc++){
    domain=_domain_from_rank(proc,gl);
    for_ijk(domain,is,js,ks,ie,je,ke){
          if (proc==rank){
            for (flux=0; flux<nf; flux++) U[flux]=np[_ai(gl,i,j,k)].bs->U[flux];
            if (proc!=0) {
              MPI_Send(U,nf,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            }
          } 
          
          if (rank==0 && proc!=0) {
            MPI_Recv(U,nf,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&MPI_Status1);
          }
          
          for (flux=0; flux<nf; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=U[flux];
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        for (flux=0; flux<nf; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->U[flux];
  }
#endif

  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        if (NODEVALID[_ai_all(gl,i,j,k)]) {
          switch (DATATYPE){
            case DATATYPE_BINARY:
              wfwrite(fluxtmp[_ai_all(gl,i,j,k)], sizeof(flux_t), 1, datafile);
            break;
            case DATATYPE_ASCII:
              for (flux=0; flux<nf; flux++)
                wfprintf(datafile, "%18.16E\n",fluxtmp[_ai_all(gl,i,j,k)][flux]);
            break;
            default:
              fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
          }
        }
  }


#ifdef EMFIELD
  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_EMFIELD, NODEVALID);
#ifdef DISTMPI
  for (proc=0; proc<numproc; proc++){
    domain=_domain_from_rank(proc,gl);
    for_ijk(domain,is,js,ks,ie,je,ke){
          if (proc==rank){
            for (flux=0; flux<nfe; flux++) U[flux]=np[_ai(gl,i,j,k)].bs->Uemfield[flux];
            if (proc!=0) {
              MPI_Send(U,nfe,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            }
          } 
          
          if (rank==0 && proc!=0) {
            MPI_Recv(U,nfe,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&MPI_Status1);
          }
          
          for (flux=0; flux<nfe; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=U[flux];
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        for (flux=0; flux<nfe; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->Uemfield[flux];
  }
#endif

  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        if (NODEVALID[_ai_all(gl,i,j,k)]) {
          switch (DATATYPE){
            case DATATYPE_BINARY:
              wfwrite(fluxtmp[_ai_all(gl,i,j,k)], sizeof(fluxemfield_t), 1, datafile);
            break;
            case DATATYPE_ASCII:
              for (flux=0; flux<nfe; flux++)
                wfprintf(datafile, "%18.16E\n",fluxtmp[_ai_all(gl,i,j,k)][flux]);
            break;
            default:
              fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
          }
        }
  }
#endif



#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_FLUID, NODEVALID);
#ifdef DISTMPI
  for (proc=0; proc<numproc; proc++){
    domain=_domain_from_rank(proc,gl);
    for_ijk(domain,is,js,ks,ie,je,ke){
          if (proc==rank){
            for (flux=0; flux<nf; flux++) U[flux]=np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux];
            if (proc!=0) {
              MPI_Send(U,nf,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            }
          } 
          
          if (rank==0 && proc!=0) {
            MPI_Recv(U,nf,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&MPI_Status1);
          }
          
          for (flux=0; flux<nf; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=U[flux];
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        for (flux=0; flux<nf; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux];
  }
#endif

  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        if (NODEVALID[_ai_all(gl,i,j,k)]) {
          switch (DATATYPE){
            case DATATYPE_BINARY:
              wfwrite(fluxtmp[_ai_all(gl,i,j,k)], sizeof(flux_t), 1, datafile);
            break;
            case DATATYPE_ASCII:
              for (flux=0; flux<nf; flux++)
                wfprintf(datafile, "%18.16E\n",fluxtmp[_ai_all(gl,i,j,k)][flux]);
            break;
            default:
              fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
          }
        }
  }
#endif







#if defined(UNSTEADY) && defined(_AVERAGEDRATES)
  find_NODEVALID_on_domain_all(np, gl, TYPELEVEL_FLUID, NODEVALID);
#ifdef DISTMPI
  for (proc=0; proc<numproc; proc++){
    domain=_domain_from_rank(proc,gl);
    for_ijk(domain,is,js,ks,ie,je,ke){
          if (proc==rank){
            for (flux=0; flux<numaveragedrates; flux++) averagedrates[flux]=np[_ai(gl,i,j,k)].bs->averagedrates[flux];
            if (proc!=0) {
              MPI_Send(averagedrates,numaveragedrates,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            }
          } 
          
          if (rank==0 && proc!=0) {
            MPI_Recv(averagedrates,numaveragedrates,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,&MPI_Status1);
          }
          
          for (flux=0; flux<numaveragedrates; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=averagedrates[flux];
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        for (flux=0; flux<numaveragedrates; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->averagedrates[flux];
  }
#endif

  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
        if (NODEVALID[_ai_all(gl,i,j,k)]) {
          switch (DATATYPE){
            case DATATYPE_BINARY:
              wfwrite(fluxtmp[_ai_all(gl,i,j,k)], numaveragedrates*sizeof(double), 1, datafile);
            break;
            case DATATYPE_ASCII:
              for (flux=0; flux<numaveragedrates; flux++)
                wfprintf(datafile, "%18.16E\n",fluxtmp[_ai_all(gl,i,j,k)][flux]);
            break;
            default:
              fatal_error("DATATYPE must be either DATATYPE_ASCII or DATATYPE_BINARY.");
          }
        }
  }
#endif










#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  wfclose(datafile);
  wfprintf(stdout,"done.\n");
  free(fluxtmp);
  free(NODEVALID);
}


void write_data_file_mpi(char *filename, np_t *np, gl_t *gl){
  FILE *datafile;
  long i,j,k;
  flux_t *fluxtmp;
#ifdef EMFIELD
  double effiter_U_emfield,effiter_R_emfield;
#endif
  long flux;
  double effiter_U,effiter_R;
#ifdef DISTMPI
  zone_t domain;
  int rank,fd;
#endif
#ifdef DISTMPI
  long offset=0;
#endif
#ifdef EMFIELD
  assert(nf>=nfe);
#endif
  fluxtmp=(flux_t *)malloc((gl->domain_lim_all.ie-gl->domain_lim_all.is+1)
#ifdef _2DL
                          *(gl->domain_lim_all.je-gl->domain_lim_all.js+1)
#endif
#ifdef _3DL
                          *(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1)
#endif
                          *sizeof(flux_t));
 
  effiter_U=gl->effiter_U;
  effiter_R=gl->effiter_R;
#ifdef EMFIELD
  effiter_U_emfield=gl->effiter_U_emfield;
  effiter_R_emfield=gl->effiter_R_emfield;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Allreduce(&gl->effiter_U, &effiter_U, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&gl->effiter_R, &effiter_R, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef EMFIELD
  MPI_Allreduce(&gl->effiter_U_emfield, &effiter_U_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&gl->effiter_R_emfield, &effiter_R_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif


  datafile = wfopen(filename, "w");
  wfprintf(stdout,"Writing to CFDWARP binary (MPI)");
  wfprintf(datafile,"WARPMPBFORMAT010");
  wfprintf(stdout," data file %s..",filename);
  wfprintf(datafile," windowis=%ld windowie=%ld iter=%ld effiter_U=%E effiter_R=%E CFL=%E",
           gl->window.is,gl->window.ie,gl->iter,effiter_U,effiter_R,gl->CFL);
  wfprintf(datafile," nd=%ld ns=%ld nf=%ld",nd,ns,nf);
  wfprintf(datafile," is=%ld ie=%ld",gl->domain_all.is,gl->domain_all.ie);
#ifdef _2DL
  wfprintf(datafile," js=%ld je=%ld",gl->domain_all.js,gl->domain_all.je);
#endif
#ifdef _3DL
  wfprintf(datafile," ks=%ld ke=%ld",gl->domain_all.ks,gl->domain_all.ke);
#endif

#if defined(UNSTEADY)
  wfprintf(datafile," time=%E",gl->time);
#else
  wfprintf(datafile," time=%E",0.0); 
#endif

#ifdef UNSTEADY
  wfprintf(datafile," dt=%E",gl->dt);
#else
  wfprintf(datafile," dt=%E",dt_steady); 
#endif

#ifdef EMFIELD
  wfprintf(datafile," Lc=%E effiter_U_emfield=%E effiter_R_emfield=%E",gl->Lc,effiter_U_emfield,effiter_R_emfield);
#endif

  wfprintf(datafile,"\n");
#ifdef DISTMPI
  if(rank==0) {
    fseek(datafile,0,SEEK_END);
    offset=ftell(datafile);
  }
  MPI_Bcast(&offset,1,MPI_LONG,0,MPI_COMM_WORLD);
#endif


#ifdef DISTMPI
  wfclose(datafile);
  do {datafile = fopen(filename, "rb+");} while(datafile==NULL);
  domain=_domain_from_rank(rank,gl);
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    for (flux=0; flux<nf; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->U[flux];
  }
#endif

#ifdef DISTMPI
  fd=fileno(datafile);
  for_ijk(domain,is,js,ks,ie,je,ke){
    for (flux=0; flux<nf; flux++) fluxtmp[_ai(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->U[flux];
    pwrite(fd,fluxtmp[_ai(gl,i,j,k)],sizeof(flux_t),_ai_mpidatafile(gl,i,j,k)*sizeof(flux_t)+offset);
  }
  offset+=((gl->domain_all.ie-gl->domain_all.is+1)*(gl->domain_all.je-gl->domain_all.js+1)*(gl->domain_all.ke-gl->domain_all.ks+1)*sizeof(flux_t));
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    wfwrite(fluxtmp[_ai_all(gl,i,j,k)], sizeof(flux_t), 1, datafile);
  }
#endif


#ifdef EMFIELD
#ifndef DISTMPI
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    for (flux=0; flux<nfe; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->Uemfield[flux];
  }
#endif

#ifdef DISTMPI
  for_ijk(domain,is,js,ks,ie,je,ke){
    for (flux=0; flux<nfe; flux++) fluxtmp[_ai(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->Uemfield[flux];
    pwrite(fd,fluxtmp[_ai(gl,i,j,k)],sizeof(fluxemfield_t),_ai_mpidatafile(gl,i,j,k)*sizeof(fluxemfield_t)+offset);
  }
  offset+=((gl->domain_all.ie-gl->domain_all.is+1)*(gl->domain_all.je-gl->domain_all.js+1)*(gl->domain_all.ke-gl->domain_all.ks+1)*sizeof(fluxemfield_t));
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    wfwrite(fluxtmp[_ai_all(gl,i,j,k)], sizeof(fluxemfield_t), 1, datafile);
  }
#endif
#endif



#ifdef _RESTIME_STORAGE_TRAPEZOIDAL
#ifndef DISTMPI
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    for (flux=0; flux<nf; flux++) fluxtmp[_ai_all(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux];
  }
#endif

#ifdef DISTMPI
  for_ijk(domain,is,js,ks,ie,je,ke){
    for (flux=0; flux<nf; flux++) fluxtmp[_ai(gl,i,j,k)][flux]=np[_ai(gl,i,j,k)].bs->trapezoidalm1[flux];
    pwrite(fd,fluxtmp[_ai(gl,i,j,k)],sizeof(flux_t),_ai_mpidatafile(gl,i,j,k)*sizeof(flux_t)+offset);
  }
#else
  for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
    wfwrite(fluxtmp[_ai_all(gl,i,j,k)], sizeof(flux_t), 1, datafile);
  }
#endif
#endif

#ifdef DISTMPI
  fclose(datafile);
#else
  wfclose(datafile);
#endif
  wfprintf(stdout,"done.\n");
  free(fluxtmp);
}




void write_data_file(np_t *np, gl_t *gl){
  char *tmp_filename,*tmp_filename2;
#ifdef DISTMPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  tmp_filename=(char *)malloc((strlen(gl->output_filename)+10)*sizeof(char));
  tmp_filename2=(char *)malloc((strlen(gl->output_filename)+10)*sizeof(char));
  strcpy(tmp_filename,gl->output_filename);
  SOAP_strins(".wbak",&tmp_filename,strlen(tmp_filename));
  strcpy(tmp_filename2,gl->output_filename);
  SOAP_strins(".wbak2",&tmp_filename2,strlen(tmp_filename2));
#ifdef DISTMPI
  if (rank==0) {
#endif
    rename(tmp_filename,tmp_filename2);
    rename(gl->output_filename,tmp_filename);
    if (!is_file_writable(gl->output_filename)) fatal_error("Insufficient permissions to write data file %s.",gl->output_filename);
#ifdef DISTMPI
  }
#endif
  if (gl->OUTPUTASCII) {
    write_data_file_binary_ascii(gl->output_filename, np, gl, DATATYPE_ASCII);
  }
  else if (gl->OUTPUTINTERPOLATION){
    write_data_file_interpolation(gl->output_filename, np, gl);
  }
  else if (gl->OUTPUTBINARYMPI){
    write_data_file_mpi(gl->output_filename, np, gl);
  }
  else {
    write_data_file_binary_ascii(gl->output_filename, np, gl, DATATYPE_BINARY);
  }
  free(tmp_filename);
  free(tmp_filename2);
}



void read_data_file(input_t input, np_t *np, gl_t *gl){
#ifdef UNSTEADY
  long i,j,k;
  long flux;
#endif
 
  long spec;
  gl->nsinit=ns;
  for (spec=0; spec<ns; spec++) gl->initspecies[spec]=spec;
  
  if (input.READDATAFILE) {
    if (input.ASCII) {
      read_data_file_binary_ascii(input.name, np, gl, 0, DATATYPE_ASCII);
    }
    else if (input.INTERPOLATION){
      read_data_file_interpolation(input, np, gl);
    }
    else if (input.BINARYMPI) {
      read_data_file_mpi(input.name, np, gl, 0);
    }
    else {
      read_data_file_binary_ascii(input.name, np, gl, 0, DATATYPE_BINARY);
    }
    gl->iter=max(gl->iter,1);
    gl->INIT_FLUID_READ=TRUE;
    gl->INIT_EMFIELD_READ=TRUE;
  }
#ifdef UNSTEADY
    for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          for (flux=0; flux<nf; flux++){
            (np)[_ai(gl,i,j,k)].bs->Um1[flux]=(np)[_ai(gl,i,j,k)].bs->U[flux];
          }
#ifdef EMFIELD
          for (flux=0; flux<nfe; flux++){
            (np)[_ai(gl,i,j,k)].bs->Uemfieldm1[flux]=(np)[_ai(gl,i,j,k)].bs->Uemfield[flux];
          }
#endif
    }
      if (input.MM1) read_data_file_mpi(input.name_m1, np, gl, 1);
      if (input.M1) read_data_file_binary_ascii(input.name_m1, np, gl, 1, DATATYPE_BINARY);
#if _RESTIME_BW > 2  
      for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
            for (flux=0; flux<nf; flux++){
              (np)[_ai(gl,i,j,k)].bs->Um2[flux]=(np)[_ai(gl,i,j,k)].bs->Um1[flux];
            }
      }
      if (input.MM2) read_data_file_mpi(input.name_m2, np, gl, 2);
      if (input.M2) read_data_file_binary_ascii(input.name_m2, np, gl, 2, DATATYPE_BINARY);
#endif
#if _RESTIME_BW > 3 
      for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
            for (flux=0; flux<nf; flux++){
              (np)[_ai(gl,i,j,k)].bs->Um3[flux]=(np)[_ai(gl,i,j,k)].bs->Um2[flux];
            }
      }
      if (input.MM3) read_data_file_mpi(input.name_m3, np, gl, 3);
      if (input.M3) read_data_file_binary_ascii(input.name_m3, np, gl, 3, DATATYPE_BINARY);
#endif
#endif

}


void find_interpolation_weight(np_t *np, gl_t *gl, long l, dim_t x_file, dim_t dx1_file, 
#ifdef _2DL 
                               dim_t dx2_file,
#endif
#ifdef _3DL
                               dim_t dx3_file,
#endif
                               double radiusmax2, double *thisweight){
  double distance;
  EXM_mat_t mat1,mat2,mat3,mat1inv;
  long dim;
  distance=0.0;
  for (dim=0; dim<nd; dim++) distance+=sqr(_x(np[l],dim)-x_file[dim]);
  *thisweight=0.0;
  if (distance<radiusmax2) {
    EXM_init_matrix(&mat1, nd, nd);
    for (dim=0; dim<nd; dim++){
      mat1.cont[EXM_aim(mat1.glm,dim,0)]=dx1_file[dim];
#ifdef _2DL
      mat1.cont[EXM_aim(mat1.glm,dim,1)]=dx2_file[dim];
#endif
#ifdef _3DL
      mat1.cont[EXM_aim(mat1.glm,dim,2)]=dx3_file[dim];
#endif
    }
    EXM_init_matrix(&mat1inv, nd, nd);
    EXM_invert_matrix_analytical(mat1, &mat1inv);
    EXM_init_matrix(&mat2, nd, 1);
    for (dim=0; dim<nd; dim++){
      mat2.cont[EXM_aim(mat2.glm,dim,0)]=_x(np[l],dim)-x_file[dim];
    }
    EXM_init_matrix(&mat3, nd, 1);
    EXM_multiply_matrices(mat1inv, mat2, &mat3);

    *thisweight=0.0;
    //for (dim=0; dim<nd; dim++) thisweight=max(thisweight,fabs(mat3.cont[EXM_aim(mat3.glm,dim,0)]));
    for (dim=0; dim<nd; dim++) *thisweight+=fabs(pow(fabs(mat3.cont[EXM_aim(mat3.glm,dim,0)]),3.0));
    *thisweight=fabs(pow(*thisweight,1.0/3.0));
    *thisweight=max(1e-16,max(0.0003-(*thisweight)*0.00001,1.0-(*thisweight)));
    EXM_free_matrix(&mat1);
    EXM_free_matrix(&mat1inv);
    EXM_free_matrix(&mat2);
    EXM_free_matrix(&mat3);
  }
}


bool is_interpolation_occurring_in_zone(np_t *np, gl_t *gl, int TYPELEVEL, dim_t x_file, double radiusmax2, zone_t zone, long *i, long *j, long *k){
  long dim;
  double distance;
  bool FOUND;
  FOUND=FALSE;
//  fprintf(stderr,"zone.is=%ld zone.ie=%ld\n",zone.is,zone.ie);
  //if (zone.ie>gl->domain_all.ie) fatal_error("problem here..");
  //if (zone.is<gl->domain_all.is) fatal_error("problem here..");
  *i=zone.is-1;
  do {
    (*i)++;
    *j=zone.js-1;
    do {
      (*j)++;
#ifdef _3DL
      *k=zone.ks-1;
      do {
        (*k)++;
#endif
        if (is_node_valid(np[_ai(gl,*i,*j,*k)],TYPELEVEL)){
          distance=0.0;
          for (dim=0; dim<nd; dim++) distance+=sqr(_x(np[_ai(gl,*i,*j,*k)],dim)-x_file[dim]);
          if (distance<radiusmax2) {
          //if (zone.is==zone.ie) fprintf(stderr,"[[%E %E %ld,%ld,%ld %E %E]]\n",distance,radiusmax2,*i,*j,*k,_x(np[_ai(gl,*i,*j,*k)],0),x_file[0]);
            FOUND=TRUE;
          }
        } 
#ifdef _3DL
      } while(!FOUND && *k<zone.ke);
#endif
    } while(!FOUND && *j<zone.je);
  } while(!FOUND && *i<zone.ie);
  return(FOUND);
}


bool find_interpolation_zone(np_t *np, gl_t *gl, int TYPELEVEL, dim_t x_file, double radiusmax2, zone_t *zone){
  long i,j,k,offset,imem;
  bool FOUNDWITHIN,FOUNDLEFT,FOUNDRIGHT,FOUND;
  zone_t istationzone,domain_eff;

  domain_eff=_zone_intersection(gl->domain_all,gl->domain_lim);
  FOUNDWITHIN=FALSE;
  FOUNDLEFT=FALSE;
  FOUNDRIGHT=FALSE;
  istationzone=domain_eff;
  if (is_interpolation_occurring_in_zone(np,gl,TYPELEVEL,x_file,radiusmax2,*zone,&i,&j,&k)){
    istationzone.is=i;
    istationzone.ie=i;
    FOUNDWITHIN=TRUE;
  } else {
    offset=0;
    do {
      offset++;
      if (zone->ie+offset<=domain_eff.ie) {
        istationzone.ie=zone->ie+offset;
        istationzone.is=zone->ie+offset;
        if (is_interpolation_occurring_in_zone(np,gl,TYPELEVEL,x_file,radiusmax2,istationzone,&i,&j,&k)){
          FOUNDLEFT=TRUE;
        }
      }
      if (zone->is-offset>=domain_eff.is && !FOUNDLEFT) {
        istationzone.ie=zone->is-offset;
        istationzone.is=zone->is-offset;
        if (is_interpolation_occurring_in_zone(np,gl,TYPELEVEL,x_file,radiusmax2,istationzone,&i,&j,&k)){
          FOUNDRIGHT=TRUE;
        }
      }
    } while (!FOUNDLEFT && !FOUNDRIGHT && offset<=(domain_eff.ie-domain_eff.is+1) );
  } 

  if (FOUNDRIGHT || FOUNDLEFT || FOUNDWITHIN){
   imem=istationzone.is;

   if (FOUNDRIGHT){
    zone->ie=imem;
   } else {
    FOUND=TRUE;
    while (istationzone.ie<domain_eff.ie && FOUND) {
      istationzone.ie++;
      istationzone.is++;
      FOUND=is_interpolation_occurring_in_zone(np,gl,TYPELEVEL,x_file,radiusmax2,istationzone,&i,&j,&k);
    } 
    if (!FOUND) zone->ie=istationzone.ie-1; else zone->ie=domain_eff.ie; 
   } 

   if (FOUNDLEFT){
    zone->is=imem;
   } else {
    istationzone.is=imem;
    istationzone.ie=imem;
    FOUND=TRUE;
    while (istationzone.is>domain_eff.is && FOUND) {
      istationzone.is--;
      istationzone.ie--;
      FOUND=is_interpolation_occurring_in_zone(np,gl,TYPELEVEL,x_file,radiusmax2,istationzone,&i,&j,&k);
    } 
    if (!FOUND) zone->is=istationzone.is+1; else zone->is=domain_eff.is;
   } 
  } 
  return(FOUNDRIGHT || FOUNDLEFT || FOUNDWITHIN);
}


bool is_data_point_in_domain(dim_t x_file, dim_t xmin, dim_t xmax, double radiusmax2){
  bool INDOMAIN;
  long dim;
  INDOMAIN=TRUE;
  for (dim=0; dim<nd; dim++) {
    if (INDOMAIN) {
      if (x_file[dim]<xmin[dim] && sqr(x_file[dim]-xmin[dim])>radiusmax2) INDOMAIN=FALSE;
    }
  }
  for (dim=0; dim<nd; dim++) {
    if (INDOMAIN) {
      if (x_file[dim]>xmax[dim] && sqr(x_file[dim]-xmax[dim])>radiusmax2) INDOMAIN=FALSE;
    }
  }
  return(INDOMAIN);
}


void find_interpolation_map(char* initvar_str_userinput, char* initvar_str_file, long numvars, long* numvars_file, mapvar_t map[], int INITVARTYPE){
  long cnt, cntr, cntp, numvars_userinput; 
  bool FOUND_FILE, FOUND_USERINPUT, FOUND_WDEFAULT, FOUND_NAME;
  char **initvar_names_file = NULL;
  char **initvar_names_userinput = NULL;
  initvarname_t *initvar_names;
  *numvars_file = 0;

  initvar_names = (initvarname_t *) malloc(sizeof(initvarname_t) * numvars * sizeof(char));
  switch (INITVARTYPE) {
    case INITVARTYPE_FLUID:
      find_default_initvar_name(initvar_names); break;
    case INITVARTYPE_EMFIELD:
#ifdef EMFIELD
      find_default_initvar_name_emfield(initvar_names); break;  
#endif 
    default: fatal_error("INITVARTYPE must be either INITVARTYPE_FLUID or INITVARTYPE_EMFIELD.");
    }

  find_words_from_string(initvar_str_file, " ", &initvar_names_file, numvars_file);
  find_words_from_string(initvar_str_userinput, ",= ", &initvar_names_userinput, &numvars_userinput);

  /* loop through the current set of variable names */
  for (cnt=0; cnt<numvars; cnt++) {
    map[cnt].newvalue=0;
    map[cnt].newindex=-1;
    FOUND_FILE = FALSE;
    FOUND_USERINPUT = FALSE;
    FOUND_WDEFAULT = FALSE;


    /* 1. loop to match current and read variable names, then map */
    for (cntr=0; cntr<*numvars_file; cntr++) {
      if (strcmp(initvar_names[cnt], initvar_names_file[cntr])==0){
        FOUND_FILE=TRUE;
        map[cnt].newindex=cntr;
      } 
    }

    /* 2a. loop to match input and current variables, then map */
    for (cntr=0; cntr<numvars_userinput; cntr+=2) {
      if (strcmp(initvar_names[cnt],initvar_names_userinput[cntr])==0) {
        if (FOUND_USERINPUT) fatal_error("Problem reading string after -imap. Variable %s can only be declared once.",initvar_names[cnt]);
        FOUND_USERINPUT=TRUE;
        if (sscanf(initvar_names_userinput[cntr+1],"%lg",&map[cnt].newvalue)==1) {
          // found a double value to be assigned to the variable
          map[cnt].newindex=-1;
        } else {
          /* 2b. a variable name is expected because no double value is assigned  */
          FOUND_NAME=FALSE;
          for (cntp=0; cntp<*numvars_file; cntp++) {
            if (strcmp(initvar_names_file[cntp],initvar_names_userinput[cntr+1])==0) {
              FOUND_NAME=TRUE;
              map[cnt].newindex=cntp;
            } 
          }
          if (!FOUND_NAME) fatal_error("Problem reading assigned value to %s in -imap input. Not a valid double or datafile variable name.",initvar_names[cnt]);      
        }   
      }
    }

    /* 3. check if w_default is present and assign its value to species not found in interpolation file */
    if (!FOUND_USERINPUT && !FOUND_FILE) {
      for (cntr=0; cntr<numvars_userinput; cntr+=2) {
        if (strcmp("w_default",initvar_names_userinput[cntr])==0 && strncmp(initvar_names[cnt],"w_",2)==0) {
          FOUND_WDEFAULT = TRUE;
          if (sscanf(initvar_names_userinput[cntr+1],"%lg",&map[cnt].newvalue)!=1) fatal_error("Problem reading assigned w_default value in -imap input. Not a valid double.");
        }
      }
    }

    /* check if a variable within memory is not present in read variables and was not given a value by user
      (it wasn't found in input, nor in the interpolation file variables, and isn't a species assigned to w_default) */
    if (!FOUND_USERINPUT && !FOUND_FILE && !FOUND_WDEFAULT) fatal_error("Variable %s was not found in the data file and was not given a value after the -imap flag.",initvar_names[cnt]);  
      
    if (map[cnt].newindex==-1) wfprintf(stdout,"%s=%lg..", initvar_names[cnt], map[cnt].newvalue); 
  } 

  for (cnt=0; cnt<*numvars_file; cnt++) free(initvar_names_file[cnt]); 
  for (cnt=0; cnt<numvars_userinput; cnt++) free(initvar_names_userinput[cnt]);
  free(initvar_names_file); 
  free(initvar_names_userinput); 
  free(initvar_names);
}


double _map_interpolation_var(initvar_t initvars_file, mapvar_t map[], long count){
  double tmp;
  tmp = map[count].newvalue;
  if (map[count].newindex!=-1) tmp = initvars_file[map[count].newindex];
  return(tmp);
}


void read_data_file_interpolation(input_t input, np_t *np, gl_t *gl){  
  FILE *datafile;
  char data_format_str[100], initvar_fluid_str_file[500], initvar_emfield_str_file[500];
  long i,j,k,l_file,cnt,dim,cntzone;
  long numsubzone, numflux_read,numspec_read,numdim_read,numnodes,numvars_fluid_file;
  double tmp_dt,tmp_time;
  double *weight,*radiusmax2_file,thisweight;
  zone_t *subzone;
  dim_t *xmin,*xmax;
  initvar_t *initvar;
  //initvar_t *initvar_file;
  zone_t zone;
  long numsubzone_desired;
#ifdef EMFIELD
  bool ISREAD_EMFIELD = FALSE;
  long numvars_emfield_file;
  initvar_emfield_t *initvar_emfield;
 // initvar_emfield_t *initvar_emfield_file;
#endif
  bool FORMAT001;
  dim_t *dx1_file,*x_file;
#ifdef _2DL
  dim_t *dx2_file;
#endif
#ifdef _3DL
  dim_t *dx3_file;
#endif
  int cnterror;
#ifdef OPENMPTHREADS
  omp_lock_t *nodelock;
#endif
#ifdef DISTMPI
  int rank,numproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  weight=(double *)malloc(sizeof(double)*(gl->domain_lim.ie+4) 
#ifdef _2DL 
    *(gl->domain_lim.je+4)
#endif
#ifdef _3DL
    *(gl->domain_lim.ke+4)
#endif
    );

#ifdef OPENMPTHREADS
  nodelock=(omp_lock_t *)malloc(sizeof(double)*(gl->domain_lim.ie+4) 
#ifdef _2DL 
    *(gl->domain_lim.je+4)
#endif
#ifdef _3DL
    *(gl->domain_lim.ke+4)
#endif
    );
#endif

  datafile = fopen(input.name, "r");
  if (datafile==NULL)
    fatal_error("Having problems opening interpolation datafile %s.",input.name);

/* first do the fluid properties */

  mapvar_t *map_fluid = (mapvar_t*)malloc(numinitvar * sizeof(mapvar_t));
#ifdef EMFIELD  
  mapvar_t *map_emfield = (mapvar_t*)malloc(numinitvar_emfield * sizeof(mapvar_t));
#endif  
  initvar=(initvar_t *)malloc(sizeof(initvar_t)*(gl->domain_lim.ie+4) 
#ifdef _2DL 
    *(gl->domain_lim.je+4)
#endif
#ifdef _3DL
    *(gl->domain_lim.ke+4)
#endif
    );


  for (cnt=0; cnt<16; cnt++){
    if (fscanf(datafile,"%c",&(data_format_str[cnt]))!=1) {
      fatal_error("Problem with fscanf in read_data_file_interpolation().");
    }
  }
  data_format_str[16]=EOS;
  wfprintf(stdout,"Reading interpolation data file %s..",input.name);
  FORMAT001=FALSE;  
  if (strcmp("WARPINTFORMAT001",data_format_str)==0) {
//    wfprintf(stdout,"in CFDWARP format 001..\n");
    FORMAT001=TRUE;
  }

  if (FORMAT001) {
      if (fscanf(datafile," numnodes=%ld nf=%ld nd=%ld ns=%ld windowis=%ld windowie=%ld iter=%ld effiter_U=%lg effiter_R=%lg CFL=%lg time=%lg dt=%lg vars_fluid=\"%[^\"]\" vars_emfield=\"%[^\"]\"%*[^\n]",
             &numnodes,&numflux_read,&numdim_read,&numspec_read,&(gl->window.is),&(gl->window.ie),
             &(gl->iter),&(gl->effiter_U),&(gl->effiter_R),&(gl->CFL),&(tmp_time),&tmp_dt,initvar_fluid_str_file,initvar_emfield_str_file)!=14) fatal_error("Problem reading interpolation data file.");

        find_interpolation_map(input.interpolationvarsmap, initvar_fluid_str_file, numinitvar, &numvars_fluid_file, map_fluid, INITVARTYPE_FLUID);
#ifdef EMFIELD
        find_interpolation_map(input.interpolationvarsmap, initvar_emfield_str_file, numinitvar_emfield, &numvars_emfield_file, map_emfield, INITVARTYPE_EMFIELD);
        if (strcmp(initvar_emfield_str_file,"NONE")!=0) ISREAD_EMFIELD = TRUE;
#endif            

#ifdef UNSTEADY
      gl->time=tmp_time;
      gl->dt=tmp_dt;
#endif

    fgetc(datafile);
    if (numdim_read!=nd) fatal_error("Number of dimensions read (%ld) does not equal current number of dimensions (%ld).",numdim_read,nd);
   // if (numspec_read!=ns) fatal_error("Number of species read (%ld) does not equal current number of species (%ld).",numspec_read,ns);
   // if (numflux_read!=nf) fatal_error("Number of fluxes read (%ld) does not equal current number of fluxes (%ld).",numflux_read,nf);
  } else {
    fatal_error("Interpolation file format unknown.");
  }

  /* read data and store in ram */
  //initvar_file=(initvar_t *)malloc(numnodes*sizeof(initvar_t));
  double(*initvar_file)[numvars_fluid_file] = malloc(numnodes*numvars_fluid_file*sizeof(double));
  x_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
  dx1_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
#ifdef _2DL
  dx2_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
#endif
#ifdef _3DL
  dx3_file=(dim_t *)malloc(numnodes*sizeof(dim_t));
#endif
  radiusmax2_file=(double *)malloc(numnodes*sizeof(double));
  for (l_file=0; l_file<numnodes; l_file++){
    cnterror=0;
    if (fread(initvar_file[l_file], numvars_fluid_file*sizeof(double), 1, datafile)!=1) cnterror++;
    if (fread(x_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (fread(dx1_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
#ifdef _2DL
    if (fread(dx2_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
#endif
#ifdef _3DL
    if (fread(dx3_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
#endif
    if (cnterror>0) fatal_error("Could not read all data properly.");
    radiusmax2_file[l_file]=0.0e0;
    for (dim=0; dim<nd; dim++) 
      radiusmax2_file[l_file]+=sqr(fabs(dx1_file[l_file][dim])
#ifdef _2DL
        +fabs(dx2_file[l_file][dim])
#endif
#ifdef _3DL
        +fabs(dx3_file[l_file][dim])
#endif
      );
    radiusmax2_file[l_file]*=1.1;

  }

  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        weight[_ai(gl,i,j,k)]=0.0e0;
#ifdef OPENMPTHREADS
        omp_init_lock(&(nodelock[_ai(gl,i,j,k)]));
#endif
        for (cnt=0; cnt<numinitvar; cnt++) (initvar[_ai(gl,i,j,k)])[cnt]=0.0;
  }

  zone=_zone_intersection(gl->domain_all,gl->domain_lim);

  subzone=(zone_t *)malloc(sizeof(zone_t));
  find_subzones_in_zone_given_zonelength(SUBZONE_DESIRED_WIDTH, zone, &numsubzone, &subzone);

#ifdef OPENMPTHREADS
  numsubzone_desired=MIN_NUMSUBZONE_PER_THREAD*omp_get_max_threads();
#else
  numsubzone_desired=MIN_NUMSUBZONE_PER_THREAD;
#endif
  if (numsubzone<numsubzone_desired)
    find_subzones_in_zone_given_numsubzone(zone, numsubzone_desired, &numsubzone, &subzone);

  xmin=(dim_t *)malloc(numsubzone*sizeof(dim_t));
  xmax=(dim_t *)malloc(numsubzone*sizeof(dim_t));

  for (cntzone=0; cntzone<numsubzone; cntzone++){
    for (dim=0; dim<nd; dim++){
      xmin[cntzone][dim]=1e99;
      xmax[cntzone][dim]=-1e99;
    }
    for_ijk(subzone[cntzone],is,js,ks,ie,je,ke){
          if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)){
            for (dim=0; dim<nd; dim++){
              xmin[cntzone][dim]=min(xmin[cntzone][dim],_x(np[_ai(gl,i,j,k)],dim));
              xmax[cntzone][dim]=max(xmax[cntzone][dim],_x(np[_ai(gl,i,j,k)],dim));
            }
          }
    }
  }

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  wfprintf(stdout,"Fluid/%ld",numsubzone*numproc);
#else
  wfprintf(stdout,"Fluid/%ld",numsubzone);
#endif

#if defined(OPENMPTHREADS) 
#pragma omp parallel for private(l_file,cntzone,dim,zone,i,j,k,cnt,thisweight) schedule(dynamic) 
#endif
  for (cntzone=0; cntzone<numsubzone; cntzone++){
    for (l_file=0; l_file<numnodes; l_file++){
      if (is_data_point_in_domain(x_file[l_file],xmin[cntzone],xmax[cntzone],radiusmax2_file[l_file])){
        zone=subzone[cntzone];
        if (find_interpolation_zone(np,gl,TYPELEVEL_FLUID,x_file[l_file],radiusmax2_file[l_file],&zone)){
          for_jik(zone,is,js,ks,ie,je,ke){
                if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)){
                  find_interpolation_weight(np,gl,_ai(gl,i,j,k),x_file[l_file],dx1_file[l_file],
#ifdef _2DL
                    dx2_file[l_file],
#endif
#ifdef _3DL
                    dx3_file[l_file],
#endif
                    radiusmax2_file[l_file],&thisweight);
#ifdef OPENMPTHREADS 
                  omp_set_lock(&(nodelock[_ai(gl,i,j,k)]));
#endif
                  if (thisweight>1e-99) {
                    weight[_ai(gl,i,j,k)]+=thisweight;
                    for (cnt=0; cnt<numinitvar; cnt++) 
                      initvar[_ai(gl,i,j,k)][cnt]+=thisweight*_map_interpolation_var(initvar_file[l_file],map_fluid,cnt);
                  }
#ifdef OPENMPTHREADS 
                  omp_unset_lock(&(nodelock[_ai(gl,i,j,k)]));
#endif
                }
          }
        }
      } 
    }
    fprintf(stdout,".");
    fflush(stdout);
//    if (mod(cntzone,numsubzone/100+1)==0) wfprintf(stdout,".");
  }

#ifdef OPENMPTHREADS
  #pragma omp parallel for private(i,j,k,cnt) schedule(static) 
#endif
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
	  if (weight[_ai(gl,i,j,k)]>1e-99 && is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID)) {
            for (cnt=0; cnt<numinitvar; cnt++) initvar[_ai(gl,i,j,k)][cnt]=initvar[_ai(gl,i,j,k)][cnt]/weight[_ai(gl,i,j,k)];
            init_node_fluid(np,_ai(gl,i,j,k), gl, defaultinitvartypefluid, initvar[_ai(gl,i,j,k)]);
            np[_ai(gl,i,j,k)].INIT_FLUID=TRUE;
          }
  }

  free(initvar);
  free(initvar_file);
  free(map_fluid);

/* second do the emfield properties */
#ifdef EMFIELD

  initvar_emfield=(initvar_emfield_t *)malloc(sizeof(initvar_emfield_t)*(gl->domain_lim.ie+4) 
#ifdef _2DL 
    *(gl->domain_lim.je+4)
#endif
#ifdef _3DL
    *(gl->domain_lim.ke+4)
#endif
    );

if (ISREAD_EMFIELD) {
  for (cnt=0; cnt<16; cnt++){
    if (fscanf(datafile,"%c",&(data_format_str[cnt]))!=1){
      fatal_error("Problem with fscanf in emfield part of read_data_file_interpolation().");
    }
  }
  data_format_str[16]=EOS;
  FORMAT001=FALSE;  
  if (strcmp("WARPINTFORMAT001",data_format_str)==0) {
    FORMAT001=TRUE;
  }

  if (FORMAT001) {
    if (fscanf(datafile," numnodes_emfield=%ld nfe=%ld nd=%ld Lc=%lg effiter_U_emfield=%lg effiter_R_emfield=%lg%*[^\n]",
             &numnodes,&numflux_read,&numdim_read,&(gl->Lc),&(gl->effiter_U_emfield),&(gl->effiter_R_emfield))!=6){
      fatal_error("Problem reading emfield preambule in interpolating file.");
    }
    fgetc(datafile);
    if (numdim_read!=nd) fatal_error("Number of dimensions read (%ld) does not equal current number of dimensions (%ld).",numdim_read,nd);
    //if (numflux_read!=nfe) fatal_error("Number of fluxes read (%ld) does not equal current number of emfield fluxes (%ld).",numflux_read,nfe);
    gl->Lc=1.0e0;

  } else {
    fatal_error("Interpolation file format unknown for EMfield variables.");
  }
}
  /* read data and store in ram */
  double(*initvar_emfield_file)[numvars_emfield_file] = malloc(numnodes*numvars_emfield_file*sizeof(double));
  //initvar_emfield_file=(initvar_emfield_t *)malloc(numnodes*sizeof(initvar_emfield_t));
if (ISREAD_EMFIELD) {
  x_file=(dim_t *)realloc(x_file,numnodes*sizeof(dim_t));
  dx1_file=(dim_t *)realloc(dx1_file,numnodes*sizeof(dim_t));
#ifdef _2DL
  dx2_file=(dim_t *)realloc(dx2_file,numnodes*sizeof(dim_t));
#endif
#ifdef _3DL
  dx3_file=(dim_t *)realloc(dx3_file,numnodes*sizeof(dim_t));
#endif
  radiusmax2_file=(double *)realloc(radiusmax2_file,numnodes*sizeof(double));

  for (l_file=0; l_file<numnodes; l_file++){
    cnterror=0;
    if (fread(initvar_emfield_file[l_file], numvars_emfield_file*sizeof(double), 1, datafile)!=1) cnterror++;
    if (fread(x_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
    if (fread(dx1_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
#ifdef _2DL
    if (fread(dx2_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
#endif
#ifdef _3DL
    if (fread(dx3_file[l_file], sizeof(dim_t), 1, datafile)!=1) cnterror++;
#endif
    if (cnterror>0) fatal_error("Could not read all data properly.");
    radiusmax2_file[l_file]=0.0e0;
    for (dim=0; dim<nd; dim++) 
      radiusmax2_file[l_file]+=sqr(fabs(dx1_file[l_file][dim])
#ifdef _2DL
        +fabs(dx2_file[l_file][dim])
#endif
#ifdef _3DL
        +fabs(dx3_file[l_file][dim])
#endif
      );
    radiusmax2_file[l_file]*=1.1;
  }
}
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        weight[_ai(gl,i,j,k)]=0.0e0;
        for (cnt=0; cnt<numinitvar_emfield; cnt++) (initvar_emfield[_ai(gl,i,j,k)])[cnt]=0.0;
  }

  for (cntzone=0; cntzone<numsubzone; cntzone++){
    for (dim=0; dim<nd; dim++){
      xmin[cntzone][dim]=1e99;
      xmax[cntzone][dim]=-1e99;
    }
    for_ijk(subzone[cntzone],is,js,ks,ie,je,ke){
          if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
            for (dim=0; dim<nd; dim++){
              xmin[cntzone][dim]=min(xmin[cntzone][dim],_x(np[_ai(gl,i,j,k)],dim));
              xmax[cntzone][dim]=max(xmax[cntzone][dim],_x(np[_ai(gl,i,j,k)],dim));
            }
          }
    }
  }

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  wfprintf(stdout,"EMfield/%ld",numsubzone*numproc);
#else
  wfprintf(stdout,"EMfield/%ld",numsubzone);
#endif

#if defined(OPENMPTHREADS) //&& !defined(DISTMPI)
#pragma omp parallel for private(l_file,cntzone,cnt,thisweight,dim,zone,i,j,k) schedule(dynamic) 
#endif
  for (cntzone=0; cntzone<numsubzone; cntzone++){
    for (l_file=0; l_file<numnodes; l_file++){
      if (is_data_point_in_domain(x_file[l_file],xmin[cntzone],xmax[cntzone],radiusmax2_file[l_file])){
        zone=subzone[cntzone];
        if (find_interpolation_zone(np,gl,TYPELEVEL_EMFIELD,x_file[l_file],radiusmax2_file[l_file],&zone)){
          for_jik(zone,is,js,ks,ie,je,ke){
                if (is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
                  find_interpolation_weight(np,gl,_ai(gl,i,j,k),x_file[l_file],dx1_file[l_file],
#ifdef _2DL
                                            dx2_file[l_file],
#endif
#ifdef _3DL
                                            dx3_file[l_file],
#endif
                                            radiusmax2_file[l_file],&thisweight);
#ifdef OPENMPTHREADS 
                  omp_set_lock(&(nodelock[_ai(gl,i,j,k)]));
#endif
                  if (thisweight>1e-99) {
                    weight[_ai(gl,i,j,k)]+=thisweight;
                    for (cnt=0; cnt<numinitvar_emfield; cnt++) 
                      initvar_emfield[_ai(gl,i,j,k)][cnt]+=thisweight*_map_interpolation_var(initvar_emfield_file[l_file],map_emfield,cnt);
                  }
#ifdef OPENMPTHREADS 
                  omp_unset_lock(&(nodelock[_ai(gl,i,j,k)]));
#endif
                }
          }
        }
      }
    }
//    if (mod(cntzone,numsubzone/100+1)==0) wfprintf(stdout,".");
    fprintf(stdout,".");
    fflush(stdout);
  }

#ifdef OPENMPTHREADS
#pragma omp parallel for private(i,j,k,cnt) schedule(static) 
#endif
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        if (weight[_ai(gl,i,j,k)]>1e-99 && is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)) {
          for (cnt=0; cnt<numinitvar_emfield; cnt++) initvar_emfield[_ai(gl,i,j,k)][cnt]=initvar_emfield[_ai(gl,i,j,k)][cnt]/weight[_ai(gl,i,j,k)];
          init_node_emfield(np[_ai(gl,i,j,k)], gl, defaultinitvartypeemfield, initvar_emfield[_ai(gl,i,j,k)]);
          np[_ai(gl,i,j,k)].INIT_EMFIELD=TRUE;
        }
  }
  free(initvar_emfield);
  free(map_emfield);
  free(initvar_emfield_file);
#endif //EMFIELD
  free(subzone);
  free(xmin);
  free(xmax);
  fclose(datafile);
  free(weight);
#ifdef OPENMPTHREADS
  for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
        omp_destroy_lock(&(nodelock[_ai(gl,i,j,k)]));
  }
  free(nodelock);
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank!=0) {
    gl->effiter_U=0.0;
    gl->effiter_R=0.0;
    #ifdef EMFIELD
    gl->effiter_U_emfield=0.0;
    gl->effiter_R_emfield=0.0;   
    #endif
  }
#endif
  wfprintf(stdout,"done;\n");

  free(x_file);
  free(dx1_file);
#ifdef _2DL
  free(dx2_file);
#endif
#ifdef _3DL
  free(dx3_file);
#endif
  free(radiusmax2_file);
}


void write_data_file_interpolation(char *filename, np_t *np, gl_t *gl){
  FILE *datafile;
  long i,j,k,cnt;
  initvarname_t *initvar_names;
  dim_t dx1,x;
#ifdef _2DL
  dim_t dx2;
#endif
#ifdef _3DL
  dim_t dx3;
#endif
  double tmp_time, tmp_dt;
  long numnodes,dim;
  int TYPELEVEL,pass,passmax;
  bool *NODEVALID;
  initvar_t initvar;
  double effiter_U,effiter_R;
#ifdef EMFIELD
  double effiter_U_emfield,effiter_R_emfield;
  initvar_emfield_t initvar_emfield;
#endif
#ifdef DISTMPI
  int rank;
  MPI_Status MPI_Status1;
#endif

  /* nodes may be suspended. Hence, ensure that appropriate nodes are
     resumed. */
  resume_nodes_in_zone(np,gl,gl->domain);

  NODEVALID=(bool *)malloc(sizeof(bool)*(gl->domain_lim_all.ie-gl->domain_lim_all.is+1) 
#ifdef _2DL 
    *(gl->domain_lim_all.je-gl->domain_lim_all.js+1)
#endif
#ifdef _3DL
    *(gl->domain_lim_all.ke-gl->domain_lim_all.ks+1)
#endif
  );

  effiter_U=gl->effiter_U;
  effiter_R=gl->effiter_R;
#ifdef EMFIELD
  effiter_U_emfield=gl->effiter_U_emfield;
  effiter_R_emfield=gl->effiter_R_emfield;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Allreduce(&gl->effiter_U, &effiter_U, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&gl->effiter_R, &effiter_R, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef EMFIELD
  MPI_Allreduce(&gl->effiter_U_emfield, &effiter_U_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&gl->effiter_R_emfield, &effiter_R_emfield, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif
  datafile = wfopen(filename, "w");
  wfprintf(stdout,"Writing to CFDWARP interpolation data file %s..",filename);

#ifdef EMFIELD
  passmax=2;
#else
  passmax=1;
#endif
  for (pass=1; pass<=passmax; pass++){
    if (pass==1){
      TYPELEVEL=TYPELEVEL_FLUID;
    } else {
#ifdef EMFIELD
      TYPELEVEL=TYPELEVEL_EMFIELD;
#endif
    }
#ifdef DISTMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    find_NODEVALID_on_domain_all(np, gl, TYPELEVEL, NODEVALID);

    numnodes=0;
    for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
	  if (NODEVALID[_ai_all(gl,i,j,k)]) {
            numnodes++;
	  }
    }

    if (pass==1){
#ifdef UNSTEADY
      tmp_time=gl->time;
      tmp_dt=gl->dt;
#else
      tmp_time=0.0;
      tmp_dt=dt_steady;
#endif
      wfprintf(datafile,"WARPINTFORMAT001 numnodes=%ld nf=%ld nd=%ld ns=%ld windowis=%ld windowie=%ld iter=%ld effiter_U=%E effiter_R=%E CFL=%E time=%E dt=%E ",numnodes,nf,nd, ns,gl->window.is,gl->window.ie,gl->iter,effiter_U,effiter_R,gl->CFL,tmp_time,tmp_dt);
      wfprintf(datafile,"vars_fluid=\"");
      initvar_names = (initvarname_t *) malloc(sizeof(initvarname_t) * numinitvar * sizeof(char));   
      find_default_initvar_name(initvar_names); 
      for (cnt=0; cnt<numinitvar; cnt++) {
        wfprintf(datafile,"%s ", initvar_names[cnt]);   
      }
      wfprintf(datafile,"\"");   
      wfprintf(datafile," vars_emfield=\"");
#ifdef EMFIELD           
      initvar_names = (initvarname_t *) malloc(sizeof(initvarname_t) * numinitvar_emfield * sizeof(char));   
      find_default_initvar_name_emfield(initvar_names);            
      for (cnt=0; cnt<numinitvar_emfield; cnt++) {
        wfprintf(datafile,"%s ", initvar_names[cnt]);   
      }
      wfprintf(datafile,"\"\n");
#else
      wfprintf(datafile,"NONE\"\n");
#endif
    } else {
#ifdef EMFIELD
      wfprintf(datafile,"WARPINTFORMAT001 numnodes_emfield=%ld nfe=%ld nd=%ld Lc=%E effiter_U_emfield=%E effiter_R_emfield=%E\n",numnodes,nfe,nd,gl->Lc,effiter_U_emfield,effiter_R_emfield);
#endif
    }

    for_ijk(gl->domain_all,is,js,ks,ie,je,ke){
#ifdef DISTMPI
          if (pass==1){
            if (_node_rank(gl,i,j,k)==rank) {
	      if (NODEVALID[_ai_all(gl,i,j,k)]) {
	        find_default_initvar(np, gl, _ai(gl,i,j,k), initvar); 
              } else {
                for (cnt=0; cnt<numinitvar; cnt++) initvar[cnt]=0.0;
	      }
              if (rank!=0) {
                MPI_Ssend(initvar,numinitvar,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
              }
            }
            if (rank==0 && _node_rank(gl,i,j,k)!=0){
              MPI_Recv(initvar,numinitvar,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
            }
          } else {
#ifdef EMFIELD
            if (_node_rank(gl,i,j,k)==rank) {
	      if (NODEVALID[_ai_all(gl,i,j,k)]) {
	        find_default_initvar_emfield(np, gl, _ai(gl,i,j,k),initvar_emfield); 
              } else {
                for (cnt=0; cnt<numinitvar_emfield; cnt++) initvar_emfield[cnt]=0.0;
	      }
              if (rank!=0) {
                MPI_Ssend(initvar_emfield,numinitvar_emfield,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
              }
            }
            if (rank==0 && _node_rank(gl,i,j,k)!=0){
              MPI_Recv(initvar_emfield,numinitvar_emfield,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
            }
#endif
          }
#else
          if (pass==1){
            if (NODEVALID[_ai_all(gl,i,j,k)]) {
	      find_default_initvar(np, gl, _ai(gl,i,j,k), initvar); 
            } else {
              for (cnt=0; cnt<numinitvar; cnt++) initvar[cnt]=0.0;
            }
          } else {
#ifdef EMFIELD
            if (NODEVALID[_ai_all(gl,i,j,k)]) {
              find_default_initvar_emfield(np, gl, _ai(gl,i,j,k), initvar_emfield); 
            } else {
              for (cnt=0; cnt<numinitvar_emfield; cnt++) initvar_emfield[cnt]=0.0;
	    }
#endif
          }
#endif

          if (NODEVALID[_ai_all(gl,i,j,k)]) {
#ifdef DISTMPI
            if (_node_rank(gl,i,j,k)==rank) {
#endif
              for (dim=0; dim<nd; dim++) x[dim]=_x(np[_ai(gl,i,j,k)],dim);
              for (dim=0; dim<nd; dim++){
	        if ((i<gl->domain_all.ie && NODEVALID[_ai_all(gl,i+1,j,k)]) && (i>gl->domain_all.is && NODEVALID[_ai_all(gl,i-1,j,k)])) {
	          dx1[dim]=0.5*(np[_ai(gl,i+1,j,k)].bs->x[dim]-np[_ai(gl,i-1,j,k)].bs->x[dim]);
                } else {
                  if (i<gl->domain_all.ie && NODEVALID[_ai_all(gl,i+1,j,k)]) {
                    dx1[dim]=(np[_ai(gl,i+1,j,k)].bs->x[dim]-np[_ai(gl,i,j,k)].bs->x[dim]);
                  } else {
                    if (i>gl->domain_all.is && NODEVALID[_ai_all(gl,i-1,j,k)]) {
	              dx1[dim]=(np[_ai(gl,i,j,k)].bs->x[dim]-np[_ai(gl,i-1,j,k)].bs->x[dim]);
		    } else {
                      fatal_error("Couldn't find adjacent valid node along i needed for interpolation.");
		    }
                  }
	        }
#ifdef _2DL
                if ((j<gl->domain_all.je && NODEVALID[_ai_all(gl,i,j+1,k)]) && (j>gl->domain_all.js && NODEVALID[_ai_all(gl,i,j-1,k)])) {
                  dx2[dim]=0.5*(np[_ai(gl,i,j+1,k)].bs->x[dim]-np[_ai(gl,i,j-1,k)].bs->x[dim]);
                } else {
                  if (j<gl->domain_all.je && NODEVALID[_ai_all(gl,i,j+1,k)]) {
	            dx2[dim]=(np[_ai(gl,i,j+1,k)].bs->x[dim]-np[_ai(gl,i,j,k)].bs->x[dim]);
                  } else {
                    if (j>gl->domain_all.js && NODEVALID[_ai_all(gl,i,j-1,k)]) {
                      dx2[dim]=(np[_ai(gl,i,j,k)].bs->x[dim]-np[_ai(gl,i,j-1,k)].bs->x[dim]);
                    } else {
                      fatal_error("Couldn't find adjacent valid node along j needed for interpolation.");
		    }
                  }
                }
#endif
#ifdef _3DL
	        if ((k<gl->domain_all.ke && NODEVALID[_ai_all(gl,i,j,k+1)]) && (k>gl->domain_all.ks && NODEVALID[_ai_all(gl,i,j,k-1)])) {
                  dx3[dim]=0.5*(np[_ai(gl,i,j,k+1)].bs->x[dim]-np[_ai(gl,i,j,k-1)].bs->x[dim]);
                } else {
                  if (k<gl->domain_all.ke && NODEVALID[_ai_all(gl,i,j,k+1)]) {
                    dx3[dim]=(np[_ai(gl,i,j,k+1)].bs->x[dim]-np[_ai(gl,i,j,k)].bs->x[dim]);
                  } else {
                    if (k>gl->domain_all.ks && NODEVALID[_ai_all(gl,i,j,k-1)]) {
                      dx3[dim]=(np[_ai(gl,i,j,k)].bs->x[dim]-np[_ai(gl,i,j,k-1)].bs->x[dim]);
                    } else {
                      fatal_error("Couldn't find adjacent valid node along k needed for interpolation.");
                    }
                  }
                }
#endif
	      }
#ifdef DISTMPI
	    }
            if (rank!=0 && _node_rank(gl,i,j,k)==rank) MPI_Ssend(x,nd,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            if (rank==0 && _node_rank(gl,i,j,k)!=0) MPI_Recv(x,nd,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
            if (rank!=0 && _node_rank(gl,i,j,k)==rank) MPI_Ssend(dx1,nd,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            if (rank==0 && _node_rank(gl,i,j,k)!=0) MPI_Recv(dx1,nd,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
#ifdef _2DL  
            if (rank!=0 && _node_rank(gl,i,j,k)==rank) MPI_Ssend(dx2,nd,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            if (rank==0 && _node_rank(gl,i,j,k)!=0) MPI_Recv(dx2,nd,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
#endif
#ifdef _3DL  
            if (rank!=0 && _node_rank(gl,i,j,k)==rank) MPI_Ssend(dx3,nd,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            if (rank==0 && _node_rank(gl,i,j,k)!=0) MPI_Recv(dx3,nd,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
#endif

#endif
            if (pass==1) {
              wfwrite(initvar, sizeof(initvar_t), 1, datafile);
            } else {
#ifdef EMFIELD
              wfwrite(initvar_emfield, sizeof(initvar_emfield_t), 1, datafile);
#endif
            }
            wfwrite(x, sizeof(dim_t), 1, datafile);
            wfwrite(dx1, sizeof(dim_t), 1, datafile);
#ifdef _2DL
            wfwrite(dx2, sizeof(dim_t), 1, datafile);
#endif
#ifdef _3DL
            wfwrite(dx3, sizeof(dim_t), 1, datafile);
#endif
          } //end if nodevalid
    } // for_ijk
  }//pass

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  wfclose(datafile);
  wfprintf(stdout,"done.\n");

  free(NODEVALID);
  free(initvar_names);
}

