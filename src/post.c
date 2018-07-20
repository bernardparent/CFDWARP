#include <src/post.h>
#include <cycle/_cycle.h>
#include <src/bdry.h>


#define FSTAR_METRICS_INTERFACE 0
#define FSTAR_METRICS_EXTRAPOLATED 1
#define FSTAR_METRICS FSTAR_METRICS_EXTRAPOLATED


typedef char dimname_t[3];


static void set_dimname(dimname_t dimname, char dim1, char dim2, char dim3){
  dimname[0]=dim1;
  dimname[1]=dim2;
  dimname[2]=dim3;
}


void write_post_file_nodplot(np_t *np, gl_t *gl, zone_t zone, char *filename, bool GRIDONLY){
  long i,j,k,cnt;
  FILE *postfile;
  char varname[100];
  double val[numpostvar+nd];
#ifdef DISTMPI
  int rank;
  MPI_Status MPI_Status1;
#endif
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  postfile = wfopen(filename, "w");
  wfprintf(postfile, "%d\n",numpostvar+nd);
  for (cnt=0; cnt<numpostvar+nd; cnt++){
    find_post_proc_var_name(cnt, varname);
    wfprintf(postfile, "%s\n",varname);
  }

#ifdef _3D
  wfprintf(postfile, "%12ld%12ld%12ld%12ld%12ld%12ld\n",
           zone.is,zone.ie,zone.js,zone.je,zone.ks,zone.ke);
#endif
#ifdef _2D
  wfprintf(postfile, "%12ld%12ld%12ld%12ld%12d%12d\n",
           zone.is,zone.ie,zone.js,zone.je,1,1);
#endif

  wfprintf(stdout,"Writing to postfile %s for NODPLOT use..",filename);

  for3DL(k,zone.ks,zone.ke)
    for1DL(i,zone.is,zone.ie)
      for2DL(j,zone.js,zone.je) 
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          for (cnt=0; cnt<numpostvar+nd; cnt++){
            find_post_proc_var_value(np, _ai(gl,i,j,k), gl, cnt, &(val[cnt]));
          }
#ifdef DISTMPI
          if (rank!=0) MPI_Ssend(val,numpostvar,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        }
        if (rank==0) {
          if (_node_rank(gl,i,j,k)!=0)
            MPI_Recv(val,numpostvar,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
#endif

          for (cnt=0; cnt<numpostvar+nd; cnt++) wfprintf(postfile, "%.10E  ",val[cnt]);
          wfprintf(postfile, "\n");
#ifdef DISTMPI
        }
#endif
      end2DL
    end1DL
  end3DL
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  wfclose(postfile);
  wfprintf(stdout,"done.\n");
}



void write_post_file_gnuplot(np_t *np, gl_t *gl, zone_t zone, char *filename, bool GRIDONLY){
  long i,j,k,dim,cnt;
  FILE *postfile;
  dimname_t dimname;
  char varname[100];
  int numval,cntval;
  double val[numpostvar+nd+2];
#ifdef DISTMPI
  int rank;
  MPI_Status MPI_Status1;
#endif

  numval=nd;
  if (!GRIDONLY) numval+=numpostvar;

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  set_dimname(dimname,'X','Y','Z');
  postfile = wfopen(filename, "w");

  wfprintf(postfile, "# ");
  for (dim=0; dim<nd; dim++) wfprintf(postfile, " %c,",dimname[dim]);
  if (!GRIDONLY) {
    for (cnt=0; cnt<numpostvar; cnt++){
      find_post_proc_var_name(cnt, varname);
      wfprintf(postfile, " %s,",varname);
    }
  }
  wfprintf(postfile,"\n");

  wfprintf(stdout,"Writing to postfile %s for GNUPLOT use..",filename);

#ifdef DISTMPI
  if (rank==0 || gl->DISTDOMAIN) {
#endif
  for3DL(k,zone.ks,zone.ke)
    for2DL(j,zone.js,zone.je)
      for1DL(i,zone.is,zone.ie)
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          cntval=0;
          for (dim=0; dim<nd; dim++) {
            val[cntval]=np[_ai(gl,i,j,k)].bs->x[dim];
            cntval++;
          }
          for (cnt=0; cnt<numpostvar; cnt++){
            if (!GRIDONLY)
              find_post_proc_var_value(np, _ai(gl,i,j,k), gl, cnt, &(val[cntval]));
            else val[cntval]=0.0;
            cntval++;
          }
#ifdef DISTMPI
          if (rank!=0) MPI_Ssend(val,numval,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        }
        if (rank==0) {
          if (_node_rank(gl,i,j,k)!=0)
            MPI_Recv(val,numval,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
#endif
          for (cntval=0; cntval<numval; cntval++) wfprintf(postfile, "%.10E  ",val[cntval]);
          wfprintf(postfile, "\n");
#ifdef DISTMPI
        }
#endif
      end1DL
      if (zone.is!=zone.ie) {
        wfprintf(postfile, "\n");
      }
    end2DL
    if (zone.js!=zone.je && zone.is==zone.ie) {
      wfprintf(postfile, "\n");
    }
  end3DL
#ifdef DISTMPI
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  wfclose(postfile);
  wfprintf(stdout,"done.\n");
}



void write_post_file_tecplot(np_t *np, gl_t *gl, zone_t zone, char *filename, bool GRIDONLY, bool NOHEADER){
  long i,j,k,dim,cnt;
  FILE *postfile;
  dimname_t dimname;
#ifdef _2DL
  dimname_t axisname;
#endif
  char varname[100];
  int numval,cntval;
  double val[numpostvar+nd+2];
#ifdef DISTMPI
  int rank;
  MPI_Status MPI_Status1;
#endif

  numval=nd;
  if (!GRIDONLY) numval+=numpostvar;
  if (NOHEADER) numval++;

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  set_dimname(dimname,'X','Y','Z');
  postfile = wfopen(filename, "w");
  if (!NOHEADER){
    wfprintf(postfile, "VARIABLES=");
    for (dim=0; dim<nd; dim++) wfprintf(postfile, " \"%c\",",dimname[dim]);
    wfprintf(postfile, "\n");
    if (!GRIDONLY) {
      for (cnt=0; cnt<numpostvar; cnt++){
        find_post_proc_var_name(cnt, varname);
        wfprintf(postfile, "\"%s\",\n",varname);
      }
    }
#ifdef _2D
    set_dimname(axisname,'I','J','K');
    /* tecplot doesn't like I=1 in 2D. Only
       J=1 is accepted. Go figure. */
    if (zone.ie-zone.is==0) set_dimname(axisname,'J','I','K');
    wfprintf(postfile, "ZONE %c=%ld, %c=%ld F=POINT\n", axisname[0],zone.ie-zone.is+1,axisname[1],zone.je-zone.js+1);
#endif
#ifdef _3D
    set_dimname(axisname,'I','J','K');
    /* tecplot doesn't like I=1 or J=1 in 3D. Only
       K=1 is accepted.  */
    if (zone.ie-zone.is==0) set_dimname(axisname,'K','I','J');
    if (zone.je-zone.js==0) set_dimname(axisname,'I','K','J');
    wfprintf(postfile, "ZONE %c=%ld, %c=%ld, %c=%ld F=POINT\n",
             axisname[0],zone.ie-zone.is+1,axisname[1],zone.je-zone.js+1,axisname[2],
             zone.ke-zone.ks+1);
#endif
  }
  wfprintf(stdout,"Writing to postfile %s for TECPLOT use..",filename);

#ifdef DISTMPI
  if (rank==0 || gl->DISTDOMAIN) {
#endif
  for3DL(k,zone.ks,zone.ke)
    for2DL(j,zone.js,zone.je)
      for1DL(i,zone.is,zone.ie)
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          cntval=0;
          for (dim=0; dim<nd; dim++) {
            val[cntval]=np[_ai(gl,i,j,k)].bs->x[dim];
            cntval++;
          }
          if (NOHEADER){
            val[cntval]=0.0; /* time variable */
            cntval++;
          }
          for (cnt=0; cnt<numpostvar; cnt++){
            if (!GRIDONLY)
              find_post_proc_var_value(np, _ai(gl,i,j,k), gl, cnt, &(val[cntval]));
            else val[cntval]=0.0;
            cntval++;
          }
#ifdef DISTMPI
          if (rank!=0) MPI_Ssend(val,numval,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        }
        if (rank==0) {
          if (_node_rank(gl,i,j,k)!=0)
            MPI_Recv(val,numval,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
#endif
          for (cntval=0; cntval<numval; cntval++) wfprintf(postfile, "%.10E  ",val[cntval]);
          wfprintf(postfile, "\n");
#ifdef DISTMPI
        }
#endif
      end1DL
    end2DL
  end3DL
#ifdef DISTMPI
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  wfclose(postfile);
  wfprintf(stdout,"done.\n");
}


#ifdef _FLUID_MULTISPECIES
static double _w_post(np_t np, long spec) {
  if (spec>=ns) fatal_error("User specified species %ld is greater or equal to ns (%ld).",spec+1,ns);
  if (spec<0) fatal_error("User specified species %ld is smaller than 0.",spec+1);
  return(_w(np,spec));
}
#endif




bool _CELLVALID(np_t *np, gl_t *gl, long i, long j, long k){
  bool CELLVALID;
  CELLVALID=FALSE;
  if (   is_node_valid(np[_ai(gl,i+0,j+0,k+0)],TYPELEVEL_FLUID)
      && is_node_valid(np[_ai(gl,i+1,j+0,k+0)],TYPELEVEL_FLUID)
#ifdef _2DL
      && is_node_valid(np[_ai(gl,i+1,j+1,k+0)],TYPELEVEL_FLUID)
      && is_node_valid(np[_ai(gl,i+0,j+1,k+0)],TYPELEVEL_FLUID)
#endif
#ifdef _3DL
      && is_node_valid(np[_ai(gl,i+0,j+0,k+1)],TYPELEVEL_FLUID)
      && is_node_valid(np[_ai(gl,i+1,j+0,k+1)],TYPELEVEL_FLUID)
      && is_node_valid(np[_ai(gl,i+1,j+1,k+1)],TYPELEVEL_FLUID)
      && is_node_valid(np[_ai(gl,i+0,j+1,k+1)],TYPELEVEL_FLUID)
#endif
     ) CELLVALID=TRUE;

#ifdef EMFIELD
  if (   is_node_valid(np[_ai(gl,i+0,j+0,k+0)],TYPELEVEL_EMFIELD)
      && is_node_valid(np[_ai(gl,i+1,j+0,k+0)],TYPELEVEL_EMFIELD)
#ifdef _2DL
      && is_node_valid(np[_ai(gl,i+1,j+1,k+0)],TYPELEVEL_EMFIELD)
      && is_node_valid(np[_ai(gl,i+0,j+1,k+0)],TYPELEVEL_EMFIELD)
#endif
#ifdef _3DL
      && is_node_valid(np[_ai(gl,i+0,j+0,k+1)],TYPELEVEL_EMFIELD)
      && is_node_valid(np[_ai(gl,i+1,j+0,k+1)],TYPELEVEL_EMFIELD)
      && is_node_valid(np[_ai(gl,i+1,j+1,k+1)],TYPELEVEL_EMFIELD)
      && is_node_valid(np[_ai(gl,i+0,j+1,k+1)],TYPELEVEL_EMFIELD)
#endif
     ) CELLVALID=TRUE;
#endif
  return(CELLVALID);
}



void write_post_file_vtk(np_t *np, gl_t *gl, zone_t zone, char *filename, bool GRIDONLY){
  long ii,cnt,i,j,k,dim,postvar;
  FILE *postfile;
  int cntval;
  char varname[100],varname_tmp[100];
  long numcells,cntpoint,numpoints,varnd;
  long *l;
  char *strmatch;
  bool VECTORVAR;
  int *CELLVALID;
#ifdef _2D
  EXM_gl2D_t gl2d;
#endif
#ifdef _3D
  EXM_gl3D_t gl3d;
#endif
  double val[numpostvar+2*nd+2];
#ifdef DISTMPI
  int rank;
  MPI_Status MPI_Status1;
#endif

#ifdef _2D
  gl2d.is=zone.is;
  gl2d.ie=zone.ie;
  gl2d.js=zone.js;
  gl2d.je=zone.je;
  ii=(zone.ie-zone.is)*(zone.je-zone.js+1)+(zone.je-zone.js)+1;
#endif

#ifdef _3D
  gl3d.is=zone.is;
  gl3d.ie=zone.ie;
  gl3d.js=zone.js;
  gl3d.je=zone.je;
  gl3d.ks=zone.ks;
  gl3d.ke=zone.ke;
  ii=(zone.ie-zone.is)*(zone.je-zone.js+1)+(zone.je-zone.js);
  ii=ii*(zone.ke-zone.ks+1)+(zone.ke-zone.ks)+1;
#endif

  l=(long *)malloc(sizeof(long)*(ii));
  CELLVALID=(int *)malloc(sizeof(int)*ii);
  numpoints=(zone.ie-zone.is+1)*(zone.je-zone.js+1)if3D(*(zone.ke-zone.ks+1));

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  postfile = wfopen(filename, "w");
  wfprintf(postfile, "# vtk DataFile Version 2.0\n");
  wfprintf(postfile, "Data generated by CFDWARP\n");
  wfprintf(postfile, "ASCII\n");
  wfprintf(postfile, "\n");
  wfprintf(postfile, "DATASET UNSTRUCTURED_GRID\n");
  wfprintf(postfile, "POINTS %ld float\n", numpoints);
  wfprintf(stdout,"Writing to VTK postfile %s for PARAVIEW use..",filename);


/* here, output the grid x,y,z to the postfile */
  cntpoint=0;
  numcells=0;
  for3DL(k,zone.ks,zone.ke)
    for2DL(j,zone.js,zone.je)
      for1DL(i,zone.is,zone.ie)
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          cntval=0;
          for (dim=0; dim<nd; dim++) {
            val[cntval]=np[_ai(gl,i,j,k)].bs->x[dim];
            cntval++;
          }
          CELLVALID[if2D(EXM_ai2(gl2d,i,j))if3D(EXM_ai3(gl3d,i,j,k))]=_CELLVALID(np,gl,i,j,k);
#ifdef DISTMPI
          if (rank!=0) MPI_Ssend(val,nd,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
          if (rank!=0) MPI_Ssend(&(CELLVALID[if2D(EXM_ai2(gl2d,i,j))if3D(EXM_ai3(gl3d,i,j,k))]),1,MPI_INT,0,0,MPI_COMM_WORLD);
        }
        if (rank==0) {
          if (_node_rank(gl,i,j,k)!=0){
            MPI_Recv(val,nd,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
            MPI_Recv(&(CELLVALID[if2D(EXM_ai2(gl2d,i,j))if3D(EXM_ai3(gl3d,i,j,k))]),1,MPI_INT,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
          }
#endif
          for (cntval=0; cntval<nd; cntval++) wfprintf(postfile, "%.10E  ",val[cntval]);
          for (cntval=nd; cntval<3; cntval++) wfprintf(postfile, "0.0  ");
          wfprintf(postfile, "\n");
          if (CELLVALID[if2D(EXM_ai2(gl2d,i,j))if3D(EXM_ai3(gl3d,i,j,k))]) numcells++;
          l[if2D(EXM_ai2(gl2d,i,j))if3D(EXM_ai3(gl3d,i,j,k))]=cntpoint;
          cntpoint++;
#ifdef DISTMPI
        }
#endif
      end1DL
    end2DL
  end3DL


  cntpoint=0;
/* here, setup the cells */
  wfprintf(postfile, "\n");
  wfprintf(postfile, "CELLS %ld ",numcells);
  wfprintf(postfile, "%ld\n",if2D(5)if3D(9)*numcells);
  
  for3DL(k,zone.ks,zone.ke-1)
    for2DL(j,zone.js,zone.je-1)
      for1DL(i,zone.is,zone.ie-1)
         if (CELLVALID[if2D(EXM_ai2(gl2d,i,j))if3D(EXM_ai3(gl3d,i,j,k))]) {
#ifdef _2D
          wfprintf(postfile, "4 %ld %ld %ld %ld \n",l[EXM_ai2(gl2d,i,j)],l[EXM_ai2(gl2d,i+1,j)],l[EXM_ai2(gl2d,i+1,j+1)],l[EXM_ai2(gl2d,i,j+1)]);
#endif          
#ifdef _3D
          wfprintf(postfile, "8 %ld %ld %ld %ld %ld %ld %ld %ld \n",l[EXM_ai3(gl3d,i,j,k+1)],l[EXM_ai3(gl3d,i+1,j,k+1)],l[EXM_ai3(gl3d,i+1,j,k)],l[EXM_ai3(gl3d,i,j,k)],l[EXM_ai3(gl3d,i,j+1,k+1)],l[EXM_ai3(gl3d,i+1,j+1,k+1)],l[EXM_ai3(gl3d,i+1,j+1,k)],l[EXM_ai3(gl3d,i,j+1,k)]);
#endif          
         }
      end1DL
    end2DL
  end3DL
  wfprintf(postfile, "\n");
  wfprintf(postfile, "CELL_TYPES %ld\n",numcells);
  for (cnt=0; cnt<numcells; cnt++){
    wfprintf(postfile, if2D("9\n")if3D("12\n"));
  }


  if (!GRIDONLY){
    /* now, output the flow properties to the postfile */
    wfprintf(postfile, "\n");
    wfprintf(postfile, "POINT_DATA %ld\n",numpoints);

    postvar=0;
    do{

      find_post_proc_var_name(postvar, varname);
      if (strstr(varname,"[0]")==NULL) {
        VECTORVAR=FALSE;
        wfprintf(postfile, "SCALARS %s float \n",varname);
        varnd=1;
      } else {
        VECTORVAR=TRUE;
        strmatch=strstr(varname,"[0]");
        strmatch[0]='\0'; /* end of string character to get rid of the [0] within the varname */
        wfprintf(postfile, "VECTORS %s float \n",varname);
        varnd=1;
        if (postvar+1<numpostvar+nd) {
          find_post_proc_var_name(postvar+1, varname_tmp);
          if (strstr(varname_tmp,"[1]")!=NULL) {
            varnd=2;
            if (postvar+2<numpostvar+nd) {
              find_post_proc_var_name(postvar+2, varname_tmp);
              if (strstr(varname_tmp,"[2]")!=NULL) varnd=3;
            }
          }
        }
      }
      if (!VECTORVAR) wfprintf(postfile, "LOOKUP_TABLE default \n");
        for3DL(k,zone.ks,zone.ke)
          for2DL(j,zone.js,zone.je)
            for1DL(i,zone.is,zone.ie)
#ifdef DISTMPI
              if (_node_rank(gl,i,j,k)==rank) {
#endif

                for (dim=0; dim<varnd; dim++)
                  find_post_proc_var_value(np, _ai(gl,i,j,k), gl, postvar+dim, &(val[dim]));

#ifdef DISTMPI
                if (rank!=0) MPI_Ssend(val,varnd,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
              }
              if (rank==0) {
                if (_node_rank(gl,i,j,k)!=0)
                  MPI_Recv(val,varnd,MPI_DOUBLE,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
#endif
                if (VECTORVAR) {
                  for (dim=0; dim<varnd; dim++) wfprintf(postfile, "%.10E  ",val[dim]);
                  for (dim=varnd; dim<3; dim++) wfprintf(postfile, "0.0  ");
                } else {
                  wfprintf(postfile, "%.10E  ",val[0]);
                }
                wfprintf(postfile, "\n");
#ifdef DISTMPI
              }
#endif
            end1DL
          end2DL
        end3DL

      postvar=postvar+varnd;
    } while (postvar<numpostvar+nd);
  }
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  wfclose(postfile);
  wfprintf(stdout,"done.\n");
  free(l);
  free(CELLVALID);
}




void write_post_file(np_t *np, gl_t *gl, zone_t zone, char *filename, char *postprocessor, bool GRIDONLY){
  bool FOUND;

    zone.is=min(gl->domain_all.ie,max(gl->domain_all.is,zone.is));
    zone.js=min(gl->domain_all.je,max(gl->domain_all.js,zone.js));
    zone.ks=min(gl->domain_all.ke,max(gl->domain_all.ks,zone.ks));
    zone.ie=max(gl->domain_all.is,min(gl->domain_all.ie,zone.ie));
    zone.je=max(gl->domain_all.js,min(gl->domain_all.je,zone.je));
    zone.ke=max(gl->domain_all.ks,min(gl->domain_all.ke,zone.ke));
    FOUND=FALSE;
    if (strcmp(postprocessor,"gnuplot")==0) {
      write_post_file_gnuplot(np, gl, zone, filename, GRIDONLY);
      FOUND=TRUE;
    }
    if (strcmp(postprocessor,"vtk")==0) {
      write_post_file_vtk(np, gl, zone, filename, GRIDONLY);
      FOUND=TRUE;
    }
    if (strcmp(postprocessor,"nodplot")==0) {
      write_post_file_nodplot(np, gl, zone, filename, GRIDONLY);
      FOUND=TRUE;
    }
    if (strcmp(postprocessor,"tecplot")==0) {
      write_post_file_tecplot(np, gl, zone, filename, GRIDONLY, FALSE);
      FOUND=TRUE;
    }
    if (strcmp(postprocessor,"tecplotnoheader")==0) {
      write_post_file_tecplot(np, gl, zone, filename, GRIDONLY, TRUE);
      FOUND=TRUE;
    }
    
    if (!FOUND) write_post_file_tecplot(np, gl, zone, filename, GRIDONLY, FALSE);

}


static void interpolate_between_nodes(np_t *np, long lA, long lB, gl_t *gl, double factA, double factB, long *initnum, initvar_t initvar){
  long cnt;
  initvar_t initvarA,initvarB;
  find_default_initvar(np, gl, lA, initvarA);
  find_default_initvar(np, gl, lB, initvarB);
  for (cnt=0; cnt<numinitvar; cnt++){
    initvar[cnt]=factA*initvarA[cnt]+factB*initvarB[cnt];
  }
  *initnum=defaultinitvartypefluid;
}


void create_domain_of_cuts_along_x(np_t *np, gl_t *gl, np_t **npcut, gl_t *glcut, double *xcut, long numcut){
  long icut,i,j,k,dim;
  double fact,x;
  zone_t zone;
  dim_t xl,xc;
  int FOUND;
  initvar_t initvar;
  long initnum;
  long nodetype;

#ifdef DISTMPI
  int rank;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank!=0) fatal_error("Can't create x-cut when number of processes is more than 1."); 
#endif

  /* need to initialize xl and xc to avoid compiler warning */
  nodetype=0;
  for (dim=0; dim<nd; dim++) {
    xl[dim]=0.0;
    xc[dim]=0.0;
  }

  wfprintf(stdout,"Creating x-station..");
  if (numcut<1) fatal_error("Number of x-stations must be greater or equal to 1.");
  *glcut=*gl;
  glcut->domain=gl->domain_all;
  glcut->domain.is=1;
  glcut->domain.ie=numcut;

  init_data_structure_and_create_nodes(npcut, glcut, glcut->domain, glcut->domain);
  update_node_type(*npcut,glcut,TYPELEVEL_FLUID,NODETYPE_UNUSED,glcut->domain);

  /* first, find x,y,z on the new domain */
  for1DL(icut,glcut->domain.is,glcut->domain.ie)
    x=xcut[icut-1];
    wfprintf(stdout,"%E..",x);
    for2DL(j,gl->domain_all.js-1,gl->domain_all.je+1)
      for3DL(k,gl->domain_all.ks-1,gl->domain_all.ke+1)
        FOUND=FALSE;
        for1DL(i,gl->domain.is+1,gl->domain.ie)
          if (_x(np[_ai(gl,i+0,j,k)],0)>=x && _x(np[_ai(gl,i-1,j,k)],0)<x) {
            FOUND=TRUE;
            for (dim=0; dim<nd; dim++) {
              xl[dim]=_x(np[_ai(gl,i-1,j,k)],dim);
              xc[dim]=_x(np[_ai(gl,i+0,j,k)],dim);
            }
          }
        end1DL
        if (FOUND) {
          (*npcut)[_ai(glcut,icut,j,k)].bs->x[0]=x;
          fact=(x-xl[0])/(xc[0]-xl[0]);
          for (dim=1; dim<nd; dim++)
            (*npcut)[_ai(glcut,icut,j,k)].bs->x[dim]=(1.0-fact)*xl[dim]+fact*xc[dim];
        } else {
          fatal_error("Could not interpolate at j=%ld k=%ld for xcut=%E m..",j,k,x);
        }
      end3DL
    end2DL
  end1DL

  /* now, for the special case of a domain with only one node along i, find x,y,z on i-1 and i+1 */
  for2DL(j,glcut->domain.js-1,glcut->domain.je+1)
    for3DL(k,glcut->domain.ks-1,glcut->domain.ke+1)
      for (dim=0; dim<nd; dim++){
        (*npcut)[_ai(glcut,glcut->domain.is-1,j,k)].bs->x[dim]=
           (*npcut)[_ai(glcut,glcut->domain.is,j,k)].bs->x[dim];
        (*npcut)[_ai(glcut,glcut->domain.ie+1,j,k)].bs->x[dim]=
           (*npcut)[_ai(glcut,glcut->domain.ie,j,k)].bs->x[dim];
      }
      (*npcut)[_ai(glcut,glcut->domain.ie+1,j,k)].bs->x[0]+=1.0;
      (*npcut)[_ai(glcut,glcut->domain.is-1,j,k)].bs->x[0]-=1.0;
    end3DL
  end2DL

  zone=gl->domain;
  /* now, find U on the new domain */
  for1DL(icut,glcut->domain.is,glcut->domain.ie)
    x=xcut[icut-1];
    /* find the lower and high i in *np that englobes x */
    zone.ie=gl->domain.is;
    zone.is=gl->domain.ie;
    for2DL(j,gl->domain.js,gl->domain.je)
      for3DL(k,gl->domain.ks,gl->domain.ke)
        for1DL(i,gl->domain.is,gl->domain.ie)
          if (_x(np[_ai(gl,i,j,k)],0)>=x) zone.is=min(zone.is,i);
          if (_x(np[_ai(gl,i,j,k)],0)<=x) zone.ie=max(zone.ie,i);
        end1DL
      end3DL
    end2DL
    zone.is=max(gl->domain.is,zone.is-1);
    zone.ie=min(gl->domain.ie,zone.ie+1);
    //resume_nodes_only_in_zone_and_update_bdry_nodes(np,gl,zone);

    for2DL(j,gl->domain_all.js,gl->domain_all.je)
      for3DL(k,gl->domain_all.ks,gl->domain_all.ke)
        FOUND=FALSE;
        for1DL(i,gl->domain_all.is+1,gl->domain_all.ie)
          if (_x(np[_ai(gl,i+0,j,k)],0)>=x && _x(np[_ai(gl,i-1,j,k)],0)<x) {
            FOUND=TRUE;
            nodetype=min(_node_type(np[_ai(gl,i+0,j,k)],TYPELEVEL_FLUID),
                         _node_type(np[_ai(gl,i-1,j,k)],TYPELEVEL_FLUID));
            if (nodetype>=NODETYPE_INNER) {
              fact=(x-_x(np[_ai(gl,i-1,j,k)],0))/(_x(np[_ai(gl,i+0,j,k)],0)-_x(np[_ai(gl,i-1,j,k)],0));
              interpolate_between_nodes(np,_ai(gl,i,j,k),_ai(gl,i-1,j,k),gl,
                                      fact,1.0e0-fact,&initnum,initvar);
            }
          }
        end1DL

        if (FOUND) {
          update_node_type((*npcut),glcut,TYPELEVEL_FLUID,nodetype, _zone_from_point(icut,j,k));
          if (nodetype>=NODETYPE_INNER)
            init_node_fluid((*npcut),_ai(glcut,icut,j,k),glcut,initnum,initvar);
        } else {
          fatal_error("Could not interpolate at j=%ld k=%ld for xcut=%E m..",j,k,x);
        }
      end3DL
    end2DL
  end1DL
  wfprintf(stdout,"done.\n");
}



#ifdef _2DL
void create_domain_of_cuts_along_y(np_t *np, gl_t *gl, np_t **npcut, gl_t *glcut, double *ycut, long numcut){
  long jcut,i,j,k,dim;
  double fact,y;
  zone_t zone;
  dim_t xl,xc;
  int FOUND;
  initvar_t initvar;
  long initnum;
  long nodetype;

#ifdef DISTMPI
  int rank;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank!=0) fatal_error("Can't create y-cut when number of processes is more than 1."); 
#endif

  /* need to initialize xl and xc to avoid compiler warning */
  nodetype=0;
  for (dim=0; dim<nd; dim++) {
    xl[dim]=0.0;
    xc[dim]=0.0;
  }

  wfprintf(stdout,"Creating y-station..");
  if (numcut<1) fatal_error("Number of y-stations must be greater or equal to 1.");
  *glcut=*gl;
  glcut->domain=gl->domain_all;
  glcut->domain.js=1;
  glcut->domain.je=numcut;

  init_data_structure_and_create_nodes(npcut, glcut, glcut->domain, glcut->domain);
  update_node_type(*npcut,glcut,TYPELEVEL_FLUID,NODETYPE_UNUSED,glcut->domain);

  /* first, find x,y,z on the new domain */
  for2DL(jcut,glcut->domain.js,glcut->domain.je)
    y=ycut[jcut-1];
    wfprintf(stdout,"%E..",y);
    for1DL(i,gl->domain_all.is-1,gl->domain_all.ie+1)
      for3DL(k,gl->domain_all.ks-1,gl->domain_all.ke+1)
        FOUND=FALSE;
        for2DL(j,gl->domain.js+1,gl->domain.je)
          if (_x(np[_ai(gl,i,j+0,k)],1)>=y && _x(np[_ai(gl,i,j-1,k)],1)<y) {
            FOUND=TRUE;
            for (dim=0; dim<nd; dim++) {
              xl[dim]=_x(np[_ai(gl,i,j-1,k)],dim);
              xc[dim]=_x(np[_ai(gl,i,j+0,k)],dim);
            }
          }
        end1DL
        if (FOUND) {
          (*npcut)[_ai(glcut,i,jcut,k)].bs->x[1]=y;
          fact=(y-xl[1])/(xc[1]-xl[1]);
          for (dim=0; dim<nd; dim++) {
            if (dim!=1) (*npcut)[_ai(glcut,i,jcut,k)].bs->x[dim]=(1.0-fact)*xl[dim]+fact*xc[dim];
          }
        } else {
          fatal_error("Could not interpolate at i=%ld k=%ld for ycut=%E m..",i,k,y);
        }
      end2DL
    end1DL
  end3DL

  /* now, for the special case of a domain with only one node along j, find x,y,z on j-1 and j+1 */
  for1DL(i,glcut->domain.is-1,glcut->domain.ie+1)
    for3DL(k,glcut->domain.ks-1,glcut->domain.ke+1)
      for (dim=0; dim<nd; dim++){
        (*npcut)[_ai(glcut,i,glcut->domain.js-1,k)].bs->x[dim]=
           (*npcut)[_ai(glcut,i,glcut->domain.js,k)].bs->x[dim];
        (*npcut)[_ai(glcut,i,glcut->domain.je+1,k)].bs->x[dim]=
           (*npcut)[_ai(glcut,i,glcut->domain.je,k)].bs->x[dim];
      }
      (*npcut)[_ai(glcut,i,glcut->domain.je+1,k)].bs->x[1]+=1.0;
      (*npcut)[_ai(glcut,i,glcut->domain.js-1,k)].bs->x[1]-=1.0;
    end3DL
  end1DL

  zone=gl->domain;
  /* now, find U on the new domain */
  for2DL(jcut,glcut->domain.js,glcut->domain.je)
    y=ycut[jcut-1];
    /* find the lower and high j in *np that englobes y */
    zone.je=gl->domain.js;
    zone.js=gl->domain.je;
    for1DL(i,gl->domain.is,gl->domain.ie)
      for2DL(j,gl->domain.js,gl->domain.je)
        for3DL(k,gl->domain.ks,gl->domain.ke)
          if (_x(np[_ai(gl,i,j,k)],1)>=y) zone.js=min(zone.js,j);
          if (_x(np[_ai(gl,i,j,k)],1)<=y) zone.je=max(zone.je,j);
        end3DL
      end2DL
    end1DL
    zone.js=max(gl->domain.js,zone.js-1);
    zone.je=min(gl->domain.je,zone.je+1);
    //resume_nodes_only_in_zone_and_update_bdry_nodes(np,gl,zone);

    for1DL(i,gl->domain_all.is,gl->domain_all.ie)
      for3DL(k,gl->domain_all.ks,gl->domain_all.ke)
        FOUND=FALSE;
        for2DL(j,gl->domain_all.js+1,gl->domain_all.je)
          if (_x(np[_ai(gl,i,j+0,k)],1)>=y && _x(np[_ai(gl,i,j-1,k)],1)<y) {
            FOUND=TRUE;
            nodetype=min(_node_type(np[_ai(gl,i,j+0,k)],TYPELEVEL_FLUID),
                         _node_type(np[_ai(gl,i,j-1,k)],TYPELEVEL_FLUID));
            if (nodetype>=NODETYPE_INNER) {
              fact=(y-_x(np[_ai(gl,i,j-1,k)],1))/(_x(np[_ai(gl,i,j+0,k)],1)-_x(np[_ai(gl,i,j-1,k)],1));
              interpolate_between_nodes(np,_ai(gl,i,j,k),_ai(gl,i,j-1,k),gl,
                                        fact,1.0e0-fact,&initnum,initvar);
            }
          }
        end2DL

        if (FOUND) {
          update_node_type((*npcut),glcut,TYPELEVEL_FLUID,nodetype, _zone_from_point(i,jcut,k));
          if (nodetype>=NODETYPE_INNER)
            init_node_fluid((*npcut),_ai(glcut,i,jcut,k),glcut,initnum,initvar);
        } else {
          fatal_error("Could not interpolate at i=%ld k=%ld for ycut=%E m..",i,k,y);
        }
      end3DL
    end1DL
  end2DL
  wfprintf(stdout,"done.\n");
}
#endif



#ifdef _3D
void create_domain_of_cuts_along_z(np_t *np, gl_t *gl, np_t **npcut, gl_t *glcut, double *zcut, long numcut){
  long kcut,i,j,k,dim;
  double fact,z;
  zone_t zone;
  dim_t xl,xc;
  int FOUND;
  initvar_t initvar;
  long initnum;
  long nodetype;

#ifdef DISTMPI
  int rank;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank!=0) fatal_error("Can't create z-cut when number of processes is more than 1."); 
#endif

  /* need to initialize xl and xc to avoid compiler warning */
  nodetype=0;
  for (dim=0; dim<nd; dim++) {
    xl[dim]=0.0;
    xc[dim]=0.0;
  }

  wfprintf(stdout,"Creating z-station..");
  if (numcut<1) fatal_error("Number of z-stations must be greater or equal to 1.");
  *glcut=*gl;
  glcut->domain=gl->domain_all;
  glcut->domain.ks=1;
  glcut->domain.ke=numcut;

  init_data_structure_and_create_nodes(npcut, glcut, glcut->domain, glcut->domain);
  update_node_type(*npcut,glcut,TYPELEVEL_FLUID,NODETYPE_UNUSED,glcut->domain);

  /* first, find x,y,z on the new domain */
  for3DL(kcut,glcut->domain.ks,glcut->domain.ke)
    z=zcut[kcut-1];
    wfprintf(stdout,"%E..",z);
    for1DL(i,gl->domain_all.is-1,gl->domain_all.ie+1)
      for2DL(j,gl->domain_all.js-1,gl->domain_all.je+1)
        FOUND=FALSE;
        for3DL(k,gl->domain.ks+1,gl->domain.ke)
          if (_x(np[_ai(gl,i,j,k+0)],2)>=z && _x(np[_ai(gl,i,j,k-1)],2)<z) {
            FOUND=TRUE;
            for (dim=0; dim<nd; dim++) {
              xl[dim]=_x(np[_ai(gl,i,j,k-1)],dim);
              xc[dim]=_x(np[_ai(gl,i,j,k+0)],dim);
            }
          }
        end1DL
        if (FOUND) {
          (*npcut)[_ai(glcut,i,j,kcut)].bs->x[2]=z;
          fact=(z-xl[2])/(xc[2]-xl[2]);
          for (dim=0; dim<nd-1; dim++)
            (*npcut)[_ai(glcut,i,j,kcut)].bs->x[dim]=(1.0-fact)*xl[dim]+fact*xc[dim];
        } else {
          fatal_error("Could not interpolate at i=%ld j=%ld for zcut=%E m..",i,j,z);
        }
      end2DL
    end1DL
  end3DL

  /* now, for the special case of a domain with only one node along k, find x,y,z on k-1 and k+1 */
  for1DL(i,glcut->domain.is-1,glcut->domain.ie+1)
    for2DL(j,glcut->domain.js-1,glcut->domain.je+1)
      for (dim=0; dim<nd; dim++){
        (*npcut)[_ai(glcut,i,j,glcut->domain.ks-1)].bs->x[dim]=
           (*npcut)[_ai(glcut,i,j,glcut->domain.ks)].bs->x[dim];
        (*npcut)[_ai(glcut,i,j,glcut->domain.ke+1)].bs->x[dim]=
           (*npcut)[_ai(glcut,i,j,glcut->domain.ke)].bs->x[dim];
      }
      (*npcut)[_ai(glcut,i,j,glcut->domain.ke+1)].bs->x[2]+=1.0;
      (*npcut)[_ai(glcut,i,j,glcut->domain.ks-1)].bs->x[2]-=1.0;
    end2DL
  end1DL

  zone=gl->domain;
  /* now, find U on the new domain */
  for3DL(kcut,glcut->domain.ks,glcut->domain.ke)
    z=zcut[kcut-1];
    /* find the lower and high k in *np that englobes z */
    zone.ke=gl->domain.ks;
    zone.ks=gl->domain.ke;
    for1DL(i,gl->domain.is,gl->domain.ie)
      for2DL(j,gl->domain.js,gl->domain.je)
        for3DL(k,gl->domain.ks,gl->domain.ke)
          if (_x(np[_ai(gl,i,j,k)],2)>=z) zone.ks=min(zone.ks,k);
          if (_x(np[_ai(gl,i,j,k)],2)<=z) zone.ke=max(zone.ke,k);
        end3DL
      end2DL
    end1DL
    zone.ks=max(gl->domain.ks,zone.ks-1);
    zone.ke=min(gl->domain.ke,zone.ke+1);
    //resume_nodes_only_in_zone_and_update_bdry_nodes(np,gl,zone);

    for1DL(i,gl->domain_all.is,gl->domain_all.ie)
      for2DL(j,gl->domain_all.js,gl->domain_all.je)
        FOUND=FALSE;
        for3DL(k,gl->domain_all.ks+1,gl->domain_all.ke)
          if (_x(np[_ai(gl,i,j,k+0)],2)>=z && _x(np[_ai(gl,i,j,k-1)],2)<z) {
            FOUND=TRUE;
            nodetype=min(_node_type(np[_ai(gl,i,j,k+0)],TYPELEVEL_FLUID),
                         _node_type(np[_ai(gl,i,j,k-1)],TYPELEVEL_FLUID));
            if (nodetype>=NODETYPE_INNER) {
              fact=(z-_x(np[_ai(gl,i,j,k-1)],2))/(_x(np[_ai(gl,i,j,k+0)],2)-_x(np[_ai(gl,i,j,k-1)],2));
              interpolate_between_nodes(np,_ai(gl,i,j,k),_ai(gl,i,j,k-1),gl,
                                        fact,1.0e0-fact,&initnum,initvar);
            }
          }
        end1DL

        if (FOUND) {
          update_node_type((*npcut),glcut,TYPELEVEL_FLUID,nodetype, _zone_from_point(i,j,kcut));
          if (nodetype>=NODETYPE_INNER)
            init_node_fluid((*npcut),_ai(glcut,i,j,kcut),glcut,initnum,initvar);
        } else {
          fatal_error("Could not interpolate at i=%ld j=%ld for zcut=%E m..",i,j,z);
        }
      end2DL
    end1DL
  end3DL
  wfprintf(stdout,"done.\n");
}
#endif



#ifdef EMFIELD

static void integrate_emfield_force(np_t *np, gl_t *gl, zone_t zone, dim_t Femfield){
  long i,j,k;
  long l,dim;
  dim_t Femfield_per_volume;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

  for (dim=0; dim<nd; dim++) Femfield[dim]=0.0;

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==NODETYPE_INNER) {
            find_emfield_force(np, gl, l, Femfield_per_volume);
            for (dim=0; dim<nd; dim++) Femfield[dim]+=_Omega(np[l],gl)*Femfield_per_volume[dim];
	  }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL


  /* here sum up all the contributions from the different processes */
#ifdef DISTMPI
  for (dim=0; dim<nd; dim++) {
    MPI_Allreduce(&(Femfield[dim]), &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Femfield[dim]=tempsum;
  }
#endif
}


static void integrate_emfield_work(np_t *np, gl_t *gl, zone_t zone,
                                   double *Wemfield){
  long i,j,k;
  long l,dim;
  dim_t Femfield_per_volume;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  *Wemfield=0.0;

  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==NODETYPE_INNER) {
            find_emfield_force(np, gl, l,Femfield_per_volume);
            for (dim=0; dim<nd; dim++) {
	            *Wemfield+=_V(np[l], dim)*Femfield_per_volume[dim]*_Omega(np[l],gl);
            }
          }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL

  /* here sum up all the contributions from the different processes */
#ifdef DISTMPI
  MPI_Allreduce(Wemfield, &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *Wemfield=tempsum;
#endif
}



static void integrate_EdotJ(np_t *np, gl_t *gl, zone_t zone, double *EdotJ){
  long i,j,k;
  long l;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  *EdotJ=0.0;

  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==NODETYPE_INNER) {
	    *EdotJ+=_E_dot_J(np,gl,l)*_Omega(np[l],gl);
	  }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL

  /* here sum up all the contributions from the different processes */
#ifdef DISTMPI
  MPI_Allreduce(EdotJ, &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *EdotJ=tempsum;
#endif
}




static void integrate_Qbeam(np_t *np, gl_t *gl, zone_t zone, double *Qbeam){
  long i,j,k;
  long l;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  *Qbeam=0.0;

  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==NODETYPE_INNER) {
	    *Qbeam+=_Omega(np[l],gl)*_Qbeam(np[l],gl);
	  }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL

  /* here sum up all the contributions from the different processes */
#ifdef DISTMPI
  MPI_Allreduce(Qbeam, &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *Qbeam=tempsum;
#endif
}

#endif           






static void find_surface_area_given_metrics(np_t np, gl_t *gl, metrics_t metrics, long theta, dim_t A){
  long dim;
  for (dim=0; dim<nd; dim++)
    A[dim]=fabs(metrics.Omega*metrics.X2[theta][dim]);
}



static void integrate_area_on_bdry(np_t *np, gl_t *gl, zone_t zone,
                                      dim_t Awall, long BDRYTYPE){
  long i,j,k,dim;
  long l,theta,thetasgn;
  flux_t tmpp1h;
  metrics_t metrics;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  for (dim=0; dim<nd; dim++) Awall[dim]=0.0;

  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==BDRYTYPE) {
            if (find_bdry_direc(np, gl, l, TYPELEVEL_FLUID, &theta, &thetasgn)) {
	              if (thetasgn>0)
                    find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metrics);
                  else find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, &metrics);
                  find_surface_area_given_metrics(np[l], gl, metrics, theta, tmpp1h);
              for (dim=0; dim<nd; dim++) Awall[dim]+=tmpp1h[dim];
            }
          }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL
  /* here sum up all the contributions in Fwall_shear from the different processes */
#ifdef DISTMPI
  for (dim=0; dim<nd; dim++) {
    MPI_Allreduce(&(Awall[dim]), &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Awall[dim]=tempsum;
  }
#endif
}

#if (fluxmom >=0)
static void integrate_shear_force_on_bdry(np_t *np, gl_t *gl, zone_t zone,
                                             dim_t Fwall_shear, long BDRYTYPE){
  long i,j,k,dim;
  long theta,thetasgn;
  long l,vartheta,flux;
  flux_t Gp1,Gp0,dGp1h,tmpp1h;
  flux_t Gp0p1,Gp0m1,Gp1p1,Gp1m1;
  sqmat_t Kp1h;
  metrics_t metricsp1h;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  for (dim=0; dim<nd; dim++) Fwall_shear[dim]=0.0;
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==BDRYTYPE) {
            if (find_bdry_direc(np, gl, l, TYPELEVEL_FLUID, &theta, &thetasgn)) {
              /* find the shear forces at the wall */
              for (vartheta=0; vartheta<nd; vartheta++){
	        if (thetasgn>0)
                find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metricsp1h);
                  else find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, &metricsp1h);

                find_Kstar_interface(np,gl,_al(gl,l,theta,+0),_al(gl,l,theta,thetasgn),metricsp1h,theta,vartheta,Kp1h);

                if (theta==vartheta) {
                  find_G(np[_al(gl,l,theta,thetasgn)],gl,Gp1);
                  find_G(np[_al(gl,l,theta,+0)],gl,Gp0);
                  for (flux=0; flux<nf; flux++)
                    dGp1h[flux]=thetasgn*(Gp1[flux]-Gp0[flux]);
                } else {
                  find_G(np[_all(gl,l,theta,+0,vartheta,+1)],gl,Gp0p1);
                  find_G(np[_all(gl,l,theta,+0,vartheta,-1)],gl,Gp0m1);
                  find_G(np[_all(gl,l,theta,thetasgn,vartheta,+1)],gl,Gp1p1);
                  find_G(np[_all(gl,l,theta,thetasgn,vartheta,-1)],gl,Gp1m1);
                  for (flux=0; flux<nf; flux++)
                    dGp1h[flux]=0.25e0*(Gp0p1[flux]+Gp1p1[flux]-Gp0m1[flux]-Gp1m1[flux]);
                }
                multiply_matrix_and_vector(Kp1h,dGp1h,tmpp1h);
                for (dim=0; dim<nd; dim++) Fwall_shear[dim]+=sign(metricsp1h.Omega)*thetasgn*tmpp1h[fluxmom+dim];
              }
            }
          }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL
  /* here sum up all the contributions in Fwall_shear from the different processes */
#ifdef DISTMPI
  for (dim=0; dim<nd; dim++) {
    MPI_Allreduce(&(Fwall_shear[dim]), &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Fwall_shear[dim]=tempsum;
  }
#endif
}
#endif

#if (fluxet>=0)
static void integrate_heat_to_surface_on_bdry(np_t *np, gl_t *gl, zone_t zone,
                                              double *heat_to_surface, long BDRYTYPE){
  long i,j,k;
  long theta,thetasgn;
  long row,col,l,vartheta,flux;
  flux_t Gp1,Gp0,dGp1h,tmpp1h;
  flux_t Gp0p1,Gp0m1,Gp1p1,Gp1m1;
  sqmat_t Kp1h;
  metrics_t metricsp1h;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  *heat_to_surface=0.0;
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==BDRYTYPE) {
            if (find_bdry_direc(np, gl, l, TYPELEVEL_FLUID, &theta, &thetasgn)) {
              /* find the shear forces at the wall */
              for (vartheta=0; vartheta<nd; vartheta++){
	        if (thetasgn>0)
                find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metricsp1h);
                  else find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, &metricsp1h);

                find_Kstar_interface(np,gl,_al(gl,l,theta,+0),_al(gl,l,theta,thetasgn),metricsp1h,theta,vartheta,Kp1h);
                for (row=0; row<nf; row++){
                  for (col=0; col<nf; col++){
                    if (row!=col) Kp1h[row][col]=0.0;
                  }
                }
                if (theta==vartheta) {
                  find_G(np[_al(gl,l,theta,thetasgn)],gl,Gp1);
                  find_G(np[_al(gl,l,theta,+0)],gl,Gp0);
                  for (flux=0; flux<nf; flux++)
                    dGp1h[flux]=thetasgn*(Gp1[flux]-Gp0[flux]);
                } else {
                  find_G(np[_all(gl,l,theta,+0,vartheta,+1)],gl,Gp0p1);
                  find_G(np[_all(gl,l,theta,+0,vartheta,-1)],gl,Gp0m1);
                  find_G(np[_all(gl,l,theta,thetasgn,vartheta,+1)],gl,Gp1p1);
                  find_G(np[_all(gl,l,theta,thetasgn,vartheta,-1)],gl,Gp1m1);
                  for (flux=0; flux<nf; flux++)
                    dGp1h[flux]=0.25e0*(Gp0p1[flux]+Gp1p1[flux]-Gp0m1[flux]-Gp1m1[flux]);
                }
                multiply_matrix_and_vector(Kp1h,dGp1h,tmpp1h);
                (*heat_to_surface)+=sign(metricsp1h.Omega)*thetasgn*tmpp1h[fluxet];
              }
            }
          }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL
  /* here sum up all the contributions in heat_to_surface from the different processes */
#ifdef DISTMPI
  MPI_Allreduce(&(*heat_to_surface), &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *heat_to_surface=tempsum;
#endif
}
#endif


#if (fluxmom>=0)

static void integrate_pressure_force_on_bdry(np_t *np, gl_t *gl, zone_t zone,
                                                dim_t Fwall_P, long BDRYTYPE){
  long i,j,k,dim;
  long l,theta,thetasgn;
  flux_t tmpp1h;
  metrics_t metrics;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  for (dim=0; dim<nd; dim++) Fwall_P[dim]=0.0;

  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (_node_type(np[l], TYPELEVEL_FLUID)==BDRYTYPE) {
            if (find_bdry_direc(np, gl, l, TYPELEVEL_FLUID, &theta, &thetasgn)) {

              switch (FSTAR_METRICS){
                case FSTAR_METRICS_INTERFACE:
                  // find the metrics at the interface between the inner node and the bdry node 
	              if (thetasgn>0)
                    find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metrics);
                  else find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, &metrics);
                  find_Fstar_given_metrics(np[l], gl, metrics, theta, tmpp1h);
                break;
                case FSTAR_METRICS_EXTRAPOLATED:
                  // find the pressure forces at the wall 
                  find_Fstar(np[l], gl, theta, tmpp1h);
                  find_metrics_at_node(np, gl, l, theta, &metrics);
                break;
                default:
                  fatal_error("In integrate_pressure_force_on_bdry(), FSTAR_METRICS is invalid.");
              }
              for (dim=0; dim<nd; dim++) Fwall_P[dim]-=sign(metrics.Omega)*thetasgn*tmpp1h[fluxmom+dim];

            }
          }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL
  /* here sum up all the contributions in Fwall_shear from the different processes */
#ifdef DISTMPI
  for (dim=0; dim<nd; dim++) {
    MPI_Allreduce(&(Fwall_P[dim]), &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Fwall_P[dim]=tempsum;
  }
#endif
}



#endif

#if (fluxet>=0)
static void integrate_metotal(np_t *np, gl_t *gl, zone_t zone, double *metotal){
  long i,j,k;
  long l;
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  *metotal=0.0;
  
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (is_node_inner(np[l], TYPELEVEL_FLUID)) {
#ifdef fluxet
            *metotal+=np[l].bs->U[fluxet]*_Omega(np[l],gl);
#endif
          }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL
  /* here sum up all the contributions in Fwall_shear from the different processes */
#ifdef DISTMPI
  MPI_Allreduce(metotal, &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *metotal=tempsum;
#endif
}
#endif


static void integrate_mass(np_t *np, gl_t *gl, zone_t zone, double *mass){
  long i,j,k;
  long l;
#ifndef _FLUID_NUMBERDENSITIES
  long spec;
#endif
#ifdef DISTMPI
  int rank;
  double tempsum;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  *mass=0.0;
  
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        l=_ai(gl,i,j,k);
#ifdef DISTMPI
        if (_node_rank(gl,i,j,k)==rank) {
#endif
          if (is_node_inner(np[l], TYPELEVEL_FLUID)) {
#ifndef _FLUID_NUMBERDENSITIES
            for (spec=0; spec<ns; spec++) *mass+=np[l].bs->U[spec]*_Omega(np[l],gl);
#endif
          }
#ifdef DISTMPI
        }
#endif
      end3DL
    end2DL
  end1DL
  /* here sum up all the contributions in Fwall_shear from the different processes */
#ifdef DISTMPI
  MPI_Allreduce(mass, &tempsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *mass=tempsum;
#endif
}




void write_post_template(FILE **controlfile){
  wfprintf(*controlfile,
  "\n\n"
  "Post(\n"
  "  xstation[1]=0.1; {m}\n"
  "  xstation[2]=0.2; {m}\n"
  "  xstation[3]=0.3; {m}\n"
  "  numsteps=300;\n"
  "  qmin=1.0; {m/s, as small a value as possible}\n"
  "  Pback_min=400; {Pa}\n"
  "  Pback_max=40000; {Pa}\n"
  "  Aback=1.0; {m2"if2D("/m")"}\n"
  "  for (cnt,1,3,\n"
  "    SetXstation(xstation[cnt]);\n"
  "    Pback=_Pback_xstation(Aback, Pback_min, Pback_max, numsteps, qmin);\n"
  "    Fpot=_Fpot_xstation(Pback, numsteps, qmin);\n"
  "    mdot=_mdot_xstation();\n"
  "    Tstag=_Tstag_xstation();\n"
  "    Pstag=_Pstag_xstation(numsteps);\n"
  "    Pstar=_Pstar_xstation();\n"
  "    T=_T_xstation();\n"
  "    q=_q_xstation();\n"
  "    rho=_rho_xstation();\n"
  "    htstar=_htstar_xstation();\n"
  "    printf(\"x      = %%E m\\n\"\n"
  "           \"Pback  = %%E Pa\\n\"\n"
  "           \"Fpot   = %%E Ns/kg\\n\"\n"
  "           \"mdot   = %%E kg/"if2D("m")"s\\n\"\n"
  "           \"htstar = %%E J/kg\\n\"\n"
  "           \"Tstag  = %%E K\\n\"\n"
  "           \"Pstag  = %%E Pa\\n\"\n"
  "           \"Pstar  = %%E Pa\\n\"\n"
  "           \"T      = %%E K\\n\"\n"
  "           \"q      = %%E m/s\\n\"\n"
  "           \"rho    = %%E kg/m3\\n\\n\"\n"
  "           ,xstation[cnt],Pback,Fpot,mdot,htstar,Tstag,Pstag,Pstar,T,q,rho);\n"
  "  );\n"
  "  printf(\"\\n\");\n"
  "  for (dim,1,%ld,\n"
  "    Area[dim]=_Area(is,js,"if3D("ks,")" ie,je,"if3D("ke,")" dim,BDRY_WALLTFIXED1);\n"
#if (fluxmom>=0)
  "    Fshear[dim]=_Fshear(is,js,"if3D("ks,")" ie,je,"if3D("ke,")" dim,BDRY_WALLTFIXED1);\n"
  "    Fpressure[dim]=_Fpressure(is,js,"if3D("ks,")" ie,je,"if3D("ke,")" dim,BDRY_WALLTFIXED1);\n"  
#endif
#ifdef EMFIELD
  "    Femfield[dim]=_Femfield(is,js,"if3D("ks,")" ie,je,"if3D("ke,")" dim);\n"  
#endif
  "  );\n"  
#if (fluxmom>=0)
  "  printf(\"Fshear    = (%%+E"if2DL(",%%+E")if3DL(",%%+E")") N"if2D("/m")"\\n\",Fshear[1]"if2DL(",Fshear[2]")if3DL(",Fshear[3]")");\n"
  "  printf(\"Fpressure = (%%+E"if2DL(",%%+E")if3DL(",%%+E")") N"if2D("/m")"\\n\",Fpressure[1]"if2DL(",Fpressure[2]")if3DL(",Fpressure[3]")");\n"
#endif
#if (fluxet>=0)
  "  printf(\"Qheat     = %%+E W"if2D("/m")"\\n\",_Qheat(is,js,"if3D("ks,")" ie,je,"if3D("ke,")" BDRY_WALLTFIXED1));\n"
  "  printf(\"metotal   = %%+E J"if2D("/m")"\\n\",_metotal(is,js,"if3D("ks,")" ie,je"if3D(",ke")"));\n"
#endif
  "  printf(\"m         = %%+E J"if2D("/m")"\\n\",_m(is,js,"if3D("ks,")" ie,je"if3D(",ke")"));\n"
#ifdef EMFIELD
  "  printf(\"Femfield  = (%%+E"if2DL(",%%+E")if3DL(",%%+E")") N"if2D("/m")"\\n\",Femfield[1]"if2DL(",Femfield[2]")if3DL(",Femfield[3]")");\n"
  "  printf(\"Qbeam     = %%+E W"if2D("/m")"\\n\",_Qbeam(is,js,"if3D("ks,")" ie,je"if3D(",ke")"));\n"
  "  printf(\"EdotJ     = %%+E W"if2D("/m")"\\n\",_EdotJ(is,js,"if3D("ks,")" ie,je"if3D(",ke")"));\n"
  "  printf(\"Wemfield  = %%+E W"if2D("/m")"\\n\",_Wemfield(is,js,"if3D("ks,")" ie,je"if3D(",ke")")); {Wemfield=Femfield dot Vn}\n"
#endif
  "  printf(\"\\n\");\n"
  "  {\n"
  "  POSTGRIDONLY=FALSE;\n"
  "  WritePostFile(is,js,"if3D("ks,")" ie,je,"if3D("ke,")" \"post.01\", \"tecplot\", POSTGRIDONLY);\n"
  "  }\n"
  ");\n",nd,nd);
}



static double area_yz(np_t *np, gl_t *gl, long i, long j, long k){
  double area;
  long js,je;
#ifdef _3D
  EXM_vec3D_t A[4];
  long cnt,dim,k2,j2,ks,ke;
#endif

#ifdef _3D
    ks=gl->domain.ks;
    js=gl->domain.js;
    je=gl->domain.je;
    ke=gl->domain.ke;
    for (cnt=0; cnt<4; cnt++){
      j2=j;
      k2=k;
      if (cnt<=1) j2=j2+1;
      if (mod(cnt,2)==0) k2=k2-1;
      for (dim=0; dim<3; dim++){
        A[cnt][dim]=0.25e0*(
           _x(np[_ai(gl,i,min(max(j2,js),je),min(max(k2,ks),ke))],dim)
          +_x(np[_ai(gl,i,min(max(j2,js),je),min(max(k2+1,ks),ke))],dim)
          +_x(np[_ai(gl,i,min(max(j2-1,js),je),min(max(k2,ks),ke))],dim)
          +_x(np[_ai(gl,i,min(max(j2-1,js),je),min(max(k2+1,ks),ke))],dim));
      }
    }
    area=EXM_area_quadrilateral(A[0], A[1], A[2], A[3]);
#endif
#ifdef _2D
    js=gl->domain.js;
    je=gl->domain.je;
    area=sqrt(sqr(_x(np[_ai(gl,i,min(j+1,je),k)],0)-_x(np[_ai(gl,i,max(j-1,js),k)],0))
           +sqr(_x(np[_ai(gl,i,min(j+1,je),k)],1)-_x(np[_ai(gl,i,max(j-1,js),k)],1)))/2.0e0;
#endif
  return(area);
}


typedef struct {
  gl_t *gl;
  np_t *np;
  long i,numsteps;
  double Aback,q_min;
  zone_t zone;
} argfP_t;


static double fPstarc(void *arg, double Pback){
  double A,qc,rhoc,q_min;
  long l,i,j,k,numsteps;
  gl_t *gl;
  np_t *np;
  zone_t zone;

  gl=((argfP_t *)arg)->gl;
  np=((argfP_t *)arg)->np;
  zone=((argfP_t *)arg)->zone;
  A=-((argfP_t *)arg)->Aback;
  numsteps=((argfP_t *)arg)->numsteps;
  q_min=((argfP_t *)arg)->q_min;
  i=((argfP_t *)arg)->i;
  /* integrate i station */
  for3DL(k,zone.ks,zone.ke)
    for2DL(j,zone.js,zone.je)
      l=_ai(gl,i,j,k);
      if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  ) {
        if (reversibly_expand_q_and_rho_to_given_P(np[l], gl, Pback, &qc, &rhoc, numsteps, q_min)==0)
          A+=_rho(np[_ai(gl,i,j,k)])*_V(np[_ai(gl,i,j,k)],0)*area_yz(np,gl,i,j,k)/rhoc/qc;
      }
    end2DL
  end3DL
  return(A);
}

static void verify_post_domain_is_single_x_plane(np_t *np, gl_t gl, zone_t zone,
                                                 SOAP_codex_t *codex, char *message){
  long j,k,l;
  double x;
  bool FOUNDONE;
  if (zone.ie!=zone.is)
    SOAP_fatal_error(codex,"The post domain sent to %s must have only one i plane.",message);
#ifdef DISTMPI
  if (gl.DISTDOMAIN)
    SOAP_fatal_error(codex,"In DISTMPI mode, the post domain sent to %s must be a single x-station.",message);
#endif
  FOUNDONE=FALSE;
  x=0.0;
  for3DL(k,zone.ks,zone.ke)
    for2DL(j,zone.js,zone.je)
      l=_ai(&gl,zone.is,j,k);
      if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  ) {
        if (FOUNDONE && fabs(_x(np[l],0)-x)>1.0e-10)
	  SOAP_fatal_error(codex,"The nodes in the i plane sent to %s do not share the same x.",message);
	FOUNDONE=TRUE;
	x=_x(np[l],0);
      }
    end2DL
  end3DL
  //resume_nodes_only_in_zone_and_update_bdry_nodes(np,&gl,gl.domain_all);
}



static double _mdot_xstation (np_t *np, gl_t gl, zone_t zone, SOAP_codex_t *codex 
                    ){
  double mdot;
  long k,j,l;
#ifdef DISTMPI
  int rank;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef DISTMPI
  if (rank==0) {
#endif
    verify_post_domain_is_single_x_plane(np,gl,zone,codex,"_mdot_xstation");
    mdot=0.0;
    for3DL(k,zone.ks,zone.ke)
      for2DL(j,zone.js,zone.je)
        l=_ai(&gl,zone.is,j,k);
        if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  ){
          mdot+=_rho(np[l])*_V(np[l],0)*area_yz(np,&gl,zone.is,j,k);
        }        
      end2DL
    end3DL
#ifdef DISTMPI
  }
  MPI_Bcast(&mdot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  return(mdot);
}





/* expand the flow to the engine back pressure */
static double _Fpot_xstation (np_t *np, gl_t gl, zone_t zone, SOAP_codex_t *codex, 
                              long numsteps, double q_min, double Pback ){
  double qc,rhoc,ThrustPotential;
  long k,j,l,i;
#ifdef DISTMPI
  int rank;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef DISTMPI
  if (rank==0) {
#endif
    i=zone.is;
    verify_post_domain_is_single_x_plane(np,gl,zone,codex,"_Fpot_xstation");
    ThrustPotential=0.0;
    for3DL(k,zone.ks,zone.ke)
      for2DL(j,zone.js,zone.je)
        l=_ai(&gl,i,j,k);
        if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  ) {
          if (reversibly_expand_q_and_rho_to_given_P(np[l], &gl, Pback, &qc, &rhoc, numsteps, q_min)==0)
            ThrustPotential+=_rho(np[l])*_V(np[l],0)*area_yz(np,&gl,i,j,k)
                  /rhoc/qc*(rhoc*sqr(qc)+Pback);
        }
      end2DL
    end3DL
#ifdef DISTMPI
  }
  MPI_Bcast(&ThrustPotential,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  ThrustPotential/=_mdot_xstation(np, gl, zone, codex);
  return(ThrustPotential);
}



/* find force per mdot associated with momentum in x direction */
static double _Fx_xstation (np_t *np, gl_t gl, zone_t zone, SOAP_codex_t *codex){
  double Fx;
  long k,j,l,i;
#ifdef DISTMPI
  int rank;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef DISTMPI
  if (rank==0) {
#endif
    i=zone.is;
    verify_post_domain_is_single_x_plane(np,gl,zone,codex,"_Fx_xstation");
    Fx=0.0;
    for3DL(k,zone.ks,zone.ke)
      for2DL(j,zone.js,zone.je)
        l=_ai(&gl,i,j,k);
        if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  ) {
            Fx+=(_rho(np[l])*_V(np[l],0)*fabs(_V(np[l],0))
                +_Pstar(np[l],&gl))*area_yz(np,&gl,i,j,k);
        }
      end2DL
    end3DL
#ifdef DISTMPI
  }
  MPI_Bcast(&Fx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  Fx/=_mdot_xstation(np, gl, zone, codex);
  return(Fx);
}





/* find iteratively the engine back pressure */
static double _Pback_xstation (np_t *np, gl_t gl, zone_t zone, SOAP_codex_t *codex,
                               double Aback, long numsteps, double q_min,
                               double Pback_min, double Pback_max){
  double Pback;
  argfP_t argfP;
  long IFLAG;
  long i;
#ifdef DISTMPI
  int rank;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef DISTMPI
  if (rank==0) {
#endif
    i=zone.is;
    verify_post_domain_is_single_x_plane(np,gl,zone,codex,"_Pback_xstation");
    argfP.i=i;
    argfP.Aback=Aback;
    argfP.np=np;
    argfP.gl=&gl;
    argfP.numsteps=numsteps;
    argfP.q_min=q_min;
    argfP.zone=zone;
    /* the relative and absolute errors are on the engine exit area */
    Pback=EXM_find_root_zero_in(&fPstarc, &argfP, Pback_min, Pback_max,
                             1.0e-5, 1.0e-10, &IFLAG);
    if (IFLAG==4) wfprintf(stderr,"Problem finding Pback..\n");
#ifdef DISTMPI
  }
  MPI_Bcast(&Pback,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

  return(Pback);
}



static double _mdotreacting_xstation (np_t *np, gl_t gl, zone_t zone, SOAP_codex_t *codex,
                             long rank1,  double wstoichio1, long rank2,
			     double wstoichio2){
  double mdotreacting;
#ifdef _FLUID_MULTISPECIES
  double creacting;
  long k,j,l;
#endif
#ifdef DISTMPI
  int rank;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef DISTMPI
  if (rank==0) {
#endif
    rank1--;
    rank2--;
    verify_post_domain_is_single_x_plane(np,gl,zone,codex,"_mdotreacting_xstation");
    mdotreacting=0.0;
#ifdef _FLUID_MULTISPECIES
    for3DL(k,zone.ks,zone.ke)
      for2DL(j,zone.js,zone.je)
        l=_ai(&gl,zone.is,j,k);
        if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  ) {
          if (_w_post(np[l],rank2)>=wstoichio2) creacting=_w_post(np[l],rank1);
            else creacting=wstoichio1*_w_post(np[l],rank2)/wstoichio2;
          mdotreacting+=_rho(np[l])*_V(np[l],0)*area_yz(np,&gl,zone.is,j,k)*creacting;
        }
      end2DL
    end3DL
#endif
#ifdef DISTMPI
  }
  MPI_Bcast(&mdotreacting,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  return(mdotreacting);
}



static double _Pstag_xstation (np_t *np, gl_t gl, zone_t zone, SOAP_codex_t *codex,
                               long numsteps){
  double Pstagave;
  long k,j,l;
#ifdef DISTMPI
  int rank;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef DISTMPI
  if (rank==0) {
#endif

    verify_post_domain_is_single_x_plane(np,gl,zone,codex,"_Pstag_xstation");
    Pstagave=0.0;
    for3DL(k,zone.ks,zone.ke)
      for2DL(j,zone.js,zone.je)
        l=_ai(&gl,zone.is,j,k);
        if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  )
          Pstagave+=_rho(np[l])*_V(np[l],0)*area_yz(np,&gl,zone.is,j,k)*
                    _Pstag(np[l],&gl,numsteps);
      end2DL
    end3DL
#ifdef DISTMPI
  }
  MPI_Bcast(&Pstagave,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  Pstagave/=_mdot_xstation(np,gl,zone,codex);
  return(Pstagave);
}


static double _massavgproperty_xstation (np_t *np, gl_t gl, zone_t zone, SOAP_codex_t *codex,
                                    double(*PROPERTY)(np_t, gl_t *), char *functionname){
  double Qave;
  long k,j,l;
#ifdef DISTMPI
  int rank;
#endif

#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef DISTMPI
  if (rank==0) {
#endif
    verify_post_domain_is_single_x_plane(np,gl,zone,codex,functionname);
    Qave=0.0;
    for3DL(k,zone.ks,zone.ke)
      for2DL(j,zone.js,zone.je)
        l=_ai(&gl,zone.is,j,k);
        if (   is_node_inner(np[l],TYPELEVEL_FLUID) || is_node_bdry (np[l],TYPELEVEL_FLUID)  )
          Qave+=_rho(np[l])*_V(np[l],0)*area_yz(np,&gl,zone.is,j,k)*(*PROPERTY)(np[l],&gl);
      end2DL
    end3DL
#ifdef DISTMPI
  }
  MPI_Bcast(&Qave,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  Qave/=_mdot_xstation(np,gl,zone,codex);

  return(Qave);
}


static double _q(np_t np, gl_t *gl){
  long dim;
  double tmp;
  tmp=0.0;
  for (dim=0; dim<nd; dim++) tmp+=sqr(_V(np,dim));
  tmp=sqrt(tmp);
  return(tmp);
}


static double _rho_post(np_t np, gl_t *gl){
  double rho;
  rho=_rho(np);
  return(rho);
}


static void read_post_functions(char *functionname, char **argum,
                                char **returnstr, SOAP_codex_t *codex){
  np_t *np;
  gl_t gl;
  double Pback,Pback_min,Pback_max,Aback;
  double wstoichio1,wstoichio2;
  long rank1,rank2;
  long numsteps;
  long dim, BDRYTYPE_SURFACE;
  double q_min,mass;

  zone_t *domain_post;
  dim_t Area;
#if (fluxmom>=0)
  dim_t Fshear,Fpressure;
#endif
#if (fluxet>=0)
  double Qheat,metotal;
#endif
#ifdef EMFIELD
  dim_t Femfield;
  double Qbeam,EdotJ,Wemfield;
#endif
  int eos=EOS;
 
  np=((readcontrolarg_t *)codex->action_args)->np_post;
  gl=((readcontrolarg_t *)codex->action_args)->gl_post;
  domain_post=&(((readcontrolarg_t *)codex->action_args)->domain_post);

  read_control_functions(functionname, argum, returnstr, codex);

  if (strcmp(functionname,"_Pback_xstation")==0) {
    if (SOAP_number_argums(*argum)!=5) SOAP_fatal_error(codex,"Wrong number of arguments to _Pback_xstation.");
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%lg,%lg,%lg,%ld,%lg%n",&Aback,&Pback_min,&Pback_max,&numsteps,&q_min,&eos)!=5 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"Problem reading arguments given to function _Pback_xstation.");
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_Pback_xstation (np, gl, *domain_post, codex, Aback, numsteps, q_min, Pback_min, Pback_max));
  }

  if (strcmp(functionname,"_Fpot_xstation")==0) {
    if (SOAP_number_argums(*argum)!=3) SOAP_fatal_error(codex,"Wrong number of arguments to _Fpot_xstation.");
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%lg,%ld,%lg%n",&Pback,&numsteps,&q_min,&eos)!=3 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"Problem reading arguments given to function _Fpot_xstation.");
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_Fpot_xstation (np, gl, *domain_post, codex,
            numsteps, q_min, Pback));
  }

  if (strcmp(functionname,"_mdot_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_mdot_xstation (np, gl, *domain_post, codex ));
  }

  if (strcmp(functionname,"_Fx_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_Fx_xstation (np, gl, *domain_post, codex ));
  }

  if (strcmp(functionname,"_mdotreacting_xstation")==0) {
    if (SOAP_number_argums(*argum)!=4) SOAP_fatal_error(codex,"Wrong number of arguments to _mdotreacting_xstation.");
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld,%lg,%ld,%lg%n",&rank1,&wstoichio1,&rank2,&wstoichio2,&eos)!=4 || (*argum)[eos]!=EOS) {
      SOAP_fatal_error(codex,"Problem reading arguments given to function _mdotreacting_xstation.");
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_mdotreacting_xstation (np, gl, *domain_post, codex, rank1, wstoichio1, rank2, wstoichio2));
  }

  if (strcmp(functionname,"_Pstag_xstation")==0) {
    if (SOAP_number_argums(*argum)!=1) SOAP_fatal_error(codex,"Wrong number of arguments to _Pstag_xstation.");
    SOAP_substitute_all_argums(argum,codex);
    if (sscanf(*argum, "%ld%n",&numsteps,&eos)!=1 || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"Problem reading argument given to _Pstag_xstation().");
    }
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_Pstag_xstation (np, gl, *domain_post, codex,  numsteps));
  }

  if (strcmp(functionname,"_Tstag_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_massavgproperty_xstation(np, gl, *domain_post, codex, &(_Tstag),
               "_Tstag_xstation"));
  }

  if (strcmp(functionname,"_Pstar_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_massavgproperty_xstation(np, gl, *domain_post, codex, &(_Pstar),
               "_Pstar_xstation"));
  }

  if (strcmp(functionname,"_T_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_massavgproperty_xstation(np, gl, *domain_post, codex, &(_T),
               "_T_xstation"));
  }

  if (strcmp(functionname,"_rho_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_massavgproperty_xstation(np, gl, *domain_post, codex, &(_rho_post),
               "_rho_xstation"));
  }

  if (strcmp(functionname,"_htstar_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_massavgproperty_xstation(np, gl, *domain_post, codex, &(_htstar),
               "_htstar_xstation"));
  }

  if (strcmp(functionname,"_q_xstation")==0) {
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    sprintf(*returnstr,"%E",_massavgproperty_xstation(np, gl, *domain_post, codex, &(_q),
               "_q_xstation"));
  }


  np=*(((readcontrolarg_t *)codex->action_args)->np);
  gl=*(((readcontrolarg_t *)codex->action_args)->gl);


  if (strcmp(functionname,"_Area")==0) {
    if (SOAP_number_argums(*argum)!=2*nd+2) SOAP_fatal_error(codex,"_Area() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"], the dimension [1..3], and the boundary condition type of the surface [eg. BDRY_WALLTFIXED1, BDRY_SYMMETRICAL1, etc].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    dim=SOAP_get_argum_long(codex,*argum,2*nd);
    BDRYTYPE_SURFACE=SOAP_get_argum_long(codex,*argum,2*nd+1);
    dim--;
    if (dim<0 || dim>=nd) SOAP_fatal_error(codex,"The specified dimension is not within range when calling _Fshear().");
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_area_on_bdry(np, &gl, *domain_post, Area, BDRYTYPE_SURFACE);
    sprintf(*returnstr,"%E",Area[dim]);
  }


#if (fluxmom>=0)
  if (strcmp(functionname,"_Fshear")==0) {
    if (SOAP_number_argums(*argum)!=2*nd+2) SOAP_fatal_error(codex,"_Fshear() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"], the dimension [1..3], and the boundary condition type of the surface [eg. BDRY_WALLTFIXED1, BDRY_SYMMETRICAL1, etc].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    dim=SOAP_get_argum_long(codex,*argum,2*nd);
    BDRYTYPE_SURFACE=SOAP_get_argum_long(codex,*argum,2*nd+1);
    dim--;
    if (dim<0 || dim>=nd) SOAP_fatal_error(codex,"The specified dimension is not within range when calling _Fshear().");
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_shear_force_on_bdry(np, &gl, *domain_post, Fshear, BDRYTYPE_SURFACE);
    sprintf(*returnstr,"%E",Fshear[dim]);
  }


  if (strcmp(functionname,"_Fpressure")==0) {
    if (SOAP_number_argums(*argum)!=nd*2+2) SOAP_fatal_error(codex,"_Fpressure() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"], the dimension [1..3],  and the boundary condition type of the surface [eg. BDRY_WALLTFIXED1, BDRY_SYMMETRICAL1, etc].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    dim=SOAP_get_argum_long(codex,*argum,2*nd);
    BDRYTYPE_SURFACE=SOAP_get_argum_long(codex,*argum,2*nd+1);
    dim--;
    if (dim<0 || dim>=nd) SOAP_fatal_error(codex,"The specified dimension is not within range when calling _Fpressure().");
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_pressure_force_on_bdry(np, &gl, *domain_post, Fpressure, BDRYTYPE_SURFACE);
    sprintf(*returnstr,"%E",Fpressure[dim]);
  }
#endif

#if (fluxet>=0)
  if (strcmp(functionname,"_Qheat")==0) {
    if (SOAP_number_argums(*argum)!=2*nd+1) SOAP_fatal_error(codex,"_Qheat() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"], and the boundary condition type of the surface [eg. BDRY_WALLTFIXED1, BDRY_SYMMETRICAL1, etc].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    BDRYTYPE_SURFACE=SOAP_get_argum_long(codex,*argum,2*nd);
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_heat_to_surface_on_bdry(np, &gl, *domain_post, &Qheat, BDRYTYPE_SURFACE);
    sprintf(*returnstr,"%E",Qheat);
  }


  if (strcmp(functionname,"_metotal")==0) {
    if (SOAP_number_argums(*argum)!=2*nd) SOAP_fatal_error(codex,"_metotal() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_metotal(np, &gl, *domain_post, &metotal);
    sprintf(*returnstr,"%E",metotal);
  }
#endif

  if (strcmp(functionname,"_m")==0) {
    if (SOAP_number_argums(*argum)!=2*nd) SOAP_fatal_error(codex,"_m() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_mass(np, &gl, *domain_post, &mass);
    sprintf(*returnstr,"%E",mass);
  }


#ifdef EMFIELD
  if (strcmp(functionname,"_Femfield")==0) {
    if (SOAP_number_argums(*argum)!=1+2*nd) SOAP_fatal_error(codex,"_Femfield() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"], and the dimension [1..3].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    dim=SOAP_get_argum_long(codex,*argum,2*nd);
    dim--;
    if (dim<0 || dim>=nd) SOAP_fatal_error(codex,"The specified dimension is not within range when calling _Femfield().");
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_emfield_force(np, &gl, *domain_post, Femfield);
    sprintf(*returnstr,"%E",Femfield[dim]);
  }

  if (strcmp(functionname,"_Qbeam")==0) {
    if (SOAP_number_argums(*argum)!=2*nd) SOAP_fatal_error(codex,"_Qbeam() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_Qbeam(np, &gl, *domain_post, &Qbeam);
    sprintf(*returnstr,"%E",Qbeam);
  }

  if (strcmp(functionname,"_EdotJ")==0) {
    if (SOAP_number_argums(*argum)!=2*nd) SOAP_fatal_error(codex,"_EdotJ() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_EdotJ(np, &gl, *domain_post, &EdotJ);
    sprintf(*returnstr,"%E",EdotJ);
  }

  if (strcmp(functionname,"_Wemfield")==0) {
    if (SOAP_number_argums(*argum)!=2*nd) SOAP_fatal_error(codex,"_Wemfield() expects the following arguments: the zone limits [is,js,"if3D("ks,")" ie,je"if3D(",ke")"].");
    SOAP_substitute_all_argums(argum,codex);
    find_zone_from_argum(*argum, 0, &gl, codex, domain_post);
    *returnstr=(char *)realloc(*returnstr,40*sizeof(char));
    integrate_emfield_work(np, &gl, *domain_post, &Wemfield);
    sprintf(*returnstr,"%E",Wemfield);
  }

#endif

}


static void read_post_actions(char *action, char **argum, SOAP_codex_t *codex){
  np_t **np,**np_post;
  gl_t *gl,*gl_post;
  zone_t *domain_post;
  long cnt,numcut;
  double *xcut;
  char *filename,*postprocessor;
  bool POSTGRIDONLY;

  np=((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;
  np_post=&(((readcontrolarg_t *)codex->action_args)->np_post);
  gl_post=&(((readcontrolarg_t *)codex->action_args)->gl_post);
  domain_post=&(((readcontrolarg_t *)codex->action_args)->domain_post);
  filename=(char *)malloc(sizeof(char));
  postprocessor=(char *)malloc(sizeof(char));

  if (strcmp(action,"SetXstation")==0) {
#ifdef DISTMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    SOAP_substitute_all_argums(argum, codex);
    numcut=SOAP_number_argums(*argum);
    if (numcut>0) {
      if (*np_post!=*np) {
        /* if *np!=*np_post, then free the xcut domain first */
        free (*np_post);
      }
      xcut=(double *)malloc(numcut*sizeof(double));
      for (cnt=0; cnt<numcut; cnt++) xcut[cnt]=SOAP_get_argum_double(codex,*argum,cnt);
      create_domain_of_cuts_along_x(*np, gl, np_post, gl_post, xcut, numcut);

      *domain_post=gl_post->domain_all;
      free(xcut);
    }
    codex->ACTIONPROCESSED=TRUE;
  }



  if (strcmp(action,"WritePostFile")==0) {
#ifdef DISTMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=3+nd*2)
      SOAP_fatal_error(codex,"Not the right number of arguments in WritePostFile.");
    find_zone_from_argum(*argum, 0, gl, codex, domain_post);
    SOAP_get_argum_string(codex,&filename, *argum, 2*nd+0);
    SOAP_get_argum_string(codex,&postprocessor, *argum, 2*nd+1);
    POSTGRIDONLY=(bool)SOAP_get_argum_long(codex,*argum, 2*nd+2);
    write_post_file(*np, gl, *domain_post, filename, postprocessor, POSTGRIDONLY);
    codex->ACTIONPROCESSED=TRUE;
  }

  free(filename);
  free(postprocessor);
}


void read_post(char *argum, SOAP_codex_t *codex){
  #ifdef DISTMPI
  int rank;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank!=0) fatal_error("Post module can only be executed with one MPI process. Try -np 1.");    
  #endif
  wfprintf(stdout,"\n");
  codex->action=&read_post_actions;
  codex->function=&read_post_functions;
  ((readcontrolarg_t *)codex->action_args)->domain_post=((readcontrolarg_t *)codex->action_args)->gl->domain_all;
  ((readcontrolarg_t *)codex->action_args)->np_post=*((readcontrolarg_t *)codex->action_args)->np;
  ((readcontrolarg_t *)codex->action_args)->gl_post=*((readcontrolarg_t *)codex->action_args)->gl;
  add_bdry_types_fluid_to_codex(codex);

  SOAP_process_code(argum, codex, SOAP_VARS_KEEP_ALL);
}


