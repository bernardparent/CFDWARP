#include <src/common.h>
#include <cycle/_cycle.h>
#include <src/control.h>
#include <src/data.h>
#include <src/post.h>
#include <src/bdry.h>
#include <signal.h>
#include <stdarg.h>
#include <model/_model.h>
#include <src/version.hh>

static int SIGUSR1_CAUGHT;

#define maxstrlen 4096
#define lengthcol1 9
#define lengthcol2 14

static void signal_1 ( int signal ) {
  SIGUSR1_CAUGHT = TRUE;
}



void write_options ( FILE * outputfile ) {
  int term_width, term_height, linewidth;  
  char tmpstr[1000];
  long minlinewidth;
  minlinewidth=43;

  find_terminal_window_size ( &term_width, &term_height );
  linewidth = min ( maxlinewidth, max ( minlinewidth, term_width - 2 ) );

  write_hline ( outputfile, linewidth, 2 );
  write_options_row ( outputfile, "Flag", "Argument(s)", "Description", linewidth, lengthcol1, lengthcol2 );
  write_hline ( outputfile, linewidth, 2 );
#ifdef DISTMPI
  sprintf(tmpstr,"%d int",nd);
  write_options_row ( outputfile, "-dom", tmpstr, "Number of MPI domains in each dimension", linewidth, lengthcol1, lengthcol2 );
#endif
  write_options_row ( outputfile, "-h", "none", "List command line options", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-i", "string", "Input binary data file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-ia", "string", "Input ascii data file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-ii", "string", "Input interpolation data file", linewidth, lengthcol1, lengthcol2 );
#ifdef UNSTEADY
  write_options_row ( outputfile, "-im1", "string", "Input time level minus 1 binary data file",
                      linewidth, lengthcol1, lengthcol2 );
#if _RESTIME_BW > 2
  write_options_row ( outputfile, "-im2", "string", "Input time level minus 2 binary data file",
                      linewidth, lengthcol1, lengthcol2 );
#endif
#if _RESTIME_BW > 3
  write_options_row ( outputfile, "-im3", "string", "Input time level minus 3 binary data file",
                      linewidth, lengthcol1, lengthcol2 );
#endif
#endif
  write_options_row ( outputfile, "-l", "none", "Display licensing terms",
                      linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-ncr", "none", "No clipping report", linewidth, lengthcol1, lengthcol2 );
#if defined(EMFIELD) && defined(DISTMPI)
  write_options_row ( outputfile, "-nib", "none", 
                      "No implicit BC between MPI subdomains when updating EMField properties", linewidth, lengthcol1, lengthcol2 );
#endif
  write_options_row ( outputfile, "-nic", "none", "Set iteration count to nil", linewidth, lengthcol1, lengthcol2 );
#if defined(POSIXTHREADS) || defined(OPENMPTHREADS)
  write_options_row ( outputfile, "-nst", "none", "No short threads", linewidth, lengthcol1, lengthcol2 );
#endif
#ifdef OPENMPTHREADS
  write_options_row ( outputfile, "-nt", "int", "Set number of threads", linewidth, lengthcol1, lengthcol2 );
#endif
  write_options_row ( outputfile, "-nvr", "none", "No vars resetting at end of cycle module", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-o", "string", "Output binary data file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-oa", "string", "Output ascii data file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-oi", "string", "Output interpolation data file", linewidth, lengthcol1, lengthcol2 );
  sprintf(tmpstr,"%d int",nd+1);
  write_options_row ( outputfile, "-on", tmpstr, 
                      "Output node types nearby the node i,j,k with the bandwidth bw ("if2D("-on i j bw")if3D("-on i j k bw")")",
                      linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-op", "string", "Output post file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-opg", "string", "Output grid-only post file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-opm", "string", "Execute the post module in control file",
                      linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-or", "none", "Output maximum residual", linewidth, lengthcol1, lengthcol2 );
  sprintf(tmpstr,"%d int",nd*2);
  write_options_row ( outputfile, "-pr", tmpstr, 
                      "Post file region ("if2D("-pr is js ie je")if3D("-pr is js ks ie je ke")")", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-pt", "string", 
                      "Post-processor type: tecplot (default), tecplotnoheader, nodplot, vtk, gnuplot", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-px", "double", "Post file x-cuts (-px 0.1 0.2 0.55 ...)",
                      linewidth, lengthcol1, lengthcol2 );
#ifdef _2DL
  write_options_row ( outputfile, "-py", "double", "Post file y-cuts (-py 0.1 0.2 0.55 ...)",
                      linewidth, lengthcol1, lengthcol2 );
#endif
#ifdef _3D
  write_options_row ( outputfile, "-pz", "double", "Post file z-cuts (-pz 0.1 0.2 0.55 ...)",
                      linewidth, lengthcol1, lengthcol2 );
#endif
  write_options_row ( outputfile, "-r", "string", "Read control file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-v", "none", "Version number and other info", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-w", "string", "Write control file", linewidth, lengthcol1, lengthcol2 );
  write_options_row ( outputfile, "-q", "none", "Quit before iterating", linewidth, lengthcol1, lengthcol2 );

  write_hline ( outputfile, linewidth, 2 );

}

void fatal_error_if_options_remaining ( int argc, char **argv ) {
  char *options;
  int RET;
  options = NULL;
  RET = find_remaining_options ( argc, argv, &options );
  if ( RET == 1 ) {
    fatal_error ( "The following CFDWARP command line option could not be processed:\n%s", options );
  }
  if ( RET > 1 ) {
    fatal_error ( "The following CFDWARP command line options could not be processed:\n%s", options );
  }
  free ( options );
}

int main ( int argc, char **argv ) {
  np_t *np, *nppost;
  gl_t gl, glpost;
  zone_t postzone;
  double *xcut;
  long i, j, k, node_i, node_j, node_k, node_bw, numcut;
  char *control_filename, *post_filename, *postprocessor;
  bool INPUTBINARY, FOUNDCUT, WRITECONTROL, QUICK, OUTPUTBINARY, OUTPUTNODETYPE, OUTPUTPOSTMODULE, CYCLEMODULE,
    READCONTROL, OUTPUTPOST, OUTPUTASCII, OUTPUTINTERPOLATION, GRIDONLY, RESETITERCOUNT, OUTPUTRESIDUAL;
  input_t input;
  int *argint;
  int RET;

#ifdef DISTMPI
  int SIGUSR1_CAUGHT_ALL;
#ifdef THREADS
  int MPI_provided;
  MPI_Init_thread ( &argc, &argv, MPI_THREAD_MULTIPLE, &MPI_provided );
#else //THREADS
  MPI_Init ( &argc, &argv );
#endif //THREADS
#endif
  node_j = 0;
  node_k = 0;
  /* take care of the command line options */
  READCONTROL = FALSE;
  WRITECONTROL = FALSE;
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
  INPUTBINARY=FALSE;
  OUTPUTBINARY = FALSE;
  OUTPUTASCII = FALSE;
  OUTPUTINTERPOLATION = FALSE;
  OUTPUTPOST = FALSE;
  QUICK = FALSE;
  OUTPUTNODETYPE = FALSE;
  argint = NULL;

#if defined(UNSTEADY)
  gl.time = 0.0;
#endif

#ifdef UNSTEADY
  gl.dt = 1e99;
#endif

#ifdef EMFIELD
  gl.numsubiter_tsemf = 4;
  gl.tsemfmethod = 0;
#endif

  GRIDONLY = FALSE;
  control_filename = NULL;
  post_filename = NULL;
  input.name = NULL;
  gl.output_filename = NULL;
  gl.post_filename = NULL;
  postprocessor = ( char * ) malloc ( 10 * sizeof ( char ) );
  strcpy ( postprocessor, "tecplot" );
  if ( process_flag ( argc, argv, "-nic" ) )
    RESETITERCOUNT = TRUE;
  else
    RESETITERCOUNT = FALSE;
  if ( process_flag ( argc, argv, "-nvr" ) )
    gl.RESETRUNTIMEVARS = FALSE;
  else
    gl.RESETRUNTIMEVARS = TRUE;
  if ( process_flag ( argc, argv, "-ncr" ) )
    gl.REPORTCLIPPING = FALSE;
  else
    gl.REPORTCLIPPING = TRUE;
#if defined(POSIXTHREADS) || defined(OPENMPTHREADS)
  if ( process_flag ( argc, argv, "-nst" ) )
    gl.NOSHORTTHREADS = TRUE;
  else
    gl.NOSHORTTHREADS = FALSE;
#endif

  if ( process_flag_string ( argc, argv, "-r", &control_filename ) )
    READCONTROL = TRUE;
  if ( process_flag_string ( argc, argv, "-w", &control_filename ) )
    WRITECONTROL = TRUE;
  if ( process_flag ( argc, argv, "-h" ) ) {
    write_options ( stdout );
    free ( postprocessor ); 
    exit(EXIT_SUCCESS);
  }
  if ( process_flag ( argc, argv, "-l" ) ) {
    write_license ( stdout );
    free ( postprocessor ); 
    exit(EXIT_SUCCESS);
  }
  if (READCONTROL && WRITECONTROL) fatal_error("CFDWARP can not be called with both the '-r' and the '-w' flags.");
  if ( process_flag ( argc, argv, "-v" ) || argc==1) {
    write_modules ( stdout );
    free ( postprocessor ); 
    exit(EXIT_SUCCESS);
  }
  if (!READCONTROL && !WRITECONTROL) fatal_error("CFDWARP must be called with either the '-r' or the '-w' flag.");

  if ( process_flag_string ( argc, argv, "-i", &input.name ) ) {
    INPUTBINARY=TRUE;
    input.READDATAFILE = TRUE;
  }
  if ( process_flag_string ( argc, argv, "-ia", &input.name ) ) {
    input.ASCII = TRUE;
    input.READDATAFILE = TRUE;
  }
  if ( process_flag_string ( argc, argv, "-ii", &input.name ) ) {
    input.INTERPOLATION = TRUE;
    input.READDATAFILE = TRUE;
  }
  if (INPUTBINARY+input.ASCII+input.INTERPOLATION>1) fatal_error("CFDWARP can not input more than one data file.");

  if ( process_flag_string ( argc, argv, "-o", &gl.output_filename ) )
    OUTPUTBINARY = TRUE;
  if ( process_flag_string ( argc, argv, "-oa", &gl.output_filename ) )
    OUTPUTASCII = TRUE;
  gl.OUTPUTASCII = OUTPUTASCII;
  if ( process_flag_string ( argc, argv, "-oi", &gl.output_filename ) )
    OUTPUTINTERPOLATION = TRUE;
  gl.OUTPUTINTERPOLATION = OUTPUTINTERPOLATION;

  if ( process_flag_string ( argc, argv, "-op", &post_filename ) )
    OUTPUTPOST = TRUE;
  if ( process_flag_string ( argc, argv, "-opg", &post_filename ) ) {
    OUTPUTPOST = TRUE;
    GRIDONLY = TRUE;
  }

  RET = process_flag_int_multiple ( argc, argv, "-on", &argint );
  if ( RET ) {
    if ( RET != nd + 2 )
      fatal_error ( "After -on flag, you must supply %d integer arguments.", nd + 1 );
    node_i = argint[0];
#ifdef _2DL
    node_j = argint[1];
#endif
#ifdef _3DL
    node_k = argint[2];
#endif
    node_bw = argint[nd];
    OUTPUTNODETYPE = TRUE;
  }
  if ( process_flag ( argc, argv, "-or" ) )
    OUTPUTRESIDUAL = TRUE;
  else
    OUTPUTRESIDUAL = FALSE;
  if ( process_flag_string ( argc, argv, "-opm", &gl.post_filename ) )
    OUTPUTPOSTMODULE = TRUE;
  else
    OUTPUTPOSTMODULE = FALSE;
  if (OUTPUTBINARY+OUTPUTASCII+OUTPUTINTERPOLATION+OUTPUTPOST+OUTPUTNODETYPE+OUTPUTPOSTMODULE>1) fatal_error("CFDWARP can not output more than one file. Use only one flag within the -o family.");

  if (!OUTPUTBINARY && !OUTPUTASCII && !OUTPUTINTERPOLATION && !OUTPUTPOST && !WRITECONTROL && !OUTPUTNODETYPE && !OUTPUTRESIDUAL && !OUTPUTPOSTMODULE) 
    fatal_error("When the -r flag is specified, CFDWARP must be called with either the '-o', '-oa', '-oi', '-on', '-op', '-opg', '-opm', or '-or' flag.");

#ifdef UNSTEADY
  if ( process_flag_string ( argc, argv, "-im1", &input.name_m1 ) )
    input.M1 = TRUE;
#if _RESTIME_BW > 2
  if ( process_flag_string ( argc, argv, "-im2", &input.name_m2 ) )
    input.M2 = TRUE;
#endif
#if _RESTIME_BW > 3
  if ( process_flag_string ( argc, argv, "-im3", &input.name_m3 ) )
    input.M3 = TRUE;
#endif
#endif //UNSTEADY

#ifdef OPENMPTHREADS
  int numthread;
  if ( process_flag_int ( argc, argv, "-nt", &numthread ) ) {
    if ( numthread < 1 ) {
      fatal_error ( "The maximum number of threads must be higher than 1." );
    }
    omp_set_num_threads ( numthread );
  }
#endif
#if defined(EMFIELD) && defined(DISTMPI)
  if ( process_flag ( argc, argv, "-nib" ) )
    gl.EM_MPIBDRY_EXPLICIT = TRUE;
  else
    gl.EM_MPIBDRY_EXPLICIT = FALSE;
#endif



#ifdef DISTMPI
  gl.DISTDOMAIN_METHOD=DISTDOMAIN_METHOD_AUTO;
  RET = process_flag_int_multiple ( argc, argv, "-dom", &argint );
  if ( RET ) {
    gl.DISTDOMAIN_METHOD=DISTDOMAIN_METHOD_USERSPEC;
    if ( RET != nd + 1 )
      fatal_error ( "After -on flag, you must supply %d integer arguments.", nd );
    gl.numdomain_i = argint[0];
#ifdef _2DL
    gl.numdomain_j = argint[1];
#else
    gl.numdomain_j = 1;
#endif
#ifdef _3DL
    gl.numdomain_k = argint[2];
#else
    gl.numdomain_k = 1;
#endif
  }
#endif


  if ( process_flag ( argc, argv, "-q" ) || OUTPUTPOST || OUTPUTNODETYPE || OUTPUTRESIDUAL || OUTPUTPOSTMODULE )
    QUICK = TRUE;



  if ( READCONTROL
       && ( OUTPUTPOST || OUTPUTBINARY || OUTPUTASCII || OUTPUTINTERPOLATION || OUTPUTNODETYPE || OUTPUTPOSTMODULE
            || OUTPUTRESIDUAL ) ) {

    if ( !OUTPUTPOST )
      fatal_error_if_options_remaining ( argc, argv );
    //if (QUICK) CYCLEMODULE=FALSE; else CYCLEMODULE=TRUE;
    if (GRIDONLY) CYCLEMODULE=FALSE; else CYCLEMODULE=TRUE;
    //CYCLEMODULE = TRUE;
    read_control ( control_filename, input, CYCLEMODULE, OUTPUTPOSTMODULE, GRIDONLY, RESETITERCOUNT, &np, &gl );
    if ( !QUICK ) {
      gl.CONVERGED = FALSE;
      signal ( SIGUSR1, &signal_1 );
      do {
        SIGUSR1_CAUGHT = FALSE;
        perform_one_iteration ( np, &gl );
        gl.iter++;
#ifdef DISTMPI
        MPI_Allreduce ( &SIGUSR1_CAUGHT, &SIGUSR1_CAUGHT_ALL, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
        if ( SIGUSR1_CAUGHT_ALL > 0 )
          SIGUSR1_CAUGHT = TRUE;
        else
          SIGUSR1_CAUGHT = FALSE;
#endif
        if ( SIGUSR1_CAUGHT )
          write_data_file ( np, &gl );
      } while ( !gl.CONVERGED );
      write_data_file ( np, &gl );
      //free_clipped_variables(&gl);

    } else {

      if ( OUTPUTBINARY || OUTPUTASCII || OUTPUTINTERPOLATION )
        write_data_file ( np, &gl );
      if ( OUTPUTNODETYPE ) {
        wfprintf ( stdout, "\n\nFLUID NODE TYPES:" );
        display_node_type ( np, &gl,
                            _zone_expansion ( _zone_from_point ( node_i, node_j, node_k ), +node_bw ),
                            TYPELEVEL_FLUID );
#ifdef EMFIELD
        wfprintf ( stdout, "\n\nEMFIELD NODE TYPES:" );
        display_node_type ( np, &gl,
                            _zone_expansion ( _zone_from_point ( node_i, node_j, node_k ), +node_bw ),
                            TYPELEVEL_EMFIELD );
#endif
      }
      if ( OUTPUTRESIDUAL ) {
        wfprintf ( stdout, "Finding maximum residual..\n" );
        check_residual ( np, &gl, gl.domain );
        wfprintf ( stdout, "ximax=%E (%ld" if2DL ( ", %ld" ) if3DL ( ", %ld" ) ")\n", gl.ximax, gl.i_ximax
#ifdef _2DL
                   , gl.j_ximax
#endif
#ifdef _3DL
                   , gl.k_ximax
#endif
           );
#ifdef EMFIELD
        wfprintf ( stdout, "ximax_emfield=%E (%ld" if2DL ( ", %ld" ) if3DL ( ", %ld" ) ")\n",
                   gl.ximax_emfield, gl.i_ximax_emfield
#ifdef _2DL
                   , gl.j_ximax_emfield
#endif
#ifdef _3DL
                   , gl.k_ximax_emfield
#endif
           );
#endif
      }
      if ( OUTPUTPOST ) {
        if ( process_flag_string ( argc, argv, "-pt", &postprocessor ) ) {
          if ( !
               ( strcmp ( postprocessor, "nodplot" ) == 0 || strcmp ( postprocessor, "vtk" ) == 0
                 || strcmp ( postprocessor, "tecplot" ) == 0 || strcmp ( postprocessor, "gnuplot" ) == 0
                 || strcmp ( postprocessor, "tecplotnoheader" ) == 0 ) ) {
            fatal_error ( "Specified post-processor \"%s\" is not recognized.", postprocessor );
          }
        }

        FOUNDCUT=FALSE;

        xcut = NULL;
        RET = process_flag_double_multiple ( argc, argv, "-px", &xcut );
        if ( RET ) {
          FOUNDCUT=TRUE;
          numcut = RET - 1;
          create_domain_of_cuts_along_x ( np, &gl, &nppost, &glpost, xcut, numcut );
        } 
        free ( xcut );

#ifdef _2DL
        xcut = NULL;
        RET = process_flag_double_multiple ( argc, argv, "-py", &xcut );
        if ( RET && !FOUNDCUT) {
          FOUNDCUT=TRUE;
          numcut = RET - 1;
          create_domain_of_cuts_along_y ( np, &gl, &nppost, &glpost, xcut, numcut );
        } 
        free ( xcut );
#endif

#ifdef _3D
        xcut = NULL;
        RET = process_flag_double_multiple ( argc, argv, "-pz", &xcut );
        if ( RET && !FOUNDCUT) {
          FOUNDCUT=TRUE;
          numcut = RET - 1;
          create_domain_of_cuts_along_z ( np, &gl, &nppost, &glpost, xcut, numcut );
        } 
        free ( xcut );
#endif

        if (!FOUNDCUT) {
          glpost = gl;
          nppost = np;
        }

        postzone = glpost.domain_all;
        RET = process_flag_int_multiple ( argc, argv, "-pr", &argint );
        if ( RET ) {
          if ( RET != nd * 2 + 1 )
            fatal_error ( "After -pr flag, you must supply %d integer arguments.", nd * 2 );
#ifdef _1D
          postzone.is = argint[0];
          postzone.ie = argint[1];
#endif
#ifdef _2D
          postzone.is = argint[0];
          postzone.js = argint[1];
          postzone.ie = argint[2];
          postzone.je = argint[3];
#endif
#ifdef _3D
          postzone.is = argint[0];
          postzone.js = argint[1];
          postzone.ks = argint[2];
          postzone.ie = argint[3];
          postzone.je = argint[4];
          postzone.ke = argint[5];
#endif
//          if ( !is_zone_in_zone ( postzone, glpost.domain_all ) )
//            fatal_error ( "Zone specified after -pr flag is out of bounds." );
          if (find_zone_intersection(postzone,glpost.domain_all,&postzone))
              fatal_error ( "Zone specified after -pr flag does not intersect domain anywhere." );
        }
        fatal_error_if_options_remaining ( argc, argv );
        write_post_file ( nppost, &glpost, postzone, post_filename, postprocessor, GRIDONLY );
        if( FOUNDCUT ){
          for1DL ( i, glpost.domain_lim.is, glpost.domain_lim.ie )
            for2DL ( j, glpost.domain_lim.js, glpost.domain_lim.je )
              for3DL ( k, glpost.domain_lim.ks, glpost.domain_lim.ke )
                dispose_node ( &( nppost[_ai ( &glpost, i, j, k )] ) );
              end3DL 
            end2DL 
          end1DL 
          free ( nppost );
        }
      }
    }

    for1DL ( i, gl.domain_lim.is, gl.domain_lim.ie )
      for2DL ( j, gl.domain_lim.js, gl.domain_lim.je )
        for3DL ( k, gl.domain_lim.ks, gl.domain_lim.ke )
          dispose_node ( &( np[_ai ( &gl, i, j, k )] ) );
        end3DL 
      end2DL 
    end1DL 
    free ( np );

  }

  fatal_error_if_options_remaining ( argc, argv );

  if ( WRITECONTROL ) {
    write_control ( control_filename );
  } else if ( !GRIDONLY ) { 
    free_clipped_variables( &gl );
    free( gl.cycle.code_runtime );
    if( OUTPUTPOSTMODULE )
      SOAP_free_codex_copy ( &gl.cycle.codex );
    else
      SOAP_free_codex ( &gl.cycle.codex );
  }

  free ( postprocessor );

  free ( control_filename );
  free ( post_filename );
  free ( input.name );
  free ( gl.output_filename );
  free ( gl.post_filename );
  free ( argint );
#ifdef UNSTEADY
  free ( input.name_m1 );
#if _RESTIME_BW > 2
  free ( input.name_m2 );
#endif
#if _RESTIME_BW > 3
  free ( input.name_m3 );
#endif
#endif //UNSTEADY

#ifdef DISTMPI
  MPI_Finalize (  );
#endif
  return ( EXIT_SUCCESS );
}
