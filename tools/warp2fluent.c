// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2025-2026 Felipe Martin Rodriguez Fuentes

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

To run this tool at each time/pseudotime step for a steady or transient CFDWARP simulation, use WriteInterpolationFileZone() within the control file Cycle() with a system() call.
*/

#include "share.h"
#include <exm.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <errno.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>

#define EOS 0
#define PROFILENAMELENGTH 256

/* ===========================================================================
   -exchange : the CFDWARP side of the CFDWARP <-> FLUENT rendezvous. The control 
    file then reads, for example

       WriteInterpolationFileZone("w2f.int", ie,js+1, ie,je-1);
       system("./warp2fluent -exchange -tent");
       ReadInterpolationFileZone("f2w.int", ie,js+1, ie,je-1);

    General steps:

     1. publish the binary payload CFDWARP just wrote, as xchg/w2f.<n>.int;
     2. interpolate it onto the FLUENT face centres -- the ordinary body of
        main() below -- producing warp_output.dat;
     3. stamp xchg/w2f.<n>.ready, wait for the peer's xchg/f2w.<n>.ready, and
        stage its payload under the fixed name f2w.int that Cycle() reads.

   All paths are relative to the current directory.
   =========================================================================== */
static int    EXCHANGE=0;                   /* -exchange given                 */
static char   XDIR[400]="xchg";             /* -xchg <dir>                     */
static double XTIMEOUT=600.0;               /* -timeout <s>, was COUPLE_TIMEOUT*/
static double XPOLL=0.2;                    /* -poll <s>                       */
static long   XSEQ=0;                       /* this exchange's number n        */
static int    XDONE=0;                      /* set once step 3 has succeeded   */
static char   XSTAGED[400]="f2w.int";       /* -staged <file>, what Cycle reads*/
static int    NOUDF=0;
static int    ZONEID=15;                    /* -zone <id> -> WARP_ZONE_ID      */
static int    XINTERVAL=20;                 /* -interval <n> -> WARP_INTERVAL  */
static int    ZONEID_GIVEN=0;               /* -zone was on the command line   */
static int    XINTERVAL_GIVEN=0;            /* -interval was on the command line*/
static char   FLSPEC[500]="";               /* -fluentspecies a,b,c,...        */
static int    TRANSIENT=0;                  /* -transient given                */
static int    GLOBALDX=0;                   /* -globaldx: keep the old single  */
                                            /* mean-spacing dx1/dx2 pair       */
/* CFDWARP's own time and dt, out of the -ii header, kept for the clock check
   that -exchange -transient does against FLUENT's answer. */
static double XWARPTIME=0.0, XWARPDT=0.0;
static int    XCLOCKOK=0;                   /* the two above have been read    */
static int warp_ieq(const char *a, const char *b) {
  while (*a && *b) {
    int ca=*a, cb=*b;
    if (ca>='A' && ca<='Z') ca += 'a'-'A';
    if (cb>='A' && cb<='Z') cb += 'a'-'A';
    if (ca!=cb) return 0;
    a++; b++;
  }
  return *a=='\0' && *b=='\0';
}

/* Index of a CFDWARP species in the FLUENT mixture order given by the flag -fluentspecies */
static long fluent_species_index(const char *warpname, char *flname, long flsize) {
  char item[64];
  const char *p, *s;
  long idx=0, n;
  if (flname && flsize>0) flname[0]='\0';
  if (FLSPEC[0]=='\0') return -1;
  p = warpname;
  if (strncmp(p,"w_",2)==0) p += 2;
  s = FLSPEC;
  while (*s) {
    n = 0;
    while (*s==' ' || *s=='\t') s++;
    while (*s && *s!=',' && n<(long)sizeof(item)-1) item[n++] = *s++;
    while (n>0 && (item[n-1]==' ' || item[n-1]=='\t')) n--;
    item[n] = '\0';
    while (*s && *s!=',') s++;
    if (*s==',') s++;
    if (n>0) {
      if (warp_ieq(item,p)) {
        if (flname && flsize>0) {
          long k = n;
          if (k > flsize-1) k = flsize-1;
          memcpy(flname,item,(size_t)k);
          flname[k] = '\0';
        }
        return idx;
      }
      idx++;
    }
  }
  return -1;
}

/* Append one line to xchg/warp_exchange.log, the running record of the coupling */
static void xlog(const char *fmt, ...) {
  char name[440];
  va_list ap;
  FILE *fp;
  snprintf(name,sizeof(name),"%s/warp_exchange.log",XDIR);
  fp=fopen(name,"a");
  if (fp==NULL) return;
  fprintf(fp,"exchange %ld: ",XSEQ);
  va_start(ap,fmt);
  vfprintf(fp,fmt,ap);
  va_end(ap);
  fprintf(fp,"\n");
  fclose(fp);
}

static int xexists(const char *path) {
  return (access(path,F_OK)==0);
}

/* Remove the file the control file is about to read, so that a failed exchange
   stops CFDWARP with a message instead of letting it march on stale data. */
static void xchg_atexit(void) {
  if (!EXCHANGE || XDONE) return;
  if (unlink(XSTAGED)==0 || errno==ENOENT)
    xlog("failed before the peer answered; %s removed so that the read which "
         "follows stops CFDWARP",XSTAGED);
  else
    xlog("failed before the peer answered, and %s could not be removed: %s. "
         "CFDWARP will read stale data",XSTAGED,strerror(errno));
}

static void xchg_begin(char *input, size_t insize) {
  char name[440], dst[400];
  FILE *fp;
  struct stat st;

  if (stat(XDIR,&st)!=0 || !S_ISDIR(st.st_mode)) {
    fprintf(stderr,"\nwarp2fluent -exchange: the exchange directory %s does not exist in the\n"
                   "working directory %s. Create it, and empty it, before the run: the\n"
                   "rendezvous keeps every message it has ever sent, so a stale xchg/ makes\n"
                   "this process answer exchange 0 with the previous run's reply.\n",
            XDIR,getcwd(name,sizeof(name))?name:"?");
    exit (EXIT_FAILURE);
  }
  snprintf(name,sizeof(name),"%s/seq",XDIR);
  fp=fopen(name,"r");
  if (fp!=NULL) {
    if (fscanf(fp,"%ld",&XSEQ)!=1) XSEQ=0;
    fclose(fp);
  }
  snprintf(name,sizeof(name),"%s/w2f.%ld.ready",XDIR,XSEQ);
  if (xexists(name)) {
    fprintf(stderr,"\nwarp2fluent -exchange: %s is already there, and exchange %ld is the one\n"
                   "about to be published. The rendezvous never deletes anything, so %s holds\n"
                   "messages from an earlier attempt: the FLUENT side will answer that stale\n"
                   "request instead of waiting for this one. Empty %s and restart both codes.\n"
                   "Removing %s/seq on its own does not do it -- it only renumbers this side\n"
                   "back onto messages that are still lying there.\n",
            name,XSEQ,XDIR,XDIR,XDIR);
    xlog("refused to start: %s is a leftover from an earlier run",name);
    exit (EXIT_FAILURE);
  }
  snprintf(name,sizeof(name),"%s/f2w.%ld.ready",XDIR,XSEQ);
  if (xexists(name)) {
    fprintf(stderr,"\nwarp2fluent -exchange: %s is already there before exchange %ld has been\n"
                   "asked. A reply cannot precede its request, so %s holds messages from an\n"
                   "earlier attempt and this exchange would return the previous run's answer\n"
                   "at once. Empty %s and restart both codes.\n",
            name,XSEQ,XDIR,XDIR);
    xlog("refused to start: %s is a leftover from an earlier run",name);
    exit (EXIT_FAILURE);
  }
  if (!xexists(input)) {
    fprintf(stderr,"\nwarp2fluent -exchange: %s is not in the working directory. That file is\n"
                   "written by the WriteInterpolationFile() call which must come immediately\n"
                   "before this system() call in Cycle(), and its name must match -ii.\n",input);
    exit (EXIT_FAILURE);
  }
  if (snprintf(dst,sizeof(dst),"%s/w2f.%ld.int",XDIR,XSEQ)>=(int)sizeof(dst) || insize<sizeof(dst)) {
    fprintf(stderr,"\nwarp2fluent -exchange: the exchange directory name %s is too long\n",XDIR);
    exit (EXIT_FAILURE);
  }
  if (rename(input,dst)!=0) {
    fprintf(stderr,"\nwarp2fluent -exchange: could not publish %s as %s: %s\n",
            input,dst,strerror(errno));
    exit (EXIT_FAILURE);
  }
  snprintf(input,insize,"%s",dst);

  snprintf(name,sizeof(name),"%s/warp2fluent.%ld.log",XDIR,XSEQ);
  if (freopen(name,"w",stdout)==NULL)
    fprintf(stderr,"warp2fluent -exchange: could not open %s, leaving stdout alone\n",name);
}

/* Step 3.  Everything below runs after warp_output.dat is on disk. */
static void xchg_finish(void) {
  char tmp[460], stamp[460], want[460], payload[460], err[460], stop[460], seqf[460];
  char buf[8192];
  FILE *fp, *in, *out;
  size_t nread;
  double waited=0.0;
  struct timespec ts;
  int c;

  snprintf(err ,sizeof(err ),"%s/ERROR",XDIR);
  snprintf(stop,sizeof(stop),"%s/STOP" ,XDIR);
  snprintf(want,sizeof(want),"%s/f2w.%ld.ready",XDIR,XSEQ);

  snprintf(tmp  ,sizeof(tmp  ),"%s/.stamp.%ld",XDIR,(long)getpid());
  snprintf(stamp,sizeof(stamp),"%s/w2f.%ld.ready",XDIR,XSEQ);
  fp=fopen(tmp,"w");
  if (fp==NULL || fclose(fp)!=0 || rename(tmp,stamp)!=0) {
    fprintf(stderr,"\nwarp2fluent -exchange: could not publish the stamp %s: %s\n",
            stamp,strerror(errno));
    xlog("could not publish the stamp %s: %s",stamp,strerror(errno));
    return;                     /* atexit removes the staged file */
  }

  ts.tv_sec =(time_t)XPOLL;
  ts.tv_nsec=(long)((XPOLL-(double)ts.tv_sec)*1.0e9);
  while (!xexists(want)) {
    if (xexists(err)) {
      /* The FLUENT side found data it refuses to send. It publishes the reason
         and no payload */
      fprintf(stderr,"\n");
      in=fopen(err,"r");
      if (in!=NULL) {
        while ((c=fgetc(in))!=EOF) fputc(c,stderr);
        fclose(in);
      }
      fprintf(stderr,"\n");
      xlog("peer refused to send; see %s",err);
      return;
    }
    if (xexists(stop)) {
      /* leave the previous staged file in place so that the read
         which follows is a re-read  */
      xlog("STOP seen, leaving %s untouched",XSTAGED);
      XDONE=1;
      return;
    }
    nanosleep(&ts,NULL);
    waited+=XPOLL;
    if (waited>XTIMEOUT) {
      fprintf(stderr,"\nwarp2fluent -exchange: no %s after %.1f s. Either the FLUENT side is not\n"
                     "iterating, or its Execute At End hook list does not hold warp_exchange, or\n"
                     "it is watching a different xchg/. Raise the wait with -timeout if the peer\n"
                     "is merely slow.\n",want,waited);
      if (TRANSIENT)
        fprintf(stderr,"On a transient run there is a fourth cause, and it is the commonest:\n"
                       "warp_wait_and_load from UDF.c is not selected under Function Hooks ->\n"
                       "Adjust. The FLUENT console then says NOT answering exchange %ld, and\n"
                       "%s/w2f.%ld.installed is missing.\n",XSEQ,XDIR,XSEQ);
      xlog("TIMEOUT after %.1f s, no %s",waited,want);
      return;
    }
  }

  snprintf(payload,sizeof(payload),"%s/f2w.%ld.int",XDIR,XSEQ);
  in =fopen(payload,"rb");
  out=fopen(XSTAGED,"wb");
  if (in==NULL || out==NULL) {
    fprintf(stderr,"\nwarp2fluent -exchange: the peer stamped %s but %s could not be staged as %s: %s\n",
            want,payload,XSTAGED,strerror(errno));
    xlog("could not stage %s as %s: %s",payload,XSTAGED,strerror(errno));
    if (in !=NULL) fclose(in);
    if (out!=NULL) fclose(out);
    return;
  }
  while ((nread=fread(buf,1,sizeof(buf),in))>0) {
    if (fwrite(buf,1,nread,out)!=nread) {
      fprintf(stderr,"\nwarp2fluent -exchange: short write staging %s as %s\n",payload,XSTAGED);
      xlog("short write staging %s as %s",payload,XSTAGED);
      fclose(in); fclose(out);
      return;
    }
  }
  fclose(in);
  if (fclose(out)!=0) {
    fprintf(stderr,"\nwarp2fluent -exchange: could not close %s: %s\n",XSTAGED,strerror(errno));
    xlog("could not close %s: %s",XSTAGED,strerror(errno));
    return;
  }

  /* -transient: both interpolation-file headers carry time= and dt=; CFDWARP's came out of
     the -ii file this call published and FLUENT's is CURRENT_TIME and
     CURRENT_TIMESTEP at the end of the time step it answered from */
  if (TRANSIENT) {
    double ftime=0.0, fdt=0.0;
    char hdr[2048], *at;
    int gott=0, gotd=0;
    in=fopen(XSTAGED,"rb");
    if (in!=NULL) {
      if (fgets(hdr,sizeof(hdr),in)!=NULL) {
        at=strstr(hdr," time="); if (at!=NULL) gott=(sscanf(at+6,"%lf",&ftime)==1);
        at=strstr(hdr," dt=");   if (at!=NULL) gotd=(sscanf(at+4,"%lf",&fdt)==1);
      }
      fclose(in);
    }
    if (!XCLOCKOK)
      fprintf(stdout,"CLOCK: exchange %ld: no time/dt was read from the CFDWARP side.\n",XSEQ);
    else if (!gott || !gotd)
      fprintf(stdout,"CLOCK: exchange %ld: no time/dt in the header of %s.\n",XSEQ,XSTAGED);
    else {
      fprintf(stdout,"CLOCK: exchange %ld: CFDWARP t=%.6E dt=%.6E, FLUENT t=%.6E dt=%.6E\n",
              XSEQ,XWARPTIME,XWARPDT,ftime,fdt);
      if (XWARPDT>=1.0e98)
        /* A steady CFDWARP run writes dt=1E+99 into its own interpolation files */
        fprintf(stdout,"CLOCK: WARNING the CFDWARP side is a STEADY run: dt=%.3E is what\n"
                       "       steady CFDWARP writes, and its time never advances. The two\n"
                       "       halves then share no physical time, so FLUENT's time steps are\n"
                       "       sub-iterations of a steady coupling, not a transient one. Run\n"
                       "       CFDWARP under DualTimeStepping with the same dt as FLUENT.\n",
                XWARPDT);
      else if (fdt<=0.0 || fdt>=1.0)
        fprintf(stdout,"CLOCK: WARNING FLUENT reports dt=%.6E. In a STEADY FLUENT run\n"
                       "       CURRENT_TIMESTEP is 1 s and CURRENT_TIME is 0, so a dt of 1\n"
                       "       here means FLUENT is not transient at all.\n",fdt);
      else if (XWARPDT>0.0 && fabs(fdt-XWARPDT)>1.0e-6*fabs(XWARPDT))
        fprintf(stdout,"CLOCK: WARNING the two dt differ by %.3E s (%.2f per cent of CFDWARP's).\n"
                       "       The exchange is then not time accurate: each side is advancing\n"
                       "       a different amount of physical time per exchange.\n",
                fdt-XWARPDT,100.0*(fdt-XWARPDT)/XWARPDT);
      if (XWARPDT>0.0 && XWARPDT<1.0e98 && fabs(ftime-XWARPTIME)>0.5*XWARPDT)
        fprintf(stdout,"CLOCK: WARNING the two clocks are %.3E s apart, more than half a step.\n"
                       "       One side has taken time steps the other did not.\n",
                ftime-XWARPTIME);
    }
  }

  snprintf(seqf,sizeof(seqf),"%s/seq",XDIR);
  fp=fopen(seqf,"w");
  if (fp!=NULL) { fprintf(fp,"%ld\n",XSEQ+1); fclose(fp); }
  xlog("answered after %.1f s",waited);
  XDONE=1;
}


/* Turn a CFDWARP variable name into a valid C identifier that FLUENT can use as a profile name */
static void find_valid_c_name(const char *name, char *cname, size_t cnamesize) {
  size_t i, k;
  k = 0;
  for (i = 0; name[i] != '\0' && k + 6 < cnamesize; i++) {
    char c = name[i];
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_') {
      cname[k++] = c;
    } else if (c == '+') {
      strcpy(&cname[k], "plus");
      k += 4;
    } else if (c == '-') {
      strcpy(&cname[k], "minus");
      k += 5;
    } else if (c == ')') {
      /* dropped */
    } else {
      /* '(' and anything else */
      if (k > 0 && cname[k - 1] != '_') cname[k++] = '_';
    }
  }
  while (k > 0 && cname[k - 1] == '_') k--;
  cname[k] = '\0';
  if (cname[0] == '\0' || (cname[0] >= '0' && cname[0] <= '9')) {
    memmove(cname + 1, cname, strlen(cname) + 1);
    cname[0] = 'v';
  }
}


/*
Build the list of profile names, one per column of the interpolated data file
*/
static void find_profile_names(char **varnames, long numvars, char **cnames) {
  long i, j, n;
  for (i = 0; i < numvars; i++) {
    find_valid_c_name(varnames[i], cnames[i], PROFILENAMELENGTH);
    n = 1;
    for (j = 0; j < i; j++) {
      if (strcmp(cnames[i], cnames[j]) == 0) {
        n++;
        snprintf(cnames[i] + strlen(cnames[i]), PROFILENAMELENGTH - strlen(cnames[i]), "_%ld", n);
        j = -1;
      }
    }
  }
}


/*
This function writes the UDF that is loaded through Define -> User-Defined ->
Functions->Compiled in the FLUENT GUI.

varnames/cnames hold the fluid variables followed by the EMField ones, in the
same order as the columns of the interpolated data file.
*/
void write_udf_file(const char *udf_filename, char **varnames, char **cnames, long numvars,
                    long numdim, const char *interpolated_filename, long numfaces) {
  FILE *udf;
  long i;
  char abspath[PATH_MAX + 1];
  const char *basename_ptr, *slash;

  udf = fopen(udf_filename, "w");
  if (!udf) {
    fprintf(stderr, "Could not open %s for writing\n", udf_filename);
    return;
  }

  if (realpath(interpolated_filename, abspath) == NULL) {
    strncpy(abspath, interpolated_filename, PATH_MAX);
    abspath[PATH_MAX] = '\0';
  }
  basename_ptr = interpolated_filename;
  for (slash = interpolated_filename; *slash != '\0'; slash++) {
    if (*slash == '/' || *slash == '\\') basename_ptr = slash + 1;
  }

  fprintf(udf,
    "/* -------------------------------------------------------------------------\n"
    "   FLUENT UDF generated by CFDWARP tools/warp2fluent.c -- do not edit by hand,\n"
    "   it is rewritten every time warp2fluent runs.\n"
    "\n"
    "   HOW TO LOAD IT\n"
    "     1. put this file and %s in the directory FLUENT was\n"
    "        started from (File -> Read -> Case shows it; /file/show-configuration\n"
    "        in the text console prints the working directory).\n"
    "     2. Define -> User-Defined -> Functions -> Compiled..., pick this file,\n"
    "        leave CPP Command Name at cpp, leave Display Assembly Listing off,\n"
    "        press Compile. No editing of this file is needed.\n"
    "     3. Define -> User-Defined -> Execute on Demand -> warp_load_data, to read\n"
    "        the data file and confirm in the console how many points came in.\n"
    "     4. in the boundary condition panel of the inflow zone, pick the matching\n"
    "        entry from each drop-down list (they all end in _profile).\n"
    "\n"
    "   Re-running warp2fluent for a new time step only needs step 3 again: the UDF\n"
    "   does not have to be recompiled as long as the variable list is unchanged.\n"
    "\n"
    "   IF THE INTERPRETER REFUSES SOMETHING, in the order worth trying:\n"
    "     - a complaint about Message or CX_Message: the WARP_PRINT line below is\n"
    "       already set to printf; set it to warp_noop if printf is refused too.\n"
    "     - a complaint about warp_thread or Thread: set warp_usecache to 0. The\n"
    "       UDF then searches per face per call, which is correct but slower.\n"
    "     - anything else: load it as a compiled UDF instead. The file is written\n"
    "       so that the same source works either way; for a compiled build set\n"
    "       WARP_PRINT to Message, which is also what a parallel run needs.\n"
    "\n"
    "   NOTE ON PRESSURE: P below is the absolute static pressure in Pa. FLUENT\n"
    "   pressure boundary conditions are specified as gauge pressure relative to\n"
    "   the operating pressure, so either set the operating pressure to 0 or\n"
    "   subtract it yourself.\n"
    "\n"
    "   PROFILES DEFINED HERE (CFDWARP name -> FLUENT profile name)\n",
    basename_ptr);
  for (i = 0; i < numvars; i++) {
    fprintf(udf, "     %-20s -> %s_profile\n", varnames[i], cnames[i]);
  }
  fprintf(udf,
    "   ------------------------------------------------------------------------- */\n\n");

  fprintf(udf, "#include \"udf.h\"\n\n");

  fprintf(udf,
    "/* The data file is looked for in FLUENT's working directory first, then at\n"
    "   the absolute path it was generated at. Change these two lines only if you\n"
    "   moved the file somewhere else. On Windows write the path as\n"
    "   \"C:\\\\Users\\\\you\\\\case\\\\%s\". */\n", basename_ptr);
  fprintf(udf, "#define WARP_FILE     \"%s\"\n", basename_ptr);
  fprintf(udf, "#define WARP_FILE_ALT \"%s\"\n\n", abspath);

  /* Stamp the UDF with where and when it was generated. */
  {
    char gendir[PATH_MAX + 1];
    char gentime[64];
    char *cut;
    time_t now;
    struct tm *lt;
    strncpy(gendir, abspath, PATH_MAX);
    gendir[PATH_MAX] = '\0';
    cut = strrchr(gendir, '/');
    if (cut != NULL) { if (cut == gendir) cut[1] = '\0'; else *cut = '\0'; }
    else strcpy(gendir, ".");
    now = time(NULL);
    lt = localtime(&now);
    if (lt == NULL || strftime(gentime, sizeof(gentime), "%Y-%m-%d %H:%M:%S", lt) == 0)
      strcpy(gentime, "unknown");
    fprintf(udf, "/* Where and when this file was generated. Nothing reads these but the\n"
                 "   banner in warp_load(); they are here so that the FLUENT console can\n"
                 "   name the UDF it is actually running. */\n");
    fprintf(udf, "#define WARP_GEN_DIR  \"%s\"\n", gendir);
    fprintf(udf, "#define WARP_GEN_TIME \"%s\"\n\n", gentime);
  }

  fprintf(udf,
    "/* Where the messages of this UDF go. FLUENT's interpreter on some builds\n"
    "   refuses Message, because Message resolves to CX_Message and that symbol is\n"
    "   looked up when the file is loaded; printf is built into the interpreter and\n"
    "   works there. Set the line below to exactly one of:\n"
    "       printf      interpreted UDF, serial              (default)\n"
    "       Message     compiled UDF, and any parallel run\n"
    "       warp_noop   silent, if the interpreter refuses both\n"
    "   Nothing else in this file depends on the choice. */\n"
    "#define WARP_PRINT printf\n"
    "\n"
    "/* Old-style empty parameter list on purpose: it takes whatever WARP_PRINT is\n"
    "   called with and discards it. Unused unless WARP_PRINT is set to it. */\n"
    "void warp_noop() { }\n\n");

  fprintf(udf, "/* set to 0 if the interpreter refuses the cached face map below */\n");
  fprintf(udf, "int warp_usecache = 1;\n\n");

  fprintf(udf, "#define WARP_NPTS     %ld   /* rows in %s */\n", numfaces, basename_ptr);
  fprintf(udf, "#define WARP_NVAR     %ld\n", numvars);
  fprintf(udf, "#define WARP_MAXFACE  %ld   /* face-index bound of the cached map */\n\n",
          2 * numfaces + 64);

  fprintf(udf, "int   warp_ndata  = 0;\n");
  fprintf(udf, "int   warp_loaded = 0;\n");
  fprintf(udf,
    "/* The profile data are held as WARP_REAL and read with WARP_SCANF. Both were\n"
    "   float and \"%%e\" up to 2026-08-24, which is about 7 significant digits, and\n"
    "   that is not enough for a run coupled through this file: the pressure sent\n"
    "   back to CFDWARP then wobbles at the 1e-3 Pa level from one exchange to the\n"
    "   next, CFDWARP's ximax floors two decades above its own xiverge, and the\n"
    "   coupled run cannot terminate on its own convergence test. Put both lines\n"
    "   back to float and \"%%e\" if the interpreter refuses the double form. */\n"
    "#define WARP_REAL  double\n"
    "#define WARP_SCANF \"%%le\"\n"
    "\n");
  fprintf(udf, "WARP_REAL warp_d2     = 0.0;   /* squared distance found by warp_nearest */\n");
  fprintf(udf, "WARP_REAL warp_x[WARP_NPTS];\n");
  fprintf(udf, "WARP_REAL warp_y[WARP_NPTS];\n");
  if (numdim == 3) fprintf(udf, "WARP_REAL warp_z[WARP_NPTS];\n");
  fprintf(udf, "int   warp_sorted[WARP_NPTS];\n");
  for (i = 0; i < numvars; i++) {
    fprintf(udf, "WARP_REAL warp_%s[WARP_NPTS];   /* %s */\n", cnames[i], varnames[i]);
  }
  fprintf(udf, "\n");

  /* ---------------- reader ---------------- */
  fprintf(udf,
    "/* read the interpolated data file, then index-sort the points by x so that\n"
    "   warp_nearest below can binary-search instead of scanning all of them */\n"
    "void warp_load(void)\n"
    "{\n"
    "    FILE *fp;\n"
    "    int n, gap, i, j, t;\n"
    "    int from_alt;\n"
    "    WARP_REAL extra;\n"
    "\n"
    "    warp_ndata  = 0;\n"
    "    warp_loaded = 0;\n"
    "    from_alt    = 0;\n"
    "    fp = fopen(WARP_FILE, \"r\");\n"
    "    if (fp == NULL) {\n"
    "        fp = fopen(WARP_FILE_ALT, \"r\");\n"
    "        if (fp != NULL) from_alt = 1;\n"
    "    }\n"
    "    if (fp == NULL) {\n"
    "        WARP_PRINT(\"warp2fluent: ERROR could not open %%s in the working directory\\n\", WARP_FILE);\n"
    "        WARP_PRINT(\"warp2fluent:       nor %%s\\n\", WARP_FILE_ALT);\n"
    "        WARP_PRINT(\"warp2fluent:       this UDF was generated %%s in %%s\\n\", WARP_GEN_TIME, WARP_GEN_DIR);\n"
    "        WARP_PRINT(\"warp2fluent:       if that is not the folder you are running in, FLUENT is\\n\");\n"
    "        WARP_PRINT(\"warp2fluent:       compiling an old UDF.c and not the one you just generated\\n\");\n"
    "        return;\n"
    "    }\n"
    "    n = 0;\n"
    "    while (n < WARP_NPTS) {\n");
  if (numdim == 3) {
    fprintf(udf, "        if (fscanf(fp, WARP_SCANF \" \" WARP_SCANF \" \" WARP_SCANF, &warp_x[n], &warp_y[n], &warp_z[n]) != 3) break;\n");
  } else {
    fprintf(udf, "        if (fscanf(fp, WARP_SCANF \" \" WARP_SCANF, &warp_x[n], &warp_y[n]) != 2) break;\n");
  }
  for (i = 0; i < numvars; i++) {
    fprintf(udf, "        if (fscanf(fp, WARP_SCANF, &warp_%s[n]) != 1) break;\n", cnames[i]);
  }
  fprintf(udf,
    "        n++;\n"
    "    }\n"
    "    /* The loop above stops at WARP_NPTS, so a file with MORE rows than this\n"
    "       UDF was written for used to load silently, taking the first WARP_NPTS\n"
    "       of them -- which is the wrong face for every one of them.  That is what\n"
    "       happens when warp2fluent is re-run against a different -m mesh than the\n"
    "       one the UDF was generated from.  One more read tells the two apart. */\n"
    "    if (n == WARP_NPTS && fscanf(fp, WARP_SCANF, &extra) == 1) n = WARP_NPTS + 1;\n"
    "    fclose(fp);\n"
    "    if (n > WARP_NPTS) {\n"
    "        WARP_PRINT(\"warp2fluent: ERROR %%s holds more than the %%d rows this UDF was\\n\", WARP_FILE, WARP_NPTS);\n"
    "        WARP_PRINT(\"warp2fluent:       written for. It was made from a different FLUENT mesh:\\n\");\n"
    "        WARP_PRINT(\"warp2fluent:       re-run warp2fluent with the -m mesh of THIS zone and\\n\");\n"
    "        WARP_PRINT(\"warp2fluent:       reload both UDFs. Nothing has been loaded.\\n\");\n"
    "        warp_ndata = 0;\n"
    "        return;\n"
    "    }\n"
    "    warp_ndata = n;\n"
    "    if (n != WARP_NPTS) {\n"
    "        WARP_PRINT(\"warp2fluent: WARNING read %%d complete rows but this UDF was written for %%d.\\n\", n, WARP_NPTS);\n"
    "        WARP_PRINT(\"warp2fluent:         the data file and the UDF are out of step -- regenerate both.\\n\");\n"
    "    }\n"
    "    if (n <= 0) return;\n"
    "\n"
    "    for (i = 0; i < n; i++) warp_sorted[i] = i;\n"
    "    gap = n / 2;\n"
    "    while (gap > 0) {\n"
    "        for (i = gap; i < n; i++) {\n"
    "            t = warp_sorted[i];\n"
    "            j = i;\n"
    "            while (j >= gap && warp_x[warp_sorted[j - gap]] > warp_x[t]) {\n"
    "                warp_sorted[j] = warp_sorted[j - gap];\n"
    "                j = j - gap;\n"
    "            }\n"
    "            warp_sorted[j] = t;\n"
    "        }\n"
    "        gap = gap / 2;\n"
    "    }\n"
    "    warp_loaded = 1;\n"
    "    WARP_PRINT(\"warp2fluent: read %%d points, %%d variables per point\\n\", n, WARP_NVAR);\n"
    "    WARP_PRINT(\"warp2fluent: UDF generated %%s in %%s\\n\", WARP_GEN_TIME, WARP_GEN_DIR);\n"
    "    if (from_alt) {\n"
    "        WARP_PRINT(\"warp2fluent: WARNING %%s was not in the working directory, so the\\n\", WARP_FILE);\n"
    "        WARP_PRINT(\"warp2fluent: WARNING profiles above were read from the absolute path\\n\");\n"
    "        WARP_PRINT(\"warp2fluent: WARNING %%s\\n\", WARP_FILE_ALT);\n"
    "        WARP_PRINT(\"warp2fluent: WARNING frozen into this UDF when it was generated. That\\n\");\n"
    "        WARP_PRINT(\"warp2fluent: WARNING file belongs to another run and is not updated by\\n\");\n"
    "        WARP_PRINT(\"warp2fluent: WARNING the exchange, so the coupling is not converging on\\n\");\n"
    "        WARP_PRINT(\"warp2fluent: WARNING anything. Regenerate the UDFs in this folder.\\n\");\n"
    "    } else {\n"
    "        WARP_PRINT(\"warp2fluent: %%s read from the working directory\\n\", WARP_FILE);\n"
    "    }\n"
    "}\n\n");

  /* ---------------- nearest point ---------------- */
  if (numdim == 3) {
    fprintf(udf, "int warp_nearest(float xf, float yf, float zf)\n");
  } else {
    fprintf(udf, "int warp_nearest(float xf, float yf)\n");
  }
  fprintf(udf,
    "{\n"
    "    int lo, hi, mid, l, r, k, best;\n"
    "    float dx, dy%s, d2, bd2;\n"
    "\n"
    "    best = -1;\n"
    "    bd2  = 1.0e30;\n"
    "    warp_d2 = 0.0;\n"
    "    if (warp_ndata <= 0) return -1;\n"
    "\n"
    "    lo = 0;\n"
    "    hi = warp_ndata;\n"
    "    while (lo < hi) {\n"
    "        mid = (lo + hi) / 2;\n"
    "        if (warp_x[warp_sorted[mid]] < xf) lo = mid + 1; else hi = mid;\n"
    "    }\n"
    "    r = lo;\n"
    "    l = lo - 1;\n"
    "    while (r < warp_ndata || l >= 0) {\n"
    "        if (r < warp_ndata) {\n"
    "            k  = warp_sorted[r];\n"
    "            dx = warp_x[k] - xf;\n"
    "            if (dx * dx > bd2) {\n"
    "                r = warp_ndata;\n"
    "            } else {\n"
    "                dy = warp_y[k] - yf;\n"
    "%s"
    "                if (d2 < bd2) { bd2 = d2; best = k; }\n"
    "                r = r + 1;\n"
    "            }\n"
    "        }\n"
    "        if (l >= 0) {\n"
    "            k  = warp_sorted[l];\n"
    "            dx = warp_x[k] - xf;\n"
    "            if (dx * dx > bd2) {\n"
    "                l = -1;\n"
    "            } else {\n"
    "                dy = warp_y[k] - yf;\n"
    "%s"
    "                if (d2 < bd2) { bd2 = d2; best = k; }\n"
    "                l = l - 1;\n"
    "            }\n"
    "        }\n"
    "    }\n"
    "    warp_d2 = bd2;\n"
    "    return best;\n"
    "}\n\n",
    (numdim == 3 ? ", dz" : ""),
    (numdim == 3 ? "                dz = warp_z[k] - zf;\n"
                   "                d2 = dx * dx + dy * dy + dz * dz;\n"
                 : "                d2 = dx * dx + dy * dy;\n"),
    (numdim == 3 ? "                dz = warp_z[k] - zf;\n"
                   "                d2 = dx * dx + dy * dy + dz * dz;\n"
                 : "                d2 = dx * dx + dy * dy;\n"));

  /* ---------------- face map ---------------- */
  fprintf(udf,
    "/* The nearest-point search only depends on the geometry, so it is done once\n"
    "   per zone and reused by every profile on that zone and by every iteration.\n"
    "   Without this, each of the %ld profiles repeats the search on every call. */\n"
    "Thread *warp_thread = NULL;\n"
    "int     warp_mapok  = 0;\n"
    "int     warp_map[WARP_MAXFACE];\n\n", numvars);

  fprintf(udf,
    "void warp_prepare(Thread *t)\n"
    "{\n"
    "    face_t f;\n"
    "    real xyz[ND_ND];\n"
    "    int idx, n, over;\n"
    "    float worst;\n"
    "\n"
    "    if (!warp_loaded) warp_load();\n"
    "    if (warp_usecache && warp_loaded && t != warp_thread) {\n"
    "        warp_thread = t;\n"
    "        warp_mapok  = 1;\n"
    "        n = 0; over = 0; worst = 0.0;\n"
    "        begin_f_loop(f, t)\n"
    "        {\n"
    "            F_CENTROID(xyz, f, t);\n");
  if (numdim == 3) {
    fprintf(udf, "            idx = warp_nearest(xyz[0], xyz[1], (ND_ND > 2 ? xyz[2] : 0.0));\n");
  } else {
    fprintf(udf, "            idx = warp_nearest(xyz[0], xyz[1]);\n");
  }
  fprintf(udf,
    "            if (warp_d2 > worst) worst = warp_d2;\n"
    "            if (f < WARP_MAXFACE) warp_map[f] = idx; else over = 1;\n"
    "            n = n + 1;\n"
    "        }\n"
    "        end_f_loop(f, t)\n"
    "        if (over) {\n"
    "            warp_mapok = 0;\n"
    "            WARP_PRINT(\"warp2fluent: WARNING zone has more faces than WARP_MAXFACE, map disabled\\n\");\n"
    "        }\n"
    "        WARP_PRINT(\"warp2fluent: matched %%d faces of this zone, largest gap %%.3e m\\n\", n, sqrt(worst));\n"
    "    }\n"
    "}\n\n");

  fprintf(udf,
    "int warp_index(Thread *t, face_t f)\n"
    "{\n"
    "    real xyz[ND_ND];\n"
    "\n"
    "    if (warp_usecache && warp_mapok && t == warp_thread && f < WARP_MAXFACE) return warp_map[f];\n"
    "    F_CENTROID(xyz, f, t);\n");
  if (numdim == 3) {
    fprintf(udf, "    return warp_nearest(xyz[0], xyz[1], (ND_ND > 2 ? xyz[2] : 0.0));\n");
  } else {
    fprintf(udf, "    return warp_nearest(xyz[0], xyz[1]);\n");
  }
  fprintf(udf, "}\n\n");

  /* ---------------- on-demand loader ---------------- */
  fprintf(udf,
    "/* Execute on Demand entry: reloads the data file and drops the cached map,\n"
    "   so a new warp2fluent output can be picked up without recompiling. */\n"
    "DEFINE_ON_DEMAND(warp_load_data)\n"
    "{\n"
    "    warp_load();\n"
    "    warp_thread = NULL;\n"
    "    warp_mapok  = 0;\n"
    "}\n\n");

  /* ---------------- automatic reload during a coupled run ---------------- */
  fprintf(udf,
    "/* ---- automatic reload during a coupled run -------------------------------\n"
    "   Does by itself what pressing warp_load_data does by hand: it watches the\n"
    "   xchg/w2f.<n>.ready stamps that the CFDWARP-side exchange.sh publishes and\n"
    "   re-reads WARP_FILE the moment exchange n is on disk.  Without it FLUENT\n"
    "   would serve the profile of exchange 0 for the whole run.\n"
    "\n"
    "   It sits in THIS file, not in UDF_fluent2warp.c, because warp_load(),\n"
    "   warp_thread and warp_mapok are file-scope symbols of this file and an\n"
    "   INTERPRETED UDF cannot see another interpreted file's symbols.  The two\n"
    "   generated files coordinate through the exchange directory, never through C,\n"
    "   so each one loads on its own and both work interpreted or compiled.\n"
    "\n"
    "   SELECT BOTH Execute-At-End functions under Define -> User-Defined ->\n"
    "   Function Hooks -> Execute At End: warp_reload_at_end from this file and\n"
    "   warp_exchange from UDF_fluent2warp.c.  That list accepts several entries.\n"
    "   Set WARP_COUPLED to 0 for a standalone FLUENT run with no CFDWARP on the\n"
    "   other end; the profiles then only change when you press warp_load_data.\n"
    "   The function below is compiled either way -- see the note on it.\n"
    "\n"
    "   IN A TRANSIENT RUN warp_reload_at_end is NOT the hook that installs the\n"
    "   profile -- warp_wait_and_load below is, and it has to be selected under\n"
    "   Function Hooks -> Adjust.  See the block above it for why. */\n"
    "#define WARP_COUPLED   1\n"
    "#define WARP_XCHG      \"xchg\"   /* same directory as exchange.sh uses */\n"
    "\n"
    "/* ---- steady or transient -------------------------------------------------\n"
    "   FLUENT calls Execute At End once per ITERATION when the solver is steady\n"
    "   and once per TIME STEP when it is transient.  Nothing in the API says\n"
    "   which one is happening, so it is set here, by warp2fluent -transient, and\n"
    "   the same value is baked into UDF_fluent2warp.c.\n"
    "     WARP_TRANSIENT 0  steady.  This file behaves exactly as it always has.\n"
    "     WARP_TRANSIENT 1  transient.  WARP_INTERVAL counts TIME STEPS, and the\n"
    "                       profile is installed by warp_wait_and_load below,\n"
    "                       which blocks.  Needs a COMPILED UDF: the wait uses\n"
    "                       nanosleep(), which the interpreter does not have. */\n"
    "#define WARP_TRANSIENT %d\n"
    "#define WARP_INTERVAL  %-3d  /* time steps per exchange; WARP_TRANSIENT only */\n"
    "#define WARP_POLL_US   20000     /* how often that wait looks, microseconds  */\n"
    "#define WARP_WAIT_MAX  600.0     /* s to wait for CFDWARP before giving up   */\n"
    "\n"
    "static int warp_xchg_n = 0;   /* exchange number this file is waiting for */\n"
    "\n"
    "#if WARP_TRANSIENT\n"
    "#include <time.h>\n"
    "/* The whole point of the transient path: give the core back while waiting,\n"
    "   instead of spinning on the directory.  nanosleep rather than usleep only\n"
    "   because it needs no useconds_t; both build warning-free under FLUENT's\n"
    "   Build (checked 2026-09-02, transient-26/probe3). */\n"
    "static void warp_nap(long us)\n"
    "{\n"
    "    struct timespec ts;\n"
    "    ts.tv_sec  = us/1000000L;\n"
    "    ts.tv_nsec = (us%%1000000L)*1000L;\n"
    "    nanosleep(&ts, NULL);\n"
    "}\n"
    "\n"
    "static int warp_there(const char *path)\n"
    "{\n"
    "    FILE *fp;\n"
    "    fp = fopen(path, \"rb\");\n"
    "    if (fp == NULL) return 0;\n"
    "    fclose(fp);\n"
    "    return 1;\n"
    "}\n"
    "#endif\n"
    "\n"
    "/* ---- transient: hold the time step until CFDWARP has answered ------------\n"
    "   In a transient run Execute At End fires once per TIME STEP, AFTER the step\n"
    "   has been solved, so a profile installed there is one whole time step late\n"
    "   and the step just computed used stale inflow data.  Worse, the look is not\n"
    "   a wait: FLUENT marches physical time whether or not CFDWARP has answered,\n"
    "   so the two codes' clocks drift apart silently.\n"
    "\n"
    "   Adjust is the hook that fixes both.  It is called at the START of every\n"
    "   inner iteration, with N_TIME and CURRENT_TIME already advanced to the step\n"
    "   about to be solved (measured, transient-26/probe2), so the first Adjust\n"
    "   call of time step m is the last moment at which the profile for step m can\n"
    "   still be installed.  This one blocks there until CFDWARP's message for\n"
    "   that step is on disk, which is what makes the exchange time accurate and\n"
    "   what stops FLUENT running ahead.\n"
    "\n"
    "   Timing over one window of WARP_INTERVAL time steps, I for short:\n"
    "     start of step m, (m-1)%%I == 0 : this function waits for w2f.<n>.ready\n"
    "                                     and loads it; steps m .. m+I-1 are then\n"
    "                                     solved on exchange n\n"
    "     end of step m+I-1, m%%I == 0   : warp_exchange in UDF_fluent2warp.c\n"
    "                                     publishes f2w.<n> and CFDWARP goes on\n"
    "   Both halves compute the window from N_TIME and WARP_INTERVAL alone, so\n"
    "   they need no shared symbol and cannot disagree about whose turn it is.\n"
    "\n"
    "   SELECT IT under Define -> User-Defined -> Function Hooks -> Adjust.  In a\n"
    "   transient run that selection is not optional: without it nothing installs\n"
    "   the profile, and warp_exchange says so and refuses to answer rather than\n"
    "   feeding CFDWARP a FLUENT state computed from a profile it never read. */\n"
    "DEFINE_ADJUST(warp_wait_and_load, warp_domain)\n"
    "{\n"
    "#if WARP_COUPLED && WARP_TRANSIENT\n"
    "    static int warp_last_time = -1;\n"
    "    static int warp_gaveup = 0;\n"
    "    char path[512], stop[512];\n"
    "    double waited = 0.0;\n"
    "\n"
    "    if (warp_gaveup) return;\n"
    "    /* Adjust runs once per inner iteration; everything below is once per step */\n"
    "    if ((int)N_TIME == warp_last_time) return;\n"
    "    if ((int)N_TIME < 1) return;        /* the call before the march starts */\n"
    "    warp_last_time = (int)N_TIME;\n"
    "    if (((int)N_TIME - 1) %% WARP_INTERVAL != 0) return;\n"
    "\n"
    "    sprintf(path, \"%%s/w2f.%%d.ready\", WARP_XCHG, warp_xchg_n);\n"
    "    sprintf(stop, \"%%s/STOP\", WARP_XCHG);\n"
    "    while (!warp_there(path)) {\n"
    "        if (warp_there(stop)) {\n"
    "            warp_gaveup = 1;\n"
    "            WARP_PRINT(\"\\nwarp2fluent: STOP seen while waiting for exchange %%d;\"\n"
    "                       \" CFDWARP has finished.\\n\", warp_xchg_n);\n"
    "            WARP_PRINT(\"             this zone keeps the last profile it was sent.\\n\");\n"
    "            return;\n"
    "        }\n"
    "        warp_nap(WARP_POLL_US);\n"
    "        waited = waited + WARP_POLL_US*1.0e-6;\n"
    "        if (waited > WARP_WAIT_MAX) {\n"
    "            warp_gaveup = 1;\n"
    "            WARP_PRINT(\"\\nwarp2fluent: gave up after %%.0f s waiting for %%s.\\n\",\n"
    "                       waited, path);\n"
    "            WARP_PRINT(\"             The CFDWARP side is not running, is not writing into\\n\");\n"
    "            WARP_PRINT(\"             %%s, or died. FLUENT keeps marching from here on the\\n\", WARP_XCHG);\n"
    "            WARP_PRINT(\"             last profile it holds, so THIS RUN IS NO LONGER COUPLED.\\n\");\n"
    "            return;\n"
    "        }\n"
    "    }\n"
    "    warp_load();\n"
    "    warp_thread = NULL;\n"
    "    warp_mapok  = 0;\n"
    "    /* Say so in the exchange directory, not just in C.  This file is what\n"
    "       UDF_fluent2warp.c waits for before it answers, and it is the only\n"
    "       thing that distinguishes \"the profile for this window was installed\"\n"
    "       from \"CFDWARP has published something\": the w2f.<n>.ready stamp\n"
    "       appears whether or not this hook is selected, so waiting on that one\n"
    "       would let a run with no Adjust hook answer from a time step solved on\n"
    "       stale data. Zero bytes, so it is whole the instant it exists and\n"
    "       several ranks creating it is harmless. */\n"
    "    sprintf(path, \"%%s/w2f.%%d.installed\", WARP_XCHG, warp_xchg_n);\n"
    "    { FILE *fp; fp = fopen(path, \"wb\"); if (fp) fclose(fp);\n"
    "      else WARP_PRINT(\"warp2fluent: cannot create %%s\\n\", path); }\n"
    "    warp_xchg_n = warp_xchg_n + 1;\n"
    "    WARP_PRINT(\"warp2fluent: exchange %%d loaded from %%s at time step %%d,\"\n"
    "               \" t=%%.6E s, waited %%.2f s\\n\",\n"
    "               warp_xchg_n - 1, WARP_FILE, (int)N_TIME, (double)CURRENT_TIME, waited);\n"
    "#else\n"
    "    /* Compiled either way, for the udf_names.c reason given below, and inert\n"
    "       in a steady run.  It says so once rather than in silence: a hook that\n"
    "       is selected and quietly does nothing is exactly how a coupling looks\n"
    "       alive while it is not. */\n"
    "    static int warp_said_adj = 0;\n"
    "    if (!warp_said_adj) {\n"
    "        warp_said_adj = 1;\n"
    "        WARP_PRINT(\"warp2fluent: warp_wait_and_load is inactive (WARP_TRANSIENT %%d,\"\n"
    "                   \" WARP_COUPLED %%d).\\n\", (int)WARP_TRANSIENT, (int)WARP_COUPLED);\n"
    "        WARP_PRINT(\"             In a steady run warp_reload_at_end installs the profile;\\n\");\n"
    "        WARP_PRINT(\"             regenerate with warp2fluent -transient for a transient one.\\n\");\n"
    "    }\n"
    "#endif\n"
    "}\n"
    "\n",
    TRANSIENT, TRANSIENT ? XINTERVAL : 1);
  fprintf(udf,
    "/* WHY THE #if IS INSIDE THIS FUNCTION AND NOT AROUND IT.  Build generates\n"
    "   libudf/src/udf_names.c by scanning these sources, as text and with no\n"
    "   preprocessing, for lines whose first word begins with a DEFINE_ prefix.\n"
    "   Every name found gets an extern declaration and an entry in udf_data[],\n"
    "   so one that sits in the dead half of an #if becomes a symbol that is\n"
    "   declared and never defined.  The library still links, a shared object\n"
    "   being allowed to carry undefined symbols, and Load then fails with\n"
    "       Error: undefined symbol: <name>\n"
    "   after which NOTHING in the library is available, not even the functions\n"
    "   that did compile.  So every DEFINE_ here is unconditional and the\n"
    "   switches live in the bodies. */\n"
    "DEFINE_EXECUTE_AT_END(warp_reload_at_end)\n"
    "{\n"
    "#if WARP_COUPLED && !WARP_TRANSIENT\n"
    "    FILE *fp;\n"
    "    char path[512];\n"
    "\n"
    "    sprintf(path, \"%%s/w2f.%%d.ready\", WARP_XCHG, warp_xchg_n);\n"
    "    fp = fopen(path, \"rb\");        /* fopen rather than stat or access: those */\n"
    "    if (fp == NULL) return;        /* two are not reliably in the interpreter */\n"
    "    fclose(fp);\n"
    "\n"
    "    warp_load();\n"
    "    warp_thread = NULL;\n"
    "    warp_mapok  = 0;\n"
    "    warp_xchg_n = warp_xchg_n + 1;\n"
    "    WARP_PRINT(\"warp2fluent: exchange %%d loaded from %%s\\n\",\n"
    "               warp_xchg_n - 1, WARP_FILE);\n"
    "#else\n"
    "    /* Inactive, and it says so once rather than returning silently, since a\n"
    "       hook that is selected under Function Hooks and quietly does nothing is\n"
    "       exactly how a coupling looks alive while serving one profile for\n"
    "       ever.  Two ways to get here:\n"
    "         WARP_COUPLED 0    no CFDWARP on the other end, nothing to watch for\n"
    "         WARP_TRANSIENT 1  a transient run installs the profile from\n"
    "                           warp_wait_and_load, an Adjust hook, because this\n"
    "                           one runs after the time step has been solved.\n"
    "                           Leaving this hook selected as well is harmless --\n"
    "                           it is this branch -- but the two would otherwise\n"
    "                           share warp_xchg_n and each consume the other's\n"
    "                           exchange, so it is switched off here rather than\n"
    "                           left to the user to deselect. */\n"
    "    static int warp_said_reload = 0;\n"
    "    if (!warp_said_reload) {\n"
    "        warp_said_reload = 1;\n"
    "#if WARP_TRANSIENT\n"
    "        WARP_PRINT(\"warp2fluent: warp_reload_at_end is inactive, WARP_TRANSIENT is 1.\\n\");\n"
    "        WARP_PRINT(\"             Select warp_wait_and_load under Function Hooks -> Adjust;\\n\");\n"
    "        WARP_PRINT(\"             that is what installs the profile in a transient run.\\n\");\n"
    "#else\n"
    "        WARP_PRINT(\"warp2fluent: warp_reload_at_end is inactive, WARP_COUPLED is 0.\\n\");\n"
    "        WARP_PRINT(\"             %%s is then re-read only by warp_load_data.\\n\", WARP_FILE);\n"
    "#endif\n"
    "    }\n"
    "#endif\n"
    "}\n\n");

  /* ---------------- one profile per variable ---------------- */
  for (i = 0; i < numvars; i++) {
    fprintf(udf,
      "/* %s */\n"
      "DEFINE_PROFILE(%s_profile, thread, position)\n"
      "{\n"
      "    face_t f;\n"
      "    int idx;\n"
      "\n"
      "    warp_prepare(thread);\n"
      "    begin_f_loop(f, thread)\n"
      "    {\n"
      "        idx = warp_index(thread, f);\n"
      "        F_PROFILE(f, thread, position) = (idx >= 0 ? warp_%s[idx] : 0.0);\n"
      "    }\n"
      "    end_f_loop(f, thread)\n"
      "}\n\n",
      varnames[i], cnames[i], cnames[i]);
  }

  fclose(udf);
}


/* -tent: drop the 0.0003 weight floor, see find_interpolation_weight() */
static int TENTWEIGHT=0;


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
    if (TENTWEIGHT) {
      *thisweight=max(0.0,1.0-(*thisweight));
    } else {
      *thisweight=max(1e-16,max(0.0003-(*thisweight)*0.00001,1.0-(*thisweight)));
    }
    EXM_free_matrix(&mat1);
    EXM_free_matrix(&mat1inv);
    EXM_free_matrix(&mat2);
    EXM_free_matrix(&mat3);
  }
}


/* Write the FLUENT-side UDF that produces a CFDWARP interpolation file */
/* Mean distance from a point of the set to its nearest other point */
static double mesh_mean_spacing(const double *x, long n, long nd) {
  long l, m, dim;
  double d2, d2min, sum;
  if (n < 2) return 0.0;
  sum = 0.0;
  for (l = 0; l < n; l++) {
    d2min = -1.0;
    for (m = 0; m < n; m++) {
      if (m == l) continue;
      d2 = 0.0;
      for (dim = 0; dim < nd; dim++) d2 += sqr(x[l*nd+dim] - x[m*nd+dim]);
      if (d2min < 0.0 || d2 < d2min) d2min = d2;
    }
    sum += sqrt(d2min);
  }
  return sum/(double)n;
}


static void write_fluent2warp_udf_file(const char *udf_filename, const char *vars_fluid_str,
                                       char **varnames, long numvars, long numdim, long numspec,
                                       double xstation, double *dx1mean, double *dx2mean,
                                       double *dx3mean,
                                       long nf, long windowis, long windowie, long iter,
                                       double effiter_U, double effiter_R, double CFL,
                                       double dxwarp, double dxfluent, long numfaces,
                                       int CLAMP, double pfloor, double tfloor,
                                       double pmean) {
  FILE *udf;
  long cnt, spec, reclen, fidx;
  double hmax, dxscale;
  char pscalestr[32];
  char flname[64];
  char editnote[900];
  long nleft;

  /* validate -fluentspecies against the CFDWARP species list */
  if (FLSPEC[0] != '\0') {
    for (cnt = 0; cnt < numvars; cnt++) {
      if (strncmp(varnames[cnt], "w_", 2) != 0) continue;
      if (fluent_species_index(varnames[cnt], NULL, 0) < 0) {
        fprintf(stderr,
          "\n\n-fluentspecies does not name the CFDWARP species %s.\n"
          "The list given was \"%s\". It must contain every CFDWARP species, spelled\n"
          "as FLUENT spells it, in the order your FLUENT mixture material lists them\n"
          "(bulk species last). Nothing was written.\n\n",
          varnames[cnt], FLSPEC);
        exit (EXIT_FAILURE);
      }
    }
  }
  nleft = 0;
  strcpy(editnote, "   WHAT YOU STILL HAVE TO EDIT: ");
  if (!ZONEID_GIVEN)    { strcat(editnote, "WARP_ZONE_ID"); nleft++; }
  if (!XINTERVAL_GIVEN) { strcat(editnote, nleft ? ", WARP_INTERVAL" : "WARP_INTERVAL"); nleft++; }
  if (FLSPEC[0] == '\0') {
    strcat(editnote, nleft ? ",\n   and the FLUENT species indices in warp_fill_record()"
                           : "the FLUENT species indices in\n   warp_fill_record()");
    nleft++;
  }
  if (nleft == 0)
    strcpy(editnote,
      "   NOTHING IS LEFT TO EDIT IN THIS FILE.  The zone ID, the exchange interval\n"
      "   and the CFDWARP-to-FLUENT species mapping -- the three things this tool\n"
      "   cannot read out of the CFDWARP interpolation file -- were all given on the\n"
      "   warp2fluent command line and are baked in below.  Regenerating with the\n"
      "   same flags reproduces this file, so there is no hand edit to lose.\n");
  else
    strcat(editnote,
      ".\n   Pass -zone, -interval and -fluentspecies to warp2fluent instead and it\n"
      "   will bake them in, which survives the next regeneration; a hand edit does\n"
      "   not, because this whole block is rewritten every time the tool runs.\n");


  if (pmean > 0.0) sprintf(pscalestr, "%.3E", pmean);
  else             strcpy(pscalestr, "the interface");
  udf = fopen(udf_filename, "w");
  if (!udf) {
    fprintf(stderr, "Could not open %s for writing\n", udf_filename);
    return;
  }
  reclen = numvars + (numdim + 1) * numdim; /* vars, then x, dx1, dx2 [, dx3] */
  hmax = 0.0;
  for (cnt = 0; cnt < numdim; cnt++)
    hmax = max(hmax, fabs(dx1mean[cnt]) + fabs(dx2mean[cnt]) + fabs(dx3mean[cnt]));
  dxscale = 1.0;
  if (hmax > 0.0 && dxwarp > 0.0) dxscale = 0.5*dxwarp/hmax;

  fprintf(udf,
    "/* -------------------------------------------------------------------------\n"
    "   FLUENT UDF generated by CFDWARP tools/warp2fluent.c -- do not edit by hand,\n"
    "   it is rewritten every time warp2fluent runs.\n"
    "\n"
    "   It writes a CFDWARP interpolation file from a FLUENT boundary zone, to be\n"
    "   read back on the CFDWARP side with\n"
    "       ReadInterpolationFileZone(\"%s\", is,js, ie,je);\n"
    "   inside Cycle(). The record layout and the vars_fluid string below were taken\n"
    "   from the CFDWARP interpolation file this tool just read, so they match that\n"
    "   run by construction. Regenerate this file whenever the CFDWARP species list,\n"
    "   the number of dimensions or the fluid module changes.\n"
    "\n"
    "%s"
    "\n"
    "   HOW TO LOAD IT.  Interpreted or compiled, either works, and this file no\n"
    "   longer resolves anything against UDF.c, so the two can even be loaded in\n"
    "   the other order.  Under Function Hooks -> Execute At End select BOTH\n"
    "   warp_exchange (this file) and warp_reload_at_end (UDF.c); the list takes\n"
    "   several entries and each half of the exchange needs its own.\n"
    "     interpreted, serial : leave WARP_PRINT at printf and WARP_USE_RENAME at 0\n"
    "     compiled, or any parallel run : WARP_PRINT Message, WARP_USE_RENAME 1\n"
    "   A parallel run has to be compiled whatever these say: the PRF_CSEND and\n"
    "   PRF_CRECV block below is not something the interpreter accepts.\n"
    "%s"
    "   ------------------------------------------------------------------------- */\n"
    "\n"
    "#include \"udf.h\"\n"
    "\n"
    "/* ===================== USER SETTINGS ===================================== */\n"
    "#define WARP_ZONE_ID   %-3d                /* FLUENT zone ID of the interface */\n"
    "#define WARP_INTERVAL  %-3d                /* FLUENT %s per exchange */\n"
    "#define WARP_FILENAME  \"%s\"\n"
    "\n"
    "/* Steady or transient FLUENT run.  FLUENT calls DEFINE_EXECUTE_AT_END once\n"
    "   per ITERATION when the solver is steady and once per TIME STEP when it is\n"
    "   transient, and there is nothing in the API that says which is happening,\n"
    "   so it is set here by warp2fluent -transient and the same value is baked\n"
    "   into UDF.c.  At 0 in a transient run WARP_INTERVAL silently becomes a\n"
    "   number of time steps divided by the inner iteration count, and FLUENT\n"
    "   marches physical time whether or not CFDWARP has answered.\n"
    "     WARP_TRANSIENT 0  count iterations, look once per iteration, never wait\n"
    "     WARP_TRANSIENT 1  count time steps; the profile is installed by\n"
    "                       warp_wait_and_load in UDF.c, an Adjust hook that\n"
    "                       BLOCKS at the start of the time step, so FLUENT and\n"
    "                       CFDWARP advance one window at a time together.\n"
    "                       Select that hook under Function Hooks -> Adjust; this\n"
    "                       file refuses to answer without it. */\n"
    "#define WARP_TRANSIENT %d\n"
    "\n"
    "/* Rendezvous with a CFDWARP run going on at the same time.\n"
    "     WARP_COUPLED 0  one-shot: write WARP_FILENAME every WARP_INTERVAL\n"
    "                     %s and never wait for anybody. This is what this\n"
    "                     file did before the rendezvous was added.\n"
    "     WARP_COUPLED 1  take turns with CFDWARP through the directory below.\n"
    "   WARP_XCHG must name the same directory the CFDWARP-side exchange.sh uses,\n"
    "   as seen from the directory FLUENT was started from. */\n"
    "#define WARP_COUPLED   1\n"
    "#define WARP_XCHG      \"xchg\"\n"
    "\n"
    "/* Where the messages of this file go, exactly as in UDF.c. Set to one of:\n"
    "       printf      interpreted UDF, serial              (default)\n"
    "       Message     compiled UDF, and any parallel run\n"
    "       warp_noop   silent, if the interpreter refuses both\n"
    "   Nothing else in this file depends on the choice. */\n"
    "#define WARP_PRINT printf\n"
    "\n"
    "/* 1 publishes every file by writing a temporary name and rename()ing it, so\n"
    "   the payload appears whole or not at all.  0 writes each file under its final\n"
    "   name instead, for the interpreter, whose C library does not include\n"
    "   rename().  The default is 1, which is right for a compiled UDF and for any\n"
    "   parallel run; set it to 0 only for an interpreted serial run, together with\n"
    "   WARP_PRINT printf above.\n"
    "   At 0 the payload grows in place under the name CFDWARP will open, so a peer\n"
    "   that looks at the directory rather than waiting on the stamp -- or a stamp\n"
    "   published early, which is what the rank guard in warp_collect_and_write()\n"
    "   below exists to prevent -- can see a half-written file.  The zero-byte stamp\n"
    "   needs no rename of its own: it is whole the instant it exists.  The names\n"
    "   carry the exchange number and are never reused, so there is no stale file to\n"
    "   mistake for a fresh one. */\n"
    "#define WARP_USE_RENAME 1\n"
    "/* ========================================================================= */\n"
    "\n"
    "/* Old-style empty parameter list on purpose: it takes whatever WARP_PRINT is\n"
    "   called with and discards it. Unused unless WARP_PRINT is set to it. A name\n"
    "   of its own so that this file and UDF.c can be compiled together. */\n"
    "void warp_noop_f2w() { }\n"
    "\n",
    "fluent2warp.interpolation", editnote,
    TRANSIENT ?
      "   THIS IS THE TRANSIENT BUILD (warp2fluent -transient).  A third hook has\n"
      "   to be selected as well: warp_wait_and_load from UDF.c, under Function\n"
      "   Hooks -> Adjust.  In a transient run Execute At End fires once per TIME\n"
      "   STEP, after the step has been solved, so it is too late to install an\n"
      "   inflow profile there; the Adjust hook installs it at the start of the\n"
      "   step and blocks until CFDWARP has sent it.  Both generated files must\n"
      "   be COMPILED, not interpreted: that wait uses nanosleep().\n"
      : "",
    ZONEID, XINTERVAL,
    TRANSIENT ? "time steps" : "iterations",
    "fluent2warp.interpolation", TRANSIENT,
    TRANSIENT ? "time steps" : "iterations");

  fprintf(udf,
    "#define WARP_NUMVARS   %ld    /* doubles per record, variables part          */\n"
    "#define WARP_ND        %ld    /* CFDWARP number of dimensions                */\n"
    "#define WARP_NS        %ld    /* CFDWARP number of species                   */\n"
    "#define WARP_RECLEN    %ld    /* WARP_NUMVARS + (WARP_ND+1)*WARP_ND          */\n"
    "#define WARP_MAXFACES  %ld    /* faces in the -m mesh, x4 headroom           */\n"
    "\n"
    "/* x of the shared interface, taken from the CFDWARP interpolation file. */\n"
    "#define WARP_XSTATION  %.8E\n"
    "\n"
    "/* dx1 and dx2 are the support of the weight CFDWARP gives this file's points.\n"
    "   read_data_file_interpolation_zone() forms w = max(1e-16, max(3e-4 - 1e-5 s,\n"
    "   1 - s)) with s the offset from the node measured in units of [dx1 dx2], so\n"
    "   the support is a tent that reaches zero one dx away, and the cutoff radius\n"
    "   sqrt(1.1 sum_dim (|dx1|+|dx2|)^2) only decides which points are looked at.\n"
    "   The tent has to be wide enough that every CFDWARP node sees at least one\n"
    "   FLUENT face -- half the FLUENT spacing, %.3E m here -- and narrow enough\n"
    "   that it does not reach the neighbouring CFDWARP node, %.3E m away, or the\n"
    "   read is a smoother rather than an identity and the CFDWARP residual never\n"
    "   settles.  The values below are the CFDWARP cell size scaled to sit in the\n"
    "   middle of that window; the raw cell size is the scale 1.0 case, and it is\n"
    "   already slightly too wide.%s */\n"
    "#define WARP_DXSCALE   %.8E   /* applied below, kept visible on purpose */\n",
    numvars, numdim, numspec, reclen, 4*numfaces+16, xstation, 0.5*dxfluent, dxwarp,
    (0.5*dxfluent < 0.5*dxwarp) ? "" :
      "\n   WARNING: the FLUENT interface mesh is not finer than the CFDWARP one, so\n"
      "   no dx satisfies both ends of that window.  The read is unavoidably a\n"
      "   smoother: judge convergence on the interface movement per exchange\n"
      "   rather than on the CFDWARP residual.", dxscale);
  for (cnt = 0; cnt < numdim; cnt++)
    fprintf(udf, "#define WARP_DX1_%ld    %.8E   /* %.8E * WARP_DXSCALE */\n",
            cnt, dx1mean[cnt]*dxscale, dx1mean[cnt]);
  for (cnt = 0; cnt < numdim; cnt++)
    fprintf(udf, "#define WARP_DX2_%ld    %.8E   /* %.8E * WARP_DXSCALE */\n",
            cnt, dx2mean[cnt]*dxscale, dx2mean[cnt]);
  /* dx3 exists only in 3D.  Left at zero -- which is what this tool used to
     write into every record -- it makes the matrix [dx1 dx2 dx3] that
     find_interpolation_weight() inverts singular, and the CFDWARP side stops
     on a zero determinant at the first zone read. */
  if (numdim > 2)
    for (cnt = 0; cnt < numdim; cnt++)
      fprintf(udf, "#define WARP_DX3_%ld    %.8E   /* %.8E * WARP_DXSCALE */\n",
              cnt, dx3mean[cnt]*dxscale, dx3mean[cnt]);

  /* Unit vectors of the two support directions */
  {
    double n1, n2;
    n1 = 0.0; n2 = 0.0;
    for (cnt = 0; cnt < numdim; cnt++) { n1 += sqr(dx1mean[cnt]); n2 += sqr(dx2mean[cnt]); }
    n1 = sqrt(n1); n2 = sqrt(n2);
    for (cnt = 0; cnt < numdim; cnt++)
      fprintf(udf, "#define WARP_DX1U_%ld   %.8E\n", cnt,
              (n1 > 0.0) ? dx1mean[cnt]/n1 : (cnt == 0 ? 1.0 : 0.0));
    for (cnt = 0; cnt < numdim; cnt++)
      fprintf(udf, "#define WARP_DX2U_%ld   %.8E\n", cnt,
              (n2 > 0.0) ? dx2mean[cnt]/n2 : (cnt == 1 ? 1.0 : 0.0));
    if (numdim > 2) {
      double n3 = 0.0;
      for (cnt = 0; cnt < numdim; cnt++) n3 += sqr(dx3mean[cnt]);
      n3 = sqrt(n3);
      for (cnt = 0; cnt < numdim; cnt++)
        fprintf(udf, "#define WARP_DX3U_%ld   %.8E\n", cnt,
                (n3 > 0.0) ? dx3mean[cnt]/n3 : (cnt == 2 ? 1.0 : 0.0));
    }
    /* The three lists below are what the warp_dxNu[] arrays are built from, so
       that one declaration serves both 2D and 3D. */
    if (numdim > 2) {
      fprintf(udf, "#define WARP_DX1U      WARP_DX1U_0, WARP_DX1U_1, WARP_DX1U_2\n");
      fprintf(udf, "#define WARP_DX2U      WARP_DX2U_0, WARP_DX2U_1, WARP_DX2U_2\n");
      fprintf(udf, "#define WARP_DX3U      WARP_DX3U_0, WARP_DX3U_1, WARP_DX3U_2\n");
      fprintf(udf, "#define WARP_DX2MAG    sqrt(WARP_DX2_0*WARP_DX2_0 + WARP_DX2_1*WARP_DX2_1"
                   " + WARP_DX2_2*WARP_DX2_2)\n");
    } else {
      fprintf(udf, "#define WARP_DX1U      WARP_DX1U_0, WARP_DX1U_1\n");
      fprintf(udf, "#define WARP_DX2U      WARP_DX2U_0, WARP_DX2U_1\n");
      fprintf(udf, "#define WARP_DX2MAG    sqrt(WARP_DX2_0*WARP_DX2_0 + WARP_DX2_1*WARP_DX2_1)\n");
    }
  }
  fprintf(udf,
    "\n/* WARP_LOCALDX 1: the constants above are NOT what is sent.  Every face gets\n"
    "   its own dx1 and dx2, scaled to the local spacing of the interface mesh by\n"
    "   warp_set_local_support() below; they coincide with the constants only on a\n"
    "   uniform interface.  One pair taken from the MEAN spacing cannot serve a mesh\n"
    "   clustered at the wall: it is several times too wide at the wall, where the\n"
    "   read then averages a dozen faces, and several times too narrow in the core,\n"
    "   where read_data_file_interpolation_zone() finds no face inside the cutoff\n"
    "   radius and leaves the CFDWARP node UNTOUCHED -- silently, because a node that\n"
    "   accumulated zero weight is simply skipped, with no message on either side.\n"
    "   Set this to 0, or regenerate with warp2fluent -globaldx, to get the old\n"
    "   single-pair behaviour back for comparison. */\n"
    "#define WARP_LOCALDX   %d\n", GLOBALDX ? 0 : 1);

  /* Same stamp UDF.c carries */
  {
    char gendir[PATH_MAX + 1];
    char gentime[64];
    time_t now;
    struct tm *lt;
    if (getcwd(gendir, PATH_MAX) == NULL) strcpy(gendir, "unknown");
    now = time(NULL);
    lt = localtime(&now);
    if (lt == NULL || strftime(gentime, sizeof(gentime), "%Y-%m-%d %H:%M:%S", lt) == 0)
      strcpy(gentime, "unknown");
    fprintf(udf, "\n/* Where and when this file was generated.  Printed once, by the first\n"
                 "   exchange it writes, so that the FLUENT console names the UDF actually\n"
                 "   running rather than the one you last regenerated. */\n");
    fprintf(udf, "#define WARP_F2W_GEN_DIR  \"%s\"\n", gendir);
    fprintf(udf, "#define WARP_F2W_GEN_TIME \"%s\"\n", gentime);
  }

  /* index of P, of T and of the first species inside a record, needed by the
     validation code emitted below.  -1 means the variable is not in vars_fluid. */
  {
    long iP = -1, iT = -1, iw = -1, nw = 0;
    for (cnt = 0; cnt < numvars; cnt++) {
      if (strcmp(varnames[cnt], "P") == 0) iP = cnt;
      else if (strcmp(varnames[cnt], "T") == 0) iT = cnt;
      else if (strncmp(varnames[cnt], "w_", 2) == 0) { if (iw < 0) iw = cnt; nw++; }
    }
    fprintf(udf,
      "\n"
      "/* ---- what a usable record has to look like ------------------------------\n"
      "   A record whose species mass fractions are all zero is the single most\n"
      "   expensive thing this file can hand CFDWARP, because the complaint comes\n"
      "   from somewhere else entirely.  read_data_file_interpolation_zone() runs\n"
      "   with REFORMAT_INITVAR_SPECIES_SUM_CHECK switched off, so src/init.c does\n"
      "   not object to fractions that do not sum to 1 -- it divides each of them by\n"
      "   their sum.  With every fraction zero that is 0/0, the node is initialised\n"
      "   with NaN, and the run dies a few lines later inside the thermodynamics:\n"
      "       Couldn't find root for temperature in _T_from_w_h subroutine.\n"
      "       T=-NAN  dT=-NAN\n"
      "       CFDWARP fatal error. Exiting.\n"
      "   naming neither the file, nor the boundary, nor the species.\n"
      "\n"
      "   Flooring each fraction at a tiny positive number, which is what this file\n"
      "   did up to 2026-08-25, is not a guard: %ld fractions all equal to 1e-30 sum\n"
      "   to %ld"
      "e-30 and normalise to 1/%ld each, so CFDWARP silently receives an\n"
      "   invented equimolar mixture instead of an error.  A bad record is refused\n"
      "   here instead, at the only place that still knows which face it came from.\n"
      "\n"
      "   Rejected: any non-finite value, T below WARP_TMIN, absolute P below\n"
      "   WARP_PMIN, species fractions summing to less than WARP_WSUMMIN.\n"
      "   Warned about but sent: a species sum that is positive but not 1 -- that is\n"
      "   a FLUENT mixture with fewer species than CFDWARP, and it is renormalised\n"
      "   here so that what CFDWARP stores is what this file measured. */\n"
      "#define WARP_TMIN      1.0e-03   /* K,  below this it is not a temperature */\n"
      "#define WARP_PMIN      1.0e-03   /* Pa, absolute, after adding warp_Pop    */\n"
      "#define WARP_WSUMMIN   1.0e-06   /* below this the fractions carry nothing */\n"
      "#define WARP_WSUMTOL   1.0e-03   /* |sum-1| tolerated without a warning    */\n"
      "#define WARP_BIG       1.0e+300  /* NaN compares false against both ends   */\n"
      "\n"
      "/* ---- clamping instead of refusing ---------------------------------------\n"
      "   With WARP_CLAMP 1 a P or T below its limit is raised to the floor below\n"
      "   and SENT, instead of stopping the exchange.  It is off by default and it\n"
      "   should stay off until you know WHY the value is low, because a floor does\n"
      "   not repair the record -- it only stops this file objecting to it.  What\n"
      "   CFDWARP then stores at that node is the floor, as a measurement.\n"
      "\n"
      "   Pick the floor at the scale of the flow, not at the scale of the check.\n"
      "   WARP_PMIN is 1e-3 Pa because that is where a number stops being a\n"
      "   pressure at all; a floor of 1e-3 Pa on an interface running near %s Pa\n"
      "   punches a seven-decade hole in the profile and CFDWARP believes it.  The\n"
      "   default below is instead the smallest absolute pressure in the CFDWARP\n"
      "   interpolation file this UDF was generated from, so a clamped face at\n"
      "   least lands in the range the rest of the interface occupies.  Override\n"
      "   both with -pfloor and -tfloor, or edit them here.\n"
      "\n"
      "   A clamped sweep says so on EVERY exchange, not once.  That is deliberate:\n"
      "   the count is the only record that the profile CFDWARP converged on was\n"
      "   partly invented here.\n"
      "\n"
      "   Species get no floor from WARP_CLAMP, because there is no honest one.  A\n"
      "   tiny value on every fraction is precisely the equimolar invention\n"
      "   described above: %ld fractions at 1e-30 renormalise to 1/%ld each, and the\n"
      "   record then looks perfectly well formed.  To substitute a composition you\n"
      "   have to name it: set WARP_WCLAMP 1 and put the mass fractions in\n"
      "   WARP_WFALLBACK, in the vars_fluid order listed just below.  Left as it is,\n"
      "   a face carrying no species is still refused. */\n"
      "#define WARP_CLAMP     %d         /* 1: floor P and T instead of refusing   */\n"
      "#define WARP_PFLOOR    %.8E   /* Pa, absolute. Used only if WARP_CLAMP  */\n"
      "#define WARP_TFLOOR    %.8E   /* K.  Used only if WARP_CLAMP            */\n"
      "#define WARP_WCLAMP    0         /* 1: use WARP_WFALLBACK below            */\n"
      "/* #define WARP_WFALLBACK 0.0, 0.0, 0.0, 0.233, 0.767 */\n"
      "\n"
      "#define WARP_IP        %ld    /* index of P in a record, -1 if absent  */\n"
      "#define WARP_IT        %ld    /* index of T                            */\n"
      "#define WARP_IW        %ld    /* index of the first species fraction   */\n"
      "#define WARP_NW        %ld    /* how many species fractions there are  */\n"
      "\n"
      "static double warp_buf[WARP_MAXFACES * WARP_RECLEN];\n"
      "static int    warp_nrec = 0;\n"
      "static int    warp_nbad = 0;    /* records refused in this sweep          */\n"
      "static int    warp_nwarn = 0;   /* records renormalised in this sweep     */\n"
      "static char   warp_why[256];    /* why the first refused one was refused  */\n"
      "static int    warp_stamped = 0; /* generation banner printed once yet?     */\n"
      "static double warp_Pop  = 0.0;   /* FLUENT operating pressure, added to F_P */\n"
      "static int    warp_nclamp = 0;   /* records floored in this sweep           */\n"
      "static char   warp_clamped[256]; /* what the first floored one looked like  */\n"
      "static double warp_pmin_seen = 0.0; /* lowest absolute P this sweep saw     */\n"
      "static int    warp_nzerop = 0;   /* faces whose absolute P was exactly 0    */\n"
      "#if WARP_LOCALDX\n"
      "static const double warp_dx1u[WARP_ND] = { WARP_DX1U };\n"
      "static const double warp_dx2u[WARP_ND] = { WARP_DX2U };\n"
      "#if WARP_ND > 2\n"
      "static const double warp_dx3u[WARP_ND] = { WARP_DX3U };\n"
      "#endif\n"
      "static int    warp_dx_reported = 0;  /* the support banner, printed once    */\n"
      "#endif\n"
      "\n"
      "/* isnan() and isfinite() are not reliably available to the interpreter.  NaN\n"
      "   compares false against everything, so this one test catches NaN and both\n"
      "   infinities. */\n"
      "static int warp_finite(double v)\n"
      "{\n"
      "  if (v > -WARP_BIG && v < WARP_BIG) return 1;\n"
      "  return 0;\n"
      "}\n"
      "\n"
      "#if WARP_CLAMP || WARP_WCLAMP\n"
      "/* Remember the first floored face of the sweep and count the rest.  Reached\n"
      "   only from the WARP_CLAMP and WARP_WCLAMP branches, so it is compiled with\n"
      "   them: gcc -Wall does warn about an unused static function, and leaving\n"
      "   this one unconditional put a fresh warning in every default build.  An\n"
      "   #if is safe here for the same reason it is safe around warp_exists() --\n"
      "   this is a plain static, not a DEFINE_ hook, so Build's text scan of this\n"
      "   file never sees it and cannot declare a name the compiler then omits. */\n"
      "static void warp_note_clamp(const char *what, double was, double now, double y)\n"
      "{\n"
      "  warp_nclamp = warp_nclamp + 1;\n"
      "  if (warp_nclamp == 1)\n"
      "    sprintf(warp_clamped, \"%%s = %%E raised to %%E at y=%%E\", what, was, now, y);\n"
      "}\n"
      "#endif\n"
      "\n", (long)nw, (long)nw, (long)nw,
           pscalestr, (long)nw, (long)nw, CLAMP, pfloor, tfloor,
           iP, iT, iw, nw);
  }
  fprintf(udf,
    "/* CFDWARP reads the variables in exactly this order, and takes the names from\n"
    "   the vars_fluid string in the header. Do not reorder either one alone.\n");
  for (cnt = 0; cnt < numvars; cnt++)
    fprintf(udf, "     %2ld  %s\n", cnt, varnames[cnt]);
  fprintf(udf, "*/\n\n");

  fprintf(udf,
    "/* Returns 1 if the record is usable, 0 if it must not be sent.  On 0 it has\n"
    "   put the reason in warp_why. */\n"
    "static int warp_fill_record(double *r, face_t f, Thread *t)\n"
    "{\n"
    "  real xc[ND_ND];\n"
    "  double wsum;\n"
    "  int k;\n"
    "  F_CENTROID(xc, f, t);\n"
    "\n");
  spec = 0;
  for (cnt = 0; cnt < numvars; cnt++) {
    if (strcmp(varnames[cnt], "V[0]") == 0)
      fprintf(udf, "  r[%2ld] = F_U(f, t);                  /* %s */\n", cnt, varnames[cnt]);
    else if (strcmp(varnames[cnt], "V[1]") == 0)
      fprintf(udf, "  r[%2ld] = F_V(f, t);                  /* %s */\n", cnt, varnames[cnt]);
    else if (strcmp(varnames[cnt], "V[2]") == 0)
      fprintf(udf, "  r[%2ld] = F_W(f, t);                  /* %s */\n", cnt, varnames[cnt]);
    else if (strcmp(varnames[cnt], "P") == 0)
      fprintf(udf, "  r[%2ld] = F_P(f, t) + warp_Pop;       /* P: F_P is gauge pressure, CFDWARP wants absolute */\n", cnt);
    else if (strcmp(varnames[cnt], "T") == 0)
      fprintf(udf, "  r[%2ld] = F_T(f, t);                  /* T */\n", cnt);
    else if (strncmp(varnames[cnt], "w_", 2) == 0) {
      fidx = fluent_species_index(varnames[cnt], flname, (long)sizeof(flname));
      if (fidx >= 0)
        fprintf(udf, "  r[%2ld] = F_YI(f, t, %ld);            /* %-4s <- FLUENT \"%s\" */\n",
                cnt, fidx, varnames[cnt], flname);
      else
        fprintf(udf, "  r[%2ld] = F_YI(f, t, %ld);            /* %s: CHECK the FLUENT species index */\n",
                cnt, spec, varnames[cnt]);
      spec++;
    }
    else if (strcmp(varnames[cnt], "Tv") == 0)
      fprintf(udf, "  r[%2ld] = F_T(f, t);                  /* Tv: FLUENT has no vibrational temperature, sending T */\n", cnt);
    else if (strcmp(varnames[cnt], "Te") == 0)
      fprintf(udf, "  r[%2ld] = F_T(f, t);                  /* Te: FLUENT has no electron temperature, sending T */\n", cnt);
    else
      fprintf(udf, "  r[%2ld] = 0.0;                        /* %s: NO FLUENT COUNTERPART, edit or drop it from vars_fluid */\n",
              cnt, varnames[cnt]);
  }
  fprintf(udf,
    "\n"
    "  /* x of the file point: the CFDWARP interface station, not the FLUENT one,\n"
    "     so that the point lands on the CFDWARP boundary node. */\n"
    "  r[WARP_NUMVARS + 0] = WARP_XSTATION;\n");
  for (cnt = 1; cnt < numdim; cnt++)
    fprintf(udf, "  r[WARP_NUMVARS + %ld] = (double)xc[%ld];\n", cnt, cnt);
  for (cnt = 0; cnt < numdim; cnt++)
    fprintf(udf, "  r[WARP_NUMVARS + %ld*WARP_ND + %ld] = WARP_DX1_%ld;\n", (long)1, cnt, cnt);
  for (cnt = 0; cnt < numdim; cnt++)
    fprintf(udf, "  r[WARP_NUMVARS + %ld*WARP_ND + %ld] = WARP_DX2_%ld;\n", (long)2, cnt, cnt);
  if (numdim > 2) {
    for (cnt = 0; cnt < numdim; cnt++)
      fprintf(udf, "  r[WARP_NUMVARS + 3*WARP_ND + %ld] = WARP_DX3_%ld;\n", cnt, cnt);
  }
  fprintf(udf,
    "\n"
    "  /* ---- refuse anything CFDWARP cannot turn into a thermodynamic state ---- */\n"
    "  for (k = 0; k < WARP_NUMVARS; k++) {\n"
    "    if (!warp_finite(r[k])) {\n"
    "      sprintf(warp_why,\n"
    "              \"variable %%d of the record is not a finite number (%%E) at y=%%E\",\n"
    "              k, r[k], (double)xc[1]);\n"
    "      return 0;\n"
    "    }\n"
    "  }\n"
    "#if WARP_IT >= 0\n"
    "  if (r[WARP_IT] < WARP_TMIN) {\n"
    "#if WARP_CLAMP\n"
    "    warp_note_clamp(\"T\", r[WARP_IT], WARP_TFLOOR, (double)xc[1]);\n"
    "    r[WARP_IT] = WARP_TFLOOR;\n"
    "#else\n"
    "    sprintf(warp_why, \"T = %%E K at y=%%E, which is not a temperature\",\n"
    "            r[WARP_IT], (double)xc[1]);\n"
    "    return 0;\n"
    "#endif\n"
    "  }\n"
    "#endif\n"
    "#if WARP_IP >= 0\n"
    "  /* The lowest absolute pressure this sweep saw, floored or not.  Reported\n"
    "     every sweep: it is what says whether the interface is merely low or is\n"
    "     sitting at exactly zero, and those two have different causes. */\n"
    "  if (warp_nrec == 0 || r[WARP_IP] < warp_pmin_seen) warp_pmin_seen = r[WARP_IP];\n"
    "  if (r[WARP_IP] == 0.0) warp_nzerop = warp_nzerop + 1;\n"
    "  if (r[WARP_IP] < WARP_PMIN) {\n"
    "#if WARP_CLAMP\n"
    "    warp_note_clamp(\"P\", r[WARP_IP], WARP_PFLOOR, (double)xc[1]);\n"
    "    r[WARP_IP] = WARP_PFLOOR;\n"
    "#else\n"
    "    sprintf(warp_why,\n"
    "            \"absolute P = %%E Pa at y=%%E (F_P gauge %%E + operating %%E). \"\n"
    "            \"CFDWARP wants absolute pressure; check the FLUENT operating pressure\",\n"
    "            r[WARP_IP], (double)xc[1], r[WARP_IP] - warp_Pop, warp_Pop);\n"
    "    return 0;\n"
    "#endif\n"
    "  }\n"
    "#endif\n"
    "#if WARP_NW > 0\n"
    "  wsum = 0.0;\n"
    "  for (k = 0; k < WARP_NW; k++) wsum = wsum + r[WARP_IW + k];\n"
    "  if (wsum < WARP_WSUMMIN) {\n"
    "#if WARP_WCLAMP && defined(WARP_WFALLBACK)\n"
    "    {\n"
    "      static const double warp_wfb[] = { WARP_WFALLBACK };\n"
    "      if ((int)(sizeof(warp_wfb)/sizeof(warp_wfb[0])) != (int)WARP_NW) {\n"
    "        sprintf(warp_why,\n"
    "                \"WARP_WFALLBACK lists %%d mass fractions but vars_fluid has %%d\",\n"
    "                (int)(sizeof(warp_wfb)/sizeof(warp_wfb[0])), (int)WARP_NW);\n"
    "        return 0;\n"
    "      }\n"
    "      for (k = 0; k < WARP_NW; k++) r[WARP_IW + k] = warp_wfb[k];\n"
    "      wsum = 0.0;\n"
    "      for (k = 0; k < WARP_NW; k++) wsum = wsum + r[WARP_IW + k];\n"
    "      if (wsum < WARP_WSUMMIN) {\n"
    "        sprintf(warp_why, \"WARP_WFALLBACK itself sums to %%E\", wsum);\n"
    "        return 0;\n"
    "      }\n"
    "      warp_note_clamp(\"species sum\", 0.0, wsum, (double)xc[1]);\n"
    "    }\n"
    "#else\n"
    "    sprintf(warp_why,\n"
    "            \"the %%d species mass fractions sum to %%E at y=%%E. Either this \"\n"
    "            \"FLUENT case carries no species, or the F_YI indices in \"\n"
    "            \"warp_fill_record() do not name real FLUENT species\",\n"
    "            (int)WARP_NW, wsum, (double)xc[1]);\n"
    "    return 0;\n"
    "#endif\n"
    "  }\n"
    "  if (fabs(wsum - 1.0) > WARP_WSUMTOL) {\n"
    "    warp_nwarn = warp_nwarn + 1;\n"
    "    if (warp_nwarn == 1)\n"
    "      WARP_PRINT(\"\\nwarp: WARNING species fractions sum to %%E, not 1, at y=%%E;\"\n"
    "                 \" renormalising. FLUENT has fewer species than CFDWARP, or an\"\n"
    "                 \" F_YI index is wrong.\\n\", wsum, (double)xc[1]);\n"
    "  }\n"
    "  for (k = 0; k < WARP_NW; k++) r[WARP_IW + k] = r[WARP_IW + k]/wsum;\n"
    "#else\n"
    "  wsum = 1.0;   /* no species in vars_fluid; silences an unused warning */\n"
    "  if (wsum < 0.0) return 0;\n"
    "#endif\n"
    "  return 1;\n");
  fprintf(udf, "}\n\n");

  fprintf(udf,
    "#if WARP_COUPLED\n"
    "/* Does a file exist?  fopen rather than stat or access, because those two\n"
    "   are not reliably available to the interpreter.  This one and warp_stamp()\n"
    "   below are used by the rendezvous only, so they are compiled with it.  They\n"
    "   are plain static functions rather than DEFINE_ hooks, which is why an #if\n"
    "   around them is safe -- see the note further down. */\n"
    "static int warp_exists(const char *path)\n"
    "{\n"
    "  FILE *fp;\n"
    "  fp = fopen(path, \"rb\");\n"
    "  if (fp == NULL) return 0;\n"
    "  fclose(fp);\n"
    "  return 1;\n"
    "}\n"
    "\n"
    "/* Publish a zero-byte stamp.  With WARP_USE_RENAME 1 it goes through a\n"
    "   temporary name; with 0 it is created directly, which is equally safe for a\n"
    "   file that is empty by definition -- it is whole the instant it exists. */\n"
    "static void warp_stamp(const char *path)\n"
    "{\n"
    "  FILE *fp;\n"
    "  char tmpname[512];\n"
    "#if WARP_USE_RENAME\n"
    "  sprintf(tmpname, \"%%s.stamptmp\", path);\n"
    "#else\n"
    "  strcpy(tmpname, path);\n"
    "#endif\n"
    "  fp = fopen(tmpname, \"wb\");\n"
    "  if (!fp) { WARP_PRINT(\"\\nwarp: cannot open %%s\\n\", tmpname); return; }\n"
    "  fclose(fp);\n"
    "#if WARP_USE_RENAME\n"
    "  if (rename(tmpname, path) != 0)\n"
    "    WARP_PRINT(\"\\nwarp: could not rename %%s to %%s\\n\", tmpname, path);\n"
    "#endif\n"
    "}\n"
    "#endif   /* WARP_COUPLED */\n"
    "\n"
    "static void warp_write_file(const char *finalname)\n"
    "{\n"
    "  FILE *fp;\n"
    "  char tmpname[512];\n"
    "  int k;\n"
    "#if WARP_USE_RENAME\n"
    "  sprintf(tmpname, \"%%s.tmp\", finalname);\n"
    "#else\n"
    "  strcpy(tmpname, finalname);\n"
    "#endif\n"
    "  fp = fopen(tmpname, \"wb\");\n"
    "  if (!fp) { WARP_PRINT(\"\\nwarp: cannot open %%s\\n\", tmpname); return; }\n"
    "  /* iter, CFL, effiter, windowis and windowie are the CFDWARP run's own state.\n"
    "     They are written here only so the header is well formed; CFDWARP must not\n"
    "     adopt them on a zone read. */\n"
    "  fprintf(fp, \"WARPINTFORMAT001 numnodes=%%d nf=%%d nd=%%d ns=%%d\"\n"
    "              \" windowis=%%d windowie=%%d iter=%%d effiter_U=%%E effiter_R=%%E\"\n"
    "              \" CFL=%%E time=%%E dt=%%E vars_fluid=\\\"%s\\\" vars_emfield=\\\"NONE\\\"\\n\",\n"
    "          warp_nrec, %ld, WARP_ND, WARP_NS, %ld, %ld, %ld, %.8E, %.8E, %.8E,\n"
    "          (double)CURRENT_TIME, (double)CURRENT_TIMESTEP);\n"
    "  for (k = 0; k < warp_nrec; k++)\n"
    "    fwrite(&warp_buf[k*WARP_RECLEN], sizeof(double), WARP_RECLEN, fp);\n"
    "  fclose(fp);\n"
    "#if WARP_USE_RENAME\n"
    "  /* rename() is atomic on the same filesystem, so CFDWARP never sees a file\n"
    "     that is half written -- a plain fopen(finalname) would race. */\n"
    "  if (rename(tmpname, finalname) != 0) {\n"
    "    WARP_PRINT(\"\\nwarp: could not rename %%s to %%s\\n\", tmpname, finalname);\n"
    "    return;\n"
    "  }\n"
    "#endif\n"
    "  WARP_PRINT(\"\\nwarp: wrote %%d faces to %%s\\n\", warp_nrec, finalname);\n"
    "  /* Once per FLUENT session, name the UDF that is running. A regenerated\n"
    "     file that FLUENT never compiled is otherwise indistinguishable here. */\n"
    "  if (!warp_stamped) {\n"
    "    warp_stamped = 1;\n"
    "    WARP_PRINT(\"warp: this UDF was generated %%s in %%s\\n\",\n"
    "               WARP_F2W_GEN_TIME, WARP_F2W_GEN_DIR);\n"
    "    WARP_PRINT(\"warp:   if that is not the folder FLUENT is running in, FLUENT compiled\\n\");\n"
    "    WARP_PRINT(\"warp:   an old UDF_fluent2warp.c, not the one you just generated.\\n\");\n"
    "  }\n"
    "}\n"
    "\n", vars_fluid_str, nf, windowis, windowie, iter, effiter_U, effiter_R, CFL);

  fprintf(udf,
    "/* The decision, kept in one place because the serial and the parallel sweep\n"
    "   both reach it.  A sweep is publishable only if it refused nothing and\n"
    "   collected something: an empty file is a valid CFDWARP interpolation file\n"
    "   with numnodes=0, and the reader accepts it and changes not one node, which\n"
    "   is the same silent nothing a stock unpatched binary produces. */\n"
    "static int warp_publish(const char *finalname)\n"
    "{\n"
    "  if (warp_nbad > 0) {\n"
    "    WARP_PRINT(\"\\nwarp: REFUSING to send exchange data: %%d of %%d faces are \"\n"
    "               \"unusable.\\n      first one: %%s\\n\",\n"
    "               warp_nbad, warp_nbad + warp_nrec, warp_why);\n"
    "#if WARP_IP >= 0\n"
    "    /* The two causes of a low interface pressure do not look alike, and this\n"
    "       one line separates them without leaving FLUENT.  Exactly 0.0 on some or\n"
    "       all faces is a boundary condition handing over its default, not a flow\n"
    "       field that came out low: the commonest is a pressure-inlet whose\n"
    "       Supersonic/Initial Gauge Pressure was left at 0 while the face went\n"
    "       supersonic, in which case that field IS the static pressure FLUENT\n"
    "       applies and F_P returns it. */\n"
    "    WARP_PRINT(\"      lowest absolute P on the zone: %%E Pa; %%d of %%d faces \"\n"
    "               \"are exactly 0.\\n\",\n"
    "               warp_pmin_seen, warp_nzerop, warp_nbad + warp_nrec);\n"
    "    if (warp_nzerop > 0)\n"
    "      WARP_PRINT(\"      an exact 0 is a boundary condition default, not a \"\n"
    "                 \"computed pressure. Check\\n      Supersonic/Initial Gauge \"\n"
    "                 \"Pressure on the interface zone -- for a supersonic\\n\"\n"
    "                 \"      face that field is the static pressure, and its \"\n"
    "                 \"default is 0.\\n\");\n"
    "#endif\n"
    "    return 0;\n"
    "  }\n"
    "  if (warp_nrec == 0) {\n"
    "    if (warp_why[0] == '\\0')\n"
    "      sprintf(warp_why, \"the interface zone contributed no faces at all\");\n"
    "    WARP_PRINT(\"\\nwarp: REFUSING to send exchange data: %%s\\n\", warp_why);\n"
    "    return 0;\n"
    "  }\n"
    "  if (warp_nclamp > 0) {\n"
    "    /* Every sweep, not once. A run that clamped its way to a converged\n"
    "       interface has to leave that in the transcript on every exchange, or the\n"
    "       converged answer is indistinguishable from one nobody had to invent. */\n"
    "    WARP_PRINT(\"\\nwarp: WARNING %%d of %%d faces were altered before being sent\"\n"
    "               \" (WARP_CLAMP %%d, WARP_WCLAMP %%d).\\n      first one: %%s\\n\"\n"
    "               \"      CFDWARP receives the substituted value as a \"\n"
    "               \"measurement. This is not a fix.\\n\",\n"
    "               warp_nclamp, warp_nrec, (int)WARP_CLAMP, (int)WARP_WCLAMP,\n"
    "               warp_clamped);\n"
    "  }\n"
    "  warp_write_file(finalname);\n"
    "  return 1;\n"
    "}\n"
    "\n"
    "#if WARP_LOCALDX\n"
    "/* Give every face the support its own neighbourhood needs.\n"
    "\n"
    "   CFDWARP weights this file's points with w = max(1e-16, max(3e-4 - 1e-5 s,\n"
    "   1 - s)), s being the offset from the CFDWARP node measured in units of\n"
    "   [dx1 dx2], and it only looks at points closer than\n"
    "   sqrt(1.1 sum_dim (|dx1|+|dx2|)^2).  Both are set here, per face.\n"
    "\n"
    "   h is the distance to the SECOND nearest face centre, not the nearest: on a\n"
    "   stretched row of face centres that is the wider of the two neighbouring\n"
    "   gaps, and the wider gap is the one a CFDWARP node sitting in the middle of\n"
    "   it has to be reached across.  The nearest face centre to any point of the\n"
    "   interface lies within h/2, so every CFDWARP node gets s <= 0.5 whatever the\n"
    "   stretching, and the two faces bracketing it dominate: the next ones out sit\n"
    "   at s >= 1.5, worth 3e-4 apiece against 0.5.\n"
    "\n"
    "   O(nrec^2), a few thousand operations once per exchange, nothing against one\n"
    "   FLUENT time step, and it needs no sort -- which matters, because the records\n"
    "   arrive in face-loop order and, in parallel, partition by partition. */\n"
    "static void warp_set_local_support(void)\n"
    "{\n"
    "  int i, j, dim;\n"
    "  double d2, d2a, d2b, h, hmin, hmax;\n"
    "  double *ri;\n"
    "  double *rj;\n"
    "  if (warp_nrec < 3) return;   /* too few points to measure a spacing from */\n"
    "  hmin = 0.0; hmax = 0.0;\n"
    "  for (i = 0; i < warp_nrec; i++) {\n"
    "    ri = &warp_buf[i*WARP_RECLEN + WARP_NUMVARS];\n"
    "    d2a = WARP_BIG; d2b = WARP_BIG;\n"
    "    for (j = 0; j < warp_nrec; j++) {\n"
    "      if (j == i) continue;\n"
    "      rj = &warp_buf[j*WARP_RECLEN + WARP_NUMVARS];\n"
    "      d2 = 0.0;\n"
    "      for (dim = 0; dim < WARP_ND; dim++) d2 += (ri[dim]-rj[dim])*(ri[dim]-rj[dim]);\n"
    "      if (d2 < d2a) { d2b = d2a; d2a = d2; }\n"
    "      else if (d2 < d2b) d2b = d2;\n"
    "    }\n"
    "    if (!(d2b > 0.0) || d2b >= WARP_BIG) continue;  /* keep the constant pair */\n"
    "    h = sqrt(d2b);\n"
    "    for (dim = 0; dim < WARP_ND; dim++) {\n"
    "      warp_buf[i*WARP_RECLEN + WARP_NUMVARS + 1*WARP_ND + dim] = warp_dx1u[dim]*h;\n"
    "      warp_buf[i*WARP_RECLEN + WARP_NUMVARS + 2*WARP_ND + dim] = warp_dx2u[dim]*h;\n"
    "#if WARP_ND > 2\n"
    "      warp_buf[i*WARP_RECLEN + WARP_NUMVARS + 3*WARP_ND + dim] = warp_dx3u[dim]*h;\n"
    "#endif\n"
    "    }\n"
    "    if (hmax == 0.0) { hmin = h; hmax = h; }\n"
    "    if (h < hmin) hmin = h;\n"
    "    if (h > hmax) hmax = h;\n"
    "  }\n"
    "  if (!warp_dx_reported && hmax > 0.0) {\n"
    "    warp_dx_reported = 1;\n"
    "    WARP_PRINT(\"warp: per-face support, local spacing %%.4E to %%.4E m over %%d\"\n"
    "               \" faces; the single WARP_DX pair was %%.4E m everywhere.\\n\",\n"
    "               hmin, hmax, warp_nrec,\n"
    "               WARP_DX2MAG);\n"
    "  }\n"
    "}\n"
    "#endif\n"
    "\n"
    "/* Sweep the interface zone into warp_buf and publish it under finalname.\n"
    "   Returns 2 if THIS rank wrote the file, 1 if the call succeeded but this rank\n"
    "   had nothing to publish, 0 if the sweep produced nothing CFDWARP could use --\n"
    "   in which case warp_why says what was wrong and NOTHING is published, because\n"
    "   a bad payload costs far more to diagnose from the CFDWARP end than a missing\n"
    "   one does.\n"
    "\n"
    "   The 2-versus-1 distinction is not cosmetic: it is what keeps the stamp behind\n"
    "   the payload.  In a parallel FLUENT -- which is any -t, including -t1 -- the\n"
    "   HOST also runs DEFINE_EXECUTE_AT_END, its half of this function is compiled\n"
    "   out by the #if RP_NODE below, so it writes nothing and falls through to the\n"
    "   final return.  When that return was a plain 1 the host went straight on to\n"
    "   warp_stamp() and could publish the .ready stamp before node-0 had written the\n"
    "   payload; CFDWARP then failed the exchange with \"could not be staged as\n"
    "   f2w.int: No such file or directory\".  Only the rank that returns 2 stamps. */\n"
    "static int warp_collect_and_write(const char *finalname)\n"
    "{\n"
    "  Domain *d;\n"
    "  Thread *t;\n"
    "  face_t f;\n"
    "\n"
    "  warp_Pop = RP_Get_Real(\"operating-pressure\");\n"
    "  d = Get_Domain(1);\n"
    "  t = Lookup_Thread(d, WARP_ZONE_ID);\n"
    "  warp_nrec  = 0;\n"
    "  warp_nbad  = 0;\n"
    "  warp_nwarn = 0;\n"
    "  warp_nclamp = 0;\n"
    "  warp_nzerop = 0;\n"
    "  warp_pmin_seen = 0.0;\n"
    "  strcpy(warp_why, \"\");\n"
    "  strcpy(warp_clamped, \"\");\n"
    "\n"
    "#if !PARALLEL\n"
    "  /* serial: host and node are the same process, no message passing at all.\n"
    "     Note that a UDF whose face loop sits inside #if RP_NODE collects nothing\n"
    "     in serial and writes an empty file. */\n"
    "  if (NNULLP(t)) {\n"
    "    begin_f_loop(f, t) {\n"
    "      if (warp_nrec < WARP_MAXFACES) {\n"
    "        if (warp_fill_record(&warp_buf[warp_nrec*WARP_RECLEN], f, t)) warp_nrec++;\n"
    "        else warp_nbad++;\n"
    "      }\n"
    "    } end_f_loop(f, t)\n"
    "  } else {\n"
    "    sprintf(warp_why,\n"
    "            \"Lookup_Thread() found no zone with ID %%d. Set WARP_ZONE_ID to the \"\n"
    "            \"FLUENT zone ID of the interface\", (int)WARP_ZONE_ID);\n"
    "  }\n"
    "#if WARP_LOCALDX\n"
    "  warp_set_local_support();\n"
    "#endif\n"
    "  return warp_publish(finalname) ? 2 : 0;   /* 2 = this rank published */\n"
    "#else\n"
    "  /* parallel: every node fills its own part, node-0 collects it (only node-0\n"
    "     may talk to the host), and node-0 writes the file. Writing from the node\n"
    "     rather than the host keeps the data off the host entirely. */\n"
    "#if RP_NODE\n"
    "  {\n"
    "    int i, n[4];\n"
    "    if (NNULLP(t)) {\n"
    "      begin_f_loop(f, t) {\n"
    "        if (PRINCIPAL_FACE_P(f, t) && warp_nrec < WARP_MAXFACES) {\n"
    "          if (warp_fill_record(&warp_buf[warp_nrec*WARP_RECLEN], f, t)) warp_nrec++;\n"
    "          else warp_nbad++;\n"
    "        }\n"
    "      } end_f_loop(f, t)\n"
    "    }\n"
    "    if (I_AM_NODE_ZERO_P) {\n"
    "      compute_node_loop_not_zero(i) {\n"
    "        /* good records, refused ones, floored ones, exact zeros -- in that\n"
    "           order.  A node that refused every face must not look like a node\n"
    "           that simply holds no part of the interface, and a node that floored\n"
    "           every face must not look like one that had nothing to floor.\n"
    "           warp_pmin_seen is NOT reduced: it stays node-0's own minimum, so in\n"
    "           a parallel run read the count of exact zeros, which is reduced, and\n"
    "           treat the printed minimum as a lower bound from one partition. */\n"
    "        PRF_CRECV_INT(i, n, 4, i);\n"
    "        warp_nbad += n[1];\n"
    "        warp_nclamp += n[2];\n"
    "        warp_nzerop += n[3];\n"
    "        if (n[0] > 0 && warp_nrec + n[0] <= WARP_MAXFACES) {\n"
    "          PRF_CRECV_REAL(i, &warp_buf[warp_nrec*WARP_RECLEN], n[0]*WARP_RECLEN, i);\n"
    "          warp_nrec += n[0];\n"
    "        }\n"
    "      }\n"
    "#if WARP_LOCALDX\n"
    "      warp_set_local_support();\n"
    "#endif\n"
    "      return warp_publish(finalname) ? 2 : 0;   /* 2 = this rank published */\n"
    "    } else {\n"
    "      n[0] = warp_nrec;\n"
    "      n[1] = warp_nbad;\n"
    "      n[2] = warp_nclamp;\n"
    "      n[3] = warp_nzerop;\n"
    "      PRF_CSEND_INT(node_zero, n, 4, myid);\n"
    "      if (n[0] > 0) PRF_CSEND_REAL(node_zero, warp_buf, n[0]*WARP_RECLEN, myid);\n"
    "    }\n"
    "  }\n"
    "#endif\n"
    "  return 1;   /* the host and every node but node-0: nothing published, and\n"
    "                 deliberately NOT 2, so that warp_exchange does not stamp */\n"
    "#endif\n"
    "}\n"
    "\n");

  fprintf(udf,
    "/* ---------------- rendezvous with the CFDWARP run ------------------------\n"
    "   Every message is a NEW file carrying the exchange number, published by\n"
    "   rename() once it is whole, and nothing is ever deleted.  A single flag\n"
    "   that one side creates and the other deletes cannot tell \"the peer has not\n"
    "   answered yet\" apart from \"the peer answered and I consumed it\", so a flag\n"
    "   left over from a previous run makes one side read stale data and a flag\n"
    "   deleted twice deadlocks both.\n"
    "\n"
    "       xchg/w2f.<n>.int    payload, CFDWARP -> FLUENT\n"
    "       xchg/w2f.<n>.ready  stamp, appears only once the payload is whole AND\n"
    "                           warp2fluent has turned it into the profile file\n"
    "       xchg/w2f.<n>.installed  WARP_TRANSIENT only: written by UDF.c's\n"
    "                           warp_wait_and_load once it has actually loaded\n"
    "                           that profile. .ready says CFDWARP sent it;\n"
    "                           .installed says FLUENT is using it, and only the\n"
    "                           second one licenses an answer\n"
    "       xchg/f2w.<n>.int    payload, FLUENT -> CFDWARP, written here\n"
    "       xchg/f2w.<n>.ready  stamp, written here\n"
    "       xchg/STOP           either side may create it; both then leave\n"
    "       xchg/ERROR          this side found data it will not send; the text in\n"
    "                           it says what, and exchange.sh fails on it at once\n"
    "                           instead of waiting out COUPLE_TIMEOUT\n"
    "\n"
    "   Nothing here refers to a symbol of UDF.c.  Reloading the profiles CFDWARP\n"
    "   sends is UDF.c's own warp_reload_at_end, watching the same stamps, because\n"
    "   an interpreted UDF cannot resolve another interpreted file's symbols.  Both\n"
    "   Execute-At-End functions have to be selected under Function Hooks.\n"
    "\n"
    "   This function does NOT block.  It looks once per iteration and returns, so\n"
    "   FLUENT keeps iterating on the profile it already has while CFDWARP works.\n"
    "   A blocking wait inside DEFINE_EXECUTE_AT_END would freeze the FLUENT GUI,\n"
    "   burn a core, and needs a portable sleep the interpreter does not have.\n"
    "   The extra iterations are not wasted: they are ordinary FLUENT iterations\n"
    "   on the last profile CFDWARP sent.\n"
    "\n"
    "   THAT REASONING IS A STEADY-RUN REASONING and it does not carry over.  In a\n"
    "   transient run the extra passes are not iterations on the last profile,\n"
    "   they are TIME STEPS taken on it, and the physical time FLUENT has reached\n"
    "   is then no longer the time CFDWARP thinks it is at.  So WARP_TRANSIENT 1\n"
    "   moves the waiting into warp_wait_and_load in UDF.c, an Adjust hook, which\n"
    "   blocks at the START of the time step, where the profile can still change\n"
    "   the step about to be solved.  This function stays non-blocking either way;\n"
    "   in the transient case there is simply nothing left for it to wait for.\n"
    "\n"
    "   BOTH Execute-At-End functions below are compiled whatever WARP_COUPLED is\n"
    "   set to, and each one tests it itself.  They have to be: Build generates\n"
    "   libudf/src/udf_names.c by scanning this file, as text and with no\n"
    "   preprocessing, for lines whose first word begins with a DEFINE_ prefix,\n"
    "   and it declares every name it finds and puts it in udf_data[].  One that\n"
    "   sits in the dead half of an #if is then declared and never defined; the\n"
    "   library still links, a shared object being allowed to carry undefined\n"
    "   symbols, and Load then fails with\n"
    "       Error: undefined symbol: warp_exchange\n"
    "       Error: undefined symbol: write_cfdwarp_interpolation_file\n"
    "   leaving nothing in the library available.  Select the one that matches\n"
    "   WARP_COUPLED under Function Hooks; the other one says so and returns.\n"
    "   ---------------------------------------------------------------------- */\n"
    "\n"
    "#if WARP_COUPLED\n"
    "/* Write the refusal where the CFDWARP side can see it, then stop taking part.\n"
    "   Leaving no answer at all would work too -- exchange.sh times out -- but that\n"
    "   costs COUPLE_TIMEOUT seconds and says nothing about why.  A plain static\n"
    "   function, not a DEFINE_, so an #if around it is harmless. */\n"
    "static void warp_error(const char *why)\n"
    "{\n"
    "  FILE *fp;\n"
    "  char path[512];\n"
    "  sprintf(path, \"%%s/ERROR\", WARP_XCHG);\n"
    "  fp = fopen(path, \"w\");\n"
    "  if (fp) { fprintf(fp, \"FLUENT side refused to send: %%s\\n\", why); fclose(fp); }\n"
    "  WARP_PRINT(\"warp: wrote %%s; the CFDWARP side will stop on it.\\n\", path);\n"
    "}\n"
    "\n"
    "static int warp_n     = 0;   /* exchange number this side is working on */\n"
    "static int warp_have  = 0;   /* CFDWARP's message n is in and installed */\n"
    "#if !WARP_TRANSIENT\n"
    "static int warp_iter0 = 0;   /* N_ITER when it was installed            */\n"
    "#endif\n"
    "static int warp_done  = 0;   /* STOP seen, this side is finished        */\n"
    "static int warp_resumed = 0; /* warp_n recovered from the directory yet? */\n"
    "\n"
    "/* warp_n is a static of this library, and a library static goes back to 0\n"
    "   every time FLUENT reloads libudf -- which re-reading the case does.  The\n"
    "   CFDWARP side does not have that problem: it keeps its counter in\n"
    "   xchg/seq, on disk, because a SOAP control file cannot hold one.  Left\n"
    "   alone, the two counters therefore drift apart the first time the case is\n"
    "   re-read, and since the rendezvous never deletes a message this side then\n"
    "   answers requests it has already answered while CFDWARP consumes replies\n"
    "   it never asked for.  Recover the number from the directory instead, the\n"
    "   same way the other side recovers it from xchg/seq: the highest answered\n"
    "   exchange plus one.  Only fopen is used, so this stays inside what the\n"
    "   interpreted UDF can do -- no opendir, no scandir. */\n"
    "static void warp_resume(void)\n"
    "{\n"
    "  char path[512];\n"
    "  int k;\n"
    "  warp_resumed = 1;\n"
    "  for (k = 0; k < 1000000; k++) {\n"
    "    sprintf(path, \"%%s/f2w.%%d.ready\", WARP_XCHG, k);\n"
    "    if (!warp_exists(path)) break;\n"
    "    warp_n = k + 1;\n"
    "  }\n"
    "  if (warp_n > 0) {\n"
    "    WARP_PRINT(\"\\nwarp: %%s already holds %%d answered exchanges; resuming at %%d.\\n\",\n"
    "               WARP_XCHG, warp_n, warp_n);\n"
    "    WARP_PRINT(\"warp:   this counter is a libudf static and restarts at 0 whenever the\\n\");\n"
    "    WARP_PRINT(\"warp:   case is re-read, so it is taken from the directory, not from 0.\\n\");\n"
    "    WARP_PRINT(\"warp:   If you meant to start a new coupled run, stop FLUENT, empty %%s\\n\", WARP_XCHG);\n"
    "    WARP_PRINT(\"warp:   and restart both codes -- do not delete %%s/seq on its own.\\n\", WARP_XCHG);\n"
    "  }\n"
    "}\n"
    "#endif\n"
    "\n"
    "DEFINE_EXECUTE_AT_END(warp_exchange)\n"
    "{\n"
    "#if !WARP_COUPLED\n"
    "  /* Compiled but inactive, so that this name resolves at Load whatever\n"
    "     WARP_COUPLED is.  Says so once: a hook selected under Function Hooks\n"
    "     that returns in silence is how a coupling looks alive and is not. */\n"
    "  static int warp_said = 0;\n"
    "  if (!warp_said) {\n"
    "    warp_said = 1;\n"
    "    WARP_PRINT(\"warp: warp_exchange is inactive, WARP_COUPLED is 0.\\n\");\n"
    "    WARP_PRINT(\"      Select write_cfdwarp_interpolation_file instead, or set\\n\");\n"
    "    WARP_PRINT(\"      WARP_COUPLED to 1 and rebuild for a coupled run.\\n\");\n"
    "  }\n"
    "#else\n"
    "  char path[512];\n"
    "  int warp_rc;\n"
    "\n"
    "  if (warp_done) return;\n"
    "  if (!warp_resumed) warp_resume();\n"
    "\n"
    "  sprintf(path, \"%%s/STOP\", WARP_XCHG);\n"
    "  if (warp_exists(path)) {\n"
    "    warp_done = 1;\n"
    "    WARP_PRINT(\"\\nwarp: STOP seen after %%d exchanges, CFDWARP has finished.\\n\", warp_n);\n"
    "    WARP_PRINT(\"warp: this zone keeps the last profile CFDWARP sent; stop FLUENT from\\n\");\n"
    "    WARP_PRINT(\"      the journal or the GUI once its own residuals are down.\\n\");\n"
    "    return;\n"
    "  }\n"
    "\n"
    "#if WARP_TRANSIENT\n"
    "  /* ---- transient ---------------------------------------------------------\n"
    "     One exchange per WARP_INTERVAL time steps, I for short.  The profile for\n"
    "     this window was installed at the START of time step m with (m-1)%%I == 0,\n"
    "     by warp_wait_and_load in UDF.c, which blocked until it was there.  So the\n"
    "     answer is due here at the END of the last step of the window, m%%I == 0,\n"
    "     and there is nothing to wait for: the state being sent is FLUENT's state\n"
    "     at physical time CURRENT_TIME, computed from CFDWARP's own message.\n"
    "     Both halves derive the window from N_TIME and WARP_INTERVAL alone, so\n"
    "     they share no symbol and cannot disagree about whose turn it is. */\n"
    "  if ((int)N_TIME < 1) return;              /* the call before the march */\n"
    "  if ((int)N_TIME %% WARP_INTERVAL != 0) return;\n"
    "\n"
    "  sprintf(path, \"%%s/w2f.%%d.installed\", WARP_XCHG, warp_n);\n"
    "  if (!warp_exists(path)) {\n"
    "    /* The profile for this window was never INSTALLED -- which is not the\n"
    "       same as CFDWARP not having sent it.  warp_wait_and_load in UDF.c\n"
    "       creates this file after it has loaded the profile, so its absence at\n"
    "       the end of the window means that hook is not selected under Function\n"
    "       Hooks -> Adjust, and the window was solved on whatever FLUENT was\n"
    "       already holding.  Answering from it would hand CFDWARP a state that\n"
    "       is not a response to its own message, which is invisible from the\n"
    "       CFDWARP end and produces a plausible converged-looking wrong answer.\n"
    "       So: say so, once, and let CFDWARP time out instead. */\n"
    "    if (!warp_have) {\n"
    "      warp_have = 1;   /* reused here only to say this once */\n"
    "      WARP_PRINT(\"\\nwarp: NOT answering exchange %%d at time step %%d: %%s is not\"\n"
    "                 \" there.\\n\", warp_n, (int)N_TIME, path);\n"
    "      WARP_PRINT(\"      This is the transient build (WARP_TRANSIENT 1). That file is\\n\");\n"
    "      WARP_PRINT(\"      written by warp_wait_and_load in UDF.c once it has installed the\\n\");\n"
    "      WARP_PRINT(\"      profile, so it is missing because that hook is not selected under\\n\");\n"
    "      WARP_PRINT(\"      Function Hooks -> Adjust. Select it and restart the run.\\n\");\n"
    "    }\n"
    "    return;\n"
    "  }\n"
    "  warp_have = 0;\n"
    "\n"
    "#else\n"
    "  if (!warp_have) {\n"
    "    /* Waiting for CFDWARP.  Its exchange.sh has already run warp2fluent on\n"
    "       w2f.<n>.int, so the interpolated profile file that the DEFINE_PROFILEs\n"
    "       in UDF.c read is on disk by the time this stamp appears; UDF.c's own\n"
    "       warp_reload_at_end picks it up on this same iteration. */\n"
    "    sprintf(path, \"%%s/w2f.%%d.ready\", WARP_XCHG, warp_n);\n"
    "    if (!warp_exists(path)) return;\n"
    "    warp_have  = 1;\n"
    "    warp_iter0 = (int)N_ITER;\n"
    "    WARP_PRINT(\"\\nwarp: exchange %%d in at iteration %%d\\n\", warp_n, (int)N_ITER);\n"
    "    return;\n"
    "  }\n"
    "\n"
    "  /* Have it.  Give FLUENT WARP_INTERVAL iterations to respond to the new\n"
    "     inflow profile, then publish payload and stamp, in that order. */\n"
    "  if ((int)N_ITER - warp_iter0 < WARP_INTERVAL) return;\n"
    "#endif\n"
    "\n"
    "  sprintf(path, \"%%s/f2w.%%d.int\", WARP_XCHG, warp_n);\n"
    "  warp_rc = warp_collect_and_write(path);\n"
    "  if (!warp_rc) {\n"
    "    /* Deliberately publish no stamp: CFDWARP must not receive this. */\n"
    "    warp_error(warp_why);\n"
    "    warp_done = 1;\n"
    "    return;\n"
    "  }\n"
    "  sprintf(path, \"%%s/f2w.%%d.ready\", WARP_XCHG, warp_n);\n"
    "  /* Only the rank that actually wrote the payload may say it is there.  The\n"
    "     host writes nothing -- its half of the sweep is inside #if RP_NODE -- and\n"
    "     stamping from it raced node-0 and killed coupled runs mid-exchange. */\n"
    "  if (warp_rc == 2) warp_stamp(path);\n"
    "#if WARP_TRANSIENT\n"
    "  WARP_PRINT(\"warp: exchange %%d answered at the end of time step %%d, t=%%.6E s,\"\n"
    "             \" dt=%%.6E s\\n\", warp_n, (int)N_TIME, (double)CURRENT_TIME,\n"
    "             (double)CURRENT_TIMESTEP);\n"
    "#else\n"
    "  WARP_PRINT(\"warp: exchange %%d answered at iteration %%d\\n\", warp_n, (int)N_ITER);\n"
    "#endif\n"
    "  warp_have = 0;\n"
    "  warp_n++;\n"
    "#endif\n"
    "}\n"
    "\n"
    "/* The one-shot writer: write WARP_FILENAME every WARP_INTERVAL iterations --\n"
    "   time steps, with WARP_TRANSIENT 1 -- and never wait for anybody.  Active\n"
    "   only when WARP_COUPLED is 0, compiled either way for the reason given\n"
    "   above.  It overwrites WARP_FILENAME each time, so what is left on disk at\n"
    "   the end of a transient run is the LAST time level written, not a history;\n"
    "   if you want the history, WARP_COUPLED 0 is the wrong mode -- run the\n"
    "   rendezvous, where every time level is kept as xchg/f2w.<n>.int. */\n"
    "DEFINE_EXECUTE_AT_END(write_cfdwarp_interpolation_file)\n"
    "{\n"
    "#if WARP_COUPLED\n"
    "  static int warp_said1 = 0;\n"
    "  if (!warp_said1) {\n"
    "    warp_said1 = 1;\n"
    "    WARP_PRINT(\"warp: write_cfdwarp_interpolation_file is inactive, WARP_COUPLED is 1.\\n\");\n"
    "    WARP_PRINT(\"      Select warp_exchange instead; it writes the same file, in turn\\n\");\n"
    "    WARP_PRINT(\"      with CFDWARP, as %%s/f2w.<n>.int.\\n\", WARP_XCHG);\n"
    "  }\n"
    "#else\n"
    "#if WARP_TRANSIENT\n"
    "  if (N_TIME == 0 || N_TIME %% WARP_INTERVAL != 0) return;\n"
    "#else\n"
    "  if (N_ITER == 0 || N_ITER %% WARP_INTERVAL != 0) return;\n"
    "#endif\n"
    "  (void)warp_collect_and_write(WARP_FILENAME);\n"
    "#endif\n"
    "}\n");
  fclose(udf);
}


int main(int argc, char **argv) {
  FILE *datafile, *meshdata, *fluent_output;
  bool VALIDOPTIONS = TRUE, IS_3D = FALSE, ISREAD_EMFIELD = TRUE;
  char *options, *input, *mesh, *output;
  char **initvar_fluid_names_file = NULL;
  char **initvar_fluid_names_orig = NULL;
  long numvars_fluid_orig;
  double xstation,dx1mean[3],dx2mean[3],dx3mean[3];
  char **initvar_emfield_names_file = NULL;
  int RET;
  options = NULL;
  bool FORMAT001;
  long cnt,cnterror,dim,l_file,l_mesh;
  long numflux_read,numspec_read,numdim_read,numnodes,numfaces,windowis,windowie,iter,numvars_fluid_file;
  long numflux_emfield_read,numdim_emfield_read,numnodes_emfield,numvars_emfield_file;
  double effiter_U,effiter_R,effiter_U_emfield,effiter_R_emfield;
  double tmp_dt,tmp_time,CFL,Lc,dxwarp,dxfluent;
  double *radiusmax2_file,*radiusmax2_emfield_file,thisweight,*weight,thisweight_emfield,*weight_emfield;
  char data_format_str[100], initvar_fluid_str_file[500], initvar_emfield_str_file[500];
  char **allvarnames, **allcnames;
  long numvars_all,numoutside,l_nearest;
  double dist2,dist2min,dist2worst;
  /* -clamp / -pfloor / -tfloor. pfloor_cl and tfloor_cl stay 0 when the flag was
     not given, which is the signal to fall back to the minimum of the CFDWARP
     interpolation file computed after the interpolation loop below. */
  int CLAMP=0;
  double pfloor_cl=0.0, tfloor_cl=0.0, pfloor_use, tfloor_use, pmean_file;
  long iPfile,iTfile;
  numfaces=0;
  /* numnodes_emfield is only set by the EMField pass of the interpolation file;
     without this it stays uninitialized when the CFDWARP build has no EMField
     module and is then used as a malloc size below */
  numnodes_emfield=0;
  numvars_emfield_file=0;

  input =(char *)malloc(400*sizeof(char));
  mesh =(char *)malloc(400*sizeof(char));
  output =(char *)malloc(400*sizeof(char));

  /* -exchange is read before the required-flag check because it supplies
     defaults for all three of them: in a control-file call the payload, the
     interface mesh and the profile file always have the same names, and a
     Cycle() line reading system("./warp2fluent -exchange -tent") is one fewer
     place for a path to be wrong. Any of the three may still be given. */
  for (cnt=1; cnt<argc; cnt++) {
    if (strcmp(argv[cnt],"-exchange")==0) { EXCHANGE=1; argv[cnt]=(char *)""; }
  }
  NOUDF=EXCHANGE;
  for (cnt=1; cnt<argc; cnt++) {
    if (strcmp(argv[cnt],"-noudf")==0) { NOUDF=1; argv[cnt]=(char *)""; }
    if (strcmp(argv[cnt],"-udf")==0)   { NOUDF=0; argv[cnt]=(char *)""; }
    if (strcmp(argv[cnt],"-transient")==0) { TRANSIENT=1; argv[cnt]=(char *)""; }
    if (strcmp(argv[cnt],"-globaldx")==0) { GLOBALDX=1; argv[cnt]=(char *)""; }
  }
  /* A transient exchange is one per time step unless asked otherwise; a steady
     one is every 20 iterations, which is what it has always been. Reading
     -interval below overrides either. */
  if (TRANSIENT) XINTERVAL=1;
  for (cnt=1; cnt<argc-1; cnt++) {
    if (strcmp(argv[cnt],"-xchg")==0) {
      snprintf(XDIR,sizeof(XDIR),"%s",argv[cnt+1]);
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
    if (strcmp(argv[cnt],"-staged")==0) {
      snprintf(XSTAGED,sizeof(XSTAGED),"%s",argv[cnt+1]);
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
    if (strcmp(argv[cnt],"-timeout")==0) {
      XTIMEOUT=atof(argv[cnt+1]);
      if (XTIMEOUT<=0.0) {
        fprintf(stderr,"\n\n-timeout must be a positive number of seconds, got %s\n\n",argv[cnt+1]);
        exit (EXIT_FAILURE);
      }
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
    if (strcmp(argv[cnt],"-poll")==0) {
      XPOLL=atof(argv[cnt+1]);
      if (XPOLL<=0.0 || XPOLL>60.0) {
        fprintf(stderr,"\n\n-poll must be a positive number of seconds not larger than 60, got %s\n\n",argv[cnt+1]);
        exit (EXIT_FAILURE);
      }
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
    if (strcmp(argv[cnt],"-zone")==0) {
      ZONEID=atoi(argv[cnt+1]);
      ZONEID_GIVEN=1;
      if (ZONEID<=0) {
        fprintf(stderr,"\n\n-zone must be a positive FLUENT zone ID, got %s\n\n",argv[cnt+1]);
        exit (EXIT_FAILURE);
      }
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
    if (strcmp(argv[cnt],"-interval")==0) {
      XINTERVAL=atoi(argv[cnt+1]);
      XINTERVAL_GIVEN=1;
      if (XINTERVAL<=0) {
        fprintf(stderr,"\n\n-interval must be a positive number of FLUENT %s, got %s\n\n",
                TRANSIENT?"time steps":"iterations",argv[cnt+1]);
        exit (EXIT_FAILURE);
      }
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
    if (strcmp(argv[cnt],"-fluentspecies")==0) {
      if (snprintf(FLSPEC,sizeof(FLSPEC),"%s",argv[cnt+1])>=(int)sizeof(FLSPEC)) {
        fprintf(stderr,"\n\n-fluentspecies list is too long\n\n");
        exit (EXIT_FAILURE);
      }
      if (FLSPEC[0]=='\0') {
        fprintf(stderr,"\n\n-fluentspecies was given an empty list\n\n");
        exit (EXIT_FAILURE);
      }
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
  }
  if (NOUDF && (ZONEID_GIVEN || XINTERVAL_GIVEN || FLSPEC[0]!='\0'))
    fprintf(stdout,"NOTE: -zone/-interval/-fluentspecies only affect the generated UDF_fluent2warp.c,\n"
                   "      and no UDF is being written on this call (-noudf, implied by -exchange).\n");
  /* -transient is NOT in that note: on an -exchange call it still does something,
     namely the clock check in xchg_finish() below. */
  if (!EXCHANGE && (strcmp(XDIR,"xchg")!=0 || XTIMEOUT!=600.0 || XPOLL!=0.2 || strcmp(XSTAGED,"f2w.int")!=0))
    fprintf(stdout,"NOTE: -xchg/-staged/-timeout/-poll only mean anything with -exchange, which was\n"
                   "      not given. This call is an ordinary one-shot conversion.\n");

  if (process_flag_string(argc, argv, "-ii", &input)!=2) {
    if (EXCHANGE) strcpy(input,"w2f.int"); else VALIDOPTIONS=FALSE;
  }
  if (process_flag_string(argc, argv, "-m", &mesh)!=2) {
    if (EXCHANGE) strcpy(mesh,"fluent_bdry.mesh"); else VALIDOPTIONS=FALSE;
  }
  if (process_flag_string(argc, argv, "-o", &output)!=2) {
    if (EXCHANGE) strcpy(output,"warp_output.dat"); else VALIDOPTIONS=FALSE;
  }

  if (!VALIDOPTIONS) {
    fprintf (stderr, "\nPROBLEM WITH FLAGS\n\n\nRequired and Optional Flags:\n\n"
             "Flag      Arg                     Arg Type     Required?\n"
             "----------------------------------------------------------\n"
             "-ii       CFDWARP datafile        string       YES\n"
             "-m        FLUENT mesh             string       YES\n"
             "-o        Output file             string       YES\n"
             "-tent     (no argument)            flag         NO\n"
             "-clamp    (no argument)            flag         NO\n"
             "-pfloor   pressure floor, Pa       double       NO\n"
             "-tfloor   temperature floor, K     double       NO\n"
             "-exchange (no argument)            flag         NO\n"
             "-noudf    (no argument)            flag         NO\n"
             "-udf      (no argument)            flag         NO\n"
             "-xchg     exchange directory       string       NO\n"
             "-staged   file Cycle() reads       string       NO\n"
             "-timeout  wait for the peer, s     double       NO\n"
             "-poll     directory look, s        double       NO\n"
             "-zone     FLUENT interface zone ID  int          NO\n"
             "-interval FLUENT iters per exchange int          NO\n"
             "-fluentspecies  FLUENT species order, comma separated  NO\n"
             "-transient (no argument)            flag         NO\n"
             "----------------------------------------------------------\n"
             "-transient generates the pair of UDFs for a TRANSIENT FLUENT run. It is\n"
             "not cosmetic: FLUENT calls DEFINE_EXECUTE_AT_END once per iteration when\n"
             "the solver is steady and once per TIME STEP when it is transient, so\n"
             "without it -interval counts a unit that no longer exists, and the\n"
             "non-blocking look that is right for a steady run lets FLUENT march\n"
             "physical time on a profile CFDWARP has not sent yet. With it:\n"
             "  - WARP_INTERVAL counts TIME STEPS and defaults to 1, one exchange per\n"
             "    time step. Give -interval n for one every n.\n"
             "  - UDF.c gains warp_wait_and_load, a DEFINE_ADJUST that BLOCKS at the\n"
             "    start of the time step until CFDWARP's message for it is on disk,\n"
             "    and warp_reload_at_end goes inert. Select warp_wait_and_load under\n"
             "    Function Hooks -> Adjust; it is not optional, and UDF_fluent2warp.c\n"
             "    refuses to answer rather than send back a step solved without it.\n"
             "  - both files must be COMPILED, not interpreted: the wait uses\n"
             "    nanosleep().\n"
             "On an -exchange call -transient instead compares the physical time and\n"
             "dt in FLUENT's answer against CFDWARP's own, and says so every exchange.\n"
             "----------------------------------------------------------\n"
             "-exchange runs one full CFDWARP <-> FLUENT exchange instead of a one-shot\n"
             "conversion: it publishes the payload CFDWARP just wrote, converts it,\n"
             "stamps it, waits for FLUENT's answer and stages that answer under the\n"
             "name the control file reads. It replaces exchange.sh, so Cycle() reads\n"
             "  WriteInterpolationFile(\"w2f.int\", ie,js+1, ie,je-1);\n"
             "  system(\"./warp2fluent -exchange -tent\");\n"
             "  ReadInterpolationFileZone(\"f2w.int\", ie,js+1, ie,je-1);\n"
             "With -exchange, -ii, -m and -o default to w2f.int, fluent_bdry.mesh and\n"
             "warp_output.dat, and -xchg, -staged, -timeout and -poll default to xchg,\n"
             "f2w.int, 600 s and 0.2 s. Every path is relative to CFDWARP's working\n"
             "directory, so there is nothing to place anywhere. -exchange also implies\n"
             "-noudf, which suppresses the rewriting of UDF.c and UDF_fluent2warp.c: an\n"
             "exchange does not change the UDF FLUENT compiled, and rewriting it once\n"
             "per exchange only advances its generation stamp away from the one baked\n"
             "into libudf, which makes checkudf.sh report a mismatch that is not one.\n"
             "Give -udf to rewrite them anyway.\n"
             "----------------------------------------------------------\n"
             "-clamp generates UDF_fluent2warp.c with WARP_CLAMP 1: a face whose\n"
             "absolute P or T is below the validity limit is raised to the floor and\n"
             "SENT, rather than the exchange being refused. It does not repair the\n"
             "face, it only stops the complaint, so leave it off until you know why\n"
             "the value is low. Without -pfloor/-tfloor the floors default to the\n"
             "smallest P and T in the CFDWARP interpolation file, which keeps a\n"
             "clamped face in the range the rest of the interface is in.\n"
             "----------------------------------------------------------\n"
             "-zone, -interval and -fluentspecies fill in the three case-specific\n"
             "settings of UDF_fluent2warp.c, the only ones this tool cannot read out\n"
             "of the CFDWARP interpolation file. Give all three and the generated UDF\n"
             "needs no editing at all; leave one out and it is emitted as a default\n"
             "with a CHECK comment, as before. -fluentspecies takes the species in the\n"
             "order your FLUENT mixture material lists them, bulk species last, eg\n"
             "  ./warp2fluent -ii w2f.int -m fluent_bdry.mesh -o warp_output.dat \\\n"
             "                -zone 13 -interval 45 -fluentspecies n,o2,o,no,n2\n"
             "The names are matched against the CFDWARP species without regard to\n"
             "case, and a CFDWARP species the list does not name is a hard error, not\n"
             "a fallback: a wrong species mapping is otherwise silent all the way to\n"
             "mass fractions arriving on the wrong species.\n"
             "----------------------------------------------------------\n"
             "Eg: \n"
             "./warp2fluent -ii datafile -m fluentmesh.dat -o warp_output.dat\nWill read input CFDWARP interpolation datafile (2D or 3D) obtained with the -oi flag, the FLUENT mesh and assign the CFDWARP fluid and/or emfield variables to the FUENT face centers.\n\n");
    exit (EXIT_FAILURE);
  }
  for (cnt=1; cnt<argc; cnt++) if (strcmp(argv[cnt],"-tent")==0) TENTWEIGHT=1;
  if (TENTWEIGHT) {
    /* remove it before find_remaining_options complains about it */
    for (cnt=1; cnt<argc; cnt++) if (strcmp(argv[cnt],"-tent")==0) argv[cnt]=(char *)"";
  }
  /* -clamp, -pfloor, -tfloor. Same blank-out trick as -tent: process_flag_string
     is for required string flags, and find_remaining_options rejects anything it
     has not seen. A negative floor is refused here rather than baked into the
     UDF, where it would silently disable the check it is supposed to replace. */
  for (cnt=1; cnt<argc; cnt++) {
    if (strcmp(argv[cnt],"-clamp")==0) { CLAMP=1; argv[cnt]=(char *)""; }
  }
  for (cnt=1; cnt<argc-1; cnt++) {
    if (strcmp(argv[cnt],"-pfloor")==0) {
      pfloor_cl=atof(argv[cnt+1]);
      if (pfloor_cl<=0.0) {
        fprintf(stderr,"\n\n-pfloor must be a positive absolute pressure in Pa, got %s\n\n",argv[cnt+1]);
        exit (EXIT_FAILURE);
      }
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
    if (strcmp(argv[cnt],"-tfloor")==0) {
      tfloor_cl=atof(argv[cnt+1]);
      if (tfloor_cl<=0.0) {
        fprintf(stderr,"\n\n-tfloor must be a positive temperature in K, got %s\n\n",argv[cnt+1]);
        exit (EXIT_FAILURE);
      }
      argv[cnt]=(char *)""; argv[cnt+1]=(char *)"";
    }
  }
  if ((pfloor_cl>0.0 || tfloor_cl>0.0) && !CLAMP)
    fprintf(stdout,"NOTE: -pfloor/-tfloor were given without -clamp, so WARP_CLAMP is 0 and the\n"
                   "      floors are written into the UDF but not used. Add -clamp to enable them,\n"
                   "      or set WARP_CLAMP to 1 at the top of UDF_fluent2warp.c.\n");
  RET = find_remaining_options ( argc, argv, &options );
  if ( RET >= 1 ) {
    fprintf ( stderr, "\n\nThe following command line options could not be processed: %s\n\n", options );
    exit (EXIT_FAILURE);
  }
  else {
    /* Step 1 of the rendezvous, and it has to come before the file is opened:
       it renames the payload into the exchange directory and rewrites input to
       name it there. The handler is registered first so that every exit path
       from here on, including the exit(EXIT_FAILURE) two lines below, removes
       the staged file and stops CFDWARP rather than letting it re-read the
       previous exchange's answer. */
    if (EXCHANGE) {
      atexit(xchg_atexit);
      xchg_begin(input,400);
    }
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
    /* kept for the -transient clock check in xchg_finish(), which runs long
       after this file has been closed */
    XWARPTIME=tmp_time; XWARPDT=tmp_dt; XCLOCKOK=1;

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

  // Keep an untouched copy of the CFDWARP names: the reverse UDF has to reproduce
  // the vars_fluid string byte for byte, so it cannot use the renamed ones below.
  find_words_from_string(initvar_fluid_str_file, " ", &initvar_fluid_names_orig, &numvars_fluid_orig);

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
  double (*interpolated_vars)[numvars_fluid_file] = calloc(numfaces, sizeof(*interpolated_vars));
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
  double (*interpolated_emfield_vars)[numvars_emfield_file] = calloc(numfaces, sizeof(*interpolated_emfield_vars));
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
  // Every line carries the same number of columns: the coordinates, then all fluid variables, then all EMField variables
  numoutside=0;
  dist2worst=0.0;
  for (l_mesh = 0; l_mesh < numfaces; l_mesh++) {
    if (IS_3D) {
      fprintf(fluent_output,"%.16E %.16E %.16E ",x_mesh[l_mesh][0],x_mesh[l_mesh][1],x_mesh[l_mesh][2]);
    } else {
      fprintf(fluent_output,"%.16E %.16E ",x_mesh[l_mesh][0],x_mesh[l_mesh][1]);
    }
    if (weight[l_mesh]>1e-99){
      for (cnt=0; cnt<numvars_fluid_file; cnt++) {
        interpolated_vars[l_mesh][cnt]=interpolated_vars[l_mesh][cnt]/weight[l_mesh];
        fprintf(fluent_output,"%.16E ",interpolated_vars[l_mesh][cnt]);
      }
    } else {
      numoutside++;
      l_nearest=0;
      dist2min=-1.0;
      for (l_file=0; l_file<numnodes; l_file++){
        dist2=0.0;
        for (dim=0; dim<numdim_read; dim++) dist2+=sqr(x_mesh[l_mesh][dim]-x_file[l_file][dim]);
        if (dist2min<0.0 || dist2<dist2min) {
          dist2min=dist2;
          l_nearest=l_file;
        }
      }
      if (dist2min>dist2worst) dist2worst=dist2min;
      for (cnt=0; cnt<numvars_fluid_file; cnt++) {
        interpolated_vars[l_mesh][cnt]=initvar_fluid_file[l_nearest][cnt];
        fprintf(fluent_output,"%.16E ",interpolated_vars[l_mesh][cnt]);
      }
    }
    if (weight_emfield[l_mesh]>1e-99){
      for (cnt=0; cnt<numvars_emfield_file; cnt++) {
        interpolated_emfield_vars[l_mesh][cnt]=interpolated_emfield_vars[l_mesh][cnt]/weight_emfield[l_mesh];
        fprintf(fluent_output,"%.16E ",interpolated_emfield_vars[l_mesh][cnt]);
      }
    } else {
      l_nearest=0;
      dist2min=-1.0;
      for (l_file=0; l_file<numnodes_emfield; l_file++){
        dist2=0.0;
        for (dim=0; dim<numdim_emfield_read; dim++) dist2+=sqr(x_mesh[l_mesh][dim]-x_emfield_file[l_file][dim]);
        if (dist2min<0.0 || dist2<dist2min) {
          dist2min=dist2;
          l_nearest=l_file;
        }
      }
      for (cnt=0; cnt<numvars_emfield_file; cnt++) {
        interpolated_emfield_vars[l_mesh][cnt]=initvar_emfield_file[l_nearest][cnt];
        fprintf(fluent_output,"%.16E ",interpolated_emfield_vars[l_mesh][cnt]);
      }
    }
    fprintf(fluent_output,"\n");
  }
  fprintf ( stdout, "done.\n");
  if (numoutside>0) {
    fprintf ( stdout, "WARNING: %ld of the %ld FLUENT faces lie outside the CFDWARP domain and were\n"
                      "         given the values of the nearest CFDWARP node instead (worst gap %E m).\n"
                      "         Check that the FLUENT boundary is covered by the CFDWARP zone.\n",
              numoutside,numfaces,sqrt(dist2worst));
  }

  // Write Fluent UDF
  numvars_all=numvars_fluid_file+numvars_emfield_file;
  allvarnames=(char **)malloc(numvars_all*sizeof(char *));
  allcnames=(char **)malloc(numvars_all*sizeof(char *));
  for (cnt=0; cnt<numvars_fluid_file; cnt++) allvarnames[cnt]=initvar_fluid_names_file[cnt];
  for (cnt=0; cnt<numvars_emfield_file; cnt++) allvarnames[numvars_fluid_file+cnt]=initvar_emfield_names_file[cnt];
  for (cnt=0; cnt<numvars_all; cnt++) allcnames[cnt]=(char *)malloc(PROFILENAMELENGTH*sizeof(char));
  find_profile_names(allvarnames,numvars_all,allcnames);

  if (NOUDF) {
    fprintf ( stdout, "Not writing UDF.c: this call only converts. The UDF FLUENT compiled is\n"
                      "unchanged by an exchange, and rewriting it would advance its generation\n"
                      "stamp away from the one baked into libudf. Add -udf to write it anyway.\n");
  } else {
    fprintf ( stdout, "Writing FLUENT UDF...");
    write_udf_file("UDF.c",allvarnames,allcnames,numvars_all,numdim_read,output,numfaces);
    fprintf ( stdout, "done.\n");
  }

  // UDF that has FLUENT write a CFDWARP interpolation file
  xstation=0.0;
  for (dim=0; dim<numdim_read; dim++){ dx1mean[dim]=0.0; dx2mean[dim]=0.0; dx3mean[dim]=0.0; }
  for (l_file=0; l_file<numnodes; l_file++){
    xstation+=x_file[l_file][0];
    for (dim=0; dim<numdim_read; dim++){
      dx1mean[dim]+=dx1_file[l_file][dim];
      dx2mean[dim]+=dx2_file[l_file][dim];
      dx3mean[dim]+=dx3_file[l_file][dim];
    }
  }
  if (numnodes>0){
    xstation/=(double)numnodes;
    for (dim=0; dim<numdim_read; dim++){ dx1mean[dim]/=(double)numnodes; dx2mean[dim]/=(double)numnodes;
                                         dx3mean[dim]/=(double)numnodes; }
  }
  dxwarp=mesh_mean_spacing((const double *)x_file,numnodes,numdim_read);
  dxfluent=mesh_mean_spacing((const double *)x_mesh,numfaces,numdim_read);
  fprintf ( stdout, "Interface spacing: CFDWARP %E m, FLUENT %E m.\n",dxwarp,dxfluent);
  if (dxfluent>=dxwarp)
    fprintf ( stdout, "WARNING: the FLUENT interface mesh is not finer than the CFDWARP one, so no\n"
                      "         dx1/dx2 can both cover every CFDWARP node and stay clear of its\n"
                      "         neighbours. The zone read is then a smoother and the CFDWARP\n"
                      "         residual will not settle; judge the coupling on the interface\n"
                      "         movement per exchange instead.\n");
  iPfile=-1; iTfile=-1;
  for (cnt=0; cnt<numvars_fluid_file; cnt++) {
    if (strcmp(initvar_fluid_names_orig[cnt],"P")==0) iPfile=cnt;
    if (strcmp(initvar_fluid_names_orig[cnt],"T")==0) iTfile=cnt;
  }
  pfloor_use=1.0e-03; tfloor_use=1.0e-03; pmean_file=-1.0;
  if (iPfile>=0 && numfaces>0) {
    pfloor_use=interpolated_vars[0][iPfile];
    pmean_file=0.0;
    for (l_mesh=0; l_mesh<numfaces; l_mesh++) {
      if (interpolated_vars[l_mesh][iPfile]<pfloor_use) pfloor_use=interpolated_vars[l_mesh][iPfile];
      pmean_file+=interpolated_vars[l_mesh][iPfile];
    }
    pmean_file/=(double)numfaces;
    if (pfloor_use<=0.0) {
      fprintf ( stdout, "WARNING: the smallest pressure in the CFDWARP file is %E Pa, which cannot be\n"
                        "         used as a floor. WARP_PFLOOR falls back to %E Pa; set it with\n"
                        "         -pfloor if you mean to use -clamp at all.\n", pfloor_use, 1.0e-03);
      pfloor_use=1.0e-03;
    }
  }
  if (iTfile>=0 && numfaces>0) {
    tfloor_use=interpolated_vars[0][iTfile];
    for (l_mesh=0; l_mesh<numfaces; l_mesh++)
      if (interpolated_vars[l_mesh][iTfile]<tfloor_use) tfloor_use=interpolated_vars[l_mesh][iTfile];
    if (tfloor_use<=0.0) tfloor_use=1.0e-03;
  }
  if (pfloor_cl>0.0) pfloor_use=pfloor_cl;
  if (tfloor_cl>0.0) tfloor_use=tfloor_cl;

  if (NOUDF) {
    fprintf ( stdout, "Not writing UDF_fluent2warp.c either, same reason.\n");
  } else {
    fprintf ( stdout, "Writing FLUENT-to-CFDWARP UDF...");
    write_fluent2warp_udf_file("UDF_fluent2warp.c",initvar_fluid_str_file,initvar_fluid_names_orig,
                               numvars_fluid_file,numdim_read,numspec_read,xstation,dx1mean,dx2mean,
                               dx3mean,
                               numflux_read,windowis,windowie,iter,effiter_U,effiter_R,CFL,
                               dxwarp,dxfluent,numfaces,CLAMP,pfloor_use,tfloor_use,pmean_file);
    fprintf ( stdout, "done.\n");
  }
  if (CLAMP && !NOUDF)
    fprintf ( stdout, "WARP_CLAMP is 1: a face whose absolute P is below %E Pa or whose T is below\n"
                      "%E K will be RAISED to %E Pa / %E K and sent, not refused. CFDWARP stores\n"
                      "the floor as a measurement. The exchange reports the count every time.\n",
              1.0e-03, 1.0e-03, pfloor_use, tfloor_use);

  fprintf ( stdout, "\nFLUENT profile names (pick these from the boundary condition drop-down lists):\n");
  for (cnt=0; cnt<numvars_all; cnt++) {
    fprintf ( stdout, "  %-20s -> %s_profile\n", allvarnames[cnt], allcnames[cnt]);
  }
  fprintf ( stdout, "\nThe interpolated data file, UDF.c and UDF_fluent2warp.c were written to the\n"
                    "directory warp2fluent was RUN from, not to the one holding the executable.\n"
                    "Both absolute paths baked into UDF.c are frozen at that directory, so run\n"
                    "warp2fluent in the folder FLUENT starts from and load the UDF.c it wrote\n"
                    "there. Loading a UDF.c generated in some other folder is not caught by\n"
                    "FLUENT and is reported only by the banner warp_load() prints.\n");

  for (cnt=0; cnt<numvars_all; cnt++) free(allcnames[cnt]);
  free(allcnames);
  free(allvarnames);

  fclose(datafile);
  fclose(meshdata);
  fclose(fluent_output);

  if (EXCHANGE) xchg_finish();

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
  for (int i = 0; i < numvars_fluid_orig; i++) free(initvar_fluid_names_orig[i]);
  free(initvar_fluid_names_orig);
  for (int i = 0; i < numvars_emfield_file; i++) free(initvar_emfield_names_file[i]);
  free(initvar_emfield_names_file);

  return(EXIT_SUCCESS);
}



