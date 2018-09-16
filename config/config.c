#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <termios.h>
#include <unistd.h>


#define FONT_BOLD "\033[1;1m"
#define FONT_STANDARD "\033[0m"

#define makefileheader "../.makefile-header"
#ifndef TRUE
  #define TRUE    1
  #define FALSE   0
#endif
#define OPENMP TRUE
#define numkeys 9
#define KEY_RET 10
#define EOS 0
#define EOL '\n'

#define CO_CC 0
#define CO_GCC 1
#define CO_CCC 2
#define CO_ICC 3
#define CO_MPICC 4
#define CO_num 5


#define mod(a,b) ((((a)%(b))+(b))%(b))

#if (OPENMP==TRUE)
  #define THREADTYPE_max 3
#else
  #define THREADTYPE_max 2
#endif
static char CO_str[CO_num][20]=
  {"cc","gcc","ccc","icc","mpicc"};

typedef unsigned char bool;

#define MAX_STR_LENGTH 5000
typedef char type_str[3*MAX_STR_LENGTH];

typedef int keys_t[numkeys];




/* insert str1 into str2 at the position pos; 
   make sure *str2 has enough memory allocated */
static char *strins(char *str1, char *str2,  long pos){
  long len1,len2,i;
  len1=(long)strlen(str1);
  len2=(long)strlen(str2);

  for (i=len2; i>=pos; i--) (str2)[i+len1]=(str2)[i];
  for (i=0; i<len1; i++) (str2)[pos+i]=str1[i];
  (str2)[len2+len1]=EOS;
  return str2;
}


static char *strupr(char *str){
  long cnt;
  for (cnt=0; str[cnt] != '\0'; cnt++){
    if (str[cnt] >= 'a' && str[cnt] <= 'z'){
      str[cnt] = str[cnt]-32;
    } 
  }
  str[cnt] = '\0';
  return(str);
}


static char *strnounderscore(char *str){
  long cnt;
  for(cnt=0; str[cnt] != '\0'; cnt++){
    if (str[cnt]=='_') str[cnt]=' ';
  }
  return(str);
}


/* add line breaks without breaking words with width the maximum number of characters per line 
   and indent the number of indented characters (either negative or positive)  */
static char *strwrpind(char *str, int width, int indent){
  long cnt,cnt2,cntbreak;
  bool CONTINUE;
  static char whitespace[2];
  strcpy(whitespace," ");

  if (indent>0){
    for (cnt=0; cnt<indent; cnt++) strins(whitespace,str,0);
  }

  CONTINUE=TRUE;
  cntbreak=0;
  cnt=0;
  do {
    cnt++;
    if (str[cnt]=='\n') cntbreak=cnt;
    if (cnt-cntbreak>width-1){
      cntbreak=cnt;
      do {
        cntbreak--;
      } while(str[cntbreak]!=' ' && str[cntbreak]!='-' && cntbreak>0);
      if (cntbreak>0){
        if (str[cntbreak]=='-') {
          strins("\n",str,cntbreak+1);
          cntbreak++;
        }
        str[cntbreak]='\n';
        cnt++;
        if (indent<0){
          // do the hang indent here
          for (cnt2=0; cnt2<-indent; cnt2++){
            strins(whitespace,str,cntbreak+1);
          }
        }
      } else {
        // problem breaking line..
        CONTINUE=FALSE;
      }
    }        
  } while (CONTINUE && cnt<strlen(str)-1); 
  return str;
}


/* - get_key function: will grab one key from stdin
   - timeout is the key waiting timeout is in tenths of a second;
     if negative, wait forever
 */
static int get_key(long timeout){
  int key;
  struct termios term_orig, term_new;

  tcgetattr(0, &term_orig);
  term_new = term_orig;
  term_new.c_lflag &= ~(ICANON|ECHO);
  if (timeout>=0) {
    term_new.c_cc[VTIME] = timeout;
    term_new.c_cc[VMIN] = timeout ? 0 : 1;
  } else {
    term_new.c_cc[VTIME] = 0;
    term_new.c_cc[VMIN] = 1;
  }
  tcsetattr(0, TCSANOW, &term_new);
  key = getc(stdin);
  tcsetattr(0, TCSANOW, &term_orig);
  return(key);
}


/* - get_keys subroutine: will grab one key or a set of keys
     whether the first key (key[0]) is the escape character (27).
   - total number of keys read corresponds to x where keys[x]=-1 */
static void get_keys(keys_t keys){
  int cnt;
  for (cnt=0; cnt<numkeys; cnt++){
    keys[cnt]=-1;
  }
  keys[0]=get_key(-1);
  if (keys[0]==27) {
    keys[1]=get_key(1);
    if (keys[1]==91) {
      cnt=1;
      do {
        cnt++;
        keys[cnt]=get_key(-1);
      } while (keys[cnt]<64);
    }
    if (keys[1]==79) {
      keys[2]=get_key(-1);
    }
  }
}


/* returns 1 on success, 0 otherwise
   seems to have problem if DISTMPI is needed but DDISTMPI is in the file*/
static int find_str_in_file(FILE *file, char *str)  {
  long len, i, c;

  if ((file==NULL) || (str==NULL)) return 0;
  len = (long)strlen(str);
  i = 0;

  while(i != len){
    c=fgetc(file);
    if (c == EOF) return 0;
    if (str[i] == (char)c) {
      i++;
    } else {
      i = 0;
    }
  }

  return 1;
}


/* returns 1 on success, 0 otherwise */
static int find_str_in_file_named(char *filename, char *string){
  FILE *file;

  file = fopen(filename, "r");
  if (file == NULL) {
    return 0;
  }

  if (find_str_in_file(file, string) == 1) {
    fclose(file);
    return 1;
  }

  fclose(file);
  return 0;
}


/* returns 1 on success, 0 otherwise */
static int find_num_after_str_in_file_named(char *filename, char *string, long *number){
  FILE *file;

  file = fopen(filename, "r");
  if (file == NULL) {
    return 0;
  }

  if (find_str_in_file(file, string) == 0) {
    fclose(file);
    return 0;
  }

  if (fscanf(file, "%ld", number) != 1) {
    fclose(file);
    return 0;
  }

  fclose(file);
  return 1;
}


static void get_one_number(long *answer, bool *ESCFLAG, bool *ENTERFLAG,
                  bool *ERRORFLAG){
  keys_t keys;
  get_keys(keys);
  *ENTERFLAG=FALSE;
  *ERRORFLAG=FALSE;
  if (keys[0]==10) {
    *ENTERFLAG=TRUE;
  }
  if (!(*ENTERFLAG)) {
    if ((keys[0]>=48) && (keys[0]<=57)) {
      *answer=keys[0]-48;
    } else {
      *ERRORFLAG=TRUE;
    }
  }
}


static void get_one_letter(long *answer, bool *ESCFLAG, bool *ENTERFLAG,
                  bool *ERRORFLAG, bool *HELPFLAG){
  keys_t keys;
  get_keys(keys);
  *ENTERFLAG=FALSE;
  *ERRORFLAG=FALSE;
  *HELPFLAG=FALSE;
  *answer=-9999;
  if (keys[0]==10) {
    *ENTERFLAG=TRUE;
  }
  if (keys[0]=='?') {
    *HELPFLAG=TRUE;
  }
  if (!(*ENTERFLAG)) {
    if ((keys[0]>='a') && (keys[0]<='z')) {
      *answer=keys[0]-'a'+1;
    } else {
      if ((keys[0]>='A') && (keys[0]<='Z')) {
        *answer=keys[0]-'A'+'z'-'a'+2;
      } else {
        *ERRORFLAG=TRUE;
      }
    }
  }
}


static void get_yes_no(bool *yes, bool *ESCFLAG, bool *ENTERFLAG,
              bool *ERRORFLAG){
  keys_t keys;
  get_keys(keys);
  *ENTERFLAG=FALSE;
  *ERRORFLAG=FALSE;
  if (keys[0]==10) {
    *ENTERFLAG=TRUE;
  }
  if (!(*ENTERFLAG)) {
    if ((keys[0]=='y') || (keys[0]=='n')) {
      if (keys[0]=='y') {
        *yes=TRUE;
      } else {
        *yes=FALSE;
      }
    } else {
      *ERRORFLAG=TRUE;
    }
  }
}


static void read_line_from_config_file(FILE *configfile, type_str str1, type_str str2){
  long cnt;
  char chr;
  bool INSTRING;

  fscanf(configfile, "%s",str1);
  if (strcmp(str1,"END")==0) {
    str2[0]=EOS;
  } else {
    if (strcmp(str1,"message")==0 || strcmp(str1,"help")==0) {
      cnt=0;
      INSTRING=FALSE;
      do {
        chr=fgetc(configfile);
        if (chr!=' ') INSTRING=TRUE;
        if (INSTRING) {
          str2[cnt]=chr;
          cnt++;
        }
      } while (chr!=EOL && cnt<MAX_STR_LENGTH);
      str2[cnt-1]=EOS;
    } else {
      fscanf(configfile, "%s*[^\n]",str2);
    }
  }
}



static void read_config_file(type_str dirname, type_str current, type_str message,
                     type_str **branch, type_str **option, type_str **optioncomment, long *numoption,
                     long *numbranch, type_str help){
  FILE *configfile;
  type_str filename,str1,str2;
  bool CONTINUE;
  long cntoption,cntbranch;

  strcpy(filename,dirname);
  strcat(filename,"/.config");
  configfile=fopen(filename, "r");
  *numoption=0;
  *numbranch=0;
  CONTINUE=TRUE;
  while (CONTINUE) {
    read_line_from_config_file(configfile,str1,str2);
    if (strcmp(str1,"END")==0) CONTINUE=FALSE;
    if (strcmp(str1,"option")==0) (*numoption)++;
#ifdef EXPERIMENTAL
    if (strcmp(str1,"optionX")==0) (*numoption)++;  //experimental module
    if (strcmp(str1,"optionPX")==0) (*numoption)++;  //experimental and proprietary module
#endif
    if (strcmp(str1,"optionP")==0) (*numoption)++;  //proprietary module
    if (strcmp(str1,"branch")==0) (*numbranch)++;
  }
  fclose(configfile);

  *branch=(type_str *) malloc((*numbranch)*sizeof(type_str));
  *option=(type_str *) malloc((*numoption)*sizeof(type_str));
  *optioncomment=(type_str *) malloc((*numoption)*sizeof(type_str));
  configfile=fopen(filename, "r");
  cntoption=0;
  cntbranch=0;
  CONTINUE=TRUE;
  strcpy(help,"\0");
  while (CONTINUE) {
    read_line_from_config_file(configfile,str1,str2);
    if (strcmp(str1,"END")==0) CONTINUE=FALSE;
    if (strcmp(str1,"current")==0) strcpy(current,str2);
    if (strcmp(str1,"help")==0) strcpy(help,str2);
    if (strcmp(str1,"message")==0) strcpy(message,str2);
    if (strcmp(str1,"option")==0) {
      strcpy((*optioncomment)[cntoption],"");
      strcpy((*option)[cntoption],str2);
      cntoption++;
    }
#ifdef EXPERIMENTAL
    if (strcmp(str1,"optionX")==0) {
      strcpy((*optioncomment)[cntoption],"[EXPERIMENTAL]");
      strcpy((*option)[cntoption],str2);
      cntoption++;
    }
    if (strcmp(str1,"optionPX")==0) {
      strcpy((*optioncomment)[cntoption],"[PROPRIETARY] [EXPERIMENTAL]");
      strcpy((*option)[cntoption],str2);
      cntoption++;
    }
#endif
    if (strcmp(str1,"optionP")==0) {
      strcpy((*optioncomment)[cntoption],"[PROPRIETARY]");
      strcpy((*option)[cntoption],str2);
      cntoption++;
    }
    if (strcmp(str1,"branch")==0) {
      strcpy((*branch)[cntbranch],str2);
      cntbranch++;
    }

  }
  fclose(configfile);
}


static void find_default_option(type_str current, type_str *option, long numoption,
                       long *defaultoption){
  long cnt;
  *defaultoption=-1; /* initialize to a negative value, which is not a valid option */
  for (cnt=0; cnt<numoption; cnt++){
    if (strcmp(current,option[cnt])==0){
      *defaultoption=cnt+1;
    }
  }
}


static void process_config_files(FILE *scriptfile_applyconfig, FILE *scriptfile_removeproprietary, type_str dirname);


static void print_section_title(type_str title){
  printf("\n  %s%s%s\n\n",FONT_BOLD,title,FONT_STANDARD);  
}


static void handle_options(FILE *scriptfile_applyconfig, FILE *scriptfile_removeproprietary, type_str dirname, type_str message, type_str *option, type_str *optioncomment, long numoption, type_str help){
  long cnt,answer;
  type_str current,current_formatted,newdirname,optionfilename,linkfilename;
  bool HELPGIVEN,ERRORFLAG,ENTERFLAG,ESCFLAG,HELPFLAG;
#ifndef PROPRIETARY
  bool HELPGIVEN_PROPRIETARY;
#endif
  type_str *garb_branch,*garb_option,*garb_optioncomment,garb_message,garb_help,help_proprietary;
  long garb_numoption,letteroption;
  long garb_numbranch;
  static char strtmp[MAX_STR_LENGTH];
  sprintf(help_proprietary,"%s","This module is not included in this distribution. Please contact the CFDWARP maintainer for details on how to obtain this module.");
  if (numoption>1) {
    /* newdirname=filename+"/.active" */
    strcpy(newdirname,dirname);
    strcat(newdirname,"/.active");
    read_config_file(newdirname, current, garb_message,
                    &garb_branch, &garb_option, &garb_optioncomment, &garb_numoption,
                    &garb_numbranch,garb_help);
    print_section_title(strupr(message));
    for (cnt=0; cnt<numoption; cnt++){
      strcpy(strtmp,option[cnt]);
      letteroption='a'+cnt;
      if (letteroption>'z') letteroption='A'+letteroption-'z'-1;
      printf("  %c. %s %s\n",(char)letteroption,strnounderscore(strtmp),optioncomment[cnt]);
    }
    if (help[0]!=EOS) printf("  ?  help\n");
    strcpy(current_formatted,current);
    strnounderscore(current_formatted);
    printf("\n  %s : ",current_formatted);
    HELPGIVEN=FALSE;

#ifndef PROPRIETARY
    HELPGIVEN_PROPRIETARY=FALSE;
#endif

    do {
      get_one_letter(&answer,&ESCFLAG,&ENTERFLAG,&ERRORFLAG,&HELPFLAG);
      if (ENTERFLAG) {
        find_default_option(current,option,numoption,&answer);
      }
      if (HELPFLAG && help[0]!=EOS && !HELPGIVEN) {
        printf("?\n\n  %s\n\n  %s : ",strwrpind(help, 50, -2),current_formatted);
        HELPGIVEN=TRUE;
      }
#ifndef PROPRIETARY
      if (!HELPFLAG && (answer>0 && answer<=numoption) ) {
        if (strstr(optioncomment[answer-1],"[PROPRIETARY]")) {
          strcpy(strtmp,option[answer-1]);
          if (!HELPGIVEN_PROPRIETARY) printf("%s\n\n  %s\n\n  %s : ",strnounderscore(strtmp),strwrpind(help_proprietary, 50, -2),current_formatted);
          HELPGIVEN_PROPRIETARY=TRUE;
          answer=0;
        }
      }
#endif
    } while (!(answer>0 && answer<=numoption));
    answer=answer-1;
    strcpy(strtmp,option[answer]);
    printf("%s\n\n",strnounderscore(strtmp));
    strcpy(optionfilename,"_");
    strcat(optionfilename,option[answer]);
    strcpy(linkfilename,dirname);
    strcat(linkfilename,"/.active");

    fprintf(scriptfile_applyconfig,"rm -f %s \n",linkfilename);
    fprintf(scriptfile_applyconfig,"ln -s %s %s \n",optionfilename,linkfilename);

    for (cnt=0; cnt<numoption; cnt++){
      if (strstr(optioncomment[cnt],"[PROPRIETARY]")) {
        /* write directives to remove proprietary directory and substitute it by the _blank directory */
        strcpy(optionfilename,dirname);
        strcat(optionfilename,"/_");
        strcat(optionfilename,option[cnt]);
        fprintf(scriptfile_removeproprietary,"echo \"Removing proprietary directory %s\"\n",optionfilename);
        fprintf(scriptfile_removeproprietary,"rm -rf %s \n",optionfilename);
        fprintf(scriptfile_removeproprietary,"cp -a _proprietary %s \n",optionfilename);
        strcat(optionfilename,"/.config");
        fprintf(scriptfile_removeproprietary,"( printf \"current %s\\nEND\" ) > %s \n",option[cnt],optionfilename);
      
      }
    }
  }
}


static void handle_branches(FILE *scriptfile_applyconfig, FILE  *scriptfile_removeproprietary, type_str dirname, type_str *branch, long numbranch){
  long cnt;
  type_str newdirname;

  for (cnt=0; cnt<numbranch; cnt++){
    /* newdirname=dirname+"/"+branch[cnt] */
    strcpy(newdirname,dirname);
    strcat(newdirname,"/");
    strcat(newdirname,branch[cnt]);
    process_config_files(scriptfile_applyconfig, scriptfile_removeproprietary, newdirname);
  }
}


static void process_config_files(FILE *scriptfile_applyconfig, FILE *scriptfile_removeproprietary, type_str dirname){
  type_str current,message,help;
  type_str *branch,*option,*optioncomment;
  long numoption,numbranch;

  read_config_file(dirname, current, message,
                  &branch, &option, &optioncomment, &numoption,&numbranch, help);
  handle_options(scriptfile_applyconfig, scriptfile_removeproprietary, dirname, message, option, optioncomment, numoption, help);
  handle_branches(scriptfile_applyconfig, scriptfile_removeproprietary, dirname, branch, numbranch);
}


static void ask_yesno_question(bool answer_def, char *filestring, char *question, char *options, bool *answer){
  bool yes,ENTERFLAG,ESCFLAG,ERRORFLAG;
  if (find_str_in_file_named(makefileheader, filestring)==1) {
    answer_def=mod(answer_def+1,2);
  }
  print_section_title(question);
  if (answer_def) {
    printf("%s\n\n  yes : ",options);
  } else {
    printf("%s\n\n  no : ",options);
  }
  do {
    yes=answer_def;
    get_yes_no(&yes, &ESCFLAG, &ENTERFLAG, &ERRORFLAG);
  } while (ERRORFLAG);
  if (yes) {
    *answer=TRUE;
    printf("yes\n\n");
  } else {
    *answer=FALSE;
    printf("no\n\n");
  }
}


static void ask_number_question(long numbermin, long numbermax, long answer_def, char *filestring, char *question, char *options, long *answer){
  bool ENTERFLAG,ESCFLAG,ERRORFLAG;
  long tmp;
  if (find_num_after_str_in_file_named(
    makefileheader, filestring, &tmp)==0){
    tmp=answer_def;
  }
  print_section_title(question);
  printf("%s\n\n  %ld : ",options,tmp);
  do {
    *answer=tmp;
    get_one_number(answer, &ESCFLAG, &ENTERFLAG, &ERRORFLAG);
    if (*answer<numbermin || *answer>numbermax) ERRORFLAG=TRUE;
  } while (ERRORFLAG);
  printf("%ld\n\n",*answer);
}


static void get_makefile_vars(long *nd, long *numopt,
                       bool *DEBUG, bool *DEBUGGER, bool *PROFIL,
                       long *THREADTYPE, 
                       long *COMPILER, bool *STATIC, bool *TEST, bool *M32, bool *SANITIZE){
  bool ESCFLAG,ENTERFLAG,ERRORFLAG;
  long nd_def,cnt,COMPILER_def,THREADTYPE_def;

  COMPILER_def=CO_CC;
  for (cnt=0; cnt<CO_num; cnt++)
    if (find_str_in_file_named(makefileheader, CO_str[cnt])==1) COMPILER_def=cnt;
  printf("\n");
  print_section_title("COMPILER");
  for (cnt=0; cnt<CO_num; cnt++) printf("  %ld. %s\n",cnt,CO_str[cnt]);
  printf("\n  %s : ",CO_str[COMPILER_def]);
  do {
    *COMPILER=COMPILER_def;
    get_one_number(COMPILER, &ESCFLAG, &ENTERFLAG, &ERRORFLAG);
  } while ((*COMPILER<0) || (*COMPILER>=CO_num) || (ERRORFLAG));
  printf("%s\n\n",CO_str[*COMPILER]);

  print_section_title("THREADING");
  printf("  0. none\n");
  printf("  1. POSIX\n");
  printf("  2. POSIX MULTIZONE\n");
  if (OPENMP){
    printf("  3. OPENMP\n");
  }
  printf("\n  ");
  if (find_str_in_file_named(makefileheader, "POSIXTHREADS")==1) {
    printf("POSIX : ");
    THREADTYPE_def=1;
  } else {
    if (find_str_in_file_named(makefileheader, "ZONETHREADS")==1) {
      printf("POSIX MULTIZONE : ");
      THREADTYPE_def=2;
    } else {
      if (OPENMP && find_str_in_file_named(makefileheader, "OPENMPTHREADS")==1) {
        printf("OPENMP : ");
        THREADTYPE_def=3;
      } else {
        printf("none : ");
        THREADTYPE_def=0;
      }
    }
  }
  do {
    *THREADTYPE=THREADTYPE_def;
    get_one_number(THREADTYPE, &ESCFLAG, &ENTERFLAG, &ERRORFLAG);
  } while ((*THREADTYPE<0) || (*THREADTYPE>THREADTYPE_max) || (ERRORFLAG));
  if (*THREADTYPE==0) printf("none\n\n");
  if (*THREADTYPE==1) printf("POSIX\n\n");
  if (*THREADTYPE==2) printf("POSIX MULTIZONE\n\n");
  if (*THREADTYPE==3) printf("OPENMP\n\n");

  ask_number_question(0, 3, 3,  "-O", "COMPILER OPTIMIZATIONS",
                                         "  0. none\n"
                                         "  1. basic optimizations, valgrind-compatible\n"
                                         "  2. more optimizations, well-tested\n"
                                         "  3. further optimizations (loop unrolling, function inlining)"
                                         ,numopt);
  ask_yesno_question(FALSE,"-static",  "STATICALLY LINKED EXECUTABLE",
                                         "  n. no\n"
                                         "  y. yes",STATIC);
  ask_yesno_question(TRUE, "-DNDEBUG", "ASSERTIONS",
                                         "  n. no\n"
                                         "  y. yes",DEBUG);
  ask_yesno_question(FALSE, "-fsanitize", "SANITIZATION",
                                         "  n. no\n"
                                         "  y. yes",SANITIZE);
  ask_yesno_question(FALSE,"-g", "DEBUGGING SYMBOLS",
                                         "  n. no\n"
                                         "  y. yes",DEBUGGER);
  ask_yesno_question(FALSE,"-pg", "CPU PROFILING",
                                         "  n. no\n"
                                         "  y. yes",PROFIL);
  ask_yesno_question(FALSE, "-DTEST", "TEST MODE",
                                         "  n. no\n"
                                         "  y. yes",TEST);
  ask_yesno_question(FALSE, "-m32", "32-BIT-ARCHITECTURE EXECUTABLE",
                                         "  n. no\n"
                                         "  y. yes",M32);

  nd_def=2;
  if (find_str_in_file_named(makefileheader, "-D_3D")==1) nd_def=3;
  if (find_str_in_file_named(makefileheader, "-D_2D")==1) nd_def=2;
  print_section_title("NUMBER OF DIMENSIONS");
  printf("  2. 2D\n"
         "  3. 3D\n\n"
         "  %ldD : ",nd_def);
  do {
    *nd=nd_def;
    get_one_number(nd, &ESCFLAG, &ENTERFLAG, &ERRORFLAG);
    if ((*nd<2) || (*nd>3)) ERRORFLAG=TRUE;
  } while (ERRORFLAG);
  printf("%ldD\n\n",*nd);

}


static void output_makefileheader(FILE *scriptfile_applyconfig, long nd, long numopt,
                           bool DEBUG, bool DEBUGGER, bool PROFIL,
                           long THREADTYPE, 
                           long COMPILER, bool STATIC, bool TEST, bool M32, bool SANITIZE){
  fprintf(scriptfile_applyconfig,"(\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## This makefile header has been automatically generated by        ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## the config program config/config.c through 'make config'.       ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## It can be changed manually prior to compiling the code.         ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## Note that should another 'make config' take place, your         ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## changes will be overridden.                                     ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## CC        --> C compiler                                        ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## CFLAGS    --> compiler flags for any language                   ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## CFLAGSCFD --> compiler flags when compiling CFD code            ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## L         --> linker                                            ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## LFLAGS    --> linker flags                                      ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"MAKEFLAGS += --stop\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"MAKEFLAGS += --no-print-directory\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"SHELL = /bin/sh -e\\n\"\n");
  if (COMPILER==CO_ICC) {
    fprintf(scriptfile_applyconfig,"  printf \"CFLAGS = -O"); 
  } else {
    fprintf(scriptfile_applyconfig,"  printf \"CFLAGS = -O%ld",numopt);
  }
  if (COMPILER==CO_GCC || COMPILER==CO_ICC || COMPILER==CO_MPICC) fprintf(scriptfile_applyconfig, " -Wall -funroll-all-loops -Wfatal-errors");
  if (COMPILER==CO_CCC) fprintf(scriptfile_applyconfig, " -ieee");
  if (COMPILER==CO_MPICC) fprintf(scriptfile_applyconfig, " -DDISTMPI");
  if (PROFIL) fprintf(scriptfile_applyconfig," -pg");
  if (DEBUGGER) fprintf(scriptfile_applyconfig," -g");
  if (DEBUGGER && (COMPILER==CO_GCC || COMPILER==CO_MPICC)) fprintf(scriptfile_applyconfig," -fno-omit-frame-pointer"); 
  if (!DEBUG) fprintf(scriptfile_applyconfig," -DNDEBUG");
  if (TEST) fprintf(scriptfile_applyconfig," -DTEST");
  if (M32) fprintf(scriptfile_applyconfig," -m32 -DM32");
  if (SANITIZE) fprintf(scriptfile_applyconfig," -fsanitize=address -fsanitize=undefined -fsanitize=bool -fsanitize=integer-divide-by-zero -fsanitize=null -fsanitize=vla-bound -fsanitize=bounds -fsanitize=signed-integer-overflow");

  if (getlogin()!=NULL){    
    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);

    fprintf (scriptfile_applyconfig," -D_USERNAME=\\\"\\\\\\\"%s\\\\\\\"\\\"",getlogin());
    fprintf (scriptfile_applyconfig," -D_HOSTNAME=\\\"\\\\\\\"%s\\\\\\\"\\\"",hostname);
  } 

  if (THREADTYPE==1) fprintf(scriptfile_applyconfig," -DPOSIXTHREADS -D_REENTRANT");
  if (THREADTYPE==2) fprintf(scriptfile_applyconfig," -DZONETHREADS -D_REENTRANT");
  if (THREADTYPE==3) {
    if (COMPILER==CO_ICC){
      fprintf(scriptfile_applyconfig," -DOPENMPTHREADS -openmp");
    } else {
      fprintf(scriptfile_applyconfig," -DOPENMPTHREADS -fopenmp");
    }
  }
  fprintf(scriptfile_applyconfig,"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"CC = %s\\n\"\n",CO_str[COMPILER]);
  fprintf(scriptfile_applyconfig,"  printf \"CFLAGSCFD = -D_%ldD ",nd);
  fprintf(scriptfile_applyconfig,"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"L = %s\\n\"\n",CO_str[COMPILER]);
  fprintf(scriptfile_applyconfig,"  printf \"LFLAGS =");
  if (THREADTYPE==3) {
    if (COMPILER==CO_ICC){
      fprintf(scriptfile_applyconfig," -openmp");
    } else {
      fprintf(scriptfile_applyconfig," -fopenmp");
    }
//    if (COMPILER==CO_GCC) fprintf(scriptfile_applyconfig," -s");
  }
  if (STATIC) fprintf(scriptfile_applyconfig," -static");
  if (DEBUGGER) fprintf(scriptfile_applyconfig," -g");
  if (DEBUGGER && (COMPILER==CO_GCC || COMPILER==CO_MPICC)) fprintf(scriptfile_applyconfig," -fno-omit-frame-pointer"); 
  if (PROFIL) fprintf(scriptfile_applyconfig," -pg");
  if (!DEBUG) fprintf(scriptfile_applyconfig," -DNDEBUG");
  if (TEST) fprintf(scriptfile_applyconfig," -DTEST");
  if (M32) fprintf(scriptfile_applyconfig," -m32 -DM32");
  if (SANITIZE) fprintf(scriptfile_applyconfig," -fsanitize=address -fsanitize=undefined -fsanitize=bool -fsanitize=integer-divide-by-zero -fsanitize=null -fsanitize=vla-bound -fsanitize=bounds -fsanitize=signed-integer-overflow");
  fprintf(scriptfile_applyconfig,"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## Implicit Rules                                                  ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \".SUFFIXES: .c \\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \".c.o:\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"	\\$(CC) -c \\$(CFLAGS) \\$(CCFLAGSLOCAL)  \\$<\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"## Build Rules                                                     ##\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"#####################################################################\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"default:\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"	make all\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"showfiles:\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"	@echo \\$(SOURCES) \\$(HEADERS) Makefile\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"clean:\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"	-rm -f *.o *.a *.so *.bak *BAK *~ *%%%% #*\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"	-rm -f \\$(TARGETS)\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"	-rm -f core\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"\\n\"\n");
  fprintf(scriptfile_applyconfig,"  printf \"cleanall: clean\\n\"\n");
  fprintf(scriptfile_applyconfig,") > %s\n",makefileheader);
}


int main(){
  long nd,numopt;
  bool DEBUG,DEBUGGER,PROFIL,STATIC,TEST,M32,SANITIZE;
  long THREADTYPE;
  long COMPILER;
  FILE *scriptfile_applyconfig;
  FILE *scriptfile_removeproprietary;
  type_str dirname;

  scriptfile_applyconfig=fopen("applyconfig.sh", "w");
  scriptfile_removeproprietary=fopen("removeproprietary.sh_tmp", "w");
  fprintf(scriptfile_applyconfig,"#!/bin/sh\n");
  fprintf(scriptfile_removeproprietary,"#!/bin/sh\n");

  get_makefile_vars(&nd,&numopt,&DEBUG,&DEBUGGER,&PROFIL,
                    &THREADTYPE,&COMPILER,&STATIC,&TEST,&M32,&SANITIZE);
  strcpy(dirname,"../.");
  process_config_files(scriptfile_applyconfig,scriptfile_removeproprietary,dirname);    
  output_makefileheader(scriptfile_applyconfig,nd,numopt,DEBUG,DEBUGGER,PROFIL,
                     THREADTYPE,COMPILER,STATIC,TEST,M32,SANITIZE);

  fclose(scriptfile_applyconfig);
  fclose(scriptfile_removeproprietary);
  if (rename("removeproprietary.sh_tmp","removeproprietary.sh")!=0){
    fprintf(stderr,"\n\nProblem encountered while trying to rename removeproprietary.sh_tmp to removeproprietary.sh.\n\n");
  }
  return(EXIT_SUCCESS);
}

