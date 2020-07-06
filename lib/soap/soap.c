// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-1999 Bernard Parent
Copyright 2020      Prasanna Thoguluva Rajendran

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



#include "soap.h"
#include "printf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#ifdef DISTMPI
  #include "mpi.h"
#endif

#define EOS 0
#define pi 3.14159265358979323846
#define sqr(a)   ((a)*(a))
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))
#define rad(a)     (((a)*pi/180.0e0))
#define deg(a)     (((a)/pi*180.0e0))
#ifndef round
//#define round(a) (floor((a)+0.5e0))
#define round(a) (a<0?ceil((a)-0.5):floor((a)+0.5))
#endif
//#define longfromdouble(a)  ((long)(a+0.5))
#define longfromdouble(a) ((a)>=0?(long)((a)+0.5):(long)((a)-0.5))

#define DOUBLEFORMAT  "%11.11E"
/* the latter will give 12 significant numbers
   when performing calculations */

/* logical operators */
#define GT   21
#define GEQ  22
#define LT   23
#define LEQ  24
#define EQ   25
#define AND  26
#define OR   27
#define NOT  28
#define NEQ  29

#define maxnumlen 30
/* 30 chars are enough to store 15 significant numbers */


/* this function calls vfprintf with the same arguments
   as it is given and exits.*/
int SOAP_fatal_error(SOAP_codex_t *codex, const char *formatstr, ...){
  va_list ap;
  char *newstr;
  int retval,term_height,term_width;
  newstr=(char *)malloc(10000*sizeof(char));
  fprintf(stderr,"\n\n");
  va_start(ap, formatstr);
  vsprintf(newstr,formatstr, ap);
  va_end(ap);
  find_terminal_window_size(&term_width,&term_height);
  fprintf(stderr,"%s",strwrp(newstr,min(term_width-1,70)));
  free(newstr);

  fprintf(stderr,"\n\nSOAP fatal error ");
  if (codex->action_being_processed!=NULL) fprintf(stderr,"within %s() ",codex->action_being_processed);
  fprintf(stderr,"in the vicinity of line %ld in file %s.\n\nExiting.\n\n",
          codex->linenum,codex->filename);
  exit(EXIT_FAILURE);
  retval=EXIT_FAILURE;
  return(retval);
}



/* cut _all characters from is to ie in *str */
void SOAP_strcut(long is, long ie, char *str){
  long i;
  i=is;
  do {
    str[i]=str[i+(ie-is)+1];
    i++;
  } while (str[i-1]!=EOS);
}

/* insert str1 into str2 at the position pos */
void SOAP_strins(char *str1, char **str2,  long pos){
  long len1,len2,i;
  len1=(long)strlen(str1);
  len2=(long)strlen(*str2);
  *str2=(char *)realloc(*str2,(len2+len1+3)*sizeof(char));

  for (i=len2; i>=pos; i--) (*str2)[i+len1]=(*str2)[i];
  for (i=0; i<len1; i++) (*str2)[pos+i]=str1[i];
  (*str2)[len2+len1]=EOS;
}

void SOAP_store_file_as_string(char *filename, char **str){
  FILE *file;
  long cnt;

  file = fopen(filename, "r");
  if (file==NULL) {
    fprintf(stderr,"\nHaving problems opening file %s.\nExiting.\n\n",filename);
    exit(EXIT_FAILURE);
  }
  cnt=0;
  do {
    *str=(char *)realloc(*str,(cnt+3)*sizeof(char));
    (*str)[cnt]=fgetc(file);
    cnt++;
  } while (!feof(file));
  fclose(file);
  (*str)[cnt-1]=EOS;
}


/* returns TRUE on success, FALSE otherwise
   find string str in expr following cnts
   anchorL and anchorR represent the boundaries of the string*/
static bool find_str_in_str(char *expr, char *str, long cnts, long *anchorL, long *anchorR)  {
  long len, i, c,cnt;

  len = (long)strlen(str);
  i = 0;
  cnt=cnts;

  while(i != len){
    c=expr[cnt];
    cnt++;
    *anchorR=cnt-1;
    *anchorL=cnt-len;
    if (c == EOS) return FALSE;
    if (str[i] == (char)c) {
      i++;
    } else {
      i = 0;
    }
  }
  return(TRUE);
}

static double random_double(double minval, double maxval){
  long randinput;
  double tmp;
  randinput=random();
  tmp=(double)(randinput)/(double)(RAND_MAX)*(maxval-minval)+minval;
  return(tmp);
}

/* returns true if expr[i] is a valid operator, false
   otherwise */
static bool is_operator(char *expr, long i){
  bool tmp;
  tmp=FALSE;
  if (expr[i]!=EOS) {
    if (   expr[i]=='*' || expr[i]=='/' || expr[i]=='^'
        || expr[i]==GT || expr[i]==GEQ || expr[i]==LT
        || expr[i]==LEQ || expr[i]==EQ || expr[i]==AND
        || expr[i]==OR || expr[i]==NEQ) tmp=TRUE;
    if (    (expr[i]=='+' || expr[i]=='-')
         && (i>0 && ((expr[i-1]>='0' && expr[i-1]<='9') || expr[i-1]=='.'))
       )
       tmp=TRUE;
  }
  return(tmp);
}

static bool is_logical(SOAP_codex_t *codex, long num){
  bool tmp;

  tmp=(bool)num;
  if (num!=0 && num!=1) {
    SOAP_fatal_error(codex,"Expecting a is_logical expression (0 or 1) but got %ld.",num);
  }
  return(tmp);
}


/* replaces the NOT operator '!' in *expr  */
static void replace_NOT_operator(SOAP_codex_t *codex, char **expr){
  long cnt,anchorL,anchorR;
  char *strR;
  double res;
  char *res_str;
  bool FOUND,res_bool;
  int eos=EOS;

  res_str=(char *)malloc(maxnumlen*sizeof(char));
  strR=(char *)malloc(sizeof(char));
  FOUND=FALSE;
  do {
    cnt=0;
    anchorL=0;
    FOUND=FALSE;
    do {
      if ( (*expr)[cnt]==NOT ) {
        anchorL=cnt;
        FOUND=TRUE;
      }
      cnt++;
    } while ((*expr)[cnt]!=EOS);
    if (FOUND) {
      anchorR=anchorL+1;
      do {
        anchorR++;
      } while ((!is_operator((*expr),anchorR)) && ((*expr)[anchorR]!=EOS));
      anchorR--;

      strR=(char *)realloc(strR,(anchorR-anchorL+3)*sizeof(char));
      for (cnt=anchorL+1; cnt<=anchorR; cnt++)
        strR[cnt-anchorL-1]=(*expr)[cnt];
      strR[anchorR-anchorL]=EOS;
      if (sscanf(strR,"%lg%n",&res,&eos)!=1 || strR[eos]!=EOS) {
        SOAP_fatal_error(codex,"Cannot read expression >%s<.",strR);
      }
      res_bool=is_logical(codex,longfromdouble(res));
      if (res_bool==0) strcpy(res_str,"1");
      if (res_bool==1) strcpy(res_str,"0");
      SOAP_strcut(anchorL,anchorR,*expr);
      SOAP_strins(res_str,expr,anchorL);
    }
  } while (FOUND);

  free(strR);
  free(res_str);
}


/* evaluate a string expression in which there are no
   parentheses "10*50/40^10.0E7" */
static double evaluate_arithmetic_1(SOAP_codex_t *codex, char *expr_orig){
  long cnt,priority,anchor,anchorL,anchorR;
  char *expr;
  char *strL,*strR;
  double tmp,res,numL,numR;
  char *res_str;
  bool SINGLENUM;
  int eos=EOS;

  res_str=(char *)malloc(maxnumlen*sizeof(char));
  expr=(char *)malloc(((long)strlen(expr_orig)+3)*sizeof(char));
  strL=(char *)malloc(sizeof(char));
  strR=(char *)malloc(sizeof(char));
  strcpy(expr,expr_orig);
  //if (expr[0]=='-' || expr[0]=='+') SOAP_strins("0.0",&expr,0); //???

//  do {
//    strrep(expr, "+-", "-");
//    strrep(expr, "--", "+");
//  } while(strstr(expr,"--")!=NULL || strstr(expr,"+-")!=NULL);
  if (expr[0]=='-' && expr[1]=='-') {
    SOAP_strins("0.0",&expr,0);
  }
  replace_NOT_operator(codex,&expr);
  SINGLENUM=FALSE;
  do {
    cnt=0;
    priority=-2;
    anchor=0;
    do {
      cnt++;
      if (is_operator(expr,cnt)) {
        if ((   expr[cnt]==AND || expr[cnt]==OR) && priority<-1) {
          priority=-1;
          anchor=cnt;
        }
        if ((   expr[cnt]==GT  || expr[cnt]==GEQ  || expr[cnt]==LT
              || expr[cnt]==LEQ || expr[cnt]==EQ || expr[cnt]==NEQ )
           && priority<0) {
          priority=0;
          anchor=cnt;
        }
        if ((expr[cnt]=='-' || expr[cnt]=='+') && priority<1) {
          priority=1;
          anchor=cnt;
        }
        if ((expr[cnt]=='*' || expr[cnt]=='/') && priority<2) {
          priority=2;
          anchor=cnt;
        }
        if ((expr[cnt]=='^') && priority<3) {
          priority=3;
          anchor=cnt;
        }
      }
    } while (expr[cnt]!=EOS);

    if (anchor!=0) {
      anchorL=anchor;
      do {
        anchorL--;
      } while ((!is_operator(expr,anchorL)) && (anchorL!=0));
      if (anchorL!=0) anchorL++; //???
      anchorR=anchor;
      do {
        anchorR++;
      } while ((!is_operator(expr,anchorR)) && (expr[anchorR]!=EOS));
      anchorR--;

      strL=(char *)realloc(strL,(anchor-anchorL+3)*sizeof(char));
      strR=(char *)realloc(strR,(anchorR-anchor+3)*sizeof(char));
      for (cnt=anchorL; cnt<anchor; cnt++)
        strL[cnt-anchorL]=expr[cnt];
      strL[anchor-anchorL]=EOS;
      for (cnt=anchor+1; cnt<=anchorR; cnt++)
        strR[cnt-anchor-1]=expr[cnt];
      strR[anchorR-anchor]=EOS;
      if (sscanf(strL,"%lg%n",&numL,&eos)!=1 || strL[eos]!=EOS) {
        SOAP_fatal_error(codex,"Problem reading expression >%s<.",strL);
      }
      if (sscanf(strR,"%lg%n",&numR,&eos)!=1 || strR[eos]!=EOS) {
        SOAP_fatal_error(codex,"Problem reading expression >%s<.",strR);
      }
      res=0.0e0; /* to avoid compiler warning */
      switch (expr[anchor]) {
        case OR:  if (is_logical(codex,longfromdouble(numL)) || is_logical(codex,longfromdouble(numR))) res=1.0e0; else res=0.0e0; break;
        case AND: if (is_logical(codex,longfromdouble(numL)) && is_logical(codex,longfromdouble(numR))) res=1.0e0; else res=0.0e0; break;
        case NEQ: if (numL!=numR) res=1.0e0; else res=0.0e0; break;
        case EQ:  if (numL==numR) res=1.0e0; else res=0.0e0; break;
        case GEQ: if (numL>=numR) res=1.0e0; else res=0.0e0; break;
        case LEQ: if (numL<=numR) res=1.0e0; else res=0.0e0; break;
        case LT:  if (numL<numR) res=1.0e0; else res=0.0e0; break;
        case GT:  if (numL>numR) res=1.0e0; else res=0.0e0; break;
        case '-': res=numL-numR; break;
        case '+': res=numL+numR; break;
        case '*': res=numL*numR; break;
        case '/': res=numL/numR; break;
        case '^': res=pow(numL,numR); break;
      }
      sprintf(res_str,DOUBLEFORMAT,res);
      SOAP_strcut(anchorL,anchorR,expr);
      SOAP_strins(res_str,&expr,anchorL);
    } else {
      SINGLENUM=TRUE;
    }

  } while (!SINGLENUM);
  if (sscanf(expr,"%lg%n",&tmp,&eos)!=1 || expr[eos]!=EOS) {
    SOAP_fatal_error(codex,"Problem reading expression >%s<.",expr);
  }
  free(strL);
  free(strR);
  free(expr);

  free(res_str);
  return(tmp);
}


/* evaluate a string expression in which there are
   parentheses "(10*50)/40^10.0E7" */
double SOAP_evaluate_arithmetic(SOAP_codex_t *codex, char *expr_orig){
  char *expr;
  char *expr2;
  char *res_str;
  long cnt,cnt2,anchorL,anchorR;
  double res;
  bool STILLSOME;

  res_str=(char *)malloc(maxnumlen*sizeof(char));
  expr=(char *)malloc(((long)strlen(expr_orig)+3)*sizeof(char));
  expr2=(char *)malloc(((long)strlen(expr_orig)+3)*sizeof(char));
  strcpy(expr,expr_orig);


  /* first check if parentheses are balanced */
  cnt2=0;
  for (cnt=0; cnt<(long)strlen(expr); cnt++){
    if (expr[cnt]=='(') cnt2++;
    if (expr[cnt]==')') cnt2--;
  }
  if (cnt2!=0) {
    SOAP_fatal_error(codex,"Parentheses not balanced: >%s<.",expr);
  }

  do {
    /* find expr2, the expression in brackets which needs to be
       evaluated first from expr and replace its value
       in expr*/
    STILLSOME=FALSE;
    cnt=0;
    anchorL=0;
    do {
      if (expr[cnt]=='(') {
        anchorL=cnt+1;
        STILLSOME=TRUE;
      }
      cnt++;
    } while (expr[cnt]!=')' && expr[cnt]!=EOS);
    anchorR=cnt-1;
    expr2=(char *)realloc(expr2,((long)strlen(expr)+3)*sizeof(char));
    for (cnt=anchorL; cnt<=anchorR; cnt++){
      expr2[cnt-anchorL]=expr[cnt];
    }
    expr2[anchorR-anchorL+1]=EOS;
    res=evaluate_arithmetic_1(codex,expr2);
    if (STILLSOME) {
      SOAP_strcut(anchorL-1,anchorR+1,expr);
      sprintf(res_str,DOUBLEFORMAT,res);
      SOAP_strins(res_str,&expr,anchorL-1);
    }
  } while (STILLSOME);
  free(expr);
  free(expr2);
  free(res_str);

  return(res);
}

/* is character a valid one for variables (or functions)? */
static bool is_part_of_var(char chr){
  bool ans;
  ans=FALSE;
  if (chr>='a' && chr<='z') ans=TRUE;
  if (chr>='A' && chr<='Z') ans=TRUE;
  if (chr>='0' && chr<='9') ans=TRUE;
  if (chr=='_' || chr=='.' || chr=='[' || chr==']') ans=TRUE;
  return(ans);
}

/* returns TRUE on success, FALSE otherwise
   find word str in expr following cnts
   anchorL and anchorR represent the boundaries of the word*/
static bool find_word_in_string(char *expr, char *word, long *anchorL, long *anchorR)  {
  bool WORD,STR;
  long cnts;
  cnts=0;
  do {
    STR=find_str_in_str(expr, word, cnts, anchorL, anchorR);
    WORD=TRUE;
    if (STR) {
      if (is_part_of_var(expr[*anchorR+1])) WORD=FALSE;
      if (*anchorL!=0) {
        if (is_part_of_var(expr[*anchorL-1])) WORD=FALSE;
      }
    } else {
      WORD=FALSE;
    }
    if (!WORD) cnts++;
  } while (!WORD && expr[cnts]!=EOS);
  return WORD;
}

/* substitute the expressions contained in _all array elements
   (between [ and ]) to their corresponding values */
static void substitute_array_elements(char **name, SOAP_codex_t *codex){
  long cnt,brackets;
  long anchorL,anchorR;
  char *expr;
  bool ENDREACHED;

  expr=(char *)malloc(sizeof(char));

  anchorL=0;
  ENDREACHED=FALSE;
  while (!ENDREACHED){
    brackets=0;
    do {
      anchorL++;
      if ((*name)[anchorL]=='[') brackets++;
      if ((*name)[anchorL]==EOS) ENDREACHED=TRUE;
    } while (brackets==0 && !ENDREACHED);
    anchorL++;
    if (!ENDREACHED) {
      cnt=anchorL;
      do {
        cnt++;
        if ((*name)[cnt]=='[') brackets++;
        if ((*name)[cnt]==']') brackets--;
      } while ( brackets!=0 && (*name)[cnt]!=EOS );
      if ((*name)[cnt]==EOS) {
          SOAP_fatal_error(codex,"Missing end of array character ] "
                                "in string >%s<.",*name);
      }
      anchorR=cnt-1;
      expr=(char *)realloc(expr,(anchorR-anchorL+3)*sizeof(char));
      for (cnt=anchorL; cnt<=anchorR; cnt++) expr[cnt-anchorL]=(*name)[cnt];
      expr[anchorR-anchorL+1]=EOS;

      SOAP_substitute_expression(&expr, codex);
      SOAP_strcut(anchorL,anchorR,*name);
      SOAP_strins(expr,name,anchorL);
    }
  }

  free(expr);
}

/* substitute variables in string *expr */
static void substitute_vars(char **expr, SOAP_codex_t *codex){
  long cnt,anchorL,anchorR;
  char *varname,*varvalue;

  /* need to set anchorL and anchorR to 0 to get rid of gcc warning */
  anchorL=0;
  anchorR=0;

  varname=(char *)malloc(maxnumlen*sizeof(char));
  varvalue=(char *)malloc(maxnumlen*sizeof(char));
  /* predefined variables first */
  for (cnt=0; cnt<7; cnt++){
    switch (cnt) {
      case 0: sprintf(varname,"pi");  sprintf(varvalue,DOUBLEFORMAT,pi); break;
      case 1: sprintf(varname,"TRUE");  sprintf(varvalue,"1"); break;
      case 2: sprintf(varname,"FALSE");  sprintf(varvalue,"0"); break;
      case 3: sprintf(varname,"YES");  sprintf(varvalue,"1"); break;
      case 4: sprintf(varname,"NO");  sprintf(varvalue,"0"); break;
      case 5: sprintf(varname,"EXIT_SUCCESS");  sprintf(varvalue,"0"); break;
      case 6: sprintf(varname,"EXIT_FAILURE");  sprintf(varvalue,"1"); break;
    }
    while (find_word_in_string(*expr, varname, &anchorL, &anchorR)){
      SOAP_strcut(anchorL,anchorR,*expr);
      SOAP_strins(varvalue,expr,anchorL);
    }
  }

  /* user defined variables second */
  substitute_array_elements(expr,codex);
  cnt=0;
  if (codex->vars[0].name!=NULL) {
    do {
      while (find_word_in_string(*expr, codex->vars[cnt].name, &anchorL, &anchorR)) {
        SOAP_strcut(anchorL,anchorR,*expr);
        SOAP_strins(codex->vars[cnt].value,expr,anchorL);
      }
      cnt++;
    } while (codex->vars[cnt].name!=NULL);
  }
  free(varname);
  free(varvalue);
}

/* get the nth argument and store it in *expr;
   arguments are counted from 0. *expr must already
   have been malloc'ed*/
static void get_argum_straight_0(SOAP_codex_t *codex, char **expr, char *argum, long n,
                              long *anchorL, long *anchorR){
  long cnt,parentheses;
  bool INSTRING;

  INSTRING=FALSE;
  *anchorL=0;
  for (cnt=0; cnt<n; cnt++){
    parentheses=0;
    do {
      if (argum[*anchorL]=='"') INSTRING=!INSTRING;
      if (argum[*anchorL]=='(' && !INSTRING) parentheses++;
      if (argum[*anchorL]==')' && !INSTRING) parentheses--;
      if (argum[*anchorL]==EOS) {
        SOAP_fatal_error(codex,"Reached end of string while trying to grab argument#%ld from string "
                              ">%s<.",n+1,argum);
      }
      (*anchorL)++;
    } while (!(argum[*anchorL]==',' && parentheses==0 && !INSTRING));
    (*anchorL)++;
  }

  cnt=*anchorL;
  parentheses=0;
  INSTRING=FALSE;
  do {
    if (argum[cnt]=='"') INSTRING=!INSTRING;
    if (argum[cnt]=='(' && !INSTRING) parentheses++;
    if (argum[cnt]==')' && !INSTRING) parentheses--;
    (*expr)=(char *)realloc(*expr,(cnt-(*anchorL)+3)*sizeof(char));
    (*expr)[cnt-(*anchorL)]=argum[cnt];
    cnt++;
  } while (!(argum[cnt]==',' && parentheses==0 && !INSTRING) && argum[cnt]!=EOS);

  (*anchorR)=cnt-1;
  (*expr)[(*anchorR)-(*anchorL)+1]=EOS;
}

/* get the nth argument and store it in *expr;
   arguments are counted from 0. *expr must already
   have been malloc'ed*/
void SOAP_get_argum_straight(SOAP_codex_t *codex, char **expr, char *argum, long n){
  long anchorR,anchorL;
  get_argum_straight_0(codex,expr,argum,n,&anchorL,&anchorR);
}

/* get the nth argument and store what is in between
   the quotes of the string in *expr;
   arguments are counted from 0. *expr must already
   have been malloc'ed*/
void SOAP_get_argum_string(SOAP_codex_t *codex, char **expr, char *argum, long n){
  long anchorR,anchorL;
  get_argum_straight_0(codex,expr,argum,n,&anchorL,&anchorR);
  if ((*expr)[0]=='"') {
    SOAP_strcut(0, 0, *expr);
  } else {
    SOAP_fatal_error(codex,"String does not start with \".");
  }
  if ((*expr)[strlen(*expr)-1]=='"') {
    SOAP_strcut(strlen(*expr)-1, strlen(*expr)-1, *expr);
  } else {
    SOAP_fatal_error(codex,"String does not end with \".");
  }
}


/* get the nth argument; arguments are counted from 0*/
double SOAP_get_argum_double(SOAP_codex_t *codex, char *argum, long n){
  char *expr;
  double tmp;
  int eos = EOS;

  expr=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&expr, argum, n);

  if (sscanf(expr,"%lg%n",&tmp,&eos)!=1 || expr[eos]!=EOS){
    SOAP_fatal_error(codex,"\"%s\" is not a float.",expr);
  }
  free(expr);
  return(tmp);
}


/* get the nth argument; arguments are counted from 0*/
long SOAP_get_argum_long(SOAP_codex_t *codex, char *argum, long n){
  char *expr;
  long tmp;
  int eos = EOS;

  expr=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&expr, argum, n);

  if (sscanf(expr,"%ld%n",&tmp,&eos)!=1 || expr[eos]!=EOS){
    SOAP_fatal_error(codex,"\"%s\" is not an integer.",expr);
  }
  free(expr);
  return(tmp);
}


/* get the nth argument; arguments are counted from 0*/
long SOAP_get_argum_bool(SOAP_codex_t *codex, char *argum, bool n){
  char *expr;
  long tmp;
  int eos = EOS;

  expr=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&expr, argum, n);

  if (sscanf(expr,"%ld%n",&tmp,&eos)!=1 || expr[eos]!=EOS){
    SOAP_fatal_error(codex,"\"%s\" is not a boolean.",expr);
  }
  if (tmp<0 || tmp>1) {
    SOAP_fatal_error(codex,"\"%s\" is not a boolean.",expr);
  }
  free(expr);
  return(tmp);
}



static void functions_builtin(char *function, char **argum,
             char **returnstr, SOAP_codex_t *codex){
  double tmp,returnval;
  long functionnum,numargum,cnt;
  int eos=EOS;
  functionnum=0;
  if (strcmp(function,"rad")==0) functionnum=1;
  if (strcmp(function,"deg")==0) functionnum=2;
  if (strcmp(function,"sin")==0) functionnum=3;
  if (strcmp(function,"cos")==0) functionnum=4;
  if (strcmp(function,"tan")==0) functionnum=5;
  if (strcmp(function,"asin")==0) functionnum=6;
  if (strcmp(function,"acos")==0) functionnum=7;
  if (strcmp(function,"atan")==0) functionnum=8;
  if (strcmp(function,"sqrt")==0) functionnum=9;
  if (strcmp(function,"sqr")==0) functionnum=10;
  if (strcmp(function,"exp")==0) functionnum=11;
  if (strcmp(function,"ln")==0) functionnum=12;
  if (strcmp(function,"round")==0) functionnum=13;
  if (strcmp(function,"floor")==0) functionnum=14;
  if (strcmp(function,"abs")==0) functionnum=15;
  if (strcmp(function,"sinh")==0) functionnum=16;
  if (strcmp(function,"cosh")==0) functionnum=17;
  if (strcmp(function,"tanh")==0) functionnum=18;
  if (strcmp(function,"asinh")==0) functionnum=19;
  if (strcmp(function,"acosh")==0) functionnum=20;
  if (strcmp(function,"atanh")==0) functionnum=21;

  if (functionnum>0 && functionnum<22) {
    SOAP_substitute_all_argums(argum, codex);
    if (sscanf(*argum,"%lg%n",&tmp,&eos)!=1 || (*argum)[eos]!=EOS)
      SOAP_fatal_error(codex,"Problem evaluating expression >%s<.",*argum);
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    returnval=0.0e0; /* to avoid compiler warning only */
    switch (functionnum) {
      case 1:  returnval=rad(tmp); break;
      case 2:  returnval=deg(tmp); break;
      case 3:  returnval=sin(tmp); break;
      case 4:  returnval=cos(tmp); break;
      case 5:  returnval=tan(tmp); break;
      case 6:  returnval=asin(tmp); break;
      case 7:  returnval=acos(tmp); break;
      case 8:  returnval=atan(tmp); break;
      case 9:  returnval=sqrt(tmp); break;
      case 10:  returnval=sqr(tmp); break;
      case 11:  returnval=exp(tmp); break;
      case 12:  returnval=log(tmp); break;
      case 13:  returnval=round(tmp); break;
      case 14:  returnval=floor(tmp); break;
      case 15:  returnval=fabs(tmp); break;
      case 16:  returnval=sinh(tmp); break;
      case 17:  returnval=cosh(tmp); break;
      case 18:  returnval=tanh(tmp); break;
      case 19:  returnval=asinh(tmp); break;
      case 20:  returnval=acosh(tmp); break;
      case 21:  returnval=atanh(tmp); break;
    }
    sprintf(*returnstr,DOUBLEFORMAT,returnval);
  }

  if (strcmp(function,"min")==0) {
    SOAP_substitute_all_argums(argum, codex);
    numargum=SOAP_number_argums(*argum);
    returnval=SOAP_get_argum_double(codex,*argum,0);
    for (cnt=1; cnt<numargum; cnt++)
       returnval=min(returnval,SOAP_get_argum_double(codex,*argum,cnt));
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    sprintf(*returnstr,DOUBLEFORMAT,returnval);
  }

  if (strcmp(function,"max")==0) {
    SOAP_substitute_all_argums(argum, codex);
    numargum=SOAP_number_argums(*argum);
    returnval=SOAP_get_argum_double(codex,*argum,0);
    for (cnt=1; cnt<numargum; cnt++)
       returnval=max(returnval,SOAP_get_argum_double(codex,*argum,cnt));
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    sprintf(*returnstr,DOUBLEFORMAT,returnval);
  }

  if (strcmp(function,"random")==0) {
    SOAP_substitute_all_argums(argum, codex);
    returnval=random_double(SOAP_get_argum_double(codex,*argum,0),SOAP_get_argum_double(codex,*argum,1));
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    sprintf(*returnstr,DOUBLEFORMAT,returnval);
  }

  if (strcmp(function,"mod")==0) {
    SOAP_substitute_all_argums(argum, codex);
    returnval=mod(longfromdouble(SOAP_get_argum_double(codex,*argum,0)),longfromdouble(SOAP_get_argum_double(codex,*argum,1)));
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    sprintf(*returnstr,DOUBLEFORMAT,returnval);
  }

  if (strcmp(function,"krodelta")==0) {
    SOAP_substitute_all_argums(argum, codex);
    returnval=krodelta(longfromdouble(SOAP_get_argum_double(codex,*argum,0)),longfromdouble(SOAP_get_argum_double(codex,*argum,1)));
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    sprintf(*returnstr,DOUBLEFORMAT,returnval);
  }


  if (strcmp(function,"defined")==0) {
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    if (sscanf(*argum,"%lg%n",&tmp,&eos)==1 && (*argum)[eos]==EOS) strcpy(*returnstr,"1");
      else strcpy(*returnstr,"0");
  }

  if (strcmp(function,"spline")==0) {
    long N,n;
    double *f,*b,*x;
    double thisx;

    N=SOAP_number_argums(*argum);
    if (mod(N-1,2)!=0) SOAP_fatal_error(codex,"Number of arguments within spline must be an odd number.");
    N=(N-1)/2;
    if (N<4) SOAP_fatal_error(codex,"Number of data points supplied within spline must be at least 4.");
  
    x=(double *)malloc(N*sizeof(double));
    f=(double *)malloc(N*sizeof(double));
    b=(double *)malloc(N*sizeof(double));
  
    for (n=0; n<N; n++) {
      SOAP_substitute_argum(argum, n*2, codex);    
      x[n]=SOAP_get_argum_double(codex, *argum, n*2);
      SOAP_substitute_argum(argum, n*2+1, codex);    
      f[n]=SOAP_get_argum_double(codex, *argum, n*2+1);
    }
    /* check if data points are valid (x[n+1]>x[n]) */
    for (n=0; n<N-1; n++){
      if (x[n+1]<=x[n]) SOAP_fatal_error(codex, "Data points supplied to spline must be such that x[i+1]>x[i]."); 
    }

    SOAP_substitute_argum(argum, N*2, codex);    
    thisx=SOAP_get_argum_double(codex, *argum, N*2);
    /* check if point is out of range */
    if (thisx<x[0] || thisx>x[N-1]) SOAP_fatal_error(codex, "Ensure that x lies between x[0] and x[N].");
  
    EXM_find_spline(N, x, f, b);
    *returnstr=(char *)realloc(*returnstr,maxnumlen*sizeof(char));
    sprintf(*returnstr,DOUBLEFORMAT,EXM_f_from_spline(N, x, f, b, thisx));

    free(x);
    free(b);
    free(f);
  }

}

static void substitute_functions(char **expr, SOAP_codex_t *codex){
  long cnt,anchorLL,anchorL,anchorR,parentheses;
  bool FOUND;
  char *function,*argum,*returnstr;

  returnstr=(char *)malloc(sizeof(char));
  do {
    cnt=0;
    anchorL=0;
    FOUND=FALSE;
    do {
      cnt++;
      if ((*expr)[cnt]=='(' && is_part_of_var((*expr)[cnt-1])) {
        FOUND=TRUE;
        anchorL=cnt;
      }
    } while ((*expr)[cnt]!=EOS);


    if (FOUND) {
      cnt=anchorL;
      do {
        cnt--;
      } while (cnt>0 && is_part_of_var((*expr)[cnt-1]) );
      anchorLL=cnt;
      cnt=anchorL;
      parentheses=0;
      do {
        if ((*expr)[cnt]=='"') SOAP_fatal_error(codex,"Found a quote inside an expression.");
        if ((*expr)[cnt]=='(') parentheses++;
        if ((*expr)[cnt]==')') parentheses--;
        cnt++;
      } while (parentheses!=0);
      anchorR=cnt-1;
      function=(char *)malloc((anchorL-anchorLL+4)*sizeof(char));
      argum=(char *)malloc((anchorR-anchorL+4)*sizeof(char));
      for (cnt=anchorLL; cnt<anchorL; cnt++)
        function[cnt-anchorLL]=(*expr)[cnt];
      function[anchorL-anchorLL]=EOS;
      for (cnt=anchorL+1; cnt<anchorR; cnt++)
        argum[cnt-anchorL-1]=(*expr)[cnt];
      argum[anchorR-anchorL-1]=EOS;
      returnstr=(char *)realloc(returnstr,(3+(long)strlen(function))*sizeof(char));
      strcpy(returnstr,function);
      if (codex->FUNCTION) (codex->function)(function,&argum,&returnstr,codex);
      functions_builtin(function, &argum, &returnstr, codex);
      SOAP_strcut(anchorLL,anchorR,*expr);
      SOAP_strins(returnstr,expr,anchorLL);
      free(argum);
      free(function);
    }
  } while (FOUND);
  free(returnstr);
}


bool SOAP_is_double_a_long(double expr_double){
  double expr_double2;
  long expr_long;
  bool RET;
  expr_long=(long)expr_double;
  expr_double2=(double)expr_long;
  if (expr_double2!=expr_double) RET=FALSE; else RET=TRUE;
  return(RET);
}

/* evaluate the expr **expr and substitute it for
   the number it stands for. This includes evaluating
   the variables, the functions and finally the arithmetic */
void SOAP_substitute_expression(char **expr, SOAP_codex_t *codex){
  double expr_double;
  char *expr_str;
  expr_str=(char *)malloc(maxnumlen*sizeof(char));
  substitute_vars(expr, codex);
  substitute_functions(expr, codex);
  expr_double=SOAP_evaluate_arithmetic(codex,*expr);
  
  (*expr)[0]=EOS;
  if (!SOAP_is_double_a_long(expr_double)){
    sprintf(expr_str,DOUBLEFORMAT,expr_double);
    SOAP_strins(expr_str,expr,0);
  } else {
    sprintf(expr_str,"%ld",(long)expr_double);
    SOAP_strins(expr_str,expr,0);
  }
  free(expr_str);
}


/* in the given *expr, substitute the string comparisons,
   like "bernard"=="jeny" would be evaluated to false,
   and substituted by 0 */
void SOAP_substitute_string_arithmetic(char **expr, SOAP_codex_t *codex){
  char *leftstring;
  char *rightstring;
  long cnt,anchor,anchorL,anchorR;
  bool INSTRING;
  int strcmp_ret;
  char newexpr[2];

  if (strlen(*expr)>4) {
    anchor=0;
    if ((*expr)[0]=='"') INSTRING=TRUE; else INSTRING=FALSE;
    do {
      anchor++;
      if ((*expr)[anchor]=='"') INSTRING=!INSTRING;
      if (   ((*expr)[anchor]==EQ || (*expr)[anchor]==NEQ
          ||  (*expr)[anchor]==LT || (*expr)[anchor]==GT
          ||  (*expr)[anchor]==GEQ || (*expr)[anchor]==LEQ)
          && (*expr)[anchor-1]=='"' && (*expr)[anchor+1]=='"' && !INSTRING) {
        /* first find left string */
        anchorL=anchor-1;
        do {
          anchorL--;
        } while ((*expr)[anchorL]!='"' && anchorL>0);
        if ((*expr)[anchorL]!='"') SOAP_fatal_error(codex,"Problem in subroutine SOAP_StringArith(2).");
        leftstring=(char *)malloc(sizeof(char)*(anchor-anchorL+2));
        for (cnt=anchorL+1; cnt<anchor-1; cnt++){
          leftstring[cnt-anchorL-1]=(*expr)[cnt];
        }
        leftstring[anchor-anchorL-2]=EOS;
        /* find right string */
        anchorR=anchor+1;
        do {
          (anchorR)++;
        } while ((*expr)[anchorR]!='"' && (*expr)[anchorR]!=EOS);
        if ((*expr)[anchorR]!='"') SOAP_fatal_error(codex,"Problem in subroutine SOAP_StringArith(3).");
        rightstring=(char *)malloc(sizeof(char)*(anchorR-anchor+2));
        for (cnt=anchor+2; cnt<anchorR; cnt++){
          rightstring[cnt-anchor-2]=(*expr)[cnt];
        }
        rightstring[anchorR-anchor-2]=EOS;
        /* printf("leftstring=>%s< rightstring=>%s<\n",leftstring,rightstring); */
        /* compare the two strings and find newexpr */
        strcmp_ret=strcmp(leftstring,rightstring);
        switch ((*expr)[anchor]) {
          case EQ:  if (strcmp_ret==0) (newexpr)[0]='1'; else (newexpr)[0]='0';
	      break;
          case NEQ: if (strcmp_ret!=0) (newexpr)[0]='1'; else (newexpr)[0]='0';
	      break;
          case LT: if (strcmp_ret<0) (newexpr)[0]='1'; else (newexpr)[0]='0';
	      break;
          case GT: if (strcmp_ret>0) (newexpr)[0]='1'; else (newexpr)[0]='0';
	      break;
          case LEQ: if (strcmp_ret<=0) (newexpr)[0]='1'; else (newexpr)[0]='0';
	      break;
          case GEQ: if (strcmp_ret>=0) (newexpr)[0]='1'; else (newexpr)[0]='0';
	      break;
        }
        (newexpr)[1]=EOS;
        SOAP_strcut(anchorL,anchorR,*expr);
        SOAP_strins(newexpr, expr, anchorL);
        anchor=anchorL;
        free(leftstring);
        free(rightstring);
      }
    } while ((*expr)[anchor+1]!=EOS);
  }
}

/* evaluate the expr **expr and substitute it for
   the number it stands for. This includes evaluating
   the variables, the functions and finally the arithmetic of all
   sub-expressions outside of the strings. Sub-expressions which are
   strings are left untouched. If some sub-expressions are strings,
   then the final expression is made into a string*/
void SOAP_substitute_expression_including_strings(char **expr, SOAP_codex_t *codex){
  long anchorL, anchorR, cnt, pass;
  char *exprtmp;
  bool QUOTEFOUND,INSTRING;
  for (pass=1; pass<=2; pass++){
    if (pass==2) SOAP_substitute_string_arithmetic(expr,codex);
    anchorL=0;
    do {
      /* find the anchorL and anchorR corresponding to non-string */
      anchorL--;
      INSTRING=FALSE;
      do {
        anchorL++;
        if ((*expr)[anchorL]=='"') INSTRING=!INSTRING;
      } while (INSTRING || (*expr)[anchorL]=='"');
      anchorR=anchorL;
      if ((*expr)[anchorR]!=EOS){
        do {
          anchorR++;
        } while ((*expr)[anchorR]!='"' && (*expr)[anchorR]!=EOS);
        anchorR--;
      }
      /* substitute variables or expressions*/
      if ((*expr)[anchorL]!=EOS) {
        exprtmp=(char *)malloc(sizeof(char)*(anchorR-anchorL+3));
        for (cnt=anchorL; cnt<=anchorR; cnt++) exprtmp[cnt-anchorL]=(*expr)[cnt];
        exprtmp[anchorR-anchorL+1]=EOS;
        if (pass==1) substitute_vars(&exprtmp, codex);
        if (pass==2) SOAP_substitute_expression(&exprtmp, codex);
        SOAP_strcut(anchorL,anchorR,*expr);
        SOAP_strins(exprtmp, expr, anchorL);
        anchorL+=strlen(exprtmp);
        free(exprtmp);
      }
    } while ((*expr)[anchorL]!=EOS);
  }
  /* clean up _all quotes and, if quotes were found, add quotes at end and start of *expr */
  anchorL=0;
  QUOTEFOUND=FALSE;
  do {
    if ((*expr)[anchorL]=='"') {
      SOAP_strcut(anchorL,anchorL,*expr);
      anchorL--;
      QUOTEFOUND=TRUE;
    }
    anchorL++;
  } while((*expr)[anchorL]!=EOS);
  if (QUOTEFOUND) {
    SOAP_strins("\"", expr, 0);
    SOAP_strins("\"", expr, strlen(*expr));
  }

}


/* sub the expression located in nth argument
   with SOAP_substitute_expression; arguments are counted from 0*/
void SOAP_substitute_argum(char **argum, long n, SOAP_codex_t *codex){
  char *expr;
  long anchorL,anchorR;

  expr=(char *)malloc((long)(strlen(*argum)+3)*sizeof(char));
  get_argum_straight_0(codex,&expr, *argum, n, &anchorL, &anchorR);

  SOAP_substitute_expression_including_strings(&expr, codex);
  SOAP_strcut(anchorL,anchorR,*argum);
  SOAP_strins(expr,argum,anchorL);
  free(expr);
}

/* returns the number of arguments */
long SOAP_number_argums(char *argum){
  long cnt,commas,parentheses;
  bool INSTRING;
  if (strlen(argum)==0) {
    commas=0;
  } else {
    cnt=0;
    commas=0;
    parentheses=0;
    INSTRING=FALSE;
    do {
      if (argum[cnt]=='"') INSTRING=!INSTRING;
      if (argum[cnt]==',' && parentheses==0 && !INSTRING) commas++;
      if (argum[cnt]=='(' && !INSTRING) parentheses++;
      if (argum[cnt]==')' && !INSTRING) parentheses--;
      cnt++;
    } while(argum[cnt]!=EOS);
    commas++;
  }
  return(commas);
}

void SOAP_substitute_all_argums(char **argum, SOAP_codex_t *codex){
  long cnt,numargum;
  numargum=SOAP_number_argums(*argum);
  for (cnt=0; cnt<numargum; cnt++) SOAP_substitute_argum(argum, cnt, codex);
}

/*
static void ShowCode(char *code){
  long cnt;
  printf("code starts here -->");
  for (cnt=0; cnt<(long)strlen(code); cnt++){
    printf("%c",code[cnt]);
  }
  printf("<-- code ended there\n");
  fflush(stdout);
}*/

static void delete_character(long cnt, char *code){
  long cnt2;
  long codelength;

  codelength=(long)strlen(code);
  for (cnt2=cnt+1; cnt2<codelength; cnt2++)
    code[cnt2-1]=code[cnt2];
  code[codelength-1]=EOS;
}

static void clean_comments(SOAP_codex_t *codex, char *code){
  long cnt;
  long parentheses,anchorL,anchorR;
  bool INSTRING;
  char *lastnewline;
  cnt=0;
  parentheses=0;
  anchorL=0;
  INSTRING=FALSE;
  while (code[cnt]!=EOS) {
    if (code[cnt]=='"') INSTRING=!INSTRING;
    if (code[cnt]=='{' && !INSTRING) {
      parentheses++;
      if (parentheses==1) anchorL=cnt;
    }
    if (code[cnt]=='}' && !INSTRING) {
      parentheses--;
    
      if (parentheses==0) {
        anchorR=cnt;
        SOAP_strcut(anchorL,anchorR,code);
        cnt=cnt-(anchorR-anchorL+1);
      }
      if (parentheses==-1){
        code[cnt]=EOS;
        lastnewline=code;
        while( (code = strstr(code,"__newline"))){
          lastnewline = code++;
        }
        if (sscanf(lastnewline,"__newline(%ld)",&(codex->linenum))!=1) codex->linenum=-1;
        SOAP_fatal_error(codex,"Comment closed but not opened.");
      }
    }
    cnt++;
  }
  if (code[cnt]==EOS && INSTRING) {
    lastnewline=code;
    while( (code = strstr(code,"__newline"))){
      lastnewline = code++;
    }
    if (sscanf(lastnewline,"__newline(%ld)",&(codex->linenum))!=1) codex->linenum=-1;
    SOAP_fatal_error(codex,"String not closed properly.");
  }
  if (code[cnt]==EOS && parentheses!=0) {
    lastnewline=code;
    while( (code = strstr(code,"__newline"))){
      lastnewline = code++;
    }
    if (sscanf(lastnewline,"__newline(%ld)",&(codex->linenum))!=1) codex->linenum=-1;
    SOAP_fatal_error(codex,"Comment not closed properly.");
  }
}



/* only insert the __newline(); action in front of an action*/
void SOAP_insert_line_numbers_in_code_backward(char **code, long linenum_start){
  long cnt,linenum,linenum2,cnt2,parentheses,commentbrackets;
  char *newlinestr;
  bool INSTRING,INSTRING2,CONTINUE;
  newlinestr=(char *)malloc(sizeof(char)*100);
  cnt=0;
  linenum=linenum_start;
  sprintf(newlinestr,"__newline(%ld);",linenum);
  SOAP_strins(newlinestr, code,  cnt);
  cnt=cnt+strlen(newlinestr);
  INSTRING=FALSE;
  commentbrackets=0;
  CONTINUE=TRUE;
  do {
    if ((*code)[cnt]=='"') INSTRING=!INSTRING;
    if ((*code)[cnt]=='{') commentbrackets++;
    if ((*code)[cnt]=='}') commentbrackets--;
    if ((*code)[cnt]=='\n') linenum++;
    if ((*code)[cnt]==';' && !INSTRING && commentbrackets==0) {
    /* here, go back to previous ; or , or ( making sure parentheses is zero
       and not in a string */
      cnt2=cnt;
      linenum2=linenum;
      parentheses=0;
      INSTRING2=FALSE;
      do {
        cnt2--;
        if ((*code)[cnt2]=='\n') linenum2--;
        if ((*code)[cnt2+1]=='(' && !INSTRING2 && commentbrackets==0) parentheses--;
        if ((*code)[cnt2+1]==')' && !INSTRING2 && commentbrackets==0) parentheses++;
        if ((*code)[cnt2+1]=='"') INSTRING2=!INSTRING2;
        if ((*code)[cnt2+1]=='{') commentbrackets--;
        if ((*code)[cnt2+1]=='}') commentbrackets++;
        if (cnt2<0) {
          CONTINUE=FALSE;
          /* fprintf(stderr,"\n\nproblem inserting line numbers aroung line %ld\n\n",linenum);
          exit(EXIT_FAILURE); */
        } 
      } while (CONTINUE && (((*code)[cnt2]!='(' && (*code)[cnt2]!=';' && (*code)[cnt2]!=',' && cnt2!=0)
               || INSTRING2 || parentheses!=0 || commentbrackets!=0));
    /* then, do this */
      if (CONTINUE){
       do {
        cnt2++;
        if ((*code)[cnt2]=='\n') linenum2++;
        if ((*code)[cnt2]=='{') commentbrackets++;
        if ((*code)[cnt2-1]=='}') commentbrackets--;
        if ((*code)[cnt2]==EOS) {
          CONTINUE=FALSE;
          /* fprintf(stderr,"\n\nproblem inserting line numbers aroung line %ld\n\n",linenum);
           exit(EXIT_FAILURE); */
        } 
       } while(CONTINUE && ((*code)[cnt2]==' '  || (*code)[cnt2]=='\n'
           || (*code)[cnt2]=='\t' || (*code)[cnt2]==13 || commentbrackets!=0));
       if (CONTINUE){
        sprintf(newlinestr,"__newline(%ld);",linenum2);
        SOAP_strins(newlinestr, code,  cnt2);
        cnt=cnt+strlen(newlinestr);
        commentbrackets=0;
       }
      }
    }
    cnt++;
  } while((*code)[cnt]!=EOS && CONTINUE);
}


/* insert the __newline(); action at all newlines */
void SOAP_insert_line_numbers_in_code(char **code, long linenum_start){
  long cnt,linenum,linenumwritten,commentbrackets;
  char *newlinestr;
  bool INSTRING;
  newlinestr=(char *)malloc(sizeof(char)*1000);
  linenum=linenum_start;
  sprintf(newlinestr,"__newline(%ld);",linenum);
  SOAP_strins(newlinestr, code,  0);


  //SOAP_insert_line_numbers_in_code_backward(code, linenum_start);

  /* add a __newline() command after each ';' followed by spaces, tabs or new lines */
  cnt=0;
  INSTRING=FALSE;
  commentbrackets=0;
  linenumwritten=-1;
  do {
    if ((*code)[cnt]=='"') INSTRING=!INSTRING;
    if ((*code)[cnt]=='{') commentbrackets++;
    if ((*code)[cnt]=='}') commentbrackets--;
    if ((*code)[cnt]=='\n') linenum++;
    if ((*code)[cnt]==';' && !INSTRING && commentbrackets==0) {
      while ((*code)[cnt+1]=='\n' || (*code)[cnt+1]==' ' || (*code)[cnt+1]=='\t' || (*code)[cnt+1]=='{'){         

        if ((*code)[cnt+1]=='{') {
          commentbrackets++;
          while (commentbrackets>0 && (*code)[cnt+2]!=EOS) {
            cnt++;
            if ((*code)[cnt+1]=='{') commentbrackets++;
            if ((*code)[cnt+1]=='\n') linenum++;
            if ((*code)[cnt+1]=='}') commentbrackets--;
/*            if ((*code)[cnt+1]==EOS) {
              fprintf(stderr,"\n\nComment not closed properly. SOAP fatal error in the vicinity of line %ld.\n\nExiting.\n\n",linenum);
              exit(EXIT_FAILURE);
            }*/
          }
          // at this point, (*code)[cnt+1]='}'
        }


        if ((*code)[cnt+1]=='\n') {
          linenum++;
          sprintf(newlinestr,"__newline(%ld);",linenum);
          linenumwritten=linenum;
          if ((*code)[cnt+2]!=EOS) {
            SOAP_strins(newlinestr, code,  cnt+2);
            cnt+=strlen(newlinestr);
          }
        }
        cnt++;
      }
      if (linenum!=linenumwritten) {
        sprintf(newlinestr,"__newline(%ld);",linenum);
        if ((*code)[cnt+1]!=EOS) SOAP_strins(newlinestr, code,  cnt+1);
        cnt=cnt+strlen(newlinestr);
      }
    }
    cnt++;
  } while((*code)[cnt]!=EOS);
  free(newlinestr);
  //printf("%s",*code);
}





static void clean_code(SOAP_codex_t *codex, char *code){
  long cnt;
  bool WINDOWS,INSTRING;


  clean_comments(codex,code);
  cnt=0;
  WINDOWS=FALSE;
  INSTRING=FALSE;
  while (code[cnt]!=EOS) {
    if (code[cnt]=='"') INSTRING=!INSTRING;
    if (!INSTRING) {
      if (code[cnt]=='=' && code[cnt+1]=='=') {
        code[cnt]=EQ;
        code[cnt+1]=' ';
      }
      if (code[cnt]=='!' && code[cnt+1]=='=') {
        code[cnt]=NEQ;
        code[cnt+1]=' ';
      }
      if (code[cnt]=='>' && code[cnt+1]=='=') {
        code[cnt]=GEQ;
        code[cnt+1]=' ';
      }
      if (code[cnt]=='<' && code[cnt+1]=='=') {
        code[cnt]=LEQ;
        code[cnt+1]=' ';
      }
      if (code[cnt]=='&' && code[cnt+1]=='&') {
        code[cnt]=AND;
        code[cnt+1]=' ';
      }
      if (code[cnt]=='|' && code[cnt+1]=='|') {
        code[cnt]=OR;
        code[cnt+1]=' ';
      }
      if (code[cnt]=='<') code[cnt]=LT;
      if (code[cnt]=='>') code[cnt]=GT;
      if (code[cnt]=='!') code[cnt]=NOT;
    }
    if (code[cnt]==13) WINDOWS=TRUE;
    if ((   code[cnt]==' '
         || code[cnt]=='\n'
         || code[cnt]=='\t'
         || code[cnt]==13
        ) && (!INSTRING)
        ) delete_character(cnt,code);
    else cnt++;
  }
  if (WINDOWS)
    fprintf(stdout,"Your code contains some DOS(TM) end-of-line characters (#13). \n");
  if (INSTRING) fprintf(stdout,"String not closed properly.\n");

}



/* update the variable named *name to the value
   specified in *argum */
static void update_var(char **name, char **argum, SOAP_codex_t *codex){
  long cnt;
  bool FOUND;

  substitute_array_elements(name,codex);
  SOAP_substitute_argum(argum,0,codex);
  cnt=0;
  FOUND=FALSE;
  if ((codex->vars)[0].name!=NULL) {
    do {
     if (strcmp(*name,(codex->vars)[cnt].name)==0) {
        FOUND=TRUE;
        (codex->vars)[cnt].value=(char *)realloc((codex->vars)[cnt].value,
                           ((long)strlen(*argum)+3)*sizeof(char));
        strcpy((codex->vars)[cnt].value,*argum);
      }
      cnt++;
    } while ((codex->vars)[cnt].name!=NULL);
  }
  if (!FOUND) {
    codex->vars=(SOAP_vars_t *)realloc(codex->vars,(cnt+5)*sizeof(SOAP_vars_t));
    (codex->vars)[cnt].name=(char *)malloc(((long)strlen(*name)+3)*sizeof(char));
    (codex->vars)[cnt].value=(char *)malloc(((long)strlen(*argum)+3)*sizeof(char));
    strcpy((codex->vars)[cnt].name,*name);
    strcpy((codex->vars)[cnt].value,*argum);
    (codex->vars)[cnt+1].name=NULL;
  }
}

/* the builtin actions */

static void BA_write(char **argum, SOAP_codex_t *codex, bool NEWLINE){
  if (codex->SCREENOUTPUT) {
    SOAP_substitute_all_argums(argum,codex);
    fprintf(stdout,"%s",*argum);
    if (NEWLINE) fprintf(stdout,"\n");
    fflush(stdout);
  }
}


static void BA_printf(char **argum, SOAP_codex_t *codex){
  long cnt,numargum,cnt2;
  char **argv;
  numargum=SOAP_number_argums(*argum);
  if (numargum<1) SOAP_fatal_error(codex,"Number of arguments given to printf must be at least 1.");
  assert(numargum>0);
  for (cnt=0; cnt<numargum; cnt++) SOAP_substitute_argum(argum,cnt,codex);
  argv=(char **)malloc(numargum*sizeof(char *));
  for (cnt=0; cnt<numargum; cnt++){
    argv[cnt]=(char *)malloc(sizeof(char));
    SOAP_get_argum_straight(codex,&(argv[cnt]), *argum, cnt);
  }
  for (cnt=0; cnt<numargum; cnt++){
    cnt2=0;
    do {
      if (argv[cnt][cnt2]=='"') {
        SOAP_strcut(cnt2,cnt2,argv[cnt]);
        cnt2--;
      }
      cnt2++;
    } while (argv[cnt][cnt2]!=EOS);
  }
  /* send everything to printf, to printf to stdout */
  if (codex->SCREENOUTPUT) SOAP_printf((int)numargum,argv,stdout);
  /* free pointers */
  for (cnt=0; cnt<numargum; cnt++) free(argv[cnt]);
  free(argv);
}




static void BA_fprintf(char **argum, SOAP_codex_t *codex){
  long numargum,cnt,cnt2;
  FILE *stream;
  char **argv;
  char *filename;

  numargum=SOAP_number_argums(*argum);
  if (numargum<2) SOAP_fatal_error(codex,"Number of arguments given to fprintf must be at least 2: the filename and the string to print.");
  assert(numargum>1);
  for (cnt=0; cnt<numargum; cnt++) SOAP_substitute_argum(argum,cnt,codex);
  argv=(char **)malloc(numargum*sizeof(char *));
  for (cnt=1; cnt<numargum; cnt++){
    argv[cnt-1]=(char *)malloc(sizeof(char));
    SOAP_get_argum_straight(codex,&(argv[cnt-1]), *argum, cnt);
  }
  for (cnt=1; cnt<numargum; cnt++){
    cnt2=0;
    do {
      if (argv[cnt-1][cnt2]=='"') {
        SOAP_strcut(cnt2,cnt2,argv[cnt-1]);
        cnt2--;
      }
      cnt2++;
    } while (argv[cnt-1][cnt2]!=EOS);
  }
  filename=(char *)malloc(sizeof(char));
  SOAP_get_argum_string(codex, &filename, *argum, 0);
  if (codex->FILEOUTPUT) {
    /* here, append to the file named in the first argument */
    stream=fopen(filename,"a");
    /* send everything to printf*/
    SOAP_printf((int)numargum-1,argv,stream);
    fclose(stream);
  }
  /* free pointers */
  for (cnt=1; cnt<numargum; cnt++) free(argv[cnt-1]);
  free(filename);
  free(argv);
}


static void BA_for(char **argum, SOAP_codex_t *codex){
  char *cntstr,*loopcode,*cntstr2;
  long cnts,cnte,cnt;

  if (SOAP_number_argums(*argum)!=4)
  SOAP_fatal_error(codex,"the for() command needs 4 arguments: "
                          "the first argument is the counter variable name; "
			  "the second argument is the start of the counting (integer); "
			  "the third argument is the end of the counting (integer); "
			  "the fourth argument is the code to be executed at every count.");
  cntstr=(char *)malloc(sizeof(char));
  cntstr2=(char *)malloc(sizeof(char));
  loopcode=(char *)malloc(sizeof(char));
  SOAP_substitute_argum(argum,1,codex);
  SOAP_substitute_argum(argum,2,codex);
  cnts=SOAP_get_argum_long(codex,*argum,1);
  cnte=SOAP_get_argum_long(codex,*argum,2);
  SOAP_get_argum_straight(codex,&cntstr,*argum,0);
  SOAP_get_argum_straight(codex,&loopcode, *argum, 3);
  if (cnts<=cnte) {
    for (cnt=cnts; cnt<=cnte; cnt++){
      /* change value of cntstr to cntstr2 in variables */
      cntstr2=(char *)realloc(cntstr2,maxnumlen*sizeof(char));
      sprintf(cntstr2,"%ld",cnt);
      update_var(&cntstr, &cntstr2, codex);
      SOAP_process_code(loopcode, codex, SOAP_VARS_KEEP_ALL);
    }
  } else {
    for (cnt=cnts; cnt>=cnte; cnt--){
      /* change value of cntstr to cntstr2 in variables */
      cntstr2=(char *)realloc(cntstr2,maxnumlen*sizeof(char));
      sprintf(cntstr2,"%ld",cnt);
      update_var(&cntstr, &cntstr2, codex);
      SOAP_process_code(loopcode, codex, SOAP_VARS_KEEP_ALL);
    }
  }
  free(cntstr);
  free(cntstr2);
  free(loopcode);

}


static void BA_for_parallel(char **argum, SOAP_codex_t *codex){
  char *cntstr,*loopcode,*cntstr2;
  long cnts,cnte,cnt,cntvar,numvars,cnttmp;
  SOAP_codex_t *codexcopy;

  if (SOAP_number_argums(*argum)!=4)
  SOAP_fatal_error(codex,"the for_parallel() command needs 4 arguments: "
                          "the first argument is the counter variable name; "
			  "the second argument is the start of the counting (integer); "
			  "the third argument is the end of the counting (integer); "
			  "the fourth argument is the code to be executed at every count.");
  SOAP_substitute_argum(argum,1,codex);
  SOAP_substitute_argum(argum,2,codex);

  cnts=SOAP_get_argum_long(codex,*argum,1);
  cnte=SOAP_get_argum_long(codex,*argum,2);
  if (cnts>cnte) {
    cnttmp=cnte;
    cnte=cnts;
    cnts=cnttmp;
  }
  codexcopy=(SOAP_codex_t *)malloc((cnte-cnts+2)*sizeof(SOAP_codex_t));
#ifdef OPENMPTHREADS
#pragma omp parallel for private(cnt,cntstr2,cntstr,loopcode) schedule(dynamic) 
#endif
  for (cnt=cnts; cnt<=cnte; cnt++){
    SOAP_copy_codex(codex, &(codexcopy[cnt-cnts]));
    (codexcopy[cnt-cnts]).vars=NULL;
    SOAP_copy_all_vars(codex->vars, &((codexcopy[cnt-cnts]).vars));
    /* change value of cntstr to cntstr2 in variables */
    loopcode=(char *)malloc(sizeof(char));
    cntstr=(char *)malloc(sizeof(char));
    SOAP_get_argum_straight(&(codexcopy[cnt-cnts]),&cntstr,*argum,0);
    SOAP_get_argum_straight(&(codexcopy[cnt-cnts]),&loopcode, *argum, 3);
    cntstr2=(char *)malloc(maxnumlen*sizeof(char));
    sprintf(cntstr2,"%ld",cnt);
    update_var(&cntstr, &cntstr2, &(codexcopy[cnt-cnts]));
    SOAP_process_code(loopcode, &(codexcopy[cnt-cnts]), SOAP_VARS_KEEP_ALL);
    free(cntstr2);
    free(cntstr);
    free(loopcode);
  }
  for (cnt=cnts; cnt<=cnte; cnt++){
    SOAP_count_all_vars(&(codexcopy[cnt-cnts]), &numvars);
    for (cntvar=0; cntvar<numvars; cntvar++){
      SOAP_add_to_vars(codex, (codexcopy[cnt-cnts]).vars[cntvar].name, (codexcopy[cnt-cnts]).vars[cntvar].value);
    }
    SOAP_free_all_vars(((codexcopy[cnt-cnts]).vars));
    SOAP_free_codex(&(codexcopy[cnt-cnts]));
  }
}


static void BA_if(char **argum, SOAP_codex_t *codex){
  char *ifcode;

  ifcode=(char *)malloc(sizeof(char));
  SOAP_substitute_argum(argum,0,codex);
  if (!(SOAP_number_argums(*argum)==2 || SOAP_number_argums(*argum)==3))
    SOAP_fatal_error(codex,"Not the right number of arguments in if() command; "
                          "the first argument is the condition; "
			  "the second argument is the code to execute if the condition is true; "
			  "the third argument (not required) is the code to execute if the condition if false.");
  if (is_logical(codex,longfromdouble(SOAP_get_argum_double(codex,*argum,0)))) {
    SOAP_get_argum_straight(codex,&ifcode, *argum, 1);
    SOAP_process_code(ifcode, codex, SOAP_VARS_KEEP_ALL);
  } else {
    if (SOAP_number_argums(*argum)==3){
      SOAP_get_argum_straight(codex,&ifcode, *argum, 2);
      SOAP_process_code(ifcode, codex, SOAP_VARS_KEEP_ALL);
    }
  }
  free(ifcode);
}

static void BA_include(char **argum, SOAP_codex_t *codex){
  char *includedcode;
  char *filename;
  char *message;
  char *filename_mem;

  filename_mem=(char *)malloc(sizeof(char)*(2+strlen(codex->filename)));
  strcpy(filename_mem,codex->filename);
  if (SOAP_number_argums(*argum)!=1)
    SOAP_fatal_error(codex,"Not the right number of arguments in include() command; "
                          "the first and only argument is the name of the file to be included.");
  SOAP_substitute_argum(argum,0,codex);

  filename=(char *)malloc(sizeof(char));
  includedcode=(char *)malloc(sizeof(char));
  SOAP_get_argum_string(codex, &filename, *argum, 0);
  SOAP_store_file_as_string(filename, &includedcode);
  SOAP_insert_line_numbers_in_code(&includedcode, 1);
  message=(char *)malloc((100+strlen(filename))*sizeof(char));
  sprintf(message,"%s\n  included on line %ld of file ",filename,codex->linenum);
  SOAP_strins(message,&codex->filename,0);
  SOAP_process_code(includedcode, codex, SOAP_VARS_KEEP_ALL);

  strcpy(codex->filename,filename_mem);

  free(filename_mem);
  free(includedcode);
  free(filename);
  free(message);
}


static void BA_while(char **argum, SOAP_codex_t *codex){
  char *loopcode,*condition;
  bool CONTINUE;


  condition=(char *)malloc(sizeof(char));
  loopcode=(char *)malloc(sizeof(char));
  SOAP_get_argum_straight(codex,&loopcode, *argum, 1);

  CONTINUE=TRUE;
  do {
    SOAP_get_argum_straight(codex,&condition,*argum,0);
    SOAP_substitute_argum(&condition,0,codex);
    if (!is_logical(codex,longfromdouble(SOAP_get_argum_double(codex,condition,0)))) CONTINUE=FALSE;
    if (CONTINUE) SOAP_process_code(loopcode, codex, SOAP_VARS_KEEP_ALL);
  } while (CONTINUE);

  free(condition);
  free(loopcode);

}

static void BA_exit(char **argum, SOAP_codex_t *codex){
  long ret;
  SOAP_substitute_argum(argum,0,codex);
  ret=longfromdouble(SOAP_get_argum_double(codex,*argum,0));
#ifdef DISTMPI
  MPI_Finalize (  );
#endif
  exit(ret);
}


static void BA_system(char **argum, SOAP_codex_t *codex){
  char *expr;
  expr=(char *)malloc(sizeof(char));
  SOAP_substitute_argum(argum,0,codex);
  SOAP_get_argum_string(codex, &expr, *argum, 0);
  if (codex->SYSTEMCALL) {
    if (system(expr)==-1) fprintf(stdout,"Problem executing system command in BA_system()");
  }
}


static void BA_newline(char **argum, SOAP_codex_t *codex){
  SOAP_substitute_argum(argum,0,codex);
  codex->linenum=SOAP_get_argum_long(codex,*argum,0);
}


static void builtin_actions(char *action, char **argum, SOAP_codex_t *codex){
  if (strcmp(action,"write")==0) { 
    BA_write(argum,codex,FALSE); 
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"writeln")==0) {
    BA_write(argum,codex,TRUE);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"printf")==0) {
    BA_printf(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"fprintf")==0) {
    BA_fprintf(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"for")==0) {
    BA_for(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"for_parallel")==0) {
    BA_for_parallel(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"if")==0) {
    BA_if(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"include")==0) {
    BA_include(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"while")==0) {
    BA_while(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"exit")==0) {
    BA_exit(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"system")==0) {
    BA_system(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
  if (strcmp(action,"__newline")==0) {
    BA_newline(argum,codex);
    codex->ACTIONPROCESSED=TRUE; 
  }
}


void SOAP_add_to_vars(SOAP_codex_t *codex, char *name, char *value){
  long varnum;
  bool FOUNDMATCH;

  varnum=0;
  FOUNDMATCH=FALSE;
  while (codex->vars[varnum].name!=NULL) {
    if (strcmp(codex->vars[varnum].name,name)==0) {
      FOUNDMATCH=TRUE;
      codex->vars[varnum].value=(char *)realloc(codex->vars[varnum].value,
                                  ((long)strlen(value)+2)*sizeof(char));
      strcpy(codex->vars[varnum].value,value);
    }
    varnum++;
  }
  if (!FOUNDMATCH) {
    codex->vars=(SOAP_vars_t *)realloc(codex->vars,(varnum+2)*sizeof(SOAP_vars_t));
    codex->vars[varnum].name=(char *)malloc(((long)strlen(name)+2)*sizeof(char));
    codex->vars[varnum].value=(char *)malloc(((long)strlen(value)+2)*sizeof(char));
    strcpy(codex->vars[varnum].name,name);
    strcpy(codex->vars[varnum].value,value);
    codex->vars[varnum+1].name=NULL;
  }
}


void SOAP_add_int_to_vars(SOAP_codex_t *codex, char *name, int value){
  char valuestr[100];
  if (sprintf(valuestr,"%d",value)<0) SOAP_fatal_error(codex,"Problem converting within SOAP_add_int_to_vars(); name=%s  value=%d valuestr=%s.",name,value,valuestr);
  SOAP_add_to_vars(codex,name,valuestr);
}


bool SOAP_is_var_in_codex(SOAP_codex_t *codex, char *name){
  long varnum;
  bool FOUNDMATCH;

  varnum=0;
  FOUNDMATCH=FALSE;
  while (codex->vars[varnum].name!=NULL) {
    if (strcmp(codex->vars[varnum].name,name)==0) {
      FOUNDMATCH=TRUE;
    }
    varnum++;
  }
  return(FOUNDMATCH);
}


double SOAP_var_value(SOAP_codex_t *codex, char *name){
  long varnum;
  bool FOUNDMATCH;
  double value;
  int eos=EOS;
  varnum=0;
  FOUNDMATCH=FALSE;
  while (codex->vars[varnum].name!=NULL) {
    if (strcmp(codex->vars[varnum].name,name)==0) {
      FOUNDMATCH=TRUE;
      if (sscanf(codex->vars[varnum].value,"%lg%n",&value,&eos)!=1 || (codex->vars[varnum].value)[eos]!=EOS)
        SOAP_fatal_error(codex,"Problem evaluating expression >%s<.",codex->vars[varnum].value);
    }
    varnum++;
  }
  if (!FOUNDMATCH) {
    SOAP_fatal_error(codex,"Can't find variable match for %s.",name);
  }
  return(value);
}


void SOAP_var_value_string(SOAP_codex_t *codex, char *name, char **value){
  long varnum;
  bool FOUNDMATCH;

  varnum=0;
  FOUNDMATCH=FALSE;

  while (codex->vars[varnum].name!=NULL) {
    if (strcmp(codex->vars[varnum].name,name)==0) {
      FOUNDMATCH=TRUE;
      *value=(char *)realloc(*value,(strlen(codex->vars[varnum].value)+3)*sizeof(char));
      strcpy(*value,codex->vars[varnum].value);
    }
    varnum++;
  }
  if (!FOUNDMATCH) {
    SOAP_fatal_error(codex,"Can't find variable match for %s.",name);
  }
}



void SOAP_copy_all_vars(SOAP_vars_t *vars1, SOAP_vars_t **vars2){
  long varnum;
  varnum=0;
  *vars2=(SOAP_vars_t *)realloc(*vars2,2*sizeof(SOAP_vars_t));
  (*vars2)[varnum].name=NULL;
  while (vars1[varnum].name!=NULL) {
    *vars2=(SOAP_vars_t *)realloc(*vars2,(varnum+2)*sizeof(SOAP_vars_t));
    (*vars2)[varnum].name=(char *)malloc(((long)strlen(vars1[varnum].name)+2)*sizeof(char));
    (*vars2)[varnum].value=(char *)malloc(((long)strlen(vars1[varnum].value)+2)*sizeof(char));
    strcpy((*vars2)[varnum].name,vars1[varnum].name);
    strcpy((*vars2)[varnum].value,vars1[varnum].value);
    (*vars2)[varnum+1].name=NULL;
    varnum++;
  }  
}


void SOAP_free_all_vars(SOAP_vars_t *vars){
  long varnum;
  varnum=0;
  while (vars[varnum].name!=NULL) {
    free(vars[varnum].name);
    free(vars[varnum].value);
    varnum++;
  }
  vars[0].name=NULL;
}


void SOAP_free_codex_copy(SOAP_codex_t *codex){
  free(codex->filename);
  if (codex->action_being_processed!=NULL) free(codex->action_being_processed);
}


void SOAP_free_codex(SOAP_codex_t *codex){
  long varnum;
  varnum=0;
  while (codex->vars[varnum].name!=NULL) {
    free(codex->vars[varnum].name);
    free(codex->vars[varnum].value);
    varnum++;
  }
  codex->vars[0].name=NULL;
  free(codex->vars);
  free(codex->filename);
  if (codex->action_being_processed!=NULL) free(codex->action_being_processed);
}


void SOAP_init_codex(SOAP_codex_t *codex, const char *filename){
  codex->vars=(SOAP_vars_t *)malloc(sizeof(SOAP_vars_t));
  codex->vars[0].name=NULL;
  codex->VERBOSE=FALSE;
  codex->FUNCTION=FALSE;
  codex->ACTION=FALSE;
  codex->SCREENOUTPUT=TRUE;
  codex->FILEOUTPUT=TRUE;
  codex->SYSTEMCALL=TRUE;
  codex->ACTIONPROCESSED=FALSE;
  codex->linenum=0;
  codex->filename=(char *)malloc((2+strlen(filename))*sizeof(char));
  strcpy(codex->filename,filename);
  codex->action_being_processed=NULL;
}


void SOAP_copy_codex(SOAP_codex_t *orig, SOAP_codex_t *copy){
  copy->vars=orig->vars;
  copy->VERBOSE=orig->VERBOSE;
  copy->FUNCTION=orig->FUNCTION;
  copy->ACTION=orig->ACTION;
  copy->SCREENOUTPUT=orig->SCREENOUTPUT;
  copy->FILEOUTPUT=orig->FILEOUTPUT;
  copy->SYSTEMCALL=orig->SYSTEMCALL;
  copy->ACTIONPROCESSED=orig->ACTIONPROCESSED;
  copy->action=orig->action;
  copy->function=orig->function;
  copy->action_args=orig->action_args;
  copy->function_args=orig->function_args;
  copy->linenum=orig->linenum;
  copy->filename=(char *)malloc((strlen(orig->filename)+2)*sizeof(char));
  strcpy(copy->filename,orig->filename);
  if (orig->action_being_processed!=NULL) {
    copy->action_being_processed=(char *)malloc((strlen(orig->action_being_processed)+2)*sizeof(char));
    strcpy(copy->action_being_processed,orig->action_being_processed);
  } else {
    copy->action_being_processed=NULL;
  }
}


void SOAP_count_all_vars(SOAP_codex_t *codex, long *numvars){
  long NULL_POS;
  NULL_POS=-1;
  do {
    NULL_POS++;
  } while ((codex->vars)[NULL_POS].name!=NULL);
  *numvars=NULL_POS;
}


void SOAP_clean_added_vars(SOAP_codex_t *codex, long numvarsinit){
  long numvars,cnt;
  SOAP_count_all_vars(codex, &numvars);
  for (cnt=numvarsinit; cnt<numvars; cnt++) {
    free((codex->vars)[cnt].name);
    free((codex->vars)[cnt].value);
  }
  (codex->vars)[numvarsinit].name=NULL;
}


/* process a piece of code defined in *code
   with the actions list defined in *action */
void SOAP_process_code(char *code, SOAP_codex_t *codex, int SOAP_VARS){
  long cnt,cnt2,codelength,numvarsinit,parentheses;
  bool ASSIGN,INSTRING;
  char *action,*argum;

  SOAP_count_all_vars(codex,&numvarsinit);

  clean_code(codex,code);
  /* ShowCode(code);   */
  codelength=(long)strlen(code);
  if (codelength>0) {
    cnt=0;
    action=(char *)malloc(sizeof(char));
    argum=(char *)malloc(sizeof(char));
    do {
      /* get action string and argument string*/
      cnt2=0;
      while (code[cnt]!='(' && code[cnt]!='=' && code[cnt]!=EOS && code[cnt]!=';' && code[cnt]!='"') {
        action=(char *)realloc(action,(cnt2+2)*sizeof(char));
        action[cnt2]=code[cnt];
        cnt++;
        cnt2++;
      }
      if (code[cnt]==EOS || code[cnt]==';' || code[cnt]=='"') {
        SOAP_fatal_error(codex,"Action name not followed by '(' or '='.");
      }
      if (code[cnt]=='=') ASSIGN=TRUE; else ASSIGN=FALSE;
      action[cnt2]=EOS;
      cnt++;
      cnt2=0;
      if (ASSIGN) parentheses=0; else parentheses=1;
      INSTRING=FALSE;
      while ((code[cnt]!=';' || parentheses>0 || INSTRING) && code[cnt]!=EOS) {
        argum=(char *)realloc(argum,(cnt2+2)*sizeof(char));
        argum[cnt2]=code[cnt];
        if (code[cnt]=='"') INSTRING=!INSTRING;
        if (code[cnt]=='(') parentheses++;
        if (code[cnt]==')') parentheses--;
        cnt++;
        cnt2++;
      }
      if (!ASSIGN){
        codex->ACTIONPROCESSED=FALSE;
        codex->action_being_processed=realloc(codex->action_being_processed,sizeof(char) * (2+strlen(action))); 
        strcpy(codex->action_being_processed,action);
      }
      if (parentheses>0) {
        SOAP_fatal_error(codex,"Missing ')' .");
      }
      if (parentheses<0) {
        SOAP_fatal_error(codex,"Too many ')' .");
      }
      if (code[cnt]==EOS) {
        SOAP_fatal_error(codex,"Expecting ';' .");
      }
      cnt++;
      if (ASSIGN) argum[cnt2]=EOS; else argum[cnt2-1]=EOS;
      /* if (codex->VERBOSE) printf("action='%s'  argum='%s'\n",action,argum); */
      if (ASSIGN) {
        substitute_functions(&argum, codex);
        update_var(&action,&argum,codex);
      } else {
        if (codex->ACTION) (codex->action)(action,&argum,codex);
        builtin_actions(action,&argum,codex);
        if (!codex->ACTIONPROCESSED && codex->SCREENOUTPUT) fprintf(stdout,"%s ignored..",action);
        free(codex->action_being_processed);
        codex->action_being_processed=NULL;
      }
    } while (cnt<codelength);
    free(action);
    free(argum);
  }
  if (SOAP_VARS==SOAP_VARS_CLEAN_ADDED) {
    SOAP_clean_added_vars(codex, numvarsinit);
  }
  if (SOAP_VARS==SOAP_VARS_CLEAN_ALL) {
    fprintf(stdout,"SOAP_VARS_CLEAN_ALL not yet implemented in soap.c.\n");
  }
}



