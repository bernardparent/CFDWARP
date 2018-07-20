#ifndef SOAP_H
#define SOAP_H

#include <memdebug.h>
#include <exm.h>

#define SOAP_VARS_KEEP_ALL 0
#define SOAP_VARS_CLEAN_ADDED 1
#define SOAP_VARS_CLEAN_ALL 2

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"



  typedef struct {
    char *value;
    char *name;
  } SOAP_vars_t;

 /* the code accessories */
  typedef struct SOAP_codex_t {
    SOAP_vars_t *vars;
    bool VERBOSE,ACTION,FUNCTION,SCREENOUTPUT,FILEOUTPUT,SYSTEMCALL,ACTIONPROCESSED;
    void (*action) (char *, char **, struct SOAP_codex_t *);
    void *action_args;
    void (*function) (char *, char **, char **, struct SOAP_codex_t *);
    void *function_args;
    long linenum;     /* this is used only when an error is encountered */
    char *filename; /* this is used only when an error is encountered */
    char *action_being_processed;
  } SOAP_codex_t;

  /* cut _all characters from is to ie in *str */
  void SOAP_strcut(long is, long ie, char *str);

  /* insert str1 into str2 at the position pos */
  void SOAP_strins(char *str1, char **str2,  long pos);

  void SOAP_store_file_as_string(char *filename, char **str);

  /* evaluate a string expression in which there are
   parentheses "(10*50)/40^10.0E7" */
  double SOAP_evaluate_arithmetic(SOAP_codex_t *codex, char *expr_orig);

  /* get the nth argument; arguments are counted from 0*/
  double SOAP_get_argum_double(SOAP_codex_t *codex, char *argum, long n);

  /* get the nth argument and store it in *expr;
     arguments are counted from 0. *expr must already
     have been malloc'ed*/
  void SOAP_get_argum_straight(SOAP_codex_t *codex, char **expr, char *argum, long n);

  /* get the nth argument; arguments are counted from 0*/
  long SOAP_get_argum_long(SOAP_codex_t *codex, char *argum, long n);


  /* get the nth argument; arguments are counted from 0*/
  long SOAP_get_argum_bool(SOAP_codex_t *codex, char *argum, bool n);

  /* get the nth argument and store what is in between
     the quotes of the string in *expr;
     arguments are counted from 0. *expr must already
     have been malloc'ed*/
  void SOAP_get_argum_string(SOAP_codex_t *codex, char **expr, char *argum, long n);

  /* evaluate the expr **expr and substitute it for
   the number it stands for */
  void SOAP_substitute_expression(char **expr, SOAP_codex_t *codex);

  /* sub the expression located in nth argument
   with SOAP_substitute_expression; arguments are counted from 0*/
  void SOAP_substitute_argum(char **argum, long n, SOAP_codex_t *codex);

  /* return the number of arguments */
  long SOAP_number_argums(char *argum);

  /* sub expressions of _all the arguments  */
  void SOAP_substitute_all_argums(char **argum, SOAP_codex_t *codex);

  /* process a piece of code defined in *code
   with the actions list defined in *funct */
  void SOAP_process_code(char *code, SOAP_codex_t *codex, int SOAP_VARS);

  /* prior to SOAP_process_code, variables can be added from
     within the C program */
  void SOAP_add_to_vars(SOAP_codex_t *codex, char *name, char *value);

  /* prior to SOAP_process_code, int variables can be added from
     within the C program */
  void SOAP_add_int_to_vars(SOAP_codex_t *codex, char *name, int value);

  /* returns the value of the variable with name *name */
  double SOAP_var_value(SOAP_codex_t *codex, char *name);

  /* returns the value of the variable with name *name  in string *value */
  void SOAP_var_value_string(SOAP_codex_t *codex, char *name, char **value);

  void SOAP_free_all_vars(SOAP_vars_t *vars);

  void SOAP_copy_all_vars(SOAP_vars_t *vars_src, SOAP_vars_t **vars_dest);

  void SOAP_count_all_vars(SOAP_codex_t *codex, long *numvars);

  void SOAP_clean_added_vars(SOAP_codex_t *codex, long numvarsinit);

  bool SOAP_is_var_in_codex(SOAP_codex_t *codex, char *name);


  /* to be called just prior to executing SOAP_process_code but before
     SOAP_add_to_vars, filename being used only when an error is encountered:
     if filename is not known, just set it to ? */
  void SOAP_init_codex(SOAP_codex_t *codex, const char *filename);

  /* Makes a perfect copy of  *orig to *copy */
  void SOAP_copy_codex(SOAP_codex_t *orig, SOAP_codex_t *copy);

  /* to be called after SOAP_process_code */
  void SOAP_free_codex(SOAP_codex_t *codex);

  void SOAP_insert_line_numbers_in_code(char **code, long linenum_start);

  int SOAP_fatal_error(SOAP_codex_t *codex, const char *formatstr, ...);

  bool SOAP_is_double_a_long(double expr_double);

#endif

