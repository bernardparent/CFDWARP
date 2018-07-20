#include "model.h"
#include <src/control.h>
#include <model/thermo/_thermo.h>
#include <model/fluid/_fluid.h>
#include <model/emfield/_emfield.h>
#include <model/beam/_beam.h>


void find_clipped_variables_list(gl_t *gl, char **cliplist){
  char strtmp[300];
  long cnt;
  strcpy(*cliplist,"");
  /* sprintf(strtmp," clipnumtot(%ld)",gl->model.clipnumtot);
  SOAP_strins(strtmp, cliplist,  0);*/
  for (cnt=0; cnt<gl->model.clipnamenum; cnt++){
    sprintf(strtmp," %s(%ld)",gl->model.clipname[cnt],gl->model.clipnum[cnt]);    
    if (strstr(strtmp, "muscl") == NULL && strstr(strtmp, "bdry") == NULL) {
      SOAP_strins(strtmp, cliplist,  strlen(*cliplist));
    }
  }
}


void find_clipped_muscl_variables_list(gl_t *gl, char **cliplist){
  char strtmp[300];
  long cnt;
  strcpy(*cliplist,"");
  for (cnt=0; cnt<gl->model.clipnamenum; cnt++){
    sprintf(strtmp," %s(%ld)",gl->model.clipname[cnt],gl->model.clipnum[cnt]);
    if (strstr(strtmp, "muscl") != NULL) {
      SOAP_strins(strtmp, cliplist,  strlen(*cliplist));
    }
  }
}


void find_clipped_bdry_variables_list(gl_t *gl, char **cliplist){
  char strtmp[300];
  long cnt;
  strcpy(*cliplist,"");
  for (cnt=0; cnt<gl->model.clipnamenum; cnt++){
    sprintf(strtmp," %s(%ld)",gl->model.clipname[cnt],gl->model.clipnum[cnt]);    
    if (strstr(strtmp, "bdry") != NULL) {
      SOAP_strins(strtmp, cliplist,  strlen(*cliplist));
    }
  }
}


void add_to_clipped_variables(gl_t *gl, char *str){
  long cnt;
  bool CLIPFOUND;
  if (gl->REPORTCLIPPING){
    thread_lock_set(&(gl->model.clip_lock),THREADTYPE_ALL);
    gl->model.clipnumtot++;
    CLIPFOUND=FALSE;
    for (cnt=0; cnt<gl->model.clipnamenum; cnt++){
      if (strcmp(str,gl->model.clipname[cnt])==0){
        gl->model.clipnum[cnt]++;
        CLIPFOUND=TRUE;
      }
    }
    if (!CLIPFOUND){
      gl->model.clipname=(clipname_t *)
       realloc(gl->model.clipname,(gl->model.clipnamenum+1)*sizeof(clipname_t));
      gl->model.clipnum=(clipnum_t *)
       realloc(gl->model.clipnum,(gl->model.clipnamenum+1)*sizeof(clipnum_t));
      gl->model.clipname[gl->model.clipnamenum]=(char *)malloc((1+strlen(str))*sizeof(char));
      strcpy(gl->model.clipname[gl->model.clipnamenum],str);
      gl->model.clipnum[gl->model.clipnamenum]=1;
      gl->model.clipnamenum++;
    }
    thread_lock_unset(&(gl->model.clip_lock),THREADTYPE_ALL);
  }
}


void add_to_clipped_variables2(gl_t *gl, char *str, char *suffix){
  char catstr[1500];
  if (gl->REPORTCLIPPING){
    catstr[0]=EOS;
    strcat(catstr,str);
    strcat(catstr,suffix);
    add_to_clipped_variables(gl,catstr);
  }
}


#ifdef DISTMPI
void find_clipped_variables_list_all(gl_t *gl, char **cliplist){
  char strtmp[300];
  long cnt;
  strcpy(*cliplist,"");
  for (cnt=0; cnt<gl->model.clipnamenum_all; cnt++){
    sprintf(strtmp," %s(%ld)",gl->model.clipname_all[cnt],gl->model.clipnum_all[cnt]);    
    if (strstr(strtmp, "muscl") == NULL && strstr(strtmp, "bdry") == NULL) {
      SOAP_strins(strtmp, cliplist,  strlen(*cliplist));
    }
  }
}


void find_clipped_muscl_variables_list_all(gl_t *gl, char **cliplist){
  char strtmp[300];
  long cnt;
  strcpy(*cliplist,"");
  for (cnt=0; cnt<gl->model.clipnamenum_all; cnt++){
    sprintf(strtmp," %s(%ld)",gl->model.clipname_all[cnt],gl->model.clipnum_all[cnt]);
    if (strstr(strtmp, "muscl") != NULL) {
      SOAP_strins(strtmp, cliplist,  strlen(*cliplist));
    }
  }
}


void find_clipped_bdry_variables_list_all(gl_t *gl, char **cliplist){
  char strtmp[300];
  long cnt;
  strcpy(*cliplist,"");
  for (cnt=0; cnt<gl->model.clipnamenum_all; cnt++){
    sprintf(strtmp," %s(%ld)",gl->model.clipname_all[cnt],gl->model.clipnum_all[cnt]);    
    if (strstr(strtmp, "bdry") != NULL) {
      SOAP_strins(strtmp, cliplist,  strlen(*cliplist));
    }
  }
}


void add_to_clipped_variables_all(gl_t *gl, char *str, long clipnum){
  long cnt;
  bool CLIPFOUND;
  thread_lock_set(&(gl->model.clip_lock),THREADTYPE_ALL);
  gl->model.clipnumtot_all+=clipnum;
  CLIPFOUND=FALSE;
  for (cnt=0; cnt<gl->model.clipnamenum_all; cnt++){
    if (strcmp(str,gl->model.clipname_all[cnt])==0){
      gl->model.clipnum_all[cnt]+=clipnum;
      CLIPFOUND=TRUE;
    }
  }
  if (!CLIPFOUND){
    gl->model.clipname_all=(clipname_t *)
       realloc(gl->model.clipname_all,(gl->model.clipnamenum_all+1)*sizeof(clipname_t));
    gl->model.clipnum_all=(clipnum_t *)
       realloc(gl->model.clipnum_all,(gl->model.clipnamenum_all+1)*sizeof(clipnum_t));
    gl->model.clipname_all[gl->model.clipnamenum_all]=(char *)malloc((1+strlen(str))*sizeof(char));
    strcpy(gl->model.clipname_all[gl->model.clipnamenum_all],str);
    gl->model.clipnum_all[gl->model.clipnamenum_all]=clipnum;
    gl->model.clipnamenum_all++;
  }
  thread_lock_unset(&(gl->model.clip_lock),THREADTYPE_ALL);
}


void reset_clipped_variables_all(gl_t *gl){
  long cnt;
  thread_lock_set(&(gl->model.clip_lock),THREADTYPE_ALL);
  gl->model.clipnumtot_all=0;
  for (cnt=0; cnt<gl->model.clipnamenum_all; cnt++) {
    free(gl->model.clipname_all[cnt]);
  }
  gl->model.clipnamenum_all=0;
  gl->model.clipname_all=(clipname_t *)realloc(gl->model.clipname_all,sizeof(clipname_t));
  gl->model.clipnum_all=(clipnum_t *)realloc(gl->model.clipnum_all,sizeof(clipnum_t));
  thread_lock_unset(&(gl->model.clip_lock),THREADTYPE_ALL);
}

#endif


void init_clipped_variables(gl_t *gl){
  thread_lock_init(&(gl->model.clip_lock),THREADTYPE_ALL);
  thread_lock_set(&(gl->model.clip_lock),THREADTYPE_ALL);
  gl->model.clipnumtot=0;
  gl->model.clipnamenum=0;
  gl->model.clipname=(clipname_t *)malloc(sizeof(clipname_t));
  gl->model.clipnum=(clipnum_t *)malloc(sizeof(clipnum_t));
#ifdef DISTMPI
  gl->model.clipnumtot_all=0;
  gl->model.clipnamenum_all=0;
  gl->model.clipname_all=(clipname_t *)malloc(sizeof(clipname_t));
  gl->model.clipnum_all=(clipnum_t *)malloc(sizeof(clipnum_t));
#endif
  thread_lock_unset(&(gl->model.clip_lock),THREADTYPE_ALL);
}


void reset_clipped_variables(gl_t *gl){
  long cnt;
  thread_lock_set(&(gl->model.clip_lock),THREADTYPE_ALL);
  gl->model.clipnumtot=0;
  for (cnt=0; cnt<gl->model.clipnamenum; cnt++) {
    free(gl->model.clipname[cnt]);
  }
  gl->model.clipnamenum=0;
  gl->model.clipname=(clipname_t *)realloc(gl->model.clipname,sizeof(clipname_t));
  gl->model.clipnum=(clipnum_t *)realloc(gl->model.clipnum,sizeof(clipnum_t));
/*  gl->model.num_Pmin=0; */
  thread_lock_unset(&(gl->model.clip_lock),THREADTYPE_ALL);
}


void free_clipped_variables(gl_t *gl){
  reset_clipped_variables(gl);
  free(gl->model.clipname);
  free(gl->model.clipnum);
#ifdef DISTMPI
  reset_clipped_variables_all(gl);
  free(gl->model.clipname_all);
  free(gl->model.clipnum_all);
#endif
}


void write_model_template(FILE **controlfile){
  wfprintf(*controlfile,
    " \n"
    " \n"
    "Model(\n"
  );
  write_model_fluid_template(controlfile);
  write_model_emfield_template(controlfile);
  write_model_beam_template(controlfile);
  wfprintf(*controlfile,
    ");\n"
  );
}


void read_model_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  read_model_fluid_actions(actionname, argum, codex);
  read_model_emfield_actions(actionname, argum, codex);
  read_model_beam_actions(actionname, argum, codex);
}


void read_model(char *argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (!gl->CONTROL_READ){
    gl->MODEL_FLUID_READ=FALSE;
    gl->MODEL_EMFIELD_READ=FALSE;
    gl->MODEL_BEAM_READ=FALSE;
    init_clipped_variables(gl);
  } 
  codex->action=&read_model_actions;
  SOAP_process_code(argum, codex, SOAP_VARS_CLEAN_ADDED);
  read_model_thermo();

  if (!gl->MODEL_FLUID_READ) 
    fatal_error("The fluid module %s was not found within Model().",_FLUID_ACTIONNAME);
  if (!gl->MODEL_EMFIELD_READ) 
    fatal_error("The emfield module %s was not found within Model().",_EMFIELD_ACTIONNAME);
  if (!gl->MODEL_BEAM_READ) 
    fatal_error("The beam module %s was not found within Model().",_BEAM_ACTIONNAME);

}





void find_post_proc_var_name(long varnum, char *varname){
  long varnumx;
  if (varnum<totalpostvarfluid) {
    find_post_variable_name_fluid(varnum, varname);
  } 
  varnumx=totalpostvarfluid;
#ifdef EMFIELD
  if (varnum>=totalpostvarfluid && varnum<totalpostvarfluid+totalpostvaremfield) {
    find_post_variable_name_emfield(varnum-totalpostvarfluid, varname);  
  }
  varnumx=totalpostvarfluid+totalpostvaremfield;
#endif
  if (varnum==varnumx) sprintf(varname,"x");
  if (varnum==varnumx+1) sprintf(varname,"y");
  if (varnum==varnumx+2) sprintf(varname,"z");
}


void find_post_proc_var_value(np_t *np, long l, gl_t *gl, long varnum, double *varvalue){
  long varnumx;
  *varvalue=0.0;
  if (varnum<totalpostvarfluid) {
    if (is_node_valid(np[l],TYPELEVEL_FLUID)) {
      find_post_variable_value_fluid(np, l, gl, varnum, varvalue);
    }
  } 
  varnumx=totalpostvarfluid;
#ifdef EMFIELD
  if (varnum>=totalpostvarfluid && varnum<totalpostvarfluid+totalpostvaremfield) {
    if (is_node_valid(np[l],TYPELEVEL_EMFIELD)) {
      find_post_variable_value_emfield(np, l, gl, varnum-totalpostvarfluid, varvalue);
    }
  }
  varnumx=totalpostvarfluid+totalpostvaremfield;
#endif
  if (varnum==varnumx) *varvalue=_x(np[l],0);
  if (varnum==varnumx+1) *varvalue=_x(np[l],1);
  if (varnum==varnumx+2) *varvalue=_x(np[l],2);
}


bool is_node_bdry_symmetry_plane(np_t np){
  double RET;
  if (is_node_bdry_symmetry_plane_fluid(np) 
#ifdef EMFIELD
   && !is_node_inner(np,TYPELEVEL_EMFIELD)
#endif
  ) RET=TRUE; else RET=FALSE;
  return(RET);
}

