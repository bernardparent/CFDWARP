// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015-2016, 2020 Bernard Parent

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

#include <src/init.h>
#include <model/_model.h>
#include <src/control.h>
#include <cycle/share/cycle_share.h>

#define INITSPECIES_DEFAULT -1

static void reorder_initvar_species(gl_t *gl, initvar_t initvar, long specstart){
  long flux,spec,specinit;
  spec_t N;
  if (gl->nsinit!=ns){
    // move values up to make space for missing species
    for (flux=numinitvar-1; flux>=specstart+ns; flux--){
      assert(flux>0);
      assert(flux<numinitvar);
      assert(flux-(ns-gl->nsinit)>=0);
      //printf("flux=%ld  ns=%ld  gl->nsinit=%ld\n",flux,ns,gl->nsinit);
      initvar[flux]=initvar[flux-(ns-gl->nsinit)]; 
    } 
    // set all densities to the default
    for (spec=0; spec<ns; spec++) N[spec]=initvar[specstart+gl->nsinit-1];
    
    // set the given densities
    for (specinit=0; specinit<gl->nsinit-1; specinit++){
      assert(gl->initspecies[specinit]<ns);
      assert(gl->initspecies[specinit]>=0);
      N[gl->initspecies[specinit]]=initvar[specstart+specinit];
    }
    
    //put the densities back into initvar
    for (spec=0; spec<ns; spec++){
      initvar[specstart+spec]=N[spec];
    }
  }
}


void reformat_initvar_densities(gl_t *gl, initvar_t initvar_orig, initvar_t initvar, long specstart){
  long flux;
  
  for (flux=0; flux<numinitvar; flux++) initvar[flux]=initvar_orig[flux];
  reorder_initvar_species(gl, initvar, specstart);
}  


void reformat_initvar_species_fractions(gl_t *gl, initvar_t initvar_orig, initvar_t initvar, long specstart){
  long spec,flux;
  double sum;
  
  for (flux=0; flux<numinitvar; flux++) initvar[flux]=initvar_orig[flux];
  reorder_initvar_species(gl, initvar, specstart);
    
  sum=0.0;
  for (spec=0; spec<ns; spec++){
    sum+=initvar[specstart+spec];
  }
  if (sum<0.999 || sum>1.001){
    fatal_error("Species mass or mole fractions sum to %E and not to 1 in Init() module.",sum);
  } else {
    for (spec=0; spec<ns; spec++){
      initvar[specstart+spec]/=sum;
    }
  }
}

void verify_positivity_of_determinative_property(initvar_t initvar, long row_start, long row_end){
  long row; 

  for (row=row_start; row<=row_end; row++){
    if (initvar[row]<=0.0) {
#ifdef _RESCONV_POSITIVITY_PRESERVING
      fatal_error("Init variable #%ld must be set to a positive value when a positivity-preserving spatial discretization method is specified. It is now set to %E.", row+1,initvar[row]);
#endif
    }
  }
}


void write_init_template(FILE **controlfile){
  wfprintf(*controlfile,"\n\n");
  wfprintf(*controlfile,
  "Init(\n");
  write_init_fluid_template(controlfile);
#ifdef EMFIELD
  write_init_emfield_template(controlfile);
#endif
  wfprintf(*controlfile,
  ");\n");
}


#ifdef EMFIELD

void read_init_actions_emfield(char *actionname, char **argum, SOAP_codex_t *codex){
  np_t **np;
  gl_t *gl;
  long inittype,bdrytype,flux,i,j,k,is,ie;
  #ifdef _2DL
    long js,je;
  #endif
  #ifdef _3DL
    long ks,ke;
  #endif
  initvar_emfield_t values;
  double dummy;
  zone_t zone;
  np=((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;

  
  if (strcmp(actionname,"All")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=numinitvar_emfield+1)
      SOAP_fatal_error(codex,"Number of arguments not equal to %d in All() command part of Init().",numinitvar_emfield+1);
    inittype=SOAP_get_argum_long(codex,*argum,0);
    if (inittype<1 || inittype>totalinitvartypeemfield) 
      SOAP_fatal_error(codex,"Initial condition type %ld is out of bounds.",inittype);
    for (flux=0; flux<numinitvar_emfield; flux++) values[flux]=SOAP_get_argum_double(codex,*argum,flux+1);
    for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          if (is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)
            ) {
            init_node_emfield((*np)[_ai(gl,i,j,k)], gl, inittype, values);
            (*np)[_ai(gl,i,j,k)].INIT_EMFIELD=TRUE;
          }
    }
    codex->ACTIONPROCESSED=TRUE;
  }

  if (strcmp(actionname,"Bdry")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=numinitvar_emfield+2)
      SOAP_fatal_error(codex,"Number of arguments not equal to %d in Bdry() command part of Init().",numinitvar_emfield+2);
    bdrytype=SOAP_get_argum_long(codex,*argum,0);
    inittype=SOAP_get_argum_long(codex,*argum,1);
    if (inittype<1 || inittype>totalinitvartypeemfield) 
      SOAP_fatal_error(codex,"Initial condition type %ld is out of bounds.",inittype);
    for (flux=0; flux<numinitvar_emfield; flux++)
      values[flux]=SOAP_get_argum_double(codex,*argum,flux+2);
    for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          if (   (is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD))
              && (_node_type((*np)[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)==bdrytype)
             ){
            init_node_emfield((*np)[_ai(gl,i,j,k)], gl, inittype, values);
            (*np)[_ai(gl,i,j,k)].INIT_EMFIELD=TRUE;
          }
    }
    codex->ACTIONPROCESSED=TRUE;
  }


  if (strcmp(actionname,"Region")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=numinitvar_emfield+2*nd+1)
      SOAP_fatal_error(codex,"Number of arguments not equal to %d in Region() command part of Init().",numinitvar_emfield+2*nd+1);
    if if1D( (sscanf(*argum,"%ld,%ld,%ld,%lg",
          &is,&ie,&inittype,&dummy)!=4) )
       if2D( (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%lg",
          &is,&js,&ie,&je,&inittype,&dummy)!=6) )
       if3D( (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%ld,%lg",
          &is,&js,&ks,&ie,&je,&ke,&inittype,&dummy)!=8) )
       {
       SOAP_fatal_error(codex,"One or more argument(s) could not be read properly "
                              "in the Region() command part of Init(). Arguments: %s.",*argum);
    }
    if (inittype<1 || inittype>totalinitvartypeemfield) 
      SOAP_fatal_error(codex,"Initial condition type %ld is out of bounds.",inittype);
    for (flux=0; flux<numinitvar_emfield; flux++)
      values[flux]=SOAP_get_argum_double(codex,*argum,flux+1+nd*2);
    find_zone_from_argum(*argum, 0, gl, codex, &zone);
    for_ijk(zone,is,js,ks,ie,je,ke){
          if (is_node_in_zone(i, j, k, gl->domain_lim)){
            if (is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL_EMFIELD)){
              init_node_emfield((*np)[_ai(gl,i,j,k)], gl, inittype, values);
              (*np)[_ai(gl,i,j,k)].INIT_EMFIELD=TRUE;
            }
          }
    }
    codex->ACTIONPROCESSED=TRUE;
  }


}

#endif



void read_init_actions_fluid(char *actionname, char **argum, SOAP_codex_t *codex){
  np_t **np;
  gl_t *gl;
  long inittype,bdrytype,flux,i,j,k,is,ie;
  #ifdef _2DL
    long js,je;
  #endif
  #ifdef _3DL
    long ks,ke;
  #endif
  initvar_t values;
  double dummy;
  zone_t zone;
  long spec,nsinit,specinit;
  char *specnameinit,*specname;
  bool FOUND;

  np=((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (strcmp(actionname,"All")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=numinitvar+1-(ns-gl->nsinit))
      SOAP_fatal_error(codex,"Number of arguments not equal to %d in All() command part of Init().",numinitvar+1-(ns-gl->nsinit));
    inittype=SOAP_get_argum_long(codex,*argum,0);
    if (inittype<1 || inittype>totalinitvartypefluid) 
      SOAP_fatal_error(codex,"Initial condition type %ld is out of bounds.",inittype);
    for (flux=0; flux<numinitvar-(ns-gl->nsinit); flux++) values[flux]=SOAP_get_argum_double(codex,*argum,flux+1);
    for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          if (is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL_FLUID)
            ) {
            init_node_fluid((*np),_ai(gl,i,j,k), gl, inittype, values);
            (*np)[_ai(gl,i,j,k)].INIT_FLUID=TRUE;
          }
    }
    codex->ACTIONPROCESSED=TRUE;
  }

  if (strcmp(actionname,"Bdry")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=numinitvar+2-(ns-gl->nsinit))
      SOAP_fatal_error(codex,"Number of arguments not equal to %d in Bdry() command part of Init().",numinitvar+2-(ns-gl->nsinit));
    bdrytype=SOAP_get_argum_long(codex,*argum,0);
    inittype=SOAP_get_argum_long(codex,*argum,1);
    if (inittype<1 || inittype>totalinitvartypefluid) 
      SOAP_fatal_error(codex,"Initial condition type %ld is out of bounds.",inittype);
    for (flux=0; flux<numinitvar-(ns-gl->nsinit); flux++)
      values[flux]=SOAP_get_argum_double(codex,*argum,flux+2);
    for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          if (   (is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL_FLUID))
              && (_node_type((*np)[_ai(gl,i,j,k)],TYPELEVEL_FLUID)==bdrytype)
             ){
            init_node_fluid((*np),_ai(gl,i,j,k), gl, inittype, values);
            (*np)[_ai(gl,i,j,k)].INIT_FLUID=TRUE;
          }
    }
    codex->ACTIONPROCESSED=TRUE;
  }


  if (strcmp(actionname,"Region")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=numinitvar+2*nd+1-(ns-gl->nsinit))
      SOAP_fatal_error(codex,"Number of arguments not equal to %d in Region() command part of Init().",numinitvar+2*nd+1-(ns-gl->nsinit));
    if (if1D( (sscanf(*argum,"%ld,%ld,%ld,%lg",
          &is,&ie,&inittype,&dummy)!=4) )
       if2D( (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%lg",
          &is,&js,&ie,&je,&inittype,&dummy)!=6) )
       if3D( (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%ld,%lg",
          &is,&js,&ks,&ie,&je,&ke,&inittype,&dummy)!=8) )
       ){
       SOAP_fatal_error(codex,"One or more argument(s) could not be read properly "
                              "in Region() command part of Init(). Arguments: %s.",*argum);
    }
    if (inittype<1 || inittype>totalinitvartypefluid) 
      SOAP_fatal_error(codex,"Initial condition type %ld is out of bounds.",inittype);

    for (flux=0; flux<numinitvar-(ns-gl->nsinit); flux++)
      values[flux]=SOAP_get_argum_double(codex,*argum,flux+1+nd*2);
    find_zone_from_argum(*argum, 0, gl, codex, &zone);
    for_ijk(zone,is,js,ks,ie,je,ke){
          if (is_node_in_zone(i, j, k, gl->domain_lim)){
            if (is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL_FLUID)) {
              init_node_fluid((*np),_ai(gl,i,j,k), gl, inittype, values);
              (*np)[_ai(gl,i,j,k)].INIT_FLUID=TRUE;
            }
          }
    }
    codex->ACTIONPROCESSED=TRUE;
  }
  
  
  if (ns>1 && strcmp(actionname,"Species")==0) {
    nsinit=SOAP_number_argums(*argum);
    if (nsinit>ns) SOAP_fatal_error(codex,"Number of species within Species() can not be higher than the number of species in the chemical solver.");
    gl->nsinit=nsinit;
    specnameinit=(char *)malloc(sizeof(char));
    specname=(char *)malloc(sizeof(char));
    for (specinit=0; specinit<nsinit; specinit++) {
      SOAP_get_argum_string(codex, &specnameinit, *argum, specinit);
      FOUND=FALSE;
      for (spec=0; spec<ns; spec++){
        find_species_name(spec, &specname);
        if (strcmp(specname,specnameinit)==0){
          gl->initspecies[specinit]=spec;
          FOUND=TRUE; 
        }
      }
      if (strcmp("default",specnameinit)==0){
        gl->initspecies[specinit]=INITSPECIES_DEFAULT;
        FOUND=TRUE; 
      }
      if (!FOUND) {
        SOAP_fatal_error(codex,"Species %s invalid within Species(). Set it to one of the species within the chemical solver or to 'default'.",specnameinit);
      }
    }
    if (gl->nsinit<ns && gl->initspecies[gl->nsinit-1]!=INITSPECIES_DEFAULT) SOAP_fatal_error(codex,"The last species within Species() must be set to 'default'.");

    free(specnameinit);
    free(specname);    
    codex->ACTIONPROCESSED=TRUE;
  }

  
}



void read_init_actions(char *actionname, char **argum, SOAP_codex_t *codex){
  gl_t *gl;
#ifdef UNSTEADY
  np_t **np;
  long i,j,k,flux;
  np=((readcontrolarg_t *)codex->action_args)->np;
#endif
  long spec;
  gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_FLUID_ACTIONNAME)==0 ){
    gl->INIT_FLUID_READ=TRUE;
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE)  wfprintf(stdout,"%s..",_FLUID_ACTIONNAME);
    gl->nsinit=ns;
    for (spec=0; spec<ns; spec++) gl->initspecies[spec]=spec;
    codex->action=&read_init_actions_fluid;
    add_bdry_types_fluid_to_codex(codex);
    add_init_types_fluid_to_codex(codex);
    SOAP_process_code(*argum, codex, SOAP_VARS_CLEAN_ADDED);
    codex->action=&read_init_actions;
    gl->nsinit=ns;
    for (spec=0; spec<ns; spec++) gl->initspecies[spec]=spec;

    /* this must not be performed if Init() is called within the cycle module or if reading data files to init the domain*/ 
    if (!gl->CONTROL_READ){
      #ifdef UNSTEADY
        for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          for (flux=0; flux<nf; flux++){
            (*np)[_ai(gl,i,j,k)].bs->Um1[flux]=(*np)[_ai(gl,i,j,k)].bs->U[flux];
            #if _RESTIME_BW > 2  
            (*np)[_ai(gl,i,j,k)].bs->Um2[flux]=(*np)[_ai(gl,i,j,k)].bs->Um1[flux];
            #endif
            #if _RESTIME_BW > 3 
            (*np)[_ai(gl,i,j,k)].bs->Um3[flux]=(*np)[_ai(gl,i,j,k)].bs->Um2[flux];
            #endif
          }
        }
      #endif

      #ifdef _RESTIME_STORAGE_TRAPEZOIDAL
//       resume_nodes_only_in_zone_and_update_bdry_nodes(*np, gl, _zone_intersection(gl->domain_all,_zone_expansion(gl->domain,+hbw_res_fluid)));
        for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          #ifdef _RESTIME_STORAGE_TRAPEZOIDAL_RESIDUAL
           for (flux=0; flux<nf; flux++){
            (*np)[_ai(gl,i,j,k)].bs->trapezoidalm1[flux]=0.0;
            (*np)[_ai(gl,i,j,k)].bs->trapezoidalm1_next[flux]=0.0;
           }
          #endif
          #ifdef _RESTIME_STORAGE_TRAPEZOIDAL_MUSCLVARS
           if (is_node_valid((*np)[_ai(gl,i,j,k)],TYPELEVEL_FLUID)) find_musclvars((*np)[_ai(gl,i,j,k)],gl,(*np)[_ai(gl,i,j,k)].bs->trapezoidalm1);
          #endif
        }
      #endif
    }
    codex->ACTIONPROCESSED=TRUE;
  }
#ifdef EMFIELD
  if (strcmp(actionname,_EMFIELD_ACTIONNAME)==0){
    if (!gl->INIT_FLUID_READ) SOAP_fatal_error(codex,"The fluid module %s must come before the emfield module %s within Init().",_FLUID_ACTIONNAME,_EMFIELD_ACTIONNAME);
    gl->INIT_EMFIELD_READ=TRUE;
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_EMFIELD_ACTIONNAME);
    codex->action=&read_init_actions_emfield;
    add_bdry_types_emfield_to_codex(codex);
    add_init_types_emfield_to_codex(codex);
    SOAP_process_code(*argum, codex, SOAP_VARS_CLEAN_ADDED);
    codex->action=&read_init_actions;
    #ifdef UNSTEADY
     /* this must not be performed if Init() is called within the cycle module or if reading data files to init the domain*/ 
     if (!gl->CONTROL_READ){
       for_ijk(gl->domain_lim,is,js,ks,ie,je,ke){
          for (flux=0; flux<nfe; flux++){
            (*np)[_ai(gl,i,j,k)].bs->Uemfieldm1[flux]=(*np)[_ai(gl,i,j,k)].bs->Uemfield[flux];
          }
       }
     }
    #endif
    codex->ACTIONPROCESSED=TRUE;
  }
#endif

}


void read_init(char *argum, SOAP_codex_t *codex){
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;

  codex->action=&read_init_actions;

  SOAP_process_code(argum, codex, SOAP_VARS_CLEAN_ADDED);
  if (!gl->CONTROL_READ){
      if (!gl->INIT_FLUID_READ) SOAP_fatal_error(codex,"The fluid module %s was not found within Init().",_FLUID_ACTIONNAME);
#ifdef EMFIELD  
      if (!gl->INIT_EMFIELD_READ) SOAP_fatal_error(codex,"The emfield module %s was not found within Init().",_EMFIELD_ACTIONNAME);
#endif
  }


}

