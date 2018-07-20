#include <src/bdry.h>
#include <model/metrics/_metrics.h>
#include <src/control.h>
#define minadjnode 2



double _bdry_param(np_t *np, gl_t *gl, long l, short param, int TYPELEVEL){
  double bdryparam;
  bdryparam=0.0; //needed to prevent compiler warning
#ifdef EMFIELD
  if (TYPELEVEL==TYPELEVEL_EMFIELD){
    if (param>np[l].numbdryparam_emf) fatal_error("Emfield boundary node (%ld,%ld,%ld) needs at least %d parameters specified. Specify those within Bdry() module using Param() action.", param+1,_i(l,gl,0),_i(l,gl,1),_i(l,gl,2));
    bdryparam=np[l].bdryparam_emf[param];
  }
#endif
  if (TYPELEVEL==TYPELEVEL_FLUID || TYPELEVEL==TYPELEVEL_FLUID_WORK){
    if (param>=np[l].numbdryparam) fatal_error("Fluid boundary node (%ld,%ld,%ld) needs at least %d parameter(s) specified. Specify those within Bdry() module using Param action.", _i(l,gl,0),_i(l,gl,1),_i(l,gl,2),param+1);
    bdryparam=np[l].bdryparam[param];
  }
  return(bdryparam);
}


void write_bdry_template(FILE **controlfile){
  wfprintf(*controlfile,"\n\n");
  wfprintf(*controlfile,
  "Bdry(\n");
  write_bdry_fluid_template(controlfile);
#ifdef EMFIELD
  write_bdry_emfield_template(controlfile);
#endif
  wfprintf(*controlfile,
  ");\n");
}



void read_bdry_actions(char *action, char **argum, SOAP_codex_t *codex){
  np_t **np;
  gl_t *gl;
  zone_t zone;
  long numparam,BCtype,i,j,k,is,ie,i1,i2,k1,k2,j1,j2,l1,l2;
#ifdef _2DL 
  long js,je;
#endif
#ifdef _3DL 
  long ks,ke;
#endif


  double *bdryparam;
  long BCtypes[nd*2];
  long TYPELEVEL,param;
  int eos=EOS;

  np=((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;
  TYPELEVEL=((readcontrolarg_t *)codex->action_args)->TYPELEVEL;


#ifdef EMFIELD
  if (strcmp(action,_EMFIELD_ACTIONNAME)==0){
    gl->BDRY_EMFIELD_READ=TRUE;
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE)  wfprintf(stdout,"%s..",_EMFIELD_ACTIONNAME);
    ((readcontrolarg_t *)codex->action_args)->TYPELEVEL=TYPELEVEL_EMFIELD;
    add_bdry_types_emfield_to_codex(codex);
    SOAP_process_code(*argum, codex, SOAP_VARS_CLEAN_ADDED);
    ((readcontrolarg_t *)codex->action_args)->TYPELEVEL=TYPELEVEL;
    codex->ACTIONPROCESSED=TRUE;
  }
#endif
  if (strcmp(action,_FLUID_ACTIONNAME)==0){
    gl->BDRY_FLUID_READ=TRUE;
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE)  wfprintf(stdout,"%s..",_FLUID_ACTIONNAME);
    add_bdry_types_fluid_to_codex(codex);
    SOAP_process_code(*argum, codex, SOAP_VARS_CLEAN_ADDED);
    codex->ACTIONPROCESSED=TRUE;
  }

  if (strcmp(action,"All")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (sscanf(*argum,"%ld%n",&BCtype,&eos)!=1 || (*argum)[eos]!=EOS)
      SOAP_fatal_error(codex,"One argument could not be read properly "
                     "in All() part of Bdry(). Argument: %s.",*argum);
    for1DL(i,gl->domain_lim.is,gl->domain_lim.ie)
      for2DL(j,gl->domain_lim.js,gl->domain_lim.je)
        for3DL(k,gl->domain_lim.ks,gl->domain_lim.ke)
          if (is_node_bdry((*np)[_ai(gl,i,j,k)],TYPELEVEL))
            update_node_type((*np),gl,TYPELEVEL,BCtype,_zone_from_point(i,j,k));
        end3DL
      end2DL
    end1DL
    codex->ACTIONPROCESSED=TRUE;
  }


  if (strcmp(action,"Region")==0) {
    SOAP_substitute_all_argums(argum, codex);
    #ifdef _2D
      if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld%n",&is,&js,&ie,&je,&BCtype,&eos)!=5
    #endif
    #ifdef _3D
      if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld,%ld%n",&is,&js,&ks,&ie,&je,&ke,&BCtype,&eos)!=7
    #endif
      || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"One or more argument(s) could not be read properly "
                             "in Region() part of Bdry(). Arguments: %s .",*argum);
    }
    find_zone_from_argum(*argum, 0, gl, codex, &zone);
    for1DL(i,zone.is,zone.ie)
      for2DL(j,zone.js,zone.je)
        for3DL(k,zone.ks,zone.ke)
          if (is_node_bdry((*np)[_ai(gl,i,j,k)],TYPELEVEL))
            update_node_type((*np),gl,TYPELEVEL,BCtype,_zone_from_point(i,j,k));
        end3DL
      end2DL
    end1DL
    codex->ACTIONPROCESSED=TRUE;
  }


  if (strcmp(action,"Param")==0) {
    SOAP_substitute_all_argums(argum, codex);
    find_zone_from_argum(*argum, 0, gl, codex, &zone);
    BCtype=SOAP_get_argum_long(codex,*argum,2*nd);
    numparam=SOAP_number_argums(*argum)-2*nd-1;
    bdryparam=(double *)malloc(numparam*sizeof(double));
    for (param=0; param<numparam; param++) bdryparam[param]=SOAP_get_argum_double(codex,*argum,2*nd+param+1);
    for1DL(i,zone.is,zone.ie)
      for2DL(j,zone.js,zone.je)
        for3DL(k,zone.ks,zone.ke)
          if (_node_type((*np)[_ai(gl,i,j,k)],TYPELEVEL)==BCtype){
            if ((*np)[_ai(gl,i,j,k)].bdryparam==NULL){
              (*np)[_ai(gl,i,j,k)].bdryparam=(double *)malloc(numparam*sizeof(double));
            } else {
              (*np)[_ai(gl,i,j,k)].bdryparam=(double *)realloc((*np)[_ai(gl,i,j,k)].bdryparam,numparam*sizeof(double));
            }
            for (param=0; param<numparam; param++){
              (*np)[_ai(gl,i,j,k)].bdryparam[param]=bdryparam[param];
            }
            (*np)[_ai(gl,i,j,k)].numbdryparam=numparam;
          }
        end3DL
      end2DL
    end1DL
    codex->ACTIONPROCESSED=TRUE;
  }



  if (strcmp(action,"Link")==0) {
    SOAP_substitute_all_argums(argum, codex);
    j1=1;
    j2=1;
    k1=1;
    k2=1;
#ifdef _2D
    if (sscanf(*argum,"%ld,%ld,%ld,%ld%n",&i1,&j1,&i2,&j2,&eos)!=4
#endif
#ifdef _3D
    if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld%n",&i1,&j1,&k1,&i2,&j2,&k2,&eos)!=6
#endif
      || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"One or more argument(s) could not be read properly "
                             "in Link() part of Bdry(). Arguments: %s .",*argum);
    }
    if (!(   (is_node_bdry((*np)[_ai(gl,i1,j1,k1)],TYPELEVEL) && is_node_inner((*np)[_ai(gl,i2,j2,k2)],TYPELEVEL))
        || (is_node_bdry((*np)[_ai(gl,i2,j2,k2)],TYPELEVEL) && is_node_inner((*np)[_ai(gl,i1,j1,k1)],TYPELEVEL))
        )){
      SOAP_fatal_error(codex,"When making a link between two nodes, one must be a bdry node and the other an inner node.");
    }

    switch (TYPELEVEL){   
      case TYPELEVEL_FLUID: 
        (*np)[_ai(gl,i1,j1,k1)].link=_ai(gl,i2,j2,k2);
        (*np)[_ai(gl,i2,j2,k2)].link=_ai(gl,i1,j1,k1);
#ifdef DISTMPI
        if (is_node_bdry((*np)[_ai(gl,i1,j1,k1)],TYPELEVEL)){
          (*np)[_ai(gl,i1,j1,k1)].linkmusclvars=(double *)realloc((*np)[_ai(gl,i1,j1,k1)].linkmusclvars,(max(0,hbw_resconv_fluid-1)*nf)*sizeof(double));
        } else {
          assert(is_node_bdry((*np)[_ai(gl,i2,j2,k2)],TYPELEVEL));
          (*np)[_ai(gl,i2,j2,k2)].linkmusclvars=(double *)realloc((*np)[_ai(gl,i2,j2,k2)].linkmusclvars,(max(0,hbw_resconv_fluid-1)*nf)*sizeof(double));
        }
#endif
      break;
#ifdef EMFIELD
      case TYPELEVEL_EMFIELD: 
        (*np)[_ai(gl,i1,j1,k1)].link_emf=_ai(gl,i2,j2,k2);
        (*np)[_ai(gl,i2,j2,k2)].link_emf=_ai(gl,i1,j1,k1);
      break;
#endif
      default:
         fatal_error("TYPELEVEL must be wither TYPELEVEL_FLUID or TYPELEVEL_EMFIELD in read_bdry_actions(), Link.");
    }
    codex->ACTIONPROCESSED=TRUE;
  }


  if (strcmp(action,"Unlink")==0) {
    SOAP_substitute_all_argums(argum, codex);
    #ifdef _2D
      if (sscanf(*argum,"%ld,%ld,%ld,%ld%n",&is,&js,&ie,&je,&eos)!=4
    #endif
    #ifdef _3D
      if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld%n",&is,&js,&ks,&ie,&je,&ke,&eos)!=6
    #endif
      || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"One or more argument(s) could not be read properly "
                             "in Unlink() part of Bdry(). Arguments: %s .",*argum);
    }
    find_zone_from_argum(*argum, 0, gl, codex, &zone);
    for1DL(i,zone.is,zone.ie)
      for2DL(j,zone.js,zone.je)
        for3DL(k,zone.ks,zone.ke)
          if (is_node_link((*np)[_ai(gl,i,j,k)],TYPELEVEL)){
            l1=_ai(gl,i,j,k);
            l2=_node_link((*np)[l1],TYPELEVEL);
            switch (TYPELEVEL){   
              case TYPELEVEL_FLUID: 
                (*np)[l1].link=LINK_NONE;
                (*np)[l2].link=LINK_NONE;
#ifdef DISTMPI
                free((*np)[l1].linkmusclvars);
                free((*np)[l2].linkmusclvars);
                (*np)[l1].linkmusclvars=NULL;
                (*np)[l2].linkmusclvars=NULL;
#endif
              break;
#ifdef EMFIELD
              case TYPELEVEL_EMFIELD: 
                (*np)[l1].link_emf=LINK_NONE;
                (*np)[l2].link_emf=LINK_NONE;
              break;
#endif
              default:
                fatal_error("TYPELEVEL must be wither TYPELEVEL_FLUID or TYPELEVEL_EMFIELD in read_bdry_actions(), Unlink.");
            }
          }
        end3DL
      end2DL
    end1DL
    codex->ACTIONPROCESSED=TRUE;
  }


  if (strcmp(action,"Cut")==0) {
    SOAP_substitute_all_argums(argum, codex);
    find_zone_from_argum(*argum, 0, gl, codex, &zone);
    if (          zone.ie-zone.is<2 
         if2DL(|| zone.je-zone.js<2)
         if3DL(|| zone.ke-zone.ks<2)
       ) SOAP_fatal_error(codex,"A Cut() command must create a cut with a bandwidth of at least 3 nodes along each dimension.");
    zone=_zone_intersection(zone,gl->domain_all);
    update_node_type(*np,gl,TYPELEVEL,NODETYPE_BDRY, zone);
    update_node_type(*np,gl,TYPELEVEL,NODETYPE_UNUSED, _zone_expansion(zone,-1));
    codex->ACTIONPROCESSED=TRUE;
  }

  if (strcmp(action,"Faces")==0) {
    SOAP_substitute_all_argums(argum, codex);
#ifdef _1D
    if (sscanf(*argum,"%ld,%ld%n",&BCtypes[0],&BCtypes[1],&eos)!=2
#endif
#ifdef _2D
    if (sscanf(*argum,"%ld,%ld,%ld,%ld%n",&BCtypes[0],&BCtypes[1],&BCtypes[2],&BCtypes[3],&eos)!=4 
#endif
#ifdef _3D
    if (sscanf(*argum,"%ld,%ld,%ld,%ld,%ld,%ld%n",&BCtypes[0],&BCtypes[1],&BCtypes[2],&BCtypes[3],&BCtypes[4],&BCtypes[5],&eos)!=6
#endif
      || (*argum)[eos]!=EOS){
      SOAP_fatal_error(codex,"One or more argument(s) could not be read properly "
                             "in Faces() part of Bdry(). Arguments: %s .",*argum);
    }

    zone=gl->domain_all;
    zone.ie=zone.is;
    update_node_type(*np,gl,TYPELEVEL,BCtypes[0],zone);
    zone=gl->domain_all;
    zone.is=zone.ie;
    update_node_type(*np,gl,TYPELEVEL,BCtypes[1], zone);
#ifdef _2DL
    zone=gl->domain_all;
    zone.je=zone.js;
    update_node_type(*np,gl,TYPELEVEL,BCtypes[2], zone);
    zone=gl->domain_all;
    zone.js=zone.je;
    update_node_type(*np,gl,TYPELEVEL,BCtypes[3], zone);
#endif
#ifdef _3DL
    zone=gl->domain_all;
    zone.ke=zone.ks;
    update_node_type(*np,gl,TYPELEVEL,BCtypes[4], zone);
    zone=gl->domain_all;
    zone.ks=zone.ke;
    update_node_type(*np,gl,TYPELEVEL,BCtypes[5], zone);
#endif
    codex->ACTIONPROCESSED=TRUE;
  }
}




/* returns TRUE if found, FALSE if not found */
bool find_bdry_direc(np_t *np, gl_t *gl, long l_A, int TYPELEVEL,
                   long *theta, long *thetasgn){
  long dim;
  bool FOUND;

  FOUND=FALSE;
  *theta=0;
  *thetasgn=+1;
//  for (dim=nd-1; dim>=0; dim--){
  for (dim=0; dim<nd; dim++){
    if (is_node_inner(np[_al(gl,l_A,dim,+1)],TYPELEVEL) && is_node_inner(np[_al(gl,l_A,dim,+2)],TYPELEVEL)) {
       if (!FOUND){
         *theta=dim;
         *thetasgn=+1;
         FOUND=TRUE;
       } else {
         //only substitute previous find if distance between inner and bdry node is less (improves convergence rates)
         if (distance2_between_nodes(np,gl,l_A,_al(gl,l_A,dim,+1))<1.001*distance2_between_nodes(np,gl,l_A,_al(gl,l_A,*theta,*thetasgn))){
           *theta=dim;
           *thetasgn=+1;
         }
       }
    }
    if (is_node_inner(np[_al(gl,l_A,dim,-1)],TYPELEVEL) && is_node_inner(np[_al(gl,l_A,dim,-2)],TYPELEVEL)) {
       if (!FOUND){
         *theta=dim;
         *thetasgn=-1;
         FOUND=TRUE;
       } else {
         //only substitute previous find if distance between inner and bdry node is less (improves convergence rates)
         if (distance2_between_nodes(np,gl,l_A,_al(gl,l_A,dim,-1))<1.001*distance2_between_nodes(np,gl,l_A,_al(gl,l_A,*theta,*thetasgn))){
           *theta=dim;
           *thetasgn=-1;
         }

       }
    }
  }
  return(FOUND);
}


void find_link_direc(np_t *np, gl_t *gl, long llink, int TYPELEVEL, long *theta, long *thetasgn){
  long dim,dimsgn;
  if (llink==LINK_NONE) fatal_error("The node given must link to another in _al_link().");
  if (!is_node_inner(np[llink],TYPELEVEL)) fatal_error("The node linked to must be an inner node in _al_link().");
  /* check on which side there is a node that links to another */
  *theta=-1;
  *thetasgn=0;
  for (dim=0; dim<nd; dim++){
    for (dimsgn=-1; dimsgn<=1; dimsgn+=2){
      if (is_node_link(np[_al(gl,llink,dim,dimsgn)],TYPELEVEL) && is_node_bdry(np[_al(gl,llink,dim,dimsgn)],TYPELEVEL)){
        if (*theta!=-1) fatal_error("link node (%ld,%ld,%ld) has more than 1 neighbor node that is both a link node and a bdry node; can't determine the direction of the surface properly.",_i(llink, gl, 0),_i(llink, gl, 1),_i(llink, gl, 2));
        *theta=dim;
        *thetasgn=-dimsgn;
      }
    }
  }
  if (*theta==-1) fatal_error("link node (%ld,%ld,%ld) has no neighbor node that is both a link node and a bdry node; can't determine the direction of the surface properly.",_i(llink, gl, 0),_i(llink, gl, 1),_i(llink, gl, 2));
}


/* returns TRUE if one boundary node has multiple inner nodes adjacent to it,
   and FALSE otherwise */
bool is_multiple_bdry_direc(np_t *np, gl_t *gl, long l_A, int TYPELEVEL){
  long dim,cnt;
  bool MULTIPLEBDRY;

  cnt=0;
  for (dim=0; dim<nd; dim++){
    if (is_node_inner(np[_al(gl,l_A,dim,+1)],TYPELEVEL) && is_node_inner(np[_al(gl,l_A,dim,+2)],TYPELEVEL)) cnt++;
    if (is_node_inner(np[_al(gl,l_A,dim,-1)],TYPELEVEL) && is_node_inner(np[_al(gl,l_A,dim,-2)],TYPELEVEL)) cnt++;
  }
  if (cnt>1) MULTIPLEBDRY=TRUE; else MULTIPLEBDRY=FALSE;
  return(MULTIPLEBDRY);
}



void update_node_type(np_t *np, gl_t *gl, int TYPELEVEL, long type, zone_t zone){
  long i,j,k;
  zone=_zone_intersection(zone,gl->domain_lim);
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        if (TYPELEVEL==TYPELEVEL_FLUID) np[_ai(gl,i,j,k)].type=type;
        if (TYPELEVEL==TYPELEVEL_FLUID_WORK) np[_ai(gl,i,j,k)].type_wk=type;
#ifdef EMFIELD
        if (TYPELEVEL==TYPELEVEL_EMFIELD) np[_ai(gl,i,j,k)].type_emf=type;
#endif

      end3DL
    end2DL
  end1DL
}


void copy_base_to_work_node_type(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        np[_ai(gl,i,j,k)].type_wk=np[_ai(gl,i,j,k)].type;
      end3DL
    end2DL
  end1DL
}


void update_inner_to_bdry_node_type(np_t *np, gl_t *gl, int TYPELEVEL, long type, zone_t zone){
  long i,j,k;
  zone=_zone_intersection(zone,gl->domain_lim);
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL)) {
          if (TYPELEVEL==TYPELEVEL_FLUID) np[_ai(gl,i,j,k)].type=type;
          if (TYPELEVEL==TYPELEVEL_FLUID_WORK) np[_ai(gl,i,j,k)].type_wk=type;
#ifdef EMFIELD
          if (TYPELEVEL==TYPELEVEL_EMFIELD) np[_ai(gl,i,j,k)].type_emf=type;
#endif
        }
      end3DL
    end2DL
  end1DL
}


void update_inner_and_bdry_to_node_type(np_t *np, gl_t *gl, int TYPELEVEL, long type, zone_t zone){
  long i,j,k;
  zone=_zone_intersection(zone,gl->domain_lim);
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL) || is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL)) {
          if (TYPELEVEL==TYPELEVEL_FLUID) np[_ai(gl,i,j,k)].type=type;
          if (TYPELEVEL==TYPELEVEL_FLUID_WORK) np[_ai(gl,i,j,k)].type_wk=type;
#ifdef EMFIELD
          if (TYPELEVEL==TYPELEVEL_EMFIELD) np[_ai(gl,i,j,k)].type_emf=type;
#endif
        }
      end3DL
    end2DL
  end1DL
}




/* ???
   note that a BC type of NODETYPE_BDRYOUTFLOW is given to any node
   which is surrounded by boundary conditions. This assumes
   that hbw_bdry_fluid=2. For hbw_bdry_fluid having another value, this subroutine
   must be rewritten */
#if (hbw_bdry_fluid<2) 
 // not sure if adjust_node_type can be used with hbw_bdry_fluid>2: this needs further testing
 #error function adjust_node_type() in src/bdry.c can only be used with hbw_bdry_fluid>=2
#endif

void adjust_node_type(np_t *np, gl_t *gl, zone_t zone, int TYPELEVEL){
  long i,j,k,dim,dim2,offset,offset2,offset3,l;
  bool FLAG;

  /* limit the zone to gl->domain_lim minus 1 ????? */
  zone=_zone_intersection(zone,_zone_expansion(gl->domain_lim,-1));
  /* transform inner node to boundary node if along
     any direction, the segment is made up of only 1 node */
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        FLAG=FALSE;
        for (dim=0; dim<nd; dim++){
          l=_ai(gl,i,j,k);
          if (is_node_inner(np[l],TYPELEVEL) && is_node_bdry(np[_al(gl,l,dim,+1)],TYPELEVEL)
              && is_node_bdry(np[_al(gl,l,dim,-1)],TYPELEVEL) ){
              FLAG=TRUE;
          }
        }
        if (FLAG) update_node_type(np,gl,TYPELEVEL,NODETYPE_BDRYOUTFLOW, _zone_from_point(i,j,k));

      end3DL
    end2DL
  end1DL


  /* set to UNUSED all boundary nodes that have no inner nodes as a neighbor */
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        FLAG=FALSE;
          for (offset=-1; offset<=1; offset++){
            for (offset2=-1; offset2<=1; offset2++){
              for (offset3=-1; offset3<=1; offset3++){
                if (is_node_inner(np[_ai(gl,i+offset,j+offset2,k+offset3)],TYPELEVEL)) {
                  FLAG=TRUE;
                }
              }
            }
          }
        if ((!FLAG) && (is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL))) {
          update_node_type(np,gl,TYPELEVEL,NODETYPE_UNUSED, _zone_from_point(i,j,k));
        }
      end3DL
    end2DL
  end1DL


  /* update BC nodes located in corners (in 3D) to outflow 
     --> don't do this anymore. This leads to problems with cross-difference terms within FDS
  */
  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        FLAG=FALSE;
        for (dim=0; dim<nd; dim++){
          for (dim2=0; dim2<nd; dim2++){
            for (offset=-1; offset<=1; offset++){
              for (offset2=-1; offset2<=1; offset2++){
                if (dim!=dim2) {
                  if (is_node_inner(np[_all(gl,_ai(gl,i,j,k),dim,offset,dim2,offset2)],TYPELEVEL)) {
                    FLAG=TRUE;
                  }
                }
              }
            }
          }
        }
        if ((!FLAG) && (is_node_bdry(np[_ai(gl,i,j,k)],TYPELEVEL))) {
          //update_node_type(np,gl,TYPELEVEL,NODETYPE_UNUSED, _zone_from_point(i,j,k));
          //update_node_type(np,gl,TYPELEVEL,NODETYPE_BDRYOUTFLOW, _zone_from_point(i,j,k));
        }
      end3DL
    end2DL
  end1DL
}


/* if (nodetype>9)
                wfprintf(outfile,"%c",(int)('A'+nodetype-10));
                else wfprintf(outfile,"%ld",nodetype);
*/
char _bdry_ID(long nodetype){
  char ID[10];
  if (nodetype>9)
    sprintf(ID,"%c",(int)('a'+nodetype-10));
  else sprintf(ID,"%ld",nodetype);
  return(ID[0]);
}


void display_node_type_window(FILE *outfile, np_t *np, gl_t *gl, int TYPELEVEL,
                        long im, long jm, long km, long bw){
  long i,j,k,is,ie;
#ifdef _2DL
  long js,je;
#endif
#ifdef _3DL
  long ks,ke;
#endif
  long nodetype,link;
  int NODERESUMED;
  bool PRINT=TRUE;
#ifdef DISTMPI
  int rank;
  MPI_Status MPI_Status1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank!=0) PRINT=FALSE;
#endif

  im=min(gl->domain_all.ie,max(gl->domain_all.is,im));
  jm=min(gl->domain_all.je,max(gl->domain_all.js,jm));
  km=min(gl->domain_all.ke,max(gl->domain_all.ks,km));

#ifdef _3DL
  ke=max(gl->domain_all.ks-1,min(km+bw,gl->domain_all.ke+1));
  ks=min(gl->domain_all.ke+1,max(km-bw,gl->domain_all.ks-1));
#endif

#ifdef _2DL
  je=max(gl->domain_all.js-1,min(jm+bw,gl->domain_all.je+1));
  js=min(gl->domain_all.je+1,max(jm-bw,gl->domain_all.js-1));
#endif

  ie=max(gl->domain_all.is-1,min(im+bw,gl->domain_all.ie+1));
  is=min(gl->domain_all.ie+1,max(im-bw,gl->domain_all.is-1));


  i=im;  j=jm;  k=km;

  if (PRINT) wfprintf(outfile,"\n\n");
#ifdef _3D
  if (PRINT) wfprintf(outfile,"Cut along k=%ld \n",km);
#endif
  k=km;
#ifdef _2DL
  for (j=je; j>=js; j--){
#endif
    for (i=is; i<=ie; i++){
#ifdef DISTMPI
      if (_node_rank(gl,i,j,k)==rank) {
        NODERESUMED=(int)is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID);
        nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
        link=_node_link(np[_ai(gl,i,j,k)],TYPELEVEL);
        if (rank!=0) {
          MPI_IBsend(&NODERESUMED,1,MPI_INT,0,0,MPI_COMM_WORLD);
          MPI_IBsend(&nodetype,1,MPI_LONG,0,0,MPI_COMM_WORLD);
          MPI_IBsend(&link,1,MPI_LONG,0,0,MPI_COMM_WORLD);
        }
      }
      if (rank==0 && _node_rank(gl,i,j,k)!=0) {
        MPI_Recv(&NODERESUMED,1,MPI_INT,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
        MPI_Recv(&nodetype,1,MPI_LONG,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
        MPI_Recv(&link,1,MPI_LONG,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
      }
#else
      nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
      link=_node_link(np[_ai(gl,i,j,k)],TYPELEVEL);
      NODERESUMED=is_node_resumed(np[_ai(gl,i,j,k)]);
#endif
      if (PRINT){
        if (NODERESUMED || TYPELEVEL==TYPELEVEL_FLUID
#ifdef EMFIELD
        || TYPELEVEL==TYPELEVEL_EMFIELD
#endif
        ) {
          if (link!=LINK_NONE){
            wfprintf(outfile,"L");
          } else {
            if (nodetype>=NODETYPE_BDRY) {
              wfprintf(outfile,"%c",_bdry_ID(nodetype));
            } else {
              if (nodetype==NODETYPE_INNER) wfprintf(outfile,"+"); /* Inner flow node */
              if (nodetype==NODETYPE_UNUSED) wfprintf(outfile,"%c",CHARACTER_BLOCK); /* unused node */
            }
          }
        } else {
          wfprintf(outfile,"S");
        }
      }
    } /* for1DL */
#ifdef _2DL
    if (PRINT && j!=jm) wfprintf(outfile,"\n"); else wfprintf(outfile,"<j=%ld\n",jm);
  } /* for2DL */
#else
  if (PRINT) wfprintf(outfile,"\n");
#endif
  if (PRINT){
    for (i=is; i<im; i++) wfprintf(outfile," ");
    wfprintf(outfile,"^\n");
    for (i=is; i<im-1; i++) wfprintf(outfile," ");
    wfprintf(outfile,"i=%ld\n",im);
  }



#ifdef _3D
  if (PRINT) wfprintf(outfile,"\nCut along i=%ld \n",im);
  i=im;
  for (j=je; j>=js; j--){
    for (k=ks; k<=ke; k++){
#ifdef DISTMPI
      if (_node_rank(gl,i,j,k)==rank) {
        NODERESUMED=(int)is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID);
        nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
        if (rank!=0) {
          MPI_IBsend(&NODERESUMED,1,MPI_INT,0,0,MPI_COMM_WORLD);
          MPI_IBsend(&nodetype,1,MPI_LONG,0,0,MPI_COMM_WORLD);
        }
      }
      if (rank==0 && _node_rank(gl,i,j,k)!=0) {
        MPI_Recv(&NODERESUMED,1,MPI_INT,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
        MPI_Recv(&nodetype,1,MPI_LONG,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
      }
#else
      nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
      link=_node_link(np[_ai(gl,i,j,k)],TYPELEVEL);
      NODERESUMED=is_node_resumed(np[_ai(gl,i,j,k)]);
#endif
      if (PRINT) {
        if (NODERESUMED || TYPELEVEL==TYPELEVEL_FLUID
#ifdef EMFIELD
        || TYPELEVEL==TYPELEVEL_EMFIELD
#endif
           ) {
          if (link!=LINK_NONE){
            wfprintf(outfile,"L");
          } else {
            if (nodetype>=NODETYPE_BDRY) {
              wfprintf(outfile,"%c",_bdry_ID(nodetype));
            } else {
              if (nodetype==NODETYPE_INNER) wfprintf(outfile,"+"); /* Inner flow node */
              if (nodetype==NODETYPE_UNUSED) wfprintf(outfile,"%c",CHARACTER_BLOCK); /* unused node */
            }
          }
        } else {
          wfprintf(outfile,"S");
        }
      }
    }
    if (PRINT && j!=jm) wfprintf(outfile,"\n"); else wfprintf(outfile,"<j=%ld\n",jm);
  }
  if (PRINT){
    for (k=ks; k<km; k++) wfprintf(outfile," ");
    wfprintf(outfile,"^\n");
    for (k=ks; k<km-1; k++) wfprintf(outfile," ");
    wfprintf(outfile,"k=%ld\n",km);

    wfprintf(outfile,"\nCut along j=%ld \n",jm);
  }
  j=jm;
  for (i=ie; i>=is; i--){
    for (k=ks; k<=ke; k++){
#ifdef DISTMPI
      if (_node_rank(gl,i,j,k)==rank) {
        NODERESUMED=(int)is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID);
        nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
        if (rank!=0) {
          MPI_IBsend(&NODERESUMED,1,MPI_INT,0,0,MPI_COMM_WORLD);
          MPI_IBsend(&nodetype,1,MPI_LONG,0,0,MPI_COMM_WORLD);
        }
      }
      if (rank==0 && _node_rank(gl,i,j,k)!=0) {
        MPI_Recv(&NODERESUMED,1,MPI_INT,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
        MPI_Recv(&nodetype,1,MPI_LONG,_node_rank(gl,i,j,k),0,MPI_COMM_WORLD,&MPI_Status1);
      }
#else
      nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
      NODERESUMED=is_node_resumed(np[_ai(gl,i,j,k)]);
#endif
      if (PRINT) {
        if (NODERESUMED || TYPELEVEL==TYPELEVEL_FLUID
#ifdef EMFIELD
        || TYPELEVEL==TYPELEVEL_EMFIELD
#endif
           ) {
          if (nodetype>=NODETYPE_BDRY) {
            if (nodetype>9)
              wfprintf(outfile,"%c",(int)('A'+nodetype-10));
              else wfprintf(outfile,"%ld",nodetype);
          } else {
            if (nodetype==NODETYPE_INNER) wfprintf(outfile,"+"); /* Inner flow node */
            if (nodetype==NODETYPE_UNUSED) wfprintf(outfile,"%c",CHARACTER_BLOCK); /* unused node */
          }
        } else {
          wfprintf(outfile,"S");
        }
      }
    }
    if (PRINT && i!=im) wfprintf(outfile,"\n"); else wfprintf(outfile,"<i=%ld\n",im);
  }
  if (PRINT){
    for (k=ks; k<km; k++) wfprintf(outfile," ");
    wfprintf(outfile,"^\n");
    for (k=ks; k<km-1; k++) wfprintf(outfile," ");
    wfprintf(outfile,"k=%ld\n",km);
  }
#endif
  if (PRINT) wfprintf(outfile,"\nNote that boundary nodes numbered 10, 11, ...,  will be shown as %c, %c, ...\n\n",_bdry_ID(10),_bdry_ID(11));
#ifdef DISTMPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}


void display_node_type(np_t *np, gl_t *gl, zone_t zone, int TYPELEVEL){

  display_node_type_window(stdout, np, gl, TYPELEVEL, (zone.ie+zone.is)/2,
                     (zone.je+zone.js)/2,(zone.ke+zone.ks)/2, max((zone.ke-zone.ks)/2,max((zone.ie-zone.is)/2,(zone.je-zone.js)/2)));
}


void display_zone(FILE *outfile, zone_t zone){
  fprintf(outfile,"%ld..%ld, %ld..%ld, %ld..%ld",zone.is,zone.ie,zone.js, zone.je, zone.ks, zone.ke);
}


void display_node_type_window_local_process(FILE *outfile, np_t *np, gl_t *gl, int TYPELEVEL,
                        long im, long jm, long km, long bw){
  long i,j,k,is,ie;
#ifdef _2DL
  long js,je;
#endif
#ifdef _3DL
  long ks,ke;
#endif
  long nodetype;
  int NODERESUMED;
  bool PRINT=TRUE;
  zone_t domain_eff;

  domain_eff=_zone_intersection(gl->domain_all,gl->domain_lim);
  if (PRINT) {
    fprintf(outfile,"\n\ndomain_lim: ");
    display_zone(outfile,gl->domain_lim);
    fprintf(outfile,"\ndomain:     ");
    display_zone(outfile,gl->domain);
    fprintf(outfile,"\ndomain_all: ");
    display_zone(outfile,gl->domain_all);
  }

  
  im=min(domain_eff.ie,max(domain_eff.is,im));
  jm=min(domain_eff.je,max(domain_eff.js,jm));
  km=min(domain_eff.ke,max(domain_eff.ks,km));

#ifdef _3DL
  ke=max(domain_eff.ks-1,min(km+bw,domain_eff.ke+1));
  ks=min(domain_eff.ke+1,max(km-bw,domain_eff.ks-1));
#endif

#ifdef _2DL
  je=max(domain_eff.js-1,min(jm+bw,domain_eff.je+1));
  js=min(domain_eff.je+1,max(jm-bw,domain_eff.js-1));
#endif

  ie=max(domain_eff.is-1,min(im+bw,domain_eff.ie+1));
  is=min(domain_eff.ie+1,max(im-bw,domain_eff.is-1));


  i=im;  j=jm;  k=km;

  if (PRINT) fprintf(outfile,"\n\n");
#ifdef _3D
  if (PRINT) fprintf(outfile,"Cut along k=%ld \n",km);
#endif
  k=km;
#ifdef _2DL
  for (j=je; j>=js; j--){
#endif
    for (i=is; i<=ie; i++){
      nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
      NODERESUMED=is_node_resumed(np[_ai(gl,i,j,k)]);
      if (PRINT){
        if (NODERESUMED || TYPELEVEL==TYPELEVEL_FLUID
#ifdef EMFIELD
        || TYPELEVEL==TYPELEVEL_EMFIELD
#endif
        ) {
          if (nodetype>=NODETYPE_BDRY) {
            if (nodetype>9)
              fprintf(outfile,"%c",(int)('A'+nodetype-10));
              else fprintf(outfile,"%ld",nodetype);
          } else {
            if (nodetype==NODETYPE_INNER) fprintf(outfile,"+"); /* Inner flow node */
            if (nodetype==NODETYPE_UNUSED) fprintf(outfile,"%c",CHARACTER_BLOCK); /* unused node */
          }
        } else {
          fprintf(outfile,"S");
        }
      }
    } /* for1DL */
#ifdef _2DL
    if (PRINT && j!=jm) fprintf(outfile,"\n"); else fprintf(outfile,"<j=%ld\n",jm);
  } /* for2DL */
#else
  if (PRINT) fprintf(outfile,"\n");
#endif
  if (PRINT){
    for (i=is; i<im; i++) fprintf(outfile," ");
    fprintf(outfile,"^\n");
    for (i=is; i<im-1; i++) fprintf(outfile," ");
    fprintf(outfile,"i=%ld\n",im);
  }



#ifdef _3D
  if (PRINT) fprintf(outfile,"\nCut along i=%ld \n",im);
  i=im;
  for (j=je; j>=js; j--){
    for (k=ks; k<=ke; k++){
      nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
      NODERESUMED=is_node_resumed(np[_ai(gl,i,j,k)]);
      if (PRINT) {
        if (NODERESUMED || TYPELEVEL==TYPELEVEL_FLUID) {
          if (nodetype>=NODETYPE_BDRY) {
            if (nodetype>9) fprintf(outfile,"%c",(int)('A'+nodetype-10));
            else fprintf(outfile,"%ld",nodetype);
          } else {
            if (nodetype==NODETYPE_INNER) fprintf(outfile,"+"); /* Inner flow node */
            if (nodetype==NODETYPE_UNUSED) fprintf(outfile,"%c",CHARACTER_BLOCK); /* unused node */
          }
        } else {
          fprintf(outfile,"S");
        }
      }
    }
    if (PRINT && j!=jm) fprintf(outfile,"\n"); else fprintf(outfile,"<j=%ld\n",jm);
  }
  if (PRINT){
    for (k=ks; k<km; k++) fprintf(outfile," ");
    fprintf(outfile,"^\n");
    for (k=ks; k<km-1; k++) fprintf(outfile," ");
    fprintf(outfile,"k=%ld\n",km);

    fprintf(outfile,"\nCut along j=%ld \n",jm);
  }
  j=jm;
  for (i=ie; i>=is; i--){
    for (k=ks; k<=ke; k++){
      nodetype=_node_type(np[_ai(gl,i,j,k)],TYPELEVEL);
      NODERESUMED=is_node_resumed(np[_ai(gl,i,j,k)]);
      if (PRINT) {
        if (NODERESUMED || TYPELEVEL==TYPELEVEL_FLUID) {
          if (nodetype>=NODETYPE_BDRY) {
            if (nodetype>9)
              fprintf(outfile,"%c",(int)('A'+nodetype-10));
              else fprintf(outfile,"%ld",nodetype);
          } else {
            if (nodetype==NODETYPE_INNER) fprintf(outfile,"+"); /* Inner flow node */
            if (nodetype==NODETYPE_UNUSED) fprintf(outfile,"%c",CHARACTER_BLOCK); /* unused node */
          }
        } else {
          fprintf(outfile,"S");
        }
      }
    }
    if (PRINT && i!=im) fprintf(outfile,"\n"); else fprintf(outfile,"<i=%ld\n",im);
  }
  if (PRINT){
    for (k=ks; k<km; k++) fprintf(outfile," ");
    fprintf(outfile,"^\n");
    for (k=ks; k<km-1; k++) fprintf(outfile," ");
    fprintf(outfile,"k=%ld\n",km);
  }
#endif
  if (PRINT) fprintf(outfile,"\nNote that boundary nodes numbered 10, 11, ...,  will be shown as A, B, ...\n\n");
}




/* this will check nodes in zone and determine
   if there are at least minadjnode adjacent nodes
   everywhere, in every dimension */
static void validate_inner_nodes(np_t *np, gl_t *gl, zone_t zone, int TYPELEVEL){
  long i,j,k,dim,l,noderem,ls;

  for1DL(i,zone.is,zone.ie)
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        ls=_ai(gl,i,j,k);
        if (is_node_inner(np[ls],TYPELEVEL)){
          for (dim=0; dim<nd; dim++){
            noderem=minadjnode-1;
            l=ls;
            do {
              l=_al(gl,l,dim,+1);
              if (is_node_inner(np[l],TYPELEVEL)) noderem=noderem-1;
            } while ((noderem>0) && (is_node_inner(np[l],TYPELEVEL)));
            l=ls;
            do {
              l=_al(gl,l,dim,-1);
              if (is_node_inner(np[l],TYPELEVEL)) noderem=noderem-1;
            } while ((noderem>0) && (is_node_inner(np[l],TYPELEVEL)));
            if (noderem>0) {
              display_node_type_window(stderr, np, gl, TYPELEVEL, i, j, k, 20);
              fatal_error("Not enough adjacent inner nodes.");
            }
          }
        }
      end3DL
    end2DL
  end1DL
}




void read_bdry(char *argum, SOAP_codex_t *codex){
  np_t *np;
  gl_t *gl;
  long i,j,k;

  gl=((readcontrolarg_t *)codex->action_args)->gl;

#ifdef DISTMPI
  np_t **np_old;
  gl_t gl_all;
  zone_t window_old,domain_lim_old,domain_all_old,domain_old;
  domain_lim_old=gl->domain_lim;
  domain_all_old=gl->domain_all;
  domain_old=gl->domain;
  window_old=gl->window;
  gl_all=*gl;
  init_data_structure(&np, &gl_all, gl->domain_all, gl->domain_all);
  gl->domain_lim=gl_all.domain_lim;
  gl->domain=gl_all.domain;
  gl->domain_all=gl_all.domain_all;
  gl->domain_lim_all=gl_all.domain_lim;
  gl->window=gl_all.window;
  np_old=((readcontrolarg_t *)codex->action_args)->np;
  ((readcontrolarg_t *)codex->action_args)->np=&np;
#else
  np=*(((readcontrolarg_t *)codex->action_args)->np);
  gl->domain_lim_all=gl->domain_lim;
#endif


  gl->BDRY_FLUID_READ=FALSE;
  gl->BDRY_EMFIELD_READ=FALSE;

  for1DL(i,gl->domain_lim.is,gl->domain_lim.ie)
    for2DL(j,gl->domain_lim.js,gl->domain_lim.je)
      for3DL(k,gl->domain_lim.ks,gl->domain_lim.ke)
        np[_ai(gl,i,j,k)].link=LINK_NONE;
        np[_ai(gl,i,j,k)].numbdryparam=0;      
        np[_ai(gl,i,j,k)].bdryparam=NULL;      
#ifdef EMFIELD
        np[_ai(gl,i,j,k)].link_emf=LINK_NONE;
        np[_ai(gl,i,j,k)].numbdryparam_emf=0;      
        np[_ai(gl,i,j,k)].bdryparam_emf=NULL;      
#endif
#ifdef DISTMPI
        np[_ai(gl,i,j,k)].numlinkmusclvars=0;
        np[_ai(gl,i,j,k)].linkmusclvars=NULL;
#endif
      end3DL
    end2DL
  end1DL

  update_node_type(np,gl,TYPELEVEL_FLUID,NODETYPE_UNUSED, gl->domain_lim);
  update_node_type(np,gl,TYPELEVEL_FLUID,NODETYPE_BDRY, gl->domain_all);
  update_node_type(np,gl,TYPELEVEL_FLUID,NODETYPE_INNER,_zone_expansion(gl->domain_all,-1));
#ifdef EMFIELD
  update_node_type(np,gl,TYPELEVEL_EMFIELD,NODETYPE_UNUSED, gl->domain_lim);
  update_node_type(np,gl,TYPELEVEL_EMFIELD,NODETYPE_BDRY, gl->domain_all);
  update_node_type(np,gl,TYPELEVEL_EMFIELD,NODETYPE_INNER,_zone_expansion(gl->domain_all,-1));
#endif
  ((readcontrolarg_t *)codex->action_args)->TYPELEVEL=TYPELEVEL_FLUID;
  codex->ACTION=TRUE;
  codex->action=&read_bdry_actions;
  SOAP_process_code(argum, codex, SOAP_VARS_KEEP_ALL);

  if (!gl->BDRY_FLUID_READ) fatal_error("The fluid module %s was not found within Bdry().",_FLUID_ACTIONNAME);
  adjust_node_type(np, gl, gl->domain_lim, TYPELEVEL_FLUID);
  validate_inner_nodes(np, gl, gl->domain, TYPELEVEL_FLUID);
#ifdef EMFIELD
  if (!gl->BDRY_EMFIELD_READ) fatal_error("The emfield module %s was not found within Bdry().",_EMFIELD_ACTIONNAME);
  adjust_node_type(np, gl, gl->domain_lim, TYPELEVEL_EMFIELD);
  validate_inner_nodes(np, gl, gl->domain, TYPELEVEL_EMFIELD);
#endif



#ifdef DISTMPI
  int proc;
  zone_t zone;
  gl->domain_lim=domain_lim_old;
  gl->domain=domain_old;
  gl->domain_all=domain_all_old;
  gl->window=window_old;
  ((readcontrolarg_t *)codex->action_args)->np=np_old;


//  printf("(%ld %ld  %ld %ld)",gl->domain_lim.is,gl->domain_lim.js,gl->domain_lim.ie,gl->domain_lim.je);
  update_node_type(*(((readcontrolarg_t *)codex->action_args)->np),gl,TYPELEVEL_FLUID,NODETYPE_UNUSED, gl->domain_lim);
#ifdef EMFIELD
  update_node_type(*(((readcontrolarg_t *)codex->action_args)->np),gl,TYPELEVEL_EMFIELD,NODETYPE_UNUSED, gl->domain_lim);
#endif
 /* update_node_type(*(((readcontrolarg_t *)codex->action_args)->np),gl,TYPELEVEL_FLUID,NODETYPE_BDRY, _zone_intersection(gl->domain_all,_zone_expansion(gl->domain_lim,-1)));
*/

  for1DL (i,gl->domain_lim.is,gl->domain_lim.ie)
    for2DL (j,gl->domain_lim.js,gl->domain_lim.je)
      for3DL (k,gl->domain_lim.ks,gl->domain_lim.ke)
//        if (is_node_in_zone(i,j,k,gl->domain) || is_node_in_zone(i,j,k,_zone_expansion(gl->domain_lim,-1))){
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].type=np[_ai(&gl_all,i,j,k)].type;
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].type_wk=np[_ai(&gl_all,i,j,k)].type_wk;

        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].link=np[_ai(&gl_all,i,j,k)].link;
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].numbdryparam=np[_ai(&gl_all,i,j,k)].numbdryparam;      
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].bdryparam=np[_ai(&gl_all,i,j,k)].bdryparam;      


#ifdef EMFIELD
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].type_emf=np[_ai(&gl_all,i,j,k)].type_emf;

        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].link_emf=np[_ai(&gl_all,i,j,k)].link_emf;
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].numbdryparam_emf=np[_ai(&gl_all,i,j,k)].numbdryparam_emf;      
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].bdryparam_emf=np[_ai(&gl_all,i,j,k)].bdryparam_emf;      

#endif


#ifdef DISTMPI

        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].numlinkmusclvars=np[_ai(&gl_all,i,j,k)].numlinkmusclvars;
        (*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].linkmusclvars=np[_ai(&gl_all,i,j,k)].linkmusclvars;
        if (is_node_bdry((*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)],TYPELEVEL_FLUID)
         && is_node_link((*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)],TYPELEVEL_FLUID)) {
          assert((*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].linkmusclvars!=NULL); 
        } else {
          assert((*(((readcontrolarg_t *)codex->action_args)->np))[_ai(gl,i,j,k)].linkmusclvars==NULL); 
        }
#endif

//        }
      end3DL
    end2DL
  end1DL

  free(np);
//  wfprintf(stdout,"[MPI check of node types]");
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  for1DL (i,gl->domain_lim.is,gl->domain_lim.ie)
    for2DL (j,gl->domain_lim.js,gl->domain_lim.je)
      for3DL (k,gl->domain_lim.ks,gl->domain_lim.ke)
        if (!is_node_in_zone(i, j, k, _domain_lim_from_rank(_node_rank(gl, i, j, k), gl))) {
          fatal_error("Problem with domain_lim_from_rank at i=%ld j=%ld k=%ld.",i,j,k);
        }
        if (_node_rank(gl, i, j, k)<0 || _node_rank(gl, i, j, k)>=proc){
          fatal_error("Node rank out of bounds at i=%ld j=%ld k=%ld.",i,j,k);
        } 
      end3DL
    end2DL
  end1DL


  for1DL (i,gl->domain_all.is,gl->domain_all.ie)
    for2DL (j,gl->domain_all.js,gl->domain_all.je)
      for3DL (k,gl->domain_all.ks,gl->domain_all.ke)
        zone=_domain_from_rank(_node_rank(gl, i, j, k), gl);
        if (!is_node_in_zone(i, j, k, zone)) {
          fatal_error("Problem with domain_from_rank at rank=%d i=%ld j=%ld k=%ld [zone.is=%ld zone.ie=%ld zone.js=%ld zone.je=%ld].",_node_rank(gl,i,j,k),i,j,k,zone.is,zone.ie,zone.js,zone.je);
        }
      end3DL
    end2DL
  end1DL

  MPI_Barrier(MPI_COMM_WORLD);

/*  int rank,proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    if (rank==1) display_node_type_window_local_process(stdout, *np_old, gl, TYPELEVEL_FLUID, 10, 24, 10, 40);
  MPI_Barrier(MPI_COMM_WORLD);
  exit(1);  */
#endif

}


