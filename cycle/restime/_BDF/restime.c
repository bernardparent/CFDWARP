#include <cycle/restime/_restime.h>
#include <cycle/share/cycle_share.h>
#include <cycle/restime/_restime.h>
#include <src/control.h>

#define STENCIL_BDF1 1
#define STENCIL_BDF2 2
#define STENCIL_BDF3 3


void write_disc_restime_template(FILE **controlfile){
  wfprintf(*controlfile,
    "  %s(\n"
    "    STENCIL=STENCIL_BDF3;\n"
    "  );\n"
  ,_RESTIME_ACTIONNAME);

}

void read_disc_restime_actions_2(char *actionname, char **argum, SOAP_codex_t *codex){
}


void read_disc_restime_actions(char *actionname, char **argum, SOAP_codex_t *codex){

  long numvarsinit;
  void (*action_original) (char *, char **, struct SOAP_codex_t *);
  gl_t *gl=((readcontrolarg_t *)codex->action_args)->gl;
  if (strcmp(actionname,_RESTIME_ACTIONNAME)==0) {
    SOAP_count_all_vars(codex, &numvarsinit);
    if (((readcontrolarg_t *)codex->action_args)->VERBOSE) wfprintf(stdout,"%s..",_RESTIME_ACTIONNAME);

    SOAP_add_int_to_vars(codex,"STENCIL_BDF1",STENCIL_BDF1); 
    SOAP_add_int_to_vars(codex,"STENCIL_BDF2",STENCIL_BDF2); 
    SOAP_add_int_to_vars(codex,"STENCIL_BDF3",STENCIL_BDF3); 

    gl->DISC_RESTIME_READ=TRUE;
 
    action_original=codex->action;
    codex->action=&read_disc_restime_actions_2;
    SOAP_process_code(*argum, codex, SOAP_VARS_KEEP_ALL);
    codex->action=action_original;

    find_int_var_from_codex(codex,"STENCIL",&gl->cycle.restime.STENCIL);
    if (gl->cycle.restime.STENCIL!=STENCIL_BDF1 && gl->cycle.restime.STENCIL!=STENCIL_BDF2 && gl->cycle.restime.STENCIL!=STENCIL_BDF3)
      SOAP_fatal_error(codex,"STENCIL must be set to either STENCIL_BDF1, STENCIL_BDF2, STENCIL_BDF3.");

    SOAP_clean_added_vars(codex,numvarsinit);
    codex->ACTIONPROCESSED=TRUE;
  }

}



void add_Z_dUstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  long l;
  long flux;
  flux_t tmp,dRes,Up1h,Um1h;
  sqmat_t Z;

    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      for (flux=0; flux<nf; flux++) {
        switch (gl->cycle.restime.STENCIL){
          case STENCIL_BDF1:
            Up1h[flux]= 
              +1.0*np[l].bs->U[flux]; 
            Um1h[flux]= 
              +1.0*np[l].bs->Um1[flux]; 
          break;
          case STENCIL_BDF2:
            Up1h[flux]= 
              +1.5*np[l].bs->U[flux]-0.5*np[l].bs->Um1[flux]; 
            Um1h[flux]= 
              +1.5*np[l].bs->Um1[flux]-0.5*np[l].bs->Um2[flux]; 
          break;
          case STENCIL_BDF3:
            Up1h[flux]= 
              +(11.0/6.0)*np[l].bs->U[flux]-(7.0/6.0)*np[l].bs->Um1[flux]+(1.0/3.0)*np[l].bs->Um2[flux]; 
            Um1h[flux]= 
              +(11.0/6.0)*np[l].bs->Um1[flux]-(7.0/6.0)*np[l].bs->Um2[flux]+(1.0/3.0)*np[l].bs->Um3[flux]; 
          break;
          default:
            fatal_error("STENCIL can not be set to %ld in add_Z_dUstar_residual().",gl->cycle.restime.STENCIL);
        }
      }
        
      for (flux=0; flux<nf; flux++){
        tmp[flux]=_Omega(np[l],gl)*(Up1h[flux]-Um1h[flux])/gl->dt;
      }
      
      find_Z(np,gl,l,Z);
      multiply_matrix_and_vector(Z, tmp, dRes);
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]+=dRes[flux];
    }
}



