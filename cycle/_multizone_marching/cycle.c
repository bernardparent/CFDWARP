#include <cycle/ts/_ts.h>
#include <cycle/tsemf/_tsemf.h>
#include <cycle/res/_res.h>
#include <cycle/_cycle.h>
#include <cycle/share/cycle_share.h>
#include <src/bdry.h>
#include <model/fluid/_fluid.h>
#include <model/emfield/_emfield.h>


#define minwidth (hbw_res_fluid+1)

#ifdef DISTMPI
#error the multizone-marching cycle is not compatible with distributed MPI
#endif

#ifdef EMFIELD
#error the multizone-marching cycle is not yet compatible with emfield 
#endif

#ifdef UNSTEADY
#error the multizone-marching cycle can not be used for time-accurate simulations
#endif



void write_runtime_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n",_CYCLE_ACTIONNAME);
  wfprintf(*controlfile,
  "    if (mod(iter,50)==0,\n"
  "      WriteDataFile(outputfilename);\n"
  "    );\n"
  "    sigma1=0.3;      {exponent of maximum pseudotime step in all dimensions}\n"
  "    sigma2=0.3;      {exponent of maximum pseudotime step in all fluxes}\n"
  "    if (iter==1,\n"
  "      CFL=0.1;\n"
  "    );\n"
  "    CFL=min(CFL*1.05,0.5);\n"
  "    phi1=20;         {maximum number of gridlines per zone in each dimension}\n"
  "    phi2=6;          {minimum number of iterations before advancing the marching-window downstream boundary}\n"
  "    phi3=17;         {minimum width of the marching window}\n"
  "    varphiverge=1e4; {parameter related to the streamwise ellipticity sensor}\n"
  "    entrance=1;      {force the marching window upstream boundary to be at least as high as entrance}\n"
  "    UpdateFluid(CFL,PRECON_LOCALTIMESTEP,sigma1,sigma2,phi1,phi2,phi3,varphiverge,entrance);\n\n"
  );

  wfprintf(*controlfile,
  "    printf(\"%%6.2f %%4ld %%4ld %%6ld %%9.3f   %%E (%%4ld"
#ifdef _2DL
               ",%%4ld"
#endif
#ifdef _3DL
               ",%%4ld"
#endif
                  ") %%3ld/%%3ld  "
#ifdef EMFIELD
		  "  %%6.0f %%E (%%4ld"
#ifdef _2DL
		  ",%%4ld"
#endif
#ifdef _3DL
		  ",%%4ld"
#endif
                  ")"
#endif
		  "  %%s \\n\",\n"
  "           CFL,window.is,window.ie,iter,effiter_U,ximax,i_ximax,\n"
  "           "
#ifdef _2DL
               "j_ximax,"
#endif
#ifdef _3DL
               "k_ximax,"
#endif
                  "numzones_updated,numzones_total,"
#ifdef EMFIELD
		  "Lc,ximax_emfield,i_ximax_emfield,"
#ifdef _2DL
		  "j_ximax_emfield,"
#endif
#ifdef _3DL
		  "k_ximax_emfield,"
#endif
#endif
		  "clipinfo);\n"
  "  );\n");

}





static void resize_marching_window(np_t *np, gl_t *gl, zone_t oldwindow, zone_t newwindow){
  long i,j,k,flux;
  zone_t zone;


  zone=gl->domain;
  zone.ie=min(gl->domain.ie,newwindow.ie);
  zone.is=max(gl->domain.is,newwindow.is-hbw_res_fluid-hbw_bdry_fluid-hbw_res_fluid-hbw_bdry_fluid);
  resume_nodes_only_in_zone_and_update_bdry_nodes(np, gl, zone);
  zone=gl->domain_lim;
  zone.is=zone.ie=newwindow.ie;
  update_inner_to_bdry_node_type(np,gl, TYPELEVEL_FLUID_WORK, NODETYPE_BDRYOUTFLOW, zone);
  zone.is=newwindow.ie+1;
  zone.ie=zone.is+2;
  update_node_type(np, gl, TYPELEVEL_FLUID_WORK, NODETYPE_UNUSED, zone);
  zone.ie=min(gl->domain.ie, newwindow.ie);
  zone.is=max(gl->domain.is, newwindow.is-hbw_res_fluid-hbw_bdry_fluid-max(hbw_res_fluid,hbw_bdry_fluid));
  adjust_node_type(np, gl, zone, TYPELEVEL_FLUID_WORK);

  if (newwindow.is<oldwindow.is){
    for1DL(i,max(gl->domain.is,newwindow.is-hbw_res_fluid-hbw_bdry_fluid-max(hbw_res_fluid,hbw_bdry_fluid)),
             max(gl->domain.is,oldwindow.is-hbw_bdry_fluid-hbw_res_fluid-1))
      for2DL(j,gl->domain.js,gl->domain.je)
        for3DL(k,gl->domain.ks,gl->domain.ke)
          if (is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)) {
            assert(is_node_resumed(np[_ai(gl,i,j,k)]));
            np[_ai(gl,i,j,k)].wk->xi=0.0e0;
            for (flux=0; flux<nf; flux++) np[_ai(gl,i,j,k)].wk->Res[flux]=0.0e0;
          }
        end3DL
      end2DL
    end1DL
  }

  if (oldwindow.ie!=newwindow.ie) {
    zone=newwindow;
    zone.is=max(oldwindow.ie-hbw_bdry_fluid, /* hbw_bdry_fluid needed here for sudden opening */
            max(gl->domain.is,newwindow.is-hbw_res_fluid-hbw_bdry_fluid));
    zone.ie=newwindow.ie;
    update_bdry_nodes(np, gl, zone);
    zone=newwindow;
    zone.is=max(oldwindow.ie-hbw_bdry_fluid-hbw_res_fluid, max(gl->domain.is,newwindow.is-hbw_res_fluid-hbw_bdry_fluid));
    zone.ie=newwindow.ie;
    find_residual(np, gl, zone);
  }

  if (newwindow.is<oldwindow.is-hbw_res_fluid-hbw_bdry_fluid) {
    zone=newwindow;
    zone.is=newwindow.is;
    zone.ie=max(zone.is,oldwindow.is-hbw_res_fluid-hbw_bdry_fluid-1);
    find_residual(np, gl, zone);
  }
}




void init_cycle(char *argum, SOAP_codex_t *codexcontrol){
  np_t **np;
  gl_t *gl;
  zone_t oldwindow;

  np=((readcontrolarg_t *)codexcontrol->action_args)->np;
  gl=((readcontrolarg_t *)codexcontrol->action_args)->gl;



  gl->window.js=gl->domain.js;
  gl->window.je=gl->domain.je;
  gl->window.ks=gl->domain.ks;
  gl->window.ke=gl->domain.ke;
  if (gl->iter==1) {
    gl->window=gl->domain;
    gl->window.ie=gl->window.is+minwidth;
  }
  //gl->window.is=max(gl->cycle.entrance,gl->window.is);
  gl->window.ie=max(gl->window.is+minwidth,gl->window.ie);
  oldwindow=gl->window;
  oldwindow.is=0;
  oldwindow.ie=0;
  resize_marching_window(*np, gl, oldwindow, gl->window); 
  find_ximax(*np, gl, gl->window, IJK_UPDATE_YES);
  SOAP_copy_codex(codexcontrol,&(gl->cycle.codex));
  gl->RESIDUAL_ALTERED=FALSE;

}




static long _max_i_elliptic(np_t *np, gl_t *gl, zone_t zone){
  long i,j,k;
  bool FOUND_ELLIPTIC;

  i=zone.ie+1;
  FOUND_ELLIPTIC=FALSE;
  do{
    i--;
    for2DL(j,zone.js,zone.je)
      for3DL(k,zone.ks,zone.ke)
        if (  is_node_inner(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID_WORK)
           && (   _varphi(np,gl,_ai(gl,i,j,k),0)>gl->cycle.varphiverge
              )
	   ) FOUND_ELLIPTIC=TRUE;
      end3DL
    end2DL
  } while (!FOUND_ELLIPTIC && i>zone.is);
  return(i);
}




void runtime_actions_cycle_specific(char *actionname, char **argum, SOAP_codex_t *codex){
  long maxIelliptic;
  multizone_t multizone;
  zone_t zonelim,newwindow;
  gl_t *gl;
  np_t *np;
  bool UPDATE_ALL_ZONES;


  np=*((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (strcmp(actionname,"UpdateFluid")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=9)
      SOAP_fatal_error(codex,"Number of arguments not equal to 9 in UpdateFluid(); action.");
    gl->CFL=SOAP_get_argum_double(codex,*argum,0);
    if (gl->CFL<0.0) fatal_error("The CFL number can not be negative when calling UpdateFluid.");
    gl->PRECONDITIONER=SOAP_get_argum_long(codex,*argum,1);
    if (   gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP && gl->PRECONDITIONER!=PRECON_CONSTANTTIMESTEP
        && gl->PRECONDITIONER!=PRECON_LOCALEIGENVALUE && gl->PRECONDITIONER!=PRECON_LOCALEIGENVALUE2)
      fatal_error("Problem reading 2nd argument of UpdateFluid: The preconditioner must be set to either PRECON_LOCALTIMESTEP, PRECON_CONSTANTTIMESTEP, or PRECON_LOCALEIGENVALUE, or PRECON_LOCALEIGENVALUE2.");
    gl->sigma1=SOAP_get_argum_double(codex,*argum,2);
    gl->sigma2=SOAP_get_argum_double(codex,*argum,3);
    if (gl->sigma1<0.0 || gl->sigma1>1.0) fatal_error("The constant sigma1 must be within 0 and 1 when calling UpdateFluid.");
    if (gl->sigma2<0.0 || gl->sigma2>1.0) fatal_error("The constant sigma2 must be within 0 and 1 when calling UpdateFluid.");
    gl->cycle.phi1=SOAP_get_argum_long(codex,*argum,4);
    gl->cycle.phi2=SOAP_get_argum_long(codex,*argum,5);
    gl->cycle.phi3=SOAP_get_argum_long(codex,*argum,6);
    gl->cycle.varphiverge=SOAP_get_argum_double(codex,*argum,7);
    gl->cycle.entrance=SOAP_get_argum_long(codex,*argum,8);
    if (gl->cycle.phi1<=0) fatal_error("The constant phi1 must be positive when calling UpdateFluid.");
    if (gl->cycle.phi2<=0) fatal_error("The constant phi2 must be positive when calling UpdateFluid.");
    if (gl->cycle.phi3<=0) fatal_error("The constant phi3 must be positive when calling UpdateFluid.");
    if (gl->cycle.varphiverge<=0.0) fatal_error("The constant varphiverge must be positive when calling UpdateFluid.");
    if (gl->PRECONDITIONER==PRECON_CONSTANTTIMESTEP) gl->dtau=SOAP_get_argum_double(codex,*argum,9);


    zonelim=gl->window;
    zonelim.is=gl->domain.is;
    if (gl->RESIDUAL_ALTERED){
      UPDATE_ALL_ZONES=TRUE;
      newwindow=gl->window;
      newwindow.is=gl->domain.is;
      resize_marching_window(np,gl,gl->window,newwindow);
      gl->window=newwindow;
    } else {
      UPDATE_ALL_ZONES=FALSE;
    }


    setup_multizone(np, gl, gl->window, zonelim, gl->cycle.fluid.xiverge, gl->cycle.phi1, UPDATE_ALL_ZONES, &multizone);

    find_ximax(np, gl, gl->window, IJK_UPDATE_YES);
    add_int_to_codex(&(gl->cycle.codex), "window.is", gl->window.is);
    add_int_to_codex(&(gl->cycle.codex), "window.ie", gl->window.ie);
    add_int_to_codex(&(gl->cycle.codex), "numzones_updated", multizone.numzones_ts);
    add_int_to_codex(&(gl->cycle.codex), "numzones_total", multizone.numzones_total);
    //update_runtime_codex_xi_from_gl(gl,codex);


    update_U_with_multizone(np,gl,multizone);
    update_linked_nodes(np, gl, TYPELEVEL_FLUID_WORK);
    update_bdry_nodes_with_multizone(np,gl,multizone);
    find_residual_with_multizone(np,gl,multizone);
    free_multizone(&multizone);
    gl->RESIDUAL_ALTERED=FALSE;

    /* determine new left and right boundaries */
    newwindow=gl->window;
    newwindow.is=max(gl->domain.is,gl->window.is-hbw_res_fluid-hbw_bdry_fluid-1);
    //if (gl->window.is==gl->cycle.entrance) newwindow.is=gl->cycle.entrance;
    newwindow.is=max(newwindow.is,min(gl->cycle.entrance-1,gl->window.ie-2*hbw_res_fluid));
    do {
      newwindow.is++;
      newwindow.ie=newwindow.is;
      find_ximax(np, gl, newwindow, IJK_UPDATE_NO);
    } while(gl->ximax<gl->cycle.fluid.xiverge && newwindow.is<gl->window.ie-minwidth);
    if (newwindow.is<gl->cycle.entrance && gl->window.is>=gl->cycle.entrance)
      newwindow.is=gl->cycle.entrance;
    if (newwindow.is<gl->window.is+1 && newwindow.is>gl->window.is)
       newwindow.is=gl->window.is;
    newwindow.ie=gl->window.ie;
    if (mod(gl->iter,gl->cycle.phi2)==0){
      if (gl->window.ie-gl->window.is>=gl->cycle.phi3) maxIelliptic=_max_i_elliptic(np,gl,gl->window);
        else maxIelliptic=gl->window.is;

      if (gl->window.ie-max(gl->window.is,maxIelliptic)<gl->cycle.phi3)
        newwindow.ie=min(gl->domain.ie,gl->window.ie+1);
    }
    if (newwindow.ie!=gl->window.ie || newwindow.is!=gl->window.is)
       resize_marching_window(np,gl,gl->window,newwindow);
    gl->window=newwindow;
    find_ximax(np, gl, gl->window, IJK_UPDATE_YES);
    update_runtime_codex_xi_from_gl(gl,codex);
    if (gl->ximax<gl->cycle.fluid.xiverge && gl->window.ie==gl->domain.ie) gl->CONVERGED=TRUE;
    codex->ACTIONPROCESSED=TRUE;

  }


}





void perform_one_iteration(np_t *np, gl_t *gl){
  double ximax_max;

  process_code_runtime(np,gl,gl->cycle.code_runtime,&(gl->cycle.codex));


#ifdef DISTMPI
  MPI_Allreduce(&gl->ximax, &ximax_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  ximax_max=gl->ximax;
#endif

  if (ximax_max<=gl->cycle.fluid.xiverge && gl->window.ie==gl->domain.ie) {
      gl->CONVERGED=TRUE;
  }

}

