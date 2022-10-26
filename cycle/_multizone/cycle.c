// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2000-2002 Bernard Parent

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

#include <cycle/ts/_ts.h>
#include <cycle/tsemf/_tsemf.h>
#include <cycle/res/_res.h>
#include <cycle/_cycle.h>
#include <cycle/share/cycle_share.h>
#include <src/bdry.h>
#include <model/fluid/_fluid.h>
#include <model/emfield/_emfield.h>

#ifdef DISTMPI
#error the multizone cycle is not compatible with distributed memory MPI
#endif


#ifdef EMFIELD
//#error the multizone cycle is not yet compatible with emfield 
#endif







void write_runtime_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n",_CYCLE_ACTIONNAME);

#ifdef UNSTEADY
  wfprintf(*controlfile,
  "    tmax=1e-4;    {simulation time span in seconds}\n"
  "    dt=tmax/30;   {physical time step in seconds}\n");
#else
  wfprintf(*controlfile,
  "    if (mod(iter,50)==0,\n"
  "      WriteDataFile(outputfilename);\n"
  "    );\n");
#endif
  wfprintf(*controlfile,
  "    sigma1=0.3;   {exponent of maximum pseudotime step in all dimensions}\n"
  "    sigma2=0.0;   {exponent of maximum pseudotime step in all fluxes}\n"
  "    if (iter==1,\n"
  "      CFL=0.1;\n"
  "    );\n"
  "    CFL=min(CFL*1.05,0.5);\n"
  "    phi1=20;      {maximum number of gridlines per zone in each dimension}\n"
  "    UpdateFluid(CFL,PRECON_LOCALTIMESTEP,sigma1,sigma2,phi1"
#ifdef UNSTEADY
  ",dt"
#endif
  ");\n\n"
#ifdef EMFIELD
  "    Lc=3e-3;      {characteristic length scale in meters used to solve the emfield equation}\n"
  "    relaxEMF=0.3; {relaxation factor forced on the update of the emfield variables}\n"
  "    UpdateEMField(Lc,relaxEMF"
#ifdef UNSTEADY
  ",dt"
#endif
  "); {optional extra parameters:  tsemfmethod [TSEMF_ADI, TSEMF_DDADI, TSEMF_IMAF, etc] and numsubiter_tsemf [default: 4] if supported}\n"
#endif
  );

  wfprintf(*controlfile,
  "    printf(\"%%6.2f %%6ld %%9.3f   %%E (%%2ld %%4ld"
#ifdef _2DL
                  ",%%4ld"
#endif
#ifdef _3DL
                  ",%%4ld"
#endif
                  ") %%3ld/%%3ld"
#ifdef EMFIELD
		  "  %%E (%%4ld"
#ifdef _2DL
		  ",%%4ld"
#endif
#ifdef _3DL
		  ",%%4ld"
#endif
                  ")"
#endif
		    "%%s \\n\",\n"
  "           CFL,iter,effiter_U,ximax,flux_ximax,i_ximax,\n"
  "           "
#ifdef _2DL
                  "j_ximax,"
#endif
#ifdef _3DL
                  "k_ximax,"
#endif
  "numzones_updated,numzones_total,"
#ifdef EMFIELD
		  "ximax_emfield,i_ximax_emfield,"
#ifdef _2DL
		  "j_ximax_emfield,"
#endif
#ifdef _3DL
		  "k_ximax_emfield,"
#endif
#endif

                  "clipinfo);\n");

#ifdef UNSTEADY
  wfprintf(*controlfile,
  "    if (ximax<xiverge"
#ifdef EMFIELD
  " && ximax_emfield<xiverge_emfield"
#endif
  ",\n"
  "      IncreaseTimeLevel();\n"
  "      WriteDataFile(outputfilename\".\"round(time/dt));\n"
  "      printf(\"dt=%%E time=%%Es\\n\",dt,time);\n"
  "      if (time>tmax-0.1*dt,\n"
  "        exit(EXIT_SUCCESS);\n"
  "      );\n"
  "    );\n");
#endif


  wfprintf(*controlfile,
  "  );\n");
  
}



void runtime_actions_cycle_specific(char *actionname, char **argum, SOAP_codex_t *codex){
  multizone_t multizone;
  gl_t *gl;
  np_t *np;
  bool UPDATE_ALL_ZONES;
#if defined(UNSTEADY) && !defined(EMFIELD)
//  zone_t zone;
#endif


  np=*((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (strcmp(actionname,"UpdateFluid")==0) {
    SOAP_substitute_all_argums(argum, codex);
#ifdef UNSTEADY
    if (SOAP_number_argums(*argum)!=6) SOAP_fatal_error(codex,"Number of arguments not equal to 6 in UpdateFluid(); action. This is likely due to the Multizone() module within the control file being valid for steady-state cases, while CFDWARP is here compiled for time accurate cases.");
    gl->dt=SOAP_get_argum_double(codex,*argum,5);
    if (gl->dt<=0.0) fatal_error("The time step dt must be positive when calling UpdateFluid.");
#endif

    gl->CFL=SOAP_get_argum_double(codex,*argum,0);
    if (gl->CFL<0.0) fatal_error("The CFL number can not be negative when calling UpdateFluid.");
    gl->PRECONDITIONER=SOAP_get_argum_long(codex,*argum,1);
    if (   gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP && gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP2 
        &&  gl->PRECONDITIONER!=PRECON_CONSTANTTIMESTEP
        && gl->PRECONDITIONER!=PRECON_LOCALEIGENVALUE && gl->PRECONDITIONER!=PRECON_LOCALEIGENVALUE2)
      fatal_error("Problem reading 2nd argument of UpdateFluid: The preconditioner must be set to either PRECON_LOCALTIMESTEP, PRECON_LOCALTIMESTEP2, PRECON_CONSTANTTIMESTEP, or PRECON_LOCALEIGENVALUE, or PRECON_LOCALEIGENVALUE2.");
    gl->sigma1=SOAP_get_argum_double(codex,*argum,2);
    gl->sigma2=SOAP_get_argum_double(codex,*argum,3);
    if (gl->sigma1<0.0 || gl->sigma1>1.0) fatal_error("The constant sigma1 must be within 0 and 1 when calling UpdateFluid.");
    if (gl->sigma2<0.0 || gl->sigma2>1.0) fatal_error("The constant sigma2 must be within 0 and 1 when calling UpdateFluid.");
    gl->cycle.phi1=SOAP_get_argum_long(codex,*argum,4);
    if (gl->cycle.phi1<=0) fatal_error("The constant phi1 must be positive when calling UpdateFluid.");
    if (gl->PRECONDITIONER==PRECON_CONSTANTTIMESTEP) gl->dtau=SOAP_get_argum_double(codex,*argum,5);


#ifdef EMFIELD
    if (gl->RESIDUAL_ALTERED_EMFIELD ){
      update_bdry_nodes_emfield(np,gl,gl->domain); //????
      gl->RESIDUAL_ALTERED_EMFIELD=FALSE;
    }
    update_prim_emfield_mem_in_zone(np,gl,gl->domain);
    UPDATE_ALL_ZONES=TRUE;
    setup_multizone(np, gl, gl->domain, gl->domain, gl->cycle.fluid.xiverge, gl->cycle.phi1, UPDATE_ALL_ZONES, &multizone);
    find_residual_with_multizone(np, gl, multizone);
    find_ximax(np, gl, gl->domain, IJK_UPDATE_YES);
    UPDATE_ALL_ZONES=FALSE;
    setup_multizone(np, gl, gl->domain, gl->domain, gl->cycle.fluid.xiverge, gl->cycle.phi1, UPDATE_ALL_ZONES, &multizone);
    add_int_to_codex(&(gl->cycle.codex), "numzones_updated", multizone.numzones_ts);
    add_int_to_codex(&(gl->cycle.codex), "numzones_total", multizone.numzones_total);
    update_U_with_multizone(np, gl, multizone);
    update_linked_nodes(np, gl, TYPELEVEL_FLUID_WORK);
    update_bdry_nodes_with_multizone(np, gl, multizone);
    free_multizone(&multizone);
#endif

#ifndef EMFIELD
    find_ximax(np, gl, gl->domain, IJK_UPDATE_YES);
    if (gl->RESIDUAL_ALTERED){
      UPDATE_ALL_ZONES=TRUE;
    } else {
      UPDATE_ALL_ZONES=FALSE;
    }
    setup_multizone(np, gl, gl->domain, gl->domain, gl->cycle.fluid.xiverge, gl->cycle.phi1, UPDATE_ALL_ZONES, &multizone);
    add_int_to_codex(&(gl->cycle.codex), "numzones_updated", multizone.numzones_ts);
    add_int_to_codex(&(gl->cycle.codex), "numzones_total", multizone.numzones_total);
    update_U_with_multizone(np,gl,multizone);
    update_linked_nodes(np, gl, TYPELEVEL_FLUID_WORK);
    update_bdry_nodes_with_multizone(np,gl,multizone);
    find_residual_with_multizone(np,gl,multizone);
    free_multizone(&multizone);
    gl->RESIDUAL_ALTERED=FALSE;
#endif


    update_runtime_codex_xi_from_gl(gl,codex);
    codex->ACTIONPROCESSED=TRUE;
  }


#ifdef UNSTEADY
  if (strcmp(actionname,"IncreaseTimeLevel")==0) {
  
    increase_time_level(np, gl);
/*  #ifndef EMFIELD
    zone=_zone_intersection(gl->domain_all,_zone_expansion(gl->domain,+hbw_res_fluid));  
    resume_nodes_only_in_zone_and_update_bdry_nodes(np, gl, zone);
    update_bdry_nodes(np, gl, gl->domain);
    find_residual(np, gl, gl->domain);
  #endif
*/
  #ifndef EMFIELD
    UPDATE_ALL_ZONES=TRUE;
    setup_multizone(np, gl, gl->domain, gl->domain, gl->cycle.fluid.xiverge, gl->cycle.phi1, UPDATE_ALL_ZONES, &multizone);
    update_bdry_nodes_with_multizone(np,gl,multizone);
    find_residual_with_multizone(np,gl,multizone);
    free_multizone(&multizone);
    gl->RESIDUAL_ALTERED=FALSE;
  #endif
    codex->ACTIONPROCESSED=TRUE;
  }
#endif


#ifdef EMFIELD
  if (strcmp(actionname,"UpdateEMField")==0) {
    read_UpdateEMField_arguments(argum,codex,gl);
    if (gl->RESIDUAL_ALTERED_EMFIELD ){
      update_bdry_nodes_emfield(np,gl,gl->domain); //????
      gl->RESIDUAL_ALTERED_EMFIELD=FALSE;
    }

    update_prim_emfield_mem_in_zone(np,gl,gl->domain);
    find_residual_emfield(np, gl, gl->domain);
    find_ximax_emfield(np, gl, gl->domain);
    update_U_emfield(np, gl, gl->domain);
    update_linked_nodes(np, gl, TYPELEVEL_EMFIELD);
    update_bdry_nodes_emfield(np,gl,gl->domain);
    update_runtime_codex_xi_from_gl(gl,codex);
    codex->ACTIONPROCESSED=TRUE;
  }  
#endif


}



void init_cycle(char *argum, SOAP_codex_t *codexcontrol){
  np_t **np;
  gl_t *gl;
  zone_t zone;

  np=((readcontrolarg_t *)codexcontrol->action_args)->np;
  gl=((readcontrolarg_t *)codexcontrol->action_args)->gl;


  
#ifdef EMFIELD
  zone=_zone_intersection(gl->domain_all,_zone_expansion(gl->domain,+hbw_res_fluid));
  update_prim_emfield_mem_in_zone(*np, gl, zone);
  resume_nodes_only_in_zone_and_update_bdry_nodes(*np, gl, zone);
  SOAP_copy_codex(codexcontrol,&(gl->cycle.codex));
  update_prim_emfield_mem_in_zone(*np, gl, zone);
  find_residual_emfield(*np, gl, gl->domain);
  find_ximax_emfield(*np, gl, gl->domain);
  find_residual(*np, gl, zone);
  find_ximax(*np,gl,zone,IJK_UPDATE_YES);
  gl->RESIDUAL_ALTERED_EMFIELD=FALSE;
#else
  zone=_zone_intersection(gl->domain_all,_zone_expansion(gl->domain,+hbw_res_fluid));
  resume_nodes_only_in_zone_and_update_bdry_nodes(*np, gl, zone);
  update_bdry_nodes(*np, gl, gl->domain);
  find_residual(*np, gl, gl->domain);
  find_ximax(*np, gl, gl->domain, IJK_UPDATE_YES);
  SOAP_copy_codex(codexcontrol,&(gl->cycle.codex));
#endif
  gl->RESIDUAL_ALTERED=FALSE;
}


void perform_one_iteration(np_t *np, gl_t *gl){
#ifndef UNSTEADY
  double ximax_max;
#endif

#ifndef EMFIELD
  update_runtime_codex_xi_from_gl(gl,&(gl->cycle.codex));
#endif
  
  process_code_runtime(np,gl,gl->cycle.code_runtime,&(gl->cycle.codex));

#ifndef UNSTEADY

#ifdef DISTMPI
  MPI_Allreduce(&gl->ximax, &ximax_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  ximax_max=gl->ximax;
#endif

  if (ximax_max<=gl->cycle.fluid.xiverge
#ifdef EMFIELD
      && gl->ximax_emfield<=gl->cycle.emfield.xiverge
#endif
     ) {
      gl->CONVERGED=TRUE;
  }

#endif
}

