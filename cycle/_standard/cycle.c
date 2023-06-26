// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 1998-2016 Bernard Parent

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



void write_runtime_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n",_CYCLE_ACTIONNAME);
  wfprintf(*controlfile,
#ifdef UNSTEADY
  "    tmax=1e-4;    {simulation time span in seconds}\n"
  "    dt=tmax/30;   {physical time step in seconds}\n"
#else 
  "    if (mod(iter,50)==0,\n"
  "      WriteDataFile(outputfilename);\n"
  "    );\n"
#endif
  "    sigma1=0.3;   {exponent of maximum pseudotime step in all dimensions}\n"
  "    sigma2=0.0;   {exponent of maximum pseudotime step in all fluxes}\n"
  "    if (iter==1,\n"
  "      CFL=0.01;\n"
  "    );\n"
  "    CFL=min(CFL*1.05,0.5);\n"
  "    UpdateFluid(CFL,PRECON_LOCALTIMESTEP,sigma1,sigma2"
#ifdef UNSTEADY
    ",dt"
#endif
  ");\n"
#ifdef EMFIELD
  "    Lc=1e0;      {characteristic length scale in meters used to solve the emfield equation}\n"
  "    relaxEMF=0.8; {relaxation factor forced on the update of the emfield variables}\n"
  "    UpdateEMField(Lc,relaxEMF"
#ifdef UNSTEADY
    ",dt"
#endif
  ",TSEMF_SOR,500); {optional extra parameters:  tsemfmethod [TSEMF_ADI, TSEMF_DDADI, TSEMF_IMAF, etc] and numsubiter_tsemf [default: 4] if supported}\n"
#endif
  "    printf(\"%%6.2f %%6ld %%9.3f   %%E (%%2ld %%4ld"
#ifdef _2DL
                  ",%%4ld"
#endif
#ifdef _3DL
                  ",%%4ld"
#endif
                  ")"
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
                  "  %%s\\n\",\n"
  "           CFL,iter,effiter_U,ximax,flux_ximax,i_ximax,\n"
  "           "
#ifdef _2DL
                  "j_ximax,"
#endif
#ifdef _3DL
                  "k_ximax,"
#endif
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
  "    );\n"
  );
#endif

  wfprintf(*controlfile,
  "  );\n");

}



void runtime_actions_cycle_specific(char *actionname, char **argum, SOAP_codex_t *codex){
#ifdef DISTMPI
  double delta_effiter_U,effiter_U_mem,delta_effiter_U_max;
#endif
  gl_t *gl;
  np_t *np;


  np=*((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (strcmp(actionname,"UpdateFluid")==0) {
    SOAP_substitute_all_argums(argum, codex);
#ifdef UNSTEADY
    if (SOAP_number_argums(*argum)!=5) SOAP_fatal_error(codex,"Number of arguments not equal to 5 in UpdateFluid(); action. This is likely due to the Standard() module within the control file being valid for steady-state cases, while CFDWARP is here compiled for time accurate cases.");
    gl->dt=SOAP_get_argum_double(codex,*argum,4);
    if (gl->dt<=0.0) fatal_error("The time step dt must be positive when calling UpdateFluid.");
#endif
    gl->CFL=SOAP_get_argum_double(codex,*argum,0);
    if (gl->CFL<0.0) fatal_error("The CFL number can not be negative when calling UpdateFluid.");
    gl->PRECONDITIONER=SOAP_get_argum_long(codex,*argum,1);
    if (   gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP && gl->PRECONDITIONER!=PRECON_LOCALTIMESTEP2 
        && gl->PRECONDITIONER!=PRECON_CONSTANTTIMESTEP
        && gl->PRECONDITIONER!=PRECON_LOCALEIGENVALUE && gl->PRECONDITIONER!=PRECON_LOCALEIGENVALUE2)
      fatal_error("Problem reading 2nd argument of UpdateFluid: The preconditioner must be set to either PRECON_LOCALTIMESTEP, PRECON_LOCALTIMESTEP2, PRECON_CONSTANTTIMESTEP, or PRECON_LOCALEIGENVALUE, or PRECON_LOCALEIGENVALUE2.");
    gl->sigma1=SOAP_get_argum_double(codex,*argum,2);
    gl->sigma2=SOAP_get_argum_double(codex,*argum,3);
    if (gl->sigma1<0.0 || gl->sigma1>1.0) fatal_error("The constant sigma1 must be within 0 and 1 when calling UpdateFluid.");
    if (gl->sigma2<0.0 || gl->sigma2>1.0) fatal_error("The constant sigma2 must be within 0 and 1 when calling UpdateFluid.");
    if (gl->PRECONDITIONER==PRECON_CONSTANTTIMESTEP) {
      gl->dtau=SOAP_get_argum_double(codex,*argum,4);
      if (gl->dtau<=0.0) fatal_error("The time step dtau must be positive when calling UpdateFluid.");
    }
#ifdef EMFIELD
    if (gl->RESIDUAL_ALTERED_EMFIELD){
      update_bdry_nodes_emfield(np,gl,gl->domain);  
      gl->RESIDUAL_ALTERED_EMFIELD=FALSE;
    }
    //update_prim_emfield_mem_in_zone(np,gl,_zone_intersection(gl->domain_all,_zone_expansion(gl->domain,hbw_mem_emfield)));
#endif
    if (gl->RESIDUAL_ALTERED){
      update_bdry_nodes(np,gl,gl->domain); 
      gl->RESIDUAL_ALTERED=FALSE;
    }
    find_residual(np, gl, gl->domain);
    find_ximax(np, gl, gl->domain, IJK_UPDATE_YES);
#ifdef DISTMPI
    effiter_U_mem=gl->effiter_U;
#endif
    update_U(np, gl, gl->domain);
    update_linked_nodes(np, gl, TYPELEVEL_FLUID);
#ifdef DISTMPI
    exchange_U(np, gl);
    delta_effiter_U=gl->effiter_U-effiter_U_mem;
    MPI_Allreduce(&delta_effiter_U, &delta_effiter_U_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    gl->effiter_U=effiter_U_mem+delta_effiter_U_max;
#endif
    update_bdry_nodes(np, gl, gl->domain);
    update_runtime_codex_xi_from_gl(gl,codex);
    codex->ACTIONPROCESSED=TRUE;
  }







#ifdef UNSTEADY
  if (strcmp(actionname,"IncreaseTimeLevel")==0) {
    increase_time_level(np, gl);
    codex->ACTIONPROCESSED=TRUE;
  }
#endif




#ifdef EMFIELD
  if (strcmp(actionname,"UpdateEMField")==0) {
    read_UpdateEMField_arguments(argum,codex,gl);
    if (gl->RESIDUAL_ALTERED_EMFIELD ){
      update_bdry_nodes_emfield(np,gl,gl->domain); 
      gl->RESIDUAL_ALTERED_EMFIELD=FALSE;
    }
    update_prim_emfield_mem_in_zone(np,gl,_zone_intersection(gl->domain_all,_zone_expansion(gl->domain,max(hbw_resplasma_fluid,hbw_mem_emfield)))); 
    find_residual_emfield(np, gl, gl->domain);
    find_ximax_emfield(np, gl, gl->domain);
    update_U_emfield(np, gl, gl->domain);
    update_linked_nodes(np, gl, TYPELEVEL_EMFIELD);
#ifdef DISTMPI
    exchange_U_emfield(np, gl);
#endif
    update_bdry_nodes_emfield(np,gl,gl->domain); 
    
    update_runtime_codex_xi_from_gl(gl,codex);
    codex->ACTIONPROCESSED=TRUE;
  }  


/*
  if (strcmp(actionname,"UpdateFluidEMField")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=6)
      SOAP_fatal_error(codex,"Number of arguments not equal to 6 in UpdateFluidEMField(); action.");
    gl->CFL=SOAP_get_argum_double(codex,*argum,0);
    gl->PRECONDITIONING=SOAP_get_argum_bool(codex,*argum,1);
    gl->sigma1=SOAP_get_argum_double(codex,*argum,2);
    gl->sigma2=SOAP_get_argum_double(codex,*argum,3);
    gl->Lc=SOAP_get_argum_double(codex,*argum,4);
    gl->relaxEMF=SOAP_get_argum_double(codex,*argum,5);
    gl->LOCAL_TIME_STEPPING=TRUE;

    if (gl->RESIDUAL_ALTERED_EMFIELD ){
      update_bdry_nodes_emfield(np,gl,gl->domain); 
      gl->RESIDUAL_ALTERED_EMFIELD=FALSE;
    }
    update_prim_emfield_mem_in_zone(np,gl,_zone_intersection(gl->domain_all,_zone_expansion(gl->domain,hbw_mem_emfield)));
    if (gl->RESIDUAL_ALTERED){
      update_bdry_nodes(np,gl,gl->domain); 
      gl->RESIDUAL_ALTERED=FALSE;
    }
    find_residual_emfield(np, gl, gl->domain);
    find_residual(np, gl, gl->domain);
    find_ximax_emfield(np, gl, gl->domain);
    find_ximax(np, gl, gl->domain, IJK_UPDATE_YES);
    update_U_emfield(np, gl, gl->domain);
    update_U(np, gl, gl->domain);
    update_bdry_nodes_emfield(np,gl,gl->domain);
    update_bdry_nodes(np, gl, gl->domain);



    update_runtime_codex_xi_from_gl(gl,codex);
    codex->ACTIONPROCESSED=TRUE;

  }
*/

#endif


}



void init_cycle(char *argum, SOAP_codex_t *codexcontrol){
  np_t **np;
  gl_t *gl;
  np=((readcontrolarg_t *)codexcontrol->action_args)->np;
  gl=((readcontrolarg_t *)codexcontrol->action_args)->gl;
#ifdef EMFIELD
  //update_prim_emfield_mem_in_zone(*np, gl, _zone_intersection(gl->domain_all,gl->domain_lim));
#endif

  resume_nodes_only_in_zone(*np, gl, _zone_intersection(gl->domain_all,_zone_expansion(gl->domain,+max(hbw_res_fluid,hbw_bdry_fluid))));
  SOAP_copy_codex(codexcontrol,&(gl->cycle.codex));
#ifdef EMFIELD
  update_prim_emfield_mem_in_zone(*np, gl, _zone_intersection(gl->domain_all,_zone_expansion(gl->domain,hbw_mem_emfield)));
#endif
  update_bdry_nodes(*np, gl, gl->domain);
  update_linked_nodes(*np, gl, TYPELEVEL_FLUID);
#ifdef EMFIELD
  update_prim_emfield_mem_in_zone(*np, gl, _zone_intersection(gl->domain_all,_zone_expansion(gl->domain,max(hbw_resplasma_fluid,hbw_mem_emfield)))); 
//  printf("\n%ld,%ld  %ld,%ld\n",gl->domain.is,gl->domain.ie,gl->domain_lim.is,gl->domain_lim.ie);
 // update_prim_emfield_mem_in_zone(*np, gl, gl->domain);
  find_residual_emfield(*np, gl, gl->domain);
  find_ximax_emfield(*np, gl, gl->domain);
#endif
  find_residual(*np, gl, gl->domain);
  find_ximax(*np,gl,gl->domain,IJK_UPDATE_YES);
  gl->RESIDUAL_ALTERED=FALSE;
#ifdef EMFIELD
  gl->RESIDUAL_ALTERED_EMFIELD=FALSE;
#endif
}




void perform_one_iteration(np_t *np, gl_t *gl){
#ifndef UNSTEADY
  double ximax_max;
 #ifdef EMFIELD
  double ximax_max_emfield;
 #endif
#endif


  process_code_runtime(np,gl,gl->cycle.code_runtime,&(gl->cycle.codex));


#ifndef UNSTEADY
#ifdef DISTMPI
  MPI_Allreduce(&gl->ximax, &ximax_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
 #ifdef EMFIELD
  MPI_Allreduce(&gl->ximax_emfield, &ximax_max_emfield, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
 #endif
#else
  ximax_max=gl->ximax;
 #ifdef EMFIELD
  ximax_max_emfield=gl->ximax_emfield;
 #endif
#endif

  if (ximax_max<=gl->cycle.fluid.xiverge
#ifdef EMFIELD
      && ximax_max_emfield<=gl->cycle.emfield.xiverge
#endif
     ) {
      gl->CONVERGED=TRUE;
  } 
#endif

}

