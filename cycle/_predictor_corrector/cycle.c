// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2017-2018 Bernard Parent

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
#include <cycle/restime/_restime.h>
#include <src/bdry.h>

#ifdef EMFIELD
#error The Predictor Corrector cycle can not be used with EMFIELD
#endif

#ifndef _RESTIME_PREDICTOR_CORRECTOR
#error The Predictor Corrector cycle can only be used with a predictor corrector time discretization method (RK2, RK3, etc)
#endif

void write_runtime_template(FILE **controlfile){
  wfprintf(*controlfile,
  "  %s(\n",_CYCLE_ACTIONNAME);
  wfprintf(*controlfile,
  "    tmax=1e-4;    {simulation time span in seconds}\n"
  "    dt=tmax/30;   {physical time step in seconds}\n"
  );
  wfprintf(*controlfile,
  "    if (time>tmax-0.1*dt,\n"
  "      WriteDataFile(outputfilename\".\"round(time/dt));\n"
  "      exit(EXIT_SUCCESS);\n"
  "    );\n"
  "    if (mod(iter,50)==1,\n"
  "      WriteDataFile(outputfilename\".\"round(time/dt));\n"
  "    );\n"
  );
  wfprintf(*controlfile,
  "    UpdateFluid(dt);\n"
  "    printf(\"%%E %%6ld %%9.3f   %%E (%%4ld"
#ifdef _2DL
                  ",%%4ld"
#endif
#ifdef _3DL
                  ",%%4ld"
#endif
                  ")"
                  "  %%s\\n\",\n"
  "           time,iter,effiter_U,ximax,i_ximax,\n"
  "           "
#ifdef _2DL
                  "j_ximax,"
#endif
#ifdef _3DL
                  "k_ximax,"
#endif
                  "clipinfo);\n");


  wfprintf(*controlfile,
  "  );\n");

}



void runtime_actions_cycle_specific(char *actionname, char **argum, SOAP_codex_t *codex){
#ifdef DISTMPI
  double delta_effiter_U,effiter_U_mem,delta_effiter_U_max;
#endif
  long l,flux,i,j,k;
  gl_t *gl;
  np_t *np;


  np=*((readcontrolarg_t *)codex->action_args)->np;
  gl=((readcontrolarg_t *)codex->action_args)->gl;

  if (strcmp(actionname,"UpdateFluid")==0) {
    SOAP_substitute_all_argums(argum, codex);
    if (SOAP_number_argums(*argum)!=1) SOAP_fatal_error(codex,"Number of arguments not equal to 1 in UpdateFluid(); action.");
    gl->dt=SOAP_get_argum_double(codex,*argum,0);
    if (gl->dt<=0.0) fatal_error("The time step dt must be positive when calling UpdateFluid.");
    for1DL(i,gl->domain_lim.is,gl->domain_lim.ie)
      for2DL(j,gl->domain_lim.js,gl->domain_lim.je)
        for3DL(k,gl->domain_lim.ks,gl->domain_lim.ke)
          l=_ai(gl,i,j,k);
          if ((is_node_valid(np[_ai(gl,i,j,k)],TYPELEVEL_FLUID))){
            for (flux=0; flux<nf; flux++){
              np[l].bs->Um1[flux]=np[l].bs->U[flux];
            }
          }
        end3DL
      end2DL
    end1DL

    for (gl->subiter_pc=0; gl->subiter_pc<gl->numsubiter_pc; gl->subiter_pc++){
      update_bdry_nodes(np,gl,gl->domain); 
      find_residual(np, gl, gl->domain);
      find_ximax(np, gl, gl->domain, IJK_UPDATE_YES);
#ifdef DISTMPI
      effiter_U_mem=gl->effiter_U;
#endif
      update_U_predictor_corrector(np, gl, gl->domain);
      update_linked_nodes(np, gl, TYPELEVEL_FLUID);
#ifdef DISTMPI
      exchange_U(np, gl);
      delta_effiter_U=gl->effiter_U-effiter_U_mem;
      MPI_Allreduce(&delta_effiter_U, &delta_effiter_U_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      gl->effiter_U=effiter_U_mem+delta_effiter_U_max;
#endif
    }
    gl->time+=gl->dt;
    update_runtime_codex_xi_from_gl(gl,codex);
    add_double_to_codex(codex,"time",gl->time);  
    codex->ACTIONPROCESSED=TRUE;
  }


}




void init_cycle(char *argum, SOAP_codex_t *codexcontrol){
  np_t **np;
  gl_t *gl;

  np=((readcontrolarg_t *)codexcontrol->action_args)->np;
  gl=((readcontrolarg_t *)codexcontrol->action_args)->gl;

  resume_nodes_only_in_zone(*np, gl, _zone_intersection(gl->domain_all,_zone_expansion(gl->domain,+max(hbw_res_fluid,hbw_bdry_fluid))));
  SOAP_copy_codex(codexcontrol,&(gl->cycle.codex));
  update_bdry_nodes(*np, gl, gl->domain);
  update_linked_nodes(*np, gl, TYPELEVEL_FLUID);
  find_residual(*np, gl, gl->domain);
  find_ximax(*np,gl,gl->domain,IJK_UPDATE_YES); 
  gl->RESIDUAL_ALTERED=FALSE;
}




void perform_one_iteration(np_t *np, gl_t *gl){


  process_code_runtime(np,gl,gl->cycle.code_runtime,&(gl->cycle.codex));


}

