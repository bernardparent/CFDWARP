// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2016-2018 Bernard Parent

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


#include <cycle/ressource/_ressource.h>
#include <cycle/share/cycle_share.h>




void add_Sstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  long l;
  long flux;
  flux_t S;
  if (theta==0) {
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      find_Sstar(np,gl,l,S);
#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]-=(1.0-gl->cycle.restime.weightm1_trapezoidal_default)*S[flux];
      for (flux=0; flux<nf; flux++) np[l].bs->trapezoidalm1_next[flux]-=gl->cycle.restime.weightm1_trapezoidal_default*S[flux];
#else
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]-=S[flux];
#endif
    }
  }
}


