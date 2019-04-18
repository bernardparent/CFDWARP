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


