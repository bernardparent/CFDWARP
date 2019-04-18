#include <cycle/ressource/_ressource.h>
#include <cycle/share/cycle_share.h>

#define EIGENVALUES_ENFORCED_POSITIVE TRUE

#define numiter 5


void add_Sstar_residual(long theta, long ls, long le, np_t *np, gl_t *gl){
  long l;
  long flux,cnt,row,col;
  metrics_t metrics;
  jacvars_t jacvarsp0,jacvarsm1;
  flux_t S,LUp0,LSp0,LUm1,fluxtmp,fluxtmp2,fluxtmp3,Splus,Sminus;
  sqmat_t Lambdaminus,Lambdaplus,Yminus,Yplus,Zminus,Zplus,Linvm1,Linvp0,Lp0,Lm1;

  if (theta==0) {
    for (l=ls; l!=_l_plus_one(le,gl,theta); l=_l_plus_one(l,gl,theta)){
      find_Sstar(np,gl,l,S);

      find_metrics_at_node(np, gl, l, theta, &metrics);
      find_jacvars(np[l], gl, metrics, theta, &jacvarsp0);
#ifdef UNSTEADY
      find_jacvars_from_U(np[l].bs->Um1, metrics, gl, theta, &jacvarsm1);
#else
      jacvarsm1=jacvarsp0;
#endif


      /* need to find LSp0,  LUp0, LUm1, Linvm1, Linvp0, Lp0, Lm1 */

      find_L_from_jacvars(jacvarsp0, metrics, Lp0);
      find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
      find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0);

      find_L_from_jacvars(jacvarsm1, metrics, Lm1);
      find_Linv_from_jacvars(jacvarsm1, metrics, Linvm1);
      find_LUstar_from_jacvars(jacvarsm1, metrics, LUm1);
      
      multiply_matrix_and_vector(Lp0,S,LSp0);

      for (row=0; row<nf; row++){
        for (col=0; col<nf; col++){
          Lambdaminus[row][col]=0.0;
          Lambdaplus[row][col]=0.0;
          Yminus[row][col]=0.0;
          Yplus[row][col]=0.0;
          Zminus[row][col]=0.0;
          Zplus[row][col]=0.0;
        }
      }

      for (flux=0; flux<nf; flux++) {
        Yminus[flux][flux]=0.0;
        Yplus[flux][flux]=0.0;
        Zplus[flux][flux]=max(0.0,
          +LSp0[flux]/notzero(LUp0[flux],1e-99) );
        Zminus[flux][flux]=min(0.0,
          +LSp0[flux]/notzero(LUp0[flux],1e-99) );
      }

  
      for (cnt=0; cnt<numiter; cnt++){

        multiply_matrix_and_vector(Zplus,LUp0,fluxtmp2);
        multiply_matrix_and_vector(Linvp0,fluxtmp2,fluxtmp);
        multiply_matrix_and_vector(Lm1,fluxtmp,fluxtmp2);

        multiply_matrix_and_vector(Yminus,LUm1,fluxtmp3);
        multiply_matrix_and_vector(Linvm1,fluxtmp3,fluxtmp);
        multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp3);

        for (flux=0; flux<nf; flux++) {
          Yminus[flux][flux]=min(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUm1[flux],1e-99));
          Zplus[flux][flux]=max(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp0[flux],1e-99));
        }
        for (flux=0; flux<nf; flux++) {
          Yplus[flux][flux]=max(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUm1[flux],1e-99));
          Zminus[flux][flux]=min(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp0[flux],1e-99));
        }
    
      }

      if (EIGENVALUES_ENFORCED_POSITIVE ){
        for (flux=0; flux<nf; flux++) {
          Lambdaplus[flux][flux]=Yplus[flux][flux]+Zplus[flux][flux];
          Lambdaminus[flux][flux]=Zminus[flux][flux]+Yminus[flux][flux];
        }
      } else {
        for (flux=0; flux<nf; flux++) {
          Lambdaplus[flux][flux]=Yplus[flux][flux]+Yminus[flux][flux];
          Lambdaminus[flux][flux]=Zplus[flux][flux]+Zminus[flux][flux];
        }
      }

      multiply_diagonal_matrix_and_vector(Lambdaminus,LUp0,fluxtmp);
      multiply_matrix_and_vector(Linvp0,fluxtmp,Sminus);

      multiply_diagonal_matrix_and_vector(Lambdaplus,LUm1,fluxtmp);
      multiply_matrix_and_vector(Linvm1,fluxtmp,Splus);

#ifdef _RESSOURCE_LAMBDAMINUS_STORAGE
      for (flux=0; flux<nf; flux++) np[l].bs->Lambda_S_minus[flux]=Lambdaminus[flux][flux];
#endif
#ifdef _RESTIME_TRAPEZOIDAL_RESIDUAL
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]-=(1.0-gl->cycle.restime.weightm1_trapezoidal_default)*(Splus[flux]+Sminus[flux]);
      for (flux=0; flux<nf; flux++) np[l].bs->trapezoidalm1_next[flux]-=gl->cycle.restime.weightm1_trapezoidal_default*(Splus[flux]+Sminus[flux]);
#else
      for (flux=0; flux<nf; flux++) np[l].wk->Res[flux]-=(Splus[flux]+Sminus[flux]);
#endif

    }
  }
}


