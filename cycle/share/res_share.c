// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2015-2018 Bernard Parent

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

#include <cycle/share/res_share.h>
#include <cycle/restime/_restime.h>
#include <src/bdry.h>
#include <model/_model.h>


#define CWENO_CENTRAL_WEIGHT 100.0
#define CWENO_IS_EXP 2
#define AOWENO_IS_EXP 2
#define WENO_EPSILON 1e-18




void find_musclvarscycle_offset(np_t *np, gl_t *gl, long l, long theta, long offset, musclvarscycle_t musclvars){
  long cnt;
#ifdef DISTMPI
  long flux;
  long lbdry;
  //long llink;
#else
  long llink,lbdry;
#endif


  if (!is_node_valid(np[_al(gl,l,theta,+0)],TYPELEVEL_FLUID_WORK)) fatal_error("must start with a valid node in find_musclvarscycle_offset");
  cnt=0;
  do {
    if (offset>0) cnt++; else cnt--;
    assert(is_node_in_domain_lim(l, gl));
  } while(is_node_valid(np[_al(gl,l,theta,cnt)],TYPELEVEL_FLUID_WORK) && labs(cnt)<=labs(offset));
  if (offset>0) cnt=min(offset,cnt-1); else cnt=max(cnt+1,offset);
  if (!is_node_valid(np[_al(gl,l,theta,cnt)],TYPELEVEL_FLUID_WORK)) fatal_error("problem finding musclvars in find_musclvarscycle_offset");
  if (cnt!=offset && is_node_link(np[_al(gl,l,theta,cnt)],TYPELEVEL_FLUID_WORK) ){
    /* need to find musclvars labs(offset-cnt) nodes away from linked node */
#ifdef DISTMPI
    lbdry=_al(gl,l,theta,cnt);
    assert(is_node_bdry(np[lbdry],TYPELEVEL_FLUID_WORK));
    assert(is_node_link(np[lbdry],TYPELEVEL_FLUID_WORK));
    assert(np[_al(gl,l,theta,cnt)].linkmusclvars!=NULL);
    if (np[lbdry].numlinkmusclvars<=(labs(offset-cnt)-1)){
      fatal_error("Problem in find_musclvarscycle_offset: numlinkmusclvars=%d must be greater to (labs(offset-cnt)-1)=%ld",np[lbdry].numlinkmusclvars,(labs(offset-cnt)-1));
    }
    for (flux=0; flux<nmc; flux++) musclvars[flux]=np[lbdry].linkmusclvars[flux+(labs(offset-cnt)-1)*nmc]; 
    assert(is_node_in_zone(_i(lbdry,gl,0),_i(lbdry,gl,1),_i(lbdry,gl,2),gl->domain_lim));
    //llink=_node_link(np[lbdry],0,TYPELEVEL_FLUID_WORK);
    //find_musclvarscycle(np[_al_link(np,gl, llink,lbdry,labs(offset-cnt),TYPELEVEL_FLUID_WORK)], gl, musclvars);
#else
    lbdry=_al(gl,l,theta,cnt);
    assert(is_node_bdry(np[lbdry],TYPELEVEL_FLUID_WORK));
    llink=_node_link(np[lbdry],0,TYPELEVEL_FLUID_WORK);
    find_musclvarscycle(np[_al_link(np,gl, llink,lbdry, labs(offset-cnt),TYPELEVEL_FLUID_WORK)], gl, musclvars);
#endif
  } else {
    find_musclvarscycle(np[_al(gl,l,theta,cnt)], gl, musclvars);
  }
  //printf("%E ",musclvars[nf]);
}


/* basis for the Legendre polynomials for the domain x=-0.5...0.5 
   DS Balsara et al, "An efficient class of WENO schemes with adaptive order", JCP 326:780-804, 2016   
*/
double _P_Legendre(int order, double x){
  switch (order){
    case 0:
      return(1.0);
    break;
    case 1:
      return(x);
    break;
    case 2:
      return(x*x-1.0/12.0);
    break;
    case 3:
      return(x*x*x-3.0/20.0*x);
    break;
    case 4:
      return(x*x*x*x-3.0/14.0*x*x+3.0/560.0);
    break;
    case 5:
      return(x*x*x*x*x-5.0/18.0*x*x*x+5.0/336.0*x);
    break;
    case 6:
      return(x*x*x*x*x*x-15.0/44.0*x*x*x*x+5.0/176.0*x*x-5.0/14784.0);
    break;
    case 7:
      return(x*x*x*x*x*x*x-21.0/52.0*x*x*x*x*x+105.0/2288.0*x*x*x-35.0/27456.0*x);
    break;
    case 8:
      return(x*x*x*x*x*x*x*x-7.0/15.0*x*x*x*x*x*x+7.0/104.0*x*x*x*x-7.0/2288.0*x*x+7.0/329472.0);
    break;
    default:
      fatal_error("_P_Legendre can not be called with order %ld.",order);
      return(1.0);
  }
}


/* evaluare u at location x with x a non-dimensional distance varying within the cell from -0.5 to 0.5 */
double _f_CWENO3(double um2, double um1, double up0, double up1, double up2, double x){
  double S1,S2,S3,ux,ux2,IS1,IS2,IS3,ret,gamma1,gamma2,gamma3,eps,w1,w2,w3,wtil1,wtil2,wtil3;

  ux=-2.0*um1+um2/2.0+3.0*up0/2.0;
  ux2=(um2-2.0*um1+up0)/2.0;
  S1=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x);
  IS1=sqr(ux)+13.0/3.0*sqr(ux2);

  ux=(up1-um1)/2.0;
  ux2=(um1-2.0*up0+up1)/2.0;
  S2=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x);
  IS2=sqr(ux)+13.0/3.0*sqr(ux2);

  ux=-3.0*up0/2.0+2.0*up1-up2/2.0;
  ux2=(up0-2.0*up1+up2)/2.0;
  S3=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x);
  IS3=sqr(ux)+13.0/3.0*sqr(ux2);

  gamma1=1.0;
  gamma2=CWENO_CENTRAL_WEIGHT;
  gamma3=1.0;

  eps=WENO_EPSILON;

  wtil1=gamma1/powint(eps+IS1,CWENO_IS_EXP);
  wtil2=gamma2/powint(eps+IS2,CWENO_IS_EXP);
  wtil3=gamma3/powint(eps+IS3,CWENO_IS_EXP);
  w1=wtil1/(wtil1+wtil2+wtil3);
  w2=wtil2/(wtil1+wtil2+wtil3);
  w3=wtil3/(wtil1+wtil2+wtil3);
  ret=w1*S1+w2*S2+w3*S3;

  return(ret);
}


/* evaluare u at location x with x a non-dimensional distance varying within the cell from -0.5 to 0.5 */
double _f_CWENO5(double um4, double um3, double um2, double um1, double up0, double up1, 
                 double up2, double up3, double up4, double x){
  double S1,S2,S3,S4,S5,ux,ux2,ux3,ux4,IS1,IS2,IS3,IS4,IS5,ret,
         gamma1,gamma2,gamma3,gamma4,gamma5,eps,w1,w2,w3,w4,w5,wtil1,wtil2,wtil3,wtil4,wtil5;

  ux = (-462.0*um1 + 336.0*um2 - 146.0*um3 + 27.0*um4 + 245.0*up0)/120.0;
  ux2 = (-240.0*um1 + 262.0*um2 - 128.0*um3 + 25.0*um4 + 81.0*up0)/56.0;
  ux3 = (-18.0*um1 + 24.0*um2 - 14.0*um3 + 3.0*um4 + 5.0*up0)/12.0;
  ux4 = (-4.0*um1 + 6.0*um2 - 4.0*um3 + um4 + up0)/24.0;
  S1=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x);
  IS1=sqr(ux+ux3/10.0)+13.0/3.0*sqr(ux2+123.0/455.0*ux4)+781.0/20.0*sqr(ux3)+1421461.0/2275.0*sqr(ux4);

  ux = (-192.0*um1 + 66.0*um2 - 11.0*um3 + 110.0*up0 + 27.0*up1)/120.0;
  ux2 = (10.0*um1 + 12.0*um2 - 3.0*um3 - 44.0*up0 + 25.0*up1)/56.0;
  ux3 = (12.0*um1 - 6.0*um2 + um3 - 10.0*up0 + 3.0*up1)/12.0;
  ux4 = (6.0*um1 - 4.0*um2 + um3 - 4.0*up0 + up1)/24.0;
  S2=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x);
  IS2=sqr(ux+ux3/10.0)+13.0/3.0*sqr(ux2+123.0/455.0*ux4)+781.0/20.0*sqr(ux3)+1421461.0/2275.0*sqr(ux4);
  
  ux = (-82.0*um1 + 11.0*um2 + 82.0*up1 - 11.0*up2)/120.0;
  ux2 = (40.0*um1 - 3.0*um2 - 74.0*up0 + 40.0*up1 - 3.0*up2)/56.0;
  ux3 = (2.0*um1 - um2 - 2.0*up1 + up2)/12.0;
  ux4 = (-4.0*um1 + um2 + 6.0*up0 - 4.0*up1 + up2)/24.0;
  S3=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x);
  IS3=sqr(ux+ux3/10.0)+13.0/3.0*sqr(ux2+123.0/455.0*ux4)+781.0/20.0*sqr(ux3)+1421461.0/2275.0*sqr(ux4);

  ux = (-27.0*um1 - 110.0*up0 + 192.0*up1 - 66.0*up2 + 11.0*up3)/120.0;
  ux2 = (25.0*um1 - 44.0*up0 + 10.0*up1 + 12.0*up2 - 3.0*up3)/56.0;
  ux3 = (-3.0*um1 + 10.0*up0 - 12.0*up1 + 6.0*up2 - up3)/12.0;
  ux4 = (um1 - 4.0*up0 + 6.0*up1 - 4.0*up2 + up3)/24.0;
  S4=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x);
  IS4=sqr(ux+ux3/10.0)+13.0/3.0*sqr(ux2+123.0/455.0*ux4)+781.0/20.0*sqr(ux3)+1421461.0/2275.0*sqr(ux4);

  ux = (-245.0*up0 + 462.0*up1 - 336.0*up2 + 146.0*up3 - 27.0*up4)/120.0;
  ux2 = (81.0*up0 - 240.0*up1 + 262.0*up2 - 128.0*up3 + 25.0*up4)/56.0;
  ux3 = (-5.0*up0 + 18.0*up1 - 24.0*up2 + 14.0*up3 - 3.0*up4)/12.0;
  ux4 = (up0 - 4.0*up1 + 6.0*up2 - 4.0*up3 + up4)/24.0;
  S5=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x);
  IS5=sqr(ux+ux3/10.0)+13.0/3.0*sqr(ux2+123.0/455.0*ux4)+781.0/20.0*sqr(ux3)+1421461.0/2275.0*sqr(ux4);


  gamma1=1.0;
  gamma2=1.0;
  gamma3=CWENO_CENTRAL_WEIGHT;
  gamma4=1.0;
  gamma5=1.0;

  eps=WENO_EPSILON;

  wtil1=gamma1/powint(eps+IS1,CWENO_IS_EXP);
  wtil2=gamma2/powint(eps+IS2,CWENO_IS_EXP);
  wtil3=gamma3/powint(eps+IS3,CWENO_IS_EXP);
  wtil4=gamma4/powint(eps+IS4,CWENO_IS_EXP);
  wtil5=gamma5/powint(eps+IS5,CWENO_IS_EXP);
  w1=wtil1/(wtil1+wtil2+wtil3+wtil4+wtil5);
  w2=wtil2/(wtil1+wtil2+wtil3+wtil4+wtil5);
  w3=wtil3/(wtil1+wtil2+wtil3+wtil4+wtil5);
  w4=wtil4/(wtil1+wtil2+wtil3+wtil4+wtil5);
  w5=wtil5/(wtil1+wtil2+wtil3+wtil4+wtil5);
  ret=w1*S1+w2*S2+w3*S3+w4*S4+w5*S5;

  return(ret);
}


/* evaluare u at location x with x a non-dimensional distance varying within the cell from -0.5 to 0.5 */
double _f_CWENO6(double um5, double um4, double um3, double um2, double um1, double up0, double up1, 
                    double up2, double up3, double up4, double up5, double x){
  double S1,S2,S3,S4,S5,S6,ux,ux2,ux3,ux4,ux5,IS1,IS2,IS3,IS4,IS5,IS6,ret,
         gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,eps,w1,w2,w3,w4,w5,w6,wtil1,wtil2,wtil3,wtil4,wtil5,wtil6;

  ux = (-23719.0*um1 + 22742.0*um2 - 14762.0*um3 + 5449.0*um4 - 863.0*um5 + 11153.0*up0)/5040.0;
  ux2 = (-350.0*um1 + 482.0*um2 - 348.0*um3 + 135.0*um4 - 22.0*um5 + 103.0*up0)/56.0;
  ux3 = (-317.0*um1 + 526.0*um2 - 436.0*um3 + 182.0*um4 - 31.0*um5 + 76.0*up0)/108.0;
  ux4 = (-14.0*um1 + 26.0*um2 - 24.0*um3 + 11.0*um4 - 2.0*um5 + 3.0*up0)/24.0;
  ux5 = (-5.0*um1 + 10.0*um2 - 10.0*um3 + 5.0*um4 - um5 + up0)/120.0;
  S1=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x)+ux5*_P_Legendre(5,x);
  IS1=sqr(ux+ux3/10.0+ux5/126.0)
     +13.0/3.0*sqr(ux2+123.0/455.0*ux4)
     +781.0/20.0*sqr(ux3+26045.0/49203.0*ux5)
     +1421461.0/2275.0*sqr(ux4)
     +21520059541.0/1377684.0*sqr(ux5);

  ux = (-10774.0*um1 + 5482.0*um2 - 1817.0*um3 + 271.0*um4 + 5975.0*up0 + 863.0*up1)/5040.0;
  ux2 = (-20.0*um1 + 42.0*um2 - 18.0*um3 + 3.0*um4 - 29.0*up0 + 22.0*up1)/56.0;
  ux3 = (148.0*um1 - 94.0*um2 + 29.0*um3 - 4.0*um4 - 110.0*up0 + 31.0*up1)/108.0;
  ux4 = (16.0*um1 - 14.0*um2 + 6.0*um3 - um4 - 9.0*up0 + 2.0*up1)/24.0;
  ux5 = (10.0*um1 - 10.0*um2 + 5.0*um3 - um4 - 5.0*up0 + up1)/120.0;
  S2=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x)+ux5*_P_Legendre(5,x);
  IS2=sqr(ux+ux3/10.0+ux5/126.0)
     +13.0/3.0*sqr(ux2+123.0/455.0*ux4)
     +781.0/20.0*sqr(ux3+26045.0/49203.0*ux5)
     +1421461.0/2275.0*sqr(ux4)
     +21520059541.0/1377684.0*sqr(ux5);

  ux = (-5354.0*um1 + 1417.0*um2 - 191.0*um3 + 1910.0*up0 + 2489.0*up1 - 271.0*up2)/5040.0;
  ux2 = (40.0*um1 - 3.0*um2 - 74.0*up0 + 40.0*up1 - 3.0*up2)/56.0;
  ux3 = (68.0*um1 - 34.0*um2 + 5.0*um3 - 50.0*up0 + 7.0*up1 + 4.0*up2)/108.0;
  ux4 = (-4.0*um1 + um2 + 6.0*up0 - 4.0*up1 + up2)/24.0;
  ux5 = (-10.0*um1 + 5.0*um2 - um3 + 10.0*up0 - 5.0*up1 + up2)/120.0;
  S3=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x)+ux5*_P_Legendre(5,x);
  IS3=sqr(ux+ux3/10.0+ux5/126.0)
     +13.0/3.0*sqr(ux2+123.0/455.0*ux4)
     +781.0/20.0*sqr(ux3+26045.0/49203.0*ux5)
     +1421461.0/2275.0*sqr(ux4)
     +21520059541.0/1377684.0*sqr(ux5);

  ux = (-2489.0*um1 + 271.0*um2 - 1910.0*up0 + 5354.0*up1 - 1417.0*up2 + 191.0*up3)/5040.0;
  ux2 = (40.0*um1 - 3.0*um2 - 74.0*up0 + 40.0*up1 - 3.0*up2)/56.0;
  ux3 = (-7.0*um1 - 4.0*um2 + 50.0*up0 - 68.0*up1 + 34.0*up2 - 5.0*up3)/108.0;
  ux4 = (-4.0*um1 + um2 + 6.0*up0 - 4.0*up1 + up2)/24.0;
  ux5 = (5.0*um1 - um2 - 10.0*up0 + 10.0*up1 - 5.0*up2 + up3)/120.0;
  S4=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x)+ux5*_P_Legendre(5,x);
  IS4=sqr(ux+ux3/10.0+ux5/126.0)
     +13.0/3.0*sqr(ux2+123.0/455.0*ux4)
     +781.0/20.0*sqr(ux3+26045.0/49203.0*ux5)
     +1421461.0/2275.0*sqr(ux4)
     +21520059541.0/1377684.0*sqr(ux5);

  ux = (-863.0*um1 - 5975.0*up0 + 10774.0*up1 - 5482.0*up2 + 1817.0*up3 - 271.0*up4)/5040.0;
  ux2 = (22.0*um1 - 29.0*up0 - 20.0*up1 + 42.0*up2 - 18.0*up3 + 3.0*up4)/56.0;
  ux3 = (-31.0*um1 + 110.0*up0 - 148.0*up1 + 94.0*up2 - 29.0*up3 + 4.0*up4)/108.0;
  ux4 = (2.0*um1 - 9.0*up0 + 16.0*up1 - 14.0*up2 + 6.0*up3 - up4)/24.0;
  ux5 = ( -um1 + 5.0*up0 - 10.0*up1 + 10.0*up2 - 5.0*up3 + up4)/120.0;
  S5=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x)+ux5*_P_Legendre(5,x);
  IS5=sqr(ux+ux3/10.0+ux5/126.0)
     +13.0/3.0*sqr(ux2+123.0/455.0*ux4)
     +781.0/20.0*sqr(ux3+26045.0/49203.0*ux5)
     +1421461.0/2275.0*sqr(ux4)
     +21520059541.0/1377684.0*sqr(ux5);

  ux = (-11153.0*up0 + 23719.0*up1 - 22742.0*up2 + 14762.0*up3 - 5449.0*up4 + 863.0*up5)/5040.0;
  ux2 = (103.0*up0 - 350.0*up1 + 482.0*up2 - 348.0*up3 + 135.0*up4 - 22.0*up5)/56.0;
  ux3 = (-76.0*up0 + 317.0*up1 - 526.0*up2 + 436.0*up3 - 182.0*up4 + 31.0*up5)/108.0;
  ux4 = (3.0*up0 - 14.0*up1 + 26.0*up2 - 24.0*up3 + 11.0*up4 - 2.0*up5)/24.0;
  ux5 = ( -up0 + 5.0*up1 - 10.0*up2 + 10.0*up3 - 5.0*up4 + up5)/120.0;
  S6=up0+ux*_P_Legendre(1,x)+ux2*_P_Legendre(2,x)+ux3*_P_Legendre(3,x)+ux4*_P_Legendre(4,x)+ux5*_P_Legendre(5,x);
  IS6=sqr(ux+ux3/10.0+ux5/126.0)
     +13.0/3.0*sqr(ux2+123.0/455.0*ux4)
     +781.0/20.0*sqr(ux3+26045.0/49203.0*ux5)
     +1421461.0/2275.0*sqr(ux4)
     +21520059541.0/1377684.0*sqr(ux5);

  gamma1=1.0;
  gamma2=1.0;
  gamma3=CWENO_CENTRAL_WEIGHT/2.0;
  gamma4=CWENO_CENTRAL_WEIGHT/2.0;
  gamma5=1.0;
  gamma6=1.0;

  eps=WENO_EPSILON;

  wtil1=gamma1/powint(eps+IS1,CWENO_IS_EXP);
  wtil2=gamma2/powint(eps+IS2,CWENO_IS_EXP);
  wtil3=gamma3/powint(eps+IS3,CWENO_IS_EXP);
  wtil4=gamma4/powint(eps+IS4,CWENO_IS_EXP);
  wtil5=gamma5/powint(eps+IS5,CWENO_IS_EXP);
  wtil6=gamma6/powint(eps+IS6,CWENO_IS_EXP);
  w1=wtil1/(wtil1+wtil2+wtil3+wtil4+wtil5+wtil6);
  w2=wtil2/(wtil1+wtil2+wtil3+wtil4+wtil5+wtil6);
  w3=wtil3/(wtil1+wtil2+wtil3+wtil4+wtil5+wtil6);
  w4=wtil4/(wtil1+wtil2+wtil3+wtil4+wtil5+wtil6);
  w5=wtil5/(wtil1+wtil2+wtil3+wtil4+wtil5+wtil6);
  w6=wtil6/(wtil1+wtil2+wtil3+wtil4+wtil5+wtil6);

  ret=w1*S1+w2*S2+w3*S3+w4*S4+w5*S5+w6*S6;

  return(ret);
}




void find_central_polynomial_CWENO5(double um2, double um1, double up0, double up1, double up2, double x, double *IS53, double *S53){
  double ux53_1,ux53_2,ux53_3,ux53_4;
  ux53_1 = (-82.0*um1 + 11.0*um2 + 82.0*up1 - 11.0*up2)/120.0;
  ux53_2 = (40.0*um1 - 3.0*um2 - 74.0*up0 + 40.0*up1 - 3.0*up2)/56.0;
  ux53_3 = (2.0*um1 - um2 - 2.0*up1 + up2)/12.0;
  ux53_4 = (-4.0*um1 + um2 + 6.0*up0 - 4.0*up1 + up2)/24.0;
  *IS53=sqr(ux53_1+ux53_3/10.0)+13.0/3.0*sqr(ux53_2+123.0/455.0*ux53_4)+781.0/20.0*sqr(ux53_3)+1421461.0/2275.0*sqr(ux53_4);
  *S53=up0+ux53_1*_P_Legendre(1,x)+ux53_2*_P_Legendre(2,x)+ux53_3*_P_Legendre(3,x)+ux53_4*_P_Legendre(4,x);
}


void find_central_polynomial_CWENO7(double um3, double um2, double um1, double up0, double up1, double up2, double up3, double x, double *IS74, double *S74){
  double ux74_1,ux74_2,ux74_3,ux74_4,ux74_5,ux74_6;
  ux74_1 = (-7843.0*um1 + 1688.0*um2 - 191.0*um3 + 7843.0*up1 - 1688.0*up2 + 191.0*up3)/10080.0;
  ux74_2 = (8385.0*um1 - 1014.0*um2 + 79.0*um3 - 14900.0*up0 + 8385.0*up1 - 1014.0*up2 + 79.0*up3)/10080.0;
  ux74_3 = (61.0*um1 - 38.0*um2 + 5.0*um3 - 61.0*up1 + 38.0*up2 - 5.0*up3)/216.0;
  ux74_4 = (-459.0*um1 + 144.0*um2 - 13.0*um3 + 656.0*up0 - 459.0*up1 + 144.0*up2 - 13.0*up3)/1584.0;
  ux74_5 = (-5.0*um1 + 4.0*um2 - um3 + 5.0*up1 - 4.0*up2 + up3)/240.0;
  ux74_6 = (15.0*um1 - 6.0*um2 + um3 - 20.0*up0 + 15.0*up1 - 6.0*up2 + up3)/720.0;
  *IS74 = sqr(ux74_1 + ux74_3/10.0 + ux74_5/126.0)
       + 13.0/3.0 * sqr(ux74_2 + 123.0/455.0*ux74_4 +85.0/2002.0*ux74_6) 
       + 781.0/20.0 * sqr(ux74_3+26045.0/49203.0*ux74_5)
       + 1421461.0/2275.0 * sqr(ux74_4+ 81596225.0/93816426.0*ux74_6)
       + 21520059541.0/1377684.0*sqr(ux74_5)
       + 15510384942580921.0/27582029244.0*sqr(ux74_6);
  *S74=up0+ux74_1*_P_Legendre(1,x)+ux74_2*_P_Legendre(2,x)+ux74_3*_P_Legendre(3,x)+ux74_4*_P_Legendre(4,x)
         +ux74_5*_P_Legendre(5,x)+ux74_6*_P_Legendre(6,x);
}


void find_central_polynomial_CWENO9(double um4, double um3, double um2, double um1, double up0, double up1, double up2, double up3, double up4, double x, double *IS95, double *S95){
  double ux95_1,ux95_2,ux95_3,ux95_4,ux95_5,ux95_6,ux95_7,ux95_8;

  ux95_1=(-505538.0*um1 + 136238.0*um2 - 26442.0*um3 + 2497.0*um4 + 505538.0*up1 - 136238.0*up2 + 26442.0*up3 - 2497.0*up4)/604800.0;
  ux95_2=(1205324.0*um1 - 183100.0*um2 + 24500.0*um3 - 1759.0*um4 - 2089930.0*up0 + 1205324.0*up1 - 183100.0*up2 + 24500.0*up3 - 1759.0*up4)/1330560.0;
  ux95_3=(34414.0*um1 - 24294.0*um2 + 5446.0*um3 - 541.0*um4 - 34414.0*up1 + 24294.0*up2 - 5446.0*up3 + 541.0*up4)/95040.0;
  ux95_4=(-186496.0*um1 + 66572.0*um2 - 10240.0*um3 + 773.0*um4 + 258782.0*up0 - 186496.0*up1 + 66572.0*up2 - 10240.0*up3 + 773.0*up4)/494208.0;
  ux95_5=(-526.0*um1 + 474.0*um2 - 166.0*um3 + 19.0*um4 + 526.0*up1 - 474.0*up2 + 166.0*up3 - 19.0*up4)/12480.0;
  ux95_6=(1852.0*um1 - 836.0*um2 + 196.0*um3 - 17.0*um4 - 2390.0*up0 + 1852.0*up1 - 836.0*up2 + 196.0*up3 - 17.0*up4)/43200.0;
  ux95_7=(14.0*um1 - 14.0*um2 + 6.0*um3 - um4 - 14.0*up1 + 14.0*up2 - 6.0*up3 + up4)/10080.0;
  ux95_8=(-56.0*um1 + 28.0*um2 - 8.0*um3 + um4 + 70.0*up0 - 56.0*up1 + 28.0*up2 - 8.0*up3 + up4)/40320.0;
  *IS95 = sqr(ux95_1 + ux95_3/10.0 + ux95_5/126.0 + ux95_7/1716.0)
          +13.0/3.0*sqr(ux95_2 + 123.0/455.0*ux95_4 + 85.0/2002.0*ux95_6 + 29.0/5577.0*ux95_8)
          +781.0/20.0*sqr(ux95_3 + 26045.0/49203.0*ux95_5 + 8395.0/60918.0*ux95_7)
          +1421461.0/2275.0*sqr(ux95_4 + 81596225.0/93816426.0*ux95_6 + 618438835.0/1829420307.0*ux95_8)
          +21520059541.0/1377684.0*sqr(ux95_5 + 722379670131.0/559521548066.0*ux95_7)
          +15510384942580921.0/27582029244.0*sqr(ux95_6 + 5423630339859998294.0/3024525063803279595.0*ux95_8)
          +12210527897166191835083.0/443141066068272.0*sqr(ux95_7)
          +75509368098103789336083731407561.0/42818201328263029226415.0*sqr(ux95_8);
  *S95=up0+ux95_1*_P_Legendre(1,x)+ux95_2*_P_Legendre(2,x)+ux95_3*_P_Legendre(3,x)+ux95_4*_P_Legendre(4,x)
          +ux95_5*_P_Legendre(5,x)+ux95_6*_P_Legendre(6,x)+ux95_7*_P_Legendre(7,x)+ux95_8*_P_Legendre(8,x);
}



/* evaluare u at location x with x a non-dimensional distance varying within the cell from -0.5 to 0.5 */
double _f_AOWENO(double um2, double um1, double up0, double up1, double up2, double x, double ISho, double Sho, double gammalo, double gammahi, int AOWENO_TYPE){
  double ret,gamma31,gamma32,gamma33,gammaho,eps,who,w31,w32,w33,wtil31,wtil32,wtil33,wtilho,
         ux31_1,ux31_2,IS31,ux32_1,ux32_2,IS32,
         ux33_1,ux33_2,IS33,tau,sum,S31,S32,S33;

  ux31_1=-2.0*um1+um2/2.0+3.0*up0/2.0;
  ux31_2=(um2-2.0*um1+up0)/2.0;
  IS31=sqr(ux31_1)+13.0/3.0*sqr(ux31_2);
  S31=up0+ux31_1*_P_Legendre(1,x)+ux31_2*_P_Legendre(2,x);
  
  ux32_1=(up1-um1)/2.0;
  ux32_2=(um1-2.0*up0+up1)/2.0;
  IS32=sqr(ux32_1)+13.0/3.0*sqr(ux32_2);
  S32=up0+ux32_1*_P_Legendre(1,x)+ux32_2*_P_Legendre(2,x);

  ux33_1=-3.0*up0/2.0+2.0*up1-up2/2.0;
  ux33_2=(up0-2.0*up1+up2)/2.0;
  IS33=sqr(ux33_1)+13.0/3.0*sqr(ux33_2);
  S33=up0+ux33_1*_P_Legendre(1,x)+ux33_2*_P_Legendre(2,x);

  tau=1.0/3.0*(fabs(ISho-IS31)+fabs(ISho-IS32)+fabs(ISho-IS33));

  eps=WENO_EPSILON;

  gammaho=gammahi;

  gamma31=(1.0-gammahi)*(1.0-gammalo)/2.0;
  gamma32=(1.0-gammahi)*gammalo;
  gamma33=(1.0-gammahi)*(1.0-gammalo)/2.0;

  switch (AOWENO_TYPE){ 
    case AOWENO_TYPE_COMPRESSIVE:
      assert((ISho+eps)>0.0);
      assert((IS31+eps)>0.0);
      assert((IS32+eps)>0.0);
      assert((IS33+eps)>0.0);
      wtilho=gammaho*(1.0+powint(tau/(ISho+eps),AOWENO_IS_EXP));
      wtil31=gamma31*(1.0+powint(tau/(IS31+eps),AOWENO_IS_EXP));
      wtil32=gamma32*(1.0+powint(tau/(IS32+eps),AOWENO_IS_EXP));
      wtil33=gamma33*(1.0+powint(tau/(IS33+eps),AOWENO_IS_EXP));
    break;
    case AOWENO_TYPE_DIFFUSIVE:
      assert((ISho+eps)>0.0);
      assert((IS31+eps)>0.0);
      assert((IS32+eps)>0.0);
      assert((IS33+eps)>0.0);
      //??? note: the following weights may be varied
      wtilho=gammaho/powint(ISho+eps,AOWENO_IS_EXP);
      wtil31=gamma31/powint(IS31+eps,AOWENO_IS_EXP);
      wtil32=gamma32/powint(IS32+eps,AOWENO_IS_EXP);
      wtil33=gamma33/powint(IS33+eps,AOWENO_IS_EXP);
    break;
    default:
      wtilho=0.0;
      wtil31=0.0;
      wtil32=0.0;
      wtil33=0.0; // to avoid compiler warning
      fatal_error("AOWENO_TYPE can not be set to %d",AOWENO_TYPE);
  }

  sum=wtilho+wtil31+wtil32+wtil33;
  who=wtilho/sum;
  w31=wtil31/sum;
  w32=wtil32/sum;
  w33=wtil33/sum;

  ret=who/gammaho*(Sho-gamma31*S31-gamma32*S32-gamma33*S33)+w31*S31+w32*S32+w33*S33;
  return(ret);
}


double _f_AOWENO5(double um2, double um1, double up0, double up1, double up2, double x, double gammalo, double gammahi, int AOWENO_TYPE){
  double IS53,S53,ret;
  find_central_polynomial_CWENO5(um2, um1, up0, up1, up2, x, &IS53, &S53);
  ret=_f_AOWENO(um2, um1, up0, up1, up2, x, IS53, S53, gammalo, gammahi, AOWENO_TYPE);
  return(ret);
}


double _f_AOWENO7(double um3, double um2, double um1, double up0, double up1, double up2, double up3, double x, double gammalo, double gammahi, int AOWENO_TYPE){
  double IS74,S74,ret;
  find_central_polynomial_CWENO7(um3, um2, um1, up0, up1, up2, up3, x, &IS74, &S74);
  ret=_f_AOWENO(um2, um1, up0, up1, up2, x, IS74, S74, gammalo, gammahi, AOWENO_TYPE);
  return(ret);
}


double _f_AOWENO9(double um4, double um3, double um2, double um1, double up0, double up1, double up2, double up3, double up4, double x, double gammalo, double gammahi, int AOWENO_TYPE){
  double IS95,S95,ret;
  find_central_polynomial_CWENO9(um4, um3, um2, um1, up0, up1, up2, up3, up4, x, &IS95, &S95);
  ret=_f_AOWENO(um2, um1, up0, up1, up2, x, IS95, S95, gammalo, gammahi, AOWENO_TYPE);
  return(ret);
}



double _f_AOWENO7_old(double um3, double um2, double um1, double up0, double up1, double up2, double up3, double x, double gammalo, double gammahi, int AOWENO_TYPE){
  double IS74,S74,IS53,S53,sigma,eps,ret53,ret74,wtil74,wtil53,w74,w53,ret;
  find_central_polynomial_CWENO5(um2, um1, up0, up1, up2, x, &IS53, &S53);
  find_central_polynomial_CWENO7(um3, um2, um1, up0, up1, up2, up3, x, &IS74, &S74);
  ret53=_f_AOWENO(um2, um1, up0, up1, up2, x, IS53, S53, gammalo, gammahi, AOWENO_TYPE);
  ret74=_f_AOWENO(um2, um1, up0, up1, up2, x, IS74, S74, gammalo, gammahi, AOWENO_TYPE);
  sigma=fabs(IS74-IS53);
  eps=1.0e-49;
  wtil74=gammahi*(1.0+sigma/(IS74+eps));
  wtil53=gammahi*(1.0+sigma/(IS53+eps));
  w74=wtil74/(wtil74+wtil53);
  w53=wtil53/(wtil74+wtil53);
  ret=w74/gammahi*(ret74-(1.0-gammahi)*ret53)+w53*ret53;
  return(ret);
}



double _f_AOWENO9_old(double um4, double um3, double um2, double um1, double up0, double up1, double up2, double up3, double up4, double x, double gammalo, double gammahi, int AOWENO_TYPE){
  double IS74,S74,IS95,S95,sigma,eps,ret95,ret74,wtil74,wtil95,w74,w95,ret;
  find_central_polynomial_CWENO9(um4, um3, um2, um1, up0, up1, up2, up3, up4, x, &IS95, &S95);
  find_central_polynomial_CWENO7(um3, um2, um1, up0, up1, up2, up3, x, &IS74, &S74);
  ret95=_f_AOWENO(um2, um1, up0, up1, up2, x, IS95, S95, gammalo, gammahi, AOWENO_TYPE);
  ret74=_f_AOWENO7(um3, um2, um1, up0, up1, up2, up3, x, gammalo, gammahi, AOWENO_TYPE);
  //ret74=_f_AOWENO(um2, um1, up0, up1, up2, x, IS74, S74, gammalo, gammahi, AOWENO_TYPE);
  sigma=fabs(IS95-IS74);
  eps=1.0e-49;
  wtil95=gammahi*(1.0+sigma/(IS95+eps));
  wtil74=gammahi*(1.0+sigma/(IS74+eps));
  w95=wtil95/(wtil95+wtil74);
  w74=wtil74/(wtil95+wtil74);
  ret=w95/gammahi*(ret95-(1.0-gammahi)*ret74)+w74*ret74;
  return(ret);
}

/*double _f_AOWENO9(double um4, double um3, double um2, double um1, double up0, double up1, double up2, double up3, double up4, double x){
  double IS95,S95,IS74,S74,IS53,S53;
  find_central_polynomial_CWENO9(um4, um3, um2, um1, up0, up1, up2, up3, up4, x, &IS95, &S95);
  find_central_polynomial_CWENO7(um3, um2, um1, up0, up1, up2, up3, x, &IS74, &S74);
  find_central_polynomial_CWENO5(um2, um1, up0, up1, up2, x, &IS53, &S53);

  return(_f_AOWENO(um2, um1, up0, up1, up2, x, IS74, S74));
}*/



/* Xu-Dong Liu, Stanley Osher, Tony Chan. "Weighted Essentially Non-Oscillatory Schemes", Journal
   of Computational Physics, 1994, Vol. 115, Pages 200-212  */
double _f_WENO3(double um1, double up0, double up1){
  double ret;
  double u0,u1,IS0,IS1;
  double gamma0,gamma1,eps,wtil0,wtil1,w0,w1;

  u0=-1.0/2.0*um1+3.0/2.0*up0;
  u1=1.0/2.0*up0+1.0/2.0*up1;

  /* Indicator of Smoothness */
  IS0=sqr(up0-um1);
  IS1=sqr(up1-up0);

  /* optimal weights */
  gamma0=1.0/4.0;
  gamma1=3.0/4.0;

  eps=WENO_EPSILON;
  wtil0=gamma0/sqr(eps+IS0);
  wtil1=gamma1/sqr(eps+IS1);
  w0=wtil0/(wtil0+wtil1);
  w1=wtil1/(wtil0+wtil1);
  ret=w0*u0+w1*u1;
  return(ret);
}


/* G Jiang and CW Shu. "Efficient Implementation of Weighted ENO Schemes", Journal of Computational Physics 126:202-228, 1996*/
double _f_WENO5(double um2, double um1, double up0, double up1, double up2){
  double u0,u1,u2,IS0,IS1,IS2,ret;
  double gamma0,gamma1,gamma2,eps,wtil0,wtil1,wtil2,w0,w1,w2;
  u0=3.0/8.0*um2-5.0/4.0*um1+15.0/8.0*up0;
  u1=-1.0/8.0*um1+3.0/4.0*up0+3.0/8.0*up1;
  u2=3.0/8.0*up0+3.0/4.0*up1-1.0/8.0*up2;

  /* Indicator of Smoothness */
  IS0=1.0/3.0*(4.0*um2*um2-19.0*um2*um1+25.0*um1*um1+11.0*um2*up0-31.0*um1*up0+10.0*up0*up0);
  IS1=1.0/3.0*(4.0*um1*um1-13.0*um1*up0+13.0*up0*up0+5.0*um1*up1-13.0*up0*up1+4.0*up1*up1);
  IS2=1.0/3.0*(10.0*up0*up0-31.0*up0*up1+25.0*up1*up1+11.0*up0*up2-19.0*up1*up2+4.0*up2*up2);

  /* optimal weights */
  gamma0=1.0/16.0;
  gamma1=5.0/8.0;
  gamma2=5.0/16.0;

  eps=WENO_EPSILON;
  wtil0=gamma0/sqr(eps+IS0);
  wtil1=gamma1/sqr(eps+IS1);
  wtil2=gamma2/sqr(eps+IS2);
  w0=wtil0/(wtil0+wtil1+wtil2);
  w1=wtil1/(wtil0+wtil1+wtil2);
  w2=wtil2/(wtil0+wtil1+wtil2);
  ret=w0*u0+w1*u1+w2*u2;
  return(ret);
}

/* see D.S. Balsara and C.-W. Shu, “Monotonicity Preserving weighted essentially non-oscillatory schemes
with increasingly high order of accuracy,” J.Comput.Phys., vol. 160, pp. 405–452, 2000.
  and Yiqing Shen and Gecheng Zha, "A Robust Seventh-order WENO Scheme and Its Applications", AIAA Paper 2008-0757 */
double _f_WENO7(double um3, double um2, double um1, double up0, double up1, double up2, double up3){
  double ret;
  double u0,u1,u2,u3,IS0,IS1,IS2,IS3,w0,w1,w2,w3;
  double gamma0,gamma1,gamma2,gamma3,eps,wtil0,wtil1,wtil2,wtil3;

  u0=-1.0/4.0*um3+13.0/12.0*um2-23.0/12.0*um1+25.0/12.0*up0;
  u1=1.0/12.0*um2-5.0/12.0*um1+13.0/12.0*up0+1.0/4.0*up1;
  u2=-1.0/12.0*um1+7.0/12.0*up0+7.0/12.0*up1-1.0/12.0*up2;
  u3=1.0/4.0*up0+13.0/12.0*up1-5.0/12.0*up2+1.0/12.0*up3;

  /* Indicator of Smoothness */
  IS0=um3*(547.0*um3-3882.0*um2+4642.0*um1-1854.0*up0)
       +um2*(7043.0*um2-17246.0*um1+7042.0*up0)
       +um1*(11003.0*um1-9402.0*up0)+2107.0*up0*up0;
  IS1=um2*(267.0*um2-1642.0*um1+1602.0*up0-494.0*up1)
       +um1*(2843.0*um1-5966.0*up0+1922.0*up1)
       +up0*(3443.0*up0-2522.0*up1)+547.0*up1*up1;
  IS2=um1*(547.0*um1-2522.0*up0+1922.0*up1-494.0*up2)
       +up0*(3443.0*up0-5966.0*up1+1602.0*up2)
       +up1*(2843.0*up1-1642.0*up2)+267.0*up2*up2;
  IS3=up0*(2107.0*up0-9402.0*up1+7042.0*up2-1854.0*up3)
       +up1*(11003.0*up1-17246.0*up2+4642.0*up3)
       +up2*(7043.0*up2-3882.0*up3)+547.0*up3*up3;

  /* optimal weights */
  gamma0=1.0/35.0;
  gamma1=12.0/35.0;
  gamma2=18.0/35.0;
  gamma3=4.0/35.0;
  eps=WENO_EPSILON;

  wtil0=gamma0/sqr(eps+IS0);
  wtil1=gamma1/sqr(eps+IS1);
  wtil2=gamma2/sqr(eps+IS2);
  wtil3=gamma3/sqr(eps+IS3);

  w0=wtil0/(wtil0+wtil1+wtil2+wtil3);
  w1=wtil1/(wtil0+wtil1+wtil2+wtil3);
  w2=wtil2/(wtil0+wtil1+wtil2+wtil3);
  w3=wtil3/(wtil0+wtil1+wtil2+wtil3);

  ret=w0*u0+w1*u1+w2*u2+w3*u3;
  return(ret);
}


/* see D.S. Balsara and C.-W. Shu, “Monotonicity Preserving weighted essentially non-oscillatory schemes
with increasingly high order of accuracy,” J.Comput.Phys., vol. 160, pp. 405–452, 2000. */
double _f_WENO9(double um4, double um3, double um2, double um1, double up0, double up1, double up2, double up3, double up4){
  double ret;
  double u0,u1,u2,u3,u4,IS0,IS1,IS2,IS3,IS4,w0,w1,w2,w3,w4;
  double gamma0,gamma1,gamma2,gamma3,gamma4,eps,wtil0,wtil1,wtil2,wtil3,wtil4;

  u0= 1.0/5.0 *um4 - 21.0/20.0*um3 + 137.0/60.0*um2 -163.0/60.0*um1 + 137.0/60.0*up0;
  u1=-1.0/20.0*um3 + 17.0/60.0*um2 -  43.0/60.0*um1 + 77.0/60.0*up0 +   1.0/5.0 *up1;
  u2= 1.0/30.0*um2 - 13.0/60.0*um1 +  47.0/60.0*up0 +  9.0/20.0*up1 -   1.0/20.0*up2;
  u3=-1.0/20.0*um1 +  9.0/20.0*up0 +  47.0/60.0*up1 - 13.0/60.0*up2 +   1.0/30.0*up3;
  u4= 1.0/5.0 *up0 + 77.0/60.0*up1 -  43.0/60.0*up2 + 17.0/60.0*up3 -   1.0/20.0*up4;

  /* Indicator of Smoothness */
  IS0=um4*(  22658.0*um4- 208501.0*um3+ 364863.0*um2-288007.0*um1+86329.0*up0)
       +um3*( 482963.0*um3-1704396.0*um2+1358458.0*um1-411487.0*up0)
       +um2*(1521393.0*um2-  2462076*um1+ 758823.0*up0)
       +um1*(1020563.0*um1- 649501.0*up0)
       +107918.0*up0*up0;
  IS1=um3*(  6908.0*um3- 60871.0*um2+ 99213.0*um1-70237.0*up0+18079.0*up1)
       +um2*(138563.0*um2-464976.0*um1+337018.0*up0-88297.0*up1)
       +um1*(406293.0*um1-611976.0*up0+165153.0*up1)
       +up0*(242723.0*up0-140251.0*up1)
       +22658.0*up1*up1;
  IS2=um2*(  6908.0*um2- 51001.0*um1+ 67923.0*up0-38947.0*up1+8209.0*up2)
       +um1*(104963.0*um1-299076.0*up0+179098.0*up1-38947.0*up2)
       +up0*(231153.0*up0-299076.0*up1+ 67923.0*up2)
       +up1*(104963.0*up1- 51001.0*up2)
       +6908.0*up2*up2;
  IS3=um1*(  22658.0*um1-140251.0*up0+165153.0*up1-88297.0*up2+18079.0*up3)
       +up0*( 242723.0*up0-611976.0*up1+337018.0*up2-70237.0*up3)
       +up1*( 406293.0*up1-464976.0*up2+ 99213.0*up3)
       +up2*( 138563.0*up2- 60871.0*up3)
       +6908.0*up3*up3;
  IS4=up0*( 107918.0*up0- 649501.0*up1+ 758823.0*up2-411487.0*up3+86329.0*up4)
       +up1*(1020563.0*up1-2462076.0*up2+1358458.0*up3-288007.0*up4)
       +up2*(1521393.0*up2-1704396.0*up3+ 364863.0*up4)
       +up3*( 482963.0*up3- 208501.0*up4)
       +22658.0*up4*up4;

  /* optimal weights */
  gamma0=1.0/126.0;
  gamma1=10.0/63.0;
  gamma2=10.0/21.0;
  gamma3=20.0/63.0;
  gamma4=5.0/126.0;
  eps=WENO_EPSILON;

  wtil0=gamma0/sqr(eps+IS0);
  wtil1=gamma1/sqr(eps+IS1);
  wtil2=gamma2/sqr(eps+IS2);
  wtil3=gamma3/sqr(eps+IS3);
  wtil4=gamma4/sqr(eps+IS4);

  w0=wtil0/(wtil0+wtil1+wtil2+wtil3+wtil4);
  w1=wtil1/(wtil0+wtil1+wtil2+wtil3+wtil4);
  w2=wtil2/(wtil0+wtil1+wtil2+wtil3+wtil4);
  w3=wtil3/(wtil0+wtil1+wtil2+wtil3+wtil4);
  w4=wtil4/(wtil0+wtil1+wtil2+wtil3+wtil4);

  ret=w0*u0+w1*u1+w2*u2+w3*u3+w4*u4;
  return(ret);
}

/* dudesired is the additional u to add to the first order u
   duinterface is the change in u at the interface */
static double _f_positive_coefficients(double u1o, double dudesired, double duinterface){
  double ret;
  if (duinterface>0.0){
    ret=u1o+min(duinterface,dudesired);
  } else {
    ret=u1o+max(duinterface,dudesired);
  } 
  return(ret);
}


static double _f_TVD(double u_desired, double  u_1o, double deltau, double deltau_upw){
  double ret;
  ret=_f_positive_coefficients(u_1o,u_desired-u_1o,deltau);
  if (sign(deltau)!=sign(deltau_upw)) ret=u_1o;
  return(ret);
}


double _f_TVD2(double um1, double up0, double up1, int LIMITER){
  double ret;
  ret=up0+0.5*(up0-um1)*_limiter_TVD((up1-up0)/notzero(up0-um1,1e-99),LIMITER);
  return(ret);
}



double _f_TVD3(double um1, double up0, double up1){
  double ret,u3o;
  u3o=-1.0/8.0*um1+3.0/4.0*up0+3.0/8.0*up1;
  ret=_f_TVD(u3o,up0,up1-up0,up0-um1);
  return(ret);
}


double _f_TVD4(double um2, double um1, double up0, double up1){
  double ret,u4o;
  u4o=1.0/12.0*um2-5.0/12.0*um1+13.0/12.0*up0+1.0/4.0*up1;
  ret=_f_TVD(u4o,up0,up1-up0,up0-um1);
  return(ret);
}


double _f_TVD5(double um2, double um1, double up0, double up1, double up2){
  double ret,u5o;
  u5o=1.0/30.0*um2 - 13.0/60.0*um1 +  47.0/60.0*up0 +  9.0/20.0*up1 -   1.0/20.0*up2;
  ret=_f_TVD(u5o,up0,up1-up0,up0-um1);
  return(ret);
}


double _limiter_TVD(double r, int LIMITER){
	double phizero,beta;

   phizero=0.0;


   switch (LIMITER){

     /* first order accurate */
     case LIMITER_FIRSTORDER:
		   phizero=0.0;
     break;

     /* standard minmod limiter, Roe 1986 */
     case LIMITER_MINMOD:
		   phizero=max(0.0,min(1.0,r));
     break;

     /* most compressive minmod-type limiter */
     case LIMITER_MINMOD2:
		   phizero=max(0.0,min(1.0,2.0*r));
     break;

     /* superbee limiter */
     case LIMITER_SUPERBEE:
       phizero=max(0.0,max(min(1.0,2.0*r),min(2.0,r)));
     break;

     /* SMART limiter (Gaskell & Lau, 1988)*/
     case LIMITER_SMART:
       phizero=max(0.0,min(2.0*r,min(4.0,0.25+0.75*r)));
     break;

     case LIMITER_VANLEER:
       phizero=(r+fabs(r))/(1.0+fabs(r));
     break;

     /* Venkatakrishnan limiter, 1993 */
     case LIMITER_VENKATAKRISHNAN:
       phizero=0.5*(r+1.0)*min(
                    4.0*r*(3.0*r+1.0)/(11.0*r*r+4.0*r+1.0), 
                    4.0*(r+3.0)/(r*r+4.0*r+11.0)
                   );
     break;

     /* Koren third-order limiter, Koren 1993 */
     case LIMITER_KOREN:
       phizero=max(0.0,min(2.0,min(2.0*r,(2.0+r)/3.0)));
     break;

     /* Van Albada limiter, symmetric, 1982 */
     case LIMITER_VANALBADA:
       phizero=(r*r+r)/(r*r+1.0);
     break;
     
     /* Ospre limiter, symmetric (Waterson & Deconinck, 1995) */
     case LIMITER_OSPRE:
       phizero=1.5*(r*r+r)/(r*r+r+1.0);
     break;

     /* Osher limiter, symmetric (Chatkravathy and Osher, 1983) */
     case LIMITER_OSHER:
       beta=1.5; /* 1<beta<2 */
       phizero=max(0.0,min(r,beta));
     break;

     /* Sweby limiter, symmetric (Sweby, 1984) */
     case LIMITER_SWEBY:
       beta=1.5; /* 1<beta<2 */
       phizero=max3(0.0,min(beta*r,1.0),min(r,beta));
     break;
  
     default:
       fatal_error("Limiter %d is not allowed when calling _limiter_TVD()",LIMITER);
       phizero=0.0;
 
  }

	return(phizero);
}




void find_Lambda_minus_dtau_FVS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, sqmat_t lambdaminus){
  sqmat_t lambda,lambdap;
  jacvars_t jacvars;
  metrics_t metrics;
  long row,col;

  find_metrics_at_node(np, gl, l, theta, &metrics);
  find_jacvars(np[l],gl,metrics,theta,&jacvars);

  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvars, metrics, EIGENVALCOND, lambdap);
  find_Lambda_from_jacvars(jacvars, metrics, lambda);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaminus[row][col]=0.5*(lambda[row][col]-lambdap[row][col]);
    }
  }

}

/* used to determine the pseudotime step */
void find_Lambda_plus_dtau_FVS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, sqmat_t lambdaplus){
  sqmat_t lambda,lambdap;
  jacvars_t jacvars;
  metrics_t metrics;
  long row,col;

  find_metrics_at_node(np, gl, l, theta, &metrics);
  find_jacvars(np[l],gl,metrics,theta,&jacvars);

  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvars, metrics, EIGENVALCOND, lambdap);
  find_Lambda_from_jacvars(jacvars, metrics, lambda);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaplus[row][col]=0.5*(lambda[row][col]+lambdap[row][col]);
    }
  }

}


void find_jacvars_at_interface(np_t *np, gl_t *gl, long lL, long lR, long theta, int AVERAGING, jacvars_t *jacvars){
  metrics_t metrics;
  jacvars_t jacvarsL,jacvarsR;

  find_metrics_at_interface(np,gl,lL,lR,theta,&metrics);
  find_jacvars(np[lL], gl, metrics, theta, &jacvarsL);
  find_jacvars(np[lR], gl, metrics, theta, &jacvarsR);
  find_jacvars_at_interface_from_jacvars(jacvarsL, jacvarsR, gl, theta, metrics, AVERAGING, jacvars);
}


void find_jacvars_at_interface_from_jacvars(jacvars_t jacvarsL, jacvars_t jacvarsR, gl_t *gl, long theta, metrics_t metrics, int AVERAGING, jacvars_t *jacvars){
  switch (AVERAGING){
    case AVERAGING_ROE:
      find_jacvars_at_interface_Roe_average(jacvarsL, jacvarsR, gl, theta, jacvars);
    break;
    case AVERAGING_ARITH:
      find_jacvars_at_interface_arith_average(jacvarsL, jacvarsR, gl, theta, jacvars);
    break;
    default:
      fatal_error("Interface averaging can not be set to %d in find_jacvars_at_interface_from_jacvars.",AVERAGING);
  }
}


void find_Lambda_minus_dtau_FDS_from_jacvars(gl_t *gl, jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t lambdaminus){
  sqmat_t lambda,lambdap;
  long row,col;

  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvars, metrics, EIGENVALCOND, lambdap);
  find_Lambda_from_jacvars(jacvars, metrics, lambda);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaminus[row][col]=0.5*(lambda[row][col]-lambdap[row][col]);
    }
  }

}


/* used to determine the pseudotime step */
void find_Lambda_minus_dtau_FDS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, int AVERAGING, sqmat_t lambdaminus){
  jacvars_t jacvarsm1h;
  metrics_t metrics;

  find_jacvars_at_interface(np,gl,_al(gl,l,theta,-1),_al(gl,l,theta,+0),theta,AVERAGING,&jacvarsm1h);
  find_metrics_at_interface(np, gl, _al(gl,l,theta,-1), _al(gl,l,theta,+0), theta, &metrics);

  find_Lambda_minus_dtau_FDS_from_jacvars(gl, jacvarsm1h, metrics, EIGENVALCOND, lambdaminus);
}



void find_Lambda_plus_dtau_FDS_from_jacvars(gl_t *gl, jacvars_t jacvars, metrics_t metrics, int EIGENVALCOND, sqmat_t lambdaplus){
  sqmat_t lambda,lambdap;
  long row,col;

  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvars, metrics, EIGENVALCOND, lambdap);
  find_Lambda_from_jacvars(jacvars, metrics, lambda);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaplus[row][col]=0.5*(lambda[row][col]+lambdap[row][col]);
    }
  }
}



/* used to determine the pseudotime step */
void find_Lambda_plus_dtau_FDS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, int AVERAGING, sqmat_t lambdaplus){
  jacvars_t jacvarsp1h;
  metrics_t metrics;

  find_jacvars_at_interface(np,gl,_al(gl,l,theta,+0),_al(gl,l,theta,+1),theta,AVERAGING,&jacvarsp1h);
  find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metrics);

  find_Lambda_plus_dtau_FDS_from_jacvars(gl, jacvarsp1h, metrics, EIGENVALCOND, lambdaplus);
}


/* used to determine the pseudotime step */
void find_Lambda_plus_minus_dtau_FDS(np_t *np, gl_t *gl, long l, long theta, int EIGENVALCOND, int AVERAGING, sqmat_t lambdaplus, sqmat_t lambdaminus){
   sqmat_t lambda,lambdap;
  jacvars_t jacvarsp1h;
  metrics_t metrics;
  long row,col;


  find_jacvars_at_interface(np,gl,_al(gl,l,theta,+0),_al(gl,l,theta,+1),theta,AVERAGING,&jacvarsp1h);
  find_metrics_at_interface(np, gl, _al(gl,l,theta,+0), _al(gl,l,theta,+1), theta, &metrics);


  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvarsp1h, metrics, EIGENVALCOND, lambdap);
  find_Lambda_from_jacvars(jacvarsp1h, metrics, lambda);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaplus[row][col]=0.5*(lambda[row][col]+lambdap[row][col]);
      lambdaminus[row][col]=0.5*(lambda[row][col]-lambdap[row][col]);
    }
  }
}


/* Find Lambda minus and plus for a FDS flux discretization where F_i+F_{i+1} is NOT replaced by A_{i+1/2}*(U_i+U_{i+1}) */ 
static void find_Lambda_minus_plus_FDSplus_muscl(np_t *np, gl_t *gl, long lp0, long lp1, jacvars_t jacvarsp0, jacvars_t jacvarsp1h, jacvars_t jacvarsp1, jacvars_t jacvarsp0_RE, jacvars_t jacvarsp1h_RE, jacvars_t jacvarsp1_RE, metrics_t metrics, sqmat_t lambdaminus, sqmat_t lambdaplus, long theta, long numiter, int EIGENVALCOND){
  sqmat_t Yminus,Yplus,Zplus,Zminus,Lp0,Linvp0,Linvp1,Lp1,Lp1h_RE,Linvp1h_RE,lambdapp1h_RE,lambdap1h_RE;
  long row,col,flux,cnt;
  flux_t Up0,fluxtmp,fluxtmp2,fluxtmp3,LUp0,LUp1,Up1,Up0_RE,Up1_RE,Fp0_RE,Fp1_RE,Vp0,Vp1;
#ifdef _RESTIME_CDF
  sqmat_t lambdaxt_p1h,lambdaxt_p1,lambdaxt_p0;
  flux_t deltaxt_p1h;
  flux_t LUp0_RE,LUp1_RE;
#endif


  find_L_from_jacvars(jacvarsp1h_RE, metrics, Lp1h_RE);
  find_Linv_from_jacvars(jacvarsp1h_RE, metrics, Linvp1h_RE);
  find_Lambda_from_jacvars(jacvarsp1h_RE, metrics, lambdap1h_RE);
  
  set_matrix_to_zero(lambdapp1h_RE);
  for (flux=0; flux<nf; flux++) lambdapp1h_RE[flux][flux]=fabs(lambdap1h_RE[flux][flux]);
  find_Ustar_from_jacvars(jacvarsp0_RE,  metrics, Up0_RE);
  find_Ustar_from_jacvars(jacvarsp1_RE,  metrics, Up1_RE);

  find_L_from_jacvars(jacvarsp0, metrics, Lp0);
  find_Ustar_from_jacvars(jacvarsp0,  metrics, Up0);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0);

  find_L_from_jacvars(jacvarsp1, metrics, Lp1);
  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_Linv_from_jacvars(jacvarsp1, metrics, Linvp1);
  find_Ustar_from_jacvars(jacvarsp1,  metrics, Up1);
  find_LUstar_from_jacvars(jacvarsp1, metrics, LUp1);


  find_Fstar_from_jacvars(jacvarsp0_RE,  metrics, Fp0_RE);
  find_Fstar_from_jacvars(jacvarsp1_RE,  metrics, Fp1_RE);
#ifdef _RESTIME_CDF
  find_LUstar_from_jacvars(jacvarsp0_RE, metrics, LUp0_RE);
  find_LUstar_from_jacvars(jacvarsp1_RE, metrics, LUp1_RE);
  find_Lambdaxt_interface(np, gl, lp0, lp1, theta, jacvarsp1h_RE, metrics, lambdaxt_p1h);
  find_Lambdaxt(np, gl, lp0, jacvarsp0_RE, metrics, lambdaxt_p0);
  find_Lambdaxt(np, gl, lp1, jacvarsp1_RE, metrics, lambdaxt_p1);

  find_Deltaxt_interface(np, gl, lp0, lp1, theta, jacvarsp0, jacvarsp1, metrics, deltaxt_p1h);
//  find_Deltaxt_interface(np, gl, lp0, lp1, theta, jacvarsp1h_RE, jacvarsp1h_RE, metrics, deltaxt_p1h);
  
  for (flux=0; flux<nf; flux++) lambdapp1h_RE[flux][flux]+=-fabs(lambdaxt_p1h[flux][flux])*deltaxt_p1h[flux];

#endif

  multiply_matrix_and_vector(Lp1h_RE,Up0_RE,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp2[flux]=(lambdapp1h_RE[flux][flux])*fluxtmp[flux];
  multiply_matrix_and_vector(Linvp1h_RE,fluxtmp2,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp[flux]+=Fp0_RE[flux];
  multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);

  multiply_matrix_and_vector(Lp1h_RE,Up1_RE,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp3[flux]=(-lambdapp1h_RE[flux][flux])*fluxtmp[flux];
  multiply_matrix_and_vector(Linvp1h_RE,fluxtmp3,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp[flux]+=Fp1_RE[flux];
  multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);

  for (flux=0; flux<nf; flux++) {
    Vp0[flux]=
#ifdef _RESTIME_CDF
//       -0.5*lambdaxt_p1h[flux][flux]*deltaxt_p1h[flux]
       -0.5*lambdaxt_p0[flux][flux]*deltaxt_p1h[flux]
#endif
       +0.5*fluxtmp2[flux]/notzero(LUp0[flux],1e-99);
    Vp1[flux]=
#ifdef _RESTIME_CDF
//       -0.5*lambdaxt_p1h[flux][flux]*deltaxt_p1h[flux]
       -0.5*lambdaxt_p1[flux][flux]*deltaxt_p1h[flux]
#endif
       +0.5*fluxtmp3[flux]/notzero(LUp1[flux],1e-99);
  }

#ifdef _RESCONV_INCLUDES_DIFFUSION
  for (flux=0; flux<nf; flux++) fluxtmp[flux]=np[lp0].wk->Fp1h_diffusion[theta][flux];
  multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);
  for (flux=0; flux<nf; flux++) Vp0[flux]+=fluxtmp2[flux]/LUp0[flux];
//  multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp2);
//  for (flux=0; flux<nf; flux++) Vp1[flux]+=0.5*fluxtmp2[flux]/LUp1[flux];
#endif

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaminus[row][col]=0.0;
      lambdaplus[row][col]=0.0;
      Yminus[row][col]=0.0;
      Yplus[row][col]=0.0;
      Zminus[row][col]=0.0;
      Zplus[row][col]=0.0;
    }
  }

  for (flux=0; flux<nf; flux++) {
    Yminus[flux][flux]=min(0.0, Vp0[flux] );
    Yplus[flux][flux]=max(0.0, Vp0[flux] );
    Zminus[flux][flux]=min(0.0, Vp1[flux] );
    Zplus[flux][flux]=max(0.0, Vp1[flux] );
  }
  
  for (cnt=0; cnt<numiter; cnt++){
    multiply_matrix_and_vector(Zplus,LUp1,fluxtmp2);
    multiply_matrix_and_vector(Linvp1,fluxtmp2,fluxtmp);
    multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);
    multiply_matrix_and_vector(Yminus,LUp0,fluxtmp3);
    multiply_matrix_and_vector(Linvp0,fluxtmp3,fluxtmp);
    multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);
    for (flux=0; flux<nf; flux++) {
      Yminus[flux][flux]=min(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1e-99));
      Zplus[flux][flux]=max(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1e-99));
    }
    for (flux=0; flux<nf; flux++) {
      Yplus[flux][flux]=max(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1e-99));
      Zminus[flux][flux]=min(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1e-99));
    }
  }

  for (flux=0; flux<nf; flux++) {
    lambdaplus[flux][flux]=Yplus[flux][flux]+Zplus[flux][flux];///LUp0[flux]*LUp1[flux];
    lambdaminus[flux][flux]=Zminus[flux][flux]+Yminus[flux][flux];///LUp1[flux]*LUp0[flux];
  }

  condition_Lambda_plus_minus(np, gl, lp0, theta, jacvarsp0, jacvarsp1, metrics, EIGENVALCOND, lambdaplus,lambdaminus);


}




void filter_Fstar_interface_positivity_preserving_PARENT(np_t *np, gl_t *gl, long lp0, long theta, metrics_t metrics, long numiter, int EIGENVALCOND, flux_t Fint, flux_t Fpositive, sqmat_t lambdaminus, sqmat_t lambdaplus){

  flux_t LUp0,LUp1,Up0,Up1,Vp1,Vp0;
  sqmat_t Lp0,Lp1,Linvp0,Linvp1;
  flux_t fluxtmp2,fluxtmp,fluxtmp3,Fp0,Fp1;
  sqmat_t Yminus,Yplus,Zminus,Zplus;
  long flux,cnt,lp1;
  jacvars_t jacvarsp0, jacvarsp1;

  lp1=_al(gl,lp0,theta,+1);

  find_jacvars(np[lp0],gl,metrics,theta,&jacvarsp0);
  find_jacvars(np[lp1],gl,metrics,theta,&jacvarsp1);

  find_L_from_jacvars(jacvarsp0, metrics, Lp0);
  find_L_from_jacvars(jacvarsp1, metrics, Lp1);
  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_Linv_from_jacvars(jacvarsp1, metrics, Linvp1);

  find_Ustar_from_jacvars(jacvarsp0,  metrics, Up0);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0);
  find_Ustar_from_jacvars(jacvarsp1,  metrics, Up1);
  find_LUstar_from_jacvars(jacvarsp1, metrics, LUp1);

  multiply_matrix_and_vector(Lp0,Fint,fluxtmp2);
  for (flux=0; flux<nf; flux++) Vp0[flux]=0.5*fluxtmp2[flux]/LUp0[flux];

  multiply_matrix_and_vector(Lp1,Fint,fluxtmp2);
  for (flux=0; flux<nf; flux++) Vp1[flux]=0.5*fluxtmp2[flux]/LUp1[flux];

  set_matrix_to_zero(lambdaminus);
  set_matrix_to_zero(lambdaplus);
  set_matrix_to_zero(Yminus);
  set_matrix_to_zero(Yplus);
  set_matrix_to_zero(Zminus);
  set_matrix_to_zero(Zplus);
    
  for (flux=0; flux<nf; flux++) {
    Yminus[flux][flux]=min(0.0, Vp0[flux] );
    Yplus[flux][flux]=max(0.0, Vp0[flux] );
    Zminus[flux][flux]=min(0.0, Vp1[flux] );
    Zplus[flux][flux]=max(0.0, Vp1[flux] );
  }

  for (cnt=0; cnt<numiter; cnt++){
    multiply_matrix_and_vector(Zplus,LUp1,fluxtmp2);
    multiply_matrix_and_vector(Linvp1,fluxtmp2,fluxtmp);
    multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);
    multiply_matrix_and_vector(Yminus,LUp0,fluxtmp3);
    multiply_matrix_and_vector(Linvp0,fluxtmp3,fluxtmp);
    multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);
    for (flux=0; flux<nf; flux++) {
      Yminus[flux][flux]=min(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1e-99));
      Zplus[flux][flux]=max(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1e-99));
    }
    for (flux=0; flux<nf; flux++) {
      Yplus[flux][flux]=max(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1e-99));
      Zminus[flux][flux]=min(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1e-99));
    }   
  }

  for (flux=0; flux<nf; flux++) {
    lambdaplus[flux][flux]=Yplus[flux][flux]+Zplus[flux][flux];///LUp0[flux]*LUp1[flux];
    lambdaminus[flux][flux]=Zminus[flux][flux]+Yminus[flux][flux];///LUp1[flux]*LUp0[flux];
  }
  
  condition_Lambda_plus_minus(np, gl, lp0, theta, jacvarsp0, jacvarsp1, metrics, EIGENVALCOND, lambdaplus,lambdaminus);

  multiply_diagonal_matrix_and_vector(lambdaplus,LUp0,fluxtmp2);
  multiply_matrix_and_vector(Linvp0,fluxtmp2,Fp0);

  multiply_diagonal_matrix_and_vector(lambdaminus,LUp1,fluxtmp2);
  multiply_matrix_and_vector(Linvp1,fluxtmp2,Fp1);

  for (flux=0; flux<nf; flux++) Fpositive[flux]=Fp0[flux]+Fp1[flux];
}


void filter_Fstar_interface_positivity_preserving_TEST(np_t *np, gl_t *gl, long lp0, long theta, metrics_t metrics, long numiter, int EIGENVALCOND, flux_t Fint, flux_t Fpositive, sqmat_t lambdaminus, sqmat_t lambdaplus){
  flux_t Fpositive_PARENT;
  long flux;
  filter_Fstar_interface_positivity_preserving_PARENT(np, gl, lp0, theta, metrics, numiter, EIGENVALCOND, Fint, Fpositive_PARENT, lambdaminus, lambdaplus);
  for (flux=0; flux<nf; flux++) {
    if (flux==fluxet) Fpositive[flux]=Fint[flux];
      else Fpositive[flux]=Fpositive_PARENT[flux];
  }
}


void find_Fstar_interface_FDSplus_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h, long theta, 
                     flux_t musclvarsm1h, flux_t musclvarsp1h, metrics_t metrics,  long numiter, 
                     int EIGENVALCOND, int AVERAGING, flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h){
  flux_t Fm1h,Fp1h,fluxtmp,fluxtmp2;
  /* flux_t Ustar; */
  sqmat_t R;
  /* sqmat_t L; */
  long flux;
  jacvars_t jacvarsm1h,jacvarsp1h,jacvarsp0;
  jacvars_t jacvarsm1h_RE,jacvarsp1h_RE,jacvarsp0_RE;

  find_jacvars_from_musclvars(musclvarsm1h, metrics, gl, theta, &jacvarsm1h_RE);
  find_jacvars_from_musclvars(musclvarsp1h, metrics, gl, theta, &jacvarsp1h_RE);
  
  find_jacvars_at_interface_from_jacvars(jacvarsm1h_RE, jacvarsp1h_RE, gl, theta, metrics, AVERAGING, &jacvarsp0_RE);
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_jacvars_at_interface_from_jacvars(jacvarsm1h, jacvarsp1h, gl, theta, metrics, AVERAGING, &jacvarsp0);

  find_Lambda_minus_plus_FDSplus_muscl(np,gl,lm1h,lp1h,jacvarsm1h,jacvarsp0,jacvarsp1h, jacvarsm1h_RE,jacvarsp0_RE,jacvarsp1h_RE,metrics, lambdaminusp1h, lambdaplusm1h, theta, numiter, EIGENVALCOND);

  find_Linv_from_jacvars(jacvarsm1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsm1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaplusm1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fm1h);

  find_Linv_from_jacvars(jacvarsp1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsp1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaminusp1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fp1h);

  for (flux=0; flux<nf; flux++) Fint[flux]=Fm1h[flux]+Fp1h[flux];
}


void find_Fstar_interface_FDS_muscl_with_CDF(np_t *np, gl_t *gl, long lm1h, long lp1h,  long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, int AVERAGING, bool RESTRAINED, flux_t Fint){

  sqmat_t Lm1h,Linvm1h,Linvp1h,Lp1h,Lp0_RE,Linvp0_RE,lambdapp0_RE,lambdap0_RE;
  long flux;
  flux_t Um1h,fluxtmp,fluxtmp2,fluxtmp3,LUm1h,LUp1h,Up1h,Um1h_RE,Up1h_RE,Fm1h_RE,Fp1h_RE;
#ifdef _RESTIME_CDF
  sqmat_t lambdaxt_p0,lambdaxt_p1h,lambdaxt_m1h;
  flux_t deltaxt_p0;
  flux_t LUm1h_RE,LUp1h_RE;
#endif
  flux_t Fm1h,Fp1h;
  sqmat_t R;
  jacvars_t jacvarsm1h,jacvarsp1h,jacvarsp0;
  jacvars_t jacvarsm1h_RE,jacvarsp1h_RE,jacvarsp0_RE;
  sqmat_t lambdaminusp1h, lambdaplusm1h;

  if (RESTRAINED) fatal_error("find_Fstar_interface_FDS_muscl_with_CDF can not be used for now in FDS RESTRAINED mode.");
  find_jacvars_from_musclvars(musclvarsm1h, metrics, gl, theta, &jacvarsm1h_RE);
  find_jacvars_from_musclvars(musclvarsp1h, metrics, gl, theta, &jacvarsp1h_RE);
  
  find_jacvars_at_interface_from_jacvars(jacvarsm1h_RE, jacvarsp1h_RE, gl, theta, metrics, AVERAGING, &jacvarsp0_RE);
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_jacvars_at_interface_from_jacvars(jacvarsm1h, jacvarsp1h, gl, theta, metrics, AVERAGING, &jacvarsp0);

  find_L_from_jacvars(jacvarsp0_RE, metrics, Lp0_RE);
  find_Linv_from_jacvars(jacvarsp0_RE, metrics, Linvp0_RE);
  find_Lambda_from_jacvars(jacvarsp0_RE, metrics, lambdap0_RE);
  
  set_matrix_to_zero(lambdapp0_RE);
  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvarsp0_RE, metrics, EIGENVALCOND, lambdapp0_RE);
  find_Ustar_from_jacvars(jacvarsm1h_RE,  metrics, Um1h_RE);
  find_Ustar_from_jacvars(jacvarsp1h_RE,  metrics, Up1h_RE);

  find_L_from_jacvars(jacvarsm1h, metrics, Lm1h);
  find_Ustar_from_jacvars(jacvarsm1h,  metrics, Um1h);
  find_LUstar_from_jacvars(jacvarsm1h, metrics, LUm1h);

  find_L_from_jacvars(jacvarsp1h, metrics, Lp1h);
  find_Linv_from_jacvars(jacvarsm1h, metrics, Linvm1h);
  find_Linv_from_jacvars(jacvarsp1h, metrics, Linvp1h);
  find_Ustar_from_jacvars(jacvarsp1h,  metrics, Up1h);
  find_LUstar_from_jacvars(jacvarsp1h, metrics, LUp1h);

  find_Fstar_from_jacvars(jacvarsm1h_RE,  metrics, Fm1h_RE);
  find_Fstar_from_jacvars(jacvarsp1h_RE,  metrics, Fp1h_RE);
#ifdef _RESTIME_CDF
  find_LUstar_from_jacvars(jacvarsm1h_RE, metrics, LUm1h_RE);
  find_LUstar_from_jacvars(jacvarsp1h_RE, metrics, LUp1h_RE);
  find_Lambdaxt_interface(np, gl, lm1h, lp1h, theta, jacvarsp0_RE, metrics, lambdaxt_p0);
  find_Lambdaxt(np, gl, lm1h, jacvarsm1h_RE, metrics, lambdaxt_m1h);
  find_Lambdaxt(np, gl, lp1h, jacvarsp1h_RE, metrics, lambdaxt_p1h);
  find_Deltaxt_interface(np, gl, lm1h, lp1h, theta, jacvarsm1h, jacvarsp1h, metrics, deltaxt_p0);
  for (flux=0; flux<nf; flux++) lambdapp0_RE[flux][flux]+=-fabs(lambdaxt_p0[flux][flux])*deltaxt_p0[flux];
#endif

  multiply_matrix_and_vector(Lp0_RE,Um1h_RE,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp2[flux]=(lambdapp0_RE[flux][flux])*fluxtmp[flux];
  multiply_matrix_and_vector(Linvp0_RE,fluxtmp2,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp[flux]+=Fm1h_RE[flux];
  multiply_matrix_and_vector(Lm1h,fluxtmp,fluxtmp2);

  multiply_matrix_and_vector(Lp0_RE,Up1h_RE,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp3[flux]=(-lambdapp0_RE[flux][flux])*fluxtmp[flux];
  multiply_matrix_and_vector(Linvp0_RE,fluxtmp3,fluxtmp);
  for (flux=0; flux<nf; flux++) fluxtmp[flux]+=Fp1h_RE[flux];
  multiply_matrix_and_vector(Lp1h,fluxtmp,fluxtmp3);

  for (flux=0; flux<nf; flux++) {
    lambdaplusm1h[flux][flux]=
#ifdef _RESTIME_CDF
       -0.5*lambdaxt_m1h[flux][flux]*deltaxt_p0[flux]
#endif
       +0.5*fluxtmp2[flux]/notzero(LUm1h[flux],1e-99);
    lambdaminusp1h[flux][flux]=
#ifdef _RESTIME_CDF
       -0.5*lambdaxt_p1h[flux][flux]*deltaxt_p0[flux]
#endif
       +0.5*fluxtmp3[flux]/notzero(LUp1h[flux],1e-99);
  }

  find_Linv_from_jacvars(jacvarsm1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsm1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaplusm1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fm1h);

  find_Linv_from_jacvars(jacvarsp1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsp1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaminusp1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fp1h);

  for (flux=0; flux<nf; flux++) Fint[flux]=Fm1h[flux]+Fp1h[flux];
}


void find_Fstar_interface_FDS_muscl_without_CDF(gl_t *gl, long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, int AVERAGING, bool RESTRAINED, flux_t Fint){
  flux_t Fm1h,Fp1h,mattmp,Um1h,Up1h,fluxtmp,alphap0;
  jacvars_t jacvarsm1h,jacvarsp1h,jacvarsp0;
  sqmat_t Linv,lambdap,L;
  long flux;

  find_jacvars_from_musclvars(musclvarsm1h, metrics, gl, theta, &jacvarsm1h);
  find_jacvars_from_musclvars(musclvarsp1h, metrics, gl, theta, &jacvarsp1h);
  find_jacvars_at_interface_from_jacvars(jacvarsm1h, jacvarsp1h, gl, theta, metrics, AVERAGING, &jacvarsp0);
  
  if (RESTRAINED) {
    find_L_restrained_from_jacvars(jacvarsp0, metrics, L);
    find_Linv_restrained_from_jacvars(jacvarsp0, metrics, Linv);
    find_conditioned_Lambda_absolute_restrained_from_jacvars(gl, jacvarsp0, metrics, EIGENVALCOND, lambdap);
  } else { 
    find_L_from_jacvars(jacvarsp0, metrics, L);
    find_Linv_from_jacvars(jacvarsp0, metrics, Linv);
    find_conditioned_Lambda_absolute_from_jacvars(gl, jacvarsp0, metrics, EIGENVALCOND, lambdap);
  }
  find_Ustar_from_musclvars(musclvarsm1h, metrics, gl, Um1h);
  find_Ustar_from_musclvars(musclvarsp1h, metrics, gl, Up1h);
  for (flux=0; flux<nf; flux++) mattmp[flux]=Up1h[flux]-Um1h[flux];
  multiply_matrix_and_vector(L,mattmp,alphap0);

  for (flux=0; flux<nf; flux++) fluxtmp[flux]=-lambdap[flux][flux]*alphap0[flux];
  multiply_matrix_and_vector(Linv,fluxtmp,Fint);

  find_Fstar_from_jacvars(jacvarsm1h, metrics, Fm1h);
  find_Fstar_from_jacvars(jacvarsp1h, metrics, Fp1h);

  for (flux=0; flux<nf; flux++){
    Fint[flux]=0.5*(Fint[flux]+Fm1h[flux]+Fp1h[flux]);
  }
}


void find_Fstar_interface_FDS_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h,  long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, int AVERAGING, flux_t Fint){
  bool RESTRAINED;
  RESTRAINED=FALSE;
#ifdef _RESTIME_CDF
  find_Fstar_interface_FDS_muscl_with_CDF(np, gl, lm1h, lp1h,  theta, musclvarsm1h, musclvarsp1h,
                                          metrics, EIGENVALCOND, AVERAGING, RESTRAINED, Fint);
#else
  find_Fstar_interface_FDS_muscl_without_CDF(gl, theta, musclvarsm1h, musclvarsp1h,
                                             metrics, EIGENVALCOND, AVERAGING, RESTRAINED, Fint);
#endif
}





void find_Fstar_interface_FDSR_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h,  long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, int AVERAGING, flux_t Fint){
  bool RESTRAINED;
  RESTRAINED=TRUE;
#ifdef _RESTIME_CDF
  find_Fstar_interface_FDS_muscl_with_CDF(np, gl, lm1h, lp1h,  theta, musclvarsm1h, musclvarsp1h,
                                          metrics, EIGENVALCOND, AVERAGING, RESTRAINED, Fint);
#else
  find_Fstar_interface_FDS_muscl_without_CDF(gl, theta, musclvarsm1h, musclvarsp1h,
                                             metrics, EIGENVALCOND, AVERAGING, RESTRAINED, Fint);
#endif
}




void find_Fstar_interface_FVS_muscl_without_CDF(gl_t *gl, long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, flux_t Fint){
  flux_t Fplusm1h,Fminusp1h,Um1h,Up1h,fluxtmp1,fluxtmp2;
  jacvars_t jacvarsp1h,jacvarsm1h;
  sqmat_t Linvp1h,Linvm1h,Lambdap1h,Lambdam1h,Lambdaabsm1h,Lambdaabsp1h,Lp1h,Lm1h,Lambdaplusm1h,Lambdaminusp1h;
  long flux;

  find_jacvars_from_musclvars(musclvarsm1h, metrics, gl, theta, &jacvarsm1h);
  find_jacvars_from_musclvars(musclvarsp1h, metrics, gl, theta, &jacvarsp1h);
  
  find_L_from_jacvars(jacvarsm1h, metrics, Lm1h);
  find_L_from_jacvars(jacvarsp1h, metrics, Lp1h);
  find_Linv_from_jacvars(jacvarsm1h, metrics, Linvm1h);
  find_Linv_from_jacvars(jacvarsp1h, metrics, Linvp1h);
  find_Ustar_from_musclvars(musclvarsm1h, metrics, gl, Um1h);
  find_Ustar_from_musclvars(musclvarsp1h, metrics, gl, Up1h);
  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvarsp1h, metrics, EIGENVALCOND, Lambdaabsp1h);
  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvarsm1h, metrics, EIGENVALCOND, Lambdaabsm1h);
  find_Lambda_from_jacvars(jacvarsp1h, metrics, Lambdap1h);
  find_Lambda_from_jacvars(jacvarsm1h, metrics, Lambdam1h);
  set_matrix_to_zero(Lambdaplusm1h);
  set_matrix_to_zero(Lambdaminusp1h);
  for (flux=0; flux<nf; flux++){
    Lambdaplusm1h[flux][flux]=0.5*(Lambdam1h[flux][flux]+Lambdaabsm1h[flux][flux]);
    Lambdaminusp1h[flux][flux]=0.5*(Lambdap1h[flux][flux]-Lambdaabsp1h[flux][flux]);
  }
  multiply_matrix_and_vector(Lm1h,Um1h,fluxtmp1);
  multiply_matrix_and_vector(Lambdaplusm1h,fluxtmp1,fluxtmp2);
  multiply_matrix_and_vector(Linvm1h,fluxtmp2,Fplusm1h);

  multiply_matrix_and_vector(Lp1h,Up1h,fluxtmp1);
  multiply_matrix_and_vector(Lambdaminusp1h,fluxtmp1,fluxtmp2);
  multiply_matrix_and_vector(Linvp1h,fluxtmp2,Fminusp1h);

  for (flux=0; flux<nf; flux++){
    Fint[flux]=Fplusm1h[flux]+Fminusp1h[flux];
  }
}


void find_Fstar_interface_FVS_muscl_with_CDF(np_t *np, gl_t *gl, long lm1h, long lp1h, long theta, 
                     flux_t musclvarsm1h, flux_t musclvarsp1h, metrics_t metrics,  
                     int EIGENVALCOND, int AVERAGING, flux_t Fint){
  flux_t Fm1h,Fp1h,fluxtmp,fluxtmp2;
  sqmat_t R;
  long flux;
  sqmat_t lambdaminusp1h, lambdaplusm1h;
  jacvars_t jacvarsm1h,jacvarsp1h,jacvarsp0;
  jacvars_t jacvarsm1h_RE,jacvarsp1h_RE,jacvarsp0_RE;
  sqmat_t Lm1h,Linvm1h,Linvp1h,Lp1h,Linvp1h_RE,Linvm1h_RE,
          Lambdaabsp1h_RE,Lambdaabsm1h_RE,Lambdap1h_RE,Lambdam1h_RE;
  flux_t Um1h,fluxtmp3,LUm1h,LUp1h,Up1h,LUm1h_RE,LUp1h_RE;
#ifdef _RESTIME_CDF
  sqmat_t lambdaxt_m1h,lambdaxt_p1h;
  flux_t Deltaxt_p0;
#endif

  find_jacvars_from_musclvars(musclvarsm1h, metrics, gl, theta, &jacvarsm1h_RE);
  find_jacvars_from_musclvars(musclvarsp1h, metrics, gl, theta, &jacvarsp1h_RE);
  find_jacvars_at_interface_from_jacvars(jacvarsm1h_RE, jacvarsp1h_RE, gl, theta, metrics, AVERAGING, &jacvarsp0_RE);
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_jacvars_at_interface_from_jacvars(jacvarsm1h, jacvarsp1h, gl, theta, metrics, AVERAGING, &jacvarsp0);

#ifdef _RESTIME_CDF
  find_Lambdaxt(np, gl, lm1h, jacvarsm1h_RE, metrics, lambdaxt_m1h);
  find_Lambdaxt(np, gl, lp1h, jacvarsp1h_RE, metrics, lambdaxt_p1h);
  find_Deltaxt_interface(np, gl, lm1h, lp1h, theta, jacvarsm1h, jacvarsp1h,  metrics, Deltaxt_p0);
#endif

  find_Linv_from_jacvars(jacvarsm1h_RE, metrics, Linvm1h_RE);
  find_Linv_from_jacvars(jacvarsp1h_RE, metrics, Linvp1h_RE);
  find_Lambda_from_jacvars(jacvarsm1h_RE, metrics, Lambdam1h_RE);
  find_Lambda_from_jacvars(jacvarsp1h_RE, metrics, Lambdap1h_RE);
  set_matrix_to_zero(Lambdaabsm1h_RE);
  set_matrix_to_zero(Lambdaabsp1h_RE);
  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvarsm1h_RE, metrics, EIGENVALCOND, Lambdaabsm1h_RE);
  find_conditioned_Lambda_absolute_from_jacvars(gl, jacvarsp1h_RE, metrics, EIGENVALCOND, Lambdaabsp1h_RE);
  find_LUstar_from_jacvars(jacvarsm1h_RE,  metrics, LUm1h_RE);
  find_LUstar_from_jacvars(jacvarsp1h_RE,  metrics, LUp1h_RE);

  find_L_from_jacvars(jacvarsm1h, metrics, Lm1h);
  find_Ustar_from_jacvars(jacvarsm1h,  metrics, Um1h);
  find_LUstar_from_jacvars(jacvarsm1h, metrics, LUm1h);

  find_L_from_jacvars(jacvarsp1h, metrics, Lp1h);
  find_Linv_from_jacvars(jacvarsm1h, metrics, Linvm1h);
  find_Linv_from_jacvars(jacvarsp1h, metrics, Linvp1h);
  find_Ustar_from_jacvars(jacvarsp1h,  metrics, Up1h);
  find_LUstar_from_jacvars(jacvarsp1h, metrics, LUp1h);

  for (flux=0; flux<nf; flux++) fluxtmp2[flux]=(Lambdam1h_RE[flux][flux]+Lambdaabsm1h_RE[flux][flux])*LUm1h_RE[flux];
#ifdef _RESTIME_CDF
  for (flux=0; flux<nf; flux++) fluxtmp2[flux]+=-(lambdaxt_m1h[flux][flux]+fabs(lambdaxt_m1h[flux][flux]))*Deltaxt_p0[flux]*LUm1h_RE[flux];
#endif

  multiply_matrix_and_vector(Linvm1h_RE,fluxtmp2,fluxtmp);
  multiply_matrix_and_vector(Lm1h,fluxtmp,fluxtmp2);

  for (flux=0; flux<nf; flux++) fluxtmp3[flux]=(Lambdap1h_RE[flux][flux]-Lambdaabsp1h_RE[flux][flux])*LUp1h_RE[flux];
#ifdef _RESTIME_CDF
  for (flux=0; flux<nf; flux++) fluxtmp3[flux]+=-(lambdaxt_p1h[flux][flux]-fabs(lambdaxt_p1h[flux][flux]))*Deltaxt_p0[flux]*LUp1h_RE[flux];
#endif

  multiply_matrix_and_vector(Linvp1h_RE,fluxtmp3,fluxtmp);
  multiply_matrix_and_vector(Lp1h,fluxtmp,fluxtmp3);

  for (flux=0; flux<nf; flux++){
    lambdaplusm1h[flux][flux]=0.5*fluxtmp2[flux]/notzero(LUm1h[flux],1e-99);
    lambdaminusp1h[flux][flux]=0.5*fluxtmp3[flux]/notzero(LUp1h[flux],1e-99);
  }
 
  find_Linv_from_jacvars(jacvarsm1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsm1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaplusm1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fm1h);

  find_Linv_from_jacvars(jacvarsp1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsp1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaminusp1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fp1h);

  for (flux=0; flux<nf; flux++) Fint[flux]=Fm1h[flux]+Fp1h[flux];

}


void find_Fstar_interface_FVS_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h,  long theta, flux_t musclvarsm1h, flux_t musclvarsp1h,
                     metrics_t metrics, int EIGENVALCOND, int AVERAGING, flux_t Fint){
#ifdef _RESTIME_CDF
  find_Fstar_interface_FVS_muscl_with_CDF(np, gl, lm1h, lp1h,  theta, musclvarsm1h, musclvarsp1h,
                                          metrics, EIGENVALCOND, AVERAGING, Fint);
#else
  find_Fstar_interface_FVS_muscl_without_CDF(gl, theta, musclvarsm1h, musclvarsp1h,
                                             metrics, EIGENVALCOND, Fint);
#endif

}


static void find_Lambda_minus_plus_FVSplus_muscl(np_t *np, gl_t *gl, long lp0, long lp1, jacvars_t jacvarsp0, jacvars_t jacvarsp1h, jacvars_t jacvarsp1, jacvars_t jacvarsp0_RE, jacvars_t jacvarsp1h_RE, jacvars_t jacvarsp1_RE, metrics_t metrics, long theta, long numiter, int EIGENVALCOND, sqmat_t lambdaminus, sqmat_t lambdaplus){
  sqmat_t Yminus,Yplus,Zplus,Zminus,Lp0,Linvp0,Linvp1,Lp1,Linvp1_RE,Linvp0_RE,
          Lambdaabsp1_RE,Lambdaabsp0_RE,Lambdap1_RE,Lambdap0_RE;
  long row,col,flux,cnt;
  flux_t Up0,fluxtmp,fluxtmp2,fluxtmp3,LUp0,LUp1,Up1,LUp0_RE,LUp1_RE;
#ifdef _RESTIME_CDF
  sqmat_t lambdaxt_p0,lambdaxt_p1;
  flux_t Deltaxt_p1h;
#endif

#ifdef _RESTIME_CDF
  find_Lambdaxt(np, gl, lp0, jacvarsp0_RE, metrics, lambdaxt_p0);
  find_Lambdaxt(np, gl, lp1, jacvarsp1_RE, metrics, lambdaxt_p1);
  find_Deltaxt_interface(np, gl, lp0, lp1, theta, jacvarsp0, jacvarsp1,  metrics, Deltaxt_p1h);
#endif

  find_Linv_from_jacvars(jacvarsp0_RE, metrics, Linvp0_RE);
  find_Linv_from_jacvars(jacvarsp1_RE, metrics, Linvp1_RE);
  find_Lambda_from_jacvars(jacvarsp0_RE, metrics, Lambdap0_RE);
  find_Lambda_from_jacvars(jacvarsp1_RE, metrics, Lambdap1_RE);
  set_matrix_to_zero(Lambdaabsp0_RE);
  set_matrix_to_zero(Lambdaabsp1_RE);
  for (flux=0; flux<nf; flux++) {
    Lambdaabsp0_RE[flux][flux]=fabs(Lambdap0_RE[flux][flux]);
    Lambdaabsp1_RE[flux][flux]=fabs(Lambdap1_RE[flux][flux]);
  }
  find_LUstar_from_jacvars(jacvarsp0_RE,  metrics, LUp0_RE);
  find_LUstar_from_jacvars(jacvarsp1_RE,  metrics, LUp1_RE);

  find_L_from_jacvars(jacvarsp0, metrics, Lp0);
  find_Ustar_from_jacvars(jacvarsp0,  metrics, Up0);
  find_LUstar_from_jacvars(jacvarsp0, metrics, LUp0);

  find_L_from_jacvars(jacvarsp1, metrics, Lp1);
  find_Linv_from_jacvars(jacvarsp0, metrics, Linvp0);
  find_Linv_from_jacvars(jacvarsp1, metrics, Linvp1);
  find_Ustar_from_jacvars(jacvarsp1,  metrics, Up1);
  find_LUstar_from_jacvars(jacvarsp1, metrics, LUp1);

  for (flux=0; flux<nf; flux++) fluxtmp2[flux]=(Lambdap0_RE[flux][flux]+Lambdaabsp0_RE[flux][flux])*LUp0_RE[flux];
#ifdef _RESTIME_CDF
  for (flux=0; flux<nf; flux++) fluxtmp2[flux]+=-(lambdaxt_p0[flux][flux]+fabs(lambdaxt_p0[flux][flux]))*Deltaxt_p1h[flux]*LUp0_RE[flux];
#endif

  multiply_matrix_and_vector(Linvp0_RE,fluxtmp2,fluxtmp);
  multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);

  for (flux=0; flux<nf; flux++) fluxtmp3[flux]=(Lambdap1_RE[flux][flux]-Lambdaabsp1_RE[flux][flux])*LUp1_RE[flux];
#ifdef _RESTIME_CDF
  for (flux=0; flux<nf; flux++) fluxtmp3[flux]+=-(lambdaxt_p1[flux][flux]-fabs(lambdaxt_p1[flux][flux]))*Deltaxt_p1h[flux]*LUp1_RE[flux];
#endif

  multiply_matrix_and_vector(Linvp1_RE,fluxtmp3,fluxtmp);
  multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);

  for (row=0; row<nf; row++){
    for (col=0; col<nf; col++){
      lambdaminus[row][col]=0.0;
      lambdaplus[row][col]=0.0;
      Yminus[row][col]=0.0;
      Yplus[row][col]=0.0;
      Zminus[row][col]=0.0;
      Zplus[row][col]=0.0;
    }
  }

  for (flux=0; flux<nf; flux++) {
    Yminus[flux][flux]=min(0.0,
       +0.5*fluxtmp2[flux]/notzero(LUp0[flux],1e-99) );
    Yplus[flux][flux]=max(0.0,
       +0.5*fluxtmp2[flux]/notzero(LUp0[flux],1e-99) );
    Zplus[flux][flux]=max(0.0,
       +0.5*fluxtmp3[flux]/notzero(LUp1[flux],1e-99) );
    Zminus[flux][flux]=min(0.0,
       +0.5*fluxtmp3[flux]/notzero(LUp1[flux],1e-99) );
  }

  for (cnt=0; cnt<numiter; cnt++){

    multiply_matrix_and_vector(Zplus,LUp1,fluxtmp2);
    multiply_matrix_and_vector(Linvp1,fluxtmp2,fluxtmp);
    multiply_matrix_and_vector(Lp0,fluxtmp,fluxtmp2);

    multiply_matrix_and_vector(Yminus,LUp0,fluxtmp3);
    multiply_matrix_and_vector(Linvp0,fluxtmp3,fluxtmp);
    multiply_matrix_and_vector(Lp1,fluxtmp,fluxtmp3);
    for (flux=0; flux<nf; flux++) {
      Yminus[flux][flux]=min(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1e-99));
      Zplus[flux][flux]=max(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1e-99));
    }
    for (flux=0; flux<nf; flux++) {
      Yplus[flux][flux]=max(0.0, Yplus[flux][flux]+fluxtmp2[flux]/notzero(LUp0[flux],1e-99));
      Zminus[flux][flux]=min(0.0, Zminus[flux][flux]+fluxtmp3[flux]/notzero(LUp1[flux],1e-99));
    }
  
  }

  for (flux=0; flux<nf; flux++) {
    lambdaplus[flux][flux]=Yplus[flux][flux]+Zplus[flux][flux];
    lambdaminus[flux][flux]=Zminus[flux][flux]+Yminus[flux][flux];
  }
  
  condition_Lambda_plus_minus(np, gl, lp0, theta, jacvarsp0, jacvarsp1,metrics, EIGENVALCOND, lambdaplus,lambdaminus);

}


void find_Fstar_interface_FVSplus_muscl(np_t *np, gl_t *gl, long lm1h, long lp1h, long theta, 
                     flux_t musclvarsm1h, flux_t musclvarsp1h, metrics_t metrics,  long numiter, 
                     int EIGENVALCOND, int AVERAGING, flux_t Fint, sqmat_t lambdaminusp1h, sqmat_t lambdaplusm1h){
  flux_t Fm1h,Fp1h,fluxtmp,fluxtmp2;
  sqmat_t R;
  long flux;
  jacvars_t jacvarsm1h,jacvarsp1h,jacvarsp0;
  jacvars_t jacvarsm1h_RE,jacvarsp1h_RE,jacvarsp0_RE;

  find_jacvars_from_musclvars(musclvarsm1h, metrics, gl, theta, &jacvarsm1h_RE);
  find_jacvars_from_musclvars(musclvarsp1h, metrics, gl, theta, &jacvarsp1h_RE);
  find_jacvars_at_interface_from_jacvars(jacvarsm1h_RE, jacvarsp1h_RE, gl, theta, metrics, AVERAGING, &jacvarsp0_RE);
  find_jacvars(np[lm1h],gl,metrics,theta,&jacvarsm1h);
  find_jacvars(np[lp1h],gl,metrics,theta,&jacvarsp1h);
  find_jacvars_at_interface_from_jacvars(jacvarsm1h, jacvarsp1h, gl, theta, metrics, AVERAGING, &jacvarsp0);

  find_Lambda_minus_plus_FVSplus_muscl(np,gl,lm1h,lp1h,jacvarsm1h,jacvarsp0,jacvarsp1h, jacvarsm1h_RE,jacvarsp0_RE,jacvarsp1h_RE,metrics, theta, numiter, EIGENVALCOND, lambdaminusp1h, lambdaplusm1h);

  find_Linv_from_jacvars(jacvarsm1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsm1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaplusm1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fm1h);

  find_Linv_from_jacvars(jacvarsp1h, metrics, R);
  find_LUstar_from_jacvars(jacvarsp1h, metrics, fluxtmp);
  multiply_diagonal_matrix_and_vector(lambdaminusp1h,fluxtmp,fluxtmp2);
  multiply_matrix_and_vector(R,fluxtmp2,Fp1h);

  for (flux=0; flux<nf; flux++) Fint[flux]=Fm1h[flux]+Fp1h[flux];

}


