// SPDX-License-Identifier: BSD-2-Clause
/*
Copyright 2016 Bernard Parent

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

#include <model/share/chem_share.h>
#include <model/chem/_chem.h>
#include <model/_model.h>
#include <model/thermo/_thermo.h>
#include <model/share/model_share.h>
#include <stdarg.h>




double _f_extrapol(int ACCURACY, ...) {
   double ret,val1,val2,val3;
   va_list valist;

   /* initialize valist for num number of arguments */
   va_start(valist, ACCURACY);
   /* access all the arguments assigned to valist */
   val1=va_arg(valist, double);
   val2=val1;
   val3=val2;
   if (ACCURACY>1) val2=va_arg(valist, double);
   if (ACCURACY>2) val3=va_arg(valist, double);
   if (ACCURACY>3) fatal_error("Order of accuracy can not be greater than 3 in _f_extrapol");	
   /* clean memory reserved for valist */
   va_end(valist);
   ret=0.0;
   switch (ACCURACY){
     case ACCURACY_FIRSTORDER:
       ret=val1;
     break;
     case ACCURACY_SECONDORDER:
       ret=2.0*val1-1.0*val2;
     break;
     case ACCURACY_THIRDORDER:
       ret=3.0*val1-3.0*val2+1.0*val3;
     break;
     default:
       fatal_error("Order of accuracy can not be equal to %d in _f_extrapol",ACCURACY);
   }

   return(ret);
}



double _f_symmetry(int ACCURACY, ...) {
   double ret,val1,val2,val3;
   va_list valist;

   /* initialize valist for num number of arguments */
   va_start(valist, ACCURACY);
   /* access all the arguments assigned to valist */
   val1=va_arg(valist, double);
   val2=val1;
   val3=val2;
   if (ACCURACY>1) val2=va_arg(valist, double);
   if (ACCURACY>2) val3=va_arg(valist, double);
   if (ACCURACY>3) fatal_error("Order of accuracy can not be greater than 3 in _f_symmetry");	
   /* clean memory reserved for valist */
   va_end(valist);
   ret=0.0;
   switch (ACCURACY){
     case ACCURACY_FIRSTORDER:
       ret=val1;
     break;
     case ACCURACY_SECONDORDER:
       ret=4.0/3.0*val1-1.0/3.0*val2;
     break;
     case ACCURACY_THIRDORDER:
       ret=32.0/21.0*val1-4.0/7.0*val2+1.0/21.0*val3;
     break;
     default:
       fatal_error("Order of accuracy can not be equal to %d in _f_symmetry",ACCURACY);
   }

   return(ret);
}



#if defined(_AVERAGEDRATES) && defined(UNSTEADY)
double _averaged_rate(np_t np, gl_t *gl, averagedrates_id_type id, double rate){
  long id_map,cntfluid,cntchem;
  double newrate;
  bool FOUND;
  FOUND=FALSE;
  id_map=0; //to avoid compiler warning
  cntfluid=0;
#ifdef _AVERAGEDRATES_FLUID
  for (cntfluid=0; cntfluid<numaveragedrates_fluid; cntfluid++){
    if (strcmp(averagedrates_fluid_id[cntfluid],id)==0){
      if (!FOUND) FOUND=TRUE;
        else fatal_error("Problem in function _averaged_rate: two averagedrates IDs %ld can not be the same.",id);
      id_map=cntfluid;
    }
  }
#endif
#ifdef _AVERAGEDRATES_CHEM
  for (cntchem=0; cntchem<numaveragedrates_chem; cntchem++){
    if (strcmp(averagedrates_chem_id[cntchem],id)==0){
      if (!FOUND) FOUND=TRUE;
        else fatal_error("Problem in function _averaged_rate: two averagedrates IDs %ld can not be the same.",id);
      id_map=cntchem+cntfluid;
    }
  }
#endif
  assert(id_map>=0);
  assert(id_map<numaveragedrates);  

  if (!FOUND) {
    for (cntfluid=0; cntfluid<numaveragedrates_fluid; cntfluid++){
      printf("%s\n",averagedrates_fluid_id[cntfluid]);
    }
    for (cntchem=0; cntchem<numaveragedrates_chem; cntchem++){
      printf("%s\n",averagedrates_chem_id[cntchem]);
    }
    fatal_error("Problem finding id '%s' in averagedrates_fluid_id[] or averagedrates_chem_id[].",id);
  }
  switch (gl->AVERAGEDRATES){
    case AVERAGEDRATES_ADD:
      if (gl->averagedrates_time<0.0) fatal_error("Averaged rates not initialized properly. Need to call SetAveragedRatesToZero() first.");
      np.bs->averagedrates[id_map]=(np.bs->averagedrates[id_map]*gl->averagedrates_time+gl->dt*rate)/(gl->averagedrates_time+gl->dt);
      newrate=rate;
    break;
    case AVERAGEDRATES_SET:
      np.bs->averagedrates[id_map]=rate;
      newrate=rate;
    break;
    case AVERAGEDRATES_ON:
      if (gl->averagedrates_time<0.0) fatal_error("Averaged rates not initialized properly. Need to call SetAveragedRatesToZero() first.");
      newrate=np.bs->averagedrates[id_map];
    break;
    case AVERAGEDRATES_OFF:
      newrate=rate;
    break;
    default:
      newrate=0.0;
      fatal_error("gl->AVERAGEDRATES can not be set to %ld.",gl->AVERAGEDRATES);
  }
  return(newrate);
}
#else
double _averaged_rate(np_t np, gl_t *gl, averagedrates_id_type id, double rate){
  return(rate);
}
#endif

