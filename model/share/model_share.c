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


