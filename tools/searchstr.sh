#!/bin/sh
###       search ()
###        {
   if [ $# = 1 ]
   then
   {
	echo \'$1\'
       for i in `find . -path './dev' -prune -o -print 2> /dev/null`
       do
       {
#           fgrep -i "$1" $i > /dev/null 2>&1;
           fgrep "$1" $i > /dev/null 2>&1;

           if [ $? = 0 ]
           then
           {
               echo $i
           } fi;
       } done;
   }
   else
   {
       echo "Wrong number of arguments!"
   } fi;
###        }
