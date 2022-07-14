#!/bin/bash

# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright 2016 Bernard Parent
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of 
#    conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list 
#    of conditions and the following disclaimer in the documentation and/or other 
#    materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#check that the number of arguments is correct; if not, give some help

reorder=0;
reorderstart=1;
newdigits=1;

while [ $# -gt 0 ]
do
    case "$1" in
        -prefix) prefix="$2"; shift;;
	-suffix)  suffix="$2"; shift;;
	-newprefix)  newprefix="$2"; shift;;
	-newsuffix)  newsuffix="$2"; shift;;
	-digits)  digits="$2"; shift;;
	-newdigits)  newdigits="$2"; shift;;
	-reorder)  reorder=1;;
	-start)  is="$2"; shift;;
	-end)  ie="$2"; shift;;
	--)	shift; break;;
	*)  break;;	# terminate while loop
    esac
    shift
done

if [ -z "$prefix" ]; then 
  argerror="1"; 
fi

if [ -z "$suffix" ]; then 
  argerror="1"; 
fi

if [ -z "$digits" ]; then 
  argerror="1"; 
fi

if [ -z "$is" ]; then 
  argerror="1"; 
fi

if [ -z "$ie" ]; then 
  argerror="1"; 
fi

if [ -z "$newprefix" ]; then 
  newprefix=$prefix; 
fi

if [ -z "$newsuffix" ]; then 
  newsuffix=$suffix; 
fi


if [ -n "$argerror" ]; then
  echo "
Flag           Arg                                      Required?
_________________________________________________________________

-prefix        file prefix [string]                     Y
-suffix        file suffix [string]                     Y
-newprefix     new file prefix [string]                 N
-newsuffix     new file suffix [string]                 N
-start         counter start [int]                      Y
-end           counter end [int]                        Y
-reorder       reorder the new files so that they       N
               numbered 1, 2, 3, etc, no matter the 
               numbers found in original files 
-digits        number of digits the number has in       N 
               original files
-newdigits     number of digits the number has in       N 
               new files

Eg: 
${0##*/} -prefix 'vlc.' -suffix '.jpg' -start 1 -end 10 -digits 3" 

  exit 1
fi


i=$is;
thread=$reorderstart;

while [ $i -le $ie ]
 do

   size=${#i}
   iold="$i";
   for ((m=1; m<=digits-$size; m++)); do
     iold="0""$iold"
   done

   if [ ! -f "$prefix$iold$suffix" ]
    then
      echo File $prefix$iold$suffix not found. Skipping.
   else
    if [ $reorder -eq 0 ];
    then
      inew="$i";
    else
      inew="$thread";
    fi
    size=${#inew}
    for ((m=1; m<=newdigits-$size; m++)); do
      inew="0""$inew"
    done

    echo Renaming $prefix$iold$suffix to $newprefix$inew$newsuffix ..
    mv $prefix$iold$suffix $newprefix$inew$newsuffix
    # increments $thread
    thread=$(( thread+1 ))	 
   fi
   # increments $i
   i=$(( i+1 ))	 

done    

echo Done.




