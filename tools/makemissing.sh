#!/bin/sh

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


while [ $# -gt 0 ]
do
    case "$1" in
        -file) filename="$2"; shift;;
	-start)  is="$2"; shift;;
	-end)  ie="$2"; shift;;
	--)	shift; break;;
	*)  break;;	# terminate while loop
    esac
    shift
done

if [ -z "$filename" ]; then 
  argerror="1"; 
fi

if [ -z "$is" ]; then 
  argerror="1"; 
fi

if [ -z "$ie" ]; then 
  argerror="1"; 
fi

if [ -n "$argerror" ]; then
  echo "
Flag         Arg                             Required?
-----------------------------------------------------------
-file        file prefix [string]            Y
-start       counter start [int]             Y
-end         counter end [int]               Y

Eg: 
${0##*/} -file data. -start 10 -end 30 
If data.11 is not present, will ln -s data.10 data.11
If data.12 is not present, will ln -s data.11 data.12
..and so on until data.30
" 
  exit 1
fi


#make sure the first file is present
    if [ ! -f "$filename$is" ]
    then
      echo File $filename$is not found. Exiting.
      exit 1
    fi


i=$(( $is+1 ))	 
while [ $i -le $ie ]
 do

    if [ ! -f "$filename$i" ]
    then
      im1=$(( $i-1 ))
      echo ln -s $filename$im1 $filename$i
      ln -s $filename$im1 $filename$i
    fi

    # increments $i
    i=$(( i+1 ))	 
done  
wait  
