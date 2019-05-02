#!/bin/bash

# SPDX-License-Identifier: BSD-2-Clause

# Copyright 2018 Bernard Parent

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#check that the number of arguments is correct; if not, give some help


#set default counter increment
increm=1;

#set default number of threads to the number of hyperthreads obtained from nproc
#maxthreads=`grep -c ^processor /proc/cpuinfo`;
maxthreads=1;

mindigits=0;

while [ $# -gt 0 ]
do
    case "$1" in
	-output)  output="$2"; shift;;
	-eqn)  equation="$2"; shift;;
	-data)  datafilename="$2"; shift;;
	-post)  postfilename="$2"; shift;;
	-macro)  macrofilename="$2"; shift;;
	-threads)  maxthreads="$2"; shift;;
	-digits)  mindigits="$2"; shift;;
	-start)  is="$2"; shift;;
	-end)  ie="$2"; shift;;
	-increm)  increm="$2"; shift;;
	--)	shift; break;;
	*)  break;;	# terminate while loop
    esac
    shift
done

if [ -z "$output" ]; then 
  argerror="1"; 
fi


if [ -z "$postfilename" ]; then 
  argerror="1"; 
fi


if [ -z "$datafilename" ]; then 
  argerror="1"; 
fi


if [ -z "$equation" ] && [ -z "$macrofilename" ]; then 
  argerror="1"; 
fi

if [ ! -z "$equation" ] && [ ! -z "$macrofilename" ]; then 
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
Flag         Arg                                 Required?
---------------------------------------------------------------
-output      output text file [string]           Y
-data        data file prefix [string]           Y
-post        post file prefix [string]           Y
-start       counter start [int]                 Y
-end         counter end [int]                   Y
-eqn         equation [string]                   Y*
-macro       tecplot macro file [string]         Y*
-digits      minimum number of digits of post    N
             file counter
-increm      counter increments [int]            N
             (i=start, start+increm, ..)
-threads     max. number of integrations made in N 
             parallel

* either one of these flags must be given but not both at the same time

Eg: 
$0 -post post. -data out.01. -start 10 -end 30 -eqn '0.5*{rho}*({V[0]}*{V[0]}+{V[1]}*{V[1]}+{V[2]}*{V[2]})' -output out.txt
will average the property specified by the equation string within the post files post.10, post.11, ..., post.30
and output them to the output text file out.txt as a function of time."
 
  exit 1
fi



# make sure that the sed, calc executable and tecplot's tec360 are available
command -v tec360 >/dev/null 2>&1 || { echo "integrate.sh requires tec360 but it is not installed.  Aborting." >&2; exit 1; }
command -v sed >/dev/null 2>&1 || { echo "integrate.sh requires sed but it is not installed.  Aborting." >&2; exit 1; }
command -v calc >/dev/null 2>&1 || { echo "integrate.sh requires calc but it is not installed.  Aborting." >&2; exit 1; }
command -v grep >/dev/null 2>&1 || { echo "integrate.sh requires grep but it is not installed.  Aborting." >&2; exit 1; }



i=$is;
thread=0;

while [ $i -le $ie ]
 do
    size=${#i}
    iw="$i"
    for ((m=1; m<=mindigits-$size; m++)); do
      iw="0""$iw"
    done





    if [ ! -f "$postfilename$i" ]
    then
      echo File $postfilename$i not found. Skipping.
    else


      if [ -z "$macrofilename" ]; then 
        echo "#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454
\$!VarSet |MFBD| = './'
\$!ALTERDATA 
  EQUATION = '{newvar}=$equation'
\$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer3'
  COMMAND = 'Integrate VariableOption=\'Average\' ScalarVar=|NUMVARS| XVariable=1 YVariable=2 ZVariable=3'
\$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer3'
  COMMAND = 'SaveIntegrationResults FileName=\'./__integrate.$iw.txt\''
\$!RemoveVar |MFBD|"  > __integrate.$iw.mcr
        tecplotexec="tec360 -b -p __integrate.$iw.mcr $postfilename$i"
      else
        echo "#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454
\$!VarSet |MFBD| = './'
\$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer3'
  COMMAND = 'SaveIntegrationResults FileName=\'./__integrate.$iw.txt\''
\$!RemoveVar |MFBD|"  > __integrate.$iw.mcr
        tecplotexec="tec360 -b -p $macrofilename -p __integrate.$iw.mcr $postfilename$i"     
      fi

      echo $tecplotexec 
      $tecplotexec 1> /dev/null  &  
      thread=$(( thread+1 ))
    fi
    # increments $i
    i=$(( i+$increm ))	 

    # wait for the threads to complete
    if [ $(( thread % $maxthreads )) -eq 0 ];
    then
      wait
    fi


done  
wait  

rm -f $output
echo "Writing to output file $output.."
i=$is;
time='0';
prop='0';
while [ $i -le $ie ]
 do
    size=${#i}
    iw="$i"
    for ((m=1; m<=mindigits-$size; m++)); do
      iw="0""$iw"
    done
    if [  -f "$postfilename$i" ]
    then
      # find the time
      if [ ! -f "$datafilename$i" ]
      then
        echo File $datafilename$i not found. Exiting.
        exit 1
      fi

      if grep -q dataformat007 "$datafilename$i"; then
          echo WARP data file $datafilename$i has data format 7 but requires data format 8 or higher 
          exit 1
      fi

      oldtime=$time;
      time=`head -1 $datafilename$i | grep -o -P '(?<=time=).*(?= dt)'`
      oldprop=$prop;
      prop=`tail -1 __integrate.$iw.txt | sed 's/.*Total://'`
      if [ $oldtime != '0' ]
      then
        dpropdtime=`calc "($prop-$oldprop)/($time-$oldtime)"`
      else
        dpropdtime='0';
      fi
      echo "$time  $prop  $dpropdtime" >> $output
    fi
    # increments $i
    i=$(( i+$increm ))	 
 done
rm -f __integrate.*.txt
rm -f __integrate.*.mcr
echo "Done."
