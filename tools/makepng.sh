#!/bin/bash

# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright 2016,2022 Bernard Parent
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


#set default counter increment
increm=1;


#set the minimum number of digits used to write png files
numdigits=1


#set default image width (in pixels)
imagewidth=800

#set default font size (in points)
fontsize=22

#set the default modulo to 1 (don't skip any post file)
modulo=1

#set default number of maximum threads to the number of hyperthreads obtained from nproc
#maxthreads=`grep -c ^processor /proc/cpuinfo`;

#set maxthreads to 1 here because 2 tec360 processes can't usually be run at the same time because of the license restrictions
maxthreads=1;

reorder=0;
reorderstart=1;

while [ $# -gt 0 ]
do
    case "$1" in
        -data) datafilename="$2"; shift;;
	-post)  postfilename="$2"; shift;;
	-style)  stylefilename="$2"; shift;;
	-png)  pngfilename="$2"; shift;;
	-macro)  macrofilename="$2"; shift;;
	-digits)  numdigits="$2"; shift;;
	-start)  is="$2"; shift;;
	-end)  ie="$2"; shift;;
	-increm)  increm="$2"; shift;;
	-imagewidth)  imagewidth="$2"; shift;;
	-fontsize)  fontsize="$2"; shift;;
	-modulo)  modulo="$2"; shift;;
	-reorderstart)  reorderstart="$2"; shift;;
	-reorder)  reorder=1;;
	-threads)  maxthreads="$2"; shift;;
	--)	shift; break;;
	*)  break;;	# terminate while loop
    esac
    shift
done

if [ -z "$datafilename" ]; then 
  argerror="1"; 
fi

if [ -z "$postfilename" ]; then 
  argerror="1"; 
fi

if [ -z "$pngfilename" ]; then 
  argerror="1"; 
fi

if [ -z "$stylefilename" ]; then 
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
Flag           Arg                                      Required?
_________________________________________________________________

-data          data file prefix [string]                Y
-post          post file prefix [string]                Y
-png           png file prefix [string]                 Y
-start         counter start [int]                      Y
-end           counter end [int]                        Y
-style         tecplot style file [string]              Y
-digits        minimum number of digits used when       N 
               writing png files
-macro         tecplot macro file [string]              N
-modulo        modulo used to skip files read [int]     N
               default: 1
-fontsize      font size in points [int]                N
-imagewidth    width of the image in pixels [int]       N
-increm        i increments [int]                       N
               (i=is, is+increm, is+2*increm, ..)
-threads       max. number of png files made in         N 
               parallel (make sure maxthreads doesn't
               exceed the number of tec360 processes
               your license allows you to run)
-reorder       reorder the png files so that they are   N
               named 1, 2, 3, etc, no matter the 
               values given to -start and -increm 
-reorderstart  counter start when reordering [int]      N

Eg: 
${0##*/} -data data. -post post. -start 10 -end 12 -style rho.sty -png rho.
will read in data.10,data.11,data.12 and post.10,post.11,post.12 and make png files rho.1.png, rho.2.png, etc according to the tecplot style rho.sty " 

  exit 1
fi




# make sure that convert, and tec360 are installed on the system
command -v tec360 >/dev/null 2>&1 || { echo "makepng requires tec360 but it is not installed.  Aborting." >&2; exit 1; }
command -v convert >/dev/null 2>&1 || { echo "makepng requires convert but it is not installed.  Aborting." >&2; exit 1; }


# make sure the style file  exists
if [ -n "$stylefilename" ]
then
  if [ -f "$stylefilename" ]
   then
     echo Style file $stylefilename found.
   else
     echo Style file $stylefilename not found. Exiting.
     exit 1
  fi
fi

if [ -n "$macrofilename" ]
then
  if [ -f "$macrofilename" ]
   then
     echo Macro file $macrofilename found.
   else
     echo Macro file $macrofilename not found. Exiting.
     exit 1
  fi
fi





i=$is;
thread=$reorderstart;
cnt=1;

while [ $i -le $ie ]
 do

   if [ $reorder -eq 0 ];
    then
      size=${#i}
      iw="$i";
    else
      size=${#thread}
      iw="$thread";
   fi
   for ((m=1; m<=numdigits-$size; m++)); do
     iw="0""$iw"
   done

   if [ ! -f "$postfilename$i" ]
    then
      echo File $postfilename$i not found. Skipping.
   else


    #create the macro file that tecplot will read 
    echo "#!MC 1410" > __makepng.$iw.mcr

    if [ -n "$stylefilename" ]
    then
      echo "\$!ReadStyleSheet  \"./$stylefilename\"" >> __makepng.$iw.mcr
    fi

    echo "  IncludePlotStyle = Yes
  IncludeText = Yes
  IncludeGeom = Yes
  IncludeAuxData = Yes
  IncludeStreamPositions = Yes
  IncludeContourLevels = Yes
  Merge = No
  IncludeFrameSizeAndPosition = No
\$!PrintSetup Palette = Color
\$!ExportSetup ImageWidth = "$imagewidth"
\$!ExportSetup ExportFName = './$pngfilename$iw.png'
\$!Export
  ExportRegion = AllFrames" >> __makepng.$iw.mcr



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

    time=`head -1 $datafilename$i | grep -o -P '(?<=time=).*(?= dt)'`
   

    # create the png using tecplot
    cnt=$(( cnt+1 ))
    if [ $(( cnt % $modulo)) -eq 0 ]; then
      echo Creating the file $pngfilename$iw.png from $postfilename$i [$time s]  
      if [ -z "$macrofilename" ]; then 
        tecplotexec="tec360 --osmesa -b -p __makepng.$iw.mcr $postfilename$i"
      else
        tecplotexec="tec360 --osmesa -b -p $macrofilename -p __makepng.$iw.mcr $postfilename$i"     
      fi
      echo $tecplotexec
      convertexec="convert $pngfilename$iw.png -font Times-Roman -pointsize $fontsize -draw \"gravity north fill white  text 0,15 '$time s' fill black  text 5,10 '$time s'\" $pngfilename$iw.png"
      echo $convertexec
      if [ $maxthreads -eq 1 ]; then
        eval $tecplotexec 1> /dev/null  
        eval $convertexec 1> /dev/null  &
      else
        eval $tecplotexec 1> /dev/null  && eval $convertexec 1> /dev/null  &
      fi
      thread=$(( thread+1 ))
    else 
      echo Skipping the creation of a png file from $postfilename$i..
    fi
   fi
    # increments $i
    i=$(( i+$increm ))	 

    # wait for the threads to complete
    if [ ! $maxthreads -eq 1 ];
    then
      if [ $(( thread % $maxthreads )) -eq 0 ];
      then
        wait
      fi
    fi

done    

#wait for the background processes to finish 
wait

#check if all png files were created successfully by tec360
echo Verifying that all png files were created successfully..
i=$is;
thread=$reorderstart;
cnt=1;

while [ $i -le $ie ]
do

 if [ ! -f "$postfilename$i" ]
   then
   garbage=2;
 else
   if [ $reorder -eq 0 ];
    then
      size=${#i}
      iw="$i";
    else
      size=${#thread}
      iw="$thread";
   fi

   for ((m=1; m<=numdigits-$size; m++)); do
     iw="0""$iw"
   done

   cnt=$(( cnt+1 ))
   if [ $(( cnt % $modulo)) -eq 0 ]; then
     if [ ! -f "$pngfilename$iw.png" ]
      then
       echo File $pngfilename$iw.png not created properly. This is probably due to some license error when running tec360. Exiting.
       exit 1
     fi
     thread=$(( thread+1 ))
   fi
 fi
   # increments $i
   i=$(( i+$increm ))	 
done    
echo Done.


# cleanup
rm -f __makepng.*.mcr


