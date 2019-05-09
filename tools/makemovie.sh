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


#set default counter increment
increm=1;

#set default frames per second
fps=24;

#set default number of maximum threads to the number of hyperthreads obtained from nproc
#maxthreads=`grep -c ^processor /proc/cpuinfo`;

#set maxthreads to 1 here because 2 tec360 processes can't usually be run at the same time because of the license restrictions
maxthreads=1;


#set default image width (in pixels)
imagewidth=800

#set default font size (in points)
fontsize=22

#set the modulo to 1 (don't skip any post file)
modulo=1


while [ $# -gt 0 ]
do
    case "$1" in
        -data) datafilename="$2"; shift;;
	-post)  postfilename="$2"; shift;;
	-style)  stylefilename="$2"; shift;;
	-movie)  moviefilename="$2"; shift;;
	-macro)  macrofilename="$2"; shift;;
	-imagewidth)  imagewidth="$2"; shift;;
	-fontsize)  fontsize="$2"; shift;;
	-modulo)  modulo="$2"; shift;;
	-start)  is="$2"; shift;;
	-end)  ie="$2"; shift;;
	-increm)  increm="$2"; shift;;
	-fps)  fps="$2"; shift;;
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

if [ -z "$moviefilename" ]; then 
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
Flag         Arg                                      Required?
--------------------------------------------------------------------
-data        data file prefix [string]                Y
-post        post file prefix [string]                Y
-start       counter start [int]                      Y
-end         counter end [int]                        Y
-movie       movie file to create [string]            Y
-fps         frames per second [int]                  N
-style       tecplot style file [string]              N
-macro       tecplot macro file [string]              N
-modulo      modulo used to skip files read [int]     N
             default: 1
-fontsize    font size in points [int]                N
-imagewidth  width of the image in pixels [int]       N
-increm      i increments [int]                       N
             (i=is, is+increm, is+2*increm, ..)
-threads     max. number of png files made in         N 
             parallel (make sure maxthreads doesn't
             exceed the number of tec360 processes
             your license allows you to run)

Eg: 
$0 -data data. -post post. -start 10 -end 12 -style rho.sty -movie rho.mp4
will read in data.10,data.11,data.12 and post.10,post.11,post.12 and make a movie rho.mp4 according to the tecplot style rho.sty" 

  exit 1
fi




# make sure that ffmpeg, convert, calc, and tec360 are installed on the system
command -v ffmpeg >/dev/null 2>&1 || { echo "makemovie requires ffmpeg but it is not installed.  Aborting." >&2; exit 1; }
command -v convert >/dev/null 2>&1 || { echo "makemovie requires convert but it is not installed.  Aborting." >&2; exit 1; }
command -v makepng >/dev/null 2>&1 || { echo "makemovie requires makepng but it is not installed.  Aborting." >&2; exit 1; }
command -v calc >/dev/null 2>&1 || { echo "makemovie requires calc but it is not installed.  Aborting." >&2; exit 1; }



#make sure there are no leftover __makemovie*png files from the previous run or those would end up being included later on by ffmpeg
rm -f __makemovie.*.png

if [ -z "$macrofilename" ]; then 
  macroarg="";
else
  macroarg="-macro $macrofilename";
fi

if [ -z "$stylefilename" ]; then 
  stylearg="";
else
  stylearg="-style $stylefilename";
fi


echo makepng -data $datafilename -post $postfilename  -png __makemovie. $stylearg $macroarg -digits 8 -start $is -end $ie -increm $increm -threads $maxthreads -modulo $modulo -imagewidth $imagewidth -fontsize $fontsize
 
makepng -data $datafilename -post $postfilename  -png __makemovie. $stylearg $macroarg -digits 8 -start $is -end $ie -increm $increm -threads $maxthreads -modulo $modulo -imagewidth $imagewidth -fontsize $fontsize

if [ $? != 0 ]; then
  echo Problem running makepng. Exiting.
  exit 1
fi


#make the movie from the png files at 24 frames per second
echo Making the movie $moviefilename from the png files..
if [[ $(ffmpeg -filters 2> /dev/null | grep setpts) ]]; then
  fact=`calc "24/$fps"`
  #ffmpeg -y -pattern_type glob  -i "__makemovie.*.png"  -pix_fmt yuv420p -r 24 -vf scale="trunc(in_w/2)*2:trunc(in_h/2)*2" -filter:v "setpts=$fact*PTS"  $moviefilename
  ffmpeg -y -pattern_type glob  -i "__makemovie.*.png"   -r 24 -filter:v "setpts=$fact*PTS"  $moviefilename

else
  #ffmpeg  -y   -i __makemovie.%8d.png  -pix_fmt yuv420p -r $fps    $moviefilename  
  fact=`calc "round(100/$fps)"`
  convert -quality 100 -delay $fact  __makemovie.*.png $moviefilename
fi




# cleanup
rm -f __makemovie.*.png

#the following can be used to change the frames per second of an existing video
#ffmpeg -y  -i old.mp4 -r $newfps -filter:v "setpts=$fact*PTS" new.mp4

