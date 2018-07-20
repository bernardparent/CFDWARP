#!/bin/bash



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
        -warp) warpexecutable="$2"; shift;;
        -data) datafilename="$2"; shift;;
	-post)  postfilename="$2"; shift;;
	-threads)  maxthreads="$2"; shift;;
	-digits)  mindigits="$2"; shift;;
	-start)  is="$2"; shift;;
	-end)  ie="$2"; shift;;
	-increm)  increm="$2"; shift;;
	-im)  im="$2"; shift;;
	--)	shift; break;;
	*)  break;;	# terminate while loop
    esac
    shift
done


if [ -z "$warpexecutable" ]; then 
  argerror="1"; 
fi


if [ -z "$datafilename" ]; then 
  argerror="1"; 
fi

if [ -z "$postfilename" ]; then 
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
-warp        warp executable [string]            Y
-data        data file prefix [string]           Y
-post        post file prefix [string]           Y
-start       counter start [int]                 Y
-end         counter end [int]                   Y
-digits      minimum number of digits of post    N
             file counter
-increm      counter increments [int]            N
             (i=start, start+increm, ..)
-threads     max. number of post files made in   N 
             parallel
-im          read data at [int] previous time    N
             levels 


Eg: 
$0 -warp './warp -r riemann.wrp' -data data. -post post. -start 10 -end 30 
will do:
./warp -r riemann.wrp -i out.01.10 -op post.10
./warp -r riemann.wrp -i out.01.11 -op post.11
..
./warp -r riemann.wrp -i out.01.30 -op post.30
and preplot the post files post.10, post.11, ..., post.30" 
  exit 1
fi



# make sure that the warp executable and tecplot's preplot are available
command -v $warpexecutable >/dev/null 2>&1 || { echo "makepost requires $warpexecutable but it is not available.  Aborting." >&2; exit 1; }
command -v preplot >/dev/null 2>&1 || { echo "makepost requires tecplot's preplot but it is not installed.  Aborting." >&2; exit 1; }



i=$is;
thread=0;

while [ $i -le $ie ]
 do
    size=${#i}
    iw="$i"
    for ((m=1; m<=mindigits-$size; m++)); do
      iw="0""$iw"
    done


    im1=$(( i-1 ))	 
    im2=$(( i-2 ))	 
    im3=$(( i-3 ))	 
    if [ -z "$im" ]; then 
      imarg="";
    else
      if [ $im -gt 2 ]; then
        imarg="-im1 $datafilename$im1 -im2 $datafilename$im2 -im3 $datafilename$im3";
      else
        if [ $im -gt 1 ]; then
          imarg="-im1 $datafilename$im1 -im2 $datafilename$im2";
        else
          if [ $im -gt 0 ]; then
            imarg="-im1 $datafilename$im1";
          fi  
        fi 
      fi 
    fi



    if [ ! -f "$datafilename$i" ]
    then
      echo File $datafilename$i not found. Skipping.
    else
      echo $warpexecutable -i $datafilename$i $imarg -op $postfilename$iw 
      $warpexecutable -i $datafilename$i $imarg -op $postfilename$iw 1> /dev/null  && preplot $postfilename$iw $postfilename$iw.plt 1> /dev/null  && rm -f $postfilename$iw >/dev/null 2>&1 && mv $postfilename$iw.plt $postfilename$iw >/dev/null 2>&1 &
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
