#!/bin/bash


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
$0 -prefix 'vlc.' -suffix '.jpg' -start 1 -end 10 -digits 3" 

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




