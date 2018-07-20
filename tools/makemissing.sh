#!/bin/sh



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
$0 -file data. -start 10 -end 30 
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
