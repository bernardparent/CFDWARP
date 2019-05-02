#!/bin/bash


# SPDX-License-Identifier: BSD-2-Clause

# Copyright 2019 Prasanna T. Rajendran

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


increm=1;

while [ $# -gt 0 ]
do
    case "$1" in
	-inputfileprefix)  fileprefix="$2"; shift;;
	-output)  output="$2"; shift;;
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

if [ -z "$fileprefix" ]; then 
  argerror="1"; 
fi

if [ -z "$is" ]; then 
  argerror="1"; 
fi

if [ -z "$ie" ]; then 
  argerror="1"; 
fi

declare -a avg_array;

if [ -n "$argerror" ]; then
  echo "
Flag         		Arg                                 Required?
---------------------------------------------------------------
-inputfileprefix  	input file prefix [string]           Y
-output      		output text file [string]            Y
-start       		counter start [int]                  Y
-end         		counter end [int]                    Y
-increm      		counter increments [int]             N
             		(i=start, start+increm, ..)

Eg: 
$0 -inputfileprefix Cf_x. -start 10 -end 30 -output avg_out.txt
will average the values in all columns in the specified files (Cf_x.10,Cf_x.11...Cf_x.30) and output the average to the file avg_out.txt ."
 
  exit 1
fi

command -v sed >/dev/null 2>&1 || { echo "average.sh requires sed but it is not installed.  Aborting." >&2; exit 1; }
command -v calc >/dev/null 2>&1 || { echo "average.sh requires calc but it is not installed.  Aborting." >&2; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "average.sh requires awk but it is not installed.  Aborting." >&2; exit 1; }
command -v sort >/dev/null 2>&1 || { echo "average.sh requires sort but it is not installed.  Aborting." >&2; exit 1; }
command -v tail >/dev/null 2>&1 || { echo "average.sh requires tail but it is not installed.  Aborting." >&2; exit 1; }

if [ ! -f "$fileprefix$is" ]
then
  echo "$fileprefix$is. Make sure atleast $fileprefix$is exists. The number of rows and columns are obtained from $fileprefix$is and checked with the rest of the files."; exit 1;
fi

numrows=$(sed -n '$=' $fileprefix$is);
numcols=$(awk '{print NF}' $fileprefix$is | sort -nu | tail -n 1);

numfiles=0;
givennumfiles=0;
numelements=$(awk '{c+=NF} END {print c+0}' $fileprefix$is);

for ((filenum=$is; filenum<=$ie; filenum=filenum+$increm))
do
	givennumfiles=$((givennumfiles+1));
done

echo "Checking for missing files."
for ((filenum=$is; filenum<=$ie; filenum=filenum+$increm))
do
	if [ ! -f "$fileprefix$filenum" ]
    then
      echo File $fileprefix$filenum not found. Skipping.
    else
	  numfiles=$((numfiles+1));
	  if [ $numrows != $(sed -n '$=' $fileprefix$filenum) ] || [ $numcols != $(awk '{print NF}' $fileprefix$filenum | sort -nu | tail -n 1) ] || [ $numelements != $(awk '{c+=NF} END {print c+0}' $fileprefix$filenum) ]
	  then
	    echo "Mismatch in number of rows/columns or missing value at a row/column in file $fileprefix$filenum, when comparing with the first file $fileprefix$is. Exiting."; exit 1;
	  fi
    fi
done
echo "$numfiles/$givennumfiles files found in total."

for ((colnum=1; colnum<=$numcols; colnum++))
do
	for ((i=1; i<=$numrows; i++))
	do
		avg_array[$i]=0.0;
	done

	for ((filenum=$is; filenum<=$ie; filenum=filenum+$increm))
	do
		if [ -f "$fileprefix$filenum" ]
		  then
		  for ((i=1; i<=$numrows; i++))
		  do
			avg_array[$i]=$(calc "$(awk -v iloop="$i" -v numcol_loop="$colnum" 'NR == iloop {print $numcol_loop}' $fileprefix$filenum) + ${avg_array[$i]}");
		  done
		  echo "Processing column number $colnum in $fileprefix$filenum";
		fi
	done

	for ((i=1; i<=$numrows; i++))
	do
		avg_array[$i]=$(calc "${avg_array[$i]}/($numfiles)");
	done
	printf "%s\n" "${avg_array[@]}" > $output$colnum
	file_names="$file_names $output$colnum"
done

paste $file_names > $output
rm -f $file_names
echo "Averaged over $numfiles files. Results outputted to $output. Done."
