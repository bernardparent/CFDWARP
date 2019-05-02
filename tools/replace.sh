#!/bin/bash

# SPDX-License-Identifier: BSD-2-Clause

# Copyright 2018 Bernard Parent

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Description: Finds and replaces a string in certain files (.tex, Makefiles, .c) in all directories recursively 


stringfind="";
stringreplace="";
mode="rpl";
filename="";
dryrun="";
errormessage="";
nogit="-path *.git -prune -o";
strtype="any";

while [ $# -gt 0 ]
do
  case "$1" in
    -filename) filename="$2";  shift;;
    -mode) mode="$2";  shift;;
    -find) stringfind="$2"; shift;;
    -replace) stringreplace="$2"; shift;;
    -git) nogit="";;
    -w) strtype="word";;
    -dryrun) dryrun="TEST";;
	--)	shift; break;;
	*) echo  "Argument $1 not recognized.";
          argerror="1";  
          break;;	# terminate while loop
  esac
  shift
done



if [ -z "$mode" ]; then 
  echo "-mode flag not found."
  argerror="1"; 
fi

if [ -z "$filename" ]; then 
  echo "-filename flag not found."
  argerror="1"; 
fi

if [ -z "$stringfind" ]; then 
  echo "-find flag not found."
  argerror="1"; 
fi

if [ -z "$stringreplace" ]; then 
  echo "-replace flag not found."
  argerror="1"; 
fi


if [ -n "$argerror" ]; then
  echo "
FLAG           ARGUMENT                                      Required?
------------------------------------------------------------------------
-filename      '*.c' or '*Makefile' or '*.wrp' or '.config'  Y 
-find          string that will be replaced                  Y
-replace       string that will replace what was found       Y
-mode          rpl or perl (default is rpl)                  N
-git           No argument. If set, will include files       N
               within .git directories.
-w             No argument. If set, will only replace the    N
               string if it is a word.
-dryrun        No argument. If set, will only list the       N 
               files that could be altered without 
               actually doing the find and replace


EXAMPLES

$0 -find 'John' -replace 'Mikey' -filename '*.txt' 
will search within all *.txt files recursively for the string John and replace it by Mikey using rpl.

$0 -git -w -find 'John' -replace 'Mikey' -filename '*Makefile' -mode perl 
will search within all Makefiles recursively for the word John and replace it by Mikey using perl.


SPECIAL CHARACTERS

  ! must be escaped as \!
  ' must be written as '\''  or as '\"'\"'
" 

  exit 1
fi

# check if there are files matching $filename
FILE=`find . $nogit -type f -name "$filename" -print -quit`
if [ -n "$FILE" ]; then
  echo "The following files match $filename:"
  find . $nogit -type f -name "$filename" -print
#  find . -path *.git -prune -o -type f -name "$filename" -print
#  echo "find . $nogit -type f -name \"$filename\" -print"
#  exit 1
  echo " "
else
  echo "ERROR: there are no files matching $filename"
  exit 1
fi

stringfind_orig=$stringfind
stringreplace_orig=$stringreplace

if [[ "$mode" = "perl" ]]; then
  #escape the \
  stringfind=${stringfind//\\/\\\\}
  stringreplace=${stringreplace//\\/\\\\}
  #escape the $ 
  stringfind=${stringfind//$/\\$}
  stringreplace=${stringreplace//$/\\$}
  #escape the / 
  stringfind=${stringfind//\//\\\/}
  stringreplace=${stringreplace//\//\\\/}
fi


if [ -n "$dryrun" ]; then
  
  echo 'Only the above listed files will be subject to the find & replace.'
  echo 'Will search for '$stringfind_orig' and replace with '$stringreplace_orig
  echo $stringfind_orig' is escaped as '$stringfind
  echo $stringreplace_orig' is escaped as '$stringreplace
  exit 1
fi

case "$mode" in
  perl)  
    case "$strtype" in
      any)
        echo "Within files named $filename, replacing string.."
        perl -i.substrperlold -pe 's/'"$stringfind"'/'"$stringreplace"'/g;' \
           `find . $nogit -type f -name "$filename" -print` 
        rm -f `find . $nogit -type f -name '*.substrperlold' -print`
      ;;
      word)
        echo "Within files named $filename, replacing word.."
        perl -i.substrperlold -pe 's/\b'"$stringfind"'\b/'"$stringreplace"'/g;' \
           `find . $nogit -type f -name "$filename" -print` 
        rm -f `find . $nogit -type f -name '*.substrperlold' -print`
      ;;
      *) 
        echo  "strtype can not be set to $strtype."
      break;;
    esac	
  ;;
  rpl) 
    case "$strtype" in
      any)
        echo "Within files named $filename,"
        rpl -b "$stringfind" "$stringreplace" \
          `find . $nogit -type f -name "$filename" -print`
      ;;
      word)
        echo "Within files named $filename,"
        rpl -w -b "$stringfind" "$stringreplace" \
          `find . $nogit -type f -name "$filename" -print`
      ;;
      *) 
        echo  "strtype can not be set to $strtype."
      break;;
    esac	
  ;;
  *) 
    echo  "-mode can not be followed by $mode but should rather be followed by rpl or perl.";
  break;;	
esac



