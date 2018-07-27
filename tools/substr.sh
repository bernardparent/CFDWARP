#!/bin/sh

# Description: Substitutes a string in certain files (.tex, Makefiles, .c) in all directories recursively 
# Author: Bernard Parent
# Date: July 2018

stringtofind="";
stringtoreplace="";
mode="";
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
    -find) stringtofind="$2"; shift;;
    -replace) stringtoreplace="$2"; shift;;
    -git) nogit="";;
    -word) strtype="word";;
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

if [ -z "$stringtofind" ]; then 
  echo "-find flag not found."
  argerror="1"; 
fi

if [ -z "$stringtoreplace" ]; then 
  echo "-replace flag not found."
  argerror="1"; 
fi


if [ -n "$argerror" ]; then
  echo "
Flag           Arg                                           Required?
------------------------------------------------------------------------
-filename      '*.c' or '*Makefile' or '*.wrp' or '.config'  Y 
-find          string that should be found                   Y
-replace       string that will replace what was found       Y
-mode          rpl or perl                                   Y
-git           No argument. If set, will include files       N
               within .git directories.
-word          No argument. If set, will only replace the    N
               string if it is a word.
-dryrun        No argument. If set, will only list the      N 
               files that could be altered without 
               actually doing the find and replace

Eg: 
$0 -find 'Bernie' -replace 'Mikey' -filename '*.txt' -mode perl
will search within all *.txt files recursively for the string Bernie and replace it by Mikey using perl.

When using perl:
  "'\\'" must be escaped as "'\\''\\'" 
  ! must be escaped as \!
  \$ must be escaped as \\\$
  / must be escaped as \\/
  ' must be written as '\''  or '\"'\"'

When using rpl:
  ! must be escaped as \!
  ' must be written as '\''  

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

if [ -n "$dryrun" ]; then
  
  echo 'Only the above listed files will be subject to the find & replace.'
  exit 1
fi

case "$mode" in
  perl)  
    case "$strtype" in
      any)
        echo "Within files named $filename, replacing string.."
        perl -i.substrperlold -pe 's/'"$stringtofind"'/'"$stringtoreplace"'/g;' \
           `find . $nogit -type f -name "$filename" -print` 
        rm -f `find . $nogit -type f -name '*.substrperlold' -print`
      ;;
      word)
        echo "Within files named $filename, replacing word.."
        perl -i.substrperlold -pe 's/\b'"$stringtofind"'\b/'"$stringtoreplace"'/g;' \
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
        rpl -b "$stringtofind" "$stringtoreplace" \
          `find . $nogit -type f -name "$filename" -print`
      ;;
      word)
        echo "Within files named $filename,"
        rpl -w -b "$stringtofind" "$stringtoreplace" \
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



