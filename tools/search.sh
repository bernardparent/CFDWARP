#!/bin/sh

# Description: Search for a string in all files in all directories recursively 
# Author: Bernard Parent
# Date: July 2018

wordflag="";
caseflag="";
stringtofind="";
errormessage="";

while [ $# -gt 0 ]
do
  case "$1" in
    -w) wordflag="-w"; shift;;
    -i) caseflag="-i"; shift;;
	 *) if [ -z "$stringtofind" ] 
        then 
          stringtofind=$1; 
          shift; 
        else  
          errormessage="Argument $1 not recognized.";
          break; 
        fi; 	
  esac
done

if [ -z "$errormessage" ] && [ ! -z "$stringtofind" ]
then {
  echo "Searching for string" \'$stringtofind\' "...";
  for i in `find . -path './dev' -prune -o -print 2> /dev/null`
  do {
    fgrep $caseflag $wordflag "$stringtofind" $i > /dev/null 2>&1;
    if [ $? = 0 ]
    then {
      echo $i
    } fi;
  } done;
} else {
  echo ""
  echo "Wrong number of arguments. $errormessage"
  echo ""
  echo "FLAGS:"
  echo "-w specifies that the string must be delimited by word boundaries."
  echo "-i specifies that the string searched is not case sensitive."  
  echo ""
  echo "Eg: search.sh -w -i stringtosearch"
} fi;

