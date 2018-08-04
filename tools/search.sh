#!/bin/sh

# Description: Search for a string in all files in all directories recursively 
# Author: Bernard Parent
# Date: July 2018

wordflag="";
caseflag="";
stringtofind="";
errormessage="";
nogit="-path *.git -prune -o";
filename='*';

while [ $# -gt 0 ]
do
  case "$1" in
    -w) wordflag="-w";;
    -i) caseflag="-i";;
    -filename) filename="$2";  shift;;
    -git) nogit="";;
	--)	shift; break;;
	 *) if [ -z "$stringtofind" ] 
        then 
          stringtofind=$1;          
        else  
          errormessage="Argument $1 not recognized.";
          break; 
        fi; 	
      ;;
  esac
  shift
done


if [ -z "$stringtofind" ]; then
  errormessage="$errormessage A string to be searched was not provided within the arguments."
fi

if [ -z "$errormessage" ] && [ ! -z "$stringtofind" ]
then {
  echo "Searching for string '$stringtofind' ...";
#  for i in `find . -path './dev' -prune -o -print 2> /dev/null`
  for i in `find . $nogit -type f -name "$filename" -print 2> /dev/null`
  do {
    fgrep $caseflag $wordflag "$stringtofind" $i > /dev/null 2>&1;
    if [ $? = 0 ]
    then {
      echo $i
    } fi;
  } done;
} else {
  echo "
Wrong number of arguments. $errormessage


FLAG           ARGUMENT                                      Required?
------------------------------------------------------------------------
-filename      '*.c' or '*Makefile' or '*.wrp' or '.config'  N 
-git           No argument. If set, will include files       N
               within .git directories.
-w             No argument. If set, will only replace the    N
               string if it is a word.
-i             No argument. If set, will do a case           N
               insensitive search.


EXAMPLES
 
$0 -w 'John' -filename '*.c'
will search within all *.c files recursively for the word John and write out the name of the files that contain such a word

$0 -i -git 'John' 
will search within all files including those within .git directories recursively for the string John (case insensitive) and write out the name of the files that contain such a string


SPECIAL CHARACTERS

  ! must be escaped as \!
  ' must be written as '\''  or as '\"'\"'
";
} fi;

