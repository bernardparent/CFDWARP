#!/bin/sh

# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright 2018 Bernard Parent
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



# Description: Search for a string in all files in all directories recursively 

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

