#!/bin/sh
if [ $# = 2 ]
then
{
  echo "replacing string "$1" with "$2
  perl -i.perlold -pe 's/'$1'/'$2'/g;' \
     `find . -type f -name '*.c' -print` \
     `find . -type f -name '*.h' -print` \
     `find . -type f -name '*.hh' -print` 
  rm -f `find . -type f -name '*.perlold' -print`
}
else
{
  echo "Wrong number of arguments!"
  echo "argument 2 replaces all argument 1 string instances in all files"
  echo "\s to escape a whitespace in first string"
  echo "\( and \) to escape ( and ) characters in first string"
  echo "\[ and \] to escape [ and ] characters in first string"
  echo "\, to escape , in first string"
} fi;
