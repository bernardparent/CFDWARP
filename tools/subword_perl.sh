#!/bin/sh
if [ $# = 2 ]
then
{
  perl -i.perlold -pe 's/\b'$1'\b/'$2'/g;' \
     `find . -type f -name '*.c' -print` \
     `find . -type f -name '*.h' -print` \
     `find . -type f -name '*.hh' -print` 
  rm -f `find . -type f -name '*.perlold' -print`
}
else
{
  echo "Wrong number of arguments!"
  echo "argument 2 replaces all argument 1 word instances in all files"
} fi;
