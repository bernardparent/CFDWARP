#!/bin/sh
if [ $# = 2 ]
then
{
  perl -i.perlold -pe 's/'$1'/'$2'/g;' \
     `find . -type f -name '*.wrp' -print` 
  rm -f `find . -type f -name '*.perlold' -print`
}
else
{
  echo "Wrong number of arguments!"
  echo "argument 2 replaces all argument 1 string instances in all files"
} fi;
