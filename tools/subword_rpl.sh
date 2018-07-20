#!/bin/sh
if [ $# = 2 ]
then
{
  rpl -w -v -b "$1" "$2" \
     `find . -path ./.git -prune -o -type f -name '*.c' -print` \
     `find . -path ./.git -prune -o -type f -name '*.h' -print` \
     `find . -path ./.git -prune -o -type f -name '*.hh' -print` 
}
else
{
  echo "Wrong number of arguments!"
  echo ""
  echo "$ subword.sh 'Bernard' 'Bern'"
  echo ""
  echo "This will replace Bernard word instances by Bern in all .c, .h, .hh files recursively (except .git directory)."
} fi;
