#!/bin/sh
if [ $# = 2 ]
then
{
  rpl -v -b "$1" "$2" \
     `find . -path ./.git -prune -o -type f -name '*.c' -print` \
     `find . -path ./.git -prune -o -type f -name '*.h' -print` \
     `find . -path ./.git -prune -o -type f -name '*.hh' -print` 
}
else
{
  echo "Wrong number of arguments!"
  echo ""
  echo "$ substr.sh '\$Bernard and \$Michel' 'Bernard and Michel'"
  echo ""
  echo "This will replace \$Bernard and \$Michel string instances by Bernard and Michel in all .c, .h, .hh files recursively (except .git directory)."
} fi;
