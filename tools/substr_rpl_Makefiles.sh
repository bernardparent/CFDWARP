#!/bin/sh
if [ $# = 2 ]
then
{
  rpl -v -b "$1" "$2" \
     `find . -path ./.git -prune -o -type f -name '*Makefile' -print`
}
else
{
  echo "Wrong number of arguments!"
  echo ""
  echo "$ substr-Makefile.sh '\$Bernard and \$Michel' 'Bernard and Michel'"
  echo ""
  echo "This will replace \$Bernard and \$Michel string instances by Bernard and Michel in all Makefile files recursively (except .git directory)."
} fi;
