#!/bin/sh
if [ $# = 2 ]
then
{
  rpl -v -b "$1" "$2" \
     `find . -path ./.git -prune -o -type f -name '.config' -print`
}
else
{
  echo "Wrong number of arguments!"
  echo ""
  echo "$ substr-texfiles.sh '\$Bernard and \$Michel' 'Bernard and Michel'"
  echo ""
  echo "This will replace \$Bernard and \$Michel string instances by Bernard and Michel in all .config files recursively (except .git directory)."
} fi;
