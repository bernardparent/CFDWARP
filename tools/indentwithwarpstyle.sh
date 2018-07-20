#!/bin/sh
if [ $# = 1 ]
then
{
  indent -nut -nv -gnu -l110 -br -brf -brs -ce -cdw -sob -npsl -nbfda -prs  $1
}
else
{
  echo "Wrong number of arguments!"
  echo ""
  echo "To change the format of test.c to the WARP standard format:"
  echo "$ indentwithwarpstyle.sh test.c"
  echo ""
  echo "To change the format of all */*.c files to the WARP standard format:"
  echo "$ indentwithwarpstyle.sh '*/*.c'"
} fi;
