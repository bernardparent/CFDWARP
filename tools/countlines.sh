#!/bin/sh
if [ $# = 1 ]
then
{
  find . -name "$1" | xargs wc -l
}
else
{
  echo "Wrong number of arguments!"
  echo ""
  echo "$ countlines.sh '*.c'"
  echo ""
  echo "This will count the lines within all *.c recursively."
} fi;
