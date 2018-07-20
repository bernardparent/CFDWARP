#!/bin/sh


filecheck=".makefile-header-default"
warpopen="../warpopen"
tmpfile="UHFGHnhKJLJGHGGHKJkljk_tmpfile_from_gitnewversion_open.txt"


if [ -d "$warpopen" ]; then
  echo "Found $warpopen directory. [OK]";
else
  echo "ERROR: The directory $warpopen does not exist. Exiting.";
  exit 1
fi

if [ -d "$warpopen/.git" ]; then
  echo "Found $warpopen/.git directory. [OK]";
else
  echo "ERROR: The directory $warpopen/.git does not exist. Exiting.";
  exit 1
fi

touch "$warpopen/$tmpfile"
if [ -f "$tmpfile" ]; then
  echo "ERROR: The current directory is $warpopen, and not the main warp directory. Exiting.";
  rm -f "$warpopen/$tmpfile"
  exit 1
else
  echo "The current directory is not $warpopen. [OK]";
fi

if [ -f "$filecheck" ]; then
  echo "The current directory is the warp main directory. [OK]";
else
  echo "ERROR: gitnewversion_open.sh must be run from the warp main directory. Exiting.";
  exit 1
fi

if [ $# -eq 1 ]; then
  echo "$1 specified as new version to commit. [OK]";
else 
  echo "ERROR: gitnewversion_open.sh needs one argument: the new version string. Exiting.";
fi

if [ -n "$(git show-ref --tags $1)" ]; then
  echo "Version $1 found on main warp git. [OK]"; 
else
  echo "ERROR: version $1 not yet committed. Please commit version $1 to main warp origin before committing it to warpopen." ;
  exit 1
fi


if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: Changes or untracked files reported by git on main warp. Please commit changes to main warp origin before committing them to warpopen.";
  exit 1
else
  echo "No changes or untracked files reported by git on main warp. [OK]";
fi




rm -rf "$warpopen"/*
rm -f "$warpopen"/.*
cp -a * "$warpopen"
cp .* "$warpopen"
cd "$warpopen"
make config
cd config
chmod u+x removeproprietary.sh
./removeproprietary.sh
chmod u-x removeproprietary.sh
cd ..
./tools/gitnewversion.sh $1
echo '[done]'


