#!/bin/sh


filecheck=".makefile-header-default"
warpopen="../warpopen"
tmpfile="UHFGHnhKJLJGHGGHKJkljk_tmpfile_from_gitnewversion_open.txt"

command -v git >/dev/null 2>&1 || { echo "gitnewversion_open.sh requires git but it is not installed.  Aborting." >&2; exit 1; }


if [ -f "$filecheck" ]; then
  echo "Checking that the current directory is the warp main directory. [OK]";
else
  echo "ERROR: gitnewversion_open.sh must be run from the warp main directory. Exiting.";
  exit 1
fi

if [ $# -eq 1 ]; then
  echo "Checking arguments: $1 specified as new version to commit. [OK]";
else 
  echo "ERROR: gitnewversion_open.sh needs one argument: the new version string. Exiting.";
  exit 1
fi

if [ -n "$(git show-ref --tags $1)" ]; then
  echo "Checking that version $1 is found on main warp git. [OK]"; 
else
  echo "ERROR: version $1 not yet committed. Please commit version $1 to main warp origin before committing it to warpopen." ;
  exit 1
fi

latesttag=$(git describe --tags `git rev-list --tags --max-count=1`);

if [ "$1" = "$latesttag" ]; then
  echo "Checking that the latest tag is $1 on the main warp git. [OK]";
else
  echo "ERROR: The tag given $1 is not set to the latest tag $latesttag on the main warp git. Exiting.";
  exit 1
fi


if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: Changes or untracked files reported by git on main warp. Please commit changes to main warp origin before committing them to warpopen.";
  exit 1
else
  echo "Checking that there is no changes or untracked files reported by git on main warp. [OK]";
fi


if [ -d "$warpopen" ]; then
  echo "Checking that the $warpopen directory exists. [OK]";
else
  echo "The directory $warpopen does not exist. Cloning it from github.";
  git clone https://bernardparent@github.com/bernardparent/CFDWARP "$warpopen"
  if [ -d "$warpopen" ]; then
    echo "Checking that the  $warpopen directory has been created properly by git. [OK]";
  else
    echo "ERROR: The directory $warpopen does not exist. Exiting.";
    exit 1
  fi
fi


touch "$warpopen/$tmpfile"
if [ -f "$tmpfile" ]; then
  echo "ERROR: The current directory is $warpopen, and not the main warp directory. Exiting.";
  rm -f "$warpopen/$tmpfile"
  exit 1
else
  echo "Checking that the current directory is not $warpopen. [OK]";
  rm -f "$warpopen/$tmpfile"
fi


if [ -d "$warpopen/.git" ]; then
  echo "Checking that the $warpopen/.git directory exists. [OK]";
else
  echo "ERROR: The directory $warpopen/.git does not exist. Exiting.";
  exit 1
fi

echo "Pulling latest warpopen from github..";
cd $warpopen
git pull
if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: Changes or untracked files reported by git on warpopen. Can not proceed. Exiting.";
  exit 1
else
  echo "Checking that there is no changes or untracked files reported by git on warpopen. [OK]";
fi
cd -

echo -n "Create new version $1 on github? (y/n)"
read answer

if [ "$answer" != "${answer#[Yy]}" ] ;then
    echo Yes
else
    echo No
    exit 1
fi

rm -rf "$warpopen"/* 
rm -f "$warpopen"/.*
cp -a * "$warpopen"
cp .* "$warpopen"
cd "$warpopen"
cd config
chmod u+x removeproprietary.sh
./removeproprietary.sh
chmod u-x removeproprietary.sh
cd ..


if [ -f "$filecheck" ]; then
  if git show-ref --tags $1 ; then
    echo ERROR: Version $1 already committed. Exiting.
    exit 1
  fi
  if git ls-remote --exit-code --tags origin $1 ; then
    echo ERROR: Version $1 already committed on github. Exiting.
    exit 1
  fi
  echo 'copy to origin using git..'
  git add -A .
  git commit -a -m "$1" 
  git tag -d $1 > /dev/null 2>&1
  git tag -a $1 -m "$1"
  git push --tags origin master
  echo '[done]'
else
 echo "ERROR: couldn't find $filecheck in warpopen directory. Exiting."
 exit 1
fi


echo '[done]'


