#!/bin/sh


filecheck=".makefile-header-default"
warppublic="../warp_public"
tmpfile="UHFGHnhKJLJGHGGHKJkljk_tmpfile_from_gitnewversion_public.txt"

command -v git >/dev/null 2>&1 || { echo "ERROR: gitnewversion_public.sh requires git but it is not installed.  Aborting." >&2; exit 1; }


if [ -f "$filecheck" ]; then
  echo "Checking that the current directory is the warp main directory. [OK]";
else
  echo "ERROR: gitnewversion_public.sh must be run from the warp main directory. Exiting.";
  exit 1
fi

if [ $# -eq 1 ]; then
  echo "Checking arguments: $1 specified as new version to commit. [OK]";
else 
  echo "ERROR: gitnewversion_public.sh needs one argument: the new version string. Exiting.";
  exit 1
fi

if [ -n "$(git show-ref --tags $1)" ]; then
  echo "Checking that version $1 is found on main warp git. [OK]"; 
else
  echo "ERROR: version $1 not yet committed. Please commit version $1 to private warp before committing it to the public warp." ;
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
  echo "ERROR: Changes or untracked files reported by git on main warp. Please commit changes to the private warp before committing them to the public warp.";
  exit 1
else
  echo "Checking that there is no changes or untracked files reported by git on main warp. [OK]";
fi


if [ -d "$warppublic" ]; then
  echo "Checking that the $warppublic directory exists. [OK]";
else
  echo "The directory $warppublic does not exist. Cloning it from github.";
  git clone https://bernardparent@github.com/bernardparent/CFDWARP "$warppublic"
  if [ -d "$warppublic" ]; then
    echo "Checking that the  $warppublic directory has been created properly by git. [OK]";
  else
    echo "ERROR: The directory $warppublic does not exist. Exiting.";
    exit 1
  fi
fi


touch "$warppublic/$tmpfile"
if [ -f "$tmpfile" ]; then
  echo "ERROR: The current directory is $warppublic, and not the main warp directory. Exiting.";
  rm -f "$warppublic/$tmpfile"
  exit 1
else
  echo "Checking that the current directory is not $warppublic. [OK]";
  rm -f "$warppublic/$tmpfile"
fi


if [ -d "$warppublic/.git" ]; then
  echo "Checking that the $warppublic/.git directory exists. [OK]";
else
  echo "ERROR: The directory $warppublic/.git does not exist. Exiting.";
  exit 1
fi

echo "Pulling latest public warp from github..";
cd $warppublic
git pull
if [ -n "$(git status --porcelain)" ]; then
#  echo "ERROR: Changes or untracked files reported by git on $warppublic. Can not proceed. Exiting.";
#  exit 1
  echo "Changes or untracked files from previous attempt. This may not be a source of concern. "
else
  echo "Checking that there is no changes or untracked files reported by git on $warppublic. [OK]";
fi
cd -


rm -rf "$warppublic"/* 
rm -f "$warppublic"/.*
cp -a * "$warppublic"
cp .* "$warppublic"
cd "$warppublic"
cd config
chmod u+x removeproprietary.sh
./removeproprietary.sh
chmod u-x removeproprietary.sh
cd ..

echo ' '
echo 'Calling git status to check what changes will be pushed..'
echo ' '
git status
echo ' '
echo -n "Add these changes in new version $1 on PUBLIC GITHUB? (y/N)"
read answer

if [ "$answer" != "${answer#[Yy]}" ] ;then
    echo Yes
else
    echo No
    exit 1
fi


if [ -f "$filecheck" ]; then
  if git show-ref --tags $1 ; then
    echo ERROR: Version $1 already committed. Exiting.
    exit 1
  fi
  if git ls-remote --exit-code --tags origin $1 ; then
    echo ERROR: Version $1 already committed on github. Exiting.
    exit 1
  fi
  echo 'Committing and pushing files to github..'
  git add -A .
  git commit -a -m "$1" 
  git tag -d $1 > /dev/null 2>&1
  git tag -a $1 -m "$1"
  git push --tags origin master
else
 echo "ERROR: couldn't find $filecheck in $warppublic directory. Exiting."
 exit 1
fi


echo '[done]'


