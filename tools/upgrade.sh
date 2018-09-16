#!/bin/sh
command -v git >/dev/null 2>&1 || { echo "upgrade.sh requires git but it is not installed.  Aborting." >&2; exit 1; }

filecheck=".makefile-header-default"
if [ -f "$filecheck" ]; then

  echo 
  echo -n "All changes made to files within CFDWARP will be lost.\n\nAre you sure you want to upgrade? (y/N)"
  read answer

  if [ "$answer" != "${answer#[Yy]}" ] ;then
    echo Yes
  else
    echo No
    exit 1
  fi
 
  git checkout master
  git reset --hard
  git pull --tags origin master
else
 echo "upgrade.sh must be run from the CFDWARP main directory."
fi

