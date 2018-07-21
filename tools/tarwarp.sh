#!/bin/sh

if [ $# -eq 1 ]; then
  echo "Checking arguments: $1 specified as new version of warp to use when creating tarball. [OK]";
else 
  echo "ERROR: tarwarp.sh needs one argument: the version string. Exiting.";
  exit 1
fi

if [ -d "warp" ]; then
  echo "Checking that the warp directory exists. [OK]";
else
  echo "ERROR: Couldn't find warp directory. Exiting.";
  exit 1
fi

if [ -f "warp/.makefile-header-default" ]; then
  echo "Checking that the file warp/.makefile-header-default exists. [OK]";
else
  echo "ERROR: the file warp/.makefile-header-default could not be found. Exiting.";
  exit 1
fi



echo 'copying warp directory to warp.'$1'...'
cp -a warp warp.$1
cd warp.$1
echo 'cleaning warp'.$1' directory..'
make cleanall
make cleanbin
cd ..
echo 'creating tarball..'
tar -cpPvzf warp.$1.tgz warp.$1
echo 'cleaning up'
rm -rf warp.$1
echo 'done'
