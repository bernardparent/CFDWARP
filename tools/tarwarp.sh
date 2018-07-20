#!/bin/sh
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
