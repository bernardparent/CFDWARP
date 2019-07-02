#!/bin/sh

# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright 1999 Bernard Parent
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of 
#    conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list 
#    of conditions and the following disclaimer in the documentation and/or other 
#    materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



if [ $# -eq 1 ]; then
  echo "Checking arguments: $1 specified as new version of CFDWARP to use when creating tarball. [OK]";
else 
  echo "ERROR: tarCFDWARP.sh needs one argument: the version string. Exiting.";
  exit 1
fi

if [ -d "CFDWARP" ]; then
  echo "Checking that the CFDWARP directory exists. [OK]";
else
  echo "ERROR: Couldn't find CFDWARP directory. Exiting.";
  exit 1
fi

if [ -f "CFDWARP/.makefile-header-default" ]; then
  echo "Checking that the file CFDWARP/.makefile-header-default exists. [OK]";
else
  echo "ERROR: the file CFDWARP/.makefile-header-default could not be found. Exiting.";
  exit 1
fi



echo 'copying CFDWARP directory to CFDWARP.'$1'...'
cp -a CFDWARP CFDWARP.$1
cd CFDWARP.$1
echo 'cleaning CFDWARP'.$1' directory..'
make cleanall
make cleanbin
cd ..
echo 'creating tarball..'
tar -cpPvzf CFDWARP.$1.tgz CFDWARP.$1
echo 'cleaning up'
rm -rf CFDWARP.$1
echo 'done'
