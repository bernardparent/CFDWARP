#!/bin/sh

# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright 2016-2018 Bernard Parent
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


command -v git >/dev/null 2>&1 || { echo "gitnewversion requires git but it is not installed.  Aborting." >&2; exit 1; }

filecheck=".makefile-header-default"
if [ -f "$filecheck" ]; then
 if [ $# -eq 1 ]; then
  if git show-ref --tags $1 ; then
    echo Version $1 already committed
    exit 
  fi
  if git ls-remote --exit-code --tags origin $1 ; then
    echo Version $1 already committed on origin
    exit 
  fi
  make cleanall
  make cleanbin
  echo ' '
  echo 'Calling git status to check what changes will be pushed..'
  echo ' '
  git status
  echo ' '
  echo -n "Add these changes in new version $1 on PRIVATE SERVER? (y/N)"
  read answer

  if [ "$answer" != "${answer#[Yy]}" ] ;then
    echo Yes
  else
    echo No
    exit 1
  fi


  echo 'setting new version in src/version.hh'
  echo '#define VERSION "'"$1"'"' > src/version.hh
  echo 'Committing and pushing files to private server..'
  git add -A .
  git commit -a -m "$1" 
  git tag -d $1 > /dev/null 2>&1
  git tag -a $1 -m "$1"
  git push --tags origin master


  tools/gitnewversion_public.sh $1

 else
  echo "gitnewversion.sh needs one argument: the new version string"
 fi
else
 echo "gitnewversion.sh must be run from the warp main directory."
fi


