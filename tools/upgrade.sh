#!/bin/sh

# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright 2018 Bernard Parent
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

