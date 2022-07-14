#!/bin/sh

# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright 2016 Bernard Parent
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



if [ $# = 3 ]
then
{
  rm -f /tmp/_tecplotpostfiles_$USER.sh
  touch /tmp/_tecplotpostfiles_$USER.sh
  echo "#!/bin/sh" > /tmp/_tecplotpostfiles_$USER.sh
  echo "tec360 -mesa \\" >> /tmp/_tecplotpostfiles_$USER.sh
  chmod u+x /tmp/_tecplotpostfiles_$USER.sh

  for i in `seq $2 $3`;
   do
     echo $1.$i \\ >> /tmp/_tecplotpostfiles_$USER.sh
   done    
  /tmp/_tecplotpostfiles_$USER.sh
}   else
{
  echo "${0##*/} must be called with 3 arguments:"
  echo " "
  echo "${0##*/}  post 1 20"
  echo " "
  echo "where post is the name of the tecplot post file excluding the time step suffices"
  echo "  1 is first time step suffix (integer)"
  echo "  20 is the last time step suffix (integer)"      
} fi;

