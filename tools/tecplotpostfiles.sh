#!/bin/sh
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
  echo "tecplotpostfiles.sh must be called with 3 arguments:"
  echo " "
  echo "./tecplotpostfiles.sh  post 1 20"
  echo " "
  echo "where post is the name of the tecplot post file excluding the time step suffices"
  echo "  1 is first time step suffix (integer)"
  echo "  20 is the last time step suffix (integer)"      
} fi;

