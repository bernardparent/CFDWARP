#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454
$!VarSet |MFBD| = '/home/parent/pnu/src/warp/bin'
$!ALTERDATA 
  EQUATION = '{rhoudiff}=abs({rho}[1]*{V[0]}[1]-{rho}[2]*{V[0]}[2])'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer3'
  COMMAND = 'Integrate VariableOption=\'Average\' ScalarVar=11 XVariable=1 YVariable=2 ZVariable=3'
$!RemoveVar |MFBD|
