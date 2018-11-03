#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454
$!VarSet |MFBD| = '/home/parent/warp/bin'
$!ALTERDATA 
  EQUATION = '{k}=0.5*{rho}*({V[0]}*{V[0]}+{V[1]}*{V[1]}+{V[2]}*{V[2]})'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer3'
  COMMAND = 'Integrate VariableOption=\'Average\' ScalarVar=14 XVariable=1 YVariable=2 ZVariable=3 IRange={MIN =1 MAX = -1 SKIP = 1} JRange={MIN =1 MAX = -1 SKIP = 1} KRange={MIN =1 MAX = -1 SKIP = 1}'
$!RemoveVar |MFBD|
