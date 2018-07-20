#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454
$!VarSet |MFBD| = '/home/parent/pnu/src/warp/bin'
$!ALTERDATA 
  EQUATION = '{Pdiff}=abs({P}[1]-{P}[2])'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer3'
  COMMAND = 'Integrate VariableOption=\'Average\' ScalarVar=11 XVariable=1 YVariable=2 ZVariable=3'
$!DRAWGRAPHICS TRUE
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer3'
  COMMAND = 'Integrate VariableOption=\'Average\' ScalarVar=11 XVariable=1 YVariable=2 ZVariable=3'
$!RemoveVar |MFBD|
