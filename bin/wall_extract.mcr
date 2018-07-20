#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454
$!VarSet |MFBD| = '/home/parent/pnu/pub/SECOND/testcases/separation'
$!ALTERDATA 
  EQUATION = '{tauxy}=2e-5*ddy({V[0]})'
$!ALTERDATA 
  EQUATION = '{Cf}={tauxy}/(0.5*0.0631*715*715)'
$!ALTERDATA 
  EQUATION = '{qy}=0.03*ddy({T})'
$!DUPLICATEZONE 
  SOURCEZONE = 1
  JRANGE
    {
    MAX = 1
    }
$!RemoveVar |MFBD|
