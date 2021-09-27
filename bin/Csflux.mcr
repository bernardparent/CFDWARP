#!MC 1410
$!AlterData 
  Equation = '{rhoVCs[0]}={V[0]}*{N}*({chi_Cs}+{chi_Cs+})*2.2e-25'
$!AlterData 
  Equation = '{rhoVCs[1]}={V[1]}*{N}*({chi_Cs}+{chi_Cs+})*2.2e-25'
$!ExtendedCommand 
  CommandProcessorID = 'CFDAnalyzer4'
  Command = 'Integrate VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=28 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'IPlanes\' IntegrateBy=\'Zones\' IRange={MIN =0 MAX = 0 SKIP = 1} JRange={MIN =70 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'
