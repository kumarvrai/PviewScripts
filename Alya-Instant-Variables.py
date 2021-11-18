#### import the simple module from the paraview
import os
import time
import operator
import numpy as np
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

case = sys.argv[1]
model = sys.argv[2]
# INITIALIZE VARIABLES
if(model=='INS'):
 IF = './'+case+'.ensi.case'
 OF = './InsDataExtend.vtm'
 print("--|| ALYA :: READING ALYA ARRAYS FROM",IF)
 startTime = time.time()
 case = OpenDataFile(IF)
 case.UpdatePipeline()
 caseVarNames = case.PointArrays
 # get animation scene
 AS1 = GetAnimationScene()
 AS1.GoToLast()
 print("--|| ALYA :: CALCULATING FLUCT DATA FROM LAST FRAME")
elif(model=='PAVG'):
 IF = './'+case+'.ensi.case'
 OF = './PAvgDataExtend.vtm'
 print("--|| ALYA :: READING ALYA ARRAYS FROM",IF)
 startTime = time.time()
 case = OpenDataFile(IF)
 case.UpdatePipeline()
 caseVarNames = case.PointArrays
 # get animation scene
 AS1 = GetAnimationScene()
 AS1.GoToLast()
 print("--|| ALYA :: CALCULATING PRATIAL AVGD DATA FROM LAST FRAME")
elif(model=='FAVG'):
 varNames = ['AVPRE_average', 'AVVEL_average', 'AVVE2_average', 'AVVXY_average']
 IF = './AvgData_3D.vtm'
 OF = './FAvgDataExtend.vtm'
 print("--|| ALYA :: READING ALYA ARRAYS FROM",IF)
 startTime = time.time()
 case = OpenDataFile(IF)
 case.PointArrayStatus = varNames
 case.UpdatePipeline()
else:
  raise ValueError('--|| ALYA ERROR :: MODEL CAN BE INS/PAVG/FAVG.')


indU = int([i for i, s in enumerate(caseVarNames) if 'AVVEL' in s][0]);
indP = int([i for i, s in enumerate(caseVarNames) if 'AVPRE' in s][0]);
indXX = int([i for i, s in enumerate(caseVarNames) if 'AVVE2' in s][0]);
indXY = int([i for i, s in enumerate(caseVarNames) if 'AVVXY' in s][0]);

print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

if(model=="INS"):
  indu = int([i for i, s in enumerate(caseVarNames) if 'VELOC' in s][0]);
  indp = int([i for i, s in enumerate(caseVarNames) if 'PRESS' in s][0]);
  print("--|| ALYA :: TEMPORAL AVERAGING  ALYA-AVERAGED ARRAYS")
  startTime = time.time()
  case2 = TemporalStatistics(Input=case)
  
  # Properties modified on temporalStatistics1
  case2.ComputeMinimum = 0
  case2.ComputeMaximum = 0
  case2.ComputeStandardDeviation = 0
  case2.UpdatePipeline()

  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA :: APPEND DATASETS")
  startTime = time.time()

  APND1 = AppendAttributes(Input=[case,case2])
  APND1.UpdatePipeline()

  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  # CALCULATE Fluctuations
  print("--|| ALYA :: CALCULATING VELOC FLUCTUATIONS")
  startTime = time.time()
  
  CAL1 = Calculator(Input=APND1)
  CAL1.ResultArrayName = "VFLUC"
  CAL1.Function = " VELOC - AVVEL_average " 
  CAL1.UpdatePipeline()
  
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print("--|| ALYA :: CALCULATING PRESSURE FLUCTUATIONS")
  startTime = time.time()

  CAL1 = Calculator(Input=CAL1)
  CAL1.ResultArrayName = "PFLUC"
  CAL1.Function = " PRESS - AVPRE_average " 
  CAL1.UpdatePipeline()

  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print "--|| ALYA :: APPEND CASE SPECIFIC VARIABLES"
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=CAL1)
  
  PF1.Script = \
  """
  for var in caseVarNames:
    output.PointData.append(inputs[0].PointData[var],var)
  output.PointData.append(inputs[0].PointData["VFLUC"],"VFLUC")
  output.PointData.append(inputs[0].PointData["PFLUC"],"PFLUC")
  
  """
  
  PF1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  # CALCULATE RStresses
  CAL1 = Calculator(Input=PF1)
  CAL1.ResultArrayName = "RS_II"
  CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXX],caseVarNames[indU],'X',caseVarNames[indU],'X',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Y',\
                   caseVarNames[indU],'Z',caseVarNames[indU],'Z')
  CAL1.UpdatePipeline() 

  CAL1 = Calculator(Input=CAL1)
  CAL1.ResultArrayName = "RS_IJ"
  CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXY],caseVarNames[indU],'X',caseVarNames[indU],'Y',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Z',\
                    caseVarNames[indU],'X',caseVarNames[indU],'Z')
  CAL1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print("--|| ALYA :: CALCULATING VELGR, QINST")
  startTime = time.time()
  CAL1 = GradientOfUnstructuredDataSet(Input=CAL1)
  CAL1.ScalarArray = ['POINTS', 'VELOC']
  CAL1.ComputeGradient = 1
  CAL1.ResultArrayName = 'VELGR'
  CAL1.ComputeVorticity = 0
  CAL1.VorticityArrayName = 'OMEG'
  CAL1.ComputeQCriterion = 1
  CAL1.QCriterionArrayName = 'QINST'
  CAL1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  # CALCULATE LAMBDA2
  print("--|| ALYA :: CALCULATING LAMBDA2 - INSTANT")
  startTime = time.time()
  CAL1 = PythonCalculator(Input=CAL1)
  CAL1.ArrayName = "VELAM"
  CAL1.Expression = "eigenvalue(strain(%s)**2 + (VELGR - strain(%s))**2)"% (caseVarNames[indu],caseVarNames[indu])
  CAL1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
elif(model=='PAVG'):
  print("--|| ALYA :: TEMPORAL AVERAGING  ALYA-AVERAGED ARRAYS")
  startTime = time.time()
  case2 = TemporalStatistics(Input=case)
  
  # Properties modified on temporalStatistics1
  case2.ComputeMinimum = 0
  case2.ComputeMaximum = 0
  case2.ComputeStandardDeviation = 0
  case2.UpdatePipeline()

  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA :: APPEND DATASETS")
  startTime = time.time()

  APND1 = AppendAttributes(Input=[case,case2])
  APND1.UpdatePipeline()

  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  # CALCULATE Fluctuations
  print("--|| ALYA :: CALCULATING AVERAGED FLUCTUATIONS")
  startTime = time.time()
  
  CAL1 = Calculator(Input=APND1)
  CAL1.ResultArrayName = "VFLUC"
  CAL1.Function = " AVVEL - AVVEL_average " 
  CAL1.UpdatePipeline()
  
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print("--|| ALYA :: CALCULATING PRESSURE FLUCTUATIONS")
  startTime = time.time()

  CAL1 = Calculator(Input=CAL1)
  CAL1.ResultArrayName = "PFLUC"
  CAL1.Function = " AVPRE - AVPRE_average " 
  CAL1.UpdatePipeline()

  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print "--|| ALYA :: APPEND CASE SPECIFIC VARIABLES"
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=CAL1)
  
  PF1.Script = \
  """
  for var in caseVarNames:
    output.PointData.append(inputs[0].PointData[var],var)
  output.PointData.append(inputs[0].PointData["VFLUC"],"VFLUC")
  output.PointData.append(inputs[0].PointData["PFLUC"],"PFLUC")
  
  """
  
  PF1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  # CALCULATE RStresses
  CAL1 = Calculator(Input=PF1)
  CAL1.ResultArrayName = "RS_II"
  CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXX],caseVarNames[indU],'X',caseVarNames[indU],'X',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Y',\
                   caseVarNames[indU],'Z',caseVarNames[indU],'Z')
  CAL1.UpdatePipeline() 
  CAL1 = Calculator(Input=CAL1)
  CAL1.ResultArrayName = "RS_IJ"
  CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXY],caseVarNames[indU],'X',caseVarNames[indU],'Y',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Z',\
                    caseVarNames[indU],'X',caseVarNames[indU],'Z')
  CAL1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
elif(model=="FAVG"):
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  # CALCULATE RStresses
  CAL1 = Calculator(Input=case)
  CAL1.ResultArrayName = "RS_II"
  CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXX],caseVarNames[indU],'X',caseVarNames[indU],'X',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Y',\
                   caseVarNames[indU],'Z',caseVarNames[indU],'Z')
  CAL1.UpdatePipeline() 
  CAL1 = Calculator(Input=CAL1)
  CAL1.ResultArrayName = "RS_IJ"
  CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXY],caseVarNames[indU],'X',caseVarNames[indU],'Y',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Z',\
                    caseVarNames[indU],'X',caseVarNames[indU],'Z')
  CAL1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
else:
  raise ValueError('--|| ALYA ERROR :: MODEL CAN BE INS/AVG..ADD MORE LATER')

# GRADIENT CALC
print("--|| ALYA :: CALCULATING PRESS GRADIENT")
CAL1 = GradientOfUnstructuredDataSet(Input=CAL1)
startTime = time.time()
CAL1.ScalarArray = ['POINTS', caseVarNames[indP]]
CAL1.ComputeGradient = 1
CAL1.ResultArrayName = 'AVPGR'
CAL1.ComputeVorticity = 0
CAL1.VorticityArrayName = 'OMEGA'
CAL1.ComputeQCriterion = 0
CAL1.QCriterionArrayName = 'QCRIT'
CAL1.UpdatePipeline()
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

print("--|| ALYA :: CALCULATING AVVEL GRADIENT, Q AND VORTICITY")
startTime = time.time()
CAL1 = GradientOfUnstructuredDataSet(Input=CAL1)
CAL1.ScalarArray = ['POINTS', caseVarNames[indU]]
CAL1.ComputeGradient = 1
CAL1.ResultArrayName = 'AVVGR'
CAL1.ComputeVorticity = 1
CAL1.VorticityArrayName = 'OMEGA'
CAL1.ComputeQCriterion = 1
CAL1.QCriterionArrayName = 'QCRIT'
CAL1.UpdatePipeline()
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
# CALCULATE LAMBDA2
print("--|| ALYA :: CALCULATING LAMBDA")
startTime = time.time()
CAL1 = PythonCalculator(Input=CAL1)
CAL1.ArrayName = "LAMDA"
CAL1.Expression = "eigenvalue(strain(%s)**2 + (AVVGR - strain(%s))**2)"% (caseVarNames[indU],caseVarNames[indU])
CAL1.UpdatePipeline()
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

#SAVE THE DATA
print("--|| ALYA :: SAVING DATA")
startTime = time.time()
SaveData(OF, proxy=CAL1,Writetimestepsasfileseries=0)
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
#  SaveData(savePath, proxy=src, Writealltimestepsasfileseries=0,
#  DataMode='Binary', HeaderType='UInt64', EncodeAppendedData=0,
#  CompressorType='None')

print('--||ALYA :: SLICES HAVE BEEN SAVED IN',OF)
