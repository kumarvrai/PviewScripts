#### import the simple module from the paraview
import os, argparse
import glob
import time
import operator
import numpy as np
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *

#pxm  = servermanager.ProxyManager()
#pxm.GetVersion()
#print("--|| NEK :: USING PARAVIEW VERSION",pxm)

#------------------------------------#
def convert_to_float(frac_str):
    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac
#------------------------------------#


caseName	= sys.argv[1]
codeName	= sys.argv[2]
fileType	= sys.argv[3]
qVal		= float(sys.argv[4])
ifTSeries		= int(sys.argv[5])

zDec = 6; xDec = 6
xc = 1.0; yc = 0.0; 
casePath = os.getcwd()


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

if('NEK' in codeName):
  print("--|| NEK :: READING NEK5000 ARRAYS")
  startTime = time.time()
  
  fileName1 = caseName+'.nek5000'
  fileName2 = 'qcrit'+caseName+'.nek5000'
  
  case1 = OpenDataFile(fileName1)
  case1.PointArrays = ['velocity','pressure','temperature']
  case1.UpdatePipeline()
  
  ## create a new 'Programmable Filter and change names'
  print("--|| NEK: CHANGING VARNAMES 1 USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  case1 = ProgrammableFilter(Input=case1)
  case1.Script = \
  """
  import numpy as np
  varNames0 = ['velocity','pressure','temperature']
  varNames1 = ['VELOC','PRESS','LAMDA']
  for (i,var) in enumerate(varNames0):
   outName = varNames1[i]
   avg = (inputs[0].PointData[var])
   output.PointData.append(avg,outName)
  """
  case1.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  case2 = OpenDataFile(fileName2)
  case2.PointArrays = ['pressure']
  case2.UpdatePipeline()
  
  ## create a new 'Programmable Filter and change names'
  print("--|| NEK: CHANGING VARNAMES 2 USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  case2 = ProgrammableFilter(Input=case2)
  case2.Script = \
  """
  import numpy as np
  varNames0 = ['pressure']
  varNames1 = ['QCRIT']
  for (i,var) in enumerate(varNames0):
   outName = varNames1[i]
   avg = (inputs[0].PointData[var])
   output.PointData.append(avg,outName)
  """
  case2.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  ##case3 = OpenDataFile(fileName3)
  ##case3.PointArrays = ['velocity']
  ##case3.UpdatePipeline()
  ##
  #### create a new 'Programmable Filter and change names'
  ##print("--|| NEK: CHANGING VARNAMES 3 USING A PROGRAMMABLE FILTER")
  ##startTime = time.time()
  ##case3 = ProgrammableFilter(Input=case3)
  ##case3.Script = \
  ##"""
  ##import numpy as np
  ##varNames0 = ['velocity']
  ##varNames1 = ['AVVXY']
  ##for (i,var) in enumerate(varNames0):
  ## outName = varNames1[i]
  ## avg = (inputs[0].PointData[var])
  ## output.PointData.append(avg,outName)
  ##"""
  ##case3.UpdatePipeline()
  ##print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  print("--|| NEK: APPEND DATASETS")
  startTime = time.time()
  #case = AppendAttributes(Input=[case1,case2,case3])
  case = AppendAttributes(Input=[case1,case2])
  case.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
elif('OFOAM' in codeName):
  print("--|| OFOAM :: READING OPENFOAM ARRAYS")
  startTime = time.time()
  
  fileName = caseName+'.foam'
  
  case = OpenDataFile(fileName)
  if('INS' in fileType):
    case.CellArrays = ['p','U']
    if('PAR' in fileType):
      case.CaseType = ['Decomposed Case']
    case.UpdatePipeline()
    print("--|| NEK: CHANGING VARNAMES USING A PROGRAMMABLE FILTER")
    startTime = time.time()
    case = ProgrammableFilter(Input=case)
    case.Script = \
    """
    import numpy as np
    varNames0 = ['p','U']
    varNames1 = ['AVPRE','AVVEL']
    for (i,var) in enumerate(varNames0):
     outName = varNames1[i]
     avg = (inputs[0].CellData[var])
     output.CellData.append(avg,outName)
    """
    case.UpdatePipeline()
    print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  elif('AVG' in fileType):
    case.CellArrays = ['pMean','UMean','pPrime2Mean','UPrime2Mean']
    if('PAR' in fileType):
      case.caseType = ['Decomposed Case']
    case.UpdatePipeline()
    print("--|| NEK: CHANGING VARNAMES USING A PROGRAMMABLE FILTER")
    startTime = time.time()
    case = ProgrammableFilter(Input=case)
    case.Script = \
    """
    import numpy as np
    varNames0 = ['pMean','UMean','pPrime2Mean','UPrime2Mean']
    varNames1 = ['AVPRE','AVVEL','AVPR2','RESTR']
    for (i,var) in enumerate(varNames0):
     outName = varNames1[i]
     avg = (inputs[0].CellData[var])
     output.CellData.append(avg,outName)
    """
    case.UpdatePipeline()
    print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  else:
    raise ValueError('--|| ALYA ERROR :: FILETYPE NOT RECONIZED.')
  case.UpdatePipeline()

  print("--|| NEK :: CONVERT CELL DATA TO POINT DATA.")
  startTime = time.time()
  case = CellDatatoPointData(Input=case)
  case.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
elif("ALYA" in codeName):
  print("--|| INFO :: READING ALYA ARRAYS")
  startTime = time.time()
  fileName = caseName+'.ensi.case'
  case = OpenDataFile(fileName)
  if(fileType == "LES"):
    case.PointArrays = ['TURBU','YPLUS','AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
  elif(fileType == "DNS"):
    case.PointArrays = ['AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
  else:
    case.PointArrays = ['YPLUS', 'AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
  case.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
elif("SOD" in codeName):
  print("--|| INFO :: READING SOD2D ARRAYS")
  startTime = time.time()
  if('AVG' in fileType):
    fileName = 'results_AVG_'+caseName+'.hdf'
    case = OpenDataFile(fileName)
    case.UpdatePipeline()
    print("--||SOD :: LOADED VARIABLES",case.PointData.keys())
    ## create a new 'Programmable Filter and change names'
    print("--|| SOD: CHANGING VARNAMES USING A PROGRAMMABLE FILTER")
    startTime = time.time()
    case = ProgrammableFilter(Input=case)
    case.Script = \
    """
    import numpy as np
    varNames0 = inputs[0].PointData.keys()
    if("avrho" in varNames0):
      rho = inputs[0].PointData["avrho"]
    #---------------------------------#
    for (i,var) in enumerate(varNames0):
     avg = inputs[0].PointData[var]
     outName = var.upper()
     if(outName in ["AVVEL","AVVE2","AVVEX"]):
       if("avrho" in varNames0):
         avg = avg/rho
     if("AVVEX" in outName):
         outName = "AVVXY"
     output.PointData.append(avg,outName)
    """
  elif('INS' in fileType):
    files = './results_'+caseName+'.hdf'
    fileList = glob.glob(files, recursive=True)
    fileList.sort(key=os.path.getmtime, reverse=False)
    fileName = []
    for (i,file) in enumerate(fileList):
      print(file)
      fileName.append(file)
    print(fileName)
    case = OpenDataFile(fileName)
    case.UpdatePipeline()
    print("--||SOD :: LOADED VARIABLES",case.PointData.keys())
    ## create a new 'Programmable Filter and change names'
    print("--|| SOD: CHANGING VARNAMES USING A PROGRAMMABLE FILTER")
    startTime = time.time()
    case = ProgrammableFilter(Input=case)
    case.Script = \
    """
    import numpy as np
    varNames0 = ['u','pr','qcrit']
    varNames0 = ['curlU', 'eta', 'mu_fluid', 'mue', 'mut', 'pr', 'qcrit', 'u']
    varNames1 = ['CRULU', 'ETAFL', 'MUFLD', 'MUENT', 'MUTUR', 'PRESS', 'QCRIT', 'VELOC']
    #---------------------------------#
    for (i,var) in enumerate(varNames0):
     outName = varNames1[i]
     avg = inputs[0].PointData[var]
     output.PointData.append(avg,outName)
    """
  case.UpdatePipeline()
  print("--|| SOD :: DONE. TIME =",time.time()-startTime,'sec')
elif("PVD" in codeName):
  print("--|| INFO :: READING PVD ARRAYS")
  startTime = time.time()
  fileName = caseName+'.pvd'
  case = OpenDataFile(fileName)
  case.UpdatePipeline()
  print("--|| PVD :: DONE. TIME =",time.time()-startTime,'sec')
elif("VTM" in codeName):
  print("--|| INFO :: READING PVD ARRAYS")
  startTime = time.time()
  fileName = caseName+'.vtm'
  case = OpenDataFile(fileName)
  case.UpdatePipeline()
  print("--|| VTM :: DONE. TIME =",time.time()-startTime,'sec')
else:      
  raise ValueError('--|| ALYA ERROR :: CODENAME NOT RECONIZED.')

if('CLEAN' in codeName):
  print("--|| ALYA :: APPLYING CLEAN TO GRID FILTER")
  startTime = time.time()
  case = CleantoGrid(Input=case)
  case.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

Ntotal = int(case.GetDataInformation().GetNumberOfPoints())
print("----|| ALYA :: WORKING WITH ",Ntotal," TOTAL NUMBER OF POINTS")

(xmin,xmax,ymin,ymax,zmin,zmax) =  case.GetDataInformation().GetBounds()
print("----|| INFO: BOX SIZE = %.2f %.2f %.2f %.2f %.2f %.2f"%(xmin,xmax,ymin,ymax,zmin,zmax))

caseVarNames = case.PointData.keys()
if("SKIP" not in codeName):
  indU = int([i for i, s in enumerate(caseVarNames) if 'VELOC' in s][0]);
  indP = int([i for i, s in enumerate(caseVarNames) if 'PRESS' in s][0]);
  # GRADIENT CALC
  print("--|| NEK :: CALCULATING PRESS GRADIENT")
  startTime = time.time()
  case = GradientOfUnstructuredDataSet(Input=case)
  case.ScalarArray = ['POINTS', caseVarNames[indP]]
  case.ComputeGradient = 1
  case.ResultArrayName = 'AVPGR'
  case.ComputeVorticity = 0
  case.VorticityArrayName = 'OMEGA'
  case.ComputeQCriterion = 0
  case.QCriterionArrayName = 'QCRIT'
  case.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| NEK :: CALCULATING AVVEL GRADIENT, Q AND VORTICITY")
  startTime = time.time()
  case = GradientOfUnstructuredDataSet(Input=case)
  case.ScalarArray = ['POINTS', caseVarNames[indU]]
  case.ComputeGradient = 1
  case.ResultArrayName = 'AVVGR'
  case.ComputeVorticity = 1
  case.VorticityArrayName = 'OMEGA'
  if('SOD' not in codeName):
    case.ComputeQCriterion = 1
  else:  
    case.ComputeQCriterion = 0
  case.QCriterionArrayName = 'QCRIT'
  case.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
if("CPCF" in fileType):
  # CALCULATE CPCF
  print("--|| NEK :: CALCULATING CPCF")
  case_clcd = ExtractSurface(Input=case)
  case_clcd.UpdatePipeline()
  case_clcd = GenerateSurfaceNormals(Input=case_clcd)
  case_clcd.UpdatePipeline()
  SetActiveSource(case_clcd)            
  Show(case_clcd)
  SetActiveView(GetActiveView())
  QuerySelect(QueryString='(mag(AVVEL) == 0)', FieldType='POINT', InsideOut=0)
  case_clcd = ExtractSelection(Input=case_clcd)
  case_clcd.UpdatePipeline()
  #case_clcd = Clip(Input=case_clcd)
  #case_clcd.ClipType = 'Box'
  #case_clcd.ClipType.Position = [-1, -1, 0]
  #case_clcd.ClipType.Length = [3, 3, 1]
  #case_clcd.Invert = 1
  #case_clcd.UpdatePipeline()
  case_clcd = Calculator(Input=case_clcd)
  case_clcd.ResultArrayName = "AVGCF"
  #case_clcd.Function = "%s*sqrt(dot((""AVVGR_0""*iHat+""AVVGR_1""*jHat),-Normals)^2+dot((""AVVGR_3""*iHat+""AVVGR_4""*jHat),-Normals)^2)"%visc
  case_clcd.Function = "2.0*%s*(dot((""AVVGR_0""*iHat+""AVVGR_1""*jHat),-Normals)+dot((""AVVGR_3""*iHat+""AVVGR_4""*jHat),-Normals))"%visc
  case_clcd.UpdatePipeline()
  case_clcd = Calculator(Input=case_clcd)
  case_clcd.ResultArrayName = "AVGCP"
  case_clcd.Function = "2.0*AVPRE"
  case_clcd.UpdatePipeline()
  savePath = casePath+"/NekBndData_1D.csv"
  #SaveData(savePath, proxy=case_clcd)
  SaveData(savePath, proxy=case_clcd,ChooseArraysToWrite=1,
                     PointDataArrays=['AVGCF', 'AVGCP'],
                     UseScientificNotation=1)
  print("----|| NEK: 1D CPCF FILE WRITTEN AS: ",savePath)
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')

## CALCULATE LAMBDA2
#print("--|| NEK :: CALCULATING LAMBDA")
#startTime = time.time()
#case = PythonCalculator(Input=case)
#case.ArrayName = "LAMDA"
#case.Expression = "eigenvalue(strain(%s)**2 + (AVVGR - strain(%s))**2)"% (caseVarNames[indU],caseVarNames[indU])
#case.UpdatePipeline()
#print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')

# create a new 'Contour'
cal = Calculator(Input=case)
cal.ResultArrayName = "AVVMG"
cal.Function = "mag(VELOC)"
cal.UpdatePipeline()

# create a new 'Contour'
c1 = Contour(registrationName='Contour', Input=cal)
c1.ContourBy = ['POINTS', 'AVVMG']
c1.Isosurfaces = [1e-6]
c1.PointMergeMethod = 'Uniform Binning'
c1.UpdatePipeline()

########### SAVE FILES ###################
if("SOD" in codeName):
  savePath = casePath+"/ContourData_CRM.pvd"
  SaveData(savePath, proxy=c1)
elif("NEK" in codeName):
  savePath = casePath+"/ContourData_CRM.vtm"
  SaveData(savePath, proxy=c1,Writetimestepsasfileseries=0)
#SaveData('/Users/Gohan/Downloads/test.vtm', proxy=contour1, ChooseArraysToWrite=1,
#    PointDataArrays=['Normals', 'pressure'],
#    Writetimestepsasfileseries=1)
print("----|| NEK: SURFACE FILE WRITTEN ")
##############################

s1 = Slice(Input=case)
s1.SliceType = 'Plane'
s1.SliceOffsetValues = [0.0]
s1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
s1.SliceType.Normal = [0.0, 0.0, 1.0]
s1.UpdatePipeline()

########### SAVE FILES ###################
if("SOD" in codeName):
  savePath = casePath+"/SliceData_zmid.pvd"
  SaveData(savePath, proxy=s1, WriteTimeSteps=ifTSeries)
elif("NEK" in codeName):
  savePath = casePath+"/SliceData_zmid.vtm"
  SaveData(savePath, proxy=s1,Writetimestepsasfileseries=0)
#SaveData('/Users/Gohan/Downloads/test.vtm', proxy=contour1, ChooseArraysToWrite=1,
#    PointDataArrays=['Normals', 'pressure'],
#    Writetimestepsasfileseries=1)
print("----|| NEK: SLICED FILE WRITTEN ")
##############################

# create a new 'Contour'
c2 = Contour(registrationName='Contour', Input=case)
c2.ContourBy = ['POINTS', 'QCRIT']
c2.Isosurfaces = [qVal]
c2.PointMergeMethod = 'Uniform Binning'
c2.UpdatePipeline()

########### SAVE FILES ###################
if("SOD" in codeName):
  savePath = casePath+"/ContourData_QCrit.pvd"
  SaveData(savePath, proxy=c2, WriteTimeSteps=ifTSeries)
elif("NEK" in codeName):
  savePath = casePath+"/ContourData_QCrit.vtm"
  SaveData(savePath, proxy=c2,Writetimestepsasfileseries=0)
print("----|| NEK: QCRIT FILE WRITTEN ")
################################################ 

