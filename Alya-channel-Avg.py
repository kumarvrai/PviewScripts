#### import the simple module from the paraview
import os
import time
import operator
import numpy as np
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *

#pxm  = servermanager.ProxyManager()
#pxm.GetVersion()
#print("--|| NEK :: USING PARAVIEW VERSION",pxm)


caseName	= sys.argv[1]
codeName	= sys.argv[2]
fileType	= sys.argv[3]
avgDim		= sys.argv[4]
geomType	= sys.argv[5]
rotDeg		= float(sys.argv[6])

zDec = 6; xDec = 6
xc = 1.0; yc = 0.0; rot=np.radians(rotDeg);
casePath = os.getcwd()

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

if('NEK' in codeName):
  print("--|| NEK :: READING NEK5000 ARRAYS")
  startTime = time.time()
  
  fileName1 = 'avg'+caseName+'.nek5000'
  fileName2 = 'rms'+caseName+'.nek5000'
  fileName3 = 'rm2'+caseName+'.nek5000'
  
  case1 = OpenDataFile(fileName1)
  case1.PointArrays = ['velocity','pressure']
  case1.UpdatePipeline()
  
  ## create a new 'Programmable Filter and change names'
  print("--|| NEK: CHANGING VARNAMES 1 USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  case1 = ProgrammableFilter(Input=case1)
  case1.Script = \
  """
  import numpy as np
  varNames0 = ['velocity','pressure']
  varNames1 = ['AVVEL','AVPRE']
  for (i,var) in enumerate(varNames0):
   outName = varNames1[i]
   avg = (inputs[0].PointData[var])
   output.PointData.append(avg,outName)
  """
  case1.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  case2 = OpenDataFile(fileName2)
  case2.PointArrays = ['velocity','pressure']
  case2.UpdatePipeline()
  
  ## create a new 'Programmable Filter and change names'
  print("--|| NEK: CHANGING VARNAMES 2 USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  case2 = ProgrammableFilter(Input=case2)
  case2.Script = \
  """
  import numpy as np
  varNames0 = ['velocity','pressure']
  varNames1 = ['AVVE2','AVPR2']
  for (i,var) in enumerate(varNames0):
   outName = varNames1[i]
   avg = (inputs[0].PointData[var])
   output.PointData.append(avg,outName)
  """
  case2.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  case3 = OpenDataFile(fileName3)
  case3.PointArrays = ['velocity']
  case3.UpdatePipeline()
  
  ## create a new 'Programmable Filter and change names'
  print("--|| NEK: CHANGING VARNAMES 3 USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  case3 = ProgrammableFilter(Input=case3)
  case3.Script = \
  """
  import numpy as np
  varNames0 = ['velocity']
  varNames1 = ['AVVXY']
  for (i,var) in enumerate(varNames0):
   outName = varNames1[i]
   avg = (inputs[0].PointData[var])
   output.PointData.append(avg,outName)
  """
  case3.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  print("--|| NEK: APPEND DATASETS")
  startTime = time.time()
  case = AppendAttributes(Input=[case1,case2,case3])
  case.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  if("ROT" in geomType):
    print("--|| NEK: ROTATE VARIABLES")
    case = Calculator(Input=case)
    case.ResultArrayName = "AVVEL"
    case.Function = "(AVVEL_X*cos(%f) - AVVEL_Y*sin(%f))*iHat \
                     + (AVVEL_X*sin(%f) + AVVEL_Y*cos(%f))*jHat + AVVEL_Z*kHat" \
                    % (rot,rot,rot,rot)
    case.UpdatePipeline()
    case = Calculator(Input=case)
    case.ResultArrayName = "AVVE2_ROT"
    case.Function = "(AVVE2_X*cos(%f)^2 + AVVE2_Y*sin(%f)^2 - 2*AVVXY_X*sin(%f)*cos(%f))*iHat \
                     +(AVVE2_X*sin(%f)^2 + AVVE2_Y*cos(%f)^2 + 2*AVVXY_X*sin(%f)*cos(%f))*jHat \
                     + AVVE2_Z*kHat" \
                    % (rot,rot,rot,rot,rot,rot,rot,rot)
    case.UpdatePipeline()
    case = Calculator(Input=case)
    case.ResultArrayName = "AVVXY_ROT"
    case.Function = "((AVVE2_X-AVVE2_Y)*sin(%f)*cos(%f)+AVVXY_X*(cos(%f)^2-sin(%f)^2))*iHat \
                     +(AVVXY_Z*sin(%f) + AVVXY_Y*cos(%f))*jHat \
                     +(AVVXY_Z*cos(%f) - AVVXY_Y*sin(%f))*kHat" \
                    % (rot,rot,rot,rot,rot,rot,rot,rot)
    case.UpdatePipeline()
    case = Calculator(Input=case)
    case.CoordinateResults = 1
    case.ResultArrayName = "result"
    case.Function = "((coordsX-%f)*cos(%f)-(coordsY-%f)*sin(%f) + %f)*iHat \
                     +((coordsX-%f)*sin(%f)+(coordsY-%f)*cos(%f) + %f)*jHat \
                     +coordsZ*kHat" \
                    % (xc,rot,yc,rot,xc,xc,rot,yc,rot,yc)
    case.UpdatePipeline()
    startTime = time.time()
    case = ProgrammableFilter(Input=case)
    case.Script = \
    """
    import numpy as np
    varNames0 = ['AVVEL','AVPRE','AVVE2_ROT','AVVXY_ROT','AVPR2']
    varNames1 = ['AVVEL','AVPRE','AVVE2','AVVXY','AVPR2']
    for (i,var) in enumerate(varNames0):
     outName = varNames1[i]
     avg = (inputs[0].PointData[var])
     output.PointData.append(avg,outName)
    """
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
  case.UpdatePipeline()
  print("--|| SOD :: DONE. TIME =",time.time()-startTime,'sec')
elif("PVD" in codeName):
  print("--|| INFO :: READING PVD ARRAYS")
  startTime = time.time()
  fileName = caseName+'.pvd'
  case = OpenDataFile(fileName)
  case.UpdatePipeline()
  print("--|| PVD :: DONE. TIME =",time.time()-startTime,'sec')
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

print("--|| NEK :: TEMPORAL AVERAGING.")
startTime = time.time()
case = TemporalStatistics(Input=case)

# Properties modified on temporalStatistics1
case.ComputeMinimum = 0
case.ComputeMaximum = 0
case.ComputeStandardDeviation = 0
case.UpdatePipeline()
print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
## create a new 'Programmable Filter and change names'
print("--|| NEK: CHANGING VARIABLE NAMES.")
startTime = time.time()
case = ProgrammableFilter(Input=case)
case.Script = \
"""
import numpy as np
varNames = inputs[0].PointData.keys()
print("----|| Alya :: ALL 3D ARRAYS --> ",varNames)
for var in varNames:
 outName = str(var[0:5])
 avg = (inputs[0].PointData[var])
 output.PointData.append(avg,outName)
"""
case.UpdatePipeline()
print("--|| NEK: DONE. TIME =",time.time()-startTime,'sec')

Ntotal = int(case.GetDataInformation().GetNumberOfPoints())
print("----|| ALYA :: WORKING WITH ",Ntotal," TOTAL NUMBER OF POINTS")

caseVarNames = case.PointData.keys()
indU = int([i for i, s in enumerate(caseVarNames) if 'AVVEL' in s][0]);
indP = int([i for i, s in enumerate(caseVarNames) if 'AVPRE' in s][0]);
if(codeName in str(["NEK","ALYA","SOD"])):
  indXX = int([i for i, s in enumerate(caseVarNames) if 'AVVE2' in s][0]);
  indXY = int([i for i, s in enumerate(caseVarNames) if 'AVVXY' in s][0]);
  print("--|| NEK :: CALCULATING R-STRESSES")
  startTime = time.time()
  # CALCULATE RStresses
  case = Calculator(Input=case)
  case.ResultArrayName = "RS_II"
  case.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXX],caseVarNames[indU],'X',caseVarNames[indU],'X',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Y',\
                   caseVarNames[indU],'Z',caseVarNames[indU],'Z')
  case.UpdatePipeline() 
  case = Calculator(Input=case)
  case.ResultArrayName = "RS_IJ"
  case.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                  % (caseVarNames[indXY],caseVarNames[indU],'X',caseVarNames[indU],'Y',\
                    caseVarNames[indU],'Y',caseVarNames[indU],'Z',\
                    caseVarNames[indU],'X',caseVarNames[indU],'Z')
  case.UpdatePipeline()
  print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
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
case.ComputeQCriterion = 1
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
  #QuerySelect(QueryString='(mag(avvel) == 0)', 
  #            FieldType='POINT', InsideOut=0)
  #generateSurfaceNormals1 = FindSource('GenerateSurfaceNormals1')            
  #SetActiveSource(case_clcd)            
  #case_clcd = ExtractSelection(registrationName='ExtractSelection1',Input=case_clcd)
  #case_clcd.UpdatePipeline()
  case_clcd = Clip(Input=case_clcd)
  case_clcd.ClipType = 'Box'
  case_clcd.ClipType.Position = [-1, -1, 0]
  case_clcd.ClipType.Length = [3, 3, 1]
  case_clcd.Invert = 1
  case_clcd.UpdatePipeline()
  case_clcd = Calculator(Input=case_clcd)
  case_clcd.ResultArrayName = "AVGCF"
  case_clcd.Function = "(1/500)*sqrt(dot((""AVVGR_0""*iHat+""AVVGR_1""*jHat),-Normals)^2+dot((""AVVGR_3""*iHat+""AVVGR_4""*jHat),-Normals)^2)"
  case_clcd.UpdatePipeline()
  case_clcd = Calculator(Input=case_clcd)
  case_clcd.ResultArrayName = "AVGCP"
  case_clcd.Function = "AVPRE"
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

########### 3D STATISTICS ###################
if('3D' in avgDim):
 # Save a 3D time averaged file
 savePath = casePath+"/AvgData_3D.pvd"
 savePath = casePath+"/AvgData_2D.vtm"
 #savePath = casePath+"/AvgData_3D.csv"
 SaveData(savePath, proxy=case)
 print("----|| NEK: 3D STATISTICS FILE WRITTEN ")
 #slice1 = Slice(Input=case)
 #slice1.SliceType = 'Plane'
 #slice1.SliceOffsetValues = [0.0]
 #slice1.SliceType.Origin = [3.0, 1.0, 1.5]
 #slice1.SliceType.Normal = [0.0, 0.0, 1.0]
 #slice1.UpdatePipeline()
 #slice1 = Slice(Input=slice1)
 #slice1.SliceType = 'Plane'
 #slice1.SliceOffsetValues = [0.0]
 #slice1.SliceType.Origin = [3.0, 1.0, 1.5]
 #slice1.SliceType.Normal = [1.0, 0.0, 0.0]
 #slice1.UpdatePipeline()
 #savePath = casePath+"/PySodAvgData_1D.csv"
 #SaveData(savePath, proxy=slice1)
 #print("----|| NEK: TAVG 1D FILE WRITTEN ")
################################################ 

print("--|| NEK :: EVALUATING DIMENSIONS FOR SPANWISE AVERAGE")
startTime = time.time()

(xmin,xmax,ymin,ymax,zmin,zmax) =  case.GetDataInformation().GetBounds()
print("----|| INFO: BOX SIZE = %.2f %.2f %.2f %.2f %.2f %.2f"%(xmin,xmax,ymin,ymax,zmin,zmax))

if("DFUSER" in geomType):
  
  slice1 = Slice(Input=case)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  ## init the 'Plane' selected for 'SliceType'
  slice1.SliceType.Origin = [0.0, 0.0, 0.0]
  ## Properties modified on slice1.SliceType
  slice1.SliceType.Normal = [1.0, 0.0, 0.0]
  slice1.UpdatePipeline()
  #------------------------------#
  slice1 = Clip(Input=slice1)
  slice1.ClipType = 'Plane'
  slice1.ClipType.Origin = [0.0, 0.0, 0.0]
  slice1.ClipType.Normal = [0.0, 1.0, 0.0]
  slice1.UpdatePipeline()

  Nplane = int(slice1.GetDataInformation().GetNumberOfPoints())
  print("----|| ALYA :: WORKING WITH ",Nplane," PLANAR POINTS")
  
  N = int(Ntotal/Nplane)
  thMax = 2.0*np.pi; thMin = 0.0;
  thMid = np.pi/2
  zpos = np.arange(N)*(thMax-thMin)/(N-1)
  print("----|| ALYA: WORKING WITH %d THETA-PLANES" % (len(zpos)))
  print("----|| ALYA: DELTA-THETA = %f" % ((thMax-thMin)/(N-1)))
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

else:
  slice1 = Slice(Input=case)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
  slice1.SliceType.Normal = [0.0, 0.0, 1.0]
  slice1.UpdatePipeline()
  
  Nplane = int(slice1.GetDataInformation().GetNumberOfPoints())
  print("----|| INFO :: WORKING WITH ",Nplane," PLANAR POINTS")
  
  N = int(Ntotal/Nplane)
  zmid = (zmin+zmax)/2
  zpos = np.around(np.asarray(np.arange(N)*(zmax-zmin)/(N-1),dtype=np.double),decimals=zDec)
  delta_z = (zmax-zmin)/(N-1)
  print("----|| ALYA: WORKING WITH %d Z-PLANES" % (len(zpos)))
  print("----|| ALYA: DELTA-Z = %f" % (delta_z))
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

########### PERFORM AVERAGING ################
print("--|| NEK :: CREATING TRANSFORMATIONS")
startTime = time.time()
resample_transforms=list();
data=list();

for i in range(N):
	# create a new 'Transform'
	transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
	# Properties modified on transform1.Transform
	if("DFUSER" in geomType):
	  transform1.Transform.Rotate = [0.0, 0.0, zpos[i]-thMid]
	else:
	  transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
	try:  
	 resampleWithDataset1 = ResampleWithDataset(Input=case,Source=transform1)
	except: 
	 resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=case,DestinationMesh=transform1)
	resample_transforms.append(resampleWithDataset1)
print("--|| NEK: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
HideAll()


## create a new 'Programmable Filter'
print("--|| NEK: AVERAGING USING A PROGRAMMABLE FILTER")
startTime = time.time()

PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
### first input is the grid
### the rest of them are data to be averaged
PF1.Script = \
"""
#from vtk.numpy_interface import dataset_adapter as dsa
#from vtk.numpy_interface import algorithms as algs
import numpy as np

varFull = []
varFull = inputs[0].PointData.keys()
print("----|| Alya - WORKING ON ARRAYS::",varFull)
N=len(inputs)-1;
print("--|| NEK: AVERAGING %d DATA-PLANES" % (N))

for varName in varFull:
   varName0=varName[0:5]
   avg = 0.0*(inputs[0].PointData[varName])
   for i in range(N):
       d = inputs[i+1].PointData[varName]
       avg = avg + d
   avg = avg/N
   output.PointData.append(avg,varName0)
"""
PF1.UpdatePipeline()
print("--|| NEK: SPANWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec')

if('2D' in avgDim):
  #if("DFUSER" in geomType):
  #  # Convert the plane (x,y,z) to (z,r,th) plane
  #  PF1 = Calculator(Input=PF1)
  #  PF1.ResultArrayName = "result"
  #  PF1.CoordinateResults = 1
  #  PF1.Function = "coordsZ*iHat + sqrt(coordsX^2+coordsY^2)*jHat"
  #  PF1.UpdatePipeline()
  if(codeName in str(["NEK","ALYA","SOD"])):
    savePath = casePath+"/AvgData_2D.vtm"
  #elif(codeName in str(["SOD"])):
  #  savePath = casePath+"/AvgData_2D.pvd"
  else:
    savePath = casePath+"/AvgData_2D.vtp"
  savePath = casePath+"/AvgData_2D.pvd"
  SaveData(savePath, proxy=PF1)
  savePath = casePath+"/AvgData_2D.csv"
  SaveData(savePath, proxy=PF1)
  print("----|| NEK: 2D STATISTICS FILE WRITTEN AS: ",savePath)

  ########### STREAMWISE AVERAGING ################
if('1D' in avgDim):
  print("--|| ALYA :: GENERATING SLICE FOR STREAMWISE AVERAGE")
  slice1 = Slice(Input=PF1)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  ## init the 'Plane' selected for 'SliceType'
  slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
  ## Properties modified on slice1.SliceType
  if("DFUSER" in geomType):
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]
  else:
    slice1.SliceType.Normal = [1.0, 0.0, 0.0]
  slice1.UpdatePipeline()
  
  Ny = int(slice1.GetDataInformation().GetNumberOfPoints())
  print("----|| ALYA :: WORKING WITH ",Ny," PLANAR POINTS")

  N = int(Nplane/Ny)
  if("DFUSER" in geomType):
    xmid = (zmin+zmax)/2
    xpos = np.around(np.asarray(np.arange(N)*(zmax-zmin)/(N-1),dtype=np.double),decimals=zDec)
  else:  
    xmid = (xmin+xmax)/2
    xpos = np.around(np.asarray(np.arange(N)*(xmax-xmin)/(N-1),dtype=np.double),decimals=xDec)
  print("----|| ALYA: WORKING WITH %d X-PLANES" % (len(xpos)))
  print("----|| ALYA: DELTA-X = %f" % ((np.amax(xpos)-np.amin(xpos))/(N-1)))
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| NEK: CREATING X TRANSFORMATIONS")
  resample_transforms=list();
  data=list();
  
  startTime = time.time()
  for i in range(N):
    # create a new 'Transform'
    transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
    ## Properties modified on transform1.Transform
    if("DFUSER" in geomType):
      transform1.Transform.Translate = [0.0, 0.0, xpos[i]-xmid]  
    else:  
      transform1.Transform.Translate = [xpos[i]-xmid, 0.0, 0.0]  
    try:
      resampleWithDataset1 = ResampleWithDataset(Input=PF1,Source=transform1)
    except:
      resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=PF1,DestinationMesh=transform1)
    resample_transforms.append(resampleWithDataset1)
    transform1.UpdatePipeline()
  print("--|| NEK: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
  
  ## create a new 'Programmable Filter'
  print("--|| NEK: X AVERAGING USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
  PF1.Script = \
  """
  import numpy as np
  varFull = inputs[0].PointData.keys()
  print("----|| INFO : WORKING ON ",varFull)
  N=len(inputs)-1;
  for varName in varFull:
     avg = 0.0*(inputs[0].PointData[varName])
     for i in range(N):
         d = inputs[i+1].PointData[varName]
         avg = avg + d
     avg = avg/N
     output.PointData.append(avg,varName)
  """
  PF1.UpdatePipeline()
  print("--|| NEK: STREAMWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec')
  #
  #### write
  print("--|| NEK: SAVING THE AVERAGED FILES")
  startTime = time.time()
  #savePath = casePath+"/AvgData_1D.vtm"
  #SaveData(savePath, proxy=PF1)
  savePath = casePath+"/AvgData_1D.csv"
  SaveData(savePath, proxy=PF1)
  print("----|| NEK: 1D STATISTICS FILE WRITTEN AS: ",savePath)
  print("--|| NEK: FILE SAVED. TIME =",time.time()-startTime,'sec')
