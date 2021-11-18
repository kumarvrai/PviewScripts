#### import the simple module from the paraview
import os
import time
import operator
import numpy as np
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *

caseName	= sys.argv[1]
model		= sys.argv[2]
nu		= float(sys.argv[3])
dim		= sys.argv[4]
geom		= sys.argv[5]

casePath = os.getcwd()

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

print "--|| ALYA :: READING ALYA-AVERAGED ARRAYS"
startTime = time.time()

#case = EnSightReader(CaseFileName=os.path.join(caseName+'.ensi.case'))

fileName = caseName+'.ensi.case'

case = OpenDataFile(fileName)

fileName = caseName+'.pvd'

if(model == "LES"):
 case.PointArrays = ['TURBU','YPLUS','AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
elif(model == "DNS"):
 case.PointArrays = ['AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
else:
 case.PointArrays = ['YPLUS', 'AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']

case.UpdatePipeline()
print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'


print "--|| ALYA :: TEMPORAL AVERAGING  ALYA-AVERAGED ARRAYS"
startTime = time.time()
nacaensicase = TemporalStatistics(Input=case)

# Properties modified on temporalStatistics1
nacaensicase.ComputeMinimum = 0
nacaensicase.ComputeMaximum = 0
nacaensicase.ComputeStandardDeviation = 0
nacaensicase.UpdatePipeline()

print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'

caseVarNames = case.PointArrays
caseVarNames = [s + '_average' for s in caseVarNames]
indU = int([i for i, s in enumerate(caseVarNames) if 'AVVEL' in s][0]);
indP = int([i for i, s in enumerate(caseVarNames) if 'AVPRE' in s][0]);
indXX = int([i for i, s in enumerate(caseVarNames) if 'AVVE2' in s][0]);
indXY = int([i for i, s in enumerate(caseVarNames) if 'AVVXY' in s][0]);
print "--|| ALYA :: CALCULATING R-STRESSES"
startTime = time.time()
# CALCULATE RStresses
CAL1 = Calculator(Input=nacaensicase)
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
# GRADIENT CALC
print("--|| ALYA :: CALCULATING PRESS GRADIENT")
startTime = time.time()
CAL1 = GradientOfUnstructuredDataSet(Input=CAL1)
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

########### 3D STATISTICS ###################
if('3D' in dim):
 ## create a new 'Programmable Filter and change names'
 print "--|| ALYA: CHANGING VARIABLE NAMES USING A PROGRAMMABLE FILTER"
 startTime = time.time()
 CAL1 = ProgrammableFilter(Input=CAL1)
 CAL1.Script = \
 """
 import numpy as np
 varNames = inputs[0].PointData.keys()
 print "----|| Alya :: ALL 3D ARRAYS --> ",varNames
 for var in varNames:
  outName = str(var[0:5])
  avg = (inputs[0].PointData[var])
  output.PointData.append(avg,outName)
 """
 CAL1.UpdatePipeline()
 print "--|| ALYA: DONE. TIME =",time.time()-startTime,'sec'
 # Save a 3D time averaged file
 savePath = casePath+"/AvgData_3D.vtm"
 SaveData(savePath, proxy=CAL1)
 print "----|| ALYA: 3D STATISTICS FILE WRITTEN "
 #exit()
################################################ 

print "--|| ALYA :: EVALUATING DIMENSIONS FOR SPANWISE AVERAGE"
startTime = time.time()

(xmin,xmax,ymin,ymax,zmin,zmax) =  case.GetDataInformation().GetBounds()

if("DFUSER" in geom):
  if('.pvd' in fileName):
    geofile_in = '../'+caseName+'.geo.dat'
    xperFile   = '../'+caseName+'.xPerCoord.dat'
    zperFile   = '../'+caseName+'.zPerCoord.dat'
  else:
    geofile_in = caseName+'.geo.dat'
    xperFile   = caseName+'.xPerCoord.dat'
    zperFile   = caseName+'.zPerCoord.dat'

  print("--||ALYA: READING PLANAR COORDINATES")
  xCenter = np.round((xmin+xmax)/2,decimals=2)
  yCenter = np.round((ymin+ymax)/2,decimals=2)
  cf = open(zperFile,'r')
  targetCoord = np.loadtxt(cf)
  cf.close()
  xPlane = np.round(targetCoord[:,0],decimals=6); 
  yPlane = np.round(targetCoord[:,1],decimals=6);
  rPlane = np.sort(np.unique(np.round(np.sqrt((xPlane)**2+(yPlane)**2),decimals=4)))
  thPlane = np.sort(np.unique(np.round(np.arctan2(yPlane,xPlane),decimals=2)));
  print("----|| ALYA: PLANE AT (XC=%.3f,YC=%.3f) WITH %d R- AND %d TH- PLANES" % (xCenter,yCenter,len(rPlane),len(thPlane)))
  N = len(thPlane);
  thMax = np.amax(thPlane,axis=None); thMin = np.amin(thPlane,axis=None);
  thMid = (thMax+thMin)/2
  zpos = np.arange(N)*(thMax-thMin)/(N-1)
  print("----|| ALYA: GRID SIZE IN THETA = %.4f" % ((thMax-thMin)/(N-1)))
  
  resample_transforms=list();
  data=list();
  
  slice1 = Slice(Input=CAL1)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  ## init the 'Plane' selected for 'SliceType'
  slice1.SliceType.Origin = [0.0, 0.0, 0.0]
  ## Properties modified on slice1.SliceType
  slice1.SliceType.Normal = [-np.sin(thMid), np.cos(thMid), 0.0]

  print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
  
else:

  if('.pvd' in fileName):
    geofile_in = '../'+caseName+'.geo.dat'
    eperFile   = '../'+caseName+'.eper.dat'
    zperFile   = '../'+caseName+'.zper.dat'
  else:
    geofile_in = caseName+'.geo.dat'
    eperFile   = caseName+'.eper.dat'
    zperFile   = caseName+'.zper.dat'
  
  bString = 'COORDINATES'
  eString = 'END_COORDINATES'
  nNodes = int(os.popen('sed -n ''/%s/,/%s/p'' %s | wc -l ' % (bString,eString,geofile_in)).read())
  nNodes = nNodes-2
  print '----||ALYA :: WORKING WITH %i NODES IN TOTAL' % nNodes
  nNodesPlane = 0
  if(caseName == "chan"):
    nNodesPlane = int(os.popen('wc -l < %s'%eperFile).read())
  nNodesPlane = nNodesPlane + len(open('%s'%zperFile).readlines()) - 2
  print '----||ALYA :: WORKING WITH %i NODES ON A PLANE' % nNodesPlane
  
  print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'

  ########### SPANWISE AVERAGING ################
  print "--|| ALYA :: GENERATING SLICE FOR SPANWISE AVERAGE"
  startTime = time.time()
  
  ## generate z-axis slices
  if(caseName == "chan"):
    N = int(np.ceil(float(nNodes)/nNodesPlane))
  else:
    N = int(np.floor(float(nNodes)/nNodesPlane))
  #N = input("\n--|| ALYA: ENTER NUMBER OF Z-PLANES IN THIS SIMULATION..")
  print("--|| ALYA: WORKING WITH %d Z-PLANES" % (N))
  print("----|| ALYA: DELTA-Z = %f" % ((zmax-zmin)/(N-1)))
  
  resample_transforms=list();
  data=list();
  
  #slice1 = Slice(Input=nacaensicase)
  slice1 = Slice(Input=CAL1)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
  slice1.SliceType.Normal = [0.0, 0.0, 1.0]
  
  zmid = (zmin+zmax)/2
  zpos = np.arange(N)*(zmax-zmin)/(N-1)
  print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
  
print "--|| ALYA :: CREATING TRANSFORMATIONS"
startTime = time.time()

for i in range(N):
	# create a new 'Transform'
	transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
	resampleWithDataset1 = ResampleWithDataset(Input=CAL1,Source=transform1,guiName="resample{}".format(i))
	resample_transforms.append(resampleWithDataset1)
	# Properties modified on transform1.Transform
	if("DFUSER" in geom):
	  transform1.Transform.Rotate = [0.0, 0.0, zpos[i]-thMid]
	else:
	  transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
	#data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])
print "--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec'
HideAll()


## create a new 'Programmable Filter'
print "--|| ALYA: AVERAGING USING A PROGRAMMABLE FILTER"
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
print "----|| Alya - WORKING ON ARRAYS::",varFull
N=len(inputs)-1;
print "--|| ALYA: AVERAGING %d DATA-PLANES" % (N)

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
print "--|| ALYA: SPANWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec'

if('2D' in dim):
  if("DFUSER" in geom):
    # Convert the plane (x,y,z) to (z,r,th) plane
    PF1 = Calculator(Input=PF1)
    PF1.ResultArrayName = "result"
    PF1.CoordinateResults = 1
    PF1.Function = "coordsZ*iHat + sqrt(coordsX^2+coordsY^2)*jHat"
    PF1.UpdatePipeline()
  savePath = casePath+"/AvgData_2D.vtm"
  SaveData(savePath, proxy=PF1)
  savePath = casePath+"/AvgData_2D.csv"
  SaveData(savePath, proxy=PF1)
  print "----|| ALYA: 2D STATISTICS FILE WRITTEN AS: ",savePath

  ########### STREAMWISE AVERAGING ################
if('1D' in dim):
  ## read Coordinates File
  xChanCoord = np.loadtxt('xChanCoord.dat')
  yChanCoord = np.loadtxt('yChanCoord.dat')
  
  N = len(xChanCoord)
  #N = input("\n--|| ALYA: ENTER NUMBER OF Z-PLANES IN THIS SIMULATION..")
  print("--|| ALYA: WORKING WITH %d X-PLANES" % (N))
  
  resample_transforms=list();
  data=list();
  
  slice1 = Slice(Input=PF1)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  ## init the 'Plane' selected for 'SliceType'
  slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
  ## Properties modified on slice1.SliceType
  slice1.SliceType.Normal = [1.0, 0.0, 0.0]
  
  xmid = (xmin+xmax)/2
  xpos = np.arange(N)*(xmax-xmin)/(N-1)
  
  print "--|| ALYA: CREATING X TRANSFORMATIONS"
  startTime = time.time()
  for i in range(N):
  	# create a new 'Transform'
  	transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
  	#resampleWithDataset1=ResampleWithDataset(Input=GOUD1,Source=transform1,guiName="resample{}".format(i))
  	resampleWithDataset1=ResampleWithDataset(Input=PF1,Source=transform1,guiName="resample{}".format(i))
  	resample_transforms.append(resampleWithDataset1)
  	# Properties modified on transform1.Transform
  	transform1.Transform.Translate = [xpos[i]-xmid, 0.0, 0.0]
  	#data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])
  print "--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec'
  
  HideAll()
  
  ## create a new 'Programmable Filter'
  print "--|| ALYA: X AVERAGING USING A PROGRAMMABLE FILTER"
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
  ### first input is the grid
  ### the rest of them are data to be averaged
  PF1.Script = \
  """
  import numpy as np
  
  varFull = []
  varFull = inputs[0].PointData.keys()
  print "----|| Alya - WORKING ON ORIGINAL ARRAYS::",varFull
  N=len(inputs);
  
  for varName in varFull:
     varName0=str(varName[0:5])
     #varName=varName+"_average"
     avg = (inputs[0].PointData[varName])
     for i in range(1,N):
         d = inputs[i].PointData[varName]
         avg = avg + d
     avg = avg/N
     output.PointData.append(avg,varName0)
  """
  PF1.UpdatePipeline()
  print "--|| ALYA: STREAMWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec'
  
  #Show(programmableFilter1)
  ##
  #### write
  print "--|| ALYA: SAVING THE AVERAGED FILES"
  startTime = time.time()
  savePath = casePath+"/AvgData_1D.csv"
  SaveData(savePath, proxy=PF1)
  print "----|| ALYA: 1D STATISTICS FILE WRITTEN AS: ",savePath
  print "--|| ALYA: FILE SAVED. TIME =",time.time()-startTime,'sec'
