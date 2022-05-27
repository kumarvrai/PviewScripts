#### import the simple module from the paraview
import os
import time
import sys
import operator
import numpy as np
from paraview.simple import *
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs 
from scipy.interpolate import griddata

casePath = os.getcwd()

caseName = sys.argv[1]
model    = sys.argv[2]
nu       = float(sys.argv[3])
dim      = sys.argv[4]
mode     = sys.argv[5]
method     = sys.argv[6]
file_fmt     = sys.argv[7]
zDec = 6; xDec = 6

c = vtk.vtkMultiProcessController.GetGlobalController()
mpi_ranks = c.GetNumberOfProcesses()
print("--|| ALYA :: WORKING WITH %d RANKS." % mpi_ranks)
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
#
print("--|| ALYA :: READING ALYA-AVERAGED ARRAYS")
startTime = time.time()
if('ENSI' in file_fmt):
 fileName = caseName+'.ensi.case'
elif('VTK' in file_fmt): 
 fileName = caseName+'.pvd'
elif('VTM' in file_fmt): 
 fileName = 'AvgData_3D.vtm'
elif('BUDPVD' in file_fmt): 
 fileName = './'+caseName+'.pvd'
else:
 raise ValueError('--|| ALYA ERROR :: FILE_FMT NOT RECONIZED.')
case = OpenDataFile(fileName)

if('BUD' not in file_fmt):
  if("FAVG" in mode):
   if(model == "LES"):
    case.PointArrays = ['TURBU','YPLUS','AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
   elif(model == "DNS"):
    case.PointArrays = ['AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
   else:
    case.PointArrays = ['YPLUS', 'AVVEL', 'AVPRE', 'AVTAN', 'AVVE2', 'AVVXY']
  elif("SAVG" in mode):
   if(model == "LES"):
    case.PointArrays = ['TURBU','YPLUS','VELOC','PRESS']
   elif(model == "DNS"):
    case.PointArrays = ['VELOC','PRESS']
   else:
    case.PointArrays = ['TURBU','VELOC','PRESS']
else:
   if(model == "BASIC"):
    case.PointArrays = ['AVVEL', 'AVPRE', 'AVPGR', 'AVVGR', 'AVVE2', 'AVVXY', 'AVTAN', 'RS_II', 'RS_IJ']
   elif(model == "TSTEP"):
    case.PointArrays = ['VELOC', 'PRESS']
case.UpdatePipeline()
caseVarNames = case.PointData.keys()
print("--|| ALYA : LOADED VARIABLES", caseVarNames)
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

inpFile = os.getcwd()+'/'+'inputFile.py'
f = open(inpFile,'w')
f.write("caseName = '%s'\n" % caseName)
f.write("fileName = '%s'\n" % fileName)
f.write("model = '%s'\n" % model)
f.write("nu = %f\n" % nu)
f.write("dim = '%s'\n" % dim)
f.write("mode = '%s'\n" % mode)
f.write("method = '%s'\n" % method)
f.write("file_fmt = '%s'\n" % file_fmt)
f.write("xDec = %d\n" % xDec)
f.write("zDec = %d\n" % zDec)
f.close()

if('D3' in file_fmt):
    print("--|| ALYA :: APPLYING D3 FILTER")
    startTime = time.time()
    case = D3(Input=case)
    case.UpdatePipeline()
    print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
if('CLEAN' in file_fmt):
    print("--|| ALYA :: APPLYING CLEAN TO GRID FILTER")
    startTime = time.time()
    case = CleantoGrid(Input=case)
    case.UpdatePipeline()
    print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')


Ntotal = int(case.GetDataInformation().GetNumberOfPoints())
print("----|| ALYA :: WORKING WITH ",Ntotal," TOTAL NUMBER OF POINTS")

if("SAVG" in mode):
   ## CALCULATE AVVE2 and AVVXY VARIABLES ###
   print("--|| ALYA :: CALCULATE AVVE2 VARIABLES")
   startTime = time.time()
   case = Calculator(Input=case)
   case.ResultArrayName = "VE2EL"
   case.Function = "VELOC_X*VELOC_X*iHat + VELOC_Y*VELOC_Y*jHat + VELOC_Z*VELOC_Z*kHat" 
   case.UpdatePipeline()
   case = Calculator(Input=case)
   case.ResultArrayName = "VXYEL"
   case.Function = "VELOC_X*VELOC_Y*iHat + VELOC_Y*VELOC_Z*jHat + VELOC_X*VELOC_Z*kHat" 
   case.UpdatePipeline()
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

elif("FAVG" in mode):
 if('SKIP' not in file_fmt):
   print("--|| ALYA :: TEMPORAL AVERAGING  ALYA-AVERAGED ARRAYS")
   startTime = time.time()
   case = TemporalStatistics(Input=case)
   # Properties modified on temporalStatistics1
   case.ComputeMinimum = 0
   case.ComputeMaximum = 0
   case.ComputeStandardDeviation = 0
   case.UpdatePipeline()
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
   caseVarNames = [s + '_average' for s in caseVarNames]
   indU = int([i for i, s in enumerate(caseVarNames) if 'AVVEL' in s][0]);
   indP = int([i for i, s in enumerate(caseVarNames) if 'AVPRE' in s][0]);
   if any("AVVE2" in s for s in caseVarNames):
    indXX = int([i for i, s in enumerate(caseVarNames) if 'AVVE2' in s][0]);
    indXY = int([i for i, s in enumerate(caseVarNames) if 'AVVXY' in s][0]);
    print("--|| ALYA :: CALCULATING R-STRESSES")
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
    print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
   # GRADIENT CALC
   print("--|| ALYA :: CALCULATING PRESS GRADIENT")
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
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
   print("--|| ALYA :: CALCULATING AVVEL GRADIENT, Q AND VORTICITY")
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
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
   ########### 3D STATISTICS ###################
   if('3D' in dim):
    ## create a new 'Programmable Filter and change names'
    print("--|| ALYA: CHANGING VARIABLE NAMES USING A PROGRAMMABLE FILTER")
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
    print("--|| ALYA: DONE. TIME =",time.time()-startTime,'sec')
    # Save a 3D time averaged file
    if('PVD' in file_fmt): 
     savePath = casePath+"/AvgData_3D.vtk"
    else:
     savePath = casePath+"/AvgData_3D.vtm"
    SaveData(savePath, proxy=case)
    print("----|| ALYA: 3D STATISTICS FILE WRITTEN ")
   ################################################

########### SPANWISE AVERAGING ################
print("--|| ALYA :: GENERATING SLICE FOR SPANWISE AVERAGE")
startTime = time.time()

(xmin,xmax,ymin,ymax,zmin,zmax) =  case.GetDataInformation().GetBounds()
print("--|| ALYA GEOM:: XMIN=",xmin,"XMAX=",xmax,"YMIN=",ymin,"YMAX=",ymax,"ZMIN=",zmin,"ZMAX=",zmax)


slice1 = Slice(Input=case)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.UpdatePipeline()

Nplane = int(slice1.GetDataInformation().GetNumberOfPoints())
print("----|| ALYA :: WORKING WITH ",Nplane," PLANAR POINTS")

N = int(Ntotal/Nplane)
#N = 145;
zmid = (zmin+zmax)/2
zpos = np.around(np.asarray(np.arange(N)*(zmax-zmin)/(N-1),dtype=np.double),decimals=zDec)
delta_z = (zmax-zmin)/(N-1)
print("----|| ALYA: WORKING WITH %d Z-PLANES" % (len(zpos)))
print("----|| ALYA: DELTA-Z = %f" % (delta_z))
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

if("PINTERP" in method):
  print("--|| ALYA :: CREATING TRANSFORMATIONS")
  startTime = time.time()
  resample_transforms=list();
  data=list();
  for i in range(N):
    # create a new 'Transform'
    transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
    # Properties modified on transform1.Transform
    transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
    try:
      resampleWithDataset1 = ResampleWithDataset(Input=case,Source=transform1)
    except:  
      resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=case,DestinationMesh=transform1)
    resample_transforms.append(resampleWithDataset1)
    #data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])
  print("--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA: AVERAGING USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
  PF1.Script = \
"""
import os
import vtk
import sys
import numpy as np

caseName = 'naca'
fileName = './naca.pvd'
model = 'DNS'
nu = 0.000005
dim = '2D'
mode = 'SAVG'
method = 'INTERP'
file_fmt = 'BUDPVD'
xDec = 6
zDec = 6

inpFile = os.getcwd()+'/'+'inputFile.py'
exec(open(inpFile).read())
exec(open(inpFile).read(), {'__file__': inpFile})

c = vtk.vtkMultiProcessController.GetGlobalController()
rank = c.GetLocalProcessId()
t = inputs[0].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())

varFull = inputs[0].PointData.keys()
#varFull = [x for x in varFull if len(x) == 5]

if(rank==0):
  print("----|| ALYA : ARRAYS AT TIME %.3f " % (t) ,varFull)

N=len(inputs)-1;

for varName in varFull:
 if("FAVG" in mode):
  varName0 = varName[0:5]
 else:
  varName0 = "AV"+varName[0:3]
 avg = 0.0*(inputs[0].PointData[varName])
 for i in range(N):
     d = inputs[i+1].PointData[varName]
     avg = avg + d
 avg = avg/N
 output.PointData.append(avg,varName0)
"""
  PF1.UpdatePipeline()
  print("--|| ALYA: SPANWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec')

else:
  caseVarOrig = slice1.PointData.keys()
  print("--|| ALYA :: CALCULATE SPANWISE AVERAGING")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=[slice1,case])
  ### first input is the grid
  ### the rest of them are data to be averaged
  PF1.Script = \
"""
import os
import vtk
import numpy as np
from scipy.interpolate import griddata

caseName = 'naca'
fileName = './naca.pvd'
model = 'DNS'
nu = 0.000005
dim = '2D'
mode = 'SAVG'
method = 'INTERP'
file_fmt = 'BUDPVD'
xDec = 6
zDec = 6

inpFile = os.getcwd()+'/'+'inputFile.py'
try:
  exec(open(inpFile).read())
except:
  exec(open(inpFile).read(), {'__file__': inpFile})

c = vtk.vtkMultiProcessController.GetGlobalController()
rank = c.GetLocalProcessId()

t = inputs[1].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
varFull = inputs[1].PointData.keys()
#varFull = [x for x in varFull if len(x) == 5]

if(rank==0):
  print("--|| ALYA :: CALCULATING FOR",varFull," AT T=",t) 

try:
 d = dsa.WrapDataObject(inputs[1].GetBlock(0))
except: 
 d = dsa.WrapDataObject(inputs[1].VTKObject)
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=zDec)

## SLICE GEOMETRY
try:
 out = dsa.WrapDataObject(inputs[0].GetBlock(0))
except:
 out = dsa.WrapDataObject(inputs[0].VTKObject)

x_2d = np.around(np.asarray(out.Points[:,0],dtype=np.double),decimals=6)
y_2d = np.around(np.asarray(out.Points[:,1],dtype=np.double),decimals=6)
ind_2d = np.lexsort((y_2d,x_2d));

if(rank==0):
  print("----|| ALYA CHECK :: SHAPE OF SLICE", np.shape(ind_2d))
z_vec = np.unique(z)
N = len(z_vec)

for varName in varFull:
 #print("--|| ALYA :: CALCULATING FOR",varName)
 if("FAVG" in mode):
  varName0 = varName[0:5]
 else:
  varName0 = "AV"+varName[0:3]
 avg = 0.0*np.asarray(out.PointData[varName],dtype=np.double)
 #s_src = 0.0*np.asarray(out.PointData[varName],dtype=np.double)
 src = np.asarray(d.PointData[varName],dtype=np.double)
 for n in range(N):
  ind_z = np.where(z == z_vec[n])
  #print("--|| ALYA CHECK :: AT n,Z=",n,zpos[n],"SHAPE OF MASTER DATASET", np.shape(ind_z))
  x_l = x[ind_z];
  y_l = y[ind_z];
  s_src = src[ind_z]
  ind = np.lexsort((y_l,x_l)); 
  if("SORT" in method):
   avg[ind_2d] = avg[ind_2d] + s_src[ind]
  elif("INTERP" in method):
   avg = avg + griddata((x_l,y_l),s_src, (x_2d,y_2d), method='nearest')
  else:
   raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')
 avg = avg/N;
 avg = np.asarray(avg, dtype=np.float64)
 #output.ShallowCopy(inputs[1])
 output.PointData.append(avg, varName0)

"""
  PF1.UpdatePipeline() 
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
print("--|| ALYA :: WORKSPACE VARIABLES AFTER SPAN AVERAGING", PF1.PointData.keys())

 ## CALCULTE STANDARD DEVIATION ###
if("SAVG" in mode): 
  print("--|| ALYA :: CALCULATE SPANWISE STANDARD DEVIATION")
  startTime = time.time()
  if("PINTERP" in method):
    print("--|| ALYA :: CREATING TRANSFORMATIONS")
    startTime = time.time()
    resample_transforms=list();
    data=list();
    for i in range(N):
      # create a new 'Transform'
      transform1 = Transform(Input=PF1,guiName="transform{}".format(i))
      try:
        resampleWithDataset1 = ResampleWithDataset(Input=case,Source=transform1)
      except:  
        resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=case,DestinationMesh=transform1)
      resample_transforms.append(resampleWithDataset1)
      # Properties modified on transform1.Transform
      transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
      #data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])
    print("--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
    
    print("--|| ALYA: STD-CALC USING A PROGRAMMABLE FILTER")
    startTime = time.time()
    PF2 = ProgrammableFilter(Input=[PF1]+resample_transforms)
    PF2.Script = \
"""
import os
import vtk
import numpy as np
from scipy.interpolate import griddata

caseName = 'naca'
fileName = './naca.pvd'
model = 'DNS'
nu = 0.000005
dim = '2D'
mode = 'SAVG'
method = 'INTERP'
file_fmt = 'BUDPVD'
xDec = 6
zDec = 6

inpFile = os.getcwd()+'/'+'inputFile.py'
exec(open(inpFile).read(), {'__file__': inpFile})

c = vtk.vtkMultiProcessController.GetGlobalController()
rank = c.GetLocalProcessId()

varFull = ['PRESS']

if(rank==0):
  t = inputs[0].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
  print("----|| ALYA : ARRAYS AT TIME %.3f " % (t), varFull)

N=len(inputs)-1;

for varName in varFull:
  avgVarName = "AV"+varName[0:3]
  varName0 = "SD"+varName[0:3]
 a_src = inputs[0].PointData[avgVarName]
 avg = 0.0*(inputs[0].PointData[varName])
 for i in range(N):
   s_src = inputs[i+1].PointData[varName]
   avg = avg + (s_src-a_src)**2
 avg = avg/N
 avg = np.asarray(avg, dtype=np.float64)
 output.PointData.append(avg,varName0)
"""
  else:
    PF2 = ProgrammableFilter(Input=[PF1,case])
    ### first input is the grid
    ### the rest of them are data to be averaged
    PF2.Script = \
  """
import os
import vtk
import numpy as np
from scipy.interpolate import griddata

caseName = 'naca'
fileName = './naca.pvd'
model = 'DNS'
nu = 0.000005
dim = '2D'
mode = 'SAVG'
method = 'INTERP'
file_fmt = 'BUDPVD'
xDec = 6
zDec = 6

inpFile = os.getcwd()+'/'+'inputFile.py'
exec(open(inpFile).read(), {'__file__': inpFile})

c = vtk.vtkMultiProcessController.GetGlobalController()
rank = c.GetLocalProcessId()

varFull = ['PRESS']

if(rank==0):
  t = inputs[0].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
  print("--|| ALYA :: ARRAYS AT TIME %.3f " % (t), varFull) 

try:
  d = dsa.WrapDataObject(inputs[1].GetBlock(0))
except:  
  d = dsa.WrapDataObject(inputs[1].VTKObject)
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=zDec)
## SLICE GEOMETRY
try:
  out = dsa.WrapDataObject(inputs[0].GetBlock(0))
except:  
  out = dsa.WrapDataObject(inputs[0].VTKObject)
x_2d = np.around(np.asarray(out.Points[:,0],dtype=np.double),decimals=6)
y_2d = np.around(np.asarray(out.Points[:,1],dtype=np.double),decimals=6)
ind_2d = np.lexsort((y_2d,x_2d));
#print("--|| ALYA CHECK :: SHAPE OF SLICE", np.shape(ind_2d))
z_vec = np.unique(z)
N = len(z_vec)

for varName in varFull:
 avgVarName = "AV"+varName[0:3]
 varName0 = "SD"+varName[0:3]
 avg = 0.0*np.asarray(out.PointData[avgVarName],dtype=np.double)
 src = np.asarray(d.PointData[varName],dtype=np.double)
 a_src = np.asarray(out.PointData[avgVarName],dtype=np.double)
 for n in range(N):
  ind_z = np.where(z == z_vec[n])
  x_l = x[ind_z];
  y_l = y[ind_z];
  s_src = src[ind_z]
  ind = np.lexsort((y_l,x_l)); 
  if("SORT" in method):
   avg[ind_2d] = avg[ind_2d] + (s_src[ind]-a_src[ind_2d])**2
  elif("INTERP" in method):
   s_src = griddata((x_l,y_l),s_src, (x_2d,y_2d), method='nearest')
   avg = avg + (s_src-a_src)**2
  else:
   raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE STD?')
 avg = avg/N;
 avg = np.asarray(avg, dtype=np.float64)
 output.PointData.append(avg, varName0)

"""
  PF2.UpdatePipeline() 
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  print("--|| ALYA :: WORKSPACE VARIABLES AFTER STD-DEV CALC", PF1.PointData.keys())

  ### APPEND SPECIFIC VARIABLES ###
  print("--|| ALYA :: APPEND AVERAGED VARIABLES")
  startTime = time.time()
  PF1 = AppendAttributes(Input=[PF1,PF2])
  PF1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  #PF1 = ProgrammableFilter(Input=PF1)
  #PF1.Script = \
  #"""
  #import os
  #import numpy as np
  #from scipy.interpolate import griddata
  #inpFile = os.getcwd()+'/'+'inputFile.py'
  #exec(open(inpFile).read(), {'__file__': inpFile})

  #varFull = inputs[0].PointData.keys()
  #for varName in varFull:
  # if("SAVG" in mode):
  #  print("----|| INFO :: APPENDING",varName)
  #  inp0 = inputs[0].PointData[varName]
  #  output.PointData.append(inp0,varName)
  # elif("FAVG" in mode):
  #  if not "_average" in varName:
  #   inp0 = inputs[0].PointData[varName]
  #   output.PointData.append(inp0,varName)
  #"""
  #PF1.UpdatePipeline() 
  print("--|| ALYA :: WORKSPACE VARIABLES AFTER APPEND", PF1.PointData.keys())

#  # EXTRACT DATA AT SPECIFIC TIMES
#  print("--|| ALYA :: PERFORM THE SAME ANALYSIS ON TIME DATA")
#  startTime = time.time()
#  dataSet = GroupTimeSteps(Input=PF1)
#  dataSet.UpdatePipeline()
#  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  caseVarNames = PF1.PointData.keys()
  indU = int([i for i, s in enumerate(caseVarNames) if 'AVVEL' in s][0]);
  indP = int([i for i, s in enumerate(caseVarNames) if 'AVPRE' in s][0]);
  if any("AVVE2" in s for s in caseVarNames):
   indXX = int([i for i, s in enumerate(caseVarNames) if 'AVVE2' in s][0]);
   indXY = int([i for i, s in enumerate(caseVarNames) if 'AVVXY' in s][0]);
   print("--|| ALYA :: CALCULATING R-STRESSES")
   startTime = time.time()
   # CALCULATE RStresses
   PF1 = Calculator(Input=PF1)
   PF1.ResultArrayName = "RS_II"
   PF1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                   % (caseVarNames[indXX],caseVarNames[indU],'X',caseVarNames[indU],'X',\
                   caseVarNames[indU],'Y',caseVarNames[indU],'Y',\
                   caseVarNames[indU],'Z',caseVarNames[indU],'Z')
   PF1.UpdatePipeline()
   PF1 = Calculator(Input=PF1)
   PF1.ResultArrayName = "RS_IJ"
   PF1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                   % (caseVarNames[indXY],caseVarNames[indU],'X',caseVarNames[indU],'Y',\
                   caseVarNames[indU],'Y',caseVarNames[indU],'Z',\
                   caseVarNames[indU],'X',caseVarNames[indU],'Z')
   PF1.UpdatePipeline()
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
   ## APPEND SPECIFIC VARIABLES ###
   print("--|| ALYA :: APPEND AVERAGED VARIABLES")
   startTime = time.time()
   PF1 = ProgrammableFilter(Input=PF1)
   PF1.Script = \
   """
   import numpy as np
   varFull = inputs[0].PointData.keys()
   for varName in varFull:
    if not varName in ["AVVE2","AVVXY"]:
     inp0 = inputs[0].PointData[varName]
     output.PointData.append(inp0,varName)
   """
   PF1.UpdatePipeline() 
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  # GRADIENT CALC
  print("--|| ALYA :: CALCULATING PRESS GRADIENT")
  startTime = time.time()
  PF1 = GradientOfUnstructuredDataSet(Input=PF1)
  PF1.ScalarArray = ['POINTS', caseVarNames[indP]]
  PF1.ComputeGradient = 1
  PF1.ResultArrayName = 'AVPGR'
  PF1.ComputeVorticity = 0
  PF1.VorticityArrayName = 'OMEGA'
  PF1.ComputeQCriterion = 0
  PF1.QCriterionArrayName = 'QCRIT'
  PF1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  print("--|| ALYA :: CALCULATING AVVEL GRADIENT, Q AND VORTICITY")
  startTime = time.time()
  PF1 = GradientOfUnstructuredDataSet(Input=PF1)
  PF1.ScalarArray = ['POINTS', caseVarNames[indU]]
  PF1.ComputeGradient = 1
  PF1.ResultArrayName = 'AVVGR'
  PF1.ComputeVorticity = 1
  PF1.VorticityArrayName = 'OMEGA'
  PF1.ComputeQCriterion = 1
  PF1.QCriterionArrayName = 'QCRIT'
  PF1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')


if('2D' in dim):
  #### WRITE THE OUTPUT FILE #########
  print("--|| ALYA: SAVING THE AVERAGED FILES")
  startTime = time.time()
  #savePath = casePath+"/AvgData_2D_"+mode+"_"+method+".vtm"
  savePath2 = casePath+"/AvgData_2D.csv"
  if("SAVG" in mode):
    if('PVD' in file_fmt): 
      savePath = casePath+"/AvgData_2D.pvd"
    else:
      savePath = casePath+"/AvgData_2D.vtm"
    SaveData(savePath, proxy=PF1, WriteTimeSteps=1)
  else:
    if('PVD' in file_fmt): 
      savePath = casePath+"/AvgData_2D.pvd"
    else:
      savePath = casePath+"/AvgData_2D.vtm"
    if(model == "TSTEP"):
      SaveData(savePath, proxy=PF1, WriteTimeSteps=1)
    else:  
      SaveData(savePath, proxy=PF1)
      SaveData(savePath2, proxy=PF1)
  print("----|| ALYA: FINAL STATISTICS FILE WRITTEN AS: ",savePath)
  print("--|| ALYA: FILE SAVED. TIME =",time.time()-startTime,'sec')
  ########### STREAMWISE AVERAGING ################
if('1D' in dim):
  print("--|| ALYA :: GENERATING SLICE FOR STREAMWISE AVERAGE")
  startTime = time.time()
  if("1DROT" in dim):
    slice1 = Slice(Input=PF1)
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
  
    Ny = int(slice1.GetDataInformation().GetNumberOfPoints())
    print("----|| ALYA :: WORKING WITH ",Nplane," PLANAR POINTS")
    
    N = int(Nplane/Ny)
    thMax = 2.0*np.pi; thMin = 0.0;
    thMid = np.pi/2
    xpos = np.arange(N)*(thMax-thMin)/(N-1)
    print("----|| ALYA: WORKING WITH %d THETA-PLANES" % (len(xpos)))
    print("----|| ALYA: DELTA-THETA = %f" % ((thMax-thMin)/(N-1)))
  else:
    print("--|| ALYA :: GENERATING SLICE FOR STREAMWISE AVERAGE")
    slice1 = Slice(Input=PF1)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    ## init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
    ## Properties modified on slice1.SliceType
    slice1.SliceType.Normal = [1.0, 0.0, 0.0]
    slice1.UpdatePipeline()
    
    Ny = slice1.GetDataInformation().GetNumberOfPoints()
    print("----|| ALYA :: WORKING WITH ",Ny," PLANAR POINTS")

    N = int(Nplane/Ny)
    xmid = (xmin+xmax)/2
    xpos = np.around(np.asarray(np.arange(N)*(xmax-xmin)/(N-1),dtype=np.double),decimals=xDec)
    print("----|| ALYA: WORKING WITH %d X-PLANES" % (len(xpos)))
    print("----|| ALYA: DELTA-X = %f" % ((xmax-xmin)/(N-1)))
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  if("PINTERP" in method):
   print("--|| ALYA :: CREATING TRANSFORMATIONS")
   startTime = time.time()
   resample_transforms=list();
   data=list();
   
   
   print("--|| ALYA: CREATING X TRANSFORMATIONS")
   startTime = time.time()
   for i in range(N):
   	# create a new 'Transform'
   	transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
   	if("1DROT" in dim):
   	  transform1.Transform.Rotate = [0.0, 0.0, xpos[i]-thMid]
   	else:
   	  transform1.Transform.Translate = [xpos[i]-xmid, 0.0, 0.0]
   	try:
   	  resampleWithDataset1=ResampleWithDataset(Input=PF1,Source=transform1)
   	except:  
   	  resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=PF1,DestinationMesh=transform1)
   	resample_transforms.append(resampleWithDataset1)
   print("--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
   
   HideAll()
   
   ## create a new 'Programmable Filter'
   print("--|| ALYA: X AVERAGING USING A PROGRAMMABLE FILTER")
   startTime = time.time()
   PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
   PF1.Script = \
   """
   import numpy as np
   
   varFull = []
   varFull = inputs[0].PointData.keys()
   print("----|| Alya - WORKING ON ORIGINAL ARRAYS::",varFull)
   N=len(inputs);
   
   for varName in varFull:
      varName0=str(varName[0:5])
      avg = (inputs[0].PointData[varName])
      for i in range(1,N):
          d = inputs[i].PointData[varName]
          avg = avg + d
      avg = avg/N
      output.PointData.append(avg,varName0)
   """
   PF1.UpdatePipeline()
   print("--|| ALYA: STREAMWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec')
  else:
   print("--|| ALYA :: STREAMWISE AVERAGING")
   startTime = time.time()
   PF1 = ProgrammableFilter(Input=[PF1,slice1])
   PF1.Script = \
   """
   import numpy as np
   t = inputs[0].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
   varFull = inputs[0].PointData.keys()
   print("--|| ALYA :: CALCULATING FOR",varFull," AT T=",t) 
   if('ensi' in fileName):
    d = dsa.WrapDataObject(inputs[0].GetBlock(0))
   elif('pvd' in fileName):
    d = dsa.WrapDataObject(inputs[0])
   x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=xDec)
   y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
   ## SLICE GEOMETRY
   out = dsa.WrapDataObject(inputs[1].GetBlock(0))
   y_2d = np.around(np.asarray(out.Points[:,1],dtype=np.double),decimals=6)
   ind_2d = np.argsort(y_2d);
   print("--|| ALYA CHECK :: SHAPE OF SLICE", np.shape(ind_2d))
   
   for varName in varFull:
    print("--|| ALYA :: CALCULATING FOR",varName)
    varName0 = varName[0:5]
    avg = 0.0*np.asarray(out.PointData[varName],dtype=np.double)
    src = np.asarray(d.PointData[varName],dtype=np.double)
    for n in range(N):
     ind_z = np.where(x == xpos[n])
     #print("--|| ALYA CHECK :: AT n,X=",n,xpos[n],"SHAPE OF MASTER DATASET", np.shape(ind_z))
     y_l = y[ind_z];
     s_src = src[ind_z]
     if("SORT" in method):
      ind = np.argsort(y_l); 
      avg[ind_2d] = avg[ind_2d] + s_src[ind]
     elif("INTERP" in method):
      avg = avg + griddata((y_l),s_src, (y_2d), method='linear')
     else:
      raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')
    avg = avg/N;
    avg = np.asarray(avg, dtype=np.float64)
    output.ShallowCopy(inputs[1].VTKObject)
    output.PointData.append(avg, varName0)
   
   """
   PF1.UpdatePipeline() 
   print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
   
  #### write
  print("--|| ALYA: SAVING THE AVERAGED FILES")
  startTime = time.time()
  savePath = casePath+"/AvgData_1D.vtm"
  SaveData(savePath, proxy=PF1)
  savePath = casePath+"/AvgData_1D.csv"
  SaveData(savePath, proxy=PF1)
  print("----|| ALYA: 1D STATISTICS FILE WRITTEN AS: ",savePath)
  print("--|| ALYA: FILE SAVED. TIME =",time.time()-startTime,'sec')
