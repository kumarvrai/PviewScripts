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
from scipy.interpolate import interp1d

casePath = os.getcwd()

airfoil = '4412'; side = 'SS'

caseName = sys.argv[1]
model    = sys.argv[2]
nu       = float(sys.argv[3])
dim      = sys.argv[4]
mode     = sys.argv[5]
method     = sys.argv[6]
zDec = 6; xDec = 6

paraview.simple._DisableFirstRenderCameraReset()
#
#-------------READ ALYA ARRAYS -----------------#
print("--|| ALYA :: READING ALYA-AVERAGED ARRAYS")
startTime = time.time()
fileName = caseName+'.ensi.case'
case = OpenDataFile(fileName)
case.PointArrays = ['VELOC']
case.UpdatePipeline()

reader = GetActiveSource()
view = GetActiveView()
times = reader.TimestepValues

avgCase = OpenDataFile('AvgData_3D.vtm')
avgCase.PointArrayStatus = ['AVVEL']
avgCase.UpdatePipeline()
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

#-------------READ ALYA ARRAYS -----------------#
print("--|| ALYA :: CREATING PLANES TO CALCULATE STATISTICS")
startTime = time.time()

if("NPLANE" in mode):
  # INITIALIZE VARIABLES
  d = 0.1
  x_loc = [0.9]
  
  # LOAD AIRFOIL SPECIFIC FILES
  #baseDir = '/home/kvishal/projects/rrg-ugo-ac/kvishal/1.post_process/1.airfoil/3.PviewScripts/'
  baseDir = '/home/kvishal/1.post_process/1.airfoil/3.PviewScripts/'
  if('0012' in airfoil):
    fAirU = baseDir+'/1.0012-Slices/naca0012-UP.txt'
    fAirL = baseDir+'/1.0012-Slices/naca0012-DOWN.txt'
  elif('4412' in airfoil):
    fAirU = baseDir+'/2.4412-Slices/naca4412-UP.txt'
    fAirL = baseDir+'/2.4412-Slices/naca4412-DOWN.txt'
  else:
    raise ValueError('--|| ALYA ERROR :: FILE NOT PROVIDED.')
  
  
  coordAirU = np.loadtxt(fAirU, delimiter=',')
  coordAirL = np.loadtxt(fAirL, delimiter=',')
  
  # INTERPOLATION ON SUCTION SIDE
  thAirU = np.arctan2(np.diff(coordAirU[:,1]),np.diff(coordAirU[:,0]))
  F0 = interp1d(coordAirU[:,0],coordAirU[:,1])
  airLen = len(coordAirU)
  coordMid = 0.5*(coordAirU[0:airLen-1,0]+coordAirU[1:airLen,0])
  Fth0 = interp1d(coordMid,thAirU)
  
  # INTERPOLATION ON PRESSURE SIDE
  thAirL = np.arctan2(np.diff(coordAirL[:,1]),np.diff(coordAirL[:,0]))
  F1 = interp1d(coordAirL[:,0],coordAirL[:,1])
  airLen = len(coordAirL)
  coordMid = 0.5*(coordAirL[0:airLen-1,0]+coordAirL[1:airLen,0])
  Fth1 = interp1d(coordMid,thAirL)
  
  #
  # EXTRACT LOCATIONS
  if('SS' in side):
   y_loc = F0(x_loc)
   th_loc = Fth0(x_loc)
  elif('PS' in side):
   y_loc = F1(x_loc)
   th_loc = Fth1(x_loc)
  center = np.empty((0,6),dtype='float')
  for ii in range(0,len(x_loc)):
    x0 = x_loc[ii]; y0 = y_loc[ii];
    m = np.tan(th_loc[ii]);
    rhs = d**2/(1+float(1.0/m**2));
    if('SS' in side):
      if(m>0.0):
        x1 = x0 - np.sqrt(rhs);
      else:
        x1 = x0 + np.sqrt(rhs);
      y1 = y0-(x1-x0)/m;
    elif('PS' in side):
      if(m<0.0):
        x1 = x0 - np.sqrt(rhs);
      else:
        x1 = x0 + np.sqrt(rhs);
      y1 = y0-(x1-x0)/m;
    center = np.vstack((center,np.array([x0,y0,np.cos(th_loc[ii]),np.sin(th_loc[ii]),x1,y1])))

  print('----|| CREATING NORMAL PLANES ON', side,'AT X/C=', x_loc)
  caseList = []
  for n in range(len(x_loc)):
    # Generate Normal Slices on Airfoil Surface
    slice1 = Slice(Input=case)
    slice1.SliceType.Origin = [center[n,0], center[n,1], 0.1]
    slice1.SliceType.Normal = [center[n,2], center[n,3], 0.0]
    slice1.UpdatePipeline()
    # Generate Normal Slices on Airfoil Surface
    slice2 = Slice(Input=avgCase)
    slice2.SliceType.Origin = [center[n,0], center[n,1], 0.1]
    slice2.SliceType.Normal = [center[n,2], center[n,3], 0.0]
    slice2.UpdatePipeline()
    # Resample Average values
    slice3 = ResampleWithDataset(Input=slice2,Source=slice1)
    slice3.UpdatePipeline()
    # Append Attributes
    slice1 = AppendAttributes(Input=[slice1,slice3])
    slice1.UpdatePipeline()
    # create a new 'Clip'
    clip = Clip(Input=slice1)
    clip.ClipType.Origin = [center[n,4], center[n,5], 0.1]
    clip.ClipType.Normal = [-center[n,3], center[n,2], 0.0]
    if('SS' in side):
     clip.Invert = 1
    elif('PS' in side):
     clip.Invert = 0
    clip.UpdatePipeline()
    # create a new 'Clip'
    clip = Clip(Input=clip)
    clip.ClipType.Origin = [0, 0, 0.1]
    clip.ClipType.Normal = [0, 1, 0.0]
    if('SS' in side):
     clip.Invert = 0
    elif('PS' in side):
     clip.Invert = 1
    clip.UpdatePipeline()
    caseList.append(clip)
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print("--|| ALYA :: APPEND ALL SLICES TOGTHER.")
  startTime = time.time()
  sliceCase = AppendAttributes(Input=caseList)  
  sliceCase.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  # EXTRACT DATA AT SPECIFIC TIMES
  print("--|| ALYA :: EXTRACT TIME DATA")
  startTime = time.time()
  dataSet = GroupTimeSteps(Input=sliceCase)  
  dataSet.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print("--|| ALYA :: PERFORM STATISTICAL ANALYSIS AND APPEND DATA")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=[dataSet,sliceCase])
  PF1.Script = \
  """
print("----|| ALYA :: MANUPULATE DATA")
startTime = time.time()

V = dsa.WrapDataObject(inputs[0].GetBlock(0)).PointData['VELOC'].Arrays[0]
V = vtk_to_numpy(V)

N = np.size(times)
L = np.shape(V)[0]

print('----|| ALYA : TEMPORAL SIZE IS ',N, 'SPATIAL SIZE IS ',L)

fields = np.zeros([N,L,3])
count_flds = np.zeros([L,3],dtype=np.double)

print ('----|| ALYA : READING TIME: FROM {} TO {}'.format(times[0],times[-1]))
for i in range(N):
 t = times[i]
 d = dsa.WrapDataObject(inputs[0].GetBlock(i))
 d = vtk_to_numpy(d.PointData["VELOC"].Arrays[0])
 fields[i]=np.copy(d)
 
print('--|| ALYA : CALCULATING FLUCTUATIONS')
fields_avg = dsa.WrapDataObject(inputs[0].GetBlock(0))
fields_avg = vtk_to_numpy(fields_avg.PointData["AVVEL"].Arrays[0])
fields = fields - fields_avg
print("----|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
eps = 1e-4
count_flds = np.asarray(np.count_nonzero(abs(fields)>eps,axis=0),np.double)/N 

output.ShallowCopy(inputs[1].VTKObject)
output.PointData.append(count_flds,"GAMMA")

  """ 
  PF1.UpdatePipeline()
  #### write
  print("--|| ALYA: SAVING THE AVERAGED FILES")
  startTime = time.time()
  savePath = casePath+"/AvgPlaneData_2D.vtm"
  SaveData(savePath, proxy=PF1)
  print("----|| ALYA: PLANAR STATISTICS FILE WRITTEN AS: ",savePath)
  print("--|| ALYA: FILE SAVED. TIME =",time.time()-startTime,'sec')
