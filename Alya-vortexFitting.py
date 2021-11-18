import sys
import os.path
import argparse
import os
import sys
import time
import numpy as np
import vtk
from paraview.simple import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

from vortexfitting import fitting
#from vortexfitting import schemes
#from vortexfitting import detection
#from vortexfitting import output
#from vortexfitting import classes

caseName	= sys.argv[1]
mode		= sys.argv[2]

casePath = os.getcwd()

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ---- READ DATA ----#
print "--|| ALYA :: READING ALYA-AVERAGED ARRAYS"
startTime = time.time()

fileName = caseName+'.ensi.case'

case = OpenDataFile(fileName)
case.UpdatePipeline()
print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'

if("NPLANE" in mode):
 # INITIALIZE VARIABLES
 d = 0.1
 x_loc = [0.9]
  
 if('chan' in caseName):
  print('----|| CREATING NORMAL PLANES FOR CHAN','AT X/C=', x_loc)
  caseList = []
  for n in range(len(x_loc)):
    # Generate Normal Slices on Airfoil Surface
    slice1 = Slice(Input=case)
    slice1.SliceType.Origin = [x_loc[n], 0, 0.0]
    slice1.SliceType.Normal = [1.0, 0.0, 0.0]
    slice1.UpdatePipeline()
    caseList.append(slice1)
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

 else:
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

##---------------------------------------------#    
print("--|| ALYA :: VORTEX-FITTING ON A 2D PLANE")
startTime = time.time()
PF1 = ProgrammableFilter(Input=caseList)
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
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
vortex_field = np.asarray(d.PointData["SWIRL"],dtype=np.double)

# ---- PEAK DETECTION ----#
peaks = fitting.find_peaks(vortex_field, threshold=0.0, box_size=0.1)
print('--|| ALYA : VORTICES FOUND = ', len(peaks[0]))

# ---- PEAKS DIRECTION OF ROTATION ----#
vortices_counterclockwise, vortices_clockwise = fitting.direction_rotation(vorticity, peaks)

# ---- MODEL FITTING ----#
vortices = list()

vortices = fitting.get_vortices(vfield, peaks, vorticity, args.rmax, args.correlation_threshold)
print('--|| ALYA : Accepted vortices=', len(vortices))

output.ShallowCopy(inputs[0].VTKObject)
output.PointData.append(avg, varName0)

"""
PF1.UpdatePipeline() 
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')


