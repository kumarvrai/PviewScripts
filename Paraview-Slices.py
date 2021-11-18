#### import the simple module from the paraview
import os
import time
import operator
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from paraview import vtk
from paraview import numpy_support
from paraview.vtk.numpy_interface import dataset_adapter as dsa
import paraview.vtk.util.numpy_support as vnp
from vtk.util import numpy_support
from paraview.simple import *

niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'
niaHome = '/home/u/ugo/kvishal/'
#Connect("localhost", 11111)

# INITIALIZE VARIABLES
airfoil = '0012'
#IF = 'AvgData_2D.vtm'
IF = 'naca.ensi.case'
OD = './0.Results/1.SliceResults/'

tSlice = True; nSlice = True; zSlice=False
d = 0.025
x_loc = [0.005, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4]
n_loc = [0.002, 0.004, 0.006]
z_loc = [0.0, 0.05, 0.1, 0.15, 0.2]



if not os.path.exists(OD):
  os.makedirs(OD)

# LOAD AIRFOIL SPECIFIC FILES
print('--|| ALYA INITIALIZING')
if('0012' in airfoil):
  fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-UP.txt'
  fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-DOWN.txt'
  sliceLoc = niaHome+'/1.PostProc/1.Airfoil/3.PviewScripts/1.0012-Slices/'
elif('trip' in airfoil):
  fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-UP.txt'
  fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-DOWN.txt'
  sliceLoc = niaHome+'/1.PostProc/1.Airfoil/3.PviewScripts/2.4412-Slices/'
elif('rough' in airfoil):
  fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-UP.txt'
  fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-DOWN.txt'
  sliceLoc = niaHome+'/1.PostProc/1.Airfoil/3.PviewScripts/2.4412-Slices/'
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

## INTERPOLATION ON SUCTION SIDE
#thAirL = np.arctan2(np.diff(coordAirL[:,1]),np.diff(coordAirL[:,0]))
#F1 = interp1d(coordAirL[:,0],coordAirU[:,1])
#airLen = len(coordAirL)
#coordMid = 0.5*(coordAirU[airLen:2*airLen-2,0]+coordAirU[airLen+1:2*airLen-1,0])
#Fth1 = interp1d(coordMid,thAirL)

#
# EXTRACT LOCATIONS
y_loc = F0(x_loc)
th_loc = Fth0(x_loc)
center = np.empty((0,6),dtype='float')
for ii in range(0,len(x_loc)):
  x0 = x_loc[ii]; y0 = y_loc[ii];
  m = np.tan(th_loc[ii]);
  rhs = d**2/(1+float(1.0/m**2));
  if(m>0.0):
    x1 = x0 - np.sqrt(rhs);
  else:
    x1 = x0 + np.sqrt(rhs);
  y1 = y0-(x1-x0)/m;
  center = np.vstack((center,np.array([x0,y0,np.cos(th_loc[ii]),np.sin(th_loc[ii]),x1,y1])))

print('--|| ALYA :: AIRFOIL VARIABLES LOADED')

# Read the Source
source = OpenDataFile(IF)
source.UpdatePipeline()


if '.ensi.case' in IF:
  caseVarNames = source.PointArrays
elif '.vtm' in IF:
  caseVarNames = source.PointArrayStatus
else:
  raise ValueError('--|| ALYA ERROR :: FILE NOT PROVIDED.')

print('--|| ALYA ARRAYS ::',caseVarNames)

indU = int([i for i, s in enumerate(caseVarNames) if 'AVVEL' in s][0]); 
indXX = int([i for i, s in enumerate(caseVarNames) if 'AVVE2' in s][0]); 
indXY = int([i for i, s in enumerate(caseVarNames) if 'AVVXY' in s][0]); 
indu = int([i for i, s in enumerate(caseVarNames) if 'VELOC' in s][0]); 


# create a new 'Extract Surface'
ES1 = ExtractSurface(Input=source)
# rename source object
RenameSource('Airfoil-Surface', ES1)
ES1.UpdatePipeline() 

# create a new 'Clip'
CLP1 = Clip(Input=ES1)
# rename source object
RenameSource('Airfoil-Surface-Clip', CLP1)
# Properties modified on clip
CLP1.ClipType = 'Box'
# Properties modified on clip1_1.ClipType
CLP1.ClipType.Position = [0.0, 0.0, 0.0]
CLP1.ClipType.Scale = [0.5, 0.5, 1.0]
CLP1.Invert = 1
CLP1.UpdatePipeline() 

print('--|| ALYA :: INDICES OF U UXX UXY',indU,indXX,indXY)

TS = TemporalStatistics(Input=source)
# Properties modified on temporalStatistics1
TS.ComputeMinimum = 0
TS.ComputeMaximum = 0
TS.ComputeStandardDeviation = 0
TS.UpdatePipeline()
# CALCULATE Fluctuations
CAL1 = PythonCalculator(Input=[source,TS])
CAL1.ArrayName = "u"
CAL1.Expression = "inputs[0].PointData['%s'] - inputs[1].PointData['%s'+'_average']" % (caseVarNames[indu],caseVarNames[indU])
RenameSource('u', CAL1)
CAL1.UpdatePipeline()

# CALCULATE RStresses
CAL1 = Calculator(Input=CAL1)
CAL1.ResultArrayName = "RS_II"
CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                % (caseVarNames[indXX],caseVarNames[indU],'X',caseVarNames[indU],'X',\
                  caseVarNames[indU],'Y',caseVarNames[indU],'Y',\
                 caseVarNames[indU],'Z',caseVarNames[indU],'Z')
RenameSource('RS-II', CAL1)
CAL1.UpdatePipeline() 
CAL1 = Calculator(Input=CAL1)
CAL1.ResultArrayName = "RS_IJ"
CAL1.Function = "%s - %s_%s*%s_%s*iHat - %s_%s*%s_%s*jHat - %s_%s*%s_%s*kHat" \
                % (caseVarNames[indXY],caseVarNames[indU],'X',caseVarNames[indU],'Y',\
                  caseVarNames[indU],'Y',caseVarNames[indU],'Z',\
                  caseVarNames[indU],'X',caseVarNames[indU],'Z')
RenameSource('RS-IJ', CAL1)
CAL1.UpdatePipeline()


############# SLICES ##########################

nNames=[]
nNames.append(['Airfoil-Surface-Clip'])

if nSlice:
  print('--|| ALYA EXTRACTS :: ', x_loc)
  for n in range(0,len(x_loc)):
    # Generate Normal Slices on Airfoil Surface
    slice1 = Slice(Input=CAL1)
    strName = 'n-slice-%s' % (x_loc[n])
    RenameSource(strName, slice1) 
    slice1.SliceType.Origin = [center[n,0], center[n,1], 0.1]
    slice1.SliceType.Normal = [center[n,2], center[n,3], 0.0]
    slice1.UpdatePipeline()
    # create a new 'Clip'
    clip = Clip(Input=slice1)
    # rename source object
    strName = 'nSlice-3D-%s' % (x_loc[n])
    RenameSource(strName, clip) 
    # Properties modified
    clip.ClipType.Origin = [center[n,4], center[n,5], 0.1]
    clip.ClipType.Normal = [-center[n,3], center[n,2], 0.0]
    clip.Invert = 1
    clip.UpdatePipeline()
    nNames.append([strName])
    #srcDisp = Show(clip, RV1)
    # create a new 'Calculator'
    CAL1 = Calculator(Input=clip)
    # rename source object
    strName = 'nSlice-%s' % (x_loc[n])
    RenameSource(strName, CAL1)
    # Properties modified on calculator1
    CAL1.CoordinateResults = 1
    # Properties modified on calculator1
    CAL1.Function='(sqrt((coordsX-%s)^2+((coordsY-%s)^2)))*jHat+coordsZ*iHat'%(center[n,0],center[n,1])
    CAL1.UpdatePipeline()
    # Collect names for saving
    nNames.append([strName])


###### Plane Parallel Visualization #############
if tSlice:
 tNames=[]
 for n in range(len(n_loc)): 
   f = 'normPlane-' + str(n_loc[n])
   f = f+'.csv'
   print('--|| ALYA ::  READING PLANE', f)
   sliceCase = OpenDataFile(sliceLoc+f)
   RenameSource(f, sliceCase)
   sliceCase.UpdatePipeline()
   
   TtP = TableToPoints(Input=sliceCase)
   strName = 'Table-%s' % n_loc[n]
   RenameSource(strName, TtP)
   # Properties modified on tableToPoints1
   TtP.XColumn = 'x'
   TtP.YColumn = 'z'
   TtP.ZColumn = 'y'
   TtP.a2DPoints = 1
   TtP.UpdatePipeline()

   
   # create a new 'Resample With Dataset'
   RwD1 = ResampleWithDataset(Input=source,Source=TtP)
   strName = 'Resample-%s' % n_loc[n]
   RenameSource(strName, RwD1)
   RwD1.UpdatePipeline()
   
   #AA = AppendAttributes(Input=[RwD1,TtP])
   #strName = 'AppendBoth-%s' % n_loc[n]
   #RenameSource(strName, AA)
   #AA.UpdatePipeline()
 
   ## create a new 'Calculator'
   #CAL1 = Calculator(Input=AA)
   #strName = 'calcCoord-%s' % n_loc[n]
   #RenameSource(strName, CAL1)
   ## Properties modified on calculator1
   #CAL1.CoordinateResults = 1
   #CAL1.Function = 's*iHat+(coordsZ)*jHat'
   #CAL1.UpdatePipeline()
   
   # create a new 'Delaunay 2D'
   d2D1 = Delaunay2D(Input=RwD1)
   strName = 'tSlice-%s' % n_loc[n]
   RenameSource(strName, d2D1)
   d2D1.UpdatePipeline()
   tNames.append([strName])

#if zSlice:
#  for n in range(0,len(z_loc)):
#    # Generate Normal Slices on Airfoil Surface
#    slice1 = Slice(Input=clipU)
#    strName = 'z-slice-%s' % (z_loc[n])
#    RenameSource(strName, slice1) 
#    slice1.SliceType.Origin = [0.5, 0.0, z_loc[n]]
#    slice1.SliceType.Normal = [0.0, 0.0, 1.0]
#    slice1.UpdatePipeline()
#    # create a calculator
#    CAL1 = Calculator(Input=slice1)
#    CAL1.ResultArrayName = 'XYZ'
#    CAL1.Function = 'coords'
#    CAL1.UpdatePipeline()
#    # create a new 'Clip'
#    PC1 = ProgrammableFilter(Input=CAL1)
#    PC1.Script=\
#    """
#import numpy as np
#from vtk.util import numpy_support
#import vtk.numpy_interface.dataset_adapter as dsa
#
#
#def closest_node(node, nodes):
# nodes = np.asarray(nodes)
# deltas = nodes - node
# dist_2 = np.einsum('ij,ij->i', deltas, deltas)
# return np.argmin(dist_2)
#
#inData = dsa.WrapDataObject(inputs[0])
#nodes  = inData.PointData["XYZ"]
#nodes = numpy_support.vtk_to_numpy(nodes.Arrays[0])
#scalF  = np.zeros(1,len(nodes))
#
#for i in range(len(nodes)):
#  scalF[i] = closes_node(nodes[i], coordAirU)
#output.PointData.append(scalF,"D")
#    PC1.UpdatePipeline()
#    #clip = Clip(Input=slice1)
#    ## rename source object
#    #strName = 'z-slice-clip-%s' % (z_loc[n])
#    #RenameSource(strName, clip) 
#    ## Properties modified
#    #clip.UpdatePipeline()
#    ## create a new 'Calculator'
#    #CAL1 = Calculator(Input=clip)
#    ## rename source object
#    #strName = 'CoordCalc-%s' % (x_loc[n])
#    #RenameSource(strName, CAL1)
#    ## Properties modified on calculator1
#    #CAL1.CoordinateResults = 1
#    ## Properties modified on calculator1
#    #CAL1.Function='(sqrt((coordsX-%s)^2+(coordsY-%s)^2))*jHat+coordsZ*iHat'%(center[n,0],center[n,1])
#    #CAL1.UpdatePipeline()
#    ## get active view
#"""

for i in range(len(nNames)):
  name = nNames[i][0]+'.vtm'
  print('--||ALYA :: WRITING',nNames[i][0])
  savePath = OD+name
  src = FindSource(nNames[i][0])
  SaveData(savePath, proxy=src)
#  SaveData(savePath, proxy=src, Writealltimestepsasfileseries=0,
#  DataMode='Binary', HeaderType='UInt64', EncodeAppendedData=0,
#  CompressorType='None')

for i in range(len(tNames)):
  name = tNames[i][0]+'.vtk'
  print('--||ALYA :: WRITING',tNames[i][0])
  savePath = OD+name
  src = FindSource(tNames[i][0])
  SaveData(savePath, proxy=src)
#  SaveData(savePath, proxy=src, Writealltimestepsasfileseries=0,
#  DataMode='Binary', HeaderType='UInt64', EncodeAppendedData=0,
#  CompressorType='None')

print('--||ALYA :: SLICES HAVE BEEN SAVED IN',OD)
