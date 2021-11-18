#### import the simple module from the paraview
import os
import time
import operator
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#Connect("localhost", 11111)

# LOAD AIRFOIL SPECIFIC FILES
print('--|| ALYA INITIALIZING')
niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'
niaHome = '/home/u/ugo/kvishal/'

fileLoc = niaScrh+'/5.BGQ-DATA/NACA-0012/2.3D/2.LES/1.Re-5E4-alpha-8/3.WALE/1.run1/1.t-0-13/'
airfoilCoord = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012.txt'

coordAirfoil = np.loadtxt(airfoilCoord)

tanAirfoil = np.arctan2(np.diff(coordAirfoil[:,1]),np.diff(coordAirfoil[:,0]))
airLen = np.floor_divide(len(coordAirfoil),2)
F0 = interp1d(coordAirfoil[0:airLen,0],coordAirfoil[0:airLen,1])
F1 = interp1d(coordAirfoil[airLen:2*airLen-1,0],coordAirfoil[airLen:2*airLen-1,1])
Fth0 = interp1d(coordAirfoil[0:airLen,0],tanAirfoil[0:airLen])
Fth1 = interp1d(coordAirfoil[airLen:2*airLen-1,0],tanAirfoil[airLen:2*airLen-1])
#xAir,zAir = np.meshgrid(coordAirfoil[:,0],np.linspace(0,0.2,10))
#yAir,zAir = np.meshgrid(coordAirfoil[:,1],np.linspace(0,0.2,10))
#
# EXTRACT LOCATIONS
x_loc = [0.05, 0.7, 0.9]
y_loc = F0(x_loc)
th_loc = Fth0(x_loc)
d = 0.05
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
print('--|| ALYA EXTRACTS :: ', x_loc)

# create a new 'XML MultiBlock Data Reader'
case = OpenDataFile(fileLoc+'/AvgData_3D.vtm')
case.PointArrayStatus = ['AVVEL_average', 'AVPRE_average', 'AVTAN_average', 'AVVE2_average', 'AVVXY_average']
RenameSource('Alya-CaseFile', case)
# update the view to ensure updated data information
case.UpdatePipeline()

#exit()
# create a new 'Gradient Of Unstructured DataSet'
GOUD1 = GradientOfUnstructuredDataSet(Input=case)
# Properties modified on gradientOfUnstructuredDataSet1
GOUD1.ScalarArray = ['POINTS', 'AVVEL_average']
GOUD1.ComputeGradient = 0
GOUD1.ComputeVorticity = 1
GOUD1.VorticityArrayName = 'Omega'
RenameSource('GradientCalc', GOUD1)
GOUD1.UpdatePipeline()

# create a new 'Extract Surface'
ES1 = ExtractSurface(Input=GOUD1)
# rename source object
RenameSource('Airfoil-Surface', ES1)
ES1.UpdatePipeline()

# create a new 'Clip'
clip = Clip(Input=ES1)
# rename source object
RenameSource('Airfoil-Surface-Clip', clip)
# Properties modified on clip
clip.ClipType = 'Box'
# Properties modified on clip1_1.ClipType
clip.ClipType.Position = [0.5, 0.0, 0.0]
clip.ClipType.Scale = [0.05, 0.05, 1.0]
clip.UpdatePipeline()

# create a new 'Feature Edges'
FE1 = FeatureEdges(Input=ES1)
# rename source object
RenameSource('Airfoil-Edges', FE1)
FE1.UpdatePipeline()

# create a new 'Clip'
clipV = Clip(Input=GOUD1)
# rename source object
RenameSource('Airfoil-Small-Volume-Clip', clipV)
# Properties modified on clip1
clipV.ClipType = 'Box'
# Properties modified on clip1.ClipType
clipV.ClipType.Position = [0.5, 0.0, 0.0]
clipV.ClipType.Scale = [0.05, 0.025, 1.0]
clipV.UpdatePipeline()


# create a new 'Clip'
clipU = Clip(Input=clipV)
# rename source object
RenameSource('Airfoil-UpperHalf-Clip', clipU)
# Properties modified on clip1_2.ClipType
clipU.ClipType.Origin = [0.0, 0.0, 0.1]
clipU.ClipType.Normal = [0.0, 1.0, 0.0]
clipU.Invert = 0
clipU.UpdatePipeline()

# create a new 'Clip'
clipL = Clip(Input=clipV)
# rename source object
RenameSource('Airfoil-LowerHalf-Clip', clipL)
# Properties modified on clip1_2.ClipType
clipL.ClipType.Origin = [0.0, 0.0, 0.1]
clipL.ClipType.Normal = [0.0, 1.0, 0.0]
clipL.Invert = 1
clipL.UpdatePipeline()

###### VISUALIZE RESULTS IN PARAVIEW ################
layOut1 = GetLayout()
row = 2
col = 2
# split cell
for i in range(row-1):
 layOut1.SplitHorizontal(2*i, float(1.0/(row-i)))
 rv = CreateView('RenderView')
for i in range(1,row+1):
 for j in range(col):
  layOut1.SplitVertical(i+row*j, float(1.0/(col-j)))
  rv = CreateView('RenderView')
views = GetRenderViews()
# find source
src1 = FindSource('Airfoil-Surface-Clip')
# get active view
RV1 = views[0]
# uncomment following to set a specific view size
#RV1.ViewSize = [int(1920/col), int(1080/row)]
# show data in view
srcDisp = Show(src1, RV1)
srcDisp.DiffuseColor = [0.0, 0.6666666666666666, 0.0]
# find source
airfoilEdges = FindSource('Airfoil-Edges')
# show data in view
srcDisp = Show(airfoilEdges, RV1)
srcDisp.DiffuseColor = [0.0, 0.0, 0.0]
srcDisp.PointSize = 3.0
srcDisp.LineWidth = 3.0

for n in range(0,len(x_loc)):
  # Generate Normal Slices on Airfoil Surface
  slice1 = Slice(Input=clipU)
  strName = 'Slice-%s' % (x_loc[n])
  RenameSource(strName, slice1) 
  slice1.SliceType.Origin = [center[n,0], center[n,1], 0.1]
  slice1.SliceType.Normal = [center[n,2], center[n,3], 0.0]
  slice1.UpdatePipeline()
  # create a new 'Clip'
  clip = Clip(Input=slice1)
  # rename source object
  strName = 'Normal-Clip-%s' % (x_loc[n])
  RenameSource(strName, clip) 
  # Properties modified
  clip.ClipType.Origin = [center[n,4], center[n,5], 0.1]
  clip.ClipType.Normal = [-center[n,3], center[n,2], 0.0]
  clip.Invert = 1
  clip.UpdatePipeline()
  srcDisp = Show(clip, RV1)
  # create a new 'Calculator'
  CAL1 = Calculator(Input=clip)
  # rename source object
  strName = 'CoordCalc-%s' % (x_loc[n])
  RenameSource(strName, CAL1)
  # Properties modified on calculator1
  CAL1.CoordinateResults = 1
  # Properties modified on calculator1
  CAL1.Function='(sqrt((coordsX-%s)^2+(coordsY-%s)^2))*jHat+coordsZ*iHat'%(center[n,0],center[n,1])
  CAL1.UpdatePipeline()
  # get active view
  RV2 = views[n+1]
  #RV2.ViewSize = [int(1920/col), int(1080/row)]
  srcDisp = Show(CAL1, RV2)




# layout1.SplitVertical(i+1, float(1.0/col))

### Create a new 'Render View'
##renderView3 = CreateView('RenderView')
##renderView3.ViewSize = [1280, 960]
##renderView3.AxesGrid = 'GridAxes3DActor'
##renderView3.StereoType = 0
##renderView3.Background = [0.32, 0.34, 0.43]
##renderView3.OSPRayMaterialLibrary = materialLibrary1
##
### init the 'GridAxes3DActor' selected for 'AxesGrid'
##renderView3.AxesGrid.Visibility = 1
##renderView3.AxesGrid.XTitle = 'X'
##renderView3.AxesGrid.YTitle = 'Y'
##renderView3.AxesGrid.ZTitle = 'Z'
##renderView3.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
##renderView3.AxesGrid.XTitleFontFamily = 'Times'
##renderView3.AxesGrid.XTitleFontFile = ''
##renderView3.AxesGrid.XTitleFontSize = 24
##renderView3.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
##renderView3.AxesGrid.YTitleFontFamily = 'Times'
##renderView3.AxesGrid.YTitleFontFile = ''
##renderView3.AxesGrid.YTitleFontSize = 24
##renderView3.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
##renderView3.AxesGrid.ZTitleFontFamily = 'Times'
##renderView3.AxesGrid.ZTitleFontFile = ''
##renderView3.AxesGrid.ZTitleFontSize = 24
##renderView3.AxesGrid.GridColor = [0.0, 0.0, 0.0]
##renderView3.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
##renderView3.AxesGrid.XLabelFontFamily = 'Times'
##renderView3.AxesGrid.XLabelFontFile = ''
##renderView3.AxesGrid.XLabelFontSize = 20
##renderView3.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
##renderView3.AxesGrid.YLabelFontFamily = 'Times'
##renderView3.AxesGrid.YLabelFontFile = ''
##renderView3.AxesGrid.YLabelFontSize = 20
##renderView3.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
##renderView3.AxesGrid.ZLabelFontFamily = 'Times'
##renderView3.AxesGrid.ZLabelFontFile = ''
##renderView3.AxesGrid.ZLabelFontSize = 20
##
### place view in the layout
##layout1.AssignView(13, renderView3)

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
