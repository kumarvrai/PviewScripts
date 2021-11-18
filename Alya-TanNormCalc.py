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

def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.sum((array - value)**2,axis=1)).argmin()
    return idx

niaHome = '/home/kvishal/projects/rrg-ugo-ac/kvishal/0.pre_process/1.pre_process_files/'

d = dsa.WrapDataObject(inputs[0].GetBlock(0))
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)

Nx = np.shape(x)[0]

# INITIALIZE VARIABLES
airfoil = '4412';

# LOAD AIRFOIL SPECIFIC FILES
print('--|| ALYA INITIALIZING')
if('0012' in airfoil):
  fAirU = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-UP.txt'
  fAirL = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-DOWN.txt'
  sliceLoc = '/home/u/ugo/kvishal/1.PostProc/1.Airfoil/3.PviewScripts/1.0012-Slices/'
elif('4412' in airfoil):
  fAirU = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-UP.txt'
  fAirL = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-DOWN.txt'
  sliceLoc = '/home/u/ugo/kvishal/1.PostProc/1.Airfoil/3.PviewScripts/2.4412-Slices/'
else:
  raise ValueError('--|| ALYA ERROR :: FILE NOT PROVIDED.')


coordAirU = np.loadtxt(fAirU, delimiter=',')
coordAirL = np.loadtxt(fAirL, delimiter=',')

FyU = interp1d(coordAirU[:,0],coordAirU[:,1])
FyL = interp1d(coordAirL[:,0],coordAirL[:,1])

xuLin = np.linspace(np.amin(coordAirU[:,0],axis=None),np.amax(coordAirU[:,0],axis=None),1000);
yuLin =  FyU(xuLin);
xlLin = np.flip(np.linspace(np.amin(coordAirL[:,0],axis=None),np.amax(coordAirL[:,0],axis=None),1000));
ylLin =  FyL(xlLin);

coordAirU = np.column_stack((xuLin,yuLin));
coordAirL = np.column_stack((xlLin,ylLin));

# X=1 POINT 
coordF = [0.0, 0.0]
thAirF = np.arctan2(1.0,0.0)

# INTERPOLATION ON SUCTION SIDE
thAirU = np.arctan2(np.diff(coordAirU[:,1]),np.diff(coordAirU[:,0]))
F0 = interp1d(coordAirU[:,0],coordAirU[:,1])
airLen = len(coordAirU)
coordMid = 0.5*(coordAirU[0:airLen-1,0]+coordAirU[1:airLen,0])
Fth0 = interp1d(coordMid,thAirU)
coordU = 0.5*(coordAirU[0:airLen-1,:]+coordAirU[1:airLen,:]); 

# X=1 POINT 
coordC = [1.0, 0.0]
thAirC = 0.0

# INTERPOLATION ON PRESSURE SIDE
thAirL = np.arctan2(np.diff(coordAirL[:,1]),np.diff(coordAirL[:,0]))
F1 = interp1d(coordAirL[:,0],coordAirL[:,1])
airLen = len(coordAirL)
coordMid = 0.5*(coordAirL[0:airLen-1,0]+coordAirL[1:airLen,0])
Fth1 = interp1d(coordMid,thAirL)
coordL = 0.5*(coordAirL[0:airLen-1,:]+coordAirL[1:airLen,:])

coordInterp = np.vstack((coordL,coordC, coordU))
tInterp = np.hstack((thAirL,thAirC,thAirU));
indSort = np.lexsort((coordInterp[:,1],coordInterp[:,0]))
coordInterp = coordInterp[indSort,:];
tInterp = tInterp[indSort]

LthArray = np.zeros((Nx,3),dtype=np.double)
VthArray = np.zeros((Nx,3),dtype=np.double)
varArray = np.zeros((Nx,3),dtype=np.double)
for i in range(Nx):
  idx = find_nearest_index(coordInterp, (x[i], y[i]))
  th = tInterp[idx]
  #if(np.logical_and(x[i]<=1.0,y[i]<0.0)):
  #  LthArray[i,:] = (-cos(th),-sin(th),0.0);
  #else:
  LthArray[i,:] = (cos(th),sin(th),0.0);
  if(np.logical_and(x[i]>1.0,y[i]<0.0)):    
    LthArray[i,:] = (-cos(th),-sin(th),0.0);
    VthArray[i,:] = (sin(th),-cos(th),0.0);
  else:
    VthArray[i,:] = (-sin(th),cos(th),0.0);

LthArray = np.asarray(LthArray, dtype=np.float64)
VthArray = np.asarray(VthArray, dtype=np.float64)
output.ShallowCopy(inputs[0].VTKObject)
output.PointData.append(LthArray,"LAXES")
output.PointData.append(VthArray,"VAXES")

## CALCULATE VARIABLES
uArr = np.asarray(d.PointData["AVVEL"],dtype=np.double)
uuArr = np.asarray(d.PointData["AVVE2"],dtype=np.double)
uvArr = np.asarray(d.PointData["AVVXY"],dtype=np.double)

#-----AVVEL---------#
varArray = np.zeros((Nx,3),dtype=np.double)
varArray[:,0] = np.multiply(uArr[:,0],LthArray[:,0])
varArray[:,1] = np.multiply(uArr[:,1],LthArray[:,1])
varArray[:,2] = uArr[:,2]
var = np.asarray(varArray, dtype=np.float64)
output.PointData.append(var,"avvel")

#-----AVVE2---------#
varArray[:,0] = np.multiply(uuArr[:,0],np.square(LthArray[:,0])) + \
                np.multiply(uuArr[:,1],np.square(LthArray[:,1])) + \
                2.0*np.multiply(uvArr[:,1],np.multiply(LthArray[:,0],LthArray[:,1]))

varArray[:,1] = np.multiply(uuArr[:,0],np.square(VthArray[:,0])) + \
                np.multiply(uuArr[:,1],np.square(VthArray[:,1])) + \
                2.0*np.multiply(uvArr[:,1],np.multiply(VthArray[:,0],VthArray[:,1])) 
varArray[:,2] = uuArr[:,2]
var = np.asarray(varArray, dtype=np.float64)
output.PointData.append(var,"avve2")

#-----AVVXY---------#
varArray[:,0] = np.multiply(uuArr[:,0],np.multiply(LthArray[:,0],VthArray[:,0])) + \
                np.multiply(uuArr[:,1],np.multiply(LthArray[:,1],VthArray[:,1])) + \
                np.multiply(uvArr[:,1],np.multiply(LthArray[:,0],VthArray[:,1])) + \
                np.multiply(uvArr[:,1],np.multiply(LthArray[:,1],VthArray[:,0]))

varArray[:,1] = np.multiply(uvArr[:,2],VthArray[:,0]) + \
                np.multiply(uvArr[:,1],VthArray[:,1])

varArray[:,2] = np.multiply(uvArr[:,2],LthArray[:,0]) + \
                np.multiply(uvArr[:,1],LthArray[:,1])

var = np.asarray(varArray, dtype=np.float64)
output.PointData.append(var,"avvxy")

#-----AVPGR---------#
uArr = np.asarray(d.PointData["AVPGR"],dtype=np.double)
varArray[:,0] = np.multiply(uArr[:,0],LthArray[:,0])
varArray[:,1] = np.multiply(uArr[:,1],LthArray[:,1])
varArray[:,2] = uArr[:,2]
var = np.asarray(varArray, dtype=np.float64)
output.PointData.append(var,"avpgr")

#-----AVVGR---------#
uArr = np.asarray(d.PointData["AVVGR"],dtype=np.double)
uArr = np.reshape(uArr,(Nx,9))
uuArr = np.zeros((Nx,3),dtype=np.double)
uvArr = np.zeros((Nx,3),dtype=np.double)
varArray = np.zeros((Nx,9),dtype=np.double)

uuArr[:,0] = uArr[:,0]
uuArr[:,1] = uArr[:,4]
uuArr[:,2] = uArr[:,8]
#
uvArr[:,0] = uArr[:,1]
uvArr[:,1] = uArr[:,5]
uvArr[:,2] = uArr[:,2]
#
varArray[:,0] = np.multiply(uuArr[:,0],np.square(LthArray[:,0])) + \
                np.multiply(uuArr[:,1],np.square(LthArray[:,1])) + \
                2.0*np.multiply(uvArr[:,1],np.multiply(LthArray[:,0],LthArray[:,1]))

varArray[:,4] = np.multiply(uuArr[:,0],np.square(VthArray[:,0])) + \
                np.multiply(uuArr[:,1],np.square(VthArray[:,1])) + \
                2.0*np.multiply(uvArr[:,1],np.multiply(VthArray[:,0],VthArray[:,1])) 

varArray[:,8] = varArray[:,2]

varArray[:,1] = np.multiply(uuArr[:,0],np.multiply(LthArray[:,0],VthArray[:,0])) + \
                np.multiply(uuArr[:,1],np.multiply(LthArray[:,1],VthArray[:,1])) + \
                np.multiply(uvArr[:,1],np.multiply(LthArray[:,0],VthArray[:,1])) + \
                np.multiply(uvArr[:,1],np.multiply(LthArray[:,1],VthArray[:,0]))

varArray[:,3] = varArray[:,1]

varArray[:,5] = np.multiply(uvArr[:,2],VthArray[:,0]) + \
                np.multiply(uvArr[:,1],VthArray[:,1])

varArray[:,7] = varArray[:,5]

varArray[:,2] = np.multiply(uvArr[:,2],LthArray[:,0]) + \
                np.multiply(uvArr[:,1],LthArray[:,1])

varArray[:,6] = varArray[:,2]

var = np.asarray(np.reshape(varArray,(Nx,3,3)), dtype=np.float64)
output.PointData.append(var,"avvgr")
