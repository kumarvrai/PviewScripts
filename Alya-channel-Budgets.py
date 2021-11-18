#### import the simple module from the paraview
#activate_this = '/home/kvishal/.virtualenvs/AlyaPythonEnv/bin/'
#pyVer = sys.version_info[0]
#if pyVer < 3:
# execfile(activate_this, dict(__file__=activate_this))
#else:
# exec(open(activate_this).read())
#import mpi4py
#mpi4py.rc.recv_mprobe = False
#
#import pyAlya
 
import os
import time
import sys 
import numpy as np
import itertools
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

caseName = sys.argv[1]
model    = sys.argv[2]
nu       = float(sys.argv[3])
dim      = sys.argv[4]
code     = sys.argv[5]
method   = sys.argv[6]
form     = 'CON'

zDec = 6; xDec = 6;

casePath = os.getcwd()

if("ALYA" in code):
 fileName = caseName+'.ensi.case'
elif("CITY" in code):
 fileName = 'read_h5.xdmf'
elif("VTK" in code):
 fileName = caseName+'.pvd'
else:
 raise ValueError('--|| ALYA ERROR :: NO CODE SPECIFIED.')

print "--|| ALYA :: READING ARRAYS"
startTime = time.time()

if("ALYA" in code):
 case1 = OpenDataFile(fileName)
 if(model == "LES"):
   case1.PointArrays = ['VELOC', 'PRESS', 'TURBU']
 else:
   case1.PointArrays = ['VELOC', 'PRESS']
elif("CITY" in code):
 case1 = XDMFReader(registrationName=fileName, FileNames= casePath+'/'+fileName)
 if(model == "LES"):
   case1.CellArrayStatus = ['VELOC', 'PRESS', 'TURBU']
 else:
   case1.CellArrayStatus = ['VELOC', 'PRESS']
   #case1.UpdatePipeline()
   #case1 = CellDatatoPointData(Input=case1)
elif("VTK" in code):
 case1 = OpenDataFile(fileName)
 if(model == "LES"):
   case1.PointArrays = ['VELOC', 'PRESS', 'TURBU']
 else:
   case1.PointArrays = ['VELOC', 'PRESS']
else:
 raise ValueError('--|| ALYA ERROR :: NO CODE SPECIFIED.')

case1.UpdatePipeline()
print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

#print "--|| ALYA :: CLIPPING DATA"
#startTime = time.time()
#case1 = Clip(Input=case1)
#case1.ClipType.Origin = [0.0, 0.0, 0.1]
#case1.ClipType.Normal = [0.0, 1.0, 0.0]
#case1.Invert=0
#case1.UpdatePipeline()
#case1 = Clip(Input=case1)
#case1.ClipType.Origin = [1.1, 0.0, 0.1]
#case1.ClipType.Normal = [1.0, 0.0, 0.0]
#case1.UpdatePipeline()
#case1 = Clip(Input=case1)
#case1.ClipType.Origin = [-0.1, 0.0, 0.1]
#case1.ClipType.Normal = [1.0, 0.0, 0.0]
#case1.Invert=0
#case1.UpdatePipeline()
#case1 = Clip(Input=case1)
#case1.ClipType.Origin = [0.0, 0.5, 0.1]
#case1.ClipType.Normal = [0.0, 1.0, 0.0]
#case1.UpdatePipeline()
#print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

Ntotal = case1.GetDataInformation().GetNumberOfPoints()
print("----|| ALYA :: WORKING WITH ",Ntotal," TOTAL NUMBER OF POINTS")


print "--|| ALYA :: PERFORMING TEMPORAL AVERAGING"
startTime = time.time()
case2 = TemporalStatistics(Input=case1)

# Properties modified on temporalStatistics1
case2.ComputeMinimum = 0
case2.ComputeMaximum = 0
case2.ComputeStandardDeviation = 0

case2.UpdatePipeline()
print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: APPEND DATA ATTRIBUTES"
startTime = time.time()

APND1 = AppendAttributes(Input=[case1,case2])
APND1.UpdatePipeline()
 
print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE VELOCITY FLUCTUATIONS"
startTime = time.time()

CAL1 = Calculator(Input=APND1)
CAL1.ResultArrayName = "u"
CAL1.Function = "VELOC - VELOC_average"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE VRMS"
startTime = time.time()
CAL1 = Calculator(Input=CAL1)
CAL1.ResultArrayName = "VRMS"
CAL1.Function = "u_X^2*iHat+u_Y^2*jHat+u_Z^2*kHat"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE RSSIJ"
startTime = time.time()
CAL1 = Calculator(Input=CAL1)
CAL1.ResultArrayName = "RSSIJ"
CAL1.Function = "u_X*u_Y*iHat+u_Y*u_Z*jHat+u_Z*u_X*kHat"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE PRESSURE FLUCTUATIONS"
startTime = time.time()

CAL1 = Calculator(Input=CAL1)
CAL1.ResultArrayName = "p"
CAL1.Function = "PRESS - PRESS_average"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE P-RMS"
startTime = time.time()
CAL1 = Calculator(Input=CAL1)
CAL1.ResultArrayName = "PRMS"
CAL1.Function = "p^2"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: APPEND CASE SPECIFIC VARIABLES"
startTime = time.time()
PF1 = ProgrammableFilter(Input=CAL1)

PF1.Script = \
"""
import numpy as np
import itertools

inp1 = inputs[0].PointData["VELOC_average"]
output.PointData.append(inp1,"AVVEL")

inp1 = inputs[0].PointData["u"]
output.PointData.append(inp1,"u")

inp1 = inputs[0].PointData["PRESS_average"]
output.PointData.append(inp1,"AVPRE")

inp1 = inputs[0].PointData["p"]
output.PointData.append(inp1,"p")

inp1 = inputs[0].PointData["PRMS"]
output.PointData.append(inp1,"PRMS")

inp1 = inputs[0].PointData["VRMS"]
output.PointData.append(inp1,"VRMS")

inp1 = inputs[0].PointData["RSSIJ"]
output.PointData.append(inp1,"RSSIJ")

if(model == "LES"):
  inp1 = inputs[0].PointData["TURBU_average"]
  output.PointData.append(inp1,"AVTUR")

"""

PF1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE INSTANTANEOUS BUDGET TERMS"
startTime = time.time()

if(model == "LES"):
  print "----|| ALYA :: CALCULATE FLUCT GRADIENT"
  sTime = time.time()
  # CALCULATE FLUCTUATIONS GRADIENT
  PF1 = GradientOfUnstructuredDataSet(Input=PF1)
  PF1.ScalarArray = ['POINTS', 'u']
  PF1.ResultArrayName = 'grad_u'
  PF1.FasterApproximation = 1
  PF1.UpdatePipeline()
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-sTime), 'sec'

  print "----|| ALYA :: CALCULATE S_IJ LES"
  sTime = time.time()
  # CALCULATE S_IJ COMPONENTS
  count = 1
  for k in range(0,3):
    for i in range(0,3):
      PF1 = Calculator(Input=PF1)
      PF1.ResultArrayName = "S_%s%s" % (k,i)
      PF1.Function = "-AVTUR*(grad_u_%s + grad_u_%s + grad_u_%s + grad_u_%s + grad_u_%s + grad_u_%s)" % (3*k,3*k+1,3*k+2,i,i+3,i+6)
      PF1.UpdatePipeline()    
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-sTime), 'sec'
  
  print "----|| ALYA :: CALCULATE S_IJ VECTOR"
  sTime = time.time()
  # MAKE S_IJ VECTORS
  for k in range(0,3):
    PF1 = Calculator(Input=PF1)
    PF1.ResultArrayName = "S_%s" % (k)
    PF1.Function = "S_%s%s*iHat + S_%s%s*jHat + S_%s%s*kHat" % (0,k,1,k,2,k)
    PF1.UpdatePipeline()    
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-sTime), 'sec'

#CALCULATE UU, VV, WW  BUDGET TERMS
count=1
if(model == "LES"):
  varNames = ["C","P","TD","PV","VDD","EPS","VD","SGSDD","SGSD","SGSEps"]
else:
  varNames = ["C","P","TD","PV","VDD","EPS","VD"]

varAddOn = ["uu","vv","ww"]

for var in varNames:
   for i in range(0,3):
      print "----|| ALYA :: CALCULATE",var+varAddOn[i], "TERM"
      sTime = time.time()

      PF1 = PythonCalculator(Input=PF1)

      if(form == "CON"):

        if var == "C":
           PF1.Expression = '2*u[:,%d]*(divergence(AVVEL*u[:,%d]))' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "P":
           PF1.Expression = '-2.0*u[:,%d]*(divergence(AVVEL[:,%d]*u))' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "TD":
           PF1.Expression = '-2*u[:,%d]*(divergence(u[:,%d]*u))' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "PV":
           PF1.Expression = '-2*u[:,%d]*(gradient(p)[:,%d])' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "VDD":
           PF1.Expression = '2*%f*u[:,%d]*(laplacian(u[:,%d]))' % (nu,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "EPS":
           PF1.Expression = '-2*%f*dot(gradient(u[:,%d]),gradient(u[:,%d]))' % (nu,i,i)
           #PF1.Expression = 'VDD_%s - VD_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "VD":
           #PF1.Expression = '%f*laplacian(u[:,%d]*u[:,%d])' % (nu,i,i)
           PF1.Expression = 'VDD_%s - EPS_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        #elif var == "VD":
        #   PF1.Expression = '%f*(divergence(gradient(u[:,%d]*u[:,%d])))' % (nu,i,i)
        #   PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        #elif var == "EPS":
        #   PF1.Expression = 'VDD_%s - VD_%s' % (varAddOn[i],varAddOn[i])
        #   PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "SGSDD":
           PF1.Expression = '-2.0*u[:,%d]*divergence(S_%s)' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "SGSD":
           PF1.Expression = '-2.0*divergence(u[:,%d]*S_%s)' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "SGSEps":
           PF1.Expression = 'SGSDD_%s - SGSD_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif(form == "NCON"):

        if var == "C":
           PF1.Expression = '2.0*u[:,%d]*(dot(AVVEL,gradient(u[:,%d]))+u[:,%d]*divergence(AVVEL))' % (i,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "P":
           PF1.Expression = '-2.0*u[:,%d]*(dot(u,gradient(AVVEL[:,%d]))+AVVEL[:,%d]*divergence(u))' % (i,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "TD":
           PF1.Expression = '-2*(u[:,%d]*dot(u,gradient(u[:,%d]))+u[:,%d]*u[:,%d]*divergence(u))' % (i,i,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "PV":
           PF1.Expression = '-2*u[:,%d]*(gradient(p)[:,%d])' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "VDD":
           PF1.Expression = '2*%f*u[:,%d]*laplacian(u[:,%d])' % (nu,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])


        elif var == "EPS":
           PF1.Expression = '-2*%f*dot(gradient(u[:,%d]),gradient(u[:,%d]))' % (nu,i,i)
           #PF1.Expression = 'VDD_%s - VD_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "VD":
           #PF1.Expression = '%f*laplacian(u[:,%d]*u[:,%d])' % (nu,i,i)
           PF1.Expression = 'VDD_%s - EPS_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])


        elif var == "SGSDD":
           PF1.Expression = '-2.0*u[:,%d]*divergence(S_%s)' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "SGSD":
           PF1.Expression = '-2.0*divergence(u[:,%d]*S_%s)' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "SGSEps":
           PF1.Expression = 'SGSDD_%s - SGSD_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif(form == "MXED"):

        if var == "C":
           PF1.Expression = '2*u[:,%d]*(divergence(AVVEL*u[:,%d]))' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "P":
           PF1.Expression = '-2.0*u[:,%d]*(dot(u,gradient(AVVEL[:,%d]))+AVVEL[:,%d]*divergence(u))' % (i,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "TD":
           PF1.Expression = '-2*u[:,%d]*(divergence(u[:,%d]*u))' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "PV":
           PF1.Expression = '-2*u[:,%d]*(gradient(p)[:,%d])' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "VDD":
           PF1.Expression = '2*%f*u[:,%d]*laplacian(u[:,%d])' % (nu,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])


        elif var == "EPS":
           PF1.Expression = '-2*%f*dot(gradient(u[:,%d]),gradient(u[:,%d]))' % (nu,i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "VD":
           PF1.Expression = 'VDD_%s - EPS_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])


        elif var == "SGSDD":
           PF1.Expression = '-2.0*u[:,%d]*divergence(S_%s)' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "SGSD":
           PF1.Expression = '-2.0*divergence(u[:,%d]*S_%s)' % (i,i)
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])

        elif var == "SGSEps":
           PF1.Expression = 'SGSDD_%s - SGSD_%s' % (varAddOn[i],varAddOn[i])
           PF1.ArrayName = '%s_%s' % (var,varAddOn[i])
      else:
        raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE THE TERMS?')

      PF1.UpdatePipeline()
      print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-sTime), 'sec'
print "--|| ALYA :: TERMS DONE. TIME =",time.time()-startTime,'sec'

#
# PERFORM TEMPORAL AVERAGING OF BUDGET TERMS
print "--|| ALYA :: TEMPORAL AVERAGING OF BUDGET TERMS"
startTime = time.time()
case = TemporalStatistics(Input=PF1)
#
# Properties modified on temporalStatistics1
case.ComputeMinimum = 0
case.ComputeMaximum = 0
case.ComputeStandardDeviation = 0
case.UpdatePipeline()
print "--|| ALYA :: AVERAGING DONE. TIME =",time.time()-startTime,'sec'

if('3D' in dim):
 # Save a 3D time averaged file
 savePath = casePath+"/BudgetData_3D.vtm"
 SaveData(savePath, proxy=case)
 print "----|| ALYA: 3D BUDGET FILE WRITTEN "

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

Nplane = slice1.GetDataInformation().GetNumberOfPoints()
print("----|| ALYA :: WORKING WITH ",Nplane," PLANAR POINTS")

N = int(Ntotal/Nplane)
zmid = (zmin+zmax)/2
zpos = np.around(np.asarray(np.arange(N)*(zmax-zmin)/(N-1),dtype=np.double),decimals=zDec)
print("----|| ALYA: WORKING WITH %d Z-PLANES" % (len(zpos)))
print("----|| ALYA: DELTA-Z = %f" % ((zmax-zmin)/(N-1)))
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

if("PINTERP" in method):
  print("--|| ALYA :: CREATING TRANSFORMATIONS")
  startTime = time.time()
  resample_transforms=list();
  data=list();
  for i in range(N):
    # create a new 'Transform'
    transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
    resampleWithDataset1 = ResampleWithDataset(Input=case,Source=transform1,guiName="resample{}".format(i))
    resample_transforms.append(resampleWithDataset1)
    # Properties modified on transform1.Transform
    transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
    #data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])
  print("--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA: AVERAGING USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  ## create a new 'Programmable Filter'
  PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
  ### the rest of them are data to be averaged
  PF1.Script = \
  """
  import numpy as np

  varFull = inputs[0].PointData.keys()
  N=len(inputs);
  print "----|| ALYA WORKING ON ARRAYS ",varNames
  # Boudget Variable
  for varName in varFull:
   print("--|| ALYA :: CALCULATING FOR",varName)
   varName0 = varName.replace('_average','')
   avg = inputs[0].PointData[varName]
   for i in range(1,N):
       d = inputs[i].PointData[varName]
       avg = avg + d
   avg = avg/N
   output.PointData.append(avg,varName0)
  
  """
  PF1.UpdatePipeline()
  print "--|| ALYA: SPANWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec'

else:
  print("--|| ALYA :: SPANWISE AVERAGING")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=[case,slice1])
  ### first input is the grid
  ### the rest of them are data to be averaged
  PF1.Script = \
"""
import numpy as np
t = inputs[0].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
varFull = inputs[0].PointData.keys()
#print("--|| ALYA :: CALCULATING FOR",varFull," AT T=",t) 

#-----------------3D GEOMETRY--------------------#
if("ALYA" in code):
 d = dsa.WrapDataObject(inputs[0].GetBlock(0))
elif("CITY" in code):
 d = inputs[0]
else:
 raise ValueError('--|| ALYA ERROR :: NO CODE SPECIFIED.')
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=zDec)

#-----------------SLICE GEOMETRY--------------------#
if("ALYA" in code):
 out = dsa.WrapDataObject(inputs[1].GetBlock(0))
elif("CITY" in code):
 out = inputs[1]
else:
 raise ValueError('--|| ALYA ERROR :: NO CODE SPECIFIED.')

x_2d = np.around(np.asarray(out.Points[:,0],dtype=np.double),decimals=6)
y_2d = np.around(np.asarray(out.Points[:,1],dtype=np.double),decimals=6)
ind_2d = np.lexsort((y_2d,x_2d));
print("--|| ALYA CHECK :: SHAPE OF SLICE", np.shape(ind_2d))

#-----------------CALCULATE AVERAGE--------------------#
for varName in varFull:
 print("--|| ALYA :: CALCULATING FOR",varName)
 varName0 = varName.replace('_average','')
 avg = 0.0*np.asarray(out.PointData[varName],dtype=np.double)
 # s_src = 0.0*np.asarray(out.PointData[varName],dtype=np.double)
 src = np.asarray(d.PointData[varName],dtype=np.double)
 for n in range(N):
  ind_z = np.where(z == zpos[n])
  # print("--|| ALYA CHECK :: AT n,Z=",n,zpos[n],"SHAPE OF MASTER DATASET", np.shape(ind_z))
  x_l = x[ind_z];
  y_l = y[ind_z];
  s_src = src[ind_z]
  if("SORT" in method):
   ind = np.lexsort((y_l,x_l)); 
   avg[ind_2d] = avg[ind_2d] + s_src[ind]
  elif("INTERP" in method):
   avg = avg + griddata((x_l,y_l),s_src, (x_2d,y_2d), method='nearest')
  else:
   raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')
 avg = avg/N;
 avg = np.asarray(avg, dtype=np.float64)
 output.ShallowCopy(inputs[1].VTKObject)
 output.PointData.append(avg, varName0)

"""
  PF1.UpdatePipeline() 
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
  if not "_average" in varName:
   inp0 = inputs[0].PointData[varName]
   output.PointData.append(inp0,varName)
"""
PF1.UpdatePipeline() 
PF1 = ProgrammableFilter(Input=PF1)
PF1.Script = \
"""
import numpy as np
varFull = inputs[0].PointData.keys()
for varName in varFull:
 if varName not in ["S_0","S_1","S_2","S_00","S_01","S_02","S_10","S_11","S_12","S_20","S_21","S_22"]:
   inp0 = inputs[0].PointData[varName]
   output.PointData.append(inp0,varName)
"""
PF1.UpdatePipeline() 
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
## APPEND SPECIFIC VARIABLES ###
print("--|| ALYA :: UNAPPEND VARIABLES")
startTime = time.time()
PF1 = ProgrammableFilter(Input=PF1)
PF1.Script = \
"""
import numpy as np
varFull = inputs[0].PointData.keys()
for varName in varFull:
 if varName not in ["u","p","grad_u"]:
   inp0 = inputs[0].PointData[varName]
   output.PointData.append(inp0,varName)
"""
PF1.UpdatePipeline() 
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

# GRADIENT CALC
print("--|| ALYA :: CALCULATING PRESS GRADIENT")
startTime = time.time()
PF1 = GradientOfUnstructuredDataSet(Input=PF1)
PF1.FasterApproximation = 1
PF1.ScalarArray = ['POINTS', 'AVPRE']
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
PF1.FasterApproximation = 1
PF1.ScalarArray = ['POINTS', 'AVVEL']
PF1.ComputeGradient = 1
PF1.ResultArrayName = 'AVVGR'
PF1.ComputeVorticity = 0
PF1.VorticityArrayName = 'OMEGA'
PF1.ComputeQCriterion = 0
PF1.QCriterionArrayName = 'QCRIT'
PF1.UpdatePipeline()
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
#
#### write Files
if('2D' in dim):
  print "--|| ALYA :: CALCULATE VRMS"
  startTime = time.time()
  PF1 = Calculator(Input=PF1)
  PF1.ResultArrayName = "VRMS"
  PF1.Function = "sqrt(VRMS_X)*iHat+sqrt(VRMS_Y)*jHat+sqrt(VRMS_Z)*kHat"
  PF1.UpdatePipeline()
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

  print "--|| ALYA :: CALCULATE PRMS"
  startTime = time.time()
  PF1 = Calculator(Input=PF1)
  PF1.ResultArrayName = "PRMS"
  PF1.Function = "sqrt(PRMS)"
  PF1.UpdatePipeline()

  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'
  savePath = casePath+"/BudgetData_2D.vtm"
  SaveData(savePath, proxy=PF1)
  savePath = casePath+"/BudgetData_2D.csv"
  SaveData(savePath, proxy=PF1)
  print "--|| ALYA :: 2D BUDGETS FILE WRITTEN "

########### STREAMWISE AVERAGING ################
## read Coordinates File
if('1D' in dim):
  (xmin,xmax,ymin,ymax,zmin,zmax) =  PF1.GetDataInformation().GetBounds()
  print("--|| ALYA GEOM:: XMIN=",xmin,"XMAX=",xmax,"YMIN=",ymin,"YMAX=",ymax,"ZMIN=",zmin,"ZMAX=",zmax)

  slice1 = Slice(Input=PF1)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  ## init the 'Plane' selected for 'SliceType'
  slice1.SliceType.Origin = [xmin, (ymin+ymax)/2, (zmin+zmax)/2]
  ## Properties modified on slice1.SliceType
  slice1.SliceType.Normal = [-1.0, 0.0, 0.0]
  slice1.UpdatePipeline()
  Ny = slice1.GetDataInformation().GetNumberOfPoints()
  print("----|| ALYA :: WORKING WITH ",Ny," PLANAR POINTS")

  N = int(Nplane/Ny)
  xmid = (xmin+xmax)/2
  xpos = np.around(np.asarray(np.arange(N)*(xmax-xmin)/(N-1),dtype=np.double),decimals=xDec)
  print("----|| ALYA: WORKING WITH %d X-PLANES" % (len(xpos)))
  print("----|| ALYA: DELTA-X = %f" % ((xmax-xmin)/(N-1)))

  print("--|| ALYA :: GENERATING SLICE FOR STREAMWISE AVERAGE")
  slice1 = Slice(Input=PF1)
  slice1.SliceType = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  ## init the 'Plane' selected for 'SliceType'
  slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
  ## Properties modified on slice1.SliceType
  slice1.SliceType.Normal = [1.0, 0.0, 0.0]
  slice1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  
  print("--|| ALYA :: CREATING TRANSFORMATIONS")
  startTime = time.time()
  resample_transforms=list();
  data=list();
  for i in range(N):
    # create a new 'Transform'
    transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
    resampleWithDataset1 = ResampleWithDataset(Input=PF1,Source=transform1,guiName="resample{}".format(i))
    resample_transforms.append(resampleWithDataset1)
    # Properties modified on transform1.Transform
    transform1.Transform.Translate = [xpos[i]-xmid, 0.0, 0.0]
    #data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])
  print("--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA: AVERAGING USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  ## create a new 'Programmable Filter'
  PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
  ### the rest of them are data to be averaged
  PF1.Script = \
  """
  import numpy as np

  varFull = inputs[0].PointData.keys()
  N=len(inputs);
  print "----|| ALYA WORKING ON ARRAYS ",varNames
  # Boudget Variable
  for varName in varFull:
   print("--|| ALYA :: CALCULATING FOR",varName)
   avg = inputs[0].PointData[varName]
   for i in range(1,N):
       #if("AVVEL" in varName):
       # print(i,np.amax(avg,axis=None))
       d = inputs[i].PointData[varName]
       avg = avg + d
   avg = avg/N
   output.PointData.append(avg,varName)
  
  """
  PF1.UpdatePipeline()
  print "--|| ALYA: DONE. TIME =",time.time()-startTime,'sec'

  print "--|| ALYA :: CALCULATE VRMS"
  startTime = time.time()
  PF1 = Calculator(Input=PF1)
  PF1.ResultArrayName = "VRMS"
  PF1.Function = "sqrt(VRMS_X)*iHat+sqrt(VRMS_Y)*jHat+sqrt(VRMS_Z)*kHat"
  PF1.UpdatePipeline()
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

  print "--|| ALYA :: CALCULATE PRMS"
  startTime = time.time()
  PF1 = Calculator(Input=PF1)
  PF1.ResultArrayName = "PRMS"
  PF1.Function = "sqrt(PRMS)"
  PF1.UpdatePipeline()

  savePath = casePath+"/BudgetData_1D.csv"
  SaveData(savePath, proxy=PF1)
  print "--|| ALYA: 1D BUDGETS FILE WRITTEN"
