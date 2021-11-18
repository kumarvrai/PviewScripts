#### import the simple module from the paraview
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
method     = sys.argv[6]

zDec = 4; xDec = 4

casePath = os.getcwd()

if("ALYA" in code):
 fileName = caseName+'.ensi.case'
elif("CITY" in code):
 fileName = 'read_h5.xdmf'
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
else:
 raise ValueError('--|| ALYA ERROR :: NO CODE SPECIFIED.')

case1.UpdatePipeline()
print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

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
CAL1.AttributeType = 'Cell Data'
CAL1.ResultArrayName = "u"
CAL1.Function = "VELOC - VELOC_average"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE VRMS"
startTime = time.time()
CAL1 = Calculator(Input=CAL1)
CAL1.AttributeType = 'Cell Data'
CAL1.ResultArrayName = "VRMS"
CAL1.Function = "u_X^2*iHat+u_Y^2*jHat+u_Z^2*kHat"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE RSSIJ"
startTime = time.time()
CAL1 = Calculator(Input=CAL1)
CAL1.AttributeType = 'Cell Data'
CAL1.ResultArrayName = "RSSIJ"
CAL1.Function = "u_X*u_Y*iHat+u_Y*u_Z*jHat+u_Z*u_X*kHat"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE PRESSURE FLUCTUATIONS"
startTime = time.time()

CAL1 = Calculator(Input=CAL1)
CAL1.AttributeType = 'Cell Data'
CAL1.ResultArrayName = "p"
CAL1.Function = "PRESS - PRESS_average"
CAL1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE P-RMS"
startTime = time.time()
CAL1 = Calculator(Input=CAL1)
CAL1.AttributeType = 'Cell Data'
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

inp1 = inputs[0].CellData["VELOC_average"]
output.CellData.append(inp1,"AVVEL")

inp1 = inputs[0].CellData["u"]
output.CellData.append(inp1,"u")

inp1 = inputs[0].CellData["PRESS_average"]
output.CellData.append(inp1,"AVPRE")

inp1 = inputs[0].CellData["p"]
output.CellData.append(inp1,"p")

inp1 = inputs[0].CellData["PRMS"]
output.CellData.append(inp1,"PRMS")

inp1 = inputs[0].CellData["VRMS"]
output.CellData.append(inp1,"VRMS")

inp1 = inputs[0].CellData["RSSIJ"]
output.CellData.append(inp1,"RSSIJ")

if(model == "LES"):
  inp1 = inputs[0].CellData["TURBU_average"]
  output.CellData.append(inp1,"AVTUR")

"""

PF1.UpdatePipeline()

print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

print "--|| ALYA :: CALCULATE INSTANTANEOUS BUDGET TERMS"
startTime = time.time()

if(model == "LES"):

  # CALCULATE FLUCTUATIONS GRADIENT
  GOUD1 = GradientOfUnstructuredDataSet(Input=PF1)
  GOUD1.ScalarArray = ['POINTS', 'u']
  GOUD1.ResultArrayName = 'grad_u'

  # CALCULATE S_IJ COMPONENTS
  count = 1
  for k in range(0,3):
    for i in range(0,3):
      if count==1: 
         pc = Calculator(Input=GOUD1)
         pc.AttributeType = 'Cell Data'
      else:
         strn = 'Calculator%d' % (count-1)
         pc_old = FindSource(strn) 
         pc = Calculator(Input=pc_old)
         pc.AttributeType = 'Cell Data'
      count = count+1
      pc.ResultArrayName = "S_%s%s" % (k,i)
      pc.Function = "-AVTUR*(grad_u_%s + grad_u_%s + grad_u_%s + grad_u_%s + grad_u_%s + grad_u_%s)" % (3*k,3*k+1,3*k+2,i,i+3,i+6)
  
  # MAKE S_IJ VECTORS
  for k in range(0,3):
    strn = 'Calculator%d' % (count-1)
    pc_old = FindSource(strn)
    pc = Calculator(Input=pc_old)
    pc.AttributeType = 'Cell Data'
    count = count+1
    pc.ResultArrayName = "S_%s" % (k)
    pc.Function = "S_%s%s*iHat + S_%s%s*jHat + S_%s%s*kHat" % (0,k,1,k,2,k)

#CALCULATE UU, VV, WW  BUDGET TERMS

count=1
if(model == "LES"):
  varNames = ["C","P","TD","PV","VDD","VD","EPS", "SGSDD", "SGSD", "SGSEps"]
else:
  varNames = ["C","P","TD","PV","VDD","VD","EPS"]

varAddOn = ["uu","vv","ww"]

for var in varNames:
   for i in range(0,3):
      if count==1: 
         if(model == "LES"):
           pc = PythonCalculator(Input=pc)
           pc.ArrayAssociation = 'Cell Data'
	 else:
	   pc = PythonCalculator(Input=PF1)
           pc.ArrayAssociation = 'Cell Data'
      else:
         strn = 'PythonCalculator%d' % (count-1)
         pc_old = FindSource(strn) 
         pc = PythonCalculator(Input=pc_old)
         pc.ArrayAssociation = 'Cell Data'
      count = count+1
      if var == "C":
         pc.Expression = '2*u[:,%d]*(divergence(AVVEL*u[:,%d]))' % (i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "P":
         pc.Expression = '-2.0*u[:,%d]*(divergence(AVVEL[:,%d]*u))' % (i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "TD":
         pc.Expression = '-2*u[:,%d]*(divergence(u[:,%d]*u))' % (i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "PV":
         pc.Expression = '-2*u[:,%d]*(gradient(p)[:,%d])' % (i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "VDD":
         pc.Expression = '2*%f*u[:,%d]*(divergence(gradient(u[:,%d])))' % (nu,i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "VD":
         pc.Expression = '%f*(divergence(gradient(u[:,%d]*u[:,%d])))' % (nu,i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "EPS":
         pc.Expression = 'VDD_%s - VD_%s' % (varAddOn[i],varAddOn[i])
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "SGSDD":
         pc.Expression = '-2.0*u[:,%d]*divergence(S_%s)' % (i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "SGSD":
         pc.Expression = '-2.0*divergence(u[:,%d]*S_%s)' % (i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "SGSEps":
         pc.Expression = 'SGSDD_%s - SGSD_%s' % (varAddOn[i],varAddOn[i])
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])


strn = 'PythonCalculator%d' % (len(varNames)*len(varAddOn))
case = FindSource(strn) 
#
case.UpdatePipeline()
print "--|| ALYA :: TERMS DONE. TIME =",time.time()-startTime,'sec'

#
# PERFORM TEMPORAL AVERAGING OF BUDGET TERMS
print "--|| ALYA :: TEMPORAL AVERAGING OF BUDGET TERMS"
startTime = time.time()
case = TemporalStatistics(Input=case)
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
    resampleWithDataset1.PassCellArrays = 1
    resample_transforms.append(resampleWithDataset1)
    # Properties modified on transform1.Transform
    transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
    #data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).CellData['AVVEL'])
  print("--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA: AVERAGING USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  ## create a new 'Programmable Filter'
  PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
  ### the rest of them are data to be averaged
  PF1.Script = \
  """
  import numpy as np

  varFull = inputs[0].CellData.keys()
  N=len(inputs);
  print "----|| ALYA WORKING ON ARRAYS ",varNames
  # Boudget Variable
  for varName in varFull:
   print("--|| ALYA :: CALCULATING FOR",varName)
   varName0 = varName.replace('_average','')
   avg = inputs[0].CellData[varName]
   for i in range(1,N):
       d = inputs[i].CellData[varName]
       avg = avg + d
   avg = avg/N
   output.CellData.append(avg,varName0)
  
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
varFull = inputs[0].CellData.keys()
print("--|| ALYA :: CALCULATING FOR",varFull," AT T=",t) 

#-----------------3D GEOMETRY--------------------#
if("ALYA" in code):
 d = dsa.WrapDataObject(inputs[0].GetBlock(0))
elif("CITY" in code):
 d = dsa.WrapDataObject(inputs[0])
else:
 raise ValueError('--|| ALYA ERROR :: NO CODE SPECIFIED.')
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=zDec)

#-----------------SLICE GEOMETRY--------------------#
if("ALYA" in code):
 out = dsa.WrapDataObject(inputs[1].GetBlock(0))
elif("CITY" in code):
 out = dsa.WrapDataObject(inputs[1].GetBlock(0))
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
 avg = 0.0*np.asarray(out.CellData[varName],dtype=np.double)
 # s_src = 0.0*np.asarray(out.CellData[varName],dtype=np.double)
 src = np.asarray(d.CellData[varName],dtype=np.double)
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
 output.CellData.append(avg, varName0)

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
varFull = inputs[0].CellData.keys()
for varName in varFull:
 if varName not in ["u","p"]:
   inp0 = inputs[0].CellData[varName]
   output.CellData.append(inp0,varName)
"""
PF1.UpdatePipeline() 
print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

# GRADIENT CALC
print("--|| ALYA :: CALCULATING PRESS GRADIENT")
startTime = time.time()
PF1 = GradientOfUnstructuredDataSet(Input=PF1)
PF1.ScalarArray = 'AVPRE'
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
PF1.ScalarArray = 'AVVEL'
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
  PF1.AttributeType = 'Cell Data'
  PF1.ResultArrayName = "VRMS"
  PF1.Function = "sqrt(VRMS_X)*iHat+sqrt(VRMS_Y)*jHat+sqrt(VRMS_Z)*kHat"
  PF1.UpdatePipeline()
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

  print "--|| ALYA :: CALCULATE PRMS"
  startTime = time.time()
  PF1 = Calculator(Input=PF1)
  PF1.AttributeType = 'Cell Data'
  PF1.ResultArrayName = "PRMS"
  PF1.Function = "sqrt(PRMS)"
  PF1.UpdatePipeline()

  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

  PF1=CellDatatoPointData(Input=PF1)
  PF1.UpdatePipeline()

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
  print("--|| ALYA :: GENERATING SLICE FOR STREAMWISE AVERAGE")
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
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA :: CREATING TRANSFORMATIONS")
  startTime = time.time()
  resample_transforms=list();
  data=list();
  for i in range(N):
    # create a new 'Transform'
    transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
    resampleWithDataset1 = ResampleWithDataset(Input=PF1,Source=transform1,guiName="resample{}".format(i))
    resampleWithDataset1.PassCellArrays = 1
    resample_transforms.append(resampleWithDataset1)
    # Properties modified on transform1.Transform
    transform1.Transform.Translate = [xpos[i]-xmid, 0.0, 0.0]
    #data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).CellData['AVVEL'])
  print("--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA: AVERAGING USING A PROGRAMMABLE FILTER")
  startTime = time.time()
  ## create a new 'Programmable Filter'
  PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
  ### the rest of them are data to be averaged
  PF1.Script = \
  """
  import numpy as np

  varFull = inputs[0].CellData.keys()
  N=len(inputs);
  print "----|| ALYA WORKING ON ARRAYS ",varNames
  # Boudget Variable
  for varName in varFull:
   print("--|| ALYA :: CALCULATING FOR",varName)
   avg = inputs[0].CellData[varName]
   for i in range(1,N):
       d = inputs[i].CellData[varName]
       avg = avg + d
   avg = avg/N
   output.CellData.append(avg,varName)
  
  """
  PF1.UpdatePipeline()
  print "--|| ALYA: SPANWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec'

  print "--|| ALYA :: CALCULATE VRMS"
  startTime = time.time()
  PF1 = Calculator(Input=PF1)
  PF1.AttributeType = 'Cell Data'
  PF1.ResultArrayName = "VRMS"
  PF1.Function = "sqrt(VRMS_X)*iHat+sqrt(VRMS_Y)*jHat+sqrt(VRMS_Z)*kHat"
  PF1.UpdatePipeline()
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

  print "--|| ALYA :: CALCULATE PRMS"
  startTime = time.time()
  PF1 = Calculator(Input=PF1)
  PF1.AttributeType = 'Cell Data'
  PF1.ResultArrayName = "PRMS"
  PF1.Function = "sqrt(PRMS)"
  PF1.UpdatePipeline()
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

  print "--|| ALYA :: CELL TO POINT DATA"
  startTime = time.time()
  PF1=CellDatatoPointData(Input=PF1)
  PF1.UpdatePipeline()
  print "----|| ALYA :: DONE. TIME TAKEN =", (time.time()-startTime), 'sec'

  savePath = casePath+"/BudgetData_1D_cell.csv"
  SaveData(savePath, proxy=PF1,FieldAssociation='Points')
  print "--|| ALYA: 1D BUDGETS FILE WRITTEN"
