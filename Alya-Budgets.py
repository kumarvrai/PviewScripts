#### import the simple module from the paraview
import os
import time
import sys 
import numpy as np
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

caseName = sys.argv[1]
nu = float(sys.argv[2])

geofile_in= os.path.join('./'+caseName+'.geo.dat')
bString = 'COORDINATES'
eString = 'END_COORDINATES'
nNodes = int(os.popen('sed -n ''/%s/,/%s/p'' %s | wc -l ' % (bString,eString,geofile_in)).read())
nNodes = nNodes-2
print '--|| ALYA :: WORKING WITH %i NODES IN TOTAL' % nNodes
nNodesPlane = 0
if(caseName == "chan"):
  nNodesPlane = int(os.popen('wc -l < %s.eper.dat'%caseName).read())
nNodesPlane = nNodesPlane + len(open('%s.zper.dat'%caseName).readlines()) - 2
print '--|| ALYA :: WORKING WITH %i NODES ON A PLANE' % nNodesPlane

## generate z-axis slices
if(caseName == "chan"):
  N = int(np.ceil(float(nNodes)/nNodesPlane))
else:
  N = int(np.floor(float(nNodes)/nNodesPlane))
#N = input("\n--|| ALYA: ENTER NUMBER OF Z-PLANES IN THIS SIMULATION..")
print("--|| ALYA: WORKING WITH %d Z-PLANES" % (N))

file0 = os.path.join('./'+caseName+'.ensi.case')
file1 = './AvgData.vtm'

case1 = EnSightReader(CaseFileName=file0)
case1.PointArrays = ['VELOC', 'PRESS']
print "--|| Alya - READING INSTANTANEOUS ARRAYS::", case1.PointArrays
case1.UpdatePipeline()

case2 = XMLMultiBlockDataReader(FileName=file1)
print "--|| Alya - READING AVERAGED ARRAYS::", case2.PointArrayStatus
case2.UpdatePipeline()

startTime = time.time()
print "--|| ALYA: APPEND DATA"
PF1 = ProgrammableFilter(Input=[case1, case2])
### first input is the grid
### the rest of them are data to be averaged
PF1.Script = \
"""
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
#from vtk.numpy_interface import algorithms as algs

def make_tensor(xx,yy,zz, xy, yz, xz):
   t = np.vstack([xx,yy,zz,xy, yz, xz]).transpose().view(dsa.VTKArray)
   t.DataSet = xx.DataSet
   t.Association = xx.Association
   return t

inp0 = inputs[0].PointData["VELOC"]
output.PointData.append(inp0,"V")

inp1 = inputs[0].PointData["AVVEL_average_avg"]
output.PointData.append(inp1,"U")

output.PointData.append(inp0-inp1,"u")
output.PointData.append((inp0-inp1)**2,"u2")

inp0 = inputs[0].PointData["PRESS"]
output.PointData.append(inp0,"PRESS")

inp1 = inputs[0].PointData["AVPRE_average_avg"]
output.PointData.append(inp1,"P")

output.PointData.append(inp0-inp1,"p")
output.PointData.append((inp0-inp1)**2,"p2")

"""

print "--|| ALYA: APPENDING DONE. TIME =",time.time()-startTime,'sec'

#calculate TKE budget terms
startTime = time.time()
print "--|| ALYA: CALCULATE BUDGET TERMS"


varNames = ["C","P","TD","PV","VDD","VD","EPS"]
varAddOn = ["uu","vv","ww"]
count = 1
for var in varNames:
   for i in range(0,3):
      if count==1: 
         pc = PythonCalculator(Input=PF1)
      else:
         strn = 'PythonCalculator%d' % (count-1)
         pc_old = FindSource(strn) 
         pc = PythonCalculator(Input=pc_old)
      count = count+1
      if var == "C":
         pc.Expression = '2*u[:,%d]*(divergence(U*u[:,%d]))' % (i,i)
         pc.ArrayName = '%s_%s' % (var,varAddOn[i])

      elif var == "P":
         pc.Expression = '-2*u[:,%d]*(divergence(U[:,%d]*u))' % (i,i)
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

print "--|| ALYA: TERMS DONE. TIME =",time.time()-startTime,'sec'

strn = 'PythonCalculator%d' % (len(varNames)*len(varAddOn))
case = FindSource(strn) 

case.UpdatePipeline()

nacaensicase = TemporalStatistics(Input=case)

# Properties modified on temporalStatistics1
nacaensicase.ComputeMinimum = 0
nacaensicase.ComputeMaximum = 0
nacaensicase.ComputeStandardDeviation = 0


# Perform spanwise Averaging
(xmin,xmax,ymin,ymax,zmin,zmax) =  case.GetDataInformation().GetBounds()

## generate z-axis slices
resample_transforms=list();
data=list();

slice1 = Slice(Input=nacaensicase)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
## init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
## Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

zmid = (zmin+zmax)/2
zpos = np.arange(N)*(zmax-zmin)/(N-1)

startTime = time.time()
print "--|| ALYA: TRANSFORMATION FOR AVG ON SLICE"

for i in range(N):
	# create a new 'Transform'
	transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
	resampleWithDataset1 = ResampleWithDataset(Input=nacaensicase,Source=transform1,guiName="resample{}".format(i))
	resample_transforms.append(resampleWithDataset1)
	# Properties modified on transform1.Transform
	transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
	#data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])

print "--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec'
HideAll()

startTime = time.time()
print "--|| ALYA: PERFORM SPATIAL AVERAGING"
## create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
### the rest of them are data to be averaged
programmableFilter1.Script = \
"""
import numpy as np

case="BUDGETS"
varFull = []
N=len(inputs);
if(case=="BUDGETS"):
   varNames = ["C","P","TD","PV","VDD","VD","EPS"]
   varAddOn = ["uu","vv","ww"]

print "----|| ALYA WORKING ON ARRAYS ",varNames

extraVarNames = ["u", "p"]
# RMS Variables
for var1 in extraVarNames:
   varName0=var1+"_rms"
   varName=var1+"2_average"
   avg = (inputs[0].PointData[varName])
   for i in range(1,N):
       d = inputs[i].PointData[varName]
       avg = avg + d
   avg = sqrt(abs(avg/N))
   output.PointData.append(avg,varName0)

# Boudget Variable
for var1 in varNames:
   for var2 in varAddOn:
      varName0=var1+"_"+var2
      varName=var1+"_"+var2+"_average"
      avg = (inputs[0].PointData[varName])
      for i in range(1,N):
          d = inputs[i].PointData[varName]
          avg = avg + d
      avg = avg/N
      output.PointData.append(avg,varName0)

"""
programmableFilter1.UpdatePipeline()
print "--|| ALYA: AVERAGING DONE. TIME =",time.time()-startTime,'sec'
##
#### write Files
print "--|| ALYA: WRITE BUDGET FILES"
startTime = time.time()
casePath = os.getcwd()
savePath = casePath+"/BudgetData.csv"
SaveData(savePath, proxy=programmableFilter1)
savePath = casePath+"/BudgetData.vtm"
SaveData(savePath, proxy=programmableFilter1)
print("----|| ALYA: BUDGET FILES WRITTEN TO: ",casePath)
print "--|| ALYA: FILES WRITTEN. TIME =",time.time()-startTime,'sec'

