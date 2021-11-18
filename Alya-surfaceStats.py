#### import the simple module from the paraview
import os
import time
import operator
import numpy as np
from paraview import vtk
from paraview import numpy_support
from paraview.simple import *
from vtk.numpy_interface import dataset_adapter as da
from vtk.numpy_interface.algorithms import sqrt as sqrt

caseName	= sys.argv[1]

casePath = os.getcwd()

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

print "--|| ALYA :: READING ALYA-GEOMETRY"
startTime = time.time()
s_path = 
r_path = 
## SMOOTH CASE
print "----|| READING :: SMOOTH ALYA-GEOMETRY"
fileName = smth_path+caseName+'.ensi.case'
s_case = OpenDataFile(fileName)
s_case.PointArrays = []
s_case.UpdatePipeline()
## ROUGH CASE
print "----|| READING :: ROUGH ALYA-GEOMETRY"
fileName = r_path+caseName+'.ensi.case'
r_case = OpenDataFile(fileName)
r_case.PointArrays = []
r_case.UpdatePipeline()
print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'

pts1 = inputs[1].Points
pts0 = inputs[0].Points
output.PointData.append(abs(pts0-pts1),'h')

F = inputs[0].GetBlock(0).GetPoints().GetData()
G = inputs[1].GetBlock(0).GetPointData().GetArray('s_local')

n = F.GetNumberOfTuples()
print("Points=",n)
ind=[]; x = []; y = []; z = [];
for j in range(0, n):
  f = F.GetTuple(j)
    g = G.GetTuple(j)
      ind.append(j);
        x.append(f[0]); y.append(f[1]);
	  z.append(g[0]);

	  indSort = np.lexsort((y,x));
	  x = x[indSort];
	  #print(shape(arr))
	  #vtk_arr = da.VTKArray(arr)
	  #output.PointData.append(vtk_arr,"A_IJ")
	  inp0 = inputs[1].PointData["s_local"]
	  output.PointData.append(inp0,'s_local')

