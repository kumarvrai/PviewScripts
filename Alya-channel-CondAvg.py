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

caseName = sys.argv[1]
task_label = sys.argv[2]
airfoil = sys.argv[3]
dim      = sys.argv[4]   #What is the dimension of the output file?
mode    = sys.argv[5]    #How is the averaging performed? [SORT/INTERP/PINTERP]
nCondTop       = 9
nCondBottom    = 8

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

print "--|| ALYA :: READING ALYA-AVERAGED ARRAYS"
startTime = time.time()

fileName = 'AvgData_3D.vtm'

fullcase = OpenDataFile(fileName)
fullcase.UpdatePipeline()
print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
(xmin,xmax,ymin,ymax,zmin,zmax) =  fullcase.GetDataInformation().GetBounds()
Ntotal = int(fullcase.GetDataInformation().GetNumberOfPoints())

if('3D' in dim):
 print "--|| ALYA :: CREATING SLICE"
 startTime = time.time()
 slicecase = Slice(Input=fullcase)
 slicecase.SliceType = 'Plane'
 slicecase.SliceOffsetValues = [0.0]
 slicecase.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
 slicecase.SliceType.Normal = [0.0, 0.0, 1.0]
 slicecase.UpdatePipeline()
 print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
 Nplane = int(slicecase.GetDataInformation().GetNumberOfPoints())
 N = int(Ntotal/Nplane)
 print('--|| ALYA :: NUMBER OF Z-PLANES',N) 
 zpos = np.arange(N)*(zmax-zmin)/(N-1)

 if("TDA" in task_label):
   nCondTop = N
   print "--|| ALYA :: PHASE AVERAGING BEGINS"
   startTime = time.time()
   PF1 = ProgrammableFilter(Input=[fullcase])
   PF1.Script = \
"""
import numpy as np
N0 = N-(N%nCondTop)
nCondFreqTop = N0//nCondTop;
print('--|| ALYA :: PERIODICITY ON SUCTION SIDE WITH',nCondTop,'SEGMENTS AND FREQ ',nCondFreqTop) 
nRepMat = np.tile(np.arange(N0),N+1);
#print(nRepMat)
varFull = inputs[0].PointData.keys()
d = dsa.WrapDataObject(inputs[0].GetBlock(0))
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=6)
for varName in varFull:
 print("--|| ALYA :: CALCULATING FOR",varName) 
 src = np.asarray(d.PointData[varName],dtype=np.double)
 avg = 0.0*np.asarray(d.PointData[varName],dtype=np.double)
 for n in range(N):
  ind_z_main = np.where(z == np.around(zpos[n],decimals=6))
  x_2d = x[ind_z_main];
  y_2d = y[ind_z_main]
  avg_2d = avg[ind_z_main]
  ind_2d = np.lexsort((y_2d,x_2d)); 
  ind_cond = np.arange(n,N0+n,nCondFreqTop);
  z_cond_loc = zpos[nRepMat[ind_cond]];
  #print(z_cond_loc)
  if("AVVEL" in varName):
   print("FOR Z=",zpos[n],"AVERAGING",z_cond_loc)
  M = len(z_cond_loc)
  ind_count = np.asarray(0*ind_2d,dtype=int);
  for m in range(M):
   ind_z = np.where(z == np.around(z_cond_loc[m],decimals=6))
   x_l = x[ind_z];
   y_l = y[ind_z];
   s_src = src[ind_z]
   if("SORT" in mode):
    ind = np.lexsort((y_l,x_l)); 
    ind_u = np.logical_and(np.equal(x_l[ind],x_2d[ind_2d]),np.equal(y_l[ind],y_2d[ind_2d]))
    ind_v = np.where(ind_u == False)
    s_src[ind_v] = 0.0
    avg_2d[ind_2d] = avg_2d[ind_2d] + s_src[ind]
    ind_count = ind_count + ind_u.astype(int)
   elif("INTERP" in mode):
    avg_2d = avg_2d + griddata((x_l,y_l),s_src, (x_2d,y_2d), method='nearest')
   else:
    raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')
  #############
  if("SORT" in mode):
   if("AVVEL" in varName): 
    print("--|| PLANE",n,"MAX",np.amax(ind_count,axis=None),"MIN",np.amin(ind_count,axis=None),\
    "AVERAGE",np.mean(ind_count,axis=None)) 
   arr_type = len(avg_2d.shape) 
   #avg[ind_z_main] = avg_2d/M
   if(arr_type == 1):
    avg[ind_z_main] = np.divide(avg_2d,ind_count);
   elif(arr_type == 2):
    num_rows,num_cols = avg_2d.shape
    avg[ind_z_main] = np.divide(avg_2d,np.tile(ind_count,(num_cols,1)).T);
   elif(arr_type == 3):
    num_rows, num_cols, num_3d = avg_2d.shape
    avg[ind_z_main] = np.divide(avg_2d,np.tile(ind_count,(num_cols,num_3d,1)).T);
  elif("INTERP" in mode):
   avg[ind_z_main]=avg_2d/M;
  else:
   raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')

 avg = np.asarray(avg, dtype=np.float64)
 output.ShallowCopy(inputs[0].VTKObject)
 output.PointData.append(avg, "DA_"+varName)

"""
   PF1.UpdatePipeline()
   print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
   ########### CLIP DATA IN A BOX [AIRFOIL ONLY]##############
   print "--|| ALYA :: CLIPPING DATA IN A BOX"
   startTime = time.time()
   #fullcase = Clip(Input=PF1)
   #fullcase.ClipType = 'Box'
   #if "4412" in airfoil:
   # print("--|| ALYA :: WORKING WITH",airfoil)
   # fullcase.ClipType.Position = [0.5, 0.03, 0.0]
   # fullcase.ClipType.Scale = [0.026, 0.005, 1.0]
   #elif "0012" in airfoil:
   # print("--|| ALYA :: WORKING WITH",airfoil)
   # fullcase.ClipType.Position = [0.5, 0.0, 0.1]
   # fullcase.ClipType.Scale = [0.026, 0.005, 1.0]
   #else:
   # raise Exception("--|| ERROR :: AIRFOIL NOT DEFINED")
   #fullcase.UpdatePipeline()
   #print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
   print "--|| ALYA :: SAVING CLIPPED DATA"
   startTime = time.time()
   savePath = casePath+"/AvgData_3D_TDA_box.vtm"
   SaveData(savePath, proxy=PF1)
   print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
  
 else:
   # ########### CLIP DATA EACH FOR SS AND PS##############
   print "--|| ALYA :: CLIPPING FOR SS"
   startTime = time.time()
   caseSS = Clip(Input=fullcase)
   caseSS.ClipType.Normal = [0.0, 1.0, 0.0]
   caseSS.ClipType.Origin = [0.0, 0.0, 0.1]
   caseSS.Invert = 0;
   caseSS.UpdatePipeline()
   print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
   print "--|| ALYA :: CLIPPING FOR PS"
   startTime = time.time()
   casePS = Clip(Input=fullcase)
   casePS.ClipType.Normal = [0.0, 1.0, 0.0]
   casePS.ClipType.Origin = [0.0, 0.0, 0.1]
   casePS.Invert = 1;
   casePS.UpdatePipeline()
   print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'

   print "--|| ALYA :: PHASE AVERAGING USING",mode
   print "--|| ALYA :: PHASE AVERAGING FOR SUCTION SIDE"
   startTime = time.time()
   PF1 = ProgrammableFilter(Input=[caseSS])
   PF1.Script = \
"""
import numpy as np
N0 = N-(N%nCondTop)
nCondFreqTop = N0//nCondTop;
print('--|| ALYA :: PERIODICITY ON SUCTION SIDE WITH',nCondTop,'SEGMENTS AND FREQ ',nCondFreqTop) 
nRepMat = np.tile(np.arange(N0),N+1);
#print(nRepMat)
varFull = inputs[0].PointData.keys()
d = dsa.WrapDataObject(inputs[0].GetBlock(0))
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=6)
for varName in varFull:
 print("--|| ALYA :: CALCULATING FOR",varName) 
 src = np.asarray(d.PointData[varName],dtype=np.double)
 avg = 0.0*np.asarray(d.PointData[varName],dtype=np.double)
 for n in range(N):
  ind_z_main = np.where(z == np.around(zpos[n],decimals=6))
  x_2d = x[ind_z_main];
  y_2d = y[ind_z_main]
  avg_2d = avg[ind_z_main]
  ind_2d = np.lexsort((y_2d,x_2d)); 
  ind_cond = np.arange(n,N0+n,nCondFreqTop);
  z_cond_loc = zpos[nRepMat[ind_cond]];
  if("AVVEL" in varName):
   print("FOR Z=",zpos[n],"AVERAGING",z_cond_loc)
  M = len(z_cond_loc)
  ind_count = np.asarray(0*ind_2d,dtype=int);
  for m in range(M):
   ind_z = np.where(z == np.around(z_cond_loc[m],decimals=6))
   x_l = x[ind_z];
   y_l = y[ind_z];
   s_src = src[ind_z]
   if("SORT" in mode):
    ind = np.lexsort((y_l,x_l)); 
    ind_u = np.logical_and(np.equal(x_l[ind],x_2d[ind_2d]),np.equal(y_l[ind],y_2d[ind_2d]))
    ind_v = np.where(ind_u == False)
    s_src[ind_v] = 0.0
    avg_2d[ind_2d] = avg_2d[ind_2d] + s_src[ind]
    ind_count = ind_count + ind_u.astype(int)
   elif("INTERP" in mode):
    avg_2d = avg_2d + griddata((x_l,y_l),s_src, (x_2d,y_2d), method='nearest')
   else:
    raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')
  #############
  if("SORT" in mode):
   if("AVVEL" in varName): 
    print("--|| PLANE",n,"MAX",np.amax(ind_count,axis=None),"MIN",np.amin(ind_count,axis=None),\
    "AVERAGE",np.mean(ind_count,axis=None)) 
   arr_type = len(avg_2d.shape) 
   #avg[ind_z_main] = avg_2d/M
   if(arr_type == 1):
    avg[ind_z_main] = np.divide(avg_2d,ind_count);
   elif(arr_type == 2):
    num_rows,num_cols = avg_2d.shape
    avg[ind_z_main] = np.divide(avg_2d,np.tile(ind_count,(num_cols,1)).T);
   elif(arr_type == 3):
    num_rows, num_cols, num_3d = avg_2d.shape
    avg[ind_z_main] = np.divide(avg_2d,np.tile(ind_count,(num_cols,num_3d,1)).T);
  elif("INTERP" in mode):
   avg[ind_z_main]=avg_2d/M;
  else:
   raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')

 avg = np.asarray(avg, dtype=np.float64)
 output.ShallowCopy(inputs[0].VTKObject)
 output.PointData.append(avg, varName)

"""
   PF1.UpdatePipeline()
   print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
   print "--|| ALYA :: PHASE AVERAGING FOR PRESSURE SIDE"
   startTime = time.time()
   PF2 = ProgrammableFilter(Input=[casePS])
   PF2.Script = \
"""
import numpy as np

N0 = N-(N%nCondBottom)
nCondFreqBottom = N0//nCondBottom;
print('--|| ALYA :: PERIODICITY ON PRESSURE SIDE WITH',nCondBottom,'SEGMENTS AND FREQ ', nCondFreqBottom) 
nRepMat = np.tile(np.arange(N0),N+1);
varFull = inputs[0].PointData.keys()
d = dsa.WrapDataObject(inputs[0].GetBlock(0))
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=6)
for varName in varFull:
 print("--|| ALYA :: CALCULATING FOR",varName) 
 src = np.asarray(d.PointData[varName],dtype=np.double)
 avg = 0.0*np.asarray(d.PointData[varName],dtype=np.double)
 for n in range(N):
  ind_z_main = np.where(z == np.around(zpos[n],decimals=6))
  x_2d = x[ind_z_main];
  y_2d = y[ind_z_main]
  avg_2d = avg[ind_z_main]
  ind_2d = np.lexsort((y_2d,x_2d)); 
  ind_cond = np.arange(n,N0+n,nCondFreqBottom);
  z_cond_loc = zpos[nRepMat[ind_cond]];
  M = len(z_cond_loc)
  ind_count = np.asarray(0*ind_2d,dtype=int);
  for m in range(M):
   ind_z = np.where(z == np.around(z_cond_loc[m],decimals=6))
   x_l = x[ind_z];
   y_l = y[ind_z];
   s_src = src[ind_z]
   if("SORT" in mode):
    ind = np.lexsort((y_l,x_l)); 
    ind_u = np.logical_and(np.equal(x_l[ind],x_2d[ind_2d]),np.equal(y_l[ind],y_2d[ind_2d]))
    ind_v = np.where(ind_u == False)
    s_src[ind_v] = 0.0
    avg_2d[ind_2d] = avg_2d[ind_2d] + s_src[ind]
    ind_count = ind_count + ind_u.astype(int)
   elif("INTERP" in mode):
    avg_2d = avg_2d + griddata((x_l,y_l),s_src, (x_2d,y_2d), method='nearest')
   else:
    raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')
  #############
  if("SORT" in mode):
   if("AVVEL" in varName): 
    print("--|| PLANE",n,"MAX",np.amax(ind_count,axis=None),"MIN",np.amin(ind_count,axis=None),\
    "AVERAGE",np.mean(ind_count,axis=None)) 
   arr_type = len(avg_2d.shape) 
   #avg[ind_z_main] = avg_2d/M
   if(arr_type == 1):
    avg[ind_z_main] = np.divide(avg_2d,ind_count);
   elif(arr_type == 2):
    num_rows,num_cols = avg_2d.shape
    avg[ind_z_main] = np.divide(avg_2d,np.tile(ind_count,(num_cols,1)).T);
   elif(arr_type == 3):
    num_rows, num_cols, num_3d = avg_2d.shape
    avg[ind_z_main] = np.divide(avg_2d,np.tile(ind_count,(num_cols,num_3d,1)).T);
  elif("INTERP" in mode):
   avg[ind_z_main]=avg_2d/M;
  else:
   raise ValueError('--|| ALYA ERROR :: HOW DO YOU WANT TO CALCULATE AVERGAES?')

 avg = np.asarray(avg, dtype=np.float64)
 output.ShallowCopy(inputs[0].VTKObject)
 output.PointData.append(avg, varName)

"""
   PF2.UpdatePipeline()
   print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
   ### GROUP BOTH SIDES TOGTHER ###
   print "--|| ALYA :: GROUP DATASETS"
   startTime = time.time()
   GD1 = GroupDatasets(Input=[PF1,PF2])
   GD1.UpdatePipeline()
   print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
   ### SAVE-FILE ###
   savePath = casePath+"/AvgData_3D_COND.vtm"
   SaveData(savePath, proxy=GD1)
   print "--|| ALYA: 3D STATISTICS FILE WRITTEN AS: ",savePath
   print "--|| ALYA: FILE SAVED. TIME =",time.time()-startTime,'sec'

if('2D' in dim):
  geomName = 'RLE2'
  loc_str = ['PSS','VSS','PPS','VPS']
  #loc_str = ['PSS']
  
  if("RLE1" in geomName):
    peak_loc_ss = [0, 0.000755680000000000, 0.00138890000000000, 0.00214460000000000, 0.00277780000000000, 0.00353350000000000, 0.0174220000000000, 0.0180560000000000, 0.0188110000000000, 0.0194440000000000, 0.0202000000000000, 0.0208330000000000, 0.0215890000000000, 0.0222220000000000, 0.0229780000000000, 0.0236110000000000, 0.0243670000000000, 0.0250000000000000, 0.0257560000000000, 0.0263890000000000, 0.0271450000000000, 0.0680560000000000, 0.0688110000000000, 0.0694440000000000, 0.0702000000000000, 0.0708330000000000, 0.0715890000000000, 0.0882560000000000, 0.0888890000000000, 0.0896450000000000, 0.104920000000000, 0.105560000000000, 0.106310000000000, 0.106940000000000, 0.107700000000000, 0.108330000000000, 0.109090000000000, 0.109720000000000, 0.110480000000000, 0.111110000000000, 0.111870000000000, 0.112500000000000, 0.129170000000000, 0.129920000000000, 0.130560000000000, 0.172980000000000, 0.173610000000000, 0.174370000000000, 0.175000000000000, 0.175760000000000, 0.176390000000000, 0.177140000000000, 0.177780000000000, 0.178530000000000, 0.179170000000000, 0.195830000000000, 0.196590000000000, 0.197220000000000, 0.197980000000000, 0.198610000000000, 0.199370000000000, 0.200000000000000]
    valley_loc_ss =[0.00833330000000000, 0.00908900000000000, 0.00972220000000000, 0.0104780000000000, 0.0111110000000000, 0.0118670000000000, 0.0333330000000000, 0.0340890000000000, 0.0347220000000000, 0.0354780000000000, 0.0361110000000000, 0.0368670000000000, 0.0375000000000000, 0.0382560000000000, 0.0388890000000000, 0.0938110000000000, 0.0944440000000000, 0.0952000000000000, 0.0958330000000000, 0.0965890000000000, 0.0972220000000000, 0.0979780000000000, 0.0986110000000000, 0.0993670000000000, 0.100000000000000, 0.100760000000000, 0.101390000000000, 0.116030000000000, 0.116670000000000, 0.117420000000000, 0.118060000000000, 0.118810000000000, 0.119440000000000, 0.120200000000000, 0.120830000000000, 0.121590000000000, 0.122220000000000, 0.122980000000000, 0.123610000000000, 0.138260000000000, 0.138890000000000, 0.139640000000000, 0.140280000000000, 0.183330000000000, 0.184090000000000, 0.184720000000000, 0.185480000000000, 0.186110000000000, 0.186870000000000, 0.187500000000000, 0.188260000000000, 0.188890000000000, 0.189640000000000, 0.190280000000000, 0.191030000000000, 0.191670000000000]
    peak_loc_ps = [0.200000000000000, 0.198890000000000, 0.198610000000000, 0.197480000000000, 0.197220000000000, 0.179170000000000, 0.178040000000000, 0.177780000000000, 0.176650000000000, 0.176390000000000, 0.175260000000000, 0.175000000000000, 0.173870000000000, 0.173610000000000, 0.172480000000000, 0.172220000000000, 0.171090000000000, 0.153040000000000, 0.152780000000000, 0.151650000000000, 0.151390000000000, 0.131940000000000, 0.129430000000000, 0.129170000000000, 0.128040000000000, 0.127780000000000, 0.126650000000000, 0.126390000000000, 0.102780000000000, 0.101650000000000, 0.101390000000000, 0.100260000000000, 0.100000000000000, 0.0988720000000000, 0.0986110000000000, 0.0974830000000000, 0.0972220000000000, 0.0960940000000000, 0.0958330000000000, 0.0780380000000000, 0.0777780000000000, 0.0766500000000000, 0.0763890000000000, 0.0752610000000000, 0.0750000000000000, 0.0738720000000000, 0.0736110000000000, 0.0513890000000000, 0.0502610000000000, 0.0500000000000000, 0.0488720000000000, 0.0486110000000000, 0.0474830000000000, 0.0472220000000000, 0.0308160000000000, 0.0305560000000000, 0.0294270000000000, 0.0291670000000000, 0.0280380000000000, 0.0277780000000000, 0.0266500000000000, 0.0263890000000000, 0.0252610000000000, 0.0250000000000000, 0.00416670000000000, 0.00303840000000000, 0.00277780000000000, 0.00164950000000000, 0.00138890000000000, 0.000266120000000000, 0]
    valley_loc_ps = [0.191930000000000, 0.191670000000000, 0.190540000000000, 0.190280000000000, 0.189150000000000, 0.188890000000000, 0.187760000000000, 0.187500000000000, 0.186370000000000, 0.186110000000000, 0.184980000000000, 0.184720000000000, 0.162760000000000, 0.162500000000000, 0.161370000000000, 0.161110000000000, 0.159980000000000, 0.159720000000000, 0.158590000000000, 0.158330000000000, 0.139150000000000, 0.138890000000000, 0.137760000000000, 0.137500000000000, 0.109980000000000, 0.109720000000000, 0.0919270000000000, 0.0916670000000000, 0.0905380000000000, 0.0902780000000000, 0.0891500000000000, 0.0888890000000000, 0.0877610000000000, 0.0875000000000000, 0.0863720000000000, 0.0861110000000000, 0.0849830000000000, 0.0847220000000000, 0.0835940000000000, 0.0833330000000000, 0.0419270000000000, 0.0416670000000000, 0.0405380000000000, 0.0402780000000000, 0.0391500000000000, 0.0388890000000000, 0.0377610000000000, 0.0375000000000000, 0.0363720000000000, 0.0361110000000000, 0.0349830000000000, 0.0347220000000000, 0.00998290000000000, 0.00972220000000000, 0.00859400000000000, 0.00833330000000000]
  
  elif("RLE2" in geomName):
    peak_loc_ss = [0.0500000000000000,0.0502960000000000,0.0513890000000000,0.0516850000000000,0.0527780000000000,0.0530740000000000,0.0541670000000000,0.0544630000000000,0.0555560000000000,0.0558520000000000,0.0569440000000000,0.0572410000000000,0.0583330000000000,0.0586300000000000,0.0597220000000000,0.0600190000000000,0.0611110000000000,0.0614070000000000,0.0625000000000000,0.0627960000000000,0.0638890000000000,0.0641850000000000,0.0652780000000000,0.0655740000000000,0.0666670000000000,0.0669630000000000,0.0680560000000000,0.143060000000000,0.143350000000000,0.144440000000000,0.144740000000000,0.145830000000000,0.146130000000000,0.161410000000000,0.162500000000000,0.162800000000000,0.163890000000000,0.164190000000000,0.165280000000000,0.165570000000000,0.166670000000000,0.166960000000000,0.168060000000000,0.168350000000000,0.169440000000000,0.169740000000000,0.170830000000000,0.171130000000000,0.172220000000000,0.172520000000000,0.173610000000000,0.173910000000000,0.175000000000000,0.175300000000000,0.176390000000000,0.176690000000000,0.177780000000000,0.178070000000000,0.179170000000000]
    valley_loc_ss =[0,0.000296380000000000,0.00416670000000000,0.00446310000000000,0.00555560000000000,0.00585190000000000,0.00694440000000000,0.00724080000000000,0.00833330000000000,0.00862970000000000,0.00972220000000000,0.0100190000000000,0.0111110000000000,0.0114070000000000,0.0125000000000000,0.0127960000000000,0.0138890000000000,0.0141850000000000,0.0152780000000000,0.0975190000000000,0.0986110000000000,0.0989070000000000,0.100000000000000,0.100300000000000,0.101390000000000,0.101690000000000,0.102780000000000,0.103070000000000,0.104170000000000,0.104460000000000,0.105560000000000,0.127780000000000,0.128070000000000,0.129170000000000,0.129460000000000,0.130560000000000,0.130850000000000,0.131940000000000,0.132240000000000,0.133330000000000,0.133630000000000,0.191670000000000,0.191960000000000,0.193060000000000,0.193350000000000,0.194440000000000,0.194740000000000,0.195830000000000,0.196130000000000,0.197220000000000,0.197520000000000,0.198610000000000,0.198910000000000,0.200000000000000]
    peak_loc_ps = [0.105560000000000,0.105820000000000,0.106940000000000,0.107210000000000,0.108330000000000,0.108590000000000,0.109720000000000,0.109980000000000,0.111110000000000,0.111370000000000,0.112500000000000,0.112760000000000,0.113890000000000,0.114150000000000,0.115280000000000,0.115540000000000,0.136110000000000,0.136370000000000,0.137500000000000,0.137760000000000,0.138890000000000,0.139150000000000,0.140280000000000,0.140540000000000,0.141670000000000,0.141930000000000,0.143060000000000,0.143320000000000,0.144440000000000,0.144710000000000,0.145830000000000,0.146090000000000,0.147220000000000,0.147480000000000,0.148610000000000,0.148870000000000,0.150000000000000,0.150260000000000,0.162500000000000,0.162760000000000,0.163890000000000,0.164150000000000,0.165280000000000,0.165540000000000,0.166670000000000,0.166930000000000,0.168060000000000,0.168320000000000,0.169440000000000,0.169710000000000,0.170830000000000,0.171090000000000,0.172220000000000,0.172480000000000,0.173610000000000,0.173870000000000,0.175000000000000,0.175260000000000,0.176390000000000]
    valley_loc_ps = [0.00998290000000000,0.0111110000000000,0.0113720000000000,0.0125000000000000,0.0127610000000000,0.0138890000000000,0.0141500000000000,0.0152780000000000,0.0155380000000000,0.0166670000000000,0.0169270000000000,0.0180560000000000,0.0183160000000000,0.0194440000000000,0.0197050000000000,0.0208330000000000,0.0210940000000000,0.0377610000000000,0.0388890000000000,0.0391500000000000,0.0402780000000000,0.0405380000000000,0.0416670000000000,0.0419270000000000,0.0430560000000000,0.0433160000000000,0.0444440000000000,0.0447050000000000,0.0458330000000000,0.0460940000000000,0.0472220000000000,0.0474830000000000,0.0486110000000000,0.0488720000000000,0.0500000000000000,0.0597220000000000,0.0599830000000000,0.0611110000000000,0.0613720000000000,0.0625000000000000,0.0627610000000000,0.0638890000000000,0.0641500000000000,0.0652780000000000,0.0655380000000000,0.0666670000000000,0.0669270000000000,0.0680560000000000,0.0683160000000000,0.0694440000000000,0.0697050000000000,0.0708330000000000,0.0710940000000000,0.0722220000000000,0.0724830000000000]
  else:
    raise Exception("--|| ERROR :: GEOMETRY DOESN'T EXIST!")
  
  if "CONDAVG" in task_label:
   for sideLabel in loc_str:
     print "--|| ALYA :: GENERATING SLICE FOR SPANWISE AVERAGE"
     startTime = time.time()
     case = Clip(Input=fullcase)
     case.ClipType.Normal = [0.0, 1.0, 0.0]
     case.ClipType.Origin = [0.0, 0.0, 0.1]
     if('SS' in sideLabel):
      case.Invert = 0;
     elif('PS' in sideLabel):
      case.Invert = 1;
     else:
      raise Exception("--|| ERROR :: NO SIDE CHOSEN")
     case.UpdatePipeline()
     print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
     
     ########### PEAK AVERAGING SUCTION SIDE################
     (xmin,xmax,ymin,ymax,zmin,zmax) =  case.GetDataInformation().GetBounds()
     print('--|| CLIPPED :: XMIN=',xmin,'XMAX=',xmax,'YMIN=',ymin,'YMAX=',ymax,'ZMIN=',zmin,'ZMAX=',zmax) 
     print "--|| ALYA :: GENERATING SLICE FOR %s AVERAGE" % (sideLabel)
     startTime = time.time()
     
     if('PSS' in sideLabel):
      zpos = np.asarray(peak_loc_ss)
     elif('VSS' in sideLabel):
      zpos = np.asarray(valley_loc_ss)
     elif('PPS' in sideLabel):
      zpos = np.asarray(peak_loc_ps)
     elif('VPS' in sideLabel):
      zpos = np.asarray(valley_loc_ps)
     else:
      raise Exception("--|| ERROR :: DIMENSIONS COULDN''T BE DETERMINED")
   
     N = len(zpos)
     print("--|| ALYA: WORKING WITH %d Z-PLANES" % (N))
     
     resample_transforms=list();
     data=list();
     
     slice1 = Slice(Input=case)
     slice1.SliceType = 'Plane'
     slice1.SliceOffsetValues = [0.0]
     slice1.SliceType.Origin = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
     slice1.SliceType.Normal = [0.0, 0.0, 1.0]
     
     zmid = (zmin+zmax)/2
     print "--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
       
     print "--|| ALYA :: CREATING TRANSFORMATIONS"
     startTime = time.time()
     
     for i in range(N):
     	# create a new 'Transform'
     	transform1 = Transform(Input=slice1,guiName="transform{}".format(i))
     	resampleWithDataset1 = ResampleWithDataset(Input=case,Source=transform1,guiName="resample{}".format(i))
     	resample_transforms.append(resampleWithDataset1)
     	# Properties modified on transform1.Transform
     	transform1.Transform.Translate = [0.0, 0.0, zpos[i]-zmid]
     	#data.append(dsa.WrapDataObject(servermanager.Fetch(resampleWithDataset1)).PointData['AVVEL'])
     print "--|| ALYA: TRANSFORMATION DONE. TIME =",time.time()-startTime,'sec'
     HideAll()
     
     ## create a new 'Programmable Filter'
     print "--|| ALYA: AVERAGING USING A PROGRAMMABLE FILTER"
     startTime = time.time()
     
     PF1 = ProgrammableFilter(Input=[slice1]+resample_transforms)
     ### first input is the grid
     ### the rest of them are data to be averaged
     PF1.Script = \
   """
   import numpy as np
   
   varFull = []
   varFull = inputs[0].PointData.keys()
   print "----|| Alya - WORKING ON ARRAYS::",varFull
   N=len(inputs)-1;
   print "--|| ALYA: AVERAGING %d DATA-PLANES" % (N)
   
   for varName in varFull:
      varName0=varName[0:5]
      avg = 0.0*(inputs[0].PointData[varName])
      for i in range(N):
          d = inputs[i+1].PointData[varName]
          avg = avg + d
      avg = avg/N
      output.PointData.append(avg,varName0)
   """
     PF1.UpdatePipeline()
     print "--|| ALYA: SPANWISE AVERAGING DONE. TIME =",time.time()-startTime,'sec'
     savePath = casePath+"/AvgData_2D_"+sideLabel+".vtm"
     SaveData(savePath, proxy=PF1)
     savePath = casePath+"/AvgData_2D_"+sideLabel+".csv"
     SaveData(savePath, proxy=PF1)
     print "----|| ALYA: 2D STATISTICS FILE WRITTEN AS: ",savePath
  elif "RGHSTAT" in task_label:
   # LOAD AIRFOIL SPECIFIC FILES
   niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'
   niaHome = '/home/u/ugo/kvishal/'
   if('0012' in airfoil):
     fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-UP.txt'
     fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-DOWN.txt'
   elif('4412' in airfoil):  
     fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-UP.txt'
     fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-DOWN.txt'
   else:
     sys.exit("--||Alya (ERROR): AIRFOIL NOT IN THE LIST")
   

   print "--|| ALYA :: CALCULATING ROUGHNESS STATISTICS AT PEAK VALLEY LOCATION"
   print "--|| ALYA :: CLIPPING ONLY AIRFOIL SURFACE"
   ES1 = ExtractSurface(Input=fullcase)
   ES1.UpdatePipeline()
   CLIP1 = Clip(Input=ES1)
   CLIP1.ClipType = 'Box'
   CLIP1.ClipType.Scale = [0.5, 0.5, 1.0]
   CLIP1.UpdatePipeline()

   for sideLabel in loc_str:
     print "----|| ALYA :: CLIPPING SS/PS VOLUME"
     startTime = time.time()
     CLIP2 = Clip(Input=CLIP1)
     CLIP2.ClipType.Normal = [0.0, 1.0, 0.0]
     CLIP2.ClipType.Origin = [0.0, 0.0, 0.1]
     if('SS' in sideLabel):
      CLIP2.Invert = 0;
      coordAir = np.loadtxt(fAirU, delimiter=',')
     elif('PS' in sideLabel):
      CLIP2.Invert = 1;
      coordAir = np.loadtxt(fAirL, delimiter=',')
     else:
      raise Exception("--|| ERROR :: NO SIDE CHOSEN")
     CLIP2.UpdatePipeline()
     F1 = interp1d(coordAir[:,0],coordAir[:,1])
     print "----|| ALYA :: DONE. TIME =",time.time()-startTime,'sec'
     if('PSS' in sideLabel):
      zpos = np.asarray(peak_loc_ss)
     elif('VSS' in sideLabel):
      zpos = np.asarray(valley_loc_ss)
     elif('PPS' in sideLabel):
      zpos = np.asarray(peak_loc_ps)
     elif('VPS' in sideLabel):
      zpos = np.asarray(valley_loc_ps)
     else:
      raise Exception("--|| ERROR :: DIMENSIONS COULDN''T BE DETERMINED")
   
     N = len(zpos)
     print("----|| ALYA: WORKING WITH %d Z-PLANES" % (N))
     print "--|| ALYA :: CALCULATE ROUGHNESS-STATS USING PFILTER"
     startTime = time.time()
     PF2 = ProgrammableFilter(Input=[CLIP2])
     PF2.Script = \
    """
    import numpy as np
    def find_nearest(array, value):
     array = np.asarray(array)
     idx = (np.abs(array - value)).argmin()
     return array[idx]
    
    # LOAD PLANAR DATA
    d = dsa.WrapDataObject(inputs[0].GetBlock(0))
    x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
    y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
    z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=6)
    zpos = np.around(zpos,decimals=6)
    # FOR EVERY ZPLANE: INTERPOLATE SMOOTH PROFILE ON ROUGH GRID
    z_loc = []; r_nrm = [];
    for n in range(N):
     ind_z_main = np.where(z == find_nearest(z,zpos[n]))
     x_2d = x[ind_z_main];
     y_2d = y[ind_z_main]
     ind_2d = np.argsort(x_2d); 
     x_2d = x_2d[ind_2d]; 
     y_2d = y_2d[ind_2d];
     y_2d_smth = F1(x_2d);

     # FOR EVERY ZPLANE: CALCULATE ROUGHNESS STATISTICS
     z_loc.append(zpos[n]);
     header_str = 'Z_LOC,  '
     r_nrm.append(np.linalg.norm((y_2d), axis=None))
     header_str += 'R_NORM'

    ## FOR EVERY ZPLANE: APPEND THE ROUGHNESS STATISTICS
    savePath = casePath+"/AvgData_2D_"+sideLabel+"_RSTAT.csv"
    np.savetxt(savePath, np.c_[z_loc,r_nrm], delimiter=',',header = header_str)
    
    """
     PF2.UpdatePipeline()
     print "--|| DONE. TIME =",time.time()-startTime,'sec'
     #savePath = casePath+"/AvgData_2D_"+sideLabel+"_RSTAT.csv"
     #SaveData(savePath, proxy=PF2)
     #print "----|| ALYA: 2D ROUGH-STAT FILE WRITTEN AS: ",savePath
  else:
   raise Exception("--|| ERROR :: LABEL DOESN'T EXIST!")
