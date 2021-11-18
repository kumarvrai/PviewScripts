#### import the simple module from the paraview
import os
import time
import operator
import numpy as np
import glob
from scipy.interpolate import interp1d
#from scipy.interpolate import interp2d
from scipy.interpolate import griddata
from paraview import vtk
from paraview import numpy_support
from paraview.vtk.numpy_interface import dataset_adapter as dsa
import paraview.vtk.util.numpy_support as vnp
from vtk.util import numpy_support
from paraview.simple import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtk.numpy_interface import algorithms as algs
import matplotlib.pyplot as plt

ls = ['-','--','-.',':','-o','-v','-^','-s','-p','-d','-*','-+','-x']
plt.close("all")
lw = 2.0;
SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
casePath = os.getcwd() 

mode     = sys.argv[1]
caseName = sys.argv[2] 
airfoil  = sys.argv[3] ; 

sides = ['SS','PS']
d = 0.1

niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'
niaHome = '/home/u/ugo/kvishal/'
# LOAD AIRFOIL SPECIFIC FILES
if('0012' in airfoil):
  fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-UP.txt'
  fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-DOWN.txt'
elif('4412' in airfoil):  
  fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-UP.txt'
  fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-DOWN.txt'
else:
  raise ValueError('--|| ALYA ERROR :: FILE NOT PROVIDED.')

print('--|| ALYA :: READING AIRFOIL DATA.')
stime = time.time()

coordAirU = np.loadtxt(fAirU, delimiter=',')
coordAirL = np.loadtxt(fAirL, delimiter=',')

# INTERPOLATION ON SUCTION SIDE
thAirU = np.arctan2(np.diff(coordAirU[:,1]),np.diff(coordAirU[:,0]))
F0 = interp1d(coordAirU[:,0],coordAirU[:,1])
airLen = len(coordAirU)
coordMid = 0.5*(coordAirU[0:airLen-1,0]+coordAirU[1:airLen,0])
Fth0 = interp1d(thAirU,coordMid)
dydx0 = Fth0(0.0).tolist()
print('----|| ALYA :: SUCTION SIDE AT X =',dydx0,' dy/dx = 0')
Fth0 = interp1d(coordMid,thAirU)

# INTERPOLATION ON PRESSURE SIDE
thAirL = np.arctan2(np.diff(coordAirL[:,1]),np.diff(coordAirL[:,0]))
F1 = interp1d(coordAirL[:,0],coordAirL[:,1])
airLen = len(coordAirL)
coordMid = 0.5*(coordAirL[0:airLen-1,0]+coordAirL[1:airLen,0])
Fth1 = interp1d(thAirL,coordMid)
dydx0 = Fth1(0.0).tolist()
print('----|| ALYA :: PRESSURE SIDE AT X =',dydx0,' dy/dx = 0')
Fth1 = interp1d(coordMid,thAirL)

print('--|| ALYA :: DONE. TIME TAKEN',time.time()-stime)

if 'WIPE' in mode:
 baseDir = os.getcwd()
 if not os.path.exists('0.InsDataSlicesFinal'):
     os.makedirs('0.InsDataSlicesFinal')
 savePath = baseDir+"/0.InsDataSlicesFinal/InsDataSlices.csv"
 f = open(savePath,'w')
 fileList = []
 caseDir = ['/1.t-0-8/','/2.t-8-16/','/3.t-16-24/']
 caseDir = ['/1.t-0-8/']
 caseDir = ['/1.t-0-85/']

 zones = []; times = []; timesShift = [];
 for (count,case) in enumerate(caseDir):
  shft = len(times)
  fileList0 = sorted(glob.glob(baseDir+case+'/0.InsDataSlices/InsData*.csv'))
  fileList.append(fileList0)
  for name in fileList0:
   splitStr = str(name).split(".")
   splitStr2 = str(splitStr[-3]).split("_")
   zones.append(int(splitStr2[-1])) 
   times.append(int(splitStr[-2])) 
   timesShift.append(int(splitStr[-2])+shft)
 fileList = np.transpose(fileList) 

 zones = np.asarray(zones,dtype=int);zones = np.transpose(zones);  
 times = np.asarray(times,dtype=int);times = np.transpose(times);  
 indSort = np.lexsort((zones,times))
 zones = zones[indSort]; times = times[indSort];
 fileList = fileList[indSort];
 L = len(np.unique(zones))
 N = len(np.unique(timesShift))
 x_loc = np.linspace(0.3,0.98,L/2)
 y_loc = np.asarray(np.concatenate((F0(x_loc),F1(x_loc)),axis=None),dtype=float)
 th_loc = np.asarray(np.concatenate((Fth0(x_loc),Fth1(x_loc)),axis=None),dtype=float)
 x_loc = np.asarray(np.concatenate((x_loc,x_loc),axis=None),dtype=float)

 pltInd = [0,2,4,5,7,9]

 fname = (baseDir+case+'/0.InsDataSlices/InsDataSlice_%s.%s.csv') % (0,0)
 data=np.around(np.loadtxt(fname,dtype=float,comments='#',delimiter=',',skiprows=1),decimals=6)
 z1 = data[:,-1];
 y1 = data[:,-2];
 x1 = data[:,-3];
 r1 = np.sqrt((x1-0.3)**2+(y1-F0(0.3))**2);
 xI = np.sort(np.unique(r1)); nr = len(xI);
 rI = np.sort(np.unique(r1));
 yI = np.sort(np.unique(z1)); nz = len(yI)
 f.write("# %s \n" % '  '.join(str(elm) for elm in xI));
 yI,xI = np.meshgrid(yI,xI);
 xI = np.reshape(xI,(np.size(xI),1)); yI = np.reshape(yI,(np.size(yI),1));

 bfEvents = np.empty((0,nr),dtype=float);
 #bfEvents = [];
 for ii in range(L):
  x0 = x_loc[ii]; y0 = y_loc[ii];
  t0 = th_loc[ii]
  m = np.tan(th_loc[ii]);
  rhs = d**2/(1+float(1.0/m**2));
  if(y0<0):
    if(m>0.0):
      x1 = x0 - np.sqrt(rhs);
    else:
      x1 = x0 + np.sqrt(rhs);
    y1 = y0-(x1-x0)/m;
  elif(y0>=0):
    if(m<0.0):
      x1 = x0 - np.sqrt(rhs);
    else:
      x1 = x0 + np.sqrt(rhs);
    y1 = y0-(x1-x0)/m;
  uArr = np.empty((0,nr),dtype=float)
  vArr = np.empty((0,nr),dtype=float)
  wArr = np.empty((0,nr),dtype=float)
  pArr = np.empty((0,nr),dtype=float)
  f.write("# ZONE %s \n" % (ii));
  f.write("# %s %s \n" % (x0,y0));

  print('--|| ALYA :: WORKING ON ZONE ',ii)
  stime = time.time()
  for jj in range(N):
   ### READ DATA AND INTERPOLATE ON A STRUCTURED GRID ###
   fname = (baseDir+case+'/0.InsDataSlices/InsDataSlice_%s.%s.csv') % (ii,jj)
   data=np.around(np.loadtxt(fname,dtype=float,comments='#',delimiter=',',skiprows=1),decimals=6)
   z1 = data[:,-1]; y1 = data[:,-2]; x1 = data[:,-3];
   r1 = np.sqrt((x1-x0)**2+(y1-y0)**2);
   
   ### Interpolate Ut
   zO = griddata((r1,z1),(data[:,0]*np.cos(t0)+data[:,1]*np.sin(t0)), (xI,yI), method='linear')
   zO = np.reshape(zO,(nz,nr))
   uArr = np.vstack((uArr,zO))
   ### Interpolate Vn
   zO = griddata((r1,z1),(-data[:,0]*np.sin(t0)+data[:,1]*np.cos(t0)), (xI,yI), method='linear')
   zO = np.reshape(zO,(nz,nr))
   vArr = np.vstack((vArr,zO))
   ### Interpolate W
   zO = griddata((r1,z1),data[:,2], (xI,yI), method='linear')
   zO = np.reshape(zO,(nz,nr))
   wArr = np.vstack((wArr,zO))
   ### Interpolate P
   zO = griddata((r1,z1),data[:,3], (xI,yI), method='linear')
   zO = np.reshape(zO,(nz,nr))
   pArr = np.vstack((pArr,zO))

  print('--|| ALYA :: DONE. TIME TAKEN',time.time()-stime)
  num_rows, num_cols = uArr.shape
  uArr = np.transpose(uArr)
  vArr = np.transpose(vArr)
  wArr = np.transpose(wArr)
  pArr = np.transpose(pArr)
  for wr in range(nr):
   f.write('#%s \n' % wr);
   f.write("%s \n" % '  '.join(str(elm) for elm in uArr[wr,:]));
   f.write("%s \n" % '  '.join(str(elm) for elm in vArr[wr,:]));
   f.write("%s \n" % '  '.join(str(elm) for elm in wArr[wr,:]));
   f.write("%s \n" % '  '.join(str(elm) for elm in pArr[wr,:]));
  ## BACKFLOW EVENTS CALC##
  bfEvents = np.vstack((bfEvents,np.asarray(np.sum(uArr > 0, axis=1),dtype=float)/num_rows))
  #print(np.shape(np.transpose(rI)),np.shape(np.sum(uArr < 0, axis=1)/len(uArr)))
  #bfEvents.append([[np.transpose(rI)],[np.sum(uArr < 0, axis=1)/len(uArr)]])
 f.close()
 ## BACKFLOW EVENTS ##
 print('--|| ALYA :: PLOTTING RESULTS.')
 plt_row,plt_clm = (2,int(len(pltInd)/2))
 fig,axs = plt.subplots(plt_row,plt_clm)
 count=0
 for row in range(plt_row):
  for clmn in range(plt_clm):
   if(count in pltInd):
    if(y_loc[ii] < 0.0):
     l1, = axs[row,clmn].plot(bfEvents[ii,:],rI,linewidth=lw)
    else:
     l1, = axs[row,clmn].plot(bfEvents[ii,:],rI,linewidth=lw)
   count=count+1  
 for ax in axs.flat:
     ax.set(xlabel=r'$x/c$', ylabel=r'$r/c$')
 # Hide x labels and tick labels for top plots and y ticks for right plots.
 for ax in axs.flat:
     ax.label_outer()
 
 savePath = casePath+'/backFlow-'
 savePath = savePath+airfoil
 savePath = savePath+'.png'
 plt.savefig(savePath, format='png',\
             dpi=600, facecolor='w', edgecolor='w',\
 	    orientation='portrait',transparent=True,\
 	    bbox_inches=None,pad_inches=0.1,\
 	    papertype='a4',frameon=None, metadata=None)
 ## QUADRANT ANALYSIS ##
 ## PDF ANALYSIS ##

elif 'DUMP' in mode:
 x_loc = np.linspace(0.3,0.98,5)
 
 if not os.path.exists('0.InsDataSlices'):
     os.makedirs('0.InsDataSlices')
 savePath = casePath+"/0.InsDataSlices/InsDataSlice_.csv" 
 
 # INITIALIZE VARIABLES
 fileName = caseName+'.ensi.case' 
 
 print('--|| ALYA :: READING FILE ',fileName)
 stime = time.time()
 case = OpenDataFile(fileName) 
 case.PointArrays = ['VELOC', 'PRESS'] 
 case.UpdatePipeline() 
 print('--|| ALYA :: DONE. TIME TAKEN',time.time()-stime)
 
 reader = GetActiveSource()
 view = GetActiveView()
 times = reader.TimestepValues
 
 nSlice = True;
 
 view = GetActiveView()
 times = reader.TimestepValues
 
 nSlice = True;
 
 
 
 
 #
 # EXTRACT LOCATIONS
 listPlanes = []
 #listPlanes = np.empty((0,2*len(x_loc)+1))
 
 ES = ExtractSurface(Input=case)
 ES.UpdatePipeline()
 clip = Clip(Input=ES)
 clip.ClipType='Box'
 clip.ClipType.Scale = [0.1,0.1,1]
 clip.UpdatePipeline()
 #listPlanes.append(clip)
 
 count=0
 for side in sides:
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
  
  print('--|| ALYA :: AIRFOIL VARIABLES LOADED', side)
  
  ########### GENERATE SLICES ###################### 
  if nSlice:
    print('--|| ALYA EXTRACTS :: ', x_loc)
    for n in range(0,len(x_loc)):
      # Generate Normal Slices on Airfoil Surface
      slice1 = Slice(Input=case)
      strName = '1.slice-%s' % (n)
      RenameSource(strName, slice1) 
      slice1.SliceType.Origin = [center[n,0], center[n,1], 0.1]
      slice1.SliceType.Normal = [center[n,2], center[n,3], 0.0]
      slice1.UpdatePipeline()
      ## Generate Surface Normals
      #slice1 = GenerateSurfaceNormals(Input=slice1)
      #slice1.UpdatePipeline()
      # Clip Surface 1
      clip = Clip(Input=slice1)
      # rename source object
      strName = 'sliceClip'
      RenameSource(strName, clip) 
      # Properties modified
      clip.ClipType.Origin = [center[n,4], center[n,5], 0.1]
      clip.ClipType.Normal = [-center[n,3], center[n,2], 0.0]
      if('SS' in side):
       clip.Invert = 1
      elif('PS' in side):
       clip.Invert = 0
      clip.UpdatePipeline()
      # Clip Surface 2
      clip = Clip(Input=clip)
      # rename source object
      strName = '2.sliceClip-%s' % (count)
      RenameSource(strName, clip) 
      # Properties modified
      clip.ClipType.Origin = [center[n,0], center[n,1], 0.1]
      clip.ClipType.Normal = [-center[n,3], center[n,2], 0.0]
      if('SS' in side):
       clip.Invert = 0
      elif('PS' in side):
       clip.Invert = 1
      clip.UpdatePipeline()
 #     print('--|| ALYA :: GROUPING TIMESET ', n)
 #     stime = time.time()
 #     dataSet = GroupTimeSteps(Input=clip)
 #     dataSet.UpdatePipeline()
 #     print('--|| ALYA :: DONE. TIME TAKEN ', time.time()-stime)
 #     print('--|| ALYA :: APPLY PROGRAMMABLE FILTER ', n)
 #     stime = time.time()
 #     PF1 = ProgrammableFilter(Input=dataSet) 
 #     PF1.Script = \
 # """ 
 # print("----|| ALYA :: MANUPULATE DATA",n) 
 # P = dsa.WrapDataObject(inputs[0].GetBlock(0)).PointData['PRESS'].Arrays[0] 
 # P = vtk_to_numpy(P) 
 # N = np.size(times) 
 # L = np.size(P) 
 # print('----|| ALYA : TEMPORAL SIZE IS ',N, 'SPATIAL SIZE IS ',L) 
 # field1 = np.zeros([N,L])
 # field2 = np.zeros([N,L,3])  
 # print ('----|| ALYA : READING TIME: FROM {} TO {}'.format(times[0],times[-1])) 
 # for i in range(N): 
 #   t = times[i] 
 #   d = dsa.WrapDataObject(inputs[0].GetBlock(i)) 
 #   d = d.PointData["PRESS"].Arrays[0] 
 #   field = vtk_to_numpy(d) 
 # field1_avg = np.average(field1, axis=0) 
 # field2_avg = np.average(field1, axis=0) 
 #
 # prms = field1 - field1_avg 
 # urms = field1 - field1_avg 
 #     """ 
 #     PF1.UpdatePipeline()
 #     print('--|| ALYA :: DONE. TIME TAKEN ', time.time()-stime)
      listPlanes.append(clip)
      count = count+1
 
 print('--|| ALYA :: GROUPING DATASETS.')
 stime=time.time()
 DS = GroupDatasets(Input = listPlanes)
 DS.UpdatePipeline()
 print('--|| ALYA :: DONE. TIME TAKEN',time.time()-stime)
 print('--|| ALYA :: SAVING DATASETS.')
 stime=time.time()
 if('.vtm' in savePath):
  SaveData(savePath, proxy=DS, Writetimestepsasfileseries=1)
 elif('.csv' in savePath): 
  SaveData(savePath, proxy=DS, WriteTimeSteps=1)
 else:
   raise ValueError('--|| ALYA ERROR :: SAVE OR NOT?')
 print('--|| ALYA :: DONE. TIME TAKEN',time.time()-stime)
else: 
 raise ValueError('--|| ALYA ERROR :: WIPE OR DUMP?')




