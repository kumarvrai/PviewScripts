import numpy as np
from scipy.interpolate import interp1d

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    indMin = np.argmin(np.sqrt(dist_2),axis=None)
    return indMin, np.amin(np.sqrt(dist_2),axis=None)

airfoil = '4412';
side = 'ss';

# LOAD AIRFOIL SPECIFIC FILES
print('--|| ALYA INITIALIZING')
baseDir = '/home/kvishal/1.post_process/1.airfoil/3.PviewScripts/'
if('0012' in airfoil):
  fAirU = baseDir+'/1.0012-Slices/naca0012-UP.txt'
  fAirL = baseDir+'/1.0012-Slices/naca0012-DOWN.txt'
  sliceLoc = baseDir+'/1.0012-Slices/'
elif('4412' in airfoil):
  fAirU = baseDir+'/2.4412-Slices/naca4412-UP.txt'
  fAirL = baseDir+'/2.4412-Slices/naca4412-DOWN.txt'
  sliceLoc = baseDir+'/2.4412-Slices/'
else:
  raise ValueError('--|| ALYA ERROR :: FILE NOT PROVIDED.')

if('ss' in side):
  coordAir = np.loadtxt(fAirU, delimiter=',');
else:
  coordAir = np.loadtxt(fAirL, delimiter=',')

# INTERPOLATION ON SUCTION SIDE
thAir = np.arctan2(np.diff(coordAir[:,1]),np.diff(coordAir[:,0]))
airLen = len(coordAir)
coordMid = 0.5*(coordAir[0:airLen-1,0]+coordAir[1:airLen,0])
Fth = interp1d(coordMid,thAir)


d = dsa.WrapDataObject(inputs[0].GetBlock(0))
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
pdist = np.zeros(np.shape(x));
tvect = np.zeros((np.shape(x)[0],2));
nvect = np.zeros((np.shape(x)[0],2));

for i in range(len(x)):
  node = (x[i],y[i]);
  indMin,pdist[i] = closest_node(node, coordAir);
  theta = Fth(x[indMin])
  tvect[i,:] = (cos(theta), sin(theta))
  nvect[i,:] = (-sin(theta), cos(theta))

pdist=pdist-np.amin(pdist,axis=None)
output.ShallowCopy(inputs[0].VTKObject)  
output.PointData.append(pdist,"PDIST")
output.PointData.append(tvect,"TVECT")
output.PointData.append(nvect,"NVECT")  
