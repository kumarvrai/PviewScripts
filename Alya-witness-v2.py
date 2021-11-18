import os
import time
import sys

os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import matplotlib.pyplot as plt
import numpy as np


ls = ['-','--','-.',':','.','-,','-o','-v','-^','-s','-p','-d','-*','-+','-x']
caseName	= sys.argv[1]
airfoil 	= sys.argv[2]

# LOAD AIRFOIL SPECIFIC FILES
niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'
niaHome = '/home/u/ugo/kvishal/'


print('--|| ALYA INITIALIZING')

if('0012' in airfoil):
  fileLoc = niaScrh+'/5.BGQ-DATA/NACA-0012/2.3D/2.LES/1.Re-5E4-alpha-8/3.WALE/1.run1/1.t-0-13/'
  fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-UP.txt'
  fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-DOWN.txt'
elif('4412' in airfoil):  
  fileLoc = niaScrh+'/2.NACA4412/2.ROUGH_CASE/1.trip/2.run2/2.t-2-12/'
  fAirU = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-UP.txt'
  fAirL = niaHome+'/0.PreProc/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-DOWN.txt'
else:
  sys.exit("--||Alya (ERROR): AIRFOIL NOT IN THE LIST")

print('--|| ALYA :: READING AIRFOIL DATA.')
stime = time.time()
coordAirU = np.loadtxt(fAirU, delimiter=',')
coordAirL = np.loadtxt(fAirL, delimiter=',')
print('--|| ALYA :: DONE. TIME=',time.time()-stime)

print('--|| ALYA :: READING WITNESS POINTS.')
stime = time.time()
witPoints = np.loadtxt('./witness-nastin.dat')
print('--|| ALYA :: DONE. TIME=',time.time()-stime)

nWit = len(witPoints)
z = np.unique(witPoints[:,2]);
nz = nWit/len(z)
#indP = range(nz); #Indices of points
indP = [0,1,2,3,4,5]
print('--|| ALYA :: Z LOCATIONS ',z)
print('--|| ALYA :: PROCESSING', nWit, 'WITNESS POINTS')
fname = caseName+'.nsi.wit'

data = np.empty((0,5),dtype='float')
for line in open(fname):
 li=line.strip()
 if not li.startswith("#"):
  li = np.array(map(float, line.split()));
  data = np.vstack((data,li))

print('--|| ALYA :: CALCULATING WIT POINT CORRELATIONS.')
stime = time.time()
uu = np.zeros((len(z),nz),dtype=float)
uv = np.zeros((len(z),nz),dtype=float)
pp = np.zeros((len(z),nz),dtype=float)
up = np.zeros((len(z),nz),dtype=float)

for j in indP:
 u0 = data[np.where(data[:,0]==(j+1)),1]; u0 = u0-np.mean(u0,axis=1);
 v0 = data[np.where(data[:,0]==(j+1)),2]; v0 = v0-np.mean(v0,axis=1);
 ip = 3
 p0 = data[np.where(data[:,0]==(j+1)),ip]; p0 = p0-np.mean(p0,axis=1);
 for i in range(len(z)):
  pair = i*nz+j+1;
  #print(j+1,pair);
  # Calculate uu
  u1 = data[np.where(data[:,0]==pair),1]; u1 = u1-np.mean(u1,axis=1);
  prod = u0*u1
  uu[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(u0*u0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
  # Calculate uv
  u1 = data[np.where(data[:,0]==pair),2]; u1 = u1-np.mean(u1,axis=1);
  prod = u0*u1
  uv[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(u0*u0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
  # Calculate pp
  u1 = data[np.where(data[:,0]==pair),ip]; u1 = u1-np.mean(u1,axis=1);
  prod = p0*u1
  pp[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(p0*p0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
  # Calculate up
  u1 = data[np.where(data[:,0]==pair),ip]; u1 = u1-np.mean(u1,axis=1);
  prod = u0*u1
  up[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(u0*u0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
print('--|| ALYA :: DONE. TIME=',time.time()-stime)
