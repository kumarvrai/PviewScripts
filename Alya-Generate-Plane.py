import os
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import matplotlib.tri as tri
#from dolointerpolation import MultilinearInterpolator
import csv
#from scipy.interpolate import interp2d
niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'
niaHome = '/home/kvishal/projects/rrg-ugo-ac/kvishal/0.pre_process/1.pre_process_files/'

# LOAD AIRFOIL SPECIFIC FILES
airfoil = '4412'; 
side = 'SS';
surf='tang';
#dVec = [0.001,0.002,0.004,0.006]
#dVec = [0.001,0.002,0.0035,0.004,0.005,0.006,0.008,0.01]
#mode = 'yplus';
dVec = [5.5,11,30] #YPLUS VALUES

mode = 'yplus';
#mode = 'd99';
#dVec = [0.1] #DELTA_99 VALUES

print('--|| ALYA INITIALIZING')

if('0012' in airfoil):
  fAirU = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-UP.txt'
  fAirL = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/3.sphereRough-NACA0012/naca0012-DOWN.txt'
  nu = 2.0e-5
elif('4412' in airfoil):  
  fAirU = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-UP.txt'
  fAirL = niaHome+'/1.GridGen/b.Airfoil/3.RoughAirfoil/4.sphereRough-NACA4412/naca4412-DOWN.txt'
  #-------HEADER----------------#
  #  'x', 'Ue', 'Ve', 'u_t', 'd_99', 'dStar', 'dTheta', 'H', 'beta', 'Re_t', 'Re_th' 
  #fBLcase ='./2.4412-Slices/1.SLE/1.SLE-5/BL_data_SLE-5-' 
  #fBLcase ='./2.4412-Slices/1.SLE/2.SLE-10/BL_data_SLE-10-' 
  #fBLcase ='./2.4412-Slices/1.SLE/3.SLE-15/BL_data_SLE-15-'

  #fBLcase ='./2.4412-Slices/2.TBL/1.TBL-5/2.Medium/BL_data_TBL-5-' 
  #fBLcase ='./2.4412-Slices/2.TBL/1.TBL-5/3.Fine/BL_data_TBL-5-' 

  #fBLcase ='./2.4412-Slices/3.RLE1/1.AoA-5/BL_data_RLE1-5-' 
  #fBLcase ='./2.4412-Slices/3.RLE1/2.AoA-10/BL_data_RLE1-10-' 
  #fBLcase ='./2.4412-Slices/3.RLE1/3.AoA-15/BL_data_RLE1-15-' 

  #fBLcase ='./2.4412-Slices/4.RLE2/1.AoA-5/BL_data_RLE2-5-' 
  #fBLcase ='./2.4412-Slices/4.RLE2/2.AoA-10/BL_data_RLE2-10-' 
  fBLcase ='./2.4412-Slices/4.RLE2/3.AoA-15/BL_data_RLE2-15-' 
  fBLcase = fBLcase+side+'.csv' 
  nu = 5.0e-6

print('--|| ALYA :: READING AIRFOIL DATA.')
stime = time.time()
coordAirU = np.loadtxt(fAirU, delimiter=',')
coordAirL = np.loadtxt(fAirL, delimiter=',')

thAirU = np.arctan2(np.diff(coordAirU[:,1]),np.diff(coordAirU[:,0]))
F0 = interp1d(coordAirU[:,0],coordAirU[:,1])
airLen = len(coordAirU)
coordMid = 0.5*(coordAirU[0:airLen-1,0]+coordAirU[1:airLen,0])
#print(airLen,len(thAirU),len(coordMid))
Fth0 = interp1d(coordMid,thAirU)

thAirL = np.arctan2(np.diff(coordAirL[:,1]),np.diff(coordAirL[:,0]))
F1 = interp1d(coordAirL[:,0],coordAirL[:,1])
airLen = len(coordAirL)
coordMid = 0.5*(coordAirL[0:airLen-1,0]+coordAirL[1:airLen,0])
#print(airLen,len(thAirL),len(coordMid))
Fth1 = interp1d(coordMid,thAirL)

blData = np.loadtxt(fBLcase, delimiter=',',skiprows=1)

F_Ue = interp1d(blData[:,0],blData[:,1])
F_d99 = interp1d(blData[:,0],blData[:,4])
F_Ke = interp1d(blData[:,0],blData[:,11])

if('yplus' in mode):
  F_intrp = interp1d(blData[:,0],blData[:,3])
elif('d99' in mode):
  F_intrp = interp1d(blData[:,0],blData[:,4])

if 'tang' in surf:
  # Define Arrays for coordinates
  #npts = 1001; npts_z = len(np.unique(zz));
  npts = 1001; npts_z = 145;
  print('--||ALYA UNIQUE Z =',npts_z);
  xlinV = np.zeros((npts,len(dVec)), dtype=float)
  ylinV = np.zeros((npts,len(dVec)), dtype=float)
  xAir  = np.zeros((npts*npts_z,1), dtype=float)
  yAir  = np.zeros((npts*npts_z,1), dtype=float)
  zAir  = np.zeros((npts*npts_z,1), dtype=float)
  sAir  = np.zeros((npts*npts_z,1), dtype=float)
  UeAir  = np.zeros((npts*npts_z,1), dtype=float)
  d99Air  = np.zeros((npts*npts_z,1), dtype=float)
  KeAir  = np.zeros((npts*npts_z,1), dtype=float)
  tVecX  = np.zeros((npts*npts_z,1), dtype=float)
  tVecY  = np.zeros((npts*npts_z,1), dtype=float)
  nVecX  = np.zeros((npts*npts_z,1), dtype=float)
  nVecY  = np.zeros((npts*npts_z,1), dtype=float)
  xAirV = np.zeros((npts*npts_z,len(dVec)), dtype=float)
  yAirV = np.zeros((npts*npts_z,len(dVec)), dtype=float)
  wAirV = np.zeros((npts*npts_z,len(dVec)), dtype=float)
  stime = time.time()
  
  
  for ik in range(len(dVec)):
    print('--|| ALYA :: CREATING A SURFACE FOR INTERPOLATION FOR D=',dVec[ik])
    d = dVec[ik]
    xlin = np.linspace(0.0,0.99,npts);
    zlin = np.linspace(0,0.2,npts_z);
    if('SS' in side):
      ylin = F0(xlin)
      th_loc = Fth0(xlin)
    elif('PS' in side):
      ylin = F1(xlin)
      th_loc = Fth1(xlin)
    slin = 0.0*xlin
    Uelin = 0.0*xlin
    d99lin = 0.0*xlin
    Kelin = 0.0*xlin
    tVX = 0.0*xlin; tVY = 0.0*xlin 
    nVX = 0.0*xlin; nVY = 0.0*xlin 
    for i in range(len(xlin)):
     x0 = xlin[i];y0 = ylin[i];
     m = np.tan(th_loc[i]);
     Uelin[i] = F_Ue(x0)
     d99lin[i] = F_d99(x0)
     Kelin[i] = F_Ke(x0)
     if('yplus' in mode):
      utau = F_intrp(x0)
      d = dVec[ik]*nu/utau;
     elif('d99' in mode):
      d = dVec[ik]*F_intrp(x0);

     if('SS' in side):
       if(m>0.0):
         rhs = d**2/(1+float(1.0/m**2));
         x1 = x0 - np.sqrt(rhs);
         y1 = y0-(x1-x0)/m;
       elif(m<0.0):
         rhs = d**2/(1+float(1.0/m**2));
         x1 = x0 + np.sqrt(rhs);
         y1 = y0-(x1-x0)/m;
       else:
         x1 = x0;
         y1 = y0+d;
     elif('PS' in side):
       if(m<0.0):
         rhs = d**2/(1+float(1.0/m**2));
         x1 = x0 - np.sqrt(rhs);
         y1 = y0-(x1-x0)/m;
       elif(m>0.0):
         rhs = d**2/(1+float(1.0/m**2));
         x1 = x0 + np.sqrt(rhs);
         y1 = y0-(x1-x0)/m;
       else:
         x1 = x0;
         y1 = y0-d;
     xlin[i] = x1; ylin[i] = y1
     if(i>0):
      #if('ynorm' in mode):
      # slin[i] = np.trapz(np.sqrt(1+np.tan(th_loc[0:i])**2),xlin[0:i]);
      #elif('yplus' in mode):
      slin[i] = slin[i-1]+np.sqrt((xlin[i-1]-x1)**2+(ylin[i-1]-y1)**2);
     tVX[i] = np.cos(th_loc[i]); tVY[i] = np.sin(th_loc[i]);
     nVX[i] = -np.sin(th_loc[i]); nVY[i] = np.cos(th_loc[i]);
  
    print('----|| ALYA :: BOX DIM X',min(xlin),max(xlin))
    print('----|| ALYA :: BOX DIM Y',min(ylin),max(ylin))
    print('----|| ALYA :: BOX DIM Z',min(zlin),max(zlin))

    count = 0;
    for k in range(npts_z):
       for i in range(npts):
        xAir[count] = xlin[i];
        yAir[count] = ylin[i];
        zAir[count] = zlin[k];
        sAir[count] = slin[i];
        tVecX[count] = tVX[i];
        tVecY[count] = tVY[i];
        nVecX[count] = nVX[i];
        nVecY[count] = nVY[i];
        UeAir[count] = Uelin[i];
        d99Air[count] = d99lin[i];
        KeAir[count] = Kelin[i];
        count=count+1;

    print('--|| ALYA :: DONE. TIME=', time.time()-stime)
  
    xAir = np.ravel(xAir)
    yAir = np.ravel(yAir)
    zAir = np.ravel(zAir)
    sAir = np.ravel(sAir)
    tVecX = np.ravel(tVecX)
    tVecY = np.ravel(tVecY)
    nVecX = np.ravel(nVecX)
    nVecY = np.ravel(nVecY)
    UeAir = np.ravel(UeAir)
    d99Air = np.ravel(d99Air)
    KeAir = np.ravel(KeAir)
    
    #### Write CSV Files
    print('--|| ALYA :: WRITING CSV FILES')
    locStr = mode+'Plane-'+side+'-'
    locStr = locStr+np.array2string(np.asarray(dVec[ik], dtype=float),separator='-')
    with open(locStr+'.csv','w') as csvfile:
     spamwriter = csv.writer(csvfile, delimiter=',')
     spamwriter.writerow(['x','y','z','s','tVX','tVY','nVX','nVY','Ue','d99','Ke'])
     spamwriter.writerows(np.transpose([xAir,yAir,zAir,sAir,tVecX,tVecY,nVecX,nVecY,UeAir,d99Air,KeAir]))
    print('--|| ALYA :: DONE. TIME=', time.time()-stime)

exit()
print('--|| ALYA :: PLOTTING RESULTS')
stime = time.time()
rows,cols = (len(dVec)), 1
fig, ax1 = plt.subplots(nrows=rows, ncols=cols);
levels = np.linspace(-0.05,0.05,256);
ind = 0
if(rows==1 or cols==1):
 if(rows != 1):
  count = rows
 else:
  count = cols
 for i in range(count):
   cf = ax1[i].tricontourf(xAirV[:,ind], yAirV[:,ind], wAirV[:,ind], levels, cmap='RdBu',extend='both')
   ind = ind + 1
else:
 for i in range(rows):
   for j in range(cols):
       cf = ax1[i,j].tricontourf(xAirV[:,ind], yAirV[:,ind], wAirV[:,ind], levels, cmap='RdBu',extend='both')
       ind = ind + 1
    #ax1[1].triplot(triang, 'ko-')

for ax in ax1.flat:
    ax.set(xlabel='s/c')

fig.text(0.02, 0.55, 'z/c', va='center', fontsize=MEDIUM_SIZE, rotation='vertical')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in ax1.flat:
    ax.label_outer()
plt.tight_layout();
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.215, 0.035, 0.7])
cbar = fig.colorbar(cf, cax=cbar_ax,extendrect=True);
cbar.set_ticks([min(levels),0,max(levels)])
cbar.set_ticklabels([min(levels),0,max(levels)])
#fig.text(0.5, 0.03, '$s/c$', ha='center')
#fig.text(0.04, 0.5, '$z/c$', va='center', rotation='vertical')
savePath = './normContour-'
varStr = var.replace(':', '-')
locStr = np.array2string(np.asarray(dVec, dtype=float),separator='-')
savePath = savePath+varStr
savePath = savePath+locStr
savePath = savePath+'.png'
plt.savefig(savePath, format='png',\
            dpi=600, facecolor='w', edgecolor='w',\
	    orientation='portrait',transparent=True,\
	    bbox_inches=None, pad_inches=0.1,\
	    frameon=None, metadata=None)

fig1, ax2 = plt.subplots();
ax2.plot(coordAirU[:,0],coordAirU[:,1],'k',linewidth=1)
ax2.plot(coordAirL[:,0],coordAirL[:,1],'k',linewidth=1)
for k in range(len(dVec)):
 ax2.plot(xlinV[:,k],ylinV[:,k],linewidth=1,markersize=0)
plt.xlabel('x/c')
plt.ylabel('y/c')
frame1 = plt.gca()
#frame1.axes.get_yaxis().set_ticks([])
ax2.set_aspect(2);
savePath = './normLoc-'
savePath = savePath+airfoil
savePath = savePath+'.png'
plt.savefig(savePath, format='png',\
            dpi=600, facecolor='w', edgecolor='w',\
	    orientation='portrait',transparent=True,\
	    bbox_inches=None, pad_inches=0.1,\
	    frameon=None, metadata=None)
#plt.show();

