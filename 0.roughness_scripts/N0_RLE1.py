"""
Created on Sun Feb 17 17:53:23 2019

@author: dpastrana
"""
import time
import os
import sys
import math as m
from mpl_toolkits import mplot3d 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plt
from matplotlib import cm
import numpy as np
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import random as rand
from scipy.stats import norm
from scipy.stats import skewnorm

plt.close('all')


def unwrap_cyl(x, y, xc0, yc0):
    r = m.sqrt((x-xc0)**2 + (y-yc0)**2)
    theta = m.atan2(y-yc0, x-xc0)
    return r, theta


def wrap_cyl(r, theta,xc0, yc0):
    x = r * m.cos(theta) + xc0
    y = r * m.sin(theta) + yc0
    return x, y

def genSphereLoc(xMin,xMax,zMin,zMax,delx,delz,h):
    airCoordU = RGH_SCRPTS+'naca0012-UP.txt'
    airCoordD = RGH_SCRPTS+'naca0012-DOWN.txt'
    aircup = np.loadtxt(airCoordU)
    aircdo = np.loadtxt(airCoordD)
    # Upper Surface
    dsA = np.sqrt(1.0 + np.array(np.diff(aircup[:,1])/np.diff(aircup[:,0])))
    ds0 = np.zeros((len(dsA)+1,1),dtype='float')
    xds = np.zeros((len(dsA)+1,1),dtype='float')
    for i in range(1,len(ds0)):
       xds[i] = aircup[i,0]
       ds0[i] = np.trapz(dsA[0:i],x=np.array(aircup[0:i,0]))
       #print(xds[i],ds0[i])
    Fds = interp1d(ds0[:,0],xds[:,0])
    nx = np.floor_divide((xMax-xMin),delx) + 1
    xLocU = Fds(np.array([0.5*h, 2.5*h]))
    xLocU = np.append(xLocU,Fds(np.linspace(5.0*h,xMax,nx)))
    # Lower Surface
    dsA = np.sqrt(abs(1.0 + np.array(np.diff(aircdo[:,1])/np.diff(aircdo[:,0]))))
    ds0 = np.zeros((len(dsA)+1,1),dtype='float')
    xds = np.zeros((len(dsA)+1,1),dtype='float')
    for i in range(1,len(ds0)):
       xds[i] = aircdo[i,0]
       ds0[i] = np.trapz(dsA[0:i],x=np.array(aircdo[0:i,0]))
    Fds = interp1d(-ds0[:,0],xds[:,0])
    nx = np.floor_divide((xMax-xMin),delx) + 1
    xLocD = Fds(-np.array([0.1*h, 1.5*h]))
    xLocD = np.append(xLocD,Fds(-np.linspace(5.0*h,xMax,nx)))
    # z-direction
    nz = np.floor_divide((zMax-zMin),delz) + 1
    zLoc = np.linspace(zMin,zMax,nz)
    return xLocU, xLocD, zLoc

fileType = sys.argv[1]
RGH_SCRPTS='/home/kvishal/1.post_process/0.alya_pv_scripts/0.roughness_scripts/input_files/'
###########################################
if fileType == 'ALYA':
   geofile_in = './naca_old.msh'
   geofile_out= './naca.msh'
   bString = '$Nodes'
   eString = '$EndNodes'
elif fileType == 'ALYA':
   geofile_in = './naca.geo.dat'
   geofile_out= './naca_r.geo.dat'
   bString = 'COORDINATES'
   eString = 'END_COORDINATES'
airfoilCoord = RGH_SCRPTS+'naca0012.txt'
expDistRough = RGH_SCRPTS+'expDistRough.txt'

f = open(geofile_in, 'r')
g = open(geofile_out, 'w')
db = open('./out.o', 'w')

nNodes = int(os.popen('sed -n ''/%s/,/%s/p'' %s | wc -l ' % (bString,eString,geofile_in)).read())
print '--||ALYA :: WORKING WITH %i NODES IN TOTAL' % nNodes

lineB = True
inCoords = False

#Load Files for generating interpolants
print('--||ALYA :: READING FILES FOR INTERPOLATION')
startTime = time.time()

coordAirfoil = np.loadtxt(airfoilCoord)
expDist = np.loadtxt(expDistRough)

tanAirfoil = np.arctan2(np.diff(coordAirfoil[:,1]),np.diff(coordAirfoil[:,0]))
airLen = np.floor_divide(len(coordAirfoil),2)
F0 = interp1d(coordAirfoil[0:airLen,0],coordAirfoil[0:airLen,1])
F1 = interp1d(coordAirfoil[airLen:2*airLen-1,0],coordAirfoil[airLen:2*airLen-1,1])
Fth0 = interp1d(coordAirfoil[0:airLen,0],tanAirfoil[0:airLen])
Fth1 = interp1d(coordAirfoil[airLen:2*airLen-1,0],tanAirfoil[airLen:2*airLen-1])
F  = interp2d(coordAirfoil[0:2*airLen-1,0],coordAirfoil[0:2*airLen-1,1],tanAirfoil)
xAir,zAir = np.meshgrid(coordAirfoil[:,0],np.linspace(0,0.2,10))
yAir,zAir = np.meshgrid(coordAirfoil[:,1],np.linspace(0,0.2,10))
Fexp = interp1d(expDist[:,0],expDist[:,1])
print('----||DONE. TIME =',time.time()-startTime,'sec')

#Define some Parameters
zdom = 0.2
r0 = 0.5
D = 0.08
h = 0.008
H = 1.0


fact = 0.0
b_radio = 10.0*h
center_y =  (h**2.0 + (0.5*D)**2.0)/(2.0*h) - h
rad = center_y+h

# Generate roughness locations locations (x,y)
vecLoc = genSphereLoc(0.0, 0.2, 0.0, 0.2, 5.0*h, 5.0*h, h)
XposU = vecLoc[0]
XposD = vecLoc[1]
factRoughU = Fexp(XposU)
factRoughU = factRoughU/np.amax(factRoughU)
factRoughD = Fexp(XposD)
factRoughD = factRoughD/np.amax(factRoughD)
    
#Create center, radius vectors
center = np.empty((0,3),dtype='float')
rS = []
# Roughness on Upper surface
for i in range(0,len(XposU)):
  if(XposU[i]<0.01 or XposU[i]>=0.1):
     Zpos = np.linspace(0.0,0.2,rand.randint(4,8))
  else:
     Zpos = np.linspace(0.0,0.2,rand.randint(8,12))
  for k in range(0,len(Zpos)):
    if(Zpos[k]==0 or Zpos[k]==0.2):
       xLoc = XposU[i]
       zLoc = Zpos[k]
    else:
       if(XposU[i]>h):
          xLoc = XposU[i]+ h*rand.uniform(-0.5,0.5)
       else:
          xLoc = XposU[i]
       zLoc = Zpos[k] + h*rand.uniform(-0.5,0.5) 
    yLoc = F0(xLoc)
    center = np.vstack((center,np.array([xLoc,yLoc,zLoc])))
    rS.append(factRoughU[i]*h)

# Roughness on Lower surface
for i in range(0,len(XposD)):
  if(XposD[i]<0.01 or XposD[i]>=0.1):
     Zpos = np.linspace(0.0,0.2,rand.randint(4,8))
  else:
     Zpos = np.linspace(0.0,0.2,rand.randint(8,12))
  for k in range(0,len(Zpos)):
    if(Zpos[k]==0 or Zpos[k]==0.2):
       xLoc = XposD[i]
       zLoc = Zpos[k]
    else:
       if(XposD[i]>h):
          xLoc = XposD[i]+ h*rand.uniform(-0.5,0.5)
       else:
          xLoc = XposD[i]
       zLoc = Zpos[k] + h*rand.uniform(-0.5,0.5) 
    yLoc = F1(xLoc)
    center = np.vstack((center,np.array([xLoc,yLoc,zLoc])))
    rS.append(factRoughD[i]*h)

#Create distance vectors for smoothing
Roughness = len(center)
aS = []
bS = []
cS = []
thetaS = []
thetaAirVec = []
for n in range(0,Roughness):
  aS.append(3.0*rS[n]) #smoothing distance in x-direction
  bS.append(3.0*rS[n]) #smoothing distance in y-direction
  cS.append(2.0*rS[n]) #smoothing distance in z-direction
  thetaS.append(0.0*rand.randint(-45,45)) #add random orientation
  #thetaAirVec.append(F(center[n,0],center[n,1]))
  if(center[n,1]<0.0):
    thetaAirVec.append(Fth1(center[n,0]))
  else:
    thetaAirVec.append(Fth0(center[n,0]))

#### USE SPECIFIC VALUES OF ROUGHNESS HEIGHTS #####
print('--||ALYA ::  ROUGHNESS INFORMATION')
print('----||INFO ::  USING RANDOM ROUGHNESS HEIGHTS')
print('----||INFO :: USING N= ',Roughness,'  ROUGHNESS ELEMENTS')

#for n in range(0,Roughness):
#  print '----||R', n+1,'x0=',center[n,0],'y0=',center[n,1],'z0=',center[n,2],'h = ',rS[n],'th=',thetaAirVec[n]
#exit()

ii,jj = 0,0
done = False

xx,yy = [],[] #Full x-y plane visualization 
xx1,yy1 = [],[] #Arrow x-y plane visualization 
xx0,yy0 = [],[]
xx00,yy00 = [],[] #Arrow y-z plane visualization 
xx11,yy11 = [],[]
nMax = rS.index(max(rS))

print('--||ALYA :: RUNNING GRID MODIFICATION SCRIPT')

#Find a box to apply the grid modification
xBoxL = np.amin(center[:,0]) - 5.0*h
xBoxR = np.amax(center[:,0]) + 5.0*h
yBoxL = np.amin(center[:,1]) - 5.0*h
yBoxR = np.amax(center[:,1]) + 5.0*h

#Read each node, change coordinate if inside roughness
for line in f:
    if inCoords:
        ii = ii + 1
        #perDone = np.floor_divide(ii*100,nNodes)
        #if(perDone==jj):
        #    print '--||ALYA :: %.1f %% COMPLETE' % perDone
        #    jj = jj + 1
        if eString in line:
            print('----||DONE ::  TIME ',time.time()-startTime,'sec')
            inCoords = False
            lineB = True
            done = True
        else:
            line_elem = list(map(float, line.split()))
            x = line_elem[1]
            y = line_elem[2]
            z = line_elem[3]


            if  x >= xBoxL and x <= xBoxR and y>=yBoxL and y<=yBoxR:
                for n in range(0,Roughness):
                    aSphere = aS[n]
                    bSphere = bS[n]
                    cSphere = cS[n]
                    rSphere = rS[n]
                    tS      = m.pi*thetaS[n]/180.0
                    xc0 = center[n,0]
                    yc0 = center[n,1]
                    zc0 = center[n,2]
#                    #if zc0==0.0 or zc0==0.2:
#                    #  rSphere = h
                    thetaAirfoil = thetaAirVec[n]
#                    thetaAirfoil = thetaAirfoil + tS
                    s,t = unwrap_cyl(x, y, xc0, yc0)
#                
                    hx = s*m.cos(t-thetaAirfoil)
                    hy = s*m.sin(t-thetaAirfoil)
                    rad_loc = (x-xc0)**2.0/(aSphere)**2+(y-yc0)**2.0/(bSphere)**2+(z-zc0)**2.0/(cSphere)**2
                    if rad_loc <= 1.0:
                         # Write a plane for visualization
                        if z == 0.0:
                           xx0.append(x)
                           yy0.append(y)
                        hy_sphere = hy + rSphere*m.sqrt(abs(1.0-rad_loc))
                        
			# smoothing factors for x-y plane
                        inloc = (s - rSphere/4.0)
                        if inloc<=0.0:
                           fact = 0.0
                        else:
                           fact = 0.5/rSphere**2

			# smoothing factors for z direction
                        fact_z = (abs(z-zc0)/rSphere)
			#if n==25:
			#   db.write("[z-zc0] = %.6f \t [z-factor] = %.6f\n"%(z-zc0,fact_z))
			   #jj = jj + 1

                        hy_new = m.exp(-fact*inloc**2)*hy_sphere + (1.0-m.exp(-fact*inloc**2))*hy
                        #hy_new = m.exp(-fact_z)*hy_new + (1.0-m.exp(-fact_z))*hy

                        hy = ((H-hy)/H)*hy_new + (1.0 - (H-hy)/H)*hy

                        x = wrap_cyl(m.sqrt(hx**2+hy**2), m.atan2(hy,hx)+thetaAirfoil, xc0, yc0)[0]
                        y = wrap_cyl(m.sqrt(hx**2+hy**2), m.atan2(hy,hx)+thetaAirfoil, xc0, yc0)[1]
			#if(m.isnan(x) or m.isnan(y)):
			#   print('--|| ERROR :: NaN Values FOR R=',n)
			#   print('--|| INFO :: R=',n,xc0,yc0,zc0,thetaAirfoil)
			#   exit()

                        # Write a plane for visualization
                        if z == 0.0:
                           xx1.append(x)
                           yy1.append(y)

                g.write("%i %.17f %.17f %.17f\n"%(ii,x,y,z))
                #o.write("%i %.17f %.17f %.17f\n"%(ii,x,y,z))
                
#                # Write a plane for visualization
#                if z == 0.0:
#                    xx.append(x)
#                    yy.append(y)

            else:
                g.write("%i %.17f %.17f %.17f\n"%(ii,x,y,z))
                #o.write("%i %.17f %.17f %.17f\n"%(ii,line_elem[1],line_elem[2],z))


    if lineB:
        g.write(line)
        if bString in line:
            if fileType == 'msh':
               line = next(f) 
               g.write(line)
            if done == False:
                inCoords = True
                lineB = False
                print('--||ALYA :: START WRITING MODIFIED GEO FILE')
                startTime = time.time()

f.close()
g.close()
db.close()


###########################################################
###        PLOT FIGURES FOR ROUGHNESS VISUALIZATION  ######
#Plot a x-y plane

#plt.figure()
#plt.plot(xx[::1], yy[::1], 'b', marker=".", linewidth=0,  markersize=4)
#plt.plot(xx0[::1], yy0[::1], 'r', marker="o", linewidth=0,  markersize=4)
#plt.plot(center[:,0], center[:,1], 'k', marker="o", linewidth=0,  markersize=4)
#for i in range(0,len(xx0)):
#    plt.arrow(np.asarray(xx0[i]), np.asarray(yy0[i]), np.asarray(xx1[i])-np.asarray(xx0[i]), np.asarray(yy1[i])-np.asarray(yy0[i]), width=h/2000.0, head_width=h/500.0, head_length=h/500.0, fc='k', ec='k')
#plt.xlabel('$x$')
#plt.ylabel('$y$')
#plt.title('Un-wrapped cylinder view at $z=0.5$')
##plt.xlim(center[Roughness-7][0]-10.0*h,center[Roughness-7][0]+10.0*h )
##plt.ylim(center[Roughness-7][1]-10.0*h,center[Roughness-7][1]+10.0*h )
#plt.savefig('1-xy.eps', format='eps')
#plt.show()

#Plot a x-y plane

#plt.figure()
##plt.plot(xx[::1], yy[::1], 'b', marker=".", linewidth=0,  markersize=4)
#plt.plot(xx00[::1], yy00[::1], 'r', marker="o", linewidth=0,  markersize=4)
#for i in range(0,len(xx00)):
#    plt.arrow(np.asarray(xx00[i]), np.asarray(yy00[i]), np.asarray(xx11[i])-np.asarray(xx00[i]), np.asarray(yy11[i])-np.asarray(yy00[i]), width=h/2000.0, head_width=h/500.0, head_length=h/500.0, fc='k', ec='k')
#plt.xlabel('$y$')
#plt.ylabel('$z$')
#plt.title('Un-wrapped cylinder view at $x=$Max. Rough')
#plt.savefig('11-xy.eps', format='eps')
#plt.show()
#
##Plot a full Boxed x-y plane
#plt.figure()
#plt.plot(xx[::1], yy[::1], 'b', marker=".", linewidth=0,  markersize=4)
#plt.plot(center[:,0], center[:,1], 'k', marker="o", linewidth=0,  markersize=4)
#plt.xlabel('$x$')
#plt.ylabel('$y$')
#plt.title('Un-wrapped cylinder view at $z=0.5$')
#plt.savefig('2-xy.eps', format='eps')
#plt.show()
#
