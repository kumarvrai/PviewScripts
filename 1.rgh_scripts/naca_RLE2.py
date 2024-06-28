"""
Created on Sun Feb 17 17:53:23 2019

@author: dpastrana
"""
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

###########################################
fileType = sys.argv[1]
foilType = sys.argv[2]
RGH_SCRPTS='/home/kvishal/1.post_process/0.alya_pv_scripts/0.roughness_scripts/input_files/'
RGH_SCRPTS=os.path.expanduser('~')+'/1.post_process/0.alya_pv_scripts/0.roughness_scripts/input_files/'
if fileType == 'GMSH':
   geofile_in = 'naca.msh'
   geofile_out= 'naca_r.msh'
   bString = '$Nodes'
   eString = '$EndNodes'
elif fileType == 'ALYA':
   geofile_in = 'naca.geo.dat'
   geofile_out= 'naca_r.geo.dat'
   bString = 'COORDINATES'
   eString = 'END_COORDINATES'
if('4412' in foilType):    
  airfoilCoord = RGH_SCRPTS+'/naca4412.txt'
elif('0012' in foilType):  
  airfoilCoord = RGH_SCRPTS+'/naca0012.txt'
else:
  raise ValueError('--|| ALYA ERROR :: AIRFOIL TYPE NOT RECONIZED.')

coordAirfoil = np.loadtxt(airfoilCoord)

tanAirfoil = np.arctan2(np.diff(coordAirfoil[:,1]),np.diff(coordAirfoil[:,0]))
airLen = int(len(coordAirfoil)/2)
F0 = interp1d(coordAirfoil[0:airLen,0],coordAirfoil[0:airLen,1])
F1 = interp1d(coordAirfoil[airLen:2*airLen-1,0],coordAirfoil[airLen:2*airLen-1,1])
Fth0 = interp1d(coordAirfoil[0:airLen,0],tanAirfoil[0:airLen])
Fth1 = interp1d(coordAirfoil[airLen:2*airLen-1,0],tanAirfoil[airLen:2*airLen-1])
F = interp2d(coordAirfoil[0:2*airLen-1,0],coordAirfoil[0:2*airLen-1,1],tanAirfoil)
xAir,zAir = np.meshgrid(coordAirfoil[:,0],np.linspace(0,0.2,10))
yAir,zAir = np.meshgrid(coordAirfoil[:,1],np.linspace(0,0.2,10))
###########################################3
# Print information about the GEO file
nNodes = int(os.popen('sed -n ''/%s/,/%s/p'' %s | wc -l ' % (bString,eString,geofile_in)).read())
nNodes = nNodes-2
print('--||ALYA :: WORKING WITH %i NODES IN TOTAL' % nNodes)

#Load Files
f = open(geofile_in, 'r')
g = open(geofile_out, 'w')

lineB = True
inCoords = False

#PARAMETERS
iCenter = 21
xc0 = coordAirfoil[iCenter][0]
yc0 = 0.0
x_min = 0.0; x_max = 0.2;
z_min = 0.0; z_max = 0.2;
zdom = z_max-z_min

z_dom_fac = int((z_max-z_min)/0.2)


r0 = 0.5
D = 0.08
#h = 0.0049
#h = 0.013
#D = 0.1
h = 0.008
H = 1.0


fact = 0.0
b_radio = 10.0*h
center_y =  (h**2.0 + (0.5*D)**2.0)/(2.0*h) - h
rad = center_y+h

Roughness = -1
delta_z = 0.05 #distance in z between spheres
delta_x = 0.15  #distance in x between spheres
initial_pos = 17.45760312372209

###########################################3
# Choose roughness start and end locations (x,y)
xIniRough = 0.0
xFinRough = 0.01
nRough = 3
Xpos = np.linspace(xIniRough,xFinRough,nRough)
factRough = np.linspace(0.3,0.75,nRough)
xIniRough = 0.015
xFinRough = 0.06
nRough = 3 
Xpos = np.append(Xpos, np.linspace(xIniRough,xFinRough,nRough), axis=0)
factRough = np.append(factRough,np.linspace(0.3,1.0,nRough),axis=0)
xIniRough = 0.065
xFinRough = 0.1
nRough = 3 
Xpos = np.append(Xpos, np.linspace(xIniRough,xFinRough,nRough), axis=0)
factRough = np.append(factRough,np.linspace(1.0,0.75,nRough),axis=0)
xIniRough = 0.105
xFinRough = 0.18
nRough = 3 
Xpos = np.append(Xpos, np.linspace(xIniRough,xFinRough,nRough), axis=0)
factRough = np.append(factRough,np.linspace(0.75,0.3,nRough),axis=0)
#   Ypos.append(F1(Xpos[n]))
center = []
#center = np.meshgrid(Xpos,Ypos, Zpos)
#for z in Zpos:
#  for y in Ypos:
#    for x in Xpos:
#       center.append((x,y,z))
    
#p1 = skewnorm.pdf(Xpos, 100.0, loc=0.008, scale=1)
#p1 = (np.exp(-500.0*abs(Xpos-0.01)))
#p1 = 0.5*p1/np.amax(p1)
#p3 = (np.exp(-50.0*abs(Xpos-0.03)))
#p3 = p3/np.amax(p3)
#fact = np.power(p1 + p3,2)
#plt.plot(Xpos,p1,Xpos,p3,Xpos,np.power(p1 + p3,2))
#plt.show()
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#cset = ax.plot_surface(xAir,zAir,yAir,color='w')
#ax.clabel(cset, fontsize=9, inline=1)
#ax.plot3D(coordAirfoil[:,0],coordAirfoil[:,1],0.2,'k')
print('----||INFO. XLOC SS/PS SIDE =',Xpos)
nzR = np.floor_divide(zdom,h)
rS = []
# Roughness on Upper surface
for i in range(0,len(Xpos)):
  if(Xpos[i]<0.01 or Xpos[i]>=0.1):
     Zpos = np.linspace(z_min, z_max, z_dom_fac*rand.randint(4,8))
  else:
     Zpos = np.linspace(z_min, z_max, z_dom_fac*rand.randint(8,12))
  for k in range(0,len(Zpos)):
    if(Zpos[k]<=(np.amin(Zpos,axis=None)) or Zpos[k]>=(np.amax(Zpos,axis=None))):
       xLoc = Xpos[i]
    else:
       xLoc = Xpos[i]*(1.0+0.10*rand.randint(-1,1)/rand.randint(1,3)) 
    yLoc = F0(xLoc)
    center.append((xLoc,yLoc,Zpos[k]))
    rS.append(factRough[i] * h)
#    ax.scatter(xLoc,Zpos[k],yLoc,s=20, c='k')

# Roughness on Lower surface
for i in range(1,len(Xpos)):
  if(Xpos[i]<0.01 or Xpos[i]>=0.1):
     Zpos = np.linspace(z_min, z_max, z_dom_fac*rand.randint(4,8))
  else:
     Zpos = np.linspace(z_min, z_max, z_dom_fac*rand.randint(8,12))
  for k in range(0,len(Zpos)):
    if(Zpos[k]<=(np.amin(Zpos,axis=None)) or Zpos[k]>=(np.amax(Zpos,axis=None))):
       xLoc = Xpos[i]
    else:
       xLoc = Xpos[i]*(1.0+0.10*rand.randint(-1,1)/rand.randint(1,3)) 
    yLoc = F1(xLoc)
    center.append((xLoc,yLoc,Zpos[k]))
    rS.append(factRough[i] * h)

Roughness = len(center)
##############################################

aS = []
bS = []
cS = []
thetaS = []
thetaAirVec = []
for n in range(0,Roughness):
  aS.append(3.0*h)
  bS.append(3.0*h)
  cS.append(2.0*h)
  thetaS.append(0.0*rand.randint(-75,75))
  if(center[n][1]<0.0):
    thetaAirVec.append(Fth1(center[n][0]))
  else:
    thetaAirVec.append(Fth0(center[n][0]))

#### USE SPECIFIC VALUES OF ROUGHNESS HEIGHTS #####
print('--||ALYA ::  ROUGHNESS INFORMATION')
print('----||ROUGH ::  USING RANDOM ROUGHNESS HEIGHTS')
print('----||ROUGH :: USING N= ',Roughness,'  ROUGHNESS ELEMENTS')

rS = np.array(rS)
print('----||INFO :: ROUGHNESS PEAKS LOCATED AT ',len(np.unique(rS,axis=None)),'X LOCATIONS')

print('----||INFO :: MAX PEAK ROUGHNESS HEIGHT =',np.amax(rS,axis=None))
print('----||INFO :: AVG PEAK ROUGHNESS HEIGHT =',np.mean(np.unique(rS,axis=None),axis=None))
for n in range(0,Roughness):
  print('----|| R -',n+1,'[x0,y0,z0]=',np.array(center[n]),' h=',np.array(rS[n]), 'THETA=', thetaAirVec[n])
#exit()


ii = 0
done = False

xx = []
yy = []
xx1 = []
yy1 = []
xx0 = []
yy0 = []
tt = []
ss = []

centerArr = np.array(center,dtype=float);
#Find a box to apply the grid modification
xBoxL = np.amin(centerArr[:,0],axis=None) - 5.0*h
xBoxR = np.amax(centerArr[:,0],axis=None) + 5.0*h
yBoxL = np.amin(centerArr[:,1],axis=None) - 5.0*h
yBoxR = np.amax(centerArr[:,1],axis=None) + 5.0*h

#print('ALYA --|| SCRIPTING BOX X IN',xBoxL,xBoxR,'Y IN',yBoxL,yBoxR) 
#print('ALYA --|| INITIALIZE THE SCRIPT')

for line in f:
    if inCoords:
        ii = ii + 1
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
                    xc0 = center[n][0]
                    yc0 = center[n][1]
                    zc0 = center[n][2]
                    if zc0==np.amin(Zpos,axis=None) or zc0==np.amax(Zpos,axis=None):
                      rSphere = h
                    thetaAirfoil = thetaAirVec[n]
                    thetaAirfoil = thetaAirfoil + tS


                    s,t = unwrap_cyl(x, y, xc0, yc0)
                
                    hx = s*m.cos(t-thetaAirfoil)
                    hy = s*m.sin(t-thetaAirfoil)

                    rad_loc = (x-xc0)**2.0/(aSphere)**2 + (y-yc0)**2.0/(bSphere)**2 + (z-zc0)**2.0/(cSphere)**2
                    if rad_loc <= 1.0:
                         # Write a plane for visualization
                        if z == 0.0:
                           xx0.append(x)
                           yy0.append(y)
                        
                        x_old = x; y_old = y;
                        distPoint = m.sqrt((y-yc0)**2)
                        hy_sphere = hy + rSphere*m.sqrt(abs(1.0-rad_loc))
                        
                        inloc = (s - rSphere/4.0)
                        if inloc<=0.0:
                           fact = 0.0
                        else:
                           fact = 0.5/rSphere**2

			# smoothing factors for z direction
                        fact_z = (abs(z-zc0)/rSphere)

                        hy_new = m.exp(-fact*inloc**2)*hy_sphere + (1.0-m.exp(-fact*inloc**2))*hy
                        hy_new = m.exp(-fact_z)*hy_new + (1.0-m.exp(-fact_z))*hy

                        hy = ((H-hy)/H)*hy_new + (1.0 - (H-hy)/H)*hy

                        x = wrap_cyl(m.sqrt(hx**2+hy**2), m.atan2(hy,hx)+thetaAirfoil, xc0, yc0)[0]
                        y = wrap_cyl(m.sqrt(hx**2+hy**2), m.atan2(hy,hx)+thetaAirfoil, xc0, yc0)[1]

                        #if(m.sqrt((x-x_old)**2+(y-y_old)**2 > 1.1*h)):
			#   print('--|| ERROR :: HEIGHT GREATER THAN PRESCRIBED')
			#   exit()
                        #if(m.isnan(x) or m.isnan(y)):
                        #   print('--|| ERROR :: NaN Values FOR R=',n)
                        #   print('--|| INFO :: R=',n,xc0,yc0,zc0,thetaAirfoil)
                        #   exit()

                         # Write a plane for visualization
                        if z == 0.0:
                           xx1.append(x)
                           yy1.append(y)

                g.write("%i %.17f %.17f %.17f\n"%(ii,x,y,z))
                
                # Write a plane for visualization
                if z == 0.0:
                    xx.append(x)
                    yy.append(y)

            else:
                g.write("%i %.17f %.17f %.17f\n"%(ii,line_elem[1],line_elem[2],z))

    if lineB:
        g.write(line)
        if bString in line:
            if fileType == 'GMSH':
               line = next(f) 
               g.write(line)
            if done == False:
                inCoords = True
                lineB = False
                print('--||ALYA :: START WRITING MODIFIED GEO FILE')
                startTime = time.time()



f.close()
g.close()


############################################################
############################################################
####        PLOT FIGURES FOR ROUGHNESS VISUALIZATION  ######
#plt.figure()
#plt.plot(xx[::1], yy[::1], 'b', marker=".", linewidth=0,  markersize=4)
#plt.plot(xx0[::1], yy0[::1], 'r', marker="o", linewidth=0,  markersize=4)
##plt.plot(xx0[::1], yy0[::1], 'r', marker=".", linewidth=0,  markersize=1)
#for i in range(0,len(xx0)):
#    plt.arrow(np.asarray(xx0[i]), np.asarray(yy0[i]), np.asarray(xx1[i])-np.asarray(xx0[i]), np.asarray(yy1[i])-np.asarray(yy0[i]), width=h/2000.0, head_width=h/500.0, head_length=h/500.0, fc='k', ec='k')
##plt.quiver(np.asarray(xx0[::1]), np.asarray(yy0[::1]), np.asarray(xx1[::1])-np.asarray(xx0[::1]), np.asarray(yy1[::1])-np.asarray(yy0[::1]))
#plt.xlabel('$x$')
#plt.ylabel('$y$')
#plt.title('Un-wrapped cylinder view at $z=0.5$')
##plt.xlim(center[Roughness-7][0]-10.0*h,center[Roughness-7][0]+10.0*h )
##plt.ylim(center[Roughness-7][1]-10.0*h,center[Roughness-7][1]+10.0*h )
#plt.savefig('1.eps', format='eps')
##fig.savefig('1.svg', format='svg', dpi=1200)
#plt.show()
##
#plt.figure()
#plt.plot(xx[::1], yy[::1], 'b', marker=".", linewidth=0,  markersize=4)
#plt.xlabel('$x$')
#plt.ylabel('$y$')
#plt.title('Un-wrapped cylinder view at $z=0.5$')
##plt.savefig('2.eps', format='eps')
#fig.savefig('2.svg', format='svg', dpi=1200)
#plt.show()
##
