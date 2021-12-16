import numpy as np
import os
import sys
from scipy.interpolate import griddata
#import matplotlib.pyplot as plt
import time

fileName        =       sys.argv[1]
airfoil         =	str(sys.argv[2])


geofile		= 	fileName+"_r.geo.dat"
coordFile	=	fileName+"Mesh.alya"

##########################################################################
def getCoordFunction(geofile,coordFile):
    f		=	open(geofile,'r')
    g		=	open(coordFile,'w')
    line	= 	True
    inCoord 	= 	False
    coordArray	=	np.empty([1,3])
    i		=	0 
    
    start_time = time.time()
    print("--||Alya: WRITING TARGET COORDINATES FROM GEO FILE")
    while line:
    	line = f.readline()
    	if inCoord:
    		if 'END' in line:
    			inCoord = False
    			line = False
    		else:
    			g.write(line)
    	if line:
    		if 'COORD' in line:
    			inCoord = True
    
    f.close()
    g.close()
    print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))

def getPerodicNodes(fileName,coordVec,perLabel,perPlane1,perPlane2):
  print("--||Alya: FINDING", perLabel, "PERIODIC NODES B/W", perPlane1,perPlane2)

  if perLabel=="z":
  	coordPlane1	=	coordVec[np.where(coordVec[:,3]==perPlane1)[0]]
  	coordPlane2	=	coordVec[np.where(coordVec[:,3]==perPlane2)[0]]
  elif perLabel=="x":
  	coordPlane1	=	coordVec[np.where(coordVec[:,1]==perPlane1)[0]]
  	coordPlane2	=	coordVec[np.where(coordVec[:,1]==perPlane2)[0]]
  
  shape1	= coordPlane1.shape
  shape2	= coordPlane2.shape
  
  if len(coordPlane1) == len(coordPlane2):
    print("--||Alya: START FINDING PERIODIC NODES")
    print("----||Alya: Shape-Plane1 = ", shape1)
    print("----||Alya: Shape-Plane2 = ", shape2)
  else:
    print("--||Alya: START FINDING PERIODIC NODES")
    print("----||Alya: Shape-Plane1 = ", shape1)
    print("----||Alya: Shape-Plane2 = ", shape2)
    print('--||Alya: Nodes-Plane1 =', (len(coordPlane1)))
    print('--||Alya: Nodes-Plane2 =', (len(coordPlane2)))
    sys.exit("--||Alya: ERROR: PLANES HAVE DIFFERENT NUMBER OF NODES")
  
  pairs	=	0
  
  if perLabel=="z":
    ind = np.lexsort((coordPlane1[:,2],coordPlane1[:,1]))
    coordPlane1 = coordPlane1[ind,:]
    ind = np.lexsort((coordPlane2[:,2],coordPlane2[:,1]))
    coordPlane2 = coordPlane2[ind,:]
    for i in range(0,len(coordPlane1)):
      if np.logical_and(coordPlane2[i,1] == coordPlane1[i,1],coordPlane2[i,2] == coordPlane1[i,2]):
        pairs = pairs + 1
      else:	
        sys.exit("--||Alya: ERROR: COULDN'T FIND PERIODIC NODE!")
    print('--||Alya:        END FINDING PERIODIC NODES')
  elif perLabel=="x":
          for i in range(0,len(coordPlane1)-1):
          	j=coordPlane2[np.where((coordPlane2[:,2] == coordPlane1[i,2]) & (coordPlane2[:,3] == coordPlane1[i,3])),0]
          	pairs = pairs +1
          print('--||Alya:        END FINDING PERIODIC NODES')
  print("--||Alya:        TOTAL", perLabel, "PERIODIC NODES = ", pairs)

################################################################################

# Read GEO file and Write point coordinates
#if(os.path.isfile(coordFile)):
#  print('--||Alya: COORD FILE ALREADY EXISTS')
#else:  
getCoordFunction(geofile,coordFile)

# Read point coordinates
start_time = time.time()
print("--||Alya: READING COORDINATES")
cf		=	open(coordFile,'r')
targetCoord = np.loadtxt(cf)
cf.close()
print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))

# Find periodic nodes
getPerodicNodes(fileName,np.around(targetCoord,decimals=6),'z', np.amin(targetCoord[:,3]), np.amax(targetCoord[:,3]))
