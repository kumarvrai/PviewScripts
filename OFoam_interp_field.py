import numpy as np
import os
import sys
from scipy.interpolate import griddata
#import matplotlib.pyplot as plt
import time

fileName        =     str(sys.argv[1])
airfoil         =     str(sys.argv[2])
varInterp       =     str(sys.argv[3])

tst = 0;
for air in ["0012", "4412", "2D", "DFUSER", "PIPE"]:
 if air in airfoil:
   tst = 1

if(tst==1):
  print ("--||Alya - WORKING WITH %s %s " % (fileName, airfoil))
else:
  sys.exit("--||Alya (ERROR): PROVIDED AIRFOIL IS NOT CODED!")

if varInterp in ["NONE", "VELOC", "TEMP", "BOTH"] :
   print ("--||Alya - INTERPOLATION FOR %s " % (varInterp))
else:
   sys.exit("--||Alya (ERROR): PROVIDED VARIABLE DOESN'T EXIST!")

coordFile     =     "cellCenters.txt"

def interpFunction(airfoil,targetCoord,var):
    if airfoil == "0012":
         caseName     =     "1.NACA0012/1.Re50000-AoA8/"
    elif airfoil == "0012pit":
         caseName     =     "1.NACA0012/2.Re50000-AoA0/"
    elif airfoil == "0012rod0":
         caseName     =     "1.NACA0012/3.Re50000-Rod/1.aoa_0/"
    elif airfoil == "0012rod8":
         caseName     =     "1.NACA0012/3.Re50000-Rod/1.aoa_0/"
    elif airfoil == "4412s1":
         caseName     =     "2.NACA4412/1.SMOOTH/1.Re-200K-AoA-5/"
    elif airfoil == "4412s2":
         caseName     =     "2.NACA4412/1.SMOOTH/2.Re-50K-AoA-8/"
    elif airfoil == "4412t":
         caseName     =     "2.NACA4412/2.TRIP/1.Re2e5/"
    elif airfoil == "4412twm":
         caseName     =     "2.NACA4412/2.TRIP/2.Re1e6/"
    elif airfoil == "4412r1":
         caseName     =     "2.NACA4412/3.ICE-ROUGH/1.R1/1.Re2e5/"
    elif airfoil == "4412r1wm":
         caseName     =     "2.NACA4412/3.ICE-ROUGH/1.R1/2.Re1e6/"
    elif airfoil == "4412r2":
         caseName     =     "2.NACA4412/3.ICE-ROUGH/1.R1/"
    elif airfoil == "4412r2wm":
         caseName     =     "2.NACA4412/3.ICE-ROUGH/2.R2/2.Re1e6/"
    elif airfoil == "2D":
         caseName     =     "2.BL_LAMINAR/"
    elif airfoil == "DFUSER":
         caseName     =     "6.Diffuser/1.Re19000/1.nek_data/"
    elif airfoil == "PIPE":
         caseName     =     "7.Pipe/"
    else:
         sys.exit("--||Alya (ERROR): PROVIDED AIRFOIL DOESN'T EXIST!")
         
    srcFolder     =     "/scratch/kvishal/1.Alya/0.PrePost/0.PreProc/0.InterpSource/"+caseName+"/"     
    srcCoordFile     =     srcFolder + "nacaRun00Mesh.dat"

    if var == "VELOC":
         srcVelFile     =     srcFolder + "nacaRun00.VELOC.dat"
         targetVelFile     =     "./VELOC.dat"
    elif var == "TEMP":
         srcVelFile     =     srcFolder + "nacaRun00.TEMP.dat"
         targetVelFile     =     "./TEMP.dat"
    else:
         sys.exit("--||Alya (ERROR): PROVIDED VARIABLE DOESN'T EXIST!")
    #     OPEN THE FILES FOR READING AND WRITING
    sf1               =       open(srcCoordFile,'r')
    sf2               =       open(srcVelFile,'r')
    tf2               =       open(targetVelFile,'w')
    
    #     READ SOURCE DATA AND TARGET COORDINATES
    start_time = time.time()
    print("--||Alya: READING SOURCE COORDINATES")
    srcCoord     =     np.loadtxt(sf1)
    sf1.close()
    print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))
    print("--||Alya: READING SOURCE VELOCITY/TEMPERATURE")
    start_time = time.time()
    srcVel          =     np.array(np.loadtxt(sf2))
    sf2.close()
    
    if airfoil=="2D":
       xSrc = np.array(srcCoord[:,1])
       ySrc = np.array(srcCoord[:,2])
       zSrc = np.array(0.0*srcCoord[:,1])
    else:
       xSrc = np.array(srcCoord[:,1])
       ySrc = np.array(srcCoord[:,2])
       zSrc = np.array(srcCoord[:,3])
    print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))
    
    
    #     PERFORM INTERPOLATION ON THE VELOCITIES
    if len(srcCoord) == len(srcVel) and var == "VELOC":
         if airfoil=="2D":
          xTrg = np.array(targetCoord[:,1])
          yTrg = np.array(targetCoord[:,2])
          zTrg = np.array(0.0*targetCoord[:,1])
         else:
          xTrg = np.array(targetCoord[:,0])
          yTrg = np.array(targetCoord[:,1])
          zTrg = np.array(targetCoord[:,2])
         print("--||Alya: NODES IN SOURCE FILES : (", (len(srcVel)), ")")
         print("--||Alya: NODES IN TARGET FILES : (", (len(targetCoord)), ")")
         fac_z = (np.amax(zSrc,axis=None)-np.amin(zSrc,axis=None))/ \
                (np.amax(zTrg,axis=None)-np.amin(zTrg,axis=None))
         print("--||Alya: STRETCHING THE SOURCE GRID",fac_z)

         start_time = time.time()
         print("--||Alya: INTERPOLATING VELOC-X" )
         tvx = griddata((xSrc,ySrc,zSrc/fac_z),srcVel[:,1], (xTrg,yTrg,zTrg), method='nearest')
         print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))
    
         start_time = time.time()
         print("--||Alya: INTERPOLATING VELOC-Y" )
         tvy = griddata((xSrc,ySrc,zSrc/fac_z),srcVel[:,2], (xTrg,yTrg,zTrg), method='nearest')
         print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))
    
         start_time = time.time()
         print("--||Alya: INTERPOLATING VELOC-Z" )
         tvz = griddata((xSrc,ySrc,zSrc/fac_z),srcVel[:,3], (xTrg,yTrg,zTrg), method='nearest')
         print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))
    
         start_time = time.time()

         print("--||Alya: WRITING VELOC.dat FILE" )
         tf2.write("velocityProfile\n")
         tf2.write("List<vector>\n")
         tf2.write("%i\n" % (len(targetCoord)))
         tf2.write("(\n")

         for i in range(0,len(targetCoord)):
              tf2.write("(%.8f\t%.8f\t%.8f)\n" % (tvx[i],tvy[i],tvz[i]))
         tf2.write(");\n")
         print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))

    elif len(srcCoord) == len(srcVel) and var == "TEMP":
         xTrg = np.array(targetCoord[:,1])
         yTrg = np.array(targetCoord[:,2])
         zTrg = np.array(targetCoord[:,3])
         print("--||Alya: NODES IN SOURCE FILES : (", (len(srcVel)), ")")
         print("--||Alya: NODES IN TARGET FILES : (", (len(targetCoord)), ")")
         start_time = time.time()
         print("--||Alya: INTERPOLATING TEMPERATURE" )

         tvx = griddata((xSrc,ySrc,zSrc),srcVel[:,1], (xTrg,yTrg,zTrg), method='nearest')
         print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))
    
         start_time = time.time()

         print("--||Alya: WRITING TEMP.dat FILE" )

         for i in range(0,len(targetCoord)):
              tf2.write("%i\t%.8f\n" % (i+1,tvx[i]))
         print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))
    else:
         sys.exit("--||Alya (ERROR): NODES IN COORD AND VEL FILES ARE THE SAME")
    tf2.close()
    print("")
    print("--||Alya: VELOC.dat/TEMP.dat file has been written!")


start_time = time.time()
print("--||Alya: READING TARGET COORDINATES")
cf          =     open(coordFile,'r')
targetCoord = np.loadtxt(cf)
cf.close()
print("----||DONE. Time Taken -- %s seconds" % (time.time() - start_time))

if varInterp == "VELOC":
   interpFunction(airfoil,targetCoord,"VELOC")
elif varInterp == "TEMP":
   interpFunction(airfoil,targetCoord,"TEMP")
elif varInterp == "BOTH":
   interpFunction(airfoil,targetCoord,"VELOC")
   interpFunction(airfoil,targetCoord,"TEMP")
elif varInterp == "NONE":
   print("--||Alya: SKIPING INTERPOLATION ROUTINE")
else:
   sys.exit("--||Alya (ERROR): PROVIDED VARIABLE DOESN'T EXIST!")

