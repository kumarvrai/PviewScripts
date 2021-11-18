

#### import modules
import numpy as np
from paraview.simple import *
import os
import sys 
from ToutanBatailleMaterialProp  import *
from importBoundaryData import *

#  pvpython ../POSTPROCESS/CalcLowMachAverages.py  "chan32T2" "LMchan003cold" 0 1 $nSamplePoints 298

casePath	=sys.argv[1]
caseName	=sys.argv[2]
outputName  =sys.argv[3]
yWall	   =float(sys.argv[4])
yMid		=float(sys.argv[5])
coordFil	=sys.argv[6]
TW		  =float(sys.argv[7])

if len(sys.argv) > 8:
	P		  =float(sys.argv[8])
else:
	P		  = 101325.0



yWallApprox=yWall+0.0001*(yMid-yWall)
#----------------------------------------------------------------# 

	
	
def getPtsData( _Obj, _T=-1.0, _Prop=None, _fname=None ): 
	from paraview.numpy_support import vtk_to_numpy

	if(_Prop==None):
	  print _Obj.GetPointDataInformation().keys()
	  sys.exit()

	_Obj.SMProxy.UpdatePipeline( _T )
	_Obj.UpdatePipelineInformation()
	B = servermanager.Fetch( _Obj )

	if( B.GetClassName()== 'vtkMultiBlockDataSet' ): B = extract_block(B)[0]

	Data = vtk_to_numpy( B.GetPointData().GetArray(_Prop) )

	if not _fname==None:
	  time = B.GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
	  print "  |_Time: %f" % ( time )
	  nout = "%s_%s" % (_fname, _Prop)
	  Fout = open(nout+".dat", "w")
	  Fout.close()
	  print "  |_'%s' " % ( Fout.name )

	return Data
	

####################################################################
####################################################################
####################################################################	


def getCellData( _Obj, _T=-1.0, _Prop=None, _fname=None ): 
	from paraview.numpy_support import vtk_to_numpy

	if(_Prop==None):
	  print _Obj.GetCellDataInformation().keys()
	  sys.exit()

	_Obj.SMProxy.UpdatePipeline( _T )
	_Obj.UpdatePipelineInformation()
	B = servermanager.Fetch( _Obj )

	if( B.GetClassName()== 'vtkMultiBlockDataSet' ): B = extract_block(B)[0]

	Data = vtk_to_numpy( B.GetCellData().GetArray(_Prop) )

	if not _fname==None:
	  time = B.GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
	  print "  |_Time: %f" % ( time )
	  nout = "%s_%s" % (_fname, _Prop)
	  Fout = open(nout+".dat", "w")
	  Fout.close()
	  print "  |_'%s' " % ( Fout.name )

	return Data
   

####################################################################
####################################################################
####################################################################   

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'EnSight Reader'
channEnsicase = EnSightReader(CaseFileName=os.path.join(casePath, caseName+'.ensi.case'))
channEnsicase.PointArrays = ['AVVEL', 'AVVE2','VELOC', 'CONDU', 'DENSI', 'ENTHA', 'SPECI', 'AVTEM', 'AVTE2', 'TEMPE', 'TURBU', 'VISCO']

# create a new 'Calculator'
calculator0 = Calculator(Input=channEnsicase)
calculator1 = Calculator(Input=calculator0)
calculator2 = Calculator(Input=calculator1)

calcCoord = ['X','Y','Z']
calcObj   = [calculator0,calculator1,calculator2]


for i in range(len(calcObj)):
	calc = calcObj[i]
	calc.ResultArrayName = "ReStress{}".format(i)
	calc.Function = "AVVE2_{0:}-AVVEL_{0:}*AVVEL_{0:}".format(calcCoord[i])
   

# temperature
calculator3 = Calculator(Input=calculator2)
calculator3.ResultArrayName = "Trms"
calculator3.Function = "sqrt(AVTE2-AVTEM*AVTEM)"

# mass flux
calculator4 = Calculator(Input=calculator3)
calculator4.ResultArrayName = "mFlux"
calculator4.Function = "DENSI*VELOC"

# kinematic viscosity
calculator5 = Calculator(Input=calculator4)
calculator5.ResultArrayName = "nu"
calculator5.Function = "VISCO/DENSI"

# temperature gradient
gradientOfUnstructuredDataSet1 = GradientOfUnstructuredDataSet(Input=calculator5)
gradientOfUnstructuredDataSet1.ScalarArray = ['POINTS', 'AVTEM']
gradientOfUnstructuredDataSet1.ResultArrayName = 'gradAVT'

# velocity gradient
gradientOfUnstructuredDataSet2 = GradientOfUnstructuredDataSet(Input=gradientOfUnstructuredDataSet1)
gradientOfUnstructuredDataSet2.ScalarArray = ['POINTS', 'AVVEL']
gradientOfUnstructuredDataSet2.ResultArrayName = 'gradAVV'


# create a new 'Temporal Statistics'
temporalStatistics1 = TemporalStatistics(Input=gradientOfUnstructuredDataSet2,ComputeMinimum = 0,ComputeMaximum = 0,ComputeStandardDeviation = 0)


#############################
# get node coordinates
#############################
coordData = np.loadtxt(os.path.join(casePath,coordFil))
yall= coordData[:,2]
yuni = sorted(np.unique(yall))

y=[]
yLoc=[]

ymin = min([yMid,yWall])
ymax = max([yMid,yWall])

tol=1e-7

for yi in yuni:
	if (yi>=ymin) and (yi<=ymax):
		if yi-ymin < tol:
			ytemp = yi + tol
		elif ymax-yi < tol:
			ytemp = yi - tol
		else:
			ytemp = yi

		y.append(ytemp)
		yLoc.append(abs(ytemp-yWall))

correctOrder=np.argsort(yLoc)
y 		= list(np.array(y)[correctOrder])
yLoc 	= list(np.array(yLoc)[correctOrder])



#############################
# sample data 
#############################
variableName=["AVVEL", "ReStress0", "ReStress1", "ReStress2", "CONDU", "DENSI", "ENTHA", "SPECI", "AVTEM", "TURBU", "VISCO", 'Trms', 'mFlux', 'nu', 'gradAVT','gradAVV']
variableName=[s+"_average" for s in variableName]

dataDict = {}

for var in variableName:
	dataDict[var]=[]


for j in range(len(y)):
	# create a new 'Slice'

	slice1 = Slice(Input=temporalStatistics1)
	slice1.SliceType = 'Plane'
	slice1.SliceOffsetValues = [0.0]
	slice1.SliceType.Normal = [0.0, 1.0, 0.0]
	slice1.SliceType.Origin = [3.14159265359, y[j] , 1.57079632679]

	#integrate over slice
	intVar  = IntegrateVariables(Input=slice1)
	
	area = getCellData( intVar , -1, "Area")
	
	prop = {}
	intProp={}
	for var in variableName:
		#print var
		intProp[var]	= getPtsData( intVar , -1, var)[0]
		prop[var]		= intProp[var]/area[0]
		dataDict[var].append(prop[var])


#################################################################
#################################################################
##########	Calculte scaling variables	#####################
#################################################################
f=open(os.path.join(casePath,caseName+'.ensi.case'),'r')

timeFound=False
tMin=np.inf
tMax=0

for line in f:
	#print line
	data = line.split()

	if timeFound:
		#print data[0]
		#print data[-1]
		tMin=min(tMin,float(data[0]))
		tMax=max(tMax,float(data[-1]))

	if len(data) >= 2 and data[0] == 'time' and data[1] == 'values:':
		#print data
		timeFound=True

f.close()

print 'tMin=', tMin
print 'tMax=', tMax



muW	 	= muTB(TW)
rhoW	= rhoTB(TW,P)
lambdaW = lambdaTB(TW)
CpW	 	= cpTB(TW)


print('{:<30}{:20.12f}'.format("muW= ", muW))
print('{:<30}{:20.12f}'.format("rhoW= ", rhoW))
print('{:<30}{:20.12f}'.format("lambdaW= ", lambdaW))
print('{:<30}{:20.12f}'.format("CpW= ", CpW))

halfHeight = abs(yMid-yWall)

tauW	= muW * abs(dataDict['gradAVV_average'][0][1])
uTau	= np.sqrt(tauW/rhoW)
Re		= (halfHeight*uTau*rhoW)/muW

QW		= -1.0*lambdaW*dataDict['gradAVT_average'][0][1]* (y[0]-yWall)/abs(y[0]-yWall)
TTau	= QW/(rhoW*CpW*uTau)


print('{:<30}{:20.12f}'.format("tauW= ",tauW))
print('{:<30}{:20.12f}'.format("utau= ",uTau))
print('{:<30}{:20.12f}'.format("ReTau= ",Re))
print('{:<30}{:20.12f}'.format("QW= ",QW))
print('{:<30}{:20.12f}'.format("TTau= ",TTau))
print('{:<30}{:20.12f}'.format("Bq= ",TTau/TW))


#################################################################
#################################################################
################	Calculte flux	############################
#################################################################


slice1 = Slice(Input=temporalStatistics1)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Normal = [1.0, 0.0, 0.0]
slice1.SliceType.Origin = [3.14159265359, 1.0 , 1.57079632679]

#integrate over slice
intVar  = IntegrateVariables(Input=slice1)

intProp=[]
for var in variableName:
	#print var
	intProp.append(getPtsData( intVar , -1, var)[0])

area = getCellData( intVar , -1, "Area")

volumeFlowRate = intProp[variableName.index('AVVEL_average')]
massFlowRate   = intProp[variableName.index('mFlux_average')]
uMean		   = intProp[variableName.index('AVVEL_average')][0]/area[0]
rhoMean	   =intProp[variableName.index('DENSI_average')]
muMean		=intProp[variableName.index('VISCO_average')]
ReBulk		=(halfHeight*uMean*rhoMean)/muMean
ReCenter	  =(halfHeight*dataDict['AVVEL_average'][-1][0]*dataDict['DENSI_average'][-1])/dataDict['VISCO_average'][-1]
ReCenter_nu   =(halfHeight*dataDict['AVVEL_average'][-1][0])/dataDict['nu_average'][-1]



print('{:<30}{:20.12f}'.format("uMean= ",uMean))
print('{:<30}{:20.12f}'.format("ReBulk= ",ReBulk))
print('{:<30}{:20.12f}'.format("ReCenter= ",ReCenter))
print('{:<30}{:20.12f}'.format("ReCenter_nu= ",ReCenter_nu))

#################################################################
#################################################################





f = open(outputName+'_1D.dat', 'w')
f.write("#  1: y\n")
f.write("#  2: y+\n")
f.write("#  3: Umean\n")
f.write("#  4: Urms\n")
f.write("#  5: Vrms\n")
f.write("#  6: Wrms\n")
f.write("#  7: Umean+\n")
f.write("#  8: Urms+\n")
f.write("#  9: Vrms+\n")
f.write("# 10: Wrms+\n")
f.write("# 11: Tmean\n")
f.write("# 12: Trms\n")
f.write("# 13: Tmean+\n")
f.write("# 14: Trms+\n")
f.write("# 15: mu_t\n")
f.write("# 16: rho\n")
f.write("# 17: lambda\n")
f.write("# 18: gradT\n")
f.write("# 19: Vmean\n")
f.write("# 20: Vmean+\n")
f.write("# y  y+  Umean  Urms  Vrms  Wrms  Umean+  Urms+  Vrms+  Wrms+  Tmean  Trms  Tmean+  Trms+ mu_t rho lambda gradT Vmean Vmean+\n")

for j in range(len(y)):
	#print y[j],ReStressXX[j]
	line=''
	line+= '{} {} '.format(yLoc[j], yLoc[j]*uTau*rhoW/muW)
	line+= '{} '.format(dataDict['AVVEL_average'][j][0]) 
	line+= '{} {} {} '.format(np.sqrt(dataDict['ReStress0_average'][j]), np.sqrt(dataDict['ReStress1_average'][j]), np.sqrt(dataDict['ReStress2_average'][j]))
	line+= '{} '.format(dataDict['AVVEL_average'][j][0]/uTau)
	line+= '{} {} {} '.format(np.sqrt(dataDict['ReStress0_average'][j])/uTau, np.sqrt(dataDict['ReStress1_average'][j])/uTau, np.sqrt(dataDict['ReStress2_average'][j])/uTau)
 
	line+= '{} '.format(dataDict['AVTEM_average'][j]) 
	line+= '{} '.format(dataDict['Trms_average'][j]) 
	line+= '{} '.format((TW-dataDict['AVTEM_average'][j])/TTau)  
	line+= '{} '.format(dataDict['Trms_average'][j]/np.fabs(TTau))  
	line+= '{} '.format(dataDict['TURBU_average'][j])  
	line+= '{} '.format(dataDict['DENSI_average'][j])  
	line+= '{} '.format(dataDict['CONDU_average'][j])  
	line+= '{} '.format(dataDict['gradAVT_average'][j][1])  
	line+= '{} '.format(dataDict['AVVEL_average'][j][1])  
	line+= '{} '.format(dataDict['AVVEL_average'][j][1]/uTau)  

	f.write(line+ "\n")

f.close()



f = open(outputName+'_0D.dat', 'w')
f.write("# mu_wall  rho_wall  y_wall  y_mid  halfHeight  tau_wall  u_tau  Re_tau  u_mean  Re_bulk t_min t_max\n")
line=''
line=line+ str(muW) + "  " + str(rhoW) + "  " + str(yWall)+ "  " + str(yMid) + "  " + str(halfHeight) + "  " + str(tauW) + "  " + str(uTau) + "  " + str(Re) + "  " + str(uMean) + "  " + str(ReBulk)  + "  " + str(tMin)  + "  " + str(tMax) 


f.write(line+ "\n")
f.close()



