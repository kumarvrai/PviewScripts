# Paraview Python script for integrating on a surface
# The script calculates average pressure and temperature on selected surface

# Start
print "Please wait. Integrating..."
srf=ExtractSurface() # Extract surface (if the selected object is not a surface)
nrm=GenerateSurfaceNormals(srf) # Generate surface normal vectors
# Calculate additional variables
# Mass flux
calcMassFlux=Calculator(nrm)
calcMassFlux.Function="AVVEL_average_avg.Normals"
calcMassFlux.ResultArrayName="MassFlux"
# Velocity magnitude
calcVelMag=Calculator(calcMassFlux) # Create calculator
calcVelMag.Function="((AVVEL_average_avg.iHat)^2+(AVVEL_average_avg.jHat)^2+(AVVEL_average_avg.kHat)^2)^0.5" # Function
calcVelMag.ResultArrayName="VelMag" # Result variable name
# Dynamic pressure
calcPrsDyn=Calculator(calcVelMag) # Create calculator
calcPrsDyn.Function="0.5*(VelMag^2)" # Function
calcPrsDyn.ResultArrayName="PrsDyn" # Result variable name
# Dynamic pressure multiplied by mass flux
calcPrsDynFlux=Calculator(calcPrsDyn) # Create calculator
calcPrsDynFlux.Function="PrsDyn*MassFlux" # Result variable name
calcPrsDynFlux.ResultArrayName="PrsDynFlux" # Function
## Temperature multiplied by mass flux
#calcTempFlux=Calculator(calcPrsDynFlux) # Create calculator
#calcTempFlux.Function="Temperature*MassFlux" # Function
#calcTempFlux.ResultArrayName="TempFlux" # Result variable name
# Integrate variables on selected surface with respect to area
integr=IntegrateVariables(calcPrsDynFlux) # Create integrator
vtkIntegr=servermanager.Fetch(integr) # Get native VTK integrator
# Get cell integrated data
integrRes=vtkIntegr.GetCellData() # Cell data array
integrArr=integrRes.GetArray("Area") # Area column
integrArea=integrArr.GetTuple(0)[0] # Area
# Get point integrated data
integrRes=vtkIntegr.GetPointData() # Point data array
integrArr=integrRes.GetArray("MassFlux") # Mass flux column
integrMassFlow=integrArr.GetTuple(0)[0] # Mass flow
integrArr=integrRes.GetArray("VelMag") # Velocity magnitude column
integrVelMag=integrArr.GetTuple(0)[0] # Velocity magnitude area integral
integrArr=integrRes.GetArray("AVPRE_average_avg") # Static pressure column
integrPrsStatArea=integrArr.GetTuple(0)[0] # Static pressure area integral
integrArr=integrRes.GetArray("PrsDynFlux") # Dynamic pressure flux column
integrPrsDynFlow=integrArr.GetTuple(0)[0] # Dynamic pressure flux area integral
##integrArr=integrRes.GetArray("TempFlux") # Temperature flux column
##integrTempFlow=integrArr.GetTuple(0)[0] # Temperature flux area integral
# Calculate average variable values
velMagAve=integrVelMag/integrArea # Area-averaged velocity magnitude [m/s]
prsStatAve=integrPrsStatArea/integrArea # Area-averaged static pressure [Pa]
prsDynAve=abs(integrPrsDynFlow/integrMassFlow) # Massflow-averaged dynamic pressure [Pa]
#tempAve=abs(integrTempFlow/integrMassFlow) # Massflow-averaged temparature [K]
#Find Bulk Flow rate
avgDatavtm = FindSource('AvgData.vtm')
plotOverLine1 = PlotOverLine(avgDatavtm)
plotOverLine1.Source='High Resolution Line Source'
plotOverLine1.Source.Point1 = [0.0, 1.0, 4.0]
plotOverLine1.Source.Point2 = [0.0, 4.0, 4.0]
integrBulkRate = IntegrateVariables(plotOverLine1)
vtkIntegr=servermanager.Fetch(integrBulkRate)
integrRes=vtkIntegr.GetPointData() # Point data array
integrArr=integrRes.GetArray("AVVEL_average_avg") # Mass flux column
integrBulkRate=integrArr.GetTuple(0)[0] # Mass flow
# Print results
print "Integrals on selected surface"
print "============================="
print "Bulk Velocity: ",integrBulkRate/3.0," [kg/s]"
print "Mass flow: ",integrMassFlow," [kg/s]"
print "Area-averaged velocity magnitude: ",velMagAve," [m/s]"
print "Area-averaged static pressure: ",prsStatAve," [Pa]"
print "Massflow-averaged dynamic pressure: ",prsDynAve," [Pa]"
print "Average total pressure: ",prsStatAve+prsDynAve," [Pa]"
#print "Massflow-averaged temperature: ",tempAve-273.15," [C]"
# Delete objects
Delete(integrBulkRate)
Delete(plotOverLine1)
Delete(integr) # Integrator
#Delete(calcTempFlux) # Calculator for temperature flux
Delete(calcPrsDynFlux) # Calculator for dynamic pressure flux
Delete(calcPrsDyn) # Calculator for dynamic pressure
Delete(calcVelMag) # Calculator for velocity magnitude
Delete(calcMassFlux) # Calculator for mass flux
Delete(nrm) # Surface normals
Delete(srf) # Extracted surface
