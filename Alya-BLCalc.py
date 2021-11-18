import numpy as np
from scipy.interpolate import griddata
from paraview.vtk.numpy_interface import dataset_adapter as dsa
d = dsa.WrapDataObject(inputs[0].GetBlock(0))
x = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
y = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
z = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=6)
nVec = -np.asarray(d.PointData["Normals"],dtype=np.double)
tVec = np.zeros(np.shape(nVec))
tVec[:,0] = np.multiply(nVec[:,1],np.sign(y))
tVec[:,1] = np.multiply(-nVec[:,0],np.sign(y))
#---------------------------------------------#
d = dsa.WrapDataObject(inputs[1].GetBlock(0))
xS = np.around(np.asarray(d.Points[:,0],dtype=np.double),decimals=6)
yS = np.around(np.asarray(d.Points[:,1],dtype=np.double),decimals=6)
zS = np.around(np.asarray(d.Points[:,2],dtype=np.double),decimals=6)
uS = np.asarray(d.PointData["VELOC"],dtype=np.double)[:,0]
vS = np.asarray(d.PointData["VELOC"],dtype=np.double)[:,1]
pS = np.asarray(d.PointData["PRESS"],dtype=np.double)

#---------------------------------------------#
nbl=0.9;
nu = 5e-6;
N = np.shape(x)[0]; NInt = 100;
print('--|| WORKING ON',N,'POINTS');
xNew = np.zeros(np.shape(x))
yNew = np.zeros(np.shape(x))
zNew = np.zeros(np.shape(x))
uEdge = np.zeros(np.shape(x))
delta99 = np.zeros(np.shape(x))
delta_star = np.zeros(np.shape(x))
delta_theta = np.zeros(np.shape(x))

for i in range(N):
  if(np.logical_and(x[i]>0.02,x[i]<0.98)):
    eps = np.linspace(0,5.0*(10.0*np.sqrt(nu*x[i])+0.01),NInt);
    xNew = x[i]+eps*nVec[i,0]
    yNew = y[i]+eps*nVec[i,1]
    uInterp = griddata((xS,yS), uS, (xNew,yNew),method='nearest')
    vInterp = griddata((xS,yS), vS, (xNew,yNew),method='nearest')
    pInterp = griddata((xS,yS), pS, (xNew,yNew),method='nearest')
    uT = uInterp*tVec[i,0] + vInterp*tVec[i,1]
    vT = uInterp*nVec[i,0] + vInterp*nVec[i,1]
    # Griffin Method to find edge
    P0 = pInterp + 0.5*1.0*(np.square(uT) + np.square(vT));
    P0ref = np.amax(P0);
    uI = np.multiply(np.sign(uT),np.sqrt(abs(2*(P0ref-pInterp)-np.square(vT))));
    indEdge = np.argmin(abs(np.divide(uT,uI+1e-6)-nbl));
    uEdge[i] = uT[indEdge]; ue = uEdge[i];
    delta99[i] = eps[indEdge]
    delta_star[i] = np.trapz((1-uT[0:indEdge]/ue),eps[0:indEdge])
    delta_theta[i] = np.trapz(np.multiply(uT[0:indEdge]/ue,(1-uT[0:indEdge]/ue)),eps[0:indEdge])

#-----------------------------------#
uEdge = np.asarray(uEdge, dtype=np.float64)
#output.ShallowCopy(inputs[0].VTKObject)
output.PointData.append(nVec, "NVEC")
output.PointData.append(tVec, "TVEC")
output.PointData.append(uEdge, "UEDGE")
output.PointData.append(delta99, "DEL99")
output.PointData.append(delta_star, "DSTAR")
output.PointData.append(delta_theta, "DTHTA")
