#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

avgStr = sys.argv[1]
Retau = float(sys.argv[2]) #Re_theta value
fileName = sys.argv[3]

strIdx = 100;

SOD_DIR='~/0.alya_pv_scripts/0.dns_data/blayer/'
decimals=3

plt.close("all");
plt.style.use('classic');
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 24};
plt.rc('font', **font);

mpl.rcParams['xtick.major.size'] = 5;
mpl.rcParams['xtick.major.width'] = 2;
mpl.rcParams['xtick.minor.size'] = 2.5;
mpl.rcParams['xtick.minor.width'] = 2;

mpl.rcParams['ytick.major.size'] = 5;
mpl.rcParams['ytick.major.width'] = 2;
mpl.rcParams['ytick.minor.size'] = 2.5;
mpl.rcParams['ytick.minor.width'] = 2;

#inputs
indices=['Points:0','AVVEL:0','AVPRE','AVVE2:0','AVVXY:0','AVRHO','AVMUE','AVVGR:0','AVVTW:0'];
index_var = np.zeros(len(indices),dtype=int);

fid = open('AvgData_3D.csv'); 
var = fid.readline();
var = var.replace('"','')
var = var.split(",")
for i in range(len(indices)):
  try:
    ind = var.index(indices[i])
  except:
    ind = -1
  index_var[i] = ind;

iX  = index_var[0];
iU  = index_var[1];
iP  = index_var[2];
iUU = index_var[3];
iUV = index_var[4];
iR  = index_var[5];
iMU = index_var[6];
idV = index_var[7];
idT = index_var[8];

#inputs
print("--||SOD: READING DNS DATA")
startTime = time.time()
#dns = np.loadtxt(SOD_DIR+'Re950_DNS.dat', delimiter=','); 
fname = SOD_DIR+'Re{}_DNS.dat'.format(int(Retau))
dns = np.loadtxt(fname, delimiter=','); 
print("--||SOD :: DONE. TIME =",time.time()-startTime,'sec')
print("--||SOD: READING CSV DATA")
startTime = time.time()
les = np.loadtxt('AvgData_3D.csv', delimiter=',', skiprows=1);
print("--||SOD :: DONE. TIME =",time.time()-startTime,'sec')

Uref=1.0;
delta =1.0;
rho = 1;
Re = np.exp((1.0/0.88)*np.log(Retau/0.09));
# mu = 5.923894374253487e-05
mu = (rho*2.0*delta*Uref)/Re; 
utau = (Retau*mu)/(delta*rho); 
print('--||SOD :: Input mu ',mu);
print('--||SOD :: Theory Reb ',Re);
print('--||SOD :: Theory utau ',utau);
print('--||SOD :: Theory PG  ',utau*utau*rho/delta);

# projection to cartesian mesh
x   = les[:,iX];
y   = les[:,iX+1];
z   = les[:,iX+2];
print('--||SOD :: Total Grids = {}'.format(len(x)))

x = np.around(x,decimals);
y = np.around(y,decimals);
z = np.around(z,decimals);

xyz = np.column_stack((x,y,z))
xyz,indx = np.unique(xyz,return_index=True, axis=0)
nxyz = len(indx)
print('--||SOD :: Unique Grids = {}'.format(nxyz))

les = les[indx,:];

x   = les[:,iX];
y   = les[:,iX+1];
z   = les[:,iX+2];
u   = les[:,iU];
v   = les[:,iU+1];
w   = les[:,iU+2];
r   = les[:,iR];
if(index_var[5]<0):
  r = 0.0*r+1.0;
p   = les[:,iP];
uu  = les[:,iUU];
vv  = les[:,iUU+1];
ww  = les[:,iUU+2];
uv  = les[:,iUV];
uw  = les[:,iUV+1];
vw  = les[:,iUV+2]
mue = les[:,iMU];
dudx = les[:,idV];
dudy = les[:,idV+1];
dudz = les[:,idV+2];
dvdx = les[:,idV+3];
dvdy = les[:,idV+4];
dvdz = les[:,idV+5];
dwdx = les[:,idV+6];
dwdy = les[:,idV+7];
dwdz = les[:,idV+8];
avtw = les[:,idT];
#avtw = np.sqrt(np.square(les[:,idT])+np.square(les[:,idT+1])+np.square(les[:,idT+2]));


xlin = np.unique(np.around(x,decimals));
ylin = np.unique(np.around(y,decimals));
zlin = np.unique(np.around(z,decimals));

indSort = np.lexsort((z,x,y))

print('--||SOD :: XDOMAIN = %.2f -- %.2f' % (np.amin(xlin),np.amax(xlin)))
print('--||SOD :: YDOMAIN = %.2f -- %.2f' % (np.amin(ylin),np.amax(ylin)))
print('--||SOD :: ZDOMAIN = %.2f -- %.2f' % (np.amin(zlin),np.amax(zlin)))

nx = len(xlin)
ny = len(ylin)
nz = len(zlin)

#nx = 127
#ny = 127
#nz = 127
print('--||SOD :: Mesh size : nx[{}] ny[{}] nz[{}]'.format(nx,ny,nz))

#x_3d,y_3d, z_3d = np.meshgrid(xlin,ylin,zlin, indexing='ij');
if(nx*ny*nz == nxyz):
  print('--||SOD :: Reshaping the input data')
  x_3d    = (   x[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  y_3d    = (   y[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  z_3d    = (   z[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  u_3d    = (   u[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  v_3d    = (   v[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  w_3d    = (   w[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  r_3d    = (   r[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  p_3d    = (   p[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  uu_3d   = (  uu[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  vv_3d   = (  vv[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  ww_3d   = (  ww[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  uv_3d   = (  uv[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  uw_3d   = (  uw[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  vw_3d   = (  vw[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  mue_3d  = ( mue[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dudx_3d = (dudx[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dudy_3d = (dudy[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dudz_3d = (dudz[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dvdx_3d = (dvdx[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dvdy_3d = (dvdy[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dvdz_3d = (dvdz[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dwdx_3d = (dwdx[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dwdy_3d = (dwdy[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  dwdz_3d = (dwdz[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  avtw_3d = (avtw[indSort].reshape((ny,nx,nz))).transpose(1, 0, 2)
  print('----||SOD :: Input data shape', np.shape(x_3d))
else:  
  x_3d,y_3d,z_3d = np.meshgrid(xlin,ylin,zlin, indexing='ij');
  print('--||SOD :: ijk mesh generated.')
  
  u_3d = griddata((x,y,z),u,(x_3d,y_3d,z_3d),method='nearest');
  v_3d = griddata((x,y,z),v,(x_3d,y_3d,z_3d),method='nearest');
  w_3d = griddata((x,y,z),w,(x_3d,y_3d,z_3d),method='nearest');
  r_3d = griddata((x,y,z),r,(x_3d,y_3d,z_3d),method='nearest');
  p_3d = griddata((x,y,z),p,(x_3d,y_3d,z_3d),method='nearest');
  
  uu_3d = griddata((x,y,z),uu,(x_3d,y_3d,z_3d),method='nearest');
  vv_3d = griddata((x,y,z),vv,(x_3d,y_3d,z_3d),method='nearest');
  ww_3d = griddata((x,y,z),ww,(x_3d,y_3d,z_3d),method='nearest');
  uv_3d = griddata((x,y,z),uv,(x_3d,y_3d,z_3d),method='nearest');
  uw_3d = griddata((x,y,z),uw,(x_3d,y_3d,z_3d),method='nearest');
  vw_3d = griddata((x,y,z),vw,(x_3d,y_3d,z_3d),method='nearest');
  mue_3d = griddata((x,y,z),mue,(x_3d,y_3d,z_3d),method='nearest');
  
  dudx_3d = griddata((x,y,z),dudx,(x_3d,y_3d,z_3d),method='nearest');
  dudy_3d = griddata((x,y,z),dudy,(x_3d,y_3d,z_3d),method='nearest');
  dudz_3d = griddata((x,y,z),dudz,(x_3d,y_3d,z_3d),method='nearest');
  dvdx_3d = griddata((x,y,z),dvdx,(x_3d,y_3d,z_3d),method='nearest');
  dvdy_3d = griddata((x,y,z),dvdy,(x_3d,y_3d,z_3d),method='nearest');
  dvdz_3d = griddata((x,y,z),dvdz,(x_3d,y_3d,z_3d),method='nearest');
  dwdx_3d = griddata((x,y,z),dwdx,(x_3d,y_3d,z_3d),method='nearest');
  dwdy_3d = griddata((x,y,z),dwdy,(x_3d,y_3d,z_3d),method='nearest');
  dwdz_3d = griddata((x,y,z),dwdz,(x_3d,y_3d,z_3d),method='nearest');
  avtw_3d = griddata((x,y,z),avtw,(x_3d,y_3d,z_3d),method='nearest');
  print('--||SOD :: Data projected to the ijk mesh.')

#z-average
u_3d    =    np.mean(u_3d, axis=2)
v_3d    =    np.mean(v_3d, axis=2)
w_3d    =    np.mean(w_3d, axis=2)
r_3d    =    np.mean(r_3d, axis=2)
p_3d    =    np.mean(p_3d, axis=2)
uu_3d   =   np.mean(uu_3d, axis=2)
vv_3d   =   np.mean(vv_3d, axis=2)
ww_3d   =   np.mean(ww_3d, axis=2)
uv_3d   =   np.mean(uv_3d, axis=2)
uw_3d   =   np.mean(uw_3d, axis=2)
vw_3d   =   np.mean(vw_3d, axis=2)
mue_3d  =  np.mean(mue_3d, axis=2)
dudx_3d = np.mean(dudx_3d, axis=2)  
dudy_3d = np.mean(dudy_3d, axis=2)  
dudz_3d = np.mean(dudz_3d, axis=2)  
dvdx_3d = np.mean(dvdx_3d, axis=2)  
dvdy_3d = np.mean(dvdy_3d, axis=2)  
dvdz_3d = np.mean(dvdz_3d, axis=2)  
dwdx_3d = np.mean(dwdx_3d, axis=2)  
dwdy_3d = np.mean(dwdy_3d, axis=2)  
dwdz_3d = np.mean(dwdz_3d, axis=2)  
avtw_3d = np.mean(avtw_3d, axis=2)
#x-average
if("X" in avgStr):
  bl_u    = np.mean(u_3d,    axis=0)
  bl_v    = np.mean(v_3d,    axis=0)
  bl_w    = np.mean(w_3d,    axis=0)
  bl_r    = np.mean(r_3d,    axis=0)
  bl_p    = np.mean(p_3d,    axis=0)
  bl_u2   = np.mean(uu_3d,   axis=0)
  bl_v2   = np.mean(vv_3d,   axis=0)
  bl_w2   = np.mean(ww_3d,   axis=0)
  bl_ux   = np.mean(uv_3d,   axis=0)
  bl_vx   = np.mean(uw_3d,   axis=0)
  bl_wx   = np.mean(vw_3d,   axis=0)
  bl_uu   = np.mean(uu_3d,   axis=0)-np.multiply(bl_u,bl_u)
  bl_vv   = np.mean(vv_3d,   axis=0)-np.multiply(bl_v,bl_v)
  bl_ww   = np.mean(ww_3d,   axis=0)-np.multiply(bl_w,bl_w)
  bl_uv   = np.mean(uv_3d,   axis=0)-np.multiply(bl_u,bl_v)
  bl_uw   = np.mean(uw_3d,   axis=0)-np.multiply(bl_u,bl_w)
  bl_vw   = np.mean(vw_3d,   axis=0)-np.multiply(bl_v,bl_w)
  bl_mue  = np.mean(mue_3d,  axis=0)
  bl_dudx = np.mean(dudx_3d, axis=0)  
  bl_dudy = np.mean(dudy_3d, axis=0)  
  bl_dudz = np.mean(dudz_3d, axis=0)  
  bl_dvdx = np.mean(dvdx_3d, axis=0)  
  bl_dvdy = np.mean(dvdy_3d, axis=0)  
  bl_dvdz = np.mean(dvdz_3d, axis=0)  
  bl_dwdx = np.mean(dwdx_3d, axis=0)  
  bl_dwdy = np.mean(dwdy_3d, axis=0)  
  bl_dwdz = np.mean(dwdz_3d, axis=0)  
  bl_avtw = np.mean(avtw_3d, axis=0)

  bl_u    = 1.0-bl_u;
  bl_u2   = bl_u2+2*bl_u-1;
  bl_uu   = bl_u2-np.multiply(bl_u,bl_u)

  #midline_size = int(ny);
  midline_size = np.argmin(abs(ylin-2.5));
  d99 = 1.0;
else:  
  #--blayer calculations
  nbl=0.95;
  N = np.shape(xlin)[0];
  print('--|| SOD :: Working with',N,'X points');
  iEdge = np.zeros(np.shape(xlin),dtype=int)
  uEdge = np.zeros(np.shape(xlin))
  pEdge = np.zeros(np.shape(xlin))
  delta99 = np.zeros(np.shape(xlin))
  delta_star = np.zeros(np.shape(xlin))
  delta_theta = np.zeros(np.shape(xlin))

  f = open('SodAvgData_2D.csv', 'w');
  f.write('Points:0,Points:1,Points:2,YPLUS,');
  f.write('AVVEL:0,AVVEL:1,AVVEL:2,AVPRE,AVRHO,AVMUE,');
  f.write('AVVE2:0,AVVE2:1,AVVE2:2,');
  f.write('AVVXY:0,AVVXY:1,AVVXY:2,');
  f.write('RS_II:0,RS_II:1,RS_II:2,');
  f.write('RS_IJ:0,RS_IJ:1,RS_IJ:2,');
  f.write('AVVGR:0,AVVGR:1,AVVGR:2,');
  f.write('AVVGR:3,AVVGR:4,AVVGR:5,');
  f.write('AVVGR:6,AVVGR:7,AVVGR:8,');
  f.write('Delta:0,Delta:1,Delta:2 \n');
  for i in range(N):
    eps   = ylin;
    uT    = u_3d[i,:];
    vT    = v_3d[i,:];
    wT    = w_3d[i,:];
    pT    = p_3d[i,:];
    rT    =    r_3d[i,:];
    u2T   =   uu_3d[i,:];
    v2T   =   vv_3d[i,:];
    w2T   =   ww_3d[i,:];
    uxT   =   uv_3d[i,:];
    vxT   =   uw_3d[i,:];
    wxT   =   vw_3d[i,:];
    uuT   =   uu_3d[i,:]-np.multiply(uT,uT);
    vvT   =   vv_3d[i,:]-np.multiply(vT,vT);
    wwT   =   ww_3d[i,:]-np.multiply(wT,wT);
    uvT   =   uv_3d[i,:]-np.multiply(uT,vT);
    uwT   =   uw_3d[i,:]-np.multiply(uT,wT);
    vwT   =   vw_3d[i,:]-np.multiply(vT,wT);
    mueT  =  mue_3d[i,:];
    dudxT = dudx_3d[i,:]; 
    dudyT = dudy_3d[i,:]; 
    dudzT = dudz_3d[i,:]; 
    dvdxT = dvdx_3d[i,:]; 
    dvdyT = dvdy_3d[i,:]; 
    dvdzT = dvdz_3d[i,:]; 
    dwdxT = dwdx_3d[i,:]; 
    dwdyT = dwdy_3d[i,:]; 
    dwdzT = dwdz_3d[i,:]; 
    avtwT = avtw_3d[i,:];
    #---------------------------------------------#
    # Griffin Method to find edge
    P0 = pT + 0.5*1.0*(np.square(uT) + np.square(vT));
    P0ref = np.amax(P0);
    uI = np.multiply(np.sign(uT),np.sqrt(abs(2*(P0ref-pT)-np.square(vT))));
    indEdge = np.argmin(abs(np.divide(uT,uI+1e-6)-nbl));
    iEdge[i] = indEdge;
    uEdge[i] = uT[indEdge];
    pEdge[i] = pT[indEdge]
    ue = uEdge[i];
    delta99[i] = eps[indEdge]
    delta_star[i] = np.trapz((1-uT[0:indEdge]/ue),eps[0:indEdge])
    delta_theta[i] = np.trapz(np.multiply(uT[0:indEdge]/ue,(1-uT[0:indEdge]/ue)),eps[0:indEdge])
    #print('idx ue d99 dstar dtht',indEdge,uEdge[i],delta99[i],delta_star[i],delta_theta[i])

    avgmu = mueT[0];
    tw = np.abs(avgmu*(uT[1])/ylin[1]);
    utau = np.sqrt(tw/rho);
    for j in range(len(ylin)):
        yplus = (np.absolute(ylin[j])*utau/avgmu);
        f.write('{},{},{},{},'.format(xlin[i],ylin[j],0.0,yplus));
        f.write('{},{},{},{},'.format(uT[j],vT[j],wT[j],pT[j]));
        f.write('{},{},'.format(rT[j],mueT[j]));
        f.write('{},{},{},'.format(u2T[j],v2T[j],w2T[j]));
        f.write('{},{},{},'.format(uxT[j],vxT[j],wxT[j]));
        f.write('{},{},{},'.format(uuT[j],vvT[j],wwT[j]));
        f.write('{},{},{},'.format(uvT[j],uwT[j],vwT[j]));
        f.write('{},{},{},'.format(dudxT[j],dudyT[j],dudzT[j]));
        f.write('{},{},{},'.format(dvdxT[j],dvdyT[j],dvdzT[j]));
        f.write('{},{},{},'.format(dwdxT[j],dwdyT[j],dwdzT[j]));
        f.write('{},{},{} \n'.format(delta99[i],delta_star[i],delta_theta[i]));
  f.close()

  Re99 = delta99/avgmu;
  ReStar = delta_star/avgmu;
  ReTheta = delta_theta/avgmu;

  #indX=int(len(xlin)/2);
  indX = np.argmin(abs(ReStar-Retau));
  print('--||SOD :: XLOC', xlin[indX]);
  print('--||SOD :: Re_d*', ReStar[indX]);
  print('--||SOD :: Re_th', ReTheta[indX]);
  d99 = delta99[indX];
  #midline_size = int(iEdge[indX]);
  midline_size = np.argmin(abs(ylin/d99-4.0));
  indY = range(midline_size);

  bl_u    = u_3d[indX,:];
  bl_v    = v_3d[indX,:];
  bl_w    = w_3d[indX,:];
  bl_r    = r_3d[indX,:];
  bl_p    = p_3d[indX,:];
  bl_u2   = uu_3d[indX,:]
  bl_v2   = vv_3d[indX,:]
  bl_w2   = ww_3d[indX,:]
  bl_ux   = uv_3d[indX,:]
  bl_vx   = uw_3d[indX,:]
  bl_wx   = vw_3d[indX,:]
  bl_uu   = uu_3d[indX,:]-np.multiply(bl_u,bl_u)
  bl_vv   = vv_3d[indX,:]-np.multiply(bl_v,bl_v)
  bl_ww   = ww_3d[indX,:]-np.multiply(bl_w,bl_w)
  bl_uv   = uv_3d[indX,:]-np.multiply(bl_u,bl_v)
  bl_uw   = uw_3d[indX,:]-np.multiply(bl_u,bl_w)
  bl_vw   = vw_3d[indX,:]-np.multiply(bl_v,bl_w)
  bl_mue  = mue_3d[indX,:];
  bl_dudx = dudx_3d[indX,:];  
  bl_dudy = dudy_3d[indX,:];  
  bl_dudz = dudz_3d[indX,:];  
  bl_dvdx = dvdx_3d[indX,:];  
  bl_dvdy = dvdy_3d[indX,:];  
  bl_dvdz = dvdz_3d[indX,:];  
  bl_dwdx = dwdx_3d[indX,:];  
  bl_dwdy = dwdy_3d[indX,:];  
  bl_dwdz = dwdz_3d[indX,:];  
  bl_avtw = avtw_3d[indX,:];

# get the averaged profiles
size_avg = nx*nz;

print('--||SOD :: Averages in x-z done.')

print('--||SOD :: mid_line',midline_size )

dx = xlin[1]-xlin[0];
dy = ylin[1]-ylin[0];
dz = zlin[1]-zlin[0];
#--if wall_model-------#
avgmu = mu;
tw = np.amax(np.abs(bl_avtw),axis=None)
#if(tw==0):
#  #--if no wall_model
#  print("--||SOD :: Calculating utau for WRLES")
#  #avgmu = 0.5*(bl_mue[0]+bl_mue[-1]);
#  avgmu = bl_mue[0];
#  print("----||SOD :: using FDiff")
#  #tw = np.abs(mu*(bl_u[1])/dy);
#  tw = np.abs(avgmu*(bl_u[1])/dy);
#else:
#  print("--||SOD :: Using utau for WMLES")
if(tw==0):
  #--if no wall_model
  print("--||SOD :: Calculating utau for WRLES")
  avgmu = bl_mue[0];
  file_2 = 'surf_code_1-'+fileName+'.dat';
  data_2 = np.loadtxt(file_2,delimiter=",",skiprows=1)
  area      = data_2[strIdx:,2]
  Fvx       = data_2[strIdx:,6];
  if(len(area)<100):
    print("----||SOD :: using FDiff")
    #tw = np.abs(mu*(bl_u[1])/dy);
    tw = np.abs(avgmu*(bl_u[1])/dy);
  else:  
    print("----||SOD :: using logfile")
    tw      = np.mean(abs(Fvx/area),axis=None)
  print('--||SOD :: MURAT WALL', avgmu/mu)
else:
  print("--||SOD :: Calculating utau for WMLES")

utauDNS  = utau;
utau = np.sqrt(tw/rho);
ubulk = np.trapz(bl_u,ylin)/(np.amax(ylin)-np.amin(ylin));

Re_sim = (2.0*d99*ubulk/avgmu);
Re_tau_sim = d99*utau/avgmu;
#Re_tau_sim = 0.09*np.exp(0.88*np.log(Re_sim));
#utauDNS = (Re_tau_sim*avgmu);
ubulkDNS = np.trapz(utauDNS*dns[:,2],dns[:,0])/(np.amax(dns[:,0])-np.amin(dns[:,0]));

err_utau = 100*(utau/ubulk-utauDNS/ubulkDNS)/(utauDNS/ubulkDNS)

print('--||SOD :: DNS Data Ub ',ubulkDNS)
print('--||SOD :: DNS Data utau ',utauDNS/ubulkDNS)
print('--||SOD :: Simulation mu ',avgmu)
print('--||SOD :: Simulation Ub ',ubulk)
print('--||SOD :: Simulation Reb ',Re_sim)
print('--||SOD :: Simulation utau ',utau/ubulk)
print('--||SOD :: Simulation Retau ',Re_tau_sim)
print('--||SOD :: Error tw %',err_utau)

bl_ustar   = np.zeros(midline_size);
bl_ystar   = np.zeros(midline_size);
bl_uustar  = np.zeros(midline_size);
bl_vvstar  = np.zeros(midline_size);
bl_wwstar  = np.zeros(midline_size);
bl_uvstar  = np.zeros(midline_size);
bl_uwstar  = np.zeros(midline_size);
bl_vwstar  = np.zeros(midline_size);

for j in range(midline_size):
    bl_ystar[j]  = ylin[j]*utau/avgmu;
    #bl_ystar[j]  = ylin[j]*utau/mu;
    bl_ustar[j]  = bl_u[j]/utau;
    bl_uustar[j] = np.sqrt(np.abs(bl_uu[j]))/(utau);
    bl_vvstar[j] = np.sqrt(np.abs(bl_vv[j]))/(utau);
    bl_wwstar[j] = np.sqrt(np.abs(bl_ww[j]))/(utau);
    bl_uvstar[j] = np.sqrt(np.abs(bl_uv[j]))/(utau);
    bl_uwstar[j] = np.sqrt(np.abs(bl_uw[j]))/(utau);
    bl_vwstar[j] = np.sqrt(np.abs(bl_vw[j]))/(utau);

print('--||SOD :: Bl averages done.')
print('--||SOD :: dx+ :', (rho*dx*utau/avgmu))
print('--||SOD :: dy+ :', np.amin(np.diff(bl_ystar)), np.amax(np.diff(bl_ystar)))
print('--||SOD :: dz+ :', (rho*dz*utau/avgmu))

# print the data

#f = open('SodAvgData_{}_{}_{}.csv'.format(len(xlin),len(ylin),len(zlin)), 'w');
f = open('SodAvgData_1D.csv', 'w');


f.write('Points:0,Points:1,Points:2,YPLUS,');
f.write('AVVEL:0,AVVEL:1,AVVEL:2,AVPRE,AVRHO,AVMUE,');
f.write('AVVE2:0,AVVE2:1,AVVE2:2,');
f.write('AVVXY:0,AVVXY:1,AVVXY:2,');
f.write('RS_II:0,RS_II:1,RS_II:2,');
f.write('RS_IJ:0,RS_IJ:1,RS_IJ:2,');
f.write('AVVGR:0,AVVGR:1,AVVGR:2,');
f.write('AVVGR:3,AVVGR:4,AVVGR:5,');
f.write('AVVGR:6,AVVGR:7,AVVGR:8 \n');
        
for j in range(len(ylin)):
    yplus = (np.absolute(1.0-ylin[j])*utau/avgmu);
    f.write('{},{},{},{},'.format(0.0,ylin[j],0.0,yplus));
    f.write('{},{},{},{},'.format(bl_u[j],bl_v[j],bl_w[j],bl_p[j]));
    f.write('{},{},'.format(bl_r[j],bl_mue[j]));
    f.write('{},{},{},'.format(bl_u2[j],bl_v2[j],bl_w2[j]));
    f.write('{},{},{},'.format(bl_ux[j],bl_vx[j],bl_wx[j]));
    f.write('{},{},{},'.format(bl_uu[j],bl_vv[j],bl_ww[j]));
    f.write('{},{},{},'.format(bl_uv[j],bl_uw[j],bl_vw[j]));
    f.write('{},{},{},'.format(bl_dudx[j],bl_dudy[j],bl_dudz[j]));
    f.write('{},{},{},'.format(bl_dvdx[j],bl_dvdy[j],bl_dvdz[j]));
    f.write('{},{},{} \n'.format(bl_dwdx[j],bl_dwdy[j],bl_dwdz[j]));
    #print('{},{},{}'.format(ylin[j],bl_ystar[j],bl_ustar[j]));
f.close()

# bl plots
lbl = r'$N_x={} \; N_y={} \;  N_z={}$'.format(len(xlin),len(ylin),len(zlin))

# U + 
fig=plt.figure(1, figsize=(8, 6), dpi=300)

plt.plot(bl_ystar[1:],bl_ustar[1:],'b',linewidth=3.0,label=lbl)
plt.plot(dns[:,1],dns[:,2],'k--',linewidth=3.0,label='DNS')
plt.axis([0.1, 2*Retau, 0, 1.2*np.amax(bl_ustar[1:],axis=None) ])
plt.xscale('log')
plt.ylabel(r'$U^+$')
plt.xlabel(r'$y^+$')
plt.tight_layout()
ax = fig.add_subplot(111)
legend = ax.legend(loc='upper left', fontsize=11)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
#plt.show()
plt.savefig('U+.png')


# Urms + 
fig=plt.figure(2, figsize=(8, 6), dpi=300)

plt.plot(bl_ystar[1:],bl_uustar[1:],'b',linewidth=3.0,label=r'$uu^+$')
plt.plot(bl_ystar[1:],bl_vvstar[1:],'r',linewidth=3.0,label=r'$vv^+$')
plt.plot(bl_ystar[1:],bl_wwstar[1:],'g',linewidth=3.0,label=r'$ww^+$')
plt.plot(dns[1:,1],dns[1:,3],'k--',linewidth=3.0,label='DNS')
plt.plot(dns[:,1],dns[:,4],'k--',linewidth=3.0)
plt.plot(dns[:,1],dns[:,5],'k--',linewidth=3.0)
plt.axis([0.1, 2*Retau, 0, 1.2*np.amax(bl_uustar[1:],axis=None)])
plt.xscale('log')
plt.ylabel(r'$U^{\prime +}$')
plt.xlabel(r'$y^+$')
plt.tight_layout()
ax = fig.add_subplot(111)
legend = ax.legend(loc='upper right', fontsize=11)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
#plt.show()
plt.savefig('Urms+.png')


# UVrms + 
fig=plt.figure(3, figsize=(8, 6), dpi=300)

plt.plot(bl_ystar[1:],bl_uvstar[1:],'b',linewidth=3.0,label=lbl)
plt.plot(dns[:,1],np.sqrt(abs(dns[:,6])),'k--',linewidth=3.0,label='DNS')
plt.axis([0.1, 2*Retau, 0, 1.2*np.amax(bl_uvstar[1:],axis=None)])
plt.xscale('log')
plt.ylabel(r'$U^{\prime}V^{\prime +}$')
plt.xlabel(r'$y^+$')
plt.tight_layout()
ax = fig.add_subplot(111)
legend = ax.legend(loc='upper left', fontsize=11)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
plt.show()
plt.savefig('UVrms+.png')

# U  
fig=plt.figure(4, figsize=(8, 6), dpi=300)

if('X' in avgStr):
  plt.plot(ylin,(bl_u)/ubulk,'b',linewidth=3.0,label=lbl)
else:
  plt.plot(ylin[indY]/delta99[indX],bl_u[indY]/ubulk,'b',linewidth=3.0,label=lbl)
plt.plot(dns[:,0],utauDNS*dns[:,2]/ubulkDNS,'k--',linewidth=3.0,label='DNS')
plt.ylabel(r'$U/U_b$')
plt.xlabel(r'$y/\delta_{99}$')
plt.tight_layout()
ax = fig.add_subplot(111)
legend = ax.legend(loc='upper left', fontsize=11)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
#plt.show()
plt.savefig('U.png')

# Blayer  
if('X' not in avgStr):
  fig=plt.figure(5, figsize=(8, 6), dpi=300)
  
  plt.plot(xlin,delta99,'o',linewidth=1.0,label=lbl)
  plt.axis([np.amin(xlin), np.amax(xlin), 0, 30.0])
  plt.ylabel(r'$\delta_{99}/\delta_0^*$')
  plt.xlabel(r'$x/\delta_0^*$')
  plt.tight_layout()
  ax = fig.add_subplot(111)
  #legend = ax.legend(loc='upper left', fontsize=11)
  for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
  #plt.show()
  plt.savefig('bLayer.png')
