#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata


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
indices=['Points:0','AVVEL:0','AVPRE','AVVE2:0','AVVXY:0','AVVGR:0'];
index_var = np.zeros(len(indices),dtype=int);

SOD_DIR='/scratch/u/ugo/kvishal/research/0.Alya/4.CHAN/'
dns = np.loadtxt(SOD_DIR+'Re950_DNS.dat', delimiter=','); 

fid = open('AvgData_3D.csv'); 
var = fid.readline();
var = var.replace('"','')
var = var.split(",")
for i in range(len(indices)):
  ind = var.index(indices[i])
  index_var[i] = ind;

iX  = index_var[0];
iU  = index_var[1];
iP  = index_var[2];
iUU = index_var[3];
iUV = index_var[4];
#iR  = index_var[5];
#iMU = index_var[6];
idV = index_var[5];

#inputs
les = np.loadtxt('AvgData_3D.csv', delimiter=',', skiprows=1);
Uref=1.0;
delta =1.0;
rho = 1;
Retau=950.0;
Re = np.exp((1.0/0.88)*np.log(Retau/0.09));
# mu = 5.923894374253487e-05
mu = (rho*2.0*delta*Uref)/Re; 
utau = (Retau*mu)/(delta*rho); 
print(' Re ',Re);
print(' utau teor ',utau);
print(' mu ',mu);
print(' Gp  ',utau*utau*rho/delta);

# projection to cartesian mesh

x   = les[:,iX];
y   = les[:,iX+1];
z   = les[:,iX+2];
u   = les[:,iU];
v   = les[:,iU+1];
w   = les[:,iU+2];
#r   = les[:,iR];
r   = 0.0*y+1.0;
p   = les[:,iP];
uu  = les[:,iUU];
vv  = les[:,iUU+1];
ww  = les[:,iUU+2];
uv  = les[:,iUV];
uw  = les[:,iUV+1];
vw  = les[:,iUV+2]
#mue = les[:,iMU];
mue = 0.0*y+5.3566e-05;
dudx = les[:,idV];
dudy = les[:,idV+1];
dudz = les[:,idV+2];
dvdx = les[:,idV+3];
dvdy = les[:,idV+4];
dvdz = les[:,idV+5];
dwdx = les[:,idV+6];
dwdy = les[:,idV+7];
dwdz = les[:,idV+8];


xlin = np.unique(x);
ylin = np.unique(y);
zlin = np.unique(z);
print('XDOMAIN = %.2f -- %.2f' % (np.amin(xlin),np.amax(xlin)))
print('YDOMAIN = %.2f -- %.2f' % (np.amin(ylin),np.amax(ylin)))
print('ZDOMAIN = %.2f -- %.2f' % (np.amin(zlin),np.amax(zlin)))

print('Mesh size : nx[{}] ny[{}] nz[{}]'.format(len(xlin),len(ylin),len(zlin)))

x_3d,y_3d, z_3d = np.meshgrid(xlin,ylin,zlin, indexing='ij');


print('ijk mesh generated.')

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
print('Data projected to the ijk mesh.')

# get the averaged profiles

bl_u   = np.zeros(len(ylin));
size_avg = len(xlin)*len(zlin);

bl_u   = np.zeros(len(ylin));  #U
bl_v   = np.zeros(len(ylin));
bl_w   = np.zeros(len(ylin));
bl_p   = np.zeros(len(ylin));
bl_r   = np.zeros(len(ylin));
bl_uu  = np.zeros(len(ylin));  #RS_II
bl_vv  = np.zeros(len(ylin));
bl_ww  = np.zeros(len(ylin));
bl_uv  = np.zeros(len(ylin));  #RS_IJ
bl_uw  = np.zeros(len(ylin));
bl_vw  = np.zeros(len(ylin));
bl_u2  = np.zeros(len(ylin));  #AVVE2
bl_v2  = np.zeros(len(ylin));
bl_w2  = np.zeros(len(ylin));
bl_ux  = np.zeros(len(ylin));  #AVVXY
bl_vx  = np.zeros(len(ylin));
bl_wx  = np.zeros(len(ylin));
bl_dudx  = np.zeros(len(ylin));  #AVVGR
bl_dudy  = np.zeros(len(ylin));
bl_dudz  = np.zeros(len(ylin));
bl_dvdx  = np.zeros(len(ylin));  #AVVGR
bl_dvdy  = np.zeros(len(ylin));
bl_dvdz  = np.zeros(len(ylin));
bl_dwdx  = np.zeros(len(ylin));  #AVVGR
bl_dwdy  = np.zeros(len(ylin));
bl_dwdz  = np.zeros(len(ylin));
bl_mue = np.zeros(len(ylin));

#u_3d = u_3d/r_3d;
#v_3d = v_3d/r_3d;
#w_3d = w_3d/r_3d;

#uu_3d = uu_3d/r_3d;
#vv_3d = vv_3d/r_3d;
#ww_3d = ww_3d/r_3d;
#uv_3d = uv_3d/r_3d;
#uw_3d = uw_3d/r_3d;
#vw_3d = vw_3d/r_3d;

for j in range(len(ylin)):
    for i in range(len(xlin)):
        for k in range(len(zlin)):
            bl_u[j]      = bl_u[j]  + u_3d[i,j,k];
            bl_v[j]      = bl_v[j]  + v_3d[i,j,k];
            bl_w[j]      = bl_w[j]  + w_3d[i,j,k];
            bl_p[j]      = bl_p[j]  + p_3d[i,j,k];
            bl_r[j]      = bl_r[j]  + r_3d[i,j,k];
            bl_mue[j]    = bl_mue[j]  + mue_3d[i,j,k];
            bl_uu[j]     = bl_uu[j] + (uu_3d[i,j,k]-u_3d[i,j,k]*u_3d[i,j,k]);
            bl_vv[j]     = bl_vv[j] + (vv_3d[i,j,k]-v_3d[i,j,k]*v_3d[i,j,k]);
            bl_ww[j]     = bl_ww[j] + (ww_3d[i,j,k]-w_3d[i,j,k]*w_3d[i,j,k]);
            bl_uv[j]     = bl_uv[j] + (uv_3d[i,j,k]-u_3d[i,j,k]*v_3d[i,j,k]);
            bl_uw[j]     = bl_uw[j] + (uw_3d[i,j,k]-u_3d[i,j,k]*w_3d[i,j,k]);
            bl_vw[j]     = bl_vw[j] + (vw_3d[i,j,k]-v_3d[i,j,k]*w_3d[i,j,k]);
            bl_u2[j]     = bl_u2[j] + uu_3d[i,j,k];
            bl_v2[j]     = bl_v2[j] + vv_3d[i,j,k];
            bl_w2[j]     = bl_w2[j] + ww_3d[i,j,k];
            bl_ux[j]     = bl_ux[j] + uv_3d[i,j,k];
            bl_vx[j]     = bl_vx[j] + uw_3d[i,j,k];
            bl_wx[j]     = bl_wx[j] + vw_3d[i,j,k];
            bl_dudx[j]   = bl_dudx[j]  + dudx_3d[i,j,k];
            bl_dudy[j]   = bl_dudy[j]  + dudy_3d[i,j,k];
            bl_dudz[j]   = bl_dudz[j]  + dudz_3d[i,j,k];
            bl_dvdx[j]   = bl_dvdx[j]  + dvdx_3d[i,j,k];
            bl_dvdy[j]   = bl_dvdy[j]  + dvdy_3d[i,j,k];
            bl_dvdz[j]   = bl_dvdz[j]  + dvdz_3d[i,j,k];
            bl_dwdx[j]   = bl_dwdx[j]  + dwdx_3d[i,j,k];
            bl_dwdy[j]   = bl_dwdy[j]  + dwdy_3d[i,j,k];
            bl_dwdz[j]   = bl_dwdz[j]  + dwdz_3d[i,j,k];
    bl_u[j]  = bl_u[j]/size_avg;
    bl_v[j]  = bl_v[j]/size_avg;
    bl_w[j]  = bl_w[j]/size_avg;
    bl_p[j]  = bl_p[j]/size_avg;
    bl_r[j]  = bl_r[j]/size_avg;
    bl_mue[j]  = bl_mue[j]/size_avg;
    bl_uu[j] = bl_uu[j]/size_avg;
    bl_vv[j] = bl_vv[j]/size_avg;
    bl_ww[j] = bl_ww[j]/size_avg;
    bl_uv[j] = bl_uv[j]/size_avg;
    bl_uw[j] = bl_uw[j]/size_avg;
    bl_vw[j] = bl_vw[j]/size_avg;
    bl_u2[j] = bl_u2[j]/size_avg;
    bl_v2[j] = bl_v2[j]/size_avg;
    bl_w2[j] = bl_w2[j]/size_avg;
    bl_ux[j] = bl_ux[j]/size_avg;
    bl_wx[j] = bl_vx[j]/size_avg;
    bl_wx[j] = bl_wx[j]/size_avg;
    bl_dudx[j] = bl_dudx[j]/size_avg;
    bl_dudy[j] = bl_dudy[j]/size_avg;
    bl_dudz[j] = bl_dudz[j]/size_avg;
    bl_dvdx[j] = bl_dvdx[j]/size_avg;
    bl_dvdy[j] = bl_dvdy[j]/size_avg;
    bl_dvdz[j] = bl_dvdz[j]/size_avg;
    bl_dwdx[j] = bl_dwdx[j]/size_avg;
    bl_dwdy[j] = bl_dwdy[j]/size_avg;
    bl_dwdz[j] = bl_dwdz[j]/size_avg;

print('Averages in x-z done.')

midline_size = int(0.5*len(ylin));
line_size    = len(ylin)-1;

print('mid_line',midline_size )

#for j in range(len(ylin)):
#    print('{},{},{}'.format(j,ylin[j],2-ylin[line_size-j]));

# for j in range(midline_size):
#     bl_u[j]  = 0.5*(bl_u[j] + bl_u[line_size-j]);
#     bl_v[j]  = 0.5*(-bl_v[j] + bl_v[line_size-j]);
#     bl_w[j]  = 0.5*(bl_w[j] + bl_w[line_size-j]);
#     bl_uu[j] = 0.5*(bl_uu[j]+ bl_uu[line_size-j]);
    # bl_vv[j] = 0.5*(bl_vv[j]+ bl_vv[line_size-j]);
    # bl_ww[j] = 0.5*(bl_ww[j]+ bl_ww[line_size-j]);

dx = xlin[1]-xlin[0];
dy = ylin[1]-ylin[0];
dz = zlin[1]-zlin[0];
avgmu = 0.5*(bl_mue[0]+bl_mue[-1]);
tw = np.abs(avgmu*(bl_u[1])/dy);
#tw = np.abs(mu*(bl_u[1])/dy);

print('Delta y:',dy)
print('tw ',tw)
print('theory tw ',utau**2*rho)
print('err tw %',100*np.abs(tw-utau**2*rho)/(utau**2*rho))
utauDNS=utau
utau = np.sqrt(tw/rho);
print('utau ',utau)
print('Re_tau ',rho*utau*delta/mu)

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
    bl_uustar[j] = np.sqrt(np.abs(bl_uu[j]))/utau;
    bl_vvstar[j] = np.sqrt(np.abs(bl_vv[j]))/utau;
    bl_wwstar[j] = np.sqrt(np.abs(bl_ww[j]))/utau;
    bl_uvstar[j] = np.sqrt(np.abs(bl_uv[j]))/utau;
    bl_uwstar[j] = np.sqrt(np.abs(bl_uw[j]))/utau;
    bl_vwstar[j] = np.sqrt(np.abs(bl_vw[j]))/utau;

print('Bl averages done.')
print('dx+ :', (rho*dx*utau/avgmu))
print('dy+ :', np.amin(np.diff(bl_ystar)), np.amax(np.diff(bl_ystar)))
print('dz+ :', (rho*dz*utau/avgmu))

# print the data

#f = open('SodAvgData_{}_{}_{}.csv'.format(len(xlin),len(ylin),len(zlin)), 'w');
f = open('AlyaAvgData_1D.csv', 'w');


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

#plots
lbl = r'$N_x={} \; N_y={} \;  N_z={}$'.format(len(xlin),len(ylin),len(zlin))

# U + 
fig=plt.figure(1, figsize=(8, 6), dpi=300)

plt.plot(bl_ystar[1:],bl_ustar[1:],'b',linewidth=3.0,label=lbl)
plt.plot(dns[:,1],dns[:,2],'k--',linewidth=3.0,label='DNS')
plt.axis([0.1, 2000, 0, 30 ])
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

plt.plot(bl_ystar[1:],bl_uustar[1:],'b',linewidth=3.0,label=lbl)
plt.plot(dns[1:,1],dns[1:,3],'k--',linewidth=3.0,label='DNS')
plt.plot(bl_ystar[1:],bl_vvstar[1:],'r',linewidth=3.0)
plt.plot(dns[:,1],dns[:,4],'k--',linewidth=3.0)
plt.plot(bl_ystar[1:],bl_wwstar[1:],'g',linewidth=3.0)
plt.plot(dns[:,1],dns[:,5],'k--',linewidth=3.0)
plt.axis([0.1, 2000, 0, 3.5])
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
fig=plt.figure(5, figsize=(8, 6), dpi=300)

plt.plot(bl_ystar[1:],-bl_uvstar[1:],'b',linewidth=3.0,label=lbl)
plt.plot(dns[:,1],dns[:,5],'k--',linewidth=3.0,label='DNS')
plt.axis([0.1, 2000, 0, 2])
plt.xscale('log')
plt.ylabel(r'$U^{\prime}V^{\prime +}$')
plt.xlabel(r'$y^+$')
plt.tight_layout()
ax = fig.add_subplot(111)
legend = ax.legend(loc='upper right', fontsize=11)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
plt.show()
plt.savefig('UVrms+.png')

# U  
fig=plt.figure(3, figsize=(8, 6), dpi=300)

plt.plot(ylin,bl_u,'b',linewidth=3.0,label=lbl)
plt.plot(dns[:,0],utauDNS*dns[:,2],'k--',linewidth=3.0,label='DNS')
plt.ylabel(r'$U$')
plt.xlabel(r'$y$')
plt.tight_layout()
ax = fig.add_subplot(111)
legend = ax.legend(loc='upper left', fontsize=11)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
#plt.show()
plt.savefig('U.png')
