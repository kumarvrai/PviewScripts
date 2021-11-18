import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

plt.close("all")

sphere = np.loadtxt('/Users/oriol/Dropbox/2019_lph_ijhff_old/data_sphere/symmetry_with_vorticity_plane00.csv', delimiter=',', skiprows=1)

x  = sphere[:,27];
y  = sphere[:,28];
wz = sphere[:,26];

npt  = 1001;
xlin = np.linspace(-1,3,npt);
ylin = np.linspace(-1.5,1.5,npt);

x_2d,y_2d = np.meshgrid(xlin,ylin);

wz_2d = griddata((x,y),wz,(x_2d,y_2d),method='linear');

for i in range(npt):
    for j in range(npt):
        r = np.sqrt(x_2d[i,j]**2 + y_2d[i,j]**2)
        if(r<= 0.5):
            wz_2d[i,j] = 0.0
        


levels = np.linspace(-5,5,256);

fig, ax1 = plt.subplots();
cf = ax1.contourf(x_2d, y_2d, wz_2d,levels, cmap='viridis',extend='both');

ax1.set_aspect('equal');
plt.tight_layout();
fig.colorbar(cf, ax=ax1,ticks=[-5, 0, 5],extendrect=True,extendfrac=0.0, fraction=0.0345, pad=0.05);

plt.show();

