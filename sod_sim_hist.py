import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys

fileName = sys.argv[1]
data = np.loadtxt('analysis_'+fileName+'.dat',skiprows=0)
data_2 = np.loadtxt('surf_code_1-'+fileName+'.dat',delimiter=",",skiprows=1)

time    = data[:,0];
Ek      = data[:,1];
eps_S   = data[:,2];
eps_D   = data[:,3];
eps_T   = data[:,4];
Mmax    = data[:,5];

time_2    = data_2[:,1];
area      = data_2[:,2]
Fpx       = data_2[:,3];
Fpy       = data_2[:,4];
Fvx       = data_2[:,6];
Fvy       = data_2[:,7];

utau      = np.sqrt(abs(Fvx/area))

ylab = [r'$E_k$',r'$eps_{S}$',r'$eps_{D}$',r'$eps_{T}$',r'$M_{max}$','']

#fig, axs = plt.subplots(3, 2, sharex=True, dpi=300)
fig = plt.figure(dpi=300)
grid = GridSpec(3, 2, figure=fig)

axs = fig.add_subplot(grid[0,0])
axs.plot(time,Ek,linewidth=1.5,label='Ek')
plt.xlabel(r'$t$')
plt.ylabel(ylab[0])
axs = fig.add_subplot(grid[0,1])
axs.plot(time,eps_S,linewidth=1.5,label=r'$eps_S$')
plt.xlabel(r'$t$')
plt.ylabel(ylab[1])
axs = fig.add_subplot(grid[1,0])
axs.plot(time,eps_D,linewidth=1.5,label=r'$eps_D$')
plt.xlabel(r'$t$')
plt.ylabel(ylab[2])
axs = fig.add_subplot(grid[1,1])
#axs[1,1].plot(time,eps_T,linewidth=1.5,label=r'$eps_T$')
axs.plot(time,Mmax,linewidth=1.5,label=r'$M_{max}$')
plt.xlabel(r'$t$')
plt.ylabel(ylab[4])
axs = fig.add_subplot(grid[2,:])
axs.plot(time,utau,linewidth=1.5,label=r'$u_{\tau}$')
#axs[2,1].plot(time,Fy,linewidth=1.5,label=r'$F_y$')

plt.xlabel(r'$t$')
plt.ylabel(r'$u_{\tau}$')

plt.tight_layout()
plt.savefig('sod_analysis.png')
#plt.show()
