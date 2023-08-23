import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys

fileName = sys.argv[1]
data = np.loadtxt('analysis_'+fileName+'.dat',skiprows=0)
data_2 = np.loadtxt('surf_code_1-'+fileName+'.dat',delimiter=",",skiprows=1)
data_3 = np.loadtxt('monitor_'+fileName+'.dat',skiprows=0)

#visco   = 5.3566e-05;
visco   = 1.0;
lw = 1.0

startInd = 100;
time    = data[startInd::,0];
Ek      = data[startInd::,1];
eps_S   = data[startInd::,2];
eps_D   = data[startInd::,3];
eps_T   = data[startInd::,4];
Mmax    = data[startInd::,5];

time_2    = data_2[startInd::,1];
area      = data_2[startInd::,2]
Fpx       = data_2[startInd::,3];
Fpy       = data_2[startInd::,4];
Fvx       = data_2[startInd::,6];
Fvy       = data_2[startInd::,7];

time_3    = data_3[startInd::,0];
maxRho    = data_3[startInd::,1]
maxMue    = data_3[startInd::,2]/visco;
maxVel    = data_3[startInd::,3];

utau      = np.sqrt(abs(Fvx/area))

#ylab = [r'$E_k$',r'$\epsilon_{S}$',r'$\epsilon_{D}$',r'$\epsilon_{T}$',r'$M_{max}$','']
ylab = [r'$E_k$',r'$eps_{T}$',r'$\max(M)$',r'$\max(U)$',r'$\max(\mu_e)/\mu_f$',r'$\max(\rho)$']

#fig, axs = plt.subplots(3, 2, sharex=True, dpi=300)
fig = plt.figure(dpi=300)
grid = GridSpec(4, 2, figure=fig)

#------------------#
axs = fig.add_subplot(grid[0,0])
axs.plot(time,Ek,linewidth=lw,label='Ek')
plt.xlabel(r'$t$')
plt.ylabel(ylab[0])
plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
#plt.gca().set_xticklabels([])
#------------------#
#------------------#
axs = fig.add_subplot(grid[0,1])
axs.plot(time_3,maxVel,linewidth=lw,label='')
plt.xlabel(r'$t$')
plt.ylabel(ylab[3])
#plt.gca().set_xticklabels([])
#------------------#
axs = fig.add_subplot(grid[1,0])
#axs.plot(time,eps_S,linewidth=lw,label=r'$eps_S$')
axs.plot(time,eps_T,linewidth=lw,label=r'$eps_T$')
plt.xlabel(r'$t$')
plt.ylabel(ylab[1])
plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
#plt.gca().set_xticklabels([])
#------------------#
#------------------#
axs = fig.add_subplot(grid[1,1])
axs.plot(time_3,maxMue,linewidth=lw,label='')
plt.xlabel(r'$t$')
plt.ylabel(ylab[4])
plt.gca().set_ylim(bottom=0)
#plt.gca().set_xticklabels([])
#------------------#
axs = fig.add_subplot(grid[2,0])
#axs.plot(time,eps_D,linewidth=lw,label=r'$eps_D$')
axs.plot(time,Mmax,linewidth=lw,label=r'$M_{max}$')
plt.xlabel(r'$t$')
plt.ylabel(ylab[2])
plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
#------------------#
#------------------#
axs = fig.add_subplot(grid[2,1])
axs.plot(time_3,maxRho,linewidth=lw,label='')
plt.xlabel(r'$t$')
plt.ylabel(ylab[5])
#------------------#
axs = fig.add_subplot(grid[3,:])
axs.plot(time_2,utau,linewidth=lw,label=r'$u_{\tau}$')
#axs[2,1].plot(time,Fy,linewidth=lw,label=r'$F_y$')
plt.xlabel(r'$t$')
plt.ylabel(r'$u_{\tau}$')
plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
#------------------#

fig.tight_layout()
plt.savefig('sod_analysis.png')
#plt.show()
