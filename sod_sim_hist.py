import numpy as np
import matplotlib.pyplot as plt
import sys

fileName = sys.argv[1]
data = np.loadtxt('analysis_'+fileName+'.dat', skiprows=0)

time    = data[:,0];
Ek      = data[:,1];
eps_S   = data[:,2];
eps_D   = data[:,3];
eps_T   = data[:,4];
Mmax    = data[:,5];
ylab = [r'$E_k$',r'$eps_{S}$',r'$eps_{D}$',r'$eps_{T}$',r'$M_{max}$','']

fig, axs = plt.subplots(3, 2, sharex=True, dpi=300)

axs[0,0].plot(time,Ek,linewidth=1.5,label='Ek')
axs[0,1].plot(time,eps_S,linewidth=1.5,label=r'$eps_S$')
axs[1,0].plot(time,eps_D,linewidth=1.5,label=r'$eps_D$')
axs[1,1].plot(time,eps_T,linewidth=1.5,label=r'$eps_T$')
axs[2,0].plot(time,Mmax,linewidth=1.5,label=r'$M_{max}$')

plt.xlabel(r'$t$')

for i,ax in enumerate(axs.flat):
    ax.set(ylabel=ylab[i])
    if(i in [4,5]):
      ax.set(xlabel=r'$t$')
plt.tight_layout()
plt.savefig('sod_analysis.png')
#plt.show()
