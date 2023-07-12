import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('t_hist.txt', skiprows=0)

indTime = 1;
nplots = [2, 3];
ylab = [r'$U_b$',r'$(u^2+v^2+w^2)/V$']

fig, axs = plt.subplots(len(nplots), 1, sharex=True, dpi=300)

for i,ind in enumerate(nplots):
  axs[i].plot(data[:,indTime],data[:,ind],linewidth=1.5)
  axs[i].set_xlabel(r'$t$')
  axs[i].set_ylabel(ylab[i])

plt.tight_layout()
plt.savefig('nrs_analysis.png')
#plt.show()
