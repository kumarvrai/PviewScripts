import os
import sys
import numpy as np
import matplotlib.pyplot as plt

case = sys.argv[1]
#visco    = convert_to_float(sys.argv[2]);
#Lz       = convert_to_float(sys.argv[3]);
#strtTime = convert_to_float(sys.argv[4]);
#aoa      = convert_to_float(sys.argv[5]);

os.system('grep "t-min-max" out.o | sed "/<</d" > t_hist_minmax.txt')
print('--||NRS :: LOADING UVWP DATA')
if('RANS' in case):
  col_to_read=np.arange(1,17,1, dtype=int)
  nplots = [1, 3, 5, 7, 11, 13];
  ylab = [r'$u/U_o$',r'$v/U_o$',r'$w/U_o$',r'$p/\rho U_o^2$', r'$k/U_o^2$', r'$\omega c/U_o$']
else:
  col_to_read=np.arange(1,10,1, dtype=int)
  nplots = [1, 3, 5, 7];
  ylab = [r'$max(U)$',r'$max(V)$',r'$max(W)$',r'$max(P)$']

minmax_data = np.loadtxt('t_hist_minmax.txt', skiprows=100,usecols = col_to_read)

indTime = 0;
time = minmax_data[:,indTime]

print('--||NRS :: PLOTTING DATA')
fig, axs = plt.subplots(int(len(nplots)/2), 2, sharex=True, dpi=300)
ax=axs.flatten()

for i,ind in enumerate(nplots):
  ax[i].plot(time,minmax_data[:,ind],linewidth=1.5)
  ax[i].plot(time,minmax_data[:,ind+1],linewidth=1.5)
  ax[i].set_xlabel(r'$t$')
  ax[i].set_ylabel(ylab[i])

plt.tight_layout()
plt.savefig('nrs_analysis.png')
#plt.show()

if('CLCD' in case):
  os.system('grep "t-cd-cl" out.o | sed "/<</d" > t_hist_minmax.txt')
  print('--||NRS :: LOADING CLCD DATA')
  col_to_read=np.arange(1,8,1, dtype=int)
  nplots = [3, 6];
  ylab = [r'$C_d$',r'$C_l$']
  
  minmax_data = np.loadtxt('t_hist_minmax.txt', skiprows=100,usecols = col_to_read)
  
  indTime = 0;
  time = minmax_data[:,indTime]
  
  print('--||NRS :: PLOTTING DATA')
  fig, axs = plt.subplots(len(nplots), 1, sharex=True, dpi=300)
  ax=axs.flatten()
  
  for i,ind in enumerate(nplots):
    ax[i].plot(time,minmax_data[:,ind],linewidth=1.5)
    ax[i].set_xlabel(r'$t$')
    ax[i].set_ylabel(ylab[i])
  
  plt.tight_layout()
  plt.savefig('nrs_clcd.png')
