import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
import os

#------------------------------------#
def convert_to_float(frac_str):
    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac
#------------------------------------#

fileName = sys.argv[1]
visco    = convert_to_float(sys.argv[2]);
Lz       = convert_to_float(sys.argv[3]);

#-------------SIM CONSTS-------------#
Cp        = 1004.0
Prt       = 0.71
vo        = 1.0
M         = 0.2
delta     = 1.0
rho0      = 1.0
gamma_gas = 1.40
Re        =  50000.0
mul       = (rho0*delta*vo)/Re
Rgas      = Cp*(gamma_gas-1.0)/gamma_gas
to        = vo*vo/(gamma_gas*Rgas*M*M)
po        = rho0*Rgas*to

#-------------------------------------#
file_1 = 'analysis_'+fileName+'.dat';
file_2 = 'surf_code_1-'+fileName+'.dat';

data = np.loadtxt(file_1,skiprows=0)
data_2 = np.loadtxt(file_2,delimiter=",",skiprows=1)

print('--||NEK :: Data Shape ',np.shape(data))
chkCnd = (np.shape(data)[1] > 6)

lw = 1.0

startInd = 10;
time    = data[startInd::,0];
Ek      = data[startInd::,1];
eps_S   = data[startInd::,2];
eps_D   = data[startInd::,3];
eps_T   = data[startInd::,4];
Mmax    = data[startInd::,5];
if(chkCnd):
  maxRho    = data[startInd::,6]
  maxMue    = data[startInd::,7]/visco;
  maxU    = data[startInd::,8];
  maxV    = data[startInd::,9];
  maxW    = data[startInd::,10];
  maxP    = data[startInd::,11];

time_2    = data_2[startInd::,1];
area      = data_2[startInd::,2]
Fpx       = data_2[startInd::,3];
Fpy       = data_2[startInd::,4];
Fpz       = data_2[startInd::,5];
Fvx       = data_2[startInd::,6];
Fvy       = data_2[startInd::,7];
Fvz       = data_2[startInd::,8];


if('pipe' in fileName):
 print("--||SOD :: Calculating utau for pipe")
 utau      = np.sqrt(abs(Fvz/area))
 print("--||SOD :: time-average utau =%f",np.mean(utau,axis=None))
elif('channel' in fileName):
 print("--||SOD :: Calculating utau for channel")
 utau      = np.sqrt(abs(Fvx/area))
 print("--||SOD :: time-average utau = ",np.mean(utau,axis=None))
elif('naca' in fileName):
 print("--||SOD :: Calculating Cl/Cd for NACA")
 print("----||Fvx,Fpx,Fvy,Fpy = ",np.mean(Fvx,axis=None),np.mean(Fpx,axis=None),np.mean(Fvy,axis=None),np.mean(Fpy,axis=None))
 cdp     = 2.0*(Fpx)/Lz
 cdv     = 2.0*(Fvx)/Lz
 cd      = 2.0*(Fvx+Fpx)/Lz
 cl      = 2.0*(Fvy+Fpy)/Lz
 print("--||SOD :: time-average Cdp = ",np.mean(cdp,axis=None))
 print("--||SOD :: time-average Cdv = ",np.mean(cdv,axis=None))
 print("--||SOD :: time-average Cd  = ",np.mean(cd,axis=None))
 print("--||SOD :: time-average Cl  = ",np.mean(cl,axis=None))
else:
 print("--||SOD :: Error. fileName handling not defined!")

if(chkCnd):
  #ylab = [r'$E_k$',r'$\epsilon_{S}$',r'$\epsilon_{D}$',r'$\epsilon_{T}$',r'$M_{max}$','']
  ylab = [r'$E_k$',r'$eps_{T}$',r'$\max(M)$',r'$\max(U)$',r'$\max(\mu_e)/\mu_f$',r'$\max(\rho)$']
  
  #fig, axs = plt.subplots(3, 2, sharex=True, dpi=300)
  fig = plt.figure(dpi=300)
  grid = GridSpec(4, 3, figure=fig)
  
  #------------------#
  axs = fig.add_subplot(grid[0,0])
  axs.plot(time,Ek,linewidth=lw,label='Ek')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[0])
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #plt.gca().set_xticklabels([])
  #------------------#
  axs = fig.add_subplot(grid[1,0])
  #axs.plot(time,eps_S,linewidth=lw,label=r'$eps_S$')
  axs.plot(time,eps_S,linewidth=lw,label=r'$eps_T$')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$eps_S$')
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #plt.gca().set_xticklabels([])
  #------------------#
  axs = fig.add_subplot(grid[2,0])
  #axs.plot(time,eps_S,linewidth=lw,label=r'$eps_S$')
  axs.plot(time,eps_D,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$eps_D$')
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #plt.gca().set_xticklabels([])
  #------------------#
  axs = fig.add_subplot(grid[3,0])
  #axs.plot(time,eps_S,linewidth=lw,label=r'$eps_S$')
  axs.plot(time,eps_T,linewidth=lw,label=r'$eps_T$')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$eps_T$')
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #plt.gca().set_xticklabels([])
  #------------------#
  #------------------#
  axs = fig.add_subplot(grid[0,1])
  #axs.plot(time,eps_D,linewidth=lw,label=r'$eps_D$')
  axs.plot(time,Mmax,linewidth=lw,label=r'$M_{max}$')
  plt.ylabel(r'$M_{max}$')
  plt.xlabel(r'$t$')
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #------------------#
  #------------------#
  axs = fig.add_subplot(grid[1,1])
  axs.plot(time,maxRho,linewidth=lw,label='')
  plt.ylabel(r'$\rho_{max}$')
  plt.xlabel(r'$t$')
  #------------------#
  axs = fig.add_subplot(grid[2,1])
  axs.plot(time,maxP/po,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$max(p)/P_o$')
  plt.title(r'$P_o=%.2f$'%po)
  #------------------#
  if('naca' in fileName):
    axs = fig.add_subplot(grid[3,1])
    axs.plot(time_2,cl,linewidth=lw,label=r'$c_{l}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{l}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))

    axs = fig.add_subplot(grid[3,2])
    axs.plot(time_2,cd,linewidth=lw,label=r'$c_{d}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{d}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  else:
    axs = fig.add_subplot(grid[3,1:2])
    axs.plot(time_2,utau,linewidth=lw,label=r'$u_{\tau}$')
    #axs[2,1].plot(time,Fy,linewidth=lw,label=r'$F_y$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$u_{\tau}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #------------------#
  #------------------#
  axs = fig.add_subplot(grid[0,2])
  axs.plot(time,maxU,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$max(u)/U_o$')
  #------------------#
  axs = fig.add_subplot(grid[1,2])
  axs.plot(time,maxV,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$max(v)/U_o$')
  #------------------#
  axs = fig.add_subplot(grid[2,2])
  axs.plot(time,maxW,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$max(w)/U_o$')
  #------------------#
  
  fig.tight_layout()
  plt.savefig('sod_analysis.png')
  #plt.show()
else:
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
  if('naca' in fileName):
    axs = fig.add_subplot(grid[2,0])
    axs.plot(time_2,cl,linewidth=lw,label=r'$c_{l}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{l}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))

    axs = fig.add_subplot(grid[2,1])
    axs.plot(time_2,cd,linewidth=lw,label=r'$c_{d}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{d}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  else:
    axs = fig.add_subplot(grid[2,:])
    axs.plot(time_2,utau,linewidth=1.5,label=r'$u_{\tau}$')
    #axs[2,1].plot(time,Fy,linewidth=1.5,label=r'$F_y$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$u_{\tau}$')

  plt.tight_layout()
  plt.savefig('sod_analysis.png')
  #plt.show()
