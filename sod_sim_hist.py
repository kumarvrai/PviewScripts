import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
import os
from scipy.interpolate import griddata
#import seaborn as sns
import scipy.signal as signal
from scipy.interpolate import interp1d

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

def calcSpectra(f, nWindows, iPeriodic):
  n = len(f)
  n = n//nWindows
  print('--|| FFT INFO : WINDOW SIGNAL LENGTH={}'.format(n))
  
  wss = 0.
  if iPeriodic:
    print('--||INFO :: SIGNAL PERIODIC, NO WINDOWING')
  else:
    wss = sum(np.hanning(n))
  
  for iWindow in range(0, nWindows, 1):
    #print('data window: {}'.format(iWindow))
    #print('start index: {}, end index: {}'.format(iWindow*n, (iWindow+1)*n-1))
    fx = f[iWindow*n:(iWindow+1)*n]
  
    if iPeriodic:
      fk = np.fft.fft(fx)
    else:
      wx = np.hanning(n)*fx
      fk = np.fft.fft(wx)
    
    if (n % 2) == 0:
      #Mind the odd-ball wavenumber, not considered here
      afk = 2*(np.abs(fk[0:(n//2+1)])/n)
      pfk = 2*(np.abs(fk[0:(n//2+1)])/n)**2
    else:
      afk = 2*(np.abs(fk[0:((n-1)//2+1)])/n)
      pfk = 2*(np.abs(fk[0:((n-1)//2+1)])/n)**2
  
    if iWindow==0:
      afks = afk
      pfks = pfk
    else:
      afks = afks + afk
      pfks = pfks + pfk
  
  afks = afks/nWindows
  pfks = pfks/nWindows
  pfks = pfks/(wss/n)
  
  afks[0] = 0.5*afks[0]
  pfks[0] = 0.5*pfks[0]
  
  return afks, pfks
#------------------------------------#

fileName = sys.argv[1]
visco    = convert_to_float(sys.argv[2]);
Lz       = convert_to_float(sys.argv[3]);
strtTime = convert_to_float(sys.argv[4]);
aoa      = convert_to_float(sys.argv[5]);

calc_psd=1;
mode='FFT'
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
#po        = 1.0

#-------------------------------------#
file_1 = 'analysis_'+fileName+'.dat';
file_2 = 'surf_code_1-'+fileName+'.dat';

data = np.loadtxt(file_1,skiprows=0)
data_2 = np.loadtxt(file_2,delimiter=",",skiprows=1)

print('--||NEK :: Data Shape ',np.shape(data))
chkCnd = (np.shape(data)[1] > 6)

lw = 1.0


startInd = np.argmin(abs(data[:,0]-strtTime),axis=None);
print('--||NEK :: Plotting Shape ',np.shape(data[startInd::,0]))
time    = data[startInd::,0];
Ek      = data[startInd::,1];
eps_S   = data[startInd::,2];
eps_D   = data[startInd::,3];
eps_T   = data[startInd::,4];
Mmax    = data[startInd::,5];
if(chkCnd):
  maxRho  = data[startInd::,6]
  maxMue  = data[startInd::,7]/visco;
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

Fpxx      = Fpx*np.cos(aoa*np.pi/180.0)+Fpy*np.sin(aoa*np.pi/180.0)
Fpyy      = -Fpx*np.sin(aoa*np.pi/180.0)+Fpy*np.cos(aoa*np.pi/180.0)
Fvxx      = Fvx*np.cos(aoa*np.pi/180.0)+Fvy*np.sin(aoa*np.pi/180.0)
Fvyy      = -Fvx*np.sin(aoa*np.pi/180.0)+Fvy*np.cos(aoa*np.pi/180.0)

Fpx = Fpxx
Fpy = Fpyy
Fvx = Fvxx
Fvy = Fvyy


if('pipe' in fileName):
 print("--||SOD :: Calculating utau for pipe")
 utau      = np.sqrt(abs(Fvz/area))
 print("--||SOD :: time-average utau = ",np.mean(utau,axis=None))
elif('channel' in fileName):
 print("--||SOD :: Calculating utau for channel")
 utau      = np.sqrt(abs(Fvx/area))
 print("--||SOD :: time-average utau = ",np.mean(utau,axis=None))
elif('naca' in fileName):
 print("--||SOD :: Calculating Cl/Cd for NACA")
 print("--||SOD :: Calculating for AoA = ",aoa)
 #print("----||Fvx,Fpx,Fvy,Fpy = ",np.mean(Fvx,axis=None),np.mean(Fpx,axis=None),np.mean(Fvy,axis=None),np.mean(Fpy,axis=None))
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
  ##------------------#
  #axs = fig.add_subplot(grid[1,0])
  ##axs.plot(time,eps_S,linewidth=lw,label=r'$eps_S$')
  #axs.plot(time,eps_S,linewidth=lw,label=r'$eps_T$')
  #plt.xlabel(r'$t$')
  #plt.ylabel(r'$eps_S$')
  #plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  ##plt.gca().set_xticklabels([])
  ##------------------#
  axs = fig.add_subplot(grid[1,0])
  #axs.plot(time,eps_S,linewidth=lw,label=r'$eps_S$')
  axs.plot(time,eps_D,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$eps_D$')
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #plt.gca().set_xticklabels([])
  #------------------#
  axs = fig.add_subplot(grid[2,0])
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
    po = 5;
    print("--||SOD :: BEST-FIT PORDER = ",po)
    pFit = np.poly1d(np.polyfit(time_2, cl, po))
    axs = fig.add_subplot(grid[3,0])
    axs.plot(time_2,cl,linewidth=lw,label=r'$c_{l}$')
    plt.plot(time_2,pFit(time_2),'k--',linewidth=lw)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{l}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))

    pFit = np.poly1d(np.polyfit(time_2, cd, po))
    axs = fig.add_subplot(grid[3,1])
    axs.plot(time_2,cd,linewidth=lw,label=r'$c_{d}$')
    plt.plot(time_2, pFit(time_2),'k--',linewidth=lw)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{d}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  else:
    axs = fig.add_subplot(grid[3,:-1])
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
  axs = fig.add_subplot(grid[3,2])
  axs.plot(time,maxMue,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(r'$max(\mu_e)/\mu_f$')
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
    po = 5;
    print("--||SOD :: BEST-FIT PORDER = ",po)
    pFit = np.poly1d(np.polyfit(time_2, cl, po))
    axs = fig.add_subplot(grid[2,0])
    axs.plot(time_2,cl,linewidth=lw,label=r'$c_{l}$')
    plt.plot(time_2,pFit(time_2),'k--',linewidth=lw)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{l}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))

    pFit = np.poly1d(np.polyfit(time_2, cd, po))
    axs = fig.add_subplot(grid[2,1])
    axs.plot(time_2,cd,linewidth=lw,label=r'$c_{d}$')
    plt.plot(time_2, pFit(time_2),'k--',linewidth=lw)
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

if(calc_psd):
  print('--|| NEK :: CALCULATE POWER SPECTRA (PSD).')
  if('naca' in fileName):
    scale=np.sin(aoa*np.pi/180.0);
    cd_data = cd
  else:  
    scale=1.0;
    cd_data = utau**2
  indx = range(0,len(cd_data),1)
  tSignal=time_2[indx]
  ySignal=cd_data[indx]
  print('--|| INFO :: Total samples',len(ySignal))
  freq_min = 0.01; freq_max = 1.0;
  if('LOMB' in mode):
    nout = 10000;
    f = np.linspace(freq_min, freq_max, nout)
    pgram = signal.lombscargle(tSignal, ySignal, f)
  elif('FFT' in mode):
    delta_t = np.amax(np.diff(tSignal),axis=None)
    s_rate = 1.0/delta_t
    print('----|| INFO :: SAMPLING FREQ = %.2f'% s_rate)
    Nt = int((np.amax(tSignal,axis=None)-np.amin(tSignal,axis=None))/delta_t)+1
    print('----|| INFO :: INTERPOLATING ON %d POINTS'%(Nt))
    print('--|| ALYA :: INTERPOLATE CONSTANT TIME.')
    fInterp = interp1d(tSignal,ySignal);
    tSignal = np.linspace(np.amin(tSignal,axis=None),np.amax(tSignal,axis=None),Nt);
    ySignal = fInterp(tSignal)
    hann_size = 1    
    amp,power = calcSpectra(ySignal, hann_size, 0)
    pgram=amp
    #pgram = np.fft.fft(ySignal)
    N = len(pgram); print('--|| INFO :: FFT LENGTH = %d'%N) 
    n = np.arange(N); 
    T = float(N)/s_rate; 
    f = 0.5*n/T
    #plt.stem(f, np.abs(pgram), 'b', \
    # markerfmt=" ", basefmt="-b")
    #plt.ylabel(r'$P_{u^\prime u^\prime}$')
    #ax.set(xlim=(0, 10))
  else:
    raise ValueError('ALYA ERROR :: METHOD TO CALC PSD NOT SPECIFIED!')
  
  #---plot PSD----#
  #indx=np.where(np.logical_and(f>0.001,f<1))
  fig=plt.figure(6, figsize=(8, 6))
  ax=plt.subplot(111)
  ax.plot(scale*f, pgram, 'k',linewidth=1.0)
  
  plt.xlabel(r'$St$')
  plt.ylabel(r'$PSD(C_f)$')
  plt.tight_layout()
  ax.set_yscale("log")
  ax.set_xscale("log")
  #ax.set_xlim(left=2e-3,right=8e-1)
  #ax.set_ylim(bottom=1e-4,top=1e1)
  for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
  #fig.patch.set_alpha(0.5)  
  #plt.savefig('fft_hist.png', transparent=False)
  plt.savefig('sod_fft_hist.png')
