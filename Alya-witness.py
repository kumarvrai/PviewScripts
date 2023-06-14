import os
import time
import sys

#os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'
import contextlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import CFDlib as cfd
from os.path import expanduser
from scipy.interpolate import interp1d

#######################################################################
@contextlib.contextmanager
def plot_kde_as_log(base=10, support_threshold=1e-4):
    """Context manager to render density estimates on a logarithmic scale.
    Usage:

        with plot_kde_as_log():
          sns.jointplot(x='x', y='y', data=df, kind='kde')
    """
    #old_stats = sns.distributions._has_statsmodels
    #old_univar = sns.distributions._scipy_univariate_kde
    #old_bivar = sns.distributions._scipy_bivariate_kde
    old_call = sns._statistics.KDE.__call__

    #sns.distributions._has_statsmodels = False
    def log_clip_fn(v):
      #global n_total
      #v *= n_total; 
      #v /= np.amax(v,axis=None) 
      v = np.log10(np.clip(v, support_threshold, np.inf))
      v -= np.log10(support_threshold)
      v /= np.log10(base)
      return v
    def new_call(*args, **kwargs):
      density, support = old_call(*args, **kwargs)
      density = log_clip_fn(density)
      return density, support
    #def new_univar(*args, **kwargs):
    #  x, y = old_univar(*args, **kwargs)
    #  y = log_clip_fn(y)
    #  return x, y
    #def new_bivar(*args, **kwargs):
    #  x, y, z = old_bivar(*args, **kwargs)
    #  z = log_clip_fn(z)
    #  return x, y, z
    
    #sns.distributions._scipy_univariate_kde = new_univar
    #sns.distributions._scipy_bivariate_kde = new_bivar
    sns._statistics.KDE.__call__ = new_call

    try:
      yield
    finally:
      sns._statistics.KDE.__call__ = old_call
      #sns.distributions._has_statsmodels = old_stats
      #sns.distributions._scipy_univariate_kde = old_univar
      #sns.distributions._scipy_bivariate_kde = old_bivar

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

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    indMin = np.argmin(np.sqrt(dist_2),axis=None)
    return indMin, np.amin(np.sqrt(dist_2),axis=None)


#######################################################################

home = expanduser("~")
casePath = os.getcwd()  


ls = ['-','--','-.',':','-o','-v','-^','-s','-p','-d','-*','-+','-x']
linestyles = [
     ('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]     
mpl.rcParams['mathtext.fontset'] = 'cm' 
mpl.rcParams['mathtext.rm'] = 'serif'

caseName	= sys.argv[1]
airfoil 	= sys.argv[2]
side 		= sys.argv[3];
mode 		= sys.argv[4]

#indP = [2,3,4]
indP = [0,1,2,3,4]
#indP = [0,2,6]
#indP = [8,10,12,14]
#indP = [12,23,28]
#indP = [28,29,30]
#indP = [4,5,6,7]
#indP = [0,1,2,3,4,5,6]
#indP = [11,12,13]
#indP = [22,23]
#indP = [5,6,7,8,9]
#indP = [2,3,4,5,6,7,8,9,10,11]
#indP = (np.arange(23,30)-1)
#indP = (np.arange(1,7)-1)
#indP = np.array([1,3,5,7,8,10])-1
#indP = np.array([2,6,10])-1
#indP = np.array([2,4,6,8])-1
#indP = [581,584] #SLE-5 (0.6,0.9) to compare 
#indP = [288,290,291] #TBL-5 (0.6,0.9) to compare 
#indP = [480,481,482,483] # RLE2 - TE-SP (0.2,0.3,0.4,0.5) location
#indP = [484,485,486,487] # RLE2 - TE-SP (0.6,0.7,0.8,0.9) location
#indP = [424,464]   #RLE1-5 (FP,SP) at x/c=0.2
#indP = [202,442]  #RLE1-5 (FP,SP) at x/c=0.6
#indP = [235,427] #RLE1-5 (FP,SP) at x/c=0.9
#indP = [236,468] #RLE1-5 TE (FP,SP) at x/c=0.95
#indP = [235,467] #RLE1-5 TE (FP,SP) at x/c=0.9
#indP = [204,96]  #RLE1-10 LE (FP,SP) at x/c=0.2

#indP = [440,384] #RLE2-5 (FP,SP) at x/c=0.2
#indP = [450,346] #RLE2-5 (FP,SP) at x/c=0.6
#indP = [555,331] #RLE2-5 (FP,SP) at x/c=0.9
#indP = [580,292] #RLE2-5 TE (FP,SP) at x/c=0.95
#indP = [869,437]  #RLE2-10 TE (FP,SP) at x/c=0.7
#indP = [481,865]  #RLE2-10 LE (FP,SP) at x/c=0.3
#indP = [480,360]  #RLE2-10 LE (FP,SP) at x/c=0.2

#-------FP-SP----------#
# Matching with beta plots
# RLE1-5 F @ z/c=0.04; S @ z/c=0.021
#indP = [232,234,235,120,122,123] #RLE1-5 [F,F,F,S,S,S] 

x_loc = [1.2, 1.5]

# LOAD AIRFOIL SPECIFIC FILES
print('--|| ALYA INITIALIZING')
plt.close("all")
lw = 1.2;
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

if('0012' in airfoil):
  fAirU = home+'/1.post_process/1.airfoil/3.PviewScripts/1.0012-Slices/naca0012-UP.txt'
  fAirL = home+'/1.post_process/1.airfoil/3.PviewScripts/1.0012-Slices/naca0012-DOWN.txt'
elif('4412' in airfoil):  
  fAirU = home+'/1.post_process/1.airfoil/3.PviewScripts/2.4412-Slices/naca4412-UP.txt'
  fAirL = home+'/1.post_process/1.airfoil/3.PviewScripts/2.4412-Slices/naca4412-DOWN.txt'
else:
  sys.exit("--||Alya (ERROR): AIRFOIL NOT IN THE LIST")

if('ss' in side):
  coordAir = np.loadtxt(fAirU, delimiter=',');
else:
  coordAir = np.loadtxt(fAirL, delimiter=',')

# INTERPOLATION ON SUCTION SIDE
thAir = np.arctan2(np.diff(coordAir[:,1]),np.diff(coordAir[:,0]))
airLen = len(coordAir)
coordMid = 0.5*(coordAir[0:airLen-1,0]+coordAir[1:airLen,0])
Fth = interp1d(coordMid,thAir)

print('--|| ALYA :: READING AIRFOIL DATA.')
stime = time.time()
coordAirU = np.loadtxt(fAirU, delimiter=',')
coordAirL = np.loadtxt(fAirL, delimiter=',')
print('--|| ALYA :: DONE. TIME=',time.time()-stime)

print('--|| ALYA :: READING WITNESS POINTS.')
stime = time.time()
witPoints = np.loadtxt('./witness-nastin.dat')
print('--|| ALYA :: DONE. TIME=',time.time()-stime)

nWit = len(witPoints)
z = np.unique(witPoints[:,-1]);
nz = int(nWit/len(z))
indPLabel = range(len(indP));
print('--|| ALYA :: PROCESSING', nWit, 'TOTAL WITNESS POINTS')
print('--|| ALYA :: NUMBER OF WIT POINTS IN A PLANE ',nz,'WITH',len(z),'PLANES')

n = len(indP)
color = plt.cm.rainbow(np.linspace(0.0,1.0,n)) # This returns RGBA; convert:
hexcolor = map(lambda rgb:'#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255),
               tuple(color[:,0:-1]))
#mpl.rcParams['axes.color_cycle'] = hexcolor

if("SHOWPTS" in mode):
  print('--|| ALYA :: PLOTTING AIRFOIL WITNESS POINTS ONLY.')
  fig = plt.figure()
  grid = plt.GridSpec(1, 3, wspace=0.3, hspace=0.1)
  axs = plt.subplot(grid[0, 0:]);
  plt.plot(coordAirU[:,0],coordAirU[:,1],'k',linewidth=2)
  plt.plot(coordAirL[:,0],coordAirL[:,1],'k',linewidth=2)
  leg = []
  #for n in range(nz):
   #lab = r'$P_{'+str(n+1)+'}$';
  for n in range(len(indP)):
   pt = indP[n]
   lab = r'$P_{'+str(n+1)+'}$';
   #if(n<11):
   # sgnX = -1
   # sgnY = 1 
   #else:
   sgnX = 1
   sgnY = 1 
   #if n in [0,1,2,11,12,13]:
   # sgnY = 0
   leg.append(lab)
   plt.plot(witPoints[pt,0], witPoints[pt,1],'ko',linewidth=0.5)
   plt.text(witPoints[pt,0]+(sgnX)*0.025, witPoints[pt,1]+(sgnY)*0.025, lab, fontsize=MEDIUM_SIZE)
  #axs.legend()
  #fig.legend(handles,     # The line objects
  #             labels,
  #             loc="upper right",   # Position of legend	i
  #             bbox_to_anchor=(0.91,0.52),
  #             ncol = 2, fancybox=True,shadow=True
  #            )
  axs.set_anchor('W')
  plt.xlabel(r'$x/c$');
  plt.yticks([])
  plt.xticks([0,0.5,1,1.3])
  axs.set_aspect(1.5);
  axs.set(xlim=(-0.01, 1.3))
  axs.spines['top'].set_visible(False)
  axs.spines['right'].set_visible(False)
  #axs.spines['bottom'].set_visible(False)
  axs.spines['left'].set_visible(False)
  plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True,
        labeltop=False) # labels along the bottom edge are off
else:  
  
  if(any(string in mode for string in ["PSD","THIS"])):
    if("UU" in mode):
      calcVar = 'U'
    elif("VV" in mode):
      calcVar = 'V'
    elif("WW" in mode):
      calcVar = 'W'
    elif("PP" in mode):
      calcVar = 'P'
    else:
      raise ValueError('ALYA ERROR :: CALCVAR NOT SPECIFIED!')
    print('--|| INFO :: PSD CALCULATIONS ON %s'% calcVar)
    print('--|| ALYA :: READING WITNESS DATA.')
    stime = time.time()
    if('U' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./%s-VELOX.wit.bin'%(caseName))
    elif('V' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./%s-VELOY.wit.bin'%(caseName))
    elif('W' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./%s-VELOZ.wit.bin'%(caseName))
    elif('P' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./%s-PRESS.wit.bin'%(caseName))
    else:
      raise ValueError('ALYA ERROR :: VARIABLE NOT SPECIFIED!')
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)

    print('----|| INFO :: %d TIME VALUES RANGE %.2f-%.2f'%(len(wit_data),time_data[0], time_data[-1]))
  

    import scipy.signal as signal
    print('--|| ALYA :: ARRANGING ARRAYS MONOTONICALLY.')
    stime = time.time()
    tM = np.maximum.accumulate(time_data)
    t, ind_unq = np.unique(tM, return_index=True)
    print('----|| INFO :: UNIQUE TIME ARRAY SIZE=',np.size(t))
    totalMarkPts = 1000;
    markFreq = int(np.amax((np.floor(len(t)/totalMarkPts),1),axis=None))
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)
    fig, axs = plt.subplots(len(indP), 2, sharex=True, figsize=(11.69,8.27))
    #ax = plt.subplot(len(indP), 1, i+1)
    ax = plt.gca()
    for (i,j) in enumerate(indP):
      print('----|| INFO :: PROCESSING POINT %d, x,y,z='%j, witPoints[j,:])
      #lab = r'$P_{'+str(j+1)+'}$';
      #lab = r'$x/c=%.3f, \, y/c=%.2f$'%(witPoints[j,0],witPoints[j,1]);
      lab = r'$x/c=%.1f$'%(witPoints[j,0]);
      if('ZAVG' in mode):
        print('----|| INFO :: PERFORMING SPANWISE AVERAGING')
        z_per_index = np.arange(j,nWit,nz)
        #for idx in z_per_index:
          #print('----|| X=%.2f, Y=%.3f, Z=%.3f '% (witPoints[idx,0],witPoints[idx,1],witPoints[idx,2]));
        dataPoint = np.mean(wit_data[:,z_per_index],axis=1)
      else:
        dataPoint = wit_data[:,j]; #Extract data
      print('----|| INFO :: WITNESS ARRAY SIZE=',np.shape(dataPoint))
      backFlow = 100.0*float(len(np.where(dataPoint<0.0)[0]))/len(dataPoint)
      print('----|| INFO :: PERCENT BACKFLOW=', backFlow)
      dataPoint = dataPoint[ind_unq];
      print('----|| INFO :: UNIQUE WITNESS ARRAY SIZE=',np.shape(dataPoint))
      tSignal = t;
      ySignal = dataPoint-np.mean(dataPoint,axis=None)
      if("THIS" in mode):
        print('--|| ALYA :: CALCULATE TIME HISTORY (THIS).')
        N = 1000;
        y1 = dataPoint - np.convolve(dataPoint, np.ones(N)/N, mode='same')
        #ax = plt.subplot(211)
        plt.plot(tSignal, dataPoint, '.', markersize=2)
        plt.plot(tSignal, np.mean(dataPoint,axis=None)*np.ones(len(tSignal)), 'k'+ls[i], linewidth=1.75)
        plt.plot(tSignal, 0.0*dataPoint, 'k--', linewidth=1.5)
        plt.tight_layout()
        if(i==len(indP)-1):
         plt.xlabel(r'$tU_o/c$');
         ylab = r'$%s/U_o$'% calcVar.lower()
         plt.ylabel(ylab);
        plt.tight_layout()
        plt.suptitle(lab)
      elif("PSD" in mode):
        print('--|| ALYA :: CALCULATE POWER SPECTRA (PSD).')
        if('LOMB' in mode):
          nout = 10000;
          freq_min = 0.01; freq_max = 1000.0;
          f = np.linspace(freq_min, freq_max, nout)
          pgram = signal.lombscargle(tSignal, ySignal, f)
        elif('FFT' in mode):
          delta_t = np.amax(np.diff(tSignal),axis=None)
          s_rate = 1.0/delta_t
          print('----|| INFO :: SAMPLING FREQ = %.2f'% s_rate)
          if('INTERP' in mode):
            Nt = int((np.amax(tSignal,axis=None)-np.amin(tSignal,axis=None))/delta_t)+1
            print('----|| INFO :: INTERPOLATING ON %d POINTS'%(Nt))
            print('--|| ALYA :: INTERPOLATE CONSTANT TIME.')
            fInterp = interp1d(tSignal,ySignal);
            tSignal = np.linspace(np.amin(tSignal,axis=None),np.amax(tSignal,axis=None),Nt);
            ySignal = fInterp(tSignal)
            print('--|| ALYA :: DONE. TIME=',time.time()-stime)
          hann_size = 1    
          amp,pgram = calcSpectra(ySignal, hann_size, 0)
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
        #---Subplot 1----#
        ax1 = axs[i,0]; ax2 = axs[i,1];
        ax1.plot(f, pgram, markevery=markFreq,linewidth=1.0,label=lab)
        if('U' in calcVar):
          ax1.set_ylabel(r'$E_{u^\prime u^\prime}$')
        elif('V' in calcVar):
          ax1.set_ylabel(r'$E_{v^\prime v^\prime}$')
        elif('W' in calcVar):
          ax1.set_ylabel(r'$E_{w^\prime w^\prime}$')
        elif('P' in calcVar):
          ax1.set_ylabel(r'$E_{p^\prime p^\prime}$')
        ax1.set_yscale("log")
        ax1.set_xscale("log")
        ax1.set_xlim(0.01, 100)
        ax1.set_ylim(bottom=1e-12)
        #ax1.legend(loc="lower left",fontsize='small',ncol=2)
        #---Subplot 2----#
        ax2.plot(f, np.multiply(f,pgram), markevery=markFreq,linewidth=1.0,label=lab)
        if('U' in calcVar):
          ax2.set_ylabel(r'$f E_{u^\prime u^\prime}$')
        elif('V' in calcVar):
          ax2.set_ylabel(r'$f E_{v^\prime v^\prime}$')
        elif('W' in calcVar):
          ax2.set_ylabel(r'$f E_{w^\prime w^\prime}$')
        elif('P' in calcVar):
          ax2.set_ylabel(r'$f E_{p^\prime p^\prime}$')
        #plt.yscale("log")
        ax2.set_xscale("log")
        ax2.set_xlim(0.01, 100)
        #plt.gca().set_xlim(right=1000);
        #ax.annotate(lab, xy=(0.05, 1.1), xycoords='axes fraction', fontsize=MEDIUM_SIZE,
        #                horizontalalignment='left', verticalalignment='top')
        #ax.legend(loc="lower left")
        #fig.legend(handles,     # The line objects
        #             labels,
        #             loc="upper right",   # Position of legend	i
        #             bbox_to_anchor=(0.91,0.52),
        #             ncol = 2, fancybox=True,shadow=True
        #            )
        ax1.plot(f[np.where(f>100)],np.power(f[np.where(f>100)],-5.0/3),'k--')
        if(i==len(indP)):
         ax1.set_xlabel(r'$fc/U_o$');
         ax2.set_xlabel(r'$fc/U_o$');
        #ax1.ticklabel_format(axis='both', style='sci')
        ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        ## Hide x labels and tick labels for top plots and y ticks for right plots.
        #for ax in axs.flat:
        #  ax.label_outer()
      else:
        raise ValueError('ALYA ERROR :: CHECK MODE THAT YOU SPECIFIED!')

      ##--------TWEAK THE SUBPLOTS-----------#
      #if(i==0):
      # # Hide the right and top spines
      # ax.spines['top'].set_visible(False)
      # ax.spines['bottom'].set_visible(True)
      # ax.tick_params(
      #  axis='x',          # changes apply to the x-axis
      #  which='both',      # both major and minor ticks are affected
      #  bottom=True,      # ticks along the bottom edge are off
      #  top=False,         # ticks along the top edge are off
      #  labelbottom=True) # labels along the bottom edge are off
      #elif(i==len(indP)-1):
      # ax.spines['top'].set_visible(False)
      # ax.spines['bottom'].set_visible(True)
      # # Only show ticks on the left and bottom spines
      # ax.tick_params(
      #  axis='x',          # changes apply to the x-axis
      #  which='both',      # both major and minor ticks are affected
      #  bottom=True,      # ticks along the bottom edge are off
      #  top=False,         # ticks along the top edge are off
      #  labelbottom=True) # labels along the bottom edge are off
      #else:
      # ax.spines['bottom'].set_visible(False)
      # ax.spines['top'].set_visible(False)
      # ax.tick_params(
      #  axis='x',          # changes apply to the x-axis
      #  which='both',      # both major and minor ticks are affected
      #  bottom=False,      # ticks along the bottom edge are off
      #  top=False,         # ticks along the top edge are off
      #  labelbottom=False) # labels along the bottom edge are off
      # 
      #ax.spines['right'].set_visible(False)
      #ax.tick_params(
      #  axis='y',          # changes apply to the x-axis
      #  which='minor',      # both major and minor ticks are affected
      #  left=False,      # ticks along the bottom edge are off
      #  right=False)         # ticks along the top edge are off
      #ax.tick_params(
      #  axis='y',          # changes apply to the x-axis
      #  which='major',      # both major and minor ticks are affected
      #  right=False)         # ticks along the top edge are off

      ## plt.locator_params(axis='y', nbins=3)
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)

     
  elif("TPCORR" in mode):
    if("UU" in mode):
      calcVar = 'UU'
    elif("VV" in mode):
      calcVar = 'VV'
    elif("WW" in mode):
      calcVar = 'WW'
    elif("PP" in mode):
      calcVar = 'PP'
    elif("UV" in mode):
      calcVar = 'UV'
    else:
      raise ValueError('ALYA ERROR :: TPC VARIABLE NOT SPECIFIED!')
      
    print('--|| INFO :: TWO-POINT CORRELATION ON %s'% calcVar)
    print('--|| ALYA :: READING WITNESS DATA.')
    stime = time.time()
    if('UU' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./%s-VELOX.wit.bin'%(caseName))
      wit_data_1 = wit_data_0;
    elif('VV' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./%s-VELOY.wit.bin'%(caseName))
      wit_data_1 = wit_data_0;
    elif('WW' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./%s-VELOZ.wit.bin'%(caseName))
      wit_data_1 = wit_data_0;
    elif('PP' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./%s-PRESS.wit.bin'%(caseName))
      wit_data_1 = wit_data_0;
    elif('UV' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./%s-VELOX.wit.bin'%(caseName))
      time_data, wit_data_1 = cfd.ExportReadProbe('./%s-VELOY.wit.bin'%(caseName))
    else:
      raise ValueError('ALYA ERROR :: VARIABLE NOT SPECIFIED!')
    print('----|| INFO :: %d TIME VALUES RANGE %.2f-%.2f'%(len(time_data),time_data[0], time_data[-1]))
    print('----|| INFO :: LOADED %d COLUMNS'%(len(wit_data_0[0,:])))

    print('--|| ALYA :: DONE. TIME=',time.time()-stime)
    markFreq = int(np.amax((np.floor(len(z)/20),1),axis=None))
    print('--|| ALYA :: MARKER FREQ=',markFreq)
    print('--|| ALYA :: CALCULATING WIT POINT CORRELATIONS.')
    stime = time.time()
    uu = np.zeros((len(z),len(indP)),dtype=float)
    for j in range(len(indP)):
     pair1 = indP[j]
     u0 = wit_data_0[:,pair1]; u0 = u0-np.mean(u0);
     if('RAVG' in mode):
       #Ns = int(np.ceil(float(len(time_data))/10));
       Ns = 2000;
       u0 = u0 - np.convolve(u0, np.ones(Ns)/Ns, mode='same')
     for i in range(len(z)):
      pair2 = i*nz+pair1;
      # print(pair1,pair2);
      # print('----|| X=%.2f, Y=%.3f, Z=%.3f '% (witPoints[pair2,0],witPoints[pair2,1],witPoints[pair2,2]));
      # Calculate the correlation b/w 2 signals
      u1 = wit_data_1[:,pair2]; u1 = u1-np.mean(u1);
      if('RAVG' in mode):
        u1 = u1 - np.convolve(u1, np.ones(Ns)/Ns, mode='same')
      prod = np.multiply(u0,u1)
      uu[i,j] = np.mean(prod)/(np.sqrt(np.mean(u0*u0))*np.sqrt(np.mean(u1*u1)))
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)

    print('--|| ALYA :: PLOTTING RESULTS.')
    
    fig = plt.figure()
    grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.8)
    
    axs = plt.subplot(grid[0, 0:]);
    for i in range(len(indP)):
      lab = r"$P_{"+str(i+1)+"}$";
      print("--|| ALYA :: PROCESSING POINT",lab);
      print('----|| X=%.2f, Y=%.3f'% (witPoints[i,0],witPoints[i,1]))
      lname, lcode = linestyles[j][0:]
      if('POLY' in mode):
        order = 15
        uu_avg = uu[:,i]
        uu_avg = np.polyval(np.polyfit(z,uu_avg,order),z)
      else:
        uu_avg = uu[:,i]
      l1, = plt.plot(z, uu_avg, ls[i], linewidth=lw, color='black',markersize=4,
            markevery=markFreq, markerfacecolor='black', markeredgecolor='black', label=lab)
    plt.axhline(y=0.0, color='g', linestyle='-', linewidth=0.5)
    plt.xlabel(r'$z/c$');
    labStr = r'$\mathcal{R}_{\langle %s^\prime %s^\prime\rangle}$'%(mode[7].lower(),mode[8].lower())
    plt.ylabel(labStr)
    plt.xticks([0,0.05,0.1])
    #plt.yticks([0,0.5,1])
    axs.set(xlim=(0,np.amax(z)))
    
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=True,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=True,
        labelright=False) # labels along the bottom edge are off
    axs.legend(loc="upper right",   # Position of legend
               bbox_to_anchor=(0.98,0.98),
               ncol = 2, frameon=True, edgecolor='black'
              )
  elif("BACKFLOW" in mode):
    import pandas as pd
    import seaborn as sns
    global n_total
    dcolor = ['Reds', 'Blues_r', 'Greys_r', 'Blues', 'Greens']
    lcolor = ['r', 'b', 'k', 'g']
    bwidth = 2.0; barScale=0; # band width(controls smoothing), log-scaling
    print('--|| INFO :: BANDWIDTH=%.1f, LOGSCALE=%d'%(bwidth,barScale))
    print('--|| INFO :: BACK-FLOW CORRELATIONS')
    print('--|| ALYA :: READING WITNESS DATA.')
    stime = time.time()
    time_data, wit_data_0 = cfd.ExportReadProbe('./%s-VELOX.wit.bin'%(caseName))
    time_data, wit_data_1 = cfd.ExportReadProbe('./%s-VELOY.wit.bin'%(caseName))
    #rows,clms = (len(indP)//2,2)
    #fig, axs = plt.subplots(nrows=rows,ncols=clms)
    #axs = axs.reshape(-1)
    for j,i in enumerate(indP):
      x_wit,y_wit,z_wit = witPoints[i,:]
      if(np.logical_and(x_wit>0.0,x_wit<1.0)):
        indMin,pdist = closest_node((x_wit,y_wit), coordAir);
        theta = Fth(coordAir[indMin,0])
      else:
        theta=0.0
      lab = "$P_{"+str(i+1)+"}$";
      print("--|| ALYA :: POINT (%.2f,%.2f,%.3f) AT d=%.3f"% (x_wit,y_wit,z_wit,pdist));
      u0 = wit_data_0[:,i]; u_avg = np.mean(u0,axis=None); 
      v0 = wit_data_1[:,i]; v_avg = np.mean(v0,axis=None);
      #print("--|| ALYA :: TH=%.2f, U=%.2f, V=%.2f"% (theta*180.0/np.pi,u_avg,v_avg));
      ut = (u0*np.cos(theta)+v0*np.sin(theta)); u_avg = np.mean(ut,axis=None);
      vn = (-u0*np.sin(theta)+v0*np.cos(theta)); v_avg = np.mean(vn,axis=None);
      # print("--|| ALYA :: Ut=%.2f, Vn=%.2f"% (u_avg,v_avg));
      ut_p = ut-np.mean(ut,axis=None)
      vn_p = vn-np.mean(vn,axis=None)
      u_rms = np.sqrt(np.mean(np.square(ut_p)));
      v_rms = np.sqrt(np.mean(np.square(vn_p)));
      n_total = len(ut)
      print("--|| ALYA :: SIGNAL LENGTH=%d"% (n_total));
      BF1 = 100.0*float(len(np.where(ut<0.0)[0]))/n_total
      BF2 = 100.0*float(len(np.where(abs(vn/v_rms)>10.0)[0]))/n_total
      print("--|| ALYA :: BACKFLOW EVENTS=%.3f, EXTREME EVENTS=%.3f"% (BF1, BF2));
      Q1 = 100.0*float(len(np.where(np.logical_and(ut_p>0.0,vn_p>0.0))[0]))/n_total
      Q2 = 100.0*float(len(np.where(np.logical_and(ut_p<0.0,vn_p>0.0))[0]))/n_total
      Q3 = 100.0*float(len(np.where(np.logical_and(ut_p<0.0,vn_p<0.0))[0]))/n_total
      Q4 = 100.0*float(len(np.where(np.logical_and(ut_p>0.0,vn_p<0.0))[0]))/n_total
      print("--|| ALYA :: (Q1,Q2,Q3,Q4)=(%.0f,%.0f,%.0f,%.0f)"% (Q1, Q2, Q3, Q4));
      if('FLUC' in mode):
        df = pd.DataFrame(np.transpose((ut_p/u_rms,vn_p/v_rms)), columns=['$u_t/u_t^{rms}$','$v_n/v_n^{rms}$'])
        if(barScale==1):
          ticks = [0,1,2,3,4,5,6];
          clevels = np.logspace(-6, 0, num=12, endpoint=False);
        else:
          ticks = [0.05,0.15,0.25]
      else:
        df = pd.DataFrame(np.transpose((ut/u_rms,vn/v_rms)), columns=['$u_t/u_t^{rms}$','$v_n/v_n^{rms}$'])
        if(barScale==1):
          ticks = [0,1,2,3,4,6];
          clevels = np.logspace(-6, 0, num=12, endpoint=False);
        else:
          ticks = [0.1,0.2,0.3,0.5];
      plt_cbar = True
      if(j==0):
        plt_shade = True
      else:
        plt_shade = False
      if(j==0):
        if(barScale==1):
          with plot_kde_as_log():
            kdeplot = sns.jointplot(data=df, x="$u_t/u_t^{rms}$", y="$v_n/v_n^{rms}$", cmap=dcolor[j], \
                kind="kde", shade=plt_shade, cbar=plt_cbar, cbar_kws={"ticks":ticks}, label=lab,\
                marginal_kws={"color":lcolor[j], "shade":plt_shade, "bw_adjust":bwidth})
        else:
          kdeplot = sns.jointplot(data=df, x="$u_t/u_t^{rms}$", y="$v_n/v_n^{rms}$", cmap=dcolor[j], \
                kind="kde", shade=plt_shade, cbar=plt_cbar, cbar_kws={"ticks":ticks}, label=lab,\
                marginal_kws={"color":lcolor[j], "shade":plt_shade, "bw_adjust":bwidth})
        if(plt_cbar):
          plt.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.15)  # shrink fig so cbar is visible
          # get the current positions of the joint ax and the ax for the marginal x
          pos_joint_ax = kdeplot.ax_joint.get_position()
          pos_marg_x_ax = kdeplot.ax_marg_x.get_position()
          # reposition the joint ax so it has the same width as the marginal x ax
          kdeplot.ax_joint.set_position([pos_joint_ax.x0,pos_joint_ax.y0,pos_marg_x_ax.width,pos_joint_ax.height])
          # reposition the colorbar using new x positions and y positions of the joint ax
          kdeplot.fig.axes[-1].set_position([.87, 0.05, .05, pos_joint_ax.height/3])
          #kdeplot.fig.axes[-1].set_position([.87, 0.0*pos_joint_ax.y0, .05, pos_joint_ax.height/3])
      else:
        kdeplot.x = df['$u_t/u_t^{rms}$'];
        kdeplot.y = df['$v_n/v_n^{rms}$'];
        if(barScale==1):
          with plot_kde_as_log():
            kdeplot.plot_joint(sns.kdeplot, cmap=dcolor[j], shade=False, \
              cbar= plt_cbar, cbar_kws={"ticks":ticks}, label=lab)
        else:
          kdeplot.plot_joint(sns.kdeplot, cmap=dcolor[j], shade=False, \
              cbar= plt_cbar, cbar_kws={"ticks":ticks}, label=lab)
        if(plt_cbar):
          ## reposition the joint ax so it has the same width as the marginal x ax
          kdeplot.ax_joint.set_position([pos_joint_ax.x0,pos_joint_ax.y0,pos_marg_x_ax.width,pos_joint_ax.height])
          # reposition the colorbar using new x positions and y positions of the joint ax
          kdeplot.fig.axes[-1].set_position([.87, 0.75, .05, pos_joint_ax.height/3])
        if(barScale==1):
          with plot_kde_as_log():
            kdeplot.plot_marginals(sns.kdeplot, color=lcolor[j], bw_adjust=bwidth, shade=plt_shade)
        else:
          kdeplot.plot_marginals(sns.kdeplot, color=lcolor[j], bw_adjust=bwidth, shade=plt_shade)
      if(barScale==0):
        kdeplot.ax_marg_x.set_xlim(-3, 3)
        kdeplot.ax_marg_y.set_ylim(-3, 3)
      else:
        kdeplot.ax_marg_x.set_xlim(-6, 6)
        kdeplot.ax_marg_y.set_ylim(-6, 6)
      kdeplot.fig.suptitle('$x/c=%.1f$'%x_wit)
      ##g = sns.JointGrid(data=df, x="$u_t/u_t^{rms}$", y="$v_n/v_n^{rms}$", space=0)
      ##g.plot_joint(sns.kdeplot, cmap=dcolor[j], shade=plt_shade, cbar=plt_cbar, label=lab)
      ##sns.kdeplot(df["$u_t/u_t^{rms}$"],color=lcolor[j], shade=plt_shade, bw=0.1, ax=g.ax_marg_x)
      ##sns.kdeplot(df["$v_n/v_n^{rms}$"], color=lcolor[j], shade=plt_shade, bw=0.1, vertical=True, ax=g.ax_marg_y)
      
  elif("WITCHECK" in mode):
    fname = 'wit_single_point_stat_info.csv'
    x_loc = [0.2, 0.9, 0.95]
    if('SAVE' in mode):
      print('--|| ALYA :: READING WITNESS DATA.')
      stime = time.time()
      time_data, wit_data_0 = cfd.ExportReadProbe('./%s-VELOX.wit.bin'%(caseName))
      time_data, wit_data_1 = cfd.ExportReadProbe('./%s-VELOY.wit.bin'%(caseName))
      PDIST=np.zeros((nWit,1),dtype=float)
      BEVNT=np.zeros((nWit,1),dtype=float)
      EEVNT=np.zeros((nWit,1),dtype=float)
      QEVNT=np.zeros((nWit,4),dtype=float)
      QEDRS=np.zeros((nWit,1),dtype=float)
      QED24=np.zeros((nWit,1),dtype=float)
      INTMT=np.zeros((nWit,1),dtype=float)
      for i in range(nWit):
        x_wit,y_wit,z_wit = witPoints[i,:]
        if(np.logical_and(x_wit>0.0,x_wit<1.0)):
          indMin,PDIST[i,0] = closest_node((x_wit,y_wit), coordAir);
          theta = Fth(coordAir[indMin,0])
        else:
          theta=0.0
        u0 = wit_data_0[:,i]; u_avg = np.mean(u0,axis=None); 
        v0 = wit_data_1[:,i]; v_avg = np.mean(v0,axis=None);
        ut = (u0*np.cos(theta)+v0*np.sin(theta)); u_avg = np.mean(ut,axis=None);
        vn = (-u0*np.sin(theta)+v0*np.cos(theta)); v_avg = np.mean(vn,axis=None);
        ut_p = ut-np.mean(ut,axis=None)
        vn_p = vn-np.mean(vn,axis=None)
        u_rms = np.sqrt(np.mean(np.square(ut_p)));
        v_rms = np.sqrt(np.mean(np.square(vn_p)));
        n_total = len(ut)
        BEVNT[i,0] = 100.0*float(len(np.where(ut<0.0)[0]))/n_total
        EEVNT[i,0] = 100.0*float(len(np.where(abs(vn/v_rms)>10.0)[0]))/n_total
        QEVNT[i,0] = 100.0*float(len(np.where(np.logical_and(ut_p>0.0,vn_p>0.0))[0]))/n_total
        QEVNT[i,1] = 100.0*float(len(np.where(np.logical_and(ut_p<0.0,vn_p>0.0))[0]))/n_total
        QEVNT[i,2] = 100.0*float(len(np.where(np.logical_and(ut_p<0.0,vn_p<0.0))[0]))/n_total
        QEVNT[i,3] = 100.0*float(len(np.where(np.logical_and(ut_p>0.0,vn_p<0.0))[0]))/n_total
        QEDRS[i,0] = (QEVNT[i,1]+QEVNT[i,3])-(QEVNT[i,0]+QEVNT[i,2])
        QED24[i,0] = (QEVNT[i,1]-QEVNT[i,3])
        # Don't know the classical definition of intermittancy.
        INTMT[i,0] = 100.0*float(len(np.where(abs(ut_p)>(0.01*u_rms))[0]))/n_total
      header_str='XWIT, YWIT, ZWIT, PDIST, BEVNT, EEVNT, '
      header_str+='QEVNT:1, QEVNT:2, QEVNT:3, QEVNT:4, QEDRS, QED24, INTMT'
      np.savetxt(fname,np.c_[witPoints,PDIST,BEVNT,EEVNT,QEVNT,QEDRS,QED24,INTMT],\
        		fmt='%.2E', delimiter = ', ',header = header_str)
    if('PLOT' in mode):
      plotData = np.loadtxt(fname, comments='#',delimiter = ', ')
      #ax = plt.figure(figsize=(11.69, 8.27))
      for i in range(len(x_loc)):
        lab = r'$x/c=%.2f$'%(x_loc[i]);
        ind_loc = np.where(witPoints[:,0]==x_loc[i])[0];
        xPlot = witPoints[ind_loc,2];
        #---subplot I----#
        yPlot = plotData[ind_loc,6];
        ax = plt.subplot(211)
        plt.plot(xPlot, yPlot, linewidth=1.5, label=lab)
        plt.ylabel(r'$\%E_{Q_{2,4}-Q_{1,3}}$');
        #---subplot II----#
        yPlot = plotData[ind_loc,7];
        ax = plt.subplot(212)
        plt.plot(xPlot, yPlot, linewidth=1.5, label=lab)
      plt.xlabel(r'$z/c$');
      plt.ylabel(r'$\%E_{Q_{2}-Q_{4}}$');
      ax.legend(loc = (0.11, 2.45),prop={'size': 10},ncol=len(x_loc))


#plt.show()
plt.tight_layout()
savePath = casePath+'/plot_'+mode.lower()+'_'
if('PSD' in mode):
  listP = [str(pt) for pt in indP]
  listStr = "_".join(listP)
  savePath = savePath+listStr+'_'
elif('THIS' in mode):
  listP = [str(pt) for pt in indP]
  listStr = "_".join(listP)
  savePath = savePath+listStr+'_'
  savePath = savePath+calcVar.lower()+'_'
#elif('TPCORR' in mode):
#  savePath = savePath+'_'
elif('BACKFLOW' in mode):
  listP = [str(pt) for pt in indP]
  listStr = "_".join(listP)
  savePath = savePath+listStr+'_JPROB_'
if('WITCHECKSAVE' not in mode):  
  savePath = savePath+airfoil
  savePath = savePath+'.png'
  print('----|| INFO :: SAVING FIGURE AS ',savePath)
  #plt.savefig(savePath, format='png',\
  #            dpi=600, facecolor='w', edgecolor='w',\
  #	    orientation='portrait',transparent=True,\
  #	    bbox_inches=None,pad_inches=0.1,\
  #	    papertype='a4',frameon=None, metadata=None)
  plt.savefig(savePath, format='png', dpi=600)
