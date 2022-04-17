import os
import time
import sys

#os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import CFDlib as cfd
from os.path import expanduser
from scipy.interpolate import interp1d

#######################################################################
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

caseName	= sys.argv[1]
airfoil 	= sys.argv[2]
mode 		= sys.argv[3]

#indP = [6,7]
#indP = [288]
#indP = [4,5,6,7]
#indP = [0,1,2,3,4,5,6]
#indP = [11,12,13]
#indP = [22,23]
#indP = [5,6,7,8,9]
#indP = [2,3,4,5,6,7,8,9,10,11]
#indP = (np.arange(23,30)-1)
indP = [22,23,24]

# LOAD AIRFOIL SPECIFIC FILES
print('--|| ALYA INITIALIZING')
plt.close("all")
lw = 2;
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
z = np.unique(witPoints[:,2]);
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
  for n in range(nz):
   lab = r'$P_{'+str(n+1)+'}$';
   
   #if(n<11):
   # sgnX = -1
   # sgnY = 1 
   #else:
   sgnX = 1
   sgnY = -1 
   #if n in [0,1,2,11,12,13]:
   # sgnY = 0
   leg.append(lab)
   plt.plot(witPoints[n,0], witPoints[n,1],'o',linewidth=0.5)
   plt.text(witPoints[n,0]+(sgnX)*0.025, witPoints[n,1]+(sgnY)*0.025, lab, fontsize=MEDIUM_SIZE)
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
    calcVar = 'V'
    print('--|| INFO :: PSD CALCULATIONS ON %s'% calcVar)
    print('--|| ALYA :: READING WITNESS DATA.')
    stime = time.time()
    if('U' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./naca-VELOX.wit.bin')
    elif('V' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./naca-VELOY.wit.bin')
    elif('W' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./naca-VELOZ.wit.bin')
    elif('P' in calcVar):
      time_data, wit_data = cfd.ExportReadProbe('./naca-PRESS.wit.bin')
    else:
      raise ValueError('ALYA ERROR :: VARIABLE NOT SPECIFIED!')
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)

    print('----|| INFO :: %d TIME VALUES RANGE %.2f-%.2f'%(len(wit_data),time_data[0], time_data[-1]))
  

    import scipy.signal as signal
    print('--|| ALYA :: ARRANGING ARRAYS MONOTONICALLY.')
    stime = time.time()
    tM = np.maximum.accumulate(time_data)
    t, ind = np.unique(tM, return_index=True)
    print('----|| INFO :: UNIQUE TIME ARRAY SIZE=',np.size(t))
    markFreq = int(np.amax((np.floor(len(t)/20),1),axis=None))
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)
    for (i,j) in enumerate(indP):
      print('----|| INFO :: PROCESSING POINT %d, x,y,z='%j, witPoints[j,:])
      #lab = r'$P_{'+str(j+1)+'}$';
      lab = r'$x/c=%.1f, \, y/c=%.2f$'%(witPoints[j,0],witPoints[j,1]);
      if('ZAVG' in mode):
        print('----|| INFO :: PERFORMING SPANWISE AVERAGING')
        z_per_index = np.arange(j,nWit,nz)
        dataPoint = np.mean(wit_data[:,z_per_index],axis=1)
      else:
        dataPoint = wit_data[:,j]; #Extract data
      print('----|| INFO :: WITNESS ARRAY SIZE=',np.shape(dataPoint))
      dataPoint = dataPoint[ind];
      print('----|| INFO :: UNIQUE WITNESS ARRAY SIZE=',np.shape(dataPoint))
      tSignal = t;
      ySignal = dataPoint-np.mean(dataPoint,axis=None)
      #ax = plt.subplot(len(indP), 1, i+1)
      ax = plt.subplot(1, 1, 1)
      if("THIS" in mode):
        print('--|| ALYA :: CALCULATE TIME HISTORY (THIS).')
        N = 1000;
        y1 = dataPoint - np.convolve(dataPoint, np.ones(N)/N, mode='same')

        plt.plot(tSignal, dataPoint, 'k',markevery=markFreq,linewidth=0.5,label=lab)
        #plt.plot(t, y1, 'r--',markevery=markFreq,linewidth=0.5,label=lab)
        ax.annotate(lab, xy=(0.05, 1.2), xycoords='axes fraction', fontsize=MEDIUM_SIZE,
                        horizontalalignment='left', verticalalignment='top')
        if(i==len(indP)-1):
         plt.xlabel(r'$tU_o/c$');
         #ylab = r'$\left(u-\overline{U}\right)/U_o$'
         ylab = r'$u/U_o$'
         ax.annotate(ylab, xy=(0.02, 0.55),xycoords='figure fraction',fontsize=MEDIUM_SIZE,rotation=90)
      elif("PSD" in mode):
        print('--|| ALYA :: CALCULATE POWER SPECTRA (PSD).')
        if('LOMB' in mode):
          nout = 10000;
          freq_min = 0.01; freq_max = 1000.0;
          f = np.linspace(freq_min, freq_max, nout)
          pgram = signal.lombscargle(tSignal, ySignal, f)
        elif('FFT' in mode):
          print('--|| ALYA :: INTERPOLATE CONSTANT TIME.')
          fInterp = interp1d(tSignal,ySignal);
          delta_t = np.amax(np.diff(tSignal),axis=None)
          s_rate = 1.0/delta_t
          print('----|| INFO :: SAMPLING FREQ = %.2f'% s_rate)
          Nt = int((np.amax(tSignal,axis=None)-np.amin(tSignal,axis=None))/delta_t)+1
          print('----|| INFO :: INTERPOLATING ON %d POINTS'%(Nt))
          tSignal = np.linspace(np.amin(tSignal,axis=None),np.amax(tSignal,axis=None),Nt);
          ySignal = fInterp(tSignal)
          print('--|| ALYA :: DONE. TIME=',time.time()-stime)
          amp,pgram = calcSpectra(ySignal, 10, 0)
          #pgram = np.fft.fft(ySignal)
          N = len(pgram); print('--|| INFO :: FFT LENGTH = %d'%N) 
          n = np.arange(N); 
          T = float(N)/s_rate; 
          f = n/T
          #plt.stem(f, np.abs(pgram), 'b', \
          # markerfmt=" ", basefmt="-b")
          #plt.ylabel(r'$P_{u^\prime u^\prime}$')
          #ax.set(xlim=(0, 10))
        else:
          raise ValueError('ALYA ERROR :: METHOD TO CALC PSD NOT SPECIFIED!')
        pgram = 10**i*(pgram/np.amax(pgram,axis=None))
        plt.plot(f, pgram, markevery=markFreq,linewidth=0.5,label=lab)
        if('U' in calcVar):
          plt.ylabel(r'$E_{u^\prime u^\prime}$')
        elif('V' in calcVar):
          plt.ylabel(r'$E_{v^\prime v^\prime}$')
        elif('W' in calcVar):
          plt.ylabel(r'$E_{w^\prime w^\prime}$')
        elif('P' in calcVar):
          plt.ylabel(r'$E_{p^\prime p^\prime}$')
        plt.yscale("log")
        plt.xscale("log")
        ax.set(xlim=(np.amin(f), 1000))
        #plt.gca().set_xlim(right=1000);
        #ax.annotate(lab, xy=(0.05, 1.1), xycoords='axes fraction', fontsize=MEDIUM_SIZE,
        #                horizontalalignment='left', verticalalignment='top')
        ax.legend()
        if(i==0):
         plt.plot(f[np.where(f>100)],10.0*np.power(f[np.where(f>100)],-5.0/3),'k--')
        elif(i==len(indP)-1):
         plt.xlabel(r'$fc/U_o$');
      else:
        raise ValueError('ALYA ERROR :: CHECK MODE THAT YOU SPECIFIED!')

      #--------TWEAK THE SUBPLOTS-----------#
      if(i==0):
       # Hide the right and top spines
       ax.spines['top'].set_visible(False)
       ax.spines['bottom'].set_visible(True)
       ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True) # labels along the bottom edge are off
      elif(i==len(indP)-1):
       ax.spines['top'].set_visible(False)
       ax.spines['bottom'].set_visible(True)
       # Only show ticks on the left and bottom spines
       ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True) # labels along the bottom edge are off
      else:
       ax.spines['bottom'].set_visible(False)
       ax.spines['top'].set_visible(False)
       ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
       
      ax.spines['right'].set_visible(False)
      ax.tick_params(
        axis='y',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False)         # ticks along the top edge are off
      ax.tick_params(
        axis='y',          # changes apply to the x-axis
	which='major',      # both major and minor ticks are affected
	right=False)         # ticks along the top edge are off

      # plt.locator_params(axis='y', nbins=3)
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)

     
  elif("TPCORR" in mode):
    calcVar = 'UU'
    print('--|| INFO :: TWO-POINT CORRELATION ON %s'% calcVar)
    print('--|| ALYA :: READING WITNESS DATA.')
    stime = time.time()
    if('UU' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./naca-VELOX.wit.bin')
      wit_data_1 = wit_data_0;
    elif('VV' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./naca-VELOY.wit.bin')
      wit_data_1 = wit_data_0;
    elif('WW' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./naca-VELOZ.wit.bin')
      wit_data_1 = wit_data_0;
    elif('PP' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./naca-PRESS.wit.bin')
      wit_data_1 = wit_data_0;
    elif('UV' in calcVar):
      time_data, wit_data_0 = cfd.ExportReadProbe('./naca-VELOX.wit.bin')
      time_data, wit_data_1 = cfd.ExportReadProbe('./naca-VELOY.wit.bin')
    else:
      raise ValueError('ALYA ERROR :: VARIABLE NOT SPECIFIED!')
    print('----|| INFO :: %d TIME VALUES RANGE %.2f-%.2f'%(len(time_data),time_data[0], time_data[-1]))

    print('--|| ALYA :: DONE. TIME=',time.time()-stime)
    markFreq = int(np.amax((np.floor(len(z)/20),1),axis=None))
    print('--|| ALYA :: MARKER FREQ=',markFreq)
    print('--|| ALYA :: CALCULATING WIT POINT CORRELATIONS.')
    stime = time.time()
    uu = np.zeros((len(z),len(indP)),dtype=float)
    for j in indPLabel:
     u0 = wit_data_0[:,j+1]; u0 = u0-np.mean(u0);
     for i in range(len(z)):
      pair = i*nz+j+1;
      #print(j+1,pair);
      # Calculate the correlation b/w 2 signals
      u1 = wit_data_1[:,pair]; u1 = u1-np.mean(u1);
      prod = np.multiply(u0,u1)
      uu[i,j] = np.mean(prod)/(np.sqrt(np.mean(u0*u0))*np.sqrt(np.mean(u1*u1)))
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)

    print('--|| ALYA :: PLOTTING RESULTS.')
    
    fig = plt.figure()
    grid = plt.GridSpec(1, 1, wspace=0.3, hspace=0.8)
    
    axs = plt.subplot(grid[0, 0:]);
    for i in indPLabel:
      lab = r"$P_{"+str(indP[i]+1)+"}$";
      print("--|| ALYA :: PROCESSING POINT",lab);
      lname, lcode = linestyles[j][0:]
      l1, = plt.plot(z, uu[:,i], ls[i], linewidth=lw, color='black', 
            markevery=markFreq, markerfacecolor='black', markeredgecolor='black', label=lab)
    plt.xlabel(r'$z/c$');
    plt.ylabel(r'$\mathcal{R}_{u^\prime u^\prime}$')
    plt.xticks([0,0.05,0.1])
    #plt.yticks([0,0.5,1])
    axs.set(xlim=(0, 0.1))
    
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=True,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=True,
        labelright=False) # labels along the bottom edge are off
    axs.legend(loc="upper right",   # Position of legend	i
               bbox_to_anchor=(0.92,0.92),
               ncol = 2, fancybox=True
              )
#plt.show()
savePath = casePath+'/plot_'+mode.lower()+'_'
if('PSD' in mode):
  listP = [str(pt) for pt in indP]
  listStr = "_".join(listP)
  savePath = savePath+listStr+'_'
elif('TPCORR' in mode):
  savePath = savePath+calcVar.lower()+'_'
savePath = savePath+airfoil
savePath = savePath+'.png'
print('----|| INFO :: SAVING FIGURE AS ',savePath)
#plt.savefig(savePath, format='png',\
#            dpi=600, facecolor='w', edgecolor='w',\
#	    orientation='portrait',transparent=True,\
#	    bbox_inches=None,pad_inches=0.1,\
#	    papertype='a4',frameon=None, metadata=None)
plt.savefig(savePath, format='png', dpi=600)
