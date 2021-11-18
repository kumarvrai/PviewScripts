import os
import time
import sys

#os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#######################################################################
def calcSpectra(f, nWindows, iPeriodic):
  n = len(f)
  n = n/nWindows
  print(n)
  
  wss = 0.
  if iPeriodic:
    print('no window')
  else:
    wss = sum(np.hanning(n))
  
  for iWindow in range(0, nWindows, 1):
    print('data window: {}'.format(iWindow))
    print('start index: {}, end index: {}'.format(iWindow*n, (iWindow+1)*n-1))
    fx = f[iWindow*n:(iWindow+1)*n]
  
    if iPeriodic:
      fk = np.fft.fft(fx)
    else:
      wx = np.hanning(n)*fx
      fk = np.fft.fft(wx)
    
    if (n % 2) == 0:
      #Mind the odd-ball wavenumber, not considered here
      afk = 2*(np.abs(fk[0:(n/2+1)])/n)
      pfk = 2*(np.abs(fk[0:(n/2+1)])/n)**2
    else:
      afk = 2*(np.abs(fk[0:((n-1)/2+1)])/n)
      pfk = 2*(np.abs(fk[0:((n-1)/2+1)])/n)**2
  
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
  
  return pfks
#######################################################################

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
skipInd 	= int(sys.argv[4]);

skipL = 0;
#indP = [6,7]
#indP = [288]
#indP = [4,5,6,7]
indP = [0,1,2,3,4,5,6]
#indP = [11,12,13]
#indP = [22,23]
#indP = [5,6,7,8,9]
#indP = [2,3,4,5,6,7,8,9,10,11]
#indP = [24,25,26,27,28,29]

# LOAD AIRFOIL SPECIFIC FILES
niaScrh = '/home/kvishal/projects/rrg-ugo-ac/kvishal/research/0.Alya/'
niaHome = '/home/u/ugo/kvishal/'


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
  fileLoc = niaScrh+'/5.BGQ-DATA/NACA-0012/2.3D/2.LES/1.Re-5E4-alpha-8/3.WALE/1.run1/1.t-0-13/'
  fAirU = '/home/kvishal/1.post_process/1.airfoil/3.PviewScripts/1.0012-Slices/naca0012-UP.txt'
  fAirL = '/home/kvishal/1.post_process/1.airfoil/3.PviewScripts/1.0012-Slices/naca0012-DOWN.txt'
elif('4412' in airfoil):  
  fileLoc = niaScrh+'/2.NACA4412/2.ROUGH_CASE/1.trip/1.WRLES/2.Medium/'
  fAirU = '/home/kvishal/1.post_process/1.airfoil/3.PviewScripts/2.4412-Slices/naca4412-UP.txt'
  fAirL = '/home/kvishal/1.post_process/1.airfoil/3.PviewScripts/2.4412-Slices/naca4412-DOWN.txt'
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
nz = nWit/len(z)
indPLabel = range(len(indP));
print('--|| ALYA :: PROCESSING', nWit, 'TOTAL WITNESS POINTS')
print('--|| ALYA :: NUMBER OF WIT POINTS IN A PLANE ',nz,'WITH',len(z),'PLANES')

n = len(indP)
color = plt.cm.rainbow(np.linspace(0.0,1.0,n)) # This returns RGBA; convert:
hexcolor = map(lambda rgb:'#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255),
               tuple(color[:,0:-1]))
mpl.rcParams['axes.color_cycle'] = hexcolor

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
  fname = caseName+'.nsi.wit'
  nLinesWit	= len(open('%s'%fname).readlines())
  if(skipL):
   skipLines = nLinesWit-30000*(nWit+2)
   print('--|| ALYA :: SKIPPING ',skipLines,'LINES OUT OF', nLinesWit,'LINES IN FILE')
  else:
   skipLines = 0
  #data = np.loadtxt(fname, dtype=float, comments='#', skiprows=skipLines)
  print('--|| ALYA :: SKIPPING EVERY %d DATA POINTS.' % skipInd)
  data = np.empty((0,5),dtype='float')
  
  print('--|| ALYA :: READING WIT POINT TIME.')
  stime = time.time()
  t = [];
  for line in open(fname):
    li=line.strip()
    if 'Time' in li:
      t.append([float(li.split()[3])])
  t = np.asarray(t,dtype='float')   
  if(skipInd):
    ind = np.arange(0, t.size, skipInd)
    t = t[ind]
  print('--|| ALYA :: DONE. TIME=',time.time()-stime)
  print('--|| ALYA :: DATA-SET LENGTH = %d' % t.size)
  
  if(any(string in mode for string in ["PSD","THIS"])):
    import scipy.signal as signal
    indSpec = 1; print('----|| INFO :: USING VARIABLE INDEX ',indSpec)
    print('--|| ALYA :: READING WIT POINT DATA.')
    stime = time.time()
    if(skipInd):
     ii = 1;  #Count iterations
     jj = 1 #Count witness points
     for line in open(fname):
       li=line.strip()
       if not li.startswith("#"):
         if (ii-1) % skipInd == 0: #Skip every skipInd iteration number
           #if any(jj==(pp+1) for pp in indP):
           li = np.array(map(float, line.split()));
           data = np.vstack((data,li))
         jj+=1;
         if (jj-1) % nWit == 0:
           ii+=1;
           jj=1
    else: 
     data = np.loadtxt(fname, comments='#')
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)
    print('--|| ALYA :: ARRANGING ARRAYS MONOTONICALLY.')
    stime = time.time()
    tM = np.maximum.accumulate(t)
    print('----|| INFO :: TIME ARRAY SIZE=',np.size(t))
    t, ind = np.unique(tM, return_index=True)
    print('----|| INFO :: UNIQUE TIME ARRAY SIZE=',np.size(t))
    t.reshape(t.size,1)
    markFreq = int(np.amax((np.floor(len(t)/20),1),axis=None))
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)
    for (i,j) in enumerate(indP):
      print('----|| INFO :: PROCESSING -->',j)
      #lab = r"$x/c="+str(witPoints[i,0])+",\; y/c="+str(witPoints[i,1])+"$";
      lab = r'$P_{'+str(j+1)+'}$';
      dataPoint = data[np.where(data[:,0]==(j+1)),:]; #Extract data
      dataPoint = np.reshape(dataPoint, (dataPoint.shape[1],dataPoint.shape[2]))
      print('----|| INFO :: WITNESS ARRAY SIZE=',np.shape(dataPoint))
      dataPoint = dataPoint[ind,:];
      print('----|| INFO :: UNIQUE WITNESS ARRAY SIZE=',np.shape(dataPoint))
      ySignal = dataPoint[:,indSpec]-np.mean(dataPoint[:,indSpec],axis=None)
      ax = plt.subplot(len(indP), 1, i+1)
      if("THIS" in mode):
        print('--|| ALYA :: CALCULATE TIME HISTORY (THIS).')
	N = 1000;
	y1 = dataPoint[:,indSpec] - np.convolve(dataPoint[:,indSpec], np.ones(N)/N, mode='same')

        plt.plot(t, dataPoint[:,indSpec], 'k',markevery=markFreq,linewidth=0.5,label=lab)
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
        nout = 10000
        f = np.linspace(0.01, 2000, nout)
        pgram = signal.lombscargle(t, ySignal, f)
        plt.plot(f, pgram, 'k',markevery=markFreq,linewidth=lw,label=lab)
        ax.annotate(lab, xy=(0.05, 1.1), xycoords='axes fraction', fontsize=MEDIUM_SIZE,
                        horizontalalignment='left', verticalalignment='top')
        plt.ylabel(r'$E_{u''u''}$')
        plt.yscale("log")
        plt.xscale("log")
        ax.set(xlim=(1, 2000))
        #ax.legend()
        if(i==0):
         plt.plot(f[np.where(f>100)],10000.0*np.power(f[np.where(f>100)],-5.0/3),'b--')
        elif(i==len(indP)-1):
         plt.xlabel(r'$fc/U_o$');
      else:
        raise ValueError('ALYA ERROR :: CHECK MODE THAT YOU SPECIFIED!')

      #--------TWEAK THE SUBPLOTS-----------#
      if(i==0):
       # Hide the right and top spines
       ax.spines['top'].set_visible(False)
       ax.spines['bottom'].set_visible(False)
       ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
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

      plt.locator_params(axis='y', nbins=3)

    plt.show()
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)

     
  elif("TPCORR" in mode):
    markFreq = int(np.amax((np.floor(len(z)/20),1),axis=None))
    print('--|| ALYA :: MARKER FREQ=',markFreq)
    print('--|| ALYA :: READING WIT POINT DATA.')
    stime = time.time()
    if(skipInd):
     cntGlob = 0
     count = 1;
     for line in open(fname):
       li=line.strip()
       if not li.startswith("#"):
         cntGlob+=1;
         if (cntGlob-1) % nWit == 0:
          count+=1;
         if count % skipInd == 0:
          li = np.array(map(float, line.split()));
          data = np.vstack((data,li))
    else:	  
     data = np.loadtxt(fname, comments='#')
    print('--|| ALYA :: DONE. TIME=',time.time()-stime)
    print('--|| ALYA :: CALCULATING WIT POINT CORRELATIONS.')
    stime = time.time()
    uu = np.zeros((len(z),len(indP)),dtype=float)
    uv = np.zeros((len(z),len(indP)),dtype=float)
    pp = np.zeros((len(z),len(indP)),dtype=float)
    up = np.zeros((len(z),len(indP)),dtype=float)
    for j in indPLabel:
     u0 = data[np.where(data[:,0]==(j+1)),1]; u0 = u0-np.mean(u0,axis=1);
     v0 = data[np.where(data[:,0]==(j+1)),2]; v0 = v0-np.mean(v0,axis=1);
     ip = 3
     p0 = data[np.where(data[:,0]==(j+1)),ip]; p0 = p0-np.mean(p0,axis=1);
     for i in range(len(z)):
      pair = i*nz+j+1;
      #print(j+1,pair);
      # Calculate uu
      u1 = data[np.where(data[:,0]==pair),1]; u1 = u1-np.mean(u1,axis=1);
      prod = np.multiply(u0,u1)
      uu[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(u0*u0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
      # Calculate uv
      u1 = data[np.where(data[:,0]==pair),2]; u1 = u1-np.mean(u1,axis=1);
      prod = np.multiply(u0,u1)
      uv[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(u0*u0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
      # Calculate pp
      u1 = data[np.where(data[:,0]==pair),3]; u1 = u1-np.mean(u1,axis=1);
      prod = np.multiply(p0,u1)
      pp[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(p0*p0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
      # Calculate up
      u1 = data[np.where(data[:,0]==pair),3]; u1 = u1-np.mean(u1,axis=1);
      prod = np.multiply(u0,u1)
      up[i,j] = np.mean(prod,axis=1)/(np.sqrt(np.mean(u0*u0,axis=1))*np.sqrt(np.mean(u1*u1,axis=1)))
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
    
    #axs = plt.subplot(grid[1, 0:]);
    #for i in indPLabel:
    #  lab = r"$P_{"+str(i)+"}$";
    #  l2, = plt.plot(z,uv[:,i],'k'+ls[i],linewidth=lw,label=lab)
    #plt.xlabel(r'$z/c$');
    #plt.ylabel(r'$\mathcal{R}_{u \prime v\prime}$')
    #plt.xticks([0,0.05,0.1])
    #axs.set(xlim=(0, 0.1))
    #  #plt.yticks([-1,0,1])
    #
    #axs = plt.subplot(grid[2, 0:]);
    #for i in indPLabel:
    #  lab = r"$P_{"+str(i)+"}$";
    #  l3, = plt.plot(z,pp[:,i],'k'+ls[i],linewidth=lw,label=lab)
    #plt.xlabel(r'$z/c$');
    #plt.ylabel(r'$\mathcal{R}_{p\prime p\prime}$')
    #axs.set(xlim=(0, 0.1))
    #plt.xticks([0,0.05,0.1])
    #  #plt.yticks([0,0.5,1])
    #
    #axs = plt.subplot(grid[3, 0:]);
    #for i in indPLabel:
    #  lab = r"$P_{"+str(i)+"}$";
    #  l4, = plt.plot(z,up[:,i],'k'+ls[i],linewidth=lw,label=lab)
    #plt.xlabel(r'$z/c$');
    #plt.ylabel(r'$\mathcal{R}_{u\prime p \prime}$')
    #axs.set(xlim=(0, 0.1))
    #plt.xticks([0,0.05,0.1])
      #plt.yticks([0,0.5,1])
    handles, labels = axs.get_legend_handles_labels()  
    #frame1 = plt.gca()
    
    #axs = plt.subplot(grid[1, 0:]);
    #plt.plot(coordAirU[:,0],coordAirU[:,1],'r',linewidth=2)
    #plt.plot(coordAirL[:,0],coordAirL[:,1],'r',linewidth=2)
    #leg = []
    #for i,n in enumerate(indP):
    # lab = r'$P_{'+str(i)+'}$';
    # print("--|| ALYA :: PROCESSING POINT",lab);
    # leg.append(lab)
    # plt.plot(witPoints[n,0], witPoints[n,1],'ko',linewidth=1)
    # plt.text(witPoints[n,0], witPoints[n,1]+0.025, lab, fontsize=SMALL_SIZE)
    #axs.set_anchor('W')
    #plt.xlabel(r'$x/c$');
    #plt.yticks([])
    #plt.xticks([0,0.5,1,1.3])
    #axs.set_aspect(1.5);
    #axs.set(xlim=(-0.001, 1.3))
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    ##axs.spines['bottom'].set_visible(False)
    #axs.spines['left'].set_visible(False)
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=True,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=True,
        labelright=False) # labels along the bottom edge are off
    fig.legend(handles,     # The line objects
               labels,
               loc="upper right",   # Position of legend	i
               bbox_to_anchor=(0.92,0.92),
               ncol = 2, fancybox=True
              )
#plt.show()
savePath = casePath+'/plot_'+mode.lower()+'_'
savePath = savePath+airfoil
savePath = savePath+'.png'
print('----|| INFO :: SAVING FIGURE AS ',savePath)
#plt.savefig(savePath, format='png',\
#            dpi=600, facecolor='w', edgecolor='w',\
#	    orientation='portrait',transparent=True,\
#	    bbox_inches=None,pad_inches=0.1,\
#	    papertype='a4',frameon=None, metadata=None)
plt.savefig(savePath, format='png', dpi=600)
