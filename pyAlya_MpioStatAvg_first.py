#!/bin/env python
#
# HiFiTurb database computation.
#
# Last rev: 19/02/2021
from __future__ import print_function, division

# Please do not delete this part otherwise it will not work
# you have been warned after a long weekend of debugging
import mpi4py
mpi4py.rc.recv_mprobe = False

import sys
import numpy as np
import pyAlya


# Parameters
rho, mu = 1.0, float(sys.argv[5])
lam     = 0.01

BASEDIR        = ''
ALT_BASEDIR    = ''
CASESTR        = str(sys.argv[1])
VARLIST        = ['AVPRE', 'AVVEL', 'AVVE2', 'AVVXY', 'AVTAN']
#START, DT, END = 2,1,415
START, DT, END = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])

FILE_FMT = str(sys.argv[6])
COMM = int(sys.argv[7])

# Partial runs
RUN_FIRST_LOOP  = True
RUN_SECOND_LOOP = False

# In case of restart, load the previous data
listOfInstants = [ii for ii in range(START,END+DT,DT)]


mesh = pyAlya.Mesh.read(CASESTR,basedir=BASEDIR,alt_basedir=ALT_BASEDIR,fmt=FILE_FMT,read_commu=True if COMM == 1 else False,read_massm=False)

pyAlya.pprint(0,'Run (%d instants)...' % len(listOfInstants),flush=True)

## Accumulate the statistics (auxiliar according to Table 5)
stats = pyAlya.Field(xyz  = pyAlya.truncate(mesh.xyz,6),
					# Here are the mandatory defined in Table 2 of HiFiTurb (GA 814837)
					# Level 1 - averaged Navier-Stokes equations
					AVPRE = mesh.newArray(),        # Averaged pressure
					AVVEL = mesh.newArray(ndim=3),  # Averaged velocity
					AVVE2 = mesh.newArray(ndim=3),        # Averaged V_ii
					AVVXY = mesh.newArray(ndim=3),  # Averaged V_ij
					AVTAN = mesh.newArray(ndim=3),  # Averaged V_ij
					RS_II = mesh.newArray(ndim=3),  # Averaged normal RS
					RS_IJ = mesh.newArray(ndim=3),  # Averaged shear RS
					AVPGR = mesh.newArray(ndim=3),  # Averaged gradient of pressure
					AVVGR = mesh.newArray(ndim=9),  # Averaged gradient of velocity
				   )

time = 0
for instant in listOfInstants:
	if instant%100 == 0: pyAlya.pprint(1,'First loop, instant %d...'%instant,flush=True)
	# Read field
	fields,header = pyAlya.Field.read(CASESTR,VARLIST,instant,mesh.xyz,basedir=BASEDIR,fmt=FILE_FMT)

	# Compute time-weighted average 
	dt   = header.time - time  # weight
	time = header.time         # sum_weights

	# Accumulate the statistics
	stats['AVPRE'] += pyAlya.stats.addS1(stats['AVPRE'],fields['AVPRE'],w=1. if instant == START else dt/time)
	stats['AVVEL'] += pyAlya.stats.addS1(stats['AVVEL'],fields['AVVEL'],w=1. if instant == START else dt/time)
	stats['AVVE2'] += pyAlya.stats.addS1(stats['AVVE2'],fields['AVVE2'],w=1. if instant == START else dt/time)
	stats['AVVXY'] += pyAlya.stats.addS1(stats['AVVXY'],fields['AVVXY'],w=1. if instant == START else dt/time)
	stats['AVTAN'] += pyAlya.stats.addS1(stats['AVTAN'],fields['AVTAN'],w=1. if instant == START else dt/time)

# Gradients of averaged velocity and pressure
stats['AVPGR'] = mesh.gradient(stats['AVPRE'])
stats['AVVGR'] = mesh.gradient(stats['AVVEL'])
# calculate RS stresses from Alya files
stats['RS_II'][:,0] = stats['AVVE2'][:,0]-stats['AVVEL'][:,0]**2   #uu
stats['RS_II'][:,1] = stats['AVVE2'][:,1]-stats['AVVEL'][:,1]**2   #vv
stats['RS_II'][:,2] = stats['AVVE2'][:,2]-stats['AVVEL'][:,2]**2   #ww
stats['RS_IJ'][:,0] = stats['AVVXY'][:,0]-stats['AVVEL'][:,0]*stats['AVVEL'][:,1]  #uv
stats['RS_IJ'][:,1] = stats['AVVXY'][:,1]-stats['AVVEL'][:,1]*stats['AVVEL'][:,2]  #vw
stats['RS_IJ'][:,2] = stats['AVVXY'][:,2]-stats['AVVEL'][:,0]*stats['AVVEL'][:,2]  #uw


stats.write(CASESTR,0,0.,basedir=ALT_BASEDIR,fmt=FILE_FMT)
pyAlya.cr_info()
