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

import numpy as np
import pyAlya
import os
import sys


# Parameters
rho, mu = 1.0, float(sys.argv[5])
lam     = 0.01

BASEDIR        = ''
ALT_BASEDIR    = ''
CASESTR        = sys.argv[1]
VARLIST        = ['VELOC','PRESS']
START, DT, END = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])

FILE_FMT = sys.argv[6]

# In case of restart, load the previous data
listOfInstants = [ii for ii in range(START,END+DT,DT)]


## Create the subdomain mesh
mesh = pyAlya.Mesh.read(CASESTR,basedir=BASEDIR,alt_basedir=ALT_BASEDIR,fmt=FILE_FMT,read_commu=True if FILE_FMT == 'mpio' else False,read_massm=False)


pyAlya.pprint(0,'Run (%d instants)...' % len(listOfInstants),flush=True)

field = pyAlya.Field(xyz  = pyAlya.truncate(mesh.xyz,6),
					OMEGA = mesh.newArray(ndim=3),  #Vorticity 
					QCRIT = mesh.newArray(),  	#Q-criterion 
					OCRIT = mesh.newArray(),        #Omega-criterion 
					LMDA2 = mesh.newArray(),        #Lambda 2 
					RORTX = mesh.newArray(),        #Rortex 
					OMRTX = mesh.newArray(),        #Omega-Rortex
					#SWIRL = mesh.newArray(),        #Swirl-Strength
				   )

time=0
for instant in listOfInstants:
  if instant%100 == 0: pyAlya.pprint(1,'First loop, instant %d...'%instant,flush=True)
  # Read field
  fields,header = pyAlya.Field.read(CASESTR,VARLIST,instant,mesh.xyz,basedir=BASEDIR,fmt=FILE_FMT)
  
  # Compute time-weighted average 
  dt   = header.time - time  # weight
  time = header.time         # sum_weights
  
  # Compute gradient of velocity
  #gradv = mesh.smooth(mesh.gradient(fields['VELOC']),iters=6)
  gradv = mesh.gradient(fields['VELOC'])
  
  field['OMEGA'] = pyAlya.postproc.vorticity(gradv)
  field['QCRIT'] = pyAlya.postproc.QCriterion(gradv)
  field['OCRIT'] = pyAlya.postproc.OmegaCriterion(gradv,epsilon=0.001,modified=False)
  field['LMDA2'] = pyAlya.postproc.Lambda2Criterion(gradv)
  field['RORTX'] = pyAlya.postproc.RortexCriterion(gradv)
  field['OMRTX'] = pyAlya.postproc.OmegaRortexCriterion(gradv,epsilon=0.001,modified=False)
  #field['SWIRL'] = pyAlya.postproc.SwirlStrength(gradv,normalized=False)
  
  
  #Write the fields
  pyAlya.pprint(1,'Writing MPIO...',flush=True)
  field.write(CASESTR,instant,time,basedir=ALT_BASEDIR,fmt=FILE_FMT)

pyAlya.cr_info()

