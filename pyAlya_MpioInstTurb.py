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
COMM = sys.argv[7]
PLOTVAR = sys.argv[8]

#Create list of instances
if(START==END):
  listOfInstants = [END]
else:
  listOfInstants = [ii for ii in range(START,END+DT,DT)]


## Create the subdomain mesh
mesh = pyAlya.Mesh.read(CASESTR,basedir=BASEDIR,alt_basedir=ALT_BASEDIR,fmt=FILE_FMT,read_commu=True if COMM == 1 else False,read_massm=False)


pyAlya.pprint(0,'Run (%d instants)...' % len(listOfInstants),flush=True)

time=0
for instant in listOfInstants:
  pyAlya.pprint(1,'First loop, instant %d...'%instant,flush=True)
  # Read field
  field, header = pyAlya.Field.read(CASESTR,VARLIST,instant,mesh.xyz,basedir=BASEDIR,fmt=FILE_FMT)
  
  # Compute and smooth the gradient of velocity
  gradv = mesh.smooth(mesh.gradient(field['VELOC']),iters=3)
  
  
  # Compute Vorticity, Q and Omega from the gradient
  if('ALL' not in PLOTVAR):
    if('GRADV' in PLOTVAR):
      # Store gradients
      field['GRADV'] = mesh.newArray(ndim=6)
      field['GRADV'][:,0] = gradv[:,0] # XX
      field['GRADV'][:,1] = gradv[:,4] # YY
      field['GRADV'][:,2] = gradv[:,8] # ZZ
      field['GRADV'][:,3] = gradv[:,1] # XY
      field['GRADV'][:,4] = gradv[:,2] # XZ
      field['GRADV'][:,5] = gradv[:,5] # YZ
    if('VOR' in PLOTVAR):
      field['VORTI'] = pyAlya.postproc.vorticity(gradv)
    if('QCR' in PLOTVAR):
      field['QCRIT'] = pyAlya.postproc.QCriterion(gradv)
    if('LAM' in PLOTVAR):
      field['LAMB2'] = pyAlya.postproc.Lambda2Criterion(gradv)
    if('OMG' in PLOTVAR):
      field['OMEGA'] = pyAlya.postproc.OmegaCriterion(gradv,epsilon=0.001,modified=False)
    if('ROR' in PLOTVAR):
      field['RORTX'] = pyAlya.postproc.RortexCriterion(gradv)
    if('ORX' in PLOTVAR):
      field['OMERX'] = pyAlya.postproc.OmegaRortexCriterion(gradv,epsilon=0.001,modified=False)
  else:
    field['VORTI'] = pyAlya.postproc.vorticity(gradv)
    field['QCRIT'] = pyAlya.postproc.QCriterion(gradv)
    field['LAMB2'] = pyAlya.postproc.Lambda2Criterion(gradv)
    field['OMEGA'] = pyAlya.postproc.OmegaCriterion(gradv,epsilon=0.001,modified=False)
    field['RORTX'] = pyAlya.postproc.RortexCriterion(gradv)
    field['OMERX'] = pyAlya.postproc.OmegaRortexCriterion(gradv,epsilon=0.001,modified=False)
  
  #Write the fields
  pyAlya.pprint(1,'Writing MPIO...',flush=True)
  field.write(CASESTR,instant,header.time,basedir=BASEDIR,fmt=FILE_FMT,exclude_vars=VARLIST)

pyAlya.cr_info()
