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
VARLIST        = ['AVVEL','RS_II','RS_IJ']
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
  #gradv = mesh.smooth(mesh.gradient(field['AVVEL']),iters=3)
  gradv = mesh.gradient(field['AVVEL'])
  rstvw = field['RS_II'][:,1] - field['RS_II'][:,2]  
  rstuv = -field['RS_IJ'][:,0]
  vorti = pyAlya.postproc.vorticity(gradv)


  field['AVVGR'] = gradv
  field['GAVOR'] = mesh.gradient(vorti[:,0]) #gradient of x-vorticity
  ggrsi = mesh.gradient(mesh.gradient(rstvw))
  field['GGRSI'] = ggrsi[:,5]
  ggrsi = mesh.gradient(mesh.gradient(rstuv))
  field['GGRSJ'] = ggrsi[:,4] - ggrsi[:,8]

  #Write the fields
  pyAlya.pprint(1,'Writing MPIO...',flush=True)
  field.write(CASESTR,instant,header.time,basedir=BASEDIR,fmt=FILE_FMT,exclude_vars=VARLIST)

pyAlya.cr_info()
