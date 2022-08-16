#!/usr/bin/python
activate_this = '/home/u/ugo/kvishal/softwares/coolVenv/bin/activate_this.py'
import os
import glob
import sys
import time
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

pyVer = sys.version_info[0]
if pyVer < 3:
  execfile(activate_this, dict(__file__=activate_this))
else:
  exec(open(activate_this).read(), {'__file__': activate_this})
import modred as mr
#####################[IOData]###########################
CODE=sys.argv[5]
MODE=sys.argv[6]
print('--|| ALYA : WORKING WITH ',CODE)
# working directory and data filename
niaHome = '/home/u/ugo/kvishal/'
niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'
# input directory
ID = os.getcwd()+'/'
# input .foam file
if('ALYA' in CODE):
 IF = ID+sys.argv[1]+'.ensi.case'
elif('NEK' in CODE):
 IF = ID+sys.argv[1]+'.nek5000'
elif('PVD' in CODE):
 IF = ID+sys.argv[1]+'.pvd'
elif('VTM' in CODE):
 fNames = glob.glob("*.vtm")
 IF = [ID+ff for ff in fNames]
else:
 raise Exception('--|| ALYA : PROVIDE A VALID CODE FILE ')
print('--|| ALYA : READING FILE ',IF)
# output directory
OD = os.getcwd()+'/'
#start time
t0=0.0
# end time
tf=200.0
# name of field to be decomposed
fieldname = sys.argv[2]
DIM = int(sys.argv[3])

print('--|| ALYA : READING VARIABLE ',fieldname)
if('TKENORM' in fieldname.upper()):
 varName_code = ['PRESS','VELOC']
 varName_calc = []
 wtsMat = np.array([0, 1, 1, 1])
 arrayList = varName_code
elif('UNORM' in fieldname.upper()):
 varName_code = ['PRESS','VELOC']
 varName_calc = []
 wtsMat = np.array([0, 1, 0, 0])
 arrayList = varName_code
elif('VNORM' in fieldname.upper()):
 varName_code = ['PRESS','VELOC']
 varName_calc = []
 wtsMat = np.array([0, 0, 1, 0])
 arrayList = varName_code
elif('PNORM' in fieldname.upper()):
 varName_code = ['PRESS','VELOC']
 varName_calc = []
 wtsMat = np.array([1, 0, 0, 0])
 arrayList = varName_code
elif('VORTINORM' in fieldname.upper()):
 varName_code = ['PRESS','VELOC']
 varName_calc = ['VORTI']
 wtsMat = np.array([0, 0, 1])
 arrayList = []
elif('BERNORM' in fieldname.upper()):
 varName_code = ['PRESS','VELOC']
 varName_calc = []
 wtsMat = np.array([1, 0.5, 0.5, 0.5])
 arrayList = varName_code
elif('COMPND' in fieldname.upper()):
 varName_code = ['PRESS','VELOC']
 varName_calc = ['RS_II','RS_IJ']
 wtsMat = np.array([0, 1, 1, 1, 0, 0, 0, 0, 0, 0])
 arrayList = varName_code
elif('UVPOS' in fieldname.upper()):
 varName_code = ['VELOC']
 varName_calc = ['UVPOS']
 wtsMat = np.array([0, 0, 0, 1])
 arrayList = varName_code
else:
 raise Exception('--|| ALYA : PROVIDED VARAIABLE IDENTIFIER MISSING ')

print('--|| ALYA : DATASET IS',DIM,'DIMENSION')
# vector or not
field_is_vector = int(sys.argv[4])
# prefix of spatial modes' name
prefix="Re40_"
# True= read data from case files; 
# False= read data from fields.npz
read_fields_from_file = True
calculate_pod_slice   = False
box_clip              = True
# Variable interpolation form= %(NAME)s is available
fields_filename="""{}fields_{}.npz""".format(OD,fieldname)
geom_filename="""{}Geometry.vtu""".format(OD)
times_filename="""{}times.npy""".format(OD)

do_POD = False; do_DMD=False;
if('POD' in MODE):
  do_POD = True
elif('DMD' in MODE):  
  do_DMD = True
else:
 raise Exception('--|| ALYA : PROVIDED MODE IDENTIFIER MISSING ')
#####################[POD]###########################
# accuracy determines how many modes will be calculated.
accuracy=0.9999
# spatial mode to be calculated, M=0 means determined by program according to accuracy
M_POD = 20
#
subtractAvg = True
#
useVolWeight = True
#
output_correlation_matrix=True
POD_cm_filename="""{}POD_correlation_matrix.csv""".format(OD)

output_POD_temporal_modes=True
POD_tm_filename="""{}POD_temporal_coefficients.csv""".format(OD)

output_POD_spatial_modes=True
POD_sm_filename="""{}POD_spatial_modes.vtm""".format(OD)

doReconstruction = False
doFullRecon = False
# reconstruct output numbers
MR_POD = 2
ReconTime= (t0+tf)/2.0
POD_reconstruction_filename="""{}POD_reconstruction.vtu""".format(OD)

#####################[DMD]###########################
# setting of DMD decomposition
M_DMD=2

output_DMD_info=True
DMD_info_filename="""{}DMD_info.csv""".format(OD)

output_DMD_build_coeffs=True
DMD_build_coeffs_filename="""{}DMD_build_coeffs.csv""".format(OD)

output_DMD_spatial_modes=True
DMD_sm_filename="""{}DMD_spatial_modes.vtm""".format(OD)
