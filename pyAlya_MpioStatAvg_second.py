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
import sys
import os


# Parameters
rho, mu = 1.0, float(sys.argv[5])
lam     = 0.01

BASEDIR        = ''
ALT_BASEDIR    = ''
CASESTR        = sys.argv[1]
VARLIST        = ['VELOC','PRESS','AVPRE', 'AVVEL', 'AVVE2', 'AVVXY','AVTAN']
#START, DT, END = 2,1,415
START, DT, END = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])

FILE_FMT = 'mpio'
SAVE_MPIO      = True
COMPUTE_EARSM  = False

# Partial runs
RUN_FIRST_LOOP  = True
RUN_SECOND_LOOP = False

# In case of restart, load the previous data
listOfInstants = [ii for ii in range(START,END+DT,DT)]


## Create the subdomain mesh
#mesh = pyAlya.Mesh.read(CASESTR,basedir=BASEDIR,alt_basedir=ALT_BASEDIR,fmt=FILE_FMT,read_commu=True if FILE_FMT == 'mpio' else False,read_massm=False)
mesh = pyAlya.Mesh.read(CASESTR,basedir=BASEDIR,alt_basedir=ALT_BASEDIR,fmt=FILE_FMT,read_commu=False,read_massm=False)

pyAlya.pprint(0,'Run (%d instants)...' % len(listOfInstants),flush=True)

## Accumulate the statistics (auxiliar according to Table 5)
stats = pyAlya.Field(xyz  = pyAlya.truncate(mesh.xyz,6),
					# Here are the mandatory defined in Table 2 of HiFiTurb (GA 814837)
					# Level 1 - averaged Navier-Stokes equations
					AVPRE = mesh.newArray(),        # Averaged pressure
					AVVEL = mesh.newArray(ndim=3),  # Averaged velocity
					AVTAN = mesh.newArray(ndim=3),  # Averaged velocity
					AVVE2 = mesh.newArray(ndim=3),  # Averaged velocity
					AVVXY = mesh.newArray(ndim=3),  # Averaged velocity
					AVTEM = mesh.newArray(),        # Averaged temperature
					GRAVP = mesh.newArray(ndim=3),  # Averaged gradient of pressure
					GRAVV = mesh.newArray(ndim=9),  # Averaged gradient of velocity
					AVHFL = mesh.newArray(ndim=3),	# Averaged heat flux
					AVSTR = mesh.newArray(ndim=9),  # Averaged strain rate
					AVROT = mesh.newArray(ndim=9),  # Averaged rotation rate
					AVSHE = mesh.newArray(ndim=9),  # Averaged shear stresses
					RESTR = mesh.newArray(ndim=9),  # Reynolds stresses
					AVSTF = mesh.newArray(ndim=9),  # Averaged strain rate
					AVRTF = mesh.newArray(ndim=9),  # Averaged rotation rate
					AVTHF = mesh.newArray(ndim=3),	# Averaged turbulent heat flux
					# Level 1 - additional quantities
					AVPF2 = mesh.newArray(),		# Pressure autocorrelation
					AVTF2 = mesh.newArray(),		# Temperature autocorrelation
					TAYMS = mesh.newArray(),        # Taylor microscale
					KOLLS = mesh.newArray(),        # Kolmogorov lenghtscale
					KOLTS = mesh.newArray(),        # Kolmogorov timescale
					# Level 2 - Reynolds stress equations budget terms
					CONVE = mesh.newArray(ndim=9),  # Convection
					PRODU = mesh.newArray(ndim=9),  # Production
					DIFF1 = mesh.newArray(ndim=9),  # Turbulent diffusion 1
					DIFF2 = mesh.newArray(ndim=9),  # Turbulent diffusion 2
					DIFF3 = mesh.newArray(ndim=9),  # Molecular diffusion
					PSTRA = mesh.newArray(ndim=9),  # Pressure strain
					DISSI = mesh.newArray(ndim=9),  # Dissipation
					# Level 2 - Reynolds stress equations - separate terms
					AVPVE = mesh.newArray(ndim=3),  # Pressure velocity correlation
					AVVE3 = mesh.newArray(ndim=27)  # Triple velocity correlation
				   )

# Only accumulate pressure, velocity and temperature (if available) to
# obtain the fluctuations on the next loop
time = 0
for instant in listOfInstants:
	if instant%100 == 0: pyAlya.pprint(1,'First loop, instant %d...'%instant,flush=True)
	# Read field
	fields,header = pyAlya.Field.read(CASESTR,VARLIST,instant,mesh.xyz,basedir=BASEDIR,fmt=FILE_FMT)

	# Compute time-weighted average 
	dt   = header.time - time  # weight
	time = header.time         # sum_weights

	# Accumulate the statistics (Welford's online algorithm)
	stats['AVPRE'] += pyAlya.stats.addS1(stats['AVPRE'],fields['AVPRE'],w=1. if instant == START else dt/time)
	stats['AVVEL'] += pyAlya.stats.addS1(stats['AVVEL'],fields['AVVEL'],w=1. if instant == START else dt/time)
	stats['AVVE2'] += pyAlya.stats.addS1(stats['AVVE2'],fields['AVVE2'],w=1. if instant == START else dt/time)
	stats['AVVXY'] += pyAlya.stats.addS1(stats['AVVXY'],fields['AVVXY'],w=1. if instant == START else dt/time)
	stats['AVTAN'] += pyAlya.stats.addS1(stats['AVTAN'],fields['AVTAN'],w=1. if instant == START else dt/time)

# Gradients of averaged velocity and pressure
# only computed once
stats['GRAVP'] = mesh.gradient(stats['AVPRE'])
stats['GRAVV'] = mesh.gradient(stats['AVVEL'])
# calculate RS stresses from Alya files
stats['RESTR'] = 0.0*stats['GRAVV']
stats['RESTR'][:,0] = stats['AVVE2'][:,0]-stats['AVVEL'][:,0]**2 #uu
stats['RESTR'][:,4] = stats['AVVE2'][:,1]-stats['AVVEL'][:,1]**2 #vv
stats['RESTR'][:,8] = stats['AVVE2'][:,2]-stats['AVVEL'][:,2]**2 #ww
stats['RESTR'][:,1] = stats['AVVXY'][:,0]-stats['AVVEL'][:,0]*stats['AVVEL'][:,1]  #uv
stats['RESTR'][:,2] = stats['AVVXY'][:,2]-stats['AVVEL'][:,1]*stats['AVVEL'][:,2]  #vw
stats['RESTR'][:,5] = stats['AVVXY'][:,1]-stats['AVVEL'][:,0]*stats['AVVEL'][:,2]  #uw
stats['RESTR'][:,3] = stats['RESTR'][:,1]
stats['RESTR'][:,6] = stats['RESTR'][:,2]
stats['RESTR'][:,7] = stats['RESTR'][:,5]


stats.write(CASESTR,0,0.,basedir=ALT_BASEDIR,fmt='mpio',exclude_vars=[
	'AVTEM','AVHFL','AVSTR','AVROT','AVSHE','AVSTF','AVRTF','AVTHF',
	'AVPF2','AVTF2','TAYMS','KOLLS','KOLTS','CONVE','PRODU','DIFF1','DIFF2',
	'DIFF3','PSTRA','DISSI','AVPVE','AVVE3'
])

exit(0) # Stop the run so as not to overload the GPFS

## Do a second loop in time
# This time compute all the necessary magnitudes and accumulate 
# them as needed 
time = 0
for instant in listOfInstants:
	if instant%100 == 0: pyAlya.pprint(1,'Second loop, instant %d...'%instant,flush=True)
	
	# Read field
	fields,tStep = pyAlya.Field.read(CASESTR,VARLIST,instant,mesh.xyz,basedir=BASEDIR,fmt=FILE_FMT)

	# Compute time-weighted average 
	dt   = header.time - time  # weight
	time = header.time         # sum_weights

	# Postprocess fields
	fields['GRADV'] = mesh.gradient(fields['VELOC'])
	fields['STRAI'] = pyAlya.stats.strainTensor(fields['GRADV'])
	fields['ROTAT'] = pyAlya.stats.vorticityTensor(fields['GRADV'])
	fields['SHEAR'] = 2.*mu*fields['STRAI']

	# Fluctuations
	fields['PFLUC'] = pyAlya.math.linopScaf(1.,fields['PRESS'],-1.,stats['AVPRE'])  # p' = p - <p>
	fields['VFLUC'] = pyAlya.math.linopArrf(1.,fields['VELOC'],-1.,stats['AVVEL'])  # u' = u - <u>
	fields['PFLU2'] = fields['PFLUC']*fields['PFLUC']
	fields['PVCOR'] = pyAlya.math.scaVecProd(fields['PFLUC'],fields['VFLUC'])
	fields['VELO3'] = rho*pyAlya.stats.tripleCorrelation(fields['VFLUC'],fields['VFLUC'],fields['VFLUC'])

	fields['GRVFL'] = pyAlya.math.linopArrf(1.,fields['GRADV'],-1.,stats['GRAVV'])
	fields['STRAF'] = pyAlya.stats.strainTensor(fields['GRVFL'])
	fields['ROTAF'] = pyAlya.stats.vorticityTensor(fields['GRVFL'])

	# Budgets
	fields['PSTRA'] = pyAlya.stats.pressureStrainBudget(fields['PFLUC'],fields['STRAF'])
	fields['DISSI'] = pyAlya.stats.dissipationBudget(mu,fields['GRVFL'])

	# Accumulate statistics
	#stats['AVPVE'] += pyAlya.stats.addS1(stats['AVPVE'],fields['PVCOR'],w=1. if instant == START else dt/time)
	#stats['AVPF2'] += pyAlya.stats.addS1(stats['AVPF2'],fields['PFLU2'],w=1. if instant == START else dt/time)
	#stats['AVVE3'] += pyAlya.stats.addS1(stats['AVVE3'],fields['VELO3'],w=1. if instant == START else dt/time)
	#stats['RESTR'] += pyAlya.stats.addS1(stats['RESTR'],fields['RESTR'],w=1. if instant == START else dt/time)
	#stats['AVSTR'] += pyAlya.stats.addS1(stats['AVSTR'],fields['STRAI'],w=1. if instant == START else dt/time)
	#stats['AVROT'] += pyAlya.stats.addS1(stats['AVROT'],fields['ROTAT'],w=1. if instant == START else dt/time)
	#stats['AVSHE'] += pyAlya.stats.addS1(stats['AVSHE'],fields['SHEAR'],w=1. if instant == START else dt/time)
	#stats['AVSTF'] += pyAlya.stats.addS1(stats['AVSTF'],fields['STRAF'],w=1. if instant == START else dt/time)
	#stats['AVRTF'] += pyAlya.stats.addS1(stats['AVRTF'],fields['ROTAF'],w=1. if instant == START else dt/time)

	#stats['PSTRA'] += pyAlya.stats.addS1(stats['PSTRA'],fields['PSTRA'],w=1. if instant == START else dt/time)
	#stats['DISSI'] += pyAlya.stats.addS1(stats['DISSI'],fields['DISSI'],w=1. if instant == START else dt/time)

	stats['AVPVE'] += fields['PVCOR']
	stats['AVPF2'] += fields['PFLU2']
	stats['AVVE3'] += fields['VELO3']
	stats['AVSTR'] += fields['STRAI']
	stats['AVROT'] += fields['ROTAT']
	stats['AVSHE'] += fields['SHEAR']
	stats['AVSTF'] += fields['STRAF']
	stats['AVRTF'] += fields['ROTAF']
	stats['PSTRA'] += fields['PSTRA']
	stats['DISSI'] += fields['DISSI']

stats['AVPVE'] = stats['AVPVE']/len(listOfInstants) 
stats['AVPF2'] = stats['AVPF2']/len(listOfInstants) 
stats['AVVE3'] = stats['AVVE3']/len(listOfInstants) 
stats['AVSTR'] = stats['AVSTR']/len(listOfInstants) 
stats['AVROT'] = stats['AVROT']/len(listOfInstants) 
stats['AVSHE'] = stats['AVSHE']/len(listOfInstants) 
stats['AVSTF'] = stats['AVSTF']/len(listOfInstants) 
stats['AVRTF'] = stats['AVRTF']/len(listOfInstants) 
stats['PSTRA'] = stats['PSTRA']/len(listOfInstants) 
stats['DISSI'] = stats['DISSI']/len(listOfInstants)

# Compute TKE and dissipation
k     = pyAlya.stats.TKE(stats['RESTR'])
dissi = 0.5*pyAlya.math.trace(stats['DISSI']) # e = 1/2*e_ii = 2*mu*<S'_ij S'_ij>

# Compute Taylor microscale and Kolmogorov length and time scales
stats['TAYMS'] = pyAlya.stats.taylorMicroscale(mu/rho,k,dissi)
stats['KOLLS'] = pyAlya.stats.kolmogorovLengthScale(mu/rho,dissi)
stats['KOLTS'] = pyAlya.stats.kolmogorovTimeScale(mu/rho,dissi)

# Finish budgets
stats['CONVE'] = pyAlya.stats.convectionBudget(stats['AVVEL'],mesh.gradient(stats['RESTR']))
stats['PRODU'] = pyAlya.stats.productionBudget(stats['RESTR'],stats['GRAVV'])
stats['DIFF1'] = pyAlya.stats.turbulentDiffusion1Budget(rho,stats['AVVE3'],mesh)
stats['DIFF2'] = pyAlya.stats.turbulentDiffusion2Budget(stats['AVPVE'],mesh)
stats['DIFF3'] = pyAlya.stats.molecularDiffusionBudget(mu,stats['RESTR'],mesh)

prod  = 0.5*pyAlya.math.trace(stats['PRODU'])

## Write MPIO if requested
if SAVE_MPIO:
	pyAlya.pprint(1,'Writing MPIO...',flush=True)
	stats.write(CASESTR,0,0.,basedir=ALT_BASEDIR,fmt='mpio',exclude_vars=[
		'AVTEM','AVHFL','AVSTR','AVROT','AVSHE','AVSTF','AVRTF','AVTHF',
		'AVPF2','AVTF2','AVPVE','AVVE3'])

exit(0) # Stop the run so as not to overload the GPFS

## Compute EARSM stats
if COMPUTE_EARSM:
	II_S = pyAlya.math.doubleDot(stats['AVSTR'],pyAlya.math.transpose(stats['AVSTR'])) # SijSji
	II_O = pyAlya.math.doubleDot(stats['AVROT'],pyAlya.math.transpose(stats['AVROT'])) # OijOji
	stats['s']    = np.sqrt(pyAlya.math.linopScaf(1.,II_S,-1.,II_O))
	stats['SIGM'] = stats['s']*k/dissi 
	stats['SIGM'][np.isnan(stats['SIGM'])] = 0. # Assign any possible NaN to 0 due to dividing by 0
	stats['r']    = -II_O/stats['s']/stats['s']
	stats['r'][np.isnan(stats['r'])] = 0. # Assign any possible NaN to 0 due to dividing by 0
	Sad  = pyAlya.math.scaTensProd(1./stats['s'],stats['AVSTR']) # Adimensionalization of Sij
	Oad  = pyAlya.math.scaTensProd(1./stats['s'],stats['AVROT']) # Adimensionalization of Oij
	stats['aij']  = pyAlya.math.scaTensProd(1./k,stats['RESTR']) - 2./3.*pyAlya.math.identity(stats['RESTR'])
	stats['aij'][np.isnan(stats['aij'])] = 0. # Assign any possible NaN to 0 due to dividing by 0
	# Invariants
	stats['III_S'] = pyAlya.math.tripleDot(Sad,Sad,Sad)
	stats['IV']    = pyAlya.math.tripleDot(Sad,Oad,Oad)
	stats['V']     = pyAlya.math.quatrupleDot(Sad,Sad,Oad,Oad) + 0.5*stats['r']*(1.-stats['r'])
	stats['II_a']  = pyAlya.math.doubleDot(stats['aij'],pyAlya.math.transpose(stats['aij'])) # aijaji
	stats['III_a'] = pyAlya.math.tripleDot(stats['aij'],stats['aij'],stats['aij'])
	# Values from the new EARSM
	OijOij = pyAlya.math.matmul(Oad,Oad)
	OijSij = pyAlya.math.matmul(Oad,Sad)
	SijOij = pyAlya.math.matmul(Sad,Oad)
	stats['Tij1']  = Sad
	stats['Tij2']  = SijOij - OijSij
	stats['Tij3']  = OijOij + 1./3.*pyAlya.math.scaTensProd(stats['r'],pyAlya.math.identity(Oad))
	stats['Tij4']  = pyAlya.math.matmul(Sad,OijOij) + pyAlya.math.matmul(OijOij,Sad) \
	               - 2./3.*pyAlya.math.scaTensProd(stats['IV'],pyAlya.math.identity(Oad)) \
	               + pyAlya.math.scaTensProd(stats['r'],Sad)
	stats['Tij5']  = pyAlya.math.matmul(OijSij,OijOij) - pyAlya.math.matmul(OijOij,SijOij) \
	               - 0.5*pyAlya.math.scaTensProd(stats['r'],pyAlya.math.linopArrf(1.,SijOij,-1.,OijSij))
	# Computation of A_kl = Tij^l T_ji^k
	Tij = [stats['Tij1'],stats['Tij2'],stats['Tij3'],stats['Tij4'],stats['Tij5']]
	A   = np.zeros((len(stats),25),dtype=np.double) # a 5x5 matrix at every grid point
	for k in range(5):
		for l in range(5):
			A[:,5*k+l] = pyAlya.math.doubleDot(Tij[l],pyAlya.math.transpose(Tij[k]))
	# Regularization Bij = Aij + lam*delta_ij
	B    = A + lam*pyAlya.math.identity(A)
	Binv = pyAlya.math.inverse(B)
	# Compute Iat vector Iat_k = a_ijT_ji^k
	Iat = np.zeros((len(stats),5),dtype=np.double)
	for ii in range(5):
		Iat[:,ii] = pyAlya.math.doubleDot(stats['aij'],pyAlya.math.transpose(Tij[ii]))
	# Compute b_k = (B^-1)_kl a_ij T_ji^l = Binv_kl Iat_k
	beta = pyAlya.math.tensVecProd(Binv,Iat)
	stats['beta1'] = beta[:,0]
	stats['beta2'] = beta[:,1]
	stats['beta3'] = beta[:,2]
	stats['beta4'] = beta[:,3]
	stats['beta5'] = beta[:,4]
	# Compute PsK
	stats['PsK'] = prod/stats['s']/k
	stats['PsK'][np.isnan(stats['PsK'])] = 0. # Assign any possible NaN to 0 due to dividing by 0
	# Compute EARSM extra fields
	stats['LAVEL'] = mesh.divergence(fields['GRAVV']) # lap(U) = div(grad(U))
	stats['SADOT'] = mesh.newArray(ndim=9)            # Sadim dot
	for ii in range(9):
		gSijs = mesh.gradient(fields['AVSTR'][:,ii])
		stats['SADOT'][:,ii] = pyAlya.math.dot(fields['AVVEL'],gSijs) # U·grad(Sij/s)



## Write h5 database
pyAlya.pprint(1,'Writing database...',flush=True)
# We will make use of the parallel capabilities of
# the h5 io so that every node will write its own 
# part of the domain and the master will write the
# masterfile. First we need the total number of nodes
nnodG, nelG = mesh.nnodG,mesh.nelG
pyAlya.pprint(0,'nnodG=',nnodG,'nelG=',nelG,flush=True) # Master prints nnod

# All nodes write the database files
# parallel HDF5 requires all nodes to write the file
# otherwise a deadlock is produced
writer = pyAlya.io.HiFiTurbDB_Writer(nnodG)

# Nodes
writer.writeDataset('Nodes',
			 mesh.filter_bc(stats.x),
			 mesh.filter_bc(stats.y),
			 mesh.filter_bc(stats.z)	
)

# Create lists with the statistics
writer.writeDataset('Inputs',
			 # Basic averages
			 mesh.filter_bc(stats['AVPRE'][:]),   # p
			 mesh.filter_bc(stats['AVVEL'][:,0]), # u
			 mesh.filter_bc(stats['AVVEL'][:,1]), # v
			 mesh.filter_bc(stats['AVVEL'][:,2]), # w
			 # Shear stress
			 mesh.filter_bc(stats['AVSHE'][:,0]), # tau_11
			 mesh.filter_bc(stats['AVSHE'][:,1]), # tau_12
			 mesh.filter_bc(stats['AVSHE'][:,4]), # tau_22
			 mesh.filter_bc(stats['AVSHE'][:,2]), # tau_13
			 mesh.filter_bc(stats['AVSHE'][:,5]), # tau_23
			 mesh.filter_bc(stats['AVSHE'][:,8]), # tau_33
			 # Reynolds stress tensor
			 mesh.filter_bc(stats['RESTR'][:,0]), # r_11
			 mesh.filter_bc(stats['RESTR'][:,1]), # r_12
			 mesh.filter_bc(stats['RESTR'][:,4]), # r_22
			 mesh.filter_bc(stats['RESTR'][:,2]), # r_13
			 mesh.filter_bc(stats['RESTR'][:,5]), # r_23
			 mesh.filter_bc(stats['RESTR'][:,8]), # r_33
		     # Gradients
			 mesh.filter_bc(stats['GRAVP'][:,0]), # px
			 mesh.filter_bc(stats['GRAVV'][:,0]), # ux
			 mesh.filter_bc(stats['GRAVV'][:,3]), # vx
			 mesh.filter_bc(stats['GRAVV'][:,6]), # wx
			 mesh.filter_bc(stats['GRAVP'][:,1]), # py
			 mesh.filter_bc(stats['GRAVV'][:,1]), # uy
			 mesh.filter_bc(stats['GRAVV'][:,4]), # vy
			 mesh.filter_bc(stats['GRAVV'][:,5]), # wy
			 mesh.filter_bc(stats['GRAVP'][:,2]), # pz
			 mesh.filter_bc(stats['GRAVV'][:,2]), # uz
			 mesh.filter_bc(stats['GRAVV'][:,5]), # vz
			 mesh.filter_bc(stats['GRAVV'][:,8])  # wz
)

# Additional quantities
writer.writeDataset('AdditionalQuantities',
			 # Pressure autocorrelation
			 mesh.filter_bc(stats['AVPF2'][:]),   # pp
			 # Taylor microscale
			 mesh.filter_bc(stats['TAYMS'][:]),   # Tm
			 # Kolmogorov length scale
			 mesh.filter_bc(stats['KOLLS'][:]),   # Kl
			 # Kolmogorov time scale
			 mesh.filter_bc(stats['KOLTS'][:]),   # Kt	
)

# Dissipation
writer.writeDataset('Dissipation',
			 mesh.filter_bc(stats['DISSI'][:,0]), # d_11
			 mesh.filter_bc(stats['DISSI'][:,1]), # d_12
			 mesh.filter_bc(stats['DISSI'][:,4]), # d_22
			 mesh.filter_bc(stats['DISSI'][:,2]), # d_13
			 mesh.filter_bc(stats['DISSI'][:,5]), # d_23
			 mesh.filter_bc(stats['DISSI'][:,8])  # d_33
)

# Convection
writer.writeDataset('Convection',
			 mesh.filter_bc(stats['CONVE'][:,0]), # c_11
			 mesh.filter_bc(stats['CONVE'][:,1]), # c_12
			 mesh.filter_bc(stats['CONVE'][:,4]), # c_22
			 mesh.filter_bc(stats['CONVE'][:,2]), # c_13
			 mesh.filter_bc(stats['CONVE'][:,5]), # c_23
			 mesh.filter_bc(stats['CONVE'][:,8])  # c_33
)

# Production
writer.writeDataset('Production',
			 mesh.filter_bc(stats['PRODU'][:,0]), # p_11
			 mesh.filter_bc(stats['PRODU'][:,1]), # p_12
			 mesh.filter_bc(stats['PRODU'][:,4]), # p_22
			 mesh.filter_bc(stats['PRODU'][:,2]), # p_13
			 mesh.filter_bc(stats['PRODU'][:,5]), # p_23
			 mesh.filter_bc(stats['PRODU'][:,8])  # p_33
)

# Molecular Diffusion
writer.writeDataset('MolecularDiffusion',
			 mesh.filter_bc(stats['DIFF3'][:,0]), # md_11
			 mesh.filter_bc(stats['DIFF3'][:,1]), # md_12
			 mesh.filter_bc(stats['DIFF3'][:,4]), # md_22
			 mesh.filter_bc(stats['DIFF3'][:,2]), # md_13
			 mesh.filter_bc(stats['DIFF3'][:,5]), # md_23
			 mesh.filter_bc(stats['DIFF3'][:,8])  # md_33
)

# Pressure Strain
writer.writeDataset('PressureStrain',
			 mesh.filter_bc(stats['PSTRA'][:,0]), # ps_11
			 mesh.filter_bc(stats['PSTRA'][:,1]), # ps_12
			 mesh.filter_bc(stats['PSTRA'][:,4]), # ps_22
			 mesh.filter_bc(stats['PSTRA'][:,2]), # ps_13
			 mesh.filter_bc(stats['PSTRA'][:,5]), # ps_23
			 mesh.filter_bc(stats['PSTRA'][:,8])  # ps_33
)

# Pressure Velocity
writer.writeDataset('PressureVelocity',
			 mesh.filter_bc(stats['AVPVE'][:,0]), # pv_1
			 mesh.filter_bc(stats['AVPVE'][:,1]), # pv_2
			 mesh.filter_bc(stats['AVPVE'][:,2])  # pv_3
)

# Turbulent Diffusion 01
writer.writeDataset('TurbulentDiffusion01',
			 mesh.filter_bc(stats['DIFF1'][:,0]), # td01_11
			 mesh.filter_bc(stats['DIFF1'][:,1]), # td01_12
			 mesh.filter_bc(stats['DIFF1'][:,4]), # td01_22
			 mesh.filter_bc(stats['DIFF1'][:,2]), # td01_13
			 mesh.filter_bc(stats['DIFF1'][:,5]), # td01_23
			 mesh.filter_bc(stats['DIFF1'][:,8])  # td01_33
)

# Turbulent Diffusion 02
writer.writeDataset('TurbulentDiffusion02',
			 mesh.filter_bc(stats['DIFF2'][:,0]), # td02_11
			 mesh.filter_bc(stats['DIFF2'][:,1]), # td02_12
			 mesh.filter_bc(stats['DIFF2'][:,4]), # td02_22
			 mesh.filter_bc(stats['DIFF2'][:,2]), # td02_13
			 mesh.filter_bc(stats['DIFF2'][:,5]), # td02_23
			 mesh.filter_bc(stats['DIFF2'][:,8])  # td02_33
)

# Triple Correlation
writer.writeDataset('TripleCorrelation',
			 mesh.filter_bc(stats['AVVE3'][:,0]),  # t_111
			 mesh.filter_bc(stats['AVVE3'][:,3]),  # t_121
			 mesh.filter_bc(stats['AVVE3'][:,12]), # t_221
			 mesh.filter_bc(stats['AVVE3'][:,6]),  # t_131
			 mesh.filter_bc(stats['AVVE3'][:,15]), # t_231
			 mesh.filter_bc(stats['AVVE3'][:,24]), # t_331
			 mesh.filter_bc(stats['AVVE3'][:,13]), # t_222
			 mesh.filter_bc(stats['AVVE3'][:,16]), # t_232
			 mesh.filter_bc(stats['AVVE3'][:,25]), # t_332
			 mesh.filter_bc(stats['AVVE3'][:,26])  # t_333
)

# EARSM statistics
if COMPUTE_EARSM:
	writer.writeDataset('EARSM',
						mesh.filter_bc(stats['PsK']),
						mesh.filter_bc(stats['aij'][:,0]), # a11
						mesh.filter_bc(stats['aij'][:,1]), # a12
						mesh.filter_bc(stats['aij'][:,4]), # a22
						mesh.filter_bc(stats['aij'][:,2]), # a13
						mesh.filter_bc(stats['aij'][:,5]), # a23
						mesh.filter_bc(stats['aij'][:,8]), # a33
						mesh.filter_bc(stats['SIGM']),
						mesh.filter_bc(stats['r']),
						mesh.filter_bc(stats['III_S']),
						mesh.filter_bc(stats['IV']),
						mesh.filter_bc(stats['V']),
						mesh.filter_bc(stats['II_a']),
						mesh.filter_bc(stats['III_a']),
						mesh.filter_bc(stats['Tij1'][:,0]),
						mesh.filter_bc(stats['Tij1'][:,1]),
						mesh.filter_bc(stats['Tij1'][:,4]),
						mesh.filter_bc(stats['Tij1'][:,2]),
						mesh.filter_bc(stats['Tij1'][:,5]),
						mesh.filter_bc(stats['Tij1'][:,8]),
						mesh.filter_bc(stats['Tij2'][:,0]),
						mesh.filter_bc(stats['Tij2'][:,1]),
						mesh.filter_bc(stats['Tij2'][:,4]),
						mesh.filter_bc(stats['Tij2'][:,2]),
						mesh.filter_bc(stats['Tij2'][:,5]),
						mesh.filter_bc(stats['Tij2'][:,8]),
						mesh.filter_bc(stats['Tij3'][:,0]),
						mesh.filter_bc(stats['Tij3'][:,1]),
						mesh.filter_bc(stats['Tij3'][:,4]),
						mesh.filter_bc(stats['Tij3'][:,2]),
						mesh.filter_bc(stats['Tij3'][:,5]),
						mesh.filter_bc(stats['Tij3'][:,8]),
						mesh.filter_bc(stats['Tij4'][:,0]),
						mesh.filter_bc(stats['Tij4'][:,1]),
						mesh.filter_bc(stats['Tij4'][:,4]),
						mesh.filter_bc(stats['Tij4'][:,2]),
						mesh.filter_bc(stats['Tij4'][:,5]),
						mesh.filter_bc(stats['Tij4'][:,8]),
						mesh.filter_bc(stats['Tij5'][:,0]),
						mesh.filter_bc(stats['Tij5'][:,1]),
						mesh.filter_bc(stats['Tij5'][:,4]),
						mesh.filter_bc(stats['Tij5'][:,2]),
						mesh.filter_bc(stats['Tij5'][:,5]),
						mesh.filter_bc(stats['Tij5'][:,8]),
						mesh.filter_bc(stats['beta1']),
						mesh.filter_bc(stats['beta2']),
						mesh.filter_bc(stats['beta3']),
						mesh.filter_bc(stats['beta4']),
						mesh.filter_bc(stats['beta5']),
	)

	writer.writeDataset('EARSM_Extra',
				# Laplacian of velocity
				mesh.filter_bc(stats['LAVEL'][:,0]),  # lap(U)
				mesh.filter_bc(stats['LAVEL'][:,1]),  # lap(V)
				mesh.filter_bc(stats['LAVEL'][:,2]),  # lap(W)
				# Convective derivative of Sij/s
				mesh.filter_bc(stats['SADOT'][:,0]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,1]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,2]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,3]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,4]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,5]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,6]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,7]),  # U·grad(Sij/s)
				mesh.filter_bc(stats['SADOT'][:,8]),  # U·grad(Sij/s)
	)


if mesh.comm.rank == 0:
	print('Writing master file',flush=True)
	# Create the dictionary for the master file
	writer.createGroup('01_Info',{
			'n_nodes'   : writer.createDataset('n_nodes',  (1,),'i',nnodG,            ret=True),
			'n_elems'   : writer.createDataset('n_elems',  (1,),'i',nelG,             ret=True),
			'density'   : writer.createDataset('density',  (1,),'f',rho,              ret=True),
			'viscosity' : writer.createDataset('viscosity',(1,),'f',mu,               ret=True),
			'model'     : writer.createDataset('model',    (1,),'s','incompressible', ret=True),
		},
	ret=False)

	writer.createGroup('02_Entries',{
			'Inputs' : writer.createExternalLink('Inputs','./Inputs.h5',ret=True),
			'output' : writer.createGroup('01_Output',{
					'AdditionalQuantities' : writer.createExternalLink('AdditionalQuantities','./AdditionalQuantities.h5',ret=True),
					'Convection'           : writer.createExternalLink('Convection','./Convection.h5',ret=True),
					'Production'           : writer.createExternalLink('Production','./Production.h5',ret=True),
					'TurbulentDiffusion01' : writer.createExternalLink('TurbulentDiffusion01','./TurbulentDiffusion01.h5',ret=True),
					'TurbulentDiffusion02' : writer.createExternalLink('TurbulentDiffusion02','./TurbulentDiffusion02.h5',ret=True),
					'MolecularDiffusion'   : writer.createExternalLink('MolecularDiffusion','./MolecularDiffusion.h5',ret=True),
					'PressureStrain'       : writer.createExternalLink('PressureStrain','./PressureStrain.h5',ret=True),
					'Dissipation'          : writer.createExternalLink('Dissipation','./Dissipation.h5',ret=True),
					'TripleCorrelation'    : writer.createExternalLink('TripleCorrelation','./TripleCorrelation.h5',ret=True),
					'PressureVelocity'     : writer.createExternalLink('PressureVelocity','./PressureVelocity.h5',ret=True),
				},
			ret=True),
			'earsm' : writer.createGroup('02_EARSM',{
					'EARSM'       : writer.createExternalLink('EARSM'      ,'./EARSM.h5'      ,ret=True),
					'EARSM_Extra' : writer.createExternalLink('EARSM_Extra','./EARSM_Extra.h5',ret=True),
				},
			ret=True)
		},
	ret=False)

	writer.createGroup('03_Nodes',{
			'Nodes' : writer.createExternalLink('Nodes','./Nodes.h5',ret=True),
		},
	ret=False)

	# Master writes the master file
	writer.writeMaster('Statistics.h5')

pyAlya.cr_info()
