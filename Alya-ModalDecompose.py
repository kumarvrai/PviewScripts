initializeVars = '/home/u/ugo/kvishal/1.post_process/0.alya_pv_scripts/initializeVariables.py'
import sys
from paraview.simple import *
pyVer = sys.version_info[0]
if pyVer < 3:
  execfile(initializeVars)
else:
  #exec(open(initializeVars).read(), {'__file__': initializeVars})
  exec(open(initializeVars).read())
# LOAD SNAPSHOTS
print('--|| ALYA :: INITIALIZING PYTHON VERSION:',pyVer)
niaHome = '/home/u/ugo/kvishal/'
niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'

if any(idstr in MODE for idstr in ["POD","DMD"]):
  print("--|| ALYA :: READING ALYA ARRAYS", varName_code)
  startTime = time.time()
  case = OpenDataFile(IF)
  if('ALYA' in CODE):
    try:
     case.PointArrays = varName_code
    except:
     print('--|| ALYA: LOADING DEFAULT VARIABLES')
    case.UpdatePipeline()
  if('NEK' in CODE):
    case.PointArrays = ['velocity','pressure']
    case.UpdatePipeline()
    ## create a new 'Programmable Filter and change names'
    print("--|| NEK: CHANGING VARNAMES 1 USING A PROGRAMMABLE FILTER")
    startTime = time.time()
    case = ProgrammableFilter(Input=case)
    case.Script = \
    """
    import numpy as np
    varNames0 = ['velocity','pressure']
    varNames1 = ['VELOC','PRESS']
    for (i,var) in enumerate(varNames0):
     outName = varNames1[i]
     avg = (inputs[0].PointData[var])
     output.PointData.append(avg,outName)
    """
    case.UpdatePipeline()
    print("--|| NEK :: DONE. TIME =",time.time()-startTime,'sec')
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  reader = GetActiveSource()
  view = GetActiveView()
  times = reader.TimestepValues
  print("--|| ALYA :: TOTAL TEMPORAL STEPS =",len(times))
  
  if calculate_pod_slice:
    print("--|| ALYA :: APPLYING CLEAN TO GRID FILTER.")
    startTime = time.time()
    case = CleantoGrid(Input=case)
    #slice.SliceType.Origin = [0.0, 0.0, 0.0]
    case.UpdatePipeline()
    print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
    if(DIM==2):
      case = Calculator(Input=case)
      case.ResultArrayName = "VELOC"
      case.Function = " VELOC_X*iHat+VELOC_Y*jHat+0*VELOC_X*kHat"
      case.UpdatePipeline()
  if box_clip:
    print("--|| ALYA :: APPLYING BOX-CLIP FILTER.")
    startTime = time.time()
    case = Clip(Input=case)
    case.ClipType = 'Box'
    case.ClipType.Position = [0.0, 0.0, 0]
    case.ClipType.Length = [0.5, 0.3, 1.0]
    #slice.SliceType.Origin = [0.0, 0.0, 0.0]
    case.UpdatePipeline()
    print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
    #print("--|| ALYA :: TAKING DATA SLICE")
    #startTime = time.time()
    #case = Slice(Input=case)
    #case.SliceType = 'Plane'
    ##slice.SliceType.Origin = [0.0, 0.0, 0.0]
    #case.SliceType.Normal = [0.0, 0.0, 1.0]
    #case.UpdatePipeline()
    #print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  # TEMPORAL AVERAGING
  #print("--|| ALYA :: TEMPORAL AVERAGING OF CODE ARRAYS")
  #startTime = time.time()
  #case2 = TemporalStatistics(Input=case)
  #case2.ComputeMinimum = 0
  #case2.ComputeMaximum = 0
  #case2.ComputeStandardDeviation = 0
  #case2.UpdatePipeline()
  #print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  ## APPEND DATASETS
  #print("--|| ALYA :: APPEND DATASETS")
  #startTime = time.time()
  #case = AppendAttributes(Input=[case,case2])
  #case.UpdatePipeline()
  #print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  # CALCULATE EXTRA VARIABLES
  print("--|| ALYA :: CALCULATING EXTRA VARIABLES")
  startTime = time.time()
  
  if("VORTI" in varName_calc):
    case = GradientOfUnstructuredDataSet(Input=case)
    case.ScalarArray = ['POINTS', 'VELOC']
    case.ComputeGradient=0
    case.ComputeQCriterion=0
    case.ComputeVorticity=1
    case.VorticityArrayName='VORTI'
    case.UpdatePipeline()
    case = Calculator(Input=case)
    case.ResultArrayName = "VORTI"
    case.Function = " VORTI_Z"
    case.UpdatePipeline()
    arrayList.append('VORTI')
  if("RS_II" in varName_calc):
    case = Calculator(Input=case)
    case.ResultArrayName = "RS_II"
    case.Function = " (VELOC_X*VELOC_X - VELOC_average_X*VELOC_X -\
    VELOC_average_X*VELOC_X + VELOC_average_X*VELOC_average_X)*iHat + \
    (VELOC_Y*VELOC_Y - VELOC_average_Y*VELOC_Y -\
    VELOC_average_Y*VELOC_Y + VELOC_average_Y*VELOC_average_Y)*jHat + \
    (VELOC_Z*VELOC_Z - VELOC_average_Z*VELOC_Z -\
    VELOC_average_Z*VELOC_Z + VELOC_average_Z*VELOC_average_Z)*kHat "
    case.UpdatePipeline()
    arrayList.append('RS_II')
  if("RS_IJ" in varName_calc):
    case = Calculator(Input=case)
    case.ResultArrayName = "RS_IJ"
    case.Function = " (VELOC_X*VELOC_Y - VELOC_average_X*VELOC_Y -\
    VELOC_average_Y*VELOC_X + VELOC_average_X*VELOC_average_Y)*iHat + \
    (VELOC_Z*VELOC_Y - VELOC_average_Z*VELOC_Y -\
    VELOC_average_Y*VELOC_Z + VELOC_average_Z*VELOC_average_Y)*jHat + \
    (VELOC_Z*VELOC_X - VELOC_average_Z*VELOC_X -\
    VELOC_average_X*VELOC_Z + VELOC_average_Z*VELOC_average_X)*kHat "
    case.UpdatePipeline()
    arrayList.append('RS_IJ')
  if("UVPOS" in varName_calc):
    case = Calculator(Input=case)
    case.ResultArrayName = "UVPOS"
    case.Function = " (VELOC_X*VELOC_Y - VELOC_average_X*VELOC_Y -\
    VELOC_average_Y*VELOC_X + VELOC_average_X*VELOC_average_Y)"
    case.UpdatePipeline()
    case = Calculator(Input=case)
    case.ResultArrayName = "UVPOS"
    case.Function = " sqrt((UVPOS + abs(UVPOS))/2)"
    case.UpdatePipeline()
    arrayList.append('UVPOS')
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA :: APPEND CASE SPECIFIC VARIABLES")
  startTime = time.time()
  case = ProgrammableFilter(Input=case)
  case.Script = \
  """
  for var in arrayList:
   output.PointData.append(inputs[0].PointData[var],var)
  
  """
  case.UpdatePipeline()
  # POINT DATA TO CELL DATA 
  PDtCD1 = PointDatatoCellData(Input=case)
  case.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA :: CALCULATE CELL VOLUME")
  startTime = time.time()
  case = PythonCalculator(Input=PDtCD1)
  #case = PythonCalculator(Input=case)
  # Properties modified on pythonCalculator1
  if(DIM==3):
   if('pvd' in IF): 
    case.Expression = 'area(inputs[0])'
   else:
    case.Expression = 'volume(inputs[0])'
  elif(DIM==2):
   case.Expression = 'area(inputs[0])'
  else:
   raise ValueError("--|| ERROR: DIMESION NOT ACCURATE")
   
  case.ArrayName = 'Volume'
  case.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  # EXTRACT DATA AT SPECIFIC TIMES
  print("--|| ALYA :: EXTRACT TIME DATA")
  case = GroupTimeSteps(Input=case)  
  case.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  print("--|| ALYA :: PERFORM MODAL ANALYSIS AND APPEND DATA")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=case)
  PF1.Script = \
"""
np.set_printoptions(threshold=5000)
inp0 = dsa.WrapDataObject(inputs[0].GetBlock(0))
#inp0 = inputs[0].GetBlock(0)

print("----|| ALYA :: MANUPULATE DATA")
print("----|| INFO :: DIRECTORY OF OBJECT ", dir(inp0))
print("----|| INFO :: BLOCK 0 CELL ARRAYS ", inp0.CellData.keys())
print("----|| INFO :: BLOCK 0 POINT ARRAYS ", inp0.PointData.keys())
startTime = time.time()

if('pvd' in IF):
  V = dsa.WrapDataObject(inputs[0].GetBlock(0)).CellData['Volume']
else: 
  V = dsa.WrapDataObject(inputs[0].GetBlock(0)).CellData['Volume'].Arrays[0]

V = vtk_to_numpy(V)

Vtotal = sum(V)
V /= Vtotal

N = np.size(times)
L = np.size(V)

print('----|| ALYA : TEMPORAL SIZE IS ',N, 'SPATIAL SIZE IS ',L)
if field_is_vector == True:
 fields = np.zeros([N,L,len(wtsMat)])
else:
 fields = np.zeros([N,L])
print ('----|| ALYA : READING TIME: FROM {} TO {}'.format(times[0],times[-1]))
for i in range(N):
 field = np.empty([L,0],dtype='float')
 for fieldname in arrayList:
  t = times[i]
  d = dsa.WrapDataObject(inputs[0].GetBlock(i))
  if('pvd' in IF): 
    d = vtk_to_numpy(d.CellData[fieldname])
  else:  
    d = vtk_to_numpy(d.CellData[fieldname].Arrays[0])
  field = np.column_stack([field,d])
 if field_is_vector == False:
  field=field.T
 fields[i]=np.copy(field)
 
if subtractAvg:
 fields_avg = np.average(fields, axis=0)
 fields = fields - fields_avg
 print('--|| ALYA : CALCULATING FLUCTUATIONS')
print("----|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

## POD SECTION IMPLEMENTATION
if do_POD:
 print("----|| ALYA :: IMPLEMENT POD CALCULATIONS")
 startTime = time.time()
 # if field is a vector, reshape the fields and corresponding volument weight
 if field_is_vector:
     shp_vec = fields.shape
     shp_flat = (fields.shape[0],fields.shape[1]*fields.shape[2])
     fields = fields.reshape(shp_flat)
     V = np.dot(np.tile(wtsMat,(shp_vec[2],1)).T,np.tile(V,(shp_vec[2],1)))
     V = V.T.reshape(shp_flat[1])

 POD_res = mr.compute_POD_arrays_snaps_method(
    fields.T, inner_product_weights=V)

 modes = POD_res.modes; eigen_vals = POD_res.eigvals
 eigen_vecs = POD_res.eigvecs; correlation_mat = POD_res.correlation_array

 print(np.shape(modes),np.shape(eigen_vals),np.shape(eigen_vecs),np.shape(correlation_mat))

 print("----|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

 # if field is a vector, reshape the output matrix
 if field_is_vector:
   fields = fields.reshape(shp_vec)
   modes = np.asarray(modes).T.reshape((modes.shape[1],shp_vec[1],shp_vec[2]))
   V = np.asarray(V).T.reshape([shp_vec[1],shp_vec[2]])
 else:
   modes = np.asarray(modes).T

 # WRITE OR APPEND THE DATA REQUIRED INTO PARAVIEW
 if output_POD_spatial_modes:
   print("----|| ALYA :: APPENDING SPATIAL MODE")
   sclArr = 0*dsa.WrapDataObject(inputs[0].GetBlock(0)).CellData['Volume']
   vecArr = 0*dsa.WrapDataObject(inputs[0].GetBlock(0)).CellData[arrayList[1]]
   if field_is_vector:
     ### ADD WEIGTH VECTORS FOR PRESSURE
     ##if('pvd' in IF): 
     ##  sclArr = V[:,0]
     ##else:
     ##  sclArr.Arrays[0] = V[:,0]
     ##output.CellData.append(sclArr,"WPRES")
     ### ADD WEIGTH VECTORS FOR VELOCITY
     ##if('pvd' in IF): 
     ##  vecArr = V[:,1:DIM]
     ##else:
     ##  vecArr.Arrays[0] = V[:,1:DIM]
     ##output.CellData.append(vecArr,"WVELO")
     ##if("RS_II" in varName_calc):
     ## # ADD WEIGTH VECTORS FOR RS_IJ
     ##  if('pvd' in IF): 
     ##    vecArr = V[:,4:6]
     ##  else:
     ##    vecArr.Arrays[0] = V[:,4:6]
     ##  output.CellData.append(vecArr,"WRSII")
     ##if("RS_IJ" in varName_calc):
     ## # ADD WEIGTH VECTORS FOR RS_IJ
     ##  if('pvd' in IF): 
     ##    vecArr = V[:,7:9]
     ##  else:
     ##    vecArr.Arrays[0] = V[:,7:9]
     ##  output.CellData.append(vecArr,"WRSIJ")

     # ADD AVERAGES FOR ALL VARIABLES 
     if subtractAvg:
       if('pvd' in IF): 
         sclArr = fields_avg[:,0]
         vecArr = fields_avg[:,1:DIM]
       else:
         sclArr.Arrays[0] = fields_avg[:,0]
         vecArr.Arrays[0] = fields_avg[:,1:DIM]
       output.CellData.append(sclArr,"AVPRE")
       output.CellData.append(vecArr,"AVVEL")
       if("RS_II" in varName_calc):
         if('pvd' in IF): 
           vecArr = fields_avg[:,4:6]
         else:
           vecArr.Arrays[0] = fields_avg[:,4:6]
         output.CellData.append(vecArr,"ARSII")
       if("RS_IJ" in varName_calc):
         if('pvd' in IF): 
           vecArr = fields_avg[:,7:9]
         else:
           vecArr.Arrays[0] = fields_avg[:,7:9]
         output.CellData.append(vecArr,"ARSIJ")
     
     if(M_POD>np.shape(modes)[0]):
       raise ValueError("--|| ERROR: YOU HAVE MORE MODES THAN SNAPS")

     # ADD MODES FOR ALL VARIABLES 
     for i in range(M_POD):
       # ADD PRESSURE MODE
       if('pvd' in IF): 
         sclArr = modes[i,:,0]
       else:
         sclArr.Arrays[0] = modes[i,:,0]
       output.CellData.append(sclArr,'POD_mode_{}_{}'.format("PRESS",i))
       # ADD VELOCITY MODE
       if('pvd' in IF): 
         vecArr = modes[i,:,1:DIM]
       else:
         vecArr.Arrays[0] = modes[i,:,1:DIM]
       output.CellData.append(vecArr,'POD_mode_{}_{}'.format("VELOC",i))
       if("RS_II" in varName_calc):
         # ADD RS_II MODE
         if('pvd' in IF): 
           vecArr = modes[i,:,4:6]
         else:
           vecArr.Arrays[0] = modes[i,:,4:6]
         output.CellData.append(vecArr,'POD_mode_{}_{}'.format("RS-II",i))
       if("RS_IJ" in varName_calc):
         # ADD RS_IJ MODE
         if('pvd' in IF): 
           vecArr = modes[i,:,7:9]
         else:
           vecArr.Arrays[0] = modes[i,:,7:9]
         output.CellData.append(vecArr,'POD_mode_{}_{}'.format("RS-IJ",i))
   else:    
     # ADD SCALARS
     if('pvd' in IF): 
      sclArr = V
     else:
      sclArr.Arrays[0] = V
     output.CellData.append(sclArr,"Volume")
     if subtractAvg:
       if('pvd' in IF): 
        sclArr = fields_avg
       else:
        sclArr.Arrays[0] = fields_avg
       output.CellData.append(sclArr,varName.upper())
     for i in range(M_POD):
       if('pvd' in IF): 
        sclArr = modes[i]
       else:
        sclArr.Arrays[0] = modes[i]
       output.CellData.append(sclArr,'POD_mode_{}_{}'.format(varName0,i))
 if output_correlation_matrix:
   print("----|| ALYA :: WRITING FIELD CORRELATION MATRIX IN",POD_cm_filename)
   np.savetxt(POD_cm_filename,correlation_mat,delimiter=',')

 if output_POD_temporal_modes:
   print("----|| ALYA :: WRITING POD TEMPORAL MODES IN",POD_tm_filename)
   singular_vals = eigen_vals**0.5
   POD_mode_energy_normalized = eigen_vals/correlation_mat.trace()
   cumsum_POD_mode_energy_normalized = np.cumsum(POD_mode_energy_normalized)
   header_str = 'TEMPORAL MODES\\n'
   header_str += 'Modes =,  '
   header_str += np.array2string(np.asarray((range(N-1)), dtype=int),separator=', ')
   header_str += '\\n'
   header_str += 'SV =,  '
   header_str += np.array2string(np.asarray(singular_vals, dtype=float),separator=', ')
   header_str += '\\n'
   header_str += 'EV =,  '
   header_str += np.array2string(np.asarray(eigen_vals, dtype=float),separator=', ')
   header_str += '\\n'
   header_str += 'NEnergy =,  '
   header_str += np.array2string(np.asarray(POD_mode_energy_normalized, dtype=float),separator=', ')
   header_str += '\\n'
   header_str += 'CEnergy =,  '
   header_str += np.array2string(np.asarray(cumsum_POD_mode_energy_normalized, dtype=float),separator=', ')
   header_str += '\\n'
   header_str += 'TIME, EIGENVEC\\n'
   np.savetxt(POD_tm_filename,np.c_[times,eigen_vecs],delimiter = ', ',header = header_str)

 #--------Time stamp reconstruction of modes----------------#
 if doReconstruction:
   print("----|| ALYA :: RECONSTRUCTING WITH {} MODES AT TIME = {}".format(MR_POD,ReconTime))
   # reconstruct from first MR_POD POD modes
   ReconN = np.searchsorted(times,ReconTime)
   print("----|| ALYA :: RECONSTRUCTING AT TIME {} RATHER THAN AT {}".format(times[ReconN],ReconTime))
   if(DIM==3):
    recon_field = np.einsum("ijk...,i,i->jk",modes[:MR_POD],eigen_vals[:MR_POD]**0.5,np.asarray(eigen_vecs)[ReconN,:MR_POD])
   elif(DIM==2):
    recon_field = np.einsum("ij...,i,i->j",modes[:MR_POD],eigen_vals[:MR_POD]**0.5,np.asarray(eigen_vecs)[ReconN,:MR_POD])
   else: 
    raise ValueError("--|| ERROR: DIMESION NOT ACCURATE")

   if subtractAvg:
     recon_field = recon_field + fields_avg
   if field_is_vector:
    if('pvd' in IF): 
     vecArr = recon_field
    else:
     vecArr.Arrays[0] = recon_field
    output.CellData.append(vecArr,'POD_{}_RECON_{}_{}'.format(MR_POD,fieldname,ReconTime))
   else: 
    if('pvd' in IF): 
     sclArr = recon_field
    else:
     sclArr.Arrays[0] = recon_field
    output.CellData.append(sclArr,'POD_{}_RECON_{}_{}'.format(MR_POD,fieldname,ReconTime))
 #--------Full reconstruction of modes----------------#
 if doFullRecon:
   print("----|| ALYA :: FULL RECONSTRUCTION WITH {} MODES ".format(MR_POD))
   # reconstruct from first MR_POD POD modes
   if(DIM==3):
    recon_field = np.einsum("ijk...,i,i->jk",modes[:MR_POD],eigen_vals[:MR_POD]**0.5,np.asarray(eigen_vecs)[:,:MR_POD])
   elif(DIM==2):
    recon_field = np.einsum("ij...,i,i->j",modes[:MR_POD],eigen_vals[:MR_POD]**0.5,np.asarray(eigen_vecs)[:,:MR_POD])
   else: 
    raise ValueError("--|| ERROR: DIMESION NOT ACCURATE")

   #if subtractAvg:
   #  recon_field = recon_field + fields_avg
   if field_is_vector:
    if('pvd' in IF): 
     vecArr = recon_field
    else:
     vecArr.Arrays[0] = recon_field
    output.CellData.append(vecArr,'POD_{}_FRECON_{}'.format(MR_POD,fieldname))
   else: 
    if('pvd' in IF): 
     sclArr = recon_field
    else:
     sclArr.Arrays[0] = recon_field
    output.CellData.append(sclArr,'POD_{}_FRECON_{}'.format(MR_POD,fieldname))
#-----------------------------#     
#-------DMD part--------------#     
#-----------------------------#     
if do_DMD:
 print("----|| ALYA :: IMPLEMENT DMD CALCULATIONS")
 startTime = time.time()
 # if field is a vector, reshape the fields and corresponding volume weight
 if field_is_vector:
   shp_vec = fields.shape
   shp_flat = (fields.shape[0],fields.shape[1]*fields.shape[2])
   fields = fields.reshape(shp_flat)
   V = np.dot(np.tile(wtsMat,(shp_vec[2],1)).T,np.tile(V,(shp_vec[2],1)))
   V = V.T.reshape(shp_flat[1])

 # Compute DMD modes using method of snapshot
 DMD_res = mr.compute_DMD_arrays_snaps_method(fields.T,inner_product_weights=V)
 print("----|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

 eigen_vals = DMD_res.eigvals
 spectral_coeffs = DMD_res.spectral_coeffs
 exact_modes = DMD_res.exact_modes; 
 proj_modes = DMD_res.proj_modes; 
 adjoint_modes = DMD_res.adjoint_modes; 
 proj_coeffs = DMD_res.proj_coeffs
 adv_proj_coeffs = DMD_res.adv_proj_coeffs
 R_low_order_eigvecs = DMD_res.R_low_order_eigvecs
 L_low_order_eigvecs = DMD_res.L_low_order_eigvecs
 correlation_array_eigvals = DMD_res.correlation_array_eigvals
 correlation_array_eigvecs = DMD_res.correlation_array_eigvecs
 correlation_array = DMD_res.correlation_array
 cross_correlation_array = DMD_res.cross_correlation_array


 # if field is a vector, reshape the fields, V and output matrix
 modes = np.real(proj_modes)

 if field_is_vector:
   fields = fields.reshape(shp_vec)
   modes = np.asarray(modes).T.reshape((modes.shape[1],shp_vec[1],shp_vec[2]))
   V = np.asarray(V).T.reshape([shp_vec[1],shp_vec[2]])
 else:
   modes = np.asarray(modes).T

 #---Calculate additional information----------# 
 mode_amp = np.linalg.lstsq(np.einsum('ij,jk->ik', modes.T ,np.diag(eigen_vals)), \
            fields[0,:].T)  #using first snapshot

 mode_amp = abs(mode_amp[0])
 #mode_norms = np.linalg.norm(proj_modes,axis=0)
 #print(np.shape(spectral_coeffs),spectral_coeffs)
 #print(np.shape(proj_coeffs),proj_coeffs)

 ### sorting modes based on norm of the modes
 #eorder = np.argsort(mode_norms)[::-1]
 ### re-order the outputs
 #eigen_vals = eigen_vals[eorder]
 #mode_norms = mode_norms[eorder]
 #proj_coeffs = proj_coeffs[:,eorder]
 ###build the DMD_modes
 ##DMD_modes = np.einsum('ijk,il->ljk', fields,build_coeffs[:,:M_DMD])

 if output_DMD_info:
   print("----|| ALYA :: WRITING DMD INFO TO ",DMD_info_filename)
   # output modes info
   header_str = 'DMD modes info\\n'
   header_str += 'eig_real, eig_imag, growth_rate, frequency, spectral_coeffs, mode_amp\\n'
   header_str += r'AU, AU, 1/s, Hz, AU, AU'
   dt = np.average(np.array(times[1:])-np.array(times[:-1])) #time step
   np.savetxt(DMD_info_filename, \
               np.c_[ np.real(eigen_vals), \
                   np.imag(eigen_vals), \
                   np.log(np.abs(eigen_vals))/dt, \
                   np.angle(eigen_vals)/dt, \
                   spectral_coeffs,\
                   mode_amp], \
               delimiter = ',', \
               header = header_str)

 if output_DMD_build_coeffs:
   print("----|| ALYA :: WRITING DMD BUILD COEFF TO ",DMD_build_coeffs_filename)
   np.savetxt(DMD_build_coeffs_filename, proj_coeffs,header = 'proj_coeffs',delimiter = ',')
   with open(DMD_build_coeffs_filename, "ab") as f:
     np.savetxt(f, adv_proj_coeffs,header = 'adv_proj_coeffs',delimiter = ',')

 if output_DMD_spatial_modes:
   print("----|| ALYA :: APPENDING SPATIAL MODE")
   sclArr = 0*dsa.WrapDataObject(inputs[0].GetBlock(0)).CellData['Volume']
   if field_is_vector:
     vecArr = 0*dsa.WrapDataObject(inputs[0].GetBlock(0)).CellData[arrayList[1]]
     # ADD WEIGTH VECTORS FOR PRESSURE
     if('pvd' in IF): 
       sclArr = V[:,0]
     else:
       sclArr.Arrays[0] = V[:,0]
     output.CellData.append(sclArr,"WPRES")
     # ADD WEIGTH VECTORS FOR VELOCITY
     if('pvd' in IF): 
       vecArr = V[:,1:DIM]
     else:
       vecArr.Arrays[0] = V[:,1:DIM]
     output.CellData.append(vecArr,"WVELO")
     if("RS_II" in varName_calc):
      # ADD WEIGTH VECTORS FOR RS_IJ
       if('pvd' in IF): 
         vecArr = V[:,4:6]
       else:
         vecArr.Arrays[0] = V[:,4:6]
       output.CellData.append(vecArr,"WRSII")
     if("RS_IJ" in varName_calc):
      # ADD WEIGTH VECTORS FOR RS_IJ
       if('pvd' in IF): 
         vecArr = V[:,7:9]
       else:
         vecArr.Arrays[0] = V[:,7:9]
       output.CellData.append(vecArr,"WRSIJ")

     # ADD AVERAGES FOR ALL VARIABLES 
     if subtractAvg:
       if('pvd' in IF): 
         sclArr = fields_avg[:,0]
         vecArr = fields_avg[:,1:DIM]
       else:
         sclArr.Arrays[0] = fields_avg[:,0]
         vecArr.Arrays[0] = fields_avg[:,1:DIM]
       output.CellData.append(sclArr,"AVPRE")
       output.CellData.append(vecArr,"AVVEL")
       if("RS_II" in varName_calc):
         if('pvd' in IF): 
           vecArr = fields_avg[:,4:6]
         else:
           vecArr.Arrays[0] = fields_avg[:,4:6]
         output.CellData.append(vecArr,"ARSII")
       if("RS_IJ" in varName_calc):
         if('pvd' in IF): 
           vecArr = fields_avg[:,7:9]
         else:
           vecArr.Arrays[0] = fields_avg[:,7:9]
         output.CellData.append(vecArr,"ARSIJ")
     for i in range(M_DMD):
       # ADD PRESSURE MODE
       if('pvd' in IF): 
         sclArr = modes[i,:,0]
       else:
         sclArr.Arrays[0] = modes[i,:,0]
       output.CellData.append(sclArr,'DMD_mode_{}_{}'.format("PRESS",i))
       # ADD VELOCITY MODE
       if('pvd' in IF): 
         vecArr = modes[i,:,1:DIM]
       else:
         vecArr.Arrays[0] = modes[i,:,1:DIM]
       output.CellData.append(vecArr,'DMD_mode_{}_{}'.format("VELOC",i))
       if("RS_II" in varName_calc):
         # ADD RS_II MODE
         if('pvd' in IF): 
           vecArr = modes[i,:,4:6]
         else:
           vecArr.Arrays[0] = modes[i,:,4:6]
         output.CellData.append(vecArr,'DMD_mode_{}_{}'.format("RS-II",i))
       if("RS_IJ" in varName_calc):
         # ADD RS_IJ MODE
         if('pvd' in IF): 
           vecArr = modes[i,:,7:9]
         else:
           vecArr.Arrays[0] = modes[i,:,7:9]
         output.CellData.append(vecArr,'DMD_mode_{}_{}'.format("RS-IJ",i))
   else:    
     # ADD SCALARS
     if('pvd' in IF): 
      sclArr = V
     else:
      sclArr.Arrays[0] = V
     output.CellData.append(sclArr,"Volume")
     if subtractAvg:
       if('pvd' in IF): 
        sclArr = fields_avg
       else:
        sclArr.Arrays[0] = fields_avg
       output.CellData.append(sclArr,"AV"+str(arrayList[0:3]))
     for i in range(M_DMD):
       if('pvd' in IF): 
        sclArr = modes[i]
       else:
        sclArr.Arrays[0] = modes[i]
       output.CellData.append(sclArr,'DMD_mode_{}_{}'.format(arrayList[0],i))
#   print("----|| ALYA :: APPEND DMD SPATIAL MODES TO ",DMD_sm_filename)
#   #output to xml vtk unstructured grid file
#   ugcd = geom.GetCellData()
#   ugcd.Reset()
#   ugcd.CopyAllOff()
#   for i in range(ugcd.GetNumberOfArrays()):
#       ugcd.RemoveArray(0)
#   #import pi
#   from numpy import pi
#
#   for i in range(M_DMD):
#       ugcd.AddArray(dsa.numpyTovtkDataArray(np.abs(DMD_modes[i]),prefix+'_DMD_mode_abs_{}_{}'.format(    fieldname,i)))
#       ugcd.AddArray(dsa.numpyTovtkDataArray(np.angle(DMD_modes[i])*180/pi,prefix+'_DMD_mode_angle_{}_    {}'.format(fieldname,i)))
#
#   ugw = vtk.vtkXMLUnstructuredGridWriter()
#   ugw.SetInputDataObject(geom)
#   ugw.SetFileName(DMD_sm_filename)
#   ugw.Write()


""" 
  
  PF1.UpdatePipeline()
  # create a new 'Extract Block'
  EB1 = ExtractBlock(Input=PF1)
  # Properties modified on extractBlock1
  EB1.BlockIndices = [1]
  EB1.UpdatePipeline()
  
  C2P = CellDatatoPointData(Input=EB1)
  C2P.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  
  
  #print("----|| ALYA :: WRITING FILE", POD_sm_filename)
  if(do_POD):
   savePath = POD_sm_filename
   SaveData(savePath, proxy=C2P)
  if(do_DMD):
   savePath = DMD_sm_filename
   SaveData(savePath, proxy=C2P)
if('RECON' in MODE):
  print("--|| ALYA :: READING POD ARRAYS")
  startTime = time.time()
  case = OpenDataFile(casename+'.vtm')
  case.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  print("--|| ALYA :: MERGING BLOCKS")
  startTime = time.time()
  case = MergeBlocks(Input=case)
  case.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  print("--|| ALYA :: PERFORM MODAL RECONSTRUCTION AND WRITE DATA")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=case)
  PF1.Script = \
"""
import os 
import vtk
import numpy as np
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from csv import reader

#--Load eigenvectors----#
casePath = os.getcwd()
fileName = casePath+'/POD_temporal_coefficients.csv'

headerStr = [];
with open(fileName, "r") as csv_file:
  csv_reader = reader(csv_file)
  count = 0
  for row in csv_reader:
    if(row[0].startswith('# Modes')):
      indModes = count
    if(row[0].startswith('# SV')):
      indSV = count
    if(row[0].startswith('# EV')):
      indEV = count
    if(row[0].startswith('# NEnergy')):
      indEnergy = count
    if(row[0].startswith('#')):
      headerStr.append(row)
    count =count +1;  

#-----Extract modes array -------#
headerArray = list(np.concatenate(headerStr[indModes:indSV]).flat)
headerArray = list(filter(None, headerArray))

headerArray = [i.replace('#', '') for i in headerArray]
headerArray = [i.replace('[', '') for i in headerArray]
headerArray = [i.replace(']', '') for i in headerArray]

headerArray = headerArray[1:]

headerArray = [int(i) for i in headerArray]
ModesArray = np.array(headerArray, dtype=int)

print('--||INFO: TOTAL MODES = %d'%len(ModesArray))
#-----Extract eigenvalues array -------#
headerArray = list(np.concatenate(headerStr[indEV:indEnergy]).flat)
headerArray = list(filter(None, headerArray))

headerArray = [i.replace('#', '') for i in headerArray]
headerArray = [i.replace('[', '') for i in headerArray]
headerArray = [i.replace(']', '') for i in headerArray]

headerArray = headerArray[1:]

headerArray = [float(i) for i in headerArray]
eigenVal = np.array(headerArray, dtype=int)

#-------------------------------------#
#-------------------------------------#

data = np.loadtxt(fileName, dtype=float, comments='#', delimiter=',')
timeArray = data[:,0]; 
nt = len(timeArray); 
eigenVec = data[:,1:];

#-------------------------------------#
#-------------------------------------#
var_recon = ['VELOC','PRESS']

outFolderName = 'PodData_2D_{}'.format(MR_POD)
if not os.path.exists(outFolderName):
  os.makedirs(outFolderName)

d = dsa.WrapDataObject(inputs[0].VTKObject)
input = self.GetInputDataObject(0,0)
data_to_write = vtk.vtkUnstructuredGrid()
data_to_write.CopyStructure(input)
dataSet = self.GetOutputDataObject(0)
dataSet.CopyStructure(input)
ugw = vtk.vtkXMLUnstructuredGridWriter()
ugw.SetInputData(data_to_write)
#ugw.SetNumberOfTimeSteps(nt)
#ugw.WriteNextTime(n)
#ugw.Stop()

#-----------------------#
f = open(outFolderName+'.pvd', "w")
f.write('<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">\\n')
f.write('\\t<Collection>\\n')

for n in range(nt):
 for pod_var in var_recon:
  varName = 'POD_mode_'+pod_var+'_'+str(0);
  avvel = d.PointData["AV"+pod_var[0:3]]
  avgFlowRecon = 0.0*d.PointData[varName]
  avgFlow = 0.0*d.PointData[varName]
  for m in range(MR_POD):
    varName = 'POD_mode_'+pod_var+'_'+str(m);
    avgFlow = avgFlow+(d.PointData[varName]*eigenVal[m]**0.5*eigenVec[n,m])
  avgFlowRecon = avvel+avgFlow
  vtk_array = numpy_to_vtk(num_array=avgFlowRecon, deep=True, array_type=vtk.VTK_FLOAT)
  vtk_array.SetName(pod_var)
  dataSet.GetPointData().AddArray(vtk_array)
 data_to_write.ShallowCopy(dataSet)
 fileWriteName = casePath+'/'+outFolderName+'/PodData_2D_{}.vtu'.format(n)
 ugw.SetFileName(fileWriteName)
 ugw.Write()
 f.write('\\t\\t<DataSet timestep="%f" part="0" file="%s"/>\\n' \
        % (timeArray[n],outFolderName+'/PodData_2D_{}.vtu'.format(n)))
f.write('\\t</Collection>\\n')
f.write('</VTKFile>')
f.close()
""" 
  PF1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
