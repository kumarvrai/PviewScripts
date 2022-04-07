initializeVars = '/home/u/ugo/kvishal/1.post_process/0.alya_pv_scripts/initializeVariables.py'
import sys
from paraview.simple import *
pyVer = sys.version_info[0]
if pyVer < 3:
  execfile(initializeVars)
else:
  exec(open(initializeVars).read())
# LOAD SNAPSHOTS
print('--|| ALYA :: INITIALIZING PYTHON VERSION:',pyVer)
niaHome = '/home/u/ugo/kvishal/'
niaScrh = '/scratch/u/ugo/kvishal/research/0.Alya/'

if read_fields_from_file:
  print("--|| ALYA :: READING ALYA ARRAYS", arrayList)
  startTime = time.time()
  case = OpenDataFile(ID+IF)
  case.PointArrays = arrayList
  case.UpdatePipeline()
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
    #print("--|| ALYA :: TAKING DATA SLICE")
    #startTime = time.time()
    #case = Slice(Input=case)
    #case.SliceType = 'Plane'
    ##slice.SliceType.Origin = [0.0, 0.0, 0.0]
    #case.SliceType.Normal = [0.0, 0.0, 1.0]
    #case.UpdatePipeline()
    #print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')
  # TEMPORAL AVERAGING
  print("--|| ALYA :: TEMPORAL AVERAGING OF CODE ARRAYS")
  startTime = time.time()
  case2 = TemporalStatistics(Input=case)
  case2.ComputeMinimum = 0
  case2.ComputeMaximum = 0
  case2.ComputeStandardDeviation = 0
  case2.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  # APPEND DATASETS
  print("--|| ALYA :: APPEND DATASETS")
  startTime = time.time()
  case = AppendAttributes(Input=[case,case2])
  case.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  # CALCULATE EXTRA VARIABLES
  print("--|| ALYA :: CALCULATING EXTRA VARIABLES")
  startTime = time.time()

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
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print "--|| ALYA :: APPEND CASE SPECIFIC VARIABLES"
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
  PC1 = PythonCalculator(Input=PDtCD1)
  # Properties modified on pythonCalculator1
  if(DIM==3):
   if('pvd' in IF): 
    PC1.Expression = 'area(inputs[0])'
   else:
    PC1.Expression = 'volume(inputs[0])'
  elif(DIM==2):
   PC1.Expression = 'area(inputs[0])'
  else:
   raise ValueError("--|| ERROR: DIMESION NOT ACCURATE")
   
  PC1.ArrayName = 'Volume'
  PC1.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  # EXTRACT DATA AT SPECIFIC TIMES
  print("--|| ALYA :: EXTRACT TIME DATA")
  dataSet = GroupTimeSteps(Input=PC1)  
  dataSet.UpdatePipeline()
  print("--|| ALYA :: DONE. TIME =",time.time()-startTime,'sec')

  print("--|| ALYA :: PERFORM MODAL ANALYSIS AND APPEND DATA")
  startTime = time.time()
  PF1 = ProgrammableFilter(Input=dataSet)
  PF1.Script = \
  """
print("----|| ALYA :: MANUPULATE DATA")
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
    fields.T,list(mr.range(M)), inner_product_weights=V)

 modes = POD_res.modes; eigen_vals = POD_res.eigvals
 eigen_vecs = POD_res.eigvecs; correlation_mat = POD_res.correlation_array

 #print(np.shape(modes),np.shape(eigen_vals),np.shape(eigen_vecs),np.shape(correlation_mat))

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
     # ADD WEIGTH VECTORS FOR PRESSURE
     if('pvd' in IF): 
       sclArr = V[:,0]
     else:
       sclArr.Arrays[0] = V[:,0]
     output.CellData.append(sclArr,"WPRES")
     # ADD WEIGTH VECTORS FOR VELOCITY
     if('pvd' in IF): 
       vecArr = V[:,1:3]
     else:
       vecArr.Arrays[0] = V[:,1:3]
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
         vecArr = fields_avg[:,1:3]
       else:
         sclArr.Arrays[0] = fields_avg[:,0]
         vecArr.Arrays[0] = fields_avg[:,1:3]
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
     for i in range(M):
       # ADD PRESSURE MODE
       if('pvd' in IF): 
         sclArr = modes[i,:,0]
       else:
         sclArr.Arrays[0] = modes[i,:,0]
       output.CellData.append(sclArr,'POD_mode_{}_{}'.format("PRESS",i))
       # ADD VELOCITY MODE
       if('pvd' in IF): 
         vecArr = modes[i,:,1:3]
       else:
         vecArr.Arrays[0] = modes[i,:,1:3]
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
     for i in range(M):
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

 if doReconstruction:
   print("----|| ALYA :: RECONSTRUCTING WITH {} MODES AT TIME = {}".format(MR,ReconTime))
   # reconstruct from first MR POD modes
   ReconN = np.searchsorted(times,ReconTime)
   print("----|| ALYA :: RECONSTRUCTING AT TIME {} RATHER THAN AT {}".format(times[ReconN],ReconTime))
   if(DIM==3):
    recon_field = np.einsum("ijk...,i,i->jk",modes[:MR],eigen_vals[:MR]**0.5,np.asarray(eigen_vecs)[ReconN,:MR])
   elif(DIM==2):
    recon_field = np.einsum("ij...,i,i->j",modes[:MR],eigen_vals[:MR]**0.5,np.asarray(eigen_vecs)[ReconN,:MR])
   else: 
    raise ValueError("--|| ERROR: DIMESION NOT ACCURATE")

   if subtractAvg:
     recon_field = recon_field + fields_avg
   if field_is_vector:
    if('pvd' in IF): 
     vecArr = recon_field
    else:
     vecArr.Arrays[0] = recon_field
    output.CellData.append(vecArr,'POD_{}_RECON_{}_{}'.format(MR,fieldname,ReconTime))
   else: 
    if('pvd' in IF): 
     sclArr = recon_field
    else:
     sclArr.Arrays[0] = recon_field
    output.CellData.append(sclArr,'POD_{}_RECON_{}_{}'.format(MR,fieldname,ReconTime))

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
  savePath = POD_sm_filename
  SaveData(savePath, proxy=C2P)
#else:
#  print('--|| ALYA :: READING FILE',ID+IF)
#  ofr = vtk.vtkEnSightGoldReader()
#  ofr.SetFileName(ID+IF)
#  ofr.SetPointArrayStatus(fieldname,1)
#  ofr.Update()
## Create random data
#num_vecs = 30
#vecs = np.random.random((100, num_vecs))
#
## Compute POD
#num_modes = 5
#POD_res = mr.compute_POD_arrays_snaps_method(
#    vecs, list(mr.range(num_modes)))
#modes = POD_res.modes
#eigvals = POD_res.eigvals

