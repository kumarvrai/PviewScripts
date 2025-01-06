import numpy as np
import pyAlya

BASEDIR = './'
CASESTR = 'AVG_naca'
VARLIST = ['avvel', 'avmueff', 'avpre']

pyAlya.pprint(0, 'Initializing!', flush=True)
mesh  = pyAlya.MeshSOD2D.read(CASESTR, basedir=BASEDIR)
pyAlya.pprint(0, 'Mesh read', flush=True)
field = pyAlya.FieldSOD2D.read(CASESTR, VARLIST, 1, mesh.xyz, basedir=BASEDIR)
pyAlya.pprint(0, 'Field read', flush=True)

field['gradv'] = mesh.gradient(field['avvel'])
pyAlya.pprint(0, 'Gradient computed', flush=True)

bcmesh, mask = mesh.extract_bc(1)
pyAlya.pprint(0, 'Mesh extracted', flush=True)
fieldbc = field.selectMask(mask)
pyAlya.pprint(0, 'Field extracted', flush=True)

fieldbc['norm']  = bcmesh.computeNormals()
pyAlya.pprint(0, 'Normals computed', flush=True)
fieldbc['dudotn'] = pyAlya.math.tensVecProd(fieldbc['gradv'], fieldbc['norm'])
pyAlya.pprint(0, 'Tensor-vector product computed', flush=True)
fieldbc['tang']  = np.transpose(np.array([fieldbc['norm'][:,1], -fieldbc['norm'][:,0], fieldbc['norm'][:,2]]))
pyAlya.pprint(0, 'Tangent computed', flush=True)
dummy            = pyAlya.math.dot(fieldbc['dudotn'], fieldbc['tang']) 
fieldbc['AVGCF']    = 2*fieldbc['avmueff']*dummy
pyAlya.pprint(0, 'Cf computed', flush=True)
fieldbc['AVGCP']    = 2*fieldbc['avpre']

bcmesh.write('output')
fieldbc.write('output')
