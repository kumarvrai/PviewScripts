# trace generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
src = FindSource('vortices_0.dat')

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(registrationName='TableToPoints1', Input=src)

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'xcenter'
tableToPoints1.YColumn = 'ycenter'
tableToPoints1.ZColumn = 'zcenter'
tableToPoints1.a2DPoints = 1
tableToPoints1.UpdatePipeline()

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=tableToPoints1)
# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
# rename source object
RenameSource('MAKE-SPHERES', glyph1)
glyph1.ScaleArray = ['POINTS', 'radius']
glyph1.ScaleFactor = 2.0
glyph1.GlyphMode = 'All Points'
glyph1.UpdatePipeline()

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=tableToPoints1)

# rename source object
RenameSource('COORD-CALC', calculator1)
# Properties modified on calculator1
calculator1.CoordinateResults = 1
calculator1.Function = 'coordsX*kHat+coordsY*jHat+coordsZ*iHat'
calculator1.UpdatePipeline()


# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=glyph1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.0, 1.0014573223888874, 1.551501527428627]
# rename source object
RenameSource('MAKE-CIRCLES', slice1)
slice1.UpdatePipeline()

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
# change solid color
slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
# Properties modified on slice1Display
slice1Display.LineWidth = 3.0
# Properties modified on slice1Display
slice1Display.RenderLinesAsTubes = 1

