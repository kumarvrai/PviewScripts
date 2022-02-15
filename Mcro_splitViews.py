import numpy as np
from paraview.simple import *
 
 def splitView(ROWS,COLS):
 # Usage: splitView('ROWS, COLS)
  layOut1 = GetLayout()
  # for i in range(ROWS-1):
  #  view = GetActiveView()
  #  layOut1.SplitViewHorizontal(view, float(1.0/(ROWS-i)))
  #  rv = CreateView('RenderView')

  for i in range(COLS-1):
    view1 = GetActiveView()
    num = float('{0:.5f}'.format(1.0/(COLS-i)))
    layOut1.SplitViewVertical(view1, num)
    rv = CreateView('RenderView')
	       

 splitView(1,3)
