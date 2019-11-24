import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np



def load_snapshot_cholla( nSnap, inDir, single_file=False  ):
  gridFileName = inDir + 'grid_{0:03}.h5'.format(nSnap)
  if single_file: gridFileName = inDir + '{0}.h5'.format(nSnap)
  
  outDir = { 'gas':{} }
  
  data_grid = h5.File( gridFileName, 'r' )
  fields_data = data_grid.keys()
  if single_file: 
    for key in data_grid.attrs.keys(): outDir[key] = data_grid.attrs[key][0]
  else: 
    for key in data_grid.attrs.keys(): outDir[key] = data_grid.attrs[key]
  fields_grid = fields_data
  for field in fields_grid:
    if field not in fields_data: continue
    outDir['gas'][field] = data_grid[field]

  return outDir
