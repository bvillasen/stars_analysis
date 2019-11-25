import os, sys
from os import listdir
from os.path import isfile, join
from data_compress_grid import compress_grid
from tools import create_directory
import numpy as np
import time


dataDir = '/data/groups/comp-astro/bruno/data/1024_cool_uv_50Mpc/'
inDir = dataDir + 'out_files/'
outDir = dataDir + 'snapshots/'

# dataDir = '/raid/bruno/data/stars/polytrope/'
# inDir = dataDir + 'output_files/'
# outDir = dataDir + 'snapshots/grav_0./'
# 



def split_name( file_name):
  nSap, name, nBox = file_name.split('.')
  return [int(nSap), int(nBox)]


print( 'Input Dir: ' + inDir )
print( 'Output Dir: ' + outDir )
create_directory( outDir )
print("")

name_base = 'h5'
dataFiles = [f for f in listdir(inDir) if (isfile(join(inDir, f)) and (f.find('.h5.') > 0 ) and ( f.find('_particles') < 0) ) ]
dataFiles = np.sort( dataFiles )
nFiles = len( dataFiles )

files_names = np.array([ split_name( file_name ) for file_name in dataFiles ])
snaps, boxes = files_names.T
snapshots_all = np.unique( snaps )
boxes = np.unique( boxes )
snapshots_all.sort()
nSnapshots = len( snapshots_all )
nBoxes = len( boxes )

print( "Number of snapshots: {0}".format(nSnapshots) )
print( "Number of files per snapshot: {0}".format(nBoxes) )


#Set wich snapshots to compress
snapshots_to_compress = range( nSnapshots )
print( "\nNumber of snapshots to compres: {0}".format(len(snapshots_to_compress)) )

#available Hydro Fields:
#[ density, momentum_x, momentum_y, momentum_z, Energy, GasEnergy ]
# hydro_fields = 'all'
hydro_fields = ['density', 'momentum_x', 'momentum_y', 'momentum_z', 'Energy' ]
print( "\nHydro fields: {0}".format(hydro_fields))


precision = np.float64
# precision = np.float32
# precision = np.float16
print( "\nPrecision: {0}".format( precision ))

print( "\nCompressing Snapshots..." )
for nSnap in snapshots_to_compress:
  start = time.time()
  out_base_name = 'grid_' 
  compress_grid( nSnap, nBoxes, name_base, out_base_name, inDir, outDir, hydro_fields,  precision=precision )
  end = time.time()
  print( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) )

