import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

stars_dir = os.path.dirname(os.getcwd()) + '/'
loadDataDirectory = stars_dir + "load_data/"
toolsDirectory = stars_dir + "tools/"
analysisDirectory = stars_dir + "analysis/"
figuresDirectory = stars_dir + "figures/"
sys.path.extend([ loadDataDirectory, toolsDirectory, analysisDirectory ] )

from load_data_cholla import load_snapshot_cholla

dataDir = '/raid/bruno/data/stars/polytrope/'
inputDir = dataDir + 'snapshots/' 

n_snapshots = 11

for n_snapshot in range( n_snapshots ):
  data = load_snapshot_cholla( n_snapshot, inputDir )

  dens = data['gas']['density'][...]
  nz, ny, nx = dens.shape
  indx_slice = nz/2
  dens_slice = dens[indx_slice, :, : ]

  fig, ax = plt.subplots()

  fig.set_size_inches( 10, 8 )
  im = ax.imshow( dens_slice )
  plt.colorbar(im)
  plt.savefig( figuresDirectory + 'dens_cholla_{0}.png'.format(n_snapshot))
# 

