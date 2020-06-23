import sys, os
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import h5py as h5

stars_dir = os.path.dirname(os.getcwd()) + '/'
subDirectories = [x[0] for x in os.walk(stars_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_cholla, load_snapshot_cholla_distributed
from tools import create_directory


dataDir = '/data/groups/comp-astro/bruno/'
inputDir  = dataDir + 'polytrope/output_files/'
outputDir = dataDir + 'polytrope/figures/'
create_directory( outputDir )



n_snapshot = 0
print('nSnapshot: {0}'.format( n_snapshot ) )

data_single = load_snapshot_cholla( n_snapshot, inputDir, single_file=True )
density_single = data_single['gas']['density'][...]
potential_single = data_single['gas']['potential'][...]


data_type = 'hydro'
fields = ['density', 'potential' ]
data_snapshot = load_snapshot_cholla_distributed( n_snapshot, data_type, fields, inputDir, subgrid = None, show_progess = True, precision = np.float64)

dx = data_snapshot['header']['dx']
density_mpi = data_snapshot[data_type]['density']
potential_mpi = data_snapshot[data_type]['potential']

diff_pot = ( potential_mpi - potential_single ) / potential_single
diff_dens = ( density_mpi - density_single ) / density_single


print(diff_pot.min(), diff_pot.max())
print( diff_dens.min(), diff_dens.max())