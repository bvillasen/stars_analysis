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
from polytrope_functions import Solve_Polytropic_Star
from load_data_cholla import load_snapshot_cholla
from tools import create_directory

dataDir = '/data/groups/comp-astro/bruno/stars/polytrope/'
inputDir = dataDir + 'output_files/grav_0./'
outputDir = figuresDirectory + 'grav_0./'
create_directory( outputDir )

n_poly = 1.5
psi_max = 5
n_points = 100000

M_star = 1.98847e33 #Msun in gr;  //Solar Mass
R_star = 6.957e+10  #Solar radius in cm
dens_min = 1e-6

R_vals, density_vals = Solve_Polytropic_Star( M_star, R_star, n_poly, psi_max, n_points, dens_min=1e-6)

# Get the polytropic density profile
r_poly = np.concatenate([ -R_vals[::-1], R_vals ]) / R_star
dens_poly = np.concatenate([ density_vals[::-1], density_vals ])

#Plot the density profile

 
 
 

n_snapshots = 11
n_snapshot = 0
for n_snapshot in range( n_snapshots ):
  print 'nSnapshot: {0}'.format( nSnapshot )
  data = load_snapshot_cholla( n_snapshot, inputDir, single_file=True )
  t = data['t']
  gamma = data['gamma']

  dens = data['gas']['density'][...]
  vx = data['gas']['momentum_x'][...] / dens
  vy = data['gas']['momentum_y'][...] / dens
  vz = data['gas']['momentum_z'][...] / dens
  E = data['gas']['Energy'][...]
  v2 = vx*vx + vy*vy + vz*vz 
  U = E - 0.5*dens*v2
  p = (gamma-1) * U
  nz, ny, nx = dens.shape
  Lbox = 20e+10
  x_cholla =  ( ( np.linspace( 0, nx-1, nx) + 0.5 ) *( Lbox / nx ) - Lbox/2 ) / R_star
  indx_slice = nz/2
  dens_slice = dens[indx_slice, :, : ]
  dens_profile = dens[nz/2, ny/2, :]

  pressure_profile = p[nz/2, ny/2, :]


  fig, ax = plt.subplots()
  fig.set_size_inches( 10, 8 )
  ax.plot( r_poly, dens_poly, label='n=1.5' )
  ax.plot( x_cholla, dens_profile, label='Cholla' )
  ax.legend()
  text = 't={0:.0f}  s'.format(t)
  plt.text(0.05, 0.95, text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, fontsize=15)
  fig.savefig( outputDir + 'polytrope_1D_{0}.png'.format(n_snapshot)) 

  fig, ax = plt.subplots()
  fig.set_size_inches( 10, 8 )
  # ax.plot( r_poly, dens_poly, label='n=1.5' )
  ax.plot( x_cholla, pressure_profile, label='Cholla' )
  ax.legend()
  text = 't={0:.0f}  s'.format(t)
  plt.text(0.05, 0.95, text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, fontsize=15)
  fig.savefig( outputDir + 'polytrope_pressure_1D_{0}.png'.format(n_snapshot)) 

    # fig, ax = plt.subplots()
    # fig.set_size_inches( 10, 8 )
    # im = ax.imshow( dens_slice )
    # plt.colorbar(im)
    # plt.savefig( outputDir + 'polytrope_2D_{0}.png'.format(n_snapshot))
  # # 

