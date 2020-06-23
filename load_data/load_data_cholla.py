import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np

from domain_decomposition import get_domain_block




def load_snapshot_cholla_distributed( n_snapshot, data_type, fields, inputDir, subgrid = None, show_progess = True, precision = np.float64):

  base_name = '.h5.'
  fileName = '{0}{1}0'.format( n_snapshot, base_name )

  #Load The first file to get parameters
  fileName = inputDir + fileName

  file = h5.File( fileName, 'r' )
  dims   = file.attrs['dims']
  nprocs = file.attrs['nprocs']
  domain = file.attrs['domain']

  header = {}
  h_keys = file.attrs.keys()
  for key in h_keys:
    data_key = file.attrs[key]
    if data_key.shape == (1,): header[key] = data_key[0]
    else: header[key] = data_key 
  file.close()

  grid_size = dims
  proc_grid = nprocs
  box_size  = domain
  domain = get_domain_block( proc_grid, box_size, grid_size )
  if subgrid == None: subgrid = [ [0,grid_size[0]], [0,grid_size[1]], [0,grid_size[2]], ]
  data_snapshot = load_snapshot_data_distributed( n_snapshot, inputDir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  data_snapshot['header'] = header
  return data_snapshot



def select_procid( proc_id, subgrid, domain, ids, ax ):
  domain_l, domain_r = domain
  subgrid_l, subgrid_r = subgrid
  if domain_l <= subgrid_l and domain_r > subgrid_l:
    ids.append(proc_id)
  if domain_l >= subgrid_l and domain_r <= subgrid_r:
    ids.append(proc_id)
  if domain_l < subgrid_r and domain_r >= subgrid_r:
    ids.append(proc_id)


def select_ids_to_load( subgrid, domain, proc_grid ):
  subgrid_x, subgrid_y, subgrid_z = subgrid
  nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2]
  ids_x, ids_y, ids_z = [], [], []
  for proc_id in range(nprocs):
    domain_local = domain[proc_id]
    domain_x = domain_local['grid']['x']
    domain_y = domain_local['grid']['y']
    domain_z = domain_local['grid']['z']
    select_procid( proc_id, subgrid_x, domain_x, ids_x, 'x' )
    select_procid( proc_id, subgrid_y, domain_y, ids_y, 'y' )
    select_procid( proc_id, subgrid_z, domain_z, ids_z, 'z' )
  set_x = set(ids_x)
  set_y = set(ids_y)
  set_z = set(ids_z)
  set_ids = (set_x.intersection(set_y)).intersection(set_z )
  return list(set_ids)
def load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid, domain, precision, proc_grid, show_progess=True, kernel_types=[] ):
  # Find the ids to load 
  ids_to_load = select_ids_to_load( subgrid, domain, proc_grid )

  #Find the boundaries of the volume to load
  domains = { 'x':{'l':[], 'r':[]}, 'y':{'l':[], 'r':[]}, 'z':{'l':[], 'r':[]}, }
  for id in ids_to_load:
    for ax in domains.keys():
      d_l, d_r = domain[id]['grid'][ax]
      domains[ax]['l'].append(d_l)
      domains[ax]['r'].append(d_r)
  boundaries = {}
  for ax in domains.keys():
    boundaries[ax] = [ min(domains[ax]['l']),  max(domains[ax]['r']) ]

  # Get the size of the volume to load
  nx = boundaries['x'][1] - boundaries['x'][0]    
  ny = boundaries['y'][1] - boundaries['y'][0]    
  nz = boundaries['z'][1] - boundaries['z'][0]    

  dims_all = [ nx, ny, nz ]

  data_out = {}
  data_out[data_type] = {}
  
  if type(fields) != list: fields = [fields]
  print( "Loading Snapshot: {0}  ->  {1}".format(nSnap, fields))

  if data_type == 'sph_kernel':
    
    added_header = True
    n_to_load = len(ids_to_load)
    
    for kernel_type in kernel_types:
      
      data_out[data_type][kernel_type] = {}
      
      for field in fields:
        data_all = np.zeros( dims_all, dtype=precision )
        for i, nBox in enumerate(ids_to_load):
          name_base = 'h5'
          inFileName = '{0}.{1}.{2}'.format(nSnap, name_base, nBox)
          
          inFile = h5.File( inDir + inFileName, 'r')
          head = inFile.attrs
          if added_header == False:
            for h_key in head.keys():
              if h_key in ['dims', 'dims_local', 'offset', 'bounds', 'domain', 'dx', ]: continue
              data_out[h_key] = head[h_key][0]
              if h_key == 'current_z': print (' current_z: {0}'.format( data_out[h_key]))
            added_header = True
            
          if show_progess:
            terminalString  = '\r Loading File: {0}/{1}   {2}'.format(i, n_to_load, field)
            sys.stdout. write(terminalString)
            sys.stdout.flush() 

          procStart_x, procStart_y, procStart_z = head['offset']
          procEnd_x, procEnd_y, procEnd_z = head['offset'] + head['dims_local']
          # Substract the offsets
          procStart_x -= boundaries['x'][0]
          procEnd_x   -= boundaries['x'][0]
          procStart_y -= boundaries['y'][0]
          procEnd_y   -= boundaries['y'][0]
          procStart_z -= boundaries['z'][0]
          procEnd_z   -= boundaries['z'][0]
          data_local = inFile[kernel_type][field][...]

          data_all[ procStart_x:procEnd_x, procStart_y:procEnd_y, procStart_z:procEnd_z] = data_local

        # Trim off the excess data on the boundaries:
        trim_x_l = subgrid[0][0] - boundaries['x'][0]
        trim_x_r = boundaries['x'][1] - subgrid[0][1]  
        trim_y_l = subgrid[1][0] - boundaries['y'][0]
        trim_y_r = boundaries['y'][1] - subgrid[1][1]  
        trim_z_l = subgrid[2][0] - boundaries['z'][0]
        trim_z_r = boundaries['z'][1] - subgrid[2][1]  
        data_output = data_all[trim_x_l:nx-trim_x_r, trim_y_l:ny-trim_y_r, trim_z_l:nz-trim_z_r,  ]
        data_out[data_type][kernel_type][field] = data_output
        if show_progess: print("")

    
  else:

    added_header = True
    n_to_load = len(ids_to_load)
    
    for field in fields:
      data_all = np.zeros( dims_all, dtype=precision )
      for i, nBox in enumerate(ids_to_load):
        name_base = 'h5'
        if data_type == 'particles': inFileName = '{0}_particles.{1}.{2}'.format(nSnap, name_base, nBox)
        if data_type == 'hydro':     inFileName = '{0}.{1}.{2}'.format(nSnap, name_base, nBox)
        if data_type == 'fft':       inFileName ='{0}_data_fft.{1}.{2}'.format(nSnap, name_base, nBox)
        
        inFile = h5.File( inDir + inFileName, 'r')
        head = inFile.attrs
        if added_header == False:
          for h_key in head.keys():
            if h_key in ['dims', 'dims_local', 'offset', 'bounds', 'domain', 'dx', ]: continue
            data_out[h_key] = head[h_key][0]
            if h_key == 'current_z': print( ' current_z: {0}'.format( data_out[h_key]))
          added_header = True
          
        if show_progess:
          terminalString  = '\r Loading File: {0}/{1}   {2}'.format(i, n_to_load, field)
          sys.stdout. write(terminalString)
          sys.stdout.flush() 

        procStart_x, procStart_y, procStart_z = head['offset']
        procEnd_x, procEnd_y, procEnd_z = head['offset'] + head['dims_local']
        # Substract the offsets
        procStart_x -= boundaries['x'][0]
        procEnd_x   -= boundaries['x'][0]
        procStart_y -= boundaries['y'][0]
        procEnd_y   -= boundaries['y'][0]
        procStart_z -= boundaries['z'][0]
        procEnd_z   -= boundaries['z'][0]
        data_local = inFile[field][...]
        data_all[ procStart_x:procEnd_x, procStart_y:procEnd_y, procStart_z:procEnd_z] = data_local

      # Trim off the excess data on the boundaries:
      trim_x_l = subgrid[0][0] - boundaries['x'][0]
      trim_x_r = boundaries['x'][1] - subgrid[0][1]  
      trim_y_l = subgrid[1][0] - boundaries['y'][0]
      trim_y_r = boundaries['y'][1] - subgrid[1][1]  
      trim_z_l = subgrid[2][0] - boundaries['z'][0]
      trim_z_r = boundaries['z'][1] - subgrid[2][1]  
      data_output = data_all[trim_x_l:nx-trim_x_r, trim_y_l:ny-trim_y_r, trim_z_l:nz-trim_z_r,  ]
      data_out[data_type][field] = data_output
      if show_progess: print("")
    
  return data_out



def load_snapshot_cholla( nSnap, inDir, single_file=False  ):
  gridFileName = inDir + 'grid_{0:03}.h5'.format(nSnap)
  if single_file: gridFileName = inDir + '{0}.h5'.format(nSnap)
  
  outDir = { 'gas':{} }
  
  data_grid = h5.File( gridFileName, 'r' )
  fields_data = data_grid.keys()
  if single_file: 
    for key in data_grid.attrs.keys(): 
      data_key = data_grid.attrs[key]
      if data_key.shape == (1,) :  
        outDir[key] = data_grid.attrs[key][0]
      else: outDir[key] = data_grid.attrs[key][...]
  else: 
    for key in data_grid.attrs.keys(): outDir[key] = data_grid.attrs[key]
  fields_grid = fields_data
  for field in fields_grid:
    if field not in fields_data: continue
    outDir['gas'][field] = data_grid[field]

  return outDir
