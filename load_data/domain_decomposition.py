

def get_domain_block( proc_grid, box_size, grid_size ):
  np_x, np_y, np_z = proc_grid
  Lx, Ly, Lz = box_size
  nx_g, ny_g, nz_g = grid_size
  dx, dy, dz = Lx/np_x, Ly/np_y, Lz/np_z
  nx_l, ny_l, nz_l = nx_g/np_x, ny_g/np_z, nz_g/np_z,

  nprocs = np_x * np_y * np_z
  domain = {}
  domain['global'] = {}
  domain['global']['dx'] = dx
  domain['global']['dy'] = dy
  domain['global']['dz'] = dz
  for k in range(np_z):
    for j in range(np_y):
      for i in range(np_x):
        pId = i + j*np_x + k*np_x*np_y
        domain[pId] = { 'box':{}, 'grid':{} }
        xMin, xMax = i*dx, (i+1)*dx
        yMin, yMax = j*dy, (j+1)*dy
        zMin, zMax = k*dz, (k+1)*dz
        domain[pId]['box']['x'] = [xMin, xMax]
        domain[pId]['box']['y'] = [yMin, yMax]
        domain[pId]['box']['z'] = [zMin, zMax]
        domain[pId]['box']['dx'] = dx
        domain[pId]['box']['dy'] = dy
        domain[pId]['box']['dz'] = dz
        domain[pId]['box']['center_x'] = ( xMin + xMax )/2.
        domain[pId]['box']['center_y'] = ( yMin + yMax )/2.
        domain[pId]['box']['center_z'] = ( zMin + zMax )/2.
        gxMin, gxMax = i*nx_l, (i+1)*nx_l
        gyMin, gyMax = j*ny_l, (j+1)*ny_l
        gzMin, gzMax = k*ny_l, (k+1)*ny_l
        domain[pId]['grid']['x'] = [gxMin, gxMax]
        domain[pId]['grid']['y'] = [gyMin, gyMax]
        domain[pId]['grid']['z'] = [gzMin, gzMax]
  return domain

def get_domain_parent( domain, domain_parent):
  for i in domain.keys():
    proc_domain = domain[i]
    center_x = proc_domain['box']['center_x']
    center_y = proc_domain['box']['center_y']
    center_z = proc_domain['box']['center_z']
    for parent_id in domain_parent.keys():
      xMin, xMax = domain_parent[parent_id]['box']['x']
      yMin, yMax = domain_parent[parent_id]['box']['y']
      zMin, zMax = domain_parent[parent_id]['box']['z']
      if ( center_x < xMin ) or ( center_x > xMax ): continue
      if ( center_y < yMin ) or ( center_y > yMax ): continue
      if ( center_z < zMin ) or ( center_z > zMax ): continue
      proc_domain['parent_id'] = parent_id
