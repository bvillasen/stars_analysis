import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
# import matplotlib

G = 6.67259e-8 # gravitational constant, cgs


def binary_search( val, N, data,  indx_l, indx_r ):
  n = indx_r - indx_l
  indx = indx_l+ n/2
  if val >= data[N-1]: return indx_r
  if val <= data[0]:  return indx_l
  if indx_r == indx_l+1: return indx_l
  if data[indx] <= val: 
    indx_r = indx_l + n
    indx_l = indx
  else:
    indx_l = indx_l
    indx_r = indx
  return binary_search( val, N, data, indx_l, indx_r)

  
  



def RK4_step( f, t, x, dt, **kargs ):
  k1 = dt * f( t, x, kargs=kargs )
  k2 = dt * f( t + 0.5*dt, x+0.5*k1, kargs=kargs )
  k3 = dt * f( t + 0.5*dt, x+0.5*k2, kargs=kargs )
  k4 = dt * f( t + dt, x+k3, kargs=kargs )
  return x + 1./6*( k1 + 2*k2 + 2*k3 + k4)


def f( x, coords, kargs=None):
  n = kargs['n']
  theta, y = coords
  if theta < 0: theta = 0
  dydt = np.array([y, - (theta ** n) - 2. * y / x])
  return dydt


def Solve_Polytropic( ns, psi_max, n_points ):

  xi = np.linspace(1.e-11, psi_max, n_points)
  y = np.zeros([n_points, 2])
  y[0] = np.array([ 1, 0 ])

  for i in range(n_points-1):
    dx = xi[i+1] - xi[i]
    coords = y[i]
    x = xi[i]
    coords_new = RK4_step( f, x, coords, dx, n=ns )
    # if coords_new[0] < 0: coords_new[0] = 0
    y[i+1] = coords_new
  theta = y[:,0]
  theta_deriv = y[:,1]
  return xi, theta, theta_deriv
  
  
  
def Solve_Polytropic_Star( M_star, R_star, n_poly, psi_max, n_points, dens_min=1e-6):

  psi_vals, theta_vals, theta_deriv = Solve_Polytropic( n_poly, psi_max, n_points)

  for i in range(n_points): 
    if theta_vals[i] * theta_vals[i+1] < 0: 
      root_indx = i
      break
  psi_root = psi_vals[root_indx] 
  theta_deriv_root = theta_deriv[root_indx]

  # Convert tho physical values
  dens_avrg = ( 3 * M_star ) / ( 4 * np.pi * R_star**3 )
  beta = - ( psi_root / 3 / theta_deriv_root )
  dens_central = beta * dens_avrg
  pressure_central = G * M_star * M_star / ( R_star** 4 ) / ( 4 * np.pi *( n_poly+1 ) * theta_deriv_root * theta_deriv_root )
  K = pressure_central * dens_central**( -(n_poly+1)/n_poly  )
  alpha = np.sqrt( (n_poly + 1) * K / ( 4 * np.pi * G ) ) *  dens_central**( (1-n_poly)/(2*n_poly) )

  theta_min = (dens_min/ dens_central)**(1./n_poly) 
  theta_vals[theta_vals<dens_min] = theta_min
  R_vals = alpha * psi_vals;
  density_vals = dens_central * theta_vals**n_poly ;
  return R_vals, density_vals
