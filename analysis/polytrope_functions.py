import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
# import matplotlib

G = 6.67259e-8 # gravitational constant, cgs


def binary_search( val, N, data,  indx_l, indx_r ):
  print indx_l, indx_r
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
  # print x
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
    # print coords_new[0]
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
  print( " Root Theta: {0}  ->  {1}  {2}".format( psi_root, theta_vals[root_indx], theta_vals[root_indx+1] ))

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

# 
# fig, ax = plt.subplots()
# ax.set_xlim(np.array([0,10]))
# ax.set_ylim(np.array([0,1]))
# 
# ax.plot(xi, solution,label=r'$n={0:.1f}$'.format(ns))
# ax.legend(loc=0)
# plt.savefig('polytropes_1D.png')
# 
# 
# def interpolate( x, x_all, y_all ):
#   n = len(x_all)
#   # print n
#   # for i in range(n):
#   #   if ( x < x_all[i] ): break
#   #   indx = i
#   if x < x_all[0]: return y_all[0]
#   if x > x_all[n-1]: return y_all[n-1]
#   indx = binary_search( x, n, x_all, 0, n-1)
#   if x < x_all[indx]: print "interpolation error"
#   if x > x_all[indx+1]: print "interpolation error"
#   # if indx == n-1:
#     # if x > x_all[indx]: print "radius outside interpolation range"
#   # else:
#   #   if x > x_all[indx+1]: print "radius outside interpolation range"
#   x_0 = x_all[indx]
#   x_1 = x_all[indx + 1]
#   dx = x_1 - x_0
#   y_0 = y_all[indx]
#   y_1 = y_all[indx + 1]
#   if x<x_0:print( 'Error in the interpolation')
#   y = y_0 + ( y_1 - y_0 )/dx*(x-x_0)
#   return y
# 
# 
# n_cells = 256
# L = 20.0
# dx = L/n_cells
# dy = L/n_cells
# dens_all = np.zeros([ n_cells, n_cells ])
# x_min = -10.0
# x_max = x_min + L
# y_min = -10.0
# y_max = y_min + L
# 
# for j in range(n_cells):
#   for i in range(n_cells):
#     x = x_min + (i + 0.5)*dx
#     y = y_min + (j + 0.5)*dy
#     r = np.sqrt( x*x + y*y )
#     dens = interpolate( r, xi, solution ) 
#     dens_all[j,i] = dens
# 
# fig, ax = plt.subplots()
# ax.imshow( dens_all )
# plt.savefig('polytropes_2D.png')

# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# def func(coords, xi, n):
#   theta, y = coords
#   dydt = [y, - (theta ** n) - 2. * y / xi]
#   return dydt
# 
# 
# xi = np.linspace(1.e-11, 10., 1000)
# ns = np.array([ 3. ])
# 
# lims = np.array([2.4494, 3.1415, 3.65375, 6.89685, xi[-1]])
# labels = [r'$n=0$',r'$n=1$',r'$n=3/2$',r'$n=3$',r'$n=5$']
# 
# sol = np.empty((len(ns), len(xi)))
# fig, ax = plt.subplots()
# ax.set_xlim(np.array([0,xi[-1]]))
# ax.set_ylim(np.array([0,1]))
# 
# # for i in range(0,len(ns)):
# ns = 3
# lim = 6.8968
# coords0 = np.array([1., 0.])
# sol = odeint(func, coords0, xi, args = (ns,))[:,0]
# solgt0 = sol[xi<lim]
# ax.plot(xi[xi<lim], solgt0,label=r'$n={0:.1f}$'.format(ns))
# ax.legend(loc=0)
# plt.savefig('polytropes.png')