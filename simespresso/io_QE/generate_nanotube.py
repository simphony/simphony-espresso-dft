__author__ = 'jeremy'
import scipy.optimize
import numpy as np

def generate_polyhedral_nanotube(m,n):
    initial_value = 0
#    scipy.optimize.minimize(fun, x0, args=(), method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None)
    x = scipy.optimize.newton(test_eqn,initial_value,args=(m,n))
    print('x {} f(x) {}'.format(x,test_eqn(x,m,n)))
    x = scipy.optimize.newton(phi_transcedental,initial_value,args=(m,n))
    print('x {} f(x) {}'.format(x,phi_transcedental(x,m,n)))

def test_eqn(x,m,n):
    print('x {} m {} n {}'.format(x,m,n))
    y = m*x**2+n*x+1
    return y

def phi_transcedental(phi,m,n):
    print('phi {} m {} n {}'.format(phi,m,n))
    ksi = (n*phi - np.pi)/m
    x = (n**2 - m**2)*(np.sin(ksi+phi)**2)-n*(n+2*m)*(np.sin(ksi))**2+m*(2*n+m)*(np.sin(phi))**2
    return x

def generate_nanotube_xyz(type = 'armchair'):
    #coords from http://turin.nss.udel.edu/research/tubegenonline.html
    if type == 'armchair':
    #armchair
        xyz_33 = [[2.064591,     0.000000,    -1.232141],
           [1.575573,     1.334205,    -1.232141],
           [1.032295,     1.787988,   -0.000000],
           [-0.367669,     2.031589,    -0.000000],
           [-1.032295,     1.787988,    -1.232141],
           [-1.943242,     0.697384,    -1.232141],
           [-2.064591,     0.000000,    -0.000000],
           [-1.575573,    -1.334205,    -0.000000],
           [-1.032295,    -1.787988,    -1.232141],
           [0.367669,   -2.031589,    -1.232141],
           [1.032295,    -1.787988,    -0.000000],
           [1.943242,    -0.697384,    -0.000000]]
        return xyz_33
    elif type == 'zigzag':
   #n = 3, m = 0  zigzag
        xyz_30= [[1.216023,     0.000000,    -0.735207],
          [1.216023,     0.000000,    -2.156207],
          [0.608011,     1.053107,     1.421000],
          [0.608011,     1.053107,     0.000000],
          [-0.608011,     1.053107,    -0.735207],
          [-0.608011,     1.053107,    -2.156207],
          [-1.216023,    -0.000000,     1.421000],
          [-1.216023,     0.000000,     0.000000],
          [-0.608011,    -1.053107,    -0.735207],
          [-0.608011,    -1.053107,    -2.156207],
          [0.608011,    -1.053107,     1.421000],
          [0.608011,    -1.053107,    -0.000000]]
        return xyz_30
    else:
           #n = 0, m = 3
        xyz_03 = [[1.216023,     0.000000,    -2.156207],
          [0.608011,     1.053107,    -1.421000],
          [-0.608011,     1.053107,    -2.156207],
          [-1.216023,    -0.000000,    -1.421000],
          [-0.608011,    -1.053107,    -2.156207],
          [ 0.608011,    -1.053107,    -1.421000],
          [ 1.216023,     0.000000,     0.735207],
          [ 0.608011,     1.053107,     0.000000],
          [-0.608011,     1.053107,     0.735207],
          [-1.216023,     0.000000,     0.000000],
          [-0.608011,    -1.053107,     0.735207],
          [ 0.608011,    -1.053107,    -0.000000]]
        return xyz_03


