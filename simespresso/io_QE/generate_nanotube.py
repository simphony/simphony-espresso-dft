import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt


def generate_polyhedral_nanotube(m, n, sigma=1.44, n_atoms_per_helix=10):
    initial_value = np.pi*(2*n+m)/(2*(n**2+n*m+m**2))
    if m == 0:
        initial_value = 0.1
#    scipy.optimize.minimize(
# fun, x0, args=(), method=None, jac=None, hess=None, hessp=None, bounds=None,
# constraints=(), tol=None, callback=None, options=None)
#    x = scipy.optimize.newton(test_eqn,initial_value,args=(m,n))
#    print('x {} f(x) {}'.format(x,test_eqn(x,m,n)))
    R0 = conventional_r0(sigma, m, n)
#    theta1 = conventional_chiral_angle(m,n)
    theta2 = conventional_chiral_angle2(m, n)
    phi = scipy.optimize.newton(phi_newton, initial_value, args=(m, n))
    print('theta {} r0 {} phi {}'.format(theta2, R0, phi))
#    print('phi {} f(phi) {} ok? {}'.format(
# phi,phi_newton(phi,m,n),check_bounds(m,n,phi)))
#    theta = theta_direct(phi,m,n)
#    print('thetadirect {} f(theta) {}'.format(theta,theta_direct(phi,m,n)))
#   theta = scipy.optimize.newton(theta_newton,theta,args=(phi,m,n))
#    print('thetanewton {} f(theta) {}'.format(
# theta,theta_newton(theta,phi,m,n)))
    n_helices = np.pi*2/phi
    print('n_helices '+str(n_helices))
    n_helices = int(np.round(n_helices))
#    n_helices = 2
    print('n_helices '+str(n_helices))
    current_phi = 0
    positions = []
    for i in range(0, n_helices):
        z0 = 0
        positions = positions + generate_helix(
            theta2, current_phi, phi, sigma, R0, n_atoms_per_helix, z0)
        current_phi = current_phi + np.pi*2 / (n_helices)

    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]
    # print('xs:'+str(xs))
    filewrite(positions)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs)
    plt.show()


def filewrite(positions):
    with open('nanotube.txt', 'w') as f:
        for p in positions:
            f.write(str(p[0])+','+str(p[1])+','+str(p[2])+'\n')


def generate_helix(theta, initial_phi, phi, sigma, R0, n_atoms, z0):
    print('new helix, theta {} phi0 {} phi {} sigma {} R0 {} n {} z0'.format(
        theta, initial_phi, phi, sigma, R0, n_atoms, z0))
    positions = []
    current_phi = initial_phi
    current_z = z0
    for i in range(n_atoms):
        current_x = R0*np.cos(current_phi)
        current_y = R0*np.sin(current_phi)
        positions.append([current_x, current_y, current_z])
        current_phi = (current_phi + phi) % ((np.pi)*2)
        current_z = current_z + sigma*np.sin(theta)
        # print('x {} y {} z {} phi {}'.format(current_x,
        # current_y,current_z,current_phi))
    return positions


def test_eqn(x, m, n):
    print('x {} m {} n {}'.format(x, m, n))
    y = m*x**2+n*x+1
    return y


def conventional_chiral_angle(m, n):
    print('chiral angle1: m {} n {}'.format(m, n))
    costheta = (2*n+m)/(np.sqrt(n**2+n*m+m**2))
    print('costheta {}'.format(costheta))
    theta = np.arccos(costheta)
#    print('theta conv:'+str(theta))
    return theta


def conventional_chiral_angle2(m, n):
    print('chiral angle2: m {} n {}'.format(m, n))
    tantheta = m*np.sqrt(3)/(2*n+m)
    print('tantheta {}'.format(tantheta))
    theta = np.arctan(tantheta)
#    print('theta conv:'+str(theta))
    return theta


def conventional_r0(sigma, m, n):
    print('sigma {} m {} n {}'.format(sigma, m, n))
    r0 = sigma*np.sqrt(3*(n**2 + n*m + m**2))/(2*np.pi)
    return r0


def check_bounds(m, n, phi):
    ok = np.pi / (n+m) < phi and phi < np.pi/n
    return ok


def theta_direct(phi, m, n):
    print('calculating theta directly: phi {} m {} n {}'.format(phi, m, n))
    if m == 0:
        ksi = 0
    else:
        ksi = (n*phi - np.pi)/m
    costheta2 = (n*(n+2*m) * (np.sin(phi))**2) / \
                ((n+m)**2 * (np.sin(phi))**2 - m**2*(np.sin(ksi+phi))**2)
    theta = np.arccos(np.sqrt(costheta2))
    print('theta calculated '+str(theta))
    return theta


def theta_newton(theta, phi, m, n):
    print('calculating theta using newton: phi {} m {} n {}'.format(
        phi, m, n))
    if m == 0:
        ksi = 0
    else:
        ksi = (n*phi - np.pi)/m
    costheta2 = (n*(n+2*m) * (np.sin(phi))**2)/((n+m)**2*(np.sin(phi))**2 -
                                                m**2*(np.sin(ksi+phi))**2)
    x = np.cos(theta)**2 - costheta2
    print('THETA RESIDUAL '+str(x))
    return x



def phi_newton(phi, m, n):
    print('phi {} m {} n {}'.format(phi, m, n))
    if m == 0:
        ksi = 0
    else:
        ksi = (n*phi - np.pi)/m
    x = (n**2 - m**2) * (np.sin(ksi+phi)**2) - n * (n+2*m) * \
        (np.sin(ksi))**2 + m*(2*n+m)*(np.sin(phi))**2
    print('phi residual '+str(x))
    return x


def read_xyz(filename):
    with open(filename, 'r') as f:
        atomsline = f.readline()
        n_atoms = int(atomsline)
        print('atomsline:'+str(atomsline))
        commentline = f.readline()
        print('commentline:'+str(commentline))
        print('n_atoms:'+str(n_atoms))
        positions_list = []
        for n in range(n_atoms):
            line = f.readline()
            print(line)
            vals = line.split()
            atomtype = vals[0]
            pos = [float(vals[1]), float(vals[2]), float(vals[3])]
            print('atom:' + str(atomtype) + 'pos:' + str(pos))
            positions_list.append(pos)
        return positions_list


def show_xyz(filename):
    positions = read_xyz(filename)
    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs)
    plt.show()


def generate_nanotube_xyz(kind='armchair'):
    # coords from http://turin.nss.udel.edu/research/tubegenonline.html
    if kind == 'armchair':  # armchair
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
    elif kind == 'zigzag':
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
                  [0.608011,    -1.053107,    -1.421000],
                  [1.216023,     0.000000,     0.735207],
                  [0.608011,     1.053107,     0.000000],
                  [-0.608011,     1.053107,     0.735207],
                  [-1.216023,     0.000000,     0.000000],
                  [-0.608011,    -1.053107,     0.735207],
                  [0.608011,    -1.053107,    -0.000000]]
        return xyz_03


if __name__ == "__main__":
    generate_polyhedral_nanotube(0,3)
