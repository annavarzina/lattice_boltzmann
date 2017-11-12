"""
        Module 'd2q5LB' contains a class D2Q5, which represents the logic of
        lattice-boltzmann method on D2Q5 lattice.
"""

import grid as gd
import numpy as np
import matplotlib.pyplot as mat_plot

from timeit import default_timer as timer

#---- Global variables
# TODO divide LB variables from solution variables
Q = 5  # directions
D = 2  # dimensions
conc = 1.  # default density
flux = gd.LatticeVelocity(0, 0) # flux
diff = 1. # default diffusion value
tau = 1.  # default tau
w = []  # weights
e = gd.LatticeVelocity(0, 0)  # lattice velocities
c = 1.  # propagation speed
c_s2 = 1. / 3.  # squared speed of sound
force = gd.LatticeVelocity(0, 0)  # external force


def initialize_weights():
    return np.array([1. / 3.,
                     1. / 6.,
                     1. / 6.,
                     1. / 6.,
                     1. / 6.])


def initialize_normal_velocities():
    e_x = np.array([0, 1, 0, -1, 0])
    e_y = np.array([0, 0, 1, 0, -1])
    return gd.LatticeVelocity(e_x, e_y)


class D2Q5:
    def __init__(self,
                 grid=None,
                 t_final=1000,
                 time_step=1,
                 relaxation_time=None,
                 concentration=None,
                 diffusion = None,
                 velocity = None,
                 lat_speed=None,
                 ext_force=None,
                 boundary=None,
                 solid_cells=None):
        global tau, conc, flux, c, c_s2, force, w, e, Q, D, diff
        w = initialize_weights()
        e = initialize_normal_velocities()
        self.grid = (grid, gd.Grid(10, 10))[grid is None]
        self.solid = (solid_cells, np.zeros((Q, self.grid.width, self.grid.height) ))[solid_cells is None]
        self.tau = (relaxation_time, tau)[relaxation_time is None]
        self.c = (lat_speed, self.grid.dx / time_step)[lat_speed is None]
        self.source = (ext_force, gd.LatticeVelocity(0, 0))[ext_force is None]
        self.iterations = int(t_final / time_step)
        self.f = np.zeros((Q, self.grid.height, self.grid.width))
        self.diff = (diffusion, diff)[diffusion is None]
        # solution initialization
        self.conc = (concentration, conc)[concentration is None]
        print(self.conc)
        self.u = (velocity, gd.LatticeVelocity(np.zeros((self.grid.height, self.grid.width)),
                                               np.zeros((self.grid.height, self.grid.width))))[velocity is None]
        self.flux = gd.LatticeVelocity(np.zeros((self.grid.height, self.grid.width)),
                                       np.zeros((self.grid.height, self.grid.width)))
        self.bc = (boundary, gd.default_poiseuille_boundaries(self.grid))[boundary is None]
        # TODO function run()
        self.run()

    def run(self):
        start = timer()
        for i in range(0, self.iterations):
            self.collision()
            self.streaming()
            self.boundary_conditions()
            self.macroscopic()
            if (self.source.x != 0) | (self.source.y != 0):
                self.source()
            np.set_printoptions(precision=3)
            print('Iteration %s' % (i + 1))
            print('Concentration:')
            print(self.conc)
            print('Flux X:')
            print(self.flux.x)
            print('Flux Y:')
            print(self.flux.y)
            print('Liquid mass is %s [mol]' %np.sum(self.conc))
        end = timer()
        print(str(end - start) + " sec.")

    def macroscopic(self):
        """
        Calculates concentration and flux
        """
        omega = 1./2./self.tau
        self.conc = self.f.sum(axis=0)
        self.flux.x = omega * self.conc * self.u.x + (1- omega) * (self.f[1, :, :] - self.f[3, :, :])
        self.flux.y = omega * self.conc * self.u.y + (1- omega) * (self.f[2, :, :] - self.f[4, :, :])
        # todo add obstacles!

    def source(self):
        """
        Adds an external force
        """
        omega = 1. / tau
        self.flux.x += (self.source.x * self.conc) / (self.conc ** 2 * omega)
        self.flux.y += (self.source.y * self.conc) / (self.conc ** 2 * omega)

    def collision(self):
        """
        Collision step
        """
        f_eq = self.compute_f_equilibrium()
        omega = 1. / tau # omega # TODO rename a to omega
        self.f = (1. - omega) * self.f + omega * f_eq

    def streaming(self):
        """
        Streaming step of LB
        """
        for i in range(Q):
            self.f[i] = np.roll(self.f[i], e.x[i], axis=1)  # roll along x
            self.f[i] = np.roll(self.f[i], e.y[i], axis=0)  # roll along y

    def boundary_conditions(self):
        bc_bounceback = [b for b in self.bc if b.type is 'bounceback']
        # TODO Zoe He flux boundary
        # bc_zoe_he_velocity = [b for b in self.bc if b.type is 'zoe_he_flux']
        bc_zoe_he_conc = [b for b in self.bc if b.type is 'zoe_he_conc']
        bc_neumann = [b for b in self.bc if b.type is 'neumann']

        if len(bc_bounceback) != 0:
            temp_f = np.copy(self.f[:, :, :])
            rows, columns = self.define_row_columns(bc_bounceback, 0)
            # rows = [bt.face / 4 / self.grid.width for bt in bc_bounceback]
            # columns = [(bt.face/ 4) % self.grid.width for bt in bc_bounceback]
            # i = [0, 3, 4, 1, 2]  # bounceback order
            self.f[1, rows, columns] = temp_f[3, rows, columns]
            self.f[2, rows, columns] = temp_f[4, rows, columns]
            self.f[3, rows, columns] = temp_f[1, rows, columns]
            self.f[4, rows, columns] = temp_f[2, rows, columns]

        if len(bc_zoe_he_conc) != 0:
            bc_bottom, bc_top, bc_left, bc_right = self.define_boundaries(bc_zoe_he_conc)
            if len(bc_top) != 0:
                rows, columns = self.define_row_columns(bc_top, 1)
                conc_bc = [bt.value for bt in bc_top]
                conc_bc = [x * 2. / 3. for x in conc_bc]
                self.f[4, rows, columns] = conc_bc - self.f[2, rows, columns]
            if len(bc_bottom) != 0:
                rows, columns = self.define_row_columns(bc_bottom, 0)
                conc_bc = [bt.value for bt in bc_bottom]
                conc_bc = [x * 2. / 3. for x in conc_bc]
                self.f[2, rows, columns] = conc_bc - self.f[4, rows, columns]
            if len(bc_right) != 0:
                rows, columns = self.define_row_columns(bc_right, 3)
                conc_bc = [bt.value for bt in bc_right]
                conc_bc = [x * 2. / 3. for x in conc_bc]
                self.f[3, rows, columns] = conc_bc - self.f[1, rows, columns]
            if len(bc_left) != 0:
                rows, columns = self.define_row_columns(bc_left, 2)
                conc_bc = [bt.value for bt in bc_left]
                conc_bc = [x * 2. / 3. for x in conc_bc]
                self.f[1, rows, columns] = conc_bc - self.f[3, rows, columns]

        if len(bc_neumann) != 0:
            bc_bottom, bc_top, bc_left, bc_right = self.define_boundaries(bc_neumann)
            # todo check the total sum of velocities on boundaries, check is there are pressure boundaries
            if len(bc_top) != 0:
                rows, columns = self.define_row_columns(bc_top, 1)
                u0 = [bt.value for bt in bc_top]
                rho0 = ((self.f[0, rows, columns] + self.f[1, rows, columns] +
                         self.f[3, rows, columns]) + 2 *
                        self.f[2, rows, columns]) / (1. + np.array(u0))
                ru = rho0 * np.array(u0)
                self.f[4, rows, columns] = self.f[2, rows, columns] - (2. / 3.) * ru
            if len(bc_bottom) != 0:
                rows, columns = self.define_row_columns(bc_bottom, 0)
                u0 = [bb.value for bb in bc_bottom]
                rho0 = ((self.f[0, rows, columns] + self.f[1, rows, columns] +
                         self.f[3, rows, columns]) + 2 *
                        self.f[4, rows, columns]) / (1. - np.array(u0))
                ru = rho0 * np.array(u0)
                self.f[2, rows, columns] = self.f[4, rows, columns] + (2. / 3.) * ru
            if len(bc_right) != 0:
                rows, columns = self.define_row_columns(bc_right, 3)
                u0 = [br.value for br in bc_right]
                rho0 = ((self.f[0, rows, columns] + self.f[2, rows, columns] +
                         self.f[4, rows, columns]) + 2 *
                        self.f[1, rows, columns]) / (1. + np.array(u0))
                ru = rho0 * np.array(u0)
                self.f[3, rows, columns] = self.f[1, rows, columns] - (2. / 3.) * ru
            if len(bc_left) != 0:
                rows, columns = self.define_row_columns(bc_left, 2)
                u0 = [bl.value for bl in bc_left]
                rho0 = ((self.f[0, rows, columns] + self.f[2, rows, columns] +
                         self.f[4, rows, columns]) + 2 *
                        self.f[3, rows, columns]) / (1. - np.array(u0))
                ru = rho0 * np.array(u0)
                self.f[1, rows, columns] = self.f[3, rows, columns] + (2. / 3.) * ru

    @staticmethod
    def define_boundaries(bc):
        bc_bottom = [bn for bn in bc if bn.face % 4 == 0]
        bc_top = [bn for bn in bc if bn.face % 4 == 1]
        bc_left = [bn for bn in bc if bn.face % 4 == 2]
        bc_right = [bn for bn in bc if bn.face % 4 == 3]
        return bc_bottom, bc_top, bc_left, bc_right

    def define_row_columns(self, bc, d):
        rows = [(bt.face - d) / 4 / self.grid.width for bt in bc]
        columns = [((bt.face - d) / 4) % self.grid.width for bt in bc]
        return rows, columns

    def compute_f_equilibrium(self):
        """
        Equilibrium distribution function:
                                 e u
        f_eq = weight*rho*(1 + 3 --- )
                                 c^2
        :return: equilibrium distribution function
        """
        f_eq = np.zeros((Q, self.grid.height, self.grid.width))
        f_eq[0,None, None] = w[0, None, None] * self.conc
        f_eq[1,None, None] = w[1, None, None] * self.conc * (1. + 3. * self.u.x)
        f_eq[2,None, None] = w[2, None, None] * self.conc * (1. + 3. * self.u.y)
        f_eq[3,None, None] = w[3, None, None] * self.conc * (1. - 3. * self.u.x)
        f_eq[4,None, None] = w[4, None, None] * self.conc * (1. - 3. * self.u.y)
        return f_eq


def plotting(solution):
    l_x = solution.grid.width
    l_y = solution.grid.height
    x = np.arange(0, l_x, 1)
    y = np.arange(0, l_y, 1)
    x_mat, y_mat = np.meshgrid(x, y)
    u = (force.y / (2 * 1. / 6.)) * (((l_x - 1) / 2.) ** 2 - (x_mat - (l_x - 1) / 2.) ** 2)  # analytical solution
    mat_plot.figure(1)
    mat_plot.plot(solution.u.y[0, :], 'r-', u[0, :], 'ro')

    mat_plot.show()

def concentration_boundaries(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_conc', value = 0.))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_conc', value = 0.))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_conc', value = 0.))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_conc', value = 0.))
    return boundaries


def bb_boundaries(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'bounceback'))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'bounceback'))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'bounceback'))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'bounceback'))
    return boundaries

if __name__ == "__main__":
    n = 3 #5  # width  - lx
    m = 3 #5  # height - ly
    t = 5  # final time
    c_init = np.zeros((n, m))
    c_init[1, 1] = 1
    #c_init[2, 2] = 1
    solid = np.zeros((Q, n, m))
    userGrid = gd.Grid(n, m)
    bc = bb_boundaries(userGrid)
    #bc = concentration_boundaries(userGrid)
    lat_bol = D2Q5(grid=userGrid,
                   t_final=t,
                   concentration=c_init,
                   boundary=bc)
    #print(lat_bol.__dict__)
    '''
    n = 31  # width  - lx
    m = 51  # height - ly
    dx = 1  # cell size
    dt = 1  # time step
    t_f = 3000  # final time
    gravity = gd.LatticeVelocity(0, 1e-8)  # acceleration by gravity
    solid = np.zeros((Q, n, m))

    userGrid = gd.Grid(n, m)
    lat_bol = D2Q5(grid=userGrid, t_final=t_f,
                   ext_force=gravity)
    plotting(lat_bol)
    '''