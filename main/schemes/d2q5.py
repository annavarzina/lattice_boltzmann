"""
        Module 'd2q5LB' contains a class D2Q5, which represents the logic of
        lattice-boltzmann method on D2Q5 lattice.
"""

import grid as gd
import numpy as np
import matplotlib.pyplot as mat_plot

from timeit import default_timer as timer

Q = 5  # directions
D = 2  # dimensions
rho = 1.  # default density
tau = 1.  # default tau
w = []  # weights
e = gd.LatticeVelocity(0, 0)  # lattice velocities
c = 1.  # propagation speed
c_s2 = 1. / 3.  # squared speed of sound
force = gd.LatticeVelocity(0, 0)  # external force


def initialize_weights():
    return np.array([1. / 3., 1. / 6., 1. / 6., 1. / 6., 1. / 6.])


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
                 density=None,
                 lat_speed=None,
                 ext_force=None,
                 boundary=None,
                 solid_cells=None):
        global tau, rho, c, c_s2, force, w, e, Q, D
        self.grid = (grid, gd.Grid(10, 10))[grid is None]
        self.solid = (solid_cells,
                      np.zeros((Q, self.grid.width, self.grid.height)
                               ))[solid_cells is None]
        tau = (relaxation_time, tau)[relaxation_time is None]
        rho = (density, rho)[density is None]
        c = (lat_speed, self.grid.dx / time_step)[lat_speed is None]
        force = (ext_force, gd.LatticeVelocity(0, 0))[ext_force is None]
        w = initialize_weights()
        e = initialize_normal_velocities()
        # c_s = c / 3.
        iterations = int(t_final / time_step)
        self.f, self.u, self.rho = self.initialize(self.grid.width,
                                                   self.grid.height)
        self.bc = (boundary, gd.default_poiseuille_boundaries(self.grid))[boundary is None]
        # self.bc = (bd.Boundaries(self.grid, boundary),
        #            bd.Boundaries(self.grid, bd.default_poiseuille_boundaries(self.grid))
        #            )[boundary is None]
        start = timer()
        for i in range(0, iterations):
            self.macroscopic()
            self.external_force()
            self.collision()
            self.streaming()
            self.boundary_conditions()
        end = timer()
        print(str(end - start) + " sec.")

    @staticmethod
    def initialize(width, height):
        """
        Initialization of u, rho and f
        """
        f_0 = np.zeros((Q, height, width))
        for i in range(Q):
            f_0[i, :, :] = w[i] * rho
        rho_0 = np.zeros((Q, height, width))
        u_0 = gd.LatticeVelocity(np.zeros((height, width)),
                                 np.zeros((height, width)))
        return f_0, u_0, rho_0

    def macroscopic(self):
        """
        Calculates rho, ux, uy
        """
        temp_x = self.f * e.x[:, None, None]
        temp_y = self.f * e.y[:, None, None]

        self.rho = self.f.sum(axis=0)  # numpy function sum
        self.u.x = 1. / (self.rho ** 2) * temp_x.sum(axis=0)
        self.u.y = 1. / (self.rho ** 2) * temp_y.sum(axis=0)
        # todo add obstacles!

    def external_force(self):
        """
        Adds an external force
        """
        a = 1. / tau
        self.u.x += (force.x * self.rho) / (self.rho ** 2 * a)
        self.u.y += (force.y * self.rho) / (self.rho ** 2 * a)

    def collision(self):
        """
        Collision step
        """
        f_eq = self.compute_f_equilibrium()
        a = 1. / tau
        self.f = (1. - a) * self.f + a * f_eq

    def streaming(self):
        """
        Streaming step of LB
        """
        for i in range(Q):
            self.f[i] = np.roll(self.f[i], e.x[i], axis=1)  # roll along x
            self.f[i] = np.roll(self.f[i], e.y[i], axis=0)  # roll along y

    def boundary_conditions(self):
        bc_bounceback = [b for b in self.bc if b.type is 'bounceback']
        # bc_zoe_he_velocity = [b for b in self.bc if b.type is 'zoe_he_velocity']
        # bc_zoe_he_pressure = [b for b in self.bc if b.type is 'zoe_he_pressure']
        bc_neumann = [b for b in self.bc if b.type is 'neumann']
        if len(bc_bounceback) != 0:
            pass
            # todo bounceback type and zoe he bc
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
                                 e u      (e u)^2       u u
        f_eq = weight*rho*(1 + 3 --- + 9/2 ------ - 3/2 ---)
                                 c^2        c^4         c^2

        e u = e_x u_x + e_y u_y
        :return: equilibrium distribution function
        """
        xx = 1. * e.x[:, None, None] * e.x[:, None, None] * self.u.x * self.u.x
        xy = 1. * e.x[:, None, None] * e.y[:, None, None] * self.u.x * self.u.y
        yy = 1. * e.y[:, None, None] * e.y[:, None, None] * self.u.y * self.u.y
        term1 = (3. / c ** 2) * (e.x[:, None, None] * self.u.x +
                                 e.y[:, None, None] * self.u.y)
        term2 = 9. / (2. * c ** 4) * (xx + yy + 2. * xy)
        term3 = 3. / (2. * c ** 2) * (1. * self.u.x * self.u.x +
                                      1. * self.u.y * self.u.y)
        f_eq = w[:, None, None] * self.rho * (1 + term1 + term2 -
                                              term3)
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


if __name__ == "__main__":
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