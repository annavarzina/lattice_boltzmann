"""
Module 'd2q9LB' contains a class D2Q9, which represents logics of
lattice-boltzmann method.
"""

# todo plot function and use tkinter to visualize
import inspect, os
import grid as gd
import numpy as np
import scipy.integrate as it
import plot_and_save as ps
from timeit import default_timer as timer

Q = 9                                       # directions
D = 2                                       # dimensions
e = gd.LatticeVelocity(0, 0)                # lattice velocities
w = []                                      # weights
c = 1.                                      # propagation speed
c_s2 = 1. / 3.                              # squared speed of sound
dt = 1
conc = 1.  # default density
flux = gd.LatticeVelocity(0, 0) # flux
diff = 1. # default diffusion value


class D2Q9:
    def __init__(self,
                 grid=None,
                 iterations=1000,
                 eps=None,
                 check_convergence=False,
                 solid_cells=None,
                 relaxation_time=None,
                 concentration=None,
                 diffusion = None,
                 velocity = None,
                 ext_force=None,
                 boundary=None,
                 path=None):

        global Q, D, e, w, c, c_s2, dt, conc, flux, diff
        w = self.initialize_weights()
        e = self.initialize_normal_velocities()
        self.iterations = iterations
        self.grid = (grid, gd.Grid(10, 10))[grid is None]
        self.solid = (solid_cells, np.zeros((self.grid.height, self.grid.width)))[solid_cells is None]
        self.bc = (boundary, gd.default_poiseuille_boundaries(self.grid))[boundary is None]
        self.tau = (relaxation_time, 1)[relaxation_time is None]  # short 'if-else' statement
        self.source = (ext_force, gd.LatticeVelocity(0., 0.))[ext_force is None]
        self.f = np.zeros((Q, self.grid.height, self.grid.width))
        self.diff = (diffusion, diff)[diffusion is None]
        # solution initialization
        self.conc = (concentration, conc)[concentration is None]
        print(self.conc)
        self.u = (velocity, gd.LatticeVelocity(np.zeros((self.grid.height, self.grid.width)),
                                               np.zeros((self.grid.height, self.grid.width))))[velocity is None]
        self.flux = gd.LatticeVelocity(np.zeros((self.grid.height, self.grid.width)),
                                       np.zeros((self.grid.height, self.grid.width)))
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
            #'''
            print('Iteration %s' % (i + 1))
            print('Concentration:')
            print(self.conc)
            #'''
            #'''
            print('Flux X:')
            print(self.flux.x)
            print('Flux Y:')
            print(self.flux.y)
            #'''
            print('Liquid mass is %s [mol]' % np.sum(self.conc))
        end = timer()
        print(str(end - start) + " sec.")

    # ============================== Initialization functions ================================= #

    def initialize(self,width, height):
        f = np.zeros((Q, height, width))
        for i in range(Q):
            f[i, :, :] = w[i] * self.initial_density
        rho = np.ones((Q, height, width)) * self.initial_density
        u = gd.LatticeVelocity(np.zeros((height, width)), np.zeros((height, width)))
        return f, u, rho

    @staticmethod
    def initialize_normal_velocities():
        e_x = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
        e_y = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
        return gd.LatticeVelocity(e_x, e_y)

    @staticmethod
    def initialize_weights():
        return np.array([4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9.,
                         1. / 36., 1. / 36., 1. / 36., 1. / 36.])

    @staticmethod
    def define_boundaries(bc):
        bc_bottom = [bn for bn in bc if bn.position == 0]
        bc_top = [bn for bn in bc if bn.position == 1]
        bc_left = [bn for bn in bc if bn.position == 2]
        bc_right = [bn for bn in bc if bn.position == 3]
        return bc_bottom, bc_top, bc_left, bc_right

    def define_row_columns(self, bc, d):
        rows = [(bt.face - d) / 4 / self.grid.width for bt in bc]
        columns = [((bt.face - d) / 4) % self.grid.width for bt in bc]
        return rows, columns

    # ============================== Lattice Boltzmann main2 functions ================================= #

    def macroscopic(self):
        """ rho, ux, uy """
        temp_x = self.f * e.x[:, None, None]
        temp_y = self.f * e.y[:, None, None]
        self.conc = self.f.sum(axis=0)
        self.flux.x = (self.conc / (1e-30 + self.conc) ** 2) * temp_x.sum(axis=0)
        self.flux.y = (self.conc / (1e-30 + self.conc) ** 2) * temp_y.sum(axis=0)
        # account for obstacles
        self.conc *= 1. * (self.solid <= 0)
        self.flux.x *= 1. * (self.solid <= 0)
        self.flux.y *= 1. * (self.solid <= 0)
        self.check_cfl()

    def collision(self):
        f_eq = self.compute_f_equilibrium()
        a = dt / self.tau
        temp_f = (1. - a) * self.f + a * f_eq
        f_bb = np.copy(self.f[:, :, :])
        i = [0, 3, 4, 1, 2, 7, 8, 5, 6]  # bounceback order
        f_bb[i, None, None] = self.f[:, None, None]
        self.f = temp_f * 1. * (self.solid <= 0) + f_bb * 1. * (self.solid > 0)

    def streaming(self):
        # todo add obstacles (not interior cell)
        # todo add bounceback?
        for i in range(Q):
            self.f[i] = np.roll(self.f[i], int(self.grid.dx) * e.x[i], axis=1)  # roll along x
            self.f[i] = np.roll(self.f[i], int(self.grid.dx) * e.y[i], axis=0)  # roll along y

    def source(self):
        self.flux.x += (self.applied_force.x * self.conc * self.tau) / ((self.conc + 1e-30) ** 2)
        self.flux.y += (self.applied_force.y * self.conc * self.tau) / ((self.conc + 1e-30) ** 2)

    def compute_f_equilibrium(self):
        """
            Equilibrium distribution function:
                                     e u      (e u)^2       u u
            f_eq = weight*rho*(1 + 3 --- + 9/2 ------ - 3/2 ---)
                                     c^2        c^4         c^2
            e u = e_x u_x + e_y u_y
        """
        xx = 1. * e.x[:, None, None] * e.x[:, None, None] * self.u.x * self.u.x
        xy = 1. * e.x[:, None, None] * e.y[:, None, None] * self.u.x * self.u.y
        yy = 1. * e.y[:, None, None] * e.y[:, None, None] * self.u.y * self.u.y
        term1 = (1. / c_s2) * (e.x[:, None, None] * self.u.x + e.y[:, None, None] * self.u.y)
        term2 = 0.5 / (c_s2 ** 2) * (xx + yy + 2. * xy)
        term3 = (0.5 / c_s2) * (1. * self.u.x * self.u.x + 1. * self.u.y * self.u.y)

        f_eq = w[:, None, None] * self.conc * (1 + term1 + term2 - term3)
        return f_eq

    def boundary_conditions(self):
        bc_bounceback = [b for b in self.bc if b.type is 'bounceback']
        bc_zoe_he_flux = [b for b in self.bc if b.type is 'zoe_he_flux']
        bc_zoe_he_conc = [b for b in self.bc if b.type is 'zoe_he_conc']

        if len(bc_bounceback) != 0:
            temp_f = np.copy(self.f[:, :, :])
            bc_bottom, bc_top, bc_left, bc_right = self.define_boundaries(bc_bounceback)
            temp_f = np.copy(self.f[:, :, :])
            if len(bc_top) != 0:
                rows, columns = self.define_row_columns(bc_top, 1)
                self.f[4, rows, columns] = temp_f[2, rows, columns]
                self.f[7, rows, columns] = temp_f[5, rows, columns]
                self.f[8, rows, columns] = temp_f[6, rows, columns]
            if len(bc_bottom) != 0:
                rows, columns = self.define_row_columns(bc_bottom, 0)
                self.f[2, rows, columns] = temp_f[4, rows, columns]
                self.f[5, rows, columns] = temp_f[7, rows, columns]
                self.f[6, rows, columns] = temp_f[8, rows, columns]
            if len(bc_right) != 0:
                rows, columns = self.define_row_columns(bc_right, 3)
                self.f[3, rows, columns] = temp_f[1, rows, columns]
                self.f[6, rows, columns] = temp_f[8, rows, columns]
                self.f[7, rows, columns] = temp_f[5, rows, columns]
            if len(bc_left) != 0:
                rows, columns = self.define_row_columns(bc_left, 2)
                self.f[1, rows, columns] = temp_f[3, rows, columns]
                self.f[5, rows, columns] = temp_f[7, rows, columns]
                self.f[8, rows, columns] = temp_f[6, rows, columns]
            pass
            '''
            rows, columns = self.define_row_columns(bc_bounceback, 0)
            # rows = [bt.face / 4 / self.grid.width for bt in bc_bounceback]
            # columns = [(bt.face/ 4) % self.grid.width for bt in bc_bounceback]
            # i = [0, 3, 4, 1, 2, 7, 8, 5, 6]  # bounceback order
            self.f[1, rows, columns] = temp_f[3, rows, columns]
            self.f[2, rows, columns] = temp_f[4, rows, columns]
            self.f[3, rows, columns] = temp_f[1, rows, columns]
            self.f[4, rows, columns] = temp_f[2, rows, columns]
            self.f[5, rows, columns] = temp_f[7, rows, columns]
            self.f[6, rows, columns] = temp_f[8, rows, columns]
            self.f[7, rows, columns] = temp_f[5, rows, columns]
            self.f[8, rows, columns] = temp_f[6, rows, columns]
            '''

        if len(bc_zoe_he_flux):
            bc_bottom, bc_top, bc_left, bc_right = self.define_boundaries(bc_zoe_he_flux)
            if len(bc_top) != 0:
                rows, columns = self.define_row_columns(bc_top, 1)
                flux0_x = [bt.value.x for bt in bc_top]
                flux0_y = [bt.value.y for bt in bc_top]
                conc0 = ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                        (self.f[2, rows, columns] + self.f[5, rows, columns] + self.f[6, rows, columns])
                        ) / (1. + np.array(flux0_y))
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f24 = (2. / 3.) * ru_y / (c ** 2)
                self.f[4, rows, columns] = self.f[2, rows, columns] - f24
                self.f[7, rows, columns] = self.f[5, rows, columns] - (1. / 2.) * ru_y + \
                    (1. / 2.) * (self.f[1, rows, columns] - self.f[3, rows, columns]) - \
                    (1. / 2.) * ru_x + f24
                self.f[8, rows, columns] = self.f[6, rows, columns] - (1. / 2.) * ru_y + \
                    (1. / 2.) * (self.f[3, rows, columns] - self.f[1, rows, columns]) + \
                    (1. / 2.) * ru_x + f24
            if len(bc_bottom) != 0:
                rows, columns = self.define_row_columns(bc_bottom, 0)
                flux0_x = [bt.value.x for bt in bc_bottom]
                flux0_y = [bt.value.y for bt in bc_bottom]
                conc0 = ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                        (self.f[4, rows, columns] + self.f[7, rows, columns] + self.f[8, rows, columns])
                        ) / (1. - np.array(flux0_y))
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f24 = (2. / 3.) * ru_y / (c ** 2)
                self.f[2, rows, columns] = self.f[4, rows, columns] + f24
                self.f[5, rows, columns] = self.f[7, rows, columns] + (1. / 2.) * ru_y - \
                    (1. / 2.) * (self.f[1, rows, columns] - self.f[3, rows, columns]) + \
                    (1. / 2.) * ru_x - f24
                self.f[6, rows, columns] = self.f[8, rows, columns] + (1. / 2.) * ru_y - \
                    (1. / 2.) * (self.f[3, rows, columns] - self.f[1, rows, columns]) - \
                    (1. / 2.) * ru_x - f24
            if len(bc_right) != 0:
                rows, columns = self.define_row_columns(bc_right, 3)
                flux0_x = [bt.value.x for bt in bc_right]
                flux0_y = [bt.value.y for bt in bc_right]
                conc0 = ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]) + 2. *
                        (self.f[1, rows, columns] + self.f[5, rows, columns] + self.f[8, rows, columns])
                        ) / (1. + np.array(flux0_x))
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f13 = (2. / 3.) * ru_x / (c ** 2)
                self.f[3, rows, columns] = self.f[1, rows, columns] - f13
                self.f[7, rows, columns] = self.f[5, rows, columns] - (1. / 2.) * ru_x + \
                    (1. / 2.) * (self.f[2, rows, columns] - self.f[4, rows, columns]) - \
                    (1. / 2.) * ru_y + f13
                self.f[6, rows, columns] = self.f[8, rows, columns] - (1. / 2.) * ru_x + \
                    (1. / 2.) * (self.f[4, rows, columns] - self.f[2, rows, columns]) + \
                    (1. / 2.) * ru_y + f13
            if len(bc_left) != 0:
                rows, columns = self.define_row_columns(bc_left, 2)
                flux0_x = [bt.value.x for bt in bc_left]
                flux0_y = [bt.value.y for bt in bc_left]
                conc0 = ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]
                         ) + 2. * (self.f[3, rows, columns] + self.f[7, rows, columns] + self.f[6, rows, columns])
                        ) / (1. - np.array(flux0_x))
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f13 = (2. / 3.) * ru_x / (c ** 2)
                self.f[1, rows, columns] = self.f[3, rows, columns] + f13
                self.f[5, rows, columns] = self.f[7, rows, columns] + (1. / 2.) * ru_x - \
                    (1. / 2.) * (self.f[2, rows, columns] - self.f[4, rows, columns]) + \
                    (1. / 2.) * ru_y - f13
                self.f[8, rows, columns] = self.f[6, rows, columns] + (1. / 2.) * ru_x - \
                    (1. / 2.) * (self.f[4, rows, columns] - self.f[2, rows, columns]) - \
                    (1. / 2.) * ru_y - f13

        if len(bc_zoe_he_conc):
            bc_bottom, bc_top, bc_left, bc_right = self.define_boundaries(bc_zoe_he_conc)
            if len(bc_top) != 0:
                rows, columns = self.define_row_columns(bc_top, 1)
                conc0 = [bt.value for bt in bc_top]
                # u0_x = 3. / 2. * (self.f[1, rows, columns] - self.f[3, rows, columns]) / np.array(rho0)
                flux0_x = np.zeros(len(conc0))
                flux0_y = -1 + ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                             (self.f[2, rows, columns] + self.f[5, rows, columns] + self.f[6, rows, columns])
                             ) / np.array(conc0)
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f24 = (2. / 3.) * ru_y / (c ** 2)

                self.f[4, rows, columns] = self.f[2, rows, columns] - f24
                self.f[7, rows, columns] = self.f[5, rows, columns] - (1. / 2.) * ru_y + \
                    (1. / 2.) * (self.f[1, rows, columns] - self.f[3, rows, columns]) - \
                    (1. / 2.) * ru_x + f24
                self.f[8, rows, columns] = self.f[6, rows, columns] - (1. / 2.) * ru_y + \
                    (1. / 2.) * (self.f[3, rows, columns] - self.f[1, rows, columns]) + \
                    (1. / 2.) * ru_x + f24
            if len(bc_bottom) != 0:
                rows, columns = self.define_row_columns(bc_bottom, 0)
                conc0 = [bt.value for bt in bc_bottom]
                # u0_x = 3. / 2. * (self.f[1, rows, columns] - self.f[3, rows, columns]) / np.array(rho0)
                flux0_x = np.zeros(len(conc0))
                flux0_y = 1 - ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                            (self.f[4, rows, columns] + self.f[7, rows, columns] + self.f[8, rows, columns])
                            ) / np.array(conc0)
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f24 = (2. / 3.) * ru_y / (c ** 2)

                self.f[2, rows, columns] = self.f[4, rows, columns] + f24
                self.f[5, rows, columns] = self.f[7, rows, columns] + (1. / 2.) * ru_y - \
                    (1. / 2.) * (self.f[1, rows, columns] - self.f[3, rows, columns]) + \
                    (1. / 2.) * ru_x - f24
                self.f[6, rows, columns] = self.f[8, rows, columns] + (1. / 2.) * ru_y - \
                    (1. / 2.) * (self.f[3, rows, columns] - self.f[1, rows, columns]) - \
                    (1. / 2.) * ru_x - f24
            if len(bc_right) != 0:
                rows, columns = self.define_row_columns(bc_right, 3)
                conc0 = [bt.value for bt in bc_right]
                flux0_x = -1. + ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]) + 2. *
                              (self.f[1, rows, columns] + self.f[5, rows, columns] + self.f[8, rows, columns])
                              ) / np.array(conc0)
                # u0_y = 3. / 2. * (self.f[4, rows, columns] - self.f[2, rows, columns]) / np.array(rho0)
                flux0_y = np.zeros(len(conc0))
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f13 = (2. / 3.) * ru_x / (c ** 2)

                self.f[3, rows, columns] = self.f[1, rows, columns] - f13
                self.f[7, rows, columns] = self.f[5, rows, columns] - (1. / 2.) * ru_x + \
                    (1. / 2.) * (self.f[2, rows, columns] - self.f[4, rows, columns]) - \
                    (1. / 2.) * ru_y + f13
                self.f[6, rows, columns] = self.f[8, rows, columns] - (1. / 2.) * ru_x + \
                    (1. / 2.) * (self.f[4, rows, columns] - self.f[2, rows, columns]) + \
                    (1. / 2.) * ru_y + f13
            if len(bc_left) != 0:
                rows, columns = self.define_row_columns(bc_left, 2)
                conc0 = [bt.value for bt in bc_left]
                flux0_x = 1 - ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]) + 2. *
                            (self.f[3, rows, columns] + self.f[6, rows, columns] + self.f[7, rows, columns])
                            ) / np.array(conc0)
                # u0_y = 3. / 2. * (self.f[4, rows, columns] - self.f[2, rows, columns]) / np.array(rho0)
                flux0_y = np.zeros(len(conc0))
                ru_x = conc0 * np.array(flux0_x)
                ru_y = conc0 * np.array(flux0_y)
                f13 = (2. / 3.) * ru_x / (c ** 2)

                self.f[1, rows, columns] = self.f[3, rows, columns] + f13
                self.f[5, rows, columns] = self.f[7, rows, columns] + (1. / 2.) * ru_x - \
                    (1. / 2.) * (self.f[2, rows, columns] - self.f[4, rows, columns]) + \
                    (1. / 2.) * ru_y - f13
                self.f[8, rows, columns] = self.f[6, rows, columns] + (1. / 2.) * ru_x - \
                    (1. / 2.) * (self.f[4, rows, columns] - self.f[2, rows, columns]) - \
                    (1. / 2.) * ru_y - f13

    # ============================== Additional functions ================================= #

    def check_cfl(self):
        u = max(np.average(abs(self.u.x)), np.average(abs(self.u.y)))
        c_num = u * dt / self.grid.dx  # Courant number
        assert c_num <= 1, "Do not satisfy CFL condition c = " + str(c_num)

    def converge(self, ux, uy, rho, eps):
        ux -= self.u.x
        print 'ux ' + str(np.average(np.abs(ux)))
        if np.average(np.abs(ux)) < eps:
            uy -= self.u.y
            print 'uy ' + str(np.average(np.abs(uy)))
            if np.average(np.abs(uy)) < eps:
                rho -= self.rho
                return True
                # print 'rho ' + str(np.average(np.abs(rho)))
                # if np.average(np.abs(rho)) < eps:
                #     return True
        return False


def concentration_boundaries(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_conc', value=0.))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_conc', value=0.))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_conc', value=0.))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_conc', value=0.))
    return boundaries


def zero_flux_boundaries(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_flux', value=gd.LatticeVelocity(0, 0)))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_flux', value=gd.LatticeVelocity(0, 0)))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_flux', value=gd.LatticeVelocity(0, 0)))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_flux', value=gd.LatticeVelocity(0, 0)))
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
    n = 5 # width  - lx
    m = 5 # height - ly
    t = 5 # final time
    c_init = np.zeros((n, m))
    c_init[2, 2] = 1

    solid = np.zeros((Q, n, m))
    userGrid = gd.Grid(n, m)
    bc = bb_boundaries(userGrid)
    #bc = zero_flux_boundaries(userGrid)
    # bc = concentration_boundaries(userGrid)
    lat_bol = D2Q9(grid=userGrid,
                   iterations=t,
                   concentration=c_init,
                   boundary=bc)
    '''
    n = 30              # width  - lx
    m = 60              # height - ly
    t = 3000            # final time
    gravity = gd.LatticeVelocity(0, 1e-8)  # acceleration by gravity
    solid = np.zeros((Q, n, m))

    userGrid = gd.Grid(n, m)
    lat_bol = D2Q9(grid=userGrid,
                   iterations=t,
                   ext_force=gravity,
                   plot_velocity=False)
    ps.plot_poiseuille_solution_0x(lat_bol)
    '''

