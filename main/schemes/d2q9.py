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


class D2Q9:
    def __init__(self,
                 grid=None,
                 iterations=1000,
                 solid_cells=None,
                 relaxation_time=None,
                 density=None,
                 ext_force=None,
                 boundary=None,
                 reynolds=None,
                 plot_velocity=True,
                 path=None):

        global Q, D, e, w, c, c_s2, dt
        w = self.initialize_weights()
        e = self.initialize_normal_velocities()
        self.grid = (grid, gd.Grid(10, 10))[grid is None]
        self.solid = (solid_cells, np.zeros((self.grid.height, self.grid.width)))[solid_cells is None]
        self.bc = (boundary, gd.default_poiseuille_boundaries(self.grid))[boundary is None]
        self.tau = (relaxation_time, 1)[relaxation_time is None]  # short 'if-else' statement
        self.initial_density = (density, 1)[density is None]
        self.applied_force = (ext_force, gd.LatticeVelocity(0., 0.))[ext_force is None]
        self.viscosity = (1. / 3.) * (self.tau - 0.5)
        self.reynolds = (reynolds, 0)[reynolds is None]
        self.f, self.u, self.rho = self.initialize(self.grid.width, self.grid.height)
        self.stream = np.zeros((self.grid.height, self.grid.width))
        self.vorticity = np.zeros((self.grid.height, self.grid.width))
        if path is None:
            path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        title = ''
        elapsed_time=0
        start = timer()
        for i in range(0, iterations):
            self.compute_macroscopic()
            if (self.applied_force.x != 0) | (self.applied_force.y != 0):
                self.external_force()
            self.collision_step()
            self.streaming_step()
            self.boundary_conditions()

            if (i % 10 == 0) & (i > 0):
                print(i)
                if i % 100 == 0:
                    title = "Re_" + str(self.reynolds) + "_vis_" + str(self.viscosity) + \
                             "_grid_" + str(self.grid.width) + 'x' + str(self.grid.height) + str(i)
                    ps.plot_streamlines(self, title=title, path=path, save=True)
                if plot_velocity:
                    ps.plot_ux(self)
                    ps.plot_uy(self)

        elapsed_time=timer() - start
        print(str(elapsed_time) + " sec.")
        # self.velocity_stream_and_vorticity()
        ps.flag = False
        # ps.save_data(title)

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

    def compute_macroscopic(self):
        """ rho, ux, uy """
        temp_x = self.f * e.x[:, None, None]
        temp_y = self.f * e.y[:, None, None]
        self.rho = self.f.sum(axis=0)
        self.u.x = (self.rho / (1e-30 + self.rho) ** 2) * temp_x.sum(axis=0)
        self.u.y = (self.rho / (1e-30 + self.rho) ** 2) * temp_y.sum(axis=0)
        # account for obstacles
        self.rho *= 1. * (self.solid <= 0)
        self.u.x *= 1. * (self.solid <= 0)
        self.u.y *= 1. * (self.solid <= 0)
        self.check_cfl()

    def collision_step(self):
        f_eq = self.compute_f_equilibrium()
        a = dt / self.tau
        temp_f = (1. - a) * self.f + a * f_eq
        f_bb = np.copy(self.f[:, :, :])
        i = [0, 3, 4, 1, 2, 7, 8, 5, 6]  # bounceback order
        f_bb[i, None, None] = self.f[:, None, None]
        self.f = temp_f * 1. * (self.solid <= 0) + f_bb * 1. * (self.solid > 0)

    def streaming_step(self):
        # todo add obstacles (not interior cell)
        # todo add bounceback?
        for i in range(Q):
            self.f[i] = np.roll(self.f[i], int(self.grid.dx) * e.x[i], axis=1)  # roll along x
            self.f[i] = np.roll(self.f[i], int(self.grid.dx) * e.y[i], axis=0)  # roll along y

    def external_force(self):
        self.u.x += (self.applied_force.x * self.rho * self.tau) / ((self.rho + 1e-30) ** 2)
        self.u.y += (self.applied_force.y * self.rho * self.tau) / ((self.rho + 1e-30) ** 2)

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

        f_eq = w[:, None, None] * self.rho * (1 + term1 + term2 - term3)
        return f_eq

    def boundary_conditions(self):
        bc_bounceback = [b for b in self.bc if b.type is 'bounceback']
        bc_zoe_he_velocity = [b for b in self.bc if b.type is 'zoe_he_velocity']
        bc_zoe_he_pressure = [b for b in self.bc if b.type is 'zoe_he_pressure']

        if len(bc_bounceback) != 0:
            temp_f = np.copy(self.f[:, :, :])
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

        if len(bc_zoe_he_velocity):
            bc_bottom, bc_top, bc_left, bc_right = self.define_boundaries(bc_zoe_he_velocity)
            if len(bc_top) != 0:
                rows, columns = self.define_row_columns(bc_top, 1)
                u0_x = [bt.value.x for bt in bc_top]
                u0_y = [bt.value.y for bt in bc_top]
                rho0 = ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                        (self.f[2, rows, columns] + self.f[5, rows, columns] + self.f[6, rows, columns])
                        ) / (1. + np.array(u0_y))
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
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
                u0_x = [bt.value.x for bt in bc_bottom]
                u0_y = [bt.value.y for bt in bc_bottom]
                rho0 = ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                        (self.f[4, rows, columns] + self.f[7, rows, columns] + self.f[8, rows, columns])
                        ) / (1. - np.array(u0_y))
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
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
                u0_x = [bt.value.x for bt in bc_right]
                u0_y = [bt.value.y for bt in bc_right]
                rho0 = ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]) + 2. *
                        (self.f[1, rows, columns] + self.f[5, rows, columns] + self.f[8, rows, columns])
                        ) / (1. + np.array(u0_x))
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
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
                u0_x = [bt.value.x for bt in bc_left]
                u0_y = [bt.value.y for bt in bc_left]
                rho0 = ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]
                         ) + 2. * (self.f[3, rows, columns] + self.f[7, rows, columns] + self.f[6, rows, columns])
                        ) / (1. - np.array(u0_x))
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
                f13 = (2. / 3.) * ru_x / (c ** 2)
                self.f[1, rows, columns] = self.f[3, rows, columns] + f13
                self.f[5, rows, columns] = self.f[7, rows, columns] + (1. / 2.) * ru_x - \
                    (1. / 2.) * (self.f[2, rows, columns] - self.f[4, rows, columns]) + \
                    (1. / 2.) * ru_y - f13
                self.f[8, rows, columns] = self.f[6, rows, columns] + (1. / 2.) * ru_x - \
                    (1. / 2.) * (self.f[4, rows, columns] - self.f[2, rows, columns]) - \
                    (1. / 2.) * ru_y - f13

        if len(bc_zoe_he_pressure):
            bc_bottom, bc_top, bc_left, bc_right = self.define_boundaries(bc_zoe_he_pressure)
            if len(bc_top) != 0:
                rows, columns = self.define_row_columns(bc_top, 1)
                rho0 = [bt.value for bt in bc_top]
                # u0_x = 3. / 2. * (self.f[1, rows, columns] - self.f[3, rows, columns]) / np.array(rho0)
                u0_x = np.zeros(len(rho0))
                u0_y = -1 + ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                             (self.f[2, rows, columns] + self.f[5, rows, columns] + self.f[6, rows, columns])
                             ) / np.array(rho0)
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
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
                rho0 = [bt.value for bt in bc_bottom]
                # u0_x = 3. / 2. * (self.f[1, rows, columns] - self.f[3, rows, columns]) / np.array(rho0)
                u0_x = np.zeros(len(rho0))
                u0_y = 1 - ((self.f[0, rows, columns] + self.f[1, rows, columns] + self.f[3, rows, columns]) + 2 *
                            (self.f[4, rows, columns] + self.f[7, rows, columns] + self.f[8, rows, columns])
                            ) / np.array(rho0)
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
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
                rho0 = [bt.value for bt in bc_right]
                u0_x = -1. + ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]) + 2. *
                              (self.f[1, rows, columns] + self.f[5, rows, columns] + self.f[8, rows, columns])
                              ) / np.array(rho0)
                # u0_y = 3. / 2. * (self.f[4, rows, columns] - self.f[2, rows, columns]) / np.array(rho0)
                u0_y = np.zeros(len(rho0))
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
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
                rho0 = [bt.value for bt in bc_left]
                u0_x = 1 - ((self.f[0, rows, columns] + self.f[2, rows, columns] + self.f[4, rows, columns]) + 2. *
                            (self.f[3, rows, columns] + self.f[6, rows, columns] + self.f[7, rows, columns])
                            ) / np.array(rho0)
                # u0_y = 3. / 2. * (self.f[4, rows, columns] - self.f[2, rows, columns]) / np.array(rho0)
                u0_y = np.zeros(len(rho0))
                ru_x = rho0 * np.array(u0_x)
                ru_y = rho0 * np.array(u0_y)
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

    def velocity_stream_and_vorticity(self):
        # todo return vorticity and stream, save them
        u = np.array(self.u.x)
        v = np.array(self.u.y)
        int_u_x = it.simps(u, axis=1)
        int_u_y = it.simps(u, axis=0)
        int_v_x = it.simps(v, axis=1)
        int_v_y = it.simps(v, axis=0)
        self.stream = int_u_y[np.newaxis, :] - int_v_x[:, np.newaxis]  # stream function (psi)
        self.vorticity = int_u_x[:, np.newaxis] + int_v_y[np.newaxis, :]  # vorticity


if __name__ == "__main__":
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
