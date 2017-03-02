# coding=utf-8
import unittest
import main.schemes.d2q9 as lb
import main.schemes.grid as gd
import numpy as np


class PipeFlowTest(unittest.TestCase):
    @staticmethod
    def testPressurePressurePipe():
        n = 50  # length  - 0x
        m = 300  # height - 0y
        v = 0.1  # viscosity
        p_left = 2. / 1.  # left pressure
        p_right = 1. / 1.  # right pressure
        t_f = 500  # final time
        tau = v * 3. + 0.5
        user_grid = gd.Grid(width=m, height=n)
        bc = pipe_pp_bc(grid=user_grid, p_right=p_right, p_left=p_left)
        lat_bol = lb.D2Q9(grid=user_grid,
                          iterations=t_f,
                          boundary=bc,
                          relaxation_time=tau,
                          plot_velocity=False)
        # ps.plot_streamlines(lat_bol, save=False)
        # ps.close_figures()

    @staticmethod
    def testVelocityPressurePipe():
        n = 50  # length  - 0x
        m = 300  # height - 0y
        v = 0.1  # viscosity
        v_left = 1. / 10.  # left velocity
        p_right = 1. / 1.  # right pressure
        t_f = 500  # final time
        tau = v * 3. + 0.5
        user_grid = gd.Grid(width=m, height=n)
        # boundary conditions
        bc = pipe_vp_bc(grid=user_grid, p_right=p_right, v_left=v_left)
        lat_bol = lb.D2Q9(grid=user_grid,
                          iterations=t_f,
                          boundary=bc,
                          relaxation_time=tau,
                          plot_velocity=False)

    @staticmethod
    def testVelocityVelocityPipe():
        n = 50  # length  - 0x
        m = 300  # height - 0y
        v = 0.1  # viscosity
        v_left = 1. / 10.  # left velocity
        v_right = 1. / 10.  # right pressure
        t_f = 500  # final time
        tau = v * 3. + 0.5
        user_grid = gd.Grid(width=m, height=n)
        # boundary conditions
        bc = pipe_vv_bc(grid=user_grid, v_right=v_right, v_left=v_left)
        lat_bol = lb.D2Q9(grid=user_grid,
                          iterations=t_f,
                          boundary=bc,
                          relaxation_time=tau,
                          plot_velocity=False)


class PoiseuilleFlowTest(unittest.TestCase):
    @staticmethod
    def testPoiseuilleFlowPeriodicZeroVelocities():
        n = 100  # width  - 0x
        m = 30  # height - 0y
        v = 0.1  # viscosity
        tau = v * 3. + 0.5  # relaxation time
        t_f = 500  # final time
        g = 1e-5
        gravity = gd.LatticeVelocity(g, 0)  # acceleration by gravity

        userGrid = gd.Grid(n, m)
        bc = horizontal_poiseuille_boundaries_zero_v(grid=userGrid)
        lat_bol = lb.D2Q9(grid=userGrid,
                          iterations=t_f,
                          relaxation_time=tau,
                          ext_force=gravity,
                          boundary=bc,
                          plot_velocity=False)

    @staticmethod
    def testPoiseuilleFlowPeriodicNonZeroTopVelocities():
        n = 100  # width  - 0x
        m = 30  # height - 0y
        v = 0.1  # viscosity
        tau = v * 3. + 0.5  # relaxation time
        t_f = 500  # final time
        v_top = 1e-3
        g = 1e-5
        gravity = gd.LatticeVelocity(g, 0)  # acceleration by gravity

        userGrid = gd.Grid(n, m)
        bc = horizontal_poiseuille_boundaries_top_vel(grid=userGrid, vel=v_top)
        lat_bol = lb.D2Q9(grid=userGrid,
                          iterations=t_f,
                          relaxation_time=tau,
                          ext_force=gravity,
                          boundary=bc,
                          plot_velocity=False)

    @staticmethod
    def testPoiseuilleFlowPeriodicBounceback():
        n = 100  # width  - 0x
        m = 30  # height - 0y
        v = 0.1  # viscosity
        tau = v * 3. + 0.5  # relaxation time
        t_f = 500  # final time
        g = 1e-5
        gravity = gd.LatticeVelocity(g, 0)  # acceleration by gravity

        userGrid = gd.Grid(n, m)
        bc = horizontal_poiseuille_boundaries_bb(grid=userGrid)
        lat_bol = lb.D2Q9(grid=userGrid,
                          iterations=t_f,
                          relaxation_time=tau,
                          ext_force=gravity,
                          boundary=bc,
                          plot_velocity=False)


class DrivenCavityTest(unittest.TestCase):

    def testDrivenCavity(self):
        n = 100  # width  - 0x
        m = 100  # length - 0y
        v = 0.1  # viscosity
        tau = v * 3. + 0.5  # relaxation time
        t_f = 500  # final time
        u_x = 0.1
        user_grid = gd.Grid(n, m)
        reynolds = n * u_x / v
        self.assertEqual(reynolds, n)
        bc = driven_cavity_boundaries(user_grid, u_x)
        lat_bol = lb.D2Q9(grid=user_grid,
                          iterations=t_f,
                          boundary=bc,
                          relaxation_time=tau,
                          plot_velocity=False)


class FlowPastCylinderTest(unittest.TestCase):
    @staticmethod
    def testPastCylinder():
        n = 20  # width  - 0x
        m = 40  # length - 0y
        rho = 1.  # density
        v = 0.03
        tau = v * 3. + 0.5  # relaxation time
        t_f = 500  # final time
        u = 0.3
        reynolds = str(m * u / v)
        solid = np.zeros((m, n))
        # circle
        radius = 6
        x0 = 20
        y0 = m / 2
        for i in range(0, n):
            for j in range(0, m):
                if (j - y0) ** 2 + (i - x0) ** 2 <= radius ** 2:
                    solid[j, i] = 1
        userGrid = gd.Grid(n, m)
        bc = flow_past_cylinder(userGrid, u, rho)
        lat_bol = lb.D2Q9(grid=userGrid,
                          iterations=t_f,
                          relaxation_time=tau,
                          boundary=bc,
                          solid_cells=solid,
                          plot_velocity=False)


if __name__ == "__main__":
    unittest.main()


# #### Boundary Conditions #### #


def pipe_pp_bc(grid, p_right, p_left):
    boundaries = []
    u = gd.LatticeVelocity(0, 0)  # zero velocity
    assert isinstance(p_right, float), "Incorrect type"
    assert isinstance(p_left, float), "Incorrect type"

    # for i in range(0, grid.width):
    #     boundaries.append(bd.Boundary(top_faces[i], 'zoe_he_pressure',value=p_left))
    #     boundaries.append(bd.Boundary(bottom_faces[i],  'zoe_he_pressure',value=p_right))
    # for i in range(0, grid.height):
    #     boundaries.append(bd.Boundary(left_faces[i], 'zoe_he_velocity', value=lv.LatticeVelocity(0, 0)))
    #     boundaries.append(bd.Boundary(right_faces[i], 'zoe_he_velocity', value=lv.LatticeVelocity(0, 0)))

    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u))
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_pressure', value=p_left))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_pressure', value=p_right))
    return boundaries


def pipe_vp_bc(grid, p_right, v_left):
    boundaries = []
    u = gd.LatticeVelocity(0, 0)  # zero velocity
    assert isinstance(p_right, float), "Incorrect type"
    assert isinstance(v_left, float), "Incorrect type"

    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u))
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(v_left, 0)))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_pressure', value=p_right))
    return boundaries


def pipe_vv_bc(grid, v_right, v_left):
    boundaries = []
    u = gd.LatticeVelocity(0, 0)  # zero velocity
    assert isinstance(v_right, float), "Incorrect type"
    assert isinstance(v_left, float), "Incorrect type"

    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u))
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(v_left, 0)))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(v_right, 0)))
    return boundaries


def horizontal_poiseuille_boundaries_zero_v(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(0, 0)))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(0, 0)))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'periodic'))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'periodic'))
    return boundaries


def horizontal_poiseuille_boundaries_top_vel(grid, vel):
    boundaries = []
    v0 = gd.LatticeVelocity(0, 0)
    v1 = gd.LatticeVelocity(vel, 0)
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=v0))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=v1))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'periodic'))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'periodic'))
    return boundaries


def horizontal_poiseuille_boundaries_bb(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'bounceback'))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'bounceback'))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'periodic'))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'periodic'))
    return boundaries


def driven_cavity_boundaries(grid, velocity):
    """
    top - slip boundary (x-direction velocity)
    bottom, left, right - zero velocity
    """
    boundaries = []
    v0 = gd.LatticeVelocity(0, 0)
    v_top = gd.LatticeVelocity(velocity, 0)
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=v0))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=v_top))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=v0))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_velocity', value=v0))
    # ======================================================================= #
    # for i in range(0, grid.width):
    #     boundaries.append(bd.Boundary(bottom_faces[i], 'bounceback'))
    #     boundaries.append(bd.Boundary(top_faces[i], 'zoe_he_velocity', value=v_top))
    # for i in range(0, grid.height):
    #     boundaries.append(bd.Boundary(left_faces[i], 'bounceback'))
    #     boundaries.append(bd.Boundary(right_faces[i], 'bounceback'))
    return boundaries


def flow_past_cylinder(grid, velocity, pressure):
    boundaries = []
    u_hor = gd.LatticeVelocity(velocity, 0)
    u_0 = gd.LatticeVelocity(0, 0)
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u_0))
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u_0))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=u_hor))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_pressure', value=pressure))
    return boundaries