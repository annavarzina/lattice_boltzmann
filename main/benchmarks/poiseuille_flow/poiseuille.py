

import main.schemes.d2q9 as lb
import main.schemes.grid as gd
import main.schemes.plot_and_save as ps
import os, inspect
import profile

def get_path():
    # return os.path.realpath(__file__)
    return os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

def horizontal_poiseuille_boundaries1(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=lv.LatticeVelocity(0, 0)))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=lv.LatticeVelocity(0, 0)))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'periodic'))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'periodic'))
    return boundaries


def horizontal_poiseuille_boundaries2(grid):
    boundaries = []
    v0 = gd.LatticeVelocity(0, 0)
    v1 = gd.LatticeVelocity(1e-6, 0)
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=v0))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=v1))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'periodic'))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'periodic'))
    return boundaries


def horizontal_poiseuille_boundaries3(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'bounceback'))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'bounceback'))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'periodic'))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'periodic'))
    return boundaries


if __name__ == "__main__":
    n = 60     # width  - 0x
    m = 30     # height - 0y
    rho = 1.   # density
    v = 0.1    # viscosity
    tau = v * 3. + 0.5   # relaxation time
    t = 3000  # final time
    gravity = gd.LatticeVelocity(1e-5, 0)  # acceleration by gravity

    userGrid = gd.Grid(n, m)
    bc = horizontal_poiseuille_boundaries2(userGrid)
    lat_bol = lb.D2Q9(grid=userGrid,
                      iterations=t,
                      relaxation_time=tau,
                      ext_force=gravity,
                      boundary=bc)

    path = get_path()
    ps.save_all(lat_bol,path=path)
    ps.plot_poiseuille_solution_0y(lat_bol)

    # profile.run('print \
    #              lb.D2Q9(grid=userGrid, \
    #                      iterations=t, \
    #                      relaxation_time=tau, \
    #                      ext_force=gravity, \
    #                      boundary=bc); \
    #              print')