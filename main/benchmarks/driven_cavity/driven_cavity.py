# coding=utf-8
import main.schemes.d2q9 as lb
import main.schemes.grid as gd
import main.schemes.plot_and_save as ps
import os, inspect


def driven_cavity_boundaries(grid, value):
    boundaries = []
    v0 = gd.LatticeVelocity(0,0)
    for i in range(0, grid.width):
        # boundaries.append(Boundary(top_faces[i], 'zoe_he_velocity', value=value))
        # boundaries.append(Boundary(bottom_faces[i], 'zoe_he_velocity', value=v0))
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=v0))
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=value))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=v0))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_velocity', value=v0))
    # ======================================================================= #
    # for i in range(0, grid.width):
    #     boundaries.append(Boundary(bottom_faces[i], 'bounceback'))
    #     boundaries.append(Boundary(top_faces[i], 'zoe_he_velocity', value=value))
    # for i in range(0, grid.height):
    #     boundaries.append(Boundary(left_faces[i], 'bounceback'))
    #     boundaries.append(Boundary(right_faces[i], 'bounceback'))
    return boundaries


if __name__ == "__main__":
    # initial parameters
    n = 500   # width  - 0x
    m = 500   # length - 0y
    rho = 1.   # density
    viscosity = 0.125
    tau = viscosity * 3. + 0.5
    t = 60000  # final time
    u_x = 0.25
    re = n * u_x / viscosity
    print(re)
    velocity = gd.LatticeVelocity(u_x, 0)
    userGrid = gd.Grid(n, m)
    bc = driven_cavity_boundaries(userGrid, velocity)
    path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    lat_bol = lb.D2Q9(grid=userGrid,
                      iterations=t,
                      boundary=bc,
                      relaxation_time=tau,
                      reynolds=re,
                      path=path)
    print("Viscosity = " + str(lat_bol.viscosity))
    print("Reynolds number = " + str(u_x * rho * n / lat_bol.viscosity))

    # print(path)
    ps.save_all(lat_bol,path=path)
    ps.save_data(lat_bol, lat_bol.stream, 'stream',path=path)
    ps.save_data(lat_bol, lat_bol.vorticity, 'vorticity',path=path)
