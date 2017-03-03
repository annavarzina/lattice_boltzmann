# coding=utf-8
import main.schemes.d2q9 as lb
import main.schemes.grid as gd
import numpy as np
import main.schemes.plot_and_save as ps
import os, inspect

def flow_past_cylinder_bc(grid, velocity, pressure):
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


if __name__ == "__main__":
    n = 500                    # width  - 0x
    m = 50                     # height - 0y
    rho = 1.                    # density
    v = 0.03                    # viscosity
    tau = v * 3. + 0.5          # relaxation time
    t = 1000                    # final time
    u = 0.3                     # inflow velocity
    print('Re = ' + str(m*u/v))
    solid = np.zeros((m, n))    # circle shape
    radius = m/8                 # radius of circle
    x0 = n/7                    # (x0,y0) - center of the circle
    y0 = m/2
    for i in range(0, n):
        for j in range(0,m):
            if (j - y0) ** 2 + (i - x0) ** 2 <= radius ** 2:
                solid[j, i] = 1
    # solid[int(m/2)-20:int(m/2)+20, int(n/10)-20:int(n/10)+20] = 1  # square
    path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    userGrid = gd.Grid(n, m)
    bc = flow_past_cylinder_bc(userGrid, u, rho)
    lat_bol = lb.D2Q9(grid=userGrid,
                      iterations=t,
                      relaxation_time=tau,
                      boundary=bc,
                      solid_cells=solid,
                      path=path)
    ps.save_all(lat_bol, path=path)
    ps.plot_streamlines(lat_bol, path=path)
    # print(isinstance(gravity, lv.LatticeVelocity))

    # profile.run('print \
    #              lb.D2Q9(grid=userGrid, \
    #              iterations=t, \
    #              relaxation_time=tau, \
    #              boundary=bc, \
    #              solid_cells=solid); print')
