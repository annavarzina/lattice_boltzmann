# coding=utf-8
import main.schemes.d2q9 as lb
import main.schemes.grid as gd
import main.schemes.plot_and_save as ps
import os, inspect


def pipe_pp_bc(grid, p_right, p_left):
    boundaries = []
    u = gd.LatticeVelocity(0, 0)  # zero velocity
    assert isinstance(p_right, float), "Incorrect type"
    assert isinstance(p_left, float), "Incorrect type"
    for i in range(0, grid.width):
        boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u))
        boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u))
    for i in range(0, grid.height):
        boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_pressure', value=p_left))
        boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_pressure', value=p_right))
    return boundaries

if __name__ == "__main__":
    n = 30  # length  - 0x
    m = 200  # height - 0y
    v = 0.1  # viscosity
    p_left = 2. / 1.  # left pressure
    p_right = 1. / 1.  # right pressure
    t_f = 500  # final time
    tau = v * 3. + 0.5
    user_grid = gd.Grid(width=m, height=n)
    # boundary conditions
    bc = pipe_pp_bc(grid=user_grid, p_right=p_right, p_left=p_left)
    lat_bol = lb.D2Q9(grid=user_grid,
                      iterations=t_f,
                      boundary=bc,
                      relaxation_time=tau)


    path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    print(path)
    ps.save_all(lat_bol,path=path)
    ps.plot_streamlines(lat_bol, path=path)