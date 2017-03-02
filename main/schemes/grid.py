import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


class Grid:
    def __init__(self, width, height):
        self.width = width    # 0x length (N)
        self.height = height  # 0y length (M)
        self.dx = 1           # lattice unit is constant
        self.bottom_faces = self.set_bottom_faces()
        self.top_faces = self.set_top_faces()
        self.left_faces = self.set_left_faces()
        self.right_faces = self.set_right_faces()

    def set_bottom_faces(self):
        n = self.width
        faces = []
        # direction = 0 for bottom
        for i in range(0, n):
            faces.append(4 * i)  # cell_no = n * j + i
        return faces

    def set_top_faces(self):
        n = self.width
        m = self.height
        faces = []
        # direction = 1 for top
        d = 1
        for i in range(0, n):
            faces.append(d + 4 * (n * (m - 1) + i))  # cell_no = n * j + i
        return faces

    def set_left_faces(self):
        n = self.width
        m = self.height
        faces = []
        # direction = 2 for left
        d = 2
        for j in range(0, m):
            faces.append(d + 4 * n * j)  # cell_no = n * j + i
        return faces

    def set_right_faces(self):
        n = self.width
        m = self.height
        faces = []
        # direction = 3 for right
        d = 3
        for j in range(0, m):
            faces.append(d + 4 * (n * j + n - 1))
        return faces

    def plot_grid(self):
        """
        Plot the grid with cell number
        :return:
        """
        board = np.zeros((self.width, self.height))
        cmap = colors.ListedColormap(['white', 'red'])
        bounds = [0, 0.5, 1]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        ind_array_x = np.arange(0.5, self.width + 0.5, 1.0)
        ind_array_y = np.arange(0.5, self.height + 0.5, 1.0)
        x, y = np.meshgrid(ind_array_x, ind_array_y)
        plt.matshow(board, cmap=cmap, norm=norm, interpolation='none', vmin=0, vmax=1)
        for i, (x_val, y_val) in enumerate(zip(x.flatten(), y.flatten())):
            # c = str(int(x_val)) + ' ' + str(int(y_val))
            c = str(i)
            plt.text(x_val, y_val, c, va='center', ha='center')
        plt.xlim(0, self.width)
        plt.ylim(0, self.height)
        plt.xticks(np.arange(self.width))
        plt.yticks(np.arange(self.height))
        plt.grid()
        plt.show()


class Boundary:
    def __init__(self, face_no, boundary_type=None, value=None):
        """
        Class Boundary represents a single boundary on a cell.
        :param face_no: type int (see class Face in 'grid' module)
        :param boundary_type: type of boundary
        :param value: type float (velocity or pressure value)
        """
        self.face = face_no
        self.type = boundary_type
        self.value = value  # velocity or pressure
        self.position = self.face % 4


class LatticeVelocity(object):
    def __init__(self, x, y):
        """
        Velocity vector class
        :param x: velocity value in x direction
        :param y: velocity value in y direction
        """
        self.x = x
        self.y = y


def default_poiseuille_boundaries(grid):
    boundaries = []
    for i in range(0, grid.width):
        boundaries.append(Boundary(grid.bottom_faces[i], 'periodic'))
        boundaries.append(Boundary(grid.top_faces[i], 'periodic'))
    for i in range(0, grid.height):
        boundaries.append(Boundary(grid.left_faces[i], 'bounceback'))
        boundaries.append(Boundary(grid.right_faces[i], 'bounceback'))
    return boundaries


if __name__ == "__main__":
    gr = Grid(5, 10)
    gr.plot_grid()