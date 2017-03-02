"""
Plot functions
Save functions
"""

import os, inspect
import csv
import numpy as np
import matplotlib.pyplot as plt

# ==================== plot functions ==================== #

def plot_ux(solution):
    plt.figure(1)
    plt.imshow(solution.u.x, origin='lower')
    # mat_plot.imshow(self.u.x)
    plt.colorbar()
    plt.pause(.001)
    plt.clf()


def plot_uy(solution):
    plt.figure(2)
    plt.imshow(solution.u.y, origin='lower')
    plt.colorbar()
    plt.pause(.001)
    plt.clf()


def plot_streamlines(solution, title=None, path=None):
    """ Print and plot streamlines """
    x = np.arange(0, solution.grid.width, solution.grid.dx)
    y = np.arange(0, solution.grid.height, solution.grid.dx)
    fig, ax = plt.subplots()
    plt.xticks(np.arange(0, solution.grid.width, solution.grid.width / 5))
    plt.yticks(np.arange(0, solution.grid.height, solution.grid.height / 5))
    title = 'Streamlines' + str(title)
    plt.title(title)
    strm = ax.streamplot(x, y, solution.u.x, solution.u.y,
                         color=solution.u.x,
                         linewidth=0.5,
                         density=4)
    plt.colorbar(strm.lines)
    if path is None:
        path = get_path()
    path += '\\fig'
    filename = 'figure_' + str(title) + '.png'
    filename = os.path.join(path, filename)
    plt.savefig(filename)
    plt.pause(.001)
    plt.show()


def plot_poiseuille_solution_0y(solution):
    l_y = solution.grid.height
    d = 1
    y = np.arange(0, l_y, d)
    ux = (solution.applied_force.x / (2 * solution.viscosity)) * y * (l_y - 1 - y)
    fig, ax1 = plt.subplots()
    ax1.plot(solution.u.x, 'r-', ux, 'ro')
    ax1.set_xlabel('y')
    ax1.set_ylabel('u')
    ax1.set_title('velocity profile along 0y')
    plt.tight_layout()
    plt.show()


def plot_poiseuille_solution_0x(solution):
    l_x = solution.grid.width
    d = 1
    x = np.arange(0, l_x, d)
    uy = (solution.applied_force.y / (2 * solution.viscosity)) * x * (l_x - 1 - x)
    fig, ax1 = plt.subplots()
    ax1.plot(solution.u.y[0, :], 'r-', uy, 'ro')
    ax1.set_xlabel('x')
    ax1.set_ylabel('u')
    ax1.set_title('velocity profile along 0x')
    plt.tight_layout()
    plt.show()


# ==================== save functions ==================== #


def get_path():
    # return os.path.realpath(__file__)
    return os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))


def save_velocity_x(solution, title=None, path=None):
    filename = 'velocity_x_' + str(title) + '.csv'
    if path is None:
        path = get_path()
    path = path + '\output'
    filename = os.path.join(path, filename)
    csvout = csv.writer(open(filename, "wb"))
    csvout.writerow((np.arange(0, 1, 1./solution.grid.width)))
    [csvout.writerow(r) for r in solution.u.x]


def save_velocity_y(solution, title=None, path=None):
    filename = 'velocity_y_' + str(title) + '.csv'
    if path is None:
        path = get_path()
    path = path + '\output'
    filename = os.path.join(path, filename)
    csvout = csv.writer(open(filename, "wb"))
    csvout.writerow((np.arange(0, solution.grid.width)))  # np.arange(0, 1, 1. / solution.grid.width
    [csvout.writerow(r) for r in solution.u.y]


def save_density(solution, title=None, path=None):
    filename = 'density_' + str(title) + '.csv'
    if path is None:
        path = get_path()
    path = path + '\output'
    filename = os.path.join(path, filename)
    csvout = csv.writer(open(filename, "wb"))
    csvout.writerow((np.arange(0, solution.grid.width)))  # np.arange(0, 1, 1. / solution.grid.width
    [csvout.writerow(r) for r in solution.rho]


def save_all(solution, title=None, path=None):
    save_velocity_x(solution, title, path=path)
    save_velocity_y(solution, title, path=path)
    save_density(solution, title, path=path)

    # fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    # ax1.plot(...)
    # ax2.plot(...)
    # plt.tight_layout()
    # plt.show()