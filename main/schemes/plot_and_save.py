"""
Plot functions
Save functions
"""

import os, inspect
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches
from numpy.matlib import rand
from pylab import *
def get_path():
    # return os.path.realpath(__file__)
    return os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))


# ==================== plot functions ==================== #


def plot_ux(solution, title=None, path=None, save=False):
    plt.figure(1)
    plt.imshow(solution.u.x, origin='lower')
    plt.colorbar()
    y = np.nonzero(solution.solid)[0]
    x = np.nonzero(solution.solid)[1]
    plt.scatter(x,y,color='black', s=100)
    plt.pause(.001)
    if save:
        if path is None:
            path = get_path()
        path += '\\fig'
        filename = 'ux' + str(title) + '.png'
        filename = os.path.join(path, filename)
        plt.savefig(filename)
    plt.clf()


def plot_uy(solution, title=None, path=None, save=False):
    plt.figure(2)
    plt.imshow(solution.u.y, origin='lower')
    plt.colorbar()
    y = np.nonzero(solution.solid)[0]
    x = np.nonzero(solution.solid)[1]
    plt.scatter(x,y,color='black', s=100)
    plt.pause(.001)
    if save:
        if path is None:
            path = get_path()
        path += '\\fig'
        filename = 'uy' + str(title) + '.png'
        filename = os.path.join(path, filename)
        plt.savefig(filename)
    plt.clf()


def plot_streamlines(solution, title=None, path=None, save=True):
    """ Print and plot streamlines """
    x = np.arange(0, solution.grid.width, solution.grid.dx)
    y = np.arange(0, solution.grid.height, solution.grid.dx)
    plt.figure(3)
    plt.clf()
    ax = plt.subplot()
    plt.xticks(np.arange(0, solution.grid.width, solution.grid.width / 5))
    plt.yticks(np.arange(0, solution.grid.height, solution.grid.height / 5))
    title = 'Streamlines' + str(title)
    plt.title(title)
    strm = ax.streamplot(x, y, solution.u.x, solution.u.y,
                         color=solution.u.x,
                         linewidth=0.5,
                         density=4)
    plt.colorbar(strm.lines)
    plt.pause(.001)
    if save:
        if path is None:
            path = get_path()
        path += '\\fig'
        filename = 'figure_' + str(title) + '.png'
        filename = os.path.join(path, filename)
        plt.savefig(filename)
    # plt.close(fig)


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


def close_figures():
    plt.hold(False)


# ==================== save functions ==================== #


def save_data(solution, a, parameter_name=None, title=None, path=None):
    filename = str(parameter_name) + '_' + str(title) + '.csv'
    if path is None:
        path = get_path()
    path += '\output'
    filename = os.path.join(path, filename)
    csvout = csv.writer(open(filename, "wb"))
    csvout.writerow((np.arange(0, 1, 1./solution.grid.width)))
    [csvout.writerow(r) for r in a]


def save_velocity_x(solution, title=None, path=None):
    save_data(solution, solution.u.x, 'velocity_x', title, path)


def save_velocity_y(solution, title=None, path=None):
    save_data(solution, solution.u.y, 'velocity_y', title, path)


def save_density(solution, title=None, path=None):
    save_data(solution, solution.rho, 'density', title, path)


def save_all(solution, title=None, path=None):
    save_velocity_x(solution, title, path=path)
    save_velocity_y(solution, title, path=path)
    save_density(solution, title, path=path)

    # fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    # ax1.plot(...)
    # ax2.plot(...)
    # plt.tight_layout()
    # plt.show()


# ==================== video functions ==================== #

# def lb_movie(solution, path):
#     pass
    # plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
    #
    # fig = plt.figure(1)
    # filename = 'demo.mp4'
    # filename = os.path.join(path, filename)
    #
    # ax = fig.add_subplot(111)
    # ax.set_aspect('equal')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    # im = ax.imshow(solution.u.x, origin='lower')
    # # im.set_clim([0, 1])
    # # fig.set_size_inches([5, 5])
    #
    # tight_layout()
    #
    # def animate(i):
    #     lb algorythm (i)
    #     return  value
    #
    # # legend(loc=0)
    # ani = animation.FuncAnimation(fig, animate, 300, interval=30)
    # writer = animation.FFMpegWriter()
    #
    # ani.save(filename, writer=writer, dpi=100)