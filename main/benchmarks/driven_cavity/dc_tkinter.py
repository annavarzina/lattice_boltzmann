from Tkinter import *
import main.schemes.d2q9 as lb
import main.schemes.grid as gd
import main.schemes.plot_and_save as ps
import os, inspect
import numpy as np


class EntryWindow:
    def __init__(self):
        self.master = Tk()
        self.errormsg = Label(text='', fg='red')
        self.errormsg.pack()
        row1 = Frame(self.master)
        row1.pack(side=TOP, fill=X, padx=5, pady=5)
        row2 = Frame(self.master)
        row2.pack(side=TOP, fill=X, padx=5, pady=5)
        row3 = Frame(self.master)
        row3.pack(side=TOP, fill=X, padx=5, pady=5)
        self.velocity_entry = Entry(row1,
                                    textvariable=StringVar(self.master, value='0.1'))
        self.viscosity_entry = Entry(row2,
                                     textvariable=StringVar(self.master, value='0.1'))
        self.length_entry = Entry(row3,
                                  textvariable=StringVar(self.master, value='100'))
        Label(row1, text='Velocity', width=15).pack(side=LEFT)
        Label(row2, text='Viscosity', width=15).pack(side=LEFT)
        Label(row3,text='Length', width=15).pack(side=LEFT)
        self.reynolds_label = Label(text='')
        self.velocity_entry = Entry(row1,
                                    textvariable=StringVar(self.master, value='0.1'),
                                    validate="focusout",
                                    validatecommand=(self.master.register(self.validate_number), '%P', 'velocity'),
                                    invalidcommand=(self.master.register(self.invalid_number), '%P', 'velocity'))
        self.velocity_entry.pack(side=RIGHT)
        self.viscosity_entry = Entry(row2,
                                     textvariable=StringVar(self.master, value='0.1'),
                                     validate="focusout",
                                     validatecommand=(self.master.register(self.validate_number), '%P', 'viscosity'),
                                     invalidcommand=(self.master.register(self.invalid_number), '%P', 'viscosity'))
        self.viscosity_entry.pack(side=RIGHT)
        self.length_entry = Entry(row3,
                                  textvariable=StringVar(self.master, value='100'),
                                  validate="focusout",
                                  validatecommand=(self.master.register(self.validate_number), '%P', 'length'),
                                  invalidcommand=(self.master.register(self.invalid_number), '%P', 'length'))
        self.length_entry.pack(side=RIGHT)
        self.reynolds_label.pack()
        Button(self.master, text="OK", command=self.create_window).pack()
        mainloop()

    def validate_number(self, p, parameter):
        self.errormsg.config(text='')
        if parameter == 'velocity':
            x = False
            try:
                x = (float(p) <= 0.2)
                ren = float(self.length_entry.get()) * \
                      float(self.velocity_entry.get()) / \
                      float(self.viscosity_entry.get())
                self.reynolds_label.config(text='Re: ' + str(ren))
                # print str(ren) + ' vel'
            except ValueError:
                pass
            return x  # always return True or False
        if parameter == 'viscosity':
            x = False
            try:
                x = (float(p) >= 0.03)
                ren = float(self.length_entry.get()) * \
                     float(self.velocity_entry.get()) / \
                     float(self.viscosity_entry.get())
                self.reynolds_label.config(text='Re: ' + str(ren))
                # print str(ren) + ' vis'
            except ValueError:
                pass
            return x  # always return True or False
        if parameter == 'length':
            x = False
            try:
                x = (float(p) <= 500)
                ren = float(self.length_entry.get()) * \
                     float(self.velocity_entry.get()) / \
                     float(self.viscosity_entry.get())
                # print str(ren) + ' len'
                self.reynolds_label.config(text='Re: ' + str(ren))
            except ValueError:
                pass
            return x  # always return True or False

    def invalid_number(self, p, parameter):
        if parameter == 'velocity':
            try:
                if float(p) > 0.2:
                    self.errormsg.config(text='Velocity is < 0.2')
            except ValueError:
                self.errormsg.config(text='Not a number')
        if parameter == 'viscosity':
            try:
                if float(p) < 0.03:
                    self.errormsg.config(text='Viscosity is < 0.03')
            except ValueError:
                self.errormsg.config(text='Not a number')
        if parameter == 'length':
            try:
                if float(p) > 500:
                    self.errormsg.config(text='Length is > 500')
            except ValueError:
                self.errormsg.config(text='Not a number')

    def create_window(self):
        entries = (('Velocity', float(self.velocity_entry.get())),
                   ('Viscosity', float(self.viscosity_entry.get())),
                   ('Length', int(self.length_entry.get())))
        self.master.destroy()
        new_window = MainWindow(entries)


class MainWindow(object):
    def __init__(self, entries):

        self.start_driven_cavity(entries)

    def start_driven_cavity(self, entries):
        print entries
        n = entries[2][1]  # width  - 0x
        m = entries[2][1]  # length - 0y
        v = entries[1][1]  # viscosity
        tau = self.calc_tau(v)
        t_f = 301  # final time
        u_x = entries[0][1]
        reynolds = n * u_x / v
        velocity = gd.LatticeVelocity(u_x, 0)
        user_grid = gd.Grid(n, m)
        # boundary conditions
        bc = self.driven_cavity_boundaries(user_grid, velocity)

        lat_bol = lb.D2Q9(grid=user_grid,
                          iterations=t_f,
                          boundary=bc,
                          relaxation_time=tau,
                          reynolds=reynolds)

    def calc_tau(self, viscosity):
        tau = viscosity * 3. + 0.5
        return tau

    def driven_cavity_boundaries(self,grid, value):
        boundaries = []
        v0 = gd.LatticeVelocity(0, 0)
        for i in range(0, grid.width):
            boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=v0))
            boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=value))
        for i in range(0, grid.height):
            boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=v0))
            boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_velocity', value=v0))
        return boundaries

if __name__ == "__main__":
    kv = EntryWindow()

